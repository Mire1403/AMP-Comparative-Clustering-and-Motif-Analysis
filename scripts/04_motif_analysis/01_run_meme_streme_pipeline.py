from __future__ import annotations

from pathlib import Path
import csv
import shutil
import subprocess
import sys

from scipy.stats import ks_2samp

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


print("\nMOTIF DISCOVERY PIPELINE (STREME + MEME) using precomputed 10x backgrounds\n")

# =====================================================
# PATH CONFIG
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

RESULTS_CLUSTER_DIR = PROJECT_ROOT / "results" / "02_clustering"
RESULTS_BG_DIR = PROJECT_ROOT / "results" / "03_background"
RESULTS_MOTIF_DIR = PROJECT_ROOT / "results" / "04_motif_discovery"

RESULTS_MOTIF_DIR.mkdir(parents=True, exist_ok=True)

# Inputs
AMP_CDHIT = RESULTS_CLUSTER_DIR / "AMP_MASTER_cdhit.fasta"
AMP_MMSEQ = RESULTS_CLUSTER_DIR / "AMP_MASTER_mmseq_rep_seq.fasta"

NONAMP10_CDHIT = RESULTS_BG_DIR / "nonamp_cdhit_10x_len_kr_matched_progressive_random.fasta"
NONAMP10_MMSEQ = RESULTS_BG_DIR / "nonamp_mmseq_10x_len_kr_matched_progressive_random.fasta"

# =====================================================
# CONFIG
# =====================================================

ALPHA = 0.05

# STREME parameters
MINW = "6"
MAXW = "20"
STREME_NMOTIFS = "45"
STREME_PVT = "0.01"

# MEME parameters
MEME_NMOTIFS = "45"
MEME_MAXSIZE = "10000000"

# Runtime behavior
CLEAN_OUTPUT_DIRS_BEFORE_RUN = True

# =====================================================
# HELPERS
# =====================================================

def die(msg: str) -> None:
    print(f"❌ {msg}")
    sys.exit(1)


def check_file(path: Path) -> None:
    if not path.exists():
        die(f"File not found: {path}")


def get_binary(name: str) -> str:
    """
    Find an executable in PATH in a repo-friendly way.
    """
    path = shutil.which(name)
    if path is None:
        die(
            f"Required executable '{name}' not found in PATH.\n"
            f"Activate the correct environment or install MEME Suite.\n"
            f"Example checks:\n"
            f"  which {name}\n"
            f"  {name} --version"
        )
    return path


def check_binaries() -> tuple[str, str]:
    streme_bin = get_binary("streme")
    meme_bin = get_binary("meme")

    print("Detected binaries:")
    print(f"- STREME: {streme_bin}")
    print(f"- MEME:   {meme_bin}")

    for label, binary in [("STREME", streme_bin), ("MEME", meme_bin)]:
        try:
            result = subprocess.run(
                [binary, "--version"],
                capture_output=True,
                text=True,
                check=False,
            )
            version_text = (result.stdout or result.stderr).strip().splitlines()
            version_line = version_text[0] if version_text else "version info not available"
            print(f"  {label} version: {version_line}")
        except Exception:
            print(f"  {label} version: could not retrieve")

    print()
    return streme_bin, meme_bin


def maybe_reset_dir(path: Path) -> None:
    if path.exists() and CLEAN_OUTPUT_DIRS_BEFORE_RUN:
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def run_command(cmd: list[str], label: str) -> None:
    print("\nRunning:", " ".join(str(c) for c in cmd))
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        die(f"{label} failed with return code {e.returncode}.")
    except FileNotFoundError:
        die(f"Executable not found while running {label}: {cmd[0]}")


def read_fasta_lengths(path: Path) -> list[int]:
    lengths: list[int] = []
    seq_len = 0

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_len:
                    lengths.append(seq_len)
                seq_len = 0
            else:
                seq_len += len(line)

        if seq_len:
            lengths.append(seq_len)

    return lengths


def percent_kr_seq(seq: str) -> float:
    return (seq.count("K") + seq.count("R")) / len(seq) if seq else 0.0


def read_fasta_kr(path: Path) -> list[float]:
    krs: list[float] = []
    seq_parts: list[str] = []

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if seq_parts:
                    seq = "".join(seq_parts).upper()
                    krs.append(percent_kr_seq(seq))
                seq_parts = []
            else:
                seq_parts.append(line)

        if seq_parts:
            seq = "".join(seq_parts).upper()
            krs.append(percent_kr_seq(seq))

    return krs


def write_validation_csv(
    out_csv: Path,
    label: str,
    amp_fasta: Path,
    nonamp_fasta: Path,
    ks_len_pvalue: float,
    ks_kr_pvalue: float,
) -> None:
    with open(out_csv, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "label",
                "amp_fasta",
                "nonamp_fasta",
                "ks_length_pvalue",
                "ks_kr_pvalue",
                "alpha",
                "length_match_ok",
                "kr_match_ok",
            ],
            delimiter=";",
        )
        writer.writeheader()
        writer.writerow({
            "label": label,
            "amp_fasta": str(amp_fasta),
            "nonamp_fasta": str(nonamp_fasta),
            "ks_length_pvalue": ks_len_pvalue,
            "ks_kr_pvalue": ks_kr_pvalue,
            "alpha": ALPHA,
            "length_match_ok": ks_len_pvalue >= ALPHA,
            "kr_match_ok": ks_kr_pvalue >= ALPHA,
        })


def validate_background(amp_fasta: Path, nonamp_fasta: Path, label: str) -> None:
    """
    Quick sanity check: KS test for length and KR distributions.
    Saves plots and a small CSV summary.
    """
    amp_len = read_fasta_lengths(amp_fasta)
    bg_len = read_fasta_lengths(nonamp_fasta)

    amp_kr = read_fasta_kr(amp_fasta)
    bg_kr = read_fasta_kr(nonamp_fasta)

    ks_len = ks_2samp(amp_len, bg_len)
    ks_kr = ks_2samp(amp_kr, bg_kr)

    print(f"\nValidation ({label})")
    print(f"- KS length p-value: {ks_len.pvalue}")
    print(f"- KS KR p-value:     {ks_kr.pvalue}")

    if ks_len.pvalue < ALPHA or ks_kr.pvalue < ALPHA:
        print(
            "⚠️ Warning: KS suggests distributions differ (alpha=0.05). "
            "Length may be matched but KR may still differ."
        )

    outdir = RESULTS_MOTIF_DIR / f"validation_{label}"
    maybe_reset_dir(outdir)

    write_validation_csv(
        out_csv=outdir / "validation_summary.csv",
        label=label,
        amp_fasta=amp_fasta,
        nonamp_fasta=nonamp_fasta,
        ks_len_pvalue=ks_len.pvalue,
        ks_kr_pvalue=ks_kr.pvalue,
    )

    plt.figure()
    plt.hist(amp_len, bins=30, density=True, alpha=0.5, label="AMP")
    plt.hist(bg_len, bins=30, density=True, alpha=0.5, label="NONAMP")
    plt.legend()
    plt.xlabel("Length")
    plt.ylabel("Density")
    plt.title(f"{label}: Length (AMP vs NONAMP)")
    plt.tight_layout()
    plt.savefig(outdir / "length.png", dpi=300)
    plt.close()

    plt.figure()
    plt.hist(amp_kr, bins=30, density=True, alpha=0.5, label="AMP")
    plt.hist(bg_kr, bins=30, density=True, alpha=0.5, label="NONAMP")
    plt.legend()
    plt.xlabel("K+R proportion")
    plt.ylabel("Density")
    plt.title(f"{label}: KR (AMP vs NONAMP)")
    plt.tight_layout()
    plt.savefig(outdir / "kr.png", dpi=300)
    plt.close()


# =====================================================
# MOTIF DISCOVERY
# =====================================================

def run_motif(
    label: str,
    amp_fasta: Path,
    nonamp_fasta: Path,
    streme_bin: str,
    meme_bin: str,
) -> None:
    validate_background(amp_fasta, nonamp_fasta, label)

    # STREME
    streme_out = RESULTS_MOTIF_DIR / f"streme_{label}"
    maybe_reset_dir(streme_out)

    run_command(
        [
            streme_bin,
            "--p", str(amp_fasta),
            "--n", str(nonamp_fasta),
            "--protein",
            "--minw", MINW,
            "--maxw", MAXW,
            "--nmotifs", STREME_NMOTIFS,
            "--pvt", STREME_PVT,
            "--oc", str(streme_out),
        ],
        label=f"STREME ({label})",
    )
    print("✅ STREME completed:", streme_out)

    # MEME
    meme_out = RESULTS_MOTIF_DIR / f"meme_{label}"
    maybe_reset_dir(meme_out)

    run_command(
        [
            meme_bin,
            str(amp_fasta),
            "-neg", str(nonamp_fasta),
            "-mod", "zoops",
            "-objfun", "de",
            "-nmotifs", MEME_NMOTIFS,
            "-minw", MINW,
            "-maxw", MAXW,
            "-maxsize", MEME_MAXSIZE,
            "-oc", str(meme_out),
        ],
        label=f"MEME ({label})",
    )
    print("✅ MEME completed:", meme_out)


def main() -> None:
    for p in [AMP_CDHIT, AMP_MMSEQ, NONAMP10_CDHIT, NONAMP10_MMSEQ]:
        check_file(p)

    streme_bin, meme_bin = check_binaries()

    run_motif("cdhit", AMP_CDHIT, NONAMP10_CDHIT, streme_bin, meme_bin)
    run_motif("mmseq", AMP_MMSEQ, NONAMP10_MMSEQ, streme_bin, meme_bin)

    print("\n✅ MOTIF DISCOVERY PIPELINE COMPLETED.\n")


if __name__ == "__main__":
    main()