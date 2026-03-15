#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.stats import ks_2samp

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


print("\nMOTIF DISCOVERY PIPELINE (STREME + MEME) using precomputed 10x backgrounds\n")

# =====================================================
# PATH CONFIG
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

RESULTS_CLUSTER_DIR = PROJECT_ROOT / "results" / "clustering"
RESULTS_BG_DIR = PROJECT_ROOT / "results" / "background_generation"
RESULTS_MOTIF_DIR = PROJECT_ROOT / "results" / "motif_discovery"
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
STREME_NMOTIFS = "25"
STREME_PVT = "0.01"

# MEME parameters
MEME_NMOTIFS = "20"
MEME_MAXSIZE = "10000000"


# =====================================================
# HELPERS
# =====================================================

def die(msg: str) -> None:
    print(f"❌ {msg}")
    sys.exit(1)

def check_file(path: Path) -> None:
    if not path.exists():
        die(f"File not found: {path}")

def run_command(cmd: list[str]) -> None:
    print("\nRunning:", " ".join(str(c) for c in cmd))
    subprocess.run(cmd, check=True)

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

def percent_KR_seq(seq: str) -> float:
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
                    s = "".join(seq_parts).upper()
                    krs.append(percent_KR_seq(s))
                seq_parts = []
            else:
                seq_parts.append(line)
        if seq_parts:
            s = "".join(seq_parts).upper()
            krs.append(percent_KR_seq(s))
    return krs

def validate_background(amp_fasta: Path, nonamp_fasta: Path, label: str) -> None:
    """Quick sanity check: KS test for length + KR distributions."""
    amp_len = read_fasta_lengths(amp_fasta)
    bg_len = read_fasta_lengths(nonamp_fasta)

    amp_kr = read_fasta_kr(amp_fasta)
    bg_kr = read_fasta_kr(nonamp_fasta)

    ks_len = ks_2samp(amp_len, bg_len)
    ks_kr = ks_2samp(amp_kr, bg_kr)

    print(f"\nValidation ({label})")
    print(f"- KS length p-value: {ks_len.pvalue}")
    print(f"- KS KR p-value:     {ks_kr.pvalue}")

    # No mato el pipeline si falla: solo aviso (porque ya lo validaste en el paso anterior)
    if ks_len.pvalue < ALPHA or ks_kr.pvalue < ALPHA:
        print("⚠️ Warning: KS suggests distributions differ (alpha=0.05). Check your background step if needed.")

    # Save quick plots (optional, light)
    outdir = RESULTS_MOTIF_DIR / f"validation_{label}"
    outdir.mkdir(parents=True, exist_ok=True)

    plt.figure()
    plt.hist(amp_len, bins=30, density=True, alpha=0.5, label="AMP")
    plt.hist(bg_len, bins=30, density=True, alpha=0.5, label="NONAMP")
    plt.legend()
    plt.xlabel("Length")
    plt.ylabel("Density")
    plt.title(f"{label}: Length (AMP vs NONAMP)")
    plt.tight_layout()
    plt.savefig(outdir / "length.png")
    plt.close()

    plt.figure()
    plt.hist(amp_kr, bins=30, density=True, alpha=0.5, label="AMP")
    plt.hist(bg_kr, bins=30, density=True, alpha=0.5, label="NONAMP")
    plt.legend()
    plt.xlabel("K+R proportion")
    plt.ylabel("Density")
    plt.title(f"{label}: KR (AMP vs NONAMP)")
    plt.tight_layout()
    plt.savefig(outdir / "kr.png")
    plt.close()


# =====================================================
# MOTIF DISCOVERY
# =====================================================

def run_motif(label: str, amp_fasta: Path, nonamp_fasta: Path) -> None:
    validate_background(amp_fasta, nonamp_fasta, label)

    # STREME
    streme_out = RESULTS_MOTIF_DIR / f"streme_{label}"
    streme_out.mkdir(parents=True, exist_ok=True)

    run_command([
        "streme",
        "--p", str(amp_fasta),
        "--n", str(nonamp_fasta),
        "--protein",
        "--minw", MINW,
        "--maxw", MAXW,
        "--nmotifs", STREME_NMOTIFS,
        "--pvt", STREME_PVT,
        "--oc", str(streme_out),
    ])
    print("✅ STREME completed:", streme_out)

    # MEME
    meme_out = RESULTS_MOTIF_DIR / f"meme_{label}"
    meme_out.mkdir(parents=True, exist_ok=True)

    run_command([
        "meme",
        str(amp_fasta),
        "-neg", str(nonamp_fasta),
        "-mod", "zoops",
        "-objfun", "de",
        "-nmotifs", MEME_NMOTIFS,
        "-minw", MINW,
        "-maxw", MAXW,
        "-maxsize", MEME_MAXSIZE,
        "-oc", str(meme_out),
    ])
    print("✅ MEME completed:", meme_out)


def main() -> None:
    for p in [AMP_CDHIT, AMP_MMSEQ, NONAMP10_CDHIT, NONAMP10_MMSEQ]:
        check_file(p)

    run_motif("cdhit", AMP_CDHIT, NONAMP10_CDHIT)
    run_motif("mmseq", AMP_MMSEQ, NONAMP10_MMSEQ)

    print("\n✅ MOTIF DISCOVERY PIPELINE COMPLETED.\n")


if __name__ == "__main__":
    main()