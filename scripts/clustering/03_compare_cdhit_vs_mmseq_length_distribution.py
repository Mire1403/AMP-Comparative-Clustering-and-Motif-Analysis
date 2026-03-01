from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

print("\n=== COMPARING CD-HIT VS MMSEQS CLUSTER REPRESENTATIVES ===\n")

# =====================================================
# PATH CONFIG
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

RESULTS_CLUSTER_DIR = PROJECT_ROOT / "results" / "clustering"
OUTPUT_DIR = RESULTS_CLUSTER_DIR / "comparison"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

CDHIT_FASTA = RESULTS_CLUSTER_DIR / "AMP_MASTER_cdhit.fasta"

# Default expected MMseq representative FASTA
MMSEQ_REP_FASTA = RESULTS_CLUSTER_DIR / "AMP_MASTER_mmseq_rep_seq.fasta"

# Fallback: search pattern if name differs
MMSEQ_REP_GLOB = "AMP_MASTER_mmseq*_rep_seq.fasta"

STATS_CSV = OUTPUT_DIR / "cluster_representatives_length_stats.csv"


# =====================================================
# HELPERS
# =====================================================

def die(msg: str) -> None:
    print(f" {msg}")
    sys.exit(1)


def resolve_mmseq_rep_fasta() -> Path:
    """Return the MMseq rep fasta path, trying a fallback search if needed."""
    if MMSEQ_REP_FASTA.exists():
        return MMSEQ_REP_FASTA

    candidates = sorted(RESULTS_CLUSTER_DIR.glob(MMSEQ_REP_GLOB))
    if candidates:
        print(f" Expected MMseq rep fasta not found: {MMSEQ_REP_FASTA.name}")
        print(f" Using fallback detected file: {candidates[0].name}")
        return candidates[0]

    die(
        "MMseq representative FASTA not found.\n"
        f"- Expected: {MMSEQ_REP_FASTA}\n"
        f"- Also tried glob: {RESULTS_CLUSTER_DIR / MMSEQ_REP_GLOB}"
    )


def read_fasta_lengths(fasta_path: Path) -> list[int]:
    if not fasta_path.exists():
        die(f"File not found: {fasta_path}")

    lengths: list[int] = []
    seq_len = 0

    with open(fasta_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_len > 0:
                    lengths.append(seq_len)
                seq_len = 0
            else:
                seq_len += len(line)

        if seq_len > 0:
            lengths.append(seq_len)

    return lengths


def compute_stats(lengths: list[int]) -> dict:
    if not lengths:
        return {"count": 0, "min": np.nan, "max": np.nan, "mean": np.nan, "std": np.nan}

    arr = np.array(lengths, dtype=float)
    return {
        "count": int(arr.size),
        "min": int(np.min(arr)),
        "max": int(np.max(arr)),
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr)),
    }


def plot_histogram(lengths: list[int], title: str, output_path: Path) -> None:
    if not lengths:
        print(f"⚠️ No data to plot for {title}")
        return

    plt.figure()
    plt.hist(lengths, bins=30)
    plt.xlabel("Sequence Length (aa)")
    plt.ylabel("Frequency")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


# =====================================================
# MAIN
# =====================================================

def main() -> None:
    if not CDHIT_FASTA.exists():
        die(f"CD-HIT FASTA not found: {CDHIT_FASTA}")

    mmseq_fasta = resolve_mmseq_rep_fasta()

    print("CD-HIT reps:", CDHIT_FASTA)
    print("MMseq reps :", mmseq_fasta)
    print("Output dir :", OUTPUT_DIR, "\n")

    cdhit_lengths = read_fasta_lengths(CDHIT_FASTA)
    mmseq_lengths = read_fasta_lengths(mmseq_fasta)

    cdhit_stats = compute_stats(cdhit_lengths)
    mmseq_stats = compute_stats(mmseq_lengths)

    print("CD-HIT STATS")
    print("------------------------------")
    for k, v in cdhit_stats.items():
        print(f"{k}: {v}")

    print("\nMMseqs STATS")
    print("------------------------------")
    for k, v in mmseq_stats.items():
        print(f"{k}: {v}")

    # Save stats table
    stats_df = pd.DataFrame([
        {"Method": "cdhit", **cdhit_stats, "fasta": str(CDHIT_FASTA)},
        {"Method": "mmseq", **mmseq_stats, "fasta": str(mmseq_fasta)},
    ])
    stats_df.to_csv(STATS_CSV, index=False)

    # Plots
    plot_histogram(
        cdhit_lengths,
        "Length Distribution - CD-HIT representatives",
        OUTPUT_DIR / "cdhit_length_hist.png",
    )
    plot_histogram(
        mmseq_lengths,
        "Length Distribution - MMseq representatives",
        OUTPUT_DIR / "mmseq_length_hist.png",
    )

    print("\n Saved:")
    print("-", OUTPUT_DIR / "cdhit_length_hist.png")
    print("-", OUTPUT_DIR / "mmseq_length_hist.png")
    print("-", STATS_CSV)
    print("\n Comparison completed successfully.\n")


if __name__ == "__main__":
    main()