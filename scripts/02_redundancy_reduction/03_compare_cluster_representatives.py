from __future__ import annotations

from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

print("\n=== CLUSTER REPRESENTATIVE ANALYSIS (4 CONDITIONS) ===\n")

# =====================================================
# PATH CONFIG
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

RESULTS_CLUSTER_DIR = PROJECT_ROOT / "results" / "clustering"
OUT_DIR = PROJECT_ROOT / "results" / "statistics" / "clustering"
OUT_DIR.mkdir(parents=True, exist_ok=True)

STATS_CSV = OUT_DIR / "cluster_representatives_length_stats.csv"

# Patterns (robust to small naming differences)
FILES = {
    "cdhit_amp": "AMP_MASTER_cdhit*.fasta",
    "mmseq_amp": "AMP_MASTER_mmseq*_rep_seq.fasta",
    "cdhit_nonamp": "NONAMP_cdhit*.fasta",
    "mmseq_nonamp": "NONAMP_mmseq*_rep_seq.fasta",
}

# =====================================================
# HELPERS
# =====================================================

def die(msg: str) -> None:
    print(f"❌ {msg}")
    sys.exit(1)


def resolve_one(glob_pat: str) -> Path:
    hits = sorted(RESULTS_CLUSTER_DIR.glob(glob_pat))
    if not hits:
        die(f"No file found for pattern: {glob_pat} in {RESULTS_CLUSTER_DIR}")
    if len(hits) > 1:
        print(f"⚠️ Multiple matches for '{glob_pat}', using: {hits[0].name}")
    return hits[0]


def read_fasta_lengths(path: Path) -> list[int]:
    if not path.exists():
        die(f"File not found: {path}")

    lengths: list[int] = []
    seq_len = 0

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
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
        return {
            "count": 0, "min": np.nan, "max": np.nan,
            "mean": np.nan, "std": np.nan,
            "p50": np.nan, "p90": np.nan
        }

    arr = np.asarray(lengths, dtype=float)
    return {
        "count": int(arr.size),
        "min": int(arr.min()),
        "max": int(arr.max()),
        "mean": float(arr.mean()),
        "std": float(arr.std()),
        "p50": float(np.percentile(arr, 50)),
        "p90": float(np.percentile(arr, 90)),
    }


def plot_hist(lengths: list[int], title: str, out_path: Path) -> None:
    if not lengths:
        print(f"⚠️ No data to plot for: {title}")
        return
    plt.figure()
    plt.hist(lengths, bins=30)
    plt.xlabel("Sequence length (aa)")
    plt.ylabel("Frequency")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


# =====================================================
# MAIN
# =====================================================

def main() -> None:
    rows = []

    for label, pattern in FILES.items():
        path = resolve_one(pattern)
        lengths = read_fasta_lengths(path)
        stats = compute_stats(lengths)

        rows.append({
            "condition": label,
            "fasta": str(path),
            **stats
        })

        plot_hist(
            lengths,
            f"{label} — representative length distribution",
            OUT_DIR / f"{label}_length_hist.png",
        )

        print(f"✅ {label}: n={stats['count']} -> {path.name}")

    df = pd.DataFrame(rows).sort_values("condition")
    df.to_csv(STATS_CSV, index=False)

    print("\n✅ Saved:")
    print("-", STATS_CSV)
    print("-", OUT_DIR, "(histograms)")

    print("\nDone.\n")


if __name__ == "__main__":
    main()