from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# =====================================================
# PATH CONFIG (RELATIVE TO REPO)
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

RESULTS_CLUSTER_DIR = PROJECT_ROOT / "results" / "clustering"
RESULTS_CLUSTER_DIR.mkdir(parents=True, exist_ok=True)

CDHIT_FASTA = RESULTS_CLUSTER_DIR / "AMP_MASTER_cdhit80.fasta"
MMSEQ_FASTA = RESULTS_CLUSTER_DIR / "AMP_MASTER_mmseq80_rep_seq.fasta"

OUTPUT_DIR = RESULTS_CLUSTER_DIR / "comparison"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# =====================================================
# FUNCTIONS
# =====================================================

def read_fasta_lengths(fasta_path):

    if not fasta_path.exists():
        print(f"File not found: {fasta_path}")
        return []

    lengths = []
    seq = ""

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq:
                    lengths.append(len(seq))
                    seq = ""
            else:
                seq += line

        if seq:
            lengths.append(len(seq))

    return lengths


def compute_stats(lengths):

    if not lengths:
        return {}

    return {
        "count": len(lengths),
        "min": min(lengths),
        "max": max(lengths),
        "mean": round(np.mean(lengths), 2),
        "std": round(np.std(lengths), 2),
    }


def plot_histogram(lengths, title, output_path):

    if not lengths:
        print(f"No data to plot for {title}")
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

if __name__ == "__main__":

    print("\n=== COMPARING CD-HIT VS MMSEQS CLUSTER REPRESENTATIVES ===\n")

    # Read FASTA files
    cdhit_lengths = read_fasta_lengths(CDHIT_FASTA)
    mmseq_lengths = read_fasta_lengths(MMSEQ_FASTA)

    # Compute statistics
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

    # Plot histograms
    plot_histogram(
        cdhit_lengths,
        "Length Distribution - CD-HIT",
        OUTPUT_DIR / "cdhit_length_hist.png"
    )

    plot_histogram(
        mmseq_lengths,
        "Length Distribution - MMseqs2",
        OUTPUT_DIR / "mmseq_length_hist.png"
    )

    print("\nHistograms saved in:", OUTPUT_DIR)
    print("\nComparison completed successfully.\n")