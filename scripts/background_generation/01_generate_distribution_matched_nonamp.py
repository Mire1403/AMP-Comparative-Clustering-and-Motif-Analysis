from pathlib import Path
import random
from collections import defaultdict
import numpy as np
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt

print("\nRANDOM PROGRESSIVE PIPELINE — MEMORY SAFE\n")

# =====================================================
# PATH CONFIG (RELATIVE TO REPO)
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

RESULTS_CLUSTER_DIR = PROJECT_ROOT / "results" / "clustering"
RESULTS_BG_DIR = PROJECT_ROOT / "results" / "background_generation"

RESULTS_BG_DIR.mkdir(parents=True, exist_ok=True)

AMP_FASTA_CDHIT = RESULTS_CLUSTER_DIR / "AMP_MASTER_cdhit80.fasta"
AMP_FASTA_MMSEQ = RESULTS_CLUSTER_DIR / "AMP_MASTER_mmseq80_rep_seq.fasta"

NONAMP_FASTA_CDHIT = RESULTS_CLUSTER_DIR / "non_amp" / "NONAMP_cdhit80.fasta"
NONAMP_FASTA_MMSEQ = RESULTS_CLUSTER_DIR / "non_amp" / "NONAMP_mmseq80_rep_seq.fasta"

OUTPUT_DIR = RESULTS_BG_DIR

MULTIPLIER = 50
MIN_LENGTH = 5
random.seed(42)

# =====================================================
# BASIC FUNCTIONS
# =====================================================

def percent_KR(seq):
    return (seq.count("K") + seq.count("R")) / len(seq)

def read_fasta(path):
    sequences = []
    with open(path) as f:
        seq = ""
        for line in f:
            if line.startswith(">"):
                if seq:
                    sequences.append(seq)
                seq = ""
            else:
                seq += line.strip()
        if seq:
            sequences.append(seq)

    # 🔥 Remove sequences < 5 aa
    initial_count = len(sequences)
    sequences = [s for s in sequences if len(s) >= MIN_LENGTH]
    removed = initial_count - len(sequences)

    print(f"{path.name} → removed (<5 aa): {removed}")

    return sequences

def compute_amp_stats(sequences):
    length_counts = defaultdict(int)
    kr_by_length = defaultdict(list)

    for seq in sequences:
        L = len(seq)
        length_counts[L] += 1
        kr_by_length[L].append(percent_KR(seq))

    kr_stats = {}
    for L in kr_by_length:
        kr_stats[L] = (np.mean(kr_by_length[L]), np.std(kr_by_length[L]))

    return length_counts, kr_stats

# =====================================================
# RANDOM PROGRESSIVE SLIDING WINDOW
# =====================================================

def progressive_random_sampling(nonamp_sequences, target_counts, kr_stats):

    selected_counts = defaultdict(int)
    selected_fragments = []

    lengths = list(target_counts.keys())
    total_required = sum(target_counts.values())

    random.shuffle(nonamp_sequences)

    for seq in nonamp_sequences:

        if len(selected_fragments) >= total_required:
            break

        seq_len = len(seq)
        random.shuffle(lengths)

        for L in lengths:

            if selected_counts[L] >= target_counts[L]:
                continue

            if seq_len < L:
                continue

            mean, sd = kr_stats[L]
            lower = mean - sd
            upper = mean + sd

            positions = list(range(seq_len - L + 1))
            random.shuffle(positions)

            for i in positions:

                if selected_counts[L] >= target_counts[L]:
                    break

                if len(selected_fragments) >= total_required:
                    break

                frag = seq[i:i+L]
                kr = percent_KR(frag)

                if lower <= kr <= upper:
                    selected_fragments.append(frag)
                    selected_counts[L] += 1

    print("Fragments generated:", len(selected_fragments))
    print("Expected:", total_required)

    return selected_fragments

# =====================================================
# VALIDATION
# =====================================================

def validate_and_plot(amp_seqs, nonamp_seqs, label):

    print(f"\n=== VALIDATION FOR {label.upper()} ===")

    amp_lengths = [len(s) for s in amp_seqs]
    nonamp_lengths = [len(s) for s in nonamp_seqs]

    amp_kr = [percent_KR(s) for s in amp_seqs]
    nonamp_kr = [percent_KR(s) for s in nonamp_seqs]

    ks_len = ks_2samp(amp_lengths, nonamp_lengths)
    ks_kr = ks_2samp(amp_kr, nonamp_kr)

    alpha = 0.05

    def interpret(p):
        return (
            "NOT significantly different (distributions statistically similar)"
            if p >= alpha
            else "Significantly different"
        )

    results_file = OUTPUT_DIR / f"{label}_validation_report.txt"

    with open(results_file, "w") as f:
        f.write("KOLMOGOROV–SMIRNOV TEST RESULTS\n\n")
        f.write("Length Distribution:\n")
        f.write(f"KS statistic: {ks_len.statistic}\n")
        f.write(f"p-value: {ks_len.pvalue}\n")
        f.write(f"Interpretation (alpha=0.05): {interpret(ks_len.pvalue)}\n\n")

        f.write("K+R Composition Distribution:\n")
        f.write(f"KS statistic: {ks_kr.statistic}\n")
        f.write(f"p-value: {ks_kr.pvalue}\n")
        f.write(f"Interpretation (alpha=0.05): {interpret(ks_kr.pvalue)}\n")

    print("Validation report saved:", results_file)

    # LENGTH
    plt.figure()
    plt.hist(amp_lengths, bins=30, density=True, alpha=0.5, label="AMP")
    plt.hist(nonamp_lengths, bins=30, density=True, alpha=0.5, label="NONAMP")
    plt.legend()
    plt.xlabel("Peptide Length")
    plt.ylabel("Density")
    plt.title(f"{label.upper()} — Length Distribution (AMP vs NONAMP)")
    plt.savefig(OUTPUT_DIR / f"{label}_length_superposed.png")
    plt.close()

    # KR
    plt.figure()
    plt.hist(amp_kr, bins=30, density=True, alpha=0.5, label="AMP")
    plt.hist(nonamp_kr, bins=30, density=True, alpha=0.5, label="NONAMP")
    plt.legend()
    plt.xlabel("K+R Proportion")
    plt.ylabel("Density")
    plt.title(f"{label.upper()} — KR Distribution (AMP vs NONAMP)")
    plt.savefig(OUTPUT_DIR / f"{label}_kr_superposed.png")
    plt.close()

    print("Superposed histograms saved.\n")

# =====================================================
# CDHIT80
# =====================================================

print("\n=== CDHIT80 ===")

amp_cdhit = read_fasta(AMP_FASTA_CDHIT)
nonamp_cdhit = read_fasta(NONAMP_FASTA_CDHIT)

length_counts_cdhit, kr_stats_cdhit = compute_amp_stats(amp_cdhit)

target_counts_cdhit = {L: count * MULTIPLIER for L, count in length_counts_cdhit.items()}

generated_cdhit = progressive_random_sampling(
    nonamp_cdhit,
    target_counts_cdhit,
    kr_stats_cdhit
)

with open(OUTPUT_DIR / "nonamp_cdhit80_50x_random_progressive.fasta", "w") as f:
    for i, frag in enumerate(generated_cdhit):
        f.write(f">NONAMP_frag_{i}\n{frag}\n")

validate_and_plot(amp_cdhit, generated_cdhit, "cdhit80")

# =====================================================
# MMSEQ80
# =====================================================

print("\n=== MMSEQ80 ===")

amp_mmseq = read_fasta(AMP_FASTA_MMSEQ)
nonamp_mmseq = read_fasta(NONAMP_FASTA_MMSEQ)

length_counts_mmseq, kr_stats_mmseq = compute_amp_stats(amp_mmseq)

target_counts_mmseq = {L: count * MULTIPLIER for L, count in length_counts_mmseq.items()}

generated_mmseq = progressive_random_sampling(
    nonamp_mmseq,
    target_counts_mmseq,
    kr_stats_mmseq
)

with open(OUTPUT_DIR / "nonamp_mmseq80_50x_random_progressive.fasta", "w") as f:
    for i, frag in enumerate(generated_mmseq):
        f.write(f">NONAMP_frag_{i}\n{frag}\n")

validate_and_plot(amp_mmseq, generated_mmseq, "mmseq80")

print("\nPIPELINE COMPLETE — RANDOM PROGRESSIVE VERSION.\n")