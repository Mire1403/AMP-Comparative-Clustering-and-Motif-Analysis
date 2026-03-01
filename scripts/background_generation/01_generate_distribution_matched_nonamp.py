#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import sys
import random
from collections import defaultdict

import numpy as np
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt

print("\nRANDOM PROGRESSIVE PIPELINE — MEMORY SAFE\n")

# =====================================================
# PATH CONFIG
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

RESULTS_CLUSTER_DIR = PROJECT_ROOT / "results" / "clustering"
RESULTS_BG_DIR = PROJECT_ROOT / "results" / "background_generation"
RESULTS_BG_DIR.mkdir(parents=True, exist_ok=True)

# Inputs produced by your UPDATED clustering scripts
AMP_FASTA_CDHIT = RESULTS_CLUSTER_DIR / "AMP_MASTER_cdhit.fasta"
AMP_FASTA_MMSEQ = RESULTS_CLUSTER_DIR / "AMP_MASTER_mmseq_rep_seq.fasta"

NONAMP_FASTA_CDHIT = RESULTS_CLUSTER_DIR / "NONAMP_cdhit.fasta"
NONAMP_FASTA_MMSEQ = RESULTS_CLUSTER_DIR / "NONAMP_mmseq_rep_seq.fasta"

OUTPUT_DIR = RESULTS_BG_DIR

# Config
MULTIPLIER = 50
MIN_LENGTH = 5
RANDOM_SEED = 42
random.seed(RANDOM_SEED)

# =====================================================
# BASIC FUNCTIONS
# =====================================================

def die(msg: str) -> None:
    print(f"❌ {msg}")
    sys.exit(1)

def percent_KR(seq: str) -> float:
    return (seq.count("K") + seq.count("R")) / len(seq)

def read_fasta(path: Path) -> list[str]:
    if not path.exists():
        die(f"FASTA not found: {path}")

    sequences: list[str] = []
    seq_parts: list[str] = []

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_parts:
                    sequences.append("".join(seq_parts).upper())
                seq_parts = []
            else:
                seq_parts.append(line)

        if seq_parts:
            sequences.append("".join(seq_parts).upper())

    # Remove sequences < MIN_LENGTH
    initial_count = len(sequences)
    sequences = [s for s in sequences if len(s) >= MIN_LENGTH]
    removed = initial_count - len(sequences)
    print(f"{path.name} → removed (<{MIN_LENGTH} aa): {removed}")

    return sequences

def compute_amp_stats(sequences: list[str]):
    length_counts = defaultdict(int)
    kr_by_length = defaultdict(list)

    for seq in sequences:
        L = len(seq)
        length_counts[L] += 1
        kr_by_length[L].append(percent_KR(seq))

    kr_stats = {}
    for L, vals in kr_by_length.items():
        kr_stats[L] = (float(np.mean(vals)), float(np.std(vals)))

    return length_counts, kr_stats

# =====================================================
# RANDOM PROGRESSIVE SLIDING WINDOW
# =====================================================

def progressive_random_sampling(nonamp_sequences, target_counts, kr_stats):
    selected_counts = defaultdict(int)
    selected_fragments: list[str] = []

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
        return "NOT significantly different" if p >= alpha else "Significantly different"

    results_file = OUTPUT_DIR / f"{label}_validation_report.txt"

    with open(results_file, "w", encoding="utf-8") as f:
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

    # Length plot
    plt.figure()
    plt.hist(amp_lengths, bins=30, density=True, alpha=0.5, label="AMP")
    plt.hist(nonamp_lengths, bins=30, density=True, alpha=0.5, label="NONAMP")
    plt.legend()
    plt.xlabel("Peptide Length")
    plt.ylabel("Density")
    plt.title(f"{label.upper()} — Length Distribution (AMP vs NONAMP)")
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / f"{label}_length_superposed.png")
    plt.close()

    # KR plot
    plt.figure()
    plt.hist(amp_kr, bins=30, density=True, alpha=0.5, label="AMP")
    plt.hist(nonamp_kr, bins=30, density=True, alpha=0.5, label="NONAMP")
    plt.legend()
    plt.xlabel("K+R Proportion")
    plt.ylabel("Density")
    plt.title(f"{label.upper()} — KR Distribution (AMP vs NONAMP)")
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / f"{label}_kr_superposed.png")
    plt.close()

    print("Superposed histograms saved.\n")

# =====================================================
# RUN
# =====================================================

def main():
    print("\n=== CDHIT ===")
    amp_cdhit = read_fasta(AMP_FASTA_CDHIT)
    nonamp_cdhit = read_fasta(NONAMP_FASTA_CDHIT)

    length_counts_cdhit, kr_stats_cdhit = compute_amp_stats(amp_cdhit)
    target_counts_cdhit = {L: count * MULTIPLIER for L, count in length_counts_cdhit.items()}

    generated_cdhit = progressive_random_sampling(nonamp_cdhit, target_counts_cdhit, kr_stats_cdhit)

    out_cdhit = OUTPUT_DIR / f"nonamp_cdhit_{MULTIPLIER}x_random_progressive.fasta"
    with open(out_cdhit, "w", encoding="utf-8") as f:
        for i, frag in enumerate(generated_cdhit, start=1):
            f.write(f">NONAMP_frag_{i}\n{frag}\n")

    validate_and_plot(amp_cdhit, generated_cdhit, "cdhit")

    print("\n=== MMSEQ ===")
    amp_mmseq = read_fasta(AMP_FASTA_MMSEQ)
    nonamp_mmseq = read_fasta(NONAMP_FASTA_MMSEQ)

    length_counts_mmseq, kr_stats_mmseq = compute_amp_stats(amp_mmseq)
    target_counts_mmseq = {L: count * MULTIPLIER for L, count in length_counts_mmseq.items()}

    generated_mmseq = progressive_random_sampling(nonamp_mmseq, target_counts_mmseq, kr_stats_mmseq)

    out_mmseq = OUTPUT_DIR / f"nonamp_mmseq_{MULTIPLIER}x_random_progressive.fasta"
    with open(out_mmseq, "w", encoding="utf-8") as f:
        for i, frag in enumerate(generated_mmseq, start=1):
            f.write(f">NONAMP_frag_{i}\n{frag}\n")

    validate_and_plot(amp_mmseq, generated_mmseq, "mmseq")

    print("\n✅ PIPELINE COMPLETE — RANDOM PROGRESSIVE VERSION.\n")

if __name__ == "__main__":
    main()