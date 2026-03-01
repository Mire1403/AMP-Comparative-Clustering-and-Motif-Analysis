#!/usr/bin/env python3

from pathlib import Path
import random
from collections import defaultdict
import numpy as np
from scipy.stats import ks_2samp
import subprocess
import sys

print("\nSMART 10x MOTIF PIPELINE\n")

# =====================================================
# PATH CONFIG (RELATIVE TO REPO)
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

RESULTS_CLUSTER_DIR = PROJECT_ROOT / "results" / "clustering"
RESULTS_BACKGROUND_DIR = PROJECT_ROOT / "results" / "background_generation"
RESULTS_MOTIF_DIR = PROJECT_ROOT / "results" / "motif_discovery"

RESULTS_MOTIF_DIR.mkdir(parents=True, exist_ok=True)

AMP_CDHIT = RESULTS_CLUSTER_DIR / "AMP_MASTER_cdhit80.fasta"
AMP_MMSEQ = RESULTS_CLUSTER_DIR / "AMP_MASTER_mmseq80_rep_seq.fasta"

NONAMP50_CDHIT = RESULTS_BACKGROUND_DIR / "nonamp_cdhit80_50x_random_progressive.fasta"
NONAMP50_MMSEQ = RESULTS_BACKGROUND_DIR / "nonamp_mmseq80_50x_random_progressive.fasta"

MULTIPLIER = 10
ALPHA = 0.05
random.seed(42)

# =====================================================
# UTILITIES
# =====================================================

def read_fasta(path):
    if not path.exists():
        print(f"File not found: {path}")
        sys.exit(1)

    seqs = []
    with open(path) as f:
        seq = ""
        for line in f:
            if line.startswith(">"):
                if seq:
                    seqs.append(seq)
                seq = ""
            else:
                seq += line.strip()
        if seq:
            seqs.append(seq)
    return seqs


def get_length_distribution(seqs):
    dist = defaultdict(int)
    for s in seqs:
        dist[len(s)] += 1
    return dist


def group_by_length(seqs):
    groups = defaultdict(list)
    for s in seqs:
        groups[len(s)].append(s)
    return groups


def run_command(cmd):
    print("\nRunning:", " ".join(str(c) for c in cmd))
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print("❌ ERROR DETECTED. Stopping pipeline.")
        sys.exit(1)

# =====================================================
# CORE GENERATOR
# =====================================================

def generate_10x(amp_fasta, nonamp50_fasta, label):

    print(f"\n=== GENERATING 10x FOR {label.upper()} ===")

    amp_seqs = read_fasta(amp_fasta)
    nonamp50_seqs = read_fasta(nonamp50_fasta)

    amp_dist = get_length_distribution(amp_seqs)
    nonamp_groups = group_by_length(nonamp50_seqs)

    target_counts = {L: amp_dist[L] * MULTIPLIER for L in amp_dist}

    selected = []

    for L in target_counts:

        available = nonamp_groups.get(L, [])
        required = target_counts[L]

        if len(available) < required:
            print(f"❌ Not enough sequences for length {L}")
            sys.exit(1)

        sampled = random.sample(available, required)
        selected.extend(sampled)

    print("Expected:", sum(target_counts.values()))
    print("Generated:", len(selected))

    # VALIDATION
    amp_lengths = [len(s) for s in amp_seqs]
    nonamp_lengths = [len(s) for s in selected]

    ks_len = ks_2samp(amp_lengths, nonamp_lengths)

    print("KS p-value:", ks_len.pvalue)

    if ks_len.pvalue < ALPHA:
        print("❌ Distributions significantly different. Stopping.")
        sys.exit(1)
    else:
        print("✅ Length distributions validated.")

    output_fasta = RESULTS_MOTIF_DIR / f"nonamp_{label}_10x.fasta"

    with open(output_fasta, "w") as f:
        for i, s in enumerate(selected):
            f.write(f">NONAMP_{i}\n{s}\n")

    return output_fasta

# =====================================================
# RUN PIPELINE FOR ONE DATASET
# =====================================================

def run_motif_pipeline(label, amp_fasta, nonamp50_fasta):

    nonamp10 = generate_10x(amp_fasta, nonamp50_fasta, label)

    # STREME
    print(f"\n🚀 Running STREME ({label})")

    streme_out = RESULTS_MOTIF_DIR / f"streme_{label}"

    run_command([
        "streme",
        "--p", str(amp_fasta),
        "--n", str(nonamp10),
        "--protein",
        "--minw", "6",
        "--maxw", "20",
        "--nmotifs", "15",
        "--pvt", "0.01",
        "--oc", str(streme_out)
    ])

    print("✅ STREME completed.")

    # MEME
    print(f"\n🚀 Running MEME ({label})")

    meme_out = RESULTS_MOTIF_DIR / f"meme_{label}"

    run_command([
        "meme",
        str(amp_fasta),
        "-neg", str(nonamp10),
        "-mod", "zoops",
        "-objfun", "de",
        "-nmotifs", "10",
        "-minw", "6",
        "-maxw", "20",
        "-maxsize", "10000000",
        "-oc", str(meme_out)
    ])

    print("✅ MEME completed.")

# =====================================================
# RUN BOTH
# =====================================================

if __name__ == "__main__":

    run_motif_pipeline("cdhit80", AMP_CDHIT, NONAMP50_CDHIT)
    run_motif_pipeline("mmseq80", AMP_MMSEQ, NONAMP50_MMSEQ)

    print("\n🎉 FULL 10x MOTIF PIPELINE COMPLETED SUCCESSFULLY.\n")