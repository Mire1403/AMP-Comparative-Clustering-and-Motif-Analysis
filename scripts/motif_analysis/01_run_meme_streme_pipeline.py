from __future__ import annotations

from pathlib import Path
import random
from collections import defaultdict
import subprocess
import sys

import numpy as np
from scipy.stats import ks_2samp

print("\nSMART 10x MOTIF PIPELINE\n")

# =====================================================
# PATH CONFIG 
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

RESULTS_CLUSTER_DIR = PROJECT_ROOT / "results" / "clustering"
RESULTS_BACKGROUND_DIR = PROJECT_ROOT / "results" / "background_generation"
RESULTS_MOTIF_DIR = PROJECT_ROOT / "results" / "motif_discovery"
RESULTS_MOTIF_DIR.mkdir(parents=True, exist_ok=True)

# Inputs produced by your UPDATED clustering/background scripts
AMP_CDHIT = RESULTS_CLUSTER_DIR / "AMP_MASTER_cdhit.fasta"
AMP_MMSEQ = RESULTS_CLUSTER_DIR / "AMP_MASTER_mmseq_rep_seq.fasta"

NONAMP50_CDHIT = RESULTS_BACKGROUND_DIR / "nonamp_cdhit_50x_random_progressive.fasta"
NONAMP50_MMSEQ = RESULTS_BACKGROUND_DIR / "nonamp_mmseq_50x_random_progressive.fasta"

# =====================================================
# CONFIG
# =====================================================

MULTIPLIER = 10
ALPHA = 0.05
RANDOM_SEED = 42
random.seed(RANDOM_SEED)

# MEME/STREME parameters (same logic as you had)
MINW = "6"
MAXW = "20"
STREME_NMOTIFS = "15"
STREME_PVT = "0.01"
MEME_NMOTIFS = "10"
MEME_MAXSIZE = "10000000"

# =====================================================
# UTILITIES
# =====================================================

def die(msg: str) -> None:
    print(f"{msg}")
    sys.exit(1)

def check_file(path: Path) -> None:
    if not path.exists():
        die(f"File not found: {path}")

def read_fasta(path: Path) -> list[str]:
    check_file(path)
    seqs: list[str] = []
    seq_parts: list[str] = []

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_parts:
                    seqs.append("".join(seq_parts).upper())
                seq_parts = []
            else:
                seq_parts.append(line)

        if seq_parts:
            seqs.append("".join(seq_parts).upper())

    return seqs

def get_length_distribution(seqs: list[str]) -> dict[int, int]:
    dist = defaultdict(int)
    for s in seqs:
        dist[len(s)] += 1
    return dict(dist)

def group_by_length(seqs: list[str]) -> dict[int, list[str]]:
    groups = defaultdict(list)
    for s in seqs:
        groups[len(s)].append(s)
    return dict(groups)

def run_command(cmd: list[str]) -> None:
    print("\nRunning:", " ".join(str(c) for c in cmd))
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError:
        die("Command failed. Stopping pipeline.")

# =====================================================
# CORE GENERATOR
# =====================================================

def generate_10x(amp_fasta: Path, nonamp50_fasta: Path, label: str) -> Path:
    """
    Creates a 10x non-AMP FASTA by sampling from a pre-generated 50x background,
    matching length distribution to AMP (x10).
    """
    print(f"\n=== GENERATING {MULTIPLIER}x FOR {label.upper()} ===")

    amp_seqs = read_fasta(amp_fasta)
    nonamp50_seqs = read_fasta(nonamp50_fasta)

    amp_dist = get_length_distribution(amp_seqs)
    nonamp_groups = group_by_length(nonamp50_seqs)

    target_counts = {L: amp_dist[L] * MULTIPLIER for L in amp_dist}

    selected: list[str] = []

    for L, required in sorted(target_counts.items()):
        available = nonamp_groups.get(L, [])
        if len(available) < required:
            die(f"Not enough nonAMP sequences for length {L}: have {len(available)}, need {required}")

        sampled = random.sample(available, required)
        selected.extend(sampled)

    print("Expected:", sum(target_counts.values()))
    print("Generated:", len(selected))

    # VALIDATION: KS test on lengths
    amp_lengths = [len(s) for s in amp_seqs]
    nonamp_lengths = [len(s) for s in selected]
    ks_len = ks_2samp(amp_lengths, nonamp_lengths)

    print("KS p-value (length):", ks_len.pvalue)

    if ks_len.pvalue < ALPHA:
        die("Length distributions significantly different (KS test).")

    print(" Length distributions validated.")

    output_fasta = RESULTS_MOTIF_DIR / f"nonamp_{label}_{MULTIPLIER}x.fasta"
    with open(output_fasta, "w", encoding="utf-8") as f:
        for i, s in enumerate(selected, start=1):
            f.write(f">NONAMP_{label}_{i}\n{s}\n")

    print("Saved:", output_fasta)
    return output_fasta

# =====================================================
# MOTIF DISCOVERY (STREME + MEME)
# =====================================================

def run_motif_pipeline(label: str, amp_fasta: Path, nonamp50_fasta: Path) -> None:
    nonamp10 = generate_10x(amp_fasta, nonamp50_fasta, label)

    # STREME
    print(f"\n Running STREME ({label})")
    streme_out = RESULTS_MOTIF_DIR / f"streme_{label}"
    run_command([
        "streme",
        "--p", str(amp_fasta),
        "--n", str(nonamp10),
        "--protein",
        "--minw", MINW,
        "--maxw", MAXW,
        "--nmotifs", STREME_NMOTIFS,
        "--pvt", STREME_PVT,
        "--oc", str(streme_out),
    ])
    print(" STREME completed:", streme_out)

    # MEME
    print(f"\n Running MEME ({label})")
    meme_out = RESULTS_MOTIF_DIR / f"meme_{label}"
    run_command([
        "meme",
        str(amp_fasta),
        "-neg", str(nonamp10),
        "-mod", "zoops",
        "-objfun", "de",
        "-nmotifs", MEME_NMOTIFS,
        "-minw", MINW,
        "-maxw", MAXW,
        "-maxsize", MEME_MAXSIZE,
        "-oc", str(meme_out),
    ])
    print(" MEME completed:", meme_out)

# =====================================================
# MAIN
# =====================================================

def main():
    check_file(AMP_CDHIT)
    check_file(AMP_MMSEQ)
    check_file(NONAMP50_CDHIT)
    check_file(NONAMP50_MMSEQ)

    run_motif_pipeline("cdhit", AMP_CDHIT, NONAMP50_CDHIT)
    run_motif_pipeline("mmseq", AMP_MMSEQ, NONAMP50_MMSEQ)

    print("\n FULL 10x MOTIF PIPELINE COMPLETED SUCCESSFULLY.\n")

if __name__ == "__main__":
    main()