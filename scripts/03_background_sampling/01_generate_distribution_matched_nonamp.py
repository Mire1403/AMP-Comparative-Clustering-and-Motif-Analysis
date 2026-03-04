from __future__ import annotations

from pathlib import Path
import sys
import random
from collections import defaultdict
from datetime import datetime

import pandas as pd
import numpy as np
from scipy.stats import ks_2samp

import matplotlib
matplotlib.use("Agg")  # ✅ avoid Qt/Wayland issues in WSL
import matplotlib.pyplot as plt

print("\nBACKGROUND SAMPLING PIPELINE — progressive random (len + KR-matched)\n")

# =====================================================
# PATH CONFIG
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

RESULTS_CLUSTER_DIR = PROJECT_ROOT / "results" / "clustering"

# Big outputs (FASTA)
RESULTS_BG_DIR = PROJECT_ROOT / "results" / "background_generation"
RESULTS_BG_DIR.mkdir(parents=True, exist_ok=True)

# Lightweight outputs (plots, reports, manifests, csv summary)
STATS_DIR = PROJECT_ROOT / "results" / "statistics" / "background_sampling"
STATS_DIR.mkdir(parents=True, exist_ok=True)

# Inputs produced by clustering scripts
AMP_FASTA_CDHIT = RESULTS_CLUSTER_DIR / "AMP_MASTER_cdhit.fasta"
AMP_FASTA_MMSEQ = RESULTS_CLUSTER_DIR / "AMP_MASTER_mmseq_rep_seq.fasta"

NONAMP_FASTA_CDHIT = RESULTS_CLUSTER_DIR / "NONAMP_cdhit.fasta"
NONAMP_FASTA_MMSEQ = RESULTS_CLUSTER_DIR / "NONAMP_mmseq_rep_seq.fasta"

# =====================================================
# CONFIG
# =====================================================

MULTIPLIER = 10
MIN_LENGTH = 5
RANDOM_SEED = 42
ALPHA = 0.05

random.seed(RANDOM_SEED)

SUMMARY_CSV = STATS_DIR / "background_sampling_summary.csv"

# =====================================================
# HELPERS
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

    initial_count = len(sequences)
    sequences = [s for s in sequences if len(s) >= MIN_LENGTH]
    removed = initial_count - len(sequences)
    print(f"{path.name} → removed (<{MIN_LENGTH} aa): {removed}")

    return sequences

def compute_amp_stats(sequences: list[str]):
    """Return length histogram + per-length (mean, std) of KR proportion."""
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

def write_manifest(outdir: Path, label: str, inputs: dict[str, Path]) -> Path:
    manifest = outdir / f"{label}_background_manifest.txt"
    with open(manifest, "w", encoding="utf-8") as f:
        f.write("BACKGROUND SAMPLING MANIFEST\n")
        f.write(f"timestamp_utc={datetime.utcnow().isoformat()}Z\n")
        f.write(f"label={label}\n")
        f.write(f"MULTIPLIER={MULTIPLIER}\n")
        f.write(f"MIN_LENGTH={MIN_LENGTH}\n")
        f.write(f"RANDOM_SEED={RANDOM_SEED}\n")
        f.write(f"ALPHA={ALPHA}\n")
        for k, v in inputs.items():
            f.write(f"{k}={v}\n")
    return manifest

# =====================================================
# PROGRESSIVE RANDOM SAMPLING (SLIDING WINDOW)
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

    if len(selected_fragments) < total_required:
        print(
            f"⚠️ Generated {len(selected_fragments)}/{total_required} fragments "
            f"({total_required - len(selected_fragments)} missing). "
            f"Tip: lower MULTIPLIER or widen KR window."
        )

    return selected_fragments, total_required

# =====================================================
# VALIDATION
# =====================================================

def validate_and_plot(amp_seqs, nonamp_seqs, label: str):
    amp_lengths = [len(s) for s in amp_seqs]
    nonamp_lengths = [len(s) for s in nonamp_seqs]

    amp_kr = [percent_KR(s) for s in amp_seqs]
    nonamp_kr = [percent_KR(s) for s in nonamp_seqs]

    ks_len = ks_2samp(amp_lengths, nonamp_lengths)
    ks_kr = ks_2samp(amp_kr, nonamp_kr)

    def interpret(p):
        return "NOT significantly different" if p >= ALPHA else "Significantly different"

    results_file = STATS_DIR / f"{label}_validation_report.txt"
    with open(results_file, "w", encoding="utf-8") as f:
        f.write("KOLMOGOROV–SMIRNOV TEST RESULTS\n\n")
        f.write("Length Distribution:\n")
        f.write(f"KS statistic: {ks_len.statistic}\n")
        f.write(f"p-value: {ks_len.pvalue}\n")
        f.write(f"Interpretation (alpha={ALPHA}): {interpret(ks_len.pvalue)}\n\n")
        f.write("K+R Composition Distribution:\n")
        f.write(f"KS statistic: {ks_kr.statistic}\n")
        f.write(f"p-value: {ks_kr.pvalue}\n")
        f.write(f"Interpretation (alpha={ALPHA}): {interpret(ks_kr.pvalue)}\n")

    # Plots (lightweight)
    plt.figure()
    plt.hist(amp_lengths, bins=30, density=True, alpha=0.5, label="AMP")
    plt.hist(nonamp_lengths, bins=30, density=True, alpha=0.5, label="NONAMP")
    plt.legend()
    plt.xlabel("Peptide Length")
    plt.ylabel("Density")
    plt.title(f"{label.upper()} — Length (AMP vs NONAMP)")
    plt.tight_layout()
    plt.savefig(STATS_DIR / f"{label}_length_superposed.png")
    plt.close()

    plt.figure()
    plt.hist(amp_kr, bins=30, density=True, alpha=0.5, label="AMP")
    plt.hist(nonamp_kr, bins=30, density=True, alpha=0.5, label="NONAMP")
    plt.legend()
    plt.xlabel("K+R Proportion")
    plt.ylabel("Density")
    plt.title(f"{label.upper()} — KR (AMP vs NONAMP)")
    plt.tight_layout()
    plt.savefig(STATS_DIR / f"{label}_kr_superposed.png")
    plt.close()

    return {
        "ks_len_stat": float(ks_len.statistic),
        "ks_len_p": float(ks_len.pvalue),
        "ks_kr_stat": float(ks_kr.statistic),
        "ks_kr_p": float(ks_kr.pvalue),
        "len_match_ok": bool(ks_len.pvalue >= ALPHA),
        "kr_match_ok": bool(ks_kr.pvalue >= ALPHA),
    }

# =====================================================
# RUN
# =====================================================

def run_one(label: str, amp_fasta: Path, nonamp_fasta: Path) -> dict:
    inputs = {"AMP_FASTA": amp_fasta, "NONAMP_FASTA": nonamp_fasta}
    manifest = write_manifest(STATS_DIR, label, inputs)
    print("Manifest saved:", manifest.name)

    amp = read_fasta(amp_fasta)
    nonamp = read_fasta(nonamp_fasta)

    length_counts, kr_stats = compute_amp_stats(amp)
    target_counts = {L: count * MULTIPLIER for L, count in length_counts.items()}

    generated, expected = progressive_random_sampling(nonamp, target_counts, kr_stats)

    # BIG OUTPUT (do not track on GitHub)
    out_fasta = RESULTS_BG_DIR / f"nonamp_{label}_{MULTIPLIER}x_len_kr_matched_progressive_random.fasta"
    with open(out_fasta, "w", encoding="utf-8") as f:
        for i, frag in enumerate(generated, start=1):
            f.write(f">NONAMP_{label}_frag_{i:08d}\n{frag}\n")

    print(f"{label.upper()} background: {len(generated)}/{expected} fragments")
    val = validate_and_plot(amp, generated, label)

    return {
        "label": label,
        "multiplier": MULTIPLIER,
        "expected": int(expected),
        "generated": int(len(generated)),
        "coverage": float(len(generated) / expected) if expected else np.nan,
        "background_fasta": str(out_fasta),
        **val,
    }

def main() -> None:
    rows = []

    print("\n=== CDHIT ===")
    rows.append(run_one("cdhit", AMP_FASTA_CDHIT, NONAMP_FASTA_CDHIT))

    print("\n=== MMSEQ ===")
    rows.append(run_one("mmseq", AMP_FASTA_MMSEQ, NONAMP_FASTA_MMSEQ))

    pd.DataFrame(rows).to_csv(SUMMARY_CSV, index=False)

    print("\n✅ PIPELINE COMPLETE — BACKGROUND SAMPLING (len + KR matched).")
    print("Summary:", SUMMARY_CSV)

if __name__ == "__main__":
    main()