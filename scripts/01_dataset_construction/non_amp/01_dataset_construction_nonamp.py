from __future__ import annotations

from pathlib import Path
import sys
import numpy as np
import pandas as pd

# =====================================================
# PATH CONFIG (RELATIVE TO REPO)
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[3]
DATA_RAW_DIR = PROJECT_ROOT / "data" / "raw" / "non_AMPs"
DATA_INTERMEDIATE_DIR = PROJECT_ROOT / "data" / "intermediate"

DATA_INTERMEDIATE_DIR.mkdir(parents=True, exist_ok=True)

INPUT_FASTA = DATA_RAW_DIR / "uniprot_nonamp_raw.fasta"
OUTPUT_FASTA = DATA_INTERMEDIATE_DIR / "nonamp_clean_min5.fasta"
REPORT_CSV = DATA_INTERMEDIATE_DIR / "nonamp_clean_min5_report.csv"

MIN_LEN = 6
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

# =====================================================
# FASTA HELPERS
# =====================================================

def read_fasta(path: Path) -> list[tuple[str, str]]:
    records = []
    header = None
    seq_parts = []

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts).upper()))
                header = line
                seq_parts = []
            else:
                seq_parts.append(line)

        if header is not None:
            records.append((header, "".join(seq_parts).upper()))

    return records


def write_fasta(records: list[tuple[str, str]], out_path: Path, wrap: int = 80) -> None:
    with open(out_path, "w", encoding="utf-8") as f:
        for header, seq in records:
            f.write(f"{header}\n")
            for i in range(0, len(seq), wrap):
                f.write(seq[i:i+wrap] + "\n")


def compute_stats(records: list[tuple[str, str]]) -> dict:
    lengths = [len(seq) for _, seq in records]

    if not lengths:
        return {"count": 0, "min": 0, "max": 0, "mean": 0.0, "std": 0.0}

    return {
        "count": len(lengths),
        "min": int(min(lengths)),
        "max": int(max(lengths)),
        "mean": float(round(np.mean(lengths), 2)),
        "std": float(round(np.std(lengths), 2)),
    }


def is_valid_sequence(seq: str) -> bool:
    return bool(seq) and all(aa in VALID_AA for aa in seq)

# =====================================================
# MAIN
# =====================================================

def main() -> None:

    if not INPUT_FASTA.exists():
        print(f"❌ Input FASTA not found: {INPUT_FASTA}")
        sys.exit(1)

    records = read_fasta(INPUT_FASTA)

    initial_stats = compute_stats(records)

    print("\n===== INITIAL STATS =====")
    for k, v in initial_stats.items():
        print(f"{k}: {v}")

    filtered = []
    removed_short = 0
    removed_non_natural = 0

    seen_sequences = set()

    for header, seq in records:

        if len(seq) < MIN_LEN:
            removed_short += 1
            continue

        if not is_valid_sequence(seq):
            removed_non_natural += 1
            continue

        if seq in seen_sequences:
            continue

        seen_sequences.add(seq)
        filtered.append((header, seq))

    final_stats = compute_stats(filtered)

    print("\n===== FILTER SUMMARY =====")
    print(f"Removed length < {MIN_LEN}: {removed_short}")
    print(f"Removed non-natural AA: {removed_non_natural}")
    print(f"Remaining unique sequences: {len(filtered)}")

    print("\n===== FINAL STATS =====")
    for k, v in final_stats.items():
        print(f"{k}: {v}")

    write_fasta(filtered, OUTPUT_FASTA, wrap=80)

    report = {
        "input_fasta": str(INPUT_FASTA),
        "output_fasta": str(OUTPUT_FASTA),
        "min_len": MIN_LEN,
        "initial_count": initial_stats["count"],
        "final_count_unique": final_stats["count"],
        "removed_short": removed_short,
        "removed_non_natural": removed_non_natural,
        "initial_mean_len": initial_stats["mean"],
        "final_mean_len": final_stats["mean"],
        "initial_min_len": initial_stats["min"],
        "final_min_len": final_stats["min"],
        "initial_max_len": initial_stats["max"],
        "final_max_len": final_stats["max"],
    }

    pd.DataFrame([report]).to_csv(REPORT_CSV, index=False)

    print("\n✅ Clean FASTA saved at:")
    print(OUTPUT_FASTA)
    print("🧾 Report saved at:")
    print(REPORT_CSV)
    print("\n✅ Non-AMP cleaning completed successfully.\n")


if __name__ == "__main__":
    main()