from pathlib import Path
import numpy as np

# =====================================================
# PATH CONFIG (RELATIVE TO REPO)
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[3]
DATA_RAW_DIR = PROJECT_ROOT / "data" / "raw" / "non_AMP"
DATA_INTERMEDIATE_DIR = PROJECT_ROOT / "data" / "intermediate"

DATA_INTERMEDIATE_DIR.mkdir(parents=True, exist_ok=True)

INPUT_FASTA = DATA_RAW_DIR / "uniprot_nonamp_raw.fasta"
OUTPUT_FASTA = DATA_INTERMEDIATE_DIR / "nonamp_clean_min5.fasta"

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

# =====================================================
# FUNCTIONS
# =====================================================

def read_fasta(fasta_path):

    sequences = []
    header = None
    seq = ""

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq:
                    sequences.append((header, seq))
                header = line
                seq = ""
            else:
                seq += line.upper()

        if seq:
            sequences.append((header, seq))

    return sequences


def compute_stats(sequences):

    lengths = [len(seq) for _, seq in sequences]

    if not lengths:
        return {"count": 0, "min": 0, "max": 0, "mean": 0, "std": 0}

    return {
        "count": len(lengths),
        "min": min(lengths),
        "max": max(lengths),
        "mean": round(np.mean(lengths), 2),
        "std": round(np.std(lengths), 2),
    }


def is_valid_sequence(seq):
    return all(aa in VALID_AA for aa in seq)

# =====================================================
# MAIN
# =====================================================

if __name__ == "__main__":

    if not INPUT_FASTA.exists():
        print(f"Input FASTA not found: {INPUT_FASTA}")
        exit()

    sequences = read_fasta(INPUT_FASTA)

    print("\n==============================")
    print("INITIAL STATS")
    print("==============================")

    initial_stats = compute_stats(sequences)
    for k, v in initial_stats.items():
        print(f"{k}: {v}")

    # -------------------------------------------------
    # FILTERING
    # -------------------------------------------------

    filtered_sequences = []

    removed_short = 0
    removed_non_natural = 0

    for header, seq in sequences:

        length = len(seq)

        if length <= 5:
            removed_short += 1
            continue

        if not is_valid_sequence(seq):
            removed_non_natural += 1
            continue

        filtered_sequences.append((header, seq))

    # -------------------------------------------------
    # FINAL STATS
    # -------------------------------------------------

    print("\n==============================")
    print("FILTERING SUMMARY")
    print("==============================")
    print("Removed length ≤ 5:", removed_short)
    print("Removed non-natural AA:", removed_non_natural)
    print("Remaining:", len(filtered_sequences))

    final_stats = compute_stats(filtered_sequences)

    print("\n==============================")
    print("FINAL STATS")
    print("==============================")
    for k, v in final_stats.items():
        print(f"{k}: {v}")

    # -------------------------------------------------
    # SAVE CLEAN FASTA
    # -------------------------------------------------

    with open(OUTPUT_FASTA, "w") as f:
        for header, seq in filtered_sequences:
            f.write(f"{header}\n{seq}\n")

    print("\nClean FASTA saved at:")
    print(OUTPUT_FASTA)
    print("\nStep 01 (non-AMP cleaning) completed successfully.\n")