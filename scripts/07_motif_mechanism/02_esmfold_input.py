import pandas as pd
import os

# ==============================
# CONFIGURATION
# ==============================

INPUT_TABLE = "results/statistics/08_fmap_input/motif_selected_sequences.csv"
OUTPUT_FASTA = "results/statistics/09_esmfold_input/esm_input.fasta"

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

# ==============================
# HELPERS
# ==============================

def clean_sequence(seq):
    """
    Replace non-standard amino acids with 'L' to avoid ESMFold errors.
    """
    return "".join([aa if aa in VALID_AA else "L" for aa in seq])

# ==============================
# MAIN
# ==============================

def main():

    print("Loading sequence table...")
    df = pd.read_csv(INPUT_TABLE)

    os.makedirs(os.path.dirname(OUTPUT_FASTA), exist_ok=True)

    print("Writing FASTA file for ESMFold...")

    with open(OUTPUT_FASTA, "w") as f:
        for _, row in df.iterrows():

            peptide_id = row["Peptide_ID"]
            sequence = clean_sequence(row["Sequence"])

            f.write(f">{peptide_id}\n{sequence}\n\n")

    print(f"\nFASTA file created: {OUTPUT_FASTA}")
    print(f"Total sequences: {len(df)}")

# ==============================
# RUN
# ==============================

if __name__ == "__main__":
    main()