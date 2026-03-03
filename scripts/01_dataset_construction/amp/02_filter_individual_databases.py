from pathlib import Path
import pandas as pd
import re

# =====================================================
# PATH CONFIG (RELATIVE TO REPO)
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_INTERMEDIATE_DIR = PROJECT_ROOT / "data" / "intermediate"

DB_NAMES = ["CAMP", "DBAASP", "dbAMP3", "DRAMP"]

MIN_LENGTH = 5

UNWANTED_NAME_PATTERNS = [
    "synthetic",
    "designed",
    "construct",
    "analog",
    "mutant",
    "fragment",
    "truncated"
]

# =====================================================
# FILTERING FUNCTION
# =====================================================

def filter_database(db_name):

    # 🔁 Ahora leemos los archivos generados por Script 01
    file_path = DATA_INTERMEDIATE_DIR / f"{db_name}_standardized.xlsx"

    if not file_path.exists():
        print(f"No {db_name}_standardized.xlsx found in {DATA_INTERMEDIATE_DIR}")
        return

    df = pd.read_excel(file_path)

    total_initial = len(df)

    # -------------------------------------------------
    # Remove sequences < 5 aa
    # -------------------------------------------------
    df["Length"] = df["Sequence"].astype(str).str.len()
    removed_short = (df["Length"] < MIN_LENGTH).sum()
    df = df[df["Length"] >= MIN_LENGTH]

    # -------------------------------------------------
    # Remove synthetic in Validation/Source
    # -------------------------------------------------
    if "Validation/Source" in df.columns:
        mask_synth_val = df["Validation/Source"].astype(str).str.contains(
            "synthetic", case=False, na=False
        )
        removed_synth_validation = mask_synth_val.sum()
        df = df[~mask_synth_val]
    else:
        removed_synth_validation = 0

    # -------------------------------------------------
    # Remove unwanted patterns in Protein_name
    # -------------------------------------------------
    if "Protein_name" in df.columns:
        pattern = "|".join(UNWANTED_NAME_PATTERNS)
        mask_unwanted_name = df["Protein_name"].astype(str).str.contains(
            pattern, case=False, na=False
        )
        removed_unwanted_name = mask_unwanted_name.sum()
        df = df[~mask_unwanted_name]
    else:
        removed_unwanted_name = 0

    # -------------------------------------------------
    # Remove 100% identical sequences
    # -------------------------------------------------
    before_dedup = len(df)
    df = df.drop_duplicates(subset=["Sequence"])
    removed_duplicates = before_dedup - len(df)

    total_final = len(df)

    # -------------------------------------------------
    # Save filtered file (not final yet — MIC comes next)
    # -------------------------------------------------
    output_file = DATA_INTERMEDIATE_DIR / f"{db_name}_structural_filtered.xlsx"
    df.drop(columns=["Length"], errors="ignore").to_excel(output_file, index=False)

    # -------------------------------------------------
    # Print stats
    # -------------------------------------------------
    print("\n========================================")
    print(f"STRUCTURAL FILTERING STATS - {db_name}")
    print("========================================")
    print("Initial:", total_initial)
    print("Removed short (<5 aa):", removed_short)
    print("Removed synthetic (Validation):", removed_synth_validation)
    print("Removed unwanted name patterns:", removed_unwanted_name)
    print("Removed duplicates (100% identical):", removed_duplicates)
    print("Final:", total_final)
    print("Saved:", output_file)
    print("========================================\n")


# =====================================================
# RUN ALL DATABASES
# =====================================================

if __name__ == "__main__":

    for db_name in DB_NAMES:
        filter_database(db_name)

    print("Step 02 (structural filtering) completed successfully.")