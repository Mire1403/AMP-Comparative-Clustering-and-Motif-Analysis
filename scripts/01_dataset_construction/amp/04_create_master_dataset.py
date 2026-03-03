from pathlib import Path
import pandas as pd

# =====================================================
# PATH CONFIG (RELATIVE TO REPO)
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_INTERMEDIATE_DIR = PROJECT_ROOT / "data" / "intermediate"

DB_NAMES = ["CAMP", "DBAASP", "dbAMP3", "DRAMP"]

OUTPUT_FILE = DATA_INTERMEDIATE_DIR / "DB_MASTER.xlsx"

# =====================================================
# MAIN
# =====================================================

def create_master():

    all_dfs = []

    for db_name in DB_NAMES:

        file_path = DATA_INTERMEDIATE_DIR / f"{db_name}_activity_filtered.xlsx"

        if not file_path.exists():
            print(f"❌ No file found: {file_path.name}")
            continue

        print(f"Reading {file_path.name}")

        df = pd.read_excel(file_path)

        df["Source_DB"] = db_name  # trazabilidad

        all_dfs.append(df)

    if not all_dfs:
        print("No data loaded.")
        return

    master_df = pd.concat(all_dfs, ignore_index=True)

    print("\n========================================")
    print("DB MASTER CREATED")
    print("========================================")
    print("Total sequences:", len(master_df))
    print("Breakdown by DB:")
    print(master_df["Source_DB"].value_counts())
    print("========================================\n")

    master_df.to_excel(OUTPUT_FILE, index=False)

    print("Saved at:", OUTPUT_FILE)


# =====================================================
# RUN
# =====================================================

if __name__ == "__main__":
    create_master()