from pathlib import Path
import pandas as pd

# =====================================================
# PATH CONFIGURATION (RELATIVE TO REPO)
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_RAW_DIR = PROJECT_ROOT / "data" / "raw"
DATA_INTERMEDIATE_DIR = PROJECT_ROOT / "data" / "intermediate"

DATA_INTERMEDIATE_DIR.mkdir(parents=True, exist_ok=True)

DB_PATHS = {
    "CAMP": DATA_RAW_DIR / "CAMP" / "CAMP.txt",
    "DBAASP": DATA_RAW_DIR / "DBAASP" / "peptides.csv",
    "dbAMP3": DATA_RAW_DIR / "dbAMP3" / "dbAMP3_pepinfo.xlsx",
    "DRAMP": DATA_RAW_DIR / "DRAMP" / "natural_amps.txt",
}

# =====================================================
# STANDARD STRUCTURE
# =====================================================

STANDARD_COLUMNS = [
    "DB", "Complexity", "Protein_name", "Organism", "Taxonomy",
    "Sequence", "Activity", "Validation/Source", "Modifications",
    "Target Group", "Target Object", "Uniprot", "PBD", "Swissprot_entry"
]

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY-")

# =====================================================
# COLUMN MAPPINGS
# =====================================================

DRAMP_MAPPING = {
    "Sequence": "Sequence",
    "Name": "Protein_name",
    "Source(organism)": "Organism",
    "Activity": "Activity",
    "Swiss_Prot_Entry": "Swissprot_entry"
}

CAMP_MAPPING = {
    "Seqence": "Sequence",
    "Title": "Protein_name",
    "Source_Organism": "Organism",
    "Activity": "Activity",
    "Taxonomy": "Taxonomy",
    "Validation": "Validation/Source",
    "Modification": "Modifications",
    "Target": "Target Group",
}

DBAASP_MAPPING = {
    "COMPLEXITY": "Complexity",
    "NAME": "Protein_name",
    "SEQUENCE": "Sequence",
    "TARGET GROUP": "Target Group",
    "TARGET OBJECT": "Target Object",
    "SYNTHESIS TYPE": "Validation/Source",
}

dbAMP3_MAPPING = {
    "Seq": "Sequence",
    "Name": "Protein_name",
    "Tax": "Taxonomy",
    "Source": "Organism",
    "Uniprot": "Uniprot",
    "PDB": "PBD",
    "Targets": "Target Group",
}

MAPPINGS = {
    "CAMP": CAMP_MAPPING,
    "DBAASP": DBAASP_MAPPING,
    "dbAMP3": dbAMP3_MAPPING,
    "DRAMP": DRAMP_MAPPING
}

# =====================================================
# CLEANING
# =====================================================

def clean_sequences(df):

    df["Sequence"] = df["Sequence"].astype(str).str.upper().str.strip()

    def has_non_natural(seq):
        return any(aa not in VALID_AA for aa in seq)

    df["non_natural"] = df["Sequence"].apply(has_non_natural)
    df = df[~df["non_natural"]]

    df["Length"] = df["Sequence"].str.len()
    df = df[df["Length"] <= 80]

    return df.drop(columns=["non_natural"])


# =====================================================
# PROCESSING
# =====================================================

def process_database(db_name, input_path):

    if not input_path.exists():
        print(f"❌ File not found for {db_name}: {input_path}")
        return

    print(f"\nProcessing {db_name}...")

    # Read file depending on extension
    if input_path.suffix.lower() == ".csv":
        df = pd.read_csv(input_path)
    elif input_path.suffix.lower() in [".xlsx", ".xls"]:
        df = pd.read_excel(input_path)
    elif input_path.suffix.lower() == ".txt":
        df = pd.read_csv(input_path, sep="\t", encoding="latin-1", engine="python", on_bad_lines="skip")
    else:
        print(f"Unsupported format for {db_name}")
        return

    standardized_df = pd.DataFrame(columns=STANDARD_COLUMNS)
    mapping = MAPPINGS[db_name]

    for orig, std in mapping.items():
        if orig in df.columns:
            standardized_df[std] = df[orig]

    standardized_df["DB"] = db_name

    standardized_df = clean_sequences(standardized_df)

    output_file = DATA_INTERMEDIATE_DIR / f"{db_name}_standardized.xlsx"
    standardized_df.to_excel(output_file, index=False)

    print(f"Saved standardized file: {output_file}")


# =====================================================
# RUN
# =====================================================

if __name__ == "__main__":

    for db_name, path in DB_PATHS.items():
        process_database(db_name, path)

    print("\n✅ Step 01 completed successfully.\n")