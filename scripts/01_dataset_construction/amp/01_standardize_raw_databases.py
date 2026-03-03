from __future__ import annotations

from pathlib import Path
import pandas as pd

from config import DB_PATHS, DATA_INTERMEDIATE_DIR

# =====================================================
# STANDARD STRUCTURE
# =====================================================

STANDARD_COLUMNS = [
    "source_db",
    "complexity",
    "protein_name",
    "organism",
    "taxonomy",
    "sequence",
    "activity",
    "validation_source",
    "modifications",
    "target_group",
    "target_object",
    "uniprot",
    "pdb",
    "swissprot_entry",
]

MAPPINGS = {
    "DRAMP": {
        "Sequence": "sequence",
        "Name": "protein_name",
        "Source": "organism",  # (si tu DRAMP real es Source(organism), cámbialo aquí)
        "Activity": "activity",
        "Swiss_Prot_Entry": "swissprot_entry",
    },
    "CAMP": {
        "Seqence": "sequence",  # typo típico en CAMP
        "Title": "protein_name",
        "Source_Organism": "organism",
        "Activity": "activity",
        "Taxonomy": "taxonomy",
        "Validation": "validation_source",
        "Modifications": "modifications",
        "Target": "target_group",
    },
    "DBAASP": {
        "COMPLEXITY": "complexity",
        "NAME": "protein_name",
        "SEQUENCE": "sequence",
        "TARGET GROUP": "target_group",
        "TARGET OBJECT": "target_object",
        "SYNTHESIS TYPE": "validation_source",
    },
    "dbAMP3": {
        "Seq": "sequence",
        "Name": "protein_name",
        "Tax": "taxonomy",
        "Source": "organism",
        "Uniprot": "uniprot",
        "PDB": "pdb",
        "Targets": "target_group",
    },
}

TEXT_COLS = [
    "source_db",
    "protein_name",
    "organism",
    "taxonomy",
    "sequence",
    "activity",
    "validation_source",
    "modifications",
    "target_group",
    "target_object",
    "uniprot",
    "pdb",
    "swissprot_entry",
    "complexity",
]

# =====================================================
# IO
# =====================================================

def _read_any(path: Path) -> pd.DataFrame:
    suf = path.suffix.lower()
    if suf == ".csv":
        return pd.read_csv(path)
    if suf in {".xlsx", ".xls"}:
        # requires openpyxl for .xlsx
        return pd.read_excel(path)
    if suf == ".txt":
        return pd.read_csv(
            path,
            sep="\t",
            encoding="latin-1",
            engine="python",
            on_bad_lines="skip",
        )
    raise ValueError(f"Unsupported file format: {path}")

def _normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Make column matching robust across weird exports."""
    cleaned = []
    for c in df.columns:
        s = str(c).strip()
        s = s.strip('"').strip("'")    # remove surrounding quotes
        s = s.replace("’", "'")        # normalize fancy apostrophes
        s = s.strip().strip("'")       # remove trailing apostrophes
        cleaned.append(s)
    df = df.copy()
    df.columns = cleaned
    return df

def _force_text_types(df: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure parquet-safe types:
    Arrow hates object columns mixing int/str. Convert to pandas 'string'.
    """
    df = df.copy()
    for c in TEXT_COLS:
        if c in df.columns:
            df[c] = df[c].astype("string")
    return df

def _write_table(df: pd.DataFrame, out_parquet: Path) -> Path:
    """
    Prefer parquet. If it fails (engine or types), fallback to CSV.
    Returns the file actually written.
    """
    try:
        df.to_parquet(out_parquet, index=False)
        return out_parquet
    except Exception as e:
        fallback = out_parquet.with_suffix(".csv")
        print(f"⚠️  Could not write parquet ({type(e).__name__}: {e}). Falling back to CSV: {fallback.name}")
        df.to_csv(fallback, index=False)
        return fallback

# =====================================================
# STANDARDIZE
# =====================================================

def standardize_one(db_name: str, input_path: Path) -> Path:
    if db_name not in MAPPINGS:
        raise KeyError(f"No mapping defined for db_name='{db_name}'")

    if not input_path.exists():
        raise FileNotFoundError(f"Missing {db_name} input: {input_path}")

    print(f"\n=== {db_name} ===")
    print(f"Reading: {input_path}")

    df = _read_any(input_path)
    df = _normalize_columns(df)
    print(f"Raw rows: {len(df)}")

    out = pd.DataFrame(columns=STANDARD_COLUMNS)
    mapping = MAPPINGS[db_name]

    missing = []
    for orig, std in mapping.items():
        if orig in df.columns:
            out[std] = df[orig]
        else:
            missing.append(orig)

    if missing:
        print(f"⚠️  Missing columns in {db_name} input (skipped): {missing}")

    out["source_db"] = db_name

    # Ensure sequence exists
    if "sequence" not in out.columns:
        out["sequence"] = pd.NA

    # Normalize sequence
    out["sequence"] = (
        out["sequence"]
        .astype("string")
        .str.upper()
        .str.strip()
    )

    # Drop empty sequences
    before = len(out)
    out = out[out["sequence"].notna() & (out["sequence"] != "")]
    dropped = before - len(out)
    if dropped:
        print(f"Dropped empty/NA sequences: {dropped}")

    # Dedup within DB
    before = len(out)
    out = out.drop_duplicates(subset=["sequence"])
    print(f"Dedup within {db_name}: removed {before - len(out)}")

    # Force safe types BEFORE saving (fixes dbAMP3 parquet error)
    out = _force_text_types(out)

    DATA_INTERMEDIATE_DIR.mkdir(parents=True, exist_ok=True)
    out_parquet = DATA_INTERMEDIATE_DIR / f"{db_name}_01_standardized.parquet"

    written = _write_table(out, out_parquet)

    print(f"Saved: {written}")
    print(f"Final rows: {len(out)}")
    return written

def main() -> None:
    for db_name, path in DB_PATHS.items():
        standardize_one(db_name, path)

    print("\n Step 01 completed successfully.\n")

if __name__ == "__main__":
    main()