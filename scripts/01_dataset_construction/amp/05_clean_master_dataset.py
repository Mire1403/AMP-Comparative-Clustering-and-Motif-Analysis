from __future__ import annotations

from pathlib import Path
import pandas as pd
import re


# =====================================================
# PATH CONFIG (RELATIVE TO REPO)
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_INTERMEDIATE_DIR = PROJECT_ROOT / "data" / "intermediate"
DATA_FINAL_DIR = PROJECT_ROOT / "data" / "final"
DATA_FINAL_DIR.mkdir(parents=True, exist_ok=True)

MASTER_FILE = DATA_INTERMEDIATE_DIR / "DB_MASTER.xlsx"
OUTPUT_FILE = DATA_FINAL_DIR / "DB_MASTER_CLEAN.xlsx"


# =====================================================
# NORMALIZATION FUNCTIONS
# =====================================================

def normalize_organism(text) -> str:
    if pd.isna(text):
        return ""

    text = str(text)

    # remove content in parentheses
    text = re.sub(r"\(.*?\)", "", text)

    # replace &&
    text = text.replace("&&", ";")

    # collapse whitespace
    text = re.sub(r"\s+", " ", text).strip()

    return text.title()


def normalize_taxonomy(text) -> str:
    if pd.isna(text):
        return ""

    text = str(text)

    # common replacements
    text = text.replace("&&", ";")
    text = text.replace(",", ";")

    # remove spaces
    text = re.sub(r"\s+", "", text)

    parts = [p for p in text.split(";") if p]
    parts = sorted(set(parts))

    return ";".join(parts)


def concat_unique(series: pd.Series) -> str:
    values = series.dropna().astype(str)
    values = [v.strip() for v in values if v.strip() != ""]
    if not values:
        return ""
    return " | ".join(sorted(set(values)))


# =====================================================
# MAIN
# =====================================================

def clean_master() -> None:
    if not MASTER_FILE.exists():
        print(f"❌ Not found: {MASTER_FILE}")
        return

    df = pd.read_excel(MASTER_FILE)
    print("Initial sequences:", len(df))

    # Normalize fields if present
    if "Organism" in df.columns:
        df["Organism"] = df["Organism"].apply(normalize_organism)

    if "Taxonomy" in df.columns:
        df["Taxonomy"] = df["Taxonomy"].apply(normalize_taxonomy)

    # Group by sequence: merge text columns by unique values
    grouped = df.groupby("Sequence", as_index=False).agg(concat_unique)

    print("After removing inter-database duplicates:", len(grouped))

    # Stats
    if "Source_DB" in grouped.columns:
        print("\nBreakdown by Source_DB:")
        print(grouped["Source_DB"].value_counts())

    grouped.to_excel(OUTPUT_FILE, index=False)

    print("\n✅ Saved cleaned master at:")
    print(OUTPUT_FILE)


if __name__ == "__main__":
    clean_master()