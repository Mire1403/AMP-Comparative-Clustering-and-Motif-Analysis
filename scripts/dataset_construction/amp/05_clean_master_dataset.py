from pathlib import Path
import pandas as pd
import re

# =====================================================
# CONFIG
# =====================================================

BASE_DIR = Path("/mnt/c/Users/mirei/Desktop/PAES/DBs/3")
MASTER_FILE = BASE_DIR / "DB_MASTER.xlsx"
OUTPUT_FILE = BASE_DIR / "DB_MASTER_CLEAN.xlsx"

# =====================================================
# NORMALIZATION FUNCTIONS
# =====================================================

def normalize_organism(text):
    if pd.isna(text):
        return ""

    text = str(text)

    # eliminar contenido entre paréntesis
    text = re.sub(r"\(.*?\)", "", text)

    # reemplazar &&
    text = text.replace("&&", ";")

    # quitar espacios extra
    text = re.sub(r"\s+", " ", text).strip()

    return text.title()


def normalize_taxonomy(text):
    if pd.isna(text):
        return ""

    text = str(text)

    # reemplazos comunes
    text = text.replace("&&", ";")
    text = text.replace(",", ";")

    # quitar espacios
    text = re.sub(r"\s+", "", text)

    parts = text.split(";")
    parts = list(set([p for p in parts if p != ""]))

    return ";".join(sorted(parts))


def concat_unique(series):
    values = series.dropna().astype(str)
    values = [v for v in values if v.strip() != ""]
    if not values:
        return ""
    return " | ".join(sorted(set(values)))

# =====================================================
# MAIN CLEANING
# =====================================================

def clean_master():

    if not MASTER_FILE.exists():
        print("DB_MASTER.xlsx not found")
        return

    df = pd.read_excel(MASTER_FILE)

    print("Initial sequences:", len(df))

    # -----------------------------------------
    # Normalizar organism y taxonomy
    # -----------------------------------------

    if "Organism" in df.columns:
        df["Organism"] = df["Organism"].apply(normalize_organism)

    if "Taxonomy" in df.columns:
        df["Taxonomy"] = df["Taxonomy"].apply(normalize_taxonomy)

    # -----------------------------------------
    # Agrupar por secuencia
    # -----------------------------------------

    grouped = df.groupby("Sequence", as_index=False).agg(concat_unique)

    print("After removing inter-database duplicates:", len(grouped))

    # -----------------------------------------
    # Estadísticas
    # -----------------------------------------

    print("\nBreakdown by Source_DB:")
    if "Source_DB" in grouped.columns:
        print(grouped["Source_DB"].value_counts())

    # -----------------------------------------
    # Guardar resultado
    # -----------------------------------------

    grouped.to_excel(OUTPUT_FILE, index=False)

    print("\nSaved cleaned master at:")
    print(OUTPUT_FILE)


# =====================================================
# RUN
# =====================================================

if __name__ == "__main__":
    clean_master()