from __future__ import annotations
import re
from pathlib import Path
import pandas as pd

from config import (
    DATA_INTERMEDIATE_DIR,
    MIN_LENGTH,
    MAX_LENGTH,
    VALID_AA,
    UNWANTED_NAME_PATTERNS,
)

def is_natural(seq: str) -> bool:
    if not seq:
        return False
    return all(aa in VALID_AA for aa in seq)

def filter_one(db_name: str) -> Path:
    in_file = DATA_INTERMEDIATE_DIR / f"{db_name}_01_standardized.parquet"
    if not in_file.exists():
        raise FileNotFoundError(in_file)

    df = pd.read_parquet(in_file)
    n0 = len(df)

    # length
    df["length"] = df["sequence"].astype(str).str.len()
    df = df[(df["length"] >= MIN_LENGTH) & (df["length"] <= MAX_LENGTH)]
    n_len = len(df)

    # natural AA only
    df = df[df["sequence"].astype(str).apply(is_natural)]
    n_nat = len(df)

    # remove "synthetic" in validation_source
    if "validation_source" in df.columns:
        mask = df["validation_source"].astype(str).str.contains("synthetic", case=False, na=False)
        df = df[~mask]
    n_syn = len(df)

    # unwanted protein_name patterns
    if "protein_name" in df.columns:
        pattern = "|".join(map(re.escape, UNWANTED_NAME_PATTERNS))
        mask = df["protein_name"].astype(str).str.contains(pattern, case=False, na=False)
        df = df[~mask]
    n_name = len(df)

    # dedup within DB
    before = len(df)
    df = df.drop_duplicates(subset=["sequence"])
    n_dedup = len(df)

    out_file = DATA_INTERMEDIATE_DIR / f"{db_name}_02_structural_filtered.parquet"
    df.to_parquet(out_file, index=False)

    print(
        f"{db_name}: start={n0} -> len={n_len} -> natural={n_nat} -> no_synth={n_syn} -> name_ok={n_name} -> dedup={n_dedup}"
    )
    return out_file

def main() -> None:
    for db_name in ["CAMP", "DBAASP", "dbAMP3", "DRAMP"]:
        filter_one(db_name)
    print("✅ step02 done")

if __name__ == "__main__":
    main()