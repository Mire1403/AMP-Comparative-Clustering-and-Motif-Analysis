from __future__ import annotations

from pathlib import Path
import re
import numpy as np
import pandas as pd

from config import DATA_INTERMEDIATE_DIR, MIC_THRESHOLD_UGML

DBS = ["CAMP", "DBAASP", "dbAMP3", "DRAMP"]

AA_WEIGHTS = {
    "A": 89.1, "R": 174.2, "N": 132.1, "D": 133.1,
    "C": 121.2, "E": 147.1, "Q": 146.2, "G": 75.1,
    "H": 155.2, "I": 131.2, "L": 131.2, "K": 146.2,
    "M": 149.2, "F": 165.2, "P": 115.1, "S": 105.1,
    "T": 119.1, "W": 204.2, "Y": 181.2, "V": 117.1
}

UNIT_PATTERNS = {
    "uM": r"(?:µm|μm|um|micromolar|µmol\/l|μmol\/l|umol\/l|micromol\/l)",
    "ug/mL": r"(?:µg\/ml|μg\/ml|ug\/ml|micrograms?\/ml)",
    "mg/mL": r"(?:mg\/ml)",
    "ng/mL": r"(?:ng\/ml|nanograms?\/ml)",
}

MIC_KEYWORD = r"(?:\bmic\b|\bmic50\b|\bmic90\b)"

REPORT_FILE = DATA_INTERMEDIATE_DIR / "step03_mic_filter_report.csv"

def calculate_mw(sequence: str) -> float:
    seq = str(sequence).upper().strip()
    return float(sum(AA_WEIGHTS.get(aa, 0.0) for aa in seq))


def _norm_text(x) -> str:
    if pd.isna(x):
        return ""
    return str(x).replace("≤", "<=").replace("≥", ">=").lower()


def detect_unit(text: str) -> str | None:
    for unit, pat in UNIT_PATTERNS.items():
        if re.search(pat, text, flags=re.I):
            return unit
    return None


def extract_mic(text: str) -> list[tuple[float, str | None]]:
    """
    Extract MIC values as (value, unit). Supports ranges '2-4' (returns both).
    """
    text = _norm_text(text)
    if not re.search(MIC_KEYWORD, text, flags=re.I):
        return []

    pat = re.compile(
        rf"{MIC_KEYWORD}[^0-9<>=]*"
        rf"(?:<=|>=|=)?\s*"
        rf"([0-9]+(?:\.[0-9]+)?)(?:\s*-\s*([0-9]+(?:\.[0-9]+)?))?"
        rf"(?:\s*([a-zµμ\/]+))?",
        flags=re.I
    )

    out: list[tuple[float, str | None]] = []
    for m in pat.finditer(text):
        v1 = float(m.group(1))
        v2 = m.group(2)
        token = (m.group(3) or "").strip()
        unit = None

        if token:
            for u, upat in UNIT_PATTERNS.items():
                if re.search(upat, token, flags=re.I):
                    unit = u
                    break
            if unit is None and token in {"um", "µm", "μm"}:
                unit = "uM"

        if unit is None:
            unit = detect_unit(text)

        out.append((v1, unit))
        if v2:
            out.append((float(v2), unit))

    return out


def to_ugml(value: float, unit: str | None, mw: float) -> float | None:
    if unit == "ug/mL":
        return value
    if unit == "mg/mL":
        return value * 1000.0
    if unit == "ng/mL":
        return value / 1000.0
    if unit == "uM":
        return (value * mw) / 1000.0 if mw > 0 else None
    return None


def filter_one(db_name: str) -> Path:
    in_file = DATA_INTERMEDIATE_DIR / f"{db_name}_02_structural_filtered.parquet"
    if not in_file.exists():
        print(f"⚠️ {db_name}: missing input -> {in_file.name} (skipping)")
        return in_file
    df = pd.read_parquet(in_file).copy()
    if "sequence" not in df.columns:
        print(f"⚠️ {db_name}: no 'sequence' column (skipping)")
        return in_file
    n0 = len(df)
    if n0 == 0:
        raise ValueError(f"Empty input file: {in_file}")

    # Build MIC search text without creating "nan"
    combined = (
        df.get("activity", pd.Series([""] * n0)).astype("string").fillna("") + " " +
        df.get("target_group", pd.Series([""] * n0)).astype("string").fillna("") + " " +
        df.get("target_object", pd.Series([""] * n0)).astype("string").fillna("")
    )

    mw = df["sequence"].astype("string").fillna("").apply(calculate_mw)

    min_mic = []
    unit_detected = []

    for txt, m in zip(combined.tolist(), mw.tolist()):
        vals = []
        units = set()

        for v, u in extract_mic(txt):
            cv = to_ugml(v, u, m)
            if cv is not None:
                vals.append(cv)
                if u:
                    units.add(u)

        min_mic.append(float(min(vals)) if vals else np.nan)
        unit_detected.append(";".join(sorted(units)) if units else "")

    df["min_mic_ugml"] = min_mic
    df["mic_unit_detected"] = unit_detected

    has = df["min_mic_ugml"].notna()
    rm = has & (df["min_mic_ugml"] >= MIC_THRESHOLD_UGML)

    df_out = df.loc[~rm].copy()
    out_file = DATA_INTERMEDIATE_DIR / f"{db_name}_03_mic_filtered.parquet"
    df_out.to_parquet(out_file, index=False)

    print(
        f"{db_name}: start={n0} "
        f"mic_parsed={int(has.sum())} "
        f"removed_mic={int(rm.sum())} "
        f"final={len(df_out)} "
        f"-> {out_file.name}"
    )
    return out_file


def main() -> None:
    for db in DBS:
        filter_one(db)
    print("✅ step03 done")


if __name__ == "__main__":
    main()