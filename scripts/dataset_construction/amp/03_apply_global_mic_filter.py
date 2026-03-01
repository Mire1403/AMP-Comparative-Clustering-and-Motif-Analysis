from __future__ import annotations

from pathlib import Path
import re
import pandas as pd
import numpy as np

# =====================================================
# PATH CONFIG (RELATIVE TO REPO)
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_INTERMEDIATE_DIR = PROJECT_ROOT / "data" / "intermediate"

DB_NAMES = ["CAMP", "DBAASP", "dbAMP3", "DRAMP"]

THRESHOLD_UGML = 65.0  # µg/mL

# Approximate residue masses (Da) for quick MW estimate (good enough for unit conversion)
AA_WEIGHTS = {
    "A": 89.1, "R": 174.2, "N": 132.1, "D": 133.1,
    "C": 121.2, "E": 147.1, "Q": 146.2, "G": 75.1,
    "H": 155.2, "I": 131.2, "L": 131.2, "K": 146.2,
    "M": 149.2, "F": 165.2, "P": 115.1, "S": 105.1,
    "T": 119.1, "W": 204.2, "Y": 181.2, "V": 117.1
}

# Unit patterns (robust)
UNIT_PATTERNS = {
    "uM": r"(?:µm|μm|um|micromolar|µmol\/l|μmol\/l|umol\/l|micromol\/l)",
    "ug/mL": r"(?:µg\/ml|μg\/ml|ug\/ml|micrograms?\/ml)",
    "mg/mL": r"(?:mg\/ml)",
    "ng/mL": r"(?:ng\/ml|nanograms?\/ml)",
}

# MIC keyword (include MIC50/MIC90 etc., but treat them as MIC values too)
MIC_KEYWORD = r"(?:\bmic\b|\bmic50\b|\bmic90\b)"


# =====================================================
# HELPERS
# =====================================================

def calculate_mw(sequence: str) -> float:
    seq = str(sequence).upper().strip()
    return float(sum(AA_WEIGHTS.get(aa, 0.0) for aa in seq))


def _normalize_text(x) -> str:
    if pd.isna(x):
        return ""
    return str(x).replace("≤", "<=").replace("≥", ">=").lower()


def detect_unit_near(text: str) -> str | None:
    """Detect a unit anywhere in text (fallback)."""
    for unit, pat in UNIT_PATTERNS.items():
        if re.search(pat, text, flags=re.I):
            return unit
    return None


def extract_mic_measurements(text: str) -> list[tuple[float, str | None, str]]:
    """
    Extract MIC measurements as (value, unit, raw_snippet).

    Tries to capture:
    - 'MIC 2 ug/ml'
    - 'MIC=2-4 ug/ml' -> takes both numbers
    - 'MIC <= 1 uM'
    - 'MIC: 0.5 μg/ml'
    """
    text = _normalize_text(text)
    if not re.search(MIC_KEYWORD, text, flags=re.I):
        return []

    # Capture patterns where MIC keyword appears before number(s) and optional unit nearby
    # Examples matched:
    # mic 2 ug/ml
    # mic: 2-4 ug/ml
    # mic <= 1 um
    pattern = re.compile(
        rf"({MIC_KEYWORD}[^0-9<>=]*"
        rf"(?:<=|>=|=)?\s*"
        rf"([0-9]+(?:\.[0-9]+)?)(?:\s*-\s*([0-9]+(?:\.[0-9]+)?))?"
        rf"(?:\s*([a-zµμ\/]+))?)",
        flags=re.I
    )

    measurements: list[tuple[float, str | None, str]] = []

    for m in pattern.finditer(text):
        snippet = m.group(1)
        v1 = float(m.group(2))
        v2 = m.group(3)
        unit_token = (m.group(4) or "").strip()

        # Try detect unit from the token first, else from the snippet, else from whole text
        unit = None
        if unit_token:
            for u, upat in UNIT_PATTERNS.items():
                if re.search(upat, unit_token, flags=re.I):
                    unit = u
                    break
            # also allow direct tokens like "um"
            if unit is None and unit_token in {"um", "µm", "μm"}:
                unit = "uM"

        if unit is None:
            unit = detect_unit_near(snippet)
        if unit is None:
            unit = detect_unit_near(text)

        measurements.append((v1, unit, snippet))
        if v2 is not None:
            measurements.append((float(v2), unit, snippet))

    return measurements


def convert_to_ugml(value: float, unit: str | None, mw: float) -> float | None:
    if unit == "ug/mL":
        return value
    if unit == "mg/mL":
        return value * 1000.0
    if unit == "ng/mL":
        return value / 1000.0
    if unit == "uM":
        if mw <= 0:
            return None
        # uM * (g/mol) -> ug/mL : (value * mw) / 1000
        return (value * mw) / 1000.0
    return None


# =====================================================
# MAIN
# =====================================================

def process_database(db_name: str) -> None:
    input_file = DATA_INTERMEDIATE_DIR / f"{db_name}_structural_filtered.xlsx"
    if not input_file.exists():
        print(f"⚠️  Missing: {input_file.name}")
        return

    df = pd.read_excel(input_file)

    total_initial = len(df)
    if total_initial == 0:
        print(f"⚠️  Empty file: {input_file.name}")
        return

    # Build combined text (same idea as you had)
    combined = (
        df.get("Activity", "").astype(str).fillna("") + " " +
        df.get("Target Group", "").astype(str).fillna("") + " " +
        df.get("Target Object", "").astype(str).fillna("")
    )

    # MW per sequence (vectorized-ish)
    seqs = df.get("Sequence", "").astype(str).fillna("")
    mw_series = seqs.apply(calculate_mw)

    # Extract MICs per row
    min_mic_ugml = []
    unit_detected = []
    raw_snippets = []

    for text, mw in zip(combined.tolist(), mw_series.tolist()):
        meas = extract_mic_measurements(text)

        conv_vals = []
        units_here = set()
        snippets_here = []

        for v, u, snip in meas:
            cv = convert_to_ugml(v, u, mw)
            if cv is not None:
                conv_vals.append(cv)
                if u:
                    units_here.add(u)
                snippets_here.append(snip)

        if conv_vals:
            min_mic_ugml.append(float(min(conv_vals)))
            unit_detected.append(";".join(sorted(units_here)) if units_here else "")
            # store a short joined snippet (avoid huge fields)
            raw_snippets.append(" | ".join(snippets_here)[:500])
        else:
            min_mic_ugml.append(np.nan)
            unit_detected.append("")
            raw_snippets.append("")

    df["Min_MIC_ug_mL"] = min_mic_ugml
    df["MIC_unit_detected"] = unit_detected
    df["MIC_values_raw"] = raw_snippets

    # Filtering rule: remove if we have MIC and min >= threshold
    has_mic = df["Min_MIC_ug_mL"].notna()
    to_remove = has_mic & (df["Min_MIC_ug_mL"] >= THRESHOLD_UGML)

    removed = int(to_remove.sum())
    df_filtered = df.loc[~to_remove].copy()

    print("\n========================================")
    print(f"GLOBAL MIC FILTER - {db_name}")
    print("========================================")
    print("Initial:", total_initial)
    print("Rows with MIC parsed:", int(has_mic.sum()))
    print(f"Removed (min MIC ≥ {THRESHOLD_UGML:g} ug/mL):", removed)
    print("Final:", len(df_filtered))
    print("========================================\n")

    output_file = DATA_INTERMEDIATE_DIR / f"{db_name}_activity_filtered.xlsx"
    df_filtered.to_excel(output_file, index=False)
    print("Saved:", output_file)


if __name__ == "__main__":
    for db_name in DB_NAMES:
        process_database(db_name)

    print("\n✅ Step 03 (activity filtering - MIC) completed successfully.\n")