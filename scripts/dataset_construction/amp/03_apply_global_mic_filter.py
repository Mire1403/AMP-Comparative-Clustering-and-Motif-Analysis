from pathlib import Path
import pandas as pd
import re

# =====================================================
# PATH CONFIG (RELATIVE TO REPO)
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_INTERMEDIATE_DIR = PROJECT_ROOT / "data" / "intermediate"

DB_NAMES = ["CAMP", "DBAASP", "dbAMP3", "DRAMP"]

THRESHOLD = 65  # µg/mL

AA_WEIGHTS = {
    'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1,
    'C': 121.2, 'E': 147.1, 'Q': 146.2, 'G': 75.1,
    'H': 155.2, 'I': 131.2, 'L': 131.2, 'K': 146.2,
    'M': 149.2, 'F': 165.2, 'P': 115.1, 'S': 105.1,
    'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1
}

UNIT_PATTERNS = {
    "uM": r"(µm|μm|um|micromolar)",
    "ug/mL": r"(µg\/ml|μg\/ml|ug\/ml|microgram\/ml|micrograms\/ml)",
    "mg/mL": r"(mg\/ml)",
    "ng/mL": r"(ng\/ml|nanogram\/ml|nanograms\/ml)"
}

# =====================================================
# FUNCTIONS
# =====================================================

def calculate_mw(sequence):
    return sum(AA_WEIGHTS.get(aa, 0) for aa in sequence)

def extract_all_mic_values(text):

    if pd.isna(text):
        return []

    text = str(text).lower()

    if "mic" not in text:
        return []

    pattern = r"mic[^0-9]*([0-9]+\.?[0-9]*)"
    values = re.findall(pattern, text)

    return [float(v) for v in values]

def detect_unit(text):

    text = str(text).lower()

    for unit_name, pattern in UNIT_PATTERNS.items():
        if re.search(pattern, text):
            return unit_name

    return None

def convert_to_ugml(value, unit, mw):

    if unit == "ug/mL":
        return value
    elif unit == "mg/mL":
        return value * 1000
    elif unit == "ng/mL":
        return value / 1000
    elif unit == "uM":
        if mw == 0:
            return None
        return (value * mw) / 1000
    else:
        return None

# =====================================================
# MAIN PROCESS
# =====================================================

def process_database(db_name):

    input_file = DATA_INTERMEDIATE_DIR / f"{db_name}_structural_filtered.xlsx"

    if not input_file.exists():
        print(f"No {db_name}_structural_filtered.xlsx found.")
        return

    df = pd.read_excel(input_file)

    total_initial = len(df)
    removed = 0

    for idx, row in df.iterrows():

        combined_text = " ".join([
            str(row.get("Activity", "")),
            str(row.get("Target Group", "")),
            str(row.get("Target Object", ""))
        ])

        mic_values = extract_all_mic_values(combined_text)

        if not mic_values:
            continue

        unit = detect_unit(combined_text)
        mw = calculate_mw(row["Sequence"])

        converted_values = []

        for v in mic_values:
            conv = convert_to_ugml(v, unit, mw)
            if conv is not None:
                converted_values.append(conv)

        if not converted_values:
            continue

        min_mic = min(converted_values)

        if min_mic >= THRESHOLD:
            df.drop(index=idx, inplace=True)
            removed += 1

    total_final = len(df)

    print("\n========================================")
    print(f"GLOBAL MIC FILTER - {db_name}")
    print("========================================")
    print("Initial:", total_initial)
    print("Removed (min MIC ≥ 65 ug/mL):", removed)
    print("Final:", total_final)
    print("========================================\n")

    # 🔁 Guardamos nuevo archivo coherente con pipeline
    output_file = DATA_INTERMEDIATE_DIR / f"{db_name}_activity_filtered.xlsx"
    df.to_excel(output_file, index=False)


# =====================================================
# RUN
# =====================================================

if __name__ == "__main__":

    for db_name in DB_NAMES:
        process_database(db_name)

    print("Step 03 (activity filtering - MIC) completed successfully.")