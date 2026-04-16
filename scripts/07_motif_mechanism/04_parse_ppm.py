import pandas as pd
import re
import os

# =========================
# CONFIG
# =========================

BASE_DIR = "results/statistics/10_ppm_output"

DATAPAR1 = f"{BASE_DIR}/datapar1"
DATAPAR2 = f"{BASE_DIR}/datapar2"
DATASUB1 = f"{BASE_DIR}/datasub1"

OUT_CSV = f"{BASE_DIR}/ppm_results_clean.csv"
OUT_XLSX = f"{BASE_DIR}/ppm_results_clean.xlsx"

# =========================
# NORMALIZE ID
# =========================

def normalize_id(filename):
    """
    Convert:
    15.pdb, 015.pdb, 00015.pdb, P0015.pdb → P0015
    """
    name = filename.replace(".pdb", "")
    name = name.replace("P", "")

    num = int(name)  # elimina ceros a la izquierda
    return f"P{num:04d}"

# =========================
# PARSE DATAPAR1
# =========================

def parse_datapar1(file):

    rows = []

    with open(file) as f:
        for line in f:
            parts = line.strip().split(";")
            if len(parts) < 7:
                continue

            fname = parts[0].strip()
            pid = normalize_id(fname)

            rows.append({
                "Peptide_ID": pid,
                "Depth": float(parts[1]),
                "Tilt": float(parts[2]),
                "Angle": float(parts[3]),
                "Shift": float(parts[4]),
                "DeltaG1": float(parts[5]),
                "DeltaG2": float(parts[6]),
            })

    df = pd.DataFrame(rows)

    # eliminar duplicados
    df = df.drop_duplicates(subset=["Peptide_ID"])

    return df

# =========================
# PARSE DATAPAR2
# =========================

def parse_datapar2(file):

    rows = []

    with open(file) as f:
        for line in f:
            parts = line.strip().split(";")
            if len(parts) < 6:
                continue

            fname = parts[0].strip()
            pid = normalize_id(fname)

            rows.append({
                "Peptide_ID": pid,
                "Param1": float(parts[1]),
                "Param2": float(parts[2]),
                "Param3": float(parts[3]),
                "Energy1": float(parts[4]),
                "Energy2": float(parts[5]),
            })

    df = pd.DataFrame(rows)
    df = df.drop_duplicates(subset=["Peptide_ID"])

    return df

# =========================
# PARSE DATASUB1
# =========================

def parse_datasub1(file):

    rows = []

    with open(file) as f:
        for line in f:

            parts = line.strip().split(";")
            if len(parts) < 4:
                continue

            fname = parts[0].strip()
            pid = normalize_id(fname)

            # extraer rango (ej: 3-21)
            m = re.search(r"\((\s*\d+)-(\s*\d+)\)", parts[3])
            if m:
                start = int(m.group(1))
                end = int(m.group(2))
            else:
                start, end = None, None

            rows.append({
                "Peptide_ID": pid,
                "Chain": parts[1].strip(),
                "Segment_ID": int(parts[2]),
                "TM_Start": start,
                "TM_End": end
            })

    df = pd.DataFrame(rows)
    df = df.drop_duplicates(subset=["Peptide_ID"])

    return df

# =========================
# MAIN
# =========================

if __name__ == "__main__":

    print("Parsing PPM files...")

    df1 = parse_datapar1(DATAPAR1)
    df2 = parse_datapar2(DATAPAR2)
    df3 = parse_datasub1(DATASUB1)

    print("Merging datasets...")

    df = df1.merge(df2, on="Peptide_ID", how="outer")
    df = df.merge(df3, on="Peptide_ID", how="outer")

    print("Total unique peptides:", len(df))

    # =========================
    # SAVE
    # =========================

    os.makedirs(BASE_DIR, exist_ok=True)

    df.to_csv(OUT_CSV, index=False)
    df.to_excel(OUT_XLSX, index=False)

    print("\n✅ PPM clean dataset generated:")
    print(f"CSV: {OUT_CSV}")
    print(f"Excel: {OUT_XLSX}")