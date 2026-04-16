import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_FILE = PROJECT_ROOT / "results" / "statistics" / "11_motif_activity_analysis" / "comparison" / "fmap_ppm_comparison.xlsx"
OUTPUT_FILE = PROJECT_ROOT / "results" / "statistics" / "11_motif_activity_analysis" / "comparison" / "motif_like_analysis.xlsx"

df = pd.read_excel(INPUT_FILE, sheet_name="Merged")

# =========================
# 🔥 CHECK MOTIF COLUMN
# =========================
if "Motif" not in df.columns:
    raise ValueError("❌ No 'Motif' column found in merged data")

# =========================
# 1. GLOBAL DISTRIBUTION
# =========================
fmap_mech = df["Mechanism_FMAP"].value_counts(normalize=True)
ppm_mech = df["Mechanism_PPM"].value_counts(normalize=True)

global_dist = pd.DataFrame({
    "FMAP": fmap_mech,
    "PPM": ppm_mech
}).fillna(0)

# =========================
# 2. CROSS MODEL
# =========================
cross_mech = pd.crosstab(
    df["Mechanism_FMAP"],
    df["Mechanism_PPM"]
)

cross_mech_norm = pd.crosstab(
    df["Mechanism_FMAP"],
    df["Mechanism_PPM"],
    normalize="index"
)

# =========================
# 🔥🔥 3. MOTIF ANALYSIS (CLAVE)
# =========================

# Conteo absoluto
motif_fmap = pd.crosstab(
    df["Motif"],
    df["Mechanism_FMAP"]
)

motif_ppm = pd.crosstab(
    df["Motif"],
    df["Mechanism_PPM"]
)

# Normalizado (MUY IMPORTANTE 🔥)
motif_fmap_norm = pd.crosstab(
    df["Motif"],
    df["Mechanism_FMAP"],
    normalize="index"
)

motif_ppm_norm = pd.crosstab(
    df["Motif"],
    df["Mechanism_PPM"],
    normalize="index"
)

# Tamaño del motif
motif_size = df["Motif"].value_counts()

# Añadir tamaño
motif_fmap["Total"] = motif_size
motif_ppm["Total"] = motif_size
motif_fmap_norm["Total"] = motif_size
motif_ppm_norm["Total"] = motif_size

# Ordenar por tamaño
motif_fmap = motif_fmap.sort_values("Total", ascending=False)
motif_ppm = motif_ppm.loc[motif_fmap.index]
motif_fmap_norm = motif_fmap_norm.loc[motif_fmap.index]
motif_ppm_norm = motif_ppm_norm.loc[motif_fmap.index]

# =========================
# 🔥 4. HIGH CONFIDENCE POR MOTIF
# =========================
high = df[df["Activity_Agreement"] == "Exact"]

motif_high = pd.crosstab(
    high["Motif"],
    high["Mechanism_PPM"],
    normalize="index"
)

# =========================
# 🔥 5. CONFLICTOS POR MOTIF
# =========================
conflicts = df[df["Mechanism_Agreement"] == "Different"]

motif_conflicts = pd.crosstab(
    conflicts["Motif"],
    conflicts["Mechanism_PPM"]
)

# =========================
# SAVE
# =========================
with pd.ExcelWriter(OUTPUT_FILE, engine="openpyxl") as writer:
    global_dist.to_excel(writer, sheet_name="Global_Distribution")
    cross_mech.to_excel(writer, sheet_name="Mechanism_Cross")
    cross_mech_norm.to_excel(writer, sheet_name="Mechanism_Cross_Normalized")

    # 🔥 MOTIF
    motif_fmap.to_excel(writer, sheet_name="Motif_FMAP_Counts")
    motif_ppm.to_excel(writer, sheet_name="Motif_PPM_Counts")
    motif_fmap_norm.to_excel(writer, sheet_name="Motif_FMAP_Normalized")
    motif_ppm_norm.to_excel(writer, sheet_name="Motif_PPM_Normalized")

    motif_high.to_excel(writer, sheet_name="Motif_HighConfidence")
    motif_conflicts.to_excel(writer, sheet_name="Motif_Conflicts")

print(OUTPUT_FILE)