import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

FMAP_FILE = PROJECT_ROOT / "results" / "11_mechanism_analysis" / "1_motif_activity_analysis" / "fmap" / "fmap_best_peptides.csv"
PPM_FILE = PROJECT_ROOT / "results" / "11_mechanism_analysis" / "1_motif_activity_analysis" / "ppm" / "ppm_processed.csv"

OUTPUT_FILE = PROJECT_ROOT / "results" / "11_mechanism_analysis" / "1_motif_activity_analysis" / "comparison" / "fmap_ppm_comparison.xlsx"
OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

# =========================
# LOAD
# =========================
fmap = pd.read_csv(FMAP_FILE)
ppm = pd.read_csv(PPM_FILE)

# =========================
# DETECT MOTIF COLUMN
# =========================
if "Motif_ID" in fmap.columns:
    motif_col = "Motif_ID"
elif "Motif" in fmap.columns:
    motif_col = "Motif"
else:
    raise ValueError("❌ No motif column found (Motif or Motif_ID)")

print(f"Using motif column: {motif_col}")

# =========================
# RENAME
# =========================
fmap = fmap.rename(columns={
    motif_col: "Motif",
    "FMAP_Class": "Class_FMAP",
    "Mechanism": "Mechanism_FMAP",
    "Best_Energy": "Best_Energy_FMAP",
    "Best_Angle": "Best_Angle_FMAP",
    "Best_Depth": "Best_Depth_FMAP"
})

ppm = ppm.rename(columns={
    "PPM_Class": "Class_PPM",
    "Mechanism": "Mechanism_PPM",
    "Best_Energy": "Best_Energy_PPM",
    "Best_Angle": "Best_Angle_PPM",
    "Best_Depth": "Best_Depth_PPM"
})

# =========================
# SELECT
# =========================
fmap = fmap[[
    "Peptide_ID",
    "Motif",
    "Class_FMAP",
    "Mechanism_FMAP",
    "Best_Energy_FMAP",
    "Best_Angle_FMAP",
    "Best_Depth_FMAP"
]]

ppm = ppm[[
    "Peptide_ID",
    "Class_PPM",
    "Mechanism_PPM",
    "Best_Energy_PPM",
    "Best_Angle_PPM",
    "Best_Depth_PPM"
]]

# =========================
# MERGE
# =========================
merged = pd.merge(fmap, ppm, on="Peptide_ID", how="inner")

# =========================
# ACTIVITY AGREEMENT
# =========================
order = {
    "No_interaction": 0,
    "Weak_interaction": 1,
    "Weak_insertion": 2,
    "Moderate_insertion": 3,
    "Strong_insertion": 4
}

merged["Class_FMAP_num"] = merged["Class_FMAP"].map(order)
merged["Class_PPM_num"] = merged["Class_PPM"].map(order)

merged["Class_diff"] = (merged["Class_FMAP_num"] - merged["Class_PPM_num"]).abs()

def classify_agreement(diff):
    if pd.isna(diff):
        return "Unknown"
    if diff == 0:
        return "Exact"
    if diff == 1:
        return "Close"
    return "Conflict"

merged["Activity_Agreement"] = merged["Class_diff"].apply(classify_agreement)

# =========================
# MECHANISM AGREEMENT
# =========================
merged["Mechanism_Agreement"] = merged.apply(
    lambda r: "Match" if r["Mechanism_FMAP"] == r["Mechanism_PPM"] else "Different",
    axis=1
)

# =========================
# METRICS
# =========================
total = len(merged)

metrics = pd.DataFrame({
    "Metric": [
        "Activity Exact",
        "Activity Close",
        "Activity Conflict",
        "Mechanism Match"
    ],
    "Value": [
        (merged["Activity_Agreement"] == "Exact").sum() / total,
        (merged["Activity_Agreement"] == "Close").sum() / total,
        (merged["Activity_Agreement"] == "Conflict").sum() / total,
        (merged["Mechanism_Agreement"] == "Match").sum() / total
    ]
})

# =========================
# 🔥 MOTIF ANALYSIS
# =========================
motif_counts = merged["Motif"].value_counts()

motif_mech = pd.crosstab(
    merged["Motif"],
    merged["Mechanism_PPM"]
)

motif_mech_norm = pd.crosstab(
    merged["Motif"],
    merged["Mechanism_PPM"],
    normalize="index"
)

motif_mech["Total"] = motif_counts
motif_mech_norm["Total"] = motif_counts

motif_mech = motif_mech.sort_values("Total", ascending=False)
motif_mech_norm = motif_mech_norm.loc[motif_mech.index]

# =========================
# SAVE
# =========================
with pd.ExcelWriter(OUTPUT_FILE, engine="openpyxl") as writer:
    merged.to_excel(writer, sheet_name="Merged", index=False)
    metrics.to_excel(writer, sheet_name="Metrics", index=False)
    motif_mech.to_excel(writer, sheet_name="Motif_Counts")
    motif_mech_norm.to_excel(writer, sheet_name="Motif_Normalized")

# =========================
# PRINT
# =========================
print("\nSaved to:")
print(OUTPUT_FILE)

print("\n=== ACTIVITY AGREEMENT ===")
print(merged["Activity_Agreement"].value_counts())

print("\n=== MECHANISM AGREEMENT ===")
print(merged["Mechanism_Agreement"].value_counts())