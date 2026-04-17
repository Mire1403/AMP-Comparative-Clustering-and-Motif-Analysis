import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_FILE = PROJECT_ROOT / "results" / "10_motif_mechanism_output" / "01_fmap_output" / "fmap_final_results.csv"
OUTPUT_DIR = PROJECT_ROOT / "results" / "11_mechanism_analysis" / "1_motif_activity_analysis" / "fmap"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT_FILE)

if "Sequence_x" in df.columns:
    df["Sequence"] = df["Sequence_x"]
elif "Sequence_y" in df.columns:
    df["Sequence"] = df["Sequence_y"]
else:
    raise ValueError("No sequence column found")

NUMERIC_COLUMNS = [
    "Membrane_Binding_Energy",
    "Tilt_Angle",
    "Depth_Thickness",
    "Helix_Start",
    "Helix_End",
    "Motif_Start",
    "Motif_End"
]

for col in NUMERIC_COLUMNS:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

df = df.dropna(subset=["Membrane_Binding_Energy", "Tilt_Angle"])

df["Helix_Length"] = df["Helix_End"] - df["Helix_Start"] + 1

def compute_overlap(row):
    if pd.isna(row["Helix_Start"]) or pd.isna(row["Helix_End"]) or pd.isna(row["Motif_Start"]) or pd.isna(row["Motif_End"]):
        return 0
    start = max(row["Helix_Start"], row["Motif_Start"])
    end = min(row["Helix_End"], row["Motif_End"])
    return max(0, end - start + 1)

df["Motif_Overlap"] = df.apply(compute_overlap, axis=1)

def is_real_helix(val):
    if pd.isna(val):
        return 0
    if isinstance(val, str):
        return 1 if val.upper() == "FALSE" else 0
    return 1 if val is False else 0

df["Real_Helix"] = df["Is_False_Helix"].apply(is_real_helix)

df = df.sort_values(
    by=["Peptide_ID", "Motif_Overlap", "Real_Helix", "Membrane_Binding_Energy"],
    ascending=[True, False, False, True]
)

best = df.groupby("Peptide_ID", as_index=False).first()

best["Best_Energy"] = best["Membrane_Binding_Energy"]
best["Best_Angle"] = best["Tilt_Angle"]
best["Best_Depth"] = best["Depth_Thickness"]

def classify_common(row):
    e = row["Best_Energy"]

    if pd.isna(e):
        return "Unknown"
    if e < -12:
        return "Strong_insertion"
    if e < -8:
        return "Moderate_insertion"
    if e < -5:
        return "Weak_insertion"
    if e < -2:
        return "Weak_interaction"
    return "No_interaction"

def classify_mechanism(row):
    a = row["Best_Angle"]
    d = row["Best_Depth"]

    if pd.isna(a) or pd.isna(d):
        return "Unknown"

    if d >= 5 and a >= 75:
        return "Barrel_stave"

    if d >= 3 and a >= 50:
        return "Toroidal_pore"

    if d < 3:
        return "Carpet"

    return "Other"

best["FMAP_Class"] = best.apply(classify_common, axis=1)
best["Mechanism"] = best.apply(classify_mechanism, axis=1)

out = OUTPUT_DIR / "fmap_best_peptides.csv"
best.to_csv(out, index=False)

plt.figure()
plt.scatter(best["Best_Angle"], best["Best_Energy"])
plt.xlabel("Angle")
plt.ylabel("Energy")
plt.savefig(OUTPUT_DIR / "energy_vs_angle.png")
plt.close()

plt.figure()
plt.hist(best["Best_Energy"], bins=30)
plt.savefig(OUTPUT_DIR / "energy_distribution.png")
plt.close()

print(out)
print("\nClass distribution:")
print(best["FMAP_Class"].value_counts())
print("\nMechanism distribution:")
print(best["Mechanism"].value_counts())