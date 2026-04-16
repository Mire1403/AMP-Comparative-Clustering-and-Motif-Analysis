import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_FILE = PROJECT_ROOT / "results" / "statistics" / "10_ppm_output" / "ppm_results_clean.csv"
OUTPUT_DIR = PROJECT_ROOT / "results" / "statistics" / "11_motif_activity_analysis" / "ppm"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT_FILE)

cols = [
    "Depth","Tilt","Angle","Shift",
    "DeltaG1","DeltaG2",
    "Param1","Param2","Param3",
    "Energy1","Energy2"
]

for c in cols:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")

df["Best_Energy"] = df[["Energy1","Energy2","DeltaG1","DeltaG2"]].min(axis=1)

df["Best_Angle"] = df["Angle"].fillna(df["Param3"])
df["Best_Depth"] = df["Depth"].fillna(df["Param1"])

df = df.dropna(subset=["Best_Energy", "Best_Angle"])

# =========================
# ACTIVITY (MISMO QUE FMAP)
# =========================
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

# =========================
# MECHANISM (MISMO QUE FMAP)
# =========================
def classify_mechanism(row):
    a = row["Best_Angle"]
    d = row["Best_Depth"]

    if pd.isna(a) or pd.isna(d):
        return "Unknown"

    if d >= 6 and a >= 75:
        return "Barrel_stave"

    if d >= 2.5 and a >= 50:
        return "Toroidal_pore"

    if d < 2.5:
        return "Carpet"

    return "Other"

df["PPM_Class"] = df.apply(classify_common, axis=1)
df["Mechanism"] = df.apply(classify_mechanism, axis=1)

out = OUTPUT_DIR / "ppm_processed.csv"
df.to_csv(out, index=False)

# =========================
# PLOTS
# =========================
plt.figure()
plt.scatter(df["Best_Angle"], df["Best_Energy"])
plt.xlabel("Angle")
plt.ylabel("Energy")
plt.savefig(OUTPUT_DIR / "ppm_energy_vs_angle.png")
plt.close()

plt.figure()
plt.hist(df["Best_Energy"], bins=30)
plt.savefig(OUTPUT_DIR / "ppm_energy_hist.png")
plt.close()

if "Best_Depth" in df.columns:
    plt.figure()
    plt.scatter(df["Best_Depth"], df["Best_Energy"])
    plt.xlabel("Depth")
    plt.ylabel("Energy")
    plt.savefig(OUTPUT_DIR / "ppm_depth_vs_energy.png")
    plt.close()

print(out)
print("\nClass distribution:")
print(df["PPM_Class"].value_counts())

print("\nMechanism distribution:")
print(df["Mechanism"].value_counts())

print("\nSanity check:")
print("Total:", len(df))
print("With depth:", df["Best_Depth"].notna().sum())