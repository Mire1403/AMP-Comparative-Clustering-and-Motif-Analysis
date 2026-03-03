
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

print("\nFIMO ENRICHMENT + FDR PIPELINE \n")

# =====================================================
# REPO ROOT (robust)
# =====================================================

def find_repo_root(start: Path) -> Path:
    markers = [".git", "README.md", "requirements.txt", "pyproject.toml", "environment.yml"]
    start = start.resolve()
    for parent in [start] + list(start.parents):
        if any((parent / m).exists() for m in markers):
            return parent
    # fallback: assume scripts/<category>/<file>.py
    return start.parents[2]

PROJECT_ROOT = find_repo_root(Path(__file__).parent)

# =====================================================
# PATH CONFIG
# =====================================================

FIMO_BASE = PROJECT_ROOT / "results" / "motif_scanning"
OUT_DIR = PROJECT_ROOT / "results" / "statistics" / "fimo_enrichment"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# =====================================================
# CONFIG
# =====================================================

ALPHA_FISHER = 0.05
ALPHA_FDR = 0.05
PSEUDOCOUNT = 1  # enrichment smoothing

# Expected locations produced by your FIMO script:
# results/motif_scanning/{cdhit|mmseq}/{meme|streme}/{fimo_amps|fimo_nonamps}/fimo.tsv
STRUCTURE = {
    ("cdhit", "meme"): {
        "amps": FIMO_BASE / "cdhit" / "meme" / "fimo_amps" / "fimo.tsv",
        "nonamps": FIMO_BASE / "cdhit" / "meme" / "fimo_nonamps" / "fimo.tsv",
    },
    ("cdhit", "streme"): {
        "amps": FIMO_BASE / "cdhit" / "streme" / "fimo_amps" / "fimo.tsv",
        "nonamps": FIMO_BASE / "cdhit" / "streme" / "fimo_nonamps" / "fimo.tsv",
    },
    ("mmseq", "meme"): {
        "amps": FIMO_BASE / "mmseq" / "meme" / "fimo_amps" / "fimo.tsv",
        "nonamps": FIMO_BASE / "mmseq" / "meme" / "fimo_nonamps" / "fimo.tsv",
    },
    ("mmseq", "streme"): {
        "amps": FIMO_BASE / "mmseq" / "streme" / "fimo_amps" / "fimo.tsv",
        "nonamps": FIMO_BASE / "mmseq" / "streme" / "fimo_nonamps" / "fimo.tsv",
    },
}

# =====================================================
# HELPERS
# =====================================================

def die(msg: str) -> None:
    print(f"{msg}")
    sys.exit(1)

def read_fimo_tsv(path: Path) -> pd.DataFrame:
    """Read FIMO TSV (skip comments). Return empty DF if missing/invalid."""
    if not path.exists():
        return pd.DataFrame()
    try:
        df = pd.read_csv(path, sep="\t", comment="#")
    except Exception:
        return pd.DataFrame()

    needed = {"motif_id", "sequence_name"}
    if not needed.issubset(df.columns):
        return pd.DataFrame()
    return df

def safe_fisher(table):
    """table = [[a,b],[c,d]]"""
    try:
        oddsratio, pvalue = fisher_exact(table)
        return oddsratio, pvalue
    except Exception:
        return None, None

def save_excel(df: pd.DataFrame, path: Path, sheet_name: str = "Results") -> None:
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name=sheet_name)

# =====================================================
# 1) BUILD ENRICHMENT TABLE
# =====================================================

all_results: list[dict] = []
any_loaded = False

for (clustering, tool), files in STRUCTURE.items():
    amps_path = files["amps"]
    nonamps_path = files["nonamps"]

    df_amp = read_fimo_tsv(amps_path)
    df_non = read_fimo_tsv(nonamps_path)

    if df_amp.empty or df_non.empty:
        print(f"Skipping {clustering}-{tool} (missing/invalid FIMO TSV).")
        continue

    any_loaded = True

    motifs = set(df_amp["motif_id"]).union(set(df_non["motif_id"]))
    total_amp_seq = int(df_amp["sequence_name"].nunique())
    total_non_seq = int(df_non["sequence_name"].nunique())

    for motif in motifs:
        amp_seq = int(df_amp.loc[df_amp["motif_id"] == motif, "sequence_name"].nunique())
        non_seq = int(df_non.loc[df_non["motif_id"] == motif, "sequence_name"].nunique())

        table = [
            [amp_seq, total_amp_seq - amp_seq],
            [non_seq, total_non_seq - non_seq],
        ]
        _, pvalue = safe_fisher(table)

        enrichment = (amp_seq + PSEUDOCOUNT) / (non_seq + PSEUDOCOUNT)

        all_results.append({
            "Clustering": clustering,
            "Tool": tool,
            "Motif": str(motif),
            "AMP_sequences_with_hit": amp_seq,
            "nonAMP_sequences_with_hit": non_seq,
            "Total_AMP_sequences": total_amp_seq,
            "Total_nonAMP_sequences": total_non_seq,
            "Enrichment_ratio": float(enrichment),
            "Fisher_pvalue": float(pvalue) if pvalue is not None else np.nan,
        })

if not any_loaded:
    die("No valid FIMO TSV files found. Run FIMO first and check paths.")

results_df = pd.DataFrame(all_results)
results_df["Enrichment_ratio"] = pd.to_numeric(results_df["Enrichment_ratio"], errors="coerce")
results_df["Fisher_pvalue"] = pd.to_numeric(results_df["Fisher_pvalue"], errors="coerce")

# =====================================================
# 2) MULTIPLE TEST CORRECTION (BH-FDR)
# =====================================================

mask_valid_p = results_df["Fisher_pvalue"].notna()
results_df["FDR_corrected"] = np.nan

if mask_valid_p.any():
    results_df.loc[mask_valid_p, "FDR_corrected"] = multipletests(
        results_df.loc[mask_valid_p, "Fisher_pvalue"].values,
        method="fdr_bh"
    )[1]

results_df["Significant_Fisher"] = (results_df["Fisher_pvalue"] < ALPHA_FISHER) & (results_df["Enrichment_ratio"] > 1)
results_df["Significant_FDR"] = (results_df["FDR_corrected"] < ALPHA_FDR) & (results_df["Enrichment_ratio"] > 1)

results_df = results_df.sort_values(
    by=["Significant_FDR", "Enrichment_ratio", "FDR_corrected"],
    ascending=[False, False, True],
)

# =====================================================
# 3) SAVE TABLES
# =====================================================

csv_all = OUT_DIR / "fimo_enrichment_all.csv"
xlsx_all = OUT_DIR / "fimo_enrichment_all.xlsx"

results_df.to_csv(csv_all, index=False)
save_excel(results_df, xlsx_all, sheet_name="All")

df_sig = results_df[results_df["Significant_FDR"]].copy()

csv_sig = OUT_DIR / "fimo_enrichment_significant_fdr.csv"
xlsx_sig = OUT_DIR / "fimo_enrichment_significant_fdr.xlsx"

df_sig.to_csv(csv_sig, index=False)
save_excel(df_sig, xlsx_sig, sheet_name="Significant_FDR")

print("\nSaved:")
print("-", csv_all)
print("-", xlsx_all)
print("-", csv_sig)
print("-", xlsx_sig)

# =====================================================
# 4) SUMMARIES
# =====================================================

summary_by_tool = results_df.groupby("Tool")["Significant_FDR"].sum().sort_values(ascending=False)
summary_by_clustering = results_df.groupby("Clustering")["Significant_FDR"].sum().sort_values(ascending=False)

summary_combo = (
    results_df.groupby(["Clustering", "Tool"])
    .agg(
        Total_motifs=("Motif", "count"),
        Significant_FDR_n=("Significant_FDR", "sum"),
        Mean_ER=("Enrichment_ratio", "mean"),
        Median_ER=("Enrichment_ratio", "median"),
    )
    .reset_index()
)

summary_combo_path = OUT_DIR / "summary_by_clustering_tool.csv"
summary_combo.to_csv(summary_combo_path, index=False)

print("\n=== Significant motifs (FDR) by Tool ===")
print(summary_by_tool)

print("\n=== Significant motifs (FDR) by Clustering ===")
print(summary_by_clustering)

print("\nSaved summary:", summary_combo_path)

# =====================================================
# 5) PLOTS
# =====================================================

plt.figure()
vals = results_df["Enrichment_ratio"].dropna().values
plt.hist(vals, bins=30)
plt.xlabel("Enrichment Ratio (AMP / nonAMP)")
plt.ylabel("Count")
plt.title("Distribution of Motif Enrichment Ratios (All motifs)")
plt.tight_layout()
plt.savefig(OUT_DIR / "enrichment_ratio_distribution_all.png")
plt.close()

top_n = 30
df_top = df_sig.sort_values("Enrichment_ratio", ascending=False).head(top_n)

if not df_top.empty:
    plt.figure(figsize=(10, 6))
    plt.bar(range(len(df_top)), df_top["Enrichment_ratio"].values)
    plt.xticks(range(len(df_top)), df_top["Motif"].astype(str).values, rotation=90)
    plt.xlabel("Motif")
    plt.ylabel("Enrichment Ratio")
    plt.title(f"Top {top_n} enriched motifs (FDR < {ALPHA_FDR})")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "top_enriched_motifs_fdr.png")
    plt.close()

if not df_sig.empty:
    data_mmseq = df_sig.loc[df_sig["Clustering"].str.contains("mmseq", case=False, na=False), "Enrichment_ratio"].dropna()
    data_cdhit = df_sig.loc[df_sig["Clustering"].str.contains("cdhit", case=False, na=False), "Enrichment_ratio"].dropna()

    plt.figure(figsize=(8, 6))
    plt.boxplot([data_mmseq.values, data_cdhit.values], labels=["MMseq", "CD-HIT"])
    plt.ylabel("Enrichment Ratio")
    plt.title("ER Distribution (Significant motifs only) — by Clustering")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "boxplot_er_significant_by_clustering.png")
    plt.close()

if not df_sig.empty:
    data_meme = df_sig.loc[df_sig["Tool"].str.contains("meme", case=False, na=False), "Enrichment_ratio"].dropna()
    data_streme = df_sig.loc[df_sig["Tool"].str.contains("streme", case=False, na=False), "Enrichment_ratio"].dropna()

    plt.figure(figsize=(8, 6))
    plt.boxplot([data_meme.values, data_streme.values], labels=["MEME", "STREME"])
    plt.ylabel("Enrichment Ratio")
    plt.title("ER Distribution (Significant motifs only) — by Tool")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "boxplot_er_significant_by_tool.png")
    plt.close()

print("\nPlots saved in:", OUT_DIR)
print("\n FIMO enrichment + FDR analysis completed.\n")