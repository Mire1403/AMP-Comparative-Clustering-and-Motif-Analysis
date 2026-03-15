from __future__ import annotations

from pathlib import Path
import sys
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

print("\nFIMO ENRICHMENT + FDR PIPELINE\n")

# =====================================================
# REPO ROOT
# =====================================================

def find_repo_root(start: Path) -> Path:
    markers = [".git", "README.md", "requirements.txt", "pyproject.toml", "environment.yml"]
    start = start.resolve()
    for parent in [start] + list(start.parents):
        if any((parent / m).exists() for m in markers):
            return parent
    return start.parents[2]

PROJECT_ROOT = find_repo_root(Path(__file__).parent)

# =====================================================
# PATH CONFIG
# =====================================================

FIMO_BASE = PROJECT_ROOT / "results" / "05_motif_scanning"
OUT_DIR = PROJECT_ROOT / "results" / "statistics" / "03_fimo_enrichment"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CLUSTER_DIR = PROJECT_ROOT / "results" / "02_clustering"
BG_DIR = PROJECT_ROOT / "results" / "03_background"

FASTA_TOTALS = {
    "cdhit_amps": CLUSTER_DIR / "AMP_MASTER_cdhit.fasta",
    "mmseq_amps": CLUSTER_DIR / "AMP_MASTER_mmseq_rep_seq.fasta",
    "cdhit_nonamps": BG_DIR / "nonamp_cdhit_10x_len_kr_matched_progressive_random.fasta",
    "mmseq_nonamps": BG_DIR / "nonamp_mmseq_10x_len_kr_matched_progressive_random.fasta",
}

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
# CONFIG
# =====================================================

ALPHA_FISHER = 0.05
ALPHA_FDR = 0.05
MIN_AMP_HITS = 1

# warning heuristics
WARN_IF_NONAMP_ZERO_FRAC_GT = 0.50
WARN_IF_SHARED_MOTIFS_LT = 3

# =====================================================
# HELPERS
# =====================================================

def die(msg: str) -> None:
    print(msg)
    sys.exit(1)

def clean_motif_name(x: str) -> str:
    """
    Removes numeric prefixes typical of STREME:
    '24-DEGEMTEEEKK' -> 'DEGEMTEEEKK'
    '8_ALKAIAKAAKKLL' -> 'ALKAIAKAAKKLL'
    """
    x = str(x).strip()
    x = re.sub(r"^\d+[-_]", "", x)
    return x.strip()

def count_fasta_sequences(fasta_path: Path) -> int:
    if not fasta_path.exists():
        die(f"FASTA not found (needed to compute totals): {fasta_path}")
    n = 0
    with open(fasta_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n

def safe_read_fimo(path: Path) -> tuple[pd.DataFrame, str]:
    """
    Robust FIMO reader.
    Returns:
      dataframe with columns: motif_name, sequence_name
      status string
    """
    if not path.exists():
        return pd.DataFrame(columns=["motif_name", "sequence_name"]), f"missing_file: {path}"

    try:
        df = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
    except Exception as e:
        return pd.DataFrame(columns=["motif_name", "sequence_name"]), f"read_error: {e}"

    if df.empty:
        return pd.DataFrame(columns=["motif_name", "sequence_name"]), "empty_table"

    if "sequence_name" not in df.columns:
        return pd.DataFrame(columns=["motif_name", "sequence_name"]), "missing_sequence_name"

    motif_col = None

    # Prefer motif_id when informative; otherwise motif_alt_id
    if "motif_id" in df.columns:
        motif_id_values = df["motif_id"].dropna().astype(str)
        informative = motif_id_values.str.contains(r"[A-Za-z]").any()
        if informative:
            motif_col = "motif_id"

    if motif_col is None and "motif_alt_id" in df.columns:
        motif_alt_values = df["motif_alt_id"].dropna().astype(str)
        if len(motif_alt_values) > 0:
            motif_col = "motif_alt_id"

    if motif_col is None:
        return pd.DataFrame(columns=["motif_name", "sequence_name"]), "missing_motif_columns"

    out = df[[motif_col, "sequence_name"]].copy()
    out.columns = ["motif_raw", "sequence_name"]
    out["motif_name"] = out["motif_raw"].astype(str).map(clean_motif_name)
    out["sequence_name"] = out["sequence_name"].astype(str).str.strip()
    out = out.dropna(subset=["motif_name", "sequence_name"])

    return out[["motif_name", "sequence_name"]], f"ok_using_{motif_col}"

def safe_fisher(amp_hit: int, amp_total: int, non_hit: int, non_total: int):
    table = [
        [amp_hit, amp_total - amp_hit],
        [non_hit, non_total - non_hit],
    ]
    try:
        oddsratio, pvalue = fisher_exact(table)
        return float(oddsratio), float(pvalue)
    except Exception:
        return np.nan, np.nan

def save_excel(df: pd.DataFrame, path: Path, sheet_name: str = "Results") -> None:
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name=sheet_name)

def save_pair(df: pd.DataFrame, stem: Path, sheet_name: str = "Results") -> None:
    df.to_csv(str(stem) + ".csv", index=False)
    save_excel(df, Path(str(stem) + ".xlsx"), sheet_name=sheet_name)

# =====================================================
# PRECOMPUTE TRUE TOTALS FROM FASTA
# =====================================================

TOTALS = {k: count_fasta_sequences(v) for k, v in FASTA_TOTALS.items()}

print("FASTA totals (true denominators):")
for k, v in TOTALS.items():
    print(f"  {k}: {v:,}  ({FASTA_TOTALS[k].name})")

# =====================================================
# 1) BUILD ENRICHMENT TABLE — per pipeline + diagnostics
# =====================================================

all_results: list[pd.DataFrame] = []
diagnostic_rows: list[dict] = []
any_loaded = False

for (clustering, tool), files in STRUCTURE.items():
    print(f"\n--- {clustering}-{tool} ---")
    print("AMP FIMO   :", files["amps"])
    print("nonAMP FIMO:", files["nonamps"])

    df_amp, amp_status = safe_read_fimo(files["amps"])
    df_non, non_status = safe_read_fimo(files["nonamps"])

    print("AMP status   :", amp_status)
    print("nonAMP status:", non_status)
    print("AMP rows     :", len(df_amp))
    print("nonAMP rows  :", len(df_non))

    # AMP is mandatory
    if df_amp.empty:
        print(f"  Skipping {clustering}-{tool}: AMP FIMO missing/invalid/empty.")
        diagnostic_rows.append({
            "Clustering": clustering,
            "Tool": tool,
            "AMP_status": amp_status,
            "nonAMP_status": non_status,
            "AMP_rows": len(df_amp),
            "nonAMP_rows": len(df_non),
            "AMP_motifs": 0,
            "nonAMP_motifs": 0,
            "Shared_motifs": 0,
            "Frac_nonAMP_zero": np.nan,
            "Suspicious": True,
            "Comment": "AMP FIMO unavailable",
        })
        continue

    any_loaded = True

    if clustering == "cdhit":
        total_amp = TOTALS["cdhit_amps"]
        total_non = TOTALS["cdhit_nonamps"]
    else:
        total_amp = TOTALS["mmseq_amps"]
        total_non = TOTALS["mmseq_nonamps"]

    amp_motifs = set(df_amp["motif_name"].dropna().astype(str))
    non_motifs = set(df_non["motif_name"].dropna().astype(str)) if not df_non.empty else set()
    shared_motifs = amp_motifs & non_motifs

    print("AMP motifs unique   :", len(amp_motifs))
    print("nonAMP motifs unique:", len(non_motifs))
    print("Shared motifs       :", len(shared_motifs))

    if len(shared_motifs) > 0:
        print("Shared examples     :", list(sorted(shared_motifs))[:10])
    else:
        print("Shared examples     : []")

    motifs = sorted(amp_motifs | non_motifs)
    pipe_results = []

    for motif in motifs:
        amp_hit = int(df_amp.loc[df_amp["motif_name"] == motif, "sequence_name"].nunique())
        non_hit = int(df_non.loc[df_non["motif_name"] == motif, "sequence_name"].nunique()) if not df_non.empty else 0

        amp_hit = min(amp_hit, total_amp)
        non_hit = min(non_hit, total_non)

        if amp_hit < MIN_AMP_HITS:
            continue

        amp_rate = amp_hit / total_amp
        non_rate = non_hit / total_non

        er = np.inf if non_rate == 0 else amp_rate / non_rate
        oddsratio, pvalue = safe_fisher(amp_hit, total_amp, non_hit, total_non)

        pipe_results.append({
            "Clustering": clustering,
            "Tool": tool,
            "Motif": motif,
            "AMP_sequences_with_hit": amp_hit,
            "nonAMP_sequences_with_hit": non_hit,
            "Total_AMP_sequences": total_amp,
            "Total_nonAMP_sequences": total_non,
            "AMP_hit_rate": round(amp_rate, 6),
            "nonAMP_hit_rate": round(non_rate, 6),
            "Enrichment_ratio": er,
            "ER_is_infinite": bool(non_rate == 0 and amp_hit > 0),
            "Fisher_oddsratio": oddsratio,
            "Fisher_pvalue": pvalue,
        })

    if not pipe_results:
        print(f"  No motifs passed filters for {clustering}-{tool}.")
        diagnostic_rows.append({
            "Clustering": clustering,
            "Tool": tool,
            "AMP_status": amp_status,
            "nonAMP_status": non_status,
            "AMP_rows": len(df_amp),
            "nonAMP_rows": len(df_non),
            "AMP_motifs": len(amp_motifs),
            "nonAMP_motifs": len(non_motifs),
            "Shared_motifs": len(shared_motifs),
            "Frac_nonAMP_zero": np.nan,
            "Suspicious": True,
            "Comment": "No motifs passed filters",
        })
        continue

    pipe_df = pd.DataFrame(pipe_results)

    valid_p = pipe_df["Fisher_pvalue"].notna()
    pipe_df["FDR_corrected"] = np.nan
    if valid_p.any():
        pipe_df.loc[valid_p, "FDR_corrected"] = multipletests(
            pipe_df.loc[valid_p, "Fisher_pvalue"].values,
            method="fdr_bh"
        )[1]

    pipe_df["Significant_Fisher"] = (
        (pipe_df["Fisher_pvalue"] < ALPHA_FISHER) &
        (pipe_df["Enrichment_ratio"] > 1)
    )
    pipe_df["Significant_FDR"] = (
        (pipe_df["FDR_corrected"] < ALPHA_FDR) &
        (pipe_df["Enrichment_ratio"] > 1)
    )

    pipe_df = pipe_df.sort_values(
        ["Significant_FDR", "ER_is_infinite", "Enrichment_ratio", "FDR_corrected", "AMP_sequences_with_hit"],
        ascending=[False, False, False, True, False],
    )

    frac_nonamp_zero = (pipe_df["nonAMP_sequences_with_hit"] == 0).mean()
    suspicious = False
    comments = []

    if len(shared_motifs) < WARN_IF_SHARED_MOTIFS_LT:
        suspicious = True
        comments.append("very low shared motif names between AMP and nonAMP")

    if frac_nonamp_zero > WARN_IF_NONAMP_ZERO_FRAC_GT:
        suspicious = True
        comments.append("high fraction of motifs with nonAMP=0")

    if df_non.empty:
        suspicious = True
        comments.append("nonAMP FIMO empty or unparsable")

    diagnostic_rows.append({
        "Clustering": clustering,
        "Tool": tool,
        "AMP_status": amp_status,
        "nonAMP_status": non_status,
        "AMP_rows": len(df_amp),
        "nonAMP_rows": len(df_non),
        "AMP_motifs": len(amp_motifs),
        "nonAMP_motifs": len(non_motifs),
        "Shared_motifs": len(shared_motifs),
        "Frac_nonAMP_zero": round(float(frac_nonamp_zero), 4),
        "Suspicious": suspicious,
        "Comment": "; ".join(comments) if comments else "ok",
    })

    if suspicious:
        print("  WARNING:", "; ".join(comments))
        print("  This pipeline may have a parsing/path/background problem.")

    tag = f"{clustering}_{tool}"
    save_pair(pipe_df, OUT_DIR / f"fimo_enrichment_{tag}", sheet_name=tag)

    sig_pipe_df = pipe_df[pipe_df["Significant_FDR"]].copy()
    save_pair(sig_pipe_df, OUT_DIR / f"fimo_enrichment_{tag}_significant", sheet_name=f"{tag}_sig")

    print(
        f"  {clustering}-{tool}: {len(pipe_df)} motifs, "
        f"{int(pipe_df['Significant_FDR'].sum())} significant, "
        f"{pipe_df['ER_is_infinite'].sum()} with infinite ER"
    )

if not any_loaded:
    die("No valid AMP FIMO TSV files found. Run FIMO first and check paths.")

diagnostic_df = pd.DataFrame(diagnostic_rows)
save_pair(diagnostic_df, OUT_DIR / "fimo_enrichment_diagnostics", sheet_name="Diagnostics")

# =====================================================
# 2) GLOBAL TABLE — FDR across all pipelines combined
# =====================================================

if not all_results and "pipe_df" in locals():
    die("No enrichment tables were produced.")
elif not all_results:
    # collect saved per-loop tables if any_loaded but nothing appended
    produced = list(OUT_DIR.glob("fimo_enrichment_*_*.csv"))
    if not produced:
        die("No enrichment tables were produced.")
else:
    pass

# rebuild from newly produced pipe_dfs in memory
per_pipeline_csvs = sorted([
    p for p in OUT_DIR.glob("fimo_enrichment_*_*.csv")
    if "significant" not in p.name and "all" not in p.name and "diagnostics" not in p.name
])

if not per_pipeline_csvs:
    die("No per-pipeline enrichment CSVs found to build global table.")

all_results = [pd.read_csv(p) for p in per_pipeline_csvs]
results_df = pd.concat(all_results, ignore_index=True)

valid_p = results_df["Fisher_pvalue"].notna()
results_df["FDR_corrected_global"] = np.nan
if valid_p.any():
    results_df.loc[valid_p, "FDR_corrected_global"] = multipletests(
        results_df.loc[valid_p, "Fisher_pvalue"].values,
        method="fdr_bh"
    )[1]

results_df["Significant_FDR_global"] = (
    (results_df["FDR_corrected_global"] < ALPHA_FDR) &
    (results_df["Enrichment_ratio"] > 1)
)

results_df = results_df.sort_values(
    ["Significant_FDR", "ER_is_infinite", "Enrichment_ratio", "FDR_corrected", "AMP_sequences_with_hit"],
    ascending=[False, False, False, True, False],
)

save_pair(results_df, OUT_DIR / "fimo_enrichment_all", sheet_name="All")

df_sig = results_df[results_df["Significant_FDR"]].copy()
save_pair(df_sig, OUT_DIR / "fimo_enrichment_significant_fdr", sheet_name="Significant_FDR")

print(
    f"\nGlobal table: {len(results_df)} motifs total, "
    f"{int(results_df['Significant_FDR'].sum())} significant (per-pipeline FDR), "
    f"{int(results_df['ER_is_infinite'].sum())} with infinite ER."
)

# =====================================================
# 3) SUMMARY TABLE
# =====================================================

summary_combo = (
    results_df.groupby(["Clustering", "Tool"])
    .agg(
        Total_motifs=("Motif", "count"),
        Significant_FDR_n=("Significant_FDR", "sum"),
        Mean_ER=("Enrichment_ratio", lambda x: x.replace(np.inf, np.nan).mean()),
        Median_ER=("Enrichment_ratio", lambda x: x.replace(np.inf, np.nan).median()),
        N_infinite_ER=("ER_is_infinite", "sum"),
        Median_AMP_hit_rate=("AMP_hit_rate", "median"),
        Median_nonAMP_hit_rate=("nonAMP_hit_rate", "median"),
    )
    .reset_index()
)

summary_combo["Significant_FDR_%"] = (
    100 * summary_combo["Significant_FDR_n"] / summary_combo["Total_motifs"]
).round(1)

save_pair(summary_combo, OUT_DIR / "summary_by_clustering_tool", sheet_name="Summary")
print("\nSummary by pipeline:")
print(summary_combo.to_string(index=False))

# =====================================================
# 4) PLOTS
# =====================================================

finite_er = results_df["Enrichment_ratio"].replace(np.inf, np.nan).dropna()
if not finite_er.empty:
    plt.figure(figsize=(8, 5))
    plt.hist(finite_er.values, bins=40, edgecolor="white")
    plt.xlabel("Enrichment Ratio (AMP rate / nonAMP rate)")
    plt.ylabel("Count")
    plt.title("Distribution of Motif Enrichment Ratios\n(infinite ER excluded)")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "enrichment_ratio_distribution_all.png", dpi=200)
    plt.close()

top_n = 30
df_top = df_sig.copy()
df_top["ER_plot"] = df_top["Enrichment_ratio"].replace(np.inf, np.nan)
df_top = df_top.dropna(subset=["ER_plot"]).sort_values("ER_plot", ascending=False).head(top_n)

if not df_top.empty:
    labels = [f"{m}\n({cl}-{tl})" for m, cl, tl in zip(df_top["Motif"], df_top["Clustering"], df_top["Tool"])]
    plt.figure(figsize=(12, 6))
    plt.bar(range(len(df_top)), df_top["ER_plot"].values)
    plt.xticks(range(len(df_top)), labels, rotation=90, fontsize=7)
    plt.ylabel("Enrichment Ratio")
    plt.title(f"Top {top_n} enriched motifs (finite ER only)")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "top_enriched_motifs_fdr.png", dpi=200)
    plt.close()

if not df_sig.empty:
    groups = {
        cl: df_sig.loc[df_sig["Clustering"] == cl, "AMP_hit_rate"].dropna().values
        for cl in df_sig["Clustering"].unique()
    }
    groups = {k: v for k, v in groups.items() if len(v) > 0}
    if groups:
        plt.figure(figsize=(7, 5))
        plt.boxplot(list(groups.values()), labels=list(groups.keys()), showfliers=False)
        plt.ylabel("AMP hit rate")
        plt.title("AMP hit-rate distribution — significant motifs by clustering")
        plt.tight_layout()
        plt.savefig(OUT_DIR / "boxplot_amp_hit_rate_significant_by_clustering.png", dpi=200)
        plt.close()

if not df_sig.empty:
    groups = {
        tl: df_sig.loc[df_sig["Tool"] == tl, "AMP_hit_rate"].dropna().values
        for tl in df_sig["Tool"].unique()
    }
    groups = {k: v for k, v in groups.items() if len(v) > 0}
    if groups:
        plt.figure(figsize=(7, 5))
        plt.boxplot(list(groups.values()), labels=list(groups.keys()), showfliers=False)
        plt.ylabel("AMP hit rate")
        plt.title("AMP hit-rate distribution — significant motifs by tool")
        plt.tight_layout()
        plt.savefig(OUT_DIR / "boxplot_amp_hit_rate_significant_by_tool.png", dpi=200)
        plt.close()

# =====================================================
# 5) PRINT SAVED FILES
# =====================================================

print("\nAll saved files:")
for f in sorted(OUT_DIR.iterdir()):
    print(" -", f.name)

print("\nFIMO enrichment + FDR analysis completed.\n")
print("=" * 60)
print("CHECK FIRST:")
print("=" * 60)
print("""
1. fimo_enrichment_diagnostics.csv
   → sanity checks for each pipeline

2. fimo_enrichment_all.csv
   → complete enrichment table

3. fimo_enrichment_significant_fdr.csv
   → significant motifs only

4. summary_by_clustering_tool.csv
   → quick overview of suspicious infinite-ER patterns
""")