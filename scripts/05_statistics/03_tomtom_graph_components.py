"""
Graph-based motif family analysis from Tomtom similarities.
Expected inputs:
- results/statistics/fimo_enrichment/fimo_enrichment_all.csv
- results/motif_similarity/tomtom_all_vs_all/tomtom.tsv

Outputs:
- results/statistics/tomtom_graph/component_analysis_all_vs_all.csv
- results/statistics/tomtom_graph/global_summary.txt
"""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import networkx as nx
from statsmodels.stats.multitest import multipletests

# =====================================================
# CONFIG
# =====================================================

ALPHA_FDR = 0.05
ALPHA_TOMTOM_Q = 0.05
TOMTOM_SUBDIR = "tomtom_all_vs_all"

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

ENRICH_CSV = PROJECT_ROOT / "results" / "statistics" / "fimo_enrichment" / "fimo_enrichment_all.csv"
TOMTOM_TSV = PROJECT_ROOT / "results" / "motif_similarity" / TOMTOM_SUBDIR / "tomtom.tsv"

OUT_DIR = PROJECT_ROOT / "results" / "statistics" / "tomtom_graph"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_COMPONENTS = OUT_DIR / "component_analysis_all_vs_all.csv"
OUT_SUMMARY = OUT_DIR / "global_summary.txt"

# =====================================================
# HELPERS
# =====================================================

def die(msg: str) -> None:
    print(f" {msg}")
    sys.exit(1)

def load_enrichment(enrich_csv: Path) -> pd.DataFrame:
    if not enrich_csv.exists():
        die(f"Enrichment CSV not found: {enrich_csv}")

    df = pd.read_csv(enrich_csv)

    required = {"Motif", "Tool", "Clustering", "Enrichment_ratio", "Fisher_pvalue"}
    missing = required - set(df.columns)
    if missing:
        die(f"Enrichment CSV missing columns: {sorted(missing)}")

    df["Motif"] = df["Motif"].astype(str).str.strip()
    df["Tool"] = df["Tool"].astype(str).str.strip().str.lower()
    df["Clustering"] = df["Clustering"].astype(str).str.strip().str.lower()
    df["Enrichment_ratio"] = pd.to_numeric(df["Enrichment_ratio"], errors="coerce")
    df["Fisher_pvalue"] = pd.to_numeric(df["Fisher_pvalue"], errors="coerce")
    df = df.dropna(subset=["Motif", "Tool", "Clustering", "Enrichment_ratio", "Fisher_pvalue"])

    # ensure FDR_corrected exists
    if "FDR_corrected" not in df.columns or df["FDR_corrected"].isna().all():
        df = df.copy()
        df["FDR_corrected"] = np.nan
        for (cl, tl), sub in df.groupby(["Clustering", "Tool"]):
            m = sub["Fisher_pvalue"].notna()
            if m.any():
                qvals = multipletests(sub.loc[m, "Fisher_pvalue"].values, method="fdr_bh")[1]
                df.loc[sub.index[m], "FDR_corrected"] = qvals
    else:
        df["FDR_corrected"] = pd.to_numeric(df["FDR_corrected"], errors="coerce")

    df["Significant_FDR"] = (df["FDR_corrected"] < ALPHA_FDR) & (df["Enrichment_ratio"] > 1)
    return df

def load_tomtom(tomtom_tsv: Path) -> pd.DataFrame:
    if not tomtom_tsv.exists():
        die(f"Tomtom TSV not found: {tomtom_tsv}")

    tt = pd.read_csv(tomtom_tsv, sep="\t", comment="#")

    required = {"Query_ID", "Target_ID", "q-value"}
    missing = required - set(tt.columns)
    if missing:
        die(f"Tomtom TSV missing columns: {sorted(missing)}")

    tt["Query_ID"] = tt["Query_ID"].astype(str).str.strip()
    tt["Target_ID"] = tt["Target_ID"].astype(str).str.strip()
    tt["q-value"] = pd.to_numeric(tt["q-value"], errors="coerce")

    tt = tt.dropna(subset=["Query_ID", "Target_ID", "q-value"])
    tt = tt[tt["Query_ID"] != tt["Target_ID"]]
    tt = tt[tt["q-value"] < ALPHA_TOMTOM_Q].copy()
    return tt

# =====================================================
# MAIN
# =====================================================

def main():
    print("\nTOMTOM GRAPH COMPONENTS \n")
    print("Project root:", PROJECT_ROOT)
    print("Enrichment :", ENRICH_CSV)
    print("Tomtom TSV :", TOMTOM_TSV)
    print("Output dir :", OUT_DIR, "\n")

    df = load_enrichment(ENRICH_CSV)
    df_sig = df[df["Significant_FDR"]].copy()

    print("Total motifs in enrichment table:", len(df))
    print("Total significant motifs (FDR & ER>1):", len(df_sig))
    if df_sig.empty:
        die("No significant motifs found. Check enrichment outputs / thresholds.")

    sig_set = set(df_sig["Motif"].tolist())

    tt = load_tomtom(TOMTOM_TSV)
    tt = tt[tt["Query_ID"].isin(sig_set) & tt["Target_ID"].isin(sig_set)]
    print(f"Significant Tomtom edges (q < {ALPHA_TOMTOM_Q} among significant motifs): {len(tt)}")

    G = nx.Graph()
    for motif in sig_set:
        G.add_node(motif)
    for _, row in tt.iterrows():
        G.add_edge(row["Query_ID"], row["Target_ID"])

    components = list(nx.connected_components(G))
    sizes = [len(c) for c in components]
    n_singletons = sum(1 for s in sizes if s == 1)
    n_families = sum(1 for s in sizes if s >= 2)

    print("\nNumber of connected components (incl. singletons):", len(components))
    print("Families (size >= 2):", n_families)
    print("Singletons (size = 1):", n_singletons)

    component_summary = []
    for i, comp in enumerate(components, start=1):
        comp_list = sorted(list(comp))
        sub_df = df_sig[df_sig["Motif"].isin(comp_list)]

        tools = sorted(set(sub_df["Tool"]))
        clusterings = sorted(set(sub_df["Clustering"]))

        component_summary.append({
            "Component_ID": i,
            "Size": len(comp_list),
            "Motifs": ", ".join(comp_list),
            "Tools_present": ", ".join(tools),
            "Clusterings_present": ", ".join(clusterings),
            "Exclusive_MEME": set(tools) == {"meme"},
            "Exclusive_STREME": set(tools) == {"streme"},
            "Exclusive_CDHIT": set(clusterings) == {"cdhit"},
            "Exclusive_MMSEQ": set(clusterings) == {"mmseq"},
            "Mean_ER": float(sub_df["Enrichment_ratio"].mean()) if not sub_df.empty else np.nan,
        })

    component_df = pd.DataFrame(component_summary).sort_values(by=["Size", "Mean_ER"], ascending=[False, False])
    component_df.to_csv(OUT_COMPONENTS, index=False)

    with open(OUT_SUMMARY, "w", encoding="utf-8") as f:
        f.write("TOMTOM GRAPH COMPONENT SUMMARY\n\n")
        f.write(f"Enrichment input: {ENRICH_CSV}\n")
        f.write(f"Tomtom input    : {TOMTOM_TSV}\n\n")
        f.write(f"Significant motifs (FDR<{ALPHA_FDR}, ER>1): {len(df_sig)}\n")
        f.write(f"Tomtom edges kept (q<{ALPHA_TOMTOM_Q}): {len(tt)}\n\n")
        f.write(f"Components total (incl. singletons): {len(components)}\n")
        f.write(f"Families (size>=2): {n_families}\n")
        f.write(f"Singletons: {n_singletons}\n\n")
        f.write("Exclusivity counts (components):\n")
        f.write(f"- Exclusive MEME   : {int(component_df['Exclusive_MEME'].sum())}\n")
        f.write(f"- Exclusive STREME : {int(component_df['Exclusive_STREME'].sum())}\n")
        f.write(f"- Exclusive CD-HIT : {int(component_df['Exclusive_CDHIT'].sum())}\n")
        f.write(f"- Exclusive MMseq  : {int(component_df['Exclusive_MMSEQ'].sum())}\n")

    print("\nSaved:")
    print("-", OUT_COMPONENTS)
    print("-", OUT_SUMMARY)
    print("\n Done.\n")

if __name__ == "__main__":
    main()