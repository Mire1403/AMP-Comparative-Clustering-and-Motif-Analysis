"""
FIMO results: reporting + robustness

INPUT (preferred):
- results/statistics/fimo_enrichment/fimo_enrichment_all.csv

Fallback:
- results/statistics/fimo_enrichment/fimo_enrichment_all.xlsx

OUTPUT:
- results/statistics/fimo_reporting/
  - summary_all_motifs.csv
  - summary_significant_only.csv
  - fimo_significant_only.csv
  - boxplot_ER_significant_by_clustering.png
  - boxplot_ER_significant_by_tool.png
  - robustness_consensus_similarity.csv
  - robustness_consensus_similarity.txt
"""
from __future__ import annotations

from pathlib import Path
import sys
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests

# =====================================================
# CONFIG
# =====================================================

ALPHA_FDR = 0.05
ER_HIGH = 2.0
ER_EXTREME = 3.0
MAX_MISMATCH = 2  # hamming threshold

plt.rcParams["figure.dpi"] = 150
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["font.size"] = 12

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

ENRICH_DIR = PROJECT_ROOT / "results" / "statistics" / "fimo_enrichment"
OUT_DIR = PROJECT_ROOT / "results" / "statistics" / "fimo_reporting"
OUT_DIR.mkdir(parents=True, exist_ok=True)

MOTIF_DIR = PROJECT_ROOT / "results" / "motif_discovery"
MOTIF_FILES = {
    ("cdhit", "meme"):   MOTIF_DIR / "meme_cdhit" / "meme.txt",
    ("mmseq", "meme"):   MOTIF_DIR / "meme_mmseq" / "meme.txt",
    ("cdhit", "streme"): MOTIF_DIR / "streme_cdhit" / "streme.txt",
    ("mmseq", "streme"): MOTIF_DIR / "streme_mmseq" / "streme.txt",
}

# =====================================================
# IO
# =====================================================

def die(msg: str) -> None:
    print(f"{msg}")
    sys.exit(1)

def load_enrichment_table(enrich_dir: Path) -> pd.DataFrame:
    csv_path = enrich_dir / "fimo_enrichment_all.csv"
    xlsx_path = enrich_dir / "fimo_enrichment_all.xlsx"

    if csv_path.exists():
        return pd.read_csv(csv_path)
    if xlsx_path.exists():
        return pd.read_excel(xlsx_path, engine="openpyxl")

    die(f"Could not find enrichment table at {csv_path} or {xlsx_path}")

def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    rename_map = {}
    if "FDR" in df.columns and "FDR_corrected" not in df.columns:
        rename_map["FDR"] = "FDR_corrected"
    df = df.rename(columns=rename_map)

    for col in ["Clustering", "Tool", "Motif"]:
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip()

    df["Clustering"] = df["Clustering"].str.lower()
    df["Tool"] = df["Tool"].str.lower()

    for col in ["Enrichment_ratio", "Fisher_pvalue", "FDR_corrected"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["Clustering", "Tool", "Motif", "Enrichment_ratio", "Fisher_pvalue"])
    return df

def ensure_fdr(df: pd.DataFrame) -> pd.DataFrame:
    if "FDR_corrected" in df.columns and df["FDR_corrected"].notna().any():
        return df

    df = df.copy()
    df["FDR_corrected"] = np.nan
    for (cl, tl), sub in df.groupby(["Clustering", "Tool"], dropna=False):
        m = sub["Fisher_pvalue"].notna()
        if m.any():
            qvals = multipletests(sub.loc[m, "Fisher_pvalue"].values, method="fdr_bh")[1]
            df.loc[sub.index[m], "FDR_corrected"] = qvals
    return df

# =====================================================
# REPORTING
# =====================================================

def build_summaries(df: pd.DataFrame):
    df = df.copy()
    df["Significant_FDR"] = (df["FDR_corrected"] < ALPHA_FDR) & (df["Enrichment_ratio"] > 1)
    df["is_high"] = df["Significant_FDR"] & (df["Enrichment_ratio"] > ER_HIGH)
    df["is_extreme"] = df["Significant_FDR"] & (df["Enrichment_ratio"] > ER_EXTREME)

    summary_all = []
    for (cluster, tool), group in df.groupby(["Clustering", "Tool"]):
        total = len(group)
        sig = int(group["Significant_FDR"].sum())
        high = int(group["is_high"].sum())
        extreme = int(group["is_extreme"].sum())
        summary_all.append({
            "Clustering": cluster,
            "Tool": tool,
            "Total_motifs": total,
            "Significant_FDR_n": sig,
            "Significant_FDR_%": (sig / total * 100) if total else 0.0,
            "Highly_enriched_%": (high / total * 100) if total else 0.0,
            "Extremely_enriched_%": (extreme / total * 100) if total else 0.0,
            "Mean_ER": float(group["Enrichment_ratio"].mean()),
            "Median_ER": float(group["Enrichment_ratio"].median()),
        })

    summary_all_df = pd.DataFrame(summary_all).sort_values(
        by=["Significant_FDR_%", "Mean_ER"],
        ascending=[False, False]
    )

    df_sig = df[df["Significant_FDR"]].copy()

    summary_sig = []
    for (cluster, tool), group in df_sig.groupby(["Clustering", "Tool"]):
        total = len(group)
        high = int((group["Enrichment_ratio"] > ER_HIGH).sum())
        extreme = int((group["Enrichment_ratio"] > ER_EXTREME).sum())
        summary_sig.append({
            "Clustering": cluster,
            "Tool": tool,
            "Total_significant_motifs": total,
            "Highly_enriched_%": (high / total * 100) if total else 0.0,
            "Extremely_enriched_%": (extreme / total * 100) if total else 0.0,
            "Mean_ER": float(group["Enrichment_ratio"].mean()) if total else np.nan,
            "Median_ER": float(group["Enrichment_ratio"].median()) if total else np.nan,
        })

    summary_sig_df = pd.DataFrame(summary_sig).sort_values(
        by=["Total_significant_motifs", "Mean_ER"],
        ascending=[False, False]
    )

    return summary_all_df, df_sig, summary_sig_df

def plot_boxplots(df_sig: pd.DataFrame, out_dir: Path) -> None:
    if df_sig.empty:
        print("No significant motifs to plot.")
        return

    # clustering boxplot (cdhit vs mmseq)
    data = []
    labels = []
    for key, label in [("mmseq", "MMseq"), ("cdhit", "CD-HIT")]:
        v = df_sig.loc[df_sig["Clustering"] == key, "Enrichment_ratio"].dropna().values
        if len(v):
            data.append(v)
            labels.append(label)

    if data:
        plt.figure(figsize=(9, 7))
        plt.boxplot(data, labels=labels)
        plt.ylabel("Enrichment Ratio")
        plt.title("ER Distribution (Significant motifs only) — by Clustering")
        plt.tight_layout()
        plt.savefig(out_dir / "boxplot_ER_significant_by_clustering.png")
        plt.close()

    # tool boxplot (meme vs streme)
    data = []
    labels = []
    for key, label in [("meme", "MEME"), ("streme", "STREME")]:
        v = df_sig.loc[df_sig["Tool"] == key, "Enrichment_ratio"].dropna().values
        if len(v):
            data.append(v)
            labels.append(label)

    if data:
        plt.figure(figsize=(9, 7))
        plt.boxplot(data, labels=labels)
        plt.ylabel("Enrichment Ratio")
        plt.title("ER Distribution (Significant motifs only) — by Tool")
        plt.tight_layout()
        plt.savefig(out_dir / "boxplot_ER_significant_by_tool.png")
        plt.close()

# =====================================================
# ROBUSTNESS (CONSENSUS)
# =====================================================

def infer_alphabet(lines: list[str]) -> str | None:
    for ln in lines[:200]:
        if ln.startswith("ALPHABET="):
            alph = ln.split("=", 1)[1].strip().replace(" ", "")
            return alph or None
    for i, ln in enumerate(lines[:500]):
        if "Background letter frequencies" in ln:
            for ln2 in lines[i:i+10]:
                toks = ln2.strip().split()
                syms = []
                for t in toks[::2]:
                    if len(t) == 1 and t.isalpha():
                        syms.append(t)
                if len(syms) >= 4:
                    return "".join(syms)
    return None

def extract_consensus(motif_file: Path) -> dict[str, str]:
    consensus_dict: dict[str, str] = {}
    if not motif_file.exists():
        return consensus_dict

    lines = motif_file.read_text(errors="ignore").splitlines(True)
    alphabet = infer_alphabet(lines) or "ACDEFGHIKLMNPQRSTVWY"

    i = 0
    while i < len(lines):
        if lines[i].startswith("MOTIF"):
            parts = lines[i].strip().split()
            if len(parts) < 2:
                i += 1
                continue
            motif_id = parts[1]
            i += 1

            while i < len(lines) and "letter-probability matrix" not in lines[i]:
                i += 1
            if i >= len(lines):
                break

            width_match = re.search(r"w=\s*(\d+)", lines[i])
            if not width_match:
                i += 1
                continue
            width = int(width_match.group(1))
            i += 1

            matrix = []
            for _ in range(width):
                if i >= len(lines):
                    break
                row = lines[i].strip().split()
                try:
                    probs = list(map(float, row[:len(alphabet)]))
                except Exception:
                    probs = []
                if len(probs) == len(alphabet):
                    matrix.append(probs)
                i += 1

            if len(matrix) != width:
                continue

            cons = []
            for position in matrix:
                cons.append(alphabet[int(np.argmax(position))] if len(position) else "X")
            consensus_dict[motif_id] = "".join(cons)
        else:
            i += 1

    return consensus_dict

def hamming(s1: str, s2: str) -> int | None:
    if len(s1) != len(s2):
        return None
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def count_similar(set1: list[str], set2: list[str], max_mismatch: int) -> int:
    matches = 0
    for s1 in set1:
        for s2 in set2:
            d = hamming(s1, s2)
            if d is not None and d <= max_mismatch:
                matches += 1
                break
    return matches

def build_group_consensus(df_sig: pd.DataFrame, clustering: str, tool: str, cons_dict: dict[str, str]) -> list[str]:
    subset = df_sig[(df_sig["Clustering"] == clustering) & (df_sig["Tool"] == tool)]
    out = []
    for motif in subset["Motif"].astype(str).tolist():
        if motif in cons_dict:
            out.append(cons_dict[motif])
    return out

def run_robustness(df_sig: pd.DataFrame, out_dir: Path) -> None:
    if df_sig.empty:
        print("No significant motifs -> skipping robustness.")
        return

    cons_maps = {k: extract_consensus(p) for k, p in MOTIF_FILES.items()}
    if any(len(v) == 0 for v in cons_maps.values()):
        print("⚠️  Some motif files missing/empty. Robustness will be partial.")

    meme_cdhit = build_group_consensus(df_sig, "cdhit", "meme", cons_maps.get(("cdhit", "meme"), {}))
    meme_mmseq = build_group_consensus(df_sig, "mmseq", "meme", cons_maps.get(("mmseq", "meme"), {}))
    streme_cdhit = build_group_consensus(df_sig, "cdhit", "streme", cons_maps.get(("cdhit", "streme"), {}))
    streme_mmseq = build_group_consensus(df_sig, "mmseq", "streme", cons_maps.get(("mmseq", "streme"), {}))

    results = []

    def add_result(name: str, a: list[str], b: list[str]):
        if len(a) == 0 or len(b) == 0:
            results.append({
                "Comparison": name,
                "SetA_n": len(a),
                "SetB_n": len(b),
                "Similar_A_to_B_n": np.nan,
                "Similar_B_to_A_n": np.nan,
                "Max_mismatch": MAX_MISMATCH,
            })
            return
        sim_a = count_similar(a, b, MAX_MISMATCH)
        sim_b = count_similar(b, a, MAX_MISMATCH)
        results.append({
            "Comparison": name,
            "SetA_n": len(a),
            "SetB_n": len(b),
            "Similar_A_to_B_n": sim_a,
            "Similar_B_to_A_n": sim_b,
            "Max_mismatch": MAX_MISMATCH,
        })

    add_result("MEME: CD-HIT vs MMseq", meme_cdhit, meme_mmseq)
    add_result("STREME: CD-HIT vs MMseq", streme_cdhit, streme_mmseq)
    add_result("CD-HIT: MEME vs STREME", meme_cdhit, streme_cdhit)
    add_result("MMseq: MEME vs STREME", meme_mmseq, streme_mmseq)

    rob_df = pd.DataFrame(results)
    rob_csv = out_dir / "robustness_consensus_similarity.csv"
    rob_txt = out_dir / "robustness_consensus_similarity.txt"
    rob_df.to_csv(rob_csv, index=False)

    with open(rob_txt, "w", encoding="utf-8") as f:
        f.write("ROBUSTNESS (consensus-level similarity)\n")
        f.write(f"Rule: similar if Hamming distance <= {MAX_MISMATCH} (same length)\n\n")
        f.write(rob_df.to_string(index=False))
        f.write("\n")

    print("\n=== ROBUSTNESS (consensus-level) ===")
    print(rob_df.to_string(index=False))
    print("\nSaved robustness:")
    print("-", rob_csv)
    print("-", rob_txt)

# =====================================================
# MAIN
# =====================================================

def main():
    print("\nFIMO REPORTING + ROBUSTNESS \n")
    print("Project root:", PROJECT_ROOT)
    print("Enrichment dir:", ENRICH_DIR)
    print("Output dir:", OUT_DIR)

    df = load_enrichment_table(ENRICH_DIR)
    df = normalize_columns(df)
    df = ensure_fdr(df)

    print("\nLoaded rows:", df.shape[0])

    summary_all_df, df_sig, summary_sig_df = build_summaries(df)

    summary_all_path = OUT_DIR / "summary_all_motifs.csv"
    summary_sig_path = OUT_DIR / "summary_significant_only.csv"
    sig_path = OUT_DIR / "fimo_significant_only.csv"

    summary_all_df.to_csv(summary_all_path, index=False)
    summary_sig_df.to_csv(summary_sig_path, index=False)
    df_sig.to_csv(sig_path, index=False)

    print("\nSaved reporting tables:")
    print("-", summary_all_path)
    print("-", summary_sig_path)
    print("-", sig_path)

    if not df.empty:
        pct_by_clustering = df.groupby("Clustering")["Significant_FDR"].mean() * 100
        pct_by_tool = df.groupby("Tool")["Significant_FDR"].mean() * 100
        print("\n% Significant (FDR) by Clustering:")
        print(pct_by_clustering.sort_values(ascending=False))
        print("\n% Significant (FDR) by Tool:")
        print(pct_by_tool.sort_values(ascending=False))

    plot_boxplots(df_sig, OUT_DIR)
    print("\nSaved plots:")
    print("-", OUT_DIR / "boxplot_ER_significant_by_clustering.png")
    print("-", OUT_DIR / "boxplot_ER_significant_by_tool.png")

    run_robustness(df_sig, OUT_DIR)

    print("\n Done.\n")

if __name__ == "__main__":
    main()