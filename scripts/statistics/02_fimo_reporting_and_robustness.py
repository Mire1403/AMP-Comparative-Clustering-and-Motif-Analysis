"""
FIMO results: reporting + robustness

INPUT (preferred):
- results/statistics/fimo_enrichment/fimo_enrichment_all.csv
  (from your enrichment+FDR pipeline)

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
MAX_MISMATCH = 2  # Hamming distance threshold for "similar" motifs (consensus-level)

# Figures
plt.rcParams["figure.dpi"] = 150
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["font.size"] = 12


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

ENRICH_DIR = PROJECT_ROOT / "results" / "statistics" / "fimo_enrichment"
OUT_DIR = PROJECT_ROOT / "results" / "statistics" / "fimo_reporting"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Motif discovery files (optional for robustness)
MOTIF_DIR = PROJECT_ROOT / "results" / "motif_discovery"
MOTIF_FILES = {
    ("cdhit80", "meme"):   MOTIF_DIR / "meme_cdhit80" / "meme.txt",
    ("mmseq80", "meme"):   MOTIF_DIR / "meme_mmseq80" / "meme.txt",
    ("cdhit80", "streme"): MOTIF_DIR / "streme_cdhit80" / "streme.txt",
    ("mmseq80", "streme"): MOTIF_DIR / "streme_mmseq80" / "streme.txt",
}


# =====================================================
# IO HELPERS
# =====================================================

def load_enrichment_table(enrich_dir: Path) -> pd.DataFrame:
    csv_path = enrich_dir / "fimo_enrichment_all.csv"
    xlsx_path = enrich_dir / "fimo_enrichment_all.xlsx"

    if csv_path.exists():
        df = pd.read_csv(csv_path)
        return df

    if xlsx_path.exists():
        df = pd.read_excel(xlsx_path, engine="openpyxl")
        return df

    print("❌ Could not find enrichment table.")
    print("Expected:", csv_path, "or", xlsx_path)
    sys.exit(1)


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    # Standardize naming if needed
    rename_map = {}
    if "FDR" in df.columns and "FDR_corrected" not in df.columns:
        rename_map["FDR"] = "FDR_corrected"
    df = df.rename(columns=rename_map)

    # Clean strings
    for col in ["Clustering", "Tool", "Motif"]:
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip()

    # Lower-case keys for grouping
    if "Clustering" in df.columns:
        df["Clustering"] = df["Clustering"].str.lower()
    if "Tool" in df.columns:
        df["Tool"] = df["Tool"].str.lower()

    # Numeric conversions
    for col in ["Enrichment_ratio", "Fisher_pvalue", "FDR_corrected"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["Clustering", "Tool", "Motif", "Enrichment_ratio", "Fisher_pvalue"])
    return df


def ensure_fdr(df: pd.DataFrame, alpha_fdr: float) -> pd.DataFrame:
    """
    Ensure df has FDR_corrected.
    If missing, compute BH-FDR per (Clustering, Tool) group.
    """
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

def build_summaries(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Returns:
    - summary_all: metrics over all motifs
    - df_sig: significant only (FDR < ALPHA_FDR and ER > 1)
    - summary_sig: metrics over significant only
    """
    df = df.copy()

    df["Significant_FDR"] = (df["FDR_corrected"] < ALPHA_FDR) & (df["Enrichment_ratio"] > 1)

    # thresholds (defined on significant status)
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

    # metrics over significant only
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

    # by clustering
    clust_vals = {}
    for key in ["mmseq80", "cdhit80"]:
        v = df_sig.loc[df_sig["Clustering"] == key, "Enrichment_ratio"].dropna().values
        if len(v):
            clust_vals[key] = v

    if len(clust_vals) >= 1:
        plt.figure(figsize=(9, 7))
        labels = []
        data = []
        for k in ["mmseq80", "cdhit80"]:
            if k in clust_vals:
                labels.append("MMseq" if "mmseq" in k else "CD-HIT")
                data.append(clust_vals[k])

        plt.boxplot(data, labels=labels)
        plt.ylabel("Enrichment Ratio")
        plt.title("ER Distribution (Significant motifs only) — by Clustering")
        plt.tight_layout()
        plt.savefig(out_dir / "boxplot_ER_significant_by_clustering.png")
        plt.close()

    # by tool
    tool_vals = {}
    for key in ["meme", "streme"]:
        v = df_sig.loc[df_sig["Tool"] == key, "Enrichment_ratio"].dropna().values
        if len(v):
            tool_vals[key] = v

    if len(tool_vals) >= 1:
        plt.figure(figsize=(9, 7))
        labels = []
        data = []
        for k in ["meme", "streme"]:
            if k in tool_vals:
                labels.append(k.upper())
                data.append(tool_vals[k])

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
    """
    Try to infer alphabet order from MEME/STREME motif file.
    Returns a string of symbols in the matrix column order (best effort).
    """
    # Common MEME header line: "ALPHABET= ACGT" or "ALPHABET= ACDE..."
    for ln in lines[:200]:
        if ln.startswith("ALPHABET="):
            alph = ln.split("=", 1)[1].strip().replace(" ", "")
            if alph:
                return alph

    # Another clue: "Background letter frequencies" often lists letters with frequencies.
    # Example: "A 0.07 C 0.01 D 0.05 ..."
    for i, ln in enumerate(lines[:500]):
        if "Background letter frequencies" in ln:
            # scan next few lines for tokens like "A 0.0"
            for ln2 in lines[i:i+10]:
                toks = ln2.strip().split()
                # collect alternating symbol/freq pairs
                syms = []
                for t in toks[::2]:
                    if len(t) == 1 and t.isalpha():
                        syms.append(t)
                if len(syms) >= 4:
                    return "".join(syms)

    return None


def extract_consensus(motif_file: Path) -> dict[str, str]:
    """
    Parse MEME/STREME motif file and build a simple consensus per motif
    by taking argmax per position.

    NOTE: The alphabet order is inferred from file when possible.
    If not found, falls back to a standard protein alphabet.
    """
    consensus_dict: dict[str, str] = {}

    if not motif_file.exists():
        return consensus_dict

    lines = motif_file.read_text(errors="ignore").splitlines(True)
    alphabet = infer_alphabet(lines)

    if alphabet is None:
        # Fallback: protein alphabet (common), but we warn
        alphabet = "ACDEFGHIKLMNPQRSTVWY"
        print(f"⚠️  Alphabet not found in {motif_file.name}. Falling back to protein alphabet order.")

    i = 0
    while i < len(lines):
        if lines[i].startswith("MOTIF"):
            parts = lines[i].strip().split()
            if len(parts) < 2:
                i += 1
                continue
            motif_id = parts[1]
            i += 1

            # Skip until matrix header
            while i < len(lines) and "letter-probability matrix" not in lines[i]:
                i += 1
            if i >= len(lines):
                break

            # Extract width w=
            width_match = re.search(r"w=\s*(\d+)", lines[i])
            if not width_match:
                i += 1
                continue
            width = int(width_match.group(1))
            i += 1

            matrix = []
            # read width rows
            for _ in range(width):
                if i >= len(lines):
                    break
                row = lines[i].strip().split()
                # Some files may include extra columns; keep only len(alphabet)
                try:
                    probs = list(map(float, row[:len(alphabet)]))
                except Exception:
                    probs = []
                if len(probs) == len(alphabet):
                    matrix.append(probs)
                i += 1

            if len(matrix) != width:
                continue

            # Build consensus
            cons = []
            for position in matrix:
                max_index = int(np.argmax(position))
                if max_index < len(alphabet):
                    cons.append(alphabet[max_index])
                else:
                    cons.append("X")
            consensus_dict[motif_id] = "".join(cons)
        else:
            i += 1

    return consensus_dict


def hamming(s1: str, s2: str) -> int | None:
    if len(s1) != len(s2):
        return None
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def count_similar(set1: list[str], set2: list[str], max_mismatch: int) -> int:
    """
    Count how many sequences in set1 have at least one partner in set2
    with Hamming distance <= max_mismatch (same length).
    """
    matches = 0
    for s1 in set1:
        found = False
        for s2 in set2:
            d = hamming(s1, s2)
            if d is not None and d <= max_mismatch:
                found = True
                break
        if found:
            matches += 1
    return matches


def build_group_consensus(df_sig: pd.DataFrame, clustering: str, tool: str, cons_dict: dict[str, str]) -> list[str]:
    subset = df_sig[(df_sig["Clustering"] == clustering) & (df_sig["Tool"] == tool)]
    seqs = []
    for motif in subset["Motif"].astype(str).tolist():
        if motif in cons_dict:
            seqs.append(cons_dict[motif])
    return seqs


def run_robustness(df_sig: pd.DataFrame, out_dir: Path) -> None:
    """
    Robustness comparisons at consensus level:
    - MEME: CD-HIT vs MMseq
    - STREME: CD-HIT vs MMseq
    - CD-HIT: MEME vs STREME
    - MMseq: MEME vs STREME
    """
    if df_sig.empty:
        print("No significant motifs -> skipping robustness.")
        return

    # load consensus dicts (optional)
    cons_maps = {}
    missing_any = False
    for key, path in MOTIF_FILES.items():
        cons = extract_consensus(path)
        cons_maps[key] = cons
        if not cons:
            missing_any = True

    if missing_any:
        print("⚠️  Some motif files missing/empty. Robustness will be partial (only where files exist).")

    # build sets
    meme_cdhit = build_group_consensus(df_sig, "cdhit80", "meme", cons_maps.get(("cdhit80", "meme"), {}))
    meme_mmseq = build_group_consensus(df_sig, "mmseq80", "meme", cons_maps.get(("mmseq80", "meme"), {}))

    streme_cdhit = build_group_consensus(df_sig, "cdhit80", "streme", cons_maps.get(("cdhit80", "streme"), {}))
    streme_mmseq = build_group_consensus(df_sig, "mmseq80", "streme", cons_maps.get(("mmseq80", "streme"), {}))

    results = []

    def add_result(name: str, a: list[str], b: list[str]):
        if len(a) == 0 or len(b) == 0:
            results.append({
                "Comparison": name,
                "SetA_n": len(a),
                "SetB_n": len(b),
                "Similar_A_to_B_n": np.nan,
                "Similar_B_to_A_n": np.nan,
                "Max_mismatch": MAX_MISMATCH
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
            "Max_mismatch": MAX_MISMATCH
        })

    add_result("MEME: CD-HIT vs MMseq", meme_cdhit, meme_mmseq)
    add_result("STREME: CD-HIT vs MMseq", streme_cdhit, streme_mmseq)
    add_result("CD-HIT: MEME vs STREME", meme_cdhit, streme_cdhit)
    add_result("MMseq: MEME vs STREME", meme_mmseq, streme_mmseq)

    rob_df = pd.DataFrame(results)
    rob_csv = out_dir / "robustness_consensus_similarity.csv"
    rob_txt = out_dir / "robustness_consensus_similarity.txt"
    rob_df.to_csv(rob_csv, index=False)

    # also write a readable txt
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
    print("\nFIMO REPORTING + ROBUSTNESS (repo-friendly)\n")
    print("Project root:", PROJECT_ROOT)
    print("Enrichment dir:", ENRICH_DIR)
    print("Output dir:", OUT_DIR)

    df = load_enrichment_table(ENRICH_DIR)
    df = normalize_columns(df)
    df = ensure_fdr(df, ALPHA_FDR)

    print("\nLoaded rows:", df.shape[0])

    summary_all_df, df_sig, summary_sig_df = build_summaries(df)

    # Save outputs
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

    # Basic comparisons
    if not df.empty:
        pct_by_clustering = df.groupby("Clustering")["Significant_FDR"].mean() * 100
        pct_by_tool = df.groupby("Tool")["Significant_FDR"].mean() * 100
        print("\n% Significant (FDR) by Clustering:")
        print(pct_by_clustering.sort_values(ascending=False))
        print("\n% Significant (FDR) by Tool:")
        print(pct_by_tool.sort_values(ascending=False))

    # Plots
    plot_boxplots(df_sig, OUT_DIR)
    print("\nSaved plots:")
    print("-", OUT_DIR / "boxplot_ER_significant_by_clustering.png")
    print("-", OUT_DIR / "boxplot_ER_significant_by_tool.png")

    # Robustness
    run_robustness(df_sig, OUT_DIR)

    print("\n✅ Done.\n")


if __name__ == "__main__":
    main()