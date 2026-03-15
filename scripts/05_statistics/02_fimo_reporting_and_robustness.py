from __future__ import annotations

from pathlib import Path
import sys
import re
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =====================================================
# CONFIG
# =====================================================

ALPHA_FDR = 0.05
ER_MIN = 2.0
TOP_N = 25
MAX_MISMATCH = 2
USE_LOG10_ER_FOR_PLOTS = True
NONAMP_BACKGROUND_FACTOR = 10

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

ENRICH_DIR = PROJECT_ROOT / "results" / "statistics" / "03_fimo_enrichment"
OUT_DIR = PROJECT_ROOT / "results" / "statistics" / "04_fimo_reporting"
OUT_DIR.mkdir(parents=True, exist_ok=True)

MOTIF_DIR = PROJECT_ROOT / "results" / "motif_discovery"
MOTIF_FILES = {
    ("cdhit", "meme"):   MOTIF_DIR / "meme_cdhit" / "meme.txt",
    ("mmseq", "meme"):   MOTIF_DIR / "meme_mmseq" / "meme.txt",
    ("cdhit", "streme"): MOTIF_DIR / "streme_cdhit" / "streme.txt",
    ("mmseq", "streme"): MOTIF_DIR / "streme_mmseq" / "streme.txt",
}

FIMO_DIRS = {
    ("cdhit", "meme"):   PROJECT_ROOT / "results" / "05_motif_scanning" / "cdhit" / "meme" / "fimo_amps",
    ("mmseq", "meme"):   PROJECT_ROOT / "results" / "05_motif_scanning" / "mmseq" / "meme" / "fimo_amps",
    ("cdhit", "streme"): PROJECT_ROOT / "results" / "05_motif_scanning" / "cdhit" / "streme" / "fimo_amps",
    ("mmseq", "streme"): PROJECT_ROOT / "results" / "05_motif_scanning" / "mmseq" / "streme" / "fimo_amps",
}

# =====================================================
# IO / HELPERS
# =====================================================

def die(msg: str) -> None:
    print(msg)
    sys.exit(1)

def clean_motif_name(x: str) -> str:
    """
    Remove numeric prefixes typically added by STREME:
    '24-DEGEMTEEEKK' -> 'DEGEMTEEEKK'
    """
    x = str(x).strip()
    x = re.sub(r"^\d+[-_]", "", x)
    return x.strip()

def load_enrichment_table(enrich_dir: Path) -> pd.DataFrame:
    csv_path = enrich_dir / "fimo_enrichment_all.csv"
    xlsx_path = enrich_dir / "fimo_enrichment_all.xlsx"

    if csv_path.exists():
        print(f"Using enrichment table: {csv_path}")
        return pd.read_csv(csv_path)

    if xlsx_path.exists():
        print(f"Using enrichment table: {xlsx_path}")
        return pd.read_excel(xlsx_path, engine="openpyxl")

    # fallback: pick first file containing 'enrichment'
    csv_candidates = sorted(enrich_dir.glob("*enrichment*.csv"))
    xlsx_candidates = sorted(enrich_dir.glob("*enrichment*.xlsx"))

    if csv_candidates:
        print(f"Using fallback enrichment CSV: {csv_candidates[0]}")
        return pd.read_csv(csv_candidates[0])

    if xlsx_candidates:
        print(f"Using fallback enrichment XLSX: {xlsx_candidates[0]}")
        return pd.read_excel(xlsx_candidates[0], engine="openpyxl")

    die(f"Could not find enrichment table in {enrich_dir}")
    return pd.DataFrame()

def load_diagnostics(enrich_dir: Path) -> pd.DataFrame:
    path = enrich_dir / "fimo_enrichment_diagnostics.csv"
    if path.exists():
        try:
            return pd.read_csv(path)
        except Exception as e:
            print(f"Warning: could not read diagnostics file {path}: {e}")
    return pd.DataFrame()

def choose_motif_column(df: pd.DataFrame) -> str | None:
    """
    Prefer motif_id when it looks informative (contains letters),
    otherwise fallback to motif_alt_id, then motif.
    """
    if "motif_id" in df.columns:
        vals = df["motif_id"].dropna().astype(str)
        if not vals.empty and vals.str.contains(r"[A-Za-z]").any():
            return "motif_id"

    if "motif_alt_id" in df.columns:
        vals = df["motif_alt_id"].dropna().astype(str)
        if not vals.empty:
            return "motif_alt_id"

    if "motif" in df.columns:
        return "motif"

    return None

# =====================================================
# RAW FIMO HIT COUNTS
# =====================================================

def load_fimo_hit_counts(fimo_dirs: dict) -> pd.DataFrame:
    """
    Reads raw FIMO hit tables for AMP scans and returns per motif:
      - Total_hits
      - AMP_sequences_hit_raw
      - AMP_total_sequences
      - AMP_hit_frequency_pct
    """
    records = []

    for (clustering, tool), fimo_dir in fimo_dirs.items():
        fimo_tsv = fimo_dir / "fimo.tsv"

        if not fimo_tsv.exists():
            print(f"  Warning: fimo.tsv not found at {fimo_tsv} — skipping hit counts for {clustering}-{tool}")
            continue

        try:
            fimo = pd.read_csv(fimo_tsv, sep="\t", comment="#", low_memory=False)
        except Exception as e:
            print(f"  Warning: could not read {fimo_tsv}: {e}")
            continue

        if fimo.empty:
            print(f"  Warning: empty fimo.tsv at {fimo_tsv}")
            continue

        fimo.columns = [c.strip().lower().replace("-", "_") for c in fimo.columns]

        motif_col = choose_motif_column(fimo)
        seq_col = next((c for c in ["sequence_name", "sequence"] if c in fimo.columns), None)

        if motif_col is None or seq_col is None:
            print(f"  Warning: could not find motif/sequence columns in {fimo_tsv} — skipping.")
            continue

        fimo["_motif_clean"] = fimo[motif_col].astype(str).map(clean_motif_name)
        fimo[seq_col] = fimo[seq_col].astype(str).str.strip()
        total_amp_seqs = fimo[seq_col].nunique()

        for motif_id, grp in fimo.groupby("_motif_clean", dropna=True):
            amp_seqs_hit = grp[seq_col].nunique()
            records.append({
                "Clustering": clustering,
                "Tool": tool,
                "Motif": motif_id,
                "Total_hits": len(grp),
                "AMP_sequences_hit_raw": amp_seqs_hit,
                "AMP_total_sequences": total_amp_seqs,
                "AMP_hit_frequency_pct": round(100.0 * amp_seqs_hit / total_amp_seqs, 2)
                                         if total_amp_seqs > 0 else np.nan,
            })

    if not records:
        return pd.DataFrame(columns=[
            "Clustering", "Tool", "Motif",
            "Total_hits", "AMP_sequences_hit_raw",
            "AMP_total_sequences", "AMP_hit_frequency_pct",
        ])

    return pd.DataFrame(records)

# =====================================================
# NORMALIZE COLUMNS
# =====================================================

def normalize_columns(df: pd.DataFrame, hit_counts: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    if "FDR" in df.columns and "FDR_corrected" not in df.columns:
        df = df.rename(columns={"FDR": "FDR_corrected"})

    required = [
        "Clustering", "Tool", "Motif",
        "AMP_sequences_with_hit", "nonAMP_sequences_with_hit",
        "Enrichment_ratio", "Fisher_pvalue", "FDR_corrected",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        die(f"Missing required columns in enrichment table: {missing}")

    for col in ["Clustering", "Tool", "Motif"]:
        df[col] = df[col].astype(str).str.strip()

    df["Clustering"] = df["Clustering"].str.lower()
    df["Tool"] = df["Tool"].str.lower()

    df["Motif_raw"] = df["Motif"].astype(str)
    df["Motif"] = df["Motif"].astype(str).map(clean_motif_name)

    numeric_cols = [
        "AMP_sequences_with_hit",
        "nonAMP_sequences_with_hit",
        "Enrichment_ratio",
        "Fisher_pvalue",
        "FDR_corrected",
    ]
    if "Fisher_oddsratio" in df.columns:
        numeric_cols.append("Fisher_oddsratio")
    if "Total_AMP_sequences" in df.columns:
        numeric_cols.append("Total_AMP_sequences")
    if "Total_nonAMP_sequences" in df.columns:
        numeric_cols.append("Total_nonAMP_sequences")

    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["ER_is_infinite"] = (
        df["nonAMP_sequences_with_hit"].fillna(0).eq(0) &
        df["AMP_sequences_with_hit"].fillna(0).gt(0)
    )

    df["ER_corrected"] = df["Enrichment_ratio"].copy()

    if "Total_AMP_sequences" in df.columns and "Total_nonAMP_sequences" in df.columns:
        df["AMP_hit_rate"] = df["AMP_sequences_with_hit"] / df["Total_AMP_sequences"]
        df["nonAMP_hit_rate"] = df["nonAMP_sequences_with_hit"] / df["Total_nonAMP_sequences"]
        df["ER_corrected"] = df["AMP_hit_rate"] / df["nonAMP_hit_rate"].replace(0, np.nan)
    else:
        df["AMP_hit_rate"] = np.nan
        df["nonAMP_hit_rate"] = np.nan

    if "Fisher_oddsratio" in df.columns:
        df["OddsRatio_is_infinite"] = np.isinf(df["Fisher_oddsratio"])
        df["Fisher_oddsratio_plot"] = df["Fisher_oddsratio"].replace([np.inf, -np.inf], np.nan)
    else:
        df["OddsRatio_is_infinite"] = False
        df["Fisher_oddsratio_plot"] = np.nan

    df["Significant_FDR"] = (df["FDR_corrected"] < ALPHA_FDR) & (df["Enrichment_ratio"] > ER_MIN)
    df["MinusLog10_FDR"] = -np.log10(df["FDR_corrected"].clip(lower=1e-300))
    df["AMP_nonAMP_hit_diff"] = df["AMP_sequences_with_hit"] - df["nonAMP_sequences_with_hit"]

    if USE_LOG10_ER_FOR_PLOTS:
        finite_er = df["Enrichment_ratio"].replace([np.inf, -np.inf], np.nan)
        df["log10_ER"] = np.log10(finite_er.clip(lower=1e-12))
    else:
        df["log10_ER"] = df["Enrichment_ratio"].replace([np.inf, -np.inf], np.nan)

    df = df.dropna(subset=["Clustering", "Tool", "Motif", "FDR_corrected"])

    if not hit_counts.empty:
        hc = hit_counts.copy()
        hc["Clustering"] = hc["Clustering"].str.lower()
        hc["Tool"] = hc["Tool"].str.lower()
        hc["Motif"] = hc["Motif"].astype(str).map(clean_motif_name)

        df = df.merge(
            hc[[
                "Clustering", "Tool", "Motif",
                "Total_hits", "AMP_total_sequences", "AMP_hit_frequency_pct"
            ]],
            on=["Clustering", "Tool", "Motif"],
            how="left",
        )
    else:
        df["Total_hits"] = np.nan
        df["AMP_total_sequences"] = np.nan
        df["AMP_hit_frequency_pct"] = np.nan

    return df

# =====================================================
# MOTIF CATEGORIES
# =====================================================

def classify_motif(seq: str) -> str:
    seq = str(seq)
    L = len(seq)

    if L == 0:
        return "Other"

    frac_c = seq.count("C") / L
    frac_g = seq.count("G") / L
    frac_kr = (seq.count("K") + seq.count("R")) / L
    frac_hydrophobic = sum(seq.count(x) for x in "AILMFWVY") / L

    if frac_c >= 0.25:
        return "Cys-rich"
    if frac_g >= 0.30:
        return "Gly-rich"
    if frac_kr >= 0.35:
        return "Cationic-rich"
    if frac_hydrophobic >= 0.60:
        return "Hydrophobic/signal-like"
    return "Other"

# =====================================================
# SUMMARIES
# =====================================================

def build_summaries(df: pd.DataFrame):
    df = df.copy()
    df["Motif_class"] = df["Motif"].map(classify_motif)

    summary_all = (
        df.groupby(["Clustering", "Tool"], dropna=False)
        .agg(
            Total_motifs=("Motif", "count"),
            Significant_FDR_n=("Significant_FDR", "sum"),
            Mean_ER=("Enrichment_ratio", lambda x: x.replace([np.inf, -np.inf], np.nan).mean()),
            Median_ER=("Enrichment_ratio", lambda x: x.replace([np.inf, -np.inf], np.nan).median()),
            N_infinite_ER=("ER_is_infinite", "sum"),
            Mean_log10_ER=("log10_ER", "mean"),
        )
        .reset_index()
    )

    summary_all["Significant_FDR_%"] = (
        100 * summary_all["Significant_FDR_n"] / summary_all["Total_motifs"]
    )
    summary_all["Frac_infinite_ER"] = (
        summary_all["N_infinite_ER"] / summary_all["Total_motifs"]
    ).round(3)

    df_sig = df[df["Significant_FDR"]].copy()

    agg_sig = dict(
        Total_significant_motifs=("Motif", "count"),
        Mean_ER=("Enrichment_ratio", lambda x: x.replace([np.inf, -np.inf], np.nan).mean()),
        Median_ER=("Enrichment_ratio", lambda x: x.replace([np.inf, -np.inf], np.nan).median()),
        N_infinite_ER=("ER_is_infinite", "sum"),
        Mean_hits_AMP=("AMP_sequences_with_hit", "mean"),
        Mean_hits_nonAMP=("nonAMP_sequences_with_hit", "mean"),
    )

    if "Total_hits" in df_sig.columns:
        agg_sig["Mean_Total_hits"] = ("Total_hits", "mean")

    if "AMP_hit_frequency_pct" in df_sig.columns:
        agg_sig["Mean_AMP_hit_frequency_pct"] = ("AMP_hit_frequency_pct", "mean")
        agg_sig["Median_AMP_hit_frequency_pct"] = ("AMP_hit_frequency_pct", "median")

    summary_sig = (
        df_sig.groupby(["Clustering", "Tool"], dropna=False)
        .agg(**agg_sig)
        .reset_index()
    )

    class_summary = (
        df_sig.groupby(["Clustering", "Tool", "Motif_class"], dropna=False)
        .size()
        .reset_index(name="Count")
        .sort_values(["Clustering", "Tool", "Count"], ascending=[True, True, False])
    )

    return summary_all, df_sig, summary_sig, class_summary

# =====================================================
# TOP TABLES
# =====================================================

def build_top_tables(df_sig: pd.DataFrame):
    if df_sig.empty:
        return pd.DataFrame(), pd.DataFrame()

    cols = [
        "Clustering", "Tool", "Motif", "Motif_raw", "Motif_class",
        "AMP_sequences_with_hit", "nonAMP_sequences_with_hit",
        "Enrichment_ratio", "ER_corrected", "ER_is_infinite",
        "FDR_corrected", "MinusLog10_FDR",
        "OddsRatio_is_infinite",
        "Total_hits", "AMP_total_sequences", "AMP_hit_frequency_pct",
    ]
    cols = [c for c in cols if c in df_sig.columns]

    sort_cols = ["ER_is_infinite", "AMP_sequences_with_hit", "Enrichment_ratio", "FDR_corrected"]
    sort_asc = [False, False, False, True]

    top_global = (
        df_sig.sort_values(by=sort_cols, ascending=sort_asc)[cols]
        .head(TOP_N)
        .copy()
    )

    top_per_combo = (
        df_sig.sort_values(by=sort_cols, ascending=sort_asc)
        .groupby(["Clustering", "Tool"], as_index=False, group_keys=False)
        .head(TOP_N)[cols]
        .copy()
    )

    return top_global, top_per_combo

# =====================================================
# OVERLAP / ROBUSTNESS
# =====================================================

def infer_alphabet(lines: list[str]) -> str | None:
    for ln in lines[:200]:
        if ln.startswith("ALPHABET="):
            alph = ln.split("=", 1)[1].strip().replace(" ", "")
            return alph or None
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

            motif_id = clean_motif_name(parts[1])
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
                cons.append(alphabet[int(np.argmax(position))] if position else "X")
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

def build_overlap_tables(df_sig: pd.DataFrame, out_dir: Path):
    if df_sig.empty:
        return pd.DataFrame(), pd.DataFrame()

    combo_to_set = {}
    for (cl, tl), sub in df_sig.groupby(["Clustering", "Tool"]):
        combo_to_set[(cl, tl)] = set(sub["Motif"].astype(str))

    exact_rows = []
    combos = sorted(combo_to_set.keys())

    for a, b in combinations(combos, 2):
        set_a = combo_to_set[a]
        set_b = combo_to_set[b]
        inter = set_a & set_b
        union = set_a | set_b
        jaccard = len(inter) / len(union) if union else np.nan

        exact_rows.append({
            "Comparison_A": f"{a[0]}_{a[1]}",
            "Comparison_B": f"{b[0]}_{b[1]}",
            "SetA_n": len(set_a),
            "SetB_n": len(set_b),
            "Intersection_n": len(inter),
            "Jaccard": jaccard,
            "Shared_motifs": ";".join(sorted(inter)) if inter else "",
        })

    exact_df = pd.DataFrame(exact_rows)
    exact_df.to_csv(out_dir / "motif_overlap_exact.csv", index=False)

    cons_maps = {k: extract_consensus(p) for k, p in MOTIF_FILES.items()}
    sim_rows = []

    for a, b in combinations(combos, 2):
        motifs_a = df_sig[
            (df_sig["Clustering"] == a[0]) & (df_sig["Tool"] == a[1])
        ]["Motif"].astype(str).tolist()

        motifs_b = df_sig[
            (df_sig["Clustering"] == b[0]) & (df_sig["Tool"] == b[1])
        ]["Motif"].astype(str).tolist()

        cons_a = [x for x in [cons_maps.get(a, {}).get(m) for m in motifs_a] if x]
        cons_b = [x for x in [cons_maps.get(b, {}).get(m) for m in motifs_b] if x]

        sim_rows.append({
            "Comparison_A": f"{a[0]}_{a[1]}",
            "Comparison_B": f"{b[0]}_{b[1]}",
            "SetA_n": len(cons_a),
            "SetB_n": len(cons_b),
            "Similar_A_to_B_n": count_similar(cons_a, cons_b, MAX_MISMATCH) if cons_a and cons_b else np.nan,
            "Similar_B_to_A_n": count_similar(cons_b, cons_a, MAX_MISMATCH) if cons_a and cons_b else np.nan,
            "Max_mismatch": MAX_MISMATCH,
        })

    sim_df = pd.DataFrame(sim_rows)
    sim_df.to_csv(out_dir / "motif_overlap_consensus_similarity.csv", index=False)

    return exact_df, sim_df

# =====================================================
# PLOTS
# =====================================================

def plot_counts_bar(summary_all: pd.DataFrame, out_dir: Path) -> None:
    if summary_all.empty:
        return

    xlabels = [f"{r.Clustering}\n{r.Tool}" for _, r in summary_all.iterrows()]
    vals = summary_all["Significant_FDR_%"].values
    abs_n = summary_all["Significant_FDR_n"].values

    plt.figure(figsize=(8, 6))
    bars = plt.bar(range(len(vals)), vals)

    plt.xticks(range(len(vals)), xlabels)
    plt.ylabel("% significant motifs")
    plt.title("Significant enriched motifs by pipeline")

    y_offset = max(vals) * 0.02 if len(vals) and np.isfinite(vals).any() else 1.0

    for bar, n in zip(bars, abs_n):
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + y_offset,
            f"n={int(n)}",
            ha="center",
            va="bottom",
            fontsize=10
        )

    plt.ylim(0, max(vals) * 1.18 if len(vals) and np.isfinite(vals).any() else 100)
    plt.tight_layout()
    plt.savefig(out_dir / "significant_motifs_by_pipeline_pct.png")
    plt.close()

def plot_top_motifs(df_sig: pd.DataFrame, out_dir: Path) -> None:
    if df_sig.empty:
        return

    top = df_sig.copy()

    top["ER_is_infinite"] = (
        top["nonAMP_sequences_with_hit"].fillna(0).eq(0) &
        top["AMP_sequences_with_hit"].fillna(0).gt(0)
    )

    finite_mask = ~top["ER_is_infinite"]
    finite_log_vals = top.loc[finite_mask, "log10_ER"].replace([np.inf, -np.inf], np.nan).dropna()

    finite_max = finite_log_vals.max() if not finite_log_vals.empty else 1.0
    infinite_display_height = finite_max + 0.35

    top["plot_height"] = top["log10_ER"]
    top.loc[top["ER_is_infinite"], "plot_height"] = infinite_display_height

    top = top.sort_values(
        by=["ER_is_infinite", "AMP_sequences_with_hit", "Enrichment_ratio", "FDR_corrected"],
        ascending=[False, False, False, True]
    ).head(TOP_N).copy()

    labels = [f"{m}\n({cl}-{tl})" for m, cl, tl in zip(top["Motif"], top["Clustering"], top["Tool"])]

    plt.figure(figsize=(14, 8))
    bars = plt.bar(range(len(top)), top["plot_height"].values)

    for bar, is_inf in zip(bars, top["ER_is_infinite"]):
        if is_inf:
            bar.set_hatch("//")

    plt.xticks(range(len(top)), labels, rotation=90)
    plt.ylabel("log10(Enrichment Ratio)")
    plt.title(f"Top {TOP_N} significant motifs")

    if not finite_log_vals.empty:
        plt.axhline(finite_max, linestyle="--", linewidth=1)

    for i, (_, row) in enumerate(top.iterrows()):
        if row["ER_is_infinite"]:
            plt.text(
                i,
                row["plot_height"] + 0.04,
                "∞",
                ha="center",
                va="bottom",
                fontsize=12,
                fontweight="bold"
            )

    plt.ylim(0, top["plot_height"].max() * 1.12 if not top.empty else 1)
    plt.tight_layout()
    plt.savefig(out_dir / "top_significant_motifs_log10ER_fixed.png")
    plt.close()

def plot_scatter_er_fdr(df_sig: pd.DataFrame, out_dir: Path) -> None:
    if df_sig.empty:
        return

    plot_df = df_sig.replace([np.inf, -np.inf], np.nan).dropna(subset=["log10_ER", "MinusLog10_FDR"])
    if plot_df.empty:
        return

    plt.figure(figsize=(8, 6))
    for (cl, tl), sub in plot_df.groupby(["Clustering", "Tool"]):
        plt.scatter(sub["log10_ER"], sub["MinusLog10_FDR"], label=f"{cl}-{tl}", alpha=0.75)

    plt.xlabel("log10(Enrichment Ratio)")
    plt.ylabel("-log10(FDR)")
    plt.title("Significant motifs: effect size vs significance")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / "scatter_log10ER_vs_minuslog10FDR.png")
    plt.close()

def plot_class_distribution(class_summary: pd.DataFrame, out_dir: Path) -> None:
    if class_summary.empty:
        return

    pivot = class_summary.pivot_table(
        index=["Clustering", "Tool"],
        columns="Motif_class",
        values="Count",
        fill_value=0,
    )

    pivot.plot(kind="bar", stacked=True, figsize=(10, 6))
    plt.ylabel("Number of significant motifs")
    plt.title("Motif class composition among significant motifs")
    plt.tight_layout()
    plt.savefig(out_dir / "significant_motif_classes_stacked.png")
    plt.close()

def plot_hit_frequency(df_sig: pd.DataFrame, out_dir: Path) -> None:
    if "AMP_hit_frequency_pct" not in df_sig.columns:
        return

    top = (
        df_sig.dropna(subset=["AMP_hit_frequency_pct"])
        .sort_values("AMP_hit_frequency_pct", ascending=False)
        .head(TOP_N)
        .copy()
    )

    if top.empty:
        return

    labels = [f"{m}\n({cl}-{tl})" for m, cl, tl in zip(top["Motif"], top["Clustering"], top["Tool"])]

    plt.figure(figsize=(12, 7))
    plt.bar(range(len(top)), top["AMP_hit_frequency_pct"].values)
    plt.xticks(range(len(top)), labels, rotation=90)
    plt.ylabel("% AMP sequences with ≥1 hit")
    plt.title(f"Top {TOP_N} significant motifs — AMP hit frequency")
    plt.tight_layout()
    plt.savefig(out_dir / "top_motifs_amp_hit_frequency_pct.png")
    plt.close()

def plot_hit_frequency_boxplot(df_sig: pd.DataFrame, out_dir: Path) -> None:
    if df_sig.empty or "AMP_hit_frequency_pct" not in df_sig.columns:
        return

    plot_df = df_sig.dropna(subset=["AMP_hit_frequency_pct"]).copy()
    if plot_df.empty:
        return

    plot_df["Pipeline"] = plot_df["Clustering"] + "-" + plot_df["Tool"]

    pipeline_order = ["cdhit-meme", "mmseq-meme", "cdhit-streme", "mmseq-streme"]
    available = [p for p in pipeline_order if p in plot_df["Pipeline"].unique()]
    data = [plot_df.loc[plot_df["Pipeline"] == p, "AMP_hit_frequency_pct"].values for p in available]

    if not data:
        return

    plt.figure(figsize=(9, 6))
    plt.boxplot(data, labels=available, showfliers=True)
    plt.ylabel("% AMP sequences with ≥1 hit")
    plt.title("Distribution of AMP hit frequency among significant motifs")
    plt.tight_layout()
    plt.savefig(out_dir / "significant_motifs_amp_hit_frequency_boxplot.png")
    plt.close()

# =====================================================
# MAIN
# =====================================================

def main():
    print("\nFIMO REPORTING + INTERPRETATION\n")
    print("Project root:", PROJECT_ROOT)
    print("Enrichment dir:", ENRICH_DIR)
    print("Output dir:", OUT_DIR)
    print(f"Background factor (nonAMP/AMP): {NONAMP_BACKGROUND_FACTOR}×")

    diag = load_diagnostics(ENRICH_DIR)
    if not diag.empty:
        suspicious = diag[diag["Suspicious"] == True].copy()
        if not suspicious.empty:
            print("\nWARNING: suspicious pipelines detected in enrichment diagnostics:")
            cols = [c for c in ["Clustering", "Tool", "Shared_motifs", "Frac_nonAMP_zero", "Comment"] if c in suspicious.columns]
            print(suspicious[cols].to_string(index=False))
            print("\nInterpret reporting outputs with caution.\n")

    df = load_enrichment_table(ENRICH_DIR)

    print("\nLoading raw FIMO hit counts...")
    hit_counts = load_fimo_hit_counts(FIMO_DIRS)
    if hit_counts.empty:
        print("  No fimo.tsv files found — Total_hits and AMP_hit_frequency_pct will be NaN.")
    else:
        print(f"  Loaded hit counts for {len(hit_counts)} motif entries across all pipelines.")

    df = normalize_columns(df, hit_counts)
    print("\nLoaded rows:", df.shape[0])

    cleaned_all_path = OUT_DIR / "fimo_enrichment_all_cleaned.csv"
    df.to_csv(cleaned_all_path, index=False)

    summary_all, df_sig, summary_sig, class_summary = build_summaries(df)
    top_global, top_per_combo = build_top_tables(df_sig)

    paths = {
        "summary_all_motifs": OUT_DIR / "summary_all_motifs.csv",
        "summary_significant_only": OUT_DIR / "summary_significant_only.csv",
        "fimo_significant_only_cleaned": OUT_DIR / "fimo_significant_only_cleaned.csv",
        "significant_motif_class_summary": OUT_DIR / "significant_motif_class_summary.csv",
        "top_global_significant_motifs": OUT_DIR / "top_global_significant_motifs.csv",
        "top_significant_motifs_by_pipeline": OUT_DIR / "top_significant_motifs_by_pipeline.csv",
    }

    summary_all.to_csv(paths["summary_all_motifs"], index=False)
    summary_sig.to_csv(paths["summary_significant_only"], index=False)
    df_sig.to_csv(paths["fimo_significant_only_cleaned"], index=False)
    class_summary.to_csv(paths["significant_motif_class_summary"], index=False)
    top_global.to_csv(paths["top_global_significant_motifs"], index=False)
    top_per_combo.to_csv(paths["top_significant_motifs_by_pipeline"], index=False)

    print("\nSaved reporting tables:")
    for p in [cleaned_all_path] + list(paths.values()):
        print("-", p)

    exact_overlap_df, sim_overlap_df = build_overlap_tables(df_sig, OUT_DIR)
    if not exact_overlap_df.empty:
        print("\nExact overlap table saved.")
    if not sim_overlap_df.empty:
        print("Consensus-similarity overlap table saved.")

    plot_counts_bar(summary_all, OUT_DIR)
    plot_top_motifs(df_sig, OUT_DIR)
    plot_scatter_er_fdr(df_sig, OUT_DIR)
    plot_class_distribution(class_summary, OUT_DIR)
    plot_hit_frequency(df_sig, OUT_DIR)
    plot_hit_frequency_boxplot(df_sig, OUT_DIR)

    print("\nSaved plots:")
    for name in [
        "significant_motifs_by_pipeline_pct.png",
        "top_significant_motifs_log10ER_fixed.png",
        "scatter_log10ER_vs_minuslog10FDR.png",
        "significant_motif_classes_stacked.png",
        "top_motifs_amp_hit_frequency_pct.png",
        "significant_motifs_amp_hit_frequency_boxplot.png",
    ]:
        print("-", OUT_DIR / name)

    if not df_sig.empty:
        print("\nTop significant motifs preview:")
        preview_cols = [
            "Clustering", "Tool", "Motif",
            "Enrichment_ratio", "ER_corrected", "ER_is_infinite", "FDR_corrected",
            "AMP_sequences_with_hit", "nonAMP_sequences_with_hit",
            "Total_hits", "AMP_hit_frequency_pct",
        ]
        preview_cols = [c for c in preview_cols if c in top_global.columns]
        print(top_global[preview_cols].head(10).to_string(index=False))

    inf_n = int(df["ER_is_infinite"].sum()) if "ER_is_infinite" in df.columns else 0
    print(f"\nRows with infinite / background-zero ER: {inf_n}")
    print(f"  → Motif present in AMPs but in 0 of the {NONAMP_BACKGROUND_FACTOR}× non-AMP background.")
    print("  → Use together with hit counts and FDR, not alone.")
    print("\nDone.\n")

if __name__ == "__main__":
    main()