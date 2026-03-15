from __future__ import annotations

from pathlib import Path
import sys
import re
import csv
from itertools import combinations
from difflib import SequenceMatcher

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =====================================================
# CONFIG
# =====================================================

MAX_HAMMING = 2
MIN_SIM_RATIO = 0.70
MIN_SHARED_SUBSTRING = 5
MIN_PIPELINES_ROBUST = 2

# More conservative grouping:
USE_SEQUENCE_RATIO_FOR_EQUAL_LENGTH_ONLY = True
ALLOW_SUBSTRING_RULE = True

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

IN_FILE = PROJECT_ROOT / "results" / "statistics" / "04_fimo_reporting" / "fimo_significant_only_cleaned.csv"
OUT_DIR = PROJECT_ROOT / "results" / "statistics" / "05_motif_family_analysis"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# =====================================================
# HELPERS
# =====================================================

def die(msg: str) -> None:
    print(msg)
    sys.exit(1)

def sanitize_token(x: str) -> str:
    x = str(x)
    x = x.replace("\t", " ").replace("\n", " ").replace("\r", " ")
    x = re.sub(r"\s+", " ", x).strip()
    return x

def clean_motif(x: str) -> str:
    x = sanitize_token(x)
    x = re.sub(r"^\d+[-_]", "", x)
    return x

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

def hamming(s1: str, s2: str) -> int | None:
    if len(s1) != len(s2):
        return None
    return sum(a != b for a, b in zip(s1, s2))

def longest_common_substring_len(a: str, b: str) -> int:
    m = [[0] * (1 + len(b)) for _ in range(1 + len(a))]
    longest = 0
    for i in range(1, len(a) + 1):
        for j in range(1, len(b) + 1):
            if a[i - 1] == b[j - 1]:
                m[i][j] = m[i - 1][j - 1] + 1
                longest = max(longest, m[i][j])
    return longest

def are_similar(a: str, b: str) -> tuple[bool, str]:
    """
    Conservative similarity rules:
    1) small Hamming distance for equal lengths
    2) one motif contains the other, if reasonably long
    3) long shared substring
    4) SequenceMatcher ratio, but optionally only for equal lengths
    """
    a = str(a)
    b = str(b)

    d = hamming(a, b)
    if d is not None and d <= MAX_HAMMING:
        return True, f"hamming<={MAX_HAMMING}"

    if ALLOW_SUBSTRING_RULE and (a in b or b in a):
        if min(len(a), len(b)) >= MIN_SHARED_SUBSTRING:
            return True, "substring"

    lcs = longest_common_substring_len(a, b)
    if lcs >= MIN_SHARED_SUBSTRING:
        return True, f"shared_substring>={MIN_SHARED_SUBSTRING}"

    if USE_SEQUENCE_RATIO_FOR_EQUAL_LENGTH_ONLY and len(a) != len(b):
        return False, ""

    sim = SequenceMatcher(None, a, b).ratio()
    if sim >= MIN_SIM_RATIO:
        return True, f"ratio>={MIN_SIM_RATIO:.2f}"

    return False, ""

def join_unique_sorted(values) -> str:
    vals = [sanitize_token(v) for v in values if pd.notna(v) and sanitize_token(v)]
    return ";".join(sorted(set(vals)))

def sanitize_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].map(lambda x: sanitize_token(x) if pd.notna(x) else x)
    return df

def save_csv(df: pd.DataFrame, path: Path) -> None:
    sanitize_dataframe(df).to_csv(
        path, index=False, sep=";", encoding="utf-8-sig", quoting=csv.QUOTE_MINIMAL
    )

def save_xlsx(df: pd.DataFrame, path: Path, sheet_name: str = "Results") -> None:
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        sanitize_dataframe(df).to_excel(writer, index=False, sheet_name=sheet_name)

def finite_mean(x: pd.Series) -> float:
    x = pd.to_numeric(x, errors="coerce").replace([np.inf, -np.inf], np.nan)
    return x.mean()

def finite_median(x: pd.Series) -> float:
    x = pd.to_numeric(x, errors="coerce").replace([np.inf, -np.inf], np.nan)
    return x.median()

# =====================================================
# UNION-FIND
# =====================================================

class UnionFind:
    def __init__(self, items):
        self.parent = {x: x for x in items}

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra

# =====================================================
# PLOTS
# =====================================================

CLASS_COLORS = {
    "Cys-rich": "#e07b54",
    "Gly-rich": "#6ab187",
    "Cationic-rich": "#5b8db8",
    "Hydrophobic/signal-like": "#a97fc4",
    "Other": "#aaaaaa",
}

def _class_color(cls: str) -> str:
    return CLASS_COLORS.get(str(cls), "#aaaaaa")

def plot_family_er_by_class(family_summary: pd.DataFrame, out_dir: Path) -> None:
    df = family_summary.dropna(subset=["median_ER"]).copy()
    if df.empty:
        return

    df = df.sort_values(
        ["n_pipelines", "median_ER", "n_unique_motifs"],
        ascending=[False, False, False]
    ).head(20)
    df = df.sort_values("median_ER", ascending=True)

    colors = [_class_color(c) for c in df["Representative_class"]]

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.barh(range(len(df)), df["median_ER"].values, color=colors)

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["Representative_motif"].astype(str).tolist(), fontsize=9)
    ax.set_xlabel("Median Enrichment Ratio")
    ax.set_title("Top robust motif families — colored by class")

    from matplotlib.patches import Patch
    handles = [Patch(color=v, label=k) for k, v in CLASS_COLORS.items()]
    ax.legend(handles=handles, loc="lower right", fontsize=8)

    plt.tight_layout()
    plt.savefig(out_dir / "family_er_by_class.png")
    plt.close()

def plot_pipeline_heatmap(presence: pd.DataFrame, out_dir: Path) -> None:
    if presence.empty:
        return

    df = presence.set_index("Family_ID") if "Family_ID" in presence.columns else presence
    if df.empty:
        return

    fig, ax = plt.subplots(figsize=(max(6, len(df.columns) * 1.5), max(4, len(df) * 0.35)))
    im = ax.imshow(df.values, aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)

    ax.set_xticks(range(len(df.columns)))
    ax.set_xticklabels(df.columns.tolist(), rotation=35, ha="right", fontsize=9)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df.index.tolist(), fontsize=7)
    ax.set_title("Motif family presence across pipelines")

    plt.colorbar(im, ax=ax, label="Present (1) / Absent (0)")
    plt.tight_layout()
    plt.savefig(out_dir / "family_pipeline_heatmap.png")
    plt.close()

def plot_er_vs_hit_freq(family_summary: pd.DataFrame, out_dir: Path) -> None:
    if "mean_AMP_hit_frequency_pct" not in family_summary.columns:
        return

    df = family_summary.dropna(subset=["median_ER", "mean_AMP_hit_frequency_pct"]).copy()
    if df.empty:
        return

    colors = [_class_color(c) for c in df["Representative_class"]]
    sizes = (df["n_unique_motifs"].fillna(1) * 40).clip(lower=30)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.scatter(
        df["median_ER"],
        df["mean_AMP_hit_frequency_pct"],
        c=colors,
        s=sizes,
        alpha=0.75,
        edgecolors="white",
        linewidths=0.5,
    )

    top = df.nlargest(8, "median_ER")
    for _, row in top.iterrows():
        ax.annotate(
            str(row["Representative_motif"]),
            (row["median_ER"], row["mean_AMP_hit_frequency_pct"]),
            fontsize=7,
            ha="left",
            va="bottom",
            xytext=(3, 3),
            textcoords="offset points",
        )

    ax.set_xlabel("Median Enrichment Ratio")
    ax.set_ylabel("Mean AMP hit frequency (%)")
    ax.set_title("Enrichment vs AMP hit frequency per family\n(size = n unique motifs)")

    from matplotlib.patches import Patch
    handles = [Patch(color=v, label=k) for k, v in CLASS_COLORS.items()]
    ax.legend(handles=handles, fontsize=8, loc="upper left")

    plt.tight_layout()
    plt.savefig(out_dir / "family_er_vs_hit_frequency.png")
    plt.close()

def plot_bubble_fdr_er_hits(df: pd.DataFrame, out_dir: Path) -> None:
    if "Total_hits" not in df.columns or "FDR_corrected" not in df.columns:
        return

    plot_df = df.dropna(subset=["Enrichment_ratio", "FDR_corrected"]).copy()
    if plot_df.empty:
        return

    plot_df["ER_for_plot"] = pd.to_numeric(plot_df["Enrichment_ratio"], errors="coerce").replace([np.inf, -np.inf], np.nan)
    plot_df = plot_df.dropna(subset=["ER_for_plot"])
    if plot_df.empty:
        return

    plot_df["log10_ER"] = np.log10(plot_df["ER_for_plot"].clip(lower=1e-12))
    plot_df["minus_log10_FDR"] = -np.log10(pd.to_numeric(plot_df["FDR_corrected"], errors="coerce").clip(lower=1e-300))

    sizes = (plot_df["Total_hits"].fillna(1).clip(lower=1) * 2).clip(upper=300)
    colors = [_class_color(c) for c in plot_df.get("Motif_class", ["Other"] * len(plot_df))]

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.scatter(
        plot_df["log10_ER"],
        plot_df["minus_log10_FDR"],
        s=sizes,
        c=colors,
        alpha=0.6,
        edgecolors="white",
        linewidths=0.4,
    )

    ax.set_xlabel("log10(Enrichment Ratio)")
    ax.set_ylabel("-log10(FDR)")
    ax.set_title("Motif enrichment landscape\n(bubble size = Total FIMO hits)")

    from matplotlib.patches import Patch
    handles = [Patch(color=v, label=k) for k, v in CLASS_COLORS.items()]
    ax.legend(handles=handles, fontsize=8)

    plt.tight_layout()
    plt.savefig(out_dir / "bubble_er_fdr_hits.png")
    plt.close()

# =====================================================
# MAIN
# =====================================================

def main():
    print("\nMOTIF FAMILY ANALYSIS\n")
    print("Input:", IN_FILE)
    print("Output:", OUT_DIR)

    if not IN_FILE.exists():
        die(f"Input file not found: {IN_FILE}")

    df = pd.read_csv(IN_FILE)

    needed = [
        "Clustering", "Tool", "Motif",
        "Enrichment_ratio", "FDR_corrected",
        "AMP_sequences_with_hit", "nonAMP_sequences_with_hit",
    ]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        die(f"Missing required columns: {missing}")

    has_hits = "Total_hits" in df.columns
    has_freq = "AMP_hit_frequency_pct" in df.columns
    has_inf_flag = "ER_is_infinite" in df.columns

    df = df.copy()
    df["Clustering"] = df["Clustering"].astype(str).map(sanitize_token).str.lower()
    df["Tool"] = df["Tool"].astype(str).map(sanitize_token).str.lower()
    df["Motif"] = df["Motif"].astype(str).map(clean_motif)
    df["Pipeline"] = df["Clustering"] + "_" + df["Tool"]
    df["Motif_class"] = df["Motif"].map(classify_motif)

    numeric_cols = [
        "Enrichment_ratio", "FDR_corrected",
        "AMP_sequences_with_hit", "nonAMP_sequences_with_hit",
    ]
    if has_hits:
        numeric_cols.append("Total_hits")
    if has_freq:
        numeric_cols.append("AMP_hit_frequency_pct")
    if has_inf_flag:
        numeric_cols.append("ER_is_infinite")

    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    unique_motifs = sorted(df["Motif"].dropna().unique().tolist())
    print("Unique significant motifs:", len(unique_motifs))

    # --------------------------------------------------
    # Similarity graph + union-find families
    # --------------------------------------------------
    uf = UnionFind(unique_motifs)
    edge_rows = []

    for a, b in combinations(unique_motifs, 2):
        ok, rule = are_similar(a, b)
        if ok:
            uf.union(a, b)
            edge_rows.append({
                "Motif_A": a,
                "Motif_B": b,
                "Rule": rule,
            })

    edges_df = (
        pd.DataFrame(edge_rows, columns=["Motif_A", "Motif_B", "Rule"])
        if edge_rows
        else pd.DataFrame(columns=["Motif_A", "Motif_B", "Rule"])
    )

    # --------------------------------------------------
    # Assign family IDs
    # --------------------------------------------------
    root_to_family: dict[str, str] = {}
    family_ids = []
    counter = 1

    for motif in unique_motifs:
        root = uf.find(motif)
        if root not in root_to_family:
            root_to_family[root] = f"FAM_{counter:03d}"
            counter += 1
        family_ids.append((motif, root_to_family[root]))

    fam_df = pd.DataFrame(family_ids, columns=["Motif", "Family_ID"])
    df = df.merge(fam_df, on="Motif", how="left")

    # --------------------------------------------------
    # Representative motif per family
    # --------------------------------------------------
    rep_cols = [
        "Family_ID", "Motif", "Motif_class", "Clustering", "Tool",
        "Enrichment_ratio", "FDR_corrected",
        "AMP_sequences_with_hit", "nonAMP_sequences_with_hit",
    ]
    if has_hits:
        rep_cols.append("Total_hits")
    if has_freq:
        rep_cols.append("AMP_hit_frequency_pct")
    if has_inf_flag:
        rep_cols.append("ER_is_infinite")

    sort_cols = ["Family_ID", "FDR_corrected", "Enrichment_ratio", "AMP_sequences_with_hit"]

    rep_df = (
        df.sort_values(sort_cols, ascending=[True, True, False, False])
        .groupby("Family_ID", as_index=False)
        .first()[rep_cols]
        .rename(columns={
            "Motif": "Representative_motif",
            "Motif_class": "Representative_class",
            "Clustering": "Representative_clustering",
            "Tool": "Representative_tool",
            "Enrichment_ratio": "Representative_ER",
            "FDR_corrected": "Representative_FDR",
            "AMP_sequences_with_hit": "Representative_AMP_hits",
            "nonAMP_sequences_with_hit": "Representative_nonAMP_hits",
            "Total_hits": "Representative_Total_hits",
            "AMP_hit_frequency_pct": "Representative_AMP_hit_freq_pct",
            "ER_is_infinite": "Representative_ER_is_infinite",
        })
    )

    # --------------------------------------------------
    # Family-level aggregations
    # --------------------------------------------------
    family_members = (
        df.groupby("Family_ID")["Motif"]
        .apply(join_unique_sorted)
        .reset_index(name="Family_members")
    )

    family_classes = (
        df.groupby("Family_ID")["Motif_class"]
        .apply(join_unique_sorted)
        .reset_index(name="Classes_present")
    )

    family_pipelines = (
        df.groupby("Family_ID")["Pipeline"]
        .apply(join_unique_sorted)
        .reset_index(name="Pipelines_present")
    )

    agg_dict = dict(
        n_rows=("Motif", "size"),
        n_unique_motifs=("Motif", "nunique"),
        n_pipelines=("Pipeline", "nunique"),
        mean_ER=("Enrichment_ratio", finite_mean),
        median_ER=("Enrichment_ratio", finite_median),
        min_FDR=("FDR_corrected", "min"),
        max_AMP_hits=("AMP_sequences_with_hit", "max"),
        min_nonAMP_hits=("nonAMP_sequences_with_hit", "min"),
    )

    if has_inf_flag:
        agg_dict["n_infinite_ER"] = ("ER_is_infinite", "sum")

    if has_hits:
        agg_dict["mean_Total_hits"] = ("Total_hits", "mean")
        agg_dict["max_Total_hits"] = ("Total_hits", "max")

    if has_freq:
        agg_dict["mean_AMP_hit_frequency_pct"] = ("AMP_hit_frequency_pct", "mean")
        agg_dict["max_AMP_hit_frequency_pct"] = ("AMP_hit_frequency_pct", "max")
        agg_dict["median_AMP_hit_frequency_pct"] = ("AMP_hit_frequency_pct", "median")

    family_stats = df.groupby("Family_ID").agg(**agg_dict).reset_index()

    family_summary = (
        family_stats
        .merge(rep_df, on="Family_ID", how="left")
        .merge(family_members, on="Family_ID", how="left")
        .merge(family_classes, on="Family_ID", how="left")
        .merge(family_pipelines, on="Family_ID", how="left")
        .sort_values(
            ["n_pipelines", "n_unique_motifs", "min_FDR", "median_ER"],
            ascending=[False, False, True, False]
        )
    )

    robust_families = (
        family_summary[family_summary["n_pipelines"] >= MIN_PIPELINES_ROBUST]
        .copy()
        .sort_values(["n_pipelines", "min_FDR", "median_ER"], ascending=[False, True, False])
    )

    presence = (
        df.assign(present=1)
        .pivot_table(
            index="Family_ID",
            columns="Pipeline",
            values="present",
            aggfunc="max",
            fill_value=0,
        )
        .reset_index()
    )

    memory_cols = [
        "Family_ID",
        "Representative_motif",
        "Representative_class",
        "n_unique_motifs",
        "n_pipelines",
        "Pipelines_present",
        "median_ER",
        "min_FDR",
        "Family_members",
    ]
    if has_hits:
        memory_cols.append("mean_Total_hits")
    if has_freq:
        memory_cols.append("mean_AMP_hit_frequency_pct")
    if has_inf_flag:
        memory_cols.append("n_infinite_ER")

    memory_table = robust_families[[c for c in memory_cols if c in robust_families.columns]].copy()

    # --------------------------------------------------
    # Save outputs
    # --------------------------------------------------
    outputs = {
        "motif_family_edges": edges_df,
        "motif_family_assignments": df,
        "motif_family_summary": family_summary,
        "motif_family_summary_robust": robust_families,
        "motif_family_pipeline_presence": presence,
        "motif_family_summary_for_memory": memory_table,
    }

    print("\nSaved:")
    for name, table in outputs.items():
        csv_path = OUT_DIR / f"{name}.csv"
        xlsx_path = OUT_DIR / f"{name}.xlsx"
        save_csv(table, csv_path)
        save_xlsx(table, xlsx_path)
        print("-", csv_path)
        print("-", xlsx_path)

    # --------------------------------------------------
    # Plots
    # --------------------------------------------------
    print("\nGenerating plots...")
    plot_family_er_by_class(robust_families, OUT_DIR)
    plot_pipeline_heatmap(presence, OUT_DIR)
    plot_er_vs_hit_freq(family_summary, OUT_DIR)
    plot_bubble_fdr_er_hits(df, OUT_DIR)

    print("- family_er_by_class.png")
    print("- family_pipeline_heatmap.png")
    if has_freq:
        print("- family_er_vs_hit_frequency.png")
    if has_hits:
        print("- bubble_er_fdr_hits.png")

    # --------------------------------------------------
    # Console summary
    # --------------------------------------------------
    print("\nTop robust families:")
    if not robust_families.empty:
        preview = [
            "Family_ID",
            "Representative_motif",
            "Representative_class",
            "n_unique_motifs",
            "n_pipelines",
            "Pipelines_present",
            "min_FDR",
        ]
        if has_freq:
            preview.append("mean_AMP_hit_frequency_pct")
        if has_inf_flag:
            preview.append("n_infinite_ER")

        print(robust_families[[c for c in preview if c in robust_families.columns]].head(15).to_string(index=False))
    else:
        print("No robust families found with current threshold.")

    print("\nDone.\n")

if __name__ == "__main__":
    main()