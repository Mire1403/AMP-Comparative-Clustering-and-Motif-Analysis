from __future__ import annotations

from pathlib import Path
import sys
import re
import csv

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# =====================================================
# CONFIG
# =====================================================

TOP_N_FAMILIES = 15

plt.rcParams["figure.dpi"] = 150
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["font.size"] = 11

CLASS_COLORS = {
    "Cys-rich": "#e07b54",
    "Gly-rich": "#6ab187",
    "Cationic-rich": "#5b8db8",
    "Hydrophobic/signal-like": "#a97fc4",
    "Other": "#aaaaaa",
}

def _class_color(cls: str) -> str:
    return CLASS_COLORS.get(str(cls).strip(), "#aaaaaa")

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

IN_FILE = PROJECT_ROOT / "results" / "07_motif_families" / "01_motif_family_analysis_tomtom" / "motif_family_summary_for_memory_tomtom.csv"
IN_FULL = PROJECT_ROOT / "results" / "07_motif_families" / "01_motif_family_analysis_tomtom" / "motif_family_summary_tomtom.csv"
IN_PRES = PROJECT_ROOT / "results" / "07_motif_families" / "01_motif_family_analysis_tomtom" / "motif_family_pipeline_presence_tomtom.csv"

OUT_DIR = PROJECT_ROOT / "results" / "07_motif_families" / "02_final_reporting"
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

def sanitize_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].map(lambda x: sanitize_token(x) if pd.notna(x) else x)
    return df

def save_csv(df: pd.DataFrame, path: Path) -> None:
    sanitize_dataframe(df).to_csv(
        path,
        index=False,
        sep=";",
        encoding="utf-8-sig",
        quoting=csv.QUOTE_MINIMAL
    )

def save_xlsx(df: pd.DataFrame, path: Path) -> None:
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        sanitize_dataframe(df).to_excel(writer, index=False, sheet_name="Results")

def save_parquet(df: pd.DataFrame, path: Path) -> None:
    try:
        sanitize_dataframe(df).to_parquet(path, index=False)
    except Exception as e:
        print(f"  Warning: could not write parquet ({e}). Install pyarrow if needed.")

def _legend_handles():
    return [mpatches.Patch(color=v, label=k) for k, v in CLASS_COLORS.items()]

def safe_log10(x: pd.Series, floor: float = 1e-300) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    return -np.log10(x.clip(lower=floor))

def finite_numeric(x: pd.Series) -> pd.Series:
    return pd.to_numeric(x, errors="coerce").replace([np.inf, -np.inf], np.nan)

# =====================================================
# FIGURES
# =====================================================

def fig_class_distribution(df: pd.DataFrame, out_dir: Path) -> Path | None:
    class_dist = (
        df.groupby("Representative_class", dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values("count", ascending=False)
    )

    if class_dist.empty:
        return None

    colors = [_class_color(c) for c in class_dist["Representative_class"]]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(class_dist["Representative_class"], class_dist["count"], color=colors)
    ax.set_ylabel("Number of robust motif families")
    ax.set_xlabel("Motif class")
    ax.set_title("Distribution of robust motif family classes")
    plt.xticks(rotation=30, ha="right")
    ax.legend(handles=_legend_handles(), fontsize=8, loc="upper right")
    plt.tight_layout()

    out = out_dir / "motif_class_distribution.png"
    plt.savefig(out)
    plt.close()
    return out

def fig_top_families_er(df: pd.DataFrame, out_dir: Path, top_n: int = TOP_N_FAMILIES) -> Path | None:
    top = df.dropna(subset=["median_ER"]).copy()
    if top.empty:
        return None

    top = top.sort_values(
        ["n_pipelines", "median_ER", "n_unique_motifs"],
        ascending=[False, False, False]
    ).head(top_n)
    top = top.sort_values("median_ER", ascending=True)

    colors = [_class_color(c) for c in top["Representative_class"]]

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.barh(range(len(top)), top["median_ER"].values, color=colors)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top["Representative_motif"].astype(str).tolist(), fontsize=9)
    ax.set_xlabel("Median Enrichment Ratio")
    ax.set_title(f"Top {top_n} robust motif families across pipelines")
    ax.legend(handles=_legend_handles(), fontsize=8, loc="lower right")
    plt.tight_layout()

    out = out_dir / "robust_motif_families_top.png"
    plt.savefig(out)
    plt.close()
    return out

def fig_pipeline_heatmap(presence_path: Path, out_dir: Path) -> Path | None:
    if not presence_path.exists():
        print(f"  Skipping heatmap — {presence_path} not found.")
        return None

    pres = pd.read_csv(presence_path, sep=";", encoding="utf-8-sig")
    if pres.empty or "Family_ID" not in pres.columns:
        return None

    pres = sanitize_dataframe(pres).set_index("Family_ID")
    if pres.empty:
        return None

    pipeline_cols = [c for c in pres.columns if c != "Family_ID"]
    pres = pres[pipeline_cols]

    fig, ax = plt.subplots(figsize=(max(6, len(pipeline_cols) * 1.5), max(4, len(pres) * 0.30)))
    im = ax.imshow(pres.values.astype(float), aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)

    ax.set_xticks(range(len(pipeline_cols)))
    ax.set_xticklabels(pipeline_cols, rotation=35, ha="right", fontsize=9)
    ax.set_yticks(range(len(pres)))
    ax.set_yticklabels(pres.index.tolist(), fontsize=7)
    ax.set_title("Motif family presence across pipelines")
    plt.colorbar(im, ax=ax, label="Present (1) / Absent (0)")
    plt.tight_layout()

    out = out_dir / "family_pipeline_heatmap.png"
    plt.savefig(out)
    plt.close()
    return out

def fig_er_vs_hit_freq(df: pd.DataFrame, out_dir: Path) -> Path | None:
    freq_col = next((c for c in ["mean_AMP_hit_frequency_pct", "AMP_hit_frequency_pct"] if c in df.columns), None)
    if freq_col is None:
        return None

    plot_df = df.dropna(subset=["median_ER", freq_col]).copy()
    if plot_df.empty:
        return None

    colors = [_class_color(c) for c in plot_df["Representative_class"]]
    sizes = (plot_df["n_unique_motifs"].fillna(1) * 40).clip(lower=30)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.scatter(
        plot_df["median_ER"],
        plot_df[freq_col],
        c=colors,
        s=sizes,
        alpha=0.75,
        edgecolors="white",
        linewidths=0.5
    )

    for _, row in plot_df.nlargest(8, "median_ER").iterrows():
        ax.annotate(
            str(row["Representative_motif"]),
            (row["median_ER"], row[freq_col]),
            fontsize=7,
            ha="left",
            va="bottom",
            xytext=(3, 3),
            textcoords="offset points"
        )

    ax.set_xlabel("Median Enrichment Ratio")
    ax.set_ylabel("Mean AMP hit frequency (%)")
    ax.set_title("Enrichment vs AMP hit frequency\n(bubble size = n unique motifs per family)")
    ax.legend(handles=_legend_handles(), fontsize=8, loc="upper left")
    plt.tight_layout()

    out = out_dir / "family_er_vs_hit_frequency.png"
    plt.savefig(out)
    plt.close()
    return out

def fig_bubble_fdr_er(df: pd.DataFrame, out_dir: Path) -> Path | None:
    if "median_ER" not in df.columns or "min_FDR" not in df.columns:
        return None

    plot_df = df.copy()
    plot_df["median_ER"] = finite_numeric(plot_df["median_ER"])
    plot_df["min_FDR"] = pd.to_numeric(plot_df["min_FDR"], errors="coerce")
    plot_df = plot_df.dropna(subset=["median_ER", "min_FDR"])
    if plot_df.empty:
        return None

    plot_df["log10_ER"] = np.log10(plot_df["median_ER"].clip(lower=1e-12))
    plot_df["minus_log10_FDR"] = -np.log10(plot_df["min_FDR"].clip(lower=1e-300))

    size_col = "mean_Total_hits" if "mean_Total_hits" in plot_df.columns else "n_unique_motifs"
    sizes = (plot_df[size_col].fillna(1).clip(lower=1) * 3).clip(upper=400)
    colors = [_class_color(c) for c in plot_df["Representative_class"]]

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.scatter(
        plot_df["log10_ER"],
        plot_df["minus_log10_FDR"],
        s=sizes,
        c=colors,
        alpha=0.65,
        edgecolors="white",
        linewidths=0.4
    )

    for _, row in plot_df.nlargest(8, "minus_log10_FDR").iterrows():
        ax.annotate(
            str(row["Representative_motif"]),
            (row["log10_ER"], row["minus_log10_FDR"]),
            fontsize=7,
            ha="left",
            va="bottom",
            xytext=(3, 3),
            textcoords="offset points"
        )

    size_label = "Total FIMO hits" if size_col == "mean_Total_hits" else "n unique motifs"
    ax.set_xlabel("log10(Median Enrichment Ratio)")
    ax.set_ylabel("-log10(min FDR)")
    ax.set_title(f"Motif family enrichment landscape\n(bubble size = {size_label})")
    ax.legend(handles=_legend_handles(), fontsize=8)
    plt.tight_layout()

    out = out_dir / "family_bubble_er_fdr.png"
    plt.savefig(out)
    plt.close()
    return out

def fig_n_pipelines_bar(df: pd.DataFrame, out_dir: Path) -> Path | None:
    if "n_pipelines" not in df.columns or df.empty:
        return None

    counts = df["n_pipelines"].value_counts().sort_index()
    if counts.empty:
        return None

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.bar(counts.index.astype(str), counts.values, color="#5b8db8")
    ax.set_xlabel("Number of pipelines")
    ax.set_ylabel("Number of families")
    ax.set_title("Robustness: how many pipelines share each family")
    plt.tight_layout()

    out = out_dir / "family_robustness_npipelines.png"
    plt.savefig(out)
    plt.close()
    return out

# =====================================================
# MAIN
# =====================================================

def main() -> None:
    print("\nFINAL MOTIF REPORTING\n")
    print("Input (memory):", IN_FILE)
    print("Input (full):  ", IN_FULL)
    print("Output:", OUT_DIR)

    if not IN_FILE.exists():
        die(f"Input file not found: {IN_FILE}")

    df = pd.read_csv(IN_FILE, sep=";", encoding="utf-8-sig")

    required = [
        "Family_ID", "Representative_motif", "Representative_class",
        "n_unique_motifs", "n_pipelines", "Pipelines_present",
        "median_ER", "min_FDR", "Family_members",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        die(f"Missing required columns: {missing}")

    df = sanitize_dataframe(df)

    for col in ["Representative_motif", "Representative_class", "Pipelines_present", "Family_members"]:
        df[col] = df[col].astype(str).str.strip()

    numeric_cols = ["n_unique_motifs", "n_pipelines", "median_ER", "min_FDR"]
    for col in ["mean_Total_hits", "mean_AMP_hit_frequency_pct", "n_infinite_ER"]:
        if col in df.columns:
            numeric_cols.append(col)

    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["median_ER"] = finite_numeric(df["median_ER"])

    # Optional full family summary
    df_full = None
    if IN_FULL.exists():
        df_full = pd.read_csv(IN_FULL, sep=";", encoding="utf-8-sig")
        df_full = sanitize_dataframe(df_full)
        for col in ["median_ER", "min_FDR", "n_unique_motifs", "mean_Total_hits"]:
            if col in df_full.columns:
                df_full[col] = pd.to_numeric(df_full[col], errors="coerce")
        if "median_ER" in df_full.columns:
            df_full["median_ER"] = finite_numeric(df_full["median_ER"])

    final_table = df.sort_values(
        ["n_pipelines", "min_FDR", "median_ER", "n_unique_motifs"],
        ascending=[False, True, False, False]
    ).copy()

    class_dist = (
        final_table.groupby("Representative_class", dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values("count", ascending=False)
    )

    summary_stats = pd.DataFrame([{
        "n_robust_families": len(final_table),
        "n_classes": final_table["Representative_class"].nunique(dropna=True),
        "mean_n_pipelines": final_table["n_pipelines"].mean(),
        "median_n_pipelines": final_table["n_pipelines"].median(),
        "mean_median_ER": finite_numeric(final_table["median_ER"]).mean(),
        "median_median_ER": finite_numeric(final_table["median_ER"]).median(),
    }])

    top_shortlist = final_table.head(TOP_N_FAMILIES).copy()

    tables = {
        "motif_class_distribution": class_dist,
        "robust_motif_families_table": final_table,
        "robust_motif_families_top_shortlist": top_shortlist,
        "robust_motif_families_summary_stats": summary_stats,
    }

    print("\nSaved tables:")
    for name, table in tables.items():
        save_csv(table, OUT_DIR / f"{name}.csv")
        save_xlsx(table, OUT_DIR / f"{name}.xlsx")
        save_parquet(table, OUT_DIR / f"{name}.parquet")
        print("-", OUT_DIR / f"{name}.csv")

    print("\nGenerating figures...")
    saved_figs = []

    for fig_path in [
        fig_class_distribution(final_table, OUT_DIR),
        fig_top_families_er(final_table, OUT_DIR),
        fig_n_pipelines_bar(final_table, OUT_DIR),
        fig_pipeline_heatmap(IN_PRES, OUT_DIR),
        fig_er_vs_hit_freq(final_table, OUT_DIR),
        fig_bubble_fdr_er(df_full if df_full is not None else final_table, OUT_DIR),
    ]:
        if fig_path:
            saved_figs.append(fig_path)

    for p in saved_figs:
        print("-", p)

    print("\nMotif class distribution:")
    print(class_dist.to_string(index=False))

    print("\nTop robust motif families:")
    preview_cols = [
        "Family_ID",
        "Representative_motif",
        "Representative_class",
        "n_pipelines",
        "median_ER",
        "Pipelines_present",
    ]
    for extra in ["mean_AMP_hit_frequency_pct", "n_infinite_ER"]:
        if extra in final_table.columns:
            preview_cols.append(extra)

    print(
        final_table[[c for c in preview_cols if c in final_table.columns]]
        .head(10)
        .to_string(index=False)
    )

    print("\nDone.\n")

if __name__ == "__main__":
    main()