from __future__ import annotations

from pathlib import Path
import sys
import re
import csv
import shutil
import subprocess
from typing import Iterable

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# =====================================================
# CONFIG
# =====================================================
TOMTOM_ENV = "amp_meme"
TOMTOM_QVALUE_THRESHOLD = 0.00001
TOMTOM_MIN_OVERLAP = 6
MIN_COVERAGE = 0.8
REQUIRE_RECIPROCAL_HIT = True
TOMTOM_DIST = "pearson"
TOMTOM_MINIMAL = True

MIN_PIPELINES_ROBUST = 2
SUSPICIOUS_FAMILY_SIZE = 8

CLASS_COLORS = {
    "Cys-rich": "#e07b54",
    "Gly-rich": "#6ab187",
    "Cationic-rich": "#5b8db8",
    "Hydrophobic/signal-like": "#a97fc4",
    "Other": "#aaaaaa",
}


# =====================================================
# REPO ROOT
# =====================================================

def find_repo_root(start: Path) -> Path:
    markers = [".git", "README.md"]
    start = start.resolve()
    for parent in [start] + list(start.parents):
        if any((parent / m).exists() for m in markers):
            return parent
    return start.parents[2]


PROJECT_ROOT = find_repo_root(Path(__file__).parent)

IN_FILE = PROJECT_ROOT / "results" / "06_motif_statistics" / "02_fimo_reporting" / "fimo_significant_only_cleaned.csv"

MOTIF_FILES = {
    "cdhit_meme": PROJECT_ROOT / "results" / "04_motif_discovery" / "meme_cdhit" / "meme.txt",
    "cdhit_streme": PROJECT_ROOT / "results" / "04_motif_discovery" / "streme_cdhit" / "streme.txt",
    "mmseq_meme": PROJECT_ROOT / "results" / "04_motif_discovery" / "meme_mmseq" / "meme.txt",
    "mmseq_streme": PROJECT_ROOT / "results" / "04_motif_discovery" / "streme_mmseq" / "streme.txt",
}

OUT_DIR = PROJECT_ROOT / "results" / "07_motif_families" / "01_motif_family_analysis_tomtom"
OUT_DIR.mkdir(parents=True, exist_ok=True)

COMBINED_MOTIF_FILE = OUT_DIR / "combined_motifs_for_tomtom.meme"
TOMTOM_OUT_DIR = OUT_DIR / "tomtom_all_vs_all"
TOMTOM_TSV = TOMTOM_OUT_DIR / "tomtom.tsv"


# =====================================================
# HELPERS
# =====================================================

def die(msg: str) -> None:
    print(msg)
    sys.exit(1)


def get_binary(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        die(
            f"Required executable '{name}' not found in PATH.\n"
            f"Activate the correct environment and check:\n"
            f"  which {name}\n"
            f"  {name} --version"
        )
    return path


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


def join_unique_sorted(values: Iterable) -> str:
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


def finite_mean(x: pd.Series) -> float:
    x = pd.to_numeric(x, errors="coerce").replace([np.inf, -np.inf], np.nan)
    return x.mean()


def finite_median(x: pd.Series) -> float:
    x = pd.to_numeric(x, errors="coerce").replace([np.inf, -np.inf], np.nan)
    return x.median()


def _class_color(cls: str) -> str:
    return CLASS_COLORS.get(str(cls), "#aaaaaa")


def canonical_pair(a: str, b: str) -> tuple[str, str]:
    return (a, b) if a <= b else (b, a)


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
# MEME PARSING / COMBINING
# =====================================================

def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def split_meme_header_and_blocks(text: str) -> tuple[list[str], list[list[str]]]:
    lines = text.splitlines()
    header = []
    blocks = []

    i = 0
    while i < len(lines):
        if lines[i].startswith("MOTIF "):
            break
        header.append(lines[i])
        i += 1

    current = []
    while i < len(lines):
        line = lines[i]
        if line.startswith("MOTIF "):
            if current:
                blocks.append(current)
            current = [line]
        else:
            if current:
                current.append(line)
        i += 1

    if current:
        blocks.append(current)

    return header, blocks


def parse_motif_name_from_block(block: list[str]) -> str:
    first = block[0]
    parts = first.split()
    if len(parts) < 2:
        raise ValueError(f"Cannot parse MOTIF line: {first}")
    return sanitize_token(parts[1])


def parse_motif_width_from_block(block: list[str]) -> int | None:
    for line in block:
        if "letter-probability matrix" in line and "w=" in line:
            m = re.search(r"\bw=\s*(\d+)", line)
            if m:
                return int(m.group(1))
    return None


def build_combined_meme_file(motif_files: dict[str, Path], output_path: Path) -> pd.DataFrame:
    all_blocks = []
    metadata_rows = []
    final_header = None

    for pipeline, path in motif_files.items():
        if not path.exists():
            die(f"Motif file not found: {path}")

        text = read_text(path)
        header, blocks = split_meme_header_and_blocks(text)

        if final_header is None:
            final_header = header

        for block in blocks:
            original_id = parse_motif_name_from_block(block)
            width = parse_motif_width_from_block(block)
            tomtom_id = f"{pipeline}__{original_id}"

            new_block = block.copy()
            parts = new_block[0].split()
            parts[1] = tomtom_id
            new_block[0] = " ".join(parts)

            all_blocks.append(new_block)
            metadata_rows.append({
                "tomtom_id": tomtom_id,
                "source_pipeline": pipeline,
                "original_motif_id": original_id,
                "motif_width": width,
            })

    if final_header is None:
        die("Could not parse MEME header from motif files.")

    with output_path.open("w", encoding="utf-8") as f:
        for line in final_header:
            f.write(line + "\n")
        if final_header and final_header[-1].strip() != "":
            f.write("\n")

        for block in all_blocks:
            for line in block:
                f.write(line + "\n")
            f.write("\n")

    return pd.DataFrame(metadata_rows)


# =====================================================
# TOMTOM
# =====================================================

def check_tomtom() -> None:
    cmd = ["conda", "run", "-n", TOMTOM_ENV, "tomtom", "--version"]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        die("Tomtom not working in env 'amp_meme'")

    print("TOMTOM OK:")
    print(result.stdout or result.stderr)


def run_tomtom(query_file: Path, target_file: Path, out_dir: Path) -> None:
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "conda", "run", "-n", "amp_meme", "tomtom",
        "-dist", TOMTOM_DIST,
        "-min-overlap", str(TOMTOM_MIN_OVERLAP),
        "-text",
        str(query_file),
        str(target_file)
    ]

    print("\nRunning Tomtom:")
    print(" ".join(cmd))

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(result.stderr)
        die(f"Tomtom failed with return code {result.returncode}.")

    TOMTOM_TSV.write_text(result.stdout, encoding="utf-8")


def load_tomtom_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        die(f"Tomtom TSV not found: {path}")

    df = pd.read_csv(path, sep="\t", comment="#")
    if df.empty:
        return df

    df.columns = [sanitize_token(c) for c in df.columns]
    return df


def detect_tomtom_columns(df: pd.DataFrame) -> dict[str, str]:
    cols = {c.lower(): c for c in df.columns}

    def pick(options: list[str]) -> str | None:
        for opt in options:
            if opt.lower() in cols:
                return cols[opt.lower()]
        return None

    mapping = {
        "query_id": pick(["Query_ID", "Query ID"]),
        "target_id": pick(["Target_ID", "Target ID"]),
        "pvalue": pick(["p-value", "pvalue"]),
        "evalue": pick(["E-value", "evalue"]),
        "qvalue": pick(["q-value", "qvalue"]),
        "overlap": pick(["Overlap", "overlap"]),
        "offset": pick(["Offset", "offset"]),
        "orientation": pick(["Orientation", "orientation"]),
    }

    needed = ["query_id", "target_id", "qvalue"]
    missing = [k for k in needed if mapping[k] is None]
    if missing:
        die(f"Could not detect required Tomtom columns: {missing}. Found columns: {list(df.columns)}")

    return mapping


def build_tomtom_pairs_detailed(
    tomtom_df: pd.DataFrame,
    motif_meta: pd.DataFrame,
    qvalue_threshold: float,
    min_overlap: int,
    min_coverage: float,
) -> pd.DataFrame:
    if tomtom_df.empty:
        return pd.DataFrame(columns=[
            "Query_ID", "Target_ID",
            "q_value", "E_value", "p_value",
            "Overlap", "Offset", "Orientation",
            "Query_width", "Target_width",
            "Coverage_min_width",
            "Query_pipeline", "Target_pipeline",
            "Query_original", "Target_original"
        ])

    colmap = detect_tomtom_columns(tomtom_df)
    meta = motif_meta.set_index("tomtom_id").to_dict(orient="index")

    rows = []

    for _, row in tomtom_df.iterrows():
        qid = sanitize_token(row[colmap["query_id"]])
        tid = sanitize_token(row[colmap["target_id"]])

        if qid == tid:
            continue

        qval = pd.to_numeric(row[colmap["qvalue"]], errors="coerce")
        if pd.isna(qval) or qval > qvalue_threshold:
            continue

        overlap = pd.to_numeric(row[colmap["overlap"]], errors="coerce") if colmap["overlap"] else np.nan
        if pd.notna(overlap) and overlap < min_overlap:
            continue

        qmeta = meta.get(qid, {})
        tmeta = meta.get(tid, {})

        qwidth = pd.to_numeric(qmeta.get("motif_width"), errors="coerce")
        twidth = pd.to_numeric(tmeta.get("motif_width"), errors="coerce")
        min_width = np.nanmin([qwidth, twidth]) if not (pd.isna(qwidth) and pd.isna(twidth)) else np.nan

        if pd.notna(overlap) and pd.notna(min_width) and min_width > 0:
            coverage = overlap / min_width
        else:
            coverage = np.nan

        if pd.notna(coverage) and coverage < min_coverage:
            continue

        rows.append({
            "Query_ID": qid,
            "Target_ID": tid,
            "q_value": qval,
            "E_value": pd.to_numeric(row[colmap["evalue"]], errors="coerce") if colmap["evalue"] else np.nan,
            "p_value": pd.to_numeric(row[colmap["pvalue"]], errors="coerce") if colmap["pvalue"] else np.nan,
            "Overlap": overlap,
            "Offset": pd.to_numeric(row[colmap["offset"]], errors="coerce") if colmap["offset"] else np.nan,
            "Orientation": sanitize_token(row[colmap["orientation"]]) if colmap["orientation"] else "",
            "Query_width": qwidth,
            "Target_width": twidth,
            "Coverage_min_width": coverage,
            "Query_pipeline": qmeta.get("source_pipeline", ""),
            "Target_pipeline": tmeta.get("source_pipeline", ""),
            "Query_original": qmeta.get("original_motif_id", ""),
            "Target_original": tmeta.get("original_motif_id", ""),
        })

    detailed = pd.DataFrame(rows)
    if detailed.empty:
        return pd.DataFrame(columns=[
            "Query_ID", "Target_ID",
            "q_value", "E_value", "p_value",
            "Overlap", "Offset", "Orientation",
            "Query_width", "Target_width",
            "Coverage_min_width",
            "Query_pipeline", "Target_pipeline",
            "Query_original", "Target_original"
        ])

    return detailed.sort_values(["q_value", "E_value", "p_value"], ascending=[True, True, True])


def build_reciprocal_edges(
    directional_df: pd.DataFrame,
    require_reciprocal: bool,
) -> pd.DataFrame:
    if directional_df.empty:
        return pd.DataFrame(columns=[
            "Motif_A", "Motif_B",
            "best_q_value", "best_E_value", "best_p_value",
            "best_overlap", "mean_overlap",
            "best_coverage", "mean_coverage",
            "reciprocal"
        ])

    records = directional_df.to_dict(orient="records")
    dir_map = {(r["Query_ID"], r["Target_ID"]): r for r in records}

    rows = []
    seen = set()

    for (a, b), rab in dir_map.items():
        pair = canonical_pair(a, b)
        if pair in seen:
            continue
        seen.add(pair)

        rba = dir_map.get((b, a))
        reciprocal = rba is not None

        if require_reciprocal and not reciprocal:
            continue

        overlaps = [x["Overlap"] for x in [rab, rba] if x is not None and pd.notna(x["Overlap"])]
        coverages = [x["Coverage_min_width"] for x in [rab, rba] if x is not None and pd.notna(x["Coverage_min_width"])]
        qvals = [x["q_value"] for x in [rab, rba] if x is not None and pd.notna(x["q_value"])]
        evals = [x["E_value"] for x in [rab, rba] if x is not None and pd.notna(x["E_value"])]
        pvals = [x["p_value"] for x in [rab, rba] if x is not None and pd.notna(x["p_value"])]

        rows.append({
            "Motif_A": pair[0],
            "Motif_B": pair[1],
            "best_q_value": np.nanmin(qvals) if qvals else np.nan,
            "best_E_value": np.nanmin(evals) if evals else np.nan,
            "best_p_value": np.nanmin(pvals) if pvals else np.nan,
            "best_overlap": np.nanmax(overlaps) if overlaps else np.nan,
            "mean_overlap": np.nanmean(overlaps) if overlaps else np.nan,
            "best_coverage": np.nanmax(coverages) if coverages else np.nan,
            "mean_coverage": np.nanmean(coverages) if coverages else np.nan,
            "reciprocal": reciprocal,
        })

    edges = pd.DataFrame(rows)
    if edges.empty:
        return pd.DataFrame(columns=[
            "Motif_A", "Motif_B",
            "best_q_value", "best_E_value", "best_p_value",
            "best_overlap", "mean_overlap",
            "best_coverage", "mean_coverage",
            "reciprocal"
        ])

    return edges.sort_values(["best_q_value", "best_E_value", "best_p_value"], ascending=[True, True, True])


# =====================================================
# PLOTS
# =====================================================

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
    plt.savefig(out_dir / "family_er_by_class.png", dpi=300)
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
    plt.savefig(out_dir / "family_pipeline_heatmap.png", dpi=300)
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
    plt.savefig(out_dir / "family_er_vs_hit_frequency.png", dpi=300)
    plt.close()


def plot_bubble_fdr_er_hits(df: pd.DataFrame, out_dir: Path) -> None:
    fdr_col = "FDR_for_reporting" if "FDR_for_reporting" in df.columns else "FDR_corrected"

    if "Total_hits" not in df.columns or fdr_col not in df.columns:
        return

    plot_df = df.dropna(subset=["Enrichment_ratio", fdr_col]).copy()
    if plot_df.empty:
        return

    plot_df["ER_for_plot"] = pd.to_numeric(plot_df["Enrichment_ratio"], errors="coerce").replace([np.inf, -np.inf], np.nan)
    plot_df = plot_df.dropna(subset=["ER_for_plot"])
    if plot_df.empty:
        return

    plot_df["log10_ER"] = np.log10(plot_df["ER_for_plot"].clip(lower=1e-12))
    plot_df["minus_log10_FDR"] = -np.log10(pd.to_numeric(plot_df[fdr_col], errors="coerce").clip(lower=1e-300))

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
    plt.savefig(out_dir / "bubble_er_fdr_hits.png", dpi=300)
    plt.close()


# =====================================================
# MAIN
# =====================================================

def main() -> None:
    print("\nMOTIF FAMILY ANALYSIS — TOMTOM-BASED\n")
    print("Input FIMO summary:", IN_FILE)
    print("Output directory:", OUT_DIR)

    check_tomtom()

    if not IN_FILE.exists():
        die(f"Input file not found: {IN_FILE}")

    missing_motif_files = [str(p) for p in MOTIF_FILES.values() if not p.exists()]
    if missing_motif_files:
        die("Missing motif discovery files:\n" + "\n".join(missing_motif_files))

    # --------------------------------------------------
    # Load FIMO reporting table
    # --------------------------------------------------
    df = pd.read_csv(IN_FILE, sep=";", encoding="utf-8-sig")

    if "FDR_for_reporting" in df.columns and "FDR_corrected" not in df.columns:
        df["FDR_corrected"] = pd.to_numeric(df["FDR_for_reporting"], errors="coerce")

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
    if "FDR_for_reporting" in df.columns:
        numeric_cols.append("FDR_for_reporting")

    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    unique_motifs_in_fimo = sorted(df["Motif"].dropna().unique().tolist())
    print("Unique significant motifs in FIMO table:", len(unique_motifs_in_fimo))

    # --------------------------------------------------
    # Build combined MEME file and motif metadata
    # --------------------------------------------------
    motif_meta = build_combined_meme_file(MOTIF_FILES, COMBINED_MOTIF_FILE)
    motif_meta["clean_original_motif_id"] = motif_meta["original_motif_id"].map(clean_motif)

    print("Combined motif file:", COMBINED_MOTIF_FILE)
    print("Total motifs in combined MEME:", len(motif_meta))

    # --------------------------------------------------
    # Run Tomtom all-vs-all
    # --------------------------------------------------
    run_tomtom(COMBINED_MOTIF_FILE, COMBINED_MOTIF_FILE, TOMTOM_OUT_DIR)
    print("Tomtom TSV:", TOMTOM_TSV)

    tomtom_df = load_tomtom_tsv(TOMTOM_TSV)
    print("Tomtom rows loaded:", len(tomtom_df))

    # --------------------------------------------------
    # Detailed directional pairs
    # --------------------------------------------------
    directional_pairs = build_tomtom_pairs_detailed(
        tomtom_df=tomtom_df,
        motif_meta=motif_meta,
        qvalue_threshold=TOMTOM_QVALUE_THRESHOLD,
        min_overlap=TOMTOM_MIN_OVERLAP,
        min_coverage=MIN_COVERAGE,
    )

    print("Directional Tomtom hits kept after filtering:", len(directional_pairs))

    # --------------------------------------------------
    # Keep only motifs present in FIMO significant table
    # --------------------------------------------------
    valid_clean_motifs = set(unique_motifs_in_fimo)
    motif_map = motif_meta[
        ["tomtom_id", "source_pipeline", "original_motif_id", "clean_original_motif_id", "motif_width"]
    ].copy()
    motif_map = motif_map[motif_map["clean_original_motif_id"].isin(valid_clean_motifs)].copy()
    tomtom_ids_for_fimo = set(motif_map["tomtom_id"])

    if directional_pairs.empty:
        print("Warning: no directional Tomtom hits passed the filters.")
    else:
        directional_pairs = directional_pairs[
            directional_pairs["Query_ID"].isin(tomtom_ids_for_fimo) &
            directional_pairs["Target_ID"].isin(tomtom_ids_for_fimo)
        ].copy()

    # --------------------------------------------------
    # Build reciprocal / undirected edges
    # --------------------------------------------------
    reciprocal_edges = build_reciprocal_edges(
        directional_df=directional_pairs,
        require_reciprocal=REQUIRE_RECIPROCAL_HIT,
    )

    print("Undirected Tomtom edges kept:", len(reciprocal_edges))

    # --------------------------------------------------
    # Convert Tomtom IDs to cleaned motif names for family assignment
    # --------------------------------------------------
    id_to_clean = dict(zip(motif_map["tomtom_id"], motif_map["clean_original_motif_id"]))

    converted_edges = []
    seen_clean_pairs = set()

    for _, row in reciprocal_edges.iterrows():
        a_clean = id_to_clean.get(row["Motif_A"])
        b_clean = id_to_clean.get(row["Motif_B"])

        if not a_clean or not b_clean or a_clean == b_clean:
            continue

        pair = canonical_pair(a_clean, b_clean)
        if pair in seen_clean_pairs:
            continue
        seen_clean_pairs.add(pair)

        converted_edges.append({
            "Motif_A": pair[0],
            "Motif_B": pair[1],
            "best_q_value": row["best_q_value"],
            "best_E_value": row["best_E_value"],
            "best_p_value": row["best_p_value"],
            "best_overlap": row["best_overlap"],
            "mean_overlap": row["mean_overlap"],
            "best_coverage": row["best_coverage"],
            "mean_coverage": row["mean_coverage"],
            "reciprocal": row["reciprocal"],
            "Rule": (
                f"tomtom_q<={TOMTOM_QVALUE_THRESHOLD};"
                f"overlap>={TOMTOM_MIN_OVERLAP};"
                f"coverage>={MIN_COVERAGE};"
                f"reciprocal={REQUIRE_RECIPROCAL_HIT}"
            ),
        })

    clean_edges_df = pd.DataFrame(converted_edges)
    if clean_edges_df.empty:
        clean_edges_df = pd.DataFrame(columns=[
            "Motif_A", "Motif_B",
            "best_q_value", "best_E_value", "best_p_value",
            "best_overlap", "mean_overlap",
            "best_coverage", "mean_coverage",
            "reciprocal", "Rule"
        ])

    # --------------------------------------------------
    # Similarity graph + union-find families
    # --------------------------------------------------
    unique_motifs = sorted(df["Motif"].dropna().unique().tolist())
    uf = UnionFind(unique_motifs)

    for _, row in clean_edges_df.iterrows():
        a = row["Motif_A"]
        b = row["Motif_B"]
        if a in uf.parent and b in uf.parent:
            uf.union(a, b)

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
    # Component size diagnostics
    # --------------------------------------------------
    family_sizes = (
        df.groupby("Family_ID")["Motif"]
        .nunique()
        .reset_index(name="n_unique_motifs")
        .sort_values("n_unique_motifs", ascending=False)
    )

    suspicious_families = (
        family_sizes[family_sizes["n_unique_motifs"] >= SUSPICIOUS_FAMILY_SIZE]
        .copy()
    )
    if not suspicious_families.empty:
        suspicious_families = suspicious_families.merge(
            df.groupby("Family_ID")["Motif"].apply(join_unique_sorted).reset_index(name="Family_members"),
            on="Family_ID",
            how="left",
        )

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
        "tomtom_raw_results": tomtom_df,
        "tomtom_directional_pairs_filtered": directional_pairs,
        "tomtom_edges_filtered_reciprocal": clean_edges_df,
        "motif_family_assignments_tomtom": df,
        "motif_family_summary_tomtom": family_summary,
        "motif_family_summary_robust_tomtom": robust_families,
        "motif_family_pipeline_presence_tomtom": presence,
        "motif_family_summary_for_memory_tomtom": memory_table,
        "suspicious_large_families_tomtom": suspicious_families,
    }

    print("\nSaved:")
    for name, table in outputs.items():
        csv_path = OUT_DIR / f"{name}.csv"
        save_csv(table, csv_path)
        print("-", csv_path)

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

    if not suspicious_families.empty:
        print("\nWarning: suspiciously large families detected:")
        print(suspicious_families.head(10).to_string(index=False))

    print("\nDone.\n")


if __name__ == "__main__":
    main()