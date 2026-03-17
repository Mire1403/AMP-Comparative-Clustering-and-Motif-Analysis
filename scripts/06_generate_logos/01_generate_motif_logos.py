from __future__ import annotations

from pathlib import Path
import sys
import re
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

try:
    import logomaker
except ImportError:
    print("ERROR: logomaker is not installed. Run: pip install logomaker")
    sys.exit(1)

# =====================================================
# CONFIG
# =====================================================

DPI              = 300
MIN_PROB         = 0.05   # minimum probability to show an alternative letter
                          # = 1/20 uniform background (IC > 0 criterion)
TOP_K            = 3      # max alternative letters per position in consensus line
LOGO_HEIGHT      = 1.8    # inches — height of the logomaker panel
CONSENSUS_HEIGHT = 0.7    # inches — height of the consensus text panel
LOGO_COLOR_SCHEME = "chemistry"

# Amino acid colors for consensus line (matching chemistry scheme)
AA_COLORS = {
    "A": "#5b8db8", "V": "#5b8db8", "I": "#5b8db8", "L": "#5b8db8",
    "M": "#5b8db8", "F": "#5b8db8", "W": "#5b8db8", "P": "#5b8db8",
    "S": "#6ab187", "T": "#6ab187", "C": "#6ab187", "Y": "#6ab187",
    "N": "#6ab187", "Q": "#6ab187", "G": "#6ab187",
    "K": "#e07b54", "R": "#e07b54", "H": "#e07b54",
    "D": "#a97fc4", "E": "#a97fc4",
    "B": "#aaaaaa", "Z": "#aaaaaa", "X": "#aaaaaa", "J": "#aaaaaa",
}
LEGEND_GROUPS = [
    ("Hydrophobic",       "#5b8db8"),
    ("Polar uncharged",   "#6ab187"),
    ("Positive (K/R/H)",  "#e07b54"),
    ("Negative (D/E)",    "#a97fc4"),
    ("Other / ambiguous", "#aaaaaa"),
]

def aa_color(aa: str) -> str:
    return AA_COLORS.get(aa.upper(), "#aaaaaa")

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

IN_SIG    = PROJECT_ROOT / "results" / "statistics" / "04_fimo_reporting" / "fimo_significant_only_cleaned.csv"
IN_FAMILY = PROJECT_ROOT / "results" / "statistics" / "05b_motif_family_analysis_tomtom" / "motif_family_assignments_tomtom.csv"

OUT_DIR        = PROJECT_ROOT / "results" / "statistics" / "07b_motif_logos"
OUT_DIR_ALL    = OUT_DIR / "all"
OUT_DIR_ROBUST = OUT_DIR / "robust"
for _d in [OUT_DIR, OUT_DIR_ALL, OUT_DIR_ROBUST]:
    _d.mkdir(parents=True, exist_ok=True)

MIN_PIPELINES_ROBUST = 2

MOTIF_FILES = {
    ("cdhit", "meme"):   PROJECT_ROOT / "results" / "04_motif_discovery" / "meme_cdhit"   / "meme.txt",
    ("mmseq", "meme"):   PROJECT_ROOT / "results" / "04_motif_discovery" / "meme_mmseq"   / "meme.txt",
    ("cdhit", "streme"): PROJECT_ROOT / "results" / "04_motif_discovery" / "streme_cdhit" / "streme.txt",
    ("mmseq", "streme"): PROJECT_ROOT / "results" / "04_motif_discovery" / "streme_mmseq" / "streme.txt",
}

# =====================================================
# HELPERS
# =====================================================

def die(msg: str) -> None:
    print(f"ERROR: {msg}")
    sys.exit(1)

def clean_motif_name(x: str) -> str:
    """Remove numeric STREME prefixes (e.g. '25-LKKAVK' -> 'LKKAVK')."""
    x = str(x).strip()
    x = re.sub(r"^\d+[-_]", "", x)
    return x.strip()

def safe_filename(x: str) -> str:
    x = str(x).strip()
    x = re.sub(r"[^\w\-\.]+", "_", x)
    return x[:100]

def infer_alphabet(lines: List[str]) -> str:
    for ln in lines[:300]:
        if ln.startswith("ALPHABET="):
            return ln.split("=", 1)[1].strip().replace(" ", "")
    return "ACDEFGHIKLMNPQRSTVWY"

# =====================================================
# PARSE MEME / STREME MOTIF FILES
# =====================================================

def parse_meme_motifs(txt_path: Path) -> Dict[str, pd.DataFrame]:
    """
    Parse a MEME or STREME output file.
    Returns: cleaned_motif_name -> letter-probability matrix (DataFrame).
    """
    if not txt_path.exists():
        print(f"  Warning: motif file not found: {txt_path}")
        return {}

    lines    = txt_path.read_text(errors="ignore").splitlines(True)
    alphabet = infer_alphabet(lines)
    result: Dict[str, pd.DataFrame] = {}

    i = 0
    while i < len(lines):
        if lines[i].startswith("MOTIF"):
            parts = lines[i].strip().split()
            if len(parts) < 2:
                i += 1
                continue
            motif_name = clean_motif_name(parts[1])
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
                row_vals = lines[i].strip().split()
                try:
                    probs = list(map(float, row_vals[:len(alphabet)]))
                except Exception:
                    probs = []
                if len(probs) == len(alphabet):
                    matrix.append(probs)
                i += 1

            if len(matrix) != width:
                continue

            result[motif_name] = pd.DataFrame(matrix, columns=list(alphabet))
        else:
            i += 1

    return result

# =====================================================
# CONSENSUS LINE BUILDER
# =====================================================

def build_consensus_positions(df: pd.DataFrame) -> List[List[Tuple[str, float]]]:
    """
    For each position return [(aa, prob), ...] above MIN_PROB,
    limited to TOP_K, sorted by probability descending.
    Always returns at least one letter per position.
    """
    positions = []
    for pos in range(len(df)):
        row     = df.iloc[pos]
        entries = [(aa, float(p)) for aa, p in row.items() if float(p) >= MIN_PROB]
        entries.sort(key=lambda x: x[1], reverse=True)
        if not entries:
            entries = [(row.idxmax(), float(row.max()))]
        positions.append(entries[:TOP_K])
    return positions


def consensus_string(positions: List[List[Tuple[str, float]]]) -> str:
    """Build a readable string, e.g. 'G-L/K/R-G-L/K-L-L-G/P'."""
    return "-".join("/".join(aa for aa, _ in pos) for pos in positions)

# =====================================================
# HYBRID FIGURE: LOGOMAKER + CONSENSUS LINE
# =====================================================

def prob_to_ic(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert a letter-probability matrix to information content (bits).
    IC_i = sum_aa [ p(aa,i) * log2(p(aa,i) / bg(aa)) ]
    Background: uniform over the 20 standard amino acids (0.05 each).
    Values are clipped to >= 0 to avoid floating-point negatives.
    """
    n_aa  = len(df.columns)
    bg    = 1.0 / n_aa
    ic_df = df.copy().astype(float)
    for col in ic_df.columns:
        p = ic_df[col].clip(lower=1e-10)
        ic_df[col] = (p * np.log2(p / bg)).clip(lower=0)
    return ic_df


def draw_hybrid_logo(
    df               : pd.DataFrame,
    motif_name       : str,
    clustering       : str,
    tool             : str,
    family_id        : str,
    er               : Optional[float],
    fdr              : Optional[float],
    amp_pct          : Optional[float],
    is_representative: bool,
    out_png          : Path,
) -> None:
    """
    Two-panel figure — identical layout for all motifs:
      Top panel    : logomaker information-content logo (bits, publication standard)
      Bottom panel : coloured consensus line  G - L/K/R - G - L/K - L - L - G/P
      Bottom        : colour legend (outside panels, never overlapping)

    The only difference between representative and non-representative motifs
    is the '★ REPRESENTATIVE' tag in the title and the gold border.
    """
    n_pos   = len(df)
    fig_w   = max(4.0, n_pos * 0.45 + 1.5)
    # Extra bottom space for legend row
    legend_h = 0.45
    total_h  = LOGO_HEIGHT + CONSENSUS_HEIGHT + 0.55 + legend_h

    fig = plt.figure(figsize=(fig_w, total_h))

    # Three rows: logo | consensus | legend
    gs = fig.add_gridspec(
        3, 1,
        height_ratios = [LOGO_HEIGHT, CONSENSUS_HEIGHT, legend_h],
        hspace        = 0.05,
        top           = 0.88,   # leave room for suptitle
        bottom        = 0.02,
    )
    ax_logo   = fig.add_subplot(gs[0])
    ax_cons   = fig.add_subplot(gs[1])
    ax_legend = fig.add_subplot(gs[2])
    ax_legend.axis("off")

    # ── Title — same format for all, ★ REPRESENTATIVE only for rep ────
    rep_tag    = "  \u2605 REPRESENTATIVE" if is_representative else ""
    er_str     = "ER = \u221e" if (er is None or (isinstance(er, float) and np.isinf(er))) \
                 else f"ER = {er:.1f}"
    fdr_str    = f"FDR = {fdr:.2e}" if fdr is not None else ""
    pct_str    = f"AMP freq. = {amp_pct:.1f}%" if amp_pct is not None else ""
    stats_line = "  |  ".join(s for s in [er_str, fdr_str, pct_str] if s)

    # Title colour: gold only for representative
    title_color = "#b8860b" if is_representative else "#222222"

    fig.suptitle(
        f"{family_id}  \u00b7  {motif_name}  [{clustering}-{tool}]{rep_tag}\n{stats_line}",
        fontsize   = 8.5,
        color      = title_color,
        fontweight = "bold" if is_representative else "normal",
        x=0.02, ha="left", y=0.99, va="top",
    )

    # Gold border only for representative
    if is_representative:
        fig.add_artist(plt.Rectangle(
            (0, 0), 1, 1,
            transform=fig.transFigure,
            fill=False, edgecolor="#b8860b", linewidth=2.0, clip_on=False,
        ))

    # ── Top panel: information-content logo ───────────────────────────
    ic_df = prob_to_ic(df)

    logo = logomaker.Logo(
        ic_df,
        ax           = ax_logo,
        color_scheme = LOGO_COLOR_SCHEME,
        stack_order  = "big_on_top",
        show_spines  = False,
    )
    logo.style_spines(visible=False)
    logo.style_spines(spines=["left"], visible=True)
    ax_logo.set_ylabel("bits", fontsize=6.5, labelpad=2)
    ax_logo.tick_params(axis="y", labelsize=6)
    ax_logo.set_xticks([])

    # ── Middle panel: consensus line ───────────────────────────────────
    ax_cons.set_xlim(-0.5, n_pos - 0.5)
    ax_cons.set_ylim(0, 1)
    ax_cons.axis("off")

    positions = build_consensus_positions(df)

    for pos, letters in enumerate(positions):
        token = "/".join(aa for aa, _ in letters)
        color = aa_color(letters[0][0])
        fs    = 8.5 if len(token) <= 2 else 7.0
        ax_cons.text(pos, 0.62, token,
                     ha="center", va="center",
                     fontsize=fs,
                     fontweight="bold" if len(letters) == 1 else "normal",
                     color=color)
        ax_cons.text(pos, 0.10, str(pos + 1),
                     ha="center", va="bottom", fontsize=5.5, color="#cccccc")

    # Dashes between positions
    for pos in range(n_pos - 1):
        ax_cons.plot([pos + 0.38, pos + 0.62], [0.62, 0.62],
                     color="#e0e0e0", lw=0.7, zorder=0)

    # Italic consensus string below the tokens
    cons_str = consensus_string(positions)
    ax_cons.text(0.5, 0.02, cons_str,
                 ha="center", va="bottom", transform=ax_cons.transAxes,
                 fontsize=5.5, color="#aaaaaa", style="italic")

    # ── Bottom row: colour legend (never overlaps panels) ─────────────
    handles = [mpatches.Patch(color=c, label=lbl) for lbl, c in LEGEND_GROUPS]
    ax_legend.legend(
        handles   = handles,
        loc       = "center",
        fontsize  = 6,
        frameon   = False,
        ncol      = len(LEGEND_GROUPS),
        handlelength = 0.8,
        handletextpad= 0.4,
        columnspacing= 1.0,
    )

    plt.savefig(out_png, dpi=DPI, bbox_inches="tight")
    plt.close()

# =====================================================
# MAIN
# =====================================================

def main() -> None:
    print("\n" + "=" * 60)
    print("  MOTIF HYBRID LOGOS — logomaker + consensus line")
    print("=" * 60 + "\n")
    print(f"Significant motifs : {IN_SIG}")
    print(f"Family assignments : {IN_FAMILY}")
    print(f"Output directory   : {OUT_DIR}\n")

    if not IN_SIG.exists():
        die(f"Significant motifs file not found: {IN_SIG}")
    if not IN_FAMILY.exists():
        die(f"Family assignments file not found: {IN_FAMILY}")

    sig = pd.read_csv(IN_SIG)
    fam = pd.read_csv(IN_FAMILY, sep=";", encoding="utf-8-sig")

    for df_ref in [sig, fam]:
        for col in ["Clustering", "Tool"]:
            if col in df_ref.columns:
                df_ref[col] = df_ref[col].astype(str).str.strip().str.lower()
        if "Motif" in df_ref.columns:
            df_ref["Motif"] = df_ref["Motif"].astype(str).map(clean_motif_name)

    # ── Representatives and robust families ───────────────────────────
    rep_set: set         = set()
    robust_families: set = set()
    summary_path = IN_FAMILY.parent / "motif_family_summary_tomtom.csv"
    if summary_path.exists():
        summary = pd.read_csv(summary_path, sep=";", encoding="utf-8-sig")
        for _, row in summary.iterrows():
            rep_cl    = str(row.get("Representative_clustering", "")).strip().lower()
            rep_tool  = str(row.get("Representative_tool",       "")).strip().lower()
            rep_motif = clean_motif_name(row.get("Representative_motif", ""))
            rep_set.add((rep_cl, rep_tool, rep_motif))
            try:
                if int(row.get("n_pipelines", 1)) >= MIN_PIPELINES_ROBUST:
                    robust_families.add(str(row.get("Family_ID", "")).strip())
            except (ValueError, TypeError):
                pass
        print(f"Representative motifs loaded : {len(rep_set)}")
        print(f"Robust families (>={MIN_PIPELINES_ROBUST} pipelines): {len(robust_families)}")
    else:
        print("Warning: motif_family_summary_tomtom.csv not found.")

    # ── Merge family IDs ──────────────────────────────────────────────
    merge_cols = [c for c in ["Clustering", "Tool", "Motif", "Family_ID"] if c in fam.columns]
    if "Family_ID" in fam.columns:
        sig = sig.merge(fam[merge_cols].drop_duplicates(),
                        on=["Clustering", "Tool", "Motif"], how="left")
    else:
        sig["Family_ID"] = "FAM_unknown"
    sig["Family_ID"] = sig["Family_ID"].fillna("FAM_unknown").astype(str)

    # ── Parse motif files ─────────────────────────────────────────────
    parsed: Dict[Tuple[str, str], Dict[str, pd.DataFrame]] = {}
    for key, path in MOTIF_FILES.items():
        print(f"Parsing {path.name} ...")
        parsed[key] = parse_meme_motifs(path)
        print(f"  -> {len(parsed[key])} motifs loaded")

    er_col  = next((c for c in ["Enrichment_ratio", "ER_corrected"] if c in sig.columns), None)
    fdr_col = next((c for c in ["FDR_corrected", "FDR"]            if c in sig.columns), None)
    pct_col = "AMP_hit_frequency_pct" if "AMP_hit_frequency_pct" in sig.columns else None

    result_rows = []
    total       = len(sig)
    generated   = 0

    print(f"\nGenerating hybrid logos for {total} significant motifs...\n")

    for _, row in sig.iterrows():
        motif      = clean_motif_name(row["Motif"])
        clustering = str(row["Clustering"]).strip().lower()
        tool       = str(row["Tool"]).strip().lower()
        family_id  = str(row["Family_ID"]).strip()
        source_key = (clustering, tool)
        is_robust  = family_id in robust_families
        is_rep     = (clustering, tool, motif) in rep_set

        er      = float(row[er_col])  if er_col  and pd.notna(row.get(er_col))  else None
        fdr     = float(row[fdr_col]) if fdr_col and pd.notna(row.get(fdr_col)) else None
        amp_pct = float(row[pct_col]) if pct_col and pd.notna(row.get(pct_col)) else None

        prefix   = "00_REP__" if is_rep else "01_____"
        filename = f"{prefix}{safe_filename(clustering)}_{safe_filename(tool)}__{safe_filename(motif)}.png"

        fam_dir_all = OUT_DIR_ALL / safe_filename(family_id)
        fam_dir_all.mkdir(parents=True, exist_ok=True)
        out_png_all = fam_dir_all / filename

        out_png_robust = None
        if is_robust:
            fam_dir_rob = OUT_DIR_ROBUST / safe_filename(family_id)
            fam_dir_rob.mkdir(parents=True, exist_ok=True)
            out_png_robust = fam_dir_rob / filename

        motif_df = parsed.get(source_key, {}).get(motif)
        found    = motif_df is not None

        if found:
            draw_hybrid_logo(
                df=motif_df, motif_name=motif,
                clustering=clustering, tool=tool,
                family_id=family_id, er=er, fdr=fdr, amp_pct=amp_pct,
                is_representative=is_rep, out_png=out_png_all,
            )
            if out_png_robust is not None:
                draw_hybrid_logo(
                    df=motif_df, motif_name=motif,
                    clustering=clustering, tool=tool,
                    family_id=family_id, er=er, fdr=fdr, amp_pct=amp_pct,
                    is_representative=is_rep, out_png=out_png_robust,
                )
            generated += 1

        cons_str = consensus_string(build_consensus_positions(motif_df)) if found else ""

        result_rows.append({
            "Family_ID":           family_id,
            "Motif":               motif,
            "Clustering":          clustering,
            "Tool":                tool,
            "Is_representative":   is_rep,
            "Is_robust_family":    is_robust,
            "Consensus_string":    cons_str,
            "Diagram_found":       found,
            "Diagram_file_all":    str(out_png_all.relative_to(OUT_DIR)) if found else "",
            "Diagram_file_robust": str(out_png_robust.relative_to(OUT_DIR))
                                   if (found and out_png_robust) else "",
        })

        rob_tag = " [ROBUST]" if is_robust else ""
        rep_tag = " [\u2605 REP]"  if is_rep    else ""
        status  = "OK  " if found else "MISS"
        print(f"  {status}  {family_id:<10}  {clustering}-{tool:<12}  {motif}{rep_tag}{rob_tag}")

    # ── Summary table ─────────────────────────────────────────────────
    out_table = pd.DataFrame(result_rows)
    out_csv   = OUT_DIR / "logo_summary.csv"
    out_xlsx  = OUT_DIR / "logo_summary.xlsx"
    out_table.to_csv(out_csv,   index=False, sep=";", encoding="utf-8-sig")
    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        out_table.to_excel(writer, index=False, sheet_name="logos")

    robust_n = out_table[out_table["Is_robust_family"] & out_table["Diagram_found"]].shape[0]

    print("\n" + "=" * 60)
    print(f"  Generated  : {generated} / {total} logos")
    print(f"  all/       : {generated} (all significant motifs)")
    print(f"  robust/    : {robust_n} (robust families, >={MIN_PIPELINES_ROBUST} pipelines)")
    print(f"  Layout     : logomaker probability logo + consensus line")
    print(f"  Output     : {OUT_DIR}")
    print(f"  Summary    : {out_csv}")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    main()