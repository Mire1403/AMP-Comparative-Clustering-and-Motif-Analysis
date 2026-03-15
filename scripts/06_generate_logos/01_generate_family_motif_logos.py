from __future__ import annotations

from pathlib import Path
import sys
import re
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import logomaker
except ImportError:
    print("ERROR: logomaker is not installed. Run: pip install logomaker")
    sys.exit(1)

# =====================================================
# CONFIG
# =====================================================

FIG_WIDTH = 8
FIG_HEIGHT = 2.4
DPI = 300
COLOR_SCHEME = "chemistry"   # good for proteins
TOP_N = None                 # None = all families

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

IN_FAMILY = PROJECT_ROOT / "results" / "statistics" / "05_motif_family_analysis" / "motif_family_summary_for_memory.csv"
OUT_DIR = PROJECT_ROOT / "results" / "statistics" / "06_final_reporting" / "motif_logos"
OUT_DIR.mkdir(parents=True, exist_ok=True)

MOTIF_FILES = {
    ("cdhit", "meme"):   PROJECT_ROOT / "results" / "04_motif_discovery" / "meme_cdhit" / "meme.txt",
    ("mmseq", "meme"):   PROJECT_ROOT / "results" / "04_motif_discovery" / "meme_mmseq" / "meme.txt",
    ("cdhit", "streme"): PROJECT_ROOT / "results" / "04_motif_discovery" / "streme_cdhit" / "streme.txt",
    ("mmseq", "streme"): PROJECT_ROOT / "results" / "04_motif_discovery" / "streme_mmseq" / "streme.txt",
}

# =====================================================
# HELPERS
# =====================================================

def die(msg: str) -> None:
    print(msg)
    sys.exit(1)

def clean_motif_name(x: str) -> str:
    x = str(x).strip()
    x = re.sub(r"^\d+[-_]", "", x)
    return x.strip()

def safe_filename(x: str) -> str:
    x = str(x).strip()
    x = re.sub(r"[^\w\-\.]+", "_", x)
    return x[:150]

def infer_alphabet(lines: List[str]) -> str:
    for ln in lines[:300]:
        if ln.startswith("ALPHABET="):
            return ln.split("=", 1)[1].strip().replace(" ", "")
    # protein fallback
    return "ACDEFGHIKLMNPQRSTVWY"

# =====================================================
# PARSE MEME/STREME TXT
# =====================================================

def parse_meme_motifs(txt_path: Path) -> Dict[str, pd.DataFrame]:
    """
    Returns dict:
      cleaned_motif_name -> probability matrix DataFrame
    """
    if not txt_path.exists():
        print(f"Warning: motif file not found: {txt_path}")
        return {}

    lines = txt_path.read_text(errors="ignore").splitlines(True)
    alphabet = infer_alphabet(lines)
    motif_to_df: Dict[str, pd.DataFrame] = {}

    i = 0
    while i < len(lines):
        if lines[i].startswith("MOTIF"):
            parts = lines[i].strip().split()
            if len(parts) < 2:
                i += 1
                continue

            motif_name_raw = parts[1]
            motif_name = clean_motif_name(motif_name_raw)
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

            df = pd.DataFrame(matrix, columns=list(alphabet))
            motif_to_df[motif_name] = df
        else:
            i += 1

    return motif_to_df

# =====================================================
# LOGO DRAWING
# =====================================================

def draw_logo(df: pd.DataFrame, title: str, out_png: Path) -> None:
    fig, ax = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))

    logo = logomaker.Logo(
        df,
        ax=ax,
        color_scheme=COLOR_SCHEME,
        stack_order="big_on_top"
    )
    logo.style_spines(visible=False)
    logo.style_spines(spines=["left", "bottom"], visible=True)

    ax.set_xlabel("Position")
    ax.set_ylabel("Probability")
    ax.set_title(title)

    plt.tight_layout()
    plt.savefig(out_png, dpi=DPI)
    plt.close()

# =====================================================
# MAIN
# =====================================================

def main() -> None:
    print("\nGENERATE MOTIF LOGOS FROM MEME/STREME\n")
    print("Family summary:", IN_FAMILY)
    print("Output dir:", OUT_DIR)

    if not IN_FAMILY.exists():
        die(f"Input family summary not found: {IN_FAMILY}")

    fam = pd.read_csv(IN_FAMILY, sep=";", encoding="utf-8-sig")

    required = ["Family_ID", "Representative_motif"]
    missing = [c for c in required if c not in fam.columns]
    if missing:
        die(f"Missing required columns in family summary: {missing}")

    if TOP_N is not None:
        fam = fam.head(TOP_N).copy()

    # Parse all motif files once
    parsed_sources: Dict[Tuple[str, str], Dict[str, pd.DataFrame]] = {}
    for key, path in MOTIF_FILES.items():
        print(f"Parsing {path} ...")
        parsed_sources[key] = parse_meme_motifs(path)
        print(f"  Loaded {len(parsed_sources[key])} motifs")

    result_rows = []

    for _, row in fam.iterrows():
        family_id = str(row["Family_ID"]).strip()
        rep_motif = clean_motif_name(row["Representative_motif"])

        found = False
        found_source = ""
        out_png = OUT_DIR / f"{family_id}__{safe_filename(rep_motif)}.png"

        # Prefer the pipeline listed in Pipelines_present if available,
        # but fall back to searching all files
        candidate_order = list(parsed_sources.keys())

        if "Pipelines_present" in fam.columns and pd.notna(row.get("Pipelines_present", np.nan)):
            pp = str(row["Pipelines_present"])
            preferred = []
            for token in [t.strip().lower() for t in pp.split(";") if t.strip()]:
                # token format expected: cdhit_meme, mmseq_streme...
                parts = token.split("_")
                if len(parts) == 2:
                    tup = (parts[0], parts[1])
                    if tup in parsed_sources:
                        preferred.append(tup)
            candidate_order = preferred + [x for x in candidate_order if x not in preferred]

        for source_key in candidate_order:
            motif_dict = parsed_sources[source_key]
            if rep_motif in motif_dict:
                title = f"{family_id} | {rep_motif} | {source_key[0]}-{source_key[1]}"
                draw_logo(motif_dict[rep_motif], title, out_png)
                found = True
                found_source = f"{source_key[0]}_{source_key[1]}"
                break

        result_rows.append({
            "Family_ID": family_id,
            "Representative_motif": rep_motif,
            "Logo_found": found,
            "Source_used": found_source,
            "Logo_file": out_png.name if found else "",
        })

        if found:
            print(f"OK  {family_id}: {rep_motif} -> {out_png.name}")
        else:
            print(f"MISS {family_id}: {rep_motif} not found in MEME/STREME files")

    out_table = pd.DataFrame(result_rows)
    out_csv = OUT_DIR / "motif_logo_generation_summary.csv"
    out_xlsx = OUT_DIR / "motif_logo_generation_summary.xlsx"

    out_table.to_csv(out_csv, index=False, sep=";", encoding="utf-8-sig")
    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        out_table.to_excel(writer, index=False, sheet_name="logos")

    print("\nSaved:")
    print("-", out_csv)
    print("-", out_xlsx)
    print("-", OUT_DIR)
    print("\nDone.\n")

if __name__ == "__main__":
    main()