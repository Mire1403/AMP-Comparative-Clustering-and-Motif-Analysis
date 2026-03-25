from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# =========================================================
# CONFIG
# =========================================================

INPUT_FASTA = Path(
    "results/statistics/08_colabfold_input/colabfold_input_selected.fasta"
)

OUTDIR = Path("results/statistics/09_colabfold_output")
BESTDIR = OUTDIR / "best_pdbs"
SUMMARY_CSV = OUTDIR / "best_pdb_summary.csv"

# Extra arguments for colabfold_batch if needed
COLABFOLD_EXTRA_ARGS: List[str] = [
    # "--num-recycle", "6",
]


# =========================================================
# HELPERS
# =========================================================

def ensure_exists(path: Path, kind: str = "path") -> None:
    if not path.exists():
        raise FileNotFoundError(f"{kind.capitalize()} not found: {path}")


def run_colabfold(input_fasta: Path, outdir: Path, extra_args: List[str]) -> None:
    cmd = ["colabfold_batch", str(input_fasta), str(outdir), *extra_args]

    print("\nRUNNING COLABFOLD\n")
    print("Command:")
    print(" ".join(cmd))
    print()

    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError as e:
        raise RuntimeError(
            "colabfold_batch was not found in PATH. Activate your ColabFold environment first."
        ) from e
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"colabfold_batch failed with exit code {e.returncode}") from e


def parse_seq_ids_from_fasta(path: Path) -> List[str]:
    seq_ids: List[str] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip()
                seq_id = header.split("|")[0].strip()
                seq_ids.append(seq_id)
    return seq_ids


def collect_all_pdbs(outdir: Path) -> List[Path]:
    return sorted(outdir.rglob("*.pdb"))


def collect_all_jsons(outdir: Path) -> List[Path]:
    return sorted(outdir.rglob("*.json"))


def safe_read_json(path: Path) -> Optional[dict]:
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return None


def score_pdb_filename(pdb_path: Path) -> tuple:
    """
    Lower tuple = better candidate.

    Last-resort fallback based on filename only.
    """
    name = pdb_path.name.lower()

    rank_score = 99
    if "rank_001" in name:
        rank_score = 1
    elif "rank_1" in name:
        rank_score = 2
    elif "ranked_0" in name:
        rank_score = 3
    elif "best" in name:
        rank_score = 4
    elif "rank" in name:
        rank_score = 10

    relax_score = 1
    if "relaxed" in name:
        relax_score = 0
    elif "unrelaxed" in name:
        relax_score = 1

    return (rank_score, relax_score, len(name), name)


# =========================================================
# RANKING-BASED MODEL SELECTION
# =========================================================

def find_seq_pdbs(seq_id: str, pdb_paths: List[Path]) -> List[Path]:
    return [p for p in pdb_paths if p.name.startswith(seq_id)]


def find_seq_jsons(seq_id: str, json_paths: List[Path]) -> List[Path]:
    return [p for p in json_paths if p.name.startswith(seq_id)]


def choose_best_from_ranking_debug(
    seq_id: str,
    seq_pdbs: List[Path],
    seq_jsons: List[Path],
) -> Tuple[Optional[Path], Dict[str, str]]:
    """
    Prefer ranking_debug.json when available.

    ranking_debug.json usually contains:
      - "order": list of model names in best-to-worst order
      - score mappings like "plddts", "iptm+ptm", etc.

    We try to map the top-ranked model to the corresponding PDB.
    """
    ranking_files = [p for p in seq_jsons if "ranking_debug" in p.name.lower()]
    if not ranking_files:
        return None, {}

    ranking_data = safe_read_json(ranking_files[0])
    if not ranking_data or "order" not in ranking_data:
        return None, {}

    order = ranking_data.get("order", [])
    if not order:
        return None, {}

    top_model = str(order[0])

    # Try to find a PDB whose filename contains the top model identifier
    candidates = [p for p in seq_pdbs if top_model in p.name]
    if candidates:
        best = sorted(candidates, key=score_pdb_filename)[0]
        return best, {
            "selection_method": "ranking_debug.json",
            "top_model_key": top_model,
            "ranking_debug_file": str(ranking_files[0]),
        }

    # Fallback: ranked_0 / rank_001 if present
    ranked_candidates = [
        p for p in seq_pdbs
        if "ranked_0" in p.name.lower() or "rank_001" in p.name.lower()
    ]
    if ranked_candidates:
        best = sorted(ranked_candidates, key=score_pdb_filename)[0]
        return best, {
            "selection_method": "ranking_debug.json_fallback_to_ranked0",
            "top_model_key": top_model,
            "ranking_debug_file": str(ranking_files[0]),
        }

    return None, {}


def extract_numeric_score(score_data: dict) -> Optional[float]:
    """
    Prefer multimer-style ranking metrics if present, otherwise pTM, otherwise pLDDT.
    """
    preferred_keys = [
        "iptm+ptm",
        "iptm",
        "ptm",
        "mean_plddt",
        "plddt",
    ]

    for key in preferred_keys:
        if key in score_data:
            value = score_data[key]
            if isinstance(value, (int, float)):
                return float(value)
    return None


def choose_best_from_score_jsons(
    seq_id: str,
    seq_pdbs: List[Path],
    seq_jsons: List[Path],
) -> Tuple[Optional[Path], Dict[str, str]]:
    """
    If ranking_debug.json is missing, inspect score JSONs written by ColabFold.
    We choose the model with the highest available confidence score and then
    map that model identifier back to a PDB filename.
    """
    score_candidates = []
    for jpath in seq_jsons:
        jname = jpath.name.lower()
        if "ranking_debug" in jname:
            continue

        data = safe_read_json(jpath)
        if not isinstance(data, dict):
            continue

        score = extract_numeric_score(data)
        if score is None:
            continue

        score_candidates.append((jpath, data, score))

    if not score_candidates:
        return None, {}

    # Highest score first
    score_candidates.sort(key=lambda x: x[2], reverse=True)
    best_json, best_data, best_score = score_candidates[0]

    # Try to match by stem fragments first
    stem = best_json.stem
    matching = [
        p for p in seq_pdbs
        if stem in p.stem or p.stem in stem
    ]

    if not matching:
        # Try looser overlap on model-like tokens
        tokens = [tok for tok in stem.replace("-", "_").split("_") if tok]
        matching = [
            p for p in seq_pdbs
            if any(tok in p.stem for tok in tokens)
        ]

    if matching:
        best = sorted(matching, key=score_pdb_filename)[0]
        return best, {
            "selection_method": "scores_json",
            "score_json_file": str(best_json),
            "score_value": str(best_score),
        }

    return None, {}


def choose_best_pdb_for_seq(
    seq_id: str,
    pdb_paths: List[Path],
    json_paths: List[Path],
) -> Tuple[Optional[Path], Dict[str, str]]:
    seq_pdbs = find_seq_pdbs(seq_id, pdb_paths)
    seq_jsons = find_seq_jsons(seq_id, json_paths)

    if not seq_pdbs:
        return None, {"selection_method": "no_pdb_found"}

    # 1) Best: official ranking_debug.json
    best, meta = choose_best_from_ranking_debug(seq_id, seq_pdbs, seq_jsons)
    if best is not None:
        return best, meta

    # 2) Next: score JSONs
    best, meta = choose_best_from_score_jsons(seq_id, seq_pdbs, seq_jsons)
    if best is not None:
        return best, meta

    # 3) Last resort: filename heuristic
    best = sorted(seq_pdbs, key=score_pdb_filename)[0]
    return best, {"selection_method": "filename_fallback"}


# =========================================================
# COPY AND SUMMARIZE
# =========================================================

def copy_best_pdbs(
    seq_ids: List[str],
    outdir: Path,
    bestdir: Path,
) -> List[Dict[str, str]]:
    bestdir.mkdir(parents=True, exist_ok=True)

    all_pdbs = collect_all_pdbs(outdir)
    all_jsons = collect_all_jsons(outdir)

    print(f"Total PDB files found: {len(all_pdbs)}")
    print(f"Total JSON files found: {len(all_jsons)}")

    rows: List[Dict[str, str]] = []

    for seq_id in seq_ids:
        best_pdb, meta = choose_best_pdb_for_seq(seq_id, all_pdbs, all_jsons)

        if best_pdb is None:
            rows.append(
                {
                    "seq_id": seq_id,
                    "best_pdb_filename": "",
                    "best_pdb_source_path": "",
                    "best_pdb_copied_path": "",
                    "selection_method": meta.get("selection_method", "not_found"),
                    "score_json_file": "",
                    "score_value": "",
                    "ranking_debug_file": "",
                    "top_model_key": "",
                    "status": "not_found",
                }
            )
            continue

        target_name = f"{seq_id}_best.pdb"
        target_path = bestdir / target_name
        shutil.copy2(best_pdb, target_path)

        rows.append(
            {
                "seq_id": seq_id,
                "best_pdb_filename": best_pdb.name,
                "best_pdb_source_path": str(best_pdb),
                "best_pdb_copied_path": str(target_path),
                "selection_method": meta.get("selection_method", ""),
                "score_json_file": meta.get("score_json_file", ""),
                "score_value": meta.get("score_value", ""),
                "ranking_debug_file": meta.get("ranking_debug_file", ""),
                "top_model_key": meta.get("top_model_key", ""),
                "status": "copied",
            }
        )

    return rows


def write_summary_csv(rows: List[Dict[str, str]], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "seq_id",
        "best_pdb_filename",
        "best_pdb_source_path",
        "best_pdb_copied_path",
        "selection_method",
        "score_json_file",
        "score_value",
        "ranking_debug_file",
        "top_model_key",
        "status",
    ]

    with open(out_csv, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


# =========================================================
# MAIN
# =========================================================

def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    BESTDIR.mkdir(parents=True, exist_ok=True)

    ensure_exists(INPUT_FASTA, "input FASTA")

    seq_ids = parse_seq_ids_from_fasta(INPUT_FASTA)
    print(f"Sequences in input FASTA: {len(seq_ids)}")

    if not seq_ids:
        raise ValueError(f"No sequences found in FASTA: {INPUT_FASTA}")

    run_colabfold(INPUT_FASTA, OUTDIR, COLABFOLD_EXTRA_ARGS)

    rows = copy_best_pdbs(seq_ids, OUTDIR, BESTDIR)
    write_summary_csv(rows, SUMMARY_CSV)

    n_copied = sum(1 for r in rows if r["status"] == "copied")
    n_missing = sum(1 for r in rows if r["status"] == "not_found")

    print("\nSUMMARY")
    print(f"Best PDBs copied: {n_copied}")
    print(f"Missing best PDBs: {n_missing}")

    print("\nOutputs:")
    print(f"  - Full ColabFold output: {OUTDIR}")
    print(f"  - Best PDBs: {BESTDIR}")
    print(f"  - Summary CSV: {SUMMARY_CSV}")
    print("\nDone.")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        sys.exit(1)