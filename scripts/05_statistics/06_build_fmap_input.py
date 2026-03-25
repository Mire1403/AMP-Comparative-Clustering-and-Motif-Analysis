from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd


# =========================================================
# CONFIG / RUTAS
# =========================================================

FAMILIES_CSV = Path(
    "results/statistics/05b_motif_family_analysis_tomtom/motif_family_summary_for_memory_tomtom.csv"
)

AMP_FASTAS = [
    Path("data/final/AMP_MASTER.fasta"),
]

OUTDIR = Path("results/statistics/08_fmap_input")

FMAP_FASTA_OUT = OUTDIR / "fmap_input_selected.fasta"
MAPPING_CSV_OUT = OUTDIR / "fmap_input_selected_mapping.csv"
MANUAL_XLSX_OUT = OUTDIR / "fmap_input_selected_manual.xlsx"
EXCLUDED_CSV_OUT = OUTDIR / "fmap_excluded_too_long_or_not_selected.csv"

FMAP_MODEL = "peptide in membrane"
FMAP_ENVIRONMENT = "Gram-negative bacteria IM"
FMAP_TEMPERATURE_K = 300
FMAP_PH = 7

MIN_FMAP_LENGTH = 14
MAX_FMAP_LENGTH = 35

MAX_MOTIFS_PER_FAMILY = 2
SEQS_PER_MOTIF = 2

HYDROPHOBIC_AA = set("AILMFWVYCG")
BASIC_AA = set("KRH")
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


# =========================================================
# FASTA IO
# =========================================================

def read_fasta(path: Path) -> Iterable[Tuple[str, str]]:
    header = None
    seq_parts: List[str] = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts).upper()
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line)

        if header is not None:
            yield header, "".join(seq_parts).upper()


# =========================================================
# MOTIF HANDLING
# =========================================================

def parse_family_members(cell: str) -> List[str]:
    if not cell or str(cell).strip() == "":
        return []
    return [x.strip().upper() for x in str(cell).split(";") if x.strip()]


def motif_to_regex(motif: str) -> str:
    """
    Convert motif to regex for searching inside peptide sequences.

    Rules:
    J -> [IL]
    X -> .
    B -> [DN]
    Z -> [EQ]
    """
    motif = motif.strip().upper()
    parts: List[str] = []

    for aa in motif:
        if aa == "J":
            parts.append("[IL]")
        elif aa == "X":
            parts.append(".")
        elif aa == "B":
            parts.append("[DN]")
        elif aa == "Z":
            parts.append("[EQ]")
        else:
            parts.append(re.escape(aa))

    return "".join(parts)


def clean_sequence_for_fmap(seq: str) -> Tuple[str, str, bool]:
    """
    Keep only standard amino acids for FMAP output.

    Replace:
      J -> L
      B -> D
      Z -> E
      U -> C
      O -> K
      X -> A

    Any remaining non-standard character is converted to A.

    Returns:
      cleaned_seq
      cleaning_notes
      had_nonstandard_before_cleaning
    """
    seq = seq.strip().upper()

    replacements = {
        "J": "L",
        "B": "D",
        "Z": "E",
        "U": "C",
        "O": "K",
        "X": "A",
    }

    cleaned: List[str] = []
    notes: List[str] = []
    had_nonstandard = False

    for aa in seq:
        original = aa
        if aa not in STANDARD_AA:
            had_nonstandard = True

        aa = replacements.get(aa, aa)
        if aa not in STANDARD_AA:
            aa = "A"
            if original not in replacements:
                notes.append(f"{original}->A")
        else:
            if original != aa:
                notes.append(f"{original}->{aa}")

        cleaned.append(aa)

    unique_notes = ";".join(sorted(set(notes))) if notes else ""
    return "".join(cleaned), unique_notes, had_nonstandard


def count_motif_matches(pattern: re.Pattern, sequence: str) -> int:
    return len(list(pattern.finditer(sequence)))


def compute_fmap_priority(seq: str, motif_match_count: int) -> float:
    """
    Simple priority score to select better FMAP candidates.

    Prioritizes:
      - lengths in a useful range for FMAP
      - higher motif match count
      - moderate hydrophobic/basic content
    Penalizes:
      - very short sequences
    """
    length = len(seq)
    hydrophobic_fraction = sum(aa in HYDROPHOBIC_AA for aa in seq) / length
    basic_fraction = sum(aa in BASIC_AA for aa in seq) / length

    score = 0.0

    if 16 <= length <= 30:
        score += 2.0
    elif 14 <= length < 16:
        score += 1.0
    elif 31 <= length <= 35:
        score += 0.5

    score += min(motif_match_count, 3) * 1.0
    score += hydrophobic_fraction * 2.0
    score += basic_fraction * 1.0

    if length < 14:
        score -= 3.0

    return round(score, 4)


# =========================================================
# LOAD FAMILIES
# =========================================================

def select_representative_motifs(
    rep_motif: str,
    all_motifs: List[str],
    max_motifs: int,
) -> List[str]:
    """
    Select up to max_motifs motifs per family.
    Strategy:
      1. include representative motif first if present
      2. fill with remaining motifs, prioritizing longer motifs
      3. tie-break alphabetically
    """
    unique_motifs = []
    seen = set()
    for m in all_motifs:
        if m and m not in seen:
            unique_motifs.append(m)
            seen.add(m)

    selected: List[str] = []

    if rep_motif and rep_motif in seen:
        selected.append(rep_motif)

    remaining = [m for m in unique_motifs if m != rep_motif]
    remaining = sorted(remaining, key=lambda x: (-len(x), x))

    for motif in remaining:
        if len(selected) >= max_motifs:
            break
        selected.append(motif)

    return selected[:max_motifs]


def load_families_table(families_path: Path) -> Dict[str, Dict]:
    """
    Expected columns in the input file:
    - Family_ID
    - Representative_motif
    - Representative_class (optional)
    - Family_members

    Input file is semicolon-separated.
    """
    families: Dict[str, Dict] = {}

    with open(families_path, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter=";")
        fieldnames = reader.fieldnames or []

        required = {"Family_ID", "Representative_motif", "Family_members"}
        missing = required - set(fieldnames)
        if missing:
            raise ValueError(
                f"Missing required columns: {sorted(missing)}. Found: {fieldnames}"
            )

        for row in reader:
            family_id = row["Family_ID"].strip()
            rep_motif = row["Representative_motif"].strip().upper()
            rep_class = row.get("Representative_class", "").strip()

            motifs: List[str] = []
            if rep_motif:
                motifs.append(rep_motif)

            for m in parse_family_members(row.get("Family_members", "")):
                if m not in motifs:
                    motifs.append(m)

            selected_motifs = select_representative_motifs(
                rep_motif=rep_motif,
                all_motifs=motifs,
                max_motifs=MAX_MOTIFS_PER_FAMILY,
            )

            patterns = [(m, re.compile(motif_to_regex(m))) for m in selected_motifs]

            families[family_id] = {
                "Family_ID": family_id,
                "Representative_motif": rep_motif,
                "Representative_class": rep_class,
                "All_motifs": motifs,
                "Selected_motifs": selected_motifs,
                "Patterns": patterns,
            }

    return families


# =========================================================
# SEARCH CANDIDATES PER MOTIF
# =========================================================

def collect_motif_candidates(
    families: Dict[str, Dict],
    fasta_paths: List[Path],
) -> List[Dict]:
    """
    Returns one row per family-motif-sequence candidate.
    Only records matches for the selected motifs of each family.
    """
    candidates: List[Dict] = []

    for fasta_path in fasta_paths:
        db_name = fasta_path.stem
        print(f"Searching in FASTA: {fasta_path}")

        for seq_header, sequence in read_fasta(fasta_path):
            fmap_seq, cleaning_notes, had_nonstandard = clean_sequence_for_fmap(sequence)
            fmap_len = len(fmap_seq)
            eligible = MIN_FMAP_LENGTH <= fmap_len <= MAX_FMAP_LENGTH

            for _, fam in families.items():
                for motif, pattern in fam["Patterns"]:
                    n_matches = count_motif_matches(pattern, sequence)
                    if n_matches > 0:
                        candidates.append(
                            {
                                "Family_ID": fam["Family_ID"],
                                "Representative_class": fam["Representative_class"],
                                "Representative_motif": fam["Representative_motif"],
                                "Selected_motif": motif,
                                "DB_name": db_name,
                                "Sequence_header": seq_header,
                                "Original_sequence": sequence,
                                "FMAP_sequence": fmap_seq,
                                "FMAP_sequence_length": fmap_len,
                                "Eligible_for_FMAP": eligible,
                                "Motif_match_count": n_matches,
                                "had_nonstandard_before_cleaning": had_nonstandard,
                                "cleaning_notes": cleaning_notes,
                                "FMAP_priority_score": compute_fmap_priority(
                                    fmap_seq, n_matches
                                ),
                            }
                        )

    return candidates


def select_sequences_per_motif(
    candidates: List[Dict],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Keep up to SEQS_PER_MOTIF sequences per (Family_ID, Selected_motif),
    but do not allow the same sequence to be selected more than once across motifs.

    Global uniqueness key:
      - Sequence_header
      - FMAP_sequence

    Ranking:
      1. higher FMAP priority score first
      2. more motif matches first
      3. shorter sequence first
      4. alphabetical header
    """
    if not candidates:
        return pd.DataFrame(), pd.DataFrame()

    df = pd.DataFrame(candidates)

    df_eligible = df[df["Eligible_for_FMAP"]].copy()
    df_not_eligible = df[~df["Eligible_for_FMAP"]].copy()

    # Remove exact duplicated candidate rows
    df_eligible = df_eligible.drop_duplicates(
        subset=["Family_ID", "Selected_motif", "FMAP_sequence", "Sequence_header"]
    ).copy()

    # Sort globally so the best candidates are seen first within each motif
    df_eligible = df_eligible.sort_values(
        by=[
            "Family_ID",
            "Selected_motif",
            "FMAP_priority_score",
            "Motif_match_count",
            "FMAP_sequence_length",
            "Sequence_header",
        ],
        ascending=[True, True, False, False, True, True],
    ).copy()

    selected_rows: List[Dict] = []
    excluded_rows: List[Dict] = []

    used_sequences = set()
    motif_counts: Dict[Tuple[str, str], int] = {}

    for _, row in df_eligible.iterrows():
        key = (row["Family_ID"], row["Selected_motif"])
        seq_key = (row["Sequence_header"], row["FMAP_sequence"])

        current_count = motif_counts.get(key, 0)

        if current_count >= SEQS_PER_MOTIF:
            row_dict = row.to_dict()
            row_dict["exclusion_reason"] = f"beyond_top_{SEQS_PER_MOTIF}_per_motif"
            excluded_rows.append(row_dict)
            continue

        if seq_key in used_sequences:
            row_dict = row.to_dict()
            row_dict["exclusion_reason"] = "sequence_already_selected_for_another_motif"
            excluded_rows.append(row_dict)
            continue

        row_dict = row.to_dict()
        row_dict["motif_seq_rank"] = current_count + 1
        selected_rows.append(row_dict)

        used_sequences.add(seq_key)
        motif_counts[key] = current_count + 1

    selected_df = pd.DataFrame(selected_rows)
    excluded_selected_df = pd.DataFrame(excluded_rows)

    if not df_not_eligible.empty:
        df_not_eligible = df_not_eligible.copy()
        df_not_eligible["exclusion_reason"] = (
            f"sequence_length_not_in_[{MIN_FMAP_LENGTH},{MAX_FMAP_LENGTH}]"
        )

    excluded_frames = []
    if not df_not_eligible.empty:
        excluded_frames.append(df_not_eligible)
    if not excluded_selected_df.empty:
        excluded_frames.append(excluded_selected_df)

    excluded_df = (
        pd.concat(excluded_frames, ignore_index=True)
        if excluded_frames
        else pd.DataFrame()
    )

    return selected_df, excluded_df


# =========================================================
# BUILD MANUAL TABLE
# =========================================================

def build_manual_table(selected_df: pd.DataFrame) -> pd.DataFrame:
    """
    Manual FMAP input table.
    One row per selected sequence per motif.
    """
    if selected_df.empty:
        return pd.DataFrame()

    selected_df = selected_df.copy().reset_index(drop=True)
    selected_df["seq_id"] = [f"seq_{i:05d}" for i in range(1, len(selected_df) + 1)]

    rows = []
    for _, row in selected_df.iterrows():
        seq = row["FMAP_sequence"]
        seq_id = row["seq_id"]
        seq_len = len(seq)

        rows.append(
            {
                "seq_id": seq_id,
                "Family_ID": row["Family_ID"],
                "Representative_class": row["Representative_class"],
                "Representative_motif": row["Representative_motif"],
                "Selected_motif": row["Selected_motif"],
                "motif_seq_rank": row["motif_seq_rank"],
                "motif_match_count": row["Motif_match_count"],
                "fmap_priority_score": row["FMAP_priority_score"],
                "FMAP_sequence": seq,
                "sequence_length": seq_len,
                "DB_name": row["DB_name"],
                "Sequence_header": row["Sequence_header"],
                "Original_sequence": row["Original_sequence"],
                "had_nonstandard_before_cleaning": row["had_nonstandard_before_cleaning"],
                "cleaning_notes": row["cleaning_notes"],
                "copy_ready_fasta": f">{seq_id}\n{seq}",

                # fixed FMAP config
                "model": FMAP_MODEL,
                "environment": FMAP_ENVIRONMENT,
                "temperature_K": FMAP_TEMPERATURE_K,
                "pH": FMAP_PH,

                # FMAP results to fill manually
                "binding_energy_kcal_mol": "",
                "has_helix": "",
                "helix_n": "",
                "helix_start": "",
                "helix_end": "",
                "helix_length": "",
                "helical_fraction": "",
                "stability_water_kcal_mol": "",
                "stability_bound_kcal_mol": "",
                "deltaGtransfer_kcal_mol": "",
                "tilt_angle_deg": "",
                "depth_A": "",

                # workflow
                "output_pdb_downloaded": "",
                "fmap_status": "",
                "fmap_notes": "",
            }
        )

    df_manual = pd.DataFrame(rows)
    df_manual = df_manual.sort_values(
        by=["Family_ID", "Selected_motif", "motif_seq_rank", "seq_id"]
    ).reset_index(drop=True)

    return df_manual


# =========================================================
# OUTPUTS
# =========================================================

def write_mapping_csv(selected_df: pd.DataFrame, out_path: Path) -> None:
    selected_df.to_csv(out_path, index=False)


def write_excluded_csv(excluded_df: pd.DataFrame, out_path: Path) -> None:
    excluded_df.to_csv(out_path, index=False)


def write_selected_fasta(df_manual: pd.DataFrame, out_path: Path) -> None:
    with open(out_path, "w", encoding="utf-8") as f:
        for _, row in df_manual.iterrows():
            header = (
                f"{row['seq_id']}"
                f"|family={row['Family_ID']}"
                f"|motif={row['Selected_motif']}"
                f"|rank={row['motif_seq_rank']}"
            )
            f.write(f">{header}\n{row['FMAP_sequence']}\n")


def write_manual_excel(
    df_manual: pd.DataFrame,
    selected_df: pd.DataFrame,
    excluded_df: pd.DataFrame,
    families: Dict[str, Dict],
    out_path: Path,
) -> None:
    family_summary_rows = []
    for fam_id, fam in families.items():
        family_summary_rows.append(
            {
                "Family_ID": fam_id,
                "Representative_motif": fam["Representative_motif"],
                "Representative_class": fam["Representative_class"],
                "n_all_motifs": len(fam["All_motifs"]),
                "n_selected_motifs": len(fam["Selected_motifs"]),
                "selected_motifs": ";".join(fam["Selected_motifs"]),
            }
        )
    df_family_summary = pd.DataFrame(family_summary_rows).sort_values("Family_ID")

    if out_path.exists():
        out_path.unlink()

    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        df_manual.to_excel(writer, sheet_name="FMAP_manual_input", index=False)
        selected_df.to_excel(writer, sheet_name="Selected_mapping", index=False)
        excluded_df.to_excel(writer, sheet_name="Excluded", index=False)
        df_family_summary.to_excel(writer, sheet_name="Family_motif_summary", index=False)

        ws = writer.sheets["FMAP_manual_input"]
        ws.freeze_panes = "A2"


# =========================================================
# MAIN
# =========================================================

def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    print("\nBUILD FMAP INPUT (priority-based selection for FMAP)\n")
    print(f"Families CSV: {FAMILIES_CSV}")
    print(f"N FASTA DBs:  {len(AMP_FASTAS)}")
    for p in AMP_FASTAS:
        print(f"  - {p}")
    print(f"Output dir:   {OUTDIR}")
    print(f"FMAP length window: [{MIN_FMAP_LENGTH}, {MAX_FMAP_LENGTH}]")
    print(f"Max motifs/family: {MAX_MOTIFS_PER_FAMILY}")
    print(f"Seqs per motif: {SEQS_PER_MOTIF}")
    print("Unique selected sequences across motifs: True\n")

    if not FAMILIES_CSV.exists():
        raise FileNotFoundError(f"Families CSV not found: {FAMILIES_CSV}")

    missing_fastas = [p for p in AMP_FASTAS if not p.exists()]
    if missing_fastas:
        raise FileNotFoundError(
            "The following FASTA files were not found:\n"
            + "\n".join(f"  - {p}" for p in missing_fastas)
        )

    families = load_families_table(FAMILIES_CSV)
    print(f"Loaded families: {len(families)}")

    candidates = collect_motif_candidates(families, AMP_FASTAS)
    print(f"Total family-motif-sequence candidates: {len(candidates)}")

    if not candidates:
        raise ValueError("No motif candidates found. Check your families and FASTA files.")

    selected_df, excluded_df = select_sequences_per_motif(candidates)

    if selected_df.empty:
        raise ValueError(
            "No valid sequences remained after applying FMAP filters and ranking."
        )

    df_manual = build_manual_table(selected_df)

    write_mapping_csv(selected_df, MAPPING_CSV_OUT)
    write_selected_fasta(df_manual, FMAP_FASTA_OUT)
    write_manual_excel(df_manual, selected_df, excluded_df, families, MANUAL_XLSX_OUT)
    write_excluded_csv(excluded_df, EXCLUDED_CSV_OUT)

    n_selected_families = selected_df["Family_ID"].nunique()
    n_selected_motifs = selected_df[["Family_ID", "Selected_motif"]].drop_duplicates().shape[0]
    n_unique_sequences = selected_df[["Sequence_header", "FMAP_sequence"]].drop_duplicates().shape[0]

    print(f"Selected families: {n_selected_families}")
    print(f"Selected motifs: {n_selected_motifs}")
    print(f"Selected sequence rows for FMAP: {len(df_manual)}")
    print(f"Unique selected sequences: {n_unique_sequences}")
    print(f"Excluded rows: {len(excluded_df)}")

    print("\nOutputs:")
    print(f"  - {FMAP_FASTA_OUT}")
    print(f"  - {MAPPING_CSV_OUT}")
    print(f"  - {MANUAL_XLSX_OUT}")
    print(f"  - {EXCLUDED_CSV_OUT}")
    print("\nDone.")


if __name__ == "__main__":
    main()