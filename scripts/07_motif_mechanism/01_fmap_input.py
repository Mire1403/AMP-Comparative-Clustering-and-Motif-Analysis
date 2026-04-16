import pandas as pd
import os
import re

# ==============================
# CONFIGURATION
# ==============================

MOTIF_FILE = "results/statistics/05b_motif_family_analysis_tomtom/motif_family_summary_tomtom.csv"
AMP_FASTA = "data/final/AMP_MASTER.fasta"
OUTDIR = "results/statistics/08_fmap_input"

SEQS_PER_MOTIF = 10
MAX_LEN = 35
MIN_IDENTITY_DIFF = 0.4

TEMP = 300
PH = 7.0
SYS = "bact"

# ==============================
# HELPERS
# ==============================

def read_fasta(path):
    header, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    yield header, "".join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header:
            yield header, "".join(seq)

def motif_to_regex(motif):
    motif = motif.replace("J", "[LI]")
    motif = motif.replace("X", ".")
    return motif

def sequence_identity(s1, s2):
    length = min(len(s1), len(s2))
    matches = sum(a == b for a, b in zip(s1[:length], s2[:length]))
    return matches / length

def is_diverse(seq, selected):
    for s in selected:
        if sequence_identity(seq, s) > (1 - MIN_IDENTITY_DIFF):
            return False
    return True

def clean_sequence(seq):
    valid = set("ACDEFGHIKLMNPQRSTVWY")
    return "".join([aa if aa in valid else "L" for aa in seq])

# ==============================
# LOAD MOTIFS
# ==============================

def load_motifs(path):
    df = pd.read_csv(path, sep=None, engine="python")

    df.columns = (
        df.columns
        .str.strip()
        .str.replace(" ", "_")
        .str.replace("\ufeff", "")
    )

    print("\nDetected columns:")
    print(df.columns.tolist())

    family_col = [c for c in df.columns if "Family" in c][0]
    fdr_col = [c for c in df.columns if "FDR" in c][0]
    er_col = [c for c in df.columns if "ER" in c and "mean" in c.lower()][0]
    members_col = [c for c in df.columns if "Family_members" in c][0]

    return df, family_col, fdr_col, er_col, members_col

# ==============================
# MAIN
# ==============================

def main():

    print("Loading motif table...")
    motifs_df, family_col, fdr_col, er_col, members_col = load_motifs(MOTIF_FILE)

    # Filtering (more permissive ER threshold)
    motifs_df = motifs_df[
        (motifs_df[fdr_col] <= 1e-3) &
        (motifs_df[er_col] >= 2)
    ]

    print(f"\nFamilies after filtering: {len(motifs_df)}")

    print("\nLoading AMP database...")
    amp_seqs = list(read_fasta(AMP_FASTA))

    results = []
    fmap_entries = []

    seq_counter = 1

    print("\nSelecting sequences...")

    for _, row in motifs_df.iterrows():

        family = row[family_col]
        motifs = str(row[members_col]).split(";")

        for motif in motifs:

            motif = motif.strip()

            if len(motif) < 6:
                continue

            # Skip overly degenerate motifs
            if motif.count("X") > len(motif) * 0.5:
                continue

            pattern = motif_to_regex(motif)
            selected = []

            for header, seq in amp_seqs:

                if len(seq) > MAX_LEN:
                    continue

                match = re.search(pattern, seq)
                if not match:
                    continue

                seq = clean_sequence(seq)

                if not is_diverse(seq, selected):
                    continue

                peptide_id = f"P{seq_counter:04d}"

                selected.append(seq)

                results.append({
                    "Peptide_ID": peptide_id,
                    "Family_ID": family,
                    "Motif": motif,
                    "Motif_Length": len(motif),
                    "Motif_Start": match.start(),
                    "Motif_End": match.end(),
                    "Sequence": seq,
                    "Header": header
                })

                fmap_entries.append((peptide_id, seq))
                seq_counter += 1

                if len(selected) == SEQS_PER_MOTIF:
                    break

    # ==============================
    # SAVE METADATA TABLE
    # ==============================

    os.makedirs(OUTDIR, exist_ok=True)

    df = pd.DataFrame(results)
    df.to_csv(f"{OUTDIR}/motif_selected_sequences.csv", index=False)

    # ==============================
    # GENERATE FMAP INPUT FILE
    # ==============================

    fmap_path = f"{OUTDIR}/fmap_input.txt"

    with open(fmap_path, "w") as f:

        f.write("  1\n")  # Required header for FMAP

        for peptide_id, seq in fmap_entries:
            f.write(f"{TEMP} {PH} {SYS} {peptide_id}\n")
            f.write(f"{seq}\n")
            f.write("  0    0.0\n")

    print(f"\nFMAP input file created: {fmap_path}")
    print(f"Total peptides: {len(fmap_entries)}")

# ==============================
# RUN
# ==============================

if __name__ == "__main__":
    main()