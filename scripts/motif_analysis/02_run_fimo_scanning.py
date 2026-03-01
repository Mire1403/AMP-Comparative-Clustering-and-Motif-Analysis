#!/usr/bin/env python3

from pathlib import Path
import subprocess
import sys

print("\nRUNNING FULL FIMO PIPELINE\n")

# =====================================================
# PATH CONFIG (RELATIVE TO REPO)
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

RESULTS_CLUSTER_DIR = PROJECT_ROOT / "results" / "clustering"
RESULTS_BACKGROUND_DIR = PROJECT_ROOT / "results" / "background_generation"
RESULTS_MOTIF_DIR = PROJECT_ROOT / "results" / "motif_discovery"
RESULTS_FIMO_DIR = PROJECT_ROOT / "results" / "motif_scanning"

RESULTS_FIMO_DIR.mkdir(parents=True, exist_ok=True)

THRESH = "1e-4"

# =====================================================
# INPUT FILES
# =====================================================

# Motif files (output of MEME/STREME)
paths = {
    "cdhit_meme": RESULTS_MOTIF_DIR / "meme_cdhit80" / "meme.txt",
    "mmseq_meme": RESULTS_MOTIF_DIR / "meme_mmseq80" / "meme.txt",
    "cdhit_streme": RESULTS_MOTIF_DIR / "streme_cdhit80" / "streme.txt",
    "mmseq_streme": RESULTS_MOTIF_DIR / "streme_mmseq80" / "streme.txt",
}

# Sequence files
seqs = {
    "cdhit_amps": RESULTS_CLUSTER_DIR / "AMP_MASTER_cdhit80.fasta",
    "mmseq_amps": RESULTS_CLUSTER_DIR / "AMP_MASTER_mmseq80_rep_seq.fasta",
    "cdhit_nonamp": RESULTS_MOTIF_DIR / "nonamp_cdhit80_10x.fasta",
    "mmseq_nonamp": RESULTS_MOTIF_DIR / "nonamp_mmseq80_10x.fasta",
}

# =====================================================
# UTILITIES
# =====================================================

def check_file(path):
    if not path.exists():
        print(f"❌ File not found: {path}")
        sys.exit(1)


def run_fimo(motif_file, seq_file, output_dir):

    check_file(motif_file)
    check_file(seq_file)

    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "fimo",
        "--thresh", THRESH,
        "--oc", str(output_dir),
        str(motif_file),
        str(seq_file)
    ]

    print("\n===================================")
    print("Running FIMO")
    print("Motifs:", motif_file)
    print("Sequences:", seq_file)
    print("Output:", output_dir)
    print("===================================\n")

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError:
        print("❌ FIMO failed. Stopping pipeline.")
        sys.exit(1)

# =====================================================
# MAIN
# =====================================================

def main():

    combinations = [

        # CD-HIT + MEME
        (paths["cdhit_meme"], seqs["cdhit_amps"], RESULTS_FIMO_DIR / "cdhit80" / "meme" / "fimo_amps"),
        (paths["cdhit_meme"], seqs["cdhit_nonamp"], RESULTS_FIMO_DIR / "cdhit80" / "meme" / "fimo_nonamps"),

        # CD-HIT + STREME
        (paths["cdhit_streme"], seqs["cdhit_amps"], RESULTS_FIMO_DIR / "cdhit80" / "streme" / "fimo_amps"),
        (paths["cdhit_streme"], seqs["cdhit_nonamp"], RESULTS_FIMO_DIR / "cdhit80" / "streme" / "fimo_nonamps"),

        # MMseq + MEME
        (paths["mmseq_meme"], seqs["mmseq_amps"], RESULTS_FIMO_DIR / "mmseq80" / "meme" / "fimo_amps"),
        (paths["mmseq_meme"], seqs["mmseq_nonamp"], RESULTS_FIMO_DIR / "mmseq80" / "meme" / "fimo_nonamps"),

        # MMseq + STREME
        (paths["mmseq_streme"], seqs["mmseq_amps"], RESULTS_FIMO_DIR / "mmseq80" / "streme" / "fimo_amps"),
        (paths["mmseq_streme"], seqs["mmseq_nonamp"], RESULTS_FIMO_DIR / "mmseq80" / "streme" / "fimo_nonamps"),
    ]

    for motif, seq, outdir in combinations:
        run_fimo(motif, seq, outdir)

    print("\n✅ All FIMO runs completed successfully.\n")


if __name__ == "__main__":
    main()