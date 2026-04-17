from __future__ import annotations

from pathlib import Path
import subprocess
import sys

print("\nRUNNING FULL FIMO PIPELINE \n")

# =====================================================
# PATH CONFIG
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

RESULTS_CLUSTER_DIR = PROJECT_ROOT / "results" / "02_clustering"
RESULTS_MOTIF_DIR   = PROJECT_ROOT / "results" / "04_motif_discovery"
RESULTS_BG_DIR      = PROJECT_ROOT / "results" / "03_background" 
RESULTS_FIMO_DIR    = PROJECT_ROOT / "results" / "05_motif_scanning"
RESULTS_FIMO_DIR.mkdir(parents=True, exist_ok=True)

THRESH = "1e-4"

# =====================================================
# INPUT FILES (coherent with script 01 outputs)
# =====================================================

# Motif files (outputs of MEME/STREME)
paths = {
    "cdhit_meme":   RESULTS_MOTIF_DIR / "meme_cdhit"   / "meme.txt",
    "mmseq_meme":   RESULTS_MOTIF_DIR / "meme_mmseq"   / "meme.txt",
    "cdhit_streme": RESULTS_MOTIF_DIR / "streme_cdhit" / "streme.txt",
    "mmseq_streme": RESULTS_MOTIF_DIR / "streme_mmseq" / "streme.txt",
}

# Sequence files
seqs = {
    "cdhit_amps":  RESULTS_CLUSTER_DIR / "AMP_MASTER_cdhit.fasta",
    "mmseq_amps":  RESULTS_CLUSTER_DIR / "AMP_MASTER_mmseq_rep_seq.fasta",

    "cdhit_nonamp": RESULTS_BG_DIR / "nonamp_cdhit_10x_len_kr_matched_progressive_random.fasta",
    "mmseq_nonamp": RESULTS_BG_DIR / "nonamp_mmseq_10x_len_kr_matched_progressive_random.fasta",
}

# =====================================================
# UTILITIES
# =====================================================

def die(msg: str) -> None:
    print(f" {msg}")
    sys.exit(1)

def check_file(path: Path) -> None:
    if not path.exists():
        die(f"File not found: {path}")

def run_fimo(motif_file: Path, seq_file: Path, output_dir: Path) -> None:
    check_file(motif_file)
    check_file(seq_file)

    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "fimo",
        "--thresh", THRESH,
        "--oc", str(output_dir),
        str(motif_file),
        str(seq_file),
    ]

    print("\n===================================")
    print("Running FIMO")
    print("Motifs   :", motif_file)
    print("Sequences:", seq_file)
    print("Output   :", output_dir)
    print("===================================\n")

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError:
        die("FIMO failed. Stopping pipeline.")

# =====================================================
# MAIN
# =====================================================

def main():
    combinations = [
        # CD-HIT + MEME
        (paths["cdhit_meme"],   seqs["cdhit_amps"],    RESULTS_FIMO_DIR / "cdhit" / "meme"   / "fimo_amps"),
        (paths["cdhit_meme"],   seqs["cdhit_nonamp"],  RESULTS_FIMO_DIR / "cdhit" / "meme"   / "fimo_nonamps"),

        # CD-HIT + STREME
        (paths["cdhit_streme"], seqs["cdhit_amps"],    RESULTS_FIMO_DIR / "cdhit" / "streme" / "fimo_amps"),
        (paths["cdhit_streme"], seqs["cdhit_nonamp"],  RESULTS_FIMO_DIR / "cdhit" / "streme" / "fimo_nonamps"),

        # MMseq + MEME
        (paths["mmseq_meme"],   seqs["mmseq_amps"],    RESULTS_FIMO_DIR / "mmseq" / "meme"   / "fimo_amps"),
        (paths["mmseq_meme"],   seqs["mmseq_nonamp"],  RESULTS_FIMO_DIR / "mmseq" / "meme"   / "fimo_nonamps"),

        # MMseq + STREME
        (paths["mmseq_streme"], seqs["mmseq_amps"],    RESULTS_FIMO_DIR / "mmseq" / "streme" / "fimo_amps"),
        (paths["mmseq_streme"], seqs["mmseq_nonamp"],  RESULTS_FIMO_DIR / "mmseq" / "streme" / "fimo_nonamps"),
    ]

    for motif, seq, outdir in combinations:
        run_fimo(motif, seq, outdir)

    print("\n All FIMO runs completed successfully.\n")

if __name__ == "__main__":
    main()