from pathlib import Path
import subprocess

# =====================================================
# PATH CONFIG (RELATIVE TO REPO)
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

DATA_INTERMEDIATE_DIR = PROJECT_ROOT / "data" / "intermediate"
RESULTS_CLUSTER_DIR = PROJECT_ROOT / "results" / "clustering" / "non_amp"

RESULTS_CLUSTER_DIR.mkdir(parents=True, exist_ok=True)

INPUT_FASTA = DATA_INTERMEDIATE_DIR / "nonamp_clean_min5.fasta"

CDHIT_OUTPUT = RESULTS_CLUSTER_DIR / "NONAMP_cdhit80.fasta"

MMSEQ_TMP = RESULTS_CLUSTER_DIR / "mmseq_tmp"
MMSEQ_OUTPUT = RESULTS_CLUSTER_DIR / "NONAMP_mmseq80"

IDENTITY_THRESHOLD = 0.8

# =====================================================
# 1️⃣ RUN CD-HIT
# =====================================================

def run_cdhit():

    if not INPUT_FASTA.exists():
        print(f"Input FASTA not found: {INPUT_FASTA}")
        return

    cmd = [
        "cd-hit",
        "-i", str(INPUT_FASTA),
        "-o", str(CDHIT_OUTPUT),
        "-c", str(IDENTITY_THRESHOLD),
        "-n", "5",
        "-aS", "1.0",
        "-T", "4",
        "-M", "16000"
    ]

    subprocess.run(cmd, check=True)

    print("CD-HIT clustering finished.")
    print("Output:", CDHIT_OUTPUT)


# =====================================================
# 2️⃣ RUN MMSEQS
# =====================================================

def run_mmseq():

    cmd = [
        "mmseqs", "easy-cluster",
        str(INPUT_FASTA),
        str(MMSEQ_OUTPUT),
        str(MMSEQ_TMP),
        "--min-seq-id", str(IDENTITY_THRESHOLD),
        "-c", str(IDENTITY_THRESHOLD),
        "--cov-mode", "1"
    ]

    subprocess.run(cmd, check=True)

    print("MMseqs clustering finished.")
    print("Output prefix:", MMSEQ_OUTPUT)


# =====================================================
# MAIN
# =====================================================

if __name__ == "__main__":

    run_cdhit()
    run_mmseq()

    print("\nNon-AMP clustering pipeline completed successfully.\n")