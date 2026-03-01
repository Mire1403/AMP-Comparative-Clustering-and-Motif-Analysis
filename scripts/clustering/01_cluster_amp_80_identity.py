from pathlib import Path
import pandas as pd
import subprocess

# =====================================================
# PATH CONFIG (RELATIVE TO REPO)
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

DATA_FINAL_DIR = PROJECT_ROOT / "data" / "final"
RESULTS_CLUSTER_DIR = PROJECT_ROOT / "results" / "clustering"

RESULTS_CLUSTER_DIR.mkdir(parents=True, exist_ok=True)

MASTER_FILE = DATA_FINAL_DIR / "DB_MASTER_CLEAN.xlsx"

FASTA_FILE = RESULTS_CLUSTER_DIR / "AMP_MASTER.fasta"

CDHIT_OUTPUT = RESULTS_CLUSTER_DIR / "AMP_MASTER_cdhit80.fasta"
MMSEQ_TMP = RESULTS_CLUSTER_DIR / "mmseq_tmp"
MMSEQ_OUTPUT = RESULTS_CLUSTER_DIR / "AMP_MASTER_mmseq80"

IDENTITY_THRESHOLD = 0.8

# =====================================================
# 1️⃣ CONVERT TO FASTA
# =====================================================

def convert_to_fasta():

    if not MASTER_FILE.exists():
        print(f"Master file not found: {MASTER_FILE}")
        return

    df = pd.read_excel(MASTER_FILE)

    sequences = df["Sequence"].dropna().unique()

    with open(FASTA_FILE, "w") as f:
        for i, seq in enumerate(sequences, start=1):
            f.write(f">AMP_{i}\n")
            f.write(f"{seq}\n")

    print("FASTA created:", FASTA_FILE)
    print("Total sequences:", len(sequences))


# =====================================================
# 2️⃣ RUN CD-HIT
# =====================================================

def run_cdhit():

    cmd = [
        "cd-hit",
        "-i", str(FASTA_FILE),
        "-o", str(CDHIT_OUTPUT),
        "-c", str(IDENTITY_THRESHOLD),
        "-n", "5",
        "-aS", "1.0",
        "-T", "4",
        "-M", "16000"
    ]

    subprocess.run(cmd)

    print("CD-HIT clustering finished.")


# =====================================================
# 3️⃣ RUN MMSEQS
# =====================================================

def run_mmseq():

    cmd = [
        "mmseqs", "easy-cluster",
        str(FASTA_FILE),
        str(MMSEQ_OUTPUT),
        str(MMSEQ_TMP),
        "--min-seq-id", str(IDENTITY_THRESHOLD),
        "-c", str(IDENTITY_THRESHOLD),
        "--cov-mode", "1"
    ]

    subprocess.run(cmd)

    print("MMseqs clustering finished.")


# =====================================================
# MAIN
# =====================================================

if __name__ == "__main__":

    convert_to_fasta()
    run_cdhit()
    run_mmseq()

    print("\nClustering pipeline completed successfully.\n")