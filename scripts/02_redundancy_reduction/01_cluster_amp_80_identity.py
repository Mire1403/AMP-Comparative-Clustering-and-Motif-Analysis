from __future__ import annotations

from pathlib import Path
import subprocess
import pandas as pd

# =====================================================
# PATH CONFIG (RELATIVE TO REPO)
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

DATA_FINAL_DIR = PROJECT_ROOT / "data" / "final"
RESULTS_CLUSTER_DIR = PROJECT_ROOT / "results" / "clustering"
RESULTS_CLUSTER_DIR.mkdir(parents=True, exist_ok=True)

MASTER_PARQUET = DATA_FINAL_DIR / "DB_MASTER_CLEAN.parquet"

FASTA_FILE = RESULTS_CLUSTER_DIR / "AMP_MASTER.fasta"

CDHIT_OUTPUT = RESULTS_CLUSTER_DIR / "AMP_MASTER_cdhit.fasta"
MMSEQ_TMP = RESULTS_CLUSTER_DIR / "mmseq_tmp_amp"
MMSEQ_OUTPUT = RESULTS_CLUSTER_DIR / "AMP_MASTER_mmseq"

IDENTITY_THRESHOLD = 0.8

# =====================================================
# 1) CONVERT TO FASTA
# =====================================================

def convert_to_fasta() -> None:
    if not MASTER_PARQUET.exists():
        raise FileNotFoundError(f"Master file not found: {MASTER_PARQUET}")

    df = pd.read_parquet(MASTER_PARQUET)

    if "sequence" not in df.columns:
        raise KeyError("Expected column 'sequence' not found in DB_MASTER_CLEAN.parquet")

    sequences = (
        df["sequence"]
        .astype("string")
        .dropna()
        .str.upper()
        .str.strip()
    )
    sequences = sequences[sequences != ""].unique()

    with open(FASTA_FILE, "w", encoding="utf-8") as f:
        for i, seq in enumerate(sequences, start=1):
            f.write(f">AMP_{i:06d}\n{seq}\n")

    print(f"✅ FASTA created: {FASTA_FILE}")
    print(f"   Total sequences: {len(sequences)}")


# =====================================================
# 2) RUN CD-HIT
# =====================================================

def run_cdhit() -> None:
    cmd = [
        "cd-hit",
        "-i", str(FASTA_FILE),
        "-o", str(CDHIT_OUTPUT),
        "-c", str(IDENTITY_THRESHOLD),
        "-n", "5",
        "-aS", "1.0",
        "-T", "4",
        "-M", "16000",
    ]
    subprocess.run(cmd, check=True)
    print(f"✅ CD-HIT done: {CDHIT_OUTPUT}")


# =====================================================
# 3) RUN MMSEQS
# =====================================================

def run_mmseq() -> None:
    cmd = [
        "mmseqs", "easy-cluster",
        str(FASTA_FILE),
        str(MMSEQ_OUTPUT),
        str(MMSEQ_TMP),
        "--min-seq-id", str(IDENTITY_THRESHOLD),
        "-c", str(IDENTITY_THRESHOLD),
        "--cov-mode", "1",
    ]
    subprocess.run(cmd, check=True)
    print(f"✅ MMseqs done: {MMSEQ_OUTPUT}*")


if __name__ == "__main__":
    convert_to_fasta()
    run_cdhit()
    run_mmseq()
    print("\n✅ AMP clustering pipeline completed.\n")