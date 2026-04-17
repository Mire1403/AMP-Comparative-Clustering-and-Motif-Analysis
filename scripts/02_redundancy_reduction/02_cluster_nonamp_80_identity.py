from __future__ import annotations

from pathlib import Path
import subprocess

PROJECT_ROOT = Path(__file__).resolve().parents[2]

DATA_INTERMEDIATE_DIR = PROJECT_ROOT / "data" / "intermediate"
RESULTS_CLUSTER_DIR = PROJECT_ROOT / "results" / "02_clustering"
RESULTS_CLUSTER_DIR.mkdir(parents=True, exist_ok=True)

INPUT_FASTA = DATA_INTERMEDIATE_DIR / "nonamp_clean_min5.fasta"

CDHIT_OUTPUT = RESULTS_CLUSTER_DIR / "NONAMP_cdhit.fasta"

MMSEQ_TMP = RESULTS_CLUSTER_DIR / "mmseq_tmp_nonamp"
MMSEQ_OUTPUT = RESULTS_CLUSTER_DIR / "NONAMP_mmseq"

IDENTITY_THRESHOLD = 0.8


def run_cdhit() -> None:
    if not INPUT_FASTA.exists():
        raise FileNotFoundError(f"Input FASTA not found: {INPUT_FASTA}")

    cmd = [
        "cd-hit",
        "-i", str(INPUT_FASTA),
        "-o", str(CDHIT_OUTPUT),
        "-c", str(IDENTITY_THRESHOLD),
        "-n", "5",
        "-aS", "1.0",
        "-T", "4",
        "-M", "16000",
    ]
    subprocess.run(cmd, check=True)
    print(f"✅ CD-HIT done: {CDHIT_OUTPUT}")


def run_mmseq() -> None:
    if not INPUT_FASTA.exists():
        raise FileNotFoundError(f"Input FASTA not found: {INPUT_FASTA}")

    cmd = [
        "mmseqs", "easy-cluster",
        str(INPUT_FASTA),
        str(MMSEQ_OUTPUT),
        str(MMSEQ_TMP),
        "--min-seq-id", str(IDENTITY_THRESHOLD),
        "-c", str(IDENTITY_THRESHOLD),
        "--cov-mode", "1",
    ]
    subprocess.run(cmd, check=True)
    print(f"✅ MMseqs done: {MMSEQ_OUTPUT}*")


if __name__ == "__main__":
    run_cdhit()
    run_mmseq()
    print("\n✅ Non-AMP clustering pipeline completed.\n")