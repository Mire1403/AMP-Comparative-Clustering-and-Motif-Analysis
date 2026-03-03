from __future__ import annotations
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[3]

DATA_RAW_DIR = PROJECT_ROOT / "data" / "raw"
DATA_INTERMEDIATE_DIR = PROJECT_ROOT / "data" / "intermediate"
DATA_FINAL_DIR = PROJECT_ROOT / "data" / "final"

RESULTS_DIR = PROJECT_ROOT / "results" / "amp_db_build"

DATA_INTERMEDIATE_DIR.mkdir(parents=True, exist_ok=True)
DATA_FINAL_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

DB_PATHS = {
    "CAMP": DATA_RAW_DIR / "CAMP" / "CAMP_2026-02-16 15-21-23.txt",
    "DBAASP": DATA_RAW_DIR / "DBAASP" / "peptides.csv",
    "dbAMP3": DATA_RAW_DIR / "dbAMP3" / "dbAMP3_pepinfo.xlsx",
    "DRAMP": DATA_RAW_DIR / "DRAMP" / "natural_amps.txt",
}

MIN_LENGTH = 5
MAX_LENGTH = 80

# MIC filter
MIC_THRESHOLD_UGML = 65.0

# Only natural 20 AA
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

UNWANTED_NAME_PATTERNS = [
    "synthetic",
    "designed",
    "construct",
    "analog",
    "mutant",
    "fragment",
    "truncated",
]