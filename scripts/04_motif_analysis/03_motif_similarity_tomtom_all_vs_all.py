"""
Build a global MEME-format motif file and run Tomtom all-vs-all.

Pipeline:
1) Concatenate MEME + STREME motifs (CD-HIT & MMseq)
2) Prefix motif IDs to preserve provenance
3) Run Tomtom all-vs-all
4) Save outputs in results/motif_similarity/

Requires:
- MEME Suite installed
- 'tomtom' available in PATH
"""

from pathlib import Path
import re
import subprocess
import shutil
import sys


# =====================================================
# REPO ROOT
# =====================================================

def find_repo_root(start: Path) -> Path:
    markers = [".git", "README.md", "requirements.txt", "pyproject.toml", "environment.yml"]
    start = start.resolve()
    for parent in [start] + list(start.parents):
        if any((parent / m).exists() for m in markers):
            return parent
    return start.parents[2]


PROJECT_ROOT = find_repo_root(Path(__file__).parent)

MOTIF_DIR = PROJECT_ROOT / "results" / "motif_discovery"
SIM_DIR = PROJECT_ROOT / "results" / "motif_similarity"
GLOBAL_DIR = SIM_DIR / "global_motifs"
TOMTOM_DIR = SIM_DIR / "tomtom_all_vs_all"

GLOBAL_DIR.mkdir(parents=True, exist_ok=True)
SIM_DIR.mkdir(parents=True, exist_ok=True)

GLOBAL_FILE = GLOBAL_DIR / "all_motifs_global.meme"


# =====================================================
# INPUTS
# =====================================================

SOURCES = {
    "meme_cdhit":   MOTIF_DIR / "meme_cdhit" / "meme.txt",
    "meme_mmseq":   MOTIF_DIR / "meme_mmseq" / "meme.txt",
    "streme_cdhit": MOTIF_DIR / "streme_cdhit" / "streme.txt",
    "streme_mmseq": MOTIF_DIR / "streme_mmseq" / "streme.txt",
}


# =====================================================
# PARSING
# =====================================================

def extract_header(text: str) -> str | None:
    m = re.search(r"(MEME version.*?)(?=\nMOTIF\s)", text, flags=re.S)
    return m.group(1).strip() if m else None


def extract_motif_blocks(text: str) -> list[str]:
    blocks = re.findall(
        r"(MOTIF\s+.*?)(?=\nMOTIF\s|\Z)",
        text,
        flags=re.S
    )
    blocks = [b.strip() for b in blocks if "letter-probability matrix" in b]
    return blocks


def prefix_motif_ids(block: str, prefix: str) -> str:
    lines = block.splitlines()
    if lines and lines[0].startswith("MOTIF"):
        parts = lines[0].split()
        if len(parts) >= 2:
            parts[1] = f"{prefix}__{parts[1]}"
            lines[0] = " ".join(parts)
    return "\n".join(lines)


# =====================================================
# TOMTOM
# =====================================================

def check_tomtom():
    try:
        subprocess.run(["tomtom", "-h"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except FileNotFoundError:
        print("❌ Tomtom not found in PATH. Install MEME Suite.")
        sys.exit(1)


def run_tomtom_all_vs_all(motif_file: Path, output_dir: Path):

    if output_dir.exists():
        shutil.rmtree(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "tomtom",
        "-oc", str(output_dir),
        "-dist", "pearson",
        "-thresh", "0.1",
        str(motif_file),
        str(motif_file)
    ]

    print("\nRunning Tomtom all-vs-all:")
    print(" ".join(cmd))

    subprocess.run(cmd, check=True)

    if not (output_dir / "tomtom.tsv").exists():
        print("⚠️  Tomtom finished but tomtom.tsv not found.")


# =====================================================
# MAIN
# =====================================================

def main():

    print("\nBUILDING GLOBAL MOTIF FILE + RUNNING TOMTOM\n")

    header = None
    all_blocks = []
    total_motifs = 0

    for prefix, path in SOURCES.items():

        if not path.exists():
            print(f"⚠️  Missing file: {path}")
            continue

        text = path.read_text(errors="ignore")

        if header is None:
            header = extract_header(text)

        blocks = extract_motif_blocks(text)
        blocks = [prefix_motif_ids(b, prefix) for b in blocks]

        all_blocks.extend(blocks)
        total_motifs += len(blocks)

        print(f"{prefix}: {len(blocks)} motifs loaded")

    if header is None or not all_blocks:
        print("❌ Failed to build global motif file.")
        sys.exit(1)

    GLOBAL_FILE.write_text(header + "\n\n" + "\n\n".join(all_blocks) + "\n")

    print(f"\n✅ Global motif file created:")
    print(f"   {GLOBAL_FILE}")
    print(f"   Total motifs: {total_motifs}")

    # Run Tomtom
    check_tomtom()
    run_tomtom_all_vs_all(GLOBAL_FILE, TOMTOM_DIR)

    print("\n Tomtom all-vs-all completed.")
    print(f"   Results in: {TOMTOM_DIR}\n")


if __name__ == "__main__":
    main()