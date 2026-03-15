# AMP Comparative Clustering & Motif Analysis

Reproducible bioinformatics pipeline to benchmark two redundancy-reduction strategies (**CD-HIT** vs **MMseqs2**) for antimicrobial peptide (**AMP**) motif discovery, enrichment analysis, and motif-family robustness.

The project evaluates whether the choice of clustering strategy systematically affects the motifs recovered downstream and their statistical support.

## Research question

Does the choice of redundancy-reduction strategy (CD-HIT vs MMseqs2) systematically influence the identity, enrichment, and robustness of antimicrobial peptide motifs discovered downstream?

## Overview

Starting from a curated AMP collection and a matched non-AMP background, the pipeline:

1. integrates four public AMP resources into a unified database,
2. prepares a cleaned non-AMP sequence pool,
3. applies independent redundancy reduction with **CD-HIT** and **MMseqs2**,
4. generates matched 10× non-AMP background sets,
5. performs discriminative motif discovery with **MEME** and **STREME**,
6. scans AMP and non-AMP sequence sets with **FIMO**,
7. computes motif enrichment using Fisher’s exact test with Benjamini–Hochberg FDR correction,
8. groups significant motifs into families using sequence similarity and identifies robust families reproduced across multiple pipelines.

The final output is a ranked set of **robust AMP motif families** reproducible across multiple pipeline configurations.

---

## Main outputs

The most relevant outputs for interpretation are located in:

- `results/statistics/03_fimo_enrichment/` — motif enrichment tables
- `results/statistics/04_fimo_reporting/` — cleaned significant motif summaries and plots
- `results/statistics/05_motif_family_analysis/` — motif family assignments and robustness analysis
- `results/statistics/06_final_reporting/` — final ranked motif families and sequence logos

---

## Repository structure

```text
.
├── data/
│   ├── raw/                              # Original source files (not tracked)
│   │   ├── CAMP/
│   │   ├── DBAASP/
│   │   ├── DRAMP/
│   │   ├── dbAMP3/
│   │   └── non_AMPs/
│   ├── intermediate/                     # Processed per-step files (parquet tracked)
│   └── final/
│       ├── AMP_MASTER.fasta              # Final curated AMP FASTA
│       └── DB_MASTER_CLEAN.parquet       # Final curated AMP metadata table
│
├── scripts/
│   ├── 01_dataset_construction/
│   │   ├── amp/
│   │   └── non_amp/
│   ├── 02_redundancy_reduction/
│   ├── 03_background_sampling/
│   ├── 04_motif_analysis/
│   ├── 05_statistics/
│   └── 06_generate_logos/
│
├── results/
│   ├── amp_db_build/
│   ├── clustering/
│   ├── background_generation/            # Not tracked
│   ├── 04_motif_discovery/               # MEME/STREME reports
│   ├── 05_motif_scanning/                # FIMO outputs (not tracked)
│   └── statistics/
│       ├── clustering/
│       ├── background_sampling/
│       ├── 03_fimo_enrichment/
│       ├── 04_fimo_reporting/
│       ├── 05_motif_family_analysis/
│       └── 06_final_reporting/
│
├── environment.yml
├── .gitignore
└── README.md
```

---

## Pipeline

### Stage 1 — Dataset construction

#### **AMP database.**
Scripts in ``scripts/01_dataset_construction/amp/`` integrate four public AMP databases (**CAMP**, **DBAASP**, **dbAMP3**, **DRAMP**) into a standardised master dataset.

Main processing steps:
- harmonisation of heterogeneous source formats and column names,
- sequence normalisation to uppercase,
- removal of missing or duplicated sequences,
- filtering by peptide length and natural amino acid alphabet,
- removal of synthetic/designed peptides by name pattern,
- MIC-based activity filtering,
- final deduplication with metadata merging across sources.

Final outputs:
- ``data/final/DB_MASTER_CLEAN.parquet``
- ``data/final/AMP_MASTER.fasta``

### **Non-AMP background.** 
``scripts/01_dataset_construction/non_amp/01_dataset_construction_nonamp.py`` filters a UniProt-derived non-AMP FASTA by:
- minimum peptide length,
- valid amino acid alphabet,
- sequence uniqueness.

This produces the cleaned non-AMP pool used for background sampling.

### Stage 2 — Redundancy reduction 
Scripts in ``scripts/02_redundancy_reduction/`` cluster AMP and non-AMP sequence sets independently at 80% sequence identity using:
- **CD-HIT**
- **MMseqs2**

This yields representative sequence sets for both clustering strategies.

A comparison script generates descriptive statistics and validation plots for the resulting representative sets.

### Stage 3 — Background sampling
`scripts/03_background_sampling/` generates a **10× non-AMP background** for each clustering condition.
- peptide length distribution,
- per-length K+R composition.

Post-generation validation is performed using two-sample Kolmogorov–Smirnov tests and diagnostic plots.

### Stage 4 — Motif analysis 
Scripts in ``scripts/04_motif_analysis/`` perform motif discovery and scanning.

#### **Motif discovery**
``01_run_meme_streme_pipeline.py`` runs:
- **STREME** (discriminative motif discovery),
- **MEME** (ZOOPS model, discriminative setting),
using the matched non-AMP background as the negative set.

#### **Motif scanning**
``02_run_fimo_scanning.py`` scans the discovered motif models with **FIMO** across:
- AMP representatives,
- non-AMP background representatives.

This covers all combinations of:
- 2 clustering strategies,
- 2 motif discovery tools,
- 2 sequence sets.

### Stage 5 — Statistics and motif-family analysis
Scripts in ``scripts/05_statistics/`` perform enrichment analysis, reporting, family-level robustness analysis, and logo generation.

#### **01 — FIMO enrichment**
``01_fimo_enrichment_fdr.py``
- counts motif-positive sequences,
- computes enrichment using Fisher’s exact test,
- applies Benjamini–Hochberg FDR correction,
- exports full and significant motif tables.

#### **02 — FIMO reporting and robustness**
``02_fimo_reporting_and_robustness.py``
- cleans and standardises enrichment outputs,
- classifies motifs by amino acid composition,
- compares significant motifs across pipelines,
- produces summary plots and overlap analyses.

#### **03 — Motif family analysis**
``03_motif_family_analysis.py``
- groups significant motifs into families using sequence-similarity criteria,
- identifies families reproduced across multiple pipelines,
- summarises robust motif families.

#### **04 — Final motif reporting**
``04_final_motif_reporting.py``
- creates final ranked tables and summary figures,
- exports the main robust-family deliverables.

#### **05 — Motif logo generation**
``01_generate_family_motif_logos.py``
- extracts motif probability matrices from MEME/STREME outputs,
- generates Logomaker sequence logos for motifs belonging to robust families,
- exports family-organised logo panels and summary tables.

Sequence logos are generated for all motifs belonging to robust motif families using Logomaker.
Each family receives its own directory containing motif logos and summary tables.

---

## Execution order

```bash
# Stage 1 — AMP database
python scripts/01_dataset_construction/amp/01_standardize_raw_databases.py
python scripts/01_dataset_construction/amp/02_filter_individual_databases.py
python scripts/01_dataset_construction/amp/03_apply_global_mic_filter.py
python scripts/01_dataset_construction/amp/04_create_master_dataset.py
python scripts/01_dataset_construction/amp/05_clean_master_dataset.py

# Stage 1 — Non-AMP background
python scripts/01_dataset_construction/non_amp/01_dataset_construction_nonamp.py

# Stage 2 — Redundancy reduction
python scripts/02_redundancy_reduction/01_cluster_amp_80_identity.py
python scripts/02_redundancy_reduction/02_cluster_nonamp_80_identity.py
python scripts/02_redundancy_reduction/03_compare_cluster_representatives.py

# Stage 3 — Background sampling
python scripts/03_background_sampling/01_generate_distribution_matched_nonamp.py

# Stage 4 — Motif analysis
python scripts/04_motif_analysis/01_run_meme_streme_pipeline.py
python scripts/04_motif_analysis/02_run_fimo_scanning.py

# Stage 5 — Statistics
python scripts/05_statistics/01_fimo_enrichment_fdr.py
python scripts/05_statistics/02_fimo_reporting_and_robustness.py
python scripts/05_statistics/03_motif_family_analysis.py
python scripts/05_statistics/04_final_motif_reporting.py
python scripts/06_generate_logos/01_generate_family_motif_logos.py
```

All scripts resolve paths relative to the repository root and can therefore be run from any working directory.

---

## Installation
Create the conda environment
```bash
conda env create -f environment.yml
conda activate amp_motif
```
Core Python dependencies include:
- pandas
- numpy
- scipy
- statsmodels
- matplotlib
- openpyxl
- pyarrow
- logomaker

External tools required on `PATH`:
- CD-HIT
- MMseqs2
- MEME Suite (meme, streme, fimo)

---

## Data sources

The AMP master dataset is built from four public resources:
- CAMP
- DBAASP
- DRAMP
- dbAMP3

The non-AMP background is derived from manually downloaded reviewed UniProt entries.

Raw source files are excluded from version control because of size and/or licensing restrictions. Intermediate parquet files are tracked where possible to facilitate partial reproducibility without redistributing raw source data.

---

## Reproducibility scope
The full pipeline is reproducible from raw data when the source files are available locally.

From Stage 2 onward, tracked intermediate parquet files and lightweight result files allow partial reproduction without redistributing raw datasets or heavy temporary outputs.

---

## What is tracked in git
Tracked:
- all scripts in scripts/
- selected processed parquet files
- final curated AMP dataset
- lightweight statistical outputs (csv, png, parquet)
- MEME/STREME report files (meme.txt, streme.txt, .html)
- validation plots

Not tracked:
- raw source datasets
- generated non-AMP background FASTAs
- full FIMO scan outputs
- temporary MMseqs2 databases
- local Excel exports and verbose auxiliary files

See .gitignore for the exact tracking policy.

---

## Author

Mireia Rivas Bermúdez  
BSc Biotechnology — Universitat Autònoma de Barcelona