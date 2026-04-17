# AMP Comparative Clustering & Motif Analysis

Reproducible bioinformatics pipeline to benchmark two redundancy-reduction strategies (**CD-HIT** vs **MMseqs2**) for antimicrobial peptide (**AMP**) motif discovery, enrichment analysis, motif-family robustness and membrane interaction mechanisms (**FMAP, PPM3**)

The project evaluates whether the choice of clustering strategy systematically affects motif discovery outcomes, and explores how identified motifs relate to membrane interaction mechanisms.

---

## Research question

How do redundancy-reduction strategies (CD-HIT vs MMseqs2) impact motif discovery in antimicrobial peptides, and can the resulting motifs be linked to specific membrane interaction mechanisms?

---

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
9. links motifs to membrane interaction mechanisms using structural modelling and membrane insertion analysis (FMAP + PPM).


### Conceptual output
The pipeline produces a multi-layer mapping:

motif → enrichment → family → structure → membrane interaction → mechanism

This enables moving from purely statistical motif discovery to **functional interpretation of AMP mechanisms**.

---

## Main outputs

Key results are located in:

- `results/07_motif_families/` — motif family assignments and robustness (Tomtom-based)
- `results/08_motif_logos/` — motif logos and consensus representations
- `results/09_motif_mechanism_input/` — selected peptide subsets (FMAP / ColabFold input)
- `results/10_motif_mechanism_output/` — FMAP + PPM parsed outputs
- `results/11_mechanism_analysis/` — FMAP vs PPM comparison and motif–mechanism analysis
- `results/statistics/11_motif_activity_analysis/` — final comparison tables and motif-level analysis

---

## Repository structure

```text
.
├── data/
│   ├── raw/                        # Not tracked
│   ├── intermediate/                     
│   └── final/
│       ├── AMP_MASTER.fasta              # Final curated AMP FASTA
│       └── DB_MASTER_CLEAN.parquet       # Final curated AMP metadata table
│
├── scripts/
│   ├── 01_dataset_construction/
│   ├── 02_redundancy_reduction/
│   ├── 03_background_sampling/
│   ├── 04_motif_analysis/
│   ├── 05_statistics/
│   ├── 06_generate_logos/
│   └── 07_motif_mechanism/
│
├── results/
│   ├── 01_amp_db_build/
│   ├── 02_clustering/
│   ├── 04_motif_discovery/
│   ├── 05_motif_scanning/
│   ├── 07_motif_families/
│   ├── 08_motif_logos/
│   ├── 09_motif_mechanism_input/
│   ├── 10_motif_mechanism_output/
│   ├── 11_mechanism_analysis/
│
├── envs/│   
├── .gitignore
└── README.md
```

---

## Pipeline

### Stage 1 — Dataset construction

#### **AMP database.**
Scripts in `scripts/01_dataset_construction/amp/` integrate four public AMP databases (**CAMP**, **DBAASP**, **dbAMP3**, **DRAMP**) into a standardised master dataset.

Main processing steps:
- harmonisation of heterogeneous source formats and column names,
- sequence normalisation to uppercase,
- removal of missing or duplicated sequences,
- filtering by peptide length and natural amino acid alphabet,
- removal of synthetic/designed peptides by name pattern,
- MIC-based activity filtering,
- final deduplication with metadata merging across sources.

Final outputs:
- `data/final/DB_MASTER_CLEAN.parquet`
- `data/final/AMP_MASTER.fasta`

#### **Non-AMP background.** 
`scripts/01_dataset_construction/non_amp/01_dataset_construction_nonamp.py` filters a UniProt-derived non-AMP FASTA by:
- minimum peptide length,
- valid amino acid alphabet,
- sequence uniqueness.

This produces the cleaned non-AMP pool used for background sampling.


### Stage 2 — Redundancy reduction 
Scripts in `scripts/02_redundancy_reduction/` cluster AMP and non-AMP sequence sets independently at 80% sequence identity using:
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
Scripts in `scripts/04_motif_analysis/` perform motif discovery and scanning.

#### **Motif discovery**
`01_run_meme_streme_pipeline.py` runs:
- **STREME** (discriminative motif discovery),
- **MEME** (ZOOPS model, discriminative setting),
using the matched non-AMP background as the negative set.

This step must be run within the MEME environment (`amp_meme`).

#### **Motif scanning**
`02_run_fimo_scanning.py` scans the discovered motif models with **FIMO** across:
- AMP representatives,
- non-AMP background representatives.

This covers all combinations of:
- 2 clustering strategies,
- 2 motif discovery tools,
- 2 sequence sets.


### Stage 5 — Statistics and motif-family analysis
Scripts in `scripts/05_statistics/` perform enrichment analysis, reporting, and family-level robustness analysis.

#### **01 — FIMO enrichment**
`01_fimo_enrichment_fdr.py`
- counts motif-positive sequences,
- computes enrichment using Fisher’s exact test,
- applies Benjamini–Hochberg FDR correction,
- exports full and significant motif tables.

#### **02 — FIMO reporting and robustness**
`02_fimo_reporting_and_robustness.py`
- cleans and standardises enrichment outputs,
- classifies motifs by amino acid composition,
- compares significant motifs across pipelines,
- produces summary plots and overlap analyses.

#### **03 — Motif family analysis**
`03_motif_family_analysis.py`
- groups significant motifs into families using sequence-similarity criteria,
- identifies families reproduced across multiple pipelines,
- summarises robust motif families.

#### **04 — Final motif reporting**
`04_final_motif_reporting.py`
- creates final ranked tables and summary figures,
- exports the main robust-family deliverables.

### Stage 6 — Motif logos

`01_generate_family_motif_logos.py`
- extracts motif probability matrices from MEME/STREME outputs,
- generates Logomaker sequence logos for motifs belonging to robust families,
- exports family-organised logo panels and summary tables.

Sequence logos are generated for all motifs belonging to robust motif families using Logomaker.
Each family receives its own directory containing motif logos and summary tables.


### Stage 7 — Motif–Mechanism analysis

Scripts in `scripts/07_motif_mechanism/` link motif-containing peptides to membrane interaction behaviour.

#### 7.1 Peptide selection

- motif-associated peptides are extracted from AMP dataset
- diversity filtering applied (sequence identity threshold)
- subsets generated for structural analysis

#### 7.2 Structure prediction
- Peptides associated with selected motifs are modelled using **ESMFold / ColabFold**
- FASTA inputs are generated

#### 7.3 Membrane interaction modelling
Two complementary approaches are used:

- **FMAP**
  - predicts membrane binding energetics
  - helix detection
  - estimates depth, and tilt

- **PPM**
  - calculates positioning of peptides in membranes
  - provides independent estimates of depth, tilt, and ΔG

#### 7.4 Comparative analysis
- FMAP and PPM outputs are harmonised
- activity classes are defined (interaction → insertion scale)
- agreement between models is evaluated

#### 7.5 Motif-level analysis 
- mechanisms are aggregated **per motif**
- frequency of mechanisms per motif is computed
- motif → mechanism relationships are identified

Outputs include:
- FMAP vs PPM comparison tables
- agreement metrics
- motif-level mechanism distributions
- high-confidence motif–mechanism associations

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

# Stage 6 — Motif logos
python scripts/06_generate_logos/01_generate_family_motif_logos.py

# Stage 7 — Mechanism analysis
python scripts/07_motif_mechanism/01_prepare_fmap_input.py
python scripts/07_motif_mechanism/02_parse_fmap_output.py
python scripts/07_motif_mechanism/03_parse_ppm_output.py
python scripts/07_motif_mechanism/04_fmap_analysis.py
python scripts/07_motif_mechanism/05_ppm_analysis.py
python scripts/07_motif_mechanism/06_compare_fmap_ppm.py
python scripts/07_motif_mechanism/07_motif_level_analysis.py

```
All scripts resolve paths relative to the repository root and can therefore be run from any working directory.

---

## Installation
### 1. Core Python environment (data processing and statistics)

```bash
conda env create -f envs/environment.yml
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

This environment is used for dataset construction, statistical analysis, and reporting.

### 2. Motif discovery environment (MEME Suite)
Motif discovery requires MEME Suite, which has specific Python and system dependencies.
A separate environment is provided for full reproducibility.

```bash
conda env create -f envs/environment_meme.yml
conda activate amp_meme
```
Check installation:

```bash
meme --version
streme --version
fimo --version
tomtom --version
```

### 3. Optional structural modelling environment (ColabFold)

Selected motif-associated peptides can be modelled with ColabFold in a separate optional environment.

```bash
conda env create -f envs/environment_colabfold.yml
conda activate colabfold
```
Check installation:

```bash
which colabfold_search
which colabfold_batch
colabfold_batch --help
```
This environment is optional and is only required for exploratory downstream peptide structure prediction.

### External tools required

The following tools must be installed and available in the corresponding environments:

- CD-HIT
- MMseqs2
- MEME Suite (`meme`, `streme`, `fimo`, `tomtom`)
- FMAP
- PPM 3.0
- ColabFold
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
The full pipeline is reproducible from raw data when the source files and external tools are available locally.

From Stage 2 onward, tracked intermediate parquet files and lightweight result files allow partial reproduction without redistributing raw datasets or heavy temporary outputs.

The motif–mechanism analysis (Stage 7) is reproducible given access to external tools such as FMAP, PPM, and ColabFold. These steps may require additional setup or web-based execution and should be considered semi-reproducible depending on the computational environment.

---

## What is tracked in git
Tracked:
- all scripts in scripts/
- selected processed parquet files
- final curated AMP dataset
- lightweight statistical outputs (csv, png, parquet)
- MEME/STREME report files (meme.txt, streme.txt, .html)
- validation plots
- final CSV/XLSX/PNG outputs from mechanism analysis
- motif-level comparison tables
- environment files in `envs/`

Not tracked:
- raw source datasets
- generated non-AMP background FASTAs
- full FIMO scan outputs
- temporary MMseqs2 databases
- local Excel exports and verbose auxiliary files
- FMAP raw outputs
- PPM raw outputs
- ColabFold structures

See .gitignore for the exact tracking policy.

---

## Author

Mireia Rivas Bermúdez  
BSc Biotechnology — Universitat Autònoma de Barcelona