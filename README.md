# AMP Comparative Clustering and Motif Analysis

Comparative computational analysis of antimicrobial peptide (AMP) clustering strategies and motif enrichment patterns.

---

## Overview

This project investigates the impact of sequence clustering methodologies on downstream motif discovery and enrichment analysis in antimicrobial peptides (AMPs).

Two clustering approaches (CD-HIT and MMseqs2) are systematically compared to evaluate differences in redundancy reduction and motif detection sensitivity. Motif enrichment across functional AMP categories is assessed using statistical testing (Chi-square and Fisher’s exact test).

---

## Objectives

- Compare sequence clustering approaches (CD-HIT vs MMseqs2)
- Identify enriched sequence motifs using MEME Suite
- Perform statistical enrichment analysis (Chi-square and Fisher tests)
- Evaluate motif distribution across functional AMP categories

---

## Tools and Technologies

- Python (NumPy, Pandas, SciPy, Biopython)
- CD-HIT
- MMseqs2
- MEME Suite
- Conda (Bioconda environment)

---

## Project Structure
```
scripts/ # Python analysis scripts
data/ # Example input data
results/ # Selected output figures and summaries
docs/ # Workflow diagrams and supporting materials
```

---

## Reproducibility

To recreate the computational environment:
```
conda env create -f environment.yml
conda activate amp-analysis
```

---

## Author

Mireia Rivas Bermúdez  
BSc Biotechnology