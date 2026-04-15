# Quick Start Guide for Test Datasets

## Location
All test data files are in: `/sessions/festive-admiring-sagan/mnt/bioinfor-claw/tests/data/`

## Running the Generation Script

To regenerate all datasets from scratch:

```bash
cd /sessions/festive-admiring-sagan/mnt/bioinfor-claw/
python3 generate_test_data.py
```

This will recreate all 32 files in 1-2 minutes with consistent, reproducible data.

## Loading Data in Python

```python
import pandas as pd
import numpy as np

# RNA-seq
rnaseq_counts = pd.read_csv("tests/data/rnaseq_counts.tsv", sep="\t", index_col=0)
rnaseq_meta = pd.read_csv("tests/data/rnaseq_metadata.tsv", sep="\t")

# Single-cell RNA-seq
scrna_counts = pd.read_csv("tests/data/scrna_counts.tsv", sep="\t", index_col=0)
scrna_meta = pd.read_csv("tests/data/scrna_metadata.tsv", sep="\t")

# Methylation
methylation = pd.read_csv("tests/data/methylation_beta.tsv", sep="\t", index_col=0)

# Proteomics
proteomics = pd.read_csv("tests/data/proteomics_intensities.tsv", sep="\t", index_col=0)

# CRISPR screen
crispr = pd.read_csv("tests/data/crispr_counts.tsv", sep="\t", index_col=[0,1])

# Gene lists
gene_list_a = open("tests/data/gene_list_A.txt").read().strip().split("\n")

# Ranked genes
ranked = pd.read_csv("tests/data/ranked_gene_list.tsv", sep="\t")
```

## Loading Data in R

```r
# RNA-seq
rnaseq_counts <- read.delim("tests/data/rnaseq_counts.tsv", row.names=1)
rnaseq_meta <- read.delim("tests/data/rnaseq_metadata.tsv")

# Single-cell
scrna_counts <- read.delim("tests/data/scrna_counts.tsv", row.names=1)
scrna_meta <- read.delim("tests/data/scrna_metadata.tsv")

# Gene lists
gene_list_a <- readLines("tests/data/gene_list_A.txt")

# Ranked genes
ranked <- read.delim("tests/data/ranked_gene_list.tsv")
```

## Key Files by Analysis Type

### Differential Expression
- `rnaseq_counts.tsv` + `rnaseq_metadata.tsv`
- Known signal: genes 1-50 upregulated 2-4x in treatment

### Peak Analysis
- `atac_peaks.narrowPeak` (2000 peaks)
- `atac_peaks2.narrowPeak` (1800 peaks, ~60% overlap)

### DNA Methylation
- `methylation_beta.tsv` + `methylation_metadata.tsv`
- Known signal: first 200 CpGs shifted +0.2-0.4 in treatment

### Proteomics
- `proteomics_intensities.tsv` + `proteomics_metadata.tsv`
- Known signal: first 100 proteins +1.5 log2FC in treatment

### Single-cell Analysis
- `scrna_counts.tsv` + `scrna_metadata.tsv`
- 3 cell types with 50 marker genes each
- ~5% cells with high MT content for QC testing

### CRISPR Screen
- `crispr_counts.tsv`
- 30 essential genes depleted (0.1-0.3x)
- 30 oncogenes enriched (2-5x)

### Gene Enrichment
- `ranked_gene_list.tsv` + `gene_list_A/B/C.txt`
- ~50 upregulated, ~50 downregulated genes

### Survival Analysis
- `survival_data.tsv`
- 200 patients with time-to-event and covariates
- Hazards: Stage IV and high expression → shorter survival

### Machine Learning
- `ml_features.tsv` + `ml_labels.tsv`
- First 20 features discriminative
- 100 Class_A, 100 Class_B samples

### Clustering
- `expression_matrix.tsv` + `sample_metadata.tsv`
- 3 clusters with 50 marker genes each

### Visualization
- `de_results.tsv`: volcano plot data
- `heatmap_matrix.tsv`: heatmap data
- `boxviolin_data.tsv`: box/violin plot data
- `scatter_data.tsv`: scatter plot data
- `bar_data.tsv`: bar plot data
- `survival_plot.tsv`: survival curves data

## File Formats

All files are tab-separated (TSV) except:
- `atac_peaks.narrowPeak`: standard narrowPeak format
- `atac_peaks2.narrowPeak`: standard narrowPeak format
- `gene_list_*.txt`: one gene per line
- `api_test_*.txt`: gene lists
- `pubmed_test_query.txt`: plain text

## Data Characteristics

| Dataset | Rows | Columns | Size | Type |
|---------|------|---------|------|------|
| RNA-seq counts | 500 | 12 | 22 KB | NegBinomial |
| Methylation | 3000 | 12 | 761 KB | Beta |
| Proteomics | 1500 | 12 | 272 KB | Log2 |
| Single-cell | 500 | 1000 | 993 KB | Poisson |
| CRISPR | 4046 | 5 | 159 KB | NegBinomial |
| Survival | 200 | 8 | 9.2 KB | Weibull |
| ML Features | 200 | 200 | 771 KB | Normal |
| Expression | 50 | 500 | 482 KB | Normal |

## Sample IDs and Naming

### RNA-seq/Methylation/Proteomics samples:
- Control_1, Control_2, ..., Control_6
- Treatment_1, Treatment_2, ..., Treatment_6

### Single-cell barcodes:
- CELL_0001, CELL_0002, ..., CELL_0500

### CRISPR sgRNAs:
- sg_0001_0, sg_0001_1, ... (gene sgRNAs)
- sg_ess_00_0, ... (essential genes)
- NonTargeting_001, ... (controls)

### Gene names:
- TP53, BRCA1, EGFR, MYC, KRAS, ... (known genes)
- Gene_001, Gene_002, ... (numbered genes)
- MT-CO1, MT-CO2, ... (mitochondrial genes)

## Documentation

- **README.md**: Comprehensive dataset descriptions
- **MANIFEST.txt**: Complete file inventory
- **VALIDATION_REPORT.txt**: Data quality validation
- **QUICKSTART.md**: This file

## Tips for Testing

1. **Check signal detection**: First 50 RNA-seq genes should be significantly DE
2. **Test missing data handling**: Proteomics has ~20% missing values
3. **Verify clustering**: Expression matrix has 3 known clusters
4. **Test survival analysis**: Stage IV/high expression should have short survival
5. **Validate cell type identification**: Single-cell data has 3 known types
6. **Test CRISPR analysis**: Essential genes should be depleted, oncogenes enriched

---

Generated: 2026-04-09
All data reproducible with seed: numpy.random.seed(42)
