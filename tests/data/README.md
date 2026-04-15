# Bioinformatics Skills Test Data

This directory contains comprehensive synthetic test datasets for all 49 bioinformatics skills.

## Dataset Summary

### 1. RNA-seq Data
- **rnaseq_counts.tsv**: 500 genes × 12 samples
  - Gene names: TP53, BRCA1, EGFR, MYC, KRAS, Gene_006..Gene_500
  - Samples: Control_1-6, Treatment_1-6
  - Counts: NegBinomial-distributed
  - First 50 genes: differentially expressed (2-4x upregulated in treatment)

- **rnaseq_metadata.tsv**: 12 samples
  - Columns: sample_id, group, batch
  - Groups: control (6), treatment (6)

### 2. ATAC-seq / ChIP-seq Peaks
- **atac_peaks.narrowPeak**: 2000 peaks
  - Format: narrowPeak (10 columns: chrom, start, end, name, score, strand, signalValue, pValue, qValue, peak)
  - Chromosomes: chr1-chr22, chrX
  - Peak widths: 200-2000 bp

- **atac_peaks2.narrowPeak**: 1800 peaks
  - ~60% overlap with atac_peaks.narrowPeak for differential analysis testing

### 3. DNA Methylation
- **methylation_beta.tsv**: 3000 CpGs × 12 samples
  - CpG IDs: cg00000001..cg00003000
  - Beta values: 0-1 range (Beta distribution)
  - First 200 CpGs: differentially methylated (treatment +0.2-0.4)
  - Columns: CpG_ID, chr, position, [12 samples]

- **methylation_metadata.tsv**: 12 samples
  - Same structure as rnaseq_metadata.tsv

### 4. Proteomics
- **proteomics_intensities.tsv**: 1500 proteins × 12 samples
  - Log2 intensities: 15-35 range
  - ~20% missing values per sample
  - First 100 proteins: differentially abundant (+1.5 log2FC in treatment)
  - ~5% of rows entirely missing (low-abundance proteins)

- **proteomics_metadata.tsv**: 12 samples
  - Same structure as rnaseq_metadata.tsv

### 5. Single-cell RNA-seq
- **scrna_counts.tsv**: 500 cells × 1000 genes
  - Cell barcodes: CELL_0001..CELL_0500
  - 20 MT genes (MT-CO1..MT-CO20)
  - 3 cell types with distinct marker genes:
    - TypeA: 167 cells, high expression in genes 21-70
    - TypeB: 167 cells, high expression in genes 100-150
    - TypeC: 166 cells, high expression in genes 200-250
  - ~5% cells with high MT percentage (>30%)

- **scrna_metadata.tsv**: 500 cells
  - Columns: cell_id, true_celltype, batch

### 6. CRISPR Screen
- **crispr_counts.tsv**: 4000 sgRNAs × 5 samples
  - sgRNA × gene mapping: ~1000 genes × 4 sgRNAs/gene + controls
  - Samples: Plasmid_lib, Control_1-2, Treatment_1-2
  - First 30 essential genes: depleted in treatment (0.1-0.3x)
  - 30 enriched oncogenes: enriched in treatment (2-5x)
  - ~30% CV between replicates

### 7. Gene Lists
- **gene_list_A.txt**: 150 cancer-related genes
- **gene_list_B.txt**: 120 genes (~70% overlap with A)
- **gene_list_C.txt**: 100 genes (~30% overlap with A, ~20% with B)
- **ranked_gene_list.tsv**: 500 genes ranked by log2FC
  - Columns: gene, log2FC, pvalue, padj
  - ~50 upregulated (log2FC>1, padj<0.05)
  - ~50 downregulated (log2FC<-1, padj<0.05)

### 8. Survival / Clinical Data
- **survival_data.tsv**: 200 patients × 8 columns
  - Columns: sample_id, time_days, event, age, stage, grade, treatment, gene_expression
  - Time: Weibull-distributed, right-censored at 1825 days
  - Events: ~60%
  - High expression & stage IV patients: shorter survival

### 9. ML Classification
- **ml_features.tsv**: 200 samples × 200 features
  - First 20 features: discriminative (Class_A higher expression)
  - Remaining: random noise

- **ml_labels.tsv**: 200 samples
  - 100 Class_A, 100 Class_B

### 10. Expression Matrix
- **expression_matrix.tsv**: 50 samples × 500 genes
  - 3 clusters with 50 marker genes each:
    - Cluster 1 (samples 1-17): genes 1-50
    - Cluster 2 (samples 18-34): genes 100-150
    - Cluster 3 (samples 35-50): genes 200-250

- **sample_metadata.tsv**: 50 samples
  - Columns: sample_id, tissue_type, batch, treatment

### 11. Plot Test Data
- **de_results.tsv**: 1000 genes
  - DE volcano plot data
  - ~80 upregulated, ~80 downregulated (padj<0.05)

- **heatmap_matrix.tsv**: 80 genes × 24 samples
  - 3 sample groups × 8 samples with distinct patterns

- **boxviolin_data.tsv**: 120 observations
  - 4 groups with different means/spreads

- **scatter_data.tsv**: 100 points
  - x and y correlated (r~0.7)

- **bar_data.tsv**: 10 genes × 3 conditions × 5 replicates
  - Expression by condition

- **survival_plot.tsv**: 180 patients
  - 3 groups: GroupA (best), GroupB (medium), GroupC (poor) survival

### 12. API Skills Parameters
- **api_test_genes.txt**: 10 well-known genes (BRCA1, TP53, EGFR, etc.)
- **api_test_gene_list.txt**: 50 cancer-related genes
- **crispr_library_genes.txt**: 20 genes for library design
- **pubmed_test_query.txt**: "CRISPR base editing cancer 2023"

## Data Characteristics

### Quality Control Features
- Single-cell data includes cells with varying mitochondrial gene content
- Proteomics data includes realistic missing values and low-abundance proteins
- CRISPR screen includes technical replicates with realistic variance

### Biological Realism
- RNA-seq: Negative binomial counts with known DE signal
- ATAC-seq: Overlapping peak sets for differential analysis
- Methylation: Beta-distributed values with differential methylation
- Survival: Weibull-distributed times with meaningful covariates
- CRISPR: Essential gene depletion and oncogene enrichment patterns

### Integration Features
- Consistent sample IDs across related datasets
- Metadata files for all experiments
- Named features (genes, CpGs, etc.) for real-world compatibility
- Multiple clustering structures for unsupervised analysis testing

## Running Tests

### Full Test Suite
Run all 49 skill tests (excluding API-requiring skills):
```bash
python tests/test_skills.py
```

### Test Single Skill
```bash
python tests/test_skills.py --skill plot-volcano
```

### Test Skill Set
```bash
python tests/test_skills.py --skill-set multiomics-data-analysis
```

### Include API-Requiring Skills
```bash
python tests/test_skills.py --run-api
```

### Verbose Output
```bash
python tests/test_skills.py --verbose
```

### Dry Run (preview tests)
```bash
python tests/test_skills.py --dry-run
```

### Using Bash Wrapper
```bash
./tests/run_all_tests.sh
./tests/run_all_tests.sh --skill-set bioinformatics-plot-generator
./tests/run_all_tests.sh --verbose
./tests/run_all_tests.sh --run-api
```

## Test Output

Tests are executed with the following command pattern:
```bash
python tests/test_skills.py
```

Results include:
- Progress counter: [N/49] Testing skill-name...
- Per-test status: [PASS], [FAIL], [SKIP], or [API]
- Per-test runtime in seconds
- Summary table with all results
- TSV report saved to: tests/results/test_report_YYYYMMDD_HHMMSS.tsv

### Expected Test Coverage

**49 skills organized by category:**
- 5 multiomics skills
- 6 CRISPR design & analysis skills
- 7 gene list analysis skills
- 7 gene-centered analysis skills
- 3 machine learning skills
- 5 bioinformatics plotting skills
- 5 protein structure analysis skills
- 4 paper search & digest skills
- 3 lab search & tracking skills
- 3 public dataset access skills

**Tests without API access:**
- 22 runnable tests (all data provided locally)
- 27 tests requiring external API access (marked as SKIP by default, enable with --run-api)

## Data Import Examples

All files are tab-separated text files that can be imported into R, Python, or other analysis tools:

```python
import pandas as pd
df = pd.read_csv("rnaseq_counts.tsv", sep="\t", index_col=0)
```

```r
df <- read.delim("rnaseq_counts.tsv", row.names=1)
```

---
Generated: 2026-04-09
Total files: 30+ test data files
Total skills: 49
