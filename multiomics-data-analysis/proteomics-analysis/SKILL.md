---
name: proteomics-analysis
description: Analyze mass spectrometry proteomics data (TMT, LFQ, or DIA-NN output). Normalizes protein intensities, filters low-quality measurements, identifies differentially expressed proteins, and generates publication-quality visualizations.
---

# Proteomics Data Analysis

## Purpose

Comprehensive analysis of quantitative proteomics data from **mass spectrometry** (TMT, Label-Free Quantification, or DIA). This skill:
- **Performs QC analysis** (sample correlation, missing value patterns, CV assessment)
- **Normalizes intensities** (median, quantile, log2-median, variance-stabilizing)
- **Imputes missing values** (minimum probability, k-NN approximate, or zero)
- **Applies batch correction** (mean-centering if batch effect present)
- **Identifies differentially expressed proteins** (Welch t-test with FDR correction)
- **Generates visualizations** (volcano plots, heatmaps, correlation matrices)

Suitable for comparing protein abundance across treatment groups, disease states, or developmental stages.

## Reuse policy

This skill is designed for:
- Quality control and normalization of proteomics datasets
- Detection of differentially expressed proteins between groups
- Batch effect correction when sample batches vary
- Publication-quality volcano plots and heatmaps
- Integration with downstream protein enrichment or pathway analysis

This skill requires **protein intensity matrix** (protein × sample) and **sample metadata** with group assignments.

## Inputs

### Required
- `--input` — protein intensity matrix (TSV file)
  - Rows = protein IDs or symbols (e.g., TP53, EGFR, MYC)
  - Columns = sample names
  - Values = protein abundances (intensities, LFQ values, TMT ratios, etc.)
  - Values can be raw or log2-transformed; script auto-detects

- `--metadata` — sample metadata table (TSV file)
  - Must contain columns: `sample_id`, and `group`
  - Can optionally include `batch` column for batch correction
  - `sample_id` must match column names in `--input`

### Optional
- `--mode` — analysis mode
  Choices: `qc`, `normalize`, `differential`, `all`
  Default: `all`

- `--quant-type` — quantification platform (for context)
  Choices: `tmt`, `lfq`, `dia`
  Default: `lfq`

- `--normalization` — normalization method
  Choices: `median`, `quantile`, `vsn_approx`, `log2_median`, `none`
  Default: `median`

- `--imputation` — missing value imputation strategy
  Choices: `minprob`, `knn_approx`, `zero`, `none`
  Default: `minprob`

- `--group-col` — column name in metadata for grouping
  Default: `group`

- `--ref-group` — reference/control group for fold-change direction
  Example: `control`, `WT`, `normal`

- `--fdr-cutoff` — FDR-adjusted p-value threshold
  Default: `0.05`

- `--fc-cutoff` — log2 fold change threshold
  Default: `1.0`

- `--min-valid-values` — fraction of valid (non-missing) values required per protein
  Default: `0.7` (at least 70% of samples must have a value)

- `--batch-col` — column name for batch (if present in metadata)
  Example: `batch`, `plate`, `date`
  If provided, batch correction is applied

- `--protein-col` — protein ID column name in input matrix
  Default: first column is used as protein identifier

- `--outdir` — output directory
  Default: `./proteomics_output`

## Input file formats

### Protein intensity matrix (--input)

```
Protein      sample1  sample2  sample3  sample4
TP53         245.5    312.1    1542.3   1891.2
EGFR         512.4    498.2    145.3    128.9
MYC          1234.5   1156.2   98.3     87.1
```

### Metadata (--metadata)

```
sample_id    group   batch
sample1      control batch1
sample2      control batch1
sample3      treat   batch2
sample4      treat   batch2
```

## Outputs

### QC mode (--mode qc or all)
- `qc_summary.txt` — per-sample QC metrics:
  - Number of proteins detected, CV, mean intensity

- `sample_correlation_heatmap.png` — Pearson correlation matrix of samples
- `missing_value_heatmap.png` — heatmap showing missing data patterns (red = missing)
- `cv_distribution.png` — histogram of coefficient of variation per sample

### Normalization mode (--mode normalize or all)
- `normalized_proteins.tsv` — normalized intensity matrix
- `normalization_summary.txt` — summary of normalization applied

### Differential expression (--mode differential or all)
- `de_results.tsv` — all tested proteins with columns:
  - `protein`, `log2FC`, `mean_ref`, `mean_treat`, `pvalue`, `fdr`, `significant`

- `significant_proteins.tsv` — subset passing `--fdr-cutoff` AND `--fc-cutoff` threshold

- `de_volcano_plot.png` — volcano plot (log2FC vs -log10 FDR)
- `de_heatmap_top50.png` — heatmap of top 50 significant proteins (log-normalized, Z-scored)

### General
- `analysis_summary.txt` — overall analysis summary
- All PNG plots at 300 DPI

## Execution policy

Construct the command from required and optional parameters.

### Command template

```
python scripts/proteomics_analysis.py \
  --input <PROTEIN_MATRIX> \
  --metadata <METADATA_FILE> \
  [--mode <MODE>] \
  [--quant-type <TYPE>] \
  [--normalization <METHOD>] \
  [--imputation <METHOD>] \
  [--group-col <COLUMN>] \
  [--ref-group <GROUP>] \
  [--fdr-cutoff <FLOAT>] \
  [--fc-cutoff <FLOAT>] \
  [--min-valid-values <FLOAT>] \
  [--batch-col <COLUMN>] \
  --outdir <OUTDIR>
```

### Example commands

#### Basic differential expression (all defaults)
```
python scripts/proteomics_analysis.py \
  --input protein_intensities.tsv \
  --metadata metadata.tsv \
  --ref-group control \
  --outdir results/
```

#### QC only
```
python scripts/proteomics_analysis.py \
  --input protein_lfq.tsv \
  --metadata meta.tsv \
  --ref-group WT \
  --mode qc \
  --outdir results/
```

#### TMT data with batch correction and imputation
```
python scripts/proteomics_analysis.py \
  --input tmt_ratios.tsv \
  --metadata metadata.tsv \
  --quant-type tmt \
  --normalization quantile \
  --imputation minprob \
  --batch-col batch \
  --ref-group UT \
  --mode all \
  --outdir results/
```

#### Stricter DE with LFQ data
```
python scripts/proteomics_analysis.py \
  --input lfq_values.tsv \
  --metadata samples.tsv \
  --quant-type lfq \
  --normalization log2_median \
  --imputation knn_approx \
  --ref-group normal \
  --fdr-cutoff 0.01 \
  --fc-cutoff 1.5 \
  --min-valid-values 0.8 \
  --outdir results/
```

## Parameter decision guide

| User intent | Parameter to set |
|---|---|
| "Check data quality first" | `--mode qc` |
| "Just normalize my data" | `--mode normalize` |
| "Find differentially expressed proteins" | `--mode differential` or `--mode all` |
| "Do complete analysis" | `--mode all` (default) |
| TMT data | `--quant-type tmt --normalization quantile` |
| Label-free (LFQ) | `--quant-type lfq --normalization log2_median` |
| DIA data | `--quant-type dia --imputation knn_approx` |
| Has batch effects | `--batch-col batch` (include in metadata) |
| Strict significance | `--fdr-cutoff 0.01 --fc-cutoff 1.5` |
| Lenient/discovery | `--fdr-cutoff 0.1 --fc-cutoff 0.5` |
| Lots of missing data | `--imputation knn_approx --min-valid-values 0.5` |
| Minimal missing data | `--imputation minprob --min-valid-values 0.8` |
| Skip imputation | `--imputation none --min-valid-values 0.9` |

## Failure conditions

Fail clearly if:
- `--input`, `--metadata`, or required files are missing
- Protein intensity matrix cannot be parsed or is empty
- Metadata cannot be parsed or is missing required columns
- Sample IDs in metadata do not match column names in intensity matrix
- `--ref-group` label is not found in metadata group column
- Less than 2 samples per group (t-test invalid)
- All proteins filtered out by `--min-valid-values` threshold
- Invalid `--normalization` or `--imputation` choices provided
- Intensity values cannot be converted to numeric format

## Agent trigger examples

**Trigger this skill when the user asks:**
- "Analyze my mass spectrometry proteomics data"
- "Find differentially expressed proteins between samples"
- "Normalize my TMT data"
- "QC check on my LFQ protein quantification"
- "Which proteins are upregulated in treated vs control?"
- "Make a volcano plot from my proteomics data"
- "Correct batch effects in my protein abundances"

**Do NOT trigger this skill when the user asks:**
- "Identify proteins from spectra" → use peptide identification skill (e.g., MaxQuant)
- "Align sequences" → use sequence alignment tools
- "RNA-seq analysis" → use rnaseq-differential-expression
- "Protein-protein interactions" → use network analysis skill
- "Pathway enrichment" → use separate enrichment skill

## Notes

- **Statistical testing**: Uses Welch's t-test (unequal variance assumption) with Benjamini-Hochberg FDR correction.
- **Log2 transformation**: Script auto-detects raw vs log-scale (if max > 100, assumes raw and applies log2).
- **Normalization methods**:
  - `median`: subtract column median
  - `quantile`: rank-based quantile normalization
  - `log2_median`: log2-transform then subtract column median
  - `vsn_approx`: variance-stabilizing via asinh transform
  - `none`: no normalization (useful if pre-normalized)
- **Imputation methods**:
  - `minprob`: fill with low random values (mean - 1.8 SD)
  - `knn_approx`: fill with k=5 nearest neighbor mean
  - `zero`: replace with 0
  - `none`: leave missing as-is
- **Fold change direction**: positive = higher in treatment; negative = lower in treatment.
- **Recommended thresholds**: FDR < 0.05, log2FC > 1.0 (publication-ready); adjust for discovery.
- All plots saved at 300 DPI in PNG format.
- Batch correction (if `--batch-col` provided) uses simple mean-centering approach.
