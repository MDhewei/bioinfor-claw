---
name: drug-sensitivity-for-gene
description: Correlate gene expression or dependency scores with drug sensitivity data from PRISM (DepMap). Identifies drugs whose sensitivity correlates with the gene's expression or essentiality, and generates correlation plots.
---

# Drug Sensitivity for Gene

## Purpose
Identify drugs whose cellular sensitivity correlates with a gene's expression or essentiality score across a panel of cancer cell lines. This enables:
- Discovery of synthetic lethal interactions (high dependency genes)
- Identification of drugs effective in cell lines with high/low gene expression
- Biomarker-driven drug selection
- Mechanistic understanding of drug efficacy across genetic contexts

## Use When
- You want to find drugs whose efficacy correlates with a specific gene's expression
- You're screening for synthetic lethal interactions (genes whose knockout increases drug sensitivity)
- You need to identify patient stratification biomarkers for drug response
- You're prioritizing gene targets based on downstream drug opportunities
- You want to understand the mechanistic connection between a gene and drug response

## Do Not Use When
- You need to predict drug response in a specific patient (use pharmacogenomics databases instead)
- You're looking for drug-drug interactions
- You want cell-cycle phase-specific sensitivity analysis
- You need epigenetic or post-translational modification effects
- You lack internet access (requires DepMap/PRISM API or local data)

## Expected Inputs & Outputs

### Inputs
- **--gene** (required): Gene symbol (e.g., TP53, BRCA1, KRAS)
- **--omics-type**: Type of gene measurement (choices: expression, dependency, copy_number; default: expression)
- **--correlation-method**: Statistical correlation method (choices: pearson, spearman; default: spearman)
- **--fdr-cutoff**: False Discovery Rate threshold (default: 0.05)
- **--top-n**: Number of top correlated drugs to visualize (default: 20)
- **--min-cell-lines**: Minimum cell lines with data for both gene and drug (default: 50)
- **--depmap-dir**: Path to local DepMap data directory (optional)
- **--prism-file**: Path to PRISM drug sensitivity TSV/HDF5 (optional)
- **--expression-file**: Path to DepMap expression matrix (optional)
- **--outdir**: Output directory

### Outputs
- `drug_correlation.tsv`: Results table (drug_name, moa, target, r, pvalue, fdr, n_cell_lines)
- `scatter_plot_top4.png`: 2x2 grid of scatter plots for top 4 drugs (gene expression vs sensitivity)
- `waterfall_plot.png`: Ranked bar chart of top N drug correlations (red=positive, blue=negative)
- `correlation_summary.txt`: Summary statistics and parameters used

## Procedure

### Step 1: Obtain Data
If local files not provided:
1. Download DepMap expression data (OmicsExpressionProteinCodingGenesExpectedCountsLogTPMPlusOne)
2. Download PRISM repurposing drug sensitivity screen data
3. Cache data locally to avoid repeated downloads

### Step 2: Extract Gene Vector
1. Retrieve gene expression or dependency scores across cell lines
2. Filter for cell lines with complete data (no missing values)
3. Normalize expression if needed (already log-transformed in DepMap)

### Step 3: Compute Correlations
1. For each drug in PRISM:
   - Identify cell lines with data for both gene and drug
   - Require minimum of --min-cell-lines with both measurements
   - Compute Pearson or Spearman correlation
   - Calculate p-value
   - Store: r, p-value, n_cell_lines, drug metadata (MOA, target)

### Step 4: Multiple Testing Correction
1. Apply Benjamini-Hochberg FDR correction across all drug correlations
2. Filter results by --fdr-cutoff (default 0.05)
3. Sort by absolute correlation value

### Step 5: Visualization
1. **Scatter plots**: Top 4 drugs, gene expression (x) vs drug sensitivity (y)
   - Add regression line
   - Color by cell line lineage if available
   - Include R and p-value
2. **Waterfall plot**: Ranked bar chart of top N drugs
   - Red bars: positive correlation (high expression = high sensitivity)
   - Blue bars: negative correlation (high expression = low sensitivity)
   - Sorted by effect size

## Key Execution Patterns

### Basic usage (gene expression vs all drugs):
```bash
python scripts/drug_sensitivity_for_gene.py \
  --gene TP53 \
  --outdir results/tp53_drugs
```

### CRISPR dependency correlation:
```bash
python scripts/drug_sensitivity_for_gene.py \
  --gene BRCA1 \
  --omics-type dependency \
  --correlation-method spearman \
  --fdr-cutoff 0.01 \
  --top-n 30 \
  --outdir results/brca1_synthetic_lethal
```

### Copy number and local data:
```bash
python scripts/drug_sensitivity_for_gene.py \
  --gene MYC \
  --omics-type copy_number \
  --prism-file /data/prism_repurposing.csv \
  --expression-file /data/depmap_expression.csv \
  --min-cell-lines 100 \
  --outdir results/myc_cn_drugs
```

### High-stringency filtering:
```bash
python scripts/drug_sensitivity_for_gene.py \
  --gene EGFR \
  --fdr-cutoff 0.001 \
  --min-cell-lines 200 \
  --top-n 10 \
  --outdir results/egfr_high_confidence
```

## Parameter Decision Guide

| Scenario | --omics-type | --correlation-method | --fdr-cutoff | --min-cell-lines | --top-n |
|----------|--------------|----------------------|--------------|------------------|---------|
| Expression biomarker | expression | spearman | 0.05 | 50 | 20 |
| Synthetic lethal | dependency | pearson | 0.01 | 100 | 30 |
| Copy number driven | copy_number | pearson | 0.05 | 75 | 15 |
| Publication quality | expression | spearman | 0.001 | 150 | 10 |
| Discovery screen | expression | spearman | 0.1 | 30 | 50 |

## Failure Modes

| Error | Cause | Solution |
|-------|-------|----------|
| "Gene not found in expression matrix" | Gene symbol not in DepMap | Verify HGNC symbol; check case sensitivity |
| "No drugs with sufficient cell lines" | Most drug screens lack complete data | Lower --min-cell-lines to 30 or 20 |
| "Download timeout" | Network issue or DepMap API slow | Provide local --prism-file and --expression-file |
| "Memory error with large matrix" | Full DepMap matrix too large | Subset to specific genes or cancer types beforehand |
| "No significant correlations" | Gene not functionally related to any drug | Check gene's biological role; try lower --fdr-cutoff |
| "All correlations positive/negative" | Data quality issue or biased scaling | Verify gene vector normalization; check for batch effects |

## Technical Notes

- Correlation is computed only on cell lines with complete data for both gene and drug
- Spearman correlation recommended for non-normal distributions (typical for sensitivity data)
- FDR correction uses Benjamini-Hochberg procedure (controls false positive rate)
- Drug metadata (MOA, target) sourced from PRISM manifest
- Results are robust to cancer type differences (uses cross-cancer cell line panel)
- Positive correlation: high gene activity -> high drug sensitivity
- Negative correlation: high gene activity -> resistance to drug (protective/therapeutic antagonism)
- Minimum 50 cell lines recommended for reliable correlation estimates
- All plots use matplotlib; vector-compatible output (PNG 300 DPI)

