# plot-heatmap

---
name: plot-heatmap
description: Generate publication-quality heatmaps from numeric matrices with hierarchical clustering, custom scaling, and flexible annotations.
commands:
  - name: generate-heatmap
    description: Create a heatmap from a matrix file
    args:
      - --input (required): Path to matrix file (TSV/CSV)
      - --output (required): Output image path (PNG/SVG/PDF)
---

## Overview

**plot-heatmap** is a bioinformatics tool for generating publication-ready heatmaps from numeric data matrices. It supports:

- **Hierarchical clustering** with multiple linkage methods and distance metrics
- **Flexible scaling**: z-score or min-max normalization across rows, columns, or both
- **Advanced visualization**: dendrograms, cell annotations, custom colormaps, and missing-value handling
- **Input flexibility**: TSV/CSV formats with automatic or manual row label detection
- **Publication quality**: high-DPI output (PNG/SVG/PDF) with full layout control

Ideal for genomics (gene expression heatmaps), proteomics, metabolomics, and any tabular numeric data requiring visual clustering and pattern discovery.

---

## When to Use

### Use When:
- Analyzing expression matrices (RNA-seq, microarray, qPCR)
- Visualizing correlation or distance matrices
- Exploring high-dimensional tabular data with clustering
- Creating publication-quality figures for scientific manuscripts
- Comparing sample groups across features (genes, proteins, metabolites)
- Need reproducible, script-based heatmap generation with version control

### Do Not Use When:
- Data is non-numeric or contains many categorical features (use box plots, violin plots instead)
- Matrix has fewer than 3 rows or columns (use tables instead)
- You need interactive exploration (consider Plotly, Shiny, or JavaScript libraries)
- Data is not rectangular/tabular (use specialized network or 3D visualization tools)
- Real-time updates are required (this generates static images)

---

## Input/Output Specification

### Input Files

#### Matrix File (--input)
- **Format**: TSV (tab-separated) or CSV (comma-separated)
- **Structure**: Rows = features (genes, proteins, metabolites), Columns = samples
- **Data Type**: Numeric values (floats or integers)
- **Missing Values**: NaN, NA, empty cells (rendered with `--na-color`)
- **Row Labels**: First column (auto-detected) or specified via `--index-col`
- **Example**:
  ```
  GeneID    Sample1    Sample2    Sample3    Sample4
  GENE001   12.5       18.3       9.2        14.1
  GENE002   5.1        NA         6.3        7.8
  GENE003   22.4       25.1       20.9       23.5
  ```

#### Annotation Files (Optional)
- **Column Annotation File** (--col-annotation): TSV with sample names in first column, annotation groups in subsequent columns
  ```
  Sample    Condition    Batch
  Sample1   Control      Batch1
  Sample2   Treatment    Batch1
  Sample3   Control      Batch2
  ```
- **Row Annotation File** (--row-annotation): Similar structure with feature names in first column

### Output Files

#### Primary Output (--output)
- **Format**: PNG (default), SVG, or PDF
- **Resolution**: Configurable via `--dpi` (default 300 dpi for publication quality)
- **Contents**: Heatmap with optional dendrograms, annotations, colorbar, title

#### Optional Outputs
- **--output-svg**: Save additional SVG version (vector graphics, scalable)
- **--output-matrix**: Export processed (scaled and reordered) matrix as TSV for downstream analysis

---

## Execution Examples

### Example 1: Basic Usage (RNA-seq Expression)
```bash
python scripts/plot_heatmap.py \
  --input expression_matrix.tsv \
  --output heatmap_basic.png
```
- Plots top 100 rows by variance
- Row-wise z-score scaling
- Hierarchical clustering (Ward linkage, Euclidean distance)
- Auto-sized figure

### Example 2: Large Expression Matrix with Variance Filtering
```bash
python scripts/plot_heatmap.py \
  --input rnaseq_all_genes.tsv \
  --output top_genes_heatmap.png \
  --select-by-variance 50 \
  --title "Top 50 Differentially Expressed Genes" \
  --fig-width 8 --fig-height 10 \
  --dpi 300
```
- Selects top 50 genes by variance instead of 100
- Custom figure dimensions
- Publication-quality 300 DPI output

### Example 3: Correlation Matrix with Custom Coloring
```bash
python scripts/plot_heatmap.py \
  --input correlation_matrix.tsv \
  --output corr_heatmap.png \
  --scale none \
  --cmap coolwarm \
  --vmin -1 --vmax 1 \
  --center 0 \
  --title "Pearson Correlation Matrix" \
  --no-cluster-rows --no-cluster-cols
```
- No scaling (raw correlation values)
- Diverging colormap centered at 0
- No clustering (preserve matrix order)

### Example 4: Heatmap with Row and Column Dendrograms
```bash
python scripts/plot_heatmap.py \
  --input proteomics_data.csv \
  --output proteomics_clustered.png \
  --scale both \
  --scale-method minmax \
  --linkage complete \
  --distance-metric correlation \
  --dendrogram-ratio 0.2 \
  --show-values \
  --value-fmt ".2f"
```
- Both row-wise and column-wise min-max scaling
- Complete linkage clustering with correlation distance
- Cell values annotated with 2 decimal places

### Example 5: With Sample Annotations and Filtered Features
```bash
python scripts/plot_heatmap.py \
  --input metabolomics.tsv \
  --output metabolites_annotated.png \
  --min-variance 0.5 \
  --sample-cols "S1,S2,S3,S4,S5,S6" \
  --col-annotation sample_metadata.tsv \
  --cluster-rows --cluster-cols \
  --colorbar-label "Log2(Abundance)" \
  --output-matrix processed_metabolites.tsv
```
- Filters features with variance < 0.5
- Selects specific sample subset
- Adds sample condition annotation
- Exports reordered and scaled matrix

### Example 6: High-Resolution Multi-Panel Figure
```bash
python scripts/plot_heatmap.py \
  --input expression.tsv \
  --output heatmap_hires.png \
  --cell-width 0.15 --cell-height 0.15 \
  --dpi 600 \
  --cmap RdBu_r \
  --clip-zscore 2.5 \
  --font-family "Helvetica" \
  --base-fontsize 11 \
  --grid-color white --grid-linewidth 0.5
```
- Fixed cell dimensions for reproducible sizing
- Ultra-high resolution (600 DPI)
- Visible grid for clarity
- Professional font choice

---

## Parameter Decision Guide

| Goal | Key Parameters | Example |
|------|-----------------|---------|
| **Reduce dimensionality** | `--select-by-variance`, `--min-variance` | `--select-by-variance 50 --min-variance 0.1` |
| **Custom scaling** | `--scale` (row/col/both), `--scale-method` (zscore/minmax) | `--scale both --scale-method minmax` |
| **No clustering** | `--no-cluster-rows --no-cluster-cols` | Use for matrices where order is meaningful |
| **Custom distance metric** | `--distance-metric` (euclidean/correlation/cosine) | `--distance-metric correlation` for gene co-expression |
| **Diverging colormap** | `--cmap` (RdBu_r, RdYlBu, seismic, coolwarm), `--center` | `--cmap RdBu_r --center 0` for z-scores |
| **Hide labels** | `--no-row-labels --no-col-labels` | For large matrices (>100 rows) |
| **Add cell values** | `--show-values --value-fmt` | `--show-values --value-fmt ".2f"` |
| **Control figure size** | `--fig-width`, `--fig-height` OR `--cell-width`, `--cell-height` | `--cell-width 0.2 --cell-height 0.25` |
| **Annotation bars** | `--col-annotation`, `--row-annotation` | Paired TSV files with sample/feature metadata |
| **Publication quality** | `--dpi 300`, `--font-family Arial`, `--base-fontsize 10` | Default 300 DPI PNG + SVG backup |
| **Subset samples** | `--sample-cols` | `--sample-cols "WT1,WT2,WT3,Mutant1,Mutant2"` |
| **Save processed data** | `--output-matrix` | Export clustered/scaled matrix for downstream analysis |

---

## Technical Details

### Scaling Methods

#### Z-score Normalization (default)
- Centers each row (or column) by subtracting mean
- Divides by standard deviation
- Clips extreme values (default ±3 sigma) to prevent color saturation
- Suitable for comparing features with different absolute ranges
- Recommended for: gene expression, proteomics

#### Min-Max Scaling
- Scales values to [0, 1] range
- Formula: `(x - min) / (max - min)`
- Preserves relative differences within each row/column
- Recommended for: abundance data, abundance ratios

### Clustering

- **Method**: Hierarchical agglomerative clustering via SciPy
- **Linkage**: Ward (default), complete, average, single
- **Distance metrics**: Euclidean (default), correlation, cosine
- **Correlation distance**: 1 - Pearson correlation coefficient
- **Handling missing values**: Forward-fill mean values before clustering
- **Dendrogram**: Optional colored dendrograms for row/column clusters (set `--row-cluster-cutoff` or `--col-cluster-cutoff` to color by N clusters)

### Colormaps

Common choices:
- **Diverging** (centered): RdBu_r, RdYlBu, coolwarm, seismic (good for z-scores)
- **Sequential**: viridis, plasma, YlOrRd, Blues (good for abundance)
- **Qualitative**: Set1, tab10 (for categorical data)

All matplotlib colormaps supported.

### Missing Value Handling
- Cells with NaN/NA displayed as `--na-color` (default lightgray)
- Variance calculations skip NaN
- Clustering: NaN filled with row/column mean before distance computation
- Cell value annotation skips NaN cells

---

## Common Workflows

### 1. RNA-seq Differential Expression Heatmap
```bash
# Filter to top 50 DE genes, row-wise z-score
python scripts/plot_heatmap.py \
  --input deseq2_results_matrix.tsv \
  --output de_genes.png \
  --select-by-variance 50 \
  --scale row \
  --col-annotation sample_conditions.tsv \
  --title "Top 50 DE Genes"
```

### 2. Correlation Matrix Between Samples
```bash
# Pre-computed correlation; preserve order
python scripts/plot_heatmap.py \
  --input sample_correlations.tsv \
  --output correlation.png \
  --scale none \
  --no-cluster-rows --no-cluster-cols \
  --cmap coolwarm --vmin -1 --vmax 1 --center 0
```

### 3. Clustering-Focused Exploration
```bash
# Emphasize clustering with large dendrograms
python scripts/plot_heatmap.py \
  --input expression.tsv \
  --output clustered.png \
  --dendrogram-ratio 0.25 \
  --row-cluster-cutoff 5 \
  --col-cluster-cutoff 3
```

### 4. Reproducible Publication Figure
```bash
# Fixed dimensions, high DPI, minimal labels
python scripts/plot_heatmap.py \
  --input data.tsv \
  --output figure.png \
  --cell-width 0.15 --cell-height 0.12 \
  --dpi 600 \
  --row-label-size 8 --col-label-size 8 \
  --output-svg figure.svg
```

---

## Troubleshooting

| Issue | Cause | Solution |
|-------|-------|----------|
| Labels overlap | Too many rows/columns | Use `--no-row-labels`, `--row-label-size 6` |
| Wrong row order detected | TSV has unnamed index column | Use `--index-col GeneID` (column name) |
| Colors look washed out | vmin/vmax too wide | Set `--vmin`, `--vmax` explicitly |
| Clustering disabled | scipy not installed | Install: `pip install scipy>=1.8` |
| All white heatmap | All NaN values | Check data format, try sample file first |
| Memory error on large files | >50k rows | Use `--select-by-variance 1000` to pre-filter |

---

## Dependencies

- **numpy** >= 1.22: Numerical computing
- **pandas** >= 1.4: Data I/O and manipulation
- **matplotlib** >= 3.5: Plotting and rendering
- **scipy** >= 1.8: Hierarchical clustering (optional; skipped if unavailable)

Install all with: `pip install -r requirements.txt`

---

## Version History

- **v1.0** (2025): Initial release with full clustering, scaling, and annotation support
