---
name: single-cell-basic-analysis
description: Basic single-cell RNA-seq analysis pipeline. Accepts a raw count matrix (cells × genes or genes × cells), performs QC filtering, normalization, highly variable gene selection, dimensionality reduction (PCA + UMAP), Leiden-like clustering, and marker gene identification. Generates UMAP plots and a marker gene heatmap.
---

## Purpose

Perform end-to-end basic single-cell RNA-seq analysis on raw count matrices. This skill automates the standard preprocessing, quality control, dimensionality reduction, clustering, and marker gene discovery workflow commonly used in single-cell genomics. All computations use pure NumPy/Pandas with no external C libraries, making it suitable for exploratory analysis and educational purposes.

## Use When

- **Exploratory scRNA-seq analysis**: You have a raw count matrix and want to quickly identify cell types and populations
- **Proof-of-concept studies**: Validating a new dataset before more detailed analysis
- **Educational/prototyping**: Learning single-cell analysis without complex dependencies (scanpy/AnnData)
- **Moderate dataset sizes**: Data fits in memory (typically <100k cells on standard compute nodes)
- **Standard preprocessing pipeline**: You need default-safe QC, normalization, and clustering steps

## Do NOT Use When

- **Large datasets (>500k cells)**: Performance will degrade; consider dask or on-disk approaches
- **Complex batch effects**: This pipeline does not include batch correction; pre-correct externally
- **Strict reproducibility required**: Stochastic algorithms (UMAP, clustering) have seeds but may differ from scanpy's implementations
- **Precise statistical testing**: Marker gene Wilcoxon test is simplified; use more rigorous methods if needed for publication
- **Multi-omics integration**: For CITE-seq, scATAC, etc., use specialized tools

## Expected Inputs

- **Count matrix** (required): TSV or CSV file with genes as rows and cells as columns (or vice versa)
  - Format: first row is cell barcodes, first column is gene names
  - Numeric counts (UMI or raw read counts)
  - Orientation auto-detected: if #columns >> #rows, assumed genes×cells and transposed

- **Cell metadata** (optional): TSV file with columns [cell_barcode, annotation_col1, annotation_col2, ...]
  - Used to color UMAP via `--color-by` parameter
  - Must have matching cell barcodes with count matrix

## Expected Outputs

All files written to `--outdir/`:

- **cell_metadata.tsv**: Per-cell QC metrics and results
  - Columns: barcode, n_genes, total_counts, pct_mito, cluster_id, UMAP1, UMAP2

- **marker_genes.tsv**: Ranked marker genes per cluster
  - Columns: cluster, gene, log2fc, pvalue, fdr, mean_cluster, mean_other

- **hvg_list.tsv**: Highly variable gene selection details
  - Columns: gene, mean, variance, dispersion, normalized_dispersion

- **qc_violin.png**: Violin plots of QC metrics (n_genes, total_counts, pct_mito)

- **umap_clusters.png**: UMAP projection colored by cluster assignment

- **umap_metadata.png**: UMAP colored by user-provided metadata (if --color-by given)

- **marker_heatmap.png**: Heatmap of top 5 marker genes × clusters

- **umap_top_markers.png**: Grid of UMAP plots showing expression of top markers per cluster

All images saved as 300 DPI PNG.

## Procedure

### Step 1: Data Loading and Orientation Detection
The script auto-detects count matrix orientation using a heuristic: if #columns >> #rows (e.g., 30000 genes, 5000 cells in genes×cells format), it transposes to cells×genes. If already cells×genes, no transpose occurs.

### Step 2: Quality Control Filtering
For each cell, compute:
- **n_genes**: count of genes with expression > 0
- **total_counts**: sum of all counts in that cell
- **pct_mito**: percentage of counts from mitochondrial genes (prefix MT- by default)

Filter criteria:
- Keep cells where: `min_genes ≤ n_genes ≤ max_genes AND pct_mito ≤ max_mito_pct`
- Keep genes expressed in ≥ `min_cells` cells
- Report before/after counts and plots

### Step 3: Normalization
1. **Library size normalization**: Divide each cell by its total counts, multiply by 10,000 (CPM-like)
2. **Log transformation**: Apply log1p(x) = log(x + 1)

### Step 4: Highly Variable Gene Selection
1. Compute per-gene mean and variance across cells
2. Compute dispersion = variance / (mean + epsilon)
3. Bin genes by mean expression (10 bins)
4. Normalize dispersion within each bin (z-score)
5. Select top `n_highly_variable` genes by normalized dispersion

### Step 5: Dimensionality Reduction (PCA)
1. Scale HVG submatrix: subtract mean, divide by std per gene
2. SVD decomposition via NumPy (no sklearn dependency)
3. Keep top `n_pcs` principal components
4. PCA matrix: cells × n_pcs

### Step 6: UMAP
1. Compute Euclidean k-NN graph (k = `n_neighbors`) in PCA space
2. Simplistic UMAP gradient descent (200 iterations):
   - Attractive forces between k-NN pairs
   - Repulsive forces between non-neighbors
   - Initialize with scaled PCA coordinates
3. Output 2D embedding for visualization

### Step 7: Clustering
1. Build shared nearest neighbor (SNN) graph from k-NN
2. Greedy label propagation community detection (~100 iterations)
3. Modulate edge weights by `--resolution`: higher resolution → more clusters
4. Assign each cell a cluster ID

### Step 8: Marker Gene Discovery
For each cluster vs. all other cells:
1. Compute per-gene mean expression in cluster and outside
2. Wilcoxon rank-sum test (non-parametric, pure NumPy)
3. Score = log2(mean_in / mean_out + 1) × -log10(pvalue)
4. Apply Benjamini-Hochberg FDR correction
5. Report top `n_top_markers` per cluster

### Step 9: Visualization and Output
Generate publication-quality PNG plots:
- QC violin plots (before and after filtering)
- UMAP colored by cluster
- UMAP colored by metadata (if provided)
- Marker heatmap and expression grids
- Write TSV tables with all results

## Key Execution Patterns

### Basic usage (minimal parameters):
```bash
python scripts/single_cell_basic_analysis.py \
  --input my_counts.tsv \
  --outdir ./results
```
Runs with default QC, clustering, and marker detection.

### With metadata coloring:
```bash
python scripts/single_cell_basic_analysis.py \
  --input my_counts.tsv \
  --metadata cell_types.tsv \
  --color-by cell_type \
  --outdir ./results
```

### Relaxed QC for low-quality samples:
```bash
python scripts/single_cell_basic_analysis.py \
  --input my_counts.tsv \
  --min-genes 100 \
  --max-genes 8000 \
  --max-mito-pct 30.0 \
  --outdir ./results
```

### More sensitive clustering:
```bash
python scripts/single_cell_basic_analysis.py \
  --input my_counts.tsv \
  --resolution 1.5 \
  --n-neighbors 30 \
  --outdir ./results
```

### Transpose input (genes × cells format):
```bash
python scripts/single_cell_basic_analysis.py \
  --input my_counts_genes_by_cells.tsv \
  --transpose \
  --outdir ./results
```

## Parameter Decision Guide

| Parameter | Typical Range | When to Adjust | Example Scenario |
|-----------|---------------|----------------|------------------|
| `--min-genes` | 100–500 | Lower if many small cells (neurons); raise if data is clean | Noisy PBMC: 200; Low-quality tissue: 100 |
| `--max-genes` | 5000–10000 | Raise if expecting large cells; lower to catch doublets | scRNA-seq: 6000; Sorted population: 8000 |
| `--max-mito-pct` | 10–30 | Lower for fresh cells; higher for old/stressed samples | Fresh: 15%; Cultured: 25% |
| `--min-cells` | 1–10 | Raise to remove rare genes and noise | Standard: 3; Conservative: 10 |
| `--n-highly-variable` | 1000–5000 | Reduce if memory-constrained; increase for rare cell types | Standard: 2000; Large data: 5000 |
| `--n-pcs` | 20–50 | Increase for complex data; decrease for simple | Multi-tissue: 50; Pure population: 20 |
| `--n-neighbors` | 10–30 | Lower for small datasets, higher for heterogeneous | Small dataset (n<5k): 10; Large: 30 |
| `--resolution` | 0.3–2.0 | Increase to discover finer subclusters | Coarse clustering: 0.3; Fine: 1.5 |
| `--n-top-markers` | 5–20 | More for comprehensive reporting | Quick QC: 5; Publication: 10 |

## Failure Modes & Troubleshooting

| Problem | Likely Cause | Solution |
|---------|--------------|----------|
| "Matrix is all zeros after filtering" | QC filters too strict | Relax `--min-genes`, `--max-mito-pct` |
| "UMAP has dense blob, no structure" | Too few HVGs or low-quality data | Increase `--n-highly-variable`, check input counts |
| "Clustering too coarse/fine" | Wrong resolution value | Adjust `--resolution` (0.3 for coarse, 1.5 for fine) |
| "Memory error on large data" | Data too big for RAM | Subset cells or genes, consider external tools |
| "Metadata column not found" | Column name mismatch | Verify column in metadata TSV matches --color-by |
| "Inconsistent cell barcodes between count and metadata" | Metadata barcodes don't match counts columns | Ensure exact barcode matches; check for whitespace |
| "Very few marker genes found" | Low biological signal or poor clustering | Check data quality; increase clustering resolution |

## Advanced Options

- **Stratification by batch**: Pre-process data separately per batch, then integrate externally
- **Custom gene filtering**: Provide a gene list file (future extension)
- **Protein markers**: Use metadata annotations to validate clusters post-hoc

## References

- Luecken & Theis (2019): "Current best practices in single-cell RNA-seq analysis"
- Wolf et al. (2018): "SCANPY, large-scale single-cell gene expression data analysis"
- McInnes & Healy (2018): "UMAP: Uniform Manifold Approximation and Projection"
