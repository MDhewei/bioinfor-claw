---
name: dimensionality-reduction
description: Apply PCA, UMAP, or t-SNE dimensionality reduction to omics data matrices (gene expression, proteomics, methylation, etc.). Projects samples into 2D/3D space, colors by metadata, and generates publication-quality scatter plots with explained variance and loadings.
---

# Dimensionality Reduction for Omics Data

## Purpose

Reduce high-dimensional omics data (thousands of features) to 2D/3D visualizations while preserving biological structure and variance. Enables discovery of sample clusters, batch effects, and outliers. Supports PCA (linear, interpretable), UMAP (non-linear, local structure), and t-SNE (non-linear, local structure).

## Use When

- Visualizing sample relationships in RNA-seq, proteomics, metabolomics, or methylation data
- Exploring batch effects and confounding variables
- Validating clustering or cell-type predictions
- Identifying outlier samples
- Presenting high-dimensional data in publications
- Investigating how metadata annotations correlate with transcriptional space

## Do Not Use When

- Sample count < 3 (insufficient for meaningful projections)
- Feature count < 2 (nothing to reduce)
- Data is already well-characterized in 2D (e.g., flow cytometry, imaging)
- Inference on unseen test data is required (PCA only; t-SNE/UMAP do not generalize)

## Expected Inputs

### Required
- **--input**: TSV file with numeric expression matrix (rows = features, columns = samples, OR transpose). First column should be feature IDs (gene names, protein IDs, etc.). No missing values allowed.

### Optional
- **--metadata**: TSV file with sample annotations. First column = sample_id (matching input column names), additional columns for categorical/continuous variables.
- **--color-by**: Metadata column name to color scatter points (e.g., "cell_type", "batch", "treatment").
- **--shape-by**: Metadata column name for point shapes (categorical only).
- **--label-samples**: Flag to label each point with sample ID (useful for n_samples < 50).

## Expected Outputs

- **pca_projection.tsv**: PC1, PC2, PC3... coordinates for each sample
- **pca_variance.tsv**: Explained variance ratio and cumulative variance per PC
- **pca_scree_plot.png**: Elbow plot showing variance explained
- **pca_2d.png**: PC1 vs PC2 scatter plot, colored by metadata
- **pca_loadings.png**: Heatmap of top gene loadings for PC1/PC2
- **pca_3d.png**: Static 3D projection (if --plot-3d flag)
- **umap_projection.tsv**: UMAP embedding coordinates (if --method umap/all)
- **umap_2d.png**: UMAP scatter plot
- **tsne_projection.tsv**: t-SNE embedding coordinates (if --method tsne/all)
- **tsne_2d.png**: t-SNE scatter plot

## Procedure

1. **Load and validate data**: Read TSV, infer orientation (if n_rows >> n_cols, assume features × samples and transpose), check for missing values.
2. **Filter by variance** (if --filter-variance < 1.0): Keep only top X% most variable features.
3. **Normalize**: Apply scaling (standard = z-score, minmax = [0,1], none = no scaling).
4. **PCA**:
   - Center data
   - Compute SVD: X = U Σ V^T
   - Extract principal components (columns of U), explained variance (Σ^2 / n_samples)
   - Project onto first 50 PCs (or --n-components if specified)
   - Generate scree plot, 2D scatter, and loadings heatmap
5. **UMAP** (optional):
   - Reduce to 50 PCA dimensions first
   - Compute k-NN graph in PCA space (k = --n-neighbors)
   - Optimize 2D embedding via gradient descent (200 iterations)
   - Apply UMAP layout algorithm
6. **t-SNE** (optional):
   - Reduce to 50 PCA dimensions first
   - Compute pairwise affinities (Gaussian kernel)
   - Optimize 2D embedding via gradient descent (200 iterations, perplexity = --perplexity)
   - Apply Barnes-Hut approximation for efficiency
7. **Generate plots**: All at 300 DPI PNG, publication-ready formatting.

## Key Execution Patterns

```bash
# PCA only, colored by sample type
python dimensionality_reduction.py \
  --input expression_matrix.tsv \
  --metadata sample_metadata.tsv \
  --method pca \
  --color-by cell_type \
  --outdir results/

# All three methods, auto-transpose, top 5000 variable genes
python dimensionality_reduction.py \
  --input counts.tsv \
  --method all \
  --filter-variance 0.95 \
  --color-by treatment \
  --metadata metadata.tsv \
  --outdir results/

# UMAP with custom k-neighbors and min_dist
python dimensionality_reduction.py \
  --input data.tsv \
  --method umap \
  --n-neighbors 30 \
  --min-dist 0.05 \
  --color-by batch \
  --outdir results/

# t-SNE with custom perplexity for large dataset
python dimensionality_reduction.py \
  --input data.tsv \
  --method tsne \
  --perplexity 50 \
  --color-by cell_type \
  --label-samples \
  --outdir results/

# PCA with loadings analysis and 3D projection
python dimensionality_reduction.py \
  --input data.tsv \
  --method pca \
  --plot-loadings \
  --n-loadings 25 \
  --plot-3d \
  --color-by condition \
  --outdir results/
```

## Parameter Decision Guide

| Scenario | Method | Key Params | Rationale |
|----------|--------|-----------|-----------|
| Interpretability critical; explore variance structure | PCA | `--n-components 50`, `--plot-loadings` | Linear, explainable feature contributions |
| Non-linear structure, batch effects, clustering | UMAP | `--n-neighbors 15`, `--min-dist 0.1` | Preserves both local & global structure |
| Publication visualization; smooth embedding | t-SNE | `--perplexity 30` | Excellent for visual discovery (non-generalizable) |
| Exploring top variance genes | PCA | `--filter-variance 0.99` | Reduces noise, faster computation |
| Single small dataset (n_samples < 100) | t-SNE | Higher perplexity (40-50) | Handles small cohorts well |
| Large dataset (n_samples > 1000) | UMAP | `--n-neighbors 30` | Scalable, fast |
| Batch correction assessment | PCA/UMAP | Color by batch variable | Identifies systematic effects |
| Feature space exploration | PCA | `--plot-loadings`, `--n-loadings 30` | Understand gene driving variation |

## Failure Modes

| Error | Cause | Solution |
|-------|-------|----------|
| "Input matrix contains NaN values" | Missing data in expression matrix | Impute (e.g., zero-fill for counts) or filter features |
| "n_samples too low (< 3)" | Insufficient samples | Need minimum 3 samples for meaningful projections |
| "Metadata sample_id not found in input columns" | Mismatch between metadata and input | Check column names; ensure sample IDs align exactly |
| "UMAP crashed / diverged" | Too few neighbors for high-dimensional data | Reduce `--n-neighbors` or pre-filter by variance |
| "t-SNE produced uniform blob" | Perplexity too high for dataset size | Reduce `--perplexity` (rule: perplexity < n_samples / 3) |
| "Plots generated but colors are blank" | Metadata column is continuous/numeric | For scatter colors, use categorical columns; use `--color-by` numeric column for gradient |
| "PCA loadings plot unreadable" | Too many features; top genes overlap | Reduce `--n-loadings` or increase figure size |
| "Memory error on large matrix" | Dataset exceeds available RAM | Filter by variance (`--filter-variance 0.95`), or subset samples |

