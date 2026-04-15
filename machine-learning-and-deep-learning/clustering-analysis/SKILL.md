---
name: clustering-analysis
description: Cluster samples or features in omics data using K-means, hierarchical agglomerative clustering, DBSCAN, or consensus clustering. Evaluates cluster quality, assigns cluster labels, and generates heatmaps, dendrograms, and silhouette plots.
---

# Clustering Analysis for Omics Data

## Purpose

Partition samples (or features) into groups based on expression similarity. Enables unsupervised discovery of cell types, disease subtypes, co-expressed gene modules. Provides cluster quality metrics and consensus-based robustness assessment.

## Use When

- Discovering cell types or disease subtypes in single-cell or bulk RNA-seq
- Finding co-expressed gene modules
- Validating cell population predictions
- Generating sample labels for downstream analysis
- Assessing clustering robustness via consensus methods
- Determining optimal number of clusters

## Do Not Use When

- Sample count < 10 (unreliable clustering)
- Ground truth labels exist; use supervised classification instead
- Data is categorical (use appropriate clustering for categorical data)
- Real-time clustering needed (consensus clustering is slow)

## Expected Inputs

### Required
- **--input**: TSV numeric matrix (rows = features, columns = samples OR features to cluster). First column = feature/sample IDs. No missing values.

### Optional
- **--metadata**: Sample metadata TSV for heatmap annotations
- **--color-by**: Metadata column for heatmap sidebar colors
- **--consensus-n**: Number of consensus iterations (default: 100; increase for more robust assessment)

## Expected Outputs

- **cluster_labels.tsv**: sample_id, cluster_id, silhouette_score per sample
- **cluster_summary.tsv**: cluster_id, n_samples, centroid_distance_to_global_centroid
- **silhouette_plot.png**: Per-sample silhouette coefficient plot, sorted by cluster
- **cluster_heatmap.png**: Hierarchical clustered heatmap of expression, samples ordered by cluster
- **dendrogram.png**: Dendrogram for hierarchical method
- **consensus_matrix.png**: Heatmap of consensus co-clustering frequencies (if consensus method)
- **elbow_plot.png**: Silhouette scores for k in --k-range (if auto k-selection)

## Procedure

1. **Load and validate data**: Read TSV, check for missing values. Optionally transpose (if --transpose flag).
2. **Normalize**: Apply scaling (standard = z-score, minmax = [0,1], none = no scaling).
3. **Auto k-selection** (if --n-clusters 0):
   - Try each k in --k-range
   - Compute mean silhouette score for each k
   - Select k with highest silhouette score
   - Plot elbow/silhouette curve
4. **Clustering** (apply selected method(s)):
   - **K-means**: Lloyd's algorithm with k-means++ initialization, 100 iterations, 10 restarts, keep best inertia
   - **Hierarchical**: Build distance matrix, apply linkage (ward/complete/average/single), cut dendrogram at --n-clusters
   - **DBSCAN**: Compute pairwise distances, identify core points (eps-neighbors >= min-samples), expand clusters
   - **Consensus**: Run hierarchical clustering --consensus-n times on subsampled data (--consensus-subsample fraction), build consensus matrix, cluster consensus matrix
5. **Compute cluster quality metrics**:
   - Silhouette score per sample: s_i = (b_i - a_i) / max(a_i, b_i), where a_i = mean intra-cluster distance, b_i = min mean inter-cluster distance
   - Mean silhouette per cluster
   - Inertia (K-means only)
   - Cophenetic correlation (hierarchical only)
   - Consensus score (consensus clustering only)
6. **Generate visualizations**:
   - Silhouette plot: horizontal bars per sample, sorted by cluster
   - Heatmap: samples (columns) ordered by cluster assignment, rows = features, metadata sidebar
   - Dendrogram: hierarchical clustering tree
   - Elbow plot: silhouette vs k
   - Consensus matrix: heatmap of co-clustering frequencies

## Key Execution Patterns

```bash
# Hierarchical clustering, auto k-selection by silhouette
python clustering_analysis.py \
  --input expression_matrix.tsv \
  --method hierarchical \
  --n-clusters 0 \
  --k-range 2,10 \
  --outdir results/

# K-means with fixed k=5, colored heatmap by cell type
python clustering_analysis.py \
  --input expression_matrix.tsv \
  --method kmeans \
  --n-clusters 5 \
  --metadata sample_metadata.tsv \
  --color-by cell_type \
  --outdir results/

# Consensus clustering for robustness (100 iterations)
python clustering_analysis.py \
  --input expression_matrix.tsv \
  --method consensus \
  --n-clusters 5 \
  --consensus-n 100 \
  --consensus-subsample 0.8 \
  --outdir results/

# All methods for comparison
python clustering_analysis.py \
  --input expression_matrix.tsv \
  --method all \
  --n-clusters 4 \
  --linkage ward \
  --metadata metadata.tsv \
  --color-by treatment \
  --outdir results/

# Cluster genes (features) instead of samples
python clustering_analysis.py \
  --input expression_matrix.tsv \
  --method hierarchical \
  --transpose \
  --n-clusters 10 \
  --linkage complete \
  --outdir results/

# DBSCAN with custom epsilon
python clustering_analysis.py \
  --input expression_matrix.tsv \
  --method dbscan \
  --eps 0.5 \
  --min-samples 5 \
  --distance-metric euclidean \
  --outdir results/
```

## Parameter Decision Guide

| Scenario | Method | Key Params | Rationale |
|----------|--------|-----------|-----------|
| Know desired k; want interpretable results | K-means | `--n-clusters K`, `--distance-metric euclidean` | Fast, deterministic, interpretable centroids |
| Unknown k; want hierarchical exploration | Hierarchical | `--n-clusters 0`, `--linkage ward` | Dendrogram enables visual k selection |
| Variable cluster sizes; unknown k | DBSCAN | `--eps` auto-tuned, `--min-samples 5` | Finds density-based clusters of arbitrary shape/size |
| Want robustness assessment | Consensus | `--consensus-n 200`, `--consensus-subsample 0.8` | Evaluates stability across resamples |
| Gene co-expression modules | Hierarchical | `--transpose`, `--linkage average`, `--distance-metric correlation` | Correlation captures co-expression; average linkage less prone to chaining |
| Single-cell data (large n) | K-means | `--n-clusters` high (20-30), `--distance-metric euclidean` | Scalable to 100k cells |
| Small sample size (n < 50) | Hierarchical | `--linkage complete` | Complete linkage robust to outliers |
| Outlier-rich data | DBSCAN | `--eps`, `--min-samples` tuned | Labels outliers as noise (-1) |
| Exploratory; want all perspectives | all | All defaults; compare outputs | Identify consensus across methods |

## Failure Modes

| Error | Cause | Solution |
|-------|-------|----------|
| "K-means failed to converge" | Incompatible data scale or initialization | Increase max_iterations (code default: 100) or rescale data |
| "Silhouette score is negative" | Clusters overlap significantly | Increase --n-clusters or try different method |
| "Dendrogram is unreadable (many leaves)" | Too many samples for dendrogram | Use heatmap instead; subset samples for visualization |
| "DBSCAN found only 1 cluster + noise" | Epsilon too large | Decrease --eps or compute optimal eps via k-distance plot |
| "Consensus matrix is blank/uniform" | All subsamples produce identical clusters | --consensus-subsample fraction too high; reduce to 0.7 or lower |
| "Memory error on large matrix" | Dataset exceeds RAM | Pre-filter by variance or subset samples |
| "Hierarchical clustering tree is a star" | Linkage inappropriate for data (e.g., single linkage, chaining effect) | Switch to --linkage ward or complete |
| "Heatmap colors don't match metadata" | Metadata column mismatch | Verify --color-by column exists and matches sample IDs |
| "No silhouette plot generated" | DBSCAN produced noise points only | Tune --eps and --min-samples; try k-means instead |

