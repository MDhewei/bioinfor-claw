#!/usr/bin/env python3
"""
Single-cell RNA-seq basic analysis pipeline.
Pure NumPy/Pandas implementation without scanpy/anndata dependency.
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

try:
    from scipy import stats as _scipy_stats
    _HAVE_SCIPY = True
except ImportError:
    _scipy_stats = None
    _HAVE_SCIPY = False


# ---------------------------------------------------------------------------
# Pure-numpy fallbacks for scipy.stats
# ---------------------------------------------------------------------------
def _erf_approx(x):
    """Vectorised Abramowitz & Stegun erf approximation."""
    x = np.asarray(x, float)
    t = 1.0 / (1.0 + 0.3275911 * np.abs(x))
    poly = (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t
              - 0.284496736) * t + 0.254829592) * t
    y = 1.0 - poly * np.exp(-(x ** 2))
    return np.sign(x) * y


def _log_gamma(x):
    """log Γ(x) via Lanczos."""
    if x < 0.5:
        return np.log(np.pi / np.sin(np.pi * x)) - _log_gamma(1 - x)
    x -= 1
    a = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
         771.32342877765313, -176.61502916214059, 12.507343278686905,
         -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]
    t = x + 7.5
    return (0.5 * np.log(2 * np.pi) + (x + 0.5) * np.log(t) - t
            + np.log(a[0] + sum(a[i] / (x + i) for i in range(1, 9))))


def _betainc_cf(a, b, x):
    """Regularised incomplete beta I_x(a,b) via Lentz continued fraction."""
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    lbeta = _log_gamma(a) + _log_gamma(b) - _log_gamma(a + b)
    front = np.exp(a * np.log(x) + b * np.log(1 - x) - lbeta) / a
    f = 1.0
    C = 1.0
    D = 1.0 - (a + b) * x / (a + 1)
    if abs(D) < 1e-30:
        D = 1e-30
    D = 1.0 / D
    f = D
    for m in range(1, 100):
        for step in (1, 2):
            if step == 1:
                num = m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m))
            else:
                num = -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1))
            D = 1.0 + num * D
            if abs(D) < 1e-30:
                D = 1e-30
            D = 1.0 / D
            C = 1.0 + num / C
            if abs(C) < 1e-30:
                C = 1e-30
            delta = C * D
            f *= delta
            if abs(delta - 1.0) < 1e-10:
                break
    return min(front * f, 1.0)


def _t_sf(t_val, df):
    """P(T > t_val) for Student-t with df degrees of freedom."""
    if df <= 0:
        return 0.5
    x = float(df) / (float(df) + float(t_val) ** 2)
    return 0.5 * _betainc_cf(df / 2.0, 0.5, x)


def _gammainc_lower(a, x):
    """Regularised lower incomplete gamma P(a, x) = 1 - Q(a, x)."""
    x = np.asarray(x, float)
    scalar = x.ndim == 0
    x = np.atleast_1d(x)
    result = np.zeros_like(x)
    for i, xi in enumerate(x):
        if xi <= 0:
            continue
        term = float(xi) ** a * np.exp(-xi) / max(a, 1e-30)
        s = term
        for n_iter in range(1, 300):
            term *= xi / (a + n_iter)
            s += term
            if abs(term) < 1e-12 * abs(s):
                break
        result[i] = min(s * np.exp(-_log_gamma(a)), 1.0)
    return float(result[0]) if scalar else result


def _welch_ttest(a, b):
    """Two-sample Welch t-test; returns (t_stat, p_value)."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return np.nan, np.nan
    va = np.var(a, ddof=1)
    vb = np.var(b, ddof=1)
    se = np.sqrt(va / na + vb / nb)
    if se == 0:
        return 0.0, 1.0
    t = float((a.mean() - b.mean()) / se)
    df = (va / na + vb / nb) ** 2 / ((va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1))
    p = 2 * _t_sf(abs(t), df)
    return t, float(min(p, 1.0))


def _mannwhitneyu_fallback(a, b, alternative='two-sided'):
    """Mann-Whitney U; returns (U, p_value)."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    na, nb = len(a), len(b)
    if na == 0 or nb == 0:
        return np.nan, np.nan
    combined = np.concatenate([a, b])
    order = np.argsort(combined, kind='stable')
    ranks = np.empty(len(combined))
    ranks[order] = np.arange(1, len(combined) + 1)
    # handle ties with midranks
    sorted_c = combined[order]
    i = 0
    while i < len(sorted_c):
        j = i
        while j < len(sorted_c) - 1 and sorted_c[j + 1] == sorted_c[j]:
            j += 1
        if j > i:
            midrank = (ranks[order[i]] + ranks[order[j]]) / 2
            for k in range(i, j + 1):
                ranks[order[k]] = midrank
        i = j + 1
    U = np.sum(ranks[:na]) - na * (na + 1) / 2.0
    mu = na * nb / 2.0
    sigma = np.sqrt(na * nb * (na + nb + 1) / 12.0)
    z = (U - mu) / (sigma + 1e-15)
    p = float(2 * (1 - 0.5 * (1 + _erf_approx(abs(z) / np.sqrt(2)))))
    return float(U), min(p, 1.0)


class _FallbackStats:
    """Drop-in replacements for scipy.stats functions used in this script."""

    @staticmethod
    def ttest_ind(a, b, equal_var=True, nan_policy='propagate'):
        return _welch_ttest(a, b)

    @staticmethod
    def mannwhitneyu(a, b, alternative='two-sided'):
        return _mannwhitneyu_fallback(a, b, alternative)

    @staticmethod
    def ranksums(a, b):
        return _mannwhitneyu_fallback(a, b, 'two-sided')

    @staticmethod
    def spearmanr(a, b=None):
        if b is None:
            a, b = a[:, 0], a[:, 1]
        a = np.asarray(a, float)
        b = np.asarray(b, float)

        def _rank(arr):
            order = np.argsort(arr)
            r = np.empty(len(arr), float)
            r[order] = np.arange(1, len(arr) + 1)
            return r

        ra = _rank(a) - (len(a) + 1) / 2.0
        rb = _rank(b) - (len(b) + 1) / 2.0
        denom = np.sqrt((ra ** 2).sum() * (rb ** 2).sum()) + 1e-15
        r = float((ra * rb).sum() / denom)
        n = len(a)
        t = r * np.sqrt(max(n - 2, 1)) / np.sqrt(max(1 - r ** 2, 1e-15))
        p = 2 * _t_sf(abs(t), n - 2)
        return r, float(min(p, 1.0))

    @staticmethod
    def pearsonr(a, b):
        a = np.asarray(a, float) - np.mean(a)
        b = np.asarray(b, float) - np.mean(b)
        denom = np.sqrt((a ** 2).sum() * (b ** 2).sum()) + 1e-15
        r = float((a * b).sum() / denom)
        n = len(a)
        t = r * np.sqrt(max(n - 2, 1)) / np.sqrt(max(1 - r ** 2, 1e-15))
        p = 2 * _t_sf(abs(t), n - 2)
        return r, float(min(p, 1.0))

    class norm:
        @staticmethod
        def cdf(x):
            x = np.asarray(x, float)
            return 0.5 * (1.0 + _erf_approx(x / np.sqrt(2)))

    class chi2:
        @staticmethod
        def cdf(x, df):
            return _gammainc_lower(df / 2.0, np.asarray(x, float) / 2.0)


stats = _scipy_stats if _HAVE_SCIPY else _FallbackStats()

import warnings

warnings.filterwarnings('ignore')

# Set matplotlib backend
plt.switch_backend('Agg')


class SingleCellAnalysis:
    """Single-cell RNA-seq analysis pipeline."""

    def __init__(self, args):
        self.args = args
        self.X = None  # cells × genes matrix
        self.gene_names = None
        self.cell_barcodes = None
        self.cell_meta = None
        self.hvg_mask = None
        self.pca_matrix = None
        self.umap_coords = None
        self.clusters = None

    def log_msg(self, msg):
        """Print timestamped message."""
        print(f"[INFO] {msg}")

    def load_data(self):
        """Load count matrix and detect orientation."""
        self.log_msg(f"Loading count matrix from {self.args.input}")

        # Load with gene/cell names
        if self.args.input.endswith('.csv'):
            df = pd.read_csv(self.args.input, index_col=0)
        else:
            df = pd.read_csv(self.args.input, sep='\t', index_col=0)

        n_rows, n_cols = df.shape
        self.log_msg(f"Loaded matrix: {n_rows} rows × {n_cols} columns")

        # Auto-detect orientation: if cols >> rows, likely genes × cells
        if n_cols > n_rows * 2 and not self.args.transpose:
            self.log_msg("Auto-detected genes × cells orientation; transposing...")
            df = df.T
        elif self.args.transpose:
            self.log_msg("Transposing per user request...")
            df = df.T

        self.X = df.values.astype(np.float32)
        self.gene_names = df.columns.values
        self.cell_barcodes = df.index.values

        self.log_msg(f"Final orientation: {self.X.shape[0]} cells × {self.X.shape[1]} genes")

    def load_metadata(self):
        """Load optional cell metadata."""
        if self.args.metadata:
            self.log_msg(f"Loading metadata from {self.args.metadata}")
            self.cell_meta = pd.read_csv(self.args.metadata, sep='\t', index_col=0)
            # Ensure barcodes match
            common_barcodes = np.intersect1d(self.cell_barcodes, self.cell_meta.index)
            self.log_msg(f"Matched {len(common_barcodes)} cells between count and metadata")

    def qc_filtering(self):
        """Compute QC metrics and filter cells/genes."""
        self.log_msg("Computing QC metrics...")

        n_genes = (self.X > 0).sum(axis=1)  # genes per cell
        total_counts = self.X.sum(axis=1)   # total counts per cell

        # Mitochondrial percentage
        mito_genes = [g for g in self.gene_names if g.startswith(self.args.mito_prefix)]
        if len(mito_genes) > 0:
            mito_mask = np.isin(self.gene_names, mito_genes)
            pct_mito = (self.X[:, mito_mask].sum(axis=1) / total_counts * 100).astype(np.float32)
        else:
            pct_mito = np.zeros(self.X.shape[0])

        # Store original metrics before filtering
        self.qc_before = {
            'n_genes': n_genes,
            'total_counts': total_counts,
            'pct_mito': pct_mito
        }

        # Filter cells
        cell_mask = (
            (n_genes >= self.args.min_genes) &
            (n_genes <= self.args.max_genes) &
            (pct_mito <= self.args.max_mito_pct)
        )

        n_cells_before = self.X.shape[0]
        self.X = self.X[cell_mask, :]
        self.cell_barcodes = self.cell_barcodes[cell_mask]
        n_genes = n_genes[cell_mask]
        total_counts = total_counts[cell_mask]
        pct_mito = pct_mito[cell_mask]

        self.log_msg(f"Cells: {n_cells_before} → {self.X.shape[0]} (filtered {n_cells_before - self.X.shape[0]})")

        # Filter genes
        gene_mask = (self.X > 0).sum(axis=0) >= self.args.min_cells
        n_genes_before = len(self.gene_names)
        self.X = self.X[:, gene_mask]
        self.gene_names = self.gene_names[gene_mask]

        self.log_msg(f"Genes: {n_genes_before} → {self.X.shape[1]} (filtered {n_genes_before - self.X.shape[1]})")

        # Store filtered QC metrics
        self.qc_metrics = pd.DataFrame({
            'barcode': self.cell_barcodes,
            'n_genes': n_genes,
            'total_counts': total_counts.astype(int),
            'pct_mito': pct_mito
        })

    def normalize(self):
        """Library size normalization and log transform."""
        self.log_msg("Normalizing (CPM + log1p)...")

        total_counts = self.X.sum(axis=1, keepdims=True)
        self.X = self.X / total_counts * 10000  # CPM-like
        self.X = np.log1p(self.X)

    def highly_variable_genes(self):
        """Select highly variable genes."""
        self.log_msg(f"Selecting {self.args.n_highly_variable} highly variable genes...")

        mean_per_gene = self.X.mean(axis=0)
        var_per_gene = self.X.var(axis=0)

        # Dispersion
        dispersion = var_per_gene / (mean_per_gene + 1e-10)

        # Normalize by expression bin
        n_bins = 10
        gene_idx = np.argsort(mean_per_gene)
        bin_size = len(gene_idx) // n_bins
        norm_disp = np.zeros_like(dispersion)

        for i in range(n_bins):
            start_idx = i * bin_size
            end_idx = (i + 1) * bin_size if i < n_bins - 1 else len(gene_idx)
            bin_genes = gene_idx[start_idx:end_idx]
            bin_disp = dispersion[bin_genes]
            bin_mean = bin_disp.mean()
            bin_std = bin_disp.std() + 1e-10
            norm_disp[bin_genes] = (bin_disp - bin_mean) / bin_std

        # Select top HVGs
        hvg_idx = np.argsort(-norm_disp)[:self.args.n_highly_variable]
        self.hvg_mask = np.zeros(len(self.gene_names), dtype=bool)
        self.hvg_mask[hvg_idx] = True

        # Save HVG list
        hvg_df = pd.DataFrame({
            'gene': self.gene_names,
            'mean': mean_per_gene,
            'variance': var_per_gene,
            'dispersion': dispersion,
            'normalized_dispersion': norm_disp
        }).iloc[hvg_idx].reset_index(drop=True)

        self.hvg_list = hvg_df

        self.log_msg(f"Selected {self.hvg_mask.sum()} HVGs")

    def pca(self):
        """PCA on HVG matrix."""
        self.log_msg(f"Running PCA with {self.args.n_pcs} components...")

        X_hvg = self.X[:, self.hvg_mask]

        # Scale
        mean = X_hvg.mean(axis=0)
        std = X_hvg.std(axis=0) + 1e-10
        X_scaled = (X_hvg - mean) / std

        # SVD
        U, S, Vt = np.linalg.svd(X_scaled, full_matrices=False)

        # Keep top PCs
        n_pcs = min(self.args.n_pcs, len(S))
        self.pca_matrix = U[:, :n_pcs]
        self.pca_var_explained = (S[:n_pcs] ** 2) / (X_scaled.shape[0] - 1)

        self.log_msg(f"PCA: {self.pca_matrix.shape[0]} cells × {self.pca_matrix.shape[1]} PCs")

    def knn_graph(self, X, k):
        """Build k-NN graph using Euclidean distance (pure NumPy)."""
        n = X.shape[0]
        # Pairwise distances
        distances = np.zeros((n, n))
        for i in range(n):
            diffs = X - X[i:i+1]
            distances[i] = np.sqrt((diffs ** 2).sum(axis=1))

        # k-NN indices
        knn_idx = np.argsort(distances, axis=1)[:, 1:k+1]  # Skip self
        return knn_idx

    def umap(self):
        """Simplified UMAP embedding."""
        self.log_msg(f"Running UMAP with {self.args.n_neighbors} neighbors...")

        # Build k-NN graph
        knn_idx = self.knn_graph(self.pca_matrix, self.args.n_neighbors)

        # Initialize with scaled PCA
        embedding = self.pca_matrix[:, :2].copy()
        embedding = (embedding - embedding.min(axis=0)) / (embedding.max(axis=0) - embedding.min(axis=0) + 1e-10)
        embedding = embedding * 2 - 1

        # Gradient descent (simplified)
        n_epochs = 200
        learning_rate = 0.1

        for epoch in range(n_epochs):
            lr = learning_rate * (1 - epoch / n_epochs)
            grad = np.zeros_like(embedding)

            # Attractive forces (k-NN pairs)
            for i in range(embedding.shape[0]):
                for j in knn_idx[i]:
                    diff = embedding[i] - embedding[j]
                    dist = np.linalg.norm(diff) + 1e-10
                    grad[i] += diff / dist

            # Repulsive forces (random non-neighbors)
            for i in range(embedding.shape[0]):
                # Sample ~5 negative samples per point
                neg_samples = np.random.choice(embedding.shape[0], size=min(5, embedding.shape[0]), replace=False)
                for j in neg_samples:
                    if j == i or j in knn_idx[i]:
                        continue
                    diff = embedding[i] - embedding[j]
                    dist = np.linalg.norm(diff) + 1e-10
                    if dist < 1.0:
                        grad[i] -= diff / dist * 0.1

            # Update
            embedding -= lr * grad / (np.linalg.norm(grad, axis=1, keepdims=True) + 1e-10)

        self.umap_coords = embedding
        self.log_msg(f"UMAP: {self.umap_coords.shape}")

    def snn_graph(self, knn_idx):
        """Build shared nearest neighbor graph."""
        n = knn_idx.shape[0]
        snn_graph = np.zeros((n, n))

        for i in range(n):
            for j in knn_idx[i]:
                # Count shared neighbors
                shared = np.intersect1d(knn_idx[i], knn_idx[j])
                snn_graph[i, j] = len(shared)
                snn_graph[j, i] = len(shared)

        return snn_graph

    def clustering(self):
        """Greedy label propagation clustering."""
        self.log_msg(f"Clustering with resolution={self.args.resolution}...")

        # Build SNN graph
        knn_idx = self.knn_graph(self.pca_matrix, self.args.n_neighbors)
        snn_graph = self.snn_graph(knn_idx)

        # Modulate by resolution
        snn_graph = snn_graph ** self.args.resolution

        # Initialize clusters
        clusters = np.arange(self.X.shape[0])

        # Label propagation
        for iteration in range(100):
            new_clusters = clusters.copy()
            for i in range(self.X.shape[0]):
                neighbors = np.where(snn_graph[i] > 0)[0]
                if len(neighbors) > 0:
                    neighbor_labels = clusters[neighbors]
                    # Vote: most common neighbor label
                    unique, counts = np.unique(neighbor_labels, return_counts=True)
                    new_clusters[i] = unique[np.argmax(counts)]

            # Check convergence
            if np.all(new_clusters == clusters):
                break
            clusters = new_clusters

        # Relabel clusters 0, 1, 2, ...
        unique_clusters = np.unique(clusters)
        cluster_map = {old: new for new, old in enumerate(unique_clusters)}
        self.clusters = np.array([cluster_map[c] for c in clusters])

        self.log_msg(f"Found {len(np.unique(self.clusters))} clusters")

    def marker_genes(self):
        """Identify marker genes per cluster."""
        self.log_msg("Identifying marker genes...")

        unique_clusters = np.unique(self.clusters)
        marker_results = []

        for cluster_id in unique_clusters:
            cluster_mask = self.clusters == cluster_id
            n_cluster = cluster_mask.sum()

            for gene_idx in range(self.X.shape[1]):
                gene_expr = self.X[:, gene_idx]
                cluster_expr = gene_expr[cluster_mask]
                other_expr = gene_expr[~cluster_mask]

                mean_cluster = cluster_expr.mean()
                mean_other = other_expr.mean()

                if mean_cluster == 0 and mean_other == 0:
                    continue

                # Wilcoxon rank-sum test
                if len(other_expr) > 0:
                    stat, pval = stats.ranksums(cluster_expr, other_expr)
                else:
                    stat, pval = 0, 1.0

                log2fc = np.log2(mean_cluster / (mean_other + 1e-10) + 1e-10)

                marker_results.append({
                    'cluster': cluster_id,
                    'gene': self.gene_names[gene_idx],
                    'log2fc': log2fc,
                    'pvalue': pval,
                    'mean_cluster': mean_cluster,
                    'mean_other': mean_other
                })

        marker_df = pd.DataFrame(marker_results)

        # BH FDR correction
        marker_df['fdr'] = self._fdr_correct(marker_df['pvalue'].values)

        # Sort by score
        marker_df['score'] = -np.log10(marker_df['pvalue'] + 1e-300) * marker_df['log2fc']
        marker_df = marker_df.sort_values(['cluster', 'score'], ascending=[True, False])

        self.marker_genes_df = marker_df

    def _fdr_correct(self, pvals):
        """Benjamini-Hochberg FDR correction."""
        ranked_p = np.argsort(pvals)
        fdr = np.zeros_like(pvals)

        n = len(pvals)
        for i, idx in enumerate(ranked_p):
            fdr[idx] = pvals[idx] * n / (i + 1)

        return np.minimum(fdr, 1.0)

    def plot_qc_violin(self):
        """Plot QC metrics before/after filtering."""
        self.log_msg("Plotting QC metrics...")

        fig, axes = plt.subplots(1, 3, figsize=(14, 4))

        # Before filtering
        for ax, (key, data) in zip(axes, self.qc_before.items()):
            ax.hist(data, bins=50, alpha=0.5, label='Before', color='red')

        # After filtering
        axes[0].hist(self.qc_metrics['n_genes'], bins=50, alpha=0.5, label='After', color='blue')
        axes[1].hist(self.qc_metrics['total_counts'], bins=50, alpha=0.5, label='After', color='blue')
        axes[2].hist(self.qc_metrics['pct_mito'], bins=50, alpha=0.5, label='After', color='blue')

        axes[0].set_xlabel('n_genes')
        axes[1].set_xlabel('total_counts')
        axes[2].set_xlabel('pct_mito')

        for ax in axes:
            ax.set_ylabel('Cells')
            ax.legend()

        plt.tight_layout()
        outfile = os.path.join(self.args.outdir, 'qc_violin.png')
        plt.savefig(outfile, dpi=300, bbox_inches='tight')
        self.log_msg(f"Saved {outfile}")
        plt.close()

    def plot_umap_clusters(self):
        """Plot UMAP colored by clusters."""
        self.log_msg("Plotting UMAP clusters...")

        fig, ax = plt.subplots(figsize=(8, 8))

        # Color palette
        n_clusters = len(np.unique(self.clusters))
        colors = plt.cm.tab20(np.linspace(0, 1, n_clusters))

        for cluster_id in np.unique(self.clusters):
            mask = self.clusters == cluster_id
            ax.scatter(self.umap_coords[mask, 0], self.umap_coords[mask, 1],
                      c=[colors[int(cluster_id)]], label=f'C{int(cluster_id)}', s=20, alpha=0.7)

        ax.set_xlabel('UMAP1')
        ax.set_ylabel('UMAP2')
        ax.set_title('UMAP: Clusters')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)

        plt.tight_layout()
        outfile = os.path.join(self.args.outdir, 'umap_clusters.png')
        plt.savefig(outfile, dpi=300, bbox_inches='tight')
        self.log_msg(f"Saved {outfile}")
        plt.close()

    def plot_umap_metadata(self):
        """Plot UMAP colored by metadata."""
        if not self.args.color_by or self.cell_meta is None:
            return

        self.log_msg(f"Plotting UMAP colored by {self.args.color_by}...")

        if self.args.color_by not in self.cell_meta.columns:
            self.log_msg(f"Warning: {self.args.color_by} not in metadata")
            return

        # Match cells
        meta_col = self.cell_meta.loc[self.cell_barcodes, self.args.color_by]

        fig, ax = plt.subplots(figsize=(8, 8))

        unique_vals = meta_col.unique()
        colors = plt.cm.tab20(np.linspace(0, 1, len(unique_vals)))

        for i, val in enumerate(unique_vals):
            mask = meta_col == val
            ax.scatter(self.umap_coords[mask, 0], self.umap_coords[mask, 1],
                      c=[colors[i]], label=str(val), s=20, alpha=0.7)

        ax.set_xlabel('UMAP1')
        ax.set_ylabel('UMAP2')
        ax.set_title(f'UMAP: {self.args.color_by}')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)

        plt.tight_layout()
        outfile = os.path.join(self.args.outdir, 'umap_metadata.png')
        plt.savefig(outfile, dpi=300, bbox_inches='tight')
        self.log_msg(f"Saved {outfile}")
        plt.close()

    def plot_marker_heatmap(self):
        """Plot heatmap of top marker genes."""
        self.log_msg("Plotting marker heatmap...")

        unique_clusters = np.unique(self.clusters)

        # Get top 5 markers per cluster
        top_markers = []
        for cluster_id in unique_clusters:
            cluster_markers = self.marker_genes_df[self.marker_genes_df['cluster'] == cluster_id]
            top_genes = cluster_markers.head(5)['gene'].values
            top_markers.extend(top_genes)

        top_markers = np.unique(top_markers)
        top_markers = top_markers[:min(len(top_markers), 40)]

        # Compute mean per cluster
        heatmap_data = np.zeros((len(top_markers), len(unique_clusters)))

        for gene_idx, gene in enumerate(top_markers):
            gene_expression = self.X[:, self.gene_names == gene].ravel()
            for cluster_idx, cluster_id in enumerate(unique_clusters):
                mask = self.clusters == cluster_id
                heatmap_data[gene_idx, cluster_idx] = gene_expression[mask].mean()

        # Plot
        fig, ax = plt.subplots(figsize=(len(unique_clusters) + 2, len(top_markers) * 0.3))

        im = ax.imshow(heatmap_data, cmap='RdYlBu_r', aspect='auto')

        ax.set_xticks(range(len(unique_clusters)))
        ax.set_xticklabels([f'C{int(c)}' for c in unique_clusters])
        ax.set_yticks(range(len(top_markers)))
        ax.set_yticklabels(top_markers, fontsize=8)

        ax.set_xlabel('Cluster')
        ax.set_ylabel('Gene')
        ax.set_title('Top Marker Genes')

        plt.colorbar(im, ax=ax, label='Mean log1p expression')
        plt.tight_layout()

        outfile = os.path.join(self.args.outdir, 'marker_heatmap.png')
        plt.savefig(outfile, dpi=300, bbox_inches='tight')
        self.log_msg(f"Saved {outfile}")
        plt.close()

    def plot_umap_top_markers(self):
        """Plot top marker genes on UMAP grid."""
        self.log_msg("Plotting UMAP with top marker expressions...")

        unique_clusters = np.unique(self.clusters)
        n_clusters_to_plot = min(len(unique_clusters), 8)

        fig, axes = plt.subplots(2, 4, figsize=(16, 8))
        axes = axes.ravel()

        for plot_idx, cluster_id in enumerate(unique_clusters[:n_clusters_to_plot]):
            ax = axes[plot_idx]

            # Get top marker for this cluster
            cluster_markers = self.marker_genes_df[self.marker_genes_df['cluster'] == cluster_id]
            top_gene = cluster_markers.iloc[0]['gene']

            # Get expression
            gene_idx = np.where(self.gene_names == top_gene)[0]
            if len(gene_idx) == 0:
                continue

            gene_expr = self.X[:, gene_idx[0]]

            # Plot
            scatter = ax.scatter(self.umap_coords[:, 0], self.umap_coords[:, 1],
                               c=gene_expr, cmap='viridis', s=20, alpha=0.7)
            ax.set_title(f'C{int(cluster_id)}: {top_gene}', fontsize=10)
            ax.set_xlabel('UMAP1')
            ax.set_ylabel('UMAP2')
            plt.colorbar(scatter, ax=ax)

        # Hide unused subplots
        for idx in range(n_clusters_to_plot, len(axes)):
            axes[idx].set_visible(False)

        plt.tight_layout()
        outfile = os.path.join(self.args.outdir, 'umap_top_markers.png')
        plt.savefig(outfile, dpi=300, bbox_inches='tight')
        self.log_msg(f"Saved {outfile}")
        plt.close()

    def save_outputs(self):
        """Save TSV outputs."""
        self.log_msg("Saving output files...")

        # Cell metadata
        cell_meta = pd.DataFrame({
            'barcode': self.cell_barcodes,
            'n_genes': self.qc_metrics['n_genes'].values,
            'total_counts': self.qc_metrics['total_counts'].values,
            'pct_mito': self.qc_metrics['pct_mito'].values,
            'cluster_id': self.clusters,
            'UMAP1': self.umap_coords[:, 0],
            'UMAP2': self.umap_coords[:, 1]
        })

        outfile = os.path.join(self.args.outdir, 'cell_metadata.tsv')
        cell_meta.to_csv(outfile, sep='\t', index=False)
        self.log_msg(f"Saved {outfile}")

        # Marker genes
        marker_out = self.marker_genes_df[['cluster', 'gene', 'log2fc', 'pvalue', 'fdr', 'mean_cluster', 'mean_other']]
        outfile = os.path.join(self.args.outdir, 'marker_genes.tsv')
        marker_out.to_csv(outfile, sep='\t', index=False)
        self.log_msg(f"Saved {outfile}")

        # HVG list
        outfile = os.path.join(self.args.outdir, 'hvg_list.tsv')
        self.hvg_list.to_csv(outfile, sep='\t', index=False)
        self.log_msg(f"Saved {outfile}")

    def run(self):
        """Execute full pipeline."""
        self.log_msg("Starting single-cell analysis pipeline...")

        self.load_data()
        self.load_metadata()
        self.qc_filtering()
        self.normalize()
        self.highly_variable_genes()
        self.pca()
        self.umap()
        self.clustering()
        self.marker_genes()

        # Plots
        self.plot_qc_violin()
        self.plot_umap_clusters()
        self.plot_umap_metadata()
        self.plot_marker_heatmap()
        self.plot_umap_top_markers()

        # Outputs
        self.save_outputs()

        self.log_msg("Pipeline complete!")


def main():
    parser = argparse.ArgumentParser(
        description='Single-cell RNA-seq basic analysis pipeline'
    )
    parser.add_argument('--input', required=True, help='Count matrix (TSV/CSV)')
    parser.add_argument('--transpose', action='store_true', help='Transpose input')
    parser.add_argument('--metadata', help='Cell metadata TSV')
    parser.add_argument('--color-by', help='Metadata column to color UMAP')
    parser.add_argument('--min-genes', type=int, default=200, help='Min genes per cell')
    parser.add_argument('--max-genes', type=int, default=6000, help='Max genes per cell')
    parser.add_argument('--min-cells', type=int, default=3, help='Min cells per gene')
    parser.add_argument('--max-mito-pct', type=float, default=20.0, help='Max mito %')
    parser.add_argument('--mito-prefix', default='MT-', help='Mito gene prefix')
    parser.add_argument('--n-highly-variable', type=int, default=2000, help='Number of HVGs')
    parser.add_argument('--n-pcs', type=int, default=30, help='Number of PCs')
    parser.add_argument('--n-neighbors', type=int, default=15, help='UMAP neighbors')
    parser.add_argument('--resolution', type=float, default=0.5, help='Clustering resolution')
    parser.add_argument('--n-top-markers', type=int, default=10, help='Top markers per cluster')
    parser.add_argument('--outdir', required=True, help='Output directory')

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    analysis = SingleCellAnalysis(args)
    analysis.run()


if __name__ == '__main__':
    main()
