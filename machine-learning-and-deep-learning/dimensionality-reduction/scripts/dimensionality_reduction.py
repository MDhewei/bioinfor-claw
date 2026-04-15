#!/usr/bin/env python3
"""
Dimensionality Reduction for Omics Data
Supports PCA, UMAP, t-SNE with publication-quality visualizations
"""

import argparse
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm
try:
    from scipy.spatial.distance import pdist, squareform
    from scipy.spatial import cKDTree
    _HAVE_SCIPY = True
except ImportError:
    _HAVE_SCIPY = False

    def pdist(X, metric='euclidean'):
        X = np.asarray(X, float)
        if metric == 'euclidean':
            D_sq = np.sum(X ** 2, axis=1)
            D = np.sqrt(np.maximum(D_sq[:, None] + D_sq[None, :] - 2 * X @ X.T, 0))
        else:
            diff = X[:, None, :] - X[None, :, :]
            D = np.sqrt((diff ** 2).sum(-1))
        idx = np.triu_indices(len(X), k=1)
        return D[idx]

    def squareform(v):
        m = len(v)
        n = int(np.round((1 + np.sqrt(1 + 8 * m)) / 2))
        D = np.zeros((n, n))
        idx = np.triu_indices(n, k=1)
        D[idx] = v
        D += D.T
        return D

    class cKDTree:
        """Brute-force KD-tree replacement."""
        def __init__(self, data):
            self.data = np.asarray(data, float)

        def query(self, X, k=1):
            X = np.asarray(X, float)
            D_sq = (np.sum(self.data ** 2, axis=1)[None, :]
                    + np.sum(X ** 2, axis=1)[:, None]
                    - 2 * X @ self.data.T)
            D = np.sqrt(np.maximum(D_sq, 0))
            idx = np.argsort(D, axis=1)[:, :k]
            dists = np.take_along_axis(D, idx, axis=1)
            return dists, idx

import warnings
warnings.filterwarnings('ignore')

def load_and_validate_data(input_file):
    """Load TSV matrix and infer orientation"""
    df = pd.read_csv(input_file, sep='\t', index_col=0)

    # Check for NaN values
    if df.isna().any().any():
        raise ValueError("Input matrix contains NaN values. Please handle missing data first.")

    # Auto-detect orientation: if features >> samples, likely needs transpose
    n_rows, n_cols = df.shape
    if n_rows > n_cols and n_rows > 100:
        print(f"Detected features x samples ({n_rows} x {n_cols}), transposing to samples x features...")
        df = df.T

    return df

def load_metadata(metadata_file):
    """Load sample metadata"""
    meta = pd.read_csv(metadata_file, sep='\t', index_col=0)
    return meta

def filter_variance(X, feature_names, variance_threshold=1.0):
    """Keep top X% most variable features"""
    if variance_threshold >= 1.0:
        return X, feature_names

    variances = np.var(X, axis=0)
    n_keep = max(1, int(len(feature_names) * variance_threshold))
    top_idx = np.argsort(variances)[-n_keep:]
    top_idx = np.sort(top_idx)

    print(f"Keeping top {n_keep} most variable features out of {len(feature_names)}")
    return X[:, top_idx], feature_names[top_idx]

def scale_data(X, method='standard'):
    """Scale data: standard (z-score), minmax, or none"""
    if method == 'none':
        return X

    if method == 'standard':
        mean = np.mean(X, axis=0)
        std = np.std(X, axis=0)
        std[std == 0] = 1.0  # Avoid division by zero
        return (X - mean) / std

    elif method == 'minmax':
        min_val = np.min(X, axis=0)
        max_val = np.max(X, axis=0)
        range_val = max_val - min_val
        range_val[range_val == 0] = 1.0
        return (X - min_val) / range_val

    return X

def pca(X, n_components=50):
    """PCA using SVD"""
    # Center data
    X_centered = X - np.mean(X, axis=0)

    # SVD
    U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)

    # Explained variance
    var_explained = (S ** 2) / (X.shape[0] - 1)
    var_explained_ratio = var_explained / np.sum(var_explained)
    cumsum_var = np.cumsum(var_explained_ratio)

    # Project onto first n_components
    n_comp = min(n_components, len(S))
    components = Vt[:n_comp, :].T
    scores = U[:, :n_comp]

    return scores, components, var_explained_ratio[:n_comp], cumsum_var[:n_comp]

def umap_embed(X_pca, n_neighbors=15, min_dist=0.1, n_iters=200, random_seed=42):
    """Simplified UMAP embedding using k-NN and gradient descent"""
    np.random.seed(random_seed)
    n_samples = X_pca.shape[0]

    # Build k-NN graph
    tree = cKDTree(X_pca)
    distances, indices = tree.query(X_pca, k=n_neighbors + 1)
    distances = distances[:, 1:]  # Exclude self
    indices = indices[:, 1:]

    # Initialize 2D embedding randomly
    embedding = np.random.randn(n_samples, 2) * 0.01

    # Simple gradient descent optimization
    learning_rate = 1.0
    for iteration in range(n_iters):
        # For each neighbor pair, apply attractive force
        for i in range(n_samples):
            for j_idx in range(n_neighbors):
                j = indices[i, j_idx]

                # Compute distance in 2D
                delta = embedding[j] - embedding[i]
                dist_2d = np.linalg.norm(delta) + 1e-6

                # Attractive force
                grad = learning_rate * delta / (dist_2d + 1e-6)
                embedding[i] += grad * 0.1
                embedding[j] -= grad * 0.1

        # Repulsive forces with random negatives
        for i in range(min(n_samples, 100)):
            neg_samples = np.random.choice(n_samples, size=min(5, n_samples), replace=False)
            for j in neg_samples:
                if i != j:
                    delta = embedding[j] - embedding[i]
                    dist_2d = np.linalg.norm(delta) + 1e-6
                    grad = learning_rate * delta / (dist_2d ** 2 + 1e-6)
                    embedding[i] -= grad * 0.01

        learning_rate *= 0.95

    return embedding

def tsne_embed(X_pca, perplexity=30, n_iters=200, random_seed=42):
    """Simplified t-SNE using perplexity and gradient descent"""
    np.random.seed(random_seed)
    n_samples = X_pca.shape[0]

    # Compute pairwise distances
    distances = pdist(X_pca, metric='euclidean')
    D = squareform(distances)

    # Compute Gaussian affinities with perplexity
    P = np.zeros((n_samples, n_samples))
    for i in range(n_samples):
        # Binary search for sigma
        sigma = 1.0
        for _ in range(50):
            exp_d = np.exp(-D[i] / (2 * sigma ** 2))
            exp_d[i] = 0
            P_i = exp_d / np.sum(exp_d)
            entropy = -np.sum(P_i * np.log(P_i + 1e-10))
            perplexity_i = np.exp(entropy)

            if perplexity_i < perplexity:
                sigma *= 1.1
            else:
                sigma *= 0.9

        P[i] = exp_d / np.sum(exp_d)

    # Symmetrize
    P = (P + P.T) / 2
    P = np.maximum(P, 1e-12)

    # Initialize 2D embedding
    embedding = np.random.randn(n_samples, 2) * 0.01

    # Gradient descent with momentum
    learning_rate = 200.0
    momentum = 0.9
    velocity = np.zeros_like(embedding)

    for iteration in range(n_iters):
        # Compute low-dimensional affinities (student t-distribution)
        D_2d = pdist(embedding, metric='euclidean')
        D_2d = squareform(D_2d)

        Q = (1 + D_2d ** 2) ** (-1)
        np.fill_diagonal(Q, 0)
        Q = Q / np.sum(Q)
        Q = np.maximum(Q, 1e-12)

        # Compute gradient (symmetric KL divergence)
        PQ = P - Q
        grad = np.zeros_like(embedding)
        for i in range(n_samples):
            for j in range(n_samples):
                if i != j:
                    delta = embedding[i] - embedding[j]
                    grad[i] += 4 * PQ[i, j] * delta * (1 + D_2d[i, j] ** 2) ** (-1)

        # Update with momentum
        velocity = momentum * velocity - learning_rate * grad
        embedding += velocity

        if (iteration + 1) % 50 == 0:
            learning_rate *= 0.95

    return embedding

def plot_2d_scatter(scores, metadata=None, color_by=None, shape_by=None, label_samples=False,
                    title="2D Projection", xlabel="PC1", ylabel="PC2", outfile="projection.png"):
    """Generate 2D scatter plot"""
    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)

    sample_ids = None
    if metadata is not None and color_by is not None:
        color_col = metadata[color_by].astype(str)
        colors = plt.cm.tab20(np.linspace(0, 1, len(color_col.unique())))
        color_map = {val: colors[i] for i, val in enumerate(color_col.unique())}
        colors = np.array([color_map[val] for val in color_col])
        sample_ids = metadata.index
    else:
        colors = 'steelblue'

    ax.scatter(scores[:, 0], scores[:, 1], c=colors, alpha=0.7, s=100, edgecolors='black', linewidth=0.5)

    if label_samples and sample_ids is not None:
        for i, txt in enumerate(sample_ids):
            ax.annotate(txt, (scores[i, 0], scores[i, 1]), fontsize=8, alpha=0.7)

    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)

    # Legend
    if metadata is not None and color_by is not None:
        legend_elements = [mpatches.Patch(facecolor=color_map[val], label=val)
                          for val in sorted(color_col.unique())]
        ax.legend(handles=legend_elements, loc='best', fontsize=9)

    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    print(f"Saved {outfile}")
    plt.close()

def plot_scree(var_explained, outfile="scree_plot.png"):
    """Plot PCA scree plot"""
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)

    x = np.arange(1, len(var_explained) + 1)
    ax.plot(x, np.cumsum(var_explained), 'bo-', linewidth=2, markersize=8, label='Cumulative')
    ax.bar(x, var_explained, alpha=0.5, label='Per PC')

    ax.set_xlabel('Principal Component', fontsize=12)
    ax.set_ylabel('Explained Variance Ratio', fontsize=12)
    ax.set_title('PCA Scree Plot', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    print(f"Saved {outfile}")
    plt.close()

def plot_loadings(components, feature_names, n_loadings=20, outfile="loadings.png"):
    """Plot PCA loadings heatmap"""
    # Get top n_loadings features for PC1 and PC2
    pc1_idx = np.argsort(np.abs(components[:, 0]))[-n_loadings:]
    pc2_idx = np.argsort(np.abs(components[:, 1]))[-n_loadings:]
    top_idx = np.unique(np.concatenate([pc1_idx, pc2_idx]))
    top_idx = np.sort(top_idx)

    loadings_subset = components[top_idx, :2]
    feature_subset = feature_names[top_idx]

    fig, ax = plt.subplots(figsize=(10, max(6, len(top_idx) * 0.15)), dpi=300)

    im = ax.imshow(loadings_subset, cmap='RdBu_r', aspect='auto', vmin=-1, vmax=1)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['PC1', 'PC2'])
    ax.set_yticks(range(len(feature_subset)))
    ax.set_yticklabels(feature_subset, fontsize=8)
    ax.set_title('Top Feature Loadings', fontsize=14, fontweight='bold')

    plt.colorbar(im, ax=ax, label='Loading')
    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    print(f"Saved {outfile}")
    plt.close()

def plot_3d_projection(scores, metadata=None, color_by=None, outfile="projection_3d.png"):
    """Generate static 3D projection as PNG"""
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure(figsize=(10, 8), dpi=300)
    ax = fig.add_subplot(111, projection='3d')

    color_col = metadata[color_by].astype(str) if metadata is not None and color_by is not None else None

    if color_col is not None:
        colors = plt.cm.tab20(np.linspace(0, 1, len(color_col.unique())))
        color_map = {val: colors[i] for i, val in enumerate(color_col.unique())}
        colors = np.array([color_map[val] for val in color_col])
    else:
        colors = 'steelblue'

    ax.scatter(scores[:, 0], scores[:, 1], scores[:, 2], c=colors, alpha=0.7, s=100, edgecolors='black', linewidth=0.5)

    ax.set_xlabel('PC1', fontsize=10)
    ax.set_ylabel('PC2', fontsize=10)
    ax.set_zlabel('PC3', fontsize=10)
    ax.set_title('3D PCA Projection', fontsize=12, fontweight='bold')

    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    print(f"Saved {outfile}")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Dimensionality Reduction for Omics Data')
    parser.add_argument('--input', required=True, help='Input TSV matrix')
    parser.add_argument('--metadata', help='Sample metadata TSV')
    parser.add_argument('--method', choices=['pca', 'umap', 'tsne', 'all'], default='pca')
    parser.add_argument('--transpose', action='store_true', help='Transpose matrix')
    parser.add_argument('--n-components', type=int, default=50)
    parser.add_argument('--color-by', help='Metadata column for coloring')
    parser.add_argument('--shape-by', help='Metadata column for shapes')
    parser.add_argument('--label-samples', action='store_true')
    parser.add_argument('--scale', choices=['standard', 'minmax', 'none'], default='standard')
    parser.add_argument('--filter-variance', type=float, default=1.0)
    parser.add_argument('--n-neighbors', type=int, default=15)
    parser.add_argument('--min-dist', type=float, default=0.1)
    parser.add_argument('--perplexity', type=float, default=30)
    parser.add_argument('--random-seed', type=int, default=42)
    parser.add_argument('--outdir', default='results')
    parser.add_argument('--plot-loadings', action='store_true')
    parser.add_argument('--n-loadings', type=int, default=20)
    parser.add_argument('--plot-3d', action='store_true')

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Load data
    df = load_and_validate_data(args.input)
    if args.transpose:
        df = df.T

    X = df.values
    feature_names = np.array(df.columns)
    sample_ids = np.array(df.index)

    # Load metadata
    metadata = None
    if args.metadata:
        metadata = load_metadata(args.metadata)

    # Filter variance
    X, feature_names = filter_variance(X, feature_names, args.filter_variance)

    # Scale
    X = scale_data(X, args.scale)

    # PCA
    if args.method in ['pca', 'all']:
        print("Running PCA...")
        pca_scores, pca_components, var_exp, cumsum_var = pca(X, args.n_components)

        # Save projection
        pca_df = pd.DataFrame(pca_scores[:, :3], index=sample_ids,
                             columns=['PC1', 'PC2', 'PC3'])
        pca_df.to_csv(f'{args.outdir}/pca_projection.tsv', sep='\t')

        # Save variance
        var_df = pd.DataFrame({'pc': range(1, len(var_exp) + 1),
                               'explained_variance': var_exp,
                               'cumulative_variance': cumsum_var})
        var_df.to_csv(f'{args.outdir}/pca_variance.tsv', sep='\t', index=False)

        # Plots
        plot_scree(var_exp, f'{args.outdir}/pca_scree_plot.png')
        plot_2d_scatter(pca_scores, metadata, args.color_by, args.shape_by, args.label_samples,
                       'PCA Projection', 'PC1', 'PC2', f'{args.outdir}/pca_2d.png')

        if args.plot_loadings:
            plot_loadings(pca_components, feature_names, args.n_loadings, f'{args.outdir}/pca_loadings.png')

        if args.plot_3d and pca_scores.shape[1] >= 3:
            plot_3d_projection(pca_scores, metadata, args.color_by, f'{args.outdir}/pca_3d.png')

    # UMAP
    if args.method in ['umap', 'all']:
        print("Running UMAP...")
        pca_scores, _, _, _ = pca(X, 50)
        umap_scores = umap_embed(pca_scores, args.n_neighbors, args.min_dist, random_seed=args.random_seed)

        umap_df = pd.DataFrame(umap_scores, index=sample_ids, columns=['UMAP1', 'UMAP2'])
        umap_df.to_csv(f'{args.outdir}/umap_projection.tsv', sep='\t')

        plot_2d_scatter(umap_scores, metadata, args.color_by, args.shape_by, args.label_samples,
                       'UMAP Projection', 'UMAP1', 'UMAP2', f'{args.outdir}/umap_2d.png')

    # t-SNE
    if args.method in ['tsne', 'all']:
        print("Running t-SNE...")
        pca_scores, _, _, _ = pca(X, 50)
        tsne_scores = tsne_embed(pca_scores, args.perplexity, random_seed=args.random_seed)

        tsne_df = pd.DataFrame(tsne_scores, index=sample_ids, columns=['tSNE1', 'tSNE2'])
        tsne_df.to_csv(f'{args.outdir}/tsne_projection.tsv', sep='\t')

        plot_2d_scatter(tsne_scores, metadata, args.color_by, args.shape_by, args.label_samples,
                       't-SNE Projection', 'tSNE1', 'tSNE2', f'{args.outdir}/tsne_2d.png')

    print("Done!")

if __name__ == '__main__':
    main()
