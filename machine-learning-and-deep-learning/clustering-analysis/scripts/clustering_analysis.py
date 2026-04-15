#!/usr/bin/env python3
"""
Clustering Analysis for Omics Data
Supports K-means, Hierarchical, DBSCAN, and Consensus clustering
"""

import argparse
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
try:
    from scipy.spatial.distance import pdist, squareform
    from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
    _HAVE_SCIPY = True
except ImportError:
    _HAVE_SCIPY = False

    def pdist(X, metric='euclidean'):
        """Pure-numpy condensed pairwise distance."""
        X = np.asarray(X, float)
        n = len(X)
        if metric == 'euclidean':
            diff = X[:, None, :] - X[None, :, :]
            D = np.sqrt((diff ** 2).sum(-1))
        elif metric == 'correlation':
            Xc = X - X.mean(1, keepdims=True)
            norms = np.linalg.norm(Xc, axis=1, keepdims=True) + 1e-10
            Xn = Xc / norms
            D = 1.0 - Xn @ Xn.T
            np.fill_diagonal(D, 0)
        elif metric == 'cosine':
            norms = np.linalg.norm(X, axis=1, keepdims=True) + 1e-10
            Xn = X / norms
            D = 1.0 - Xn @ Xn.T
            np.fill_diagonal(D, 0)
        else:
            diff = X[:, None, :] - X[None, :, :]
            D = np.sqrt((diff ** 2).sum(-1))
        idx = np.triu_indices(n, k=1)
        return D[idx]

    def squareform(condensed):
        """Convert condensed to square distance matrix."""
        m = len(condensed)
        n = int(np.round((1 + np.sqrt(1 + 8 * m)) / 2))
        D = np.zeros((n, n))
        idx = np.triu_indices(n, k=1)
        D[idx] = condensed
        D += D.T
        return D

    def linkage(distances, method='ward'):
        """Agglomerative hierarchical clustering; returns linkage matrix Z."""
        D = squareform(distances) if distances.ndim == 1 else np.asarray(distances, float)
        n = len(D)
        Z = []
        cluster_members = {i: [i] for i in range(n)}
        active = list(range(n))
        # Current distance matrix (modifiable)
        cur_D = D.copy()
        next_id = n

        for _ in range(n - 1):
            active_arr = np.array(active)
            # Find minimum off-diagonal distance
            min_d = np.inf
            mi, mj = -1, -1
            for ii in range(len(active_arr)):
                for jj in range(ii + 1, len(active_arr)):
                    i, j = active_arr[ii], active_arr[jj]
                    if cur_D[i, j] < min_d:
                        min_d = cur_D[i, j]
                        mi, mj = i, j

            if mi < 0:
                break

            ni = len(cluster_members[mi])
            nj = len(cluster_members[mj])
            Z.append([float(mi), float(mj), float(min_d), float(ni + nj)])

            # New cluster ID in expanded distance matrix
            new_id = next_id
            next_id += 1
            cluster_members[new_id] = cluster_members[mi] + cluster_members[mj]
            # Expand matrix
            new_row = np.zeros(len(cur_D) + 1)
            new_col_size = len(cur_D) + 1
            new_D = np.zeros((new_col_size, new_col_size))
            new_D[:len(cur_D), :len(cur_D)] = cur_D

            for k in active:
                dk_i = cur_D[k, mi] if max(k, mi) < len(cur_D) else np.inf
                dk_j = cur_D[k, mj] if max(k, mj) < len(cur_D) else np.inf
                if method == 'single':
                    d_new = min(dk_i, dk_j)
                elif method == 'complete':
                    d_new = max(dk_i, dk_j)
                elif method == 'average':
                    d_new = (ni * dk_i + nj * dk_j) / (ni + nj)
                else:  # ward
                    nk = len(cluster_members[k]) if k in cluster_members else 1
                    d_new = np.sqrt(max(((nk + ni) * dk_i ** 2 + (nk + nj) * dk_j ** 2
                                         - nk * min_d ** 2) / (nk + ni + nj), 0))
                new_D[len(cur_D), k] = d_new
                new_D[k, len(cur_D)] = d_new

            cur_D = new_D
            active.remove(mi)
            active.remove(mj)
            active.append(len(cur_D) - 1)  # index of new_id in cur_D

        return np.array(Z)

    def fcluster(Z, t, criterion='maxclust'):
        """Extract flat clusters from linkage matrix."""
        n = len(Z) + 1
        if criterion == 'maxclust':
            k = int(t)
            if k >= n:
                return np.arange(1, n + 1, dtype=int)
            heights = Z[:, 2]
            sorted_h = np.sort(heights)[::-1]
            cut = sorted_h[k - 1] if k <= len(sorted_h) else 0.0
            # Assign via union-find
            parent = list(range(2 * n))

            def find(x):
                while parent[x] != x:
                    parent[x] = parent[parent[x]]
                    x = parent[x]
                return x

            for row_idx, row in enumerate(Z):
                a, b, h = int(row[0]), int(row[1]), row[2]
                new_node = n + row_idx
                if h > cut:
                    parent[find(a)] = find(new_node)
                    parent[find(b)] = find(new_node)
                else:
                    parent[find(a)] = find(new_node)
                    parent[find(b)] = find(new_node)
            roots = {find(i) for i in range(n)}
            root_map = {r: idx + 1 for idx, r in enumerate(sorted(roots))}
            return np.array([root_map[find(i)] for i in range(n)], dtype=int)
        elif criterion == 'distance':
            cut = float(t)
            parent = list(range(2 * n))

            def find(x):
                while parent[x] != x:
                    parent[x] = parent[parent[x]]
                    x = parent[x]
                return x

            for row_idx, row in enumerate(Z):
                a, b, h = int(row[0]), int(row[1]), row[2]
                new_node = n + row_idx
                parent[find(a)] = find(new_node)
                parent[find(b)] = find(new_node)
                if h >= cut:
                    break
            roots = {find(i) for i in range(n)}
            root_map = {r: idx + 1 for idx, r in enumerate(sorted(roots))}
            return np.array([root_map[find(i)] for i in range(n)], dtype=int)
        else:
            return np.ones(n, dtype=int)

    def dendrogram(Z, ax=None, no_labels=False, **kwargs):
        """Stub dendrogram - draws a simplified tree if ax provided."""
        if ax is None:
            return {}
        n = len(Z) + 1
        ax.plot([0, n], [0, 0], 'k-', linewidth=0.5)
        ax.set_xlim(0, n)
        ax.set_title('Dendrogram (scipy unavailable - simplified)')
        return {'icoord': [], 'dcoord': [], 'leaves': list(range(n))}

import warnings
warnings.filterwarnings('ignore')

def load_and_validate_data(input_file, transpose=False):
    """Load TSV matrix"""
    df = pd.read_csv(input_file, sep='\t', index_col=0)

    if df.isna().any().any():
        raise ValueError("Input matrix contains NaN values.")

    if transpose:
        df = df.T

    return df

def load_metadata(metadata_file):
    """Load sample metadata"""
    meta = pd.read_csv(metadata_file, sep='\t', index_col=0)
    return meta

def scale_data(X, method='standard'):
    """Scale data"""
    if method == 'none':
        return X

    if method == 'standard':
        mean = np.mean(X, axis=0)
        std = np.std(X, axis=0)
        std[std == 0] = 1.0
        return (X - mean) / std

    elif method == 'minmax':
        min_val = np.min(X, axis=0)
        max_val = np.max(X, axis=0)
        range_val = max_val - min_val
        range_val[range_val == 0] = 1.0
        return (X - min_val) / range_val

    return X

def compute_pairwise_distances(X, metric='euclidean'):
    """Compute pairwise distance matrix"""
    if metric == 'euclidean':
        distances = pdist(X, metric='euclidean')
    elif metric == 'correlation':
        distances = pdist(X, metric='correlation')
    elif metric == 'cosine':
        distances = pdist(X, metric='cosine')
    else:
        distances = pdist(X, metric=metric)

    return squareform(distances)

def kmeans_plusplus_init(X, k):
    """K-means++ initialization"""
    np.random.seed(42)
    n_samples = X.shape[0]

    # Choose first centroid randomly
    centroids = [X[np.random.randint(n_samples)]]

    for _ in range(1, k):
        # Compute distances to nearest centroid
        distances = np.array([np.min([np.linalg.norm(x - c) for c in centroids]) for x in X])
        # Choose next centroid with probability proportional to distance squared
        probs = distances ** 2
        probs /= np.sum(probs)
        cumsum_probs = np.cumsum(probs)
        r = np.random.random()
        next_centroid_idx = np.searchsorted(cumsum_probs, r)
        centroids.append(X[next_centroid_idx])

    return np.array(centroids)

def kmeans(X, k, max_iter=100, n_restarts=10):
    """K-means clustering with k-means++ initialization"""
    best_inertia = np.inf
    best_labels = None
    best_centroids = None

    for restart in range(n_restarts):
        centroids = kmeans_plusplus_init(X, k)

        for iteration in range(max_iter):
            # Assign samples to nearest centroid
            distances = np.array([[np.linalg.norm(x - c) for c in centroids] for x in X])
            labels = np.argmin(distances, axis=1)

            # Update centroids
            new_centroids = np.array([X[labels == i].mean(axis=0) if np.sum(labels == i) > 0 else centroids[i]
                                     for i in range(k)])

            # Check convergence
            if np.allclose(centroids, new_centroids):
                break

            centroids = new_centroids

        # Compute inertia
        distances = np.array([[np.linalg.norm(x - c) for c in centroids] for x in X])
        labels = np.argmin(distances, axis=1)
        inertia = np.sum([np.min(distances[i]) ** 2 for i in range(len(X))])

        if inertia < best_inertia:
            best_inertia = inertia
            best_labels = labels
            best_centroids = centroids

    return best_labels, best_centroids

def hierarchical_clustering(X, linkage_method='ward', distance_metric='euclidean', n_clusters=None):
    """Hierarchical agglomerative clustering"""
    distances = pdist(X, metric=distance_metric)

    if linkage_method == 'ward' and distance_metric != 'euclidean':
        print(f"Ward linkage requires euclidean distance; using {distance_metric}")

    Z = linkage(distances, method=linkage_method)

    if n_clusters is not None:
        labels = fcluster(Z, n_clusters, criterion='maxclust') - 1
    else:
        labels = fcluster(Z, 0.5 * np.max(Z[:, 2]), criterion='distance') - 1

    return labels, Z

def dbscan_clustering(X, eps=0.5, min_samples=5):
    """DBSCAN clustering"""
    D = compute_pairwise_distances(X, metric='euclidean')

    labels = np.full(X.shape[0], -1)  # -1 = noise
    cluster_id = 0

    for i in range(X.shape[0]):
        if labels[i] != -1:
            continue

        # Find neighbors
        neighbors = np.where(D[i] <= eps)[0]

        if len(neighbors) < min_samples:
            continue  # Mark as noise for now

        # Start cluster
        labels[i] = cluster_id
        queue = list(neighbors)

        while queue:
            j = queue.pop(0)
            if labels[j] == -1:
                labels[j] = cluster_id
                neighbors_j = np.where(D[j] <= eps)[0]
                if len(neighbors_j) >= min_samples:
                    queue.extend(neighbors_j)

        cluster_id += 1

    return labels

def consensus_clustering(X, k, n_iterations=100, subsample_frac=0.8):
    """Consensus clustering for robustness assessment"""
    n_samples = X.shape[0]
    consensus_matrix = np.zeros((n_samples, n_samples))

    for iteration in range(n_iterations):
        # Subsample
        sample_idx = np.random.choice(n_samples, size=int(n_samples * subsample_frac), replace=False)
        X_sub = X[sample_idx]

        # Cluster subsampled data
        labels_sub, _ = kmeans(X_sub, k)

        # Update consensus matrix
        for i_idx, i in enumerate(sample_idx):
            for j_idx, j in enumerate(sample_idx):
                if labels_sub[i_idx] == labels_sub[j_idx]:
                    consensus_matrix[i, j] += 1

    # Normalize
    consensus_matrix /= n_iterations

    return consensus_matrix

def silhouette_score(X, labels):
    """Compute silhouette score per sample"""
    D = compute_pairwise_distances(X, metric='euclidean')

    silhouette_scores = np.zeros(len(labels))

    for i in range(len(labels)):
        # Intra-cluster distance
        same_cluster = labels == labels[i]
        if np.sum(same_cluster) > 1:
            a_i = np.mean(D[i, same_cluster])
        else:
            a_i = 0

        # Inter-cluster distance
        other_clusters = np.unique(labels[labels != labels[i]])
        if len(other_clusters) > 0:
            b_i = np.min([np.mean(D[i, labels == c]) for c in other_clusters])
        else:
            b_i = 0

        if max(a_i, b_i) > 0:
            silhouette_scores[i] = (b_i - a_i) / max(a_i, b_i)

    return silhouette_scores

def plot_silhouette(X, labels, outfile="silhouette.png"):
    """Plot silhouette scores"""
    silhouette_scores = silhouette_score(X, labels)

    fig, ax = plt.subplots(figsize=(10, max(6, len(labels) * 0.02)), dpi=300)

    # Sort by cluster first, then by silhouette score within cluster
    sorted_idx = np.lexsort((silhouette_scores, labels))
    sorted_labels = labels[sorted_idx]
    sorted_scores = silhouette_scores[sorted_idx]

    y_pos = np.arange(len(sorted_scores))
    n_clusters = len(np.unique(sorted_labels))
    cmap_colors = plt.cm.tab20(np.linspace(0, 1, max(n_clusters, 2)))
    # Build a list of colors, one per sample
    cluster_colors = [cmap_colors[int(lbl) % len(cmap_colors)] for lbl in sorted_labels]

    ax.barh(y_pos, sorted_scores, color=cluster_colors, alpha=0.7)

    ax.set_xlabel('Silhouette Score', fontsize=12)
    ax.set_title('Silhouette Plot', fontsize=14, fontweight='bold')
    ax.axvline(x=np.mean(silhouette_scores), color='red', linestyle='--', label=f'Mean: {np.mean(silhouette_scores):.2f}')
    ax.legend()

    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    print(f"Saved {outfile}")
    plt.close()

def plot_heatmap(X, labels, metadata=None, color_by=None, outfile="heatmap.png"):
    """Plot clustered heatmap"""
    # Sort by cluster
    sorted_idx = np.argsort(labels)
    X_sorted = X[sorted_idx]
    labels_sorted = labels[sorted_idx]

    # Normalize columns for visualization
    X_norm = (X_sorted - X_sorted.mean(axis=0)) / (X_sorted.std(axis=0) + 1e-6)
    X_norm = np.clip(X_norm, -3, 3)

    fig, ax = plt.subplots(figsize=(15, max(6, len(X_norm) * 0.05)), dpi=300)

    im = ax.imshow(X_norm, cmap='RdBu_r', aspect='auto', interpolation='nearest')

    # Draw cluster boundaries
    cluster_boundaries = np.where(np.diff(labels_sorted) != 0)[0]
    for boundary in cluster_boundaries:
        ax.axhline(y=boundary, color='black', linewidth=1)

    ax.set_xlabel('Features', fontsize=12)
    ax.set_ylabel('Samples', fontsize=12)
    ax.set_title('Hierarchical Clustered Heatmap', fontsize=14, fontweight='bold')

    plt.colorbar(im, ax=ax, label='Normalized Expression')
    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    print(f"Saved {outfile}")
    plt.close()

def plot_dendrogram(Z, outfile="dendrogram.png"):
    """Plot hierarchical dendrogram"""
    fig, ax = plt.subplots(figsize=(15, 8), dpi=300)

    dendrogram(Z, ax=ax, no_labels=True)

    ax.set_xlabel('Sample Index', fontsize=12)
    ax.set_ylabel('Distance', fontsize=12)
    ax.set_title('Hierarchical Clustering Dendrogram', fontsize=14, fontweight='bold')

    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    print(f"Saved {outfile}")
    plt.close()

def plot_elbow(X, k_range, outfile="elbow.png"):
    """Plot silhouette scores for k selection"""
    silhouette_scores = []

    for k in range(k_range[0], k_range[1] + 1):
        labels, _ = kmeans(X, k)
        sil_score = np.mean(silhouette_score(X, labels))
        silhouette_scores.append(sil_score)

    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)

    ax.plot(range(k_range[0], k_range[1] + 1), silhouette_scores, 'bo-', linewidth=2, markersize=8)
    ax.set_xlabel('Number of Clusters (k)', fontsize=12)
    ax.set_ylabel('Mean Silhouette Score', fontsize=12)
    ax.set_title('Silhouette Score for k Selection', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    print(f"Saved {outfile}")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Clustering Analysis for Omics Data')
    parser.add_argument('--input', required=True, help='Input TSV matrix')
    parser.add_argument('--method', choices=['kmeans', 'hierarchical', 'dbscan', 'consensus', 'all'],
                       default='hierarchical')
    parser.add_argument('--n-clusters', type=int, default=0, help='0 = auto-select')
    parser.add_argument('--k-range', default='2,10', help='k range for auto-selection (min,max)')
    parser.add_argument('--linkage', choices=['ward', 'complete', 'average', 'single'], default='ward')
    parser.add_argument('--distance-metric', choices=['euclidean', 'correlation', 'cosine'],
                       default='euclidean')
    parser.add_argument('--eps', type=float, default=0.5, help='DBSCAN epsilon')
    parser.add_argument('--min-samples', type=int, default=5, help='DBSCAN min_samples')
    parser.add_argument('--scale', choices=['standard', 'minmax', 'none'], default='standard')
    parser.add_argument('--transpose', action='store_true', help='Cluster features instead of samples')
    parser.add_argument('--metadata', help='Sample metadata TSV')
    parser.add_argument('--color-by', help='Metadata column for heatmap colors')
    parser.add_argument('--consensus-n', type=int, default=100)
    parser.add_argument('--consensus-subsample', type=float, default=0.8)
    parser.add_argument('--outdir', default='results')

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Load data
    df = load_and_validate_data(args.input, args.transpose)
    X = df.values
    sample_ids = np.array(df.index)

    # Load metadata
    metadata = None
    if args.metadata:
        metadata = load_metadata(args.metadata)

    # Scale
    X = scale_data(X, args.scale)

    # Parse k_range
    k_range = list(map(int, args.k_range.split(',')))

    # Auto k-selection
    n_clusters = args.n_clusters
    if n_clusters == 0:
        print("Auto-selecting optimal k...")
        plot_elbow(X, k_range, f'{args.outdir}/elbow_plot.png')
        silhouette_scores = []
        for k in range(k_range[0], k_range[1] + 1):
            labels, _ = kmeans(X, k)
            sil = np.mean(silhouette_score(X, labels))
            silhouette_scores.append(sil)
        n_clusters = k_range[0] + np.argmax(silhouette_scores)
        print(f"Selected k={n_clusters}")

    # Clustering
    if args.method in ['kmeans', 'all']:
        print("Running K-means...")
        labels, centroids = kmeans(X, n_clusters)

        # Save labels
        sil_scores = silhouette_score(X, labels)
        out_df = pd.DataFrame({
            'cluster': labels,
            'silhouette_score': sil_scores
        }, index=sample_ids)
        out_df.to_csv(f'{args.outdir}/kmeans_labels.tsv', sep='\t')

        # Plots
        plot_silhouette(X, labels, f'{args.outdir}/kmeans_silhouette.png')
        plot_heatmap(X, labels, metadata, args.color_by, f'{args.outdir}/kmeans_heatmap.png')

    if args.method in ['hierarchical', 'all']:
        print("Running Hierarchical Clustering...")
        labels, Z = hierarchical_clustering(X, args.linkage, args.distance_metric, n_clusters)

        sil_scores = silhouette_score(X, labels)
        out_df = pd.DataFrame({
            'cluster': labels,
            'silhouette_score': sil_scores
        }, index=sample_ids)
        out_df.to_csv(f'{args.outdir}/hierarchical_labels.tsv', sep='\t')

        plot_silhouette(X, labels, f'{args.outdir}/hierarchical_silhouette.png')
        plot_heatmap(X, labels, metadata, args.color_by, f'{args.outdir}/hierarchical_heatmap.png')
        plot_dendrogram(Z, f'{args.outdir}/dendrogram.png')

    if args.method in ['dbscan', 'all']:
        print("Running DBSCAN...")
        labels = dbscan_clustering(X, args.eps, args.min_samples)

        sil_scores = silhouette_score(X, labels)
        out_df = pd.DataFrame({
            'cluster': labels,
            'silhouette_score': sil_scores
        }, index=sample_ids)
        out_df.to_csv(f'{args.outdir}/dbscan_labels.tsv', sep='\t')

        plot_silhouette(X, labels, f'{args.outdir}/dbscan_silhouette.png')
        plot_heatmap(X, labels, metadata, args.color_by, f'{args.outdir}/dbscan_heatmap.png')

    if args.method in ['consensus', 'all']:
        print("Running Consensus Clustering...")
        consensus_mat = consensus_clustering(X, n_clusters, args.consensus_n, args.consensus_subsample)

        # Cluster consensus matrix
        labels, _ = hierarchical_clustering(consensus_mat, 'average', 'euclidean', n_clusters)

        sil_scores = silhouette_score(X, labels)
        out_df = pd.DataFrame({
            'cluster': labels,
            'silhouette_score': sil_scores
        }, index=sample_ids)
        out_df.to_csv(f'{args.outdir}/consensus_labels.tsv', sep='\t')

        # Plot consensus matrix
        fig, ax = plt.subplots(figsize=(12, 10), dpi=300)
        im = ax.imshow(consensus_mat, cmap='viridis', aspect='auto')
        ax.set_title('Consensus Matrix', fontsize=14, fontweight='bold')
        plt.colorbar(im, ax=ax, label='Co-clustering Frequency')
        plt.tight_layout()
        plt.savefig(f'{args.outdir}/consensus_matrix.png', dpi=300, bbox_inches='tight')
        print(f"Saved {args.outdir}/consensus_matrix.png")
        plt.close()

        plot_silhouette(X, labels, f'{args.outdir}/consensus_silhouette.png')
        plot_heatmap(X, labels, metadata, args.color_by, f'{args.outdir}/consensus_heatmap.png')

    print("Done!")

if __name__ == '__main__':
    main()
