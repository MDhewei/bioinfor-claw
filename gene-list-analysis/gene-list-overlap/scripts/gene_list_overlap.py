#!/usr/bin/env python3
"""
Gene List Overlap & Enrichment Analysis
Compare 2-6 gene lists for overlap, compute Fisher's exact test for enrichment,
generate Venn diagrams, UpSet plots, and Jaccard heatmap.
"""

import argparse
import sys
import os
from pathlib import Path
from collections import defaultdict
import math
from itertools import combinations

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), '..', '..', '..', '_shared'))
from plot_style import init_style
import matplotlib.patches as patches
from matplotlib.patches import Circle, Wedge
try:
    from scipy.cluster.hierarchy import dendrogram, linkage
    from scipy.spatial.distance import pdist, squareform
    _HAVE_SCIPY = True
except ImportError:
    _HAVE_SCIPY = False

    def pdist(X, metric='euclidean'):
        X = np.asarray(X, float)
        D_sq = np.sum(X ** 2, axis=1)
        D = np.sqrt(np.maximum(D_sq[:, None] + D_sq[None, :] - 2 * X @ X.T, 0))
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

    def linkage(distances, method='ward'):
        D = squareform(distances) if np.asarray(distances).ndim == 1 else np.asarray(distances, float)
        n = len(D)
        Z = []
        members = {i: [i] for i in range(n)}
        active = list(range(n))
        cur_D = D.copy()
        for _ in range(n - 1):
            aa = np.array(active)
            min_d, mi, mj = np.inf, -1, -1
            for ii in range(len(aa)):
                for jj in range(ii + 1, len(aa)):
                    if cur_D[aa[ii], aa[jj]] < min_d:
                        min_d = cur_D[aa[ii], aa[jj]]
                        mi, mj = aa[ii], aa[jj]
            if mi < 0:
                break
            ni, nj = len(members[mi]), len(members[mj])
            Z.append([float(mi), float(mj), float(min_d), float(ni + nj)])
            new_id = n + len(Z) - 1
            members[new_id] = members[mi] + members[mj]
            active.remove(mi); active.remove(mj)
            # Update distances (average linkage)
            for k in active:
                dk_i = cur_D[k, mi]; dk_j = cur_D[k, mj]
                nk = len(members[k]) if k in members else 1
                if method == 'single':
                    d_new = min(dk_i, dk_j)
                elif method == 'complete':
                    d_new = max(dk_i, dk_j)
                elif method == 'ward':
                    d_new = np.sqrt(max(((nk+ni)*dk_i**2 + (nk+nj)*dk_j**2 - nk*min_d**2)/(nk+ni+nj), 0))
                else:
                    d_new = (ni*dk_i + nj*dk_j)/(ni+nj)
                cur_D[k, mi] = d_new; cur_D[mi, k] = d_new
            active.append(mi)
        return np.array(Z) if Z else np.zeros((0, 4))

    def dendrogram(Z, no_plot=False, **kwargs):
        n = len(Z) + 1
        leaves = list(range(n))
        if not no_plot:
            ax = kwargs.get('ax', None)
            if ax:
                ax.set_title('Dendrogram (simplified)')
        return {'icoord': [], 'dcoord': [], 'leaves': leaves, 'ivl': [str(i) for i in leaves]}


def load_gene_list(source):
    """Load genes from file or comma-separated string."""
    if os.path.isfile(source):
        with open(source, 'r') as f:
            genes = [line.strip() for line in f if line.strip()]
    else:
        genes = [g.strip() for g in source.split(',') if g.strip()]
    return set(genes)


def parse_lists_input(lists_str):
    """Parse --lists argument into dict of {label: gene_set}."""
    lists_dict = {}
    items = lists_str.split(',')
    current_label = None
    genes_buffer = []

    for item in items:
        if ':' in item:
            # Save previous label if any
            if current_label:
                source = ','.join(genes_buffer)
                lists_dict[current_label] = load_gene_list(source)
                genes_buffer = []
            # Parse new label:source
            parts = item.split(':', 1)
            current_label = parts[0].strip()
            genes_buffer.append(parts[1].strip())
        else:
            genes_buffer.append(item.strip())

    if current_label:
        source = ','.join(genes_buffer)
        lists_dict[current_label] = load_gene_list(source)

    return lists_dict


def fisher_exact_test(a, b, c, d):
    """
    Compute Fisher's exact test (hypergeometric).
    Inputs: 2x2 contingency table
        In B    Not in B
    In A        a       b
    Not in A    c       d

    Returns: p-value (one-tailed, hypergeometric)
    """
    if a + b == 0 or c + d == 0:
        return 1.0
    if a + c == 0 or b + d == 0:
        return 1.0

    # Use hypergeometric: P(X >= k) where k = a
    # P(X = k) = C(K, k) * C(N-K, n-k) / C(N, n)
    # where N = a+b+c+d, K = a+c, n = a+b, k = a

    def comb(n, k):
        if k > n or k < 0:
            return 0
        if k == 0 or k == n:
            return 1
        k = min(k, n - k)
        result = 1
        for i in range(k):
            result = result * (n - i) // (i + 1)
        return result

    N = a + b + c + d
    K = a + c
    n = a + b

    p_value = 0.0
    for k in range(max(0, n + K - N), min(n, K) + 1):
        p = (comb(K, k) * comb(N - K, n - k)) / comb(N, n)
        if k >= a:
            p_value += p

    return min(p_value, 1.0)


def odds_ratio(a, b, c, d):
    """Compute odds ratio."""
    if b == 0 or c == 0:
        return float('inf') if a * d > 0 else 0
    return (a * d) / (b * c)


def benjamini_hochberg_fdr(p_values):
    """Apply Benjamini-Hochberg FDR correction."""
    n = len(p_values)
    if n == 0:
        return []

    sorted_indices = np.argsort(p_values)
    sorted_pvals = np.array(p_values)[sorted_indices]

    # Compute adjusted p-values
    adjusted = np.zeros(n)
    for i, idx in enumerate(sorted_indices):
        adjusted[idx] = min(sorted_pvals[i] * n / (i + 1), 1.0)

    # Monotonicity correction
    for i in range(n - 2, -1, -1):
        adjusted[i] = min(adjusted[i], adjusted[i + 1])

    return adjusted.tolist()


def compute_pairwise_overlap(lists_dict, background_size):
    """Compute pairwise overlap statistics."""
    labels = list(lists_dict.keys())
    results = []

    p_values = []

    for i, label_a in enumerate(labels):
        for label_b in labels[i+1:]:
            genes_a = lists_dict[label_a]
            genes_b = lists_dict[label_b]

            intersection = len(genes_a & genes_b)
            union = len(genes_a | genes_b)
            jaccard = intersection / union if union > 0 else 0
            overlap_coeff = intersection / min(len(genes_a), len(genes_b)) if min(len(genes_a), len(genes_b)) > 0 else 0

            # Fisher's exact test
            a = intersection
            b = len(genes_a) - a
            c = len(genes_b) - a
            d = background_size - len(genes_a) - len(genes_b) + a

            p_val = fisher_exact_test(a, b, c, d)
            p_values.append(p_val)
            or_val = odds_ratio(a, b, c, d)

            results.append({
                'List_A': label_a,
                'List_B': label_b,
                'Size_A': len(genes_a),
                'Size_B': len(genes_b),
                'Intersection': intersection,
                'Union': union,
                'Jaccard': round(jaccard, 4),
                'OverlapCoefficient': round(overlap_coeff, 4),
                'PValue': p_val,
                'OddsRatio': round(or_val, 4) if or_val != float('inf') else float('inf'),
                'FDR': 0.0,  # Will fill after BH correction
            })

    # Apply FDR correction
    if p_values:
        fdr_values = benjamini_hochberg_fdr(p_values)
        for i, result in enumerate(results):
            result['FDR'] = round(fdr_values[i], 4)

    return results


def compute_jaccard_matrix(lists_dict):
    """Compute NxN Jaccard index matrix."""
    labels = list(lists_dict.keys())
    n = len(labels)
    matrix = np.zeros((n, n))

    for i, label_a in enumerate(labels):
        for j, label_b in enumerate(labels):
            if i == j:
                matrix[i, j] = 1.0
            else:
                genes_a = lists_dict[label_a]
                genes_b = lists_dict[label_b]
                intersection = len(genes_a & genes_b)
                union = len(genes_a | genes_b)
                jaccard = intersection / union if union > 0 else 0
                matrix[i, j] = jaccard

    return matrix, labels


def get_unique_genes(lists_dict):
    """Get genes unique to each list."""
    results = []
    union_all = set().union(*lists_dict.values())

    for label, genes in lists_dict.items():
        unique = genes - (union_all - genes)
        results.append({
            'List': label,
            'UniqueCount': len(unique),
            'UniqueGenes': ';'.join(sorted(unique)) if unique else '',
        })

    return results


def get_all_intersections(lists_dict):
    """Get all pairwise and global intersections."""
    labels = list(lists_dict.keys())
    results = []

    # All pairwise
    for i, label_a in enumerate(labels):
        for label_b in labels[i+1:]:
            intersection = lists_dict[label_a] & lists_dict[label_b]
            results.append({
                'SetCombination': f'{label_a}_{label_b}',
                'Count': len(intersection),
                'Genes': ';'.join(sorted(intersection)) if intersection else '',
            })

    # Global intersection (all lists)
    global_intersection = set.intersection(*lists_dict.values()) if lists_dict else set()
    results.append({
        'SetCombination': 'All_' + '_'.join(labels),
        'Count': len(global_intersection),
        'Genes': ';'.join(sorted(global_intersection)) if global_intersection else '',
    })

    return results


def plot_venn_diagram(lists_dict, outdir):
    """Generate Venn diagram for 2-4 lists."""
    labels = list(lists_dict.keys())
    n_lists = len(labels)

    if n_lists < 2 or n_lists > 4:
        print(f"Skipping Venn diagram: only works for 2-4 lists (have {n_lists})", file=sys.stderr)
        return

    fig, ax = plt.subplots(figsize=(10, 10))

    if n_lists == 2:
        # Two circles
        circle1 = Circle((0.35, 0.5), 0.25, color='red', alpha=0.5, label=labels[0])
        circle2 = Circle((0.65, 0.5), 0.25, color='blue', alpha=0.5, label=labels[1])
        ax.add_patch(circle1)
        ax.add_patch(circle2)

        genes1 = lists_dict[labels[0]]
        genes2 = lists_dict[labels[1]]
        intersection = genes1 & genes2
        only1 = genes1 - genes2
        only2 = genes2 - genes1

        ax.text(0.25, 0.5, f"{len(only1)}", ha='center', va='center', fontsize=14, fontweight='bold')
        ax.text(0.5, 0.5, f"{len(intersection)}", ha='center', va='center', fontsize=14, fontweight='bold')
        ax.text(0.75, 0.5, f"{len(only2)}", ha='center', va='center', fontsize=14, fontweight='bold')

    elif n_lists == 3:
        # Three circles
        circle1 = Circle((0.3, 0.6), 0.25, color='red', alpha=0.5, label=labels[0])
        circle2 = Circle((0.7, 0.6), 0.25, color='blue', alpha=0.5, label=labels[1])
        circle3 = Circle((0.5, 0.3), 0.25, color='green', alpha=0.5, label=labels[2])
        ax.add_patch(circle1)
        ax.add_patch(circle2)
        ax.add_patch(circle3)

        g1, g2, g3 = lists_dict[labels[0]], lists_dict[labels[1]], lists_dict[labels[2]]
        i12 = g1 & g2 & g3
        i1_only = len(g1 - g2 - g3)
        i2_only = len(g2 - g1 - g3)
        i3_only = len(g3 - g1 - g2)
        i12_only = len((g1 & g2) - g3)
        i13_only = len((g1 & g3) - g2)
        i23_only = len((g2 & g3) - g1)

        ax.text(0.15, 0.6, f"{i1_only}", ha='center', va='center', fontsize=12, fontweight='bold')
        ax.text(0.85, 0.6, f"{i2_only}", ha='center', va='center', fontsize=12, fontweight='bold')
        ax.text(0.5, 0.1, f"{i3_only}", ha='center', va='center', fontsize=12, fontweight='bold')
        ax.text(0.5, 0.65, f"{i12_only}", ha='center', va='center', fontsize=12, fontweight='bold')
        ax.text(0.35, 0.35, f"{i13_only}", ha='center', va='center', fontsize=12, fontweight='bold')
        ax.text(0.65, 0.35, f"{i23_only}", ha='center', va='center', fontsize=12, fontweight='bold')
        ax.text(0.5, 0.45, f"{len(i12)}", ha='center', va='center', fontsize=12, fontweight='bold')

    elif n_lists == 4:
        # Four ellipses (approximate)
        circle1 = Circle((0.25, 0.25), 0.2, color='red', alpha=0.4, label=labels[0])
        circle2 = Circle((0.75, 0.25), 0.2, color='blue', alpha=0.4, label=labels[1])
        circle3 = Circle((0.25, 0.75), 0.2, color='green', alpha=0.4, label=labels[2])
        circle4 = Circle((0.75, 0.75), 0.2, color='yellow', alpha=0.4, label=labels[3])
        ax.add_patch(circle1)
        ax.add_patch(circle2)
        ax.add_patch(circle3)
        ax.add_patch(circle4)

        ax.text(0.1, 0.25, f"{len(lists_dict[labels[0]])}", ha='center', va='center', fontsize=10)
        ax.text(0.9, 0.25, f"{len(lists_dict[labels[1]])}", ha='center', va='center', fontsize=10)
        ax.text(0.1, 0.9, f"{len(lists_dict[labels[2]])}", ha='center', va='center', fontsize=10)
        ax.text(0.9, 0.9, f"{len(lists_dict[labels[3]])}", ha='center', va='center', fontsize=10)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.legend(loc='upper right', fontsize=12)
    ax.set_title('Gene List Overlap Venn Diagram', fontsize=14, fontweight='bold', pad=20)

    outpath = os.path.join(outdir, 'venn_diagram.png')
    plt.savefig(outpath, dpi=300, bbox_inches='tight')
    print(f"Saved Venn diagram to {outpath}")
    plt.close()


def plot_upset_plot(lists_dict, outdir):
    """Generate UpSet plot for >= 3 lists."""
    labels = list(lists_dict.keys())
    n_lists = len(labels)

    if n_lists < 3:
        print(f"UpSet plot requires >= 3 lists (have {n_lists})", file=sys.stderr)
        return

    # Create binary matrix: rows = genes, columns = lists
    all_genes = set().union(*lists_dict.values())
    matrix = pd.DataFrame(0, index=sorted(all_genes), columns=labels)
    for label, genes in lists_dict.items():
        matrix[label] = [1 if gene in genes else 0 for gene in matrix.index]

    # Get intersection sizes for all combinations
    intersections = []
    for r in range(1, n_lists + 1):
        for combo in combinations(labels, r):
            # Find genes in all lists in combo and not in others
            mask = matrix[list(combo)].sum(axis=1) == len(combo)
            other_lists = [l for l in labels if l not in combo]
            if other_lists:
                mask = mask & (matrix[other_lists].sum(axis=1) == 0)
            count = mask.sum()
            if count > 0:
                intersections.append({
                    'combination': ','.join(combo),
                    'count': count,
                })

    # Sort by count descending
    intersections.sort(key=lambda x: x['count'], reverse=True)

    # Plot
    fig, ax = plt.subplots(figsize=(max(12, len(intersections) * 0.5), 8))

    combos = [x['combination'] for x in intersections]
    counts = [x['count'] for x in intersections]

    bars = ax.bar(range(len(combos)), counts, color='steelblue', alpha=0.8, edgecolor='black')

    # Color code by list involvement
    for i, bar in enumerate(bars):
        combo_lists = combos[i].split(',')
        if len(combo_lists) == n_lists:
            bar.set_color('darkgreen')
        elif len(combo_lists) == 1:
            bar.set_color('lightcoral')
        else:
            bar.set_color('steelblue')

    ax.set_xticks(range(len(combos)))
    ax.set_xticklabels(combos, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('Intersection Size', fontsize=11)
    ax.set_xlabel('Set Combination', fontsize=11)
    ax.set_title('UpSet Plot: Gene List Intersections', fontsize=13, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    outpath = os.path.join(outdir, 'upset_plot.png')
    plt.savefig(outpath, dpi=300, bbox_inches='tight')
    print(f"Saved UpSet plot to {outpath}")
    plt.close()


def plot_jaccard_heatmap(matrix, labels, outdir):
    """Plot clustered Jaccard heatmap."""
    fig, ax = plt.subplots(figsize=(10, 9))

    # Compute distances and linkage for clustering
    distances = pdist(matrix, metric='euclidean')
    linkage_matrix = linkage(distances, method='ward')

    # Reorder matrix by clustering
    dendro = dendrogram(linkage_matrix, no_plot=True)
    order = dendro['leaves']

    matrix_ordered = matrix[np.ix_(order, order)]
    labels_ordered = [labels[i] for i in order]

    # Plot heatmap
    im = ax.imshow(matrix_ordered, cmap='RdYlGn', vmin=0, vmax=1, aspect='auto')

    ax.set_xticks(range(len(labels_ordered)))
    ax.set_yticks(range(len(labels_ordered)))
    ax.set_xticklabels(labels_ordered, rotation=45, ha='right')
    ax.set_yticklabels(labels_ordered)

    # Add values to heatmap
    for i in range(len(labels_ordered)):
        for j in range(len(labels_ordered)):
            text = ax.text(j, i, f'{matrix_ordered[i, j]:.2f}',
                          ha="center", va="center", color="black", fontsize=9)

    ax.set_title('Jaccard Index Heatmap (Clustered)', fontsize=13, fontweight='bold')
    cbar = plt.colorbar(im, ax=ax, label='Jaccard Index')

    plt.tight_layout()
    outpath = os.path.join(outdir, 'jaccard_heatmap.png')
    plt.savefig(outpath, dpi=300, bbox_inches='tight')
    print(f"Saved Jaccard heatmap to {outpath}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Gene List Overlap & Enrichment Analysis')
    parser.add_argument('--lists', required=True, action='append',
                       help='Input lists format: label:path_or_genes (e.g., "A:genes1.txt,B:GENE1,GENE2"). Can use multiple times.')
    parser.add_argument('--background', default='20000', help='Background universe size or file')
    parser.add_argument('--outdir', default='./overlap_output', help='Output directory')
    parser.add_argument('--min-overlap', type=int, default=1, help='Minimum overlap to report')
    parser.add_argument('--output-unique', action='store_true', help='Save unique genes per list')
    parser.add_argument('--output-shared', action='store_true', help='Save all pairwise/global overlaps')
    parser.add_argument('--plot-venn', action='store_true', help='Generate Venn diagram (2-4 lists only)')
    parser.add_argument('--plot-upset', action='store_true', help='Generate UpSet plot (3+ lists)')

    args = parser.parse_args()

    init_style(
        font_family=getattr(args, 'font_family', None),
        font_size=getattr(args, 'font_size', None),
    )

    os.makedirs(args.outdir, exist_ok=True)

    # Determine background size
    if os.path.isfile(args.background):
        with open(args.background, 'r') as f:
            background_genes = set(line.strip() for line in f if line.strip())
        background_size = len(background_genes)
    else:
        background_size = int(args.background)

    # Parse input lists
    print("Loading gene lists...")
    lists_dict = {}
    for list_arg in args.lists:
        list_items = parse_lists_input(list_arg)
        lists_dict.update(list_items)

    if len(lists_dict) < 2:
        print("Error: Need at least 2 gene lists", file=sys.stderr)
        return 1

    print(f"Loaded {len(lists_dict)} lists:")
    for label, genes in lists_dict.items():
        print(f"  {label}: {len(genes)} genes")

    # Compute pairwise overlap
    print("Computing pairwise overlaps...")
    overlap_results = compute_pairwise_overlap(lists_dict, background_size)
    overlap_df = pd.DataFrame(overlap_results)
    overlap_df = overlap_df[overlap_df['Intersection'] >= args.min_overlap]
    overlap_df.to_csv(os.path.join(args.outdir, 'overlap_summary.tsv'), sep='\t', index=False)
    print(f"Saved {len(overlap_df)} pairwise comparisons to overlap_summary.tsv")

    # Compute Jaccard matrix
    print("Computing Jaccard matrix...")
    jaccard_matrix, labels_ordered = compute_jaccard_matrix(lists_dict)
    jaccard_df = pd.DataFrame(jaccard_matrix, index=labels_ordered, columns=labels_ordered)
    jaccard_df.to_csv(os.path.join(args.outdir, 'jaccard_matrix.tsv'), sep='\t')
    print("Saved Jaccard matrix to jaccard_matrix.tsv")

    # Optional: unique genes
    if args.output_unique:
        print("Computing unique genes...")
        unique_results = get_unique_genes(lists_dict)
        unique_df = pd.DataFrame(unique_results)
        unique_df.to_csv(os.path.join(args.outdir, 'unique_per_list.tsv'), sep='\t', index=False)
        print("Saved unique genes to unique_per_list.tsv")

    # Optional: all intersections
    if args.output_shared:
        print("Computing all intersections...")
        intersections = get_all_intersections(lists_dict)
        intersections_df = pd.DataFrame(intersections)
        intersections_df.to_csv(os.path.join(args.outdir, 'all_intersections.tsv'), sep='\t', index=False)
        print("Saved all intersections to all_intersections.tsv")

    # Visualizations
    print("Generating visualizations...")

    if args.plot_venn:
        plot_venn_diagram(lists_dict, args.outdir)

    if args.plot_upset:
        plot_upset_plot(lists_dict, args.outdir)

    plot_jaccard_heatmap(jaccard_matrix, list(lists_dict.keys()), args.outdir)

    print(f"\nCompleted! Results saved to {args.outdir}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
