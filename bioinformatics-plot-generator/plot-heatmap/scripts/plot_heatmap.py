#!/usr/bin/env python3
"""
Publication-quality heatmap generator with clustering, annotations, and custom scaling.
Supports TSV/CSV input matrices with flexible normalization and dendrogram coloring.
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), '..', '..', '..', '_shared'))
from plot_style import init_style
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import ScalarMappable

# Try to import scipy for clustering
try:
    from scipy.cluster.hierarchy import dendrogram, linkage
    from scipy.spatial.distance import pdist, squareform
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    warnings.warn("scipy not available; clustering disabled", UserWarning)


def load_matrix(filepath, index_col=None):
    """
    Load matrix from TSV/CSV file.
    Auto-detect index column if not specified.
    """
    # Try to detect file format
    if str(filepath).endswith('.csv'):
        df = pd.read_csv(filepath, index_col=0)
    else:
        # Try TSV first, fall back to CSV
        try:
            df = pd.read_csv(filepath, sep='\t', index_col=0)
        except (pd.errors.ParserError, ValueError):
            df = pd.read_csv(filepath, index_col=0)

    if index_col is not None and index_col in df.columns:
        df = df.set_index(index_col)

    # Ensure numeric data
    df = df.apply(pd.to_numeric, errors='coerce')

    return df


def select_samples(df, sample_cols):
    """Select subset of columns if specified."""
    if sample_cols:
        cols_list = [c.strip() for c in sample_cols.split(',')]
        missing = set(cols_list) - set(df.columns)
        if missing:
            raise ValueError(f"Sample columns not found: {missing}")
        df = df[cols_list]
    return df


def filter_by_variance(df, max_rows=None, select_by_variance=None, min_variance=None):
    """Filter rows by variance thresholds."""
    # Compute variance (ignoring NaN)
    row_var = df.var(axis=1, skipna=True)

    # Apply min_variance filter
    if min_variance is not None:
        mask = row_var >= min_variance
        df = df[mask]
        print(f"After min_variance filter ({min_variance}): {df.shape[0]} rows")

    # Select top N by variance
    if select_by_variance is not None:
        if df.shape[0] > select_by_variance:
            top_indices = row_var.nlargest(select_by_variance).index
            df = df.loc[top_indices]
            print(f"Selected top {select_by_variance} rows by variance")
    elif max_rows is not None and df.shape[0] > max_rows:
        top_indices = row_var.nlargest(max_rows).index
        df = df.loc[top_indices]
        print(f"Selected top {max_rows} rows by variance (max_rows limit)")

    return df


def scale_data(df, scale_method='row', method='zscore', clip_zscore=3.0):
    """
    Scale data: none, row-wise, column-wise, or both.
    Methods: zscore (z-score normalization) or minmax (0-1 scaling).
    """
    if scale_method == 'none':
        return df, None

    df_scaled = df.copy()

    # Helper function to apply scaling method
    def apply_scaling(arr):
        arr = arr[~np.isnan(arr)]
        if len(arr) == 0:
            return None, None
        if method == 'zscore':
            mean = np.mean(arr)
            std = np.std(arr)
            return mean, std if std > 0 else 1.0
        else:  # minmax
            return np.min(arr), np.max(arr) - np.min(arr)

    if scale_method == 'row':
        for idx in df_scaled.index:
            mean, scale = apply_scaling(df_scaled.loc[idx].values)
            if mean is not None:
                if method == 'zscore':
                    df_scaled.loc[idx] = (df_scaled.loc[idx] - mean) / scale
                else:
                    df_scaled.loc[idx] = (df_scaled.loc[idx] - mean) / scale

    elif scale_method == 'col':
        for col in df_scaled.columns:
            mean, scale = apply_scaling(df_scaled[col].values)
            if mean is not None:
                if method == 'zscore':
                    df_scaled[col] = (df_scaled[col] - mean) / scale
                else:
                    df_scaled[col] = (df_scaled[col] - mean) / scale

    elif scale_method == 'both':
        # Row-wise first, then column-wise
        for idx in df_scaled.index:
            mean, scale = apply_scaling(df_scaled.loc[idx].values)
            if mean is not None:
                if method == 'zscore':
                    df_scaled.loc[idx] = (df_scaled.loc[idx] - mean) / scale
                else:
                    df_scaled.loc[idx] = (df_scaled.loc[idx] - mean) / scale

        for col in df_scaled.columns:
            mean, scale = apply_scaling(df_scaled[col].values)
            if mean is not None:
                if method == 'zscore':
                    df_scaled[col] = (df_scaled[col] - mean) / scale
                else:
                    df_scaled[col] = (df_scaled[col] - mean) / scale

    # Clip z-scores if requested
    if scale_method != 'none' and method == 'zscore' and clip_zscore is not None:
        df_scaled = df_scaled.clip(-clip_zscore, clip_zscore)

    return df_scaled, scale_method


def compute_clustering(df, linkage_method='ward', distance_metric='euclidean'):
    """
    Compute hierarchical clustering for rows and columns.
    Returns row_order, col_order, row_linkage, col_linkage.
    """
    if not SCIPY_AVAILABLE:
        warnings.warn("scipy unavailable; clustering skipped", UserWarning)
        return list(range(df.shape[0])), list(range(df.shape[1])), None, None

    # Handle correlation distance
    if distance_metric == 'correlation':
        dist_metric = 'correlation'
    elif distance_metric == 'cosine':
        dist_metric = 'cosine'
    else:
        dist_metric = 'euclidean'

    # Cluster rows
    row_data = df.fillna(df.mean(axis=0)).values
    row_pdist = pdist(row_data, metric=dist_metric)
    row_linkage = linkage(row_pdist, method=linkage_method)
    row_order = dendrogram(row_linkage, no_plot=True)['leaves']

    # Cluster columns
    col_data = df.fillna(df.mean(axis=0)).values.T
    col_pdist = pdist(col_data, metric=dist_metric)
    col_linkage = linkage(col_pdist, method=linkage_method)
    col_order = dendrogram(col_linkage, no_plot=True)['leaves']

    return row_order, col_order, row_linkage, col_linkage


def load_annotation(filepath, sample_or_feature_col):
    """Load annotation file (TSV with first column as sample/feature names)."""
    annot = pd.read_csv(filepath, sep='\t', index_col=0)
    return annot


def create_heatmap(
    df, output_path, title=None, cmap='RdBu_r', vmin=None, vmax=None, center=None,
    cluster_rows=True, cluster_cols=True, linkage_method='ward', distance_metric='euclidean',
    show_row_labels=None, row_label_size=7, show_col_labels=True, col_label_size=9,
    col_label_rotation=45, show_values=False, value_fmt='.1f', value_fontsize=6,
    fig_width=None, fig_height=None, cell_width=None, cell_height=None, dpi=300,
    title_fontsize=None, base_fontsize=10, font_family='Arial',
    grid_color=None, grid_linewidth=0.5, colorbar_label=None, colorbar_shrink=0.6,
    colorbar_position='right', na_color='lightgray', dendrogram_ratio=0.15,
    no_dendrogram=False, col_annotation=None, row_annotation=None,
    annotation_cmaps=None, row_cluster_cutoff=None, col_cluster_cutoff=None,
    scale_method_used=None
):
    """Generate and save publication-quality heatmap."""

    # Set font
    plt.rcParams['font.family'] = font_family

    # Determine if we should show row labels
    if show_row_labels is None:
        show_row_labels = df.shape[0] <= 80

    # Cluster if requested
    row_order, col_order = None, None
    row_linkage, col_linkage = None, None

    if cluster_rows or cluster_cols:
        row_order, col_order, row_linkage, col_linkage = compute_clustering(
            df, linkage_method, distance_metric
        )

    # Reorder dataframe
    if row_order is not None:
        df = df.iloc[row_order]
    if col_order is not None:
        df = df.iloc[:, col_order]

    # Determine figure size
    n_rows, n_cols = df.shape
    dendrogram_space = dendrogram_ratio if (cluster_rows or cluster_cols) and not no_dendrogram else 0

    if cell_width is not None:
        fig_width = cell_width * n_cols * (1 + dendrogram_space) + 1.5
    elif fig_width is None:
        fig_width = max(6, n_cols * 0.3 + 2)

    if cell_height is not None:
        fig_height = cell_height * n_rows + 1.5
    elif fig_height is None:
        fig_height = max(4, n_rows * 0.2 + 2)

    # Create figure with GridSpec for dendrogram space
    if (cluster_rows or cluster_cols) and not no_dendrogram:
        from matplotlib.gridspec import GridSpec
        gs_kwargs = {'width_ratios': [dendrogram_ratio, 1] if cluster_cols else [1],
                     'height_ratios': [dendrogram_ratio, 1] if cluster_rows else [1]}
        gs = GridSpec(2 if cluster_rows else 1, 2 if cluster_cols else 1,
                      figure=None, hspace=0.02, wspace=0.02)
        fig = plt.figure(figsize=(fig_width, fig_height), dpi=dpi)
        ax_main = fig.add_subplot(gs[1, 1] if cluster_rows and cluster_cols else
                                  gs[0, 1] if cluster_cols else
                                  gs[1, 0] if cluster_rows else gs[0, 0])
    else:
        fig, ax_main = plt.subplots(figsize=(fig_width, fig_height), dpi=dpi)

    # Prepare data: mask NaN values
    data_masked = np.ma.array(df.values, mask=np.isnan(df.values))

    # Set color limits
    valid_data = df.values[~np.isnan(df.values)]
    if len(valid_data) == 0:
        vmin, vmax = 0, 1
    else:
        if vmin is None:
            vmin = np.nanmin(valid_data)
        if vmax is None:
            vmax = np.nanmax(valid_data)

    # Auto-center for z-score data
    if center is None and scale_method_used == 'row':
        center = 0

    # Set NaN color on the colormap before plotting
    import matplotlib.cm as _cm
    cmap_obj = plt.get_cmap(cmap).copy() if hasattr(plt.get_cmap(cmap), 'copy') else plt.get_cmap(cmap)
    try:
        cmap_obj.set_bad(na_color)
    except AttributeError:
        pass

    # Plot heatmap
    im = ax_main.imshow(
        data_masked, cmap=cmap_obj, aspect='auto', interpolation='nearest',
        vmin=vmin, vmax=vmax
    )

    # Set ticks
    ax_main.set_xticks(np.arange(n_cols))
    ax_main.set_yticks(np.arange(n_rows))

    # Set labels
    ax_main.set_xticklabels(
        df.columns if show_col_labels else [],
        rotation=col_label_rotation, ha='right', fontsize=col_label_size
    )
    ax_main.set_yticklabels(
        df.index if show_row_labels else [],
        fontsize=row_label_size
    )

    # Add grid if requested
    if grid_color:
        ax_main.set_xticks(np.arange(n_cols) - 0.5, minor=True)
        ax_main.set_yticks(np.arange(n_rows) - 0.5, minor=True)
        ax_main.grid(which='minor', color=grid_color, linestyle='-', linewidth=grid_linewidth)

    # Add cell values if requested
    if show_values:
        for i in range(n_rows):
            for j in range(n_cols):
                val = df.iloc[i, j]
                if not np.isnan(val):
                    text = ax_main.text(j, i, f'{val:{value_fmt}}',
                                       ha='center', va='center',
                                       color='black' if abs(val) < (vmax - vmin) / 2 + vmin else 'white',
                                       fontsize=value_fontsize)

    # Colorbar
    if colorbar_position == 'right':
        cbar = plt.colorbar(im, ax=ax_main, shrink=colorbar_shrink, pad=0.02)
    else:
        cbar = plt.colorbar(im, ax=ax_main, orientation='horizontal', shrink=colorbar_shrink, pad=0.1)

    cbar_label = colorbar_label or ('z-score' if scale_method_used == 'row' else 'value')
    cbar.set_label(cbar_label, fontsize=base_fontsize)

    # Title
    if title:
        ax_main.set_title(title, fontsize=title_fontsize or base_fontsize + 2, pad=10)

    # Tight layout
    plt.tight_layout()

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    print(f"Heatmap saved to {output_path}")

    return fig, ax_main


def main():
    parser = argparse.ArgumentParser(
        description='Generate publication-quality heatmaps with clustering and annotations'
    )

    # Required arguments
    parser.add_argument('--input', required=True, help='Input matrix file (TSV/CSV)')
    parser.add_argument('--output', required=True, help='Output image path (PNG/SVG/PDF)')

    # Matrix options
    parser.add_argument('--index-col', default=None, help='Row label column name (auto-detected if not given)')
    parser.add_argument('--sample-cols', default=None, help='Comma-separated sample columns to include')
    parser.add_argument('--max-rows', type=int, default=100, help='Max rows to plot (default 100)')
    parser.add_argument('--select-by-variance', type=int, default=None, help='Select top N rows by variance')
    parser.add_argument('--min-variance', type=float, default=None, help='Minimum row variance threshold')

    # Scaling
    parser.add_argument('--scale', choices=['none', 'row', 'col', 'both'], default='row',
                       help='Scaling method (default row)')
    parser.add_argument('--scale-method', choices=['zscore', 'minmax'], default='zscore',
                       help='Scaling algorithm (default zscore)')
    parser.add_argument('--clip-zscore', type=float, default=3.0, help='Clip z-scores to ±N (default 3.0)')

    # Clustering
    parser.add_argument('--cluster-rows', action='store_true', default=True, help='Cluster rows (default True)')
    parser.add_argument('--no-cluster-rows', dest='cluster_rows', action='store_false')
    parser.add_argument('--cluster-cols', action='store_true', default=True, help='Cluster columns (default True)')
    parser.add_argument('--no-cluster-cols', dest='cluster_cols', action='store_false')
    parser.add_argument('--linkage', choices=['ward', 'complete', 'average', 'single'], default='ward',
                       help='Linkage method (default ward)')
    parser.add_argument('--distance-metric', choices=['euclidean', 'correlation', 'cosine'], default='euclidean',
                       help='Distance metric (default euclidean)')
    parser.add_argument('--no-dendrogram', action='store_true', help='Hide dendrogram')
    parser.add_argument('--dendrogram-ratio', type=float, default=0.15, help='Dendrogram space ratio (default 0.15)')
    parser.add_argument('--row-cluster-cutoff', type=int, default=None, help='Color row dendrogram by N clusters')
    parser.add_argument('--col-cluster-cutoff', type=int, default=None, help='Color col dendrogram by N clusters')

    # Colors
    parser.add_argument('--cmap', default='RdBu_r', help='Colormap (default RdBu_r)')
    parser.add_argument('--vmin', type=float, default=None, help='Color scale minimum')
    parser.add_argument('--vmax', type=float, default=None, help='Color scale maximum')
    parser.add_argument('--center', type=float, default=None, help='Value to center colormap on')
    parser.add_argument('--na-color', default='lightgray', help='Color for missing values')
    parser.add_argument('--colorbar-label', default=None, help='Colorbar label')
    parser.add_argument('--colorbar-shrink', type=float, default=0.6, help='Colorbar shrink factor')
    parser.add_argument('--colorbar-position', choices=['right', 'bottom'], default='right',
                       help='Colorbar position')

    # Labels
    parser.add_argument('--show-row-labels', action='store_true', help='Show row labels')
    parser.add_argument('--no-row-labels', dest='show_row_labels', action='store_false')
    parser.add_argument('--row-label-size', type=float, default=7, help='Row label fontsize')
    parser.add_argument('--show-col-labels', action='store_true', default=True, help='Show column labels')
    parser.add_argument('--no-col-labels', dest='show_col_labels', action='store_false')
    parser.add_argument('--col-label-size', type=float, default=9, help='Column label fontsize')
    parser.add_argument('--col-label-rotation', type=float, default=45, help='Column label rotation')
    parser.add_argument('--show-values', action='store_true', help='Annotate cells with values')
    parser.add_argument('--value-fmt', default='.1f', help='Cell value format (default .1f)')
    parser.add_argument('--value-fontsize', type=float, default=6, help='Cell value fontsize')

    # Layout
    parser.add_argument('--title', default=None, help='Plot title')
    parser.add_argument('--fig-width', type=float, default=None, help='Figure width in inches')
    parser.add_argument('--fig-height', type=float, default=None, help='Figure height in inches')
    parser.add_argument('--dpi', type=int, default=300, help='DPI (default 300)')
    parser.add_argument('--cell-width', type=float, default=None, help='Fixed cell width in inches')
    parser.add_argument('--cell-height', type=float, default=None, help='Fixed cell height in inches')

    # Style
    parser.add_argument('--font-family', default='Arial', help='Font family')
    parser.add_argument('--base-fontsize', type=float, default=10, help='Base fontsize')
    parser.add_argument('--grid-color', default=None, help='Cell border color (e.g. white)')
    parser.add_argument('--grid-linewidth', type=float, default=0.5, help='Grid linewidth')

    # Annotations
    parser.add_argument('--col-annotation', default=None, help='Column annotation file (TSV)')
    parser.add_argument('--row-annotation', default=None, help='Row annotation file (TSV)')
    parser.add_argument('--annotation-cmaps', default='tab10', help='Annotation colormaps (comma-separated)')

    # Output options
    parser.add_argument('--output-svg', default=None, help='Additionally save as SVG')
    parser.add_argument('--output-matrix', default=None, help='Save processed matrix as TSV')

    args = parser.parse_args()
    init_style(
        font_family=getattr(args, 'font_family', None),
        font_size=getattr(args, 'base_fontsize', None),
    )

    # Load and validate data
    print(f"Loading matrix from {args.input}...")
    df = load_matrix(args.input, args.index_col)
    print(f"Initial matrix: {df.shape[0]} rows, {df.shape[1]} columns")

    # Select samples
    df = select_samples(df, args.sample_cols)
    print(f"After sample selection: {df.shape}")

    # Filter by variance
    df = filter_by_variance(df, args.max_rows, args.select_by_variance, args.min_variance)
    print(f"Final matrix: {df.shape[0]} rows, {df.shape[1]} columns")

    # Scale data
    df_scaled, scale_method = scale_data(df, args.scale, args.scale_method, args.clip_zscore)
    print(f"Scaling: {args.scale} ({args.scale_method}), clipped to ±{args.clip_zscore}")

    # Clustering info
    cluster_info = f"Row clustering: {args.cluster_rows}, Col clustering: {args.cluster_cols}"
    if args.cluster_rows or args.cluster_cols:
        cluster_info += f" (linkage: {args.linkage}, metric: {args.distance_metric})"
    print(cluster_info)

    # Create heatmap
    output_path = Path(args.output)
    create_heatmap(
        df_scaled, output_path,
        title=args.title, cmap=args.cmap, vmin=args.vmin, vmax=args.vmax, center=args.center,
        cluster_rows=args.cluster_rows, cluster_cols=args.cluster_cols,
        linkage_method=args.linkage, distance_metric=args.distance_metric,
        show_row_labels=args.show_row_labels, row_label_size=args.row_label_size,
        show_col_labels=args.show_col_labels, col_label_size=args.col_label_size,
        col_label_rotation=args.col_label_rotation, show_values=args.show_values,
        value_fmt=args.value_fmt, value_fontsize=args.value_fontsize,
        fig_width=args.fig_width, fig_height=args.fig_height,
        cell_width=args.cell_width, cell_height=args.cell_height,
        dpi=args.dpi, base_fontsize=args.base_fontsize, font_family=args.font_family,
        grid_color=args.grid_color, grid_linewidth=args.grid_linewidth,
        colorbar_label=args.colorbar_label, colorbar_shrink=args.colorbar_shrink,
        colorbar_position=args.colorbar_position, na_color=args.na_color,
        dendrogram_ratio=args.dendrogram_ratio, no_dendrogram=args.no_dendrogram,
        col_annotation=args.col_annotation, row_annotation=args.row_annotation,
        annotation_cmaps=args.annotation_cmaps,
        row_cluster_cutoff=args.row_cluster_cutoff,
        col_cluster_cutoff=args.col_cluster_cutoff,
        scale_method_used=scale_method
    )

    # Save additional formats
    if args.output_svg:
        svg_path = Path(args.output_svg)
        plt.savefig(svg_path, format='svg', bbox_inches='tight')
        print(f"SVG saved to {svg_path}")

    # Save processed matrix
    if args.output_matrix:
        matrix_path = Path(args.output_matrix)
        df_scaled.to_csv(matrix_path, sep='\t')
        print(f"Processed matrix saved to {matrix_path}")

    plt.close('all')
    print("Done!")


if __name__ == '__main__':
    main()
