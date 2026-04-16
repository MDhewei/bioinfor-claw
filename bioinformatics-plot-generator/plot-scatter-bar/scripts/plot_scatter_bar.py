#!/usr/bin/env python3
"""
Bioinformatics plot generator: scatter, bar, MA, correlation matrix, and bubble plots.

Supports multiple plot types with publication-quality output using matplotlib.
No seaborn dependency; all visualizations built from numpy, pandas, and matplotlib.
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
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch
try:
    from scipy import stats as _scipy_stats
    _HAVE_SCIPY = True
except ImportError:
    _scipy_stats = None
    _HAVE_SCIPY = False


def _linregress(x, y):
    """Pure-numpy linear regression. Returns (slope, intercept, r, p, se)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    xm, ym = x.mean(), y.mean()
    ssxx = np.sum((x - xm) ** 2)
    ssxy = np.sum((x - xm) * (y - ym))
    slope = ssxy / (ssxx + 1e-300)
    intercept = ym - slope * xm
    y_pred = slope * x + intercept
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - ym) ** 2)
    r_value = np.sqrt(max(0.0, 1.0 - ss_res / (ss_tot + 1e-300))) * np.sign(slope)
    mse = ss_res / max(1, n - 2)
    se = np.sqrt(mse / (ssxx + 1e-300))
    # t-based p-value approximation
    t_stat = slope / (se + 1e-300)
    p_value = _t_sf_approx(abs(t_stat), n - 2) * 2.0
    return slope, intercept, r_value, min(p_value, 1.0), se


def _t_sf_approx(t, df):
    """Survival function of t-distribution (rough Wilson-Hilferty approx)."""
    import math
    if df <= 0:
        return 0.5
    x = df / (df + t * t)
    z = (x ** (1.0 / 3.0) - (1.0 - 2.0 / (9.0 * df))) / np.sqrt(2.0 / (9.0 * df))
    return 0.5 * (1.0 - math.erf(z / math.sqrt(2.0)))


def _t_ppf_approx(p, df):
    """Approximate t-distribution quantile via normal approximation."""
    import math
    # Normal quantile via Newton on erfinv approximation
    # Use iterative approach for |p-0.5| < 0.5
    if df > 200:
        z = _norm_ppf(p)
        return z
    # Cornish-Fisher expansion
    z = _norm_ppf(p)
    g1 = (z ** 3 + z) / (4.0 * df)
    g2 = (5 * z ** 5 + 16 * z ** 3 + 3 * z) / (96.0 * df ** 2)
    return z + g1 + g2


def _norm_ppf(p):
    """Rational approximation for the normal quantile (Beasley-Springer-Moro)."""
    import math
    p = float(p)
    if p <= 0.0:
        return -float('inf')
    if p >= 1.0:
        return float('inf')
    if p < 0.5:
        t = math.sqrt(-2.0 * math.log(p))
    else:
        t = math.sqrt(-2.0 * math.log(1.0 - p))
    # Rational approx
    c = [2.515517, 0.802853, 0.010328]
    d = [1.432788, 0.189269, 0.001308]
    num = c[0] + c[1] * t + c[2] * t * t
    den = 1.0 + d[0] * t + d[1] * t * t + d[2] * t ** 3
    z = t - num / den
    return z if p >= 0.5 else -z


class _FallbackStats:
    @staticmethod
    def linregress(x, y):
        return _linregress(x, y)

    class t:
        @staticmethod
        def ppf(p, df):
            return _t_ppf_approx(p, df)

    @staticmethod
    def pearsonr(x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        xm = x - x.mean(); ym = y - y.mean()
        r = np.sum(xm * ym) / (np.sqrt(np.sum(xm**2) * np.sum(ym**2)) + 1e-300)
        r = max(-1.0, min(1.0, r))
        n = len(x)
        t_stat = r * np.sqrt(n - 2) / (np.sqrt(1 - r**2) + 1e-300)
        p = _t_sf_approx(abs(t_stat), n - 2) * 2.0
        return r, min(p, 1.0)

    @staticmethod
    def spearmanr(x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        rx = np.argsort(np.argsort(x)).astype(float)
        ry = np.argsort(np.argsort(y)).astype(float)
        return _FallbackStats.pearsonr(rx, ry)


stats = _scipy_stats if _HAVE_SCIPY else _FallbackStats()

warnings.filterwarnings('ignore', category=UserWarning)


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def load_data(input_file):
    """Load data from CSV/TSV file."""
    if input_file.endswith(('.csv', '.tsv', '.txt')):
        sep = '\t' if input_file.endswith(('.tsv', '.txt')) else ','
        return pd.read_csv(input_file, sep=sep, index_col=False)
    else:
        return pd.read_csv(input_file, sep=None, engine='python', index_col=False)


def parse_color_list(color_str):
    """Parse comma-separated hex colors."""
    return [c.strip() for c in color_str.split(',')]


def parse_feature_list(feature_str, file_path=None):
    """Parse feature list from comma-separated string or file."""
    if file_path and Path(file_path).exists():
        with open(file_path) as f:
            return [line.strip() for line in f if line.strip()]
    return [f.strip() for f in feature_str.split(',') if f.strip()]


def get_color_for_category(cat, unique_cats, palette='Set2', colors=None):
    """Get color for categorical value."""
    if colors:
        idx = min(list(unique_cats).index(cat) if cat in unique_cats else 0, len(colors) - 1)
        return colors[idx]

    palette_map = {
        'Set2': ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3'],
        'tab10': ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'],
        'Pastel1': ['#fbb4ae', '#b3cde3', '#ccebc5', '#decbe4', '#fed9a6', '#ffffcc', '#e5d8bd', '#fddaec', '#f2f2f2']
    }

    palette_colors = palette_map.get(palette, palette_map['Set2'])
    idx = list(unique_cats).index(cat) % len(palette_colors)
    return palette_colors[idx]


def loess_smoothing(x, y, frac=0.3, num_points=100):
    """
    Simple LOESS (Local Polynomial Regression) smoothing using tricubic kernel.

    Args:
        x, y: arrays
        frac: fraction of data for local window (default 0.3)
        num_points: number of points for output curve

    Returns:
        x_smooth, y_smooth: smoothed curve
    """
    x = np.asarray(x).flatten()
    y = np.asarray(y).flatten()

    # Sort by x
    sort_idx = np.argsort(x)
    x_sorted = x[sort_idx]
    y_sorted = y[sort_idx]

    # Generate interpolation points
    x_smooth = np.linspace(x_sorted.min(), x_sorted.max(), num_points)
    y_smooth = []

    window_size = max(3, int(len(x_sorted) * frac))

    for x_i in x_smooth:
        # Find nearest neighbors
        distances = np.abs(x_sorted - x_i)
        nearest_idx = np.argsort(distances)[:window_size]

        # Tricubic kernel weights
        d_max = distances[nearest_idx[-1]]
        if d_max > 0:
            d_normalized = distances[nearest_idx] / d_max
            weights = (1 - d_normalized**3)**3
        else:
            weights = np.ones(window_size)

        # Weighted least squares
        X_local = np.column_stack([np.ones(window_size), x_sorted[nearest_idx]])
        W = np.diag(weights)
        beta = np.linalg.lstsq(X_local.T @ W @ X_local, X_local.T @ W @ y_sorted[nearest_idx], rcond=None)[0]
        y_smooth.append(beta[0] + beta[1] * x_i)

    return x_smooth, np.array(y_smooth)


def regression_line_and_ci(x, y, alpha=0.05):
    """
    Compute regression line and 95% confidence interval band.

    Returns:
        slope, intercept, x_sorted, y_pred, y_upper, y_lower
    """
    x = np.asarray(x).flatten()
    y = np.asarray(y).flatten()

    # Linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

    # Predictions
    x_sorted = np.linspace(x.min(), x.max(), 100)
    y_pred = slope * x_sorted + intercept

    # 95% CI: t * SE
    n = len(x)
    t_val = stats.t.ppf(1 - alpha/2, n - 2)
    residuals = y - (slope * x + intercept)
    mse = np.sum(residuals**2) / (n - 2)
    se = np.sqrt(mse * (1/n + (x_sorted - x.mean())**2 / np.sum((x - x.mean())**2)))

    y_upper = y_pred + t_val * se
    y_lower = y_pred - t_val * se

    return slope, intercept, x_sorted, y_pred, y_upper, y_lower


def set_style(ax, args):
    """Apply shared style parameters to axis."""
    ax.set_facecolor('white')
    ax.grid(False)

    if args.grid != 'none':
        ax.grid(True, which='major' if args.grid == 'major' else 'both',
                alpha=args.grid_alpha, linestyle='--', linewidth=0.5)

    if args.spine_style == 'minimal':
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.setp(ax.get_xticklabels(), fontsize=args.tick_size, fontfamily=args.font_family)
    plt.setp(ax.get_yticklabels(), fontsize=args.tick_size, fontfamily=args.font_family)

    if ax.get_xlabel():
        ax.xaxis.label.set_fontsize(args.axis_label_size)
        ax.xaxis.label.set_fontfamily(args.font_family)
    if ax.get_ylabel():
        ax.yaxis.label.set_fontsize(args.axis_label_size)
        ax.yaxis.label.set_fontfamily(args.font_family)

    if ax.get_title():
        ax.title.set_fontsize(args.title_size)
        ax.title.set_fontfamily(args.font_family)


# ============================================================================
# PLOT FUNCTIONS
# ============================================================================

def plot_scatter(data, args):
    """
    Scatter plot with optional regression line, highlights, and marginal histograms.
    """
    print(f"Creating scatter plot: {args.x_col} vs {args.y_col}")

    # Extract data
    x = data[args.x_col].values.astype(float)
    y = data[args.y_col].values.astype(float)

    # Remove NaN
    valid_idx = ~(np.isnan(x) | np.isnan(y))
    x = x[valid_idx]
    y = y[valid_idx]

    # Log transform if requested
    if args.log_x:
        x = np.log10(x)
    if args.log_y:
        y = np.log10(y)

    # Setup figure
    if args.marginal_hist:
        fig = plt.figure(figsize=(args.fig_width, args.fig_height))
        gs = fig.add_gridspec(2, 2, width_ratios=[3, 1], height_ratios=[1, 3],
                              hspace=0.05, wspace=0.05)
        ax = fig.add_subplot(gs[1, 0])
        ax_hx = fig.add_subplot(gs[0, 0], sharex=ax)
        ax_hy = fig.add_subplot(gs[1, 1], sharey=ax)
        ax_hx.set_xticklabels([])
        ax_hy.set_yticklabels([])
    else:
        fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height), dpi=args.dpi)

    # Determine colors
    colors = None
    colorbar_obj = None
    if args.color_col:
        color_data = data.loc[valid_idx, args.color_col].values
        if pd.api.types.is_numeric_dtype(color_data):
            # Numeric colormap
            norm = Normalize(vmin=np.nanmin(color_data), vmax=np.nanmax(color_data))
            cmap = plt.get_cmap(args.colormap)
            colors = cmap(norm(color_data))
            if args.show_colorbar:
                sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
                sm.set_array([])
                cbar = plt.colorbar(sm, ax=ax, pad=0.02)
                cbar.set_label(args.color_col, fontsize=args.axis_label_size, fontfamily=args.font_family)
        else:
            # Categorical colors
            unique_cats = color_data.unique()
            colors = [get_color_for_category(cat, unique_cats, args.palette) for cat in color_data]

    # Determine sizes
    sizes = np.ones_like(x) * args.point_size
    if args.size_col and args.size_col in data.columns:
        size_data = data.loc[valid_idx, args.size_col].values.astype(float)
        sizes = args.point_size + (size_data - size_data.min()) / (size_data.max() - size_data.min() + 1e-10) * 50

    # Plot points
    if colors is None:
        ax.scatter(x, y, s=sizes, alpha=args.alpha, color=args.color_ns, edgecolors='none')
    else:
        ax.scatter(x, y, s=sizes, alpha=args.alpha, c=colors, edgecolors='none')

    # Highlight specific features
    if args.highlight:
        highlight_list = parse_feature_list(args.highlight)
        if args.label_col and args.label_col in data.columns:
            label_data = data.loc[valid_idx, args.label_col].values
            for hl in highlight_list:
                hl_idx = np.where(label_data == hl)[0]
                if len(hl_idx) > 0:
                    ax.scatter(x[hl_idx], y[hl_idx], s=sizes[hl_idx]*1.5,
                              color=args.highlight_color, edgecolors='black', linewidths=1.5, zorder=5)

    # Regression line
    if args.show_regression:
        slope, intercept, x_reg, y_pred, y_upper, y_lower = regression_line_and_ci(x, y)
        ax.plot(x_reg, y_pred, color=args.regression_color, linewidth=2, zorder=4)
        if args.regression_ci:
            ax.fill_between(x_reg, y_lower, y_upper, alpha=0.2, color=args.regression_color, zorder=3)

    # Identity line
    if args.identity_line:
        lims = [np.min([ax.get_xlim(), ax.get_ylim()]), np.max([ax.get_xlim(), ax.get_ylim()])]
        ax.plot(lims, lims, 'k--', alpha=0.3, zorder=2)

    # Labels for top N points
    if args.label_top_n > 0 and args.label_col and args.label_col in data.columns:
        label_data = data.loc[valid_idx, args.label_col].values
        distances = np.sqrt((x - x.mean())**2 + (y - y.mean())**2)
        top_idx = np.argsort(distances)[-args.label_top_n:]
        for idx in top_idx:
            ax.annotate(label_data[idx], (x[idx], y[idx]), fontsize=8, alpha=0.7,
                       xytext=(5, 5), textcoords='offset points')

    # Labels
    ax.set_xlabel(args.xlabel or args.x_col, fontsize=args.axis_label_size, fontfamily=args.font_family)
    ax.set_ylabel(args.ylabel or args.y_col, fontsize=args.axis_label_size, fontfamily=args.font_family)
    if args.title:
        ax.set_title(args.title, fontsize=args.title_size, fontfamily=args.font_family, pad=20)

    # Limits
    if args.xlim:
        ax.set_xlim(args.xlim)
    if args.ylim:
        ax.set_ylim(args.ylim)

    # Marginal histograms
    if args.marginal_hist:
        ax_hx.hist(x, bins=30, color=args.regression_color, alpha=0.5, edgecolor='none')
        ax_hy.hist(y, bins=30, orientation='horizontal', color=args.regression_color, alpha=0.5, edgecolor='none')
        ax_hx.set_ylabel('Count', fontsize=args.axis_label_size)
        ax_hy.set_xlabel('Count', fontsize=args.axis_label_size)
        set_style(ax_hx, args)
        set_style(ax_hy, args)

    set_style(ax, args)

    print(f"  Points: {len(x)}")
    if args.show_regression:
        print(f"  Regression: y = {slope:.4f}*x + {intercept:.4f}")

    return fig


def plot_bar(data, args):
    """
    Bar plot with optional error bars, grouping, stacking, and statistical annotations.
    """
    print(f"Creating bar plot: {args.x_col} vs {args.y_col}")

    # Prepare data
    plot_data = data[[args.x_col, args.y_col]].copy()
    if args.color_col:
        plot_data[args.color_col] = data[args.color_col]
    if args.error_col:
        plot_data[args.error_col] = data[args.error_col]

    plot_data = plot_data.dropna(subset=[args.x_col, args.y_col])

    # Aggregate if needed
    if args.error_col and args.error_col in plot_data.columns:
        # Pre-computed errors
        if args.color_col:
            agg_data = plot_data.groupby([args.x_col, args.color_col]).agg({
                args.y_col: 'mean',
                args.error_col: 'first'
            }).reset_index()
        else:
            agg_data = plot_data.groupby(args.x_col).agg({
                args.y_col: 'mean',
                args.error_col: 'first'
            }).reset_index()
    else:
        # Compute from raw data
        if args.color_col:
            agg_data = plot_data.groupby([args.x_col, args.color_col])[args.y_col].agg(['mean', 'sem', 'std']).reset_index()
            agg_data.columns = [args.x_col, args.color_col, 'mean', 'sem', 'std']
            if args.error_type == 'sem':
                agg_data[args.error_col] = agg_data['sem']
            elif args.error_type == 'sd':
                agg_data[args.error_col] = agg_data['std']
            else:  # ci95
                agg_data[args.error_col] = 1.96 * agg_data['sem']
        else:
            agg_data = plot_data.groupby(args.x_col)[args.y_col].agg(['mean', 'sem', 'std']).reset_index()
            agg_data.columns = [args.x_col, 'mean', 'sem', 'std']
            if args.error_type == 'sem':
                agg_data[args.error_col] = agg_data['sem']
            elif args.error_type == 'sd':
                agg_data[args.error_col] = agg_data['std']
            else:  # ci95
                agg_data[args.error_col] = 1.96 * agg_data['sem']
            # Restore y_col name for downstream access
            agg_data[args.y_col] = agg_data['mean']

    # Ensure grouped case also has y_col accessible
    if args.color_col and 'mean' in agg_data.columns and args.y_col not in agg_data.columns:
        agg_data[args.y_col] = agg_data['mean']

    # Sort if requested
    if args.sort_by == 'value':
        agg_data = agg_data.sort_values(args.y_col, ascending=False)
    elif args.sort_by == 'name':
        agg_data = agg_data.sort_values(args.x_col)

    # Figure setup
    fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height), dpi=args.dpi)

    if args.color_col:
        # Grouped or stacked bars
        groups = agg_data[args.color_col].unique()
        if args.group_order:
            groups = [g for g in args.group_order.split(',') if g in groups]

        x_cats = agg_data[args.x_col].unique()
        if args.group_order:
            x_cats = sorted(x_cats)

        x_pos = np.arange(len(x_cats))
        bar_width = args.bar_width / len(groups)

        for i, group in enumerate(groups):
            group_data = agg_data[agg_data[args.color_col] == group]
            group_data = group_data.set_index(args.x_col).reindex(x_cats).reset_index()

            color = get_color_for_category(group, groups, args.palette,
                                          parse_color_list(args.bar_colors) if args.bar_colors else None)

            offset = (i - len(groups)/2 + 0.5) * bar_width
            bars = ax.bar(x_pos + offset, group_data[args.y_col].values, bar_width,
                         label=group, color=color, alpha=0.8, edgecolor='black', linewidth=0.5)

            if args.error_col in group_data.columns:
                ax.errorbar(x_pos + offset, group_data[args.y_col].values,
                           yerr=group_data[args.error_col].values,
                           fmt='none', color='black', capsize=args.capsize, linewidth=1)

            if args.show_values:
                for j, (pos, val) in enumerate(zip(x_pos + offset, group_data[args.y_col].values)):
                    ax.text(pos, val, f'{val:{args.value_fmt}}', ha='center', va='bottom', fontsize=8)

        ax.set_xticks(x_pos)
        ax.set_xticklabels(x_cats, rotation=45, ha='right')
        ax.legend(title=args.color_col, fontsize=args.legend_size, title_fontsize=args.axis_label_size,
                 loc='best', frameon=True, edgecolor='black')
    else:
        # Single bar plot
        x_cats = agg_data[args.x_col].values
        x_pos = np.arange(len(x_cats))

        bars = ax.bar(x_pos, agg_data[args.y_col].values, args.bar_width,
                     color=args.color_ns, alpha=0.7, edgecolor='black', linewidth=0.5)

        if args.error_col in agg_data.columns:
            ax.errorbar(x_pos, agg_data[args.y_col].values,
                       yerr=agg_data[args.error_col].values,
                       fmt='none', color='black', capsize=args.capsize, linewidth=1)

        if args.show_values:
            for pos, val in zip(x_pos, agg_data[args.y_col].values):
                ax.text(pos, val, f'{val:{args.value_fmt}}', ha='center', va='bottom', fontsize=8)

        ax.set_xticks(x_pos)
        ax.set_xticklabels(x_cats, rotation=45, ha='right')

    # Overlay points if requested
    if args.show_points and args.color_col not in data.columns:
        for i, cat in enumerate(agg_data[args.x_col].unique()):
            cat_data = plot_data[plot_data[args.x_col] == cat][args.y_col].values
            jitter = np.random.normal(i, args.point_jitter, len(cat_data))
            ax.scatter(jitter, cat_data, alpha=0.4, s=20, color='black', zorder=3)

    ax.set_xlabel(args.xlabel or args.x_col, fontsize=args.axis_label_size, fontfamily=args.font_family)
    ax.set_ylabel(args.ylabel or args.y_col, fontsize=args.axis_label_size, fontfamily=args.font_family)
    if args.title:
        ax.set_title(args.title, fontsize=args.title_size, fontfamily=args.font_family, pad=20)

    if args.ylim:
        ax.set_ylim(args.ylim)

    set_style(ax, args)
    fig.tight_layout()

    print(f"  Groups: {len(agg_data)}")

    return fig


def plot_ma(data, args):
    """
    MA plot: M (log2FC) vs A (mean expression).
    """
    print(f"Creating MA plot: {args.x_col} vs {args.y_col}")

    # Extract data
    a = data[args.x_col].values.astype(float)
    m = data[args.y_col].values.astype(float)

    # Auto-log A if needed
    if np.max(a) > 100:
        a = np.log2(a + 1)

    p_values = data[args.p_col].values.astype(float) if args.p_col and args.p_col in data.columns else np.ones_like(m)

    # Remove NaN
    valid_idx = ~(np.isnan(a) | np.isnan(m))
    a = a[valid_idx]
    m = m[valid_idx]
    p_values = p_values[valid_idx]

    # Classify points
    sig_idx = (p_values < args.p_cutoff) & (np.abs(m) > args.fc_cutoff)
    up_idx = sig_idx & (m > 0)
    down_idx = sig_idx & (m < 0)
    ns_idx = ~sig_idx

    # Figure
    fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height), dpi=args.dpi)

    # Plot points
    ax.scatter(a[ns_idx], m[ns_idx], s=args.point_size, alpha=args.alpha,
              color=args.color_ns, edgecolors='none', label='Not significant')
    ax.scatter(a[up_idx], m[up_idx], s=args.point_size, alpha=args.alpha,
              color=args.color_up, edgecolors='none', label=f'Up (p<{args.p_cutoff})')
    ax.scatter(a[down_idx], m[down_idx], s=args.point_size, alpha=args.alpha,
              color=args.color_down, edgecolors='none', label=f'Down (p<{args.p_cutoff})')

    # LOESS line
    if args.loess_line:
        x_smooth, y_smooth = loess_smoothing(a, m, frac=0.3)
        ax.plot(x_smooth, y_smooth, color=args.loess_color, linewidth=2, zorder=4)

    # Zero line
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1, alpha=0.5)

    # Significance lines
    ax.axhline(y=args.fc_cutoff, color='red', linestyle='--', linewidth=1, alpha=0.3)
    ax.axhline(y=-args.fc_cutoff, color='red', linestyle='--', linewidth=1, alpha=0.3)

    # Annotations
    if args.annotate_top_n > 0 and args.feature_col and args.feature_col in data.columns:
        feature_data = data.loc[valid_idx, args.feature_col].values
        sig_features = np.where(sig_idx)[0]
        if len(sig_features) > 0:
            # Top by p-value and absolute FC
            scores = -np.log10(p_values[sig_features] + 1e-300) * np.abs(m[sig_features])
            top_idx = sig_features[np.argsort(scores)[-args.annotate_top_n:]]

            for idx in top_idx:
                ax.annotate(feature_data[idx], (a[idx], m[idx]), fontsize=8, alpha=0.8,
                           xytext=(5, 5), textcoords='offset points')

    ax.set_xlabel(args.xlabel, fontsize=args.axis_label_size, fontfamily=args.font_family)
    ax.set_ylabel(args.ylabel, fontsize=args.axis_label_size, fontfamily=args.font_family)
    if args.title:
        ax.set_title(args.title, fontsize=args.title_size, fontfamily=args.font_family, pad=20)

    ax.legend(fontsize=args.legend_size, loc='best', frameon=True, edgecolor='black')
    set_style(ax, args)
    fig.tight_layout()

    print(f"  Points: {len(a)}")
    print(f"  Significant: {np.sum(sig_idx)} (up: {np.sum(up_idx)}, down: {np.sum(down_idx)})")

    return fig


def plot_corrmat(data, args):
    """
    Correlation matrix heatmap with optional clustering.
    """
    print(f"Creating correlation matrix plot")

    # Extract numeric columns
    numeric_data = data.select_dtypes(include=[np.number])
    if args.index_col and args.index_col in data.columns:
        numeric_data.index = data[args.index_col]

    # Compute correlation
    if args.method == 'spearman':
        corr = numeric_data.corr(method='spearman')
    else:
        corr = numeric_data.corr(method='pearson')

    # Cluster if requested
    if args.cluster:
        try:
            from scipy.cluster.hierarchy import dendrogram, linkage
            from scipy.spatial.distance import pdist, squareform

            # Distance matrix
            dist = 1 - corr.abs()
            Z = linkage(pdist(dist, metric='euclidean'), method='average')
            dendro = dendrogram(Z, no_plot=True)
            order = dendro['leaves']
            corr = corr.iloc[order, order]
        except:
            # Fallback: sort by first PC
            pass

    # Figure
    fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height), dpi=args.dpi)

    # Heatmap
    im = ax.imshow(corr.values, cmap=args.cmap, vmin=args.vmin, vmax=args.vmax, aspect='auto')

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label('Correlation', fontsize=args.axis_label_size)

    # Labels
    ax.set_xticks(np.arange(len(corr.columns)))
    ax.set_yticks(np.arange(len(corr.index)))
    ax.set_xticklabels(corr.columns, rotation=args.label_rotation, ha='right', fontsize=args.tick_size)
    ax.set_yticklabels(corr.index, fontsize=args.tick_size)

    # Values
    if args.show_values:
        for i in range(len(corr)):
            for j in range(len(corr.columns)):
                if not (args.mask_upper and i < j):
                    text = ax.text(j, i, f'{corr.iloc[i, j]:{args.value_fmt}}',
                                  ha='center', va='center', color='black', fontsize=args.value_fontsize)

    # Mask upper triangle
    if args.mask_upper:
        for i in range(len(corr)):
            for j in range(i + 1, len(corr.columns)):
                ax.add_patch(mpatches.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                               fill=True, facecolor='white', edgecolor='none'))

    if args.title:
        ax.set_title(args.title, fontsize=args.title_size, fontfamily=args.font_family, pad=20)

    fig.tight_layout()

    print(f"  Samples: {len(corr)}")
    print(f"  Method: {args.method}")

    return fig


def plot_bubble(data, args):
    """
    Bubble plot with size and optional color encoding.
    """
    print(f"Creating bubble plot: {args.x_col} vs {args.y_col}, size by {args.size_col}")

    # Extract data
    x = data[args.x_col].values.astype(float)
    y = data[args.y_col].values.astype(float)
    size_data = data[args.size_col].values.astype(float)

    # Normalize sizes
    size_norm = (size_data - size_data.min()) / (size_data.max() - size_data.min() + 1e-10)
    sizes = args.min_size + size_norm * (args.max_size - args.min_size)

    # Color data
    colors = None
    if args.color_col and args.color_col in data.columns:
        color_data = data[args.color_col].values
        if pd.api.types.is_numeric_dtype(color_data):
            norm = Normalize(vmin=np.nanmin(color_data), vmax=np.nanmax(color_data))
            cmap = plt.get_cmap(args.colormap)
            colors = cmap(norm(color_data))
        else:
            unique_cats = np.unique(color_data)
            colors = [get_color_for_category(cat, unique_cats, args.palette) for cat in color_data]

    # Remove NaN
    valid_idx = ~(np.isnan(x) | np.isnan(y) | np.isnan(sizes))
    x = x[valid_idx]
    y = y[valid_idx]
    sizes = sizes[valid_idx]
    if colors is not None:
        colors = colors[valid_idx]

    # Figure
    fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height), dpi=args.dpi)

    # Plot bubbles
    if colors is None:
        ax.scatter(x, y, s=sizes, alpha=args.alpha, color='#999999', edgecolors='black', linewidths=0.5)
    else:
        ax.scatter(x, y, s=sizes, alpha=args.alpha, c=colors, edgecolors='black', linewidths=0.5)

    # Labels
    if args.label_all and args.label_col and args.label_col in data.columns:
        label_data = data.loc[valid_idx, args.label_col].values
        for i, (xi, yi, label) in enumerate(zip(x, y, label_data)):
            ax.annotate(label, (xi, yi), fontsize=7, alpha=0.7,
                       xytext=(3, 3), textcoords='offset points')

    ax.set_xlabel(args.xlabel or args.x_col, fontsize=args.axis_label_size, fontfamily=args.font_family)
    ax.set_ylabel(args.ylabel or args.y_col, fontsize=args.axis_label_size, fontfamily=args.font_family)
    if args.title:
        ax.set_title(args.title, fontsize=args.title_size, fontfamily=args.font_family, pad=20)

    # Size legend
    if args.show_size_legend:
        legend_sizes = [args.min_size, (args.min_size + args.max_size) / 2, args.max_size]
        legend_labels = [f'{size_data.min():.1f}', f'{(size_data.min() + size_data.max())/2:.1f}', f'{size_data.max():.1f}']
        legend_bubbles = [ax.scatter([], [], s=s, alpha=args.alpha, color='#999999', edgecolors='black', linewidths=0.5)
                         for s in legend_sizes]
        leg1 = ax.legend(legend_bubbles, legend_labels, scatterpoints=1, title=args.size_col,
                        fontsize=args.legend_size, loc='upper left', frameon=True, edgecolor='black')

    set_style(ax, args)
    fig.tight_layout()

    print(f"  Bubbles: {len(x)}")

    return fig


# ============================================================================
# MAIN FUNCTION
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Generate publication-quality scatter, bar, MA, correlation, and bubble plots',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Input/output
    parser.add_argument('input', help='Input CSV/TSV file')
    parser.add_argument('--plot-type', required=True, choices=['scatter', 'bar', 'ma', 'corrmat', 'bubble'],
                       help='Plot type to generate')
    parser.add_argument('--output', required=True, help='Output file (PNG/SVG/PDF)')
    parser.add_argument('--output-svg', help='Additionally save as SVG')
    parser.add_argument('--output-table', help='Save processed data as TSV')

    # Scatter plot arguments
    parser.add_argument('--x-col', help='X-axis column (required for scatter, bar, bubble)')
    parser.add_argument('--y-col', help='Y-axis column (required for scatter, bar, bubble, ma)')
    parser.add_argument('--color-col', help='Column to color points by')
    parser.add_argument('--size-col', help='Column to scale point size by')
    parser.add_argument('--label-col', help='Column for point labels')
    parser.add_argument('--label-top-n', type=int, default=0, help='Label top N points by distance from origin')
    parser.add_argument('--highlight', help='Comma-separated feature names to highlight')
    parser.add_argument('--highlight-color', default='#D55E00', help='Highlight color')
    parser.add_argument('--corr-method', choices=['pearson', 'spearman', 'none'], default='pearson',
                       help='Correlation method for scatter')
    parser.add_argument('--show-regression', type=bool, default=True, help='Show regression line')
    parser.add_argument('--regression-ci', type=bool, default=True, help='Show 95% confidence interval')
    parser.add_argument('--regression-color', default='#2166AC', help='Regression line color')
    parser.add_argument('--point-size', type=float, default=25, help='Point size')
    parser.add_argument('--alpha', type=float, default=0.7, help='Point alpha transparency')
    parser.add_argument('--color-ns', default='#999999', help='Color for non-significant points')
    parser.add_argument('--colormap', default='viridis', help='Colormap for numeric colors')
    parser.add_argument('--show-colorbar', type=bool, default=True, help='Show colorbar')
    parser.add_argument('--identity-line', type=bool, default=False, help='Draw y=x diagonal line')
    parser.add_argument('--xlabel', help='X-axis label')
    parser.add_argument('--ylabel', help='Y-axis label')
    parser.add_argument('--title', help='Plot title')
    parser.add_argument('--xlim', type=float, nargs=2, help='X-axis limits')
    parser.add_argument('--ylim', type=float, nargs=2, help='Y-axis limits')
    parser.add_argument('--log-x', type=bool, default=False, help='Log10-scale x axis')
    parser.add_argument('--log-y', type=bool, default=False, help='Log10-scale y axis')
    parser.add_argument('--marginal-hist', type=bool, default=False, help='Show marginal histograms')

    # Bar plot arguments
    parser.add_argument('--error-col', help='Column with pre-computed error bars')
    parser.add_argument('--error-type', choices=['sem', 'sd', 'ci95'], default='sem',
                       help='Error type if computing from data')
    parser.add_argument('--group-order', help='Comma-separated x-axis order')
    parser.add_argument('--color-order', help='Comma-separated color group order')
    parser.add_argument('--palette', choices=['Set2', 'tab10', 'Pastel1'], default='Set2',
                       help='Color palette')
    parser.add_argument('--bar-colors', help='Comma-separated hex colors for bars')
    parser.add_argument('--orientation', choices=['vertical', 'horizontal'], default='vertical',
                       help='Bar orientation')
    parser.add_argument('--bar-width', type=float, default=0.7, help='Bar width')
    parser.add_argument('--show-points', type=bool, default=False, help='Overlay individual data points')
    parser.add_argument('--point-jitter', type=float, default=0.05, help='Point jitter amount')
    parser.add_argument('--capsize', type=float, default=4, help='Error bar cap size')
    parser.add_argument('--show-values', type=bool, default=False, help='Show bar height as text')
    parser.add_argument('--value-fmt', default='.2f', help='Format for displayed values')
    parser.add_argument('--sort-by', choices=['value', 'name', 'none'], default='none',
                       help='Sort bars by value or name')
    parser.add_argument('--stacked', type=bool, default=False, help='Stacked bar chart')
    parser.add_argument('--stats', choices=['none', 'pairwise'], default='none',
                       help='Statistical annotations')

    # MA plot arguments
    parser.add_argument('--p-col', help='P-value column for MA plot')
    parser.add_argument('--feature-col', help='Gene/feature name column')
    parser.add_argument('--p-cutoff', type=float, default=0.05, help='P-value cutoff')
    parser.add_argument('--fc-cutoff', type=float, default=1.0, help='Fold-change cutoff (log2)')
    parser.add_argument('--color-up', default='#D55E00', help='Color for upregulated')
    parser.add_argument('--color-down', default='#0072B2', help='Color for downregulated')
    parser.add_argument('--annotate-top-n', type=int, default=10, help='Annotate top N DE genes')
    parser.add_argument('--loess-line', type=bool, default=True, help='Draw LOESS smoothing line')
    parser.add_argument('--loess-color', default='black', help='LOESS line color')

    # Correlation matrix arguments
    parser.add_argument('--index-col', help='Row label column for corrmat')
    parser.add_argument('--method', choices=['pearson', 'spearman'], default='pearson',
                       help='Correlation method')
    parser.add_argument('--cluster', type=bool, default=True, help='Cluster correlation matrix')
    parser.add_argument('--cmap', default='RdBu_r', help='Colormap for corrmat')
    parser.add_argument('--vmin', type=float, default=-1, help='Colormap min value')
    parser.add_argument('--vmax', type=float, default=1, help='Colormap max value')
    parser.add_argument('--mask-upper', type=bool, default=False, help='Mask upper triangle')
    parser.add_argument('--label-rotation', type=float, default=45, help='Label rotation angle')

    # Bubble plot arguments
    parser.add_argument('--min-size', type=float, default=50, help='Min bubble area')
    parser.add_argument('--max-size', type=float, default=2000, help='Max bubble area')
    parser.add_argument('--show-size-legend', type=bool, default=True, help='Show bubble size legend')
    parser.add_argument('--show-color-legend', type=bool, default=True, help='Show color legend')
    parser.add_argument('--label-all', type=bool, default=False, help='Label all bubbles')

    # Style arguments
    parser.add_argument('--font-family', default='Arial', help='Font family')
    parser.add_argument('--base-fontsize', type=float, default=11, help='Base font size')
    parser.add_argument('--axis-label-size', type=float, default=12, help='Axis label size')
    parser.add_argument('--tick-size', type=float, default=10, help='Tick label size')
    parser.add_argument('--legend-size', type=float, default=10, help='Legend font size')
    parser.add_argument('--title-size', type=float, default=13, help='Title font size')
    parser.add_argument('--fig-width', type=float, default=8, help='Figure width (inches)')
    parser.add_argument('--fig-height', type=float, default=6, help='Figure height (inches)')
    parser.add_argument('--dpi', type=int, default=300, help='Output DPI')
    parser.add_argument('--spine-style', choices=['all', 'minimal'], default='minimal',
                       help='Spine visibility')
    parser.add_argument('--grid', choices=['none', 'major', 'both'], default='none',
                       help='Grid display')
    parser.add_argument('--grid-alpha', type=float, default=0.3, help='Grid transparency')

    args = parser.parse_args()
    init_style(
        font_family=getattr(args, 'font_family', None),
        font_size=getattr(args, 'base_fontsize', None),
    )

    # Load data
    try:
        data = load_data(args.input)
        print(f"Loaded {len(data)} rows × {len(data.columns)} columns from {args.input}")
    except Exception as e:
        print(f"Error loading input file: {e}", file=sys.stderr)
        sys.exit(1)

    # Generate plot
    try:
        if args.plot_type == 'scatter':
            fig = plot_scatter(data, args)
        elif args.plot_type == 'bar':
            fig = plot_bar(data, args)
        elif args.plot_type == 'ma':
            fig = plot_ma(data, args)
        elif args.plot_type == 'corrmat':
            fig = plot_corrmat(data, args)
        elif args.plot_type == 'bubble':
            fig = plot_bubble(data, args)
        else:
            print(f"Unknown plot type: {args.plot_type}", file=sys.stderr)
            sys.exit(1)

        # Save figure
        fig.savefig(args.output, dpi=args.dpi, bbox_inches='tight', facecolor='white')
        print(f"Plot saved to {args.output}")

        if args.output_svg:
            fig.savefig(args.output_svg, format='svg', bbox_inches='tight', facecolor='white')
            print(f"SVG saved to {args.output_svg}")

        if args.output_table:
            data.to_csv(args.output_table, sep='\t', index=False)
            print(f"Data saved to {args.output_table}")

        plt.close(fig)

    except Exception as e:
        print(f"Error generating plot: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
