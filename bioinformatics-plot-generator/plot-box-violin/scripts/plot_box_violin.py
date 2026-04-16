#!/usr/bin/env python3
"""
Plot Box, Violin, and Raincloud plots with statistical annotations.
Publication-quality visualization for comparing distributions across groups.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), '..', '..', '..', '_shared'))
from plot_style import init_style
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle

try:
    from scipy import stats as _scipy_stats
    _HAVE_SCIPY = True
except ImportError:
    _scipy_stats = None
    _HAVE_SCIPY = False


class _FallbackStats:
    """Pure-numpy fallbacks for scipy.stats functions."""

    @staticmethod
    def mannwhitneyu(x, y):
        """Mann-Whitney U test (two-sided approximation)."""
        x, y = np.asarray(x, float), np.asarray(y, float)
        nx, ny = len(x), len(y)
        combined = np.concatenate([x, y])
        ranks = np.argsort(np.argsort(combined)) + 1.0
        u1 = np.sum(ranks[:nx]) - nx * (nx + 1) / 2.0
        u2 = nx * ny - u1
        u = min(u1, u2)
        # Normal approximation
        mu = nx * ny / 2.0
        sigma = np.sqrt(nx * ny * (nx + ny + 1) / 12.0)
        if sigma == 0:
            return u, 1.0
        z = (u - mu) / sigma
        # two-sided p-value via erfc
        p = 2.0 * (1.0 - 0.5 * (1.0 + np.sign(z) * (1.0 - np.exp(-2.0 / np.pi * z ** 2))))
        p = min(max(p, 0.0), 1.0)
        return u, p

    @staticmethod
    def ttest_ind(x, y):
        """Welch's t-test."""
        x, y = np.asarray(x, float), np.asarray(y, float)
        nx, ny = len(x), len(y)
        mx, my = np.mean(x), np.mean(y)
        vx, vy = np.var(x, ddof=1), np.var(y, ddof=1)
        t = (mx - my) / np.sqrt(vx / nx + vy / ny + 1e-300)
        df = (vx / nx + vy / ny) ** 2 / (
            (vx / nx) ** 2 / (nx - 1) + (vy / ny) ** 2 / (ny - 1) + 1e-300
        )
        # p-value approximation using incomplete beta
        x_beta = df / (df + t ** 2)
        # Simple approximation: use erfc-based normal CDF for large df
        p = min(1.0, 2.0 * (1.0 - _norm_cdf(abs(t))))
        return t, p

    @staticmethod
    def kruskal(*groups):
        """Kruskal-Wallis H-test approximation."""
        all_data = np.concatenate([np.asarray(g, float) for g in groups])
        n = len(all_data)
        ranks = np.argsort(np.argsort(all_data)) + 1.0
        h = 12.0 / (n * (n + 1))
        offset = 0
        h_stat = 0.0
        for g in groups:
            ng = len(g)
            rg = np.sum(ranks[offset:offset + ng])
            h_stat += rg ** 2 / ng
            offset += ng
        h_stat *= h
        h_stat -= 3 * (n + 1)
        k = len(groups)
        # chi2 p-value approximation (df = k-1)
        p = _chi2_sf(h_stat, k - 1)
        return h_stat, p


def _norm_cdf(z):
    """Standard normal CDF via erf."""
    import math
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def _chi2_sf(x, df):
    """Chi-squared survival function approximation."""
    if x <= 0:
        return 1.0
    # Regularised incomplete gamma via series (rough approximation)
    # Use Wilson-Hilferty normal approximation
    if df == 0:
        return 0.0
    z = (x / df) ** (1.0 / 3.0) - (1.0 - 2.0 / (9.0 * df))
    z /= np.sqrt(2.0 / (9.0 * df))
    return 1.0 - _norm_cdf(z)


# Use scipy if available, else fallback
stats = _scipy_stats if _HAVE_SCIPY else _FallbackStats()


class PlotGenerator:
    """Generate publication-quality distribution plots with statistics."""

    def __init__(self, args):
        """Initialize plot generator with command-line arguments."""
        self.args = args
        self.data = None
        self.groups = None
        self.values = None
        self.group_names = None
        self.fig = None
        self.ax = None
        self.group_colors = None
        self.stats_results = {}

    def load_data(self):
        """Load and validate input data.

        Supports two input shapes:

        1. LONG format (original) — specify --value-col and --group-col.
           One row per observation; one numeric value column, one group column.

        2. WIDE format — each column is a group; each row is one observation.
           Activated by:
             a) --wide-cols "colA,colB,colC"  (explicit subset of columns)
             b) --wide-cols auto              (use all numeric columns)
             c) Neither --value-col nor --group-col given, but the file
                has 2+ numeric columns (auto-detected wide format).
           Different groups may have different sample sizes — NaN rows in
           each column are dropped independently.
        """
        input_path = Path(self.args.input)
        if not input_path.exists():
            raise FileNotFoundError(f"Input file not found: {self.args.input}")

        # Load TSV or CSV
        if str(input_path).endswith('.tsv'):
            self.data = pd.read_csv(input_path, sep='\t')
        else:
            self.data = pd.read_csv(input_path)

        # Decide layout: long vs wide
        has_long_args = bool(self.args.value_col) and bool(self.args.group_col)
        wide_mode = None  # None | 'explicit' | 'auto-all' | 'auto-detected'

        if self.args.wide_cols:
            wide_mode = 'auto-all' if self.args.wide_cols.strip().lower() == 'auto' else 'explicit'
        elif not has_long_args:
            # Auto-detect: if no long-format args given, fall back to treating
            # every numeric column as a group.
            numeric_cols = [c for c in self.data.columns
                            if pd.api.types.is_numeric_dtype(self.data[c])]
            if len(numeric_cols) >= 2:
                wide_mode = 'auto-detected'
                print(f"[auto] No --value-col/--group-col provided; detected "
                      f"{len(numeric_cols)} numeric columns — treating as wide format.")
            else:
                raise ValueError(
                    "No --value-col/--group-col given and only "
                    f"{len(numeric_cols)} numeric column(s) found. "
                    "Provide --value-col + --group-col (long format) or "
                    "--wide-cols 'colA,colB,colC' (wide format)."
                )

        if wide_mode:
            # ── WIDE FORMAT ─────────────────────────────────────────────────
            if wide_mode == 'explicit':
                wanted = [c.strip() for c in self.args.wide_cols.split(',') if c.strip()]
                missing = [c for c in wanted if c not in self.data.columns]
                if missing:
                    raise ValueError(f"--wide-cols: columns not found in data: {missing}")
                cols_to_use = wanted
            else:
                cols_to_use = [c for c in self.data.columns
                               if pd.api.types.is_numeric_dtype(self.data[c])]
                if len(cols_to_use) < 2:
                    raise ValueError(
                        f"--wide-cols auto: found only {len(cols_to_use)} numeric "
                        f"column(s); need at least 2 to compare."
                    )

            # Melt: stack each column into (value, group) pairs, dropping NaNs per column.
            values_list = []
            groups_list = []
            for col in cols_to_use:
                col_vals = pd.to_numeric(self.data[col], errors='coerce').dropna().values
                values_list.append(col_vals)
                groups_list.append(np.full(len(col_vals), col, dtype=object))

            self.values = np.concatenate(values_list) if values_list else np.array([])
            self.groups = np.concatenate(groups_list) if groups_list else np.array([])

            if self.args.group_order:
                self.group_names = self.args.group_order.split(',')
                # Validate: every requested group should be one of the selected cols.
                bad = [g for g in self.group_names if g not in cols_to_use]
                if bad:
                    raise ValueError(f"--group-order: unknown group(s): {bad}. "
                                     f"Available: {cols_to_use}")
            else:
                self.group_names = list(cols_to_use)  # preserve file column order

            # Back-fill value_col / group_col labels so axis labels still work.
            if not self.args.value_col:
                self.args.value_col = 'Value'
            if not self.args.group_col:
                self.args.group_col = 'Group'

            print(f"Loaded WIDE data: {len(self.values)} observations across "
                  f"{len(self.group_names)} groups ({', '.join(self.group_names)})")
            for grp in self.group_names:
                n = int(np.sum(self.groups == grp))
                print(f"  {grp}: n={n}")
            return

        # ── LONG FORMAT (original path) ─────────────────────────────────────
        if self.args.value_col not in self.data.columns:
            raise ValueError(f"Value column '{self.args.value_col}' not found in data")
        if self.args.group_col not in self.data.columns:
            raise ValueError(f"Group column '{self.args.group_col}' not found in data")

        # Extract data
        self.values = self.data[self.args.value_col].astype(float)
        self.groups = self.data[self.args.group_col].astype(str)

        # Remove NaN values
        mask = self.values.notna() & self.groups.notna()
        self.values = self.values[mask].values
        self.groups = self.groups[mask].values

        # Get unique groups
        if self.args.group_order:
            self.group_names = self.args.group_order.split(',')
        else:
            self.group_names = sorted(np.unique(self.groups))

        print(f"Loaded data: {len(self.values)} observations across {len(self.group_names)} groups")
        for grp in self.group_names:
            n = np.sum(self.groups == grp)
            print(f"  {grp}: n={n}")

    def setup_colors(self):
        """Setup color palette for groups."""
        palettes = {
            'tab10': plt.cm.tab10,
            'Set2': plt.cm.Set2,
            'Set3': plt.cm.Set3,
            'Pastel1': plt.cm.Pastel1,
        }

        if self.args.group_colors:
            # Custom hex colors
            colors = self.args.group_colors.split(',')
            self.group_colors = colors[:len(self.group_names)]
        elif self.args.palette in palettes:
            # Named palette
            n_groups = len(self.group_names)
            cmap = palettes[self.args.palette]
            self.group_colors = [cmap(i % cmap.N) for i in range(n_groups)]
        else:
            self.group_colors = ['#1f77b4'] * len(self.group_names)

    def prepare_plot_data(self):
        """Prepare data for plotting."""
        self.plot_data = {}
        for grp in self.group_names:
            self.plot_data[grp] = self.values[self.groups == grp]

    def create_figure(self):
        """Create figure with specified dimensions and style."""
        fig_width = self.args.fig_width
        fig_height = self.args.fig_height
        dpi = self.args.dpi

        self.fig, self.ax = plt.subplots(figsize=(fig_width, fig_height), dpi=dpi)

        # Set font
        plt.rcParams['font.family'] = self.args.font_family
        plt.rcParams['font.size'] = self.args.base_fontsize

        # Set background
        if self.args.background == 'dark':
            self.fig.patch.set_facecolor('#2b2b2b')
            self.ax.set_facecolor('#3b3b3b')
        elif self.args.background == 'light':
            self.ax.set_facecolor('#f5f5f5')
        else:
            self.ax.set_facecolor('white')

        # Set spines
        if self.args.spine_style == 'none':
            self.ax.spines['top'].set_visible(False)
            self.ax.spines['right'].set_visible(False)
            self.ax.spines['left'].set_visible(False)
            self.ax.spines['bottom'].set_visible(False)
        elif self.args.spine_style == 'minimal':
            self.ax.spines['top'].set_visible(False)
            self.ax.spines['right'].set_visible(False)

    def add_grid(self):
        """Add grid based on orientation and settings."""
        if self.args.grid == 'none':
            return

        kwargs = {'alpha': self.args.grid_alpha, 'linestyle': self.args.grid_linestyle}

        if self.args.orientation == 'vertical':
            if self.args.grid in ('y', 'both'):
                self.ax.yaxis.grid(True, **kwargs)
            if self.args.grid in ('x', 'both'):
                self.ax.xaxis.grid(True, **kwargs)
        else:  # horizontal
            if self.args.grid in ('x', 'both'):
                self.ax.xaxis.grid(True, **kwargs)
            if self.args.grid in ('y', 'both'):
                self.ax.yaxis.grid(True, **kwargs)

    def draw_boxplot(self, positions):
        """Draw boxplot."""
        data_list = [self.plot_data[grp] for grp in self.group_names]

        bp = self.ax.boxplot(
            data_list,
            positions=positions,
            widths=self.args.width,
            patch_artist=True,
            notch=self.args.notch,
            showfliers=self.args.show_fliers,
            vert=(self.args.orientation == 'vertical'),
        )

        # Style box elements
        for patch, color in zip(bp['boxes'], self.group_colors):
            patch.set_facecolor(color)
            patch.set_alpha(self.args.color_alpha)
            patch.set_linewidth(self.args.box_linewidth)

        for median in bp['medians']:
            median.set_color(self.args.median_color)
            median.set_linewidth(self.args.median_linewidth)

        for flier in bp['fliers']:
            flier.set_markersize(self.args.flier_size)
            flier.set_marker(self.args.flier_symbol)

        for whisker in bp['whiskers']:
            whisker.set_linewidth(self.args.box_linewidth)

        for cap in bp['caps']:
            cap.set_linewidth(self.args.box_linewidth)

        return bp

    def draw_violin(self, positions):
        """Draw violin plot."""
        data_list = [self.plot_data[grp] for grp in self.group_names]

        parts = self.ax.violinplot(
            data_list,
            positions=positions,
            widths=self.args.width,
            vert=(self.args.orientation == 'vertical'),
            showmeans=False,
            showmedians=False,
            showextrema=False,
        )

        # Style violin bodies
        for pc, color in zip(parts['bodies'], self.group_colors):
            pc.set_facecolor(color)
            pc.set_alpha(self.args.color_alpha)
            pc.set_edgecolor('black')
            pc.set_linewidth(self.args.box_linewidth)

        # Style lines
        for partname in ('cmeans', 'cmedians', 'cbars', 'cmaxes', 'cmins'):
            if partname in parts:
                vp = parts[partname]
                vp.set_edgecolor('black')
                vp.set_linewidth(self.args.box_linewidth)

        return parts

    def draw_raincloud(self, positions):
        """Draw raincloud plot (half-violin + strip + boxplot)."""
        for pos, grp, color in zip(positions, self.group_names, self.group_colors):
            data = self.plot_data[grp]

            # Draw half-violin
            self._draw_half_violin(data, pos, color)

            # Draw data points (strip plot with jitter)
            if self.args.show_points:
                self._draw_strip_points(data, pos, color)

            # Draw small boxplot in center
            self._draw_mini_boxplot(data, pos)

    def _draw_half_violin(self, data, pos, color):
        """Draw half-violin (right side by default)."""
        parts = self.ax.violinplot(
            [data],
            positions=[pos],
            widths=self.args.width,
            vert=(self.args.orientation == 'vertical'),
            showmeans=False,
            showmedians=False,
            showextrema=False,
        )

        for pc in parts['bodies']:
            pc.set_facecolor(color)
            pc.set_alpha(self.args.color_alpha)
            pc.set_edgecolor('black')
            pc.set_linewidth(self.args.box_linewidth)

    def _draw_strip_points(self, data, pos, color):
        """Draw individual data points with optional jitter."""
        if self.args.point_style == 'strip':
            jitter = np.random.normal(0, self.args.jitter, len(data))
        elif self.args.point_style == 'jitter':
            jitter = np.random.uniform(-self.args.jitter, self.args.jitter, len(data))
        else:  # swarm - simplified swarm layout
            jitter = np.random.normal(0, self.args.jitter * 0.5, len(data))

        point_positions = pos + jitter

        if self.args.point_color == 'auto':
            point_color = color
        else:
            point_color = self.args.point_color

        if self.args.orientation == 'vertical':
            self.ax.scatter(
                point_positions, data,
                s=self.args.point_size, alpha=self.args.point_alpha,
                color=point_color, edgecolor=self.args.point_edge_color,
                zorder=3
            )
        else:
            self.ax.scatter(
                data, point_positions,
                s=self.args.point_size, alpha=self.args.point_alpha,
                color=point_color, edgecolor=self.args.point_edge_color,
                zorder=3
            )

    def _draw_mini_boxplot(self, data, pos):
        """Draw small boxplot in center of raincloud."""
        bp = self.ax.boxplot(
            [data],
            positions=[pos],
            widths=self.args.width * 0.3,
            patch_artist=True,
            vert=(self.args.orientation == 'vertical'),
        )

        for patch in bp['boxes']:
            patch.set_facecolor('white')
            patch.set_alpha(0.8)

    def draw_data_points(self, positions):
        """Draw individual data points overlay."""
        for pos, grp, color in zip(positions, self.group_names, self.group_colors):
            data = self.plot_data[grp]

            if self.args.point_style == 'strip':
                jitter = np.random.normal(0, self.args.jitter, len(data))
            elif self.args.point_style == 'jitter':
                jitter = np.random.uniform(-self.args.jitter, self.args.jitter, len(data))
            else:  # swarm
                jitter = np.random.normal(0, self.args.jitter * 0.5, len(data))

            point_positions = pos + jitter

            if self.args.point_color == 'auto':
                point_color = color
            else:
                point_color = self.args.point_color

            if self.args.orientation == 'vertical':
                self.ax.scatter(
                    point_positions, data,
                    s=self.args.point_size, alpha=self.args.point_alpha,
                    color=point_color, edgecolor=self.args.point_edge_color,
                    zorder=3
                )
            else:
                self.ax.scatter(
                    data, point_positions,
                    s=self.args.point_size, alpha=self.args.point_alpha,
                    color=point_color, edgecolor=self.args.point_edge_color,
                    zorder=3
                )

    def draw_mean_line(self, positions):
        """Draw horizontal/vertical line at grand mean."""
        grand_mean = np.mean(self.values)
        if self.args.orientation == 'vertical':
            self.ax.axhline(grand_mean, color='red', linestyle='--', alpha=0.5, linewidth=1)
        else:
            self.ax.axvline(grand_mean, color='red', linestyle='--', alpha=0.5, linewidth=1)

    def add_sample_sizes(self, positions):
        """Add sample size labels per group."""
        for pos, grp in zip(positions, self.group_names):
            n = len(self.plot_data[grp])
            if self.args.orientation == 'vertical':
                y_pos = self.ax.get_ylim()[0] if self.args.n_position == 'bottom' else self.ax.get_ylim()[1]
                self.ax.text(pos, y_pos, f'n={n}', ha='center', va='bottom' if self.args.n_position == 'bottom' else 'top',
                             fontsize=self.args.tick_size)
            else:
                x_pos = self.ax.get_xlim()[0] if self.args.n_position == 'bottom' else self.ax.get_xlim()[1]
                self.ax.text(x_pos, pos, f'n={n}', ha='left' if self.args.n_position == 'bottom' else 'right', va='center',
                             fontsize=self.args.tick_size)

    def add_mean_labels(self, positions):
        """Add mean value labels above each group."""
        for pos, grp in zip(positions, self.group_names):
            mean_val = np.mean(self.plot_data[grp])
            if self.args.orientation == 'vertical':
                self.ax.text(pos, mean_val, f'{mean_val:.2f}', ha='center', va='bottom',
                             fontsize=self.args.tick_size)
            else:
                self.ax.text(mean_val, pos, f'{mean_val:.2f}', ha='left', va='center',
                             fontsize=self.args.tick_size)

    def perform_statistics(self):
        """Perform statistical tests."""
        n_groups = len(self.group_names)

        if self.args.stats == 'none':
            return {}

        results = {
            'groups': self.group_names,
            'sample_sizes': {grp: len(self.plot_data[grp]) for grp in self.group_names},
            'means': {grp: np.mean(self.plot_data[grp]) for grp in self.group_names},
            'medians': {grp: np.median(self.plot_data[grp]) for grp in self.group_names},
            'stds': {grp: np.std(self.plot_data[grp]) for grp in self.group_names},
        }

        # Determine tests to run
        if self.args.stats == 'auto':
            if n_groups == 2:
                pairwise_comparisons = [(self.group_names[0], self.group_names[1])]
            else:
                # Kruskal-Wallis test
                data_list = [self.plot_data[grp] for grp in self.group_names]
                h_stat, h_pval = stats.kruskal(*data_list)
                results['kruskal_wallis'] = {'statistic': h_stat, 'p_value': h_pval}
                pairwise_comparisons = list(self._get_all_pairs())
        elif self.args.stats == 'all_pairs':
            pairwise_comparisons = list(self._get_all_pairs())
        elif self.args.stats == 'vs_first':
            pairwise_comparisons = [(self.group_names[0], grp) for grp in self.group_names[1:]]
        elif self.args.stats == 'vs_last':
            pairwise_comparisons = [(self.group_names[-1], grp) for grp in self.group_names[:-1]]
        else:
            pairwise_comparisons = []

        # Perform pairwise tests
        p_values = []
        for grp1, grp2 in pairwise_comparisons:
            data1 = self.plot_data[grp1]
            data2 = self.plot_data[grp2]

            if self.args.test == 'auto':
                test_type = 'mannwhitney'
            else:
                test_type = self.args.test

            if test_type == 'mannwhitney':
                stat, pval = stats.mannwhitneyu(data1, data2)
            elif test_type == 'ttest':
                stat, pval = stats.ttest_ind(data1, data2)
            elif test_type == 'kruskal':
                stat, pval = stats.kruskal(data1, data2)

            p_values.append(pval)
            results.setdefault('pairwise', []).append({
                'pair': (grp1, grp2),
                'statistic': stat,
                'p_value': pval,
            })

        # Multiple testing correction
        if p_values:
            if self.args.correction == 'bonferroni':
                corrected_pvals = np.minimum(np.array(p_values) * len(p_values), 1.0)
            elif self.args.correction == 'fdr_bh':
                corrected_pvals = self._fdr_correction(np.array(p_values))
            else:
                corrected_pvals = np.array(p_values)

            for i, result in enumerate(results.get('pairwise', [])):
                result['p_value_corrected'] = corrected_pvals[i]
                result['significant'] = corrected_pvals[i] < self.args.alpha_level

        self.stats_results = results
        return results

    def _get_all_pairs(self):
        """Get all pairwise combinations."""
        for i, grp1 in enumerate(self.group_names):
            for grp2 in self.group_names[i+1:]:
                yield (grp1, grp2)

    def _fdr_correction(self, p_values):
        """Benjamini-Hochberg FDR correction."""
        n = len(p_values)
        sorted_idx = np.argsort(p_values)
        sorted_pvals = p_values[sorted_idx]

        # Calculate adjusted p-values
        adjusted = np.zeros(n)
        for i, idx in enumerate(sorted_idx):
            adjusted[idx] = sorted_pvals[i] * n / (i + 1)

        return np.minimum(adjusted, 1.0)

    def add_statistical_annotations(self, positions):
        """Add statistical annotations to plot."""
        if not self.stats_results.get('pairwise'):
            return

        y_max = self.ax.get_ylim()[1] if self.args.orientation == 'vertical' else self.ax.get_xlim()[1]
        y_start = y_max * 0.95

        for i, result in enumerate(self.stats_results['pairwise']):
            grp1, grp2 = result['pair']
            significant = result.get('significant', False)
            pval = result.get('p_value_corrected', result.get('p_value'))

            if significant or self.args.show_ns:
                pos1 = positions[self.group_names.index(grp1)]
                pos2 = positions[self.group_names.index(grp2)]

                if self.args.stat_annotation_style == 'bracket':
                    self._draw_bracket(pos1, pos2, y_start + i * (y_max * 0.05), pval, significant)
                else:
                    self._draw_text_annotation(pos1, pos2, y_start + i * (y_max * 0.05), pval, significant)

    def _draw_bracket(self, x1, x2, y, pval, significant):
        """Draw bracket annotation."""
        if self.args.orientation == 'vertical':
            # Vertical lines
            self.ax.plot([x1, x1], [y, y + (self.ax.get_ylim()[1] * 0.02)], 'k-', linewidth=self.args.bracket_linewidth)
            self.ax.plot([x2, x2], [y, y + (self.ax.get_ylim()[1] * 0.02)], 'k-', linewidth=self.args.bracket_linewidth)
            # Horizontal line
            self.ax.plot([x1, x2], [y + (self.ax.get_ylim()[1] * 0.02), y + (self.ax.get_ylim()[1] * 0.02)],
                         'k-', linewidth=self.args.bracket_linewidth)
            # Significance marker
            if significant:
                if pval < 0.001:
                    marker = '***'
                elif pval < 0.01:
                    marker = '**'
                else:
                    marker = '*'
            else:
                marker = 'ns' if self.args.show_ns else ''

            if marker:
                self.ax.text((x1 + x2) / 2, y + (self.ax.get_ylim()[1] * 0.03), marker, ha='center', va='bottom',
                             fontsize=self.args.tick_size)

    def _draw_text_annotation(self, x1, x2, y, pval, significant):
        """Draw text annotation."""
        if significant:
            if pval < 0.001:
                marker = '***'
            elif pval < 0.01:
                marker = '**'
            else:
                marker = '*'
        else:
            marker = 'ns' if self.args.show_ns else ''

        if marker:
            self.ax.text((x1 + x2) / 2, y, marker, ha='center', va='bottom', fontsize=self.args.tick_size)

    def configure_axes(self, positions):
        """Configure axis labels, limits, and formatting."""
        # Set labels
        if self.args.orientation == 'vertical':
            self.ax.set_xlabel(self.args.xlabel or self.args.group_col, fontsize=self.args.axis_label_size)
            self.ax.set_ylabel(self.args.ylabel or self.args.value_col, fontsize=self.args.axis_label_size)
            self.ax.set_xticks(positions)
            self.ax.set_xticklabels(self.group_names, fontsize=self.args.tick_size)

            if self.args.ylim:
                ymin, ymax = map(float, self.args.ylim.split(','))
                self.ax.set_ylim(ymin, ymax)

            if self.args.y_log:
                self.ax.set_yscale('log')
        else:
            self.ax.set_ylabel(self.args.xlabel or self.args.group_col, fontsize=self.args.axis_label_size)
            self.ax.set_xlabel(self.args.ylabel or self.args.value_col, fontsize=self.args.axis_label_size)
            self.ax.set_yticks(positions)
            self.ax.set_yticklabels(self.group_names, fontsize=self.args.tick_size)

            if self.args.xlim:
                xmin, xmax = map(float, self.args.xlim.split(','))
                self.ax.set_xlim(xmin, xmax)

            if self.args.y_log:
                self.ax.set_xscale('log')

        # Set title
        if self.args.title:
            self.ax.set_title(self.args.title, fontsize=self.args.axis_label_size, pad=10)

    def generate_plot(self):
        """Generate the main plot."""
        self.load_data()
        self.setup_colors()
        self.prepare_plot_data()
        self.create_figure()

        # Positions for groups
        positions = np.arange(len(self.group_names))

        # Draw plots based on type
        if self.args.plot_type == 'box':
            self.draw_boxplot(positions)
        elif self.args.plot_type == 'violin':
            self.draw_violin(positions)
        elif self.args.plot_type == 'both':
            self.draw_violin(positions)
            self.draw_boxplot(positions)
        elif self.args.plot_type == 'raincloud':
            self.draw_raincloud(positions)

        # Add overlays
        if self.args.show_points and self.args.plot_type != 'raincloud':
            self.draw_data_points(positions)

        if self.args.show_mean:
            self.draw_mean_line(positions)

        if self.args.add_mean_line:
            self.draw_mean_line(positions)

        if self.args.show_n:
            self.add_sample_sizes(positions)

        if self.args.show_mean_label:
            self.add_mean_labels(positions)

        # Add grid and configure axes
        self.add_grid()
        self.configure_axes(positions)

        # Perform statistics
        self.perform_statistics()

        # Add statistical annotations
        if self.args.stats != 'none':
            self.add_statistical_annotations(positions)

        plt.tight_layout()

    def save_plot(self):
        """Save plot to file(s)."""
        output_path = Path(self.args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Save main format
        self.fig.savefig(output_path, dpi=self.args.dpi, bbox_inches='tight')
        print(f"Saved plot: {output_path}")

        # Save SVG if requested
        if self.args.output_svg:
            svg_path = output_path.with_suffix('.svg')
            self.fig.savefig(svg_path, format='svg', bbox_inches='tight')
            print(f"Saved SVG: {svg_path}")

    def save_statistics(self):
        """Save statistics table."""
        if not self.args.output_stats:
            return

        stats_path = Path(self.args.output_stats)
        stats_path.parent.mkdir(parents=True, exist_ok=True)

        # Build statistics dataframe
        if self.stats_results.get('pairwise'):
            rows = []
            for result in self.stats_results['pairwise']:
                grp1, grp2 = result['pair']
                rows.append({
                    'Group1': grp1,
                    'Group2': grp2,
                    'Statistic': result.get('statistic', ''),
                    'P-value': result.get('p_value', ''),
                    'P-value (corrected)': result.get('p_value_corrected', ''),
                    'Significant': result.get('significant', False),
                })
            stats_df = pd.DataFrame(rows)
            stats_df.to_csv(stats_path, sep='\t', index=False)
            print(f"Saved statistics: {stats_path}")

    def print_summary(self):
        """Print summary statistics."""
        print("\nSummary Statistics:")
        print("=" * 60)
        for grp in self.group_names:
            data = self.plot_data[grp]
            print(f"\n{grp}:")
            print(f"  n = {len(data)}")
            print(f"  mean = {np.mean(data):.4f}")
            print(f"  median = {np.median(data):.4f}")
            print(f"  std = {np.std(data):.4f}")
            print(f"  min = {np.min(data):.4f}")
            print(f"  max = {np.max(data):.4f}")

        if self.stats_results.get('pairwise'):
            print("\nStatistical Tests:")
            print("-" * 60)
            for result in self.stats_results['pairwise']:
                grp1, grp2 = result['pair']
                pval = result.get('p_value_corrected', result.get('p_value'))
                sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
                print(f"  {grp1} vs {grp2}: p={pval:.4f} {sig}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Generate publication-quality boxplots, violin plots, and raincloud plots',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Required arguments
    parser.add_argument('--input', required=True, help='Input TSV/CSV file')
    parser.add_argument('--output', required=True, help='Output file (PNG/SVG/PDF)')

    # Data shape — two mutually compatible modes:
    #   LONG: specify both --value-col and --group-col.
    #   WIDE: specify --wide-cols "colA,colB,colC" (or "auto" for all numeric
    #         columns). In wide format every row is one observation and every
    #         selected column is one group.
    # If none of --value-col / --group-col / --wide-cols is supplied, the
    # script auto-detects wide format when the file has >=2 numeric columns.
    parser.add_argument('--value-col', default=None,
                        help='LONG format: column name for numeric values. '
                             'Omit together with --group-col if using --wide-cols.')
    parser.add_argument('--group-col', default=None,
                        help='LONG format: column name for grouping.')
    parser.add_argument('--wide-cols', default=None,
                        help='WIDE format: comma-separated column names (each becomes '
                             'a group), or "auto" to use every numeric column.')

    # Plot style
    parser.add_argument('--plot-type', default='box', choices=['box', 'violin', 'both', 'raincloud'],
                        help='Type of plot')
    parser.add_argument('--orientation', default='vertical', choices=['vertical', 'horizontal'],
                        help='Plot orientation')
    parser.add_argument('--width', type=float, default=0.6, help='Box/violin width')
    parser.add_argument('--violin-bw', default='scott', help='Violin bandwidth')

    # Groups and colors
    parser.add_argument('--group-order', help='Comma-separated group order')
    parser.add_argument('--group-colors', help='Comma-separated hex colors')
    parser.add_argument('--palette', default='Set2',
                        choices=['tab10', 'Set2', 'Set3', 'Pastel1', 'custom'],
                        help='Color palette')
    parser.add_argument('--color-alpha', type=float, default=0.75, help='Fill transparency')

    # Data points
    parser.add_argument('--show-points', action='store_true', help='Show individual data points')
    parser.add_argument('--point-style', default='strip', choices=['strip', 'swarm', 'jitter'],
                        help='Point style')
    parser.add_argument('--jitter', type=float, default=0.05, help='Jitter amount')
    parser.add_argument('--point-size', type=float, default=4, help='Point size')
    parser.add_argument('--point-alpha', type=float, default=0.7, help='Point transparency')
    parser.add_argument('--point-color', default='auto', help='Point color (hex or auto)')
    parser.add_argument('--point-edge-color', default='none',
                        choices=['none', 'black', 'white'], help='Point edge color')

    # Box options
    parser.add_argument('--show-mean', action='store_true', help='Show mean marker')
    parser.add_argument('--notch', action='store_true', help='Notched boxplot')
    parser.add_argument('--flier-size', type=float, default=4, help='Outlier size')
    parser.add_argument('--flier-symbol', default='o', help='Outlier marker')
    parser.add_argument('--show-fliers', action='store_true', default=True, help='Show outliers')
    parser.add_argument('--box-linewidth', type=float, default=1.2, help='Box line width')
    parser.add_argument('--median-color', default='black', help='Median line color')
    parser.add_argument('--median-linewidth', type=float, default=2.0, help='Median line width')

    # Violin options
    parser.add_argument('--violin-inner', default='quartile',
                        choices=['box', 'quartile', 'point', 'None'],
                        help='Violin inner plot')
    parser.add_argument('--violin-scale', default='width',
                        choices=['width', 'area', 'count'], help='Violin scale')
    parser.add_argument('--half-violin', default='right',
                        choices=['left', 'right', 'full'], help='Half violin direction')

    # Statistics
    parser.add_argument('--stats', default='auto',
                        choices=['none', 'auto', 'all_pairs', 'vs_first', 'vs_last'],
                        help='Statistical tests to perform')
    parser.add_argument('--test', default='auto',
                        choices=['mannwhitney', 'ttest', 'kruskal', 'auto'],
                        help='Statistical test type')
    parser.add_argument('--correction', default='bonferroni',
                        choices=['bonferroni', 'fdr_bh', 'none'],
                        help='Multiple testing correction')
    parser.add_argument('--alpha-level', type=float, default=0.05, help='Significance threshold')
    parser.add_argument('--show-ns', action='store_true', default=True,
                        help='Show "ns" for non-significant pairs')
    parser.add_argument('--stat-annotation-style', default='bracket',
                        choices=['bracket', 'text'], help='Annotation style')
    parser.add_argument('--bracket-linewidth', type=float, default=1.2, help='Bracket line width')

    # Axes
    parser.add_argument('--ylim', help='Y-axis limits (min,max)')
    parser.add_argument('--xlim', help='X-axis limits (min,max)')
    parser.add_argument('--xlabel', help='X-axis label')
    parser.add_argument('--ylabel', help='Y-axis label')
    parser.add_argument('--title', help='Plot title')
    parser.add_argument('--y-log', action='store_true', help='Log-scale value axis')
    parser.add_argument('--add-mean-line', action='store_true', help='Add grand mean line')

    # Summary stats
    parser.add_argument('--show-n', action='store_true', default=True, help='Show sample size')
    parser.add_argument('--n-position', default='bottom', choices=['top', 'bottom'],
                        help='Sample size position')
    parser.add_argument('--show-mean-label', action='store_true', help='Show mean values')

    # Style
    parser.add_argument('--font-family', default='Arial', help='Font family')
    parser.add_argument('--base-fontsize', type=float, default=11, help='Base font size')
    parser.add_argument('--axis-label-size', type=float, default=12, help='Axis label size')
    parser.add_argument('--tick-size', type=float, default=10, help='Tick label size')
    parser.add_argument('--legend-size', type=float, default=10, help='Legend size')
    parser.add_argument('--fig-width', type=float, default=8, help='Figure width')
    parser.add_argument('--fig-height', type=float, default=6, help='Figure height')
    parser.add_argument('--dpi', type=int, default=300, help='DPI')
    parser.add_argument('--spine-style', default='minimal',
                        choices=['all', 'minimal', 'none'],
                        help='Spine style')
    parser.add_argument('--grid', default='y',
                        choices=['none', 'x', 'y', 'both'],
                        help='Grid lines')
    parser.add_argument('--grid-alpha', type=float, default=0.3, help='Grid transparency')
    parser.add_argument('--grid-linestyle', default='--',
                        choices=['--', ':', '-'],
                        help='Grid line style')
    parser.add_argument('--background', default='white',
                        choices=['white', 'light', 'dark'],
                        help='Background style')

    # Output
    parser.add_argument('--output-svg', action='store_true', help='Also save SVG')
    parser.add_argument('--output-stats', help='Save statistics table')

    args = parser.parse_args()
    init_style(
        font_family=getattr(args, 'font_family', None),
        font_size=getattr(args, 'base_fontsize', None),
    )

    try:
        generator = PlotGenerator(args)
        generator.generate_plot()
        generator.print_summary()
        generator.save_plot()
        generator.save_statistics()
        print("\nPlot generation complete!")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
