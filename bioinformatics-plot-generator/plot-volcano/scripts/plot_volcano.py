#!/usr/bin/env python3
"""
Volcano Plot Generator for Differential Expression Analysis

A publication-quality volcano plot generator from DE result tables using only
matplotlib (no seaborn dependency). Supports extensive customization of colors,
thresholds, labels, and styling while maintaining clean aesthetics.
"""

import argparse
import sys
from pathlib import Path
from typing import Tuple, List, Set, Dict, Optional
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), '..', '..', '..', '_shared'))
from plot_style import init_style
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# Optional: try to import adjustText for smart label positioning
try:
    from adjustText import adjust_text
    HAS_ADJUSTTEXT = True
except ImportError:
    HAS_ADJUSTTEXT = False


class VolcanoPlotGenerator:
    """Generate publication-quality volcano plots from DE results."""

    def __init__(self, args):
        """Initialize with command-line arguments."""
        self.args = args
        self.df = None
        self.x_col = args.x_col
        self.p_col = args.p_col
        self.feature_col = None
        self.color_col = args.color_col

        # Store processed data
        self.x_data = None
        self.y_data = None
        self.features = None
        self.color_data = None
        self.categories = None
        self.labels_to_draw = None
        self.label_colors = None

    def validate_input(self) -> bool:
        """Validate input arguments and file."""
        if not Path(self.args.input).exists():
            raise FileNotFoundError(f"Input file not found: {self.args.input}")

        if self.args.fc_cutoff <= 0:
            raise ValueError(f"--fc-cutoff must be > 0, got {self.args.fc_cutoff}")
        if self.args.p_cutoff <= 0 or self.args.p_cutoff > 1:
            raise ValueError(f"--p-cutoff must be in (0, 1], got {self.args.p_cutoff}")
        if self.args.point_size <= 0:
            raise ValueError(f"--point-size must be > 0")

        return True

    def load_data(self):
        """Load and validate input file."""
        # Detect file format
        file_path = self.args.input
        if file_path.endswith('.csv'):
            self.df = pd.read_csv(file_path)
        elif file_path.endswith('.tsv') or file_path.endswith('.txt'):
            self.df = pd.read_csv(file_path, sep='\t')
        else:
            # Try auto-detect
            try:
                self.df = pd.read_csv(file_path, sep='\t')
            except:
                self.df = pd.read_csv(file_path)

        print(f"Loaded {len(self.df)} features from {self.args.input}")

        # Validate required columns
        if self.x_col not in self.df.columns:
            raise ValueError(
                f"log2FC column '{self.x_col}' not found in data. "
                f"Available columns: {list(self.df.columns)}"
            )
        if self.p_col not in self.df.columns:
            raise ValueError(
                f"p-value column '{self.p_col}' not found in data. "
                f"Available columns: {list(self.df.columns)}"
            )

        # Auto-detect feature column if not provided
        if self.args.feature_col:
            self.feature_col = self.args.feature_col
        else:
            self.feature_col = self._auto_detect_feature_col()

        # Validate feature column
        if self.feature_col and self.feature_col not in self.df.columns:
            raise ValueError(
                f"Feature column '{self.feature_col}' not found in data. "
                f"Available columns: {list(self.df.columns)}"
            )

        # Validate color column if provided
        if self.color_col and self.color_col not in self.df.columns:
            raise ValueError(
                f"Color column '{self.color_col}' not found in data. "
                f"Available columns: {list(self.df.columns)}"
            )

    def _auto_detect_feature_col(self) -> Optional[str]:
        """Auto-detect feature/gene name column from common column names."""
        common_names = ['gene', 'Gene', 'symbol', 'Symbol', 'id', 'ID',
                       'feature', 'Feature', 'name', 'Name']

        for name in common_names:
            if name in self.df.columns:
                print(f"Auto-detected feature column: {name}")
                return name

        # If no common name found, use index if it's not default
        if self.df.index.name and self.df.index.name != 'index':
            print(f"Using index as feature names: {self.df.index.name}")
            return self.df.index.name

        print("Warning: No feature column detected. Points will not be labeled.")
        return None

    def prepare_data(self):
        """Prepare and transform data for plotting."""
        # Extract columns
        self.x_data = self.df[self.x_col].values.astype(float)
        p_values = self.df[self.p_col].values.astype(float)

        # Handle log10 transformation
        if self.args.p_is_log:
            self.y_data = p_values  # Already -log10(p)
        else:
            # Convert p-values to -log10(p), handling zeros
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                self.y_data = -np.log10(np.maximum(p_values, 1e-300))

        # Extract feature names
        if self.feature_col:
            self.features = self.df[self.feature_col].values.astype(str)
        else:
            self.features = np.array([f"Feature_{i}" for i in range(len(self.df))])

        # Extract color data if provided
        if self.color_col:
            self.color_data = self.df[self.color_col].values.astype(float)

        # Categorize points: up, down, ns
        self.categories = self._categorize_points()

        print(f"Summary: {len(self.df)} total features")
        print(f"  - Upregulated: {np.sum(self.categories == 'up')}")
        print(f"  - Downregulated: {np.sum(self.categories == 'down')}")
        print(f"  - Non-significant: {np.sum(self.categories == 'ns')}")

    def _categorize_points(self) -> np.ndarray:
        """Categorize points as up, down, or non-significant."""
        categories = np.full(len(self.x_data), 'ns', dtype=object)

        p_threshold = self._get_p_threshold()

        # Upregulated: log2FC >= fc_cutoff AND p < p_cutoff
        up_mask = (self.x_data >= self.args.fc_cutoff) & (self._get_p_values() < p_threshold)
        categories[up_mask] = 'up'

        # Downregulated: log2FC <= -fc_cutoff_neg AND p < p_cutoff
        down_mask = (self.x_data <= -self.args.fc_cutoff_neg) & (self._get_p_values() < p_threshold)
        categories[down_mask] = 'down'

        return categories

    def _get_p_values(self) -> np.ndarray:
        """Get p-values, converting from -log10 if needed."""
        if self.args.p_is_log:
            return 10 ** (-self.y_data)
        else:
            return self.df[self.p_col].values.astype(float)

    def _get_p_threshold(self) -> float:
        """Get p-value threshold (un-log'd if necessary)."""
        return self.args.p_cutoff

    def determine_labels(self):
        """Determine which points to label based on label mode."""
        self.labels_to_draw = {}
        self.label_colors = {}

        # Score for ranking: -log10(p) * |log2FC|
        scores = self.y_data * np.abs(self.x_data)

        # Initialize label dictionary by category
        label_dict = {'up': {}, 'down': {}, 'other': {}}

        if self.args.label_mode in ['top', 'top_and_highlight']:
            # Top upregulated
            up_mask = self.categories == 'up'
            if np.any(up_mask):
                up_indices = np.where(up_mask)[0]
                up_scores = scores[up_indices]
                top_up_idx = up_indices[np.argsort(-up_scores)[:self.args.top_up_n]]
                for idx in top_up_idx:
                    label_dict['up'][idx] = self.features[idx]

            # Top downregulated
            down_mask = self.categories == 'down'
            if np.any(down_mask):
                down_indices = np.where(down_mask)[0]
                down_scores = scores[down_indices]
                top_down_idx = down_indices[np.argsort(-down_scores)[:self.args.top_down_n]]
                for idx in top_down_idx:
                    label_dict['down'][idx] = self.features[idx]

            # Top by score overall
            if self.args.annotate_top_n > 0:
                top_overall_idx = np.argsort(-scores)[:self.args.annotate_top_n]
                for idx in top_overall_idx:
                    if idx not in label_dict['up'] and idx not in label_dict['down']:
                        label_dict['other'][idx] = self.features[idx]

        if self.args.label_mode in ['highlight', 'top_and_highlight']:
            # Highlight groups
            highlight_up_set = self._load_highlight_set(self.args.highlight_up)
            highlight_down_set = self._load_highlight_set(self.args.highlight_down)
            highlight_other_set = self._load_highlight_set(self.args.highlight_other)

            for idx, feature in enumerate(self.features):
                if feature in highlight_up_set:
                    label_dict['up'][idx] = feature
                elif feature in highlight_down_set:
                    label_dict['down'][idx] = feature
                elif feature in highlight_other_set:
                    label_dict['other'][idx] = feature

        # Flatten and enforce max_labels
        all_labels = {}
        for cat, labels in label_dict.items():
            all_labels.update(labels)

        if len(all_labels) > self.args.max_labels:
            # Keep top by score
            all_indices = list(all_labels.keys())
            all_scores = [scores[idx] for idx in all_indices]
            top_idx_list = [all_indices[i] for i in np.argsort(-np.array(all_scores))[:self.args.max_labels]]
            all_labels = {idx: all_labels[idx] for idx in top_idx_list}

        self.labels_to_draw = all_labels

        # Assign label colors
        for idx, label in self.labels_to_draw.items():
            if self.categories[idx] == 'up' or label in self._load_highlight_set(self.args.highlight_up):
                self.label_colors[idx] = self.args.color_highlight_up
            elif self.categories[idx] == 'down' or label in self._load_highlight_set(self.args.highlight_down):
                self.label_colors[idx] = self.args.color_highlight_down
            else:
                self.label_colors[idx] = self.args.highlight_other_color

        print(f"Labeling {len(self.labels_to_draw)} features")

    def _load_highlight_set(self, highlight_arg: Optional[str]) -> Set[str]:
        """Load highlight set from comma-separated string or file."""
        if not highlight_arg:
            return set()

        # Try as file first
        try:
            with open(highlight_arg) as f:
                return set(line.strip() for line in f if line.strip())
        except FileNotFoundError:
            pass

        # Parse as comma-separated
        return set(s.strip() for s in highlight_arg.split(',') if s.strip())

    def determine_point_highlighting(self) -> Dict[str, np.ndarray]:
        """Determine which points to highlight and their properties."""
        highlight_up_set = self._load_highlight_set(self.args.highlight_up)
        highlight_down_set = self._load_highlight_set(self.args.highlight_down)
        highlight_other_set = self._load_highlight_set(self.args.highlight_other)

        highlight_mask = {
            'up': np.array([f in highlight_up_set for f in self.features]),
            'down': np.array([f in highlight_down_set for f in self.features]),
            'other': np.array([f in highlight_other_set for f in self.features])
        }

        return highlight_mask

    def create_plot(self):
        """Create the volcano plot figure."""
        # Set fonts BEFORE creating axes so tick/label artists pick up the
        # correct family at construction time.
        #
        # Common proprietary font names (Helvetica, Arial, Times) are often
        # missing on Linux servers. matplotlib silently falls back to
        # DejaVu Sans when the requested font isn't installed, which is
        # why --font-family Helvetica seemed to "do nothing". To make the
        # flag actually work we map the requested name to the first
        # available substitute from a small alias list, and seed
        # font.sans-serif with that list so matplotlib's own fallback
        # chain also lands on something reasonable.
        import matplotlib.font_manager as _fm
        _installed = {f.name for f in _fm.fontManager.ttflist}
        _alias = {
            'helvetica':        ['Helvetica', 'Nimbus Sans', 'Liberation Sans',
                                 'Arial', 'DejaVu Sans'],
            'arial':            ['Arial', 'Liberation Sans', 'Nimbus Sans',
                                 'Helvetica', 'DejaVu Sans'],
            'times':            ['Times New Roman', 'Times', 'Liberation Serif',
                                 'Nimbus Roman', 'DejaVu Serif'],
            'times new roman':  ['Times New Roman', 'Times', 'Liberation Serif',
                                 'Nimbus Roman', 'DejaVu Serif'],
            'courier':          ['Courier New', 'Courier', 'Liberation Mono',
                                 'Nimbus Mono PS', 'DejaVu Sans Mono'],
        }
        _requested = (self.args.font_family or 'Arial').strip()
        _chain = _alias.get(_requested.lower(), [_requested, 'DejaVu Sans'])
        _resolved = next((f for f in _chain if f in _installed), 'DejaVu Sans')
        if _resolved.lower() != _requested.lower():
            print(f"[font] '{_requested}' not installed; "
                  f"using '{_resolved}' (nearest available).",
                  file=sys.stderr)
        plt.rcParams['font.family'] = _resolved
        plt.rcParams['font.sans-serif'] = _chain + ['DejaVu Sans']
        plt.rcParams['font.size'] = self.args.base_fontsize

        fig, ax = plt.subplots(figsize=(self.args.fig_width, self.args.fig_height),
                               dpi=self.args.dpi)

        # Set background color
        ax.set_facecolor(self.args.background_color)
        fig.patch.set_facecolor('white')

        # Plot points by category
        self._plot_points(ax)

        # Plot threshold lines
        self._plot_thresholds(ax)

        # Add quadrant counts
        if self.args.show_quadrant_counts:
            self._add_quadrant_counts(ax)

        # Add labels
        if self.labels_to_draw:
            self._add_labels(ax)

        # Handle color bar
        if self.color_col:
            self._add_colorbar(ax)

        # Set axes limits and labels
        if self.args.xlim:
            ax.set_xlim(self.args.xlim)
        if self.args.ylim:
            ax.set_ylim(self.args.ylim)

        xlabel = self.args.xlabel or self.x_col
        ylabel = self.args.ylabel or f"-log10({self.p_col})"

        ax.set_xlabel(xlabel, fontsize=self.args.axis_label_size)
        ax.set_ylabel(ylabel, fontsize=self.args.axis_label_size)

        if self.args.title:
            ax.set_title(self.args.title, fontsize=self.args.axis_label_size)

        # Style axes
        ax.tick_params(labelsize=self.args.tick_size)
        for spine in ax.spines.values():
            spine.set_linewidth(self.args.spine_width)

        # Hide selected spines (e.g. publication-style "top-right off")
        _hide_raw = (self.args.hide_spines or '').strip().lower()
        if _hide_raw:
            if _hide_raw in ('all', 'box'):
                _to_hide = ['top', 'right', 'bottom', 'left']
            elif _hide_raw in ('top-right', 'topright', 'tr'):
                _to_hide = ['top', 'right']
            else:
                _to_hide = [s.strip() for s in _hide_raw.split(',') if s.strip()]
            for _name in _to_hide:
                if _name in ax.spines:
                    ax.spines[_name].set_visible(False)
            # Also drop the corresponding ticks so the gap is clean
            if 'top' in _to_hide and 'bottom' not in _to_hide:
                ax.tick_params(top=False)
            if 'right' in _to_hide and 'left' not in _to_hide:
                ax.tick_params(right=False)

        # Add grid
        if self.args.grid == 'major':
            ax.grid(True, which='major', alpha=self.args.grid_alpha, linestyle=':')
        elif self.args.grid == 'minor':
            ax.grid(True, which='minor', alpha=self.args.grid_alpha, linestyle=':')

        plt.tight_layout()
        return fig, ax

    def _plot_points(self, ax):
        """Plot data points with appropriate colors and sizes."""
        highlight_mask = self.determine_point_highlighting()

        # Determine colors
        if self.color_col:
            # Use continuous colormap
            self._plot_colored_points(ax, highlight_mask)
        else:
            # Use categorical colors
            self._plot_categorical_points(ax, highlight_mask)

    def _plot_colored_points(self, ax, highlight_mask):
        """Plot points colored by continuous color column."""
        norm = Normalize(vmin=np.nanmin(self.color_data), vmax=np.nanmax(self.color_data))
        cmap = plt.get_cmap(self.args.colormap)

        # Plot non-highlighted points first
        for cat in ['up', 'down', 'ns']:
            cat_mask = self.categories == cat
            non_hl = cat_mask & ~(highlight_mask['up'] | highlight_mask['down'] | highlight_mask['other'])

            if np.any(non_hl):
                colors = cmap(norm(self.color_data[non_hl]))
                ax.scatter(self.x_data[non_hl], self.y_data[non_hl],
                          c=self.color_data[non_hl], cmap=self.args.colormap, norm=norm,
                          s=self.args.point_size, alpha=self.args.alpha_ns if cat == 'ns' else self.args.alpha,
                          edgecolors='none', rasterized=True)

        # Plot highlighted points on top
        for hl_type in ['up', 'down', 'other']:
            if np.any(highlight_mask[hl_type]):
                colors = cmap(norm(self.color_data[highlight_mask[hl_type]]))
                ax.scatter(self.x_data[highlight_mask[hl_type]], self.y_data[highlight_mask[hl_type]],
                          c=self.color_data[highlight_mask[hl_type]], cmap=self.args.colormap, norm=norm,
                          s=self.args.point_size_highlight, alpha=self.args.alpha,
                          edgecolors='black', linewidths=self.args.edge_width, rasterized=True)

        # Add colorbar
        if self.args.show_colorbar:
            sm = ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax)
            cbar.set_label(self.color_col, fontsize=self.args.legend_size)

    def _plot_categorical_points(self, ax, highlight_mask):
        """Plot points with categorical colors."""
        # Plot non-significant first (background)
        ns_mask = self.categories == 'ns'
        non_hl_ns = ns_mask & ~(highlight_mask['up'] | highlight_mask['down'] | highlight_mask['other'])
        if np.any(non_hl_ns):
            ax.scatter(self.x_data[non_hl_ns], self.y_data[non_hl_ns],
                      c=self.args.color_ns, s=self.args.point_size,
                      alpha=self.args.alpha_ns, edgecolors='none', rasterized=True)

        # Plot downregulated
        down_mask = self.categories == 'down'
        non_hl_down = down_mask & ~(highlight_mask['up'] | highlight_mask['down'] | highlight_mask['other'])
        if np.any(non_hl_down):
            ax.scatter(self.x_data[non_hl_down], self.y_data[non_hl_down],
                      c=self.args.color_down, s=self.args.point_size,
                      alpha=self.args.alpha, edgecolors='none', rasterized=True)

        # Plot upregulated
        up_mask = self.categories == 'up'
        non_hl_up = up_mask & ~(highlight_mask['up'] | highlight_mask['down'] | highlight_mask['other'])
        if np.any(non_hl_up):
            ax.scatter(self.x_data[non_hl_up], self.y_data[non_hl_up],
                      c=self.args.color_up, s=self.args.point_size,
                      alpha=self.args.alpha, edgecolors='none', rasterized=True)

        # Plot highlighted points on top with larger size and edge
        if np.any(highlight_mask['up']):
            ax.scatter(self.x_data[highlight_mask['up']], self.y_data[highlight_mask['up']],
                      c=self.args.color_highlight_up, s=self.args.point_size_highlight,
                      alpha=self.args.alpha, edgecolors='black', linewidths=self.args.edge_width,
                      rasterized=True)

        if np.any(highlight_mask['down']):
            ax.scatter(self.x_data[highlight_mask['down']], self.y_data[highlight_mask['down']],
                      c=self.args.color_highlight_down, s=self.args.point_size_highlight,
                      alpha=self.args.alpha, edgecolors='black', linewidths=self.args.edge_width,
                      rasterized=True)

        if np.any(highlight_mask['other']):
            ax.scatter(self.x_data[highlight_mask['other']], self.y_data[highlight_mask['other']],
                      c=self.args.highlight_other_color, s=self.args.point_size_highlight,
                      alpha=self.args.alpha, edgecolors='black', linewidths=self.args.edge_width,
                      rasterized=True)

    def _plot_thresholds(self, ax):
        """Plot threshold lines for FC and p-value."""
        p_threshold = self._get_p_threshold()
        y_threshold = -np.log10(p_threshold) if not self.args.p_is_log else p_threshold

        if self.args.show_fc_line:
            # Left FC cutoff
            ax.axvline(-self.args.fc_cutoff_neg, color=self.args.threshold_color,
                      linestyle=self.args.threshold_linestyle, linewidth=self.args.threshold_linewidth,
                      zorder=1, alpha=0.7)
            # Right FC cutoff
            ax.axvline(self.args.fc_cutoff, color=self.args.threshold_color,
                      linestyle=self.args.threshold_linestyle, linewidth=self.args.threshold_linewidth,
                      zorder=1, alpha=0.7)

        if self.args.show_p_line:
            ax.axhline(y_threshold, color=self.args.threshold_color,
                      linestyle=self.args.threshold_linestyle, linewidth=self.args.threshold_linewidth,
                      zorder=1, alpha=0.7)

    def _add_quadrant_counts(self, ax):
        """Add counts of features in each quadrant."""
        p_threshold = self._get_p_threshold()
        y_threshold = -np.log10(p_threshold) if not self.args.p_is_log else p_threshold

        # Get axis limits
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        # Count in each quadrant
        quad_up_sig = np.sum((self.x_data >= self.args.fc_cutoff) & (self.y_data >= y_threshold))
        quad_down_sig = np.sum((self.x_data <= -self.args.fc_cutoff_neg) & (self.y_data >= y_threshold))

        # Top-right
        if quad_up_sig > 0:
            ax.text(xlim[1] * 0.95, ylim[1] * 0.95, f"n={quad_up_sig}",
                   ha='right', va='top', fontsize=self.args.quadrant_fontsize,
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

        # Top-left
        if quad_down_sig > 0:
            ax.text(xlim[0] * 0.95, ylim[1] * 0.95, f"n={quad_down_sig}",
                   ha='left', va='top', fontsize=self.args.quadrant_fontsize,
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

    def _add_labels(self, ax):
        """Add text labels to selected points."""
        texts = []

        for idx, label in self.labels_to_draw.items():
            x, y = self.x_data[idx], self.y_data[idx]
            color = self.label_colors[idx]

            # Determine fontweight
            fontweight = 'bold'  # Default
            if label not in self._load_highlight_set(self.args.highlight_up) and \
               label not in self._load_highlight_set(self.args.highlight_down) and \
               label not in self._load_highlight_set(self.args.highlight_other):
                fontweight = 'normal'

            if self.args.label_fontweight:
                fontweight = self.args.label_fontweight

            # Add annotation
            if self.args.annotation_arrow:
                xytext = (x + 0.1, y + 0.1)
                ax.annotate(label, xy=(x, y), xytext=xytext,
                           fontsize=self.args.annotation_size, color=color,
                           fontweight=fontweight, alpha=0.9,
                           arrowprops=dict(arrowstyle='->', color=color, alpha=0.5))
            else:
                text_obj = ax.text(x, y, label, fontsize=self.args.annotation_size,
                                  color=color, fontweight=fontweight, alpha=0.9)
                texts.append(text_obj)

        # Use adjustText if available and requested
        if texts and self.args.adjust_labels and HAS_ADJUSTTEXT:
            try:
                adjust_text(texts, arrowprops=dict(arrowstyle='->', color='gray', alpha=0.3))
            except:
                pass  # Fall back to default positioning

    def _add_colorbar(self, ax):
        """Add colorbar for continuous color column."""
        sm = ScalarMappable(cmap=plt.get_cmap(self.args.colormap),
                           norm=Normalize(vmin=np.nanmin(self.color_data),
                                        vmax=np.nanmax(self.color_data)))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label(self.color_col, fontsize=self.args.legend_size)

    def save_plot(self, fig):
        """Save plot to file."""
        output_path = Path(self.args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        fig.savefig(str(output_path), dpi=self.args.dpi, bbox_inches='tight')
        print(f"Saved plot to: {output_path}")

        # Save SVG if requested
        if self.args.output_svg:
            svg_path = output_path.with_suffix('.svg')
            fig.savefig(str(svg_path), format='svg', bbox_inches='tight')
            print(f"Saved SVG to: {svg_path}")

    def save_annotated_table(self):
        """Save annotated table with category and label information."""
        if not self.args.output_table:
            return

        output_df = self.df.copy()
        output_df['category'] = self.categories
        output_df['labeled'] = False

        for idx in self.labels_to_draw:
            output_df.loc[idx, 'labeled'] = True

        output_path = Path(self.args.output_table)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_df.to_csv(str(output_path), sep='\t', index=False)
        print(f"Saved annotated table to: {output_path}")

    def run(self):
        """Execute the full pipeline."""
        self.validate_input()
        self.load_data()
        self.prepare_data()
        self.determine_labels()
        fig, ax = self.create_plot()
        self.save_plot(fig)
        self.save_annotated_table()
        plt.close(fig)


def main():
    """Parse arguments and run volcano plot generator."""
    parser = argparse.ArgumentParser(
        description='Generate publication-quality volcano plots from DE results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with auto-detected columns
  python plot_volcano.py -i de_results.tsv -x logFC -p p.adj -o volcano.png

  # With custom thresholds and highlighting
  python plot_volcano.py -i results.tsv -x logFC -p p.adj \\
    --fc-cutoff 1.5 --p-cutoff 0.01 \\
    --highlight-up interesting_genes.txt \\
    -o volcano_custom.png

  # With color coding by expression level
  python plot_volcano.py -i results.tsv -x logFC -p p.adj \\
    --color-col baseMean --colormap viridis \\
    -o volcano_colored.png
        """
    )

    # Required arguments
    parser.add_argument('-i', '--input', required=True, help='Input TSV/CSV file')
    parser.add_argument('-x', '--x-col', required=True, help='log2FC column name')
    parser.add_argument('-p', '--p-col', required=True, help='p-value or adjusted p-value column')

    # Column options
    parser.add_argument('--feature-col', help='Gene/feature name column (auto-detected if omitted)')
    parser.add_argument('--color-col', help='Numeric column to color points by (e.g., baseMean)')

    # Thresholds
    parser.add_argument('--fc-cutoff', type=float, default=1.0, help='|log2FC| cutoff (default 1.0)')
    parser.add_argument('--p-cutoff', type=float, default=0.05, help='Significance cutoff (default 0.05)')
    parser.add_argument('--fc-cutoff-neg', type=float, default=None,
                       help='Separate left-side FC cutoff (default: same as fc-cutoff)')
    parser.add_argument('--p-is-log', action='store_true', help='Flag if p-col is already -log10(p)')

    # Labeling
    parser.add_argument('--label-mode', choices=['none', 'top', 'highlight', 'top_and_highlight'],
                       default='top_and_highlight', help='Labeling strategy')
    parser.add_argument('--annotate-top-n', type=int, default=0, help='Label top N by score')
    parser.add_argument('--top-up-n', type=int, default=10, help='Label top N upregulated')
    parser.add_argument('--top-down-n', type=int, default=10, help='Label top N downregulated')
    parser.add_argument('--highlight-up', help='Comma-separated genes or file to highlight red')
    parser.add_argument('--highlight-down', help='Comma-separated genes or file to highlight blue')
    parser.add_argument('--highlight-other', help='Comma-separated genes or file with custom color')
    parser.add_argument('--highlight-other-color', default='#FF6B35', help='Hex color for highlight-other')
    parser.add_argument('--max-labels', type=int, default=20, help='Max total labels')
    parser.add_argument('--adjust-labels', action='store_true', help='Use adjustText for smart positioning')
    parser.add_argument('--annotation-arrow', action='store_true', help='Draw arrow from label to point')
    parser.add_argument('--annotation-size', type=float, default=8, help='Font size for labels')
    parser.add_argument('--label-fontweight', choices=['normal', 'bold'],
                       help='Font weight for labels (auto-selected if omitted)')

    # Colors
    parser.add_argument('--color-up', default='#D55E00', help='Color for upregulated')
    parser.add_argument('--color-down', default='#0072B2', help='Color for downregulated')
    parser.add_argument('--color-ns', default='#BDBDBD', help='Color for non-significant')
    parser.add_argument('--color-highlight-up', default='#B2182B', help='Color for highlight-up')
    parser.add_argument('--color-highlight-down', default='#2166AC', help='Color for highlight-down')
    parser.add_argument('--colormap', default='viridis', help='Colormap for --color-col')
    parser.add_argument('--show-colorbar', action='store_true', help='Show colorbar when color-col is active')

    # Points
    parser.add_argument('--point-size', type=float, default=20, help='Point size')
    parser.add_argument('--point-size-highlight', type=float, default=45, help='Size for highlighted points')
    parser.add_argument('--alpha', type=float, default=0.75, help='Transparency')
    parser.add_argument('--alpha-ns', type=float, default=0.4, help='Transparency for non-significant')
    parser.add_argument('--edge-width', type=float, default=0.4, help='Edge linewidth for highlights')

    # Threshold lines
    parser.add_argument('--threshold-linestyle', choices=['dashed', 'dotted', 'solid'],
                       default='dashed', help='Threshold line style')
    parser.add_argument('--threshold-linewidth', type=float, default=0.9, help='Threshold line width')
    parser.add_argument('--threshold-color', default='gray', help='Threshold line color')
    parser.add_argument('--show-fc-line', action='store_true', default=True,
                       help='Show FC threshold lines (default True)')
    parser.add_argument('--show-p-line', action='store_true', default=True,
                       help='Show p-value threshold line (default True)')

    # Axes and layout
    parser.add_argument('--xlim', help='X-axis limits (min,max)')
    parser.add_argument('--ylim', help='Y-axis limits (min,max)')
    parser.add_argument('--xlabel', help='Custom x-axis label')
    parser.add_argument('--ylabel', help='Custom y-axis label')
    parser.add_argument('--title', help='Plot title')
    parser.add_argument('--fig-width', type=float, default=8, help='Figure width')
    parser.add_argument('--fig-height', type=float, default=6, help='Figure height')
    parser.add_argument('--dpi', type=int, default=300, help='DPI for output')

    # Quadrant annotations
    parser.add_argument('--show-quadrant-counts', action='store_true', default=True,
                       help='Show n= counts in quadrants')
    parser.add_argument('--quadrant-fontsize', type=float, default=9, help='Quadrant label font size')

    # Style
    parser.add_argument('--font-family', default='Arial', help='Font family')
    parser.add_argument('--base-fontsize', type=float, default=11, help='Base font size')
    parser.add_argument('--axis-label-size', type=float, default=12, help='Axis label font size')
    parser.add_argument('--tick-size', type=float, default=10, help='Tick font size')
    parser.add_argument('--legend-size', type=float, default=10, help='Legend font size')
    parser.add_argument('--spine-width', type=float, default=1.0, help='Spine linewidth')
    parser.add_argument('--hide-spines', default='',
                        help='Comma-separated spines to hide. '
                             'Choices: top,right,bottom,left. '
                             'Shortcuts: "top-right" hides top+right; '
                             '"all" hides all four; "" (default) hides none. '
                             'Example: --hide-spines top,right')
    parser.add_argument('--grid', choices=['none', 'major', 'minor'], default='none', help='Grid style')
    parser.add_argument('--grid-alpha', type=float, default=0.3, help='Grid transparency')
    parser.add_argument('--background-color', default='white', help='Axes background color')

    # Output
    parser.add_argument('-o', '--output', default='volcano.png', help='Output file path')
    parser.add_argument('--output-svg', action='store_true', help='Also save SVG')
    parser.add_argument('--output-table', help='Save annotated table as TSV')

    args = parser.parse_args()
    init_style(
        font_family=getattr(args, 'font_family', None),
        font_size=getattr(args, 'base_fontsize', None),
    )

    # Set default for fc-cutoff-neg if not provided
    if args.fc_cutoff_neg is None:
        args.fc_cutoff_neg = args.fc_cutoff

    # Parse xlim and ylim
    if args.xlim:
        args.xlim = tuple(map(float, args.xlim.split(',')))
    if args.ylim:
        args.ylim = tuple(map(float, args.ylim.split(',')))

    try:
        generator = VolcanoPlotGenerator(args)
        generator.run()
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
