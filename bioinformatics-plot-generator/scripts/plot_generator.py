#!/usr/bin/env python3

import argparse
import os
import textwrap
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), '..', '_shared'))
from plot_style import init_style
from matplotlib.ticker import MaxNLocator
from scipy import stats

try:
    from adjustText import adjust_text
    HAS_ADJUSTTEXT = True
except ImportError:
    HAS_ADJUSTTEXT = False


# =============================
# Utilities
# =============================
def read_table(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine="python")


def ensure_parent_dir(path: str) -> None:
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)


def coerce_numeric(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def wrap_title(title: Optional[str], width: int = 60) -> Optional[str]:
    if not title:
        return title
    return "\n".join(textwrap.wrap(title, width=width))


def parse_list_or_file(x: Optional[str]) -> List[str]:
    if x is None or str(x).strip() == "":
        return []
    if os.path.exists(x):
        vals = []
        with open(x) as f:
            for line in f:
                line = line.strip()
                if line:
                    vals.append(line)
        return vals
    return [i.strip() for i in str(x).split(",") if i.strip()]


def parse_group_order(x: Optional[str]) -> Optional[List[str]]:
    if x is None or str(x).strip() == "":
        return None
    return [i.strip() for i in str(x).split(",") if i.strip()]


def parse_xlim_ylim(x: Optional[str]) -> Optional[Tuple[float, float]]:
    if x is None or str(x).strip() == "":
        return None
    parts = [i.strip() for i in str(x).split(",")]
    if len(parts) != 2:
        raise ValueError("Axis limits must be in format: min,max")
    return float(parts[0]), float(parts[1])


def pick_feature_col(df: pd.DataFrame, feature_col: Optional[str]) -> str:
    if feature_col:
        if feature_col not in df.columns:
            raise ValueError(f"Feature column '{feature_col}' not found.")
        return feature_col

    candidates = ["gene", "symbol", "feature", "id", "name"]
    lower_map = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c in lower_map:
            return lower_map[c]

    return df.columns[0]


def format_pvalue(p: float) -> str:
    if p < 1e-4:
        return f"p = {p:.1e}"
    return f"p = {p:.4f}"


def significance_stars(p: float) -> str:
    if p < 1e-4:
        return "****"
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def save_figure(fig: plt.Figure, output: str, dpi: int = 300) -> None:
    ensure_parent_dir(output)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def apply_style(
    font_family: str = "Arial",
    base_fontsize: float = 11,
    title_size: float = 13,
    axis_label_size: float = 12,
    tick_size: float = 10,
    legend_size: float = 10,
):
    plt.rcParams.update({
        "font.family": font_family,
        "font.size": base_fontsize,
        "axes.labelsize": axis_label_size,
        "axes.titlesize": title_size,
        "xtick.labelsize": tick_size,
        "ytick.labelsize": tick_size,
        "legend.fontsize": legend_size,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 1.0,
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
    })


def get_figsize(plot_type: str, fig_width: Optional[float], fig_height: Optional[float]):
    defaults = {
        "volcano": (6.8, 5.8),
        "heatmap": (7.0, 6.0),
        "boxplot": (6.6, 5.6),
        "violin": (6.6, 5.6),
        "scatter": (6.2, 5.6),
        "correlation": (6.2, 5.6),
    }
    w, h = defaults.get(plot_type, (6.5, 5.5))
    return fig_width or w, fig_height or h


def add_stat_bracket(ax, x1, x2, y, h, text, lw=1.2):
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=lw, c="black")
    ax.text((x1 + x2) / 2, y + h, text, ha="center", va="bottom", fontsize=10)


# =============================
# Volcano
# =============================
def make_volcano(
    df: pd.DataFrame,
    feature_col: str,
    x_col: str,
    p_col: str,
    output: str,
    title: Optional[str],
    fc_cutoff: float,
    p_cutoff: float,
    annotate_top_n: int,
    top_up_n: int,
    top_down_n: int,
    highlight_up: List[str],
    highlight_down: List[str],
    label_mode: str,
    max_labels: int,
    adjust_labels_flag: bool,
    annotation_arrow: bool,
    annotation_size: float,
    point_size: float,
    alpha: float,
    fig_width: Optional[float],
    fig_height: Optional[float],
    xlim: Optional[Tuple[float, float]],
    ylim: Optional[Tuple[float, float]],
    dpi: int,
):
    dat = df[[feature_col, x_col, p_col]].copy()
    dat[x_col] = coerce_numeric(dat[x_col])
    dat[p_col] = coerce_numeric(dat[p_col])
    dat = dat.replace([np.inf, -np.inf], np.nan).dropna()

    dat = dat[dat[p_col] > 0].copy()
    if len(dat) == 0:
        raise ValueError("No valid rows remain after removing missing/non-positive p-values.")

    dat[p_col] = dat[p_col].clip(lower=np.nextafter(0, 1))
    dat["neglog10p"] = -np.log10(dat[p_col])
    dat["_feature_upper"] = dat[feature_col].astype(str).str.upper()

    dat["category"] = "nonsig"
    dat.loc[(dat[x_col] >= fc_cutoff) & (dat[p_col] < p_cutoff), "category"] = "up"
    dat.loc[(dat[x_col] <= -fc_cutoff) & (dat[p_col] < p_cutoff), "category"] = "down"

    highlight_up = set([x.upper() for x in highlight_up])
    highlight_down = set([x.upper() for x in highlight_down])

    dat["is_highlight_up"] = dat["_feature_upper"].isin(highlight_up)
    dat["is_highlight_down"] = dat["_feature_upper"].isin(highlight_down)
    dat["label_me"] = False

    dat["score"] = dat["neglog10p"] * dat[x_col].abs()

    if label_mode in ["highlight", "top_and_highlight"]:
        dat.loc[dat["is_highlight_up"] | dat["is_highlight_down"], "label_me"] = True

    if label_mode in ["top", "top_and_highlight"]:
        if annotate_top_n > 0:
            idx = dat.sort_values("score", ascending=False).head(annotate_top_n).index
            dat.loc[idx, "label_me"] = True

        if top_up_n > 0:
            up_idx = dat[dat["category"] == "up"].sort_values("score", ascending=False).head(top_up_n).index
            dat.loc[up_idx, "label_me"] = True

        if top_down_n > 0:
            down_idx = dat[dat["category"] == "down"].sort_values("score", ascending=False).head(top_down_n).index
            dat.loc[down_idx, "label_me"] = True

    label_df = dat[dat["label_me"]].copy()
    if len(label_df) > max_labels:
        keep_idx = label_df.sort_values(
            by=["is_highlight_up", "is_highlight_down", "score"],
            ascending=[False, False, False]
        ).head(max_labels).index
        dat["label_me"] = False
        dat.loc[keep_idx, "label_me"] = True
        label_df = dat.loc[keep_idx].copy()

    fig_w, fig_h = get_figsize("volcano", fig_width, fig_height)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    palette = {
        "nonsig": "#BDBDBD",
        "up": "#D55E00",
        "down": "#0072B2",
    }

    for cat in ["nonsig", "down", "up"]:
        sub = dat[dat["category"] == cat]
        ax.scatter(
            sub[x_col],
            sub["neglog10p"],
            s=point_size,
            alpha=alpha,
            c=palette[cat],
            edgecolors="none",
            label=cat,
            zorder=2,
        )

    sub_hu = dat[dat["is_highlight_up"]]
    if len(sub_hu) > 0:
        ax.scatter(
            sub_hu[x_col], sub_hu["neglog10p"],
            s=point_size * 2.3,
            c="#B2182B",
            edgecolors="black",
            linewidths=0.5,
            zorder=4,
            label="highlight_up",
        )

    sub_hd = dat[dat["is_highlight_down"]]
    if len(sub_hd) > 0:
        ax.scatter(
            sub_hd[x_col], sub_hd["neglog10p"],
            s=point_size * 2.3,
            c="#2166AC",
            edgecolors="black",
            linewidths=0.5,
            zorder=4,
            label="highlight_down",
        )

    ax.axvline(fc_cutoff, linestyle="--", linewidth=1, color="gray")
    ax.axvline(-fc_cutoff, linestyle="--", linewidth=1, color="gray")
    ax.axhline(-np.log10(p_cutoff), linestyle="--", linewidth=1, color="gray")

    texts = []
    if len(label_df) > 0:
        for _, row in label_df.iterrows():
            texts.append(
                ax.text(
                    row[x_col],
                    row["neglog10p"],
                    str(row[feature_col]),
                    fontsize=annotation_size,
                    fontweight="bold" if (row["is_highlight_up"] or row["is_highlight_down"]) else "normal",
                    zorder=5,
                )
            )

    if adjust_labels_flag and HAS_ADJUSTTEXT and texts:
        arrowprops = dict(arrowstyle="-", lw=0.5, color="black") if annotation_arrow else None
        adjust_text(texts, ax=ax, arrowprops=arrowprops)

    ax.set_xlabel(x_col)
    ax.set_ylabel(f"-log10({p_col})")
    ax.set_title(wrap_title(title) or "Volcano plot")

    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    ax.legend(frameon=False, loc="upper right")
    save_figure(fig, output, dpi=dpi)


# =============================
# Heatmap
# =============================
def zscore_rows(mat: pd.DataFrame) -> pd.DataFrame:
    arr = mat.values.astype(float)
    mean = arr.mean(axis=1, keepdims=True)
    std = arr.std(axis=1, keepdims=True)
    std[std == 0] = 1.0
    out = (arr - mean) / std
    return pd.DataFrame(out, index=mat.index, columns=mat.columns)


def make_heatmap(
    df: pd.DataFrame,
    output: str,
    index_col: Optional[str],
    title: Optional[str],
    scale_rows: bool,
    cmap: str,
    vmin: Optional[float],
    vmax: Optional[float],
    fig_width: Optional[float],
    fig_height: Optional[float],
    dpi: int,
):
    dat = df.copy()
    if index_col is not None:
        if index_col not in dat.columns:
            raise ValueError(f"Index column '{index_col}' not found.")
        dat = dat.set_index(index_col)

    dat = dat.select_dtypes(include=[np.number])
    if dat.shape[1] == 0:
        raise ValueError("No numeric columns available for heatmap.")

    if scale_rows:
        dat = zscore_rows(dat)

    if fig_width is None:
        fig_width = max(6, min(12, 0.45 * dat.shape[1] + 3))
    if fig_height is None:
        fig_height = max(5, min(14, 0.22 * dat.shape[0] + 2))

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    im = ax.imshow(dat.values, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(im, ax=ax, shrink=0.85)
    cbar.ax.set_ylabel("value", rotation=270, labelpad=15)

    ax.set_xticks(range(dat.shape[1]))
    ax.set_xticklabels(dat.columns, rotation=90)

    if dat.shape[0] <= 80:
        ax.set_yticks(range(dat.shape[0]))
        ax.set_yticklabels(dat.index)
    else:
        ax.set_yticks([])

    ax.set_title(wrap_title(title) or "Heatmap")
    save_figure(fig, output, dpi=dpi)


# =============================
# Boxplot / Violin
# =============================
def two_group_test(
    dat: pd.DataFrame,
    value_col: str,
    group_col: str,
    group_order: List[str],
    test_method: str = "mannwhitney",
):
    if len(group_order) != 2:
        return None

    g1, g2 = group_order
    x = coerce_numeric(dat.loc[dat[group_col] == g1, value_col]).dropna()
    y = coerce_numeric(dat.loc[dat[group_col] == g2, value_col]).dropna()

    if len(x) < 2 or len(y) < 2:
        return None

    if test_method == "ttest":
        stat, p = stats.ttest_ind(x, y, equal_var=False, nan_policy="omit")
        method = "Welch t-test"
    else:
        stat, p = stats.mannwhitneyu(x, y, alternative="two-sided")
        method = "Mann–Whitney U"

    return {"group1": g1, "group2": g2, "p": p, "method": method}


def make_box_or_violin(
    df: pd.DataFrame,
    value_col: str,
    group_col: str,
    output: str,
    kind: str,
    title: Optional[str],
    group_order: Optional[List[str]],
    show_points: bool,
    point_size: float,
    point_alpha: float,
    test_method: str,
    fig_width: Optional[float],
    fig_height: Optional[float],
    dpi: int,
):
    dat = df[[value_col, group_col]].copy()
    dat[value_col] = coerce_numeric(dat[value_col])
    dat[group_col] = dat[group_col].astype(str)
    dat = dat.dropna()

    if len(dat) == 0:
        raise ValueError("No valid rows remain for box/violin plot.")

    if group_order is None:
        group_order = list(pd.unique(dat[group_col]))
    else:
        missing = [g for g in group_order if g not in set(dat[group_col])]
        if missing:
            raise ValueError(f"Groups not found in data: {missing}")

    arrays = [dat.loc[dat[group_col] == g, value_col].values for g in group_order]

    fig_w, fig_h = get_figsize(kind, fig_width, fig_height)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    if kind == "violin":
        parts = ax.violinplot(arrays, showmeans=False, showmedians=True, showextrema=True)
        for pc in parts["bodies"]:
            pc.set_alpha(0.6)
    else:
        ax.boxplot(
            arrays,
            tick_labels=group_order,
            widths=0.6,
            medianprops={"linewidth": 1.8},
            whiskerprops={"linewidth": 1.2},
            capprops={"linewidth": 1.2},
            boxprops={"linewidth": 1.2},
        )

    if show_points:
        rng = np.random.default_rng(42)
        for i, g in enumerate(group_order, start=1):
            vals = dat.loc[dat[group_col] == g, value_col].values
            xj = rng.normal(loc=i, scale=0.05, size=len(vals))
            ax.scatter(xj, vals, s=point_size, alpha=point_alpha, zorder=3)

    test_res = two_group_test(dat, value_col, group_col, group_order, test_method=test_method)
    if test_res is not None:
        ymax = np.nanmax(dat[value_col].values)
        ymin = np.nanmin(dat[value_col].values)
        yr = ymax - ymin if ymax > ymin else 1.0
        y = ymax + 0.08 * yr
        h = 0.04 * yr
        label = f"{significance_stars(test_res['p'])}\n{format_pvalue(test_res['p'])}"
        add_stat_bracket(ax, 1, 2, y, h, label)
        ax.text(
            0.02, 0.98,
            test_res["method"],
            transform=ax.transAxes,
            ha="left", va="top", fontsize=10
        )

    ax.set_xlabel(group_col)
    ax.set_ylabel(value_col)
    ax.set_title(wrap_title(title) or ("Violin plot" if kind == "violin" else "Boxplot"))
    save_figure(fig, output, dpi=dpi)


# =============================
# Scatter / Correlation
# =============================
def make_scatter(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    output: str,
    title: Optional[str],
    add_regression: bool,
    corr_method: str,
    point_size: float,
    alpha: float,
    fig_width: Optional[float],
    fig_height: Optional[float],
    dpi: int,
):
    dat = df[[x_col, y_col]].copy()
    dat[x_col] = coerce_numeric(dat[x_col])
    dat[y_col] = coerce_numeric(dat[y_col])
    dat = dat.replace([np.inf, -np.inf], np.nan).dropna()

    if len(dat) < 3:
        raise ValueError("Scatter plot requires at least 3 valid points.")

    fig_w, fig_h = get_figsize("scatter", fig_width, fig_height)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    ax.scatter(dat[x_col], dat[y_col], s=point_size, alpha=alpha, zorder=3)

    if corr_method == "spearman":
        r, p = stats.spearmanr(dat[x_col], dat[y_col])
        method_label = "Spearman"
    else:
        r, p = stats.pearsonr(dat[x_col], dat[y_col])
        method_label = "Pearson"

    if add_regression:
        slope, intercept, _, _, _ = stats.linregress(dat[x_col], dat[y_col])
        xfit = np.linspace(dat[x_col].min(), dat[x_col].max(), 200)
        yfit = slope * xfit + intercept
        ax.plot(xfit, yfit, linewidth=1.5, zorder=2)

    txt = f"{method_label} r = {r:.3f}\n{format_pvalue(p)}\nn = {len(dat)}"
    ax.text(
        0.03, 0.97, txt,
        transform=ax.transAxes,
        ha="left", va="top", fontsize=10,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.75)
    )

    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_title(wrap_title(title) or "Scatter plot")
    ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
    save_figure(fig, output, dpi=dpi)


# =============================
# CLI
# =============================
def build_parser():
    parser = argparse.ArgumentParser(
        description="Publication-quality bioinformatics plot generator"
    )

    parser.add_argument("--input", required=True, help="Input table")
    parser.add_argument(
        "--plot-type",
        required=True,
        choices=["volcano", "heatmap", "boxplot", "violin", "scatter", "correlation"],
        help="Plot type",
    )
    parser.add_argument("--output", required=True, help="Output plot path")
    parser.add_argument("--title", default=None, help="Plot title")

    # General style
    parser.add_argument("--font-family", default="Arial")
    parser.add_argument("--base-fontsize", type=float, default=11)
    parser.add_argument("--title-size", type=float, default=13)
    parser.add_argument("--axis-label-size", type=float, default=12)
    parser.add_argument("--tick-size", type=float, default=10)
    parser.add_argument("--legend-size", type=float, default=10)
    parser.add_argument("--annotation-size", type=float, default=9)
    parser.add_argument("--fig-width", type=float, default=None)
    parser.add_argument("--fig-height", type=float, default=None)
    parser.add_argument("--dpi", type=int, default=300)

    # Common columns
    parser.add_argument("--feature-col", default=None)
    parser.add_argument("--x-col", default=None)
    parser.add_argument("--y-col", default=None)
    parser.add_argument("--p-col", default=None)
    parser.add_argument("--value-col", default=None)
    parser.add_argument("--group-col", default=None)
    parser.add_argument("--index-col", default=None)

    # Volcano
    parser.add_argument("--fc-cutoff", type=float, default=1.0)
    parser.add_argument("--p-cutoff", type=float, default=0.05)
    parser.add_argument("--annotate-top-n", type=int, default=0)
    parser.add_argument("--top-up-n", type=int, default=0)
    parser.add_argument("--top-down-n", type=int, default=0)
    parser.add_argument("--highlight-up", default=None, help="Comma list or file")
    parser.add_argument("--highlight-down", default=None, help="Comma list or file")
    parser.add_argument(
        "--label-mode",
        choices=["none", "top", "highlight", "top_and_highlight"],
        default="top_and_highlight",
    )
    parser.add_argument("--max-labels", type=int, default=20)
    parser.add_argument("--adjust-labels", action="store_true")
    parser.add_argument("--annotation-arrow", action="store_true")
    parser.add_argument("--point-size", type=float, default=18)
    parser.add_argument("--alpha", type=float, default=0.75)
    parser.add_argument("--xlim", default=None, help="min,max")
    parser.add_argument("--ylim", default=None, help="min,max")

    # Heatmap
    parser.add_argument("--scale-rows", action="store_true")
    parser.add_argument("--cmap", default="coolwarm")
    parser.add_argument("--vmin", type=float, default=None)
    parser.add_argument("--vmax", type=float, default=None)

    # Box/Violin
    parser.add_argument("--group-order", default=None, help="A,B,C")
    parser.add_argument("--show-points", action="store_true")
    parser.add_argument("--point-alpha", type=float, default=0.65)
    parser.add_argument(
        "--test-method",
        choices=["mannwhitney", "ttest"],
        default="mannwhitney",
    )

    # Scatter
    parser.add_argument(
        "--corr-method",
        choices=["pearson", "spearman"],
        default="pearson",
    )
    parser.add_argument("--no-regression", action="store_true")

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    init_style(
        font_family=getattr(args, 'font_family', None),
        font_size=getattr(args, 'base_fontsize', None),
    )

    apply_style(
        font_family=args.font_family,
        base_fontsize=args.base_fontsize,
        title_size=args.title_size,
        axis_label_size=args.axis_label_size,
        tick_size=args.tick_size,
        legend_size=args.legend_size,
    )

    df = read_table(args.input)

    if args.plot_type == "volcano":
        if args.x_col is None or args.p_col is None:
            raise ValueError("Volcano plot requires --x-col and --p-col")
        feature_col = pick_feature_col(df, args.feature_col)

        make_volcano(
            df=df,
            feature_col=feature_col,
            x_col=args.x_col,
            p_col=args.p_col,
            output=args.output,
            title=args.title,
            fc_cutoff=args.fc_cutoff,
            p_cutoff=args.p_cutoff,
            annotate_top_n=args.annotate_top_n,
            top_up_n=args.top_up_n,
            top_down_n=args.top_down_n,
            highlight_up=parse_list_or_file(args.highlight_up),
            highlight_down=parse_list_or_file(args.highlight_down),
            label_mode=args.label_mode,
            max_labels=args.max_labels,
            adjust_labels_flag=args.adjust_labels,
            annotation_arrow=args.annotation_arrow,
            annotation_size=args.annotation_size,
            point_size=args.point_size,
            alpha=args.alpha,
            fig_width=args.fig_width,
            fig_height=args.fig_height,
            xlim=parse_xlim_ylim(args.xlim),
            ylim=parse_xlim_ylim(args.ylim),
            dpi=args.dpi,
        )

    elif args.plot_type == "heatmap":
        make_heatmap(
            df=df,
            output=args.output,
            index_col=args.index_col,
            title=args.title,
            scale_rows=args.scale_rows,
            cmap=args.cmap,
            vmin=args.vmin,
            vmax=args.vmax,
            fig_width=args.fig_width,
            fig_height=args.fig_height,
            dpi=args.dpi,
        )

    elif args.plot_type in ["boxplot", "violin"]:
        if args.value_col is None or args.group_col is None:
            raise ValueError(f"{args.plot_type} requires --value-col and --group-col")
        make_box_or_violin(
            df=df,
            value_col=args.value_col,
            group_col=args.group_col,
            output=args.output,
            kind=args.plot_type,
            title=args.title,
            group_order=parse_group_order(args.group_order),
            show_points=args.show_points,
            point_size=args.point_size,
            point_alpha=args.point_alpha,
            test_method=args.test_method,
            fig_width=args.fig_width,
            fig_height=args.fig_height,
            dpi=args.dpi,
        )

    elif args.plot_type in ["scatter", "correlation"]:
        if args.x_col is None or args.y_col is None:
            raise ValueError(f"{args.plot_type} requires --x-col and --y-col")
        make_scatter(
            df=df,
            x_col=args.x_col,
            y_col=args.y_col,
            output=args.output,
            title=args.title,
            add_regression=not args.no_regression,
            corr_method=args.corr_method,
            point_size=args.point_size,
            alpha=args.alpha,
            fig_width=args.fig_width,
            fig_height=args.fig_height,
            dpi=args.dpi,
        )

    print(f"Saved plot to: {args.output}")


if __name__ == "__main__":
    main()