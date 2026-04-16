#!/usr/bin/env python3
"""
Plot Survival Skill
Generate Kaplan-Meier survival curves with statistical testing and annotations.
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict

import matplotlib.pyplot as plt
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), '..', '..', '..', '_shared'))
from plot_style import init_style
import numpy as np
import pandas as pd

try:
    from scipy import stats as _scipy_stats
    _HAVE_SCIPY = True
except ImportError:
    _scipy_stats = None
    _HAVE_SCIPY = False


class _FallbackStats:
    """Minimal scipy.stats replacements using pure numpy."""
    class chi2:
        @staticmethod
        def cdf(x, df):
            return _gammainc_lower(df / 2.0, np.asarray(x, float) / 2.0)
    class norm:
        @staticmethod
        def cdf(x):
            return 0.5 * (1.0 + np.sign(x) * _erf_vec(np.abs(x) / np.sqrt(2)))


def _erf_vec(x):
    x = np.asarray(x, float)
    t = 1.0 / (1.0 + 0.3275911 * x)
    return 1.0 - (((((1.061405429*t - 1.453152027)*t) + 1.421413741)*t
                    - 0.284496736)*t + 0.254829592)*t*np.exp(-x*x)


def _log_gamma(x):
    if x < 0.5:
        return np.log(np.pi / np.sin(np.pi * x)) - _log_gamma(1 - x)
    x -= 1
    a = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
         771.32342877765313, -176.61502916214059, 12.507343278686905,
         -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]
    t = x + 7.5
    return 0.5*np.log(2*np.pi) + (x+0.5)*np.log(t) - t + np.log(
        a[0] + sum(a[i]/(x+i) for i in range(1, 9)))


def _gammainc_lower(a, x):
    x = np.asarray(x, float)
    scalar = x.ndim == 0
    x = np.atleast_1d(x)
    result = np.zeros_like(x)
    for idx in range(len(x)):
        xi = x[idx]
        if xi <= 0:
            continue
        term = xi**a * np.exp(-xi) / max(a, 1e-30)
        s = term
        for n_iter in range(1, 200):
            term *= xi / (a + n_iter)
            s += term
            if abs(term) < 1e-10 * abs(s):
                break
        result[idx] = min(s * np.exp(-_log_gamma(a)), 1.0)
    return float(result[0]) if scalar else result


stats = _scipy_stats if _HAVE_SCIPY else _FallbackStats()


def kaplan_meier_estimator(time, event, group=None):
    """
    Compute Kaplan-Meier survival curves.
    Returns: (times, survival_prob, ci_lower, ci_upper, median_survival)
    """
    if group is not None:
        data = pd.DataFrame({"time": time, "event": event, "group": group})
    else:
        data = pd.DataFrame({"time": time, "event": event})

    # Sort by time
    data = data.sort_values("time").reset_index(drop=True)

    # Get unique event times
    event_times = data[data["event"] == 1]["time"].unique()
    event_times = np.sort(event_times)

    s_t = []
    ci_lower = []
    ci_upper = []
    at_risk_times = []

    cumulative_survival = 1.0
    cumulative_variance = 0.0

    for t in event_times:
        # Number at risk at time t
        n_at_risk = len(data[data["time"] >= t])
        # Number of events at time t
        d_t = len(data[(data["time"] == t) & (data["event"] == 1)])

        # Update survival probability (S(t) = S(t-1) * (1 - d_t / n_at_risk))
        cumulative_survival *= (1.0 - d_t / n_at_risk)

        # Greenwood's formula for variance (guard against n_at_risk == d_t)
        denom = n_at_risk * (n_at_risk - d_t)
        if denom > 0:
            cumulative_variance += d_t / denom

        # Log-log CI
        if cumulative_survival > 0:
            se_log_log_s = np.sqrt(cumulative_variance) / np.log(cumulative_survival)
            log_log_s = np.log(-np.log(max(cumulative_survival, 1e-10)))
            ci_lower_loglog = log_log_s - 1.96 * se_log_log_s
            ci_upper_loglog = log_log_s + 1.96 * se_log_log_s

            lower = np.exp(-np.exp(ci_upper_loglog))
            upper = np.exp(-np.exp(ci_lower_loglog))
        else:
            lower, upper = 0, 0

        s_t.append(cumulative_survival)
        ci_lower.append(max(0, lower))
        ci_upper.append(min(1, upper))
        at_risk_times.append(t)

    # Find median survival time
    median_survival = np.inf
    for i, t in enumerate(at_risk_times):
        if s_t[i] <= 0.5:
            median_survival = t
            break

    return np.array(at_risk_times), np.array(s_t), np.array(ci_lower), np.array(ci_upper), median_survival


def logrank_test(time, event, group):
    """
    Perform log-rank test between groups.
    Returns: (test_statistic, p_value, chi2_stat)
    """
    data = pd.DataFrame({"time": time, "event": event, "group": group})
    data = data.sort_values("time").reset_index(drop=True)

    groups = data["group"].unique()
    event_times = data[data["event"] == 1]["time"].unique()
    event_times = np.sort(event_times)

    # Build contingency table
    observed = defaultdict(float)
    expected = defaultdict(float)

    for t in event_times:
        at_risk_t = data[data["time"] >= t]
        d_t = len(at_risk_t[at_risk_t["event"] == 1])
        n_t = len(at_risk_t)

        for g in groups:
            d_gt = len(at_risk_t[(at_risk_t["group"] == g) & (at_risk_t["event"] == 1)])
            n_gt = len(at_risk_t[at_risk_t["group"] == g])

            observed[g] += d_gt
            expected[g] += (d_t * n_gt) / n_t if n_t > 0 else 0

    # Compute chi-squared statistic
    chi2_stat = 0
    for g in groups:
        if expected[g] > 0:
            chi2_stat += (observed[g] - expected[g]) ** 2 / expected[g]

    df = len(groups) - 1
    p_value = 1 - stats.chi2.cdf(chi2_stat, df)

    return chi2_stat, p_value


def wilcoxon_test(time, event, group):
    """
    Perform Wilcoxon (Peto-Wilcoxon) test between groups.
    Returns: (test_statistic, p_value)
    """
    data = pd.DataFrame({"time": time, "event": event, "group": group})
    data = data.sort_values("time").reset_index(drop=True)

    groups = data["group"].unique()
    event_times = data[data["event"] == 1]["time"].unique()
    event_times = np.sort(event_times)

    observed = defaultdict(float)
    expected = defaultdict(float)
    variance = defaultdict(float)

    for t in event_times:
        at_risk_t = data[data["time"] >= t]
        d_t = len(at_risk_t[at_risk_t["event"] == 1])
        n_t = len(at_risk_t)

        for g in groups:
            d_gt = len(at_risk_t[(at_risk_t["group"] == g) & (at_risk_t["event"] == 1)])
            n_gt = len(at_risk_t[at_risk_t["group"] == g])

            weight = n_gt  # Weight by number at risk
            observed[g] += d_gt * weight
            exp_gt = (d_t * n_gt) / n_t if n_t > 0 else 0
            expected[g] += exp_gt * weight

            if n_t > 1:
                var_gt = (d_t * (n_t - d_t) * n_gt * (n_t - n_gt)) / (n_t**2 * (n_t - 1))
                variance[g] += (weight ** 2) * var_gt

    # Compute test statistic
    test_stat = 0
    total_var = 0
    for g in groups:
        test_stat += (observed[g] - expected[g]) ** 2
        total_var += variance[g]

    if total_var > 0:
        test_stat /= total_var

    p_value = 1 - stats.chi2.cdf(test_stat, len(groups) - 1)

    return test_stat, p_value


def get_at_risk_counts(time, event, group, event_times):
    """
    Get number at risk for each group at given time points.
    """
    data = pd.DataFrame({"time": time, "event": event, "group": group})
    groups = sorted(data["group"].unique())

    at_risk = defaultdict(list)
    for t in event_times:
        for g in groups:
            n = len(data[(data["group"] == g) & (data["time"] >= t)])
            at_risk[g].append(n)

    return groups, at_risk


def main():
    parser = argparse.ArgumentParser(description="Plot Kaplan-Meier survival curves")
    parser.add_argument("--input", required=True, help="Input CSV/TSV file")
    parser.add_argument("--time-col", required=True, help="Time-to-event column name")
    parser.add_argument("--event-col", required=True, help="Event indicator column name")
    parser.add_argument("--group-col", default="", help="Grouping variable column")
    parser.add_argument("--group-order", default="", help="Comma-separated group order")
    parser.add_argument("--group-colors", default="", help="Comma-separated hex colors")
    parser.add_argument("--split-by-median", action="store_true", help="Split continuous variable at median")
    parser.add_argument("--split-col", default="", help="Column to split by median")
    parser.add_argument("--output", required=True, help="Output PNG file")
    parser.add_argument("--output-svg", action="store_true", help="Also save SVG")
    parser.add_argument("--title", default="Survival Curve", help="Plot title")
    parser.add_argument("--xlabel", default="Time (days)", help="X-axis label")
    parser.add_argument("--ylabel", default="Survival probability", help="Y-axis label")
    parser.add_argument("--show-ci", action="store_true", help="Show confidence intervals")
    parser.add_argument("--ci-alpha", type=float, default=0.2, help="CI transparency")
    parser.add_argument("--show-at-risk", action="store_true", help="Show at-risk table")
    parser.add_argument("--show-median", action="store_true", help="Annotate median survival")
    parser.add_argument("--show-pvalue", action="store_true", help="Show p-value")
    parser.add_argument("--test", choices=["logrank", "wilcoxon", "both"], default="logrank")
    parser.add_argument("--pvalue-location", choices=["top_left", "top_right", "bottom_left", "bottom_right"],
                       default="top_right")
    parser.add_argument("--time-unit", default="days", help="Time unit for labels")
    parser.add_argument("--xlim", type=float, default=0, help="Max time to show (0=auto)")
    parser.add_argument("--ylim-min", type=float, default=0, help="Min y-axis value")
    parser.add_argument("--fig-width", type=float, default=8, help="Figure width")
    parser.add_argument("--fig-height", type=float, default=6, help="Figure height")
    parser.add_argument("--dpi", type=int, default=300, help="Output DPI")
    parser.add_argument("--palette", choices=["Set1", "tab10", "custom"], default="Set1")
    parser.add_argument("--font-family", default="sans-serif", help="Font family")
    parser.add_argument("--linewidth", type=float, default=2.0, help="Line width")
    parser.add_argument("--censoring-marks", action="store_true", default=True, help="Show censoring marks")
    parser.add_argument("--censoring-symbol", default="+", help="Censoring marker symbol")

    args = parser.parse_args()
    init_style(
        font_family=getattr(args, 'font_family', None),
        font_size=getattr(args, 'base_fontsize', None),
    )

    # Load data
    print(f"Loading data from {args.input}...")
    input_path = Path(args.input)
    if input_path.suffix == ".tsv":
        df = pd.read_csv(args.input, sep="\t")
    else:
        df = pd.read_csv(args.input)

    # Validate columns
    if args.time_col not in df.columns:
        print(f"Error: Column '{args.time_col}' not found", file=sys.stderr)
        sys.exit(1)
    if args.event_col not in df.columns:
        print(f"Error: Column '{args.event_col}' not found", file=sys.stderr)
        sys.exit(1)

    # Handle group splitting
    if args.split_by_median:
        if not args.split_col:
            print("Error: --split-col required with --split-by-median", file=sys.stderr)
            sys.exit(1)
        if args.split_col not in df.columns:
            print(f"Error: Column '{args.split_col}' not found", file=sys.stderr)
            sys.exit(1)
        median_val = df[args.split_col].median()
        df["_split_group"] = df[args.split_col].apply(lambda x: f"High (>{median_val:.2f})" if x > median_val else f"Low (≤{median_val:.2f})")
        group_col = "_split_group"
    elif args.group_col:
        if args.group_col not in df.columns:
            print(f"Error: Column '{args.group_col}' not found", file=sys.stderr)
            sys.exit(1)
        group_col = args.group_col
    else:
        group_col = None

    # Extract data
    time = df[args.time_col].values
    event = df[args.event_col].values

    if group_col:
        group = df[group_col].values
    else:
        group = np.ones(len(time))

    # Remove missing values
    mask = ~(pd.isna(time) | pd.isna(event))
    time = time[mask]
    event = event[mask]
    group = group[mask]

    print(f"Loaded {len(time)} records, {sum(event)} events")

    # Compute KM curves per group
    unique_groups = np.unique(group)
    if args.group_order:
        ordered = [g.strip() for g in args.group_order.split(",")]
        unique_groups = [g for g in ordered if g in unique_groups]

    km_data = {}
    medians = {}

    for g in unique_groups:
        mask_g = group == g
        times, surv, ci_lower, ci_upper, median = kaplan_meier_estimator(time[mask_g], event[mask_g])
        km_data[g] = {
            "times": times,
            "surv": surv,
            "ci_lower": ci_lower,
            "ci_upper": ci_upper,
            "n": sum(mask_g),
            "n_events": sum(event[mask_g]),
        }
        medians[g] = median

    # Statistical tests
    p_value_lr = None
    p_value_wilcox = None

    if group_col and len(unique_groups) > 1:
        if args.test in ["logrank", "both"]:
            _, p_value_lr = logrank_test(time, event, group)
        if args.test in ["wilcoxon", "both"]:
            _, p_value_wilcox = wilcoxon_test(time, event, group)

    # Set up colors
    if args.palette == "custom" and args.group_colors:
        colors = args.group_colors.split(",")
        color_map = {g: colors[i % len(colors)] for i, g in enumerate(unique_groups)}
    elif args.palette == "Set1":
        cmap = plt.cm.Set1
        color_map = {g: cmap(i % 9) for i, g in enumerate(unique_groups)}
    else:  # tab10
        cmap = plt.cm.tab10
        color_map = {g: cmap(i % 10) for i, g in enumerate(unique_groups)}

    # Create plot
    fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height), dpi=100)
    plt.rcParams["font.family"] = args.font_family

    # Plot KM curves
    max_time = 0
    for g in unique_groups:
        data = km_data[g]
        times = data["times"]
        surv = data["surv"]
        ci_lower = data["ci_lower"]
        ci_upper = data["ci_upper"]
        max_time = max(max_time, np.max(times))

        # Add starting point (t=0, S=1)
        times_plot = np.concatenate([[0], times])
        surv_plot = np.concatenate([[1], surv])
        ci_lower_plot = np.concatenate([[1], ci_lower])
        ci_upper_plot = np.concatenate([[1], ci_upper])

        ax.step(times_plot, surv_plot, where="post", label=f"{g} (n={data['n']})",
                linewidth=args.linewidth, color=color_map[g])

        # Confidence intervals
        if args.show_ci:
            ax.fill_between(times_plot, ci_lower_plot, ci_upper_plot, where=(times_plot <= times_plot[-1]),
                           step="post", alpha=args.ci_alpha, color=color_map[g])

        # Censoring marks
        if args.censoring_marks:
            censor_mask = np.concatenate([[False], data["surv"] > 0])
            # Find censored points (difference in survival = 0)
            censor_times = times[np.concatenate([[False], np.diff(surv) == 0])]
            # Actually, we need to find censored from raw data
            mask_g = group == g
            censor_raw = time[mask_g][event[mask_g] == 0]
            ax.scatter(censor_raw, np.interp(censor_raw, times_plot, surv_plot, left=1, right=surv_plot[-1]),
                      marker=args.censoring_symbol, s=50, color=color_map[g], alpha=0.7)

    # X-axis limit
    if args.xlim > 0:
        ax.set_xlim(0, args.xlim)
    else:
        ax.set_xlim(0, max_time * 1.05)

    ax.set_ylim(args.ylim_min, 1.05)
    ax.set_xlabel(args.xlabel, fontsize=12)
    ax.set_ylabel(args.ylabel, fontsize=12)
    ax.set_title(args.title, fontsize=14, fontweight="bold")
    ax.legend(loc="best", fontsize=10)
    ax.grid(True, alpha=0.3)

    # Add p-value
    if args.show_pvalue and (p_value_lr is not None or p_value_wilcox is not None):
        pval_text = ""
        if p_value_lr is not None:
            pval_text += f"Log-rank p = {p_value_lr:.4f}\n"
        if p_value_wilcox is not None:
            pval_text += f"Wilcoxon p = {p_value_wilcox:.4f}"

        if args.pvalue_location == "top_left":
            ax.text(0.02, 0.98, pval_text, transform=ax.transAxes,
                   fontsize=10, verticalalignment="top", bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
        elif args.pvalue_location == "top_right":
            ax.text(0.98, 0.98, pval_text, transform=ax.transAxes,
                   fontsize=10, verticalalignment="top", horizontalalignment="right",
                   bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
        elif args.pvalue_location == "bottom_left":
            ax.text(0.02, 0.02, pval_text, transform=ax.transAxes,
                   fontsize=10, bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
        else:  # bottom_right
            ax.text(0.98, 0.02, pval_text, transform=ax.transAxes,
                   fontsize=10, horizontalalignment="right",
                   bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

    # Add median survival lines
    if args.show_median:
        for g in unique_groups:
            if medians[g] != np.inf:
                ax.axvline(medians[g], linestyle="--", alpha=0.5, color=color_map[g], linewidth=1)
                ax.axhline(0.5, linestyle="--", alpha=0.3, color="gray", linewidth=0.5)

    # At-risk table
    if args.show_at_risk and group_col:
        # Sample time points for at-risk table
        all_times = np.concatenate([km_data[g]["times"] for g in unique_groups])
        event_times = np.sort(np.unique(all_times))
        step = max(1, len(event_times) // 5)  # Show ~5 time points
        table_times = event_times[::step]

        groups_for_table, at_risk = get_at_risk_counts(time, event, group, table_times)

        # Add text below x-axis
        y_table = -0.25
        for i, t in enumerate(table_times):
            ax.text(t, y_table, f"{t:.0f}", ha="center", fontsize=8, transform=ax.get_xaxis_transform())
            for j, g in enumerate(groups_for_table):
                ax.text(t, y_table - 0.05 - j * 0.03, f"{at_risk[g][i]}", ha="center", fontsize=7,
                       transform=ax.get_xaxis_transform())

    plt.tight_layout()
    plt.savefig(args.output, dpi=args.dpi, bbox_inches="tight")
    print(f"Saved PNG to {args.output}")

    if args.output_svg:
        svg_path = Path(args.output).with_suffix(".svg")
        plt.savefig(svg_path, format="svg", bbox_inches="tight")
        print(f"Saved SVG to {svg_path}")

    plt.close()

    # Write statistics file
    stats_data = []
    for g in unique_groups:
        stats_data.append({
            "group": g,
            "n": km_data[g]["n"],
            "n_events": km_data[g]["n_events"],
            "median_survival": medians[g],
            "ci_lower": "N/A",
            "ci_upper": "N/A",
            "pvalue": p_value_lr if p_value_lr is not None else "N/A",
        })

    stats_df = pd.DataFrame(stats_data)
    stats_file = Path(args.output).with_stem(Path(args.output).stem + "_stats").with_suffix(".tsv")
    stats_df.to_csv(stats_file, sep="\t", index=False)
    print(f"Saved statistics to {stats_file}")

    # Print summary
    print("\n" + "="*50)
    print(f"Survival curve: {args.title}")
    print(f"Total subjects: {len(time)}")
    print(f"Total events: {sum(event)}")
    if p_value_lr:
        print(f"Log-rank p-value: {p_value_lr:.6f}")
    if p_value_wilcox:
        print(f"Wilcoxon p-value: {p_value_wilcox:.6f}")
    print("="*50)


if __name__ == "__main__":
    main()
