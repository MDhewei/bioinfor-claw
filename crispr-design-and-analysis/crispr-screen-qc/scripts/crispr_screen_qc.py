#!/usr/bin/env python3
"""
Quality control analysis for pooled CRISPR screens.

Accepts sgRNA read count matrix and computes library representation, read depth,
guide-level QC metrics, replicate correlations, Gini index, and generates report.
"""

import argparse
import os
import sys
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import sys as _sys, os as _os
try:
    _sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), *(['..'] * 3), '_shared'))
    from plot_style import init_style
except ImportError:
    def init_style(**kw): pass  # graceful fallback if _shared not available
import numpy as np
import pandas as pd
try:
    from scipy.stats import pearsonr, spearmanr
except ImportError:
    def _t_sf_simple(t, df):
        if df <= 0:
            return 0.5
        x = float(df) / (float(df) + float(t) ** 2)
        a, b = df / 2.0, 0.5
        if x <= 0:
            return 0.0
        if x >= 1:
            return 1.0
        lbeta = sum(_lgam(v) for v in [a, b]) - _lgam(a + b)
        front = np.exp(a * np.log(x) + b * np.log(1 - x) - lbeta) / a
        f = 1.0
        C, D = 1.0, 1.0 - (a + b) * x / (a + 1)
        if abs(D) < 1e-30:
            D = 1e-30
        D = 1.0 / D
        f = D
        for m in range(1, 100):
            for step in (1, 2):
                num = (m * (b - m) * x / ((a + 2*m - 1) * (a + 2*m)) if step == 1
                       else -(a + m) * (a + b + m) * x / ((a + 2*m) * (a + 2*m + 1)))
                D = 1.0 / max(abs(1.0 + num * D), 1e-30) * (1 if 1 + num * D >= 0 else -1)
                C = 1.0 + num / max(abs(C), 1e-30) * (1 if C >= 0 else -1)
                delta = C * D
                f *= delta
                if abs(delta - 1.0) < 1e-10:
                    break
        return min(front * f, 1.0)

    def _lgam(x):
        if x < 0.5:
            return np.log(np.pi / np.sin(np.pi * x)) - _lgam(1 - x)
        x -= 1
        a = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
             771.32342877765313, -176.61502916214059, 12.507343278686905,
             -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]
        t = x + 7.5
        return 0.5*np.log(2*np.pi) + (x+0.5)*np.log(t) - t + np.log(
            a[0] + sum(a[i]/(x+i) for i in range(1, 9)))

    def pearsonr(a, b):
        a = np.asarray(a, float) - np.mean(a)
        b = np.asarray(b, float) - np.mean(b)
        r = float((a * b).sum() / (np.sqrt((a**2).sum() * (b**2).sum()) + 1e-15))
        n = len(a)
        t = r * np.sqrt(max(n - 2, 1)) / np.sqrt(max(1 - r**2, 1e-15))
        return r, min(2 * _t_sf_simple(abs(t), n - 2), 1.0)

    def spearmanr(a, b=None):
        if b is None:
            a, b = np.asarray(a)[:, 0], np.asarray(a)[:, 1]
        def _rank(arr):
            order = np.argsort(arr)
            r = np.empty(len(arr), float)
            r[order] = np.arange(1, len(arr) + 1)
            return r
        ra = _rank(np.asarray(a, float))
        rb = _rank(np.asarray(b, float))
        return pearsonr(ra, rb)


def calculate_gini(counts: np.ndarray) -> float:
    """
    Calculate Gini index from count array.
    Gini = 0 (uniform) to 1 (maximally skewed).
    """
    counts = np.asarray(counts, dtype=float)
    counts = counts[counts > 0]  # Remove zeros

    if len(counts) == 0:
        return np.nan

    sorted_counts = np.sort(counts)
    n = len(sorted_counts)
    numerator = 2 * np.sum(np.arange(1, n + 1) * sorted_counts)
    denominator = n * np.sum(sorted_counts)

    gini = (numerator / denominator) - (n + 1) / n
    return max(0.0, gini)  # Clamp to [0, 1)


def normalize_to_cpm(counts: np.ndarray) -> np.ndarray:
    """Normalize counts to CPM (counts per million)."""
    total = np.sum(counts)
    if total == 0:
        return counts
    return (counts / total) * 1e6


def normalize_log2_cpm(counts: np.ndarray, pseudocount: float = 1.0) -> np.ndarray:
    """Normalize to log2(CPM + pseudocount)."""
    cpm = normalize_to_cpm(counts)
    return np.log2(cpm + pseudocount)


def load_count_matrix(filepath: str, sgrna_col: str, gene_col: str, samples: List[str] = None) -> Tuple[pd.DataFrame, List[str]]:
    """Load count matrix TSV."""
    df = pd.read_csv(filepath, sep="\t")

    # Validate columns
    if sgrna_col not in df.columns:
        raise ValueError(f"Column '{sgrna_col}' not found in {filepath}")
    if gene_col not in df.columns:
        raise ValueError(f"Column '{gene_col}' not found in {filepath}")

    # Identify sample columns (all except sgrna_col and gene_col)
    sample_cols = [c for c in df.columns if c not in [sgrna_col, gene_col]]

    if samples:
        sample_cols = [c for c in sample_cols if c in samples]

    if not sample_cols:
        raise ValueError("No valid sample columns found")

    # Filter to selected samples
    df = df[[sgrna_col, gene_col] + sample_cols].copy()

    # Remove all-zero rows
    count_cols = sample_cols
    df = df.loc[(df[count_cols] > 0).any(axis=1)]

    return df, sample_cols


def compute_sample_qc(counts: pd.Series, min_reads: int = 30) -> Dict:
    """Compute QC metrics for a single sample."""
    counts_arr = counts.values.astype(float)

    total_reads = np.sum(counts_arr)
    num_guides = len(counts_arr)
    num_represented = np.sum(counts_arr >= min_reads)
    representation = num_represented / num_guides if num_guides > 0 else 0.0

    median_reads = np.median(counts_arr[counts_arr > 0]) if np.any(counts_arr > 0) else 0.0
    mean_reads = np.mean(counts_arr[counts_arr > 0]) if np.any(counts_arr > 0) else 0.0
    min_reads_val = np.min(counts_arr[counts_arr > 0]) if np.any(counts_arr > 0) else 0.0
    max_reads_val = np.max(counts_arr)

    num_missing = np.sum(counts_arr == 0)
    missing_fraction = num_missing / num_guides if num_guides > 0 else 0.0

    gini = calculate_gini(counts_arr)

    return {
        "total_reads": int(total_reads),
        "num_guides": num_guides,
        "represented_guides": num_represented,
        "representation_fraction": representation,
        "median_reads": median_reads,
        "mean_reads": mean_reads,
        "min_reads": min_reads_val,
        "max_reads": max_reads_val,
        "missing_guides": num_missing,
        "missing_fraction": missing_fraction,
        "gini_index": gini,
    }


def compute_control_stats(
    df: pd.DataFrame, gene_col: str, sample_cols: List[str], control_prefix: str
) -> Dict:
    """Compute statistics for control vs targeting guides."""
    control_guides = df[df[gene_col].str.startswith(control_prefix, na=False)]
    targeting_guides = df[~df[gene_col].str.startswith(control_prefix, na=False)]

    results = {}

    for sample in sample_cols:
        ctrl_counts = control_guides[sample].values.astype(float)
        target_counts = targeting_guides[sample].values.astype(float)

        ctrl_mean = np.mean(ctrl_counts[ctrl_counts > 0]) if np.any(ctrl_counts > 0) else 0.0
        target_mean = np.mean(target_counts[target_counts > 0]) if np.any(target_counts > 0) else 0.0

        ctrl_cv = np.std(ctrl_counts) / ctrl_mean if ctrl_mean > 0 else np.nan
        target_cv = np.std(target_counts) / target_mean if target_mean > 0 else np.nan

        results[sample] = {
            "control_mean": ctrl_mean,
            "targeting_mean": target_mean,
            "control_cv": ctrl_cv,
            "targeting_cv": target_cv,
            "num_controls": len(ctrl_counts),
            "num_targeting": len(target_counts),
        }

    return results


def compute_guide_cv(df: pd.DataFrame, sgrna_col: str, sample_cols: List[str]) -> pd.DataFrame:
    """Compute per-guide CV and basic stats across samples."""
    counts = df[sample_cols].values.astype(float)

    guide_stats = pd.DataFrame({
        sgrna_col: df[sgrna_col],
        "mean_reads": np.mean(counts, axis=1),
        "cv": np.std(counts, axis=1) / (np.mean(counts, axis=1) + 1e-6),
        "min_reads": np.min(counts, axis=1),
        "max_reads": np.max(counts, axis=1),
    })

    return guide_stats


def compute_correlations(df: pd.DataFrame, sample_cols: List[str], log_norm: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    """Compute Pearson and Spearman correlations between all sample pairs."""
    counts = df[sample_cols].values.astype(float)

    if log_norm:
        counts = np.log2(normalize_to_cpm(counts.T).T + 1.0)

    num_samples = len(sample_cols)
    pearson_corr = np.ones((num_samples, num_samples))
    spearman_corr = np.ones((num_samples, num_samples))

    for i in range(num_samples):
        for j in range(i + 1, num_samples):
            p_corr, _ = pearsonr(counts[:, i], counts[:, j])
            s_corr, _ = spearmanr(counts[:, i], counts[:, j])

            pearson_corr[i, j] = p_corr
            pearson_corr[j, i] = p_corr
            spearman_corr[i, j] = s_corr
            spearman_corr[j, i] = s_corr

    return pearson_corr, spearman_corr


def plot_distributions(df: pd.DataFrame, sample_cols: List[str], outdir: str):
    """Plot read count distributions per sample."""
    counts = df[sample_cols].values.astype(float)

    fig, ax = plt.subplots(figsize=(12, 6))

    data_to_plot = [counts[:, i] for i in range(len(sample_cols))]
    bp = ax.boxplot(data_to_plot, labels=sample_cols, patch_artist=True)

    for patch in bp["boxes"]:
        patch.set_facecolor("lightblue")

    ax.set_ylabel("Read Count")
    ax.set_title("Read Count Distribution per Sample")
    ax.set_yscale("log")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "read_distributions.png"), dpi=300, bbox_inches="tight")
    plt.close()


def plot_gini(gini_dict: Dict[str, float], outdir: str):
    """Plot Gini index per sample."""
    samples = list(gini_dict.keys())
    ginis = list(gini_dict.values())

    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar(samples, ginis, color="coral", edgecolor="black")

    ax.axhline(0.25, color="green", linestyle="--", label="Good (<0.25)")
    ax.axhline(0.35, color="red", linestyle="--", label="Poor (>0.35)")

    ax.set_ylabel("Gini Index")
    ax.set_title("Gini Index (Guide Representation Evenness)")
    ax.set_ylim(0, max(ginis) * 1.1 if ginis else 1.0)
    ax.legend()
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "gini_index.png"), dpi=300, bbox_inches="tight")
    plt.close()


def plot_cdf(df: pd.DataFrame, sample_cols: List[str], outdir: str):
    """Plot CDF of reads per guide per sample."""
    fig, ax = plt.subplots(figsize=(12, 6))

    for sample in sample_cols:
        counts = np.sort(df[sample].values.astype(float))
        cdf = np.arange(1, len(counts) + 1) / len(counts)
        ax.plot(counts, cdf, label=sample, linewidth=2)

    ax.set_xlabel("Reads per Guide")
    ax.set_ylabel("Cumulative Fraction")
    ax.set_title("CDF of Read Counts per Guide")
    ax.set_xscale("log")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "cdf_reads.png"), dpi=300, bbox_inches="tight")
    plt.close()


def plot_correlation_heatmap(corr: np.ndarray, sample_cols: List[str], outdir: str, name: str = "pearson"):
    """Plot correlation heatmap."""
    fig, ax = plt.subplots(figsize=(10, 10))

    im = ax.imshow(corr, cmap="coolwarm", vmin=-1, vmax=1)
    ax.set_xticks(np.arange(len(sample_cols)))
    ax.set_yticks(np.arange(len(sample_cols)))
    ax.set_xticklabels(sample_cols, rotation=45, ha="right")
    ax.set_yticklabels(sample_cols)

    # Add text annotations
    for i in range(len(sample_cols)):
        for j in range(len(sample_cols)):
            text = ax.text(j, i, f"{corr[i, j]:.2f}", ha="center", va="center", color="black", fontsize=9)

    ax.set_title(f"{name.capitalize()} Correlation Heatmap")
    plt.colorbar(im, ax=ax)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"correlation_{name}.png"), dpi=300, bbox_inches="tight")
    plt.close()


def plot_control_vs_targeting(
    df: pd.DataFrame, gene_col: str, sample_cols: List[str], control_prefix: str, outdir: str
):
    """Plot control vs targeting guide distributions."""
    fig, axes = plt.subplots(1, len(sample_cols), figsize=(4 * len(sample_cols), 5))

    if len(sample_cols) == 1:
        axes = [axes]

    for idx, sample in enumerate(sample_cols):
        ctrl_counts = df[df[gene_col].str.startswith(control_prefix, na=False)][sample].values.astype(float)
        target_counts = df[~df[gene_col].str.startswith(control_prefix, na=False)][sample].values.astype(float)

        axes[idx].hist([target_counts[target_counts > 0], ctrl_counts[ctrl_counts > 0]],
                       bins=20, label=["Targeting", "Control"], color=["steelblue", "coral"])
        axes[idx].set_xlabel("Reads")
        axes[idx].set_ylabel("Count")
        axes[idx].set_title(f"Control vs Targeting\n({sample})")
        axes[idx].set_yscale("log")
        axes[idx].legend()

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "control_vs_targeting.png"), dpi=300, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="CRISPR screen QC analysis")
    parser.add_argument("--counts", type=str, required=True, help="TSV with counts matrix")
    parser.add_argument("--sgrna-col", type=str, default="sgRNA", help="sgRNA ID column")
    parser.add_argument("--gene-col", type=str, default="Gene", help="Gene column")
    parser.add_argument("--outdir", type=str, default="qc_results")
    parser.add_argument("--min-reads", type=int, default=30)
    parser.add_argument("--min-representation", type=float, default=0.8)
    parser.add_argument("--n-sgrnas-per-gene", type=int, default=4)
    parser.add_argument("--control-prefix", type=str, default="non")
    parser.add_argument("--samples", type=str, help="Comma-separated sample names to include")
    parser.add_argument("--log2-norm", action="store_true", help="Apply log2 normalization")

    args = parser.parse_args()

    init_style(
        font_family=getattr(args, 'font_family', None),
        font_size=getattr(args, 'font_size', None),
    )

    # Create output directory
    os.makedirs(args.outdir, exist_ok=True)

    # Parse samples
    sample_list = args.samples.split(",") if args.samples else None

    # Load data
    try:
        df, sample_cols = load_count_matrix(args.counts, args.sgrna_col, args.gene_col, sample_list)
    except Exception as e:
        print(f"Error loading count matrix: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Loaded {len(df)} guides, {len(sample_cols)} samples")

    # Compute sample-level QC
    sample_qc = {}
    for sample in sample_cols:
        sample_qc[sample] = compute_sample_qc(df[sample], args.min_reads)

    # Compute control stats
    control_stats = compute_control_stats(df, args.gene_col, sample_cols, args.control_prefix)

    # Save QC summary
    qc_summary = pd.DataFrame([
        {
            "sample": sample,
            "total_reads": sample_qc[sample]["total_reads"],
            "num_guides": sample_qc[sample]["num_guides"],
            "represented_guides": sample_qc[sample]["represented_guides"],
            "representation_fraction": round(sample_qc[sample]["representation_fraction"], 3),
            "median_reads": round(sample_qc[sample]["median_reads"], 1),
            "mean_reads": round(sample_qc[sample]["mean_reads"], 1),
            "gini_index": round(sample_qc[sample]["gini_index"], 3),
            "missing_guides": sample_qc[sample]["missing_guides"],
            "missing_fraction": round(sample_qc[sample]["missing_fraction"], 3),
            "control_mean": round(control_stats[sample]["control_mean"], 1),
            "targeting_mean": round(control_stats[sample]["targeting_mean"], 1),
        }
        for sample in sample_cols
    ])

    qc_summary.to_csv(os.path.join(args.outdir, "qc_summary.tsv"), sep="\t", index=False)
    print(f"Saved QC summary to {args.outdir}/qc_summary.tsv")
    print(qc_summary)

    # Compute and save guide-level QC
    guide_qc = compute_guide_cv(df, args.sgrna_col, sample_cols)
    guide_qc.to_csv(os.path.join(args.outdir, "guide_qc.tsv"), sep="\t", index=False)
    print(f"Saved guide-level QC to {args.outdir}/guide_qc.tsv")

    # Compute correlations
    pearson_corr, spearman_corr = compute_correlations(df, sample_cols, args.log2_norm)

    # Plots
    plot_distributions(df, sample_cols, args.outdir)
    print(f"Saved distribution plot to {args.outdir}/read_distributions.png")

    gini_dict = {sample: sample_qc[sample]["gini_index"] for sample in sample_cols}
    plot_gini(gini_dict, args.outdir)
    print(f"Saved Gini plot to {args.outdir}/gini_index.png")

    plot_cdf(df, sample_cols, args.outdir)
    print(f"Saved CDF plot to {args.outdir}/cdf_reads.png")

    plot_correlation_heatmap(pearson_corr, sample_cols, args.outdir, "pearson")
    print(f"Saved correlation heatmap to {args.outdir}/correlation_pearson.png")

    plot_control_vs_targeting(df, args.gene_col, sample_cols, args.control_prefix, args.outdir)
    print(f"Saved control vs targeting plot to {args.outdir}/control_vs_targeting.png")

    # Summary report
    print("\n=== QC Summary ===")
    for sample in sample_cols:
        rep_frac = sample_qc[sample]["representation_fraction"]
        gini = sample_qc[sample]["gini_index"]
        status = "PASS" if rep_frac >= args.min_representation and gini < 0.35 else "WARN"
        print(f"{sample}: {status} (representation={rep_frac:.1%}, gini={gini:.3f})")


if __name__ == "__main__":
    main()
