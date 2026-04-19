#!/usr/bin/env python3
"""
RNA-seq Differential Expression Analysis

Performs differential expression analysis between two groups using
a pseudo-bulk approach with statistical testing. Works with raw
count matrices (genes × samples) in TSV or CSV format.

Statistical methods supported:
  - DESeq2-style (via pydeseq2) — recommended for count data
  - t-test + log2FC                — fast fallback for normalized data
  - Mann-Whitney U                 — non-parametric alternative

Outputs:
  - Full DE results table (TSV)
  - Significant genes table (TSV)
  - Volcano plot (PNG + PDF)
  - MA plot (PNG + PDF)
  - Top-genes heatmap (PNG + PDF)
  - Summary statistics (TSV + TXT)

Dependencies: pandas, numpy, scipy, matplotlib, seaborn
Optional:     pydeseq2 (for DESeq2-style analysis)
"""

import argparse
import os
import sys
import warnings
from typing import Dict, List, Optional, Tuple

import matplotlib
import sys as _sys, os as _os
try:
    _sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), *(['..'] * 3), '_shared'))
    from plot_style import init_style
except ImportError:
    def init_style(**kw): pass  # graceful fallback if _shared not available
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

try:
    from scipy import stats as _scipy_stats
    _HAVE_SCIPY = True
except ImportError:
    _scipy_stats = None
    _HAVE_SCIPY = False


# ---------------------------------------------------------------------------
# Pure-numpy fallbacks for scipy.stats
# ---------------------------------------------------------------------------
def _erf_approx(x):
    """Vectorised Abramowitz & Stegun erf approximation."""
    x = np.asarray(x, float)
    t = 1.0 / (1.0 + 0.3275911 * np.abs(x))
    poly = (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t
              - 0.284496736) * t + 0.254829592) * t
    y = 1.0 - poly * np.exp(-(x ** 2))
    return np.sign(x) * y


def _log_gamma(x):
    """log Γ(x) via Lanczos."""
    if x < 0.5:
        return np.log(np.pi / np.sin(np.pi * x)) - _log_gamma(1 - x)
    x -= 1
    a = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
         771.32342877765313, -176.61502916214059, 12.507343278686905,
         -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]
    t = x + 7.5
    return (0.5 * np.log(2 * np.pi) + (x + 0.5) * np.log(t) - t
            + np.log(a[0] + sum(a[i] / (x + i) for i in range(1, 9))))


def _betainc_cf(a, b, x):
    """Regularised incomplete beta I_x(a,b) via Lentz continued fraction."""
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    lbeta = _log_gamma(a) + _log_gamma(b) - _log_gamma(a + b)
    front = np.exp(a * np.log(x) + b * np.log(1 - x) - lbeta) / a
    f = 1.0
    C = 1.0
    D = 1.0 - (a + b) * x / (a + 1)
    if abs(D) < 1e-30:
        D = 1e-30
    D = 1.0 / D
    f = D
    for m in range(1, 100):
        for step in (1, 2):
            if step == 1:
                num = m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m))
            else:
                num = -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1))
            D = 1.0 + num * D
            if abs(D) < 1e-30:
                D = 1e-30
            D = 1.0 / D
            C = 1.0 + num / C
            if abs(C) < 1e-30:
                C = 1e-30
            delta = C * D
            f *= delta
            if abs(delta - 1.0) < 1e-10:
                break
    return min(front * f, 1.0)


def _t_sf(t_val, df):
    """P(T > t_val) for Student-t with df degrees of freedom."""
    if df <= 0:
        return 0.5
    x = float(df) / (float(df) + float(t_val) ** 2)
    return 0.5 * _betainc_cf(df / 2.0, 0.5, x)


def _gammainc_lower(a, x):
    """Regularised lower incomplete gamma P(a, x) = 1 - Q(a, x)."""
    x = np.asarray(x, float)
    scalar = x.ndim == 0
    x = np.atleast_1d(x)
    result = np.zeros_like(x)
    for i, xi in enumerate(x):
        if xi <= 0:
            continue
        term = float(xi) ** a * np.exp(-xi) / max(a, 1e-30)
        s = term
        for n_iter in range(1, 300):
            term *= xi / (a + n_iter)
            s += term
            if abs(term) < 1e-12 * abs(s):
                break
        result[i] = min(s * np.exp(-_log_gamma(a)), 1.0)
    return float(result[0]) if scalar else result


def _welch_ttest(a, b):
    """Two-sample Welch t-test; returns (t_stat, p_value)."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return np.nan, np.nan
    va = np.var(a, ddof=1)
    vb = np.var(b, ddof=1)
    se = np.sqrt(va / na + vb / nb)
    if se == 0:
        return 0.0, 1.0
    t = float((a.mean() - b.mean()) / se)
    df = (va / na + vb / nb) ** 2 / ((va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1))
    p = 2 * _t_sf(abs(t), df)
    return t, float(min(p, 1.0))


def _mannwhitneyu_fallback(a, b, alternative='two-sided'):
    """Mann-Whitney U; returns (U, p_value)."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    na, nb = len(a), len(b)
    if na == 0 or nb == 0:
        return np.nan, np.nan
    combined = np.concatenate([a, b])
    order = np.argsort(combined, kind='stable')
    ranks = np.empty(len(combined))
    ranks[order] = np.arange(1, len(combined) + 1)
    # handle ties with midranks
    sorted_c = combined[order]
    i = 0
    while i < len(sorted_c):
        j = i
        while j < len(sorted_c) - 1 and sorted_c[j + 1] == sorted_c[j]:
            j += 1
        if j > i:
            midrank = (ranks[order[i]] + ranks[order[j]]) / 2
            for k in range(i, j + 1):
                ranks[order[k]] = midrank
        i = j + 1
    U = np.sum(ranks[:na]) - na * (na + 1) / 2.0
    mu = na * nb / 2.0
    sigma = np.sqrt(na * nb * (na + nb + 1) / 12.0)
    z = (U - mu) / (sigma + 1e-15)
    p = float(2 * (1 - 0.5 * (1 + _erf_approx(abs(z) / np.sqrt(2)))))
    return float(U), min(p, 1.0)


class _FallbackStats:
    """Drop-in replacements for scipy.stats functions used in this script."""

    @staticmethod
    def ttest_ind(a, b, equal_var=True, nan_policy='propagate'):
        return _welch_ttest(a, b)

    @staticmethod
    def mannwhitneyu(a, b, alternative='two-sided'):
        return _mannwhitneyu_fallback(a, b, alternative)

    @staticmethod
    def ranksums(a, b):
        return _mannwhitneyu_fallback(a, b, 'two-sided')

    @staticmethod
    def spearmanr(a, b=None):
        if b is None:
            a, b = a[:, 0], a[:, 1]
        a = np.asarray(a, float)
        b = np.asarray(b, float)

        def _rank(arr):
            order = np.argsort(arr)
            r = np.empty(len(arr), float)
            r[order] = np.arange(1, len(arr) + 1)
            return r

        ra = _rank(a) - (len(a) + 1) / 2.0
        rb = _rank(b) - (len(b) + 1) / 2.0
        denom = np.sqrt((ra ** 2).sum() * (rb ** 2).sum()) + 1e-15
        r = float((ra * rb).sum() / denom)
        n = len(a)
        t = r * np.sqrt(max(n - 2, 1)) / np.sqrt(max(1 - r ** 2, 1e-15))
        p = 2 * _t_sf(abs(t), n - 2)
        return r, float(min(p, 1.0))

    @staticmethod
    def pearsonr(a, b):
        a = np.asarray(a, float) - np.mean(a)
        b = np.asarray(b, float) - np.mean(b)
        denom = np.sqrt((a ** 2).sum() * (b ** 2).sum()) + 1e-15
        r = float((a * b).sum() / denom)
        n = len(a)
        t = r * np.sqrt(max(n - 2, 1)) / np.sqrt(max(1 - r ** 2, 1e-15))
        p = 2 * _t_sf(abs(t), n - 2)
        return r, float(min(p, 1.0))

    class norm:
        @staticmethod
        def cdf(x):
            x = np.asarray(x, float)
            return 0.5 * (1.0 + _erf_approx(x / np.sqrt(2)))

    class chi2:
        @staticmethod
        def cdf(x, df):
            return _gammainc_lower(df / 2.0, np.asarray(x, float) / 2.0)


stats = _scipy_stats if _HAVE_SCIPY else _FallbackStats()


warnings.filterwarnings("ignore", category=RuntimeWarning)

# =========================================================
# Utilities
# =========================================================

def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def load_count_matrix(path: str) -> pd.DataFrame:
    """
    Load a count/expression matrix (genes × samples).
    Accepts TSV or CSV. First column is treated as gene identifiers.
    """
    sep = "\t" if path.endswith(".tsv") or path.endswith(".txt") else ","
    df = pd.read_csv(path, sep=sep, index_col=0)
    df.index.name = "gene"
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(how="all")
    return df


def load_sample_metadata(path: str) -> pd.DataFrame:
    """
    Load sample metadata table.
    Must contain 'sample' (or 'sample_id') and 'group' columns.
    """
    sep = "\t" if path.endswith(".tsv") or path.endswith(".txt") else ","
    meta = pd.read_csv(path, sep=sep)
    meta.columns = meta.columns.str.lower()
    # Accept 'sample_id' as alias for 'sample'
    if "sample" not in meta.columns and "sample_id" in meta.columns:
        meta = meta.rename(columns={"sample_id": "sample"})
    required = {"sample", "group"}
    missing = required - set(meta.columns)
    if missing:
        raise ValueError(
            f"Metadata file must contain columns: {required}. "
            f"Missing: {missing}. Found: {list(meta.columns)}"
        )
    return meta


def filter_low_expression(
    counts: pd.DataFrame,
    min_count: int = 10,
    min_samples: int = 2,
) -> pd.DataFrame:
    """Remove genes with low counts across samples."""
    mask = (counts >= min_count).sum(axis=1) >= min_samples
    filtered = counts[mask]
    return filtered


def log2_normalize(counts: pd.DataFrame) -> pd.DataFrame:
    """CPM normalization followed by log2(CPM + 1)."""
    lib_sizes = counts.sum(axis=0)
    cpm = counts.div(lib_sizes, axis=1) * 1e6
    return np.log2(cpm + 1)


# =========================================================
# Statistical testing
# =========================================================

def run_ttest(
    counts: pd.DataFrame,
    group_a_samples: List[str],
    group_b_samples: List[str],
    use_log2norm: bool = True,
) -> pd.DataFrame:
    """
    Differential expression via Welch's t-test on log2-normalized counts.
    Returns DataFrame with: gene, log2FC, mean_a, mean_b, pvalue, padj.
    """
    expr = log2_normalize(counts) if use_log2norm else counts.copy()

    a = expr[group_a_samples]
    b = expr[group_b_samples]

    rows = []
    for gene in expr.index:
        va = a.loc[gene].dropna().values
        vb = b.loc[gene].dropna().values
        if len(va) < 2 or len(vb) < 2:
            continue
        _, pval = stats.ttest_ind(vb, va, equal_var=False)
        mean_a  = float(np.mean(va))
        mean_b  = float(np.mean(vb))
        log2fc  = mean_b - mean_a  # B vs A (condition vs control)
        rows.append(
            {
                "gene":   gene,
                "log2FC": log2fc,
                "mean_control": mean_a,
                "mean_condition": mean_b,
                "pvalue": float(pval),
            }
        )

    result = pd.DataFrame(rows)
    if result.empty:
        return result

    # Benjamini-Hochberg FDR
    result = result.sort_values("pvalue").reset_index(drop=True)
    n = len(result)
    result["padj"] = np.minimum(
        1.0,
        result["pvalue"] * n / (np.arange(1, n + 1)),
    )
    # Enforce monotonicity (step-up)
    result["padj"] = result["padj"][::-1].cummin()[::-1]
    result = result.sort_values("log2FC", ascending=False).reset_index(drop=True)
    return result


def run_mannwhitney(
    counts: pd.DataFrame,
    group_a_samples: List[str],
    group_b_samples: List[str],
    use_log2norm: bool = True,
) -> pd.DataFrame:
    """Differential expression via Mann-Whitney U test."""
    expr = log2_normalize(counts) if use_log2norm else counts.copy()

    a = expr[group_a_samples]
    b = expr[group_b_samples]

    rows = []
    for gene in expr.index:
        va = a.loc[gene].dropna().values
        vb = b.loc[gene].dropna().values
        if len(va) < 2 or len(vb) < 2:
            continue
        _, pval = stats.mannwhitneyu(vb, va, alternative="two-sided")
        mean_a = float(np.mean(va))
        mean_b = float(np.mean(vb))
        rows.append(
            {
                "gene":          gene,
                "log2FC":        mean_b - mean_a,
                "mean_control":  mean_a,
                "mean_condition": mean_b,
                "pvalue":        float(pval),
            }
        )

    result = pd.DataFrame(rows)
    if result.empty:
        return result

    result = result.sort_values("pvalue").reset_index(drop=True)
    n = len(result)
    result["padj"] = np.minimum(
        1.0,
        result["pvalue"] * n / (np.arange(1, n + 1)),
    )
    result["padj"] = result["padj"][::-1].cummin()[::-1]
    result = result.sort_values("log2FC", ascending=False).reset_index(drop=True)
    return result


def run_pydeseq2(
    counts: pd.DataFrame,
    group_a_samples: List[str],
    group_b_samples: List[str],
    group_a_label: str,
    group_b_label: str,
) -> pd.DataFrame:
    """
    Differential expression via PyDESeq2 (DESeq2 Python port).
    Requires integer raw counts. Falls back to t-test if unavailable.
    """
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError:
        print("[WARN] pydeseq2 not installed. Falling back to t-test.")
        return run_ttest(counts, group_a_samples, group_b_samples)

    # PyDESeq2 requires integer count matrix (samples × genes)
    sub = counts[group_a_samples + group_b_samples].copy()
    sub = sub.fillna(0).astype(int)
    sub_T = sub.T  # samples × genes

    meta = pd.DataFrame(
        {"condition": [group_a_label] * len(group_a_samples) + [group_b_label] * len(group_b_samples)},
        index=group_a_samples + group_b_samples,
    )

    try:
        dds = DeseqDataSet(
            counts=sub_T,
            metadata=meta,
            design_factors="condition",
            quiet=True,
        )
        dds.deseq2()
        ds = DeseqStats(dds, contrast=["condition", group_b_label, group_a_label], quiet=True)
        ds.summary()
        res = ds.results_df.reset_index()
        res = res.rename(columns={
            "index":    "gene",
            "log2FoldChange": "log2FC",
            "pvalue":   "pvalue",
            "padj":     "padj",
            "baseMean": "base_mean",
            "lfcSE":    "lfcSE",
            "stat":     "stat",
        })
        res["mean_control"]   = None
        res["mean_condition"] = None
        return res[["gene", "log2FC", "mean_control", "mean_condition", "pvalue", "padj",
                    "base_mean", "lfcSE", "stat"]].sort_values("padj").reset_index(drop=True)

    except Exception as ex:
        print(f"[WARN] PyDESeq2 failed ({ex}). Falling back to t-test.")
        return run_ttest(counts, group_a_samples, group_b_samples)


# =========================================================
# Plots
# =========================================================

def plot_volcano(
    result: pd.DataFrame,
    outdir: str,
    prefix: str,
    fc_thresh: float,
    padj_thresh: float,
    top_n_label: int = 15,
) -> None:
    """Volcano plot: log2FC vs -log10(padj)."""
    df = result.copy().dropna(subset=["log2FC", "padj"])
    df["neg_log10_padj"] = -np.log10(df["padj"].clip(lower=1e-300))

    up   = (df["log2FC"] >= fc_thresh)  & (df["padj"] <= padj_thresh)
    down = (df["log2FC"] <= -fc_thresh) & (df["padj"] <= padj_thresh)
    ns   = ~(up | down)

    fig, ax = plt.subplots(figsize=(7, 6), dpi=300)
    ax.scatter(df.loc[ns,   "log2FC"], df.loc[ns,   "neg_log10_padj"], s=8, color="#AAAAAA", alpha=0.5, linewidths=0)
    ax.scatter(df.loc[up,   "log2FC"], df.loc[up,   "neg_log10_padj"], s=10, color="#D62728", alpha=0.8, linewidths=0, label=f"Up ({up.sum()})")
    ax.scatter(df.loc[down, "log2FC"], df.loc[down, "neg_log10_padj"], s=10, color="#1F77B4", alpha=0.8, linewidths=0, label=f"Down ({down.sum()})")

    # Threshold lines
    ax.axvline( fc_thresh, color="gray", lw=0.8, ls="--", alpha=0.7)
    ax.axvline(-fc_thresh, color="gray", lw=0.8, ls="--", alpha=0.7)
    ax.axhline(-np.log10(padj_thresh), color="gray", lw=0.8, ls="--", alpha=0.7)

    # Label top genes by significance
    top = df[up | down].nlargest(top_n_label, "neg_log10_padj")
    for _, row in top.iterrows():
        ax.text(row["log2FC"], row["neg_log10_padj"] + 0.15,
                row["gene"], fontsize=6.5, ha="center", va="bottom",
                color="black", alpha=0.85)

    ax.set_xlabel("log₂ Fold Change", fontsize=12)
    ax.set_ylabel("-log₁₀(adjusted p-value)", fontsize=12)
    ax.set_title(f"Volcano Plot  ·  {prefix}", fontsize=13)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(fontsize=9, loc="upper left")
    plt.tight_layout()

    base = os.path.join(outdir, f"{prefix}.volcano")
    fig.savefig(f"{base}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{base}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_ma(
    result: pd.DataFrame,
    outdir: str,
    prefix: str,
    fc_thresh: float,
    padj_thresh: float,
) -> None:
    """MA plot: mean expression (A) vs log2FC (M)."""
    df = result.copy().dropna(subset=["log2FC", "padj"])
    # A = average of mean_control and mean_condition; fall back to base_mean
    if "mean_control" in df.columns and df["mean_control"].notna().any():
        df["A"] = (df["mean_control"].fillna(0) + df["mean_condition"].fillna(0)) / 2
    elif "base_mean" in df.columns:
        df["A"] = df["base_mean"]
    else:
        df["A"] = 0

    sig = (df["padj"] <= padj_thresh) & (df["log2FC"].abs() >= fc_thresh)

    fig, ax = plt.subplots(figsize=(7, 5), dpi=300)
    ax.scatter(df.loc[~sig, "A"], df.loc[~sig, "log2FC"], s=5, color="#AAAAAA", alpha=0.4, linewidths=0)
    ax.scatter(df.loc[sig,  "A"], df.loc[sig,  "log2FC"], s=8, color="#D62728", alpha=0.8, linewidths=0, label=f"Significant ({sig.sum()})")
    ax.axhline(0, color="black", lw=0.8)
    ax.axhline( fc_thresh, color="gray", lw=0.7, ls="--", alpha=0.7)
    ax.axhline(-fc_thresh, color="gray", lw=0.7, ls="--", alpha=0.7)

    ax.set_xlabel("Mean expression (log₂ normalized)", fontsize=11)
    ax.set_ylabel("log₂ Fold Change", fontsize=11)
    ax.set_title(f"MA Plot  ·  {prefix}", fontsize=13)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(fontsize=9)
    plt.tight_layout()

    base = os.path.join(outdir, f"{prefix}.ma_plot")
    fig.savefig(f"{base}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{base}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_heatmap(
    counts: pd.DataFrame,
    sig_genes: List[str],
    group_a_samples: List[str],
    group_b_samples: List[str],
    group_a_label: str,
    group_b_label: str,
    outdir: str,
    prefix: str,
    top_n: int = 50,
) -> None:
    """Heatmap of top significant genes across samples."""
    genes_to_plot = sig_genes[:top_n]
    if not genes_to_plot:
        return

    avail = [g for g in genes_to_plot if g in counts.index]
    if not avail:
        return

    expr = log2_normalize(counts)
    sub = expr.loc[avail, group_a_samples + group_b_samples]

    # Z-score across samples per gene
    z = sub.subtract(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, 1), axis=0)

    col_colors = (
        [plt.cm.Set1(0)] * len(group_a_samples) +
        [plt.cm.Set1(1)] * len(group_b_samples)
    )

    h = max(5, len(avail) * 0.22)
    w = max(6, len(sub.columns) * 0.35)

    fig, ax = plt.subplots(figsize=(w, h), dpi=300)
    sns.heatmap(
        z, ax=ax,
        cmap="RdBu_r", center=0,
        xticklabels=True, yticklabels=True,
        linewidths=0.0, linecolor="none",
        vmin=-2.5, vmax=2.5,
        cbar_kws={"label": "Z-score", "shrink": 0.6},
    )
    ax.set_title(f"Top {len(avail)} DE genes  ·  {prefix}", fontsize=12, pad=8)
    ax.set_xlabel("")
    ax.set_ylabel("Gene", fontsize=10)
    ax.tick_params(axis="x", labelsize=7, rotation=45)
    ax.tick_params(axis="y", labelsize=7)

    # Group label bar above heatmap
    import matplotlib.patches as mpatches
    patch_a = mpatches.Patch(color=plt.cm.Set1(0), label=group_a_label)
    patch_b = mpatches.Patch(color=plt.cm.Set1(1), label=group_b_label)
    ax.legend(handles=[patch_a, patch_b], loc="upper right",
              bbox_to_anchor=(1.18, 1.02), fontsize=8)

    plt.tight_layout()

    base = os.path.join(outdir, f"{prefix}.top_genes_heatmap")
    fig.savefig(f"{base}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{base}.pdf", bbox_inches="tight")
    plt.close(fig)


# =========================================================
# Main
# =========================================================

def main():
    parser = argparse.ArgumentParser(
        description="RNA-seq differential expression analysis.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Inputs
    parser.add_argument("--counts", required=True,
        help="Count/expression matrix (TSV or CSV). Rows=genes, Cols=samples.")
    parser.add_argument("--metadata", required=True,
        help="Sample metadata table (TSV or CSV). Must have 'sample' and 'group' columns.")
    parser.add_argument("--group-a", required=True,
        help="Label for the reference/control group (must match values in metadata 'group' column).")
    parser.add_argument("--group-b", required=True,
        help="Label for the condition/treatment group (must match values in metadata 'group' column).")

    # Method
    parser.add_argument(
        "--method",
        choices=["deseq2", "ttest", "mannwhitney"],
        default="deseq2",
        help="Statistical method. deseq2 requires pydeseq2 and integer raw counts.",
    )

    # Filtering
    parser.add_argument("--min-count", type=int, default=10,
        help="Minimum count in at least --min-samples samples to keep a gene.")
    parser.add_argument("--min-samples", type=int, default=2,
        help="Minimum number of samples a gene must pass --min-count in.")

    # Thresholds
    parser.add_argument("--fc-thresh", type=float, default=1.0,
        help="log2 fold-change threshold for significance (|log2FC| >= this).")
    parser.add_argument("--padj-thresh", type=float, default=0.05,
        help="Adjusted p-value threshold for significance.")

    # Output
    parser.add_argument("--prefix", default=None,
        help="Output file prefix. Defaults to '<group_b>_vs_<group_a>'.")
    parser.add_argument("--top-n-heatmap", type=int, default=50,
        help="Number of top significant genes to show in heatmap.")
    parser.add_argument("--outdir", required=True,
        help="Output directory.")

    args = parser.parse_args()

    init_style(
        font_family=getattr(args, 'font_family', None),
        font_size=getattr(args, 'font_size', None),
    )
    ensure_dir(args.outdir)

    prefix = args.prefix or f"{args.group_b}_vs_{args.group_a}"

    # ----------------------------------------------------------
    # Load data
    # ----------------------------------------------------------
    print(f"[INFO] Loading count matrix: {args.counts}")
    counts = load_count_matrix(args.counts)
    print(f"[INFO] Matrix dimensions: {counts.shape[0]} genes × {counts.shape[1]} samples")

    print(f"[INFO] Loading metadata: {args.metadata}")
    meta = load_sample_metadata(args.metadata)

    # Resolve sample groups
    group_a_samples = meta[meta["group"] == args.group_a]["sample"].tolist()
    group_b_samples = meta[meta["group"] == args.group_b]["sample"].tolist()

    if not group_a_samples:
        raise ValueError(f"No samples found for group '{args.group_a}' in metadata.")
    if not group_b_samples:
        raise ValueError(f"No samples found for group '{args.group_b}' in metadata.")

    # Verify samples exist in count matrix
    missing_a = [s for s in group_a_samples if s not in counts.columns]
    missing_b = [s for s in group_b_samples if s not in counts.columns]
    if missing_a or missing_b:
        raise ValueError(
            f"Samples in metadata not found in count matrix.\n"
            f"  Missing from group A: {missing_a}\n"
            f"  Missing from group B: {missing_b}"
        )

    print(f"[INFO] Group A ({args.group_a}): {len(group_a_samples)} samples")
    print(f"[INFO] Group B ({args.group_b}): {len(group_b_samples)} samples")

    # ----------------------------------------------------------
    # Filter low-expression genes
    # ----------------------------------------------------------
    all_samples = group_a_samples + group_b_samples
    counts_sub = counts[all_samples]
    counts_filtered = filter_low_expression(counts_sub, args.min_count, args.min_samples)
    print(f"[INFO] Genes after low-expression filter: {len(counts_filtered)} / {len(counts_sub)}")

    if len(counts_filtered) == 0:
        raise ValueError("All genes were filtered out. Lower --min-count or --min-samples.")

    # ----------------------------------------------------------
    # Differential expression
    # ----------------------------------------------------------
    print(f"[INFO] Running DE analysis: {args.method}  ({args.group_b} vs {args.group_a})")

    if args.method == "deseq2":
        result = run_pydeseq2(
            counts_filtered, group_a_samples, group_b_samples,
            args.group_a, args.group_b,
        )
    elif args.method == "ttest":
        result = run_ttest(counts_filtered, group_a_samples, group_b_samples)
    else:
        result = run_mannwhitney(counts_filtered, group_a_samples, group_b_samples)

    if result.empty:
        raise ValueError("DE analysis returned no results. Check input data and groups.")

    print(f"[INFO] Total genes tested: {len(result)}")

    # ----------------------------------------------------------
    # Significant genes
    # ----------------------------------------------------------
    sig = result[
        (result["padj"].notna()) &
        (result["padj"] <= args.padj_thresh) &
        (result["log2FC"].abs() >= args.fc_thresh)
    ].copy()

    n_up   = int((sig["log2FC"] >= args.fc_thresh).sum())
    n_down = int((sig["log2FC"] <= -args.fc_thresh).sum())
    print(f"[INFO] Significant: {len(sig)} genes  ({n_up} up, {n_down} down)")

    # ----------------------------------------------------------
    # Save tables
    # ----------------------------------------------------------
    full_path = os.path.join(args.outdir, f"{prefix}.de_results.tsv")
    result.to_csv(full_path, sep="\t", index=False)
    print(f"[INFO] Full DE results: {full_path}")

    sig_path = os.path.join(args.outdir, f"{prefix}.significant_genes.tsv")
    sig.to_csv(sig_path, sep="\t", index=False)
    print(f"[INFO] Significant genes: {sig_path}")

    # ----------------------------------------------------------
    # Plots
    # ----------------------------------------------------------
    print(f"[INFO] Generating volcano plot")
    plot_volcano(result, args.outdir, prefix, args.fc_thresh, args.padj_thresh)

    print(f"[INFO] Generating MA plot")
    plot_ma(result, args.outdir, prefix, args.fc_thresh, args.padj_thresh)

    print(f"[INFO] Generating heatmap")
    sig_genes_ordered = sig.sort_values("padj")["gene"].tolist()
    plot_heatmap(
        counts_filtered, sig_genes_ordered,
        group_a_samples, group_b_samples,
        args.group_a, args.group_b,
        args.outdir, prefix,
        top_n=args.top_n_heatmap,
    )

    # ----------------------------------------------------------
    # Summary
    # ----------------------------------------------------------
    summary = {
        "prefix":          prefix,
        "group_a":         args.group_a,
        "group_b":         args.group_b,
        "n_samples_a":     len(group_a_samples),
        "n_samples_b":     len(group_b_samples),
        "method":          args.method,
        "genes_tested":    len(result),
        "fc_threshold":    args.fc_thresh,
        "padj_threshold":  args.padj_thresh,
        "n_significant":   len(sig),
        "n_upregulated":   n_up,
        "n_downregulated": n_down,
    }

    pd.DataFrame([summary]).to_csv(
        os.path.join(args.outdir, "de_summary.tsv"), sep="\t", index=False
    )

    with open(os.path.join(args.outdir, "summary.txt"), "w", encoding="utf-8") as f:
        for k, v in summary.items():
            f.write(f"{k}: {v}\n")

    print(f"[DONE] Results written to: {args.outdir}")


if __name__ == "__main__":
    main()
