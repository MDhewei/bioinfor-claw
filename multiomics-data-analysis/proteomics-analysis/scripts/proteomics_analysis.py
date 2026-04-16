#!/usr/bin/env python3
"""
Proteomics Data Analysis

Performs comprehensive analysis of quantitative proteomics data:
  - QC metrics (sample correlation, missing values, CV)
  - Normalization (median, quantile, log2-median, VSN-approx)
  - Imputation (minprob, k-NN, zero)
  - Batch correction (mean-centering)
  - Differential expression (Welch t-test with FDR)
  - Visualization (volcano plots, heatmaps, correlation matrices)

Input: Protein intensity matrix (protein × sample) + metadata
Output: Normalized/imputed matrix, DE results, QC plots, volcano plots

Dependencies: pandas, numpy, scipy, matplotlib
"""

import argparse
import os
import sys
import warnings
from typing import Dict, List, Optional, Tuple

import matplotlib
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), '..', '..', '..', '_shared'))
from plot_style import init_style
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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


def load_protein_matrix(path: str) -> pd.DataFrame:
    """Load protein intensity matrix (protein × sample)."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Protein matrix not found: {path}")

    sep = "\t" if path.endswith(".tsv") or path.endswith(".txt") else ","
    df = pd.read_csv(path, sep=sep, index_col=0)
    df.index.name = "protein"
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(how="all")

    return df


def load_metadata(path: str, group_col: str = "group") -> pd.DataFrame:
    """Load sample metadata."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Metadata file not found: {path}")

    sep = "\t" if path.endswith(".tsv") or path.endswith(".txt") else ","
    meta = pd.read_csv(path, sep=sep)

    required = {"sample_id", group_col.lower()}
    cols_lower = {col.lower(): col for col in meta.columns}
    missing = required - set(cols_lower.keys())

    if missing:
        raise ValueError(
            f"Metadata must contain columns: sample_id, {group_col.lower()}. "
            f"Missing: {missing}. Found: {list(meta.columns)}"
        )

    # Rename columns to lowercase
    meta.columns = [col.lower() for col in meta.columns]
    return meta


def detect_log_scale(data: pd.DataFrame) -> bool:
    """Auto-detect if data is already log-transformed."""
    max_val = data.values.max()
    return max_val <= 20 or (data.values < 1).sum() > len(data) * 0.1


def apply_log2_transform(data: pd.DataFrame) -> pd.DataFrame:
    """Apply log2 transformation if needed."""
    if detect_log_scale(data):
        return data.copy()
    return np.log2(data.replace(0, np.nan) + 1)


def filter_proteins(data: pd.DataFrame, min_valid_fraction: float = 0.7) -> pd.DataFrame:
    """Filter proteins with insufficient valid measurements."""
    valid_count = (~data.isna()).sum(axis=1)
    threshold = min_valid_fraction * len(data.columns)
    mask = valid_count >= threshold
    return data[mask].copy()


def normalize_median(data: pd.DataFrame) -> pd.DataFrame:
    """Median normalization."""
    medians = data.median(axis=0)
    return data.subtract(medians, axis=1)


def normalize_quantile(data: pd.DataFrame) -> pd.DataFrame:
    """Quantile normalization (rank-based)."""
    rank_data = data.rank(axis=0, method="average")
    sorted_data = np.sort(data.values, axis=0)
    mean_sorted = np.mean(sorted_data, axis=1)

    normalized = np.empty_like(data.values, dtype=float)
    for i in range(data.shape[1]):
        idx = np.argsort(rank_data.iloc[:, i].values).argsort()
        normalized[:, i] = mean_sorted[idx]

    return pd.DataFrame(normalized, index=data.index, columns=data.columns)


def normalize_log2_median(data: pd.DataFrame) -> pd.DataFrame:
    """Log2-transform then median normalize."""
    log_data = apply_log2_transform(data)
    return normalize_median(log_data)


def normalize_vsn_approx(data: pd.DataFrame) -> pd.DataFrame:
    """Variance-stabilizing normalization (asinh approximation)."""
    # Simple VSN approximation using asinh transform
    return np.arcsinh(data / 2) * 2


def normalize_none(data: pd.DataFrame) -> pd.DataFrame:
    """No normalization."""
    return data.copy()


def impute_minprob(data: pd.DataFrame) -> pd.DataFrame:
    """Impute missing values with low random values (minprob strategy)."""
    imputed = data.copy()
    for col in imputed.columns:
        valid = imputed[col].dropna()
        if len(valid) > 0:
            mean_val = valid.mean()
            std_val = valid.std()
            if np.isnan(std_val) or std_val == 0:
                std_val = 1.0
            # Fill missing with mean - 1.8*SD + small noise
            fill_value = mean_val - 1.8 * std_val
            mask = imputed[col].isna()
            imputed.loc[mask, col] = fill_value + np.random.normal(0, 0.1, mask.sum())

    return imputed


def impute_knn_approx(data: pd.DataFrame, k: int = 5) -> pd.DataFrame:
    """Impute missing values using k-NN approximation."""
    imputed = data.copy()

    # Compute correlation between proteins based on non-missing values
    for idx, row in imputed.iterrows():
        missing_cols = row.isna()
        if not missing_cols.any():
            continue

        # Find k nearest proteins
        distances = []
        for other_idx, other_row in imputed.iterrows():
            if idx == other_idx:
                continue
            # Correlation on shared non-missing values
            shared = ~(row.isna() | other_row.isna())
            if shared.sum() > 0:
                corr = row[shared].corr(other_row[shared])
                if not np.isnan(corr):
                    distances.append((1 - corr, other_idx))

        distances.sort()
        neighbors = [x[1] for x in distances[:k]]

        # Impute using neighbor means
        for col in imputed.columns[missing_cols]:
            neighbor_vals = [imputed.loc[n, col] for n in neighbors if not pd.isna(imputed.loc[n, col])]
            if neighbor_vals:
                imputed.loc[idx, col] = np.mean(neighbor_vals)

    return imputed


def impute_zero(data: pd.DataFrame) -> pd.DataFrame:
    """Impute missing values with zero."""
    return data.fillna(0)


def impute_none(data: pd.DataFrame) -> pd.DataFrame:
    """No imputation."""
    return data.copy()


def apply_batch_correction(data: pd.DataFrame, meta: pd.DataFrame, batch_col: str) -> pd.DataFrame:
    """Apply simple batch correction (mean-centering per batch)."""
    corrected = data.copy()

    if batch_col not in meta.columns:
        print(f"Warning: batch column '{batch_col}' not found in metadata. Skipping batch correction.")
        return corrected

    for batch in meta[batch_col].unique():
        batch_samples = meta[meta[batch_col] == batch]["sample_id"].tolist()
        batch_samples = [s for s in batch_samples if s in corrected.columns]

        if len(batch_samples) > 0:
            batch_mean = corrected[batch_samples].mean(axis=1)
            corrected[batch_samples] = corrected[batch_samples].subtract(batch_mean, axis=0)

    return corrected


def compute_cv(data: pd.DataFrame) -> pd.Series:
    """Compute coefficient of variation per sample."""
    cv = (data.std(axis=0) / data.mean(axis=0).abs()) * 100
    return cv


def compute_qc_metrics(data: pd.DataFrame, meta: pd.DataFrame) -> Dict:
    """Compute QC metrics per sample."""
    metrics = {}
    for col in data.columns:
        valid_data = data[col].dropna()
        if len(valid_data) > 0:
            metrics[col] = {
                "n_proteins": len(valid_data),
                "mean_intensity": float(valid_data.mean()),
                "cv": float(compute_cv(data[[col]])[col]) if valid_data.mean() != 0 else 0,
            }
    return metrics


def plot_sample_correlation(data: pd.DataFrame, outdir: str) -> None:
    """Plot sample correlation heatmap."""
    corr = data.corr(method="pearson")

    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
    im = ax.imshow(corr.values, cmap="RdYlBu_r", vmin=0, vmax=1)
    ax.set_xticks(range(len(corr.columns)))
    ax.set_yticks(range(len(corr.columns)))
    ax.set_xticklabels(corr.columns, rotation=45, ha="right")
    ax.set_yticklabels(corr.columns)
    ax.set_title("Sample Correlation Matrix")

    for i in range(len(corr.columns)):
        for j in range(len(corr.columns)):
            text = ax.text(j, i, f"{corr.iloc[i, j]:.2f}",
                         ha="center", va="center", color="black", fontsize=8)

    plt.colorbar(im, ax=ax, label="Pearson Correlation")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "sample_correlation_heatmap.png"), dpi=300, bbox_inches="tight")
    plt.close()


def plot_missing_values(raw_data: pd.DataFrame, outdir: str) -> None:
    """Plot missing value heatmap."""
    missing = raw_data.isna().astype(int)

    fig, ax = plt.subplots(figsize=(12, len(raw_data) / 20 + 2), dpi=300)
    im = ax.imshow(missing.values, cmap="Reds", aspect="auto")
    ax.set_xticks(range(len(missing.columns)))
    ax.set_xticklabels(missing.columns, rotation=45, ha="right")
    ax.set_ylabel("Proteins")
    ax.set_title("Missing Value Heatmap (Red = Missing)")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "missing_value_heatmap.png"), dpi=300, bbox_inches="tight")
    plt.close()


def plot_cv_distribution(data: pd.DataFrame, outdir: str) -> None:
    """Plot CV distribution per sample."""
    cv = compute_cv(data)

    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    ax.bar(range(len(cv)), cv.values, color="#4ECDC4", edgecolor="black", alpha=0.7)
    ax.set_xticks(range(len(cv)))
    ax.set_xticklabels(cv.index, rotation=45, ha="right")
    ax.set_ylabel("Coefficient of Variation (%)")
    ax.set_title("CV Distribution per Sample")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "cv_distribution.png"), dpi=300, bbox_inches="tight")
    plt.close()


def welch_ttest(group1: np.ndarray, group2: np.ndarray) -> Tuple[float, float]:
    """Compute Welch's t-test (unequal variance)."""
    g1_valid = group1[~np.isnan(group1)]
    g2_valid = group2[~np.isnan(group2)]

    if len(g1_valid) < 2 or len(g2_valid) < 2:
        return np.nan, 1.0

    t_stat, p_val = stats.ttest_ind(g1_valid, g2_valid, equal_var=False, nan_policy="omit")
    return float(t_stat), float(p_val)


def benjamini_hochberg_correction(p_values: np.ndarray) -> np.ndarray:
    """Apply Benjamini-Hochberg FDR correction."""
    n = len(p_values)
    sorted_idx = np.argsort(p_values)
    sorted_p = p_values[sorted_idx]

    adjusted = np.ones_like(sorted_p)
    for i, p in enumerate(sorted_p):
        adjusted[i] = p * n / (i + 1)

    for i in range(len(adjusted) - 2, -1, -1):
        adjusted[i] = min(adjusted[i], adjusted[i + 1])

    fdr = np.empty_like(adjusted)
    fdr[sorted_idx] = adjusted
    fdr = np.clip(fdr, 0, 1)

    return fdr


def identify_de_proteins(
    data: pd.DataFrame,
    meta: pd.DataFrame,
    ref_group: str,
    group_col: str = "group",
    fdr_cutoff: float = 0.05,
    fc_cutoff: float = 1.0,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Identify differentially expressed proteins."""
    ref_samples = meta[meta[group_col] == ref_group]["sample_id"].tolist()
    ref_samples = [s for s in ref_samples if s in data.columns]

    if len(ref_samples) < 2:
        raise ValueError(f"Reference group '{ref_group}' has < 2 samples")

    treat_samples = meta[meta[group_col] != ref_group]["sample_id"].tolist()
    treat_samples = [s for s in treat_samples if s in data.columns]

    if len(treat_samples) < 2:
        raise ValueError("Treatment group has < 2 samples")

    ref_data = data[ref_samples]
    treat_data = data[treat_samples]

    # Compute statistics
    results = []
    for protein_id in data.index:
        ref_vals = ref_data.loc[protein_id].values
        treat_vals = treat_data.loc[protein_id].values

        ref_mean = np.nanmean(ref_vals)
        treat_mean = np.nanmean(treat_vals)
        log2fc = treat_mean - ref_mean

        t_stat, p_val = welch_ttest(ref_vals, treat_vals)

        results.append({
            "protein": protein_id,
            "log2FC": log2fc,
            "mean_ref": ref_mean,
            "mean_treat": treat_mean,
            "pvalue": p_val,
        })

    de_df = pd.DataFrame(results)

    # Apply FDR correction
    de_df["fdr"] = benjamini_hochberg_correction(de_df["pvalue"].values)
    de_df["significant"] = (de_df["fdr"] <= fdr_cutoff) & (np.abs(de_df["log2FC"]) >= fc_cutoff)

    significant = de_df[de_df["significant"]].copy()

    return de_df, significant


def plot_de_volcano(de_df: pd.DataFrame, outdir: str) -> None:
    """Plot volcano plot for DE proteins."""
    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)

    log10_fdr = -np.log10(de_df["fdr"].replace(0, 1e-300))
    log2fc = de_df["log2FC"]

    colors = np.where(
        de_df["significant"],
        "#FF6B6B",  # Significant: red
        "#CCCCCC"   # Non-significant: gray
    )

    ax.scatter(log2fc, log10_fdr, c=colors, alpha=0.6, s=20, edgecolor="none")
    ax.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=1, label="FDR = 0.05")
    ax.axvline(1.0, color="blue", linestyle="--", linewidth=1, label="log2FC = 1.0")
    ax.axvline(-1.0, color="blue", linestyle="--", linewidth=1)

    ax.set_xlabel("log2(Fold Change)")
    ax.set_ylabel("-log10(FDR)")
    ax.set_title("Volcano Plot: Differential Protein Expression")
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "de_volcano_plot.png"), dpi=300, bbox_inches="tight")
    plt.close()


def plot_de_heatmap(data: pd.DataFrame, de_df: pd.DataFrame, meta: pd.DataFrame, outdir: str, top_n: int = 50) -> None:
    """Plot heatmap of top DE proteins."""
    if len(de_df) == 0:
        return

    top_proteins = de_df.nlargest(top_n, "pvalue")["protein"].tolist()
    top_proteins = [p for p in top_proteins if p in data.index][:top_n]

    if len(top_proteins) == 0:
        return

    heatmap_data = data.loc[top_proteins].copy()

    # Log transform and Z-score
    heatmap_data = np.log2(heatmap_data.replace(0, np.nan) + 1)
    z_scores = (heatmap_data - heatmap_data.mean(axis=1, keepdims=True)) / heatmap_data.std(axis=1, keepdims=True)

    fig, ax = plt.subplots(figsize=(12, len(top_proteins) / 2 + 2), dpi=300)
    im = ax.imshow(z_scores.fillna(0).values, cmap="RdBu_r", aspect="auto", vmin=-3, vmax=3)
    ax.set_xticks(range(len(z_scores.columns)))
    ax.set_xticklabels(z_scores.columns, rotation=45, ha="right")
    ax.set_yticks(range(0, len(z_scores.index), max(1, len(z_scores.index) // 20)))
    ax.set_yticklabels(z_scores.index[ax.get_yticks().astype(int).tolist()], fontsize=8)
    ax.set_title(f"Top {top_n} DE Proteins (Z-scored)")
    plt.colorbar(im, ax=ax, label="Z-score")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "de_heatmap_top50.png"), dpi=300, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Proteomics data analysis (QC, normalization, DE)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", required=True, help="Protein intensity matrix (TSV)")
    parser.add_argument("--metadata", required=True, help="Sample metadata (TSV)")
    parser.add_argument("--mode", choices=["qc", "normalize", "differential", "all"], default="all",
                       help="Analysis mode")
    parser.add_argument("--quant-type", choices=["tmt", "lfq", "dia"], default="lfq",
                       help="Quantification platform")
    parser.add_argument("--normalization", choices=["median", "quantile", "vsn_approx", "log2_median", "none"],
                       default="median", help="Normalization method")
    parser.add_argument("--imputation", choices=["minprob", "knn_approx", "zero", "none"], default="minprob",
                       help="Imputation method")
    parser.add_argument("--group-col", default="group", help="Group column in metadata")
    parser.add_argument("--ref-group", default="", help="Reference group for fold-change direction")
    parser.add_argument("--fdr-cutoff", type=float, default=0.05, help="FDR threshold")
    parser.add_argument("--fc-cutoff", type=float, default=1.0, help="log2 FC threshold")
    parser.add_argument("--min-valid-values", type=float, default=0.7, help="Min fraction of valid values per protein")
    parser.add_argument("--batch-col", default="", help="Batch column in metadata")
    parser.add_argument("--outdir", default="./proteomics_output", help="Output directory")

    args = parser.parse_args()

    init_style(
        font_family=getattr(args, 'font_family', None),
        font_size=getattr(args, 'font_size', None),
    )

    ensure_dir(args.outdir)

    # Load data
    print(f"Loading protein matrix from {args.input}...")
    raw_data = load_protein_matrix(args.input)
    print(f"  Loaded {len(raw_data)} proteins, {len(raw_data.columns)} samples")

    print(f"Loading metadata from {args.metadata}...")
    meta = load_metadata(args.metadata, args.group_col)
    print(f"  Loaded {len(meta)} samples")

    # QC analysis (on raw data)
    if args.mode in ["qc", "all"]:
        print("\nRunning QC analysis...")
        metrics = compute_qc_metrics(raw_data, meta)

        qc_file = os.path.join(args.outdir, "qc_summary.txt")
        with open(qc_file, "w") as f:
            f.write("QC Summary\n")
            f.write("=" * 50 + "\n")
            for sample_id, m in metrics.items():
                f.write(f"\nSample: {sample_id}\n")
                for k, v in m.items():
                    f.write(f"  {k}: {v}\n")

        plot_sample_correlation(raw_data.fillna(raw_data.mean()), args.outdir)
        plot_missing_values(raw_data, args.outdir)
        plot_cv_distribution(raw_data.fillna(raw_data.mean()), args.outdir)
        print(f"  QC plots saved")

    # Normalize and impute
    data = raw_data.copy()
    data = filter_proteins(data, args.min_valid_values)
    print(f"After filtering: {len(data)} proteins")

    if args.mode in ["normalize", "differential", "all"]:
        print(f"\nApplying normalization ({args.normalization})...")
        norm_func = {
            "median": normalize_median,
            "quantile": normalize_quantile,
            "log2_median": normalize_log2_median,
            "vsn_approx": normalize_vsn_approx,
            "none": normalize_none,
        }.get(args.normalization, normalize_median)

        data = norm_func(data)

        print(f"Applying imputation ({args.imputation})...")
        imp_func = {
            "minprob": impute_minprob,
            "knn_approx": impute_knn_approx,
            "zero": impute_zero,
            "none": impute_none,
        }.get(args.imputation, impute_minprob)

        data = imp_func(data)

        if args.batch_col:
            print(f"Applying batch correction (batch_col={args.batch_col})...")
            data = apply_batch_correction(data, meta, args.batch_col)

        if args.mode in ["normalize", "all"]:
            norm_file = os.path.join(args.outdir, "normalized_proteins.tsv")
            data.to_csv(norm_file, sep="\t")
            print(f"  Normalized data saved to {norm_file}")

    # Differential expression
    if args.mode in ["differential", "all"]:
        if not args.ref_group:
            print("Error: --ref-group required for differential expression analysis")
            sys.exit(1)

        print(f"\nRunning differential expression analysis (ref_group={args.ref_group})...")
        try:
            de_df, significant_de = identify_de_proteins(
                data, meta, args.ref_group, args.group_col,
                args.fdr_cutoff, args.fc_cutoff
            )

            de_file = os.path.join(args.outdir, "de_results.tsv")
            de_df.to_csv(de_file, sep="\t", index=False)
            print(f"  DE results saved to {de_file}")
            print(f"  Significant proteins: {len(significant_de)}")

            sig_file = os.path.join(args.outdir, "significant_proteins.tsv")
            significant_de.to_csv(sig_file, sep="\t", index=False)

            plot_de_volcano(de_df, args.outdir)
            plot_de_heatmap(raw_data, de_df, meta, args.outdir)
        except Exception as e:
            print(f"  Error in DE analysis: {e}")

    # Summary
    summary_file = os.path.join(args.outdir, "analysis_summary.txt")
    with open(summary_file, "w") as f:
        f.write("Proteomics Analysis Summary\n")
        f.write("=" * 60 + "\n")
        f.write(f"Input matrix: {args.input}\n")
        f.write(f"Proteins analyzed: {len(data)}\n")
        f.write(f"Samples: {len(data.columns)}\n")
        f.write(f"Quantification type: {args.quant_type}\n")
        f.write(f"Normalization: {args.normalization}\n")
        f.write(f"Imputation: {args.imputation}\n")
        f.write(f"Mode: {args.mode}\n")
        if args.ref_group:
            f.write(f"Reference group: {args.ref_group}\n")
            f.write(f"FDR cutoff: {args.fdr_cutoff}\n")
            f.write(f"log2FC cutoff: {args.fc_cutoff}\n")

    print(f"\nAnalysis complete. Results saved to {args.outdir}")


if __name__ == "__main__":
    main()
