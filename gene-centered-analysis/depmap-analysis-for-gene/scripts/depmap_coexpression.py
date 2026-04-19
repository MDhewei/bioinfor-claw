#!/usr/bin/env python3
"""DepMap co-expression analysis for a single gene.

Computes Pearson/Spearman correlations between the query gene's expression
vector and all other genes across DepMap cell lines. Produces:
  - Ranked correlation table with FDR correction (TSV)
  - Top co-expressed genes horizontal bar plot (PNG + PDF)
  - Co-expression network visualization (PNG + PDF)

Data: real DepMap expression matrix (OmicsExpressionProteinCodingGenesTPMLogp1.csv
or similar). NOT patient/tumor samples — these are cancer cell lines.
For patient-level co-expression, use coexpression-for-gene (TCGA/GTEx).
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from typing import List, Optional

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), *(['..'] * 3), '_shared'))
    from plot_style import init_style, save_fig
except ImportError:
    def init_style(**kw): pass
    def save_fig(fig, path, close=True, **kw):
        fig.savefig(path, dpi=300, bbox_inches='tight')
        if close: import matplotlib.pyplot as _p; _p.close(fig)


# ─── scipy-free fallback implementations ─────────────────────────────────────
# These allow the script to run even when scipy is not installed.
# Uses Numerical Recipes continued-fraction method for the regularised
# incomplete beta function, which gives the two-tailed p-value for Pearson r.

def _betacf(a, b, x):
    """Continued fraction for incomplete beta (Numerical Recipes)."""
    MAXIT, EPS, FPMIN = 200, 3e-14, 1e-30
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < FPMIN: d = FPMIN
    d = 1.0 / d
    h = d
    for m in range(1, MAXIT + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN: d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN: c = FPMIN
        d = 1.0 / d; h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN: d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN: c = FPMIN
        d = 1.0 / d; delta = d * c; h *= delta
        if abs(delta - 1.0) < EPS: break
    return h


def _betainc(a, b, x):
    """Regularised incomplete beta function I_x(a, b)."""
    if x <= 0: return 0.0
    if x >= 1: return 1.0
    ln_pre = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
              + a * math.log(x) + b * math.log(1 - x))
    front = math.exp(ln_pre)
    if x < (a + 1) / (a + b + 2):
        return front * _betacf(a, b, x) / a
    else:
        return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def _pearsonr_np(x, y):
    """Pearson correlation with two-tailed p-value (numpy only)."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    n = len(x)
    if n < 3:
        return np.nan, 1.0
    xm = x - x.mean()
    ym = y - y.mean()
    denom = math.sqrt(np.dot(xm, xm) * np.dot(ym, ym))
    if denom < 1e-16:
        return 0.0, 1.0
    r = float(max(-1.0, min(1.0, np.dot(xm, ym) / denom)))
    if abs(r) == 1.0:
        return r, 0.0
    t2 = r * r * (n - 2) / (1.0 - r * r)
    df = n - 2
    p = _betainc(df / 2.0, 0.5, df / (df + t2))
    return r, p


def _rankdata_np(a):
    """Rank data array (average ties), numpy only."""
    a = np.asarray(a, dtype=np.float64)
    sorter = np.argsort(a, kind='mergesort')
    inv = np.empty_like(sorter)
    inv[sorter] = np.arange(len(a))
    a_sorted = a[sorter]
    obs = np.concatenate(([True], a_sorted[1:] != a_sorted[:-1]))
    dense = np.cumsum(obs)[inv]
    count = np.bincount(dense)
    cumcount = np.concatenate(([0], np.cumsum(count)))
    ranks = np.empty(len(a), dtype=np.float64)
    for i in range(len(a)):
        d = dense[i]
        ranks[i] = 0.5 * (cumcount[d] + cumcount[d] + count[d] + 1)
    return ranks


def _spearmanr_np(x, y):
    """Spearman rank correlation with two-tailed p-value (numpy only)."""
    return _pearsonr_np(_rankdata_np(x), _rankdata_np(y))


# Try to import scipy; fall back to numpy implementations
try:
    from scipy.stats import pearsonr as _scipy_pearsonr, spearmanr as _scipy_spearmanr
    from scipy.stats import rankdata as _scipy_rankdata
    _pearsonr = _scipy_pearsonr
    _rankdata = _scipy_rankdata
    def _spearmanr(x, y):
        res = _scipy_spearmanr(x, y)
        return res.correlation if hasattr(res, 'correlation') else res[0], \
               res.pvalue if hasattr(res, 'pvalue') else res[1]
    print("[INFO] Using scipy for correlation functions")
except ImportError:
    _pearsonr = _pearsonr_np
    _spearmanr = _spearmanr_np
    _rankdata = _rankdata_np
    print("[INFO] scipy not available — using numpy fallback for correlations")


# ─── Utilities ────────────────────────────────────────────────────────────────

def normalize_gene_symbol(gene: str) -> str:
    return re.sub(r"\s+", "", gene.strip()).upper()


def find_gene_column(columns: List[str], gene_symbol: str) -> Optional[str]:
    if gene_symbol in columns:
        return gene_symbol
    pat = re.compile(rf"^{re.escape(gene_symbol)}(\s|\(|\[|$)", re.IGNORECASE)
    matches = [c for c in columns if pat.search(str(c))]
    return matches[0] if matches else None


def extract_symbol(col: str) -> str:
    m = re.match(r"^([A-Za-z0-9_.-]+)", str(col))
    return m.group(1).upper() if m else str(col).upper()


# ─── Core analysis ────────────────────────────────────────────────────────────

def compute_correlations(
    matrix: pd.DataFrame,
    gene: str,
    method: str = "pearson",
) -> pd.DataFrame:
    """Correlate `gene` against all other genes. Returns full ranked DataFrame
    with gene_symbol, correlation, p-value, FDR."""

    col = find_gene_column(list(matrix.columns), gene)
    if col is None:
        raise ValueError(f"Gene {gene} not found in expression matrix columns.")

    target = matrix[col].dropna()
    other_cols = [c for c in matrix.columns if c != col]
    common_idx = target.index

    results = []
    corr_fn = _pearsonr if method == "pearson" else _spearmanr

    for c in other_cols:
        vec = matrix[c].reindex(common_idx).dropna()
        shared = target.index.intersection(vec.index)
        if len(shared) < 10:
            continue
        t = target.loc[shared].values
        v = vec.loc[shared].values
        if np.std(v) < 1e-6:
            continue
        try:
            r, p = corr_fn(t, v)
            if np.isnan(r):
                continue
            results.append({
                "gene_column": c,
                "gene_symbol": extract_symbol(c),
                "correlation": r,
                "pvalue": p,
                "n_samples": len(shared),
            })
        except Exception:
            continue

    df = pd.DataFrame(results)
    if len(df) == 0:
        return df

    # FDR (Benjamini-Hochberg with monotonicity enforcement)
    pvals = df["pvalue"].values
    n = len(pvals)
    sort_idx = np.argsort(pvals)
    sorted_pvals = pvals[sort_idx]
    # BH adjusted p-values: p_adj[i] = p[i] * n / (i+1)
    fdr_sorted = sorted_pvals * n / np.arange(1, n + 1)
    # Enforce monotonicity: walk backwards, each value must be <= the next
    for i in range(n - 2, -1, -1):
        fdr_sorted[i] = min(fdr_sorted[i], fdr_sorted[i + 1])
    fdr_sorted = np.minimum(fdr_sorted, 1.0)
    # Map back to original order
    fdr = np.empty(n)
    fdr[sort_idx] = fdr_sorted
    df["fdr"] = fdr
    df = df.sort_values("correlation", key=abs, ascending=False).reset_index(drop=True)
    return df


# ─── Helper: select top positive AND top negative hits ───────────────────────

def _select_top_pos_neg(df: pd.DataFrame, top_n: int) -> pd.DataFrame:
    """Select top_n/2 most positive and top_n/2 most negative correlations.
    If one side has fewer genes, the other side fills the gap."""
    pos = df[df["correlation"] > 0].sort_values("correlation", ascending=False)
    neg = df[df["correlation"] < 0].sort_values("correlation", ascending=True)
    half = top_n // 2
    # If one side is short, give the surplus to the other
    n_pos = min(len(pos), half)
    n_neg = min(len(neg), half)
    if n_pos < half:
        n_neg = min(len(neg), top_n - n_pos)
    elif n_neg < half:
        n_pos = min(len(pos), top_n - n_neg)
    selected = pd.concat([pos.head(n_pos), neg.head(n_neg)], ignore_index=True)
    return selected


# ─── Plotting: horizontal bar chart ──────────────────────────────────────────

def plot_top_coexpressed(
    df: pd.DataFrame,
    gene: str,
    top_n: int = 30,
    output_prefix: str = "coexpression_barplot",
):
    if len(df) < 1:
        print("[warn] No results to plot for bar chart")
        return
    top = _select_top_pos_neg(df, top_n).sort_values("correlation")

    fig_h = max(5, len(top) * 0.28 + 1.5)
    fig, ax = plt.subplots(figsize=(9, fig_h))

    colors = ["#d62728" if r > 0 else "#1f77b4" for r in top["correlation"]]
    ax.barh(range(len(top)), top["correlation"], color=colors,
            edgecolor="black", linewidth=0.4)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top["gene_symbol"], fontsize=9)
    ax.set_xlabel("Correlation coefficient")
    ax.set_title(f"Top {len(top)} genes co-expressed with {gene} (DepMap)")
    ax.axvline(0, color="black", linewidth=0.8)
    ax.grid(axis="x", alpha=0.3, linestyle="--")

    save_fig(fig, f"{output_prefix}.png")
    save_fig(fig, f"{output_prefix}.pdf", close=False)
    plt.close(fig)


# ─── Plotting: network ───────────────────────────────────────────────────────

def plot_network(
    df: pd.DataFrame,
    gene: str,
    top_n: int = 30,
    output_prefix: str = "coexpression_network",
):
    if len(df) < 2:
        print("[warn] Too few results for network plot")
        return
    top = _select_top_pos_neg(df, top_n).copy()

    np.random.seed(42)
    n = len(top)
    positions = {g: np.random.randn(2) for g in top["gene_symbol"]}

    # Simple spring layout
    gene_list = list(top["gene_symbol"])
    for _ in range(15):
        forces = {g: np.zeros(2) for g in gene_list}
        for i, g1 in enumerate(gene_list):
            for j, g2 in enumerate(gene_list):
                if i >= j:
                    continue
                d = positions[g2] - positions[g1]
                dist = np.linalg.norm(d) + 0.1
                forces[g1] -= d / (dist ** 2)
                forces[g2] += d / (dist ** 2)
        for g in positions:
            positions[g] += 0.01 * forces[g]

    all_pos = np.array([positions[g] for g in gene_list])
    span = all_pos.max(axis=0) - all_pos.min(axis=0) + 0.1
    all_pos = (all_pos - all_pos.min(axis=0)) / span
    for i, g in enumerate(gene_list):
        positions[g] = all_pos[i]

    # Edges: connect genes with similar correlation to query
    corr_vals = dict(zip(top["gene_symbol"], top["correlation"]))
    edges = []
    for i, g1 in enumerate(gene_list):
        for j, g2 in enumerate(gene_list):
            if i >= j:
                continue
            diff = abs(corr_vals[g1] - corr_vals[g2])
            if diff < 0.3:
                edges.append((g1, g2, 1 - diff))

    fig, ax = plt.subplots(figsize=(11, 9))

    # Draw edges
    for g1, g2, w in edges:
        xs = [positions[g1][0], positions[g2][0]]
        ys = [positions[g1][1], positions[g2][1]]
        ax.plot(xs, ys, "gray", alpha=0.25, linewidth=w * 1.5, zorder=1)

    # Draw query gene at center
    center = np.array([0.5, 0.5])
    ax.scatter(*center, s=600, c="#FFD700", edgecolors="black", linewidth=1.5, zorder=4)
    ax.text(*center, gene, ha="center", va="center", fontsize=9, fontweight="bold", zorder=5)

    # Draw co-expressed gene nodes
    for g in gene_list:
        r = corr_vals[g]
        pos = positions[g]
        color = "#d62728" if r > 0 else "#1f77b4"
        size = abs(r) * 400 + 50
        alpha = min(0.95, abs(r) + 0.4)
        ax.scatter(pos[0], pos[1], s=size, c=color, alpha=alpha,
                   edgecolors="black", linewidth=0.8, zorder=2)
        ax.text(pos[0], pos[1] + 0.025, g, ha="center", va="bottom",
                fontsize=7, fontweight="bold", zorder=3)

    # Edges from query to each gene
    for g in gene_list:
        r = corr_vals[g]
        pos = positions[g]
        lw = abs(r) * 2.5
        color = "#d62728" if r > 0 else "#1f77b4"
        ax.plot([center[0], pos[0]], [center[1], pos[1]],
                color=color, alpha=0.35, linewidth=lw, zorder=0)

    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, 1.1)
    ax.set_title(f"Co-expression network: {gene} (DepMap cell lines, top {n})")
    ax.axis("off")

    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(facecolor="#d62728", label="Positive correlation"),
        Patch(facecolor="#1f77b4", label="Negative correlation"),
        Patch(facecolor="#FFD700", label=f"Query: {gene}"),
    ], loc="upper right", fontsize=9, framealpha=0.9)

    save_fig(fig, f"{output_prefix}.png")
    save_fig(fig, f"{output_prefix}.pdf", close=False)
    plt.close(fig)


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="DepMap co-expression analysis for a single gene. "
                    "Correlates expression across cancer cell lines."
    )
    parser.add_argument("--gene", required=True, help="Gene symbol (e.g. TP53)")
    parser.add_argument("--expression-file", required=True,
                        help="DepMap expression CSV (genes as columns, cell lines as rows)")
    parser.add_argument("--metadata-file", default=None,
                        help="(optional) DepMap Model.csv for cell line annotations")
    parser.add_argument("--method", choices=["pearson", "spearman"],
                        default="pearson", help="Correlation method (default: pearson)")
    parser.add_argument("--top-n", type=int, default=30,
                        help="Number of top co-expressed genes to plot (default: 30)")
    parser.add_argument("--fdr-cutoff", type=float, default=0.01,
                        help="FDR threshold for significance (default: 0.01)")
    parser.add_argument("--min-corr", type=float, default=0.2,
                        help="Minimum |correlation| for significant genes (default: 0.2)")
    parser.add_argument("--network-top-n", type=int, default=30,
                        help="Top N genes for network visualization (default: 30)")
    parser.add_argument("--outdir", default=".", help="Output directory")
    parser.add_argument("--font-family", default=None, help="Font family override")
    parser.add_argument("--font-size", type=float, default=None, help="Base font size")
    args = parser.parse_args()

    init_style(font_family=args.font_family, font_size=args.font_size)

    os.makedirs(args.outdir, exist_ok=True)
    gene = normalize_gene_symbol(args.gene)

    # Load expression matrix
    print(f"[INFO] Loading expression matrix: {args.expression_file}")
    expr = pd.read_csv(args.expression_file, low_memory=False, index_col=0)
    # Drop non-numeric columns (metadata like SequencingID, ModelID, etc.)
    num_cols = expr.select_dtypes(include=[np.number]).columns
    dropped = len(expr.columns) - len(num_cols)
    if dropped > 0:
        print(f"[INFO] Dropped {dropped} non-numeric metadata columns")
        expr = expr[num_cols]
    print(f"[INFO] Matrix: {expr.shape[0]} cell lines × {expr.shape[1]} genes")

    # Compute correlations
    print(f"[INFO] Computing {args.method} correlations for {gene}…")
    corr_df = compute_correlations(expr, gene, method=args.method)

    if len(corr_df) == 0:
        print(f"[ERROR] No valid correlations computed for {gene}.", file=sys.stderr)
        sys.exit(1)

    # Full table
    full_path = os.path.join(args.outdir, f"{gene}.depmap_coexpression_full.tsv")
    corr_df.to_csv(full_path, sep="\t", index=False)
    print(f"Saved: {full_path} ({len(corr_df)} genes)")

    # Significant genes: |correlation| >= min_corr AND FDR <= fdr_cutoff
    sig = corr_df[
        (corr_df["fdr"] <= args.fdr_cutoff) &
        (corr_df["correlation"].abs() >= args.min_corr)
    ]
    sig_path = os.path.join(args.outdir, f"{gene}.depmap_coexpression_sig.tsv")
    sig.to_csv(sig_path, sep="\t", index=False)
    print(f"Saved: {sig_path} ({len(sig)} significant genes, |r| ≥ {args.min_corr} & FDR ≤ {args.fdr_cutoff})")

    # Bar plot
    plot_top_coexpressed(
        corr_df, gene, top_n=args.top_n,
        output_prefix=os.path.join(args.outdir, f"{gene}.depmap_coexpression_barplot"),
    )

    # Network
    plot_network(
        corr_df, gene, top_n=args.network_top_n,
        output_prefix=os.path.join(args.outdir, f"{gene}.depmap_coexpression_network"),
    )

    # Summary
    print(f"\n=== DepMap Co-expression Summary for {gene} ===")
    print(f"Total genes correlated : {len(corr_df)}")
    print(f"Significant (|r|≥{args.min_corr}, FDR≤{args.fdr_cutoff}): {len(sig)}")
    if len(corr_df) > 0:
        top = corr_df.iloc[0]
        print(f"Top co-expressed       : {top['gene_symbol']} (r={top['correlation']:.3f}, FDR={top['fdr']:.2e})")
    print(f"Method                 : {args.method}")
    print(f"Output dir             : {args.outdir}")


if __name__ == "__main__":
    main()
