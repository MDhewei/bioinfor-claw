#!/usr/bin/env python3
"""DepMap co-essentiality analysis for a single gene.

Computes Pearson/Spearman correlations between the query gene's CRISPR
dependency score and all other genes' dependency scores across DepMap cell
lines. Genes that are co-essential tend to work in the same pathway or
complex — they are "needed together" in the same cell lines.

Produces:
  - Ranked correlation table with FDR correction (TSV)
  - Top co-essential genes horizontal bar plot (PNG + PDF)
  - Co-essentiality network visualization (PNG + PDF)

Data: real DepMap CRISPR gene effect CSV (CRISPRGeneEffect.csv).
NOT patient/tumor samples — these are cancer cell lines.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from typing import List, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', '_shared'))
from plot_style import init_style, save_fig


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
    """Correlate `gene` dependency scores against all other genes.
    Returns full ranked DataFrame with gene_symbol, correlation, p-value, FDR."""
    from scipy.stats import pearsonr, spearmanr

    col = find_gene_column(list(matrix.columns), gene)
    if col is None:
        raise ValueError(f"Gene {gene} not found in essentiality matrix columns.")

    target = matrix[col].dropna()
    other_cols = [c for c in matrix.columns if c != col]

    results = []
    corr_fn = pearsonr if method == "pearson" else spearmanr

    for c in other_cols:
        vec = matrix[c].reindex(target.index).dropna()
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

    # FDR (Benjamini-Hochberg)
    from scipy.stats import rankdata
    pvals = df["pvalue"].values
    n = len(pvals)
    ranked = rankdata(pvals)
    fdr = np.minimum((pvals * n) / ranked, 1.0)
    df["fdr"] = fdr
    df = df.sort_values("correlation", key=abs, ascending=False).reset_index(drop=True)
    return df


# ─── Plotting: horizontal bar chart ──────────────────────────────────────────

def plot_top_coessential(
    df: pd.DataFrame,
    gene: str,
    top_n: int = 30,
    output_prefix: str = "coessentiality_barplot",
):
    if len(df) < 1:
        print("[warn] No results to plot for bar chart")
        return
    top = df.head(top_n).copy().sort_values("correlation")

    fig_h = max(5, len(top) * 0.28 + 1.5)
    fig, ax = plt.subplots(figsize=(9, fig_h))

    colors = ["#d62728" if r > 0 else "#1f77b4" for r in top["correlation"]]
    ax.barh(range(len(top)), top["correlation"], color=colors,
            edgecolor="black", linewidth=0.4)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top["gene_symbol"], fontsize=9)
    ax.set_xlabel("Correlation coefficient (dependency scores)")
    ax.set_title(f"Top {len(top)} genes co-essential with {gene} (DepMap)")
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
    output_prefix: str = "coessentiality_network",
):
    if len(df) < 2:
        print("[warn] Too few results for network plot")
        return
    top = df.head(top_n).copy()

    np.random.seed(42)
    n = len(top)
    gene_list = list(top["gene_symbol"])
    positions = {g: np.random.randn(2) for g in gene_list}

    # Simple spring layout
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

    for g1, g2, w in edges:
        xs = [positions[g1][0], positions[g2][0]]
        ys = [positions[g1][1], positions[g2][1]]
        ax.plot(xs, ys, "gray", alpha=0.25, linewidth=w * 1.5, zorder=1)

    # Query gene at center
    center = np.array([0.5, 0.5])
    ax.scatter(*center, s=600, c="#FFD700", edgecolors="black", linewidth=1.5, zorder=4)
    ax.text(*center, gene, ha="center", va="center", fontsize=9, fontweight="bold", zorder=5)

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

    for g in gene_list:
        r = corr_vals[g]
        pos = positions[g]
        lw = abs(r) * 2.5
        color = "#d62728" if r > 0 else "#1f77b4"
        ax.plot([center[0], pos[0]], [center[1], pos[1]],
                color=color, alpha=0.35, linewidth=lw, zorder=0)

    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, 1.1)
    ax.set_title(f"Co-essentiality network: {gene} (DepMap CRISPR, top {n})")
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
        description="DepMap co-essentiality analysis for a single gene. "
                    "Correlates CRISPR dependency scores across cell lines."
    )
    parser.add_argument("--gene", required=True, help="Gene symbol (e.g. TP53)")
    parser.add_argument("--essentiality-file", required=True,
                        help="DepMap CRISPR gene effect CSV (genes as columns, cell lines as rows)")
    parser.add_argument("--metadata-file", default=None,
                        help="(optional) DepMap Model.csv for cell line annotations")
    parser.add_argument("--method", choices=["pearson", "spearman"],
                        default="pearson", help="Correlation method (default: pearson)")
    parser.add_argument("--top-n", type=int, default=30,
                        help="Number of top co-essential genes to plot (default: 30)")
    parser.add_argument("--fdr-cutoff", type=float, default=0.05,
                        help="FDR threshold for reporting (default: 0.05)")
    parser.add_argument("--network-top-n", type=int, default=30,
                        help="Top N genes for network visualization (default: 30)")
    parser.add_argument("--outdir", default=".", help="Output directory")
    parser.add_argument("--font-family", default=None, help="Font family override")
    parser.add_argument("--font-size", type=float, default=None, help="Base font size")
    args = parser.parse_args()

    init_style(font_family=args.font_family, font_size=args.font_size)

    os.makedirs(args.outdir, exist_ok=True)
    gene = normalize_gene_symbol(args.gene)

    # Load essentiality matrix
    print(f"[INFO] Loading essentiality matrix: {args.essentiality_file}")
    dep = pd.read_csv(args.essentiality_file, low_memory=False, index_col=0)
    print(f"[INFO] Matrix: {dep.shape[0]} cell lines × {dep.shape[1]} genes")

    # Compute correlations
    print(f"[INFO] Computing {args.method} correlations for {gene}…")
    corr_df = compute_correlations(dep, gene, method=args.method)

    if len(corr_df) == 0:
        print(f"[ERROR] No valid correlations computed for {gene}.", file=sys.stderr)
        sys.exit(1)

    # Full table
    full_path = os.path.join(args.outdir, f"{gene}.depmap_coessentiality_full.tsv")
    corr_df.to_csv(full_path, sep="\t", index=False)
    print(f"Saved: {full_path} ({len(corr_df)} genes)")

    # FDR-filtered
    sig = corr_df[corr_df["fdr"] <= args.fdr_cutoff]
    sig_path = os.path.join(args.outdir, f"{gene}.depmap_coessentiality_sig.tsv")
    sig.to_csv(sig_path, sep="\t", index=False)
    print(f"Saved: {sig_path} ({len(sig)} significant genes, FDR ≤ {args.fdr_cutoff})")

    # Bar plot
    plot_top_coessential(
        corr_df, gene, top_n=args.top_n,
        output_prefix=os.path.join(args.outdir, f"{gene}.depmap_coessentiality_barplot"),
    )

    # Network
    plot_network(
        corr_df, gene, top_n=args.network_top_n,
        output_prefix=os.path.join(args.outdir, f"{gene}.depmap_coessentiality_network"),
    )

    # Summary
    print(f"\n=== DepMap Co-essentiality Summary for {gene} ===")
    print(f"Total genes correlated : {len(corr_df)}")
    print(f"Significant (FDR≤{args.fdr_cutoff})  : {len(sig)}")
    if len(corr_df) > 0:
        top = corr_df.iloc[0]
        print(f"Top co-essential       : {top['gene_symbol']} (r={top['correlation']:.3f}, FDR={top['fdr']:.2e})")
    print(f"Method                 : {args.method}")
    print(f"Output dir             : {args.outdir}")


if __name__ == "__main__":
    main()
