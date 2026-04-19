#!/usr/bin/env python3
"""
Drug Sensitivity for Gene: Correlate gene expression/dependency with drug sensitivity.

Identifies drugs whose sensitivity correlates with a gene's expression or essentiality
across cell lines using DepMap and PRISM data.
"""

import argparse
import json
import sys
import os
from typing import Dict, Tuple, Optional, List
import urllib.request
import urllib.error

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys as _sys, os as _os
try:
    _sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), *(['..'] * 3), '_shared'))
    from plot_style import init_style
except ImportError:
    def init_style(**kw): pass  # graceful fallback if _shared not available
from scipy.stats import pearsonr, spearmanr


def download_depmap_expression() -> pd.DataFrame:
    """
    Download DepMap expression data (log TPM).
    Returns DataFrame: rows=genes, columns=cell lines, values=log expression
    """
    print("Downloading DepMap expression data...")
    url = "https://figshare.com/ndownloader/files/34990036"

    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            # This is a gzip file, need to handle
            import gzip
            data = gzip.GzipFile(fileobj=response).read().decode()
            # Parse as TSV
            from io import StringIO
            df = pd.read_csv(StringIO(data), sep="\t", index_col=0)
            print(f"Loaded expression matrix: {df.shape[0]} genes x {df.shape[1]} cell lines")
            return df
    except Exception as e:
        print(f"Warning: Failed to download expression data: {e}", file=sys.stderr)
        return pd.DataFrame()


def download_prism_sensitivity() -> pd.DataFrame:
    """
    Download PRISM drug sensitivity data.
    Returns DataFrame: rows=cell lines, columns=drugs, values=log fold change (sensitivity)
    """
    print("Downloading PRISM drug sensitivity data...")
    url = "https://figshare.com/ndownloader/files/20237710"

    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            import gzip
            data = gzip.GzipFile(fileobj=response).read().decode()
            from io import StringIO
            df = pd.read_csv(StringIO(data), sep="\t", index_col=0)
            print(f"Loaded PRISM data: {df.shape[0]} cell lines x {df.shape[1]} drugs")
            return df
    except Exception as e:
        print(f"Warning: Failed to download PRISM data: {e}", file=sys.stderr)
        return pd.DataFrame()


def load_local_expression(filepath: str) -> pd.DataFrame:
    """Load expression matrix from local file."""
    try:
        df = pd.read_csv(filepath, sep="\t", index_col=0)
        print(f"Loaded expression matrix: {df.shape[0]} genes x {df.shape[1]} cell lines")
        return df
    except Exception as e:
        print(f"Error loading expression file: {e}", file=sys.stderr)
        return pd.DataFrame()


def load_local_prism(filepath: str) -> pd.DataFrame:
    """Load PRISM data from local file."""
    try:
        df = pd.read_csv(filepath, sep="\t", index_col=0)
        print(f"Loaded PRISM data: {df.shape[0]} cell lines x {df.shape[1]} drugs")
        return df
    except Exception as e:
        print(f"Error loading PRISM file: {e}", file=sys.stderr)
        return pd.DataFrame()


def extract_gene_vector(
    expression_df: pd.DataFrame,
    gene_symbol: str,
    omics_type: str = "expression"
) -> Optional[pd.Series]:
    """
    Extract gene expression or dependency vector.
    """
    # Check for exact match first
    if gene_symbol in expression_df.index:
        return expression_df.loc[gene_symbol]

    # Try case-insensitive match
    matches = [idx for idx in expression_df.index if idx.upper() == gene_symbol.upper()]
    if matches:
        return expression_df.loc[matches[0]]

    # Try partial match
    matches = [idx for idx in expression_df.index if gene_symbol in idx]
    if matches:
        print(f"Warning: Using {matches[0]} instead of {gene_symbol}")
        return expression_df.loc[matches[0]]

    return None


def compute_correlations(
    gene_vector: pd.Series,
    prism_df: pd.DataFrame,
    method: str = "spearman",
    min_cell_lines: int = 50
) -> pd.DataFrame:
    """
    Compute correlation between gene vector and each drug sensitivity profile.
    """
    results = []

    for drug_col in prism_df.columns:
        drug_data = prism_df[drug_col]

        # Find common cell lines with data for both gene and drug
        common_indices = gene_vector.index.intersection(drug_data.index)
        common_indices = common_indices[~(gene_vector[common_indices].isna() | drug_data[common_indices].isna())]

        if len(common_indices) < min_cell_lines:
            continue

        gene_vals = gene_vector[common_indices].values
        drug_vals = drug_data[common_indices].values

        # Compute correlation
        if method == "pearson":
            r, pval = pearsonr(gene_vals, drug_vals)
        else:  # spearman
            r, pval = spearmanr(gene_vals, drug_vals)

        results.append({
            "drug": drug_col,
            "r": r,
            "pvalue": pval,
            "n_cell_lines": len(common_indices)
        })

    results_df = pd.DataFrame(results)

    # Apply FDR correction (Benjamini-Hochberg)
    if len(results_df) > 0:
        from scipy.stats import rankdata
        pvals = results_df["pvalue"].values
        n = len(pvals)
        ranked = rankdata(pvals)
        fdr = (pvals * n) / ranked
        fdr = np.minimum(fdr, 1.0)
        results_df["fdr"] = fdr

    return results_df


def plot_scatter_grid(
    gene_vector: pd.Series,
    prism_df: pd.DataFrame,
    results_df: pd.DataFrame,
    top_n: int = 4,
    output_file: str = "scatter_plot_top4.png"
):
    """
    Create 2x2 grid of scatter plots for top drugs.
    """
    if len(results_df) < 1:
        print("Warning: No results to plot")
        return

    top_results = results_df.nlargest(top_n, lambda x: abs(x["r"]))

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    for idx, (ax, (_, row)) in enumerate(zip(axes, top_results.iterrows())):
        drug = row["drug"]
        if drug not in prism_df.columns:
            continue

        drug_data = prism_df[drug]
        common_idx = gene_vector.index.intersection(drug_data.index)
        common_idx = common_idx[~(gene_vector[common_idx].isna() | drug_data[common_idx].isna())]

        gene_vals = gene_vector[common_idx].values
        drug_vals = drug_data[common_idx].values

        # Scatter plot
        ax.scatter(gene_vals, drug_vals, alpha=0.6, s=30, edgecolors="black", linewidth=0.5)

        # Regression line
        z = np.polyfit(gene_vals, drug_vals, 1)
        p = np.poly1d(z)
        x_line = np.linspace(gene_vals.min(), gene_vals.max(), 100)
        ax.plot(x_line, p(x_line), "r-", linewidth=2, alpha=0.8)

        ax.set_xlabel("Gene Expression (log TPM)", fontsize=10)
        ax.set_ylabel("Drug Sensitivity", fontsize=10)
        ax.set_title(f"{drug}\nR={row['r']:.3f}, p={row['pvalue']:.2e}", fontsize=11, fontweight="bold")
        ax.grid(alpha=0.3, linestyle="--")

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_file}")


def plot_waterfall(
    results_df: pd.DataFrame,
    top_n: int = 20,
    output_file: str = "waterfall_plot.png"
):
    """
    Create waterfall/bar plot of top drug correlations.
    """
    if len(results_df) < 1:
        print("Warning: No results to plot")
        return

    top_results = results_df.nlargest(top_n, lambda x: abs(x["r"])).copy()
    top_results = top_results.sort_values("r")

    fig, ax = plt.subplots(figsize=(12, 8))

    colors = ["#d62728" if r > 0 else "#1f77b4" for r in top_results["r"]]
    bars = ax.barh(range(len(top_results)), top_results["r"], color=colors, edgecolor="black", linewidth=1)

    ax.set_yticks(range(len(top_results)))
    ax.set_yticklabels(top_results["drug"], fontsize=9)
    ax.set_xlabel("Pearson/Spearman Correlation", fontsize=12, fontweight="bold")
    ax.set_title(f"Top {top_n} Drug-Gene Correlations", fontsize=14, fontweight="bold")
    ax.axvline(0, color="black", linewidth=1)
    ax.grid(axis="x", alpha=0.3, linestyle="--")

    # Add value labels
    for i, (idx, row) in enumerate(top_results.iterrows()):
        r = row["r"]
        ax.text(r, i, f"  {r:.3f}", va="center", fontsize=8)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Correlate gene expression/dependency with drug sensitivity."
    )
    parser.add_argument("--gene", required=True, help="Gene symbol")
    parser.add_argument("--omics-type", choices=["expression", "dependency", "copy_number"],
                       default="expression", help="Type of gene measurement")
    parser.add_argument("--correlation-method", choices=["pearson", "spearman"],
                       default="spearman", help="Correlation method")
    parser.add_argument("--fdr-cutoff", type=float, default=0.05, help="FDR threshold")
    parser.add_argument("--top-n", type=int, default=20, help="Top N drugs to visualize")
    parser.add_argument("--min-cell-lines", type=int, default=50,
                       help="Minimum cell lines with data")
    parser.add_argument("--depmap-dir", default=None, help="Local DepMap data directory")
    parser.add_argument("--prism-file", default=None, help="Local PRISM file")
    parser.add_argument("--expression-file", default=None, help="Local expression file")
    parser.add_argument("--outdir", default=".", help="Output directory")

    args = parser.parse_args()

    init_style(
        font_family=getattr(args, 'font_family', None),
        font_size=getattr(args, 'font_size', None),
    )

    os.makedirs(args.outdir, exist_ok=True)

    # Load data
    if args.expression_file:
        expression_df = load_local_expression(args.expression_file)
    else:
        expression_df = download_depmap_expression()

    if expression_df.empty:
        print("Error: Could not load expression data", file=sys.stderr)
        return 1

    if args.prism_file:
        prism_df = load_local_prism(args.prism_file)
    else:
        prism_df = download_prism_sensitivity()

    if prism_df.empty:
        print("Error: Could not load PRISM data", file=sys.stderr)
        return 1

    # Extract gene vector
    print(f"Extracting {args.gene} vector...")
    gene_vector = extract_gene_vector(expression_df, args.gene, args.omics_type)

    if gene_vector is None:
        print(f"Error: Gene {args.gene} not found in expression matrix", file=sys.stderr)
        return 1

    print(f"Gene vector: {len(gene_vector)} cell lines")

    # Compute correlations
    print("Computing correlations...")
    results_df = compute_correlations(
        gene_vector, prism_df,
        method=args.correlation_method,
        min_cell_lines=args.min_cell_lines
    )

    if results_df.empty:
        print("Error: No correlations computed", file=sys.stderr)
        return 1

    print(f"Computed {len(results_df)} drug correlations")

    # Filter by FDR
    fdr_results = results_df[results_df["fdr"] <= args.fdr_cutoff].copy()
    print(f"Significant correlations (FDR <= {args.fdr_cutoff}): {len(fdr_results)}")

    # Sort by absolute correlation
    fdr_results = fdr_results.sort_values("r", key=abs, ascending=False)

    # Save results
    results_file = os.path.join(args.outdir, "drug_correlation.tsv")
    fdr_results.to_csv(results_file, sep="\t", index=False)
    print(f"Saved: {results_file}")

    # Generate plots
    if len(fdr_results) > 0:
        scatter_file = os.path.join(args.outdir, "scatter_plot_top4.png")
        plot_scatter_grid(gene_vector, prism_df, fdr_results, output_file=scatter_file)

        waterfall_file = os.path.join(args.outdir, "waterfall_plot.png")
        plot_waterfall(fdr_results, top_n=args.top_n, output_file=waterfall_file)

    # Summary file
    summary_file = os.path.join(args.outdir, "correlation_summary.txt")
    with open(summary_file, "w") as f:
        f.write(f"Drug Sensitivity Analysis for {args.gene}\n")
        f.write(f"=" * 60 + "\n\n")
        f.write(f"Gene: {args.gene}\n")
        f.write(f"Omics type: {args.omics_type}\n")
        f.write(f"Correlation method: {args.correlation_method}\n")
        f.write(f"FDR cutoff: {args.fdr_cutoff}\n")
        f.write(f"Minimum cell lines: {args.min_cell_lines}\n\n")
        f.write(f"Total drugs tested: {len(results_df)}\n")
        f.write(f"Significant correlations (FDR <= {args.fdr_cutoff}): {len(fdr_results)}\n")
        if len(fdr_results) > 0:
            f.write(f"Strongest positive correlation: {fdr_results.iloc[0]['drug']} (R={fdr_results.iloc[0]['r']:.3f})\n")
            f.write(f"Strongest negative correlation: {fdr_results.iloc[-1]['drug']} (R={fdr_results.iloc[-1]['r']:.3f})\n")
        f.write("\n\nTop 20 Correlations:\n")
        f.write(fdr_results.head(20).to_string(index=False))

    print(f"Saved: {summary_file}")
    print(f"\nAnalysis complete! Results saved to {args.outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
