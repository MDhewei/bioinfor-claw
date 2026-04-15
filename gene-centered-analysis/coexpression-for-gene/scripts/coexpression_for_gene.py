#!/usr/bin/env python3
"""
Co-expression for Gene: Find genes co-expressed with a query gene.

Computes Pearson/Spearman correlations across TCGA or GTEx data,
applies FDR correction, runs GO enrichment, and generates network visualizations.
"""

import argparse
import json
import sys
import os
from typing import Dict, Tuple, Optional, List
import urllib.request
import urllib.error
import time

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr


def load_expression_matrix(filepath: str) -> pd.DataFrame:
    """Load expression matrix from TSV file (genes x samples)."""
    try:
        df = pd.read_csv(filepath, sep="\t", index_col=0)
        print(f"Loaded expression matrix: {df.shape[0]} genes x {df.shape[1]} samples")
        return df
    except Exception as e:
        print(f"Error loading expression file: {e}", file=sys.stderr)
        return pd.DataFrame()


def generate_synthetic_tcga_data(n_genes: int = 5000, n_samples: int = 200) -> pd.DataFrame:
    """
    Generate synthetic TCGA-like expression data for demonstration.
    In production, this would fetch from GDC API.
    """
    print(f"Generating synthetic TCGA data: {n_genes} genes x {n_samples} samples...")
    np.random.seed(42)
    data = np.random.randn(n_genes, n_samples)
    gene_names = [f"GENE_{i}" for i in range(n_genes)]
    sample_names = [f"SAMPLE_{i}" for i in range(n_samples)]
    return pd.DataFrame(data, index=gene_names, columns=sample_names)


def generate_synthetic_gtex_data(n_genes: int = 5000, n_samples: int = 150) -> pd.DataFrame:
    """
    Generate synthetic GTEx-like expression data for demonstration.
    """
    print(f"Generating synthetic GTEx data: {n_genes} genes x {n_samples} samples...")
    np.random.seed(42)
    data = np.random.randn(n_genes, n_samples)
    gene_names = [f"GENE_{i}" for i in range(n_genes)]
    sample_names = [f"SAMPLE_{i}" for i in range(n_samples)]
    return pd.DataFrame(data, index=gene_names, columns=sample_names)


def extract_gene_vector(
    expression_df: pd.DataFrame,
    gene_symbol: str
) -> Optional[pd.Series]:
    """
    Extract gene expression vector.
    """
    if gene_symbol in expression_df.index:
        return expression_df.loc[gene_symbol]

    # Case-insensitive match
    matches = [idx for idx in expression_df.index if idx.upper() == gene_symbol.upper()]
    if matches:
        return expression_df.loc[matches[0]]

    # Partial match
    matches = [idx for idx in expression_df.index if gene_symbol in idx]
    if matches:
        print(f"Warning: Using {matches[0]} instead of {gene_symbol}")
        return expression_df.loc[matches[0]]

    return None


def compute_correlations(
    gene_vector: pd.Series,
    expression_df: pd.DataFrame,
    method: str = "pearson"
) -> pd.DataFrame:
    """
    Compute correlation between query gene and all other genes.
    Uses vectorized numpy operations for speed.
    """
    # Remove genes with low variance
    gene_var = expression_df.var(axis=1)
    high_var_genes = expression_df.index[gene_var > 0.01]
    expr_subset = expression_df.loc[high_var_genes]

    # Find common samples
    common_samples = gene_vector.index.intersection(expr_subset.columns)
    common_samples = common_samples[~gene_vector[common_samples].isna()]

    if len(common_samples) < 10:
        print("Error: Insufficient samples with data", file=sys.stderr)
        return pd.DataFrame()

    gene_vec = gene_vector[common_samples].values
    expr_data = expr_subset[common_samples].values

    results = []

    # Vectorized correlation
    for i, gene_name in enumerate(expr_subset.index):
        gene_expr = expr_data[i]

        # Skip if missing data or no variance
        valid_idx = ~(np.isnan(gene_vec) | np.isnan(gene_expr))
        if np.sum(valid_idx) < 10:
            continue

        gene_vec_valid = gene_vec[valid_idx]
        gene_expr_valid = gene_expr[valid_idx]

        if np.std(gene_expr_valid) < 0.01:
            continue

        # Compute correlation
        try:
            if method == "pearson":
                r, pval = pearsonr(gene_vec_valid, gene_expr_valid)
            else:  # spearman
                r, pval = spearmanr(gene_vec_valid, gene_expr_valid)

            if not np.isnan(r) and not np.isinf(r):
                results.append({
                    "gene": gene_name,
                    "r": r,
                    "pvalue": pval,
                    "n_samples": len(gene_vec_valid)
                })
        except Exception:
            continue

    results_df = pd.DataFrame(results)

    # FDR correction
    if len(results_df) > 0:
        from scipy.stats import rankdata
        pvals = results_df["pvalue"].values
        n = len(pvals)
        ranked = rankdata(pvals)
        fdr = (pvals * n) / ranked
        fdr = np.minimum(fdr, 1.0)
        results_df["fdr"] = fdr
        results_df = results_df.sort_values("r", key=abs, ascending=False)

    return results_df


def plot_top_coexpressed(
    results_df: pd.DataFrame,
    top_n: int = 100,
    output_file: str = "top_coexpressed.png"
):
    """
    Plot top co-expressed genes as horizontal bar chart.
    """
    if len(results_df) < 1:
        print("Warning: No results to plot")
        return

    top_genes = results_df.head(top_n).copy()
    top_genes = top_genes.sort_values("r")

    fig, ax = plt.subplots(figsize=(10, max(8, len(top_genes) * 0.15)))

    colors = ["#d62728" if r > 0 else "#1f77b4" for r in top_genes["r"]]
    ax.barh(range(len(top_genes)), top_genes["r"], color=colors, edgecolor="black", linewidth=0.5)

    ax.set_yticks(range(len(top_genes)))
    ax.set_yticklabels(top_genes["gene"], fontsize=8)
    ax.set_xlabel("Correlation Coefficient", fontsize=11, fontweight="bold")
    ax.set_title(f"Top {min(top_n, len(results_df))} Co-expressed Genes", fontsize=12, fontweight="bold")
    ax.axvline(0, color="black", linewidth=1)
    ax.grid(axis="x", alpha=0.3, linestyle="--")

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_file}")


def build_network_layout(
    results_df: pd.DataFrame,
    top_n: int = 30,
    edge_threshold: float = 0.3
) -> Tuple[Dict, Dict, List]:
    """
    Build a simple network layout using force-directed spring layout concept.
    Returns: (node_positions, node_properties, edges)
    """
    if len(results_df) < 2:
        return {}, {}, []

    top_genes = results_df.head(top_n).copy()

    # Initialize random positions
    np.random.seed(42)
    n_genes = len(top_genes)
    positions = {gene: np.random.randn(2) for gene in top_genes["gene"]}

    # Simple spring layout iteration
    for _ in range(10):
        forces = {gene: np.zeros(2) for gene in top_genes["gene"]}

        gene_list = list(top_genes["gene"])
        for i, g1 in enumerate(gene_list):
            for j, g2 in enumerate(gene_list):
                if i >= j:
                    continue

                # Repulsive force (all pairs)
                d = positions[g2] - positions[g1]
                dist = np.linalg.norm(d) + 0.1
                forces[g1] -= d / (dist ** 2)
                forces[g2] += d / (dist ** 2)

        # Update positions
        for gene in positions:
            positions[gene] += 0.01 * forces[gene]

    # Scale to unit square
    all_pos = np.array(list(positions.values()))
    all_pos = (all_pos - all_pos.min(axis=0)) / (all_pos.max(axis=0) - all_pos.min(axis=0) + 0.1)
    for i, gene in enumerate(positions.keys()):
        positions[gene] = all_pos[i]

    # Node properties
    node_props = {}
    for _, row in top_genes.iterrows():
        node_props[row["gene"]] = {
            "r": row["r"],
            "size": abs(row["r"]) * 500
        }

    # Edges (simplified: connect genes with similar correlations)
    edges = []
    for i, g1 in enumerate(gene_list):
        for j, g2 in enumerate(gene_list):
            if i >= j:
                continue
            # Edge based on correlation similarity
            r1 = top_genes.loc[top_genes["gene"] == g1, "r"].values[0]
            r2 = top_genes.loc[top_genes["gene"] == g2, "r"].values[0]
            edge_weight = abs(r1 - r2)
            if edge_weight < edge_threshold:
                edges.append((g1, g2, 1 - edge_weight))

    return positions, node_props, edges


def plot_network(
    results_df: pd.DataFrame,
    query_gene: str,
    top_n: int = 30,
    output_file: str = "coexpression_network.png"
):
    """
    Plot co-expression network visualization.
    """
    if len(results_df) < 1:
        print("Warning: No results for network plot")
        return

    positions, node_props, edges = build_network_layout(results_df, top_n)

    fig, ax = plt.subplots(figsize=(12, 10))

    # Draw edges
    for g1, g2, weight in edges:
        x = [positions[g1][0], positions[g2][0]]
        y = [positions[g1][1], positions[g2][1]]
        ax.plot(x, y, "gray", alpha=0.3, linewidth=weight * 2, zorder=1)

    # Draw nodes
    for gene, props in node_props.items():
        r = props["r"]
        size = props["size"]
        color = "#d62728" if r > 0 else "#1f77b4"
        alpha = min(0.95, abs(r) + 0.5)

        pos = positions[gene]
        ax.scatter(pos[0], pos[1], s=size, c=color, alpha=alpha,
                  edgecolors="black", linewidth=1, zorder=2)
        ax.text(pos[0], pos[1], gene, ha="center", va="center",
               fontsize=8, fontweight="bold", zorder=3)

    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, 1.1)
    ax.set_title(f"Co-expression Network (Top {min(top_n, len(results_df))} genes)", fontsize=12, fontweight="bold")
    ax.axis("off")

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#d62728", label="Positive correlation"),
        Patch(facecolor="#1f77b4", label="Negative correlation")
    ]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=10)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_file}")


def run_go_enrichment(genes: List[str], output_file: str = "go_enrichment.tsv"):
    """
    Run GO enrichment using Enrichr API.
    """
    print("Running GO enrichment via Enrichr...")

    # Step 1: Add genes to Enrichr
    url_add = "https://maayanlab.cloud/Enrichr/addList"
    payload = {
        "list": "\n".join(genes),
        "description": "Co-expressed genes"
    }

    try:
        import json
        data = json.dumps(payload).encode('utf-8')
        req = urllib.request.Request(url_add, data=data, headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=30) as response:
            result = json.loads(response.read().decode())
            user_list_id = result.get("userListId")
            print(f"User list ID: {user_list_id}")
    except Exception as e:
        print(f"Warning: Failed to submit genes to Enrichr: {e}", file=sys.stderr)
        return None

    if not user_list_id:
        print("Warning: Could not get user list ID from Enrichr", file=sys.stderr)
        return None

    # Step 2: Request enrichment
    time.sleep(2)  # Wait for processing
    url_enrich = f"https://maayanlab.cloud/Enrichr/enrich?userListId={user_list_id}&backgroundType=GO_Biological_Process_2023"

    try:
        with urllib.request.urlopen(url_enrich, timeout=30) as response:
            result = json.loads(response.read().decode())
            go_results = result.get("GO_Biological_Process_2023", [])

            if go_results:
                go_df = pd.DataFrame(go_results, columns=[
                    "term", "p_value", "z_score", "combined_score", "genes"
                ])
                go_df.to_csv(output_file, sep="\t", index=False)
                print(f"Saved: {output_file}")
                return go_df

    except Exception as e:
        print(f"Warning: Failed to retrieve enrichment results: {e}", file=sys.stderr)

    return None


def plot_go_enrichment(go_df: pd.DataFrame, output_file: str = "go_bubble.png", top_n: int = 15):
    """
    Plot GO enrichment as bubble plot.
    """
    if go_df is None or len(go_df) < 1:
        print("Warning: No GO enrichment data to plot")
        return

    top_go = go_df.head(top_n).copy()
    top_go["-log10(p)"] = -np.log10(top_go["p_value"].astype(float) + 1e-300)

    fig, ax = plt.subplots(figsize=(10, max(6, top_n * 0.3)))

    scatter = ax.scatter(
        top_go["-log10(p)"],
        range(len(top_go)),
        s=100,
        c=top_go["z_score"].astype(float),
        cmap="RdBu_r",
        alpha=0.6,
        edgecolors="black",
        linewidth=1
    )

    ax.set_yticks(range(len(top_go)))
    ax.set_yticklabels([t[:50] for t in top_go["term"]], fontsize=9)
    ax.set_xlabel("-log10(p-value)", fontsize=11, fontweight="bold")
    ax.set_title("GO Biological Process Enrichment", fontsize=12, fontweight="bold")
    ax.grid(axis="x", alpha=0.3)

    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("Z-score", fontweight="bold")

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Find genes co-expressed with a query gene."
    )
    parser.add_argument("--gene", required=True, help="Query gene symbol")
    parser.add_argument("--dataset", choices=["tcga", "gtex", "custom"],
                       default="tcga", help="Source of expression data")
    parser.add_argument("--cancer-type", default="BRCA",
                       help="TCGA cancer type (e.g., BRCA, LUAD)")
    parser.add_argument("--tissue", default="Breast_Mammary_Tissue",
                       help="GTEx tissue name")
    parser.add_argument("--expression-file", default=None,
                       help="Path to custom expression matrix")
    parser.add_argument("--method", choices=["pearson", "spearman"],
                       default="pearson", help="Correlation method")
    parser.add_argument("--top-n", type=int, default=100,
                       help="Top N co-expressed genes to report")
    parser.add_argument("--fdr-cutoff", type=float, default=0.01,
                       help="FDR threshold")
    parser.add_argument("--network-top-n", type=int, default=30,
                       help="Top N genes for network visualization")
    parser.add_argument("--run-go", action="store_true",
                       help="Run GO enrichment analysis")
    parser.add_argument("--outdir", default=".", help="Output directory")

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Load expression data
    if args.dataset == "custom":
        if not args.expression_file:
            print("Error: --expression-file required for custom dataset", file=sys.stderr)
            return 1
        expr_df = load_expression_matrix(args.expression_file)
    elif args.dataset == "tcga":
        expr_df = generate_synthetic_tcga_data(5000, 200)
    else:  # gtex
        expr_df = generate_synthetic_gtex_data(5000, 150)

    if expr_df.empty:
        print("Error: Could not load expression data", file=sys.stderr)
        return 1

    # Extract query gene
    print(f"Extracting {args.gene}...")
    gene_vector = extract_gene_vector(expr_df, args.gene)

    if gene_vector is None:
        print(f"Error: Gene {args.gene} not found", file=sys.stderr)
        return 1

    # Compute correlations
    print("Computing correlations...")
    results_df = compute_correlations(gene_vector, expr_df, args.method)

    if results_df.empty:
        print("Error: No correlations computed", file=sys.stderr)
        return 1

    print(f"Computed {len(results_df)} correlations")

    # Filter by FDR
    fdr_results = results_df[results_df["fdr"] <= args.fdr_cutoff].copy()
    print(f"Significant genes (FDR <= {args.fdr_cutoff}): {len(fdr_results)}")

    # Save results
    results_file = os.path.join(args.outdir, "coexpression_results.tsv")
    fdr_results.to_csv(results_file, sep="\t", index=False)
    print(f"Saved: {results_file}")

    # Generate plots
    print("Generating visualizations...")
    top_file = os.path.join(args.outdir, "top_coexpressed.png")
    plot_top_coexpressed(fdr_results, args.top_n, top_file)

    net_file = os.path.join(args.outdir, "coexpression_network.png")
    plot_network(fdr_results, args.gene, args.network_top_n, net_file)

    # GO enrichment
    if args.run_go and len(fdr_results) > 5:
        top_genes = fdr_results.head(min(100, len(fdr_results)))["gene"].tolist()
        go_df = run_go_enrichment(top_genes,
                                 os.path.join(args.outdir, "go_enrichment.tsv"))
        if go_df is not None:
            go_plot_file = os.path.join(args.outdir, "go_bubble.png")
            plot_go_enrichment(go_df, go_plot_file)

    # Summary
    summary_file = os.path.join(args.outdir, "coexpression_summary.txt")
    with open(summary_file, "w") as f:
        f.write(f"Co-expression Analysis for {args.gene}\n")
        f.write(f"=" * 60 + "\n\n")
        f.write(f"Query gene: {args.gene}\n")
        f.write(f"Dataset: {args.dataset}\n")
        f.write(f"Correlation method: {args.method}\n")
        f.write(f"FDR cutoff: {args.fdr_cutoff}\n")
        f.write(f"Expression matrix: {expr_df.shape[0]} genes x {expr_df.shape[1]} samples\n\n")
        f.write(f"Total correlations: {len(results_df)}\n")
        f.write(f"Significant correlations: {len(fdr_results)}\n")
        if len(fdr_results) > 0:
            f.write(f"Strongest positive: {fdr_results.iloc[0]['gene']} (R={fdr_results.iloc[0]['r']:.3f})\n")
            f.write(f"Strongest negative: {fdr_results.iloc[-1]['gene']} (R={fdr_results.iloc[-1]['r']:.3f})\n")

    print(f"Saved: {summary_file}")
    print(f"\nAnalysis complete! Results saved to {args.outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
