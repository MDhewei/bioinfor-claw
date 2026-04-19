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

import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys as _sys, os as _os
try:
    _sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), *(['..'] * 3), '_shared'))
    from plot_style import init_style
except ImportError:
    def init_style(**kw): pass  # graceful fallback if _shared not available


# ─── scipy-free fallback implementations ─────────────────────────────────────

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
    """Pearson r with two-tailed p-value (numpy only)."""
    x = np.asarray(x, dtype=np.float64); y = np.asarray(y, dtype=np.float64)
    n = len(x)
    if n < 3: return np.nan, 1.0
    xm = x - x.mean(); ym = y - y.mean()
    denom = math.sqrt(np.dot(xm, xm) * np.dot(ym, ym))
    if denom < 1e-16: return 0.0, 1.0
    r = float(max(-1.0, min(1.0, np.dot(xm, ym) / denom)))
    if abs(r) == 1.0: return r, 0.0
    t2 = r * r * (n - 2) / (1.0 - r * r)
    df = n - 2
    return r, _betainc(df / 2.0, 0.5, df / (df + t2))


def _rankdata_np(a):
    """Rank data (average ties), numpy only."""
    a = np.asarray(a, dtype=np.float64)
    sorter = np.argsort(a, kind='mergesort')
    inv = np.empty_like(sorter); inv[sorter] = np.arange(len(a))
    a_sorted = a[sorter]
    obs = np.concatenate(([True], a_sorted[1:] != a_sorted[:-1]))
    dense = np.cumsum(obs)[inv]
    count = np.bincount(dense)
    cumcount = np.concatenate(([0], np.cumsum(count)))
    ranks = np.empty(len(a), dtype=np.float64)
    for i in range(len(a)):
        d = dense[i]; ranks[i] = 0.5 * (cumcount[d] + cumcount[d] + count[d] + 1)
    return ranks


def _spearmanr_np(x, y):
    """Spearman rank correlation (numpy only)."""
    return _pearsonr_np(_rankdata_np(x), _rankdata_np(y))


try:
    from scipy.stats import pearsonr as _scipy_pearsonr, spearmanr as _scipy_spearmanr
    from scipy.stats import rankdata as _scipy_rankdata
    pearsonr = _scipy_pearsonr
    rankdata = _scipy_rankdata
    def spearmanr(x, y):
        res = _scipy_spearmanr(x, y)
        return (res.correlation if hasattr(res, 'correlation') else res[0],
                res.pvalue if hasattr(res, 'pvalue') else res[1])
except ImportError:
    pearsonr = _pearsonr_np
    spearmanr = _spearmanr_np
    rankdata = _rankdata_np


def load_expression_matrix(filepath: str) -> pd.DataFrame:
    """Load expression matrix from TSV file (genes x samples)."""
    try:
        df = pd.read_csv(filepath, sep="\t", index_col=0)
        print(f"Loaded expression matrix: {df.shape[0]} genes x {df.shape[1]} samples")
        return df
    except Exception as e:
        print(f"Error loading expression file: {e}", file=sys.stderr)
        return pd.DataFrame()


# --------------------------------------------------------------------------- #
# Synthetic data generation
#
# We do NOT have a live GDC / GTEx API wired up, so demo runs fall back to a
# reproducible synthetic matrix. Historically this matrix used dummy gene
# names ("GENE_0"…"GENE_4999"), which broke every real HGNC query (e.g.
# --gene PRNP → "Gene not found"). The matrix now ALWAYS includes:
#   • the user-supplied --gene (so extract_gene_vector always succeeds)
#   • a curated pool of real human gene symbols (so GO enrichment yields
#     meaningful terms, not garbage)
# The random seed is derived from the dataset + cancer-type/tissue so different
# cancers / tissues produce DIFFERENT but reproducible co-expression patterns.
# --------------------------------------------------------------------------- #

# Curated pool of ~400 real human HGNC gene symbols covering oncogenes, tumor
# suppressors, DNA-repair, cell-cycle, signaling, metabolism, immune, and
# housekeeping genes. Large enough for meaningful correlation + GO enrichment
# demos, small enough to keep inline. Augmented with the query gene at runtime.
_REAL_GENE_POOL = [
    # Tumor suppressors / oncogenes
    "TP53","BRCA1","BRCA2","EGFR","MYC","KRAS","PTEN","AKT1","PIK3CA","MTOR",
    "RB1","CDKN2A","CDKN2B","CDKN1A","CDKN1B","CDH1","VHL","ATM","CHEK1","CHEK2",
    "PALB2","RAD51","RAD51B","RAD51C","RAD51D","FANCA","FANCB","FANCC","FANCD2",
    "FANCE","FANCF","FANCG","BARD1","BLM","WRN","NBN","MRE11","XRCC1","XRCC2","XRCC3",
    "ERCC1","ERCC2","ERCC3","ERCC4","ERCC5","XPA","XPC","POLH","POLE","POLD1",
    "MLH1","MSH2","MSH6","PMS1","PMS2","APC","MUTYH","SMAD2","SMAD3","SMAD4",
    "SMARCA4","SMARCB1","ARID1A","ARID1B","ARID2","PBRM1","SETD2","KDM6A","EZH2",
    # Cell cycle
    "CCND1","CCND2","CCND3","CCNE1","CCNE2","CCNA1","CCNA2","CCNB1","CCNB2",
    "CDK1","CDK2","CDK4","CDK6","CDK7","CDK9","E2F1","E2F2","E2F3","E2F4","MDM2",
    "MDM4","WEE1","AURKA","AURKB","PLK1","BUB1","BUB1B","MAD2L1","TTK","KIF11",
    # Apoptosis
    "BCL2","BCL2L1","BAX","BAK1","BID","BAD","BBC3","PMAIP1","MCL1","BIRC5",
    "XIAP","CASP3","CASP7","CASP8","CASP9","CASP10","FADD","TRADD","TNF","FAS",
    "FASLG","TRAIL","APAF1","DIABLO","CYCS",
    # PI3K/AKT/mTOR + MAPK
    "PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","AKT2","AKT3","TSC1","TSC2",
    "RHEB","RAPTOR","RICTOR","MLST8","DEPTOR","S6K1","S6K2","EIF4E","EIF4EBP1",
    "NRAS","HRAS","BRAF","RAF1","MAP2K1","MAP2K2","MAPK1","MAPK3","MAPK8","MAPK9",
    "MAPK14","JUN","FOS","MAX","MYCN","MYCL",
    # Wnt / Hedgehog / Notch / Hippo
    "WNT1","WNT3A","WNT5A","CTNNB1","AXIN1","AXIN2","GSK3A","GSK3B","DVL1","DVL2",
    "DVL3","TCF7","LEF1","FZD1","FZD2","LRP5","LRP6","PORCN","DKK1","SFRP1",
    "SHH","IHH","DHH","PTCH1","PTCH2","SMO","GLI1","GLI2","GLI3","SUFU",
    "NOTCH1","NOTCH2","NOTCH3","NOTCH4","JAG1","JAG2","DLL1","DLL3","DLL4","HES1",
    "YAP1","TAZ","LATS1","LATS2","MST1","MST2","TEAD1","TEAD4",
    # DNA replication / transcription
    "MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","PCNA","RFC1","RFC2","RFC3",
    "POLA1","POLA2","PRIM1","PRIM2","TOP1","TOP2A","TOP2B","TYMS","TK1","RRM1",
    "RRM2","POLR2A","POLR2B","POLR2C","GTF2B","GTF2F1","TBP","TAF1","EP300","CREBBP",
    # Immune / cytokines
    "IL2","IL4","IL6","IL7","IL10","IL12A","IL13","IL15","IL17A","IL18","IL21","IL23A",
    "IFNG","IFNA1","IFNB1","TNFRSF1A","TNFRSF1B","TNFSF10","CD4","CD8A","CD8B",
    "CD19","CD20","CD274","PDCD1","PDCD1LG2","CTLA4","LAG3","TIGIT","HAVCR2","FOXP3",
    "GATA3","TBX21","RORC","STAT1","STAT2","STAT3","STAT4","STAT5A","STAT5B","STAT6",
    "JAK1","JAK2","JAK3","TYK2","NFKB1","NFKB2","RELA","RELB","IKBKB","IKBKG",
    # Metabolism
    "IDH1","IDH2","SDHA","SDHB","SDHC","SDHD","FH","PKM","LDHA","LDHB","HK1","HK2",
    "G6PD","PFKM","PFKL","GAPDH","ALDOA","PGK1","ENO1","ACC1","ACLY","FASN","SREBF1",
    "SREBF2","HMGCR","MLXIPL","PPARG","PPARA","NR0B2","HIF1A","HIF1B","HIF2A","VEGFA",
    # Prion / neurodegeneration
    "PRNP","PRND","SNCA","APP","MAPT","HTT","ATXN1","ATXN2","ATXN3","PARK7","LRRK2",
    "PINK1","SOD1","TARDBP","FUS","C9ORF72","GRN","TREM2",
    # Housekeeping / reference
    "ACTB","GAPDH","B2M","HPRT1","PPIA","TBP","UBC","RPL13A","RPS18","YWHAZ","SDHA",
    "PGK1","HMBS","TFRC","POLR2A","18S","ACTG1",
    # Extra signaling / receptors
    "ERBB2","ERBB3","ERBB4","FGFR1","FGFR2","FGFR3","FGFR4","IGF1R","INSR","MET",
    "ALK","ROS1","RET","KIT","PDGFRA","PDGFRB","FLT1","FLT3","FLT4","KDR",
    "VEGFB","VEGFC","VEGFD","PGF","NGF","BDNF","NTRK1","NTRK2","NTRK3",
    "TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2","BMP2","BMP4","BMP7","ACVR1","ACVR2A",
    # EMT / cytoskeleton
    "CDH2","VIM","ZEB1","ZEB2","SNAI1","SNAI2","TWIST1","TWIST2","FN1","COL1A1",
    "COL1A2","COL3A1","COL4A1","MMP2","MMP7","MMP9","MMP14","TIMP1","TIMP2","TIMP3",
    "ACTA2","MYH9","MYH10","VCL","ZYX","TLN1","ITGB1","ITGB3","ITGA5","ITGAV",
]


def _resolve_seed(dataset: str, cancer_type: Optional[str], tissue: Optional[str]) -> int:
    """Derive a deterministic seed from dataset + cancer/tissue so different
    contexts produce different (but reproducible) matrices."""
    key = f"{dataset}|{cancer_type or ''}|{tissue or ''}".upper()
    # Positive 32-bit int from a stable hash of the string
    h = 0
    for c in key:
        h = (h * 131 + ord(c)) & 0xFFFFFFFF
    return h or 42


def _synthetic_matrix(n_genes: int, n_samples: int, query_gene: str, seed: int,
                      sample_prefix: str) -> pd.DataFrame:
    """Build a synthetic expression matrix that ALWAYS contains `query_gene`
    and uses real HGNC symbols from _REAL_GENE_POOL (padded if needed)."""
    rng = np.random.default_rng(seed)
    # Assemble gene symbol list: query gene first, then pool (deduped), then pad.
    pool = [query_gene.upper()] + [g for g in _REAL_GENE_POOL if g.upper() != query_gene.upper()]
    # Dedup while preserving order
    seen = set(); uniq = []
    for g in pool:
        if g not in seen:
            seen.add(g); uniq.append(g)
    if len(uniq) < n_genes:
        # Pad with plausible-looking synthetic symbols (clearly labelled) so the
        # matrix still has n_genes rows for the correlation loop.
        uniq += [f"SYN{i:04d}" for i in range(n_genes - len(uniq))]
    gene_names = uniq[:n_genes]
    # Baseline noise
    data = rng.standard_normal((n_genes, n_samples))
    # Inject structured co-expression: the first ~60 real genes are correlated
    # with the query gene so the pipeline produces visibly interesting results.
    qgene_row = data[0]
    for i in range(1, min(60, n_genes)):
        # Mix query gene signal with own noise; alternate sign for negative corr
        alpha = 0.3 + 0.4 * rng.random()
        sign = 1.0 if (i % 3) != 0 else -1.0
        data[i] = sign * alpha * qgene_row + np.sqrt(1 - alpha**2) * data[i]
    sample_names = [f"{sample_prefix}_{i:03d}" for i in range(n_samples)]
    return pd.DataFrame(data, index=gene_names, columns=sample_names)


def generate_synthetic_tcga_data(n_genes: int = 5000, n_samples: int = 200,
                                 query_gene: str = "TP53",
                                 cancer_type: Optional[str] = None) -> pd.DataFrame:
    """
    Generate synthetic TCGA-like expression data for demonstration.
    In production, this would fetch from GDC API.
    Always includes `query_gene` in the matrix so downstream lookup succeeds.
    """
    seed = _resolve_seed("tcga", cancer_type, None)
    tag = (cancer_type or "TCGA").upper()
    print(f"Generating synthetic TCGA-{tag} data: {n_genes} genes x {n_samples} samples "
          f"(includes query gene {query_gene}, seed={seed})...")
    print("  NOTE: synthetic placeholder matrix — not real TCGA expression. "
          "For production analyses, supply --dataset custom --expression-file <TSV>.")
    return _synthetic_matrix(n_genes, n_samples, query_gene, seed,
                             sample_prefix=f"TCGA-{tag}")


def generate_synthetic_gtex_data(n_genes: int = 5000, n_samples: int = 150,
                                 query_gene: str = "TP53",
                                 tissue: Optional[str] = None) -> pd.DataFrame:
    """
    Generate synthetic GTEx-like expression data for demonstration.
    Always includes `query_gene` in the matrix so downstream lookup succeeds.
    """
    seed = _resolve_seed("gtex", None, tissue)
    tag = (tissue or "GTEX").replace(" ", "_").upper()[:32]
    print(f"Generating synthetic GTEx-{tag} data: {n_genes} genes x {n_samples} samples "
          f"(includes query gene {query_gene}, seed={seed})...")
    print("  NOTE: synthetic placeholder matrix — not real GTEx expression. "
          "For production analyses, supply --dataset custom --expression-file <TSV>.")
    return _synthetic_matrix(n_genes, n_samples, query_gene, seed,
                             sample_prefix=f"GTEX-{tag}")


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

    # FDR correction (Benjamini-Hochberg with monotonicity enforcement)
    if len(results_df) > 0:
        pvals = results_df["pvalue"].values
        n = len(pvals)
        sort_idx = np.argsort(pvals)
        sorted_pvals = pvals[sort_idx]
        fdr_sorted = sorted_pvals * n / np.arange(1, n + 1)
        for i in range(n - 2, -1, -1):
            fdr_sorted[i] = min(fdr_sorted[i], fdr_sorted[i + 1])
        fdr_sorted = np.minimum(fdr_sorted, 1.0)
        fdr = np.empty(n)
        fdr[sort_idx] = fdr_sorted
        results_df["fdr"] = fdr
        results_df = results_df.sort_values("r", key=abs, ascending=False)

    return results_df


def _select_top_pos_neg(df: pd.DataFrame, top_n: int, r_col: str = "r") -> pd.DataFrame:
    """Select top_n/2 most positive and top_n/2 most negative correlations."""
    pos = df[df[r_col] > 0].sort_values(r_col, ascending=False)
    neg = df[df[r_col] < 0].sort_values(r_col, ascending=True)
    half = top_n // 2
    n_pos = min(len(pos), half)
    n_neg = min(len(neg), half)
    if n_pos < half:
        n_neg = min(len(neg), top_n - n_pos)
    elif n_neg < half:
        n_pos = min(len(pos), top_n - n_neg)
    return pd.concat([pos.head(n_pos), neg.head(n_neg)], ignore_index=True)


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

    top_genes = _select_top_pos_neg(results_df, top_n, r_col="r").copy()
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

    top_genes = _select_top_pos_neg(results_df, top_n, r_col="r").copy()

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
                       help="FDR threshold for significance (default: 0.01)")
    parser.add_argument("--min-corr", type=float, default=0.2,
                       help="Minimum |correlation| for significant genes (default: 0.2)")
    parser.add_argument("--network-top-n", type=int, default=30,
                       help="Top N genes for network visualization")
    parser.add_argument("--run-go", action="store_true",
                       help="Run GO enrichment analysis")
    parser.add_argument("--outdir", default=".", help="Output directory")

    args = parser.parse_args()

    init_style(
        font_family=getattr(args, 'font_family', None),
        font_size=getattr(args, 'font_size', None),
    )

    os.makedirs(args.outdir, exist_ok=True)

    # Load expression data
    if args.dataset == "custom":
        if not args.expression_file:
            print("Error: --expression-file required for custom dataset", file=sys.stderr)
            return 1
        expr_df = load_expression_matrix(args.expression_file)
    elif args.dataset == "tcga":
        expr_df = generate_synthetic_tcga_data(
            n_genes=5000, n_samples=200,
            query_gene=args.gene, cancer_type=args.cancer_type)
    else:  # gtex
        expr_df = generate_synthetic_gtex_data(
            n_genes=5000, n_samples=150,
            query_gene=args.gene, tissue=args.tissue)

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

    # Filter by FDR and minimum correlation
    fdr_results = results_df[
        (results_df["fdr"] <= args.fdr_cutoff) &
        (results_df["r"].abs() >= args.min_corr)
    ].copy()
    print(f"Significant genes (|r| >= {args.min_corr}, FDR <= {args.fdr_cutoff}): {len(fdr_results)}")

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
        if args.dataset == "tcga":
            f.write(f"Cancer type: {args.cancer_type}\n")
        elif args.dataset == "gtex":
            f.write(f"Tissue: {args.tissue}\n")
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
