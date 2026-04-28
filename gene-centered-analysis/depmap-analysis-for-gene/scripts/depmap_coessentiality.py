#!/usr/bin/env python3
"""DepMap co-essentiality analysis for a single gene.

Computes Pearson correlations between the query gene's CRISPR dependency
score and all other genes' dependency scores across DepMap cell lines.
Genes that are co-essential tend to work in the same pathway or complex.

DATA LOADING: Streams the essentiality matrix directly from the DepMap API
in two passes — never loads the full multi-GB file into memory.

Produces:
  - Ranked correlation table with FDR correction (TSV)
  - Top co-essential genes horizontal bar plot (PNG + PDF)
  - Co-essentiality network visualization (PNG + PDF)
"""

from __future__ import annotations

import argparse
import csv
import io
import math
import os
import re
import sys
from typing import Dict, List, Optional, Tuple
from urllib.request import urlopen

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
        if close: plt.close(fig)


# ═══════════════════════════════════════════════════════════
# DepMap streaming infrastructure
# ═══════════════════════════════════════════════════════════
DEPMAP_RELEASE = "DepMap Public 26Q1"
DEPMAP_INDEX = "https://depmap.org/portal/api/download/files"
ESSENTIALITY_FILE = "CRISPRGeneEffect.csv"

_file_url_cache: Dict[str, str] = {}


def depmap_file_url(filename: str) -> str:
    if filename in _file_url_cache:
        return _file_url_cache[filename]
    index_data = urlopen(DEPMAP_INDEX, timeout=120).read().decode("utf-8")
    for row in csv.DictReader(io.StringIO(index_data)):
        if row.get("release") == DEPMAP_RELEASE and row.get("filename") == filename:
            _file_url_cache[filename] = row["url"]
            return row["url"]
    raise RuntimeError(f"Missing DepMap file: {filename}")


def stream_csv_url(url: str):
    with urlopen(url, timeout=600) as response:
        text = io.TextIOWrapper(response, encoding="utf-8", newline="")
        reader = csv.reader(text)
        for row in reader:
            yield row


# ═══════════════════════════════════════════════════════════
# scipy-free statistics
# ═══════════════════════════════════════════════════════════
def _betacf(a, b, x):
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
    if x <= 0: return 0.0
    if x >= 1: return 1.0
    ln_pre = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
              + a * math.log(x) + b * math.log(1 - x))
    front = math.exp(ln_pre)
    if x < (a + 1) / (a + b + 2):
        return front * _betacf(a, b, x) / a
    else:
        return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def pearson_r_and_p(n, sum_x, sum_y, sum_xy, sum_x2, sum_y2):
    if n < 3:
        return 0.0, 1.0
    denom_x = n * sum_x2 - sum_x * sum_x
    denom_y = n * sum_y2 - sum_y * sum_y
    if denom_x <= 0 or denom_y <= 0:
        return 0.0, 1.0
    r = (n * sum_xy - sum_x * sum_y) / math.sqrt(denom_x * denom_y)
    r = max(-1.0, min(1.0, r))
    if abs(r) == 1.0:
        return r, 0.0
    t2 = r * r * (n - 2) / (1.0 - r * r)
    df = n - 2
    p = _betainc(df / 2.0, 0.5, df / (df + t2))
    return r, p


# ═══════════════════════════════════════════════════════════
# Two-pass streaming correlation
# ═══════════════════════════════════════════════════════════
def find_gene_col_index(header, gene):
    for i, h in enumerate(header):
        if h == gene or h.startswith(f"{gene} "):
            return i
    raise ValueError(f"Gene {gene} not found in header. First 10: {header[:10]}")


def extract_symbol(col):
    m = re.match(r"^([A-Za-z0-9_.-]+)", str(col))
    return m.group(1).upper() if m else str(col).upper()


def streaming_correlations(url, gene, method="pearson"):
    # Pass 1: extract target gene column
    print(f"[INFO] Pass 1: extracting {gene} dependency scores...")
    rows_iter = stream_csv_url(url)
    header = next(rows_iter)
    gene_col_idx = find_gene_col_index(header, gene)
    model_col = header.index("ModelID") if "ModelID" in header else 0

    target_values = {}
    for row in rows_iter:
        try:
            target_values[row[model_col]] = float(row[gene_col_idx])
        except (ValueError, IndexError):
            continue

    n_models = len(target_values)
    print(f"[INFO] Got {n_models} models for {gene}")
    if n_models < 10:
        raise ValueError(f"Too few models ({n_models}) for {gene}")

    # Pass 2: accumulate running stats
    print(f"[INFO] Pass 2: computing correlations against all genes...")
    rows_iter = stream_csv_url(url)
    header = next(rows_iter)

    gene_columns = []
    for i, h in enumerate(header):
        if i == model_col or i == gene_col_idx:
            continue
        if h and (h[0].isalpha() or h[0].isdigit()):
            gene_columns.append((i, h))

    n_genes = len(gene_columns)
    print(f"[INFO] Will correlate against {n_genes} gene columns")

    acc_n = np.zeros(n_genes, dtype=np.int32)
    acc_sx = np.zeros(n_genes, dtype=np.float64)
    acc_sy = np.zeros(n_genes, dtype=np.float64)
    acc_sxy = np.zeros(n_genes, dtype=np.float64)
    acc_sx2 = np.zeros(n_genes, dtype=np.float64)
    acc_sy2 = np.zeros(n_genes, dtype=np.float64)

    rows_processed = 0
    for row in rows_iter:
        model_id = row[model_col]
        if model_id not in target_values:
            continue
        x = target_values[model_id]
        for j, (col_idx, _) in enumerate(gene_columns):
            try:
                y = float(row[col_idx])
            except (ValueError, IndexError):
                continue
            acc_n[j] += 1
            acc_sx[j] += x
            acc_sy[j] += y
            acc_sxy[j] += x * y
            acc_sx2[j] += x * x
            acc_sy2[j] += y * y
        rows_processed += 1
        if rows_processed % 500 == 0:
            print(f"[INFO] Processed {rows_processed}/{n_models} cell lines...")

    print(f"[INFO] Computing final correlations for {n_genes} genes...")
    results = []
    for j, (col_idx, col_name) in enumerate(gene_columns):
        n = int(acc_n[j])
        if n < 10:
            continue
        r, p = pearson_r_and_p(n, acc_sx[j], acc_sy[j], acc_sxy[j], acc_sx2[j], acc_sy2[j])
        if math.isnan(r):
            continue
        results.append({
            "gene_column": col_name,
            "gene_symbol": extract_symbol(col_name),
            "correlation": r,
            "pvalue": p,
            "n_samples": n,
        })

    df = pd.DataFrame(results)
    if len(df) == 0:
        return df

    # FDR (Benjamini-Hochberg)
    pvals = df["pvalue"].values
    n = len(pvals)
    sort_idx = np.argsort(pvals)
    sorted_pvals = pvals[sort_idx]
    fdr_sorted = sorted_pvals * n / np.arange(1, n + 1)
    for i in range(n - 2, -1, -1):
        fdr_sorted[i] = min(fdr_sorted[i], fdr_sorted[i + 1])
    fdr_sorted = np.minimum(fdr_sorted, 1.0)
    fdr = np.empty(n)
    fdr[sort_idx] = fdr_sorted
    df["fdr"] = fdr
    df = df.sort_values("correlation", key=abs, ascending=False).reset_index(drop=True)
    return df


# ═══════════════════════════════════════════════════════════
# Plotting
# ═══════════════════════════════════════════════════════════
def _select_top_pos_neg(df, top_n):
    pos = df[df["correlation"] > 0].sort_values("correlation", ascending=False)
    neg = df[df["correlation"] < 0].sort_values("correlation", ascending=True)
    half = top_n // 2
    n_pos = min(len(pos), half)
    n_neg = min(len(neg), half)
    if n_pos < half:
        n_neg = min(len(neg), top_n - n_pos)
    elif n_neg < half:
        n_pos = min(len(pos), top_n - n_neg)
    return pd.concat([pos.head(n_pos), neg.head(n_neg)], ignore_index=True)


def plot_top_coessential(df, gene, top_n=30, output_prefix="coessentiality_barplot"):
    if len(df) < 1:
        return
    top = _select_top_pos_neg(df, top_n).sort_values("correlation")
    fig_h = max(5, len(top) * 0.28 + 1.5)
    fig, ax = plt.subplots(figsize=(9, fig_h))
    colors = ["#d62728" if r > 0 else "#1f77b4" for r in top["correlation"]]
    ax.barh(range(len(top)), top["correlation"], color=colors, edgecolor="black", linewidth=0.4)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top["gene_symbol"], fontsize=9)
    ax.set_xlabel("Correlation coefficient (CRISPR gene effect)")
    ax.set_title(f"Top {len(top)} genes co-essential with {gene} (DepMap)")
    ax.axvline(0, color="black", linewidth=0.8)
    ax.grid(axis="x", alpha=0.3, linestyle="--")
    save_fig(fig, f"{output_prefix}.png")
    save_fig(fig, f"{output_prefix}.pdf", close=False)
    plt.close(fig)


def plot_network(df, gene, top_n=30, output_prefix="coessentiality_network"):
    if len(df) < 2:
        return
    top = _select_top_pos_neg(df, top_n).copy()
    np.random.seed(42)
    gene_list = list(top["gene_symbol"])
    positions = {g: np.random.randn(2) for g in gene_list}
    for _ in range(15):
        forces = {g: np.zeros(2) for g in gene_list}
        for i, g1 in enumerate(gene_list):
            for j, g2 in enumerate(gene_list):
                if i >= j: continue
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
            if i >= j: continue
            if abs(corr_vals[g1] - corr_vals[g2]) < 0.3:
                edges.append((g1, g2, 1 - abs(corr_vals[g1] - corr_vals[g2])))
    fig, ax = plt.subplots(figsize=(11, 9))
    for g1, g2, w in edges:
        ax.plot([positions[g1][0], positions[g2][0]],
                [positions[g1][1], positions[g2][1]],
                "gray", alpha=0.25, linewidth=w * 1.5, zorder=1)
    center = np.array([0.5, 0.5])
    ax.scatter(*center, s=600, c="#FFD700", edgecolors="black", linewidth=1.5, zorder=4)
    ax.text(*center, gene, ha="center", va="center", fontsize=9, fontweight="bold", zorder=5)
    for g in gene_list:
        r = corr_vals[g]
        pos = positions[g]
        color = "#d62728" if r > 0 else "#1f77b4"
        ax.scatter(pos[0], pos[1], s=abs(r) * 400 + 50, c=color,
                   alpha=min(0.95, abs(r) + 0.4), edgecolors="black", linewidth=0.8, zorder=2)
        ax.text(pos[0], pos[1] + 0.025, g, ha="center", va="bottom", fontsize=7, fontweight="bold", zorder=3)
        ax.plot([center[0], pos[0]], [center[1], pos[1]],
                color=color, alpha=0.35, linewidth=abs(r) * 2.5, zorder=0)
    ax.set_xlim(-0.1, 1.1); ax.set_ylim(-0.1, 1.1)
    ax.set_title(f"Co-essentiality network: {gene} (DepMap, top {len(gene_list)})")
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


# ═══════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════
def main():
    parser = argparse.ArgumentParser(
        description="DepMap co-essentiality analysis for a single gene. "
                    "Streams data from DepMap API — no full dataset download needed."
    )
    parser.add_argument("--gene", required=True, help="Gene symbol (e.g. TP53)")
    parser.add_argument("--method", choices=["pearson", "spearman"],
                        default="pearson", help="Correlation method (default: pearson)")
    parser.add_argument("--top-n", type=int, default=30)
    parser.add_argument("--fdr-cutoff", type=float, default=0.01)
    parser.add_argument("--min-corr", type=float, default=0.2)
    parser.add_argument("--network-top-n", type=int, default=30)
    parser.add_argument("--outdir", default=".", help="Output directory")
    parser.add_argument("--font-family", default=None)
    parser.add_argument("--font-size", type=float, default=None)
    # Legacy flags — accepted but IGNORED
    parser.add_argument("--essentiality-file", help=argparse.SUPPRESS)
    parser.add_argument("--metadata-file", help=argparse.SUPPRESS)
    args = parser.parse_args()

    init_style(font_family=args.font_family, font_size=args.font_size)
    os.makedirs(args.outdir, exist_ok=True)
    gene = re.sub(r"\s+", "", args.gene.strip()).upper()

    if args.method == "spearman":
        print("[WARN] Spearman requires ranking — using Pearson for streaming mode.")
        method = "pearson"
    else:
        method = args.method

    url = depmap_file_url(ESSENTIALITY_FILE)
    print(f"[INFO] Streaming co-essentiality for {gene} from DepMap API...")
    corr_df = streaming_correlations(url, gene, method=method)

    if len(corr_df) == 0:
        print(f"[ERROR] No valid correlations computed for {gene}.", file=sys.stderr)
        sys.exit(1)

    full_path = os.path.join(args.outdir, f"{gene}.depmap_coessentiality_full.tsv")
    corr_df.to_csv(full_path, sep="\t", index=False)

    sig = corr_df[
        (corr_df["fdr"] <= args.fdr_cutoff) &
        (corr_df["correlation"].abs() >= args.min_corr)
    ]
    sig_path = os.path.join(args.outdir, f"{gene}.depmap_coessentiality_sig.tsv")
    sig.to_csv(sig_path, sep="\t", index=False)

    plot_top_coessential(
        corr_df, gene, top_n=args.top_n,
        output_prefix=os.path.join(args.outdir, f"{gene}.depmap_coessentiality_barplot"),
    )
    plot_network(
        corr_df, gene, top_n=args.network_top_n,
        output_prefix=os.path.join(args.outdir, f"{gene}.depmap_coessentiality_network"),
    )

    print(f"\n[RESULTS] === DepMap Co-essentiality: {gene} ===")
    print(f"[RESULTS] Total genes correlated: {len(corr_df)}")
    print(f"[RESULTS] Significant (|r|>={args.min_corr}, FDR<={args.fdr_cutoff}): {len(sig)}")
    if len(corr_df) > 0:
        for _, row in corr_df.head(5).iterrows():
            print(f"[RESULTS]   {row['gene_symbol']}: r={row['correlation']:.3f}, FDR={row['fdr']:.2e}")
    print(f"[RESULTS] === END ===")
    print(f"[DONE] Results written to: {args.outdir}")


if __name__ == "__main__":
    main()
