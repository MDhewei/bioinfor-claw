#!/usr/bin/env python3
"""
crispr_screen_analysis.py
─────────────────────────
CRISPR pooled screen analysis using MAGeCK.

Modes
-----
test  : MAGeCK RRA — pairwise treatment vs control
mle   : MAGeCK MLE — multi-condition / time-course
count : MAGeCK count — FASTQ → count table
all   : count then test (full pipeline from FASTQs)

If MAGeCK is not installed the script falls back to a pure-Python
implementation of the RRA algorithm for the 'test' mode so that the
workflow always produces usable output.
"""

import argparse
import os
import sys
import subprocess
import shutil
import json
import warnings
from pathlib import Path

import math
import numpy as np
import pandas as pd
import matplotlib
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), '..', '..', '..', '_shared'))
from plot_style import init_style
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# scipy is optional — use it when available, fall back to pure-Python otherwise
try:
    from scipy import stats as _scipy_stats
    _HAVE_SCIPY = True
except ImportError:
    _HAVE_SCIPY = False

warnings.filterwarnings("ignore")

# ── colours ───────────────────────────────────────────────────────────────────
NEG_COLOR  = "#2166AC"   # blue  – negative selection
POS_COLOR  = "#D6604D"   # red   – positive selection
GRAY_COLOR = "#AAAAAA"
HIT_ALPHA  = 0.85
DOT_ALPHA  = 0.35
FIG_DPI    = 150

plt.rcParams.update({
    "figure.dpi":      FIG_DPI,
    "font.family":     "sans-serif",
    "axes.spines.top": False,
    "axes.spines.right": False,
})


# ─────────────────────────────────────────────────────────────────────────────
#  MAGeCK detection
# ─────────────────────────────────────────────────────────────────────────────

def _mageck_available() -> bool:
    return shutil.which("mageck") is not None


def _run(cmd: list, label: str = "") -> str:
    """Run a shell command, stream stdout/stderr, return stdout."""
    print(f"\n[RUN] {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stderr[-2000:])
        raise RuntimeError(f"Command failed ({label}): {result.returncode}")
    return result.stdout


# ─────────────────────────────────────────────────────────────────────────────
#  Python-native RRA (fallback when MAGeCK is absent)
# ─────────────────────────────────────────────────────────────────────────────

def _beta_cdf(x: float, a: float, b: float) -> float:
    """Regularised incomplete beta function via continued fraction (pure Python)."""
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    # Use scipy when available for accuracy; otherwise use a simple approximation
    if _HAVE_SCIPY:
        return float(_scipy_stats.beta.cdf(x, a, b))
    # Pure-Python regularised incomplete beta via numerical integration (Simpson)
    n_steps = 200
    h = x / n_steps
    def integrand(t):
        if t <= 0:
            return 0.0
        return (t ** (a - 1)) * ((1 - t) ** (b - 1))
    xs = [i * h for i in range(n_steps + 1)]
    ys = [integrand(xi) for xi in xs]
    integral = h / 3 * (ys[0] + ys[-1] + 4 * sum(ys[i] for i in range(1, n_steps, 2))
                                        + 2 * sum(ys[i] for i in range(2, n_steps, 2)))
    # Normalise by beta function B(a, b) via log-gamma
    log_beta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    return min(1.0, integral / math.exp(log_beta))


def _rra_score(ranks: np.ndarray, n_total: int, alpha: float = 0.1) -> float:
    """
    Robust Rank Aggregation score.
    ranks  : sorted ascending rank fractions (0..1) for guides of a gene
    n_total: total number of guides in library
    alpha  : RRA threshold parameter
    """
    k = len(ranks)
    if k == 0:
        return 1.0
    scores = []
    for i, r in enumerate(sorted(ranks)):
        p = _beta_cdf(r, i + 1, k - i)
        scores.append(p)
    return float(min(scores))


def _lfc_per_gene(df: pd.DataFrame, treat_cols: list, ctrl_cols: list,
                  pseudocount: float = 1.0, method: str = "median") -> pd.Series:
    """Median (or alpha-median) log2 fold-change per guide, then gene-level."""
    treat = df[treat_cols].mean(axis=1) + pseudocount
    ctrl  = df[ctrl_cols].mean(axis=1)  + pseudocount
    df = df.copy()
    df["lfc"] = np.log2(treat / ctrl)
    if method == "alphamedian":
        gene_lfc = df.groupby("gene")["lfc"].apply(
            lambda x: np.mean(sorted(x)[len(x)//5: -len(x)//5 or None])
        )
    else:
        gene_lfc = df.groupby("gene")["lfc"].median()
    return gene_lfc


def _python_rra_test(count_table: str, treatment: list, control: list,
                     norm_method: str, gene_lfc_method: str, prefix: str,
                     outdir: Path) -> pd.DataFrame:
    """
    Pure-Python RRA-style analysis.
    Returns gene_summary DataFrame with columns:
        gene, num_sgrna, neg.lfc, neg.score, neg.fdr, pos.lfc, pos.score, pos.fdr
    """
    print("[INFO] MAGeCK not found — running Python RRA fallback.")
    df = pd.read_csv(count_table, sep="\t")

    # Detect column structure
    if "sgRNA" in df.columns:
        df = df.rename(columns={"sgRNA": "sgrna"})
    if "Gene" in df.columns:
        df = df.rename(columns={"Gene": "gene"})
    # lowercase first two columns as fallback
    cols = df.columns.tolist()
    if cols[0].lower() not in ("sgrna", "id", "sgrna_id"):
        df = df.rename(columns={cols[0]: "sgrna", cols[1]: "gene"})

    all_samples = treatment + control

    # ── Normalise ──────────────────────────────────────────────────────────
    count_cols = df[all_samples]
    if norm_method == "median":
        size_factors = count_cols.sum(axis=0) / count_cols.sum(axis=0).median()
        df[all_samples] = count_cols.div(size_factors, axis=1)
    elif norm_method == "total":
        df[all_samples] = count_cols.div(count_cols.sum(axis=0), axis=1) * 1e6

    # ── sgRNA-level LFC ───────────────────────────────────────────────────
    eps = 1.0
    treat_mean = df[treatment].mean(axis=1) + eps
    ctrl_mean  = df[control].mean(axis=1)   + eps
    df["lfc"]  = np.log2(treat_mean / ctrl_mean)

    n_total = len(df)

    # ── Rank guides (neg = most depleted, pos = most enriched) ───────────
    df["rank_neg"] = df["lfc"].rank(method="average") / n_total    # low rank = depleted
    df["rank_pos"] = df["lfc"].rank(method="average", ascending=False) / n_total

    # ── Gene-level RRA ────────────────────────────────────────────────────
    gene_lfc = _lfc_per_gene(df, treatment, control, method=gene_lfc_method)

    records = []
    for gene, grp in df.groupby("gene"):
        n = len(grp)
        neg_score = _rra_score(grp["rank_neg"].values, n_total)
        pos_score = _rra_score(grp["rank_pos"].values, n_total)
        records.append({
            "gene":      gene,
            "num_sgrna": n,
            "neg.lfc":   gene_lfc.get(gene, 0),
            "neg.score": neg_score,
            "pos.lfc":   gene_lfc.get(gene, 0),
            "pos.score": pos_score,
        })

    gene_df = pd.DataFrame(records)

    # ── BH FDR ───────────────────────────────────────────────────────────
    def _bh_fdr(pvals):
        """Benjamini-Hochberg FDR correction (pure Python / numpy)."""
        pv = np.array(pvals, dtype=float)
        n = len(pv)
        order = np.argsort(pv)
        ranks = np.empty(n, dtype=int)
        ranks[order] = np.arange(1, n + 1)
        fdr = np.minimum(1.0, pv * n / ranks)
        # Enforce monotonicity from right
        for i in range(n - 2, -1, -1):
            fdr[order[i]] = min(fdr[order[i]], fdr[order[i + 1]])
        return fdr

    for col in ("neg.score", "pos.score"):
        gene_df[col.replace("score", "fdr")] = _bh_fdr(gene_df[col].values)

    gene_df = gene_df.sort_values("neg.score")

    out_path = outdir / f"{prefix}.gene_summary.txt"
    gene_df.to_csv(out_path, sep="\t", index=False)
    print(f"[OK] Gene summary saved: {out_path}")

    # sgRNA summary
    sgrna_out = outdir / f"{prefix}.sgrna_summary.txt"
    df[["sgrna", "gene", "lfc", "rank_neg", "rank_pos"]].to_csv(
        sgrna_out, sep="\t", index=False)

    return gene_df, df


# ─────────────────────────────────────────────────────────────────────────────
#  MAGeCK wrappers
# ─────────────────────────────────────────────────────────────────────────────

def run_mageck_count(fastqs: list, library: str, treatment: list, control: list,
                     prefix: str, outdir: Path, norm_method: str) -> str:
    """Run mageck count. Returns path to generated count table."""
    # MAGeCK count requires sample labels and corresponding fastqs
    all_labels = treatment + control
    if len(all_labels) != len(fastqs):
        raise ValueError(
            f"Number of FASTQ files ({len(fastqs)}) must equal number of "
            f"sample labels ({len(all_labels)}). "
            f"Labels: {all_labels}"
        )
    sample_str = ",".join(all_labels)
    fastq_str  = " ".join(fastqs)

    cmd = [
        "mageck", "count",
        "-l", library,
        "-n", str(outdir / prefix),
        "--sample-label", sample_str,
        "--fastq"
    ] + fastqs + [
        "--norm-method", norm_method,
    ]
    _run(cmd, "mageck count")
    count_path = str(outdir / f"{prefix}.count.txt")
    print(f"[OK] Count table: {count_path}")
    return count_path


def run_mageck_test(count_table: str, treatment: list, control: list,
                    prefix: str, outdir: Path, norm_method: str,
                    gene_lfc_method: str) -> pd.DataFrame:
    """Run mageck test (RRA). Returns gene_summary DataFrame."""
    cmd = [
        "mageck", "test",
        "-k", count_table,
        "-t", ",".join(treatment),
        "-c", ",".join(control),
        "-n", str(outdir / prefix),
        "--norm-method", norm_method,
        "--gene-lfc-method", gene_lfc_method,
    ]
    _run(cmd, "mageck test")
    gene_path = outdir / f"{prefix}.gene_summary.txt"
    gene_df = pd.read_csv(gene_path, sep="\t")
    sgrna_df = pd.read_csv(outdir / f"{prefix}.sgrna_summary.txt", sep="\t")
    return gene_df, sgrna_df


def run_mageck_mle(count_table: str, design_matrix: str,
                   prefix: str, outdir: Path, norm_method: str) -> pd.DataFrame:
    """Run mageck mle (beta score). Returns gene_summary DataFrame."""
    cmd = [
        "mageck", "mle",
        "-k", count_table,
        "-d", design_matrix,
        "-n", str(outdir / prefix),
        "--norm-method", norm_method,
    ]
    _run(cmd, "mageck mle")
    gene_path = outdir / f"{prefix}.gene_summary.txt"
    gene_df = pd.read_csv(gene_path, sep="\t")
    return gene_df, None


# ─────────────────────────────────────────────────────────────────────────────
#  QC plot
# ─────────────────────────────────────────────────────────────────────────────

def plot_qc(count_table: str, treatment: list, control: list,
            prefix: str, outdir: Path):
    """Read count distribution and Gini index per sample."""
    df = pd.read_csv(count_table, sep="\t")
    cols = df.columns.tolist()
    if cols[0].lower() not in ("sgrna", "id"):
        df = df.rename(columns={cols[0]: "sgrna", cols[1]: "gene"})

    all_samples = treatment + control
    sample_cols = [c for c in all_samples if c in df.columns]
    if not sample_cols:
        print("[WARN] Could not match sample names to count table columns — skipping QC plot.")
        return

    def gini(arr):
        arr = np.sort(arr.astype(float) + 1)
        n = len(arr)
        idx = np.arange(1, n + 1)
        return (2 * (idx * arr).sum()) / (n * arr.sum()) - (n + 1) / n

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    colors = plt.cm.tab10(np.linspace(0, 1, len(sample_cols)))

    # Log-scale count distribution
    ax = axes[0]
    for col, color in zip(sample_cols, colors):
        vals = df[col].values.astype(float) + 1
        vals = vals[vals > 0]
        ax.hist(np.log10(vals), bins=50, alpha=0.55, color=color, label=col)
    ax.set_xlabel("log10(read count + 1)")
    ax.set_ylabel("Number of sgRNAs")
    ax.set_title("Read count distribution")
    ax.legend(fontsize=7, ncol=2)

    # Gini index
    ax = axes[1]
    gini_vals = [gini(df[c].values) for c in sample_cols]
    bars = ax.bar(range(len(sample_cols)), gini_vals, color=colors)
    ax.set_xticks(range(len(sample_cols)))
    ax.set_xticklabels(sample_cols, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Gini index")
    ax.set_title("Library evenness (lower = more even)")
    ax.axhline(0.1, ls="--", color="gray", lw=1, label="threshold 0.1")
    ax.legend(fontsize=8)

    fig.tight_layout()
    out = outdir / f"{prefix}_qc.png"
    fig.savefig(out, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] QC plot: {out}")


# ─────────────────────────────────────────────────────────────────────────────
#  Volcano plot
# ─────────────────────────────────────────────────────────────────────────────

def plot_volcano(gene_df: pd.DataFrame, fdr_cutoff: float, lfc_cutoff: float,
                 top_n: int, prefix: str, outdir: Path, mode: str = "test"):
    """Gene-level volcano: LFC vs -log10(FDR) for both neg and pos arms."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    for ax, arm, color, label in [
        (axes[0], "neg", NEG_COLOR, "Negative selection (dropout)"),
        (axes[1], "pos", POS_COLOR, "Positive selection (enrichment)"),
    ]:
        if mode == "mle":
            # MLE output columns differ per condition — use first beta column
            beta_cols  = [c for c in gene_df.columns if c.endswith("|beta")]
            fdr_cols   = [c for c in gene_df.columns if c.endswith("|fdr")]
            if not beta_cols:
                ax.text(0.5, 0.5, "MLE: no beta column found",
                        ha="center", transform=ax.transAxes)
                continue
            lfc_col = beta_cols[0]
            fdr_col = fdr_cols[0] if fdr_cols else None
            xvals = gene_df[lfc_col]
            yvals = -np.log10(gene_df[fdr_col].clip(1e-300)) if fdr_col else np.zeros(len(gene_df))
            hit_mask = (gene_df[fdr_col] < fdr_cutoff) if fdr_col else pd.Series([False] * len(gene_df))
        else:
            lfc_col = f"{arm}.lfc"
            fdr_col = f"{arm}.fdr"
            if lfc_col not in gene_df.columns:
                continue
            xvals = gene_df[lfc_col]
            yvals = -np.log10(gene_df[fdr_col].clip(1e-300))
            hit_mask = (gene_df[fdr_col] < fdr_cutoff) & (xvals.abs() >= lfc_cutoff)

        # Background
        ax.scatter(xvals[~hit_mask], yvals[~hit_mask],
                   c=GRAY_COLOR, s=18, alpha=DOT_ALPHA, linewidths=0, rasterized=True)
        # Hits
        ax.scatter(xvals[hit_mask], yvals[hit_mask],
                   c=color, s=30, alpha=HIT_ALPHA, linewidths=0, zorder=3)

        # Labels for top hits
        hit_genes = gene_df[hit_mask].copy()
        if arm == "neg":
            hit_genes = hit_genes.nsmallest(top_n, lfc_col)
        else:
            hit_genes = hit_genes.nlargest(top_n, lfc_col)

        for _, row in hit_genes.iterrows():
            ax.annotate(
                row["gene"],
                xy=(row[lfc_col], -np.log10(max(row[fdr_col], 1e-300))),
                fontsize=6.5, color=color,
                xytext=(4, 0), textcoords="offset points",
            )

        # Threshold lines
        ax.axhline(-np.log10(fdr_cutoff), ls="--", lw=0.8, color="black", alpha=0.5)
        if arm != "mle":
            ax.axvline(-lfc_cutoff if arm == "neg" else lfc_cutoff,
                       ls="--", lw=0.8, color="black", alpha=0.5)

        ax.set_xlabel(f"log2 Fold Change ({arm})")
        ax.set_ylabel(f"-log10(FDR)")
        ax.set_title(label)
        n_hits = hit_mask.sum()
        ax.text(0.02, 0.97, f"Hits (FDR<{fdr_cutoff}): {n_hits}",
                transform=ax.transAxes, va="top", fontsize=8)

    fig.suptitle(f"{prefix} — Volcano Plot", fontweight="bold")
    fig.tight_layout()
    out = outdir / f"{prefix}_volcano.png"
    fig.savefig(out, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Volcano plot: {out}")


# ─────────────────────────────────────────────────────────────────────────────
#  Rank plot
# ─────────────────────────────────────────────────────────────────────────────

def plot_rank(gene_df: pd.DataFrame, fdr_cutoff: float, top_n: int,
              prefix: str, outdir: Path):
    """Rank plot for negative and positive selection."""
    for arm, color, title, lfc_col, fdr_col, ascending in [
        ("neg", NEG_COLOR, "Negative selection rank", "neg.lfc", "neg.fdr", True),
        ("pos", POS_COLOR, "Positive selection rank", "pos.lfc", "pos.fdr", False),
    ]:
        if lfc_col not in gene_df.columns:
            continue

        df = gene_df.sort_values(lfc_col, ascending=ascending).reset_index(drop=True)
        df["rank"] = df.index + 1
        hit_mask = df[fdr_col] < fdr_cutoff

        fig, ax = plt.subplots(figsize=(9, 5))
        ax.scatter(df["rank"][~hit_mask], df[lfc_col][~hit_mask],
                   c=GRAY_COLOR, s=14, alpha=DOT_ALPHA, linewidths=0, rasterized=True)
        ax.scatter(df["rank"][hit_mask], df[lfc_col][hit_mask],
                   c=color, s=22, alpha=HIT_ALPHA, linewidths=0, zorder=3)

        # Label top genes
        top = df[hit_mask].head(top_n)
        for _, row in top.iterrows():
            ax.annotate(row["gene"], xy=(row["rank"], row[lfc_col]),
                        fontsize=6, color=color,
                        xytext=(3, 2), textcoords="offset points")

        ax.set_xlabel("Gene rank")
        ax.set_ylabel("log2 Fold Change")
        ax.set_title(f"{prefix} — {title}")
        ax.axhline(0, ls="--", lw=0.7, color="black", alpha=0.4)

        patch = mpatches.Patch(color=color, label=f"FDR < {fdr_cutoff}: {hit_mask.sum()} genes")
        ax.legend(handles=[patch], fontsize=8)

        fig.tight_layout()
        out = outdir / f"{prefix}_rank_{arm}.png"
        fig.savefig(out, dpi=FIG_DPI, bbox_inches="tight")
        plt.close(fig)
        print(f"[OK] Rank plot ({arm}): {out}")


# ─────────────────────────────────────────────────────────────────────────────
#  Hit summary bar chart
# ─────────────────────────────────────────────────────────────────────────────

def plot_hit_summary(gene_df: pd.DataFrame, fdr_cutoff: float, top_n: int,
                     prefix: str, outdir: Path):
    """Horizontal bar chart of top hits, colored by neg/pos selection."""
    if "neg.fdr" not in gene_df.columns:
        return

    neg_hits = gene_df[gene_df["neg.fdr"] < fdr_cutoff].nsmallest(top_n, "neg.lfc")
    pos_hits = gene_df[gene_df["pos.fdr"] < fdr_cutoff].nlargest(top_n, "pos.lfc")

    records = []
    for _, row in neg_hits.iterrows():
        records.append({"gene": row["gene"], "lfc": row["neg.lfc"], "arm": "neg"})
    for _, row in pos_hits.iterrows():
        records.append({"gene": row["gene"], "lfc": row["pos.lfc"], "arm": "pos"})

    if not records:
        print("[INFO] No hits to plot in summary chart.")
        return

    df = pd.DataFrame(records).sort_values("lfc")
    colors = [NEG_COLOR if a == "neg" else POS_COLOR for a in df["arm"]]

    fig, ax = plt.subplots(figsize=(7, max(4, len(df) * 0.32)))
    bars = ax.barh(range(len(df)), df["lfc"], color=colors, alpha=0.85, edgecolor="none")
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["gene"], fontsize=8)
    ax.set_xlabel("log2 Fold Change")
    ax.set_title(f"{prefix} — Top screen hits (FDR < {fdr_cutoff})")
    ax.axvline(0, color="black", lw=0.7)

    neg_patch = mpatches.Patch(color=NEG_COLOR, label="Negative selection")
    pos_patch = mpatches.Patch(color=POS_COLOR, label="Positive selection")
    ax.legend(handles=[neg_patch, pos_patch], fontsize=8, loc="lower right")

    fig.tight_layout()
    out = outdir / f"{prefix}_hit_summary.png"
    fig.savefig(out, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Hit summary: {out}")


# ─────────────────────────────────────────────────────────────────────────────
#  sgRNA scatter
# ─────────────────────────────────────────────────────────────────────────────

def plot_sgrna_scatter(sgrna_df: pd.DataFrame, gene_df: pd.DataFrame,
                       fdr_cutoff: float, top_n: int, prefix: str, outdir: Path):
    """Scatter plot of per-guide LFC, colored by hit status."""
    if sgrna_df is None or "lfc" not in sgrna_df.columns:
        return

    # Determine hit genes
    hit_genes_neg = set()
    hit_genes_pos = set()
    if "neg.fdr" in gene_df.columns:
        hit_genes_neg = set(gene_df[gene_df["neg.fdr"] < fdr_cutoff]["gene"])
    if "pos.fdr" in gene_df.columns:
        hit_genes_pos = set(gene_df[gene_df["pos.fdr"] < fdr_cutoff]["gene"])

    df = sgrna_df.copy()
    df["color"] = GRAY_COLOR
    df.loc[df["gene"].isin(hit_genes_neg), "color"] = NEG_COLOR
    df.loc[df["gene"].isin(hit_genes_pos), "color"] = POS_COLOR

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.scatter(range(len(df)), df["lfc"].sort_values().values,
               c=df.loc[df["lfc"].sort_values().index, "color"].values,
               s=8, alpha=0.5, linewidths=0, rasterized=True)
    ax.axhline(0, color="black", lw=0.7, ls="--")
    ax.set_xlabel("sgRNA rank")
    ax.set_ylabel("log2 Fold Change")
    ax.set_title(f"{prefix} — sgRNA-level LFC")

    neg_patch = mpatches.Patch(color=NEG_COLOR, label=f"Negative hits ({len(hit_genes_neg)} genes)")
    pos_patch = mpatches.Patch(color=POS_COLOR, label=f"Positive hits ({len(hit_genes_pos)} genes)")
    gray_patch = mpatches.Patch(color=GRAY_COLOR, label="Non-significant")
    ax.legend(handles=[neg_patch, pos_patch, gray_patch], fontsize=8)

    fig.tight_layout()
    out = outdir / f"{prefix}_sgrna_scatter.png"
    fig.savefig(out, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] sgRNA scatter: {out}")


# ─────────────────────────────────────────────────────────────────────────────
#  Save hit table
# ─────────────────────────────────────────────────────────────────────────────

def save_hits(gene_df: pd.DataFrame, fdr_cutoff: float, prefix: str, outdir: Path):
    """Export filtered hit table (FDR < cutoff in either arm)."""
    if "neg.fdr" not in gene_df.columns and "pos.fdr" not in gene_df.columns:
        return

    mask = pd.Series([False] * len(gene_df), index=gene_df.index)
    if "neg.fdr" in gene_df.columns:
        mask |= gene_df["neg.fdr"] < fdr_cutoff
    if "pos.fdr" in gene_df.columns:
        mask |= gene_df["pos.fdr"] < fdr_cutoff

    hits = gene_df[mask].copy()
    out = outdir / f"{prefix}_hits.tsv"
    hits.to_csv(out, sep="\t", index=False)
    print(f"[OK] Hit table ({len(hits)} genes, FDR<{fdr_cutoff}): {out}")
    return hits


# ─────────────────────────────────────────────────────────────────────────────
#  Main
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="CRISPR screen analysis via MAGeCK",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--count-table",      help="Count table TSV (sgRNA × samples)")
    p.add_argument("--treatment",        help="Treatment sample columns, comma-separated")
    p.add_argument("--control",          help="Control sample columns, comma-separated")
    p.add_argument("--mode",             default="test",
                   choices=["test", "mle", "count", "all"],
                   help="Analysis mode (default: test)")
    p.add_argument("--fastq",            nargs="+",
                   help="FASTQ files (for count/all mode)")
    p.add_argument("--library",          help="sgRNA library TSV (for count/all mode)")
    p.add_argument("--design-matrix",    help="Design matrix TSV (for mle mode)")
    p.add_argument("--norm-method",      default="median",
                   choices=["median", "total", "control", "none"],
                   help="Normalisation method (default: median)")
    p.add_argument("--fdr-cutoff",       type=float, default=0.25,
                   help="FDR significance threshold (default: 0.25)")
    p.add_argument("--lfc-cutoff",       type=float, default=1.0,
                   help="|LFC| threshold for volcano labelling (default: 1.0)")
    p.add_argument("--gene-lfc-method",  default="median",
                   choices=["median", "alphamedian"],
                   help="Gene-level LFC method (default: median)")
    p.add_argument("--top-n",            type=int, default=20,
                   help="Top N hits to label in plots (default: 20)")
    p.add_argument("--prefix",           default="screen",
                   help="Output file prefix (default: screen)")
    p.add_argument("--outdir",           default="./crispr_screen_results",
                   help="Output directory (default: ./crispr_screen_results/)")
    return p.parse_args()


def main():
    args = parse_args()
    init_style(
        font_family=getattr(args, 'font_family', None),
        font_size=getattr(args, 'font_size', None),
    )
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    treatment = [t.strip() for t in args.treatment.split(",")] if args.treatment else []
    control   = [c.strip() for c in args.control.split(",")]   if args.control   else []

    use_mageck = _mageck_available()
    print(f"[INFO] MAGeCK {'detected' if use_mageck else 'NOT found — using Python RRA fallback'}")
    print(f"[INFO] Mode: {args.mode}")
    print(f"[INFO] Output: {outdir}/")

    gene_df  = None
    sgrna_df = None
    count_table = args.count_table

    # ── Step 1: count (FASTQ → count table) ──────────────────────────────────
    if args.mode in ("count", "all"):
        if not use_mageck:
            print("[ERROR] mageck count requires MAGeCK to be installed.")
            print("  Install via: conda install -c bioconda mageck")
            sys.exit(1)
        if not args.fastq or not args.library:
            print("[ERROR] --fastq and --library are required for count/all mode.")
            sys.exit(1)
        count_table = run_mageck_count(
            args.fastq, args.library, treatment, control,
            args.prefix, outdir, args.norm_method,
        )
        if args.mode == "count":
            print(f"\n[DONE] Count table written to: {count_table}")
            return

    # ── Step 2: test or mle ───────────────────────────────────────────────────
    if not count_table:
        print("[ERROR] --count-table is required for test/mle mode.")
        sys.exit(1)

    if args.mode in ("test", "all"):
        if not treatment or not control:
            print("[ERROR] --treatment and --control are required for test mode.")
            sys.exit(1)
        if use_mageck:
            gene_df, sgrna_df = run_mageck_test(
                count_table, treatment, control,
                args.prefix, outdir, args.norm_method, args.gene_lfc_method,
            )
        else:
            gene_df, sgrna_df = _python_rra_test(
                count_table, treatment, control,
                args.norm_method, args.gene_lfc_method, args.prefix, outdir,
            )

    elif args.mode == "mle":
        if not args.design_matrix:
            print("[ERROR] --design-matrix is required for mle mode.")
            sys.exit(1)
        if use_mageck:
            gene_df, sgrna_df = run_mageck_mle(
                count_table, args.design_matrix,
                args.prefix, outdir, args.norm_method,
            )
        else:
            print("[ERROR] MLE mode requires MAGeCK. Install: conda install -c bioconda mageck")
            sys.exit(1)

    if gene_df is None:
        print("[ERROR] No gene summary produced.")
        sys.exit(1)

    # ── Plots ──────────────────────────────────────────────────────────────────
    print("\n[INFO] Generating plots ...")

    # QC (only if we have the count table locally)
    if count_table and os.path.isfile(count_table) and treatment and control:
        try:
            plot_qc(count_table, treatment, control, args.prefix, outdir)
        except Exception as e:
            print(f"[WARN] QC plot failed: {e}")

    mode_label = "mle" if args.mode == "mle" else "test"
    plot_volcano(gene_df, args.fdr_cutoff, args.lfc_cutoff, args.top_n,
                 args.prefix, outdir, mode=mode_label)

    if args.mode in ("test", "all"):
        plot_rank(gene_df, args.fdr_cutoff, args.top_n, args.prefix, outdir)
        plot_hit_summary(gene_df, args.fdr_cutoff, args.top_n, args.prefix, outdir)
        if sgrna_df is not None:
            plot_sgrna_scatter(sgrna_df, gene_df, args.fdr_cutoff,
                               args.top_n, args.prefix, outdir)

    hits = save_hits(gene_df, args.fdr_cutoff, args.prefix, outdir)

    # ── Save metadata ─────────────────────────────────────────────────────────
    meta = {
        "mode":           args.mode,
        "treatment":      treatment,
        "control":        control,
        "norm_method":    args.norm_method,
        "fdr_cutoff":     args.fdr_cutoff,
        "lfc_cutoff":     args.lfc_cutoff,
        "mageck_used":    use_mageck,
        "n_genes":        int(len(gene_df)),
        "n_hits":         int(len(hits)) if hits is not None else 0,
    }
    with open(outdir / f"{args.prefix}_meta.json", "w") as f:
        json.dump(meta, f, indent=2)

    # ── Summary ───────────────────────────────────────────────────────────────
    print("\n" + "═" * 56)
    print(f"  CRISPR screen analysis complete")
    print("═" * 56)
    print(f"  Mode          : {args.mode}")
    print(f"  MAGeCK used   : {use_mageck}")
    print(f"  Genes tested  : {len(gene_df)}")
    if hits is not None:
        print(f"  Hits (FDR<{args.fdr_cutoff}): {len(hits)}")
    print(f"  Output dir    : {outdir}/")
    print("═" * 56 + "\n")


if __name__ == "__main__":
    main()
