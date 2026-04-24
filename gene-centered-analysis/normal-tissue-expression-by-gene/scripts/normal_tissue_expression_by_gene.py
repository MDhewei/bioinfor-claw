#!/usr/bin/env python3

import argparse
import math
import os
import re
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

import numpy as np
import matplotlib.pyplot as plt
import sys as _sys, os as _os
try:
    _sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), *(['..'] * 3), '_shared'))
    from plot_style import init_style
except ImportError:
    def init_style(**kw): pass  # graceful fallback if _shared not available
import pandas as pd
import requests


GTEX_API_BASE = "https://gtexportal.org/api/v2"
SESSION = requests.Session()
SESSION.headers.update(
    {
        "User-Agent": "normal-tissue-expression-for-gene/0.5",
        "Accept": "application/json",
    }
)


@dataclass
class GeneResult:
    query_gene: str
    resolved_gene_symbol: str
    gencode_id: Optional[str]
    tissue_df: pd.DataFrame
    expression_category: str


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def safe_get_json(
    url: str,
    params: Optional[Dict[str, Any]] = None,
    timeout: int = 60,
) -> Dict[str, Any]:
    r = SESSION.get(url, params=params, timeout=timeout)
    r.raise_for_status()
    return r.json()


def first_non_null(d: Dict[str, Any], keys: List[str]) -> Any:
    for k in keys:
        if k in d and d[k] is not None:
            return d[k]
    return None


def normalize_gene_symbol(gene: str) -> str:
    return re.sub(r"\s+", "", gene.strip()).upper()


def search_gene_in_gtex(gene_symbol: str) -> Dict[str, Optional[str]]:
    """
    Resolve a gene symbol using GTEx reference/gene.
    """
    url = f"{GTEX_API_BASE}/reference/gene"

    candidate_params = [
        {"geneId": gene_symbol, "format": "json"},
        {"geneId": gene_symbol},
    ]

    last_err = None
    items: List[Dict[str, Any]] = []

    for params in candidate_params:
        try:
            data = safe_get_json(url, params=params)
            items = data.get("data", []) or []
            if items:
                break
        except Exception as e:
            last_err = e

    if not items:
        raise RuntimeError(
            f"Could not resolve gene {gene_symbol} via GTEx reference/gene. Last error: {last_err}"
        )

    chosen = None
    for item in items:
        sym = first_non_null(item, ["geneSymbol", "symbol", "gene_name", "geneName"])
        if sym and str(sym).upper() == gene_symbol.upper():
            chosen = item
            break

    if chosen is None:
        chosen = items[0]

    resolved_symbol = first_non_null(
        chosen, ["geneSymbol", "symbol", "gene_name", "geneName"]
    ) or gene_symbol

    gencode_id = first_non_null(
        chosen, ["gencodeId", "geneId", "gencode_id"]
    )

    return {
        "resolved_gene_symbol": str(resolved_symbol),
        "gencode_id": str(gencode_id) if gencode_id is not None else None,
    }


def _normalize_median_items(items: List[Dict[str, Any]]) -> pd.DataFrame:
    rows = []
    for item in items:
        tissue = first_non_null(
            item,
            [
                "tissueSiteDetailId",
                "tissueSiteDetail",
                "tissueSiteDetailAbbr",
                "tissueSite",
                "tissue",
            ],
        )
        expr = first_non_null(
            item,
            [
                "median",
                "medianExpression",
                "expression",
                "value",
            ],
        )

        if tissue is None or expr is None:
            continue

        try:
            expr_value = float(expr)
        except Exception:
            continue

        rows.append({"tissue": str(tissue), "expression": expr_value})

    if not rows:
        return pd.DataFrame(columns=["tissue", "expression", "log2_tpm_plus1"])

    df = pd.DataFrame(rows)
    df = (
        df.groupby("tissue", as_index=False)["expression"]
        .mean()
        .sort_values("expression", ascending=False)
        .reset_index(drop=True)
    )
    # Add log2(TPM+1) transform — standard for GTEx visualisation
    df["log2_tpm_plus1"] = np.log2(df["expression"] + 1)
    return df


def fetch_median_tissue_expression(
    gencode_id: Optional[str],
    gene_symbol: str,
) -> pd.DataFrame:
    """
    Fetch GTEx median tissue expression.

    Preferred route:
      /expression/medianGeneExpression?gencodeIds=...&datasetId=gtex_v8

    If versioned GENCODE ID fails, also try the unversioned ID.
    """
    if not gencode_id:
        raise RuntimeError(f"No gencode_id available for gene {gene_symbol}")

    url = f"{GTEX_API_BASE}/expression/medianGeneExpression"

    candidate_gencode_ids = [gencode_id]
    if "." in gencode_id:
        candidate_gencode_ids.append(gencode_id.split(".")[0])

    last_http_err = None
    empty_attempts = []

    for gid in candidate_gencode_ids:
        candidate_params = [
            {"gencodeIds": gid, "datasetId": "gtex_v8", "page": 0, "itemsPerPage": 100000},
            {"gencodeId": gid, "datasetId": "gtex_v8", "page": 0, "itemsPerPage": 100000},
            {"gencodeIds": gid, "page": 0, "itemsPerPage": 100000},
            {"gencodeId": gid, "page": 0, "itemsPerPage": 100000},
        ]

        for params in candidate_params:
            try:
                data = safe_get_json(url, params=params)
                items = data.get("data", []) or []
                if not items:
                    empty_attempts.append(params)
                    continue

                df = _normalize_median_items(items)
                if not df.empty:
                    return df

                empty_attempts.append(params)

            except Exception as e:
                last_http_err = e
                continue

    if last_http_err is not None:
        raise RuntimeError(
            f"Could not retrieve GTEx median tissue expression for gene={gene_symbol}, "
            f"gencode_id={gencode_id}. Last HTTP error: {last_http_err}"
        )

    raise RuntimeError(
        f"GTEx returned no median tissue-expression rows for gene={gene_symbol}, "
        f"gencode_id={gencode_id}. Tried parameter sets: {empty_attempts}"
    )


def classify_expression_pattern(
    df: pd.DataFrame,
    expr_detect_threshold: float = 1.0,
    tissue_specific_max_threshold: float = 5.0,
    tissue_specific_fold_threshold: float = 5.0,
    tissue_specific_detect_frac_max: float = 0.30,
    universal_detect_frac_min: float = 0.80,
    universal_median_min: float = 1.0,
) -> str:
    """Classify expression pattern using log2(TPM+1) values."""
    if df.empty:
        return "no-data"

    # Use log2(TPM+1) for classification
    col = "log2_tpm_plus1" if "log2_tpm_plus1" in df.columns else "expression"
    expr = df[col].astype(float)

    # Convert thresholds to log2 scale if using log2 column
    if col == "log2_tpm_plus1":
        detect_thr = np.log2(expr_detect_threshold + 1)
        specific_max_thr = np.log2(tissue_specific_max_threshold + 1)
        univ_median_min = np.log2(universal_median_min + 1)
    else:
        detect_thr = expr_detect_threshold
        specific_max_thr = tissue_specific_max_threshold
        univ_median_min = universal_median_min

    max_expr = float(expr.max())
    median_expr = float(expr.median())
    n_total = len(expr)
    n_detected = int((expr >= detect_thr).sum())
    frac_detected = n_detected / n_total if n_total > 0 else 0.0

    if max_expr < detect_thr:
        return "non-expressed"

    fold_vs_median = max_expr / (median_expr + 1e-6)

    if (
        max_expr >= specific_max_thr
        and fold_vs_median >= tissue_specific_fold_threshold
        and frac_detected <= tissue_specific_detect_frac_max
    ):
        return "tissue-specific"

    if frac_detected >= universal_detect_frac_min and median_expr >= univ_median_min:
        return "universally-expressed"

    return "mixed"


def write_summary(
    result: GeneResult,
    out_txt: str,
    top_n: int = 10,
    expr_detect_threshold: float = 1.0,
) -> None:
    df = result.tissue_df.copy()

    lines = [
        f"Query gene: {result.query_gene}",
        f"Resolved gene symbol: {result.resolved_gene_symbol}",
        f"GENCODE ID: {result.gencode_id}",
        f"Expression category: {result.expression_category}",
        f"Number of tissues: {len(df)}",
    ]

    if not df.empty:
        raw_expr = df["expression"].astype(float)
        log2_expr = df["log2_tpm_plus1"].astype(float) if "log2_tpm_plus1" in df.columns else np.log2(raw_expr + 1)
        log2_detect_thr = np.log2(expr_detect_threshold + 1)
        n_detected = int((log2_expr >= log2_detect_thr).sum())
        frac_detected = n_detected / len(df)

        lines.append(f"Max expression: {raw_expr.max():.4f} TPM  ({log2_expr.max():.4f} log2(TPM+1))")
        lines.append(f"Median expression across tissues: {raw_expr.median():.4f} TPM  ({log2_expr.median():.4f} log2(TPM+1))")
        lines.append(
            f"Tissues above threshold ({expr_detect_threshold} TPM): "
            f"{n_detected}/{len(df)} ({frac_detected:.2%})"
        )
        lines.append("")
        lines.append(f"Top {min(top_n, len(df))} tissues (sorted by expression):")
        for i, row in enumerate(df.head(top_n).itertuples(index=False), start=1):
            log2_val = np.log2(row.expression + 1)
            lines.append(f"{i}. {row.tissue}: {row.expression:.4f} TPM  ({log2_val:.4f} log2(TPM+1))")

    with open(out_txt, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def make_barplot(
    df: pd.DataFrame,
    out_png: str,
    gene_symbol: str,
    expr_detect_threshold: float = 1.0,
) -> None:
    """
    Publication-quality horizontal barplot for all tissues.
    Uses log2(TPM+1) on the x-axis (standard for GTEx visualisation).
    """
    if df.empty:
        return

    plot_df = df.copy()
    # Ensure log2(TPM+1) column exists
    if "log2_tpm_plus1" not in plot_df.columns:
        plot_df["log2_tpm_plus1"] = np.log2(plot_df["expression"] + 1)
    plot_df = plot_df.sort_values("log2_tpm_plus1", ascending=True)

    log2_threshold = np.log2(expr_detect_threshold + 1)

    fig_w = 10
    fig_h = max(6, 0.28 * len(plot_df) + 1.8)

    plt.figure(figsize=(fig_w, fig_h), dpi=300)
    ax = plt.gca()

    ax.barh(
        plot_df["tissue"],
        plot_df["log2_tpm_plus1"],
        edgecolor="black",
        linewidth=0.4,
    )

    ax.axvline(
        x=log2_threshold,
        linestyle="--",
        linewidth=1.2,
        color="red",
        alpha=0.9,
        label=f"Detection threshold = {expr_detect_threshold:g} TPM",
    )

    ax.set_xlabel("Median expression — log2(TPM+1)", fontsize=12)
    ax.set_ylabel("Tissue", fontsize=12)
    ax.set_title(f"{gene_symbol} normal tissue expression (GTEx)", fontsize=14, pad=12)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)

    ax.tick_params(axis="x", labelsize=10, width=0.8)
    ax.tick_params(axis="y", labelsize=9, width=0.8)

    ax.grid(axis="x", linestyle=":", linewidth=0.5, alpha=0.5)
    ax.set_axisbelow(True)

    ax.legend(frameon=False, fontsize=10, loc="lower right")

    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")

    out_pdf = os.path.splitext(out_png)[0] + ".pdf"
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Retrieve GTEx normal tissue expression for a single gene."
    )
    parser.add_argument("--gene", required=True, help="Gene symbol, e.g. TP53")
    parser.add_argument(
        "--top-n",
        type=int,
        default=10,
        help="Top tissues to summarize in summary.txt",
    )
    parser.add_argument(
        "--expr-detect-threshold",
        type=float,
        default=1.0,
        help="Expression threshold for calling a tissue detected",
    )
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    init_style(
        font_family=getattr(args, 'font_family', None),
        font_size=getattr(args, 'font_size', None),
    )

    gene = normalize_gene_symbol(args.gene)
    ensure_dir(args.outdir)

    print(f"[INFO] Resolving gene: {gene}")
    gene_info = search_gene_in_gtex(gene)

    resolved_gene_symbol = gene_info["resolved_gene_symbol"]
    gencode_id = gene_info["gencode_id"]

    print(f"[INFO] Resolved gene symbol: {resolved_gene_symbol}")
    print(f"[INFO] GENCODE ID: {gencode_id}")

    print("[INFO] Fetching normal tissue expression from GTEx ...")
    tissue_df = fetch_median_tissue_expression(
        gencode_id=gencode_id,
        gene_symbol=resolved_gene_symbol,
    )

    expression_category = classify_expression_pattern(
        tissue_df,
        expr_detect_threshold=args.expr_detect_threshold,
    )

    tissue_df = tissue_df.copy()
    tissue_df.insert(0, "gene_symbol", resolved_gene_symbol)
    tissue_df.insert(1, "gencode_id", gencode_id)
    tissue_df["rank_within_gene"] = range(1, len(tissue_df) + 1)

    result = GeneResult(
        query_gene=gene,
        resolved_gene_symbol=resolved_gene_symbol,
        gencode_id=gencode_id,
        tissue_df=tissue_df,
        expression_category=expression_category,
    )

    tsv_path = os.path.join(args.outdir, f"{gene}.gtex_tissues.tsv")
    summary_path = os.path.join(args.outdir, "summary.txt")
    plot_path = os.path.join(args.outdir, f"{gene}.gtex_all_tissues.png")
    summary_tsv_path = os.path.join(args.outdir, "gene_expression_summary.tsv")

    tissue_df.to_csv(tsv_path, sep="\t", index=False)
    write_summary(
        result,
        summary_path,
        top_n=args.top_n,
        expr_detect_threshold=args.expr_detect_threshold,
    )
    # Ensure log2 column is present for plotting
    if "log2_tpm_plus1" not in tissue_df.columns:
        tissue_df["log2_tpm_plus1"] = np.log2(tissue_df["expression"].astype(float) + 1)

    make_barplot(
        tissue_df[["tissue", "expression", "log2_tpm_plus1"]].copy(),
        plot_path,
        gene_symbol=resolved_gene_symbol,
        expr_detect_threshold=args.expr_detect_threshold,
    )

    expr = tissue_df["expression"].astype(float)
    log2_expr = tissue_df["log2_tpm_plus1"].astype(float)
    log2_detect_thr = np.log2(args.expr_detect_threshold + 1)
    n_detected = int((log2_expr >= log2_detect_thr).sum())
    frac_detected = n_detected / len(tissue_df) if len(tissue_df) > 0 else 0.0

    summary_df = pd.DataFrame(
        [
            {
                "query_gene": gene,
                "resolved_gene_symbol": resolved_gene_symbol,
                "gencode_id": gencode_id,
                "expression_category": expression_category,
                "max_expression_TPM": float(expr.max()) if len(expr) else None,
                "max_expression_log2_TPM_plus1": float(log2_expr.max()) if len(log2_expr) else None,
                "median_expression_TPM": float(expr.median()) if len(expr) else None,
                "median_expression_log2_TPM_plus1": float(log2_expr.median()) if len(log2_expr) else None,
                "num_tissues": int(len(tissue_df)),
                "num_detected_tissues": n_detected,
                "fraction_detected_tissues": frac_detected,
            }
        ]
    )
    summary_df.to_csv(summary_tsv_path, sep="\t", index=False)

    # Print key results to stdout for agent consumption
    print(f"\n[RESULTS] Gene: {resolved_gene_symbol} | Category: {expression_category}")
    print(f"[RESULTS] Tissues: {len(tissue_df)} total, {n_detected} detected (threshold >= {args.expr_detect_threshold} TPM)")
    if len(expr) > 0:
        print(f"[RESULTS] Max expression: {float(expr.max()):.2f} TPM | Median: {float(expr.median()):.2f} TPM")
        # Top 5 tissues by expression
        top5 = tissue_df.nlargest(5, 'expression')
        for _, row in top5.iterrows():
            print(f"[RESULTS]   {row['tissue']}: {float(row['expression']):.2f} TPM")
    print(f"[DONE] Results written to: {os.path.dirname(tsv_path)}")


if __name__ == "__main__":
    main()