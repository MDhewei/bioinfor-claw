#!/usr/bin/env python3

import argparse
import os
import re
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

import matplotlib.pyplot as plt
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), '..', '..', '..', '_shared'))
from plot_style import init_style
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
        return pd.DataFrame(columns=["tissue", "expression"])

    df = pd.DataFrame(rows)
    df = (
        df.groupby("tissue", as_index=False)["expression"]
        .mean()
        .sort_values("expression", ascending=False)
        .reset_index(drop=True)
    )
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
    if df.empty:
        return "no-data"

    expr = df["expression"].astype(float)
    max_expr = float(expr.max())
    median_expr = float(expr.median())
    n_total = len(expr)
    n_detected = int((expr >= expr_detect_threshold).sum())
    frac_detected = n_detected / n_total if n_total > 0 else 0.0

    if max_expr < expr_detect_threshold:
        return "non-expressed"

    fold_vs_median = max_expr / (median_expr + 1e-6)

    if (
        max_expr >= tissue_specific_max_threshold
        and fold_vs_median >= tissue_specific_fold_threshold
        and frac_detected <= tissue_specific_detect_frac_max
    ):
        return "tissue-specific"

    if frac_detected >= universal_detect_frac_min and median_expr >= universal_median_min:
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
        expr = df["expression"].astype(float)
        n_detected = int((expr >= expr_detect_threshold).sum())
        frac_detected = n_detected / len(df)

        lines.append(f"Max expression: {expr.max():.4f}")
        lines.append(f"Median expression across tissues: {expr.median():.4f}")
        lines.append(
            f"Tissues above threshold ({expr_detect_threshold}): "
            f"{n_detected}/{len(df)} ({frac_detected:.2%})"
        )
        lines.append("")
        lines.append(f"Top {min(top_n, len(df))} tissues:")
        for i, row in enumerate(df.head(top_n).itertuples(index=False), start=1):
            lines.append(f"{i}. {row.tissue}: {row.expression:.4f}")

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
    Because expression is on the x-axis, the threshold is drawn as a vertical line.
    """
    if df.empty:
        return

    plot_df = df.copy().sort_values("expression", ascending=True)

    fig_w = 10
    fig_h = max(6, 0.28 * len(plot_df) + 1.8)

    plt.figure(figsize=(fig_w, fig_h), dpi=300)
    ax = plt.gca()

    ax.barh(
        plot_df["tissue"],
        plot_df["expression"],
        edgecolor="black",
        linewidth=0.4,
    )

    ax.axvline(
        x=expr_detect_threshold,
        linestyle="--",
        linewidth=1.2,
        color="red",
        alpha=0.9,
        label=f"Detection threshold = {expr_detect_threshold:g}",
    )

    ax.set_xlabel("Median expression", fontsize=12)
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
    make_barplot(
        tissue_df[["tissue", "expression"]].copy(),
        plot_path,
        gene_symbol=resolved_gene_symbol,
        expr_detect_threshold=args.expr_detect_threshold,
    )

    expr = tissue_df["expression"].astype(float)
    n_detected = int((expr >= args.expr_detect_threshold).sum())
    frac_detected = n_detected / len(tissue_df) if len(tissue_df) > 0 else 0.0

    summary_df = pd.DataFrame(
        [
            {
                "query_gene": gene,
                "resolved_gene_symbol": resolved_gene_symbol,
                "gencode_id": gencode_id,
                "expression_category": expression_category,
                "max_expression": float(expr.max()) if len(expr) else None,
                "median_expression_across_tissues": float(expr.median()) if len(expr) else None,
                "num_tissues": int(len(tissue_df)),
                "num_detected_tissues": n_detected,
                "fraction_detected_tissues": frac_detected,
            }
        ]
    )
    summary_df.to_csv(summary_tsv_path, sep="\t", index=False)

    print(f"[DONE] Tissue table: {tsv_path}")
    print(f"[DONE] Summary text: {summary_path}")
    print(f"[DONE] Summary table: {summary_tsv_path}")
    print(f"[DONE] Barplot PNG: {plot_path}")
    print(f"[DONE] Barplot PDF: {os.path.splitext(plot_path)[0] + '.pdf'}")
    print(f"[DONE] Expression category: {expression_category}")


if __name__ == "__main__":
    main()