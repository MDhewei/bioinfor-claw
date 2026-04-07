#!/usr/bin/env python3

import argparse
import os
import re
import textwrap
from pathlib import Path
from typing import Dict, Tuple

import pandas as pd
import numpy as np
import gseapy as gp
from gseapy.plot import gseaplot


LIBRARY_MAP: Dict[str, Dict[str, str]] = {
    "HALLMARK": {
        "human": "MSigDB_Hallmark_2020",
        "mouse": "MSigDB_Hallmark_2020",
    },
    "KEGG": {
        "human": "KEGG_2021_Human",
        "mouse": "KEGG_2021_Mouse",
    },
    "REACTOME": {
        "human": "Reactome_2022",
        "mouse": "Reactome_2022",
    },
    "GO_BP": {
        "human": "GO_Biological_Process_2023",
        "mouse": "GO_Biological_Process_2023",
    },
    "GO_MF": {
        "human": "GO_Molecular_Function_2023",
        "mouse": "GO_Molecular_Function_2023",
    },
    "GO_CC": {
        "human": "GO_Cellular_Component_2023",
        "mouse": "GO_Cellular_Component_2023",
    },
}


ORGANISM_MAP = {
    "human": "human",
    "homo sapiens": "human",
    "hs": "human",
    "hsapiens": "human",
    "mouse": "mouse",
    "mus musculus": "mouse",
    "mm": "mouse",
}


def ensure_parent_dir(path: str) -> None:
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)


def normalize_organism(org: str) -> str:
    key = str(org).strip().lower()
    if key not in ORGANISM_MAP:
        raise ValueError(
            f"Unsupported organism '{org}'. Supported examples: human, mouse, hsapiens, mm"
        )
    return ORGANISM_MAP[key]


def resolve_library(lib: str, organism: str) -> str:
    lib = lib.strip().upper()
    if lib not in LIBRARY_MAP:
        raise ValueError(f"Unsupported library: {lib}")
    if organism not in LIBRARY_MAP[lib]:
        raise ValueError(f"No library configured for {lib} and organism {organism}")
    return LIBRARY_MAP[lib][organism]


def normalize_gene_symbol(gene: str) -> str:
    gene = str(gene).strip()
    gene = re.sub(r"\s+", "", gene)
    return gene.upper()


def detect_gene_col(df: pd.DataFrame, user_col: str = None) -> str:
    if user_col:
        if user_col not in df.columns:
            raise ValueError(f"Gene column '{user_col}' not found")
        return user_col

    candidates = ["gene", "Gene", "symbol", "Symbol", "GeneSymbol", "gene_symbol"]
    for c in candidates:
        if c in df.columns:
            return c

    lower_map = {c.lower(): c for c in df.columns}
    for c in ["gene", "symbol", "genesymbol", "gene_symbol"]:
        if c in lower_map:
            return lower_map[c]

    raise ValueError("Could not detect gene column automatically")


def detect_score_col(df: pd.DataFrame, user_col: str = None) -> str:
    if user_col:
        if user_col not in df.columns:
            raise ValueError(f"Score column '{user_col}' not found")
        return user_col

    candidates = [
        "score", "Score", "rank", "Rank",
        "log2FC", "logFC", "stat", "Statistic", "t", "T",
        "wald", "Wald", "correlation", "Correlation"
    ]
    for c in candidates:
        if c in df.columns:
            return c

    lower_map = {c.lower(): c for c in df.columns}
    for c in ["score", "rank", "log2fc", "logfc", "stat", "wald", "correlation"]:
        if c in lower_map:
            return lower_map[c]

    raise ValueError("Could not detect score column automatically")


def read_ranked_gene_table(path: str, gene_col: str = None, score_col: str = None) -> Tuple[pd.DataFrame, str, str]:
    df = pd.read_csv(path, sep=None, engine="python")
    if len(df) == 0:
        raise ValueError("Input ranked gene file is empty")

    gcol = detect_gene_col(df, gene_col)
    scol = detect_score_col(df, score_col)

    out = df[[gcol, scol]].copy()
    out[gcol] = out[gcol].map(normalize_gene_symbol)
    out[scol] = pd.to_numeric(out[scol], errors="coerce")
    out = out.dropna(subset=[gcol, scol]).copy()
    out = out[out[gcol] != ""].copy()

    # keep highest absolute score for duplicated genes
    out["_abs_score"] = out[scol].abs()
    out = out.sort_values("_abs_score", ascending=False).drop_duplicates(subset=[gcol], keep="first")
    out = out.drop(columns="_abs_score")
    out = out.sort_values(scol, ascending=False).reset_index(drop=True)

    if len(out) == 0:
        raise ValueError("No valid ranked genes remain after cleaning")

    return out, gcol, scol


def find_result_column(df: pd.DataFrame, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    lower_map = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower_map:
            return lower_map[c.lower()]
    return None


def clean_prerank_results(df: pd.DataFrame) -> pd.DataFrame:
    if df is None or len(df) == 0:
        return pd.DataFrame()

    out = df.copy()

    term_col = find_result_column(out, ["Term", "Name", "pathway"])
    nes_col = find_result_column(out, ["NES"])
    fdr_col = find_result_column(out, ["FDR q-val", "FDR q-value", "FDR", "fdr"])
    es_col = find_result_column(out, ["ES"])
    pval_col = find_result_column(out, ["NOM p-val", "P-value", "pval"])

    rename_map = {}
    if term_col and term_col != "Term":
        rename_map[term_col] = "Term"
    if nes_col and nes_col != "NES":
        rename_map[nes_col] = "NES"
    if fdr_col and fdr_col != "FDR q-val":
        rename_map[fdr_col] = "FDR q-val"
    if es_col and es_col != "ES":
        rename_map[es_col] = "ES"
    if pval_col and pval_col != "NOM p-val":
        rename_map[pval_col] = "NOM p-val"

    out = out.rename(columns=rename_map)

    required = ["Term", "NES", "FDR q-val"]
    missing = [c for c in required if c not in out.columns]
    if missing:
        raise ValueError(f"Missing expected GSEA result columns: {missing}")

    out["NES"] = pd.to_numeric(out["NES"], errors="coerce")
    out["FDR q-val"] = pd.to_numeric(out["FDR q-val"], errors="coerce")
    if "ES" in out.columns:
        out["ES"] = pd.to_numeric(out["ES"], errors="coerce")
    if "NOM p-val" in out.columns:
        out["NOM p-val"] = pd.to_numeric(out["NOM p-val"], errors="coerce")

    out = out.dropna(subset=["Term", "NES", "FDR q-val"]).copy()
    out["minus_log10_fdr"] = -np.log10(out["FDR q-val"].clip(lower=np.nextafter(0, 1)))
    out["Direction"] = np.where(out["NES"] >= 0, "positive", "negative")
    out = out.sort_values(["FDR q-val", "NES"], ascending=[True, False]).reset_index(drop=True)
    return out


def run_prerank(
    rnk_df: pd.DataFrame,
    gene_col: str,
    score_col: str,
    gene_sets: str,
    permutation_num: int,
    min_size: int,
    max_size: int,
    seed: int = 123,
):
    prerank_input = rnk_df[[gene_col, score_col]].copy()

    pre_res = gp.prerank(
        rnk=prerank_input,
        gene_sets=gene_sets,
        permutation_num=permutation_num,
        min_size=min_size,
        max_size=max_size,
        seed=seed,
        outdir=None,
        no_plot=True,
        verbose=False,
    )
    return pre_res


def select_top_terms(df: pd.DataFrame, fdr_cutoff: float, top_n: int) -> pd.DataFrame:
    if len(df) == 0:
        return df

    sig = df[df["FDR q-val"] <= fdr_cutoff].copy()
    if len(sig) == 0:
        sig = df.copy()

    pos = sig[sig["NES"] > 0].sort_values(
        ["FDR q-val", "NES"], ascending=[True, False]
    ).head(top_n)

    neg = sig[sig["NES"] < 0].sort_values(
        ["FDR q-val", "NES"], ascending=[True, True]
    ).head(top_n)

    out = pd.concat([pos, neg], axis=0).drop_duplicates(subset=["Term"]).reset_index(drop=True)
    return out


def make_classic_gsea_plots(
    pre_res,
    res_df: pd.DataFrame,
    output_prefix: str,
    top_n: int = 5,
    fdr_cutoff: float = 0.25,
):
    if len(res_df) == 0:
        raise ValueError("No GSEA results available for plotting")

    sig = res_df[res_df["FDR q-val"] <= fdr_cutoff].copy()
    if len(sig) == 0:
        sig = res_df.copy()

    pos = sig[sig["NES"] > 0].sort_values(
        ["FDR q-val", "NES"], ascending=[True, False]
    ).head(top_n)

    neg = sig[sig["NES"] < 0].sort_values(
        ["FDR q-val", "NES"], ascending=[True, True]
    ).head(top_n)

    selected = pd.concat([pos, neg], axis=0).drop_duplicates(subset=["Term"]).reset_index(drop=True)

    plot_paths = []

    # gseapy stores ranking as a Series/DataFrame-like object
    rank_metric = pre_res.ranking

    for _, row in selected.iterrows():
        term = row["Term"]
        if term not in pre_res.results:
            continue

        safe_term = re.sub(r"[^A-Za-z0-9._-]+", "_", str(term))[:120]
        out_png = f"{output_prefix}.{safe_term}.gsea.png"
        out_pdf = f"{output_prefix}.{safe_term}.gsea.pdf"

        ensure_parent_dir(out_png)

        term_result = pre_res.results[term]

        gseaplot(
            rank_metric=rank_metric,
            term=term,
            ofname=out_png,
            **term_result,
        )

        gseaplot(
            rank_metric=rank_metric,
            term=term,
            ofname=out_pdf,
            **term_result,
        )

        plot_paths.append((term, out_png, out_pdf))

    return plot_paths


def save_plot_input_table(df: pd.DataFrame, output: str) -> None:
    cols = ["Term", "NES", "ES", "NOM p-val", "FDR q-val", "minus_log10_fdr", "Direction"]
    existing = [c for c in cols if c in df.columns]
    ensure_parent_dir(output)
    df[existing].to_csv(output, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(
        description="Preranked GSEA with classic enrichment plots"
    )
    parser.add_argument("--input", required=True, help="Input ranked gene file")
    parser.add_argument("--output-prefix", required=True, help="Output file prefix")
    parser.add_argument("--gene-col", default=None)
    parser.add_argument("--score-col", default=None)
    parser.add_argument("--organism", default="human")
    parser.add_argument("--library", default="HALLMARK",
                        choices=["HALLMARK", "KEGG", "REACTOME", "GO_BP", "GO_MF", "GO_CC"])
    parser.add_argument("--permutation-num", type=int, default=1000)
    parser.add_argument("--min-size", type=int, default=5)
    parser.add_argument("--max-size", type=int, default=500)
    parser.add_argument("--top-n", type=int, default=5)
    parser.add_argument("--fdr-cutoff", type=float, default=0.25)
    args = parser.parse_args()

    organism = normalize_organism(args.organism)
    gene_set = resolve_library(args.library, organism)

    ranked_df, gene_col, score_col = read_ranked_gene_table(
        args.input,
        gene_col=args.gene_col,
        score_col=args.score_col,
    )

    print(f"[INFO] Ranked genes: {len(ranked_df)}")
    print(f"[INFO] Gene column: {gene_col}")
    print(f"[INFO] Score column: {score_col}")
    print(f"[INFO] Library: {gene_set}")

    pre_res = run_prerank(
        rnk_df=ranked_df,
        gene_col=gene_col,
        score_col=score_col,
        gene_sets=gene_set,
        permutation_num=args.permutation_num,
        min_size=args.min_size,
        max_size=args.max_size,
    )

    res_df = clean_prerank_results(pre_res.res2d)

    if len(res_df) == 0:
        raise RuntimeError("No valid GSEA results were returned")

    full_tsv = f"{args.output_prefix}.{args.library}.full.tsv"
    sig_tsv = f"{args.output_prefix}.{args.library}.significant.tsv"
    plot_tsv = f"{args.output_prefix}.{args.library}.plot_input.tsv"

    ensure_parent_dir(full_tsv)
    res_df.to_csv(full_tsv, sep="\t", index=False)

    sig_df = res_df[res_df["FDR q-val"] <= args.fdr_cutoff].copy()
    if len(sig_df) == 0:
        sig_df = res_df.copy()
    sig_df.to_csv(sig_tsv, sep="\t", index=False)

    plot_df = select_top_terms(res_df, fdr_cutoff=args.fdr_cutoff, top_n=args.top_n)
    save_plot_input_table(plot_df, plot_tsv)

    classic_prefix = f"{args.output_prefix}.{args.library}"
    plot_paths = make_classic_gsea_plots(
        pre_res=pre_res,
        res_df=res_df,
        output_prefix=classic_prefix,
        top_n=args.top_n,
        fdr_cutoff=args.fdr_cutoff,
    )

    print(f"[INFO] Saved full result table: {full_tsv}")
    print(f"[INFO] Saved significant result table: {sig_tsv}")
    print(f"[INFO] Saved plot input table: {plot_tsv}")
    for term, png, pdf in plot_paths:
        print(f"[INFO] Saved classic GSEA plot for {term}: {png}")
        print(f"[INFO] Saved classic GSEA plot for {term}: {pdf}")

    if len(plot_paths) == 0:
        raise RuntimeError("No classic GSEA plots were generated.")


if __name__ == "__main__":
    main()