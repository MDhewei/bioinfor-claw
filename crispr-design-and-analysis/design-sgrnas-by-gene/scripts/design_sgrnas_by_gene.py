#!/usr/bin/env python3

import argparse
import os
from pathlib import Path
from typing import List, Optional

import pandas as pd


RAW_URLS = {
    "human": "https://raw.githubusercontent.com/MDhewei/GuidePro/master/Downloads/Genome-wide-sgRNA-Selection-human.csv",
    "monkey": "https://raw.githubusercontent.com/MDhewei/GuidePro/master/Downloads/Genome-wide-sgRNA-Selection-monkey.csv",
    "mouse": "https://raw.githubusercontent.com/MDhewei/GuidePro/master/Downloads/Genome-wide-sgRNA-Selection-mouse.csv",
}


def read_gene_list(gene: Optional[str], gene_file: Optional[str]) -> List[str]:
    genes: List[str] = []
    if gene:
        genes.extend([x.strip() for x in gene.split(",") if x.strip()])
    if gene_file:
        with open(gene_file) as f:
            genes.extend([x.strip() for x in f if x.strip()])
    genes = list(dict.fromkeys(genes))
    if not genes:
        raise ValueError("Please provide --gene or --gene-file")
    return genes


def ensure_dir(path: str) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


def detect_gene_col(df: pd.DataFrame, user_col: Optional[str] = None) -> str:
    if user_col:
        if user_col not in df.columns:
            raise ValueError(f"Gene column '{user_col}' not found")
        return user_col

    candidates = [
        "gene", "Gene", "symbol", "Symbol", "GeneSymbol", "gene_symbol",
        "Gene symbol", "Gene Symbol"
    ]
    for c in candidates:
        if c in df.columns:
            return c

    lower_map = {c.lower(): c for c in df.columns}
    for c in ["gene", "symbol", "genesymbol", "gene_symbol"]:
        if c in lower_map:
            return lower_map[c]

    raise ValueError(
        "Could not detect gene column automatically. Use --gene-col to specify it."
    )


def detect_score_col(df: pd.DataFrame, user_col: Optional[str] = None) -> Optional[str]:
    if user_col:
        if user_col not in df.columns:
            raise ValueError(f"Score column '{user_col}' not found")
        return user_col

    candidates = [
        "GuidePro score", "GuidePro_score", "GuideProScore",
        "score", "Score", "sgRNA_score", "sgRNA score"
    ]
    for c in candidates:
        if c in df.columns:
            return c

    lower_map = {c.lower(): c for c in df.columns}
    for c in ["guidepro score", "guidepro_score", "guideproscore", "score"]:
        if c in lower_map:
            return lower_map[c]

    return None


def detect_sgrna_col(df: pd.DataFrame, user_col: Optional[str] = None) -> Optional[str]:
    if user_col:
        if user_col not in df.columns:
            raise ValueError(f"sgRNA column '{user_col}' not found")
        return user_col

    candidates = [
        "sgRNA", "sgrna", "guide", "guide_seq", "Guide sequence",
        "sequence", "Spacer", "spacer", "sgRNA sequence"
    ]
    for c in candidates:
        if c in df.columns:
            return c

    lower_map = {c.lower(): c for c in df.columns}
    for c in ["sgrna", "guide", "guide_seq", "sequence", "spacer"]:
        if c in lower_map:
            return lower_map[c]

    return None


def get_cache_file(cache_dir: str, genome: str) -> str:
    return os.path.join(cache_dir, f"guidepro_{genome}.csv")


def download_if_needed(genome: str, cache_dir: str, refresh: bool = False) -> str:
    if genome not in RAW_URLS:
        raise ValueError(f"Unsupported genome: {genome}")

    ensure_dir(cache_dir)
    cache_file = get_cache_file(cache_dir, genome)

    if os.path.exists(cache_file) and not refresh:
        return cache_file

    url = RAW_URLS[genome]
    print(f"[INFO] Downloading {genome} GuidePro reference...")
    df = pd.read_csv(url)
    df.to_csv(cache_file, index=False)
    print(f"[INFO] Saved cache to: {cache_file}")
    return cache_file


def filter_guides(
    df: pd.DataFrame,
    genes: List[str],
    gene_col: str,
    score_col: Optional[str],
    sgrna_col: Optional[str],
    top_n_per_gene: int,
    drop_duplicates: bool,
) -> pd.DataFrame:
    genes_upper = {g.upper() for g in genes}
    out = df.copy()

    out["_gene_upper"] = out[gene_col].astype(str).str.upper()
    out = out[out["_gene_upper"].isin(genes_upper)].copy()

    if len(out) == 0:
        return out

    if score_col is not None:
        out[score_col] = pd.to_numeric(out[score_col], errors="coerce")

    if drop_duplicates and sgrna_col is not None:
        out = out.drop_duplicates(subset=[gene_col, sgrna_col], keep="first").copy()

    if score_col is not None:
        out = out.sort_values([gene_col, score_col], ascending=[True, False]).copy()
    else:
        out = out.sort_values([gene_col]).copy()

    if top_n_per_gene is not None and top_n_per_gene > 0:
        out = (
            out.groupby(gene_col, dropna=False, group_keys=False)
            .head(top_n_per_gene)
            .copy()
        )

    out["query_source"] = "GuidePro"
    return out


def write_summary(
    out_df: pd.DataFrame,
    genes: List[str],
    gene_col: str,
    summary_file: Optional[str],
) -> None:
    if summary_file is None:
        return

    requested = len(genes)
    found_genes = out_df[gene_col].astype(str).nunique() if len(out_df) else 0
    returned_rows = len(out_df)

    summary = pd.DataFrame([{
        "requested_genes": requested,
        "genes_with_results": found_genes,
        "returned_rows": returned_rows,
    }])

    Path(summary_file).parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_file, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(
        description="Retrieve sgRNAs by gene or gene list from GuidePro genome-wide reference tables"
    )
    parser.add_argument("--gene", default=None, help="Comma-separated genes, e.g. TP53,EGFR,EP300")
    parser.add_argument("--gene-file", default=None, help="One gene per line")
    parser.add_argument("--genome", choices=["human", "monkey", "mouse"], default="human")
    parser.add_argument("--cache-dir", default="./guidepro_cache")
    parser.add_argument("--refresh-cache", action="store_true")
    parser.add_argument("--gene-col", default=None, help="Optional gene column name")
    parser.add_argument("--score-col", default=None, help="Optional score column name")
    parser.add_argument("--sgrna-col", default=None, help="Optional sgRNA column name")
    parser.add_argument("--top-n-per-gene", type=int, default=10)
    parser.add_argument("--drop-duplicates", action="store_true")
    parser.add_argument("--output", required=True, help="Output TSV")
    parser.add_argument("--summary", default=None, help="Optional summary TSV")
    args = parser.parse_args()

    genes = read_gene_list(args.gene, args.gene_file)
    cache_file = download_if_needed(
        genome=args.genome,
        cache_dir=args.cache_dir,
        refresh=args.refresh_cache,
    )

    print(f"[INFO] Reading reference table: {cache_file}")
    df = pd.read_csv(cache_file)

    gene_col = detect_gene_col(df, args.gene_col)
    score_col = detect_score_col(df, args.score_col)
    sgrna_col = detect_sgrna_col(df, args.sgrna_col)

    print(f"[INFO] Gene column: {gene_col}")
    print(f"[INFO] Score column: {score_col if score_col else 'not detected'}")
    print(f"[INFO] sgRNA column: {sgrna_col if sgrna_col else 'not detected'}")

    out = filter_guides(
        df=df,
        genes=genes,
        gene_col=gene_col,
        score_col=score_col,
        sgrna_col=sgrna_col,
        top_n_per_gene=args.top_n_per_gene,
        drop_duplicates=args.drop_duplicates,
    )

    if len(out) == 0:
        raise RuntimeError("No matching genes were found in the GuidePro reference table")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)
    write_summary(out, genes, gene_col, args.summary)

    print(f"[INFO] Requested genes: {len(genes)}")
    print(f"[INFO] Returned rows: {len(out)}")
    print(f"[INFO] Saved results to: {args.output}")
    if args.summary:
        print(f"[INFO] Saved summary to: {args.summary}")


if __name__ == "__main__":
    main()
