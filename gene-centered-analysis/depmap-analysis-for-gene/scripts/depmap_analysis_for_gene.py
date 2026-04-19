#!/usr/bin/env python3

import argparse
import os
import re
from dataclasses import dataclass
from typing import List, Optional, Set, Tuple

import matplotlib.pyplot as plt
import sys as _sys, os as _os
try:
    _sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), *(['..'] * 3), '_shared'))
    from plot_style import init_style
except ImportError:
    def init_style(**kw): pass  # graceful fallback if _shared not available
import pandas as pd


# =========================================================
# Utilities
# =========================================================
VALID_MODULES = {
    "expression",
    "mutation",
    "copy_number",
    "essentiality",
    "coexpression",
    "coessentiality",
    "full",
}


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def safe_read_csv(path: str, nrows: Optional[int] = None) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False, nrows=nrows)


def first_existing_column(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def normalize_gene_symbol(gene: str) -> str:
    return re.sub(r"\s+", "", gene.strip()).upper()


def find_gene_column(columns: List[str], gene_symbol: str) -> Optional[str]:
    """
    Match common DepMap matrix column styles:
      TP53
      TP53 (7157)
      TP53 [ENSG...]
      TP53 something...
    """
    if gene_symbol in columns:
        return gene_symbol

    pat = re.compile(rf"^{re.escape(gene_symbol)}(\s|\(|\[|$)", re.IGNORECASE)
    matches = [c for c in columns if pat.search(str(c))]
    if matches:
        return matches[0]

    return None


def infer_depmap_id_column(df: pd.DataFrame, context: str) -> str:
    """
    Try standard names first, then fall back to the first column if it looks like ACH-xxxxx.
    """
    id_col = first_existing_column(df, ["DepMap_ID", "ModelID", "depmap_id"])
    if id_col is not None:
        return id_col

    if len(df.columns) == 0:
        raise ValueError(f"{context}: file has no columns")

    first_col = df.columns[0]
    vals = df[first_col].astype(str)
    ach_like = vals.str.match(r"^ACH-\d+").sum()

    if first_col == "Unnamed: 0" or ach_like > 0:
        return first_col

    raise ValueError(
        f"{context}: could not find DepMap ID column. "
        f"First 20 columns: {list(df.columns[:20])}"
    )


def extract_symbol_from_gene_column(col: str) -> str:
    m = re.match(r"^([A-Za-z0-9_.-]+)", str(col))
    return m.group(1).upper() if m else str(col).upper()


def parse_modules(modules_arg: Optional[str]) -> Set[str]:
    if not modules_arg:
        return {"full"}

    modules = {
        x.strip().lower()
        for x in re.split(r"[,\s;]+", modules_arg)
        if x.strip()
    }

    unknown = modules - VALID_MODULES
    if unknown:
        raise ValueError(
            f"Unknown module(s): {sorted(unknown)}. "
            f"Valid modules: {sorted(VALID_MODULES)}"
        )

    if "full" in modules:
        return {"expression", "mutation", "copy_number", "essentiality", "coexpression", "coessentiality"}

    return modules


def required_inputs_for_modules(modules: Set[str]) -> Set[str]:
    req = set()
    if "expression" in modules or "coexpression" in modules:
        req.add("expression")
    if "mutation" in modules:
        req.add("mutations")
    if "copy_number" in modules:
        req.add("copy_number")
    if "essentiality" in modules or "coessentiality" in modules:
        req.add("essentiality")
    if modules:
        req.add("metadata")
    return req


# =========================================================
# Data classes
# =========================================================
@dataclass
class DepMapPaths:
    expression: Optional[str] = None
    mutations: Optional[str] = None
    copy_number: Optional[str] = None
    essentiality: Optional[str] = None
    metadata: Optional[str] = None


# =========================================================
# Provider
# =========================================================
class DepMapProvider:
    def __init__(self, paths: DepMapPaths):
        self.paths = paths
        self.expr_df: Optional[pd.DataFrame] = None
        self.cn_df: Optional[pd.DataFrame] = None
        self.dep_df: Optional[pd.DataFrame] = None
        self.mut_df: Optional[pd.DataFrame] = None
        self.meta_df: Optional[pd.DataFrame] = None

    def _load_matrix(self, path: str, context: str) -> pd.DataFrame:
        df = safe_read_csv(path)
        id_col = infer_depmap_id_column(df, context)
        df = df.set_index(id_col)

        num_df = df.apply(pd.to_numeric, errors="coerce")
        keep_cols = num_df.notna().any(axis=0)
        num_df = num_df.loc[:, keep_cols]

        if num_df.shape[1] == 0:
            raise ValueError(
                f"{context}: no numeric gene columns found after parsing. "
                f"Original columns (first 20): {list(df.columns[:20])}"
            )

        return num_df

    def load_expression(self) -> pd.DataFrame:
        if self.paths.expression is None:
            raise ValueError("Expression file path not provided.")
        if self.expr_df is None:
            self.expr_df = self._load_matrix(self.paths.expression, "expression file")
        return self.expr_df

    def load_copy_number(self) -> pd.DataFrame:
        if self.paths.copy_number is None:
            raise ValueError("Copy number file path not provided.")
        if self.cn_df is None:
            self.cn_df = self._load_matrix(self.paths.copy_number, "copy number file")
        return self.cn_df

    def load_essentiality(self) -> pd.DataFrame:
        if self.paths.essentiality is None:
            raise ValueError("Essentiality file path not provided.")
        if self.dep_df is None:
            self.dep_df = self._load_matrix(self.paths.essentiality, "essentiality file")
        return self.dep_df

    def load_mutations(self) -> pd.DataFrame:
        if self.paths.mutations is None:
            raise ValueError("Mutation file path not provided.")
        if self.mut_df is None:
            self.mut_df = safe_read_csv(self.paths.mutations)
        return self.mut_df

    def load_metadata(self) -> pd.DataFrame:
        if self.paths.metadata is None:
            raise ValueError("Metadata file path not provided.")
        if self.meta_df is None:
            df = safe_read_csv(self.paths.metadata)
            id_col = infer_depmap_id_column(df, "metadata file")
            self.meta_df = df.set_index(id_col)
        return self.meta_df

    def join_metadata(self, df: pd.DataFrame) -> pd.DataFrame:
        meta = self.load_metadata()
        return df.join(meta, how="left")

    def get_lineage_column(self, df: pd.DataFrame) -> Optional[str]:
        return first_existing_column(
            df,
            ["OncotreeLineage", "lineage", "Lineage", "primary_disease", "PrimaryDisease"],
        )

    def get_cellline_name_column(self, df: pd.DataFrame) -> Optional[str]:
        return first_existing_column(
            df,
            ["stripped_cell_line_name", "CellLineName", "cell_line_name", "CCLE_Name", "ModelName"],
        )

    def get_expression_vector(self, gene: str) -> pd.Series:
        df = self.load_expression()
        col = find_gene_column(list(df.columns), gene)
        if col is None:
            raise ValueError(f"Gene {gene} not found in expression matrix.")
        return df[col].rename("expression")

    def get_copy_number_vector(self, gene: str) -> pd.Series:
        df = self.load_copy_number()
        col = find_gene_column(list(df.columns), gene)
        if col is None:
            raise ValueError(f"Gene {gene} not found in copy number matrix.")
        return df[col].rename("copy_number")

    def get_essentiality_vector(self, gene: str) -> pd.Series:
        df = self.load_essentiality()
        col = find_gene_column(list(df.columns), gene)
        if col is None:
            raise ValueError(f"Gene {gene} not found in essentiality matrix.")
        return df[col].rename("essentiality")

    def get_mutation_rows(self, gene: str) -> pd.DataFrame:
        df = self.load_mutations()
        gene_col = first_existing_column(df, ["HugoSymbol", "gene", "Gene", "gene_symbol"])
        if gene_col is None:
            raise ValueError("Could not find gene symbol column in mutation file.")
        return df[df[gene_col].astype(str).str.upper() == gene.upper()].copy()


# =========================================================
# Analysis functions
# =========================================================
def summarize_vector_by_cellline_and_lineage(
    s: pd.Series,
    provider: DepMapProvider,
    value_name: str,
    top_n: int,
    ascending: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    df = s.to_frame()
    merged = provider.join_metadata(df)

    lineage_col = provider.get_lineage_column(merged)
    cellline_col = provider.get_cellline_name_column(merged)

    merged["cell_line"] = merged[cellline_col] if cellline_col else merged.index

    top_cells = merged.sort_values(value_name, ascending=ascending).copy()
    keep_cols = ["cell_line", value_name]
    if lineage_col is not None:
        top_cells["lineage"] = top_cells[lineage_col]
        keep_cols.append("lineage")

    top_cells = (
        top_cells[keep_cols]
        .head(top_n)
        .reset_index()
        .rename(columns={"index": "depmap_id"})
    )

    if lineage_col is not None:
        lineage_summary = (
            merged.groupby(lineage_col, dropna=False)[value_name]
            .agg(["median", "mean", "count"])
            .reset_index()
            .rename(columns={lineage_col: "lineage"})
            .sort_values("median", ascending=ascending)
            .reset_index(drop=True)
        )
    else:
        lineage_summary = pd.DataFrame(columns=["lineage", "median", "mean", "count"])

    return top_cells, lineage_summary


def summarize_mutations(
    mut_df: pd.DataFrame,
    provider: DepMapProvider,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if mut_df.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    meta = provider.load_metadata()
    depmap_col = first_existing_column(mut_df, ["DepMap_ID", "ModelID", "depmap_id"])
    protein_col = first_existing_column(mut_df, ["ProteinChange", "protein_change", "HGVSp_Short", "HGVSp"])
    variant_col = first_existing_column(mut_df, ["VariantInfo", "variant_info", "VariantClassification", "variant_classification"])

    out = mut_df.copy()
    if depmap_col and depmap_col in out.columns:
        out = out.merge(meta, left_on=depmap_col, right_index=True, how="left")

    if variant_col:
        mut_type_summary = out[variant_col].fillna("NA").value_counts().reset_index()
        mut_type_summary.columns = ["mutation_type", "count"]
    else:
        mut_type_summary = pd.DataFrame(columns=["mutation_type", "count"])

    if protein_col:
        protein_summary = out[protein_col].fillna("NA").value_counts().reset_index()
        protein_summary.columns = ["protein_change", "count"]
    else:
        protein_summary = pd.DataFrame(columns=["protein_change", "count"])

    return out, mut_type_summary, protein_summary


def classify_copy_number(value: float, amp_threshold: float = 1.0, del_threshold: float = -1.0) -> str:
    if pd.isna(value):
        return "NA"
    if value >= amp_threshold:
        return "amplified"
    if value <= del_threshold:
        return "deleted"
    return "neutral"


def summarize_copy_number(
    s: pd.Series,
    provider: DepMapProvider,
    top_n: int,
    amp_threshold: float = 1.0,
    del_threshold: float = -1.0,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    df = s.to_frame()
    df["copy_number_status"] = df["copy_number"].apply(
        lambda x: classify_copy_number(x, amp_threshold=amp_threshold, del_threshold=del_threshold)
    )

    merged = provider.join_metadata(df)
    lineage_col = provider.get_lineage_column(merged)
    cellline_col = provider.get_cellline_name_column(merged)
    merged["cell_line"] = merged[cellline_col] if cellline_col else merged.index

    top_amp = (
        merged.sort_values("copy_number", ascending=False)
        [["cell_line", "copy_number", "copy_number_status"]]
        .head(top_n)
        .reset_index()
        .rename(columns={"index": "depmap_id"})
    )
    top_del = (
        merged.sort_values("copy_number", ascending=True)
        [["cell_line", "copy_number", "copy_number_status"]]
        .head(top_n)
        .reset_index()
        .rename(columns={"index": "depmap_id"})
    )

    if lineage_col:
        amp_lineages = merged.sort_values("copy_number", ascending=False)[lineage_col].head(top_n).values
        del_lineages = merged.sort_values("copy_number", ascending=True)[lineage_col].head(top_n).values
        top_amp["lineage"] = amp_lineages
        top_del["lineage"] = del_lineages

        lineage_summary = (
            merged.groupby(lineage_col, dropna=False)["copy_number"]
            .agg(["median", "mean", "count"])
            .reset_index()
            .rename(columns={lineage_col: "lineage"})
            .sort_values("median", ascending=False)
            .reset_index(drop=True)
        )
    else:
        lineage_summary = pd.DataFrame(columns=["lineage", "median", "mean", "count"])

    return top_amp, top_del, lineage_summary


def compute_correlations_to_all_genes(
    matrix_df: pd.DataFrame,
    gene: str,
    method: str = "pearson",
    top_n: int = 20,
) -> pd.DataFrame:
    col = find_gene_column(list(matrix_df.columns), gene)
    if col is None:
        raise ValueError(f"Gene {gene} not found in matrix for correlation.")

    target = matrix_df[col]
    corr = matrix_df.corrwith(target, method=method).dropna()

    out = corr.reset_index()
    out.columns = ["gene_column", "correlation"]
    out["gene_symbol"] = out["gene_column"].apply(extract_symbol_from_gene_column)
    out = out[out["gene_symbol"] != gene.upper()].copy()
    out = out.sort_values("correlation", ascending=False).reset_index(drop=True)
    return out.head(top_n)


# =========================================================
# Plotting
# =========================================================
def make_horizontal_barplot(
    df: pd.DataFrame,
    label_col: str,
    value_col: str,
    out_png: str,
    title: str,
    xlabel: str,
    top_n: int = 20,
    ascending: bool = False,
    threshold: Optional[float] = None,
):
    if df.empty:
        return

    plot_df = df.head(top_n).copy().sort_values(value_col, ascending=ascending)

    fig_w = 10
    fig_h = max(5, 0.35 * len(plot_df) + 1.5)

    plt.figure(figsize=(fig_w, fig_h), dpi=300)
    ax = plt.gca()

    ax.barh(plot_df[label_col], plot_df[value_col], edgecolor="black", linewidth=0.4)

    if threshold is not None:
        ax.axvline(
            x=threshold,
            linestyle="--",
            linewidth=1.2,
            color="red",
            alpha=0.9,
            label=f"Threshold = {threshold:g}",
        )
        ax.legend(frameon=False, fontsize=10, loc="best")

    ax.set_title(title, fontsize=14, pad=10)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel("Cell line", fontsize=12)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", linestyle=":", linewidth=0.5, alpha=0.5)
    ax.set_axisbelow(True)

    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    pdf_path = os.path.splitext(out_png)[0] + ".pdf"
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()


# =========================================================
# Module runners
# =========================================================
def run_expression(gene: str, provider: DepMapProvider, args) -> dict:
    print(f"[INFO] Expression analysis for {gene}")
    expr_vec = provider.get_expression_vector(gene)
    expr_top_cells, expr_lineage_summary = summarize_vector_by_cellline_and_lineage(
        expr_vec, provider, value_name="expression", top_n=args.top_n, ascending=False
    )

    expr_top_cells.to_csv(os.path.join(args.outdir, f"{gene}.expression_top_cell_lines.tsv"), sep="\t", index=False)
    expr_lineage_summary.to_csv(os.path.join(args.outdir, f"{gene}.expression_lineage_summary.tsv"), sep="\t", index=False)

    make_horizontal_barplot(
        expr_top_cells,
        label_col="cell_line",
        value_col="expression",
        out_png=os.path.join(args.outdir, f"{gene}.expression_barplot.png"),
        title=f"{gene} expression in DepMap cell lines",
        xlabel="Expression",
        top_n=min(args.top_n, len(expr_top_cells)),
        ascending=False,
    )

    expr_detected_frac = float((expr_vec > 1.0).sum() / len(expr_vec)) if len(expr_vec) > 0 else None
    return {
        "expression_max": float(expr_vec.max()) if len(expr_vec) else None,
        "expression_median": float(expr_vec.median()) if len(expr_vec) else None,
        "expression_detected_fraction_gt1": expr_detected_frac,
    }


def run_mutation(gene: str, provider: DepMapProvider, args) -> dict:
    print(f"[INFO] Mutation analysis for {gene}")
    mut_rows = provider.get_mutation_rows(gene)
    mut_full, mut_type_summary, mut_protein_summary = summarize_mutations(mut_rows, provider)

    mut_full.to_csv(os.path.join(args.outdir, f"{gene}.mutations.tsv"), sep="\t", index=False)
    mut_type_summary.to_csv(os.path.join(args.outdir, f"{gene}.mutation_type_summary.tsv"), sep="\t", index=False)
    mut_protein_summary.to_csv(os.path.join(args.outdir, f"{gene}.mutation_protein_summary.tsv"), sep="\t", index=False)

    return {
        "mutation_rows": int(len(mut_rows)),
    }


def run_copy_number(gene: str, provider: DepMapProvider, args) -> dict:
    print(f"[INFO] Copy number analysis for {gene}")
    cn_vec = provider.get_copy_number_vector(gene)
    cn_df = cn_vec.to_frame()
    cn_df["copy_number_status"] = cn_df["copy_number"].apply(
        lambda x: classify_copy_number(x, amp_threshold=args.amp_threshold, del_threshold=args.del_threshold)
    )
    cn_df.to_csv(os.path.join(args.outdir, f"{gene}.copy_number.tsv"), sep="\t")

    cn_top_amp, cn_top_del, cn_lineage_summary = summarize_copy_number(
        cn_vec, provider, top_n=args.top_n, amp_threshold=args.amp_threshold, del_threshold=args.del_threshold
    )
    cn_top_amp.to_csv(os.path.join(args.outdir, f"{gene}.copy_number_top_amplified.tsv"), sep="\t", index=False)
    cn_top_del.to_csv(os.path.join(args.outdir, f"{gene}.copy_number_top_deleted.tsv"), sep="\t", index=False)
    cn_lineage_summary.to_csv(os.path.join(args.outdir, f"{gene}.copy_number_lineage_summary.tsv"), sep="\t", index=False)

    make_horizontal_barplot(
        cn_top_amp,
        label_col="cell_line",
        value_col="copy_number",
        out_png=os.path.join(args.outdir, f"{gene}.copy_number_amplified_barplot.png"),
        title=f"{gene} top amplified cell lines",
        xlabel="Copy number",
        top_n=min(args.top_n, len(cn_top_amp)),
        ascending=False,
        threshold=args.amp_threshold,
    )
    make_horizontal_barplot(
        cn_top_del,
        label_col="cell_line",
        value_col="copy_number",
        out_png=os.path.join(args.outdir, f"{gene}.copy_number_deleted_barplot.png"),
        title=f"{gene} top deleted cell lines",
        xlabel="Copy number",
        top_n=min(args.top_n, len(cn_top_del)),
        ascending=True,
        threshold=args.del_threshold,
    )

    cn_amp_frac = float((cn_vec >= args.amp_threshold).sum() / len(cn_vec)) if len(cn_vec) > 0 else None
    cn_del_frac = float((cn_vec <= args.del_threshold).sum() / len(cn_vec)) if len(cn_vec) > 0 else None

    return {
        "copy_number_max": float(cn_vec.max()) if len(cn_vec) else None,
        "copy_number_min": float(cn_vec.min()) if len(cn_vec) else None,
        "copy_number_amplified_fraction": cn_amp_frac,
        "copy_number_deleted_fraction": cn_del_frac,
    }


def run_essentiality(gene: str, provider: DepMapProvider, args) -> dict:
    print(f"[INFO] Essentiality analysis for {gene}")
    dep_vec = provider.get_essentiality_vector(gene)
    dep_df = dep_vec.to_frame()
    dep_df.to_csv(os.path.join(args.outdir, f"{gene}.essentiality.tsv"), sep="\t")

    dep_top_cells, dep_lineage_summary = summarize_vector_by_cellline_and_lineage(
        dep_vec, provider, value_name="essentiality", top_n=args.top_n, ascending=True
    )
    dep_top_cells.to_csv(os.path.join(args.outdir, f"{gene}.essentiality_top_cell_lines.tsv"), sep="\t", index=False)
    dep_lineage_summary.to_csv(os.path.join(args.outdir, f"{gene}.essentiality_lineage_summary.tsv"), sep="\t", index=False)

    make_horizontal_barplot(
        dep_top_cells,
        label_col="cell_line",
        value_col="essentiality",
        out_png=os.path.join(args.outdir, f"{gene}.essentiality_barplot.png"),
        title=f"{gene} most dependent cell lines",
        xlabel="Gene effect / dependency",
        top_n=min(args.top_n, len(dep_top_cells)),
        ascending=True,
        threshold=args.essential_threshold,
    )

    essential_frac = float((dep_vec < args.essential_threshold).sum() / len(dep_vec)) if len(dep_vec) > 0 else None
    return {
        "essentiality_min": float(dep_vec.min()) if len(dep_vec) else None,
        "essentiality_median": float(dep_vec.median()) if len(dep_vec) else None,
        "essential_fraction_lt_threshold": essential_frac,
    }


def run_coexpression(gene: str, provider: DepMapProvider, args) -> dict:
    print(f"[INFO] Co-expression analysis for {gene}")
    expr_matrix = provider.load_expression()
    coexpr_df = compute_correlations_to_all_genes(
        expr_matrix, gene=gene, method=args.corr_method, top_n=args.top_n
    )
    coexpr_df.to_csv(os.path.join(args.outdir, f"{gene}.coexpression.tsv"), sep="\t", index=False)

    return {
        "top_coexpressed_gene": coexpr_df.iloc[0]["gene_symbol"] if len(coexpr_df) else None,
        "top_coexpressed_corr": float(coexpr_df.iloc[0]["correlation"]) if len(coexpr_df) else None,
    }


def run_coessentiality(gene: str, provider: DepMapProvider, args) -> dict:
    print(f"[INFO] Co-essentiality analysis for {gene}")
    dep_matrix = provider.load_essentiality()
    coess_df = compute_correlations_to_all_genes(
        dep_matrix, gene=gene, method=args.corr_method, top_n=args.top_n
    )
    coess_df.to_csv(os.path.join(args.outdir, f"{gene}.coessentiality.tsv"), sep="\t", index=False)

    return {
        "top_coessential_gene": coess_df.iloc[0]["gene_symbol"] if len(coess_df) else None,
        "top_coessential_corr": float(coess_df.iloc[0]["correlation"]) if len(coess_df) else None,
    }


# =========================================================
# Main
# =========================================================
def main():
    parser = argparse.ArgumentParser(
        description="DepMap single-gene analysis: expression, mutation, copy number, essentiality, co-expression, co-essentiality."
    )
    parser.add_argument("--gene", required=True, help="Gene symbol, e.g. TP53")
    parser.add_argument("--expression-file", help="DepMap expression CSV")
    parser.add_argument("--mutation-file", help="DepMap mutations CSV")
    parser.add_argument("--copy-number-file", help="DepMap copy number CSV")
    parser.add_argument("--essentiality-file", help="DepMap CRISPR gene effect CSV")
    parser.add_argument("--metadata-file", help="DepMap model metadata CSV")
    parser.add_argument(
        "--modules",
        default="full",
        help="Comma-separated modules to run. Choices: expression,mutation,copy_number,essentiality,coexpression,coessentiality,full. Default: full",
    )
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--top-n", type=int, default=20, help="Top rows to report")
    parser.add_argument("--corr-method", choices=["pearson", "spearman"], default="pearson")
    parser.add_argument("--amp-threshold", type=float, default=1.0)
    parser.add_argument("--del-threshold", type=float, default=-1.0)
    parser.add_argument("--essential-threshold", type=float, default=-0.5)
    args = parser.parse_args()

    init_style(
        font_family=getattr(args, 'font_family', None),
        font_size=getattr(args, 'font_size', None),
    )

    modules = parse_modules(args.modules)
    gene = normalize_gene_symbol(args.gene)
    ensure_dir(args.outdir)

    required = required_inputs_for_modules(modules)

    missing_cli = []
    if "expression" in required and not args.expression_file:
        missing_cli.append("--expression-file")
    if "mutations" in required and not args.mutation_file:
        missing_cli.append("--mutation-file")
    if "copy_number" in required and not args.copy_number_file:
        missing_cli.append("--copy-number-file")
    if "essentiality" in required and not args.essentiality_file:
        missing_cli.append("--essentiality-file")
    if "metadata" in required and not args.metadata_file:
        missing_cli.append("--metadata-file")

    if missing_cli:
        raise ValueError(
            "Missing required argument(s) for requested module(s): "
            + ", ".join(missing_cli)
        )

    provider = DepMapProvider(
        DepMapPaths(
            expression=args.expression_file,
            mutations=args.mutation_file,
            copy_number=args.copy_number_file,
            essentiality=args.essentiality_file,
            metadata=args.metadata_file,
        )
    )

    summary_data = {
        "query_gene": gene,
        "modules_run": ",".join(sorted(modules)),
    }

    if "expression" in modules:
        summary_data.update(run_expression(gene, provider, args))

    if "mutation" in modules:
        summary_data.update(run_mutation(gene, provider, args))

    if "copy_number" in modules:
        summary_data.update(run_copy_number(gene, provider, args))

    if "essentiality" in modules:
        summary_data.update(run_essentiality(gene, provider, args))

    if "coexpression" in modules:
        summary_data.update(run_coexpression(gene, provider, args))

    if "coessentiality" in modules:
        summary_data.update(run_coessentiality(gene, provider, args))

    summary = pd.DataFrame([summary_data])
    summary.to_csv(os.path.join(args.outdir, "gene_depmap_summary.tsv"), sep="\t", index=False)

    with open(os.path.join(args.outdir, "summary.txt"), "w", encoding="utf-8") as f:
        for key, value in summary_data.items():
            f.write(f"{key}: {value}\n")

    print(f"[DONE] Results written to: {args.outdir}")


if __name__ == "__main__":
    main()