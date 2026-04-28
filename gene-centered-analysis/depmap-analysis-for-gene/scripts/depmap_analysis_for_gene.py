#!/usr/bin/env python3
"""
DepMap single-gene analysis: expression, mutation, copy number, essentiality.

This script streams DepMap data directly from the DepMap download API and
extracts ONLY the target gene's column, so it never needs to download or
load multi-GB whole-genome matrices into memory.

Co-expression and co-essentiality require the full matrix and are handled
by the standalone scripts depmap_coexpression.py / depmap_coessentiality.py.
"""

import argparse
import csv
import gzip
import io
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from urllib.request import urlopen, Request

import matplotlib.pyplot as plt
import sys as _sys, os as _os
try:
    _sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), *(['..'] * 3), '_shared'))
    from plot_style import init_style
except ImportError:
    def init_style(**kw): pass

import pandas as pd


# =========================================================
# DepMap API constants
# =========================================================
DEPMAP_RELEASE = "DepMap Public 26Q1"
DEPMAP_INDEX = "https://depmap.org/portal/api/download/files"

DATASET_FILES = {
    "expression": "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
    "copy_number": "PortalOmicsCNGeneLog2.csv",
    "mutations": "OmicsSomaticMutations.csv",
    "essentiality": "CRISPRGeneEffect.csv",
    "metadata": "Model.csv",
}

VALID_MODULES = {
    "expression", "mutation", "copy_number", "essentiality",
    "coexpression", "coessentiality", "full",
}


# =========================================================
# Utilities
# =========================================================
def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def normalize_gene_symbol(gene: str) -> str:
    return re.sub(r"\s+", "", gene.strip()).upper()


def parse_modules(modules_arg: Optional[str]) -> Set[str]:
    if not modules_arg:
        return {"full"}
    modules = {x.strip().lower() for x in re.split(r"[,\s;]+", modules_arg) if x.strip()}
    unknown = modules - VALID_MODULES
    if unknown:
        raise ValueError(f"Unknown module(s): {sorted(unknown)}. Valid: {sorted(VALID_MODULES)}")
    if "full" in modules:
        return {"expression", "mutation", "copy_number", "essentiality"}
    return modules


def required_datasets_for_modules(modules: Set[str]) -> Set[str]:
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
# Streaming DepMap data (gene-specific, no full download)
# =========================================================
_file_url_cache: Dict[str, str] = {}


def depmap_file_url(filename: str) -> str:
    """Resolve a DepMap filename to its download URL via the file index API."""
    if filename in _file_url_cache:
        return _file_url_cache[filename]

    index_data = urlopen(DEPMAP_INDEX, timeout=120).read().decode("utf-8")
    rows = list(csv.DictReader(io.StringIO(index_data)))
    for row in rows:
        if row.get("release") == DEPMAP_RELEASE and row.get("filename") == filename:
            _file_url_cache[filename] = row["url"]
            return row["url"]
    raise RuntimeError(f"Missing DepMap file: {filename} in release {DEPMAP_RELEASE}")


def stream_csv_url(url: str):
    """Stream a CSV from a URL row-by-row without loading the whole file."""
    with urlopen(url, timeout=300) as response:
        text = io.TextIOWrapper(response, encoding="utf-8", newline="")
        reader = csv.reader(text)
        for row in reader:
            yield row


def fetch_gene_column(filename: str, gene_symbol: str, cache_dir: Path) -> pd.DataFrame:
    """
    Stream a DepMap matrix CSV and extract only the target gene column.
    Returns a DataFrame with columns [ModelID, <gene_symbol>].
    Caches result locally for repeated calls in the same session.
    """
    cache = cache_dir / f"DepMap_{filename}_{gene_symbol}.csv"
    if cache.exists():
        return pd.read_csv(cache)

    print(f"[INFO] Streaming {filename} for {gene_symbol} (gene-only, not full download)...")
    url = depmap_file_url(filename)
    rows_iter = stream_csv_url(url)
    header = next(rows_iter)

    # Find the gene column (e.g. "TP53 (7157)" or "TP53")
    gene_col = None
    for i, h in enumerate(header):
        if h == gene_symbol or h.startswith(f"{gene_symbol} "):
            gene_col = i
            break
    if gene_col is None:
        raise ValueError(f"Gene {gene_symbol} not found in {filename}. First 10 columns: {header[:10]}")

    model_col = header.index("ModelID") if "ModelID" in header else 0
    default_col = header.index("IsDefaultEntryForModel") if "IsDefaultEntryForModel" in header else None

    rows = []
    for row in rows_iter:
        if default_col is not None and row[default_col] != "Yes":
            continue
        try:
            rows.append({"ModelID": row[model_col], gene_symbol: float(row[gene_col])})
        except (ValueError, IndexError):
            continue

    df = pd.DataFrame(rows)
    cache_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(cache, index=False)
    print(f"[INFO] Got {len(df)} models for {gene_symbol} from {filename}")
    return df


def fetch_gene_mutations(gene_symbol: str, cache_dir: Path) -> pd.DataFrame:
    """Stream the mutations CSV and extract only rows for the target gene."""
    cache = cache_dir / f"DepMap_{gene_symbol}_mutations.csv"
    if cache.exists():
        return pd.read_csv(cache)

    print(f"[INFO] Streaming OmicsSomaticMutations.csv for {gene_symbol}...")
    url = depmap_file_url(DATASET_FILES["mutations"])
    rows = []
    with urlopen(url, timeout=300) as response:
        text = io.TextIOWrapper(response, encoding="utf-8", newline="")
        for row in csv.DictReader(text):
            if row.get("HugoSymbol") != gene_symbol:
                continue
            if row.get("IsDefaultEntryForModel") not in ("Yes", "True", "true", "1"):
                continue
            rows.append({
                "ModelID": row.get("ModelID", ""),
                "Chrom": row.get("Chrom", ""),
                "Pos": row.get("Pos", ""),
                "ProteinChange": row.get("ProteinChange", ""),
                "Consequence": row.get("MolecularConsequence", row.get("VariantClassification", "")),
                "Impact": row.get("VepImpact", ""),
                "AF": row.get("AF", ""),
                "Hotspot": row.get("Hotspot", ""),
                "LikelyLoF": row.get("LikelyLoF", ""),
            })

    df = pd.DataFrame(rows)
    cache_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(cache, index=False)
    print(f"[INFO] Got {len(df)} mutation records for {gene_symbol}")
    return df


def fetch_model_metadata(cache_dir: Path) -> pd.DataFrame:
    """Fetch and cache DepMap model metadata."""
    cache = cache_dir / "DepMap_Model_metadata.csv"
    if cache.exists():
        return pd.read_csv(cache)

    print("[INFO] Fetching DepMap model metadata...")
    url = depmap_file_url(DATASET_FILES["metadata"])
    rows = []
    with urlopen(url, timeout=120) as response:
        text = io.TextIOWrapper(response, encoding="utf-8", newline="")
        for row in csv.DictReader(text):
            rows.append({
                "ModelID": row["ModelID"],
                "CellLineName": row.get("CellLineName", ""),
                "Lineage": row.get("OncotreeLineage") or "Unknown",
                "PrimaryDisease": row.get("OncotreePrimaryDisease") or "Unknown",
                "Subtype": row.get("OncotreeSubtype") or "Unknown",
            })
    df = pd.DataFrame(rows)
    cache_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(cache, index=False)
    return df


# =========================================================
# Data provider (streams gene-specific data)
# =========================================================
@dataclass
class DepMapGeneData:
    gene: str
    cache_dir: Path
    metadata: Optional[pd.DataFrame] = None
    _expr: Optional[pd.DataFrame] = None
    _cn: Optional[pd.DataFrame] = None
    _ess: Optional[pd.DataFrame] = None
    _mut: Optional[pd.DataFrame] = None

    def get_metadata(self) -> pd.DataFrame:
        if self.metadata is None:
            self.metadata = fetch_model_metadata(self.cache_dir)
        return self.metadata

    def get_expression(self) -> pd.Series:
        if self._expr is None:
            self._expr = fetch_gene_column(DATASET_FILES["expression"], self.gene, self.cache_dir)
        return self._expr.set_index("ModelID")[self.gene].rename("expression")

    def get_copy_number(self) -> pd.Series:
        if self._cn is None:
            self._cn = fetch_gene_column(DATASET_FILES["copy_number"], self.gene, self.cache_dir)
        return self._cn.set_index("ModelID")[self.gene].rename("copy_number")

    def get_essentiality(self) -> pd.Series:
        if self._ess is None:
            self._ess = fetch_gene_column(DATASET_FILES["essentiality"], self.gene, self.cache_dir)
        return self._ess.set_index("ModelID")[self.gene].rename("essentiality")

    def get_mutations(self) -> pd.DataFrame:
        if self._mut is None:
            self._mut = fetch_gene_mutations(self.gene, self.cache_dir)
        return self._mut

    def join_metadata(self, s: pd.Series) -> pd.DataFrame:
        meta = self.get_metadata().set_index("ModelID")
        return s.to_frame().join(meta, how="left")


# =========================================================
# Analysis helpers
# =========================================================
def summarize_vector(
    s: pd.Series, data: DepMapGeneData, value_name: str, top_n: int, ascending: bool = False
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    merged = data.join_metadata(s)
    merged["cell_line"] = merged["CellLineName"].fillna(merged.index.to_series())
    top_cells = (
        merged.sort_values(value_name, ascending=ascending)
        [["cell_line", value_name, "Lineage"]]
        .head(top_n)
        .reset_index()
        .rename(columns={"index": "depmap_id", "Lineage": "lineage"})
    )
    lineage_summary = (
        merged.groupby("Lineage", dropna=False)[value_name]
        .agg(["median", "mean", "count"])
        .reset_index()
        .rename(columns={"Lineage": "lineage"})
        .sort_values("median", ascending=ascending)
        .reset_index(drop=True)
    )
    return top_cells, lineage_summary


def classify_copy_number(value: float, amp_threshold: float = 1.0, del_threshold: float = -1.0) -> str:
    if pd.isna(value):
        return "NA"
    if value >= amp_threshold:
        return "amplified"
    if value <= del_threshold:
        return "deleted"
    return "neutral"


# =========================================================
# Plotting
# =========================================================
def make_horizontal_barplot(
    df: pd.DataFrame, label_col: str, value_col: str, out_png: str,
    title: str, xlabel: str, top_n: int = 20, ascending: bool = False,
    threshold: Optional[float] = None,
):
    if df.empty:
        return
    plot_df = df.head(top_n).copy().sort_values(value_col, ascending=ascending)
    fig_h = max(5, 0.35 * len(plot_df) + 1.5)
    plt.figure(figsize=(10, fig_h), dpi=300)
    ax = plt.gca()
    ax.barh(plot_df[label_col], plot_df[value_col], edgecolor="black", linewidth=0.4)
    if threshold is not None:
        ax.axvline(x=threshold, linestyle="--", linewidth=1.2, color="red", alpha=0.9,
                   label=f"Threshold = {threshold:g}")
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
    plt.savefig(os.path.splitext(out_png)[0] + ".pdf", bbox_inches="tight")
    plt.close()


# =========================================================
# Module runners
# =========================================================
def run_expression(gene: str, data: DepMapGeneData, args) -> dict:
    print(f"[INFO] Expression analysis for {gene}")
    expr_vec = data.get_expression()
    top_cells, lineage_summary = summarize_vector(
        expr_vec, data, value_name="expression", top_n=args.top_n, ascending=False
    )
    top_cells.to_csv(os.path.join(args.outdir, f"{gene}.expression_top_cell_lines.tsv"), sep="\t", index=False)
    lineage_summary.to_csv(os.path.join(args.outdir, f"{gene}.expression_lineage_summary.tsv"), sep="\t", index=False)
    make_horizontal_barplot(
        top_cells, label_col="cell_line", value_col="expression",
        out_png=os.path.join(args.outdir, f"{gene}.expression_barplot.png"),
        title=f"{gene} expression in DepMap cell lines", xlabel="Expression (TPM log2+1)",
        top_n=min(args.top_n, len(top_cells)), ascending=False,
    )
    frac = float((expr_vec > 1.0).sum() / len(expr_vec)) if len(expr_vec) > 0 else None
    return {
        "expression_models": len(expr_vec),
        "expression_max": float(expr_vec.max()) if len(expr_vec) else None,
        "expression_median": float(expr_vec.median()) if len(expr_vec) else None,
        "expression_detected_fraction_gt1": frac,
    }


def run_mutation(gene: str, data: DepMapGeneData, args) -> dict:
    print(f"[INFO] Mutation analysis for {gene}")
    mut_df = data.get_mutations()
    mut_df.to_csv(os.path.join(args.outdir, f"{gene}.mutations.tsv"), sep="\t", index=False)

    if not mut_df.empty:
        # Annotate with metadata
        meta = data.get_metadata().set_index("ModelID")
        mut_annotated = mut_df.merge(meta, left_on="ModelID", right_index=True, how="left")

        if "Consequence" in mut_df.columns:
            type_summary = mut_df["Consequence"].fillna("NA").value_counts().reset_index()
            type_summary.columns = ["mutation_type", "count"]
            type_summary.to_csv(os.path.join(args.outdir, f"{gene}.mutation_type_summary.tsv"), sep="\t", index=False)

        if "ProteinChange" in mut_df.columns:
            protein_summary = mut_df["ProteinChange"].fillna("NA").value_counts().reset_index()
            protein_summary.columns = ["protein_change", "count"]
            protein_summary.to_csv(os.path.join(args.outdir, f"{gene}.mutation_protein_summary.tsv"), sep="\t", index=False)

    return {
        "mutation_records": len(mut_df),
        "mutated_models": mut_df["ModelID"].nunique() if len(mut_df) else 0,
    }


def run_copy_number(gene: str, data: DepMapGeneData, args) -> dict:
    print(f"[INFO] Copy number analysis for {gene}")
    cn_vec = data.get_copy_number()
    cn_df = cn_vec.to_frame()
    cn_df["copy_number_status"] = cn_df["copy_number"].apply(
        lambda x: classify_copy_number(x, args.amp_threshold, args.del_threshold)
    )
    cn_df.to_csv(os.path.join(args.outdir, f"{gene}.copy_number.tsv"), sep="\t")

    top_cells, lineage_summary = summarize_vector(
        cn_vec, data, value_name="copy_number", top_n=args.top_n, ascending=False
    )

    # Top amplified
    top_amp = cn_df.nlargest(args.top_n, "copy_number").reset_index().rename(columns={"index": "depmap_id"})
    top_del = cn_df.nsmallest(args.top_n, "copy_number").reset_index().rename(columns={"index": "depmap_id"})
    top_amp.to_csv(os.path.join(args.outdir, f"{gene}.copy_number_top_amplified.tsv"), sep="\t", index=False)
    top_del.to_csv(os.path.join(args.outdir, f"{gene}.copy_number_top_deleted.tsv"), sep="\t", index=False)
    lineage_summary.to_csv(os.path.join(args.outdir, f"{gene}.copy_number_lineage_summary.tsv"), sep="\t", index=False)

    make_horizontal_barplot(
        top_cells, label_col="cell_line", value_col="copy_number",
        out_png=os.path.join(args.outdir, f"{gene}.copy_number_amplified_barplot.png"),
        title=f"{gene} top amplified cell lines", xlabel="Copy number (log2)",
        top_n=min(args.top_n, len(top_cells)), ascending=False, threshold=args.amp_threshold,
    )

    cn_amp_frac = float((cn_vec >= args.amp_threshold).sum() / len(cn_vec)) if len(cn_vec) > 0 else None
    cn_del_frac = float((cn_vec <= args.del_threshold).sum() / len(cn_vec)) if len(cn_vec) > 0 else None

    return {
        "copy_number_models": len(cn_vec),
        "copy_number_max": float(cn_vec.max()) if len(cn_vec) else None,
        "copy_number_min": float(cn_vec.min()) if len(cn_vec) else None,
        "copy_number_amplified_fraction": cn_amp_frac,
        "copy_number_deleted_fraction": cn_del_frac,
    }


def run_essentiality(gene: str, data: DepMapGeneData, args) -> dict:
    print(f"[INFO] Essentiality analysis for {gene}")
    dep_vec = data.get_essentiality()
    dep_vec.to_frame().to_csv(os.path.join(args.outdir, f"{gene}.essentiality.tsv"), sep="\t")

    top_cells, lineage_summary = summarize_vector(
        dep_vec, data, value_name="essentiality", top_n=args.top_n, ascending=True
    )
    top_cells.to_csv(os.path.join(args.outdir, f"{gene}.essentiality_top_cell_lines.tsv"), sep="\t", index=False)
    lineage_summary.to_csv(os.path.join(args.outdir, f"{gene}.essentiality_lineage_summary.tsv"), sep="\t", index=False)

    make_horizontal_barplot(
        top_cells, label_col="cell_line", value_col="essentiality",
        out_png=os.path.join(args.outdir, f"{gene}.essentiality_barplot.png"),
        title=f"{gene} most dependent cell lines", xlabel="Gene effect (CRISPR)",
        top_n=min(args.top_n, len(top_cells)), ascending=True, threshold=args.essential_threshold,
    )

    essential_frac = float((dep_vec < args.essential_threshold).sum() / len(dep_vec)) if len(dep_vec) > 0 else None
    return {
        "essentiality_models": len(dep_vec),
        "essentiality_min": float(dep_vec.min()) if len(dep_vec) else None,
        "essentiality_median": float(dep_vec.median()) if len(dep_vec) else None,
        "essential_fraction_lt_threshold": essential_frac,
    }


# =========================================================
# Main
# =========================================================
def main():
    parser = argparse.ArgumentParser(
        description="DepMap single-gene analysis: expression, mutation, copy number, essentiality. "
                    "Streams data directly from DepMap API — no full dataset download needed."
    )
    parser.add_argument("--gene", required=True, help="Gene symbol, e.g. TP53")
    parser.add_argument(
        "--modules", default="full",
        help="Comma-separated modules: expression,mutation,copy_number,essentiality,full. Default: full",
    )
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--top-n", type=int, default=20, help="Top rows to report")
    parser.add_argument("--amp-threshold", type=float, default=1.0)
    parser.add_argument("--del-threshold", type=float, default=-1.0)
    parser.add_argument("--essential-threshold", type=float, default=-0.5)

    # Legacy flags — accepted but IGNORED (data is now streamed directly)
    parser.add_argument("--expression-file", help=argparse.SUPPRESS)
    parser.add_argument("--mutation-file", help=argparse.SUPPRESS)
    parser.add_argument("--copy-number-file", help=argparse.SUPPRESS)
    parser.add_argument("--essentiality-file", help=argparse.SUPPRESS)
    parser.add_argument("--metadata-file", help=argparse.SUPPRESS)

    args = parser.parse_args()

    init_style(
        font_family=getattr(args, 'font_family', None),
        font_size=getattr(args, 'font_size', None),
    )

    modules = parse_modules(args.modules)
    gene = normalize_gene_symbol(args.gene)
    ensure_dir(args.outdir)

    if "coexpression" in modules or "coessentiality" in modules:
        print("[ERROR] Co-expression and co-essentiality require the full matrix.")
        print("[ERROR] Use the standalone scripts: depmap_coexpression.py or depmap_coessentiality.py")
        raise SystemExit(1)

    # Use a cache directory inside outdir for gene-specific extracts
    cache_dir = Path(args.outdir) / ".depmap_cache"
    data = DepMapGeneData(gene=gene, cache_dir=cache_dir)

    summary_data = {"query_gene": gene, "modules_run": ",".join(sorted(modules))}

    if "expression" in modules:
        summary_data.update(run_expression(gene, data, args))

    if "mutation" in modules:
        summary_data.update(run_mutation(gene, data, args))

    if "copy_number" in modules:
        summary_data.update(run_copy_number(gene, data, args))

    if "essentiality" in modules:
        summary_data.update(run_essentiality(gene, data, args))

    summary = pd.DataFrame([summary_data])
    summary.to_csv(os.path.join(args.outdir, "gene_depmap_summary.tsv"), sep="\t", index=False)

    with open(os.path.join(args.outdir, "summary.txt"), "w", encoding="utf-8") as f:
        for key, value in summary_data.items():
            f.write(f"{key}: {value}\n")

    # Clean up cache
    import shutil
    if cache_dir.exists():
        shutil.rmtree(cache_dir, ignore_errors=True)

    # Print key results to stdout for agent consumption
    print(f"\n[RESULTS] === DepMap Analysis: {gene} | Modules: {', '.join(sorted(modules))} ===")
    for key, value in summary_data.items():
        if key in ('query_gene', 'modules_run'):
            continue
        if value is not None and value != '' and str(value) != 'nan':
            print(f"[RESULTS] {key}: {value}")
    print(f"[RESULTS] === END ===")
    print(f"[DONE] Results written to: {args.outdir}")


if __name__ == "__main__":
    main()
