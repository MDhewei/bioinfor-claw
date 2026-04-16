#!/usr/bin/env python3

import argparse
import os
import re
import sys
import textwrap
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), '..', '..', '..', '_shared'))
from plot_style import init_style
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import MaxNLocator
import gseapy as gp


LIBRARY_MAP: Dict[str, Dict[str, str]] = {
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
    "KEGG": {
        "human": "KEGG_2021_Human",
        "mouse": "KEGG_2021_Mouse",
    },
    "REACTOME": {
        "human": "Reactome_2022",
        "mouse": "Reactome_2022",
    },
}


ORGANISM_MAP: Dict[str, str] = {
    "human": "human",
    "homo sapiens": "human",
    "hs": "human",
    "hsapiens": "human",
    "mouse": "mouse",
    "mus musculus": "mouse",
    "mm": "mouse",
}


DISPLAY_NAME_MAP: Dict[str, str] = {
    "GO_BP": "GO Biological Process",
    "GO_MF": "GO Molecular Function",
    "GO_CC": "GO Cellular Component",
    "KEGG": "KEGG Pathway",
    "REACTOME": "Reactome Pathway",
}


def normalize_organism(org: str) -> str:
    key = str(org).strip().lower()
    if key not in ORGANISM_MAP:
        raise ValueError(
            f"Unsupported organism '{org}'. Supported examples: human, mouse, hsapiens, mm"
        )
    return ORGANISM_MAP[key]


def ensure_parent_dir(path: str) -> None:
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)


def normalize_gene_symbol(gene: str) -> str:
    gene = str(gene).strip()
    gene = re.sub(r"\s+", "", gene)
    return gene.upper()


def read_gene_list(path: str) -> List[str]:
    genes: List[str] = []
    with open(path) as f:
        for line in f:
            x = line.strip()
            if not x:
                continue
            x = re.split(r"[\t,; ]+", x)[0]
            if x:
                genes.append(normalize_gene_symbol(x))

    genes = [g for g in genes if g]
    genes = list(dict.fromkeys(genes))

    if not genes:
        raise ValueError("Input gene list is empty after cleaning")

    return genes


def wrap_term(term: str, width: int = 45) -> str:
    return "\n".join(textwrap.wrap(str(term), width=width))


def parse_overlap_to_counts(overlap: str) -> Tuple[float, float]:
    try:
        a, b = str(overlap).split("/")
        return int(a), int(b)
    except Exception:
        return np.nan, np.nan


def clean_enrichr_results(df: pd.DataFrame) -> pd.DataFrame:
    if df is None or len(df) == 0:
        return pd.DataFrame()

    dat = df.copy()

    required = {"Term", "Adjusted P-value", "Overlap"}
    missing = required - set(dat.columns)
    if missing:
        raise ValueError(f"Missing expected enrichment columns: {sorted(missing)}")

    dat["Adjusted P-value"] = pd.to_numeric(dat["Adjusted P-value"], errors="coerce")
    dat = dat.dropna(subset=["Adjusted P-value"]).copy()

    hit_bg = dat["Overlap"].map(parse_overlap_to_counts)
    dat["Hit_Count"] = [x[0] for x in hit_bg]
    dat["Bg_Count"] = [x[1] for x in hit_bg]
    dat["Gene_Ratio"] = dat["Hit_Count"] / dat["Bg_Count"]
    dat["minus_log10_fdr"] = -np.log10(dat["Adjusted P-value"].clip(lower=np.nextafter(0, 1)))

    dat["Term_clean"] = dat["Term"].astype(str).str.replace(r"\s*\(GO:\d+\)\s*$", "", regex=True)
    dat["Term_wrapped"] = dat["Term_clean"].map(lambda x: wrap_term(x, width=45))

    dat = dat.sort_values(
        ["Adjusted P-value", "Hit_Count", "Gene_Ratio", "Term_clean"],
        ascending=[True, False, False, True],
    ).reset_index(drop=True)

    return dat


def resolve_library(category: str, organism: str) -> str:
    if category not in LIBRARY_MAP:
        raise ValueError(f"Unsupported category: {category}")
    if organism not in LIBRARY_MAP[category]:
        raise ValueError(f"No library configured for category={category}, organism={organism}")
    return LIBRARY_MAP[category][organism]


def run_enrichr(genes: List[str], gene_set: str, organism: str) -> pd.DataFrame:
    enr = gp.enrichr(
        gene_list=genes,
        gene_sets=gene_set,
        organism=organism,
        outdir=None,
        no_plot=True,
    )

    if enr.results is None or len(enr.results) == 0:
        return pd.DataFrame()

    return clean_enrichr_results(enr.results)


def choose_terms_for_plot(df: pd.DataFrame, top_n: int, sig_cutoff: float) -> pd.DataFrame:
    if len(df) == 0:
        return df

    sig = df[df["Adjusted P-value"] <= sig_cutoff].copy()
    if len(sig) == 0:
        sig = df.copy()

    return sig.head(top_n).copy()


def publication_style(font_family: str, font_size: float) -> None:
    plt.rcParams.update({
        "font.family": font_family,
        "font.size": font_size,
        "axes.labelsize": font_size + 1,
        "axes.titlesize": font_size + 2,
        "xtick.labelsize": font_size,
        "ytick.labelsize": font_size,
        "legend.fontsize": font_size - 0.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 1.0,
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
        "savefig.bbox": "tight",
    })


def make_bubble_plot(
    df: pd.DataFrame,
    output: str,
    title: str,
    top_n: int = 15,
    sig_cutoff: float = 0.05,
    fig_width: float = 8.5,
    fig_height: float = 6.5,
    dpi: int = 300,
    font_family: str = "Arial",
    font_size: float = 10.0,
    cmap: str = "viridis_r",
    min_bubble: float = 80.0,
    max_bubble: float = 650.0,
):
    if len(df) == 0:
        raise ValueError("No enrichment results available for plotting")

    dat = choose_terms_for_plot(df, top_n=top_n, sig_cutoff=sig_cutoff)
    if len(dat) == 0:
        raise ValueError("No enrichment terms available after filtering")

    publication_style(font_family=font_family, font_size=font_size)

    fig_height = max(fig_height, 0.38 * len(dat) + 1.8)
    dat = dat.iloc[::-1].copy()
    y_pos = np.arange(len(dat))

    hit_min = dat["Hit_Count"].min()
    hit_max = dat["Hit_Count"].max()
    if hit_min == hit_max:
        bubble_sizes = np.full(len(dat), (min_bubble + max_bubble) / 2)
    else:
        bubble_sizes = min_bubble + (
            (dat["Hit_Count"] - hit_min) / (hit_max - hit_min)
        ) * (max_bubble - min_bubble)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    norm = Normalize(
        vmin=float(dat["minus_log10_fdr"].min()),
        vmax=float(dat["minus_log10_fdr"].max()),
    )

    sc = ax.scatter(
        dat["Gene_Ratio"],
        y_pos,
        s=bubble_sizes,
        c=dat["minus_log10_fdr"],
        cmap=cmap,
        norm=norm,
        alpha=0.9,
        edgecolors="black",
        linewidths=0.4,
        zorder=3,
    )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(dat["Term_wrapped"])
    ax.set_xlabel("Gene Ratio")
    ax.set_title(title)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.grid(axis="x", linestyle="--", alpha=0.3, zorder=1)

    cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.88)
    cbar.set_label("-log10(FDR)")

    unique_sizes = sorted(dat["Hit_Count"].dropna().astype(int).unique())
    if len(unique_sizes) > 0:
        if len(unique_sizes) > 3:
            size_labels = [
                int(np.percentile(unique_sizes, 25)),
                int(np.percentile(unique_sizes, 50)),
                int(np.percentile(unique_sizes, 75)),
            ]
            size_labels = sorted(set(size_labels))
        else:
            size_labels = unique_sizes

        handles = []
        labels = []
        for s in size_labels:
            if hit_min == hit_max:
                area = (min_bubble + max_bubble) / 2
            else:
                area = min_bubble + ((s - hit_min) / (hit_max - hit_min)) * (max_bubble - min_bubble)
            h = ax.scatter([], [], s=area, facecolor="lightgray", edgecolor="black", linewidth=0.4)
            handles.append(h)
            labels.append(str(s))

        leg = ax.legend(
            handles,
            labels,
            title="Hit count",
            loc="lower right",
            frameon=False,
            borderpad=0.3,
            labelspacing=0.8,
            handletextpad=0.8,
        )
        ax.add_artist(leg)

    ensure_parent_dir(output)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi)
    plt.close(fig)


def make_summary_plot_table(df: pd.DataFrame, output: str) -> None:
    out = df.copy()
    cols = [
        "Term_clean",
        "Adjusted P-value",
        "minus_log10_fdr",
        "Overlap",
        "Hit_Count",
        "Bg_Count",
        "Gene_Ratio",
        "Combined Score",
        "Genes",
    ]
    existing = [c for c in cols if c in out.columns]
    out = out[existing].copy()
    ensure_parent_dir(output)
    out.to_csv(output, sep="\t", index=False)


def parse_panel(panel: str) -> List[str]:
    panel = panel.strip().upper()
    if panel == "ALL":
        return ["GO_BP", "GO_MF", "GO_CC", "KEGG", "REACTOME"]
    if panel == "GO":
        return ["GO_BP", "GO_MF", "GO_CC"]
    if panel in {"GO_BP", "GO_MF", "GO_CC", "KEGG", "REACTOME"}:
        return [panel]
    raise ValueError("Unsupported --library value")


def main():
    parser = argparse.ArgumentParser(
        description="GO, KEGG, and Reactome enrichment with publication-quality bubble plots"
    )
    parser.add_argument("--input", required=True, help="Input gene list file")
    parser.add_argument("--output-prefix", required=True, help="Output file prefix")
    parser.add_argument("--organism", default="human", help="Organism, e.g. human or mouse")
    parser.add_argument(
        "--library",
        default="ALL",
        choices=["GO_BP", "GO_MF", "GO_CC", "GO", "KEGG", "REACTOME", "ALL"],
    )
    parser.add_argument("--top-n", type=int, default=15)
    parser.add_argument("--sig-cutoff", type=float, default=0.05)
    parser.add_argument("--fig-width", type=float, default=8.5)
    parser.add_argument("--fig-height", type=float, default=6.5)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--font-family", default="Arial")
    parser.add_argument("--font-size", type=float, default=10.0)
    parser.add_argument("--cmap", default="viridis_r")
    args = parser.parse_args()

    init_style(
        font_family=getattr(args, 'font_family', None),
        font_size=getattr(args, 'font_size', None),
    )

    organism = normalize_organism(args.organism)
    genes = read_gene_list(args.input)
    categories = parse_panel(args.library)

    generated_any = False

    for category in categories:
        gene_set = resolve_library(category, organism)
        display_name = DISPLAY_NAME_MAP[category]
        print(f"[INFO] Running enrichment for {display_name}: {gene_set}")

        try:
            res = run_enrichr(genes, gene_set, organism=organism)
        except Exception as e:
            print(f"[ERROR] Failed enrichment for {category}: {e}", file=sys.stderr)
            continue

        out_tsv = f"{args.output_prefix}.{category}.tsv"
        out_plot_tsv = f"{args.output_prefix}.{category}.plot_input.tsv"
        out_png = f"{args.output_prefix}.{category}.bubble.png"
        out_pdf = f"{args.output_prefix}.{category}.bubble.pdf"

        if len(res) == 0:
            print(f"[WARN] No enrichment results for {category}")
            continue

        ensure_parent_dir(out_tsv)
        res.to_csv(out_tsv, sep="\t", index=False)

        plot_df = choose_terms_for_plot(res, top_n=args.top_n, sig_cutoff=args.sig_cutoff)
        if len(plot_df) == 0:
            print(f"[WARN] No plottable terms for {category}")
            continue

        make_summary_plot_table(plot_df, out_plot_tsv)

        title = f"{display_name} enrichment"

        try:
            make_bubble_plot(
                plot_df,
                output=out_png,
                title=title,
                top_n=args.top_n,
                sig_cutoff=args.sig_cutoff,
                fig_width=args.fig_width,
                fig_height=args.fig_height,
                dpi=args.dpi,
                font_family=args.font_family,
                font_size=args.font_size,
                cmap=args.cmap,
            )
            make_bubble_plot(
                plot_df,
                output=out_pdf,
                title=title,
                top_n=args.top_n,
                sig_cutoff=args.sig_cutoff,
                fig_width=args.fig_width,
                fig_height=args.fig_height,
                dpi=args.dpi,
                font_family=args.font_family,
                font_size=args.font_size,
                cmap=args.cmap,
            )
        except Exception as e:
            print(f"[ERROR] Failed plotting for {category}: {e}", file=sys.stderr)
            continue

        print(f"[INFO] Saved full result table: {out_tsv}")
        print(f"[INFO] Saved plot input table: {out_plot_tsv}")
        print(f"[INFO] Saved bubble plot PNG: {out_png}")
        print(f"[INFO] Saved bubble plot PDF: {out_pdf}")
        generated_any = True

    if not generated_any:
        raise RuntimeError("No enrichment result or plot was successfully generated.")


if __name__ == "__main__":
    main()