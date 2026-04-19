#!/usr/bin/env python3

import argparse
import csv
import json
import os
import re
import sys
import textwrap
import urllib.request
import urllib.error
import urllib.parse
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys as _sys, os as _os
try:
    _sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), *(['..'] * 3), '_shared'))
    from plot_style import init_style
except ImportError:
    def init_style(**kw): pass  # graceful fallback if _shared not available
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import MaxNLocator

try:
    import gseapy as gp
except ImportError:
    gp = None  # g:Profiler is used as primary; gseapy/Enrichr as fallback


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

# g:Profiler source names (maps our category names → g:Profiler source IDs)
GPROFILER_SOURCE_MAP: Dict[str, str] = {
    "GO_BP": "GO:BP",
    "GO_MF": "GO:MF",
    "GO_CC": "GO:CC",
    "KEGG": "KEGG",
    "REACTOME": "REAC",
}

GPROFILER_ORGANISM_MAP: Dict[str, str] = {
    "human": "hsapiens",
    "mouse": "mmusculus",
}

GPROFILER_API = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"


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

    # Drop Enrichr-specific columns that are not meaningful for g:Profiler
    drop_cols = {"Old P-value", "Old Adjusted P-value", "Odds Ratio", "Combined Score"}
    dat = dat.drop(columns=[c for c in drop_cols if c in dat.columns])

    return dat


def resolve_library(category: str, organism: str) -> str:
    if category not in LIBRARY_MAP:
        raise ValueError(f"Unsupported category: {category}")
    if organism not in LIBRARY_MAP[category]:
        raise ValueError(f"No library configured for category={category}, organism={organism}")
    return LIBRARY_MAP[category][organism]


def _gprofiler_df_to_enrichr_format(result_df: pd.DataFrame,
                                      gp_source: str) -> pd.DataFrame:
    """Convert a gprofiler-official result DataFrame to our enrichr-like format."""
    rows = []
    for _, r in result_df.iterrows():
        term_name = str(r.get("name", ""))
        term_id = str(r.get("native", ""))
        p_value = float(r.get("p_value", 1.0))
        intersection_size = int(r.get("intersection_size", 0))
        term_size = int(r.get("term_size", 1))

        # Build display term (include GO ID for GO terms)
        if term_id.startswith("GO:"):
            display = f"{term_name} ({term_id})"
        else:
            display = term_name

        # Get intersection genes if available
        intersections = r.get("intersections", [])
        if isinstance(intersections, list):
            # Filter out empty strings / None
            gene_list = [str(g) for g in intersections if g and str(g).strip()]
            gene_str = ";".join(gene_list)
        else:
            gene_str = str(intersections) if pd.notna(intersections) else ""

        rows.append({
            "Gene_set": f"g:Profiler_{gp_source}",
            "Term": display,
            "Overlap": f"{intersection_size}/{term_size}",
            "P-value": p_value,
            "Adjusted P-value": p_value,  # g:Profiler returns already-corrected p-values
            "Genes": gene_str,
        })

    if not rows:
        return pd.DataFrame()
    return clean_enrichr_results(pd.DataFrame(rows))


def run_gprofiler(genes: List[str], category: str, organism: str,
                   sig_cutoff: float = 0.05) -> pd.DataFrame:
    """
    Run enrichment via g:Profiler.
    Uses official GO consortium annotations with automatic gene symbol mapping.
    Tries gprofiler-official package first, then raw REST API fallback.
    Returns DataFrame in the same format as clean_enrichr_results().
    """
    gp_organism = GPROFILER_ORGANISM_MAP.get(organism)
    gp_source = GPROFILER_SOURCE_MAP.get(category)
    if not gp_organism or not gp_source:
        raise ValueError(f"g:Profiler does not support organism={organism} or category={category}")

    # ── Method 1: gprofiler-official package (most reliable) ──────────
    try:
        from gprofiler import GProfiler
        print(f"[INFO] Using gprofiler-official package")
        gp_inst = GProfiler(return_dataframe=True)
        result_df = gp_inst.profile(
            organism=gp_organism,
            query=genes,
            sources=[gp_source],
            user_threshold=sig_cutoff,
            significance_threshold_method="fdr",
            no_evidences=False,  # include intersecting gene names
        )

        if result_df is None or len(result_df) == 0:
            print(f"[INFO] g:Profiler returned 0 results")
            return pd.DataFrame()

        print(f"[INFO] g:Profiler returned {len(result_df)} terms")
        print(f"[DEBUG] Result columns: {list(result_df.columns)}")
        if len(result_df) > 0:
            top = result_df.iloc[0]
            print(f"[DEBUG] Top term: {top.get('name', '?')} "
                  f"p_value={top.get('p_value', '?')} "
                  f"intersection_size={top.get('intersection_size', '?')} "
                  f"term_size={top.get('term_size', '?')}")

        return _gprofiler_df_to_enrichr_format(result_df, gp_source)

    except ImportError:
        print(f"[INFO] gprofiler-official not installed; trying REST API…")
    except Exception as e:
        print(f"[WARN] gprofiler-official failed: {e}; trying REST API…")

    # ── Method 2: Raw REST API fallback ───────────────────────────────
    body = {
        "organism": gp_organism,
        "query": genes,
        "sources": [gp_source],
        "user_threshold": sig_cutoff,
        "significance_threshold_method": "fdr",
        "no_evidences": False,
    }

    data = json.dumps(body).encode()
    req = urllib.request.Request(
        GPROFILER_API, data=data, method="POST",
        headers={"Content-Type": "application/json", "Accept": "application/json"})

    with urllib.request.urlopen(req, timeout=60) as resp:
        raw = json.loads(resp.read().decode())

    gp_results = raw.get("result", [])
    if not gp_results:
        print(f"[INFO] g:Profiler REST API returned 0 results")
        return pd.DataFrame()

    # Debug: dump first result's keys and a few field values
    first = gp_results[0]
    if isinstance(first, dict):
        print(f"[DEBUG] REST API result keys: {sorted(first.keys())}")
        print(f"[DEBUG] First result sample: name={first.get('name')}, "
              f"p_value={first.get('p_value')}, "
              f"intersection_size={first.get('intersection_size')}, "
              f"term_size={first.get('term_size')}")
    else:
        print(f"[DEBUG] REST API result type: {type(first).__name__}, value: {first!r:.200s}")
        raise ValueError(f"Unexpected g:Profiler response format: results are {type(first).__name__}, not dict")

    # Convert list of dicts → DataFrame
    result_df = pd.DataFrame(gp_results)
    print(f"[INFO] g:Profiler REST API returned {len(result_df)} terms")
    return _gprofiler_df_to_enrichr_format(result_df, gp_source)


def run_enrichr(genes: List[str], gene_set: str, organism: str) -> pd.DataFrame:
    if gp is None:
        raise ImportError("gseapy is not installed; cannot use Enrichr backend")
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
    """Plot already-selected terms as a bubble plot. df should be pre-filtered."""
    if len(df) == 0:
        raise ValueError("No enrichment results available for plotting")

    dat = df.copy()
    if len(dat) == 0:
        raise ValueError("No enrichment terms available after filtering")

    # Check how many terms pass significance threshold
    n_sig = int((dat["Adjusted P-value"] <= sig_cutoff).sum())

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

    # Anchor color scale: vmin=0, vmax extends to at least -log10(sig_cutoff)
    sig_line = -np.log10(sig_cutoff)  # e.g. 1.3 for 0.05
    data_max = float(dat["minus_log10_fdr"].max())
    norm = Normalize(
        vmin=0.0,
        vmax=max(data_max * 1.05, sig_line * 1.3),
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

    # Add subtitle annotation when nothing is significant
    if n_sig == 0:
        ax.set_title(f"{title}\n(none significant at FDR ≤ {sig_cutoff})",
                      fontsize=font_size + 2)
    else:
        ax.set_title(title)

    ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.grid(axis="x", linestyle="--", alpha=0.3, zorder=1)

    cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.88)
    cbar.set_label("-log10(FDR)")

    # Draw significance threshold line on colorbar
    if 0 < sig_line < norm.vmax:
        cbar.ax.axhline(y=sig_line, color="red", linewidth=1.0, linestyle="--")
        cbar.ax.text(
            1.05, sig_line, f" FDR={sig_cutoff}",
            transform=cbar.ax.get_yaxis_transform(),
            va="center", ha="left", fontsize=font_size - 1.5, color="red",
        )

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
        "Genes",
    ]
    existing = [c for c in cols if c in out.columns]
    out = out[existing].copy()
    ensure_parent_dir(output)
    out.to_csv(output, sep="\t", index=False, quoting=csv.QUOTE_NONNUMERIC)


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
    parser.add_argument("--backend", default="auto",
                        choices=["auto", "gprofiler", "enrichr"],
                        help="Enrichment backend: auto (g:Profiler first, Enrichr fallback), "
                             "gprofiler, or enrichr")
    args = parser.parse_args()

    init_style(
        font_family=getattr(args, 'font_family', None),
        font_size=getattr(args, 'font_size', None),
    )

    organism = normalize_organism(args.organism)
    genes = read_gene_list(args.input)
    print(f"[INFO] Loaded {len(genes)} genes from {args.input}")
    categories = parse_panel(args.library)

    generated_any = False

    for category in categories:
        display_name = DISPLAY_NAME_MAP[category]
        res = pd.DataFrame()

        # Try g:Profiler first (official GO annotations, gene alias mapping)
        if args.backend in ("auto", "gprofiler") and category in GPROFILER_SOURCE_MAP:
            try:
                print(f"[INFO] Running g:Profiler for {display_name} "
                      f"(source={GPROFILER_SOURCE_MAP[category]})…")
                res = run_gprofiler(genes, category, organism,
                                    sig_cutoff=args.sig_cutoff)
                if len(res) > 0:
                    print(f"[INFO] g:Profiler returned {len(res)} terms")
                else:
                    print(f"[INFO] g:Profiler returned 0 terms")
            except Exception as e:
                print(f"[WARN] g:Profiler failed for {category}: {e}")
                res = pd.DataFrame()

        # Fallback to Enrichr if g:Profiler failed or returned nothing
        if len(res) == 0 and args.backend in ("auto", "enrichr"):
            gene_set = resolve_library(category, organism)
            print(f"[INFO] Running Enrichr for {display_name}: {gene_set}")
            try:
                if gp is None:
                    raise ImportError("gseapy not installed; cannot use Enrichr backend")
                res = run_enrichr(genes, gene_set, organism=organism)
            except Exception as e:
                print(f"[ERROR] Enrichr failed for {category}: {e}", file=sys.stderr)
                continue

        out_tsv = f"{args.output_prefix}.{category}.tsv"
        out_plot_tsv = f"{args.output_prefix}.{category}.plot_input.tsv"
        out_png = f"{args.output_prefix}.{category}.bubble.png"
        out_pdf = f"{args.output_prefix}.{category}.bubble.pdf"

        if len(res) == 0:
            print(f"[WARN] No enrichment results for {category}")
            continue

        ensure_parent_dir(out_tsv)
        # Drop internal/Enrichr-specific columns from output TSV
        drop_cols = {"Term_wrapped", "Old P-value", "Old Adjusted P-value",
                     "Odds Ratio", "Combined Score"}
        save_cols = [c for c in res.columns if c not in drop_cols]
        # Use QUOTE_NONNUMERIC to prevent Excel date-corruption of
        # the Overlap column (e.g. "2/23" → "23-Feb")
        res[save_cols].to_csv(out_tsv, sep="\t", index=False,
                              quoting=csv.QUOTE_NONNUMERIC)

        n_total = len(res)
        n_sig = int((res["Adjusted P-value"] <= args.sig_cutoff).sum())
        print(f"[INFO] {display_name}: {n_total} terms total, {n_sig} significant at FDR ≤ {args.sig_cutoff}")
        if n_sig > 0:
            best = res.iloc[0]
            print(f"[INFO]   Most significant: {best['Term_clean']} (FDR={best['Adjusted P-value']:.2e})")

        plot_df = choose_terms_for_plot(res, top_n=args.top_n, sig_cutoff=args.sig_cutoff)
        if len(plot_df) == 0:
            print(f"[WARN] No plottable terms for {category}")
            continue

        print(f"[INFO]   Plotting top {len(plot_df)} terms (FDR range: "
              f"{plot_df['Adjusted P-value'].min():.2e} – {plot_df['Adjusted P-value'].max():.2e})")

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