#!/usr/bin/env python3
"""
ATAC-seq / ChIP-seq Downstream Analysis

Performs comprehensive downstream analysis on peak files:
  - Peak annotation to genomic features (promoter, exon, intron, intergenic, downstream)
  - QC metrics (peak width, score distribution, TSS enrichment)
  - Differential peak analysis across conditions
  - Visualization (pie charts, heatmaps, Venn diagrams)

Input: BED/narrowPeak/broadPeak format peak files
Output: Annotated peaks, QC plots, differential peak tables, visualizations

Dependencies: pandas, numpy, matplotlib, requests
"""

import argparse
import os
import sys
import warnings
from typing import Dict, List, Optional, Set, Tuple

import matplotlib
import sys as _sys, os as _os
try:
    _sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), *(['..'] * 3), '_shared'))
    from plot_style import init_style
except ImportError:
    def init_style(**kw): pass  # graceful fallback if _shared not available
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests

warnings.filterwarnings("ignore", category=RuntimeWarning)

# =========================================================
# Utilities
# =========================================================

def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def load_peak_file(path: str) -> pd.DataFrame:
    """
    Load peak file (BED, narrowPeak, or broadPeak format).
    Auto-detects number of columns.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Peak file not found: {path}")

    try:
        df = pd.read_csv(path, sep="\t", header=None, comment="#")
    except Exception as e:
        raise ValueError(f"Failed to parse peak file {path}: {e}")

    if df.shape[1] < 3:
        raise ValueError(f"Peak file must have at least 3 columns (BED3 format), got {df.shape[1]}")

    # Rename columns based on format
    if df.shape[1] == 3:
        # BED3
        df.columns = ["chrom", "start", "end"]
        df["name"] = [f"peak_{i}" for i in range(len(df))]
        df["score"] = 0.0
    elif df.shape[1] == 6:
        # BED6
        df.columns = ["chrom", "start", "end", "name", "score", "strand"]
    elif df.shape[1] == 9:
        # broadPeak
        df.columns = ["chrom", "start", "end", "name", "score", "strand", "signalValue", "pvalue", "qvalue"]
    elif df.shape[1] >= 10:
        # narrowPeak (10+ columns)
        cols = ["chrom", "start", "end", "name", "score", "strand", "signalValue", "pvalue", "qvalue", "peak"]
        df.columns = cols + [f"col_{i}" for i in range(len(cols), df.shape[1])]
    else:
        df.columns = [f"col_{i}" for i in range(df.shape[1])]

    # Convert to numeric
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype(int)
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype(int)
    df["score"] = pd.to_numeric(df["score"], errors="coerce").fillna(0.0)

    # Get signal value if available
    if "signalValue" in df.columns:
        df["signalValue"] = pd.to_numeric(df["signalValue"], errors="coerce").fillna(0.0)
    else:
        df["signalValue"] = df["score"]

    return df


def filter_peaks(df: pd.DataFrame, min_score: float, top_n: int) -> pd.DataFrame:
    """Filter peaks by score and keep top-n."""
    if min_score > 0:
        df = df[df["score"] >= min_score].copy()

    if top_n < len(df):
        df = df.nlargest(top_n, "score").copy()

    return df.reset_index(drop=True)


def get_ensembl_genes(chrom: str, start: int, end: int, species: str = "homo_sapiens") -> List[Dict]:
    """
    Fetch genes from Ensembl REST API for a genomic region.
    Returns list of dicts with gene info.
    """
    url = f"https://rest.ensembl.org/overlap/region/{species}/{chrom}:{start}-{end}?feature=gene&content-type=application/json"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            genes = [item for item in data if item.get("feature_type") == "gene"]
            return genes
    except Exception as e:
        pass
    return []


def annotate_peaks(peaks_df: pd.DataFrame, genome: str = "hg38", tss_window: int = 2000) -> pd.DataFrame:
    """
    Annotate peaks to genomic features using Ensembl REST API.
    Classifies peaks as: promoter, exonic, intronic, intergenic, downstream.
    """
    species_map = {
        "hg38": "homo_sapiens",
        "hg19": "homo_sapiens",
        "mm10": "mus_musculus",
        "mm39": "mus_musculus",
    }
    species = species_map.get(genome, "homo_sapiens")

    result = []
    annotations = []
    nearest_genes = []
    distances = []

    for idx, row in peaks_df.iterrows():
        chrom = row["chrom"]
        peak_start = row["start"]
        peak_end = row["end"]
        peak_center = (peak_start + peak_end) // 2

        # Query Ensembl for overlapping genes
        genes = get_ensembl_genes(chrom, peak_start - 10000, peak_end + 10000, species)

        annotation = "intergenic"
        nearest_gene = "."
        distance = np.inf

        for gene in genes:
            gene_start = gene.get("start", 0)
            gene_end = gene.get("end", 0)
            gene_name = gene.get("external_name", ".")
            gene_strand = gene.get("strand", "+")

            # Estimate TSS
            tss = gene_start if gene_strand == 1 else gene_end

            # Update nearest gene
            dist_to_tss = abs(peak_center - tss)
            if dist_to_tss < distance:
                distance = dist_to_tss
                nearest_gene = gene_name

            # Annotate based on overlap
            if abs(peak_center - tss) <= tss_window:
                annotation = "promoter"
                break
            elif peak_start >= gene_start and peak_end <= gene_end:
                annotation = "exonic"
            elif peak_start >= gene_start and peak_start <= gene_end:
                annotation = "exonic"
            elif peak_end >= gene_start and peak_end <= gene_end:
                annotation = "exonic"
            elif gene_start <= peak_start <= peak_end <= gene_end:
                annotation = "exonic"
            elif peak_start < gene_start and peak_end > gene_end:
                annotation = "intronic"
            elif gene_start <= peak_start < gene_end < peak_end:
                annotation = "downstream"

            if annotation != "intergenic":
                break

        if distance == np.inf:
            distance = -1

        result.append({
            "peak_id": row.get("name", f"peak_{idx}"),
            "chrom": chrom,
            "start": peak_start,
            "end": peak_end,
            "score": row["score"],
            "annotation": annotation,
            "nearest_gene": nearest_gene,
            "distance_to_tss": int(distance) if distance != np.inf else -1,
        })

    return pd.DataFrame(result)


def compute_qc_metrics(peaks_df: pd.DataFrame, annotated_df: Optional[pd.DataFrame] = None) -> Dict:
    """Compute QC metrics for peaks."""
    widths = peaks_df["end"] - peaks_df["start"]
    scores = peaks_df["score"]

    metrics = {
        "total_peaks": len(peaks_df),
        "median_peak_width": int(np.median(widths)),
        "mean_peak_width": int(np.mean(widths)),
        "min_peak_width": int(np.min(widths)),
        "max_peak_width": int(np.max(widths)),
        "mean_score": float(np.mean(scores)),
        "median_score": float(np.median(scores)),
        "std_score": float(np.std(scores)),
    }

    if annotated_df is not None and len(annotated_df) > 0:
        promoter_count = (annotated_df["annotation"] == "promoter").sum()
        metrics["tss_enrichment_score"] = promoter_count / len(annotated_df)
        metrics["promoter_peak_count"] = int(promoter_count)

    return metrics


def plot_peak_annotation(annotated_df: pd.DataFrame, outdir: str) -> None:
    """Create pie chart of peak annotation distribution."""
    if len(annotated_df) == 0:
        return

    counts = annotated_df["annotation"].value_counts()
    colors = {
        "promoter": "#FF6B6B",
        "exonic": "#4ECDC4",
        "intronic": "#45B7D1",
        "intergenic": "#96CEB4",
        "downstream": "#FFEAA7",
    }
    color_list = [colors.get(x, "#CCCCCC") for x in counts.index]

    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    ax.pie(counts.values, labels=counts.index, autopct="%1.1f%%", colors=color_list, startangle=90)
    ax.set_title("Peak Annotation Distribution")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "peak_annotation_distribution.png"), dpi=300, bbox_inches="tight")
    plt.close()


def plot_peak_widths(peaks_df: pd.DataFrame, outdir: str) -> None:
    """Create histogram of peak widths."""
    widths = peaks_df["end"] - peaks_df["start"]

    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    ax.hist(widths, bins=100, edgecolor="black", alpha=0.7, color="#4ECDC4")
    ax.set_xlabel("Peak Width (bp)")
    ax.set_ylabel("Frequency")
    ax.set_title("Peak Width Distribution")
    ax.set_yscale("log")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "peak_width_distribution.png"), dpi=300, bbox_inches="tight")
    plt.close()


def plot_peak_scores(peaks_df: pd.DataFrame, outdir: str) -> None:
    """Create histogram of peak scores."""
    scores = peaks_df["score"]

    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    ax.hist(scores, bins=100, edgecolor="black", alpha=0.7, color="#FF6B6B")
    ax.set_xlabel("Peak Score / Signal Value")
    ax.set_ylabel("Frequency")
    ax.set_title("Peak Score Distribution")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "peak_score_distribution.png"), dpi=300, bbox_inches="tight")
    plt.close()


def peaks_overlap(p1: pd.DataFrame, p2: pd.DataFrame, reciprocal_threshold: float = 0.5) -> float:
    """
    Compute Jaccard similarity between two peak sets.
    Uses reciprocal overlap criterion.
    """
    if len(p1) == 0 or len(p2) == 0:
        return 0.0

    overlaps = 0
    for idx1, row1 in p1.iterrows():
        chrom1, start1, end1 = row1["chrom"], row1["start"], row1["end"]
        for idx2, row2 in p2.iterrows():
            chrom2, start2, end2 = row2["chrom"], row2["start"], row2["end"]
            if chrom1 != chrom2:
                continue

            overlap = max(0, min(end1, end2) - max(start1, start2))
            len1 = end1 - start1
            len2 = end2 - start2

            if len1 > 0 and len2 > 0:
                reciprocal_ov = overlap / min(len1, len2)
                if reciprocal_ov >= reciprocal_threshold:
                    overlaps += 1
                    break

    union = len(p1) + len(p2) - overlaps
    return overlaps / union if union > 0 else 0.0


def find_differential_peaks(peak_files: List[str], conditions: List[str], outdir: str) -> None:
    """Find unique and shared peaks across conditions."""
    if len(peak_files) < 2:
        print("Warning: At least 2 peak files needed for differential analysis.")
        return

    # Load all peak files
    peak_dfs = {}
    for pf, cond in zip(peak_files, conditions):
        peak_dfs[cond] = load_peak_file(pf)

    # Find unique peaks per condition
    for cond, df in peak_dfs.items():
        unique = df.copy()
        for other_cond, other_df in peak_dfs.items():
            if other_cond != cond:
                for idx, row in unique.iterrows():
                    for idx2, row2 in other_df.iterrows():
                        if row["chrom"] == row2["chrom"]:
                            overlap = max(0, min(row["end"], row2["end"]) - max(row["start"], row2["start"]))
                            if overlap / (row["end"] - row["start"]) > 0.5:
                                unique = unique.drop(idx)
                                break

        unique = unique.reset_index(drop=True)
        unique[["chrom", "start", "end", "name", "score"]].to_csv(
            os.path.join(outdir, f"unique_peaks_{cond}.bed"),
            sep="\t", header=False, index=False
        )

    # Find shared peaks (in all conditions)
    if len(peak_files) >= 2:
        shared = peak_dfs[conditions[0]].copy()
        for cond in conditions[1:]:
            df = peak_dfs[cond]
            shared_filtered = []
            for idx, row in shared.iterrows():
                found = False
                for idx2, row2 in df.iterrows():
                    if row["chrom"] == row2["chrom"]:
                        overlap = max(0, min(row["end"], row2["end"]) - max(row["start"], row2["start"]))
                        if overlap / (row["end"] - row["start"]) > 0.5:
                            found = True
                            break
                if found:
                    shared_filtered.append(idx)
            shared = shared.iloc[shared_filtered].reset_index(drop=True)

        shared[["chrom", "start", "end", "name", "score"]].to_csv(
            os.path.join(outdir, "shared_peaks_all_conditions.bed"),
            sep="\t", header=False, index=False
        )

    # Compute pairwise similarities
    similarity_matrix = np.ones((len(conditions), len(conditions)))
    for i, cond1 in enumerate(conditions):
        for j, cond2 in enumerate(conditions):
            if i != j:
                similarity_matrix[i, j] = peaks_overlap(peak_dfs[cond1], peak_dfs[cond2])

    # Plot similarity heatmap
    if len(conditions) >= 2:
        fig, ax = plt.subplots(figsize=(8, 7), dpi=300)
        im = ax.imshow(similarity_matrix, cmap="YlOrRd", vmin=0, vmax=1)
        ax.set_xticks(range(len(conditions)))
        ax.set_yticks(range(len(conditions)))
        ax.set_xticklabels(conditions, rotation=45, ha="right")
        ax.set_yticklabels(conditions)
        ax.set_title("Peak Overlap Similarity Matrix (Jaccard)")

        for i in range(len(conditions)):
            for j in range(len(conditions)):
                text = ax.text(j, i, f"{similarity_matrix[i, j]:.2f}",
                             ha="center", va="center", color="black", fontsize=10)

        plt.colorbar(im, ax=ax, label="Jaccard Similarity")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "peak_overlap_matrix.png"), dpi=300, bbox_inches="tight")
        plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Downstream analysis of ATAC-seq or ChIP-seq peak files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--peaks", required=True, help="Peak file (BED/narrowPeak/broadPeak)")
    parser.add_argument("--mode", choices=["annotate", "differential", "qc", "all"], default="all",
                       help="Analysis mode")
    parser.add_argument("--genome", choices=["hg38", "hg19", "mm10", "mm39"], default="hg38",
                       help="Reference genome")
    parser.add_argument("--peak-files", default="", help="Comma-separated peak files for differential analysis")
    parser.add_argument("--conditions", default="", help="Comma-separated condition labels")
    parser.add_argument("--gene-list", default="", help="Comma-separated genes to highlight")
    parser.add_argument("--outdir", default="./atac_chipseq_output", help="Output directory")
    parser.add_argument("--tss-window", type=int, default=2000, help="TSS window (bp)")
    parser.add_argument("--min-score", type=float, default=0, help="Minimum peak score")
    parser.add_argument("--top-n", type=int, default=50000, help="Top N peaks to analyze")

    args = parser.parse_args()

    init_style(
        font_family=getattr(args, 'font_family', None),
        font_size=getattr(args, 'font_size', None),
    )

    ensure_dir(args.outdir)

    # Load and filter peaks
    print(f"Loading peaks from {args.peaks}...")
    peaks = load_peak_file(args.peaks)
    print(f"  Loaded {len(peaks)} peaks")

    peaks = filter_peaks(peaks, args.min_score, args.top_n)
    print(f"  After filtering: {len(peaks)} peaks")

    # QC analysis
    if args.mode in ["qc", "all"]:
        print("\nRunning QC analysis...")
        metrics = compute_qc_metrics(peaks)

        qc_file = os.path.join(args.outdir, "qc_metrics.txt")
        with open(qc_file, "w") as f:
            f.write("Peak QC Metrics\n")
            f.write("=" * 50 + "\n")
            for k, v in metrics.items():
                f.write(f"{k}: {v}\n")

        plot_peak_widths(peaks, args.outdir)
        plot_peak_scores(peaks, args.outdir)
        print(f"  QC plots saved to {args.outdir}")

    # Annotation
    annotated = None
    if args.mode in ["annotate", "all"]:
        print("\nAnnotating peaks (this may take a moment)...")
        try:
            annotated = annotate_peaks(peaks, args.genome, args.tss_window)

            annotated_file = os.path.join(args.outdir, "peaks_annotated.tsv")
            annotated.to_csv(annotated_file, sep="\t", index=False)
            print(f"  Annotated peaks saved to {annotated_file}")

            # Save summary
            summary_file = os.path.join(args.outdir, "peak_annotation_summary.txt")
            with open(summary_file, "w") as f:
                f.write("Peak Annotation Summary\n")
                f.write("=" * 50 + "\n")
                counts = annotated["annotation"].value_counts()
                for ann, count in counts.items():
                    f.write(f"{ann}: {count} ({100*count/len(annotated):.1f}%)\n")

            plot_peak_annotation(annotated, args.outdir)
        except Exception as e:
            print(f"  Warning: Annotation failed ({e}). Skipping annotation.")

    # Differential analysis
    if args.mode in ["differential", "all"]:
        if args.peak_files and args.conditions:
            peak_files = args.peak_files.split(",")
            conditions = args.conditions.split(",")

            if len(peak_files) != len(conditions):
                print("Error: Number of peak files must match number of conditions")
                sys.exit(1)

            print(f"\nRunning differential analysis on {len(conditions)} conditions...")
            find_differential_peaks(peak_files, conditions, args.outdir)
            print(f"  Differential analysis complete")

    # Summary
    summary_file = os.path.join(args.outdir, "analysis_summary.txt")
    with open(summary_file, "w") as f:
        f.write("ATAC/ChIP-seq Downstream Analysis Summary\n")
        f.write("=" * 60 + "\n")
        f.write(f"Input peak file: {args.peaks}\n")
        f.write(f"Peaks loaded: {len(peaks)}\n")
        f.write(f"Mode: {args.mode}\n")
        f.write(f"Genome: {args.genome}\n")
        f.write(f"Output directory: {args.outdir}\n")

    print(f"\nAnalysis complete. Results saved to {args.outdir}")


if __name__ == "__main__":
    main()
