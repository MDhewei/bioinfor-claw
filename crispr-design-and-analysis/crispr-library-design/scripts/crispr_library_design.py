#!/usr/bin/env python3
"""
Design a pooled CRISPR sgRNA library for target genes.

Designs sgRNAs for each gene, selects top guides by efficiency,
includes non-targeting controls, and outputs oligo order file.
"""

import argparse
import os
import re
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

import matplotlib.pyplot as plt
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), '..', '..', '..', '_shared'))
from plot_style import init_style
import numpy as np
import pandas as pd
import requests


IUPAC_CODES = {
    "N": "ACGT",
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
}

# Simplified scoring based on nucleotide preferences (Doench 2016)
DOENCH_WEIGHTS = {
    1: {"T": 0.30, "G": -0.28, "C": -0.31, "A": 0.00},
    2: {"A": -0.22, "G": 0.00, "C": -0.26, "T": 0.00},
}


@dataclass
class SgRNADesign:
    """Single sgRNA design."""
    gene: str
    spacer: str
    pam: str
    strand: str
    chromosome: str
    position: int
    gc_content: float
    efficiency_score: float
    rank: int


def reverse_complement(seq: str) -> str:
    """Reverse complement DNA sequence."""
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(complement.get(b, "N") for b in reversed(seq))


def iupac_to_regex(pam: str) -> str:
    """Convert IUPAC PAM to regex pattern."""
    pattern = ""
    for base in pam.upper():
        if base in IUPAC_CODES:
            pattern += "[" + IUPAC_CODES[base] + "]"
        else:
            pattern += base
    return pattern


def calculate_gc(seq: str) -> float:
    """Calculate GC content."""
    if len(seq) == 0:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq)


def check_problematic_sequences(seq: str) -> bool:
    """Check for problematic sequences (poly-T, BsmBI sites, etc.)."""
    # Poly-T
    if "TTTT" in seq:
        return True
    # BsmBI site (CACCG)
    if "CACCG" in seq:
        return True
    # Excessive homopolymer
    for base in "ACGT":
        if base * 5 in seq:
            return True
    return False


def score_sgrna_doench(spacer: str, gc_content: float) -> float:
    """
    Score sgRNA using simplified Doench 2016 model.
    """
    score = 0.5  # Baseline

    # GC content penalty
    gc_penalty = abs(gc_content - 0.5)
    score -= gc_penalty * 0.3

    # Position-specific nucleotide weights (simplified)
    for pos in [0, 1]:
        if pos < len(spacer) and pos in DOENCH_WEIGHTS:
            base = spacer[pos]
            score += DOENCH_WEIGHTS[pos].get(base, 0.0)

    # Seed region quality (17-20)
    seed = spacer[-4:].upper()
    seed_gc = calculate_gc(seed)
    if 0.25 <= seed_gc <= 0.75:
        score += 0.2
    else:
        score -= 0.1

    return max(0.0, min(1.0, score))


def fetch_gene_sequence_ensembl(
    gene: str, species: str = "human"
) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """Fetch gene coding sequence from Ensembl."""
    species_code = "homo_sapiens" if species.lower() == "human" else "mus_musculus"

    try:
        # Lookup gene
        lookup_url = f"https://rest.ensembl.org/lookup/symbol/{species_code}/{gene}"
        resp = requests.get(lookup_url, headers={"Content-Type": "application/json"}, timeout=10)
        resp.raise_for_status()
        data = resp.json()

        chrom = data.get("seq_region_name")
        start = data.get("start")
        end = data.get("end")

        # Fetch sequence
        seq_url = f"https://rest.ensembl.org/sequence/region/{species_code}/{chrom}:{start}-{end}:1"
        seq_resp = requests.get(seq_url, headers={"Content-Type": "application/json"}, timeout=10)
        seq_resp.raise_for_status()
        sequence = seq_resp.json().get("seq", "").upper()

        return sequence, chrom, str(start)
    except Exception as e:
        print(f"Warning: Could not fetch {gene} from Ensembl: {e}", file=sys.stderr)
        return None, None, None


def find_pams_in_sequence(sequence: str, pam: str, spacer_len: int = 20) -> List[Tuple[int, str, str]]:
    """Find PAM sites in sequence."""
    results = []
    pam_pattern = iupac_to_regex(pam)

    # Forward strand
    for match in re.finditer(pam_pattern, sequence):
        pam_pos = match.start()
        if pam_pos >= spacer_len:
            spacer = sequence[pam_pos - spacer_len : pam_pos]
            if "N" not in spacer and not check_problematic_sequences(spacer):
                results.append((pam_pos - spacer_len, "+", spacer))

    # Reverse strand
    rc_seq = reverse_complement(sequence)
    for match in re.finditer(pam_pattern, rc_seq):
        pam_pos = match.start()
        if pam_pos >= spacer_len:
            spacer = rc_seq[pam_pos - spacer_len : pam_pos]
            if "N" not in spacer and not check_problematic_sequences(spacer):
                genomic_pos = len(sequence) - pam_pos - spacer_len
                results.append((genomic_pos, "-", spacer))

    return results


def design_sgrnas_for_gene(
    gene: str,
    sequence: str,
    chromosome: str,
    start_pos: str,
    pam: str = "NGG",
    spacer_len: int = 20,
    gc_min: float = 0.3,
    gc_max: float = 0.8,
) -> List[SgRNADesign]:
    """Design sgRNAs for a single gene."""
    designs = []

    pam_sites = find_pams_in_sequence(sequence, pam, spacer_len)

    for pos, strand, spacer in pam_sites:
        gc = calculate_gc(spacer)

        # Filter GC content
        if not (gc_min <= gc <= gc_max):
            continue

        score = score_sgrna_doench(spacer, gc)

        design = SgRNADesign(
            gene=gene,
            spacer=spacer,
            pam=pam,
            strand=strand,
            chromosome=chromosome,
            position=int(start_pos) + pos if start_pos else 0,
            gc_content=gc,
            efficiency_score=score,
            rank=0,
        )
        designs.append(design)

    # Sort by score and rank
    designs.sort(key=lambda d: d.efficiency_score, reverse=True)
    for i, d in enumerate(designs, 1):
        d.rank = i

    return designs


def generate_non_targeting_controls(
    n_controls: int, spacer_len: int = 20, gc_min: float = 0.3, gc_max: float = 0.8
) -> List[SgRNADesign]:
    """Generate non-targeting control sgRNAs."""
    controls = []
    bases = "ACGT"
    np.random.seed(42)

    attempts = 0
    max_attempts = n_controls * 10

    while len(controls) < n_controls and attempts < max_attempts:
        # Generate random spacer
        spacer = "".join(np.random.choice(list(bases), spacer_len))

        # Filter
        if check_problematic_sequences(spacer):
            attempts += 1
            continue

        gc = calculate_gc(spacer)
        if not (gc_min <= gc <= gc_max):
            attempts += 1
            continue

        score = score_sgrna_doench(spacer, gc)

        control = SgRNADesign(
            gene=f"non-targeting_{len(controls):04d}",
            spacer=spacer,
            pam="NGG",
            strand="+",
            chromosome="control",
            position=0,
            gc_content=gc,
            efficiency_score=score,
            rank=len(controls) + 1,
        )
        controls.append(control)
        attempts += 1

    if len(controls) < n_controls:
        print(f"Warning: Only generated {len(controls)}/{n_controls} controls", file=sys.stderr)

    return controls


def load_gene_list(genes_input: str) -> List[str]:
    """Parse gene list from comma-separated string or file."""
    if "," in genes_input:
        return [g.strip() for g in genes_input.split(",") if g.strip()]
    elif os.path.isfile(genes_input):
        with open(genes_input) as f:
            return [line.strip() for line in f if line.strip()]
    else:
        return [genes_input]


def save_library(
    all_designs: List[SgRNADesign],
    outdir: str,
    library_name: str,
    vector_prefix: str,
    vector_suffix: str,
    output_format: str = "tsv",
):
    """Save library design to file(s)."""
    os.makedirs(outdir, exist_ok=True)

    # Prepare records
    records = []
    for i, design in enumerate(all_designs, 1):
        full_oligo = vector_prefix + design.spacer + vector_suffix
        records.append({
            "guide_id": f"{design.gene}_guide{design.rank}",
            "gene": design.gene,
            "spacer_sequence": design.spacer,
            "full_oligo": full_oligo,
            "strand": design.strand,
            "chromosome": design.chromosome,
            "position": design.position,
            "gc_content": round(design.gc_content, 3),
            "efficiency_score": round(design.efficiency_score, 3),
            "guide_rank": design.rank,
            "library": library_name,
        })

    df = pd.DataFrame(records)

    # Save in requested format(s)
    if output_format in ["tsv", "both"]:
        tsv_path = os.path.join(outdir, f"{library_name}_library.tsv")
        df.to_csv(tsv_path, sep="\t", index=False)
        print(f"Saved library (TSV) to {tsv_path}")

    if output_format in ["csv", "both"]:
        csv_path = os.path.join(outdir, f"{library_name}_library.csv")
        df.to_csv(csv_path, sep=",", index=False)
        print(f"Saved library (CSV) to {csv_path}")

    if output_format in ["fasta", "both"]:
        fasta_path = os.path.join(outdir, f"{library_name}_library.fasta")
        with open(fasta_path, "w") as f:
            for _, row in df.iterrows():
                f.write(f">{row['guide_id']}\n{row['spacer_sequence']}\n")
        print(f"Saved library (FASTA) to {fasta_path}")

    return df


def generate_plots(all_designs: List[SgRNADesign], outdir: str, library_name: str):
    """Generate summary plots."""
    gcs = [d.gc_content for d in all_designs]
    scores = [d.efficiency_score for d in all_designs]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # GC content histogram
    axes[0, 0].hist(gcs, bins=30, color="steelblue", edgecolor="black")
    axes[0, 0].axvline(0.5, color="red", linestyle="--", label="GC=50%")
    axes[0, 0].set_xlabel("GC Content")
    axes[0, 0].set_ylabel("Count")
    axes[0, 0].set_title("GC Content Distribution")
    axes[0, 0].legend()

    # Efficiency score histogram
    axes[0, 1].hist(scores, bins=30, color="coral", edgecolor="black")
    axes[0, 1].set_xlabel("Efficiency Score")
    axes[0, 1].set_ylabel("Count")
    axes[0, 1].set_title("Guide Efficiency Score Distribution")

    # GC vs Score scatter
    axes[1, 0].scatter(gcs, scores, alpha=0.3, s=10)
    axes[1, 0].set_xlabel("GC Content")
    axes[1, 0].set_ylabel("Efficiency Score")
    axes[1, 0].set_title("GC Content vs Efficiency")

    # Summary statistics
    ax_stats = axes[1, 1]
    ax_stats.axis("off")
    stats_text = (
        f"Total Guides: {len(all_designs)}\n"
        f"Mean GC: {np.mean(gcs):.2%}\n"
        f"Mean Score: {np.mean(scores):.3f}\n"
        f"GC Range: {np.min(gcs):.2%} - {np.max(gcs):.2%}\n"
        f"Score Range: {np.min(scores):.3f} - {np.max(scores):.3f}"
    )
    ax_stats.text(0.1, 0.5, stats_text, fontsize=12, verticalalignment="center",
                  family="monospace", bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    plt.tight_layout()
    plot_path = os.path.join(outdir, f"{library_name}_summary.png")
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    print(f"Saved plot to {plot_path}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Design pooled CRISPR sgRNA library")
    parser.add_argument("--genes", type=str, required=True, help="Comma-separated gene symbols or file")
    parser.add_argument("--species", type=str, default="human", choices=["human", "mouse"])
    parser.add_argument("--n-guides-per-gene", type=int, default=4)
    parser.add_argument("--n-controls", type=int, default=500)
    parser.add_argument("--guide-len", type=int, default=20)
    parser.add_argument("--pam", type=str, default="NGG")
    parser.add_argument("--library-name", type=str, default="CRISPRLibrary")
    parser.add_argument("--vector-prefix", type=str, default="ACCG")
    parser.add_argument("--vector-suffix", type=str, default="GTTTTAGAGCTA")
    parser.add_argument("--outdir", type=str, default="library_design")
    parser.add_argument("--output-format", type=str, default="tsv", choices=["csv", "tsv", "fasta"])
    parser.add_argument("--genome", type=str, default="hg38")
    parser.add_argument("--targeting-exon-only", action="store_true")
    parser.add_argument("--no-t4-stretch", action="store_true")
    parser.add_argument("--gc-min", type=float, default=0.3)
    parser.add_argument("--gc-max", type=float, default=0.8)

    args = parser.parse_args()

    init_style(
        font_family=getattr(args, 'font_family', None),
        font_size=getattr(args, 'font_size', None),
    )

    # Parse gene list
    genes = load_gene_list(args.genes)
    print(f"Loaded {len(genes)} genes")

    all_designs = []
    successful_genes = 0

    # Design sgRNAs for each gene
    for gene in genes:
        print(f"Designing sgRNAs for {gene}...", end=" ")
        sequence, chromosome, start_pos = fetch_gene_sequence_ensembl(gene, args.species)

        if not sequence:
            print(f"SKIP (no sequence)")
            continue

        designs = design_sgrnas_for_gene(
            gene=gene,
            sequence=sequence,
            chromosome=chromosome or "unknown",
            start_pos=start_pos or "0",
            pam=args.pam,
            spacer_len=args.guide_len,
            gc_min=args.gc_min,
            gc_max=args.gc_max,
        )

        if not designs:
            print(f"SKIP (no valid designs)")
            continue

        # Select top N guides per gene
        top_designs = designs[: args.n_guides_per_gene]
        all_designs.extend(top_designs)
        successful_genes += 1
        print(f"OK ({len(top_designs)} guides)")

    if not all_designs:
        print("Error: No designs generated for any gene", file=sys.stderr)
        sys.exit(1)

    print(f"Successfully designed sgRNAs for {successful_genes}/{len(genes)} genes")

    # Generate non-targeting controls
    print(f"Generating {args.n_controls} non-targeting controls...", end=" ")
    controls = generate_non_targeting_controls(
        args.n_controls, args.guide_len, args.gc_min, args.gc_max
    )
    all_designs.extend(controls)
    print(f"OK ({len(controls)} controls)")

    # Save library
    df = save_library(
        all_designs,
        args.outdir,
        args.library_name,
        args.vector_prefix,
        args.vector_suffix,
        args.output_format,
    )

    # Generate plots
    generate_plots(all_designs, args.outdir, args.library_name)

    # Summary report
    print("\n=== Library Summary ===")
    print(f"Total guides: {len(all_designs)}")
    print(f"Targeting guides: {len([d for d in all_designs if not d.gene.startswith('non-targeting')])}")
    print(f"Non-targeting controls: {len(controls)}")
    print(f"Mean GC content: {np.mean([d.gc_content for d in all_designs]):.2%}")
    print(f"Mean efficiency score: {np.mean([d.efficiency_score for d in all_designs]):.3f}")

    # Per-gene guide counts
    gene_counts = {}
    for design in all_designs:
        if not design.gene.startswith("non-targeting"):
            gene_counts[design.gene] = gene_counts.get(design.gene, 0) + 1

    print(f"\nGuides per gene: min={min(gene_counts.values()) if gene_counts else 0}, "
          f"max={max(gene_counts.values()) if gene_counts else 0}, "
          f"mean={np.mean(list(gene_counts.values())) if gene_counts else 0:.1f}")


if __name__ == "__main__":
    main()
