#!/usr/bin/env python3
"""
Design pegRNAs and nicking sgRNAs for prime editing.

Given a target gene and desired edit (SNV, insertion, deletion), identifies
PE-compatible protospacers, designs RT templates and PBS sequences, and scores
guides by predicted efficiency.
"""

import argparse
import json
import os
import re
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import sys as _sys, os as _os
try:
    _sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), *(['..'] * 3), '_shared'))
    from plot_style import init_style
except ImportError:
    def init_style(**kw): pass  # graceful fallback if _shared not available
import numpy as np
import pandas as pd
import requests


# Nearest-neighbor Tm lookup table (simplified, all dinucleotides)
TM_LOOKUP = {
    "AA": 9.6, "TT": 9.6, "AT": 9.6, "TA": 9.6,
    "CA": 13.8, "GT": 13.8, "CT": 10.4, "AG": 10.4,
    "GA": 13.8, "TC": 13.8, "GC": 14.4, "CG": 14.4,
    "GG": 13.8, "CC": 13.8, "AC": 12.4, "TG": 12.4,
    "GTA": 9.4, "TAC": 9.4, "CAC": 12.9, "GTG": 12.9,
}

IUPAC_CODES = {
    "N": "ACGT",
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
}


@dataclass
class PeguideDesign:
    """Single pegRNA design."""
    spacer: str
    strand: str
    pam: str
    pbs: str
    rt_template: str
    pbs_tm: float
    edit_distance_from_spacer: int
    gc_content: float
    efficiency_score: float
    nick_sgrna: Optional[str] = None
    nick_pam: Optional[str] = None


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


def calculate_tm(seq: str) -> float:
    """Simple nearest-neighbor Tm calculation."""
    if len(seq) < 2:
        return 0.0
    tm = 0.0
    for i in range(len(seq) - 1):
        dinuc = seq[i : i + 2].upper()
        tm += TM_LOOKUP.get(dinuc, 12.0)
    # Basic formula: Tm = 4(G+C) + 2(A+T) for short oligos
    gc = seq.count("G") + seq.count("C")
    at = seq.count("A") + seq.count("T")
    tm_simple = 4 * gc + 2 * at
    return (tm_simple + tm) / 2.0


def fetch_sequence_ensembl(
    gene: str, species: str = "human", flank: int = 500
) -> Tuple[Optional[str], Optional[str]]:
    """Fetch gene sequence from Ensembl REST API."""
    species_code = "homo_sapiens" if species.lower() == "human" else "mus_musculus"
    lookup_url = f"https://rest.ensembl.org/lookup/symbol/{species_code}/{gene}"

    try:
        resp = requests.get(lookup_url, headers={"Content-Type": "application/json"}, timeout=10)
        resp.raise_for_status()
        data = resp.json()
        chrom = data.get("seq_region_name")
        start = max(1, data.get("start", 1) - flank)
        end = data.get("end", 1) + flank

        seq_url = f"https://rest.ensembl.org/sequence/region/{species_code}/{chrom}:{start}-{end}:1"
        seq_resp = requests.get(seq_url, headers={"Content-Type": "application/json"}, timeout=10)
        seq_resp.raise_for_status()
        sequence = seq_resp.json().get("seq", "").upper()
        return sequence, chrom
    except Exception as e:
        print(f"Error fetching from Ensembl: {e}", file=sys.stderr)
        return None, None


def find_pams(sequence: str, pam: str, spacer_len: int = 20) -> List[Tuple[int, str, str]]:
    """Find PAM sites and return (position, strand, full_protospacer)."""
    results = []
    pam_pattern = iupac_to_regex(pam)

    # Forward strand
    for match in re.finditer(pam_pattern, sequence):
        pam_pos = match.start()
        if pam_pos >= spacer_len:
            spacer = sequence[pam_pos - spacer_len : pam_pos]
            if "N" not in spacer:
                results.append((pam_pos - spacer_len, "+", spacer))

    # Reverse strand
    rc_seq = reverse_complement(sequence)
    for match in re.finditer(pam_pattern, rc_seq):
        pam_pos = match.start()
        if pam_pos >= spacer_len:
            spacer = rc_seq[pam_pos - spacer_len : pam_pos]
            if "N" not in spacer:
                genomic_pos = len(sequence) - pam_pos - spacer_len
                results.append((genomic_pos, "-", spacer))

    return results


def design_pbs(spacer: str, pbs_min: int, pbs_max: int) -> Tuple[str, float]:
    """Design PBS (reverse complement of 3' end) and calculate Tm."""
    best_pbs = ""
    best_tm = float("inf")
    best_dist = pbs_max - pbs_min

    for pbs_len in range(pbs_min, pbs_max + 1):
        pbs_target = spacer[-pbs_len:]
        pbs = reverse_complement(pbs_target)
        tm = calculate_tm(pbs)
        dist_from_opt = abs(tm - 33.0)

        if dist_from_opt < best_dist:
            best_pbs = pbs
            best_tm = tm
            best_dist = dist_from_opt

    return best_pbs, best_tm


def design_rt_template(
    wt_seq: str, edit_seq: str, spacer: str, edit_pos: int
) -> str:
    """Design RT template that encodes the edit."""
    # RT template should include edit context
    # Simple approach: extract region around edit site, incorporate edit_seq
    rt_len = min(30, len(edit_seq) + 10)
    rt_template = edit_seq[:rt_len]
    return rt_template


def score_guide(spacer: str, pbs_tm: float, rt_template: str) -> float:
    """Score guide by efficiency heuristics."""
    score = 0.0

    # GC content: optimal 40-70%
    gc = (spacer.count("G") + spacer.count("C")) / len(spacer)
    gc_score = 1.0 if 0.4 <= gc <= 0.7 else 0.5
    score += gc_score * 0.3

    # PBS Tm: optimal 30-37°C
    tm_dist = abs(pbs_tm - 33.0)
    tm_score = 1.0 if tm_dist < 5.0 else max(0.3, 1.0 - tm_dist / 20.0)
    score += tm_score * 0.3

    # No poly-T or other problematic sequences
    if "TTTT" not in spacer and "AAAA" not in spacer:
        score += 0.2
    else:
        score -= 0.1

    # Seed region quality (17-20bp): prefer GC-rich
    seed = spacer[-17:]
    seed_gc = (seed.count("G") + seed.count("C")) / len(seed)
    seed_score = seed_gc / 1.0
    score += min(seed_score, 1.0) * 0.2

    return min(max(score, 0.0), 1.0)


def design_pegguides(
    target_seq: str,
    edit_type: str,
    wt_seq: str,
    edit_seq: str,
    edit_pos: int,
    pam: str = "NGG",
    spacer_len: int = 20,
    pbs_range: Tuple[int, int] = (8, 17),
    rt_range: Tuple[int, int] = (10, 40),
    nick_range: Tuple[int, int] = (40, 100),
    editor: str = "PE3",
    top_n: int = 20,
) -> List[PeguideDesign]:
    """Design pegRNAs and nicking sgRNAs."""
    guides = []

    # Find candidate protospacers near edit site
    pam_sites = find_pams(target_seq, pam, spacer_len)

    for spacer_start, strand, spacer in pam_sites:
        # Filter: must be close to edit site (within 30bp)
        dist_from_edit = abs(spacer_start + spacer_len // 2 - edit_pos)
        if dist_from_edit > 30:
            continue

        # Filter: no poly-T
        if "TTTT" in spacer:
            continue

        # Design PBS
        pbs, pbs_tm = design_pbs(spacer, pbs_range[0], pbs_range[1])

        # Design RT template
        rt_template = design_rt_template(wt_seq, edit_seq, spacer, edit_pos)

        # Score
        gc = (spacer.count("G") + spacer.count("C")) / len(spacer)
        eff_score = score_guide(spacer, pbs_tm, rt_template)

        # For PE3, find nicking sgRNA
        nick_sgrna = None
        nick_pam = None
        if editor in ["PE3", "PE3b", "PE-max", "PEmax-MLH1dn"]:
            # Find nicking sgRNA on opposite strand
            opposite_seq = reverse_complement(target_seq)
            nick_sites = find_pams(opposite_seq, pam, spacer_len)
            for nick_start, _, nick_spacer in nick_sites:
                nick_dist = abs(nick_start - spacer_start)
                if nick_range[0] <= nick_dist <= nick_range[1]:
                    nick_sgrna = nick_spacer
                    nick_pam = "NGG"
                    break

        design = PeguideDesign(
            spacer=spacer,
            strand=strand,
            pam=pam,
            pbs=pbs,
            rt_template=rt_template,
            pbs_tm=pbs_tm,
            edit_distance_from_spacer=dist_from_edit,
            gc_content=gc,
            efficiency_score=eff_score,
            nick_sgrna=nick_sgrna,
            nick_pam=nick_pam,
        )
        guides.append(design)

    # Sort by efficiency score and return top N
    guides.sort(key=lambda g: g.efficiency_score, reverse=True)
    return guides[:top_n]


def save_results(guides: List[PeguideDesign], outdir: str):
    """Save guides to TSV and generate plots."""
    os.makedirs(outdir, exist_ok=True)

    # Save TSV
    records = []
    for i, g in enumerate(guides, 1):
        records.append({
            "rank": i,
            "spacer": g.spacer,
            "strand": g.strand,
            "pam": g.pam,
            "pbs": g.pbs,
            "rt_template": g.rt_template,
            "pbs_tm": round(g.pbs_tm, 2),
            "edit_distance": g.edit_distance_from_spacer,
            "gc_content": round(g.gc_content, 3),
            "efficiency_score": round(g.efficiency_score, 3),
            "nick_sgrna": g.nick_sgrna or "N/A",
            "nick_pam": g.nick_pam or "N/A",
        })

    df = pd.DataFrame(records)
    df.to_csv(os.path.join(outdir, "pegguides.tsv"), sep="\t", index=False)
    print(f"Saved {len(guides)} guides to {outdir}/pegguides.tsv")

    # Plot efficiency scores
    if guides:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # Efficiency score distribution
        scores = [g.efficiency_score for g in guides]
        axes[0, 0].hist(scores, bins=10, color="steelblue", edgecolor="black")
        axes[0, 0].set_xlabel("Efficiency Score")
        axes[0, 0].set_ylabel("Count")
        axes[0, 0].set_title("Guide Efficiency Score Distribution")

        # PBS Tm distribution
        tms = [g.pbs_tm for g in guides]
        axes[0, 1].hist(tms, bins=10, color="coral", edgecolor="black")
        axes[0, 1].axvline(33.0, color="red", linestyle="--", label="Optimal (33°C)")
        axes[0, 1].set_xlabel("PBS Tm (°C)")
        axes[0, 1].set_ylabel("Count")
        axes[0, 1].set_title("PBS Tm Distribution")
        axes[0, 1].legend()

        # GC content distribution
        gcs = [g.gc_content for g in guides]
        axes[1, 0].hist(gcs, bins=10, color="lightgreen", edgecolor="black")
        axes[1, 0].axvspan(0.4, 0.7, alpha=0.2, color="green", label="Optimal")
        axes[1, 0].set_xlabel("GC Content")
        axes[1, 0].set_ylabel("Count")
        axes[1, 0].set_title("GC Content Distribution")
        axes[1, 0].legend()

        # Edit distance distribution
        dists = [g.edit_distance_from_spacer for g in guides]
        axes[1, 1].hist(dists, bins=10, color="plum", edgecolor="black")
        axes[1, 1].set_xlabel("Distance from Edit Site (bp)")
        axes[1, 1].set_ylabel("Count")
        axes[1, 1].set_title("Spacer-to-Edit Distance")

        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "pegguide_distributions.png"), dpi=300, bbox_inches="tight")
        print(f"Saved plot to {outdir}/pegguide_distributions.png")
        plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Design pegRNAs and nicking sgRNAs for prime editing."
    )
    parser.add_argument("--gene", type=str, help="Gene symbol (e.g., TP53)")
    parser.add_argument("--species", type=str, default="human", choices=["human", "mouse"])
    parser.add_argument("--edit-type", type=str, required=True, choices=["snv", "insertion", "deletion"])
    parser.add_argument("--target-sequence", type=str, help="Target sequence (20bp+)")
    parser.add_argument("--wt-seq", type=str, help="Wild-type sequence around edit site")
    parser.add_argument("--edit-seq", type=str, help="Desired edited sequence (same length as wt-seq for SNV)")
    parser.add_argument("--position", type=int, help="Genomic position of edit (1-indexed)")
    parser.add_argument("--editor", type=str, default="PE3",
                       choices=["PE2", "PE3", "PE3b", "PE-max", "PEmax-MLH1dn"])
    parser.add_argument("--pam", type=str, default="NGG", help="PAM sequence")
    parser.add_argument("--spacer-len", type=int, default=20)
    parser.add_argument("--pbs-length-range", type=str, default="8,17", help="min,max PBS length")
    parser.add_argument("--rt-length-range", type=str, default="10,40", help="min,max RT length")
    parser.add_argument("--nick-distance-range", type=str, default="40,100", help="min,max nick distance (PE3)")
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--outdir", type=str, default="pegguide_results")
    parser.add_argument("--genome", type=str, default="hg38", choices=["hg38", "hg19", "mm10"])

    args = parser.parse_args()

    init_style(
        font_family=getattr(args, 'font_family', None),
        font_size=getattr(args, 'font_size', None),
    )

    # Resolve target sequence
    target_seq = args.target_sequence
    edit_pos = args.position or 0

    if not target_seq:
        if args.gene:
            target_seq, _ = fetch_sequence_ensembl(args.gene, args.species)
            if not target_seq:
                print("Error: Could not fetch sequence from Ensembl", file=sys.stderr)
                sys.exit(1)
            edit_pos = len(target_seq) // 2  # Default to middle
        else:
            print("Error: Provide --target-sequence or --gene", file=sys.stderr)
            sys.exit(1)

    if not args.wt_seq or not args.edit_seq:
        print("Error: Provide --wt-seq and --edit-seq", file=sys.stderr)
        sys.exit(1)

    target_seq = target_seq.upper()

    # Parse ranges
    pbs_min, pbs_max = map(int, args.pbs_length_range.split(","))
    rt_min, rt_max = map(int, args.rt_length_range.split(","))
    nick_min, nick_max = map(int, args.nick_distance_range.split(","))

    # Design pegguides
    guides = design_pegguides(
        target_seq=target_seq,
        edit_type=args.edit_type,
        wt_seq=args.wt_seq,
        edit_seq=args.edit_seq,
        edit_pos=edit_pos,
        pam=args.pam,
        spacer_len=args.spacer_len,
        pbs_range=(pbs_min, pbs_max),
        rt_range=(rt_min, rt_max),
        nick_range=(nick_min, nick_max),
        editor=args.editor,
        top_n=args.top_n,
    )

    if not guides:
        print("No candidate guides found. Try relaxing filter parameters.", file=sys.stderr)
        sys.exit(1)

    print(f"Designed {len(guides)} pegRNA candidates")
    save_results(guides, args.outdir)


if __name__ == "__main__":
    main()
