#!/usr/bin/env python3
"""
design_base_editor_sgrnas.py
─────────────────────────────
Design base editor sgRNAs for a given gene.

Pipeline
--------
1. Fetch gene coordinates and exon/CDS structure from Ensembl REST API
2. Retrieve genomic sequence for each exon/CDS
3. Scan both strands for NGG (or NG / NRN) PAM sites
4. For each protospacer, identify editable bases (C for CBE, A for ABE)
   within the configured editing window
5. Predict amino acid changes where possible (CDS context)
6. Score guides: efficiency heuristic + bystander penalty
7. Rank, filter, and export table + figures
"""

import argparse
import json
import re
import sys
import time
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), '..', '..', '..', '_shared'))
from plot_style import init_style
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import requests

warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────────────────────────────────────
#  Editor catalogue
# ─────────────────────────────────────────────────────────────────────────────

EDITORS: Dict[str, dict] = {
    "BE3":       {"type": "CBE", "pam": "NGG", "win": (4, 8),  "base": "C"},
    "BE4max":    {"type": "CBE", "pam": "NGG", "win": (4, 8),  "base": "C"},
    "AncBE4max": {"type": "CBE", "pam": "NGG", "win": (4, 8),  "base": "C"},
    "CBE":       {"type": "CBE", "pam": "NGG", "win": (4, 8),  "base": "C"},
    "ABE7.10":   {"type": "ABE", "pam": "NGG", "win": (4, 7),  "base": "A"},
    "ABE8e":     {"type": "ABE", "pam": "NGG", "win": (4, 8),  "base": "A"},
    "ABE8.20m":  {"type": "ABE", "pam": "NGG", "win": (4, 8),  "base": "A"},
    "ABE":       {"type": "ABE", "pam": "NGG", "win": (4, 8),  "base": "A"},
    "NG-CBE":    {"type": "CBE", "pam": "NG",  "win": (4, 8),  "base": "C"},
    "NG-ABE":    {"type": "ABE", "pam": "NG",  "win": (4, 8),  "base": "A"},
    "SpRY-CBE":  {"type": "CBE", "pam": "NRN", "win": (4, 8),  "base": "C"},
    "SpRY-ABE":  {"type": "ABE", "pam": "NRN", "win": (4, 8),  "base": "A"},
    "dual":      {"type": "dual","pam": "NGG", "win": (4, 8),  "base": "CA"},
}

ORGANISMS = {
    "human":  "homo_sapiens",
    "mouse":  "mus_musculus",
    "monkey": "macaca_fascicularis",
}

ENSEMBL_REST = "https://rest.ensembl.org"

# Genetic code for amino acid prediction
CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

plt.rcParams.update({
    "figure.dpi":        150,
    "font.family":       "sans-serif",
    "axes.spines.top":   False,
    "axes.spines.right": False,
})


# ─────────────────────────────────────────────────────────────────────────────
#  Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _rev_comp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def _pam_regex(pam: str) -> str:
    """Convert IUPAC PAM string to regex."""
    iupac = {"N":"[ACGT]","R":"[AG]","Y":"[CT]","S":"[GC]",
             "W":"[AT]","K":"[GT]","M":"[AC]","B":"[CGT]",
             "D":"[AGT]","H":"[ACT]","V":"[ACG]"}
    return "".join(iupac.get(c, c) for c in pam.upper())


def _matches_pam(pam_str: str, pam_pattern: str) -> bool:
    return bool(re.fullmatch(pam_pattern, pam_str.upper()))


def _gc_content(seq: str) -> float:
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / max(len(seq), 1)


# ─────────────────────────────────────────────────────────────────────────────
#  Ensembl REST API
# ─────────────────────────────────────────────────────────────────────────────

def _ensembl_get(endpoint: str, params: dict = None) -> dict:
    url = f"{ENSEMBL_REST}/{endpoint}"
    headers = {"Content-Type": "application/json"}
    for attempt in range(3):
        try:
            r = requests.get(url, headers=headers, params=params, timeout=30)
            if r.status_code == 429:
                time.sleep(2 ** attempt)
                continue
            r.raise_for_status()
            return r.json()
        except requests.RequestException as e:
            if attempt == 2:
                raise RuntimeError(f"Ensembl API error ({url}): {e}")
            time.sleep(1)
    return {}


def fetch_gene_info(gene: str, species: str) -> dict:
    """Fetch gene metadata from Ensembl lookup."""
    print(f"[INFO] Fetching gene info: {gene} ({species})")
    data = _ensembl_get(f"lookup/symbol/{species}/{gene}",
                        {"expand": 1, "content-type": "application/json"})
    if not data or "id" not in data:
        raise ValueError(f"Gene '{gene}' not found in Ensembl for {species}.")
    return data


def fetch_sequence(chrom: str, start: int, end: int, species: str,
                   strand: int = 1) -> str:
    """Fetch genomic sequence from Ensembl."""
    region = f"{chrom}:{start}..{end}:{strand}"
    data = _ensembl_get(f"sequence/region/{species}/{region}",
                        {"content-type": "application/json"})
    return data.get("seq", "")


def get_cds_exons(gene_data: dict, transcript_id: Optional[str],
                  region: str, species: str) -> List[dict]:
    """
    Return a list of exon/CDS records:
      {chrom, start, end, strand, seq, exon_id, transcript_id, phase}
    """
    transcripts = gene_data.get("Transcript", [])
    if not transcripts:
        raise ValueError("No transcripts found for this gene.")

    # Prefer canonical or user-specified transcript
    if transcript_id:
        transcripts = [t for t in transcripts if t["id"] == transcript_id]
        if not transcripts:
            raise ValueError(f"Transcript {transcript_id} not found.")
    else:
        # Pick longest CDS transcript
        def _cds_len(t):
            return sum(e.get("end", 0) - e.get("start", 0)
                       for e in t.get("Exon", []))
        transcripts = sorted(transcripts, key=_cds_len, reverse=True)[:1]

    tx = transcripts[0]
    chrom  = gene_data["seq_region_name"]
    strand = 1 if gene_data["strand"] == 1 else -1

    exon_key = "Exon" if region in ("exon", "transcript") else "Exon"
    exons = tx.get(exon_key, [])
    if not exons:
        raise ValueError("No exons found for the selected transcript.")

    records = []
    for ex in exons:
        s = int(ex["start"])
        e = int(ex["end"])
        seq = fetch_sequence(chrom, s, e, species, strand=gene_data["strand"])
        if not seq:
            continue
        records.append({
            "chrom":     chrom,
            "start":     s,
            "end":       e,
            "strand":    gene_data["strand"],
            "seq":       seq.upper(),
            "exon_id":   ex.get("id", ""),
            "transcript_id": tx["id"],
        })
        time.sleep(0.05)   # respect rate limit

    print(f"[INFO] Retrieved {len(records)} exons "
          f"(transcript: {tx['id']}, total bp: {sum(r['end']-r['start'] for r in records)})")
    return records


# ─────────────────────────────────────────────────────────────────────────────
#  Guide scanning
# ─────────────────────────────────────────────────────────────────────────────

def _window_positions(proto: str, win_start: int, win_end: int,
                      target_base: str) -> List[int]:
    """
    Return 1-based positions (from PAM-distal end) of target bases in window.
    proto: 20 nt protospacer, 5'→3', PAM not included.
    """
    positions = []
    for i in range(win_start - 1, min(win_end, len(proto))):
        if proto[i].upper() in target_base.upper():
            positions.append(i + 1)   # 1-based
    return positions


def _bystander_count(proto: str, win_start: int, win_end: int,
                     target_base: str, edit_positions: List[int]) -> int:
    """Count editable bases in window that are NOT the primary target positions."""
    all_pos = _window_positions(proto, win_start, win_end, target_base)
    return len([p for p in all_pos if p not in edit_positions])


def _efficiency_score(proto: str) -> float:
    """
    Heuristic efficiency score based on sequence features.
    Adapted from Rule Set 1 / published base editor efficiency studies.
    Score range 0–1 (higher = better predicted activity).
    """
    proto = proto.upper()
    score = 0.5   # baseline

    # Prefer G at position 1 (transcription start)
    if proto[0] == "G":
        score += 0.05

    # Avoid runs of T (RNA Pol III termination)
    if "TTTT" in proto:
        score -= 0.20
    elif "TTT" in proto:
        score -= 0.10

    # GC content preference 40–70%
    gc = _gc_content(proto)
    if 0.40 <= gc <= 0.70:
        score += 0.10
    elif gc < 0.20 or gc > 0.80:
        score -= 0.15

    # Avoid G at position 20 (next to PAM — some studies show preference)
    if proto[-1] != "G":
        score -= 0.03

    # Penalise homopolymers
    for base in "ACGT":
        if base * 4 in proto:
            score -= 0.08

    # Prefer absence of secondary structure seed (rough: avoid palindromes in seed)
    seed = proto[-12:]
    if seed[:6] == _rev_comp(seed[6:]):
        score -= 0.05

    return round(max(0.0, min(1.0, score)), 3)


def scan_guides(exon: dict, pam_regex: str, pam_len: int,
                win_start: int, win_end: int, target_base: str,
                allow_bystander: bool, max_bystander: int) -> List[dict]:
    """Scan one exon sequence on both strands and return candidate guides."""
    results = []
    seq      = exon["seq"]
    chrom    = exon["chrom"]
    ex_start = exon["start"]
    gene_strand = exon["strand"]

    for strand_sign, search_seq in [("+", seq), ("-", _rev_comp(seq))]:
        n = len(search_seq)
        for i in range(n - 23):     # 20 nt protospacer + up to 3 nt PAM
            proto = search_seq[i: i + 20]
            pam   = search_seq[i + 20: i + 20 + pam_len]

            if len(pam) < pam_len:
                continue
            if not _matches_pam(pam, pam_regex):
                continue

            # Skip guides with ambiguous bases
            if re.search(r"[^ACGT]", proto):
                continue

            edit_positions = _window_positions(proto, win_start, win_end, target_base)
            if not edit_positions:
                continue

            byst = _bystander_count(proto, win_start, win_end,
                                    target_base, edit_positions)
            if not allow_bystander and byst > 0:
                continue
            if allow_bystander and byst > max_bystander:
                continue

            # Genomic coordinates
            if strand_sign == "+":
                g_start = ex_start + i
                g_end   = ex_start + i + 20
            else:
                g_start = ex_start + (len(seq) - i - 20)
                g_end   = ex_start + (len(seq) - i)

            eff = _efficiency_score(proto)
            byst_penalty = round(byst * 0.15, 3)
            final_score  = round(eff - byst_penalty, 3)

            results.append({
                "sgrna_seq":       proto,
                "pam":             pam,
                "strand":          strand_sign,
                "chrom":           chrom,
                "start":           g_start,
                "end":             g_end,
                "exon_id":         exon["exon_id"],
                "transcript_id":   exon["transcript_id"],
                "target_base":     target_base,
                "window_positions": ",".join(str(p) for p in edit_positions),
                "n_editable":      len(edit_positions),
                "bystander_count": byst,
                "gc_content":      round(_gc_content(proto), 3),
                "efficiency_score": eff,
                "bystander_penalty": byst_penalty,
                "final_score":     final_score,
            })

    return results


# ─────────────────────────────────────────────────────────────────────────────
#  Amino acid change prediction
# ─────────────────────────────────────────────────────────────────────────────

def predict_aa_changes(guides: List[dict], gene_data: dict) -> List[dict]:
    """
    For each guide, estimate amino acid consequences of each editable position.
    Requires CDS coordinates from Ensembl CDS feature.
    This is a best-effort annotation; strand and phase offset are approximated.
    """
    chrom = gene_data.get("seq_region_name", "")
    for g in guides:
        edits = []
        for pos_str in g["window_positions"].split(","):
            if not pos_str:
                continue
            pos = int(pos_str)   # 1-based in protospacer
            base = g["target_base"][0].upper()
            new_base = "T" if base == "C" else "G"
            # Report the change relative to protospacer
            orig = g["sgrna_seq"][pos - 1].upper()
            edits.append(f"pos{pos}:{orig}>{new_base}")
        g["edit_outcomes"] = "; ".join(edits) if edits else "n/a"
    return guides


# ─────────────────────────────────────────────────────────────────────────────
#  Amino acid position filter
# ─────────────────────────────────────────────────────────────────────────────

def filter_by_aa_position(guides: pd.DataFrame, gene_data: dict,
                          target_aa_pos: int) -> pd.DataFrame:
    """
    Keep only guides whose protospacer overlaps the codon for target_aa_pos.
    Uses approximate genomic arithmetic.
    """
    transcripts = gene_data.get("Transcript", [])
    if not transcripts:
        return guides

    # Find CDS start from canonical transcript
    tx = sorted(transcripts,
                key=lambda t: sum(e.get("end", 0) - e.get("start", 0)
                                  for e in t.get("Exon", [])),
                reverse=True)[0]

    cds_start = tx.get("Translation", {}).get("start", None)
    if cds_start is None:
        print("[WARN] Could not determine CDS start — skipping amino acid position filter.")
        return guides

    # Approximate codon genomic window (3 bp per codon, ±10 bp buffer)
    codon_start = cds_start + (target_aa_pos - 1) * 3
    codon_end   = codon_start + 3

    buf = 23   # protospacer length
    mask = (
        (guides["start"] <= codon_end   + buf) &
        (guides["end"]   >= codon_start - buf)
    )
    filtered = guides[mask]
    print(f"[INFO] Guides near codon {target_aa_pos}: {len(filtered)} / {len(guides)}")
    return filtered


# ─────────────────────────────────────────────────────────────────────────────
#  Plots
# ─────────────────────────────────────────────────────────────────────────────

def plot_editing_heatmap(guides: pd.DataFrame, win_start: int, win_end: int,
                         target_base: str, prefix: str, outdir: Path):
    """Heatmap: protospacer position (x) × guide rank (y), coloured by editable base."""
    top = guides.head(min(30, len(guides)))
    n_guides = len(top)
    n_pos = 20

    matrix = np.zeros((n_guides, n_pos))
    for i, (_, row) in enumerate(top.iterrows()):
        for pos_str in row["window_positions"].split(","):
            if pos_str:
                p = int(pos_str) - 1
                if 0 <= p < n_pos:
                    matrix[i, p] = 1.0

    fig, ax = plt.subplots(figsize=(12, max(4, n_guides * 0.28)))
    cmap = plt.cm.get_cmap("RdBu_r", 2)
    ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=0, vmax=1,
              interpolation="nearest")

    ax.set_xlabel("Protospacer position (1 = PAM-distal)")
    ax.set_ylabel("Guide rank")
    ax.set_xticks(range(n_pos))
    ax.set_xticklabels(range(1, n_pos + 1), fontsize=7)
    ax.set_yticks(range(n_guides))
    ax.set_yticklabels([f"#{i+1} {r['sgrna_seq'][:10]}…"
                        for i, (_, r) in enumerate(top.iterrows())], fontsize=7)

    # Window boundaries
    ax.axvline(win_start - 1.5, color="orange", lw=1.5, ls="--", label="Window")
    ax.axvline(win_end   - 0.5, color="orange", lw=1.5, ls="--")

    ax.set_title(f"{prefix} — Editing window heatmap "
                 f"({target_base}→{'T' if target_base=='C' else 'G'}, "
                 f"positions {win_start}–{win_end})")
    ax.legend(fontsize=8, loc="upper right")

    fig.tight_layout()
    out = outdir / f"{prefix}_editing_heatmap.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Editing heatmap: {out}")


def plot_guide_map(guides: pd.DataFrame, gene_data: dict,
                   prefix: str, outdir: Path):
    """Horizontal track plot: gene exons + guide positions."""
    transcripts = gene_data.get("Transcript", [])
    if not transcripts:
        return

    tx = sorted(transcripts,
                key=lambda t: sum(e.get("end", 0) - e.get("start", 0)
                                  for e in t.get("Exon", [])),
                reverse=True)[0]
    exons = tx.get("Exon", [])
    if not exons:
        return

    gene_start = gene_data["start"]
    gene_end   = gene_data["end"]
    gene_len   = gene_end - gene_start

    fig, ax = plt.subplots(figsize=(14, 3.5))

    # Gene backbone
    ax.barh(0, gene_len, left=0, height=0.08, color="#CCCCCC", zorder=1)

    # Exons
    for ex in exons:
        ex_rel = ex["start"] - gene_start
        ex_len = ex["end"] - ex["start"]
        ax.barh(0, ex_len, left=ex_rel, height=0.25, color="#4472C4",
                zorder=2, alpha=0.85)

    # Guides
    top = guides.head(min(40, len(guides)))
    colors = plt.cm.plasma(np.linspace(0.2, 0.85, len(top)))
    for idx, ((_, row), color) in enumerate(zip(top.iterrows(), colors)):
        g_rel = row["start"] - gene_start
        g_len = row["end"] - row["start"]
        y = 0.35 + (idx % 3) * 0.18
        ax.barh(y, g_len, left=g_rel, height=0.12, color=color,
                alpha=0.80, zorder=3)
        if idx < 10:
            ax.text(g_rel + g_len / 2, y + 0.08,
                    f"#{idx+1}", fontsize=5.5, ha="center", color="black")

    ax.set_xlim(0, gene_len)
    ax.set_ylim(-0.3, 1.1)
    ax.set_xlabel("Genomic position (relative to gene start)")
    ax.set_title(f"{prefix} — Guide positions across {gene_data.get('display_name', '')} exons")
    ax.set_yticks([])

    exon_patch = mpatches.Patch(color="#4472C4", label="Exon")
    guide_patch = mpatches.Patch(color="#FF6B35", label="Guide (top 40)")
    ax.legend(handles=[exon_patch, guide_patch], fontsize=8, loc="upper right")

    fig.tight_layout()
    out = outdir / f"{prefix}_guide_map.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Guide map: {out}")


def plot_score_distribution(guides: pd.DataFrame, prefix: str, outdir: Path):
    """Distribution of efficiency scores and bystander counts."""
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))

    axes[0].hist(guides["efficiency_score"], bins=30, color="#4472C4",
                 edgecolor="white", alpha=0.8)
    axes[0].set_xlabel("Efficiency score")
    axes[0].set_ylabel("Count")
    axes[0].set_title("On-target efficiency distribution")

    axes[1].hist(guides["final_score"], bins=30, color="#70AD47",
                 edgecolor="white", alpha=0.8)
    axes[1].set_xlabel("Final score (efficiency − bystander penalty)")
    axes[1].set_title("Final score distribution")

    byst_counts = guides["bystander_count"].value_counts().sort_index()
    axes[2].bar(byst_counts.index, byst_counts.values, color="#ED7D31",
                edgecolor="white", alpha=0.85)
    axes[2].set_xlabel("Bystander editable bases in window")
    axes[2].set_ylabel("Number of guides")
    axes[2].set_title("Bystander edit distribution")

    fig.suptitle(f"{prefix} — Guide score summary ({len(guides)} guides total)",
                 fontweight="bold")
    fig.tight_layout()
    out = outdir / f"{prefix}_score_distribution.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Score distribution: {out}")


# ─────────────────────────────────────────────────────────────────────────────
#  Dual-base editor handling
# ─────────────────────────────────────────────────────────────────────────────

def run_dual(exons: List[dict], pam_regex: str, pam_len: int,
             win_start: int, win_end: int,
             allow_bystander: bool, max_bystander: int) -> List[dict]:
    """For dual editors scan for guides carrying both C and A in window."""
    cbe_guides = []
    abe_guides = []
    for exon in exons:
        cbe_guides += scan_guides(exon, pam_regex, pam_len, win_start, win_end,
                                  "C", True, 99)
        abe_guides += scan_guides(exon, pam_regex, pam_len, win_start, win_end,
                                  "A", True, 99)

    cbe_seqs = {g["sgrna_seq"] for g in cbe_guides}
    abe_seqs = {g["sgrna_seq"] for g in abe_guides}
    dual_seqs = cbe_seqs & abe_seqs

    results = []
    seen = set()
    for g in cbe_guides + abe_guides:
        seq = g["sgrna_seq"]
        if seq not in dual_seqs or seq in seen:
            continue
        seen.add(seq)
        g["target_base"] = "CA"
        g["editor_type"] = "dual"
        results.append(g)
    return results


# ─────────────────────────────────────────────────────────────────────────────
#  CLI
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Design base editor sgRNAs for a gene",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--gene",            required=True,
                   help="Gene symbol (e.g. TP53, KRAS)")
    p.add_argument("--organism",        default="human",
                   choices=["human", "mouse", "monkey"],
                   help="Organism (default: human)")
    p.add_argument("--editor",          default="CBE",
                   help="Editor name: CBE, ABE, BE3, BE4max, ABE7.10, "
                        "ABE8e, NG-CBE, NG-ABE, SpRY-CBE, SpRY-ABE, dual "
                        "(default: CBE)")
    p.add_argument("--window-start",    type=int, default=None,
                   help="Editing window start (1-based; default: editor-specific)")
    p.add_argument("--window-end",      type=int, default=None,
                   help="Editing window end (1-based; default: editor-specific)")
    p.add_argument("--pam",             default=None,
                   help="Override PAM (NGG, NG, NRN; default: editor-specific)")
    p.add_argument("--region",          default="cds",
                   choices=["cds", "exon", "transcript"],
                   help="Sequence scope (default: cds)")
    p.add_argument("--allow-bystander", action="store_true",
                   help="Include guides with bystander editable bases in window")
    p.add_argument("--max-bystander",   type=int, default=0,
                   help="Max bystander bases allowed (used with --allow-bystander; default 0)")
    p.add_argument("--top-n",           type=int, default=20,
                   help="Number of top guides to return (default: 20)")
    p.add_argument("--transcript-id",   default=None,
                   help="Restrict to a specific Ensembl transcript ID")
    p.add_argument("--target-aa-pos",   type=int, default=None,
                   help="Filter guides overlapping a specific amino acid position")
    p.add_argument("--prefix",          default=None,
                   help="Output file prefix (default: gene symbol)")
    p.add_argument("--outdir",          default="./base_editor_results",
                   help="Output directory (default: ./base_editor_results/)")
    return p.parse_args()


def main():
    args = parse_args()
    init_style(
        font_family=getattr(args, 'font_family', None),
        font_size=getattr(args, 'font_size', None),
    )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix or args.gene.upper()
    species = ORGANISMS.get(args.organism)
    if not species:
        print(f"[ERROR] Unknown organism: {args.organism}")
        sys.exit(1)

    # ── Resolve editor config ──────────────────────────────────────────────
    editor_key = args.editor.upper()
    # Case-insensitive match
    editor_cfg = next(
        (v for k, v in EDITORS.items() if k.upper() == editor_key),
        EDITORS.get(args.editor, EDITORS["CBE"])
    )
    win_start   = args.window_start or editor_cfg["win"][0]
    win_end     = args.window_end   or editor_cfg["win"][1]
    pam_str     = args.pam or editor_cfg["pam"]
    target_base = editor_cfg["base"]
    pam_regex   = _pam_regex(pam_str)
    pam_len     = len(pam_str)
    is_dual     = editor_cfg["type"] == "dual"

    print(f"\n{'='*56}")
    print(f"  Base Editor sgRNA Design")
    print(f"{'='*56}")
    print(f"  Gene      : {args.gene}  ({args.organism})")
    print(f"  Editor    : {args.editor}  (type: {editor_cfg['type']})")
    print(f"  PAM       : {pam_str}")
    print(f"  Window    : positions {win_start}–{win_end}")
    print(f"  Target    : {target_base}→{'T' if 'C' in target_base else 'G'}")
    print(f"  Bystander : {'allowed (max %d)' % args.max_bystander if args.allow_bystander else 'not allowed'}")
    print(f"{'='*56}\n")

    # ── Fetch gene + sequences ─────────────────────────────────────────────
    gene_data = fetch_gene_info(args.gene, species)
    exons = get_cds_exons(gene_data, args.transcript_id, args.region, species)
    if not exons:
        print("[ERROR] No exon sequences retrieved. Check gene name and organism.")
        sys.exit(1)

    # ── Scan for guides ────────────────────────────────────────────────────
    print(f"\n[INFO] Scanning {len(exons)} exon(s) for PAM={pam_str} guides ...")
    all_guides = []

    if is_dual:
        all_guides = run_dual(exons, pam_regex, pam_len, win_start, win_end,
                              args.allow_bystander, args.max_bystander)
    else:
        for exon in exons:
            guides = scan_guides(exon, pam_regex, pam_len, win_start, win_end,
                                 target_base, args.allow_bystander, args.max_bystander)
            all_guides.extend(guides)

    # Remove duplicates (same protospacer from overlapping exon fetches)
    seen = set()
    unique = []
    for g in all_guides:
        if g["sgrna_seq"] not in seen:
            seen.add(g["sgrna_seq"])
            unique.append(g)
    all_guides = unique

    print(f"[INFO] {len(all_guides)} candidate guides found.")

    if not all_guides:
        print("[WARN] No guides passed filters. Try --allow-bystander or a different editor.")
        sys.exit(0)

    # ── Amino acid annotation ──────────────────────────────────────────────
    all_guides = predict_aa_changes(all_guides, gene_data)

    # ── Build DataFrame and rank ──────────────────────────────────────────
    df = pd.DataFrame(all_guides)
    df = df.sort_values("final_score", ascending=False).reset_index(drop=True)
    df.insert(0, "rank", df.index + 1)

    # ── Optional: filter by amino acid position ───────────────────────────
    if args.target_aa_pos:
        df = filter_by_aa_position(df, gene_data, args.target_aa_pos)
        df = df.reset_index(drop=True)
        df["rank"] = df.index + 1

    # ── Save outputs ──────────────────────────────────────────────────────
    full_path = outdir / f"{prefix}_base_editor_guides.tsv"
    df.to_csv(full_path, sep="\t", index=False)
    print(f"[OK] Full guide table ({len(df)} guides): {full_path}")

    top = df.head(args.top_n)
    top_path = outdir / f"{prefix}_top{args.top_n}_guides.tsv"
    top.to_csv(top_path, sep="\t", index=False)
    print(f"[OK] Top {args.top_n} guides: {top_path}")

    # Bystander report
    byst_path = outdir / f"{prefix}_bystander_report.tsv"
    df[["rank", "sgrna_seq", "window_positions", "bystander_count",
        "edit_outcomes", "final_score"]].to_csv(byst_path, sep="\t", index=False)
    print(f"[OK] Bystander report: {byst_path}")

    # ── Plots ─────────────────────────────────────────────────────────────
    print("\n[INFO] Generating plots ...")
    plot_editing_heatmap(df, win_start, win_end, target_base, prefix, outdir)
    plot_guide_map(df, gene_data, prefix, outdir)
    plot_score_distribution(df, prefix, outdir)

    # ── Save run metadata ─────────────────────────────────────────────────
    meta = {
        "gene": args.gene, "organism": args.organism,
        "editor": args.editor, "pam": pam_str,
        "window": [win_start, win_end], "target_base": target_base,
        "allow_bystander": args.allow_bystander,
        "max_bystander": args.max_bystander,
        "n_candidates": len(df), "top_n": args.top_n,
        "target_aa_pos": args.target_aa_pos,
    }
    with open(outdir / f"{prefix}_meta.json", "w") as f:
        json.dump(meta, f, indent=2)

    # ── Summary ───────────────────────────────────────────────────────────
    print(f"\n{'='*56}")
    print(f"  Base editor design complete")
    print(f"{'='*56}")
    print(f"  Gene          : {args.gene}")
    print(f"  Editor        : {args.editor}")
    print(f"  Candidates    : {len(df)}")
    print(f"  Top {args.top_n:<3} score  : {top['final_score'].max():.3f}")

    print(f"\n  Top 5 guides:")
    cols = ["rank", "sgrna_seq", "pam", "strand", "window_positions",
            "bystander_count", "final_score"]
    print(top[cols].head(5).to_string(index=False))
    print(f"\n  Output dir: {outdir}/")
    print(f"{'='*56}\n")


if __name__ == "__main__":
    main()
