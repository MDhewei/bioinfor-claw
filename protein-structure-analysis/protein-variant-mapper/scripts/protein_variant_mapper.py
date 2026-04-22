#!/usr/bin/env python3
"""
Protein Variant / Mutation Mapper

Maps a list of amino-acid substitution variants onto a protein 3D structure and produces:

  - Variant summary table: position, wild-type residue, mutant residue,
    SASA (solvent accessibility), secondary structure, B-factor,
    proximity to known pockets (if pocket_residues.tsv provided)
  - Interactive HTML 3D viewer: each variant coloured by impact class
    (pathogenic / benign / uncertain / custom), overlaid on the full protein
  - Linear protein map (PNG + PDF): variant positions along the sequence

Variant input format  (--variants):
  Comma-separated substitution strings, e.g.  "A123V,G45S,R280*,L12P"
  - First character = wild-type amino acid (1-letter)
  - Middle digits   = sequence position (must match PDB resseq)
  - Last character  = mutant amino acid (1-letter, or * for stop / frameshift)

Optionally fetch known UniProt natural variants to auto-classify.

Dependencies:
  biopython>=1.81, matplotlib>=3.7, numpy>=1.24, pandas>=2.0,
  requests>=2.28, py3Dmol>=2.0
"""

import argparse
import os
import re
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from Bio import PDB
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

# =========================================================
# Constants
# =========================================================

RCSB_FILE     = "https://files.rcsb.org/download"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api"
ALPHAFOLD_FILE = "https://alphafold.ebi.ac.uk/files"
UNIPROT_API   = "https://rest.uniprot.org/uniprotkb"
HEADERS = {"User-Agent": "protein-variant-mapper/0.1 (bioinfor-claw)"}

# Impact class colours
IMPACT_COLORS = {
    "pathogenic":  "#D62728",   # red
    "likely_pathogenic": "#FF7F0E",  # orange
    "benign":      "#2CA02C",   # green
    "likely_benign": "#98DF8A", # light green
    "uncertain":   "#9467BD",   # purple
    "custom":      "#1F77B4",   # blue (user-supplied, no classification)
    "stop":        "#000000",   # black — stop codon / frameshift
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O", "UNK": "X",
}


# =========================================================
# Utilities
# =========================================================

def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def save_fig(fig: plt.Figure, base: str) -> None:
    fig.savefig(f"{base}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{base}.pdf", bbox_inches="tight")
    plt.close(fig)


def _get(url: str, params: Optional[dict] = None, timeout: int = 60) -> requests.Response:
    r = requests.get(url, params=params, headers=HEADERS, timeout=timeout)
    r.raise_for_status()
    return r


# =========================================================
# Variant parsing
# =========================================================

VARIANT_RE = re.compile(r"^([A-Za-z])(\d+)([A-Za-z*=])$")


def parse_variants(variant_str: str) -> List[Dict]:
    """
    Parse a comma-separated variant string into a list of dicts.

    Each dict has keys: raw, wt_aa, pos, mut_aa, impact_class.
    """
    variants = []
    for token in variant_str.split(","):
        token = token.strip()
        if not token:
            continue
        m = VARIANT_RE.match(token)
        if not m:
            print(f"[WARN] Could not parse variant '{token}' — skipping. "
                  f"Expected format: A123V")
            continue
        wt, pos, mut = m.group(1).upper(), int(m.group(2)), m.group(3).upper()
        impact = "stop" if mut in ("*", "X") else "custom"
        variants.append({
            "raw":          token,
            "wt_aa":        wt,
            "pos":          pos,
            "mut_aa":       mut,
            "impact_class": impact,
        })
    return variants


def load_variants_file(path: str) -> List[Dict]:
    """
    Load variants from a file. Supports two formats:
      1. One variant per line (e.g. A123V)
      2. TSV with at minimum columns: variant, impact_class
    """
    with open(path) as f:
        first_line = f.readline().strip()

    # Detect TSV header
    if "\t" in first_line and not VARIANT_RE.match(first_line):
        df = pd.read_csv(path, sep="\t")
        results = []
        for _, row in df.iterrows():
            parsed = parse_variants(str(row["variant"]))
            if parsed:
                v = parsed[0]
                if "impact_class" in df.columns:
                    v["impact_class"] = str(row["impact_class"]).lower()
                results.append(v)
        return results

    # One variant per line
    with open(path) as f:
        lines = [l.strip() for l in f if l.strip()]
    return parse_variants(",".join(lines))


# =========================================================
# Structure fetching and parsing
# =========================================================

def fetch_pdb_structure(pdb_id: str, outdir: str) -> str:
    pdb_id = pdb_id.upper()
    out_path = os.path.join(outdir, f"{pdb_id}.pdb")
    if os.path.exists(out_path):
        return out_path
    r = _get(f"{RCSB_FILE}/{pdb_id}.pdb")
    with open(out_path, "w") as f:
        f.write(r.text)
    print(f"[INFO] Downloaded PDB: {pdb_id}")
    return out_path


def fetch_alphafold_structure(uniprot_acc: str, outdir: str) -> str:
    uniprot_acc = uniprot_acc.upper()
    out_path = os.path.join(outdir, f"AF_{uniprot_acc}.pdb")
    if os.path.exists(out_path):
        return out_path
    r = _get(f"{ALPHAFOLD_API}/prediction/{uniprot_acc}")
    entries = r.json()
    if not entries:
        raise ValueError(f"No AlphaFold entry for {uniprot_acc}")
    r2 = _get(entries[0]["pdbUrl"])
    with open(out_path, "w") as f:
        f.write(r2.text)
    print(f"[INFO] Downloaded AlphaFold: {uniprot_acc}")
    return out_path


def load_structure(pdb_path: str, struct_id: str = "STR") -> PDB.Structure.Structure:
    parser = PDBParser(QUIET=True)
    return parser.get_structure(struct_id, pdb_path)


def compute_sasa(structure: PDB.Structure.Structure, probe_radius: float = 1.4) -> Dict:
    """Return {(chain, resseq): sasa} mapping."""
    sr = ShrakeRupley(probe_radius=probe_radius)
    sr.compute(structure, level="R")
    sasa_map = {}
    for chain in structure[0].get_chains():
        for res in chain.get_residues():
            hetflag, resseq, _ = res.get_id()
            if hetflag.strip():
                continue
            sasa_map[(chain.get_id(), resseq)] = getattr(res, "sasa", 0.0)
    return sasa_map


def parse_secondary_structure(pdb_path: str) -> Dict:
    """
    Parse HELIX / SHEET records from PDB to assign secondary structure
    labels {(chain, resseq): 'H'|'E'|'C'}.
    Falls back to 'C' (coil) for all unassigned residues.
    """
    ss_map: Dict[Tuple, str] = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("HELIX"):
                chain = line[19]
                start = int(line[21:25])
                end   = int(line[33:37])
                for r in range(start, end + 1):
                    ss_map[(chain, r)] = "H"
            elif line.startswith("SHEET"):
                chain = line[21]
                start = int(line[22:26])
                end   = int(line[33:37])
                for r in range(start, end + 1):
                    ss_map[(chain, r)] = "E"
    return ss_map


# =========================================================
# UniProt natural variants
# =========================================================

def fetch_uniprot_variants(uniprot_acc: str) -> List[Dict]:
    """
    Fetch natural variants from UniProt and extract clinical significance
    to pre-classify variants as pathogenic/benign/uncertain.
    """
    url = f"{UNIPROT_API}/{uniprot_acc}.json"
    try:
        r = _get(url, timeout=30)
        data = r.json()
    except Exception as e:
        print(f"[WARN] Could not fetch UniProt variants: {e}")
        return []

    variants = []
    features = data.get("features", [])
    for feat in features:
        if feat.get("type") not in ("Natural variant", "Mutagenesis"):
            continue
        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end   = loc.get("end",   {}).get("value")
        if start is None:
            continue
        desc = feat.get("description", "")
        # Extract mut_aa from description like "A -> V"
        aa_match = re.search(r"([A-Z])\s*->\s*([A-Z*])", desc)
        wt_aa = aa_match.group(1) if aa_match else "?"
        mut_aa = aa_match.group(2) if aa_match else "?"

        # Classify based on keywords
        desc_lower = desc.lower()
        if any(w in desc_lower for w in ["pathogenic", "disease", "associated with"]):
            impact = "pathogenic"
        elif any(w in desc_lower for w in ["benign", "polymorphism", "no effect"]):
            impact = "benign"
        else:
            impact = "uncertain"

        variants.append({
            "pos":          int(start),
            "wt_aa":        wt_aa,
            "mut_aa":       mut_aa,
            "impact_class": impact,
            "description":  desc,
        })

    print(f"[INFO] Fetched {len(variants)} natural variants from UniProt.")
    return variants


# =========================================================
# Pocket proximity
# =========================================================

def load_pocket_residues(pocket_residues_path: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(pocket_residues_path, sep="\t")
        print(f"[INFO] Loaded {len(df)} pocket residues from {pocket_residues_path}")
        return df
    except Exception as e:
        print(f"[WARN] Could not load pocket residues: {e}")
        return pd.DataFrame()


def nearest_pocket_distance(
    pos: int,
    chain: str,
    structure: PDB.Structure.Structure,
    pocket_df: pd.DataFrame,
) -> float:
    """
    Return the minimum distance (Å) between this residue's Cα and
    any pocket residue's Cα. Returns np.inf if pocket_df is empty or
    residue has no Cα.
    """
    if pocket_df.empty:
        return np.inf

    try:
        res = structure[0][chain][(" ", pos, " ")]
    except KeyError:
        return np.inf

    if "CA" not in res:
        return np.inf

    coord = res["CA"].get_vector().get_array()
    pocket_coords = pocket_df[["x", "y", "z"]].values
    if len(pocket_coords) == 0:
        return np.inf

    dists = np.linalg.norm(pocket_coords - coord, axis=1)
    return float(dists.min())


# =========================================================
# Variant annotation
# =========================================================

def annotate_variants(
    variants: List[Dict],
    structure: PDB.Structure.Structure,
    sasa_map: Dict,
    ss_map: Dict,
    pdb_path: str,
    pocket_df: pd.DataFrame,
    chain_id: Optional[str],
) -> pd.DataFrame:
    """
    For each variant, find the corresponding PDB residue and annotate
    with: chain, resname, sasa, solvent_class, secondary_structure, b_factor,
    pocket_dist_A.
    """
    model = structure[0]

    # Build resseq → (chain, residue) index
    residue_index: Dict[Tuple, PDB.Residue.Residue] = {}
    for chain in model.get_chains():
        cid = chain.get_id()
        if chain_id and cid != chain_id:
            continue
        for res in chain.get_residues():
            hetflag, resseq, _ = res.get_id()
            if hetflag.strip():
                continue
            residue_index[(cid, resseq)] = res

    ss_map_parsed = parse_secondary_structure(pdb_path)

    rows = []
    for v in variants:
        pos = v["pos"]
        # Find matching residue — try specified chain or all chains
        res_obj = None
        found_chain = None
        if chain_id:
            res_obj = residue_index.get((chain_id, pos))
            found_chain = chain_id
        else:
            for (cid, rseq), r in residue_index.items():
                if rseq == pos:
                    res_obj = r
                    found_chain = cid
                    break

        if res_obj is None:
            print(f"[WARN] Residue {pos} not found in structure — variant {v['raw']} skipped.")
            continue

        resname_3 = res_obj.get_resname().strip()
        resname_1 = THREE_TO_ONE.get(resname_3, "X")
        bfactor = float(np.mean([a.get_bfactor() for a in res_obj.get_atoms()]))

        sasa = sasa_map.get((found_chain, pos), np.nan)
        # Solvent classification: buried <25, intermediate 25-100, exposed >100
        if np.isnan(sasa):
            solvent_class = "unknown"
        elif sasa < 25:
            solvent_class = "buried"
        elif sasa < 100:
            solvent_class = "intermediate"
        else:
            solvent_class = "exposed"

        ss_label = ss_map_parsed.get((found_chain, pos), "C")
        ss_full = {"H": "Helix", "E": "Strand", "C": "Coil"}.get(ss_label, "Coil")

        pocket_dist = nearest_pocket_distance(pos, found_chain, structure, pocket_df)
        near_pocket = pocket_dist < 8.0

        # Verify WT aa matches structure
        structure_aa = resname_1
        wt_match = (structure_aa == v["wt_aa"]) if v["wt_aa"] != "?" else None

        rows.append({
            "variant":              v["raw"],
            "wt_aa":                v["wt_aa"],
            "pos":                  pos,
            "mut_aa":               v["mut_aa"],
            "impact_class":         v["impact_class"],
            "chain":                found_chain,
            "resname_3letter":      resname_3,
            "structure_aa_1letter": structure_aa,
            "wt_matches_structure": wt_match,
            "sasa_A2":              round(sasa, 2) if not np.isnan(sasa) else np.nan,
            "solvent_exposure":     solvent_class,
            "secondary_structure":  ss_full,
            "b_factor":             round(bfactor, 2),
            "nearest_pocket_dist_A": round(pocket_dist, 2) if not np.isinf(pocket_dist) else np.nan,
            "near_pocket_lt8A":     near_pocket,
        })

    return pd.DataFrame(rows)


# =========================================================
# Visualisations
# =========================================================

def plot_linear_variant_map(
    variant_df: pd.DataFrame,
    total_residues: int,
    label: str,
    outdir: str,
    prefix: str,
) -> None:
    """
    Linear protein schematic with coloured lollipops at variant positions.
    """
    fig, ax = plt.subplots(figsize=(14, 3.5), dpi=300)

    # Backbone
    ax.axhline(0, color="#CCCCCC", lw=6, solid_capstyle="round", zorder=1)

    # Lollipops
    for _, row in variant_df.iterrows():
        color = IMPACT_COLORS.get(row["impact_class"], "#888888")
        pos_x = row["pos"]
        ax.plot([pos_x, pos_x], [0, 1], color=color, lw=1.2, zorder=2)
        ax.scatter(pos_x, 1, color=color, s=60, zorder=3, edgecolors="white", linewidths=0.5)
        ax.text(
            pos_x, 1.12,
            row["variant"],
            ha="center", va="bottom",
            fontsize=7.5,
            color=color,
            rotation=60,
        )

    ax.set_xlim(-total_residues * 0.02, total_residues * 1.05)
    ax.set_ylim(-0.5, 2.0)
    ax.set_xlabel("Residue position", fontsize=11)
    ax.set_yticks([])
    ax.set_title(f"Variant map  ·  {label}", fontsize=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    # Legend
    handles = [
        mpatches.Patch(color=c, label=lbl.replace("_", " ").capitalize())
        for lbl, c in IMPACT_COLORS.items()
        if lbl in variant_df["impact_class"].values
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=9, framealpha=0.8)

    plt.tight_layout()
    save_fig(fig, os.path.join(outdir, f"{prefix}.variant_linear_map"))
    print("[INFO] Linear variant map saved.")


def render_variant_viewer(
    pdb_path: str,
    variant_df: pd.DataFrame,
    label: str,
    outdir: str,
    prefix: str,
) -> str:
    """
    Build a self-contained HTML 3D viewer with variants shown as
    coloured spheres overlaid on the full protein cartoon.
    """
    try:
        import py3Dmol  # type: ignore
    except ImportError:
        print("[WARN] py3Dmol not installed — skipping 3D viewer.")
        return ""

    with open(pdb_path) as f:
        pdb_data = f.read()

    # Pre-compute Cα coordinates from BioPython for reliable label placement
    parser = PDBParser(QUIET=True)
    structure_obj = parser.get_structure("var", pdb_path)
    ca_coords: Dict[Tuple[str, int], Tuple[float, float, float]] = {}
    for model in structure_obj:
        for chain in model:
            for residue in chain:
                if "CA" in residue:
                    ca = residue["CA"].get_vector()
                    ca_coords[(chain.id, residue.id[1])] = (
                        float(ca[0]), float(ca[1]), float(ca[2]),
                    )
        break  # first model only

    view = py3Dmol.view(width=900, height=650)
    view.addModel(pdb_data, "pdb")

    # Full-protein cartoon background
    view.setStyle({}, {"cartoon": {"color": "lightgray", "opacity": 0.65}})

    # Each variant: sphere at Cα + stick for side chain + label
    legend_items: List[str] = []
    seen_classes = set()

    for _, row in variant_df.iterrows():
        color = IMPACT_COLORS.get(row["impact_class"], "#888888")
        sel = {"chain": row["chain"], "resi": str(int(row["pos"]))}

        # Coloured sphere at Cα
        view.addStyle(
            {**sel, "atom": "CA"},
            {"sphere": {"color": color, "radius": 0.9}},
        )
        # Sticks for full residue
        view.addStyle(sel, {"stick": {"color": color, "radius": 0.18}})

        # Label with explicit x/y/z coordinates (py3Dmol position selector
        # with chain/resi/atom is unreliable; explicit coords always work)
        coord_key = (row["chain"], int(row["pos"]))
        if coord_key in ca_coords:
            x, y, z = ca_coords[coord_key]
            view.addLabel(
                row["variant"],
                {
                    "position":         {"x": x, "y": y + 1.5, "z": z},
                    "font":             "Arial",
                    "fontSize":         12,
                    "fontColor":        "black",
                    "fontOpacity":      1.0,
                    "backgroundOpacity": 0.85,
                    "backgroundColor":  "white",
                    "borderThickness":  1.0,
                    "borderColor":      color,
                    "borderOpacity":    1.0,
                    "inFront":          True,
                    "showBackground":   True,
                },
            )

        if row["impact_class"] not in seen_classes:
            seen_classes.add(row["impact_class"])
            legend_items.append(
                f'<span style="display:inline-flex;align-items:center;margin:0 8px;">'
                f'<span style="width:12px;height:12px;border-radius:50%;background:{color};'
                f'display:inline-block;margin-right:4px;border:1px solid #888;"></span>'
                f'{row["impact_class"].replace("_"," ").capitalize()}</span>'
            )

    view.zoomTo()
    view.setBackgroundColor("white")

    viewer_html = view._make_html()
    legend_html = (
        f'<div style="font-family:sans-serif;font-size:13px;padding:8px 4px;'
        f'background:#f5f5f5;border-top:1px solid #ddd;text-align:center;">'
        f'<strong>{label}</strong> &nbsp;|&nbsp; '
        f'{len(variant_df)} variant(s) &nbsp;|&nbsp; '
        + "".join(legend_items) +
        f'</div>'
    )
    full_html = viewer_html.replace("</body>", legend_html + "\n</body>")

    out_path = os.path.join(outdir, f"{prefix}.variant_map.html")
    with open(out_path, "w") as f:
        f.write(full_html)

    print(f"[INFO] 3D variant viewer saved → {out_path}")
    return out_path


# =========================================================
# CLI
# =========================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Map amino-acid variants onto a protein 3D structure.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Map variants onto a PDB structure
  python protein_variant_mapper.py --pdb-id 1TIM --variants "A23V,G45S,K132E" --outdir results/

  # Map variants onto AlphaFold structure, classify by UniProt
  python protein_variant_mapper.py --uniprot P00533 --variants "T790M,L858R,G719S" \\
      --fetch-uniprot-variants --outdir results/

  # Load variants from file, specify chain
  python protein_variant_mapper.py --pdb-id 2HHB --variants-file variants.tsv \\
      --chain A --outdir results/

  # Include pocket proximity if you already ran the visualizer pocket module
  python protein_variant_mapper.py --pdb-id 1TIM --variants "A23V,G45S" \\
      --pocket-residues results/1TIM.pocket_residues.tsv --outdir results/
""",
    )

    # Structure source (mutually exclusive)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--pdb-id",    metavar="ID",   help="RCSB PDB ID")
    src.add_argument("--uniprot",   metavar="ACC",  help="UniProt accession → AlphaFold structure")
    src.add_argument("--pdb-file",  metavar="PATH", help="Local PDB file")

    # Variant input (mutually exclusive)
    var = p.add_mutually_exclusive_group(required=True)
    var.add_argument("--variants",      metavar="STR",  help='Comma-separated variants, e.g. "A123V,G45S"')
    var.add_argument("--variants-file", metavar="PATH", help="File with one variant per line, or TSV with variant/impact_class columns")

    p.add_argument("--chain",                 metavar="C",    default=None, help="Restrict to this chain (default: all)")
    p.add_argument("--fetch-uniprot-variants",action="store_true",          help="Also fetch UniProt natural variants and annotate with disease association")
    p.add_argument("--uniprot-acc",           metavar="ACC",  default=None, help="UniProt accession for variant fetching (if different from --uniprot)")
    p.add_argument("--pocket-residues",       metavar="PATH", default=None, help="Pocket residues TSV from protein-structure-visualizer pocket module")
    p.add_argument("--probe-radius",          metavar="F",    type=float, default=1.4, help="SASA probe radius in Å (default: 1.4)")
    p.add_argument("--outdir",                default="results", help="Output directory (default: results)")
    p.add_argument("--prefix",                default=None,      help="Output file prefix (default: auto)")

    return p.parse_args()


def main() -> None:
    args = parse_args()
    ensure_dir(args.outdir)

    # ── Structure ─────────────────────────────────────────────────────────
    if args.pdb_id:
        pdb_path = fetch_pdb_structure(args.pdb_id, args.outdir)
        label = args.pdb_id.upper()
        uniprot_acc = args.uniprot_acc
    elif args.uniprot:
        pdb_path = fetch_alphafold_structure(args.uniprot, args.outdir)
        label = f"AF-{args.uniprot.upper()}"
        uniprot_acc = args.uniprot
    else:
        pdb_path = args.pdb_file
        label = os.path.splitext(os.path.basename(pdb_path))[0]
        uniprot_acc = args.uniprot_acc

    prefix = args.prefix or f"{label}_variants"

    structure = load_structure(pdb_path, label)

    # ── Variants ──────────────────────────────────────────────────────────
    if args.variants:
        variants = parse_variants(args.variants)
    else:
        variants = load_variants_file(args.variants_file)

    print(f"[INFO] {len(variants)} variants to map.")

    # Optionally merge UniProt natural variants for classification
    if args.fetch_uniprot_variants and uniprot_acc:
        uniprot_vars = fetch_uniprot_variants(uniprot_acc)
        uniprot_pos_map = {(v["wt_aa"], v["pos"]): v["impact_class"] for v in uniprot_vars}
        for v in variants:
            key = (v["wt_aa"], v["pos"])
            if key in uniprot_pos_map and v["impact_class"] == "custom":
                v["impact_class"] = uniprot_pos_map[key]

    # ── SASA ──────────────────────────────────────────────────────────────
    print("[INFO] Computing per-residue SASA...")
    sasa_map = compute_sasa(structure, probe_radius=args.probe_radius)

    # ── Secondary structure ───────────────────────────────────────────────
    ss_map = parse_secondary_structure(pdb_path)

    # ── Pocket proximity ──────────────────────────────────────────────────
    pocket_df = pd.DataFrame()
    if args.pocket_residues:
        pocket_df = load_pocket_residues(args.pocket_residues)

    # ── Annotate variants ─────────────────────────────────────────────────
    variant_df = annotate_variants(
        variants, structure, sasa_map, ss_map, pdb_path, pocket_df, args.chain
    )

    if variant_df.empty:
        print("[ERROR] No variants could be mapped to the structure. Check positions and chain.")
        return

    # ── Total sequence length (approximate from model) ────────────────────
    all_resseqs = []
    for chain in structure[0].get_chains():
        for res in chain.get_residues():
            hetflag, resseq, _ = res.get_id()
            if not hetflag.strip():
                all_resseqs.append(resseq)
    total_residues = max(all_resseqs) if all_resseqs else variant_df["pos"].max()

    # ── Save table ────────────────────────────────────────────────────────
    summary_path = os.path.join(args.outdir, f"{prefix}.variant_summary.tsv")
    variant_df.to_csv(summary_path, sep="\t", index=False)

    # ── Plots ─────────────────────────────────────────────────────────────
    plot_linear_variant_map(variant_df, total_residues, label, args.outdir, prefix)
    render_variant_viewer(pdb_path, variant_df, label, args.outdir, prefix)

    # ── Print summary ─────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print(f"  Structure: {label}")
    print(f"  Variants mapped: {len(variant_df)}")
    for cls, grp in variant_df.groupby("impact_class"):
        print(f"    {cls:20s}: {len(grp)}")
    print(f"\n  Buried variants:      {(variant_df['solvent_exposure']=='buried').sum()}")
    print(f"  Exposed variants:     {(variant_df['solvent_exposure']=='exposed').sum()}")
    if not pocket_df.empty:
        print(f"  Near pocket (<8 Å):  {variant_df['near_pocket_lt8A'].sum()}")
    print("=" * 60)


if __name__ == "__main__":
    main()
