#!/usr/bin/env python3
"""
Protein Structure Alignment

Fetches two protein structures (PDB, AlphaFold, or local PDB file),
superimposes them by Cα atoms using Biopython's Superimposer, and produces:

  - Global RMSD and alignment statistics (TSV)
  - Per-residue Cα distance after superimposition (TSV + plot)
  - Self-contained HTML 3D viewer showing both aligned structures (py3Dmol)
  - Superimposed structure saved as PDB file

Typical use-cases:
  - Wild-type vs mutant comparison
  - Apo vs holo (ligand-bound) conformation
  - Homolog comparison across species
  - AlphaFold predicted vs experimentally determined

Dependencies:
  biopython>=1.81, matplotlib>=3.7, numpy>=1.24, pandas>=2.0,
  requests>=2.28, py3Dmol>=2.0
"""

import argparse
import io
import os
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from Bio import PDB
from Bio.PDB import MMCIFParser, PDBParser, Superimposer

# =========================================================
# Constants
# =========================================================

RCSB_FILE     = "https://files.rcsb.org/download"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api"
ALPHAFOLD_FILE = "https://alphafold.ebi.ac.uk/files"
HEADERS = {"User-Agent": "protein-structure-alignment/0.1 (bioinfor-claw)"}

STRUCT_COLORS = ["#1f77b4", "#d62728"]   # blue = structure 1, red = structure 2


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
# Structure fetching
# =========================================================

def fetch_pdb_structure(pdb_id: str, outdir: str) -> str:
    pdb_id = pdb_id.upper()
    out_path = os.path.join(outdir, f"{pdb_id}.pdb")
    if os.path.exists(out_path):
        return out_path
    url = f"{RCSB_FILE}/{pdb_id}.pdb"
    r = _get(url)
    with open(out_path, "w") as f:
        f.write(r.text)
    print(f"[INFO] Downloaded PDB: {pdb_id} → {out_path}")
    return out_path


def fetch_alphafold_structure(uniprot_acc: str, outdir: str) -> str:
    uniprot_acc = uniprot_acc.upper()
    out_path = os.path.join(outdir, f"AF_{uniprot_acc}.pdb")
    if os.path.exists(out_path):
        return out_path
    # Get latest model version
    meta_url = f"{ALPHAFOLD_API}/prediction/{uniprot_acc}"
    r = _get(meta_url)
    entries = r.json()
    if not entries:
        raise ValueError(f"No AlphaFold entry for {uniprot_acc}")
    pdb_url = entries[0]["pdbUrl"]
    r2 = _get(pdb_url)
    with open(out_path, "w") as f:
        f.write(r2.text)
    print(f"[INFO] Downloaded AlphaFold: {uniprot_acc} → {out_path}")
    return out_path


def load_structure(pdb_path: str, struct_id: str = "STR") -> PDB.Structure.Structure:
    parser = PDBParser(QUIET=True)
    return parser.get_structure(struct_id, pdb_path)


# =========================================================
# Alignment core
# =========================================================

def get_ca_atoms(
    structure: PDB.Structure.Structure,
    chain_id: Optional[str] = None,
    res_start: Optional[int] = None,
    res_end: Optional[int] = None,
) -> List[PDB.Atom.Atom]:
    """Extract Cα atoms from a structure, optionally filtered by chain and residue range."""
    model = structure[0]
    ca_atoms = []
    chains = [model[chain_id]] if chain_id and chain_id in model else model.get_chains()
    for chain in chains:
        for res in chain.get_residues():
            hetflag, resseq, _ = res.get_id()
            if hetflag.strip():
                continue
            if res_start is not None and resseq < res_start:
                continue
            if res_end is not None and resseq > res_end:
                continue
            if "CA" in res:
                ca_atoms.append(res["CA"])
    return ca_atoms


def pair_ca_atoms(
    ca1: List[PDB.Atom.Atom],
    ca2: List[PDB.Atom.Atom],
) -> Tuple[List[PDB.Atom.Atom], List[PDB.Atom.Atom]]:
    """
    Pair Cα atoms by sequence position (resseq).
    Returns only residues present in both structures at the same sequence number.
    Pairing is done by (chain_id, resseq) key.
    """
    def key(atom):
        res = atom.get_parent()
        chain = res.get_parent().get_id()
        _, resseq, _ = res.get_id()
        return (chain, resseq)

    map1 = {key(a): a for a in ca1}
    map2 = {key(a): a for a in ca2}
    common = sorted(set(map1) & set(map2))
    if not common:
        raise ValueError(
            "No common (chain, resseq) Cα pairs found. "
            "Try specifying --chain1/--chain2 or --res-start/--res-end."
        )
    return [map1[k] for k in common], [map2[k] for k in common]


def superimpose(
    struct1: PDB.Structure.Structure,
    struct2: PDB.Structure.Structure,
    chain1: Optional[str] = None,
    chain2: Optional[str] = None,
    res_start: Optional[int] = None,
    res_end: Optional[int] = None,
) -> Tuple[float, pd.DataFrame, PDB.Structure.Structure]:
    """
    Superimpose struct2 onto struct1 using paired Cα atoms.

    Returns:
        rmsd        — global RMSD after superimposition (Å)
        per_res_df  — per-residue Cα distance DataFrame
        struct2     — struct2 with all atoms moved into struct1 frame
    """
    ca1_all = get_ca_atoms(struct1, chain1, res_start, res_end)
    ca2_all = get_ca_atoms(struct2, chain2, res_start, res_end)

    paired1, paired2 = pair_ca_atoms(ca1_all, ca2_all)
    n_paired = len(paired1)
    print(f"[INFO] Paired {n_paired} Cα atoms for superimposition.")

    sup = Superimposer()
    sup.set_atoms(paired1, paired2)
    # Apply rotation/translation to ALL atoms in struct2
    all_atoms2 = list(struct2.get_atoms())
    sup.apply(all_atoms2)

    rmsd = sup.rms
    print(f"[INFO] Global RMSD = {rmsd:.3f} Å over {n_paired} Cα pairs.")

    # Per-residue Cα distances after superimposition
    rows = []
    for a1, a2 in zip(paired1, paired2):
        res1 = a1.get_parent()
        _, resseq, _ = res1.get_id()
        chain = res1.get_parent().get_id()
        dist = float(np.linalg.norm(a1.get_vector().get_array() - a2.get_vector().get_array()))
        rows.append({
            "chain":   chain,
            "resseq":  resseq,
            "resname": res1.get_resname(),
            "ca_dist_A": dist,
        })

    per_res_df = pd.DataFrame(rows)
    return rmsd, per_res_df, struct2


# =========================================================
# Outputs
# =========================================================

def save_superimposed_pdb(
    struct2: PDB.Structure.Structure,
    outdir: str,
    prefix: str,
) -> str:
    """Write the superimposed (rotated/translated) structure 2 to PDB."""
    io_obj = PDB.PDBIO()
    io_obj.set_structure(struct2)
    out_path = os.path.join(outdir, f"{prefix}.superimposed.pdb")
    io_obj.save(out_path)
    print(f"[INFO] Superimposed structure saved → {out_path}")
    return out_path


def plot_per_residue_rmsd(
    per_res_df: pd.DataFrame,
    rmsd: float,
    label1: str,
    label2: str,
    outdir: str,
    prefix: str,
) -> None:
    fig, ax = plt.subplots(figsize=(14, 4), dpi=300)

    x = np.arange(len(per_res_df))
    dists = per_res_df["ca_dist_A"].values

    # Colour by magnitude: <1Å green, 1-3Å orange, >3Å red
    colors = []
    for d in dists:
        if d < 1.0:
            colors.append("#2ECC40")
        elif d < 3.0:
            colors.append("#FF851B")
        else:
            colors.append("#D62728")

    ax.bar(x, dists, color=colors, width=1.0, linewidth=0)
    ax.axhline(1.0, color="#2ECC40", lw=0.8, ls="--", alpha=0.7, label="1.0 Å")
    ax.axhline(3.0, color="#D62728", lw=0.8, ls="--", alpha=0.7, label="3.0 Å")

    # X tick labels — every 20 residues
    step = max(1, len(per_res_df) // 20)
    tick_idx = x[::step]
    tick_labels = per_res_df["resseq"].values[::step]
    ax.set_xticks(tick_idx)
    ax.set_xticklabels(tick_labels, fontsize=8, rotation=45, ha="right")

    ax.set_xlabel("Residue number", fontsize=11)
    ax.set_ylabel("Cα distance (Å)", fontsize=11)
    ax.set_title(
        f"Per-residue Cα distance after superimposition\n"
        f"{label1}  vs  {label2}  |  Global RMSD = {rmsd:.3f} Å",
        fontsize=12,
    )
    ax.legend(fontsize=9, title="Thresholds", title_fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    save_fig(fig, os.path.join(outdir, f"{prefix}.per_residue_rmsd"))
    print(f"[INFO] Per-residue RMSD plot saved.")


def render_superimposed_viewer(
    pdb_path1: str,
    pdb_path2: str,
    label1: str,
    label2: str,
    per_res_df: pd.DataFrame,
    rmsd: float,
    outdir: str,
    prefix: str,
) -> str:
    """
    Build a self-contained HTML 3D viewer showing both structures
    superimposed in two colours.

    Structure 1 = blue cartoon.
    Structure 2 = red cartoon.
    Regions with per-residue Cα distance > 3 Å are highlighted with
    orange sticks on both structures to mark flexible/divergent regions.
    """
    try:
        import py3Dmol  # type: ignore
    except ImportError:
        print("[WARN] py3Dmol not installed — skipping 3D viewer.")
        return ""

    with open(pdb_path1) as f:
        pdb1 = f.read()
    with open(pdb_path2) as f:
        pdb2 = f.read()

    view = py3Dmol.view(width=900, height=650)

    # Model 0 = structure 1 (blue)
    view.addModel(pdb1, "pdb")
    view.setStyle({"model": 0}, {"cartoon": {"color": STRUCT_COLORS[0], "opacity": 0.85}})

    # Model 1 = structure 2 (red)
    view.addModel(pdb2, "pdb")
    view.setStyle({"model": 1}, {"cartoon": {"color": STRUCT_COLORS[1], "opacity": 0.85}})

    # Highlight divergent residues (Cα dist > 3 Å) with orange sticks
    divergent = per_res_df[per_res_df["ca_dist_A"] > 3.0]
    if not divergent.empty:
        for chain_id, grp in divergent.groupby("chain"):
            resi_list = [str(r) for r in sorted(grp["resseq"])]
            resi_str = ",".join(resi_list)
            for model_idx in [0, 1]:
                view.addStyle(
                    {"model": model_idx, "chain": chain_id, "resi": resi_str},
                    {"stick": {"color": "#FF851B", "radius": 0.22}},
                )

    view.zoomTo()
    view.setBackgroundColor("white")

    viewer_html = view._make_html()

    # Legend footer
    n_div = len(divergent)
    pct_div = 100 * n_div / max(len(per_res_df), 1)
    legend_html = (
        f'<div style="font-family:sans-serif;font-size:13px;padding:8px 6px;'
        f'background:#f5f5f5;border-top:1px solid #ddd;text-align:center;">'
        f'<span style="color:{STRUCT_COLORS[0]};font-weight:bold;">■</span> {label1} &nbsp;|&nbsp; '
        f'<span style="color:{STRUCT_COLORS[1]};font-weight:bold;">■</span> {label2} &nbsp;|&nbsp; '
        f'Global RMSD = <strong>{rmsd:.3f} Å</strong> &nbsp;|&nbsp; '
        f'<span style="color:#FF851B;font-weight:bold;">■</span> Divergent &gt;3 Å '
        f'({n_div} residues, {pct_div:.1f}%)'
        f'</div>'
    )

    full_html = viewer_html.replace("</body>", legend_html + "\n</body>")
    out_path = os.path.join(outdir, f"{prefix}.superimposed.html")
    with open(out_path, "w") as f:
        f.write(full_html)

    print(f"[INFO] 3D viewer saved → {out_path}")
    return out_path


# =========================================================
# CLI
# =========================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Superimpose two protein structures and visualize alignment.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Two PDB structures
  python protein_structure_alignment.py --pdb1 1AKE --pdb2 4AKE --outdir results/

  # AlphaFold vs experimental
  python protein_structure_alignment.py --pdb1 1TIM --uniprot2 P00940 --outdir results/

  # Specify chains and residue range
  python protein_structure_alignment.py --pdb1 1AKE --pdb2 4AKE \\
      --chain1 A --chain2 A --res-start 1 --res-end 180 --outdir results/

  # Local PDB files
  python protein_structure_alignment.py --file1 wt.pdb --file2 mutant.pdb --outdir results/
""",
    )

    # Structure 1
    g1 = p.add_mutually_exclusive_group(required=True)
    g1.add_argument("--pdb1",    metavar="ID",   help="RCSB PDB ID for structure 1")
    g1.add_argument("--uniprot1",metavar="ACC",  help="UniProt accession → AlphaFold for structure 1")
    g1.add_argument("--file1",   metavar="PATH", help="Local PDB file for structure 1")

    # Structure 2
    g2 = p.add_mutually_exclusive_group(required=True)
    g2.add_argument("--pdb2",    metavar="ID",   help="RCSB PDB ID for structure 2")
    g2.add_argument("--uniprot2",metavar="ACC",  help="UniProt accession → AlphaFold for structure 2")
    g2.add_argument("--file2",   metavar="PATH", help="Local PDB file for structure 2")

    # Selection
    p.add_argument("--chain1",    metavar="C",  default=None, help="Chain ID to use in structure 1")
    p.add_argument("--chain2",    metavar="C",  default=None, help="Chain ID to use in structure 2")
    p.add_argument("--res-start", metavar="N",  type=int, default=None, help="First residue to include")
    p.add_argument("--res-end",   metavar="N",  type=int, default=None, help="Last residue to include")

    # Output
    p.add_argument("--outdir", default="results", help="Output directory (default: results)")
    p.add_argument("--prefix", default=None,      help="Output file prefix (default: auto)")

    return p.parse_args()


def main() -> None:
    args = parse_args()
    ensure_dir(args.outdir)

    # ── Fetch / load structure 1 ──────────────────────────────────────────
    if args.pdb1:
        pdb_path1 = fetch_pdb_structure(args.pdb1, args.outdir)
        label1 = args.pdb1.upper()
    elif args.uniprot1:
        pdb_path1 = fetch_alphafold_structure(args.uniprot1, args.outdir)
        label1 = f"AF-{args.uniprot1.upper()}"
    else:
        pdb_path1 = args.file1
        label1 = os.path.splitext(os.path.basename(pdb_path1))[0]

    # ── Fetch / load structure 2 ──────────────────────────────────────────
    if args.pdb2:
        pdb_path2 = fetch_pdb_structure(args.pdb2, args.outdir)
        label2 = args.pdb2.upper()
    elif args.uniprot2:
        pdb_path2 = fetch_alphafold_structure(args.uniprot2, args.outdir)
        label2 = f"AF-{args.uniprot2.upper()}"
    else:
        pdb_path2 = args.file2
        label2 = os.path.splitext(os.path.basename(pdb_path2))[0]

    prefix = args.prefix or f"{label1}_vs_{label2}"

    print(f"[INFO] Structure 1: {label1}  ({pdb_path1})")
    print(f"[INFO] Structure 2: {label2}  ({pdb_path2})")

    struct1 = load_structure(pdb_path1, label1)
    struct2 = load_structure(pdb_path2, label2)

    # ── Superimpose ───────────────────────────────────────────────────────
    rmsd, per_res_df, struct2_aligned = superimpose(
        struct1, struct2,
        chain1=args.chain1,
        chain2=args.chain2,
        res_start=args.res_start,
        res_end=args.res_end,
    )

    # ── Save aligned structure ────────────────────────────────────────────
    aligned_pdb_path = save_superimposed_pdb(struct2_aligned, args.outdir, prefix)

    # ── Per-residue table ─────────────────────────────────────────────────
    per_res_path = os.path.join(args.outdir, f"{prefix}.per_residue_rmsd.tsv")
    per_res_df.to_csv(per_res_path, sep="\t", index=False)

    # ── Summary table ─────────────────────────────────────────────────────
    n_paired = len(per_res_df)
    n_close  = (per_res_df["ca_dist_A"] < 1.0).sum()
    n_medium = ((per_res_df["ca_dist_A"] >= 1.0) & (per_res_df["ca_dist_A"] < 3.0)).sum()
    n_far    = (per_res_df["ca_dist_A"] >= 3.0).sum()

    summary = {
        "label1":            label1,
        "label2":            label2,
        "n_ca_pairs":        n_paired,
        "global_rmsd_A":     round(rmsd, 4),
        "n_close_lt1A":      int(n_close),
        "n_medium_1to3A":    int(n_medium),
        "n_divergent_gt3A":  int(n_far),
        "pct_divergent":     round(100 * n_far / max(n_paired, 1), 2),
        "max_ca_dist_A":     round(float(per_res_df["ca_dist_A"].max()), 4),
        "median_ca_dist_A":  round(float(per_res_df["ca_dist_A"].median()), 4),
    }
    summary_df = pd.DataFrame([summary])
    summary_path = os.path.join(args.outdir, f"{prefix}.alignment_summary.tsv")
    summary_df.to_csv(summary_path, sep="\t", index=False)

    # ── Plots ─────────────────────────────────────────────────────────────
    plot_per_residue_rmsd(per_res_df, rmsd, label1, label2, args.outdir, prefix)

    # ── 3D viewer ─────────────────────────────────────────────────────────
    render_superimposed_viewer(
        pdb_path1, aligned_pdb_path,
        label1, label2,
        per_res_df, rmsd,
        args.outdir, prefix,
    )

    # ── Final summary ─────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print(f"  Alignment: {label1}  vs  {label2}")
    print(f"  Paired Cα:   {n_paired}")
    print(f"  Global RMSD: {rmsd:.3f} Å")
    print(f"  < 1 Å:       {n_close} ({100*n_close/max(n_paired,1):.1f}%)")
    print(f"  1–3 Å:       {n_medium} ({100*n_medium/max(n_paired,1):.1f}%)")
    print(f"  > 3 Å:       {n_far} ({100*n_far/max(n_paired,1):.1f}%)")
    print("=" * 60)


if __name__ == "__main__":
    main()
