#!/usr/bin/env python3
"""
Protein Structure Visualizer

Fetches a PDB or AlphaFold structure and produces:
  - An interactive HTML 3D viewer (full protein or zoomed region)
  - A residue contact map (PNG + PDF)
  - A B-factor / pLDDT plot (PNG + PDF)
  - A secondary structure timeline (PNG + PDF)
  - Pocket / cavity detection (solvent-exposure + geometry based)
  - A protein-protein interaction (PPI) network via STRING API (PNG + PDF + TSV)

All structure data is fetched live from RCSB PDB or AlphaFold EBI.
PPI data is fetched from the STRING database REST API.

Supported modules:
  view         — interactive HTML 3D viewer (full protein or zoomed residue range)
  contact_map  — Cα distance contact map
  bfactor      — per-residue B-factor / pLDDT plot
  secondary    — secondary structure composition and timeline
  pocket       — solvent-accessible cavity detection
  ppi          — protein-protein interaction network (STRING)
  all          — run all modules

Dependencies:
  biopython>=1.81, matplotlib>=3.7, numpy>=1.24, pandas>=2.0,
  requests>=2.28, networkx>=3.0, py3Dmol>=2.0
"""

import argparse
import io
import json
import os
import textwrap
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests

# Biopython — PDB parsing
from Bio import PDB
from Bio.PDB import DSSP, MMCIFParser, PDBParser, Select
from Bio.PDB.SASA import ShrakeRupley

# =========================================================
# Constants
# =========================================================

RCSB_FILE      = "https://files.rcsb.org/download"
ALPHAFOLD_API  = "https://alphafold.ebi.ac.uk/api"
ALPHAFOLD_FILE = "https://alphafold.ebi.ac.uk/files"
STRING_API     = "https://string-db.org/api"

HEADERS = {"User-Agent": "protein-structure-visualizer/0.1 (bioinfor-claw)"}

SS_COLORS = {"H": "#D62728", "E": "#1F77B4", "C": "#AAAAAA", "-": "#AAAAAA"}
SS_LABELS = {"H": "Helix", "E": "Strand", "C": "Coil / Loop", "-": "Coil / Loop"}

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
    """Download a PDB file from RCSB. Returns local file path."""
    pdb_id = pdb_id.upper()
    path = os.path.join(outdir, f"{pdb_id}.pdb")
    if os.path.exists(path):
        return path
    url = f"{RCSB_FILE}/{pdb_id}.pdb"
    r = _get(url)
    with open(path, "w") as f:
        f.write(r.text)
    print(f"[INFO] PDB structure saved: {path}")
    return path


def fetch_alphafold_structure(uniprot_acc: str, outdir: str) -> str:
    """Download AlphaFold2 PDB file from EBI. Returns local file path."""
    # Get latest version
    url = f"{ALPHAFOLD_API}/prediction/{uniprot_acc}"
    r = _get(url)
    entries = r.json()
    if not entries:
        raise ValueError(f"No AlphaFold entry for {uniprot_acc}")
    entry = entries[0]
    pdb_url = entry.get("pdbUrl")
    if not pdb_url:
        raise ValueError(f"No PDB URL in AlphaFold entry for {uniprot_acc}")

    path = os.path.join(outdir, f"AF_{uniprot_acc}.pdb")
    if not os.path.exists(path):
        r2 = _get(pdb_url)
        with open(path, "w") as f:
            f.write(r2.text)
        print(f"[INFO] AlphaFold structure saved: {path}")
    return path


def parse_structure(pdb_path: str, structure_id: str = "STR") -> PDB.Structure.Structure:
    """Parse PDB file using Biopython. Returns Structure object."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(structure_id, pdb_path)
    return structure


def get_residue_data(structure: PDB.Structure.Structure, chain_id: Optional[str] = None) -> pd.DataFrame:
    """
    Extract per-residue data: residue number, name, chain, Cα coords, B-factor.
    Returns a DataFrame suitable for most downstream analyses.
    """
    rows = []
    model = structure[0]
    chains = [model[chain_id]] if chain_id and chain_id in model else model.get_chains()

    for chain in chains:
        for residue in chain.get_residues():
            hetflag, resseq, icode = residue.get_id()
            if hetflag.strip():   # skip HET atoms (ligands, water)
                continue
            res_name = residue.get_resname()
            ca = residue["CA"] if "CA" in residue else None
            if ca is None:
                continue
            coord = ca.get_vector().get_array()
            bfactor = ca.get_bfactor()
            rows.append({
                "chain":    chain.get_id(),
                "resseq":   resseq,
                "resname":  res_name,
                "x":        coord[0],
                "y":        coord[1],
                "z":        coord[2],
                "bfactor":  bfactor,
            })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["chain", "resseq"]).reset_index(drop=True)
        df["index"] = range(len(df))
    return df


# =========================================================
# Module 1: Interactive HTML viewer (py3Dmol)
# =========================================================

def render_html_viewer(
    pdb_path: str,
    structure_id: str,
    outdir: str,
    prefix: str,
    zoom_chain: Optional[str] = None,
    zoom_start: Optional[int] = None,
    zoom_end:   Optional[int] = None,
    style: str = "cartoon",
    color_scheme: str = "spectrum",
    is_alphafold: bool = False,
) -> str:
    """
    Generate an interactive 3D HTML viewer using py3Dmol embedded in a standalone HTML page.
    Works in any browser without installation.
    Returns the path to the HTML file.
    """
    try:
        import py3Dmol
    except ImportError:
        raise ImportError("py3Dmol is required for HTML visualization. Install with: pip install py3Dmol")

    with open(pdb_path, "r") as f:
        pdb_data = f.read()

    # Build viewer
    view = py3Dmol.view(width=900, height=600)
    view.addModel(pdb_data, "pdb")

    # Base style
    if color_scheme == "spectrum":
        view.setStyle({"cartoon": {"color": "spectrum"}})
    elif color_scheme == "chain":
        view.setStyle({"cartoon": {"colorscheme": "chain"}})
    elif color_scheme == "bfactor" and is_alphafold:
        # AlphaFold pLDDT colouring: >90 blue, 70-90 cyan, 50-70 yellow, <50 orange
        view.setStyle({"cartoon": {"colorscheme": {"prop": "b", "gradient": "roygb", "min": 50, "max": 90}}})
    else:
        view.setStyle({"cartoon": {"color": "spectrum"}})

    title_parts = [structure_id]

    # Zoom to region if requested
    if zoom_chain and zoom_start and zoom_end:
        sel = {"chain": zoom_chain, "resi": f"{zoom_start}-{zoom_end}"}
        # Highlight selection
        view.addStyle(sel, {"stick": {"colorscheme": "yellowCarbon", "radius": 0.3}})
        view.addStyle(sel, {"cartoon": {"color": "#FF6600"}})
        view.zoomTo(sel)
        title_parts.append(f"Chain {zoom_chain} res {zoom_start}–{zoom_end}")
    elif zoom_chain:
        view.zoomTo({"chain": zoom_chain})
        title_parts.append(f"Chain {zoom_chain}")
    else:
        view.zoomTo()

    view.setBackgroundColor("white")

    # Add surface (semi-transparent) if showing a pocket region
    if zoom_start and zoom_end and zoom_chain:
        view.addSurface(
            py3Dmol.VDW,
            {"opacity": 0.3, "color": "lightblue"},
            {"chain": zoom_chain, "resi": f"{zoom_start}-{zoom_end}"}
        )

    title = " · ".join(title_parts)

    # Embed in self-contained HTML
    viewer_html = view._make_html()
    html_content = f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>Protein Viewer — {title}</title>
  <style>
    body {{ font-family: Arial, sans-serif; background: #f5f5f5; margin: 0; padding: 20px; }}
    h2 {{ color: #333; margin-bottom: 4px; }}
    .subtitle {{ color: #666; font-size: 13px; margin-bottom: 16px; }}
    .viewer-wrap {{ background: white; border: 1px solid #ddd; border-radius: 6px;
                   box-shadow: 0 2px 8px rgba(0,0,0,0.08); display: inline-block; padding: 12px; }}
    .legend {{ margin-top: 12px; font-size: 12px; color: #555; }}
  </style>
</head>
<body>
  <h2>🧬 {title}</h2>
  <div class="subtitle">
    Interactive 3D viewer · Drag to rotate · Scroll to zoom · Right-click for options
    {"· <b>AlphaFold2 model</b> — colour encodes pLDDT confidence (blue = high, orange = low)" if is_alphafold else ""}
  </div>
  <div class="viewer-wrap">
    {viewer_html}
  </div>
  <div class="legend">
    Generated by <b>bioinfor-claw / protein-structure-visualizer</b>
  </div>
</body>
</html>"""

    suffix = "zoomed" if (zoom_start or zoom_chain) else "full"
    out_path = os.path.join(outdir, f"{prefix}.{suffix}_view.html")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html_content)

    print(f"[INFO] HTML viewer saved: {out_path}")
    return out_path


# =========================================================
# Module 2: Contact map
# =========================================================

def plot_contact_map(
    res_df: pd.DataFrame,
    outdir: str,
    prefix: str,
    threshold_A: float = 8.0,
    chain_id: Optional[str] = None,
    zoom_start: Optional[int] = None,
    zoom_end:   Optional[int] = None,
) -> None:
    """
    Compute and plot a Cα–Cα distance contact map.
    Contacts are defined as residue pairs with Cα distance ≤ threshold_A.
    """
    df = res_df.copy()
    if chain_id:
        df = df[df["chain"] == chain_id]
    if zoom_start and zoom_end:
        df = df[(df["resseq"] >= zoom_start) & (df["resseq"] <= zoom_end)]
    if df.empty:
        print("[WARN] No residues for contact map after filtering.")
        return

    coords = df[["x", "y", "z"]].values
    n = len(coords)
    dist = np.sqrt(((coords[:, None, :] - coords[None, :, :]) ** 2).sum(axis=-1))
    contact = (dist <= threshold_A).astype(float)
    np.fill_diagonal(contact, 0)

    labels = df["resseq"].tolist()
    tick_step = max(1, n // 10)

    fig, ax = plt.subplots(figsize=(7, 6), dpi=300)
    im = ax.imshow(contact, cmap="binary", origin="lower", aspect="equal")
    ax.set_xticks(range(0, n, tick_step))
    ax.set_yticks(range(0, n, tick_step))
    ax.set_xticklabels(labels[::tick_step], fontsize=7, rotation=45)
    ax.set_yticklabels(labels[::tick_step], fontsize=7)
    ax.set_xlabel("Residue number", fontsize=11)
    ax.set_ylabel("Residue number", fontsize=11)

    region = f" (chain {chain_id}" + (f", res {zoom_start}–{zoom_end})" if zoom_start else ")") if chain_id else ""
    ax.set_title(f"Contact Map{region}  ·  {prefix}\n(Cα–Cα ≤ {threshold_A} Å)", fontsize=12)
    plt.tight_layout()
    save_fig(fig, os.path.join(outdir, f"{prefix}.contact_map"))
    print(f"[INFO] Contact map saved.")

    # Save distance matrix as TSV
    dist_df = pd.DataFrame(dist, index=labels, columns=labels)
    dist_df.to_csv(os.path.join(outdir, f"{prefix}.ca_distances.tsv"), sep="\t")


# =========================================================
# Module 3: B-factor / pLDDT plot
# =========================================================

def plot_bfactor(
    res_df: pd.DataFrame,
    outdir: str,
    prefix: str,
    is_alphafold: bool = False,
    chain_id: Optional[str] = None,
    zoom_start: Optional[int] = None,
    zoom_end:   Optional[int] = None,
) -> None:
    """Per-residue B-factor (crystal) or pLDDT (AlphaFold) line plot."""
    df = res_df.copy()
    if chain_id:
        df = df[df["chain"] == chain_id]
    if zoom_start and zoom_end:
        df = df[(df["resseq"] >= zoom_start) & (df["resseq"] <= zoom_end)]
    if df.empty:
        print("[WARN] No residues for B-factor plot after filtering.")
        return

    label = "pLDDT (AlphaFold confidence)" if is_alphafold else "B-factor (Å²)"
    title = f"{'pLDDT' if is_alphafold else 'B-factor'} per Residue  ·  {prefix}"

    fig, ax = plt.subplots(figsize=(max(8, len(df) * 0.04), 4), dpi=300)

    if is_alphafold:
        # Colour bands for pLDDT confidence tiers
        ax.axhspan(90, 100, color="#0053D6", alpha=0.08, label="Very high (>90)")
        ax.axhspan(70,  90, color="#65CBF3", alpha=0.08, label="Confident (70–90)")
        ax.axhspan(50,  70, color="#FFDB13", alpha=0.08, label="Low (50–70)")
        ax.axhspan(0,   50, color="#FF7D45", alpha=0.08, label="Very low (<50)")

        colors = []
        for b in df["bfactor"]:
            if b >= 90:   colors.append("#0053D6")
            elif b >= 70: colors.append("#65CBF3")
            elif b >= 50: colors.append("#FFDB13")
            else:         colors.append("#FF7D45")

        ax.scatter(df["resseq"], df["bfactor"], c=colors, s=4, zorder=3, linewidths=0)
        ax.plot(df["resseq"], df["bfactor"], color="#333333", lw=0.6, zorder=2, alpha=0.6)
        ax.set_ylim(0, 100)
        ax.legend(fontsize=8, loc="lower right", ncol=2)
    else:
        ax.plot(df["resseq"], df["bfactor"], color="#4E79A7", lw=1.0)
        ax.fill_between(df["resseq"], df["bfactor"], alpha=0.2, color="#4E79A7")

    ax.set_xlabel("Residue number", fontsize=11)
    ax.set_ylabel(label, fontsize=11)
    ax.set_title(title, fontsize=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    save_fig(fig, os.path.join(outdir, f"{prefix}.bfactor"))
    print(f"[INFO] B-factor/pLDDT plot saved.")

    # Save per-residue data
    df[["chain", "resseq", "resname", "bfactor"]].to_csv(
        os.path.join(outdir, f"{prefix}.bfactor.tsv"), sep="\t", index=False
    )


# =========================================================
# Module 4: Secondary structure timeline
# =========================================================

def plot_secondary_structure(
    structure: PDB.Structure.Structure,
    pdb_path: str,
    outdir: str,
    prefix: str,
    chain_id: Optional[str] = None,
) -> None:
    """
    Secondary structure assignment (DSSP) and composition bar chart.
    Falls back to a simplified helix/sheet detection from PDB HELIX/SHEET records
    if DSSP binary is unavailable.
    """
    model = structure[0]
    ss_data: List[Dict] = []

    # Try DSSP
    dssp_success = False
    try:
        dssp = DSSP(model, pdb_path, dssp="mkdssp")
        for key in dssp.property_dict:
            chain, (hetflag, resseq, icode) = key
            if chain_id and chain != chain_id:
                continue
            ss = dssp[key][2]
            ss_data.append({"chain": chain, "resseq": resseq, "ss": ss if ss in "HBEGITS-" else "-"})
        dssp_success = True
    except Exception:
        pass

    if not dssp_success:
        # Fallback: parse HELIX / SHEET records from PDB file
        helix_ranges: List[Tuple] = []
        sheet_ranges: List[Tuple] = []
        with open(pdb_path) as f:
            for line in f:
                if line.startswith("HELIX"):
                    try:
                        ch = line[19]
                        s = int(line[21:25])
                        e = int(line[33:37])
                        helix_ranges.append((ch, s, e))
                    except Exception:
                        pass
                elif line.startswith("SHEET"):
                    try:
                        ch = line[21]
                        s = int(line[22:26])
                        e = int(line[33:37])
                        sheet_ranges.append((ch, s, e))
                    except Exception:
                        pass

        for chain in model.get_chains():
            cid = chain.get_id()
            if chain_id and cid != chain_id:
                continue
            for res in chain.get_residues():
                hetflag, resseq, _ = res.get_id()
                if hetflag.strip():
                    continue
                ss = "-"
                for (ch, s, e) in helix_ranges:
                    if ch == cid and s <= resseq <= e:
                        ss = "H"
                        break
                if ss == "-":
                    for (ch, s, e) in sheet_ranges:
                        if ch == cid and s <= resseq <= e:
                            ss = "E"
                            break
                ss_data.append({"chain": cid, "resseq": resseq, "ss": ss})

    if not ss_data:
        print("[WARN] No secondary structure data retrieved.")
        return

    ss_df = pd.DataFrame(ss_data).sort_values(["chain", "resseq"]).reset_index(drop=True)
    # Simplify DSSP codes to H / E / C
    ss_df["ss_simple"] = ss_df["ss"].map(
        lambda s: "H" if s == "H" else ("E" if s in "BE" else "C")
    )

    # --- Timeline strip ---
    chains = ss_df["chain"].unique()
    n_chains = len(chains)
    fig, axes = plt.subplots(n_chains + 1, 1,
                              figsize=(max(10, len(ss_df) * 0.04 + 2), 2.5 + n_chains * 0.9),
                              dpi=300,
                              gridspec_kw={"height_ratios": [1] * n_chains + [2]})
    if n_chains == 1 and not isinstance(axes, np.ndarray):
        axes = [axes, axes]
    elif not isinstance(axes, np.ndarray):
        axes = list(axes)

    for idx, cid in enumerate(chains):
        ax = axes[idx]
        sub = ss_df[ss_df["chain"] == cid]
        for _, row in sub.iterrows():
            color = SS_COLORS.get(row["ss_simple"], "#AAAAAA")
            ax.barh(0, 1, left=row["resseq"], height=0.7,
                    color=color, linewidth=0)
        ax.set_xlim(sub["resseq"].min() - 1, sub["resseq"].max() + 1)
        ax.set_ylim(-0.5, 0.5)
        ax.set_yticks([0])
        ax.set_yticklabels([f"Chain {cid}"], fontsize=8)
        ax.set_xlabel("Residue number" if idx == n_chains - 1 else "", fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)

    # --- Composition bar ---
    ax_comp = axes[n_chains]
    counts = ss_df["ss_simple"].value_counts()
    total = len(ss_df)
    comp_data = {k: counts.get(k, 0) / total * 100 for k in ["H", "E", "C"]}
    bars = ax_comp.bar(["Helix (H)", "Strand (E)", "Coil (C)"],
                       [comp_data["H"], comp_data["E"], comp_data["C"]],
                       color=["#D62728", "#1F77B4", "#AAAAAA"],
                       edgecolor="white", width=0.5)
    for bar, val in zip(bars, comp_data.values()):
        ax_comp.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                     f"{val:.1f}%", ha="center", va="bottom", fontsize=9)
    ax_comp.set_ylabel("% of residues", fontsize=10)
    ax_comp.set_title("Secondary Structure Composition", fontsize=11)
    ax_comp.spines["top"].set_visible(False)
    ax_comp.spines["right"].set_visible(False)
    ax_comp.set_ylim(0, max(comp_data.values()) * 1.2 + 5)

    fig.suptitle(f"Secondary Structure  ·  {prefix}", fontsize=13, y=1.01)
    plt.tight_layout()
    save_fig(fig, os.path.join(outdir, f"{prefix}.secondary_structure"))
    ss_df.to_csv(os.path.join(outdir, f"{prefix}.secondary_structure.tsv"), sep="\t", index=False)
    print(f"[INFO] Secondary structure plot saved.")


# =========================================================
# Module 5: Pocket / cavity detection
# =========================================================

def detect_pockets(
    structure: PDB.Structure.Structure,
    outdir: str,
    prefix: str,
    probe_radius: float = 1.4,
    top_n: int = 5,
) -> pd.DataFrame:
    """
    Detect putative binding pockets using solvent-accessible surface area (SASA).
    Low-SASA residues buried in clusters are candidate pocket-lining residues.

    Strategy:
      1. Compute SASA for all residues (Shrake-Rupley algorithm via Biopython)
      2. Identify buried residues (SASA < 25 Å²)
      3. Cluster buried residues by spatial proximity (8 Å cutoff)
      4. Report top clusters as candidate pockets

    Returns DataFrame of pocket-lining residues.
    """
    model = structure[0]

    # Compute SASA
    sr = ShrakeRupley(probe_radius=probe_radius)
    sr.compute(structure, level="R")  # per-residue

    residues_data = []
    for chain in model.get_chains():
        for res in chain.get_residues():
            hetflag, resseq, _ = res.get_id()
            if hetflag.strip():
                continue
            if "CA" not in res:
                continue
            sasa = res.sasa if hasattr(res, "sasa") else 0.0
            coord = res["CA"].get_vector().get_array()
            residues_data.append({
                "chain":   chain.get_id(),
                "resseq":  resseq,
                "resname": res.get_resname(),
                "sasa":    float(sasa),
                "x":       coord[0],
                "y":       coord[1],
                "z":       coord[2],
            })

    res_df = pd.DataFrame(residues_data)
    if res_df.empty:
        print("[WARN] No residue data for pocket detection.")
        return pd.DataFrame(), pd.DataFrame()

    # Buried residues
    buried = res_df[res_df["sasa"] < 25.0].copy()
    if buried.empty:
        print("[WARN] No buried residues found (SASA < 25 Å²). Try increasing probe_radius.")
        return pd.DataFrame(), pd.DataFrame()

    # Spatial clustering — simple greedy approach
    coords = buried[["x", "y", "z"]].values
    assigned = np.full(len(coords), -1, dtype=int)
    cluster_id = 0

    for i in range(len(coords)):
        if assigned[i] != -1:
            continue
        assigned[i] = cluster_id
        for j in range(i + 1, len(coords)):
            if assigned[j] == -1:
                dist = np.linalg.norm(coords[i] - coords[j])
                if dist <= 8.0:
                    assigned[j] = cluster_id
        cluster_id += 1

    buried = buried.copy()
    buried["pocket_id"] = assigned

    pocket_summary = (
        buried.groupby("pocket_id")
        .agg(
            n_residues=("resseq", "count"),
            mean_sasa=("sasa", "mean"),
            centroid_x=("x", "mean"),
            centroid_y=("y", "mean"),
            centroid_z=("z", "mean"),
            chains=("chain", lambda x: ";".join(sorted(x.unique()))),
            residues=("resseq", lambda x: ";".join(str(r) for r in sorted(x))),
        )
        .reset_index()
        .sort_values("n_residues", ascending=False)
        .reset_index(drop=True)
    )
    pocket_summary["pocket_rank"] = range(1, len(pocket_summary) + 1)
    top_pockets = pocket_summary.head(top_n)

    # Save tables
    pocket_summary.to_csv(os.path.join(outdir, f"{prefix}.pockets.tsv"), sep="\t", index=False)
    buried.to_csv(os.path.join(outdir, f"{prefix}.pocket_residues.tsv"), sep="\t", index=False)
    res_df.to_csv(os.path.join(outdir, f"{prefix}.sasa_per_residue.tsv"), sep="\t", index=False)

    # Plot: SASA distribution + top pocket highlight
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4), dpi=300)

    # SASA distribution
    ax1.hist(res_df["sasa"], bins=40, color="#4E79A7", edgecolor="white", linewidth=0.4)
    ax1.axvline(25, color="#D62728", lw=1.5, ls="--", label="Buried threshold (25 Å²)")
    ax1.set_xlabel("SASA per residue (Å²)", fontsize=11)
    ax1.set_ylabel("Number of residues", fontsize=11)
    ax1.set_title("SASA Distribution", fontsize=12)
    ax1.legend(fontsize=9)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Top pockets bar
    palette = plt.cm.Set2(np.linspace(0, 0.7, len(top_pockets)))
    ax2.barh(
        [f"Pocket {r}" for r in top_pockets["pocket_rank"]],
        top_pockets["n_residues"],
        color=palette, edgecolor="white"
    )
    ax2.set_xlabel("Number of buried residues", fontsize=11)
    ax2.set_title(f"Top {len(top_pockets)} Candidate Pockets", fontsize=12)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.invert_yaxis()

    fig.suptitle(f"Pocket Detection  ·  {prefix}", fontsize=13)
    plt.tight_layout()
    save_fig(fig, os.path.join(outdir, f"{prefix}.pockets"))
    print(f"[INFO] Pocket analysis saved. Top pocket: {top_pockets.iloc[0]['n_residues']} residues.")

    return top_pockets, buried


# =========================================================
# Module 5b: Pocket surface 3D viewer
# =========================================================

# Distinct colours for up to 7 pockets (cycles if more)
_POCKET_COLORS = [
    "#FF4136", "#2ECC40", "#0074D9", "#FF851B",
    "#B10DC9", "#FFDC00", "#7FDBFF",
]


def _resi_string(reslist: List[int]) -> str:
    """Compress a list of residue numbers into a py3Dmol range string.

    e.g. [45, 46, 47, 80, 100] → '45-47,80,100'
    """
    sorted_res = sorted(set(reslist))
    ranges: List[str] = []
    start = end = sorted_res[0]
    for r in sorted_res[1:]:
        if r == end + 1:
            end = r
        else:
            ranges.append(f"{start}-{end}" if start != end else str(start))
            start = end = r
    ranges.append(f"{start}-{end}" if start != end else str(start))
    return ",".join(ranges)


def render_pocket_surface_viewer(
    pdb_path: str,
    pocket_summary: pd.DataFrame,
    buried_df: pd.DataFrame,
    prefix: str,
    outdir: str,
) -> str:
    """
    Build a self-contained HTML 3D viewer that highlights pocket-lining
    residues as coloured VDW surfaces on the full protein structure,
    with per-residue labels showing amino-acid type and sequence number.

    Layout:
      - Full protein: light-gray semi-transparent cartoon (backbone visible)
      - Faint white VDW surface over whole protein (shape context)
      - Each pocket: opaque coloured VDW surface + stick highlights
      - Per-residue label on Cα: "<3-letter code><resseq>" e.g. "HIS57"
        (white background, pocket-colour border, black text, shown in front)
      - Inline colour legend in the page footer

    Labels are anchored to Cα atoms so each residue gets exactly one label
    regardless of how many atoms it contains.

    Returns the path to the saved HTML file.
    """
    try:
        import py3Dmol  # type: ignore
    except ImportError:
        print("[WARN] py3Dmol not installed — skipping pocket surface viewer.")
        return ""

    with open(pdb_path, "r") as fh:
        pdb_data = fh.read()

    view = py3Dmol.view(width=900, height=650)
    view.addModel(pdb_data, "pdb")

    # Full-protein cartoon — drawn first so it is the background layer.
    # No global surface is added; surfaces are only drawn for pocket residues
    # so the backbone cartoon stays clearly visible everywhere else.
    view.setStyle({}, {"cartoon": {"color": "lightgray", "opacity": 0.70}})

    # Per-pocket coloured surfaces + labels
    legend_items: List[str] = []
    for _, row in pocket_summary.iterrows():
        rank = int(row["pocket_rank"])
        color = _POCKET_COLORS[(rank - 1) % len(_POCKET_COLORS)]
        pocket_id = int(row["pocket_id"])

        pocket_res = buried_df[buried_df["pocket_id"] == pocket_id]
        if pocket_res.empty:
            continue

        for chain_id, chain_grp in pocket_res.groupby("chain"):
            resi_str = _resi_string(chain_grp["resseq"].tolist())
            sel = {"chain": chain_id, "resi": resi_str}

            # ── coloured VDW surface (slightly transparent so sticks show) ─
            view.addSurface(py3Dmol.VDW, {"opacity": 0.70, "color": color}, sel)

            # ── stick representation for side-chain detail ────────────────
            view.addStyle(sel, {"stick": {"color": color, "radius": 0.16}})

            # ── residue labels anchored to Cα atoms ───────────────────────
            # addResLabels uses the PDB resname+resi automatically,
            # producing labels like "HIS57", "GLY102", etc.
            # We restrict to atom "CA" so only one label appears per residue.
            label_sel = {"chain": chain_id, "resi": resi_str, "atom": "CA"}
            view.addResLabels(
                label_sel,
                {
                    "font":             "Arial",
                    "fontSize":         11,
                    "fontColor":        "black",
                    "fontOpacity":      1.0,
                    "backgroundOpacity": 0.80,
                    "backgroundColor":  "white",
                    "borderThickness":  1.0,
                    "borderColor":      color,
                    "borderOpacity":    1.0,
                    "inFront":          True,
                    "showBackground":   True,
                },
            )

        legend_items.append(
            f'<span style="display:inline-flex;align-items:center;margin:0 10px;">'
            f'<span style="width:14px;height:14px;border-radius:50%;background:{color};'
            f'display:inline-block;margin-right:5px;border:1px solid #888;"></span>'
            f'Pocket {rank}&nbsp;({int(row["n_residues"])} res)</span>'
        )

    view.zoomTo()
    view.setBackgroundColor("white")

    # ── Assemble the page with an inline legend below the viewer ──────────
    viewer_html = view._make_html()

    legend_html = (
        '<div style="font-family:sans-serif;font-size:13px;padding:8px 4px;'
        'background:#f5f5f5;border-top:1px solid #ddd;text-align:center;">'
        "<strong>Pocket legend</strong>&nbsp;(labels: amino-acid type + residue number on C\u03b1):&nbsp;"
        + "".join(legend_items)
        + "</div>"
    )

    # Insert legend just before </body>
    full_html = viewer_html.replace("</body>", legend_html + "\n</body>")

    out_path = os.path.join(outdir, f"{prefix}.pocket_surface.html")
    with open(out_path, "w") as fh:
        fh.write(full_html)

    print(f"[INFO] Pocket surface viewer saved → {out_path}")
    return out_path


# =========================================================
# Module 6: Evolutionary conservation
# =========================================================

EBI_HMMER   = "https://www.ebi.ac.uk/Tools/hmmer/search/jackhmmer"
UNIPROT_API = "https://rest.uniprot.org/uniprotkb"

# Conservation colour gradient: variable (blue) → conserved (maroon)
# 9 grades as used by ConSurf (1 = most variable, 9 = most conserved)
CONSERVATION_GRADIENT = [
    "#009EC8",  # grade 1 — most variable (cyan-blue)
    "#00ADEF",
    "#53C8F5",
    "#ABE0F9",
    "#F0F0F0",  # grade 5 — intermediate (white/gray)
    "#FAC8AF",
    "#F5946A",
    "#EF5B2B",
    "#8B0000",  # grade 9 — most conserved (dark red)
]


def _fetch_uniprot_sequence(uniprot_acc: str) -> str:
    """Return the canonical FASTA sequence for a UniProt accession."""
    url = f"{UNIPROT_API}/{uniprot_acc}.fasta"
    r = requests.get(url, headers=HEADERS, timeout=30)
    r.raise_for_status()
    lines = r.text.strip().splitlines()
    return "".join(l for l in lines if not l.startswith(">"))


def _extract_pdb_sequence(structure: PDB.Structure.Structure, chain_id: str) -> Tuple[List[int], str]:
    """Extract (resseq_list, one-letter sequence) from a PDB chain."""
    three2one = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }
    resseqs: List[int] = []
    seq: List[str] = []
    chain = structure[0][chain_id]
    for res in chain.get_residues():
        hetflag, resseq, _ = res.get_id()
        if hetflag.strip():
            continue
        aa = three2one.get(res.get_resname(), "X")
        resseqs.append(resseq)
        seq.append(aa)
    return resseqs, "".join(seq)


def _submit_jackhmmer(sequence: str, n_iterations: int = 1, evalue: float = 1e-4) -> str:
    """
    Submit a jackhmmer search to EBI HMMER and return the job UUID.
    Searches against UniRef90 by default.
    """
    payload = {
        "seqdb":      "uniref90",
        "seq":        f">query\n{sequence}",
        "iterations": n_iterations,
        "E":          evalue,
        "incE":       evalue,
    }
    r = requests.post(EBI_HMMER, data=payload, headers=HEADERS, timeout=60,
                      allow_redirects=False)
    # EBI returns 303 redirect with UUID in Location header
    location = r.headers.get("Location", "")
    if not location:
        raise RuntimeError(f"EBI HMMER submission failed (status {r.status_code}).")
    uuid = location.rstrip("/").split("/")[-1]
    print(f"[INFO] EBI HMMER job submitted: {uuid}")
    return uuid


def _poll_jackhmmer(uuid: str, poll_interval: int = 5, max_wait: int = 300) -> str:
    """Poll EBI HMMER until the job completes; return Stockholm MSA text."""
    import time
    status_url = f"https://www.ebi.ac.uk/Tools/hmmer/results/{uuid}/score"
    sto_url    = f"https://www.ebi.ac.uk/Tools/hmmer/download/{uuid}/score?format=sto&filename=msa.sto"

    waited = 0
    while waited < max_wait:
        try:
            r = requests.get(status_url, headers=HEADERS, timeout=30)
            if r.status_code == 200:
                break
        except Exception:
            pass
        time.sleep(poll_interval)
        waited += poll_interval
        print(f"[INFO] Waiting for HMMER results… ({waited}s)")

    if waited >= max_wait:
        raise RuntimeError("EBI HMMER job timed out after {max_wait}s.")

    r = requests.get(sto_url, headers=HEADERS, timeout=60)
    r.raise_for_status()
    return r.text


def _parse_stockholm_msa(sto_text: str) -> Dict[str, str]:
    """
    Parse a Stockholm-format MSA.
    Returns {seq_name: aligned_sequence} for all sequences.
    """
    seqs: Dict[str, List[str]] = {}
    for line in sto_text.splitlines():
        if line.startswith("#") or line.startswith("//") or not line.strip():
            continue
        parts = line.split()
        if len(parts) == 2:
            name, seq = parts
            seqs.setdefault(name, []).append(seq)
    return {name: "".join(frags) for name, frags in seqs.items()}


def _compute_conservation_entropy(msa: Dict[str, str]) -> List[float]:
    """
    Compute per-column conservation scores as 1 - normalised Shannon entropy.
    Returns a list of scores in [0, 1] (1 = fully conserved, 0 = maximally variable).
    Gaps are treated as a 21st character.
    """
    sequences = list(msa.values())
    if not sequences:
        return []
    ncols = len(sequences[0])
    scores: List[float] = []
    for col in range(ncols):
        chars = [s[col].upper() for s in sequences if col < len(s)]
        total = len(chars)
        if total == 0:
            scores.append(0.0)
            continue
        counts: Dict[str, int] = {}
        for c in chars:
            counts[c] = counts.get(c, 0) + 1
        entropy = 0.0
        for count in counts.values():
            p = count / total
            if p > 0:
                entropy -= p * np.log2(p)
        max_entropy = np.log2(21)  # 20 aa + gap
        scores.append(1.0 - entropy / max_entropy if max_entropy > 0 else 1.0)
    return scores


def _score_to_grade(score: float) -> int:
    """Map a conservation score [0, 1] to a ConSurf-style grade 1–9."""
    return max(1, min(9, int(score * 8) + 1))


def fetch_conservation_scores(
    structure: PDB.Structure.Structure,
    chain_id: str,
    outdir: str,
    prefix: str,
    uniprot_acc: Optional[str] = None,
    n_iterations: int = 1,
) -> pd.DataFrame:
    """
    Compute per-residue evolutionary conservation scores using EBI HMMER jackhmmer.

    Strategy:
      1. If uniprot_acc provided, use the canonical UniProt sequence
         (more reliable than PDB SEQRES for MSA search).
         Otherwise extract the sequence directly from the PDB structure.
      2. Submit to EBI HMMER jackhmmer against UniRef90
      3. Parse Stockholm MSA → compute column-wise conservation entropy
      4. Map scores back onto PDB resseq numbers

    Returns a DataFrame with columns: chain, resseq, resname, conservation_score,
    conservation_grade (1–9), conservation_color.
    """
    resseqs, pdb_seq = _extract_pdb_sequence(structure, chain_id)

    if uniprot_acc:
        try:
            query_seq = _fetch_uniprot_sequence(uniprot_acc)
            print(f"[INFO] Using UniProt sequence ({len(query_seq)} aa) for MSA search.")
        except Exception as e:
            print(f"[WARN] UniProt fetch failed ({e}); using PDB sequence instead.")
            query_seq = pdb_seq
    else:
        query_seq = pdb_seq
        print(f"[INFO] Using PDB chain {chain_id} sequence ({len(pdb_seq)} aa) for MSA search.")

    # Submit and poll HMMER
    try:
        uuid = _submit_jackhmmer(query_seq, n_iterations=n_iterations)
        sto_text = _poll_jackhmmer(uuid)
    except Exception as e:
        print(f"[WARN] EBI HMMER failed: {e}. Returning uniform conservation scores.")
        # Return a uniform "unknown" score so the module doesn't crash
        rows = [
            {"chain": chain_id, "resseq": r, "resname": "UNK",
             "conservation_score": 0.5, "conservation_grade": 5,
             "conservation_color": CONSERVATION_GRADIENT[4]}
            for r in resseqs
        ]
        return pd.DataFrame(rows)

    # Parse MSA
    msa = _parse_stockholm_msa(sto_text)
    if not msa:
        print("[WARN] Empty MSA returned — conservation scores unavailable.")
        return pd.DataFrame()

    print(f"[INFO] MSA contains {len(msa)} sequences.")

    # Per-column conservation
    col_scores = _compute_conservation_entropy(msa)

    # The query sequence is the first entry in the MSA; align column scores
    # back to PDB resseq positions, skipping gap columns in the query.
    query_name = next(iter(msa))
    query_aligned = msa[query_name]

    # Map aligned positions → PDB resseq
    pdb_pos_iter = iter(range(len(resseqs)))
    conservation_per_residue: List[float] = []
    col_idx = 0
    for aa_aligned, score in zip(query_aligned, col_scores):
        if aa_aligned == "-":
            col_idx += 1
            continue
        conservation_per_residue.append(score)
        col_idx += 1

    # Truncate or pad to match PDB residue count
    n = min(len(resseqs), len(conservation_per_residue))

    # Build per-residue table
    chain_obj = structure[0][chain_id]
    res_list = [res for res in chain_obj.get_residues()
                if not res.get_id()[0].strip()]

    rows = []
    for i in range(n):
        score = conservation_per_residue[i]
        grade = _score_to_grade(score)
        color = CONSERVATION_GRADIENT[grade - 1]
        resname = res_list[i].get_resname() if i < len(res_list) else "UNK"
        rows.append({
            "chain":              chain_id,
            "resseq":             resseqs[i],
            "resname":            resname,
            "conservation_score": round(score, 4),
            "conservation_grade": grade,
            "conservation_color": color,
        })

    df = pd.DataFrame(rows)
    out_path = os.path.join(outdir, f"{prefix}.conservation_scores.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"[INFO] Conservation scores saved → {out_path}")
    return df


def plot_conservation(
    cons_df: pd.DataFrame,
    prefix: str,
    outdir: str,
) -> None:
    """Bar plot of per-residue conservation score, coloured by grade."""
    if cons_df.empty:
        return

    fig, ax = plt.subplots(figsize=(14, 3.5), dpi=300)
    x = np.arange(len(cons_df))
    bar_colors = cons_df["conservation_color"].tolist()
    ax.bar(x, cons_df["conservation_score"], color=bar_colors, width=1.0, linewidth=0)

    step = max(1, len(cons_df) // 20)
    ax.set_xticks(x[::step])
    ax.set_xticklabels(cons_df["resseq"].values[::step], fontsize=8, rotation=45, ha="right")
    ax.set_ylabel("Conservation score\n(1 = conserved, 0 = variable)", fontsize=10)
    ax.set_xlabel("Residue number", fontsize=11)
    ax.set_ylim(0, 1.05)
    ax.set_title(f"Evolutionary conservation  ·  {prefix}", fontsize=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Gradient legend
    from matplotlib.colors import LinearSegmentedColormap
    import matplotlib.cm as cm
    cmap = LinearSegmentedColormap.from_list("cons", CONSERVATION_GRADIENT, N=256)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(1, 9))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation="vertical", pad=0.01, shrink=0.85)
    cbar.set_label("ConSurf-style grade\n(1=variable, 9=conserved)", fontsize=9)
    cbar.set_ticks([1, 3, 5, 7, 9])

    plt.tight_layout()
    save_fig(fig, os.path.join(outdir, f"{prefix}.conservation"))
    print("[INFO] Conservation plot saved.")


def render_conservation_viewer(
    pdb_path: str,
    cons_df: pd.DataFrame,
    prefix: str,
    outdir: str,
) -> str:
    """
    Build a self-contained HTML 3D viewer coloured by conservation grade.
    Variable residues = cyan-blue; conserved residues = dark red.
    """
    if cons_df.empty:
        return ""

    try:
        import py3Dmol  # type: ignore
    except ImportError:
        print("[WARN] py3Dmol not installed — skipping conservation viewer.")
        return ""

    with open(pdb_path) as fh:
        pdb_data = fh.read()

    view = py3Dmol.view(width=900, height=650)
    view.addModel(pdb_data, "pdb")

    # Default everything to gray first, then colour per residue
    view.setStyle({}, {"cartoon": {"color": "#DDDDDD", "opacity": 0.3}})

    chain_id = cons_df["chain"].iloc[0]
    for _, row in cons_df.iterrows():
        color = row["conservation_color"]
        sel = {"chain": chain_id, "resi": str(int(row["resseq"]))}
        view.setStyle(sel, {"cartoon": {"color": color, "opacity": 0.95}})

    view.zoomTo()
    view.setBackgroundColor("white")

    viewer_html = view._make_html()

    # Gradient legend bar in HTML
    grad_css = ", ".join(CONSERVATION_GRADIENT)
    legend_html = (
        '<div style="font-family:sans-serif;font-size:12px;padding:8px 16px;'
        'background:#f5f5f5;border-top:1px solid #ddd;display:flex;'
        'align-items:center;gap:12px;">'
        '<span style="white-space:nowrap;">Variable (grade 1)</span>'
        f'<div style="flex:1;height:14px;border-radius:4px;'
        f'background:linear-gradient(to right,{grad_css});'
        f'border:1px solid #ccc;"></div>'
        '<span style="white-space:nowrap;">Conserved (grade 9)</span>'
        '</div>'
    )

    full_html = viewer_html.replace("</body>", legend_html + "\n</body>")
    out_path = os.path.join(outdir, f"{prefix}.conservation.html")
    with open(out_path, "w") as fh:
        fh.write(full_html)

    print(f"[INFO] Conservation viewer saved → {out_path}")
    return out_path


# =========================================================
# Module 7: Protein-protein interaction network (STRING)
# =========================================================

def fetch_ppi_string(
    gene_symbol: str,
    species: int = 9606,
    score_threshold: int = 400,
    limit: int = 50,
) -> pd.DataFrame:
    """
    Fetch interaction partners from the STRING database REST API.
    Returns DataFrame with interaction data.
    """
    # Map gene symbol to STRING ID
    url_map = f"{STRING_API}/json/get_string_ids"
    r = requests.post(url_map, data={
        "identifiers": gene_symbol,
        "species":     species,
        "limit":       1,
        "echo_query":  1,
    }, headers=HEADERS, timeout=30)
    r.raise_for_status()
    ids = r.json()
    if not ids:
        raise ValueError(f"Gene '{gene_symbol}' not found in STRING (species {species}).")

    string_id = ids[0]["stringId"]
    print(f"[INFO] STRING ID: {string_id}")

    # Fetch interactions
    url_int = f"{STRING_API}/json/network"
    r2 = requests.post(url_int, data={
        "identifiers":      string_id,
        "species":          species,
        "required_score":   score_threshold,
        "limit":            limit,
        "add_nodes":        limit,
        "network_type":     "functional",
    }, headers=HEADERS, timeout=60)
    r2.raise_for_status()
    interactions = r2.json()

    if not interactions:
        print(f"[WARN] No interactions found for {gene_symbol} above score {score_threshold}.")
        return pd.DataFrame()

    rows = []
    for item in interactions:
        rows.append({
            "proteinA":      item.get("preferredName_A", item.get("stringId_A", "")),
            "proteinB":      item.get("preferredName_B", item.get("stringId_B", "")),
            "score":         item.get("score", 0),
            "score_coexp":   item.get("escore", 0),
            "score_exp":     item.get("ascore", 0),
            "score_db":      item.get("dscore", 0),
            "score_textmine": item.get("tscore", 0),
        })

    return pd.DataFrame(rows).drop_duplicates().reset_index(drop=True)


def plot_ppi_network(
    ppi_df: pd.DataFrame,
    gene_symbol: str,
    outdir: str,
    prefix: str,
    top_n: int = 30,
) -> None:
    """Draw PPI network with networkx + matplotlib."""
    try:
        import networkx as nx
    except ImportError:
        raise ImportError("networkx is required for PPI network plots. Install with: pip install networkx")

    if ppi_df.empty:
        return

    plot_df = ppi_df.head(top_n)
    G = nx.Graph()

    for _, row in plot_df.iterrows():
        a, b = row["proteinA"], row["proteinB"]
        G.add_edge(a, b, weight=float(row["score"]))

    if len(G.nodes()) == 0:
        return

    # Node sizes: degree-proportional; query gene highlighted
    degrees = dict(G.degree())
    max_deg = max(degrees.values()) if degrees else 1
    node_sizes = [
        800 if n == gene_symbol.upper() else 200 + 400 * (degrees[n] / max_deg)
        for n in G.nodes()
    ]
    node_colors = [
        "#D62728" if n == gene_symbol.upper() else "#4E79A7"
        for n in G.nodes()
    ]

    # Edge weights → widths
    weights = [G[u][v].get("weight", 400) / 1000 * 3 for u, v in G.edges()]

    pos = nx.spring_layout(G, seed=42, k=2.0 / max(1, len(G.nodes()) ** 0.5))

    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
    nx.draw_networkx_edges(G, pos, ax=ax, width=weights, alpha=0.4, edge_color="#AAAAAA")
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=node_sizes,
                           node_color=node_colors, alpha=0.9)
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=7, font_color="white",
                            font_weight="bold",
                            labels={n: n for n in G.nodes()})

    legend_handles = [
        mpatches.Patch(color="#D62728", label=f"Query: {gene_symbol}"),
        mpatches.Patch(color="#4E79A7", label="Interaction partners"),
    ]
    ax.legend(handles=legend_handles, fontsize=9, loc="upper left")
    ax.set_title(f"PPI Network  ·  {gene_symbol}  (STRING score ≥ {ppi_df['score'].min():.0f})",
                 fontsize=13, pad=10)
    ax.axis("off")
    plt.tight_layout()
    save_fig(fig, os.path.join(outdir, f"{prefix}.ppi_network"))
    print(f"[INFO] PPI network saved. Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")


# =========================================================
# Main
# =========================================================

def main():
    parser = argparse.ArgumentParser(
        description="Protein structure visualization: 3D viewer, contact map, pockets, PPI.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Structure source (mutually exclusive)
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--pdb-id", help="RCSB PDB ID to fetch (e.g. 4ZJH, 1TUP)")
    src.add_argument("--uniprot", help="UniProt accession for AlphaFold structure (e.g. P04637)")
    src.add_argument("--pdb-file", help="Local PDB file path")

    # Module selection
    parser.add_argument(
        "--modules", default="all",
        help=(
            "Comma-separated list: view, contact_map, bfactor, secondary, pocket, ppi. "
            "Use 'all' to run everything."
        ),
    )

    # Zoom / region
    parser.add_argument("--chain", default=None,
        help="Chain ID to focus on (e.g. A). Default: first chain.")
    parser.add_argument("--zoom-start", type=int, default=None,
        help="Residue number start for zoomed view / contact map region.")
    parser.add_argument("--zoom-end", type=int, default=None,
        help="Residue number end for zoomed view / contact map region.")

    # Viewer style
    parser.add_argument(
        "--color-scheme",
        choices=["spectrum", "chain", "bfactor"],
        default="spectrum",
        help="3D viewer colour scheme. Use 'bfactor' for AlphaFold pLDDT colouring.",
    )

    # Contact map
    parser.add_argument("--contact-threshold", type=float, default=8.0,
        help="Cα–Cα distance threshold (Å) for contact map.")

    # Pocket detection
    parser.add_argument("--top-pockets", type=int, default=5,
        help="Number of top candidate pockets to report.")
    parser.add_argument("--probe-radius", type=float, default=1.4,
        help="Probe radius (Å) for SASA calculation.")

    # Conservation
    parser.add_argument("--uniprot-for-conservation", metavar="ACC", default=None,
        help="UniProt accession used to fetch the canonical sequence for MSA-based "
             "conservation scoring (conservation module). If omitted, the PDB sequence "
             "is used directly.")
    parser.add_argument("--hmmer-iterations", type=int, default=1,
        help="Number of jackhmmer iterations for conservation MSA (default: 1).")

    # PPI
    parser.add_argument("--gene", default=None,
        help="Gene symbol for PPI search (required for ppi module).")
    parser.add_argument("--ppi-species", type=int, default=9606,
        help="NCBI taxon ID for STRING PPI search. Default: 9606 (human).")
    parser.add_argument("--ppi-score", type=int, default=400,
        help="Minimum combined STRING score (0–1000). Default: 400 (medium confidence).")
    parser.add_argument("--ppi-limit", type=int, default=50,
        help="Maximum number of STRING interaction partners to retrieve.")

    parser.add_argument("--outdir", required=True, help="Output directory.")

    args = parser.parse_args()
    ensure_dir(args.outdir)

    # Parse modules
    all_modules = {"view", "contact_map", "bfactor", "secondary", "pocket", "ppi", "conservation"}
    if args.modules.strip().lower() == "all":
        modules = all_modules.copy()
    else:
        modules = {m.strip().lower() for m in args.modules.split(",")}
        invalid = modules - all_modules
        if invalid:
            raise ValueError(f"Unknown module(s): {invalid}. Choose from: {all_modules}")

    # ----------------------------------------------------------
    # Determine prefix and structure source
    # ----------------------------------------------------------
    is_alphafold = False

    if args.pdb_id:
        prefix = args.pdb_id.upper()
        pdb_path = fetch_pdb_structure(args.pdb_id, args.outdir)
    elif args.uniprot:
        prefix = f"AF_{args.uniprot}"
        pdb_path = fetch_alphafold_structure(args.uniprot, args.outdir)
        is_alphafold = True
        if args.color_scheme == "spectrum":
            args.color_scheme = "bfactor"   # default to pLDDT colours for AF
    else:
        pdb_path = args.pdb_file
        prefix = os.path.splitext(os.path.basename(pdb_path))[0]

    # Parse structure
    structure_id = prefix
    print(f"[INFO] Parsing structure: {pdb_path}")
    structure = parse_structure(pdb_path, structure_id)

    # Determine chain
    model = structure[0]
    available_chains = [c.get_id() for c in model.get_chains()]
    chain_id = args.chain if args.chain and args.chain in available_chains else available_chains[0]
    if args.chain and args.chain not in available_chains:
        print(f"[WARN] Chain '{args.chain}' not found. Using chain '{chain_id}'.")
    print(f"[INFO] Using chain: {chain_id}  |  Available: {available_chains}")

    # Extract residue data (used by multiple modules)
    res_df = get_residue_data(structure, chain_id=None)   # all chains for pocket
    print(f"[INFO] Total residues: {len(res_df)}")

    # ----------------------------------------------------------
    # Run modules
    # ----------------------------------------------------------
    summary: Dict = {
        "structure_id": prefix,
        "pdb_path":     pdb_path,
        "is_alphafold": is_alphafold,
        "chain_used":   chain_id,
        "n_residues":   len(res_df),
        "modules":      args.modules,
    }

    if "view" in modules:
        print(f"[INFO] Generating HTML 3D viewer")
        render_html_viewer(
            pdb_path, prefix, args.outdir, prefix,
            zoom_chain=chain_id if (args.zoom_start or args.zoom_end) else (chain_id if args.chain else None),
            zoom_start=args.zoom_start,
            zoom_end=args.zoom_end,
            color_scheme=args.color_scheme,
            is_alphafold=is_alphafold,
        )

    if "contact_map" in modules:
        print(f"[INFO] Computing contact map")
        plot_contact_map(
            res_df, args.outdir, prefix,
            threshold_A=args.contact_threshold,
            chain_id=chain_id,
            zoom_start=args.zoom_start,
            zoom_end=args.zoom_end,
        )

    if "bfactor" in modules:
        print(f"[INFO] Plotting B-factor / pLDDT")
        plot_bfactor(
            res_df, args.outdir, prefix,
            is_alphafold=is_alphafold,
            chain_id=chain_id,
            zoom_start=args.zoom_start,
            zoom_end=args.zoom_end,
        )

    if "secondary" in modules:
        print(f"[INFO] Analysing secondary structure")
        plot_secondary_structure(
            structure, pdb_path, args.outdir, prefix,
            chain_id=chain_id if args.chain else None,
        )

    if "pocket" in modules:
        print(f"[INFO] Detecting binding pockets")
        pocket_df, buried_df = detect_pockets(
            structure, args.outdir, prefix,
            probe_radius=args.probe_radius,
            top_n=args.top_pockets,
        )
        if not pocket_df.empty:
            summary["n_pockets_found"] = len(pocket_df)
            summary["top_pocket_residues"] = int(pocket_df.iloc[0]["n_residues"])
            render_pocket_surface_viewer(
                pdb_path, pocket_df, buried_df, prefix, args.outdir
            )

    if "conservation" in modules:
        print(f"[INFO] Computing evolutionary conservation (EBI HMMER jackhmmer)")
        uniprot_for_cons = getattr(args, "uniprot_for_conservation", None)
        # If the user ran with --uniprot, reuse that accession automatically
        if not uniprot_for_cons and args.uniprot:
            uniprot_for_cons = args.uniprot
        try:
            cons_df = fetch_conservation_scores(
                structure, chain_id, args.outdir, prefix,
                uniprot_acc=uniprot_for_cons,
                n_iterations=args.hmmer_iterations,
            )
            if not cons_df.empty:
                plot_conservation(cons_df, prefix, args.outdir)
                render_conservation_viewer(pdb_path, cons_df, prefix, args.outdir)
                summary["conservation_residues"] = len(cons_df)
                summary["mean_conservation"] = round(
                    float(cons_df["conservation_score"].mean()), 4
                )
        except Exception as ex:
            print(f"[WARN] Conservation module failed: {ex}")

    if "ppi" in modules:
        gene = args.gene
        if not gene:
            print("[WARN] --gene is required for the ppi module. Skipping PPI.")
        else:
            print(f"[INFO] Fetching PPI from STRING for {gene}")
            try:
                ppi_df = fetch_ppi_string(
                    gene, args.ppi_species, args.ppi_score, args.ppi_limit
                )
                if not ppi_df.empty:
                    ppi_df.to_csv(
                        os.path.join(args.outdir, f"{prefix}.ppi_interactions.tsv"),
                        sep="\t", index=False
                    )
                    plot_ppi_network(ppi_df, gene, args.outdir, prefix)
                    summary["n_ppi_partners"] = ppi_df[["proteinA", "proteinB"]].values.ravel()
                    unique_partners = set(ppi_df["proteinA"].tolist() + ppi_df["proteinB"].tolist())
                    unique_partners.discard(gene.upper())
                    summary["n_ppi_partners"] = len(unique_partners)
            except Exception as ex:
                print(f"[WARN] PPI fetch failed: {ex}")

    # ----------------------------------------------------------
    # Summary
    # ----------------------------------------------------------
    pd.DataFrame([summary]).to_csv(
        os.path.join(args.outdir, "structure_summary.tsv"), sep="\t", index=False
    )
    with open(os.path.join(args.outdir, "summary.txt"), "w") as f:
        for k, v in summary.items():
            f.write(f"{k}: {v}\n")

    print(f"[DONE] Results written to: {args.outdir}")


if __name__ == "__main__":
    main()
