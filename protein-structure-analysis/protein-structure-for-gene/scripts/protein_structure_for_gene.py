#!/usr/bin/env python3
"""
Protein Structure Analysis for a Single Gene

For a given gene symbol, this script:
  1. Resolves the gene → UniProt canonical accession
  2. Fetches domain/feature annotations from UniProt
  3. Fetches experimental PDB structures from RCSB
  4. Fetches the AlphaFold2 predicted structure entry from EBI
  5. Draws a publication-quality domain map (linear schematic)
  6. Writes structured TSV outputs and a plain-text summary

Data sources (all public, no authentication required):
  - UniProt REST API   https://rest.uniprot.org
  - RCSB PDB REST API https://search.rcsb.org / https://data.rcsb.org
  - AlphaFold DB API  https://alphafold.ebi.ac.uk/api

Dependencies: matplotlib, pandas, requests
"""

import argparse
import os
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import requests

# =========================================================
# Constants
# =========================================================

UNIPROT_BASE   = "https://rest.uniprot.org"
RCSB_SEARCH    = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_GRAPHQL   = "https://data.rcsb.org/graphql"
ALPHAFOLD_BASE = "https://alphafold.ebi.ac.uk/api"

FEATURE_COLORS = {
    "Domain":           "#4E79A7",
    "Region":           "#F28E2B",
    "Motif":            "#59A14F",
    "Binding site":     "#E15759",
    "Active site":      "#B07AA1",
    "Signal":           "#76B7B2",
    "Transmembrane":    "#FF9DA7",
    "Coiled coil":      "#9C755F",
    "Compositional bias": "#BAB0AC",
    "Other":            "#D3D3D3",
}

DEFAULT_FEATURE_TYPES = {
    "Domain", "Region", "Motif", "Binding site", "Active site",
    "Signal", "Transmembrane", "Coiled coil",
}

HEADERS = {"User-Agent": "protein-structure-for-gene/0.1 (bioinfor-claw)"}

# =========================================================
# Utilities
# =========================================================

def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def _get(url: str, params: Optional[dict] = None, timeout: int = 30) -> requests.Response:
    r = requests.get(url, params=params, headers=HEADERS, timeout=timeout)
    r.raise_for_status()
    return r


def _post(url: str, payload: dict, timeout: int = 30) -> requests.Response:
    r = requests.post(url, json=payload, headers={**HEADERS, "Content-Type": "application/json"}, timeout=timeout)
    r.raise_for_status()
    return r


# =========================================================
# UniProt
# =========================================================

def resolve_uniprot(gene_symbol: str, organism: str = "human") -> Tuple[str, str, int]:
    """
    Map gene symbol → canonical UniProt accession.

    Returns (accession, entry_name, sequence_length).
    Prefers reviewed (Swiss-Prot) entries for human / mouse / rat.
    """
    taxon_map = {
        "human": "9606",
        "mouse": "10090",
        "rat":   "10116",
    }
    taxon = taxon_map.get(organism.lower(), organism)

    query = f"gene_exact:{gene_symbol} AND organism_id:{taxon} AND reviewed:true"
    params = {
        "query": query,
        "fields": "accession,entry_name,gene_names,protein_name,length,organism_name",
        "format": "json",
        "size": 5,
    }
    r = _get(f"{UNIPROT_BASE}/uniprotkb/search", params=params)
    results = r.json().get("results", [])

    if not results:
        # Relax to unreviewed
        query2 = f"gene_exact:{gene_symbol} AND organism_id:{taxon}"
        params["query"] = query2
        r2 = _get(f"{UNIPROT_BASE}/uniprotkb/search", params=params)
        results = r2.json().get("results", [])

    if not results:
        raise ValueError(
            f"Could not resolve gene '{gene_symbol}' (organism: {organism}) in UniProt."
        )

    entry = results[0]
    accession   = entry["primaryAccession"]
    entry_name  = entry.get("uniProtkbId", accession)
    seq_len     = entry.get("sequence", {}).get("length", 0)
    return accession, entry_name, int(seq_len)


def fetch_uniprot_features(accession: str) -> List[Dict]:
    """
    Fetch annotated features from UniProt for a given accession.
    Returns list of dicts with keys: type, description, start, end.
    """
    r = _get(f"{UNIPROT_BASE}/uniprotkb/{accession}", params={"format": "json"})
    data = r.json()
    features_raw = data.get("features", [])

    features = []
    for f in features_raw:
        ftype = f.get("type", "Other")
        desc  = f.get("description", "")
        loc   = f.get("location", {})
        start = loc.get("start", {}).get("value")
        end   = loc.get("end", {}).get("value")
        if start is None or end is None:
            continue
        features.append(
            {
                "type":        ftype,
                "description": desc or ftype,
                "start":       int(start),
                "end":         int(end),
                "length":      int(end) - int(start) + 1,
            }
        )

    return features


def fetch_uniprot_summary(accession: str) -> Dict:
    """Fetch basic protein metadata from UniProt."""
    r = _get(f"{UNIPROT_BASE}/uniprotkb/{accession}", params={"format": "json"})
    data = r.json()

    protein_names = data.get("proteinDescription", {})
    rec_name = protein_names.get("recommendedName", {})
    full_name = rec_name.get("fullName", {}).get("value", "")

    genes = data.get("genes", [])
    primary_gene = ""
    if genes:
        gnames = genes[0].get("geneName", {})
        primary_gene = gnames.get("value", "")

    organism = data.get("organism", {}).get("scientificName", "")
    seq_len  = data.get("sequence", {}).get("length", 0)
    reviewed = data.get("entryType", "") == "UniProtKB reviewed (Swiss-Prot)"

    return {
        "accession":   accession,
        "protein_name": full_name,
        "gene_name":   primary_gene,
        "organism":    organism,
        "seq_length":  seq_len,
        "reviewed":    reviewed,
        "uniprot_url": f"https://www.uniprot.org/uniprot/{accession}",
    }


# =========================================================
# RCSB PDB
# =========================================================

def fetch_pdb_structures(accession: str, max_structures: int = 20) -> pd.DataFrame:
    """
    Search RCSB PDB for structures linked to a UniProt accession.
    Returns DataFrame with PDB-level metadata.
    """
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator":  "exact_match",
                "value":     accession,
            },
        },
        "request_options": {
            "paginate": {"start": 0, "rows": max_structures},
            "sort":     [{"sort_by": "score", "direction": "desc"}],
        },
        "return_type": "entry",
    }

    try:
        r = _post(RCSB_SEARCH, query, timeout=30)
        result_set = r.json().get("result_set", [])
        pdb_ids = [item["identifier"] for item in result_set if "identifier" in item]
    except Exception:
        return pd.DataFrame()

    if not pdb_ids:
        return pd.DataFrame()

    # Fetch metadata via GraphQL
    rows = []
    for pdb_id in pdb_ids:
        try:
            gql = {
                "query": """
                query($id: String!) {
                  entry(entry_id: $id) {
                    rcsb_id
                    struct { title }
                    rcsb_entry_info {
                      resolution_combined
                      experimental_method
                      deposited_atom_count
                      polymer_entity_count_protein
                    }
                    rcsb_entry_container_identifiers { pubmed_id }
                    audit_author { name }
                    pdbx_audit_revision_history { revision_date }
                  }
                }
                """,
                "variables": {"id": pdb_id.upper()},
            }
            gr = _post(RCSB_GRAPHQL, gql, timeout=20)
            e = gr.json().get("data", {}).get("entry", {}) or {}

            ei     = e.get("rcsb_entry_info") or {}
            auths  = e.get("audit_author") or []
            revs   = e.get("pdbx_audit_revision_history") or []
            conts  = e.get("rcsb_entry_container_identifiers") or {}

            author_str = "; ".join(a.get("name", "") for a in auths[:3])
            if len(auths) > 3:
                author_str += f" et al. (+{len(auths)-3})"

            dep_date = ""
            if revs:
                dep_date = revs[0].get("revision_date", "")

            res = ei.get("resolution_combined")
            resolution = round(float(res[0]), 2) if res else None

            rows.append(
                {
                    "pdb_id":         pdb_id.upper(),
                    "title":          (e.get("struct") or {}).get("title", ""),
                    "method":         ei.get("experimental_method", ""),
                    "resolution_A":   resolution,
                    "n_protein_chains": ei.get("polymer_entity_count_protein"),
                    "n_atoms":        ei.get("deposited_atom_count"),
                    "pubmed_id":      conts.get("pubmed_id"),
                    "authors":        author_str,
                    "deposition_date": dep_date,
                    "rcsb_url":       f"https://www.rcsb.org/structure/{pdb_id.upper()}",
                }
            )
        except Exception:
            continue

    return pd.DataFrame(rows)


# =========================================================
# AlphaFold
# =========================================================

def fetch_alphafold_entry(accession: str) -> Optional[Dict]:
    """
    Fetch AlphaFold2 prediction entry from EBI for a UniProt accession.
    Returns dict with metadata, or None if not available.
    """
    url = f"{ALPHAFOLD_BASE}/prediction/{accession}"
    try:
        r = _get(url, timeout=20)
        entries = r.json()
        if not entries:
            return None
        e = entries[0]
        return {
            "alphafold_id":      e.get("entryId"),
            "uniprot_accession": e.get("uniprotAccession"),
            "gene":              e.get("gene"),
            "organism":          e.get("organismScientificName"),
            "seq_length":        e.get("sequenceLength"),
            "model_created":     e.get("modelCreatedDate"),
            "latest_version":    e.get("latestVersion"),
            "pdb_url":           e.get("pdbUrl"),
            "cif_url":           e.get("cifUrl"),
            "pae_image_url":     e.get("paeImageUrl"),
            "alphafold_page":    f"https://alphafold.ebi.ac.uk/entry/{accession}",
        }
    except requests.HTTPError as e:
        if e.response is not None and e.response.status_code == 404:
            return None
        raise


# =========================================================
# Domain map figure
# =========================================================

def _feature_color(ftype: str) -> str:
    for key in FEATURE_COLORS:
        if ftype.lower().startswith(key.lower()):
            return FEATURE_COLORS[key]
    return FEATURE_COLORS["Other"]


def plot_domain_map(
    features: List[Dict],
    seq_length: int,
    gene: str,
    accession: str,
    outdir: str,
    feature_types: Optional[set] = None,
) -> str:
    """
    Draw a linear domain/feature map for the protein.
    Saves PNG + PDF. Returns path to PNG.
    """
    if feature_types is None:
        feature_types = DEFAULT_FEATURE_TYPES

    # Filter to displayable features
    plot_feats = [
        f for f in features
        if any(f["type"].lower().startswith(t.lower()) for t in feature_types)
    ]

    if not plot_feats:
        # Draw an empty backbone if nothing to show
        plot_feats = []

    # Assign rows to avoid overlapping features
    rows: List[List[Dict]] = []
    for feat in sorted(plot_feats, key=lambda x: x["start"]):
        placed = False
        for row in rows:
            if all(feat["start"] > existing["end"] + 5 for existing in row):
                row.append(feat)
                placed = True
                break
        if not placed:
            rows.append([feat])

    n_rows = max(len(rows), 1)
    fig_h  = max(2.5, 1.2 + n_rows * 0.55)
    fig, ax = plt.subplots(figsize=(12, fig_h), dpi=300)

    # Backbone
    ax.barh(0, seq_length, left=1, height=0.18, color="#CCCCCC", zorder=2)

    # Features
    legend_handles: Dict[str, mpatches.Patch] = {}
    for row_idx, row in enumerate(rows):
        y = -(row_idx * 0.55)
        for feat in row:
            color = _feature_color(feat["type"])
            ax.barh(
                y, feat["length"], left=feat["start"],
                height=0.40, color=color, edgecolor="white", linewidth=0.4,
                zorder=3, alpha=0.92,
            )
            # Label if wide enough
            if feat["length"] > seq_length * 0.04:
                ax.text(
                    feat["start"] + feat["length"] / 2, y,
                    feat["description"][:22],
                    ha="center", va="center", fontsize=6.5, color="white",
                    fontweight="bold", zorder=4, clip_on=True,
                )
            if feat["type"] not in legend_handles:
                legend_handles[feat["type"]] = mpatches.Patch(
                    color=color, label=feat["type"]
                )

    # Axis styling
    ax.set_xlim(0, seq_length + 10)
    ax.set_ylim(-(n_rows * 0.55) - 0.4, 0.6)
    ax.set_xlabel("Amino acid position", fontsize=11)
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title(
        f"{gene} ({accession})  ·  {seq_length} aa",
        fontsize=13, fontweight="bold", pad=10,
    )

    if legend_handles:
        ax.legend(
            handles=list(legend_handles.values()),
            loc="upper right", fontsize=8,
            framealpha=0.9, ncol=min(4, len(legend_handles)),
        )

    plt.tight_layout()
    out_base = os.path.join(outdir, f"{gene}.{accession}.domain_map")
    fig.savefig(f"{out_base}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{out_base}.pdf", bbox_inches="tight")
    plt.close(fig)
    return f"{out_base}.png"


# =========================================================
# Main
# =========================================================

def main():
    parser = argparse.ArgumentParser(
        description="Protein structure and domain analysis for a single gene.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--gene", required=True, help="Gene symbol, e.g. EGFR, TP53, BRCA1")
    parser.add_argument(
        "--organism",
        default="human",
        help="Organism: human, mouse, rat, or NCBI taxon ID",
    )
    parser.add_argument(
        "--modules",
        default="all",
        help=(
            "Comma-separated list of modules to run: "
            "uniprot, pdb, alphafold, domain_map. "
            "Use 'all' to run everything."
        ),
    )
    parser.add_argument(
        "--max-pdb",
        type=int,
        default=20,
        help="Maximum number of PDB structures to retrieve",
    )
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    ensure_dir(args.outdir)

    # Parse modules
    all_modules = {"uniprot", "pdb", "alphafold", "domain_map"}
    if args.modules.strip().lower() == "all":
        modules = all_modules
    else:
        modules = {m.strip().lower() for m in args.modules.split(",")}
        invalid = modules - all_modules
        if invalid:
            raise ValueError(f"Unknown module(s): {invalid}. Choose from: {all_modules}")

    # domain_map requires uniprot
    if "domain_map" in modules:
        modules.add("uniprot")

    summary: Dict = {
        "query_gene": args.gene,
        "organism":   args.organism,
        "modules":    args.modules,
    }

    # ----------------------------------------------------------
    # 1. UniProt resolution + features
    # ----------------------------------------------------------
    accession = None
    seq_length = 0
    features: List[Dict] = []

    if "uniprot" in modules:
        print(f"[INFO] Resolving {args.gene} in UniProt ({args.organism})")
        accession, entry_name, seq_length = resolve_uniprot(args.gene, args.organism)
        print(f"[INFO] UniProt accession: {accession}  ({entry_name},  {seq_length} aa)")

        uniprot_meta = fetch_uniprot_summary(accession)
        summary.update(uniprot_meta)

        print(f"[INFO] Fetching domain/feature annotations from UniProt")
        features = fetch_uniprot_features(accession)
        print(f"[INFO] Features retrieved: {len(features)}")

        feat_df = pd.DataFrame(features)
        feat_path = os.path.join(args.outdir, f"{args.gene}.{accession}.features.tsv")
        feat_df.to_csv(feat_path, sep="\t", index=False)
        print(f"[INFO] Features written: {feat_path}")

        summary["n_features"] = len(features)
        summary["feature_types"] = ";".join(sorted(feat_df["type"].unique())) if not feat_df.empty else ""

    # ----------------------------------------------------------
    # 2. PDB structures
    # ----------------------------------------------------------
    if "pdb" in modules:
        if accession is None:
            print(f"[INFO] Resolving {args.gene} in UniProt ({args.organism}) for PDB search")
            accession, entry_name, seq_length = resolve_uniprot(args.gene, args.organism)

        print(f"[INFO] Fetching PDB structures for {accession} from RCSB")
        pdb_df = fetch_pdb_structures(accession, max_structures=args.max_pdb)
        print(f"[INFO] PDB structures found: {len(pdb_df)}")

        pdb_path = os.path.join(args.outdir, f"{args.gene}.{accession}.pdb_structures.tsv")
        pdb_df.to_csv(pdb_path, sep="\t", index=False)
        print(f"[INFO] PDB table written: {pdb_path}")

        summary["n_pdb_structures"] = len(pdb_df)
        if not pdb_df.empty:
            methods = pdb_df["method"].dropna().unique().tolist()
            summary["pdb_methods"] = ";".join(str(m) for m in methods)
            res_vals = pdb_df["resolution_A"].dropna()
            summary["best_resolution_A"] = float(res_vals.min()) if len(res_vals) else None
            summary["pdb_ids"] = ";".join(pdb_df["pdb_id"].tolist())

    # ----------------------------------------------------------
    # 3. AlphaFold
    # ----------------------------------------------------------
    if "alphafold" in modules:
        if accession is None:
            print(f"[INFO] Resolving {args.gene} in UniProt ({args.organism}) for AlphaFold lookup")
            accession, entry_name, seq_length = resolve_uniprot(args.gene, args.organism)

        print(f"[INFO] Fetching AlphaFold2 entry for {accession}")
        af_entry = fetch_alphafold_entry(accession)

        if af_entry:
            af_df = pd.DataFrame([af_entry])
            af_path = os.path.join(args.outdir, f"{args.gene}.{accession}.alphafold.tsv")
            af_df.to_csv(af_path, sep="\t", index=False)
            print(f"[INFO] AlphaFold entry written: {af_path}")
            summary["alphafold_available"] = True
            summary["alphafold_version"]   = af_entry.get("latest_version")
            summary["alphafold_page"]      = af_entry.get("alphafold_page")
            summary["alphafold_pdb_url"]   = af_entry.get("pdb_url")
            summary["alphafold_cif_url"]   = af_entry.get("cif_url")
        else:
            print(f"[WARN] No AlphaFold entry found for {accession}")
            summary["alphafold_available"] = False

    # ----------------------------------------------------------
    # 4. Domain map
    # ----------------------------------------------------------
    if "domain_map" in modules and features:
        print(f"[INFO] Rendering domain map")
        plot_domain_map(
            features, seq_length,
            gene=args.gene,
            accession=accession,
            outdir=args.outdir,
        )
        print(f"[INFO] Domain map saved to {args.outdir}")
        summary["domain_map_generated"] = True
    elif "domain_map" in modules:
        summary["domain_map_generated"] = False

    # ----------------------------------------------------------
    # Summary outputs
    # ----------------------------------------------------------
    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(
        os.path.join(args.outdir, "protein_structure_summary.tsv"), sep="\t", index=False
    )

    with open(os.path.join(args.outdir, "summary.txt"), "w", encoding="utf-8") as f:
        for k, v in summary.items():
            f.write(f"{k}: {v}\n")

    print(f"[DONE] Results written to: {args.outdir}")


if __name__ == "__main__":
    main()
