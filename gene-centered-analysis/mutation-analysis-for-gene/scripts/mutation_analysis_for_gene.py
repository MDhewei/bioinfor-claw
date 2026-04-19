#!/usr/bin/env python3
"""
Mutation Analysis for Gene: Analyze somatic mutations across TCGA cohorts.

Data source: cBioPortal public API (https://www.cbioportal.org/api/).
TCGA PanCancer Atlas studies are queried by default. A local MAF file
can also be provided via --mutation-file.

Outputs: lollipop plot, mutation-frequency bar chart, mutation-type pie chart,
hotspot details TSV, and a summary TSV.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import urllib.request
import urllib.error
import urllib.parse
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

try:
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), *(['..'] * 3), '_shared'))
    from plot_style import init_style
except ImportError:
    def init_style(**kw): pass


# ─── cBioPortal API ──────────────────────────────────────────────────────────

CBIO_API = "https://www.cbioportal.org/api"

# TCGA PanCancer Atlas study IDs in cBioPortal (suffix _tcga_pan_can_atlas_2018)
TCGA_PANCAN_STUDIES = [
    "acc_tcga_pan_can_atlas_2018", "blca_tcga_pan_can_atlas_2018",
    "brca_tcga_pan_can_atlas_2018", "cesc_tcga_pan_can_atlas_2018",
    "chol_tcga_pan_can_atlas_2018", "coadread_tcga_pan_can_atlas_2018",
    "dlbc_tcga_pan_can_atlas_2018", "esca_tcga_pan_can_atlas_2018",
    "gbm_tcga_pan_can_atlas_2018", "hnsc_tcga_pan_can_atlas_2018",
    "kich_tcga_pan_can_atlas_2018", "kirc_tcga_pan_can_atlas_2018",
    "kirp_tcga_pan_can_atlas_2018", "laml_tcga_pan_can_atlas_2018",
    "lgg_tcga_pan_can_atlas_2018", "lihc_tcga_pan_can_atlas_2018",
    "luad_tcga_pan_can_atlas_2018", "lusc_tcga_pan_can_atlas_2018",
    "meso_tcga_pan_can_atlas_2018", "ov_tcga_pan_can_atlas_2018",
    "paad_tcga_pan_can_atlas_2018", "pcpg_tcga_pan_can_atlas_2018",
    "prad_tcga_pan_can_atlas_2018", "sarc_tcga_pan_can_atlas_2018",
    "skcm_tcga_pan_can_atlas_2018", "stad_tcga_pan_can_atlas_2018",
    "tgct_tcga_pan_can_atlas_2018", "thca_tcga_pan_can_atlas_2018",
    "thym_tcga_pan_can_atlas_2018", "ucec_tcga_pan_can_atlas_2018",
    "ucs_tcga_pan_can_atlas_2018", "uvm_tcga_pan_can_atlas_2018",
]

# Map short TCGA codes (e.g. BRCA) → cBioPortal study IDs
_CODE_TO_STUDY = {}
for _s in TCGA_PANCAN_STUDIES:
    _code = _s.replace("_tcga_pan_can_atlas_2018", "").upper()
    _CODE_TO_STUDY[_code] = _s


def _cbio_get(endpoint: str, params: dict | None = None, timeout: int = 30):
    """GET from cBioPortal REST API. Returns parsed JSON."""
    url = f"{CBIO_API}{endpoint}"
    if params:
        url += "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode())


def _cbio_post(endpoint: str, body: dict, timeout: int = 60,
               params: dict | None = None):
    """POST JSON to cBioPortal REST API. Returns parsed JSON."""
    url = f"{CBIO_API}{endpoint}"
    if params:
        url += "?" + urllib.parse.urlencode(params)
    data = json.dumps(body).encode()
    req = urllib.request.Request(url, data=data, method="POST",
                                headers={"Content-Type": "application/json",
                                         "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode())


def _resolve_entrez_id(gene_symbol: str) -> Optional[int]:
    """Look up Entrez Gene ID from cBioPortal."""
    try:
        data = _cbio_get(f"/genes/{gene_symbol}")
        return data.get("entrezGeneId")
    except Exception:
        return None


def fetch_cbio_mutations(gene_symbol: str,
                         cancer_types: Optional[List[str]] = None,
                         ) -> Tuple[List[Dict], Dict[str, int], Dict[str, int]]:
    """
    Fetch mutations from cBioPortal for a gene across TCGA PanCancer studies.

    Returns (mutations, mut_counts_per_study, total_profiled_per_study).
    """
    entrez_id = _resolve_entrez_id(gene_symbol)
    if entrez_id is None:
        print(f"[WARN] Could not resolve Entrez ID for {gene_symbol}; using symbol only")

    # Determine which studies to query
    if cancer_types and cancer_types != ["all"]:
        studies = []
        for ct in cancer_types:
            ct_upper = ct.upper().replace("TCGA-", "")
            if ct_upper in _CODE_TO_STUDY:
                studies.append(_CODE_TO_STUDY[ct_upper])
            else:
                # Try literal
                studies.append(ct)
    else:
        studies = list(TCGA_PANCAN_STUDIES)

    # Collect molecular-profile IDs for mutations
    profile_ids = []
    study_for_profile = {}
    for sid in studies:
        try:
            profiles = _cbio_get(f"/studies/{sid}/molecular-profiles")
            for p in profiles:
                if p.get("molecularAlterationType") == "MUTATION_EXTENDED":
                    pid = p["molecularProfileId"]
                    profile_ids.append(pid)
                    study_for_profile[pid] = sid
        except Exception as e:
            print(f"  [skip] {sid}: {e}")

    if not profile_ids:
        print("[ERROR] No mutation profiles found", file=sys.stderr)
        return [], {}, {}

    print(f"[INFO] Querying {len(profile_ids)} mutation profiles for {gene_symbol} "
          f"(projection=DETAILED)…")

    # Fetch mutations via POST /mutations/fetch
    all_mutations: List[Dict] = []
    mut_counts: Dict[str, int] = defaultdict(int)
    profiled_counts: Dict[str, int] = {}

    # Process in batches of profiles to avoid huge single requests
    batch_size = 10
    for i in range(0, len(profile_ids), batch_size):
        batch = profile_ids[i:i + batch_size]
        body = {
            "molecularProfileIds": batch,
            "entrezGeneIds": [entrez_id] if entrez_id else [],
        }
        if not entrez_id:
            # Fall back: query each profile individually by gene keyword
            for pid in batch:
                try:
                    muts = _cbio_get(f"/molecular-profiles/{pid}/mutations",
                                     params={"entrezGeneId": 0, "sampleListId": pid.replace("_mutations", "_all")})
                    # Filter by gene symbol
                    muts = [m for m in muts if m.get("gene", {}).get("hugoGeneSymbol", "").upper() == gene_symbol.upper()]
                    all_mutations.extend(muts)
                except Exception:
                    pass
            continue

        try:
            muts = _cbio_post("/mutations/fetch", body,
                              params={"projection": "DETAILED"})
            all_mutations.extend(muts)
        except Exception as e:
            print(f"  [batch error] profiles {i}-{i+len(batch)}: {e}")

    # Count profiled samples per study
    for pid in profile_ids:
        sid = study_for_profile[pid]
        try:
            sample_list_id = sid + "_all"
            samples = _cbio_get(f"/sample-lists/{sample_list_id}")
            profiled_counts[sid] = len(samples.get("sampleIds", []))
        except Exception:
            profiled_counts[sid] = 0

    # Diagnostic: check a sample mutation for available fields
    if all_mutations:
        sample = all_mutations[0]
        prot_fields = {k: v for k, v in sample.items()
                       if 'prot' in k.lower() or 'amino' in k.lower() or 'change' in k.lower()}
        print(f"[DEBUG] Sample mutation fields with 'prot/amino/change': {prot_fields}")
        n_with_prot = sum(1 for m in all_mutations
                         if m.get("proteinChange") or m.get("aminoAcidChange"))
        print(f"[INFO] {n_with_prot}/{len(all_mutations)} mutations have protein change info")

    # Parse mutations into our standard format
    mutations = []
    for m in all_mutations:
        pid = m.get("molecularProfileId", "")
        sid = study_for_profile.get(pid, "Unknown")
        cancer_code = sid.replace("_tcga_pan_can_atlas_2018", "").upper()

        # Protein change: try multiple field names (varies by API projection)
        prot_change = m.get("proteinChange") or m.get("aminoAcidChange") or ""
        if not prot_change or prot_change in ("NA", "N/A", "MUTATED"):
            prot_change = ""
        # Ensure p. prefix for consistency
        if prot_change and not prot_change.startswith("p."):
            prot_change = f"p.{prot_change}"

        mut_type = m.get("mutationType", "Other") or "Other"
        # Keyword field is sometimes different
        if mut_type == "Other":
            mut_type = m.get("mutationStatus", "Other")

        aa_changes = [prot_change] if prot_change else []

        mutations.append({
            "ssm_id": m.get("uniqueSampleKey", ""),
            "project_id": f"TCGA-{cancer_code}",
            "consequence_types": [mut_type],
            "aa_changes": aa_changes,
            "genomic_change": f"chr{m.get('chr', '')}:g.{m.get('startPosition', '')}"
        })
        mut_counts[f"TCGA-{cancer_code}"] += 1

    # Convert profiled_counts keys to TCGA-XX format
    profiled = {}
    for sid, n in profiled_counts.items():
        code = sid.replace("_tcga_pan_can_atlas_2018", "").upper()
        profiled[f"TCGA-{code}"] = n

    return mutations, dict(mut_counts), profiled


# ─── Local MAF loader ────────────────────────────────────────────────────────

def load_maf_mutations(filepath: str, gene_symbol: str) -> Tuple[List[Dict], Dict[str, int], Dict[str, int]]:
    """Load mutations from a local MAF file."""
    print(f"[INFO] Loading MAF mutations from: {filepath}")
    with open(filepath) as f:
        skip = 0
        for line in f:
            if line.startswith('#'):
                skip += 1
            else:
                break
    df = pd.read_csv(filepath, sep='\t', skiprows=skip, low_memory=False)

    gene_col = next((c for c in df.columns if c in ['Hugo_Symbol', 'HugoSymbol', 'gene']), None)
    if gene_col is None:
        print("[ERROR] Cannot find gene column in MAF", file=sys.stderr)
        return [], {}, {}

    gene_df = df[df[gene_col].str.upper() == gene_symbol.upper()].copy()
    print(f"[INFO] Found {len(gene_df)} {gene_symbol} mutations")

    mutations = []
    mut_counts: Dict[str, int] = defaultdict(int)
    prot_col = next((c for c in df.columns if c in ['HGVSp_Short', 'Protein_Change', 'AAChange']), None)
    cons_col = next((c for c in df.columns if c in ['Variant_Classification', 'Consequence']), None)
    project_col = next((c for c in df.columns if c in ['project_id', 'Tumor_Sample_Barcode']), None)

    for _, row in gene_df.iterrows():
        prot_change = str(row.get(prot_col, '')) if prot_col else ''
        cons_type = str(row.get(cons_col, '')) if cons_col else 'Other'
        project = str(row.get(project_col, 'Unknown')) if project_col else 'Unknown'
        if project_col and 'Barcode' in str(project_col):
            parts = project.split('-')
            if len(parts) >= 3 and parts[0] == 'TCGA':
                project = f"TCGA-{parts[1]}"

        aa_changes = [prot_change] if prot_change and prot_change not in ('nan', 'p.?', '') else []
        mutations.append({
            "ssm_id": f"{row.get('Chromosome', '')}:{row.get('Start_Position', '')}",
            "project_id": project,
            "consequence_types": [cons_type],
            "aa_changes": aa_changes,
            "genomic_change": ""
        })
        mut_counts[project] += 1

    return mutations, dict(mut_counts), {}


# ─── Mutation classification ─────────────────────────────────────────────────

_MUTTYPE_PRIORITY = {
    "Nonsense_Mutation": 1, "Frameshift": 2, "Frame_Shift_Del": 2,
    "Frame_Shift_Ins": 2, "Splice_Site": 3, "Splice_Region": 3,
    "Missense_Mutation": 4, "Missense": 4,
    "In_Frame_Del": 5, "In_Frame_Ins": 5, "In_Frame": 5,
    "Nonstop_Mutation": 6, "Translation_Start_Site": 7,
    "Silent": 8,
}

_MUTTYPE_NORMALIZE = {
    "missense": "Missense_Mutation", "missense_mutation": "Missense_Mutation",
    "missense_variant": "Missense_Mutation",
    "nonsense_mutation": "Nonsense_Mutation", "nonsense": "Nonsense_Mutation",
    "stop_gained": "Nonsense_Mutation",
    "frame_shift_del": "Frameshift", "frame_shift_ins": "Frameshift",
    "frameshift_variant": "Frameshift", "frameshift": "Frameshift",
    "frameshift_deletion": "Frameshift", "frameshift_insertion": "Frameshift",
    "splice_site": "Splice_Site", "splice_donor_variant": "Splice_Site",
    "splice_acceptor_variant": "Splice_Site", "splice_region": "Splice_Site",
    "in_frame_del": "In_Frame", "in_frame_ins": "In_Frame",
    "inframe_deletion": "In_Frame", "inframe_insertion": "In_Frame",
    "silent": "Silent", "synonymous_variant": "Silent",
}


def classify_mutation_type(consequence_types: List[str]) -> str:
    """Classify mutation from a list of consequence type strings."""
    best, best_pri = "Other", 99
    for ct in consequence_types:
        norm = _MUTTYPE_NORMALIZE.get(ct.lower().strip(), ct)
        pri = _MUTTYPE_PRIORITY.get(norm, 90)
        if pri < best_pri:
            best, best_pri = norm, pri
    return best


# ─── Amino-acid position extraction ──────────────────────────────────────────

def extract_amino_acid_position(aa_change: str) -> Optional[Tuple[str, int]]:
    """Parse p.V600E-style strings → ('V600E', 600)."""
    if not aa_change:
        return None
    m = re.search(r'p\.([A-Z*])(\d+)([A-Z*_].*)?', aa_change, re.IGNORECASE)
    if m:
        try:
            pos = int(m.group(2))
            ref = m.group(1).upper()
            alt = (m.group(3) or '?').lstrip('_').upper()
            if len(alt) > 5:
                alt = alt[:5]
            return (f"{ref}{pos}{alt}", pos)
        except (ValueError, IndexError):
            pass
    return None


# ─── Statistics ──────────────────────────────────────────────────────────────

def compute_statistics(
    mutations: List[Dict],
    gene_symbol: str,
    mutation_type_filter: List[str],
    profiled_counts: Dict[str, int],
) -> Tuple[pd.DataFrame, Dict, Dict, int]:
    """Compute mutation stats: hotspots, frequency, type distribution."""
    hotspots: Dict[int, List] = defaultdict(list)
    type_counts: Dict[str, int] = defaultdict(int)
    cancer_mutation_counts: Dict[str, int] = defaultdict(int)

    for mutation in mutations:
        mut_type = classify_mutation_type(mutation["consequence_types"])
        if mutation_type_filter and mut_type not in mutation_type_filter:
            continue
        type_counts[mut_type] += 1
        project = mutation["project_id"]
        cancer_mutation_counts[project] += 1

        for aa_change in mutation["aa_changes"]:
            result = extract_amino_acid_position(aa_change)
            if result:
                aa_str, position = result
                hotspots[position].append({
                    "aa_change": aa_str, "type": mut_type, "project": project
                })

    # Build summary
    rows = []
    for project in sorted(cancer_mutation_counts):
        n_muts = cancer_mutation_counts[project]
        n_profiled = profiled_counts.get(project, 0)
        freq = n_muts / n_profiled if n_profiled > 0 else 0.0

        project_hotspots = sorted(
            [p for p, ms in hotspots.items() if any(m["project"] == project for m in ms)],
            key=lambda p: len(hotspots[p]), reverse=True
        )[:5]

        rows.append({
            "cancer_type": project.replace("TCGA-", ""),
            "n_mutations": n_muts,
            "n_profiled": n_profiled,
            "frequency": round(freq, 4),
            "top_hotspots": ";".join(str(p) for p in project_hotspots),
        })

    return pd.DataFrame(rows), dict(hotspots), dict(type_counts), len(mutations)


# ─── Protein length ──────────────────────────────────────────────────────────

def fetch_uniprot_protein_length(gene_symbol: str) -> Optional[int]:
    """Fetch protein length from UniProt."""
    url = (f"https://rest.uniprot.org/uniprotkb/search?"
           f"query=gene_exact:{gene_symbol}+AND+organism_id:9606+AND+reviewed:true"
           f"&fields=length&format=json")
    try:
        with urllib.request.urlopen(url, timeout=10) as resp:
            data = json.loads(resp.read().decode())
            results = data.get("results", [])
            if results:
                return int(results[0].get("sequence", {}).get("length", 0)) or None
    except Exception:
        pass
    return None


# Well-known protein lengths (fallback when UniProt is unreachable)
_KNOWN_LENGTHS = {
    "TP53": 393, "KRAS": 189, "BRAF": 766, "EGFR": 1210,
    "BRCA1": 1863, "BRCA2": 3418, "PIK3CA": 1068, "APC": 2843,
    "PTEN": 403, "RB1": 928, "MYC": 439, "NRAS": 189,
    "HRAS": 189, "IDH1": 414, "IDH2": 452, "VHL": 213,
    "CTNNB1": 781, "SMAD4": 552, "CDKN2A": 156, "ARID1A": 2285,
    "KMT2D": 5537, "ATM": 3056, "NOTCH1": 2555, "FBXW7": 707,
    "NFE2L2": 605, "KEAP1": 624, "STK11": 433, "FGFR3": 806,
    "ERBB2": 1255, "ALK": 1620, "ROS1": 2347, "MET": 1390,
    "KIT": 976, "PDGFRA": 1089, "RET": 1114, "FGFR2": 821,
    "NF1": 2818, "NF2": 595, "PRNP": 253,
}


def get_protein_length(gene: str, user_length: Optional[int]) -> int:
    if user_length:
        return user_length
    length = fetch_uniprot_protein_length(gene)
    if length:
        return length
    if gene.upper() in _KNOWN_LENGTHS:
        print(f"[INFO] Using built-in protein length for {gene}")
        return _KNOWN_LENGTHS[gene.upper()]
    print(f"[WARN] Unknown protein length for {gene}; using 500", file=sys.stderr)
    return 500


# ─── Plotting ────────────────────────────────────────────────────────────────

TYPE_COLORS = {
    "Missense_Mutation": "#1f77b4",
    "Nonsense_Mutation": "#d62728",
    "Frameshift": "#2ca02c",
    "Splice_Site": "#ff7f0e",
    "In_Frame": "#9467bd",
    "Silent": "#bcbd22",
    "Other": "#7f7f7f",
}


def plot_lollipop(hotspots, protein_length, domains, gene, output_file):
    """Lollipop plot of protein mutations."""
    fig, (ax_main, ax_dom) = plt.subplots(
        2, 1, figsize=(14, 8), gridspec_kw={'height_ratios': [4, 1]})

    sorted_hs = sorted(hotspots.items(), key=lambda x: len(x[1]), reverse=True)
    max_count = max((len(v) for v in hotspots.values()), default=1)

    for pos, muts in sorted_hs:
        if pos < 1 or pos > protein_length:
            continue
        count = len(muts)
        primary = max(set(m["type"] for m in muts),
                      key=lambda t: sum(1 for m in muts if m["type"] == t))
        color = TYPE_COLORS.get(primary, "#7f7f7f")
        ax_main.plot([pos, pos], [0, count], color=color, lw=2, alpha=0.5, zorder=2)
        ax_main.scatter(pos, count, s=200, color=color, alpha=0.7, zorder=3,
                        edgecolors="black", linewidth=0.5)
        if count >= max(3, max_count * 0.15):
            label = muts[0]["aa_change"] if muts[0].get("aa_change") else str(pos)
            ax_main.text(pos, count + max_count * 0.03, label,
                         ha="center", fontsize=7, fontweight="bold", rotation=45)

    ax_main.set_xlim(0, protein_length + 10)
    ax_main.set_ylim(0, max_count * 1.15)
    ax_main.set_xlabel("Amino Acid Position")
    ax_main.set_ylabel("Mutation Count")
    ax_main.set_title(f"{gene} — Protein Mutation Lollipop Plot (TCGA)", fontweight="bold")
    ax_main.grid(axis="y", alpha=0.3, linestyle="--")
    legend_els = [mpatches.Patch(color=c, label=t)
                  for t, c in TYPE_COLORS.items() if t != "Other"]
    ax_main.legend(handles=legend_els, loc="upper right", fontsize=9)

    # Domain track
    ax_dom.set_xlim(0, protein_length + 10)
    ax_dom.set_ylim(0, 1)
    # Draw protein backbone
    ax_dom.barh(0.5, protein_length, left=0, height=0.15, color="#dddddd",
                edgecolor="black", linewidth=0.5)
    for d in domains:
        ax_dom.barh(0.5, d["end"] - d["start"], left=d["start"], height=0.3,
                     alpha=0.6, edgecolor="black", linewidth=1)
        ax_dom.text((d["start"] + d["end"]) / 2, 0.5, d["name"],
                     ha="center", va="center", fontsize=8)
    ax_dom.set_yticks([])
    ax_dom.set_xlabel("Amino Acid Position")
    ax_dom.set_title("Protein Domains", fontsize=11)
    for spine in ("left", "right", "top"):
        ax_dom.spines[spine].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_file}")


def plot_mutation_frequency(summary_df, gene, output_file):
    """Bar chart of mutation frequency by cancer type."""
    if summary_df.empty:
        print("[WARN] No data for frequency plot")
        return
    df = summary_df[summary_df["frequency"] > 0].sort_values("frequency", ascending=True)
    if df.empty:
        df = summary_df.sort_values("n_mutations", ascending=True)

    fig, ax = plt.subplots(figsize=(10, max(4, len(df) * 0.3)))
    colors = plt.cm.RdYlBu_r(np.linspace(0.2, 0.8, len(df)))
    ax.barh(df["cancer_type"], df["frequency"], color=colors,
            edgecolor="black", linewidth=0.5)
    ax.set_xlabel("Mutation Frequency (fraction of profiled cases)")
    ax.set_title(f"{gene} — Mutation Frequency by Cancer Type (TCGA)", fontweight="bold")
    ax.grid(axis="x", alpha=0.3, linestyle="--")
    for i, (ct, freq) in enumerate(zip(df["cancer_type"], df["frequency"])):
        ax.text(freq, i, f" {freq:.3f}", va="center", fontsize=8)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_file}")


def plot_mutation_types(type_counts, gene, output_file):
    """Pie chart of mutation type distribution."""
    if not type_counts:
        print("[WARN] No type data for pie chart")
        return
    # Filter out tiny categories
    total = sum(type_counts.values())
    filtered = {k: v for k, v in type_counts.items() if v / total > 0.01}
    other = total - sum(filtered.values())
    if other > 0:
        filtered["Other"] = filtered.get("Other", 0) + other

    types = list(filtered.keys())
    counts = list(filtered.values())
    colors = [TYPE_COLORS.get(t, "#7f7f7f") for t in types]

    fig, ax = plt.subplots(figsize=(8, 8))
    wedges, texts, autotexts = ax.pie(
        counts, labels=types, autopct="%1.1f%%", colors=colors,
        startangle=90, textprops={"fontsize": 10})
    for at in autotexts:
        at.set_color("white")
        at.set_fontweight("bold")
    ax.set_title(f"{gene} — Mutation Type Distribution (TCGA)", fontweight="bold")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_file}")


# ─── Domain loading ──────────────────────────────────────────────────────────

def load_domain_file(path: str) -> List[Dict]:
    try:
        df = pd.read_csv(path, sep="\t")
        return [{"name": r.get("domain_name", r.get("name", "")),
                 "start": int(r.get("start_aa", r.get("start", 0))),
                 "end": int(r.get("end_aa", r.get("end", 0)))}
                for _, r in df.iterrows()]
    except Exception as e:
        print(f"[WARN] Failed to load domain file: {e}", file=sys.stderr)
        return []


# ─── Main ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Analyze somatic mutations in a gene across TCGA cohorts (via cBioPortal).")
    parser.add_argument("--gene", required=True, help="Gene symbol (e.g., TP53, KRAS)")
    parser.add_argument("--mutation-file", default=None,
                        help="Local MAF file (optional; uses cBioPortal API if omitted)")
    parser.add_argument("--cancer-types", default="all",
                        help="Comma-separated TCGA codes (e.g., BRCA,LUAD) or 'all'")
    parser.add_argument("--mutation-types", default="",
                        help="Filter: comma-separated types (Missense_Mutation, Nonsense_Mutation, etc.)")
    parser.add_argument("--min-frequency", type=float, default=0.01,
                        help="Min frequency to label on plots (default: 0.01)")
    parser.add_argument("--top-hotspots", type=int, default=10,
                        help="N top hotspots to highlight (default: 10)")
    parser.add_argument("--protein-length", type=int, default=None,
                        help="Protein length in AA (auto-fetched if omitted)")
    parser.add_argument("--domain-file", default=None, help="TSV: domain_name, start_aa, end_aa")
    parser.add_argument("--outdir", default=".", help="Output directory")
    parser.add_argument("--font-family", default=None)
    parser.add_argument("--font-size", type=float, default=None)
    args = parser.parse_args()

    init_style(font_family=args.font_family, font_size=args.font_size)
    os.makedirs(args.outdir, exist_ok=True)
    gene = args.gene.upper()

    # Parse cancer types
    if args.cancer_types.lower() == "all":
        cancer_types = None
    else:
        cancer_types = [ct.strip() for ct in args.cancer_types.split(",")]

    mutation_type_filter = [t.strip() for t in args.mutation_types.split(",") if t.strip()]

    # ── Load mutations ────────────────────────────────────────────────────
    profiled_counts: Dict[str, int] = {}

    if args.mutation_file and os.path.exists(args.mutation_file):
        mutations, mut_counts, profiled_counts = load_maf_mutations(args.mutation_file, gene)
    else:
        print(f"[INFO] Fetching {gene} mutations from cBioPortal (TCGA PanCancer Atlas)…")
        try:
            mutations, mut_counts, profiled_counts = fetch_cbio_mutations(gene, cancer_types)
        except Exception as e:
            print(f"[ERROR] cBioPortal API failed: {e}", file=sys.stderr)
            print("  Hint: supply --mutation-file with a local TCGA MAF file.", file=sys.stderr)
            return 1

    if not mutations:
        print(f"[ERROR] No mutations found for {gene}.", file=sys.stderr)
        return 1

    print(f"[INFO] Found {len(mutations)} mutations across {len(mut_counts)} cancer types")

    # ── Statistics ────────────────────────────────────────────────────────
    summary_df, hotspots, type_counts, total_muts = compute_statistics(
        mutations, gene, mutation_type_filter, profiled_counts)

    # ── Protein length ────────────────────────────────────────────────────
    protein_length = get_protein_length(gene, args.protein_length)
    print(f"[INFO] Protein length: {protein_length} AA")

    # ── Domains ───────────────────────────────────────────────────────────
    domains = load_domain_file(args.domain_file) if args.domain_file else []

    # ── Outputs ───────────────────────────────────────────────────────────
    summary_file = os.path.join(args.outdir, "mutation_summary.tsv")
    summary_df.to_csv(summary_file, sep="\t", index=False)
    print(f"Saved: {summary_file}")

    plot_lollipop(hotspots, protein_length, domains, gene,
                  os.path.join(args.outdir, "lollipop_plot.png"))
    plot_mutation_frequency(summary_df, gene,
                           os.path.join(args.outdir, "mutation_frequency.png"))
    plot_mutation_types(type_counts, gene,
                       os.path.join(args.outdir, "mutation_types.png"))

    # Hotspot details
    hs_rows = []
    for pos in sorted(hotspots):
        muts = hotspots[pos]
        aa_labels = set(m["aa_change"] for m in muts if m.get("aa_change"))
        types = set(m["type"] for m in muts)
        hs_rows.append({
            "position": pos, "count": len(muts),
            "aa_changes": ";".join(sorted(aa_labels)),
            "types": ";".join(sorted(types)),
        })
    hs_df = pd.DataFrame(hs_rows).sort_values("count", ascending=False)
    hs_file = os.path.join(args.outdir, "hotspot_details.tsv")
    hs_df.to_csv(hs_file, sep="\t", index=False)
    print(f"Saved: {hs_file}")

    # ── Summary ───────────────────────────────────────────────────────────
    print(f"\n=== Mutation Analysis Summary for {gene} ===")
    print(f"Total mutations        : {total_muts}")
    print(f"Cancer types           : {len(mut_counts)}")
    print(f"Mutation types         : {dict(type_counts)}")
    if len(hs_df) > 0:
        top = hs_df.iloc[0]
        print(f"Top hotspot            : position {top['position']} ({top['count']} mutations)")
    if not summary_df.empty:
        top_ct = summary_df.sort_values("frequency", ascending=False).iloc[0]
        print(f"Highest frequency      : {top_ct['cancer_type']} ({top_ct['frequency']:.3f})")
    print(f"Output dir             : {args.outdir}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
