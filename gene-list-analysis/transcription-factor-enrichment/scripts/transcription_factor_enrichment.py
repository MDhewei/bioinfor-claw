#!/usr/bin/env python3
"""
Transcription Factor Enrichment Analysis
Identifies TFs that significantly regulate a given gene list using
ChEA3 API and DoRothEA built-in database.
"""

import argparse
import sys
import os
import json
from collections import defaultdict

import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Rectangle


# Built-in DoRothEA: compact human/mouse TF-target regulons (high confidence A+B)
DOROTHEA_HUMAN = {
    'TP53': {'targets': ['CDKN1A', 'MDM2', 'PMAIP1', 'BBC3', 'GADD45A', 'SESN2', 'ZMAT3', 'FAS', 'DDB2'], 'confidence': 'A'},
    'MYC': {'targets': ['LDHA', 'LDHB', 'PFKFB4', 'GLUT1', 'ENO1', 'FASTK', 'ODC1', 'CYC1', 'UQCRC2'], 'confidence': 'A'},
    'E2F1': {'targets': ['CCNE1', 'CCNB1', 'CCNA2', 'DHFR', 'PCNA', 'RRM1', 'RRM2', 'CDK2', 'CCNE2'], 'confidence': 'A'},
    'STAT3': {'targets': ['SOCS3', 'IL6', 'IL6R', 'ICAM1', 'BCAM', 'CREBBP', 'PIAS3', 'CISH', 'PIM1'], 'confidence': 'A'},
    'RELA': {'targets': ['NFKBIA', 'IL8', 'TNF', 'IL6', 'ICAM1', 'VCAM1', 'RELB', 'NFKB2', 'PTGS2'], 'confidence': 'A'},
    'SP1': {'targets': ['CDKN1A', 'TIMP3', 'EGFR', 'FOS', 'JUN', 'TFF1', 'EREG', 'HBEGF', 'AREG'], 'confidence': 'A'},
    'HIF1A': {'targets': ['VEGFA', 'PGK1', 'LDHA', 'SLC2A1', 'CA9', 'ENO1', 'NDRG1', 'ADM', 'LOXL2'], 'confidence': 'A'},
    'FOXO3': {'targets': ['CDKN1B', 'CDKN2D', 'FAS', 'PMAIP1', 'SOD2', 'GADD45A', 'BIM', 'CAT', 'MAPK1'], 'confidence': 'A'},
    'TP63': {'targets': ['CDKN1A', 'PERP', 'PUMA', 'BAX', 'SESN3', 'SPRR2A', 'SPRR2E', 'KRT5', 'IVL'], 'confidence': 'A'},
    'TP73': {'targets': ['CDKN1A', 'PUMA', 'NOXA', 'DDB2', 'XPC', 'GADD45A', 'SESN1', 'BAX', 'FAS'], 'confidence': 'A'},
    'NFKB1': {'targets': ['NFKBIA', 'IL8', 'TNF', 'IL6', 'ICAM1', 'VCAM1', 'RELB', 'PTGS2', 'CCL2'], 'confidence': 'A'},
    'JUN': {'targets': ['FOS', 'FOSB', 'JunB', 'CDKN1A', 'MMP9', 'MMP2', 'IL2', 'IL8', 'BIRC3'], 'confidence': 'A'},
    'FOS': {'targets': ['JUN', 'JUNB', 'CDKN1A', 'MMP9', 'TNF', 'IL2', 'IL8', 'FOSB', 'MMP1'], 'confidence': 'A'},
    'GATA3': {'targets': ['IL5', 'IL13', 'IL4', 'STAT5A', 'STAT5B', 'IFNG', 'ERBB2', 'FOXP3', 'RUNX3'], 'confidence': 'A'},
    'GATA1': {'targets': ['HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'HBG1', 'HBZ', 'ALAS2', 'KLF1'], 'confidence': 'A'},
    'TBX21': {'targets': ['IFNG', 'IL2', 'CCL1', 'CCL3', 'CCL4', 'CCL5', 'STAT4', 'IL12RB2', 'CXCL10'], 'confidence': 'A'},
    'ESR1': {'targets': ['PGR', 'PRS', 'GREB1', 'SDF1', 'FOXO3', 'p21', 'AREG', 'EREG', 'BGN'], 'confidence': 'A'},
    'AR': {'targets': ['PSA', 'TMPRSS2', 'KLK2', 'KLK3', 'NKX3', 'PCA3', 'SPOP', 'FOXA1', 'CTBP2'], 'confidence': 'A'},
    'TCF7': {'targets': ['CCND1', 'c-Myc', 'LEF1', 'AXIN2', 'DKK1', 'NODAL', 'NANOG', 'SOX2', 'OCT4'], 'confidence': 'A'},
    'CTNNB1': {'targets': ['CCND1', 'c-Myc', 'AXIN2', 'DKK1', 'NODAL', 'SOX2', 'NANOG', 'TCF4', 'LEF1'], 'confidence': 'A'},
    'PTEN': {'targets': ['PTEN', 'PHIP', 'PCTA', 'ATP7A', 'ATP7B', 'LDLR', 'SLC9A1', 'EMP1', 'GJA1'], 'confidence': 'B'},
    'AHR': {'targets': ['CYP1A1', 'CYP1A2', 'CYP1B1', 'TIPARP', 'IDO1', 'ALDH3A1', 'NFE2L2', 'IL22', 'IL17'], 'confidence': 'B'},
    'NR3C1': {'targets': ['FKBP5', 'GILZ', 'TSC22D3', 'DUSP1', 'CDKN1A', 'IL1R1', 'SGK1', 'IGFBP5', 'KLF9'], 'confidence': 'A'},
    'CEBPA': {'targets': ['C/EBPB', 'C/EBPD', 'CD14', 'CD33', 'G-CSF', 'M-CSF', 'IL6', 'TNF', 'PTGS2'], 'confidence': 'A'},
    'CEBPB': {'targets': ['CEBPA', 'CEBPD', 'IL6', 'TNF', 'IL1B', 'IL8', 'G-CSF', 'LIF', 'SOCS3'], 'confidence': 'A'},
    'RUNX1': {'targets': ['RUNX2', 'RUNX3', 'CD19', 'RAG1', 'RAG2', 'NOTCH1', 'IL2', 'IL3', 'IL7'], 'confidence': 'A'},
    'RUNX2': {'targets': ['ALP', 'BGLAP', 'IBSP', 'SPP1', 'COL10A1', 'RANKL', 'RUNX1', 'RUNX3', 'OSTERIX'], 'confidence': 'A'},
    'KLF4': {'targets': ['CDKN1A', 'CDKN1B', 'CDKN2B', 'OCT4', 'SOX2', 'NANOG', 'c-Myc', 'POSTN', 'TGFBI'], 'confidence': 'B'},
    'SOX2': {'targets': ['OCT4', 'NANOG', 'FGF4', 'HESX1', 'LEFTY1', 'LEFTY2', 'NODAL', 'GDF3', 'UTF1'], 'confidence': 'A'},
    'FOXA1': {'targets': ['ESR1', 'PGR', 'KLK2', 'KLK3', 'GATA3', 'AR', 'PSA', 'TMPRSS2', 'FOS'], 'confidence': 'A'},
    'FOXA2': {'targets': ['PDX1', 'HNF4A', 'HNF1A', 'NEUROD1', 'NEUROG3', 'INSULIN', 'GCG', 'TTR', 'APOB'], 'confidence': 'A'},
}

DOROTHEA_MOUSE = {
    'Tp53': {'targets': ['Cdkn1a', 'Mdm2', 'Pmaip1', 'Bbc3', 'Gadd45a', 'Sesn2', 'Ddb2', 'Fas', 'Zmat3'], 'confidence': 'A'},
    'Myc': {'targets': ['Ldha', 'Ldhb', 'Pfkfb4', 'Slc2a1', 'Eno1', 'Odc1', 'Cyc1', 'Pcna', 'Top2a'], 'confidence': 'A'},
    'E2f1': {'targets': ['Ccne1', 'Ccnb1', 'Ccna2', 'Dhfr', 'Pcna', 'Rrm1', 'Rrm2', 'Cdk2', 'Ccne2'], 'confidence': 'A'},
    'Stat3': {'targets': ['Socs3', 'Il6', 'Icam1', 'Bcam', 'Crebbp', 'Pias3', 'Cish', 'Pim1', 'Ctf22d3'], 'confidence': 'A'},
    'Nfkb1': {'targets': ['Nfkbia', 'Il8', 'Tnf', 'Il6', 'Icam1', 'Vcam1', 'Relb', 'Nfkb2', 'Ptgs2'], 'confidence': 'A'},
    'Rel': {'targets': ['Nfkbia', 'Il8', 'Tnf', 'Il6', 'Icam1', 'Vcam1', 'Cxcl10', 'Il12a', 'Ccl2'], 'confidence': 'A'},
    'Sp1': {'targets': ['Cdkn1a', 'Timp3', 'Egfr', 'Fos', 'Jun', 'Tff1', 'Ereg', 'Areg', 'Hbegf'], 'confidence': 'A'},
    'Hif1a': {'targets': ['Vegfa', 'Pgk1', 'Ldha', 'Slc2a1', 'Ca9', 'Eno1', 'Ndrg1', 'Adm', 'Loxl2'], 'confidence': 'A'},
    'Foxo3': {'targets': ['Cdkn1b', 'Cdkn2d', 'Fas', 'Pmaip1', 'Sod2', 'Gadd45a', 'Bim', 'Cat', 'Mapk1'], 'confidence': 'A'},
    'Jun': {'targets': ['Fos', 'Fosb', 'Junb', 'Cdkn1a', 'Mmp9', 'Mmp2', 'Il2', 'Il8', 'Birc3'], 'confidence': 'A'},
    'Fos': {'targets': ['Jun', 'Junb', 'Fosl2', 'Cdkn1a', 'Mmp9', 'Tnf', 'Il2', 'Il8', 'Mmp1'], 'confidence': 'A'},
    'Gata3': {'targets': ['Il5', 'Il13', 'Il4', 'Stat5a', 'Stat5b', 'Ifng', 'Rbpj', 'Foxp3', 'Runx3'], 'confidence': 'A'},
    'Gata1': {'targets': ['Hba1', 'Hba2', 'Hbb', 'Hbd', 'Hbe1', 'Hbg1', 'Hbz', 'Alas2', 'Klf1'], 'confidence': 'A'},
    'Tbx21': {'targets': ['Ifng', 'Il2', 'Ccl1', 'Ccl3', 'Ccl4', 'Ccl5', 'Stat4', 'Il12rb2', 'Cxcl10'], 'confidence': 'A'},
    'Ar': {'targets': ['Psa', 'Tmprss2', 'Klk2', 'Klk3', 'Nkx31', 'Pca3', 'Foxa1', 'Ctbp2', 'HOXB13'], 'confidence': 'A'},
    'Tcf7': {'targets': ['Ccnd1', 'Myc', 'Lef1', 'Axin2', 'Dkk1', 'Nodal', 'Nanog', 'Sox2', 'Oct4'], 'confidence': 'A'},
    'Ctnnb1': {'targets': ['Ccnd1', 'Myc', 'Axin2', 'Dkk1', 'Nodal', 'Sox2', 'Nanog', 'Tcf4', 'Lef1'], 'confidence': 'A'},
    'Cebpa': {'targets': ['Cebpb', 'Cebpd', 'Cd14', 'Cd33', 'Csf3', 'Csf1', 'Il6', 'Tnf', 'Ptgs2'], 'confidence': 'A'},
    'Cebpb': {'targets': ['Cebpa', 'Cebpd', 'Il6', 'Tnf', 'Il1b', 'Il8', 'Csf3', 'Lif', 'Socs3'], 'confidence': 'A'},
    'Runx1': {'targets': ['Runx2', 'Runx3', 'Cd19', 'Rag1', 'Rag2', 'Notch1', 'Il2', 'Il3', 'Il7'], 'confidence': 'A'},
    'Runx2': {'targets': ['Alpl', 'Bglap', 'Ibsp', 'Spp1', 'Col10a1', 'Rankl', 'Runx1', 'Runx3', 'Osterix'], 'confidence': 'A'},
    'Klf4': {'targets': ['Cdkn1a', 'Cdkn1b', 'Cdkn2b', 'Oct4', 'Sox2', 'Nanog', 'Myc', 'Postn', 'Tgfbi'], 'confidence': 'B'},
    'Sox2': {'targets': ['Oct4', 'Nanog', 'Fgf4', 'Hesx1', 'Lefty1', 'Lefty2', 'Nodal', 'Gdf3', 'Utf1'], 'confidence': 'A'},
    'Foxa1': {'targets': ['Esr1', 'Pgr', 'Klk2', 'Klk3', 'Gata3', 'Ar', 'Tmprss2', 'Fos', 'Cdkn1a'], 'confidence': 'A'},
    'Foxa2': {'targets': ['Pdx1', 'Hnf4a', 'Hnf1a', 'Neurod1', 'Neurog3', 'Ins', 'Gcg', 'Ttr', 'Apob'], 'confidence': 'A'},
}


def load_gene_list(genes_input):
    """Load genes from file or comma-separated string."""
    if os.path.isfile(genes_input):
        with open(genes_input, 'r') as f:
            genes = set(line.strip().upper() for line in f if line.strip())
    else:
        genes = set(g.strip().upper() for g in genes_input.split(',') if g.strip())
    return genes


def query_chea3(genes, query_name='query'):
    """Query ChEA3 API for TF enrichment."""
    try:
        url = 'https://maayanlab.cloud/chea3/api/enrich/'
        payload = {
            'query_name': query_name,
            'gene_set': list(genes),
        }
        resp = requests.post(url, json=payload, timeout=30)
        resp.raise_for_status()
        data = resp.json()

        results = {}
        # Parse response: typically returns dict of library -> list of [tf, pval, genes, rank]
        if isinstance(data, dict):
            for lib, tf_list in data.items():
                if isinstance(tf_list, list):
                    for item in tf_list[:50]:  # Top 50 per library
                        if isinstance(item, (list, tuple)) and len(item) >= 2:
                            tf = item[0].upper()
                            p_val = float(item[1]) if isinstance(item[1], (int, float)) else 1.0
                            if tf not in results:
                                results[tf] = {'p_values': [], 'sources': []}
                            results[tf]['p_values'].append(p_val)
                            results[tf]['sources'].append(lib)

        return results
    except Exception as e:
        print(f"Warning: ChEA3 API call failed: {e}", file=sys.stderr)
        return {}


def dorothea_enrichment(genes, species='human', confidence_levels='A,B'):
    """Run DoRothEA enrichment using built-in database."""
    db = DOROTHEA_HUMAN if species.lower() == 'human' else DOROTHEA_MOUSE
    confidence_set = set(confidence_levels.upper().split(','))

    results = {}

    for tf, info in db.items():
        if info['confidence'] not in confidence_set:
            continue

        targets = set(g.upper() for g in info['targets'])
        overlap = genes & targets
        results[tf.upper()] = {
            'targets': targets,
            'overlap': overlap,
            'overlap_count': len(overlap),
        }

    return results


def fisher_exact_test(a, b, c, d):
    """Fisher's exact test (hypergeometric)."""
    if a + b == 0 or c + d == 0:
        return 1.0
    if a + c == 0 or b + d == 0:
        return 1.0

    def comb(n, k):
        if k > n or k < 0:
            return 0
        if k == 0 or k == n:
            return 1
        k = min(k, n - k)
        result = 1
        for i in range(k):
            result = result * (n - i) // (i + 1)
        return result

    N = a + b + c + d
    K = a + c
    n = a + b

    p_value = 0.0
    for k in range(max(0, n + K - N), min(n, K) + 1):
        p = (comb(K, k) * comb(N - K, n - k)) / comb(N, n)
        if k >= a:
            p_value += p

    return min(p_value, 1.0)


def benjamini_hochberg_fdr(p_values):
    """Benjamini-Hochberg FDR correction."""
    n = len(p_values)
    if n == 0:
        return []

    sorted_indices = np.argsort(p_values)
    sorted_pvals = np.array(p_values)[sorted_indices]

    adjusted = np.zeros(n)
    for i, idx in enumerate(sorted_indices):
        adjusted[idx] = min(sorted_pvals[i] * n / (i + 1), 1.0)

    for i in range(n - 2, -1, -1):
        adjusted[i] = min(adjusted[i], adjusted[i + 1])

    return adjusted.tolist()


def merge_chea3_dorothea(chea3_results, dorothea_results, genes, background_size):
    """Merge ChEA3 and DoRothEA results."""
    all_tfs = set(chea3_results.keys()) | set(dorothea_results.keys())
    merged = {}

    for tf in all_tfs:
        chea3_pval = 1.0
        sources_chea3 = []
        if tf in chea3_results:
            p_vals = chea3_results[tf]['p_values']
            chea3_pval = np.mean(p_vals) if p_vals else 1.0
            sources_chea3 = chea3_results[tf]['sources']

        dorothea_pval = 1.0
        dorothea_overlap = 0
        overlap_genes = []
        targets = set()
        if tf in dorothea_results:
            targets = dorothea_results[tf]['targets']
            overlap = dorothea_results[tf]['overlap']
            dorothea_overlap = len(overlap)
            overlap_genes = sorted(list(overlap))

            # Compute p-value via Fisher's exact
            a = dorothea_overlap
            b = len(genes) - a
            c = len(targets) - a
            d = background_size - len(genes) - len(targets) + a
            dorothea_pval = fisher_exact_test(a, b, c, d)

        # Aggregate p-values: use minimum (most significant)
        combined_pval = min(chea3_pval, dorothea_pval) if chea3_pval < 1 or dorothea_pval < 1 else 1.0

        # Odds ratio
        if tf in dorothea_results:
            a = dorothea_overlap
            b = len(genes) - a
            c = len(targets) - a
            d = background_size - len(genes) - len(targets) + a
            or_val = (a * d) / (b * c) if b > 0 and c > 0 else float('inf')
        else:
            or_val = 0

        merged[tf] = {
            'pvalue': combined_pval,
            'chea3_sources': sources_chea3,
            'dorothea_overlap': dorothea_overlap,
            'overlap_genes': overlap_genes,
            'targets': targets,
            'odds_ratio': or_val,
            'from_chea3': tf in chea3_results,
            'from_dorothea': tf in dorothea_results,
        }

    return merged


def main():
    parser = argparse.ArgumentParser(description='Transcription Factor Enrichment Analysis')
    parser.add_argument('--genes', required=True, help='Gene list (file or comma-separated)')
    parser.add_argument('--species', choices=['human', 'mouse'], default='human', help='Species')
    parser.add_argument('--method', choices=['chea3', 'dorothea', 'both'], default='both')
    parser.add_argument('--confidence', default='A,B', help='DoRothEA confidence levels')
    parser.add_argument('--background', type=int, default=20000, help='Background gene universe size')
    parser.add_argument('--fdr-cutoff', type=float, default=0.05, help='FDR significance cutoff')
    parser.add_argument('--top-n', type=int, default=25, help='Top N TFs to report')
    parser.add_argument('--outdir', default='./tf_enrichment_output', help='Output directory')

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("Loading gene list...")
    genes = load_gene_list(args.genes)
    print(f"Loaded {len(genes)} genes")

    chea3_results = {}
    dorothea_results = {}

    if args.method in ['chea3', 'both']:
        print("Querying ChEA3 API...")
        chea3_results = query_chea3(genes)
        print(f"Found {len(chea3_results)} TFs from ChEA3")

    if args.method in ['dorothea', 'both']:
        print("Running DoRothEA enrichment...")
        dorothea_results = dorothea_enrichment(genes, args.species, args.confidence)
        print(f"Found {len(dorothea_results)} TFs from DoRothEA")

    print("Merging and computing statistics...")
    tf_results = merge_chea3_dorothea(chea3_results, dorothea_results, genes, args.background)

    if not tf_results:
        print("No TFs found in analysis", file=sys.stderr)
        return 1

    # Apply FDR correction
    p_values = [tf_results[tf]['pvalue'] for tf in tf_results]
    fdr_values = benjamini_hochberg_fdr(p_values)
    for i, tf in enumerate(sorted(tf_results.keys())):
        tf_results[tf]['fdr'] = fdr_values[i]

    # Filter by FDR and sort
    significant_tfs = [(tf, data) for tf, data in tf_results.items() if data['fdr'] <= args.fdr_cutoff]
    significant_tfs.sort(key=lambda x: x[1]['pvalue'])

    print(f"\nFound {len(significant_tfs)} significant TFs (FDR <= {args.fdr_cutoff})")
    print("\nTop 10 TFs:")
    for i, (tf, data) in enumerate(significant_tfs[:10]):
        print(f"  {i+1}. {tf}: p={data['pvalue']:.4e}, FDR={data['fdr']:.4f}, overlap={data['dorothea_overlap']}")

    # Save results
    output_list = []
    for tf in sorted(tf_results.keys(), key=lambda x: tf_results[x]['pvalue']):
        data = tf_results[tf]
        output_list.append({
            'TF': tf,
            'PValue': data['pvalue'],
            'FDR': data['fdr'],
            'OverlapCount': data['dorothea_overlap'],
            'OverlapGenes': ';'.join(data['overlap_genes']) if data['overlap_genes'] else '',
            'OddsRatio': round(data['odds_ratio'], 4) if data['odds_ratio'] != float('inf') else 'inf',
            'Source': ', '.join(['ChEA3'] if data['from_chea3'] else []) + (', DoRothEA' if data['from_dorothea'] else ''),
        })

    output_df = pd.DataFrame(output_list)
    output_df.to_csv(os.path.join(args.outdir, 'tf_enrichment.tsv'), sep='\t', index=False)
    print(f"\nSaved all TF results to tf_enrichment.tsv")

    # Visualizations
    top_tfs_to_plot = min(args.top_n, len(significant_tfs))
    if top_tfs_to_plot > 0:
        # Bar chart
        top_data = significant_tfs[:top_tfs_to_plot]
        tfs_plot = [tf for tf, _ in top_data]
        fdr_log = [-np.log10(max(data['fdr'], 1e-300)) for _, data in top_data]

        fig, ax = plt.subplots(figsize=(10, max(6, top_tfs_to_plot * 0.25)))
        colors = ['darkgreen' if fdr < -np.log10(0.05) else 'orange' for fdr in fdr_log]
        ax.barh(tfs_plot, fdr_log, color=colors, alpha=0.8, edgecolor='black')
        ax.set_xlabel('-log10(FDR)', fontsize=11)
        ax.set_title(f'Top {top_tfs_to_plot} Transcription Factors', fontsize=13, fontweight='bold')
        ax.axvline(-np.log10(0.05), color='red', linestyle='--', linewidth=1.5, label='FDR=0.05')
        ax.legend()
        ax.grid(axis='x', alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(args.outdir, 'top_tfs.png'), dpi=300, bbox_inches='tight')
        print("Saved top TFs bar chart to top_tfs.png")
        plt.close()

        # Network visualization (bipartite)
        top_10_tfs = significant_tfs[:min(10, len(significant_tfs))]
        all_target_genes = set()
        for tf, data in top_10_tfs:
            all_target_genes.update(data['overlap_genes'])

        if all_target_genes and len(top_10_tfs) > 0:
            fig, ax = plt.subplots(figsize=(12, 8))

            # Position TFs on left, genes on right
            tf_count = len(top_10_tfs)
            gene_count = len(all_target_genes)

            tf_positions = {tf: (0, i * 1.0) for i, (tf, _) in enumerate(top_10_tfs)}
            gene_positions = {gene: (3, i * 1.0) for i, gene in enumerate(sorted(all_target_genes))}

            # Draw edges
            for tf, data in top_10_tfs:
                for gene in data['overlap_genes']:
                    if gene in gene_positions:
                        x_vals = [tf_positions[tf][0], gene_positions[gene][0]]
                        y_vals = [tf_positions[tf][1], gene_positions[gene][1]]
                        ax.plot(x_vals, y_vals, 'k-', alpha=0.3, linewidth=0.8)

            # Draw nodes
            for tf, (x, y) in tf_positions.items():
                ax.scatter(x, y, s=500, c='red', alpha=0.7, edgecolors='darkred', linewidth=2, zorder=10)
                ax.text(x - 0.3, y, tf, ha='right', va='center', fontsize=9, fontweight='bold')

            for gene, (x, y) in gene_positions.items():
                ax.scatter(x, y, s=300, c='lightblue', alpha=0.7, edgecolors='darkblue', linewidth=1.5, zorder=10)
                ax.text(x + 0.15, y, gene, ha='left', va='center', fontsize=8)

            ax.set_xlim(-1, 4)
            ax.set_ylim(-1, max(tf_count, gene_count))
            ax.axis('off')
            ax.set_title('Top TF-Gene Regulatory Network', fontsize=13, fontweight='bold')

            plt.tight_layout()
            plt.savefig(os.path.join(args.outdir, 'tf_gene_network.png'), dpi=300, bbox_inches='tight')
            print("Saved TF-gene network to tf_gene_network.png")
            plt.close()

    # Summary report
    summary_path = os.path.join(args.outdir, 'tf_enrichment_summary.txt')
    with open(summary_path, 'w') as f:
        f.write("=== Transcription Factor Enrichment Analysis ===\n\n")
        f.write(f"Input genes: {len(genes)}\n")
        f.write(f"Total TFs tested: {len(tf_results)}\n")
        f.write(f"Significant TFs (FDR <= {args.fdr_cutoff}): {len(significant_tfs)}\n\n")
        f.write("Top 10 Enriched TFs:\n")
        for i, (tf, data) in enumerate(significant_tfs[:10]):
            f.write(f"{i+1}. {tf}\n")
            f.write(f"   p-value: {data['pvalue']:.4e}\n")
            f.write(f"   FDR: {data['fdr']:.4f}\n")
            f.write(f"   Target overlap: {data['dorothea_overlap']} genes\n")
            if data['overlap_genes']:
                f.write(f"   Overlapping genes: {', '.join(data['overlap_genes'][:5])}")
                if len(data['overlap_genes']) > 5:
                    f.write(f", +{len(data['overlap_genes']) - 5} more")
                f.write("\n")
            f.write("\n")

    print(f"Saved summary to tf_enrichment_summary.txt")
    print(f"\nCompleted! Results saved to {args.outdir}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
