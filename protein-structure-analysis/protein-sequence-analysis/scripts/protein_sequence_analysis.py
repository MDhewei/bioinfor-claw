#!/usr/bin/env python3
"""
Protein Sequence Analysis
Analyze sequences for properties, motifs, PTMs, and disorder
"""

import argparse
import sys
import os
import re
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import requests
import warnings
warnings.filterwarnings('ignore')

# Amino acid properties
AA_WEIGHTS = {
    'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.15,
    'Q': 146.15, 'E': 147.13, 'G': 75.07, 'H': 155.16, 'I': 131.18,
    'L': 131.18, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
    'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15
}

AA_CHARGE = {
    'A': 0, 'R': 1, 'N': 0, 'D': -1, 'C': 0,
    'Q': 0, 'E': -1, 'G': 0, 'H': 0.1, 'I': 0,
    'L': 0, 'K': 1, 'M': 0, 'F': 0, 'P': 0,
    'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0
}

KYTE_DOOLITTLE = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

INSTABILITY_INDEX_FACTORS = {
    ('A', 'A'): 0.0, ('A', 'R'): 0.0, ('A', 'N'): 0.0, ('A', 'D'): 0.23, ('A', 'C'): 0.0,
    ('R', 'A'): 0.0, ('R', 'N'): 0.0, ('R', 'D'): 0.22, ('R', 'K'): 0.0,
    ('D', 'A'): 0.23, ('E', 'A'): 0.37, ('L', 'A'): 0.13,
}

def fetch_uniprot_sequence(gene, species='human'):
    """Fetch protein sequence from UniProt"""
    species_id = {'human': 9606, 'mouse': 10090}
    organism_id = species_id.get(species, 9606)

    url = 'https://rest.uniprot.org/uniprotkb/search'
    params = {
        'query': f'gene_exact:{gene} AND organism_id:{organism_id} AND reviewed:true',
        'format': 'fasta',
        'size': 1
    }

    try:
        response = requests.get(url, params=params, timeout=10)
        if response.status_code == 200:
            fasta = response.text.strip()
            lines = fasta.split('\n')
            if len(lines) >= 2:
                header = lines[0]
                sequence = ''.join(lines[1:])
                accession = header.split('|')[1] if '|' in header else 'unknown'
                return sequence.upper(), accession
    except Exception as e:
        print(f"Error fetching from UniProt: {e}")

    return None, None

def fetch_uniprot_features(uniprot_id):
    """Fetch features from UniProt JSON API"""
    url = f'https://rest.uniprot.org/uniprotkb/{uniprot_id}.json'

    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            features = []
            if 'features' in data:
                for feature in data['features']:
                    ftype = feature.get('type', '')
                    location = feature.get('location', {})
                    start = location.get('start', {}).get('value')
                    end = location.get('end', {}).get('value')
                    if start and end:
                        features.append({
                            'type': ftype,
                            'start': start,
                            'end': end,
                            'description': feature.get('description', '')
                        })
            return features
    except Exception as e:
        print(f"Error fetching UniProt features: {e}")

    return []

def read_fasta(fasta_file):
    """Read FASTA file"""
    with open(fasta_file, 'r') as f:
        lines = f.readlines()

    header = lines[0].strip()
    sequence = ''.join([line.strip() for line in lines[1:]])
    return sequence.upper(), header

def validate_sequence(seq):
    """Check if sequence contains only standard amino acids"""
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    return all(aa in valid_aa for aa in seq.upper())

def compute_molecular_weight(seq):
    """Compute protein molecular weight"""
    mw = sum(AA_WEIGHTS.get(aa, 0) for aa in seq) - (len(seq) - 1) * 18.015
    return mw

def compute_isoelectric_point(seq):
    """Estimate isoelectric point via binary search"""
    def net_charge(seq, pH):
        """Compute net charge at given pH"""
        charge = 0
        for aa in seq:
            if aa == 'K':
                charge += 1.0 / (1.0 + 10 ** (3.5 - pH))
            elif aa == 'R':
                charge += 1.0 / (1.0 + 10 ** (3.1 - pH))
            elif aa == 'D':
                charge -= 1.0 / (1.0 + 10 ** (pH - 4.0))
            elif aa == 'E':
                charge -= 1.0 / (1.0 + 10 ** (pH - 4.2))
        return charge

    # Binary search for pH where net_charge ~ 0
    pH_low, pH_high = 2.0, 12.0
    for _ in range(50):
        pH_mid = (pH_low + pH_high) / 2
        if net_charge(seq, pH_mid) > 0:
            pH_low = pH_mid
        else:
            pH_high = pH_mid

    return (pH_low + pH_high) / 2

def compute_gravy(seq):
    """Compute Grand Average of Hydropathy"""
    gravy = np.mean([KYTE_DOOLITTLE.get(aa, 0) for aa in seq])
    return gravy

def compute_instability_index(seq):
    """Estimate instability index"""
    ii = 0
    for i in range(len(seq) - 1):
        pair = (seq[i], seq[i + 1])
        ii += INSTABILITY_INDEX_FACTORS.get(pair, 0)

    ii = (ii / (len(seq) - 1)) * 100
    return ii

def compute_secondary_structure_propensity(seq):
    """Estimate secondary structure propensity"""
    chou_fasman_helix = {
        'A': 1.42, 'R': 0.98, 'N': 0.67, 'D': 1.01, 'C': 0.70,
        'Q': 1.11, 'E': 1.51, 'G': 0.57, 'H': 1.00, 'I': 1.08,
        'L': 1.21, 'K': 1.16, 'M': 1.45, 'F': 1.13, 'P': 0.57,
        'S': 0.77, 'T': 0.83, 'W': 1.08, 'Y': 0.69, 'V': 1.06
    }

    helix_prob = np.mean([chou_fasman_helix.get(aa, 0) for aa in seq])
    return helix_prob

def scan_motifs(seq):
    """Scan for functional motifs"""
    motifs = {}

    # NLS patterns
    nls_pattern = r'K[KR]{2}|KR'
    for match in re.finditer(nls_pattern, seq):
        if 'NLS' not in motifs:
            motifs['NLS'] = []
        motifs['NLS'].append((match.start(), match.end(), match.group()))

    # NES pattern
    nes_pattern = r'L.{2,3}[LIVMF].{2,3}[LIVMF].L'
    for match in re.finditer(nes_pattern, seq):
        if 'NES' not in motifs:
            motifs['NES'] = []
        motifs['NES'].append((match.start(), match.end(), match.group()))

    # RGD
    rgd_pattern = r'RGD'
    for match in re.finditer(rgd_pattern, seq):
        if 'RGD' not in motifs:
            motifs['RGD'] = []
        motifs['RGD'].append((match.start(), match.end(), match.group()))

    # CAAX prenylation
    caax_pattern = r'C[AVILM]{2}[AVILM]$'
    for match in re.finditer(caax_pattern, seq):
        if 'CAAX' not in motifs:
            motifs['CAAX'] = []
        motifs['CAAX'].append((match.start(), match.end(), match.group()))

    # N-glycosylation
    nglyc_pattern = r'N[^P][ST]'
    for match in re.finditer(nglyc_pattern, seq):
        if 'N-glycosylation' not in motifs:
            motifs['N-glycosylation'] = []
        motifs['N-glycosylation'].append((match.start(), match.end(), match.group()))

    # SH2 binding
    sh2_pattern = r'Y[A-Z]{2}[LIVMF]'
    for match in re.finditer(sh2_pattern, seq):
        if 'SH2 binding' not in motifs:
            motifs['SH2 binding'] = []
        motifs['SH2 binding'].append((match.start(), match.end(), match.group()))

    # SH3 binding
    sh3_pattern = r'P.{2}P'
    for match in re.finditer(sh3_pattern, seq):
        if 'SH3 binding' not in motifs:
            motifs['SH3 binding'] = []
        motifs['SH3 binding'].append((match.start(), match.end(), match.group()))

    return motifs

def predict_disorder(seq, window_size=20):
    """Predict disordered regions using hydrophobicity + charge"""
    disorder_scores = []
    n = len(seq)

    for i in range(n):
        start = max(0, i - window_size // 2)
        end = min(n, i + window_size // 2)
        window = seq[start:end]

        # Hydrophobicity
        hydrophobicity = np.mean([KYTE_DOOLITTLE.get(aa, 0) for aa in window])

        # Charge
        charge = sum(AA_CHARGE.get(aa, 0) for aa in window) / len(window)

        # Simple disorder heuristic: low hydrophobicity + high charge
        if hydrophobicity < -0.5 and abs(charge) > 0.3:
            disorder_score = 1.0
        else:
            disorder_score = 0.0

        disorder_scores.append(disorder_score)

    return np.array(disorder_scores)

def plot_feature_map(seq, features=None, motifs=None, disorder_scores=None, outfile="feature_map.png"):
    """Plot annotated protein feature map"""
    seq_len = len(seq)

    fig, ax = plt.subplots(figsize=(15, 8), dpi=300)

    # Draw main protein bar
    ax.add_patch(mpatches.Rectangle((0, 5), seq_len, 1, facecolor='lightgray', edgecolor='black'))

    # Draw features
    track_y = 7
    if features:
        for feature in features:
            ftype = feature.get('type', '')
            start = feature.get('start', 0)
            end = feature.get('end', 0)

            if ftype == 'SIGNAL':
                ax.add_patch(mpatches.Rectangle((start, track_y), end - start, 0.8,
                                               facecolor='red', alpha=0.7, edgecolor='black'))
            elif ftype == 'TRANSMEM':
                ax.add_patch(mpatches.Rectangle((start, track_y), end - start, 0.8,
                                               facecolor='blue', alpha=0.7, edgecolor='black'))
            elif ftype == 'DOMAIN':
                ax.add_patch(mpatches.Rectangle((start, track_y), end - start, 0.8,
                                               facecolor='green', alpha=0.7, edgecolor='black'))

        track_y += 1.5

    # Draw motifs
    if motifs:
        for motif_name, matches in motifs.items():
            for start, end, _ in matches:
                ax.plot([start, end], [track_y, track_y], 'o-', markersize=5, color='purple')
        track_y += 1.5

    # Draw disorder regions
    if disorder_scores is not None:
        for i, score in enumerate(disorder_scores):
            if score > 0.5:
                ax.add_patch(mpatches.Rectangle((i, track_y - 0.4), 1, 0.8,
                                               facecolor='yellow', alpha=0.5))

    ax.set_xlim(0, seq_len)
    ax.set_ylim(0, track_y + 1)
    ax.set_xlabel('Residue Position', fontsize=12)
    ax.set_title('Protein Feature Map', fontsize=14, fontweight='bold')
    ax.set_yticks([])

    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    print(f"Saved {outfile}")
    plt.close()

def plot_aa_composition(seq, outfile="aa_composition.png"):
    """Plot amino acid composition"""
    aa_counts = {}
    for aa in seq:
        aa_counts[aa] = aa_counts.get(aa, 0) + 1

    # Sort by frequency
    sorted_aa = sorted(aa_counts.items(), key=lambda x: x[1], reverse=True)
    aa_names = [x[0] for x in sorted_aa]
    counts = [x[1] for x in sorted_aa]
    fractions = [c / len(seq) for c in counts]

    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)

    bars = ax.bar(aa_names, fractions, color='steelblue', alpha=0.7, edgecolor='black')

    ax.set_xlabel('Amino Acid', fontsize=12)
    ax.set_ylabel('Fraction', fontsize=12)
    ax.set_title('Amino Acid Composition', fontsize=14, fontweight='bold')
    ax.set_ylim(0, max(fractions) * 1.1)

    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    print(f"Saved {outfile}")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Protein Sequence Analysis')
    parser.add_argument('--gene', help='Gene symbol (e.g., TP53)')
    parser.add_argument('--species', choices=['human', 'mouse'], default='human')
    parser.add_argument('--uniprot-id', help='UniProt accession')
    parser.add_argument('--fasta', help='Local FASTA file')
    parser.add_argument('--outdir', default='results')
    parser.add_argument('--compute-properties', action='store_true')
    parser.add_argument('--find-motifs', action='store_true')
    parser.add_argument('--predict-disorder', action='store_true')
    parser.add_argument('--find-ptm-sites', action='store_true')
    parser.add_argument('--conservation-file', help='Conservation scores TSV')

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Load sequence
    sequence = None
    uniprot_id = args.uniprot_id

    if args.fasta:
        sequence, header = read_fasta(args.fasta)
        print(f"Loaded sequence from {args.fasta} ({len(sequence)} aa)")
    elif args.gene:
        sequence, uniprot_id = fetch_uniprot_sequence(args.gene, args.species)
        if sequence:
            print(f"Fetched {args.gene} from UniProt ({len(sequence)} aa, ID: {uniprot_id})")
        else:
            print(f"Could not fetch {args.gene}. Try using --fasta or --uniprot-id.")
            sys.exit(1)

    if sequence is None:
        print("No sequence provided. Use --gene, --fasta, or --uniprot-id.")
        sys.exit(1)

    if not validate_sequence(sequence):
        print("Sequence contains non-standard amino acids.")
        sys.exit(1)

    # Compute properties
    properties = {}
    if args.compute_properties:
        properties['Sequence Length'] = len(sequence)
        properties['Molecular Weight (Da)'] = compute_molecular_weight(sequence)
        properties['Isoelectric Point (pI)'] = compute_isoelectric_point(sequence)
        properties['GRAVY'] = compute_gravy(sequence)
        properties['Instability Index'] = compute_instability_index(sequence)
        properties['Helix Propensity'] = compute_secondary_structure_propensity(sequence)

        with open(f'{args.outdir}/sequence_properties.txt', 'w') as f:
            for key, val in properties.items():
                if isinstance(val, float):
                    f.write(f"{key}: {val:.3f}\n")
                else:
                    f.write(f"{key}: {val}\n")
        print(f"Saved properties to {args.outdir}/sequence_properties.txt")

    # Find motifs
    motifs = {}
    if args.find_motifs:
        motifs = scan_motifs(sequence)

        motif_rows = []
        for motif_name, matches in motifs.items():
            for start, end, seq_match in matches:
                motif_rows.append({
                    'motif': motif_name,
                    'start': start,
                    'end': end,
                    'sequence': seq_match
                })

        if motif_rows:
            motif_df = pd.DataFrame(motif_rows)
            motif_df.to_csv(f'{args.outdir}/motifs.tsv', sep='\t', index=False)
            print(f"Found {len(motif_df)} motifs. Saved to {args.outdir}/motifs.tsv")

    # Predict disorder
    disorder_scores = None
    if args.predict_disorder:
        disorder_scores = predict_disorder(sequence)

        disorder_rows = []
        for i, score in enumerate(disorder_scores):
            if score > 0.5:
                disorder_rows.append({
                    'position': i + 1,
                    'disorder_score': score
                })

        if disorder_rows:
            disorder_df = pd.DataFrame(disorder_rows)
            disorder_df.to_csv(f'{args.outdir}/disorder_regions.tsv', sep='\t', index=False)
            print(f"Predicted {len(disorder_df)} disordered positions")

    # Fetch UniProt features
    features = []
    if args.find_ptm_sites and uniprot_id:
        features = fetch_uniprot_features(uniprot_id)

        feature_rows = []
        for feature in features:
            feature_rows.append({
                'type': feature['type'],
                'start': feature['start'],
                'end': feature['end'],
                'description': feature.get('description', '')
            })

        if feature_rows:
            feature_df = pd.DataFrame(feature_rows)
            feature_df.to_csv(f'{args.outdir}/features.tsv', sep='\t', index=False)
            print(f"Fetched {len(feature_df)} UniProt features")

    # Generate plots
    plot_feature_map(sequence, features, motifs, disorder_scores,
                     f'{args.outdir}/feature_map.png')
    plot_aa_composition(sequence, f'{args.outdir}/aa_composition.png')

    print("Done!")

if __name__ == '__main__':
    main()
