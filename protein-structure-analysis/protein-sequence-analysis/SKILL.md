---
name: protein-sequence-analysis
description: Analyze protein sequences for functional features. Fetches sequence from UniProt, computes physicochemical properties, identifies conserved motifs, PTM sites, disordered regions, signal peptides, transmembrane helices, and generates annotated sequence feature maps.
---

# Protein Sequence Analysis

## Purpose

Characterize protein sequences by computing properties, predicting functional domains/motifs, and identifying post-translational modifications. Enables rapid annotation of newly discovered proteins, assessment of evolutionary constraints, and functional feature prediction.

## Use When

- Characterizing a newly sequenced protein
- Predicting functional domains, motifs, signal peptides, transmembrane topology
- Assessing intrinsic disorder / structural features
- Planning experimental design (e.g., where to place fluorescent tags)
- Comparing proteins across species
- Building feature vectors for machine learning

## Do Not Use When

- Precise 3D structure is critical (use structure prediction tools)
- Sequence length < 5 amino acids
- Intrinsically disordered regions require high-resolution disorder scoring (use full IUPred/DSSP)
- De novo domain discovery needed (use profile HMM searches)

## Expected Inputs

### Required (one of)
- **--gene**: Gene symbol (e.g., "TP53", "BRCA1"). Will fetch sequence from UniProt API.
- **--fasta**: Local FASTA file with protein sequence(s).

### Optional
- **--species**: Human (9606), mouse (10090), other; default: human
- **--uniprot-id**: Override automatic UniProt lookup with explicit accession (e.g., "P04637")
- **--conservation-file**: TSV with position, conservation_score columns
- **--compute-properties**: Flag to enable physicochemical properties
- **--find-motifs**: Flag to scan for functional motifs
- **--predict-disorder**: Flag to predict disordered regions
- **--find-ptm-sites**: Flag to identify known PTM sites

## Expected Outputs

- **sequence_properties.txt**: All computed physicochemical properties (MW, pI, GRAVY, instability, etc.)
- **features.tsv**: All annotated features (signal peptide, transmembrane, domains, active sites, binding sites) with start/end positions
- **motifs.tsv**: All identified motifs (NLS, NES, RGD, etc.) with positions and matches
- **disorder_regions.tsv**: Predicted disordered regions with disorder score (if --predict-disorder)
- **feature_map.png**: Linear protein diagram with all features annotated
- **aa_composition.png**: Bar chart of amino acid frequency

## Procedure

1. **Fetch sequence**:
   - If --fasta provided: read local FASTA file, parse sequence and header
   - Else if --uniprot-id provided: fetch from UniProt REST API
   - Else: search UniProt by --gene and --species, retrieve first high-confidence match
   - Validate: sequence is protein (only standard amino acids + gap characters)
2. **Fetch UniProt annotations**:
   - Call UniProt REST API to retrieve features (signal peptide, transmembrane helices, domains, PTM sites, active sites, binding sites)
   - Parse JSON response, extract feature type and position
3. **Compute physicochemical properties**:
   - **Molecular Weight**: Sum of standard amino acid masses (MW lookup table)
   - **Isoelectric Point (pI)**: Binary search to find pH where net charge = 0 (Henderson-Hasselbalch)
   - **GRAVY**: (Sum of hydropathy scores per AA) / sequence length (Kyte-Doolittle scale)
   - **Instability Index**: GRAVY-based formula (Guruprasad et al. 1990)
   - **Secondary structure propensity**: Fraction helix/sheet/coil based on Chou-Fasman preferences
   - **Aliphatic Index**: Related to hydrophobic AA content
   - **Aromaticity**: Fraction of aromatic amino acids (F, W, Y)
4. **Scan for motifs** (regex pattern matching):
   - **NLS** (nuclear localization signal): K[KR]{2}X{0,4}K, or KR pattern
   - **NES** (nuclear export signal): [LI]X{2}[LVI]X{2}[LVI]L
   - **RGD** (integrin binding): RGD
   - **CAAX** (prenylation): C[AVILM]{2}[AVILM]$
   - **N-glycosylation**: N[^P][ST]
   - **SH2 binding**: Y[A-Z]{2}[LIVMF]
   - **SH3 binding**: P.{2}P
   - Report all matches with position and context
5. **Predict disordered regions** (simplified heuristic):
   - Sliding window (size = 20 residues)
   - Compute window hydrophobicity (Kyte-Doolittle) and charge (K+R-D-E)
   - Disordered if: hydrophobicity < -0.5 AND |charge| > 0.3
   - Output regions with disorder score (0-1 scale)
6. **Generate feature map**:
   - Linear diagram: x-axis = residue position (1 to length), y-axis = feature tracks
   - Tracks (from top): Signal peptide, Transmembrane helices (colored regions), Domains (colored boxes), PTM sites (dots), Motifs (arrows), Disorder regions (gray shading)
   - Sequence beneath with residue numbers
   - Add conservation track if --conservation-file provided
7. **Generate AA composition plot**:
   - Bar chart: 20 amino acids (or 21 with ambiguous)
   - Y-axis = fraction, sorted by frequency

## Key Execution Patterns

```bash
# Fetch TP53, compute all properties
python protein_sequence_analysis.py \
  --gene TP53 \
  --species human \
  --compute-properties \
  --find-motifs \
  --predict-disorder \
  --find-ptm-sites \
  --outdir results/

# Local FASTA, minimal analysis
python protein_sequence_analysis.py \
  --fasta my_protein.fasta \
  --compute-properties \
  --outdir results/

# UniProt accession override with conservation
python protein_sequence_analysis.py \
  --uniprot-id P12345 \
  --compute-properties \
  --find-motifs \
  --conservation-file conservation.tsv \
  --outdir results/

# Mouse protein, full analysis
python protein_sequence_analysis.py \
  --gene Tp53 \
  --species mouse \
  --compute-properties \
  --find-motifs \
  --predict-disorder \
  --find-ptm-sites \
  --outdir results/

# Quick motif scan only
python protein_sequence_analysis.py \
  --gene BRCA1 \
  --find-motifs \
  --outdir results/
```

## Parameter Decision Guide

| Scenario | Key Flags | Rationale |
|----------|-----------|-----------|
| First look at new protein | All flags enabled | Comprehensive initial characterization |
| Planning fluorescent tag placement | `--find-motifs`, `--predict-disorder`, `--find-ptm-sites` | Avoid tags at functional sites or disordered regions |
| Comparative analysis across species | `--gene` + `--species` | UniProt lookup works for all model organisms |
| Localization prediction | `--find-motifs` | NLS/NES motifs predict nuclear trafficking |
| Structure-function study | `--predict-disorder`, `--compute-properties` | Disorder correlates with conformational flexibility |
| Publication figure | `--find-motifs`, `--predict-disorder`, `--conservation-file` | Feature map is publication-ready |
| High-throughput batch analysis | `--fasta` (local files) | Faster than API calls; no network dependency |

## Failure Modes

| Error | Cause | Solution |
|-------|-------|----------|
| "Gene symbol not found in UniProt" | Gene name not recognized or organism not in database | Try --uniprot-id with explicit accession or use --fasta |
| "UniProt API request failed" | Network issue or API temporarily unavailable | Retry; check internet connection; use local --fasta as alternative |
| "Sequence contains non-standard amino acids" | FASTA file contains 'X', 'U', 'Z', or other ambiguous characters | Filter to standard 20 amino acids or validate FASTA format |
| "Feature map is crowded/unreadable" | Too many overlapping features | Split into separate subplot tracks or filter feature types |
| "Disorder prediction all 0 or all 1" | Heuristic not tuned for sequence composition | Use full IUPred tool for detailed disorder scoring |
| "No motifs found" | Motif patterns too stringent; rare in query protein | Expected; not all proteins contain all motif types |
| "Conservation track missing in plot" | --conservation-file provided but not used | Verify TSV format: position (int), conservation_score (float) |
| "AA composition bar labels overlap" | Figure too narrow | Rotate labels 45 degrees or increase figure width (code default: 10 inches) |
| "Memory error on very large sequence" | Sequence > 10,000 aa; unusual but possible | Pre-filter regions of interest or split analysis |

