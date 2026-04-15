---
name: mutation-analysis-for-gene
description: Analyze somatic mutations in a gene across TCGA cancer cohorts. Fetches mutation data from the GDC/TCGA API, computes mutation frequencies, identifies hotspot residues, classifies mutation types, and generates lollipop plots and oncoprint-style summaries.
---

# Mutation Analysis for Gene

## Purpose
Characterize the mutational landscape of a specific gene across TCGA cancer cohorts, including:
- Frequency of mutations in each cancer type
- Identification of mutational hotspots (recurrently mutated amino acids)
- Classification of mutation consequence types (missense, nonsense, frameshift, splice site, etc.)
- Visual summaries: lollipop plots, mutation frequency bar charts, and pie charts

## Use When
- You want to understand which amino acid positions in a gene are frequently mutated
- You need to compare mutation frequencies across cancer types
- You're looking for actionable mutations (hotspots) for drug target selection
- You want to characterize the mutational spectrum of a candidate gene
- You need to generate publication-quality visualizations of protein mutation patterns

## Do Not Use When
- You need germline variant analysis (this tool is for somatic mutations only)
- You want to perform rare variant association studies (focuses on recurrent mutations)
- You need RNA-level mutation impact prediction (use VEP/SIFT for that)
- You lack internet access (requires GDC API access)

## Expected Inputs & Outputs

### Inputs
- **--gene** (required): Gene symbol (e.g., TP53, BRCA1, KRAS)
- **--cancer-types**: Comma-separated TCGA project codes (e.g., TCGA-BRCA,TCGA-LUAD) or "all" (default: all)
- **--mutation-types**: Comma-separated consequence types to include (default: all except Silent)
  - Options: Missense_Mutation, Nonsense_Mutation, Frame_Shift_Del, Frame_Shift_Ins, Splice_Site, In_Frame_Del, In_Frame_Ins, Silent
- **--min-frequency**: Minimum mutation frequency (0.0-1.0) to label on plots (default: 0.01)
- **--top-hotspots**: Number of top hotspot codons to highlight (default: 10)
- **--protein-length**: Protein length in amino acids (optional; fetches from UniProt if not provided)
- **--domain-file**: Optional TSV file with columns: domain_name, start_aa, end_aa
- **--outdir**: Output directory for results (default: current directory)

### Outputs
- `mutation_summary.tsv`: Summary table (cancer_type, n_mutations, n_cases, frequency, top_hotspots)
- `lollipop_plot.png`: Protein mutation lollipop plot with domain tracks
- `mutation_frequency.png`: Bar chart of mutation frequency by cancer type
- `mutation_types.png`: Pie chart of mutation type distribution
- `hotspot_details.tsv`: Detailed hotspot information (codon, amino_acid, position, count, types)

## Procedure

### Step 1: Query GDC Mutations API
1. Construct filter query for the specified gene symbol
2. Request mutational data from GDC with pagination (1000 per batch)
3. Extract: ssm_id, case.project, aa_change, consequence_type, genomic position

### Step 2: Compute Statistics
1. Count mutations per cancer type
2. Calculate mutation frequency (n_mutated_cases / n_total_cases)
3. Identify hotspot residues (amino acids with ≥3 mutations)
4. Classify mutation types

### Step 3: Fetch Protein Information
1. Query UniProt API for protein length and canonical sequence
2. If domain file provided, parse and validate domain coordinates

### Step 4: Generate Visualizations
1. **Lollipop plot**:
   - X-axis: amino acid position (0 to protein_length)
   - Y-axis: mutation count at each position
   - Circle markers colored by mutation type (Missense=blue, Nonsense=red, Frameshift=green, Splice=orange)
   - Domain tracks as colored rectangles below x-axis
   - Top hotspots labeled with residue and count

2. **Mutation frequency bar chart**: Cancer types sorted by frequency

3. **Mutation type pie chart**: Proportions of consequence types

## Key Execution Patterns

### Basic usage (all cancer types):
```bash
python scripts/mutation_analysis_for_gene.py \
  --gene TP53 \
  --outdir results/tp53
```

### Specific cancer types:
```bash
python scripts/mutation_analysis_for_gene.py \
  --gene KRAS \
  --cancer-types TCGA-BRCA,TCGA-LUAD,TCGA-COAD \
  --outdir results/kras_subset
```

### With protein length and domains:
```bash
python scripts/mutation_analysis_for_gene.py \
  --gene BRCA1 \
  --protein-length 1863 \
  --domain-file brca1_domains.tsv \
  --min-frequency 0.005 \
  --top-hotspots 15 \
  --outdir results/brca1
```

### Specific mutation types only:
```bash
python scripts/mutation_analysis_for_gene.py \
  --gene MYC \
  --mutation-types Missense_Mutation,Nonsense_Mutation \
  --outdir results/myc_coding
```

## Parameter Decision Guide

| Scenario | --gene | --cancer-types | --mutation-types | --protein-length | --domain-file |
|----------|--------|-----------------|------------------|------------------|---------------|
| Pan-TCGA analysis | TP53 | all | all (default) | (auto-fetch) | None |
| Cancer-specific | BRCA1 | TCGA-BRCA | all | (auto-fetch) | brca1.tsv |
| High-frequency only | KRAS | TCGA-LUAD,TCGA-COAD,TCGA-PAAD | Missense_Mutation | 189 | None |
| Coding variants | MYC | all | Missense_Mutation,Nonsense_Mutation,Frame_Shift_Del,Frame_Shift_Ins | 439 | None |

## Failure Modes

| Error | Cause | Solution |
|-------|-------|----------|
| "Gene not found in GDC" | Gene symbol not recognized or not in TCGA | Verify gene symbol (case-sensitive); try HGNC official symbol |
| "No mutations found" | Gene is rarely mutated in specified cohort | Lower --min-frequency, expand --cancer-types |
| "Connection timeout" | GDC API unreachable | Check internet connection; retry (GDC may be temporarily down) |
| "Invalid protein length" | Protein has multiple isoforms | Manually specify --protein-length using canonical isoform |
| "Domain file not found" | Path to domain file incorrect | Verify TSV path and format (tab-delimited with headers) |
| "UniProt lookup failed" | Gene symbol not recognized by UniProt | Provide explicit --protein-length |

## Technical Notes

- Mutation frequencies are computed as: (number of unique cases with mutation) / (total number of cases in project)
- Hotspots are defined as amino acid positions with ≥3 independent mutations
- Consequence types follow VEP (Variant Effect Predictor) nomenclature
- Domain tracks are optional; if provided, must be TSV with columns: domain_name, start_aa, end_aa
- All plots use matplotlib; high-resolution PNG output (300 DPI)
- Results are reproducible; same gene/cancer-type combo produces identical output

