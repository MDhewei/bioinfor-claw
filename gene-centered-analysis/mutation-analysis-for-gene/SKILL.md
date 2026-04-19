---
name: mutation-analysis-for-gene
description: "Analyze TCGA somatic mutations in a gene via cBioPortal API (PanCancer Atlas). Supports local TCGA MAF file as fallback. Computes mutation frequencies, identifies hotspot residues, classifies mutation types, and generates lollipop plots and summaries."
---

# Mutation Analysis for Gene

## Purpose
Characterize the mutational landscape of a specific gene, including:
- Frequency of mutations across cancer types or cell lines
- Identification of mutational hotspots (recurrently mutated amino acids)
- Classification of mutation consequence types (missense, nonsense, frameshift, splice site, etc.)
- Visual summaries: lollipop plots, mutation frequency bar charts, and pie charts

## Data sources (in priority order)
1. **cBioPortal API (default)** — Fetches TCGA PanCancer Atlas somatic mutations from `www.cbioportal.org`. Requires internet access.
2. **`--mutation-file` (local file)** — TCGA MAF or any TSV with Hugo_Symbol column. Use as fallback if cBioPortal is unreachable.

**Do NOT include `--mutation-file` in EXEC_ARGS unless the user explicitly provides a local MAF file.** The script queries cBioPortal by default.

## Inputs
- **--gene** (required): Gene symbol (e.g., TP53, BRCA1, KRAS)
- **--mutation-file**: Path to local mutation file (DepMap CSV or TCGA MAF). Auto-detected from cache if omitted.
- **--cancer-types**: Comma-separated TCGA project codes or "all" (default: all). Used with cBioPortal API.
- **--mutation-types**: Comma-separated consequence types to include (default: all)
- **--min-frequency**: Minimum mutation frequency to label (default: 0.01)
- **--top-hotspots**: Number of top hotspot codons to highlight (default: 10)
- **--protein-length**: Protein length in amino acids (fetches from UniProt if omitted)
- **--domain-file**: Optional TSV with columns: domain_name, start_aa, end_aa
- **--outdir**: Output directory

## Outputs
- `mutation_summary.tsv`: Summary table (cancer_type/cell_line, n_mutations, frequency, top_hotspots)
- `lollipop_plot.png`: Protein mutation lollipop plot with optional domain tracks
- `mutation_frequency.png`: Bar chart of mutation frequency by cancer type / cell line group
- `mutation_types.png`: Pie chart of mutation type distribution
- `hotspot_details.tsv`: Detailed hotspot information (position, count, types)

## Example commands

### Basic (auto-detects DepMap cache):
```bash
python scripts/mutation_analysis_for_gene.py \
  --gene KRAS \
  --outdir results/kras
```

### With explicit DepMap file:
```bash
python scripts/mutation_analysis_for_gene.py \
  --gene TP53 \
  --mutation-file web_results/depmap_cache/OmicsSomaticMutations.csv \
  --outdir results/tp53
```

### With TCGA MAF file:
```bash
python scripts/mutation_analysis_for_gene.py \
  --gene BRCA1 \
  --mutation-file combined_mutations.maf \
  --protein-length 1863 \
  --outdir results/brca1
```

## Parameter Decision Guide

| Signal in user request | Parameter to set |
|---|---|
| "mutation landscape of X" | `--gene X` (auto-detect data) |
| "mutations in DepMap cell lines" | `--gene X` (uses DepMap cache) |
| "TCGA mutations" + has MAF file | `--gene X --mutation-file <MAF>` |
| "specific cancer types only" | `--cancer-types BRCA,LUAD` (GDC API only) |
| "only missense mutations" | `--mutation-types Missense_Mutation` |

## Failure Modes

| Error | Cause | Solution |
|-------|-------|----------|
| "No mutations found" | Gene not in data source | Check gene symbol; try different data source |
| "cBioPortal API unreachable" | Network blocked | Use `--mutation-file` with local MAF data |
| "Cannot find gene column" | Non-standard file format | Ensure file has HugoSymbol or Hugo_Symbol column |
| "UniProt lookup failed" | API unreachable | Provide `--protein-length` manually |

## Technical Notes

- Supports both VEP SO terms (missense_variant) and MAF names (Missense_Mutation)
- Hotspots defined as positions with ≥3 independent mutations
- Data source: TCGA patient tumors via cBioPortal PanCancer Atlas (32 cancer types)
- Protein length auto-fetched from UniProt when possible; fallback lengths for common genes built in
