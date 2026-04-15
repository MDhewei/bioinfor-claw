---
name: design-prime-editor-sgrnas
description: Design pegRNAs (prime editing guide RNAs) and nicking sgRNAs for prime editing. Given a target gene and desired edit (SNV, small insertion, or deletion), identifies PE-compatible protospacers, designs RT templates and PBS sequences, and scores guides by predicted efficiency.
tags:
  - crispr
  - prime-editing
  - guide-design
  - pegRNA
version: 1.0.0
---

## Purpose

Design prime editor (PE) sgRNAs and pegRNAs for precise genomic edits without double-strand breaks. This tool automates identification of PE-compatible protospacers, design of reverse-transcriptase (RT) templates and poly(A) binding site (PBS) sequences, and scoring of predicted guide efficiency. Supports multiple PE variants (PE2, PE3, PE3b, PE-max) and custom PAM sequences.

## When to Use

**Use this skill when you:**
- Need to perform precise edits (SNVs, small insertions, small deletions) at a known genomic locus
- Want to leverage prime editing for edits where traditional CRISPR would require DSBs
- Need to design pegRNAs with optimized PBS and RT template sequences
- Want efficiency scoring and multi-guide candidates ranked by predicted effectiveness
- Are designing for human or mouse genomes with known genomic coordinates

**Do NOT use this skill when you:**
- Targeting large deletions (>100bp) — use standard CRISPR-Cas9 instead
- Working with organisms without available Ensembl genome annotations (hg38, hg19, mm10)
- Need wet-lab synthesis and validation — this generates computational designs only
- Require functional validation; all predictions are computational

## Expected Inputs

- **Gene symbol** (e.g., "TP53", "BRCA1") or explicit target sequence
- **Edit specification**: edit type (SNV/insertion/deletion), desired nucleotide change, position
- **Genome build**: hg38, hg19, or mm10
- **PE architecture**: PE2 (pegRNA alone), PE3 (pegRNA + nicking sgRNA)
- **Optional**: custom PAM sequence, PBS/RT length constraints

## Expected Outputs

- **TSV guide file**: spacer sequence, PBS, RT template, efficiency scores, nick sgRNA (if PE3)
- **Visualizations**: efficiency distribution, PBS Tm distribution, genomic map of guides
- **Summary**: top N designs ranked by combined score

## Procedure

1. **Fetch genomic sequence**
   - Query Ensembl REST API using gene symbol or explicit coordinates
   - Extract reference sequence spanning edit site ± 500bp

2. **Identify candidate protospacers**
   - Scan both DNA strands for NGG (or custom) PAM sequences
   - Filter: protospacers must be ≤30bp from edit site
   - Exclude protospacers with poly-T stretches (TTTT)

3. **Design PBS and RT template**
   - PBS: reverse-complement of pegRNA 3' end (8-17nt)
   - RT template: encodes desired edit + flanking context (10-40nt)
   - Calculate PBS Tm using nearest-neighbor thermodynamics

4. **Score guide efficiency**
   - GC content: optimal 40-70%
   - PBS Tm: optimal 30-37°C
   - Seed region (17-20bp): no off-target-prone patterns
   - Combined score: normalized average of subscores

5. **Design nicking sgRNAs (PE3 only)**
   - Find NGG PAMs on opposite strand
   - Distance from pegRNA: 40-100bp (optimal ~70bp)
   - Report top nicking guides alongside pegRNA

6. **Generate output**
   - Rank by efficiency score
   - Export TSV with all guide properties
   - Plot distributions and genomic context

## Key Execution Patterns

### Basic SNV design
```bash
python scripts/design_prime_editor_sgrnas.py \
  --gene TP53 \
  --species human \
  --edit-type snv \
  --position 7577121 \
  --wt-seq TCGCTAGGCGGATGCT \
  --edit-seq TCGCAAGGTGGATGCT \
  --editor PE3 \
  --top-n 20 \
  --outdir results_tp53_snv/
```

### Small insertion with custom PAM
```bash
python scripts/design_prime_editor_sgrnas.py \
  --gene BRCA1 \
  --edit-type insertion \
  --target-sequence AGGATCAAAGAGCGAGCAGAGTC \
  --editor PE2 \
  --pam NG \
  --pbs-length-range 10,15 \
  --outdir results_brca1_insert/
```

### Deletion using wild-type sequences
```bash
python scripts/design_prime_editor_sgrnas.py \
  --gene PTEN \
  --species mouse \
  --edit-type deletion \
  --wt-seq CACTGGTTATAAAGAGCGTAGCTAA \
  --edit-seq CACTGTAA \
  --genome mm10 \
  --editor PE3b \
  --nick-distance-range 50,80 \
  --outdir results_pten_del/
```

## Parameter Decision Guide

| Parameter | Guidance | Default | Notes |
|-----------|----------|---------|-------|
| `--editor` | PE2: pegRNA only (low WT off-target); PE3: add nicking sgRNA (higher efficiency) | PE3 | PE3b includes RT domain optimization |
| `--pam` | NGG for SpCas9 (most common); NNG for xCas9; NNGRRT for VQR | NGG | Custom PAM for alternative nuclease variants |
| `--pbs-length-range` | Longer PBS (14-17) = higher specificity; shorter (8-11) = higher activity | 8,17 | 10-13 often best compromise |
| `--rt-length-range` | Longer RT (25-40) encodes more context; shorter (10-15) faster folding | 10,40 | 15-25 commonly used in publications |
| `--nick-distance-range` | PE3 spacing; 60-80bp often optimal; <40bp or >100bp = reduced activity | 40,100 | Empirical from Anzalone et al. 2019 |
| `--spacer-len` | SpCas9/PE standard is 20bp; xCas9 can use 17-20bp | 20 | Fixed at 20 for PE compatibility |
| `--position` | Genomic coordinate (1-indexed); required if not providing explicit sequence | — | Must match selected genome build |

## Failure Modes

| Error | Cause | Solution |
|-------|-------|----------|
| "Gene not found" | Gene symbol not in Ensembl | Verify gene name; use explicit --target-sequence instead |
| "No candidate guides found" | All protospacers filtered out | Relax --gc-min/--gc-max; extend search radius from edit site |
| "HTTP 400 from Ensembl API" | Invalid genome build or coordinates | Check --genome flag; verify --position matches build (hg38 vs hg19 differ) |
| PBS Tm out of range | Few guides with optimal Tm | Accept 25-40°C as acceptable range or reduce --pbs-length-range |
| All poly-T guides filtered | Region has high AT content | Use --target-sequence to manually specify; consider alternative edits |

## Implementation Details

- **Sequence fetching**: Ensembl REST API (no authentication required)
- **PAM scanning**: IUPAC regex for flexibility (e.g., "NG" = [ACGT]G)
- **Tm calculation**: O(N) nearest-neighbor lookup table (no external thermodynamics library)
- **Scoring**: Weighted combination of GC%, Tm, seed quality, position heuristics
- **Plotting**: matplotlib with 300 DPI PNG output

## References

- Anzalone, A. V., et al. (2019). Search-and-replace genome editing. Nature 576, 149-157.
- Chen, P. J., et al. (2021). Enhanced prime editing systems by manipulating pegRNA-polymerase interactions. Nature Biotechnology 39, 1288-1293.

