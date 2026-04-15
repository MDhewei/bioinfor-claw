---
name: crispr-library-design
description: Design a pooled CRISPR sgRNA library for a list of target genes. Designs sgRNAs for each gene, selects top guides based on predicted efficiency and off-target scores, includes non-targeting controls, and outputs an oligo order file ready for library synthesis.
tags:
  - crispr
  - library-design
  - pooled-screen
  - guide-design
  - oligo-synthesis
version: 1.0.0
---

## Purpose

Automate design of genome-scale or focused pooled CRISPR sgRNA libraries. Given a list of target genes, this tool designs sgRNAs for each gene, applies efficiency and off-target scoring, selects top guides per gene, includes non-targeting controls, and generates an oligo order file formatted for DNA synthesis vendors (Twist Bioscience, CustomArray, etc.). Supports vector design with configurable flanking sequences for directional cloning.

## When to Use

**Use this skill when you:**
- Designing a pooled CRISPR screen library targeting specific genes
- Need sgRNA sequences in a standard format for synthesis
- Want multiple guides per gene to maximize coverage
- Designing for human or mouse genomes (hg38, hg19, mm10)
- Need to include non-targeting controls for phenotypic validation
- Working with SpCas9 or variants with different PAM sequences

**Do NOT use this skill when you:**
- Designing single guides for a few genes — use design-prime-editor-sgrnas instead
- Need wet-lab cloning/validation — this generates sequence designs only
- Target organisms without available Ensembl annotations
- Require off-target analysis against non-canonical targets (e.g., cryptic PAMs)
- Need to design guide RNAs for non-Cas9 nucleases (use specialized tools)

## Expected Inputs

- **Gene list**: comma-separated symbols or file (one gene per line)
- **Number of guides per gene**: typically 3-5 for saturation screens, 1-2 for focused screens
- **Number of non-targeting controls**: 500-1000 typical
- **Vector sequences**: 5' and 3' flanking for cloning
- **Filters**: GC content range, no poly-T restriction sites, etc.

## Expected Outputs

- **Library design TSV**: guide ID, gene, spacer, full oligo, position, strand, GC content, efficiency score, rank
- **FASTA file** (optional): spacer sequences only
- **GC distribution plot**: histogram of GC content across library
- **Efficiency distribution plot**: histogram of efficiency scores
- **Library statistics**: total guides, per-gene counts, genome coverage

## Procedure

1. **Parse gene list**
   - If comma-separated, split by comma; if file path, read one gene per line
   - Validate gene symbols against Ensembl

2. **Fetch exon sequences**
   - For each gene, query Ensembl `/lookup/symbol/{species}/{gene}`
   - Extract exon coordinates
   - Fetch coding sequence via `/sequence/region`
   - If --targeting-exon-only, restrict to exon regions

3. **Scan for PAM sites**
   - Search both forward and reverse strands
   - Find all NGG (or custom PAM) sites in exon regions
   - Extract 20bp spacer upstream of PAM

4. **Filter candidates**
   - GC content: user-specified range (default 30-80%)
   - Poly-T restriction: remove TTTT stretches if --no-t4-stretch
   - BsmBI sites: remove CACCG (common in Golden Gate cloning)
   - Sequence complexity: reject guides with excessive homopolymer runs

5. **Score each guide**
   - GC content: quadratic penalty if outside optimal 40-60%
   - Seed region (bp 17-20): penalize off-target-prone patterns
   - Position in exon: prefer 5' exon and middle positions
   - Nucleotide rules: simplified Doench 2016 model
     - Penalize GC-rich seeds
     - Favor specific nucleotides at positions 1-10, 17-20
   - Combined score: weighted average of subscores

6. **Select top N guides per gene**
   - Sort by efficiency score
   - Select top N (default 4) guides
   - Ensure diversity (ideally different exons/positions)

7. **Generate non-targeting controls**
   - Create random 20bp sequences (no PAM constraint)
   - Filter: same GC range as targeting guides
   - Check against genome: ensure no CACCG match
   - Validate: no TTTT, ≤4bp homopolymer runs
   - Generate ~500-1000 controls

8. **Add vector sequences**
   - 5' prefix: default ACCG (for BsmBI-mediated Golden Gate)
   - 3' suffix: default GTTTTAGAGCTA (for GFP expression)
   - Full oligo = prefix + spacer + suffix

9. **Output and summary**
   - TSV format: one row per guide
   - FASTA format (optional): spacer sequences
   - Summary statistics: total guides, per-gene breakdown
   - Visualizations: GC distribution, efficiency scores

## Key Execution Patterns

### Basic library design (4 guides per gene, 500 controls)
```bash
python scripts/crispr_library_design.py \
  --genes TP53,BRCA1,EGFR,MYC \
  --species human \
  --n-guides-per-gene 4 \
  --n-controls 500 \
  --outdir my_library/
```

### Genome-scale library from file
```bash
python scripts/crispr_library_design.py \
  --genes human_genes.txt \
  --species human \
  --n-guides-per-gene 4 \
  --n-controls 1000 \
  --genome hg38 \
  --outdir library_hg38/
```

### Focused library with strict GC and exon-only guides
```bash
python scripts/crispr_library_design.py \
  --genes DNA_repair_genes.txt \
  --n-guides-per-gene 3 \
  --targeting-exon-only \
  --no-t4-stretch \
  --gc-min 0.4 \
  --gc-max 0.6 \
  --outdir strict_library/
```

### Library with custom vector sequences for directional cloning
```bash
python scripts/crispr_library_design.py \
  --genes target_list.txt \
  --vector-prefix ACCGCACCG \
  --vector-suffix GTTTTAGAGCTACCTTGC \
  --library-name MyCustomLib \
  --output-format fasta \
  --outdir custom_lib/
```

### Mouse library with alternative PAM
```bash
python scripts/crispr_library_design.py \
  --genes mouse_genes.txt \
  --species mouse \
  --genome mm10 \
  --pam NG \
  --n-guides-per-gene 5 \
  --outdir mm10_library/
```

## Parameter Decision Guide

| Parameter | Guidance | Default | Notes |
|-----------|----------|---------|-------|
| `--n-guides-per-gene` | More guides = higher coverage, larger library, higher cost; 4 is standard | 4 | Use 2-3 for focused, 5-10 for saturation |
| `--n-controls` | Non-targeting controls for validation; 1-2% of total guides typical | 500 | Use 500-1000; critical for phenotypic validation |
| `--gc-min`, `--gc-max` | GC content range; 40-60% often optimal for specificity | 0.3, 0.8 | 30-80% is permissive; 40-60% is strict |
| `--pam` | NGG for SpCas9; NG for xCas9; custom for variants | NGG | Must match nuclease used in experiments |
| `--vector-prefix` | 5' flanking for cloning; ACCG (BsmBI recognition); adjust for your vector | ACCG | Verify with your cloning strategy |
| `--vector-suffix` | 3' flanking; GTTTTAGAGCTA standard for GFP polyA | GTTTTAGAGCTA | Adjust for your expression cassette |
| `--targeting-exon-only` | Restrict to exon regions (coding); off-target risk vs coverage trade-off | False | Use for focused screens; more specific |
| `--no-t4-stretch` | Remove TTTT sequences (T4 polymerase risk, secondary structure) | False | Recommended for mammalian screens |
| `--output-format` | csv/tsv for spreadsheet; fasta for sequence analysis | tsv | tsv most common for synthesis orders |

## Failure Modes

| Error | Cause | Solution |
|-------|-------|----------|
| "Gene not found in Ensembl" | Gene symbol not recognized | Check spelling; use official HGNC/MGI symbol; verify genome build |
| "No guides found for gene X" | PAM sites filtered out by GC/sequence rules | Relax --gc-min/--gc-max; disable --no-t4-stretch; use --targeting-exon-only=false |
| "Library has 0 non-targeting controls" | Random sequence generation failing | Check --gc-min/--gc-max are reasonable (0.3-0.8) |
| Mostly low efficiency scores | Scoring parameters misaligned | Verify GC range optimal for your Cas9; check --pam matches nuclease |
| Off-target warnings | Some guides may cut at secondary sites | Add off-target filtering step (separate tool); use xCas9/eSpCas9 if high specificity needed |
| FASTA output has wrong format | Sequence errors or vector issues | Check --vector-prefix/--vector-suffix; verify header format |

## Implementation Details

- **PAM scanning**: IUPAC regex pattern matching for flexibility
- **GC calculation**: (G_count + C_count) / total_length
- **Scoring**: Combination of GC%, position, seed region quality (simplified Doench 2016)
- **Ensembl API**: REST-based, no authentication; handles coordinates, sequences
- **Random control generation**: numpy.random with reproducible seed
- **Plotting**: matplotlib 300 DPI PNG

## Design Principles

1. **Multiple guides per gene**: provides redundancy, improves hit calling
2. **Non-targeting controls**: critical for assessing off-target effects and phenotype baseline
3. **Filtered PAM sites**: quality over quantity; fewer better-designed guides outperform many poor ones
4. **Vector design**: flank with restriction sites for efficient cloning (BsmBI, BbsI common)
5. **GC content**: 40-60% reduces off-target binding, improves secondary structure

## Quality Checklist

- [ ] All gene symbols validated in Ensembl
- [ ] GC content histogram shows expected distribution (40-80% most guides)
- [ ] Efficiency scores reasonable (0.5-1.0 for top guides)
- [ ] No guides with TTTT or excessive homopolymers
- [ ] Controls properly labeled and counted
- [ ] Vector flanking sequences verified against cloning protocol
- [ ] Total guide count reasonable (N_genes × N_guides_per_gene + N_controls)

## References

- Doench, J. G., et al. (2016). Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9. Nature Biotechnology 34, 184-191.
- Broad Institute sgRNA Designer: https://portals.broadinstitute.org/gppx/

