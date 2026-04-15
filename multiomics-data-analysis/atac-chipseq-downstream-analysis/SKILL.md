---
name: atac-chipseq-downstream-analysis
description: Downstream analysis of ATAC-seq or ChIP-seq peak files. Annotates peaks to genomic features, computes signal summaries, finds differential peaks between conditions, and generates QC plots.
---

# ATAC/ChIP-seq Downstream Analysis

## Purpose

Perform comprehensive downstream analysis on ATAC-seq or ChIP-seq peak files. This skill:
- **Annotates peaks** to genomic features (promoter, exonic, intronic, intergenic, downstream)
- **Computes QC metrics** (peak width distribution, score distribution, TSS enrichment)
- **Identifies differential peaks** between conditions (unique vs shared)
- **Generates visualizations** (peak annotation pie charts, heatmaps, Venn diagrams)

Suitable for both narrow peaks (punctate transcription factor binding) and broad peaks (histone marks).

## Reuse policy

This skill is designed for:
- Downstream analysis of peak files from MACS2, SEACR, or similar peak callers
- Annotation to genes and genomic features
- QC assessment of peak quality
- Comparative peak analysis across multiple conditions
- Generating publication-quality peak annotation figures

This skill requires **peak files** (BED/narrowPeak/broadPeak format) and optional **sample metadata**.

## Inputs

### Required
- `--peaks` — single peak file in BED/narrowPeak/broadPeak format
  - BED3: chrom, start, end
  - BED6: chrom, start, end, name, score, strand
  - narrowPeak (10 col): chrom, start, end, name, score, strand, signalValue, pvalue, qvalue, peak
  - broadPeak (9 col): chrom, start, end, name, score, strand, signalValue, pvalue, qvalue
  - File extension can be .bed, .narrowPeak, .broadPeak, or .gff

### Optional
- `--mode` — analysis mode
  Choices: `annotate`, `differential`, `qc`, `all`
  Default: `all`

- `--genome` — reference genome for annotation
  Choices: `hg38`, `hg19`, `mm10`, `mm39`
  Default: `hg38`

- `--peak-files` — comma-separated list of peak files for differential analysis
  Example: `peaks_condition1.bed,peaks_condition2.bed,peaks_condition3.bed`
  Required if `--mode` includes `differential`

- `--conditions` — comma-separated condition labels matching `--peak-files`
  Example: `untreated,treated_1h,treated_6h`
  Must match length of `--peak-files`

- `--gene-list` — optional comma-separated gene symbols to highlight in annotation
  Example: `TP53,EGFR,MYC`

- `--outdir` — output directory
  Default: `./atac_chipseq_output`

- `--tss-window` — base pairs around TSS to count peaks (for TSS enrichment)
  Default: `2000`

- `--min-score` — minimum peak score/signalValue to keep (filters low-confidence peaks)
  Default: `0` (no filtering)

- `--top-n` — analyze top N peaks (sorted by score descending)
  Default: `50000` (no limit)

## Inputs file formats

### Peak file (--peaks)

narrowPeak format (ENCODE standard):
```
chr1  1000   2000  peak1  100  .  1.5  10.2  5.1  500
chr1  5000   6500  peak2  200  .  2.0  15.5  7.2  450
chr2  3000   4200  peak3  150  .  1.8  12.1  6.0  600
```

BED3 format (minimal):
```
chr1  1000  2000
chr1  5000  6500
chr2  3000  4200
```

### Condition labels (--conditions)

Simple comma-separated list:
```
--conditions "ATAC_untreated,ATAC_dexamethasone_4h,ATAC_dexamethasone_24h"
```

## Outputs

### Annotation mode (--mode annotate or all)
- `peaks_annotated.tsv` — annotated peak table with columns:
  - `peak_id`, `chrom`, `start`, `end`, `score`, `annotation`, `nearest_gene`, `distance_to_tss`
  - annotation: one of {promoter, exonic, intronic, intergenic, downstream}

- `peak_annotation_summary.txt` — human-readable summary of annotation counts

- `peak_annotation_distribution.png` — pie chart of peak class distribution

### QC mode (--mode qc or all)
- `qc_metrics.txt` — summary metrics:
  - Total peaks, median peak width, mean/median score
  - TSS enrichment score, fraction of promoter peaks

- `peak_width_distribution.png` — histogram of peak widths (log scale)
- `peak_score_distribution.png` — histogram of peak signal values

### Differential mode (--mode differential or all)
- `condition_comparison.txt` — summary of overlap statistics
- `unique_peaks_<CONDITION>.bed` — condition-specific unique peaks (not overlapping >50% with others)
- `shared_peaks_all_conditions.bed` — peaks present in all conditions
- `peak_overlap_venn.png` — Venn diagram of peak overlaps (2-3 conditions)
- `peak_overlap_matrix.png` — pairwise Jaccard similarity heatmap (4+ conditions)

### General
- `analysis_summary.txt` — overall analysis summary with command used
- All PNG plots at 300 DPI

## Execution policy

Construct the command from required and optional parameters.

### Command template

```
python scripts/atac_chipseq_downstream.py \
  --peaks <PEAK_FILE> \
  [--mode <MODE>] \
  [--genome <GENOME>] \
  [--peak-files <FILE1,FILE2,...>] \
  [--conditions <COND1,COND2,...>] \
  [--gene-list <GENE1,GENE2,...>] \
  [--outdir <OUTDIR>] \
  [--tss-window <INT>] \
  [--min-score <FLOAT>] \
  [--top-n <INT>]
```

### Example commands

#### Basic annotation of single peak file
```
python scripts/atac_chipseq_downstream.py \
  --peaks atac_peaks.narrowPeak \
  --genome hg38 \
  --outdir results/
```

#### QC analysis with filtering
```
python scripts/atac_chipseq_downstream.py \
  --peaks peaks.bed \
  --mode qc \
  --min-score 10 \
  --top-n 10000 \
  --tss-window 3000 \
  --outdir results/
```

#### Differential analysis across three conditions
```
python scripts/atac_chipseq_downstream.py \
  --peaks cond1.narrowPeak \
  --mode differential \
  --genome mm10 \
  --peak-files "cond1.narrowPeak,cond2.narrowPeak,cond3.narrowPeak" \
  --conditions "untreated,treated_4h,treated_24h" \
  --gene-list "NR3C1,FKBP5" \
  --outdir results/
```

#### Full analysis with gene highlighting
```
python scripts/atac_chipseq_downstream.py \
  --peaks merged_peaks.narrowPeak \
  --mode all \
  --genome hg38 \
  --gene-list "TP53,MDM2,CDKN1A" \
  --min-score 5 \
  --outdir results/
```

## Parameter decision guide

| User intent | Parameter to set |
|---|---|
| "Annotate my ATAC peaks" | `--mode annotate` (default) |
| "Compare peak sets between conditions" | `--mode differential --peak-files ... --conditions ...` |
| "Just check peak quality" | `--mode qc` |
| "Do everything" | `--mode all` (default) |
| Broad histone mark peaks | Use as-is; script auto-detects format |
| Narrow TF binding peaks | Use as-is; script auto-detects format |
| Filter noisy low-signal peaks | `--min-score 10` (or appropriate value) |
| Focus on high-confidence peaks | `--top-n 5000` or `--top-n 10000` |
| Stricter promoter definition | `--tss-window 1000` |
| Lenient promoter definition | `--tss-window 5000` |
| Different genome | `--genome mm10` or `--genome hg19` |
| Highlight specific genes in results | `--gene-list "GENE1,GENE2,GENE3"` |

## Failure conditions

Fail clearly if:
- `--peaks` is missing or file does not exist
- Peak file cannot be parsed (unknown format or corrupted)
- `--mode differential` requested but `--peak-files` or `--conditions` missing
- `--peak-files` and `--conditions` have mismatched lengths
- `--genome` is not in {hg38, hg19, mm10, mm39} (used for gene lookup)
- Ensembl REST API is unreachable (user should check internet connection)
- All peaks are filtered out by `--min-score` or `--top-n`
- Less than 2 peak files provided for differential analysis

## Agent trigger examples

**Trigger this skill when the user asks:**
- "Annotate my ATAC peaks to genes"
- "What genes are near my ChIP-seq peaks?"
- "Compare peak sets between treated and untreated"
- "Generate QC plots for my peak file"
- "Find differential peaks between conditions"
- "Make a Venn diagram of my peak overlaps"
- "Assess peak quality and distribution"

**Do NOT trigger this skill when the user asks:**
- "Call peaks from my BAM file" → use a peak calling skill (e.g., MACS2)
- "Align my FASTQ files" → use a read alignment skill
- "Do motif enrichment on my peaks" → use a motif analysis skill
- "RNA-seq analysis" → use rnaseq-differential-expression
- "Convert peak format" → use file conversion tools

## Notes

- **Gene annotation** uses Ensembl REST API (https://rest.ensembl.org/) to fetch gene coordinates by genomic region. Requires internet connection.
- **Peak overlap** (for differential analysis) uses ≥50% reciprocal overlap as the criterion for "same peak".
- **TSS enrichment** = fraction of input peaks with center within `--tss-window` bp of any known gene TSS.
- **Promoter definition**: within 2000 bp (default) of gene TSS.
- All plots saved at 300 DPI in PNG format (PDF optional in future versions).
- If Ensembl API fails, annotation is skipped with a warning; basic peak QC always completes.
