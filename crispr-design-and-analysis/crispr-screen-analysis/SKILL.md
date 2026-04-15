---
name: crispr-screen-analysis
description: Analyze CRISPR pooled screens using MAGeCK. Accepts a raw count table or FASTQ files. Runs MAGeCK RRA (pairwise test) or MAGeCK MLE (multi-condition). Produces gene-level hit rankings, volcano plots, sgRNA rank plots, and hit summary figures for both positive and negative selection. Supports genome-wide KO screens, CRISPRi, and CRISPRa experiments.
---

# CRISPR Screen Analysis (MAGeCK)

## Purpose

Run a complete CRISPR pooled screen analysis pipeline using [MAGeCK](https://sourceforge.net/p/mageck/wiki/Home/). Starting from either a raw count matrix or FASTQ files, identify positively and negatively selected genes, and generate publication-quality figures.

## Use when

- the user has a CRISPR screen count table and wants to identify hits
- the user wants to find essential genes (negative selection) or enriched genes (positive selection)
- the user wants to compare treatment vs control conditions in a pooled screen
- the user has multiple conditions or time points (use MLE mode)
- the user has raw FASTQ files from a pooled screen and wants to go from reads to hits
- the user mentions MAGeCK, RRA, beta score, sgRNA count, or pooled screen analysis
- the user wants a volcano plot, rank plot, or gene score plot from screen data

## Do not use when

- the user wants arrayed screen analysis (not pooled)
- the user wants single-cell CRISPR screen analysis (Perturb-seq / CROP-seq)
- the user wants to design sgRNAs (use design-sgrnas-by-gene instead)
- the user has only a gene hit list without raw count data

## Supported modes

| Mode | Input | When to use |
|---|---|---|
| `test` | count table TSV | Pairwise comparison (treatment vs control); standard RRA |
| `mle` | count table TSV + design matrix | Multiple conditions, time-course, or complex designs |
| `count` | FASTQ files + sgRNA library | Convert raw reads to count table |
| `all` | FASTQ files + library | Full pipeline: count → test |

## Inputs

- `--count-table` — TSV with columns: sgRNA, gene, sample1, sample2, ... (required for `test`/`mle`)
- `--treatment` — comma-separated sample column names for the treatment group
- `--control` — comma-separated sample column names for the control group
- `--mode` — `test` (default), `mle`, `count`, or `all`
- `--fastq` — space-separated FASTQ paths (for `count`/`all` mode)
- `--library` — sgRNA library file: TSV with columns sgRNA_id, sequence, gene (for `count` mode)
- `--design-matrix` — tab-separated design matrix file (for `mle` mode)
- `--norm-method` — `median` (default), `total`, `control`, or `none`
- `--fdr-cutoff` — FDR significance threshold (default 0.25)
- `--lfc-cutoff` — |log2FC| threshold for volcano plot labeling (default 1.0)
- `--gene-lfc-method` — `median` or `alphamedian` (default `median`)
- `--top-n` — number of top hits to label in plots (default 20)
- `--prefix` — output file prefix (default `screen`)
- `--outdir` — output directory (default `./crispr_screen_results/`)

## Outputs

| File | Description |
|---|---|
| `{prefix}.gene_summary.txt` | MAGeCK gene-level results (RRA scores, p-values, FDR) |
| `{prefix}.sgrna_summary.txt` | sgRNA-level results (LFC, p-value per guide) |
| `{prefix}_volcano.png` | Gene-level volcano plot (neg/pos selection) |
| `{prefix}_rank_neg.png` | Negative selection rank plot (essential genes) |
| `{prefix}_rank_pos.png` | Positive selection rank plot |
| `{prefix}_hit_summary.png` | Top hit genes bar chart |
| `{prefix}_sgrna_scatter.png` | sgRNA-level LFC scatter (treatment vs control) |
| `{prefix}_qc.png` | Read count distribution and Gini index QC |
| `{prefix}_hits.tsv` | Filtered hit gene table (FDR < cutoff) |
| `{prefix}_count_table.txt` | (count mode) Normalized count table |

## Execution

```bash
# Pairwise test (most common)
python scripts/crispr_screen_analysis.py \
  --count-table counts.txt \
  --treatment Day14_rep1,Day14_rep2 \
  --control Day0_rep1,Day0_rep2 \
  --mode test \
  --outdir results/

# MLE with design matrix
python scripts/crispr_screen_analysis.py \
  --count-table counts.txt \
  --mode mle \
  --design-matrix design.txt \
  --outdir results/

# From FASTQ files
python scripts/crispr_screen_analysis.py \
  --fastq sample1.fastq.gz sample2.fastq.gz \
  --library sgrna_library.txt \
  --treatment sample1 \
  --control sample2 \
  --mode all \
  --outdir results/
```

## Parameter decision guide

| Signal in user request | Parameter to set |
|---|---|
| "pairwise", "treatment vs control", "two groups", "KO screen" | `--mode test` |
| "multiple conditions", "time course", "complex design", "beta score" | `--mode mle` |
| "from FASTQ", "raw reads", "sequencing data" | `--mode count` or `--mode all` |
| "essential genes", "negative selection", "dropout" | `--mode test`; focus on neg.score column |
| "enriched genes", "positive selection", "gain of function", "CRISPRa" | `--mode test`; focus on pos.score column |
| "more stringent", "high confidence hits only" | `--fdr-cutoff 0.05` |
| "permissive", "screening", "exploratory" | `--fdr-cutoff 0.25` (default) |
| "median normalization" or no mention | `--norm-method median` (default) |
| "spike-in controls", "non-targeting controls for normalization" | `--norm-method control` |
| "label top N hits" | `--top-n N` |
| replicate columns named Day0_rep1, Day0_rep2 | `--control Day0_rep1,Day0_rep2` |
| replicate columns named treated_A, treated_B | `--treatment treated_A,treated_B` |
| "strong effect size only" | `--lfc-cutoff 1.5` or `--lfc-cutoff 2.0` |

## MAGeCK installation note

MAGeCK must be installed separately. The recommended method is conda:

```bash
conda install -c bioconda mageck
```

Or via pip (Python-only implementation):

```bash
pip install mageck
```

This skill will automatically detect whether MAGeCK is available and fall back to a Python re-implementation of the core RRA algorithm if it is not installed.
