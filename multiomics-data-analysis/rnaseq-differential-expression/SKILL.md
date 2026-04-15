---
name: rnaseq-differential-expression
description: Perform RNA-seq differential expression analysis between two sample groups from a count matrix. Supports DESeq2-style (via pydeseq2), Welch's t-test, and Mann-Whitney U methods. Outputs full DE results, significant gene table, volcano plot, MA plot, and top-gene heatmap. Outputs feed directly into go-analysis-for-gene-list and gsea-for-ranked-gene-list.
---

# RNA-seq Differential Expression Analysis

## Purpose

Identify differentially expressed genes between two sample groups from a **gene count or expression matrix**.

Supported statistical methods:
- `deseq2` — DESeq2-style negative binomial model via `pydeseq2` (recommended for raw counts)
- `ttest` — Welch's t-test on log2-normalized CPM (fast; good for normalized data)
- `mannwhitney` — Mann-Whitney U test on log2-normalized CPM (non-parametric)

Default method: `deseq2`

## Reuse policy

This skill is designed for:
- Bulk RNA-seq two-group comparison (treated vs control, tumor vs normal, etc.)
- Generating a DE gene list to feed into `go-analysis-for-gene-list`
- Generating a ranked gene list to feed into `gsea-for-ranked-gene-list`
- Producing publication-quality volcano plots, MA plots, and heatmaps

This skill requires **user-provided local files** (count matrix + metadata).

## Inputs

### Required
- `--counts` — count/expression matrix (TSV or CSV)
  - Rows = genes (gene symbols or IDs)
  - Columns = sample names
  - Values = raw counts (for `deseq2`) or normalized expression (for `ttest`, `mannwhitney`)

- `--metadata` — sample metadata table (TSV or CSV)
  - Must contain columns: `sample`, `group`
  - `sample` must match column names in `--counts`
  - `group` contains the group labels used in `--group-a` and `--group-b`

- `--group-a` — reference/control group label (e.g. `control`, `normal`)
- `--group-b` — condition/treatment group label (e.g. `treated`, `tumor`)
- `--outdir` — output directory

### Optional

- `--method`
  Choices: `deseq2`, `ttest`, `mannwhitney`
  Default: `deseq2`

- `--min-count`
  Default: `10`
  Genes must have at least this count in at least `--min-samples` samples.

- `--min-samples`
  Default: `2`
  Minimum number of samples a gene must pass `--min-count` in.

- `--fc-thresh`
  Default: `1.0`
  |log2FC| threshold for calling significance.

- `--padj-thresh`
  Default: `0.05`
  Adjusted p-value (BH FDR) threshold.

- `--prefix`
  Default: `<group_b>_vs_<group_a>`
  Output file prefix.

- `--top-n-heatmap`
  Default: `50`
  Number of top significant genes to show in the heatmap.

## Input file formats

### Count matrix (--counts)

```
gene       sample1  sample2  sample3  sample4
TP53       1200     980      450      420
EGFR       320      290      1800     2100
MYC        5000     4800     900      870
```

### Metadata (--metadata)

```
sample    group
sample1   control
sample2   control
sample3   treated
sample4   treated
```

## Outputs

### Always
- `<prefix>.de_results.tsv` — full DE results for all tested genes

  Columns: `gene`, `log2FC`, `mean_control`, `mean_condition`, `pvalue`, `padj`
  (+ `base_mean`, `lfcSE`, `stat` when using `deseq2`)

- `<prefix>.significant_genes.tsv` — subset passing `--fc-thresh` and `--padj-thresh`

- `de_summary.tsv` — one-row analysis summary

- `summary.txt` — human-readable summary

### Plots
- `<prefix>.volcano.png` / `.pdf` — volcano plot (log2FC vs -log10 padj)
- `<prefix>.ma_plot.png` / `.pdf` — MA plot (mean expression vs log2FC)
- `<prefix>.top_genes_heatmap.png` / `.pdf` — Z-score heatmap of top significant genes

## Execution policy

Construct the command from the required arguments and optional parameters.

### Command template

```
python scripts/rnaseq_differential_expression.py \
  --counts <COUNTS_FILE> \
  --metadata <METADATA_FILE> \
  --group-a <CONTROL_LABEL> \
  --group-b <CONDITION_LABEL> \
  [--method <METHOD>] \
  [--min-count <N>] \
  [--min-samples <N>] \
  [--fc-thresh <FC>] \
  [--padj-thresh <PADJ>] \
  [--prefix <PREFIX>] \
  [--top-n-heatmap <N>] \
  --outdir <OUTDIR>
```

### Example commands

#### DESeq2 (default) — raw counts
```
python scripts/rnaseq_differential_expression.py \
  --counts counts.tsv \
  --metadata metadata.tsv \
  --group-a control \
  --group-b treated \
  --outdir results/
```

#### t-test on normalized data
```
python scripts/rnaseq_differential_expression.py \
  --counts normalized_expr.tsv \
  --metadata metadata.tsv \
  --group-a normal \
  --group-b tumor \
  --method ttest \
  --fc-thresh 1.5 \
  --padj-thresh 0.01 \
  --outdir results/
```

#### Mann-Whitney with custom prefix
```
python scripts/rnaseq_differential_expression.py \
  --counts counts.csv \
  --metadata meta.csv \
  --group-a WT \
  --group-b KO \
  --method mannwhitney \
  --prefix KO_vs_WT_mw \
  --outdir results/
```

## Downstream skill chaining

The outputs of this skill feed directly into other bioinfor-claw skills:

| This skill's output                     | Next skill                      | How                                  |
|-----------------------------------------|---------------------------------|--------------------------------------|
| `<prefix>.significant_genes.tsv`        | `go-analysis-for-gene-list`     | Pass significant gene list           |
| `<prefix>.de_results.tsv`               | `gsea-for-ranked-gene-list`     | Rank by log2FC or -log10(padj)×sign  |
| `<prefix>.significant_genes.tsv`        | `annotation-for-gene-list`      | Annotate significant genes           |
| Gene name from results                  | `depmap-analysis-for-gene`      | Validate top hits in DepMap          |

## Parameter decision guide

| Signal in user request | Parameter to set |
|---|---|
| Raw integer counts (from featureCounts, HTSeq, STAR) | `--method deseq2` (default) |
| Already normalized (TPM, RPKM, log2) | `--method ttest` |
| Small n (< 5 samples per group) | `--method mannwhitney` (non-parametric, no distributional assumption) |
| "fast / exploratory" DE | `--method ttest` |
| Stricter significance (e.g. drug targets) | `--fc-thresh 1.5 --padj-thresh 0.01` |
| Lenient / discovery screen | `--fc-thresh 0.5 --padj-thresh 0.1` |
| Noisy / lowly expressed data | raise `--min-count 20 --min-samples 3` |
| "heatmap with more genes" | `--top-n-heatmap 100` |
| "heatmap with top 20" | `--top-n-heatmap 20` |
| Custom output prefix | `--prefix tumor_vs_normal` |



| Data type                   | Recommended method  |
|-----------------------------|---------------------|
| Raw integer counts          | `deseq2`            |
| Log2-normalized / TPM / RPKM | `ttest`            |
| Small n (< 5 per group)     | `mannwhitney`       |
| Fast exploratory analysis   | `ttest`             |

## Failure conditions

Fail clearly if:
- `--counts`, `--metadata`, `--group-a`, `--group-b`, or `--outdir` are missing
- `--group-a` or `--group-b` labels are not found in the metadata `group` column
- Samples in metadata are not found as columns in the count matrix
- All genes are filtered out by `--min-count`/`--min-samples`
- The count matrix is empty or not parseable

## Agent trigger examples

**Trigger this skill when the user asks:**
- "Run differential expression analysis on my RNA-seq data"
- "Which genes are upregulated in treated vs control samples?"
- "Make a volcano plot from my count matrix"
- "Perform DESeq2 analysis comparing tumor and normal"
- "I have a counts TSV and a metadata file, run DE for me"
- "Find significantly changed genes between KO and WT"

**Do NOT trigger this skill when the user asks:**
- "Align my FASTQ files" → use a read alignment skill (not yet available)
- "Run GO enrichment on a gene list" → use `go-analysis-for-gene-list`
- "Perform GSEA" → use `gsea-for-ranked-gene-list`
- "Show TCGA expression for a gene" → use `tcga-expression-for-gene`
- "Analyse scRNA-seq data" → single-cell requires a different skill

## Notes

- `deseq2` requires `pydeseq2` (`pip install pydeseq2`) and **integer raw counts**.
  If `pydeseq2` is not installed, the method automatically falls back to `ttest`.
- For `ttest` and `mannwhitney`, the script applies log2(CPM + 1) normalization
  internally. Do not pre-normalize if using these methods on raw counts.
- Adjusted p-values use Benjamini-Hochberg FDR correction.
- log2FC direction: **positive = higher in group-b**, negative = lower in group-b.
- The heatmap shows Z-scores (per-gene standardized across samples), not raw counts.
