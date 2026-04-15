---
name: methylation-analysis
description: Analyze DNA methylation data from Illumina 450K/EPIC arrays or bisulfite sequencing. Computes differential methylation between groups, identifies differentially methylated regions (DMRs), annotates CpGs, and visualizes results.
---

# DNA Methylation Analysis

## Purpose

Comprehensive analysis of DNA methylation data from **Illumina 450K, EPIC arrays, or whole-genome bisulfite sequencing (WGBS)**. This skill:
- **Computes QC metrics** (beta value distributions, sample clustering)
- **Identifies differentially methylated positions (DMPs)** using t-tests with FDR correction
- **Identifies differentially methylated regions (DMRs)** by clustering adjacent DMPs
- **Annotates methylation** to genes and genomic features
- **Generates visualizations** (volcano plots, Manhattan plots, DMR heatmaps)

Suitable for comparing methylation across treatment groups, disease states, or developmental stages.

## Reuse policy

This skill is designed for:
- Detection of DMPs and DMRs between two or more sample groups
- Quality control and clustering analysis of methylation datasets
- Annotation of CpGs to genes and regulatory features
- Publication-quality volcano plots and methylation heatmaps
- Integration with downstream enrichment analyses

This skill requires **beta value matrix** (CpG × sample) and **sample metadata** with group assignments.

## Inputs

### Required
- `--input` — beta value matrix (TSV file)
  - Rows = CpG IDs (e.g., cg00000029)
  - Columns = sample names
  - Values = beta values (0 to 1, representing methylation level)
  - Can also accept M-values (log2-transformed)

- `--metadata` — sample metadata table (TSV file)
  - Must contain columns: `sample_id`, and at least one grouping column (default: `group`)
  - `sample_id` must match column names in `--input`

- `--ref-group` — reference/control group label
  - Delta beta computed as: treatment_mean - reference_mean
  - Example: `normal`, `control`, `WT`

### Optional
- `--mode` — analysis mode
  Choices: `dmp`, `dmr`, `qc`, `all`
  Default: `all`

- `--group-col` — column name in metadata for grouping
  Default: `group`

- `--fdr-cutoff` — FDR-adjusted p-value threshold for significance
  Default: `0.05`

- `--delta-beta` — minimum absolute delta beta for significance
  Default: `0.2` (common in methylation studies)

- `--outdir` — output directory
  Default: `./methylation_output`

- `--cpg-annotation` — optional TSV file with CpG annotations
  - Columns: `CpG_id`, `gene`, `feature_type` (promoter/gene_body/intergenic)

- `--array-type` — array platform for context (affects interpretation)
  Choices: `450K`, `EPIC`, `WGBS`
  Default: `450K`

## Input file formats

### Beta value matrix (--input)

```
CpG_id       sample1  sample2  sample3  sample4
cg00000029   0.234    0.256    0.845    0.901
cg00000108   0.512    0.498    0.234    0.201
cg00000109   0.145    0.167    0.812    0.798
```

### Metadata (--metadata)

```
sample_id    group
sample1      normal
sample2      normal
sample3      tumor
sample4      tumor
```

### CpG annotation (--cpg-annotation, optional)

```
CpG_id       gene         feature_type
cg00000029   TP53         promoter
cg00000108   EGFR         gene_body
cg00000109   MYC          intergenic
```

## Outputs

### QC mode (--mode qc or all)
- `qc_summary.txt` — per-sample QC metrics:
  - Beta value mean, median, SD per sample
  - Detection p-value summary (if available)

- `beta_distribution.png` — density histograms of beta values per sample
- `sample_clustering_heatmap.png` — hierarchical clustering heatmap of samples based on CpG methylation

### DMP mode (--mode dmp or all)
- `dmp_results.tsv` — all tested CpGs with columns:
  - `CpG_id`, `chr`, `position`, `mean_ref`, `mean_treat`, `delta_beta`, `pvalue`, `fdr`
  - Filtered to include significant DMPs if requested

- `significant_dmps.tsv` — subset passing `--fdr-cutoff` AND `--delta-beta` threshold

- `dmp_volcano_plot.png` — volcano plot (delta_beta vs -log10 FDR)
- `dmp_manhattan_plot.png` — Manhattan-style plot (if chromosome/position info available)

### DMR mode (--mode dmr or all)
- `dmr_results.tsv` — identified DMRs with columns:
  - `chrom`, `start`, `end`, `n_cpgs`, `mean_delta_beta`, `mean_fdr`

- `dmr_heatmap.png` — per-CpG beta value heatmap for top 20 DMRs (samples as columns, CpGs as rows)

### Annotation (if --cpg-annotation provided)
- `dmp_annotated.tsv` — DMPs annotated with gene and feature type
- `feature_distribution.png` — bar chart of DMPs by feature (promoter/gene_body/intergenic)

### General
- `analysis_summary.txt` — overall analysis summary
- All PNG plots at 300 DPI

## Execution policy

Construct the command from required and optional parameters.

### Command template

```
python scripts/methylation_analysis.py \
  --input <BETA_MATRIX> \
  --metadata <METADATA_FILE> \
  --ref-group <REF_GROUP> \
  [--mode <MODE>] \
  [--group-col <COLUMN>] \
  [--fdr-cutoff <FLOAT>] \
  [--delta-beta <FLOAT>] \
  [--cpg-annotation <ANNOTATION_FILE>] \
  [--array-type <TYPE>] \
  --outdir <OUTDIR>
```

### Example commands

#### Basic DMP analysis
```
python scripts/methylation_analysis.py \
  --input beta_values.tsv \
  --metadata metadata.tsv \
  --ref-group control \
  --outdir results/
```

#### QC only
```
python scripts/methylation_analysis.py \
  --input beta_values.tsv \
  --metadata metadata.tsv \
  --ref-group control \
  --mode qc \
  --outdir results/
```

#### Stricter thresholds with annotation
```
python scripts/methylation_analysis.py \
  --input beta_values.tsv \
  --metadata metadata.tsv \
  --ref-group normal \
  --cpg-annotation cpg_anno.tsv \
  --mode all \
  --fdr-cutoff 0.01 \
  --delta-beta 0.3 \
  --outdir results/
```

#### EPIC array with DMR focus
```
python scripts/methylation_analysis.py \
  --input beta_epic.tsv \
  --metadata meta.tsv \
  --ref-group WT \
  --array-type EPIC \
  --mode dmr \
  --fdr-cutoff 0.05 \
  --outdir results/
```

## Parameter decision guide

| User intent | Parameter to set |
|---|---|
| "Just QC my samples" | `--mode qc` |
| "Find differential methylation" | `--mode dmp` (default with --mode all) |
| "Find DMRs, not just individual CpGs" | `--mode dmr` |
| "Do full analysis" | `--mode all` (default) |
| Stricter significance | `--fdr-cutoff 0.01 --delta-beta 0.3` |
| Lenient/discovery screen | `--fdr-cutoff 0.1 --delta-beta 0.1` |
| Different group column | `--group-col disease_status` |
| EPIC array (vs 450K) | `--array-type EPIC` |
| Custom CpG annotations | `--cpg-annotation genes.tsv` |
| Two-group comparison | Use `--ref-group control` for control-vs-treatment |

## Failure conditions

Fail clearly if:
- `--input`, `--metadata`, or `--ref-group` are missing
- Beta value matrix cannot be parsed or is empty
- Metadata cannot be parsed or is missing required columns
- `--ref-group` label is not found in metadata group column
- Sample IDs in metadata do not match column names in beta matrix
- All CpGs filtered out by significance thresholds
- Less than 2 samples per group (t-test invalid)
- Beta values are outside [0, 1] range (without M-value conversion)

## Agent trigger examples

**Trigger this skill when the user asks:**
- "Find differentially methylated positions in my data"
- "DMR analysis on methylation array data"
- "What CpGs are differently methylated between tumor and normal?"
- "Annotate methylation changes to genes"
- "QC check on my EPIC methylation dataset"
- "Identify DMRs from 450K data"

**Do NOT trigger this skill when the user asks:**
- "Align bisulfite-seq reads" → use read alignment skill
- "Call methylation from BAM" → use methylation calling skill (e.g., Bismark)
- "RNA-seq analysis" → use rnaseq-differential-expression
- "ChIP-seq peak analysis" → use atac-chipseq-downstream-analysis

## Notes

- **Statistical testing**: Uses Welch's t-test (not assuming equal variance) with Benjamini-Hochberg FDR correction.
- **Delta beta direction**: positive = higher methylation in treatment; negative = lower methylation in treatment.
- **DMR identification**: Adjacent significant DMPs (within 1000 bp) are clustered into regions.
- **Beta values** range from 0 (unmethylated) to 1 (fully methylated). Script auto-detects if values are M-values (outside [0,1]).
- **Recommended thresholds**: FDR < 0.05, delta beta > 0.2 (for publication); adjust as needed for discovery.
- All plots saved at 300 DPI in PNG format.
- CpG annotation is optional; without it, only positional information is provided.
