---
name: tcga-expression-for-gene
description: Query and visualise TCGA expression data for a single gene across cancer cohorts. Supports pan-cancer bar plots, single-cohort boxplots, and tumor-vs-normal comparisons. All data is fetched live from the GDC API — no local files required.
---

# TCGA Expression for Gene

## Purpose

Retrieve and visualise **gene expression across TCGA cancer cohorts** for a single gene of interest.

Three analysis modes are available:

| Mode | What it does |
|------|-------------|
| `pan_cancer` (default) | Ranks all TCGA cohorts by median expression of the gene; produces a bar chart |
| `single_cohort` | Plots the expression distribution for one cancer type (e.g. BRCA, LUAD) |
| `tumor_vs_normal` | Side-by-side boxplot of tumour vs matched normal tissue within one cohort |

All expression data uses **UQ-FPKM** units (upper-quartile normalised FPKM) from the GDC RNA-Seq pipeline. An alternative `median_centered_log2_uqfpkm` unit is also available.

## Reuse policy

Use this skill when the user wants to:
- Know which cancer type has the **highest / lowest expression** of a gene
- Compare **tumour vs. normal** expression in a specific cohort
- Get a quick **pan-cancer expression landscape** for a gene of interest
- Find cohorts for downstream survival or mutation analysis

This skill does **not** require pre-downloaded data. Data is fetched live from the GDC REST API.

## Script location

```
gene-centered-analysis/tcga_expression_for_gene/scripts/tcga_expression_for_gene.py
```

## Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--gene` | *(required)* | Gene symbol, e.g. `TP53`, `EGFR`, `MYC` |
| `--mode` | `pan_cancer` | One of `pan_cancer`, `single_cohort`, `tumor_vs_normal` |
| `--cancer-type` | — | TCGA cohort code for single-cohort / TVN modes, e.g. `BRCA`, `LUAD`, `GBM` |
| `--top-n` | `20` | Number of top cohorts shown in pan-cancer bar chart |
| `--tsv-units` | `uqfpkm` | Expression unit: `uqfpkm` or `median_centered_log2_uqfpkm` |
| `--outdir` | *(required)* | Output directory for plots and TSV tables |

## Outputs

All outputs are written to `--outdir`:

- `{gene}_pan_cancer_expression.tsv` / `{gene}_pan_cancer_barplot.png` — Pan-cancer summary table and bar chart
- `{gene}_{cohort}_expression.tsv` / `{gene}_{cohort}_boxplot.png` — Single-cohort distribution
- `{gene}_{cohort}_tumor_vs_normal.tsv` / `{gene}_{cohort}_tvn_boxplot.png` — Tumour vs. normal comparison

## Example usage

```bash
# Pan-cancer expression landscape for TP53
python scripts/tcga_expression_for_gene.py \
  --gene TP53 \
  --mode pan_cancer \
  --top-n 20 \
  --outdir results/tp53_pan_cancer

# EGFR expression in lung adenocarcinoma (tumour only)
python scripts/tcga_expression_for_gene.py \
  --gene EGFR \
  --mode single_cohort \
  --cancer-type LUAD \
  --outdir results/egfr_luad

# MYC tumour vs. normal in breast cancer
python scripts/tcga_expression_for_gene.py \
  --gene MYC \
  --mode tumor_vs_normal \
  --cancer-type BRCA \
  --outdir results/myc_brca_tvn
```

## Parameter decision guide

**Which mode should I use?**
- "Which cancer has the most/least expression?" → `pan_cancer`
- "What is expression like in a specific cancer?" → `single_cohort`
- "Is the gene upregulated in tumour vs. normal?" → `tumor_vs_normal`

**`--cancer-type` is required for `single_cohort` and `tumor_vs_normal` modes.**
Common TCGA cohort codes: `BRCA`, `LUAD`, `LUSC`, `GBM`, `LGG`, `COAD`, `READ`, `PRAD`, `KIRC`, `UCEC`, `SKCM`, `LIHC`, `STAD`, `BLCA`, `HNSC`, `OV`, `PAAD`, `THCA`, `CESC`, `SARC`.

**`--tsv-units`:**
- `uqfpkm` — raw UQ-FPKM (better for between-cohort comparisons)
- `median_centered_log2_uqfpkm` — log2 centred (better for visualising relative expression)

## Dependencies

```
pandas
matplotlib
requests
```

Install: `pip install pandas matplotlib requests`
