---
name: gsea_for_ranked_gene_list
description: Use when the user wants to run GSEA from a ranked gene list and generate publication-quality enrichment figures and result tables.
---

# GSEA Analysis for Ranked Gene List

Use this skill to run preranked GSEA from a ranked gene list and generate publication-quality plots and tabular outputs.

## Purpose

This skill performs Gene Set Enrichment Analysis using a ranked gene list. It is intended for users who already have a ranking metric for genes, such as log2 fold change, signal-to-noise score, Wald statistic, t statistic, correlation coefficient, or another signed ranking score.

This skill is appropriate for:
- preranked GSEA
- pathway-level enrichment from ranked genes
- publication-quality enrichment plots
- pathway result tables for downstream interpretation

## Use when

Use this skill when the user:
- provides a ranked gene list
- wants GSEA rather than over-representation analysis
- wants pathway-level enrichment from a continuous ranking
- wants KEGG, Reactome, Hallmark, or GO-based GSEA
- wants publication-quality enrichment plots and result tables

Typical requests:
- “Run GSEA on this ranked gene list.”
- “Use Hallmark pathways for this ranked list.”
- “Run Reactome GSEA and generate plots.”
- “Plot the top enriched pathways from my ranked genes.”
- “Perform preranked GSEA for this differential expression result.”

## Do not use when

Do not use this skill when the user:
- provides only an unranked gene list
- wants simple GO over-representation analysis instead of GSEA
- wants sgRNA design or CRISPR analysis
- wants gene-level volcano or heatmap plotting rather than pathway enrichment
- wants a custom enrichment method not supported by the script

If the user provides only a gene list without ranking values, use a gene-list enrichment skill instead.

## Expected inputs

The skill expects:
- a ranked gene list file with at least two columns:
  - gene
  - ranking score

Optional inputs:
- organism
- gene set library
- number of permutations
- minimum and maximum gene set size
- number of top pathways to plot
- output prefix
- figure style parameters

## Supported libraries

Supported libraries in the script:
- `HALLMARK`
- `KEGG`
- `REACTOME`
- `GO_BP`
- `GO_MF`
- `GO_CC`

## Expected outputs

For each run, generate:
- full GSEA result table
- filtered significant result table
- publication-quality enrichment plot(s) for top pathways
- optional summary table

## Assumptions

This skill assumes:
- the input is a ranked gene list
- higher positive values correspond to one biological direction and lower negative values to the opposite direction
- the ranking metric is already computed by the user
- the user wants preranked GSEA

## Workflow

1. Read the ranked gene list.
2. Validate that the file contains a gene column and a ranking-score column.
3. Clean gene symbols and numeric scores.
4. Sort the ranked list by score in descending order.
5. Resolve the requested gene set library.
6. Run preranked GSEA.
7. Save the full result table.
8. Filter significant pathways.
9. Generate publication-quality enrichment plots for the top positively and negatively enriched pathways.
10. Save outputs and return file paths.

## Execution

Run:

`python scripts/gsea_analysis_for_ranked_gene_list.py --input <ranked_gene_file> --output-prefix <prefix>`

Common optional arguments:
- `--gene-col <col>`
- `--score-col <col>`
- `--library <HALLMARK|KEGG|REACTOME|GO_BP|GO_MF|GO_CC>`
- `--organism <human|mouse>`
- `--permutation-num <int>`
- `--min-size <int>`
- `--max-size <int>`
- `--top-n <int>`
- `--fdr-cutoff <float>`
- `--fig-width <float>`
- `--fig-height <float>`
- `--dpi <int>`
- `--font-family <string>`
- `--font-size <float>`

## Interaction model

Users may express requests in natural language. Translate supported requests into script arguments when possible.

Examples:
- “Use Hallmark gene sets”
- “Run Reactome GSEA”
- “Plot top 5 enriched pathways”
- “Use human genes”
- “Use FDR < 0.1”
- “Make the figure publication-ready”
- “Increase font size”
- “Save all outputs with a custom prefix”

## Output conventions

- save tables as TSV
- save figures as PNG and PDF
- preserve NES, ES, nominal p-value, and FDR when available
- clearly distinguish positively and negatively enriched pathways
- keep output deterministic and readable

## Dependencies

Install with:
`pip install -r requirements.txt`

## Failure modes

Fail clearly when:
- the input file is empty
- the gene column is missing
- the score column is missing
- ranking values are non-numeric
- no valid ranked genes remain after cleaning
- the requested library is unsupported
- GSEA returns no valid pathways
- no output file can be written

## Trigger examples

This skill should trigger for prompts like:
- “Run GSEA on this ranked gene list.”
- “Use Hallmark pathways for this ranked list.”
- “Perform preranked GSEA and plot the top pathways.”
- “Run Reactome GSEA from this differential expression ranking.”

## Non-trigger examples

This skill should not trigger for prompts like:
- “Run GO enrichment on this gene list.”
- “Make a volcano plot.”
- “Design sgRNAs for TP53.”
- “Curate a gene list from literature.”
