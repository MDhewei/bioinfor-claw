---
name: normal-tissue-expression-for-gene
description: Retrieve and summarize normal tissue expression for a single human gene using GTEx, including full tissue-level expression, publication-quality barplot output, and a simple expression-pattern classification such as universally expressed, tissue-specific, mixed, or non-expressed.
---

# Normal Tissue Expression for Gene

## Purpose

Use this skill to retrieve and summarize **normal tissue expression for a single gene**.

This skill is intended for common daily bioinformatics questions such as:
- In which normal tissues is this gene highly expressed?
- Is this gene broadly expressed across tissues or tissue-specific?
- Is this gene lowly expressed or effectively non-expressed in GTEx?
- What does the expression distribution across all normal tissues look like?

This skill is **single-gene focused** by design. It is optimized for fast, interpretable inspection of one gene at a time rather than large matrix-style gene-list analysis.

## Scope

This skill covers:
- **normal tissue expression**
- **GTEx-based tissue-level summaries**
- **publication-quality visualization**
- **basic expression-pattern classification**

This skill does **not** cover:
- cancer cell line expression
- tumor expression
- protein abundance
- gene essentiality
- pathway enrichment
- multi-gene heatmaps as the primary use case

Those should be handled by separate skills.

## Data source

Primary data source:
- **GTEx normal tissue expression**

The implementation should prefer **official GTEx API access** when available.

Recommended retrieval flow:
1. Resolve the queried gene against GTEx gene reference data
2. Retrieve median tissue expression values for that gene
3. Normalize the returned fields into a consistent table
4. Generate summary outputs and figures

## When to use

Use this skill when the user asks things like:
- “Show me normal tissue expression of TP53.”
- “Which tissues highly express EGFR?”
- “Is GATA3 tissue-specific?”
- “How broadly is MYC expressed across tissues?”
- “Plot GTEx tissue expression for this gene.”

## When not to use

Do not use this skill when the user asks for:
- cancer cell line expression
- tumor expression
- gene-list heatmaps as the main deliverable
- protein expression
- CRISPR dependency
- gene annotation unrelated to tissue expression

Use another skill for those requests.

## Input

## Required input

### `--gene`
A single gene symbol.

Examples:
- `TP53`
- `EGFR`
- `ERBB2`
- `GATA3`

## Optional input

### `--top-n`
Number of top tissues to summarize in the text summary.

Default:
- `10`

### `--expr-detect-threshold`
Expression threshold used to define whether a tissue is considered “detected” or expressed.

Default:
- `1.0`

This threshold is also shown on the barplot.

### `--outdir`
Output directory for all generated files.

---

## Parameter decision guide

| Signal in user request | Parameter to set |
|---|---|
| General tissue expression question | use defaults (`--top-n 10 --expr-detect-threshold 1.0`) |
| "top 5 tissues" or quick summary | `--top-n 5` |
| "broadly expressed" / "housekeeping gene" | use default threshold; classify result as `universally-expressed` |
| "tissue-specific" question | check output `expression_category`; lower threshold (`--expr-detect-threshold 0.5`) if gene has moderate tissue-specific expression |
| "lowly expressed" gene (e.g. receptor, non-abundant TF) | lower threshold: `--expr-detect-threshold 0.5` to catch subtle tissue bias |
| Gene expected to be very highly expressed | threshold can be raised (`--expr-detect-threshold 5.0`) to focus on strong expression |

## Why this skill is gene-level instead of gene-list-level

This skill is intentionally centered on **one gene at a time**, because the most common practical use case is:

1. the user has a gene of interest,
2. they want to quickly inspect its tissue distribution,
3. they want an interpretable classification and figure.

For this use case, a single-gene workflow is usually more useful than a matrix-first workflow.

A multi-gene normal-tissue-expression skill can be built separately if needed.

## Outputs

A successful run should produce:

- one tissue-level expression table
- one concise summary text file
- one summary TSV file
- one publication-quality barplot in PNG
- one publication-quality barplot in PDF

### Typical output files

- `<GENE>.gtex_tissues.tsv`
- `summary.txt`
- `gene_expression_summary.tsv`
- `<GENE>.gtex_all_tissues.png`
- `<GENE>.gtex_all_tissues.pdf`

## Output table structure

### Tissue-level table
The per-gene tissue table should contain:

- `gene_symbol`
- `gencode_id`
- `tissue`
- `expression`
- `rank_within_gene`

### Summary table
The summary TSV should contain:

- `query_gene`
- `resolved_gene_symbol`
- `gencode_id`
- `expression_category`
- `max_expression`
- `median_expression_across_tissues`
- `num_tissues`
- `num_detected_tissues`
- `fraction_detected_tissues`

## Core workflow

### Step 1. Parse and normalize the input gene
- Read the gene from `--gene`
- Normalize formatting
- Treat the query as a human gene symbol

### Step 2. Resolve the gene in GTEx
- Query GTEx reference gene information
- Resolve:
  - canonical gene symbol
  - GENCODE ID

If the gene cannot be resolved, fail clearly.

### Step 3. Retrieve GTEx normal tissue expression
- Query GTEx median tissue expression
- Normalize the response into a consistent table with:
  - tissue
  - expression

### Step 4. Rank and summarize tissues
- Sort tissues by expression
- Keep all tissues in the output table
- Use `top-n` only for the short text summary

### Step 5. Classify the expression pattern
Assign one of the following categories:
- `universally-expressed`
- `tissue-specific`
- `mixed`
- `non-expressed`
- `no-data`

### Step 6. Generate a publication-quality barplot
- Show **all tissues by default**
- Use a horizontal barplot
- Draw a threshold line for the expression detection threshold
- Save as both PNG and PDF

### Step 7. Write output files
Write the tissue table, summary files, and figure to `outdir`.

## Expression classification logic

This skill should provide a **simple, interpretable classification** rather than an overly complicated tissue-specificity score.

Recommended logic:

### `non-expressed`
Use when:
- the maximum tissue expression is below the detection threshold

### `tissue-specific`
Use when:
- maximum expression is clearly high
- expression in the top tissue is much higher than the median across tissues
- only a relatively small fraction of tissues are above the detection threshold

### `universally-expressed`
Use when:
- most tissues are above the detection threshold
- median expression across tissues is not low

### `mixed`
Use when:
- the gene is expressed in some tissues
- but does not strongly fit either universal or highly tissue-specific behavior

### `no-data`
Use when:
- no GTEx result is returned for the gene

## Threshold guidance

The default expression detection threshold is:

- **1.0**

This is a practical and commonly used threshold for judging whether expression is present at a meaningful level.

Interpretation:
- below threshold: low or not detected
- above threshold: detectable expression

The threshold should also be displayed visually on the barplot.

## Plotting requirements

The barplot should be **publication quality**.

### Plot style requirements
- horizontal barplot
- all tissues shown by default
- tissues sorted by expression
- high resolution
- readable tissue labels
- clean axes
- minimal unnecessary styling
- threshold line clearly visible

### Because the plot is horizontal
Expression is on the x-axis, so the threshold should be plotted as a **vertical line**, not a horizontal line.

### Recommended outputs
- PNG for quick viewing
- PDF for publication or slides

## Execution

## Entry point

```bash
python scripts/normal_tissue_expression_for_gene.py --gene "<GENE>" --outdir <output_dir>