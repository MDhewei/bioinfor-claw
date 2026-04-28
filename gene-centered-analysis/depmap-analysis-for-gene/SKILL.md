---
name: depmap-analysis-for-gene
description: "DepMap analysis for a gene across CANCER CELL LINES (NOT patient samples). Modules: expression, mutation, copy number, essentiality. Data is streamed directly from the DepMap API — no full dataset download needed. Also includes standalone co-expression (depmap_coexpression.py) and co-essentiality (depmap_coessentiality.py) scripts for cell-line correlations. For patient/tumor co-expression use coexpression-for-gene (TCGA/GTEx) instead."
---

# DepMap Analysis for Gene

## Purpose

Analyze a **single gene** in DepMap across cancer cell lines.

The main script streams data directly from the DepMap download API and extracts ONLY the target gene's column. It never downloads or loads multi-GB whole-genome matrices into memory.

Supported modules (main script):
- expression
- mutation
- copy_number
- essentiality

Co-expression and co-essentiality require full matrices and use **standalone scripts** (see below).

## Inputs

### Required
- `--gene`: Gene symbol, e.g. TP53
- `--outdir`: Output directory

### Module selection
- `--modules`: Comma-separated list. Choices: expression, mutation, copy_number, essentiality, full. Default: full.

Examples:
- `--modules expression`
- `--modules mutation,copy_number`
- `--modules full`

### Analysis options
- `--top-n`: Top rows to report (default: 20)
- `--amp-threshold`: Copy number amplification threshold (default: 1.0)
- `--del-threshold`: Copy number deletion threshold (default: -1.0)
- `--essential-threshold`: Essentiality cutoff (default: -0.5)

### IMPORTANT: No data file arguments needed

The script streams directly from DepMap URLs. Do NOT pass `--expression-file`, `--mutation-file`, `--copy-number-file`, `--essentiality-file`, or `--metadata-file`. These legacy flags are accepted but ignored.

Do NOT call `depmap-download-data` or `depmap_downloader.py` before running this script. The script handles its own data fetching.

## Command template

```bash
python scripts/depmap_analysis_for_gene.py --gene <GENE> --modules <MODULES> --outdir <OUTDIR> [--top-n N] [--amp-threshold X] [--del-threshold X] [--essential-threshold X]
```

### Example commands

Expression only:
```bash
python scripts/depmap_analysis_for_gene.py --gene TP53 --modules expression --outdir <OUTDIR>
```

Mutation + copy number:
```bash
python scripts/depmap_analysis_for_gene.py --gene BRCA1 --modules mutation,copy_number --outdir <OUTDIR>
```

Full analysis:
```bash
python scripts/depmap_analysis_for_gene.py --gene KRAS --modules full --outdir <OUTDIR>
```

## Parameter decision guide

| Signal in user request | Parameter to set |
|---|---|
| "expression across cell lines" | `--modules expression` |
| "mutation landscape" | `--modules mutation` |
| "copy number / amplification / deletion" | `--modules copy_number` |
| "essential? dependency? CRISPR screen?" | `--modules essentiality` |
| "what genes are co-expressed with X?" | **MUST use standalone script** `depmap_coexpression.py` |
| "what genes are co-essential with X?" | **MUST use standalone script** `depmap_coessentiality.py` |
| "co-expression in TCGA / patient samples" | **Wrong skill!** Use `coexpression-for-gene` |
| "full profile" or no specific module | `--modules full` |

## Outputs

### Always
- `gene_depmap_summary.tsv`
- `summary.txt`

### Expression
- `<GENE>.expression_top_cell_lines.tsv`
- `<GENE>.expression_lineage_summary.tsv`
- `<GENE>.expression_barplot.png / .pdf`

### Mutation
- `<GENE>.mutations.tsv`
- `<GENE>.mutation_type_summary.tsv`
- `<GENE>.mutation_protein_summary.tsv`

### Copy number
- `<GENE>.copy_number.tsv`
- `<GENE>.copy_number_top_amplified.tsv`
- `<GENE>.copy_number_top_deleted.tsv`
- `<GENE>.copy_number_lineage_summary.tsv`
- `<GENE>.copy_number_amplified_barplot.png / .pdf`

### Essentiality
- `<GENE>.essentiality.tsv`
- `<GENE>.essentiality_top_cell_lines.tsv`
- `<GENE>.essentiality_lineage_summary.tsv`
- `<GENE>.essentiality_barplot.png / .pdf`

## Standalone co-expression / co-essentiality scripts

These require the full expression/essentiality matrix and produce richer outputs (barplot, network, FDR table).

### DepMap co-expression
```bash
python scripts/depmap_coexpression.py \
    --gene TP53 \
    --method pearson \
    --top-n 30 \
    --fdr-cutoff 0.05
```

### DepMap co-essentiality
```bash
python scripts/depmap_coessentiality.py \
    --gene TP53 \
    --method pearson \
    --top-n 30 \
    --fdr-cutoff 0.05
```

For co-expression/co-essentiality: do NOT include `--expression-file` or `--essentiality-file` in EXEC_ARGS. The server auto-locates cached DepMap data and injects the correct paths.

## Failure conditions

Fail clearly if:
- `--gene` is missing
- `--outdir` is missing
- the gene cannot be found in a DepMap dataset
- network access to DepMap API fails
