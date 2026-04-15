---
name: design_sgrnas_by_gene
description: Use when the user wants to retrieve CRISPR sgRNA candidates for one gene or a gene list from the GuidePro genome-wide reference table for human, monkey, or mouse, and return the top guides per gene.
---

## Purpose
Retrieve CRISPR sgRNA candidates for a gene or gene list from the GuidePro genome-wide reference table, filter them by the requested genes, and return a ranked output table suitable for downstream review or ordering.

## Use when
- the user asks to design sgRNAs for one gene or multiple genes
- the user provides a gene name, gene list, or a file containing genes
- the user wants top sgRNAs for CRISPR knockout design
- the user wants a lightweight and stable retrieval workflow based on GuidePro reference tables
- the user wants human, monkey, or mouse sgRNAs supported by GuidePro

## Do not use when
- the user wants genome-wide off-target analysis
- the user wants base editor or prime editor guide design
- the user wants guides designed from genomic coordinates rather than gene identifiers
- the user wants species other than human, monkey, or mouse
- the user wants a live website automation workflow rather than reference-table retrieval
- the user wants custom scoring beyond what is available in the GuidePro reference table

## Supported inputs
The user may provide one of the following:
- a single gene
- a comma-separated gene list
- a text file with one gene per line

Supported genome values:
- human
- monkey
- mouse

Common optional preferences:
- number of guides per gene
- whether to remove duplicated sgRNAs
- whether to refresh the cached GuidePro reference file
- whether to write a summary table

## Expected outputs
- a TSV table containing sgRNA candidates for the requested genes
- optional summary table
- cached local copy of the selected GuidePro genome-wide reference file
- source annotation indicating GuidePro

## Interaction model
Users may express their needs in natural language. Translate supported requests into script arguments when possible.

Examples of supported requests:
- "design sgRNAs for TP53"
- "design sgRNAs for TP53, EGFR, and EP300"
- "use the mouse reference"
- "keep top 5 guides per gene"
- "remove duplicates"
- "refresh the cache"
- "save a summary table"

If the request is supported, convert it into the appropriate command-line arguments.
If the request is unsupported, use the closest valid default and clearly state the limitation.

## Procedure
1. Parse the user input and determine whether the request contains a single gene, multiple genes, or a gene list file.
2. Determine the requested genome. Use the user-provided genome if specified; otherwise use the default genome.
3. Check whether the corresponding GuidePro genome-wide reference file is already present in the local cache.
4. If the reference file is missing or the user requests a refresh, download the correct GuidePro reference file.
5. Load the reference table.
6. Detect the gene column automatically unless the user provides a specific column name.
7. Detect the sgRNA column and score column when available.
8. Filter the reference table for the requested genes.
9. Remove duplicate sgRNAs if requested.
10. Sort guides by score when a score column is available.
11. Keep the top N guides per gene when requested.
12. Write the result table and optional summary table.
13. Return the output path and a short summary of what was done.

## Execution
Run the main script in this skill directory:

`python scripts/design_sgrnas_by_gene.py --output <output_tsv>`

Use the following argument patterns as needed:

### Single gene
`python scripts/design_sgrnas_by_gene.py --gene TP53 --genome human --output results/tp53.tsv`

### Multiple genes
`python scripts/design_sgrnas_by_gene.py --gene TP53,EGFR,EP300 --genome human --top-n-per-gene 10 --output results/multi_gene.tsv`

### Gene list file
`python scripts/design_sgrnas_by_gene.py --gene-file examples/example_genes.txt --genome mouse --top-n-per-gene 5 --output results/example_mouse.tsv`

### Common optional arguments
- `--gene-file <file>`
- `--genome <human|monkey|mouse>`
- `--cache-dir <dir>`
- `--refresh-cache`
- `--gene-col <col>`
- `--score-col <col>`
- `--sgrna-col <col>`
- `--top-n-per-gene <int>`
- `--drop-duplicates`
- `--summary <summary_tsv>`

## Parameter decision guide

Choose parameters based on what the user needs:

| Signal in user request | Parameter to set |
|---|---|
| Single gene (e.g. "design guides for TP53") | `--gene TP53` |
| Multiple genes (e.g. "for TP53, EGFR, MYC") | `--gene TP53,EGFR,MYC` |
| Gene list from a file | `--gene-file genes.txt` |
| "mouse" or "Trp53" (mouse gene symbol) | `--genome mouse` |
| "monkey" or primate study | `--genome monkey` |
| "top 5 guides" | `--top-n-per-gene 5` |
| "top 3 guides" (for a quick screen) | `--top-n-per-gene 3` |
| "all guides" (no limit) | omit `--top-n-per-gene` |
| "remove duplicates" or guides shared across genes | `--drop-duplicates` |
| "save a summary" or downstream reporting needed | `--summary results/summary.tsv` |
| "update / refresh" the reference | `--refresh-cache` |

## Defaults
Unless the user specifies otherwise:
- use `human` as the default genome
- use local cache if the reference file already exists
- keep top 10 guides per gene
- do not remove duplicates unless requested
- do not write a summary table unless requested

## Output conventions
- output format should be TSV
- preserve original GuidePro columns whenever possible
- include a source column or equivalent annotation when applicable
- preserve ranking-relevant columns such as GuidePro score if present
- keep the workflow deterministic and lightweight

## Dependencies
Install with:
`pip install -r requirements.txt`

## Failure modes
- no gene input is provided
- the requested genome is unsupported
- the GuidePro reference file cannot be downloaded
- the gene column cannot be detected automatically
- no requested genes are found in the reference table
- the score column is not detected, so ranking falls back to original row order
- the output directory is not writable

## Examples of when this skill should trigger
- "Design sgRNAs for TP53"
- "Get top 5 sgRNAs for EGFR and ERBB2 in human"
- "Use the mouse GuidePro reference to design guides for Trp53"
- "Retrieve sgRNAs for this gene list and save a summary"

## Examples of when this skill should not trigger
- "Check off-targets for these guides"
- "Design base editor guides for TP53"
- "Plot sgRNA scores as a heatmap"
- "Find CRISPRi guides around this TSS coordinate"