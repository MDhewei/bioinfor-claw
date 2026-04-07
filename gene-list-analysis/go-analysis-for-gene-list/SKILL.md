---
name: go-analysis-for-gene-list
description: Use when the user wants to run enrichment analysis from a gene list and generate publication-quality bubble plots for Gene Ontology, KEGG, and Reactome terms.
---

# GO / KEGG / Reactome Enrichment Plotter

Use this skill to perform over-representation analysis from a user-provided gene list and generate publication-quality enrichment figures and result tables.

## Purpose

This skill is designed for routine functional enrichment analysis of a gene list. It retrieves enriched terms from supported Enrichr libraries and produces clean, publication-ready bubble plots together with tabular outputs.

Use this skill when the user wants:
- Gene Ontology enrichment
- KEGG pathway enrichment
- Reactome pathway enrichment
- bubble plots for enriched terms
- result tables suitable for downstream interpretation or figure preparation

This skill is optimized for:
- human and mouse gene lists
- unranked gene lists
- over-representation analysis
- publication-ready static figures

## Use when

Use this skill when the user:
- provides a gene list directly or through a file
- asks for GO enrichment analysis
- asks for KEGG or Reactome pathway enrichment
- asks for BP, MF, or CC enrichment
- asks for a bubble plot, dot plot, or enrichment figure
- wants top enriched terms visualized
- wants a result table and figure saved to disk

Typical requests:
- “Run GO enrichment on this gene list.”
- “Make a publication-quality bubble plot for BP terms.”
- “Do KEGG and Reactome enrichment for these genes.”
- “Show the top 15 GO terms with adjusted p-value.”
- “Generate enrichment plots for human genes.”

## Do not use when

Do not use this skill when the user:
- provides a ranked gene list and wants GSEA rather than over-representation analysis
- wants gene set enrichment from fold change values or a preranked vector
- wants genome-wide annotation downloads without enrichment analysis
- wants custom pathway databases that are not supported by the script
- wants interactive web-based plots rather than static publication figures
- wants organism support outside the currently supported set
- wants network visualization, pathway maps, or graph-based pathway diagrams

If the user wants ranked-list enrichment or GSEA, use a dedicated GSEA workflow instead of this skill.

## Supported organisms

Supported organism inputs are normalized internally.

Accepted user inputs include:
- `human`
- `homo sapiens`
- `hs`
- `hsapiens`
- `mouse`
- `mus musculus`
- `mm`

These map to:
- `human`
- `mouse`

If the user gives an unsupported organism, fail clearly and explain the allowed values.

## Supported libraries

The script supports the following categories:

- `GO_BP` for Gene Ontology Biological Process
- `GO_MF` for Gene Ontology Molecular Function
- `GO_CC` for Gene Ontology Cellular Component
- `GO` for all three GO categories above
- `KEGG` for KEGG pathway enrichment
- `REACTOME` for Reactome pathway enrichment
- `ALL` for GO_BP, GO_MF, GO_CC, KEGG, and REACTOME together

## Expected inputs

The skill expects:
- a plain-text gene list file with one gene per line, or equivalent clean gene list input
- an output prefix for saving results

Optional inputs include:
- organism
- library category
- number of top terms to plot
- significance cutoff
- figure width and height
- DPI
- font family
- font size
- colormap

## Expected outputs

For each requested category, generate:

1. a full enrichment result table  
   Example:
   - `<output_prefix>.GO_BP.tsv`
   - `<output_prefix>.KEGG.tsv`

2. a reduced plotting table containing the plotted terms  
   Example:
   - `<output_prefix>.GO_BP.plot_input.tsv`

3. a publication-quality bubble plot in PNG format  
   Example:
   - `<output_prefix>.GO_BP.bubble.png`

4. a publication-quality bubble plot in PDF format  
   Example:
   - `<output_prefix>.GO_BP.bubble.pdf`

## Assumptions

This skill assumes:
- the input is an unranked gene list
- gene symbols are already reasonably clean
- the user wants over-representation analysis
- enrichment is performed against Enrichr-supported libraries through `gseapy`
- adjusted p-value is the main significance metric for ranking and plotting

If the user does not specify otherwise, use:
- `human` as the default organism
- `ALL` as the default library set only if the user explicitly asks for multiple categories; otherwise use the requested category
- `top_n = 15`
- `sig_cutoff = 0.05`
- publication-style figure defaults from the script

## Workflow

1. Read the input gene list.
2. Clean and normalize gene symbols by trimming whitespace and converting to uppercase.
3. Normalize the organism input to a supported value.
4. Resolve the requested enrichment library or library panel.
5. Run enrichment analysis using `gseapy.enrichr`.
6. Validate the returned result table and expected columns.
7. Compute plotting fields such as:
   - hit count
   - background count
   - gene ratio
   - `-log10(FDR)`
8. Rank terms by adjusted p-value and enrichment strength.
9. Select significant terms for plotting. If no term passes the significance cutoff, fall back to the top ranked terms.
10. Generate a publication-quality bubble plot.
11. Save full results, plotting input, PNG figure, and PDF figure.
12. Return the saved output paths and summarize what was generated.

## Plot behavior

The default plot is a bubble plot where:
- x-axis = gene ratio
- y-axis = enriched term
- bubble size = hit count
- bubble color = `-log10(FDR)`

The script is designed for publication-style output:
- wrapped term labels for readability
- clean axes
- no top/right spines
- bubble size legend
- color bar for significance
- PNG and PDF export

## Execution

Run the main script:

`python scripts/go_enrichment_plotter.py --input <gene_file> --output-prefix <prefix>`

### Common examples

Run all supported libraries for a human gene list:

`python scripts/go_enrichment_plotter.py --input examples/example_genes.txt --output-prefix results/example_all --organism human --library ALL --top-n 15`

Run only GO Biological Process:

`python scripts/go_enrichment_plotter.py --input examples/example_genes.txt --output-prefix results/example_bp --organism human --library GO_BP --top-n 15`

Run only KEGG:

`python scripts/go_enrichment_plotter.py --input examples/example_genes.txt --output-prefix results/example_kegg --organism human --library KEGG --top-n 15`

Run only Reactome:

`python scripts/go_enrichment_plotter.py --input examples/example_genes.txt --output-prefix results/example_reactome --organism human --library REACTOME --top-n 15`

### Supported arguments

Required:
- `--input <gene_file>`
- `--output-prefix <prefix>`

Optional:
- `--organism <human|mouse|hs|hsapiens|homo sapiens|mm|mus musculus|mouse>`
- `--library <GO_BP|GO_MF|GO_CC|GO|KEGG|REACTOME|ALL>`
- `--top-n <int>`
- `--sig-cutoff <float>`
- `--fig-width <float>`
- `--fig-height <float>`
- `--dpi <int>`
- `--font-family <string>`
- `--font-size <float>`
- `--cmap <matplotlib_colormap>`

## Parameter interpretation

Translate natural-language user requests into script arguments where possible.

Examples:

- “Use human genes”  
  → `--organism human`

- “Run BP only”  
  → `--library GO_BP`

- “Run KEGG and Reactome”  
  → run the script twice or use separate category handling if needed by the caller

- “Show top 10 terms”  
  → `--top-n 10`

- “Use adjusted p < 0.01”  
  → `--sig-cutoff 0.01`

- “Use Arial font and 12 pt labels”  
  → `--font-family Arial --font-size 12`

- “Make the figure higher resolution”  
  → `--dpi 300` or higher if explicitly requested

## Output conventions

- Save tables as TSV.
- Save figures as both PNG and PDF.
- Preserve original enrichment fields from the enrichment result whenever possible.
- Use deterministic sorting:
  - adjusted p-value ascending
  - hit count descending
  - gene ratio descending
  - term name ascending
- If no term passes the significance threshold, still generate a plot from the top ranked terms and note that fallback behavior.

## Dependencies

Install with:

`pip install -r requirements.txt`

Expected packages include:
- `pandas`
- `numpy`
- `matplotlib`
- `gseapy`

## Failure modes

Fail clearly when:
- the input gene list is empty after cleaning
- the organism is unsupported
- the requested library is unsupported
- enrichment returns no results
- the result table is missing expected columns
- no output file can be written
- no enrichment category succeeds when multiple categories were requested

If one category fails but others succeed, continue with the successful categories and report which category failed.

## Robot usage notes

When another robot uses this skill:
- prefer explicit user-specified library and organism values
- do not assume ranked-list enrichment; this skill is for unranked lists
- do not silently switch to GSEA
- do not invent unsupported organisms or library names
- return the output file paths in the response
- mention whether the plot used significant terms only or fell back to top ranked terms

## Trigger examples

This skill should trigger for prompts like:
- “Run GO enrichment on this gene list.”
- “Make a bubble plot for BP enrichment.”
- “Do KEGG pathway enrichment for these genes.”
- “Run Reactome analysis and generate publication-ready figures.”
- “Perform GO, KEGG, and Reactome enrichment for this list.”

## Non-trigger examples

This skill should not trigger for prompts like:
- “Run GSEA on my ranked gene list.”
- “Plot a volcano plot for these DEGs.”
- “Find sgRNAs for these genes.”
- “Curate a literature-based gene list for ferroptosis.”
- “Download all GO annotations for human genes.”