---
name: annotation_for_gene_list
description: Use when the user wants comprehensive annotation for one gene or a gene list, including gene identifiers, protein information, concise functional summary, domain architecture, Gene Ontology terms, pathway annotations, disease associations, and source-backed metadata.
---

# Annotation for Gene List

Use this skill to build a structured annotation table for one gene or a gene list. It is designed for users who want a concise but information-rich summary of gene identity, protein features, biological role, domain composition, pathway context, and disease relevance.

## Purpose

This skill aggregates per-gene annotations into a clean, analysis-ready table. It is intended to help users quickly understand what a gene is, what protein it encodes, what domains it contains, what broad biological functions it performs, what pathways it is associated with, and what disease connections are known.

This skill is especially useful when the user wants:
- identifier mapping
- gene and protein names
- concise gene-function summaries
- protein domain annotation
- GO annotations
- pathway annotations
- disease associations
- a structured annotation sheet for downstream analysis or figure preparation

## Use when

Use this skill when the user:
- provides one gene or a gene list
- asks to annotate genes
- wants gene identifiers mapped
- wants protein names or UniProt IDs
- wants domain information
- wants a concise function summary
- wants GO, pathway, or disease annotations
- wants a consolidated annotation table
- wants a TSV or CSV output for further analysis

Typical requests:
- “Annotate this gene list.”
- “Give me full gene and protein annotation for these genes.”
- “Add domain information for these genes.”
- “Summarize what these genes do.”
- “Map these symbols to Ensembl, Entrez, and UniProt IDs.”
- “Give pathway and disease annotations for these genes.”

## Do not use when

Do not use this skill when the user:
- wants over-representation enrichment analysis rather than per-gene annotation
- wants preranked GSEA
- wants CRISPR sgRNA design
- wants literature curation as the primary task
- wants expression profiling or essentiality analysis as the primary task
- wants TCGA-specific survival, mutation, or expression analysis
- wants detailed wet-lab protocol information

If the user specifically wants expression, dependency, DepMap, CCLE, or TCGA analysis, use a dedicated downstream skill for those tasks instead.

## Expected inputs

The skill expects one of the following:
- a single gene
- a comma-separated gene list
- a plain-text file with one gene per line

Supported optional inputs:
- organism
- output format
- a request for concise annotation only
- a request for specific fields such as domain information or disease annotations

## Supported organisms

Supported user inputs include:
- `human`
- `homo sapiens`
- `hs`
- `hsapiens`
- `mouse`
- `mus musculus`
- `mm`

These are normalized internally.

Default organism:
- `human`

## Expected outputs

The skill generates a structured annotation table with one row per gene and one column per annotation field.

Output formats supported by the script:
- TSV
- CSV

The default behavior is to produce a TSV file.

## Default annotation schema

Use these columns unless the user requests otherwise:

- `input_gene`
- `mapping_status`
- `mapping_note`
- `approved_symbol`
- `full_gene_name`
- `aliases`
- `organism`
- `entrez_id`
- `ensembl_gene_id`
- `uniprot_accession`
- `protein_name`
- `protein_length`
- `gene_type`
- `chromosome`
- `cytoband`
- `strand`
- `genomic_location`
- `reported_function`
- `subcellular_location`
- `domain_annotation`
- `domain_coordinates`
- `key_domain_summary`
- `functional_domain_class`
- `molecular_function`
- `biological_process`
- `cellular_component`
- `reactome_pathways`
- `kegg_pathways`
- `hallmark_pathways`
- `pathway_annotation_summary`
- `disease_association`
- `protein_class`
- `druggability`
- `source_summary`

## Field intent

Use the fields with the following interpretation:

### Identifier and mapping fields
- `input_gene`: original user-provided gene symbol or token
- `mapping_status`: whether the gene matched exactly, by alias, by best match, or failed
- `mapping_note`: short explanation of mapping behavior
- `approved_symbol`: resolved preferred symbol
- `full_gene_name`: resolved full gene name
- `aliases`: known aliases or alternate symbols
- `entrez_id`: Entrez Gene ID
- `ensembl_gene_id`: Ensembl gene ID
- `uniprot_accession`: canonical reviewed UniProt accession where available

### Gene and protein identity
- `protein_name`: recommended protein name
- `protein_length`: protein sequence length
- `gene_type`: broad gene class such as protein-coding
- `chromosome`, `cytoband`, `strand`, `genomic_location`: genomic context

### Function and localization
- `reported_function`: concise one-sentence summary, not a raw long text dump
- `subcellular_location`: brief localization summary
- `molecular_function`, `biological_process`, `cellular_component`: concise GO-derived summaries

### Domain-related fields
- `domain_annotation`: semicolon-separated domain or feature list
- `domain_coordinates`: domain names with residue coordinates
- `key_domain_summary`: concise domain architecture summary
- `functional_domain_class`: broad interpreted class such as reader, writer, eraser, kinase, transcription factor, RNA-binding protein, or membrane protein

### Pathway-related fields
- `reactome_pathways`: concise Reactome pathway memberships
- `kegg_pathways`: concise KEGG pathway memberships
- `hallmark_pathways`: reserved column for future Hallmark membership support
- `pathway_annotation_summary`: short merged pathway summary

### Disease and interpretation fields
- `disease_association`: concise disease names, preferring readable full names over acronyms
- `protein_class`: broad class inferred from keywords and domain context
- `druggability`: coarse interpretive label such as potentially_druggable or challenging_but_relevant
- `source_summary`: concise list of annotation sources used

## Annotation principles

This skill should prioritize:
- concise and readable summaries over long copied text
- structured identifiers over prose
- full disease names over unexplained acronyms
- domain architecture as a first-class annotation
- clean pathway splitting by source where possible
- deterministic table output

Avoid:
- long raw paragraphs copied from source databases
- duplicated information across columns
- unexplained abbreviations when fuller wording is available
- mixing expression and dependency analysis into this skill

## Workflow

1. Read the input gene list.
2. Normalize gene tokens and remove empty entries.
3. Resolve the organism.
4. Query core gene metadata and identifier mapping.
5. Determine the best gene match and mapping status.
6. Extract identifiers, genomic location, aliases, and basic metadata.
7. Resolve the UniProt accession if available.
8. Retrieve protein-level annotations.
9. Extract protein name, length, functional text, localization, keywords, and features.
10. Convert raw function text into a concise function summary.
11. Extract domain annotations and domain coordinates.
12. Infer a broad functional domain class when possible.
13. Split pathway annotations into separate source-specific columns when available.
14. Summarize disease associations into short readable names.
15. Infer broad protein class and rough druggability label when possible.
16. Assemble the final table using the default schema.
17. Write the result table to TSV or CSV.

## Interaction model

Users may describe their needs in natural language. Translate supported requests into appropriate script arguments or output expectations.

Examples:
- “Annotate these genes”  
  → run the skill with default schema

- “Include domain information”  
  → keep domain fields in output

- “Give me concise function summaries”  
  → preserve short `reported_function` behavior

- “Map these genes to UniProt and Ensembl”  
  → ensure identifier columns are included

- “Output a CSV”  
  → use `--output-format csv`

- “Use mouse genes”  
  → set `--organism mouse`

## Execution

Run the main script:

`python scripts/annotation_for_gene_list.py --input <gene_file> --output <output_tsv>`

### Common examples

Annotate a human gene list and output TSV:

`python scripts/annotation_for_gene_list.py --input examples/example_genes.txt --output results/example_annotations.tsv --organism human`

Annotate a mouse gene list and output CSV:

`python scripts/annotation_for_gene_list.py --input examples/example_genes.txt --output results/example_annotations.csv --organism mouse --output-format csv`

### Supported arguments

Required:
- `--input <gene_file>`
- `--output <output_path>`

Optional:
- `--organism <human|mouse|hs|hsapiens|homo sapiens|mm|mus musculus>`
- `--output-format <tsv|csv>`

## Defaults

Unless the user specifies otherwise:
- use `human` as the organism
- write TSV output
- include the full default annotation schema
- keep function, disease, and pathway fields concise
- do not attempt expression, essentiality, or TCGA analysis

## Source behavior

This skill is built around source-backed gene and protein annotation retrieval.

Core source categories include:
- gene-level identifier and metadata source
- protein-level source for reviewed protein annotation
- source-specific pathway fields when available
- protein feature and disease comments when available

When a field is unavailable from the sources used, leave it blank rather than hallucinating content.

## Output conventions

- one row per gene
- one column per annotation field
- preserve input order as much as possible through row processing
- use concise summaries rather than full copied paragraphs
- prefer full disease names over acronyms where possible
- keep pathway source columns separate when possible
- keep the workflow deterministic and reproducible

## Failure modes

Fail clearly when:
- the input gene list is empty
- the organism is unsupported
- the output path is not writable

Handle gracefully when:
- a gene cannot be mapped
- a UniProt accession cannot be found
- pathway annotations are sparse
- hallmark membership is unavailable
- disease annotations are absent
- one or more fields are missing from source records

When a gene cannot be fully annotated:
- still return a row
- set `mapping_status` appropriately
- explain the issue in `mapping_note`

## Robot usage notes

When another robot uses this skill:
- do not overload the output with raw copied source text
- prefer concise summaries in `reported_function`
- prefer explicit identifiers and domain fields
- leave unsupported fields blank rather than guessing
- do not use this skill for expression, essentiality, or TCGA analysis
- do not silently convert this into literature curation or enrichment analysis
- return the saved file path in the final response

## Trigger examples

This skill should trigger for prompts like:
- “Annotate this gene list.”
- “Give me IDs and protein information for these genes.”
- “Add domain annotation for these genes.”
- “Summarize what these genes do.”
- “Map these symbols to Ensembl and UniProt.”
- “Give pathway and disease information for these genes.”

## Non-trigger examples

This skill should not trigger for prompts like:
- “Run GO enrichment on this gene list.”
- “Run GSEA on this ranked list.”
- “Design sgRNAs for these genes.”
- “Find essential genes in DepMap.”
- “Analyze TCGA survival for these genes.”
- “Curate a literature-supported ferroptosis gene list.”