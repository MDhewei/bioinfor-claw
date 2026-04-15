---
name: bioinformatics_big_papers_weekly_report
description: Use when the user wants a weekly report of recent high-impact bioinformatics papers from selected journals, ranked by journal prestige, generic bioinformatics relevance, and optional user-interest keywords, with a PDF report and tabular output.
---

# Bioinformatics Big Papers Weekly Report

Use this skill to retrieve recent papers from a configured journal scope, rank them for likely bioinformatics relevance and potential impact, summarize the top papers, and generate a polished weekly PDF report plus a tabular results file.

## Purpose

This skill is designed for weekly literature surveillance in bioinformatics and computational biology. It helps users identify “big papers” from a recent time window, especially from high-profile journals or user-selected journals, and produces a concise report suitable for reading, sharing, or archiving.

This skill is intended for users who want:
- a weekly scan of recent bioinformatics papers
- a report focused on high-profile journals
- ranking that reflects both field relevance and user interests
- one PDF report with concise paper-by-paper summaries
- a machine-readable TSV table of the selected papers

This skill is especially useful for:
- lab literature roundup
- PI or group weekly reading reports
- topic-focused surveillance, such as single-cell, spatial, CRISPR, or protein structure
- monitoring methods, resources, benchmark papers, and computational trends

## Use when

Use this skill when the user:
- wants recent important papers in bioinformatics
- wants papers from the past week or another recent date range
- wants a PDF report summarizing recent papers
- wants a ranked list of papers from selected journals
- wants to prioritize papers by topic keywords
- wants a weekly or periodic “big papers” report
- wants a literature scan from journals such as Nature, Science, Cell, Nature Methods, Nature Biotechnology, Nature Genetics, Cell Systems, or Cell Genomics

Typical requests:
- “Find the big bioinformatics papers from the past week.”
- “Make a weekly PDF report of important recent bioinformatics papers.”
- “Focus on Nature Methods and Nature Biotechnology.”
- “Rank papers by relevance to single-cell and spatial omics.”
- “Summarize the most important computational biology papers this week.”

## Do not use when

Do not use this skill when the user:
- wants a comprehensive systematic review
- wants full-text paper analysis rather than title/abstract-level summarization
- wants enrichment analysis, gene annotation, or CRISPR design
- wants only a raw list of PubMed hits without scoring or report generation
- wants a bibliography manager export as the primary output
- wants detailed extraction from PDFs of full papers
- wants only one known paper summarized rather than a weekly report

If the user wants:
- full-text deep review of specific papers
- citation network analysis
- meta-analysis
- topic modeling over a large literature corpus
- manuscript drafting based on a literature set

then use a dedicated literature-analysis workflow instead of this skill.

## What this skill does

This skill performs the following tasks:

1. retrieves recent journal articles from a configured journal scope
2. collects title, metadata, and abstract text when available
3. ranks papers using a transparent heuristic scoring system
4. optionally prioritizes papers matching user-interest keywords
5. summarizes the top papers
6. writes a structured TSV results table
7. writes a polished PDF report

## Report scope

By default, this skill is intended to cover recent papers from a curated set of high-profile journals in biology, medicine, and computational biology.

The default journal scope may include journals such as:
- Nature
- Science
- Cell
- Nature Biotechnology
- Nature Methods
- Nature Genetics
- Nature Medicine
- Nature Communications
- Nature Machine Intelligence
- Nature Computational Science
- Cell Systems
- Cell Genomics
- Science Advances

The user may also supply:
- additional journals
- a fully custom journal list
- a flag to replace the default journal list with only user-specified journals

## Expected inputs

Required:
- output prefix

Optional but commonly used:
- date range
- number of papers to keep
- user-interest keywords
- journal list
- whether to use the default journals or replace them
- whether to skip PubMed enrichment
- number of candidate papers to enrich with PubMed abstracts
- scoring weights

## Expected outputs

This skill produces:

1. a TSV file containing ranked papers and metadata  
   Example:
   - `<output_prefix>.tsv`

2. a PDF report summarizing the selected papers  
   Example:
   - `<output_prefix>.pdf`

The TSV usually includes:
- DOI
- title
- journal
- publication date
- authors
- abstract
- matched keywords
- scoring components
- summary fields

The PDF usually includes:
- report title
- coverage window
- journal scope
- scoring method
- summary table of selected papers
- one section per selected paper
- concise interpretation of why each paper matters

## Defaults

Unless the user specifies otherwise:
- date window = last 7 days
- top papers in final report = 10
- default journal whitelist is used
- user journals are appended unless the user explicitly replaces defaults
- generic bioinformatics score weight = 1.0
- user-interest keyword score weight = 2.0
- Crossref metadata is used first
- PubMed enrichment is used only for top candidate papers in fast mode
- output includes both TSV and PDF

## Scoring model

This skill uses a heuristic impact score for triage and prioritization.

The score is not a citation-based impact metric and should not be interpreted as a formal measure of scientific importance.

### Score components

The impact score is composed of:

- `journal_score`
- `generic_bioinfo_score × generic_bioinfo_weight`
- `interest_keyword_score × interest_weight`
- `citation_bonus`
- `abstract_bonus`
- `bonus_score`

### Journal score

A manual journal prestige score is assigned based on journal identity. Higher-profile journals receive higher baseline values.

### Generic bioinformatics score

This score reflects whether the paper appears likely to be relevant to bioinformatics or computational biology based on broad built-in keywords in:
- title
- abstract
- subject metadata

Examples of generic signal categories include:
- computational methods
- software and tools
- benchmarks
- resources
- atlases
- machine learning
- deep learning
- single-cell
- spatial
- genomics
- proteomics
- transcriptomics
- protein structure
- CRISPR

### Interest keyword score

This score reflects user-specified interests. If the user provides keywords such as:
- single-cell
- spatial
- foundation model
- protein structure

then papers matching those keywords receive additional score.

This component is intentionally weighted more strongly than the generic bioinformatics score by default, because it reflects the user’s current interests.

### Citation bonus

A small citation-count bonus is included when Crossref citation metadata is available.

### Abstract bonus

A small bonus is included when an abstract is available.

### Extra bonus

A small additional bonus is applied to papers whose titles suggest especially reusable or high-impact resource types, such as:
- benchmark
- atlas
- resource

## Keyword matching behavior

This skill supports two types of keyword matching:

### User-interest keywords
These come from:
- `--interest-keywords`
- `--interest-keywords-file`

These keywords are used to personalize ranking.

### Generic bioinformatics keywords
These are built into the script and are used to estimate field relevance.

### Displayed matched keywords

The report shows a matched-keywords column for transparency.

Behavior:
- if user-interest keywords are provided, the display prioritizes matched user-interest keywords
- if no user-interest keywords are provided, the display falls back to matched generic bioinformatics keywords

This prevents the keyword column from appearing blank in normal runs.

## Fast mode and performance behavior

This skill is optimized for practical weekly use and may run in a fast mode.

### Fast mode workflow

1. query recent papers from the configured journals using Crossref
2. perform a first-pass ranking using Crossref metadata
3. keep likely candidates
4. optionally enrich only the top candidate subset with PubMed abstracts
5. rerank and finalize the report

### Why this is useful

This design reduces runtime substantially compared with fetching PubMed abstracts for every paper.

### PubMed enrichment

PubMed enrichment is optional and can be skipped.

If enabled, it is typically limited to the top candidate subset, controlled by:
- `--pubmed-top-k`

### Cache behavior

A DOI-to-abstract cache may be used to avoid repeated PubMed retrieval across runs.

## Journal selection behavior

This skill supports three journal-selection modes:

### 1. Default journals only
Use the built-in journal whitelist.

### 2. Default journals plus user journals
Append user-supplied journals to the default list.

### 3. User journals only
Use:
- `--replace-default-journals`

to restrict the scan to only the user-provided journals.

## Execution

Run the main script:

`python scripts/bioinformatics_big_papers_weekly_report.py --output-prefix <prefix>`

### Common examples

Run a default weekly report:

`python scripts/bioinformatics_big_papers_weekly_report.py --output-prefix results/bioinfo_weekly_report`

Run a topic-focused weekly report:

`python scripts/bioinformatics_big_papers_weekly_report.py --output-prefix results/bioinfo_weekly_report --interest-keywords "single-cell,spatial,foundation model"`

Use only a custom journal set:

`python scripts/bioinformatics_big_papers_weekly_report.py --output-prefix results/bioinfo_weekly_report --journals "Nature Biotechnology,Nature Methods,Cell Systems,Cell Genomics" --replace-default-journals`

Use a custom journal file:

`python scripts/bioinformatics_big_papers_weekly_report.py --output-prefix results/bioinfo_weekly_report --journal-list-file journals.txt --replace-default-journals`

Run in fast Crossref-only mode:

`python scripts/bioinformatics_big_papers_weekly_report.py --output-prefix results/bioinfo_weekly_report --skip-pubmed`

Limit PubMed enrichment to top candidates:

`python scripts/bioinformatics_big_papers_weekly_report.py --output-prefix results/bioinfo_weekly_report --pubmed-top-k 15`

Adjust scoring weights:

`python scripts/bioinformatics_big_papers_weekly_report.py --output-prefix results/bioinfo_weekly_report --interest-keywords "protein structure,alphafold" --generic-bioinfo-weight 1.0 --interest-weight 3.0`

## Supported arguments

Required:
- `--output-prefix <prefix>`

Optional:
- `--date-from <YYYY-MM-DD>`
- `--date-to <YYYY-MM-DD>`
- `--top-n <int>`
- `--rows-per-journal <int>`
- `--email <email>`
- `--ncbi-api-key <key>`
- `--journals <comma-separated journal list>`
- `--journal-list-file <path>`
- `--replace-default-journals`
- `--interest-keywords <comma-separated keywords>`
- `--interest-keywords-file <path>`
- `--generic-bioinfo-weight <float>`
- `--interest-weight <float>`
- `--skip-pubmed`
- `--pubmed-top-k <int>`
- `--cache-file <path>`
- `--crossref-workers <int>`

## Output columns

The TSV output may include fields such as:
- `doi`
- `title`
- `journal`
- `pub_date`
- `authors`
- `subjects`
- `citation_count`
- `abstract`
- `matched_interest_keywords`
- `matched_generic_keywords`
- `matched_keywords_display`
- `interest_keyword_score`
- `generic_keyword_count`
- `journal_score`
- `generic_bioinfo_score`
- `citation_bonus`
- `abstract_bonus`
- `bonus_score`
- `impact_score`
- `why_it_matters`
- `major_discoveries`
- `significance`

## PDF structure

The PDF report should include:

### First page
- report title
- coverage window
- journal scope
- executive summary
- impact score method

### Summary table
- rank
- title
- journal
- date
- matched keywords
- impact score

### Per-paper sections
- title
- journal
- date
- DOI
- authors
- matched keywords
- score breakdown
- why it matters
- major discoveries
- likely impact
- abstract excerpt

## Report-writing principles

This skill should prioritize:
- concise summaries
- clear ranking transparency
- readable PDF layout
- explicit explanation of scoring
- practical usefulness for weekly literature scanning

Avoid:
- long copied abstracts as the main summary
- opaque scoring with no explanation
- blank keyword columns
- overclaiming impact
- implying that the heuristic score is a formal scientific metric

## Parameter decision guide

| Signal in user request | Parameter to set |
|---|---|
| "past week" (default) | use defaults (last 7 days) |
| "past two weeks" | `--date-from YYYY-MM-DD` (14 days ago) `--date-to YYYY-MM-DD` (today) |
| "top 5 papers" (brief report) | `--top-n 5` |
| "top 20 papers" (comprehensive) | `--top-n 20` |
| "focus on single-cell / spatial / protein structure" | `--interest-keywords "single-cell,spatial omics"` |
| "only Nature Methods and Nature Biotechnology" | `--journals "Nature Methods,Nature Biotechnology" --replace-default-journals` |
| "add Cell Genomics to the list" | `--journals "Cell Genomics"` (appended to defaults) |
| "fast mode / no PubMed" | `--skip-pubmed` |
| "richer abstracts for top papers" | `--pubmed-top-k 20` |
| "strongly prioritize user topic" | `--interest-weight 5.0` |
| "balance field relevance and user interest equally" | `--generic-bioinfo-weight 1.0 --interest-weight 1.0` |

## Interaction model

Translate natural-language user requests into supported arguments where possible.

Examples:
- “Focus on single-cell and spatial papers.”  
  → `--interest-keywords "single-cell,spatial"`

- “Only use Nature Methods and Nature Biotechnology.”  
  → `--journals "Nature Methods,Nature Biotechnology" --replace-default-journals`

- “Give me a faster version.”  
  → use `--skip-pubmed` or smaller `--pubmed-top-k`

- “Make this a broader weekly scan.”  
  → keep default journals and no user keyword restriction

- “Focus on protein structure papers.”  
  → `--interest-keywords "protein structure,alphafold"`

## Robot usage notes

When another robot uses this skill:
- preserve the PDF + TSV dual output
- explain that the impact score is heuristic
- do not present the score as a citation-based or objective importance metric
- use user-interest keywords when supplied
- display matched generic keywords if no user-interest keywords are supplied
- prefer fast mode for routine use unless deeper PubMed enrichment is needed
- return the output file paths in the final response
- do not use this skill for full-text literature review or systematic review tasks

## Failure modes

Fail clearly when:
- no journals are available after parsing defaults and user input
- no records are retrieved from the selected journals and date range
- no likely bioinformatics papers remain after filtering
- output paths cannot be written

Handle gracefully when:
- some journal queries fail
- Crossref rate limits occur
- PubMed enrichment fails for some papers
- abstracts are unavailable
- some metadata fields are missing

If some journal queries fail but others succeed:
- continue with successful journals
- still produce the report if enough papers remain

## Performance notes

This skill may be slowed by:
- Crossref rate limits
- PubMed abstract retrieval
- large journal scopes
- large per-journal row counts

Recommended performance-friendly usage:
- moderate `rows-per-journal`
- limited `pubmed-top-k`
- cache enabled
- low or moderate `crossref-workers`

## Trigger examples

This skill should trigger for prompts like:
- “Find the big bioinformatics papers from the past week.”
- “Make a weekly PDF report of important recent bioinformatics papers.”
- “Focus on recent Nature Methods and Nature Biotechnology papers.”
- “Give me a report on recent bioinformatics papers about single-cell and spatial omics.”
- “Summarize the most important recent computational biology papers.”

## Non-trigger examples

This skill should not trigger for prompts like:
- “Annotate this gene list.”
- “Run GO enrichment.”
- “Run GSEA on this ranked list.”
- “Design sgRNAs for these genes.”
- “Do TCGA analysis for this gene.”
- “Summarize this one specific paper in depth.”