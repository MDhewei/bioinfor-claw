---
name: search_big_labs_by_field
description: Discover leading PI-led research groups in a scientific field using OpenAlex publication evidence, with ranking based mainly on overall scholarly impact, PI-like authorship patterns, and recent 10-year and 5-year activity summaries.
---

# Search Big Labs by Field

## Purpose

Use this skill when the user wants to identify important research labs, principal investigators, and institutions in a scientific field.

Typical use cases:
- find major labs in CRISPR, bioinformatics, or computational biology
- build a postdoc or collaboration target list
- understand the lab landscape of a field
- compare overall impact and recent activity of candidate PIs

This skill is **publication-based** and **PI-oriented**:
- it retrieves field-relevant papers
- aggregates authors from those papers
- filters for PI-like candidates
- profiles candidates with overall impact and recent publication summaries
- de-duplicates repetitive results

The final output is a **ranked candidate lab list**, where each row is a best-effort PI-centered lab proxy.

---

## When to use

Use this skill for questions like:
- “Find the big labs in CRISPR.”
- “Who are the top PIs in computational biology?”
- “Give me strong bioinformatics labs.”
- “Show leading groups in cancer genomics.”
- “Find important labs in single-cell computational biology.”
- “What are the main labs working on protein structure prediction?”

---

## When not to use

Do not use this skill when the user wants:
- a single official lab website
- exact lab size or membership
- a faculty directory lookup
- exact email scraping
- a paper summary or weekly paper report
- gene enrichment, annotation, or CRISPR design

This skill is for **lab discovery**, not for website lookup or administrative verification.

---

## Data source

This skill uses:
- **OpenAlex works**
- **OpenAlex authors**
- **OpenAlex institutions**

Current version:
- does **not** scrape lab websites
- does **not** search the web for official lab pages
- relies on structured publication and affiliation data only

---

## What the skill does

The workflow has seven main stages.

### Step 1. Retrieve field-relevant papers
The script searches OpenAlex works using the user keyword and a configurable recent time window.

Goal:
- build a discovery corpus of papers relevant to the requested field

### Step 2. Remove unsuitable papers
The script excludes:
- unsupported work types
- very large collaboration papers above a configurable authorship threshold

Goal:
- reduce noise from consortium papers and massive coauthor lists

### Step 3. Aggregate authors
From the filtered papers, the script aggregates author-level evidence:
- relevant paper count
- citations across relevant papers
- authorship positions
- institutions seen
- matched keyword breadth
- candidate representative works

Goal:
- convert the paper set into a candidate PI pool

### Step 4. Keep PI-like candidates
The script keeps only authors with sufficient PI-like authorship evidence, mainly:
- minimum last-author count
- or minimum first/last-author count

Goal:
- reduce inflation from middle authors and frequent collaborators

### Step 5. Build PI profiles
For each surviving candidate, the script retrieves:
- overall citation count
- overall works count
- institution metadata
- recent 10-year and 5-year publication summaries
- representative papers

Goal:
- provide a richer profile for each candidate

### Step 6. Rank by overall impact
Final ranking emphasizes:
- overall scholarly impact
- total work count
- field relevance
- PI-like authorship evidence

Recent publication counts are included for interpretation, but they are **not** the main ranking signal.

Goal:
- rank “important labs overall,” not merely the most recently prolific ones

### Step 7. De-duplicate final results
The script reduces redundancy by:
- limiting rows per institution
- limiting excessive overlap of representative paper sets across rows

Goal:
- improve diversity and readability of the final candidate list

---

## Inputs

## Required arguments

### `--keyword`
The scientific field or topic to search.

Examples:
- `CRISPR`
- `bioinformatics`
- `computational biology`
- `single-cell genomics`
- `cancer genomics`
- `protein structure prediction`

### `--output`
Path to the TSV output file.

Example:
- `results/computational_biology_labs.tsv`

---

## Optional arguments

### `--top-n`
Number of final candidate labs / PIs to keep.

Default:
- `30`

Use a larger value when the user wants a broader landscape.
Use a smaller value when the user wants a short, high-confidence shortlist.

---

### `--max-works`
Maximum number of keyword-relevant works retrieved from OpenAlex during the initial discovery stage.

Default:
- `300`

Higher values:
- improve recall
- increase runtime
- may increase noise for broad keywords

Lower values:
- run faster
- may miss some important candidates

---

### `--per-page`
Pagination size for OpenAlex works retrieval.

Default:
- `100`

Usually this can be left unchanged.

---

### `--years-back`
How many recent years of works are used for the initial discovery stage.

Default:
- `5`

Meaning:
- the initial paper pool is built from this recent window

Use smaller values if:
- the user wants the current landscape only

Use larger values if:
- the field has important older foundational papers
- famous PIs may otherwise be missed

For example:
- `5` is often good for broad current fields
- `10` is often better for CRISPR and other iconic mature fields

---

### `--author-works-limit`
Maximum number of recent works retrieved per candidate PI during the profiling stage.

Default:
- `150`

This affects:
- recent paper summaries
- representative papers
- recent activity fields

Higher values:
- provide richer profiles
- increase runtime

Lower values:
- run faster
- may reduce profile completeness

---

### `--sleep-seconds`
Pause between detail requests to reduce API pressure.

Default:
- `0.2`

Increase if:
- API requests appear throttled

Decrease slightly if:
- speed matters and the API remains stable

---

### `--max-authors-per-paper`
Exclude papers with more than this many authorships.

Default:
- `20`

This is one of the most important quality-control parameters.

Purpose:
- remove very large collaboration papers
- reduce repeated institutions and paper clusters
- reduce coauthor inflation

Use stricter values for:
- broad fields like `computational biology` or `bioinformatics`

Use looser values for:
- collaboration-heavy fields like `CRISPR`

Examples:
- `20` for broad/noisy fields
- `60` for CRISPR-like fields

---

### `--min-last-author-count`
Minimum last-author count for PI-like filtering.

Default:
- `2`

Purpose:
- favor likely lab leaders rather than collaborators

Higher values:
- stricter PI filtering
- more precision
- risk missing some important PIs

Lower values:
- better recall
- more noise

---

### `--min-first-last-author-count`
Minimum combined first/last-author count for PI-like filtering.

Default:
- `3`

Purpose:
- rescue candidates who may not have many last-author papers but still show strong leadership-style authorship

This works together with `--min-last-author-count`.

---

### `--max-per-institution`
Maximum number of final rows retained per institution.

Default:
- `2`

Purpose:
- prevent one institution from dominating the final list

Lower values:
- more institution diversity

Higher values:
- more depth within dominant institutions

---

### `--max-paper-overlap`
Maximum allowed DOI overlap between selected rows.

Default:
- `2`

Purpose:
- avoid highly repetitive candidate rows
- reduce cases where multiple coauthors from one paper cluster all appear in the final list

Lower values:
- stronger diversification

Higher values:
- more tolerance for similar paper sets

---

### `--output-jsonl`
Optional JSONL output path.

Use this when a machine-readable row-by-row format is also needed.

Example:
- `results/computational_biology_labs.jsonl`

---

## Outputs

## Main output
A TSV file with one row per candidate lab / PI.

## Optional output
A JSONL file with one record per row.

## Typical columns

### Identity and affiliation
- `query_keyword`
- `lab_name`
- `pi_name`
- `institution_name`
- `institution_type`
- `city`
- `country_code`
- `institution_homepage`
- `orcid`
- `openalex_author_id`
- `institution_openalex_id`
- `institution_ror`

### Ranking and evidence
- `relevance_score`
- `overall_cited_by_count`
- `overall_works_count`
- `query_relevant_works_count`
- `query_total_citations_from_relevant_works`
- `first_last_author_count`
- `last_author_count`
- `first_author_count`
- `matched_keyword_terms`

### Topic summary
- `lab_type`
- `main_methods`
- `main_topics`
- `top_topics`
- `topic_share`
- `research_summary`

### Recent activity
- `publications_last_10y`
- `publications_last_5y`
- `publications_last_3y`
- `keyword_relevant_papers_last_10y`
- `keyword_relevant_papers_last_5y`
- `keyword_relevant_papers_last_3y`
- `publication_year_distribution_last_10y`

### Representative papers
- `representative_papers_overall_10y`
- `representative_keyword_relevant_papers`
- `recent_papers_last_5y`

---

## How to interpret the output

### `lab_name`
Usually a convenience label like:
- `<PI name> Lab`

It is not guaranteed to be the official lab title.

### `relevance_score`
A heuristic ranking score combining:
- field relevance
- overall impact
- PI-like authorship evidence

Useful for:
- sorting
- comparison
- shortlist construction

Not intended as:
- a universal bibliometric metric

### `overall_cited_by_count`
Overall OpenAlex citations for the candidate PI.

### `overall_works_count`
Overall OpenAlex works for the candidate PI.

### `publications_last_10y` and `publications_last_5y`
Recent activity summaries.
These are included for context, not as the dominant ranking driver.

### `representative_papers_overall_10y`
Top overall representative papers from the last 10 years.

### `representative_keyword_relevant_papers`
Top representative papers specifically relevant to the queried field.

### `recent_papers_last_5y`
Representative recent papers to help judge current activity.

---

## Execution

## Entry point

```bash
python scripts/search_big_labs_by_field.py --keyword "<field>" --output <output_tsv>