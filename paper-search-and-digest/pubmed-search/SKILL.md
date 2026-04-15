---
name: pubmed-search
description: Search PubMed for scientific papers matching a query. Retrieves ranked results with abstracts, journals, citation counts (via Europe PMC), and author information. Filters by date range, journal type, publication type, and study type. Saves a structured TSV and a formatted report.
---

# PubMed Search

## Purpose
Execute targeted searches of the PubMed database to find relevant papers on any biomedical topic. Returns ranked results with rich metadata (abstracts, citations, MeSH terms, authors) and generates formatted reports plus visualizations of publication trends.

## When to Use / When NOT to Use

**Use when:**
- You need to find papers on a specific research topic
- Building a comprehensive literature review from scratch
- Discovering recent publications in your field
- Filtering papers by publication type (reviews, clinical trials, meta-analyses)
- Analyzing publication trends over time
- Finding papers from specific journals or date ranges

**Do NOT use when:**
- You only need to digest a single known paper (use `paper-digest-single`)
- You want real-time trending topics (use `preprint-tracker` for recent preprints)
- Your query is so broad it returns >10,000 results (refine query with MeSH terms)
- You need full-text PDF content analysis beyond abstracts

## Expected Inputs & Outputs

**Inputs:**
- `--query STRING` — PubMed query string (required; supports MeSH, Author, Journal filters)
- `--max-results INT` — Maximum papers to retrieve (default: 50)
- `--date-from` — Start date (YYYY or YYYY/MM/DD; default: 2020/01/01)
- `--date-to` — End date (YYYY or YYYY/MM/DD; default: today)
- `--pub-type` — Filter by publication type: review, clinical_trial, meta_analysis, original
- `--journal-filter` — Comma-separated journal names to include
- `--sort-by` — Sort results: relevance (default), date, citations
- `--fetch-citations` — Flag: fetch citation counts from Europe PMC (slower, richer)
- `--abstract-words` INT — Max abstract preview length (default: 100)

**Outputs (in --outdir):**
- `search_results.tsv` — Tabular results: pmid, title, authors, journal, year, doi, abstract_preview, pub_type, citations
- `search_report.md` — Formatted markdown list of papers with full citations
- `keyword_cloud.png` — Bar chart of top 20 keywords across all abstracts
- `publication_timeline.png` — Papers per year bar chart
- `journal_distribution.png` — Pie chart of papers by journal
- `search_log.txt` — Query parameters and processing notes

## Procedure

1. **Execute PubMed eSearch**
   - Query PubMed with supplied search string, date range, date type
   - Get list of PMIDs matching query
   - Limit to --max-results

2. **Batch Fetch Metadata (eFetch)**
   - Download XML metadata for all matched PMIDs in batches of 200
   - Parse: title, authors, journal, year, volume, pages, DOI, abstract, publication types, MeSH terms

3. **Extract Abstract Keywords**
   - Tokenize all abstracts
   - Remove English stopwords
   - Calculate term frequency (TF-IDF simplified)
   - Identify top 20 most relevant keywords

4. **Fetch Citation Counts (optional)**
   - If --fetch-citations flag, query Europe PMC API
   - Get citedByCount for each paper

5. **Sort Results**
   - By relevance (default), publication date (newest first), or citation count (highest first)

6. **Generate Outputs**
   - TSV: machine-readable tabular format
   - Markdown: human-readable formatted list
   - Visualizations: publication trends, keyword cloud, journal distribution

## Key Execution Patterns

**Simple search, last 2 years:**
```bash
python pubmed_search.py --query "CRISPR gene therapy" --max-results 50 --outdir ./results
```

**Filter by publication type (reviews only):**
```bash
python pubmed_search.py --query "cancer immunotherapy" --pub-type review --max-results 100 --outdir ./reviews
```

**Clinical trials from top journals:**
```bash
python pubmed_search.py --query "CAR-T therapy" --pub-type clinical_trial \
  --journal-filter "Nature,Science,Lancet" --outdir ./trials
```

**Sort by citation count (most cited first):**
```bash
python pubmed_search.py --query "GWAS risk factors" --sort-by citations \
  --fetch-citations --max-results 100 --outdir ./cited
```

**Custom date range with MeSH filtering:**
```bash
python pubmed_search.py --query "Alzheimer's Disease[MeSH]" \
  --date-from 2022/01/01 --date-to 2024/12/31 --max-results 200 --outdir ./results
```

**Export as JSON for programmatic use:**
```bash
python pubmed_search.py --query "transcriptomics single cell" --output-format json --outdir ./data
```

## Parameter Decision Guide

| Parameter | When to Use | Example |
|-----------|------------|---------|
| `--query` | Always required; use MeSH terms for precision | `--query "Diabetes Mellitus[MeSH]"` |
| `--max-results` | Start with 50-100, increase if comprehensive review needed | `--max-results 200` |
| `--date-from` | Restrict to recent papers (e.g., last 2-3 years for fast-moving fields) | `--date-from 2022/01/01` |
| `--pub-type` | Filter to review articles for quick synthesis | `--pub-type review` |
| `--journal-filter` | When you trust only high-impact journals | `--journal-filter "Nature,Science,Cell"` |
| `--sort-by citations` | Find most influential papers in the field | `--sort-by citations --fetch-citations` |
| `--sort-by date` | Track recent developments in a fast-moving area | `--sort-by date` |
| `--fetch-citations` | When citation impact matters (adds 30-60 sec per search) | `--fetch-citations` |

## Failure Modes & Recovery

| Error | Cause | Solution |
|-------|-------|----------|
| "No results found" | Query is too specific or uses incorrect MeSH | Try simpler keywords, remove MeSH filters |
| "Query syntax error" | Malformed PubMed query string | Test query on pubmed.gov first, then copy exact syntax |
| ">10,000 results found" | Query too broad; will retrieve only first N results | Add year filter, specific MeSH terms, or author names |
| "Empty abstracts" | Some papers lack abstracts (older, review types) | Expected; filters still work on title and metadata |
| "Citation counts unavailable" | Europe PMC API timeout or down | Retry with --fetch-citations, or skip citations |
| "Memory error on large results" | Fetching too many papers at once (>500) | Reduce --max-results or split query into smaller searches |
| "Duplicate papers in TSV" | Different PMID versions of same paper | Expected; manually deduplicate if needed |

## Requirements
See `requirements.txt` for exact versions. Core dependencies:
- `requests` — HTTP API calls to PubMed and Europe PMC
- `pandas` — Data organization and TSV export
- `matplotlib` — Visualization of trends and distributions
- `xml.etree.ElementTree` — PubMed XML parsing (stdlib)
- `json` — Metadata serialization (stdlib)
- `re` — Regular expressions for text processing (stdlib)
