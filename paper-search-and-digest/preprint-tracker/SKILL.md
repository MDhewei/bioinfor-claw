---
name: preprint-tracker
description: Search and track recent preprints from bioRxiv and medRxiv on a given topic. Retrieves preprint metadata, abstracts, and author information via the bioRxiv/medRxiv API. Identifies trending topics, prolific authors, and generates a digest report.
---

# Preprint Tracker

## Purpose
Monitor emerging research trends by searching recent preprints on bioRxiv (life sciences) and medRxiv (medical/clinical research). Identifies trending topics, most active authors, and generates formatted reports to stay ahead of published literature.

## When to Use / When NOT to Use

**Use when:**
- You want to track the very latest research before peer review and formal publication
- Monitoring a field for rapid developments (e.g., new treatments, emerging pathogens)
- Identifying prolific authors and research groups in a specific area
- Building a curated digest of preprints on a topic for a journal club or group
- Finding cutting-edge methods before they appear in peer-reviewed journals

**Do NOT use when:**
- You need established, peer-reviewed evidence (use `pubmed-search` instead)
- You want to cite papers (preprints are not peer-reviewed; use published versions)
- You need papers older than 3-6 months (bioRxiv/medRxiv archive is limited)
- Your focus is on validating conclusions (preprints lack peer review filtering)

## Expected Inputs & Outputs

**Inputs:**
- `--query STRING` — Search keywords (required; searches title and abstract)
- `--server` — Source: biorxiv, medrxiv, or both (default: both)
- `--date-from STR` — Start date (YYYY-MM-DD; default: 30 days ago)
- `--date-to STR` — End date (YYYY-MM-DD; default: today)
- `--max-results INT` — Maximum preprints to retrieve (default: 100)
- `--category STR` — bioRxiv category filter (e.g., "genomics", "bioinformatics", "cancer biology")
- `--sort-by` — Sort by: date (default, newest first) or relevance
- `--output-format` — Choices: tsv, markdown (default: markdown)

**Outputs (in --outdir):**
- `preprints.tsv` — Tabular results: doi, title, authors, category, date, word_count, abstract_preview
- `trending_topics.tsv` — Top keywords with frequency and top 3 papers using each keyword
- `digest_report.md` — Formatted preprint summaries with full metadata
- `category_distribution.png` — Bar chart of preprints by bioRxiv category
- `timeline.png` — Preprints published per week (bar chart)
- `topic_frequency.png` — Bar chart of top 20 keywords
- `prolific_authors.tsv` — Authors ranked by number of preprints
- `search_log.txt` — Processing notes and summary statistics

## Procedure

1. **Query bioRxiv/medRxiv API**
   - Call `https://api.biorxiv.org/details/<server>/<date_from>/<date_to>/<cursor>/json`
   - Paginate through results (cursor: 0, 100, 200, etc.)
   - Extract: doi, title, authors, category, abstract, date, jatsxml

2. **Filter by Keywords**
   - Search title and abstract text for --query keywords (case-insensitive)
   - Apply --category filter if provided

3. **Extract Metadata**
   - Parse author list (first, last, affiliations where available)
   - Count word count in abstract
   - Parse publication date and category

4. **Analyze Trending Topics**
   - Tokenize all abstracts and titles
   - Remove stopwords
   - Calculate word frequency
   - Identify top 20 keywords

5. **Identify Prolific Authors**
   - Count papers per author (both first author and all authors)
   - Rank by publication count

6. **Generate Outputs**
   - TSV: machine-readable tabular format
   - Markdown: human-readable formatted preprint list
   - Visualizations: publication timeline, category distribution, keyword frequency
   - Author rankings for identifying key researchers

## Key Execution Patterns

**Track latest preprints in a topic:**
```bash
python preprint_tracker.py --query "CRISPR cancer" --max-results 100 --outdir ./latest
```

**Monitor bioRxiv only for computational methods:**
```bash
python preprint_tracker.py --query "deep learning genomics" --server biorxiv \
  --category bioinformatics --outdir ./compbio
```

**Track medical preprints from last week:**
```bash
python preprint_tracker.py --query "treatment COVID-19" --server medrxiv \
  --date-from 2024-03-01 --date-to 2024-03-07 --outdir ./weekly
```

**Find preprints from a specific timeframe:**
```bash
python preprint_tracker.py --query "gene therapy" \
  --date-from 2024-01-01 --date-to 2024-03-31 --outdir ./q1_2024
```

**Export as TSV for spreadsheet analysis:**
```bash
python preprint_tracker.py --query "immunotherapy cancer" --output-format tsv \
  --max-results 200 --outdir ./analysis
```

## Parameter Decision Guide

| Parameter | When to Use | Example |
|-----------|------------|---------|
| `--query` | Always required; use specific keywords for best results | `--query "CRISPR off-target effects"` |
| `--server biorxiv` | For life sciences and biology topics | `--server biorxiv` |
| `--server medrxiv` | For clinical, medical, and health topics | `--server medrxiv` |
| `--server both` | When topic spans biomedical research | `--server both` (default) |
| `--date-from` | Narrow timeframe for emerging topics | `--date-from 2024-02-01` for last month |
| `--max-results` | Use 50-100 for quick digest, 200+ for comprehensive analysis | `--max-results 200` |
| `--category` | Filter to specific bioRxiv category (genomics, bioinformatics, cancer biology, etc.) | `--category genomics` |
| `--sort-by date` | Get newest preprints first (default) | `--sort-by date` |
| `--output-format markdown` | For reading and sharing digests | `--output-format markdown` |
| `--output-format tsv` | For data analysis or spreadsheet import | `--output-format tsv` |

## Failure Modes & Recovery

| Error | Cause | Solution |
|-------|-------|----------|
| "No preprints found" | Query too specific or keywords too rare | Broaden query terms, try single keywords |
| "API timeout" | bioRxiv/medRxiv server slow or unreachable | Retry in 1-2 minutes, reduce --max-results |
| "Invalid date format" | Dates not in YYYY-MM-DD format | Use correct format, e.g., 2024-03-15 |
| "Category not recognized" | Invalid bioRxiv category string | Use valid category (genomics, bioinformatics, etc.) |
| "Memory error on large results" | Fetching too many preprints at once | Reduce --max-results or narrow date range |
| "Empty abstracts" | Some preprints lack full abstract in API | Expected; keywords still match titles |
| "Duplicate preprints" | Same preprint revised multiple times on server | Expected; manually filter if needed |

## Special Notes

**bioRxiv categories** include: genomics, bioinformatics, cancer biology, cell biology, computational biology, synthetic biology, systems biology, evolutionary biology, genetics, immunology, microbiology, neuroscience, molecular biology, plant biology, and many others.

**medRxiv categories** include: epidemiology, clinical practice, infectious diseases, oncology, cardiology, psychiatry, and others.

**API Rate Limits:** bioRxiv/medRxiv APIs are generally unrestricted, but large batch requests (>500 results) may see delays.

## Requirements
See `requirements.txt` for exact versions. Core dependencies:
- `requests` — HTTP API calls to bioRxiv/medRxiv
- `pandas` — Data organization and TSV export
- `matplotlib` — Visualization of trends and distributions
- `json` — Metadata serialization (stdlib)
- `re` — Regular expressions for text processing (stdlib)
- `datetime` — Date handling (stdlib)
