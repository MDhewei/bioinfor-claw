# Bioinformatics Paper Search & Digest Skills

Three complete, production-ready Python skills for scientific literature research and analysis.

## Overview

| Skill | Purpose | Use When |
|-------|---------|----------|
| **paper-digest-single** | Extract and structure a single paper into a markdown report | You have a specific PMID, DOI, arXiv ID, or PDF to analyze |
| **pubmed-search** | Search PubMed and retrieve ranked papers with metadata | Building a literature review or finding papers on a topic |
| **preprint-tracker** | Monitor recent preprints from bioRxiv/medRxiv | Tracking emerging research trends before peer review |

---

## Skill 1: paper-digest-single

**Location:** `paper-digest-single/`

### Purpose
Digest a single scientific paper into a structured summary with automatic extraction of key findings, methods, results, and clinical implications.

### Key Features
- Accepts: PubMed ID, DOI, arXiv ID, or local PDF file
- Fetches metadata from PubMed/CrossRef/arXiv APIs
- Rule-based NLP extraction of Background, Methods, Results, Conclusions
- Entity recognition: genes, proteins, drugs, diseases
- Statistical extraction: p-values, fold changes, sample sizes
- Journal impact assessment
- Multiple output formats: markdown, JSON, plain text

### Usage Examples

```bash
# Digest a PubMed paper
python paper-digest-single/scripts/paper_digest_single.py --pmid 35769483

# Digest by DOI with citation info
python paper-digest-single/scripts/paper_digest_single.py \
  --doi 10.1038/s41586-022-04826-7 \
  --include-citations 5

# Track topic-specific findings
python paper-digest-single/scripts/paper_digest_single.py \
  --pmid 35769483 \
  --topic-keywords "CRISPR,off-target,immune"

# Export as JSON for programmatic use
python paper-digest-single/scripts/paper_digest_single.py \
  --arxiv 2301.12345 \
  --output-format json
```

### Outputs
- `paper_digest_<ID>.md` — Structured markdown report
- `paper_metadata_<ID>.json` — Complete metadata
- `extraction_log_<ID>.txt` — Processing notes

### Requirements
```
requests>=2.28.0
pandas>=1.5.0
```

---

## Skill 2: pubmed-search

**Location:** `pubmed-search/`

### Purpose
Execute targeted searches of PubMed to find relevant papers, with rich metadata and visualizations of publication trends.

### Key Features
- Full PubMed search syntax support (MeSH terms, author names, date ranges)
- Batch metadata fetching (up to 200 papers per request)
- Optional citation count retrieval from Europe PMC
- Filtering by publication type (review, clinical trial, meta-analysis, original)
- Journal filtering and sorting
- Keyword extraction and frequency analysis
- Automated visualizations:
  - Publication timeline (papers per year)
  - Keyword cloud (top keywords)
  - Journal distribution (pie chart)
- TSV and JSON exports

### Usage Examples

```bash
# Simple search
python pubmed-search/scripts/pubmed_search.py \
  --query "CRISPR gene therapy" \
  --max-results 50

# Reviews only, from top journals
python pubmed-search/scripts/pubmed_search.py \
  --query "cancer immunotherapy" \
  --pub-type review \
  --journal-filter "Nature,Science,Lancet" \
  --max-results 100

# Most cited papers (slower, with Europe PMC API)
python pubmed-search/scripts/pubmed_search.py \
  --query "GWAS risk factors" \
  --sort-by citations \
  --fetch-citations \
  --max-results 100

# Custom date range with MeSH terms
python pubmed-search/scripts/pubmed_search.py \
  --query "Alzheimer's Disease[MeSH]" \
  --date-from 2022/01/01 \
  --date-to 2024/12/31
```

### Outputs
- `search_results.tsv` — Tabular results (machine-readable)
- `search_report.md` — Formatted markdown list
- `trending_keywords.tsv` — Top keywords with frequencies
- `publication_timeline.png` — Papers per year chart
- `keyword_cloud.png` — Keyword frequency bar chart
- `journal_distribution.png` — Papers per journal pie chart
- `search_log.txt` — Query parameters and execution notes

### Requirements
```
requests>=2.28.0
pandas>=1.5.0
matplotlib>=3.5.0
```

---

## Skill 3: preprint-tracker

**Location:** `preprint-tracker/`

### Purpose
Monitor emerging research trends by searching recent preprints from bioRxiv and medRxiv, identifying trending topics and prolific authors.

### Key Features
- Search bioRxiv (life sciences) and medRxiv (clinical/medical) independently or together
- API-based fetching of preprint metadata and abstracts
- Keyword trend analysis across all preprints
- Author prolific identification (ranking by preprint count)
- Category filtering (genomics, bioinformatics, cancer biology, etc.)
- 30-day default lookback (customizable)
- Automated visualizations:
  - Publication timeline (per week)
  - Keyword frequency (top 20)
  - Category distribution (bar chart)
- TSV and markdown exports

### Usage Examples

```bash
# Track latest preprints on a topic
python preprint-tracker/scripts/preprint_tracker.py \
  --query "CRISPR cancer" \
  --max-results 100

# bioRxiv only, computational biology
python preprint-tracker/scripts/preprint_tracker.py \
  --query "deep learning genomics" \
  --server biorxiv \
  --category bioinformatics

# Track medical preprints from specific week
python preprint-tracker/scripts/preprint_tracker.py \
  --query "COVID-19 treatment" \
  --server medrxiv \
  --date-from 2024-03-01 \
  --date-to 2024-03-07

# Comprehensive 3-month analysis
python preprint-tracker/scripts/preprint_tracker.py \
  --query "gene therapy" \
  --date-from 2024-01-01 \
  --date-to 2024-03-31 \
  --max-results 200
```

### Outputs
- `preprints.tsv` — Tabular preprint data
- `trending_topics.tsv` — Top keywords with frequencies
- `prolific_authors.tsv` — Authors ranked by preprint count
- `digest_report.md` — Formatted preprint summaries
- `timeline.png` — Weekly publication bar chart
- `topic_frequency.png` — Keyword frequency chart
- `category_distribution.png` — Category distribution bar chart
- `search_log.txt` — Processing notes and summary

### Requirements
```
requests>=2.28.0
pandas>=1.5.0
matplotlib>=3.5.0
```

---

## Installation & Setup

### 1. Install Python 3.8+
```bash
python3 --version  # Should be 3.8 or higher
```

### 2. Install skill dependencies
```bash
# For paper-digest-single
pip install -r paper-digest-single/requirements.txt

# For pubmed-search
pip install -r pubmed-search/requirements.txt

# For preprint-tracker
pip install -r preprint-tracker/requirements.txt

# Or install all at once
pip install requests pandas matplotlib
```

### 3. Verify installation
```bash
# Check help for each skill
python3 paper-digest-single/scripts/paper_digest_single.py --help
python3 pubmed-search/scripts/pubmed_search.py --help
python3 preprint-tracker/scripts/preprint_tracker.py --help
```

---

## Common Workflows

### Workflow 1: Quick Literature Review
```bash
# 1. Search for papers on your topic
python3 pubmed-search/scripts/pubmed_search.py \
  --query "your topic" \
  --max-results 50 \
  --pub-type review

# 2. Digest top papers individually
python3 paper-digest-single/scripts/paper_digest_single.py \
  --pmid 12345678
```

### Workflow 2: Track Emerging Field
```bash
# 1. Get latest preprints this week
python3 preprint-tracker/scripts/preprint_tracker.py \
  --query "your topic" \
  --date-from 2024-04-01 \
  --max-results 100

# 2. Deep dive into interesting preprints
python3 paper-digest-single/scripts/paper_digest_single.py \
  --arxiv 2301.12345
```

### Workflow 3: Comprehensive Evidence Synthesis
```bash
# 1. Search comprehensive database
python3 pubmed-search/scripts/pubmed_search.py \
  --query "disease[MeSH]" \
  --date-from 2021/01/01 \
  --max-results 200 \
  --sort-by citations \
  --fetch-citations

# 2. Export results for analysis
# TSV results are in search_results.tsv

# 3. Digest key papers
for pmid in $(cut -f1 pubmed_search/search_results.tsv | tail -10); do
  python3 paper-digest-single/scripts/paper_digest_single.py --pmid "$pmid"
done
```

---

## API Note

These skills use public APIs:
- **PubMed eSearch/eFetch**: No API key required; NCBI requests reasonable use (no automated scraping)
- **CrossRef**: No API key required
- **arXiv**: No API key required
- **bioRxiv/medRxiv**: No API key required
- **Europe PMC**: No API key required (optional citations feature)

All APIs have built-in rate limiting; scripts respect this automatically.

---

## Performance Notes

- **paper-digest-single**: ~2-5 seconds per paper (depends on API response time)
- **pubmed-search**: ~1-3 seconds per result batch; 30-60 seconds for 50 papers + visualizations
- **preprint-tracker**: ~2-5 seconds per batch; 50-100 papers typically 30-60 seconds

For large batch operations (>500 papers), consider running searches overnight.

---

## Troubleshooting

### "Network timeout" errors
- Check internet connectivity
- Retry after 1-2 minutes (API server may be temporarily busy)
- Reduce --max-results if timeouts persist

### "No results found"
- Try simpler keywords (remove MeSH terms, author names)
- Expand date range
- Check query syntax against pubmed.gov

### Import errors (requests, pandas, matplotlib)
- Install missing packages: `pip install requests pandas matplotlib`
- Check Python version: `python3 --version`

### API rate limits
- Scripts have built-in delays; just wait and retry
- Europe PMC citations feature is slower; skip with `--fetch-citations` omitted

---

## File Structure

```
paper-search-and-digest/
├── paper-digest-single/
│   ├── SKILL.md                          # Full documentation
│   ├── requirements.txt                  # Python dependencies
│   └── scripts/
│       └── paper_digest_single.py        # Main script
├── pubmed-search/
│   ├── SKILL.md                          # Full documentation
│   ├── requirements.txt                  # Python dependencies
│   └── scripts/
│       └── pubmed_search.py              # Main script
├── preprint-tracker/
│   ├── SKILL.md                          # Full documentation
│   ├── requirements.txt                  # Python dependencies
│   └── scripts/
│       └── preprint_tracker.py           # Main script
└── README.md                              # This file
```

---

## License & Attribution

These skills are designed for academic and research use. Always respect copyright and terms of service when using external APIs and databases.

---

**Last Updated:** April 9, 2026
