---
name: track-lab-publications
description: Track recent publications from a specific research lab or PI. Given a PI name and/or institution, fetches their recent papers from PubMed and Europe PMC, summarizes research themes, citation impact, and generates a formatted publication report.
---

# Track Lab Publications

## Purpose
Monitor and analyze a Principal Investigator's (PI) recent publication output, research themes, and academic impact. This skill automatically fetches publications from PubMed, extracts metadata, analyzes research keywords and trends, and generates comprehensive publication reports with citation metrics and collaboration networks.

## Use when / Not when

### Use when
- You need to **monitor a PI's recent publication activity** (within last 1-5 years)
- You want to **identify research themes and topics** a lab focuses on
- You need to **quantify citation impact** and journal diversity for a research group
- You're building a **collaborator network analysis** or funding proposal
- You want to **generate a formatted report** for grant applications or lab websites
- You need to **identify co-author patterns** and collaboration networks

### Not when
- You need **real-time trending topics** (PubMed data has 2-4 week lag)
- You require **proprietary databases** (preprints, conference abstracts)
- You need **full-text PDF analysis** (this skill uses only metadata)
- The PI has **very common name** without disambiguation (institution required)
- You need **complete citation graphs** (Europe PMC has limited historical coverage)

## Expected inputs / outputs

### Inputs
- **CSV/TSV** (optional): Pre-compiled PI and institution list
- **Command-line arguments**: PI name, institution, year range, result limits
- **API calls**: PubMed eUtils (free, no key required) and Europe PMC (free)

### Outputs
- **publications.tsv**: Tab-separated table with PMID, title, authors, journal, year, DOI, citation count, abstract
- **lab_report.md**: Formatted markdown report with PI summary, yearly trends, top journals, top keywords, collaborators
- **publications_timeline.png**: Bar chart of publication count per year
- **keyword_distribution.png**: Horizontal bar chart of top 15 keywords/MeSH terms
- **coauthor_network.png** (optional): Network graph of co-author collaborations
- **Summary to console**: PI name, total papers found, date range, top 5 journals

## Procedure

### Step 1: Query PubMed for publications
1. Construct search query: `<PI_name>[Author] AND <institution>[Affiliation]`
2. Call PubMed eSearch API with `retmax` (default 100), `mindate` and `maxdate` filters
3. Retrieve list of PMIDs matching criteria
4. If no results, retry with relaxed query (PI name only) with warning

### Step 2: Fetch publication details
1. Batch PMIDs into groups of 10 for eFetch API calls
2. Parse XML response to extract:
   - PubMed ID (PMID)
   - Title, authors (last name + initials), author affiliations
   - Journal name and impact factor (if available)
   - Publication year, DOI
   - Abstract (truncate to 500 chars for report)
   - MeSH descriptors (primary research keywords)
3. Handle missing fields gracefully (e.g., abstract may not exist)

### Step 3: Fetch citation counts (optional)
1. For each PMID, query Europe PMC API: `/search?query=PMID:<pmid>&format=json`
2. Extract `citedByCount` field
3. Handle rate limiting (max 25 requests/second)

### Step 4: Analyze research themes
1. **Keyword extraction**:
   - Combine all MeSH terms from all papers
   - Add high-frequency words from abstracts (exclude stopwords)
   - Rank by frequency across entire publication set
2. **Yearly trends**:
   - Count papers per year
   - Calculate growth rate
3. **Journal analysis**:
   - Identify top 5-10 journals by paper count
   - Track journal diversity (Shannon entropy)

### Step 5: Build co-author network (optional)
1. Extract all author names and affiliations from papers
2. Count author co-occurrence frequencies
3. Identify top collaborators (co-authors with 2+ joint papers)
4. Rank by collaboration frequency

### Step 6: Generate outputs
1. Write TSV table of all publications
2. Generate markdown report with:
   - PI profile header
   - Summary statistics
   - Year-by-year breakdown
   - Research theme summary
   - Top journals and keywords
   - Top collaborators
   - Full publication list with abstracts
3. Create matplotlib visualizations:
   - Time-series bar chart
   - Keyword frequency distribution
   - Co-author network graph (if requested)

## Key execution patterns

```bash
# Basic: Track a single PI by name and institution
python scripts/track_lab_publications.py \
  --pi-name "Jennifer Doudna" \
  --institution "UC Berkeley" \
  --years-back 5 \
  --outdir ./doudna_publications

# Track PI with citation counts and save as JSON
python scripts/track_lab_publications.py \
  --pi-name "Feng Zhang" \
  --institution "Broad Institute" \
  --fetch-citations \
  --output-format json \
  --outdir ./zhang_pubs

# Build co-author network for a PI
python scripts/track_lab_publications.py \
  --pi-name "Edith Heard" \
  --institution "Institut Curie" \
  --co-author-network \
  --max-results 200 \
  --outdir ./heard_collab

# Track publications and limit to 50 papers
python scripts/track_lab_publications.py \
  --pi-name "David Baker" \
  --years-back 3 \
  --max-results 50 \
  --outdir ./baker_recent

# Relaxed search (PI name only, no institution)
python scripts/track_lab_publications.py \
  --pi-name "Doudna" \
  --max-results 100
```

## Parameter decision guide

| Parameter | Value | When to use | Rationale |
|-----------|-------|-------------|-----------|
| `--pi-name` | Last name only (e.g., "Doudna") | PI has very unique surname | Broader search, may include unrelated authors |
| `--pi-name` | Full name (e.g., "Jennifer Doudna") | PI has common last name | Better specificity, fewer false positives |
| `--institution` | "UC Berkeley" or "UC San Francisco" | PI may work at multiple institutions | Disambiguate authors with same name |
| `--years-back` | 3 (default) | Current publication trend | Recent impact and research direction |
| `--years-back` | 5-10 | Track research evolution over time | Historical perspective on research trajectory |
| `--max-results` | 100 (default) | Most research groups | Captures 95% of active researcher output |
| `--max-results` | 200+ | Very prolific PIs (>20 papers/year) | Ensures complete publication list |
| `--fetch-citations` | True | Grant proposals or impact assessment | Demonstrates research influence |
| `--fetch-citations` | False | Quick turnaround needed | Saves 30-60 seconds API time |
| `--co-author-network` | True | Collaboration analysis or network studies | Visualize collaboration ecosystem |
| `--co-author-network` | False | Publication list only | Simpler output, faster execution |
| `--output-format` | markdown (default) | Human review or lab website | Best readability |
| `--output-format` | json | Data pipeline or downstream processing | Machine-readable, structured |
| `--output-format` | tsv | Import to spreadsheet or analysis | Excel/R-compatible format |

## Failure modes

| Failure | Cause | Solution |
|---------|-------|----------|
| **No papers found** | PI name misspelled or very common | Use full name + institution; check PubMed directly |
| **Too many results (1000+)** | PI name is common, institution too broad | Add institution, narrow years-back, or use initials |
| **Missing abstracts** | Older papers or non-indexed journals | Script handles gracefully; continues with available data |
| **API rate limit exceeded** | Too many citation fetch requests | Reduce `--max-results` or use `--fetch-citations=false` |
| **Europe PMC timeout** | Network issue or service unavailable | Retry; continue without citation counts |
| **Garbled author affiliations** | XML encoding issues from PubMed | Script sanitizes UTF-8; may lose diacritics |
| **Empty co-author network** | Few papers or no co-authors | Script skips network visualization if <5 papers |
| **Network graph too dense** | Many collaborators with many papers | Filter to top-50 edges; use `--co-author-network` sparingly |

