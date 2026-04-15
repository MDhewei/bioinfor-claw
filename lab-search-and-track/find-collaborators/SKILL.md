---
name: find-collaborators
description: Identify potential research collaborators based on topic overlap. Searches PubMed for researchers publishing on given topics, ranks them by publication output and topic relevance, and generates a ranked collaborator report with their research profiles.
---

# Find Collaborators

## Purpose
Discover and rank potential research collaborators based on specific research topics or keyword combinations. This skill systematically searches PubMed for active researchers in your research area, analyzes their publication profiles, identifies research overlap, and generates a ranked collaborator report with contact information and research summaries.

## Use when / Not when

### Use when
- You want to **find collaborators in a specific research area** (e.g., CRISPR, cancer genomics)
- You're **expanding your research network** within a specific field
- You need to **identify co-authors for grant proposals** (NIH, NSF, etc.)
- You're **analyzing the competitive landscape** for your research topic
- You want to **identify industry partners** or **academic collaborators** in your region/country
- You need to **discover rising stars** in your field (early-career researchers with high output)

### Not when
- You need **real-time trending topics** (PubMed lag: 2-4 weeks)
- You require **proprietary databases** (MEDLINE-only publications included)
- You need **researcher contact information** (affiliations only; manual lookup required)
- Topics are **highly specialized** with <10 active researchers (search may return too few)
- You need **patent analysis** or **non-publication research** tracking
- Collaborators are **outside life sciences** (PubMed does not index computer science, engineering, etc.)

## Expected inputs / outputs

### Inputs
- **Topic list** (comma-separated): "CRISPR base editing,cancer genomics,personalized medicine"
- **Command-line arguments**: year range, result limits, species/country filters
- **API calls**: PubMed eSearch/eFetch (free, no key required)

### Outputs
- **collaborators.tsv**: Tab-separated table with author name, affiliation, n_papers, topics_covered, top_journals, rank_score
- **collaborator_profiles.md**: Formatted markdown with one section per top-20 collaborator showing their full publication list and research summary
- **topic_author_heatmap.png**: Heatmap with authors on Y-axis, topics on X-axis, cell values = number of papers per author per topic
- **ranking_chart.png**: Horizontal bar chart of top-20 collaborators sorted by rank_score
- **Summary to console**: Total authors found, topic coverage stats, top 5 collaborators

## Procedure

### Step 1: Search for papers per topic
1. For each topic: construct PubMed query (default: topic as text search)
2. Call eSearch with year filters and max results per topic
3. Collect all PMIDs, deduplicate across topics
4. Track which topics each PMID maps to

### Step 2: Fetch publication details
1. Batch eFetch PMID records (10 per request)
2. Parse XML for:
   - All author names and affiliations
   - Publication year, journal, title
   - Journal impact factor (if available from XML)
3. Extract author last name + initials for deduplication

### Step 3: Build author-topic matrix
1. For each author: count papers per topic
2. Identify authors with presence in multiple topics
3. Filter by `--min-papers` threshold (default: 3 papers minimum)

### Step 4: Rank collaborators
1. Compute collaboration score per author:
   ```
   score = (num_topics_covered * 10) + log1p(total_papers)
   ```
2. Apply country/institution filters if specified
3. Sort by score, return top N

### Step 5: Augment author profiles (optional)
1. For top-30 authors: fetch their complete recent publication list
2. Extract research keywords (MeSH terms, abstracts)
3. Identify co-authors and collaboration frequency

### Step 6: Generate outputs
1. Write TSV table of all qualifying authors
2. Create markdown report with top-20 profiles
3. Generate visualizations:
   - Author-topic heatmap
   - Ranking bar chart
   - (Optional) research keyword cloud per author

## Key execution patterns

```bash
# Basic: Find collaborators in a single research area
python scripts/find_collaborators.py \
  --topics "CRISPR gene editing" \
  --years-back 5 \
  --outdir ./collaborators

# Multiple topics with overlap analysis
python scripts/find_collaborators.py \
  --topics "CRISPR base editing,epigenetics,precision medicine" \
  --years-back 3 \
  --max-per-topic 150 \
  --min-papers 5 \
  --outdir ./multi_topic_collab

# Find collaborators in a specific country
python scripts/find_collaborators.py \
  --topics "immunotherapy,cancer genomics" \
  --country-filter "Australia" \
  --years-back 2 \
  --outdir ./au_collab

# Find collaborators in human-focused research
python scripts/find_collaborators.py \
  --topics "CRISPR,genome editing" \
  --species-focus human \
  --exclude-institution "model organism" \
  --years-back 5 \
  --outdir ./human_crispr

# Large-scale collaborator discovery
python scripts/find_collaborators.py \
  --topics "machine learning,drug discovery,biotech" \
  --years-back 4 \
  --max-per-topic 300 \
  --min-papers 3 \
  --outdir ./ml_drug_discovery
```

## Parameter decision guide

| Parameter | Value | When to use | Rationale |
|-----------|-------|-------------|-----------|
| `--topics` | Single broad topic | Initial scoping | Fast search, large result set |
| `--topics` | 2-3 specific topics | Medium-specificity search | Good balance of precision and recall |
| `--topics` | 4+ topics | Highly interdisciplinary work | More comprehensive coverage |
| `--years-back` | 2 | Recent trends only | Captures active researchers in hot topics |
| `--years-back` | 3-5 (default) | Standard collaboration search | Captures established and emerging researchers |
| `--years-back` | 10+ | Established researcher tracking | Includes tenured researchers, may miss new entrants |
| `--max-per-topic` | 50-100 | Focused research areas | Fast search, top researchers by publication count |
| `--max-per-topic` | 200+ | Broad or popular topics | Comprehensive coverage, longer runtime |
| `--min-papers` | 3 (default) | Standard threshold | Most researchers meet this; ~50-100 results typical |
| `--min-papers` | 5-10 | High-output researchers only | Focus on prolific authors; ~5-50 results typical |
| `--species-focus` | human | Clinical or translational work | Filters out pure model organism research |
| `--species-focus` | mouse | Basic research, model organism labs | Filters out human-centric research |
| `--country-filter` | Set (e.g., "USA") | Regional network building | Enables local collaboration |
| `--exclude-institution` | "text" | Remove noise | E.g., exclude "company name" or research consortia |

## Failure modes

| Failure | Cause | Solution |
|---------|-------|----------|
| **Too many results (1000+)** | Topic too broad or common | Add year filters, increase `--min-papers`, split topics |
| **Too few results (<10)** | Topic too specific or niche | Reduce `--years-back`, lower `--min-papers` to 1-2 |
| **Duplicate authors** | Author name ambiguity (e.g., "Smith J") | Results include multiple people; manual review needed |
| **Poor affiliation parsing** | Garbled or non-standard affiliation strings | Extract institution manually from top results |
| **Missing recent researchers** | Lag in PubMed indexing (2-4 weeks) | Check preprints (bioRxiv) separately |
| **API rate limit exceeded** | Too many eFetch calls for large author set | Reduce `--max-per-topic` or increase delays |
| **Empty heatmap** | No author publishes in multiple topics | Topics are too specialized/non-overlapping |
| **Heatmap too dense** | Many authors × many topics | Filter to top-50 authors only |
| **Missing keywords** | MeSH terms incomplete or abstract not in XML | Script gracefully continues with available data |

