---
name: paper-digest-single
description: Digest a single scientific paper into a structured summary. Accepts a PubMed ID, DOI, arXiv ID, or local PDF path. Fetches abstract and metadata from PubMed/CrossRef APIs, extracts key findings, methods, results, and clinical implications, and saves a structured markdown report.
---

# Paper Digest Single

## Purpose
Extract and synthesize key information from a single scientific paper (PubMed, arXiv, CrossRef, or local PDF) into a structured, human-readable markdown report. Ideal for rapid literature review, evidence synthesis, and knowledge base building.

## When to Use / When NOT to Use

**Use when:**
- You need to quickly understand a single paper's main contributions
- Building a personal knowledge base from specific papers
- Preparing literature review sections that cite specific findings
- Evaluating clinical or research evidence for a specific topic
- You have a PubMed ID, DOI, arXiv ID, or PDF file ready to process

**Do NOT use when:**
- You need to compare multiple papers at once (use `pubmed-search` or `preprint-tracker` instead)
- You want to search for papers on a topic (use `pubmed-search` first)
- You need full-text PDF analysis beyond abstract-level summaries
- Processing confidential or restricted-access papers

## Expected Inputs & Outputs

**Inputs (choose ONE):**
- `--pmid INT` — PubMed ID (e.g., 35769483)
- `--doi STRING` — Digital Object Identifier (e.g., 10.1038/s41586-022-04826-7)
- `--arxiv STRING` — arXiv ID (e.g., 2301.12345)
- `--pdf PATH` — Local PDF file path

**Outputs (in --outdir):**
- `paper_digest_<identifier>.md` — Structured markdown report
- `paper_metadata.json` — Complete metadata (authors, MeSH, keywords, journal, etc.)
- `extraction_log.txt` — Processing notes and warnings

## Procedure

1. **Fetch Metadata**
   - PubMed (via eFetch XML API) or CrossRef API for DOI
   - Extract: title, authors, journal, year, abstract, MeSH terms, publication type

2. **Structured Text Extraction** (rule-based NLP)
   - Background: opening sentences (study motivation, context)
   - Methods: sentences with method keywords (RNA-seq, CRISPR, statistical tests, etc.)
   - Results: sentences containing quantitative findings (p-values, fold changes, percentages)
   - Conclusions: closing statements or sentences with conclusion keywords

3. **Entity Recognition**
   - Gene/protein names (HGNC symbols, patterns like ALL_CAPS or CamelCase ≥3 chars)
   - Drug names and treatment modalities
   - Disease/condition names

4. **Statistics Extraction**
   - Extract p-values, confidence intervals, fold changes, effect sizes
   - Identify sample sizes (n=, N=)

5. **Impact Assessment**
   - Classify journal impact (high/medium/low based on known high-impact journals)
   - Identify publication type (original research, review, meta-analysis, clinical trial)

6. **Generate Report**
   - Write markdown with sections: Citation, Abstract, Key Findings, Methods, Results, Conclusions, Entities, Keywords, Impact Notes

## Key Execution Patterns

**Single paper digest:**
```bash
python paper_digest_single.py --pmid 35769483 --outdir ./results
```

**DOI with citations:**
```bash
python paper_digest_single.py --doi 10.1038/s41586-022-04826-7 --include-citations 10 --outdir ./digest
```

**arXiv preprint:**
```bash
python paper_digest_single.py --arxiv 2301.12345 --outdir ./arxiv_digest
```

**Local PDF (metadata from filename):**
```bash
python paper_digest_single.py --pdf ./papers/my_study.pdf --outdir ./results
```

**JSON output for programmatic use:**
```bash
python paper_digest_single.py --pmid 35769483 --output-format json --outdir ./results
```

**Topic-specific extraction:**
```bash
python paper_digest_single.py --pmid 35769483 --topic-keywords "CRISPR,cancer,immunotherapy" --outdir ./results
```

## Parameter Decision Guide

| Parameter | When to Use | Example |
|-----------|------------|---------|
| `--pmid` | When you have a PubMed ID from pubmed.gov | `--pmid 35769483` |
| `--doi` | When citing a paper by DOI (most journals support this) | `--doi 10.1038/nature12373` |
| `--arxiv` | For preprints before peer review publication | `--arxiv 2301.12345` |
| `--pdf` | When you have a local copy and no ID available | `--pdf ./paper.pdf` |
| `--include-citations` | To include top cited references in output | `--include-citations 5` for top 5 |
| `--topic-keywords` | To highlight domain-specific findings | `--topic-keywords "CRISPR,off-target"` |
| `--output-format` | Use `json` for pipeline integration, `markdown` for reading | `--output-format json` |

## Failure Modes & Recovery

| Error | Cause | Solution |
|-------|-------|----------|
| "PMID not found" | Invalid PubMed ID or network issue | Verify PMID at pubmed.gov, check internet connection |
| "DOI resolution failed" | Invalid DOI format or CrossRef API timeout | Verify DOI format (10.xxxx/...), try using PMID instead |
| "arXiv ID not found" | Invalid ID format or not an arXiv paper | Verify ID at arxiv.org |
| "PDF parsing incomplete" | Complex PDF structure or scanned image | Try converting PDF to text first, or use PMID/DOI |
| "No abstract found" | Paper predates abstract requirement or restricted access | Use locally available metadata or title/keywords only |
| "Rate limit exceeded" | Too many API requests in short time | Wait 30 seconds before retrying |
| "Network timeout" | API server unreachable | Retry in 1-2 minutes, check network connectivity |

## Requirements
See `requirements.txt` for exact versions. Core dependencies:
- `requests` — HTTP API calls to PubMed, CrossRef, arXiv
- `pandas` — Data organization and export
- `xml.etree.ElementTree` — PubMed XML parsing (stdlib)
- `json` — Metadata serialization (stdlib)
- `re` — Regular expression entity extraction (stdlib)
