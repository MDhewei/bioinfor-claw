#!/usr/bin/env python3
"""
pubmed_search.py: Search PubMed for papers and generate formatted reports
Retrieves ranked results with abstracts, citations, and metadata
Generates TSV, markdown, and visualization outputs
"""

import argparse
import json
import os
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from xml.etree import ElementTree as ET
from collections import Counter

import requests
import pandas as pd
import matplotlib.pyplot as plt


class PubMedSearcher:
    """Search PubMed and format results."""

    # Common English stopwords
    STOPWORDS = {
        'the', 'a', 'an', 'and', 'or', 'but', 'in', 'on', 'at', 'to', 'for',
        'is', 'was', 'are', 'were', 'be', 'been', 'being', 'have', 'has', 'had',
        'do', 'does', 'did', 'will', 'would', 'could', 'should', 'may', 'might',
        'can', 'this', 'that', 'these', 'those', 'i', 'you', 'he', 'she', 'it',
        'we', 'they', 'what', 'which', 'who', 'when', 'where', 'why', 'how',
        'all', 'each', 'every', 'both', 'few', 'more', 'most', 'other', 'some',
        'such', 'no', 'nor', 'not', 'only', 'own', 'same', 'so', 'than', 'too',
        'very', 'by', 'from', 'up', 'about', 'of', 'with', 'as', 'into', 'out'
    }

    def __init__(self, outdir: str = './pubmed_search'):
        self.outdir = Path(outdir)
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.log_messages = []
        self.papers = []

    def log(self, msg: str, level: str = 'INFO'):
        """Add message to processing log."""
        timestamp = datetime.now().isoformat(timespec='seconds')
        log_line = f"[{timestamp}] {level}: {msg}"
        self.log_messages.append(log_line)
        print(log_line)

    def search_pubmed(self, query: str, max_results: int, date_from: str,
                     date_to: str, sort_by: str = 'relevance') -> List[int]:
        """Execute PubMed eSearch to get list of PMIDs."""
        self.log(f"Searching PubMed with query: '{query}'")

        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        params = {
            'db': 'pubmed',
            'term': query,
            'retmode': 'json',
            'retmax': min(max_results, 10000),
            'mindate': date_from,
            'maxdate': date_to,
            'datetype': 'pdat',
            'sort': 'date' if sort_by == 'date' else 'relevance'
        }

        try:
            response = requests.get(url, params=params, timeout=15)
            response.raise_for_status()
            data = response.json()

            result_count = int(data['esearchresult'].get('count', 0))
            pmid_list = data['esearchresult'].get('idlist', [])

            self.log(f"Found {result_count} total results; retrieving {len(pmid_list)} PMIDs")
            return [int(pmid) for pmid in pmid_list]

        except Exception as e:
            self.log(f"Error searching PubMed: {e}", 'ERROR')
            return []

    def fetch_metadata_batch(self, pmids: List[int]) -> List[Dict]:
        """Fetch metadata for a batch of PMIDs using eFetch."""
        if not pmids:
            return []

        self.log(f"Fetching metadata for {len(pmids)} papers...")
        papers = []

        # Batch fetch (up to 200 per request)
        batch_size = 200
        for i in range(0, len(pmids), batch_size):
            batch = pmids[i:i+batch_size]
            pmid_str = ','.join(str(p) for p in batch)

            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            params = {
                'db': 'pubmed',
                'id': pmid_str,
                'retmode': 'xml'
            }

            try:
                response = requests.get(url, params=params, timeout=30)
                response.raise_for_status()

                root = ET.fromstring(response.content)
                for article in root.findall('.//PubmedArticle'):
                    paper_data = self._parse_article(article)
                    if paper_data:
                        papers.append(paper_data)

                self.log(f"Fetched batch {i//batch_size + 1}: {len(papers)} papers so far")

            except Exception as e:
                self.log(f"Error fetching batch starting at index {i}: {e}", 'WARN')

        self.log(f"Fetched total of {len(papers)} papers with complete metadata")
        return papers

    def _parse_article(self, article: ET.Element) -> Optional[Dict]:
        """Parse a single PubmedArticle XML element."""
        paper = {
            'pmid': '',
            'title': '',
            'authors': [],
            'journal': '',
            'year': '',
            'volume': '',
            'pages': '',
            'doi': '',
            'abstract': '',
            'abstract_preview': '',
            'publication_types': [],
            'mesh_terms': [],
            'keywords': [],
            'cited_by_count': 0
        }

        # PMID
        pmid_elem = article.find(".//PMID[@Version='1']")
        if pmid_elem is not None and pmid_elem.text:
            paper['pmid'] = pmid_elem.text

        # Title
        title_elem = article.find('.//ArticleTitle')
        if title_elem is not None and title_elem.text:
            paper['title'] = title_elem.text.strip()

        # Authors
        author_list = article.find('.//AuthorList')
        if author_list is not None:
            for author in author_list.findall('Author'):
                last_name = author.findtext('LastName', '')
                initials = author.findtext('Initials', '')
                if last_name:
                    paper['authors'].append(f"{last_name} {initials}".strip())

        # Journal
        journal_elem = article.find('.//Journal/Title')
        if journal_elem is not None and journal_elem.text:
            paper['journal'] = journal_elem.text.strip()

        # Publication date
        pub_date = article.find('.//PubDate')
        if pub_date is not None:
            year = pub_date.findtext('Year', '')
            if year:
                paper['year'] = year

        # Volume and pages
        volume = article.find('.//Volume')
        if volume is not None and volume.text:
            paper['volume'] = volume.text

        pages = article.find('.//MedlinePgn')
        if pages is not None and pages.text:
            paper['pages'] = pages.text

        # DOI
        doi_elem = article.find(".//ArticleId[@IdType='doi']")
        if doi_elem is not None and doi_elem.text:
            paper['doi'] = doi_elem.text.strip()

        # Abstract
        abstract_elem = article.find('.//Abstract')
        if abstract_elem is not None:
            abstract_parts = []
            for section in abstract_elem.findall('.//AbstractText'):
                if section.text:
                    abstract_parts.append(section.text.strip())
            paper['abstract'] = ' '.join(abstract_parts)
            # Create preview (first 100 words)
            words = paper['abstract'].split()[:100]
            paper['abstract_preview'] = ' '.join(words)
            if len(paper['abstract'].split()) > 100:
                paper['abstract_preview'] += '...'

        # Publication types
        pub_type_list = article.find('.//PublicationTypeList')
        if pub_type_list is not None:
            for pub_type in pub_type_list.findall('PublicationType'):
                if pub_type.text:
                    paper['publication_types'].append(pub_type.text.strip())

        # MeSH terms
        mesh_list = article.find('.//MeshHeadingList')
        if mesh_list is not None:
            for mesh in mesh_list.findall('MeshHeading'):
                descriptor = mesh.find('DescriptorName')
                if descriptor is not None and descriptor.text:
                    paper['mesh_terms'].append(descriptor.text.strip())

        # Keywords
        keyword_list = article.find('.//KeywordList')
        if keyword_list is not None:
            for keyword in keyword_list.findall('Keyword'):
                if keyword.text:
                    paper['keywords'].append(keyword.text.strip())

        return paper if paper['pmid'] else None

    def fetch_citations(self, pmids: List[int]) -> Dict[int, int]:
        """Fetch citation counts from Europe PMC API."""
        self.log(f"Fetching citation counts for {len(pmids)} papers from Europe PMC...")
        citation_map = {}

        for pmid in pmids[:100]:  # Limit to avoid rate limiting
            try:
                url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
                params = {
                    'query': f'EXT_ID:{pmid} AND SRC:MED',
                    'resulttype': 'core',
                    'format': 'json'
                }
                response = requests.get(url, params=params, timeout=10)
                response.raise_for_status()
                data = response.json()

                if 'resultList' in data and data['resultList']['result']:
                    result = data['resultList']['result'][0]
                    citation_map[pmid] = result.get('citedByCount', 0)
                else:
                    citation_map[pmid] = 0

            except Exception as e:
                self.log(f"Warning: Could not fetch citations for PMID {pmid}", 'WARN')
                citation_map[pmid] = 0

        return citation_map

    def extract_keywords(self, papers: List[Dict], top_n: int = 20) -> List[Tuple[str, int]]:
        """Extract and rank keywords from all abstracts using TF-IDF."""
        all_text = ' '.join([p['abstract'] for p in papers])

        # Tokenize and clean
        words = re.findall(r'\b[a-z]+\b', all_text.lower())

        # Filter stopwords and short words
        filtered_words = [w for w in words if w not in self.STOPWORDS and len(w) > 3]

        # Count frequency
        word_freq = Counter(filtered_words)

        # Return top N
        return word_freq.most_common(top_n)

    def filter_papers(self, papers: List[Dict], pub_type: Optional[str] = None,
                     journal_filter: Optional[str] = None) -> List[Dict]:
        """Filter papers by publication type and journal."""
        filtered = papers

        # Filter by publication type
        if pub_type:
            type_map = {
                'review': ['Review'],
                'clinical_trial': ['Clinical Trial'],
                'meta_analysis': ['Meta-Analysis'],
                'original': ['Journal Article']
            }
            allowed_types = type_map.get(pub_type, [])
            filtered = [p for p in filtered if any(t in p['publication_types'] for t in allowed_types)]
            self.log(f"Filtered to {len(filtered)} papers of type '{pub_type}'")

        # Filter by journal
        if journal_filter:
            journals = [j.strip().lower() for j in journal_filter.split(',')]
            filtered = [p for p in filtered if any(j in p['journal'].lower() for j in journals)]
            self.log(f"Filtered to {len(filtered)} papers from specified journals")

        return filtered

    def sort_papers(self, papers: List[Dict], sort_by: str = 'relevance',
                   citation_map: Optional[Dict] = None) -> List[Dict]:
        """Sort papers by relevance, date, or citations."""
        if sort_by == 'date':
            papers_sorted = sorted(papers, key=lambda x: x.get('year', '0'), reverse=True)
        elif sort_by == 'citations':
            if citation_map:
                papers_sorted = sorted(papers,
                                      key=lambda x: citation_map.get(int(x['pmid']), 0),
                                      reverse=True)
            else:
                papers_sorted = papers
        else:  # relevance (maintain original order from eSearch)
            papers_sorted = papers

        return papers_sorted

    def generate_markdown_report(self, papers: List[Dict]) -> str:
        """Generate formatted markdown report."""
        md = f"""# PubMed Search Results

**Query Date:** {datetime.now().isoformat(timespec='seconds')}
**Total Results:** {len(papers)}

## Papers

"""
        for i, paper in enumerate(papers, 1):
            authors = ', '.join(paper['authors'][:3])
            if len(paper['authors']) > 3:
                authors += f" + {len(paper['authors']) - 3} more"

            md += f"""### {i}. {paper['title']}

- **PMID:** {paper['pmid']}
- **Authors:** {authors}
- **Journal:** {paper['journal']} ({paper['year']})
- **DOI:** {paper['doi'] if paper['doi'] else 'Not available'}
- **Publication Type:** {', '.join(paper['publication_types'])}

**Abstract:** {paper['abstract_preview']}

---

"""
        return md

    def generate_visualizations(self, papers: List[Dict], keywords: List[Tuple[str, int]]):
        """Generate publication trend and keyword visualizations."""
        # Timeline: papers per year
        year_counts = {}
        for paper in papers:
            year = paper['year']
            if year:
                year_counts[year] = year_counts.get(year, 0) + 1

        if year_counts:
            years = sorted(year_counts.keys())
            counts = [year_counts[y] for y in years]

            plt.figure(figsize=(12, 5))
            plt.bar(years, counts, color='steelblue', alpha=0.7)
            plt.xlabel('Year', fontsize=12)
            plt.ylabel('Number of Papers', fontsize=12)
            plt.title('Publication Timeline', fontsize=14, fontweight='bold')
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(self.outdir / 'publication_timeline.png', dpi=100)
            plt.close()
            self.log("Saved publication_timeline.png")

        # Keyword cloud: top keywords
        if keywords:
            kwords = [k[0] for k in keywords]
            kfreqs = [k[1] for k in keywords]

            plt.figure(figsize=(12, 6))
            plt.barh(kwords, kfreqs, color='coral', alpha=0.7)
            plt.xlabel('Frequency', fontsize=12)
            plt.title('Top Keywords in Abstracts', fontsize=14, fontweight='bold')
            plt.tight_layout()
            plt.savefig(self.outdir / 'keyword_cloud.png', dpi=100)
            plt.close()
            self.log("Saved keyword_cloud.png")

        # Journal distribution: pie chart
        journal_counts = {}
        for paper in papers:
            j = paper['journal']
            journal_counts[j] = journal_counts.get(j, 0) + 1

        # Group small journals as "Other"
        top_journals = dict(sorted(journal_counts.items(), key=lambda x: x[1], reverse=True)[:10])
        other_count = sum(v for k, v in journal_counts.items() if k not in top_journals)
        if other_count > 0:
            top_journals['Other'] = other_count

        if top_journals:
            plt.figure(figsize=(10, 8))
            plt.pie(top_journals.values(), labels=top_journals.keys(), autopct='%1.1f%%',
                   startangle=90)
            plt.title('Journal Distribution', fontsize=14, fontweight='bold')
            plt.tight_layout()
            plt.savefig(self.outdir / 'journal_distribution.png', dpi=100)
            plt.close()
            self.log("Saved journal_distribution.png")

    def save_results(self, papers: List[Dict], keywords: List[Tuple[str, int]],
                    output_format: str = 'tsv'):
        """Save results to files."""
        # Prepare DataFrame
        df_data = []
        for paper in papers:
            df_data.append({
                'PMID': paper['pmid'],
                'Title': paper['title'],
                'Authors': '; '.join(paper['authors'][:5]),
                'Journal': paper['journal'],
                'Year': paper['year'],
                'DOI': paper['doi'],
                'Abstract_Preview': paper['abstract_preview'],
                'Publication_Types': '; '.join(paper['publication_types']),
                'MeSH_Terms': '; '.join(paper['mesh_terms'][:5]),
            })

        df = pd.DataFrame(df_data)

        # Save TSV
        tsv_file = self.outdir / 'search_results.tsv'
        df.to_csv(tsv_file, sep='\t', index=False)
        self.log(f"Saved results to {tsv_file}")

        # Save JSON if requested
        if output_format == 'json':
            json_file = self.outdir / 'search_results.json'
            with open(json_file, 'w') as f:
                json.dump(df_data, f, indent=2)
            self.log(f"Saved results to {json_file}")

        # Save markdown report
        md_report = self.generate_markdown_report(papers)
        md_file = self.outdir / 'search_report.md'
        with open(md_file, 'w') as f:
            f.write(md_report)
        self.log(f"Saved markdown report to {md_file}")

        # Save keywords
        keywords_file = self.outdir / 'trending_keywords.tsv'
        keywords_df = pd.DataFrame(keywords, columns=['Keyword', 'Frequency'])
        keywords_df.to_csv(keywords_file, sep='\t', index=False)
        self.log(f"Saved keywords to {keywords_file}")

    def save_log(self):
        """Save processing log."""
        log_file = self.outdir / 'search_log.txt'
        with open(log_file, 'w') as f:
            f.write('\n'.join(self.log_messages))
        self.log(f"Saved log to {log_file}")

    def search(self, query: str, max_results: int = 50, date_from: str = '2020/01/01',
              date_to: str = '2024/12/31', pub_type: Optional[str] = None,
              journal_filter: Optional[str] = None, sort_by: str = 'relevance',
              fetch_citations: bool = False, output_format: str = 'tsv',
              abstract_words: int = 100) -> bool:
        """Execute full search pipeline."""

        # Search PubMed
        pmids = self.search_pubmed(query, max_results, date_from, date_to, sort_by)
        if not pmids:
            self.log("No results found for query", 'ERROR')
            return False

        # Fetch metadata
        papers = self.fetch_metadata_batch(pmids)
        if not papers:
            self.log("Failed to fetch paper metadata", 'ERROR')
            return False

        # Filter papers
        papers = self.filter_papers(papers, pub_type, journal_filter)

        # Fetch citations if requested
        citation_map = {}
        if fetch_citations:
            citation_map = self.fetch_citations([int(p['pmid']) for p in papers])
            # Add to papers
            for paper in papers:
                paper['cited_by_count'] = citation_map.get(int(paper['pmid']), 0)

        # Sort papers
        papers = self.sort_papers(papers, sort_by, citation_map)

        # Extract keywords
        keywords = self.extract_keywords(papers)

        # Save results
        self.save_results(papers, keywords, output_format)

        # Generate visualizations
        self.generate_visualizations(papers, keywords)

        # Save log
        self.save_log()

        # Print summary
        self.log(f"\n=== SEARCH SUMMARY ===")
        self.log(f"Query: {query}")
        self.log(f"Date range: {date_from} to {date_to}")
        self.log(f"Total results: {len(pmids)}")
        self.log(f"Papers retrieved: {len(papers)}")
        self.log(f"Top 5 keywords: {', '.join([k[0] for k in keywords[:5]])}")

        return True


def main():
    parser = argparse.ArgumentParser(
        description='Search PubMed and generate formatted reports',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python pubmed_search.py --query "CRISPR gene therapy" --max-results 50
  python pubmed_search.py --query "cancer immunotherapy" --pub-type review
  python pubmed_search.py --query "CAR-T therapy" --sort-by citations --fetch-citations
  python pubmed_search.py --query "Alzheimer's[MeSH]" --date-from 2022/01/01
        """
    )

    parser.add_argument('--query', type=str, required=True, help='PubMed search query')
    parser.add_argument('--max-results', type=int, default=50, help='Max results to retrieve')
    parser.add_argument('--date-from', default='2020/01/01', help='Start date (YYYY or YYYY/MM/DD)')
    parser.add_argument('--date-to', default='2024/12/31', help='End date (YYYY or YYYY/MM/DD)')
    parser.add_argument('--pub-type', choices=['review', 'clinical_trial', 'meta_analysis', 'original'],
                       help='Filter by publication type')
    parser.add_argument('--journal-filter', help='Comma-separated journal names')
    parser.add_argument('--sort-by', choices=['relevance', 'date', 'citations'],
                       default='relevance', help='Sort results by...')
    parser.add_argument('--fetch-citations', action='store_true', help='Fetch citation counts')
    parser.add_argument('--outdir', default='./pubmed_search', help='Output directory')
    parser.add_argument('--output-format', choices=['tsv', 'json'], default='tsv')
    parser.add_argument('--abstract-words', type=int, default=100, help='Abstract preview word count')

    args = parser.parse_args()

    searcher = PubMedSearcher(args.outdir)
    success = searcher.search(
        query=args.query,
        max_results=args.max_results,
        date_from=args.date_from,
        date_to=args.date_to,
        pub_type=args.pub_type,
        journal_filter=args.journal_filter,
        sort_by=args.sort_by,
        fetch_citations=args.fetch_citations,
        output_format=args.output_format,
        abstract_words=args.abstract_words
    )

    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
