#!/usr/bin/env python3
"""
paper_digest_single.py: Extract and digest a single scientific paper
Accepts: PubMed ID, DOI, arXiv ID, or local PDF path
Outputs: Structured markdown report and JSON metadata
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

import requests
import pandas as pd


class PaperDigester:
    """Extract and structure metadata from a single scientific paper."""

    # High-impact journal patterns
    HIGH_IMPACT_JOURNALS = {
        'Nature', 'Science', 'Cell', 'PNAS', 'JAMA', 'Lancet',
        'New England Journal of Medicine', 'Nature Medicine', 'Nature Genetics',
        'Nature Biotechnology', 'Science Translational Medicine', 'Cell Reports',
        'eLife', 'EMBO Journal', 'Nature Communications'
    }

    # Gene/protein entity patterns
    GENE_PATTERN = re.compile(r'\b[A-Z]{2,}[A-Z0-9]*(?:\d)?(?:\-[A-Z0-9]+)?\b')
    PVALUE_PATTERN = re.compile(r'[pP]\s*[<>=]+\s*(?:0\.)?0*(\d+)', re.IGNORECASE)
    FOLD_CHANGE_PATTERN = re.compile(r'([\d.]+)[\s\-]*fold(?:\s+change)?', re.IGNORECASE)

    def __init__(self, outdir: str = './paper_digest'):
        self.outdir = Path(outdir)
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.log_messages = []

    def log(self, msg: str, level: str = 'INFO'):
        """Add message to processing log."""
        timestamp = datetime.now().isoformat(timespec='seconds')
        log_line = f"[{timestamp}] {level}: {msg}"
        self.log_messages.append(log_line)
        print(log_line)

    def fetch_pubmed(self, pmid: int) -> Dict:
        """Fetch paper metadata from PubMed via eFetch."""
        self.log(f"Fetching PubMed metadata for PMID {pmid}...")

        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        params = {
            'db': 'pubmed',
            'id': str(pmid),
            'retmode': 'xml'
        }

        try:
            response = requests.get(url, params=params, timeout=10)
            response.raise_for_status()
        except requests.RequestException as e:
            self.log(f"Failed to fetch PMID {pmid}: {e}", 'ERROR')
            return {}

        try:
            root = ET.fromstring(response.content)
            article = root.find('.//PubmedArticle')
            if article is None:
                self.log(f"No article found for PMID {pmid}", 'WARN')
                return {}

            metadata = self._parse_pubmed_xml(article, pmid)
            self.log(f"Successfully extracted metadata for: {metadata.get('title', 'Unknown')}")
            return metadata

        except Exception as e:
            self.log(f"Error parsing PubMed XML: {e}", 'ERROR')
            return {}

    def _parse_pubmed_xml(self, article: ET.Element, pmid: int) -> Dict:
        """Parse PubMed XML article element."""
        metadata = {
            'identifier': str(pmid),
            'identifier_type': 'pmid',
            'title': '',
            'authors': [],
            'journal': '',
            'year': '',
            'doi': '',
            'abstract': '',
            'mesh_terms': [],
            'keywords': [],
            'publication_type': [],
            'abstract_text': ''
        }

        # Title
        title_elem = article.find('.//ArticleTitle')
        if title_elem is not None and title_elem.text:
            metadata['title'] = title_elem.text.strip()

        # Authors
        authors_section = article.find('.//AuthorList')
        if authors_section is not None:
            for author_elem in authors_section.findall('Author'):
                last_name = author_elem.findtext('LastName', '')
                initials = author_elem.findtext('Initials', '')
                if last_name:
                    metadata['authors'].append(f"{last_name} {initials}".strip())

        # Journal
        journal_elem = article.find('.//Journal/Title')
        if journal_elem is not None and journal_elem.text:
            metadata['journal'] = journal_elem.text.strip()

        # Publication date
        pub_date = article.find('.//PubDate')
        if pub_date is not None:
            year = pub_date.findtext('Year', '')
            if year:
                metadata['year'] = year

        # DOI
        doi_elem = article.find(".//ArticleId[@IdType='doi']")
        if doi_elem is not None and doi_elem.text:
            metadata['doi'] = doi_elem.text.strip()

        # Abstract
        abstract_elem = article.find('.//Abstract')
        if abstract_elem is not None:
            abstract_texts = []
            for section in abstract_elem.findall('.//AbstractText'):
                if section.text:
                    abstract_texts.append(section.text.strip())
            metadata['abstract'] = ' '.join(abstract_texts)
            metadata['abstract_text'] = metadata['abstract']

        # MeSH terms
        mesh_list = article.find('.//MeshHeadingList')
        if mesh_list is not None:
            for mesh in mesh_list.findall('MeshHeading'):
                descriptor = mesh.find('DescriptorName')
                if descriptor is not None and descriptor.text:
                    metadata['mesh_terms'].append(descriptor.text.strip())

        # Keywords
        keyword_list = article.find('.//KeywordList')
        if keyword_list is not None:
            for keyword in keyword_list.findall('Keyword'):
                if keyword.text:
                    metadata['keywords'].append(keyword.text.strip())

        # Publication types
        pub_type_list = article.find('.//PublicationTypeList')
        if pub_type_list is not None:
            for pub_type in pub_type_list.findall('PublicationType'):
                if pub_type.text:
                    metadata['publication_type'].append(pub_type.text.strip())

        return metadata

    def fetch_doi(self, doi: str) -> Dict:
        """Fetch metadata from CrossRef API or resolve to PubMed."""
        self.log(f"Fetching metadata for DOI: {doi}...")

        # First try to resolve DOI to PubMed ID
        pmid = self._doi_to_pmid(doi)
        if pmid:
            self.log(f"Resolved DOI to PMID {pmid}")
            return self.fetch_pubmed(pmid)

        # Fallback to CrossRef
        return self._fetch_crossref(doi)

    def _doi_to_pmid(self, doi: str) -> Optional[int]:
        """Try to find PMID for a given DOI."""
        try:
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            params = {
                'db': 'pubmed',
                'term': f'{doi}[DOI]',
                'retmode': 'json'
            }
            response = requests.get(url, params=params, timeout=10)
            response.raise_for_status()
            data = response.json()
            ids = data.get('esearchresult', {}).get('idlist', [])
            if ids:
                return int(ids[0])
        except Exception as e:
            self.log(f"Error resolving DOI to PMID: {e}", 'WARN')
        return None

    def _fetch_crossref(self, doi: str) -> Dict:
        """Fetch metadata from CrossRef API."""
        metadata = {
            'identifier': doi,
            'identifier_type': 'doi',
            'title': '',
            'authors': [],
            'journal': '',
            'year': '',
            'doi': doi,
            'abstract': '',
            'mesh_terms': [],
            'keywords': [],
            'publication_type': []
        }

        try:
            url = f"https://api.crossref.org/works/{doi}"
            response = requests.get(url, timeout=10)
            response.raise_for_status()
            data = response.json().get('message', {})

            metadata['title'] = data.get('title', [''])[0] if data.get('title') else ''

            if 'author' in data:
                for author in data['author']:
                    name = f"{author.get('family', '')} {author.get('given', '')}".strip()
                    if name:
                        metadata['authors'].append(name)

            metadata['journal'] = data.get('container-title', [''])[0] if data.get('container-title') else ''
            metadata['year'] = str(data.get('issued', {}).get('date-parts', [['']])[0][0])
            metadata['abstract'] = data.get('abstract', '')

            self.log(f"Successfully extracted metadata from CrossRef: {metadata['title']}")
            return metadata

        except Exception as e:
            self.log(f"Failed to fetch from CrossRef: {e}", 'ERROR')
            return metadata

    def fetch_arxiv(self, arxiv_id: str) -> Dict:
        """Fetch metadata from arXiv."""
        self.log(f"Fetching arXiv metadata for ID: {arxiv_id}...")

        metadata = {
            'identifier': arxiv_id,
            'identifier_type': 'arxiv',
            'title': '',
            'authors': [],
            'journal': 'arXiv',
            'year': '',
            'doi': '',
            'abstract': '',
            'mesh_terms': [],
            'keywords': [],
            'publication_type': ['Preprint']
        }

        try:
            # Normalize arxiv ID
            clean_id = arxiv_id.replace('/', '').replace('arxiv:', '')
            url = f"https://export.arxiv.org/api/query"
            params = {'id_list': clean_id, 'start': 0, 'max_results': 1}

            response = requests.get(url, params=params, timeout=10)
            response.raise_for_status()

            # Parse Atom XML
            root = ET.fromstring(response.content)
            namespaces = {
                'atom': 'http://www.w3.org/2005/Atom',
                'arxiv': 'http://arxiv.org/schemas/atom'
            }

            entry = root.find('atom:entry', namespaces)
            if entry is None:
                self.log(f"No arXiv entry found for {arxiv_id}", 'WARN')
                return metadata

            # Title
            title_elem = entry.find('atom:title', namespaces)
            if title_elem is not None:
                metadata['title'] = title_elem.text.strip() if title_elem.text else ''

            # Authors
            for author in entry.findall('atom:author', namespaces):
                name_elem = author.find('atom:name', namespaces)
                if name_elem is not None and name_elem.text:
                    metadata['authors'].append(name_elem.text.strip())

            # Published date
            pub_elem = entry.find('atom:published', namespaces)
            if pub_elem is not None and pub_elem.text:
                metadata['year'] = pub_elem.text[:4]

            # Abstract
            summary = entry.find('atom:summary', namespaces)
            if summary is not None and summary.text:
                metadata['abstract'] = summary.text.strip()

            # Keywords (from categories)
            for category in entry.findall('atom:category', namespaces):
                term = category.get('term', '')
                if term:
                    metadata['keywords'].append(term)

            self.log(f"Successfully extracted metadata from arXiv: {metadata['title']}")
            return metadata

        except Exception as e:
            self.log(f"Failed to fetch from arXiv: {e}", 'ERROR')
            return metadata

    def extract_sections(self, abstract: str) -> Dict[str, str]:
        """Extract structured sections from abstract text."""
        sections = {
            'background': '',
            'methods': '',
            'results': '',
            'conclusions': ''
        }

        if not abstract:
            return sections

        sentences = [s.strip() for s in re.split(r'(?<=[.!?])\s+', abstract) if s.strip()]

        if len(sentences) == 0:
            return sections

        # Background: first 1-2 sentences
        sections['background'] = sentences[0]
        if len(sentences) > 1 and len(sections['background']) < 100:
            sections['background'] += ' ' + sentences[1]

        # Methods keywords
        method_keywords = [
            'used', 'applied', 'conducted', 'performed', 'employed',
            'analyzed', 'studied', 'examined', 'investigated',
            'rna-seq', 'crispr', 'sequencing', 'pcr', 'blot',
            'statistical', 'regression', 'correlation', 'cohort'
        ]

        # Results keywords
        results_keywords = ['found', 'identified', 'showed', 'revealed', 'demonstrated',
                           'increased', 'decreased', 'associated', 'significant']

        method_sents = []
        results_sents = []

        for sent in sentences[1:-1]:
            sent_lower = sent.lower()
            if any(kw in sent_lower for kw in method_keywords):
                method_sents.append(sent)
            elif any(kw in sent_lower for kw in results_keywords) or re.search(r'\d+[%]|\d+\.\d+\s*fold', sent):
                results_sents.append(sent)

        if method_sents:
            sections['methods'] = ' '.join(method_sents[:2])

        if results_sents:
            sections['results'] = ' '.join(results_sents[:3])

        # Conclusions: last 1-2 sentences
        sections['conclusions'] = sentences[-1]
        if len(sentences) > 1:
            sections['conclusions'] = sentences[-2] + ' ' + sentences[-1]

        return sections

    def extract_entities(self, text: str, topic_keywords: Optional[str] = None) -> Dict[str, List[str]]:
        """Extract genes, drugs, diseases from text."""
        entities = {
            'genes': [],
            'drugs': [],
            'diseases': [],
            'topic_mentions': []
        }

        if not text:
            return entities

        # Gene extraction (ALL_CAPS patterns)
        gene_candidates = set(self.GENE_PATTERN.findall(text))
        # Filter: common acronyms and abbreviations to exclude
        exclude = {'AND', 'OR', 'NOT', 'THE', 'BY', 'FOR', 'WITH', 'FROM', 'WAS', 'WERE'}
        entities['genes'] = [g for g in gene_candidates if g not in exclude and len(g) >= 2][:30]

        # Topic keywords mentions
        if topic_keywords:
            topic_list = [k.strip() for k in topic_keywords.split(',')]
            for topic in topic_list:
                if topic.lower() in text.lower():
                    entities['topic_mentions'].append(topic)

        return entities

    def extract_statistics(self, text: str) -> Dict[str, List]:
        """Extract p-values, fold changes, sample sizes."""
        stats = {
            'pvalues': [],
            'fold_changes': [],
            'sample_sizes': []
        }

        if not text:
            return stats

        # P-values
        pval_matches = self.PVALUE_PATTERN.findall(text)
        stats['pvalues'] = sorted(list(set(pval_matches)))[:10]

        # Fold changes
        fc_matches = re.findall(self.FOLD_CHANGE_PATTERN, text)
        stats['fold_changes'] = sorted(list(set(fc_matches)), reverse=True)[:10]

        # Sample sizes
        n_matches = re.findall(r'[nN]\s*[=:]\s*(\d+)', text)
        stats['sample_sizes'] = [int(n) for n in n_matches]

        return stats

    def assess_impact(self, metadata: Dict) -> Dict:
        """Assess paper impact based on journal and publication type."""
        impact = {
            'journal_tier': 'Unknown',
            'is_high_impact': False,
            'publication_category': ''
        }

        journal = metadata.get('journal', '')
        pub_types = metadata.get('publication_type', [])

        # Check high-impact journals
        for hij in self.HIGH_IMPACT_JOURNALS:
            if hij.lower() in journal.lower():
                impact['journal_tier'] = 'High Impact'
                impact['is_high_impact'] = True
                break
        else:
            impact['journal_tier'] = 'Standard'

        # Categorize publication type
        for ptype in pub_types:
            if 'Review' in ptype or 'Meta' in ptype:
                impact['publication_category'] = 'Review/Meta-analysis'
                break
            elif 'Clinical Trial' in ptype:
                impact['publication_category'] = 'Clinical Trial'
                break
        else:
            impact['publication_category'] = 'Original Research'

        return impact

    def generate_markdown(self, metadata: Dict, sections: Dict, entities: Dict,
                         statistics: Dict, impact: Dict) -> str:
        """Generate structured markdown report."""
        title = metadata.get('title', 'Unknown Paper')
        authors = ', '.join(metadata.get('authors', [])[:5])
        if len(metadata.get('authors', [])) > 5:
            authors += f" + {len(metadata['authors']) - 5} more"

        md = f"""# Paper Digest: {title}

## Citation
- **Authors:** {authors}
- **Journal:** {metadata.get('journal', 'Unknown')}
- **Year:** {metadata.get('year', 'Unknown')}
- **DOI:** {metadata.get('doi', 'Not available')}
- **Identifier:** {metadata.get('identifier_type', '').upper()}: {metadata.get('identifier', '')}

## Abstract
{metadata.get('abstract', 'No abstract available')}

## Key Findings
"""

        if sections.get('results'):
            for sent in sections['results'].split('. '):
                if sent.strip():
                    md += f"- {sent.strip()}\n"
        else:
            md += "- No structured findings extracted.\n"

        md += f"""
## Methods Summary
{sections.get('methods', 'Methods information not separately extracted from abstract.')}

## Main Results
{sections.get('results', 'Results not separately extracted.')}

## Key Statistics
"""
        if statistics['pvalues']:
            md += f"- **P-values:** {', '.join(statistics['pvalues'][:5])}\n"
        if statistics['fold_changes']:
            md += f"- **Fold changes:** {', '.join(statistics['fold_changes'][:5])}\n"
        if statistics['sample_sizes']:
            md += f"- **Sample sizes (n):** {', '.join(str(n) for n in statistics['sample_sizes'][:5])}\n"

        if not any([statistics['pvalues'], statistics['fold_changes'], statistics['sample_sizes']]):
            md += "- No statistical measures extracted from abstract.\n"

        md += f"""
## Entities Mentioned
"""
        if entities['genes']:
            md += f"- **Genes/Proteins:** {', '.join(entities['genes'][:15])}\n"
        if entities['topic_mentions']:
            md += f"- **Topic Keywords:** {', '.join(entities['topic_mentions'])}\n"
        if metadata.get('mesh_terms'):
            md += f"- **MeSH Terms:** {', '.join(metadata['mesh_terms'][:10])}\n"

        md += f"""
## Conclusions
{sections.get('conclusions', 'Conclusion not separately extracted.')}

## Keywords & Classification
- **Publication Type:** {', '.join(metadata.get('publication_type', ['Unknown']))}
- **Keywords:** {', '.join(metadata.get('keywords', [])[:10])}

## Research Impact Notes
- **Journal Impact Tier:** {impact['journal_tier']}
- **Publication Category:** {impact['publication_category']}
- **High Impact Journal:** {'Yes' if impact['is_high_impact'] else 'No'}

---
*Report generated: {datetime.now().isoformat(timespec='seconds')}*
"""
        return md

    def digest(self, pmid: Optional[int] = None, doi: Optional[str] = None,
               arxiv: Optional[str] = None, pdf: Optional[str] = None,
               output_format: str = 'markdown',
               include_citations: int = 0,
               topic_keywords: Optional[str] = None) -> bool:
        """Main digest pipeline."""

        # Fetch metadata
        metadata = {}
        if pmid:
            metadata = self.fetch_pubmed(pmid)
        elif doi:
            metadata = self.fetch_doi(doi)
        elif arxiv:
            metadata = self.fetch_arxiv(arxiv)
        elif pdf:
            self.log(f"Note: PDF processing limited; metadata from filename/headers only", 'WARN')
            metadata = {'identifier': pdf, 'identifier_type': 'pdf', 'title': Path(pdf).stem}
        else:
            self.log("No input source provided", 'ERROR')
            return False

        if not metadata.get('title'):
            self.log("Failed to retrieve paper metadata", 'ERROR')
            return False

        # Extract sections from abstract
        abstract = metadata.get('abstract', '')
        sections = self.extract_sections(abstract)

        # Extract entities
        entities = self.extract_entities(abstract + ' ' + metadata.get('title', ''), topic_keywords)

        # Extract statistics
        statistics = self.extract_statistics(abstract)

        # Assess impact
        impact = self.assess_impact(metadata)

        # Save metadata as JSON
        metadata_file = self.outdir / f"paper_metadata_{metadata['identifier']}.json"
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2, default=str)
        self.log(f"Metadata saved to {metadata_file}")

        # Generate and save report
        if output_format == 'markdown' or output_format == 'md':
            md_report = self.generate_markdown(metadata, sections, entities, statistics, impact)
            report_file = self.outdir / f"paper_digest_{metadata['identifier']}.md"
            with open(report_file, 'w') as f:
                f.write(md_report)
            self.log(f"Markdown report saved to {report_file}")

        elif output_format == 'json':
            report_data = {
                'metadata': metadata,
                'sections': sections,
                'entities': entities,
                'statistics': statistics,
                'impact': impact
            }
            report_file = self.outdir / f"paper_digest_{metadata['identifier']}.json"
            with open(report_file, 'w') as f:
                json.dump(report_data, f, indent=2)
            self.log(f"JSON report saved to {report_file}")

        elif output_format == 'txt':
            txt_report = f"""PAPER DIGEST REPORT
Title: {metadata['title']}
Authors: {', '.join(metadata['authors'][:5])}
Journal: {metadata['journal']}
Year: {metadata['year']}

ABSTRACT:
{abstract}

KEY FINDINGS:
{sections['results']}

ENTITIES: {', '.join(entities['genes'][:10])}
"""
            report_file = self.outdir / f"paper_digest_{metadata['identifier']}.txt"
            with open(report_file, 'w') as f:
                f.write(txt_report)
            self.log(f"Text report saved to {report_file}")

        # Save processing log
        log_file = self.outdir / f"extraction_log_{metadata['identifier']}.txt"
        with open(log_file, 'w') as f:
            f.write('\n'.join(self.log_messages))

        return True


def main():
    parser = argparse.ArgumentParser(
        description='Digest a single scientific paper into a structured summary',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python paper_digest_single.py --pmid 35769483 --outdir ./results
  python paper_digest_single.py --doi 10.1038/s41586-022-04826-7 --output-format json
  python paper_digest_single.py --arxiv 2301.12345 --topic-keywords "CRISPR,immunotherapy"
  python paper_digest_single.py --pdf ./paper.pdf --outdir ./digest
        """
    )

    parser.add_argument('--pmid', type=int, help='PubMed ID')
    parser.add_argument('--doi', type=str, help='Digital Object Identifier')
    parser.add_argument('--arxiv', type=str, help='arXiv ID')
    parser.add_argument('--pdf', type=str, help='Local PDF file path')
    parser.add_argument('--outdir', default='./paper_digest', help='Output directory')
    parser.add_argument('--output-format', choices=['markdown', 'json', 'txt'],
                       default='markdown', help='Output format')
    parser.add_argument('--include-citations', type=int, default=0,
                       help='Include top N cited references')
    parser.add_argument('--topic-keywords', type=str,
                       help='Comma-separated topic keywords for extraction')

    args = parser.parse_args()

    if not any([args.pmid, args.doi, args.arxiv, args.pdf]):
        parser.error('At least one of --pmid, --doi, --arxiv, or --pdf must be provided')

    digester = PaperDigester(args.outdir)
    success = digester.digest(
        pmid=args.pmid,
        doi=args.doi,
        arxiv=args.arxiv,
        pdf=args.pdf,
        output_format=args.output_format,
        include_citations=args.include_citations,
        topic_keywords=args.topic_keywords
    )

    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
