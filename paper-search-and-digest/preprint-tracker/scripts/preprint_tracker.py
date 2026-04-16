#!/usr/bin/env python3
"""
preprint_tracker.py: Search and track recent preprints from bioRxiv and medRxiv
Identifies trending topics, prolific authors, and generates digest reports
"""

import argparse
import json
import sys
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import Counter
import re

import requests
import pandas as pd
import matplotlib.pyplot as plt
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), '..', '..', '..', '_shared'))
from plot_style import init_style


class PreprintTracker:
    """Search and track preprints from bioRxiv and medRxiv."""

    # Stopwords for keyword extraction
    STOPWORDS = {
        'the', 'a', 'an', 'and', 'or', 'but', 'in', 'on', 'at', 'to', 'for',
        'is', 'was', 'are', 'were', 'be', 'been', 'being', 'have', 'has', 'had',
        'do', 'does', 'did', 'will', 'would', 'could', 'should', 'may', 'might',
        'can', 'this', 'that', 'these', 'those', 'i', 'you', 'he', 'she', 'it',
        'we', 'they', 'what', 'which', 'who', 'when', 'where', 'why', 'how',
        'all', 'each', 'every', 'both', 'few', 'more', 'most', 'other', 'some',
        'such', 'no', 'nor', 'not', 'only', 'own', 'same', 'so', 'than', 'too',
        'very', 'by', 'from', 'up', 'about', 'of', 'with', 'as', 'into', 'out',
        'during', 'through', 'before', 'after', 'above', 'below', 'between'
    }

    def __init__(self, outdir: str = './preprint_tracker'):
        self.outdir = Path(outdir)
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.log_messages = []
        self.preprints = []

    def log(self, msg: str, level: str = 'INFO'):
        """Add message to processing log."""
        timestamp = datetime.now().isoformat(timespec='seconds')
        log_line = f"[{timestamp}] {level}: {msg}"
        self.log_messages.append(log_line)
        print(log_line)

    def fetch_biorxiv_preprints(self, date_from: str, date_to: str,
                                max_results: int = 100) -> List[Dict]:
        """Fetch preprints from bioRxiv API."""
        self.log(f"Fetching bioRxiv preprints from {date_from} to {date_to}...")

        preprints = []
        cursor = 0
        batch_size = 100

        while len(preprints) < max_results:
            try:
                url = f"https://api.biorxiv.org/details/biorxiv/{date_from}/{date_to}/{cursor}/json"
                response = requests.get(url, timeout=15)
                response.raise_for_status()
                data = response.json()

                if 'collection' not in data or not data['collection']:
                    break

                for item in data['collection']:
                    preprints.append({
                        'doi': item.get('doi', ''),
                        'title': item.get('title', ''),
                        'authors': self._parse_authors(item.get('authors', '')),
                        'category': item.get('category', 'Unknown'),
                        'date': item.get('date', ''),
                        'abstract': item.get('abstract', ''),
                        'server': 'bioRxiv'
                    })

                if len(preprints) >= max_results:
                    break

                cursor += batch_size
                self.log(f"Fetched {len(preprints)} bioRxiv preprints so far...")

            except Exception as e:
                self.log(f"Error fetching bioRxiv batch at cursor {cursor}: {e}", 'WARN')
                break

        self.log(f"Total bioRxiv preprints fetched: {len(preprints)}")
        return preprints[:max_results]

    def fetch_medrxiv_preprints(self, date_from: str, date_to: str,
                                max_results: int = 100) -> List[Dict]:
        """Fetch preprints from medRxiv API."""
        self.log(f"Fetching medRxiv preprints from {date_from} to {date_to}...")

        preprints = []
        cursor = 0
        batch_size = 100

        while len(preprints) < max_results:
            try:
                url = f"https://api.biorxiv.org/details/medrxiv/{date_from}/{date_to}/{cursor}/json"
                response = requests.get(url, timeout=15)
                response.raise_for_status()
                data = response.json()

                if 'collection' not in data or not data['collection']:
                    break

                for item in data['collection']:
                    preprints.append({
                        'doi': item.get('doi', ''),
                        'title': item.get('title', ''),
                        'authors': self._parse_authors(item.get('authors', '')),
                        'category': item.get('category', 'Unknown'),
                        'date': item.get('date', ''),
                        'abstract': item.get('abstract', ''),
                        'server': 'medRxiv'
                    })

                if len(preprints) >= max_results:
                    break

                cursor += batch_size
                self.log(f"Fetched {len(preprints)} medRxiv preprints so far...")

            except Exception as e:
                self.log(f"Error fetching medRxiv batch at cursor {cursor}: {e}", 'WARN')
                break

        self.log(f"Total medRxiv preprints fetched: {len(preprints)}")
        return preprints[:max_results]

    def _parse_authors(self, authors_str: str) -> List[str]:
        """Parse author string into list of names."""
        if not authors_str:
            return []

        # Split by semicolon or comma
        authors = re.split(r'[;,]', authors_str)
        # Clean and return
        return [a.strip() for a in authors if a.strip()][:10]  # Limit to first 10

    def filter_by_keywords(self, preprints: List[Dict], keywords: str) -> List[Dict]:
        """Filter preprints by keyword search in title and abstract."""
        keyword_list = [kw.strip().lower() for kw in keywords.split()]

        filtered = []
        for preprint in preprints:
            full_text = (preprint['title'] + ' ' + preprint['abstract']).lower()
            if any(kw in full_text for kw in keyword_list):
                filtered.append(preprint)

        self.log(f"Filtered to {len(filtered)} preprints matching keywords: {keywords}")
        return filtered

    def filter_by_category(self, preprints: List[Dict], category: str) -> List[Dict]:
        """Filter preprints by category."""
        filtered = [p for p in preprints if category.lower() in p['category'].lower()]
        self.log(f"Filtered to {len(filtered)} preprints in category: {category}")
        return filtered

    def extract_keywords(self, preprints: List[Dict], top_n: int = 20) -> List[Tuple[str, int]]:
        """Extract and rank keywords from all preprints."""
        all_text = ' '.join([p['title'] + ' ' + p['abstract'] for p in preprints])

        # Tokenize and clean
        words = re.findall(r'\b[a-z]+\b', all_text.lower())

        # Filter stopwords and short words
        filtered_words = [w for w in words if w not in self.STOPWORDS and len(w) > 3]

        # Count frequency
        word_freq = Counter(filtered_words)

        # Return top N
        return word_freq.most_common(top_n)

    def analyze_trends(self, preprints: List[Dict],
                      keywords: List[Tuple[str, int]]) -> Dict:
        """Analyze trending topics."""
        trends = {
            'top_keywords': keywords[:10],
            'top_keywords_with_papers': []
        }

        # Find top 3 papers for each top keyword
        for keyword, freq in keywords[:10]:
            papers_with_keyword = [
                p['title'] for p in preprints
                if keyword.lower() in (p['title'] + ' ' + p['abstract']).lower()
            ][:3]
            trends['top_keywords_with_papers'].append({
                'keyword': keyword,
                'frequency': freq,
                'top_papers': papers_with_keyword
            })

        return trends

    def identify_prolific_authors(self, preprints: List[Dict], top_n: int = 20) -> List[Tuple[str, int]]:
        """Identify authors with most preprints."""
        author_counts = Counter()

        for preprint in preprints:
            for author in preprint['authors']:
                author_counts[author] += 1

        return author_counts.most_common(top_n)

    def sort_preprints(self, preprints: List[Dict], sort_by: str = 'date') -> List[Dict]:
        """Sort preprints by date or relevance."""
        if sort_by == 'date':
            # Newest first
            return sorted(preprints, key=lambda x: x['date'], reverse=True)
        else:
            # Maintain order from API (relevance)
            return preprints

    def generate_markdown_report(self, preprints: List[Dict],
                                trends: Dict, prolific_authors: List[Tuple[str, int]]) -> str:
        """Generate formatted markdown report."""
        md = f"""# Preprint Tracker Report

**Generated:** {datetime.now().isoformat(timespec='seconds')}
**Total Preprints:** {len(preprints)}

## Trending Topics

"""
        for item in trends['top_keywords_with_papers'][:10]:
            md += f"- **{item['keyword']}** (frequency: {item['frequency']})\n"

        md += f"""

## Prolific Authors

"""
        for author, count in prolific_authors[:10]:
            md += f"- {author}: {count} preprints\n"

        md += f"""

## Preprints

"""
        for i, preprint in enumerate(preprints, 1):
            authors = ', '.join(preprint['authors'][:3])
            if len(preprint['authors']) > 3:
                authors += f" + {len(preprint['authors']) - 3} more"

            abstract_preview = preprint['abstract'][:200] + ('...' if len(preprint['abstract']) > 200 else '')

            md += f"""### {i}. {preprint['title']}

- **DOI:** {preprint['doi']}
- **Authors:** {authors}
- **Date:** {preprint['date']}
- **Category:** {preprint['category']}
- **Server:** {preprint['server']}

**Abstract:** {abstract_preview}

---

"""

        return md

    def generate_visualizations(self, preprints: List[Dict],
                               keywords: List[Tuple[str, int]]):
        """Generate visualization charts."""
        # Timeline: preprints per week
        self._plot_timeline(preprints)

        # Keyword frequency
        self._plot_keyword_frequency(keywords)

        # Category distribution
        self._plot_category_distribution(preprints)

    def _plot_timeline(self, preprints: List[Dict]):
        """Generate timeline chart."""
        if not preprints:
            return

        # Group by week
        week_counts = {}
        for preprint in preprints:
            date_str = preprint['date']
            if date_str:
                try:
                    date = datetime.strptime(date_str, '%Y-%m-%d')
                    week_key = date.strftime('%Y-W%U')
                    week_counts[week_key] = week_counts.get(week_key, 0) + 1
                except ValueError:
                    pass

        if week_counts:
            weeks = sorted(week_counts.keys())
            counts = [week_counts[w] for w in weeks]

            plt.figure(figsize=(14, 5))
            plt.bar(weeks, counts, color='steelblue', alpha=0.7)
            plt.xlabel('Week', fontsize=12)
            plt.ylabel('Number of Preprints', fontsize=12)
            plt.title('Publication Timeline (Weekly)', fontsize=14, fontweight='bold')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(self.outdir / 'timeline.png', dpi=100)
            plt.close()
            self.log("Saved timeline.png")

    def _plot_keyword_frequency(self, keywords: List[Tuple[str, int]]):
        """Generate keyword frequency chart."""
        if not keywords:
            return

        kwords = [k[0] for k in keywords[:20]]
        kfreqs = [k[1] for k in keywords[:20]]

        plt.figure(figsize=(12, 6))
        plt.barh(kwords, kfreqs, color='coral', alpha=0.7)
        plt.xlabel('Frequency', fontsize=12)
        plt.title('Top Keywords in Preprints', fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.savefig(self.outdir / 'topic_frequency.png', dpi=100)
        plt.close()
        self.log("Saved topic_frequency.png")

    def _plot_category_distribution(self, preprints: List[Dict]):
        """Generate category distribution chart."""
        category_counts = {}
        for preprint in preprints:
            cat = preprint['category']
            category_counts[cat] = category_counts.get(cat, 0) + 1

        if category_counts:
            # Top 10 categories
            top_cats = dict(sorted(category_counts.items(), key=lambda x: x[1], reverse=True)[:10])

            plt.figure(figsize=(12, 6))
            categories = list(top_cats.keys())
            counts = list(top_cats.values())
            plt.bar(categories, counts, color='lightgreen', alpha=0.7)
            plt.xlabel('Category', fontsize=12)
            plt.ylabel('Number of Preprints', fontsize=12)
            plt.title('Category Distribution', fontsize=14, fontweight='bold')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(self.outdir / 'category_distribution.png', dpi=100)
            plt.close()
            self.log("Saved category_distribution.png")

    def save_results(self, preprints: List[Dict], keywords: List[Tuple[str, int]],
                    prolific_authors: List[Tuple[str, int]], output_format: str = 'markdown'):
        """Save results to files."""
        # Prepare DataFrame for preprints
        preprint_data = []
        for preprint in preprints:
            preprint_data.append({
                'DOI': preprint['doi'],
                'Title': preprint['title'],
                'Authors': '; '.join(preprint['authors']),
                'Category': preprint['category'],
                'Date': preprint['date'],
                'Server': preprint['server'],
                'Word_Count': len(preprint['abstract'].split()),
                'Abstract_Preview': preprint['abstract'][:150] + ('...' if len(preprint['abstract']) > 150 else '')
            })

        df_preprints = pd.DataFrame(preprint_data)

        # Save preprints TSV
        tsv_file = self.outdir / 'preprints.tsv'
        df_preprints.to_csv(tsv_file, sep='\t', index=False)
        self.log(f"Saved preprints to {tsv_file}")

        # Save trending topics TSV
        keywords_data = []
        for keyword, freq in keywords:
            keywords_data.append({
                'Keyword': keyword,
                'Frequency': freq
            })

        df_keywords = pd.DataFrame(keywords_data)
        keywords_file = self.outdir / 'trending_topics.tsv'
        df_keywords.to_csv(keywords_file, sep='\t', index=False)
        self.log(f"Saved trending topics to {keywords_file}")

        # Save prolific authors TSV
        authors_data = []
        for author, count in prolific_authors:
            authors_data.append({
                'Author': author,
                'Preprint_Count': count
            })

        df_authors = pd.DataFrame(authors_data)
        authors_file = self.outdir / 'prolific_authors.tsv'
        df_authors.to_csv(authors_file, sep='\t', index=False)
        self.log(f"Saved prolific authors to {authors_file}")

        # Save markdown report
        trends = self.analyze_trends(preprints, keywords)
        md_report = self.generate_markdown_report(preprints, trends, prolific_authors)
        md_file = self.outdir / 'digest_report.md'
        with open(md_file, 'w') as f:
            f.write(md_report)
        self.log(f"Saved markdown report to {md_file}")

    def save_log(self):
        """Save processing log."""
        log_file = self.outdir / 'search_log.txt'
        with open(log_file, 'w') as f:
            f.write('\n'.join(self.log_messages))
        self.log(f"Saved log to {log_file}")

    def track(self, query: str, server: str = 'both', date_from: Optional[str] = None,
             date_to: Optional[str] = None, max_results: int = 100,
             category: Optional[str] = None, sort_by: str = 'date',
             output_format: str = 'markdown') -> bool:
        """Execute full preprint tracking pipeline."""

        # Set default dates (last 30 days)
        if date_to is None:
            date_to = datetime.now().strftime('%Y-%m-%d')
        if date_from is None:
            date_from = (datetime.now() - timedelta(days=30)).strftime('%Y-%m-%d')

        self.log(f"Tracking preprints for query: '{query}'")
        self.log(f"Date range: {date_from} to {date_to}")

        # Fetch preprints
        preprints = []

        if server in ['biorxiv', 'both']:
            biorxiv_preprints = self.fetch_biorxiv_preprints(date_from, date_to, max_results)
            preprints.extend(biorxiv_preprints)

        if server in ['medrxiv', 'both']:
            medrxiv_preprints = self.fetch_medrxiv_preprints(date_from, date_to, max_results)
            preprints.extend(medrxiv_preprints)

        if not preprints:
            self.log("No preprints found", 'ERROR')
            return False

        self.log(f"Total preprints fetched: {len(preprints)}")

        # Filter by keywords
        preprints = self.filter_by_keywords(preprints, query)

        if not preprints:
            self.log("No preprints matching query keywords", 'WARN')
            self.save_log()
            return False

        # Filter by category if provided
        if category:
            preprints = self.filter_by_category(preprints, category)

        if not preprints:
            self.log(f"No preprints found in category: {category}", 'WARN')
            self.save_log()
            return False

        # Extract keywords
        keywords = self.extract_keywords(preprints)

        # Identify prolific authors
        prolific_authors = self.identify_prolific_authors(preprints)

        # Sort preprints
        preprints = self.sort_preprints(preprints, sort_by)

        # Save results
        self.save_results(preprints, keywords, prolific_authors, output_format)

        # Generate visualizations
        self.generate_visualizations(preprints, keywords)

        # Save log
        self.save_log()

        # Print summary
        self.log(f"\n=== TRACKING SUMMARY ===")
        self.log(f"Query: {query}")
        self.log(f"Server(s): {server}")
        self.log(f"Date range: {date_from} to {date_to}")
        self.log(f"Preprints retrieved: {len(preprints)}")
        self.log(f"Top 5 keywords: {', '.join([k[0] for k in keywords[:5]])}")
        self.log(f"Top author: {prolific_authors[0][0]} ({prolific_authors[0][1]} preprints)")

        return True


def main():
    parser = argparse.ArgumentParser(
        description='Track recent preprints from bioRxiv and medRxiv',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python preprint_tracker.py --query "CRISPR cancer" --max-results 100
  python preprint_tracker.py --query "deep learning genomics" --server biorxiv --category bioinformatics
  python preprint_tracker.py --query "COVID-19 treatment" --server medrxiv --date-from 2024-02-01
        """
    )

    parser.add_argument('--query', type=str, required=True, help='Search keywords')
    parser.add_argument('--server', choices=['biorxiv', 'medrxiv', 'both'],
                       default='both', help='Preprint server(s)')
    parser.add_argument('--date-from', help='Start date (YYYY-MM-DD; default: 30 days ago)')
    parser.add_argument('--date-to', help='End date (YYYY-MM-DD; default: today)')
    parser.add_argument('--max-results', type=int, default=100, help='Max preprints to retrieve')
    parser.add_argument('--category', help='bioRxiv category filter')
    parser.add_argument('--sort-by', choices=['date', 'relevance'],
                       default='date', help='Sort order')
    parser.add_argument('--outdir', default='./preprint_tracker', help='Output directory')
    parser.add_argument('--output-format', choices=['tsv', 'markdown'],
                       default='markdown', help='Output format')

    args = parser.parse_args()

    init_style(
        font_family=getattr(args, 'font_family', None),
        font_size=getattr(args, 'font_size', None),
    )

    tracker = PreprintTracker(args.outdir)
    success = tracker.track(
        query=args.query,
        server=args.server,
        date_from=args.date_from,
        date_to=args.date_to,
        max_results=args.max_results,
        category=args.category,
        sort_by=args.sort_by,
        output_format=args.output_format
    )

    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
