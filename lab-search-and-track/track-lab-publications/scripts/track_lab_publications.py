#!/usr/bin/env python3
"""
Track Lab Publications Skill
Fetches and analyzes publications from PubMed/Europe PMC for a given PI and institution.
"""

import argparse
import json
import sys
import time
import urllib.parse
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests


class PubMedAPI:
    """Interface to PubMed eUtils API."""

    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    def __init__(self, email="bioinformatics@example.com"):
        self.email = email
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": "Mozilla/5.0"})

    def esearch(self, query, retmax=100, mindate=None, maxdate=None):
        """Search PubMed and return list of PMIDs."""
        params = {
            "db": "pubmed",
            "term": query,
            "retmax": retmax,
            "retmode": "json",
            "email": self.email,
        }
        if mindate:
            params["mindate"] = mindate
        if maxdate:
            params["maxdate"] = maxdate
            params["datetype"] = "pdat"

        url = f"{self.BASE_URL}/esearch.fcgi"
        response = self.session.get(url, params=params)
        response.raise_for_status()
        data = response.json()

        pmids = data.get("esearchresult", {}).get("idlist", [])
        total = int(data.get("esearchresult", {}).get("count", 0))
        return pmids, total

    def efetch(self, pmids):
        """Fetch full records for given PMIDs."""
        if not pmids:
            return []

        params = {
            "db": "pubmed",
            "id": ",".join(pmids),
            "retmode": "xml",
            "email": self.email,
        }

        url = f"{self.BASE_URL}/efetch.fcgi"
        response = self.session.get(url, params=params)
        response.raise_for_status()
        return response.text


class PubMedParser:
    """Parse PubMed XML responses."""

    @staticmethod
    def parse_efetch_xml(xml_text):
        """Parse eFetch XML and extract publication records."""
        try:
            root = ET.fromstring(xml_text)
        except ET.ParseError as e:
            print(f"Error parsing XML: {e}", file=sys.stderr)
            return []

        records = []
        for article in root.findall(".//PubmedArticle"):
            record = {}

            # PMID
            pmid_elem = article.find(".//PMID")
            record["pmid"] = pmid_elem.text if pmid_elem is not None else "N/A"

            # Title
            title_elem = article.find(".//ArticleTitle")
            record["title"] = title_elem.text if title_elem is not None else "N/A"

            # Year
            year_elem = article.find(".//PubDate/Year")
            if year_elem is None:
                year_elem = article.find(".//PubDate/MedlineDate")
            record["year"] = year_elem.text[:4] if year_elem is not None else "N/A"

            # Journal
            journal_elem = article.find(".//Journal/Title")
            record["journal"] = journal_elem.text if journal_elem is not None else "N/A"

            # Authors
            authors = []
            for author in article.findall(".//Author"):
                last_name = author.findtext("LastName", "")
                initials = author.findtext("Initials", "")
                if last_name:
                    authors.append(f"{last_name} {initials}".strip())
            record["authors"] = "; ".join(authors) if authors else "N/A"

            # Author affiliations (first author affiliation)
            affiliation = "N/A"
            first_author = article.find(".//Author[1]")
            if first_author is not None:
                aff_elem = first_author.find(".//Affiliation")
                if aff_elem is not None:
                    affiliation = aff_elem.text
            record["affiliation"] = affiliation

            # DOI
            doi = "N/A"
            for pid in article.findall(".//ArticleId"):
                if pid.get("IdType") == "doi":
                    doi = pid.text
                    break
            record["doi"] = doi

            # Abstract
            abstract_elem = article.find(".//Abstract/AbstractText")
            abstract = ""
            if abstract_elem is not None and abstract_elem.text:
                abstract = abstract_elem.text[:500]
            record["abstract"] = abstract

            # MeSH terms (keywords)
            mesh_terms = []
            for mesh in article.findall(".//MeshHeading/DescriptorName"):
                if mesh.text:
                    mesh_terms.append(mesh.text)
            record["mesh_terms"] = "; ".join(mesh_terms) if mesh_terms else "N/A"

            records.append(record)

        return records


class EuropePMCAPI:
    """Interface to Europe PMC API for citation counts."""

    BASE_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

    def __init__(self):
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": "Mozilla/5.0"})

    def get_citations(self, pmid):
        """Fetch citation count for a PMID from Europe PMC."""
        try:
            params = {
                "query": f"PMID:{pmid}",
                "format": "json",
                "pageSize": 1,
            }
            response = self.session.get(self.BASE_URL, params=params, timeout=5)
            response.raise_for_status()
            data = response.json()

            hits = data.get("hitCount", 0)
            if hits > 0:
                result = data.get("resultList", {}).get("result", [])
                if result:
                    return result[0].get("citedByCount", 0)
            return 0
        except Exception as e:
            print(f"  Warning: Could not fetch citation count for PMID {pmid}: {e}", file=sys.stderr)
            return None


def extract_keywords(records, top_n=15):
    """Extract top keywords from MeSH terms and abstracts."""
    all_terms = []

    # Collect MeSH terms
    for record in records:
        if record["mesh_terms"] != "N/A":
            terms = [t.strip() for t in record["mesh_terms"].split(";")]
            all_terms.extend(terms)

    # Add common words from abstracts (simple word frequency)
    stopwords = {
        "the", "a", "an", "and", "or", "in", "on", "at", "to", "for", "of",
        "is", "are", "was", "were", "be", "been", "have", "has", "do", "does",
        "did", "will", "would", "should", "could", "may", "might", "can",
        "this", "that", "these", "those", "i", "you", "he", "she", "it",
        "we", "they", "with", "by", "from", "as", "by", "using", "use",
    }

    for record in records:
        if record["abstract"] != "":
            words = record["abstract"].lower().split()
            for word in words:
                word = word.strip(".,!?;:\"'()[]{}").lower()
                if len(word) > 3 and word not in stopwords:
                    all_terms.append(word)

    # Count and rank
    term_counts = Counter(all_terms)
    return term_counts.most_common(top_n)


def build_coauthor_network(records):
    """Build co-author collaboration network."""
    coauthor_pairs = defaultdict(int)
    author_counts = Counter()

    for record in records:
        if record["authors"] == "N/A":
            continue
        authors = [a.strip() for a in record["authors"].split(";")]
        author_counts.update(authors)

        # Count pairs
        for i, a1 in enumerate(authors):
            for a2 in authors[i+1:]:
                pair = tuple(sorted([a1, a2]))
                coauthor_pairs[pair] += 1

    return author_counts, coauthor_pairs


def plot_timeline(records, output_path):
    """Create publication timeline bar chart."""
    year_counts = Counter()
    for record in records:
        year = record["year"]
        if year != "N/A":
            year_counts[int(year)] += 1

    if not year_counts:
        return

    years = sorted(year_counts.keys())
    counts = [year_counts[y] for y in years]

    fig, ax = plt.subplots(figsize=(10, 5), dpi=100)
    ax.bar(years, counts, color="steelblue", edgecolor="black", alpha=0.7)
    ax.set_xlabel("Year", fontsize=12)
    ax.set_ylabel("Number of Publications", fontsize=12)
    ax.set_title("Publications Per Year", fontsize=14, fontweight="bold")
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved timeline plot to {output_path}")


def plot_keywords(keywords, output_path):
    """Create keyword frequency distribution chart."""
    if not keywords:
        return

    terms, counts = zip(*keywords)

    fig, ax = plt.subplots(figsize=(10, 8), dpi=100)
    ax.barh(range(len(terms)), counts, color="coral", edgecolor="black", alpha=0.7)
    ax.set_yticks(range(len(terms)))
    ax.set_yticklabels(terms, fontsize=10)
    ax.set_xlabel("Frequency", fontsize=12)
    ax.set_title("Top Research Keywords/MeSH Terms", fontsize=14, fontweight="bold")
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved keyword plot to {output_path}")


def plot_coauthor_network(author_counts, coauthor_pairs, output_path):
    """Create co-author network visualization."""
    if not author_counts or len(author_counts) < 2:
        return

    # Keep only top collaborators
    top_authors = [a for a, _ in author_counts.most_common(20)]
    top_pairs = {p: c for p, c in coauthor_pairs.items() if p[0] in top_authors and p[1] in top_authors}

    if not top_pairs:
        return

    try:
        import networkx as nx

        G = nx.Graph()
        G.add_nodes_from(top_authors)
        for (a1, a2), weight in top_pairs.items():
            G.add_edge(a1, a2, weight=weight)

        fig, ax = plt.subplots(figsize=(12, 10), dpi=100)
        pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

        # Node sizes by collaboration count
        node_sizes = [author_counts.get(n, 1) * 100 for n in G.nodes()]

        # Edge widths by weight
        edges = G.edges()
        weights = [top_pairs.get((min(e), max(e)), 1) for e in edges]

        nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color="lightblue",
                               edgecolors="black", linewidths=1, ax=ax)
        nx.draw_networkx_edges(G, pos, width=[w * 0.5 for w in weights], alpha=0.5, ax=ax)
        nx.draw_networkx_labels(G, pos, font_size=8, ax=ax)

        ax.set_title("Co-Author Collaboration Network", fontsize=14, fontweight="bold")
        ax.axis("off")
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"  Saved co-author network to {output_path}")
    except ImportError:
        print("  Skipping network plot (networkx not installed)", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description="Track lab publications from PubMed")
    parser.add_argument("--pi-name", required=True, help="PI name (e.g., 'Jennifer Doudna')")
    parser.add_argument("--institution", default="", help="Institution name (optional)")
    parser.add_argument("--years-back", type=int, default=3, help="Years of publications to retrieve")
    parser.add_argument("--max-results", type=int, default=100, help="Max papers to fetch")
    parser.add_argument("--outdir", default=".", help="Output directory")
    parser.add_argument("--output-format", choices=["markdown", "tsv", "json"], default="markdown")
    parser.add_argument("--fetch-citations", action="store_true", help="Fetch citation counts")
    parser.add_argument("--co-author-network", action="store_true", help="Build co-author network")

    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Tracking publications for: {args.pi_name}")
    if args.institution:
        print(f"Institution: {args.institution}")
    print(f"Years back: {args.years_back}")
    print()

    # Build search query
    query = f"{args.pi_name}[Author]"
    if args.institution:
        query += f" AND {args.institution}[Affiliation]"

    # Calculate date range
    from datetime import datetime, timedelta
    max_date = datetime.now().strftime("%Y/%m/%d")
    min_date = (datetime.now() - timedelta(days=365*args.years_back)).strftime("%Y/%m/%d")

    # Search PubMed
    print(f"Searching PubMed with query: {query}")
    api = PubMedAPI()
    pmids, total = api.esearch(query, retmax=args.max_results, mindate=min_date.replace("/", ""))

    print(f"Found {total} papers matching query, fetching up to {len(pmids)} full records...")

    if not pmids:
        print("No papers found. Try a broader search (e.g., last name only).")
        return

    # Fetch full records
    records = []
    batch_size = 10
    for i in range(0, len(pmids), batch_size):
        batch = pmids[i:i+batch_size]
        xml_text = api.efetch(batch)
        batch_records = PubMedParser.parse_efetch_xml(xml_text)
        records.extend(batch_records)
        time.sleep(0.5)  # Be respectful to API

    print(f"Retrieved {len(records)} records")

    # Fetch citations if requested
    if args.fetch_citations:
        print("Fetching citation counts from Europe PMC...")
        epmcapi = EuropePMCAPI()
        for i, record in enumerate(records):
            if i % 10 == 0:
                print(f"  {i}/{len(records)}")
            citations = epmcapi.get_citations(record["pmid"])
            record["citations"] = citations if citations is not None else 0
            time.sleep(0.1)
    else:
        for record in records:
            record["citations"] = None

    # Analyze data
    print("Analyzing publication data...")
    keywords = extract_keywords(records, top_n=15)
    author_counts, coauthor_pairs = build_coauthor_network(records)

    # Create DataFrame
    df = pd.DataFrame(records)

    # Output formats
    print(f"Writing outputs to {outdir}...")

    if args.output_format in ["markdown", "tsv"]:
        # TSV output
        tsv_cols = ["pmid", "title", "authors", "journal", "year", "doi", "citations", "abstract"]
        df[tsv_cols].to_csv(outdir / "publications.tsv", sep="\t", index=False)
        print(f"  Saved TSV to publications.tsv")

    if args.output_format in ["markdown", "json"]:
        # JSON output
        df.to_json(outdir / "publications.json", orient="records", indent=2)
        print(f"  Saved JSON to publications.json")

    if args.output_format == "markdown":
        # Markdown report
        with open(outdir / "lab_report.md", "w") as f:
            f.write(f"# Publication Report: {args.pi_name}\n\n")
            if args.institution:
                f.write(f"**Institution:** {args.institution}\n\n")

            f.write(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

            f.write("## Summary\n\n")
            f.write(f"- **Total papers:** {len(records)}\n")
            f.write(f"- **Date range:** {records[0]['year'] if records else 'N/A'} - {records[-1]['year'] if records else 'N/A'}\n")
            f.write(f"- **Top journals:** {', '.join([j for j, _ in Counter([r['journal'] for r in records if r['journal'] != 'N/A']).most_common(5)])}\n")
            f.write(f"- **Top keywords:** {', '.join([k for k, _ in keywords[:5]])}\n\n")

            # Top collaborators
            if author_counts:
                f.write("## Top Collaborators\n\n")
                for author, count in author_counts.most_common(10):
                    f.write(f"- {author}: {count} papers\n")
                f.write("\n")

            # Publications by year
            f.write("## Publications by Year\n\n")
            year_counts = Counter([r["year"] for r in records if r["year"] != "N/A"])
            for year in sorted(year_counts.keys(), reverse=True):
                f.write(f"- {year}: {year_counts[year]} papers\n")
            f.write("\n")

            # All publications
            f.write("## Publications\n\n")
            for i, record in enumerate(records, 1):
                f.write(f"### {i}. {record['title']}\n\n")
                f.write(f"**Authors:** {record['authors']}\n\n")
                f.write(f"**Journal:** {record['journal']}\n\n")
                f.write(f"**Year:** {record['year']}\n\n")
                if record['doi'] != "N/A":
                    f.write(f"**DOI:** {record['doi']}\n\n")
                if record['citations'] is not None:
                    f.write(f"**Citations:** {record['citations']}\n\n")
                if record['abstract']:
                    f.write(f"**Abstract:** {record['abstract']}\n\n")
                f.write("---\n\n")

        print(f"  Saved markdown report to lab_report.md")

    # Generate plots
    plot_timeline(records, outdir / "publications_timeline.png")
    plot_keywords(keywords, outdir / "keyword_distribution.png")

    if args.co_author_network:
        plot_coauthor_network(author_counts, coauthor_pairs, outdir / "coauthor_network.png")

    # Print summary
    print("\n" + "="*50)
    print(f"PI: {args.pi_name}")
    print(f"Papers found: {len(records)}")
    print(f"Year range: {min([r['year'] for r in records if r['year'] != 'N/A']) if records else 'N/A'} - {max([r['year'] for r in records if r['year'] != 'N/A']) if records else 'N/A'}")
    top_journals = Counter([r['journal'] for r in records if r['journal'] != 'N/A']).most_common(5)
    print(f"Top journals: {', '.join([j for j, _ in top_journals])}")
    print("="*50)


if __name__ == "__main__":
    main()
