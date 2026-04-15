#!/usr/bin/env python3
"""
Find Collaborators Skill
Identifies potential research collaborators based on topic overlap from PubMed.
"""

import argparse
import sys
import time
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from datetime import datetime, timedelta
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

    def esearch(self, query, retmax=100, mindate=None):
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

            # Journal
            journal_elem = article.find(".//Journal/Title")
            record["journal"] = journal_elem.text if journal_elem is not None else "N/A"

            # Year
            year_elem = article.find(".//PubDate/Year")
            if year_elem is None:
                year_elem = article.find(".//PubDate/MedlineDate")
            record["year"] = year_elem.text[:4] if year_elem is not None else "N/A"

            # Authors with affiliations
            authors = []
            for author in article.findall(".//Author"):
                last_name = author.findtext("LastName", "")
                initials = author.findtext("Initials", "")
                aff = author.findtext("Affiliation", "N/A")
                if last_name:
                    authors.append({
                        "name": f"{last_name} {initials}".strip(),
                        "affiliation": aff,
                    })
            record["authors"] = authors

            # Abstract
            abstract_elem = article.find(".//Abstract/AbstractText")
            abstract = ""
            if abstract_elem is not None and abstract_elem.text:
                abstract = abstract_elem.text[:500]
            record["abstract"] = abstract

            # MeSH terms
            mesh_terms = []
            for mesh in article.findall(".//MeshHeading/DescriptorName"):
                if mesh.text:
                    mesh_terms.append(mesh.text)
            record["mesh_terms"] = "; ".join(mesh_terms) if mesh_terms else ""

            records.append(record)

        return records


def extract_keywords(records, top_n=10):
    """Extract top keywords from abstracts and MeSH terms."""
    all_terms = []

    # Collect MeSH terms
    for record in records:
        if record["mesh_terms"]:
            terms = [t.strip() for t in record["mesh_terms"].split(";")]
            all_terms.extend(terms)

    # Add common words from abstracts
    stopwords = {
        "the", "a", "an", "and", "or", "in", "on", "at", "to", "for", "of",
        "is", "are", "was", "were", "be", "been", "have", "has", "with", "by"
    }

    for record in records:
        if record["abstract"]:
            words = record["abstract"].lower().split()
            for word in words:
                word = word.strip(".,!?;:\"'()[]{}").lower()
                if len(word) > 3 and word not in stopwords:
                    all_terms.append(word)

    term_counts = Counter(all_terms)
    return term_counts.most_common(top_n)


def main():
    parser = argparse.ArgumentParser(description="Find research collaborators by topic")
    parser.add_argument("--topics", required=True, help="Comma-separated research topics")
    parser.add_argument("--years-back", type=int, default=3, help="Years back to search")
    parser.add_argument("--max-per-topic", type=int, default=100, help="Max papers per topic")
    parser.add_argument("--exclude-institution", default="", help="Exclude institution string")
    parser.add_argument("--min-papers", type=int, default=3, help="Min papers to qualify")
    parser.add_argument("--outdir", default=".", help="Output directory")
    parser.add_argument("--species-focus", choices=["human", "mouse", "both", "none"],
                       default="none", help="Filter by species")
    parser.add_argument("--country-filter", default="", help="Filter by country")

    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    topics = [t.strip() for t in args.topics.split(",")]
    print(f"Finding collaborators for topics: {', '.join(topics)}")
    print(f"Years back: {args.years_back}")
    print()

    # Calculate date range
    mindate = (datetime.now() - timedelta(days=365*args.years_back)).strftime("%Y%m%d")

    # Search for papers per topic
    api = PubMedAPI()
    all_pmids = set()
    pmid_to_topics = defaultdict(set)
    topic_pmids = {}

    for topic in topics:
        print(f"Searching for '{topic}'...")
        query = topic
        if args.species_focus == "human":
            query += " AND human"
        elif args.species_focus == "mouse":
            query += " AND mouse"

        pmids, total = api.esearch(query, retmax=args.max_per_topic, mindate=mindate)
        print(f"  Found {total} papers, fetching up to {len(pmids)}")

        topic_pmids[topic] = pmids[:args.max_per_topic]
        for pmid in pmids:
            all_pmids.add(pmid)
            pmid_to_topics[pmid].add(topic)

        time.sleep(0.5)

    all_pmids = list(all_pmids)
    print(f"\nTotal unique PMIDs: {len(all_pmids)}")

    # Fetch full records
    print("Fetching publication details...")
    all_records = []
    batch_size = 10
    for i in range(0, len(all_pmids), batch_size):
        batch = all_pmids[i:i+batch_size]
        xml_text = api.efetch(batch)
        batch_records = PubMedParser.parse_efetch_xml(xml_text)
        all_records.extend(batch_records)
        if (i // batch_size) % 10 == 0:
            print(f"  Processed {i}/{len(all_pmids)} records")
        time.sleep(0.5)

    print(f"Retrieved {len(all_records)} records")

    # Build author-topic matrix
    print("Building author profiles...")
    author_papers = defaultdict(lambda: defaultdict(list))
    author_affiliations = {}
    author_journals = defaultdict(Counter)

    for record in all_records:
        pmid = record["pmid"]
        topics_for_paper = pmid_to_topics.get(pmid, set())

        for author_info in record["authors"]:
            author_name = author_info["name"]
            affiliation = author_info["affiliation"]

            # Skip if excluded
            if args.exclude_institution and args.exclude_institution.lower() in affiliation.lower():
                continue

            # Skip if country filter doesn't match
            if args.country_filter and args.country_filter.lower() not in affiliation.lower():
                continue

            author_affiliations[author_name] = affiliation
            for topic in topics_for_paper:
                author_papers[author_name][topic].append({
                    "pmid": pmid,
                    "title": record["title"],
                    "journal": record["journal"],
                    "year": record["year"],
                    "abstract": record["abstract"],
                })

            # Track journals
            if record["journal"] != "N/A":
                author_journals[author_name][record["journal"]] += 1

    # Filter by min_papers and compute ranking
    print("Computing rankings...")
    author_scores = []

    for author_name, topics_dict in author_papers.items():
        total_papers = sum(len(pmids) for pmids in topics_dict.values())

        if total_papers < args.min_papers:
            continue

        num_topics = len(topics_dict)
        score = (num_topics * 10) + np.log1p(total_papers)

        top_journals = [j for j, _ in author_journals[author_name].most_common(3)]

        author_scores.append({
            "author_name": author_name,
            "affiliation": author_affiliations.get(author_name, "N/A"),
            "n_papers": total_papers,
            "topics_covered": num_topics,
            "topics": list(topics_dict.keys()),
            "rank_score": score,
            "top_journals": "; ".join(top_journals),
        })

    author_scores.sort(key=lambda x: x["rank_score"], reverse=True)

    print(f"Found {len(author_scores)} qualifying authors")
    print("\nTop 10 collaborators:")
    for i, author in enumerate(author_scores[:10], 1):
        print(f"  {i}. {author['author_name']} ({author['n_papers']} papers, {author['topics_covered']} topics)")

    # Create DataFrame
    df = pd.DataFrame(author_scores)

    # Output TSV
    output_cols = ["author_name", "affiliation", "n_papers", "topics_covered", "top_journals", "rank_score"]
    df[output_cols].to_csv(outdir / "collaborators.tsv", sep="\t", index=False)
    print(f"\nSaved TSV to {outdir}/collaborators.tsv")

    # Generate markdown report
    print("Generating markdown report...")
    with open(outdir / "collaborator_profiles.md", "w") as f:
        f.write(f"# Collaborator Profiles\n\n")
        f.write(f"**Topics:** {', '.join(topics)}\n\n")
        f.write(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"**Total authors:** {len(author_scores)}\n\n")

        for i, author_data in enumerate(author_scores[:20], 1):
            author_name = author_data["author_name"]
            f.write(f"## {i}. {author_name}\n\n")
            f.write(f"**Affiliation:** {author_data['affiliation']}\n\n")
            f.write(f"**Papers:** {author_data['n_papers']} | ")
            f.write(f"**Topics covered:** {author_data['topics_covered']} | ")
            f.write(f"**Score:** {author_data['rank_score']:.2f}\n\n")
            f.write(f"**Top journals:** {author_data['top_journals']}\n\n")
            f.write(f"**Topics:** {', '.join(author_data['topics'])}\n\n")

            # List papers for this author
            f.write("**Recent publications:**\n\n")
            for topic in author_data['topics']:
                papers = author_papers[author_name][topic]
                for j, paper in enumerate(papers[:5], 1):  # Show top 5 per topic
                    f.write(f"- {paper['title'][:100]}... ({paper['year']}, {paper['journal']})\n")
            f.write("\n---\n\n")

    print(f"Saved markdown report to {outdir}/collaborator_profiles.md")

    # Generate heatmap
    print("Generating heatmap...")
    top_authors = author_scores[:30]
    author_names = [a["author_name"] for a in top_authors]

    heatmap_data = np.zeros((len(author_names), len(topics)))
    for i, author_data in enumerate(top_authors):
        author_name = author_data["author_name"]
        for j, topic in enumerate(topics):
            if topic in author_papers[author_name]:
                heatmap_data[i, j] = len(author_papers[author_name][topic])

    fig, ax = plt.subplots(figsize=(10, 12), dpi=100)
    im = ax.imshow(heatmap_data, cmap="YlOrRd", aspect="auto")

    ax.set_xticks(range(len(topics)))
    ax.set_xticklabels(topics, rotation=45, ha="right")
    ax.set_yticks(range(len(author_names)))
    ax.set_yticklabels(author_names, fontsize=8)

    plt.colorbar(im, ax=ax, label="Number of papers")
    ax.set_title("Author-Topic Heatmap", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(outdir / "topic_author_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved heatmap to {outdir}/topic_author_heatmap.png")

    # Generate ranking chart
    print("Generating ranking chart...")
    top_n = min(20, len(author_scores))
    top_authors_ranked = author_scores[:top_n]
    names = [a["author_name"][:30] for a in top_authors_ranked]
    scores = [a["rank_score"] for a in top_authors_ranked]

    fig, ax = plt.subplots(figsize=(10, 10), dpi=100)
    y_pos = np.arange(len(names))
    ax.barh(y_pos, scores, color="steelblue", edgecolor="black", alpha=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel("Rank Score", fontsize=12)
    ax.set_title(f"Top {top_n} Collaborators by Score", fontsize=14, fontweight="bold")
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.3)
    plt.tight_layout()
    plt.savefig(outdir / "ranking_chart.png", dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved ranking chart to {outdir}/ranking_chart.png")

    print("\n" + "="*50)
    print(f"Topics: {', '.join(topics)}")
    print(f"Authors found: {len(author_scores)}")
    print(f"Qualifying authors (≥{args.min_papers} papers): {len(author_scores)}")
    print("="*50)


if __name__ == "__main__":
    main()
