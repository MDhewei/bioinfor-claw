#!/usr/bin/env python3
"""Run a reproducible PubMed-centered discovery pass for a gene-list curation task.

This script does not decide final biological membership. It standardizes the repeatable part:
- build a fixed query set from category and organism
- search PubMed via E-utilities
- fetch summaries for the top PMIDs
- write reproducible outputs for later human/agent review
"""

from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, Iterable, List
from xml.etree import ElementTree


PUBMED_SEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_SUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def ncbi_get(url: str, params: Dict[str, str]) -> bytes:
    request_url = f"{url}?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(request_url, timeout=30) as response:
        return response.read()


def build_queries(category: str, organism: str) -> List[str]:
    return [
        f"{category} review {organism} genes proteins",
        f"{category} comprehensive review {organism}",
        f"{category} key regulators {organism}",
        f"{category} domain family review {organism}",
    ]


def search_pubmed(term: str, email: str, retmax: int) -> List[str]:
    payload = {
        "db": "pubmed",
        "term": term,
        "retmode": "xml",
        "retmax": str(retmax),
        "sort": "relevance",
        "email": email,
    }
    root = ElementTree.fromstring(ncbi_get(PUBMED_SEARCH, payload))
    return [node.text for node in root.findall("./IdList/Id") if node.text]


def fetch_summaries(pmids: Iterable[str], email: str) -> Dict[str, Dict[str, str]]:
    pmid_list = list(dict.fromkeys(pmids))
    if not pmid_list:
        return {}
    payload = {
        "db": "pubmed",
        "id": ",".join(pmid_list),
        "retmode": "json",
        "email": email,
    }
    data = json.loads(ncbi_get(PUBMED_SUMMARY, payload).decode("utf-8"))
    result = data.get("result", {})
    return {pmid: result.get(pmid, {}) for pmid in pmid_list}


def write_jsonl(path: Path, rows: Iterable[Dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=True) + "\n")


def write_candidate_csv(path: Path) -> None:
    headers = [
        "Gene name",
        "Protein ID",
        "Organism",
        "Functional class",
        "Evidence or role",
        "UniProt accession",
        "PMID",
        "References",
        "Inclusion tier",
        "Notes",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--category", required=True, help="Requested biological category, e.g. transcription factors")
    parser.add_argument("--organism", default="Homo sapiens", help="Organism name")
    parser.add_argument("--email", required=True, help="Email for NCBI E-utilities")
    parser.add_argument("--retmax", type=int, default=20, help="PMIDs per query")
    parser.add_argument("--delay-seconds", type=float, default=0.4, help="Delay between PubMed requests")
    parser.add_argument("--output-dir", required=True, help="Directory for JSONL and CSV outputs")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    queries = build_queries(args.category, args.organism)
    search_rows = []
    all_pmids: List[str] = []
    for query in queries:
        pmids = search_pubmed(query, email=args.email, retmax=args.retmax)
        search_rows.append({"query": query, "pmids": pmids, "count": len(pmids)})
        all_pmids.extend(pmids)
        time.sleep(args.delay_seconds)

    summaries = fetch_summaries(all_pmids, email=args.email)

    write_jsonl(output_dir / "search_queries.jsonl", search_rows)
    write_jsonl(
        output_dir / "pubmed_summaries.jsonl",
        (
            {
                "pmid": pmid,
                "title": summary.get("title", ""),
                "pubdate": summary.get("pubdate", ""),
                "source": summary.get("source", ""),
                "authors": [author.get("name", "") for author in summary.get("authors", [])],
            }
            for pmid, summary in summaries.items()
        ),
    )
    write_candidate_csv(output_dir / "curated_candidates.csv")
    print(output_dir)


if __name__ == "__main__":
    main()
