#!/usr/bin/env python3
"""Normalize curated gene rows with UniProt accessions using the UniProt REST API."""

from __future__ import annotations

import argparse
import csv
import json
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, List


UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"


def uniprot_lookup(gene_name: str, organism_id: int) -> Dict[str, str]:
    query = f"gene:{gene_name} AND organism_id:{organism_id} AND reviewed:true"
    params = {
        "query": query,
        "format": "json",
        "fields": "accession,protein_name,gene_names,organism_name",
        "size": "1",
    }
    url = f"{UNIPROT_SEARCH}?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(url, timeout=30) as response:
        payload = json.loads(response.read().decode("utf-8"))
    results: List[Dict[str, object]] = payload.get("results", [])
    if not results:
        return {}
    result = results[0]
    protein = result.get("proteinDescription", {})
    recommended = protein.get("recommendedName", {}) if isinstance(protein, dict) else {}
    full_name = recommended.get("fullName", {}).get("value", "") if isinstance(recommended, dict) else ""
    return {
        "UniProt accession": result.get("primaryAccession", ""),
        "Protein ID": full_name,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Curated CSV with Gene name column")
    parser.add_argument("--organism-id", type=int, default=9606, help="NCBI organism id, default human")
    parser.add_argument("--output", required=True, help="Output CSV path")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    with input_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = reader.fieldnames or []

    for row in rows:
        gene_name = (row.get("Gene name") or "").strip()
        if not gene_name:
            continue
        lookup = uniprot_lookup(gene_name, args.organism_id)
        if lookup.get("UniProt accession") and not row.get("UniProt accession"):
            row["UniProt accession"] = lookup["UniProt accession"]
        if lookup.get("Protein ID") and not row.get("Protein ID"):
            row["Protein ID"] = lookup["Protein ID"]

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(output_path)


if __name__ == "__main__":
    main()
