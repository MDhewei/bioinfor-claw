#!/usr/bin/env python3
"""Gene-list discovery pipeline.

Given a functional category (e.g. "DNA methylation readers", "transcription
factors", "RTK receptors"), discover candidate genes/proteins from authoritative
databases AND collect supporting PubMed evidence, then write a populated
candidate CSV ready for human/agent QC.

Data sources, in order of priority:
  1. UniProtKB (reviewed entries) — primary candidate source. Free-text
     query against protein name + comments + keywords for the requested
     organism.
  2. NCBI Gene — secondary candidate source by free-text term.
  3. PubMed (E-utilities) — supporting bibliographic evidence (review
     papers, methods papers) attached to the run for the agent / human
     to consult during curation.

Output files (written to --output-dir):
  - curated_candidates.csv   ← POPULATED with merged UniProt + Gene hits.
                               Each row has gene symbol, protein name,
                               UniProt accession, organism, evidence source,
                               and an `Inclusion tier` of `Candidate` (the
                               human/agent should re-classify to Core /
                               Extended / Disputed during review).
  - uniprot_hits.jsonl       ← raw UniProt response per hit
  - ncbi_gene_hits.jsonl     ← raw NCBI Gene response per hit
  - search_queries.jsonl     ← PubMed queries that were issued
  - pubmed_summaries.jsonl   ← PubMed paper metadata for evidence

This script does NOT make final biological membership decisions. It produces
the standardized, reproducible candidate set that downstream curation
(`normalize_with_uniprot.py`, manual review) refines.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
from xml.etree import ElementTree


# ---------------------------------------------------------------------------- #
# Endpoints
# ---------------------------------------------------------------------------- #
PUBMED_SEARCH    = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_SUMMARY   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
NCBI_GENE_SEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_GENE_SUMMARY= "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
UNIPROT_SEARCH   = "https://rest.uniprot.org/uniprotkb/search"

# Map a few common organism strings to NCBI taxonomy IDs.
ORGANISM_TAXIDS: Dict[str, int] = {
    "homo sapiens": 9606,
    "human":        9606,
    "mus musculus": 10090,
    "mouse":        10090,
    "rattus norvegicus": 10116,
    "rat":          10116,
    "danio rerio":  7955,
    "zebrafish":    7955,
    "drosophila melanogaster": 7227,
    "caenorhabditis elegans":  6239,
    "saccharomyces cerevisiae": 4932,
}


# ---------------------------------------------------------------------------- #
# Low-level HTTP
# ---------------------------------------------------------------------------- #
def _http_get(url: str, params: Dict[str, str], timeout: int = 30) -> bytes:
    request_url = f"{url}?{urllib.parse.urlencode(params)}"
    req = urllib.request.Request(
        request_url,
        headers={"User-Agent": "bioinfor-claw/curate-gene-list (python urllib)"},
    )
    with urllib.request.urlopen(req, timeout=timeout) as response:
        return response.read()


# ---------------------------------------------------------------------------- #
# UniProt — primary candidate source
# ---------------------------------------------------------------------------- #
def _uniprot_query_string(category: str, organism_id: int, reviewed_only: bool) -> str:
    """Build a UniProt advanced query covering protein name + family/domain
    keywords. Free text search is wrapped in parentheses so multi-word
    categories don't get tokenised away."""
    cat = category.strip().strip('"')
    base = f'("{cat}")'
    if organism_id > 0:
        base += f" AND organism_id:{organism_id}"
    if reviewed_only:
        base += " AND reviewed:true"
    return base


def query_uniprot(category: str, organism_id: int,
                  size: int = 200, reviewed_only: bool = True) -> List[Dict]:
    """Search UniProt for proteins matching `category` in the chosen organism.
    Returns a list of dicts with accession, gene_symbol, protein_name, etc."""
    query = _uniprot_query_string(category, organism_id, reviewed_only)
    params = {
        "query":  query,
        "format": "json",
        "fields": "accession,id,gene_names,protein_name,organism_name,"
                  "cc_function,keyword,xref_pubmed",
        "size":   str(min(size, 500)),
    }
    print(f"[UniProt] query: {query} (size={params['size']})", file=sys.stderr)
    try:
        raw = _http_get(UNIPROT_SEARCH, params, timeout=45)
    except Exception as e:
        print(f"[UniProt] request failed: {e}", file=sys.stderr)
        return []
    try:
        data = json.loads(raw.decode("utf-8"))
    except Exception as e:
        print(f"[UniProt] JSON parse failed: {e}", file=sys.stderr)
        return []
    out: List[Dict] = []
    for r in data.get("results", []):
        acc = r.get("primaryAccession", "")
        # Gene name (preferred symbol). UniProt sometimes lists multiple genes.
        gene_symbol = ""
        synonyms: List[str] = []
        for g in r.get("genes", []) or []:
            gn = g.get("geneName", {}).get("value", "")
            if gn and not gene_symbol:
                gene_symbol = gn
            for syn in g.get("synonyms", []) or []:
                v = syn.get("value", "")
                if v: synonyms.append(v)
        protein_name = ""
        pdesc = r.get("proteinDescription", {}) or {}
        rec = pdesc.get("recommendedName", {}) or {}
        protein_name = (rec.get("fullName", {}) or {}).get("value", "")
        organism_name = (r.get("organism", {}) or {}).get("scientificName", "")
        # Function blurb (first ~250 chars)
        function_text = ""
        for c in r.get("comments", []) or []:
            if c.get("commentType") == "FUNCTION":
                texts = c.get("texts", []) or []
                if texts:
                    function_text = (texts[0].get("value", "") or "")[:250]
                    break
        keywords = [k.get("name", "") for k in (r.get("keywords") or []) if k.get("name")]
        # Cross-referenced PubMed IDs (selected primary refs)
        pubmed_xrefs: List[str] = []
        for x in r.get("uniProtKBCrossReferences", []) or []:
            if x.get("database") == "PubMed":
                pid = x.get("id", "")
                if pid: pubmed_xrefs.append(pid)
        out.append({
            "source":         "UniProt",
            "uniprot_accession": acc,
            "gene_symbol":    gene_symbol or (synonyms[0] if synonyms else ""),
            "synonyms":       synonyms,
            "protein_name":   protein_name,
            "organism_name":  organism_name,
            "function_text":  function_text,
            "keywords":       keywords,
            "pubmed_xrefs":   pubmed_xrefs[:5],
        })
    print(f"[UniProt] {len(out)} hits", file=sys.stderr)
    return out


# ---------------------------------------------------------------------------- #
# NCBI Gene — secondary candidate source
# ---------------------------------------------------------------------------- #
def query_ncbi_gene(category: str, organism: str, email: str,
                    retmax: int = 200, delay: float = 0.4) -> List[Dict]:
    """Search NCBI Gene database for genes matching `category` in `organism`."""
    term = f'"{category}"[All Fields] AND "{organism}"[Organism]'
    params = {
        "db":     "gene",
        "term":   term,
        "retmode":"xml",
        "retmax": str(retmax),
        "sort":   "relevance",
        "email":  email,
    }
    print(f"[NCBI Gene] search: {term}", file=sys.stderr)
    try:
        raw = _http_get(NCBI_GENE_SEARCH, params, timeout=30)
    except Exception as e:
        print(f"[NCBI Gene] esearch failed: {e}", file=sys.stderr)
        return []
    try:
        root = ElementTree.fromstring(raw)
    except Exception as e:
        print(f"[NCBI Gene] esearch XML parse failed: {e}", file=sys.stderr)
        return []
    gene_ids = [n.text for n in root.findall("./IdList/Id") if n.text]
    if not gene_ids:
        print("[NCBI Gene] no hits", file=sys.stderr)
        return []
    time.sleep(delay)
    # esummary in JSON for richer fields
    sum_params = {
        "db":     "gene",
        "id":     ",".join(gene_ids),
        "retmode":"json",
        "email":  email,
    }
    try:
        raw = _http_get(NCBI_GENE_SUMMARY, sum_params, timeout=45)
        data = json.loads(raw.decode("utf-8"))
    except Exception as e:
        print(f"[NCBI Gene] esummary failed: {e}", file=sys.stderr)
        return []
    result = data.get("result", {})
    out: List[Dict] = []
    for gid in gene_ids:
        s = result.get(gid)
        if not isinstance(s, dict):
            continue
        out.append({
            "source":      "NCBI Gene",
            "gene_id":     gid,
            "gene_symbol": s.get("name", ""),
            "description": s.get("description", "") or s.get("summary", "")[:250],
            "other_aliases": s.get("otheraliases", ""),
            "organism_name": (s.get("organism", {}) or {}).get("scientificname", ""),
            "chromosome":  s.get("chromosome", ""),
        })
    print(f"[NCBI Gene] {len(out)} hits", file=sys.stderr)
    return out


# ---------------------------------------------------------------------------- #
# PubMed — supporting bibliographic evidence
# ---------------------------------------------------------------------------- #
def build_pubmed_queries(category: str, organism: str) -> List[str]:
    return [
        f"{category} review {organism}",
        f"{category} comprehensive review {organism}",
        f"{category} family {organism}",
        f"{category} domain annotation {organism}",
    ]


def search_pubmed(term: str, email: str, retmax: int) -> List[str]:
    payload = {
        "db":      "pubmed",
        "term":    term,
        "retmode": "xml",
        "retmax":  str(retmax),
        "sort":    "relevance",
        "email":   email,
    }
    try:
        raw = _http_get(PUBMED_SEARCH, payload, timeout=30)
        root = ElementTree.fromstring(raw)
        return [n.text for n in root.findall("./IdList/Id") if n.text]
    except Exception as e:
        print(f"[PubMed] esearch failed for '{term}': {e}", file=sys.stderr)
        return []


def fetch_pubmed_summaries(pmids: Iterable[str], email: str) -> Dict[str, Dict]:
    pmid_list = list(dict.fromkeys(pmids))
    if not pmid_list:
        return {}
    # Chunk to stay within URL length limits
    out: Dict[str, Dict] = {}
    for i in range(0, len(pmid_list), 200):
        chunk = pmid_list[i:i+200]
        payload = {
            "db":      "pubmed",
            "id":      ",".join(chunk),
            "retmode": "json",
            "email":   email,
        }
        try:
            raw = _http_get(PUBMED_SUMMARY, payload, timeout=45)
            data = json.loads(raw.decode("utf-8"))
            result = data.get("result", {})
            for p in chunk:
                if isinstance(result.get(p), dict):
                    out[p] = result[p]
        except Exception as e:
            print(f"[PubMed] esummary chunk failed: {e}", file=sys.stderr)
    return out


# ---------------------------------------------------------------------------- #
# Merge UniProt + NCBI Gene candidates into populated rows
# ---------------------------------------------------------------------------- #
def merge_candidates(uniprot_hits: List[Dict],
                     gene_hits: List[Dict],
                     category: str,
                     organism: str) -> List[Dict[str, str]]:
    """Merge by gene symbol (case-insensitive), preferring UniProt fields,
    backfilling from NCBI Gene. Yields rows ready for the candidate CSV."""
    # Index by upper-cased symbol
    rows_by_sym: Dict[str, Dict[str, str]] = {}
    for h in uniprot_hits:
        sym = (h.get("gene_symbol") or "").strip()
        if not sym:
            continue
        key = sym.upper()
        rows_by_sym[key] = {
            "Gene name":        sym,
            "Protein ID":       h.get("protein_name", ""),
            "Organism":         h.get("organism_name") or organism,
            "Functional class": category,
            "Evidence or role": (h.get("function_text") or "")[:250],
            "UniProt accession":h.get("uniprot_accession", ""),
            "PMID":             ";".join(h.get("pubmed_xrefs", [])),
            "References":       "UniProt:" + (h.get("uniprot_accession") or ""),
            "Inclusion tier":   "Candidate",
            "Notes":            ("aliases=" + ",".join(h.get("synonyms") or [])
                                 if h.get("synonyms") else ""),
        }
    for g in gene_hits:
        sym = (g.get("gene_symbol") or "").strip()
        if not sym:
            continue
        key = sym.upper()
        if key in rows_by_sym:
            # Backfill description if UniProt didn't have one
            if not rows_by_sym[key]["Evidence or role"]:
                rows_by_sym[key]["Evidence or role"] = (g.get("description") or "")[:250]
            # Append aliases to Notes
            extras = g.get("other_aliases", "")
            if extras:
                cur = rows_by_sym[key]["Notes"]
                rows_by_sym[key]["Notes"] = (cur + "; " if cur else "") + "ncbi_aliases=" + extras
            # Tag refs
            rows_by_sym[key]["References"] += ";NCBI_Gene:" + (g.get("gene_id") or "")
        else:
            rows_by_sym[key] = {
                "Gene name":        sym,
                "Protein ID":       "",
                "Organism":         g.get("organism_name") or organism,
                "Functional class": category,
                "Evidence or role": (g.get("description") or "")[:250],
                "UniProt accession":"",
                "PMID":             "",
                "References":       "NCBI_Gene:" + (g.get("gene_id") or ""),
                "Inclusion tier":   "Candidate",
                "Notes":            ("ncbi_aliases=" + g.get("other_aliases", "")
                                     if g.get("other_aliases") else ""),
            }
    # Stable order: UniProt-backed (i.e. has accession) first, then alphabetical
    rows = list(rows_by_sym.values())
    rows.sort(key=lambda r: (r["UniProt accession"] == "", r["Gene name"].upper()))
    return rows


# ---------------------------------------------------------------------------- #
# Output writers
# ---------------------------------------------------------------------------- #
CSV_HEADERS = [
    "Gene name", "Protein ID", "Organism", "Functional class",
    "Evidence or role", "UniProt accession", "PMID", "References",
    "Inclusion tier", "Notes",
]


def write_jsonl(path: Path, rows: Iterable) -> None:
    with path.open("w", encoding="utf-8") as fh:
        for row in rows:
            fh.write(json.dumps(row, ensure_ascii=True) + "\n")


def write_candidate_csv(path: Path, rows: List[Dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_HEADERS)
        writer.writeheader()
        for r in rows:
            writer.writerow({h: r.get(h, "") for h in CSV_HEADERS})


# ---------------------------------------------------------------------------- #
# Main
# ---------------------------------------------------------------------------- #
def resolve_organism_id(organism: str) -> int:
    return ORGANISM_TAXIDS.get(organism.strip().lower(), 0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--category", required=True,
                        help='Functional category (e.g. "DNA methylation readers", '
                             '"transcription factors", "RTK receptors")')
    parser.add_argument("--organism", default="Homo sapiens",
                        help='Organism scientific name (default: "Homo sapiens")')
    parser.add_argument("--email", required=True,
                        help="Email for NCBI E-utilities (required by NCBI policy)")
    parser.add_argument("--retmax", type=int, default=50,
                        help="Per-query result cap for PubMed (default: 50)")
    parser.add_argument("--uniprot-size", type=int, default=200,
                        help="Max UniProt hits to retrieve (default: 200, max: 500)")
    parser.add_argument("--ncbi-gene-retmax", type=int, default=200,
                        help="Max NCBI Gene hits to retrieve (default: 200)")
    parser.add_argument("--delay-seconds", type=float, default=0.4,
                        help="Delay between NCBI requests (default: 0.4s)")
    parser.add_argument("--include-unreviewed", action="store_true",
                        help="Include UniProt TrEMBL (unreviewed) entries — noisier")
    parser.add_argument("--skip-pubmed", action="store_true",
                        help="Skip PubMed evidence collection (faster)")
    parser.add_argument("--output-dir", required=True,
                        help="Directory for JSONL and CSV outputs")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    organism_id = resolve_organism_id(args.organism)
    if organism_id == 0:
        print(f"[warn] Unknown organism '{args.organism}'; UniProt query will not "
              f"filter by organism_id. Consider passing one of: "
              f"{', '.join(sorted(ORGANISM_TAXIDS.keys()))}",
              file=sys.stderr)

    # ── 1. UniProt: primary candidate source ──
    uniprot_hits = query_uniprot(
        category=args.category,
        organism_id=organism_id,
        size=args.uniprot_size,
        reviewed_only=not args.include_unreviewed,
    )
    write_jsonl(output_dir / "uniprot_hits.jsonl", uniprot_hits)

    # ── 2. NCBI Gene: secondary candidate source ──
    time.sleep(args.delay_seconds)
    gene_hits = query_ncbi_gene(
        category=args.category,
        organism=args.organism,
        email=args.email,
        retmax=args.ncbi_gene_retmax,
        delay=args.delay_seconds,
    )
    write_jsonl(output_dir / "ncbi_gene_hits.jsonl", gene_hits)

    # ── 3. Merge + write populated candidate CSV ──
    rows = merge_candidates(uniprot_hits, gene_hits, args.category, args.organism)
    write_candidate_csv(output_dir / "curated_candidates.csv", rows)

    # ── 4. PubMed evidence (optional) ──
    pubmed_pmids: List[str] = []
    if not args.skip_pubmed:
        time.sleep(args.delay_seconds)
        queries = build_pubmed_queries(args.category, args.organism)
        search_rows = []
        for q in queries:
            pmids = search_pubmed(q, email=args.email, retmax=args.retmax)
            search_rows.append({"query": q, "pmids": pmids, "count": len(pmids)})
            pubmed_pmids.extend(pmids)
            time.sleep(args.delay_seconds)
        write_jsonl(output_dir / "search_queries.jsonl", search_rows)
        summaries = fetch_pubmed_summaries(pubmed_pmids, email=args.email)
        write_jsonl(
            output_dir / "pubmed_summaries.jsonl",
            ({
                "pmid":     pmid,
                "title":    s.get("title", ""),
                "pubdate":  s.get("pubdate", ""),
                "source":   s.get("source", ""),
                "authors":  [a.get("name", "") for a in s.get("authors", [])],
            } for pmid, s in summaries.items()),
        )

    # ── 5. Stdout summary so the agent / caller sees results immediately ──
    n_uniprot = len(uniprot_hits)
    n_gene    = len(gene_hits)
    n_rows    = len(rows)
    n_pmids   = len(set(pubmed_pmids))
    print("\n=== Discovery summary ===")
    print(f"Category         : {args.category}")
    print(f"Organism         : {args.organism} (taxid={organism_id or 'unknown'})")
    print(f"UniProt hits     : {n_uniprot}")
    print(f"NCBI Gene hits   : {n_gene}")
    print(f"Candidate rows   : {n_rows}  ← curated_candidates.csv")
    print(f"PubMed evidence  : {n_pmids} unique PMIDs")
    print(f"Output dir       : {output_dir}")
    if n_rows == 0:
        print("\n[!] Zero candidates found. Likely causes:")
        print("    • category phrasing is too specific — try a broader term "
              "(e.g. 'methyl-CpG binding' instead of 'DNA methylation readers')")
        print("    • organism not recognised — pass canonical scientific name")
        print("    • UniProt/NCBI rate-limited or transiently down — retry in a few minutes")
        print("    • try --include-unreviewed to widen the UniProt search to TrEMBL")
        sys.exit(2)


if __name__ == "__main__":
    main()
