#!/usr/bin/env python3

import argparse
import json
import os
import re
import sys
import time
from collections import defaultdict, Counter
from typing import Dict, List, Optional, Tuple, Set

import pandas as pd
import requests


OPENALEX_BASE = "https://api.openalex.org"
WORKS_URL = f"{OPENALEX_BASE}/works"
AUTHORS_URL = f"{OPENALEX_BASE}/authors"
INSTITUTIONS_URL = f"{OPENALEX_BASE}/institutions"

SESSION = requests.Session()
SESSION.headers.update(
    {
        "User-Agent": "search-big-labs-by-field/8.0",
        "Accept-Language": "en-US,en;q=0.9",
    }
)


def normalize_space(text: str) -> str:
    return re.sub(r"\s+", " ", str(text or "").strip())


def ensure_parent_dir(path: str) -> None:
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)


def safe_get(
    url: str,
    params: Optional[dict] = None,
    timeout: int = 60,
    max_retries: int = 4,
    sleep_base: float = 1.2,
) -> requests.Response:
    last_err = None
    for attempt in range(max_retries):
        try:
            r = SESSION.get(url, params=params, timeout=timeout)
            if r.status_code in {403, 429}:
                wait = min(sleep_base * (2 ** attempt), 20)
                print(f"[WARN] {r.status_code} for {url}; retrying in {wait:.1f}s", file=sys.stderr)
                time.sleep(wait)
                continue
            r.raise_for_status()
            return r
        except Exception as e:
            last_err = e
            wait = min(sleep_base * (2 ** attempt), 15)
            if attempt < max_retries - 1:
                print(f"[WARN] request failed for {url}; retrying in {wait:.1f}s: {e}", file=sys.stderr)
                time.sleep(wait)
            else:
                raise
    raise last_err


def safe_get_json(url: str, params: Optional[dict] = None, timeout: int = 60) -> dict:
    return safe_get(url, params=params, timeout=timeout).json()


def parse_keywords(keyword: str) -> List[str]:
    parts = [normalize_space(x).lower() for x in re.split(r"[;,]", keyword) if normalize_space(x)]
    if not parts and normalize_space(keyword):
        parts = [normalize_space(keyword).lower()]
    return list(dict.fromkeys(parts))


def format_openalex_id(full_id: str) -> str:
    return normalize_space(full_id).rstrip("/").split("/")[-1]


def reconstruct_abstract(inv_idx: Optional[dict]) -> str:
    if not inv_idx:
        return ""
    pairs = []
    for word, positions in inv_idx.items():
        for pos in positions:
            pairs.append((pos, word))
    pairs.sort(key=lambda x: x[0])
    return " ".join(word for _, word in pairs)


def flatten_topics(topics: List[dict], n: int = 5) -> str:
    out = []
    for t in topics[:n]:
        name = normalize_space(t.get("display_name", ""))
        score = t.get("score")
        if name:
            if isinstance(score, (int, float)):
                out.append(f"{name} ({score:.2f})")
            else:
                out.append(name)
    return "; ".join(out)


def top_key_from_counter_dict(d: Dict[str, int]) -> str:
    if not d:
        return ""
    return sorted(d.items(), key=lambda x: x[1], reverse=True)[0][0]


def score_work_relevance(work: dict, keyword_terms: List[str]) -> float:
    title = normalize_space(work.get("display_name", "")).lower()
    abstract = reconstruct_abstract(work.get("abstract_inverted_index")).lower()
    topic_text = " ".join(
        [normalize_space(t.get("display_name", "")) for t in work.get("topics", [])]
    ).lower()

    text = f"{title} {abstract} {topic_text}"
    score = 0.0

    for kw in keyword_terms:
        if kw in title:
            score += 4.0
        if kw in abstract:
            score += 3.0
        if kw in topic_text:
            score += 2.0

    cited = float(work.get("cited_by_count", 0) or 0)
    score += min(cited, 1000.0) / 100.0

    pub_year = work.get("publication_year")
    if isinstance(pub_year, int):
        current_year = pd.Timestamp.today().year
        age = max(current_year - pub_year, 0)
        if age <= 1:
            score += 2.0
        elif age <= 3:
            score += 1.0

    return round(score, 3)


def infer_lab_type(text: str) -> str:
    txt = text.lower()
    computational_hits = [
        "computational", "bioinformatics", "algorithm", "machine learning",
        "deep learning", "modeling", "software", "statistics", "data science",
        "systems biology", "network biology", "ai", "foundation model",
    ]
    experimental_hits = [
        "screen", "crispr", "microscopy", "sequencing", "proteomics",
        "cell biology", "molecular biology", "genome editing", "assay",
        "mouse", "zebrafish", "organoid", "imaging",
    ]

    comp = sum(1 for k in computational_hits if k in txt)
    exp = sum(1 for k in experimental_hits if k in txt)

    if comp > 0 and exp > 0:
        return "hybrid"
    if comp > 0:
        return "computational"
    if exp > 0:
        return "experimental"
    return "unknown"


def is_mega_collaboration(work: dict, max_authors: int = 20) -> bool:
    authorships = work.get("authorships", []) or []
    return len(authorships) > max_authors


def fetch_works_for_keyword(
    keyword: str,
    years_back: int,
    per_page: int,
    max_works: int,
    max_authors_per_paper: int,
) -> List[dict]:
    current_year = pd.Timestamp.today().year
    from_year = current_year - years_back + 1

    works = []
    page = 1

    while True:
        remaining = max_works - len(works)
        if remaining <= 0:
            break

        params = {
            "search": keyword,
            "per-page": min(per_page, remaining),
            "page": page,
            "filter": f"from_publication_date:{from_year}-01-01",
            "sort": "cited_by_count:desc",
            "select": ",".join(
                [
                    "id",
                    "doi",
                    "display_name",
                    "publication_year",
                    "publication_date",
                    "cited_by_count",
                    "authorships",
                    "topics",
                    "abstract_inverted_index",
                    "type",
                ]
            ),
        }

        data = safe_get_json(WORKS_URL, params=params)
        results = data.get("results", [])
        if not results:
            break

        for w in results:
            wt = normalize_space(w.get("type", "")).lower()
            if wt not in {"article", "preprint", "review"}:
                continue
            if is_mega_collaboration(w, max_authors=max_authors_per_paper):
                continue
            works.append(w)
            if len(works) >= max_works:
                break

        page += 1
        meta = data.get("meta", {})
        count = int(meta.get("count", 0) or 0)
        if len(results) < params["per-page"]:
            break
        if page > (count // max(per_page, 1)) + 3:
            break
        time.sleep(0.4)

    return works[:max_works]


def collect_author_candidates(works: List[dict], keyword_terms: List[str]) -> Dict[str, dict]:
    authors = {}

    for work in works:
        work_score = score_work_relevance(work, keyword_terms)
        work_title = normalize_space(work.get("display_name", ""))
        work_year = work.get("publication_year")
        work_cites = int(work.get("cited_by_count", 0) or 0)
        work_doi = normalize_space(work.get("doi", ""))
        work_id = normalize_space(work.get("id", ""))

        authorships = work.get("authorships", []) or []
        for a in authorships:
            author = a.get("author") or {}
            author_id = normalize_space(author.get("id", ""))
            if not author_id:
                continue

            author_name = normalize_space(author.get("display_name", ""))
            position = normalize_space(a.get("author_position", ""))
            institutions = a.get("institutions", []) or []

            inst_names = []
            inst_ids = []
            countries = []
            for inst in institutions:
                inst_name = normalize_space(inst.get("display_name", ""))
                inst_id = normalize_space(inst.get("id", ""))
                cc = normalize_space(inst.get("country_code", ""))
                if inst_name:
                    inst_names.append(inst_name)
                if inst_id:
                    inst_ids.append(inst_id)
                if cc:
                    countries.append(cc)

            if author_id not in authors:
                authors[author_id] = {
                    "author_id": author_id,
                    "author_name": author_name,
                    "query_relevant_works_count": 0,
                    "query_total_citations_from_relevant_works": 0,
                    "first_last_author_count": 0,
                    "last_author_count": 0,
                    "first_author_count": 0,
                    "middle_author_count": 0,
                    "relevance_score": 0.0,
                    "matched_keyword_terms": set(),
                    "institutions_seen": defaultdict(int),
                    "institution_ids_seen": defaultdict(int),
                    "countries_seen": defaultdict(int),
                    "top_works": [],
                }

            rec = authors[author_id]
            rec["query_relevant_works_count"] += 1
            rec["query_total_citations_from_relevant_works"] += work_cites
            rec["relevance_score"] += work_score

            if position in {"first", "last"}:
                rec["first_last_author_count"] += 1
            if position == "first":
                rec["first_author_count"] += 1
            elif position == "last":
                rec["last_author_count"] += 1
            else:
                rec["middle_author_count"] += 1

            title_l = work_title.lower()
            abstract_l = reconstruct_abstract(work.get("abstract_inverted_index")).lower()
            topic_l = " ".join([normalize_space(t.get("display_name", "")) for t in work.get("topics", [])]).lower()
            all_text = f"{title_l} {abstract_l} {topic_l}"
            for kw in keyword_terms:
                if kw in all_text:
                    rec["matched_keyword_terms"].add(kw)

            for n in inst_names:
                rec["institutions_seen"][n] += 1
            for i in inst_ids:
                rec["institution_ids_seen"][i] += 1
            for c in countries:
                rec["countries_seen"][c] += 1

            rec["top_works"].append(
                {
                    "work_id": work_id,
                    "doi": work_doi,
                    "title": work_title,
                    "year": work_year,
                    "cited_by_count": work_cites,
                    "position": position,
                    "score": work_score,
                }
            )

    for _, rec in authors.items():
        rec["top_works"] = sorted(
            rec["top_works"],
            key=lambda x: (x["score"], x["cited_by_count"]),
            reverse=True,
        )[:10]
        rec["relevance_score"] = round(rec["relevance_score"], 3)

    return authors


def filter_pi_like_candidates(
    author_candidates: Dict[str, dict],
    min_last_author_count: int = 2,
    min_first_last_author_count: int = 3,
) -> Dict[str, dict]:
    filtered = {}
    for k, v in author_candidates.items():
        if (
            v.get("last_author_count", 0) >= min_last_author_count
            or v.get("first_last_author_count", 0) >= min_first_last_author_count
        ):
            filtered[k] = v
    return filtered


def fetch_author_details(author_id: str) -> dict:
    params = {
        "select": ",".join(
            [
                "id",
                "display_name",
                "orcid",
                "works_count",
                "cited_by_count",
                "summary_stats",
                "ids",
                "last_known_institutions",
                "affiliations",
                "topics",
                "topic_share",
            ]
        )
    }
    url = f"{AUTHORS_URL}/{format_openalex_id(author_id)}"
    return safe_get_json(url, params=params)


def fetch_institution_details(inst_id: str) -> dict:
    params = {
        "select": ",".join(
            [
                "id",
                "display_name",
                "country_code",
                "type",
                "homepage_url",
                "ror",
                "geo",
            ]
        )
    }
    url = f"{INSTITUTIONS_URL}/{format_openalex_id(inst_id)}"
    return safe_get_json(url, params=params)


def fetch_top_works_for_author(
    author_id: str,
    years_back: int,
    max_works: int,
    max_authors_per_paper: int,
) -> List[dict]:
    author_short = format_openalex_id(author_id)
    current_year = pd.Timestamp.today().year
    from_year = current_year - years_back + 1

    works = []
    page = 1
    per_page = min(100, max_works)

    while True:
        remaining = max_works - len(works)
        if remaining <= 0:
            break

        params = {
            "filter": f"authorships.author.id:{author_short},from_publication_date:{from_year}-01-01",
            "per-page": min(per_page, remaining),
            "page": page,
            "sort": "publication_date:desc",
            "select": ",".join(
                [
                    "id",
                    "doi",
                    "display_name",
                    "publication_year",
                    "publication_date",
                    "cited_by_count",
                    "topics",
                    "abstract_inverted_index",
                    "type",
                    "authorships",
                ]
            ),
        }

        data = safe_get_json(WORKS_URL, params=params)
        results = data.get("results", [])
        if not results:
            break

        for w in results:
            wt = normalize_space(w.get("type", "")).lower()
            if wt not in {"article", "preprint", "review"}:
                continue
            if is_mega_collaboration(w, max_authors=max_authors_per_paper):
                continue
            works.append(w)
            if len(works) >= max_works:
                break

        if len(results) < params["per-page"]:
            break

        meta = data.get("meta", {})
        count = int(meta.get("count", 0) or 0)
        page += 1
        if page > (count // max(per_page, 1)) + 3:
            break
        time.sleep(0.3)

    return works[:max_works]


def filter_works_relevant_to_keyword(works: List[dict], keyword_terms: List[str]) -> List[dict]:
    filtered = []
    for w in works:
        score = score_work_relevance(w, keyword_terms)
        if score > 0:
            filtered.append((score, w))
    filtered.sort(key=lambda x: (x[0], x[1].get("cited_by_count", 0)), reverse=True)
    return [w for _, w in filtered]


def build_year_counter(works: List[dict]) -> Counter:
    c = Counter()
    for w in works:
        y = w.get("publication_year")
        if isinstance(y, int):
            c[y] += 1
    return c


def stringify_year_distribution(year_counter: Counter, top_n: int = 10) -> str:
    items = sorted(year_counter.items(), key=lambda x: x[0], reverse=True)[:top_n]
    return "; ".join([f"{year}:{count}" for year, count in items])


def stringify_works(works: List[dict], top_n: int = 5) -> str:
    out = []
    for w in works[:top_n]:
        title = normalize_space(w.get("display_name", "") or w.get("title", ""))
        year = w.get("publication_year") or w.get("year", "")
        cites = int(w.get("cited_by_count", 0) or 0)
        doi = normalize_space(w.get("doi", ""))
        if title:
            if doi:
                out.append(f"{title} [{year}; cites={cites}; doi={doi}]")
            else:
                out.append(f"{title} [{year}; cites={cites}]")
    return " | ".join(out)


def extract_main_methods_and_topics(
    topic_text: str,
    relevant_works: List[dict],
) -> Tuple[str, str]:
    text = " ".join(
        [topic_text]
        + [normalize_space(w.get("display_name", "")) for w in relevant_works[:20]]
        + [reconstruct_abstract(w.get("abstract_inverted_index")) for w in relevant_works[:10]]
    ).lower()

    methods = []
    topics = []

    method_terms = [
        "crispr", "single-cell", "single cell", "spatial", "proteomics",
        "transcriptomics", "genomics", "machine learning", "deep learning",
        "foundation model", "language model", "alphafold", "protein structure",
        "screen", "imaging", "sequencing", "computational",
    ]
    topic_terms = [
        "cancer", "immunity", "epigenetics", "development", "neuroscience",
        "evolution", "microbiome", "precision medicine", "functional genomics",
        "genome editing", "rna biology", "systems biology",
    ]

    for term in method_terms:
        if term in text:
            methods.append(term)
    for term in topic_terms:
        if term in text:
            topics.append(term)

    return "; ".join(dict.fromkeys(methods)), "; ".join(dict.fromkeys(topics))


def author_relevance_final_score(
    author_rec: dict,
    author_detail: dict,
) -> float:
    score = 0.0
    score += float(author_rec.get("relevance_score", 0.0))
    score += min(float(author_detail.get("cited_by_count", 0) or 0), 300000.0) / 3000.0
    score += min(float(author_detail.get("works_count", 0) or 0), 2000.0) / 200.0
    score += 3.0 * float(author_rec.get("last_author_count", 0))
    score += 1.5 * float(author_rec.get("first_author_count", 0))
    score += 1.0 * float(author_rec.get("first_last_author_count", 0))
    score += 0.5 * float(len(author_rec.get("matched_keyword_terms", set())))
    return round(score, 3)


def overlap_doi_set(row: dict) -> Set[str]:
    dois = set()
    for field in [
        "representative_papers_overall_10y_dois",
        "representative_keyword_relevant_papers_dois",
        "recent_papers_last_5y_dois",
    ]:
        val = normalize_space(row.get(field, ""))
        if val:
            for x in val.split(";"):
                x = normalize_space(x)
                if x:
                    dois.add(x.lower())
    return dois


def cap_rows_per_institution(rows: List[dict], max_per_institution: int = 2) -> List[dict]:
    kept = []
    counts = defaultdict(int)

    for row in sorted(rows, key=lambda x: x["relevance_score"], reverse=True):
        inst = normalize_space(row.get("institution_name", "")) or "UNKNOWN"
        if counts[inst] >= max_per_institution:
            continue
        kept.append(row)
        counts[inst] += 1
    return kept


def diversify_by_paper_overlap(
    rows: List[dict],
    top_n: int,
    max_overlap: int = 2,
) -> List[dict]:
    selected = []
    selected_doi_sets = []

    for row in sorted(rows, key=lambda x: x["relevance_score"], reverse=True):
        row_dois = overlap_doi_set(row)
        too_similar = False
        for seen in selected_doi_sets:
            overlap = len(row_dois & seen)
            if overlap > max_overlap:
                too_similar = True
                break
        if too_similar:
            continue
        selected.append(row)
        selected_doi_sets.append(row_dois)
        if len(selected) >= top_n:
            break
    return selected


def build_lab_rows(
    author_candidates: Dict[str, dict],
    keyword: str,
    top_n: int,
    sleep_seconds: float,
    author_works_limit: int,
    max_authors_per_paper: int,
) -> List[dict]:
    ranked = sorted(
        author_candidates.values(),
        key=lambda x: (
            x["relevance_score"],
            x["query_total_citations_from_relevant_works"],
            x["last_author_count"],
            x["query_relevant_works_count"],
        ),
        reverse=True,
    )[: min(top_n * 3, 120)]

    rows = []
    keyword_terms = parse_keywords(keyword)

    print(f"[INFO] building lab profiles for {len(ranked)} candidate PIs...", file=sys.stderr)

    for idx, rec in enumerate(ranked, start=1):
        author_id = rec["author_id"]
        print(f"[INFO] profiling {idx}/{len(ranked)}: {rec['author_name']}", file=sys.stderr)

        try:
            detail = fetch_author_details(author_id)
        except Exception as e:
            print(f"[WARN] failed author detail for {author_id}: {e}", file=sys.stderr)
            continue

        current_inst_name = top_key_from_counter_dict(rec["institutions_seen"])
        current_inst_id = top_key_from_counter_dict(rec["institution_ids_seen"])
        current_country = top_key_from_counter_dict(rec["countries_seen"])

        author_topics = detail.get("topics", []) or []
        topic_share = detail.get("topic_share", []) or []
        topic_text = flatten_topics(author_topics, n=5)
        topic_share_text = flatten_topics(topic_share, n=5)

        inst_homepage = ""
        inst_type = ""
        inst_ror = ""
        city = ""
        if current_inst_id:
            try:
                inst = fetch_institution_details(current_inst_id)
                inst_homepage = normalize_space(inst.get("homepage_url", ""))
                inst_type = normalize_space(inst.get("type", ""))
                inst_ror = normalize_space(inst.get("ror", ""))
                geo = inst.get("geo") or {}
                city = normalize_space(geo.get("city", ""))
                if not current_country:
                    current_country = normalize_space(inst.get("country_code", ""))
            except Exception as e:
                print(f"[WARN] failed institution detail for {current_inst_id}: {e}", file=sys.stderr)

        try:
            author_works_10y = fetch_top_works_for_author(
                author_id=author_id,
                years_back=10,
                max_works=author_works_limit,
                max_authors_per_paper=max_authors_per_paper,
            )
        except Exception as e:
            print(f"[WARN] failed recent works for {author_id}: {e}", file=sys.stderr)
            author_works_10y = []

        relevant_works_10y = filter_works_relevant_to_keyword(author_works_10y, keyword_terms)

        current_year = pd.Timestamp.today().year
        works_last_10y = [w for w in author_works_10y if isinstance(w.get("publication_year"), int) and w["publication_year"] >= current_year - 9]
        works_last_5y = [w for w in author_works_10y if isinstance(w.get("publication_year"), int) and w["publication_year"] >= current_year - 4]
        works_last_3y = [w for w in author_works_10y if isinstance(w.get("publication_year"), int) and w["publication_year"] >= current_year - 2]

        relevant_last_10y = [w for w in relevant_works_10y if isinstance(w.get("publication_year"), int) and w["publication_year"] >= current_year - 9]
        relevant_last_5y = [w for w in relevant_works_10y if isinstance(w.get("publication_year"), int) and w["publication_year"] >= current_year - 4]
        relevant_last_3y = [w for w in relevant_works_10y if isinstance(w.get("publication_year"), int) and w["publication_year"] >= current_year - 2]

        year_dist = build_year_counter(works_last_10y)

        top_cited_overall_10y = sorted(
            works_last_10y,
            key=lambda x: int(x.get("cited_by_count", 0) or 0),
            reverse=True,
        )
        top_recent_overall = sorted(
            works_last_5y,
            key=lambda x: (
                int(x.get("publication_year", 0) or 0),
                int(x.get("cited_by_count", 0) or 0),
            ),
            reverse=True,
        )
        top_relevant = sorted(
            relevant_works_10y,
            key=lambda x: (
                score_work_relevance(x, keyword_terms),
                int(x.get("cited_by_count", 0) or 0),
            ),
            reverse=True,
        )

        main_methods, main_topics = extract_main_methods_and_topics(topic_text + " " + topic_share_text, relevant_works_10y)

        lab_text = " ".join(
            [
                normalize_space(detail.get("display_name", "")),
                current_inst_name,
                topic_text,
                topic_share_text,
                main_methods,
                main_topics,
            ]
        )
        lab_type = infer_lab_type(lab_text)

        final_score = author_relevance_final_score(
            author_rec=rec,
            author_detail=detail,
        )

        summary_parts = [
            f"Relevant to {keyword}.",
            f"Overall impact: {detail.get('cited_by_count', 0)} total citations and {detail.get('works_count', 0)} total works in OpenAlex.",
            f"Recent output: {len(works_last_10y)} papers in last 10 years, {len(works_last_5y)} in last 5 years, {len(works_last_3y)} in last 3 years.",
            f"Keyword-relevant output: {len(relevant_last_10y)} papers in last 10 years and {len(relevant_last_5y)} in last 5 years.",
            f"PI-like authorship evidence: last-author={rec.get('last_author_count', 0)}, first/last-author={rec.get('first_last_author_count', 0)}.",
        ]
        if main_methods:
            summary_parts.append(f"Likely strengths include {main_methods}.")
        if main_topics:
            summary_parts.append(f"Main topical areas include {main_topics}.")
        research_summary = normalize_space(" ".join(summary_parts))

        row = {
            "query_keyword": keyword,
            "lab_name": f"{normalize_space(detail.get('display_name', ''))} Lab" if normalize_space(detail.get("display_name", "")) else "",
            "pi_name": normalize_space(detail.get("display_name", "")),
            "openalex_author_id": normalize_space(detail.get("id", "")),
            "orcid": normalize_space(detail.get("orcid", "")),
            "institution_name": current_inst_name,
            "institution_openalex_id": current_inst_id,
            "institution_ror": inst_ror,
            "institution_type": inst_type,
            "city": city,
            "country_code": current_country,
            "institution_homepage": inst_homepage,
            "works_count_total": detail.get("works_count", ""),
            "cited_by_count_total": detail.get("cited_by_count", ""),
            "overall_works_count": detail.get("works_count", ""),
            "overall_cited_by_count": detail.get("cited_by_count", ""),
            "query_relevant_works_count": rec.get("query_relevant_works_count", 0),
            "query_total_citations_from_relevant_works": rec.get("query_total_citations_from_relevant_works", 0),
            "last_author_count": rec.get("last_author_count", 0),
            "first_author_count": rec.get("first_author_count", 0),
            "first_last_author_count": rec.get("first_last_author_count", 0),
            "matched_keyword_terms": "; ".join(sorted(rec.get("matched_keyword_terms", set()))),
            "top_topics": topic_text,
            "topic_share": topic_share_text,
            "lab_type": lab_type,
            "main_methods": main_methods,
            "main_topics": main_topics,
            "research_summary": research_summary,
            "publications_last_10y": len(works_last_10y),
            "publications_last_5y": len(works_last_5y),
            "publications_last_3y": len(works_last_3y),
            "keyword_relevant_papers_last_10y": len(relevant_last_10y),
            "keyword_relevant_papers_last_5y": len(relevant_last_5y),
            "keyword_relevant_papers_last_3y": len(relevant_last_3y),
            "publication_year_distribution_last_10y": stringify_year_distribution(year_dist, top_n=10),
            "representative_papers_overall_10y": stringify_works(top_cited_overall_10y, top_n=5),
            "representative_keyword_relevant_papers": stringify_works(top_relevant, top_n=5),
            "recent_papers_last_5y": stringify_works(top_recent_overall, top_n=5),
            "representative_papers_overall_10y_dois": ";".join(
                [normalize_space(w.get("doi", "")) for w in top_cited_overall_10y[:5] if normalize_space(w.get("doi", ""))]
            ),
            "representative_keyword_relevant_papers_dois": ";".join(
                [normalize_space(w.get("doi", "")) for w in top_relevant[:5] if normalize_space(w.get("doi", ""))]
            ),
            "recent_papers_last_5y_dois": ";".join(
                [normalize_space(w.get("doi", "")) for w in top_recent_overall[:5] if normalize_space(w.get("doi", ""))]
            ),
            "relevance_score": final_score,
            "source_summary": "OpenAlex works + authors + institutions",
        }
        rows.append(row)
        time.sleep(sleep_seconds)

    return rows


def save_outputs(rows: List[dict], output_tsv: str, output_jsonl: Optional[str]) -> None:
    ensure_parent_dir(output_tsv)
    df = pd.DataFrame(rows)

    preferred = [
        "query_keyword",
        "lab_name",
        "pi_name",
        "institution_name",
        "institution_type",
        "city",
        "country_code",
        "institution_homepage",
        "orcid",
        "lab_type",
        "research_summary",
        "main_methods",
        "main_topics",
        "matched_keyword_terms",
        "top_topics",
        "topic_share",
        "overall_cited_by_count",
        "overall_works_count",
        "publications_last_10y",
        "publications_last_5y",
        "publications_last_3y",
        "keyword_relevant_papers_last_10y",
        "keyword_relevant_papers_last_5y",
        "keyword_relevant_papers_last_3y",
        "publication_year_distribution_last_10y",
        "representative_papers_overall_10y",
        "representative_keyword_relevant_papers",
        "recent_papers_last_5y",
        "query_relevant_works_count",
        "query_total_citations_from_relevant_works",
        "first_last_author_count",
        "last_author_count",
        "first_author_count",
        "works_count_total",
        "cited_by_count_total",
        "relevance_score",
        "source_summary",
        "openalex_author_id",
        "institution_openalex_id",
        "institution_ror",
    ]
    cols = [c for c in preferred if c in df.columns] + [c for c in df.columns if c not in preferred]
    df = df[cols]
    df.to_csv(output_tsv, sep="\t", index=False)
    print(f"[saved] TSV -> {output_tsv} ({len(df)} rows)")

    if output_jsonl:
        ensure_parent_dir(output_jsonl)
        with open(output_jsonl, "w", encoding="utf-8") as f:
            for row in rows:
                f.write(json.dumps(row, ensure_ascii=False) + "\n")
        print(f"[saved] JSONL -> {output_jsonl}")


def main():
    parser = argparse.ArgumentParser(
        description="Search for big labs in a specific field using OpenAlex discovery, PI-like filtering, overall-impact ranking, and institution/paper de-duplication."
    )
    parser.add_argument("--keyword", required=True, help="Field keyword, e.g. CRISPR, bioinformatics, computational biology")
    parser.add_argument("--top-n", type=int, default=30, help="Number of final labs to keep")
    parser.add_argument("--max-works", type=int, default=300, help="Maximum relevant works to retrieve from OpenAlex")
    parser.add_argument("--per-page", type=int, default=100, help="Works per OpenAlex page")
    parser.add_argument("--years-back", type=int, default=5, help="How many recent years of works to use for initial discovery")
    parser.add_argument("--author-works-limit", type=int, default=150, help="How many recent works to retrieve per candidate PI for profiling")
    parser.add_argument("--sleep-seconds", type=float, default=0.2, help="Sleep between detail requests")
    parser.add_argument("--max-authors-per-paper", type=int, default=30, help="Exclude papers with more than this many authorships")
    parser.add_argument("--min-last-author-count", type=int, default=2, help="Minimum last-author count for PI-like filtering")
    parser.add_argument("--min-first-last-author-count", type=int, default=3, help="Minimum first/last-author count for PI-like filtering")
    parser.add_argument("--max-per-institution", type=int, default=5, help="Maximum final rows to keep per institution")
    parser.add_argument("--max-paper-overlap", type=int, default=5, help="Maximum allowed DOI overlap with already selected rows")
    parser.add_argument("--output", required=True, help="Output TSV path")
    parser.add_argument("--output-jsonl", default=None, help="Optional output JSONL path")
    args = parser.parse_args()

    keyword_terms = parse_keywords(args.keyword)
    print(f"[INFO] keyword terms: {keyword_terms}", file=sys.stderr)

    works = fetch_works_for_keyword(
        keyword=args.keyword,
        years_back=args.years_back,
        per_page=args.per_page,
        max_works=args.max_works,
        max_authors_per_paper=args.max_authors_per_paper,
    )
    print(f"[INFO] works retrieved: {len(works)}", file=sys.stderr)
    if not works:
        raise RuntimeError("No works retrieved from OpenAlex.")

    author_candidates = collect_author_candidates(works, keyword_terms=keyword_terms)
    print(f"[INFO] raw author candidates: {len(author_candidates)}", file=sys.stderr)
    if not author_candidates:
        raise RuntimeError("No author candidates found after aggregation.")

    author_candidates = filter_pi_like_candidates(
        author_candidates,
        min_last_author_count=args.min_last_author_count,
        min_first_last_author_count=args.min_first_last_author_count,
    )
    print(f"[INFO] PI-like author candidates: {len(author_candidates)}", file=sys.stderr)
    if not author_candidates:
        raise RuntimeError("No PI-like author candidates remained after filtering.")

    rows = build_lab_rows(
        author_candidates=author_candidates,
        keyword=args.keyword,
        top_n=args.top_n,
        sleep_seconds=args.sleep_seconds,
        author_works_limit=args.author_works_limit,
        max_authors_per_paper=args.max_authors_per_paper,
    )
    if not rows:
        raise RuntimeError("No lab rows built.")

    rows = cap_rows_per_institution(rows, max_per_institution=args.max_per_institution)
    rows = diversify_by_paper_overlap(rows, top_n=args.top_n, max_overlap=args.max_paper_overlap)

    if not rows:
        raise RuntimeError("No rows remained after institution and paper-overlap de-duplication.")

    save_outputs(rows, output_tsv=args.output, output_jsonl=args.output_jsonl)


if __name__ == "__main__":
    main()