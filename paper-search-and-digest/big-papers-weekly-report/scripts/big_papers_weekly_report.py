#!/usr/bin/env python3

import argparse
import datetime as dt
import html
import json
import os
import re
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Optional

import pandas as pd
import requests
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import (
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)


CROSSREF_WORKS = "https://api.crossref.org/works"
PUBMED_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

SESSION = requests.Session()
SESSION.headers.update(
    {"User-Agent": "bioinfor-claw/0.1 (bioinformatics-big-papers-weekly-report-fast)"}
)

DEFAULT_JOURNALS = [
    "Nature",
    "Science",
    "Cell",
    "Nature Biotechnology",
    "Nature Methods",
    "Nature Genetics",
    "Nature Medicine",
    "Nature Machine Intelligence",
    "Nature Computational Science",
    "Cell Systems",
    "Cell Genomics",
]

JOURNAL_WEIGHTS = {
    "Nature": 10.0,
    "Science": 10.0,
    "Cell": 10.0,
    "Nature Biotechnology": 9.5,
    "Nature Methods": 9.5,
    "Nature Genetics": 9.2,
    "Nature Medicine": 8.8,
    "Nature Machine Intelligence": 8.2,
    "Nature Computational Science": 8.4,
    "Cell Systems": 8.4,
    "Cell Genomics": 8.7,
}

GENERIC_BIOINFO_KEYWORDS = {
    "bioinformatics": 4.0,
    "alphafold": 4.0,
    "crispr": 4.0,
    "foundation model": 3.5,
    "language model": 3.5,
    "llm": 3.0,
    "benchmark": 3.0,
    "resource": 2.2,
    "atlas": 2.4,
    "database": 2.5,
    "single-cell": 3.0,
    "single cell": 2.4,
    "spatial": 2.1,
    "computational": 2.5,
    "algorithm": 2.2,
    "software": 2.5,
    "methods": 1.8,
    "multimodal": 2.5,
    "genomics": 2.1,
    "proteomics": 2.1,
    "transcriptomics": 2.1,
    "machine learning": 2.8,
    "deep learning": 2.8,
    "protein structure": 2.5,
}

NEGATIVE_KEYWORDS = {
    "case report": -4.0,
    "editorial": -5.0,
    "correction": -6.0,
    "erratum": -6.0,
    "retraction": -8.0,
    "news": -3.0,
    "perspective": -1.5,
    "comment": -2.0,
}


def ensure_parent_dir(path: str) -> None:
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)


def normalize_space(text: str) -> str:
    return re.sub(r"\s+", " ", str(text or "").strip())


def strip_jats_tags(text: str) -> str:
    if not text:
        return ""
    text = html.unescape(text)
    text = re.sub(r"<[^>]+>", " ", text)
    return normalize_space(text)


def safe_get_json(url: str, params: Optional[dict] = None, timeout: int = 60) -> dict:
    r = SESSION.get(url, params=params, timeout=timeout)
    r.raise_for_status()
    return r.json()


def safe_get_json_with_retry(
    url: str,
    params: Optional[dict] = None,
    timeout: int = 60,
    max_retries: int = 5,
) -> dict:
    last_err = None
    for attempt in range(max_retries):
        try:
            r = SESSION.get(url, params=params, timeout=timeout)
            if r.status_code == 429:
                wait = min(2 ** attempt, 30)
                print(f"[WARN] 429 from {url}; retrying in {wait}s")
                time.sleep(wait)
                continue
            r.raise_for_status()
            return r.json()
        except requests.HTTPError as e:
            last_err = e
            if getattr(e.response, "status_code", None) == 429:
                wait = min(2 ** attempt, 30)
                print(f"[WARN] 429 from {url}; retrying in {wait}s")
                time.sleep(wait)
                continue
            raise
        except Exception as e:
            last_err = e
            wait = min(2 ** attempt, 10)
            print(f"[WARN] Request failed; retrying in {wait}s: {e}")
            time.sleep(wait)
    raise last_err


def safe_get_text(url: str, params: Optional[dict] = None, timeout: int = 60) -> str:
    r = SESSION.get(url, params=params, timeout=timeout)
    r.raise_for_status()
    return r.text


def today_ymd() -> str:
    return dt.date.today().strftime("%Y-%m-%d")


def days_ago_ymd(days: int) -> str:
    return (dt.date.today() - dt.timedelta(days=days)).strftime("%Y-%m-%d")


def choose_strict_publication_date(item: dict) -> str:
    # Strict publication date preference:
    # published-print > published-online > published
    for key in ["published-print", "published-online", "published"]:
        part = item.get(key)
        if isinstance(part, dict):
            date_parts = part.get("date-parts", [])
            if date_parts and date_parts[0]:
                vals = list(date_parts[0]) + [1] * (3 - len(date_parts[0]))
                return f"{vals[0]:04d}-{vals[1]:02d}-{vals[2]:02d}"
    return ""


def filter_by_publication_window(df: pd.DataFrame, date_from: str, date_to: str) -> pd.DataFrame:
    out = df.copy()
    out["pub_date"] = pd.to_datetime(out["pub_date"], errors="coerce")
    start = pd.to_datetime(date_from)
    end = pd.to_datetime(date_to)

    out = out[out["pub_date"].notna()].copy()
    out = out[(out["pub_date"] >= start) & (out["pub_date"] <= end)].copy()
    out["pub_date"] = out["pub_date"].dt.strftime("%Y-%m-%d")
    return out


def split_sentences(text: str) -> List[str]:
    text = normalize_space(text)
    if not text:
        return []
    return [x.strip() for x in re.split(r"(?<=[.!?])\s+", text) if x.strip()]


def sentence_excerpt(text: str, max_sentences: int = 4, max_chars: int = 2200) -> str:
    text = normalize_space(text)
    if not text:
        return ""
    sents = split_sentences(text)
    out = []
    total = 0
    for s in sents:
        if len(out) >= max_sentences:
            break
        if total + len(s) > max_chars:
            break
        out.append(s)
        total += len(s) + 1
    if out:
        return " ".join(out)
    return text[:max_chars].rstrip()


def load_interest_keywords(keyword_str: str, keyword_file: Optional[str]) -> List[str]:
    kws = []
    if keyword_str:
        kws.extend([x.strip() for x in keyword_str.split(",") if x.strip()])
    if keyword_file:
        with open(keyword_file) as f:
            kws.extend([x.strip() for x in f if x.strip()])
    return list(dict.fromkeys([k.lower() for k in kws if k.strip()]))


def load_journals(
    journal_str: str,
    journal_file: Optional[str],
    replace_default: bool,
) -> List[str]:
    user_journals = []
    if journal_str:
        user_journals.extend([x.strip() for x in journal_str.split(",") if x.strip()])
    if journal_file:
        with open(journal_file) as f:
            user_journals.extend([x.strip() for x in f if x.strip()])
    user_journals = list(dict.fromkeys(user_journals))

    journals = user_journals if replace_default else DEFAULT_JOURNALS + user_journals
    journals = list(dict.fromkeys([j for j in journals if j.strip()]))

    if not journals:
        raise ValueError("No journals available after parsing defaults and user input.")
    return journals


def fetch_crossref_for_journal(
    journal: str,
    date_from: str,
    date_to: str,
    rows: int = 20,
) -> List[dict]:
    params = {
        "filter": f"from-pub-date:{date_from},until-pub-date:{date_to},container-title:{journal},type:journal-article",
        "rows": rows,
        "select": ",".join(
            [
                "DOI",
                "title",
                "container-title",
                "published",
                "published-print",
                "published-online",
                "abstract",
                "author",
                "URL",
                "subject",
                "is-referenced-by-count",
                "type",
            ]
        ),
    }
    data = safe_get_json_with_retry(CROSSREF_WORKS, params=params)
    return data.get("message", {}).get("items", [])


def parse_crossref_item(item: dict) -> dict:
    title = ""
    titles = item.get("title", [])
    if isinstance(titles, list) and titles:
        title = normalize_space(titles[0])

    journal = ""
    containers = item.get("container-title", [])
    if isinstance(containers, list) and containers:
        journal = normalize_space(containers[0])

    abstract = strip_jats_tags(item.get("abstract", ""))

    authors = []
    for a in item.get("author", []) or []:
        given = normalize_space(a.get("given", ""))
        family = normalize_space(a.get("family", ""))
        name = " ".join([x for x in [given, family] if x])
        if name:
            authors.append(name)

    subjects = item.get("subject", []) or []
    if not isinstance(subjects, list):
        subjects = []

    return {
        "doi": normalize_space(item.get("DOI", "")),
        "title": title,
        "journal": journal,
        "pub_date": choose_strict_publication_date(item),
        "abstract_crossref": abstract,
        "authors": "; ".join(authors[:12]),
        "url": normalize_space(item.get("URL", "")),
        "subjects": "; ".join([normalize_space(x) for x in subjects if x]),
        "citation_count": item.get("is-referenced-by-count", 0) or 0,
        "type": normalize_space(item.get("type", "")),
    }


def fetch_crossref_records_parallel(
    journals: List[str],
    date_from: str,
    date_to: str,
    rows_per_journal: int,
    max_workers: int = 2,
) -> List[dict]:
    records = []
    seen_doi = set()

    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        future_map = {
            ex.submit(fetch_crossref_for_journal, journal, date_from, date_to, rows_per_journal): journal
            for journal in journals
        }
        for fut in as_completed(future_map):
            journal = future_map[fut]
            try:
                items = fut.result()
                print(f"[INFO] Retrieved {len(items)} records from Crossref for: {journal}")
            except Exception as e:
                print(f"[WARN] Crossref query failed for {journal}: {e}")
                continue

            for item in items:
                rec = parse_crossref_item(item)
                doi = rec["doi"].lower()
                if not doi or doi in seen_doi:
                    continue
                seen_doi.add(doi)
                records.append(rec)

    return records


def pubmed_search_by_doi(doi: str, email: str, api_key: Optional[str] = None) -> List[str]:
    params = {
        "db": "pubmed",
        "term": f"{doi}[AID]",
        "retmode": "json",
        "retmax": 3,
        "tool": "bioinfor_claw",
        "email": email,
    }
    if api_key:
        params["api_key"] = api_key
    data = safe_get_json(PUBMED_ESEARCH, params=params)
    return data.get("esearchresult", {}).get("idlist", [])


def pubmed_fetch_xml(pmids: List[str], email: str, api_key: Optional[str] = None) -> str:
    if not pmids:
        return ""
    params = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "xml",
        "tool": "bioinfor_claw",
        "email": email,
    }
    if api_key:
        params["api_key"] = api_key
    return safe_get_text(PUBMED_EFETCH, params=params)


def parse_pubmed_abstract_from_xml(xml_text: str) -> str:
    if not xml_text.strip():
        return ""
    root = ET.fromstring(xml_text)
    abstracts = []
    for abstract in root.findall(".//Abstract/AbstractText"):
        label = abstract.attrib.get("Label", "")
        txt = "".join(abstract.itertext()).strip()
        if txt:
            abstracts.append(f"{label}: {txt}" if label else txt)
    return normalize_space(" ".join(abstracts))


def fetch_pubmed_abstract_by_doi(doi: str, email: str, api_key: Optional[str] = None) -> str:
    try:
        pmids = pubmed_search_by_doi(doi, email=email, api_key=api_key)
        if not pmids:
            return ""
        time.sleep(0.34 if not api_key else 0.12)
        xml_text = pubmed_fetch_xml(pmids, email=email, api_key=api_key)
        return parse_pubmed_abstract_from_xml(xml_text)
    except Exception:
        return ""


def load_cache(cache_file: str) -> Dict[str, str]:
    if not cache_file or not os.path.exists(cache_file):
        return {}
    try:
        with open(cache_file) as f:
            data = json.load(f)
        if isinstance(data, dict):
            return {str(k): str(v) for k, v in data.items()}
    except Exception:
        pass
    return {}


def save_cache(cache: Dict[str, str], cache_file: str) -> None:
    ensure_parent_dir(cache_file)
    with open(cache_file, "w") as f:
        json.dump(cache, f, indent=2, sort_keys=True)


def enrich_abstracts_with_pubmed(
    df: pd.DataFrame,
    email: str,
    api_key: Optional[str],
    cache: Dict[str, str],
) -> tuple[pd.DataFrame, bool]:
    abstracts = []
    updated = False

    for _, row in df.iterrows():
        doi = str(row["doi"]).lower().strip()
        crossref_abs = normalize_space(row.get("abstract_crossref", "") or "")

        if doi in cache:
            abstracts.append(cache[doi])
            continue

        abstract = ""
        if doi:
            abstract = fetch_pubmed_abstract_by_doi(doi=doi, email=email, api_key=api_key)

        if not abstract:
            abstract = crossref_abs

        abstract = normalize_space(abstract)
        cache[doi] = abstract
        abstracts.append(abstract)
        updated = True

    df = df.copy()
    df["abstract"] = abstracts
    return df, updated


def generic_bioinformatics_score(title: str, abstract: str, subjects: str) -> float:
    text = f"{title} {abstract} {subjects}".lower()
    score = 0.0
    for kw, val in GENERIC_BIOINFO_KEYWORDS.items():
        if kw in text:
            score += val
    for kw, val in NEGATIVE_KEYWORDS.items():
        if kw in text:
            score += val
    return score


def match_generic_keywords(title: str, abstract: str, subjects: str) -> Dict[str, object]:
    text = f"{title} {abstract} {subjects}".lower()
    matched = []

    for kw in GENERIC_BIOINFO_KEYWORDS:
        if kw in text and kw not in matched:
            matched.append(kw)

    return {
        "matched_generic_keywords": "; ".join(matched),
        "generic_keyword_count": float(len(matched)),
    }


def journal_score(journal: str) -> float:
    return JOURNAL_WEIGHTS.get(journal, 5.0)


def match_interest_keywords(title: str, abstract: str, subjects: str, interest_keywords: List[str]) -> Dict[str, object]:
    text = f"{title} {abstract} {subjects}".lower()
    matched = []
    for kw in interest_keywords:
        if kw in text and kw not in matched:
            matched.append(kw)
    return {
        "matched_interest_keywords": "; ".join(matched),
        "interest_keyword_score": float(len(matched)),
    }


def compute_score_components(
    journal: str,
    title: str,
    abstract: str,
    subjects: str,
    citation_count: float,
    interest_keyword_score: float,
    generic_bioinfo_weight: float,
    interest_weight: float,
) -> Dict[str, float]:
    jscore = journal_score(journal)
    gscore = generic_bioinformatics_score(title, abstract, subjects)
    citation_bonus = min(float(citation_count or 0), 50.0) * 0.03
    abstract_bonus = 1.0 if abstract else 0.0

    bonus = 0.0
    tlow = title.lower()
    if "benchmark" in tlow:
        bonus += 1.5
    if "atlas" in tlow:
        bonus += 1.0
    if "resource" in tlow:
        bonus += 1.0

    total = (
        jscore
        + gscore * generic_bioinfo_weight
        + interest_keyword_score * interest_weight
        + citation_bonus
        + abstract_bonus
        + bonus
    )

    return {
        "journal_score": round(jscore, 3),
        "generic_bioinfo_score": round(gscore, 3),
        "citation_bonus": round(citation_bonus, 3),
        "abstract_bonus": round(abstract_bonus, 3),
        "bonus_score": round(bonus, 3),
        "impact_score": round(total, 3),
    }


def add_scoring_columns(
    df: pd.DataFrame,
    interest_keywords: List[str],
    generic_bioinfo_weight: float,
    interest_weight: float,
) -> pd.DataFrame:
    matched_kw_col = []
    matched_generic_kw_col = []
    interest_score_col = []
    generic_keyword_count_col = []
    journal_score_col = []
    generic_score_col = []
    citation_bonus_col = []
    abstract_bonus_col = []
    bonus_score_col = []
    impact_score_col = []

    for _, row in df.iterrows():
        title = row["title"]
        abstract = row.get("abstract", "") or row.get("abstract_crossref", "")
        subjects = row["subjects"]

        match_info = match_interest_keywords(title, abstract, subjects, interest_keywords)
        generic_match_info = match_generic_keywords(title, abstract, subjects)

        score_info = compute_score_components(
            journal=row["journal"],
            title=title,
            abstract=abstract,
            subjects=subjects,
            citation_count=row["citation_count"],
            interest_keyword_score=match_info["interest_keyword_score"],
            generic_bioinfo_weight=generic_bioinfo_weight,
            interest_weight=interest_weight,
        )

        matched_kw_col.append(match_info["matched_interest_keywords"])
        matched_generic_kw_col.append(generic_match_info["matched_generic_keywords"])
        interest_score_col.append(match_info["interest_keyword_score"])
        generic_keyword_count_col.append(generic_match_info["generic_keyword_count"])
        journal_score_col.append(score_info["journal_score"])
        generic_score_col.append(score_info["generic_bioinfo_score"])
        citation_bonus_col.append(score_info["citation_bonus"])
        abstract_bonus_col.append(score_info["abstract_bonus"])
        bonus_score_col.append(score_info["bonus_score"])
        impact_score_col.append(score_info["impact_score"])

    out = df.copy()
    out["matched_interest_keywords"] = matched_kw_col
    out["matched_generic_keywords"] = matched_generic_kw_col

    has_user_interest = len(interest_keywords) > 0
    if has_user_interest:
        out["matched_keywords_display"] = matched_kw_col
        out["matched_user_interest"] = [bool(str(x).strip()) for x in matched_kw_col]
    else:
        out["matched_keywords_display"] = matched_generic_kw_col
        out["matched_user_interest"] = [False] * len(out)

    out["interest_keyword_score"] = interest_score_col
    out["generic_keyword_count"] = generic_keyword_count_col
    out["journal_score"] = journal_score_col
    out["generic_bioinfo_score"] = generic_score_col
    out["citation_bonus"] = citation_bonus_col
    out["abstract_bonus"] = abstract_bonus_col
    out["bonus_score"] = bonus_score_col
    out["impact_score"] = impact_score_col
    return out


def impact_score_method_text(
    interest_keywords: List[str],
    generic_bioinfo_weight: float,
    interest_weight: float,
    fast_mode: bool,
    pubmed_top_k: int,
    skip_pubmed: bool,
) -> List[str]:
    kw_text = ", ".join(interest_keywords) if interest_keywords else "None provided"
    mode_text = (
        "Fast mode was used: the script first ranked papers using Crossref metadata and only enriched the top candidate set with PubMed abstracts."
        if fast_mode and not skip_pubmed
        else "PubMed enrichment was skipped; ranking used Crossref metadata only."
        if skip_pubmed
        else "Full mode was used."
    )
    return [
        "Impact score is a heuristic weekly triage score, not a citation-based impact metric.",
        (
            "Score = journal prestige score + "
            f"generic bioinformatics relevance score × {generic_bioinfo_weight:.2f} + "
            f"user-interest keyword score × {interest_weight:.2f} + "
            "small citation-count bonus + abstract-availability bonus + small bonus for benchmark/resource/atlas-style papers."
        ),
        "Only papers whose strict publication date falls inside the requested date window are retained. Publication date is taken in this priority order: published-print, then published-online, then published.",
        "Journal prestige score is manually assigned by journal title.",
        "Generic bioinformatics relevance is estimated from broad field-related keywords in title, abstract, and subject metadata.",
        f"User-interest keyword relevance is computed from matched user-supplied keywords. Keywords used in this run: {kw_text}.",
        "If user-interest keywords are supplied, only papers matching at least one user-interest keyword are retained in the report.",
        "When user-interest keywords are supplied, the report displays only matched user-interest keywords. Built-in generic keywords are displayed only when no user-interest keywords are provided.",
        "Each paper includes a matched-keywords column so the ranking can be interpreted transparently.",
        mode_text,
        f"If PubMed enrichment was used, only the top {pubmed_top_k} candidates from the first-pass ranking were enriched with PubMed abstracts.",
        "This score is intended to prioritize likely high-interest papers for manual review.",
    ]


def summarize_paper(title: str, abstract: str, journal: str, matched_keywords: str) -> Dict[str, str]:
    abs_sents = split_sentences(abstract)

    if abs_sents:
        why_it_matters = abs_sents[0]
        major_discoveries = " ".join(abs_sents[1:3]).strip() if len(abs_sents) > 1 else abs_sents[0]
    else:
        why_it_matters = (
            f"This paper appears to present a potentially high-impact contribution in "
            f"bioinformatics or computational biology published in {journal}."
        )
        major_discoveries = "Abstract unavailable from the retrieval sources used by this script."

    text = f"{title} {abstract}".lower()
    impact_bits = []

    if matched_keywords:
        impact_bits.append(f"Matched keywords: {matched_keywords}.")
    if "benchmark" in text:
        impact_bits.append("Provides a benchmark or evaluation framework.")
    if "resource" in text or "database" in text:
        impact_bits.append("Introduces a reusable community resource.")
    if "single-cell" in text or "single cell" in text:
        impact_bits.append("Likely relevant to single-cell analysis workflows.")
    if "foundation model" in text or "language model" in text or "llm" in text:
        impact_bits.append("Likely important for AI-driven bioinformatics.")
    if "atlas" in text:
        impact_bits.append("Builds an atlas-scale reference or dataset.")
    if "method" in text or "algorithm" in text or "software" in text or "tool" in text:
        impact_bits.append("Likely offers methodological or software impact.")
    if "spatial" in text:
        impact_bits.append("Potentially important for spatial omics analysis.")
    if "protein structure" in text or "alphafold" in text:
        impact_bits.append("Potentially important for structural bioinformatics.")

    if not impact_bits:
        impact_bits.append(
            "Potentially important because it appears to advance a broadly useful computational or quantitative approach."
        )

    return {
        "why_it_matters": normalize_space(why_it_matters),
        "major_discoveries": normalize_space(major_discoveries),
        "significance": normalize_space(" ".join(dict.fromkeys(impact_bits))),
    }

def is_non_primary_article(rec: dict) -> bool:
    text = " ".join([
        str(rec.get("title", "")),
        str(rec.get("type", "")),
        str(rec.get("abstract_crossref", "")),
    ]).lower()

    bad_terms = [
        "erratum",
        "correction",
        "corrigendum",
        "retraction",
        "editorial",
        "comment",
        "perspective",
        "news",
    ]
    return any(term in text for term in bad_terms)


def build_pdf_report(
    df: pd.DataFrame,
    out_pdf: str,
    title_text: str,
    date_from: str,
    date_to: str,
    interest_keywords: List[str],
    generic_bioinfo_weight: float,
    interest_weight: float,
    journals: List[str],
    fast_mode: bool,
    pubmed_top_k: int,
    skip_pubmed: bool,
) -> None:
    ensure_parent_dir(out_pdf)

    doc = SimpleDocTemplate(
        out_pdf,
        pagesize=letter,
        leftMargin=0.7 * inch,
        rightMargin=0.7 * inch,
        topMargin=0.7 * inch,
        bottomMargin=0.7 * inch,
    )

    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="TitleCenter", parent=styles["Title"], alignment=TA_CENTER, spaceAfter=12))
    styles.add(
        ParagraphStyle(
            name="SectionHead",
            parent=styles["Heading2"],
            alignment=TA_LEFT,
            textColor=colors.HexColor("#1f3b73"),
            spaceBefore=8,
            spaceAfter=6,
        )
    )
    styles.add(
        ParagraphStyle(
            name="PaperTitle",
            parent=styles["Heading3"],
            alignment=TA_LEFT,
            textColor=colors.black,
            spaceBefore=8,
            spaceAfter=4,
        )
    )
    styles.add(ParagraphStyle(name="BodySmall", parent=styles["BodyText"], fontSize=9.5, leading=12, spaceAfter=6))
    styles.add(
        ParagraphStyle(
            name="AbstractStyle",
            parent=styles["BodyText"],
            fontSize=9,
            leading=11.5,
            spaceAfter=6,
            textColor=colors.HexColor("#333333"),
        )
    )

    story = []
    story.append(Paragraph(title_text, styles["TitleCenter"]))
    story.append(Paragraph(f"Coverage window: {date_from} to {date_to}", styles["BodySmall"]))
    story.append(Paragraph(f"Journal scope: {html.escape('; '.join(journals))}", styles["BodySmall"]))
    story.append(Spacer(1, 0.12 * inch))

    story.append(Paragraph("Executive Summary", styles["SectionHead"]))
    story.append(
        Paragraph(
            f"This report summarizes {len(df)} recent papers selected from the configured journal scope and ranked for likely bioinformatics relevance, user-interest relevance, and potential field impact.",
            styles["BodySmall"],
        )
    )

    story.append(Paragraph("Impact Score Method", styles["SectionHead"]))
    for line in impact_score_method_text(
        interest_keywords, generic_bioinfo_weight, interest_weight, fast_mode, pubmed_top_k, skip_pubmed
    ):
        story.append(Paragraph("• " + html.escape(line), styles["BodySmall"]))

    story.append(Spacer(1, 0.1 * inch))
    story.append(Paragraph("Top Papers This Week", styles["SectionHead"]))

    table_data = [[
        Paragraph("<b>Rank</b>", styles["BodySmall"]),
        Paragraph("<b>Title</b>", styles["BodySmall"]),
        Paragraph("<b>Journal</b>", styles["BodySmall"]),
        Paragraph("<b>Date</b>", styles["BodySmall"]),
        Paragraph("<b>Matched keywords</b>", styles["BodySmall"]),
        Paragraph("<b>Impact score</b>", styles["BodySmall"]),
    ]]

    for i, row in df.reset_index(drop=True).iterrows():
        table_data.append(
            [
                Paragraph(str(i + 1), styles["BodySmall"]),
                Paragraph(html.escape(str(row["title"])), styles["BodySmall"]),
                Paragraph(html.escape(str(row["journal"])), styles["BodySmall"]),
                Paragraph(html.escape(str(row["pub_date"])), styles["BodySmall"]),
                Paragraph(html.escape(str(row["matched_keywords_display"] or "")), styles["BodySmall"]),
                Paragraph(f'{row["impact_score"]:.2f}', styles["BodySmall"]),
            ]
        )

    tbl = Table(
        table_data,
        colWidths=[0.42 * inch, 3.15 * inch, 1.15 * inch, 0.72 * inch, 1.35 * inch, 0.65 * inch],
        repeatRows=1,
    )
    tbl.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#dbe7ff")),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.black),
                ("GRID", (0, 0), (-1, -1), 0.35, colors.grey),
                ("FONTSIZE", (0, 0), (-1, -1), 8.2),
                ("LEADING", (0, 0), (-1, -1), 10),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                ("TOPPADDING", (0, 0), (-1, -1), 3),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
            ]
        )
    )
    story.append(tbl)
    story.append(PageBreak())

    for idx, row in df.reset_index(drop=True).iterrows():
        story.append(Paragraph(f"{idx+1}. {html.escape(row['title'])}", styles["PaperTitle"]))
        meta = (
            f"<b>Journal:</b> {html.escape(row['journal'])} &nbsp;&nbsp; "
            f"<b>Date:</b> {html.escape(row['pub_date'])} &nbsp;&nbsp; "
            f"<b>DOI:</b> {html.escape(row['doi'])}"
        )
        story.append(Paragraph(meta, styles["BodySmall"]))

        if row.get("authors"):
            story.append(Paragraph(f"<b>Authors:</b> {html.escape(str(row['authors']))}", styles["BodySmall"]))

        if row.get("matched_keywords_display"):
            if len(interest_keywords) > 0:
                label = "Matched interest keywords"
            else:
                label = "Matched generic keywords"
            story.append(
                Paragraph(
                    f"<b>{label}:</b> {html.escape(str(row['matched_keywords_display']))}",
                    styles["BodySmall"],
                )
            )

        score_line = (
            f"<b>Impact score:</b> {row['impact_score']:.2f} "
            f"(journal={row['journal_score']:.2f}, "
            f"generic={row['generic_bioinfo_score']:.2f}, "
            f"interest={row['interest_keyword_score']:.2f}, "
            f"citation bonus={row['citation_bonus']:.2f}, "
            f"abstract bonus={row['abstract_bonus']:.2f}, "
            f"extra bonus={row['bonus_score']:.2f})"
        )
        story.append(Paragraph(score_line, styles["BodySmall"]))

        story.append(Paragraph(f"<b>Why it matters:</b> {html.escape(str(row['why_it_matters']))}", styles["BodySmall"]))
        story.append(Paragraph(f"<b>Major discoveries:</b> {html.escape(str(row['major_discoveries']))}", styles["BodySmall"]))
        story.append(Paragraph(f"<b>Likely impact:</b> {html.escape(str(row['significance']))}", styles["BodySmall"]))

        if row.get("abstract"):
            abs_text = sentence_excerpt(str(row["abstract"]), max_sentences=4, max_chars=2200)
            story.append(
                Paragraph(
                    f"<b>Abstract excerpt:</b> {html.escape(abs_text)}",
                    styles["AbstractStyle"],
                )
            )

        story.append(Spacer(1, 0.12 * inch))

    doc.build(story)
  


def main():
    parser = argparse.ArgumentParser(
        description="Fast weekly report for big recent bioinformatics papers from selected journals."
    )
    parser.add_argument("--date-from", default=days_ago_ymd(7), help="Start date YYYY-MM-DD")
    parser.add_argument("--date-to", default=today_ymd(), help="End date YYYY-MM-DD")
    parser.add_argument("--top-n", type=int, default=20, help="Number of papers to include in final report")
    parser.add_argument("--rows-per-journal", type=int, default=30, help="Crossref rows per journal query")
    parser.add_argument("--email", default="your_email@example.com", help="Email for PubMed E-utilities")
    parser.add_argument("--ncbi-api-key", default=None, help="Optional NCBI API key")
    parser.add_argument("--output-prefix", required=True, help="Output prefix for TSV and PDF")

    parser.add_argument(
        "--journals",
        default="",
        help="Comma-separated journal list, e.g. 'Nature,Nature Methods,Cell Systems'",
    )
    parser.add_argument(
        "--journal-list-file",
        default=None,
        help="Optional text file with one journal per line",
    )
    parser.add_argument(
        "--replace-default-journals",
        action="store_true",
        help="Use only user-supplied journals instead of appending to default journals",
    )

    parser.add_argument(
        "--interest-keywords",
        default="",
        help="Comma-separated user-interest keywords, e.g. single-cell,spatial,foundation model",
    )
    parser.add_argument(
        "--interest-keywords-file",
        default=None,
        help="Optional text file with one interest keyword per line",
    )
    parser.add_argument(
        "--generic-bioinfo-weight",
        type=float,
        default=1.0,
        help="Weight for generic bioinformatics relevance score",
    )
    parser.add_argument(
        "--interest-weight",
        type=float,
        default=2.0,
        help="Weight for user-interest keyword relevance score",
    )

    parser.add_argument(
        "--skip-pubmed",
        action="store_true",
        help="Skip PubMed enrichment and use Crossref abstracts only",
    )
    parser.add_argument(
        "--pubmed-top-k",
        type=int,
        default=20,
        help="Only enrich top K first-pass candidates with PubMed abstracts",
    )
    parser.add_argument(
        "--cache-file",
        default="cache/pubmed_abstract_cache.json",
        help="JSON cache file for DOI -> abstract",
    )
    parser.add_argument(
        "--crossref-workers",
        type=int,
        default=2,
        help="Number of parallel workers for Crossref journal retrieval",
    )
    args = parser.parse_args()

    journals = load_journals(
        journal_str=args.journals,
        journal_file=args.journal_list_file,
        replace_default=args.replace_default_journals,
    )
    interest_keywords = load_interest_keywords(args.interest_keywords, args.interest_keywords_file)

    records = fetch_crossref_records_parallel(
        journals=journals,
        date_from=args.date_from,
        date_to=args.date_to,
        rows_per_journal=args.rows_per_journal,
        max_workers=args.crossref_workers,
    )
    if not records:
        raise RuntimeError("No records retrieved from Crossref for the selected journals and dates.")

    df = pd.DataFrame(records)
    
    df = df[~df.apply(is_non_primary_article, axis=1)].copy()

    # Strict publication-date filtering
    df = filter_by_publication_window(df, args.date_from, args.date_to)
    if len(df) == 0:
        raise RuntimeError("No papers with publication dates inside the requested date window were found.")

    # First pass: Crossref only
    df["abstract"] = df["abstract_crossref"].fillna("").map(normalize_space)
    df = add_scoring_columns(
        df,
        interest_keywords=interest_keywords,
        generic_bioinfo_weight=args.generic_bioinfo_weight,
        interest_weight=args.interest_weight,
    )

    generic_signal = [
        generic_bioinformatics_score(t, a, s)
        for t, a, s in zip(df["title"], df["abstract"], df["subjects"])
    ]
    df["generic_signal_tmp"] = generic_signal
    df["keep_candidate"] = (df["generic_signal_tmp"] >= 2.0) | (df["impact_score"] >= 11.0)

    df = df[df["keep_candidate"]].copy()

    if len(interest_keywords) > 0:
        df = df[df["matched_user_interest"]].copy()

    if len(df) == 0:
        raise RuntimeError("No papers remained after filtering. Try broader journals or broader interest keywords.")

    df = df.sort_values(
        ["impact_score", "citation_count", "journal", "pub_date"],
        ascending=[False, False, True, False],
    ).reset_index(drop=True)

    # Fast PubMed enrichment only for top K
    if not args.skip_pubmed:
        top_k = min(args.pubmed_top_k, len(df))
        cache = load_cache(args.cache_file)

        df_top = df.head(top_k).copy()
        df_rest = df.iloc[top_k:].copy()

        df_top, updated = enrich_abstracts_with_pubmed(
            df=df_top,
            email=args.email,
            api_key=args.ncbi_api_key,
            cache=cache,
        )
        if updated:
            save_cache(cache, args.cache_file)

        df_top = add_scoring_columns(
            df_top,
            interest_keywords=interest_keywords,
            generic_bioinfo_weight=args.generic_bioinfo_weight,
            interest_weight=args.interest_weight,
        )

        if len(interest_keywords) > 0:
            df_top = df_top[df_top["matched_user_interest"]].copy()

        df = pd.concat([df_top, df_rest], axis=0).reset_index(drop=True)

    if len(df) == 0:
        raise RuntimeError("No papers remained after PubMed enrichment and keyword filtering.")

    why_list = []
    disc_list = []
    sig_list = []

    for _, row in df.iterrows():
        summ = summarize_paper(
            title=row["title"],
            abstract=row.get("abstract", "") or row.get("abstract_crossref", ""),
            journal=row["journal"],
            matched_keywords=row["matched_keywords_display"],
        )
        why_list.append(summ["why_it_matters"])
        disc_list.append(summ["major_discoveries"])
        sig_list.append(summ["significance"])

    df["why_it_matters"] = why_list
    df["major_discoveries"] = disc_list
    df["significance"] = sig_list
    
    df = df.sort_values(
        ["impact_score", "citation_count", "journal", "pub_date"],
        ascending=[False, False, True, False],
    ).head(args.top_n).reset_index(drop=True)

    out_tsv = f"{args.output_prefix}.tsv"
    out_pdf = f"{args.output_prefix}.pdf"

    ensure_parent_dir(out_tsv)
    df.to_csv(out_tsv, sep="\t", index=False)

    build_pdf_report(
        df=df,
        out_pdf=out_pdf,
        title_text="Weekly Big Bioinformatics Papers Report",
        date_from=args.date_from,
        date_to=args.date_to,
        interest_keywords=interest_keywords,
        generic_bioinfo_weight=args.generic_bioinfo_weight,
        interest_weight=args.interest_weight,
        journals=journals,
        fast_mode=True,
        pubmed_top_k=args.pubmed_top_k,
        skip_pubmed=args.skip_pubmed,
    )

    print(f"[INFO] Saved paper table: {out_tsv}")
    print(f"[INFO] Saved PDF report: {out_pdf}")
    print(f"[INFO] Final report size: {len(df)} papers")
    if len(interest_keywords) > 0:
        print("[INFO] User-interest keyword filtering was applied.")
    print("[INFO] Strict publication-date filtering was applied.")
    
if __name__ == "__main__":
    main()