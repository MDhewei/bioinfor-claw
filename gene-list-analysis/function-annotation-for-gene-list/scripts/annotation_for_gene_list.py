#!/usr/bin/env python3

import argparse
import os
import re
import time
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
import requests


MYGENE_QUERY_URL = "https://mygene.info/v3/query"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"

SESSION = requests.Session()
SESSION.headers.update({
    "User-Agent": "bioinfor-claw/0.1 (annotation-for-gene-list)"
})


ORGANISM_MAP = {
    "human": {"mygene": "human", "taxid": 9606, "uniprot_organism_id": 9606, "label": "human"},
    "homo sapiens": {"mygene": "human", "taxid": 9606, "uniprot_organism_id": 9606, "label": "human"},
    "hs": {"mygene": "human", "taxid": 9606, "uniprot_organism_id": 9606, "label": "human"},
    "hsapiens": {"mygene": "human", "taxid": 9606, "uniprot_organism_id": 9606, "label": "human"},
    "mouse": {"mygene": "mouse", "taxid": 10090, "uniprot_organism_id": 10090, "label": "mouse"},
    "mus musculus": {"mygene": "mouse", "taxid": 10090, "uniprot_organism_id": 10090, "label": "mouse"},
    "mm": {"mygene": "mouse", "taxid": 10090, "uniprot_organism_id": 10090, "label": "mouse"},
}


OUTPUT_SCHEMA = [
    "input_gene",
    "mapping_status",
    "mapping_note",
    "approved_symbol",
    "full_gene_name",
    "aliases",
    "organism",
    "entrez_id",
    "ensembl_gene_id",
    "uniprot_accession",
    "protein_name",
    "protein_length",
    "gene_type",
    "chromosome",
    "cytoband",
    "strand",
    "genomic_location",
    "reported_function",
    "subcellular_location",
    "domain_annotation",
    "domain_coordinates",
    "key_domain_summary",
    "functional_domain_class",
    "molecular_function",
    "biological_process",
    "cellular_component",
    "reactome_pathways",
    "kegg_pathways",
    "hallmark_pathways",
    "pathway_annotation_summary",
    "disease_association",
    "keyword_annotation",
    "protein_class",
    "druggability",
    "source_summary",
]


def normalize_organism(org: str) -> Dict[str, Any]:
    key = str(org).strip().lower()
    if key not in ORGANISM_MAP:
        raise ValueError(
            f"Unsupported organism '{org}'. Supported examples: human, mouse, hsapiens, mm"
        )
    return ORGANISM_MAP[key]


def ensure_parent_dir(path: str) -> None:
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)


def clean_gene(g: str) -> str:
    g = str(g).strip()
    g = re.sub(r"\s+", "", g)
    return g


def read_gene_list(path: str) -> List[str]:
    genes = []
    with open(path) as f:
        for line in f:
            x = line.strip()
            if not x:
                continue
            x = re.split(r"[\t,; ]+", x)[0]
            if x:
                genes.append(clean_gene(x))
    genes = list(dict.fromkeys(genes))
    if not genes:
        raise ValueError("Input gene list is empty")
    return genes


def safe_get_json(url: str, params: Optional[dict] = None, timeout: int = 60) -> dict:
    r = SESSION.get(url, params=params, timeout=timeout)
    r.raise_for_status()
    return r.json()


def _as_list(x):
    if x is None:
        return []
    if isinstance(x, list):
        return x
    return [x]


def _join_unique(values, sep="; "):
    out = []
    seen = set()
    for v in values:
        if v is None:
            continue
        s = str(v).strip()
        if not s:
            continue
        if s not in seen:
            seen.add(s)
            out.append(s)
    return sep.join(out)


def summarize_text(text: str, max_len: int = 300) -> str:
    text = re.sub(r"\s+", " ", str(text or "").strip())
    if len(text) <= max_len:
        return text
    return text[: max_len - 3].rstrip() + "..."


def clean_sentence(text: str) -> str:
    text = str(text or "").strip()
    text = re.sub(r"\s+", " ", text)
    text = re.sub(r"\s*PubMed:\d+\.?", "", text)
    text = re.sub(r"\s*\[.*?\]", "", text)
    text = text.strip(" ;|")
    return text


def split_sentences(text: str) -> List[str]:
    text = clean_sentence(text)
    if not text:
        return []
    parts = re.split(r"(?<=[.!?])\s+", text)
    return [p.strip() for p in parts if p.strip()]


def concise_function_summary(mygene_summary: str, uniprot_function: str, max_len: int = 220) -> str:
    candidates = []

    for source_text in [mygene_summary, uniprot_function]:
        for sent in split_sentences(source_text):
            sent = clean_sentence(sent)
            if len(sent) < 25:
                continue
            if sent not in candidates:
                candidates.append(sent)

    if not candidates:
        return ""

    summary = candidates[0]

    if len(summary) < 120 and len(candidates) > 1:
        s2 = candidates[1]
        w1 = set(re.findall(r"[A-Za-z0-9\-]+", summary.lower()))
        w2 = set(re.findall(r"[A-Za-z0-9\-]+", s2.lower()))
        overlap = len(w1 & w2) / max(1, len(w2))
        if overlap < 0.6 and s2:
            summary = summary.rstrip(".") + "; " + s2[0].lower() + s2[1:]

    if len(summary) > max_len:
        summary = summary[:max_len - 3].rstrip() + "..."

    return summary


def augment_function_with_domain(summary: str, domain_class: str) -> str:
    summary = clean_sentence(summary)
    domain_class = clean_sentence(domain_class)

    if not summary and not domain_class:
        return ""

    if summary and domain_class:
        return re.sub(r"\s+", " ", f"{summary} Functional class: {domain_class}.").strip()
    if summary:
        return summary
    return f"Functional class: {domain_class}."


# -----------------------------
# MyGene.info
# -----------------------------
def query_mygene(gene: str, organism_cfg: Dict[str, Any]) -> Dict[str, Any]:
    fields = ",".join([
        "symbol",
        "name",
        "alias",
        "taxid",
        "entrezgene",
        "ensembl.gene",
        "uniprot.Swiss-Prot",
        "genomic_pos",
        "summary",
        "type_of_gene",
        "go",
        "pathway",
        "map_location",
    ])

    params = {
        "q": gene,
        "species": organism_cfg["mygene"],
        "size": 10,
        "fields": fields,
        "scopes": "symbol,alias,ensemblgene,entrezgene",
    }
    data = safe_get_json(MYGENE_QUERY_URL, params=params)
    hits = data.get("hits", [])

    if not hits:
        return {}

    gene_upper = gene.upper()

    def hit_rank(hit):
        symbol = str(hit.get("symbol", "")).upper()
        aliases = [str(x).upper() for x in _as_list(hit.get("alias"))]
        exact_symbol = 1 if symbol == gene_upper else 0
        alias_match = 1 if gene_upper in aliases else 0
        score = float(hit.get("_score", 0))
        return (exact_symbol, alias_match, score)

    return sorted(hits, key=hit_rank, reverse=True)[0]


def parse_ensembl_gene(hit: Dict[str, Any]) -> str:
    ens = hit.get("ensembl")
    if isinstance(ens, dict):
        return str(ens.get("gene", "") or "")
    if isinstance(ens, list) and ens:
        first = ens[0]
        if isinstance(first, dict):
            return str(first.get("gene", "") or "")
    return ""


def parse_uniprot_accession_from_mygene(hit: Dict[str, Any]) -> str:
    up = hit.get("uniprot")
    if not up:
        return ""
    swiss = up.get("Swiss-Prot") if isinstance(up, dict) else None
    if isinstance(swiss, list):
        return str(swiss[0]) if swiss else ""
    return str(swiss or "")


def parse_genomic_location(hit: Dict[str, Any]) -> Tuple[str, str, str]:
    gp = hit.get("genomic_pos")
    if isinstance(gp, list) and gp:
        gp = gp[0]
    if not isinstance(gp, dict):
        return "", "", ""

    chr_ = str(gp.get("chr", "") or "")
    strand = gp.get("strand", "")
    start = gp.get("start", "")
    end = gp.get("end", "")
    loc = ""
    if chr_ and start != "" and end != "":
        loc = f"chr{chr_}:{start}-{end}"
    strand_str = str(strand) if strand != "" else ""
    return chr_, strand_str, loc


def parse_go_terms(go_obj: Dict[str, Any], max_items: int = 5) -> Tuple[str, str, str]:
    if not isinstance(go_obj, dict):
        return "", "", ""

    def extract_branch(branch):
        vals = _as_list(go_obj.get(branch))
        terms = []
        for x in vals:
            if isinstance(x, dict):
                term = x.get("term")
                if term and term not in terms:
                    terms.append(term)
        return _join_unique(terms[:max_items])

    bp = extract_branch("BP")
    mf = extract_branch("MF")
    cc = extract_branch("CC")
    return mf, bp, cc


def parse_pathways_split(pathway_obj: Dict[str, Any], max_items: int = 5) -> Dict[str, str]:
    out = {
        "reactome_pathways": "",
        "kegg_pathways": "",
        "hallmark_pathways": "",
        "pathway_annotation_summary": "",
    }

    if not isinstance(pathway_obj, dict):
        return out

    collected = {
        "reactome_pathways": [],
        "kegg_pathways": [],
    }

    source_map = {
        "reactome": "reactome_pathways",
        "kegg": "kegg_pathways",
    }

    for src, target_col in source_map.items():
        vals = _as_list(pathway_obj.get(src))
        for x in vals:
            if isinstance(x, dict):
                name = str(x.get("name", "") or "").strip()
                if name and name not in collected[target_col]:
                    collected[target_col].append(name)

    out["reactome_pathways"] = _join_unique(collected["reactome_pathways"][:max_items])
    out["kegg_pathways"] = _join_unique(collected["kegg_pathways"][:max_items])

    summary_parts = []
    if out["reactome_pathways"]:
        summary_parts.append(f"Reactome: {out['reactome_pathways']}")
    if out["kegg_pathways"]:
        summary_parts.append(f"KEGG: {out['kegg_pathways']}")
    if out["hallmark_pathways"]:
        summary_parts.append(f"Hallmark: {out['hallmark_pathways']}")

    out["pathway_annotation_summary"] = " | ".join(summary_parts)
    return out


# -----------------------------
# UniProt
# -----------------------------
def search_uniprot_accession(symbol: str, organism_cfg: Dict[str, Any]) -> str:
    query = f'(gene_exact:{symbol}) AND (organism_id:{organism_cfg["uniprot_organism_id"]}) AND (reviewed:true)'
    params = {
        "query": query,
        "format": "json",
        "size": 1,
        "fields": "accession,gene_names,protein_name,organism_name",
    }
    data = safe_get_json(UNIPROT_SEARCH_URL, params=params)
    results = data.get("results", [])
    if not results:
        return ""
    return str(results[0].get("primaryAccession", "") or "")


def fetch_uniprot_entry(accession: str) -> Dict[str, Any]:
    if not accession:
        return {}
    return safe_get_json(UNIPROT_ENTRY_URL.format(accession=accession))


def parse_uniprot_protein_name(entry: Dict[str, Any]) -> str:
    pdsc = entry.get("proteinDescription", {})
    if not isinstance(pdsc, dict):
        return ""
    rec = pdsc.get("recommendedName", {})
    if not isinstance(rec, dict):
        return ""
    full = rec.get("fullName", {})
    if isinstance(full, dict):
        return str(full.get("value", "") or "")
    return ""


def parse_uniprot_length(entry: Dict[str, Any]) -> str:
    seq = entry.get("sequence", {})
    if isinstance(seq, dict):
        length = seq.get("length")
        return str(length) if length is not None else ""
    return ""


def parse_uniprot_function(entry: Dict[str, Any]) -> str:
    comments = _as_list(entry.get("comments"))
    texts = []
    for c in comments:
        if not isinstance(c, dict):
            continue
        if c.get("commentType") == "FUNCTION":
            for t in _as_list(c.get("texts")):
                if isinstance(t, dict) and t.get("value"):
                    texts.append(t["value"])
    return summarize_text(_join_unique(texts), max_len=600)


def parse_uniprot_subcellular_location(entry: Dict[str, Any], max_items: int = 5) -> str:
    comments = _as_list(entry.get("comments"))
    locs = []
    for c in comments:
        if not isinstance(c, dict):
            continue
        if c.get("commentType") == "SUBCELLULAR LOCATION":
            for sl in _as_list(c.get("subcellularLocations")):
                if not isinstance(sl, dict):
                    continue
                location = sl.get("location", {})
                if isinstance(location, dict) and location.get("value"):
                    loc = str(location["value"]).strip()
                    if loc and loc not in locs:
                        locs.append(loc)
    return _join_unique(locs[:max_items])


def parse_uniprot_disease(entry: Dict[str, Any], max_items: int = 5) -> str:
    comments = _as_list(entry.get("comments"))
    diseases = []

    for c in comments:
        if not isinstance(c, dict):
            continue
        if c.get("commentType") != "DISEASE":
            continue

        disease = c.get("disease", {})
        if not isinstance(disease, dict):
            continue

        label = ""
        if disease.get("description"):
            desc = str(disease["description"]).strip()
            desc = re.sub(r"\s+", " ", desc)
            label = desc.split(".")[0]
        elif disease.get("acronym"):
            label = str(disease["acronym"]).strip()
        elif disease.get("diseaseId"):
            label = str(disease["diseaseId"]).strip()

        if label and label not in diseases:
            diseases.append(label)

    return _join_unique(diseases[:max_items])


def parse_uniprot_keywords(entry: Dict[str, Any], max_items: int = 8) -> str:
    kws = []
    for k in _as_list(entry.get("keywords")):
        if isinstance(k, dict) and k.get("name"):
            name = str(k["name"]).strip()
            if name and name not in kws:
                kws.append(name)
    return _join_unique(kws[:max_items])


def parse_uniprot_features(entry: Dict[str, Any]) -> Tuple[str, str, str, str]:
    features = _as_list(entry.get("features"))
    domain_like_types = {
        "Domain",
        "Region",
        "Repeat",
        "Zinc finger",
        "DNA binding",
        "Motif",
        "Coiled coil",
        "Transmembrane",
        "Intramembrane",
        "Topological domain",
        "Compositional bias",
        "Active site",
        "Binding site",
    }

    annotations = []
    coords = []
    classes = []

    for feat in features:
        if not isinstance(feat, dict):
            continue
        ftype = str(feat.get("type", "") or "")
        if ftype not in domain_like_types:
            continue

        desc = str(feat.get("description", "") or ftype).strip()
        loc = feat.get("location", {})
        start = ""
        end = ""
        if isinstance(loc, dict):
            b = loc.get("start", {})
            e = loc.get("end", {})
            if isinstance(b, dict):
                start = str(b.get("value", "") or "")
            if isinstance(e, dict):
                end = str(e.get("value", "") or "")

        if desc and desc not in annotations:
            annotations.append(desc)
        if desc and (start or end):
            coords.append(f"{desc}: {start}-{end}")

        dlow = desc.lower()
        if "bromo" in dlow:
            classes.append("reader")
        if "chromo" in dlow or "pwwp" in dlow or "tudor" in dlow or "mbt" in dlow:
            classes.append("reader")
        if "set" in dlow or "methyltransferase" in dlow or "acetyltransferase" in dlow:
            classes.append("writer")
        if "deacetylase" in dlow or "demethylase" in dlow:
            classes.append("eraser")
        if "kinase" in dlow:
            classes.append("kinase")
        if "homeobox" in dlow or "zinc finger" in dlow or "dna binding" in dlow:
            classes.append("transcription factor")
        if "rrm" in dlow or "rna recognition" in dlow or "kh domain" in dlow:
            classes.append("rna-binding protein")
        if "transmembrane" in dlow:
            classes.append("membrane protein")

    functional_domain_class = _join_unique(classes)

    summary_parts = []
    if annotations:
        summary_parts.append(f"Contains {', '.join(annotations[:5])}")
    if functional_domain_class:
        summary_parts.append(f"Class: {functional_domain_class}")
    key_domain_summary = summarize_text(". ".join(summary_parts), max_len=300)

    return (
        _join_unique(annotations),
        _join_unique(coords),
        key_domain_summary,
        functional_domain_class,
    )


def infer_protein_class_and_druggability(keyword_text: str, functional_domain_class: str) -> Tuple[str, str]:
    text = " ".join([keyword_text or "", functional_domain_class or ""]).lower()

    protein_class = "other"
    druggability = ""

    if "kinase" in text:
        protein_class = "kinase"
        druggability = "potentially_druggable"
    elif "receptor" in text:
        protein_class = "receptor"
        druggability = "potentially_druggable"
    elif "transcription factor" in text:
        protein_class = "transcription factor"
        druggability = "challenging_but_relevant"
    elif "reader" in text:
        protein_class = "reader"
        druggability = "potentially_druggable"
    elif "writer" in text:
        protein_class = "writer"
        druggability = "potentially_druggable"
    elif "eraser" in text:
        protein_class = "eraser"
        druggability = "potentially_druggable"
    elif "rna-binding" in text or "rna binding" in text:
        protein_class = "rna-binding protein"
    elif "membrane protein" in text:
        protein_class = "membrane protein"
    elif "enzyme" in text or "catalytic activity" in text:
        protein_class = "enzyme"
        druggability = "potentially_druggable"
    elif "chromatin" in text:
        protein_class = "chromatin regulator"

    return protein_class, druggability


# -----------------------------
# Main annotation logic
# -----------------------------
def annotate_one_gene(input_gene: str, organism_cfg: Dict[str, Any]) -> Dict[str, Any]:
    row = {k: "" for k in OUTPUT_SCHEMA}
    row["input_gene"] = input_gene
    row["organism"] = organism_cfg["label"]

    hit = query_mygene(input_gene, organism_cfg)

    if not hit:
        row["mapping_status"] = "not_found"
        row["mapping_note"] = "No MyGene.info hit found"
        row["source_summary"] = "MyGene"
        return row

    symbol = str(hit.get("symbol", "") or "")
    aliases = _join_unique(_as_list(hit.get("alias")))
    input_upper = input_gene.upper()

    if input_upper == symbol.upper():
        mapping_status = "exact_symbol"
        mapping_note = "Matched input to approved symbol"
    elif input_upper in [str(x).upper() for x in _as_list(hit.get("alias"))]:
        mapping_status = "alias_match"
        mapping_note = "Mapped input alias to approved symbol"
    else:
        mapping_status = "best_match"
        mapping_note = "Best available MyGene match"

    row["mapping_status"] = mapping_status
    row["mapping_note"] = mapping_note
    row["approved_symbol"] = symbol
    row["full_gene_name"] = str(hit.get("name", "") or "")
    row["aliases"] = aliases
    row["entrez_id"] = str(hit.get("entrezgene", "") or "")
    row["ensembl_gene_id"] = parse_ensembl_gene(hit)
    row["gene_type"] = str(hit.get("type_of_gene", "") or "")
    row["cytoband"] = str(hit.get("map_location", "") or "")

    chr_, strand, genomic_location = parse_genomic_location(hit)
    row["chromosome"] = chr_
    row["strand"] = strand
    row["genomic_location"] = genomic_location

    mygene_summary = str(hit.get("summary", "") or "")
    mf, bp, cc = parse_go_terms(hit.get("go", {}), max_items=5)
    row["molecular_function"] = mf
    row["biological_process"] = bp
    row["cellular_component"] = cc

    pathway_info = parse_pathways_split(hit.get("pathway", {}), max_items=5)
    row["reactome_pathways"] = pathway_info["reactome_pathways"]
    row["kegg_pathways"] = pathway_info["kegg_pathways"]
    row["hallmark_pathways"] = pathway_info["hallmark_pathways"]
    row["pathway_annotation_summary"] = pathway_info["pathway_annotation_summary"]

    up_acc = parse_uniprot_accession_from_mygene(hit)
    if not up_acc:
        up_acc = search_uniprot_accession(symbol, organism_cfg)
    row["uniprot_accession"] = up_acc

    up_func = ""

    if up_acc:
        entry = fetch_uniprot_entry(up_acc)
        row["protein_name"] = parse_uniprot_protein_name(entry)
        row["protein_length"] = parse_uniprot_length(entry)
        up_func = parse_uniprot_function(entry)

        row["subcellular_location"] = parse_uniprot_subcellular_location(entry, max_items=5)
        row["disease_association"] = parse_uniprot_disease(entry, max_items=5)

        domain_annotation, domain_coordinates, key_domain_summary, functional_domain_class = parse_uniprot_features(entry)
        row["domain_annotation"] = domain_annotation
        row["domain_coordinates"] = domain_coordinates
        row["key_domain_summary"] = key_domain_summary
        row["functional_domain_class"] = functional_domain_class

        keyword_text = parse_uniprot_keywords(entry, max_items=8)
        row["keyword_annotation"] = keyword_text

        inferred_class, inferred_druggability = infer_protein_class_and_druggability(
            keyword_text,
            row["functional_domain_class"],
        )
        row["protein_class"] = inferred_class
        row["druggability"] = inferred_druggability

    row["reported_function"] = concise_function_summary(mygene_summary, up_func, max_len=220)
    row["reported_function"] = augment_function_with_domain(
        row["reported_function"],
        row["functional_domain_class"],
    )

    sources = ["MyGene"]
    if row["uniprot_accession"]:
        sources.append("UniProt")
    row["source_summary"] = _join_unique(sources)

    return row


def main():
    parser = argparse.ArgumentParser(
        description="Comprehensive annotation for one gene or a gene list"
    )
    parser.add_argument("--input", required=True, help="Input gene list file")
    parser.add_argument("--output", required=True, help="Output TSV or CSV")
    parser.add_argument("--organism", default="human", help="Supported: human, mouse")
    parser.add_argument("--output-format", choices=["tsv", "csv"], default="tsv")
    args = parser.parse_args()

    organism_cfg = normalize_organism(args.organism)
    genes = read_gene_list(args.input)

    rows = []
    for i, gene in enumerate(genes, start=1):
        print(f"[INFO] Annotating {gene} ({i}/{len(genes)})")
        try:
            rows.append(annotate_one_gene(gene, organism_cfg))
        except Exception as e:
            fallback = {k: "" for k in OUTPUT_SCHEMA}
            fallback["input_gene"] = gene
            fallback["mapping_status"] = "error"
            fallback["mapping_note"] = str(e)
            rows.append(fallback)
        time.sleep(0.1)

    out = pd.DataFrame(rows)
    for col in OUTPUT_SCHEMA:
        if col not in out.columns:
            out[col] = ""
    out = out[OUTPUT_SCHEMA].copy()

    ensure_parent_dir(args.output)
    if args.output_format == "csv":
        out.to_csv(args.output, index=False)
    else:
        out.to_csv(args.output, sep="\t", index=False)

    print(f"[INFO] Genes processed: {len(out)}")
    print(f"[INFO] Saved annotation table to: {args.output}")


if __name__ == "__main__":
    main()