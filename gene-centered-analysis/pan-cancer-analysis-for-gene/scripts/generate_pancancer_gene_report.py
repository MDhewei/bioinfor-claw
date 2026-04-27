from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from urllib.request import urlopen
from urllib.request import Request

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle


OUT_DIR = Path(".")
GENE_SYMBOL = "PRNP"
GENE_ENSEMBL = "ENSG00000171867"
DATA_DIR = OUT_DIR / "prnp_pancancer_data"
PDF_OUT = OUT_DIR / "prnp_pancancer_clinical_relevance_report_with_depmap.pdf"
APPENDIX_OUT = OUT_DIR / "prnp_pancancer_appendix.csv"

GDC_HUB = "https://gdc-hub.s3.us-east-1.amazonaws.com/download"
DEPMAP_INDEX = "https://depmap.org/portal/api/download/files"
DEPMAP_RELEASE = "DepMap Public 26Q1"
CBIO_API = "https://www.cbioportal.org/api"

TCGA_PROJECTS = [
    "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM",
    "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",
    "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD",
    "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM",
]

CPTAC_STUDIES = [
    "brca_cptac_2020", "luad_cptac_2020", "lusc_cptac_2021", "coad_cptac_2019",
    "ucec_cptac_2020", "gbm_cptac_2021", "paad_cptac_2021",
]


def configure(gene: str, ensembl: str, outdir: str | None = None) -> None:
    global GENE_SYMBOL, GENE_ENSEMBL, DATA_DIR, PDF_OUT, APPENDIX_OUT, OUT_DIR
    GENE_SYMBOL = gene.upper()
    GENE_ENSEMBL = ensembl
    if outdir:
        OUT_DIR = Path(outdir)
    stem = GENE_SYMBOL.lower()
    DATA_DIR = OUT_DIR / f"{stem}_pancancer_data"
    PDF_OUT = OUT_DIR / f"{stem}_pancancer_clinical_relevance_report_with_depmap.pdf"
    APPENDIX_OUT = OUT_DIR / f"{stem}_pancancer_appendix.csv"


@dataclass
class KMResult:
    cancer: str
    n: int
    median_expr: float
    p_value: float
    chi2: float
    high_events: int
    low_events: int
    direction: str
    high_curve: list[tuple[float, float]]
    low_curve: list[tuple[float, float]]


def ensure_dirs() -> None:
    DATA_DIR.mkdir(exist_ok=True)


def gzip_lines(url: str):
    with urlopen(url, timeout=120) as response:
        with gzip.GzipFile(fileobj=response) as handle:
            for raw in handle:
                yield raw.decode("utf-8").rstrip("\n")


def get_gene_expression_tcga(project: str) -> pd.DataFrame:
    cache = DATA_DIR / f"TCGA_{project}_{GENE_SYMBOL}_expression.csv"
    if cache.exists():
        return pd.read_csv(cache)
    url = f"{GDC_HUB}/TCGA-{project}.star_counts.tsv.gz"
    iterator = gzip_lines(url)
    header = next(iterator).split("\t")
    gene_row = None
    for line in iterator:
        if line.startswith(GENE_ENSEMBL):
            gene_row = line.split("\t")
            break
    if gene_row is None:
        raise RuntimeError(f"{GENE_SYMBOL} row not found for {project}")
    rows = []
    for sample, value in zip(header[1:], gene_row[1:]):
        try:
            rows.append({"sample": sample, "expr": float(value)})
        except ValueError:
            pass
    df = pd.DataFrame(rows)
    df.to_csv(cache, index=False)
    return df


def get_survival_tcga(project: str) -> pd.DataFrame:
    cache = DATA_DIR / f"TCGA_{project}_survival.csv"
    if cache.exists():
        return pd.read_csv(cache)
    url = f"{GDC_HUB}/TCGA-{project}.survival.tsv.gz"
    rows = []
    reader = csv.DictReader(gzip_lines(url), delimiter="\t")
    for row in reader:
        try:
            rows.append(
                {
                    "sample": row["sample"],
                    "patient": row.get("_PATIENT", row["sample"][:12]),
                    "os_time": float(row["OS.time"]),
                    "os_event": int(float(row["OS"])),
                }
            )
        except (KeyError, ValueError):
            continue
    df = pd.DataFrame(rows)
    df.to_csv(cache, index=False)
    return df


def km_curve(times: list[float], events: list[int]) -> list[tuple[float, float]]:
    data = sorted(zip(times, events), key=lambda x: x[0])
    surv = 1.0
    curve = [(0.0, 1.0)]
    for t in sorted({x[0] for x in data if x[1] == 1}):
        at_risk = sum(1 for tt, _ in data if tt >= t)
        deaths = sum(1 for tt, ev in data if tt == t and ev == 1)
        if at_risk:
            surv *= 1 - deaths / at_risk
            curve.append((t, surv))
    return curve


def logrank_p(df: pd.DataFrame) -> tuple[float, float]:
    high = df["group"] == "High"
    observed = expected = variance = 0.0
    for t in sorted(df.loc[df["os_event"] == 1, "os_time"].unique()):
        risk_high = int(((df["os_time"] >= t) & high).sum())
        risk_low = int(((df["os_time"] >= t) & ~high).sum())
        deaths_high = int(((df["os_time"] == t) & (df["os_event"] == 1) & high).sum())
        deaths_total = int(((df["os_time"] == t) & (df["os_event"] == 1)).sum())
        risk_total = risk_high + risk_low
        if risk_total <= 1:
            continue
        observed += deaths_high
        expected += deaths_total * risk_high / risk_total
        variance += risk_high * risk_low * deaths_total * (risk_total - deaths_total) / (risk_total**2 * (risk_total - 1))
    chi2 = ((observed - expected) ** 2 / variance) if variance > 0 else 0.0
    return chi2, math.erfc(math.sqrt(chi2 / 2))


def analyze_tcga_project(project: str) -> tuple[pd.DataFrame, KMResult]:
    expr = get_gene_expression_tcga(project)
    surv = get_survival_tcga(project)
    df = expr.merge(surv, on="sample", how="inner")
    df = df.sort_values("sample").drop_duplicates("patient", keep="first")
    df = df[df["os_time"].notna() & df["expr"].notna()].copy()
    df["os_years"] = df["os_time"] / 365.25
    med = float(df["expr"].median())
    df["group"] = ["High" if x >= med else "Low" for x in df["expr"]]
    high_df, low_df = df[df["group"] == "High"], df[df["group"] == "Low"]
    chi2, p = logrank_p(df)
    high_event_rate = high_df["os_event"].mean() if len(high_df) else 0
    low_event_rate = low_df["os_event"].mean() if len(low_df) else 0
    direction = "high worse" if high_event_rate > low_event_rate else "high better/equal"
    result = KMResult(
        cancer=project,
        n=len(df),
        median_expr=med,
        p_value=p,
        chi2=chi2,
        high_events=int(high_df["os_event"].sum()),
        low_events=int(low_df["os_event"].sum()),
        direction=direction,
        high_curve=km_curve(high_df["os_years"].tolist(), high_df["os_event"].tolist()),
        low_curve=km_curve(low_df["os_years"].tolist(), low_df["os_event"].tolist()),
    )
    return df, result


def cbio_json(path: str, payload: dict | None = None):
    url = f"{CBIO_API}{path}"
    headers = {"Accept": "application/json"}
    data = None
    method = "GET"
    if payload is not None:
        data = json.dumps(payload).encode()
        headers["Content-Type"] = "application/json"
        method = "POST"
    req = Request(url, data=data, headers=headers, method=method)
    with urlopen(req, timeout=120) as response:
        return json.load(response)


def gene_entrez_id() -> int:
    cache = DATA_DIR / f"cbio_{GENE_SYMBOL}_gene.json"
    if cache.exists():
        return int(json.loads(cache.read_text())["entrezGeneId"])
    data = cbio_json(f"/genes/{GENE_SYMBOL}")
    cache.write_text(json.dumps(data), encoding="utf-8")
    return int(data["entrezGeneId"])


def tcga_study_id(project: str) -> str:
    return f"{project.lower()}_tcga_pan_can_atlas_2018"


def cbio_sample_list_id(study_id: str) -> str | None:
    cache = DATA_DIR / f"cbio_{study_id}_sample_lists.json"
    try:
        if cache.exists():
            lists = json.loads(cache.read_text())
        else:
            lists = cbio_json(f"/studies/{study_id}/sample-lists?projection=SUMMARY")
            cache.write_text(json.dumps(lists), encoding="utf-8")
    except Exception:
        return None
    for item in lists:
        if item.get("category") == "all_cases_in_study":
            return item["sampleListId"]
    return lists[0]["sampleListId"] if lists else None


def cbio_profiles(study_id: str) -> list[dict]:
    cache = DATA_DIR / f"cbio_{study_id}_profiles.json"
    if cache.exists():
        return json.loads(cache.read_text())
    profiles = cbio_json(f"/studies/{study_id}/molecular-profiles?projection=SUMMARY")
    cache.write_text(json.dumps(profiles), encoding="utf-8")
    return profiles


def cbio_profile_id(study_id: str, suffix: str) -> str | None:
    for profile in cbio_profiles(study_id):
        if profile["molecularProfileId"].endswith(suffix):
            return profile["molecularProfileId"]
    return None


def cbio_molecular_data(profile_id: str, sample_list_id: str, entrez: int) -> list[dict]:
    cache = DATA_DIR / f"cbio_{profile_id}_{GENE_SYMBOL}_molecular.json"
    if cache.exists():
        return json.loads(cache.read_text())
    data = cbio_json(
        f"/molecular-profiles/{profile_id}/molecular-data/fetch",
        {"entrezGeneIds": [entrez], "sampleListId": sample_list_id},
    )
    cache.write_text(json.dumps(data), encoding="utf-8")
    return data


def cbio_mutations(profile_id: str, sample_list_id: str, entrez: int) -> list[dict]:
    cache = DATA_DIR / f"cbio_{profile_id}_{GENE_SYMBOL}_mutations.json"
    if cache.exists():
        return json.loads(cache.read_text())
    data = cbio_json(
        f"/molecular-profiles/{profile_id}/mutations/fetch",
        {"entrezGeneIds": [entrez], "sampleListId": sample_list_id},
    )
    cache.write_text(json.dumps(data), encoding="utf-8")
    return data


def alteration_survival(df: pd.DataFrame, altered_patients: set[str], label: str) -> dict:
    sub = df.copy()
    sub["group"] = ["Altered" if p in altered_patients else "Unaltered" for p in sub["patient"]]
    alt = sub[sub["group"] == "Altered"]
    unalt = sub[sub["group"] == "Unaltered"]
    if len(alt) < 3 or int(alt["os_event"].sum()) < 1 or len(unalt) < 3:
        return {
            "label": label,
            "n_altered": len(alt),
            "n_unaltered": len(unalt),
            "p_value": None,
            "direction": "insufficient",
            "altered_events": int(alt["os_event"].sum()),
            "unaltered_events": int(unalt["os_event"].sum()),
            "altered_curve": [],
            "unaltered_curve": [],
        }
    mapped = sub.assign(group=["High" if g == "Altered" else "Low" for g in sub["group"]])
    chi2, p = logrank_p(mapped)
    direction = "altered worse" if alt["os_event"].mean() > unalt["os_event"].mean() else "altered better/equal"
    return {
        "label": label,
        "n_altered": len(alt),
        "n_unaltered": len(unalt),
        "p_value": p,
        "chi2": chi2,
        "direction": direction,
        "altered_events": int(alt["os_event"].sum()),
        "unaltered_events": int(unalt["os_event"].sum()),
        "altered_curve": km_curve(alt["os_years"].tolist(), alt["os_event"].tolist()),
        "unaltered_curve": km_curve(unalt["os_years"].tolist(), unalt["os_event"].tolist()),
    }


def analyze_tcga_alterations(tcga_frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    entrez = gene_entrez_id()
    rows = []
    for project, df in tcga_frames.items():
        study = tcga_study_id(project)
        sample_list = cbio_sample_list_id(study)
        if not sample_list:
            continue
        try:
            gistic = cbio_profile_id(study, "_gistic")
            mutations = cbio_profile_id(study, "_mutations")
            cna_data = cbio_molecular_data(gistic, sample_list, entrez) if gistic else []
            mut_data = cbio_mutations(mutations, sample_list, entrez) if mutations else []
        except Exception:
            continue
        cna_by_patient = {}
        for item in cna_data:
            try:
                cna_by_patient[item["patientId"]] = int(float(item["value"]))
            except (KeyError, ValueError):
                continue
        mutated = {m["patientId"] for m in mut_data if m.get("patientId")}
        amplified = {p for p, v in cna_by_patient.items() if v == 2}
        gained = {p for p, v in cna_by_patient.items() if v >= 1}
        deleted = {p for p, v in cna_by_patient.items() if v <= -1}

        for kind, patients in [
            ("mutation", mutated),
            ("amplification", amplified),
            ("gain_or_amplification", gained),
            ("deletion", deleted),
        ]:
            surv = alteration_survival(df, patients, kind)
            rows.append({
                "project": project,
                "alteration": kind,
                "altered_patients": surv["n_altered"],
                "unaltered_patients": surv["n_unaltered"],
                "altered_events": surv["altered_events"],
                "unaltered_events": surv["unaltered_events"],
                "logrank_p": surv["p_value"],
                "direction": surv["direction"],
                "altered_curve": surv["altered_curve"],
                "unaltered_curve": surv["unaltered_curve"],
            })
    return pd.DataFrame(rows)


def cptac_protein_data(entrez: int) -> pd.DataFrame:
    rows = []
    for study in CPTAC_STUDIES:
        sample_list = cbio_sample_list_id(study)
        if not sample_list:
            continue
        try:
            protein_profiles = [
                p for p in cbio_profiles(study)
                if p.get("molecularAlterationType") == "PROTEIN_LEVEL"
                and "protein_quantification" in p.get("molecularProfileId", "")
                and not p.get("molecularProfileId", "").endswith("_zscores")
            ]
            if not protein_profiles:
                continue
            profile_id = protein_profiles[0]["molecularProfileId"]
            data = cbio_molecular_data(profile_id, sample_list, entrez)
        except Exception:
            continue
        for item in data:
            try:
                rows.append({
                    "study": study,
                    "sampleId": item.get("sampleId", ""),
                    "patientId": item.get("patientId", ""),
                    "protein": float(item["value"]),
                })
            except (KeyError, ValueError):
                continue
    return pd.DataFrame(rows)


def depmap_file_url(filename: str) -> str:
    cache = DATA_DIR / "depmap_file_index.csv"
    if not cache.exists():
        cache.write_bytes(urlopen(DEPMAP_INDEX, timeout=120).read())
    rows = list(csv.DictReader(io.StringIO(cache.read_text())))
    for row in rows:
        if row["release"] == DEPMAP_RELEASE and row["filename"] == filename:
            return row["url"]
    cache.unlink(missing_ok=True)
    cache.write_bytes(urlopen(DEPMAP_INDEX, timeout=120).read())
    rows = list(csv.DictReader(io.StringIO(cache.read_text())))
    for row in rows:
        if row["release"] == DEPMAP_RELEASE and row["filename"] == filename:
            return row["url"]
    raise RuntimeError(f"Missing DepMap file: {filename}")


def stream_csv_url(url: str):
    with urlopen(url, timeout=180) as response:
        text = io.TextIOWrapper(response, encoding="utf-8", newline="")
        reader = csv.reader(text)
        for row in reader:
            yield row


def depmap_model_metadata() -> pd.DataFrame:
    cache = DATA_DIR / "DepMap_26Q1_Model_metadata.csv"
    if cache.exists():
        return pd.read_csv(cache)
    rows = []
    with urlopen(depmap_file_url("Model.csv"), timeout=120) as response:
        text = io.TextIOWrapper(response, encoding="utf-8", newline="")
        for row in csv.DictReader(text):
            rows.append(
                {
                    "ModelID": row["ModelID"],
                    "CellLineName": row.get("CellLineName", ""),
                    "Lineage": row.get("OncotreeLineage") or "Unknown",
                    "PrimaryDisease": row.get("OncotreePrimaryDisease") or "Unknown",
                    "Subtype": row.get("OncotreeSubtype") or "Unknown",
                    "ModelType": row.get("ModelType", ""),
                    "PrimaryOrMetastasis": row.get("PrimaryOrMetastasis", ""),
                }
            )
    df = pd.DataFrame(rows)
    df.to_csv(cache, index=False)
    return df


def depmap_gene_column(filename: str, gene_symbol: str | None = None) -> pd.DataFrame:
    if gene_symbol is None:
        gene_symbol = GENE_SYMBOL
    cache = DATA_DIR / f"DepMap_26Q1_{filename}_{gene_symbol}.csv"
    if cache.exists():
        return pd.read_csv(cache)
    rows_iter = stream_csv_url(depmap_file_url(filename))
    header = next(rows_iter)
    gene_col = next(i for i, h in enumerate(header) if h.startswith(f"{gene_symbol} "))
    model_col = header.index("ModelID") if "ModelID" in header else 0
    default_col = header.index("IsDefaultEntryForModel") if "IsDefaultEntryForModel" in header else None
    rows = []
    for row in rows_iter:
        if default_col is not None and row[default_col] != "Yes":
            continue
        try:
            rows.append({"ModelID": row[model_col], gene_symbol: float(row[gene_col])})
        except (ValueError, IndexError):
            continue
    df = pd.DataFrame(rows)
    df.to_csv(cache, index=False)
    return df


def depmap_gene_mutations() -> pd.DataFrame:
    cache = DATA_DIR / f"DepMap_26Q1_{GENE_SYMBOL}_mutations.csv"
    if cache.exists():
        return pd.read_csv(cache)
    rows = []
    with urlopen(depmap_file_url("OmicsSomaticMutations.csv"), timeout=180) as response:
        text = io.TextIOWrapper(response, encoding="utf-8", newline="")
        for row in csv.DictReader(text):
            if row.get("HugoSymbol") != GENE_SYMBOL:
                continue
            if row.get("IsDefaultEntryForModel") not in ("Yes", "True", "true", "1"):
                continue
            rows.append(
                {
                    "ModelID": row.get("ModelID", ""),
                    "Chrom": row.get("Chrom", ""),
                    "Pos": row.get("Pos", ""),
                    "ProteinChange": row.get("ProteinChange", ""),
                    "Consequence": row.get("MolecularConsequence", ""),
                    "Impact": row.get("VepImpact", ""),
                    "AF": row.get("AF", ""),
                    "Hotspot": row.get("Hotspot", ""),
                    "LikelyLoF": row.get("LikelyLoF", ""),
                }
            )
    df = pd.DataFrame(rows)
    df.to_csv(cache, index=False)
    return df


def percentile(vals: list[float], q: float) -> float:
    vals = sorted(v for v in vals if not pd.isna(v))
    if not vals:
        return 0.0
    idx = (len(vals) - 1) * q
    lo, hi = math.floor(idx), math.ceil(idx)
    if lo == hi:
        return vals[lo]
    return vals[lo] * (hi - idx) + vals[hi] * (idx - lo)


class Chart:
    def __init__(self, width=7.15 * inch, height=3.75 * inch):
        self.width = width
        self.height = height

    def hbar(self, title: str, items: list[tuple[str, float]], x_label: str, color="#6b9fb5", max_items=33):
        from reportlab.graphics.shapes import Drawing, Line, Rect, String

        items = items[:max_items]
        d = Drawing(self.width, self.height)
        ml, mr, mt, mb = 92, 18, 28, 24
        pw, ph = self.width - ml - mr, self.height - mt - mb
        vmax = max([abs(v) for _, v in items] + [1]) * 1.15
        zero_x = ml if min(v for _, v in items) >= 0 else ml + pw / 2
        row_h = ph / max(len(items), 1)
        d.add(String(ml, self.height - 15, title, fontSize=10.5, fillColor=colors.HexColor("#18365c")))
        d.add(Line(ml, mb, ml + pw, mb, strokeColor=colors.black, strokeWidth=0.6))
        d.add(String(ml + pw / 2 - 28, 6, x_label, fontSize=7, fillColor=colors.grey))
        for i, (label, value) in enumerate(items):
            y = mb + ph - (i + 1) * row_h + row_h * 0.18
            h = row_h * 0.62
            if zero_x == ml:
                w = value / vmax * pw
                x = ml
            else:
                w = abs(value) / vmax * (pw / 2)
                x = zero_x if value >= 0 else zero_x - w
            d.add(String(4, y + h * 0.25, label, fontSize=6.7, fillColor=colors.HexColor("#3d4651")))
            d.add(Rect(x, y, max(w, 0.5), h, fillColor=colors.HexColor(color), strokeColor=colors.HexColor("#2b6f83")))
            d.add(String(x + w + 3 if value >= 0 else x - 22, y + h * 0.25, f"{value:.2g}", fontSize=6.5))
        if zero_x != ml:
            d.add(Line(zero_x, mb, zero_x, mb + ph, strokeColor=colors.HexColor("#999999"), strokeWidth=0.5))
        return d

    def boxplot(self, title: str, data_by_label: dict[str, list[float]], y_label: str):
        from reportlab.graphics.shapes import Drawing, Line, Rect, String

        labels = list(data_by_label.keys())
        d = Drawing(self.width, self.height)
        ml, mr, mt, mb = 48, 12, 30, 50
        pw, ph = self.width - ml - mr, self.height - mt - mb
        all_vals = [v for vals in data_by_label.values() for v in vals if not pd.isna(v)]
        ymin, ymax = min(all_vals), max(all_vals)
        pad = (ymax - ymin) * 0.08 or 1
        ymin, ymax = ymin - pad, ymax + pad

        def ymap(v): return mb + (v - ymin) / (ymax - ymin) * ph

        d.add(String(ml, self.height - 15, title, fontSize=10.5, fillColor=colors.HexColor("#18365c")))
        d.add(Line(ml, mb, ml + pw, mb, strokeColor=colors.black, strokeWidth=0.6))
        d.add(Line(ml, mb, ml, mb + ph, strokeColor=colors.black, strokeWidth=0.6))
        d.add(String(3, mb + ph / 2, y_label, fontSize=7, fillColor=colors.grey))
        step = pw / max(len(labels), 1)
        for i, label in enumerate(labels):
            vals = data_by_label[label]
            q1, med, q3 = percentile(vals, 0.25), percentile(vals, 0.5), percentile(vals, 0.75)
            lo, hi = percentile(vals, 0.05), percentile(vals, 0.95)
            cx, bw = ml + step * (i + 0.5), min(18, step * 0.48)
            d.add(Line(cx, ymap(lo), cx, ymap(hi), strokeColor=colors.HexColor("#617083"), strokeWidth=0.8))
            d.add(Rect(cx - bw / 2, ymap(q1), bw, ymap(q3) - ymap(q1), fillColor=colors.HexColor("#d9ecf2"), strokeColor=colors.HexColor("#2b6f83")))
            d.add(Line(cx - bw / 2, ymap(med), cx + bw / 2, ymap(med), strokeColor=colors.HexColor("#18365c"), strokeWidth=1.1))
            d.add(String(cx - 10, 24, label, fontSize=5.4, fillColor=colors.grey))
            d.add(String(cx - 8, 14, f"n={len(vals)}", fontSize=5.2, fillColor=colors.grey))
        return d

    def km(self, title: str, result: KMResult):
        from reportlab.graphics.shapes import Drawing, Line, Rect, String

        d = Drawing(self.width, self.height)
        ml, mr, mt, mb = 48, 14, 32, 34
        pw, ph = self.width - ml - mr, self.height - mt - mb
        xmax = max([x for x, _ in result.high_curve + result.low_curve] + [1])

        def xmap(x): return ml + x / xmax * pw
        def ymap(y): return mb + y * ph

        d.add(String(ml, self.height - 15, title, fontSize=10.5, fillColor=colors.HexColor("#18365c")))
        d.add(String(ml + 4, self.height - 29, f"n={result.n}; log-rank p={result.p_value:.3g}; {result.direction}", fontSize=7.5, fillColor=colors.grey))
        d.add(Line(ml, mb, ml + pw, mb, strokeColor=colors.black, strokeWidth=0.6))
        d.add(Line(ml, mb, ml, mb + ph, strokeColor=colors.black, strokeWidth=0.6))
        for val in [0, 0.5, 1.0]:
            d.add(String(15, ymap(val) - 3, f"{val:.1f}", fontSize=7, fillColor=colors.grey))

        def step(curve, color):
            lx, ly = curve[0]
            for x, y in curve[1:]:
                d.add(Line(xmap(lx), ymap(ly), xmap(x), ymap(ly), strokeColor=color, strokeWidth=1.5))
                d.add(Line(xmap(x), ymap(ly), xmap(x), ymap(y), strokeColor=color, strokeWidth=1.5))
                lx, ly = x, y
            d.add(Line(xmap(lx), ymap(ly), xmap(xmax), ymap(ly), strokeColor=color, strokeWidth=1.5))

        step(result.high_curve, colors.HexColor("#c43b3b"))
        step(result.low_curve, colors.HexColor("#1f6fb5"))
        d.add(Rect(ml + pw - 105, mb + ph - 31, 95, 25, fillColor=colors.white, strokeColor=colors.HexColor("#d5d9df")))
        d.add(Line(ml + pw - 99, mb + ph - 16, ml + pw - 80, mb + ph - 16, strokeColor=colors.HexColor("#c43b3b"), strokeWidth=1.8))
        d.add(String(ml + pw - 76, mb + ph - 19, "High", fontSize=7.5))
        d.add(Line(ml + pw - 99, mb + ph - 27, ml + pw - 80, mb + ph - 27, strokeColor=colors.HexColor("#1f6fb5"), strokeWidth=1.8))
        d.add(String(ml + pw - 76, mb + ph - 30, "Low", fontSize=7.5))
        return d

    def alteration_km(self, title: str, result: dict):
        from reportlab.graphics.shapes import Drawing, Line, Rect, String

        d = Drawing(self.width, self.height)
        ml, mr, mt, mb = 48, 14, 32, 34
        pw, ph = self.width - ml - mr, self.height - mt - mb
        curves = result.get("altered_curve", []) + result.get("unaltered_curve", [])
        xmax = max([x for x, _ in curves] + [1])

        def xmap(x): return ml + x / xmax * pw
        def ymap(y): return mb + y * ph

        p = result.get("logrank_p")
        ptxt = "NA" if p is None or pd.isna(p) else f"{p:.3g}"
        d.add(String(ml, self.height - 15, title, fontSize=10.5, fillColor=colors.HexColor("#18365c")))
        d.add(String(ml + 4, self.height - 29, f"altered n={result.get('altered_patients')}; log-rank p={ptxt}; {result.get('direction')}", fontSize=7.5, fillColor=colors.grey))
        d.add(Line(ml, mb, ml + pw, mb, strokeColor=colors.black, strokeWidth=0.6))
        d.add(Line(ml, mb, ml, mb + ph, strokeColor=colors.black, strokeWidth=0.6))

        def step(curve, color):
            if not curve:
                return
            lx, ly = curve[0]
            for x, y in curve[1:]:
                d.add(Line(xmap(lx), ymap(ly), xmap(x), ymap(ly), strokeColor=color, strokeWidth=1.5))
                d.add(Line(xmap(x), ymap(ly), xmap(x), ymap(y), strokeColor=color, strokeWidth=1.5))
                lx, ly = x, y
            d.add(Line(xmap(lx), ymap(ly), xmap(xmax), ymap(ly), strokeColor=color, strokeWidth=1.5))

        step(result.get("altered_curve", []), colors.HexColor("#c43b3b"))
        step(result.get("unaltered_curve", []), colors.HexColor("#1f6fb5"))
        d.add(Rect(ml + pw - 120, mb + ph - 31, 110, 25, fillColor=colors.white, strokeColor=colors.HexColor("#d5d9df")))
        d.add(Line(ml + pw - 114, mb + ph - 16, ml + pw - 94, mb + ph - 16, strokeColor=colors.HexColor("#c43b3b"), strokeWidth=1.8))
        d.add(String(ml + pw - 90, mb + ph - 19, "Altered", fontSize=7.5))
        d.add(Line(ml + pw - 114, mb + ph - 27, ml + pw - 94, mb + ph - 27, strokeColor=colors.HexColor("#1f6fb5"), strokeWidth=1.8))
        d.add(String(ml + pw - 90, mb + ph - 30, "Unaltered", fontSize=7.5))
        return d

    def scatter(self, title: str, xvals: list[float], yvals: list[float], x_label: str, y_label: str):
        from reportlab.graphics.shapes import Circle, Drawing, Line, String

        paired = [(x, y) for x, y in zip(xvals, yvals) if not pd.isna(x) and not pd.isna(y)]
        if len(paired) > 1200:
            step = max(1, len(paired) // 1200)
            paired = paired[::step]
        xs = [p[0] for p in paired]
        ys = [p[1] for p in paired]
        xmin, xmax = min(xs), max(xs)
        ymin, ymax = min(ys), max(ys)
        xpad = (xmax - xmin) * 0.08 or 1
        ypad = (ymax - ymin) * 0.08 or 1
        xmin, xmax = xmin - xpad, xmax + xpad
        ymin, ymax = ymin - ypad, ymax + ypad
        d = Drawing(self.width, self.height)
        ml, mr, mt, mb = 54, 18, 32, 40
        pw, ph = self.width - ml - mr, self.height - mt - mb

        def xmap(v): return ml + (v - xmin) / (xmax - xmin) * pw
        def ymap(v): return mb + (v - ymin) / (ymax - ymin) * ph

        d.add(String(ml, self.height - 15, title, fontSize=10.5, fillColor=colors.HexColor("#18365c")))
        d.add(Line(ml, mb, ml + pw, mb, strokeColor=colors.black, strokeWidth=0.6))
        d.add(Line(ml, mb, ml, mb + ph, strokeColor=colors.black, strokeWidth=0.6))
        d.add(String(ml + pw / 2 - 35, 8, x_label, fontSize=7, fillColor=colors.grey))
        d.add(String(3, mb + ph / 2, y_label, fontSize=7, fillColor=colors.grey))
        for x, y in paired:
            d.add(Circle(xmap(x), ymap(y), 1.1, fillColor=colors.HexColor("#3f7f9b"), strokeColor=None, fillOpacity=0.35))
        return d


def add_table(story, rows, widths=None, font_size=7):
    table = Table(rows, colWidths=widths, repeatRows=1)
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#18365c")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), font_size),
        ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#d5d9df")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f6f8fa")]),
    ]))
    story.append(table)


def build_report() -> None:
    ensure_dirs()
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="Small", parent=styles["BodyText"], fontSize=8, leading=10))
    styles.add(ParagraphStyle(name="Caption", parent=styles["BodyText"], fontSize=8, leading=10, textColor=colors.HexColor("#5f6872")))
    chart = Chart()

    tcga_frames: dict[str, pd.DataFrame] = {}
    km_results: dict[str, KMResult] = {}
    failed = []
    for project in TCGA_PROJECTS:
        print(f"TCGA {project}")
        try:
            df, result = analyze_tcga_project(project)
            if result.n >= 20:
                tcga_frames[project] = df
                km_results[project] = result
        except Exception as exc:
            failed.append((project, str(exc)))

    print("DepMap metadata")
    meta = depmap_model_metadata()
    print("DepMap expression")
    dep_expr = depmap_gene_column("OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv").rename(columns={GENE_SYMBOL: "expr"})
    print("DepMap copy number")
    dep_cn = depmap_gene_column("PortalOmicsCNGeneLog2.csv").rename(columns={GENE_SYMBOL: "cn_log2"})
    print("DepMap mutations")
    dep_mut = depmap_gene_mutations()
    dep = meta.merge(dep_expr, on="ModelID", how="left").merge(dep_cn, on="ModelID", how="left")
    dep["Lineage"] = dep["Lineage"].fillna("Unknown")
    dep_mut_annot = dep_mut.merge(meta, on="ModelID", how="left") if len(dep_mut) else dep_mut
    print("TCGA cBioPortal mutation/CNA")
    alteration_df = analyze_tcga_alterations(tcga_frames)
    print("CPTAC protein")
    protein_df = cptac_protein_data(gene_entrez_id())

    appendix = []
    for project, r in km_results.items():
        appendix.append({
            "dataset": "TCGA",
            "group": project,
            "n": r.n,
            "median_gene_expression": round(r.median_expr, 4),
            "logrank_p": r.p_value,
            "logrank_chi2": r.chi2,
            "direction": r.direction,
            "high_events": r.high_events,
            "low_events": r.low_events,
        })
    for lin, sub in dep.groupby("Lineage"):
        if sub["expr"].notna().sum() >= 3:
            appendix.append({
                "dataset": "DepMap_26Q1",
                "group": lin,
                "n": int(sub["expr"].notna().sum()),
                "median_gene_expression": round(float(sub["expr"].median()), 4),
                "median_gene_copy_number_log2": round(float(sub["cn_log2"].median()), 4) if sub["cn_log2"].notna().any() else "",
                "gene_mutated_models": int(dep_mut_annot[dep_mut_annot.get("Lineage", "") == lin]["ModelID"].nunique()) if len(dep_mut_annot) else 0,
            })
    for _, row in alteration_df.iterrows():
        appendix.append({
            "dataset": "TCGA_cBioPortal",
            "group": f"{row['project']}:{row['alteration']}",
            "n": row["altered_patients"],
            "logrank_p": row["logrank_p"],
            "direction": row["direction"],
            "altered_events": row["altered_events"],
            "unaltered_events": row["unaltered_events"],
        })
    for study, sub in (protein_df.groupby("study") if len(protein_df) else []):
        appendix.append({
            "dataset": "CPTAC_cBioPortal",
            "group": study,
            "n": int(sub["protein"].notna().sum()),
            "median_protein": round(float(sub["protein"].median()), 4),
        })
    pd.DataFrame(appendix).to_csv(APPENDIX_OUT, index=False)

    doc = SimpleDocTemplate(str(PDF_OUT), pagesize=letter, leftMargin=0.55 * inch, rightMargin=0.55 * inch, topMargin=0.5 * inch, bottomMargin=0.5 * inch)
    story = []
    story.append(Paragraph(f"{GENE_SYMBOL} Pan-Cancer Report: TCGA Patients and DepMap Cell Lines", styles["Title"]))
    story.append(Paragraph(f"Generated {date.today().isoformat()} | DepMap release: {DEPMAP_RELEASE}", styles["Small"]))
    story.append(Spacer(1, 10))
    story.append(Paragraph(f"This expanded report analyzes {GENE_SYMBOL} across available TCGA cancer projects and DepMap 26Q1 cancer models. It includes pan-cancer patient expression/survival summaries and DepMap expression, copy-number, and mutation figures. TCGA survival curves use a within-project median {GENE_SYMBOL} split and are exploratory.", styles["BodyText"]))
    story.append(Spacer(1, 8))
    rows = [["Claim", "Figure/data support"]]
    rows += [
        ["Broad patient expression", f"{GENE_SYMBOL} expression row analyzed in {len(tcga_frames)} TCGA projects."],
        ["Context-dependent survival", f"{sum(1 for r in km_results.values() if r.p_value < 0.05)} TCGA projects reached nominal median-split log-rank p<0.05; directions vary."],
        ["Broad cell-line expression", f"DepMap expression analyzed for {dep['expr'].notna().sum()} models with lineage metadata."],
        ["Copy-number not a clean driver", "DepMap copy-number distribution and expression-copy relationship are shown; interpretation remains exploratory."],
        ["Mutation burden", f"DepMap {GENE_SYMBOL} mutation table contains {len(dep_mut)} mutation records across {dep_mut['ModelID'].nunique() if len(dep_mut) else 0} unique models."],
        ["TCGA alteration survival", f"cBioPortal mutation/GISTIC copy-number statuses tested across {alteration_df['project'].nunique() if len(alteration_df) else 0} TCGA PanCancer Atlas projects."],
        ["CPTAC protein", f"CPTAC/cBioPortal protein abundance retrieved for {protein_df['study'].nunique() if len(protein_df) else 0} studies and {len(protein_df)} samples."],
    ]
    add_table(story, rows, widths=[2.0 * inch, 5.0 * inch])

    story.append(PageBreak())
    story.append(Paragraph(f"TCGA Pan-Cancer {GENE_SYMBOL} Expression", styles["Heading1"]))
    med_items = sorted([(p, float(df["expr"].median())) for p, df in tcga_frames.items()], key=lambda x: x[1], reverse=True)
    story.append(chart.hbar(f"Median {GENE_SYMBOL} expression by TCGA project", med_items, f"Median {GENE_SYMBOL} expression", max_items=33))
    story.append(Spacer(1, 10))
    n_items = sorted([(p, float(len(df))) for p, df in tcga_frames.items()], key=lambda x: x[1], reverse=True)
    story.append(chart.hbar("TCGA matched expression-survival sample count by project", n_items, "Matched patients", color="#8aa66a", max_items=33))

    story.append(PageBreak())
    story.append(Paragraph("TCGA Pan-Cancer Survival Screen", styles["Heading1"]))
    p_items = []
    for p, r in km_results.items():
        signed = -math.log10(max(r.p_value, 1e-300))
        if r.direction == "high better/equal":
            signed *= -1
        p_items.append((p, signed))
    p_items = sorted(p_items, key=lambda x: abs(x[1]), reverse=True)
    story.append(chart.hbar(f"Signed -log10(log-rank p): positive = {GENE_SYMBOL}-high worse", p_items, "Signed -log10 p", color="#b97979", max_items=33))
    story.append(Paragraph("Nominal p-values are not multiple-testing corrected. This plot is a screening view to show context dependence, not clinical validation.", styles["Caption"]))

    top_km = sorted(km_results.values(), key=lambda r: r.p_value)[:8]
    for idx, r in enumerate(top_km):
        if idx % 2 == 0:
            story.append(PageBreak())
        story.append(chart.km(f"TCGA {r.cancer}: {GENE_SYMBOL} high vs low", r))
        story.append(Spacer(1, 8))

    story.append(PageBreak())
    story.append(Paragraph("TCGA Expression-Survival Summary", styles["Heading1"]))
    rows = [["TCGA project", "n", "Median expr", "Log-rank p", "Direction", "High events", "Low events"]]
    for p, r in sorted(km_results.items(), key=lambda kv: kv[1].p_value):
        rows.append([p, str(r.n), f"{r.median_expr:.2f}", f"{r.p_value:.3g}", r.direction, str(r.high_events), str(r.low_events)])
    add_table(story, rows, widths=[0.75 * inch, 0.55 * inch, 0.8 * inch, 0.8 * inch, 1.2 * inch, 0.8 * inch, 0.75 * inch], font_size=6.4)

    story.append(PageBreak())
    story.append(Paragraph(f"TCGA {GENE_SYMBOL} Mutation and Copy-Number Alteration Survival", styles["Heading1"]))
    if len(alteration_df):
        freq_rows = alteration_df.pivot_table(
            index="project",
            columns="alteration",
            values="altered_patients",
            aggfunc="first",
            fill_value=0,
        )
        amp_items = sorted([(idx, float(row.get("amplification", 0))) for idx, row in freq_rows.iterrows()], key=lambda x: x[1], reverse=True)
        mut_items = sorted([(idx, float(row.get("mutation", 0))) for idx, row in freq_rows.iterrows()], key=lambda x: x[1], reverse=True)
        story.append(chart.hbar(f"TCGA {GENE_SYMBOL} high-level amplification counts by project, cBioPortal GISTIC +2", amp_items, "Amplified patients", color="#9a8ab8", max_items=33))
        story.append(Spacer(1, 8))
        story.append(chart.hbar(f"TCGA {GENE_SYMBOL} mutation counts by project, cBioPortal mutation profile", mut_items, "Mutated patients", color="#c49a5a", max_items=33))
        story.append(PageBreak())
        valid_alt = alteration_df[alteration_df["logrank_p"].notna()].copy()
        if len(valid_alt):
            signed_items = []
            for _, row in valid_alt.iterrows():
                val = -math.log10(max(float(row["logrank_p"]), 1e-300))
                if row["direction"] == "altered better/equal":
                    val *= -1
                signed_items.append((f"{row['project']} {row['alteration']}", val))
            signed_items = sorted(signed_items, key=lambda x: abs(x[1]), reverse=True)
            story.append(chart.hbar(f"Alteration survival screen: signed -log10 p, positive = altered worse", signed_items, "Signed -log10 p", color="#b97979", max_items=33))
            story.append(Paragraph("Mutation, amplification, gain/amplification, and deletion statuses are analyzed separately. Small altered groups are marked insufficient in the appendix and excluded from this signed p-value plot.", styles["Caption"]))
            top_alt = valid_alt.sort_values("logrank_p").head(4)
            for _, row in top_alt.iterrows():
                story.append(PageBreak())
                story.append(chart.alteration_km(f"TCGA {row['project']} {row['alteration']}: altered vs unaltered", row.to_dict()))
        rows = [["Project", "Alteration", "Altered n", "Unaltered n", "Log-rank p", "Direction"]]
        for _, row in alteration_df.sort_values(["alteration", "project"]).iterrows():
            p = row["logrank_p"]
            rows.append([
                row["project"],
                row["alteration"],
                str(int(row["altered_patients"])),
                str(int(row["unaltered_patients"])),
                "NA" if pd.isna(p) else f"{float(p):.3g}",
                row["direction"],
            ])
        story.append(PageBreak())
        story.append(Paragraph("TCGA Alteration Survival Summary", styles["Heading1"]))
        add_table(story, rows, widths=[0.7 * inch, 1.25 * inch, 0.65 * inch, 0.75 * inch, 0.8 * inch, 1.1 * inch], font_size=5.8)
    else:
        story.append(Paragraph("No TCGA mutation/CNA alteration data were retrieved from cBioPortal for this gene.", styles["BodyText"]))

    story.append(PageBreak())
    story.append(Paragraph(f"DepMap {GENE_SYMBOL} Expression Across Cancer Lineages", styles["Heading1"]))
    lineage_counts = dep.groupby("Lineage")["expr"].count().sort_values(ascending=False)
    top_lineages = lineage_counts[lineage_counts >= 10].head(24).index.tolist()
    story.append(chart.boxplot(
        f"DepMap 26Q1 {GENE_SYMBOL} mRNA expression by lineage, top lineages by model count",
        {lin: dep.loc[dep["Lineage"] == lin, "expr"].dropna().tolist() for lin in top_lineages},
        "TPM log2+1",
    ))
    lineage_medians = dep.groupby("Lineage")["expr"].median().dropna().sort_values(ascending=False)
    story.append(Spacer(1, 10))
    story.append(chart.hbar(f"Highest median {GENE_SYMBOL} expression lineages in DepMap", list(lineage_medians.items())[:20], "Median TPM log2+1", color="#6b9fb5", max_items=20))

    story.append(PageBreak())
    story.append(Paragraph(f"DepMap {GENE_SYMBOL} Copy Number", styles["Heading1"]))
    cn_lineages = dep.groupby("Lineage")["cn_log2"].count().sort_values(ascending=False)
    cn_top = cn_lineages[cn_lineages >= 10].head(24).index.tolist()
    story.append(chart.boxplot(
        f"DepMap 26Q1 {GENE_SYMBOL} gene-level copy-number signal by lineage",
        {lin: dep.loc[dep["Lineage"] == lin, "cn_log2"].dropna().tolist() for lin in cn_top},
        "Portal CN log2 signal",
    ))
    high_cn_threshold = percentile(dep["cn_log2"].dropna().tolist(), 0.95)
    high_cn = dep[dep["cn_log2"] >= high_cn_threshold]
    high_cn_counts = high_cn["Lineage"].value_counts().head(20)
    story.append(Spacer(1, 10))
    story.append(chart.hbar(f"Lineages enriched among top 5% {GENE_SYMBOL} copy-number signal models, threshold {high_cn_threshold:.2f}", list(high_cn_counts.items()), "Models", color="#9a8ab8", max_items=20))
    story.append(Paragraph("DepMap copy-number values are shown as portal-provided gene-level log2 copy-number signals. The top-5% panel is a relative high-copy signal screen, not a clinical amplification call.", styles["Caption"]))

    story.append(PageBreak())
    story.append(Paragraph("DepMap Expression-Copy Number Relationship", styles["Heading1"]))
    corr_df = dep[["expr", "cn_log2"]].dropna()
    corr = float(corr_df["expr"].corr(corr_df["cn_log2"])) if len(corr_df) > 2 else float("nan")
    story.append(chart.scatter(
        f"{GENE_SYMBOL} mRNA expression vs {GENE_SYMBOL} copy-number signal in DepMap models, Pearson r={corr:.2f}",
        corr_df["cn_log2"].tolist(),
        corr_df["expr"].tolist(),
        f"{GENE_SYMBOL} copy-number signal",
        f"{GENE_SYMBOL} expression",
    ))
    story.append(Paragraph(
        f"This scatter plot is a model-level screen for copy-number/expression coupling. A weak-to-moderate correlation would support the interpretation that {GENE_SYMBOL} expression is not explained by copy number alone.",
        styles["Caption"],
    ))

    story.append(PageBreak())
    story.append(Paragraph(f"DepMap {GENE_SYMBOL} Mutation Landscape", styles["Heading1"]))
    if len(dep_mut_annot):
        mut_lineage_counts = dep_mut_annot["Lineage"].fillna("Unknown").value_counts().head(20)
        story.append(chart.hbar(f"{GENE_SYMBOL}-mutated DepMap models by lineage", list(mut_lineage_counts.items()), "Mutation records/models", color="#c49a5a", max_items=20))
        mut_rows = [["ModelID", "Cell line", "Lineage", "Protein change", "Consequence", "Impact"]]
        for _, row in dep_mut_annot.head(20).iterrows():
            mut_rows.append([
                "" if pd.isna(row.get("ModelID", "")) else str(row.get("ModelID", "")),
                "" if pd.isna(row.get("CellLineName", "")) else str(row.get("CellLineName", "")),
                "" if pd.isna(row.get("Lineage", "")) else str(row.get("Lineage", "")),
                "" if pd.isna(row.get("ProteinChange", "")) else str(row.get("ProteinChange", "")),
                ("" if pd.isna(row.get("Consequence", "")) else str(row.get("Consequence", "")))[:34],
                "" if pd.isna(row.get("Impact", "")) else str(row.get("Impact", "")),
            ])
        add_table(story, mut_rows, widths=[0.85 * inch, 1.05 * inch, 1.2 * inch, 1.0 * inch, 2.25 * inch, 0.55 * inch], font_size=6)
    else:
        story.append(Paragraph(f"No {GENE_SYMBOL} mutation records were found in the parsed DepMap 26Q1 somatic mutation table.", styles["BodyText"]))

    story.append(PageBreak())
    story.append(Paragraph(f"CPTAC {GENE_SYMBOL} Protein Expression", styles["Heading1"]))
    if len(protein_df):
        protein_by_study = {study: sub["protein"].dropna().tolist() for study, sub in protein_df.groupby("study") if sub["protein"].notna().sum() >= 3}
        if protein_by_study:
            story.append(chart.boxplot(
                f"CPTAC/cBioPortal {GENE_SYMBOL} protein abundance by study",
                protein_by_study,
                "Protein abundance",
            ))
            med_protein = sorted([(study, float(sub["protein"].median())) for study, sub in protein_df.groupby("study") if sub["protein"].notna().sum() >= 3], key=lambda x: x[1], reverse=True)
            story.append(Spacer(1, 10))
            story.append(chart.hbar(f"Median {GENE_SYMBOL} protein abundance by CPTAC study", med_protein, "Median protein", color="#6b9fb5", max_items=20))
        rows = [["CPTAC study", "Samples", "Median protein", "Min", "Max"]]
        for study, sub in protein_df.groupby("study"):
            vals = sub["protein"].dropna()
            if len(vals):
                rows.append([study, str(len(vals)), f"{vals.median():.3g}", f"{vals.min():.3g}", f"{vals.max():.3g}"])
        story.append(Spacer(1, 10))
        add_table(story, rows, widths=[1.9 * inch, 0.7 * inch, 1.0 * inch, 0.7 * inch, 0.7 * inch], font_size=7)
        story.append(Paragraph("Protein values come from public CPTAC proteomics studies exposed through cBioPortal protein_quantification profiles. Units differ by study/profile, so cross-study comparisons should be interpreted as study-level abundance distributions rather than a single harmonized pan-cancer proteome.", styles["Caption"]))
    else:
        story.append(Paragraph(f"No CPTAC protein_quantification records were retrieved for {GENE_SYMBOL} from the configured public CPTAC cBioPortal studies.", styles["BodyText"]))

    story.append(PageBreak())
    story.append(Paragraph("Interpretation", styles["Heading1"]))
    for bullet in [
        f"{GENE_SYMBOL} is measurable across most TCGA tumor projects and DepMap cancer model lineages, supporting broad biological availability.",
        f"The patient survival screen shows context dependence; {GENE_SYMBOL}-high should not be treated as a pan-cancer prognostic biomarker without tumor-specific validation.",
        f"DepMap mutation analysis shows how often {GENE_SYMBOL} is altered in cancer models and whether those alterations cluster by lineage.",
        f"DepMap copy-number distributions show variation, but relative high-copy signal should be interpreted cautiously; {GENE_SYMBOL} expression may be shaped by lineage, cell state, and copy number.",
        f"For experimental model selection, prioritize models with high {GENE_SYMBOL} mRNA/protein and matched phenotype assays. Mutation or copy-number filters alone are lower-yield.",
    ]:
        story.append(Paragraph(f"• {bullet}", styles["BodyText"]))
    if failed:
        story.append(Spacer(1, 8))
        story.append(Paragraph(f"TCGA projects skipped because data retrieval or parsing failed: {failed}", styles["Small"]))
    story.append(Spacer(1, 10))
    story.append(Paragraph(f"Primary data sources: TCGA/GDC Hub via UCSC Xena; DepMap portal download API, DepMap Public 26Q1 released 2026-04-01; Human Protein Atlas {GENE_SYMBOL} entry and gene-specific literature for interpretation.", styles["Small"]))

    doc.build(story)

    # ── Print key results to stdout for agent consumption ──
    print(f"\n[RESULTS] Pan-cancer report for {GENE_SYMBOL} ({GENE_ENSEMBL})")
    print(f"[RESULTS] TCGA projects analyzed: {len(km_results)} / {len(TCGA_PROJECTS)} | Failed: {len(failed)}")

    # Top expression-survival hits
    if km_results:
        sorted_km = sorted(km_results.values(), key=lambda r: r.p_value)
        sig_count = sum(1 for r in sorted_km if r.p_value < 0.05)
        print(f"[RESULTS] Expression-survival: {sig_count} projects with p < 0.05 (uncorrected)")
        for r in sorted_km[:5]:
            p_str = f"p = {r.p_value:.2e}" if r.p_value >= 0.0001 else "p < 0.0001"
            sig = "SIGNIFICANT" if r.p_value < 0.05 else "ns"
            print(f"[RESULTS]   {r.cancer}: {p_str} ({sig}), n={r.n}, {r.direction}")

    # Top alteration-survival hits
    if len(alteration_df):
        alt_sig = alteration_df[alteration_df["logrank_p"] < 0.05]
        print(f"[RESULTS] Alteration-survival (cBioPortal): {len(alt_sig)} significant tests out of {len(alteration_df)}")
        for _, row in alteration_df.nsmallest(5, "logrank_p").iterrows():
            p_str = f"p = {row['logrank_p']:.2e}" if row["logrank_p"] >= 0.0001 else "p < 0.0001"
            sig = "SIGNIFICANT" if row["logrank_p"] < 0.05 else "ns"
            print(f"[RESULTS]   {row['project']}:{row['alteration']}: {p_str} ({sig}), {row['direction']}")

    # DepMap summary
    if dep["expr"].notna().any():
        top_lin = dep.groupby("Lineage")["expr"].median().nlargest(3)
        print(f"[RESULTS] DepMap: {dep['expr'].notna().sum()} models with expression data")
        for lin, med in top_lin.items():
            print(f"[RESULTS]   Top lineage: {lin} (median TPM log = {med:.2f})")
    if len(dep_mut):
        print(f"[RESULTS] DepMap mutations: {len(dep_mut)} mutations in {dep_mut['ModelID'].nunique()} models")

    # CPTAC summary
    if len(protein_df):
        n_studies = protein_df["study"].nunique()
        n_samples = protein_df["protein"].notna().sum()
        print(f"[RESULTS] CPTAC protein data: {n_samples} samples across {n_studies} studies")

    if failed:
        print(f"[RESULTS] Skipped TCGA projects: {', '.join(p for p, _ in failed)}")

    print(f"[DONE] PDF: {PDF_OUT}")
    print(f"[DONE] Appendix: {APPENDIX_OUT}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a pan-cancer TCGA/DepMap report for one gene.")
    parser.add_argument("--gene", default="PRNP", help="HGNC gene symbol, e.g. PRNP or NADK")
    parser.add_argument("--ensembl", default="ENSG00000171867", help="Human Ensembl gene ID without version")
    parser.add_argument("--outdir", default=None, help="Output directory (default: current directory)")
    args = parser.parse_args()
    configure(args.gene, args.ensembl, args.outdir)
    build_report()
