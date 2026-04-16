#!/usr/bin/env python3

import argparse
import io
import math
import os
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), '..', '..', '..', '_shared'))
from plot_style import init_style
import pandas as pd
import requests


GDC_BASE = "https://api.gdc.cancer.gov"


# =========================================================
# Utilities
# =========================================================
def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def chunked(seq: Sequence[str], n: int) -> List[List[str]]:
    return [list(seq[i:i + n]) for i in range(0, len(seq), n)]


def median_or_none(series: pd.Series) -> Optional[float]:
    s = pd.to_numeric(series, errors="coerce").dropna()
    return float(s.median()) if len(s) else None


def mean_or_none(series: pd.Series) -> Optional[float]:
    s = pd.to_numeric(series, errors="coerce").dropna()
    return float(s.mean()) if len(s) else None


def clean_none_rows(df: pd.DataFrame, col: str) -> pd.DataFrame:
    out = df.copy()
    out[col] = pd.to_numeric(out[col], errors="coerce")
    return out.dropna(subset=[col]).reset_index(drop=True)


# =========================================================
# GDC client
# =========================================================
class GDCClient:
    def __init__(self, base_url: str = GDC_BASE, timeout: int = 120):
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": "tcga-expression-for-gene/0.1"})

    def _post_json(self, endpoint: str, payload: dict, accept: str = "application/json"):
        url = f"{self.base_url}{endpoint}"
        headers = {"Content-Type": "application/json", "Accept": accept}
        r = self.session.post(url, json=payload, headers=headers, timeout=self.timeout)
        r.raise_for_status()
        return r

    def _get_json(self, endpoint: str, params: Optional[dict] = None):
        url = f"{self.base_url}{endpoint}"
        r = self.session.get(url, params=params, timeout=self.timeout)
        r.raise_for_status()
        return r.json()

    def resolve_gene(self, gene_symbol: str) -> Tuple[str, str]:
        """
        Resolve gene symbol to GDC gene_id (Ensembl without version).
        """
        payload = {
            "filters": {
                "op": "=",
                "content": {
                    "field": "symbol",
                    "value": gene_symbol
                }
            },
            "fields": "gene_id,symbol,gene_type",
            "format": "JSON",
            "size": 10
        }
        r = self._post_json("/genes", payload)
        data = r.json().get("data", {})
        hits = data.get("hits", [])
        if not hits:
            raise ValueError(f"Could not resolve gene symbol in GDC: {gene_symbol}")

        # Prefer exact symbol + protein_coding
        exact_pc = [
            h for h in hits
            if str(h.get("symbol", "")).upper() == gene_symbol.upper()
            and str(h.get("gene_type", "")).lower() == "protein_coding"
        ]
        if exact_pc:
            hit = exact_pc[0]
        else:
            exact = [h for h in hits if str(h.get("symbol", "")).upper() == gene_symbol.upper()]
            hit = exact[0] if exact else hits[0]

        gene_id = hit["gene_id"]
        symbol = hit["symbol"]
        return gene_id, symbol

    def get_tcga_projects(self) -> List[str]:
        """
        Get TCGA project IDs using wildcard support on project_id.
        GDC docs describe wildcard use in filter values.  [oai_citation:1‡GDC Docs](https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/)
        """
        payload = {
            "filters": {
                "op": "=",
                "content": {
                    "field": "project_id",
                    "value": "TCGA-*"
                }
            },
            "fields": "project_id",
            "format": "JSON",
            "size": 200
        }
        r = self._post_json("/projects", payload)
        hits = r.json().get("data", {}).get("hits", [])
        projects = sorted({h["project_id"] for h in hits if h.get("project_id")})
        if not projects:
            raise ValueError("No TCGA projects found from GDC /projects endpoint.")
        return projects

    def get_cases_for_project(
        self,
        project_id: str,
        sample_types: Optional[List[str]] = None,
    ) -> pd.DataFrame:
        """
        Fetch cases for a TCGA project, optionally restricted by sample type.
        """
        filters = {
            "op": "and",
            "content": [
                {
                    "op": "=",
                    "content": {
                        "field": "project.project_id",
                        "value": project_id
                    }
                }
            ]
        }

        if sample_types:
            filters["content"].append(
                {
                    "op": "in",
                    "content": {
                        "field": "samples.sample_type",
                        "value": sample_types
                    }
                }
            )

        payload = {
            "filters": filters,
            "fields": "case_id,submitter_id,project.project_id,samples.sample_type",
            "format": "JSON",
            "size": 10000
        }

        r = self._post_json("/cases", payload)
        hits = r.json().get("data", {}).get("hits", [])

        rows = []
        for h in hits:
            case_id = h.get("case_id") or h.get("id")
            submitter_id = h.get("submitter_id")
            proj = h.get("project", {})
            proj_id = proj.get("project_id", project_id)
            samples = h.get("samples", []) or []
            sample_type_values = sorted({s.get("sample_type") for s in samples if s.get("sample_type")})
            rows.append(
                {
                    "case_id": case_id,
                    "submitter_id": submitter_id,
                    "project_id": proj_id,
                    "sample_types": ";".join(sample_type_values) if sample_type_values else None,
                }
            )

        return pd.DataFrame(rows)

    def get_expression_values(
        self,
        case_ids: List[str],
        gene_id: str,
        tsv_units: str = "uqfpkm",
        batch_size: int = 2000,
    ) -> pd.DataFrame:
        """
        Query /gene_expression/values for one gene across many cases.
        GDC docs specify case_ids + gene_ids and tsv_units in the POST body,
        returning TSV. Supported units: uqfpkm or median_centered_log2_uqfpkm.  [oai_citation:2‡GDC Docs](https://docs.gdc.cancer.gov/API/Users_Guide/Data_Analysis/)
        """
        if not case_ids:
            return pd.DataFrame(columns=["case_id", "expression"])

        all_rows = []
        for batch in chunked(case_ids, batch_size):
            payload = {
                "case_ids": batch,
                "gene_ids": [gene_id],
                "tsv_units": tsv_units,
                "format": "tsv"
            }
            r = self._post_json("/gene_expression/values", payload, accept="text/tab-separated-values")
            text = r.text.strip()
            if not text:
                continue

            df = pd.read_csv(io.StringIO(text), sep="\t")
            if df.empty:
                continue

            # Expected format: gene_id | case1 | case2 | ...
            gene_col = df.columns[0]
            row = df.iloc[0]
            for col in df.columns[1:]:
                all_rows.append(
                    {
                        "case_id": col,
                        "expression": pd.to_numeric(row[col], errors="coerce"),
                    }
                )

        out = pd.DataFrame(all_rows)
        if out.empty:
            return pd.DataFrame(columns=["case_id", "expression"])

        out = out.dropna(subset=["expression"]).drop_duplicates(subset=["case_id"]).reset_index(drop=True)
        return out


# =========================================================
# Analysis
# =========================================================
TUMOR_SAMPLE_TYPES = [
    "Primary Tumor",
    "Recurrent Tumor",
    "Metastatic",
    "Additional Metastatic",
    "Additional - New Primary",
]

NORMAL_SAMPLE_TYPES = [
    "Solid Tissue Normal",
    "Blood Derived Normal",
    "Buccal Cell Normal",
]


def run_pan_cancer(
    client: GDCClient,
    gene_id: str,
    top_n: int,
    tsv_units: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    projects = client.get_tcga_projects()
    rows = []

    for project_id in projects:
        try:
            case_df = client.get_cases_for_project(project_id)
            if case_df.empty:
                continue

            expr_df = client.get_expression_values(case_df["case_id"].dropna().tolist(), gene_id, tsv_units=tsv_units)
            if expr_df.empty:
                continue

            merged = case_df.merge(expr_df, on="case_id", how="inner")
            if merged.empty:
                continue

            rows.append(
                {
                    "project_id": project_id,
                    "cancer_type": project_id.replace("TCGA-", ""),
                    "median_expression": median_or_none(merged["expression"]),
                    "mean_expression": mean_or_none(merged["expression"]),
                    "num_cases": int(len(merged)),
                }
            )
        except Exception:
            continue

    full_df = pd.DataFrame(rows)
    full_df = full_df.dropna(subset=["median_expression"]).sort_values("median_expression", ascending=False).reset_index(drop=True)

    if full_df.empty:
        raise ValueError("No pan-cancer expression results retrieved from GDC.")

    top_df = full_df.head(top_n).copy()
    return top_df, full_df


def run_single_cohort(
    client: GDCClient,
    gene_id: str,
    cancer_type: str,
    tsv_units: str,
) -> pd.DataFrame:
    project_id = cancer_type if cancer_type.startswith("TCGA-") else f"TCGA-{cancer_type.upper()}"
    case_df = client.get_cases_for_project(project_id)
    if case_df.empty:
        raise ValueError(f"No cases found for cohort {project_id}.")

    expr_df = client.get_expression_values(case_df["case_id"].dropna().tolist(), gene_id, tsv_units=tsv_units)
    if expr_df.empty:
        raise ValueError(f"No expression values retrieved for cohort {project_id}.")

    merged = case_df.merge(expr_df, on="case_id", how="inner")
    if merged.empty:
        raise ValueError(f"No matched expression values for cohort {project_id}.")

    merged["cancer_type"] = project_id.replace("TCGA-", "")
    return merged.reset_index(drop=True)


def run_tumor_vs_normal(
    client: GDCClient,
    gene_id: str,
    cancer_type: str,
    tsv_units: str,
) -> pd.DataFrame:
    project_id = cancer_type if cancer_type.startswith("TCGA-") else f"TCGA-{cancer_type.upper()}"

    tumor_cases = client.get_cases_for_project(project_id, sample_types=TUMOR_SAMPLE_TYPES)
    normal_cases = client.get_cases_for_project(project_id, sample_types=NORMAL_SAMPLE_TYPES)

    if tumor_cases.empty or normal_cases.empty:
        raise ValueError(f"Tumor vs normal comparison unavailable for cohort {project_id}.")

    # Remove overlapping case_ids to reduce ambiguity at case level
    tumor_ids = set(tumor_cases["case_id"].dropna())
    normal_ids = set(normal_cases["case_id"].dropna())
    overlap = tumor_ids & normal_ids
    tumor_ids = sorted(tumor_ids - overlap)
    normal_ids = sorted(normal_ids - overlap)

    if not tumor_ids or not normal_ids:
        raise ValueError(f"Tumor vs normal comparison unavailable for cohort {project_id} after removing overlapping case IDs.")

    tumor_expr = client.get_expression_values(tumor_ids, gene_id, tsv_units=tsv_units)
    normal_expr = client.get_expression_values(normal_ids, gene_id, tsv_units=tsv_units)

    if tumor_expr.empty or normal_expr.empty:
        raise ValueError(f"No tumor/normal expression values retrieved for cohort {project_id}.")

    tumor_df = tumor_cases[tumor_cases["case_id"].isin(tumor_expr["case_id"])].merge(tumor_expr, on="case_id", how="inner")
    tumor_df["group"] = "Tumor"

    normal_df = normal_cases[normal_cases["case_id"].isin(normal_expr["case_id"])].merge(normal_expr, on="case_id", how="inner")
    normal_df["group"] = "Normal"

    out = pd.concat([tumor_df, normal_df], ignore_index=True)
    out["cancer_type"] = project_id.replace("TCGA-", "")

    if out.empty:
        raise ValueError(f"Tumor vs normal comparison unavailable for cohort {project_id}.")

    return out.reset_index(drop=True)


# =========================================================
# Plotting
# =========================================================
def make_pan_cancer_barplot(df: pd.DataFrame, out_png: str, gene: str, units: str) -> None:
    if df.empty:
        return

    plot_df = df.copy().sort_values("median_expression", ascending=True)
    plt.figure(figsize=(10, max(5, 0.35 * len(plot_df) + 1.5)), dpi=300)
    ax = plt.gca()
    ax.barh(plot_df["cancer_type"], plot_df["median_expression"], edgecolor="black", linewidth=0.4)
    ax.set_title(f"{gene} TCGA pan-cancer expression", fontsize=14, pad=10)
    ax.set_xlabel(units, fontsize=12)
    ax.set_ylabel("Cancer type", fontsize=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", linestyle=":", linewidth=0.5, alpha=0.5)
    ax.set_axisbelow(True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(os.path.splitext(out_png)[0] + ".pdf", bbox_inches="tight")
    plt.close()


def make_single_cohort_boxplot(df: pd.DataFrame, out_png: str, gene: str, cohort_code: str, units: str) -> None:
    if df.empty:
        return

    plt.figure(figsize=(4, 6), dpi=300)
    ax = plt.gca()
    ax.boxplot(df["expression"].dropna().values, vert=True)
    ax.set_title(f"{gene} expression in {cohort_code}", fontsize=14, pad=10)
    ax.set_ylabel(units, fontsize=12)
    ax.set_xticks([1])
    ax.set_xticklabels([cohort_code])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(os.path.splitext(out_png)[0] + ".pdf", bbox_inches="tight")
    plt.close()


def make_tumor_vs_normal_boxplot(df: pd.DataFrame, out_png: str, gene: str, cohort_code: str, units: str) -> None:
    if df.empty:
        return

    tumor = df.loc[df["group"] == "Tumor", "expression"].dropna().values
    normal = df.loc[df["group"] == "Normal", "expression"].dropna().values

    plt.figure(figsize=(5, 6), dpi=300)
    ax = plt.gca()
    ax.boxplot([tumor, normal], labels=["Tumor", "Normal"])
    ax.set_title(f"{gene} tumor vs normal in {cohort_code}", fontsize=14, pad=10)
    ax.set_ylabel(units, fontsize=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(os.path.splitext(out_png)[0] + ".pdf", bbox_inches="tight")
    plt.close()


# =========================================================
# Main
# =========================================================
def main():
    parser = argparse.ArgumentParser(description="TCGA single-gene expression analysis via GDC API.")
    parser.add_argument("--gene", required=True, help="Gene symbol, e.g. TP53")
    parser.add_argument(
        "--mode",
        choices=["pan_cancer", "single_cohort", "tumor_vs_normal"],
        default="pan_cancer",
        help="Analysis mode",
    )
    parser.add_argument("--cancer-type", help="TCGA cohort code, e.g. BRCA, LUAD")
    parser.add_argument("--top-n", type=int, default=10, help="Top cohorts to report in pan-cancer mode")
    parser.add_argument(
        "--tsv-units",
        choices=["uqfpkm", "median_centered_log2_uqfpkm"],
        default="uqfpkm",
        help="GDC gene expression units",
    )
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    init_style(
        font_family=getattr(args, 'font_family', None),
        font_size=getattr(args, 'font_size', None),
    )

    ensure_dir(args.outdir)

    client = GDCClient()
    gene_id, resolved_symbol = client.resolve_gene(args.gene)

    summary: Dict[str, object] = {
        "query_gene": args.gene,
        "resolved_symbol": resolved_symbol,
        "gene_id": gene_id,
        "mode": args.mode,
        "tsv_units": args.tsv_units,
    }

    if args.mode == "pan_cancer":
        top_df, full_df = run_pan_cancer(client, gene_id, args.top_n, args.tsv_units)
        full_df.to_csv(os.path.join(args.outdir, f"{resolved_symbol}.pan_cancer_expression.tsv"), sep="\t", index=False)
        make_pan_cancer_barplot(
            top_df,
            os.path.join(args.outdir, f"{resolved_symbol}.pan_cancer_barplot.png"),
            resolved_symbol,
            args.tsv_units,
        )
        summary["num_cohorts"] = int(len(full_df))
        summary["top_cohort"] = full_df.iloc[0]["cancer_type"] if len(full_df) else None
        summary["top_cohort_median_expression"] = float(full_df.iloc[0]["median_expression"]) if len(full_df) else None

    elif args.mode == "single_cohort":
        if not args.cancer_type:
            raise ValueError("--cancer-type is required for single_cohort mode.")
        cohort_df = run_single_cohort(client, gene_id, args.cancer_type, args.tsv_units)
        cohort_code = cohort_df["cancer_type"].iloc[0]
        cohort_df.to_csv(os.path.join(args.outdir, f"{resolved_symbol}.{cohort_code}.expression.tsv"), sep="\t", index=False)
        make_single_cohort_boxplot(
            cohort_df,
            os.path.join(args.outdir, f"{resolved_symbol}.{cohort_code}.expression_boxplot.png"),
            resolved_symbol,
            cohort_code,
            args.tsv_units,
        )
        summary["cancer_type"] = cohort_code
        summary["num_cases"] = int(len(cohort_df))
        summary["median_expression"] = median_or_none(cohort_df["expression"])

    elif args.mode == "tumor_vs_normal":
        if not args.cancer_type:
            raise ValueError("--cancer-type is required for tumor_vs_normal mode.")
        tvn_df = run_tumor_vs_normal(client, gene_id, args.cancer_type, args.tsv_units)
        cohort_code = tvn_df["cancer_type"].iloc[0]
        tvn_df.to_csv(os.path.join(args.outdir, f"{resolved_symbol}.{cohort_code}.tumor_vs_normal.tsv"), sep="\t", index=False)
        make_tumor_vs_normal_boxplot(
            tvn_df,
            os.path.join(args.outdir, f"{resolved_symbol}.{cohort_code}.tumor_vs_normal_boxplot.png"),
            resolved_symbol,
            cohort_code,
            args.tsv_units,
        )
        summary["cancer_type"] = cohort_code
        summary["num_tumor"] = int((tvn_df["group"] == "Tumor").sum())
        summary["num_normal"] = int((tvn_df["group"] == "Normal").sum())
        summary["tumor_median_expression"] = median_or_none(tvn_df.loc[tvn_df["group"] == "Tumor", "expression"])
        summary["normal_median_expression"] = median_or_none(tvn_df.loc[tvn_df["group"] == "Normal", "expression"])

    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(os.path.join(args.outdir, "tcga_expression_summary.tsv"), sep="\t", index=False)

    with open(os.path.join(args.outdir, "summary.txt"), "w", encoding="utf-8") as f:
        for k, v in summary.items():
            f.write(f"{k}: {v}\n")

    print(f"[DONE] Results written to: {args.outdir}")


if __name__ == "__main__":
    main()