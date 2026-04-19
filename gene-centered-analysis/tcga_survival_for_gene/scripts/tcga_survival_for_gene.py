#!/usr/bin/env python3
"""
TCGA Survival Analysis for a Single Gene

Stratifies TCGA patients by gene expression (high vs low) within a given
cohort and performs Kaplan-Meier survival analysis for Overall Survival (OS)
and/or Disease-Free Survival (DFS).

Data source: GDC REST API  https://api.gdc.cancer.gov
Dependencies: lifelines, matplotlib, pandas, requests, numpy
"""

import argparse
import io
import os
from typing import Dict, List, Optional, Tuple

import matplotlib
import sys as _sys, os as _os
try:
    _sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), *(['..'] * 3), '_shared'))
    from plot_style import init_style
except ImportError:
    def init_style(**kw): pass  # graceful fallback if _shared not available
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test


GDC_BASE = "https://api.gdc.cancer.gov"

# =========================================================
# Utilities
# =========================================================

def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def chunked(seq: List, n: int) -> List[List]:
    return [seq[i : i + n] for i in range(0, len(seq), n)]


def safe_float(val) -> Optional[float]:
    try:
        f = float(val)
        return None if np.isnan(f) else f
    except (TypeError, ValueError):
        return None


# =========================================================
# GDC client
# =========================================================

class GDCClient:
    def __init__(self, base_url: str = GDC_BASE, timeout: int = 120):
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": "tcga-survival-for-gene/0.1"})

    def _post_json(self, endpoint: str, payload: dict, accept: str = "application/json"):
        url = f"{self.base_url}{endpoint}"
        headers = {"Content-Type": "application/json", "Accept": accept}
        r = self.session.post(url, json=payload, headers=headers, timeout=self.timeout)
        r.raise_for_status()
        return r

    def resolve_gene(self, gene_symbol: str) -> Tuple[str, str]:
        """Resolve gene symbol → GDC gene_id (Ensembl, no version)."""
        payload = {
            "filters": {
                "op": "=",
                "content": {"field": "symbol", "value": gene_symbol},
            },
            "fields": "gene_id,symbol,gene_type",
            "format": "JSON",
            "size": 10,
        }
        r = self._post_json("/genes", payload)
        hits = r.json().get("data", {}).get("hits", [])
        if not hits:
            raise ValueError(f"Gene symbol not found in GDC: {gene_symbol}")

        exact_pc = [
            h for h in hits
            if str(h.get("symbol", "")).upper() == gene_symbol.upper()
            and str(h.get("gene_type", "")).lower() == "protein_coding"
        ]
        hit = exact_pc[0] if exact_pc else (
            next((h for h in hits if str(h.get("symbol", "")).upper() == gene_symbol.upper()), hits[0])
        )
        return hit["gene_id"], hit["symbol"]

    def get_cases_with_clinical(self, project_id: str) -> pd.DataFrame:
        """
        Fetch all cases for a TCGA project with clinical survival fields.

        Returns a DataFrame with columns:
          case_id, submitter_id, vital_status, days_to_death,
          days_to_last_follow_up, days_to_recurrence, progression_or_recurrence,
          days_to_new_tumor_event, new_tumor_event_type
        """
        # NOTE: In the modern GDC schema, `vital_status` and `days_to_death`
        # live under the `demographic` node, NOT `diagnoses`. We request both
        # locations and prefer `demographic` first, falling back to
        # `diagnoses` for older records / schema variants.
        #
        # For DFS/PFI: recurrence data lives under `diagnoses.treatments`
        # and `follow_ups` (new_tumor_event fields), NOT under
        # `diagnoses.days_to_recurrence` which is almost always NULL.
        payload = {
            "filters": {
                "op": "=",
                "content": {"field": "project.project_id", "value": project_id},
            },
            "fields": ",".join([
                "case_id",
                "submitter_id",
                "demographic.vital_status",
                "demographic.days_to_death",
                "diagnoses.vital_status",
                "diagnoses.days_to_death",
                "diagnoses.days_to_last_follow_up",
                "diagnoses.days_to_recurrence",
                "diagnoses.progression_or_recurrence",
                # Follow-up / new tumor event fields for DFS
                "diagnoses.treatments.days_to_treatment_start",
                "diagnoses.treatments.treatment_or_therapy",
                "follow_ups.days_to_follow_up",
                "follow_ups.progression_or_recurrence",
                "follow_ups.disease_response",
                "follow_ups.molecular_tests.days_to_molecular_test",
            ]),
            "format": "JSON",
            "size": 10000,
        }
        r = self._post_json("/cases", payload)
        hits = r.json().get("data", {}).get("hits", [])

        def _first(*vals):
            """Return the first non-empty value from candidates."""
            for v in vals:
                if v is None:
                    continue
                if isinstance(v, str) and not v.strip():
                    continue
                return v
            return None

        rows = []
        for h in hits:
            case_id = h.get("case_id") or h.get("id")
            submitter_id = h.get("submitter_id")
            demo = h.get("demographic") or {}
            diags = h.get("diagnoses") or []
            # Use first diagnosis record (most patients have one)
            d = diags[0] if diags else {}
            vital = _first(demo.get("vital_status"), d.get("vital_status"))
            d_death = _first(demo.get("days_to_death"), d.get("days_to_death"))

            # --- Extract new tumor / recurrence event from follow_ups ---
            follow_ups = h.get("follow_ups") or []
            days_to_new_tumor = None
            new_tumor_type = None
            fu_progression = None

            for fu in follow_ups:
                # Check follow_up-level progression_or_recurrence
                fu_prog = str(fu.get("progression_or_recurrence") or "").strip().lower()
                fu_response = str(fu.get("disease_response") or "").strip().lower()
                fu_days = safe_float(fu.get("days_to_follow_up"))

                is_progression = (
                    fu_prog in ("yes", "yes, progression", "yes, recurrence")
                    or "progressive" in fu_response
                    or "recurrence" in fu_response
                )

                if is_progression and fu_days is not None:
                    if days_to_new_tumor is None or fu_days < days_to_new_tumor:
                        days_to_new_tumor = fu_days
                        new_tumor_type = fu_prog or fu_response or "yes"
                        fu_progression = "yes"

            # Also check diagnoses-level recurrence as fallback
            diag_recurrence = safe_float(d.get("days_to_recurrence"))
            diag_prog = d.get("progression_or_recurrence")
            if diag_recurrence is not None and (
                days_to_new_tumor is None or diag_recurrence < days_to_new_tumor
            ):
                days_to_new_tumor = diag_recurrence
                new_tumor_type = "recurrence"
                fu_progression = diag_prog or "yes"

            rows.append(
                {
                    "case_id": case_id,
                    "submitter_id": submitter_id,
                    "vital_status": vital,
                    "days_to_death": safe_float(d_death),
                    "days_to_last_follow_up": safe_float(d.get("days_to_last_follow_up")),
                    "days_to_recurrence": safe_float(d.get("days_to_recurrence")),
                    "progression_or_recurrence": fu_progression or diag_prog,
                    "days_to_new_tumor_event": days_to_new_tumor,
                    "new_tumor_event_type": new_tumor_type,
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
        """Query /gene_expression/values for one gene across many cases."""
        if not case_ids:
            return pd.DataFrame(columns=["case_id", "expression"])

        all_rows = []
        for batch in chunked(case_ids, batch_size):
            payload = {
                "case_ids": batch,
                "gene_ids": [gene_id],
                "tsv_units": tsv_units,
                "format": "tsv",
            }
            r = self._post_json(
                "/gene_expression/values", payload, accept="text/tab-separated-values"
            )
            text = r.text.strip()
            if not text:
                continue
            df = pd.read_csv(io.StringIO(text), sep="\t")
            if df.empty:
                continue
            row = df.iloc[0]
            for col in df.columns[1:]:
                val = safe_float(row[col])
                if val is not None:
                    all_rows.append({"case_id": col, "expression": val})

        if not all_rows:
            return pd.DataFrame(columns=["case_id", "expression"])

        out = pd.DataFrame(all_rows).drop_duplicates(subset=["case_id"]).reset_index(drop=True)
        return out


# =========================================================
# Survival data construction
# =========================================================

def build_os_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add OS_time and OS_event columns.

    OS_event : 1 = death observed, 0 = censored
    OS_time  : days_to_death (if dead) else days_to_last_follow_up
    """
    out = df.copy()
    events, times = [], []
    for _, row in out.iterrows():
        vs = str(row.get("vital_status") or "").lower()
        dead = vs == "dead"
        if dead and row.get("days_to_death") is not None:
            events.append(1)
            times.append(row["days_to_death"])
        elif row.get("days_to_last_follow_up") is not None:
            events.append(0)
            times.append(row["days_to_last_follow_up"])
        else:
            events.append(None)
            times.append(None)
    out["OS_event"] = events
    out["OS_time"] = times
    return out


def build_dfs_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add DFS_time and DFS_event columns.

    DFS (Disease-Free Survival) / PFI (Progression-Free Interval):
      DFS_event : 1 = recurrence/progression/death observed, 0 = censored
      DFS_time  : earliest of days_to_new_tumor_event, days_to_recurrence,
                  or days_to_death; falls back to days_to_last_follow_up
                  when censored.

    Uses new_tumor_event data from GDC follow_ups (primary) and
    diagnoses.days_to_recurrence (fallback). This produces DFS curves
    that differ from OS because recurrence/progression events are
    captured independently of death.
    """
    out = df.copy()
    events, times = [], []
    for _, row in out.iterrows():
        vs = str(row.get("vital_status") or "").lower()
        prog = str(row.get("progression_or_recurrence") or "").lower()

        # Collect candidate event times from all available sources
        d_new_tumor = row.get("days_to_new_tumor_event")
        d_recurrence = row.get("days_to_recurrence")
        d_death = row.get("days_to_death")
        d_follow_up = row.get("days_to_last_follow_up")

        # Determine if a DFS event occurred
        has_recurrence = (
            prog in ("yes", "yes, progression", "yes, recurrence")
            or (d_new_tumor is not None and pd.notna(d_new_tumor))
            or (d_recurrence is not None and pd.notna(d_recurrence))
        )
        is_dead = (vs == "dead")
        has_event = has_recurrence or is_dead

        if has_event:
            # Collect all valid event times
            candidate_times = [
                t for t in [d_new_tumor, d_recurrence, d_death]
                if t is not None and pd.notna(t)
            ]
            if candidate_times:
                events.append(1)
                times.append(min(candidate_times))  # earliest event
            elif d_follow_up is not None and pd.notna(d_follow_up):
                # Event known but no time → use follow-up as approximation
                events.append(1)
                times.append(d_follow_up)
            else:
                events.append(None)
                times.append(None)
        else:
            # Censored: no recurrence and alive
            if d_follow_up is not None and pd.notna(d_follow_up):
                events.append(0)
                times.append(d_follow_up)
            else:
                events.append(None)
                times.append(None)

    out["DFS_event"] = events
    out["DFS_time"] = times
    return out


# =========================================================
# Stratification
# =========================================================

def stratify_expression(
    df: pd.DataFrame,
    method: str,
    custom_cutoff: Optional[float],
) -> pd.DataFrame:
    """
    Add expression_group column: 'High' or 'Low'.

    method:
      median   - split at median (default)
      quartile - top 25% = High, bottom 25% = Low (middle excluded)
      custom   - split at custom_cutoff value
    """
    out = df.copy()
    expr = out["expression"].dropna()

    if method == "median":
        cutoff = float(expr.median())
        out["expression_group"] = out["expression"].apply(
            lambda x: "High" if pd.notna(x) and x >= cutoff else ("Low" if pd.notna(x) else None)
        )
        out.attrs["cutoff"] = cutoff
        out.attrs["cutoff_method"] = f"median ({cutoff:.4f})"

    elif method == "quartile":
        q25 = float(expr.quantile(0.25))
        q75 = float(expr.quantile(0.75))
        def assign(x):
            if pd.isna(x):
                return None
            if x >= q75:
                return "High"
            if x <= q25:
                return "Low"
            return None  # middle-excluded
        out["expression_group"] = out["expression"].apply(assign)
        out.attrs["cutoff"] = (q25, q75)
        out.attrs["cutoff_method"] = f"quartile (Q1={q25:.4f}, Q3={q75:.4f})"

    elif method == "custom":
        if custom_cutoff is None:
            raise ValueError("--custom-cutoff must be provided when --stratify custom is used.")
        out["expression_group"] = out["expression"].apply(
            lambda x: "High" if pd.notna(x) and x >= custom_cutoff else ("Low" if pd.notna(x) else None)
        )
        out.attrs["cutoff"] = custom_cutoff
        out.attrs["cutoff_method"] = f"custom ({custom_cutoff:.4f})"

    else:
        raise ValueError(f"Unknown stratification method: {method}")

    return out


# =========================================================
# Kaplan-Meier plotting
# =========================================================

PALETTE = {"High": "#D62728", "Low": "#1F77B4"}


def _km_ax(
    df: pd.DataFrame,
    time_col: str,
    event_col: str,
    ax: plt.Axes,
    gene: str,
    cohort: str,
    endpoint_label: str,
    cutoff_method: str,
) -> Dict:
    """
    Fit and plot KM curves for High vs Low groups on ax.
    Returns dict with n_high, n_low, logrank_p.
    """
    valid = df[[time_col, event_col, "expression_group"]].dropna()
    valid = valid[valid[time_col] > 0]
    valid = valid[valid["expression_group"].isin(["High", "Low"])]

    if valid.empty or valid["expression_group"].nunique() < 2:
        ax.text(0.5, 0.5, "Insufficient data", transform=ax.transAxes, ha="center")
        return {"n_high": 0, "n_low": 0, "logrank_p": None}

    groups = {}
    kmf = KaplanMeierFitter()
    for grp in ["High", "Low"]:
        sub = valid[valid["expression_group"] == grp]
        if sub.empty:
            continue
        kmf.fit(
            sub[time_col],
            event_observed=sub[event_col],
            label=f"{grp} (n={len(sub)})",
        )
        kmf.plot_survival_function(
            ax=ax,
            ci_show=True,
            color=PALETTE[grp],
            ci_alpha=0.12,
        )
        groups[grp] = sub

    # Log-rank test
    logrank_p = None
    if "High" in groups and "Low" in groups:
        result = logrank_test(
            groups["High"][time_col],
            groups["Low"][time_col],
            event_observed_A=groups["High"][event_col],
            event_observed_B=groups["Low"][event_col],
        )
        logrank_p = float(result.p_value)
        p_str = f"p = {logrank_p:.4f}" if logrank_p >= 0.0001 else f"p < 0.0001"
        ax.text(
            0.98, 0.98, p_str,
            transform=ax.transAxes,
            ha="right", va="top",
            fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="gray", alpha=0.8),
        )

    ax.set_title(f"{gene} – {endpoint_label} ({cohort})", fontsize=13, pad=8)
    ax.set_xlabel("Time (days)", fontsize=11)
    ax.set_ylabel("Survival Probability", fontsize=11)
    ax.set_ylim(0, 1.05)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(fontsize=9, loc="lower left")
    ax.text(
        0.02, 0.02,
        f"Stratification: {cutoff_method}",
        transform=ax.transAxes,
        fontsize=7,
        color="gray",
    )

    n_high = len(groups.get("High", pd.DataFrame()))
    n_low = len(groups.get("Low", pd.DataFrame()))
    return {"n_high": n_high, "n_low": n_low, "logrank_p": logrank_p}


def plot_km(
    df: pd.DataFrame,
    outdir: str,
    gene: str,
    cohort: str,
    modes: List[str],
    cutoff_method: str,
) -> Dict:
    """
    Plot KM curves for requested endpoints. Returns summary stats.
    Saves PNG + PDF for each endpoint.
    """
    stats = {}
    for mode in modes:
        if mode == "os":
            time_col, event_col = "OS_time", "OS_event"
            label = "Overall Survival"
            suffix = "os_km"
        else:
            time_col, event_col = "DFS_time", "DFS_event"
            label = "Disease-Free Survival"
            suffix = "dfs_km"

        if time_col not in df.columns or event_col not in df.columns:
            continue

        fig, ax = plt.subplots(figsize=(7, 5), dpi=300)
        result = _km_ax(df, time_col, event_col, ax, gene, cohort, label, cutoff_method)
        plt.tight_layout()

        base = os.path.join(outdir, f"{gene}.{cohort}.{suffix}")
        fig.savefig(f"{base}.png", dpi=300, bbox_inches="tight")
        fig.savefig(f"{base}.pdf", bbox_inches="tight")
        plt.close(fig)

        stats[mode] = result

    return stats


# =========================================================
# Main
# =========================================================

def main():
    parser = argparse.ArgumentParser(
        description="TCGA survival analysis for a single gene (KM + log-rank).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--gene", required=True, help="Gene symbol, e.g. TP53, EGFR")
    parser.add_argument(
        "--cancer-type", required=True,
        help="TCGA cohort code, e.g. BRCA, LUAD, COAD (with or without 'TCGA-' prefix)",
    )
    parser.add_argument(
        "--mode",
        choices=["os", "dfs", "both"],
        default="os",
        help="Survival endpoint(s) to analyse",
    )
    parser.add_argument(
        "--stratify",
        choices=["median", "quartile", "custom"],
        default="median",
        help="Patient stratification method based on gene expression",
    )
    parser.add_argument(
        "--custom-cutoff",
        type=float,
        default=None,
        help="Expression cutoff when --stratify custom is used",
    )
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

    # Normalise cancer type
    cohort = args.cancer_type.upper()
    if cohort.startswith("TCGA-"):
        project_id = cohort
        cohort_code = cohort.replace("TCGA-", "")
    else:
        project_id = f"TCGA-{cohort}"
        cohort_code = cohort

    # Determine requested endpoints
    if args.mode == "both":
        modes = ["os", "dfs"]
    else:
        modes = [args.mode]

    # -------------------------------------------------------
    print(f"[INFO] Resolving gene: {args.gene}")
    client = GDCClient()
    gene_id, resolved_symbol = client.resolve_gene(args.gene)
    print(f"[INFO] Resolved → {resolved_symbol} ({gene_id})")

    # -------------------------------------------------------
    print(f"[INFO] Fetching clinical data for {project_id}")
    clinical_df = client.get_cases_with_clinical(project_id)
    if clinical_df.empty:
        raise ValueError(f"No clinical cases found for {project_id}.")
    print(f"[INFO] Clinical records: {len(clinical_df)}")

    # -------------------------------------------------------
    print(f"[INFO] Fetching expression data for {resolved_symbol} in {project_id}")
    case_ids = clinical_df["case_id"].dropna().tolist()
    expr_df = client.get_expression_values(case_ids, gene_id, tsv_units=args.tsv_units)
    if expr_df.empty:
        raise ValueError(f"No expression values retrieved for {resolved_symbol} in {project_id}.")
    print(f"[INFO] Expression records: {len(expr_df)}")

    # -------------------------------------------------------
    # Merge clinical + expression
    merged = clinical_df.merge(expr_df, on="case_id", how="inner")
    if merged.empty:
        raise ValueError("No cases matched between clinical and expression data.")
    print(f"[INFO] Matched cases: {len(merged)}")

    # Build survival columns
    merged = build_os_columns(merged)
    merged = build_dfs_columns(merged)

    # Diagnostic: report DFS data availability
    n_with_new_tumor = int(merged["days_to_new_tumor_event"].notna().sum()) if "days_to_new_tumor_event" in merged.columns else 0
    n_with_recurrence = int(merged["days_to_recurrence"].notna().sum()) if "days_to_recurrence" in merged.columns else 0
    n_dfs_events = int((merged["DFS_event"] == 1).sum()) if "DFS_event" in merged.columns else 0
    n_os_events = int((merged["OS_event"] == 1).sum()) if "OS_event" in merged.columns else 0
    n_dfs_differs = 0
    if "DFS_time" in merged.columns and "OS_time" in merged.columns:
        both_valid = merged[["DFS_time", "OS_time"]].dropna()
        n_dfs_differs = int((both_valid["DFS_time"] != both_valid["OS_time"]).sum())
    print(f"[INFO] Recurrence data: {n_with_new_tumor} from follow_ups, {n_with_recurrence} from diagnoses")
    print(f"[INFO] OS events: {n_os_events}, DFS events: {n_dfs_events}")
    print(f"[INFO] Cases where DFS_time ≠ OS_time: {n_dfs_differs}")

    # Stratify
    merged = stratify_expression(merged, args.stratify, args.custom_cutoff)
    cutoff_method = merged.attrs.get("cutoff_method", args.stratify)

    # -------------------------------------------------------
    # Save per-case data table
    out_cols = [
        "case_id", "submitter_id", "vital_status",
        "days_to_death", "days_to_last_follow_up",
        "days_to_recurrence", "progression_or_recurrence",
        "days_to_new_tumor_event", "new_tumor_event_type",
        "expression", "expression_group",
        "OS_time", "OS_event",
        "DFS_time", "DFS_event",
    ]
    out_cols = [c for c in out_cols if c in merged.columns]
    data_path = os.path.join(args.outdir, f"{resolved_symbol}.{cohort_code}.survival_data.tsv")
    merged[out_cols].to_csv(data_path, sep="\t", index=False)
    print(f"[INFO] Survival data written: {data_path}")

    # -------------------------------------------------------
    # KM plots
    km_stats = plot_km(merged, args.outdir, resolved_symbol, cohort_code, modes, cutoff_method)
    print(f"[INFO] KM plots generated for modes: {modes}")

    # -------------------------------------------------------
    # Summary
    expr_vals = merged["expression"].dropna()
    summary: Dict = {
        "query_gene": args.gene,
        "resolved_symbol": resolved_symbol,
        "gene_id": gene_id,
        "cancer_type": cohort_code,
        "project_id": project_id,
        "mode": args.mode,
        "stratification": cutoff_method,
        "tsv_units": args.tsv_units,
        "total_matched_cases": int(len(merged)),
        "cases_with_expression": int(expr_vals.count()),
        "expression_median": float(expr_vals.median()) if len(expr_vals) else None,
        "expression_mean": float(expr_vals.mean()) if len(expr_vals) else None,
        "n_high_group": int((merged["expression_group"] == "High").sum()),
        "n_low_group": int((merged["expression_group"] == "Low").sum()),
    }

    for endpoint, stats in km_stats.items():
        prefix = endpoint.upper()
        summary[f"{prefix}_n_high"] = stats.get("n_high")
        summary[f"{prefix}_n_low"] = stats.get("n_low")
        summary[f"{prefix}_logrank_p"] = stats.get("logrank_p")

    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(
        os.path.join(args.outdir, "survival_summary.tsv"), sep="\t", index=False
    )

    with open(os.path.join(args.outdir, "summary.txt"), "w", encoding="utf-8") as f:
        for k, v in summary.items():
            f.write(f"{k}: {v}\n")

    print(f"[DONE] Results written to: {args.outdir}")


if __name__ == "__main__":
    main()
