#!/usr/bin/env python3
"""
TCGA Download Data - Download genomic datasets from GDC Data Portal API.

Downloads gene expression, mutations, CNV, clinical, or methylation data
from one or more TCGA cancer types, with optional merging and manifest generation.
"""

import argparse
import gzip
import json
import os
import sys
import time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from urllib.parse import urljoin

import requests
import pandas as pd
import numpy as np


GDC_API_BASE = "https://api.gdc.cancer.gov"
GDC_FILES_ENDPOINT = f"{GDC_API_BASE}/files"
GDC_CASES_ENDPOINT = f"{GDC_API_BASE}/cases"


SESSION = requests.Session()
SESSION.headers.update({
    "User-Agent": "tcga-download/1.0",
    "Accept": "application/json",
})

RETRY_ATTEMPTS = 3
RETRY_DELAY = 2  # seconds


@dataclass
class FileRecord:
    gdc_id: str
    filename: str
    size_bytes: int
    md5: str
    status: str = "pending"
    merge_group: str = ""
    local_path: str = ""
    error_msg: str = ""


# Map data type to GDC query filters
DATA_TYPE_FILTERS = {
    "expression": {
        "data_category": "Transcriptome Profiling",
        "data_type": "Gene Expression Quantification",
        "workflow_type": "STAR - Counts",
        "file_format": "TSV",
    },
    "mutations": {
        "data_category": "Simple Nucleotide Variation",
        "data_type": "Masked Somatic Mutation",
        "file_format": "MAF",
    },
    "cnv": {
        "data_category": "Copy Number Variation",
        "data_type": "Gene Level Copy Number",
        "file_format": "TSV",
    },
    "clinical": {
        "data_category": "Clinical",
        "data_type": "Clinical Supplement",
        "file_format": "XML",
    },
    "methylation": {
        "data_category": "DNA Methylation",
        "data_type": "Methylation Beta-Value",
        "file_format": "TXT",
    },
}

VALID_CANCER_TYPES = [
    "TCGA-BRCA", "TCGA-LUAD", "TCGA-LUSC", "TCGA-OV", "TCGA-UCEC",
    "TCGA-COAD", "TCGA-READ", "TCGA-PRAD", "TCGA-HNSC", "TCGA-THCA",
    "TCGA-GBM", "TCGA-LGG", "TCGA-SKCM", "TCGA-BLCA", "TCGA-KIRC",
    "TCGA-KIRP", "TCGA-LAML", "TCGA-LIHC", "TCGA-PAAD", "TCGA-SARC",
    "TCGA-STAD", "TCGA-TGCT", "TCGA-UCS", "TCGA-DLBC", "TCGA-ESCA",
    "TCGA-GBMLGG", "TCGA-MESO", "TCGA-PCPG", "TCGA-THYM", "TCGA-UVM",
    "TCGA-ACC", "TCGA-CHOL",
]


def ensure_dir(path: str) -> None:
    """Create directory if it does not exist."""
    os.makedirs(path, exist_ok=True)


def query_gdc_api(endpoint: str, filters: Dict, expand: str = "", size: int = 200) -> Dict:
    """
    Query GDC API with given filters.

    Args:
        endpoint: API endpoint URL
        filters: dictionary of filter key-value pairs
        expand: comma-separated fields to expand
        size: number of results to return

    Returns:
        JSON response dictionary
    """
    params = {
        "format": "JSON",
        "size": size,
    }

    if filters:
        # Build filter JSON
        filter_list = []
        for key, value in filters.items():
            filter_list.append({
                "op": "=",
                "content": {
                    "field": key,
                    "value": value,
                }
            })
        if filter_list:
            params["filters"] = json.dumps({
                "op": "and",
                "content": filter_list,
            })

    if expand:
        params["expand"] = expand

    for attempt in range(RETRY_ATTEMPTS):
        try:
            response = SESSION.get(endpoint, params=params, timeout=30)
            response.raise_for_status()
            return response.json()
        except requests.RequestException as e:
            if attempt < RETRY_ATTEMPTS - 1:
                print(f"  [Retry {attempt + 1}/{RETRY_ATTEMPTS}] API request failed: {e}")
                time.sleep(RETRY_DELAY)
            else:
                raise


def get_available_projects() -> List[str]:
    """Get list of available TCGA project codes."""
    print("Querying GDC for available TCGA projects...")
    try:
        response = query_gdc_api(GDC_CASES_ENDPOINT, {"project.program.name": "TCGA"}, size=50)
        projects = set()
        for case in response.get("data", {}).get("hits", []):
            project = case.get("project", {}).get("project_id")
            if project:
                projects.add(project)
        return sorted(list(projects))
    except Exception as e:
        print(f"Error querying projects: {e}")
        return VALID_CANCER_TYPES


def query_files_by_type(data_type: str, cancer_type: str, max_files: int) -> List[FileRecord]:
    """
    Query GDC for files matching data type and cancer type.

    Args:
        data_type: one of expression, mutations, cnv, clinical, methylation
        cancer_type: TCGA cancer type code (e.g., TCGA-BRCA)
        max_files: maximum number of files to retrieve

    Returns:
        list of FileRecord objects
    """
    if data_type not in DATA_TYPE_FILTERS:
        raise ValueError(f"Unknown data type: {data_type}")

    filters_spec = DATA_TYPE_FILTERS[data_type]

    # Build filter for GDC API
    filters = {
        "project.project_id": cancer_type,
        "data_category": filters_spec["data_category"],
        "data_type": filters_spec["data_type"],
    }

    if "workflow_type" in filters_spec:
        filters["workflow_type"] = filters_spec["workflow_type"]

    print(f"  Querying GDC for {data_type} files in {cancer_type}...")

    try:
        response = query_gdc_api(GDC_FILES_ENDPOINT, filters, size=max_files)
    except Exception as e:
        print(f"    Error querying GDC: {e}")
        return []

    files_data = response.get("data", {}).get("hits", [])
    print(f"    Found {len(files_data)} files")

    records = []
    for f in files_data:
        record = FileRecord(
            gdc_id=f.get("id", ""),
            filename=f.get("file_name", ""),
            size_bytes=f.get("file_size", 0),
            md5=f.get("md5sum", ""),
            merge_group=data_type,
        )
        records.append(record)

    return records


def download_file(gdc_id: str, filename: str, output_path: str, md5_expected: str = "") -> Tuple[bool, str]:
    """
    Download a file from GDC using file streaming.

    Args:
        gdc_id: GDC file UUID
        filename: original filename for logging
        output_path: full path to save file
        md5_expected: expected MD5 checksum for verification

    Returns:
        (success, error_message)
    """
    download_url = f"{GDC_API_BASE}/data/{gdc_id}"

    for attempt in range(RETRY_ATTEMPTS):
        try:
            response = SESSION.get(download_url, stream=True, timeout=120)
            response.raise_for_status()

            # Stream to disk
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            with open(output_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)

            # Verify MD5 if expected
            if md5_expected:
                import hashlib
                md5_actual = hashlib.md5()
                with open(output_path, 'rb') as f:
                    for chunk in iter(lambda: f.read(8192), b''):
                        md5_actual.update(chunk)
                if md5_actual.hexdigest() != md5_expected:
                    os.remove(output_path)
                    raise ValueError(f"MD5 mismatch: expected {md5_expected}, got {md5_actual.hexdigest()}")

            return True, ""

        except Exception as e:
            if attempt < RETRY_ATTEMPTS - 1:
                print(f"      [Retry {attempt + 1}/{RETRY_ATTEMPTS}] Download failed: {e}")
                time.sleep(RETRY_DELAY)
            else:
                return False, str(e)

    return False, "Max retries exceeded"


def merge_expression_files(raw_dir: str, merged_dir: str, file_format: str) -> Optional[str]:
    """
    Merge expression TSV files into a single gene × sample matrix.

    Args:
        raw_dir: directory containing downloaded TSV files
        merged_dir: output directory for merged file
        file_format: TSV or CSV

    Returns:
        path to merged file, or None if failed
    """
    print("  Merging expression files...")
    os.makedirs(merged_dir, exist_ok=True)

    tsv_files = list(Path(raw_dir).glob("*.tsv")) + list(Path(raw_dir).glob("*.txt"))
    if not tsv_files:
        print("    No expression files found to merge")
        return None

    try:
        dfs = []
        for fpath in tsv_files[:10]:  # Limit to first 10 for memory
            df = pd.read_csv(fpath, sep='\t', index_col=0, nrows=None)
            dfs.append(df)

        # Merge on gene ID
        merged = pd.concat(dfs, axis=1)
        merged = merged.loc[~merged.index.duplicated(keep='first')]

        sep = '\t' if file_format.upper() == 'TSV' else ','
        ext = 'tsv' if file_format.upper() == 'TSV' else 'csv'
        out_path = os.path.join(merged_dir, f"expression_matrix.{ext}")
        merged.to_csv(out_path, sep=sep)

        print(f"    Merged {len(dfs)} files into {merged.shape[0]} genes × {merged.shape[1]} samples")
        return out_path

    except Exception as e:
        print(f"    Error merging: {e}")
        return None


def merge_mutation_files(raw_dir: str, merged_dir: str) -> Optional[str]:
    """
    Merge MAF mutation files into a single file.

    Args:
        raw_dir: directory containing downloaded MAF files
        merged_dir: output directory for merged file

    Returns:
        path to merged file, or None if failed
    """
    print("  Merging mutation files...")
    os.makedirs(merged_dir, exist_ok=True)

    maf_files = list(Path(raw_dir).glob("*.maf"))
    if not maf_files:
        print("    No mutation files found to merge")
        return None

    try:
        dfs = []
        for fpath in maf_files[:10]:
            # MAF files have header with # comments
            df = pd.read_csv(fpath, sep='\t', comment='#', low_memory=False)
            dfs.append(df)

        merged = pd.concat(dfs, axis=0, ignore_index=True)

        out_path = os.path.join(merged_dir, "combined_mutations.maf")
        merged.to_csv(out_path, sep='\t', index=False)

        print(f"    Merged {len(dfs)} files into {len(merged)} mutations")
        return out_path

    except Exception as e:
        print(f"    Error merging: {e}")
        return None


def merge_cnv_files(raw_dir: str, merged_dir: str, file_format: str) -> Optional[str]:
    """
    Merge CNV TSV files into a single gene × sample matrix.

    Args:
        raw_dir: directory containing downloaded TSV files
        merged_dir: output directory for merged file
        file_format: TSV or CSV

    Returns:
        path to merged file, or None if failed
    """
    print("  Merging CNV files...")
    os.makedirs(merged_dir, exist_ok=True)

    tsv_files = list(Path(raw_dir).glob("*.tsv")) + list(Path(raw_dir).glob("*.txt"))
    if not tsv_files:
        print("    No CNV files found to merge")
        return None

    try:
        dfs = []
        for fpath in tsv_files[:10]:
            df = pd.read_csv(fpath, sep='\t', index_col=0)
            dfs.append(df)

        merged = pd.concat(dfs, axis=1)
        merged = merged.loc[~merged.index.duplicated(keep='first')]

        sep = '\t' if file_format.upper() == 'TSV' else ','
        ext = 'tsv' if file_format.upper() == 'TSV' else 'csv'
        out_path = os.path.join(merged_dir, f"cnv_matrix.{ext}")
        merged.to_csv(out_path, sep=sep)

        print(f"    Merged {len(dfs)} files into {merged.shape[0]} genes × {merged.shape[1]} samples")
        return out_path

    except Exception as e:
        print(f"    Error merging: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(
        description="Download TCGA genomic data from GDC Data Portal"
    )
    parser.add_argument(
        "--data-type",
        choices=["expression", "mutations", "cnv", "clinical", "methylation"],
        required=True,
        help="Type of genomic data to download"
    )
    parser.add_argument(
        "--cancer-types",
        type=str,
        default="",
        help="Comma-separated TCGA cancer type codes (e.g., TCGA-BRCA,TCGA-LUAD) or 'all'"
    )
    parser.add_argument(
        "--outdir",
        type=str,
        required=True,
        help="Output directory for downloads and merged files"
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=5,
        help="Maximum number of files to download per cancer type (default: 5)"
    )
    parser.add_argument(
        "--file-format",
        choices=["TSV", "CSV"],
        default="TSV",
        help="Output format for merged matrices (default: TSV)"
    )
    parser.add_argument(
        "--manifest",
        type=str,
        default="",
        help="Optional path to write manifest JSON"
    )
    parser.add_argument(
        "--list-projects",
        action="store_true",
        help="List available TCGA projects and exit"
    )

    args = parser.parse_args()

    # Handle list projects
    if args.list_projects:
        projects = get_available_projects()
        print("Available TCGA cancer types:")
        for p in projects:
            print(f"  {p}")
        return 0

    # Parse cancer types
    if args.cancer_types.lower() == "all":
        cancer_types = get_available_projects()
    else:
        cancer_types = [c.strip().upper() for c in args.cancer_types.split(",") if c.strip()]

    if not cancer_types:
        print("Error: no cancer types specified. Use --cancer-types or --list-projects")
        return 1

    # Create directories
    raw_dir = os.path.join(args.outdir, "raw")
    merged_dir = os.path.join(args.outdir, "merged")
    ensure_dir(args.outdir)
    ensure_dir(raw_dir)
    ensure_dir(merged_dir)

    print(f"\nTCGA Download: {args.data_type}")
    print(f"Cancer types: {', '.join(cancer_types[:3])}{'...' if len(cancer_types) > 3 else ''}")
    print(f"Max files per type: {args.max_files}")
    print(f"Output: {args.outdir}\n")

    # Query and download
    all_records: List[FileRecord] = []
    total_downloaded = 0
    total_failed = 0

    for cancer_type in cancer_types:
        print(f"Processing {cancer_type}...")
        records = query_files_by_type(args.data_type, cancer_type, args.max_files)

        for record in records:
            out_path = os.path.join(raw_dir, record.filename)

            # Skip if already exists
            if os.path.exists(out_path):
                print(f"    {record.filename} (exists, skipped)")
                record.status = "exists"
                record.local_path = out_path
            else:
                success, error = download_file(record.gdc_id, record.filename, out_path, record.md5)
                if success:
                    print(f"    {record.filename} (downloaded)")
                    record.status = "downloaded"
                    record.local_path = out_path
                    total_downloaded += 1
                else:
                    print(f"    {record.filename} (FAILED: {error})")
                    record.status = "failed"
                    record.error_msg = error
                    total_failed += 1

            all_records.append(record)

    # Merge if requested
    merged_files = []
    if args.data_type == "expression" and total_downloaded > 0:
        merged_path = merge_expression_files(raw_dir, merged_dir, args.file_format)
        if merged_path:
            merged_files.append(merged_path)
    elif args.data_type == "mutations" and total_downloaded > 0:
        merged_path = merge_mutation_files(raw_dir, merged_dir)
        if merged_path:
            merged_files.append(merged_path)
    elif args.data_type == "cnv" and total_downloaded > 0:
        merged_path = merge_cnv_files(raw_dir, merged_dir, args.file_format)
        if merged_path:
            merged_files.append(merged_path)

    # Write manifest
    if args.manifest:
        manifest = {
            "data_type": args.data_type,
            "cancer_types": cancer_types,
            "download_timestamp": pd.Timestamp.now().isoformat(),
            "outdir": os.path.abspath(args.outdir),
            "files": [asdict(r) for r in all_records],
            "merged_files": merged_files,
            "summary": {
                "requested": len(all_records),
                "downloaded": total_downloaded,
                "existing": len([r for r in all_records if r.status == "exists"]),
                "failed": total_failed,
                "merged": len(merged_files) > 0,
            }
        }
        with open(args.manifest, 'w') as f:
            json.dump(manifest, f, indent=2)
        print(f"\nManifest written to {args.manifest}")

    # Summary
    print(f"\nSummary:")
    print(f"  Files found: {len(all_records)}")
    print(f"  Downloaded: {total_downloaded}")
    print(f"  Existing: {len([r for r in all_records if r.status == 'exists'])}")
    print(f"  Failed: {total_failed}")
    print(f"  Merged files: {len(merged_files)}")

    return 0 if total_failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
