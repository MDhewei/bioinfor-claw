#!/usr/bin/env python3
"""
GEO Download Data - Download datasets from NCBI Gene Expression Omnibus.

Downloads expression matrices and metadata from GEO series (GSE) or datasets (GDS)
using the NCBI FTP and eUtils API.
"""

import argparse
import gzip
import json
import os
import sys
import time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.parse import urljoin
from io import StringIO

import requests
import pandas as pd
import numpy as np


NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GEO_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo"

SESSION = requests.Session()
SESSION.headers.update({
    "User-Agent": "geo-download/1.0",
    "Accept": "application/json,text/plain",
})

RETRY_ATTEMPTS = 3
RETRY_DELAY = 2  # seconds


@dataclass
class FileRecord:
    accession: str
    record_type: str  # "series" or "dataset"
    filename: str
    size_bytes: int
    status: str = "pending"
    local_path: str = ""
    error_msg: str = ""


def ensure_dir(path: str) -> None:
    """Create directory if it does not exist."""
    os.makedirs(path, exist_ok=True)


def get_accession_type(accession: str) -> Optional[str]:
    """Determine if accession is GSE, GDS, or GSM."""
    accession = accession.upper().strip()
    if accession.startswith("GSE"):
        return "GSE"
    elif accession.startswith("GDS"):
        return "GDS"
    elif accession.startswith("GSM"):
        return "GSM"
    return None


def query_eutils(db: str, query: str, rettype: str = "json") -> Dict:
    """
    Query NCBI eUtils.

    Args:
        db: database (gds, gse, etc.)
        query: search query
        rettype: return type (json, xml, etc.)

    Returns:
        parsed response
    """
    url = f"{NCBI_EUTILS_BASE}/esearch.fcgi"
    params = {
        "db": db,
        "term": query,
        "retmode": rettype,
        "retmax": 1,
    }

    for attempt in range(RETRY_ATTEMPTS):
        try:
            response = SESSION.get(url, params=params, timeout=30)
            response.raise_for_status()
            return response.json()
        except requests.RequestException as e:
            if attempt < RETRY_ATTEMPTS - 1:
                print(f"  [Retry {attempt + 1}/{RETRY_ATTEMPTS}] eUtils query failed: {e}")
                time.sleep(RETRY_DELAY)
            else:
                raise


def get_series_summary(accession: str) -> Dict:
    """
    Fetch series metadata from eUtils.

    Args:
        accession: GSE accession

    Returns:
        dict with series information
    """
    try:
        url = f"{NCBI_EUTILS_BASE}/esummary.fcgi"
        params = {
            "db": "gds",
            "term": accession,
            "retmode": "json",
        }

        response = SESSION.get(url, params=params, timeout=30)
        response.raise_for_status()
        data = response.json()

        result = {
            "accession": accession,
            "title": "",
            "summary": "",
            "organism": "",
            "platform": "",
            "samples": [],
            "n_samples": 0,
        }

        # Parse eUtils response
        if "result" in data and isinstance(data["result"], dict):
            for key, value in data["result"].items():
                if isinstance(value, dict):
                    if "title" in value:
                        result["title"] = value.get("title", "")
                    if "summary" in value:
                        result["summary"] = value.get("summary", "")
                    if "organism" in value:
                        result["organism"] = value.get("organism", "")

        return result

    except Exception as e:
        print(f"  Warning: Could not fetch series summary: {e}")
        return {
            "accession": accession,
            "title": "",
            "summary": "",
            "organism": "",
            "platform": "",
            "samples": [],
            "n_samples": 0,
        }


def construct_ftp_url(accession: str) -> str:
    """
    Construct FTP URL for GEO series_matrix.txt.gz file.

    GEO FTP structure:
    /geo/series/GSEnnn/GSExxxxx/matrix/GSExxxxx_series_matrix.txt.gz

    Where nnn = first 3 digits + "nnn"
    Example: GSE12345 -> GSE12nnn

    Args:
        accession: GSE accession (e.g., GSE12345)

    Returns:
        full HTTPS URL to matrix file
    """
    accession = accession.upper().strip()

    # Extract numeric part and build prefix
    if accession.startswith("GSE"):
        numeric_part = accession[3:]
        # Get first 2-3 digits, then append 'nnn'
        if len(numeric_part) >= 3:
            prefix = numeric_part[:len(numeric_part)-3] + "nnn"
        else:
            prefix = "nnn"

        base_path = f"series/GSE{prefix}/{accession}/matrix"
        filename = f"{accession}_series_matrix.txt.gz"
        url = f"{GEO_FTP_BASE}/{base_path}/{filename}"
        return url
    else:
        raise ValueError(f"Unsupported accession type: {accession}")


def download_file(url: str, output_path: str) -> Tuple[bool, str, int]:
    """
    Download a file with streaming.

    Args:
        url: download URL
        output_path: path to save file

    Returns:
        (success, error_message, file_size)
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    for attempt in range(RETRY_ATTEMPTS):
        try:
            response = SESSION.get(url, stream=True, timeout=120)
            response.raise_for_status()

            file_size = 0
            with open(output_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
                        file_size += len(chunk)

            return True, "", file_size

        except requests.RequestException as e:
            if attempt < RETRY_ATTEMPTS - 1:
                print(f"      [Retry {attempt + 1}/{RETRY_ATTEMPTS}] Download failed: {e}")
                time.sleep(RETRY_DELAY)
            else:
                return False, str(e), 0

    return False, "Max retries exceeded", 0


def parse_series_matrix(filepath: str) -> Tuple[pd.DataFrame, Dict, List[str]]:
    """
    Parse GEO series_matrix.txt file (optionally gzipped).

    Format:
    - Lines starting with ! are metadata
    - Section !series_matrix_table_begin to !series_matrix_table_end contains matrix
    - Columns after ID_REF are sample IDs

    Args:
        filepath: path to series_matrix.txt or series_matrix.txt.gz

    Returns:
        (expression_dataframe, metadata_dict, sample_list)
    """
    # Detect if gzipped
    if filepath.endswith('.gz'):
        f = gzip.open(filepath, 'rt', encoding='utf-8', errors='ignore')
    else:
        f = open(filepath, 'r', encoding='utf-8', errors='ignore')

    metadata = {
        "series_title": "",
        "series_summary": "",
        "series_organism": "",
        "series_overall_design": "",
        "series_pubmed_id": "",
        "series_contact_name": "",
        "series_contact_email": "",
    }
    sample_metadata = {}
    samples = []
    expression_lines = []
    in_matrix = False

    with f:
        for line in f:
            line = line.rstrip('\n')

            # Skip comments and empty lines
            if line.startswith('#') or not line.strip():
                continue

            # Parse metadata lines
            if line.startswith('!Series_'):
                parts = line.split('=', 1)
                if len(parts) == 2:
                    key = parts[0].strip().lower()
                    value = parts[1].strip().strip('"')
                    if "title" in key:
                        metadata["series_title"] = value
                    elif "summary" in key:
                        metadata["series_summary"] = value
                    elif "organism" in key:
                        metadata["series_organism"] = value
                    elif "overall_design" in key:
                        metadata["series_overall_design"] = value
                    elif "pubmed_id" in key:
                        metadata["series_pubmed_id"] = value
                    elif "contact_name" in key:
                        metadata["series_contact_name"] = value
                    elif "contact_email" in key:
                        metadata["series_contact_email"] = value

            # Parse sample metadata
            elif line.startswith('!Sample_'):
                parts = line.split('=', 1)
                if len(parts) == 2:
                    field = parts[0].replace('!Sample_', '').strip()
                    value = parts[1].strip().strip('"')
                    # Extract sample ID from first sample_title or similar
                    if "title" in field.lower() and value and value not in sample_metadata:
                        samples.append(value)
                        sample_metadata[value] = {}
                    # Store all sample attributes
                    if samples:
                        if field not in sample_metadata[samples[-1]]:
                            sample_metadata[samples[-1]][field] = []
                        sample_metadata[samples[-1]][field].append(value)

            # Mark matrix section
            elif "!series_matrix_table_begin" in line:
                in_matrix = True
                continue
            elif "!series_matrix_table_end" in line:
                in_matrix = False
                continue

            # Collect matrix lines
            if in_matrix and line and not line.startswith('!'):
                expression_lines.append(line)

    # Parse expression matrix
    if expression_lines:
        matrix_text = '\n'.join(expression_lines)
        df = pd.read_csv(StringIO(matrix_text), sep='\t', index_col=0)
    else:
        df = pd.DataFrame()

    return df, metadata, samples


def extract_sample_metadata(sample_metadata_dict: Dict, samples: List) -> pd.DataFrame:
    """
    Convert sample metadata dict to DataFrame.

    Args:
        sample_metadata_dict: dict mapping sample name to attributes
        samples: ordered list of samples

    Returns:
        DataFrame with samples as rows, attributes as columns
    """
    rows = []
    for sample in samples:
        if sample in sample_metadata_dict:
            row = {"sample_id": sample}
            for attr, values in sample_metadata_dict[sample].items():
                # Join multiple values if present
                row[attr] = "; ".join(values) if isinstance(values, list) else values
            rows.append(row)
        else:
            rows.append({"sample_id": sample})

    if rows:
        return pd.DataFrame(rows)
    else:
        return pd.DataFrame({"sample_id": samples})


def main():
    parser = argparse.ArgumentParser(
        description="Download a GEO dataset by accession"
    )
    parser.add_argument(
        "--accession",
        type=str,
        required=True,
        help="GEO accession (GSE12345 or GDS100)"
    )
    parser.add_argument(
        "--outdir",
        type=str,
        required=True,
        help="Output directory"
    )
    parser.add_argument(
        "--soft",
        action="store_true",
        help="Also download SOFT file"
    )
    parser.add_argument(
        "--matrix",
        action="store_true",
        default=True,
        help="Download series matrix (default: true)"
    )
    parser.add_argument(
        "--manifest",
        type=str,
        default="",
        help="Optional path to write manifest JSON"
    )

    args = parser.parse_args()

    accession = args.accession.upper().strip()
    acc_type = get_accession_type(accession)

    if not acc_type:
        print(f"Error: Invalid accession format: {accession}")
        print("Use GSE* (series) or GDS* (dataset)")
        return 1

    if acc_type not in ["GSE", "GDS"]:
        print(f"Error: Unsupported accession type: {acc_type}")
        print("Only GSE and GDS are supported")
        return 1

    ensure_dir(args.outdir)

    print(f"\nGEO Download: {accession}")
    print(f"Type: {acc_type}")
    print(f"Output: {args.outdir}\n")

    records: List[FileRecord] = []

    # Get series metadata
    print(f"Querying NCBI for series information...")
    series_info = get_series_summary(accession)
    if series_info.get("title"):
        print(f"  Title: {series_info['title']}")
    if series_info.get("organism"):
        print(f"  Organism: {series_info['organism']}")

    # Download series matrix
    if acc_type == "GSE":
        print(f"\nDownloading expression matrix...")
        ftp_url = construct_ftp_url(accession)
        print(f"  URL: {ftp_url}")

        filename = f"{accession}_series_matrix.txt.gz"
        output_path = os.path.join(args.outdir, filename)

        if os.path.exists(output_path):
            print(f"  {filename} (exists, skipped)")
            size = os.path.getsize(output_path)
            records.append(FileRecord(
                accession=accession,
                record_type="series_matrix",
                filename=filename,
                size_bytes=size,
                status="exists",
                local_path=output_path,
            ))
        else:
            success, error, size = download_file(ftp_url, output_path)

            if success:
                print(f"  Downloaded {filename} ({size // 1024} KB)")
                records.append(FileRecord(
                    accession=accession,
                    record_type="series_matrix",
                    filename=filename,
                    size_bytes=size,
                    status="downloaded",
                    local_path=output_path,
                ))

                # Parse matrix
                try:
                    print(f"  Parsing expression matrix...")
                    df_expr, metadata, samples = parse_series_matrix(output_path)

                    if not df_expr.empty:
                        print(f"    Expression: {df_expr.shape[0]} probes × {df_expr.shape[1]} samples")

                        # Save expression matrix
                        expr_path = os.path.join(args.outdir, "expression_matrix.tsv")
                        df_expr.to_csv(expr_path, sep='\t')
                        print(f"    Saved: expression_matrix.tsv")

                        # Save sample metadata
                        if samples:
                            df_samples = extract_sample_metadata(
                                {s: {"sample": s} for s in samples}, samples
                            )
                            sample_path = os.path.join(args.outdir, "sample_metadata.tsv")
                            df_samples.to_csv(sample_path, sep='\t', index=False)
                            print(f"    Saved: sample_metadata.tsv ({len(samples)} samples)")

                        # Save series info
                        series_info.update(metadata)
                        series_info["n_samples"] = len(samples)
                        series_path = os.path.join(args.outdir, "series_info.json")
                        with open(series_path, 'w') as f:
                            json.dump(series_info, f, indent=2)
                        print(f"    Saved: series_info.json")

                    else:
                        print(f"    Warning: No expression matrix found in file")

                except Exception as e:
                    print(f"    Error parsing: {e}")

            else:
                print(f"  Download failed: {error}")
                records.append(FileRecord(
                    accession=accession,
                    record_type="series_matrix",
                    filename=filename,
                    size_bytes=0,
                    status="failed",
                    error_msg=error,
                ))

    # Download SOFT file if requested
    if args.soft:
        print(f"\nDownloading SOFT file...")
        soft_filename = f"{accession}_family.soft.gz"
        soft_url = f"{GEO_FTP_BASE}/series/{accession[:6]}nnn/{accession}/{soft_filename}"
        soft_path = os.path.join(args.outdir, soft_filename)

        if os.path.exists(soft_path):
            print(f"  {soft_filename} (exists, skipped)")
            size = os.path.getsize(soft_path)
            records.append(FileRecord(
                accession=accession,
                record_type="soft",
                filename=soft_filename,
                size_bytes=size,
                status="exists",
                local_path=soft_path,
            ))
        else:
            success, error, size = download_file(soft_url, soft_path)
            if success:
                print(f"  Downloaded {soft_filename} ({size // (1024*1024)} MB)")
                records.append(FileRecord(
                    accession=accession,
                    record_type="soft",
                    filename=soft_filename,
                    size_bytes=size,
                    status="downloaded",
                    local_path=soft_path,
                ))
            else:
                print(f"  Download failed: {error}")
                records.append(FileRecord(
                    accession=accession,
                    record_type="soft",
                    filename=soft_filename,
                    size_bytes=0,
                    status="failed",
                    error_msg=error,
                ))

    # Write manifest
    if args.manifest:
        manifest = {
            "accession": accession,
            "type": acc_type,
            "download_timestamp": pd.Timestamp.now().isoformat(),
            "outdir": os.path.abspath(args.outdir),
            "title": series_info.get("title", ""),
            "organism": series_info.get("organism", ""),
            "files": [asdict(r) for r in records],
            "summary": {
                "successful": len([r for r in records if r.status in ["downloaded", "exists"]]),
                "failed": len([r for r in records if r.status == "failed"]),
            }
        }
        with open(args.manifest, 'w') as f:
            json.dump(manifest, f, indent=2)
        print(f"\nManifest written to {args.manifest}")

    # Summary
    print(f"\nSummary:")
    print(f"  Files processed: {len(records)}")
    successful = len([r for r in records if r.status in ["downloaded", "exists"]])
    print(f"  Successful: {successful}")
    failed = len([r for r in records if r.status == "failed"])
    if failed > 0:
        print(f"  Failed: {failed}")

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
