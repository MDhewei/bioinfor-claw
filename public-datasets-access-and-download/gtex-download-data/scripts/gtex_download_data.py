#!/usr/bin/env python3
"""
GTEx Download Data - Download gene expression data from GTEx Portal.

Downloads tissue-gene expression matrices from GTEx v8 or v10,
with optional filtering by tissue and gene.
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
from io import TextIOWrapper

import requests
import pandas as pd
import numpy as np


GTEX_API_BASE = "https://gtexportal.org/api/v2"
GTEX_PORTAL = "https://gtexportal.org"

# Static file URLs for GTEx releases
GTEX_STATIC_URLS = {
    "v8": {
        "median_tpm": "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
        "read_counts": "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz",
    },
    "v10": {
        "median_tpm": "https://storage.googleapis.com/gtex_analysis_v10_releases/single_tissue_analysis_data/expression_matrices/GTEx_Tissue_Sample_Attributes.txt",
    }
}

SESSION = requests.Session()
SESSION.headers.update({
    "User-Agent": "gtex-download/1.0",
    "Accept": "application/json",
})

RETRY_ATTEMPTS = 3
RETRY_DELAY = 2  # seconds

# Sample GTEx tissue codes and names
GTEX_TISSUES = {
    "Adrenal Gland": "ADRNL",
    "Adipose - Subcutaneous": "ADPSBQ",
    "Adipose - Visceral (Omentum)": "ADPVSC",
    "Artery - Aorta": "ARTAOR",
    "Artery - Coronary": "ARTCOR",
    "Artery - Tibial": "ARTTIB",
    "Bladder": "BLDBDR",
    "Brain - Amygdala": "BRNAMG",
    "Brain - Anterior cingulate cortex (BA24)": "BRNACG",
    "Brain - Caudate (basal ganglia)": "BRNCDT",
    "Brain - Cerebellar Hemisphere": "BRNCRB",
    "Brain - Cerebellum": "BRNCBM",
    "Brain - Cortex": "BRNCRT",
    "Brain - Frontal Cortex (BA9)": "BRNFBA",
    "Brain - Hippocampus": "BRNHIP",
    "Brain - Hypothalamus": "BRNHYP",
    "Brain - Nucleus accumbens (basal ganglia)": "BRNNAC",
    "Brain - Putamen (basal ganglia)": "BRNPTM",
    "Brain - Spinal cord (cervical c-1)": "BRNSPC",
    "Brain - Substantia nigra": "BRNSNC",
    "Breast - Mammary Tissue": "BRSTMM",
    "Cells - EBV-transformed lymphocytes": "CLSEBV",
    "Cells - Leukemia cell line (CML)": "CLSLCL",
    "Cervix - Ectocervix": "CERVEC",
    "Cervix - Endocervix": "CERVIC",
    "Colon - Sigmoid": "COLSIG",
    "Colon - Transverse": "COLTRS",
    "Esophagus - Gastroesophageal Junction": "ESGJCN",
    "Esophagus - Muscularis": "ESPMSC",
    "Esophagus - Mucosa": "ESPMCO",
    "Fallopian Tube": "FLLPBN",
    "Heart - Atrial Appendage": "HRTAA",
    "Heart - Left Ventricle": "HRTLV",
    "Kidney - Cortex": "KDNCTX",
    "Kidney - Medulla": "KDNMDL",
    "Liver": "LIVER",
    "Lung": "LUNG",
    "Muscle - Skeletal": "MSCLSK",
    "Nerve - Tibial": "NERV",
    "Ovary": "OVARY",
    "Pancreas": "PNCRS",
    "Pituitary": "PTRY",
    "Prostate": "PROST",
    "Salivary Gland": "SLVRY",
    "Skin - Not Sun Exposed (Suprapubic)": "SKNSUP",
    "Skin - Sun Exposed (Lower leg)": "SKNEXP",
    "Small Intestine - Terminal Ileum": "SNTIL",
    "Spleen": "SPLEN",
    "Stomach": "STMCH",
    "Testis": "TESTIS",
    "Thyroid": "THYRD",
    "Uterus": "UTERUS",
    "Vagina": "VAGINA",
    "Whole Blood": "WHLBLOD",
}


@dataclass
class FileRecord:
    data_type: str
    version: str
    filename: str
    size_bytes: int
    status: str = "pending"
    local_path: str = ""
    error_msg: str = ""
    matrix_shape: str = ""


def ensure_dir(path: str) -> None:
    """Create directory if it does not exist."""
    os.makedirs(path, exist_ok=True)


def query_gtex_api(endpoint: str, params: Dict = None) -> Dict:
    """
    Query GTEx Portal API.

    Args:
        endpoint: API endpoint path (e.g., '/expression/medianTranscriptExpression')
        params: query parameters

    Returns:
        JSON response dictionary
    """
    if params is None:
        params = {}

    url = f"{GTEX_API_BASE}{endpoint}"

    for attempt in range(RETRY_ATTEMPTS):
        try:
            response = SESSION.get(url, params=params, timeout=30)
            response.raise_for_status()
            return response.json()
        except requests.RequestException as e:
            if attempt < RETRY_ATTEMPTS - 1:
                print(f"  [Retry {attempt + 1}/{RETRY_ATTEMPTS}] API request failed: {e}")
                time.sleep(RETRY_DELAY)
            else:
                raise


def get_tissue_sites() -> Dict[str, str]:
    """
    Fetch GTEx tissue sites from API.

    Returns:
        dict mapping tissue name to tissue code
    """
    try:
        response = query_gtex_api("/reference/tissueSites")
        tissue_map = {}
        for tissue in response.get("data", []):
            tissue_map[tissue.get("tissueSiteDetail")] = tissue.get("tissueSiteDetailId")
        return tissue_map
    except Exception as e:
        print(f"Warning: Could not fetch tissue sites from API: {e}")
        return GTEX_TISSUES


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

        except Exception as e:
            if attempt < RETRY_ATTEMPTS - 1:
                print(f"      [Retry {attempt + 1}/{RETRY_ATTEMPTS}] Download failed: {e}")
                time.sleep(RETRY_DELAY)
            else:
                return False, str(e), 0

    return False, "Max retries exceeded", 0


def parse_gct_file(filepath: str) -> pd.DataFrame:
    """
    Parse GTEx GCT format file.

    GCT format:
    - Line 1: version (e.g., #1.2)
    - Line 2: dimensions (genes, samples)
    - Line 3: header (Name, Description, sample1, sample2, ...)
    - Lines 4+: gene data

    Args:
        filepath: path to .gct or .gct.gz file

    Returns:
        pandas DataFrame with genes as index
    """
    open_func = gzip.open if filepath.endswith('.gz') else open

    with open_func(filepath, 'rt') as f:
        # Read version and dimensions
        version = f.readline().strip()
        dims = f.readline().strip().split('\t')
        n_genes, n_samples = int(dims[0]), int(dims[1])

        # Read data
        df = pd.read_csv(f, sep='\t', index_col=0)
        # Remove Description column if present
        if 'Description' in df.columns:
            df = df.drop('Description', axis=1)

    return df


def download_median_tpm(version: str, outdir: str) -> Optional[FileRecord]:
    """
    Download median TPM expression matrix.

    Args:
        version: GTEx version (v8 or v10)
        outdir: output directory

    Returns:
        FileRecord with download status
    """
    ensure_dir(outdir)

    if version not in GTEX_STATIC_URLS:
        return FileRecord(
            data_type="median_tpm",
            version=version,
            filename="",
            size_bytes=0,
            status="failed",
            error_msg=f"Unsupported version: {version}",
        )

    url = GTEX_STATIC_URLS[version].get("median_tpm")
    if not url:
        return FileRecord(
            data_type="median_tpm",
            version=version,
            filename="",
            size_bytes=0,
            status="failed",
            error_msg=f"No median_tpm available for version {version}",
        )

    filename = f"gtex_{version}_median_tpm.gct.gz"
    output_path = os.path.join(outdir, filename)

    if os.path.exists(output_path):
        print(f"  {filename} (exists, skipped)")
        size = os.path.getsize(output_path)
        return FileRecord(
            data_type="median_tpm",
            version=version,
            filename=filename,
            size_bytes=size,
            status="exists",
            local_path=output_path,
        )

    print(f"  Downloading {filename}...")
    success, error, size = download_file(url, output_path)

    if success:
        print(f"  Downloaded {filename} ({size // (1024*1024)} MB)")
        # Parse and get dimensions
        try:
            df = parse_gct_file(output_path)
            record = FileRecord(
                data_type="median_tpm",
                version=version,
                filename=filename,
                size_bytes=size,
                status="downloaded",
                local_path=output_path,
                matrix_shape=f"{df.shape[0]} genes x {df.shape[1]} tissues",
            )
            return record
        except Exception as e:
            print(f"  Error parsing: {e}")
            return FileRecord(
                data_type="median_tpm",
                version=version,
                filename=filename,
                size_bytes=size,
                status="failed",
                local_path=output_path,
                error_msg=str(e),
            )
    else:
        print(f"  Download failed: {error}")
        return FileRecord(
            data_type="median_tpm",
            version=version,
            filename=filename,
            size_bytes=0,
            status="failed",
            error_msg=error,
        )


def filter_and_save_matrix(df: pd.DataFrame, tissues: List[str], genes: List[str],
                           output_path: str, file_format: str = "tsv") -> pd.DataFrame:
    """
    Filter matrix by tissues and genes, then save.

    Args:
        df: input dataframe
        tissues: list of tissue names to keep (None = all)
        genes: list of gene symbols to keep (None = all)
        output_path: path to save filtered matrix
        file_format: tsv or csv

    Returns:
        filtered dataframe
    """
    result = df.copy()

    # Filter tissues (columns)
    if tissues:
        valid_cols = [c for c in result.columns if any(t.lower() in c.lower() for t in tissues)]
        result = result[valid_cols]

    # Filter genes (rows)
    if genes:
        valid_rows = [r for r in result.index if any(g.lower() in r.lower() for g in genes)]
        result = result.loc[valid_rows]

    # Save
    sep = '\t' if file_format.lower() == 'tsv' else ','
    result.to_csv(output_path, sep=sep)

    return result


def query_gene_expression(genes: List[str], version: str = "v8") -> Optional[pd.DataFrame]:
    """
    Query specific gene expression via GTEx API.

    Args:
        genes: list of gene symbols or ENSG IDs
        version: GTEx version

    Returns:
        dataframe with genes as rows, tissues as columns
    """
    try:
        # Get tissue sites
        tissues = get_tissue_sites()

        # Query each gene
        all_data = []
        for gene in genes[:50]:  # Limit to first 50 to avoid too many API calls
            print(f"    Querying {gene}...")
            try:
                response = query_gtex_api("/expression/geneExpression", {
                    "geneId": gene,
                    "datasetId": "gtex_v8_sr" if version == "v8" else "gtex_v10_sr",
                })

                for tissue_expr in response.get("data", []):
                    all_data.append({
                        "gene": gene,
                        "tissue": tissue_expr.get("tissueSiteDetailId"),
                        "median_tpm": tissue_expr.get("median", 0),
                    })
            except Exception as e:
                print(f"      Error querying {gene}: {e}")
                continue

        if not all_data:
            return None

        # Pivot to genes x tissues
        df = pd.DataFrame(all_data)
        df_pivot = df.pivot_table(index="gene", columns="tissue", values="median_tpm", fill_value=0)
        return df_pivot

    except Exception as e:
        print(f"  Error querying gene expression: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(
        description="Download GTEx gene expression data"
    )
    parser.add_argument(
        "--data-type",
        choices=["median_tpm", "sample_tpm", "read_counts", "tissue_metadata", "gene_tpm"],
        required=True,
        help="Type of expression data to download"
    )
    parser.add_argument(
        "--genes",
        type=str,
        default="",
        help="Comma-separated gene symbols (for gene_tpm)"
    )
    parser.add_argument(
        "--tissues",
        type=str,
        default="",
        help="Comma-separated tissue names to filter"
    )
    parser.add_argument(
        "--version",
        choices=["v8", "v10"],
        default="v8",
        help="GTEx version (default: v8)"
    )
    parser.add_argument(
        "--outdir",
        type=str,
        required=True,
        help="Output directory"
    )
    parser.add_argument(
        "--list-tissues",
        action="store_true",
        help="List available tissues and exit"
    )
    parser.add_argument(
        "--manifest",
        type=str,
        default="",
        help="Optional path to write manifest JSON"
    )

    args = parser.parse_args()

    # Handle list tissues
    if args.list_tissues:
        print("Available GTEx tissues:")
        tissue_map = get_tissue_sites()
        for tissue in sorted(tissue_map.keys()):
            print(f"  {tissue}")
        return 0

    # Parse genes and tissues
    genes = [g.strip() for g in args.genes.split(",") if g.strip()] if args.genes else None
    tissues = [t.strip() for t in args.tissues.split(",") if t.strip()] if args.tissues else None

    ensure_dir(args.outdir)

    print(f"\nGTEx Download: {args.data_type}")
    print(f"Version: {args.version}")
    if tissues:
        print(f"Tissues: {', '.join(tissues[:3])}{'...' if len(tissues) > 3 else ''}")
    if genes:
        print(f"Genes: {', '.join(genes[:3])}{'...' if len(genes) > 3 else ''}")
    print(f"Output: {args.outdir}\n")

    records: List[FileRecord] = []

    # Download based on data type
    if args.data_type == "median_tpm":
        record = download_median_tpm(args.version, args.outdir)
        records.append(record)

        # Parse and filter if needed
        if record.status in ["downloaded", "exists"]:
            try:
                print(f"  Parsing {record.filename}...")
                df = parse_gct_file(record.local_path)
                print(f"    Shape: {df.shape[0]} genes × {df.shape[1]} tissues")

                # Filter and save
                if tissues or genes:
                    df_filtered = filter_and_save_matrix(
                        df, tissues, genes,
                        os.path.join(args.outdir, "expression_matrix_filtered.tsv"),
                        "tsv"
                    )
                    print(f"    Filtered: {df_filtered.shape[0]} genes × {df_filtered.shape[1]} tissues")
                else:
                    # Save full matrix
                    df.to_csv(os.path.join(args.outdir, "expression_matrix.tsv"), sep='\t')

            except Exception as e:
                print(f"  Error processing: {e}")
                record.status = "failed"
                record.error_msg = str(e)

    elif args.data_type == "gene_tpm":
        if not genes:
            print("Error: --genes required for gene_tpm data type")
            return 1

        df = query_gene_expression(genes, args.version)
        if df is not None:
            output_path = os.path.join(args.outdir, "gene_expression.tsv")
            df.to_csv(output_path, sep='\t')
            print(f"  Saved gene expression: {df.shape[0]} genes × {df.shape[1]} tissues")
            records.append(FileRecord(
                data_type="gene_tpm",
                version=args.version,
                filename="gene_expression.tsv",
                size_bytes=os.path.getsize(output_path),
                status="downloaded",
                local_path=output_path,
                matrix_shape=f"{df.shape[0]} genes × {df.shape[1]} tissues",
            ))
        else:
            print("  Error downloading gene expression data")
            return 1

    elif args.data_type == "tissue_metadata":
        try:
            tissues_data = get_tissue_sites()
            df = pd.DataFrame([
                {"tissue": k, "code": v} for k, v in tissues_data.items()
            ])
            output_path = os.path.join(args.outdir, "tissue_metadata.tsv")
            df.to_csv(output_path, sep='\t', index=False)
            print(f"  Saved tissue metadata: {len(df)} tissues")
            records.append(FileRecord(
                data_type="tissue_metadata",
                version=args.version,
                filename="tissue_metadata.tsv",
                size_bytes=os.path.getsize(output_path),
                status="downloaded",
                local_path=output_path,
            ))
        except Exception as e:
            print(f"  Error: {e}")
            return 1

    else:
        print(f"Data type '{args.data_type}' not yet fully implemented")
        return 1

    # Write manifest
    if args.manifest and records:
        manifest = {
            "data_type": args.data_type,
            "version": args.version,
            "download_timestamp": pd.Timestamp.now().isoformat(),
            "outdir": os.path.abspath(args.outdir),
            "genes_requested": genes,
            "tissues_requested": tissues,
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
