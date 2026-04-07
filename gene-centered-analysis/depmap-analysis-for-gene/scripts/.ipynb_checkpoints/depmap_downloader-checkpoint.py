#!/usr/bin/env python3

import argparse
import json
import os
import re
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional
from urllib.parse import urljoin

import requests


DEP_MAP_DOWNLOAD_CATALOG = "https://depmap.org/portal/download/api/downloads"
DEP_MAP_BASE = "https://depmap.org"

SESSION = requests.Session()
SESSION.headers.update(
    {
        "User-Agent": "depmap-download/0.1",
        "Accept": "application/json",
    }
)


@dataclass
class DownloadItem:
    dataset_key: str
    release_name: str
    file_name: str
    download_url: str
    local_path: str
    status: str
    note: str = ""


DEFAULT_FILE_PATTERNS = {
    "expression": [
        r"^OmicsExpressionTPMLogp1HumanProteinCodingGenes\.csv$",
        r"^OmicsExpressionProteinCodingGenesTPMLogp1.*\.csv$",
        r"^OmicsExpression.*ProteinCoding.*TPM.*logp1.*\.csv$",
    ],
    "mutations": [
        r"^OmicsSomaticMutations(?:Profile)?\.csv$",
        r"^OmicsSomaticMutations.*\.csv$",
        r"^OmicsMutations.*\.csv$",
    ],
    "copy_number": [
        r"^OmicsCNGene\.csv$",
        r"^OmicsCNGene.*\.csv$",
        r"^OmicsCopyNumber.*Gene.*\.csv$",
    ],
    "essentiality": [
        r"^CRISPRGeneEffect\.csv$",
        r"^CRISPRGeneEffect.*\.csv$",
        r"^CRISPR.*GeneEffect.*\.csv$",
        r"^Chronos.*GeneEffect.*\.csv$",
    ],
    "metadata": [
        r"^Model\.csv$",
        r"^Model.*\.csv$",
        r"^.*Model.*\.csv$",
    ],
}


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def safe_get_json(url: str, timeout: int = 120) -> Dict:
    r = SESSION.get(url, timeout=timeout)
    r.raise_for_status()
    return r.json()


def normalize_download_url(download_url: str) -> str:
    if not download_url:
        raise ValueError("Empty download URL")
    if download_url.startswith("http://") or download_url.startswith("https://"):
        return download_url
    return urljoin(DEP_MAP_BASE, download_url)


def fetch_download_catalog() -> Dict:
    return safe_get_json(DEP_MAP_DOWNLOAD_CATALOG)


def get_release_names(catalog: Dict) -> List[str]:
    releases = catalog.get("releaseData", []) or []
    names = [x.get("releaseName") for x in releases if x.get("releaseName")]
    return names


def choose_release_name(catalog: Dict, requested_release: Optional[str]) -> str:
    releases = catalog.get("releaseData", []) or []

    if requested_release:
        for rel in releases:
            if rel.get("releaseName") == requested_release:
                return requested_release
        available = ", ".join(get_release_names(catalog))
        raise ValueError(
            f"Requested release '{requested_release}' not found. Available releases include: {available}"
        )

    for rel in releases:
        if rel.get("isLatest") is True:
            return rel["releaseName"]

    names = get_release_names(catalog)
    if not names:
        raise RuntimeError("No release names found in DepMap download catalog.")
    return names[0]


def collect_all_file_records(catalog: Dict, release_name: str) -> List[Dict]:
    records: List[Dict] = []

    def walk(obj):
        if isinstance(obj, dict):
            if obj.get("releaseName") == release_name and obj.get("fileName") and obj.get("downloadUrl"):
                records.append(obj)
            for v in obj.values():
                walk(v)
        elif isinstance(obj, list):
            for x in obj:
                walk(x)

    walk(catalog)

    seen = set()
    unique_records = []
    for r in records:
        key = (r.get("releaseName"), r.get("fileName"), r.get("downloadUrl"))
        if key not in seen:
            seen.add(key)
            unique_records.append(r)

    return unique_records


def score_filename_match(file_name: str, dataset_key: str) -> int:
    patterns = DEFAULT_FILE_PATTERNS[dataset_key]
    for i, pat in enumerate(patterns):
        if re.search(pat, file_name, flags=re.IGNORECASE):
            return 100 - i

    # Soft fallback keyword scoring
    fn = file_name.lower()
    score = 0
    if dataset_key == "expression":
        if "expression" in fn:
            score += 30
        if "tpm" in fn:
            score += 20
        if "proteincoding" in fn or "protein_coding" in fn:
            score += 20
        if "logp1" in fn:
            score += 15
    elif dataset_key == "mutations":
        if "mutation" in fn:
            score += 40
        if "somatic" in fn:
            score += 20
    elif dataset_key == "copy_number":
        if "cn" in fn or "copy" in fn:
            score += 30
        if "gene" in fn:
            score += 20
    elif dataset_key == "essentiality":
        if "geneeffect" in fn:
            score += 40
        if "crispr" in fn or "chronos" in fn:
            score += 20
    elif dataset_key == "metadata":
        if fn == "model.csv":
            score += 100
        elif "model" in fn:
            score += 30

    return score


def find_best_file_record(records: List[Dict], dataset_key: str) -> Optional[Dict]:
    scored = []
    for r in records:
        file_name = r.get("fileName", "")
        score = score_filename_match(file_name, dataset_key)
        if score > 0:
            scored.append((score, r))

    if not scored:
        return None

    scored.sort(key=lambda x: x[0], reverse=True)
    return scored[0][1]


def stream_download(url: str, out_path: str, chunk_size: int = 1024 * 1024) -> None:
    ensure_dir(os.path.dirname(os.path.abspath(out_path)))
    with SESSION.get(url, stream=True, timeout=600) as r:
        r.raise_for_status()
        with open(out_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)


def build_output_filename(dataset_key: str, file_name: str) -> str:
    # Keep original filename when possible
    return file_name


def download_dataset(
    records: List[Dict],
    dataset_key: str,
    release_name: str,
    outdir: str,
    overwrite: bool = False,
) -> DownloadItem:
    record = find_best_file_record(records, dataset_key)
    if record is None:
        available = sorted({r.get("fileName", "") for r in records})
        raise RuntimeError(
            f"Could not find a matching file for dataset '{dataset_key}' in release '{release_name}'. "
            f"Available files include: {available[:20]}"
        )

    file_name = record["fileName"]
    download_url = normalize_download_url(record["downloadUrl"])
    local_path = os.path.join(outdir, build_output_filename(dataset_key, file_name))

    if os.path.exists(local_path) and not overwrite:
        return DownloadItem(
            dataset_key=dataset_key,
            release_name=release_name,
            file_name=file_name,
            download_url=download_url,
            local_path=local_path,
            status="exists",
            note="File already exists; skipped download",
        )

    stream_download(download_url, local_path)

    return DownloadItem(
        dataset_key=dataset_key,
        release_name=release_name,
        file_name=file_name,
        download_url=download_url,
        local_path=local_path,
        status="downloaded",
        note="Downloaded successfully",
    )


def main():
    parser = argparse.ArgumentParser(
        description="Download DepMap release files using the official DepMap downloads API catalog."
    )
    parser.add_argument(
        "--release",
        default=None,
        help='Release name, e.g. "DepMap Public 26Q1". If omitted, use latest release.',
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Directory to save downloaded files",
    )
    parser.add_argument("--expression", action="store_true", help="Download expression dataset")
    parser.add_argument("--mutations", action="store_true", help="Download mutations dataset")
    parser.add_argument("--copy-number", action="store_true", help="Download copy number dataset")
    parser.add_argument("--essentiality", action="store_true", help="Download CRISPR gene effect / essentiality dataset")
    parser.add_argument("--metadata", action="store_true", help="Download model metadata dataset")
    parser.add_argument("--all", action="store_true", help="Download all recommended datasets")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing local files")
    parser.add_argument("--manifest", default=None, help="Optional JSON manifest output path")
    parser.add_argument("--list-releases", action="store_true", help="List releases and exit")
    parser.add_argument("--list-files", action="store_true", help="List files in the chosen release and exit")
    args = parser.parse_args()

    ensure_dir(args.outdir)

    print("[INFO] Fetching DepMap download catalog ...")
    catalog = fetch_download_catalog()

    if args.list_releases:
        print("[INFO] Available releases:")
        for name in get_release_names(catalog):
            print(name)
        return

    release_name = choose_release_name(catalog, args.release)
    print(f"[INFO] Using release: {release_name}")

    records = collect_all_file_records(catalog, release_name)
    if not records:
        raise RuntimeError(f"No downloadable file records found for release '{release_name}'.")

    if args.list_files:
        print(f"[INFO] Files in release '{release_name}':")
        for name in sorted({r.get('fileName', '') for r in records}):
            print(name)
        return

    if args.all:
        wanted_keys = ["expression", "mutations", "copy_number", "essentiality", "metadata"]
    else:
        wanted_keys = []
        if args.expression:
            wanted_keys.append("expression")
        if args.mutations:
            wanted_keys.append("mutations")
        if args.copy_number:
            wanted_keys.append("copy_number")
        if args.essentiality:
            wanted_keys.append("essentiality")
        if args.metadata:
            wanted_keys.append("metadata")

    if not wanted_keys:
        raise ValueError(
            "No datasets selected. Use --all or one or more of "
            "--expression --mutations --copy-number --essentiality --metadata"
        )

    manifest: List[DownloadItem] = []

    for dataset_key in wanted_keys:
        print(f"[INFO] Downloading dataset: {dataset_key}")
        item = download_dataset(
            records=records,
            dataset_key=dataset_key,
            release_name=release_name,
            outdir=args.outdir,
            overwrite=args.overwrite,
        )
        print(f"[DONE] {dataset_key}: {item.local_path} ({item.status})")
        manifest.append(item)

    if args.manifest:
        manifest_path = args.manifest
        ensure_dir(os.path.dirname(os.path.abspath(manifest_path)))
        with open(manifest_path, "w", encoding="utf-8") as f:
            json.dump([asdict(x) for x in manifest], f, indent=2, ensure_ascii=False)
        print(f"[DONE] Manifest written to: {manifest_path}")


if __name__ == "__main__":
    main()