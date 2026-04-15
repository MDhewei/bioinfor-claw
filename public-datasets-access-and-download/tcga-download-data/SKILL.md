---
name: tcga-download-data
description: Download TCGA genomic data (gene expression, somatic mutations, CNV, clinical) from the GDC Data Portal API for one or more cancer types. Returns per-file manifests and merged TSV matrices.
---

# TCGA Download Data

## Purpose

Download cancer genomic datasets from the National Cancer Institute's GDC Data Portal.

This skill supports:
- querying available TCGA cancer types
- filtering by genomic data type (expression, mutations, CNV, clinical, methylation)
- downloading raw data files in batch
- merging related data files into matrices
- generating file manifests with checksums
- handling large datasets with streaming downloads

This skill is for **data acquisition only**. It does not perform biological analysis.

## Supported data types

Logical dataset types available from GDC:

- `expression` – Gene Expression Quantification (STAR RNA-seq counts, TPM)
- `mutations` – Masked Somatic Mutations (MAF format)
- `cnv` – Copy Number Variation (gene-level CN calls)
- `clinical` – Clinical metadata (demographics, diagnosis, treatment)
- `methylation` – DNA methylation (HM27 or HM450 beta values)

## When to use

Use this skill when you need to:
- download gene expression profiles from TCGA for specific cancers
- obtain somatic mutation data for mutational analysis
- retrieve copy number variation calls
- fetch clinical annotations and survival data
- prepare TCGA data for downstream genomic analysis

This skill may also be called internally by analysis or integration skills.

## When not to use

Do not use this skill when you want:
- differential expression analysis results
- mutation burden summaries or interpretation
- copy number analysis or visualization
- clinical survival prediction
- integrated multi-omics interpretation

Use a TCGA analysis skill for those tasks.

## Data source

Primary data source:
- GDC Data Portal REST API: `https://api.gdc.cancer.gov/`

The API provides:
- no authentication required for public TCGA data
- structured querying by project, data type, file format
- direct file download with MD5 verification
- project metadata and sample-level annotations

## Parameters

### `--data-type`
Type of genomic data to download.

Choices:
- `expression` – gene expression quantification
- `mutations` – somatic mutations
- `cnv` – copy number variation
- `clinical` – clinical metadata
- `methylation` – DNA methylation

Required: yes
Default: none (must be specified)

### `--cancer-types`
One or more TCGA cancer types to download.

Format:
- comma-separated TCGA project codes: `TCGA-BRCA,TCGA-LUAD`
- or special value `all` to download all available projects

Examples:
- `TCGA-BRCA` – breast cancer
- `TCGA-LUAD` – lung adenocarcinoma
- `TCGA-OV` – ovarian cancer
- `all` – all TCGA cancer types

If omitted:
- an interactive list of available cancer types will be shown

### `--outdir`
Output directory for downloaded and merged files.

Use when:
- files need to be cached for reuse
- downstream tools need predictable file locations

Example:
- `data/tcga_expression`

Recommended:
- use separate directories per data type or cancer type
- the script creates subdirectories: `raw/`, `merged/`

### `--max-files`
Maximum number of files to download per cancer type.

Use when:
- testing with a small subset first
- limiting download size for slow networks

Default: 5 (suitable for testing)

Example:
- `--max-files 50` – download up to 50 files per project

### `--file-format`
Output format for merged matrices.

Choices:
- `TSV` – tab-separated values (default)
- `CSV` – comma-separated values

Note:
- individual downloaded files retain their original format
- merged matrices use the specified format

### `--manifest`
Output path for manifest JSON describing downloaded files.

Use when:
- downstream analysis needs file locations and checksums
- reproducibility requires tracking data provenance
- multiple files are downloaded

Example:
- `data/tcga_expression/manifest.json`

Manifest includes:
- file ID, name, size, MD5 checksum
- download status
- merge status
- timestamp

## Parameter decision guide

| Goal | Data type | Cancer types | Notes |
|------|-----------|--------------|-------|
| Gene expression for 1 cancer | `expression` | `TCGA-BRCA` | merged into genes × samples matrix |
| Mutations in multiple cancers | `mutations` | `TCGA-BRCA,TCGA-LUAD` | MAF files merged by mutation locus |
| Copy number across all TCGA | `cnv` | `all` | large download; consider `--max-files` |
| Clinical trial data | `clinical` | `TCGA-BRCA` | includes survival and treatment info |
| Test run (small sample) | `expression` | `TCGA-BRCA` | add `--max-files 2` |

## Procedure

1. **Query GDC API** for available TCGA projects
2. **Parse cancer type codes** from command line or interactive list
3. **Filter files** by data_category and data_type
4. **Stream download** each file with progress tracking
5. **Save to `outdir/raw/`** preserving original filename and format
6. **Merge related files**:
   - expression: concatenate gene counts across samples → genes × samples TSV
   - mutations: combine MAF files by chromosome/locus
   - CNV: stack gene-level calls across samples
7. **Write manifest** JSON with file metadata and checksums
8. **Report**: files found, downloaded, merged, any errors

## Key execution patterns

### Download gene expression for breast cancer
```bash
python scripts/tcga_download_data.py \
  --data-type expression \
  --cancer-types TCGA-BRCA \
  --outdir data/tcga_brca_expr \
  --max-files 10
```

### Download mutations for multiple cancers, merge to single MAF
```bash
python scripts/tcga_download_data.py \
  --data-type mutations \
  --cancer-types TCGA-BRCA,TCGA-LUAD,TCGA-OV \
  --outdir data/tcga_mutations \
  --max-files 20 \
  --manifest data/tcga_mutations/manifest.json
```

### Test run with small subset
```bash
python scripts/tcga_download_data.py \
  --data-type expression \
  --cancer-types TCGA-BRCA \
  --outdir data/test_tcga \
  --max-files 2
```

### Download clinical data for multiple cancers
```bash
python scripts/tcga_download_data.py \
  --data-type clinical \
  --cancer-types all \
  --outdir data/tcga_clinical \
  --manifest data/tcga_clinical/manifest.json
```

## GDC API query patterns

### Expression quantification
```
data_category = "Transcriptome Profiling"
data_type = "Gene Expression Quantification"
workflow_type = "STAR - Counts"
```

### Somatic mutations
```
data_category = "Simple Nucleotide Variation"
data_type = "Masked Somatic Mutation"
```

### Copy number variation
```
data_category = "Copy Number Variation"
data_type = "Gene Level Copy Number"
```

### Clinical supplements
```
data_category = "Clinical"
data_type = "Clinical Supplement"
```

### DNA methylation
```
data_category = "DNA Methylation"
data_type = "Methylation Beta-Value"
```

## Valid TCGA cancer types

Common cancer type codes:
- `TCGA-BRCA` – breast invasive carcinoma
- `TCGA-LUAD` – lung adenocarcinoma
- `TCGA-LUSC` – lung squamous cell carcinoma
- `TCGA-OV` – ovarian serous cystadenocarcinoma
- `TCGA-UCEC` – uterine corpus endometrial carcinoma
- `TCGA-COAD` – colon adenocarcinoma
- `TCGA-READ` – rectum adenocarcinoma
- `TCGA-PRAD` – prostate adenocarcinoma
- `TCGA-HNSC` – head and neck squamous cell carcinoma
- `TCGA-THCA` – thyroid carcinoma
- `TCGA-GBM` – glioblastoma multiforme
- `TCGA-LGG` – brain lower grade glioma
- `TCGA-SKCM` – skin cutaneous melanoma

For a complete list, query the GDC API or add `--list-projects` flag to the script.

## Expected outputs

### Downloaded files (in `outdir/raw/`)
- Original data files as returned by GDC (various formats: TSV, MAF, TXT)
- Named with GDC file IDs for provenance tracking

### Merged matrices (in `outdir/merged/`)

**Expression:** `expression_matrix.tsv`
- Rows: gene IDs or symbols
- Columns: TCGA sample barcodes
- Values: read counts or normalized expression
- Example: 20,000 genes × 500 samples

**Mutations:** `combined_mutations.maf`
- Standard MAF format
- All samples concatenated
- Rows: unique somatic mutations

**CNV:** `cnv_matrix.tsv`
- Rows: genes
- Columns: samples
- Values: copy number estimates (continuous or integer)

**Clinical:** `clinical_metadata.tsv`
- Rows: TCGA case/sample IDs
- Columns: clinical variables (age, stage, subtype, etc.)

### Manifest (if `--manifest` specified)
JSON file with structure:
```json
{
  "data_type": "expression",
  "cancer_types": ["TCGA-BRCA"],
  "download_timestamp": "2026-04-09T12:34:56Z",
  "outdir": "data/tcga_brca_expr",
  "files": [
    {
      "gdc_id": "abc123...",
      "filename": "...",
      "size_bytes": 1234567,
      "md5": "d41d8cd98f00b204e9800998ecf8427e",
      "status": "downloaded",
      "merge_group": "expression"
    }
  ],
  "summary": {
    "requested": 10,
    "downloaded": 8,
    "failed": 0,
    "merged": true
  }
}
```

## Failure modes

Fail with clear error message if:
- the GDC API is unreachable (network error)
- invalid `--data-type` specified
- invalid cancer type code provided
- no files match the filters (check data availability)
- output directory cannot be created
- download fails after 3 retry attempts
- insufficient disk space
- MD5 checksum mismatch after download (data corruption)

Recover gracefully by:
- retrying failed downloads up to 3 times
- skipping corrupted files and continuing
- reporting failed downloads in manifest
- providing suggestions for common issues (e.g., "No expression data found for TCGA-BRCA in GDC; verify project code")

## Best practices

- use `--max-files` when first testing a new data type
- include `--manifest` for reproducibility
- check manifest for download failures before downstream analysis
- reuse `outdir` across multiple runs (script skips existing files)
- use separate directories for different data types or cancer types
- verify merged matrix dimensions match expected sample counts

## Agent guidance

When using this skill:
- infer data type from user's analysis goals (e.g., "mutation analysis" → `--data-type mutations`)
- choose cancer types based on user's study design (specific types or `all` if exploratory)
- set `--max-files` conservatively first, increase if needed
- always include `--manifest` for downstream tracking
- if download stalls or times out, try reducing `--max-files`
