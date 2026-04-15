---
name: gtex-download-data
description: Download GTEx v8/v10 gene expression data (bulk RNA-seq TPMs, read counts, or median expression per tissue) from the GTEx Portal API or static URLs. Returns tissue-gene matrices.
---

# GTEx Download Data

## Purpose

Download gene expression profiles from the Genotype-Tissue Expression project.

This skill supports:
- downloading GTEx median expression data (v8 or v10)
- querying tissue-specific expression for individual genes
- downloading sample-level TPM or read count matrices
- filtering by tissue type
- generating tissue × gene expression matrices

This skill is for **data acquisition only**. It does not perform differential expression or tissue comparison analysis.

## Supported data types

Expression data types available from GTEx:

- `median_tpm` – median TPM expression per gene per tissue (summary statistic)
- `sample_tpm` – sample-level TPM values (all replicates)
- `read_counts` – raw read counts from RNA-seq (lower-level data)
- `tissue_metadata` – tissue annotation and sample information
- `gene_tpm` – TPM values for specific genes across tissues (via API)

## When to use

Use this skill when you need to:
- obtain gene expression profiles across normal tissues (GTEx)
- compare gene expression between tissue types
- download median expression for tissue-gene integrative analysis
- retrieve sample-level expression data for statistical modeling
- access tissue metadata and sample annotations

This skill may also be called internally by analysis or integration skills.

## When not to use

Do not use this skill when you want:
- differential expression analysis between conditions
- tissue-specific regulatory network inference
- expression quantitative trait loci (eQTL) analysis
- disease-specific expression data (use TCGA for cancer)
- coexpression network analysis results

Use a GTEx analysis skill for those tasks.

## Data source

Primary data source:
- GTEx Portal: `https://gtexportal.org/`
- GTEx Portal API: `https://gtexportal.org/api/v2/`
- GTEx Analysis static files (GCP): `https://storage.googleapis.com/gtex_analysis_v*/`

The dataset provides:
- RNA-seq expression data from 930+ donors
- 54 tissue types
- no authentication required for public data
- structured query API for programmatic access
- static release files (v8 and v10)

## Parameters

### `--data-type`
Type of expression data to download.

Choices:
- `median_tpm` – tissue-level median TPM (summary, recommended for quick access)
- `sample_tpm` – sample-level TPM values
- `read_counts` – raw read counts
- `tissue_metadata` – tissue and sample information
- `gene_tpm` – expression for specific genes

Required: yes
Default: none (must be specified)

### `--genes`
Comma-separated list of gene symbols to download.

Use when:
- interested in specific genes (e.g., TP53, BRCA1)
- combining with `--data-type gene_tpm`
- reducing data volume for targeted analysis

Format:
- gene symbols: `TP53,BRCA1,EGFR`
- case-insensitive
- matched to GTEx gene annotations

If omitted with `gene_tpm`:
- all genes are returned (very large download)

Example:
- `--genes TP53,KRAS,EGFR` – three genes only

### `--tissues`
Comma-separated list of tissue types to download.

Use when:
- interested in specific tissues (e.g., liver, brain)
- reducing data volume

Format:
- tissue names as shown in GTEx (with spaces/dashes): `Liver,Brain - Cortex,Muscle - Skeletal`
- case-sensitive; must match GTEx exactly
- use `--list-tissues` to see valid names

If omitted:
- all tissues are returned

Example:
- `--tissues Liver,Kidney,Heart` – three tissues only

### `--version`
GTEx analysis version to download.

Choices:
- `v8` – GTEx Analysis v8 (default, larger dataset)
- `v10` – GTEx Analysis v10 (latest)

Default: v8

### `--outdir`
Output directory for downloaded files.

Use when:
- files need to be cached for reuse
- downstream tools need predictable file locations

Example:
- `data/gtex_expression`

Recommended:
- use separate directories for different versions or data types

### `--list-tissues`
Print list of available GTEx tissue types and exit.

Use when:
- exploring available tissues
- constructing `--tissues` parameter

### `--manifest`
Optional path to write manifest JSON with file metadata.

Use when:
- downstream analysis needs file paths
- reproducibility requires tracking data provenance

## Parameter decision guide

| Goal | Data type | Version | Tissues | Notes |
|------|-----------|---------|---------|-------|
| Quick tissue × gene summary | `median_tpm` | v8 | omitted | smallest download, typical use |
| Compare 2-3 genes across tissues | `gene_tpm` | v8 | omitted | fast, via API |
| Tissue-specific analysis | `median_tpm` | v8 | `Liver,Brain` | filter to tissues of interest |
| Raw read counts for modeling | `read_counts` | v8 | omitted | larger file, statistical testing |
| Sample-level replicates | `sample_tpm` | v8 | omitted | full dataset, memory intensive |
| Tissue metadata only | `tissue_metadata` | v8 | omitted | ~1 MB, informational |

## Procedure

1. **Validate parameters**: check data type, version, tissues, genes
2. **Query GTEx Portal API** or identify static file URL
3. **Download file** from GTEx storage or API (streaming for large files)
4. **Parse data format**:
   - GCT format: skip header lines, parse gene × sample matrix
   - TSV format: parse directly with pandas
5. **Filter by tissue** if specified (optional)
6. **Filter by gene** if specified (optional)
7. **Save to TSV** with proper header and index
8. **Write manifest** JSON with metadata
9. **Report**: file size, matrix dimensions, tissues/genes included

## Key execution patterns

### Download median TPM for all tissues (v8)
```bash
python scripts/gtex_download_data.py \
  --data-type median_tpm \
  --version v8 \
  --outdir data/gtex_v8
```

### Get median TPM for specific tissues
```bash
python scripts/gtex_download_data.py \
  --data-type median_tpm \
  --version v8 \
  --tissues "Liver,Kidney,Heart" \
  --outdir data/gtex_v8_tissue_subset
```

### Download expression for specific genes via API
```bash
python scripts/gtex_download_data.py \
  --data-type gene_tpm \
  --version v8 \
  --genes TP53,BRCA1,EGFR \
  --outdir data/gtex_target_genes
```

### List available tissues
```bash
python scripts/gtex_download_data.py \
  --list-tissues
```

### Download with manifest for downstream use
```bash
python scripts/gtex_download_data.py \
  --data-type median_tpm \
  --version v8 \
  --outdir data/gtex_v8 \
  --manifest data/gtex_v8/manifest.json
```

## GTEx Portal API patterns

### Median expression endpoint
```
GET /expression/medianTranscriptExpression
Query: geneId (ENSG...), tissueSiteDetailId (tissue code), pageSize, offset
```

### Gene expression endpoint
```
GET /expression/geneExpression
Query: geneId (ENSG), tissueSiteDetailId, datasetId
```

### Tissue sites endpoint
```
GET /reference/tissueSites
Returns: all tissue metadata (names, codes, descriptions)
```

## Available GTEx tissues (v8/v10)

Common tissue types:
- **Brain**: Brain - Amygdala, Brain - Anterior cingulate cortex, Brain - Caudate, Brain - Cerebellar Hemisphere, Brain - Cerebellum, Brain - Cortex, Brain - Frontal Cortex (BA9), Brain - Hippocampus, Brain - Hypothalamus, Brain - Nucleus accumbens (basal ganglia), Brain - Putamen (basal ganglia), Brain - Spinal cord (cervical c-1), Brain - Substantia nigra
- **Heart**: Heart - Atrial Appendage, Heart - Left Ventricle
- **Liver**: Liver
- **Kidney**: Kidney - Cortex
- **Muscle**: Muscle - Skeletal
- **Thyroid**: Thyroid
- **Blood**: Whole Blood
- **Gastrointestinal**: Colon - Sigmoid, Colon - Transverse, Esophagus - Gastroesophageal Junction, Esophagus - Muscularis, Esophagus - Mucosa, Small Intestine - Terminal Ileum, Stomach
- **Reproductive**: Ovary, Prostate, Testis, Uterus, Vagina
- **Other**: Adrenal Gland, Pancreas, Pituitary, Skin - Not Sun Exposed (Suprapubic), Skin - Sun Exposed (Lower leg), Spleen, Adipose - Subcutaneous, Adipose - Visceral (Omentum), Artery - Aorta, Artery - Coronary, Artery - Tibial, Breast - Mammary Tissue, Cells - EBV-transformed lymphocytes, Cells - Leukemia cell line (CML), Nerve - Tibial, Salivary Gland

Use `--list-tissues` for the complete current list.

## Expected outputs

### Median TPM matrix (if `--data-type median_tpm`)
File: `expression_matrix.tsv`
- Rows: gene IDs (ENSG) or symbols
- Columns: tissue types
- Values: median TPM (transcripts per million)
- Example: 55,000 genes × 54 tissues

### Sample TPM matrix (if `--data-type sample_tpm`)
File: `sample_expression_matrix.tsv`
- Rows: genes
- Columns: GTEx sample IDs
- Values: TPM
- Example: 55,000 genes × 10,000+ samples (large file)

### Read counts matrix (if `--data-type read_counts`)
File: `read_counts_matrix.tsv`
- Rows: genes
- Columns: samples
- Values: raw read counts
- Example: 55,000 genes × samples

### Tissue metadata (if `--data-type tissue_metadata`)
File: `tissue_metadata.tsv`
- Rows: tissues
- Columns: tissue properties (ID, description, sample count, etc.)

### Gene expression (if `--data-type gene_tpm` with specific genes)
File: `gene_tpm.tsv`
- Rows: tissues
- Columns: specified genes
- Values: median TPM

### Manifest (if `--manifest` specified)
JSON file with:
- version, data type, download date
- tissues included, genes included
- local file path, size
- number of samples, number of genes

## Failure modes

Fail with clear error message if:
- GTEx Portal/API is unreachable
- invalid `--data-type` specified
- invalid tissue names provided (not in GTEx)
- invalid gene symbols or IDs provided
- unsupported `--version`
- output directory cannot be created
- download fails and cannot be retried
- GCT format parsing fails (unexpected header)

Recover gracefully by:
- retrying failed downloads up to 3 times
- skipping invalid tissues/genes and continuing
- printing available tissues when `--tissues` has no matches
- providing suggestions for common issues

## Best practices

- use `--list-tissues` before specifying `--tissues` to avoid typos
- start with `--data-type median_tpm` for most analyses (smaller, sufficient for summaries)
- include `--manifest` for reproducibility and downstream tracking
- reuse cached data from `--outdir` (script skips existing files)
- use `--tissues` to reduce data volume for large analyses
- use `--genes` with `gene_tpm` API queries for targeted access
- verify tissue names exactly (case-sensitive): "Brain - Cortex" not "brain-cortex"

## Agent guidance

When using this skill:
- infer data type from user's analysis goals (e.g., "compare 3 genes across tissues" → `gene_tpm`)
- choose version based on requirements (v8 is default, larger sample set)
- use `--list-tissues` if user mentions specific tissues but format is unclear
- use tissue filtering for faster downloads if only a few tissues are needed
- use gene filtering with `gene_tpm` for fast, targeted queries
- always include `--manifest` when outputs feed downstream analysis
