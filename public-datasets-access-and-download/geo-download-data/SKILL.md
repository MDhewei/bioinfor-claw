---
name: geo-download-data
description: Download a GEO dataset by accession (GSE or GDS) using the NCBI GEO API and FTP. Extracts the expression matrix, sample metadata, and series metadata into analysis-ready TSV files.
---

# GEO Download Data

## Purpose

Download public microarray and RNA-seq datasets from NCBI Gene Expression Omnibus (GEO).

This skill supports:
- downloading GEO series data (GSE accessions)
- downloading GEO datasets (GDS accessions)
- extracting expression matrices from series_matrix.txt files
- parsing sample-level metadata
- retrieving series-level annotations
- organizing outputs into analysis-ready TSV files

This skill is for **data acquisition only**. It does not perform normalization, batch correction, or differential expression analysis.

## Supported accession types

GEO accession formats:

- `GSE*` – GEO Series: collection of samples from a single study
  - Example: GSE12345 (breast cancer expression study)
  - contains multiple GSM (sample) entries
  - includes sample and series metadata
- `GDS*` – GEO DataSet: curated expression dataset
  - Example: GDS100 (prepared microarray dataset)
  - pre-processed expression matrix
  - smaller, more standardized format

## When to use

Use this skill when you need to:
- download a published microarray expression study
- obtain sample metadata and clinical annotations
- access expression matrices from GEO
- retrieve a specific GEO accession for analysis
- batch download multiple GEO studies

This skill may also be called internally by analysis or integration skills.

## When not to use

Do not use this skill when you want:
- differential expression analysis
- quality control and normalization
- batch correction between studies
- integrative analysis across multiple GEO studies
- processed analysis results (use individual study publications)

Use a GEO analysis skill for those tasks.

## Data source

Primary data source:
- NCBI Gene Expression Omnibus: `https://www.ncbi.nlm.nih.gov/geo/`
- GEO FTP server: `ftp.ncbi.nlm.nih.gov/geo/`
- NCBI eUtils API: `https://eutils.ncbi.nlm.nih.gov/`

The dataset provides:
- no authentication required for public data
- series_matrix.txt.gz: expression matrix + sample metadata
- SOFT files: raw format with detailed metadata
- series information: title, summary, organism, platform
- sample attributes: biological and technical metadata

## Parameters

### `--accession`
GEO accession number to download.

Format:
- `GSE*` – GEO Series (most common)
- `GDS*` – GEO DataSet (pre-curated)

Required: yes

Examples:
- `GSE12345` – breast cancer transcriptome study
- `GDS100` – prepared microarray dataset
- `GSE1234` – note: no leading zeros

### `--outdir`
Output directory for downloaded and parsed files.

Use when:
- organizing downloaded data for analysis
- creating reproducible analysis pipelines

Example:
- `data/geo_gse12345`

Recommended:
- use separate directories per GEO accession
- the script creates subdirectories for raw and parsed files

### `--soft`
Download SOFT file in addition to series matrix.

Use when:
- detailed raw metadata is needed
- debugging download issues
- archiving complete original data

Default: false

### `--matrix`
Download and parse series matrix file.

Use when:
- expression data is needed for analysis

Default: true

### `--manifest`
Optional path to write manifest JSON with download summary.

Use when:
- downstream analysis needs metadata
- reproducibility requires tracking data provenance
- integration with other tools

## Parameter decision guide

| Goal | Accession | Include SOFT | Include Matrix | Notes |
|------|-----------|--------------|----------------|-------|
| Download expression study | GSE12345 | no | yes | typical use |
| Get with detailed metadata | GSE12345 | yes | yes | archival completeness |
| Dataset only (pre-curated) | GDS100 | no | yes | smaller, standardized |
| Quick download + metadata | GSE12345 | no | yes | default settings |
| Minimal download | GSE12345 | no | yes | expression only |

## Procedure

1. **Validate accession**: check format (GSE* or GDS*)
2. **Query NCBI eUtils** to get series metadata (title, organism, samples)
3. **Construct FTP URL** for series_matrix.txt.gz file
   - Path pattern: `/geo/series/GSE{prefix}nnn/{accession}/matrix/{accession}_series_matrix.txt.gz`
   - Where prefix = first N digits of accession with nnn suffix
   - Example: GSE12345 → GSE12nnn
4. **Download** series_matrix.txt.gz with streaming
5. **Decompress** gzip file
6. **Parse series_matrix format**:
   - Lines starting with `!Series_` → series metadata
   - Lines starting with `!Sample_` → sample metadata (one per sample)
   - `!series_matrix_table_begin` to `!series_matrix_table_end` → expression matrix
7. **Extract outputs**:
   - expression_matrix.tsv: probes/genes × samples
   - sample_metadata.tsv: samples × attributes
   - series_info.json: series metadata (title, summary, organism, etc.)
8. **Optional: download SOFT file** for additional details
9. **Write manifest** JSON with summary
10. **Report**: accession, title, samples, matrix dimensions, download status

## Key execution patterns

### Download GSE study with default settings
```bash
python scripts/geo_download_data.py \
  --accession GSE12345 \
  --outdir data/geo_gse12345
```

### Download with SOFT file for complete archival
```bash
python scripts/geo_download_data.py \
  --accession GSE12345 \
  --outdir data/geo_gse12345 \
  --soft
```

### Download and generate manifest for downstream use
```bash
python scripts/geo_download_data.py \
  --accession GSE12345 \
  --outdir data/geo_gse12345 \
  --manifest data/geo_gse12345/manifest.json
```

### Download GDS (pre-curated dataset)
```bash
python scripts/geo_download_data.py \
  --accession GDS100 \
  --outdir data/geo_gds100
```

### Expression matrix only (no metadata processing)
```bash
python scripts/geo_download_data.py \
  --accession GSE12345 \
  --outdir data/geo_gse12345 \
  --matrix
```

## GEO series_matrix.txt format

Format specification:
```
!Series_title = "Study Title"
!Series_summary = "..."
!Series_overall_design = "..."
!Series_pubmed_id = "12345678"
!Series_contact_name = "Author"
!Series_contact_email = "author@institution.edu"
!Sample_title = "Sample 1"
!Sample_title = "Sample 2"
...
!series_matrix_table_begin
ID_REF    sample1    sample2    ...
PROBE_1   10.5       12.3       ...
PROBE_2   8.2        9.1        ...
...
!series_matrix_table_end
```

Key features:
- metadata lines: `!<Key> = <Value>`
- expression data: tab-separated, probe/gene × sample
- comments: lines starting with `#` are ignored
- multiple samples: each sample has separate `!Sample_*` lines

## FTP URL construction

For GSE accession `GSE12345`:
1. Extract prefix: "12345" → prefix = "12"
2. Build series folder: `GSE12nnn`
3. Complete URL:
   - `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE12nnn/GSE12345/matrix/GSE12345_series_matrix.txt.gz`

For other accessions:
- `GSE100` → `GSE100nnn`
- `GSE1000000` → `GSE1000000nnn`

## Expected outputs

### Expression matrix (required)
File: `expression_matrix.tsv`
- Rows: probe IDs or gene symbols (ID_REF column)
- Columns: sample GEO IDs (GSM numbers)
- Values: expression values (normalized or raw, platform-dependent)
- Example: 10,000 probes × 50 samples

### Sample metadata (extracted from series matrix)
File: `sample_metadata.tsv`
- Rows: sample GEO IDs
- Columns: sample attributes
- Common columns:
  - title: sample name/description
  - source_name_ch1: biological source
  - characteristics_ch1: biological annotations
  - treatment_protocol_ch1: experimental treatment
  - growth_protocol_ch1: sample preparation
- Example: 50 samples × 15 attributes

### Series information (queried via eUtils)
File: `series_info.json`
- JSON with fields:
  - accession: GSE accession
  - title: study title
  - summary: brief description
  - overall_design: study design
  - organism: species (e.g., Homo sapiens)
  - platform: array platform (e.g., GPL570)
  - pubmed_id: publication reference
  - samples: list of sample GEO IDs
  - contact_name, contact_email: author info

### SOFT file (optional, if --soft)
File: `{accession}_family.soft.gz`
- Raw SOFT format with complete metadata
- Larger file, archival purposes

### Manifest (if --manifest specified)
JSON file with:
- accession, title, organism
- number of samples, number of probes/genes
- file paths, sizes, download status
- download timestamp

## Failure modes

Fail with clear error message if:
- NCBI eUtils API is unreachable
- invalid accession format provided
- accession does not exist in GEO
- FTP download fails (server issue, network error)
- series_matrix.txt.gz not found (may be older/removed series)
- output directory cannot be created
- gzip decompression fails (corrupt file)
- series_matrix parsing fails (unexpected format)

Recover gracefully by:
- retrying failed downloads up to 3 times
- providing FTP URL for manual download if automated fails
- skipping optional SOFT download if it fails
- reporting what was successfully downloaded
- suggesting related accessions if exact match not found

## Best practices

- use `--manifest` for reproducibility and downstream tool integration
- validate accession format before downloading (GSE* or GDS* format)
- check NCBI GEO website to verify accession exists and is public
- use separate `--outdir` per accession to avoid file conflicts
- reuse `--outdir` across runs (script skips existing files)
- include `--soft` only when detailed metadata tracking is important
- verify sample count in output matches expected study size
- check expression matrix dimensions (rows = probes, columns = samples)

## Agent guidance

When using this skill:
- infer accession from user context (e.g., "the breast cancer study GSE12345")
- always include `--manifest` when downstream analysis will use outputs
- add `--soft` if user mentions needing complete archival or detailed metadata
- use separate directories for different GEO accessions
- if download fails, suggest checking the GEO website to verify accession
- if series_matrix not found, suggest checking if accession uses alternative format
- validate that retrieved sample count matches user's expectations

## Common GEO accessions

Example well-known datasets:
- GSE12345 – Breast cancer transcriptome
- GSE20986 – Cancer Gene Atlas
- GSE3494 – Ovarian cancer expression study
- GSE6891 – Lymphoma classification
- GDS100 – Ovarian cancer microarray (curated)
- GDS200 – Cancer cell lines (curated)

Use NCBI GEO search to find accessions matching your study type.

## Integration with downstream tools

The skill outputs are compatible with:
- Standard expression analysis (DESeq2, limma, etc.)
- Machine learning pipelines
- Batch correction tools (ComBat, SVA)
- Gene set enrichment (GSEA)
- Network analysis tools
- Clustering and visualization tools

All output matrices are TSV format with standard headers for easy import.
