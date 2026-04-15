---
name: tcga-survival-for-gene
description: Perform Kaplan-Meier survival analysis for a single gene in a TCGA cohort. Stratifies patients by expression level (high vs low) and tests Overall Survival (OS) and/or Disease-Free Survival (DFS) with log-rank statistics. All data is fetched from the GDC API — no local data files required.
---

# TCGA Survival Analysis for Gene

## Purpose

Analyse the **survival impact of a single gene's expression** within a specific TCGA cancer cohort.

Supported endpoints:
- `os` — Overall Survival (OS)
- `dfs` — Disease-Free Survival (DFS)
- `both` — OS and DFS simultaneously

Supported stratification methods:
- `median` — split at median expression (default, maximum sample size)
- `quartile` — top 25% = High, bottom 25% = Low (middle quartile excluded)
- `custom` — split at a user-specified expression value

## Reuse policy

This skill is designed for:
- Single-gene prognostic analysis in any TCGA cohort
- Downstream survival validation after DepMap, GTEx, or TCGA expression analysis
- Hypothesis testing: "Is high expression of gene X associated with worse survival?"

This skill does **not** require local data files.
All data (clinical + expression) is fetched live from the GDC API.

## Data sources

| Data type       | GDC endpoint                    | Notes                                      |
|-----------------|---------------------------------|--------------------------------------------|
| Gene resolution | `/genes`                        | Resolves symbol → Ensembl gene ID          |
| Clinical data   | `/cases`                        | vital_status, days_to_death, follow-up, recurrence |
| Expression      | `/gene_expression/values`       | UQFPKM or log2 median-centered UQFPKM     |

## Inputs

### Required
- `--gene` — gene symbol (e.g. `TP53`, `EGFR`, `MYC`)
- `--cancer-type` — TCGA cohort code (e.g. `BRCA`, `LUAD`, `COAD`)
  - Accepted with or without `TCGA-` prefix
- `--outdir` — output directory

### Endpoint selection
- `--mode`

  Choices: `os`, `dfs`, `both`
  Default: `os`

### Stratification
- `--stratify`

  Choices: `median`, `quartile`, `custom`
  Default: `median`

- `--custom-cutoff`

  Required only when `--stratify custom` is used.
  Value must be in the same units as `--tsv-units`.

### Expression units
- `--tsv-units`

  Choices: `uqfpkm`, `median_centered_log2_uqfpkm`
  Default: `uqfpkm`

## Outputs

### Always
- `<GENE>.<COHORT>.survival_data.tsv`

  Per-case table with columns:
  `case_id`, `submitter_id`, `vital_status`, `days_to_death`,
  `days_to_last_follow_up`, `days_to_recurrence`, `progression_or_recurrence`,
  `expression`, `expression_group`, `OS_time`, `OS_event`, `DFS_time`, `DFS_event`

- `survival_summary.tsv`

  One-row summary: gene, cohort, stratification method, n per group, log-rank p-values

- `summary.txt`

  Human-readable version of the summary

### OS mode (when `--mode os` or `--mode both`)
- `<GENE>.<COHORT>.os_km.png`
- `<GENE>.<COHORT>.os_km.pdf`

### DFS mode (when `--mode dfs` or `--mode both`)
- `<GENE>.<COHORT>.dfs_km.png`
- `<GENE>.<COHORT>.dfs_km.pdf`

## Stratification behaviour

| Method   | High group           | Low group             | Notes                            |
|----------|---------------------|-----------------------|----------------------------------|
| median   | expression ≥ median | expression < median   | All cases included               |
| quartile | top 25% (Q3+)       | bottom 25% (Q1-)      | Middle 50% excluded from plot    |
| custom   | expression ≥ cutoff | expression < cutoff   | cutoff in same units as --tsv-units |

## Plot behaviour

- KM curves are rendered with shaded 95% confidence intervals
- Log-rank p-value is annotated on each plot
- Stratification method is labelled in the lower-left corner
- Both PNG (300 dpi) and PDF are saved
- High group is shown in red, Low group in blue

## Execution policy

Construct the command from gene, cohort, mode, and optional parameters.

### Command template

```
python scripts/tcga_survival_for_gene.py \
  --gene <GENE> \
  --cancer-type <CANCER_TYPE> \
  --mode <MODE> \
  --stratify <STRATIFY> \
  [--custom-cutoff <CUTOFF>] \
  [--tsv-units <UNITS>] \
  --outdir <OUTDIR>
```

### Example commands

#### OS analysis with median stratification (default)
```
python scripts/tcga_survival_for_gene.py \
  --gene TP53 \
  --cancer-type BRCA \
  --mode os \
  --outdir results/
```

#### Both endpoints with quartile stratification
```
python scripts/tcga_survival_for_gene.py \
  --gene EGFR \
  --cancer-type LUAD \
  --mode both \
  --stratify quartile \
  --outdir results/
```

#### DFS with custom cutoff
```
python scripts/tcga_survival_for_gene.py \
  --gene MYC \
  --cancer-type COAD \
  --mode dfs \
  --stratify custom \
  --custom-cutoff 500.0 \
  --outdir results/
```

#### Log2 units
```
python scripts/tcga_survival_for_gene.py \
  --gene KRAS \
  --cancer-type PAAD \
  --mode os \
  --tsv-units median_centered_log2_uqfpkm \
  --outdir results/
```

## Parameter decision guide

Choose parameters based on the clinical question and cohort context:

| Signal in user request | Parameter to set |
|---|---|
| "overall survival" only | `--mode os` (default) |
| "disease-free survival" / "recurrence-free" | `--mode dfs` |
| "both OS and DFS" | `--mode both` |
| General question, maximize sample size | `--stratify median` (default) — uses all cases |
| "top quartile vs bottom quartile" / extreme expressors | `--stratify quartile` |
| "expression cutoff = 500" (user-specified) | `--stratify custom --custom-cutoff 500.0` |
| Using log2-normalized data (not raw UQFPKM) | `--tsv-units median_centered_log2_uqfpkm` |
| Cohort with few DFS events (e.g. BRCA) | prefer `--mode os`; DFS data is sparse in many cohorts |
| Small cohort (< 100 cases) | prefer `--stratify median` (quartile loses too many cases) |
| "lung cancer" → LUAD, "breast" → BRCA, "pancreatic" → PAAD | translate to TCGA cohort code |

**Common TCGA cohort codes**: BRCA, LUAD, LUSC, COAD, STAD, PRAD, KIRC, BLCA, GBM, OV, PAAD, SARC, UCEC, THCA, LIHC, HNSC

## Failure conditions

Fail clearly if:
- `--gene` is missing or cannot be resolved in GDC
- `--cancer-type` is missing or returns no cases from GDC
- No expression values are found for the gene in the cohort
- Clinical and expression data cannot be merged (zero matching cases)
- `--stratify custom` is requested but `--custom-cutoff` is not provided
- Output files cannot be written to `--outdir`

## Agent trigger examples

**Trigger this skill when the user asks:**
- "Is high EGFR expression associated with worse survival in lung cancer?"
- "Run KM analysis for TP53 in BRCA"
- "Show me OS curves for MYC in colon cancer stratified by expression"
- "Perform survival analysis for KRAS in TCGA PAAD"
- "Does VEGFA expression predict DFS in TCGA STAD?"
- "Run both OS and DFS for CDK4 in SARC with quartile split"

**Do NOT trigger this skill when the user asks:**
- "Show me TCGA expression levels for a gene" → use `tcga-expression-for-gene`
- "Find the top genes associated with survival" → multi-gene task, out of scope
- "Calculate hazard ratios" → use a dedicated Cox regression skill (not yet available)
- "Download TCGA clinical data files" → use `depmap-download-data` or GDC portal

## Notes

- This skill fetches data live from GDC; an internet connection is required.
- The GDC API may be slow for large cohorts (BRCA, LUAD, UCEC). Allow 1–3 minutes.
- Cohorts with < 50 matched cases may yield unreliable KM estimates.
- DFS analysis depends on GDC having recurrence/progression data, which varies by cohort.
  If DFS data is sparse, the `os` endpoint is more reliable.
- Expression values are at the **case level** (not sample level). Each patient contributes one value.
