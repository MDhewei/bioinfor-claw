---
name: pan-cancer-analysis-for-gene
description: Pan-cancer analysis for a gene across all 33 TCGA cancer types plus DepMap and CPTAC. Use this skill when the user asks for pan-cancer analysis, cross-cancer comparison, multi-cancer gene report, or wants to analyze a gene across ALL cancer types at once. Generates a comprehensive PDF report with expression, survival (KM curves + log-rank), mutation, copy-number, DepMap cell-line data, and CPTAC protein data. Covers TCGA expression + survival, TCGA mutation/CNA survival, DepMap 26Q1 expression/copy-number/mutations, and CPTAC proteomics.
---

# Pan-Cancer Gene Report

## Overview

Use this skill to create a figure-backed pan-cancer report for one human gene. The script generates a multi-page PDF and CSV appendix using public TCGA/GDC Hub, cBioPortal, DepMap, and CPTAC/cBioPortal sources.

Core outputs:

- TCGA expression distributions across 33 cancer projects.
- TCGA expression-associated survival using a within-project median split.
- TCGA mutation and GISTIC copy-number alteration counts plus alteration-associated survival.
- DepMap Public 26Q1 cell-line expression, copy-number signal, expression-copy-number correlation, and mutation landscape.
- CPTAC protein abundance distributions from public cBioPortal CPTAC proteomics studies.
- A CSV appendix with summary statistics.

## Quick Start

1. Resolve the gene's human Ensembl gene ID without version, for example `PRNP -> ENSG00000171867` or `NADK -> ENSG00000008130`.
2. Run the bundled script:

```bash
python scripts/generate_pancancer_gene_report.py \
  --gene NADK \
  --ensembl ENSG00000008130 \
  --outdir ./output
```

3. The script fetches public data from TCGA/GDC Hub, DepMap, and cBioPortal (network access required).

Output files:

- `<gene>_pancancer_clinical_relevance_report_with_depmap.pdf`
- `<gene>_pancancer_appendix.csv`

## Script

The reusable script is:

```text
scripts/generate_pancancer_gene_report.py
```

Arguments:

- `--gene`: HGNC symbol, uppercase preferred.
- `--ensembl`: Ensembl gene ID without version.
- `--outdir`: Output directory (auto-created if it does not exist).
- `--keep-data`: Keep downloaded intermediate data files. By default, the script deletes all cached data after the report is built so only the PDF and CSV appendix are delivered. Pass `--keep-data` only when the user explicitly asks for the raw datasets.

The script is intentionally self-contained. Prefer running it rather than rewriting the analysis. Patch it only when the user asks for a changed report structure or additional data source.

## Interpretation Rules

Expression survival:

- Use a separate median split within each TCGA project.
- Label samples `High` when expression is `>=` that cancer-specific median.
- Label samples `Low` when expression is below the median.
- The report table's `High events` and `Low events` columns are deaths/events, not group sizes.
- Treat median-split p-values as exploratory and not multiple-testing corrected.

Alteration survival:

- TCGA mutation and copy-number statuses are fetched from cBioPortal PanCancer Atlas profiles.
- Copy-number categories use cBioPortal/GISTIC discrete calls: `2 = high-level amplification`, `1 = gain`, `0 = neutral`, `-1 = shallow deletion`, `-2 = deep deletion`.
- The script tests mutation, amplification, gain-or-amplification, and deletion separately.
- Skip or mark as insufficient when altered groups are too small for a meaningful KM/log-rank comparison.

DepMap:

- Use DepMap Public 26Q1 from the DepMap download API.
- Expression uses `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`.
- Copy number uses `PortalOmicsCNGeneLog2.csv`.
- Mutations use `OmicsSomaticMutations.csv`.
- Interpret high copy-number signals as relative screens, not clinical amplification calls.

CPTAC:

- Use public cBioPortal CPTAC `protein_quantification` profiles where the gene is present.
- Do not merge CPTAC protein values into one harmonized scale without additional normalization; present study-level distributions.

## Validation Checklist

After generating a report:

- Confirm the PDF exists and opens as a PDF.
- Confirm page ordering: expression plots, expression survival screen/KMs, expression-survival summary, TCGA alteration survival, DepMap, CPTAC, interpretation.
- Confirm the appendix has TCGA, TCGA_cBioPortal, DepMap_26Q1, and CPTAC_cBioPortal rows when data are available.
- Summarize counts for the user: TCGA projects, alteration rows, CPTAC studies/samples, and any nominal hits if relevant.

Useful QA snippet:

```bash
/Users/whe3/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 - <<'PY'
from pypdf import PdfReader
r = PdfReader("nadk_pancancer_clinical_relevance_report_with_depmap.pdf")
print("pages", len(r.pages))
for i, p in enumerate(r.pages, 1):
    txt = (p.extract_text() or "")[:240].replace("\n", " | ")
    if "Expression-Survival Summary" in txt or "Mutation and Copy-Number" in txt or "CPTAC" in txt:
        print(i, txt)
PY
```
