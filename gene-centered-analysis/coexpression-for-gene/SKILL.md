---
name: coexpression-for-gene
description: "Co-expression in PATIENT SAMPLES (TCGA tumors / GTEx normal tissues). NOT for cell lines — use depmap_coexpression.py for DepMap cell-line co-expression. Fetches real TCGA expression from cBioPortal API (PanCancer Atlas). GTEx requires user-provided expression file."
---

# Co-expression for Gene (TCGA / GTEx Patient Samples)

> **⚠️ SCOPE:** This skill analyzes co-expression across **patient/tissue samples**
> (TCGA tumors, GTEx normal tissues). For co-expression across **cancer cell lines**
> (DepMap), use `depmap-analysis-for-gene` → `depmap_coexpression.py` instead.

## Data Sources

1. **TCGA (default)** — Fetches mRNA expression (RNA-Seq V2 RSEM) from the
   cBioPortal REST API (`www.cbioportal.org`). Queries TCGA PanCancer Atlas
   studies (32 cancer types). Requires internet access.
2. **GTEx** — Requires a user-provided expression matrix file (genes × samples
   TSV). GTEx data is not available via cBioPortal.
3. **Custom** — User-provided expression matrix TSV via `--expression-file`.

## Purpose
Discover genes whose expression patterns correlate with a query gene, enabling:
- Identification of functionally related genes (co-expressed gene modules)
- Gene set enrichment and pathway analysis
- Functional annotation of uncharacterized genes
- Network-based gene prioritization for disease studies
- Tissue-specific co-expression patterns

## Use When
- You want to find genes functionally related to a query gene
- You're characterizing gene regulatory networks in specific tissues
- You need functional annotations for uncharacterized genes
- You're building disease-associated gene modules
- You want tissue-specific vs pan-cancer co-expression patterns
- You need pathway enrichment for co-expressed gene sets

## Do Not Use When
- You need causal regulatory relationships (use ChIP-seq or ATAC-seq data)
- You want direct protein-protein interactions (use STRING or BioGRID)
- You're analyzing time-series or developmental gene expression
- You need single-cell resolution (use scRNA-seq tools)
- You want cell-type specific co-expression (tissue level analysis only)
- You want co-expression across **cancer cell lines** → use `depmap-analysis-for-gene/depmap_coexpression.py`
- You want **co-essentiality** (CRISPR dependency correlations) → use `depmap-analysis-for-gene/depmap_coessentiality.py`

## Expected Inputs & Outputs

### Inputs
- **--gene** (required): Query gene symbol (e.g., TP53, BRCA1)
- **--dataset**: Source of expression data (choices: tcga, gtex, custom; default: tcga)
- **--cancer-type**: TCGA project code (default: BRCA; ignored if dataset != tcga)
- **--tissue**: GTEx tissue name (default: Breast_Mammary_Tissue; ignored if dataset != gtex)
- **--expression-file**: Path to custom expression matrix TSV (required if dataset=custom or gtex)
- **--method**: Correlation method (choices: pearson, spearman; default: pearson)
- **--top-n**: Number of top co-expressed genes to report (default: 100)
- **--fdr-cutoff**: FDR threshold for significance (default: 0.01)
- **--network-top-n**: Number of top genes to include in network visualization (default: 30)
- **--run-go**: Flag to run GO enrichment analysis using Enrichr API
- **--outdir**: Output directory

### Outputs
- `coexpression_results.tsv`: Complete results (gene, correlation, pvalue, fdr)
- `top_coexpressed.png`: Horizontal bar chart of top N co-expressed genes
- `coexpression_network.png`: Force-directed network of top 30 genes
- `go_enrichment.tsv`: GO term enrichment results (if --run-go)
- `go_bubble.png`: Bubble plot of GO enrichment (if --run-go)
- `coexpression_summary.txt`: Summary statistics and parameters

## Key Execution Patterns

### Pan-TCGA breast cancer (default):
```bash
python scripts/coexpression_for_gene.py \
  --gene BRCA1 \
  --dataset tcga \
  --cancer-type BRCA \
  --outdir results/brca1_coexpr
```

### Specific TCGA cohort with GO enrichment:
```bash
python scripts/coexpression_for_gene.py \
  --gene TP53 \
  --dataset tcga \
  --cancer-type LUAD \
  --run-go \
  --top-n 150 \
  --outdir results/tp53_luad
```

### GTEx tissue-specific (requires local file):
```bash
python scripts/coexpression_for_gene.py \
  --gene EGFR \
  --dataset gtex \
  --expression-file gtex_lung_expression.tsv \
  --tissue Lung \
  --method spearman \
  --fdr-cutoff 0.001 \
  --outdir results/egfr_lung
```

### Custom expression matrix:
```bash
python scripts/coexpression_for_gene.py \
  --gene MYC \
  --dataset custom \
  --expression-file my_expression_matrix.tsv \
  --run-go \
  --network-top-n 50 \
  --outdir results/myc_custom
```

## Parameter Decision Guide

| Scenario | --dataset | --method | --fdr-cutoff | --top-n | --run-go |
|----------|-----------|----------|--------------|---------|----------|
| Discovery (pan-TCGA) | tcga | pearson | 0.05 | 200 | yes |
| High-stringency | tcga | spearman | 0.001 | 50 | yes |
| Tissue-specific | gtex | spearman | 0.01 | 100 | yes |
| Custom cohort | custom | pearson | 0.01 | 100 | yes |
| Network analysis | tcga | pearson | 0.05 | 100 | no |
| Publication quality | tcga | spearman | 0.001 | 50 | yes |

## Failure Modes

| Error | Cause | Solution |
|-------|-------|----------|
| "Gene not found" | Gene symbol not in study | Verify HGNC symbol; try different cancer type |
| "No mRNA expression profile" | Study lacks RNA-seq data | Try a different TCGA cohort |
| "cBioPortal unreachable" | Network blocked | Use --dataset custom with local expression matrix |
| "Insufficient samples" | Dataset too small | Use larger cohort (BRCA, LUAD have 500+ samples) |
| "GO enrichment timeout" | Enrichr API slow | Rerun with --run-go, or skip |
| "No significant correlations" | Gene not co-expressed | Lower --fdr-cutoff or check gene validity |

## Technical Notes

- TCGA expression data: RNA-Seq V2 RSEM normalized values from cBioPortal PanCancer Atlas
- Correlations computed only on samples with complete data for both genes
- Gene variance threshold: requires SD > 0.01 to avoid spurious correlations
- FDR correction uses Benjamini-Hochberg procedure
- Network layout: spring/force-directed layout
- Pearson assumes linear relationships; Spearman is more robust to outliers
- GO enrichment uses Enrichr database (GO_Biological_Process_2023)
- Minimum 30 samples recommended for reliable correlation estimates
