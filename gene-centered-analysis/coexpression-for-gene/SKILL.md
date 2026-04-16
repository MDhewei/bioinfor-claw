---
name: coexpression-for-gene
description: "Co-expression in PATIENT SAMPLES (TCGA tumors / GTEx normal tissues). NOT for cell lines — use depmap_coexpression.py for DepMap cell-line co-expression. Computes Pearson/Spearman correlations, FDR correction, GO enrichment, and network visualization. Currently uses synthetic placeholder data unless --dataset custom --expression-file is supplied."
---

# Co-expression for Gene (TCGA / GTEx Patient Samples)

> **⚠️ DATA STATUS:** This skill currently uses **synthetic placeholder data**
> when `--dataset tcga` or `--dataset gtex` is selected. The synthetic matrix
> includes real HGNC gene symbols and structured correlations so the pipeline
> runs end-to-end, but the results are NOT derived from real patient expression.
> For real analyses, supply your own expression matrix via
> `--dataset custom --expression-file <TSV>` (e.g. a UCSC Xena per-cancer matrix).

> **⚠️ SCOPE:** This skill analyzes co-expression across **patient/tissue samples**
> (TCGA tumors, GTEx normal tissues). For co-expression across **cancer cell lines**
> (DepMap), use `depmap-analysis-for-gene` → `depmap_coexpression.py` instead.

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
- **--expression-file**: Path to custom expression matrix TSV (required if dataset=custom)
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

## Procedure

### Step 1: Load Expression Data
1. **TCGA path**: Query GDC API for STAR read counts for specified project
   - Fetch log-normalized TPM values for top 200-500 samples
   - Genes as rows, samples as columns
2. **GTEx path**: Query GTEx API for tissue-specific expression
   - Use `/expression/geneExpression` endpoint
   - Filter for specified tissue
3. **Custom path**: Load user-provided TSV (genes x samples, log-scale expression)

### Step 2: Extract Query Gene
1. Retrieve expression vector for query gene
2. Remove samples with missing data
3. Check for sufficient variance (skip if constant expression)

### Step 3: Compute Correlations
1. Vectorized correlation computation:
   - For each gene, compute correlation with query gene across samples
   - Pure numpy implementation for speed
   - Skip genes with low variance
2. Calculate p-values for correlations
3. Filter for genes with valid correlation (exclude NaN/inf)

### Step 4: Multiple Testing Correction
1. Apply Benjamini-Hochberg FDR correction
2. Filter by --fdr-cutoff
3. Sort by absolute correlation value (|r|)

### Step 5: Functional Enrichment (Optional)
1. If --run-go flag:
   - POST top-N gene list to Enrichr API (`addList` endpoint)
   - Request GO Biological Process enrichment
   - Parse results and compute adjusted p-values

### Step 6: Network Visualization
1. Build network graph from top-N genes
2. Compute force-directed layout (spring layout using numpy)
3. Node properties:
   - Size proportional to |correlation|
   - Color representing correlation sign (red=positive, blue=negative)
4. Edge properties:
   - Width proportional to |correlation| between connected genes
   - Only draw edges for |r| > 0.3 (to reduce clutter)

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

### GTEx tissue-specific:
```bash
python scripts/coexpression_for_gene.py \
  --gene EGFR \
  --dataset gtex \
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

### High-stringency analysis:
```bash
python scripts/coexpression_for_gene.py \
  --gene KRAS \
  --method spearman \
  --fdr-cutoff 0.001 \
  --top-n 50 \
  --run-go \
  --outdir results/kras_stringent
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
| "Gene not found in expression matrix" | Gene symbol not recognized | Verify HGNC symbol; check case sensitivity |
| "Insufficient samples" | Dataset too small | Use larger dataset (TCGA > GTEx for sample size) |
| "Download timeout" | Network issue | Use local --expression-file if available |
| "GO enrichment timeout" | Enrichr API slow | Rerun with --run-go, or skip with --top-n < 5 |
| "No significant correlations" | Gene not co-expressed with others | Lower --fdr-cutoff or check gene validity |
| "Network visualization too dense" | Too many edges | Reduce --network-top-n or increase edge threshold |
| "Memory error" | Large expression matrix | Subset genes beforehand or use GTEx (smaller) |

## Technical Notes

- Correlations computed only on samples with complete data for both genes
- Gene variance threshold: requires SD > 0.01 in log scale to avoid spurious correlations
- FDR correction uses Benjamini-Hochberg procedure (controls false positive rate)
- Network layout: spring/force-directed layout using random walk sampling
- Edges filtered to |r| > 0.3 to improve visualization clarity
- All plots use matplotlib with high DPI (300); PNG output
- Pearson assumes linear relationships; Spearman is more robust to outliers
- Results reproducible with same dataset and seed
- GO enrichment uses Enrichr database (GO_Biological_Process_2023)
- Minimum 30 samples recommended for reliable correlation estimates

