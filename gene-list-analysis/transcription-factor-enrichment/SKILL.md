---
name: transcription-factor-enrichment
description: Find transcription factors (TFs) that significantly regulate a given gene list. Uses enrichment analysis against curated TF-target databases (ENCODE ChEA3 API, DoRothEA gene sets built-in) and returns ranked TFs with enrichment statistics and network visualization.
---

# Transcription Factor Enrichment Analysis

## Purpose

Identify which transcription factors (TFs) are likely to regulate the genes in your query list. This is essential for understanding the transcriptional drivers of your biological process, predicting regulatory networks, and prioritizing TFs for functional validation (e.g., CRISPRi repression, ChIP-seq, ATAC-seq).

## When to Use

- You have a set of co-regulated genes (e.g., DEGs from a condition, GWAS hits, genes in a pathway) and want to find their putative TF regulators
- You're studying a disease and want to identify dysregulated TFs
- You want to prioritize TFs for experimental perturbation (CRISPR, drug screening)
- You're integrating multiple omics layers (RNA-seq, ChIP-seq, ATAC-seq) and need TF predictions
- You want to predict cell type-specific TFs for lineage commitment studies

## When NOT to Use

- Your gene list is very small (< 5 genes): statistical power is low; consider literature review instead
- You need detailed TF binding site sequences or motif information: use motif discovery tools (MEME, HOMER)
- You're studying prokaryotes: DoRothEA and ChEA3 are human/mouse focused
- You have tissue-specific requirements: consider supplementing with tissue-specific TF databases (e.g., GTEx)
- You want trans-acting regulatory factors (miRNAs, signaling proteins): this tool focuses on TF regulation

## Expected Inputs & Outputs

**Inputs:**
- Gene list (comma-separated or file; required)
- Species (human or mouse; default: human)
- Analysis method (ChEA3, DoRothEA, or both; default: both)
- DoRothEA confidence levels (A,B,C,D,E; default: A,B for high confidence)
- Background gene universe (default: 20,000)
- FDR cutoff for significance (default: 0.05)
- Top N TFs to report (default: 25)

**Outputs:**
- `tf_enrichment.tsv`: All significant TFs with statistics (p-value, FDR, odds ratio, source)
- `top_tfs.png`: Horizontal bar chart of top 25 TFs by -log10(FDR)
- `tf_gene_network.png`: Bipartite network showing top 10 TFs and their overlapping target genes
- `tf_enrichment_summary.txt`: Text summary with top 10 TFs
- Console output: immediate report of top TF hits

## Procedure

### Step 1: Input Validation & Gene List Prep
- Load gene list from file or inline string
- Validate against known gene symbols
- Convert to uppercase for database matching

### Step 2: ChEA3 Analysis (if enabled)
- POST request to ChEA3 API: `https://maayanlab.cloud/chea3/api/enrich/`
- Request body: `{"query_name": "query", "gene_set": [list of genes]}`
- Parse response: extract TF rankings from multiple curated libraries:
  - ENCODE TF ChIP-seq peaks
  - GTEx regulatory regions
  - ARCHS4 TF co-expression
  - Enrichr TF target libraries
- Aggregate scores across libraries (mean or harmonic mean)
- Extract p-values if available; compute FDR

### Step 3: DoRothEA Analysis (if enabled)
- Use built-in human/mouse TF-target regulons (compact dictionary)
- Filter by confidence level (A = highest, E = lowest; default: A,B)
- For each TF:
  - Count overlap between TF's known targets and input gene list
  - Build 2×2 contingency table for Fisher's exact test
  - Compute p-value, odds ratio

### Step 4: Result Integration
- If both methods enabled: merge results by TF
  - Average or rank-aggregate p-values across methods
  - Combine target overlap counts
  - Flag TFs from each source
- Apply Benjamini-Hochberg FDR correction across all TFs

### Step 5: Network Visualization
- Build bipartite graph: top 10 TFs (left) → target genes (right)
- Edge thickness proportional to overlap count
- Node size proportional to degree or specificity
- Color code TFs by FDR significance

### Step 6: Bar Chart Visualization
- Horizontal bar chart sorted by -log10(FDR)
- Color by significance threshold (FDR < 0.05 = green, 0.05–0.1 = orange, > 0.1 = gray)
- Show TF names and p-values as labels

### Step 7: Output & Reporting
- Write TSV with all TFs (sorted by p-value)
- Print top 10 to console
- Generate summary text report with interpretation guidance

## Key Execution Patterns

### Example 1: Quick enrichment of DEG list using both methods
```bash
python transcription_factor_enrichment.py \
  --genes degs_from_condition.txt \
  --species human \
  --method both \
  --confidence A,B \
  --fdr-cutoff 0.05 \
  --top-n 25 \
  --outdir ./tf_enrichment_results
```

**Expected output:** 10–40 significant TFs (depending on gene list size and coherence); likely includes known regulators of your condition.

### Example 2: Mouse data with DoRothEA only (ChEA3 API unavailable)
```bash
python transcription_factor_enrichment.py \
  --genes ~/my_genes.txt \
  --species mouse \
  --method dorothea \
  --confidence A,B,C \
  --background 20000 \
  --outdir ./tf_mouse
```

**Expected output:** DoRothEA database for mouse; slower but reliable offline analysis.

### Example 3: Very stringent (high-confidence TFs only)
```bash
python transcription_factor_enrichment.py \
  --genes TP53,BRCA1,MDM2,CDKN1A,PMAIP1,BBC3 \
  --species human \
  --method dorothea \
  --confidence A \
  --fdr-cutoff 0.01 \
  --top-n 10 \
  --outdir ./tf_strict
```

**Expected output:** Very few hits (2–5), but high confidence; TP53 should be top hit if input list is TP53-regulated.

### Example 4: Lenient analysis (exploratory, lower stringency)
```bash
python transcription_factor_enrichment.py \
  --genes large_gene_list.txt \
  --species human \
  --method both \
  --confidence A,B,C,D \
  --fdr-cutoff 0.1 \
  --background 25000 \
  --outdir ./tf_exploratory
```

**Expected output:** Larger set of hits (50+); useful for hypothesis generation; check results carefully.

## Parameter Decision Guide

| Scenario | Recommendation |
|----------|---|
| **Gene list < 10 genes** | Use high-confidence TFs only (--confidence A); stricter FDR (0.01); focus on top 5 |
| **Gene list 10–100 genes** | Default settings (--confidence A,B, --fdr 0.05); good balance |
| **Gene list > 100 genes** | Can relax confidence (A,B,C); but risk of noise; filter by effect size (odds ratio > 2) |
| **RNA-seq DEGs** | Use ChEA3 + DoRothEA; likely to find known drivers |
| **GWAS hits** | Can use either method; expect fewer/weaker hits than RNA-seq (more heterogeneous) |
| **CRISPR screen essential genes** | Both methods; often find metabolic TFs (e.g., MYC, TP53) |
| **Tissue-specific question** | Supplement DoRothEA with tissue-specific DB (GTEx in ChEA3); post-filter by tissue |
| **Publication quality** | Both methods; FDR < 0.05; top 10–25 TFs; include network and bar chart |
| **Exploratory/hypothesis generation** | Either method; FDR < 0.1; larger top-N (50+); use for candidate prioritization |
| **Offline/reproducible** | Use DoRothEA only (no API dependency) |

## Failure Modes & Solutions

| Issue | Cause | Solution |
|-------|-------|----------|
| **ChEA3 API request fails** | Network error or API down; gene identifiers not recognized | Fall back to `--method dorothea`; check internet connection; test with known genes (TP53, MYC) |
| **No significant TFs found** | Gene list too small or heterogeneous; score cutoff too high | Relax FDR threshold (0.1 or 0.2); use lower confidence TFs (add D,E); increase background size |
| **Top TF is unexpected/wrong** | Could indicate method artifact or non-TF regulatory mechanism | Check gene list for batch effects; filter by manual curation; inspect network to identify outliers |
| **High variance between ChEA3 and DoRothEA** | Methods use different data sources and algorithms | Average results; highlight consensus TFs in visualization; note discrepancies in report |
| **Memory error** | Gene list very large (> 50,000); or API response is huge | Pre-filter gene list to top 5,000 by statistical significance; reduce --top-n |
| **All p-values = 1.0** | Gene list not in database; wrong species selection; background too small | Check species matches input genes; validate gene symbols; increase background size |
| **Odds ratio = inf** | Contingency table has zero cell; TF targets entirely within input list | This is expected for very specific TFs; treat as strong signal; note in report |
| **Parsing error on gene names** | Non-ASCII characters or spaces in input | Use simple gene symbols (e.g., TP53 not "TP53 gene"); ensure file is UTF-8 encoded |

## Technical Notes

- **ChEA3 Data Sources:** ChEA3 aggregates from ENCODE, GTEx, ARCHS4, Enrichr, and literature-curated databases. Strength is comprehensive but aggregation can reduce signal.
- **DoRothEA Confidence Levels:**
  - **A:** TF-target pairs curated from literature, ChIP-seq, and reporter assays (highest confidence)
  - **B:** Interactions inferred from ATAC-seq, RNA-seq (medium-high)
  - **C–E:** Predicted from sequence and computational methods (lower confidence)
- **Built-in DoRothEA:** Script includes a compact ~50 TF × 20–30 targets dictionary for human/mouse (most common TFs). Full DoRothEA (~500 TFs) available via external database if needed.
- **Fisher's Exact Test:** Used for statistical significance of overlap (vs. random expectation). Appropriate when target count is small.
- **FDR Correction:** Benjamini-Hochberg correction applied across all TFs tested to control false discovery rate.
- **Background Universe:** Default 20,000 approximates human protein-coding genes. Affects p-value stringency. Use larger value if your gene list is selected from a specific large catalog.
- **Scalability:** Script handles up to 10,000 input genes. Beyond that, consider pre-filtering or local methods.
- **Reproducibility:** DoRothEA results are deterministic; ChEA3 may vary if API is updated. For publication, lock API version or use local DoRothEA database.

## Interpretation Tips

- **Top TFs are known regulators of your process?** Good sign; method is capturing expected biology.
- **Top TF is TP53, MYC, or another ubiquitous factor?** Could indicate cell-stress response or artifact. Check for batch effects; filter results by tissue specificity.
- **TF specificity is high (unique targets)?** Strong candidate for functional perturbation.
- **Multiple TFs in same family (e.g., STAT1, STAT2, STAT3)?** Suggests coordinated pathway activation; consider combinatorial regulation.
- **Network shows hub TFs regulating most genes?** Central regulators identified; candidates for master regulators.
- **Consensus between ChEA3 and DoRothEA?** Higher confidence in result; highlight for validation.
