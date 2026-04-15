---
name: gene-list-overlap
description: Compare 2 to 6 gene lists for overlap and enrichment. Computes pairwise Jaccard similarity, overlap coefficients, and Fisher's exact test for enrichment relative to a background. Generates Venn diagrams (up to 4 sets), UpSet plots, and an overlap matrix heatmap.
---

# Gene List Overlap & Enrichment Analysis

## Purpose

Compare multiple gene lists to identify shared members, compute statistical enrichment of pairwise overlaps, and visualize intersection patterns. This is essential for validating consistency across experimental replicates, comparing results from different omics platforms or screening approaches, and discovering which genes are consensus hits across conditions.

## When to Use

- You have 2–6 independent gene lists (e.g., RNA-seq DEGs, GWAS hits, CRISPR screen results, ChIP-seq targets) and want to find shared hits
- You need to quantify overlap significance via Fisher's exact test (vs. random expectation)
- You want publication-quality Venn or UpSet visualizations
- You're aggregating multi-experiment results and need to identify consensus genes
- You want to prioritize genes based on how many independent approaches selected them

## When NOT to Use

- You have > 6 gene lists: UpSet plots become hard to read; consider splitting comparisons
- Your lists are unbalanced (one list is 10x larger than others): statistical tests may be dominated by size effects; consider filtering
- You're looking for tissue-specific or condition-specific patterns in expression data: use dedicated tools (e.g., co-expression networks)
- Your background universe is very large or unknown: use generic default (20,000) or supply from external source

## Expected Inputs & Outputs

**Inputs:**
- 2–6 gene lists (files or inline comma-separated; format: `label:path_or_genes`)
- Background universe size or file (default: 20,000)
- Minimum overlap threshold (optional; default: 1)

**Outputs:**
- `overlap_summary.tsv`: pairwise statistics (Jaccard, overlap coefficient, Fisher's p-value, FDR, odds ratio)
- `jaccard_matrix.tsv`: N×N matrix of Jaccard indices for all pairs
- `unique_per_list.tsv`: genes appearing in only one list (if `--output-unique` flag)
- `all_intersections.tsv`: all pairwise and global intersection sets (if `--output-shared` flag)
- `venn_diagram.png`: Venn diagram for 2–4 lists (if `--plot-venn` flag)
- `upset_plot.png`: UpSet plot for ≥3 lists (if `--plot-upset` flag)
- `jaccard_heatmap.png`: clustered heatmap of Jaccard matrix

## Procedure

### Step 1: Input Parsing
- Parse label:source pairs
- Load genes from files (one per line) or inline (comma-separated)
- Validate and deduplicate within each list
- Map to background universe (e.g., all human protein-coding genes)

### Step 2: Pairwise Overlap Computation
For each pair (A, B):
- Intersection |A ∩ B|
- Union |A ∪ B|
- Jaccard index: |A ∩ B| / |A ∪ B|
- Overlap coefficient: |A ∩ B| / min(|A|, |B|)

### Step 3: Fisher's Exact Test
For each pair, compute 2×2 contingency table:
```
                In B    Not in B
In A            a       b           (a+b = |A|)
Not in A        c       d           (c+d = background - |A|)
             --------
             (a+c=|B|) (b+d=...)
```
Where:
- a = |A ∩ B|
- b = |A| - a
- c = |B| - a
- d = background - |A| - |B| + a

Compute:
- p-value (hypergeometric/Fisher's exact)
- Odds ratio: (a * d) / (b * c)
- Benjamini-Hochberg FDR across all pairwise tests

### Step 4: Global Overlap & Module Detection
- Identify genes appearing in k lists (for k = 1, 2, ..., N)
- Compute pairwise module assignments

### Step 5: Visualization

**Venn Diagram (2–4 lists):**
- Use matplotlib circles/ellipses with proportional overlap areas (approximate)
- Label set sizes and intersection sizes
- Color each region with alpha transparency

**UpSet Plot (≥3 lists):**
- Binary matrix: rows = genes, columns = list membership
- Bar chart: x-axis = set combinations, height = intersection size
- Sorted by descending intersection size
- Optional: show aggregated statistics (total unique, total shared, etc.)

**Jaccard Heatmap:**
- N×N symmetric matrix colored by Jaccard index (0–1)
- Clustered via hierarchical agglomerative clustering
- Dendrogram on left/top

### Step 6: Output
- TSV files with 4 decimal precision for statistics
- PNG plots at 300 dpi
- Summary report with key findings

## Key Execution Patterns

### Example 1: Compare two RNA-seq replicates
```bash
python gene_list_overlap.py \
  --lists "Rep1:rna_rep1_degs.txt,Rep2:rna_rep2_degs.txt" \
  --background 20000 \
  --plot-venn \
  --output-shared \
  --outdir ./rna_overlap
```

**Expected output:** Venn diagram with two overlapping circles; high Jaccard (0.7–0.9) indicates good replication; unique genes per replicate noted in output.

### Example 2: Integrate GWAS hits, CRISPR screen, and ChIP-seq targets
```bash
python gene_list_overlap.py \
  --lists "GWAS:gwas_hits.txt,CRISPR:screen_essential_genes.txt,ChIP:chipseq_tf_targets.txt,CoEssential:coessential_cluster.txt" \
  --background 20000 \
  --plot-upset \
  --output-unique \
  --outdir ./multi_omics_consensus
```

**Expected output:** UpSet plot showing all 15 pairwise combinations; identifies genes shared across all 4 approaches (consensus hits) vs. approach-specific findings.

### Example 3: Inline gene lists with minimal output
```bash
python gene_list_overlap.py \
  --lists "A:TP53,BRCA1,PTEN,MYC,RB1,E2F1,CDKN1A,PIK3CA,EGFR,KRAS" \
  --lists "B:TP53,BRCA1,PTEN,EGFR,AKT1,PIP4K2A" \
  --background 20000 \
  --outdir ./simple_overlap
```

**Expected output:** Small overlap_summary.tsv; Jaccard matrix.

### Example 4: Compare disease gene lists from different sources
```bash
python gene_list_overlap.py \
  --lists "ClinVar:clinvar_cancer.txt,COSMIC:cosmic_hotspots.txt,DisGeNET:disgenet_diabetes.txt,OpenTargets:opentargets_hits.txt,PanelApp:panelapp_diagnostic.txt,Monarch:monarch_disease_genes.txt" \
  --background 25000 \
  --min-overlap 2 \
  --plot-upset \
  --plot-venn \
  --outdir ./disease_databases_overlap
```

## Parameter Decision Guide

| Scenario | Recommendation |
|----------|---|
| **2 lists only (e.g., replicates)** | Use Venn 2-way; focus on Jaccard index and Fisher's p-value; set `--min-overlap 1` |
| **3 lists (typical multi-platform)** | Both Venn and UpSet; highlight consensus and list-specific genes |
| **4 lists** | Venn is still readable; UpSet is preferred; set `--output-shared` to inspect all combinations |
| **5–6 lists** | UpSet only (Venn becomes unreadable); inspect `overlap_summary.tsv` for top pairs |
| **Highly imbalanced sizes** | Jaccard is more informative than overlap coefficient; set `--min-overlap` to filter noise |
| **Background unknown** | Use biological default: 20,000 (approx. human protein-coding genes) |
| **Background known (e.g., from experiment)** | Supply as file or specific count; affects Fisher's p-value significance |
| **Publication figure** | Include both Venn/UpSet AND Jaccard heatmap; save at 300 dpi |

## Failure Modes & Solutions

| Issue | Cause | Solution |
|-------|-------|----------|
| **Parser error on input format** | Incorrect label:source syntax or missing file | Ensure format is exactly `label:path` or `label:gene1,gene2,gene3`; check file paths exist |
| **Empty gene list** | File empty or malformed; inline string parsed incorrectly | Validate input file (at least one gene per line); test inline parsing with a simpler example |
| **Fisher's p-value = 1.0 for all pairs** | Overlap expected by random chance; background set too small | Increase background size; check if genes are actually independent (e.g., same source) |
| **Venn diagram not displaying** | Too many lists (>4); matplotlib rendering issue | Use `--plot-venn` only for ≤4 lists; switch to UpSet for >4 |
| **UpSet plot: all bars same height** | Lists have identical membership; no variation in overlap | Check if lists are duplicates; inspect unique_per_list.tsv |
| **Heatmap: no clustering visible** | All Jaccard values very similar (all ~0.7 or ~0.3) | This is expected if lists are from similar sources; colors still informative |
| **Memory error on large lists** | Gene count > 50,000 per list; matrix operations slow | Consider filtering lists by statistical significance before input; script is tested up to 20,000 genes |
| **Odd ratio = inf** | One cell in contingency table is zero | This is mathematically expected; p-value is still valid; interpret as extreme odds ratio |

## Technical Notes

- **Fisher's Exact vs. Chi-Square:** Script uses Fisher's exact test (hypergeometric distribution) which is appropriate for small counts and skewed contingency tables
- **Benjamini-Hochberg FDR:** Applied across all pairwise comparisons to control for multiple testing
- **Jaccard Index Properties:** 0 = no overlap, 1 = identical sets. Not symmetric in overlap coefficient (|A ∩ B| / min(|A|, |B|))
- **Background Universe:** Larger background → smaller p-values. Use biologically realistic value. For mammalian protein-coding genes, ~19,000–20,000 is standard
- **Visualization Limitations:** Venn diagrams for 4+ sets are hard to interpret; UpSet is more scalable
- **Computational Complexity:** O(N²) for N lists; manageable for N ≤ 6. Beyond that, precompute pairwise summaries
