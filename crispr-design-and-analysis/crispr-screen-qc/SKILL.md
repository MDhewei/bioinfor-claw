---
name: crispr-screen-qc
description: Quality control analysis for pooled CRISPR screens. Accepts a sgRNA read count matrix (guides × samples) and computes library representation, read depth, guide-level QC metrics, replicate correlations, and Gini index. Generates a QC report with plots.
tags:
  - crispr
  - screen
  - qc
  - quality-control
  - pooled-library
version: 1.0.0
---

## Purpose

Perform comprehensive quality control on pooled CRISPR screening data. This tool analyzes raw sgRNA read counts across samples to detect issues with library representation, sequencing depth, replicate consistency, and guide saturation. Identifies under-represented guides, batch effects, and data quality issues before downstream analysis (e.g., hit calling, phenotypic analysis).

## When to Use

**Use this skill when you:**
- Have completed DNA sequencing of a pooled CRISPR screen library
- Need to validate sample quality before proceeding to hit calling
- Want to detect replicates with low coverage or high variance
- Are checking for guides with insufficient reads across all conditions
- Need to document library distribution and representation metrics
- Working with screens in mammalian cells (human/mouse)

**Do NOT use this skill when you:**
- Analyzing single guide experiments — use instead for pools of ≥100 guides
- Need to call hits or perform association tests — this is QC only
- Working with very deep sequencing (>1M reads per guide) where Gini is less informative
- Analyzing arrayed screens (one guide per well) — designed for pooled only
- Samples have severe batch effects requiring plate/cohort corrections

## Expected Inputs

- **Count matrix TSV**: rows=sgRNAs, columns=samples; includes sgRNA ID and Gene columns
- **Sample metadata** (optional): condition names, replicate numbers
- **Thresholds** (optional): minimum reads per guide, minimum library representation fraction

## Expected Outputs

- **QC summary TSV**: one row per sample with metrics (read depth, Gini, representation %, control stats)
- **Guide-level QC TSV**: per-guide mean, CV, min, max reads across replicates
- **Correlation heatmap**: Pearson/Spearman correlations between all sample pairs
- **Distribution plots**: read counts, Gini index, CDF, control vs targeting guides
- **Pass/Fail flags**: indicates samples failing defined QC thresholds

## Procedure

1. **Load count matrix**
   - Read TSV with sgRNA ID (first column), Gene annotation (second), sample counts (remaining)
   - Remove rows with all zeros
   - Filter samples if --samples specified

2. **Compute library representation**
   - Count guides with ≥ min_reads threshold
   - Calculate fraction of guides passing threshold
   - Flag if <80% (default) of library represented

3. **Compute read depth metrics**
   - Total reads per sample
   - Median reads per guide
   - Mean reads per guide
   - Min/max reads per guide

4. **Calculate Gini index**
   - Measure of count evenness: 0 (perfectly uniform) to 1 (maximally skewed)
   - Formula: Gini = (2 * sum(i * x_i)) / (n * sum(x_i)) - (n + 1) / n
     where x_i = sorted read counts, n = number of guides
   - High Gini (>0.3) indicates some guides vastly over-represented

5. **Detect missing guides**
   - Count and fraction of guides with zero reads per sample
   - Flag if >10% of library absent

6. **Compute replicate correlations**
   - Normalize to log2(CPM + 1) within each sample
   - Pearson and Spearman correlations between all sample pairs
   - Generate correlation matrix heatmap

7. **Analyze control guides**
   - Extract non-targeting controls (prefix matching --control-prefix)
   - Compare mean reads and CV vs. gene-targeting guides
   - Flag if control CV >> targeting CV (library amplification bias)

8. **Per-guide statistics**
   - Mean reads across replicates of same condition
   - Coefficient of variation (CV) across replicates
   - Min/max reads
   - Flag guides with CV > 0.8 (poor replicate consistency)

9. **Generate visualizations**
   - Violin/box plots of read distributions per sample
   - Scatter plots of replicate correlations (all pairs)
   - Bar chart of Gini index per sample
   - Cumulative distribution function (CDF) of reads per guide
   - Heatmap of Pearson correlations
   - Box plots: control vs targeting guide read distributions

## Key Execution Patterns

### Basic QC on raw counts
```bash
python scripts/crispr_screen_qc.py \
  --counts counts_raw.tsv \
  --outdir qc_results/ \
  --min-reads 30 \
  --min-representation 0.8
```

### QC with log2 normalization and custom thresholds
```bash
python scripts/crispr_screen_qc.py \
  --counts counts.tsv \
  --outdir qc_results/ \
  --min-reads 100 \
  --min-representation 0.9 \
  --log2-norm \
  --control-prefix "NonTargeting"
```

### QC on subset of samples
```bash
python scripts/crispr_screen_qc.py \
  --counts counts.tsv \
  --samples T0,T24,T48,T72 \
  --outdir qc_subset/ \
  --n-sgrnas-per-gene 4
```

### QC with custom column names
```bash
python scripts/crispr_screen_qc.py \
  --counts my_counts.tsv \
  --sgrna-col "guide_id" \
  --gene-col "target" \
  --outdir results/ \
  --control-prefix "ctrl"
```

## Parameter Decision Guide

| Parameter | Guidance | Default | Notes |
|-----------|----------|---------|-------|
| `--min-reads` | Threshold below which a guide is considered "not represented"; higher = stricter | 30 | 20-50 typical; use 100+ for very deep screens |
| `--min-representation` | Fraction of library (guides) that must pass --min-reads; lower = more tolerant | 0.8 | 0.7-0.9 typical; <0.7 suggests quality issue |
| `--n-sgrnas-per-gene` | Expected number of sgRNAs per gene for uniformity check | 4 | Affects guide-level CV interpretation |
| `--control-prefix` | String prefix identifying non-targeting control guides | "non" | Common: "NonTargeting", "ctrl", "safe" |
| `--log2-norm` | Apply log2(CPM+1) normalization before correlation/visualization | False | Recommended for variance-stabilization |
| `--sgrna-col`, `--gene-col` | Column headers in input TSV | "sgRNA", "Gene" | Verify column names in your file |

## Failure Modes

| Error | Cause | Solution |
|-------|-------|----------|
| "All samples have 0 reads" | Empty input file or wrong delimiter | Verify TSV is tab-separated; check for BOM or encoding issues |
| "Column 'sgRNA' not found" | Header mismatch | Use --sgrna-col to specify actual header name |
| Gini > 0.5 for all samples | Severe library bias or over-amplification | Check PCR cycle number during library prep; consider earlier passage |
| Correlation < 0.8 between replicates | High variance or batch effect | Check for sample swaps; verify replicate identity; consider removing outlier |
| 0% of library represented | min-reads too high or library too shallow | Lower --min-reads; check sequencing depth (aim for 100-1000x coverage) |
| Missing guide clusters | Technical dropout or selective depletion | Check for guide GC bias; may indicate true selection not QC failure |

## Implementation Details

- **Normalization**: CPM = (reads / total_reads) * 1e6; log2(CPM + 1) for log transform
- **Gini calculation**: O(N log N) sorting followed by O(N) summation
- **Correlation**: scipy.stats Pearson (parametric) and Spearman (rank-based)
- **Plotting**: matplotlib 300 DPI PNG; seaborn for heatmaps
- **Missing data**: treats 0 reads as genuine zeros, not NAs

## Quality Thresholds

| Metric | Pass | Warn | Fail |
|--------|------|------|------|
| Library Representation | >85% | 80-85% | <80% |
| Gini Index | <0.25 | 0.25-0.35 | >0.35 |
| Replicate Correlation (Pearson) | >0.9 | 0.8-0.9 | <0.8 |
| Control Guide CV | <0.3 | 0.3-0.5 | >0.5 |
| Median Reads per Guide | >100 | 50-100 | <50 |
| Missing Guides (%) | <5% | 5-10% | >10% |

## References

- Li, W., et al. (2014). MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biology 15, 554.
- Gini Index: https://en.wikipedia.org/wiki/Gini_coefficient

