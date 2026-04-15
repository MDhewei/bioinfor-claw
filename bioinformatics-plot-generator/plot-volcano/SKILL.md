# plot-volcano

Generate publication-quality volcano plots from differential expression analysis results with extensive customization of colors, thresholds, labels, and styling. Powered by matplotlib for lightweight, reproducible visualization.

## Overview

Volcano plots are standard visualizations in genomics for displaying the relationship between fold change and statistical significance in differential expression analysis. This tool transforms raw DE result tables into polished, publication-ready plots suitable for journals, presentations, and supplementary materials.

**Key features:**
- Auto-detection of common column names (gene, symbol, ID, etc.)
- Flexible significance and fold-change thresholds with asymmetric cutoffs
- Multiple labeling strategies: top genes by score, highlighted gene sets, or both
- Continuous coloring by expression metrics (baseMean, expression levels, etc.)
- Publication-quality styling with extensive customization
- Smart label positioning with adjustText integration
- Quadrant annotation showing feature counts
- Optional output table with categorization for downstream analysis

## When to Use

**Use this tool when you:**
- Have raw or preprocessed differential expression results (DESeq2, edgeR, limma, etc.)
- Need quick exploratory volcano plots during analysis
- Want publication-quality visualizations with minimal post-processing
- Need to highlight specific gene sets or top candidates
- Want to compare results across multiple conditions/tissues/cell types
- Need to integrate volcano plots into analysis pipelines or reports

**Do not use this tool when:**
- You need volcano plots for non-genomic fold-change vs. p-value data (consider adapting the tool)
- Your input data is heavily filtered or curated beyond standard thresholds
- You require interactive plots (consider plotly/ggplot2 alternatives)
- You need volcano plots integrated directly into genome browsers or complex multi-panel figures

## Input/Output Specification

### Input Format

**Input file:** TSV, CSV, or whitespace-delimited table with DE results

**Required columns:**
- **log2FC column:** Raw or normalized log2 fold change values (typically -5 to +5 range)
- **p-value column:** Raw p-values (0-1 range) OR already -log10 transformed (0-300+ range)
  - If using adjusted p-values (p.adj, padj, FDR): specify via `--p-col`
  - If pre-transformed: use `--p-is-log` flag

**Optional columns:**
- **Feature names:** gene symbols, ENSEMBL IDs, etc. Auto-detected from headers: `gene`, `Gene`, `symbol`, `Symbol`, `id`, `ID`, `feature`, `Feature`, `name`, `Name`
- **Color column:** Numeric values for continuous coloring (e.g., `baseMean`, `log2CPM`, expression levels)

**Example input (TSV):**
```
gene      logFC    p.adj       baseMean
BRCA1     -2.5     0.0001      5000
TP53       3.2     0.0001      8000
GAPDH      0.1     0.8         50000
...
```

### Output Files

**Default output:** PNG image (300 DPI)
- Configurable dimensions (default 8x6 inches)
- Supports PNG, PDF, SVG via file extension

**Optional outputs:**
- `--output-svg`: Additional SVG version (vector format for editing)
- `--output-table`: TSV table with categorization (up/down/ns) and label assignments

## Execution Examples

### Example 1: Basic Volcano Plot from DESeq2 Results
```bash
python plot_volcano.py \
  --input deseq_results.tsv \
  --x-col log2FoldChange \
  --p-col padj \
  --output volcano_basic.png
```
Auto-detects gene symbols, uses default thresholds (FC>1, p<0.05).

### Example 2: Custom Thresholds with Highlighted Gene Sets
```bash
python plot_volcano.py \
  --input deseq_results.tsv \
  --x-col log2FoldChange \
  --p-col padj \
  --fc-cutoff 1.5 \
  --p-cutoff 0.01 \
  --highlight-up my_interest_genes.txt \
  --highlight-down tumor_suppressors.txt \
  --output volcano_highlighted.png
```
Uses stricter thresholds; highlights specific gene groups in distinct colors.

### Example 3: Coloring by Expression Level
```bash
python plot_volcano.py \
  --input deseq_results.tsv \
  --x-col log2FoldChange \
  --p-col padj \
  --color-col baseMean \
  --colormap plasma \
  --show-colorbar \
  --output volcano_colored.png
```
Points colored by expression level; includes colorbar legend.

### Example 4: Raw p-values with Asymmetric Thresholds
```bash
python plot_volcano.py \
  --input deseq_results.tsv \
  --x-col logFC \
  --p-col pvalue \
  --fc-cutoff 1.0 \
  --fc-cutoff-neg 0.8 \
  --annotate-top-n 15 \
  --max-labels 25 \
  --output volcano_asymmetric.png
```
Down-regulation threshold (0.8) stricter than up-regulation (1.0); labels top 15 genes overall.

### Example 5: Publication-Ready with Custom Styling
```bash
python plot_volcano.py \
  --input de_results.tsv \
  --x-col log2FC \
  --p-col p.adj \
  --label-mode top_and_highlight \
  --top-up-n 5 \
  --top-down-n 5 \
  --highlight-up key_genes.txt \
  --font-family Helvetica \
  --title "Expression Changes in Treatment vs. Control" \
  --xlabel "log2(Fold Change)" \
  --ylabel "-log10(Adjusted P-value)" \
  --fig-width 10 --fig-height 7 \
  --output-svg \
  --output-table annotated_results.tsv \
  --output volcano_publication.png
```
Publication quality: custom fonts, titles, legends, and output formats.

## Parameter Decision Guide

| Use Case | Key Parameters | Notes |
|----------|---|---|
| **Quick exploration** | `--input`, `--x-col`, `--p-col` | Uses defaults for all other settings |
| **Stricter significance** | `--p-cutoff 0.01`, `--fc-cutoff 1.5` | More conservative feature selection |
| **Raw p-values** | `--p-col pvalue` (leave `--p-is-log` off) | Automatically transforms to -log10 |
| **Pre-transformed p-values** | `--p-col log10_p`, `--p-is-log` | Skips transformation |
| **Highlight key genes** | `--highlight-up genes.txt`, `--highlight-down genes.txt` | Separate file per group or comma-separated |
| **Top features only** | `--label-mode top`, `--top-up-n 10`, `--top-down-n 10` | Labels only highest-scoring genes |
| **Expression coloring** | `--color-col baseMean`, `--colormap viridis`, `--show-colorbar` | Color intensity shows expression level |
| **Asymmetric cutoffs** | `--fc-cutoff 1.0`, `--fc-cutoff-neg 0.8` | Different thresholds for up vs. down |
| **Smart label placement** | `--adjust-labels` | Requires adjustText; prevents overlaps |
| **Multiple formats** | `--output volcano.png`, `--output-svg` | PNG + SVG for flexibility |
| **Downstream analysis** | `--output-table results_annotated.tsv` | Re-use categorization in other tools |
| **High-quality publication** | `--fig-width 10 --fig-height 8 --dpi 300`, `--font-family Arial` | Large, print-ready output |
| **Hide top/right axis frame** | `--hide-spines top,right` (or shortcut `--hide-spines top-right`) | Open "L-shape" axis style; common in publications. Use `--hide-spines all` to remove the full box, or list any combo of `top,right,bottom,left`. |

### Labeling Strategy Selection

| Strategy | Use When | Result |
|---|---|---|
| `none` | Too crowded; only need the plot | No labels; clean background |
| `top` | Want to see most significant genes | Labels top N genes by (−log10p) × |logFC| |
| `highlight` | Have specific gene sets to focus on | Only labeled genes highlighted; others gray |
| `top_and_highlight` | Want both general interest + specific groups | Top genes + highlighted sets both labeled |

### Color Customization

**Categorical coloring (default):**
- `--color-up`: Upregulated (default orange #D55E00)
- `--color-down`: Downregulated (default blue #0072B2)
- `--color-ns`: Non-significant (default gray #BDBDBD)
- Highlighted groups override categorical colors

**Continuous coloring (when `--color-col` is set):**
- `--colormap`: Choose from matplotlib colormaps (viridis, plasma, coolwarm, etc.)
- Points colored by numeric column values (e.g., expression level)
- `--show-colorbar`: Display colorbar legend

### Threshold Line Customization

Fine-tune appearance of FC and p-value threshold lines:
- `--threshold-linestyle`: dashed (default), dotted, solid
- `--threshold-linewidth`: 0.9 (default)
- `--threshold-color`: gray (default)
- `--show-fc-line`, `--show-p-line`: Toggle visibility

## Advanced Options

### Scoring and Ranking

Features ranked by **combined score**: (−log10 p) × |log2FC|

This penalizes features with:
- Small fold changes (even if significant)
- High p-values (even if large effect)
- Favors features with both strong effect and high significance

### Feature Auto-Detection

Attempts to auto-detect feature/gene name column from:
1. Exact matches: gene, Gene, symbol, Symbol, id, ID, feature, Feature, name, Name
2. Dataframe index (if named)
3. Falls back to generic Feature_0, Feature_1, ... if not found

Override with `--feature-col` if auto-detection fails.

### Label Conflict Resolution

When too many features qualify for labeling:
- `--max-labels 20` limits total labels
- Enforces limit by keeping highest-score features
- If crowded, try `--adjust-labels` (with adjustText library) for smart positioning

### Missing Value Handling

- NaN values in FC or p-value columns are excluded automatically
- Color column NaNs are handled gracefully (ignored if using continuous colors)
- Feature names converted to strings; missing values become string "nan"

## Performance Notes

- **Speed:** Typical runtime < 2 seconds for 20,000 features
- **Memory:** ~500MB for 100,000 features
- **Output file size:** PNG ~500KB–2MB; SVG ~10–50MB (uncompressed)

For very large datasets (>100K features):
- Consider pre-filtering to significant features
- Use `--label-mode none` for faster rendering
- Export as PNG rather than SVG

## Troubleshooting

**Error: Column 'X' not found**
- Check column name spelling and case sensitivity
- List available columns with: `head -n 2 your_file.tsv`

**No labels appearing**
- Check `--label-mode none` is not set
- Ensure features meet threshold criteria
- Try increasing `--top-up-n` and `--top-down-n`

**Labels overlapping**
- Use `--adjust-labels` (requires adjustText: `pip install adjustText`)
- Increase figure size: `--fig-width 12 --fig-height 9`
- Reduce `--max-labels`

**Colors not showing for highlights**
- Verify gene names in highlight file exactly match input data
- Check that highlighted genes exceed significance thresholds
- Try listing genes as comma-separated: `--highlight-up "GENE1,GENE2,GENE3"`

**Output file missing or blank**
- Verify output directory exists or use absolute path
- Check disk space
- Ensure matplotlib backend configured (usually automatic)

## Implementation Details

**Dependencies:**
- numpy: Numerical operations and statistics
- pandas: Data I/O and manipulation
- matplotlib: Plotting and styling
- adjustText (optional): Smart label positioning

**Zero scipy dependency:** All statistical transformations (log transformations, categorization) implemented in pure numpy.

**Encoding:** Input files assumed UTF-8. Non-ASCII characters in gene names preserved when possible.

**Reproducibility:** All random elements removed; output deterministic given identical inputs.

---

*Built for genomics researchers and bioinformaticians. For questions or feature requests, please open an issue.*
