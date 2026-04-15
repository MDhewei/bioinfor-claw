---
name: bioinformatics-plot-generator
description: Route to the correct publication-quality plot sub-skill for volcano plots, heatmaps, box/violin plots, scatter plots, bar charts, MA plots, correlation matrices, and bubble charts from bioinformatics data.
---

# Bioinformatics Plot Generator

## Purpose

This is the **router skill** for all bioinformatics plotting tasks. It selects and delegates to the appropriate sub-skill based on the type of plot requested. Each sub-skill is a fully self-contained, publication-quality plotting tool with 40–70 user-configurable parameters.

## Sub-skills

| Sub-skill | Location | Plot types |
|---|---|---|
| `plot-volcano` | `plot-volcano/` | Volcano plots (DE/CRISPR/GWAS results) |
| `plot-heatmap` | `plot-heatmap/` | Heatmaps with clustering and annotations |
| `plot-box-violin` | `plot-box-violin/` | Box plots, violin plots, raincloud plots |
| `plot-scatter-bar` | `plot-scatter-bar/` | Scatter, bar, MA, correlation matrix, bubble |

## Routing guide

Use this table to pick the correct sub-skill:

| User request signal | Sub-skill to use |
|---|---|
| "volcano plot", "differential expression", "CRISPR screen hits", "-log10 p-value vs fold change" | `plot-volcano` |
| "heatmap", "expression matrix", "gene expression heatmap", "z-score heatmap", "clustered heatmap" | `plot-heatmap` |
| "boxplot", "violin plot", "box and whisker", "raincloud", "distribution comparison", "group comparison" | `plot-box-violin` |
| "scatter plot", "correlation plot", "bar chart", "bar graph", "MA plot", "correlation matrix", "bubble chart", "bubble plot" | `plot-scatter-bar` |

## Use when

- The user wants any of the supported plot types from tabular data or a numeric matrix
- The user wants a figure suitable for publication (300 DPI PNG + SVG)
- The user has a result table (differential expression, CRISPR screen, proteomics, etc.) and wants to visualize it

## Do not use when

- The user wants genome browser tracks or signal plots from BAM/bigWig/BEDGraph files
- The user wants protein 3D structure visualization
- The user wants single-cell UMAP/tSNE trajectory plots requiring Scanpy/Seurat
- The user wants interactive plots (use a Plotly skill instead)

---

## Sub-skill details

### plot-volcano

**Use for:** Volcano plots from differential expression, CRISPR screens, GWAS, proteomics, or any table with a fold-change column and a p-value column.

**Key features:**
- Symmetric or asymmetric fold-change cutoffs (`--fc-cutoff`, `--fc-cutoff-neg`)
- Color by discrete group (up/down/ns) or by continuous column with a colormap
- Three highlight layers: up-regulated, down-regulated, other (custom gene lists or files)
- Top-N auto-annotation with `adjustText` label collision avoidance
- Quadrant counts displayed on plot
- 300 DPI PNG + SVG dual output
- Annotated TSV output table with assigned group per feature

**Script:** `plot-volcano/scripts/plot_volcano.py`

**Minimal run:**
```bash
python plot-volcano/scripts/plot_volcano.py \
  --input results.tsv \
  --feature-col gene \
  --x-col log2FoldChange \
  --p-col padj \
  --output volcano.png
```

---

### plot-heatmap

**Use for:** Heatmaps from any numeric matrix — gene expression, protein abundance, methylation, pathway scores, etc.

**Key features:**
- Hierarchical clustering with dendrogram display (scipy linkage methods: ward, complete, average, single)
- Row and column annotation bars from separate TSV files with custom color palettes
- Variance-based row filtering (keep top N most variable rows)
- Z-score normalization per row or column with optional clipping
- Cell value annotation (show numbers inside cells)
- Auto figure sizing based on matrix dimensions
- Flexible colormap, vmin/vmax, missing value handling

**Script:** `plot-heatmap/scripts/plot_heatmap.py`

**Minimal run:**
```bash
python plot-heatmap/scripts/plot_heatmap.py \
  --input expression_matrix.tsv \
  --index-col gene \
  --output heatmap.png
```

---

### plot-box-violin

**Use for:** Comparing distributions across groups — box plots, violin plots, combined box+violin, or raincloud plots.

**Key features:**
- Four plot types: `box`, `violin`, `both` (violin with inner box), `raincloud`
- Jittered individual points with customizable size and alpha
- Multi-group statistical testing: auto, all_pairs, vs_first, vs_last
- Tests: Mann–Whitney U, Welch t-test, Dunn's test (pure numpy, no scipy dependency)
- Multiple testing correction: Bonferroni or Benjamini–Hochberg FDR
- Significance bracket annotations drawn above each pair
- Horizontal orientation option
- Custom group ordering and color palettes

**Script:** `plot-box-violin/scripts/plot_box_violin.py`

**Minimal run:**
```bash
python plot-box-violin/scripts/plot_box_violin.py \
  --input data.tsv \
  --value-col expression \
  --group-col condition \
  --plot-type violin \
  --output violin.png
```

---

### plot-scatter-bar

**Use for:** Five plot types in one script — scatter plots, bar charts, MA plots, correlation matrices, and bubble charts.

**Key features:**

**Scatter:**
- Linear regression line with 95% CI band (bootstrap)
- Pearson/Spearman correlation annotation
- Optional marginal histograms
- Highlight gene list with separate color/label

**Bar:**
- Grouped or stacked bar charts
- Error bars: SEM, SD, or 95% CI
- Value labels on bars
- Horizontal orientation

**MA (ratio vs. mean):**
- M-A plot (log-ratio vs. average intensity)
- LOESS smoothing line (tricubic kernel, pure numpy)
- Significant feature highlighting

**Correlation matrix:**
- Pairwise correlation heatmap (Pearson or Spearman)
- Hierarchical clustering of samples
- Correlation values displayed in cells

**Bubble:**
- Bubble chart with x/y/size/color columns
- Log-scale size normalization
- Separate size and color legends

**Script:** `plot-scatter-bar/scripts/plot_scatter_bar.py`

**Minimal runs:**
```bash
# Scatter
python plot-scatter-bar/scripts/plot_scatter_bar.py \
  --input data.tsv --plot-type scatter \
  --x-col sample1 --y-col sample2 --output scatter.png

# Bar chart
python plot-scatter-bar/scripts/plot_scatter_bar.py \
  --input counts.tsv --plot-type bar \
  --value-col count --group-col condition --output bar.png

# MA plot
python plot-scatter-bar/scripts/plot_scatter_bar.py \
  --input de_results.tsv --plot-type ma \
  --mean-col baseMean --ratio-col log2FC --output ma.png

# Correlation matrix
python plot-scatter-bar/scripts/plot_scatter_bar.py \
  --input expr_matrix.tsv --plot-type corrmat \
  --index-col gene --output corrmat.png

# Bubble chart
python plot-scatter-bar/scripts/plot_scatter_bar.py \
  --input enrichment.tsv --plot-type bubble \
  --x-col NES --y-col pathway \
  --size-col gene_count --color-col padj --output bubble.png
```

---

## Common style parameters (all sub-skills)

| Parameter | Default | Description |
|---|---|---|
| `--fig-width` | auto | Figure width in inches |
| `--fig-height` | auto | Figure height in inches |
| `--dpi` | 300 | Output resolution |
| `--font-family` | `sans-serif` | Font family |
| `--base-fontsize` | 11 | Base font size (pt) |
| `--title` | — | Plot title |
| `--xlabel` | auto | X-axis label |
| `--ylabel` | auto | Y-axis label |
| `--output` | required | Output PNG path |
| `--output-svg` | off | Also save SVG alongside PNG |

## Procedure for routing

1. Identify the plot type from the user's request.
2. Read the corresponding sub-skill SKILL.md (e.g., `plot-volcano/SKILL.md`) for the full parameter reference.
3. Identify required input columns from the user's data or description.
4. Construct the command with appropriate parameters from the parameter decision guide in the sub-skill SKILL.md.
5. Run the script and return the output path.
6. If the user requests customization beyond the defaults, consult the sub-skill's parameter decision guide and full argument list.
