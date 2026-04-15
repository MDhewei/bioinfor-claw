# plot-scatter-bar

**Version:** 1.0
**Author:** Bioinformatics Team
**Updated:** 2026-04-08

Publication-quality scatter, bar, MA, correlation matrix, and bubble plots for bioinformatics data visualization. Built with matplotlib, numpy, and pandas.

## When to Use

- **Scatter plots**: Compare two continuous variables; detect correlations and outliers
- **Bar plots**: Display means/aggregates for categorical groups; compare across conditions
- **MA plots**: Differential expression analysis; identify significantly regulated genes
- **Correlation matrices**: Visualize pairwise correlations among samples/features
- **Bubble charts**: Show three or more dimensions (x, y, size, color)

## Do Not Use When

- Data is better represented as line plots (time series, continuous trajectories)
- Displaying very large sample sizes (>10,000 points) without aggregation
- Creating 3D plots or highly specialized domain-specific visualizations
- Handling non-numeric data that cannot be encoded as position, size, or color

---

## Input/Output Specification

### Input Format
- **File types**: CSV, TSV, or plain-text table (auto-detected)
- **Columns**: Must contain numeric columns for plotting; categorical columns optional
- **Missing values**: Automatically removed (NaN, blank cells)
- **Row labels**: Optional index column for sample names

### Output Formats
- **Primary output** (`--output`): PNG (default), SVG, or PDF
- **Optional SVG** (`--output-svg`): Save additional vector format for editing
- **Optional table** (`--output-table`): Save processed/annotated data as TSV

### Plot Types

#### SCATTER
Compare two continuous variables with optional regression analysis.

**Input:**
- Two numeric columns: `--x-col` and `--y-col` (required)
- Optional: categorical/numeric column for color encoding (`--color-col`)
- Optional: numeric column to scale point size (`--size-col`)
- Optional: column with point labels (`--label-col`)

**Output:** Scatter plot with regression line, confidence interval, marginal histograms (optional)

---

#### BAR
Display aggregated values for categorical groups with error bars.

**Input:**
- Category column: `--x-col` (required)
- Value column: `--y-col` (required; auto-aggregated as mean if raw data provided)
- Optional: pre-computed error column (`--error-col`)
- Optional: grouping/coloring column (`--color-col`)

**Output:** Bar plot with optional grouped/stacked bars, error bars, overlaid points, statistical annotations

---

#### MA
Differential expression visualization: log2 fold-change (M) vs mean expression (A).

**Input:**
- Mean expression column: `--x-col` (log2 or raw; auto-log if max >100)
- log2 fold-change column: `--y-col` (required)
- Optional: p-value column for significance coloring (`--p-col`)
- Optional: gene/feature names (`--feature-col`)

**Output:** MA plot with up/down/not-significant coloring, LOESS smoothing line, fold-change thresholds, top DE gene annotations

---

#### CORRMAT
Pairwise correlation matrix between numeric columns.

**Input:**
- Numeric columns: All numeric data in input file (or subset if specified)
- Optional: row label column for sample names (`--index-col`)

**Output:** Heatmap with dendrograms (if clustering enabled), correlation values, optional masking of upper triangle

---

#### BUBBLE
Multi-dimensional plot: x, y, size, and optional color encoding.

**Input:**
- X and Y axes: `--x-col`, `--y-col` (required)
- Bubble size: `--size-col` (required)
- Optional: color encoding (`--color-col`; categorical or numeric)
- Optional: labels (`--label-col`)

**Output:** Bubble plot with size and color legends, optional label annotations

---

## Full Parameter Specification

### Common Input/Output Parameters
```
--input FILE                  Input CSV/TSV table (required)
--plot-type TYPE              One of: scatter, bar, ma, corrmat, bubble (required)
--output FILE                 Output file path: PNG/SVG/PDF (required)
--output-svg FILE             Additionally save as SVG
--output-table FILE           Save processed data as TSV
```

### Scatter Plot Parameters
```
--x-col COL                   X-axis column name (required)
--y-col COL                   Y-axis column name (required)
--color-col COL               Column to color points by (categorical → discrete, numeric → colormap)
--size-col COL                Column to scale point size by
--label-col COL               Column containing point labels
--label-top-n N               Label top N points by distance from origin (default: 0)
--highlight FEATURES          Comma-separated feature names to highlight; can be a file path
--highlight-color HEX         Color for highlighted points (default: #D55E00)
--corr-method METHOD          Correlation method: pearson, spearman, none (default: pearson)
--show-regression BOOL        Show regression line (default: True)
--regression-ci BOOL          Show 95% confidence interval band (default: True)
--regression-color HEX        Regression line color (default: #2166AC)
--point-size FLOAT            Point size in scatter units (default: 25)
--alpha FLOAT                 Transparency: 0–1 (default: 0.7)
--color-ns HEX                Color for non-significant/uncolored points (default: #999999)
--colormap NAME               Colormap name: viridis, plasma, cool, hot, etc. (default: viridis)
--show-colorbar BOOL          Show colorbar for numeric color-col (default: True)
--identity-line BOOL          Draw y=x diagonal reference line (default: False)
--xlabel STRING               X-axis label (default: column name)
--ylabel STRING               Y-axis label (default: column name)
--title STRING                Plot title
--xlim MIN MAX                X-axis limits (space-separated floats)
--ylim MIN MAX                Y-axis limits (space-separated floats)
--log-x BOOL                  Log10-scale x-axis (default: False)
--log-y BOOL                  Log10-scale y-axis (default: False)
--marginal-hist BOOL          Show marginal histograms on sides (default: False)
```

### Bar Plot Parameters
```
--x-col COL                   Category column (required)
--y-col COL                   Value column (required; auto-aggregated as mean if raw data)
--error-col COL               Pre-computed error bar column (SEM, SD, or CI95)
--error-type TYPE             Error type if computing from data: sem, sd, ci95 (default: sem)
--color-col COL               Column to group/color bars by (enables grouped bars)
--group-order CATS            Comma-separated x-axis category order
--color-order GROUPS          Comma-separated color group order
--palette NAME                Color palette: Set2, tab10, Pastel1 (default: Set2)
--bar-colors HEXES            Comma-separated hex colors (overrides palette)
--orientation ORIENT          vertical or horizontal (default: vertical)
--bar-width FLOAT             Bar width: 0–1 (default: 0.7)
--show-points BOOL            Overlay individual data points (default: False)
--point-jitter FLOAT          Point jitter: amount of random horizontal noise (default: 0.05)
--capsize FLOAT               Error bar cap width (default: 4)
--show-values BOOL            Show bar height values as text (default: False)
--value-fmt FORMAT            Format string for values: .2f, .1e, etc. (default: .2f)
--sort-by SORT                Sort bars: value, name, or none (default: none)
--stacked BOOL                Stacked bar chart (requires color-col; default: False)
--ylabel STRING               Y-axis label (default: column name)
--ylim MIN MAX                Y-axis limits
--xlim MIN MAX                X-axis limits
--title STRING                Plot title
--stats TYPE                  Statistical annotations: none or pairwise (default: none)
```

### MA Plot Parameters
```
--x-col COL                   Mean expression column (log or raw; auto-log if max >100) (required)
--y-col COL                   log2 fold-change column (required)
--p-col COL                   P-value column for significance coloring
--feature-col COL             Gene/feature name column for annotations
--p-cutoff FLOAT              Adjusted p-value significance threshold (default: 0.05)
--fc-cutoff FLOAT             log2 fold-change magnitude threshold (default: 1.0)
--color-up HEX                Color for upregulated genes (default: #D55E00)
--color-down HEX              Color for downregulated genes (default: #0072B2)
--color-ns HEX                Color for non-significant genes (default: #BDBDBD)
--annotate-top-n N            Label top N DE genes by p-value × |FC| (default: 10)
--loess-line BOOL             Draw LOESS smoothing line (default: True)
--loess-color HEX             LOESS line color (default: black)
--xlabel STRING               X-axis label (default: "Mean expression (log2)")
--ylabel STRING               Y-axis label (default: "log2 Fold Change")
--title STRING                Plot title
--point-size FLOAT            Point size (default: 12)
--alpha FLOAT                 Transparency: 0–1 (default: 0.5)
```

### Correlation Matrix Parameters
```
--index-col COL               Row label column (optional; uses column names if not provided)
--method METHOD               Correlation method: pearson or spearman (default: pearson)
--cluster BOOL                Cluster rows/columns by correlation distance (default: True)
--cmap NAME                   Colormap: RdBu_r, coolwarm, bwr, etc. (default: RdBu_r)
--vmin FLOAT                  Colormap minimum value (default: -1)
--vmax FLOAT                  Colormap maximum value (default: 1)
--show-values BOOL            Show correlation coefficients in cells (default: True)
--value-fmt FORMAT            Format for correlation values: .2f, .3f, etc. (default: .2f)
--value-fontsize FLOAT        Font size for cell values (default: 8)
--mask-upper BOOL             Mask upper triangle (default: False)
--label-rotation FLOAT        Rotation angle for x-axis labels in degrees (default: 45)
--title STRING                Plot title
```

### Bubble Plot Parameters
```
--x-col COL                   X-axis column (required)
--y-col COL                   Y-axis column (required)
--size-col COL                Bubble size column (required)
--color-col COL               Color encoding column (categorical or numeric, optional)
--label-col COL               Column for bubble labels
--label-all BOOL              Label all bubbles (default: False)
--min-size FLOAT              Minimum bubble area (default: 50)
--max-size FLOAT              Maximum bubble area (default: 2000)
--alpha FLOAT                 Transparency: 0–1 (default: 0.75)
--palette NAME                Color palette for categorical color-col: tab10, Set2, etc. (default: tab10)
--colormap NAME               Colormap for numeric color-col (default: viridis)
--show-size-legend BOOL       Show bubble size legend (default: True)
--show-color-legend BOOL      Show color legend (default: True)
--xlabel STRING               X-axis label (default: column name)
--ylabel STRING               Y-axis label (default: column name)
--title STRING                Plot title
```

### Style Parameters (All Plot Types)
```
--font-family FONT            Font family: Arial, Helvetica, DejaVu, etc. (default: Arial)
--base-fontsize FLOAT         Base font size in points (default: 11)
--axis-label-size FLOAT       Axis label font size (default: 12)
--tick-size FLOAT             Tick label font size (default: 10)
--legend-size FLOAT           Legend font size (default: 10)
--title-size FLOAT            Title font size (default: 13)
--fig-width FLOAT             Figure width in inches (default: 8)
--fig-height FLOAT            Figure height in inches (default: 6)
--dpi INT                     Output resolution in DPI (default: 300)
--spine-style STYLE           Border visibility: all or minimal (default: minimal)
--grid GRID                   Grid lines: none, major, or both (default: none)
--grid-alpha FLOAT            Grid transparency: 0–1 (default: 0.3)
```

---

## Execution Examples

### Example 1: Scatter Plot with Regression
```bash
python scripts/plot_scatter_bar.py \
  expression.csv \
  --plot-type scatter \
  --x-col expression_log2 \
  --y-col expression_log10 \
  --show-regression True \
  --regression-ci True \
  --title "Expression Correlation" \
  --xlabel "log2(Expression)" \
  --ylabel "log10(Expression)" \
  --output scatter_example.png
```

### Example 2: Bar Plot with Grouped Bars
```bash
python scripts/plot_scatter_bar.py \
  conditions.csv \
  --plot-type bar \
  --x-col genotype \
  --y-col expression_level \
  --color-col treatment \
  --error-col sem \
  --show-values True \
  --title "Expression by Genotype and Treatment" \
  --ylabel "Expression Level" \
  --palette Set2 \
  --output bar_grouped.png
```

### Example 3: MA Plot for Differential Expression
```bash
python scripts/plot_scatter_bar.py \
  deseq_results.tsv \
  --plot-type ma \
  --x-col baseMean \
  --y-col log2FoldChange \
  --p-col padj \
  --feature-col gene_name \
  --p-cutoff 0.05 \
  --fc-cutoff 1.0 \
  --annotate-top-n 15 \
  --loess-line True \
  --title "MA Plot: Treatment vs Control" \
  --output ma_plot.png
```

### Example 4: Correlation Matrix with Clustering
```bash
python scripts/plot_scatter_bar.py \
  sample_matrix.csv \
  --plot-type corrmat \
  --index-col sample_id \
  --method pearson \
  --cluster True \
  --show-values True \
  --cmap RdBu_r \
  --title "Sample-to-Sample Correlation" \
  --output corrmat.png
```

### Example 5: Bubble Plot with Color Encoding
```bash
python scripts/plot_scatter_bar.py \
  genes.csv \
  --plot-type bubble \
  --x-col log2fc \
  --y-col -log10_pvalue \
  --size-col expression_level \
  --color-col pathway \
  --min-size 50 \
  --max-size 1000 \
  --label-top-n 10 \
  --title "Gene Significance and Expression" \
  --output bubble_genes.png
```

### Example 6: Scatter with Highlights and Marginal Histograms
```bash
python scripts/plot_scatter_bar.py \
  expression_pairs.csv \
  --plot-type scatter \
  --x-col gene1_expr \
  --y-col gene2_expr \
  --highlight genes_of_interest.txt \
  --highlight-color "#FF0000" \
  --marginal-hist True \
  --show-regression True \
  --title "Gene-Gene Expression Correlation" \
  --output scatter_highlighted.png
```

---

## Parameter Decision Guide

| Scenario | Plot Type | Key Parameters |
|----------|-----------|-----------------|
| Compare two numeric variables | **scatter** | `--x-col`, `--y-col`, `--color-col` (optional), `--show-regression` |
| Show means ± SEM by category | **bar** | `--x-col` (category), `--y-col` (value), `--error-type sem` |
| Grouped bars across conditions | **bar** | `--color-col` (group), `--group-order`, `--error-col` |
| Stacked bar chart | **bar** | `--color-col`, `--stacked True` |
| Visualize DE results | **ma** | `--x-col` (baseMean), `--y-col` (log2FC), `--p-col`, `--feature-col` |
| Sample similarity matrix | **corrmat** | All numeric columns, `--cluster True`, `--method pearson` |
| Three+ dimensional data | **bubble** | `--x-col`, `--y-col`, `--size-col`, `--color-col` (optional) |
| Highlight outliers/specific points | **scatter** | `--highlight FEATURES`, `--highlight-color`, `--label-col` |
| Publication-quality vector output | Any | `--output-svg FILE.svg`, `--dpi 300` |
| Log-scale axes | **scatter** | `--log-x True` and/or `--log-y True` |

---

## Implementation Details

- **LOESS smoothing** (MA plot): Tricubic kernel weighted least squares with frac=0.3
- **Regression confidence interval** (scatter): 95% analytical band computed with t-distribution
- **Clustering** (corrmat): Scipy hierarchical clustering if available; falls back to sorting
- **Size normalization** (bubble): Min-max scaling using Dunn normalization
- **Color encoding**: Automatic colormap selection for numeric; discrete colors for categorical
- **Missing values**: Automatically removed from analysis
- **Summary output**: Each plot prints summary statistics (n points, key parameters)

---

## Dependencies
- matplotlib >= 3.5
- numpy >= 1.22
- pandas >= 1.4
- scipy >= 1.8 (optional for clustering, PCA fallback available)

## Notes
- For large datasets (>10,000 points), consider aggregation or alpha transparency to avoid overplotting
- Regression line requires at least 3 points; CI band requires at least 5 points
- LOESS smoothing is computationally intensive for >5,000 points
- Colormap choices: viridis, plasma, inferno, magma, cividis (perceptually uniform); RdBu_r, coolwarm, bwr (diverging)
