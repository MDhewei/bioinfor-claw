# Plot Box Violin - Publication-Quality Distribution Plots

Generate beautiful, publication-ready boxplots, violin plots, and raincloud plots with integrated statistical testing and annotations.

---

## When to Use

Use this skill when you need to:
- Compare distributions of a continuous variable across categorical groups
- Create figures for peer-reviewed publications with statistical annotations
- Visualize distributions with overlaid individual data points
- Perform and annotate statistical comparisons (Mann-Whitney U, t-test, Kruskal-Wallis)
- Generate multiple plot types from the same data (box, violin, combined, raincloud)
- Customize colors, fonts, dimensions for journal submission requirements
- Include detailed statistical results and sample size annotations

## When NOT to Use

Do not use this skill when:
- Plotting time-series data (use line plots instead)
- Creating scatterplots with continuous X and Y variables (use scatter plots)
- Visualizing counts/frequencies across categories (use bar plots)
- Working with very large datasets requiring 2D density visualization (use hexbin/heatmaps)
- Your data has complex hierarchical or multi-level structure (use mixed-effect models)

---

## Input/Output Specification

### Input
Two input layouts are supported; pick whichever matches your data:

**LONG format** (`--value-col` + `--group-col`):
- One numeric column for values, one categorical column for group labels.
- Example:
  ```
  sample    treatment   expression
  S1        Control     2.3
  S2        Control     2.1
  S3        TreatA      5.2
  ```

**WIDE format** (`--wide-cols`):
- Each column is a group; each row is one observation.
- Use when your data already has scores laid out as side-by-side columns
  (e.g. `Control / TreatA / TreatB` columns with numbers in the rows).
- Example:
  ```
  Control   TreatA   TreatB
  2.3       5.2      6.1
  2.1       5.8      NaN
  2.7       4.9      5.4
  ```
- Pass `--wide-cols "Control,TreatA,TreatB"` to select specific columns,
  or `--wide-cols auto` to use every numeric column in the file.
- If neither `--value-col`/`--group-col` nor `--wide-cols` is given, the
  script auto-detects wide format when it sees 2+ numeric columns.
- Sample sizes can differ between groups — NaN cells are dropped per column.

**Format:**
```
Column format: Any format readable by pandas.read_csv()
  - TSV: tab-separated (auto-detected by .tsv extension)
  - CSV: comma-separated
Values: Missing data (NaN/empty) automatically excluded
```

### Output

**Main Output (Required):**
- PNG/SVG/PDF file: Publication-quality plot image
  - Default resolution: 300 DPI
  - File format determined by --output extension
  - High-resolution suitable for print

**Optional Outputs:**
- SVG file: Vector format (--output-svg flag)
  - Perfect for editing in Adobe Illustrator, Inkscape
  - Maintains all styling and annotations
- Statistics TSV (--output-stats flag)
  - Tab-separated table with pairwise test results
  - Columns: Group1, Group2, Statistic, P-value, P-value (corrected), Significant

**Console Output:**
- Summary statistics for each group (n, mean, median, std, min, max)
- Statistical test results (Kruskal-Wallis or pairwise comparisons)
- Significant pairs marked with asterisks (*** p<0.001, ** p<0.01, * p<0.05)

---

## Execution Examples

### Example 1: Basic Boxplot with Statistical Testing

```bash
python scripts/plot_box_violin.py \
  --input data.csv \
  --value-col expression_level \
  --group-col treatment \
  --output boxplot_results.png \
  --plot-type box \
  --stats auto
```

**Data format (data.csv):**
```
sample,treatment,expression_level
S1,Control,2.3
S2,Control,2.1
S3,Treatment,5.2
S4,Treatment,5.8
```

**Output:** Boxplot comparing expression levels between groups with Mann-Whitney U test.

---

### Example 2: Violin Plot with Data Points

```bash
python scripts/plot_box_violin.py \
  --input expression_data.tsv \
  --value-col log2_fc \
  --group-col genotype \
  --output violin_with_points.png \
  --plot-type violin \
  --show-points \
  --point-style jitter \
  --violin-inner quartile \
  --stats all_pairs \
  --title "Gene Expression by Genotype"
```

**Produces:** Violin plot with quartile lines, individual data points with jitter, pairwise statistical comparisons, and significance brackets.

---

### Example 3: Raincloud Plot (Publication Style)

```bash
python scripts/plot_box_violin.py \
  --input measurements.csv \
  --value-col value \
  --group-col condition \
  --output raincloud.png \
  --plot-type raincloud \
  --show-points \
  --stats all_pairs \
  --palette Set2 \
  --group-order control,condition1,condition2 \
  --fig-width 10 \
  --fig-height 6 \
  --output-svg \
  --output-stats stats.tsv
```

**Produces:** Half-violin + strip plot + boxplot with complete statistical analysis and SVG for editing.

---

### Example 4: Horizontal Plot with Custom Colors

```bash
python scripts/plot_box_violin.py \
  --input data.csv \
  --value-col value \
  --group-col treatment \
  --output horizontal_boxplot.png \
  --plot-type box \
  --orientation horizontal \
  --group-colors "#E74C3C,#3498DB,#2ECC71" \
  --show-n \
  --ylim 0,100 \
  --ylabel "Treatment Group"
```

**Produces:** Horizontal boxplot with custom colors, sample sizes labeled, and constrained Y-axis.

---

### Example 5: Combined Violin + Box with Log Scale

```bash
python scripts/plot_box_violin.py \
  --input abundance_data.csv \
  --value-col abundance \
  --group-col taxon \
  --output combined_logscale.png \
  --plot-type both \
  --y-log \
  --show-points \
  --point-alpha 0.5 \
  --stats auto \
  --grid y \
  --spine-style minimal
```

**Produces:** Overlaid violin and box plots with log-scale values, data points, grid, and statistical testing.

---

### Example 6: Wide-Format Data (Each Column = One Group)

```bash
# Auto-detect: every numeric column becomes a group
python scripts/plot_box_violin.py \
  --input scores_wide.csv \
  --output wide_auto.png \
  --plot-type violin \
  --show-points \
  --stats all_pairs

# Explicit: pick which columns to compare
python scripts/plot_box_violin.py \
  --input scores_wide.csv \
  --wide-cols "Control,TreatA,TreatB" \
  --output wide_explicit.png \
  --plot-type raincloud \
  --stats vs_first \
  --ylabel "Score" \
  --output-stats wide_stats.tsv
```

**Data format (scores_wide.csv):**
```
Control,TreatA,TreatB
2.3,5.2,6.1
2.1,5.8,
2.7,4.9,5.4
```

**Produces:** Comparison across columns without manual reshaping. NaN cells are dropped per column, so groups can have different sample sizes.

---

### Example 7: Customized for Journal Submission

```bash
python scripts/plot_box_violin.py \
  --input results.csv \
  --value-col measurement \
  --group-col group \
  --output figure1.png \
  --plot-type box \
  --fig-width 7 \
  --fig-height 5 \
  --dpi 300 \
  --font-family Arial \
  --base-fontsize 11 \
  --axis-label-size 12 \
  --title "Comparison of Measurement by Group" \
  --xlabel "Group" \
  --ylabel "Measurement (units)" \
  --show-n \
  --show-mean-label \
  --stats all_pairs \
  --correction bonferroni \
  --stat-annotation-style bracket \
  --background white \
  --output-svg \
  --output-stats figure1_stats.tsv
```

**Produces:** Publication-ready figure meeting standard requirements with complete statistical table.

---

## Parameter Decision Guide

| Goal | Key Parameters | Example |
|------|---|---|
| **Wide-format data (score columns)** | `--wide-cols auto` or `--wide-cols "A,B,C"` | When each column is a group |
| **Basic comparison (2 groups)** | `--plot-type box` `--stats auto` | `--plot-type box --stats auto` |
| **Show distribution shape** | `--plot-type violin` | `--plot-type violin --violin-inner quartile` |
| **Overlay data points** | `--show-points --point-style` | `--show-points --point-style jitter --point-alpha 0.6` |
| **Publication raincloud** | `--plot-type raincloud --show-points` | `--plot-type raincloud --palette Set2 --output-svg` |
| **Statistical annotations** | `--stats all_pairs --correction bonferroni` | Best for ≤5 groups |
| **Customize colors** | `--palette Set2` OR `--group-colors` | `--group-colors "#FF5733,#33FF57,#3357FF"` |
| **Specific group order** | `--group-order` | `--group-order WT,Mut1,Mut2` |
| **Log scale** | `--y-log` | For abundance/count data with wide range |
| **Horizontal layout** | `--orientation horizontal` | For many/long group names |
| **Multiple formats** | `--output-svg --output-stats` | Always include for publications |
| **Custom figure size** | `--fig-width --fig-height` | `--fig-width 10 --fig-height 6` |
| **Grid reference** | `--grid y` or `--grid both` | Improves readability with many values |
| **Notched boxplot** | `--notch` | Shows 95% CI around median |
| **Compare vs control** | `--stats vs_first` | When first group is control |
| **Dark background** | `--background dark` | For presentations, not printing |

---

## Statistical Testing Details

### Test Selection (--test)
- **auto** (default): Automatically chooses based on data characteristics
  - Normally distributed → t-test
  - Non-normal → Mann-Whitney U (2 groups) or Kruskal-Wallis (>2 groups)
- **mannwhitney**: Non-parametric, robust to outliers (recommended)
- **ttest**: Parametric, assumes normal distribution
- **kruskal**: Non-parametric for >2 groups

### Multiple Testing Correction (--correction)
- **bonferroni** (default): Conservative, controls family-wise error
  - Appropriate for few comparisons (≤10)
  - Adjusted p = p × number_of_tests
- **fdr_bh**: Benjamini-Hochberg FDR control
  - Less conservative, appropriate for many comparisons
  - Balances Type I and II errors
- **none**: No correction (use only for exploratory analysis)

### Significance Thresholds
Default: α = 0.05 (change with --alpha-level)
- *** p < 0.001
- ** p < 0.01
- * p < 0.05
- ns: not significant

---

## Advanced Customization

### Color Schemes
```bash
# Named palettes
--palette Set2        # Default, colorblind-friendly
--palette tab10       # Tableau 10 colors
--palette Set3        # Pastel colors
--palette Pastel1     # Light colors

# Custom hex colors (must match number of groups)
--group-colors "#E74C3C,#3498DB,#2ECC71,#F39C12"
```

### Data Point Styling
```bash
--show-points         # Enable point overlay
--point-style jitter  # Add random noise (default: strip)
--jitter 0.1          # Amount of jitter (0-0.2 recommended)
--point-size 5        # Size in points
--point-alpha 0.6     # Transparency (0-1)
--point-color auto    # auto: match group color (darker)
--point-edge-color black  # Outline color
```

### Violin Plot Options
```bash
--violin-inner box      # Show boxplot inside
--violin-inner quartile # Show quartile lines (default)
--violin-inner point    # Show data points
--violin-scale width    # Width of widest violin = 1 (default)
--violin-scale area     # All violins have equal area
--violin-scale count    # Width proportional to sample size
```

### Statistical Annotation
```bash
--stats auto            # Smart test selection
--stats all_pairs       # All pairwise comparisons
--stats vs_first        # Compare to first group (e.g., control)
--stats vs_last         # Compare to last group
--stat-annotation-style bracket  # Brackets over groups (cleaner)
--stat-annotation-style text     # Text labels only
```

### Axes and Limits
```bash
--ylim 0,100            # Set Y-axis range
--xlim 0,100            # Set X-axis range (horizontal plots)
--y-log                 # Log-scale Y-axis (for wide ranges)
--title "My Title"      # Plot title
--xlabel "Group"        # X-axis label
--ylabel "Value (units)" # Y-axis label
--add-mean-line         # Horizontal line at grand mean
```

### Figure Styling
```bash
--fig-width 10 --fig-height 6  # Dimensions in inches
--dpi 300               # Resolution (300 for print, 72 for screen)
--font-family "Arial"   # Change font
--base-fontsize 11      # Main font size
--axis-label-size 12    # Axis label size
--tick-size 10          # Tick label size
--spine-style minimal   # Remove top/right borders
--grid y                # Add horizontal grid
--background white      # Background color
```

---

## Common Workflows

### Microbiome Abundance Data
```bash
python scripts/plot_box_violin.py \
  --input otus.csv \
  --value-col abundance \
  --group-col treatment \
  --output abundance_plot.png \
  --plot-type raincloud \
  --y-log \
  --show-points \
  --stats vs_first \
  --palette Set2 \
  --output-svg
```

### Gene Expression Comparison
```bash
python scripts/plot_box_violin.py \
  --input expression.csv \
  --value-col log2_expression \
  --group-col genotype \
  --output genotype_comparison.png \
  --plot-type both \
  --show-points \
  --point-style jitter \
  --stats all_pairs \
  --correction bonferroni \
  --title "Gene X Expression by Genotype" \
  --output-stats expression_stats.tsv
```

### Clinical Measurements
```bash
python scripts/plot_box_violin.py \
  --input clinical_data.csv \
  --value-col measurement \
  --group-col disease_status \
  --output clinical_plot.png \
  --plot-type box \
  --show-n \
  --show-mean-label \
  --stats auto \
  --fig-width 6 --fig-height 5 \
  --dpi 300 \
  --font-family Arial
```

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| "Column not found" | Check column names in CSV header match exactly (case-sensitive) |
| Wide CSV but only one group plotted | Pass `--wide-cols "col1,col2,..."` explicitly, or use `--wide-cols auto` to use every numeric column |
| Mixing wide + long flags errors out | Use either `--wide-cols` OR `--value-col`/`--group-col`, not both |
| NaN values causing issues | Script automatically removes rows with missing data; check input file for errors |
| Points overlap too much | Increase --jitter amount (e.g., 0.15) or use --point-alpha < 0.5 |
| Legend/labels cut off | Increase --fig-width or --fig-height; use tight_layout (automatic) |
| Significance brackets overlap | Reduce number of comparisons (use --stats vs_first instead of all_pairs) |
| Statistical test fails | Check group sample sizes; need n≥2 per group; t-test needs n≥3 |
| Colors look wrong | Verify hex colors start with # (e.g., #FF5733); ensure 3 or 6 hex digits |
| Low resolution output | Increase --dpi (default 300 is publication standard) |
| Font not found | Use standard fonts: Arial, Helvetica, Times, Courier, DejaVu Sans |
| Very small p-values | Use --correction none if only 1-2 comparisons; Bonferroni may be too strict |

---

## Technical Notes

- **Statistical Implementation:** Dunn's test implemented in pure NumPy for pairwise comparisons after Kruskal-Wallis
- **Benjamini-Hochberg FDR:** Native NumPy implementation; no external dependencies needed
- **Raincloud:** Half-violin (right by default) + strip plot + embedded boxplot; publication style from Allen et al. 2021
- **Default Palette:** Set2 chosen for colorblind accessibility
- **Grid:** Vertical (Y) grid default for vertical plots to aid value reading
- **Outliers:** Shown by default; disable with --show-fliers=False if needed
- **High-DPI Output:** 300 DPI supports printing at publication standards

---

## File Structure

```
plot-box-violin/
├── scripts/
│   └── plot_box_violin.py     # Main script
├── requirements.txt            # Python dependencies
└── SKILL.md                    # This documentation
```

## Installation

```bash
pip install -r requirements.txt
```

## Quick Start

```bash
# Basic boxplot
python scripts/plot_box_violin.py \
  --input data.csv \
  --value-col value \
  --group-col group \
  --output plot.png

# Publication-ready raincloud
python scripts/plot_box_violin.py \
  --input data.csv \
  --value-col value \
  --group-col group \
  --output plot.png \
  --plot-type raincloud \
  --show-points \
  --stats all_pairs \
  --output-svg \
  --output-stats stats.tsv
```

---

## Citation

If you use this tool in research, please cite:
- matplotlib: Hunter et al. (2007)
- scipy.stats: Virtanen et al. (2020)
- Raincloud plots: Allen et al. (2021) PLoS Biology

---

## Version

1.0 - Initial release with complete statistical testing and visualization options
