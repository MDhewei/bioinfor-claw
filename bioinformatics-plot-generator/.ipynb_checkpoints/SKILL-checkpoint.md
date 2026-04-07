---
name: bioinformatics_plot_generator
description: Use when the user wants a publication-ready bioinformatics plot from tabular data or a numeric matrix, including volcano plots, heatmaps, boxplots, violin plots, scatter plots, and correlation plots.
---

## Purpose
Generate publication-ready bioinformatics plots from tables or matrices using a lightweight, script-based workflow.

## Use when
- the user wants a supported plot from a result table or numeric matrix
- the user wants a figure suitable for publication or presentation
- the required columns are present or can be identified reliably
- the task is a routine plotting task rather than a highly specialized visualization

## Do not use when
- the user wants genome browser tracks or signal plots from BAM, bigWig, or BEDGraph files
- the user wants protein 3D structure visualization
- the user wants a complex multi-panel layout assembled across multiple figures
- the input does not contain the columns needed for the requested plot
- the task requires a specialized plotting ecosystem not supported by this script

## Supported plot types
- volcano
- heatmap
- boxplot
- violin
- scatter
- correlation

## Expected inputs
One of the following:
- a table with effect-size and significance columns for a volcano plot
- a numeric matrix for a heatmap
- a table with a numeric value column and a grouping column for a boxplot or violin plot
- a table with two numeric columns for a scatter or correlation plot

## Expected outputs
- a publication-ready plot file, typically PNG
- a brief summary of the plot type, columns used, and key settings applied

## Procedure
1. Determine the requested plot type from the user’s instruction. If the plot type is not explicitly provided, infer the most likely supported type from the available columns.
2. Validate that the required columns exist and contain usable numeric values where needed.
3. Clean missing, invalid, or non-numeric values relevant to the selected plot.
4. Run `scripts/plot_generator.py` with the appropriate arguments for the selected plot type.
5. Save the output plot.
6. Return the output path and summarize the settings used.

## Plot-specific capabilities

### Volcano
- significance and fold-change threshold lines
- highlight selected upregulated or downregulated features
- label highlighted features
- label top-ranked hits
- separate top up and top down labeling
- optional automatic label adjustment when `adjustText` is available
- configurable maximum number of labels
- configurable axis limits
- configurable point size, transparency, and annotation size

### Heatmap
- numeric matrix plotting
- optional row scaling by z-score
- configurable colormap
- configurable color scale limits
- automatic figure sizing when dimensions are not provided

### Boxplot and violin
- grouped comparison from value and group columns
- optional overlay of individual points
- optional two-group statistical testing
- supported tests:
  - Mann–Whitney U
  - Welch t-test
- p-value and significance annotation for two-group comparisons
- configurable group order

### Scatter and correlation
- scatter plot from two numeric columns
- Pearson or Spearman correlation
- p-value annotation
- sample size annotation
- optional regression line

## Key execution patterns

### Volcano
`python scripts/plot_generator.py --input <input_file> --plot-type volcano --feature-col <feature_col> --x-col <effect_col> --p-col <pvalue_col> --output <output_png>`

Common optional arguments:
- `--fc-cutoff <float>`
- `--p-cutoff <float>`
- `--annotate-top-n <int>`
- `--top-up-n <int>`
- `--top-down-n <int>`
- `--highlight-up <comma_list_or_file>`
- `--highlight-down <comma_list_or_file>`
- `--label-mode <none|top|highlight|top_and_highlight>`
- `--max-labels <int>`
- `--adjust-labels`
- `--annotation-arrow`
- `--xlim <min,max>`
- `--ylim <min,max>`

### Heatmap
`python scripts/plot_generator.py --input <input_file> --plot-type heatmap --index-col <rowname_col> --output <output_png>`

Common optional arguments:
- `--scale-rows`
- `--cmap <colormap>`
- `--vmin <float>`
- `--vmax <float>`

### Boxplot
`python scripts/plot_generator.py --input <input_file> --plot-type boxplot --value-col <value_col> --group-col <group_col> --output <output_png>`

Common optional arguments:
- `--group-order <comma_list>`
- `--show-points`
- `--test-method <mannwhitney|ttest>`

### Violin
`python scripts/plot_generator.py --input <input_file> --plot-type violin --value-col <value_col> --group-col <group_col> --output <output_png>`

Common optional arguments:
- `--group-order <comma_list>`
- `--show-points`
- `--test-method <mannwhitney|ttest>`

### Scatter or correlation
`python scripts/plot_generator.py --input <input_file> --plot-type scatter --x-col <x_col> --y-col <y_col> --output <output_png>`

Common optional arguments:
- `--corr-method <pearson|spearman>`
- `--no-regression`

## Style controls
The script supports customizable figure styling for publication use. These options can be used across plot types when appropriate:

- `--font-family`
- `--base-fontsize`
- `--title-size`
- `--axis-label-size`
- `--tick-size`
- `--legend-size`
- `--annotation-size`
- `--fig-width`
- `--fig-height`
- `--dpi`
- `--point-size`
- `--alpha`

## Conventions
- Prefer explicit user-provided column names when available.
- Prefer absolute file paths if the working directory is uncertain.
- Save plots as PNG unless another format is explicitly requested and supported.
- Keep labeling selective and readable rather than attempting to annotate too many features.
- Use `--adjust-labels` for crowded volcano plots when `adjustText` is installed.

## Dependencies
Expected Python packages:
- pandas
- numpy
- matplotlib
- scipy

Optional but recommended:
- adjustText

## Failure modes
- required columns are missing
- numeric columns contain invalid values
- too many missing values remain after cleaning
- matrix input contains no usable numeric columns
- too many requested labels cause unreadable output
- the Python environment is missing required packages
- `adjustText` is requested implicitly through label adjustment expectations but is not installed

## Interaction model
Users may specify plotting preferences in natural language. When the request is clear and supported, translate those preferences into script arguments.

Examples of supported user preferences include:
- font family and font sizes
- figure width, height, and dpi
- plot title
- significance thresholds
- highlighted features
- label strategy
- group order
- statistical test choice
- correlation method
- whether to show regression lines or sample points

If the user’s request is supported by the plotting script, apply it directly.
If the request is unsupported, use the closest sensible default and briefly state the limitation.