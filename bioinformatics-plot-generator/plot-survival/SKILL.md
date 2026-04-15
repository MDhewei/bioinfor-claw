---
name: plot-survival
description: Generate publication-quality Kaplan-Meier survival curves from time-to-event data. Supports multiple groups, log-rank test, Wilcoxon test, median survival annotation, confidence intervals, at-risk tables, and flexible stratification by continuous or categorical variables.
---

# Plot Survival

## Purpose
Create professional Kaplan-Meier survival curves for time-to-event data with comprehensive statistical testing and annotations. This skill computes the Kaplan-Meier estimator from raw survival data, performs log-rank and Wilcoxon tests, generates publication-quality plots with confidence intervals, median annotations, at-risk tables, and optional output in PNG and SVG formats.

## Use when / Not when

### Use when
- You need to **visualize time-to-event data** (survival, recurrence-free survival, disease-free survival, etc.)
- You want to **compare survival between groups** (treatment vs. control, responder vs. non-responder)
- You need to **test for significant group differences** (log-rank or Wilcoxon p-values)
- You're **preparing figures for a publication** (high-quality, customizable aesthetics)
- You want to **stratify by continuous variables** (e.g., high vs. low expression, split at median)
- You need **confidence intervals** and **at-risk tables** (clinical/biomedical standard)
- You're **reporting median survival** and survival probabilities at specific time points

### Not when
- Your data has **competing risks** (use competing risks packages like cmprsk)
- You have **more than 10-15 groups** (plot becomes unreadable; consider binning)
- You need **Cox regression** or **adjusted curves** (use survival R package or lifelines Python)
- You have **hierarchical data** with repeated events per individual (use frailty models)
- Data is **left-censored** or **interval-censored** (only right-censoring is supported)
- You need **permutation-based p-values** (script uses parametric chi-squared approximation)

## Expected inputs / outputs

### Inputs
- **CSV/TSV file** with columns:
  - Time-to-event (required): numeric, column name specified by `--time-col`
  - Event status (required): 1 = event occurred, 0 = censored
  - Group (required for multi-group KM): categorical variable name via `--group-col`
  - (Optional) Continuous variable for median stratification (e.g., gene expression)

### Outputs
- **PNG file** (default): High-resolution Kaplan-Meier curve with all annotations
- **SVG file** (if `--output-svg`): Scalable vector graphics for editing in Illustrator/Inkscape
- **survival_stats.tsv**: Table with group, n, n_events, median_survival, CI_lower, CI_upper, p_value
- **Console output**: Summary statistics and test results

## Procedure

### Step 1: Load and validate data
1. Read CSV/TSV into pandas DataFrame
2. Check for required columns: `--time-col`, `--event-col`, `--group-col`
3. Validate data types: time should be numeric, event should be 0/1
4. Remove rows with missing values in key columns
5. Handle median stratification if `--split-by-median` is set

### Step 2: Kaplan-Meier estimation (per group)
For each group:
1. Sort records by time-to-event
2. Identify unique event times
3. At each event time t, compute:
   - d_i = number of events at time t
   - n_i = number at risk at time t
   - S(t) = S(t-1) × (1 - d_i / n_i)  [product of conditional survival]
4. Compute Greenwood's variance:
   - Var(S(t)) = S(t)^2 × sum(d_i / (n_i × (n_i - d_i)))
5. Compute 95% CI using log-log transformation:
   - α = sqrt(Var(log(-log(S))))
   - CI = exp(-exp(log(-log(S)) ± 1.96 × α))
6. Identify median survival time (time at which S(t) = 0.5)

### Step 3: Statistical testing
**Log-rank test:**
1. Compute observed vs. expected events across all event times
2. O_i = observed events in group i at each time
3. E_i = expected events = (d_i × n_i) / N, where N = total at risk
4. Chi-squared statistic: χ² = (O_i - E_i)^2 / E_i, df = n_groups - 1
5. Compute p-value from chi-squared distribution (Wilson-Hilferty normal approx for small sample)

**Wilcoxon (Peto-Wilcoxon) test:**
1. Weight each comparison by number at risk: w_i = n_i (gives more weight to early events)
2. Otherwise same as log-rank test with weighted statistics

### Step 4: Plot generation (matplotlib)
1. Create figure with specified dimensions and DPI
2. For each group, plot step function:
   - x = event times, y = S(t)
   - Use "where=post" for step function
   - Color from palette (Set1, tab10, or custom)
   - Linewidth customizable
3. Add shaded confidence interval bands (fill_between with alpha transparency)
4. Add censoring tick marks at censored times (scatter plot with marker="|" or "+")
5. Add vertical dashed lines at median survival times (if `--show-median`)
6. Add horizontal line at S(t)=0.5 (if `--show-median`)
7. Add at-risk table below x-axis (ax.text for each time point × group, showing n at risk)
8. Annotate p-value in corner (if `--show-pvalue`)
9. Customize labels, legend, colors, fonts, axis limits

### Step 5: Output generation
1. Save figure to PNG at specified DPI (default 300)
2. Save to SVG if requested
3. Write `survival_stats.tsv` with summary statistics per group
4. Print summary to console

## Key execution patterns

```bash
# Basic: Two-group KM curve (control vs. treatment)
python scripts/plot_survival.py \
  --input survival_data.csv \
  --time-col time_months \
  --event-col event \
  --group-col treatment \
  --output survival_plot.png \
  --title "Overall Survival by Treatment Group" \
  --show-pvalue \
  --show-at-risk

# Median stratification for gene expression (high vs. low)
python scripts/plot_survival.py \
  --input tcga_data.tsv \
  --time-col os_months \
  --event-col os_event \
  --split-by-median \
  --split-col gene_expression_level \
  --output gene_expr_survival.png \
  --title "Overall Survival by Gene Expression" \
  --show-median \
  --show-ci

# Multiple groups with custom colors and log-rank test
python scripts/plot_survival.py \
  --input clinical_trial.csv \
  --time-col time_days \
  --event-col event_status \
  --group-col dose_group \
  --group-order "placebo,low_dose,high_dose" \
  --group-colors "#e41a1c,#377eb8,#4daf4a" \
  --output dose_response_survival.png \
  --test logrank \
  --show-pvalue \
  --time-unit months \
  --fig-width 10 \
  --fig-height 7

# Clinical publication-quality figure
python scripts/plot_survival.py \
  --input cancer_cohort.tsv \
  --time-col efs_months \
  --event-col efs_event \
  --group-col biomarker_status \
  --output publication_figure.png \
  --output-svg \
  --title "Event-Free Survival by Biomarker Status" \
  --xlabel "Months from diagnosis" \
  --ylabel "Event-free survival probability" \
  --show-at-risk \
  --show-median \
  --show-ci \
  --show-pvalue \
  --test both \
  --dpi 300 \
  --fig-width 9 \
  --fig-height 6

# Compare log-rank and Wilcoxon tests
python scripts/plot_survival.py \
  --input early_stopping_trial.csv \
  --time-col weeks \
  --event-col event \
  --group-col arm \
  --output early_events.png \
  --test both \
  --show-pvalue \
  --pvalue-location top_left

# Stratified analysis (multiple strata)
# Create stratified dataset, then run separately:
python scripts/plot_survival.py \
  --input subset_stratum_a.csv \
  --time-col time \
  --event-col event \
  --group-col treatment \
  --output survival_stratum_a.png
```

## Parameter decision guide

| Parameter | Value | When to use | Rationale |
|-----------|-------|-------------|-----------|
| `--time-col` | e.g., "time_months" | Standard; always required | Different datasets use different column names |
| `--event-col` | e.g., "event_status" | Standard; always required | Must specify which column codes censoring |
| `--group-col` | e.g., "treatment" | Comparing 2+ groups | Leave empty for single-group KM |
| `--split-by-median` | True | High-dimensional data (gene expression, biomarkers) | Dichotomize continuous variable at median |
| `--split-col` | "expression_level" | Used with `--split-by-median` | Specify the continuous variable to split |
| `--show-ci` | True | Publication figures | Standard in clinical/biomedical publications |
| `--show-at-risk` | True | Clinical presentations | Shows how many subjects at risk at each time |
| `--show-median` | True | When median survival meaningful | Most cancer/event-free survival analyses |
| `--show-pvalue` | True | Comparative studies | Report p-value to show statistical significance |
| `--test` | "logrank" (default) | Most time-to-event data | Standard test for overall group differences |
| `--test` | "wilcoxon" | Early events more relevant | More weight on early time points |
| `--test` | "both" | Robust reporting | Show both tests to assess consistency |
| `--time-unit` | "months" | Human-readable labels | Adjust axis label (days, months, years) |
| `--xlim` | 24 | Focus on early follow-up | Truncate x-axis to relevant time window |
| `--dpi` | 300 (default) | Publication submission | High resolution for print/PDF |
| `--dpi` | 150 | Slide presentations | Adequate for screens, smaller file size |
| `--fig-width` | 8-10 | Standard presentation width | 8 for slides, 9-10 for papers |
| `--fig-height` | 6-7 | Standard presentation height | 6 for compact figure, 7 for readability |
| `--palette` | "Set1" (default) | Most applications | Good color contrast and colorblind-friendly |
| `--palette` | "tab10" | Many groups (5+) | More distinct colors than Set1 |
| `--palette` | "custom" | Custom brand colors | Provide via `--group-colors` |
| `--linewidth` | 2.0 (default) | Standard line width | Increase to 3+ for large poster figures |
| `--censoring-marks` | True (default) | Publication figures | Shows where censoring occurred |
| `--censoring-symbol` | "+" (default) | Standard in clinical papers | Can change to "o" or other marker |
| `--pvalue-location` | "top_right" (default) | Most figures | Avoids overlapping survival curves |

## Failure modes

| Failure | Cause | Solution |
|---------|-------|----------|
| **Error: No event column** | Column name mismatch | Check exact column name in CSV, use `--event-col exact_name` |
| **Error: Group has 0 events** | All subjects censored in a group | Script skips CI for that group; check data |
| **Empty plot (no curves)** | All times missing or negative | Validate time values are numeric and positive |
| **Identical curves for all groups** | No real group difference in data | Check grouping variable and data assignment |
| **p-value = NaN** | Insufficient events for valid test | May occur with <5 total events; report but interpret cautiously |
| **At-risk table overlaps plot** | Too many time points or groups | Reduce `--xlim` or number of groups |
| **Median survival = infinity** | More than 50% censoring in group | Script handles; report as "not reached" |
| **CI band too wide** | Small sample size or few events | Normal; indicates high uncertainty |
| **Censoring marks hide curves** | Too many censoring points at same time | Reduce marker size or use transparency |
| **Colors not appearing** | Custom hex codes invalid | Use format "#RRGGBB" with capital hex digits |
| **SVG file huge (>50 MB)** | Very large dataset with many points | Reduce DPI or simplify plot (fewer groups) |
| **Figure cuts off labels** | Axis labels too long | Shorten via `--title`, use `--fig-width 12` |

