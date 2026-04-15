---
name: cox-survival-analysis
description: Run Cox proportional hazards regression on clinical and molecular data. Fits univariate and multivariate Cox models, computes hazard ratios with confidence intervals, tests the proportional hazards assumption, and generates forest plots and Schoenfeld residual plots.
---

## Purpose

Perform survival analysis using Cox proportional hazards regression on clinical and molecular datasets. This skill automates univariate and multivariate Cox model fitting, hazard ratio estimation, proportional hazards assumption testing, and publication-ready visualization. Pure NumPy/SciPy implementation enables reproducible, transparent survival analyses without external C dependencies.

## Use When

- **Clinical outcome prediction**: Associating molecular/clinical variables with time-to-event (survival, recurrence, disease progression)
- **Univariate screening**: Identifying significant predictors before multivariate modeling
- **Multivariable risk modeling**: Building comprehensive Cox models accounting for multiple covariates
- **Proportional hazards validation**: Testing model assumptions with Schoenfeld residual plots
- **Publication-ready results**: Generating forest plots and summary tables for manuscripts
- **Dataset size**: 50–10,000 samples with 5–500 covariates

## Do NOT Use When

- **Proportional hazards severely violated**: Use accelerated failure time (AFT) models instead
- **Time-varying covariates required**: This skill assumes time-constant effects; use stratified or extended Cox models
- **Small samples (n < 20)**: High uncertainty in Cox estimates; consider Bayesian approaches
- **Competing risks**: Use Fine-Gray competing risks model instead
- **High-dimensional data (p >> n)**: Use penalized Cox (ridge/lasso) or external tools

## Expected Inputs

- **Clinical/molecular data** (required): TSV file with samples as rows and features as columns
  - First column: sample IDs
  - Numeric data for continuous covariates
  - String/categorical data for factors (auto-detected and one-hot encoded)
  - Missing values: specify as NA or empty cells (will be dropped with warning)

- **Time column** (required): Column name with time-to-event (days, years, etc.)
  - Must be numeric and ≥ 0
  - Interpretation: days/months/years from baseline to event or censoring

- **Event column** (required): Binary indicator of event occurrence
  - Values: 0 (censored) or 1 (event observed)
  - Other values: will be coerced to 0/1 (0 = no event, >0 = event)

- **Covariates** (required): Comma-separated list of column names to include in model
  - Can be continuous (numeric) or categorical (string/factor)

## Expected Outputs

All files written to `--outdir/`:

- **univariate_results.tsv**: Univariate Cox model results
  - Columns: covariate, HR, CI_lower, CI_upper, pvalue, significant

- **multivariate_results.tsv**: Multivariate Cox model results
  - Columns: covariate, HR, CI_lower, CI_upper, pvalue, significant

- **ph_test.tsv**: Proportional hazards assumption test results
  - Columns: covariate, correlation, pvalue, ph_violated (TRUE if p < 0.05)

- **risk_scores.tsv**: Per-sample Cox linear predictor and risk group assignment
  - Columns: sample_id, linear_predictor, risk_group (low/high)

- **forest_plot.png**: Horizontal forest plot showing HRs, 95% CIs, p-values
  - Left: univariate estimates; Right: multivariate estimates

- **schoenfeld_residuals.png**: Scatter plots of Schoenfeld residuals vs log(time) per covariate
  - Shows violations of proportional hazards assumption

- **survival_curves.png**: Kaplan-Meier curves for low/high predicted risk groups from Cox model

## Procedure

### Step 1: Data Loading and Preparation
1. Load TSV data
2. Identify numeric vs. categorical columns
3. One-hot encode categorical variables (drop reference level to avoid multicollinearity)
4. Standardize continuous covariates (z-score normalization)
5. Perform complete-case analysis: drop samples with missing values in time, event, or covariates

### Step 2: Event/Time Validation
- Verify time ≥ 0 for all samples
- Coerce event to binary: event > 0 → 1, else 0
- Report number of events and censored samples
- Sort data by event time (ascending)

### Step 3: Cox Partial Likelihood Optimization
For each Cox model (univariate or multivariate):

1. **Partial likelihood function**:
   ```
   L(β) = ∏_i [ exp(X_i·β) / Σ_{j ∈ R(t_i)} exp(X_j·β) ]^δ_i

   where R(t_i) is the risk set (all samples still at risk at time t_i)
   and δ_i is the event indicator.
   ```

2. **Log partial likelihood** and gradient (score function):
   ```
   ℓ(β) = Σ_i δ_i [ X_i·β - log(Σ_{j ∈ R(t_i)} exp(X_j·β)) ]

   ∂ℓ/∂β = Σ_i δ_i [ X_i - Σ_{j ∈ R(t_i)} w_j·X_j ]
   ```

3. **Hessian** (observed information matrix):
   ```
   ∂²ℓ/∂β² (numerically via finite differences)
   ```

4. **Newton-Raphson iteration**:
   ```
   β_{t+1} = β_t - H^{-1} · ∇ℓ(β_t)
   ```
   - Max 20 iterations
   - Convergence: ||∇ℓ|| < 1e-6 or relative change < 1e-6

### Step 4: Inference on Coefficients
- Covariance matrix: Var(β̂) = -H^{-1}
- Standard errors: SE = sqrt(diag(Var(β̂)))
- Wald test: z = β/SE, p-value = 2·(1 - Φ(|z|)) [two-tailed normal test]
- Hazard ratios: HR = exp(β)
- 95% CIs: exp(β ± 1.96·SE)

### Step 5: Univariate Analysis
Run separate Cox model for each covariate independently:
- Fit single-covariate Cox model
- Report HR, CI, p-value
- Use to screen for significant predictors

### Step 6: Multivariate Analysis
Fit full Cox model with all specified covariates simultaneously:
- Includes all covariates from --covariates list
- Reports HR, CI, p-value adjusted for other covariates
- Assess convergence and data separability

### Step 7: Proportional Hazards Assumption Test
For significant covariates (p < 0.05 in multivariate model):

1. **Schoenfeld residuals**:
   ```
   r_ij = X_ij - E[X_j | t_i]

   E[X_j | t_i] = Σ_{k ∈ R(t_i)} w_k·X_kj
   ```

2. **Test correlation** between residuals and log(time):
   - Spearman rank correlation
   - p-value from Spearman test
   - If p < 0.05: PH assumption violated for that covariate

3. **Visual inspection**: Scatter plots with smoothed LOWESS curves

### Step 8: Risk Stratification
1. Compute linear predictor for each sample: LP = X·β̂ (from multivariate model)
2. Stratify by median LP: risk_group = "low" if LP < median, else "high"
3. Compute Kaplan-Meier curves per risk group
4. Plot with log-rank test p-value

### Step 9: Visualization and Output
- **Forest plot**: Horizontal bars showing HR and 95% CI per covariate, both univariate (left) and multivariate (right)
  - Vertical line at HR = 1 (no effect)
  - Color code: significant (p < 0.05) vs. not significant
  - Include p-values as text

- **Schoenfeld residual plots**: Scatter with LOWESS smoothing per covariate
  - If smooth line deviates from 0, PH assumption violated
  - Annotate with Spearman p-value

- **Survival curves**: Kaplan-Meier curves for low/high risk groups
  - Include number at risk table
  - Log-rank test p-value in title

## Key Execution Patterns

### Basic univariate screening:
```bash
python scripts/cox_survival_analysis.py \
  --input clinical_data.tsv \
  --time-col overall_survival_days \
  --event-col vital_status \
  --covariates age,gene_expr_score,stage \
  --mode univariate \
  --outdir ./results
```

### Complete analysis with categorical variables:
```bash
python scripts/cox_survival_analysis.py \
  --input clinical_data.tsv \
  --time-col time_to_event \
  --event-col event \
  --covariates age,tumor_stage,mutation_status,gene_score \
  --categorical tumor_stage,mutation_status \
  --mode both \
  --outdir ./results
```

### Multivariate model with continuous covariates only:
```bash
python scripts/cox_survival_analysis.py \
  --input clinical_data.tsv \
  --time-col recurrence_free_survival \
  --event-col recurrence \
  --covariates age,psa,gleason_score \
  --continuous age,psa,gleason_score \
  --mode multivariate \
  --alpha 0.01 \
  --outdir ./results
```

### With interaction terms:
```bash
python scripts/cox_survival_analysis.py \
  --input data.tsv \
  --time-col time \
  --event-col event \
  --covariates age,treatment,marker \
  --interaction-terms "treatment:marker" \
  --outdir ./results
```

### Stratified analysis (e.g., by center):
```bash
python scripts/cox_survival_analysis.py \
  --input multi_center_data.tsv \
  --time-col time \
  --event-col event \
  --covariates age,score \
  --strata center \
  --outdir ./results
```

## Parameter Decision Guide

| Parameter | Typical Values | When to Adjust | Example Scenario |
|-----------|----------------|----------------|------------------|
| `--mode` | both, univariate, multivariate | Start with univariate, then multivariate | Screening: univariate; publication: both |
| `--categorical` | comma-separated column names | Auto-detected as categorical if non-numeric | Tumor stage (1,2,3,4) → categorical; age → continuous |
| `--continuous` | all numeric by default | Specify if want z-score normalization | Gene expression continuous; explicit scaling |
| `--alpha` | 0.05, 0.01 | More stringent for large-scale analyses | Candidate screening: 0.05; genome-wide: 0.01 |
| `--reference-categories` | "col:ref_val" | Manual ref level selection | stage:1, treatment:placebo |
| `--interaction-terms` | "col1:col2" comma-separated | Add if expecting synergy | treatment:biomarker interactions |
| `--strata` | column name | Use if have natural groupings | Multi-center, multiple batches |

## Failure Modes & Troubleshooting

| Problem | Likely Cause | Solution |
|---------|--------------|----------|
| "Newton-Raphson did not converge" | Quasi-complete/complete separation in data | Check for 0 variance covariates; remove perfectly predictive variables |
| "Singular Hessian matrix" | Perfect multicollinearity among covariates | Check correlation between variables; remove highly correlated predictors |
| "No events (all censored)" | Event column has no 1s | Verify event column coding; check data quality |
| "Negative hazard ratios" | Interpretation error | HR is always > 0; if exp(β) < 0, check coefficient output |
| "Very large/small p-values" | Extreme effect sizes or overfitting | Check for separation; reduce number of covariates if n is small |
| "PH assumption violated for all covariates" | Model misspecification or time-varying effects | Consider stratification by covariate; use extended Cox |
| "Risk curves not separated" | Weak predictive signal | Investigate covariate quality; consider feature engineering |

## Advanced Options

- **Stratified Cox**: Use `--strata` to fit separate baseline hazards per stratum
- **Reference categories**: Manually set reference levels via `--reference-categories col:value`
- **Interaction terms**: Include product terms via `--interaction-terms col1:col2`
- **Feature engineering**: Pre-compute derived variables in input file (e.g., log(gene_expr))

## Statistical Notes

- **Hazard ratio interpretation**: HR = 1.5 means 50% increase in instantaneous risk per unit increase in covariate
- **Log scale**: Forest plots show HR on log scale for symmetry
- **P-values**: Two-tailed Wald tests; 95% CIs always 1.96·SE for normal approximation
- **Proportional hazards**: Can be relaxed by stratification or extended Cox; see Therneau & Grambsch (2000)

## References

- Therneau, T. M., & Grambsch, P. M. (2000). *Modeling Survival Data: Extending the Cox Model*. Springer.
- Cox, D. R. (1972). "Regression Models and Life-Tables". *Journal of the Royal Statistical Society*, 34(2), 187-220.
- Kalbfleisch, J. D., & Prentice, R. L. (2002). *The Statistical Analysis of Failure Time Data* (2nd ed.). Wiley.
