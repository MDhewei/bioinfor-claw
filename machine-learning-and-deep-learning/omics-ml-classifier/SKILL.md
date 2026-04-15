---
name: omics-ml-classifier
description: Train and evaluate machine learning classifiers on omics data (gene expression, mutation scores, copy number, or any numeric feature matrix). Supports Random Forest, Logistic Regression, SVM, and XGBoost. Outputs cross-validation metrics, test set performance, feature importances, ROC curve, confusion matrix, and optional UMAP embedding. Works with any two-class or multi-class label.
---

# Omics ML Classifier

## Purpose

Train a **supervised machine learning classifier** on omics data and evaluate its performance.

Typical use cases:
- Classify samples by cancer subtype from gene expression
- Predict treatment response from mutation or expression profiles
- Distinguish tumor from normal from multi-omics features
- Identify the most predictive genes or biomarkers for a phenotype

Supported models:
- `random_forest` — Random Forest (default; handles high-dimensional data, provides feature importance)
- `logistic` — Logistic Regression with L2 regularisation
- `svm` — Support Vector Machine (RBF kernel)
- `xgboost` — XGBoost gradient boosting (requires `xgboost` package)

Default model: `random_forest`

## Reuse policy

This skill is designed for:
- Phenotype classification from any numeric omics matrix
- Biomarker discovery (via feature importance ranking)
- Model benchmarking across multiple classifier types
- Downstream validation of DE or pathway analysis results

This skill requires **user-provided local files** (feature matrix + label file).

## Inputs

### Required
- `--features` — feature matrix (TSV or CSV)
  - Rows = samples
  - Columns = features (gene names, protein IDs, etc.)
  - Values = numeric (expression, mutation score, etc.)

- `--labels` — label file (TSV or CSV)
  - Must have columns: `sample`, `label`
  - `sample` must match row names in `--features`
  - `label` is the class to predict (e.g. `responder`, `non-responder`)

- `--outdir` — output directory

### Optional

- `--model`
  Choices: `random_forest`, `logistic`, `svm`, `xgboost`
  Default: `random_forest`

- `--select-k`
  Default: None (use all features)
  Select top-K features by ANOVA F-score before training.
  Recommended for very high-dimensional data (> 5000 features).

- `--test-size`
  Default: `0.2`
  Fraction of samples held out for final test evaluation.

- `--cv-folds`
  Default: `5`
  Number of stratified cross-validation folds.

- `--random-state`
  Default: `42`
  Random seed for reproducibility.

- `--prefix`
  Default: `<model>_classifier`
  Output file prefix.

- `--top-n-features`
  Default: `30`
  Number of top features shown in the importance plot.

- `--umap`
  Flag. Generate UMAP embedding coloured by class label.
  Requires `umap-learn` package.

## Input file formats

### Feature matrix (--features)

```
sample      GENE1    GENE2    GENE3
sample_001  12.4     0.3      8.7
sample_002  3.1      11.2     2.0
sample_003  14.8     0.1      9.1
```

### Label file (--labels)

```
sample      label
sample_001  responder
sample_002  non-responder
sample_003  responder
```

## Outputs

### Always
- `<prefix>.cv_results.tsv` — per-fold and mean CV metrics (accuracy, ROC-AUC, F1)
- `<prefix>.test_metrics.tsv` — full classification report on held-out test set
- `<prefix>.feature_importance.tsv` — feature importance or coefficient magnitudes, ranked
- `ml_summary.tsv` — one-row summary of all key metrics
- `summary.txt` — human-readable summary

### Plots
- `<prefix>.roc_curve.png` / `.pdf` — ROC curve with AUC annotation
- `<prefix>.confusion_matrix.png` / `.pdf` — normalised confusion matrix
- `<prefix>.feature_importance.png` / `.pdf` — top feature importance bar chart
- `<prefix>.umap.png` / `.pdf` — UMAP embedding (only when `--umap` is set)

## Execution policy

Construct the command from required arguments and optional parameters.

### Command template

```
python scripts/omics_ml_classifier.py \
  --features <FEATURES_FILE> \
  --labels <LABELS_FILE> \
  [--model <MODEL>] \
  [--select-k <K>] \
  [--test-size <FRAC>] \
  [--cv-folds <N>] \
  [--random-state <SEED>] \
  [--prefix <PREFIX>] \
  [--top-n-features <N>] \
  [--umap] \
  --outdir <OUTDIR>
```

### Example commands

#### Default Random Forest
```
python scripts/omics_ml_classifier.py \
  --features expression_matrix.tsv \
  --labels sample_labels.tsv \
  --outdir results/
```

#### Logistic Regression with feature selection
```
python scripts/omics_ml_classifier.py \
  --features expr.tsv \
  --labels labels.tsv \
  --model logistic \
  --select-k 200 \
  --outdir results/
```

#### XGBoost with UMAP visualisation
```
python scripts/omics_ml_classifier.py \
  --features multi_omics.tsv \
  --labels subtype_labels.tsv \
  --model xgboost \
  --umap \
  --prefix xgb_subtype \
  --outdir results/
```

#### SVM with 10-fold CV
```
python scripts/omics_ml_classifier.py \
  --features expr.tsv \
  --labels labels.tsv \
  --model svm \
  --cv-folds 10 \
  --outdir results/
```

## Downstream skill chaining

| This skill's output                   | Next skill                   | How                                         |
|---------------------------------------|------------------------------|---------------------------------------------|
| `<prefix>.feature_importance.tsv`     | `go-analysis-for-gene-list`  | Pass top features as gene list              |
| `<prefix>.feature_importance.tsv`     | `annotation-for-gene-list`   | Annotate top predictive genes               |
| Top feature genes                     | `depmap-analysis-for-gene`   | Validate biomarker in DepMap                |
| Top feature genes                     | `tcga-survival-for-gene`     | Test survival impact of top biomarker       |

## Parameter decision guide

| Signal in user request | Parameter to set |
|---|---|
| High-dimensional omics (> 1000 genes) | `--model random_forest` (default; handles dimensionality well) |
| "interpretable / linear model" | `--model logistic` |
| Small dataset (< 100 samples) | `--model svm` |
| "gradient boosting / XGBoost" | `--model xgboost` |
| Very high dimensions (> 5000 features) | add `--select-k 500` to reduce features before training |
| "top 200 genes" (feature selection) | `--select-k 200` |
| "show UMAP embedding" | add `--umap` flag (requires `umap-learn`) |
| "10-fold cross-validation" | `--cv-folds 10` |
| "small dataset, leave-one-out-like" | `--cv-folds 10` with small test-size |
| "holdout test set = 30%" | `--test-size 0.3` |
| "show top 50 features" | `--top-n-features 50` |
| "reproducible results" | `--random-state 42` (default; set explicitly if comparing runs) |



| Scenario                                      | Recommended model   |
|-----------------------------------------------|---------------------|
| High-dimensional omics (> 1000 features)      | `random_forest`     |
| Need interpretable coefficients               | `logistic`          |
| Small dataset (< 100 samples)                 | `svm`               |
| Structured/tabular multi-omics                | `xgboost`           |
| First pass / exploratory                      | `random_forest`     |

## Failure conditions

Fail clearly if:
- `--features`, `--labels`, or `--outdir` are missing
- No common samples between feature matrix and label file
- Label file is missing `sample` or `label` columns
- Feature matrix has no numeric columns after NaN filtering
- `--model xgboost` is requested but `xgboost` is not installed
- Fewer samples than `--cv-folds` (reduce `--cv-folds`)

## Agent trigger examples

**Trigger this skill when the user asks:**
- "Build a classifier to predict cancer subtype from expression data"
- "Which genes best distinguish responders from non-responders?"
- "Train a Random Forest on my omics matrix"
- "Can you predict treatment response from this gene expression data?"
- "Run ML classification with cross-validation on my dataset"
- "What are the most important features for classifying tumor vs normal?"
- "Use XGBoost to classify my samples and show me a UMAP"

**Do NOT trigger this skill when the user asks:**
- "Cluster my samples" → use unsupervised clustering (not yet available)
- "Run differential expression" → use `rnaseq-differential-expression`
- "Train a deep learning model" → neural network skill (not yet available)
- "Predict drug sensitivity scores" → DepMap analysis, use `depmap-analysis-for-gene`

## Notes

- All features are standardised (zero mean, unit variance) before training.
- Missing values (NaN) in the feature matrix are imputed with column medians.
- `random_forest` and `xgboost` provide native feature importances (Gini / gain).
  `logistic` and `svm` use absolute coefficient values as proxy importances.
- ROC-AUC is reported as weighted OvR (one-vs-rest) for multi-class problems.
- For imbalanced classes, all models use `class_weight="balanced"` (except XGBoost).
- `--select-k` applies ANOVA F-score feature selection on the **training fold only**
  to avoid data leakage.
- `--umap` requires `umap-learn` (`pip install umap-learn`); skipped gracefully if absent.
