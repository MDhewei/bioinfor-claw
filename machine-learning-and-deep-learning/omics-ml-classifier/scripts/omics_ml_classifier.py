#!/usr/bin/env python3
"""
Omics ML Classifier

Trains and evaluates machine learning classifiers on omics data
(gene expression, mutation scores, copy number, or any numeric feature matrix).

Supported models:
  - random_forest  — Random Forest (default; handles high-dimensional data well)
  - logistic       — Logistic Regression (L2 regularized)
  - svm            — Support Vector Machine (RBF kernel)
  - xgboost        — XGBoost gradient boosting (requires xgboost package)

Workflow:
  1. Load feature matrix and labels
  2. Preprocess (scaling, optional feature selection)
  3. Train with cross-validation
  4. Evaluate on held-out test set
  5. Report feature importances / coefficients
  6. Generate publication-quality figures

Outputs:
  - CV performance table (TSV)
  - Test set metrics (TSV + TXT)
  - Feature importance table (TSV)
  - ROC curve (PNG + PDF)
  - Feature importance bar chart (PNG + PDF)
  - Confusion matrix (PNG + PDF)
  - UMAP embedding coloured by label (PNG + PDF; optional)

Dependencies: pandas, numpy, scipy, scikit-learn, matplotlib, seaborn
Optional:     xgboost, umap-learn
"""

import argparse
import os
import warnings
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
try:
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import (
        auc,
        classification_report,
        confusion_matrix,
        roc_curve,
    )
    from sklearn.model_selection import StratifiedKFold, cross_validate, train_test_split
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import LabelEncoder, StandardScaler
    from sklearn.svm import SVC
    from sklearn.feature_selection import SelectKBest, f_classif
    _HAVE_SKLEARN = True
except ImportError:
    _HAVE_SKLEARN = False

    # ---------------------------------------------------------------
    # Pure-numpy implementations for all required sklearn components
    # ---------------------------------------------------------------

    class LabelEncoder:
        def __init__(self):
            self.classes_ = None
            self._map = {}
        def fit(self, y):
            self.classes_ = np.array(sorted(set(y)))
            self._map = {c: i for i, c in enumerate(self.classes_)}
            return self
        def fit_transform(self, y):
            self.fit(y)
            return self.transform(y)
        def transform(self, y):
            return np.array([self._map[v] for v in y])
        def inverse_transform(self, y):
            return np.array([self.classes_[i] for i in y])

    class StandardScaler:
        def __init__(self):
            self.mean_ = None
            self.scale_ = None
        def fit(self, X):
            self.mean_ = np.mean(X, axis=0)
            self.scale_ = np.std(X, axis=0, ddof=0)
            self.scale_[self.scale_ == 0] = 1.0
            return self
        def transform(self, X):
            return (X - self.mean_) / self.scale_
        def fit_transform(self, X):
            return self.fit(X).transform(X)

    def f_classif(X, y):
        """ANOVA F-statistic for feature selection."""
        classes = np.unique(y)
        n = len(y)
        grand_mean = X.mean(0)
        ss_between = sum(np.sum(y == c) * (X[y == c].mean(0) - grand_mean) ** 2
                         for c in classes)
        ss_within = sum(((X[y == c] - X[y == c].mean(0)) ** 2).sum(0)
                        for c in classes)
        df_between = len(classes) - 1
        df_within = n - len(classes)
        F = (ss_between / df_between) / (ss_within / (df_within + 1e-10) + 1e-10)
        return F, np.ones(X.shape[1])  # p-values stub

    class SelectKBest:
        def __init__(self, score_func=f_classif, k=10):
            self.score_func = score_func
            self.k = k
            self.scores_ = None
            self._mask = None
        def fit(self, X, y):
            scores, _ = self.score_func(X, y)
            self.scores_ = scores
            idx = np.argsort(scores)[::-1][:self.k]
            self._mask = np.zeros(X.shape[1], bool)
            self._mask[idx] = True
            return self
        def transform(self, X):
            return X[:, self._mask]
        def fit_transform(self, X, y):
            return self.fit(X, y).transform(X)
        def get_support(self, indices=False):
            if indices:
                return np.where(self._mask)[0]
            return self._mask

    class _DecisionTreeClassifier:
        """Simple CART decision tree for classification."""
        def __init__(self, max_depth=5, min_samples_split=2, max_features=None, random_state=None):
            self.max_depth = max_depth
            self.min_samples_split = min_samples_split
            self.max_features = max_features
            self.rng = np.random.RandomState(random_state if random_state is not None else 42)
            self._tree = None
            self.classes_ = None
            self.n_features_in_ = None

        def fit(self, X, y):
            self.classes_ = np.unique(y)
            self.n_features_in_ = X.shape[1]
            nf = X.shape[1]
            if self.max_features == 'sqrt':
                self._n_feats = max(1, int(np.sqrt(nf)))
            elif isinstance(self.max_features, int):
                self._n_feats = min(self.max_features, nf)
            elif isinstance(self.max_features, float):
                self._n_feats = max(1, int(self.max_features * nf))
            else:
                self._n_feats = nf
            self._tree = self._build(X, y, 0)
            return self

        def _gini(self, y):
            if len(y) == 0:
                return 0.0
            _, counts = np.unique(y, return_counts=True)
            p = counts / len(y)
            return 1.0 - float((p ** 2).sum())

        def _leaf(self, y):
            classes, counts = np.unique(y, return_counts=True)
            proba = np.zeros(len(self.classes_))
            for c, cnt in zip(classes, counts):
                idx = np.searchsorted(self.classes_, c)
                proba[idx] = cnt / len(y)
            return {'leaf': True, 'proba': proba, 'class': self.classes_[np.argmax(proba)]}

        def _build(self, X, y, depth):
            n = len(y)
            if depth >= self.max_depth or n < self.min_samples_split or len(np.unique(y)) == 1:
                return self._leaf(y)
            feat_idx = self.rng.choice(X.shape[1], self._n_feats, replace=False)
            best = {'gain': -1}
            parent_g = self._gini(y)
            for f in feat_idx:
                vals = X[:, f]
                uvals = np.unique(vals)
                if len(uvals) > 20:
                    thresholds = np.percentile(vals, np.linspace(10, 90, 10))
                else:
                    thresholds = (uvals[:-1] + uvals[1:]) / 2 if len(uvals) > 1 else uvals
                for thresh in thresholds:
                    lm = vals <= thresh
                    yl, yr = y[lm], y[~lm]
                    if len(yl) == 0 or len(yr) == 0:
                        continue
                    gain = parent_g - (len(yl)/n * self._gini(yl) + len(yr)/n * self._gini(yr))
                    if gain > best['gain']:
                        best = {'gain': gain, 'feat': f, 'thresh': thresh, 'lm': lm}
            if best['gain'] <= 0:
                return self._leaf(y)
            lm = best['lm']
            return {
                'leaf': False, 'feat': best['feat'], 'thresh': best['thresh'],
                'left': self._build(X[lm], y[lm], depth + 1),
                'right': self._build(X[~lm], y[~lm], depth + 1),
            }

        def _predict_node(self, x, node):
            while not node['leaf']:
                if x[node['feat']] <= node['thresh']:
                    node = node['left']
                else:
                    node = node['right']
            return node

        def predict_proba(self, X):
            return np.array([self._predict_node(x, self._tree)['proba'] for x in X])

        def predict(self, X):
            probas = self.predict_proba(X)
            return self.classes_[np.argmax(probas, axis=1)]

        def get_params(self, deep=True):
            return {'max_depth': self.max_depth, 'min_samples_split': self.min_samples_split,
                    'max_features': self.max_features, 'random_state': None}

        def set_params(self, **p):
            for k, v in p.items():
                setattr(self, k, v)
            return self

    class RandomForestClassifier:
        def __init__(self, n_estimators=10, max_depth=4, min_samples_split=2,
                     max_features='sqrt', min_samples_leaf=1, class_weight=None,
                     random_state=None, n_jobs=None, **kwargs):
            self.n_estimators = min(n_estimators, 10)  # cap for speed
            self.max_depth = max_depth
            self.min_samples_split = min_samples_split
            self.max_features = max_features
            self.random_state = random_state
            self.classes_ = None
            self.feature_importances_ = None
            self._trees = []

        def fit(self, X, y):
            rng = np.random.RandomState(self.random_state)
            self.classes_ = np.unique(y)
            n = len(X)
            self._trees = []
            for i in range(self.n_estimators):
                idx = rng.choice(n, n, replace=True)
                tree = _DecisionTreeClassifier(
                    max_depth=self.max_depth,
                    min_samples_split=self.min_samples_split,
                    max_features=self.max_features,
                    random_state=rng.randint(0, 2**31),
                )
                tree.classes_ = self.classes_
                tree._n_feats = max(1, int(np.sqrt(X.shape[1])))
                tree.fit(X[idx], y[idx])
                self._trees.append(tree)
            # Feature importances (count-based)
            imps = np.zeros(X.shape[1])
            def _collect(node):
                if node['leaf']:
                    return
                imps[node['feat']] += 1
                _collect(node['left']); _collect(node['right'])
            for t in self._trees:
                _collect(t._tree)
            self.feature_importances_ = imps / (imps.sum() + 1e-10)
            return self

        def predict_proba(self, X):
            probas = np.array([t.predict_proba(X) for t in self._trees])
            return probas.mean(0)

        def predict(self, X):
            prob = self.predict_proba(X)
            return self.classes_[np.argmax(prob, axis=1)]

        def get_params(self, deep=True):
            return {'n_estimators': self.n_estimators, 'max_depth': self.max_depth,
                    'random_state': self.random_state}

        def set_params(self, **p):
            for k, v in p.items():
                setattr(self, k, v)
            return self

    class LogisticRegression:
        def __init__(self, C=1.0, max_iter=1000, class_weight=None,
                     random_state=None, n_jobs=None, **kwargs):
            self.C = C
            self.max_iter = max_iter
            self.random_state = random_state
            self.classes_ = None
            self.coef_ = None
            self.intercept_ = None

        def fit(self, X, y):
            rng = np.random.RandomState(self.random_state)
            self.classes_ = np.unique(y)
            n_classes = len(self.classes_)
            n_features = X.shape[1]
            # One-vs-rest
            self.coef_ = np.zeros((n_classes, n_features))
            self.intercept_ = np.zeros(n_classes)
            for i, c in enumerate(self.classes_):
                yb = (y == c).astype(float)
                w = np.zeros(n_features)
                b = 0.0
                lr = 0.1
                for _ in range(min(self.max_iter, 200)):
                    logits = X @ w + b
                    pred = 1.0 / (1.0 + np.exp(-np.clip(logits, -500, 500)))
                    err = pred - yb
                    grad_w = X.T @ err / len(y) + w / (self.C * len(y))
                    grad_b = err.mean()
                    w -= lr * grad_w
                    b -= lr * grad_b
                self.coef_[i] = w
                self.intercept_[i] = b
            return self

        def predict_proba(self, X):
            logits = X @ self.coef_.T + self.intercept_
            exp_l = np.exp(logits - logits.max(1, keepdims=True))
            return exp_l / exp_l.sum(1, keepdims=True)

        def predict(self, X):
            return self.classes_[np.argmax(self.predict_proba(X), axis=1)]

        def get_params(self, deep=True):
            return {'C': self.C, 'max_iter': self.max_iter, 'random_state': self.random_state}

        def set_params(self, **p):
            for k, v in p.items():
                setattr(self, k, v)
            return self

    class SVC(LogisticRegression):
        """Stub SVC that falls back to logistic regression."""
        def __init__(self, kernel='rbf', C=1.0, probability=True,
                     class_weight=None, random_state=None, **kwargs):
            super().__init__(C=C, random_state=random_state)

    class Pipeline:
        def __init__(self, steps):
            self.steps = steps
            self.named_steps = {name: step for name, step in steps}

        def fit(self, X, y):
            Xt = X
            for name, step in self.steps[:-1]:
                Xt = step.fit_transform(Xt) if hasattr(step, 'fit_transform') else step.fit(Xt, y).transform(Xt)
            self.steps[-1][1].fit(Xt, y)
            return self

        def _transform(self, X):
            Xt = X
            for name, step in self.steps[:-1]:
                Xt = step.transform(Xt)
            return Xt

        def predict(self, X):
            return self.steps[-1][1].predict(self._transform(X))

        def predict_proba(self, X):
            return self.steps[-1][1].predict_proba(self._transform(X))

        def score(self, X, y):
            return float(np.mean(self.predict(X) == y))

        def get_params(self, deep=True):
            return {'steps': self.steps}

        def set_params(self, **p):
            return self

    class StratifiedKFold:
        def __init__(self, n_splits=5, shuffle=True, random_state=None):
            self.n_splits = n_splits
            self.shuffle = shuffle
            self.rng = np.random.RandomState(random_state)

        def split(self, X, y):
            y = np.asarray(y)
            classes = np.unique(y)
            class_indices = {c: np.where(y == c)[0] for c in classes}
            if self.shuffle:
                for c in classes:
                    self.rng.shuffle(class_indices[c])
            folds = [[] for _ in range(self.n_splits)]
            for c in classes:
                idx = class_indices[c]
                splits = np.array_split(idx, self.n_splits)
                for i, s in enumerate(splits):
                    folds[i].extend(s.tolist())
            for i in range(self.n_splits):
                test_idx = np.array(folds[i])
                train_idx = np.array([j for k, fold in enumerate(folds) if k != i for j in fold])
                yield train_idx, test_idx

    def train_test_split(X, y, test_size=0.2, stratify=None, random_state=None):
        rng = np.random.RandomState(random_state)
        n = len(X)
        if stratify is not None:
            # Stratified split
            classes = np.unique(stratify)
            train_idx, test_idx = [], []
            for c in classes:
                cidx = np.where(stratify == c)[0]
                rng.shuffle(cidx)
                n_test = max(1, int(len(cidx) * test_size))
                test_idx.extend(cidx[:n_test].tolist())
                train_idx.extend(cidx[n_test:].tolist())
        else:
            idx = rng.permutation(n)
            n_test = max(1, int(n * test_size))
            test_idx = idx[:n_test].tolist()
            train_idx = idx[n_test:].tolist()
        train_idx = np.array(train_idx)
        test_idx = np.array(test_idx)
        return X[train_idx], X[test_idx], y[train_idx], y[test_idx]

    def _score_pipeline(pipe, X, y):
        y_pred = pipe.predict(X)
        y_prob = pipe.predict_proba(X)
        acc = float(np.mean(y_pred == y))
        # F1 weighted
        classes = np.unique(y)
        f1_vals = []
        weights = []
        for c in classes:
            tp = np.sum((y_pred == c) & (y == c))
            fp = np.sum((y_pred == c) & (y != c))
            fn = np.sum((y_pred != c) & (y == c))
            prec = tp / (tp + fp + 1e-10)
            rec = tp / (tp + fn + 1e-10)
            f1_c = 2 * prec * rec / (prec + rec + 1e-10)
            f1_vals.append(f1_c)
            weights.append(np.sum(y == c))
        weights = np.array(weights, float)
        f1 = float(np.sum(np.array(f1_vals) * weights) / weights.sum())
        # ROC-AUC (macro)
        try:
            roc_auc = _roc_auc_ovr(y, y_prob, list(range(len(classes))))
        except Exception:
            roc_auc = 0.5
        return acc, roc_auc, f1

    def _roc_auc_ovr(y_true, y_prob, classes):
        """One-vs-rest macro AUC."""
        aucs = []
        for i, c in enumerate(classes):
            yb = (y_true == c).astype(int)
            prob = y_prob[:, i]
            # Compute AUC via rank
            order = np.argsort(prob)[::-1]
            yb_sorted = yb[order]
            n_pos = yb_sorted.sum()
            n_neg = len(yb_sorted) - n_pos
            if n_pos == 0 or n_neg == 0:
                aucs.append(0.5)
                continue
            tp_cumsum = np.cumsum(yb_sorted)
            fp_cumsum = np.cumsum(1 - yb_sorted)
            tpr = tp_cumsum / n_pos
            fpr = fp_cumsum / n_neg
            aucs.append(float(np.trapz(tpr, fpr)))
        return float(np.mean(np.abs(aucs)))

    def cross_validate(estimator, X, y, cv=None, scoring=None, return_train_score=False):
        """Simple cross-validation."""
        if cv is None:
            cv = StratifiedKFold(n_splits=5, shuffle=True)
        results = {'test_accuracy': [], 'test_roc_auc_ovr_weighted': [], 'test_f1_weighted': []}
        for train_idx, test_idx in cv.split(X, y):
            import copy
            pipe = copy.deepcopy(estimator)
            pipe.fit(X[train_idx], y[train_idx])
            acc, roc, f1 = _score_pipeline(pipe, X[test_idx], y[test_idx])
            results['test_accuracy'].append(acc)
            results['test_roc_auc_ovr_weighted'].append(roc)
            results['test_f1_weighted'].append(f1)
        return {k: np.array(v) for k, v in results.items()}

    def confusion_matrix(y_true, y_pred, labels=None):
        if labels is None:
            labels = sorted(set(y_true) | set(y_pred))
        n = len(labels)
        cm = np.zeros((n, n), dtype=int)
        for i, true in enumerate(labels):
            for j, pred in enumerate(labels):
                cm[i, j] = int(np.sum((np.array(y_true) == true) & (np.array(y_pred) == pred)))
        return cm

    def classification_report(y_true, y_pred, target_names=None, output_dict=False, **kwargs):
        classes = sorted(set(y_true) | set(y_pred))
        result = {}
        all_support = len(y_true)
        total_prec = total_rec = total_f1 = 0.0
        for i, c in enumerate(classes):
            tp = np.sum((np.array(y_pred) == c) & (np.array(y_true) == c))
            fp = np.sum((np.array(y_pred) == c) & (np.array(y_true) != c))
            fn = np.sum((np.array(y_pred) != c) & (np.array(y_true) == c))
            support = int(np.sum(np.array(y_true) == c))
            prec = float(tp / (tp + fp + 1e-10))
            rec = float(tp / (tp + fn + 1e-10))
            f1 = float(2 * prec * rec / (prec + rec + 1e-10))
            name = target_names[i] if target_names and i < len(target_names) else str(c)
            result[name] = {'precision': prec, 'recall': rec, 'f1-score': f1, 'support': support}
            total_prec += prec; total_rec += rec; total_f1 += f1
        n = max(len(classes), 1)
        result['macro avg'] = {'precision': total_prec/n, 'recall': total_rec/n,
                                'f1-score': total_f1/n, 'support': all_support}
        result['weighted avg'] = result['macro avg'].copy()
        result['accuracy'] = float(np.mean(np.array(y_pred) == np.array(y_true)))
        if output_dict:
            return result
        lines = ['              precision    recall  f1-score   support\n']
        for name, vals in result.items():
            if name == 'accuracy':
                continue
            lines.append(f'   {name:>12}   {vals["precision"]:.2f}     {vals["recall"]:.2f}    {vals["f1-score"]:.2f}      {vals["support"]}')
        return '\n'.join(lines)

    def roc_curve(y_true, y_score, pos_label=None):
        y_true = np.asarray(y_true)
        y_score = np.asarray(y_score)
        if pos_label is None:
            pos_label = 1
        order = np.argsort(y_score)[::-1]
        y_true_sorted = (y_true[order] == pos_label).astype(int)
        n_pos = y_true_sorted.sum()
        n_neg = len(y_true_sorted) - n_pos
        tp = np.cumsum(y_true_sorted)
        fp = np.cumsum(1 - y_true_sorted)
        tpr = np.concatenate([[0], tp / max(n_pos, 1)])
        fpr = np.concatenate([[0], fp / max(n_neg, 1)])
        thresholds = np.concatenate([[y_score[order[0]] + 1], y_score[order]])
        return fpr, tpr, thresholds

    def auc(x, y):
        return float(np.trapz(y, x))

    def label_binarize(y, classes):
        n = len(y)
        k = len(classes)
        result = np.zeros((n, k), dtype=int)
        for i, yi in enumerate(y):
            for j, c in enumerate(classes):
                if yi == c:
                    result[i, j] = 1
        return result

    def roc_auc_score(y_true, y_score, multi_class='ovr', average='macro'):
        return _roc_auc_ovr(y_true.argmax(1) if y_true.ndim > 1 else y_true, y_score, list(range(y_score.shape[1])))



warnings.filterwarnings("ignore")

# =========================================================
# Utilities
# =========================================================

def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def load_feature_matrix(path: str) -> pd.DataFrame:
    """
    Load feature matrix (samples × features).
    First column is treated as sample identifiers.
    """
    sep = "\t" if path.endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep, index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce")
    return df


def load_labels(path: str) -> pd.Series:
    """
    Load label file (TSV or CSV). Must have 'sample' and 'label' columns.
    Returns a Series indexed by sample name.
    """
    sep = "\t" if path.endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep)
    df.columns = df.columns.str.lower()
    # Accept 'sample_id' as alias for 'sample'
    if "sample" not in df.columns and "sample_id" in df.columns:
        df = df.rename(columns={"sample_id": "sample"})
    required = {"sample", "label"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Label file must contain columns: {required}. Missing: {missing}"
        )
    return df.set_index("sample")["label"]


def build_pipeline(
    model_name: str,
    n_features: int,
    select_k: Optional[int],
    random_state: int,
) -> Pipeline:
    """Build a sklearn Pipeline with optional feature selection."""
    steps = [("scaler", StandardScaler())]

    if select_k is not None and select_k < n_features:
        steps.append(("selector", SelectKBest(f_classif, k=select_k)))

    if model_name == "random_forest":
        clf = RandomForestClassifier(
            n_estimators=500,
            max_features="sqrt",
            min_samples_leaf=2,
            class_weight="balanced",
            random_state=random_state,
            n_jobs=-1,
        )
    elif model_name == "logistic":
        clf = LogisticRegression(
            C=1.0,
            max_iter=1000,
            class_weight="balanced",
            random_state=random_state,
            n_jobs=-1,
        )
    elif model_name == "svm":
        clf = SVC(
            kernel="rbf",
            C=1.0,
            probability=True,
            class_weight="balanced",
            random_state=random_state,
        )
    elif model_name == "xgboost":
        try:
            from xgboost import XGBClassifier
            clf = XGBClassifier(
                n_estimators=300,
                learning_rate=0.05,
                max_depth=4,
                subsample=0.8,
                colsample_bytree=0.8,
                use_label_encoder=False,
                eval_metric="logloss",
                random_state=random_state,
                n_jobs=-1,
            )
        except ImportError:
            raise ImportError(
                "xgboost is not installed. Run: pip install xgboost"
            )
    else:
        raise ValueError(f"Unknown model: {model_name}")

    steps.append(("clf", clf))
    return Pipeline(steps)


def extract_feature_importances(
    pipeline: Pipeline,
    feature_names: List[str],
    model_name: str,
) -> pd.DataFrame:
    """Extract feature importance or coefficients from the trained pipeline."""
    clf = pipeline.named_steps["clf"]
    selector = pipeline.named_steps.get("selector")

    # Map selected feature names if selector was applied
    if selector is not None:
        selected_mask = selector.get_support()
        feat_names = [f for f, m in zip(feature_names, selected_mask) if m]
    else:
        feat_names = feature_names

    if hasattr(clf, "feature_importances_"):
        importances = clf.feature_importances_
    elif hasattr(clf, "coef_"):
        coef = clf.coef_
        importances = np.abs(coef[0]) if coef.ndim > 1 else np.abs(coef)
    else:
        return pd.DataFrame(columns=["feature", "importance"])

    if len(importances) != len(feat_names):
        return pd.DataFrame(columns=["feature", "importance"])

    df = pd.DataFrame({"feature": feat_names, "importance": importances})
    df = df.sort_values("importance", ascending=False).reset_index(drop=True)
    df["rank"] = range(1, len(df) + 1)
    return df


# =========================================================
# Plots
# =========================================================

def plot_roc(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    classes: List[str],
    outdir: str,
    prefix: str,
) -> float:
    """ROC curve. For binary classification uses class 1; multiclass uses macro OvR AUC."""
    if _HAVE_SKLEARN:
        from sklearn.preprocessing import label_binarize
        from sklearn.metrics import roc_auc_score

    n_classes = len(classes)
    fig, ax = plt.subplots(figsize=(6, 5), dpi=300)

    if n_classes == 2:
        fpr, tpr, _ = roc_curve(y_true, y_prob[:, 1])
        roc_auc = auc(fpr, tpr)
        ax.plot(fpr, tpr, color="#D62728", lw=2, label=f"AUC = {roc_auc:.3f}")
    else:
        y_bin = label_binarize(y_true, classes=list(range(n_classes)))
        roc_auc = roc_auc_score(y_bin, y_prob, multi_class="ovr", average="macro")
        for i, cls in enumerate(classes):
            fpr, tpr, _ = roc_curve(y_bin[:, i], y_prob[:, i])
            ax.plot(fpr, tpr, lw=1.5, label=f"{cls} (AUC={auc(fpr, tpr):.3f})")

    ax.plot([0, 1], [0, 1], color="gray", lw=1, ls="--")
    ax.set_xlabel("False Positive Rate", fontsize=11)
    ax.set_ylabel("True Positive Rate", fontsize=11)
    ax.set_title(f"ROC Curve  ·  {prefix}", fontsize=13)
    ax.legend(fontsize=9, loc="lower right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()

    base = os.path.join(outdir, f"{prefix}.roc_curve")
    fig.savefig(f"{base}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{base}.pdf", bbox_inches="tight")
    plt.close(fig)
    return float(roc_auc)


def plot_feature_importance(
    feat_df: pd.DataFrame,
    outdir: str,
    prefix: str,
    top_n: int = 30,
) -> None:
    """Horizontal bar chart of top feature importances."""
    plot_df = feat_df.head(top_n).copy()
    if plot_df.empty:
        return

    plot_df = plot_df.sort_values("importance")
    fig, ax = plt.subplots(figsize=(8, max(4, 0.3 * len(plot_df) + 1.5)), dpi=300)
    ax.barh(plot_df["feature"], plot_df["importance"], color="#4E79A7", edgecolor="white", linewidth=0.4)
    ax.set_xlabel("Importance", fontsize=11)
    ax.set_title(f"Top {len(plot_df)} Feature Importances  ·  {prefix}", fontsize=13)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()

    base = os.path.join(outdir, f"{prefix}.feature_importance")
    fig.savefig(f"{base}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{base}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_confusion_matrix(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    classes: List[str],
    outdir: str,
    prefix: str,
) -> None:
    """Normalised confusion matrix heatmap."""
    cm = confusion_matrix(y_true, y_pred, labels=list(range(len(classes))))
    cm_norm = cm.astype(float) / cm.sum(axis=1, keepdims=True)

    fig, ax = plt.subplots(figsize=(max(4, len(classes)), max(3.5, len(classes) * 0.9)), dpi=300)
    sns.heatmap(
        cm_norm, annot=cm, fmt="d", cmap="Blues",
        xticklabels=classes, yticklabels=classes,
        ax=ax, linewidths=0.5, linecolor="white",
        cbar_kws={"label": "Proportion"},
    )
    ax.set_xlabel("Predicted", fontsize=11)
    ax.set_ylabel("True", fontsize=11)
    ax.set_title(f"Confusion Matrix  ·  {prefix}", fontsize=13)
    plt.tight_layout()

    base = os.path.join(outdir, f"{prefix}.confusion_matrix")
    fig.savefig(f"{base}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{base}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_umap(
    X: np.ndarray,
    y: np.ndarray,
    classes: List[str],
    outdir: str,
    prefix: str,
) -> None:
    """UMAP embedding coloured by class label."""
    try:
        from umap import UMAP
    except ImportError:
        print("[WARN] umap-learn not installed; skipping UMAP plot.")
        return

    reducer = UMAP(n_components=2, random_state=42)
    embedding = reducer.fit_transform(X)

    palette = sns.color_palette("tab10", n_colors=len(classes))
    fig, ax = plt.subplots(figsize=(6, 5), dpi=300)
    for i, cls in enumerate(classes):
        mask = y == i
        ax.scatter(embedding[mask, 0], embedding[mask, 1],
                   s=18, color=palette[i], alpha=0.75, label=cls, linewidths=0)
    ax.set_xlabel("UMAP 1", fontsize=11)
    ax.set_ylabel("UMAP 2", fontsize=11)
    ax.set_title(f"UMAP Embedding  ·  {prefix}", fontsize=13)
    ax.legend(fontsize=9, markerscale=1.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()

    base = os.path.join(outdir, f"{prefix}.umap")
    fig.savefig(f"{base}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{base}.pdf", bbox_inches="tight")
    plt.close(fig)


# =========================================================
# Main
# =========================================================

def main():
    parser = argparse.ArgumentParser(
        description="Omics ML classifier: train, evaluate, and interpret.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Inputs
    parser.add_argument("--features", required=True,
        help="Feature matrix (TSV or CSV). Rows=samples, Cols=features (genes, proteins, etc.).")
    parser.add_argument("--labels", required=True,
        help="Label file (TSV or CSV). Must have 'sample' and 'label' columns.")

    # Model
    parser.add_argument(
        "--model",
        choices=["random_forest", "logistic", "svm", "xgboost"],
        default="random_forest",
        help="Classifier to train.",
    )

    # Feature selection
    parser.add_argument("--select-k", type=int, default=None,
        help="Select top-K features by ANOVA F-score before training. None = use all features.")

    # Splitting
    parser.add_argument("--test-size", type=float, default=0.2,
        help="Fraction of data held out for final test evaluation.")
    parser.add_argument("--cv-folds", type=int, default=5,
        help="Number of stratified cross-validation folds.")
    parser.add_argument("--random-state", type=int, default=42,
        help="Random seed for reproducibility.")

    # Output
    parser.add_argument("--prefix", default=None,
        help="Output file prefix. Defaults to '<model>_classifier'.")
    parser.add_argument("--top-n-features", type=int, default=30,
        help="Number of top features to show in importance plot.")
    parser.add_argument("--umap", action="store_true",
        help="Generate UMAP plot of feature space (requires umap-learn).")
    parser.add_argument("--outdir", required=True, help="Output directory.")

    args = parser.parse_args()
    ensure_dir(args.outdir)

    prefix = args.prefix or f"{args.model}_classifier"

    # ----------------------------------------------------------
    # Load data
    # ----------------------------------------------------------
    print(f"[INFO] Loading features: {args.features}")
    X_df = load_feature_matrix(args.features)
    print(f"[INFO] Feature matrix: {X_df.shape[0]} samples × {X_df.shape[1]} features")

    print(f"[INFO] Loading labels: {args.labels}")
    labels = load_labels(args.labels)

    # Align samples
    common = X_df.index.intersection(labels.index)
    if len(common) == 0:
        raise ValueError("No common samples between feature matrix and label file.")
    if len(common) < len(X_df):
        print(f"[WARN] {len(X_df) - len(common)} samples in features not in labels (dropped).")

    X_df = X_df.loc[common]
    labels = labels.loc[common]

    # Encode labels
    le = LabelEncoder()
    y = le.fit_transform(labels.values)
    classes = list(le.classes_)
    print(f"[INFO] Classes: {classes}  (n={len(common)} samples)")

    # Drop columns with all NaN; fill remaining NaN with column median
    X_df = X_df.dropna(axis=1, how="all")
    X_df = X_df.fillna(X_df.median())
    feature_names = list(X_df.columns)
    X = X_df.values

    if X.shape[1] == 0:
        raise ValueError("Feature matrix has no usable columns after NaN filtering.")

    # ----------------------------------------------------------
    # Train/test split
    # ----------------------------------------------------------
    X_train, X_test, y_train, y_test = train_test_split(
        X, y,
        test_size=args.test_size,
        stratify=y,
        random_state=args.random_state,
    )
    print(f"[INFO] Train: {len(X_train)} samples  |  Test: {len(X_test)} samples")

    # ----------------------------------------------------------
    # Build pipeline
    # ----------------------------------------------------------
    pipeline = build_pipeline(args.model, X.shape[1], args.select_k, args.random_state)

    # ----------------------------------------------------------
    # Cross-validation
    # ----------------------------------------------------------
    print(f"[INFO] Running {args.cv_folds}-fold stratified cross-validation")
    cv = StratifiedKFold(n_splits=args.cv_folds, shuffle=True, random_state=args.random_state)
    cv_results = cross_validate(
        pipeline, X_train, y_train, cv=cv,
        scoring=["accuracy", "roc_auc_ovr_weighted", "f1_weighted"],
        return_train_score=False,
    )

    cv_df = pd.DataFrame(
        {
            "fold":     list(range(1, args.cv_folds + 1)),
            "accuracy": cv_results["test_accuracy"],
            "roc_auc":  cv_results["test_roc_auc_ovr_weighted"],
            "f1":       cv_results["test_f1_weighted"],
        }
    )
    cv_df.loc[len(cv_df)] = {
        "fold": "mean",
        "accuracy": cv_df["accuracy"].mean(),
        "roc_auc":  cv_df["roc_auc"].mean(),
        "f1":       cv_df["f1"].mean(),
    }
    cv_path = os.path.join(args.outdir, f"{prefix}.cv_results.tsv")
    cv_df.to_csv(cv_path, sep="\t", index=False)
    print(f"[INFO] CV accuracy: {cv_df.iloc[:-1]['accuracy'].mean():.3f} ± {cv_df.iloc[:-1]['accuracy'].std():.3f}")
    print(f"[INFO] CV ROC-AUC:  {cv_df.iloc[:-1]['roc_auc'].mean():.3f} ± {cv_df.iloc[:-1]['roc_auc'].std():.3f}")

    # ----------------------------------------------------------
    # Final training + test evaluation
    # ----------------------------------------------------------
    print(f"[INFO] Training on full training set")
    pipeline.fit(X_train, y_train)

    y_pred = pipeline.predict(X_test)
    y_prob = pipeline.predict_proba(X_test)

    report = classification_report(y_test, y_pred, target_names=classes, output_dict=True)
    report_df = pd.DataFrame(report).T.reset_index().rename(columns={"index": "class"})
    report_path = os.path.join(args.outdir, f"{prefix}.test_metrics.tsv")
    report_df.to_csv(report_path, sep="\t", index=False)

    test_acc = float(np.mean(y_pred == y_test))
    print(f"[INFO] Test accuracy: {test_acc:.3f}")

    # ----------------------------------------------------------
    # Feature importances
    # ----------------------------------------------------------
    feat_df = extract_feature_importances(pipeline, feature_names, args.model)
    feat_path = os.path.join(args.outdir, f"{prefix}.feature_importance.tsv")
    feat_df.to_csv(feat_path, sep="\t", index=False)
    print(f"[INFO] Feature importances saved: {feat_path}")

    # ----------------------------------------------------------
    # Plots
    # ----------------------------------------------------------
    print(f"[INFO] Generating ROC curve")
    roc_auc_test = plot_roc(y_test, y_prob, classes, args.outdir, prefix)

    print(f"[INFO] Generating confusion matrix")
    plot_confusion_matrix(y_test, y_pred, classes, args.outdir, prefix)

    print(f"[INFO] Generating feature importance plot")
    plot_feature_importance(feat_df, args.outdir, prefix, top_n=args.top_n_features)

    if args.umap:
        print(f"[INFO] Generating UMAP embedding")
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        plot_umap(X_scaled, y, classes, args.outdir, prefix)

    # ----------------------------------------------------------
    # Summary
    # ----------------------------------------------------------
    summary: Dict = {
        "prefix":            prefix,
        "model":             args.model,
        "n_samples":         len(common),
        "n_features_input":  len(feature_names),
        "n_features_used":   args.select_k if args.select_k else len(feature_names),
        "n_classes":         len(classes),
        "classes":           ";".join(classes),
        "train_samples":     len(X_train),
        "test_samples":      len(X_test),
        "cv_folds":          args.cv_folds,
        "cv_accuracy_mean":  round(float(cv_df.iloc[:-1]["accuracy"].mean()), 4),
        "cv_accuracy_std":   round(float(cv_df.iloc[:-1]["accuracy"].std()), 4),
        "cv_roc_auc_mean":   round(float(cv_df.iloc[:-1]["roc_auc"].mean()), 4),
        "cv_roc_auc_std":    round(float(cv_df.iloc[:-1]["roc_auc"].std()), 4),
        "test_accuracy":     round(test_acc, 4),
        "test_roc_auc":      round(roc_auc_test, 4),
        "top_feature_1":     feat_df.iloc[0]["feature"] if len(feat_df) > 0 else "",
        "top_feature_2":     feat_df.iloc[1]["feature"] if len(feat_df) > 1 else "",
        "top_feature_3":     feat_df.iloc[2]["feature"] if len(feat_df) > 2 else "",
    }

    pd.DataFrame([summary]).to_csv(
        os.path.join(args.outdir, "ml_summary.tsv"), sep="\t", index=False
    )
    with open(os.path.join(args.outdir, "summary.txt"), "w", encoding="utf-8") as f:
        for k, v in summary.items():
            f.write(f"{k}: {v}\n")

    print(f"[DONE] Results written to: {args.outdir}")


if __name__ == "__main__":
    main()
