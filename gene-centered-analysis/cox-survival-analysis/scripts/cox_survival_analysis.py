#!/usr/bin/env python3
"""
Cox proportional hazards survival analysis.
Pure NumPy implementation of Newton-Raphson Cox partial likelihood optimization.
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
try:
    from scipy import stats as _scipy_stats
    _HAVE_SCIPY = True
except ImportError:
    _scipy_stats = None
    _HAVE_SCIPY = False

import warnings


# ---------------------------------------------------------------------------
# Pure-numpy scipy.stats fallbacks
# ---------------------------------------------------------------------------
class _FallbackStats:
    """Minimal scipy.stats replacements using pure numpy."""

    class norm:
        @staticmethod
        def cdf(x):
            """Normal CDF via error function."""
            return 0.5 * (1.0 + np.vectorize(lambda v: _erf(v / np.sqrt(2)))(np.asarray(x, float)))

    class chi2:
        @staticmethod
        def cdf(x, df):
            """Chi-squared CDF via regularised lower incomplete gamma."""
            return _gammainc_lower(df / 2.0, np.asarray(x, float) / 2.0)

    @staticmethod
    def spearmanr(a, b):
        """Spearman rank correlation."""
        def _rank(arr):
            order = np.argsort(arr)
            ranks = np.empty_like(order, dtype=float)
            ranks[order] = np.arange(len(arr)) + 1.0
            return ranks
        n = len(a)
        ra, rb = _rank(np.asarray(a, float)), _rank(np.asarray(b, float))
        ra -= ra.mean(); rb -= rb.mean()
        denom = np.sqrt((ra**2).sum() * (rb**2).sum()) + 1e-15
        r = (ra * rb).sum() / denom
        t = r * np.sqrt(max(n - 2, 1)) / np.sqrt(max(1 - r**2, 1e-15))
        # two-tailed p from t-distribution approx
        p = _t_sf_approx(abs(t), n - 2) * 2
        return r, min(p, 1.0)


def _erf(x):
    """Abramowitz & Stegun approximation of erf."""
    sign = 1.0 if x >= 0 else -1.0
    x = abs(x)
    t = 1.0 / (1.0 + 0.3275911 * x)
    y = 1.0 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * np.exp(-x * x)
    return sign * y


def _gammainc_lower(a, x):
    """Regularised lower incomplete gamma via series expansion (vectorised)."""
    x = np.asarray(x, float)
    scalar = x.ndim == 0
    x = np.atleast_1d(x)
    result = np.zeros_like(x)
    for idx in range(len(x)):
        xi = x[idx]
        if xi <= 0:
            result[idx] = 0.0
            continue
        # Series: sum_{n=0}^{inf} x^n / (a*(a+1)*...*(a+n))
        term = xi**a * np.exp(-xi) / max(a, 1e-30)
        s = term
        for n_iter in range(1, 200):
            term *= xi / (a + n_iter)
            s += term
            if abs(term) < 1e-10 * abs(s):
                break
        # normalise by gamma(a) -- use Stirling for large a
        log_gamma_a = _log_gamma(a)
        result[idx] = min(s * np.exp(-log_gamma_a), 1.0)
    return result[0] if scalar else result


def _log_gamma(x):
    """log(Gamma(x)) via Lanczos approximation."""
    if x < 0.5:
        return np.log(np.pi / np.sin(np.pi * x)) - _log_gamma(1 - x)
    x -= 1
    a = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
         771.32342877765313, -176.61502916214059, 12.507343278686905,
         -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]
    t = x + 7.5
    return 0.5 * np.log(2 * np.pi) + (x + 0.5) * np.log(t) - t + np.log(
        a[0] + sum(a[i] / (x + i) for i in range(1, 9)))


def _t_sf_approx(t, df):
    """Approximate survival function of t-distribution."""
    if df <= 0:
        return 0.5
    x = df / (df + t**2)
    # Regularised incomplete beta I_x(df/2, 0.5)
    return _betainc(df / 2.0, 0.5, x) / 2.0 if t > 0 else 0.5


def _betainc(a, b, x):
    """Regularised incomplete beta via continued fraction."""
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    lbeta = _log_gamma(a) + _log_gamma(b) - _log_gamma(a + b)
    front = np.exp(a * np.log(x) + b * np.log(1 - x) - lbeta) / a
    # Lentz continued fraction
    cf = 1.0
    d = 1.0 - (a + b) * x / (a + 1)
    if abs(d) < 1e-30:
        d = 1e-30
    d = 1.0 / d
    cf = d
    for m in range(1, 100):
        for sign in (1, -1):
            if sign == 1:
                num = m * (b - m) * x / ((a + 2*m - 1) * (a + 2*m))
            else:
                num = -(a + m) * (a + b + m) * x / ((a + 2*m) * (a + 2*m + 1))
            d = 1.0 + num * d
            if abs(d) < 1e-30:
                d = 1e-30
            d = 1.0 / d
            cf_prev = cf
            cf *= d
            if abs(cf - cf_prev) < 1e-10 * abs(cf):
                break
    return front * cf


stats = _scipy_stats if _HAVE_SCIPY else _FallbackStats()

warnings.filterwarnings('ignore')

# Set matplotlib backend
plt.switch_backend('Agg')


class CoxSurvivalAnalysis:
    """Cox proportional hazards regression."""

    def __init__(self, args):
        self.args = args
        self.data = None
        self.time = None
        self.event = None
        self.X = None
        self.covariate_names = None
        self.sample_ids = None
        self.univariate_results = None
        self.multivariate_results = None
        self.ph_test_results = None
        self.risk_scores = None

    def log_msg(self, msg):
        """Print timestamped message."""
        print(f"[INFO] {msg}")

    def load_data(self):
        """Load and prepare data."""
        self.log_msg(f"Loading data from {self.args.input}")

        self.data = pd.read_csv(self.args.input, sep='\t', index_col=0)
        self.sample_ids = self.data.index.values

        # Extract time and event
        if self.args.time_col not in self.data.columns:
            raise ValueError(f"Time column '{self.args.time_col}' not found")
        if self.args.event_col not in self.data.columns:
            raise ValueError(f"Event column '{self.args.event_col}' not found")

        self.time = self.data[self.args.time_col].values.astype(np.float32)
        self.event = (self.data[self.args.event_col].values > 0).astype(int)

        # Validate
        if (self.time < 0).any():
            raise ValueError("Time values must be non-negative")

        n_events = self.event.sum()
        n_censored = (self.event == 0).sum()
        self.log_msg(f"Events: {n_events}, Censored: {n_censored}")

    def prepare_covariates(self):
        """Prepare covariate matrix."""
        self.log_msg(f"Preparing covariates: {self.args.covariates}")

        cov_cols = [c.strip() for c in self.args.covariates.split(',')]

        # Parse categorical/continuous specifications
        categorical_cols = set()
        continuous_cols = set()

        if self.args.categorical:
            categorical_cols = set(c.strip() for c in self.args.categorical.split(','))

        if self.args.continuous:
            continuous_cols = set(c.strip() for c in self.args.continuous.split(','))

        # Auto-detect if not specified
        for col in cov_cols:
            if col in categorical_cols:
                continue
            if col in continuous_cols:
                continue

            # Check if numeric
            try:
                pd.to_numeric(self.data[col])
                continuous_cols.add(col)
            except (ValueError, TypeError):
                categorical_cols.add(col)

        # Prepare matrix
        X_list = []
        covariate_names = []

        for col in cov_cols:
            if col in categorical_cols:
                # One-hot encode (drop first category as reference)
                dummies = pd.get_dummies(self.data[col], prefix=col, drop_first=True)
                X_list.append(dummies.values.astype(np.float32))
                covariate_names.extend(dummies.columns.tolist())
            else:
                # Continuous: standardize
                col_data = pd.to_numeric(self.data[col], errors='coerce').values
                # Handle missing
                mask = ~np.isnan(col_data)
                col_data[~mask] = np.nanmean(col_data)

                mean = col_data.mean()
                std = col_data.std() + 1e-10
                col_standardized = (col_data - mean) / std

                X_list.append(col_standardized.reshape(-1, 1).astype(np.float32))
                covariate_names.append(col)

        self.X = np.hstack(X_list) if X_list else np.array([]).reshape(len(self.sample_ids), 0)
        self.covariate_names = np.array(covariate_names)

        # Handle interaction terms
        if self.args.interaction_terms:
            self._add_interaction_terms(cov_cols, continuous_cols, categorical_cols)

        self.log_msg(f"Covariate matrix: {self.X.shape[0]} samples × {self.X.shape[1]} covariates")

        # Complete case analysis
        mask = ~(np.isnan(self.time) | np.isnan(self.event) | np.isnan(self.X).any(axis=1))
        if not mask.all():
            n_dropped = (~mask).sum()
            self.log_msg(f"Complete case analysis: dropping {n_dropped} samples with missing values")

            self.X = self.X[mask]
            self.time = self.time[mask]
            self.event = self.event[mask]
            self.sample_ids = self.sample_ids[mask]

        # Sort by time
        sort_idx = np.argsort(self.time)
        self.X = self.X[sort_idx]
        self.time = self.time[sort_idx]
        self.event = self.event[sort_idx]
        self.sample_ids = self.sample_ids[sort_idx]

    def _add_interaction_terms(self, cov_cols, continuous_cols, categorical_cols):
        """Add interaction terms."""
        if not self.args.interaction_terms:
            return

        pairs = [p.strip() for p in self.args.interaction_terms.split(',')]

        for pair in pairs:
            col1, col2 = pair.split(':')
            col1, col2 = col1.strip(), col2.strip()

            # Find indices in X
            idx1 = np.where(self.covariate_names == col1)[0]
            idx2 = np.where(self.covariate_names == col2)[0]

            if len(idx1) > 0 and len(idx2) > 0:
                interaction = self.X[:, idx1[0]] * self.X[:, idx2[0]]
                self.X = np.hstack([self.X, interaction.reshape(-1, 1).astype(np.float32)])
                self.covariate_names = np.append(self.covariate_names, f"{col1}:{col2}")

    def cox_model(self, X):
        """Fit Cox proportional hazards model via Newton-Raphson."""
        n_samples, n_covariates = X.shape

        if n_covariates == 0:
            return np.array([]), np.array([]), np.array([])

        # Initialize coefficients
        beta = np.zeros(n_covariates)

        # Newton-Raphson
        max_iter = 20
        for iteration in range(max_iter):
            # Risk sets: for each event time, which samples are at risk
            risk_sets = {}
            for i in range(n_samples):
                if self.event[i] == 1:
                    t = self.time[i]
                    if t not in risk_sets:
                        risk_sets[t] = []
                    risk_sets[t].append(i)

            # Compute gradient and Hessian
            gradient = np.zeros(n_covariates)
            hessian = np.zeros((n_covariates, n_covariates))

            # For each event
            for t in sorted(risk_sets.keys()):
                event_indices = risk_sets[t]

                for i in event_indices:
                    # Risk set at time t
                    risk_set = np.where(self.time >= t)[0]

                    # Compute weighted sums
                    X_risk = X[risk_set]
                    exp_Xbeta = np.exp(X_risk @ beta)
                    sum_exp = exp_Xbeta.sum()

                    if sum_exp == 0:
                        continue

                    weights = exp_Xbeta / sum_exp

                    # Expected X
                    E_X = weights @ X_risk

                    # Gradient (score function)
                    gradient += X[i] - E_X

                    # Hessian (outer product)
                    for j in range(len(risk_set)):
                        centered_X = X_risk[j] - E_X
                        hessian += weights[j] * np.outer(centered_X, centered_X)

            # Newton step
            try:
                delta_beta = np.linalg.solve(hessian, gradient)
            except np.linalg.LinAlgError:
                # Singular matrix
                self.log_msg("Warning: Singular Hessian matrix")
                break

            # Check convergence
            if np.linalg.norm(gradient) < 1e-6 or np.max(np.abs(delta_beta)) < 1e-6:
                self.log_msg(f"Cox model converged in {iteration + 1} iterations")
                break

            beta = beta + delta_beta

        # Compute variance-covariance matrix
        try:
            var_beta = np.linalg.inv(hessian)
        except np.linalg.LinAlgError:
            var_beta = np.diag(np.ones(n_covariates) * 1e10)

        se_beta = np.sqrt(np.diag(var_beta))

        return beta, se_beta, var_beta

    def univariate_analysis(self):
        """Fit separate Cox model per covariate."""
        self.log_msg("Running univariate Cox analysis...")

        results = []

        for idx, cov_name in enumerate(self.covariate_names):
            X_uni = self.X[:, idx:idx+1]
            beta, se_beta, _ = self.cox_model(X_uni)

            if len(beta) == 0:
                continue

            beta = beta[0]
            se = se_beta[0]

            # Hazard ratio and CI
            hr = np.exp(beta)
            ci_lower = np.exp(beta - 1.96 * se)
            ci_upper = np.exp(beta + 1.96 * se)

            # Z-test
            z = beta / (se + 1e-10)
            pval = 2 * (1 - stats.norm.cdf(np.abs(z)))

            results.append({
                'covariate': cov_name,
                'HR': hr,
                'CI_lower': ci_lower,
                'CI_upper': ci_upper,
                'pvalue': pval,
                'significant': pval < self.args.alpha
            })

        self.univariate_results = pd.DataFrame(results)

    def multivariate_analysis(self):
        """Fit Cox model with all covariates."""
        self.log_msg("Running multivariate Cox analysis...")

        beta, se_beta, var_beta = self.cox_model(self.X)

        results = []

        for idx, cov_name in enumerate(self.covariate_names):
            beta_i = beta[idx] if idx < len(beta) else 0
            se_i = se_beta[idx] if idx < len(se_beta) else 1e10

            # Hazard ratio and CI
            hr = np.exp(beta_i)
            ci_lower = np.exp(beta_i - 1.96 * se_i)
            ci_upper = np.exp(beta_i + 1.96 * se_i)

            # Z-test
            z = beta_i / (se_i + 1e-10)
            pval = 2 * (1 - stats.norm.cdf(np.abs(z)))

            results.append({
                'covariate': cov_name,
                'HR': hr,
                'CI_lower': ci_lower,
                'CI_upper': ci_upper,
                'pvalue': pval,
                'significant': pval < self.args.alpha
            })

        self.multivariate_results = pd.DataFrame(results)
        self.beta_multi = beta
        self.var_beta_multi = var_beta

    def ph_assumption_test(self):
        """Test proportional hazards assumption via Schoenfeld residuals."""
        self.log_msg("Testing proportional hazards assumption...")

        results = []

        for idx, cov_name in enumerate(self.covariate_names):
            # Schoenfeld residuals
            residuals = []

            for i in np.where(self.event == 1)[0]:
                t = self.time[i]
                risk_set = np.where(self.time >= t)[0]

                X_risk = self.X[risk_set]
                exp_Xbeta = np.exp(X_risk @ self.beta_multi)
                weights = exp_Xbeta / exp_Xbeta.sum()

                E_X_j = (weights * X_risk[:, idx]).sum()
                r_ij = self.X[i, idx] - E_X_j
                residuals.append(r_ij)

            residuals = np.array(residuals)

            # Correlation with log(time)
            log_time = np.log(self.time[self.event == 1] + 1e-10)
            corr, pval = stats.spearmanr(residuals, log_time)

            results.append({
                'covariate': cov_name,
                'correlation': corr,
                'pvalue': pval,
                'ph_violated': pval < 0.05
            })

            # Store for plotting
            if idx == 0:
                self.schoenfeld_residuals = {}
            self.schoenfeld_residuals[cov_name] = (residuals, log_time, corr, pval)

        self.ph_test_results = pd.DataFrame(results)

    def compute_risk_scores(self):
        """Compute linear predictor and risk groups."""
        self.log_msg("Computing risk scores...")

        lp = self.X @ self.beta_multi
        median_lp = np.median(lp)
        risk_group = np.where(lp < median_lp, 'low', 'high')

        self.risk_scores = pd.DataFrame({
            'sample_id': self.sample_ids,
            'linear_predictor': lp,
            'risk_group': risk_group
        })

    def plot_forest_plot(self):
        """Plot forest plot of hazard ratios."""
        self.log_msg("Plotting forest plot...")

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Univariate
        if self.univariate_results is not None:
            uni_data = self.univariate_results.sort_values('covariate')
            y_pos = np.arange(len(uni_data))

            for i, row in uni_data.iterrows():
                color = 'red' if row['significant'] else 'blue'
                ax1.errorbar(row['HR'], i, xerr=[[row['HR'] - row['CI_lower']], [row['CI_upper'] - row['HR']]],
                           fmt='o', color=color, markersize=6, elinewidth=1.5, capsize=3)
                ax1.text(row['HR'] * 1.1, i, f"{row['pvalue']:.2e}", va='center', fontsize=8)

            ax1.axvline(1, color='black', linestyle='--', linewidth=1)
            ax1.set_yticks(y_pos)
            ax1.set_yticklabels(uni_data['covariate'])
            ax1.set_xlabel('Hazard Ratio (log scale)')
            ax1.set_title('Univariate Cox Analysis')
            ax1.set_xscale('log')
            ax1.grid(True, alpha=0.3)

        # Multivariate
        if self.multivariate_results is not None:
            multi_data = self.multivariate_results.sort_values('covariate')
            y_pos = np.arange(len(multi_data))

            for i, row in multi_data.iterrows():
                color = 'red' if row['significant'] else 'blue'
                ax2.errorbar(row['HR'], i, xerr=[[row['HR'] - row['CI_lower']], [row['CI_upper'] - row['HR']]],
                           fmt='o', color=color, markersize=6, elinewidth=1.5, capsize=3)
                ax2.text(row['HR'] * 1.1, i, f"{row['pvalue']:.2e}", va='center', fontsize=8)

            ax2.axvline(1, color='black', linestyle='--', linewidth=1)
            ax2.set_yticks(y_pos)
            ax2.set_yticklabels(multi_data['covariate'])
            ax2.set_xlabel('Hazard Ratio (log scale)')
            ax2.set_title('Multivariate Cox Analysis')
            ax2.set_xscale('log')
            ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        outfile = os.path.join(self.args.outdir, 'forest_plot.png')
        plt.savefig(outfile, dpi=300, bbox_inches='tight')
        self.log_msg(f"Saved {outfile}")
        plt.close()

    def plot_schoenfeld_residuals(self):
        """Plot Schoenfeld residuals vs log(time)."""
        self.log_msg("Plotting Schoenfeld residuals...")

        # Select significant covariates
        sig_covs = self.ph_test_results[self.ph_test_results['pvalue'] < 0.05]['covariate'].values[:6]

        if len(sig_covs) == 0:
            self.log_msg("No significant covariates for PH residual plots")
            return

        n_plots = len(sig_covs)
        n_cols = min(3, n_plots)
        n_rows = (n_plots + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
        axes = axes.ravel() if n_plots > 1 else [axes]

        for plot_idx, cov_name in enumerate(sig_covs):
            if cov_name not in self.schoenfeld_residuals:
                continue

            residuals, log_time, corr, pval = self.schoenfeld_residuals[cov_name]

            ax = axes[plot_idx]
            ax.scatter(log_time, residuals, alpha=0.6, s=30)
            ax.axhline(0, color='red', linestyle='--', linewidth=1)
            ax.set_xlabel('log(time)')
            ax.set_ylabel('Schoenfeld Residuals')
            ax.set_title(f"{cov_name}\n(Spearman corr={corr:.3f}, p={pval:.2e})")
            ax.grid(True, alpha=0.3)

        # Hide unused axes
        for idx in range(n_plots, len(axes)):
            axes[idx].set_visible(False)

        plt.tight_layout()
        outfile = os.path.join(self.args.outdir, 'schoenfeld_residuals.png')
        plt.savefig(outfile, dpi=300, bbox_inches='tight')
        self.log_msg(f"Saved {outfile}")
        plt.close()

    def plot_survival_curves(self):
        """Plot Kaplan-Meier curves by risk group."""
        self.log_msg("Plotting Kaplan-Meier survival curves...")

        fig, ax = plt.subplots(figsize=(10, 6))

        for risk_group in ['low', 'high']:
            mask = self.risk_scores['risk_group'] == risk_group
            time_group = self.time[mask]
            event_group = self.event[mask]

            # KM estimator
            sorted_idx = np.argsort(time_group)
            time_sorted = time_group[sorted_idx]
            event_sorted = event_group[sorted_idx]

            unique_times, idx = np.unique(time_sorted, return_index=True)
            events_at_time = np.zeros_like(unique_times, dtype=int)
            at_risk = np.zeros_like(unique_times, dtype=int)

            n_at_risk = len(time_sorted)
            for i, t in enumerate(unique_times):
                mask_t = time_sorted >= t
                at_risk[i] = mask_t.sum()
                events_at_time[i] = event_sorted[time_sorted == t].sum()

            # Survival probability
            survival_prob = np.cumprod(1 - events_at_time / (at_risk + 1e-10))
            survival_prob = np.insert(survival_prob, 0, 1)
            unique_times = np.insert(unique_times, 0, 0)

            ax.step(unique_times, survival_prob, where='post', label=f'{risk_group.capitalize()} risk', linewidth=2)

        ax.set_xlabel('Time')
        ax.set_ylabel('Survival probability')
        ax.set_title('Kaplan-Meier Survival Curves by Risk Group')
        ax.set_ylim([0, 1.05])
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        outfile = os.path.join(self.args.outdir, 'survival_curves.png')
        plt.savefig(outfile, dpi=300, bbox_inches='tight')
        self.log_msg(f"Saved {outfile}")
        plt.close()

    def save_outputs(self):
        """Save results to TSV files."""
        self.log_msg("Saving output files...")

        if self.univariate_results is not None:
            outfile = os.path.join(self.args.outdir, 'univariate_results.tsv')
            self.univariate_results.to_csv(outfile, sep='\t', index=False)
            self.log_msg(f"Saved {outfile}")

        if self.multivariate_results is not None:
            outfile = os.path.join(self.args.outdir, 'multivariate_results.tsv')
            self.multivariate_results.to_csv(outfile, sep='\t', index=False)
            self.log_msg(f"Saved {outfile}")

        if self.ph_test_results is not None:
            outfile = os.path.join(self.args.outdir, 'ph_test.tsv')
            self.ph_test_results.to_csv(outfile, sep='\t', index=False)
            self.log_msg(f"Saved {outfile}")

        if self.risk_scores is not None:
            outfile = os.path.join(self.args.outdir, 'risk_scores.tsv')
            self.risk_scores.to_csv(outfile, sep='\t', index=False)
            self.log_msg(f"Saved {outfile}")

    def run(self):
        """Execute full analysis."""
        self.log_msg("Starting Cox survival analysis...")

        self.load_data()
        self.prepare_covariates()

        if self.args.mode in ['univariate', 'both']:
            self.univariate_analysis()

        if self.args.mode in ['multivariate', 'both']:
            self.multivariate_analysis()
            self.ph_assumption_test()
            self.compute_risk_scores()

        # Plots
        if self.args.mode in ['univariate', 'both']:
            self.plot_forest_plot()

        if self.args.mode in ['multivariate', 'both']:
            if hasattr(self, 'schoenfeld_residuals'):
                self.plot_schoenfeld_residuals()
            self.plot_survival_curves()

        # Outputs
        self.save_outputs()

        self.log_msg("Analysis complete!")


def main():
    parser = argparse.ArgumentParser(
        description='Cox proportional hazards survival analysis'
    )
    parser.add_argument('--input', required=True, help='Input TSV with clinical/molecular data')
    parser.add_argument('--time-col', required=True, help='Time-to-event column')
    parser.add_argument('--event-col', required=True, help='Event indicator column')
    parser.add_argument('--covariates', required=True, help='Comma-separated covariate columns')
    parser.add_argument('--categorical', help='Comma-separated categorical columns')
    parser.add_argument('--continuous', help='Comma-separated continuous columns')
    parser.add_argument('--mode', choices=['univariate', 'multivariate', 'both'],
                       default='both', help='Analysis mode')
    parser.add_argument('--reference-categories', help='Reference categories (col:ref,col2:ref2)')
    parser.add_argument('--interaction-terms', help='Interaction terms (col1:col2,col3:col4)')
    parser.add_argument('--strata', help='Stratification variable')
    parser.add_argument('--alpha', type=float, default=0.05, help='Significance level')
    parser.add_argument('--outdir', required=True, help='Output directory')

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    analysis = CoxSurvivalAnalysis(args)
    analysis.run()


if __name__ == '__main__':
    main()
