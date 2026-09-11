# SINCA PM2.5 AOAS — Inferential protocol V1

**Status:** PRESPECIFIED BEFORE SINCA SURROGATE RESULTS  
**Cohort:** `SINCA_PM25_AOAS_FROZEN_COHORT_V1` (31 stations; station-list SHA-256 `59b019cd344b1cdf9c0f6fb90b381ad12d416d72461d94d31d6ce0039ff010e6`)  
**Primary window:** 2018-01-01 through 2025-12-28 (417 complete Monday–Sunday weeks). The final three days of 2025 are not used in the weekly analysis.  
**Pollutant:** daily PM2.5 as supplied in the frozen source parquet.  

## 1. Scientific target

The applied analysis asks whether clustering of high weekly PM2.5 exceedance burden depends materially on whether the null contract is imposed on the native daily PM2.5 process and transported through thresholding/aggregation, or imposed directly on the transformed weekly index.

The primary event is the policy-relevant fixed exceedance

\[
I_{s,t}=\mathbf 1\{PM_{2.5,s,t}>50\}.
\]

The application is not used to claim that one null is universally correct. It is used to determine whether the aggregate-only null gives materially different conclusions from the native transported null in observational environmental data, and whether the discrepancy varies with seasonal event-probability concentration.

## 2. Frozen stations

All 31 stations in `FROZEN_COHORT_V1.csv` are retained. No station may be added or removed based on surrogate p-values, effect directions, null widths, or visual inspection of inferential results. Any technical input failure must be documented and versioned rather than silently altering the cohort.

## 3. Missingness and weekly observable

No stochastic imputation of exceedance events is used.

For each station and calendar month, estimate the observed fixed-threshold exceedance probability

\[
\hat p_{s,m}=\frac{\#\{t:\,m(t)=m,\ X_{s,t}>50\}}{\#\{t:\,m(t)=m,\ X_{s,t}\text{ observed}\}}.
\]

For an observed day, set \(\tilde I_{s,t}=I_{s,t}\). For a missing day, set the event contribution to its fixed month-specific expectation,

\[
\tilde I_{s,t}=\hat p_{s,m(t)}.
\]

The primary weekly burden is

\[
C_{s,w}=\sum_{t\in w}\tilde I_{s,t},
\]

for each of the 417 fixed Monday–Sunday weeks. This rule keeps the weekly time axis regular and makes missing-day handling deterministic and identical for observed and surrogate paths. The supplied cohort has at least 95% of weeks with >=6 observed days by design.

A prespecified sensitivity analysis replaces the expected-value fill by the exposure-standardized count

\[
C^{\mathrm{exp}}_{s,w}=7\,E_{s,w}/n_{s,w},
\]

where \(E_{s,w}\) and \(n_{s,w}\) are the observed exceedance and observed-day counts. This sensitivity is evaluated only when every week has at least one observed day; otherwise it is omitted for that station and the omission is reported.

## 4. Fixed weekly preprocessing

The deterministic annual cycle and smooth secular drift are removed before the clustering statistic. For each weekly series, ordinary least squares is fit to

\[
1,t,t^2,t^3,\{\sin(2\pi k d/365.2425),\cos(2\pi k d/365.2425)\}_{k=1}^3,
\]

where \(t\in[-1,1]\) spans the 417 weeks and \(d\) is the week midpoint in days from 2018-01-01. The residual series \(Z_{s,w}\) is the analyzed weekly index. The same regression design is re-fit separately to every native surrogate weekly series. The index-resolution surrogate is applied to the observed residual weekly series.

Three harmonics are fixed in advance to remove the strong deterministic seasonal cycle without changing the event definition.

## 5. Primary statistic and test direction

The primary station-level statistic is the clipped Ferro–Segers intervals estimator at \(q=0.90\), implemented identically to the synthetic benchmark. Smaller \(\hat\theta\) means stronger clustering of high weekly burden anomalies.

For a null ensemble \(\{\hat\theta_b^\star\}_{b=1}^B\), the lower-tail Monte Carlo p-value is

\[
p=\frac{1+\#\{b:\hat\theta_b^\star\le \hat\theta_{obs}\}}{B+1}.
\]

The primary rejection rule is strict \(p<0.05\).

A secondary, tie-robust order statistic is

\[
T_{\mathrm{diff}}=\frac{\frac1{n-1}\sum_{w=1}^{n-1}(Z_{w+1}-Z_w)^2}{\operatorname{Var}(Z)},
\]

with lower values indicating smoother/persistent organization. It uses the same lower-tail Monte Carlo rule.

## 6. Native-resolution surrogate contract

The native contract acts on the daily PM2.5 series before thresholding and weekly aggregation.

1. Compute observed calendar-month means and standard deviations from nonmissing daily values.
2. Standardize observed values by calendar month. Missing daily positions are set to zero **only for construction of the complete standardized sequence used in the FFT**; zero is the month-specific standardized mean.
3. The target Fourier amplitudes are those of this completed standardized daily sequence.
4. Exact-value groups are calendar month over the full 2018–2025 window. Surrogate initialization permutes observed raw PM2.5 values only among observed positions within each calendar-month group. Missing positions remain fixed at the observed month mean and never enter the exact-value multiset.
5. Each iteration standardizes the current series with the observed month-specific mean/sd, resets missing standardized positions to zero, replaces Fourier amplitudes by the target amplitudes while retaining current phases, inverse transforms, and rank-remaps scores at **observed positions only** within each calendar month to the exact sorted observed raw values. Missing positions are reset to the month mean.
6. Twelve iterations are used, matching the synthetic native contract.
7. After native surrogation, exceedances are evaluated only on originally observed days. Missing-day event contributions remain the fixed \(\hat p_{s,m}\) from Section 3. Weekly aggregation, harmonic/cubic residualization, and the statistic are then applied.

This contract preserves the exact observed PM2.5 value multiset in each calendar month and therefore preserves the exact fixed-threshold exceedance count among observed days in each calendar month.

## 7. Index-resolution surrogate contract

The index contract is IAAFT applied to the observed preprocessed weekly residual series \(Z_s\).

Each surrogate starts from a random permutation, replaces Fourier amplitudes by those of the centered observed weekly residual, inverse transforms, and rank-remaps to the exact sorted observed residual values. Thirty iterations are used. The empirical marginal is preserved exactly and the weekly Fourier-amplitude spectrum approximately.

## 8. Monte Carlo size and seeds

Primary analysis: \(B=499\) surrogates per station and contract. With strict \(p<0.05\), the attainable rejection level is 24/500 = 0.048.

All random streams use deterministic station-ID/contract namespaces from base seed `20260911`. Native and index ensembles use distinct namespaces. Seed derivation is recorded in the output manifest.

## 9. Primary station-level outputs

For every station and contract report:

- observed \(\hat\theta\);
- Monte Carlo p-value and rejection indicator;
- null median, mean, and standard deviation;
- \(\Delta_{nat}=\operatorname{median}(\theta^\star_{nat})-\theta_{obs}\);
- \(\Delta_{idx}=\operatorname{median}(\theta^\star_{idx})-\theta_{obs}\);
- center gap \(G=\operatorname{median}(\theta^\star_{idx})-\operatorname{median}(\theta^\star_{nat})\);
- width ratio \(R=sd(\theta^\star_{idx})/sd(\theta^\star_{nat})\);
- paired outcome: both / neither / index-only / native-only;
- native and index spectral-fidelity diagnostics.

## 10. Across-station primary summaries

The application has two prespecified across-station summaries.

**A. Paired asymmetry.** Count index-only and native-only stations. Because stations within a region are dependent, the primary cluster-level directional assessment uses an exact region sign-flip test over the eight represented regions, with each region contributing its signed difference `#index-only - #native-only`. The ordinary station-level exact binomial/McNemar-style calculation is reported only as a descriptive sensitivity.

**B. Seasonal concentration mechanism.** The prespecified continuous mechanism outcome is

\[
\log R_s=\log\{sd(\theta^\star_{idx,s})/sd(\theta^\star_{nat,s})\}.
\]

Its association with the frozen \(N_{eff,p,s}\) is summarized by Spearman correlation. A 95% region-cluster bootstrap interval is reported using 2,000 bootstrap resamples of the eight regions. The directional synthetic prediction is that lower \(N_{eff,p}\) is associated with stronger index-null contraction, hence a positive association between \(N_{eff,p}\) and \(\log R\).

No station is selected or excluded based on this association.

## 11. Multiplicity

Station-level native and index p-values are each adjusted across the 31 frozen stations using Benjamini–Yekutieli FDR at q=0.05, which does not require independence. Raw p-values remain visible. FDR correction is secondary to the paired null-transport comparison.

## 12. Probability-equalization sensitivity

The regulatory/fixed threshold `50` remains the primary scientific event.

A mechanistic secondary analysis defines, for station \(s\),

\[
\bar p_s=\frac1{12}\sum_{m=1}^{12}\hat p_{s,m},
\]

and month-specific empirical thresholds

\[
u_{s,m}=Q_{1-\bar p_s}(X_{s,t}:m(t)=m).
\]

The entire native/index analysis is repeated using these seasonally relative thresholds. This branch tests whether equalizing monthly event probability reduces the null-transport discrepancy; it is not presented as a replacement for the regulatory event definition.

## 13. Prespecified sensitivities

1. Exposure-standardized weekly count instead of expected-value fill, where technically defined.
2. Ferro–Segers q=0.85 in addition to the primary q=0.90.
3. Secondary normalized \(T_{diff}\) statistic.
4. Probability-equalized thresholds (Section 12).

No alternative is promoted to primary based on results.

## 14. Technical gates

Before inference the implementation must verify for every station:

- exact frozen station list and source hashes;
- 417 weekly blocks;
- exact preservation of observed calendar-month raw-value multisets by native surrogates (numerical tolerance 1e-12 after sorting);
- exact preservation of observed calendar-month fixed-threshold exceedance counts among observed days by the native surrogates;
- finite weekly residuals and finite primary statistics;
- recorded spectral errors for both surrogate contracts.

A failed gate stops that station's run and is reported; it does not silently change the cohort.

## 15. Interpretation rule

Index-only rejection is interpreted as evidence that the aggregate-level marginal-and-spectrum null yields a clustering conclusion not reproduced by the native transported contract. Native non-rejection is not interpreted as evidence of absence of additional temporal organization. The application is an inferential-contract comparison, not a causal attribution analysis of air pollution dynamics.
