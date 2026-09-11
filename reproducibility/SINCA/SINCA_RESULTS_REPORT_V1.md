# SINCA PM2.5 null-transport analysis — V1 results

## Frozen design

- Primary cohort: 31 stations, frozen before surrogate computation.
- Primary window: 2018-01-01 to 2025-12-28 (417 complete Monday–Sunday weeks).
- Scientific event: daily PM2.5 > 50 µg m^-3.
- Native contract: daily calendar-month exact-value constrained marginal-spectrum surrogate, then thresholding, weekly aggregation and identical preprocessing.
- Index contract: IAAFT on the weekly residual series.
- Primary statistic: Ferro–Segers extremal index, q = 0.90.
- B = 499 surrogates per contract and station.
- Sensitivity/intervention: month-specific probability-equalized thresholds.
- Input adapter: the complete CSV contains pollutant-specific rows in addition to a consolidated station-day row. The PM2.5 analysis retains only rows with nonmissing `pm25`; within the frozen window this yields at most one PM2.5 row per station-date.
- Primary-input SHA-256: `3c467b1cbcfaba06a885a947232c10e2a2b9f6c8e07661e0db5db1e4776ba607`.

The PM2.5-only extraction reproduces the frozen 2018–2025 station coverage, exceedance counts, and N_eff,p values exactly (up to floating-point roundoff).

## Primary fixed-threshold result

At PM2.5 > 50 µg m^-3:

- Index-level IAAFT rejects in **31/31** stations.
- Native transported-null rejects in **2/31** stations.
- Discordance: **29 index-only**, 0 native-only, 2 both, 0 neither.
- Region-level exact sign-flip test for index-only minus native-only: **p = 0.003906**.
- After Benjamini–Yekutieli correction across stations: **27 index-level discoveries**, **0 native-null discoveries**.

For the primary Ferro–Segers statistic, across stations:

- median observed theta: **0.535**
- median native-null median: **0.557**
- median index-null median: **0.993**
- median index-minus-native null-center gap: **0.422**
- median index/native null-width ratio: **1.120**

Thus the real-data discrepancy is dominated by a strong **upward shift of the index-level null center toward theta = 1**, rather than by universal null-width contraction.

## Probability-equalization intervention

After making monthly exceedance probabilities approximately uniform:

- realized N_eff,p is approximately 12 at every station;
- index-level rejections fall from **31 to 6**;
- native-null rejections change from **2 to 1**;
- discordance becomes **5 index-only**, 0 native-only, 1 both, 25 neither;
- region sign-flip p = **0.0625**;
- no station remains significant after BY correction under either contract.

The median null-center gap drops from **0.422** to **0.097**.
It decreases in **31/31** stations.
The index-level p-value increases in **30/31** stations.

This is the strongest applied mechanism result: seasonal event geometry is not merely associated with the discrepancy; an intervention that removes the seasonal concentration sharply reduces it.

## Prespecified mechanistic prediction that was not supported

The prespecified width-contraction prediction was:

N_eff,p decreases -> SD(index null) / SD(native null) decreases.

The observed fixed-threshold Spearman correlation is **rho = -0.379**, with region-bootstrap 95% interval
**[-0.624, 0.136]**.

This is not supportive of the prespecified prediction and its point estimate has the opposite sign.
Therefore the SINCA application should not be presented as validation of a universal width-contraction mechanism.

Instead, the dominant applied mechanism is **null-center displacement**. The fixed-threshold index null has a median theta near 1, while the native-null center remains close to the observed statistic.

## Robustness across statistics

For Ferro–Segers q = 0.85:

- fixed threshold: index rejects 31/31, native rejects 10/31;
- equalized threshold: index rejects 2/31, native rejects 0/31.

For normalized T_diff:

- fixed threshold: index rejects 4/31, native rejects 1/31;
- equalized threshold: index rejects 5/31, native rejects 2/31.

Hence the magnitude of null-transport distortion is strongly statistic-dependent, consistent with the kernel -> statistic -> decision hierarchy in the manuscript.

## Interpretation for the AOAS manuscript

The applied result is substantially stronger than the previous Amazon-CMIP6 illustration:

1. the cohort and inferential protocol were frozen before surrogate results were inspected;
2. the disagreement is systematic across 8 regions;
3. the native and index contracts yield radically different scientific conclusions for the policy-relevant fixed threshold;
4. probability equalization nearly removes the disagreement;
5. the real-data mechanism refines the synthetic story: the key distortion is not always null narrowing; here it is primarily a shift in the projected extremal-index null.

The manuscript should therefore distinguish:
- synthetic hard-count null contraction as one mechanism,
- real-data null-center displacement as another manifestation of the same failure of null transport,
- and the broader object of interest as the full projected-null discrepancy, not only its variance.
