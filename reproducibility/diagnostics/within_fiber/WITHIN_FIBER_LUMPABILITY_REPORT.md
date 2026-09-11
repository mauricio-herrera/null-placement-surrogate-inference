# Within-fiber strong-lumpability diagnostic

## Design

For each Gaussian AR(1) trajectory `x`, a fiber mate `x_tilde` was constructed by permuting the **exact amplitudes within each calendar-month × threshold-indicator class**. Thus, by construction, `x` and `x_tilde` have:

- identical monthly wet/dry indicator sequences;
- identical annual count sequence `O(x)=O(x_tilde)`;
- identical calendar-month empirical value multisets, means, and variances;
- generally different standardized monthly Fourier-amplitude targets.

The implemented native surrogate kernel depends on the observed trajectory through the calendar-month value multisets and the standardized monthly target Fourier amplitudes. Because the first component is exactly preserved by the fiber construction, the only intended change in the native surrogate contract is the hidden monthly spectral target.

For every `x`, three independent native-surrogate ensembles were generated with `B=249`: two from the same `x` (Monte Carlo baseline) and one from `x_tilde` (within-fiber comparison). The test used 30 trajectories at each `phi = 0.2, 0.5, 0.8`, for 90 fiber pairs per seasonal profile, with `lambda=0` and `lambda=4`.

All construction checks passed exactly: maximum indicator mismatch = 0, maximum annual-count difference = 0, and maximum calendar-month multiset difference = 0.

## Main result

The native pushforward distribution depends strongly on microstate information that is absent from the annual count observable. This directly contradicts the empirical implication of strong lumpability with respect to the threshold-and-annual-count operator.

For the Ferro–Segers projection:

- `lambda=0`: median scaled Wasserstein distance was **0.175** across fiber mates versus **0.086** between independent same-x ensembles (ratio **1.93**). The fiber distance exceeded the same-x baseline in **73.3%** of trajectories (one-sided sign-test p = **5.45e-06**).
- `lambda=4`: median scaled Wasserstein distance was **0.335** versus **0.090** (ratio **3.31**). The fiber distance exceeded baseline in **91.1%** of trajectories (p = **6.92e-17**).
- A two-sample KS comparison rejected equality of the two Ferro–Segers null distributions for **70.0%** of `lambda=4` fiber pairs, compared with only **1.1%** of same-x Monte Carlo controls.

The effect is much larger for second-order/order-sensitive projections:

| projection | lambda | median fiber/same-x distance ratio | fiber KS p<.05 | same-x KS p<.05 |
|---|---:|---:|---:|---:|
| lag-1 correlation | 0 | 5.27 | 78.9% | 5.6% |
| lag-1 correlation | 4 | 4.55 | 68.9% | 3.3% |
| raw Tdiff | 0 | 15.62 | 100.0% | 2.2% |
| raw Tdiff | 4 | 15.28 | 100.0% | 7.8% |
| low-frequency share | 0 | 5.36 | 86.7% | 3.3% |
| low-frequency share | 4 | 3.97 | 75.6% | 2.2% |

For raw `Tdiff`, **100%** of fiber-pair null distributions were distinguishable by the KS test for both seasonal profiles. The median fiber/same-x Wasserstein-distance ratio was approximately **15.6** at `lambda=0` and **15.3** at `lambda=4`.

## Mechanistic interpretation

This construction preserves substantially more than the annual count: it preserves the **entire monthly indicator chronology** and every calendar-month empirical amplitude multiset. Yet the native-surrogate pushforward changes when the amplitudes are rearranged inside wet/dry classes. For the implemented native surrogate, the only intended contract component altered by this construction is the standardized monthly target spectrum.

Therefore the annual count `O(x)` is not sufficient to determine the native transported null. Hidden monthly amplitude organization, erased by thresholding, remains active in the native surrogate contract and changes the distribution after thresholding and aggregation.

The target-spectrum separation is larger under seasonal concentration (median normalized distance 0.473 at `lambda=4` versus 0.349 at `lambda=0`). Across trajectories, target-spectrum distance is strongly associated with the transported-null distance for `Tdiff` (Spearman rho about 0.90 at `lambda=4` in the post-analysis diagnostic), but not monotonically with Ferro–Segers. This is consistent with a hierarchy in which kernel-level differences can be large while their projection onto a particular extremal statistic is mechanism-dependent.

## What this establishes, and what it does not

The experiment provides strong numerical evidence against **strong lumpability of the implemented native randomization kernel with respect to the threshold + annual-count operator**. A strict mathematical proof would require an analytic counterexample or exact kernel calculation rather than Monte Carlo estimates.

The result does **not** say that no reduced null can be defined. When strong lumpability fails, a distribution-dependent transported null can still be defined by averaging over unresolved microstates conditional on the observed aggregate. It also does not imply that every statistic or every decision rule will reveal the kernel difference: the decision-level rejection discordance is much smaller than the distributional differences, which is exactly why kernel-, statistic-, and decision-level compatibility should be separated in the manuscript.

## Implication for manuscript reframing

This result supports a `null transport` formulation with two distinct questions:

1. **Lumpability / existence:** does `O_# K_X(x,·)` depend on `x` only through `O(x)`? The within-fiber experiment says no for the implemented native kernel.
2. **Surrogate substitution:** when a unique quotient null does not exist (or when a model-dependent transported null is used), how far is an index-level IAAFT kernel from the transported native null, especially after hard thresholding and seasonal concentration?

This makes the current index/native discordance a downstream consequence of an information-loss problem rather than merely a comparison of two arbitrary surrogate procedures.
