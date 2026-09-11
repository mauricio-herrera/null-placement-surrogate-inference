# Near-zero operator-ablation and mechanism diagnostic

## Design
Gaussian AR(1) branch, lambda=4, cubic-detrended annual soft counts
S_y(tau)=sum_m expit((u_m-X_ym)/tau), approaching the hard annual count as tau -> 0.

Near-zero pooled diagnostic: n=180 trajectory-condition replicates per tau (two independent batches,
phi in {0.2,0.5,0.8}, B=99 index/native/oracle Monte Carlo surrogates per trajectory).
This is a mechanism diagnostic rather than the final confirmatory B=249 experiment.

## Near-zero gradient
|   tau |   n |   oracle_reject_rate |   native_reject_rate |   index_reject_rate |   index_only_rate |   native_only_rate |   index_only_count |   native_only_count |   paired_exact_p |   median_sd_ratio_index_oracle |   median_sd_ratio_native_oracle |   median_index_spectral_error |   median_corr_soft_hard |   median_nrmse_soft_hard |   median_raw_unique_soft |   median_null_mean_index_minus_oracle |
|------:|----:|---------------------:|---------------------:|--------------------:|------------------:|-------------------:|-------------------:|--------------------:|-----------------:|-------------------------------:|--------------------------------:|------------------------------:|------------------------:|-------------------------:|-------------------------:|--------------------------------------:|
| 0.1   | 180 |            0.0388889 |            0.0222222 |           0.0444444 |         0.0277778 |         0.00555556 |                  5 |                   1 |      0.21875     |                       0.846155 |                        1.00312  |                   0.000249937 |                0.981649 |                0.192905  |                    251   |                            0.00862764 |
| 0.05  | 180 |            0.0333333 |            0.0222222 |           0.0333333 |         0.0166667 |         0.00555556 |                  3 |                   1 |      0.625       |                       0.764146 |                        0.994145 |                   0.000451306 |                0.990526 |                0.138386  |                    251   |                            0.0159346  |
| 0.025 | 180 |            0.0388889 |            0.0277778 |           0.0777778 |         0.05      |         0          |                  9 |                   0 |      0.00390625  |                       0.691045 |                        0.975201 |                   0.000890876 |                0.995393 |                0.0967594 |                    251   |                            0.0292066  |
| 0.01  | 180 |            0.0444444 |            0.0388889 |           0.122222  |         0.0833333 |         0          |                 15 |                   0 |      6.10352e-05 |                       0.583132 |                        0.986368 |                   0.00193545  |                0.998315 |                0.0584375 |                    214   |                            0.0510765  |
| 0.005 | 180 |            0.0444444 |            0.0388889 |           0.127778  |         0.0888889 |         0          |                 16 |                   0 |      3.05176e-05 |                       0.555895 |                        0.980374 |                   0.00271418  |                0.99918  |                0.04063   |                    152   |                            0.0595791  |
| 0.002 | 180 |            0.0444444 |            0.0388889 |           0.15      |         0.111111  |         0          |                 20 |                   0 |      1.90735e-06 |                       0.537791 |                        0.979629 |                   0.00355653  |                0.999683 |                0.0252362 |                     84   |                            0.0640624  |
| 0.001 | 180 |            0.0444444 |            0.0444444 |           0.155556  |         0.111111  |         0          |                 20 |                   0 |      1.90735e-06 |                       0.538189 |                        0.9863   |                   0.00378119  |                0.999853 |                0.0171693 |                     53.5 |                            0.0694463  |
| 0     | 180 |            0.0388889 |            0.0388889 |           0.161111  |         0.127778  |         0.00555556 |                 23 |                   1 |      2.98023e-06 |                       0.532295 |                        0.990282 |                   0.00423902  |                1        |                0         |                      8   |                            0.0637029  |

## Main interpretation
The transition is gradual rather than discontinuous at tau=0. Index-resolution rejection increases from
4.44% at tau=0.1 to 7.78% at tau=0.025, 12.22% at tau=0.01, 15.56% at tau=0.001, and 16.11% at the hard limit.
Native and oracle rejection remain near nominal. Simultaneously, the median index/oracle null-SD ratio falls
from 0.846 to 0.532 and the IAAFT spectral-error floor rises from 2.50e-4 to 4.24e-3.

Across tau-level medians, spectral error and index rejection co-move strongly, while spectral error and the
index/oracle null-SD ratio move oppositely. Within a fixed tau, however, trajectory-level spectral error is only
weakly related to null-SD contraction. Therefore spectral error is best treated as a diagnostic marker of the
same increasing constraint rigidity, not yet as a demonstrated causal mediator.

## IAAFT convergence diagnostic
| condition    |   iters |   median_spec_error |   median_q90_spec_error |   median_sd_theta |   median_theta |
|:-------------|--------:|--------------------:|------------------------:|------------------:|---------------:|
| hard_count   |      10 |         0.0031886   |             0.00404316  |         0.0775868 |              1 |
| hard_count   |      30 |         0.00285942  |             0.00393243  |         0.0779187 |              1 |
| hard_count   |     100 |         0.00285942  |             0.00393243  |         0.0779187 |              1 |
| hard_count   |     300 |         0.00285942  |             0.00393243  |         0.0779187 |              1 |
| mean_cubic   |      10 |         0.000225985 |             0.000304559 |         0.0781592 |              1 |
| mean_cubic   |      30 |         0.00018705  |             0.000249757 |         0.0769267 |              1 |
| mean_cubic   |     100 |         0.00018705  |             0.000249757 |         0.0769267 |              1 |
| mean_cubic   |     300 |         0.00018705  |             0.000249757 |         0.0769267 |              1 |
| soft_tau0.01 |      10 |         0.00176501  |             0.00240438  |         0.0791019 |              1 |
| soft_tau0.01 |      30 |         0.00148447  |             0.00212046  |         0.0752514 |              1 |
| soft_tau0.01 |     100 |         0.00147565  |             0.00212046  |         0.0752514 |              1 |
| soft_tau0.01 |     300 |         0.00147565  |             0.00212046  |         0.0752514 |              1 |

Thirty iterations are already at the observed error floor: hard-count median spectral error is 0.002859 at
30, 100, and 300 iterations; soft tau=0.01 changes only from 0.001484 to 0.001476. The phenomenon is therefore
not explained by stopping IAAFT too early.

## Extremal-spacing geometry
Pooled over two independent small diagnostics (48 observed trajectories; B=99 null replicates per trajectory):
|   tau | null   |   mean_adj |   sd_adj |   mean_clusters |   sd_clusters |   mean_gap_sd |   sd_theta |   mean_theta |
|------:|:-------|-----------:|---------:|----------------:|--------------:|--------------:|-----------:|-------------:|
| 0     | index  |    2.34343 |  1.34135 |         22.6566 |       1.34135 |       8.58505 |  0.06865   |     0.972366 |
| 0     | native |    2.89394 |  1.5282  |         22.1061 |       1.5282  |       9.79023 |  0.131849  |     0.901874 |
| 0     | oracle |    2.88384 |  1.56211 |         22.1162 |       1.56211 |       9.94636 |  0.131551  |     0.891602 |
| 0.001 | index  |    2.39394 |  1.3287  |         22.6061 |       1.3287  |       8.59183 |  0.0689514 |     0.970694 |
| 0.001 | native |    2.87374 |  1.52786 |         22.1263 |       1.52786 |       9.76699 |  0.129085  |     0.903733 |
| 0.001 | oracle |    2.87374 |  1.53815 |         22.1263 |       1.53815 |       9.92473 |  0.131722  |     0.894413 |
| 0.01  | index  |    2.44949 |  1.3239  |         22.5505 |       1.3239  |       8.59466 |  0.0620066 |     0.974176 |
| 0.01  | native |    2.70707 |  1.51029 |         22.2929 |       1.51029 |       9.50888 |  0.11846   |     0.915993 |
| 0.01  | oracle |    2.81818 |  1.52668 |         22.1818 |       1.52668 |       9.62491 |  0.120688  |     0.910576 |
| 0.1   | index  |    2.4697  |  1.34118 |         22.5303 |       1.34118 |       8.59578 |  0.0649903 |     0.971147 |
| 0.1   | native |    2.44444 |  1.45437 |         22.5556 |       1.45437 |       8.79406 |  0.0810025 |     0.960138 |
| 0.1   | oracle |    2.42929 |  1.42969 |         22.5707 |       1.42969 |       8.79278 |  0.0819115 |     0.95884  |

At the hard limit the index-resolution null has fewer adjacent upper-tail exceedance pairs (median mean 2.34)
than oracle/native (~2.88/2.89), smaller inter-exceedance-gap variability (8.59 vs 9.95 oracle), a much higher
mean Ferro-Segers theta (0.972 vs 0.892 oracle), and about half the theta dispersion (0.069 vs 0.132 oracle).
At tau=0.1, these discrepancies are small. This directly identifies the projection through which the lower-tail
test over-rejects: the annual IAAFT null regularizes upper-tail spacing and suppresses clustering relative to the
generative/native null.

## Global phase-dispersion diagnostic
| condition   |   median_phase_resultant |   q90_phase_resultant |   mean_phase_resultant |   median_spec_error |   sd_theta |    sd_ntv |
|:------------|-------------------------:|----------------------:|-----------------------:|--------------------:|-----------:|----------:|
| mean        |                0.0682473 |              0.121625 |              0.0733591 |         0.000183097 |  0.0827348 | 0.0182938 |
| tau0.1      |                0.0692689 |              0.121618 |              0.0726814 |         0.000244713 |  0.0760899 | 0.0174315 |
| tau0.025    |                0.0688554 |              0.122394 |              0.0719225 |         0.000714202 |  0.0785958 | 0.0184401 |
| tau0.01     |                0.0687852 |              0.12323  |              0.0725615 |         0.00135595  |  0.0763226 | 0.0184573 |
| tau0.001    |                0.0689281 |              0.121674 |              0.0728467 |         0.00248325  |  0.0820323 | 0.0191276 |
| hard        |                0.0665099 |              0.122093 |              0.0717479 |         0.00266862  |  0.0795596 | 0.0196015 |

The median circular phase resultant stays close to ~0.068 across continuous, near-hard and hard conditions.
Thus the effect is not a simple collapse or phase locking of the entire surrogate ensemble; it is much more
specific to extremal/rank geometry.

## Current mechanism statement
The evidence supports a mechanism of post-transformation constraint rigidity. Thresholding plus annual
aggregation is a many-to-one map that compresses continuous native microstates into an annual observable with
near-lattice / low-effective-support geometry. Imposing the exact annual sample marginal and the annual
sample-spectrum constraint after that map restricts the admissible orderings of upper-rank values much more
strongly than imposing the spectral constraint at native resolution and then pushing it through the
threshold-aggregation operator. The index-resolution IAAFT null consequently under-disperses extremal spacing,
suppresses adjacent exceedance clustering, shifts the Ferro-Segers null toward theta=1, and yields excess
lower-tail rejection. Seasonal concentration further reduces effective event geometry and amplifies the effect.

This is stronger than the original descriptive explanation, but it is not yet a mathematical proof that
spectral-error magnitude itself causes the distortion.

## Theoretical reframing
Let O:X->Y be the observation/coarse-graining operator and K_X a native randomization kernel.
Strong null-transport compatibility is the one-step lumpability condition
O_# K_X(x,.) = O_# K_X(x',.) whenever O(x)=O(x').
Then a quotient kernel Kbar on Y exists. Exact index-level transport additionally requires K_Y=Kbar.

If strong lumpability fails, a P0-dependent averaged transport can be defined by
Kbar_P0(y,B)=E_P0[K_X(X,O^-1B) | O(X)=y],
which is distinct from an arbitrary index surrogate such as IAAFT.

The paper can therefore separate:
(1) lumpability/transport defect,
(2) mismatch between the chosen index surrogate and the quotient/averaged transported null,
(3) statistic-projected defect, and
(4) decision discordance.

The current experiments provide strong evidence for (2)-(4); a dedicated within-fibre construction would be
needed to demonstrate strong-lumpability failure directly.
