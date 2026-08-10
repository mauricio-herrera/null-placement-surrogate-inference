from __future__ import annotations

"""
Paper A — Stage 4 power under prespecified injected alternatives v0.50

Frozen after Stage-3 confirmatory success.

Scientific target
-----------------
Quantify the power of the native-resolution monthly gate in an eight-unit
synthetic design analogous to the model-level Amazon analysis. Power is
reported separately for two injected mechanisms; it is not treated as a
universal function of Delta-theta.

Primary eight-unit gate
-----------------------
For each synthetic unit:
    delta_i = median(theta_native_null_i) - theta_obs_i.

For a cohort of eight units:
    1) at least 6/8 units must have delta_i > 0; and
    2) a synchronized-surrogate cohort statistic must have p < 0.05.

The cohort statistic is the mean centered extremal-index deficit across units.
For surrogate b:
    T_b = mean_i[median(theta_null_i) - theta_null_{i,b}],
and T_obs = mean_i delta_i.

The Monte Carlo aggregate p-value is
    p = (1 + #{T_b >= T_obs}) / (B + 1).

Frozen power design
-------------------
- 8 independent model-level units per cohort.
- Each cohort samples 8 of the 9 family x phi generator cells without
  replacement, at the frozen strong seasonal concentration lambda=4.
- 251 annual observations from 252 simulated monthly years.
- native null only, block=50 years, q=0.90, B=500.
- 150 independent cohorts per design point.
- one shared baseline/null point plus four levels for each of two mechanisms.

Mechanisms
----------
A. persistent_regime_alignment:
   rank-preserving within-calendar-month reordering toward a common persistent
   annual latent regime (AR(1) rho=0.85). Exact monthly marginals and threshold
   counts are preserved.

B. history_feedback:
   iterative rank-copula perturbation in which recent threshold events lower the
   rank score of subsequent observations through a normalized exponential
   24-month history kernel (tau=6 months). Exact calendar-month marginals and
   threshold counts are preserved.

Run confirmatory Stage 4:
    python paperA_stage4_power_v050.py --jobs 4 \
        --out paperA_stage4_power_v050_run1

Resume:
    python paperA_stage4_power_v050.py --jobs 4 --resume \
        --out paperA_stage4_power_v050_run1

Smoke:
    python paperA_stage4_power_v050.py --smoke --jobs 1 \
        --out paperA_stage4_smoke_v050

Only --jobs, --out, --resume and --smoke are exposed. The scientific grid is
not editable from the command line.
"""

import argparse
import hashlib
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import norm

VERSION = "v0.50"
SEED_NAMESPACE = "paperA_stage4_v050"

FAMILIES = ("gaussian", "periodic", "tinnov")
PHIS = (0.2, 0.5, 0.8)
GENERATOR_CELLS = tuple((f, p) for f in FAMILIES for p in PHIS)
LAMBDA = 4.0

# One common null point and two mechanism-specific curves.
DESIGNS = (
    ("baseline", 0.0),
    ("persistent_regime_alignment", 0.2),
    ("persistent_regime_alignment", 0.4),
    ("persistent_regime_alignment", 0.6),
    ("persistent_regime_alignment", 0.8),
    ("history_feedback", 0.25),
    ("history_feedback", 0.50),
    ("history_feedback", 0.75),
    ("history_feedback", 1.00),
)

COHORT_SIZE = 8
CONFIRM_COHORTS = 150
CONFIRM_B = 500
CHUNK_SIZE = 5
SMOKE_COHORTS = 2
SMOKE_B = 20

NU_T = 5
SIM_YEARS = 252
OUTPUT_YEARS = 251
MONTHS = SIM_YEARS * 12
BLOCK_YEARS = 50
Q_EXTREME = 0.90
K_DRY = 3.0
PERIODIC_AMP = 0.65
ALPHA = 0.05
MONTHLY_IAAFT_ITER = 12
AMAZON_DELTA_REFERENCE = 0.085009

REGIME_RHO = 0.85
HISTORY_HORIZON = 24
HISTORY_TAU = 6.0
HISTORY_ITERATIONS = 4


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="paperA_stage4_power_v050_run1")
    p.add_argument("--jobs", type=int, default=max(1, min(4, os.cpu_count() or 1)))
    p.add_argument("--resume", action="store_true")
    p.add_argument("--smoke", action="store_true")
    return p.parse_args()


def seed_from(*parts) -> int:
    h = hashlib.sha256("|".join(map(str, parts)).encode()).digest()
    return int.from_bytes(h[:8], "little") % (2**32)


def logistic(x):
    return 1.0 / (1.0 + np.exp(-x))


def seasonal_probabilities(lam: float = LAMBDA, K: float = K_DRY, phase_month: int = 7):
    months = np.arange(12)
    seasonal = np.cos(2 * np.pi * (months - phase_month) / 12.0)

    def f(a):
        return logistic(a + lam * seasonal).sum() - K

    a = brentq(f, -30.0, 30.0)
    probs = logistic(a + lam * seasonal)
    neff_p = probs.sum() ** 2 / np.sum(probs**2)
    v = probs * (1 - probs)
    neff_var = v.sum() ** 2 / np.sum(v**2)
    return probs, a, neff_p, neff_var


@lru_cache(maxsize=None)
def periodic_sd(phi: float, amp: float = PERIODIC_AMP, phase_month: int = 0):
    m = np.arange(12)
    sig = np.exp(amp * np.cos(2 * np.pi * (m - phase_month) / 12.0))
    v = 1.0
    sd = np.empty(12, dtype=float)
    for k in range(20000):
        mm = k % 12
        v = phi * phi * v + sig[mm] ** 2
        if k >= 20000 - 12:
            sd[mm] = math.sqrt(v)
    return sig.copy(), sd.copy()


def simulate_gaussian(phi: float, rng: np.random.Generator, n: int = MONTHS):
    x = np.empty(n)
    x[0] = rng.normal()
    s = math.sqrt(1 - phi * phi)
    for t in range(1, n):
        x[t] = phi * x[t - 1] + s * rng.normal()
    return x


def simulate_periodic(phi: float, rng: np.random.Generator, n: int = MONTHS):
    sig, sd = periodic_sd(phi)
    burn = 2400
    xx = np.empty(n + burn)
    xx[0] = 0.0
    for t in range(1, len(xx)):
        mm = t % 12
        xx[t] = phi * xx[t - 1] + sig[mm] * rng.normal()
    return xx[burn:], sd


def simulate_tinnov(phi: float, rng: np.random.Generator, n: int = MONTHS):
    scale = math.sqrt((NU_T - 2) / NU_T)
    burn = 2000
    x = np.empty(n + burn)
    x[0] = 0.0
    for t in range(1, len(x)):
        eps = rng.standard_t(NU_T) * scale
        x[t] = phi * x[t - 1] + eps
    return x[burn:]


@lru_cache(maxsize=None)
def tinnov_reference_series(phi: float):
    rng = np.random.default_rng(seed_from("tinnov_reference_stage4", phi))
    return simulate_tinnov(phi, rng, n=300000)


def thresholds_for(family: str, phi: float, probs: np.ndarray):
    if family == "gaussian":
        return norm.ppf(probs)
    if family == "periodic":
        _, sd = periodic_sd(phi)
        return sd * norm.ppf(probs)
    if family == "tinnov":
        return np.quantile(tinnov_reference_series(float(phi)), probs)
    raise ValueError(family)


def threshold_calibration_error(family: str, phi: float, thresholds: np.ndarray, probs: np.ndarray):
    if family == "gaussian":
        return float(np.max(np.abs(norm.cdf(thresholds) - probs)))
    if family == "periodic":
        _, state_sd = periodic_sd(phi)
        err = float(np.max(np.abs(norm.cdf(thresholds / state_sd) - probs)))
        if err > 1e-10:
            raise RuntimeError(f"Periodic threshold calibration failed: {err:.3e}")
        return err
    return float("nan")


def simulate_family(family: str, phi: float, rng: np.random.Generator):
    if family == "gaussian":
        return simulate_gaussian(phi, rng)
    if family == "periodic":
        return simulate_periodic(phi, rng)[0]
    if family == "tinnov":
        return simulate_tinnov(phi, rng)
    raise ValueError(family)


def indicators(x: np.ndarray, thresholds: np.ndarray):
    month = np.arange(len(x)) % 12
    return (x < thresholds[month]).astype(np.int8)


def annual_count(ind: np.ndarray):
    seg = ind[: OUTPUT_YEARS * 12]
    return seg.reshape(OUTPUT_YEARS, 12).sum(axis=1).astype(float)


def detrend_cubic(x: np.ndarray):
    x = np.asarray(x, float)
    t = np.linspace(-1, 1, len(x))
    coef = np.polyfit(t, x, 3)
    return x - np.polyval(coef, t)


def theta_intervals(x: np.ndarray, q: float = Q_EXTREME):
    x = np.asarray(x, float)
    if not np.all(np.isfinite(x)):
        return np.nan
    u = np.quantile(x, q)
    idx = np.flatnonzero(x > u)
    if len(idx) < 3:
        return np.nan
    T = np.diff(idx).astype(float)
    nt = len(T)
    if np.max(T) <= 2:
        den = nt * np.sum(T * T)
        num = 2 * (np.sum(T) ** 2)
    else:
        S = T - 1
        den = nt * np.sum(S * (S - 1))
        num = 2 * (np.sum(S) ** 2)
    if den <= 0:
        return 1.0
    return float(np.clip(num / den, 0.0, 1.0))


def monthly_groups(n: int, block_years: int = BLOCK_YEARS):
    year = np.arange(n) // 12
    month = np.arange(n) % 12
    block = year // block_years
    gid = block * 12 + month
    return [np.flatnonzero(gid == g) for g in np.unique(gid)]


def constrained_monthly_iaaft(x: np.ndarray, rng: np.random.Generator,
                              block_years: int = BLOCK_YEARS,
                              max_iter: int = MONTHLY_IAAFT_ITER):
    x = np.asarray(x, float)
    month = np.arange(len(x)) % 12
    mu = np.array([x[month == m].mean() for m in range(12)])
    sd = np.array([x[month == m].std(ddof=0) for m in range(12)])
    sd = np.where(sd > 0, sd, 1.0)
    z = (x - mu[month]) / sd[month]

    groups = monthly_groups(len(x), block_years)
    sorted_vals = [np.sort(x[idx]) for idx in groups]
    target_amp = np.abs(np.fft.rfft(z - z.mean()))

    y = x.copy()
    for idx in groups:
        y[idx] = rng.permutation(y[idx])

    last = np.inf
    err = np.inf
    for _ in range(max_iter):
        yz = (y - mu[month]) / sd[month]
        fy = np.fft.rfft(yz - yz.mean())
        zz = np.fft.irfft(target_amp * np.exp(1j * np.angle(fy)), n=len(x)) + z.mean()

        yn = np.empty_like(y)
        for idx, vals in zip(groups, sorted_vals):
            order = np.argsort(zz[idx], kind="mergesort")
            tmp = np.empty(len(idx))
            tmp[order] = vals
            yn[idx] = tmp

        ynz = (yn - mu[month]) / sd[month]
        amp = np.abs(np.fft.rfft(ynz - ynz.mean()))
        err = np.mean((amp - target_amp) ** 2) / (np.mean(target_amp**2) + 1e-15)
        y = yn
        if abs(last - err) < 1e-8:
            break
        last = err
    return y, float(err)


def mc_lower_p(null_thetas: np.ndarray, obs: float):
    valid = np.isfinite(null_thetas)
    if not np.isfinite(obs) or valid.sum() == 0:
        return np.nan
    return float((1 + np.sum(null_thetas[valid] <= obs)) / (valid.sum() + 1))


def rank_normal_scores(vals: np.ndarray):
    vals = np.asarray(vals, float)
    order = np.argsort(vals, kind="mergesort")
    ranks = np.empty(len(vals), dtype=float)
    ranks[order] = np.arange(len(vals), dtype=float)
    u = (ranks + 0.5) / len(vals)
    return norm.ppf(u)


def rank_remap(original_vals: np.ndarray, target_score: np.ndarray):
    vals = np.sort(np.asarray(original_vals, float))
    order = np.argsort(np.asarray(target_score, float), kind="mergesort")
    out = np.empty_like(vals)
    out[order] = vals
    return out


def inject_persistent_regime_alignment(x: np.ndarray, kappa: float, rng: np.random.Generator):
    if kappa == 0:
        return x.copy()
    latent = np.empty(SIM_YEARS)
    latent[0] = rng.normal()
    s = math.sqrt(1 - REGIME_RHO**2)
    for y in range(1, SIM_YEARS):
        latent[y] = REGIME_RHO * latent[y-1] + s * rng.normal()
    latent = (latent - latent.mean()) / (latent.std(ddof=0) + 1e-15)

    out = x.copy()
    for m in range(12):
        idx = np.arange(m, len(x), 12)
        vals = x[idx]
        base = rank_normal_scores(vals)
        target = (1.0 - kappa) * base + kappa * latent
        out[idx] = rank_remap(vals, target)
    return out


def history_score(ind: np.ndarray):
    lags = np.arange(1, HISTORY_HORIZON + 1)
    w = np.exp(-lags / HISTORY_TAU)
    w /= w.sum()
    h = np.zeros(len(ind), dtype=float)
    for lag, ww in zip(lags, w):
        h[lag:] += ww * ind[:-lag]
    return h


def inject_history_feedback(x: np.ndarray, thresholds: np.ndarray, kappa: float):
    if kappa == 0:
        return x.copy()
    original = x.copy()
    base_scores = {}
    for m in range(12):
        idx = np.arange(m, len(x), 12)
        base_scores[m] = rank_normal_scores(original[idx])

    y = original.copy()
    for _ in range(HISTORY_ITERATIONS):
        ind = indicators(y, thresholds)
        h = history_score(ind)
        yn = y.copy()
        for m in range(12):
            idx = np.arange(m, len(x), 12)
            hm = h[idx]
            hz = (hm - hm.mean()) / (hm.std(ddof=0) + 1e-15)
            target = base_scores[m] - kappa * hz
            yn[idx] = rank_remap(original[idx], target)
        y = yn
    return y


def inject_mechanism(x: np.ndarray, thresholds: np.ndarray, mechanism: str,
                     kappa: float, rng: np.random.Generator):
    if mechanism == "baseline" or kappa == 0:
        return x.copy()
    if mechanism == "persistent_regime_alignment":
        return inject_persistent_regime_alignment(x, kappa, rng)
    if mechanism == "history_feedback":
        return inject_history_feedback(x, thresholds, kappa)
    raise ValueError(mechanism)


def monthly_multiset_exact(x0: np.ndarray, x1: np.ndarray):
    for m in range(12):
        a = np.sort(x0[m::12])
        b = np.sort(x1[m::12])
        if not np.array_equal(a, b):
            return False
    return True


def monthly_event_count_exact(x0: np.ndarray, x1: np.ndarray, thresholds: np.ndarray):
    i0 = indicators(x0, thresholds)
    i1 = indicators(x1, thresholds)
    for m in range(12):
        if int(i0[m::12].sum()) != int(i1[m::12].sum()):
            return False
    return True


def lagcorr(x: np.ndarray, lag: int):
    x = np.asarray(x, float)
    if lag <= 0 or len(x) <= lag:
        return np.nan
    a, b = x[:-lag], x[lag:]
    if np.std(a) == 0 or np.std(b) == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0,1])


def unit_analysis(family: str, phi: float, mechanism: str, kappa: float,
                  cohort_rep: int, unit_index: int, B: int):
    probs, intercept, neff_p, neff_var = seasonal_probabilities(LAMBDA)
    thresholds = thresholds_for(family, phi, probs)
    cal_err = threshold_calibration_error(family, phi, thresholds, probs)

    base_rng = np.random.default_rng(seed_from(SEED_NAMESPACE, "base", cohort_rep, unit_index, family, phi))
    x0 = simulate_family(family, phi, base_rng)

    inj_rng = np.random.default_rng(seed_from(SEED_NAMESPACE, "inject", mechanism, kappa, cohort_rep, unit_index))
    x = inject_mechanism(x0, thresholds, mechanism, kappa, inj_rng)

    multiset_exact = monthly_multiset_exact(x0, x)
    event_count_exact = monthly_event_count_exact(x0, x, thresholds)
    if not multiset_exact or not event_count_exact:
        raise RuntimeError("Injected mechanism failed exact calendar-month marginal/event-count conservation")

    ind = indicators(x, thresholds)
    count = annual_count(ind)
    residual = detrend_cubic(count)
    theta_obs = theta_intervals(residual)

    null_theta = np.empty(B)
    spectral_err = np.empty(B)
    exact_total = np.empty(B, dtype=bool)
    exact_group = np.empty(B, dtype=bool)
    groups = monthly_groups(len(x), BLOCK_YEARS)
    orig_group_counts = np.array([ind[idx].sum() for idx in groups], dtype=int)

    null_seed_mech = "baseline" if kappa == 0 else mechanism
    null_seed_kappa = 0.0 if kappa == 0 else kappa
    for b in range(B):
        rng = np.random.default_rng(seed_from(SEED_NAMESPACE, "native", null_seed_mech,
                                              null_seed_kappa, cohort_rep, unit_index, b))
        sx, err = constrained_monthly_iaaft(x, rng, BLOCK_YEARS)
        sind = indicators(sx, thresholds)
        scount = annual_count(sind)
        null_theta[b] = theta_intervals(detrend_cubic(scount))
        spectral_err[b] = err
        exact_total[b] = int(sind.sum()) == int(ind.sum())
        new_group_counts = np.array([sind[idx].sum() for idx in groups], dtype=int)
        exact_group[b] = bool(np.array_equal(orig_group_counts, new_group_counts))

    center = float(np.nanmedian(null_theta))
    delta = center - theta_obs
    p_native = mc_lower_p(null_theta, theta_obs)

    unit_row = {
        "mechanism": mechanism,
        "kappa": kappa,
        "cohort_replicate": cohort_rep,
        "unit_index": unit_index,
        "family": family,
        "phi": phi,
        "lambda": LAMBDA,
        "neff_p": float(neff_p),
        "neff_var": float(neff_var),
        "threshold_calibration_max_error": cal_err,
        "theta_obs": float(theta_obs),
        "theta_native_null_median": center,
        "delta_theta": float(delta),
        "p_native": float(p_native),
        "native_reject_005": bool(p_native < ALPHA),
        "delta_positive": bool(delta > 0),
        "native_null_valid_fraction": float(np.mean(np.isfinite(null_theta))),
        "native_spectral_error_median": float(np.nanmedian(spectral_err)),
        "native_count_exact_fraction": float(np.mean(exact_total)),
        "native_group_count_exact_fraction": float(np.mean(exact_group)),
        "injection_monthly_multiset_exact": bool(multiset_exact),
        "injection_monthly_event_count_exact": bool(event_count_exact),
        "annual_count_mean": float(np.mean(count)),
        "annual_count_sd": float(np.std(count, ddof=0)),
        "annual_count_ac1": lagcorr(count, 1),
        "monthly_indicator_ac1": lagcorr(ind, 1),
        "monthly_indicator_ac12": lagcorr(ind, 12),
    }
    return unit_row, null_theta


def cohort_analysis(mechanism: str, kappa: float, cohort_rep: int, B: int):
    # Same balanced random cohort composition for all mechanisms/kappa at a given replicate.
    crng = np.random.default_rng(seed_from(SEED_NAMESPACE, "composition", cohort_rep))
    chosen_idx = crng.choice(len(GENERATOR_CELLS), size=COHORT_SIZE, replace=False)
    cells = [GENERATOR_CELLS[i] for i in chosen_idx]

    unit_rows = []
    null_matrix = np.empty((COHORT_SIZE, B))
    for u, (family, phi) in enumerate(cells):
        row, null = unit_analysis(family, phi, mechanism, kappa, cohort_rep, u, B)
        unit_rows.append(row)
        null_matrix[u] = null

    obs_delta = np.array([r["delta_theta"] for r in unit_rows], float)
    centers = np.array([r["theta_native_null_median"] for r in unit_rows], float)
    theta_obs = np.array([r["theta_obs"] for r in unit_rows], float)

    T_obs = float(np.mean(centers - theta_obs))
    T_null = np.mean(centers[:, None] - null_matrix, axis=0)
    valid = np.isfinite(T_null)
    aggregate_p = float((1 + np.sum(T_null[valid] >= T_obs)) / (valid.sum() + 1))
    positive_units = int(np.sum(obs_delta > 0))
    aggregate_reject = bool(aggregate_p < ALPHA)
    sign_gate = bool(positive_units >= 6)
    gate_pass = bool(sign_gate and aggregate_reject)

    cohort_row = {
        "mechanism": mechanism,
        "kappa": kappa,
        "cohort_replicate": cohort_rep,
        "n_units": COHORT_SIZE,
        "positive_units": positive_units,
        "sign_gate_6of8": sign_gate,
        "aggregate_stat_mean_delta": T_obs,
        "aggregate_p": aggregate_p,
        "aggregate_reject_005": aggregate_reject,
        "full_gate_pass": gate_pass,
        "mean_delta_theta": float(np.mean(obs_delta)),
        "median_delta_theta": float(np.median(obs_delta)),
        "mean_theta_obs": float(np.mean(theta_obs)),
        "mean_theta_null_center": float(np.mean(centers)),
        "unit_native_reject_fraction": float(np.mean([r["native_reject_005"] for r in unit_rows])),
        "unit_positive_fraction": float(np.mean(obs_delta > 0)),
        "native_null_valid_fraction_min": float(min(r["native_null_valid_fraction"] for r in unit_rows)),
        "native_spectral_error_median": float(np.median([r["native_spectral_error_median"] for r in unit_rows])),
        "native_count_exact_fraction_min": float(min(r["native_count_exact_fraction"] for r in unit_rows)),
        "native_group_count_exact_fraction_min": float(min(r["native_group_count_exact_fraction"] for r in unit_rows)),
        "injection_multiset_exact_all": bool(all(r["injection_monthly_multiset_exact"] for r in unit_rows)),
        "injection_event_count_exact_all": bool(all(r["injection_monthly_event_count_exact"] for r in unit_rows)),
    }
    return cohort_row, unit_rows


def run_chunk(task):
    mechanism, kappa, start, stop, B = task
    cohort_rows = []
    unit_rows = []
    for rep in range(start, stop):
        c, u = cohort_analysis(mechanism, kappa, rep, B)
        cohort_rows.append(c)
        unit_rows.extend(u)
    return cohort_rows, unit_rows


def wilson_interval(k: int, n: int, alpha: float = 0.05):
    if n == 0:
        return np.nan, np.nan
    z = norm.ppf(1 - alpha/2)
    phat = k/n
    den = 1 + z*z/n
    center = (phat + z*z/(2*n))/den
    half = z * math.sqrt(phat*(1-phat)/n + z*z/(4*n*n))/den
    return center-half, center+half


def summarize(cohorts: pd.DataFrame, units: pd.DataFrame, B: int):
    rows = []
    for (mech, kap), g in cohorts.groupby(["mechanism", "kappa"], sort=False):
        n = len(g)
        passes = int(g.full_gate_pass.sum())
        agg = int(g.aggregate_reject_005.sum())
        sign = int(g.sign_gate_6of8.sum())
        lo, hi = wilson_interval(passes, n)
        ug = units[(units.mechanism == mech) & (units.kappa == kap)]
        rows.append({
            "mechanism": mech,
            "kappa": kap,
            "n_cohorts": n,
            "full_gate_power": passes/n,
            "full_gate_power_ci_low": lo,
            "full_gate_power_ci_high": hi,
            "aggregate_p_power": agg/n,
            "sign_gate_rate": sign/n,
            "single_unit_native_rejection_rate": float(ug.native_reject_005.mean()),
            "median_cohort_mean_delta": float(g.mean_delta_theta.median()),
            "mean_cohort_mean_delta": float(g.mean_delta_theta.mean()),
            "median_cohort_median_delta": float(g.median_delta_theta.median()),
            "median_positive_units": float(g.positive_units.median()),
            "median_aggregate_p": float(g.aggregate_p.median()),
            "median_annual_count_ac1": float(ug.annual_count_ac1.median()),
            "median_monthly_indicator_ac1": float(ug.monthly_indicator_ac1.median()),
            "median_monthly_indicator_ac12": float(ug.monthly_indicator_ac12.median()),
            "distance_to_amazon_delta": abs(float(g.mean_delta_theta.median()) - AMAZON_DELTA_REFERENCE),
        })
    summary = pd.DataFrame(rows)

    # Add common null point to both curves in a presentation table.
    curve_rows = []
    nullrow = summary[summary.mechanism == "baseline"].iloc[0].to_dict()
    for mech in ("persistent_regime_alignment", "history_feedback"):
        nr = dict(nullrow)
        nr["curve_mechanism"] = mech
        nr["mechanism_level"] = 0.0
        curve_rows.append(nr)
        for _, r in summary[summary.mechanism == mech].sort_values("kappa").iterrows():
            rr = r.to_dict()
            rr["curve_mechanism"] = mech
            rr["mechanism_level"] = rr["kappa"]
            curve_rows.append(rr)
    curves = pd.DataFrame(curve_rows)

    # Closest simulated delta to frozen Amazon reference, separately by mechanism.
    closest = []
    for mech in ("persistent_regime_alignment", "history_feedback"):
        g = curves[curves.curve_mechanism == mech].copy()
        j = (g.median_cohort_mean_delta - AMAZON_DELTA_REFERENCE).abs().idxmin()
        r = g.loc[j]
        closest.append({
            "mechanism": mech,
            "amazon_delta_reference": AMAZON_DELTA_REFERENCE,
            "closest_kappa": float(r.mechanism_level),
            "closest_median_cohort_mean_delta": float(r.median_cohort_mean_delta),
            "absolute_delta_distance": float(abs(r.median_cohort_mean_delta-AMAZON_DELTA_REFERENCE)),
            "gate_power_at_closest_grid_point": float(r.full_gate_power),
            "power_ci_low": float(r.full_gate_power_ci_low),
            "power_ci_high": float(r.full_gate_power_ci_high),
            "note": "descriptive nearest-grid comparison only; not a universal power-at-Delta mapping",
        })
    closest = pd.DataFrame(closest)

    # Detection thresholds by simulated mechanism level.
    thresholds = []
    for mech in ("persistent_regime_alignment", "history_feedback"):
        g = curves[curves.curve_mechanism == mech].sort_values("mechanism_level")
        for target in (0.50, 0.80):
            hit = g[g.full_gate_power >= target]
            if len(hit):
                r = hit.iloc[0]
                thresholds.append({
                    "mechanism": mech,
                    "target_power": target,
                    "first_grid_kappa_reaching_target": float(r.mechanism_level),
                    "power": float(r.full_gate_power),
                    "median_cohort_mean_delta": float(r.median_cohort_mean_delta),
                })
            else:
                thresholds.append({
                    "mechanism": mech,
                    "target_power": target,
                    "first_grid_kappa_reaching_target": np.nan,
                    "power": np.nan,
                    "median_cohort_mean_delta": np.nan,
                })
    thresholds = pd.DataFrame(thresholds)

    nullc = cohorts[cohorts.mechanism == "baseline"]
    nullu = units[units.mechanism == "baseline"]
    attainable = math.floor(ALPHA*(B+1)-1e-12)/(B+1)
    tech = {
        "B": B,
        "attainable_single_test_size_strict_p_lt_005": attainable,
        "n_null_cohorts": int(len(nullc)),
        "null_full_gate_false_positive_rate": float(nullc.full_gate_pass.mean()),
        "null_aggregate_rejection_rate": float(nullc.aggregate_reject_005.mean()),
        "null_single_unit_rejection_rate": float(nullu.native_reject_005.mean()),
        "null_full_gate_calibration_gate_le_0_10": bool(nullc.full_gate_pass.mean() <= 0.10),
        "null_single_unit_calibration_gate_le_0_065": bool(nullu.native_reject_005.mean() <= 0.065),
        "native_valid_fraction_min": float(units.native_null_valid_fraction.min()),
        "native_valid_gate": bool(units.native_null_valid_fraction.min() >= 0.99),
        "native_spectral_error_median": float(units.native_spectral_error_median.median()),
        "native_spectral_error_p99": float(units.native_spectral_error_median.quantile(.99)),
        "spectral_fidelity_gate": bool(
            units.native_spectral_error_median.median() <= .01
            and units.native_spectral_error_median.quantile(.99) <= .03
        ),
        "native_count_exact_min": float(units.native_count_exact_fraction.min()),
        "native_group_count_exact_min": float(units.native_group_count_exact_fraction.min()),
        "surrogate_exactness_gate": bool(
            units.native_count_exact_fraction.min() == 1.0
            and units.native_group_count_exact_fraction.min() == 1.0
        ),
        "injection_multiset_exact_all": bool(units.injection_monthly_multiset_exact.all()),
        "injection_event_count_exact_all": bool(units.injection_monthly_event_count_exact.all()),
    }
    tech["all_technical_gates_pass"] = bool(
        tech["null_full_gate_calibration_gate_le_0_10"]
        and tech["null_single_unit_calibration_gate_le_0_065"]
        and tech["native_valid_gate"]
        and tech["spectral_fidelity_gate"]
        and tech["surrogate_exactness_gate"]
        and tech["injection_multiset_exact_all"]
        and tech["injection_event_count_exact_all"]
    )
    return summary, curves, closest, thresholds, tech


def main():
    a = parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    chunks = out / "chunks"
    chunks.mkdir(exist_ok=True)

    if a.smoke:
        designs = (("baseline", 0.0), ("persistent_regime_alignment", 0.8), ("history_feedback", 1.0))
        ncoh = SMOKE_COHORTS
        B = SMOKE_B
        chunk_size = SMOKE_COHORTS
        purpose = "technical smoke only; no scientific interpretation"
    else:
        designs = DESIGNS
        ncoh = CONFIRM_COHORTS
        B = CONFIRM_B
        chunk_size = CHUNK_SIZE
        purpose = "frozen Stage-4 power study under prespecified injected alternatives"

    probs, intercept, neff_p, neff_var = seasonal_probabilities(LAMBDA)
    meta = {
        "version": VERSION,
        "seed_namespace": SEED_NAMESPACE,
        "purpose": purpose,
        "smoke": a.smoke,
        "cohort_size": COHORT_SIZE,
        "generator_cells": [[f,p] for f,p in GENERATOR_CELLS],
        "cohort_sampling": "8 of 9 family x phi cells without replacement per cohort; composition shared across mechanisms/levels",
        "lambda": LAMBDA,
        "seasonal_probabilities": probs.tolist(),
        "neff_p": neff_p,
        "neff_var": neff_var,
        "designs": [[m,k] for m,k in designs],
        "cohorts_per_design": ncoh,
        "B": B,
        "block_years": BLOCK_YEARS,
        "q_extreme": Q_EXTREME,
        "alpha": ALPHA,
        "sim_years": SIM_YEARS,
        "output_years": OUTPUT_YEARS,
        "primary_gate": "positive delta in >=6/8 units AND synchronized-surrogate aggregate p<0.05",
        "aggregate_statistic": "mean_i[median(theta_null_i)-theta_obs_i]",
        "aggregate_null": "mean_i[median(theta_null_i)-theta_null_i,b], synchronized by surrogate index b",
        "amazon_delta_reference": AMAZON_DELTA_REFERENCE,
        "amazon_reference_usage": "descriptive nearest-grid comparison only; never used to tune mechanisms or kappa grid",
        "regime_rho": REGIME_RHO,
        "history_horizon_months": HISTORY_HORIZON,
        "history_tau_months": HISTORY_TAU,
        "history_iterations": HISTORY_ITERATIONS,
        "monthly_iaaft_iter": MONTHLY_IAAFT_ITER,
    }
    (out / "RUN_METADATA.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    tasks = []
    for mech, kap in designs:
        for start in range(0, ncoh, chunk_size):
            stop = min(ncoh, start+chunk_size)
            fn = f"{mech}_k{kap:g}_coh{start:04d}_{stop-1:04d}"
            cpath = chunks / f"{fn}_cohorts.csv"
            upath = chunks / f"{fn}_units.csv"
            if a.resume and cpath.exists() and upath.exists():
                continue
            tasks.append((mech, float(kap), start, stop, B, str(cpath), str(upath)))

    print(f"Stage 4 {VERSION}: {len(tasks)} chunks; jobs={a.jobs}; B={B}; cohorts/design={ncoh}", flush=True)

    def save(full, result):
        mech, kap, start, stop, Bx, cpath, upath = full
        cohort_rows, unit_rows = result
        for path, rows in ((Path(cpath), cohort_rows), (Path(upath), unit_rows)):
            tmp = path.with_suffix(".tmp")
            pd.DataFrame(rows).to_csv(tmp, index=False)
            tmp.replace(path)

    comp = [((t[0],t[1],t[2],t[3],t[4]), t) for t in tasks]
    if a.jobs == 1:
        for i,(ct,full) in enumerate(comp,1):
            save(full, run_chunk(ct))
            if i % 5 == 0 or i == len(comp):
                print(f"[{i}/{len(comp)} chunks]", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=a.jobs) as ex:
            fmap = {ex.submit(run_chunk, ct): full for ct,full in comp}
            for i,fut in enumerate(as_completed(fmap),1):
                save(fmap[fut], fut.result())
                if i % 5 == 0 or i == len(fmap):
                    print(f"[{i}/{len(fmap)} chunks]", flush=True)

    cfiles = sorted(chunks.glob("*_cohorts.csv"))
    ufiles = sorted(chunks.glob("*_units.csv"))
    if not cfiles or not ufiles:
        raise RuntimeError("No Stage-4 chunk files found")
    cohorts = pd.concat((pd.read_csv(f) for f in cfiles), ignore_index=True)
    units = pd.concat((pd.read_csv(f) for f in ufiles), ignore_index=True)

    expected_c = len(designs)*ncoh
    expected_u = expected_c*COHORT_SIZE
    if len(cohorts) != expected_c:
        raise RuntimeError(f"Expected {expected_c} cohort rows, found {len(cohorts)}")
    if len(units) != expected_u:
        raise RuntimeError(f"Expected {expected_u} unit rows, found {len(units)}")
    if cohorts.duplicated(["mechanism","kappa","cohort_replicate"]).any():
        raise RuntimeError("Duplicate cohort keys")
    if units.duplicated(["mechanism","kappa","cohort_replicate","unit_index"]).any():
        raise RuntimeError("Duplicate unit keys")

    cohorts = cohorts.sort_values(["mechanism","kappa","cohort_replicate"]).reset_index(drop=True)
    units = units.sort_values(["mechanism","kappa","cohort_replicate","unit_index"]).reset_index(drop=True)
    cohorts.to_csv(out / "cohort_results.csv", index=False)
    units.to_csv(out / "unit_results.csv", index=False)

    summary, curves, closest, thresholds, tech = summarize(cohorts, units, B)
    summary.to_csv(out / "power_summary.csv", index=False)
    curves.to_csv(out / "power_curves_by_mechanism.csv", index=False)
    closest.to_csv(out / "amazon_delta_nearest_grid.csv", index=False)
    thresholds.to_csv(out / "power_detection_thresholds.csv", index=False)
    (out / "TECHNICAL_GATES.json").write_text(json.dumps(tech, indent=2), encoding="utf-8")

    if a.smoke:
        decision = {
            "technical_gates_pass": None,
            "stage4_complete": False,
            "smoke_only": True,
            "note": "Smoke output is execution-only. Power estimates and confirmatory technical gates are not interpreted with B=20 and two cohorts/design.",
        }
    else:
        decision = {
            "technical_gates_pass": bool(tech["all_technical_gates_pass"]),
            "stage4_complete": bool(tech["all_technical_gates_pass"]),
            "interpretation_rule": "Power is mechanism-specific. Amazon Delta=0.085 is compared only to nearest simulated grid points; no universal interpolation is authorized.",
            "next_step_if_pass": "audit Stage 4, then integrate Stage 3 + Stage 4 + frozen Amazon result into Paper A evidence matrix",
        }
    (out / "STAGE4_DECISION.json").write_text(json.dumps(decision, indent=2), encoding="utf-8")

    print("\n=== POWER SUMMARY ===")
    print(summary.to_string(index=False))
    print("\n=== AMAZON REFERENCE: NEAREST GRID ===")
    print(closest.to_string(index=False))
    print("\n=== TECHNICAL GATES ===")
    print(json.dumps(tech, indent=2))
    print("\n=== DECISION ===")
    print(json.dumps(decision, indent=2))
    print(f"\nOutputs: {out.resolve()}")


if __name__ == "__main__":
    main()
