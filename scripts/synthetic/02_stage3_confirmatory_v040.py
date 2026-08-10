from __future__ import annotations

"""
Paper A — Stage 3 focused confirmatory synthetic experiment v0.40

Frozen after Stage-2 mechanism audit v0.30.

Primary scientific target
-------------------------
For the annual dry-month count, test whether an index-resolution IAAFT null
rejects more often than a native-resolution null passed through the same
threshold-and-aggregation operator, in known finite-sample short-memory
monthly generators.

Primary analysis
----------------
- observable: annual count
- phase: 0 (non-overlapping calendar-year windows)
- least-conditioned native contract: full-record/month-of-year groups (252 y)
- paired exact discordance test (annual-only reject vs native-only reject)

Prespecified robustness
-----------------------
- native block = 50 y
- same 3 generator families, 3 phi values, 5 seasonal concentration levels
- maxrun as a secondary confirmatory observable
- nruns is excluded from Stage 3 confirmatory claims because Stage 2 showed
  materially poorer annual IAAFT spectral fidelity for that statistic.

Frozen confirmatory design
---------------------------
- 3 families x 3 phi x 5 lambda x 300 independent generating trajectories
- B = 249 surrogates per null
- strict Monte Carlo rule p < 0.05
- with B=249, attainable null size under exchangeability = 12/250 = 0.048

Run confirmatory:
    python paperA_stage3_confirmatory_v040.py --jobs 4 \
        --out paperA_stage3_confirmatory_v040_run1

Smoke only:
    python paperA_stage3_confirmatory_v040.py --smoke --jobs 1 \
        --out paperA_stage3_smoke_v040

The confirmatory grid cannot be changed from the command line. Only --jobs,
--out, --resume, and --smoke are exposed.
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
from scipy.stats import binomtest, fisher_exact, norm

VERSION = "v0.40"
SEED_NAMESPACE = "paperA_stage3_v040"

FAMILIES = ("gaussian", "periodic", "tinnov")
PHIS = (0.2, 0.5, 0.8)
LAMBDAS = (0.0, 0.5, 1.0, 2.0, 4.0)
BLOCKS = (50, 252)
OBSERVABLES = ("count", "maxrun")
PRIMARY_OBSERVABLE = "count"
PRIMARY_BLOCK = 252
ROBUSTNESS_BLOCK = 50
PHASE = 0

CONFIRM_REPS = 300
CONFIRM_B = 249
CHUNK_SIZE = 25
SMOKE_REPS = 2
SMOKE_B = 19

NU_T = 5
OUTPUT_YEARS = 251
SIM_YEARS = 252
MONTHS = SIM_YEARS * 12
Q_EXTREME = 0.90
K_DRY = 3.0
PERIODIC_AMP = 0.65
ALPHA = 0.05
ANNUAL_IAAFT_ITER = 30
MONTHLY_IAAFT_ITER = 12


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="paperA_stage3_confirmatory_v040_run1")
    p.add_argument("--jobs", type=int, default=max(1, min(4, os.cpu_count() or 1)))
    p.add_argument("--resume", action="store_true", help="reuse completed chunk CSVs")
    p.add_argument("--smoke", action="store_true", help="tiny technical smoke; not scientific")
    return p.parse_args()


def seed_from(*parts) -> int:
    h = hashlib.sha256("|".join(map(str, parts)).encode()).digest()
    return int.from_bytes(h[:8], "little") % (2**32)


def logistic(x):
    return 1.0 / (1.0 + np.exp(-x))


def seasonal_probabilities(lam: float, K: float = K_DRY, phase_month: int = 7):
    months = np.arange(12)
    seasonal = np.cos(2 * np.pi * (months - phase_month) / 12.0)

    def f(a):
        return logistic(a + lam * seasonal).sum() - K

    a = brentq(f, -30.0, 30.0)
    probs = logistic(a + lam * seasonal)
    neff_p = probs.sum() ** 2 / np.sum(probs**2)
    v = probs * (1 - probs)
    neff_var = v.sum() ** 2 / np.sum(v**2) if np.sum(v**2) > 0 else np.nan
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
    rng = np.random.default_rng(seed_from("tinnov_reference_stage3", phi))
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


def annual_metrics(ind: np.ndarray):
    seg = ind[: OUTPUT_YEARS * 12]
    w = seg.reshape(OUTPUT_YEARS, 12)
    count = w.sum(axis=1).astype(float)
    maxrun = np.zeros(OUTPUT_YEARS, dtype=float)
    for j, row in enumerate(w):
        best = cur = 0
        for val in row:
            if val:
                cur += 1
                best = max(best, cur)
            else:
                cur = 0
        maxrun[j] = best
    return {"count": count, "maxrun": maxrun}


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


def iaaft_1d(x: np.ndarray, rng: np.random.Generator, max_iter: int = ANNUAL_IAAFT_ITER):
    x = np.asarray(x, float)
    sorted_x = np.sort(x)
    target_amp = np.abs(np.fft.rfft(x - x.mean()))
    y = rng.permutation(x)
    last = np.inf
    err = np.inf
    for _ in range(max_iter):
        fy = np.fft.rfft(y - y.mean())
        z = np.fft.irfft(target_amp * np.exp(1j * np.angle(fy)), n=len(x)) + x.mean()
        order = np.argsort(z, kind="mergesort")
        yn = np.empty_like(y)
        yn[order] = sorted_x
        amp = np.abs(np.fft.rfft(yn - yn.mean()))
        err = np.mean((amp - target_amp) ** 2) / (np.mean(target_amp**2) + 1e-15)
        y = yn
        if abs(last - err) < 1e-9:
            break
        last = err
    return y, float(err)


def monthly_groups(n: int, block_years: int):
    year = np.arange(n) // 12
    month = np.arange(n) % 12
    if block_years >= SIM_YEARS:
        gid = month
    else:
        block = year // block_years
        gid = block * 12 + month
    return [np.flatnonzero(gid == g) for g in np.unique(gid)]


def constrained_monthly_iaaft(x: np.ndarray, rng: np.random.Generator, block_years: int,
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


def one_replicate(family: str, phi: float, lam: float, rep: int, B: int):
    rng = np.random.default_rng(seed_from(SEED_NAMESPACE, family, phi, lam, rep))
    probs, intercept, neff_p, neff_var = seasonal_probabilities(lam)
    thresholds = thresholds_for(family, phi, probs)
    cal_err = threshold_calibration_error(family, phi, thresholds, probs)

    x = simulate_family(family, phi, rng)
    ind = indicators(x, thresholds)
    obs_metrics = annual_metrics(ind)

    obs_theta = {}
    annual_nulls = {}
    annual_errs = {}
    for obs_name in OBSERVABLES:
        series = detrend_cubic(obs_metrics[obs_name])
        theta_obs = theta_intervals(series)
        vals = np.empty(B)
        errs = np.empty(B)
        for b in range(B):
            ra = np.random.default_rng(seed_from("annual_stage3", family, phi, lam, rep, obs_name, b))
            sa, ea = iaaft_1d(series, ra)
            vals[b] = theta_intervals(sa)
            errs[b] = ea
        obs_theta[obs_name] = theta_obs
        annual_nulls[obs_name] = vals
        annual_errs[obs_name] = errs

    rows = []
    for block_years in BLOCKS:
        native_vals = {obs: np.empty(B) for obs in OBSERVABLES}
        native_errs = np.empty(B)
        exact_total = np.empty(B, dtype=bool)
        exact_group = np.empty(B, dtype=bool)
        orig_groups = monthly_groups(len(x), block_years)
        orig_group_counts = np.array([ind[idx].sum() for idx in orig_groups], dtype=int)

        for b in range(B):
            rn = np.random.default_rng(seed_from("native_stage3", family, phi, lam, rep, block_years, b))
            sx, en = constrained_monthly_iaaft(x, rn, block_years)
            sind = indicators(sx, thresholds)
            smetrics = annual_metrics(sind)
            native_errs[b] = en
            exact_total[b] = int(sind.sum()) == int(ind.sum())
            new_group_counts = np.array([sind[idx].sum() for idx in orig_groups], dtype=int)
            exact_group[b] = bool(np.array_equal(orig_group_counts, new_group_counts))
            for obs_name in OBSERVABLES:
                native_vals[obs_name][b] = theta_intervals(detrend_cubic(smetrics[obs_name]))

        for obs_name in OBSERVABLES:
            theta_obs = obs_theta[obs_name]
            ann = annual_nulls[obs_name]
            nat = native_vals[obs_name]
            p_ann = mc_lower_p(ann, theta_obs)
            p_nat = mc_lower_p(nat, theta_obs)
            ann_rej = bool(np.isfinite(p_ann) and p_ann < ALPHA)
            nat_rej = bool(np.isfinite(p_nat) and p_nat < ALPHA)
            rows.append({
                "family": family,
                "phi": phi,
                "lambda": lam,
                "replicate": rep,
                "observable": obs_name,
                "block_years": block_years,
                "K": K_DRY,
                "p_min": float(probs.min()),
                "p_max": float(probs.max()),
                "neff_p": float(neff_p),
                "neff_var": float(neff_var),
                "threshold_calibration_max_error": cal_err,
                "theta_obs": float(theta_obs) if np.isfinite(theta_obs) else np.nan,
                "theta_annual_null_median": float(np.nanmedian(ann)),
                "theta_native_null_median": float(np.nanmedian(nat)),
                "null_gap": float(np.nanmedian(ann) - np.nanmedian(nat)),
                "p_annual": p_ann,
                "p_native": p_nat,
                "annual_reject_005": ann_rej,
                "native_reject_005": nat_rej,
                "annual_only_reject_005": bool(ann_rej and not nat_rej),
                "native_only_reject_005": bool(nat_rej and not ann_rej),
                "both_reject_005": bool(ann_rej and nat_rej),
                "neither_reject_005": bool((not ann_rej) and (not nat_rej)),
                "annual_null_valid_fraction": float(np.mean(np.isfinite(ann))),
                "native_null_valid_fraction": float(np.mean(np.isfinite(nat))),
                "annual_spectral_error_median": float(np.nanmedian(annual_errs[obs_name])),
                "native_spectral_error_median": float(np.nanmedian(native_errs)),
                "native_count_exact_fraction": float(np.mean(exact_total)),
                "native_group_count_exact_fraction": float(np.mean(exact_group)),
            })
    return rows


def run_chunk(task):
    family, phi, lam, start, stop, B = task
    rows = []
    for rep in range(start, stop):
        rows.extend(one_replicate(family, phi, lam, rep, B))
    return rows


def wilson_interval(k: int, n: int, alpha: float = 0.05):
    if n == 0:
        return np.nan, np.nan
    z = norm.ppf(1 - alpha / 2)
    phat = k / n
    den = 1 + z*z/n
    center = (phat + z*z/(2*n)) / den
    half = z * math.sqrt(phat*(1-phat)/n + z*z/(4*n*n)) / den
    return center-half, center+half


def paired_exact_summary(df: pd.DataFrame, observable: str, block: int):
    d = df[(df.observable == observable) & (df.block_years == block)]
    a_only = int(d.annual_only_reject_005.sum())
    n_only = int(d.native_only_reject_005.sum())
    discord = a_only + n_only
    p = binomtest(a_only, discord, p=0.5, alternative="greater").pvalue if discord else 1.0
    n = len(d)
    annual = int(d.annual_reject_005.sum())
    native = int(d.native_reject_005.sum())
    alo, ahi = wilson_interval(annual, n)
    nlo, nhi = wilson_interval(native, n)
    mlo, mhi = wilson_interval(a_only, n)
    return {
        "observable": observable,
        "block_years": block,
        "n": n,
        "annual_rejections": annual,
        "annual_rejection_rate": annual/n,
        "annual_rate_ci_low": alo,
        "annual_rate_ci_high": ahi,
        "native_rejections": native,
        "native_rejection_rate": native/n,
        "native_rate_ci_low": nlo,
        "native_rate_ci_high": nhi,
        "annual_only": a_only,
        "native_only": n_only,
        "discordant_pairs": discord,
        "paired_rate_difference": (a_only-n_only)/n,
        "paired_exact_one_sided_p": float(p),
        "mechanistic_misattribution_rate": a_only/n,
        "misattrib_ci_low": mlo,
        "misattrib_ci_high": mhi,
    }


def concentration_test(df: pd.DataFrame, observable: str, block: int):
    d = df[(df.observable == observable) & (df.block_years == block)]
    low = d[d["lambda"] == 0.0]
    high = d[d["lambda"] == 4.0]
    a = int(high.annual_only_reject_005.sum())
    b = len(high) - a
    c = int(low.annual_only_reject_005.sum())
    e = len(low) - c
    odds, p = fisher_exact([[a, b], [c, e]], alternative="greater")
    return {
        "observable": observable,
        "block_years": block,
        "lambda_high": 4.0,
        "lambda_low": 0.0,
        "n_high": len(high),
        "n_low": len(low),
        "misattrib_high": a/len(high),
        "misattrib_low": c/len(low),
        "risk_difference_high_minus_low": a/len(high)-c/len(low),
        "odds_ratio": float(odds),
        "fisher_one_sided_p": float(p),
    }


def trend_test(df: pd.DataFrame, observable: str, block: int):
    """Cochran-Armitage-style one-sided trend score using frozen lambda scores."""
    d = df[(df.observable == observable) & (df.block_years == block)]
    groups = []
    for lam in LAMBDAS:
        g = d[d["lambda"] == lam]
        groups.append((float(lam), int(g.annual_only_reject_005.sum()), len(g)))
    scores = np.array([g[0] for g in groups], float)
    x = np.array([g[1] for g in groups], float)
    n = np.array([g[2] for g in groups], float)
    p0 = x.sum()/n.sum()
    wbar = np.sum(n*scores)/np.sum(n)
    num = np.sum(scores*(x-n*p0))
    den = math.sqrt(max(1e-30, p0*(1-p0)*np.sum(n*(scores-wbar)**2)))
    z = num/den if den > 0 else 0.0
    p = float(norm.sf(z))
    return {
        "observable": observable,
        "block_years": block,
        "trend_score": "lambda",
        "z": float(z),
        "one_sided_p": p,
        "common_rate": float(p0),
    }


def summarize(df: pd.DataFrame, B: int):
    cell = (
        df.groupby(["family","phi","lambda","observable","block_years"], as_index=False)
          .agg(
              n=("replicate","size"),
              annual_rejection_rate=("annual_reject_005","mean"),
              native_rejection_rate=("native_reject_005","mean"),
              mechanistic_misattribution_rate=("annual_only_reject_005","mean"),
              native_only_rate=("native_only_reject_005","mean"),
              median_null_gap=("null_gap","median"),
              annual_valid_fraction=("annual_null_valid_fraction","mean"),
              native_valid_fraction=("native_null_valid_fraction","mean"),
              annual_spectral_error_median=("annual_spectral_error_median","median"),
              native_spectral_error_median=("native_spectral_error_median","median"),
              exact_total=("native_count_exact_fraction","mean"),
              exact_group=("native_group_count_exact_fraction","mean"),
          )
    )

    by_lambda = (
        df.groupby(["lambda","observable","block_years"], as_index=False)
          .agg(
              n=("replicate","size"),
              annual_rejection_rate=("annual_reject_005","mean"),
              native_rejection_rate=("native_reject_005","mean"),
              mechanistic_misattribution_rate=("annual_only_reject_005","mean"),
              native_only_rate=("native_only_reject_005","mean"),
              median_null_gap=("null_gap","median"),
              median_neff_p=("neff_p","median"),
          )
    )

    primary_tests = []
    concentration_tests = []
    trend_tests = []
    for obs in OBSERVABLES:
        for block in BLOCKS:
            primary_tests.append(paired_exact_summary(df, obs, block))
            concentration_tests.append(concentration_test(df, obs, block))
            trend_tests.append(trend_test(df, obs, block))

    # Technical gates, frozen before confirmatory run.
    n_reject_grid = max(0, math.ceil(ALPHA*(B+1) - 1.0 - 1e-12))
    q_mc = n_reject_grid / (B+1)
    # For B=249 and strict p<.05: 12/250 = .048.
    primary = df[(df.observable == PRIMARY_OBSERVABLE) & (df.block_years == PRIMARY_BLOCK)]
    tech = {
        "B": B,
        "attainable_mc_size_strict_p_lt_005": q_mc,
        "exact_total_min": float(df.native_count_exact_fraction.min()),
        "exact_group_min": float(df.native_group_count_exact_fraction.min()),
        "annual_valid_min": float(df.annual_null_valid_fraction.min()),
        "native_valid_min": float(df.native_null_valid_fraction.min()),
        "annual_spectral_error_median": float(df[df.observable == "count"].annual_spectral_error_median.median()),
        "annual_spectral_error_p99": float(df[df.observable == "count"].annual_spectral_error_median.quantile(.99)),
        "native_spectral_error_median": float(df[df.observable == "count"].native_spectral_error_median.median()),
        "native_spectral_error_p99": float(df[df.observable == "count"].native_spectral_error_median.quantile(.99)),
        "primary_native_rejection_rate": float(primary.native_reject_005.mean()),
        "primary_native_not_anticonservative_gate": bool(primary.native_reject_005.mean() <= 0.065),
        "count_exact_gate": bool(df.native_count_exact_fraction.min() == 1.0 and df.native_group_count_exact_fraction.min() == 1.0),
        "valid_fraction_gate": bool(primary.annual_null_valid_fraction.min() >= .99 and primary.native_null_valid_fraction.min() >= .99),
        "spectral_fidelity_gate": bool(
            df[df.observable == "count"].annual_spectral_error_median.median() <= .005
            and df[df.observable == "count"].annual_spectral_error_median.quantile(.99) <= .02
            and df[df.observable == "count"].native_spectral_error_median.median() <= .01
            and df[df.observable == "count"].native_spectral_error_median.quantile(.99) <= .03
        ),
    }
    return cell, by_lambda, pd.DataFrame(primary_tests), pd.DataFrame(concentration_tests), pd.DataFrame(trend_tests), tech


def main():
    a = parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    chunks_dir = out / "chunks"
    chunks_dir.mkdir(exist_ok=True)

    if a.smoke:
        families = ("gaussian", "periodic")
        phis = (0.5,)
        lambdas = (0.0, 4.0)
        reps = SMOKE_REPS
        B = SMOKE_B
        chunk_size = SMOKE_REPS
        purpose = "technical smoke only; no scientific interpretation"
    else:
        families = FAMILIES
        phis = PHIS
        lambdas = LAMBDAS
        reps = CONFIRM_REPS
        B = CONFIRM_B
        chunk_size = CHUNK_SIZE
        purpose = "frozen Stage-3 focused confirmatory synthetic experiment"

    meta = {
        "version": VERSION,
        "seed_namespace": SEED_NAMESPACE,
        "purpose": purpose,
        "smoke": a.smoke,
        "families": list(families),
        "phis": list(phis),
        "lambdas": list(lambdas),
        "blocks": list(BLOCKS),
        "observables": list(OBSERVABLES),
        "primary_observable": PRIMARY_OBSERVABLE,
        "primary_block": PRIMARY_BLOCK,
        "robustness_block": ROBUSTNESS_BLOCK,
        "phase": PHASE,
        "reps_per_cell": reps,
        "B": B,
        "strict_mc_rule": "p < 0.05",
        "alpha": ALPHA,
        "attainable_mc_size": max(0, math.ceil(ALPHA*(B+1) - 1.0 - 1e-12))/(B+1),
        "q_extreme": Q_EXTREME,
        "K_dry": K_DRY,
        "sim_years": SIM_YEARS,
        "output_years": OUTPUT_YEARS,
        "chunk_size": chunk_size,
        "annual_iaaft_iter": ANNUAL_IAAFT_ITER,
        "monthly_iaaft_iter": MONTHLY_IAAFT_ITER,
        "confirmatory_primary_test": "paired exact discordance: annual-only > native-only, count, block=252",
        "confirmatory_concentration_test": "one-sided Fisher exact: misattribution lambda=4 > lambda=0",
        "secondary_trend_test": "one-sided Cochran-Armitage-style score using frozen lambda values",
        "technical_native_anticonservative_tolerance": "pooled primary native rejection <= 0.065",
    }
    (out / "RUN_METADATA.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    tasks = []
    for family in families:
        for phi in phis:
            for lam in lambdas:
                for start in range(0, reps, chunk_size):
                    stop = min(reps, start + chunk_size)
                    fn = f"{family}_phi{phi:g}_lam{lam:g}_rep{start:04d}_{stop-1:04d}.csv"
                    path = chunks_dir / fn
                    if a.resume and path.exists():
                        continue
                    tasks.append((family, float(phi), float(lam), start, stop, B, str(path)))

    print(f"Stage 3 {VERSION}: {len(tasks)} chunks to run; jobs={a.jobs}; B={B}; reps/cell={reps}", flush=True)

    def save_result(task, rows):
        path = Path(task[-1])
        tmp = path.with_suffix(".tmp")
        pd.DataFrame(rows).to_csv(tmp, index=False)
        tmp.replace(path)

    # Path is carried in task only for resume/output; worker receives computational subset.
    computational_tasks = [(t[:-1], t) for t in tasks]
    if a.jobs == 1:
        for i, (ct, full) in enumerate(computational_tasks, 1):
            rows = run_chunk(ct)
            save_result(full, rows)
            if i % 5 == 0 or i == len(computational_tasks):
                print(f"[{i}/{len(computational_tasks)} chunks]", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=a.jobs) as ex:
            futmap = {ex.submit(run_chunk, ct): full for ct, full in computational_tasks}
            for i, fut in enumerate(as_completed(futmap), 1):
                full = futmap[fut]
                rows = fut.result()
                save_result(full, rows)
                if i % 5 == 0 or i == len(futmap):
                    print(f"[{i}/{len(futmap)} chunks]", flush=True)

    files = sorted(chunks_dir.glob("*.csv"))
    if not files:
        raise RuntimeError("No chunk files found")
    df = pd.concat((pd.read_csv(f) for f in files), ignore_index=True)
    df = df.sort_values(["family","phi","lambda","replicate","observable","block_years"]).reset_index(drop=True)

    expected = len(families)*len(phis)*len(lambdas)*reps*len(OBSERVABLES)*len(BLOCKS)
    if len(df) != expected:
        raise RuntimeError(f"Expected {expected} result rows, found {len(df)}. Use --resume after checking chunks.")

    # Ensure uniqueness before summarizing.
    key = ["family","phi","lambda","replicate","observable","block_years"]
    if df.duplicated(key).any():
        raise RuntimeError("Duplicate confirmatory result keys detected")

    df.to_csv(out / "replicate_results.csv", index=False)
    cell, by_lambda, paired, concentration, trend, tech = summarize(df, B)
    cell.to_csv(out / "cell_summary.csv", index=False)
    by_lambda.to_csv(out / "lambda_summary.csv", index=False)
    paired.to_csv(out / "paired_confirmatory_tests.csv", index=False)
    concentration.to_csv(out / "concentration_confirmatory_tests.csv", index=False)
    trend.to_csv(out / "concentration_trend_tests.csv", index=False)
    (out / "TECHNICAL_GATES.json").write_text(json.dumps(tech, indent=2), encoding="utf-8")

    primary_row = paired[(paired.observable == PRIMARY_OBSERVABLE) & (paired.block_years == PRIMARY_BLOCK)].iloc[0]
    conc_row = concentration[(concentration.observable == PRIMARY_OBSERVABLE) & (concentration.block_years == PRIMARY_BLOCK)].iloc[0]

    decision = {
        "technical_gates_pass": bool(
            tech["primary_native_not_anticonservative_gate"]
            and tech["count_exact_gate"]
            and tech["valid_fraction_gate"]
            and tech["spectral_fidelity_gate"]
        ),
        "H1_primary_paired_direction_pass": bool(primary_row.paired_exact_one_sided_p < ALPHA and primary_row.paired_rate_difference > 0),
        "H2_high_concentration_pass": bool(conc_row.fisher_one_sided_p < ALPHA and conc_row.risk_difference_high_minus_low > 0),
        "note": "H1 is primary. H2 is tested only as the prespecified concentration contrast after H1. maxrun and trend tests are secondary.",
    }
    decision["stage3_primary_success"] = bool(
        decision["technical_gates_pass"] and decision["H1_primary_paired_direction_pass"]
    )
    decision["stage3_full_directional_success"] = bool(
        decision["stage3_primary_success"] and decision["H2_high_concentration_pass"]
    )
    (out / "CONFIRMATORY_DECISION.json").write_text(json.dumps(decision, indent=2), encoding="utf-8")

    print("\n=== STAGE 3 SUMMARY ===")
    print(paired.to_string(index=False))
    print("\n=== CONCENTRATION CONTRAST ===")
    print(concentration.to_string(index=False))
    print("\n=== TECHNICAL GATES ===")
    print(json.dumps(tech, indent=2))
    print("\n=== FROZEN DECISION ===")
    print(json.dumps(decision, indent=2))
    print(f"\nOutputs: {out.resolve()}")


if __name__ == "__main__":
    main()
