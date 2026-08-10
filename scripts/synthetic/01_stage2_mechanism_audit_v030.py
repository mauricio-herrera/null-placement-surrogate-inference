from __future__ import annotations

"""
Paper A — Stage 2 mechanism audit v0.30
Native-resolution nulls, threshold aggregation, and mechanistic misattribution.

This stage is prespecified after the successful Stage-1 v0.21 pilot. It is NOT
an Amazon re-analysis and it is NOT the final confirmatory N=1000/B=500 run.

Frozen mechanisms examined
--------------------------
1. Native-null block conditioning: 25, 50, 100, and full-record (252 yr) blocks.
2. Aggregation phase: twelve possible starts of a 12-month non-overlapping year.
3. Sliding 12-month diagnostics (descriptive only; overlapping windows are not
   treated as an inferentially equivalent replacement for annual windows).
4. Three derived observables: annual count (primary), maximum within-window run
   length, and number of within-window runs (secondary mechanistic diagnostics).

The synthetic process is simulated for 252 years so every phase 0..11 can yield
exactly 251 non-overlapping 12-month windows without circular wrap-around.

Default mechanism-audit pilot
-----------------------------
    python paperA_stage2_mechanism_audit_v030.py --reps 12 --B 30 --jobs 4

Quick smoke
-----------
    python paperA_stage2_mechanism_audit_v030.py \
        --reps 1 --B 5 --families gaussian --phis 0.5 --lambdas 0 4 --jobs 1
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

VERSION = "v0.30"
SEED_NAMESPACE = "paperA_stage2_v030"

DEFAULT_PHIS = (0.2, 0.5, 0.8)
DEFAULT_LAMBDAS = (0.0, 0.5, 1.0, 2.0, 4.0)
DEFAULT_FAMILIES = ("gaussian", "periodic", "tinnov")
DEFAULT_BLOCKS = (25, 50, 100, 252)
DEFAULT_PHASES = tuple(range(12))
OBSERVABLES = ("count", "maxrun", "nruns")
PRIMARY_OBSERVABLE = "count"
PRIMARY_PHASE = 0
PRIMARY_BLOCK = 50

NU_T = 5
OUTPUT_YEARS = 251
SIM_YEARS = 252
MONTHS = SIM_YEARS * 12
Q_EXTREME = 0.90
K_DRY = 3.0
PERIODIC_AMP = 0.65


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="paperA_stage2_mechanism_audit_v030")
    p.add_argument("--reps", type=int, default=12)
    p.add_argument("--B", type=int, default=30)
    p.add_argument("--jobs", type=int, default=max(1, min(4, os.cpu_count() or 1)))
    p.add_argument("--families", nargs="+", default=list(DEFAULT_FAMILIES), choices=DEFAULT_FAMILIES)
    p.add_argument("--phis", nargs="+", type=float, default=list(DEFAULT_PHIS))
    p.add_argument("--lambdas", nargs="+", type=float, default=list(DEFAULT_LAMBDAS))
    p.add_argument("--blocks", nargs="+", type=int, default=list(DEFAULT_BLOCKS))
    p.add_argument("--phases", nargs="+", type=int, default=list(DEFAULT_PHASES))
    p.add_argument("--K", type=float, default=K_DRY)
    p.add_argument("--annual-iaaft-iter", type=int, default=30)
    p.add_argument("--monthly-iaaft-iter", type=int, default=12)
    return p.parse_args()


def seed_from(*parts) -> int:
    h = hashlib.sha256("|".join(map(str, parts)).encode()).digest()
    return int.from_bytes(h[:8], "little") % (2**32)


def logistic(x):
    return 1.0 / (1.0 + np.exp(-x))


def seasonal_probabilities(lam: float, K: float, phase_month: int = 7):
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
    """Periodic innovation/state SD indexed explicitly by calendar month 0..11."""
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
    # burn is divisible by 12, so calendar phase is retained.
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
    rng = np.random.default_rng(seed_from("tinnov_reference_stage2", phi))
    return simulate_tinnov(phi, rng, n=250000)


def tinnov_reference_quantiles(phi: float, probs: np.ndarray):
    return np.quantile(tinnov_reference_series(float(phi)), probs)


def thresholds_for(family: str, phi: float, probs: np.ndarray):
    if family == "gaussian":
        return norm.ppf(probs)
    if family == "periodic":
        _, sd = periodic_sd(phi)
        return sd * norm.ppf(probs)
    if family == "tinnov":
        return tinnov_reference_quantiles(phi, probs)
    raise ValueError(family)


def simulate_family(family: str, phi: float, rng: np.random.Generator, n: int = MONTHS):
    if family == "gaussian":
        return simulate_gaussian(phi, rng, n=n)
    if family == "periodic":
        x, _ = simulate_periodic(phi, rng, n=n)
        return x
    if family == "tinnov":
        return simulate_tinnov(phi, rng, n=n)
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


def indicators(x: np.ndarray, thresholds: np.ndarray):
    month = np.arange(len(x)) % 12
    return (x < thresholds[month]).astype(np.int8)


def window_metrics_from_matrix(w: np.ndarray):
    """Return count, max run of ones, and number of runs for rows of a 0/1 matrix."""
    w = np.asarray(w, dtype=np.int8)
    count = w.sum(axis=1).astype(float)
    starts = (w == 1) & np.concatenate([np.ones((len(w), 1), dtype=bool), w[:, :-1] == 0], axis=1)
    nruns = starts.sum(axis=1).astype(float)

    maxrun = np.zeros(len(w), dtype=float)
    for j, row in enumerate(w):
        best = cur = 0
        for val in row:
            if val:
                cur += 1
                if cur > best:
                    best = cur
            else:
                cur = 0
        maxrun[j] = best
    return {"count": count, "maxrun": maxrun, "nruns": nruns}


def phase_metrics(ind: np.ndarray, phase: int, years: int = OUTPUT_YEARS):
    n = years * 12
    seg = ind[phase:phase + n]
    if len(seg) != n:
        raise ValueError("Insufficient months for requested phase")
    return window_metrics_from_matrix(seg.reshape(years, 12))


def rolling_metrics(ind: np.ndarray):
    w = np.lib.stride_tricks.sliding_window_view(ind, 12)
    return window_metrics_from_matrix(w)


def detrend_cubic(x: np.ndarray):
    x = np.asarray(x, float)
    if len(x) < 5:
        return x - np.mean(x)
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


def iaaft_1d(x: np.ndarray, rng: np.random.Generator, max_iter: int = 30):
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


def constrained_monthly_iaaft(x: np.ndarray, rng: np.random.Generator, block_years: int, max_iter: int = 12):
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


def compute_one(task):
    family, phi, lam, rep, B, K, blocks, phases, annual_iter, monthly_iter = task
    rng = np.random.default_rng(seed_from(SEED_NAMESPACE, family, phi, lam, rep))
    probs, intercept, neff_p, neff_var = seasonal_probabilities(lam, K)
    thresholds = thresholds_for(family, phi, probs)
    cal_err = threshold_calibration_error(family, phi, thresholds, probs)

    x = simulate_family(family, phi, rng, n=MONTHS)
    ind = indicators(x, thresholds)

    observed = {}
    annual_null = {}
    annual_err = {}

    phase_rows = []
    for phase in phases:
        metrics = phase_metrics(ind, phase)
        for obs_name in OBSERVABLES:
            series = metrics[obs_name]
            res = detrend_cubic(series)
            theta_obs = theta_intervals(res)
            observed[(phase, obs_name)] = theta_obs

            vals = np.empty(B)
            errs = np.empty(B)
            for b in range(B):
                ra = np.random.default_rng(seed_from("annual_stage2", family, phi, lam, rep, phase, obs_name, b))
                sa, ea = iaaft_1d(res, ra, max_iter=annual_iter)
                vals[b] = theta_intervals(sa)
                errs[b] = ea
            annual_null[(phase, obs_name)] = vals
            annual_err[(phase, obs_name)] = errs
            p_ann = mc_lower_p(vals, theta_obs)
            phase_rows.append({
                "family": family, "phi": phi, "lambda": lam, "replicate": rep,
                "K": K, "p_min": float(probs.min()), "p_max": float(probs.max()),
                "neff_p": float(neff_p), "neff_var": float(neff_var),
                "threshold_calibration_max_error": cal_err,
                "phase": int(phase), "observable": obs_name,
                "theta_obs": float(theta_obs) if np.isfinite(theta_obs) else np.nan,
                "theta_annual_null_median": float(np.nanmedian(vals)),
                "delta_annual": float(np.nanmedian(vals) - theta_obs) if np.isfinite(theta_obs) else np.nan,
                "p_annual": p_ann,
                "annual_reject_005": bool(np.isfinite(p_ann) and p_ann < 0.05),
                "annual_null_valid_fraction": float(np.mean(np.isfinite(vals))),
                "annual_spectral_error_median": float(np.nanmedian(errs)),
            })

    rolling_rows = []
    rmet = rolling_metrics(ind)
    for obs_name in OBSERVABLES:
        rs = rmet[obs_name]
        theta_roll = theta_intervals(detrend_cubic(rs))
        rolling_rows.append({
            "family": family, "phi": phi, "lambda": lam, "replicate": rep,
            "observable": obs_name,
            "rolling_n_windows": len(rs),
            "theta_rolling_obs": float(theta_roll) if np.isfinite(theta_roll) else np.nan,
            "rolling_descriptive_only": True,
        })

    native_rows = []
    for block_years in blocks:
        native_vals = {(phase, obs_name): np.empty(B) for phase in phases for obs_name in OBSERVABLES}
        native_errs = np.empty(B)
        exact_total = np.empty(B, dtype=bool)
        exact_group = np.empty(B, dtype=bool)

        orig_groups = monthly_groups(len(x), block_years)
        orig_group_counts = np.array([ind[idx].sum() for idx in orig_groups], dtype=int)

        for b in range(B):
            rn = np.random.default_rng(seed_from("native_stage2", family, phi, lam, rep, block_years, b))
            sx, en = constrained_monthly_iaaft(x, rn, block_years=block_years, max_iter=monthly_iter)
            sind = indicators(sx, thresholds)
            native_errs[b] = en
            exact_total[b] = int(sind.sum()) == int(ind.sum())
            new_group_counts = np.array([sind[idx].sum() for idx in orig_groups], dtype=int)
            exact_group[b] = bool(np.array_equal(orig_group_counts, new_group_counts))

            for phase in phases:
                smet = phase_metrics(sind, phase)
                for obs_name in OBSERVABLES:
                    native_vals[(phase, obs_name)][b] = theta_intervals(detrend_cubic(smet[obs_name]))

        for phase in phases:
            for obs_name in OBSERVABLES:
                vals = native_vals[(phase, obs_name)]
                theta_obs = observed[(phase, obs_name)]
                p_nat = mc_lower_p(vals, theta_obs)
                annvals = annual_null[(phase, obs_name)]
                p_ann = mc_lower_p(annvals, theta_obs)
                ann_rej = bool(np.isfinite(p_ann) and p_ann < 0.05)
                nat_rej = bool(np.isfinite(p_nat) and p_nat < 0.05)
                native_rows.append({
                    "family": family, "phi": phi, "lambda": lam, "replicate": rep,
                    "K": K, "neff_p": float(neff_p), "neff_var": float(neff_var),
                    "phase": int(phase), "observable": obs_name,
                    "block_years": int(block_years),
                    "theta_obs": float(theta_obs) if np.isfinite(theta_obs) else np.nan,
                    "theta_annual_null_median": float(np.nanmedian(annvals)),
                    "theta_native_null_median": float(np.nanmedian(vals)),
                    "delta_annual": float(np.nanmedian(annvals) - theta_obs) if np.isfinite(theta_obs) else np.nan,
                    "delta_native": float(np.nanmedian(vals) - theta_obs) if np.isfinite(theta_obs) else np.nan,
                    "null_gap": float(np.nanmedian(annvals) - np.nanmedian(vals)),
                    "p_annual": p_ann, "p_native": p_nat,
                    "annual_reject_005": ann_rej, "native_reject_005": nat_rej,
                    "mechanistic_misattribution_005": bool(ann_rej and not nat_rej),
                    "native_null_valid_fraction": float(np.mean(np.isfinite(vals))),
                    "native_spectral_error_median": float(np.nanmedian(native_errs)),
                    "native_count_exact_fraction": float(np.mean(exact_total)),
                    "native_group_count_exact_fraction": float(np.mean(exact_group)),
                })

    return phase_rows, native_rows, rolling_rows


def summarize_outputs(phase_df: pd.DataFrame, native_df: pd.DataFrame, rolling_df: pd.DataFrame):
    phase_summary = (
        native_df[native_df["block_years"] == PRIMARY_BLOCK]
        .groupby(["family", "phi", "lambda", "observable", "phase"], as_index=False)
        .agg(
            n=("replicate", "size"),
            annual_rejection_rate=("annual_reject_005", "mean"),
            native_rejection_rate=("native_reject_005", "mean"),
            mechanistic_misattribution_rate=("mechanistic_misattribution_005", "mean"),
            median_delta_annual=("delta_annual", "median"),
            median_delta_native=("delta_native", "median"),
            median_null_gap=("null_gap", "median"),
            native_valid_fraction=("native_null_valid_fraction", "mean"),
        )
    )

    block_summary = (
        native_df[native_df["phase"] == PRIMARY_PHASE]
        .groupby(["family", "phi", "lambda", "observable", "block_years"], as_index=False)
        .agg(
            n=("replicate", "size"),
            annual_rejection_rate=("annual_reject_005", "mean"),
            native_rejection_rate=("native_reject_005", "mean"),
            mechanistic_misattribution_rate=("mechanistic_misattribution_005", "mean"),
            median_delta_native=("delta_native", "median"),
            median_null_gap=("null_gap", "median"),
            native_count_exact_fraction=("native_count_exact_fraction", "mean"),
            native_group_count_exact_fraction=("native_group_count_exact_fraction", "mean"),
            native_spectral_error=("native_spectral_error_median", "median"),
            native_valid_fraction=("native_null_valid_fraction", "mean"),
        )
    )

    overview = (
        native_df[(native_df["observable"] == PRIMARY_OBSERVABLE) & (native_df["phase"] == PRIMARY_PHASE)]
        .groupby(["lambda", "block_years"], as_index=False)
        .agg(
            n=("replicate", "size"),
            annual_rejection_rate=("annual_reject_005", "mean"),
            native_rejection_rate=("native_reject_005", "mean"),
            mechanistic_misattribution_rate=("mechanistic_misattribution_005", "mean"),
            median_null_gap=("null_gap", "median"),
        )
    )

    phase_overview = (
        native_df[(native_df["observable"] == PRIMARY_OBSERVABLE) & (native_df["block_years"] == PRIMARY_BLOCK)]
        .groupby(["lambda", "phase"], as_index=False)
        .agg(
            n=("replicate", "size"),
            annual_rejection_rate=("annual_reject_005", "mean"),
            native_rejection_rate=("native_reject_005", "mean"),
            mechanistic_misattribution_rate=("mechanistic_misattribution_005", "mean"),
            median_delta_native=("delta_native", "median"),
            median_null_gap=("null_gap", "median"),
        )
    )

    phase_dispersion = (
        native_df[(native_df["observable"] == PRIMARY_OBSERVABLE) & (native_df["block_years"] == PRIMARY_BLOCK)]
        .groupby(["family", "phi", "lambda", "replicate"], as_index=False)
        .agg(
            phase_min_delta_native=("delta_native", "min"),
            phase_max_delta_native=("delta_native", "max"),
            phase_sd_delta_native=("delta_native", "std"),
            phase_min_null_gap=("null_gap", "min"),
            phase_max_null_gap=("null_gap", "max"),
        )
    )
    phase_dispersion["phase_range_delta_native"] = (
        phase_dispersion["phase_max_delta_native"] - phase_dispersion["phase_min_delta_native"]
    )
    phase_dispersion["phase_range_null_gap"] = (
        phase_dispersion["phase_max_null_gap"] - phase_dispersion["phase_min_null_gap"]
    )

    rolling_summary = (
        rolling_df.groupby(["family", "phi", "lambda", "observable"], as_index=False)
        .agg(
            n=("replicate", "size"),
            median_theta_rolling=("theta_rolling_obs", "median"),
            valid_fraction=("theta_rolling_obs", lambda x: float(np.mean(np.isfinite(x)))),
        )
    )
    return phase_summary, block_summary, overview, phase_overview, phase_dispersion, rolling_summary


def main():
    a = parse_args()
    if any(p < 0 or p > 11 for p in a.phases):
        raise ValueError("phases must be integers 0..11")
    if 50 not in a.blocks:
        raise ValueError("Frozen mechanism audit requires block 50 in --blocks")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    tasks = []
    for family in a.families:
        for phi in a.phis:
            for lam in a.lambdas:
                for rep in range(a.reps):
                    tasks.append((family, float(phi), float(lam), rep, a.B, a.K,
                                  tuple(a.blocks), tuple(a.phases),
                                  a.annual_iaaft_iter, a.monthly_iaaft_iter))

    meta = {
        "version": VERSION,
        "seed_namespace": SEED_NAMESPACE,
        "purpose": "prespecified Stage-2 synthetic mechanism audit; not final confirmatory run",
        "families": a.families,
        "phis": a.phis,
        "lambdas": a.lambdas,
        "reps": a.reps,
        "B": a.B,
        "K": a.K,
        "output_years_per_phase": OUTPUT_YEARS,
        "sim_years": SIM_YEARS,
        "q_extreme": Q_EXTREME,
        "blocks": a.blocks,
        "phases": a.phases,
        "observables": list(OBSERVABLES),
        "primary_observable": PRIMARY_OBSERVABLE,
        "primary_phase": PRIMARY_PHASE,
        "primary_block": PRIMARY_BLOCK,
        "nu_t": NU_T,
        "periodic_amp": PERIODIC_AMP,
        "rolling_windows": "descriptive only; overlapping 12-month windows",
    }
    (out / "RUN_METADATA.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    phase_rows, native_rows, rolling_rows = [], [], []
    if a.jobs == 1:
        for i, task in enumerate(tasks, 1):
            p, n, r = compute_one(task)
            phase_rows.extend(p); native_rows.extend(n); rolling_rows.extend(r)
            if i % 5 == 0 or i == len(tasks):
                print(f"[{i}/{len(tasks)}]", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=a.jobs) as ex:
            futs = [ex.submit(compute_one, t) for t in tasks]
            for i, fut in enumerate(as_completed(futs), 1):
                p, n, r = fut.result()
                phase_rows.extend(p); native_rows.extend(n); rolling_rows.extend(r)
                if i % 5 == 0 or i == len(tasks):
                    print(f"[{i}/{len(tasks)}]", flush=True)

    phase_df = pd.DataFrame(phase_rows).sort_values(["family", "phi", "lambda", "replicate", "observable", "phase"])
    native_df = pd.DataFrame(native_rows).sort_values(["family", "phi", "lambda", "replicate", "observable", "phase", "block_years"])
    rolling_df = pd.DataFrame(rolling_rows).sort_values(["family", "phi", "lambda", "replicate", "observable"])

    phase_df.to_csv(out / "phase_observed_and_annual.csv", index=False)
    native_df.to_csv(out / "native_block_results.csv", index=False)
    rolling_df.to_csv(out / "rolling_diagnostics.csv", index=False)

    summaries = summarize_outputs(phase_df, native_df, rolling_df)
    names = [
        "phase_summary_block50.csv",
        "block_summary_phase0.csv",
        "mechanism_overview_count_phase0.csv",
        "phase_overview_count_block50.csv",
        "phase_dispersion_count_block50.csv",
        "rolling_summary.csv",
    ]
    for name, df in zip(names, summaries):
        df.to_csv(out / name, index=False)

    overview = summaries[2]
    phase_overview = summaries[3]

    print("\n=== STAGE-2 MECHANISM OVERVIEW: COUNT, PHASE 0 ===")
    print(overview.to_string(index=False))
    print("\n=== STAGE-2 PHASE OVERVIEW: COUNT, BLOCK 50 ===")
    print(phase_overview.to_string(index=False))
    print(f"\nOutputs: {out.resolve()}")
    print("\nIMPORTANT: This is a mechanism audit, not the final confirmatory run.")


if __name__ == "__main__":
    main()
