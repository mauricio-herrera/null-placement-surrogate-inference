from __future__ import annotations

"""
Paper A — Stage 5A/5B closure experiment v0.60

5A: diagnose finite-sample offset of the Ferro-Segers extremal-index deficit
    using raw (unclipped) and clipped estimates, ceiling fractions, and a
    non-clipped auxiliary exceedance-run statistic.

5B: test percentile-threshold mitigation of mechanistic misattribution.

Frozen design
-------------
- 3 generator families x 3 phi values x 150 replicates = 1350 trajectories.
- 5 threshold contracts on each evaluation trajectory:
    absolute_lambda4
    oracle_percentile
    estimated_percentile_15y
    estimated_percentile_30y
    estimated_percentile_60y
- annual/index-resolution IAAFT versus native-resolution constrained monthly
  IAAFT with block=252 (calendar-month groups over full record).
- B=149, strict p<0.05.
- count only, phase 0.

Scientific hypotheses are frozen in PROTOCOLO_STAGE5_CLOSURE_V060.md.
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
from scipy.stats import binomtest, norm

VERSION = "v0.60"
SEED_NAMESPACE = "paperA_stage5AB_v060"

FAMILIES = ("gaussian", "periodic", "tinnov")
PHIS = (0.2, 0.5, 0.8)
SCENARIOS = (
    "absolute_lambda4",
    "oracle_percentile",
    "estimated_percentile_15y",
    "estimated_percentile_30y",
    "estimated_percentile_60y",
)
REF_YEARS = {
    "estimated_percentile_15y": 15,
    "estimated_percentile_30y": 30,
    "estimated_percentile_60y": 60,
}

CONFIRM_REPS = 150
CONFIRM_B = 149
SMOKE_REPS = 2
SMOKE_B = 11
CHUNK_SIZE = 10

NU_T = 5
SIM_YEARS = 252
OUTPUT_YEARS = 251
MONTHS = SIM_YEARS * 12
Q_EXTREME = 0.90
K_DRY = 3.0
P_ORACLE = K_DRY / 12.0
LAMBDA_ABS = 4.0
PERIODIC_AMP = 0.65
ALPHA = 0.05
ANNUAL_IAAFT_ITER = 30
MONTHLY_IAAFT_ITER = 12


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="paperA_stage5AB_mitigation_v060_run1")
    p.add_argument("--jobs", type=int, default=max(1, min(4, os.cpu_count() or 1)))
    p.add_argument("--resume", action="store_true")
    p.add_argument("--smoke", action="store_true")
    return p.parse_args()


def seed_from(*parts) -> int:
    h = hashlib.sha256("|".join(map(str, parts)).encode()).digest()
    return int.from_bytes(h[:8], "little") % (2**32)


def logistic(x):
    return 1.0 / (1.0 + np.exp(-x))


def seasonal_probabilities(lam: float, K: float = K_DRY, phase_month: int = 7):
    m = np.arange(12)
    seasonal = np.cos(2 * np.pi * (m - phase_month) / 12.0)

    def f(a):
        return logistic(a + lam * seasonal).sum() - K

    a = brentq(f, -30.0, 30.0)
    probs = logistic(a + lam * seasonal)
    return probs


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


def simulate_gaussian(phi: float, rng: np.random.Generator, n: int):
    x = np.empty(n)
    x[0] = rng.normal()
    s = math.sqrt(1.0 - phi * phi)
    for t in range(1, n):
        x[t] = phi * x[t - 1] + s * rng.normal()
    return x


def simulate_periodic(phi: float, rng: np.random.Generator, n: int):
    sig, _ = periodic_sd(phi)
    burn = 2400
    xx = np.empty(n + burn)
    xx[0] = 0.0
    for t in range(1, len(xx)):
        mm = t % 12
        xx[t] = phi * xx[t - 1] + sig[mm] * rng.normal()
    return xx[burn:]


def simulate_tinnov(phi: float, rng: np.random.Generator, n: int):
    scale = math.sqrt((NU_T - 2) / NU_T)
    burn = 2000
    xx = np.empty(n + burn)
    xx[0] = 0.0
    for t in range(1, len(xx)):
        xx[t] = phi * xx[t - 1] + rng.standard_t(NU_T) * scale
    return xx[burn:]


def simulate_family(family: str, phi: float, rng: np.random.Generator, n: int):
    if family == "gaussian":
        return simulate_gaussian(phi, rng, n)
    if family == "periodic":
        return simulate_periodic(phi, rng, n)
    if family == "tinnov":
        return simulate_tinnov(phi, rng, n)
    raise ValueError(family)


@lru_cache(maxsize=None)
def tinnov_reference_sorted(phi: float):
    rng = np.random.default_rng(seed_from("tinnov_true_reference", phi))
    x = simulate_tinnov(phi, rng, 600000)
    return np.sort(x)


def thresholds_from_probs(family: str, phi: float, probs: np.ndarray):
    probs = np.asarray(probs, float)
    if family == "gaussian":
        return norm.ppf(probs)
    if family == "periodic":
        _, sd = periodic_sd(phi)
        return sd * norm.ppf(probs)
    if family == "tinnov":
        ref = tinnov_reference_sorted(float(phi))
        return np.quantile(ref, probs)
    raise ValueError(family)


def true_probs_from_thresholds(family: str, phi: float, thresholds: np.ndarray):
    thresholds = np.asarray(thresholds, float)
    if family == "gaussian":
        return norm.cdf(thresholds)
    if family == "periodic":
        _, sd = periodic_sd(phi)
        return norm.cdf(thresholds / sd)
    if family == "tinnov":
        ref = tinnov_reference_sorted(float(phi))
        return np.searchsorted(ref, thresholds, side="right") / len(ref)
    raise ValueError(family)


def estimated_monthly_thresholds(family: str, phi: float, rep: int, ref_years: int):
    rng = np.random.default_rng(seed_from(SEED_NAMESPACE, "calibration", family, phi, rep, ref_years))
    x = simulate_family(family, phi, rng, ref_years * 12)
    month = np.arange(len(x)) % 12
    return np.array([np.quantile(x[month == m], P_ORACLE) for m in range(12)], float)


def threshold_contracts(family: str, phi: float, rep: int):
    p_abs = seasonal_probabilities(LAMBDA_ABS)
    p_flat = seasonal_probabilities(0.0)
    t_abs = thresholds_from_probs(family, phi, p_abs)
    t_oracle = thresholds_from_probs(family, phi, np.full(12, P_ORACLE))
    t_lambda0 = thresholds_from_probs(family, phi, p_flat)
    d = {
        "absolute_lambda4": t_abs,
        "oracle_percentile": t_oracle,
        "estimated_percentile_15y": estimated_monthly_thresholds(family, phi, rep, 15),
        "estimated_percentile_30y": estimated_monthly_thresholds(family, phi, rep, 30),
        "estimated_percentile_60y": estimated_monthly_thresholds(family, phi, rep, 60),
    }
    return d, t_lambda0


def indicators(x: np.ndarray, thresholds: np.ndarray):
    month = np.arange(len(x)) % 12
    return (x < thresholds[month]).astype(np.int8)


def annual_count(ind: np.ndarray):
    return ind[: OUTPUT_YEARS * 12].reshape(OUTPUT_YEARS, 12).sum(axis=1).astype(float)


def detrend_cubic(x: np.ndarray):
    x = np.asarray(x, float)
    t = np.linspace(-1.0, 1.0, len(x))
    coef = np.polyfit(t, x, 3)
    return x - np.polyval(coef, t)


def theta_intervals_raw_clipped(x: np.ndarray, q: float = Q_EXTREME):
    x = np.asarray(x, float)
    if not np.all(np.isfinite(x)):
        return np.nan, np.nan
    u = np.quantile(x, q)
    idx = np.flatnonzero(x > u)
    if len(idx) < 3:
        return np.nan, np.nan
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
        raw = 1.0
    else:
        raw = float(num / den)
    clipped = float(np.clip(raw, 0.0, 1.0))
    return raw, clipped


def mean_extreme_run_length(x: np.ndarray, q: float = Q_EXTREME):
    x = np.asarray(x, float)
    if not np.all(np.isfinite(x)):
        return np.nan
    u = np.quantile(x, q)
    b = x > u
    runs = []
    cur = 0
    for val in b:
        if val:
            cur += 1
        elif cur:
            runs.append(cur)
            cur = 0
    if cur:
        runs.append(cur)
    return float(np.mean(runs)) if runs else np.nan


def iaaft_annual(x: np.ndarray, rng: np.random.Generator, max_iter: int = ANNUAL_IAAFT_ITER):
    x = np.asarray(x, float)
    sorted_x = np.sort(x)
    target = np.abs(np.fft.rfft(x - x.mean()))
    y = rng.permutation(x)
    last = np.inf
    err = np.inf
    for _ in range(max_iter):
        fy = np.fft.rfft(y - y.mean())
        phase = np.exp(1j * np.angle(fy))
        z = np.fft.irfft(target * phase, n=len(x)) + x.mean()
        order = np.argsort(z, kind="mergesort")
        yn = np.empty_like(y)
        yn[order] = sorted_x
        amp = np.abs(np.fft.rfft(yn - yn.mean()))
        err = np.mean((amp - target) ** 2) / (np.mean(target**2) + 1e-15)
        y = yn
        if abs(last - err) < 1e-9:
            break
        last = err
    return y, float(err)


def monthly_groups_full_record(n: int):
    month = np.arange(n) % 12
    return [np.flatnonzero(month == m) for m in range(12)]


def constrained_monthly_iaaft_full_record(x: np.ndarray, rng: np.random.Generator,
                                           max_iter: int = MONTHLY_IAAFT_ITER):
    x = np.asarray(x, float)
    month = np.arange(len(x)) % 12
    mu = np.array([x[month == m].mean() for m in range(12)])
    sd = np.array([x[month == m].std(ddof=0) for m in range(12)])
    sd = np.where(sd > 0, sd, 1.0)
    z = (x - mu[month]) / sd[month]
    groups = monthly_groups_full_record(len(x))
    sorted_vals = [np.sort(x[idx]) for idx in groups]
    target = np.abs(np.fft.rfft(z - z.mean()))

    y = x.copy()
    for idx in groups:
        y[idx] = rng.permutation(y[idx])

    last = np.inf
    err = np.inf
    for _ in range(max_iter):
        yz = (y - mu[month]) / sd[month]
        fy = np.fft.rfft(yz - yz.mean())
        zz = np.fft.irfft(target * np.exp(1j * np.angle(fy)), n=len(x)) + z.mean()
        yn = np.empty_like(y)
        for idx, vals in zip(groups, sorted_vals):
            order = np.argsort(zz[idx], kind="mergesort")
            tmp = np.empty(len(idx))
            tmp[order] = vals
            yn[idx] = tmp
        ynz = (yn - mu[month]) / sd[month]
        amp = np.abs(np.fft.rfft(ynz - ynz.mean()))
        err = np.mean((amp - target) ** 2) / (np.mean(target**2) + 1e-15)
        y = yn
        if abs(last - err) < 1e-8:
            break
        last = err
    return y, float(err)


def mc_lower_p(null_theta: np.ndarray, obs: float):
    valid = np.isfinite(null_theta)
    if not np.isfinite(obs) or valid.sum() == 0:
        return np.nan
    return float((1 + np.sum(null_theta[valid] <= obs)) / (valid.sum() + 1))


def one_replicate(family: str, phi: float, rep: int, B: int):
    rng = np.random.default_rng(seed_from(SEED_NAMESPACE, "evaluation", family, phi, rep))
    x = simulate_family(family, phi, rng, MONTHS)
    contracts, t_lambda0 = threshold_contracts(family, phi, rep)

    # Gate P0: oracle percentile and lambda=0 must produce identical indicators.
    ind_oracle = indicators(x, contracts["oracle_percentile"])
    ind_lambda0 = indicators(x, t_lambda0)
    p0_identical = bool(np.array_equal(ind_oracle, ind_lambda0))

    obs = {}
    for scenario, th in contracts.items():
        ind = indicators(x, th)
        annual = detrend_cubic(annual_count(ind))
        raw, clip = theta_intervals_raw_clipped(annual)
        run = mean_extreme_run_length(annual)
        true_p = true_probs_from_thresholds(family, phi, th)
        neff_p = (true_p.sum() ** 2 / np.sum(true_p**2)) if np.sum(true_p**2) > 0 else np.nan
        obs[scenario] = {
            "thresholds": th,
            "indicator": ind,
            "annual": annual,
            "theta_raw": raw,
            "theta_clip": clip,
            "runlen": run,
            "true_p_rmse": float(np.sqrt(np.mean((true_p - P_ORACLE) ** 2))),
            "true_p_range": float(true_p.max() - true_p.min()),
            "true_p_neff": float(neff_p),
            "total_events": int(ind.sum()),
            "month_counts": np.array([ind[np.arange(len(ind)) % 12 == m].sum() for m in range(12)], int),
        }

    # Annual nulls are scenario-specific.
    ann = {}
    for scenario in SCENARIOS:
        series = obs[scenario]["annual"]
        rawv = np.empty(B)
        clipv = np.empty(B)
        runv = np.empty(B)
        errv = np.empty(B)
        for b in range(B):
            ra = np.random.default_rng(seed_from(SEED_NAMESPACE, "annual", family, phi, rep, scenario, b))
            sa, ea = iaaft_annual(series, ra)
            rr, cc = theta_intervals_raw_clipped(sa)
            rawv[b], clipv[b], runv[b], errv[b] = rr, cc, mean_extreme_run_length(sa), ea
        ann[scenario] = (rawv, clipv, runv, errv)

    # Native monthly surrogates are generated once per trajectory and reused across thresholds.
    nat_raw = {s: np.empty(B) for s in SCENARIOS}
    nat_clip = {s: np.empty(B) for s in SCENARIOS}
    nat_run = {s: np.empty(B) for s in SCENARIOS}
    nat_exact_total = {s: np.empty(B, dtype=bool) for s in SCENARIOS}
    nat_exact_month = {s: np.empty(B, dtype=bool) for s in SCENARIOS}
    nat_err = np.empty(B)
    month_id = np.arange(len(x)) % 12

    for b in range(B):
        rn = np.random.default_rng(seed_from(SEED_NAMESPACE, "native", family, phi, rep, b))
        sx, en = constrained_monthly_iaaft_full_record(x, rn)
        nat_err[b] = en
        for scenario in SCENARIOS:
            th = obs[scenario]["thresholds"]
            sind = indicators(sx, th)
            sannual = detrend_cubic(annual_count(sind))
            rr, cc = theta_intervals_raw_clipped(sannual)
            nat_raw[scenario][b] = rr
            nat_clip[scenario][b] = cc
            nat_run[scenario][b] = mean_extreme_run_length(sannual)
            nat_exact_total[scenario][b] = int(sind.sum()) == obs[scenario]["total_events"]
            counts = np.array([sind[month_id == m].sum() for m in range(12)], int)
            nat_exact_month[scenario][b] = bool(np.array_equal(counts, obs[scenario]["month_counts"]))

    rows = []
    for scenario in SCENARIOS:
        ar, ac, aru, ae = ann[scenario]
        nr, nc, nru = nat_raw[scenario], nat_clip[scenario], nat_run[scenario]
        o = obs[scenario]
        p_ann = mc_lower_p(ac, o["theta_clip"])
        p_nat = mc_lower_p(nc, o["theta_clip"])
        ann_rej = bool(np.isfinite(p_ann) and p_ann < ALPHA)
        nat_rej = bool(np.isfinite(p_nat) and p_nat < ALPHA)
        rows.append({
            "family": family,
            "phi": phi,
            "replicate": rep,
            "scenario": scenario,
            "p0_oracle_lambda0_indicator_identical": p0_identical,
            "true_p_rmse_to_025": o["true_p_rmse"],
            "true_p_range": o["true_p_range"],
            "true_p_neff": o["true_p_neff"],
            "theta_obs_raw": o["theta_raw"],
            "theta_obs_clipped": o["theta_clip"],
            "annual_null_raw_median": float(np.nanmedian(ar)),
            "annual_null_clipped_median": float(np.nanmedian(ac)),
            "native_null_raw_median": float(np.nanmedian(nr)),
            "native_null_clipped_median": float(np.nanmedian(nc)),
            "annual_raw_offset": float(np.nanmedian(ar) - o["theta_raw"]),
            "annual_clipped_offset": float(np.nanmedian(ac) - o["theta_clip"]),
            "native_raw_offset": float(np.nanmedian(nr) - o["theta_raw"]),
            "native_clipped_offset": float(np.nanmedian(nc) - o["theta_clip"]),
            "obs_raw_at_or_above_one": bool(np.isfinite(o["theta_raw"]) and o["theta_raw"] >= 1.0),
            "annual_null_ceiling_fraction": float(np.mean(ar >= 1.0)),
            "native_null_ceiling_fraction": float(np.mean(nr >= 1.0)),
            "annual_all_ceiling_free": bool(np.isfinite(o["theta_raw"]) and o["theta_raw"] < 1.0 and np.all(ar < 1.0)),
            "native_all_ceiling_free": bool(np.isfinite(o["theta_raw"]) and o["theta_raw"] < 1.0 and np.all(nr < 1.0)),
            "runlen_obs": o["runlen"],
            "annual_null_runlen_median": float(np.nanmedian(aru)),
            "native_null_runlen_median": float(np.nanmedian(nru)),
            "p_annual": p_ann,
            "p_native": p_nat,
            "annual_reject_005": ann_rej,
            "native_reject_005": nat_rej,
            "mechanistic_misattribution": bool(ann_rej and not nat_rej),
            "native_only": bool(nat_rej and not ann_rej),
            "annual_valid_fraction": float(np.mean(np.isfinite(ac))),
            "native_valid_fraction": float(np.mean(np.isfinite(nc))),
            "annual_spectral_error_median": float(np.nanmedian(ae)),
            "native_spectral_error_median": float(np.nanmedian(nat_err)),
            "native_count_exact_fraction": float(np.mean(nat_exact_total[scenario])),
            "native_month_count_exact_fraction": float(np.mean(nat_exact_month[scenario])),
        })
    return rows


def run_chunk(task):
    family, phi, start, stop, B = task
    rows = []
    for rep in range(start, stop):
        rows.extend(one_replicate(family, phi, rep, B))
    return rows


def wilson(k: int, n: int):
    if n == 0:
        return np.nan, np.nan
    z = norm.ppf(0.975)
    p = k / n
    den = 1 + z*z/n
    cen = (p + z*z/(2*n))/den
    half = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n))/den
    return cen-half, cen+half


def paired_mitigation(df: pd.DataFrame, reference: str, comparison: str):
    a = df[df.scenario == reference].set_index(["family","phi","replicate"])
    b = df[df.scenario == comparison].set_index(["family","phi","replicate"])
    idx = a.index.intersection(b.index)
    ma = a.loc[idx, "mechanistic_misattribution"].astype(bool).to_numpy()
    mb = b.loc[idx, "mechanistic_misattribution"].astype(bool).to_numpy()
    ref_only = int(np.sum(ma & ~mb))
    cmp_only = int(np.sum(~ma & mb))
    disc = ref_only + cmp_only
    p = binomtest(ref_only, disc, p=0.5, alternative="greater").pvalue if disc else 1.0
    return {
        "reference": reference,
        "comparison": comparison,
        "n": len(idx),
        "reference_rate": float(ma.mean()),
        "comparison_rate": float(mb.mean()),
        "rate_reduction": float(ma.mean() - mb.mean()),
        "reference_only": ref_only,
        "comparison_only": cmp_only,
        "discordant": disc,
        "paired_exact_one_sided_p": float(p),
        "pass": bool((ma.mean() > mb.mean()) and p < 0.05),
    }


def summarize(df: pd.DataFrame, B: int):
    rows = []
    for s in SCENARIOS:
        d = df[df.scenario == s]
        m = int(d.mechanistic_misattribution.sum())
        lo, hi = wilson(m, len(d))
        rows.append({
            "scenario": s,
            "n": len(d),
            "annual_rejection_rate": float(d.annual_reject_005.mean()),
            "native_rejection_rate": float(d.native_reject_005.mean()),
            "mechanistic_misattribution_rate": float(d.mechanistic_misattribution.mean()),
            "misattrib_ci_low": lo,
            "misattrib_ci_high": hi,
            "median_true_p_rmse_to_025": float(d.true_p_rmse_to_025.median()),
            "median_true_p_range": float(d.true_p_range.median()),
            "median_true_p_neff": float(d.true_p_neff.median()),
            "median_annual_raw_offset": float(d.annual_raw_offset.median()),
            "median_annual_clipped_offset": float(d.annual_clipped_offset.median()),
            "median_native_raw_offset": float(d.native_raw_offset.median()),
            "median_native_clipped_offset": float(d.native_clipped_offset.median()),
            "median_annual_null_ceiling_fraction": float(d.annual_null_ceiling_fraction.median()),
            "median_native_null_ceiling_fraction": float(d.native_null_ceiling_fraction.median()),
            "annual_ceiling_free_trajectory_fraction": float(d.annual_all_ceiling_free.mean()),
            "native_ceiling_free_trajectory_fraction": float(d.native_all_ceiling_free.mean()),
            "median_runlen_obs": float(d.runlen_obs.median()),
            "median_annual_null_runlen": float(d.annual_null_runlen_median.median()),
            "median_native_null_runlen": float(d.native_null_runlen_median.median()),
        })
    scenario_summary = pd.DataFrame(rows)

    p1 = paired_mitigation(df, "absolute_lambda4", "oracle_percentile")
    p2 = paired_mitigation(df, "absolute_lambda4", "estimated_percentile_30y")
    paired = pd.DataFrame([p1, p2])

    attainable = math.ceil(ALPHA * (B + 1) - 1.0 - 1e-12) / (B + 1)
    oracle = df[df.scenario == "oracle_percentile"]
    tech = {
        "version": VERSION,
        "B": B,
        "attainable_mc_size_strict_p_lt_005": attainable,
        "P0_indicator_identity_fraction": float(df.p0_oracle_lambda0_indicator_identical.mean()),
        "P0_indicator_identity_gate": bool(df.p0_oracle_lambda0_indicator_identical.all()),
        "native_count_exact_min": float(df.native_count_exact_fraction.min()),
        "native_month_count_exact_min": float(df.native_month_count_exact_fraction.min()),
        "annual_valid_min": float(df.annual_valid_fraction.min()),
        "native_valid_min": float(df.native_valid_fraction.min()),
        "annual_spectral_error_median": float(df.annual_spectral_error_median.median()),
        "annual_spectral_error_p99": float(df.annual_spectral_error_median.quantile(.99)),
        "native_spectral_error_median": float(df.native_spectral_error_median.median()),
        "native_spectral_error_p99": float(df.native_spectral_error_median.quantile(.99)),
        "oracle_native_rejection_rate": float(oracle.native_reject_005.mean()),
    }
    tech["technical_pass"] = bool(
        tech["P0_indicator_identity_gate"]
        and tech["native_count_exact_min"] == 1.0
        and tech["native_month_count_exact_min"] == 1.0
        and tech["annual_valid_min"] >= .99
        and tech["native_valid_min"] >= .99
        and tech["annual_spectral_error_median"] <= .006
        and tech["annual_spectral_error_p99"] <= .025
        and tech["native_spectral_error_median"] <= .012
        and tech["native_spectral_error_p99"] <= .035
        and tech["oracle_native_rejection_rate"] <= .07
    )
    decision = {
        "technical_pass": tech["technical_pass"],
        "H_P0_consistency_pass": tech["P0_indicator_identity_gate"],
        "H_P1_oracle_mitigation_pass": bool(p1["pass"]),
        "H_P2_30y_mitigation_pass": bool(p2["pass"]),
        "stage5AB_primary_success": bool(tech["technical_pass"] and p1["pass"]),
        "stage5AB_applied_success": bool(tech["technical_pass"] and p1["pass"] and p2["pass"]),
    }
    return scenario_summary, paired, tech, decision


def main():
    a = parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    chunks = out / "chunks"
    chunks.mkdir(exist_ok=True)

    reps = SMOKE_REPS if a.smoke else CONFIRM_REPS
    B = SMOKE_B if a.smoke else CONFIRM_B
    families = ("gaussian",) if a.smoke else FAMILIES
    phis = (0.5,) if a.smoke else PHIS

    tasks = []
    for f in families:
        for p in phis:
            for start in range(0, reps, CHUNK_SIZE):
                stop = min(reps, start + CHUNK_SIZE)
                fn = chunks / f"{f}_phi{p:g}_rep{start:04d}_{stop-1:04d}.csv"
                if a.resume and fn.exists():
                    continue
                tasks.append((f, p, start, stop, B))

    def save_chunk(task, rows):
        f, p, start, stop, _ = task
        fn = chunks / f"{f}_phi{p:g}_rep{start:04d}_{stop-1:04d}.csv"
        pd.DataFrame(rows).to_csv(fn, index=False)

    if a.jobs == 1:
        for i, task in enumerate(tasks, 1):
            save_chunk(task, run_chunk(task))
            print(f"[{i}/{len(tasks)}] {task[:4]}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=a.jobs) as ex:
            fut = {ex.submit(run_chunk, t): t for t in tasks}
            for i, ftr in enumerate(as_completed(fut), 1):
                t = fut[ftr]
                save_chunk(t, ftr.result())
                print(f"[{i}/{len(tasks)}] {t[:4]}", flush=True)

    files = sorted(chunks.glob("*.csv"))
    if not files:
        raise RuntimeError("No chunk files found")
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    df.to_csv(out / "replicate_results.csv", index=False)

    summary, paired, tech, decision = summarize(df, B)
    summary.to_csv(out / "scenario_summary.csv", index=False)
    paired.to_csv(out / "paired_mitigation_tests.csv", index=False)

    offset = (
        df.groupby("scenario", as_index=False)
          .agg(
              n=("replicate", "size"),
              obs_raw_ge1=("obs_raw_at_or_above_one", "mean"),
              annual_null_ceiling_fraction=("annual_null_ceiling_fraction", "mean"),
              native_null_ceiling_fraction=("native_null_ceiling_fraction", "mean"),
              annual_raw_offset=("annual_raw_offset", "median"),
              annual_clipped_offset=("annual_clipped_offset", "median"),
              native_raw_offset=("native_raw_offset", "median"),
              native_clipped_offset=("native_clipped_offset", "median"),
              annual_ceiling_free_fraction=("annual_all_ceiling_free", "mean"),
              native_ceiling_free_fraction=("native_all_ceiling_free", "mean"),
          )
    )
    offset.to_csv(out / "estimator_offset_diagnostics.csv", index=False)

    with open(out / "TECHNICAL_GATES.json", "w") as f:
        json.dump(tech, f, indent=2)
    with open(out / "STAGE5AB_DECISION.json", "w") as f:
        json.dump(decision, f, indent=2)
    metadata = {
        "version": VERSION,
        "smoke": bool(a.smoke),
        "families": list(families),
        "phis": list(phis),
        "reps_per_cell": reps,
        "B": B,
        "scenarios": list(SCENARIOS),
        "native_block_years": 252,
        "primary_outcome": "count",
        "seed_namespace": SEED_NAMESPACE,
    }
    with open(out / "RUN_METADATA.json", "w") as f:
        json.dump(metadata, f, indent=2)

    print("\n=== STAGE 5A/5B SUMMARY ===")
    print(summary.to_string(index=False))
    print("\n=== PAIRED MITIGATION TESTS ===")
    print(paired.to_string(index=False))
    print("\nTECHNICAL:", json.dumps(tech, indent=2))
    print("\nDECISION:", json.dumps(decision, indent=2))
    if a.smoke:
        print("\nSMOKE ONLY — no scientific interpretation.")


if __name__ == "__main__":
    main()
