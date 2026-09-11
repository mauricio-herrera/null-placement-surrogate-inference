#!/usr/bin/env python3
"""Gaussian AR(1) operator-ablation control for surrogate-null transport.

Ablates the analysis operator in stages on the same Gaussian AR(1) trajectories:
  1) linear annual aggregation (annual mean),
  2) thresholding only at monthly resolution,
  3) thresholding + annual aggregation (annual counts),
  4) cubic detrending added to annual observables.

Three null constructions are compared whenever possible:
  - oracle: new trajectories from the known Gaussian AR(1) generator;
  - native: manuscript's monthly constrained surrogate, then the operator;
  - index: IAAFT applied directly after the operator.

The main cross-stage statistic is normalized total variation (NTV), which is
order-sensitive and tie-safe. Ferro-Segers is retained for the annual
continuous/full-pipeline conditions to anchor the ablation to the manuscript.
For monthly indicators, runs per event is also reported.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.special import expit
from scipy.stats import binomtest, norm, wasserstein_distance

N_YEARS_SIM = 252
N_YEARS_ANALYSIS = 251
MONTHS = 12
N_MONTHS = N_YEARS_SIM * MONTHS
N_ANALYSIS_MONTHS = N_YEARS_ANALYSIS * MONTHS
Q = 0.90


def p_profile(lam: float, K: float = 3.0) -> np.ndarray:
    m = np.arange(1, 13)
    c = np.cos(2 * np.pi * (m - 8) / 12)
    f = lambda a: expit(a + lam * c).sum() - K
    a = brentq(f, -30.0, 30.0)
    return expit(a + lam * c)


def neff(p: np.ndarray) -> float:
    p = np.asarray(p, float)
    return float(p.sum() ** 2 / np.sum(p * p))


def generate_ar1(phi: float, n: int, rng: np.random.Generator, B: int | None = None) -> np.ndarray:
    sig = math.sqrt(1.0 - phi * phi)
    if B is None:
        x = np.empty(n, dtype=float)
        x[0] = rng.normal()
        eps = rng.normal(size=n - 1)
        for t in range(1, n):
            x[t] = phi * x[t - 1] + sig * eps[t - 1]
        return x
    x = np.empty((B, n), dtype=float)
    x[:, 0] = rng.normal(size=B)
    eps = rng.normal(size=(B, n - 1))
    for t in range(1, n):
        x[:, t] = phi * x[:, t - 1] + sig * eps[:, t - 1]
    return x


def annual_mean(x: np.ndarray) -> np.ndarray:
    if x.ndim == 1:
        return x[:N_ANALYSIS_MONTHS].reshape(N_YEARS_ANALYSIS, MONTHS).mean(axis=1)
    return x[:, :N_ANALYSIS_MONTHS].reshape(x.shape[0], N_YEARS_ANALYSIS, MONTHS).mean(axis=2)


def monthly_indicator(x: np.ndarray, thresholds: np.ndarray) -> np.ndarray:
    if x.ndim == 1:
        xx = x[:N_ANALYSIS_MONTHS].reshape(N_YEARS_ANALYSIS, MONTHS)
        return (xx < thresholds[None, :]).astype(float).reshape(-1)
    xx = x[:, :N_ANALYSIS_MONTHS].reshape(x.shape[0], N_YEARS_ANALYSIS, MONTHS)
    return (xx < thresholds[None, None, :]).astype(float).reshape(x.shape[0], -1)


def annual_count(x: np.ndarray, thresholds: np.ndarray) -> np.ndarray:
    if x.ndim == 1:
        return monthly_indicator(x, thresholds).reshape(N_YEARS_ANALYSIS, MONTHS).sum(axis=1)
    return monthly_indicator(x, thresholds).reshape(x.shape[0], N_YEARS_ANALYSIS, MONTHS).sum(axis=2)


def cubic_detrend(y: np.ndarray) -> np.ndarray:
    n = y.shape[-1]
    t = np.linspace(-1.0, 1.0, n)
    X = np.column_stack([np.ones(n), t, t * t, t * t * t])
    if y.ndim == 1:
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        return y - X @ beta
    beta, *_ = np.linalg.lstsq(X, y.T, rcond=None)
    return y - (X @ beta).T


def ferro_segers_one(z: np.ndarray, q: float = Q) -> float:
    u = np.quantile(z, q)
    idx = np.flatnonzero(z > u)
    if idx.size < 2:
        return 1.0
    T = np.diff(idx).astype(float)
    nT = T.size
    if T.max() <= 2:
        den = nT * np.sum(T * T)
        raw = 1.0 if den <= 0 else 2.0 * (T.sum() ** 2) / den
    else:
        S = T - 1.0
        den = nT * np.sum(S * (S - 1.0))
        raw = 1.0 if den <= 0 else 2.0 * (S.sum() ** 2) / den
    return float(np.clip(raw, 0.0, 1.0))


def ferro_segers(arr: np.ndarray) -> np.ndarray:
    if arr.ndim == 1:
        return np.array([ferro_segers_one(arr)])
    return np.fromiter((ferro_segers_one(row) for row in arr), dtype=float, count=arr.shape[0])


def ntv_one(y: np.ndarray) -> float:
    y = np.asarray(y, float)
    sd = float(np.std(y, ddof=0))
    if sd <= 1e-15:
        return 0.0
    return float(np.mean(np.abs(np.diff(y))) / sd)


def ntv(arr: np.ndarray) -> np.ndarray:
    if arr.ndim == 1:
        return np.array([ntv_one(arr)])
    return np.fromiter((ntv_one(row) for row in arr), dtype=float, count=arr.shape[0])


def runs_per_event_one(ind: np.ndarray) -> float:
    b = np.asarray(ind) > 0.5
    n1 = int(b.sum())
    if n1 == 0:
        return 1.0
    starts = b & np.r_[True, ~b[:-1]]
    return float(starts.sum() / n1)


def runs_per_event(arr: np.ndarray) -> np.ndarray:
    if arr.ndim == 1:
        return np.array([runs_per_event_one(arr)])
    return np.fromiter((runs_per_event_one(row) for row in arr), dtype=float, count=arr.shape[0])


def spectral_error_batch(y: np.ndarray, target_amp: np.ndarray) -> np.ndarray:
    yc = y - y.mean(axis=1, keepdims=True)
    amp = np.abs(np.fft.rfft(yc, axis=1))
    den = np.mean(target_amp ** 2) + 1e-15
    return np.mean((amp - target_amp[None, :]) ** 2, axis=1) / den


def iaaft_batch(z: np.ndarray, B: int, rng: np.random.Generator, max_iter: int = 30) -> tuple[np.ndarray, np.ndarray]:
    z = np.asarray(z, float)
    n = z.size
    z_sorted = np.sort(z)
    target_amp = np.abs(np.fft.rfft(z - z.mean()))
    order0 = np.argsort(rng.random((B, n)), axis=1)
    cur = np.empty((B, n), float)
    cur[np.arange(B)[:, None], order0] = z_sorted[None, :]
    for _ in range(max_iter):
        cc = cur - cur.mean(axis=1, keepdims=True)
        F = np.fft.rfft(cc, axis=1)
        phase = np.exp(1j * np.angle(F))
        scores = np.fft.irfft(target_amp[None, :] * phase, n=n, axis=1)
        order = np.argsort(scores, axis=1)
        nxt = np.empty_like(cur)
        nxt[np.arange(B)[:, None], order] = z_sorted[None, :]
        cur = nxt
    return cur, spectral_error_batch(cur, target_amp)


def iaaft_binary_batch(z: np.ndarray, B: int, rng: np.random.Generator, max_iter: int = 30) -> tuple[np.ndarray, np.ndarray]:
    """IAAFT rank-remap specialized to a 0/1 multiset; equivalent target contract, faster."""
    z = np.asarray(z, float)
    n = z.size
    k = int(np.sum(z > 0.5))
    target_amp = np.abs(np.fft.rfft(z - z.mean()))
    if k == 0 or k == n:
        cur = np.tile(z[None, :], (B, 1))
        return cur, spectral_error_batch(cur, target_amp)
    rnd = rng.random((B, n))
    idx = np.argpartition(rnd, n - k, axis=1)[:, n - k:]
    cur = np.zeros((B, n), dtype=float)
    cur[np.arange(B)[:, None], idx] = 1.0
    for _ in range(max_iter):
        F = np.fft.rfft(cur - cur.mean(axis=1, keepdims=True), axis=1)
        phase = np.exp(1j * np.angle(F))
        scores = np.fft.irfft(target_amp[None, :] * phase, n=n, axis=1)
        idx = np.argpartition(scores, n - k, axis=1)[:, n - k:]
        cur.fill(0.0)
        cur[np.arange(B)[:, None], idx] = 1.0
    return cur, spectral_error_batch(cur, target_amp)


def native_surrogate_batch(x: np.ndarray, B: int, rng: np.random.Generator, max_iter: int = 12) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x, float)
    n = x.size
    month_idx = [np.arange(m, n, 12) for m in range(12)]
    mu = np.array([x[idx].mean() for idx in month_idx])
    sd = np.array([x[idx].std(ddof=0) for idx in month_idx])
    sd = np.where(sd > 0, sd, 1.0)
    mu_t = np.tile(mu, N_YEARS_SIM)
    sd_t = np.tile(sd, N_YEARS_SIM)
    zobs = (x - mu_t) / sd_t
    target_amp = np.abs(np.fft.rfft(zobs))

    cur = np.empty((B, n), float)
    sorted_groups = []
    for idx in month_idx:
        vals = np.sort(x[idx])
        sorted_groups.append(vals)
        ord0 = np.argsort(rng.random((B, idx.size)), axis=1)
        cur[:, idx] = vals[ord0]

    rr = np.arange(B)[:, None]
    for _ in range(max_iter):
        zs = (cur - mu_t[None, :]) / sd_t[None, :]
        F = np.fft.rfft(zs, axis=1)
        phase = np.exp(1j * np.angle(F))
        scores = np.fft.irfft(target_amp[None, :] * phase, n=n, axis=1)
        nxt = np.empty_like(cur)
        for idx, vals in zip(month_idx, sorted_groups):
            order = np.argsort(scores[:, idx], axis=1)
            block = np.empty((B, idx.size), float)
            block[rr, order] = vals[None, :]
            nxt[:, idx] = block
        cur = nxt
    zs = (cur - mu_t[None, :]) / sd_t[None, :]
    amp = np.abs(np.fft.rfft(zs, axis=1))
    den = np.mean(target_amp ** 2) + 1e-15
    err = np.mean((amp - target_amp[None, :]) ** 2, axis=1) / den
    return cur, err


def p_lower(obs: float, null: np.ndarray) -> float:
    return float((1.0 + np.sum(null <= obs)) / (null.size + 1.0))


def evaluate_stat(obs_y: np.ndarray, native_y: np.ndarray, oracle_y: np.ndarray, index_y: np.ndarray,
                  stat_name: str, statfun) -> dict:
    obs = float(statfun(obs_y)[0])
    ns = statfun(native_y)
    os = statfun(oracle_y)
    ix = statfun(index_y)
    # All ablation statistics are oriented so lower = stronger temporal clustering/smoothness.
    pn, po, pi = p_lower(obs, ns), p_lower(obs, os), p_lower(obs, ix)
    osd = float(np.std(os, ddof=1))
    scale = osd if osd > 1e-12 else 1.0
    return {
        "statistic": stat_name,
        "obs_stat": obs,
        "p_oracle": po,
        "p_native": pn,
        "p_index": pi,
        "rej_oracle": po < 0.05,
        "rej_native": pn < 0.05,
        "rej_index": pi < 0.05,
        "null_mean_oracle": float(np.mean(os)),
        "null_sd_oracle": osd,
        "null_mean_native": float(np.mean(ns)),
        "null_sd_native": float(np.std(ns, ddof=1)),
        "null_mean_index": float(np.mean(ix)),
        "null_sd_index": float(np.std(ix, ddof=1)),
        "wasserstein_native_oracle": float(wasserstein_distance(ns, os)),
        "wasserstein_index_oracle": float(wasserstein_distance(ix, os)),
        "wasserstein_native_oracle_scaled": float(wasserstein_distance(ns, os) / scale),
        "wasserstein_index_oracle_scaled": float(wasserstein_distance(ix, os) / scale),
    }


def make_conditions(x, native_x, oracle_x, u0, u4):
    conds = {}
    # Linear aggregation control.
    mean_obs = annual_mean(x); mean_nat = annual_mean(native_x); mean_ora = annual_mean(oracle_x)
    conds["mean_raw"] = (mean_obs, mean_nat, mean_ora, "annual", "linear aggregation")
    conds["mean_cubic"] = (cubic_detrend(mean_obs), cubic_detrend(mean_nat), cubic_detrend(mean_ora),
                           "annual", "linear aggregation + cubic detrend")

    # Threshold only, no annual aggregation: retain monthly indicator sequence.
    for lam, thr in [(0, u0), (4, u4)]:
        io = monthly_indicator(x, thr)
        inn = monthly_indicator(native_x, thr)
        ior = monthly_indicator(oracle_x, thr)
        conds[f"indicator_lam{lam}"] = (io, inn, ior, "binary_monthly", f"threshold only, lambda={lam}")

        # Threshold + annual aggregation, with/without detrending.
        co = annual_count(x, thr); cn = annual_count(native_x, thr); cr = annual_count(oracle_x, thr)
        conds[f"count_lam{lam}_raw"] = (co, cn, cr, "annual_count", f"threshold + aggregate, lambda={lam}")
        conds[f"count_lam{lam}_cubic"] = (cubic_detrend(co), cubic_detrend(cn), cubic_detrend(cr),
                                           "annual", f"threshold + aggregate + cubic, lambda={lam}")
    return conds


def condition_seed(base_seed: int, phi_i: int, rep: int, cond_i: int) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence([base_seed, 777, phi_i, rep, cond_i]))


def run(phi_values, reps, B, seed, outdir):
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    p0, p4 = p_profile(0.0), p_profile(4.0)
    u0, u4 = norm.ppf(p0), norm.ppf(p4)
    records = []
    start = time.time()
    base_ss = np.random.SeedSequence(seed)
    trajectory_seeds = base_ss.spawn(len(phi_values) * reps)
    si = 0

    condition_order = [
        "mean_raw", "mean_cubic", "indicator_lam0", "indicator_lam4",
        "count_lam0_raw", "count_lam0_cubic", "count_lam4_raw", "count_lam4_cubic"
    ]

    for phi_i, phi in enumerate(phi_values):
        for rep in range(reps):
            # Same observed trajectory seed structure as prior Gaussian control.
            obs_rng = np.random.default_rng(trajectory_seeds[si]); si += 1
            x = generate_ar1(phi, N_MONTHS, obs_rng)
            # Separate deterministic RNGs for null constructions.
            nat_rng = np.random.default_rng(np.random.SeedSequence([seed, 111, phi_i, rep]))
            ora_rng = np.random.default_rng(np.random.SeedSequence([seed, 222, phi_i, rep]))
            native_x, native_spec = native_surrogate_batch(x, B, nat_rng)
            oracle_x = generate_ar1(phi, N_MONTHS, ora_rng, B=B)
            conditions = make_conditions(x, native_x, oracle_x, u0, u4)

            for cond_i, cname in enumerate(condition_order):
                obs_y, nat_y, ora_y, kind, description = conditions[cname]
                ix_rng = condition_seed(seed, phi_i, rep, cond_i)
                if kind == "binary_monthly":
                    index_y, index_spec = iaaft_binary_batch(obs_y, B, ix_rng)
                else:
                    index_y, index_spec = iaaft_batch(obs_y, B, ix_rng)

                stat_specs = [("ntv", ntv)]
                if kind == "binary_monthly":
                    stat_specs.append(("runs_per_event", runs_per_event))
                if cname in {"mean_raw", "mean_cubic", "count_lam0_cubic", "count_lam4_cubic"}:
                    stat_specs.append(("ferro_segers", ferro_segers))

                for stat_name, statfun in stat_specs:
                    rec = evaluate_stat(obs_y, nat_y, ora_y, index_y, stat_name, statfun)
                    rec.update({
                        "phi": phi, "rep": rep, "condition": cname,
                        "operator_description": description, "kind": kind,
                        "native_spectral_error_median": float(np.median(native_spec)),
                        "index_spectral_error_median": float(np.median(index_spec)),
                    })
                    records.append(rec)

            if (rep + 1) % max(1, reps // 5) == 0:
                print(f"phi={phi:.2f}: {rep+1}/{reps}; elapsed={time.time()-start:.1f}s", flush=True)

    df = pd.DataFrame.from_records(records)
    df.to_csv(outdir / "operator_ablation_trajectory_results.csv", index=False)

    def aggregate(g, cond, stat, phi_label):
        A = g.rej_index.to_numpy(bool); N = g.rej_native.to_numpy(bool); O = g.rej_oracle.to_numpy(bool)
        ni = int(np.sum(A & ~N)); nn = int(np.sum(N & ~A)); disc = ni + nn
        paired_p = 1.0 if disc == 0 else float(binomtest(min(ni, nn), disc, 0.5, alternative="two-sided").pvalue)
        return {
            "condition": cond, "operator_description": g.operator_description.iloc[0],
            "statistic": stat, "phi": phi_label, "n": len(g),
            "oracle_reject_rate": float(O.mean()), "native_reject_rate": float(N.mean()), "index_reject_rate": float(A.mean()),
            "index_only_vs_native": float(np.mean(A & ~N)), "native_only_vs_index": float(np.mean(N & ~A)),
            "index_only_count": ni, "native_only_count": nn, "paired_exact_p": paired_p,
            "index_only_vs_oracle": float(np.mean(A & ~O)), "oracle_only_vs_index": float(np.mean(O & ~A)),
            "native_only_vs_oracle": float(np.mean(N & ~O)), "oracle_only_vs_native": float(np.mean(O & ~N)),
            "median_null_sd_oracle": float(g.null_sd_oracle.median()),
            "median_null_sd_native": float(g.null_sd_native.median()),
            "median_null_sd_index": float(g.null_sd_index.median()),
            "median_sd_ratio_native_oracle": float(np.median(g.null_sd_native / g.null_sd_oracle.replace(0, np.nan))),
            "median_sd_ratio_index_oracle": float(np.median(g.null_sd_index / g.null_sd_oracle.replace(0, np.nan))),
            "median_wasserstein_native_oracle_scaled": float(g.wasserstein_native_oracle_scaled.median()),
            "median_wasserstein_index_oracle_scaled": float(g.wasserstein_index_oracle_scaled.median()),
            "median_index_spectral_error": float(g.index_spectral_error_median.median()),
            "median_native_spectral_error": float(g.native_spectral_error_median.median()),
        }

    rows = []
    for (cond, stat, phi), g in df.groupby(["condition", "statistic", "phi"], sort=False):
        rows.append(aggregate(g, cond, stat, phi))
    for (cond, stat), g in df.groupby(["condition", "statistic"], sort=False):
        rows.append(aggregate(g, cond, stat, "pooled"))
    summary = pd.DataFrame(rows)
    summary.to_csv(outdir / "operator_ablation_summary.csv", index=False)

    meta = {
        "seed": seed, "reps_per_phi": reps, "B": B, "phi_values": list(phi_values),
        "years_simulated": N_YEARS_SIM, "years_analyzed": N_YEARS_ANALYSIS,
        "analysis_months": N_ANALYSIS_MONTHS,
        "iaaft_iterations": 30, "native_iterations": 12,
        "p_lambda0": p0.tolist(), "p_lambda4": p4.tolist(),
        "Neff_lambda0": neff(p0), "Neff_lambda4": neff(p4),
        "attainable_strict_alpha": math.floor(0.05 * (B + 1) - 1e-12) / (B + 1),
        "notes": [
            "All index IAAFT runs use the stated 30-iteration maximum because the manuscript source does not provide a numeric convergence tolerance.",
            "Binary IAAFT uses an argpartition implementation of the same rank-remap contract for 0/1 data.",
            "NTV = mean absolute first difference / population SD; lower tail is treated as stronger temporal smoothness/clustering.",
            "Monthly indicator conditions isolate thresholding without annual aggregation but operate at monthly rather than annual resolution; they are diagnostic, not a like-for-like replacement for the manuscript's annual test.",
        ],
        "elapsed_seconds": time.time() - start,
    }
    (outdir / "operator_ablation_metadata.json").write_text(json.dumps(meta, indent=2))
    print(summary[summary.phi.astype(str)=="pooled"].to_string(index=False), flush=True)
    print(f"Saved to {outdir}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--B", type=int, default=149)
    ap.add_argument("--seed", type=int, default=20260910)
    ap.add_argument("--phis", type=float, nargs="+", default=[0.2, 0.5, 0.8])
    ap.add_argument("--outdir", default="operator_ablation_results")
    args = ap.parse_args()
    run(args.phis, args.reps, args.B, args.seed, args.outdir)
