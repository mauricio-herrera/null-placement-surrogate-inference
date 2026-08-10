from __future__ import annotations

"""
Paper A — Stage 5C independent implementation audit v0.60

This file intentionally does NOT import Stage-3 code. It independently
reimplements:
- Stage-3 generator specification and deterministic monthly->annual mapping,
- cubic detrending,
- Ferro-Segers intervals estimator,
- annual IAAFT,
- full-record calendar-month constrained monthly IAAFT,
- Monte Carlo lower-tail p-values.

The audit has three frozen checks:
C1. theta_obs identity on a preselected 450-trajectory Stage-3 subset.
C2. exact reaggregation of Stage-3 frozen CSV headline counts.
C3. independent-surrogate rerun direction on the same 450 trajectories, B=99.

See PROTOCOLO_STAGE5_CLOSURE_V060.md and STAGE5C_SELECTION_V060.csv.
"""

import argparse
import hashlib
import io
import json
import math
import os
import zipfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import binomtest, norm

VERSION = "v0.60"
AUDIT_NAMESPACE = "paperA_stage5C_independent_v060"
STAGE3_NAMESPACE = "paperA_stage3_v040"

NU_T = 5
SIM_YEARS = 252
OUTPUT_YEARS = 251
MONTHS = SIM_YEARS * 12
Q = 0.90
K_DRY = 3.0
PERIODIC_AMP = 0.65
ALPHA = 0.05
B_AUDIT = 99
B_SMOKE = 9
ANNUAL_ITER = 40
NATIVE_ITER = 20

EXPECTED_STAGE3 = {
    "annual_rejections": 1697,
    "native_rejections": 629,
    "annual_only": 1084,
    "native_only": 16,
}


def parse_args():
    p = argparse.ArgumentParser()
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--stage3-zip")
    g.add_argument("--stage3-dir")
    p.add_argument("--selection", default="STAGE5C_SELECTION_V060.csv")
    p.add_argument("--out", default="paperA_stage5C_independent_audit_v060_run1")
    p.add_argument("--jobs", type=int, default=max(1, min(4, os.cpu_count() or 1)))
    p.add_argument("--smoke", action="store_true")
    return p.parse_args()


def sha_seed(*parts):
    digest = hashlib.sha256("|".join(map(str, parts)).encode()).digest()
    return int.from_bytes(digest[:8], "little") % (2**32)


def read_stage3_results(args):
    if args.stage3_dir:
        return pd.read_csv(Path(args.stage3_dir) / "replicate_results.csv")
    with zipfile.ZipFile(args.stage3_zip) as z:
        candidates = [n for n in z.namelist() if n.endswith("/replicate_results.csv") or n == "replicate_results.csv"]
        if len(candidates) != 1:
            raise RuntimeError(f"Expected one replicate_results.csv, found {candidates}")
        return pd.read_csv(io.BytesIO(z.read(candidates[0])))


def sigmoid(z):
    return 1.0 / (1.0 + np.exp(-z))


def event_probabilities(lam: float):
    mm = np.arange(12)
    seasonal = np.cos(2*np.pi*(mm-7)/12.0)
    root = brentq(lambda a: sigmoid(a + lam*seasonal).sum() - K_DRY, -30, 30)
    return sigmoid(root + lam*seasonal)


@lru_cache(maxsize=None)
def periodic_scales(phi: float):
    mm = np.arange(12)
    innovation_sd = np.exp(PERIODIC_AMP * np.cos(2*np.pi*mm/12.0))
    state_var = 1.0
    state_sd = np.zeros(12)
    for k in range(20000):
        m = k % 12
        state_var = phi**2 * state_var + innovation_sd[m]**2
        if k >= 19988:
            state_sd[m] = math.sqrt(state_var)
    return innovation_sd, state_sd


def gen_gaussian(phi, rng, n):
    e = rng.normal(size=n)
    x = np.empty(n)
    x[0] = e[0]
    scale = math.sqrt(1-phi*phi)
    for i in range(1, n):
        x[i] = phi*x[i-1] + scale*e[i]
    return x


def gen_periodic(phi, rng, n):
    sig, _ = periodic_scales(phi)
    burn = 2400
    # Stage-3 draws exactly one innovation for each i=1,...,n+burn-1.
    e = rng.normal(size=n+burn-1)
    x = np.zeros(n+burn)
    for i in range(1, n+burn):
        x[i] = phi*x[i-1] + sig[i % 12]*e[i-1]
    return x[burn:]


def gen_t(phi, rng, n):
    burn = 2000
    # Stage-3 draws exactly one innovation for each i=1,...,n+burn-1.
    e = rng.standard_t(NU_T, size=n+burn-1) * math.sqrt((NU_T-2)/NU_T)
    x = np.zeros(n+burn)
    for i in range(1, n+burn):
        x[i] = phi*x[i-1] + e[i-1]
    return x[burn:]


def generate(family, phi, rng, n=MONTHS):
    if family == "gaussian": return gen_gaussian(phi, rng, n)
    if family == "periodic": return gen_periodic(phi, rng, n)
    if family == "tinnov": return gen_t(phi, rng, n)
    raise ValueError(family)


@lru_cache(maxsize=None)
def stage3_t_reference(phi: float):
    rng = np.random.default_rng(sha_seed("tinnov_reference_stage3", phi))
    return np.sort(gen_t(phi, rng, 300000))


def thresholds(family, phi, lam):
    p = event_probabilities(lam)
    if family == "gaussian":
        return norm.ppf(p)
    if family == "periodic":
        _, sd = periodic_scales(phi)
        return sd * norm.ppf(p)
    ref = stage3_t_reference(float(phi))
    return np.quantile(ref, p)


def to_indicator(x, th):
    m = np.arange(x.size) % 12
    return (x < th[m]).astype(np.int8)


def annual_sum(ind):
    # independent implementation using explicit slicing, not reshape-sum helper
    y = np.empty(OUTPUT_YEARS, dtype=float)
    for k in range(OUTPUT_YEARS):
        y[k] = float(np.sum(ind[12*k:12*(k+1)]))
    return y


def cubic_residual(y):
    y = np.asarray(y, float)
    t = np.linspace(-1.0, 1.0, y.size)
    X = np.column_stack((np.ones_like(t), t, t*t, t*t*t))
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


def ferro_raw_clip(series):
    a = np.asarray(series, float)
    if not np.all(np.isfinite(a)):
        return np.nan, np.nan
    u = np.quantile(a, Q)
    loc = np.flatnonzero(a > u)
    if loc.size < 3:
        return np.nan, np.nan
    gaps = np.diff(loc).astype(float)
    n = gaps.size
    if gaps.max() <= 2:
        denominator = n * np.dot(gaps, gaps)
        numerator = 2.0 * gaps.sum()**2
    else:
        z = gaps - 1.0
        denominator = n * np.sum(z*(z-1.0))
        numerator = 2.0 * z.sum()**2
    raw = 1.0 if denominator <= 0 else numerator/denominator
    return float(raw), float(min(1.0, max(0.0, raw)))


def annual_iaaft_B(x, rng):
    """Independent IAAFT implementation with stable rank assignment."""
    x = np.asarray(x, float)
    sorted_values = np.sort(x)
    amp0 = np.abs(np.fft.rfft(x - np.mean(x)))
    y = x[rng.permutation(x.size)]
    last = None
    err = np.inf
    for _ in range(ANNUAL_ITER):
        F = np.fft.rfft(y - y.mean())
        phase = F / np.where(np.abs(F) > 0, np.abs(F), 1.0)
        candidate = np.fft.irfft(amp0 * phase, n=x.size) + x.mean()
        rank = np.argsort(candidate, kind="stable")
        next_y = np.empty_like(y)
        next_y[rank] = sorted_values
        amp = np.abs(np.fft.rfft(next_y-next_y.mean()))
        err = np.mean((amp-amp0)**2)/(np.mean(amp0**2)+1e-15)
        y = next_y
        if last is not None and abs(last-err) < 1e-10:
            break
        last = err
    return y, float(err)


def native_iaaft_B(x, rng):
    x = np.asarray(x, float)
    month = np.arange(x.size) % 12
    groups = [np.where(month == m)[0] for m in range(12)]
    means = np.array([x[g].mean() for g in groups])
    sds = np.array([x[g].std(ddof=0) for g in groups])
    sds[sds == 0] = 1.0
    standardized = (x - means[month])/sds[month]
    amp0 = np.abs(np.fft.rfft(standardized-standardized.mean()))
    sorted_by_month = [np.sort(x[g]) for g in groups]

    y = x.copy()
    for g in groups:
        y[g] = y[g][rng.permutation(g.size)]

    last = None
    err = np.inf
    for _ in range(NATIVE_ITER):
        yz = (y-means[month])/sds[month]
        F = np.fft.rfft(yz-yz.mean())
        phase = F / np.where(np.abs(F) > 0, np.abs(F), 1.0)
        target_series = np.fft.irfft(amp0*phase, n=x.size) + standardized.mean()
        next_y = np.empty_like(y)
        for g, vals in zip(groups, sorted_by_month):
            order = np.argsort(target_series[g], kind="stable")
            assigned = np.empty(g.size)
            assigned[order] = vals
            next_y[g] = assigned
        nz = (next_y-means[month])/sds[month]
        amp = np.abs(np.fft.rfft(nz-nz.mean()))
        err = np.mean((amp-amp0)**2)/(np.mean(amp0**2)+1e-15)
        y = next_y
        if last is not None and abs(last-err) < 1e-9:
            break
        last = err
    return y, float(err)


def lower_p(vals, obs):
    vals = np.asarray(vals)
    ok = np.isfinite(vals)
    return (1 + np.sum(vals[ok] <= obs))/(ok.sum()+1) if ok.any() and np.isfinite(obs) else np.nan


def audit_one(row, B):
    family = row["family"]
    phi = float(row["phi"])
    lam = float(row["lambda"])
    rep = int(row["replicate"])
    rng = np.random.default_rng(sha_seed(STAGE3_NAMESPACE, family, phi, lam, rep))
    x = generate(family, phi, rng)
    th = thresholds(family, phi, lam)
    ind = to_indicator(x, th)
    annual = cubic_residual(annual_sum(ind))
    raw_obs, clip_obs = ferro_raw_clip(annual)

    ann_theta = np.empty(B)
    ann_err = np.empty(B)
    nat_theta = np.empty(B)
    nat_err = np.empty(B)
    exact = np.empty(B, dtype=bool)
    month = np.arange(x.size)%12
    base_month_counts = np.array([ind[month==m].sum() for m in range(12)])

    for b in range(B):
        ra = np.random.default_rng(sha_seed(AUDIT_NAMESPACE, "annual", family, phi, lam, rep, b))
        sa, ea = annual_iaaft_B(annual, ra)
        ann_theta[b] = ferro_raw_clip(sa)[1]
        ann_err[b] = ea

        rn = np.random.default_rng(sha_seed(AUDIT_NAMESPACE, "native", family, phi, lam, rep, b))
        sx, en = native_iaaft_B(x, rn)
        sind = to_indicator(sx, th)
        nat_theta[b] = ferro_raw_clip(cubic_residual(annual_sum(sind)))[1]
        nat_err[b] = en
        exact[b] = np.array_equal(base_month_counts, np.array([sind[month==m].sum() for m in range(12)]))

    pa = lower_p(ann_theta, clip_obs)
    pn = lower_p(nat_theta, clip_obs)
    ar = bool(pa < ALPHA)
    nr = bool(pn < ALPHA)
    return {
        "family": family, "phi": phi, "lambda": lam, "replicate": rep,
        "theta_obs_independent": clip_obs,
        "theta_obs_raw_independent": raw_obs,
        "p_annual_independent": pa,
        "p_native_independent": pn,
        "annual_reject_independent": ar,
        "native_reject_independent": nr,
        "annual_only_independent": bool(ar and not nr),
        "native_only_independent": bool(nr and not ar),
        "annual_spectral_error_median": float(np.median(ann_err)),
        "native_spectral_error_median": float(np.median(nat_err)),
        "native_month_exact_fraction": float(exact.mean()),
    }


def main():
    a = parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    stage3 = read_stage3_results(a)
    selection = pd.read_csv(a.selection)
    if a.smoke:
        selection = selection.head(6).copy()
        B = B_SMOKE
    else:
        B = B_AUDIT

    # C2 exact recomputation from frozen CSV.
    d = stage3[(stage3.observable == "count") & (stage3.block_years == 252)]
    c2 = {
        "annual_rejections": int(d.annual_reject_005.sum()),
        "native_rejections": int(d.native_reject_005.sum()),
        "annual_only": int(d.annual_only_reject_005.sum()),
        "native_only": int(d.native_only_reject_005.sum()),
    }
    c2_pass = c2 == EXPECTED_STAGE3

    # Saved theta lookup for C1.
    saved = d[["family","phi","lambda","replicate","theta_obs"]].drop_duplicates()
    selection = selection.merge(saved, on=["family","phi","lambda","replicate"], how="left", validate="one_to_one")
    if selection.theta_obs.isna().any():
        raise RuntimeError("Selection contains trajectories absent from Stage-3 frozen CSV")

    rows = []
    records = selection.to_dict("records")
    if a.jobs == 1:
        for i, r in enumerate(records, 1):
            rows.append(audit_one(r, B))
            print(f"[{i}/{len(records)}] {r['family']} phi={r['phi']} lambda={r['lambda']} rep={r['replicate']}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=a.jobs) as ex:
            fut = {ex.submit(audit_one, r, B): r for r in records}
            for i, f in enumerate(as_completed(fut), 1):
                rows.append(f.result())
                r = fut[f]
                print(f"[{i}/{len(records)}] {r['family']} phi={r['phi']} lambda={r['lambda']} rep={r['replicate']}", flush=True)

    res = pd.DataFrame(rows)
    res = res.merge(selection[["family","phi","lambda","replicate","theta_obs"]], on=["family","phi","lambda","replicate"], how="left")
    res["theta_abs_diff"] = np.abs(res.theta_obs_independent - res.theta_obs)
    res.to_csv(out/"independent_subset_results.csv", index=False)

    max_theta_diff = float(res.theta_abs_diff.max())
    c1_pass = bool(max_theta_diff <= 1e-10)
    annual_rate = float(res.annual_reject_independent.mean())
    native_rate = float(res.native_reject_independent.mean())
    annual_only = int(res.annual_only_independent.sum())
    native_only = int(res.native_only_independent.sum())
    n = len(res)
    disc = annual_only + native_only
    paired_p = float(binomtest(annual_only, disc, 0.5, alternative="greater").pvalue) if disc else 1.0
    rate_diff = (annual_only-native_only)/n
    c3_pass = bool(native_rate <= 0.07 and annual_rate > native_rate and annual_only > native_only and rate_diff >= 0.03)

    summary = {
        "version": VERSION,
        "smoke": bool(a.smoke),
        "B": B,
        "n_selected": n,
        "C1_max_theta_abs_diff": max_theta_diff,
        "C1_deterministic_theta_identity_pass": c1_pass,
        "C2_recomputed": c2,
        "C2_expected": EXPECTED_STAGE3,
        "C2_exact_headline_reaggregation_pass": c2_pass,
        "C3_annual_rejection_rate": annual_rate,
        "C3_native_rejection_rate": native_rate,
        "C3_annual_only": annual_only,
        "C3_native_only": native_only,
        "C3_paired_rate_difference": rate_diff,
        "C3_paired_one_sided_p": paired_p,
        "C3_native_month_exact_min": float(res.native_month_exact_fraction.min()),
        "C3_annual_spectral_error_median": float(res.annual_spectral_error_median.median()),
        "C3_native_spectral_error_median": float(res.native_spectral_error_median.median()),
        "C3_independent_surrogate_direction_pass": c3_pass,
        "overall_independent_closure_pass": bool(c1_pass and c2_pass and c3_pass),
    }
    with open(out/"INDEPENDENT_AUDIT_DECISION.json","w") as f:
        json.dump(summary, f, indent=2)
    pd.DataFrame([summary]).to_csv(out/"independent_audit_summary.csv", index=False)
    print("\n=== STAGE 5C INDEPENDENT AUDIT ===")
    print(json.dumps(summary, indent=2))
    if a.smoke:
        print("\nSMOKE ONLY — C3 is not scientifically interpretable with this tiny subset/B.")


if __name__ == "__main__":
    main()
