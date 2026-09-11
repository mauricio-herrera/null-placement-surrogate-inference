#!/usr/bin/env python3
"""SINCA PM2.5 null-transport analysis — prespecified V1.

Primary branch:
  daily PM2.5 -> threshold >50 -> deterministic missing-day expectation fill ->
  417 Monday-Sunday weekly burdens -> cubic + 3-harmonic residualization ->
  Ferro-Segers q=.90.

Contracts:
  (1) native daily constrained marginal-spectrum surrogate, transported through
      the identical threshold/aggregation/preprocessing pipeline;
  (2) weekly-index IAAFT on the observed residual series.

The cohort and inferential choices are frozen in the accompanying protocol.
"""
from __future__ import annotations
import argparse, hashlib, json, math, os, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

START = pd.Timestamp("2018-01-01")
END = pd.Timestamp("2025-12-28")
N_WEEKS = 417
N_DAYS = N_WEEKS * 7
B_DEFAULT = 499
BASE_SEED_DEFAULT = 20260911
EXPECTED_COHORT_SHA = "59b019cd344b1cdf9c0f6fb90b381ad12d416d72461d94d31d6ce0039ff010e6"


def stable_int(text: str) -> int:
    return int(hashlib.sha256(text.encode("utf-8")).hexdigest()[:8], 16)


def seed_rng(base: int, *parts) -> np.random.Generator:
    ints = [int(base)] + [stable_int(str(p)) for p in parts]
    return np.random.default_rng(np.random.SeedSequence(ints))


def weekly_design(n=N_WEEKS):
    week_start = pd.date_range(START, periods=n, freq="7D")
    mid = week_start + pd.Timedelta(days=3)
    d = (mid - START).days.to_numpy(float)
    t = np.linspace(-1.0, 1.0, n)
    cols = [np.ones(n), t, t*t, t*t*t]
    for k in (1,2,3):
        phase = 2*np.pi*k*d/365.2425
        cols.extend([np.sin(phase), np.cos(phase)])
    X = np.column_stack(cols)
    pinv = np.linalg.pinv(X)
    return X, pinv

XW, PINVW = weekly_design()


def residualize_weekly(y):
    a = np.asarray(y, float)
    if a.ndim == 1:
        beta = PINVW @ a
        return a - XW @ beta
    beta = a @ PINVW.T
    return a - beta @ XW.T


def fs_one(z, q=.90):
    z=np.asarray(z,float)
    u=np.quantile(z,q)
    idx=np.flatnonzero(z>u)
    if idx.size<2:
        return 1.0
    T=np.diff(idx).astype(float); nT=T.size
    if T.max()<=2:
        den=nT*np.sum(T*T)
        raw=1.0 if den<=0 else 2*(T.sum()**2)/den
    else:
        S=T-1
        den=nT*np.sum(S*(S-1))
        raw=1.0 if den<=0 else 2*(S.sum()**2)/den
    return float(np.clip(raw,0,1))


def fs(arr,q=.90):
    a=np.asarray(arr,float)
    if a.ndim==1: return np.array([fs_one(a,q)])
    return np.fromiter((fs_one(r,q) for r in a),float,count=a.shape[0])


def tdiff_one(z):
    z=np.asarray(z,float)
    v=float(np.var(z,ddof=1))
    if v<=1e-15: return 0.0
    return float(np.mean(np.diff(z)**2)/v)


def tdiff(arr):
    a=np.asarray(arr,float)
    if a.ndim==1: return np.array([tdiff_one(a)])
    return np.fromiter((tdiff_one(r) for r in a),float,count=a.shape[0])


def p_lower(obs, null):
    return float((1 + np.sum(np.asarray(null) <= obs)) / (len(null)+1))


def iaaft_batch(z, B, rng, max_iter=30):
    z=np.asarray(z,float); n=z.size
    zs=np.sort(z)
    target=np.abs(np.fft.rfft(z-z.mean()))
    order0=np.argsort(rng.random((B,n)),axis=1)
    cur=np.empty((B,n)); cur[np.arange(B)[:,None],order0]=zs[None,:]
    for _ in range(max_iter):
        F=np.fft.rfft(cur-cur.mean(axis=1,keepdims=True),axis=1)
        phase=np.exp(1j*np.angle(F))
        scores=np.fft.irfft(target[None,:]*phase,n=n,axis=1)
        order=np.argsort(scores,axis=1)
        nxt=np.empty_like(cur); nxt[np.arange(B)[:,None],order]=zs[None,:]
        cur=nxt
    amp=np.abs(np.fft.rfft(cur-cur.mean(axis=1,keepdims=True),axis=1))
    den=np.mean(target**2)+1e-15
    err=np.mean((amp-target[None,:])**2,axis=1)/den
    return cur,err


def native_surrogate_missing(x, dates, B, rng, max_iter=12):
    """Calendar-month grouped exact-value surrogate with a fixed missing mask.

    Observed values are rank-remapped within calendar month. Missing positions
    are held at the observed calendar-month mean for FFT construction and never
    enter the exact-value multiset.
    """
    x=np.asarray(x,float)
    obs=np.isfinite(x)
    month=dates.month.to_numpy()
    n=x.size
    mu=np.empty(12); sd=np.empty(12); vals=[]; idxs=[]; miss_idx=[]
    for m in range(1,13):
        idx=np.flatnonzero((month==m)&obs)
        mis=np.flatnonzero((month==m)&(~obs))
        if idx.size<20:
            raise RuntimeError(f"calendar month {m} has only {idx.size} observed daily values")
        mu[m-1]=np.mean(x[idx]); sd[m-1]=np.std(x[idx],ddof=0)
        if sd[m-1] <= 0: sd[m-1]=1.0
        vals.append(np.sort(x[idx])); idxs.append(idx); miss_idx.append(mis)
    mut=mu[month-1]; sdt=sd[month-1]
    zobs=np.zeros(n,float); zobs[obs]=(x[obs]-mut[obs])/sdt[obs]
    target=np.abs(np.fft.rfft(zobs))
    cur=np.empty((B,n),float)
    # initialize all at month mean; replace observed positions by permutations
    cur[:,:]=mut[None,:]
    for idx,v in zip(idxs,vals):
        order=np.argsort(rng.random((B,idx.size)),axis=1)
        cur[:,idx]=v[order]
    rr=np.arange(B)[:,None]
    for _ in range(max_iter):
        z=(cur-mut[None,:])/sdt[None,:]
        z[:,~obs]=0.0
        F=np.fft.rfft(z,axis=1)
        phase=np.exp(1j*np.angle(F))
        scores=np.fft.irfft(target[None,:]*phase,n=n,axis=1)
        nxt=np.broadcast_to(mut,(B,n)).copy()
        for idx,v in zip(idxs,vals):
            order=np.argsort(scores[:,idx],axis=1)
            block=np.empty((B,idx.size)); block[rr,order]=v[None,:]
            nxt[:,idx]=block
        cur=nxt
    z=(cur-mut[None,:])/sdt[None,:]; z[:,~obs]=0.0
    amp=np.abs(np.fft.rfft(z,axis=1)); den=np.mean(target**2)+1e-15
    err=np.mean((amp-target[None,:])**2,axis=1)/den
    # exact multiset gate
    for idx,v in zip(idxs,vals):
        if not np.allclose(np.sort(cur[:,idx],axis=1),v[None,:],rtol=0,atol=1e-12):
            raise RuntimeError("native exact-value conservation gate failed")
    return cur,err


def check_native_event_count_conservation(x, surr, dates, threshold_day):
    obs=np.isfinite(x); month=dates.month.to_numpy()
    for m in range(1,13):
        idx=(month==m)&obs
        c0=int(np.sum(x[idx]>threshold_day[idx]))
        cs=np.sum(surr[:,idx]>threshold_day[idx][None,:],axis=1)
        if not np.all(cs==c0):
            raise RuntimeError(f"native threshold-count conservation gate failed for month {m}")


def fixed_month_prob(x, dates, threshold_day):
    obs=np.isfinite(x); month=dates.month.to_numpy(); p=np.empty(12)
    for m in range(1,13):
        idx=(month==m)&obs
        th = threshold_day[idx]
        if idx.sum()==0: raise RuntimeError(f"no observations in month {m}")
        p[m-1]=np.mean(x[idx] > th)
    return p


def equalized_thresholds(x, dates, fixed_threshold=50.0):
    obs=np.isfinite(x); month=dates.month.to_numpy()
    p_fixed=np.array([np.mean(x[(month==m)&obs] > fixed_threshold) for m in range(1,13)])
    pbar=float(np.mean(p_fixed))
    q=1-pbar
    uth=np.empty(12)
    for m in range(1,13):
        vals=x[(month==m)&obs]
        uth[m-1]=float(np.quantile(vals,q,method="linear"))
    return uth[month-1], pbar, p_fixed


def weekly_from_daily(x, dates, threshold_day, p_missing_month):
    x=np.asarray(x,float); obs=np.isfinite(x); month=dates.month.to_numpy()
    contrib=np.empty(N_DAYS,float)
    contrib[obs]=(x[obs] > threshold_day[obs]).astype(float)
    contrib[~obs]=p_missing_month[month[~obs]-1]
    return contrib.reshape(N_WEEKS,7).sum(axis=1)


def weekly_from_surrogates(surr, obs_mask, dates, threshold_day, p_missing_month):
    month=dates.month.to_numpy(); B=surr.shape[0]
    evt=np.zeros((B,N_DAYS),float)
    evt[:,obs_mask]=(surr[:,obs_mask] > threshold_day[obs_mask][None,:]).astype(float)
    if (~obs_mask).any():
        fill=p_missing_month[month[~obs_mask]-1]
        evt[:,~obs_mask]=fill[None,:]
    return evt.reshape(B,N_WEEKS,7).sum(axis=2)


def weekly_exposure_standardized(x, dates, threshold_day):
    obs=np.isfinite(x); evt=np.zeros(N_DAYS,float)
    evt[obs]=(x[obs] > threshold_day[obs]).astype(float)
    e=evt.reshape(N_WEEKS,7).sum(axis=1)
    n=obs.reshape(N_WEEKS,7).sum(axis=1)
    if np.any(n==0):
        return None
    return 7*e/n


def summarize_stat(obs, nat, idx, fun, name):
    o=float(fun(obs)[0]); n=fun(nat); i=fun(idx)
    pn=p_lower(o,n); pi=p_lower(o,i)
    return {
        f"obs_{name}":o,
        f"p_native_{name}":pn, f"p_index_{name}":pi,
        f"rej_native_{name}":pn<.05, f"rej_index_{name}":pi<.05,
        f"median_native_{name}":float(np.median(n)), f"median_index_{name}":float(np.median(i)),
        f"mean_native_{name}":float(np.mean(n)), f"mean_index_{name}":float(np.mean(i)),
        f"sd_native_{name}":float(np.std(n,ddof=1)), f"sd_index_{name}":float(np.std(i,ddof=1)),
        f"delta_native_{name}":float(np.median(n)-o), f"delta_index_{name}":float(np.median(i)-o),
        f"center_gap_{name}":float(np.median(i)-np.median(n)),
        f"width_ratio_{name}":float(np.std(i,ddof=1)/(np.std(n,ddof=1)+1e-15)),
    }


def by_adjust(pvals):
    p=np.asarray(pvals,float); m=len(p); c=np.sum(1/np.arange(1,m+1))
    order=np.argsort(p); ps=p[order]; adj=np.empty(m,float)
    raw=ps*m*c/np.arange(1,m+1)
    raw=np.minimum.accumulate(raw[::-1])[::-1]
    adj[order]=np.minimum(raw,1.0)
    return adj


def exact_region_signflip(df, col="paired_outcome"):
    # d_r = index-only - native-only
    ds=[]
    for region,g in df.groupby("region"):
        d=int(np.sum(g[col]=="index_only")-np.sum(g[col]=="native_only"))
        ds.append((region,d))
    dvals=np.array([d for _,d in ds],int)
    obs=int(dvals.sum())
    vals=[]
    for mask in range(1<<len(dvals)):
        signs=np.array([1 if (mask>>j)&1 else -1 for j in range(len(dvals))])
        vals.append(int(np.sum(signs*dvals)))
    vals=np.array(vals)
    p=float(np.mean(vals>=obs))
    return obs,p,ds


def region_bootstrap_spearman(df, xcol, ycol, B=2000, seed=20260911):
    regs=sorted(df.region.unique()); rng=np.random.default_rng(seed); vals=[]
    point=float(spearmanr(df[xcol],df[ycol]).statistic)
    groups={r:df[df.region==r] for r in regs}
    for _ in range(B):
        pick=rng.choice(regs,size=len(regs),replace=True)
        chunks=[]
        for j,r in enumerate(pick):
            q=groups[r].copy(); q["_boot_region"]=j; chunks.append(q)
        z=pd.concat(chunks,ignore_index=True)
        rho=spearmanr(z[xcol],z[ycol]).statistic
        if np.isfinite(rho): vals.append(float(rho))
    lo,hi=np.quantile(vals,[.025,.975]) if vals else (np.nan,np.nan)
    return point,float(lo),float(hi)


def load_data(data_path, cohort_path):
    data_path = str(data_path)
    if data_path.lower().endswith(".csv"):
        df = pd.read_csv(data_path, usecols=["station_name","date","pm25"], low_memory=False)
    else:
        try:
            df=pd.read_parquet(data_path,columns=["station_name","date","pm25"])
        except Exception as e:
            raise RuntimeError("Reading the parquet requires a pandas parquet engine (normally pyarrow). "+str(e))
    cohort=pd.read_csv(cohort_path)
    df["date"]=pd.to_datetime(df["date"])
    # The complete CSV contains pollutant-specific rows plus one consolidated
    # station-day row. For the PM2.5 analysis, retain only rows carrying a
    # PM2.5 value; there is at most one such row per station-date in the
    # frozen primary window. This is an input-shape adapter only and does not
    # alter any frozen inferential choice.
    df["pm25"]=pd.to_numeric(df["pm25"],errors="coerce")
    df=df[df["pm25"].notna()].copy()
    df=df[(df.date>=START)&(df.date<=END)&df.station_name.isin(cohort.short_code)].copy()
    if df.duplicated(["station_name","date"]).any():
        raise RuntimeError("duplicate nonmissing-PM2.5 station-date rows in primary window")
    return df,cohort


def prepare_station(df, short_code):
    dates=pd.date_range(START,END,freq="D")
    s=df[df.station_name==short_code][["date","pm25"]].set_index("date").reindex(dates)
    x=pd.to_numeric(s.pm25,errors="coerce").to_numpy(float)
    return dates,x


def run_branch(short_code, station_id, x, dates, B, seed, branch="fixed"):
    obs=np.isfinite(x)
    if branch=="fixed":
        threshold_day=np.full(N_DAYS,50.0)
        p_month=fixed_month_prob(x,dates,threshold_day)
        pbar=np.nan
    elif branch=="equalized":
        threshold_day,pbar,pfix=equalized_thresholds(x,dates,50.0)
        p_month=fixed_month_prob(x,dates,threshold_day)
    else:
        raise ValueError(branch)
    weekly=weekly_from_daily(x,dates,threshold_day,p_month)
    z=residualize_weekly(weekly)
    nat,errn=native_surrogate_missing(x,dates,B,seed_rng(seed,station_id,branch,"native"),12)
    check_native_event_count_conservation(x,nat,dates,threshold_day)
    wnat=weekly_from_surrogates(nat,obs,dates,threshold_day,p_month)
    znat=residualize_weekly(wnat)
    idx,erri=iaaft_batch(z,B,seed_rng(seed,station_id,branch,"index"),30)
    out={"branch":branch,"n_obs_days":int(obs.sum()),"missing_days":int((~obs).sum()),
         "neff_realized":float(p_month.sum()**2/np.sum(p_month*p_month)),
         "p_month_mean":float(np.mean(p_month)),"pbar_equalized":pbar,
         "native_spectral_error_median":float(np.median(errn)),"index_spectral_error_median":float(np.median(erri))}
    out.update(summarize_stat(z,znat,idx,lambda a:fs(a,.90),"fs90"))
    out.update(summarize_stat(z,znat,idx,lambda a:fs(a,.85),"fs85"))
    out.update(summarize_stat(z,znat,idx,tdiff,"tdiff"))
    rn=out["rej_native_fs90"]; ri=out["rej_index_fs90"]
    out["paired_outcome"]="both" if rn and ri else "native_only" if rn else "index_only" if ri else "neither"
    exp=weekly_exposure_standardized(x,dates,threshold_day)
    out["exposure_sensitivity_defined"]=exp is not None
    if exp is not None:
        ze=residualize_weekly(exp); out["obs_fs90_exposure"]=float(fs(ze,.90)[0])
    return out


def run(args):
    outdir=Path(args.outdir); outdir.mkdir(parents=True,exist_ok=True)
    df,cohort=load_data(args.data,args.cohort)
    ids=sorted(cohort.station_id.astype(str).tolist())
    cohort_sha=hashlib.sha256(("\n".join(ids)+"\n").encode()).hexdigest()
    if cohort_sha != EXPECTED_COHORT_SHA:
        raise RuntimeError(f"frozen cohort hash mismatch: {cohort_sha} != {EXPECTED_COHORT_SHA}")
    recs=[]; start=time.time()
    for j,row in cohort.sort_values("station_id").reset_index(drop=True).iterrows():
        dates,x=prepare_station(df,row.short_code)
        for branch in (["fixed","equalized"] if args.equalized else ["fixed"]):
            r=run_branch(row.short_code,row.station_id,x,dates,args.B,args.seed,branch)
            r.update(station_id=row.station_id,short_code=row.short_code,station_name=row.station_name,
                     region=row.region,neff_frozen=float(row.neff_p_monthly_gt50))
            recs.append(r)
        print(f"[{j+1:02d}/{len(cohort)}] {row.station_id} elapsed={time.time()-start:.1f}s",flush=True)
    res=pd.DataFrame(recs)
    # BY adjustment within each branch and contract for primary FS90
    for branch,gidx in res.groupby("branch").groups.items():
        idx=list(gidx)
        res.loc[idx,"p_native_fs90_BY"]=by_adjust(res.loc[idx,"p_native_fs90"].to_numpy(float))
        res.loc[idx,"p_index_fs90_BY"]=by_adjust(res.loc[idx,"p_index_fs90"].to_numpy(float))
    res.to_csv(outdir/"SINCA_station_results_V1.csv",index=False)
    summaries=[]
    for branch,g in res.groupby("branch"):
        obs,p,ds=exact_region_signflip(g,"paired_outcome")
        rho,lo,hi=region_bootstrap_spearman(g,"neff_frozen","log_width_ratio",2000,args.seed+17) if "log_width_ratio" in g else (np.nan,np.nan,np.nan)
        # create log width ratio explicitly if not present
        lwr=np.log(g.width_ratio_fs90.astype(float))
        gg=g.copy(); gg["log_width_ratio"]=lwr
        rho,lo,hi=region_bootstrap_spearman(gg,"neff_frozen","log_width_ratio",2000,args.seed+17)
        summaries.append({
            "branch":branch,"n_stations":len(g),
            "index_reject_count":int(g.rej_index_fs90.sum()),"native_reject_count":int(g.rej_native_fs90.sum()),
            "index_only_count":int(np.sum(g.paired_outcome=="index_only")),"native_only_count":int(np.sum(g.paired_outcome=="native_only")),
            "both_count":int(np.sum(g.paired_outcome=="both")),"neither_count":int(np.sum(g.paired_outcome=="neither")),
            "region_signflip_observed_difference":obs,"region_signflip_one_sided_p":p,
            "spearman_neff_log_width_ratio":rho,"cluster_bootstrap_95_lo":lo,"cluster_bootstrap_95_hi":hi,
            "median_width_ratio_fs90":float(g.width_ratio_fs90.median()),
            "median_center_gap_fs90":float(g.center_gap_fs90.median()),
            "BY_native_discoveries":int(np.sum(g.p_native_fs90_BY<.05)),
            "BY_index_discoveries":int(np.sum(g.p_index_fs90_BY<.05)),
            "region_differences_json":json.dumps(ds),
        })
    sm=pd.DataFrame(summaries); sm.to_csv(outdir/"SINCA_across_station_summary_V1.csv",index=False)
    meta={"B":args.B,"base_seed":args.seed,"data":str(args.data),"cohort":str(args.cohort),
          "analysis_window":[str(START.date()),str(END.date())],"n_weeks":N_WEEKS,
          "equalized_branch":bool(args.equalized),"elapsed_seconds":time.time()-start}
    (outdir/"RUN_MANIFEST_V1.json").write_text(json.dumps(meta,indent=2))
    print(sm.to_string(index=False))


def smoke_test():
    # synthetic 8-year daily PM2.5-like series with seasonality and short memory
    dates=pd.date_range(START,END,freq="D"); n=len(dates); rng=np.random.default_rng(123)
    e=rng.normal(size=n); ar=np.empty(n); ar[0]=e[0]
    for i in range(1,n): ar[i]=.6*ar[i-1]+math.sqrt(1-.6**2)*e[i]
    doy=dates.dayofyear.to_numpy(); x=25+18*np.cos(2*np.pi*(doy-190)/365.2425)+8*ar
    miss=rng.random(n)<.02; x[miss]=np.nan
    out=run_branch("XX","TEST/001",x,dates,B=9,seed=1234,branch="fixed")
    needed=["obs_fs90","p_native_fs90","p_index_fs90","paired_outcome","width_ratio_fs90"]
    assert all(k in out and np.isfinite(out[k]) if k!="paired_outcome" else k in out for k in needed)
    print(json.dumps({k:out[k] for k in needed},indent=2,default=float))


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--data",default="/mnt/data/sinca_csv_full/sinca_daily_station.csv")
    ap.add_argument("--cohort",default="/mnt/data/SINCA_PM25_AOAS_FROZEN_COHORT_V1/FROZEN_COHORT_V1.csv")
    ap.add_argument("--outdir",default="SINCA_PM25_AOAS_RESULTS_V1")
    ap.add_argument("--B",type=int,default=B_DEFAULT)
    ap.add_argument("--seed",type=int,default=BASE_SEED_DEFAULT)
    ap.add_argument("--equalized",action=argparse.BooleanOptionalAction,default=True)
    ap.add_argument("--smoke-test",action="store_true")
    a=ap.parse_args()
    if a.smoke_test: smoke_test()
    else: run(a)
