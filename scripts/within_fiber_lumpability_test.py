#!/usr/bin/env python3
"""Within-fiber empirical strong-lumpability diagnostic for native surrogate kernel.

Construct x_tilde from a Gaussian AR(1) trajectory x by permuting exact amplitudes
*within each calendar-month x wet/dry indicator class*.  Therefore x and x_tilde have:
  - identical monthly threshold indicator sequence,
  - identical annual count sequence O(x)=O(x_tilde),
  - identical calendar-month empirical value multisets, means and variances,
while generally having different standardized monthly Fourier-amplitude targets.

If the native constrained-surrogate kernel were strongly lumpable with respect to
O (or P=D o O), the pushforward native-surrogate distributions from x and x_tilde
would be identical.  We compare their projected distributions and calibrate Monte
Carlo distance against two independent ensembles generated from the same x.
"""
from __future__ import annotations
import argparse, json, math, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.special import expit
from scipy.stats import binomtest, ks_2samp, norm, wasserstein_distance

N_YEARS_SIM=252; N_YEARS_ANALYSIS=251; MONTHS=12
N_MONTHS=N_YEARS_SIM*MONTHS; N_ANALYSIS_MONTHS=N_YEARS_ANALYSIS*MONTHS; Q=.90


def p_profile(lam,K=3.0):
    m=np.arange(1,13); c=np.cos(2*np.pi*(m-8)/12)
    a=brentq(lambda a: expit(a+lam*c).sum()-K,-30,30)
    return expit(a+lam*c)

def neff(p):
    p=np.asarray(p,float); return float(p.sum()**2/np.sum(p*p))

def generate_ar1(phi,n,rng):
    sig=math.sqrt(1-phi*phi); x=np.empty(n); x[0]=rng.normal(); eps=rng.normal(size=n-1)
    for t in range(1,n): x[t]=phi*x[t-1]+sig*eps[t-1]
    return x

def reshape_months(x):
    if x.ndim==1: return x[:N_ANALYSIS_MONTHS].reshape(N_YEARS_ANALYSIS,MONTHS)
    return x[:,:N_ANALYSIS_MONTHS].reshape(x.shape[0],N_YEARS_ANALYSIS,MONTHS)

def annual_count(x,thr):
    xx=reshape_months(x)
    if x.ndim==1: return (xx<thr[None,:]).sum(axis=1).astype(float)
    return (xx<thr[None,None,:]).sum(axis=2).astype(float)

def cubic_detrend(y):
    n=y.shape[-1]; t=np.linspace(-1,1,n); X=np.column_stack([np.ones(n),t,t*t,t*t*t])
    if y.ndim==1:
        b,*_=np.linalg.lstsq(X,y,rcond=None); return y-X@b
    b,*_=np.linalg.lstsq(X,y.T,rcond=None); return y-(X@b).T

def fs_one(z):
    u=np.quantile(z,Q); idx=np.flatnonzero(z>u)
    if idx.size<2: return 1.0
    T=np.diff(idx).astype(float); nT=T.size
    if T.max()<=2:
        den=nT*np.sum(T*T); raw=1.0 if den<=0 else 2*(T.sum()**2)/den
    else:
        S=T-1; den=nT*np.sum(S*(S-1)); raw=1.0 if den<=0 else 2*(S.sum()**2)/den
    return float(np.clip(raw,0,1))

def fs(arr):
    if arr.ndim==1: return np.array([fs_one(arr)])
    return np.fromiter((fs_one(r) for r in arr),float,count=arr.shape[0])

def lag1_one(y):
    a=np.asarray(y[:-1],float); b=np.asarray(y[1:],float)
    sa=a.std(); sb=b.std()
    if sa<1e-12 or sb<1e-12: return 0.0
    return float(np.corrcoef(a,b)[0,1])

def lag1(arr):
    if arr.ndim==1: return np.array([lag1_one(arr)])
    return np.fromiter((lag1_one(r) for r in arr),float,count=arr.shape[0])

def tdiff_one(y): return float(np.mean(np.diff(y)**2))
def tdiff(arr):
    if arr.ndim==1: return np.array([tdiff_one(arr)])
    return np.fromiter((tdiff_one(r) for r in arr),float,count=arr.shape[0])

def lowfreq_one(y,frac=.10):
    z=np.asarray(y,float)-np.mean(y); p=np.abs(np.fft.rfft(z))**2
    if p.size<=2 or p[1:].sum()<=1e-15: return 0.0
    k=max(1,int(np.floor(frac*(p.size-1))))
    return float(p[1:1+k].sum()/p[1:].sum())
def lowfreq(arr):
    if arr.ndim==1: return np.array([lowfreq_one(arr)])
    return np.fromiter((lowfreq_one(r) for r in arr),float,count=arr.shape[0])

def native_surrogate_batch(x,B,rng,max_iter=12):
    n=x.size; idxs=[np.arange(m,n,12) for m in range(12)]
    mu=np.array([x[i].mean() for i in idxs]); sd=np.array([x[i].std(ddof=0) for i in idxs]); sd=np.where(sd>0,sd,1)
    mut=np.tile(mu,N_YEARS_SIM); sdt=np.tile(sd,N_YEARS_SIM); zobs=(x-mut)/sdt; target=np.abs(np.fft.rfft(zobs))
    cur=np.empty((B,n)); vals=[]
    for idx in idxs:
        v=np.sort(x[idx]); vals.append(v); o=np.argsort(rng.random((B,idx.size)),axis=1); cur[:,idx]=v[o]
    rr=np.arange(B)[:,None]
    for _ in range(max_iter):
        z=(cur-mut[None,:])/sdt[None,:]; F=np.fft.rfft(z,axis=1); phase=np.exp(1j*np.angle(F)); scores=np.fft.irfft(target[None,:]*phase,n=n,axis=1)
        nxt=np.empty_like(cur)
        for idx,v in zip(idxs,vals):
            o=np.argsort(scores[:,idx],axis=1); block=np.empty((B,idx.size)); block[rr,o]=v[None,:]; nxt[:,idx]=block
        cur=nxt
    z=(cur-mut[None,:])/sdt[None,:]; amp=np.abs(np.fft.rfft(z,axis=1)); den=np.mean(target**2)+1e-15
    err=np.mean((amp-target[None,:])**2,axis=1)/den
    return cur,err,target

def within_state_month_permutation(x,thr,rng):
    """Permute exact values within each calendar-month and same indicator class."""
    xt=np.array(x,copy=True)
    for m in range(12):
        idx=np.arange(m,x.size,12)
        status=x[idx] < thr[m]
        for val in (False,True):
            ids=idx[status==val]
            if ids.size>1:
                xt[ids]=x[rng.permutation(ids)]
    return xt

def checks(x,xt,thr):
    ind=(x < np.tile(thr,N_YEARS_SIM)); indt=(xt < np.tile(thr,N_YEARS_SIM))
    c=annual_count(x,thr); ct=annual_count(xt,thr)
    max_multiset=0.0
    for m in range(12):
        idx=np.arange(m,x.size,12)
        max_multiset=max(max_multiset,float(np.max(np.abs(np.sort(x[idx])-np.sort(xt[idx])))))
    return dict(indicator_mismatches=int(np.sum(ind!=indt)),annual_count_max_abs=float(np.max(np.abs(c-ct))),month_multiset_max_abs=max_multiset)

def target_spectral_distance(t1,t2):
    den=np.mean(t1*t1)+1e-15
    return float(np.mean((t1-t2)**2)/den)

def p_lower(obs,null): return float((1+np.sum(null<=obs))/(len(null)+1))

def compare_feature(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    pooled=np.std(np.r_[a,b],ddof=1); scale=pooled if pooled>1e-12 else 1.0
    ks=ks_2samp(a,b,method='auto')
    return dict(wass=float(wasserstein_distance(a,b)),wass_scaled=float(wasserstein_distance(a,b)/scale),ks=float(ks.statistic),ks_p=float(ks.pvalue),mean_diff=float(np.mean(b)-np.mean(a)),sd_ratio=float(np.std(b,ddof=1)/(np.std(a,ddof=1)+1e-15)))

def seed_rng(seed,*parts): return np.random.default_rng(np.random.SeedSequence([seed,*parts]))

def run(phis,reps,B,seed,lams,outdir):
    outdir=Path(outdir); outdir.mkdir(parents=True,exist_ok=True); start=time.time(); records=[]
    profiles={lam:p_profile(lam) for lam in lams}; thresholds={lam:norm.ppf(profiles[lam]) for lam in lams}
    seq=np.random.SeedSequence(seed).spawn(len(phis)*reps); si=0
    features={'ferro_segers':fs,'lag1':lag1,'tdiff_raw':tdiff,'lowfreq_cubic':lowfreq}
    for pi,phi in enumerate(phis):
        for rep in range(reps):
            x=generate_ar1(phi,N_MONTHS,np.random.default_rng(seq[si])); si+=1
            for li,lam in enumerate(lams):
                thr=thresholds[lam]
                xt=within_state_month_permutation(x,thr,seed_rng(seed,301,pi,rep,li))
                ck=checks(x,xt,thr)
                # Three independent ensembles: x-A, x-B baseline replicate, and fiber mate xt-C.
                sxA,eA,tA=native_surrogate_batch(x,B,seed_rng(seed,401,pi,rep,li,0))
                sxB,eB,tB=native_surrogate_batch(x,B,seed_rng(seed,401,pi,rep,li,1))
                stC,eC,tC=native_surrogate_batch(xt,B,seed_rng(seed,401,pi,rep,li,2))
                raw_obs=annual_count(x,thr); raw_t=annual_count(xt,thr)
                zobs=cubic_detrend(raw_obs)
                rawA=annual_count(sxA,thr); rawB=annual_count(sxB,thr); rawC=annual_count(stC,thr)
                cubA,cubB,cubC=cubic_detrend(rawA),cubic_detrend(rawB),cubic_detrend(rawC)
                feat_arrays={
                    'ferro_segers':(fs(cubA),fs(cubB),fs(cubC),float(fs(zobs)[0])),
                    'lag1':(lag1(cubA),lag1(cubB),lag1(cubC),float(lag1(zobs)[0])),
                    'tdiff_raw':(tdiff(rawA),tdiff(rawB),tdiff(rawC),float(tdiff(raw_obs)[0])),
                    'lowfreq_cubic':(lowfreq(cubA),lowfreq(cubB),lowfreq(cubC),float(lowfreq(zobs)[0])),
                }
                specdist=target_spectral_distance(tA,tC)
                for fname,(A,Bb,C,obsf) in feat_arrays.items():
                    base=compare_feature(A,Bb); fib=compare_feature(A,C)
                    pA=p_lower(obsf,A) if fname=='ferro_segers' else np.nan
                    pB=p_lower(obsf,Bb) if fname=='ferro_segers' else np.nan
                    pC=p_lower(obsf,C) if fname=='ferro_segers' else np.nan
                    records.append(dict(phi=phi,rep=rep,lambda_value=lam,Neff=neff(profiles[lam]),feature=fname,
                        obs_feature=obsf,indicator_mismatches=ck['indicator_mismatches'],annual_count_max_abs=ck['annual_count_max_abs'],month_multiset_max_abs=ck['month_multiset_max_abs'],
                        native_target_spectral_distance=specdist,
                        same_wass=base['wass'],same_wass_scaled=base['wass_scaled'],same_ks=base['ks'],same_ks_p=base['ks_p'],same_mean_diff=base['mean_diff'],same_sd_ratio=base['sd_ratio'],
                        fiber_wass=fib['wass'],fiber_wass_scaled=fib['wass_scaled'],fiber_ks=fib['ks'],fiber_ks_p=fib['ks_p'],fiber_mean_diff=fib['mean_diff'],fiber_sd_ratio=fib['sd_ratio'],
                        distance_ratio=(fib['wass_scaled']/(base['wass_scaled']+1e-15)),
                        p_xA=pA,p_xB=pB,p_xtilde=pC,abs_p_same=(abs(pA-pB) if fname=='ferro_segers' else np.nan),abs_p_fiber=(abs(pA-pC) if fname=='ferro_segers' else np.nan),
                        rej_xA=(pA<.05 if fname=='ferro_segers' else False),rej_xB=(pB<.05 if fname=='ferro_segers' else False),rej_xtilde=(pC<.05 if fname=='ferro_segers' else False),
                        median_specerr_xA=float(np.median(eA)),median_specerr_xB=float(np.median(eB)),median_specerr_xtilde=float(np.median(eC))))
            if (rep+1)%max(1,reps//5)==0: print(f'phi={phi:.2f} {rep+1}/{reps} elapsed={time.time()-start:.1f}s',flush=True)
    df=pd.DataFrame(records); df.to_csv(outdir/'within_fiber_trajectory_results.csv',index=False)
    rows=[]
    for (lam,feat,phi),g in df.groupby(['lambda_value','feature','phi'],sort=False):
        rows.append(summary_row(g,lam,feat,phi))
    for (lam,feat),g in df.groupby(['lambda_value','feature'],sort=False):
        rows.append(summary_row(g,lam,feat,'pooled'))
    sm=pd.DataFrame(rows); sm.to_csv(outdir/'within_fiber_summary.csv',index=False)
    meta=dict(seed=seed,reps_per_phi=reps,B=B,phi_values=phis,lambda_values=lams,Neff={str(l):neff(profiles[l]) for l in lams},years_simulated=N_YEARS_SIM,years_analyzed=N_YEARS_ANALYSIS,native_iterations=12,elapsed_seconds=time.time()-start,
              construction='Permute exact x values within each calendar-month x threshold-indicator class. Preserves indicator sequence and each calendar-month empirical multiset exactly; changes target standardized monthly spectrum.',
              interpretation='If fiber-pair pushforward feature distributions differ beyond same-x independent-ensemble Monte Carlo baseline, this falsifies strong lumpability for the implemented native kernel with respect to the threshold+annual-count operator (and deterministic cubic detrending).')
    (outdir/'within_fiber_metadata.json').write_text(json.dumps(meta,indent=2))
    print(sm[sm.phi.astype(str).eq('pooled')].to_string(index=False),flush=True)

def summary_row(g,lam,feat,phi):
    # Compare fiber vs same Monte Carlo distances trajectory-wise.
    diff=g.fiber_wass_scaled.to_numpy()-g.same_wass_scaled.to_numpy()
    gt=float(np.mean(diff>0)); med_diff=float(np.median(diff))
    # Sign test for fiber distance > same baseline.
    nz=diff[np.abs(diff)>1e-15]; p_sign=1.0 if nz.size==0 else float(binomtest(int(np.sum(nz>0)),int(nz.size),.5,alternative='greater').pvalue)
    out=dict(lambda_value=lam,Neff=float(g.Neff.iloc[0]),feature=feat,phi=phi,n=len(g),
        median_target_spectral_distance=float(g.native_target_spectral_distance.median()),
        median_same_wass_scaled=float(g.same_wass_scaled.median()),median_fiber_wass_scaled=float(g.fiber_wass_scaled.median()),median_distance_ratio=float(g.distance_ratio.median()),
        median_fiber_minus_same=float(med_diff),fraction_fiber_distance_gt_same=gt,sign_test_p=float(p_sign),
        median_same_ks=float(g.same_ks.median()),median_fiber_ks=float(g.fiber_ks.median()),fraction_fiber_KS_p_lt_005=float(np.mean(g.fiber_ks_p<.05)),fraction_same_KS_p_lt_005=float(np.mean(g.same_ks_p<.05)),
        median_abs_p_same=(float(g.abs_p_same.median()) if feat=='ferro_segers' else np.nan),median_abs_p_fiber=(float(g.abs_p_fiber.median()) if feat=='ferro_segers' else np.nan),
        reject_discord_same=(float(np.mean(g.rej_xA.to_numpy()!=g.rej_xB.to_numpy())) if feat=='ferro_segers' else np.nan),reject_discord_fiber=(float(np.mean(g.rej_xA.to_numpy()!=g.rej_xtilde.to_numpy())) if feat=='ferro_segers' else np.nan),
        max_indicator_mismatches=int(g.indicator_mismatches.max()),max_annual_count_abs=float(g.annual_count_max_abs.max()),max_month_multiset_abs=float(g.month_multiset_max_abs.max()),
        median_specerr_x=float(g.median_specerr_xA.median()),median_specerr_xtilde=float(g.median_specerr_xtilde.median()))
    return out

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--reps',type=int,default=40); ap.add_argument('--B',type=int,default=249); ap.add_argument('--seed',type=int,default=20260910); ap.add_argument('--phis',type=float,nargs='+',default=[.2,.5,.8]); ap.add_argument('--lams',type=float,nargs='+',default=[0.,4.]); ap.add_argument('--outdir',default='/mnt/data/within_fiber_lumpability')
    a=ap.parse_args(); run(a.phis,a.reps,a.B,a.seed,a.lams,a.outdir)
