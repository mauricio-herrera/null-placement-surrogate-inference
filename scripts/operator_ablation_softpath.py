#!/usr/bin/env python3
"""Refined Gaussian operator-ablation via a soft-threshold path.

This experiment keeps the annual output length fixed and continuously moves the
observation operator from a smooth transformation toward the hard monthly
indicator count:

    S_y(tau) = sum_m expit((u_m - X_{y,m}) / tau),  tau > 0,
    C_y      = sum_m 1{X_{y,m} < u_m},               tau = 0 (hard limit).

The same observed Gaussian AR(1) trajectories are compared under:
  oracle generative null, native monthly constrained surrogate + operator,
  and index-resolution IAAFT after the operator.

Soft path is evaluated after the manuscript's cubic detrending. Raw annual mean
and raw hard-count anchors are also included to isolate detrending.
"""
from __future__ import annotations
import argparse, json, math, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.special import expit
from scipy.stats import binomtest, norm, wasserstein_distance

N_YEARS_SIM=252; N_YEARS_ANALYSIS=251; MONTHS=12
N_MONTHS=N_YEARS_SIM*MONTHS; N_ANALYSIS_MONTHS=N_YEARS_ANALYSIS*MONTHS; Q=0.90


def p_profile(lam, K=3.0):
    m=np.arange(1,13); c=np.cos(2*np.pi*(m-8)/12)
    a=brentq(lambda a: expit(a+lam*c).sum()-K,-30,30)
    return expit(a+lam*c)

def neff(p):
    p=np.asarray(p,float); return float(p.sum()**2/np.sum(p*p))

def generate_ar1(phi,n,rng,B=None):
    sig=math.sqrt(1-phi*phi)
    if B is None:
        x=np.empty(n); x[0]=rng.normal(); eps=rng.normal(size=n-1)
        for t in range(1,n): x[t]=phi*x[t-1]+sig*eps[t-1]
        return x
    x=np.empty((B,n)); x[:,0]=rng.normal(size=B); eps=rng.normal(size=(B,n-1))
    for t in range(1,n): x[:,t]=phi*x[:,t-1]+sig*eps[:,t-1]
    return x

def reshape_months(x):
    if x.ndim==1: return x[:N_ANALYSIS_MONTHS].reshape(N_YEARS_ANALYSIS,MONTHS)
    return x[:,:N_ANALYSIS_MONTHS].reshape(x.shape[0],N_YEARS_ANALYSIS,MONTHS)

def annual_mean(x): return reshape_months(x).mean(axis=-1)
def annual_count(x,thr):
    xx=reshape_months(x)
    if x.ndim==1: return (xx<thr[None,:]).sum(axis=1).astype(float)
    return (xx<thr[None,None,:]).sum(axis=2).astype(float)
def annual_soft_count(x,thr,tau):
    xx=reshape_months(x)
    if x.ndim==1: return expit((thr[None,:]-xx)/tau).sum(axis=1)
    return expit((thr[None,None,:]-xx)/tau).sum(axis=2)
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
def ntv_one(y):
    s=float(np.std(y)); return 0.0 if s<=1e-15 else float(np.mean(np.abs(np.diff(y)))/s)
def ntv(arr):
    if arr.ndim==1: return np.array([ntv_one(arr)])
    return np.fromiter((ntv_one(r) for r in arr),float,count=arr.shape[0])

def spectral_error_batch(y,target_amp):
    amp=np.abs(np.fft.rfft(y-y.mean(axis=1,keepdims=True),axis=1)); den=np.mean(target_amp**2)+1e-15
    return np.mean((amp-target_amp[None,:])**2,axis=1)/den

def iaaft_batch(z,B,rng,max_iter=30):
    z=np.asarray(z,float); n=z.size; zs=np.sort(z); target=np.abs(np.fft.rfft(z-z.mean()))
    order0=np.argsort(rng.random((B,n)),axis=1); cur=np.empty((B,n)); cur[np.arange(B)[:,None],order0]=zs[None,:]
    for _ in range(max_iter):
        F=np.fft.rfft(cur-cur.mean(axis=1,keepdims=True),axis=1); phase=np.exp(1j*np.angle(F))
        scores=np.fft.irfft(target[None,:]*phase,n=n,axis=1); order=np.argsort(scores,axis=1)
        nxt=np.empty_like(cur); nxt[np.arange(B)[:,None],order]=zs[None,:]; cur=nxt
    return cur,spectral_error_batch(cur,target)

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
    return cur,np.mean((amp-target[None,:])**2,axis=1)/den

def p_lower(obs,null): return float((1+np.sum(null<=obs))/(len(null)+1))
def eval_stat(obs,nat,ora,idx,name,fun):
    o=float(fun(obs)[0]); n=fun(nat); r=fun(ora); i=fun(idx)
    ps=[p_lower(o,a) for a in (r,n,i)]; rsd=float(np.std(r,ddof=1)); scale=rsd if rsd>1e-12 else 1.0
    return dict(statistic=name,obs_stat=o,p_oracle=ps[0],p_native=ps[1],p_index=ps[2],
                rej_oracle=ps[0]<.05,rej_native=ps[1]<.05,rej_index=ps[2]<.05,
                null_mean_oracle=float(np.mean(r)),null_sd_oracle=rsd,
                null_mean_native=float(np.mean(n)),null_sd_native=float(np.std(n,ddof=1)),
                null_mean_index=float(np.mean(i)),null_sd_index=float(np.std(i,ddof=1)),
                wasserstein_native_oracle=float(wasserstein_distance(n,r)),wasserstein_index_oracle=float(wasserstein_distance(i,r)),
                wasserstein_native_oracle_scaled=float(wasserstein_distance(n,r)/scale),wasserstein_index_oracle_scaled=float(wasserstein_distance(i,r)/scale))

def seed_rng(seed,*parts): return np.random.default_rng(np.random.SeedSequence([seed,*parts]))

def run(phis,reps,B,seed,taus,outdir):
    outdir=Path(outdir); outdir.mkdir(parents=True,exist_ok=True); start=time.time()
    p0,p4=p_profile(0),p_profile(4); u0,u4=norm.ppf(p0),norm.ppf(p4)
    t_seeds=np.random.SeedSequence(seed).spawn(len(phis)*reps); si=0; recs=[]
    condition_labels=[]
    for lam in (0,4):
        for tau in taus: condition_labels.append((f"soft_lam{lam}_tau{tau:g}_cubic",lam,tau,"soft"))
        condition_labels.append((f"count_lam{lam}_cubic",lam,0.0,"hard"))
    anchors=[("mean_raw",None,None,"mean_raw"),("mean_cubic",None,None,"mean_cubic"),
             ("count_lam0_raw",0,0.0,"hard_raw"),("count_lam4_raw",4,0.0,"hard_raw")]
    allconds=anchors+condition_labels
    for pi,phi in enumerate(phis):
        for rep in range(reps):
            rng=np.random.default_rng(t_seeds[si]); si+=1; x=generate_ar1(phi,N_MONTHS,rng)
            nx,nspec=native_surrogate_batch(x,B,seed_rng(seed,111,pi,rep)); ox=generate_ar1(phi,N_MONTHS,seed_rng(seed,222,pi,rep),B=B)
            mean_o,mean_n,mean_r=annual_mean(x),annual_mean(nx),annual_mean(ox)
            cache={"mean_raw":(mean_o,mean_n,mean_r),"mean_cubic":(cubic_detrend(mean_o),cubic_detrend(mean_n),cubic_detrend(mean_r))}
            for lam,thr in [(0,u0),(4,u4)]:
                co,cn,cr=annual_count(x,thr),annual_count(nx,thr),annual_count(ox,thr)
                cache[f"count_lam{lam}_raw"]=(co,cn,cr); cache[f"count_lam{lam}_cubic"]=(cubic_detrend(co),cubic_detrend(cn),cubic_detrend(cr))
                for tau in taus:
                    so,sn,sr=annual_soft_count(x,thr,tau),annual_soft_count(nx,thr,tau),annual_soft_count(ox,thr,tau)
                    cache[f"soft_lam{lam}_tau{tau:g}_cubic"]=(cubic_detrend(so),cubic_detrend(sn),cubic_detrend(sr))
            for ci,(cname,lam,tau,ctype) in enumerate(allconds):
                obs,nat,ora=cache[cname]; idx,ispec=iaaft_batch(obs,B,seed_rng(seed,777,pi,rep,ci))
                stats=[("ntv",ntv)]
                if ctype!="hard_raw": stats.append(("ferro_segers",fs))
                for sname,sfun in stats:
                    d=eval_stat(obs,nat,ora,idx,sname,sfun); d.update(phi=phi,rep=rep,condition=cname,condition_type=ctype,
                        lambda_value=lam,tau=tau,native_spectral_error_median=float(np.median(nspec)),index_spectral_error_median=float(np.median(ispec)))
                    recs.append(d)
            if (rep+1)%max(1,reps//5)==0: print(f"phi={phi:.2f} {rep+1}/{reps} elapsed={time.time()-start:.1f}s",flush=True)
    df=pd.DataFrame(recs); df.to_csv(outdir/'operator_ablation_softpath_trajectory_results.csv',index=False)
    def agg(g,cond,stat,phi):
        A=g.rej_index.to_numpy(bool); N=g.rej_native.to_numpy(bool); O=g.rej_oracle.to_numpy(bool); ni=int(np.sum(A&~N)); nn=int(np.sum(N&~A)); d=ni+nn
        pp=1.0 if d==0 else float(binomtest(min(ni,nn),d,.5).pvalue)
        den=g.null_sd_oracle.replace(0,np.nan)
        return dict(condition=cond,condition_type=g.condition_type.iloc[0],lambda_value=g.lambda_value.iloc[0],tau=g.tau.iloc[0],statistic=stat,phi=phi,n=len(g),
            oracle_reject_rate=float(O.mean()),native_reject_rate=float(N.mean()),index_reject_rate=float(A.mean()),index_only_vs_native=float(np.mean(A&~N)),native_only_vs_index=float(np.mean(N&~A)),
            index_only_count=ni,native_only_count=nn,paired_exact_p=pp,index_only_vs_oracle=float(np.mean(A&~O)),oracle_only_vs_index=float(np.mean(O&~A)),
            median_null_sd_oracle=float(g.null_sd_oracle.median()),median_null_sd_native=float(g.null_sd_native.median()),median_null_sd_index=float(g.null_sd_index.median()),
            median_sd_ratio_native_oracle=float(np.nanmedian(g.null_sd_native/den)),median_sd_ratio_index_oracle=float(np.nanmedian(g.null_sd_index/den)),
            median_wasserstein_native_oracle_scaled=float(g.wasserstein_native_oracle_scaled.median()),median_wasserstein_index_oracle_scaled=float(g.wasserstein_index_oracle_scaled.median()),
            median_index_spectral_error=float(g.index_spectral_error_median.median()),median_native_spectral_error=float(g.native_spectral_error_median.median()))
    rows=[]
    for (c,s,p),g in df.groupby(['condition','statistic','phi'],sort=False): rows.append(agg(g,c,s,p))
    for (c,s),g in df.groupby(['condition','statistic'],sort=False): rows.append(agg(g,c,s,'pooled'))
    sm=pd.DataFrame(rows); sm.to_csv(outdir/'operator_ablation_softpath_summary.csv',index=False)
    meta=dict(seed=seed,reps_per_phi=reps,B=B,phi_values=phis,taus=taus,Neff_lambda0=neff(p0),Neff_lambda4=neff(p4),p_lambda0=p0.tolist(),p_lambda4=p4.tolist(),
              years_simulated=N_YEARS_SIM,years_analyzed=N_YEARS_ANALYSIS,iaaft_iterations=30,native_iterations=12,
              attainable_strict_alpha=math.floor(.05*(B+1)-1e-12)/(B+1),elapsed_seconds=time.time()-start,
              soft_operator='sum_m expit((u_m-X_ym)/tau); tau->0 gives hard count',
              note='Index IAAFT uses fixed 30 iterations because manuscript source does not specify numerical convergence tolerance.')
    (outdir/'operator_ablation_softpath_metadata.json').write_text(json.dumps(meta,indent=2))
    print(sm[(sm.phi.astype(str)=='pooled') & (sm.statistic=='ferro_segers')].to_string(index=False),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--reps',type=int,default=100); ap.add_argument('--B',type=int,default=249); ap.add_argument('--seed',type=int,default=20260910)
    ap.add_argument('--phis',type=float,nargs='+',default=[.2,.5,.8]); ap.add_argument('--taus',type=float,nargs='+',default=[2.0,1.0,.5,.25,.1]); ap.add_argument('--outdir',default='operator_ablation_softpath_results')
    a=ap.parse_args(); run(a.phis,a.reps,a.B,a.seed,a.taus,a.outdir)
