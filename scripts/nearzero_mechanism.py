#!/usr/bin/env python3
from __future__ import annotations
import argparse, math, time, json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import binomtest, spearmanr
sys.path.insert(0,'/mnt/data/operator_ablation_src')
import operator_ablation_softpath as base

TAUS_DEFAULT=[0.1,0.05,0.025,0.01,0.005,0.002,0.001]

def pearson(a,b):
    a=np.asarray(a); b=np.asarray(b)
    if np.std(a)<1e-15 or np.std(b)<1e-15: return np.nan
    return float(np.corrcoef(a,b)[0,1])

def sim_metrics(soft_raw, hard_raw, soft_det, hard_det):
    d=soft_det-hard_det
    sd=float(np.std(hard_det))
    return dict(
        corr_raw_soft_hard=pearson(soft_raw,hard_raw),
        corr_det_soft_hard=pearson(soft_det,hard_det),
        rmse_det_soft_hard=float(np.sqrt(np.mean(d*d))),
        nrmse_det_soft_hard=float(np.sqrt(np.mean(d*d))/(sd+1e-15)),
        maxabs_det_soft_hard=float(np.max(np.abs(d))),
        raw_unique_soft=int(np.unique(soft_raw).size),
        raw_unique_hard=int(np.unique(hard_raw).size),
    )

def eval_fs(obs,nat,ora,idx):
    o=float(base.fs(obs)[0]); n=base.fs(nat); r=base.fs(ora); i=base.fs(idx)
    po,pn,pi=[base.p_lower(o,a) for a in (r,n,i)]
    rsd=float(np.std(r,ddof=1)); nsd=float(np.std(n,ddof=1)); isd=float(np.std(i,ddof=1))
    return dict(obs_theta=o,p_oracle=po,p_native=pn,p_index=pi,
                rej_oracle=po<.05,rej_native=pn<.05,rej_index=pi<.05,
                null_sd_oracle=rsd,null_sd_native=nsd,null_sd_index=isd,
                sd_ratio_native_oracle=nsd/(rsd+1e-15),sd_ratio_index_oracle=isd/(rsd+1e-15),
                null_mean_oracle=float(np.mean(r)),null_mean_native=float(np.mean(n)),null_mean_index=float(np.mean(i)))

def run(phis,reps,B,seed,taus,outdir,iters=30):
    outdir=Path(outdir); outdir.mkdir(parents=True,exist_ok=True); start=time.time()
    p4=base.p_profile(4); thr=base.norm.ppf(p4)
    tseeds=np.random.SeedSequence(seed).spawn(len(phis)*reps); si=0; rec=[]
    for pi,phi in enumerate(phis):
        for rep in range(reps):
            rng=np.random.default_rng(tseeds[si]); si+=1
            x=base.generate_ar1(phi,base.N_MONTHS,rng)
            nx,nspec=base.native_surrogate_batch(x,B,base.seed_rng(seed,111,pi,rep))
            ox=base.generate_ar1(phi,base.N_MONTHS,base.seed_rng(seed,222,pi,rep),B=B)
            hard_o_raw=base.annual_count(x,thr); hard_n_raw=base.annual_count(nx,thr); hard_r_raw=base.annual_count(ox,thr)
            hard_o=base.cubic_detrend(hard_o_raw); hard_n=base.cubic_detrend(hard_n_raw); hard_r=base.cubic_detrend(hard_r_raw)
            conds=[]
            for tau in taus:
                so_raw=base.annual_soft_count(x,thr,tau); sn_raw=base.annual_soft_count(nx,thr,tau); sr_raw=base.annual_soft_count(ox,thr,tau)
                so=base.cubic_detrend(so_raw); sn=base.cubic_detrend(sn_raw); sr=base.cubic_detrend(sr_raw)
                conds.append((tau,so,sn,sr,sim_metrics(so_raw,hard_o_raw,so,hard_o)))
            conds.append((0.0,hard_o,hard_n,hard_r,sim_metrics(hard_o_raw,hard_o_raw,hard_o,hard_o)))
            for ci,(tau,obs,nat,ora,sim) in enumerate(conds):
                idx,ispec=base.iaaft_batch(obs,B,base.seed_rng(seed,777,pi,rep,ci),max_iter=iters)
                d=eval_fs(obs,nat,ora,idx)
                d.update(phi=phi,rep=rep,tau=tau,condition='hard' if tau==0 else f'soft_tau{tau:g}',
                         index_spectral_error_median=float(np.median(ispec)),
                         native_spectral_error_median=float(np.median(nspec)),**sim)
                rec.append(d)
            if (rep+1)%max(1,reps//5)==0:
                print(f'phi={phi} rep={rep+1}/{reps}, elapsed={time.time()-start:.1f}s',flush=True)
    df=pd.DataFrame(rec); df.to_csv(outdir/'nearzero_trajectory_results.csv',index=False)
    rows=[]
    for tau,g in df.groupby('tau',sort=False):
        A=g.rej_index.to_numpy(bool); N=g.rej_native.to_numpy(bool); O=g.rej_oracle.to_numpy(bool)
        ni=int(np.sum(A&~N)); nn=int(np.sum(N&~A)); disc=ni+nn
        rows.append(dict(tau=tau,n=len(g),oracle_reject_rate=float(O.mean()),native_reject_rate=float(N.mean()),index_reject_rate=float(A.mean()),
                         index_only_rate=float(np.mean(A&~N)),native_only_rate=float(np.mean(N&~A)),index_only_count=ni,native_only_count=nn,
                         paired_exact_p=1.0 if disc==0 else float(binomtest(min(ni,nn),disc,.5).pvalue),
                         median_sd_ratio_index_oracle=float(g.sd_ratio_index_oracle.median()),median_sd_ratio_native_oracle=float(g.sd_ratio_native_oracle.median()),
                         median_index_spectral_error=float(g.index_spectral_error_median.median()),median_native_spectral_error=float(g.native_spectral_error_median.median()),
                         median_corr_det_soft_hard=float(g.corr_det_soft_hard.median()),median_nrmse_det_soft_hard=float(g.nrmse_det_soft_hard.median()),
                         median_raw_unique_soft=float(g.raw_unique_soft.median()),median_raw_unique_hard=float(g.raw_unique_hard.median()),
                         median_null_mean_diff_index_oracle=float((g.null_mean_index-g.null_mean_oracle).median())))
    sm=pd.DataFrame(rows).sort_values('tau',ascending=False); sm.to_csv(outdir/'nearzero_summary.csv',index=False)
    # pooled mechanism correlations across all soft + hard trajectory-condition points
    dd=df.copy(); dd['log10_specerr']=np.log10(dd.index_spectral_error_median.clip(lower=1e-12))
    corr={}
    for y in ['sd_ratio_index_oracle','p_index','corr_det_soft_hard','nrmse_det_soft_hard']:
        rr=spearmanr(dd.log10_specerr,dd[y],nan_policy='omit')
        corr[y]={'rho':float(rr.statistic),'p':float(rr.pvalue)}
    meta=dict(seed=seed,phis=phis,reps_per_phi=reps,B=B,taus=taus,iaaft_iterations=iters,lambda_value=4,Neff=float(base.neff(p4)),elapsed_seconds=time.time()-start,
              mechanism_correlations=corr)
    (outdir/'nearzero_metadata.json').write_text(json.dumps(meta,indent=2))
    print(sm.to_string(index=False)); print(json.dumps(corr,indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--reps',type=int,default=100); ap.add_argument('--B',type=int,default=249); ap.add_argument('--seed',type=int,default=20260910)
    ap.add_argument('--phis',type=float,nargs='+',default=[.2,.5,.8]); ap.add_argument('--taus',type=float,nargs='+',default=TAUS_DEFAULT); ap.add_argument('--iters',type=int,default=30)
    ap.add_argument('--outdir',default='/mnt/data/operator_ablation_nearzero_B249')
    a=ap.parse_args(); run(a.phis,a.reps,a.B,a.seed,a.taus,a.outdir,a.iters)
