#!/usr/bin/env python3
import sys, numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0,'/mnt/data/operator_ablation_src')
import operator_ablation_softpath as base

def geom(arr):
    if arr.ndim==1: arr=arr[None,:]
    out=[]
    for z in arr:
        u=np.quantile(z,.90); ii=(z>u).astype(int); pos=np.flatnonzero(ii)
        adj=int(np.sum(ii[:-1]*ii[1:])); gaps=np.diff(pos) if len(pos)>1 else np.array([np.nan])
        clusters=int(ii[0]+np.sum((ii[1:]==1)&(ii[:-1]==0)))
        out.append((adj,clusters,float(np.nanmean(gaps)),float(np.nanstd(gaps,ddof=1)) if len(gaps)>2 else np.nan))
    return np.array(out,float)

def run(out='/mnt/data/exceedance_geometry_diagnostic',reps=20,B=99,seed=20260950):
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    p4=base.p_profile(4); thr=base.norm.ppf(p4); rows=[]
    ss=np.random.SeedSequence(seed).spawn(3*reps); si=0
    for pi,phi in enumerate([.2,.5,.8]):
      for rep in range(reps):
        x=base.generate_ar1(phi,base.N_MONTHS,np.random.default_rng(ss[si]));si+=1
        nx,_=base.native_surrogate_batch(x,B,base.seed_rng(seed,111,pi,rep))
        ox=base.generate_ar1(phi,base.N_MONTHS,base.seed_rng(seed,222,pi,rep),B=B)
        conds=[]
        for tau in [.1,.01,.001]:
            obs=base.cubic_detrend(base.annual_soft_count(x,thr,tau)); nat=base.cubic_detrend(base.annual_soft_count(nx,thr,tau)); ora=base.cubic_detrend(base.annual_soft_count(ox,thr,tau)); conds.append((tau,obs,nat,ora))
        obs=base.cubic_detrend(base.annual_count(x,thr)); nat=base.cubic_detrend(base.annual_count(nx,thr)); ora=base.cubic_detrend(base.annual_count(ox,thr)); conds.append((0.,obs,nat,ora))
        for ci,(tau,obs,nat,ora) in enumerate(conds):
            idx,_=base.iaaft_batch(obs,B,base.seed_rng(seed,777,pi,rep,ci),max_iter=30)
            for name,arr in [('oracle',ora),('native',nat),('index',idx)]:
                gg=geom(arr); th=base.fs(arr)
                rows.append(dict(phi=phi,rep=rep,tau=tau,null=name,
                    mean_adj=float(np.mean(gg[:,0])),sd_adj=float(np.std(gg[:,0],ddof=1)),
                    mean_clusters=float(np.mean(gg[:,1])),sd_clusters=float(np.std(gg[:,1],ddof=1)),
                    mean_gap_sd=float(np.nanmean(gg[:,3])),sd_theta=float(np.std(th,ddof=1)),mean_theta=float(np.mean(th))))
    df=pd.DataFrame(rows); df.to_csv(out/'exceedance_geometry_trajectory.csv',index=False)
    sm=df.groupby(['tau','null']).median(numeric_only=True).reset_index()[['tau','null','mean_adj','sd_adj','mean_clusters','sd_clusters','mean_gap_sd','sd_theta','mean_theta']]
    sm.to_csv(out/'exceedance_geometry_summary.csv',index=False); print(sm.to_string(index=False))
if __name__=='__main__': run()
