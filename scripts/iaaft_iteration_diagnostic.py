#!/usr/bin/env python3
import sys, math, json
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0,'/mnt/data/operator_ablation_src')
import operator_ablation_softpath as base

def run(outdir='/mnt/data/iaaft_iteration_diagnostic', reps=5, B=49, seed=20260931, iters=(10,30,100,300), phis=(.2,.5,.8)):
    outdir=Path(outdir); outdir.mkdir(parents=True,exist_ok=True)
    p4=base.p_profile(4); thr=base.norm.ppf(p4)
    seeds=np.random.SeedSequence(seed).spawn(len(phis)*reps); si=0; rows=[]
    for pi,phi in enumerate(phis):
        for rep in range(reps):
            x=base.generate_ar1(phi,base.N_MONTHS,np.random.default_rng(seeds[si])); si+=1
            mean=base.cubic_detrend(base.annual_mean(x))
            soft=base.cubic_detrend(base.annual_soft_count(x,thr,.01))
            hard=base.cubic_detrend(base.annual_count(x,thr))
            for ci,(name,obs) in enumerate([('mean_cubic',mean),('soft_tau0.01',soft),('hard_count',hard)]):
                for it in iters:
                    idx,err=base.iaaft_batch(obs,B,base.seed_rng(seed,777,pi,rep,ci),max_iter=it)
                    th=base.fs(idx)
                    rows.append(dict(phi=phi,rep=rep,condition=name,iters=it,
                                     median_spec_error=float(np.median(err)),q90_spec_error=float(np.quantile(err,.9)),
                                     median_theta=float(np.median(th)),sd_theta=float(np.std(th,ddof=1))))
    df=pd.DataFrame(rows); df.to_csv(outdir/'iaaft_iteration_diagnostic.csv',index=False)
    sm=df.groupby(['condition','iters']).agg(median_spec_error=('median_spec_error','median'),median_q90_spec_error=('q90_spec_error','median'),median_sd_theta=('sd_theta','median'),median_theta=('median_theta','median')).reset_index()
    sm.to_csv(outdir/'iaaft_iteration_summary.csv',index=False)
    print(sm.to_string(index=False))
if __name__=='__main__': run()
