#!/usr/bin/env python3
import sys, numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0,'/mnt/data/operator_ablation_src')
import operator_ablation_softpath as base

def phase_metrics(idx):
    z=idx-idx.mean(axis=1,keepdims=True)
    F=np.fft.rfft(z,axis=1)
    ph=np.angle(F[:,1:])
    R=np.abs(np.mean(np.exp(1j*ph),axis=0))
    return float(np.median(R)), float(np.quantile(R,.9)), float(np.mean(R))

def run(out='/mnt/data/phase_rigidity_diagnostic', reps=10,B=149,seed=20260940):
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    p4=base.p_profile(4); thr=base.norm.ppf(p4); rows=[]
    ss=np.random.SeedSequence(seed).spawn(3*reps); si=0
    for pi,phi in enumerate([.2,.5,.8]):
      for rep in range(reps):
        x=base.generate_ar1(phi,base.N_MONTHS,np.random.default_rng(ss[si])); si+=1
        conds=[('mean',base.cubic_detrend(base.annual_mean(x)))]
        for tau in [.1,.025,.01,.001]:
            conds.append((f'tau{tau}',base.cubic_detrend(base.annual_soft_count(x,thr,tau))))
        conds.append(('hard',base.cubic_detrend(base.annual_count(x,thr))))
        for ci,(name,obs) in enumerate(conds):
            idx,err=base.iaaft_batch(obs,B,base.seed_rng(seed,777,pi,rep,ci),max_iter=30)
            medR,q90R,meanR=phase_metrics(idx)
            theta=base.fs(idx)
            # diversity of adjacent-difference statistic
            ntv=base.ntv(idx)
            rows.append(dict(phi=phi,rep=rep,condition=name,median_phase_resultant=medR,q90_phase_resultant=q90R,mean_phase_resultant=meanR,
                             median_spec_error=float(np.median(err)),sd_theta=float(np.std(theta,ddof=1)),sd_ntv=float(np.std(ntv,ddof=1))))
    df=pd.DataFrame(rows); df.to_csv(out/'phase_rigidity_trajectory.csv',index=False)
    sm=df.groupby('condition').median(numeric_only=True).reset_index()[['condition','median_phase_resultant','q90_phase_resultant','mean_phase_resultant','median_spec_error','sd_theta','sd_ntv']]
    order=['mean','tau0.1','tau0.025','tau0.01','tau0.001','hard']; sm['ord']=sm.condition.map({x:i for i,x in enumerate(order)}); sm=sm.sort_values('ord').drop(columns='ord')
    sm.to_csv(out/'phase_rigidity_summary.csv',index=False); print(sm.to_string(index=False))
if __name__=='__main__': run()
