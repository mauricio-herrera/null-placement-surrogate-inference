#!/usr/bin/env python3
import sys, math, time
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,'/mnt/data/SINCA_PM25_AOAS_RESULTS_V1_FINAL')
import run_sinca_null_transport_v1_csv as r

DATA='/mnt/data/sinca_csv_full/sinca_daily_station.csv'
COHORT='/mnt/data/SINCA_PM25_AOAS_FROZEN_COHORT_V1/FROZEN_COHORT_V1.csv'
B=499; SEED=20260911

def weekly_exposure_surr(surr, obs, threshold_day):
    nobs=obs.reshape(r.N_WEEKS,7).sum(axis=1)
    if np.any(nobs==0):
        return None
    evt=np.zeros((surr.shape[0],r.N_DAYS),float)
    evt[:,obs]=(surr[:,obs]>threshold_day[obs][None,:]).astype(float)
    e=evt.reshape(surr.shape[0],r.N_WEEKS,7).sum(axis=2)
    return 7*e/nobs[None,:]

import argparse
ap=argparse.ArgumentParser(); ap.add_argument('--start',type=int,default=0); ap.add_argument('--end',type=int,default=31); ap.add_argument('--out',required=True); args=ap.parse_args()
df, cohort = r.load_data(DATA,COHORT)
cohort=cohort.sort_values('station_id').reset_index(drop=True).iloc[args.start:args.end]
rows=[]; t0=time.time()
for j,row in cohort.iterrows():
    dates,x=r.prepare_station(df,row.short_code)
    obs=np.isfinite(x)
    th=np.full(r.N_DAYS,50.0)
    wobs=r.weekly_exposure_standardized(x,dates,th)
    rec={'station_id':row.station_id,'short_code':row.short_code,'station_name':row.station_name,'region':row.region,
         'defined':wobs is not None}
    if wobs is None:
        rows.append(rec); continue
    z=r.residualize_weekly(wobs)
    nat,errn=r.native_surrogate_missing(x,dates,B,r.seed_rng(SEED,row.station_id,'fixed','native'),12)
    r.check_native_event_count_conservation(x,nat,dates,th)
    wnat=weekly_exposure_surr(nat,obs,th)
    znat=r.residualize_weekly(wnat)
    idx,erri=r.iaaft_batch(z,B,r.seed_rng(SEED,row.station_id,'fixed','index_exposure'),30)
    s=r.summarize_stat(z,znat,idx,lambda a:r.fs(a,.90),'fs90_exposure')
    rec.update(s)
    rn=rec['rej_native_fs90_exposure']; ri=rec['rej_index_fs90_exposure']
    rec['paired_outcome']='both' if rn and ri else 'native_only' if rn else 'index_only' if ri else 'neither'
    rows.append(rec)
    print(j+1,row.station_id,rec.get('p_native_fs90_exposure'),rec.get('p_index_fs90_exposure'),flush=True)
out=pd.DataFrame(rows)
out.to_csv(args.out,index=False)
g=out[out.defined==True]
print('defined',len(g),'native rej',int(g.rej_native_fs90_exposure.sum()),'index rej',int(g.rej_index_fs90_exposure.sum()))
print(g.paired_outcome.value_counts())
print('elapsed',time.time()-t0)
