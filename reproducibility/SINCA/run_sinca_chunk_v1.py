#!/usr/bin/env python3
import argparse, hashlib, time
from pathlib import Path
import pandas as pd, numpy as np
import importlib.util
base_path='/mnt/data/SINCA_PM25_AOAS_INFERENCE_V1/run_sinca_null_transport_v1_csv.py'
spec=importlib.util.spec_from_file_location('sinca_v1',base_path)
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

ap=argparse.ArgumentParser()
ap.add_argument('--data',required=True)
ap.add_argument('--cohort',required=True)
ap.add_argument('--out',required=True)
ap.add_argument('--start',type=int,required=True)
ap.add_argument('--end',type=int,required=True)
ap.add_argument('--B',type=int,default=499)
ap.add_argument('--seed',type=int,default=20260911)
a=ap.parse_args()

df,cohort=m.load_data(a.data,a.cohort)
ids=sorted(cohort.station_id.astype(str).tolist())
sha=hashlib.sha256(('\n'.join(ids)+'\n').encode()).hexdigest()
if sha != m.EXPECTED_COHORT_SHA: raise RuntimeError('frozen cohort hash mismatch')
ordered=cohort.sort_values('station_id').reset_index(drop=True)
sub=ordered.iloc[a.start:a.end]
recs=[]; t=time.time()
for j,row in sub.iterrows():
    dates,x=m.prepare_station(df,row.short_code)
    for branch in ['fixed','equalized']:
        r=m.run_branch(row.short_code,row.station_id,x,dates,a.B,a.seed,branch)
        r.update(station_id=row.station_id,short_code=row.short_code,station_name=row.station_name,
                 region=row.region,neff_frozen=float(row.neff_p_monthly_gt50))
        recs.append(r)
    print(f'{row.station_id} done {time.time()-t:.1f}s', flush=True)
pd.DataFrame(recs).to_csv(a.out,index=False)
print('wrote',a.out,'rows',len(recs))
