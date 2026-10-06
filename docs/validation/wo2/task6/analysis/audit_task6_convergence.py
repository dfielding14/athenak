#!/usr/bin/env python3
"""Audit health and direct full-precision conservation at every probe resolution."""
import argparse
import json
import math
from pathlib import Path
import sys
import numpy as np
from read_cgl_restart import read_restart
sys.path.insert(0, '/autofs/nccs-svm1_home2/dfielding/athenak-cgl-wo2-baseline/vis/python')
from athena_read import hst
p=argparse.ArgumentParser();p.add_argument('root',type=Path);a=p.parse_args()
records=[]
for work in sorted(a.root.glob('n*-*')):
    if not work.is_dir(): continue
    n=int(work.name.split('-')[0][1:])
    states=[read_restart(path) for path in sorted((work/'rst').glob('*.rst'))]
    first,last=states[0],states[-1]
    ng,nx,ny,nz=first['indcs'][:4]
    sl=(slice(0,1),slice(ng,ng+ny),slice(ng,ng+nx))
    volume=2.0**(-2*(first['loc'][:,3]-first['root_level']))/n**2
    def total(state,name):
        values=state['fields'][name][(slice(None),)+sl]
        return math.fsum(map(float, values.reshape(len(volume),-1).sum(axis=1)*volume))
    conservation={name:total(last,name)-total(first,name) for name in ('energy','mu')}
    history=hst(str(work/'smr_outflow.mhd.hst'))
    user=hst(str(work/'smr_outflow.user.hst'))
    counters={name:float(np.max(np.abs(history[name]))) for name in
      ('lf_dfloor','lf_pfloor','lf_nonfin','lf_nonpos','lf_hardbd','lf_mirror','lf_firehs','lf_hwproj')}
    counters.update({name:float(np.max(np.abs(value))) for name,value in user.items()
                     if name.startswith('amr_') or name=='bad_state'})
    frozen={name:bool(np.array_equal(first['fields'][name][(slice(None),)+sl],
                                    last['fields'][name][(slice(None),)+sl]))
            for name in ('rho','scalar','B1','B2','B3')}
    positive={name:float(np.min(last['fields'][name][(slice(None),)+sl]))
              for name in ('rho','ppar','pperp')}
    record=dict(case=work.name,time=last['time'],cycles=last['cycle'],
                conservation_delta=conservation,repair_counters=counters,
                frozen=frozen,minima=positive,max_ndiv=float(np.max(user['max_ndiv'])))
    record['pass']=(all(abs(v)<5e-12 for v in conservation.values()) and
                    all(v==0 for v in counters.values()) and all(frozen.values()) and
                    all(v>0 for v in positive.values()) and record['max_ndiv']<1e-12 and
                    record['time']==.002)
    records.append(record)
report=dict(records=records,passed=all(r['pass'] for r in records),
            max_conservation_delta=max(abs(v) for r in records for v in r['conservation_delta'].values()))
(a.root/'health-conservation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
raise SystemExit(not report['passed'])
