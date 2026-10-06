#!/usr/bin/env python3
"""Measure full-precision primitive/conserved SMR agreement and its order."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
from read_cgl_restart import read_restart

p=argparse.ArgumentParser()
p.add_argument('root',type=Path)
a=p.parse_args()
sys.path.insert(0,'/autofs/nccs-svm1_home2/dfielding/athenak-cgl-wo2-baseline/vis/python')
from athena_read import hst
records=[]
for directory in sorted(a.root.glob('n*-conserved'),key=lambda p:int(p.name.split('-')[0][1:])):
    n=int(directory.name.split('-')[0][1:]);other=a.root/f'n{n}-primitive'
    for phase,index in [('initial',0),('final',-1)]:
        x=read_restart(sorted((directory/'rst').glob('*.rst'))[index])
        y=read_restart(sorted((other/'rst').glob('*.rst'))[index])
        np.testing.assert_array_equal(x['loc'],y['loc'])
        assert x['time']==y['time'],(directory,x['time'],y['time'])
        ind=x['indcs'];ng,nx,ny,nz=ind[:4]
        active=np.zeros(next(iter(x['fields'].values())).shape[1:],dtype=bool)
        active[0:nz,ng:ng+ny,ng:ng+nx]=True
        halo=active.copy();halo[0:nz,ng-1:ng+ny+1,ng-1:ng+nx+1]=True
        ghost=halo & ~active
        fine=x['loc'][:,3]>x['root_level']
        weights=2.0**(-2*(x['loc'][:,3]-x['root_level']))/n**2
        record=dict(resolution=n,phase=phase,time=x['time'],cycle_conserved=x['cycle'],
                    cycle_primitive=y['cycle'],fields={})
        for field in x['fields']:
            delta=np.abs(x['fields'][field]-y['fields'][field])
            record['fields'][field]={
                'active_max':float(delta[:,active].max()),
                'active_volume_l1':float((delta[:,active].sum(axis=1)*weights).sum()),
                'fine_ghost_max':float(delta[fine][:,ghost].max()),
                'fine_ghost_mean':float(delta[fine][:,ghost].mean())}
        for mode,work in [('conserved',directory),('primitive',other)]:
            hist=hst(str(work/'smr_outflow.mhd.hst'))
            record[mode+'_energy_relative_drift']=float(np.max(np.abs(hist['tot-E']/hist['tot-E'][0]-1)))
            record[mode+'_bad_lf_counters']={key:float(hist[key][-1]) for key in
                ('lf_dfloor','lf_pfloor','lf_nonfin','lf_nonpos','lf_hardbd')}
        records.append(record)
orders=[]
for phase in ('initial','final'):
    rows=[r for r in records if r['phase']==phase]
    for x,y in zip(rows,rows[1:]):
        rate={'phase':phase,'coarse':x['resolution'],'fine':y['resolution'],'fields':{}}
        for field in x['fields']:
            rate['fields'][field]={}
            for metric,ex in x['fields'][field].items():
                ey=y['fields'][field][metric]
                rate['fields'][field][metric]=(float(np.log(ex/ey)/np.log(y['resolution']/x['resolution']))
                                               if ex>0 and ey>0 else None)
        orders.append(rate)
report={'measure':'primitive versus conserved prolongation at matched physical time; double restart state',
        'records':records,'orders':orders}
(a.root/'convergence.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(orders,indent=2))
