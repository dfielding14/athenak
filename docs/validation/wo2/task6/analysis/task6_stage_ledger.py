#!/usr/bin/env python3
"""Locate integral changes and active magnetic changes inside the frozen LF sweep."""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from read_cgl_restart import read_restart

p=argparse.ArgumentParser();p.add_argument('run',type=Path);a=p.parse_args()
initial=read_restart(sorted((a.run/'rst').glob('*.rst'))[0])
dx=1.0/32/2.0**(initial['loc'][:,3]-initial['root_level'])
records=[];previous=None;previous_b={}
for file in sorted((a.run/'trace').glob('*.json')):
    d=json.loads(file.read_text());gids=np.array(d['gids']);volume=dx[gids]**2
    def get(name):
        x=d['arrays'][name]
        return np.fromfile(file.parent/x['file'],dtype='=f8').reshape(x['shape'])
    bounds=d['active_kji'];sl=tuple(slice(lo,hi+1) for lo,hi in bounds)
    u=get('u')[(slice(None),slice(None))+sl]
    b=get('bcc')[(slice(None),slice(None))+sl]
    def integral(x):return math.fsum(float(v) for v in (x.reshape(len(gids),-1).sum(axis=1)*volume))
    if d['representation']=='mu':moment=u[:,5]
    else:
        rho=u[:,0];bsqr=(b*b).sum(axis=1);bm=np.sqrt(bsqr)
        internal=u[:,4]-.5*(u[:,1:4]**2).sum(axis=1)/rho-.5*bsqr
        ratio=np.exp(u[:,5]/rho-2*np.log(rho)+3*np.log(bm))
        moment=internal*ratio/(.5+ratio)/bm
    entry={k:d[k] for k in ('point','sweep','stage','representation')}
    entry.update(file=file.name,energy=integral(u[:,4]),mu=integral(moment),
                 active_face_changes={})
    for axis in (1,2,3):
        fb=[list(x) for x in bounds];fb[3-axis][1]+=1
        fs=tuple(slice(lo,hi+1) for lo,hi in fb)
        field=get('b'+str(axis))[(slice(None),)+fs]
        if axis in previous_b:
            delta=np.abs(field-previous_b[axis])
            entry['active_face_changes'][str(axis)]=dict(max_abs=float(delta.max()),values=int((delta!=0).sum()))
        previous_b[axis]=field
    if d['point'] in ('STSFluxes_end','SendFlux_end','RecvFlux_end'):
        entry['global_flux_divergence']={}
        for variable,label in [(4,'energy'),(5,'mu')]:
            boundary=np.zeros(len(gids))
            for axis in (1,2):
                fb=[list(x) for x in bounds];fb[3-axis][1]+=1
                fs=tuple(slice(lo,hi+1) for lo,hi in fb)
                f=get('flux'+str(axis))[(slice(None),variable)+fs]
                difference=np.take(f,-1,axis=4-axis)-np.take(f,0,axis=4-axis)
                boundary+=difference.reshape(len(gids),-1).sum(axis=1)*dx[gids]
            entry['global_flux_divergence'][label]=math.fsum(map(float,boundary))
    if previous:
        entry['energy_change']=entry['energy']-previous['energy']
        entry['mu_change']=entry['mu']-previous['mu']
    records.append(entry);previous=entry
(a.run/'stage-ledger.json').write_text(json.dumps(records,indent=2)+'\n')
for r in records:
    if abs(r.get('energy_change',0))>2e-14 or 'global_flux_divergence' in r:
        print(r['file'],r.get('energy_change'),r.get('mu_change'),r.get('global_flux_divergence'))
