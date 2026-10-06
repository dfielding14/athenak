#!/usr/bin/env python3
"""Check CGL A/mu consistency from independent thermodynamic identities."""
import argparse
import json
from pathlib import Path
import numpy as np

p = argparse.ArgumentParser()
p.add_argument('trace', type=Path)
p.add_argument('--output', type=Path, required=True)
p.add_argument('--tolerance', type=float, default=2e-12)
a = p.parse_args()
points = {'Init_after_c2p', 'BeginCGLLandauFluidSTSSweep_begin',
          'BeginCGLLandauFluidSTSSweep_end',
          'CGLLandauFluidPrimitiveRefresh_end', 'EndCGLLandauFluidSTSSweep_end'}
report = {'tolerance': a.tolerance, 'snapshots': [], 'pass': True}
for file in sorted(a.trace.glob('*.json')):
    d = json.loads(file.read_text())
    if d['point'] not in points:
        continue
    def array(name):
        meta = d['arrays'][name]
        return np.fromfile(file.parent/meta['file'],dtype='=f8').reshape(meta['shape'])
    bounds = d['active_kji']
    halo = tuple(slice(lo-1 if hi>lo else lo,hi+2 if hi>lo else hi+1)
                 for lo,hi in bounds)
    idx = (slice(None),slice(None))+halo
    u,w,b = (array(name)[idx] for name in ('u','w','bcc'))
    b1,b2,b3 = (array(name) for name in ('b1','b2','b3'))
    face_b = np.stack((.5*(b1[:,:,:,:-1]+b1[:,:,:,1:]),
                       .5*(b2[:,:,:-1,:]+b2[:,:,1:,:]),
                       .5*(b3[:,:-1,:,:]+b3[:,1:,:,:])),axis=1)[idx]
    rho = w[:,0]
    ppar,pperp = w[:,4],w[:,5]
    bsqr = (b*b).sum(axis=1)
    bmag = np.sqrt(bsqr)
    finite = all(np.isfinite(x).all() for x in (u,w,b,face_b))
    admissible = finite and bool((rho>0).all() and (ppar>0).all() and (pperp>0).all())
    expected_energy = pperp+.5*ppar+.5*rho*(w[:,1:4]**2).sum(axis=1)+.5*bsqr
    with np.errstate(divide='ignore',invalid='ignore'):
        expected_slot = (pperp/bmag if d['representation']=='mu' else
            rho*(np.log(pperp)-np.log(ppar)+2*np.log(rho)-3*np.log(bmag)))
    comparisons = {'density':(u[:,0],rho), 'momentum':(u[:,1:4],rho[:,None]*w[:,1:4]),
                   'energy':(u[:,4],expected_energy), 'slot':(u[:,5],expected_slot),
                   'face_to_cell_B':(b,face_b)}
    if u.shape[1]>6:
        comparisons['scalars']=(u[:,6:],rho[:,None]*w[:,6:])
    residual = {name:float(np.max(np.abs(actual-exact)/(1+np.abs(exact))))
                for name,(actual,exact) in comparisons.items()}
    passed = admissible and all(np.isfinite(v) and v<a.tolerance for v in residual.values())
    report['pass'] &= passed
    report['snapshots'].append({k:d[k] for k in ('rank','point','sweep','stage','representation')} |
        {'file':file.name,'admissible':admissible,'normalized_residual':residual,'pass':passed})
report['snapshots_checked'] = len(report['snapshots'])
report['pass'] &= bool(report['snapshots'])
report['max_normalized_residual'] = max((max(x['normalized_residual'].values())
                                        for x in report['snapshots']),default=None)
a.output.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='snapshots'},indent=2))
raise SystemExit(not report['pass'])
