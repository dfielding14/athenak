#!/usr/bin/env python3
"""Bounded signed pressure-band diagnostic on three retained snapshots.

This scratch script reads existing data only. It does not change the reusable
analyzer, physical input, acceptance criteria, or time-window definitions.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import sys

import numpy as np

p=argparse.ArgumentParser()
p.add_argument('--source',type=Path,required=True)
p.add_argument('--metrics',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
a.output.mkdir(parents=True,exist_ok=True)
script=a.source/'scripts/analyze_cgl_lf_physics_benchmark.py'
spec=importlib.util.spec_from_file_location('benchmark_analyzer',script)
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
d=json.loads(a.metrics.read_text())
def retained(path):
 return {'path':str(path.resolve()),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'bytes':path.stat().st_size}
assert retained(script)['sha256']==d['provenance']['analysis_script']['sha256']
results={'purpose':'Localize the observed global pressure-correlation reversal in wavenumber; three snapshots only, no campaign or new acceptance threshold.',
 'definitions':{
  'fields':'a=p_perp-<p_perp>, b=B^2/2-<B^2/2>; volume mean on uniform cells. Primitive float32 snapshots are read into float64 through the current benchmark reader.',
  'FFT':'A=fftn(a)/N, Bhat=fftn(b)/N on the full complex FFT grid; physical k_i=2*pi*n_i/L_i. Both members of each conjugate pair are included, exactly as Parseval requires.',
  'variance':'V_a(K)=sum_K |A|^2, V_b(K)=sum_K |Bhat|^2, in pressure^2. These are integrated band powers, not per-shell densities.',
  'covariance':'cov(K)=sum_K Re[A*conj(Bhat)], signed. Different Fourier bands are orthogonal; covariances and variances add across the disjoint bands.',
  'C':'cov(K)/sqrt(V_a(K)*V_b(K)); null if either band variance vanishes.',
  'R':'sum_K |A+Bhat|^2/[V_a(K)+V_b(K)]; null if denominator vanishes. Also verify equivalence to 1+2*cov/(V_a+V_b).',
  'forcing':'0<|k|<=3*pi (includes pure parallel modes). In this box the smallest nonzero |k| is pi, so this contains the full imposed forcing shell pi..3pi.',
  'resolved_subforcing':'3*pi<|k|<=min_i(pi*N_i/L_i)/4; terminology denotes resolved scales smaller than the forcing wavelengths, not k below forcing.',
  'high_k':'|k|>min_i(pi*N_i/L_i)/4, through every represented FFT mode; excludes DC and includes modes with high parallel wavenumber even at low k_perp.',
  'optional_union':'full_resolved combines the forcing and resolved_subforcing bands; all_nonzero combines all three. DC is separately reported as a roundoff check.',
  'bounds':'The squared physical-k comparisons use the same grids and Nyquist/4 convention as the reusable analyzer; no additional physical filter or taper.',
  'verification':'Full Fourier sums and direct real-space values must reproduce stored snapshot variances, signed covariance, C and R within 256*float64_epsilon*max(1,|reference|). This is an arithmetic/reader check, not a physical acceptance tolerance.',
  'limits':'Instantaneous three-snapshot diagnostic; not time-averaged cross spectra, not a causal forcing decomposition, not resolution convergence. Matching auto-spectra alone is not used to infer pressure balance.'},
 'provenance':{'script':retained(Path(__file__)),'analyzer':retained(script),'reader':retained(Path(module.paper.__file__)),'binary_reader':retained(Path(module.paper.bin_convert.__file__)),'metrics':retained(a.metrics),'simulation_revision':d['provenance']['simulation_revision'],'analysis_revision':d['provenance']['analysis_revision'],'command':sys.argv,'python':sys.version,'numpy':np.__version__,'host':platform.node(),'slurm_job_id':os.environ.get('SLURM_JOB_ID'),'slurm_step_id':os.environ.get('SLURM_STEP_ID')},'snapshots':[]}
for target in (10.,12.,14.):
 row=min(d['snapshots'],key=lambda r:abs(r['info']['time']-target))
 paths=row['info']['files'];path=Path(paths[0]['path'])
 fields,lengths,header,info=module.read_uniform(path,('p_perp','bcc1','bcc2','bcc3'))
 assert info['files']==paths, 'snapshot changed after retained analysis'
 perp=fields['p_perp'];mag=sum(fields[k]*fields[k] for k in ('bcc1','bcc2','bcc3'))/2
 aa=perp-np.mean(perp);bb=mag-np.mean(mag)
 fa=np.fft.fftn(aa)/aa.size;fb=np.fft.fftn(bb)/bb.size
 pa=abs(fa)**2;pb=abs(fb)**2;cross=(fa*np.conj(fb)).real;residual=abs(fa+fb)**2
 kx,ky,kz=module.grids(aa.shape,lengths);k2=kx*kx+ky*ky+kz*kz
 kf=3*np.pi;kr=min(np.pi*n/l for n,l in zip(aa.shape,lengths[::-1]))/4
 masks={'forcing':(k2>0)&(k2<=kf*kf),
        'resolved_subforcing':(k2>kf*kf)&(k2<=kr*kr),
        'high_k':k2>kr*kr,
        'full_resolved':(k2>0)&(k2<=kr*kr),
        'all_nonzero':k2>0,'DC':k2==0,'full_grid':np.ones(aa.shape,dtype=bool)}
 def band(mask):
  va=float(pa[mask].sum());vb=float(pb[mask].sum());cov=float(cross[mask].sum());vr=float(residual[mask].sum())
  return {'modes':int(mask.sum()),'variance_perp':va,'variance_magnetic':vb,'covariance':cov,
          'sum_pressure_variance':vr,'correlation':float(cov/np.sqrt(va*vb)) if va*vb>0 else None,
          'normalized_residual_variance':vr/(va+vb) if va+vb>0 else None,
          'R_from_covariance':1+2*cov/(va+vb) if va+vb>0 else None}
 bands={name:band(mask) for name,mask in masks.items()}
 stored={key.removeprefix('pressure_'):val for key,val in row['scalars'].items() if key.startswith('pressure_')}
 real=module.pressure_balance(perp,mag)
 checks={}
 for key,val in stored.items():
  fft=bands['full_grid'][key];direct=real[key]
  allowed=256*np.finfo(float).eps*max(1.,abs(val))
  assert abs(fft-val)<=allowed,(target,key,'FFT',fft,val)
  assert abs(direct-val)<=allowed,(target,key,'direct',direct,val)
  checks[key]={'stored':val,'direct_real_space':direct,'full_fft':fft,'fft_minus_stored':fft-val,'allowed_arithmetic_roundoff':allowed}
 for key in ('variance_perp','variance_magnetic','covariance','sum_pressure_variance'):
  partition=sum(bands[name][key] for name in ('forcing','resolved_subforcing','high_k'))
  assert abs(partition-bands['all_nonzero'][key])<=256*np.finfo(float).eps*max(1.,abs(partition))
 assert bands['forcing']['modes']==38, 'canonical physical forcing shell count changed'
 snapshot={'target':target,'time':info['time'],'cycle':info['cycle'],'info':info,
  'forcing_kmax':kf,'resolved_kmax':kr,'bands':bands,'verification':checks}
 results['snapshots'].append(snapshot)
 print(json.dumps({'time':info['time'],'bands':{name:{key:val for key,val in value.items() if key in ('covariance','correlation','normalized_residual_variance','variance_perp','variance_magnetic')} for name,value in bands.items() if name in ('forcing','resolved_subforcing','high_k','all_nonzero')},'verification':'passed'}),flush=True)
(a.output/'results.json').write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')
