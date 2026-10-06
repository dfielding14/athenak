from pathlib import Path
import json, subprocess, re, hashlib
import numpy as np
root=Path(__file__).parent
inp=root/'input.athinput'
s=inp.read_text()
if '<output4>' not in s:
 inp.write_text(s+'\n<output4>\nfile_type = hst\ndata_format = %24.17e\ndcycle = 1\n')
exe=root/'build/src/athena'
results=[]
for recon in ['dc','plm']:
 for vel in [10.,-10.,1.,0.]:
  for floor in [1e-14,1e-10,1e-6]:
   case=root/'runs'/f'{recon}_v{vel:g}_bf{floor:g}'; case.mkdir(parents=True,exist_ok=True)
   cmd=[str(exe),'-i',str(inp),f'mhd/reconstruct={recon}',f'problem/weak_field_velocity={vel}',f'mhd/bfloor={floor}']
   run=subprocess.run(cmd,cwd=case,capture_output=True,text=True)
   (case/'stdout.log').write_text(run.stdout+run.stderr)
   (case/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
   r={'case':case.name,'reconstruct':recon,'velocity':vel,'bfloor':floor,'returncode':run.returncode}
   if run.returncode:
    results.append(r);continue
   wfiles=sorted((case/'tab').glob('*.mhd_w.*.tab'))
   r['cycles']=[]
   for wp in wfiles:
    w=np.loadtxt(wp); assert np.isfinite(w).all() and (w[:,3]>0).all(); ratio=w[:,8]/w[:,7]
    cycle=int(re.search(r'cycle=(\d+)',wp.read_text().splitlines()[0])[1])
    r['cycles'].append({'cycle':cycle,'min_ratio':float(ratio.min()),'max_ratio':float(ratio.max()),'min_ppar':float(w[:,7].min()),'min_pperp':float(w[:,8].min())})
   r['min_ratio']=min(x['min_ratio'] for x in r['cycles']);r['max_ratio']=max(x['max_ratio'] for x in r['cycles'])
   r['positive']=all(x['min_ppar']>1e-10 and x['min_pperp']>1e-10 for x in r['cycles'])
   r['bounded']=r['positive'] and r['min_ratio']>=.5 and r['max_ratio']<=2
   upaths=sorted((case/'tab').glob('*.mhd_u_bcc.*.tab'))
   u0=np.loadtxt(upaths[0]); uf=np.loadtxt(upaths[-1]); wf=np.loadtxt(wfiles[-1]); w0=np.loadtxt(wfiles[0])
   hpath=list(case.glob('*.hst'))[0]; hist=np.loadtxt(hpath)
   t=float(hist[-1,0]); fluxdiff=.25*abs(vel)
   r['time']=t; r['energy_initial']=float(u0[:,7].mean());r['energy_final']=float(uf[:,7].mean())
   r['energy_expected']=r['energy_initial']+fluxdiff*t
   r['energy_residual']=r['energy_final']-r['energy_expected']
   r['max_boundary_w_change']=float(np.max(np.abs(wf[[0,-1],3:]-w0[[0,-1],3:])))
   r['event_log']=(case/'cgl_weak_field_transport.log').read_text()
   r['physical_digest']=hashlib.sha256(b''.join(p.read_bytes() for p in wfiles+upaths)).hexdigest()
   results.append(r)
(root/'results.json').write_text(json.dumps(results,indent=2)+'\n')
for r in results:
 print(r['case'],r['returncode'], 'range',r.get('min_ratio'),r.get('max_ratio'),'energy residual',r.get('energy_residual'),'boundary change',r.get('max_boundary_w_change'),'bounded',r.get('bounded'))

assert len(results) == 24
assert all(r['bounded'] and r['returncode'] == 0 for r in results)
assert all({c['cycle'] for c in r['cycles']} == set(range(51)) for r in results)
assert max(abs(r['energy_residual']) for r in results) < 1e-13
assert all(r['max_boundary_w_change'] == 0 for r in results)
assert all(not [l for l in r['event_log'].splitlines() if l.strip() and not l.startswith('#')] for r in results)
