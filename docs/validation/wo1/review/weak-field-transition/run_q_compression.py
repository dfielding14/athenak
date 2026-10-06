from pathlib import Path
import json, subprocess, re
import numpy as np
root=Path(__file__).parent
s=(root/'input.athinput').read_text()
s=s.replace('test_mode = weak_field_transport','test_mode = q_affine_compression\nq_affine_direction = 0\nq_affine_rate = -0.2')
s=s.replace('x1min = 0.0','x1min = -50.0').replace('x1max = 1.0','x1max = 50.0')
s=s.replace('nlim = 50','nlim = 10000').replace('tlim = 100.0','tlim = 0.2')
inp=root/'affine.athinput'; inp.write_text(s)
res=[]
for direction,label in enumerate(['parallel','perpendicular','subfloor']):
 for nx in [128,256,512]:
  case=root/'affine'/f'{label}_{nx}';case.mkdir(parents=True,exist_ok=True)
  cmd=[str(root/'build/src/athena'),'-i',str(inp),f'problem/q_affine_direction={direction}',f'mesh/nx1={nx}',f'meshblock/nx1={nx}']
  run=subprocess.run(cmd,cwd=case,capture_output=True,text=True)
  (case/'stdout.log').write_text(run.stdout+run.stderr)
  (case/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
  if run.returncode: raise RuntimeError(run.stdout+run.stderr)
  wp=sorted((case/'tab').glob('*.mhd_w.*.tab'))[-1]
  w=np.loadtxt(wp);w=w[np.abs(w[:,2])<10]
  hist=np.loadtxt(next(case.glob('*.hst')));t=float(hist[-1,0]);rho_exact=1/(1-.2*t)
  q=np.log(w[:,8]/w[:,7]);factor=[-2,1,0][direction]
  q_expected=factor*np.log(rho_exact)
  r={'direction':label,'nx':nx,'time':t,'q_expected':q_expected,'q_mean':float(q.mean()),'q_min':float(q.min()),'q_max':float(q.max()),'max_q_error':float(np.max(np.abs(q-q_expected))),'max_q_density_error':float(np.max(np.abs(q-factor*np.log(w[:,3])))),'max_density_error':float(np.max(np.abs(w[:,3]-rho_exact))),'min_pressure':float(w[:,7:9].min()),'event_log':(case/'cgl_weak_field_transport.log').read_text()}
  assert r['max_q_error']<1e-4,r
  res.append(r)
(root/'affine-results.json').write_text(json.dumps(res,indent=2)+'\n')
for r in res: print(r)
