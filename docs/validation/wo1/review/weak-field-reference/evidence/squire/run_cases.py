from pathlib import Path
import subprocess, csv, json, math, re
root=Path(__file__).resolve().parent
common='''<job>
problem_id=b4_reference
<time>
cfl_number=0.4
integrator=rk2
xorder=2
ncycle_out=1
nlim={nlim}
tlim={tlim}
<mesh>
nx1={nx}
x1min=0.0
x1max=1.0
ix1_bc=outflow
ox1_bc=outflow
nx2=1
x2min=0.0
x2max=1.0
ix2_bc=outflow
ox2_bc=outflow
nx3=1
x3min=0.0
x3max=1.0
ix3_bc=outflow
ox3_bc=outflow
<meshblock>
nx1={nx}
nx2=1
nx3=1
<hydro>
gamma=1.6666666666666667
dfloor=1.0e-12
pfloor=1.0e-12
bsqr_floor={floor}
nu_coll=0.0
firehose_limiter=false
mirror_limiter=false
<problem>
velocity={velocity}
profile_width={width}
'''

def run(name,velocity,floor,nlim=50,tlim=100.,nx=128,width=0.):
    out=root/'runs'/name
    out.mkdir(parents=True,exist_ok=True)
    (out/'input.athinput').write_text(common.format(velocity=velocity,floor=floor,nlim=nlim,tlim=tlim,nx=nx,width=width))
    cmd=[str(root/'bin/athena'),'-i','input.athinput']
    with (out/'run.log').open('w') as f:
        try: status=subprocess.run(cmd,cwd=out,stdout=f,stderr=subprocess.STDOUT,timeout=40).returncode
        except subprocess.TimeoutExpired: status='timeout'
    samples=[]
    for path in sorted(out.glob('state_*_cycle_0.csv')):
        rows=list(csv.DictReader(path.open()))
        ratios=[float(r['pperp'])/float(r['ppar']) for r in rows]
        finite=all(math.isfinite(float(v)) for r in rows for v in r.values())
        sample={'cycle':int(rows[0]['cycle']),'time':float(rows[0]['time']), 'ratio_min':min(ratios), 'ratio_max':max(ratios), 'min_pressure':min(float(r[p]) for r in rows for p in ('pperp','ppar')), 'finite':finite}
        samples.append(sample)
    fails=[s for s in samples if not s['finite'] or s['ratio_min']<.5 or s['ratio_max']>2]
    result={'case':name,'command':cmd,'cwd':str(out),'returncode':status,'velocity':velocity,'bsqr_floor':floor,'nx':nx,'profile_width':width,'requested_nlim':nlim,'requested_tlim':tlim,'completed_cycles':len(samples),'first_failure':fails[0] if fails else None,'final':samples[-1] if samples else None,'overall_ratio_min':min((s['ratio_min'] for s in samples),default=None),'overall_ratio_max':max((s['ratio_max'] for s in samples),default=None),'overall_min_pressure':min((s['min_pressure'] for s in samples),default=None)}
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    (out/'cycle_summary.json').write_text(json.dumps(samples,indent=2)+'\n')
    print(json.dumps(result),flush=True)
    return result

if __name__=='__main__':
    results=[]
    for label,floor in [('facefloor',2e-20),('magnetized',2e-28)]:
        for v in (10.,-10.):
            for duration,nlim,tlim in [('50cycles',50,100.),('t012',10000,.012)]:
                results.append(run(f'{label}_{"pos" if v>0 else "neg"}_{duration}',v,floor,nlim,tlim))
    if any(r['first_failure'] for r in results):
        for nx in (128,256,512):
            results.append(run(f'smooth_n{nx}_t012',10.,2e-28,10000,.012,nx,.04))
    (root/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
