from pathlib import Path
import json,re,math,csv
base=Path(__file__).resolve().parent
summary=[]
for run in sorted((base/'runs').iterdir()):
  cells=[]; states=[]
  for path in sorted((run/'tab').glob('*.mhd_w_bcc.*.tab')):
    lines=path.read_text().splitlines(); head=lines[0]
    time=float(re.search(r'time=(\S+)',head).group(1)); cycle=int(re.search(r'cycle=(\d+)',head).group(1))
    data=[[float(v) for v in ln.split()] for ln in lines if ln and not ln.startswith('#')]
    consfile=Path(str(path).replace('.mhd_w_bcc.','.mhd_u_bcc.'))
    cons=[[float(v) for v in ln.split()] for ln in consfile.read_text().splitlines() if ln and not ln.startswith('#')]
    ratios=[r[8]/r[7] if r[7]!=0 else math.inf for r in data]
    finite=all(math.isfinite(x) for r in data for x in r)
    pressures=[r[7] for r in data]+[r[8] for r in data]
    item=dict(cycle=cycle,time=time,min_ratio=min(ratios),max_ratio=max(ratios),min_pressure=min(pressures),finite=finite)
    states.append(item)
    for w,u,ratio in zip(data,cons,ratios):
      cells.append(dict(cycle=cycle,time=time,i=int(w[1])-3,x=w[2],rho=w[3],vx=w[4],ppar=w[7],pperp=w[8],ratio=ratio,Bx=w[9],By=w[10],Bz=w[11],A=u[8],E=u[7]))
  firstfail=next((s for s in states if not s['finite'] or s['min_ratio']<.5 or s['max_ratio']>2),None)
  unique={s['cycle']:s for s in states}
  final=states[-1] if states else {}
  log=(run/'run.log').read_text()
  item=dict(name=run.name,working_directory=str(run),command=[str(base/'build/src/athena'),'-i','input.athinput'],ncycle=final.get('cycle'),time=final.get('time'),min_ratio=final.get('min_ratio'),max_ratio=final.get('max_ratio'),min_pressure=final.get('min_pressure'),first_fail_cycle=firstfail['cycle'] if firstfail else None,first_fail_time=firstfail['time'] if firstfail else None,first_fail_max_ratio=firstfail['max_ratio'] if firstfail else None,status='completed' if 'Terminating' in log else 'inspect_log',min_pressure_all_outputs=min((s['min_pressure'] for s in states),default=None),first_five=[s for c,s in sorted(unique.items()) if c<=5],max_ratio_all_outputs=max((s['max_ratio'] for s in states),default=None))
  summary.append(item)
  (run/'summary.json').write_text(json.dumps(item,indent=2)+'\n')
  (run/'cycles.json').write_text(json.dumps(list(unique.values()),indent=2)+'\n')
  if cells:
    with (run/'cells.csv').open('w') as f:
      writer=csv.DictWriter(f,fieldnames=list(cells[0]));writer.writeheader();writer.writerows(cells)
(base/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
for x in summary:
  print(x['name'],x['ncycle'],x['time'],x['min_ratio'],x['max_ratio'],x['min_pressure'],'firstfail',x['first_fail_cycle'],x['first_fail_max_ratio'])
