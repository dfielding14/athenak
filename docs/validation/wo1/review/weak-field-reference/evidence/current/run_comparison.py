from pathlib import Path
import csv, hashlib, json, os, subprocess, re
R=Path(__file__).parent
BIN=R/'athena-instrumented-4way'
INPUT=R/'comparison.athinput'

def summarize(path):
    groups={}; stages={}
    with (path/'b4-state.csv').open() as f:
        for raw in csv.DictReader(f):
            phase=raw['phase']; row={k:(v if k=='phase' else float(v)) for k,v in raw.items()}
            c=int(row['cycle']); cell=int(row['cell']); stage=int(row['stage'])
            if phase=='coll_post': groups.setdefault(c+1,[]).append(row)
            if phase in ('coll_pre','coll_post','c2p_pre','c2p_post','flux_pre_fofc','flux_post_fofc'):
                stages[(c,stage,phase,cell)]=row
    cycles=[]
    for c,rows in sorted(groups.items()):
        ratios=[x['pperp']/x['ppar'] for x in rows]
        cycle=dict(cycle=c,time=rows[0]['time']+rows[0]['dt'],dt=rows[0]['dt'],
            min_ratio=min(ratios),max_ratio=max(ratios),min_pressure=min(min(x['ppar'],x['pperp']) for x in rows),
            sum_rho=sum(x['rho'] for x in rows),sum_E=sum(x['E'] for x in rows),sum_A=sum(x['A'] for x in rows),
            max_ratio_cell=int(rows[ratios.index(max(ratios))]['cell']))
        cycles.append(cycle)
    assert cycles, path
    changes={}
    for pre,post,cols in [('coll_pre','coll_post',['A','E','ppar','pperp']),('c2p_pre','c2p_post',['A','E']),('flux_pre_fofc','flux_post_fofc',['FrhoL','FAL','FEL','FByL'])]:
        changed=[]
        for (c,s,ph,i),a in stages.items():
            if ph!=pre or (c,s,post,i) not in stages: continue
            b=stages[(c,s,post,i)]
            ds={x:b[x]-a[x] for x in cols if b[x]!=a[x]}
            if ds: changed.append(dict(cycle=c+1,stage=s,cell=i,deltas=ds))
        changes[pre+'_'+post]=changed
    events=[]
    for line in (path/'cgl_weak_field_transport.log').read_text().splitlines():
        if line and not line.startswith('#'): events.append([int(float(x)) for x in line.split()])
    first=next((x for x in cycles if x['min_ratio']<0.5 or x['max_ratio']>2),None)
    summary=dict(ncycle=cycles[-1]['cycle'],time=cycles[-1]['time'],
        min_ratio=min(x['min_ratio'] for x in cycles),max_ratio=max(x['max_ratio'] for x in cycles),
        min_pressure=min(x['min_pressure'] for x in cycles),first_failure_ratio_outside_half_to_two=first,
        final=cycles[-1],cycles=cycles,first_five_changes=changes,event_counter_rows=events,
        event_counter_columns=['cycle','eos_dfloor','eos_efloor','eos_tfloor','eos_vceil','eos_fail','c2p_it','fofc'])
    (path/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary

def run(name,face='wo1',wall=True,v=10,floor='1e-10',stop='cycles',nx=128,width=0):
    path=R/'runs'/name;path.mkdir(parents=True,exist_ok=True)
    args=['output1/variable=mhd_w_bcc','time/ndiag=1','mhd/nu_coll=0','mhd/mirror_limiter=false',
          'mhd/firehose_limiter=false','mhd/backup_limiters=false',f'problem/weak_field_velocity={v}',
          f'mhd/bfloor={floor}',f'mesh/nx1={nx}',f'meshblock/nx1={nx}',f'problem/profile_width={width}']
    args+=['time/nlim=50','time/tlim=100'] if stop=='cycles' else ['time/nlim=10000','time/tlim=0.012']
    env=dict(os.environ)
    for k in ['B4_TRACE','B4_NO_WALL','B4_ORIGINAL_FACE']:env.pop(k,None)
    env['B4_TRACE']='1'
    if not wall:env['B4_NO_WALL']='1'
    if face=='original':env['B4_ORIGINAL_FACE']='1'
    cmd=[str(BIN),'-i',str(INPUT)]+args
    meta=dict(face=face,wall=wall,velocity=v,bfloor=floor,stop=stop,nx=nx,width=width,command=cmd,
              env={k:v for k,v in env.items() if k.startswith('B4_')},binary_sha256=hashlib.sha256(BIN.read_bytes()).hexdigest())
    (path/'command.json').write_text(json.dumps(meta,indent=2)+'\n')
    with (path/'run.log').open('w') as log:
        result=subprocess.run(cmd,cwd=path,env=env,stdout=log,stderr=subprocess.STDOUT)
    meta['exitcode']=result.returncode
    if result.returncode: raise RuntimeError(name+' failed; see '+str(path/'run.log'))
    summary=summarize(path); summary.update(meta);summary['path']=str(path)
    (path/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(name, summary['ncycle'],summary['time'],summary['min_ratio'],summary['max_ratio'],summary['min_pressure'],flush=True)
    return summary

if __name__=='__main__':
    all_runs=[]
    for stop in ['cycles','time']:
      for floor in ['1e-10','1e-14']:
       for v in [10,-10]:
        for face in ['wo1','original']:
         for wall in [True,False]:
          name=f'{face}_fh{int(wall)}_v{v:+d}_bf{floor}_{stop}'
          all_runs.append(run(name,face,wall,v,floor,stop))
          (R/'summary.json').write_text(json.dumps(all_runs,indent=2)+'\n')
    for nx in [128,256,512]:
      for wall in [True,False]:
        name=f'smooth_wo1_fh{int(wall)}_nx{nx}'
        all_runs.append(run(name,'wo1',wall,10,'1e-14','time',nx,.04))
        (R/'summary.json').write_text(json.dumps(all_runs,indent=2)+'\n')
