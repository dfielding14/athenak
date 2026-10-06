from pathlib import Path
import subprocess, time, json, hashlib
base=Path(__file__).resolve().parent
exe=base/'build/src/athena'
text=(base/'base.athinput').read_text()
results=[]
for fofc in (False, True):
  for floor in ('1.0e-10','1.0e-14'):
    for velocity in ('10.0','-10.0'):
      for stop in ('cycles50','time012'):
        name=f"{'fofc' if fofc else 'nofc'}_bf{floor}_v{velocity}_{stop}"
        out=base/'runs'/name
        out.mkdir(parents=True,exist_ok=True)
        inp=text.replace('bfloor = 1.0e-10',f'bfloor = {floor}').replace('weak_field_velocity = 10.0',f'weak_field_velocity = {velocity}')
        if fofc: inp=inp.replace('fofc = false','fofc = true')
        if stop=='time012': inp=inp.replace('nlim = 50','nlim = 1000').replace('tlim = 100.0','tlim = 0.012')
        (out/'input.athinput').write_text(inp)
        cmd=[str(exe),'-i','input.athinput']
        start=time.time()
        try:
          with (out/'run.log').open('w') as log:
            run=subprocess.run(cmd,cwd=out,stdout=log,stderr=subprocess.STDOUT,timeout=60)
          rc=run.returncode
        except subprocess.TimeoutExpired:
          rc='timeout'
        result=dict(name=name,command=cmd,returncode=rc,walltime=time.time()-start)
        results.append(result)
        print(json.dumps(result),flush=True)
        (base/'run_results.json').write_text(json.dumps(results,indent=2)+'\n')
