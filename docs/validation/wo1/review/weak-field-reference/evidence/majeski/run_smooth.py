from pathlib import Path
import subprocess,time,json
base=Path(__file__).resolve().parent
exe=base/'build/src/athena'
text=(base/'base.athinput').read_text().replace('bfloor = 1.0e-10','bfloor = 1.0e-14').replace('nlim = 50','nlim = 1000').replace('tlim = 100.0','tlim = 0.012').replace('<problem>','<problem>\nprofile_width = 0.04')
results=[]
for nx in (128,256,512):
  name=f'smooth_nx{nx}_bf1.0e-14_v10.0_time012'
  out=base/'runs'/name;out.mkdir(parents=True,exist_ok=True)
  (out/'input.athinput').write_text(text.replace('nx1 = 128',f'nx1 = {nx}'))
  cmd=[str(exe),'-i','input.athinput']
  start=time.time()
  with (out/'run.log').open('w') as log:
    run=subprocess.run(cmd,cwd=out,stdout=log,stderr=subprocess.STDOUT,timeout=60)
  result=dict(name=name,command=cmd,returncode=run.returncode,walltime=time.time()-start)
  results.append(result);print(json.dumps(result),flush=True)
(base/'smooth_results.json').write_text(json.dumps(results,indent=2)+'\n')
