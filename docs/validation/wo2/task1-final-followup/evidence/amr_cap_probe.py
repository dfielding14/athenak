from pathlib import Path
import json,os,subprocess,re,hashlib
here=Path(__file__).resolve().parent
root=here/'amr-cap-probe';root.mkdir(exist_ok=True)
binary=here/'bin/athena-cpu'
records=[]
for dimension,deck in [('2d','cgl_lf_amr_primitive_churn.athinput'),('3d','cgl_lf_amr_3d_current.athinput')]:
 for ratio in [1,.1]:
  name=f'{dimension}-ratio{ratio}-initial';out=root/name;out.mkdir()
  command=['srun','--exact','--overlap','-N','1','-n','1','-c','1','--gpus-per-node=0',str(binary),'-i',str(here/'source-fused/inputs/tests'/deck),'time/nlim=0',f'time/sts_max_dt_ratio={ratio}']
  result=subprocess.run(command,cwd=out,env=os.environ.copy(),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
  (out/'stdout.log').write_text(result.stdout)
  (out/'command.json').write_text(json.dumps(command,indent=2)+'\n')
  match=re.search(r'cycle=0 time=\S+ dt=(\S+)',result.stdout)
  hist=next(out.glob('*.mhd.hst'));row=[x for x in hist.read_text().splitlines() if x and not x.startswith('#')][0].split()
  record={'case':name,'returncode':result.returncode,'dt':float(row[1]),'sha256':hashlib.sha256(binary.read_bytes()).hexdigest()}
  records.append(record);print(record,flush=True)
  assert result.returncode==0
  (root/'results.json').write_text(json.dumps(records,indent=2)+'\n')
