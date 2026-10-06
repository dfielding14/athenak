from pathlib import Path
import json,os,subprocess,hashlib,sys
here=Path(__file__).resolve().parent
sys.path.insert(0,str(here/'source-fused/vis/python'))
import athena_read,bin_convert
import numpy as np
root=here/'amr-cap-candidate';root.mkdir(exist_ok=True)
binary=here/'bin/athena-cpu';records=[]
for dimension,deck,ratio,nlim in [('2d','cgl_lf_amr_primitive_churn.athinput',.1,6),('3d','cgl_lf_amr_3d_current.athinput',.13,3)]:
 out=root/dimension;out.mkdir()
 command=['srun','--exact','--overlap','-N','1','-n','1','-c','1','--gpus-per-node=0',str(binary),'-i',str(here/'source-fused/inputs/tests'/deck),f'time/sts_max_dt_ratio={ratio}',f'time/nlim={nlim}']
 if dimension=='3d':command+=['output2/variable=mhd_w_bcc','output2/id=state']
 result=subprocess.run(command,cwd=out,env=os.environ.copy(),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
 (out/'stdout.log').write_text(result.stdout);(out/'command.json').write_text(json.dumps(command,indent=2)+'\n')
 user=athena_read.hst(str(next(out.glob('*.user.hst'))));mhd=athena_read.hst(str(next(out.glob('*.mhd.hst'))))
 record={'dimension':dimension,'cap':ratio,'returncode':result.returncode,'time':user['time'].tolist(),'ncell':user['ncell'].tolist(),'dt':user['dt'].tolist(),
  'health':{k:float(mhd[k][-1]) for k in ('lf_dfloor','lf_pfloor','lf_nonfin','lf_nonpos','lf_hardbd','lf_hwproj')},
  'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest()}
 if dimension=='3d':
  states=[bin_convert.read_binary(str(p)) for p in sorted((out/'bin').glob('*.state.*.bin'))]
  record['states']=[]
  for state in states:
   f=state['mb_data'];margin=np.asarray(f['p_perp'],float)-np.asarray(f['eint'],float)+sum(np.asarray(f[n],float)**2 for n in ('bcc1','bcc2','bcc3'))
   record['states'].append({'cycle':state['cycle'],'n_mbs':state['n_mbs'],'min_margin':float(np.min(margin)),'wall_cells':int(np.count_nonzero(np.abs(margin)<=2e-6))})
 records.append(record);(root/'results.json').write_text(json.dumps(records,indent=2)+'\n');print(record,flush=True)
