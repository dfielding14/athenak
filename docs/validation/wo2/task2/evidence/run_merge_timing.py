"""Paired GPU timing on a separately specified eligible smooth LF fixture.

The original strict tests remain unchanged. This timing fixture explicitly uses
nonstrict, collisionless evolution, with the original rotated-decay amplitude and
analytic error criteria; every run must pass those checks and have zero repairs.
"""
from pathlib import Path
import argparse,hashlib,json,math,os,re,shutil,statistics,subprocess,sys,time
import numpy as np
here=Path(__file__).resolve().parent
repo=Path('/autofs/nccs-svm1_home2/dfielding/athenak-cgl-wo2-baseline')
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--binary',type=Path,required=True)
parser.add_argument('--manifest',type=Path,required=True)
parser.add_argument('--label',required=True)
parser.add_argument('--plan-only',action='store_true')
args=parser.parse_args()
assert re.fullmatch('[A-Za-z0-9_.-]+',args.label)
root=here/'timing'/args.label
root.mkdir(parents=True,exist_ok=False)
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def stage(flags):
 text=(repo/'inputs/unit_tests/cgl_lf_rotated_decay.athinput').read_text()
 text+='\n<output2>\nfile_type = bin\nvariable = mhd_w_bcc\ndcycle = 100000\n<output3>\nfile_type = rst\ndcycle = 100000\n'
 for key,value in flags.items():
  block,name=key.split('/',1)
  m=re.search(r'(^<'+re.escape(block)+r'>\s*\n)(.*?)(?=^<|\Z)',text,re.M|re.S)
  assert m,key
  body=m.group(2);pattern=re.compile(r'^'+re.escape(name)+r'\s*=.*$',re.M)
  body=pattern.sub(name+' = '+str(value),body) if pattern.search(body) else body+'\n'+name+' = '+str(value)+'\n'
  text=text[:m.start(2)]+body+text[m.end(2):]
 return text
base_flags={'job/basename':'merge_timing','time/nlim':96,'time/tlim':1000,
 'time/ndiag':100000,'time/cfl_number':0.3,'mhd/cgl_lf_strict_admissibility':'false',
 'mhd/nu_coll':0,'mhd/mirror_limiter':'false','mhd/firehose_limiter':'false',
 'mhd/limiter_nu_coll':0,'mhd/backup_limiters':'false',
 'mhd/cgl_lf_profile':'false','mhd/cgl_lf_profile_detail':'false',
 'mhd/cgl_lf_arithmetic':'safe','mhd/cgl_lf_diagnostics':'full','mhd/cgl_lf_sts_flux':'weighted',
 'output1/dcycle':100000,'output1/data_format':'%24.16e'}
for n in (1,2,3):base_flags['mesh/nx'+str(n)]=64;base_flags['meshblock/nx'+str(n)]=32
plans=[]
for ratio in (2,20):
 for enabled in (False,True):
  flags=base_flags|{'time/sts_max_dt_ratio':ratio,'time/sts_merge_half_sweeps':str(enabled).lower()}
  text=stage(flags)
  # Keep all scientific signal/analytic tolerances from the source fixture.
  for line in ['amp = 1.0e-4','decay_rel_tol = 3.0e-2','phase_rel_tol = 2.0e-2']:
   assert line in text
  path=root/f'ratio{ratio}-merge{enabled}.athinput';path.write_text(text)
  plans.append({'ratio':ratio,'enabled':enabled,'input':str(path),'sha256':sha(path)})
result={'classification':'separate eligible smooth GPU fixture, not original strict paper acceptance',
 'input_source':str(repo/'inputs/unit_tests/cgl_lf_rotated_decay.athinput'),
 'geometry':[64,64,64],'meshblock':[32,32,32],'cycles':96,
 'signal_and_analytic_tolerances':'original amp1e-4, decay3e-2, phase2e-2 unchanged',
 'intentional_scope':'nonstrict collisionless uniform periodic kinematic LF, 3D volume with x-directed thermal wave',
 'ratio_cases':[2,20],'interpretation':'ratio2 exercises minimum three-stage ordinary halves; ratio20 exercises larger sweeps',
 'timing':'solver reported seconds per cycle, including original LF diagnostics and sparse initial/final output; warmups excluded; off/on order alternates',
 'transaction_cost':'reported snapshot_seconds; 96 bytes per stored ghosted cell for u0/w0 snapshots; accepted-copy traffic192 bytes per stored cell',
 'plans':plans,'results':[],'failures':[],'runner_sha256':sha(__file__)}
def save(): (root/'results.json').write_text(json.dumps(result,indent=2)+'\n')
save()
if args.plan_only:print(root);raise SystemExit(0)
assert os.environ.get('SLURM_JOB_ID')
def gate():
 if (here/'PAUSE_APPLICATION_TESTS').exists():raise SystemExit('Task2 applications remain paused')
gate()
binary=root/'athena';shutil.copy2(args.binary,binary)
result.update(binary_source=str(args.binary.resolve()),binary_sha256=sha(binary),
 source_manifest=str(args.manifest.resolve()),source_provenance=json.loads(args.manifest.read_text()),
 source_manifest_sha256=sha(args.manifest),slurm_job=os.environ['SLURM_JOB_ID'])
def history(path):
 lines=path.read_text().splitlines();labels=re.findall(r'\[\d+\]=(\S+)',next(x for x in lines if '[1]=' in x))
 rows=[dict(zip(labels,map(float,x.split()))) for x in lines if x.strip() and not x.startswith('#')]
 assert rows and all(math.isfinite(v) for row in rows for v in row.values())
 return rows
for ratio in (2,20):
 for repeat in (-1,0,1,2):
  for enabled in ((False,True) if repeat%2==0 else (True,False)):
   gate()
   run=root/f'ratio{ratio}-repeat{repeat}-merge{enabled}';run.mkdir()
   plan=next(x for x in plans if x['ratio']==ratio and x['enabled']==enabled)
   command=['srun','--exact','-N1','-n1','--threads-per-core=1','--cpu-bind=threads','-c7',
            '--gpus-per-task=1','--gpu-bind=closest',str(binary),'-i',plan['input']]
   env=os.environ.copy()
   for key in list(env):
    if key.startswith('ATHENAK_CGL_LF_') or key in ('KOKKOS_TOOLS_LIBS','KOKKOS_PROFILE_LIBRARY'):env.pop(key)
   (run/'tmp').mkdir();env['TMPDIR']=str(run/'tmp');env['OMP_NUM_THREADS']='1'
   (run/'command.json').write_text(json.dumps(command,indent=2)+'\n')
   start=time.perf_counter()
   with (run/'stdout.log').open('w') as log:
    proc=subprocess.run(command,cwd=run,env=env,stdout=log,stderr=subprocess.STDOUT)
   sample=dict(ratio=ratio,repeat=repeat,enabled=enabled,timed=repeat>=0,
               returncode=proc.returncode,launch_inclusive_seconds=time.perf_counter()-start)
   try:
    assert proc.returncode==0,'application/analytic acceptance failure'
    stdout=(run/'stdout.log').read_text()
    assert 'CGL LF rotated_decay passed:' in stdout
    cycles=int(re.findall(r'^time=\S+ cycle=(\d+)',stdout,re.M)[-1]);assert cycles==96
    rows=history(next(run.glob('*.mhd.hst')))
    for name in ('lf_dfloor','lf_pfloor','lf_nonfin','lf_nonpos','lf_hwproj'):assert rows[-1][name]==0,name
    stats={k:float(v) for k,v in re.findall(r'(accepted|rejected|rejected_cfl|rejected_admissibility|accepted_stages|rejected_stages|snapshot_seconds|deferred|consumed|flushed|pending)=([+\-.0-9eE]+)',stdout)}
    if enabled:
     assert stats['accepted']>0 and stats['pending']==0
    else:assert not stats
    seconds=float(re.findall(r'^cpu time used\s*=\s*(\S+)',stdout,re.M)[-1])
    stages=(rows[-1]['lf_nstage']-rows[0]['lf_nstage'])/(64**3)
    sample.update(cycles=cycles,solver_seconds_per_cycle=seconds/cycles,stages=stages,
                  stages_per_cycle=stages/cycles,transaction=stats,history_final=rows[-1],
                  analytic_relative_error=float(re.findall(r'rotated_decay passed:.* rel_err=(\S+)',stdout)[-1]),
                  mass_relative_error=abs(rows[-1]['mass']/rows[0]['mass']-1),
                  output_hashes={str(p.relative_to(run)):sha(p) for p in run.rglob('*') if p.is_file() and p.suffix in ('.bin','.rst','.hst')})
    assert sample['mass_relative_error']<1e-10
   except (AssertionError,ValueError,KeyError,IndexError,StopIteration) as error:
    sample['error']=str(error);result['failures'].append({'run':str(run),'reason':str(error)})
   result['results'].append(sample);save();print(run.name,{key:sample.get(key) for key in ('returncode','solver_seconds_per_cycle','stages_per_cycle','transaction','error')},flush=True)
   if 'error' in sample:raise SystemExit(1)
result['medians']={}
for ratio in (2,20):
 modes={}
 for enabled in (False,True):
  samples=[x for x in result['results'] if x['ratio']==ratio and x['enabled']==enabled and x['timed']]
  assert len(samples)==3
  modes[str(enabled)]={key:statistics.median(x[key] for x in samples) for key in ('solver_seconds_per_cycle','stages_per_cycle','launch_inclusive_seconds')}
 result['medians'][str(ratio)]={'off':modes['False'],'on':modes['True'],
    'speedup_off_over_on':modes['False']['solver_seconds_per_cycle']/modes['True']['solver_seconds_per_cycle']}
save()
