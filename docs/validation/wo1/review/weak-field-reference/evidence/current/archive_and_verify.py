from pathlib import Path
import difflib, hashlib, json, os, subprocess, csv
R=Path(__file__).parent; repo='/Users/dbf75/.codex/worktrees/bf22/athenak-DF'
sha=(R/'base-head.txt').read_text().strip()
changed=['src/mhd/b4_trace.hpp','src/mhd/mhd_tasks.cpp','src/mhd/mhd_update.cpp','src/mhd/mhd_ct.cpp','src/eos/cgl_mhd.cpp','src/mhd/rsolvers/hlle_cgl.hpp','src/mhd/rsolvers/llf_mhd_singlestate.hpp','src/eos/ideal_c2p_mhd.hpp','src/eos/b4_controls.hpp','src/pgen/tests/cgl_fofc.cpp']
patch=[]
for name in changed:
 old=subprocess.run(['git','show',sha+':'+name],cwd=repo,capture_output=True)
 oldtext=old.stdout.decode() if old.returncode==0 else ''
 new=(R/'source'/name).read_text()
 patch.extend(difflib.unified_diff(oldtext.splitlines(True),new.splitlines(True),fromfile='a/'+name if old.returncode==0 else '/dev/null',tofile='b/'+name))
(R/'instrumentation-and-controls.patch').write_text(''.join(patch))
verification=[]
env={k:v for k,v in os.environ.items() if not k.startswith('B4_')}
for v in [10,-10]:
 for floor in ['1e-10','1e-14']:
  name=f'baseline_v{v:+d}_bf{floor}';d=R/'verification'/name;d.mkdir(exist_ok=True)
  cmd=[str(R/'athena-uninstrumented-7b3345fd'),'-i',str(R/'comparison.athinput'),f'problem/weak_field_velocity={v}',f'mhd/bfloor={floor}','output1/variable=mhd_w_bcc','time/nlim=50']
  with (d/'run.log').open('w') as log:subprocess.run(cmd,cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
  target=R/'runs'/f'wo1_fh1_v{v:+d}_bf{floor}_cycles'
  files=sorted((d/'tab').glob('*.tab'))
  diff=[p.name for p in files if p.read_bytes()!=(target/'tab'/p.name).read_bytes()]
  item=dict(baseline=str(d),instrumented=str(target),command=cmd,nfiles=len(files),differing_files=diff)
  verification.append(item);assert not diff,item
# Original and WO1 have identical arithmetic when every face is magnetized.
allmag=[]
for wall in [0,1]:
 for v in [10,-10]:
  for stop in ['cycles','time']:
   a=R/'runs'/f'wo1_fh{wall}_v{v:+d}_bf1e-14_{stop}'
   b=R/'runs'/f'original_fh{wall}_v{v:+d}_bf1e-14_{stop}'
   files=sorted((a/'tab').glob('*.tab'));diff=[p.name for p in files if p.read_bytes()!=(b/'tab'/p.name).read_bytes()]
   allmag.append(dict(wo1=str(a),original=str(b),nfiles=len(files),differing_files=diff));assert not diff
# Validate trace count and exact requested common times.
runs=json.loads((R/'summary.json').read_text())
for run in runs:
 p=Path(run['path']);counts={}
 with (p/'b4-state.csv').open() as f:
  for row in csv.DictReader(f):
   if row['phase']=='coll_post':counts[int(row['cycle'])]=counts.get(int(row['cycle']),0)+1
 assert len(counts)==run['ncycle'] and set(counts.values())=={run['nx']},p
 if run['stop']=='time':assert run['time']==.012,p
result=dict(baseline_checks=verification,allmag_face_checks=allmag,trace_counts_passed=True,physical_time_checks_passed=True)
(R/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
hashes={}
for name in ['athena-uninstrumented-7b3345fd','athena-instrumented-4way','comparison.athinput','instrumentation-and-controls.patch','prepare_instrumented.py','run_comparison.py','archive_and_verify.py']:
 hashes[name]=hashlib.sha256((R/name).read_bytes()).hexdigest()
(R/'sha256.json').write_text(json.dumps(hashes,indent=2)+'\n')
print('Verified 4 baselines, 8 allmag face pairs, full per-cycle trace counts, and 22 exact time endpoints.')
