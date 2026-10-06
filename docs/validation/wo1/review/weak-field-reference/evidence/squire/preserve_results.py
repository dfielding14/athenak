from pathlib import Path
import subprocess,csv,json,hashlib,difflib
root=Path(__file__).resolve().parent
gitdir='/Users/dbf75/Work/Research/AthenaK/Athena++/AthenaCGL/.git_old_dev_version_history'
sha='44a29a70f2869b44e97d725393d4b9a020a60adf'
def original(p):
 return subprocess.check_output(['git','--git-dir='+gitdir,'show',sha+':'+p])
patch=[]
for p in ['src/pgen/b4_reference.cpp','src/task_list/time_integrator.cpp']:
 old=[] if 'b4_reference' in p else original(p).decode().splitlines(keepends=True)
 patch.extend(difflib.unified_diff(old,(root/p).read_text().splitlines(keepends=True),fromfile='/dev/null' if not old else 'a/'+p,tofile='b/'+p))
(root/'experiment.patch').write_text(''.join(patch))
files=subprocess.check_output(['git','--git-dir='+gitdir,'ls-tree','-r','--name-only',sha,'src']).decode().splitlines()
changed=[]
manifest={}
for p in files+['Makefile','bin/athena','run_cases.py','src/pgen/b4_reference.cpp']:
 f=root/p
 if not f.is_file(): continue
 data=f.read_bytes()
 manifest[p]=hashlib.sha256(data).hexdigest()
 if p in files and data!=original(p): changed.append(p)
(root/'source_and_binary_sha256.json').write_text(json.dumps(manifest,indent=2)+'\n')
summary=json.loads((root/'summary.json').read_text())
for r in summary:
 r['actual_floor_note']='Exact archive keeps bsq_rms=0, so effective EOS and face floors are zero for both configured controls.'
(root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
stages=[]
for case in ['magnetized_pos_50cycles','facefloor_pos_50cycles']:
 for cyc in range(1,6):
  for stage in [1,2]:
   def load(phase): return list(csv.DictReader((root/'runs'/case/f'state_{cyc:03d}_{phase}_{stage}.csv').open()))
   before,after=load('post_c2p'),load('post_p2c')
   worst=max(after,key=lambda r:float(r['pperp'])/float(r['ppar']))
   changes=[{'i':int(a['i']),'A_before':float(b['A']),'A_after':float(a['A']),'delta_A':float(a['A'])-float(b['A']),'delta_E':float(a['E'])-float(b['E']),'ppar':float(a['ppar']),'pperp':float(a['pperp'])} for b,a in zip(before,after) if float(a['A'])!=float(b['A']) or float(a['E'])!=float(b['E'])]
   stages.append({'case':case,'cycle':cyc,'stage':stage,'maxratio':float(worst['pperp'])/float(worst['ppar']),'worst_i':int(worst['i']),'B':float(worst['By']),'A':float(worst['A']),'ppar':float(worst['ppar']),'pperp':float(worst['pperp']),'p2c_changed_cells':changes})
(root/'first5_stage_summary.json').write_text(json.dumps(stages,indent=2)+'\n')
(root/'PROVENANCE.md').write_text('''# Squire archived reference B4 experiment

Immutable archived source: `'''+sha+'''`.

Archive command: `git --git-dir='''+gitdir+''' archive '''+sha+''' | tar -x -C /tmp/cgl-b4-compare-20261006/squire`.

Configure: `python3 configure.py --prob=b4_reference --coord=cartesian --eos=cgl --flux=hlle -b --cxx=clang++-apple`.

Build: `make -j2`. Commands and output are in configure.log, build.log, rebuild.log. Runs: `python3 run_cases.py`; full per-case commands/inputs/status are in summary.json and runs/*/.

Only custom pgen and read-only stage diagnostics were added. All archived physics files retain their exact bytes. Changed tracked source files: '''+repr(changed)+'''. The two stage calls dump after C2P and after collision/P2C, first five cycles only. Snapshot stage time is beginning-of-cycle time; cycle snapshots record end-of-cycle time.

## Actual magnetic floor behavior

The exact archive initializes Mesh::bsq_rms=0. Its src/hydro/new_blockdt.cpp:208 has `if (CGL_EOS && false)` around both an inactive field timestep constraint and the RMS accumulation (lines291–299); bsq_av initializes0 at65 and is assigned to block_bsq_rms at361. Thus global bsq_rms remains zero in these runs. GetBSqrFloor() returns zero regardless of configured bsqr_floor. Both EOS and face cutoffs are zero throughout, as recorded percell in CSVs. The nominal facefloor2e-20 and magnetized2e-28 controls are consequently identical. This behavior is preserved; no floor fix was applied. The initial physical mean B² of the sharp contact is .5 but is not reflected in archived bsq_rms.

The source additionally has an EOS-versus-face floor-unit inconsistency if RMS were nonzero: EOS compares |B| against bsqr_floor*bsq_rms, faces against its square root. This experiment does not exercise that branch because actual RMS stays zero.

## Initial conditions and outputs

Sharp: rho=1, v=±10, Bx=Bz=0, By=1e-12 on incoming side and1 on the outgoing side, isotropic p=1.5−By²/2. The v<0 state is the reflected fixture, matching current AthenaK. PLM RK2 CFL.4,128 cells,one block,outflow, optional kinetic limiters and nu_coll off, no heat flux. There is no Squire archived FOFC implementation. Initialization calls unmodified PrimitiveToConserved with actual B; weak A/rho=82.89306334778564.

Smooth conditional follow-up: v=10, By=1e-12+(1−1e-12)*(1+tanh((x−.5)/.04))/2, p=1.5−By²/2; nx128/256/512, fixed t=.012. All states above actual zero floor. Raw percell everycycle data and first-five-stage dumps retain rho,momentum,E,A,pressures,B and actual floor diagnostics.

No files in the archived reference directory or production worktree were modified. Source/binary SHA256 manifest and experiment.patch are retained.
''')
print('Changed archived source:',changed)
for s in stages[:10]: print({k:v for k,v in s.items() if k!='p2c_changed_cells'},'changed',len(s['p2c_changed_cells']))
