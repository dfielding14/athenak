from pathlib import Path
import sys,os,subprocess,json
repo=Path('/Users/dbf75/.codex/worktrees/bf22/athenak-DF')
root=Path('/tmp/cgl-review-shear')
os.chdir(repo/'tst');sys.path.insert(0,str(repo/'tst'))
from test_suite.cgl import test_cgl_lf_sbox_mpicpu as shear
from test_suite.cgl import test_cgl_landau_fluid_cpu as cpu
binary=repo/'build-review-shear/src/athena'
result={}
for name,flags in [('safe',[]),('fast-profile',[]),('explicit',['mhd/cgl_heat_flux_integrator=explicit','time/sts_integrator=none'])]:
 run=root/('check-'+name);run.mkdir(exist_ok=True);os.chdir(run)
 env=os.environ.copy()
 if name=='fast-profile': env.update(ATHENAK_CGL_LF_ARITHMETIC='fast',ATHENAK_CGL_LF_PROFILE='1',ATHENAK_CGL_LF_PROFILE_DETAIL='1')
 p=subprocess.run([str(binary),'-i',str(repo/'tst/inputs/cgl_lf_sts_sbox.athinput'),*flags],capture_output=True,text=True,env=env)
 (run/'stdout.log').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stdout[-2000:]
 output=shear._read('cgl_lf_sts_sbox','mhd_w_bcc');hst=shear.testutils.athena_read.hst('cgl_lf_sts_sbox.mhd.hst')
 shear._assert_admissible(output,hst)
 shear._assert_magnetic_state(output,shear._read('cgl_lf_sts_sbox','mhd_divb'))
 result[name]={'cycle':int(output['cycle']),'time':float(output['time']),'nstage':float(hst['lf_nstage'][-1]),'hard_bound':float(hst['lf_hardbd'][-1])}
print(json.dumps(result,indent=2));(root/'serial-results.json').write_text(json.dumps(result,indent=2)+'\n')
# Run retained pytest negative test directly against the tested binary.
run=root/'pytest-entry';run.mkdir(exist_ok=True);os.chdir(run)
(run/'athena').symlink_to(binary) if not (run/'athena').exists() else None
cpu.INPUT_ROOT=str(repo/'inputs/tests')
for backup in ('false','true'):
 for integ in ('sts','explicit'):
  cpu.test_cgl_lf_strict_initial_hard_bound_is_rejected_at_sweep_entry(backup,integ)
print('four invalid-entry negative tests passed')
