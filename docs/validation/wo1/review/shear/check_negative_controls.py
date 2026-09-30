from pathlib import Path
import re,shlex,subprocess,json
repo=Path('/Users/dbf75/.codex/worktrees/bf22/athenak-DF')
root=Path('/tmp/cgl-review-shear')
build=repo/'build-review-shear'
source=(repo/'src/mhd/mhd_sts.cpp').read_text()
needle='''    peos->Collisions(w0, bcc0, u0, pdrive->sts.dt_cycle, mode,
                     0, n1m1, 0, n2m1, 0, n3m1);'''
assert source.count(needle)==1
mutant=root/'mhd_sts.no-sweep-wall.cpp'
mutant.write_text(source.replace(needle,'    // Negative control: omit the scheduled rates/wall update.'))
flags=(build/'src/CMakeFiles/athena.dir/flags.make').read_text()
args=[]
for name in ('CXX_DEFINES','CXX_INCLUDES','CXX_FLAGS'):
 args+=shlex.split(re.search(r'^'+name+r' = (.*)$',flags,re.M)[1])
obj=root/'mhd_sts.no-sweep-wall.o'
subprocess.run(['/usr/bin/c++',*args,'-I'+str(repo/'src/mhd'),'-c',str(mutant),'-o',str(obj)],check=True)
link=shlex.split((build/'src/CMakeFiles/athena.dir/link.txt').read_text())
link=[str(obj) if a=='CMakeFiles/athena.dir/mhd/mhd_sts.cpp.o' else a for a in link]
binary=root/'athena-no-sweep-wall'
link[link.index('-o')+1]=str(binary)
subprocess.run(link,cwd=build/'src',check=True)
run=root/'negative-no-wall';run.mkdir(exist_ok=True)
p=subprocess.run([str(binary),'-i',str(repo/'tst/inputs/cgl_lf_sts_sbox.athinput'),'time/ndiag=1'],cwd=run,capture_output=True,text=True)
(run/'stdout.log').write_text(p.stdout+p.stderr)
assert p.returncode!=0
assert 'wall checkpoint: sweep=post stage=21/21' in p.stdout
print(p.stdout.splitlines()[-1])
# The unmodified initial-state contract is retained for either backup policy.
for backup in ('false','true'):
 run=root/('negative-initial-'+backup);run.mkdir(exist_ok=True)
 p=subprocess.run([str(build/'src/athena'),'-i',str(repo/'inputs/tests/cgl_lf_firehose_policy.athinput'),'mhd/cgl_firehose_threshold=parallel','mhd/backup_limiters='+backup,'problem/ppar0=3.0','problem/pperp0=1.0'],cwd=run,capture_output=True,text=True)
 (run/'stdout.log').write_text(p.stdout+p.stderr)
 assert p.returncode!=0
 assert 'wall checkpoint: sweep=pre stage=0/' in p.stdout
 print(backup,p.stdout.splitlines()[-1])
