from pathlib import Path
import os,subprocess,hashlib,json
base=Path('/tmp/cgl-b4-transition-20261001')
inp=Path('/Users/dbf75/.codex/worktrees/bf22/athenak-DF/inputs/unit_tests/cgl_weak_field_transport.athinput')
exe=base/'build/src/athena'
records=[]
for name,recon,vel,floor in [('plm', 'plm',10,'1e-10'),('dc','dc',10,'1e-10'),('reverse','plm',-10,'1e-10'),('floor6','plm',10,'1e-6'),('all_magnetized','plm',10,'1e-14')]:
    work=base/name;work.mkdir(exist_ok=True)
    args=['-i',str(inp),'time/nlim=5',f'mhd/reconstruct={recon}',f'problem/weak_field_velocity={vel}',f'mhd/bfloor={floor}','output1/variable=mhd_w_bcc']
    with (work/'run.log').open('w') as f:
        p=subprocess.run([str(exe),*args],cwd=work,stdout=f,stderr=subprocess.STDOUT,env={**os.environ,'B4_TRACE':'1'})
    assert p.returncode==0,(name,p.returncode)
    records.append(dict(name=name,args=args,returncode=p.returncode))
    print(name,'complete',flush=True)
    if name=='plm':
        ref=base/'uninstrumented';ref.mkdir(exist_ok=True)
        with (ref/'run.log').open('w') as f:
            p=subprocess.run(['/Users/dbf75/.codex/worktrees/bf22/athenak-DF/tst/build/src/athena',*args],cwd=ref,stdout=f,stderr=subprocess.STDOUT)
        assert p.returncode==0
        hashes=lambda root:{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in (root/'tab').glob('*.tab')}
        assert hashes(ref)==hashes(work),'Instrumentation changed output'
        print('Instrumented/uninstrumented tab outputs byte-identical',flush=True)
(base/'trace-commands.json').write_text(json.dumps(records,indent=2)+'\n')
