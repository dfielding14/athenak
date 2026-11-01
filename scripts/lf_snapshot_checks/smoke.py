"""Run both included pilot snapshots outside the source tree and check the handoff outputs."""
import argparse
import csv
from pathlib import Path
import subprocess
import sys
import numpy as np

home=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('out',type=Path,help='new directory for the pilot smoke outputs')
base=parser.parse_args().out.resolve()
base.mkdir(parents=True,exist_ok=False)

def run(*argv):
    subprocess.run([sys.executable,*argv],cwd=home,check=True)

run('-m','unittest','test_snapshot_io')
run('scale_analysis.py')
run('test_optical.py')
for name in ('nu0','collisional'):
    dest=base/name
    run('run.py','--snapshot',str(home/'examples'/f'{name}.bin'),'--case',name,
        '--la','.25','--ell-nu','.16667','--la-source','smoke fixture, arbitrary separation',
        '--visc-source','smoke fixture, arbitrary separation','--out',str(dest),
        '--cache-dir',str(base),'--reaches','.125','.25','.5','--step-cells','.5',
        '--anchors','10','--optical-anchors','4','--optical-reach','.125')
    from run import stencils, paired_values
    from snapshot_io import load_fields
    from fluxes import evaluate_faces
    fields,meta=load_fields(home/'examples'/f'{name}.bin')
    with np.load(home/'examples'/f'{name}_reference.npz') as ref:
        pairs=ref['points'].reshape(1,1,-1,2,3)
        xyz,endpoint,ids,inv=stencils(pairs,fields.shape[:3])
        dx=np.asarray(meta['spacing'])
        result=evaluate_faces(fields,ids,dx,meta['config'],[.125,.25,.5],dx.min()/4,[.1,100])
        assert np.array_equal(ids,ref['ids'])
        np.testing.assert_allclose(result['q'][:4],ref['q'],rtol=1e-9,atol=1e-11)
        _,td,_=paired_values(fields,xyz,endpoint,inv,result['q'][:4],
                             result['bn'],result['bf'],dx)
        np.testing.assert_allclose(td.reshape(4,-1),ref['rate'],rtol=1e-9,atol=1e-10)
run('combine.py',str(base/'nu0'),str(base/'collisional'),'--out',str(base/'combined'))
with (base/'combined'/'scales.csv').open() as f:
    rows=list(csv.DictReader(f))
assert {float(r['nu_coll']) for r in rows} == {0.,20.}
assert {float(r['factor']) for r in rows} == {.5,1.,2.}
assert any(r['family']=='path scattering' for r in rows)
assert any(r['resolution_flag']=='below 4 cells' for r in rows)
print('Both pilot cases and independent signed face-flux/tendency checks passed.')
