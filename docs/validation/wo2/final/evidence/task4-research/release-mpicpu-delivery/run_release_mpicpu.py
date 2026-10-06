#!/usr/bin/env python3
"""Run the permanent CPU-MPI passive package on the exact relinked release binary."""
from pathlib import Path
import datetime
import hashlib
import json
import os
import subprocess
import sys

root = Path(__file__).resolve().parent
final = root.parent/'final-research'
out = root/'release-mpicpu-runs'/sys.argv[1]
assert out.is_dir() and out.resolve().is_relative_to(root)
assert os.environ['SLURM_JOB_ID']
assert os.environ['MPICH_GPU_SUPPORT_ENABLED'] == '0'
source = final/'source-release'
binary = Path(os.environ['ATHENAK_CGL_PASSIVE_BINARY']).resolve()
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
link = final/'release-cpu-link-manifest.json'
link_record = json.loads(link.read_text())
assert str(binary) == link_record['binary'] and sha(binary) == link_record['sha256']
assert 'libamdhip' not in link_record['dynamic_section']
assert 'libmpi_gtl_hsa' not in link_record['dynamic_section']
module = source/'tst/test_suite/cgl/test_cgl_passive_mpicpu.py'
paths = [module, source/'tst/test_suite/cgl/passive_acceptance.py',
         root/'run_release_mpicpu.sh', Path(__file__).resolve(),
         final/'release-source-manifest.json', link]
cmd = [sys.executable, '-m', 'pytest', '-q', '-s', '-x', str(module),
       '--basetemp='+str(out/'work')]
record = {'purpose': 'New four-rank CPU permanent passive acceptance, separate from prior CPU1/HIP1/HIP4 counts',
          'started': datetime.datetime.now(datetime.timezone.utc).isoformat(),
          'backend': 'cpu', 'ranks': 4, 'slurm_job': os.environ['SLURM_JOB_ID'],
          'binary': str(binary), 'binary_sha256': sha(binary),
          'source': str(source), 'command': cmd, 'launcher': os.environ['ATHENAK_CGL_PASSIVE_MPI_LAUNCHER'],
          'expected_pytest_groups': 6, 'expected_applications': 80, 'expected_checks': 65,
          'files': {str(p): sha(p) for p in paths},
          'environment': {k:v for k,v in os.environ.items()
             if k.startswith(('MPICH_', 'HSA_', 'FI_', 'ROCR_', 'HIP_', 'KOKKOS_'))
             or k in ('LD_LIBRARY_PATH','OMP_NUM_THREADS','TMPDIR','PYTHONPATH')},
          'status': 'running'}
def save():
    (out/'manifest.json').write_text(json.dumps(record, indent=2)+'\n')
save()
print('Running six permanent CPU-MPI groups; log:',out/'pytest.log',flush=True)
with (out/'pytest.log').open('w') as log:
    proc = subprocess.run(cmd, cwd=root, env=os.environ.copy(), stdout=log, stderr=subprocess.STDOUT)
record['returncode'] = proc.returncode
record['finished'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
assert sha(binary) == record['binary_sha256']
record['binary_hash_unchanged'] = True
record['groups'] = {}
for p in sorted((out/'work').glob('test_cgl_passive_mpicpu*')):
    if p.is_symlink() or not (p/'results.json').exists(): continue
    result = json.loads((p/'results.json').read_text())
    record['groups'][p.name] = {'result':str(p/'results.json'),
                              'sha256':sha(p/'results.json'),
                              'applications':len(result['records']),
                              'checks':len(result['checks'])}
record['applications'] = sum(g['applications'] for g in record['groups'].values())
record['checks'] = sum(g['checks'] for g in record['groups'].values())
record['status'] = 'passed' if proc.returncode == 0 else 'failed'
save()
if proc.returncode == 0:
    assert len(record['groups']) == 6 and record['applications'] == 80 and record['checks'] == 65
print(json.dumps({k:record[k] for k in ['status','returncode','applications','checks']},indent=2),flush=True)
raise SystemExit(proc.returncode)
