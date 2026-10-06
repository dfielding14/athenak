#!/usr/bin/env python3
"""Stage compact CPU4 follow-up evidence without altering original Task4 reports."""
from pathlib import Path
import hashlib
import json
import shutil
import sys

root=Path(__file__).resolve().parent
run=root/'release-mpicpu-runs'/sys.argv[1]
manifest=json.loads((run/'manifest.json').read_text())
assert manifest['status']=='passed' and manifest['returncode']==0
assert len(manifest['groups'])==6 and manifest['applications']==80 and manifest['checks']==65
out=root/'release-mpicpu-delivery'
assert not out.exists()
out.mkdir()
for name in ['manifest.json','pytest.log','modules.txt']:
    shutil.copy2(run/name,out/name)
for name in ['run_release_mpicpu.sh','run_release_mpicpu.py','stage_release_mpicpu.py']:
    shutil.copy2(root/name,out/name)
(out/'groups').mkdir()
for group,data in manifest['groups'].items():
    path=Path(data['result'])
    assert hashlib.sha256(path.read_bytes()).hexdigest()==data['sha256']
    shutil.copy2(path,out/'groups'/f'{group}.json')
(out/'README.md').write_text('''# Final release four-rank CPU passive acceptance

The six permanent `test_cgl_passive_mpicpu.py` groups pass on four CPU MPI ranks:
80 application launches and 65 acceptance checks. This is additional CPU-MPI
coverage; the earlier CPU1/HIP1/HIP4 result counts remain unchanged.

The exact final relinked CPU binary is recorded by path and SHA-256 in
`manifest.json`. The runner uses CPU-only Cray modules,
`MPICH_GPU_SUPPORT_ENABLED=0`, four MPI ranks and zero GPUs. The original CPU
release binary with unintended ROCm link dependencies was not used. Every
application's working directory, temporary directory and outputs reside under
WO2. These concurrent functional runs support no timing claims.

Coverage includes all five released reconstructors, exact native-isothermal
flow and timestep identity, forcing in three dimensions, floor/weak-field
conditions, independent periodic heating, linear response, thermal advection,
full-state resumed restart identity, and explicit unsupported-mode fences.
The package ran its unchanged assertions on the frozen `source-release` tests.

`pytest.log` records all six passing groups. `groups/` contains the complete
per-group result records and commands; the full-precision restart artifacts
remain at the paths in those records. `manifest.json` includes source and
runner hashes, binary identity, loaded runtime contract and completed counts.
''')
files={str(p.relative_to(out)):{'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}
       for p in sorted(out.rglob('*')) if p.is_file()}
(out/'artifact-manifest.json').write_text(json.dumps({'source_run':str(run),'files':files},indent=2)+'\n')
print(out)
