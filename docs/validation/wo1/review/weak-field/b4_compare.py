"""Run the current/pre-B4 grid comparison and retain every-cycle measurements.

Usage: python b4_compare.py REPOSITORY OUTPUT_DIRECTORY
"""

from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys

import numpy as np

root = Path(sys.argv[1]).resolve()
base = Path(sys.argv[2]).resolve() / 'b4-comparison'
base.mkdir(parents=True, exist_ok=True)
results = []
executables = [('current', root / 'tst/build/src/athena'),
               ('pre_b4', base.parent / 'pre-b4/athena')]
for version, executable in executables:
    for reconstruction in ('dc', 'plm'):
        for velocity in (10, -10):
            case = base / f'{version}-{reconstruction}-{velocity}'
            case.mkdir(exist_ok=True)
            command = [
                str(executable), '-i',
                str(root / 'inputs/unit_tests/cgl_weak_field_transport.athinput'),
                f'problem/weak_field_velocity={velocity}',
                f'mhd/reconstruct={reconstruction}', 'output1/variable=mhd_w_bcc',
            ]
            run = subprocess.run(command, cwd=case, text=True, capture_output=True)
            (case / 'run.log').write_text(run.stdout + run.stderr)
            record = {
                'version': version, 'reconstruct': reconstruction,
                'velocity': velocity, 'command': command,
                'executable_sha256': hashlib.sha256(executable.read_bytes()).hexdigest(),
                'returncode': run.returncode, 'cycles': [],
            }
            for path in sorted((case / 'tab').glob('*.tab')):
                data = np.loadtxt(path)
                cycle = int(re.search(r'cycle=(\d+)', path.read_text())[1])
                ratio = data[:, 8] / data[:, 7]
                j = ratio.argmax()
                record['cycles'].append({
                    'cycle': cycle, 'min_ratio': float(ratio.min()),
                    'max_ratio': float(ratio[j]), 'x': float(data[j, 2]),
                    'B_at_max': float(data[j, 10]),
                    'ppar_at_max': float(data[j, 7]),
                    'pperp_at_max': float(data[j, 8]),
                })
            results.append(record)
            print(version, reconstruction, velocity, 'exit', run.returncode,
                  'cycle1', record['cycles'][1], 'final', record['cycles'][-1],
                  flush=True)
(base / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
