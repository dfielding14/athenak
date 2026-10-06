#!/usr/bin/env python3
"""Final-state release equivalence only: no face traces or performance claims."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

import numpy as np

W = Path('/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2')
F = W / 'final-research'
sys.path.insert(0, str(W / 'p1-research'))
sys.path.insert(0, str(F / 'source-unfused/vis/python'))
from restart_compare import restart_hashes
from athena_read import hst


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--backend', required=True, choices=('cpu', 'hip'))
parser.add_argument('--output', type=Path)
parser.add_argument('--new-binary', type=Path)
parser.add_argument('--cases', nargs='+', help='Exact case names; default all 25')
parser.add_argument('--list-cases', action='store_true')
parser.add_argument('--node', help='Optional allocation node constraint')
args = parser.parse_args()

prior = F / 'task7-fusion/checks-5629018'
prior_results = json.loads((prior / 'results.json').read_text())
assert len(prior_results['cases']) == 23 and not prior_results['failures']
cases = []
for name, entry in prior_results['cases'].items():
    old_run = entry['runs']['fused']
    cmd = old_run['command']
    path = prior / 'inputs' / (name + '.athinput')
    assert sha(path) == old_run['input_sha256']
    cases.append(dict(name=name, input=path, ranks=int(cmd[cmd.index('-n') + 1]),
                      forced=False, input_sha256=sha(path)))
forced = F / 'release-equivalence-fixtures/forced-passive3d.athinput'
for ranks in (1, 4):
    cases.append(dict(name=f'forced-passive3d-safe-full-{ranks}rank', input=forced,
                      ranks=ranks, forced=True, input_sha256=sha(forced)))
if args.cases:
    assert set(args.cases) <= {row['name'] for row in cases}, 'unknown case name'
    cases = [row for row in cases if row['name'] in args.cases]
assert cases
if args.list_cases:
    print(json.dumps([{**row, 'input': str(row['input'])} for row in cases], indent=2))
    raise SystemExit(0)

assert args.output is not None, '--output is required for applications'
assert os.environ.get('SLURM_JOB_ID'), 'requires parent-managed allocation'
out = args.output.resolve()
assert out.is_relative_to(W) and not out.exists(), 'fresh output under WO2 required'
old_manifest_path = F / f'bin/manifest-{args.backend}.json'
old_manifest = json.loads(old_manifest_path.read_text())
old_binary = Path(old_manifest['path']).resolve()
new_binary = (args.new_binary or F / f'release-build-{args.backend}/src/athena').resolve()
assert old_binary.is_relative_to(W) and new_binary.is_relative_to(W)
assert old_binary != new_binary
assert sha(old_binary) == old_manifest['sha256'], 'old accepted binary changed'
binaries = {'accepted': old_binary, 'release': new_binary}
hashes = {key: sha(path) for key, path in binaries.items()}
required = {'MPICH_GPU_SUPPORT_ENABLED': '1' if args.backend == 'hip' else '0'}
if args.backend == 'hip':
    required.update(MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED='0',
        MPICH_OFI_NIC_POLICY='GPU', MPICH_GPU_IPC_CACHE_MAX_SIZE='1000',
        HSA_XNACK='1', MPICH_OFI_NUM_CQ_ENTRIES='131072',
        FI_MR_CACHE_MONITOR='kdreg2', FI_CXI_RX_MATCH_MODE='software')
for key, value in required.items():
    assert os.environ.get(key) == value, (key, os.environ.get(key), value)

out.mkdir(parents=True)
(out / 'inputs').mkdir()
record = {'purpose': 'Final-state equivalence of accepted fused and rebuilt release binaries',
    'not_claimed': ['new face-array coverage', 'new performance or speedup measurement'],
    'backend': args.backend, 'slurm_job': os.environ['SLURM_JOB_ID'],
    'binaries': {key: {'path': str(path), 'sha256': hashes[key]} for key, path in binaries.items()},
    'accepted_binary_manifest': str(old_manifest_path),
    'accepted_binary_manifest_sha256': sha(old_manifest_path),
    'matrix_source': str(prior / 'results.json'), 'runner_sha256': sha(Path(__file__)),
    'comparison_rule': 'Exact field/history bytes and all normalized restart bytes, including live forcing state and every diagnostic. No roundoff allowance.',
    'restart_rule': 'Always normalize only 36 known unused root-index bytes. Forced cases independently opt in to validated dormant startup RNG members and four ABI-padding bytes; record all ranges and raw hashes.',
    'required_environment': required,
    'runtime_environment': {key: value for key, value in os.environ.items()
        if key.startswith(('MPICH_', 'HSA_', 'FI_', 'ROCR_', 'HIP_', 'KOKKOS_'))
        or key in ('OMP_NUM_THREADS', 'LD_LIBRARY_PATH')},
    'expected_cases': len(cases), 'cases': {}, 'failures': []}


def save():
    (out / 'results.json').write_text(json.dumps(record, indent=2) + '\n')


def launch(key, binary, fixture, ranks, work):
    assert sha(binary) == hashes[key], f'{key} binary mutated during run'
    work.mkdir(parents=True)
    (work / 'tmp').mkdir()
    env = os.environ.copy()
    env['TMPDIR'] = str(work / 'tmp')
    for name in list(env):
        if name.startswith('ATHENAK_CGL_LF_P1_') or name in (
                'ATHENAK_CGL_LF_ARITHMETIC', 'ATHENAK_CGL_LF_DIAGNOSTICS',
                'ATHENAK_CGL_LF_STS_FLUX', 'ATHENAK_CGL_LF_PROFILE',
                'ATHENAK_CGL_LF_PROFILE_DETAIL', 'ATHENAK_CGL_LF_TASK_TRACE',
                'KOKKOS_TOOLS_LIBS', 'KOKKOS_PROFILE_LIBRARY'):
            env.pop(name, None)
    cmd = ['srun', '--exact', '-N1', '-n', str(ranks), '--ntasks-per-node', str(ranks),
           '--threads-per-core=1', '--cpu-bind=threads']
    if args.backend == 'hip':
        cmd += ['-c7', '--gpus-per-task=1', '--gpu-bind=closest']
    else:
        cmd += ['-c1', '--gpus-per-node=0', '--overlap']
    if args.node:
        cmd += ['-w', args.node]
    cmd += [str(binary), '-i', str(fixture)]
    try:
        proc = subprocess.run(cmd, cwd=work, env=env, capture_output=True,
                              text=True, timeout=900)
        text, code = proc.stdout + proc.stderr, proc.returncode
    except subprocess.TimeoutExpired as exc:
        text = ''.join(x.decode(errors='replace') if isinstance(x, bytes) else x or ''
                       for x in (exc.stdout, exc.stderr))
        text += '\nRELEASE EQUIVALENCE HARNESS: application timeout\n'
        code = 124
    (work / 'stdout.log').write_text(text)
    row = {'command': cmd, 'cwd': str(work), 'returncode': code,
           'binary_sha256': hashes[key], 'input_sha256': sha(fixture)}
    if code == 0:
        history = hst(str(next(work.glob('*.mhd.hst'))))
        row['final_time'] = float(history['time'][-1])
        row['lf_cell_stages'] = float(history['lf_nstage'][-1])
        row['repair_counters'] = {name: float(history[name][-1]) for name in
            ('lf_dfloor', 'lf_pfloor', 'lf_nonfin', 'lf_nonpos', 'lf_hardbd', 'lf_hwproj')}
    (work / 'command.json').write_text(json.dumps(row, indent=2) + '\n')
    return row


def compare(left, right, forced):
    inventory = lambda directory: {str(p.relative_to(directory)): p
        for p in directory.rglob('*') if p.is_file() and p.suffix in ('.rst', '.hst', '.bin')}
    a, b = inventory(left), inventory(right)
    result = {'failures': [], 'files': {}}
    if not a or a.keys() != b.keys():
        result['failures'].append('missing/different output inventory')
    for name in sorted(a.keys() & b.keys()):
        row = {'accepted_raw_sha256': sha(a[name]), 'release_raw_sha256': sha(b[name])}
        if name.endswith('.rst'):
            options = dict(normalize_dormant_startup_forcing=forced,
                           normalize_forcing_padding=forced)
            x, y = restart_hashes(a[name], **options), restart_hashes(b[name], **options)
            row.update(accepted_restart=x, release_restart=y,
                       equal=x['normalized_sha256'] == y['normalized_sha256'])
        else:
            row['equal'] = row['accepted_raw_sha256'] == row['release_raw_sha256']
            if name.endswith('.hst') and not row['equal']:
                x, y = hst(str(a[name])), hst(str(b[name]))
                row['columns'] = {}
                for key in x.keys() & y.keys():
                    same_shape = x[key].shape == y[key].shape
                    row['columns'][key] = {'same_shape': same_shape,
                        'equal_values': bool(same_shape and np.array_equal(x[key], y[key])),
                        'max_abs': float(np.max(np.abs(x[key] - y[key]))) if same_shape else None}
        result['files'][name] = row
        if not row['equal']:
            result['failures'].append(name)
    return result


save()
for case in cases:
    name = case['name']
    fixture = out / 'inputs' / (name + '.athinput')
    shutil.copy2(case['input'], fixture)
    assert sha(fixture) == case['input_sha256']
    entry = {'ranks': case['ranks'], 'forced': case['forced'],
             'source_input': str(case['input']), 'input_sha256': sha(fixture), 'runs': {}}
    record['cases'][name] = entry
    for key, binary in binaries.items():
        entry['runs'][key] = launch(key, binary, fixture, case['ranks'], out / name / key)
        save()
    if all(row['returncode'] == 0 for row in entry['runs'].values()):
        entry['comparison'] = compare(out / name / 'accepted', out / name / 'release', case['forced'])
        if entry['comparison']['failures']:
            record['failures'].append(name)
    else:
        record['failures'].append(name)
    save()
    print(name, 'failures', len(record['failures']), flush=True)
assert all(sha(path) == hashes[key] for key, path in binaries.items())
record['binary_hashes_unchanged'] = True
record['completed_cases'] = len(record['cases'])
save()
print('Final-state equivalence failures:', len(record['failures']), flush=True)
raise SystemExit(bool(record['failures']))
