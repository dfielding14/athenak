#!/usr/bin/env python3
"""Fixed-input byte comparisons and median CGL-LF performance; scratch only."""
import argparse
import filecmp
import hashlib
import json
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import time

ROOT = Path(__file__).resolve().parent
REPO = Path('/Users/dbf75/.codex/worktrees/bf22/athenak-DF')
STAGED = ROOT / 'inputs'
# These regions execute only in LF/STS work; heat_flux_total already includes
# heat_flux_precompute, directional fluxes, and heat-flux work diagnostics.
STS_BUCKETS = (
    'heat_flux_total', 'sweep_begin_conversion', 'sts_clear_flux',
    'sts_update_copies', 'sts_update_kernel', 'primitive_refresh', 'admissibility',
    'sweep_end_conversion', 'post_sweep_collisions', 'parabolic_init_recv',
)
SHARED_BUCKETS = (
    'parabolic_send_flux', 'parabolic_recv_flux', 'parabolic_restrict_u',
    'parabolic_send_u', 'parabolic_recv_u', 'parabolic_physical_bcs',
    'parabolic_prolongate',
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_input(path):
    blocks = {}
    for line in path.read_text().splitlines():
        line = line.split('#', 1)[0].strip()
        if line.startswith('<'):
            block = line.strip('<>')
            blocks[block] = {}
        elif '=' in line:
            key, value = line.split('=', 1)
            blocks[block][key.strip()] = value.strip()
    return {k: v for k, v in blocks.items() if not k.startswith('output')}


def stage():
    STAGED.mkdir(exist_ok=True)
    g5 = REPO / 'inputs/unit_tests/cgl_lf_oblique_decay_2d.athinput'
    g6 = REPO / 'inputs/unit_tests/cgl_lf_smr_decay_2d.athinput'
    if not g5.exists():
        g5 = Path('/tmp/cgl-wo1-g5/inputs/unit_tests/cgl_lf_oblique_decay_2d.athinput')
    if not g6.exists():
        g6 = Path('/tmp/cgl-wo1-g6/inputs/unit_tests/cgl_lf_smr_decay_2d.athinput')
    cases = [
        ('lf1d', REPO/'inputs/unit_tests/cgl_lf_quant_parallel.athinput', {}, False, True),
        ('lf2d', g5, {'meshblock': {'nx1': 32, 'nx2': 32}}, True, True),
        ('lf3d', REPO/'inputs/unit_tests/cgl_lf_rotated_decay.athinput', {
            'mesh': {'nx1': 32, 'nx2': 32, 'nx3': 32},
            'meshblock': {'nx1': 16, 'nx2': 16, 'nx3': 16},
            'time': {'sts_max_dt_ratio': 20, 'nlim': 10, 'tlim': 1},
            'problem': {'b0': 1, 'by0': 0.3, 'bz0': 0.2}}, True, True),
        ('smr', g6, {}, True, True),
        ('pure_cgl', REPO/'inputs/unit_tests/cgl_pure_paper_oblique_wave.athinput', {}, False, True),
        ('shear', REPO/'tst/inputs/cgl_lf_sts_sbox.athinput', {
            'time': {'nlim': 10, 'tlim': 0.1, 'sts_max_dt_ratio': 20},
            'problem': {'amp': 0.0001, 'bx0': 1}}, True, False),
        ('smr_varied', ROOT/'smr_varied_source.athinput', {}, True, False),
        ('smr_outflow', ROOT/'smr_varied_source.athinput', {
            'mesh': {'ix1_bc': 'outflow', 'ox1_bc': 'outflow'},
            'refined_region1': {'x1min': 0, 'x1max': 0.5},
            'refined_region2': {'x1min': 0.25, 'x1max': 0.5}}, True, False),
        ('outflow', REPO/'inputs/unit_tests/cgl_lf_field_reversal_1d.athinput', {}, False, False),
        ('density', REPO/'inputs/unit_tests/cgl_lf_density_contact.athinput', {}, False, False),
        ('low_b', REPO/'inputs/unit_tests/cgl_lf_quant_parallel.athinput', {
            'time': {'nlim': 10, 'tlim': 1},
            'problem': {'test_mode': 'low_field', 'b0': 5e-11, 'amp': 0.5}}, False, False),
        ('velocity', REPO/'inputs/unit_tests/cgl_lf_field_wave.athinput', {}, False, False),
    ]
    manifest = {}
    for name, source, overrides, mpi, required in cases:
        blocks = read_input(source)
        for block, values in overrides.items():
            blocks.setdefault(block, {}).update(values)
        blocks['job']['basename'] = name
        blocks['time']['ndiag'] = 100000
        blocks['problem']['validation_output'] = 'false'
        if name != 'pure_cgl':
            blocks['mhd']['cgl_lf_profile'] = 'true'
            blocks['mhd']['cgl_lf_profile_detail'] = 'false'
        blocks['output1'] = dict(file_type='hst', data_format='%24.16e', dcycle=1)
        # Initial, every ten cycles, and final snapshots of both representations.
        for i, variable in enumerate(('mhd_u_bcc', 'mhd_w_bcc'), 2):
            blocks[f'output{i}'] = dict(file_type='bin', variable=variable, dcycle=10)
        blocks['output4'] = dict(file_type='rst', dcycle=10)
        if name in ('smr_varied', 'smr_outflow'):
            blocks['output2']['ghost_zones'] = 'true'
            blocks['output3']['ghost_zones'] = 'true'
        text = '# Frozen scratch input for WO1 byte comparison.\n'
        text += '\n'.join(f'\n<{block}>\n' + '\n'.join(f'{k} = {v}' for k, v in values.items())
                           for block, values in blocks.items()) + '\n'
        target = STAGED/f'{name}.athinput'
        target.write_text(text)
        manifest[name] = dict(input=str(target), sha256=sha(target), source=str(source),
                              source_sha256=sha(source), mpi=mpi, required=required)
    (ROOT/'inputs.json').write_text(json.dumps(manifest, indent=2))
    print('Staged', len(manifest), 'fixed inputs at', STAGED)


def metrics(stdout, wall, run):
    cycles = int(re.findall(r'^time=\S+ cycle=(\d+)', stdout, re.M)[-1])
    solver = float(re.findall(r'^cpu time used\s*=\s*(\S+)', stdout, re.M)[-1])
    profile = {}
    for line in stdout.splitlines():
        parts = line.split()
        if len(parts) == 6 and (parts[0] in STS_BUCKETS + SHARED_BUCKETS or
                                parts[0] == 'heat_flux_precompute'):
            profile[parts[0]] = dict(zip(('mean_s', 'max_s', 'mean_calls', 'max_calls',
                                         'max_s_per_call'), map(float, parts[1:])))
    stages = int(profile.get('heat_flux_total', {}).get('mean_calls', 0))
    if stages:
        assert stages == profile['heat_flux_total']['max_calls']
        assert stages == profile['admissibility']['mean_calls']
        geometry = re.search(r'meshblocks_total=(\d+).*meshblock_cells=(\d+)x(\d+)x(\d+)', stdout)
        physical_cells = 1
        for n in geometry.groups():
            physical_cells *= int(n)
        lines = next(run.glob('*.mhd.hst')).read_text().splitlines()
        names = re.findall(r'\[\d+\]=(\S+)', next(x for x in lines if '[1]=' in x))
        last = list(map(float, next(x for x in reversed(lines) if not x.startswith('#')).split()))
        history = dict(zip(names, last))
        assert history['lf_nstage'] == stages*physical_cells
    return dict(cycles=cycles, stages=stages, wall_s_per_cycle=wall/cycles,
                solver_s_per_cycle=solver/cycles,
                sts_profile_s_per_stage=(sum(profile.get(k, {}).get('mean_s', 0)
                                             for k in STS_BUCKETS)/stages if stages else None),
                shared_transport_s_per_stage=(sum(profile.get(k, {}).get('mean_s', 0)
                                                  for k in SHARED_BUCKETS)/stages if stages else None),
                profile=profile)


def outputs(run):
    return {str(p.relative_to(run)): sha(p) for p in sorted(run.rglob('*'))
            if p.is_file() and p.suffix in ('.bin', '.hst', '.rst')}


def compare(left, right):
    first, second = outputs(left), outputs(right)
    if first.keys() != second.keys():
        return ['file inventory differs']
    return [name for name in first if not filecmp.cmp(left/name, right/name, shallow=False)]


def run(args):
    if args.repeats < 3:
        raise SystemExit('Use at least three repeats for timing medians.')
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', args.label):
        raise SystemExit('Label must be a simple directory name.')
    manifest = json.loads((ROOT/'inputs.json').read_text())
    selected = args.cases or list(manifest)
    label_root = ROOT/'runs'/args.label
    label_root.mkdir(parents=True, exist_ok=False)
    baseline = json.loads((ROOT/'runs'/args.compare/'results.json').read_text()) if args.compare else None
    result = dict(label=args.label, baseline=args.compare, inputs=manifest, repeats=args.repeats,
                  binaries={}, cases={}, failures=[])
    for kind, filename in [('serial', args.serial), ('mpi', args.mpi)]:
        if not filename:
            continue
        binary = Path(filename).resolve()
        frozen = label_root/f'athena-{kind}'
        shutil.copy2(binary, frozen)
        result['binaries'][kind] = dict(source=str(binary), sha256=sha(frozen))
        for name in selected:
            case = manifest[name]
            if sha(Path(case['input'])) != case['sha256']:
                raise SystemExit(f'Staged input changed: {name}')
            if baseline and baseline['inputs'][name]['sha256'] != case['sha256']:
                raise SystemExit(f'Baseline input differs: {name}')
            ranks_list = args.mpi_ranks if kind == 'mpi' and case['mpi'] else [0]
            if kind == 'mpi' and not case['mpi']:
                continue
            for nranks in ranks_list:
                key = f'{name}-{kind}-{nranks}'
                samples = []
                for repeat in range(args.repeats):
                    path = label_root/key/f'r{repeat}'
                    path.mkdir(parents=True)
                    command = [str(frozen), '-i', case['input']]
                    if nranks:
                        command = [args.mpiexec, '-n', str(nranks)] + command
                    start = time.perf_counter()
                    proc = subprocess.run(command, cwd=path, capture_output=True, text=True)
                    wall = time.perf_counter() - start
                    (path/'stdout.log').write_text(proc.stdout + proc.stderr)
                    if proc.returncode:
                        raise RuntimeError(f'{key} repeat {repeat} failed; see {path}/stdout.log')
                    sample = metrics(proc.stdout, wall, path)
                    sample['sha256'] = outputs(path)
                    assert any(x.endswith('.bin') for x in sample['sha256'])
                    assert any(x.endswith('.hst') for x in sample['sha256'])
                    if repeat:
                        differences = compare(path, label_root/key/'r0')
                        if differences:
                            result['failures'].append(dict(case=key, repeat=repeat,
                                                           type='within-run nondeterminism', files=differences))
                    if baseline:
                        reference = ROOT/'runs'/args.compare/key/'r0'
                        differences = compare(path, reference)
                        if differences:
                            result['failures'].append(dict(case=key, repeat=repeat,
                                                           type='baseline mismatch', files=differences))
                    samples.append(sample)
                medians = {k: (statistics.median(s[k] for s in samples)
                               if samples[0][k] is not None else None)
                           for k in ('wall_s_per_cycle', 'solver_s_per_cycle',
                                     'sts_profile_s_per_stage', 'shared_transport_s_per_stage')}
                result['cases'][key] = dict(median=medians, samples=samples)
                if baseline:
                    result['cases'][key]['candidate_over_baseline'] = {
                        k: (v/baseline['cases'][key]['median'][k]
                            if v is not None and baseline['cases'][key]['median'][k] else None)
                        for k, v in medians.items()}
                (label_root/'results.json').write_text(json.dumps(result, indent=2))
                print(key, 'cycles', samples[0]['cycles'], 'global_stages', samples[0]['stages'],
                      'median', medians, flush=True)
    print('Byte comparison failures:', len(result['failures']), flush=True)
    if result['failures']:
        raise SystemExit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    sub.add_parser('stage')
    p = sub.add_parser('run')
    p.add_argument('--label', required=True)
    p.add_argument('--serial')
    p.add_argument('--mpi')
    p.add_argument('--mpiexec', default='/opt/homebrew/bin/mpiexec')
    p.add_argument('--mpi-ranks', type=int, nargs='+', default=[1, 4])
    p.add_argument('--repeats', type=int, default=3)
    p.add_argument('--cases', nargs='+')
    p.add_argument('--compare')
    args = parser.parse_args()
    stage() if args.action == 'stage' else run(args)
