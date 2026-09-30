#!/usr/bin/env python3
"""Summarize frozen baseline/candidate medians without changing run data."""
import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
KEYS = [(name, 'serial-0') for name in ('lf1d', 'lf2d', 'lf3d', 'smr', 'pure_cgl')]
KEYS += [(name, 'mpi-4') for name in ('lf2d', 'lf3d', 'smr')]


def report(label, reference=None):
    run = ROOT/'runs'/label
    candidate = json.loads((run/'results.json').read_text())
    baseline_label = reference or candidate['baseline']
    baseline = json.loads((ROOT/'runs'/baseline_label/'results.json').read_text())
    frozen = json.loads((ROOT/'runs'/candidate['baseline']/'results.json').read_text())
    cumulative = baseline_label != candidate['baseline']
    lines = [f'# {label} versus {baseline_label}', '',
             f'{len(candidate["cases"])} configurations, {candidate["repeats"]} repeats; '
             f'{len(candidate["failures"])} byte-comparison failures.', '',
             '| Case | Solver ms/cycle before | After | After/before | Partial compute µs/stage before | After | After/before |',
             '|---|---:|---:|---:|---:|---:|---:|']
    if cumulative:
        lines[-2] += ' Cumulative solver ratio vs post-G | Cumulative compute ratio vs post-G |'
        lines[-1] += '---:|---:|'
    for name, suffix in KEYS:
        key = f'{name}-{suffix}'
        if key not in candidate['cases']:
            continue
        before = baseline['cases'][key]['median']
        after = candidate['cases'][key]['median']
        values = []
        for metric, scale in [('solver_s_per_cycle', 1e3), ('sts_profile_s_per_stage', 1e6)]:
            a, b = before[metric], after[metric]
            values.extend(['—']*3 if a is None else [f'{a*scale:.6g}', f'{b*scale:.6g}', f'{b/a:.4f}'])
        if cumulative:
            values.extend('—' if after[metric] is None else
                          f'{after[metric]/frozen["cases"][key]["median"][metric]:.4f}'
                          for metric in ('solver_s_per_cycle', 'sts_profile_s_per_stage'))
        lines.append('| '+ ' | '.join([key]+values)+' |')
    lines += ['', 'Partial compute per stage sums existing exclusive LF/STS buckets (rank mean for MPI). Whole-STS wall time is not separately available. Shared transport includes RK/initialization and is excluded here; all samples and timing categories remain in results.json.', '',
              'Baseline JSON: '+str(ROOT/'runs'/baseline_label/'results.json'),
              'Candidate JSON: '+str(run/'results.json'), '']
    if candidate['failures']:
        lines += ['Exact failed files:', '', '```json', json.dumps(candidate['failures'], indent=2), '```', '']
    text = '\n'.join(lines)
    filename = f'comparison-vs-{baseline_label}.md' if reference else 'comparison.md'
    (run/filename).write_text(text)
    print(text)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('label')
    p.add_argument('--reference')
    args = p.parse_args()
    report(args.label, args.reference)
