#!/usr/bin/env python3
"""Summarize completed release comparisons without modifying raw evidence."""
from collections import Counter
import hashlib
import json
from pathlib import Path

W = Path('/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2')
F = W / 'final-research'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


summary = {
    'purpose': 'Accepted fused to final release state equivalence',
    'claims': 'Exact field and history bytes, exact normalized restart bytes including diagnostics and live RNG state; no numerical tolerance.',
    'limits': 'Ordinary evolved-run outputs only. No resumed-run validation, new face coverage, or performance claim. Passive fixtures use PLM.',
    'normalization': 'Known 36 unused root-index bytes; forced cases additionally opt in to validated dormant seed1 startup RNG members and native RNG ABI padding. Raw hashes and precise byte ranges remain in result files.',
    'backends': {},
}
for backend in ('cpu', 'hip'):
    path = F / f'release-equivalence-{backend}-ready-5629018/results.json'
    data = json.loads(path.read_text())
    assert data.get('completed_cases') == data['expected_cases'] == 25
    assert data.get('binary_hashes_unchanged') is True
    assert not data['failures']
    files = Counter()
    ranks = Counter()
    normalized_only = []
    cases = []
    for name, case in data['cases'].items():
        assert len(case['runs']) == 2
        assert not case['comparison']['failures']
        ranks[case['ranks']] += 1
        for run in case['runs'].values():
            assert run['returncode'] == 0
            assert all(value == 0 for value in run['repair_counters'].values())
        local = Counter()
        for relative, row in case['comparison']['files'].items():
            assert row['equal']
            suffix = Path(relative).suffix
            files[suffix] += 1
            local[suffix] += 1
            if row['accepted_raw_sha256'] != row['release_raw_sha256']:
                assert suffix == '.rst'
                normalized_only.append({'case': name, 'file': relative,
                    'accepted': row['accepted_restart'], 'release': row['release_restart']})
        assert all(local[suffix] > 0 for suffix in ('.bin', '.hst', '.rst'))
        cases.append({'name': name, 'ranks': case['ranks'], 'forced': case['forced'],
            'input_sha256': case['input_sha256'], 'files': dict(local),
            'final_time': case['runs']['release']['final_time'],
            'lf_cell_stages': case['runs']['release']['lf_cell_stages']})
    summary['backends'][backend] = {
        'results': str(path), 'results_sha256': sha(path),
        'binaries': data['binaries'], 'cases': cases, 'pairs': len(cases),
        'applications': 2 * len(cases), 'rank_counts': dict(ranks),
        'file_pairs': dict(files), 'raw_restart_differences': normalized_only,
        'all_repair_counters_zero': True, 'binary_hashes_unchanged': True,
        'runtime_environment': data['runtime_environment'],
    }
summary['source_manifest'] = {
    'path': str(F / 'release-source-manifest.json'),
    'sha256': sha(F / 'release-source-manifest.json'),
}
summary['binary_manifests'] = {
    backend: {'path': str(F / f'release-bin/manifest-{backend}.json'),
              'sha256': sha(F / f'release-bin/manifest-{backend}.json')}
    for backend in ('cpu', 'hip')
}
summary['retained_failed_attempt'] = {
    'path': str(F / 'release-equivalence-cpu-5629018'),
    'status': 'Stopped after 13 pairs whose accepted applications passed and whose release processes failed at dynamic loading, before solver initialization.',
    'cause': 'The initial CPU link inherited ROCm and required libamdhip64.so.6. Unchanged objects were relinked in the CPU-only environment; the original failed-loader binary and logs are preserved.',
}
target = F / 'release-equivalence-summary.json'
target.write_text(json.dumps(summary, indent=2) + '\n')
print(target)
print(json.dumps({backend: {key: row[key] for key in
    ('pairs', 'applications', 'rank_counts', 'file_pairs')}
    for backend, row in summary['backends'].items()}, indent=2))
