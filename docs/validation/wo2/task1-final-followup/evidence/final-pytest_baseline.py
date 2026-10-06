#!/usr/bin/env python3
"""Stage the repository's test working-directory conventions under WO2."""
import argparse
import os
from pathlib import Path
import sys
import shlex

parser = argparse.ArgumentParser()
parser.add_argument('backend', choices=['hip', 'cpu'])
parser.add_argument('suite', choices=['cpu', 'gpu', 'mpi', 'mpi-gpu'])
args = parser.parse_args()
root = Path('/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2')
repo = Path(os.environ.get('WO2_TEST_SOURCE', '/autofs/nccs-svm1_home2/dfielding/athenak-cgl-wo2-baseline'))
label = os.environ.get('WO2_TEST_LABEL', args.backend)
layout = root / 'baseline' / f'{label}-{args.suite}'
run = layout / 'tst/build/src'
run.mkdir(parents=True, exist_ok=True)
for dest, source in (
    (layout / 'inputs', repo / 'inputs'),
    (layout / 'scripts', repo / 'scripts'),
    (layout / 'docs', repo / 'docs'),
    (layout / 'vis', repo / 'vis'),
    (layout / 'tst/inputs', repo / 'tst/inputs'),
    (layout / 'tst/test_suite', repo / 'tst/test_suite'),
    (run / 'inputs', repo / 'tst/inputs'),
    (run / 'athena', root / 'scripts/launch.py'),
):
    if not dest.exists():
        dest.symlink_to(source)
os.environ['WO2_BACKEND'] = args.backend
os.environ['WO2_BINARY'] = os.environ.get('WO2_TEST_BINARY', str(root / 'baseline/bin' / f'athena-{args.backend}'))
os.environ['TMPDIR'] = str(layout / 'tmp')
Path(os.environ['TMPDIR']).mkdir(exist_ok=True)
os.chdir(layout / 'tst')
sys.path.insert(0, str(repo / 'tst'))
sys.path.insert(0, str(repo / 'vis/python'))
import test_suite.testutils  # noqa: E402
import pytest  # noqa: E402
os.chdir(run)
suite = repo / 'tst/test_suite/cgl'
selection = {
    'cpu': [str(suite), '-k', os.environ.get('WO2_PYTEST_K', '_cpu')],
    'gpu': [str(suite / 'test_cgl_amr_gpu.py')],
    'mpi': [str(suite), '-k', '_mpicpu'],
    'mpi-gpu': [str(suite / 'test_cgl_amr_mpi_gpu.py')],
}[args.suite]
raise SystemExit(pytest.main([
    *selection, '-q', '-ra', '--tb=short',
    '--basetemp', str(layout / 'pytest-tmp'),
    '-o', f'cache_dir={layout / "pytest-cache"}',
    f'--junitxml={layout / "results.xml"}',
    *shlex.split(os.environ.get('WO2_PYTEST_EXTRA', '')),
]))
