# Regression test for shared, per-rank, and per-node PDF output.

import logging
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import scripts.utils.athena as athena  # noqa: F401

_REPO_DIR = Path(__file__).resolve().parents[3]
sys.path.insert(0, str((_REPO_DIR / 'vis' / 'python').resolve()))
import read_pdf  # noqa

logger = logging.getLogger('athena' + __name__[7:])

_NPROC = 2
_INPUT = 'tests/pdf_partitioned.athinput'
_BUILD_DIR = _REPO_DIR / 'tst' / 'build'
_EXE_DIR = _BUILD_DIR / 'src'
_PDF_DIR = _EXE_DIR / 'pdf_hydro_w_d'
_PDF_ND_DIR = _EXE_DIR / 'pdf_hydro_w_d_hydro_w_vx_hydro_w_vy'
_SKIP_REASON = None


def _mpi_enabled():
    config_path = _BUILD_DIR / 'config.hpp'
    if not config_path.is_file():
        return False
    return '#define MPI_PARALLEL_ENABLED 1' in config_path.read_text()


def _run_athena(extra_args):
    exe_path = str((_EXE_DIR / 'athena').resolve())
    input_path = str((_REPO_DIR / 'inputs' / _INPUT).resolve())
    cmd = ['mpiexec', '-n', str(_NPROC), exe_path, '-i', input_path] + extra_args
    logger.debug('Executing: %s (cwd=%s)', ' '.join(cmd), _EXE_DIR)
    subprocess.check_call(cmd, cwd=str(_EXE_DIR))


def _clean_pdf_dir(path):
    if path.is_dir():
        shutil.rmtree(path)


def _clean_output_dirs():
    _clean_pdf_dir(_PDF_DIR)
    _clean_pdf_dir(_PDF_ND_DIR)


def _latest_data_file(root):
    files = sorted(path for path in root.rglob('*.pdf')
                   if not path.name.endswith('.header.pdf'))
    assert files, f'no PDF outputs found under {root}'
    return files[-1]


def _assert_matching_data(reference, candidate):
    assert candidate['header']['format'] in ('dense', 'sparse_coo')
    assert candidate['header']['weight'] == reference['header']['weight']
    assert candidate['header']['ndim'] == reference['header']['ndim']
    assert candidate['header']['shape'] == reference['header']['shape']
    assert np.isclose(candidate['time'], reference['time'])
    assert np.allclose(candidate['pdf'], reference['pdf'], rtol=1.0e-13, atol=1.0e-14)


def run(**kwargs):
    global _SKIP_REASON
    _SKIP_REASON = None
    logger.debug('Running test %s', __name__)

    if not _mpi_enabled():
        _SKIP_REASON = 'MPI is not enabled in build/config.hpp'
        logger.warning(_SKIP_REASON)
        return
    if shutil.which('mpiexec') is None:
        _SKIP_REASON = 'mpiexec is not available'
        logger.warning(_SKIP_REASON)
        return

    common_args = ['time/nlim=0', 'time/tlim=0.0']

    _clean_output_dirs()
    _run_athena(['job/basename=pdf_shared'] + common_args)
    shared = read_pdf.read_pdf(str(_latest_data_file(_PDF_DIR)))
    nd_shared = read_pdf.read_pdf(str(_latest_data_file(_PDF_ND_DIR)))
    assert shared['header']['format'] == 'dense'
    assert shared['header']['weight'] == 'mass'
    assert nd_shared['header']['format'] == 'dense'
    assert nd_shared['header']['weight'] == 'variable'
    assert nd_shared['header']['ndim'] == 3
    assert nd_shared['header']['dimensions'][1]['scale'] == 'symlog'

    _clean_output_dirs()
    _run_athena(['job/basename=pdf_node'] + common_args
                + ['output1/single_file_per_node=true',
                   'output2/single_file_per_node=true'])
    node = read_pdf.read_pdf(str(_latest_data_file(_PDF_DIR / 'node_00000000')))
    nd_node = read_pdf.read_pdf(str(_latest_data_file(_PDF_ND_DIR / 'node_00000000')))
    assert node['header']['format'] == 'sparse_coo'
    assert nd_node['header']['format'] == 'sparse_coo'
    _assert_matching_data(shared, node)
    _assert_matching_data(nd_shared, nd_node)

    _clean_output_dirs()
    _run_athena(['job/basename=pdf_rank'] + common_args
                + ['output1/single_file_per_rank=true',
                   'output2/single_file_per_rank=true'])
    rank = read_pdf.read_pdf(str(_latest_data_file(_PDF_DIR / 'rank_00000000')))
    nd_rank = read_pdf.read_pdf(str(_latest_data_file(_PDF_ND_DIR / 'rank_00000000')))
    assert rank['header']['format'] == 'sparse_coo'
    assert nd_rank['header']['format'] == 'sparse_coo'
    _assert_matching_data(shared, rank)
    _assert_matching_data(nd_shared, nd_rank)


def analyze():
    if _SKIP_REASON is not None:
        logger.warning('Skipping %s: %s', __name__, _SKIP_REASON)
    return True
