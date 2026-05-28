# Regression test for shared, per-rank, and per-node binary outputs.

import logging
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import scripts.utils.athena as athena  # noqa: F401

_REPO_DIR = Path(__file__).resolve().parents[3]
sys.path.insert(0, str((_REPO_DIR / 'vis' / 'python').resolve()))
import bin_convert  # noqa

logger = logging.getLogger('athena' + __name__[7:])

_NPROC = 2
_INPUT = 'tests/binary_partitioned.athinput'
_BUILD_DIR = _REPO_DIR / 'tst' / 'build'
_EXE_DIR = _BUILD_DIR / 'src'
_BIN_DIR = _EXE_DIR / 'bin'
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


def _clean_output_dirs():
    if _BIN_DIR.is_dir():
        shutil.rmtree(_BIN_DIR)
    for path in _EXE_DIR.glob('cbin_*'):
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()


def _latest_file(root, suffix, shard=None):
    files = sorted(root.rglob('*' + suffix))
    if shard is not None:
        files = [path for path in files if shard in path.parts]
    assert files, f'no {suffix} outputs found under {root} for shard={shard}'
    return files[-1]


def _read_outputs(shard=None):
    bin_path = _latest_file(_BIN_DIR, '.bin', shard)
    cbin_path = _latest_file(_EXE_DIR, '.cbin', shard)
    if shard is None:
        return (bin_convert.read_binary(str(bin_path)),
                bin_convert.read_coarsened_binary(str(cbin_path)))
    return (bin_convert.read_all_ranks_binary(str(bin_path)),
            bin_convert.read_all_ranks_coarsened_binary(str(cbin_path)))


def _meshblock_indices(data):
    logical = {}
    for index, location in enumerate(data['mb_logical']):
        key = tuple(int(value) for value in location)
        assert key not in logical, f'duplicate meshblock logical location {key}'
        logical[key] = index
    return logical


def _assert_matching_data(reference, candidate):
    for key in ('cycle', 'n_mbs', 'Nx1', 'Nx2', 'Nx3', 'var_names'):
        assert candidate[key] == reference[key], f'metadata mismatch for {key}'
    assert np.isclose(candidate['time'], reference['time'])

    ref_indices = _meshblock_indices(reference)
    candidate_indices = _meshblock_indices(candidate)
    assert candidate_indices.keys() == ref_indices.keys()
    for location, ref_index in ref_indices.items():
        candidate_index = candidate_indices[location]
        assert np.array_equal(candidate['mb_index'][candidate_index],
                              reference['mb_index'][ref_index])
        assert np.allclose(candidate['mb_geometry'][candidate_index],
                           reference['mb_geometry'][ref_index])
        for variable in reference['var_names']:
            assert np.array_equal(candidate['mb_data'][variable][candidate_index],
                                  reference['mb_data'][variable][ref_index])


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
    _run_athena(['job/basename=binary_shared'] + common_args)
    shared_bin, shared_cbin = _read_outputs()

    _clean_output_dirs()
    _run_athena(['job/basename=binary_node'] + common_args
                + ['output1/single_file_per_node=true',
                   'output2/single_file_per_node=true'])
    node_bin, node_cbin = _read_outputs('node_00000000')
    _assert_matching_data(shared_bin, node_bin)
    _assert_matching_data(shared_cbin, node_cbin)

    _clean_output_dirs()
    _run_athena(['job/basename=binary_rank'] + common_args
                + ['output1/single_file_per_rank=true',
                   'output2/single_file_per_rank=true'])
    rank_bin, rank_cbin = _read_outputs('rank_00000000')
    _assert_matching_data(shared_bin, rank_bin)
    _assert_matching_data(shared_cbin, rank_cbin)

    _clean_output_dirs()
    _run_athena(['job/basename=binary_gid_shared'] + common_args
                + ['output1/gid=7', 'output2/gid=7'])
    gid_shared_bin, gid_shared_cbin = _read_outputs()
    assert gid_shared_bin['n_mbs'] == 1
    assert gid_shared_cbin['n_mbs'] == 1

    _clean_output_dirs()
    _run_athena(['job/basename=binary_gid_rank'] + common_args
                + ['output1/gid=7', 'output2/gid=7',
                   'output1/single_file_per_rank=true',
                   'output2/single_file_per_rank=true'])
    gid_rank_bin, gid_rank_cbin = _read_outputs('rank_00000000')
    _assert_matching_data(gid_shared_bin, gid_rank_bin)
    _assert_matching_data(gid_shared_cbin, gid_rank_cbin)


def analyze():
    if _SKIP_REASON is not None:
        logger.warning('Skipping %s: %s', __name__, _SKIP_REASON)
    return True
