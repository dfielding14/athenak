"""Opt-in GPU MPI coordination of transactional LF sweeps."""
import shlex
from test_suite.cgl import cgl_lf_sts_merge as merge
from test_suite.cgl import test_cgl_amr_mpi_gpu as amr

pytestmark = amr.pytestmark


def test_cgl_lf_sts_merge_rejection_mpi_gpu(tmp_path):
    assert amr.MPI_GPU_LAUNCHER
    prefix = tuple(shlex.split(amr.MPI_GPU_LAUNCHER))
    merge.dynamic_rejection(tmp_path, tuple(
        prefix + ("-N", "1", "-n", str(ranks), "--ntasks-per-node", str(ranks))
        for ranks in (1, 4)))
