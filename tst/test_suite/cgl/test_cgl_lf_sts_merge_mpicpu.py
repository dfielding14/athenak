"""Every rank must make the same CFL rollback decision."""
from test_suite.cgl import cgl_lf_sts_merge as merge


def test_cgl_lf_sts_merge_rejection_mpicpu(tmp_path):
    merge.dynamic_rejection(tmp_path, (("mpirun", "-np", "1"),
                                      ("mpirun", "-np", "4")))
