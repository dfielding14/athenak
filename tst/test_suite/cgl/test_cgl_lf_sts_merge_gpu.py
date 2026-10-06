"""Optional LF half-sweep transaction application regressions."""
import pytest
from test_suite.cgl import cgl_lf_sts_merge as merge


def test_cgl_lf_sts_merge_default_off_gpu(tmp_path):
    merge.default_off(tmp_path)


@pytest.mark.parametrize("reason", ["collision", "limiter", "strict", "explicit"])
def test_cgl_lf_sts_merge_fallback_gpu(tmp_path, reason):
    merge.fallback(tmp_path, reason)


def test_cgl_lf_sts_merge_barriers_gpu(tmp_path):
    merge.smooth_and_output_barriers(tmp_path)


def test_cgl_lf_sts_merge_restart_gpu(tmp_path):
    merge.restart_sync(tmp_path)


def test_cgl_lf_sts_merge_rejection_gpu(tmp_path):
    merge.dynamic_rejection(tmp_path)


def test_cgl_lf_sts_merge_wall_gpu(tmp_path):
    merge.wall_sync(tmp_path)


def test_cgl_lf_sts_merge_temporal_gpu(tmp_path):
    merge.temporal_refinement(tmp_path)
