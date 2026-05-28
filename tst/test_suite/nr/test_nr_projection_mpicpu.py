"""MPI regression tests for reduced projection output products."""

from pathlib import Path
import os
import shutil
import subprocess

import numpy as np
import projection
import test_suite.testutils as testutils


def _remove_projection_outputs():
    path = Path("proj")
    if path.exists():
        shutil.rmtree(path)


def _run_projection_stats(input_file):
    env = os.environ.copy()
    env["ATHENAK_OUTPUT_IO_STATS"] = "1"
    command = [
        "mpirun",
        "-np",
        "2",
        "./athena",
        "-i",
        f"inputs/{input_file}",
        "output1/single_file_per_node=true",
        "output2/single_file_per_node=true",
    ]
    result = subprocess.run(command, env=env, text=True, capture_output=True)
    assert result.returncode == 0, result.stdout + result.stderr
    stats = {}
    for line in (result.stdout + result.stderr).splitlines():
        if not line.startswith("[output-io] type=proj"):
            continue
        fields = {}
        for token in line.split():
            if "=" in token:
                key, value = token.split("=", 1)
                fields[key] = value
        stats[fields["id"]] = fields
    assert {"full_depth_statistics", "bounded_slab_statistics"} <= stats.keys()
    return stats


def test_projection_mpi_reduces_multi_axis_and_stddev_products():
    try:
        flags = [f"output{index}/single_file_per_node=true" for index in range(1, 10)]
        assert testutils.mpi_run("inputs/projection.athinput", flags=flags, threads=2)
        profile = projection.read_projection(
            "proj/node_00000000/projection.bounded_profile.00000.proj"
        )
        stddev = projection.read_projection(
            "proj/node_00000000/projection.volume_stddev.00000.proj"
        )

        z = profile["x"]
        expected_z = 1.0 + 0.5 * np.sin(2.0 * np.pi * (z + 0.5))
        assert np.allclose(profile["fields"]["dens"][0], 0.05 * expected_z,
                           atol=1.0e-13)
        assert np.allclose(stddev["fields"]["dens"], np.sqrt(0.125), atol=1.0e-13)
    finally:
        _remove_projection_outputs()
        testutils.cleanup()


def test_projection_mpi_preserves_bounded_profile_after_amr():
    try:
        assert testutils.mpi_run(
            "inputs/projection_amr.athinput",
            flags=[
                "output1/single_file_per_node=true",
                "output2/single_file_per_node=true",
            ],
            threads=2,
        )
        profile = projection.read_projection(
            "proj/node_00000000/projection_amr.adaptive_profile.00001.proj"
        )
        assert profile["fields"]["dens"].shape == (1, 16)
        assert np.allclose(profile["fields"]["dens"], 0.05, atol=1.0e-13)
    finally:
        _remove_projection_outputs()
        testutils.cleanup()


def test_native_amr_projection_reports_scaling_metrics_for_refinement_and_bounds():
    try:
        central = _run_projection_stats("projection_scale_central.athinput")
        _remove_projection_outputs()
        broad = _run_projection_stats("projection_scale_broad.athinput")

        central_full = central["full_depth_statistics"]
        central_slab = central["bounded_slab_statistics"]
        broad_full = broad["full_depth_statistics"]
        for metrics in (central_full, central_slab, broad_full):
            assert metrics["layout"] == "native_amr"
            assert metrics["mode"] == "node"
            assert int(metrics["shard_patches"]) > 0
            assert int(metrics["shard_payload_bytes"]) > 0
            assert int(metrics["peak_output_bytes"]) > 0
            assert int(metrics["communicated_bytes"]) > 0
            assert float(metrics["elapsed"]) >= 0.0

        assert int(central_slab["shard_staged_cells"]) < int(
            central_full["shard_staged_cells"]
        )
        assert int(central_slab["shard_payload_bytes"]) < int(
            central_full["shard_payload_bytes"]
        )
        assert int(broad_full["shard_patches"]) > int(central_full["shard_patches"])
        assert int(broad_full["shard_payload_bytes"]) > int(
            central_full["shard_payload_bytes"]
        )
    finally:
        _remove_projection_outputs()
        testutils.cleanup()
