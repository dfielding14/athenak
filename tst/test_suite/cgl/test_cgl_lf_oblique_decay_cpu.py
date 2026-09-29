"""Oblique LF decay across blocks, ranks, and refinement interfaces."""

from pathlib import Path
import numpy as np
import pytest
import test_suite.testutils as testutils

INPUT = "../../../inputs/unit_tests/cgl_lf_oblique_decay_2d.athinput"


def run_oblique_decay(axis, block_size, nranks=0):
    basename = f"cgl_oblique_{axis}_{block_size}_{nranks}"
    chi = np.sqrt(8 / np.pi) / (2 * np.pi)
    b_parallel = np.cos(np.pi / 6) if axis == "x" else np.sin(np.pi / 6)
    duration = 1 / (chi * (2 * np.pi * b_parallel) ** 2)
    flags = [f"job/basename={basename}", f"problem/rotated_axis={axis}",
             f"meshblock/nx1={block_size}", f"meshblock/nx2={block_size}",
             f"time/tlim={duration:.17g}"]
    try:
        if nranks:
            testutils.mpi_run(INPUT, flags, threads=nranks)
        else:
            testutils.run(INPUT, flags)
        history = testutils.athena_read.hst(f"{basename}.mhd.hst")
        assert len(history["time"]) > 10
        for name in ("lf_dfloor", "lf_pfloor", "lf_nonfin", "lf_nonpos", "lf_hardbd"):
            assert history[name][-1] == 0
        metrics = dict(np.loadtxt(f"{basename}.rotated_decay.csv", dtype=str,
                                  delimiter=",", skiprows=1))
        assert np.isclose(float(metrics["expected_amp"]), 1.0e-4 / np.e,
                          rtol=1.0e-13, atol=0)
        assert float(metrics["rel_err"]) < 5.0e-3
        table = np.loadtxt(sorted(Path("tab").glob(f"{basename}.{axis}slice.*.tab"))[-1])
        # Ignore block IDs and local cell indices; compare physical positions/state.
        table = table[np.argsort(table[:, 2]), 2:]
        projection = np.array([float(metrics[k]) for k in ("mean", "sin_amp", "cos_amp")])
        return history, projection, table
    finally:
        for path in Path(".").glob(f"{basename}.*"):
            path.unlink()
        for path in Path("tab").glob(f"{basename}.*"):
            path.unlink()


def assert_oblique_agreement(first, second):
    for key in ("time", "dt", "tot-E", "lf_nstage"):
        np.testing.assert_allclose(first[0][key], second[0][key], rtol=2.0e-13, atol=0)
    # The mean sums 4096 values in a different order across block/rank layouts.
    np.testing.assert_allclose(first[1][0], second[1][0], rtol=0,
                               atol=512*np.finfo(float).eps)
    np.testing.assert_allclose(first[1][1:], second[1][1:], rtol=0,
                               atol=64*np.finfo(float).eps)
    np.testing.assert_allclose(first[2], second[2], rtol=0, atol=64*np.finfo(float).eps)


@pytest.mark.parametrize("axis", ["x", "y"])
def test_cgl_lf_oblique_decay_agrees_across_blocks(axis):
    assert_oblique_agreement(run_oblique_decay(axis, 64), run_oblique_decay(axis, 32))
