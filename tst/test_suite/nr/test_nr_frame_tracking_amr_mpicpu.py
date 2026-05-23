"""MPI/AMR regression for frame-tracker global moment sampling."""

from pathlib import Path

import athena_read
import numpy as np
import test_suite.testutils as testutils


INPUT = "inputs/frame_tracking_amr.athinput"


def remove_outputs(basename: str) -> None:
    for path in Path(".").glob(f"{basename}.*.hst"):
        path.unlink()


def test_mpi_amr_agrees_with_serial() -> None:
    serial = "FrameTrackAMRSerial"
    parallel = "FrameTrackAMRMPI"
    try:
        assert testutils.run(INPUT, [f"job/basename={serial}"])
        assert testutils.mpi_run(INPUT, [f"job/basename={parallel}"], threads=2)
        serial_hist = athena_read.hst(f"{serial}.frame_tracker.hst")
        parallel_hist = athena_read.hst(f"{parallel}.frame_tracker.hst")
        for field in ("ft_vf_x1", "ft_vf_x2", "ft_weight", "ft_misses"):
            np.testing.assert_allclose(
                parallel_hist[field][-1],
                serial_hist[field][-1],
                rtol=1.0e-12,
                atol=1.0e-12,
            )
    finally:
        remove_outputs(serial)
        remove_outputs(parallel)
        testutils.cleanup()
