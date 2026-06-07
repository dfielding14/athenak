"""MPI and AMR regression tests for second-moment Ito mass-flux tracers."""

from pathlib import Path
import shutil
import sys

import numpy as np

import test_suite.testutils as testutils
from test_suite.particles.ito_restart_test_utils import make_legacy_restart


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
from read_prtcl_thermo_history import read_history  # noqa: E402


UNIFORM_INPUT = str(ROOT / "inputs/particles/ito_tracers.athinput")
AMR_INPUT = str(ROOT / "inputs/particles/ito_tracers_amr.athinput")
RUN_UNIFORM = Path("run_particles_ito_mpi")
RUN_UNIFORM_FOUR = Path("run_particles_ito_mpi_four")
RUN_UNIFORM_SERIAL = Path("run_particles_ito_mpi_serial_reference")
RUN_AMR = Path("run_particles_ito_amr_mpi")
RUN_AMR_UNINTERRUPTED = Path("run_particles_ito_amr_mpi_uninterrupted")
RUN_LEGACY_V2 = Path("run_particles_ito_legacy_v2_mpi")
RUN_LEGACY_V2_UNINTERRUPTED = Path(
    "run_particles_ito_legacy_v2_mpi_uninterrupted"
)
HIGH_TAG = 2**63 + 12345


def _final_state(path):
    data = read_history(path)
    final = data["cycle"] == np.max(data["cycle"])
    final_indices = np.flatnonzero(final)
    _, unique = np.unique(data["tag"][final], return_index=True)
    final_indices = final_indices[unique]
    order = final_indices[np.argsort(data["tag"][final_indices])]
    tags = data["tag"][order]
    state = np.column_stack(
        [data[name][order] for name in ("gid", "x1", "x2", "x3")]
    )
    return tags, state


def test_ito2_mpi_and_amr():
    """Ito-2 migrates across ranks and preserves continuous positions through AMR."""
    shutil.rmtree(RUN_UNIFORM, ignore_errors=True)
    shutil.rmtree(RUN_UNIFORM_FOUR, ignore_errors=True)
    shutil.rmtree(RUN_UNIFORM_SERIAL, ignore_errors=True)
    shutil.rmtree(RUN_AMR, ignore_errors=True)
    shutil.rmtree(RUN_AMR_UNINTERRUPTED, ignore_errors=True)
    try:
        flags = [
            "meshblock/nx1=8",
            "particles/next_tracer_tag=" + str(HIGH_TAG),
            "particles/ito_covariance_model=full_finite_step",
        ]
        assert testutils.run(
            UNIFORM_INPUT, ["-d", str(RUN_UNIFORM_SERIAL), *flags]
        )
        assert testutils.mpi_run(
            UNIFORM_INPUT, ["-d", str(RUN_UNIFORM), *flags], threads=2
        )
        assert testutils.mpi_run(
            UNIFORM_INPUT,
            ["-d", str(RUN_UNIFORM_FOUR), *flags],
            threads=4,
        )
        assert testutils.mpi_run(
            AMR_INPUT, ["-d", str(RUN_AMR), *flags], threads=2
        )

        serial_history = RUN_UNIFORM_SERIAL / (
            "prtcl_thermo_history/ito_tracers.prtcl_thermo_history.thp"
        )
        mpi_history = RUN_UNIFORM / (
            "prtcl_thermo_history/ito_tracers.prtcl_thermo_history.thp"
        )
        serial_tags, serial_state = _final_state(serial_history)
        mpi_tags, mpi_state = _final_state(mpi_history)
        four_tags, four_state = _final_state(
            RUN_UNIFORM_FOUR
            / "prtcl_thermo_history/ito_tracers.prtcl_thermo_history.thp"
        )
        assert np.min(mpi_tags) == HIGH_TAG
        np.testing.assert_array_equal(mpi_tags, serial_tags)
        np.testing.assert_array_equal(mpi_state, serial_state)
        np.testing.assert_array_equal(four_tags, serial_tags)
        np.testing.assert_array_equal(four_state, serial_state)

        history = RUN_AMR / (
            "prtcl_thermo_history/ito_tracers_amr.prtcl_thermo_history.thp"
        )
        assert history.exists()
        assert history.stat().st_size > 0
        data = read_history(history)
        moved = data["cycle"] > 0
        assert np.any(moved)
        assert np.min(data["tag"]) == HIGH_TAG
        assert np.unique(data["tag"][data["cycle"] == 0]).size == 128

        # A continuous Ito trajectory must not be snapped back to AMR cell centers.
        distance_to_center = np.full(np.count_nonzero(moved), np.inf)
        for dx in (1.0 / 32.0, 1.0 / 64.0):
            cell_coordinate = (data["x1"][moved] + 0.5) / dx - 0.5
            distance_to_center = np.minimum(
                distance_to_center, np.abs(cell_coordinate - np.rint(cell_coordinate))
            )
        assert np.any(distance_to_center > 1.0e-8)

        restart_files = sorted((RUN_AMR / "rst").glob("rank_*/*.rst"))
        assert restart_files
        sizes = np.array([path.stat().st_size for path in restart_files])
        assert np.all(sizes > 0)
        latest_rank0_restart = sorted(
            (RUN_AMR / "rst/rank_00000000").glob("*.rst")
        )[-1]
        assert testutils.run_command(
            [
                "mpirun",
                "-np",
                "2",
                "./athena",
                "-r",
                str(latest_rank0_restart),
                "-d",
                str(RUN_AMR / "restart"),
                "time/nlim=6",
                "time/tlim=0.035",
            ]
        )
        assert testutils.mpi_run(
            AMR_INPUT,
            [
                "-d",
                str(RUN_AMR_UNINTERRUPTED),
                *flags,
                "time/nlim=6",
                "time/tlim=0.035",
            ],
            threads=2,
        )
        restarted_tags, restarted_state = _final_state(
            RUN_AMR
            / "restart/prtcl_thermo_history/ito_tracers_amr.prtcl_thermo_history.thp"
        )
        continuous_tags, continuous_state = _final_state(
            RUN_AMR_UNINTERRUPTED
            / "prtcl_thermo_history/ito_tracers_amr.prtcl_thermo_history.thp"
        )
        np.testing.assert_array_equal(restarted_tags, continuous_tags)
        np.testing.assert_array_equal(restarted_state, continuous_state)
    finally:
        shutil.rmtree(RUN_UNIFORM, ignore_errors=True)
        shutil.rmtree(RUN_UNIFORM_FOUR, ignore_errors=True)
        shutil.rmtree(RUN_UNIFORM_SERIAL, ignore_errors=True)
        shutil.rmtree(RUN_AMR, ignore_errors=True)
        shutil.rmtree(RUN_AMR_UNINTERRUPTED, ignore_errors=True)


def test_bare_legacy_v2_restart_keeps_mpi_message_extents_consistent():
    """A v2 restart infers diagonal mode without changing allocated MPI extents."""
    shutil.rmtree(RUN_LEGACY_V2, ignore_errors=True)
    shutil.rmtree(RUN_LEGACY_V2_UNINTERRUPTED, ignore_errors=True)
    try:
        source = RUN_LEGACY_V2 / "source"
        resumed = RUN_LEGACY_V2 / "resumed"
        source.mkdir(parents=True)
        flags = [
            "meshblock/nx1=8",
            "particles/ito_covariance_model=published_diagonal",
        ]
        assert testutils.mpi_run(
            UNIFORM_INPUT,
            ["-d", str(source), *flags, "time/nlim=1"],
            threads=2,
        )
        latest_rank0 = sorted(
            (source / "rst/rank_00000000").glob("*.rst")
        )[-1]
        rank_files = sorted(
            (source / "rst").glob(f"rank_*/{latest_rank0.name}")
        )
        assert len(rank_files) == 2
        for restart in rank_files:
            make_legacy_restart(restart, restart, 2)

        assert testutils.run_command(
            [
                "mpirun",
                "-np",
                "2",
                "./athena",
                "-r",
                str(latest_rank0),
                "-d",
                str(resumed),
                "time/nlim=2",
            ]
        )
        assert testutils.mpi_run(
            UNIFORM_INPUT,
            [
                "-d",
                str(RUN_LEGACY_V2_UNINTERRUPTED),
                *flags,
                "time/nlim=2",
            ],
            threads=2,
        )
        resumed_tags, resumed_state = _final_state(
            resumed / "prtcl_thermo_history/ito_tracers.prtcl_thermo_history.thp"
        )
        continuous_tags, continuous_state = _final_state(
            RUN_LEGACY_V2_UNINTERRUPTED
            / "prtcl_thermo_history/ito_tracers.prtcl_thermo_history.thp"
        )
        np.testing.assert_array_equal(resumed_tags, continuous_tags)
        np.testing.assert_array_equal(resumed_state, continuous_state)
    finally:
        shutil.rmtree(RUN_LEGACY_V2, ignore_errors=True)
        shutil.rmtree(RUN_LEGACY_V2_UNINTERRUPTED, ignore_errors=True)
