"""MPI regression tests for star-particle gravity."""

import math
from pathlib import Path

import test_suite.testutils as testutils


def _last_history_row(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            rows.append([float(value) for value in line.split()])
    assert rows, f"No data rows found in {path}"
    return rows[-1]


def _clean_outputs(prefix):
    for path in Path(".").glob(f"{prefix}.*.hst"):
        path.unlink()
    for path in Path("rst").glob(f"{prefix}.*.rst"):
        path.unlink()
    for path in Path("rst_prtcl").glob(f"{prefix}.*.rst_prtcl"):
        path.unlink()


def _assert_rows_close(row_a, row_b, atol=1.0e-12):
    assert len(row_a) == len(row_b)
    for i, (value_a, value_b) in enumerate(zip(row_a, row_b)):
        assert math.isclose(value_a, value_b, rel_tol=0.0, abs_tol=atol), (
            f"history column {i} differs: {value_a} != {value_b}"
        )


def _two_x1_blocks():
    return [
        "mesh/nx1=8",
        "meshblock/nx1=4",
    ]


def test_star_gravity_binary_direct_mpicpu():
    """Direct self-gravity should agree globally when particles live on separate ranks."""
    _clean_outputs("star_gravity_binary")
    assert testutils.mpi_run(
        "inputs/particles/star_gravity_binary.athinput",
        _two_x1_blocks(),
        threads=2,
    )
    row = _last_history_row("star_gravity_binary.part.hst")
    assert int(round(row[2])) == 2
    assert math.isclose(row[3], 2.0, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(row[10], 0.0, rel_tol=0.0, abs_tol=1.0e-12)
    assert math.isclose(row[11], 0.0, rel_tol=0.0, abs_tol=1.0e-12)
    assert math.isclose(row[16], 0.5, rel_tol=0.0, abs_tol=1.0e-12)
    assert math.isclose(row[20], 0.25, rel_tol=0.0, abs_tol=1.0e-6)


def test_star_gravity_tree_mpicpu():
    """The replicated tree backend should run with particles on separate ranks."""
    _clean_outputs("star_gravity_tree")
    assert testutils.mpi_run(
        "inputs/particles/star_gravity_tree.athinput",
        _two_x1_blocks(),
        threads=2,
    )
    row = _last_history_row("star_gravity_tree.part.hst")
    assert int(round(row[2])) == 2
    assert math.isclose(row[3], 2.0, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(row[10], 0.0, rel_tol=0.0, abs_tol=1.0e-12)
    assert math.isclose(row[11], 0.0, rel_tol=0.0, abs_tol=1.0e-12)
    assert math.isclose(row[16], 0.5, rel_tol=0.0, abs_tol=1.0e-12)


def test_star_gravity_zero_particle_rank_external_mpicpu():
    """A rank with no local particles should still join gravity collectives."""
    _clean_outputs("star_gravity_external")
    assert testutils.mpi_run(
        "inputs/particles/star_gravity_external.athinput",
        _two_x1_blocks(),
        threads=2,
    )
    row = _last_history_row("star_gravity_external.part.hst")
    assert int(round(row[2])) == 1
    assert math.isclose(row[11], -1.0e-2, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(row[17], 0.49995, rel_tol=0.0, abs_tol=1.0e-14)


def test_star_gravity_formation_one_rank_mpicpu():
    """A star formed on one rank should participate in external gravity collectively."""
    _clean_outputs("star_gravity_formation")
    assert testutils.mpi_run(
        "inputs/particles/star_gravity_formation.athinput",
        _two_x1_blocks() + ["problem/x1c=0.6875"],
        threads=2,
    )
    row = _last_history_row("star_gravity_formation.part.hst")
    assert int(round(row[2])) == 1
    assert math.isclose(row[3], 0.02, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(row[4], 0.02, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(row[11], -2.0e-4, rel_tol=0.0, abs_tol=1.0e-14)


def test_star_gravity_accretion_crosses_rank_boundary_mpicpu():
    """A neighboring-cell sink should collect gas from the adjacent MPI rank."""
    _clean_outputs("star_gravity_accretion_neighbor")
    assert testutils.mpi_run(
        "inputs/particles/star_gravity_accretion.athinput",
        _two_x1_blocks()
        + [
            "job/basename=star_gravity_accretion_neighbor",
            "particles/particle_file="
            "inputs/particles/star_gravity_accretion_neighbor.tbl",
            "particles/star_accretion_radius_cells=1",
            "problem/x1c=0.4375",
        ],
        threads=2,
    )
    row = _last_history_row("star_gravity_accretion_neighbor.part.hst")
    expected_gain = 0.000703125
    assert math.isclose(row[5], expected_gain, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(row[3], 1.0 + expected_gain, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(row[10], expected_gain, rel_tol=0.0, abs_tol=1.0e-14)


def test_star_gravity_sidecar_restart_repartitions_mpicpu():
    """Particle sidecars should reload by geometry when the MPI rank count changes."""
    full = "star_gravity_restart_mpi_full"
    segment = "star_gravity_restart_mpi_segment"
    joined = "star_gravity_restart_mpi_joined"
    for prefix in (full, segment, joined):
        _clean_outputs(prefix)

    assert testutils.mpi_run(
        "inputs/particles/star_gravity_restart.athinput",
        _two_x1_blocks() + [f"job/basename={full}"],
        threads=2,
    )
    assert testutils.mpi_run(
        "inputs/particles/star_gravity_restart.athinput",
        _two_x1_blocks()
        + [
            f"job/basename={segment}",
            "time/tlim=0.0025",
            "time/nlim=5",
            f"particles/particle_restart_file=rst_prtcl/{segment}.00001.rst_prtcl",
        ],
        threads=2,
    )
    assert testutils.run_command(
        [
            "mpirun",
            "-np",
            "1",
            "./athena",
            "-r",
            f"rst/{segment}.00001.rst",
            f"job/basename={joined}",
            "time/tlim=0.005",
            "time/nlim=10",
            "particles/star_init=restart",
            f"particles/particle_restart_file=rst_prtcl/{segment}.00001.rst_prtcl",
        ]
    )

    _assert_rows_close(
        _last_history_row(f"{full}.part.hst"),
        _last_history_row(f"{joined}.part.hst"),
    )
