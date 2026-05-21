"""MPI regression tests for star particle accounting and removal."""

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


def _two_block_flags():
    return [
        "mesh/nx1=8",
        "meshblock/nx1=4",
    ]


def test_star_particle_file_initialization_mpicpu():
    """File-loaded particles should reduce to the same total count and mass under MPI."""
    for path in Path(".").glob("star_particle_file.*.hst"):
        path.unlink()
    assert testutils.mpi_run(
        "inputs/particles/star_particle_file.athinput",
        _two_block_flags(),
        threads=2,
    )
    row = _last_history_row("star_particle_file.part.hst")
    assert int(round(row[2])) == 2
    assert math.isclose(row[3], 3.0, rel_tol=0.0, abs_tol=1.0e-12)


def test_star_particle_formation_and_accretion_mpicpu():
    """Formation on one rank should update the global particle count collectively."""
    for path in Path(".").glob("star_particle_formation.*.hst"):
        path.unlink()
    assert testutils.mpi_run(
        "inputs/particles/star_particle_formation.athinput",
        _two_block_flags() + ["problem/x1c=0.6875"],
        threads=2,
    )
    row = _last_history_row("star_particle_formation.part.hst")
    assert int(round(row[2])) == 1
    assert row[4] > 0.0
    assert row[5] > 0.0
    assert math.isclose(row[3], row[4] + row[5], rel_tol=0.0, abs_tol=1.0e-12)


def test_star_particle_outflow_removal_mpicpu():
    """A particle leaving an outflow boundary should be removed from the MPI total."""
    for path in Path(".").glob("star_particle_outflow.*.hst"):
        path.unlink()
    assert testutils.mpi_run(
        "inputs/particles/star_particle_outflow.athinput",
        _two_block_flags(),
        threads=2,
    )
    row = _last_history_row("star_particle_outflow.part.hst")
    assert int(round(row[2])) == 0
    assert math.isclose(row[3], 0.0, rel_tol=0.0, abs_tol=1.0e-12)
