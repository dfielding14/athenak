"""Regression tests for star particle initialization, formation, and accretion."""

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


def test_star_particle_file_initialization():
    """File initialization should keep the requested particle count and mass."""
    for path in Path(".").glob("star_particle_file.*.hst"):
        path.unlink()
    assert testutils.run("inputs/particles/star_particle_file.athinput")
    row = _last_history_row("star_particle_file.part.hst")
    assert int(round(row[2])) == 2
    assert math.isclose(row[3], 3.0, rel_tol=0.0, abs_tol=1.0e-12)
    assert math.isclose(row[4], 0.0, rel_tol=0.0, abs_tol=1.0e-12)
    assert math.isclose(row[5], 0.0, rel_tol=0.0, abs_tol=1.0e-12)


def test_star_particle_formation_and_accretion():
    """A dense cell should create one particle, then grow it by accretion."""
    for path in Path(".").glob("star_particle_formation.*.hst"):
        path.unlink()
    assert testutils.run("inputs/particles/star_particle_formation.athinput")
    row = _last_history_row("star_particle_formation.part.hst")
    assert int(round(row[2])) == 1
    assert row[4] > 0.0
    assert row[5] > 0.0
    assert math.isclose(row[3], row[4] + row[5], rel_tol=0.0, abs_tol=1.0e-12)


def test_star_particle_outflow_removal():
    """A star particle leaving an outflow domain should be removed cleanly."""
    for path in Path(".").glob("star_particle_outflow.*.hst"):
        path.unlink()
    assert testutils.run("inputs/particles/star_particle_outflow.athinput")
    row = _last_history_row("star_particle_outflow.part.hst")
    assert int(round(row[2])) == 0
    assert math.isclose(row[3], 0.0, rel_tol=0.0, abs_tol=1.0e-12)


def test_star_particle_exact_outer_boundary_removal():
    """A star particle landing exactly on xmax should be removed cleanly."""
    for path in Path(".").glob("star_particle_outflow_exact.*.hst"):
        path.unlink()
    assert testutils.run(
        "inputs/particles/star_particle_outflow.athinput",
        [
            "job/basename=star_particle_outflow_exact",
            "particles/particle_file=inputs/particles/star_particles_outflow_exact.tbl",
        ],
    )
    row = _last_history_row("star_particle_outflow_exact.part.hst")
    assert int(round(row[2])) == 0
    assert math.isclose(row[3], 0.0, rel_tol=0.0, abs_tol=1.0e-12)


def test_star_particle_rejects_multi_meshblock_crossing():
    """A particle timestep that skips MeshBlocks should fail clearly."""
    particle_file = Path("star_particle_multi_block_crossing.tbl")
    particle_file.write_text(
        "# x y z vx vy vz mass t_create\n"
        "0.125 0.5 0.5 20.0 0.0 0.0 1.0 0.0\n"
    )
    try:
        assert not testutils.run_command(
            [
                "./athena",
                "-i",
                "inputs/particles/star_particle_outflow.athinput",
                "job/basename=star_particle_multi_block_crossing",
                "mesh/nx1=16",
                "meshblock/nx1=4",
                "problem/pressure=0.0001",
                f"particles/particle_file={particle_file}",
            ]
        )
    finally:
        particle_file.unlink(missing_ok=True)
