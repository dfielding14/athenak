"""Regression tests for star-particle gravity."""

import math
import struct
from pathlib import Path

import test_suite.testutils as testutils


_RESTART_HEADER = "@16s11i3d128s"
_RESTART_MAGIC = b"ATHKPRTCLRST1"


def _last_history_row(path):
    return _history_rows(path)[-1]


def _history_rows(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            rows.append([float(value) for value in line.split()])
    assert rows, f"No data rows found in {path}"
    return rows


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


def _read_restart_header(path):
    data = Path(path).read_bytes()[:struct.calcsize(_RESTART_HEADER)]
    unpacked = struct.unpack(_RESTART_HEADER, data)
    ints = unpacked[1:12]
    doubles = unpacked[12:15]
    basename = unpacked[15].split(b"\0", 1)[0].decode()
    return {
        "magic": unpacked[0].split(b"\0", 1)[0],
        "version": ints[0],
        "header_size": ints[1],
        "endian_marker": ints[2],
        "real_size": ints[3],
        "int_size": ints[4],
        "ncycle": ints[5],
        "total_particles": ints[6],
        "nrdata": ints[7],
        "nidata": ints[8],
        "max_tag": ints[9],
        "basename_length": ints[10],
        "time": doubles[0],
        "formed": doubles[1],
        "accreted": doubles[2],
        "basename": basename,
    }


def _run_sidecar_restart(inputfile, full_prefix, segment_prefix, joined_prefix,
                         stop_time, final_time, stop_cycle, final_cycle):
    _clean_outputs(full_prefix)
    _clean_outputs(segment_prefix)
    _clean_outputs(joined_prefix)

    assert testutils.run(inputfile, [f"job/basename={full_prefix}"])
    assert testutils.run(
        inputfile,
        [
            f"job/basename={segment_prefix}",
            f"time/tlim={stop_time}",
            f"time/nlim={stop_cycle}",
            f"particles/particle_restart_file=rst_prtcl/{segment_prefix}.00001.rst_prtcl",
        ],
    )
    assert testutils.run_command(
        [
            "./athena",
            "-r",
            f"rst/{segment_prefix}.00001.rst",
            f"job/basename={joined_prefix}",
            f"time/tlim={final_time}",
            f"time/nlim={final_cycle}",
            "particles/star_init=restart",
            f"particles/particle_restart_file=rst_prtcl/{segment_prefix}.00001.rst_prtcl",
        ]
    )


def test_star_gravity_constant_external_acceleration():
    """A single star should follow the configured constant external acceleration."""
    _clean_outputs("star_gravity_external")
    assert testutils.run("inputs/particles/star_gravity_external.athinput")
    row = _last_history_row("star_gravity_external.part.hst")
    assert int(round(row[2])) == 1
    assert math.isclose(row[11], -1.0e-2, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(row[17], 0.49995, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(row[19], 1.0, rel_tol=0.0, abs_tol=1.0e-14)


def test_star_gravity_rk4_constant_external_acceleration():
    """The RK4 pusher should match the constant-acceleration analytic update."""
    _clean_outputs("star_gravity_external_rk4")
    assert testutils.run(
        "inputs/particles/star_gravity_external.athinput",
        [
            "job/basename=star_gravity_external_rk4",
            "star_gravity/integrator=rk4",
        ],
    )
    row = _last_history_row("star_gravity_external_rk4.part.hst")
    assert int(round(row[2])) == 1
    assert math.isclose(row[11], -1.0e-2, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(row[17], 0.49995, rel_tol=0.0, abs_tol=1.0e-14)


def test_star_gravity_rk4_user_acceleration_uses_stage_state():
    """User callbacks should receive RK4 stage positions and times."""
    _clean_outputs("star_gravity_external_user")
    assert testutils.run(
        "inputs/particles/star_gravity_external.athinput",
        [
            "job/basename=star_gravity_external_user",
            "star_gravity/integrator=rk4",
            "star_gravity/external_acceleration=user",
        ],
    )
    row = _last_history_row("star_gravity_external_user.part.hst")
    assert math.isclose(row[10], 0.5*(0.01 + 0.01**3/6.0),
                        rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(row[11], 0.5*0.01**2, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(row[16], 0.5*(1.0 + 0.01**2/2.0 + 0.01**4/24.0),
                        rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(row[17], 0.5 + 0.01**3/6.0, rel_tol=0.0, abs_tol=1.0e-14)


def test_star_gravity_point_mass_external_acceleration():
    """A fixed point-mass potential should accelerate stars toward the source."""
    _clean_outputs("star_gravity_point_mass")
    assert testutils.run(
        "inputs/particles/star_gravity_external.athinput",
        [
            "job/basename=star_gravity_point_mass",
            "star_gravity/external_acceleration=point_mass",
            "star_gravity/external_mass=2.0",
            "star_gravity/external_x=0.75",
            "star_gravity/external_y=0.5",
            "star_gravity/external_z=0.5",
        ],
    )
    row = _last_history_row("star_gravity_point_mass.part.hst")
    assert int(round(row[2])) == 1
    assert row[10] > 0.32
    assert math.isclose(row[16], 0.5016, rel_tol=0.0, abs_tol=1.0e-14)
    assert row[19] > 32.0


def test_star_gravity_direct_binary_conservation():
    """The direct backend should preserve binary center of mass and total momentum."""
    _clean_outputs("star_gravity_binary")
    assert testutils.run("inputs/particles/star_gravity_binary.athinput")
    row = _last_history_row("star_gravity_binary.part.hst")
    assert int(round(row[2])) == 2
    assert math.isclose(row[3], 2.0, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(row[10], 0.0, rel_tol=0.0, abs_tol=1.0e-12)
    assert math.isclose(row[11], 0.0, rel_tol=0.0, abs_tol=1.0e-12)
    assert math.isclose(row[12], 0.0, rel_tol=0.0, abs_tol=1.0e-12)
    assert math.isclose(row[16], 0.5, rel_tol=0.0, abs_tol=1.0e-12)
    assert math.isclose(row[17], 0.5, rel_tol=0.0, abs_tol=1.0e-12)
    assert math.isclose(row[18], 0.5, rel_tol=0.0, abs_tol=1.0e-12)
    assert row[6] == 5.0e-4
    assert row[8] < 0.0
    assert math.isclose(row[20], 0.25, rel_tol=0.0, abs_tol=1.0e-6)


def test_star_gravity_tree_binary_conservation():
    """The replicated tree backend should handle the binary orbit path."""
    _clean_outputs("star_gravity_tree")
    assert testutils.run("inputs/particles/star_gravity_tree.athinput")
    row = _last_history_row("star_gravity_tree.part.hst")
    assert int(round(row[2])) == 2
    assert math.isclose(row[3], 2.0, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(row[10], 0.0, rel_tol=0.0, abs_tol=1.0e-12)
    assert math.isclose(row[11], 0.0, rel_tol=0.0, abs_tol=1.0e-12)
    assert math.isclose(row[16], 0.5, rel_tol=0.0, abs_tol=1.0e-12)
    assert math.isclose(row[20], 0.25, rel_tol=0.0, abs_tol=1.0e-6)


def test_star_gravity_tree_retains_coincident_particles():
    """Tree leaves should retain all particles when multiple stars share a position."""
    overrides = [
        "particles/particle_file=inputs/particles/star_gravity_tree_coincident.tbl",
        "time/nlim=1",
        "time/tlim=0.0005",
    ]
    _clean_outputs("star_gravity_tree_coincident_direct")
    assert testutils.run(
        "inputs/particles/star_gravity_tree.athinput",
        ["job/basename=star_gravity_tree_coincident_direct",
         "star_gravity/force_method=direct"] + overrides,
    )
    _clean_outputs("star_gravity_tree_coincident")
    assert testutils.run(
        "inputs/particles/star_gravity_tree.athinput",
        ["job/basename=star_gravity_tree_coincident"] + overrides,
    )
    _assert_rows_close(
        _last_history_row("star_gravity_tree_coincident_direct.part.hst"),
        _last_history_row("star_gravity_tree_coincident.part.hst"),
    )


def test_star_gravity_tree_can_skip_exact_diagnostics():
    """Large tree runs should be able to avoid the exact pairwise history pass."""
    _clean_outputs("star_gravity_tree_fast_diagnostics")
    assert testutils.run(
        "inputs/particles/star_gravity_tree.athinput",
        [
            "job/basename=star_gravity_tree_fast_diagnostics",
            "star_gravity/exact_diagnostics=false",
        ],
    )
    row = _last_history_row("star_gravity_tree_fast_diagnostics.part.hst")
    assert int(round(row[2])) == 2
    assert math.isclose(row[8], 0.0, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(row[20], 0.0, rel_tol=0.0, abs_tol=1.0e-14)


def test_star_gravity_pair_orbit_requires_exact_diagnostics():
    """Pair-orbit timestep control should reject disabled exact diagnostics."""
    _clean_outputs("star_gravity_tree_bad_pair_orbit")
    assert not testutils.run_command(
        [
            "./athena",
            "-i",
            "inputs/particles/star_gravity_tree.athinput",
            "job/basename=star_gravity_tree_bad_pair_orbit",
            "star_gravity/exact_diagnostics=false",
            "star_gravity/timestep_mode=pair_orbit",
        ]
    )


def test_star_gravity_formation_feels_external_acceleration():
    """A newly formed star should be accelerated during the same particle update."""
    _clean_outputs("star_gravity_formation")
    assert testutils.run("inputs/particles/star_gravity_formation.athinput")
    row = _last_history_row("star_gravity_formation.part.hst")
    assert int(round(row[2])) == 1
    assert math.isclose(row[3], 0.02, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(row[4], 0.02, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(row[11], -2.0e-4, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(row[17], 0.62495, rel_tol=0.0, abs_tol=1.0e-14)


def test_star_gravity_accretion_conserves_particle_momentum_gain():
    """Accreted gas momentum should be added to the star when requested."""
    _clean_outputs("star_gravity_accretion")
    assert testutils.run("inputs/particles/star_gravity_accretion.athinput")
    row = _last_history_row("star_gravity_accretion.part.hst")
    expected_gain = 0.00140625
    assert int(round(row[2])) == 1
    assert math.isclose(row[5], expected_gain, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(row[3], 1.0 + expected_gain, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(row[10], expected_gain, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(row[19], 0.0, rel_tol=0.0, abs_tol=1.0e-14)


def test_star_gravity_accretion_crosses_meshblock_boundary():
    """A neighboring-cell sink should cross a local MeshBlock boundary."""
    _clean_outputs("star_gravity_accretion_neighbor")
    assert testutils.run(
        "inputs/particles/star_gravity_accretion.athinput",
        [
            "job/basename=star_gravity_accretion_neighbor",
            "mesh/nx1=8",
            "meshblock/nx1=4",
            "particles/particle_file="
            "inputs/particles/star_gravity_accretion_neighbor.tbl",
            "particles/star_accretion_radius_cells=1",
            "problem/x1c=0.4375",
        ],
    )
    row = _last_history_row("star_gravity_accretion_neighbor.part.hst")
    expected_gain = 0.000703125
    assert math.isclose(row[5], expected_gain, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(row[3], 1.0 + expected_gain, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(row[10], expected_gain, rel_tol=0.0, abs_tol=1.0e-14)


def test_star_gravity_acceleration_timestep_can_recover():
    """The acceleration timestep should increase again when the force weakens."""
    _clean_outputs("star_gravity_timestep_recovery")
    assert testutils.run(
        "inputs/particles/star_gravity_external.athinput",
        [
            "job/basename=star_gravity_timestep_recovery",
            "time/nlim=3",
            "time/tlim=1.0",
            "particles/dt=0.1",
            "particles/particle_file=inputs/particles/star_gravity_outward.tbl",
            "star_gravity/timestep_mode=acceleration",
            "star_gravity/max_timestep=0.1",
            "star_gravity/external_acceleration=point_mass",
            "star_gravity/external_mass=0.001",
            "star_gravity/external_x=0.25",
            "star_gravity/external_y=0.5",
            "star_gravity/external_z=0.5",
            "star_gravity/external_softening_length=0.01",
        ],
    )
    rows = _history_rows("star_gravity_timestep_recovery.part.hst")
    unique_rows = []
    for row in rows:
        if not unique_rows or row[0] != unique_rows[-1][0]:
            unique_rows.append(row)
    assert len(unique_rows) >= 3
    assert unique_rows[-1][6] > unique_rows[-2][6]


def test_star_gravity_sidecar_restart_binary_equivalence():
    """A sidecar particle restart should reproduce the uninterrupted binary run."""
    full = "star_gravity_restart_full"
    segment = "star_gravity_restart_segment"
    joined = "star_gravity_restart_joined"
    _run_sidecar_restart(
        "inputs/particles/star_gravity_restart.athinput",
        full,
        segment,
        joined,
        "0.0025",
        "0.005",
        "5",
        "10",
    )

    header = _read_restart_header(f"rst_prtcl/{segment}.00001.rst_prtcl")
    assert header["magic"] == _RESTART_MAGIC
    assert header["version"] == 2
    assert header["header_size"] == struct.calcsize(_RESTART_HEADER)
    assert header["endian_marker"] == 0x01020304
    assert header["real_size"] == 8
    assert header["int_size"] == 4
    assert header["ncycle"] == 5
    assert header["total_particles"] == 2
    assert header["nrdata"] == 8
    assert header["nidata"] == 2
    assert header["max_tag"] == 1
    assert header["basename"] == segment
    assert header["basename_length"] == len(segment)
    assert math.isclose(header["time"], 0.0025, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(header["formed"], 0.0, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(header["accreted"], 0.0, rel_tol=0.0, abs_tol=1.0e-14)

    _assert_rows_close(
        _last_history_row(f"{full}.part.hst"),
        _last_history_row(f"{joined}.part.hst"),
    )


def test_star_gravity_sidecar_restart_preserves_formation_history():
    """Restarted formation should preserve counters and continue tag allocation."""
    full = "star_gravity_restart_formation_full"
    segment = "star_gravity_restart_formation_segment"
    joined = "star_gravity_restart_formation_joined"
    _run_sidecar_restart(
        "inputs/particles/star_gravity_restart_formation.athinput",
        full,
        segment,
        joined,
        "0.01",
        "0.02",
        "1",
        "2",
    )

    segment_header = _read_restart_header(f"rst_prtcl/{segment}.00001.rst_prtcl")
    assert segment_header["ncycle"] == 1
    assert segment_header["total_particles"] == 1
    assert segment_header["max_tag"] == 0
    assert math.isclose(segment_header["time"], 0.01, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(segment_header["formed"], 0.02, rel_tol=0.0,
                        abs_tol=1.0e-14)

    joined_header = _read_restart_header(f"rst_prtcl/{joined}.00002.rst_prtcl")
    assert joined_header["total_particles"] == 2
    assert joined_header["max_tag"] == 1
    assert math.isclose(joined_header["formed"], 0.04, rel_tol=0.0,
                        abs_tol=1.0e-14)

    full_row = _last_history_row(f"{full}.part.hst")
    joined_row = _last_history_row(f"{joined}.part.hst")
    assert int(round(full_row[2])) == 2
    assert math.isclose(full_row[4], 0.04, rel_tol=0.0, abs_tol=1.0e-14)
    _assert_rows_close(full_row, joined_row)


def test_star_gravity_sidecar_restart_rejects_time_mismatch():
    """The particle sidecar must match the mesh restart time and cycle."""
    segment = "star_gravity_restart_mismatch"
    bad = "star_gravity_restart_mismatch_bad"
    _clean_outputs(segment)
    _clean_outputs(bad)
    assert testutils.run(
        "inputs/particles/star_gravity_restart.athinput",
        [
            f"job/basename={segment}",
            "time/tlim=0.0025",
            "time/nlim=5",
            f"particles/particle_restart_file=rst_prtcl/{segment}.00001.rst_prtcl",
        ],
    )

    assert not testutils.run_command(
        [
            "./athena",
            "-r",
            f"rst/{segment}.00001.rst",
            f"job/basename={bad}",
            "particles/star_init=restart",
            f"particles/particle_restart_file=rst_prtcl/{segment}.00000.rst_prtcl",
        ]
    )


def test_star_gravity_sidecar_restart_rejects_trailing_data():
    """The loader should reject sidecars whose byte count does not match the header."""
    segment = "star_gravity_restart_trailing"
    bad = "star_gravity_restart_trailing_bad"
    corrupt_sidecar = Path("rst_prtcl") / f"{segment}.corrupt.rst_prtcl"
    _clean_outputs(segment)
    _clean_outputs(bad)
    corrupt_sidecar.unlink(missing_ok=True)
    assert testutils.run(
        "inputs/particles/star_gravity_restart.athinput",
        [
            f"job/basename={segment}",
            "time/tlim=0.0025",
            "time/nlim=5",
            f"particles/particle_restart_file=rst_prtcl/{segment}.00001.rst_prtcl",
        ],
    )
    source = Path("rst_prtcl") / f"{segment}.00001.rst_prtcl"
    corrupt_sidecar.write_bytes(source.read_bytes() + b"corrupt")
    assert not testutils.run_command(
        [
            "./athena",
            "-r",
            f"rst/{segment}.00001.rst",
            f"job/basename={bad}",
            "particles/star_init=restart",
            f"particles/particle_restart_file={corrupt_sidecar}",
        ]
    )


def test_star_gravity_sidecar_restart_rejects_duplicate_tags():
    """The loader should reject sidecars with non-unique particle identities."""
    segment = "star_gravity_restart_duplicate_tags"
    bad = "star_gravity_restart_duplicate_tags_bad"
    corrupt_sidecar = Path("rst_prtcl") / f"{segment}.corrupt.rst_prtcl"
    _clean_outputs(segment)
    _clean_outputs(bad)
    corrupt_sidecar.unlink(missing_ok=True)
    assert testutils.run(
        "inputs/particles/star_gravity_restart.athinput",
        [
            f"job/basename={segment}",
            "time/tlim=0.0025",
            "time/nlim=5",
            f"particles/particle_restart_file=rst_prtcl/{segment}.00001.rst_prtcl",
        ],
    )
    source = Path("rst_prtcl") / f"{segment}.00001.rst_prtcl"
    header = _read_restart_header(source)
    data = bytearray(source.read_bytes())
    second_tag = (header["header_size"] +
                  header["total_particles"]*header["nrdata"]*header["real_size"] +
                  (header["nidata"] + 1)*header["int_size"])
    struct.pack_into("@i", data, second_tag, 0)
    corrupt_sidecar.write_bytes(data)
    assert not testutils.run_command(
        [
            "./athena",
            "-r",
            f"rst/{segment}.00001.rst",
            f"job/basename={bad}",
            "particles/star_init=restart",
            f"particles/particle_restart_file={corrupt_sidecar}",
        ]
    )


def test_star_gravity_file_init_rejects_nonpositive_mass():
    """ASCII initialization should fail before evolution for invalid stellar mass."""
    particle_file = Path("star_gravity_invalid_mass.tbl")
    particle_file.write_text("0.5 0.5 0.5 0.0 0.0 0.0 -1.0 0.0\n")
    try:
        assert not testutils.run_command(
            [
                "./athena",
                "-i",
                "inputs/particles/star_gravity_external.athinput",
                "job/basename=star_gravity_invalid_mass",
                f"particles/particle_file={particle_file}",
            ]
        )
    finally:
        particle_file.unlink(missing_ok=True)
