"""Compact LF exchanges preserve transported and frozen fields, including ghosts."""

from pathlib import Path
import struct

import numpy as np

from test_suite import testutils
from test_suite.cgl import cgl_lf_smr_conservation as conservation
from test_suite.cgl.test_cgl_c2p_pressure_floor_cpu import _run, _storage_executable


def test_cgl_lf_compact_exchange_matches_full_exchange(tmp_path):
    # Strict checks select compact exchange without changing the configured physics.
    common = ["mesh_refinement/refinement=none", "time/evolution=dynamic",
              "mhd/rsolver=hlle", "time/nlim=4", "time/tlim=1.0",
              "problem/vx0=0.08", "problem/vy0=0.05", "problem/vz0=0.03"]
    common += [f"mesh/nx{axis}=16" for axis in (1, 2, 3)]
    common += [f"meshblock/nx{axis}=8" for axis in (1, 2, 3)]
    reference = None
    # Eight 8^3 blocks, two ghosts, seven conserved variables and three CT faces.
    payload_bytes = 8*(7*12**3 + 3*13*12**2)*np.dtype("=f8").itemsize
    for ranks in (1, 4):
        histories = []
        for strict in (False, True):
            directory = tmp_path / f"ranks{ranks}-strict{strict}"
            directory.mkdir()
            testutils.mpi_run(conservation.INPUT, ["-d", str(directory), *common,
                              f"mhd/cgl_lf_strict_admissibility={str(strict).lower()}"],
                              threads=ranks)
            paths = sorted((directory / "rst").glob(f"{conservation.BASENAME}.*.rst"))
            assert len(paths) == 2
            initial, final = [conservation.read_state(path) for path in paths]
            assert initial["cycle"] == 0 and final["cycle"] == 4
            assert initial["nmb"] == final["nmb"] == 8
            for variable in (0, 1, 2, 3, 6):
                assert np.ptp(initial["u"][:, variable]) > 0
            assert np.max(abs(final["u"] - initial["u"])) > 1e-8
            history = testutils.athena_read.hst(
                str(directory / f"{conservation.BASENAME}.mhd.hst"))
            for name in ("lf_dfloor", "lf_pfloor", "lf_nonfin", "lf_nonpos", "lf_hardbd"):
                assert np.all(history[name] == 0), name
            rhs = history["lf_nstage"][-1]/16**3
            assert rhs > 2*final["cycle"], "the compact stage>1 path was not exercised"
            assert abs(history["lf_qprwrk"][-1]) + abs(history["lf_qpewrk"][-1]) > 0
            histories.append(history)
            # The existing reader validates the layout; compare the complete payload,
            # including outer ghosts and staggered fields omitted from its active view.
            payload = paths[-1].read_bytes()[-payload_bytes:]
            if reference is None:
                reference = (final, payload)
            conservation.assert_same(reference[0], final)
            assert payload == reference[1], "native conserved/CT/ghost payload changed"
        for name in ("time", "dt", "lf_nstage"):
            np.testing.assert_array_equal(histories[0][name], histories[1][name])


def test_cgl_lf_compact_exchange_preserves_density_repairs(tmp_path):
    """Run the storage custom pgen on one and two MPI ranks, optionally on HIP."""
    executable = _storage_executable(tmp_path, mpi=True)
    source = Path(conservation.INPUT).read_text()
    source += "\n<problem>\nlf_floor_case = 3\n<mhd>\niso_sound_speed = 1\n"
    input_file = tmp_path / "boundary-repair.athinput"
    input_file.write_text(source)
    common = ["mesh_refinement/refinement=none", "time/evolution=dynamic",
              "mhd/rsolver=hlle", "time/nlim=2", "time/tlim=1", "problem/user_hist=false",
              "mesh/nx1=16", "meshblock/nx1=8", "mhd/dfloor=1", "mhd/pfloor=1",
              "mhd/lf_k_parallel=0.01"]
    common += [f"{block}/nx{axis}=4" for block in ("mesh", "meshblock") for axis in (2, 3)]
    # Two 8x4x4 blocks, two ghosts, seven conserved variables and staggered CT faces.
    block_bytes = (7*12*8*8 + 13*8*8 + 12*9*8 + 12*8*9)*8
    for passive in (False, True):
        reference = None
        for ranks in (1, 2):
            for strict in (False, True):
                directory = tmp_path / f"passive{passive}-ranks{ranks}-strict{strict}"
                directory.mkdir()
                result = _run(["mpirun", "-np", str(ranks), str(executable),
                               "-i", str(input_file), *common,
                               f"mhd/passive={str(passive).lower()}",
                               f"mhd/cgl_lf_strict_admissibility={str(strict).lower()}"],
                              cwd=directory)
                (directory / "run.log").write_text(result.stdout + result.stderr)
                assert "CGL LF finite repair recorded: case=3" in result.stdout
                paths = sorted((directory / "rst").glob(f"{conservation.BASENAME}.*.rst"))
                assert len(paths) == 2
                history = testutils.athena_read.hst(
                    str(directory / f"{conservation.BASENAME}.mhd.hst"))
                assert history["time"][-1] > history["time"][0]
                assert history["lf_dfloor"][0] > 0
                assert history["lf_nstage"][-1] > history["lf_nstage"][0]
                for name in ("lf_pfloor", "lf_nonfin", "lf_nonpos"):
                    assert np.all(history[name] == 0), name
                payloads = []
                for path in paths:
                    raw = path.read_bytes()
                    offset = len(raw) - 2*block_bytes
                    assert struct.unpack_from("=Q", raw, offset-8)[0] == block_bytes
                    payloads.append(raw[offset:])
                assert payloads[0] != payloads[1], "scalar advection must evolve the repaired state"
                # Counts can differ: full exchange refreshes unused outer ghosts earlier.
                # The complete state, including both ghost layers and CT faces, must agree.
                if reference is None:
                    reference = payloads
                assert payloads == reference, "density repair changed with exchange mode or ranks"
