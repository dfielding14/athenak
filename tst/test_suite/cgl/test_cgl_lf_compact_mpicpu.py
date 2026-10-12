"""Compact LF exchanges preserve transported and frozen fields, including ghosts."""

import numpy as np

from test_suite import testutils
from test_suite.cgl import cgl_lf_smr_conservation as conservation


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
