"""Uniform CGL states exercise all RK4 conserved and magnetic registers."""

from pathlib import Path

import numpy as np
import pytest

import test_suite.testutils as testutils


@pytest.mark.parametrize(
    ("integrator", "fofc"),
    [("rk1", False), ("rk2", False), ("rk3", False), ("rk4", False), ("rk4", True)],
    ids=["rk1", "rk2", "rk3", "rk4", "rk4_fofc"],
)
@pytest.mark.parametrize("dimension", [1, 2, 3])
def test_uniform_cgl_registers(integrator, fofc, dimension, tmp_path):
    basename = f"cgl_registers_{dimension}d_{integrator}" + ("_fofc" if fofc else "")
    files = list(Path("tab").glob(f"{basename}.*.tab"))
    files += list(Path(".").glob(f"{basename}.*.hst"))
    files += list(Path(".").glob(f"{basename}.log"))
    for path in files:
        path.unlink()
    template = Path("../../../inputs/tests/cgl_amr_primitive_uniform.athinput")
    input_file = tmp_path / "athinput"
    text = template.read_text() + """
<mhd>
nscalars = 1
<output2>
file_type = tab
variable = mhd_u_bcc
data_format = %.17e
dcycle = 1
slice_x2 = 0.5
slice_x3 = 0.5
"""
    if fofc:
        text += """
<mesh>
nghost = 3
<mhd>
fofc = true
<output3>
file_type = log
dcycle = 1
"""
    input_file.write_text(text)
    testutils.run(
        str(input_file),
        [f"job/basename={basename}", f"time/integrator={integrator}",
         "time/nlim=4", "time/tlim=1", "mesh_refinement/refinement=none",
         "mesh/nx1=8", "meshblock/nx1=4",
         f"mesh/nx2={4 if dimension > 1 else 1}",
         f"meshblock/nx2={4 if dimension > 1 else 1}",
         f"mesh/nx3={4 if dimension > 2 else 1}",
         f"meshblock/nx3={4 if dimension > 2 else 1}",
         "problem/field_amp=0", "problem/vz0=0.03"],
    )
    snapshots = sorted(Path("tab").glob(f"{basename}.mhd_u_bcc.*.tab"))
    assert len(snapshots) >= 5
    initial = testutils.athena_read.tab(str(snapshots[0]))
    assert initial["cycle"] == 0
    assert initial["r_00"].min() > 0
    for snapshot in snapshots[1:]:
        final = testutils.athena_read.tab(str(snapshot))
        for name, expected in initial.items():
            if name not in ("time", "cycle"):
                np.testing.assert_allclose(final[name], expected,
                                           rtol=2.0e-12, atol=2.0e-12,
                                           err_msg=f"{snapshot}: {name}")
    assert final["cycle"] == 4
    history = testutils.athena_read.hst(f"{basename}.user.hst")
    assert np.max(np.abs(history["max_divb"])) < 2.0e-12
    assert np.max(history["bad_state"]) == 0
    if fofc:
        events = Path(f"{basename}.log")
        if events.exists():
            rows = [line.split() for line in events.read_text().splitlines()
                    if line.strip() and not line.startswith("#")]
            assert all(int(row[-1]) == 0 for row in rows), rows
