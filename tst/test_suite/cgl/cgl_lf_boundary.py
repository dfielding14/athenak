"""Independent fixed-state and equivalent user/inflow LF boundary references."""
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "vis/python"))
import bin_convert  # noqa: E402

INPUT_NAME = "cgl_lf_boundary.athinput"
INPUT = str(ROOT / "inputs/tests" / INPUT_NAME)
BASENAME = "cgl_lf_boundary"
EXPECTED = {
    "dens": 1.0, "velx": 0.08, "vely": 0.05, "velz": 0.0,
    "eint": 1.0, "p_perp": 1.01, "s_00": 0.5,
    "bcc1": 0.7, "bcc2": 0.2, "bcc3": -0.15,
}


def flags(directory, integrator, geometry, boundary, amplitude=0.0):
    options = ["-d", str(directory), f"job/basename={BASENAME}",
               f"mhd/cgl_heat_flux_integrator={integrator}",
               f"mesh/ix1_bc={boundary}", f"problem/amp={amplitude}"]
    if integrator == "explicit":
        options.append("time/sts_integrator=none")
    if geometry != "smr":
        options.append("mesh_refinement/refinement=none")
    if geometry == "1d":
        options.extend(["mesh/nx2=1", "meshblock/nx2=1"])
    elif geometry == "mixed2d":
        options.extend(["mesh/ix2_bc=outflow", "mesh/ox2_bc=outflow"])
    elif geometry == "3d":
        for direction in (1, 2, 3):
            options.extend([f"mesh/nx{direction}=8", f"meshblock/nx{direction}=4",
                            f"mesh/ix{direction}_bc={boundary}",
                            f"mesh/ox{direction}_bc={boundary}"])
    else:
        assert geometry == "smr"
    if amplitude:
        # B normal to the 1D boundary stays constant, so the stored conserved
        # inflow and prescribed primitive callback are independent equal states.
        assert geometry == "1d"
        options.extend(["problem/guide_b2=0", "problem/guide_b3=0"])
    return options


def states(directory):
    paths = sorted((directory / "bin").glob(f"{BASENAME}.state.*.bin"))
    assert paths, "boundary fixture produced no primitive fields"
    data = [bin_convert.read_binary(str(path)) for path in paths]
    assert data[0]["cycle"] == 0
    assert data[-1]["cycle"] == 4
    for state in data:
        for values in state["mb_data"].values():
            assert np.all(np.isfinite(values))
    return data


def check_uniform(directory):
    data = states(directory)
    for state in data:
        assert set(state["var_names"]) == set(EXPECTED)
        halo = tuple(slice(1, size + 3) if size > 1 else slice(None)
                     for size in (state["nx3_mb"], state["nx2_mb"], state["nx1_mb"]))
        for block in range(state["n_mbs"]):
            for name, value in EXPECTED.items():
                np.testing.assert_allclose(
                    np.asarray(state["mb_data"][name][block])[halo],
                    np.float32(value), rtol=0.0, atol=2e-7,
                    err_msg=f"{name}: active/ghost state at cycle {state['cycle']}",
                )
    return data[-1]


def same_state(left, right):
    assert left["cycle"] == right["cycle"]
    assert left["time"] == right["time"]
    assert left["var_names"] == right["var_names"]
    assert left["n_mbs"] == right["n_mbs"]
    orders = [sorted(range(s["n_mbs"]), key=lambda m: tuple(s["mb_logical"][m]))
              for s in (left, right)]
    for lm, rm in zip(*orders):
        np.testing.assert_array_equal(left["mb_logical"][lm], right["mb_logical"][rm])
        for name in left["var_names"]:
            np.testing.assert_allclose(left["mb_data"][name][lm],
                                       right["mb_data"][name][rm], rtol=0.0, atol=2e-7)
