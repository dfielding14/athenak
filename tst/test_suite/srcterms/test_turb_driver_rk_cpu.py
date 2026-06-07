"""Regression tests for turbulence-driver task ordering and RK weighting."""

from pathlib import Path
import re
import shutil
import subprocess
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
INPUT = ROOT / "tst/inputs/turb_driver_rk.athinput"
RUN_ROOT = Path("run_turb_driver_rk_cpu")

sys.path.insert(0, str(ROOT / "scripts"))
from plot_initial_perturbations_example import read_vtk_scalars  # noqa: E402


def _read_vector(run_dir, variable, names):
    files = sorted((run_dir / "vtk").glob(f"TurbDriverRK.{variable}.*.vtk"))
    assert files, f"No {variable} VTK output found in {run_dir}"
    fields = read_vtk_scalars(files[-1])["scalars"]
    return np.stack([fields[name] for name in names])


def _run(integrator, tcorr, label):
    run_dir = RUN_ROOT / label
    shutil.rmtree(run_dir, ignore_errors=True)
    run_dir.mkdir(parents=True)
    command = [
        "./athena",
        "-i",
        str(INPUT),
        "-d",
        str(run_dir),
        f"time/integrator={integrator}",
        f"turb_driving/tcorr={tcorr}",
    ]
    result = subprocess.run(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False
    )
    assert result.returncode == 0, result.stdout + result.stderr
    force = _read_vector(run_dir, "turb_force", ("force1", "force2", "force3"))
    momentum = _read_vector(run_dir, "hydro_u", ("mom1", "mom2", "mom3"))
    return force, momentum


def test_ou_state_advances_once_per_full_step():
    """The one-step O-U state is independent of the number of RK stages."""
    shutil.rmtree(RUN_ROOT, ignore_errors=True)
    try:
        forces = {}
        for integrator in ("rk1", "rk2", "rk3"):
            forces[integrator], _ = _run(
                integrator, tcorr=0.5, label=f"ou_{integrator}"
            )

        np.testing.assert_allclose(forces["rk2"], forces["rk1"], rtol=0.0, atol=0.0)
        np.testing.assert_allclose(forces["rk3"], forces["rk1"], rtol=0.0, atol=0.0)
    finally:
        shutil.rmtree(RUN_ROOT, ignore_errors=True)


def test_forcing_uses_explicit_rk_source_weights():
    """At tiny CFL, one-step forcing agrees across RK1, RK2, and RK3."""
    shutil.rmtree(RUN_ROOT, ignore_errors=True)
    try:
        norms = {}
        for integrator in ("rk1", "rk2", "rk3"):
            _, momentum = _run(
                integrator, tcorr=0.0, label=f"source_{integrator}"
            )
            norms[integrator] = np.linalg.norm(momentum)

        reference = norms["rk1"]
        assert reference > 0.0
        assert abs(norms["rk2"] / reference - 1.0) < 2.0e-4
        assert abs(norms["rk3"] / reference - 1.0) < 2.0e-4
    finally:
        shutil.rmtree(RUN_ROOT, ignore_errors=True)


def test_forcing_tasks_follow_each_supported_rk_update():
    """Hydro, MHD, and ion-neutral insert forcing before their source task."""
    source = (ROOT / "src/srcterms/turb_driver.cpp").read_text()
    pairs = (
        ("phydro", "rkupdt", "srctrms"),
        ("pmhd", "rkupdt", "srctrms"),
        ("pionn", "n_rkupdt", "n_srctrms"),
    )
    for owner, update, source_task in pairs:
        pattern = (
            rf"{owner}->id\.{update}\s*,\s*"
            rf"pmy_pack->{owner}->id\.{source_task}"
        )
        assert re.search(pattern, source)
