"""Regression test for one-time initial Fourier perturbations in MHD."""

from pathlib import Path
import sys

import numpy as np

import test_suite.testutils as testutils

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
from plot_initial_perturbations_example import (  # noqa: E402
    density_contrast,
    density_metrics,
    read_vtk_scalars,
)


def test_run():
    """Run the example and verify amplitudes and Fourier support."""
    results = testutils.run("inputs/initial_perturbations.athinput")
    assert results, "Initial perturbation MHD smoke test failed."

    vtk_files = sorted(Path("vtk").glob("InitialPerturbations.mhd_w_bcc.*.vtk"))
    assert vtk_files, "Initial perturbation run did not write VTK diagnostics."

    data = read_vtk_scalars(vtk_files[0])
    fields = data["scalars"]
    delta = density_contrast(fields)
    metrics = density_metrics(delta, nlow=1, nhigh=4)

    assert abs(metrics["mean_delta"]) < 1.0e-7
    assert abs(metrics["rms_delta"] - 1.0e-2) < 5.0e-6
    assert metrics["zero_mode_fraction"] < 1.0e-10
    assert metrics["leakage_fraction"] < 1.0e-10

    velocity_rms = np.sqrt(
        np.mean(fields["velx"] ** 2 + fields["vely"] ** 2 + fields["velz"] ** 2)
    )
    assert abs(velocity_rms - 1.0e-3) < 5.0e-6

    magnetic_rms = np.sqrt(
        np.mean(
            (fields["bcc1"] - np.mean(fields["bcc1"])) ** 2
            + (fields["bcc2"] - np.mean(fields["bcc2"])) ** 2
            + (fields["bcc3"] - np.mean(fields["bcc3"])) ** 2
        )
    )
    assert abs(magnetic_rms - 1.0e-3) < 5.0e-6
