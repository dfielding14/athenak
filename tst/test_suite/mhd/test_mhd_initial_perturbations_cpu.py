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
    vector_decomposition_metrics,
)


def check_initial_perturbation_run(input_file, basename, nlow, nhigh, cell_dims):
    """Run one input and verify amplitudes and Fourier support."""
    results = testutils.run(input_file)
    assert results, "Initial perturbation MHD smoke test failed."

    vtk_files = sorted(Path("vtk").glob(f"{basename}.mhd_w_bcc.*.vtk"))
    assert vtk_files, "Initial perturbation run did not write VTK diagnostics."

    data = read_vtk_scalars(vtk_files[0])
    assert data["cell_dims"] == cell_dims

    fields = data["scalars"]
    delta = density_contrast(fields)
    metrics = density_metrics(delta, nlow=nlow, nhigh=nhigh)

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

    return fields


def test_3d_run():
    """Run the 32^3 example and verify amplitudes and Fourier support."""
    check_initial_perturbation_run(
        "inputs/initial_perturbations.athinput",
        "InitialPerturbations",
        nlow=1,
        nhigh=4,
        cell_dims=(32, 32, 32),
    )


def test_2d_run():
    """Run the 128^2 example and verify amplitudes and Fourier support."""
    fields = check_initial_perturbation_run(
        "inputs/initial_perturbations_2d.athinput",
        "InitialPerturbations2D",
        nlow=1,
        nhigh=8,
        cell_dims=(128, 128, 1),
    )

    velocity = np.stack([fields["velx"], fields["vely"], fields["velz"]])
    velocity_metrics = vector_decomposition_metrics(velocity)
    assert abs(velocity_metrics["amplitude_solenoidal_fraction"] - 0.75) < 2.0e-2

    divb_files = sorted(Path("vtk").glob("InitialPerturbations2D.mhd_divb.*.vtk"))
    assert divb_files, "2D initial perturbation run did not write divB diagnostics."
    divb_fields = read_vtk_scalars(divb_files[0])["scalars"]
    assert np.max(np.abs(divb_fields["divb"])) < 1.0e-10
    assert np.sqrt(np.mean(divb_fields["divb"] ** 2)) < 1.0e-11
