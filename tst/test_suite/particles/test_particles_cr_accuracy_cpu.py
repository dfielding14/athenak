"""CPU accuracy tests for cosmic-ray tracer particles."""

from pathlib import Path

import pytest
import test_suite.testutils as testutils

from test_suite.particles import cr_accuracy_utils as acu


ROOT = Path(__file__).resolve().parents[3]
INPUTS = ROOT / "inputs" / "particles" / "accuracy"


def _run(input_name, flags):
    acu.clean_particle_outputs()
    assert testutils.run(str(INPUTS / input_name), flags)
    return acu.latest_restart_summary()


def _magnetic_moment_proxy(velocity, bfield):
    bmag = acu.norm3(bfield)
    speed = acu.norm3(velocity)
    vpar = sum(velocity[i]*bfield[i] for i in range(3))/bmag
    return (speed*speed - vpar*vpar)/bmag


def _one_particle_state(summary):
    particle = acu.single_particle(summary)
    return {
        "position": acu.particle_position(particle),
        "velocity": acu.particle_velocity(particle),
        "bfield": acu.particle_bfield(particle),
    }


def test_uniform_gyro_phase_converges_with_particle_timestep_cpu():
    """Verify uniform-field Boris phase error decreases with particle dt."""
    cfl_values = [0.04, 0.02, 0.01]
    errors = []
    for cfl in cfl_values:
        summary = _run(
            "cr_uniform_gyro.athinput",
            [
                f"job/basename=cr_acc_gyro_cfl_{cfl}",
                f"particles/cfl_part={cfl}",
                "particles/interpolation=tsc",
                "time/nlim=200",
                "time/tlim=0.32",
                "output1/dt=0.32",
                "output2/dt=0.32",
            ])
        particle = acu.single_particle(summary)
        vx, vy, _ = acu.particle_velocity(particle)
        errors.append(acu.phase_error(vx, vy, summary["restart"][0]["time"]))
        assert acu.speed(acu.particle_velocity(particle)) == pytest.approx(
            (1.0**2 + 0.2**2)**0.5, rel=1.0e-12, abs=1.0e-12)

    assert errors[1] < errors[0]
    assert errors[2] < errors[1]
    slope = acu.fit_log_slope(cfl_values, errors)
    assert slope > 1.7


@pytest.mark.parametrize("interpolation", ["lin", "trilinear", "tsc"])
def test_uniform_gyro_methods_agree_in_uniform_field_cpu(interpolation):
    """Uniform B should not discriminate field-gather methods."""
    summary = _run(
        "cr_uniform_gyro.athinput",
        [
            f"job/basename=cr_acc_gyro_{interpolation}",
            f"particles/interpolation={interpolation}",
            "time/nlim=20",
            "time/tlim=0.2",
            "output1/dt=0.2",
            "output2/dt=0.2",
        ])
    particle = acu.single_particle(summary)
    assert acu.vector_error(acu.particle_bfield(particle), (0.0, 0.0, 1.0)) < 1.0e-14
    assert acu.speed(acu.particle_velocity(particle)) == pytest.approx(
        (1.0**2 + 0.2**2)**0.5, rel=1.0e-12, abs=1.0e-12)


def test_linear_field_trilinear_gather_is_exact_cpu():
    """Trilinear gather should recover the linear-cross field at a fixed particle."""
    summary = _run(
        "cr_linear_gather.athinput",
        [
            "job/basename=cr_acc_linear_trilinear",
            "particles/interpolation=trilinear",
            "time/nlim=1",
            "time/tlim=0.01",
        ])
    error = acu.field_sample_error(
        summary, "linear_cross", (0.3, -0.2, 1.0), bgrad=0.7)
    assert error < 1.0e-12


def test_manufactured_divergence_free_gather_converges_cpu():
    """Smooth manufactured-field sampling should converge with resolution."""
    nxs = [8, 16, 32]
    errors = []
    for nx in nxs:
        ppc = 1.0/(nx*nx*nx)
        summary = _run(
            "cr_manufactured_divb_free.athinput",
            [
                f"job/basename=cr_acc_manufactured_{nx}",
                f"mesh/nx1={nx}",
                f"mesh/nx2={nx}",
                f"mesh/nx3={nx}",
                f"meshblock/nx1={nx}",
                f"meshblock/nx2={nx}",
                f"meshblock/nx3={nx}",
                f"particles/ppc={ppc:.17e}",
                "particles/interpolation=trilinear",
                "time/nlim=1",
                "time/tlim=0.01",
            ])
        errors.append(acu.field_sample_error(
            summary, "sinusoidal_divb_free", (0.1, -0.2, 1.0),
            bamp=0.025, bwave=1.0))

    assert errors[1] < errors[0]
    assert errors[2] < errors[1]
    slope = -acu.fit_log_slope(nxs, errors)
    assert slope > 1.5


def test_isotropic_ensemble_moments_and_joint_spectrum_cpu():
    """Isotropic initialization should produce near-isotropic reduced moments."""
    summary = _run(
        "cr_isotropic_ensemble.athinput",
        [
            "job/basename=cr_acc_isotropic",
            "time/nlim=1",
            "time/tlim=0.02",
        ])
    expected_species = [1024, 1024]
    acu.cr_tracer_inspect.validate_expected_counts(
        summary, sum(expected_species), expected_species)
    acu.cr_tracer_inspect.validate_histograms(Path("."), expected_species, 16)
    acu.cr_tracer_inspect.validate_spectra(Path("."), expected_species, 16)
    acu.cr_tracer_inspect.validate_joint_spectra(Path("."), expected_species, 16, 8)
    moments = acu.cr_tracer_inspect.validate_moments(Path("."), expected_species)
    for row in moments:
        assert abs(row["mean_mu"]) < 0.06
        assert row["mean_mu2"] == pytest.approx(1.0/3.0, abs=0.08)


def test_smooth_field_orbit_approaches_high_resolution_reference_cpu():
    """Coupled pusher/gather error should shrink toward a finer reference run."""
    states = []
    for nx in [16, 32, 64]:
        ppc = 1.0/(nx*nx*nx)
        dt_fraction = 0.2/nx
        summary = _run(
            "cr_smooth_orbit_reference.athinput",
            [
                f"job/basename=cr_acc_smooth_orbit_{nx}",
                f"mesh/nx1={nx}",
                f"mesh/nx2={nx}",
                f"mesh/nx3={nx}",
                f"meshblock/nx1={nx}",
                f"meshblock/nx2={nx}",
                f"meshblock/nx3={nx}",
                f"particles/ppc={ppc:.17e}",
                f"particles/cfl_part={dt_fraction}",
                f"particles/subcycle_gyro_fraction={dt_fraction}",
                "time/nlim=200",
                "time/tlim=0.1",
            ])
        states.append(_one_particle_state(summary))

    err16 = acu.vector_error(states[0]["velocity"], states[2]["velocity"])
    err32 = acu.vector_error(states[1]["velocity"], states[2]["velocity"])
    assert err32 < err16
    for state in states:
        assert acu.speed(state["velocity"]) == pytest.approx(
            (0.6**2 + 0.3**2 + 0.4**2)**0.5, rel=1.0e-12, abs=1.0e-12)


def test_magnetic_mirror_invariants_are_bounded_cpu():
    """Mirror-field orbit should conserve speed and keep moment drift bounded."""
    summary = _run(
        "cr_magnetic_mirror.athinput",
        [
            "job/basename=cr_acc_mirror_invariants",
            "time/nlim=100",
            "time/tlim=0.12",
        ])
    state = _one_particle_state(summary)
    initial_velocity = (0.7, 0.0, 0.25)
    initial_bfield = acu.exact_bfield(
        "mirror", 0.05, -0.03, -0.12, (0.0, 0.0, 1.0), bgrad=0.8)
    initial_moment = _magnetic_moment_proxy(initial_velocity, initial_bfield)
    final_moment = _magnetic_moment_proxy(state["velocity"], state["bfield"])
    assert acu.speed(state["velocity"]) == pytest.approx(
        acu.speed(initial_velocity), rel=1.0e-12, abs=1.0e-12)
    assert abs(final_moment - initial_moment)/initial_moment < 0.02


def test_gradb_final_state_responds_to_gradient_strength_cpu():
    """Grad-B setup should show gradient-dependent drift without speed growth."""
    states = []
    for bgrad in [0.0, 0.5, 1.0]:
        summary = _run(
            "cr_gradb_curvature_drift.athinput",
            [
                f"job/basename=cr_acc_gradb_{bgrad}",
                f"problem/Bgrad={bgrad}",
                "particles/cfl_part=0.01",
                "particles/subcycle_gyro_fraction=0.01",
                "time/nlim=1000",
                "time/tlim=0.64",
            ])
        states.append(_one_particle_state(summary))

    speed0 = (0.8**2 + 0.2**2)**0.5
    for state in states:
        assert acu.speed(state["velocity"]) == pytest.approx(
            speed0, rel=1.0e-12, abs=1.0e-12)
    weak_delta = acu.vector_error(states[1]["velocity"], states[0]["velocity"])
    strong_delta = acu.vector_error(states[2]["velocity"], states[0]["velocity"])
    assert strong_delta > weak_delta > 0.0


def test_frozen_turbulent_field_diagnostics_cpu():
    """Production-like frozen field should keep reduced diagnostics normalized."""
    summary = _run(
        "cr_frozen_turbulent_field.athinput",
        [
            "job/basename=cr_acc_turbulent_diag",
            "time/nlim=4",
            "time/tlim=0.08",
        ])
    expected_species = [128, 128]
    acu.cr_tracer_inspect.validate_expected_counts(
        summary, sum(expected_species), expected_species)
    moments = acu.cr_tracer_inspect.validate_moments(Path("."), expected_species)
    acu.cr_tracer_inspect.validate_joint_spectra(Path("."), expected_species, 16, 16)
    for row in moments:
        assert row["mean_speed2"] == pytest.approx(1.0, rel=1.0e-12, abs=1.0e-12)
