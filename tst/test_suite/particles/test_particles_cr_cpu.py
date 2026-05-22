"""CPU regression tests for cosmic-ray tracer particles."""

import math
import shutil
import sys
from pathlib import Path

import pytest
import test_suite.testutils as testutils


ROOT = Path(__file__).resolve().parents[3]
INPUTS = ROOT / "inputs"
sys.path.insert(0, str(ROOT / "scripts" / "particles"))
import cr_tracer_inspect  # noqa: E402


EXPECTED_TOTAL = 1024
EXPECTED_SPECIES = [512, 512]


def _clean_particle_outputs():
    for directory in ("ppd", "prst", "df", "dxh", "drh", "dparh", "pmom",
                      "pspec", "pspec2", "psamp"):
        shutil.rmtree(directory, ignore_errors=True)
    for pattern in ("*.vtk", "*.pvtk", "*.xdmf"):
        for file_path in Path(".").glob(pattern):
            file_path.unlink()


def _assert_restart_counts(expected_species=EXPECTED_SPECIES):
    summary = cr_tracer_inspect.summarize_restart(Path("prst"))
    cr_tracer_inspect.validate_expected_counts(
        summary, sum(expected_species), expected_species)
    return summary


def test_serial_drift_smoke_and_restart_cpu():
    """Run a small drift case and reload its particle restart file."""
    _clean_particle_outputs()
    input_file = str(INPUTS / "particles" / "random_particle_drift.athinput")
    flags = [
        "job/basename=cr_drift_restart_src",
        "mesh/nx1=16",
        "mesh/nx2=16",
        "mesh/nx3=16",
        "meshblock/nx1=8",
        "meshblock/nx2=8",
        "meshblock/nx3=8",
        "particles/ppc=0.125",
        "particles/nspecies=2",
        "particles/check_consistency=true",
        "time/nlim=2",
        "time/tlim=0.02",
        "output1/file_type=prst",
        "output1/dt=0.01",
        "output2/file_type=ppd",
        "output2/dt=0.01",
    ]
    assert testutils.run(input_file, flags)
    summary = _assert_restart_counts()
    assert len(summary["files"]) == 1

    restart_file = summary["files"][0]
    restart_flags = [
        "job/basename=cr_drift_restart_reload",
        "mesh/nx1=16",
        "mesh/nx2=16",
        "mesh/nx3=16",
        "meshblock/nx1=8",
        "meshblock/nx2=8",
        "meshblock/nx3=8",
        "particles/ppc=0.125",
        "particles/nspecies=2",
        "particles/check_consistency=true",
        "problem/prtcl_rst_flag=1",
        f"problem/prtcl_res_file={restart_file}",
        "time/nlim=1",
        "time/tlim=0.03",
        "output1/file_type=prst",
        "output1/dt=0.01",
        "output2/file_type=ppd",
        "output2/dt=0.01",
    ]
    assert testutils.run(input_file, restart_flags)
    _assert_restart_counts()


@pytest.mark.parametrize("interpolation", ["tsc", "lin", "trilinear"])
def test_serial_boris_amr_interpolation_cpu(interpolation):
    """Run Boris AMR with both supported interpolation paths."""
    _clean_particle_outputs()
    input_file = str(INPUTS / "particles" / "cr_tracer_boris_amr.athinput")
    flags = [
        f"job/basename=cr_boris_{interpolation}",
        f"particles/interpolation={interpolation}",
        "particles/check_consistency=true",
        "time/nlim=3",
        "time/tlim=0.03",
        "output1/dt=0.02",
        "output2/dt=0.02",
        "output3/dt=0.02",
        "output4/dt=0.02",
    ]
    assert testutils.run(input_file, flags)
    _assert_restart_counts()
    cr_tracer_inspect.validate_histograms(Path("."), EXPECTED_SPECIES, 8)
    cr_tracer_inspect.validate_moments(Path("."), EXPECTED_SPECIES)


def test_serial_multispecies_ppd_cpu():
    """Verify compact position output agrees with restart species counts."""
    _clean_particle_outputs()
    input_file = str(INPUTS / "particles" / "cr_tracer_boris_amr.athinput")
    flags = [
        "job/basename=cr_ppd_counts",
        "particles/check_consistency=true",
        "time/nlim=2",
        "time/tlim=0.02",
        "output1/dt=0.01",
        "output2/dt=0.01",
    ]
    assert testutils.run(input_file, flags)
    _assert_restart_counts()
    ppd_summary = cr_tracer_inspect.summarize_ppd(Path("."))
    assert ppd_summary["count"] == EXPECTED_TOTAL
    assert dict(ppd_summary["species_counts"]) == {0: 512, 1: 512}


def test_serial_selected_particle_sample_cpu():
    """Verify deterministic selected-field sample output and inspector readback."""
    _clean_particle_outputs()
    input_file = str(INPUTS / "particles" / "cr_tracer_boris_amr.athinput")
    flags = [
        "job/basename=cr_psamp_sample",
        "particles/check_consistency=true",
        "time/nlim=2",
        "time/tlim=0.02",
    ]
    assert testutils.run(input_file, flags)
    summary = _assert_restart_counts()
    sample_summary = cr_tracer_inspect.summarize_samples(Path("."))
    expected = cr_tracer_inspect.expected_sample_counts_from_restart(
        summary, sample_species=-1, sample_stride=4, sample_offset=0)
    assert sample_summary["total"] == expected["total"]
    assert sample_summary["species_counts"] == expected["species_counts"]
    assert {"x", "y", "z", "vx", "vy", "vz", "mu", "bmag"}.issubset(
        set(sample_summary["fields"]))
    for row in sample_summary["rows"]:
        assert row["rank"] == 0
        assert row["species"] in (0, 1)
        if "mu" in row["fields"]:
            assert abs(row["fields"]["mu"]) <= 1.0 + 1.0e-12
        if "bmag" in row["fields"]:
            assert row["fields"]["bmag"] >= 0.0


@pytest.mark.parametrize("interpolation", ["tsc", "lin", "trilinear"])
def test_serial_uniform_boris_energy_cpu(interpolation):
    """Verify uniform-field Boris pushing conserves speed for each gather."""
    _clean_particle_outputs()
    input_file = str(INPUTS / "particles" / "cr_tracer_boris_uniform.athinput")
    flags = [
        f"job/basename=cr_boris_uniform_energy_{interpolation}",
        f"particles/interpolation={interpolation}",
        "particles/check_consistency_mode=counts",
        "time/nlim=4",
        "time/tlim=0.2",
    ]
    assert testutils.run(input_file, flags)
    summary = _assert_restart_counts([1])
    particle = summary["restart"][0]["particles"][0]
    vx, vy, vz = particle["reals"][3:6]
    speed = math.sqrt(vx*vx + vy*vy + vz*vz)
    assert speed == pytest.approx(1.0, rel=1.0e-12, abs=1.0e-12)


def _field_gather_error(interpolation, nx):
    _clean_particle_outputs()
    input_file = str(INPUTS / "particles" / "cr_tracer_boris_uniform.athinput")
    ppc = 1.0/(nx*nx*nx)
    particle_x1 = 0.133
    particle_x2 = -0.071
    b0x = 0.3
    b0y = -0.2
    b0z = 1.0
    bgrad = 0.7
    flags = [
        f"job/basename=cr_gather_{interpolation}_{nx}",
        f"mesh/nx1={nx}",
        f"mesh/nx2={nx}",
        f"mesh/nx3={nx}",
        f"meshblock/nx1={nx}",
        f"meshblock/nx2={nx}",
        f"meshblock/nx3={nx}",
        f"particles/ppc={ppc:.17e}",
        f"particles/interpolation={interpolation}",
        "particles/check_consistency_mode=counts",
        "problem/particle_position=fixed",
        f"problem/particle_x1={particle_x1}",
        f"problem/particle_x2={particle_x2}",
        "problem/particle_x3=0.0",
        "problem/particle_velocity=uniform",
        "problem/v0x=0.0",
        "problem/v0y=0.0",
        "problem/v0z=0.0",
        "problem/B_profile=linear_cross",
        f"problem/B0x={b0x}",
        f"problem/B0y={b0y}",
        f"problem/B0z={b0z}",
        f"problem/Bgrad={bgrad}",
        "time/nlim=1",
        "time/tlim=0.01",
        "output1/file_type=prst",
        "output1/dt=0.001",
    ]
    assert testutils.run(input_file, flags)
    summary = _assert_restart_counts([1])
    bx, by, bz = summary["restart"][0]["particles"][0]["reals"][7:10]
    exact_bx = b0x + bgrad*particle_x2
    exact_by = b0y + bgrad*particle_x1
    exact_bz = b0z
    return math.sqrt((bx - exact_bx)**2 + (by - exact_by)**2 +
                     (bz - exact_bz)**2)


def test_serial_trilinear_smooth_field_convergence_cpu():
    """Verify trilinear gather improves a smooth cross-field over legacy lin."""
    legacy_errors = [_field_gather_error("lin", nx) for nx in (8, 16, 32)]
    trilinear_errors = [
        _field_gather_error("trilinear", nx) for nx in (8, 16, 32)]
    assert legacy_errors[1] < legacy_errors[0]
    assert legacy_errors[2] < legacy_errors[1]
    assert max(trilinear_errors) < 1.0e-12
    assert max(trilinear_errors) < 1.0e-8*min(legacy_errors)


@pytest.mark.parametrize("quantity", ["p", "E", "logE", "mu", "magnetic_moment"])
def test_serial_particle_spectrum_quantities_cpu(quantity):
    """Verify reduced spectrum outputs normalize for supported quantities."""
    _clean_particle_outputs()
    input_file = str(INPUTS / "particles" / "cr_tracer_boris_uniform.athinput")
    flags = [
        f"job/basename=cr_pspec_{quantity}",
        "particles/check_consistency_mode=counts",
        "time/nlim=1",
        "time/tlim=0.05",
        "output3/file_type=pspec",
        "output3/dt=0.01",
        f"output3/quantity={quantity}",
        "output3/nbin=8",
    ]
    if quantity == "logE":
        flags += ["output3/vmin=-2.0", "output3/vmax=1.0"]
    assert testutils.run(input_file, flags)
    _assert_restart_counts([1])
    cr_tracer_inspect.validate_spectra(Path("."), [1], 8)


@pytest.mark.parametrize("quantity", ["mu_p", "mu_E", "vpar_vperp"])
def test_serial_particle_joint_spectrum_quantities_cpu(quantity):
    """Verify reduced 2D spectrum outputs normalize for supported quantities."""
    _clean_particle_outputs()
    input_file = str(INPUTS / "particles" / "cr_tracer_boris_uniform.athinput")
    flags = [
        f"job/basename=cr_pspec2_{quantity}",
        "particles/check_consistency_mode=counts",
        "time/nlim=1",
        "time/tlim=0.05",
        "output4/file_type=pspec2",
        "output4/dt=0.01",
        f"output4/quantity={quantity}",
        "output4/nbin1=8",
        "output4/nbin2=8",
    ]
    if quantity == "vpar_vperp":
        flags += ["output4/vmin1=-2.0", "output4/vmax1=2.0"]
    assert testutils.run(input_file, flags)
    _assert_restart_counts([1])
    cr_tracer_inspect.validate_joint_spectra(Path("."), [1], 8, 8)


def test_serial_uniform_boris_subcycle_phase_cpu():
    """Verify subcycling improves uniform-field Boris phase accuracy."""
    input_file = str(INPUTS / "particles" / "cr_tracer_boris_uniform.athinput")

    def run_case(basename, extra_flags):
        _clean_particle_outputs()
        flags = [
            f"job/basename={basename}",
            "particles/interpolation=tsc",
            "particles/check_consistency_mode=counts",
            "time/nlim=8",
            "time/tlim=0.4",
        ] + extra_flags
        assert testutils.run(input_file, flags)
        summary = _assert_restart_counts([1])
        restart = summary["restart"][0]
        particle = restart["particles"][0]
        vx, vy = particle["reals"][3:5]
        angle = math.atan2(vy, vx)
        target = -restart["time"]
        return abs(math.atan2(math.sin(angle - target), math.cos(angle - target)))

    base_error = run_case("cr_boris_phase_base", [])
    subcycle_error = run_case(
        "cr_boris_phase_subcycle",
        [
            "particles/subcycle=true",
            "particles/subcycle_max_steps=64",
            "particles/subcycle_gyro_fraction=0.005",
        ])
    assert subcycle_error < 0.25*base_error


def test_serial_periodic_displacement_and_motion_guard_cpu():
    """Check periodic drift displacement and the optional large-motion guard."""
    _clean_particle_outputs()
    input_file = str(INPUTS / "particles" / "random_particle_drift.athinput")
    flags = [
        "job/basename=cr_periodic_drift",
        "mesh/nx1=8",
        "mesh/nx2=8",
        "mesh/nx3=8",
        "meshblock/nx1=8",
        "meshblock/nx2=8",
        "meshblock/nx3=8",
        "particles/ppc=0.001953125",
        "particles/nspecies=1",
        "particles/check_consistency_mode=local",
        "problem/particle_position=center",
        "problem/particle_velocity=uniform",
        "problem/v0x=30.0",
        "problem/v0y=0.0",
        "problem/v0z=0.0",
        "time/nlim=1",
        "time/tlim=0.125",
        "output1/file_type=prst",
        "output1/dt=0.01",
    ]
    assert testutils.run(input_file, flags)
    summary = _assert_restart_counts([1])
    particle = summary["restart"][0]["particles"][0]
    x = particle["reals"][0]
    dx = particle["reals"][10]
    assert -0.5 <= x <= 0.5
    assert dx == pytest.approx(0.75, rel=1.0e-12, abs=1.0e-12)

    _clean_particle_outputs()
    guard_flags = [
        "job/basename=cr_motion_guard",
        "mesh/nx1=16",
        "mesh/nx2=8",
        "mesh/nx3=8",
        "meshblock/nx1=8",
        "meshblock/nx2=8",
        "meshblock/nx3=8",
        "particles/ppc=0.0009765625",
        "particles/nspecies=1",
        "particles/check_motion_bounds=true",
        "problem/particle_position=center",
        "problem/particle_velocity=uniform",
        "problem/v0x=100.0",
        "problem/v0y=0.0",
        "problem/v0z=0.0",
        "time/nlim=1",
        "time/tlim=0.125",
        "output1/file_type=prst",
        "output1/dt=0.01",
    ]
    with pytest.raises(RuntimeError):
        testutils.run(input_file, guard_flags)

    _clean_particle_outputs()
    subcycle_flags = [
        "job/basename=cr_motion_subcycle",
        "mesh/nx1=16",
        "mesh/nx2=8",
        "mesh/nx3=8",
        "meshblock/nx1=8",
        "meshblock/nx2=8",
        "meshblock/nx3=8",
        "particles/ppc=0.0009765625",
        "particles/nspecies=1",
        "particles/check_consistency_mode=local",
        "particles/check_motion_bounds=true",
        "particles/subcycle=true",
        "particles/subcycle_max_steps=256",
        "problem/particle_position=center",
        "problem/particle_velocity=uniform",
        "problem/v0x=100.0",
        "problem/v0y=0.0",
        "problem/v0z=0.0",
        "time/nlim=1",
        "time/tlim=0.125",
        "output1/file_type=prst",
        "output1/dt=0.01",
    ]
    assert testutils.run(input_file, subcycle_flags)
    summary = _assert_restart_counts([1])
    particle = summary["restart"][0]["particles"][0]
    x = particle["reals"][0]
    dx = particle["reals"][10]
    assert -0.5 <= x <= 0.5
    assert dx == pytest.approx(2.5, rel=1.0e-12, abs=1.0e-12)

    strict_flags = list(subcycle_flags)
    strict_flags[0] = "job/basename=cr_motion_subcycle_strict"
    strict_flags[strict_flags.index("particles/subcycle_max_steps=256")] = (
        "particles/subcycle_max_steps=8")
    with pytest.raises(RuntimeError):
        testutils.run(input_file, strict_flags)
