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
    for directory in ("ppd", "prst", "df", "dxh", "drh", "dparh", "pmom"):
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


@pytest.mark.parametrize("interpolation", ["tsc", "lin"])
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


def test_serial_uniform_boris_energy_cpu():
    """Verify uniform-field Boris pushing conserves particle speed."""
    _clean_particle_outputs()
    input_file = str(INPUTS / "particles" / "cr_tracer_boris_uniform.athinput")
    flags = [
        "job/basename=cr_boris_uniform_energy",
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
