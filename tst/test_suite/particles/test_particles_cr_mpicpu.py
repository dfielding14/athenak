"""MPI CPU regression tests for cosmic-ray tracer particles."""

import re
import shutil
import subprocess
import sys
from pathlib import Path

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


def _run_mpi_capture(input_file, flags, ranks=2):
    command = ["mpirun", "-np", str(ranks), "./athena", "-i", input_file] + flags
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            "MPI run failed\nSTDOUT:\n"
            + result.stdout
            + "\nSTDERR:\n"
            + result.stderr)
    return result.stdout + result.stderr


def _assert_restart_counts(expected_species=EXPECTED_SPECIES):
    summary = cr_tracer_inspect.summarize_restart(Path("prst"))
    cr_tracer_inspect.validate_expected_counts(
        summary, sum(expected_species), expected_species)
    assert sum(summary["rank_counts"].values()) == sum(expected_species)
    return summary


def test_mpi_boris_amr_migration_mpicpu():
    """Run two-rank Boris AMR and validate migrated particle state."""
    _clean_particle_outputs()
    input_file = str(INPUTS / "particles" / "cr_tracer_boris_amr.athinput")
    flags = [
        "job/basename=cr_boris_mpi",
        "particles/check_consistency=true",
        "time/nlim=3",
        "time/tlim=0.03",
    ]
    assert testutils.mpi_run(input_file, flags, threads=2)
    summary = _assert_restart_counts()
    assert set(summary["rank_counts"]) == {0, 1}
    assert all(count > 0 for count in summary["rank_counts"].values())
    cr_tracer_inspect.validate_histograms(Path("."), EXPECTED_SPECIES, 8)
    cr_tracer_inspect.validate_moments(Path("."), EXPECTED_SPECIES)
    samples = cr_tracer_inspect.summarize_samples(Path("."))
    expected = cr_tracer_inspect.expected_sample_counts_from_restart(
        summary, -1, 4, 0)
    assert samples["total"] == expected["total"]
    assert samples["species_counts"] == expected["species_counts"]
    assert samples["rank_counts"] == expected["rank_counts"]
    assert {"mu", "bmag"}.issubset(set(samples["fields"]))


def test_mpi_refine_derefine_stress_mpicpu():
    """Force 4-rank AMR creation/deletion while conserving particles."""
    _clean_particle_outputs()
    output = _run_mpi_capture(
        str(INPUTS / "particles" / "cr_tracer_boris_amr_stress.athinput"),
        [
            "job/basename=cr_boris_stress_mpi",
            "particles/check_consistency_mode=full",
            "particles/validate_amr_lookup=true",
            "mesh_refinement/particle_load_weight=0.001",
        ],
        ranks=4)
    match = re.search(r"(\d+) MeshBlocks created, (\d+) deleted by AMR", output)
    assert match is not None
    created = int(match.group(1))
    deleted = int(match.group(2))
    assert created > 0
    assert deleted > 0
    summary = _assert_restart_counts()
    assert set(summary["rank_counts"]) == {0, 1, 2, 3}
    cr_tracer_inspect.validate_histograms(Path("."), EXPECTED_SPECIES, 8)
    cr_tracer_inspect.validate_joint_spectra(Path("."), EXPECTED_SPECIES, 8, 8)
    cr_tracer_inspect.validate_moments(Path("."), EXPECTED_SPECIES)


def test_mpi_drift_restart_round_trip_mpicpu():
    """Reload a two-rank drift particle restart and revalidate it."""
    _clean_particle_outputs()
    input_file = str(INPUTS / "particles" / "random_particle_drift.athinput")
    source_flags = [
        "job/basename=cr_drift_mpi_restart_src",
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
    _run_mpi_capture(input_file, source_flags)
    summary = _assert_restart_counts()
    rank0_file = next(
        file_path for file_path in summary["files"]
        if "rank_00000000" in str(file_path))

    restart_flags = [
        "job/basename=cr_drift_mpi_restart_reload",
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
        f"problem/prtcl_res_file={rank0_file}",
        "time/nlim=1",
        "time/tlim=0.03",
        "output1/file_type=prst",
        "output1/dt=0.01",
        "output2/file_type=ppd",
        "output2/dt=0.01",
    ]
    _run_mpi_capture(input_file, restart_flags)
    _assert_restart_counts()


def test_mpi_drift_subcycle_migration_mpicpu():
    """Subcycled high-velocity drift preserves ownership across MPI ranks."""
    _clean_particle_outputs()
    input_file = str(INPUTS / "particles" / "random_particle_drift.athinput")
    flags = [
        "job/basename=cr_drift_subcycle_mpi",
        "mesh/nx1=32",
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
    _run_mpi_capture(input_file, flags)
    summary = _assert_restart_counts([2])
    assert set(summary["rank_counts"]) == {0, 1}
