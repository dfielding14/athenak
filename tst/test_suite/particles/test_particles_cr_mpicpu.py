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
    for directory in ("ppd", "prst", "df", "dxh"):
        shutil.rmtree(directory, ignore_errors=True)
    for pattern in ("*.vtk", "*.pvtk", "*.xdmf"):
        for file_path in Path(".").glob(pattern):
            file_path.unlink()


def _run_mpi_capture(input_file, flags):
    command = ["mpirun", "-np", "2", "./athena", "-i", input_file] + flags
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


def test_mpi_refine_derefine_stress_mpicpu():
    """Force AMR creation and deletion while conserving particle counts."""
    _clean_particle_outputs()
    output = _run_mpi_capture(
        str(INPUTS / "particles" / "cr_tracer_boris_amr_stress.athinput"),
        ["job/basename=cr_boris_stress_mpi"])
    match = re.search(r"(\d+) MeshBlocks created, (\d+) deleted by AMR", output)
    assert match is not None
    created = int(match.group(1))
    deleted = int(match.group(2))
    assert created > 0
    assert deleted > 0
    _assert_restart_counts()
    cr_tracer_inspect.validate_histograms(Path("."), EXPECTED_SPECIES, 8)


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
