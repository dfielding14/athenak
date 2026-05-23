"""MPI CPU accuracy tests for cosmic-ray tracer particles."""

import re
import shutil
import subprocess
from pathlib import Path

from test_suite.particles import cr_accuracy_utils as acu


ROOT = Path(__file__).resolve().parents[3]
INPUTS = ROOT / "inputs" / "particles" / "accuracy"


def _run_mpi(input_name, flags, ranks=2):
    acu.clean_particle_outputs()
    command = ["mpirun", "-np", str(ranks), "./athena",
               "-i", str(INPUTS / input_name)] + flags
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            "MPI run failed\nSTDOUT:\n"
            + result.stdout
            + "\nSTDERR:\n"
            + result.stderr)
    return result.stdout + result.stderr


def test_uniform_gyro_amr_mpi_crossing_preserves_counts_mpicpu():
    """Uniform-B AMR/MPI crossing keeps particle counts and diagnostics valid."""
    output = _run_mpi(
        "cr_uniform_gyro_amr.athinput",
        [
            "job/basename=cr_acc_gyro_amr_mpi",
            "particles/check_consistency_mode=full",
            "particles/validate_amr_lookup=true",
        ],
        ranks=4)
    match = re.search(r"(\d+) MeshBlocks created, (\d+) deleted by AMR", output)
    assert match is not None
    assert int(match.group(1)) > 0
    assert int(match.group(2)) > 0

    expected_species = [512, 512]
    summary = acu.latest_restart_summary()
    acu.cr_tracer_inspect.validate_expected_counts(
        summary, sum(expected_species), expected_species)
    assert set(summary["rank_counts"]) == {0, 1, 2, 3}
    moments = acu.cr_tracer_inspect.validate_moments(Path("."), expected_species)
    for row in moments:
        assert row["mean_speed2"] > 0.0


def test_mpi_decomposition_reduced_moments_are_consistent_mpicpu():
    """One- and two-rank runs should agree on global reduced particle moments."""
    input_name = "cr_mpi_decomposition_invariance.athinput"
    run_dirs = []
    for ranks in (1, 2):
        run_dir = Path(f"cr_acc_mpi_invariance_n{ranks}")
        shutil.rmtree(run_dir, ignore_errors=True)
        command = [
            "mpirun", "-np", str(ranks), "./athena",
            "-i", str(INPUTS / input_name),
            "-d", str(run_dir),
            f"job/basename=cr_acc_mpi_invariance_n{ranks}",
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(result.stdout + result.stderr)
        run_dirs.append(run_dir)

    expected_species = [2, 2]
    summaries = [
        acu.cr_tracer_inspect.summarize_restart(path / "prst")
        for path in run_dirs
    ]
    for summary in summaries:
        acu.cr_tracer_inspect.validate_expected_counts(
            summary, sum(expected_species), expected_species)

    moments = [
        acu.cr_tracer_inspect.validate_moments(path, expected_species)
        for path in run_dirs
    ]
    for species in range(2):
        assert moments[0][species]["count"] == moments[1][species]["count"]
        assert abs(moments[0][species]["mean_speed2"] -
                   moments[1][species]["mean_speed2"]) < 1.0e-12

    particles = [acu.particles_by_species_tag(summary) for summary in summaries]
    assert set(particles[0]) == set(particles[1])
    for key in particles[0]:
        p0 = particles[0][key]
        p1 = particles[1][key]
        assert acu.vector_error(
            acu.particle_position(p0), acu.particle_position(p1)) < 1.0e-12
        assert acu.vector_error(
            acu.particle_velocity(p0), acu.particle_velocity(p1)) < 1.0e-12
        assert acu.vector_error(
            acu.particle_bfield(p0), acu.particle_bfield(p1)) < 1.0e-12


def test_amr_boundary_stress_conserves_particles_mpicpu():
    """Smooth-field AMR refinement/derefinement keeps diagnostics valid."""
    output = _run_mpi(
        "cr_amr_boundary_convergence.athinput",
        [
            "job/basename=cr_acc_amr_boundary_mpi",
            "particles/check_consistency_mode=full",
            "particles/validate_amr_lookup=true",
        ],
        ranks=2)
    match = re.search(r"(\d+) MeshBlocks created, (\d+) deleted by AMR", output)
    assert match is not None
    assert int(match.group(1)) > 0
    assert int(match.group(2)) > 0

    expected_species = [512, 512]
    summary = acu.latest_restart_summary()
    acu.cr_tracer_inspect.validate_expected_counts(
        summary, sum(expected_species), expected_species)
    assert set(summary["rank_counts"]) == {0, 1}
    moments = acu.cr_tracer_inspect.validate_moments(Path("."), expected_species)
    for row in moments:
        assert row["count"] == expected_species[row["species"]]
        assert row["mean_speed2"] > 0.0
