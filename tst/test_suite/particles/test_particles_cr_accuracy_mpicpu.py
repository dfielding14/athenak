"""MPI CPU accuracy tests for cosmic-ray tracer particles."""

import math
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
    """One, two, four, and eight-rank runs agree for tag-defined particles."""
    input_name = "cr_mpi_decomposition_invariance.athinput"
    run_dirs = []
    for ranks in (1, 2, 4, 8):
        run_dir = Path(f"cr_acc_mpi_invariance_n{ranks}")
        shutil.rmtree(run_dir, ignore_errors=True)
        command = [
            "mpirun", "-np", str(ranks), "./athena",
            "-i", str(INPUTS / input_name),
            "-d", str(run_dir),
            f"job/basename=cr_acc_mpi_invariance_n{ranks}",
            "time/nlim=256",
            "time/tlim=0.40",
            "output1/dt=0.20",
            "output2/dt=0.20",
            "particles/subcycle_cell_fraction=0.05",
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(result.stdout + result.stderr)
        run_dirs.append(run_dir)

    expected_species = [8, 8]
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
    for compared in moments[1:]:
        for species in range(2):
            assert moments[0][species]["count"] == compared[species]["count"]
            assert abs(moments[0][species]["mean_speed2"] -
                       compared[species]["mean_speed2"]) < 1.0e-12

    particles = [acu.particles_by_species_tag(summary) for summary in summaries]
    for compared in particles[1:]:
        assert set(particles[0]) == set(compared)
        for key in particles[0]:
            p0 = particles[0][key]
            p1 = compared[key]
            assert acu.vector_error(
                acu.particle_position(p0), acu.particle_position(p1)) < 1.0e-12
            assert acu.vector_error(
                acu.particle_velocity(p0), acu.particle_velocity(p1)) < 1.0e-12
            assert acu.vector_error(
                acu.particle_bfield(p0), acu.particle_bfield(p1)) < 1.0e-12
    eight_rank_snapshots = [
        acu.cr_tracer_inspect.read_prst_file(path)
        for path in (run_dirs[-1] / "prst").rglob("*.prst")
    ]
    assert any(len(snapshot["particles"]) == 0 for snapshot in eight_rank_snapshots)


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


def test_mirror_field_ensemble_exchange_is_valid_mpicpu():
    """Eight-rank mirror-field ensemble keeps deterministic particles valid."""
    _run_mpi(
        "cr_magnetic_mirror.athinput",
        [
            "job/basename=cr_acc_mirror_ensemble_mpi",
            "meshblock/nx1=8",
            "meshblock/nx2=8",
            "meshblock/nx3=8",
            "particles/ppc=0.25",
            "particles/check_consistency_mode=full",
            "problem/particle_position=tag_random",
            "problem/particle_seed=57721",
            "time/nlim=100",
            "time/tlim=0.12",
            "output1/dt=0.12",
            "output2/dt=0.12",
        ],
        ranks=8)
    summary = acu.latest_restart_summary()
    acu.cr_tracer_inspect.validate_expected_counts(summary, 1024, [1024])
    assert set(summary["rank_counts"]) == set(range(8))
    particle_map = acu.particles_by_species_tag(summary)
    assert len(particle_map) == 1024
    for particle in particle_map.values():
        assert all(math.isfinite(value) for value in particle["reals"])


def test_deterministic_pitch_angle_evolution_mpicpu():
    """Two-rank structured-field evolution retains valid deterministic particles."""
    _run_mpi(
        "cr_pitch_angle_decorrelation.athinput",
        [
            "job/basename=cr_acc_pitch_decorrelation_mpi",
            "time/nlim=2000",
            "time/tlim=0.64",
            "output1/dt=0.64",
            "output2/dt=0.64",
            "output3/dt=0.64",
        ],
        ranks=2)
    sequence = acu.restart_sequence()
    correlation = acu.pitch_angle_correlation(sequence[0], sequence[-1])
    summary = acu.latest_restart_summary()
    acu.cr_tracer_inspect.validate_expected_counts(summary, 1024, [1024])
    assert set(summary["rank_counts"]) == {0, 1}
    assert correlation < 0.995
