"""CPU regression tests for second-moment Ito mass-flux tracers."""

from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np

import test_suite.testutils as testutils
from test_suite.particles.ito_restart_test_utils import make_legacy_restart


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
from read_prtcl_thermo_history import read_history  # noqa: E402


INPUT = str(ROOT / "inputs/particles/ito_tracers.athinput")
MHD_INPUT = str(ROOT / "inputs/particles/lagrangian_mc_thermo_mhd.athinput")
RUN_DIR = Path("run_particles_ito_cpu")
RUN_RK3 = Path("run_particles_ito_rk3_cpu")
RUN_MHD = Path("run_particles_ito_mhd_cpu")


def _cycle_displacement(data, component):
    initial = {
        int(tag): value
        for tag, value in zip(data["tag"][data["cycle"] == 0],
                              data[component][data["cycle"] == 0])
    }
    final = {
        int(tag): value
        for tag, value in zip(data["tag"][data["cycle"] == 1],
                              data[component][data["cycle"] == 1])
    }
    displacement = np.array([final[tag] - initial[tag] for tag in initial])
    return (displacement + 0.5) % 1.0 - 0.5


def _assert_uniform_flow_moments(
    run_dir, velocities=(0.5, 0.0, 0.0), cell_widths=(1.0 / 32.0, 1.0 / 8.0, 1.0)
):
    history = read_history(
        run_dir / "prtcl_thermo_history/ito_tracers.prtcl_thermo_history.thp"
    )
    dt = np.max(history["time"])
    for component, velocity, cell_width in zip(
        ("x1", "x2", "x3"), velocities, cell_widths
    ):
        displacement = _cycle_displacement(history, component)
        assert displacement.size == 4096
        if velocity == 0.0:
            assert np.max(np.abs(displacement)) == 0.0
            continue
        courant = velocity * dt / cell_width
        expected_mean = velocity * dt
        expected_variance = cell_width**2 * courant * (1.0 - courant)
        mean_tolerance = 5.0 * np.sqrt(expected_variance / displacement.size)
        assert abs(np.mean(displacement) - expected_mean) < mean_tolerance
        assert abs(np.var(displacement) - expected_variance) < 0.08 * expected_variance


def _final_state(run_dir):
    history = read_history(
        run_dir / "prtcl_thermo_history/ito_tracers.prtcl_thermo_history.thp"
    )
    final = history["cycle"] == np.max(history["cycle"])
    order = np.argsort(history["tag"][final])
    return history["tag"][final][order], np.column_stack(
        [history[name][final][order] for name in ("x1", "x2", "x3")]
    )


def test_ito2_uniform_flow_moments_restart_and_validation():
    """Ito-2 matches MC moments, restarts, and explicitly rejects Ito-3."""
    shutil.rmtree(RUN_DIR, ignore_errors=True)
    try:
        assert testutils.run(INPUT, ["-d", str(RUN_DIR)])
        _assert_uniform_flow_moments(RUN_DIR)

        restart = RUN_DIR / "rst/rank_00000000/ito_tracers.00001.rst"
        assert testutils.run_command(
            ["./athena", "-r", str(restart), "-d", str(RUN_DIR / "restart"),
             "time/nlim=2"]
        )
        mismatch = subprocess.run(
            [
                "./athena",
                "-r",
                str(restart),
                "-d",
                str(RUN_DIR / "restart_wrong_covariance"),
                "particles/ito_covariance_model=full_finite_step",
                "time/nlim=2",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        assert mismatch.returncode != 0
        assert "particle restart covariance model does not match" in (
            mismatch.stdout + mismatch.stderr
        )
        assert not testutils.run_command(
            ["./athena", "-i", INPUT, "-d", str(RUN_DIR / "ito3"),
             "particles/pusher=ito3"]
        )
        assert not testutils.run_command(
            ["./athena", "-i", INPUT, "-d", str(RUN_DIR / "order3"),
             "particles/ito_order=3"]
        )
        assert not testutils.run_command(
            ["./athena", "-i", INPUT, "-d", str(RUN_DIR / "wrong_type"),
             "particles/particle_type=cosmic_ray"]
        )
    finally:
        shutil.rmtree(RUN_DIR, ignore_errors=True)


def test_legacy_v2_and_v3_ito_restarts_infer_covariance_model():
    """Legacy restarts resume without a supplemental input file."""
    root = Path("run_particles_ito_legacy_restart")
    shutil.rmtree(root, ignore_errors=True)
    try:
        root.mkdir(parents=True)
        for version, model in (
            (2, "published_diagonal"),
            (3, "full_finite_step"),
        ):
            source = root / f"v{version}_source"
            resumed = root / f"v{version}_resumed"
            uninterrupted = root / f"v{version}_uninterrupted"
            assert testutils.run(
                INPUT,
                [
                    "-d",
                    str(source),
                    "time/nlim=1",
                    f"particles/ito_covariance_model={model}",
                ],
            )
            restart = sorted(
                (source / "rst/rank_00000000").glob("ito_tracers.*.rst")
            )[-1]
            legacy_restart = root / f"ito_tracers_v{version}.rst"
            make_legacy_restart(restart, legacy_restart, version)

            dry_run = subprocess.run(
                ["./athena", "-r", str(legacy_restart), "-n"],
                check=False,
                capture_output=True,
                text=True,
            )
            assert dry_run.returncode == 0
            assert (
                "ito_covariance_model"
                in dry_run.stdout
                and "inferred_from_restart_payload" in dry_run.stdout
            )
            incompatible_model = (
                "full_finite_step"
                if model == "published_diagonal"
                else "published_diagonal"
            )
            assert not testutils.run_command(
                [
                    "./athena",
                    "-r",
                    str(legacy_restart),
                    "-d",
                    str(root / f"v{version}_incompatible"),
                    f"particles/ito_covariance_model={incompatible_model}",
                    "time/nlim=2",
                ]
            )
            assert testutils.run_command(
                [
                    "./athena",
                    "-r",
                    str(legacy_restart),
                    "-d",
                    str(resumed),
                    "time/nlim=2",
                ]
            )
            assert testutils.run(
                INPUT,
                [
                    "-d",
                    str(uninterrupted),
                    "time/nlim=2",
                    f"particles/ito_covariance_model={model}",
                ],
            )
            resumed_tags, resumed_state = _final_state(resumed)
            uninterrupted_tags, uninterrupted_state = _final_state(uninterrupted)
            np.testing.assert_array_equal(resumed_tags, uninterrupted_tags)
            np.testing.assert_array_equal(resumed_state, uninterrupted_state)
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_ito2_rk3_flux_weights_and_mhd():
    """Final-RK flux weights preserve moments under RK3, and MHD runs Ito-2."""
    shutil.rmtree(RUN_RK3, ignore_errors=True)
    shutil.rmtree(RUN_MHD, ignore_errors=True)
    try:
        assert testutils.run(
            INPUT,
            [
                "-d",
                str(RUN_RK3),
                "time/integrator=rk3",
                "mesh/nx3=8",
                "meshblock/nx3=4",
                "problem/vz0=0.25",
            ],
        )
        _assert_uniform_flow_moments(
            RUN_RK3,
            velocities=(0.5, 0.0, 0.25),
            cell_widths=(1.0 / 32.0, 1.0 / 8.0, 1.0 / 8.0),
        )

        assert testutils.run(
            MHD_INPUT,
            [
                "-d",
                str(RUN_MHD),
                "particles/particle_type=lagrangian_ito",
                "particles/pusher=ito2",
            ],
        )
        history = read_history(
            RUN_MHD
            / "prtcl_thermo_history/"
            "lagrangian_mc_thermo_mhd.prtcl_thermo_history.thp"
        )
        for field in ("bmag", "beta", "alfven_speed", "mach"):
            assert np.all(np.isfinite(history[field]))
    finally:
        shutil.rmtree(RUN_RK3, ignore_errors=True)
        shutil.rmtree(RUN_MHD, ignore_errors=True)
