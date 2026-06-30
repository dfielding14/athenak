"""Integrated regressions for simple_TRML, frame tracking, and MC tracers."""

from pathlib import Path
from subprocess import PIPE, run
import sys

import numpy as np
import pytest


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "vis" / "python"))
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import athena_read  # noqa: E402
import bin_convert  # noqa: E402
from read_prtcl_thermo_history import read_history  # noqa: E402


@pytest.fixture(scope="module")
def simple_trml_binary(tmp_path_factory: pytest.TempPathFactory) -> Path:
    root = tmp_path_factory.mktemp("simple_trml_build")
    build = root / "build"
    configure = run(
        [
            "cmake",
            "-S",
            str(REPO_ROOT),
            "-B",
            str(build),
            "-DPROBLEM=simple_TRML",
            "-DCMAKE_BUILD_TYPE=Release",
        ],
        check=False,
        stdout=PIPE,
        stderr=PIPE,
        text=True,
    )
    assert configure.returncode == 0, configure.stdout + configure.stderr
    compile_result = run(
        ["cmake", "--build", str(build), "-j", "4"],
        check=False,
        stdout=PIPE,
        stderr=PIPE,
        text=True,
    )
    assert compile_result.returncode == 0, compile_result.stdout + compile_result.stderr
    return build / "src" / "athena"


def run_case(
    binary: Path,
    input_path: Path,
    run_dir: Path,
    overrides: list[str],
    restart: Path | None = None,
) -> None:
    run_dir.mkdir(exist_ok=True)
    command = [str(binary)]
    if restart is not None:
        command.extend(["-r", str(restart)])
    command.extend(["-i", str(input_path), "-d", str(run_dir), *overrides])
    result = run(command, check=False, stdout=PIPE, stderr=PIPE, text=True)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize("integrator", ["rk2", "rk3", "rk4"])
def test_stage_aware_cooling_matches_energy_loss(
    simple_trml_binary: Path, tmp_path: Path, integrator: str
) -> None:
    basename = f"SimpleTRML{integrator.upper()}"
    run_dir = tmp_path / integrator
    run_case(
        simple_trml_binary,
        REPO_ROOT / "tst" / "inputs" / "simple_trml_cooling.athinput",
        run_dir,
        [f"job/basename={basename}", f"time/integrator={integrator}"],
    )

    hydro = athena_read.hst(str(run_dir / f"{basename}.hydro.hst"))
    check_nan = athena_read.check_nan_flag
    athena_read.check_nan_flag = False
    try:
        # Empty phase bins intentionally report NaN extrema in this uniform test.
        user = athena_read.hst(str(run_dir / f"{basename}.user.hst"))
    finally:
        athena_read.check_nan_flag = check_nan
    assert hydro["time"].size == 2
    cooling_name = next(name for name in user if name.startswith("cooling"))
    assert np.all(np.isfinite(user[cooling_name]))
    elapsed = hydro["time"][-1] - hydro["time"][0]
    energy_loss = hydro["tot-E"][0] - hydro["tot-E"][-1]
    reported_loss = user[cooling_name][-1] * elapsed
    np.testing.assert_allclose(reported_loss, energy_loss, rtol=2.0e-12, atol=2.0e-14)


def test_x3_outflow_and_outer_reservoir_vx_modes(
    simple_trml_binary: Path, tmp_path: Path
) -> None:
    input_path = (
        REPO_ROOT
        / "inputs"
        / "hydro"
        / "TRML"
        / "TRML_with_Tracers_and_Tracking.athinput"
    )
    ng = 4
    nx = 8
    rho_hot = 1.0
    pressure_fixed = 1.0
    shear_velocity = 1.0
    single_precision_atol = 5.0e-7
    zero_gradient_input = tmp_path / "zero_gradient_vx.athinput"
    canonical_text = input_path.read_text()
    assert "ix3_bc = outflow\n" in canonical_text
    assert "ox3_bc = user\n" in canonical_text
    assert "zero_gradient_vx = false\n" in canonical_text
    zero_gradient_input.write_text(
        canonical_text.replace("zero_gradient_vx = false\n", "zero_gradient_vx = true\n")
    )

    for zero_gradient in (False, True):
        mode = "zero_gradient" if zero_gradient else "fixed"
        basename = f"TRMLBoundary{mode}"
        run_dir = tmp_path / mode
        run_case(
            simple_trml_binary,
            zero_gradient_input if zero_gradient else input_path,
            run_dir,
            [
                f"job/basename={basename}",
                "mesh/nx1=8",
                "mesh/nx2=8",
                "mesh/nx3=8",
                "meshblock/nx1=8",
                "meshblock/nx2=8",
                "meshblock/nx3=8",
                "time/nlim=1",
                "time/tlim=1.0",
                "problem/velocity=1.0",
                "problem/phase_sharpness=8",
                "initial_perturbations/nlow=1",
                "initial_perturbations/nhigh=2",
                "initial_perturbations/localization=none",
                "initial_perturbations/x3_scale=-1",
                "frame_tracking/enabled=false",
                "tracer_seed1/count_per_event=1",
                "tracer_seed2/end_time=0",
                "tracer_seed2/cadence=-1",
                "tracer_seed2/count_per_event=1",
                "tracer_seed2/slab_min=0.75",
                "output1/dt=10",
                "output2/dt=10",
                "output3/variable=hydro_u",
                "output3/dt=1.0e-20",
                "output3/ghost_zones=true",
                "output4/dt=10",
                "output5/dt=10",
                "output6/dt=10",
                "output7/dt=10",
                "output8/dt=10",
            ],
        )

        snapshots = sorted((run_dir / "bin").glob("*.bin"))
        assert snapshots
        state = bin_convert.read_binary(str(snapshots[-1]))
        fields = {
            name: np.asarray(values)[0] for name, values in state["mb_data"].items()
        }
        density = fields["dens"]
        vx = fields["mom1"] / density
        kinetic = (
            0.5
            * (fields["mom1"] ** 2 + fields["mom2"] ** 2 + fields["mom3"] ** 2)
            / density
        )
        pressure = (1.666666667 - 1.0) * (fields["ener"] - kinetic)

        interior = slice(ng, ng + nx)
        bottom = (slice(0, ng), interior, interior)
        top = (slice(ng + nx, ng + nx + ng), interior, interior)
        bottom_density = density[ng, interior, interior]
        bottom_pressure = pressure[ng, interior, interior]
        bottom_scalar = fields["r_00"][ng, interior, interior]
        bottom_vx = vx[ng, interior, interior]
        np.testing.assert_allclose(
            density[bottom],
            np.broadcast_to(bottom_density, density[bottom].shape),
            rtol=0.0,
            atol=1.0e-14,
        )
        np.testing.assert_allclose(density[top], rho_hot, rtol=0.0, atol=1.0e-14)
        np.testing.assert_allclose(
            pressure[bottom],
            np.broadcast_to(bottom_pressure, pressure[bottom].shape),
            rtol=0.0,
            atol=2.0e-14,
        )
        np.testing.assert_allclose(
            pressure[top], pressure_fixed, rtol=0.0, atol=single_precision_atol
        )
        np.testing.assert_allclose(
            fields["r_00"][bottom],
            np.broadcast_to(bottom_scalar, fields["r_00"][bottom].shape),
            rtol=0.0,
            atol=1.0e-14,
        )
        np.testing.assert_allclose(fields["r_00"][top], 0.0, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(
            vx[bottom],
            np.broadcast_to(bottom_vx, vx[bottom].shape),
            rtol=0.0,
            atol=single_precision_atol,
        )

        top_active = vx[ng + nx - 1, interior, interior]
        if zero_gradient:
            np.testing.assert_allclose(
                vx[top],
                np.broadcast_to(top_active, vx[top].shape),
                rtol=0.0,
                atol=single_precision_atol,
            )
        else:
            np.testing.assert_allclose(
                vx[top],
                0.5 * shear_velocity,
                rtol=0.0,
                atol=single_precision_atol,
            )


def combined_overrides(basename: str, nlim: int) -> list[str]:
    return [
        f"job/basename={basename}",
        "mesh/nx1=32",
        "mesh/nx2=32",
        "mesh/nx3=32",
        "meshblock/nx1=16",
        "meshblock/nx2=16",
        "meshblock/nx3=16",
        f"time/nlim={nlim}",
        "time/tlim=1.0",
        "problem/phase_sharpness=32",
        "initial_perturbations/nhigh=8",
        "frame_tracking/diagnostic_every=1000",
        "tracer_seed1/count_per_event=48",
        "tracer_seed2/end_time=0.03",
        "tracer_seed2/cadence=0.015",
        "tracer_seed2/count_per_event=16",
        "tracer_seed2/slab_min=0.9375",
        "output1/dt=1.0e-20",
        "output2/dt=1.0e-20",
        "output3/variable=hydro_u",
        "output3/dt=1.0e-20",
        "output4/dt=10",
        "output5/dt=10",
        "output6/dt=10",
        "output7/dt=10",
        "output8/dt=1.0e-20",
    ]


def load_thermo(run_dir: Path) -> dict[str, np.ndarray]:
    files = sorted((run_dir / "prtcl_thermo_history").glob("*.thp"))
    assert len(files) == 1
    return read_history(files[0])


def first_particle_rows(data: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    rows = np.asarray(
        [np.flatnonzero(data["tag"] == tag)[0] for tag in np.unique(data["tag"])]
    )
    order = np.argsort(data["tag"][rows])
    return {name: values[rows][order] for name, values in data.items()}


def assert_seed_populations(data: dict[str, np.ndarray]) -> None:
    first = first_particle_rows(data)
    initial = first["seed_id"] == 1
    injected = first["seed_id"] == 2
    assert np.count_nonzero(initial) == 48
    assert np.count_nonzero(injected) == 48
    np.testing.assert_allclose(first["time"][initial], 0.0, rtol=0.0, atol=0.0)

    # The one-shot population samples the full initial volume, not one phase or face.
    for axis in ("x1", "x2", "x3"):
        assert np.min(first[axis][initial]) < -0.25
        assert np.max(first[axis][initial]) > 0.25

    # Sixteen particles are introduced in the top active-cell layer at each of
    # three uniformly spaced schedule times. The first history sample can lag the
    # requested creation time by at most one timestep.
    expected_times = np.repeat(np.asarray([0.0, 0.015, 0.03]), 16)
    injection_order = np.argsort(first["tag"][injected])
    first_times = first["time"][injected][injection_order]
    assert np.all(first_times >= expected_times - 2.0e-14)
    assert np.all(first_times <= expected_times + 0.01)
    assert np.all(first["x3"][injected] >= 0.9375)
    assert np.all(first["x3"][injected] <= 1.0)


def assert_mc_grid_steps(data: dict[str, np.ndarray]) -> None:
    moved = False
    dx = (1.0 / 32.0, 1.0 / 32.0, 2.0 / 32.0)
    lengths = (1.0, 1.0, 2.0)
    for tag in np.unique(data["tag"]):
        rows = np.flatnonzero(data["tag"] == tag)
        rows = rows[np.argsort(data["cycle"][rows])]
        unique_rows = []
        for cycle in np.unique(data["cycle"][rows]):
            same_cycle = rows[data["cycle"][rows] == cycle]
            for axis in ("x1", "x2", "x3"):
                np.testing.assert_allclose(
                    data[axis][same_cycle], data[axis][same_cycle[0]], rtol=0.0, atol=0.0
                )
            unique_rows.append(same_cycle[0])
        rows = np.asarray(unique_rows)
        deltas = []
        for axis, spacing, length in zip(("x1", "x2", "x3"), dx, lengths):
            delta = np.abs(np.diff(data[axis][rows]))
            if axis != "x3":
                delta = np.minimum(delta, length - delta)
            assert np.all(delta <= spacing + 1.0e-12)
            deltas.append(delta > 1.0e-12)
        if deltas[0].size:
            moved_axes = np.vstack(deltas).sum(axis=0)
            assert np.all(moved_axes <= 1)
            moved = moved or bool(np.any(moved_axes == 1))
    assert moved


def final_particle_rows(data: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    final = data["cycle"] == np.max(data["cycle"])
    final_rows = np.flatnonzero(final)
    rows = np.asarray(
        [final_rows[np.flatnonzero(data["tag"][final_rows] == tag)[-1]]
         for tag in np.unique(data["tag"][final_rows])]
    )
    order = np.argsort(data["tag"][rows])
    return {name: values[rows][order] for name, values in data.items()}


def test_canonical_population_split_includes_both_endpoints(
    simple_trml_binary: Path, tmp_path: Path
) -> None:
    input_path = (
        REPO_ROOT
        / "inputs"
        / "hydro"
        / "TRML"
        / "TRML_with_Tracers_and_Tracking.athinput"
    )
    run_dir = tmp_path / "population"
    run_case(
        simple_trml_binary,
        input_path,
        run_dir,
        [
            "job/basename=CanonicalPopulation",
            "mesh/nx1=8",
            "mesh/nx2=8",
            "mesh/nx3=8",
            "meshblock/nx1=8",
            "meshblock/nx2=8",
            "meshblock/nx3=8",
            "time/nlim=1",
            "time/tlim=1.0",
            "problem/phase_sharpness=8",
            "initial_perturbations/nhigh=4",
            # Compress all 51 canonical events into the first timestep.
            "tracer_seed2/end_time=0.003",
            "tracer_seed2/cadence=0.00006",
            "tracer_seed2/slab_min=0.75",
            "output1/dt=1.0e-20",
            "output2/dt=1.0e-20",
            "output3/dt=10",
            "output4/dt=10",
            "output5/dt=10",
            "output6/dt=10",
            "output7/dt=10",
            "output8/dt=10",
        ],
    )

    first = first_particle_rows(load_thermo(run_dir))
    assert len(first["tag"]) == 3072
    np.testing.assert_array_equal(first["tag"], np.arange(3072))
    initial = first["seed_id"] == 1
    injected = first["seed_id"] == 2
    assert np.count_nonzero(initial) == 777
    assert np.count_nonzero(injected) == 45 * 51
    assert np.all(first["x3"][injected] >= 0.75)
    assert np.all(first["x3"][injected] <= 1.0)


def test_combined_run_and_restart_are_consistent(
    simple_trml_binary: Path, tmp_path: Path
) -> None:
    input_path = (
        REPO_ROOT
        / "inputs"
        / "hydro"
        / "TRML"
        / "TRML_with_Tracers_and_Tracking.athinput"
    )
    continuous = tmp_path / "continuous"
    split = tmp_path / "split"

    run_case(
        simple_trml_binary,
        input_path,
        continuous,
        combined_overrides("CombinedContinuous", 8),
    )
    run_case(
        simple_trml_binary,
        input_path,
        split,
        combined_overrides("CombinedSplit", 4),
    )
    restart_files = sorted(
        (split / "rst" / "rank_00000000").glob("CombinedSplit.*.rst")
    )
    assert restart_files
    run_case(
        simple_trml_binary,
        input_path,
        split,
        combined_overrides("CombinedSplit", 8),
        restart=restart_files[-1].resolve(),
    )

    frame = athena_read.hst(
        str(continuous / "CombinedContinuous.frame_tracker.hst")
    )
    assert np.all(np.isfinite(frame["ft_weight"]))
    assert np.max(frame["ft_weight"]) > 0.0
    assert np.max(frame["ft_misses"]) == 0.0
    assert np.any(np.abs(frame["ft_dv_x3"]) > 0.0)
    assert np.any(np.abs(frame["ft_dx_x3"]) > 0.0)

    continuous_thermo = load_thermo(continuous)
    split_thermo = load_thermo(split)
    assert len(np.unique(continuous_thermo["tag"])) == 96
    assert set(np.unique(continuous_thermo["seed_id"])) == {1, 2}
    assert_seed_populations(continuous_thermo)
    for name, values in continuous_thermo.items():
        if np.issubdtype(values.dtype, np.floating):
            assert np.all(np.isfinite(values)), name
    assert_mc_grid_steps(continuous_thermo)

    # A frame boost changes fluid velocity, not the grid-frame position of an MC tracer.
    # The one-cell/one-axis condition above rejects any extra continuous particle shift.
    final_continuous = final_particle_rows(continuous_thermo)
    final_split = final_particle_rows(split_thermo)
    np.testing.assert_array_equal(final_split["tag"], final_continuous["tag"])
    np.testing.assert_array_equal(final_split["seed_id"], final_continuous["seed_id"])
    for name in continuous_thermo:
        if name in {"time", "cycle", "tag", "seed_id", "gid"}:
            continue
        np.testing.assert_allclose(
            final_split[name], final_continuous[name], rtol=2.0e-13, atol=2.0e-14
        )

    split_frame = athena_read.hst(str(split / "CombinedSplit.frame_tracker.hst"))
    for name in ("ft_vf_x3", "ft_dx_x3", "ft_pos_x3", "ft_err_x3", "ft_dv_x3"):
        np.testing.assert_allclose(
            split_frame[name][-1], frame[name][-1], rtol=2.0e-13, atol=2.0e-14
        )

    continuous_bins = sorted((continuous / "bin").glob("*.bin"))
    split_bins = sorted((split / "bin").glob("*.bin"))
    assert continuous_bins and split_bins
    continuous_state = bin_convert.read_binary(str(continuous_bins[-1]))["mb_data"]
    split_state = bin_convert.read_binary(str(split_bins[-1]))["mb_data"]
    for field in ("dens", "mom1", "mom2", "mom3", "ener", "r_00"):
        np.testing.assert_allclose(
            split_state[field], continuous_state[field], rtol=2.0e-13, atol=2.0e-14
        )


def test_combined_amr_updates_tracker_and_particles(
    simple_trml_binary: Path, tmp_path: Path
) -> None:
    input_path = (
        REPO_ROOT
        / "inputs"
        / "hydro"
        / "TRML"
        / "TRML_with_Tracers_and_Tracking.athinput"
    )
    run_dir = tmp_path / "amr"
    overrides = combined_overrides("CombinedAMR", 4) + [
        "mesh_refinement/refinement=adaptive",
        "tracer_seed1/count_per_event=12",
        "tracer_seed2/end_time=0.0",
        "tracer_seed2/cadence=-1.0",
        "tracer_seed2/count_per_event=12",
    ]
    run_case(simple_trml_binary, input_path, run_dir, overrides)

    frame = athena_read.hst(str(run_dir / "CombinedAMR.frame_tracker.hst"))
    assert np.all(np.isfinite(frame["ft_weight"]))
    assert np.max(frame["ft_weight"]) > 0.0

    thermo = load_thermo(run_dir)
    assert len(np.unique(thermo["tag"])) == 24
    assert set(np.unique(thermo["seed_id"])) == {1, 2}
    assert np.all(np.isfinite(thermo["x1"]))
    assert np.all(np.isfinite(thermo["x2"]))
    assert np.all(np.isfinite(thermo["x3"]))
    assert np.all(thermo["gid"] >= 0)

    snapshots = sorted((run_dir / "bin").glob("*.bin"))
    assert snapshots
    initial_state = bin_convert.read_binary(str(snapshots[0]))
    final_state = bin_convert.read_binary(str(snapshots[-1]))
    initial_levels = initial_state["mb_logical"][:, 3]
    final_levels = final_state["mb_logical"][:, 3]
    assert np.max(final_levels) > np.max(initial_levels)
    assert final_levels.size > initial_levels.size
