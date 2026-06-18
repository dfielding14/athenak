"""Regression tests for fully 2D turbulence and coarsened SGS output."""

from pathlib import Path
import math
import subprocess
import sys

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
INPUT = REPO_ROOT / "tst" / "inputs" / "turb_sgs_2d.athinput"
PRODUCTION_INPUTS = (
    REPO_ROOT / "inputs" / "hydro" / "tiegan_sgs" / "mach010_12288_k16_viscous.athinput",
    REPO_ROOT / "inputs" / "hydro" / "tiegan_sgs" / "mach010_16384_k16_viscous.athinput",
)
ATHENA = Path.cwd() / "athena"
sys.path.insert(0, str(REPO_ROOT / "vis" / "python"))

from bin_convert import read_binary, read_coarsened_binary  # noqa: E402


def run_athena(output_dir, *overrides, restart=None):
    """Run the focused input in an isolated output directory."""
    output_dir.mkdir(parents=True, exist_ok=True)
    command = [str(ATHENA), "-d", str(output_dir)]
    if restart is None:
        command.extend(["-i", str(INPUT)])
    else:
        command.extend(["-r", str(restart)])
    command.extend(overrides)
    return subprocess.run(command, capture_output=True, text=True, check=False)


def require_success(result):
    """Report executable output when a regression run fails."""
    assert result.returncode == 0, result.stdout + result.stderr


def latest(path, pattern):
    """Return the lexically latest numbered output file."""
    outputs = sorted(path.glob(pattern))
    assert outputs
    return outputs[-1]


def input_parameter(path, block, name):
    """Read one scalar parameter from an AthenaK input block."""
    current_block = None
    for raw_line in path.read_text().splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if line.startswith("<") and line.endswith(">"):
            current_block = line[1:-1]
            continue
        if current_block == block and "=" in line:
            key, value = (item.strip() for item in line.split("=", 1))
            if key == name:
                return value
    raise KeyError(f"missing <{block}>/{name} in {path}")


def square_mean(values, factor):
    """Average square, non-overlapping boxes in the active x-y plane."""
    nz, ny, nx = values.shape
    return values.reshape(nz, ny // factor, factor, nx // factor, factor).mean(
        axis=(2, 4)
    )


def assemble_2d_blocks(data, names):
    """Assemble uniform-grid MeshBlocks into global two-dimensional arrays."""
    logical = np.asarray(data["mb_logical"])
    sample = np.asarray(data["mb_data"][names[0]][0])[0]
    ny_block, nx_block = sample.shape
    result = {
        name: np.empty(
            ((logical[:, 1].max() + 1) * ny_block,
             (logical[:, 0].max() + 1) * nx_block)
        )
        for name in names
    }
    for block, (lx1, lx2, _lx3, _level) in enumerate(logical):
        j0 = lx2 * ny_block
        i0 = lx1 * nx_block
        for name in names:
            result[name][j0:j0 + ny_block, i0:i0 + nx_block] = np.asarray(
                data["mb_data"][name][block]
            )[0]
    return result


def test_2d_sgs_output_matches_direct_favre_filter(tmp_path):
    """The producer writes final Favre velocities and SGS stresses on a 2D mesh."""
    run_dir = tmp_path / "run"
    require_success(run_athena(run_dir))

    state = read_binary(str(latest(run_dir / "bin", "*.state.*.bin")))
    force = read_binary(str(latest(run_dir / "bin", "*.force.*.bin")))
    sgs = read_coarsened_binary(
        str(latest(run_dir / "cbin_sgs_2", "*.sgs.*.cbin"))
    )

    assert sgs["var_names"] == ["dens", "velx", "vely", "tau_xx", "tau_xy", "tau_yy"]
    assert sgs["Nx3"] == 1
    assert sgs["nx3_mb"] == 1
    np.testing.assert_array_equal(state["mb_logical"], sgs["mb_logical"])

    for block in range(state["n_mbs"]):
        rho = np.asarray(state["mb_data"]["dens"][block])
        mx = np.asarray(state["mb_data"]["mom1"][block])
        my = np.asarray(state["mb_data"]["mom2"][block])
        rho_bar = square_mean(rho, 2)
        mx_bar = square_mean(mx, 2)
        my_bar = square_mean(my, 2)
        expected = {
            "dens": rho_bar,
            "velx": mx_bar / rho_bar,
            "vely": my_bar / rho_bar,
            "tau_xx": square_mean(mx * mx / rho, 2) - mx_bar * mx_bar / rho_bar,
            "tau_xy": square_mean(mx * my / rho, 2) - mx_bar * my_bar / rho_bar,
            "tau_yy": square_mean(my * my / rho, 2) - my_bar * my_bar / rho_bar,
        }
        for name, values in expected.items():
            np.testing.assert_allclose(
                sgs["mb_data"][name][block], values, rtol=5.0e-6, atol=5.0e-8
            )

    assert np.all(np.asarray(force["mb_data"]["force3"]) == 0.0)
    assert np.all(np.asarray(state["mb_data"]["mom3"]) == 0.0)
    history = np.loadtxt(run_dir / "turb_sgs_2d_test.hydro.hst")
    assert np.all(history[:, 5] == 0.0)
    assert np.all(history[:, 8] == 0.0)


def test_isothermal_viscosity_changes_the_evolution(tmp_path):
    """The configured shear viscosity must have a measurable dynamical effect."""
    viscous_dir = tmp_path / "viscous"
    inviscid_dir = tmp_path / "inviscid"
    require_success(run_athena(viscous_dir))
    require_success(run_athena(inviscid_dir, "hydro/viscosity=1.0e-30"))

    viscous = read_binary(str(latest(viscous_dir / "bin", "*.state.*.bin")))
    inviscid = read_binary(str(latest(inviscid_dir / "bin", "*.state.*.bin")))
    momentum_difference = max(
        np.max(
            np.abs(
                np.asarray(viscous["mb_data"][name])
                - np.asarray(inviscid["mb_data"][name])
            )
        )
        for name in ("mom1", "mom2")
    )
    assert momentum_difference > 1.0e-4


def test_2d_parabolic_spectrum_uses_active_dimensions(tmp_path):
    """Inactive x3 spacing cannot broaden a two-dimensional forcing spectrum."""
    run_dir = tmp_path / "spectrum"
    require_success(run_athena(run_dir))
    force = read_binary(str(latest(run_dir / "bin", "*.force.*.bin")))
    global_force = assemble_2d_blocks(force, ("force1", "force2"))
    power = sum(
        np.abs(np.fft.rfft2(global_force[name])) ** 2
        for name in ("force1", "force2")
    )

    peak_power = power[0, 2] + power[2, 0]
    edge_power = power[0, 1] + power[1, 0] + power[0, 3] + power[3, 0]
    assert peak_power > 0.0
    assert edge_power < peak_power * 1.0e-12


def test_sparse_annulus_has_global_isotropic_fourier_support(tmp_path):
    """Sparse annulus forcing uses only its global, angularly balanced mode set."""
    run_dir = tmp_path / "sparse"
    sparse_overrides = (
        "turb_driving/mode_sampling=sparse_annulus",
        "turb_driving/sparse_mode_count=8",
        "turb_driving/nlow=5",
        "turb_driving/nhigh=7",
        "turb_driving/npeak=6",
    )
    require_success(run_athena(run_dir, *sparse_overrides))
    force = read_binary(str(latest(run_dir / "bin", "*.force.*.bin")))
    global_force = assemble_2d_blocks(force, ("force1", "force2"))
    force_fft = {
        name: np.fft.fft2(global_force[name]) for name in ("force1", "force2")
    }
    power = sum(np.abs(values) ** 2 for values in force_fft.values())
    ny, nx = power.shape
    ky = np.fft.fftfreq(ny) * ny
    kx = np.fft.fftfreq(nx) * nx
    kx_grid, ky_grid = np.meshgrid(kx, ky)
    radius = np.sqrt(kx_grid * kx_grid + ky_grid * ky_grid)
    support = (radius > 0.0) & (power > power.max() * 1.0e-12)

    assert np.count_nonzero(support) == 16
    assert np.all((radius[support] >= 5.0) & (radius[support] <= 7.0))
    assert np.any((kx_grid * ky_grid)[support] < 0.0)
    assert np.any((kx_grid * ky_grid)[support] > 0.0)
    divergence = kx_grid * force_fft["force1"] + ky_grid * force_fft["force2"]
    relative_divergence = np.max(np.abs(divergence[support])) / np.max(
        radius[support] * np.sqrt(power[support])
    )
    assert relative_divergence < 1.0e-7
    assert np.max(
        np.abs(global_force["force1"] - np.roll(global_force["force1"], nx // 2, axis=1))
    ) > 1.0e-6

    split_dir = tmp_path / "sparse_split"
    require_success(
        run_athena(
            split_dir,
            *sparse_overrides,
            "time/nlim=1",
            "output4/file_type=rst",
            "output4/dcycle=1",
        )
    )
    restart = latest(split_dir / "rst", "*.rst")
    resume_dir = tmp_path / "sparse_resume"
    require_success(run_athena(resume_dir, "time/nlim=3", restart=restart))
    resumed_force = read_binary(str(latest(resume_dir / "bin", "*.force.*.bin")))
    for name in ("force1", "force2", "force3"):
        np.testing.assert_array_equal(
            force["mb_data"][name], resumed_force["mb_data"][name]
        )


def test_sparse_annulus_construction_does_not_scan_mode_volume(tmp_path):
    """Sparse construction fills a narrow shell without scanning a large mode volume."""
    narrow_result = run_athena(
        tmp_path / "narrow_sparse_mode",
        "time/nlim=0",
        "turb_driving/mode_sampling=sparse_annulus",
        "turb_driving/sparse_mode_count=64",
        "turb_driving/nlow=15",
        "turb_driving/nhigh=17",
        "turb_driving/npeak=16",
    )
    require_success(narrow_result)
    assert "turbulence modes = 64" in narrow_result.stdout

    result = run_athena(
        tmp_path / "large_sparse_mode",
        "time/nlim=0",
        "turb_driving/mode_sampling=sparse_annulus",
        "turb_driving/sparse_mode_count=8",
        "turb_driving/nlow=9999",
        "turb_driving/nhigh=10001",
        "turb_driving/npeak=10000",
    )
    require_success(result)
    assert "turbulence modes = 8" in result.stdout


def test_2d_turbulence_rejects_kz_modes(tmp_path):
    """A nominally 2D mesh cannot silently enable out-of-plane forcing modes."""
    result = run_athena(tmp_path / "bad_kz", "turb_driving/max_kz=1")
    assert result.returncode != 0
    assert "two-dimensional meshes require min_kz = max_kz = 0" in (
        result.stdout + result.stderr
    )


def test_isothermal_hllc_is_rejected(tmp_path):
    """HLLC is not implemented for AthenaK's isothermal hydro system."""
    result = run_athena(tmp_path / "bad_hllc", "hydro/rsolver=hllc", "time/nlim=0")
    assert result.returncode != 0
    assert "hllc cannot be used with isothermal EOS" in result.stdout + result.stderr


def test_production_inputs_have_resolved_explicit_viscosity():
    """The production pair keeps one physical viscosity with a resolved cutoff."""
    expected_cells = {12288: 7.5, 16384: 10.0}
    for path in PRODUCTION_INPUTS:
        resolution = int(input_parameter(path, "mesh", "nx1"))
        assert int(input_parameter(path, "mesh", "nx2")) == resolution
        assert int(input_parameter(path, "mesh", "nx3")) == 1
        assert input_parameter(path, "hydro", "eos") == "isothermal"
        assert input_parameter(path, "hydro", "rsolver") == "roe"

        viscosity = float(input_parameter(path, "hydro", "viscosity"))
        injection = float(input_parameter(path, "turb_driving", "dedt"))
        forcing_mode = float(input_parameter(path, "turb_driving", "npeak"))
        box_length = (
            float(input_parameter(path, "mesh", "x1max"))
            - float(input_parameter(path, "mesh", "x1min"))
        )
        forcing_wavenumber = 2.0 * math.pi * forcing_mode / box_length
        enstrophy_injection = forcing_wavenumber**2 * injection
        viscous_length = (viscosity**3 / enstrophy_injection) ** (1.0 / 6.0)
        dissipation_mode = box_length / (2.0 * math.pi * viscous_length)

        assert viscosity == 8.4e-7
        assert viscous_length * resolution >= expected_cells[resolution]
        assert 15.5 <= dissipation_mode / forcing_mode <= 16.5
