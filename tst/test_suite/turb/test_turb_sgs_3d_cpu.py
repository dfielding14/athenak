"""Regression test for 3D coarsened SGS output."""

from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest


REPO_ROOT = Path(__file__).resolve().parents[3]
INPUT = REPO_ROOT / "tst" / "inputs" / "turb_sgs_3d.athinput"
MHD_INPUT = REPO_ROOT / "tst" / "inputs" / "turb_mhd_sgs_3d.athinput"
ATHENA = Path.cwd() / "athena"
SGS_NAMES = [
    "dens",
    "velx",
    "vely",
    "velz",
    "tau_xx",
    "tau_xy",
    "tau_xz",
    "tau_yy",
    "tau_yz",
    "tau_zz",
]
sys.path.insert(0, str(REPO_ROOT / "vis" / "python"))

from bin_convert import read_binary, read_coarsened_binary  # noqa: E402


def run_athena(output_dir, *overrides, input_file=INPUT, restart=None, launcher=()):
    """Run the focused input in an isolated output directory."""
    output_dir.mkdir(parents=True, exist_ok=True)
    source = ["-r", str(restart)] if restart else ["-i", str(input_file)]
    return subprocess.run(
        [*launcher, str(ATHENA), "-d", str(output_dir), *source, *overrides],
        capture_output=True,
        text=True,
        check=False,
    )


def require_success(result):
    """Report executable output when a regression run fails."""
    assert result.returncode == 0, result.stdout + result.stderr


def latest(path, pattern):
    """Return the lexically latest numbered output file."""
    outputs = sorted(path.glob(pattern))
    assert outputs
    return outputs[-1]


def cube_mean(values, factor):
    """Average non-overlapping boxes in all active dimensions."""
    nz, ny, nx = values.shape
    fz, fy, fx = [factor if size > 1 else 1 for size in (nz, ny, nx)]
    return values.reshape(
        nz // fz, fz, ny // fy, fy, nx // fx, fx
    ).mean(axis=(1, 3, 5))


def expected_sgs_fields(state, block, factor):
    """Compute the expected coarsened Favre state directly from hydro_u."""
    rho, mx, my, mz = [np.asarray(state["mb_data"][name][block], dtype=np.float64)
                       for name in ("dens", "mom1", "mom2", "mom3")]

    rho_bar = cube_mean(rho, factor)
    mx_bar = cube_mean(mx, factor)
    my_bar = cube_mean(my, factor)
    mz_bar = cube_mean(mz, factor)

    return {
        "dens": rho_bar,
        "velx": mx_bar / rho_bar,
        "vely": my_bar / rho_bar,
        "velz": mz_bar / rho_bar,
        "tau_xx": cube_mean(mx * mx / rho, factor)
        - mx_bar * mx_bar / rho_bar,
        "tau_xy": cube_mean(mx * my / rho, factor)
        - mx_bar * my_bar / rho_bar,
        "tau_xz": cube_mean(mx * mz / rho, factor)
        - mx_bar * mz_bar / rho_bar,
        "tau_yy": cube_mean(my * my / rho, factor)
        - my_bar * my_bar / rho_bar,
        "tau_yz": cube_mean(my * mz / rho, factor)
        - my_bar * mz_bar / rho_bar,
        "tau_zz": cube_mean(mz * mz / rho, factor)
        - mz_bar * mz_bar / rho_bar,
    }


def expected_mhd_fields(state, block, factor):
    """All 59 legacy raw moments, including conserved total energy, from fine state."""
    rho, mx, my, mz, energy, bx, by, bz = [
        np.asarray(state["mb_data"][name][block], dtype=np.float64)
        for name in ("dens", "mom1", "mom2", "mom3", "ener", "bcc1", "bcc2", "bcc3")
    ]
    momentum, magnetic = (mx, my, mz), (bx, by, bz)
    pairs = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))
    fields = [rho, mx, my, mz, energy, bx, by, bz]
    fields += [momentum[i] * momentum[j] / rho for i, j in pairs]
    fields += [magnetic[i] * magnetic[j] for i, j in pairs]
    fields += [m * b / rho for m in momentum for b in magnetic]
    fields += [m * energy / rho for m in momentum]
    fields += [mi * mj**2 / rho**2 for mi in momentum for mj in momentum]
    fields += [m * b**2 / rho for m in momentum for b in magnetic]
    fields += [m * bi * bj / rho for m, bi in zip(momentum, magnetic) for bj in magnetic]
    assert len(fields) == 59
    return {f"mhd_sgs_{i + 1}": cube_mean(field, factor)
            for i, field in enumerate(fields)}


def check_sgs(run_dir, factor, mhd=False):
    """Compare every field to a float64 reference formed from the float32 fine dump."""

    state = read_binary(str(latest(run_dir / "bin", "*.state.*.bin")))
    sgs = read_coarsened_binary(
        str(latest(run_dir / f"cbin_sgs_{factor}", "*.sgs.*.cbin"))
    )

    names = [f"mhd_sgs_{i}" for i in range(1, 60)] if mhd else SGS_NAMES
    assert sgs["var_names"] == names
    assert sgs["n_mbs"] == state["n_mbs"]
    for axis in (1, 2, 3):
        active_factor = factor if state[f"Nx{axis}"] > 1 else 1
        assert sgs[f"Nx{axis}"] == state[f"Nx{axis}"] // active_factor
        assert sgs[f"nx{axis}_mb"] == state[f"nx{axis}_mb"] // active_factor
    assert sgs["cycle"] == state["cycle"]
    assert sgs["time"] == pytest.approx(state["time"], rel=5.0e-6)
    np.testing.assert_array_equal(state["mb_logical"], sgs["mb_logical"])

    for block in range(state["n_mbs"]):
        reference = expected_mhd_fields if mhd else expected_sgs_fields
        expected = reference(state, block, factor)
        fine_scale = reference(state, block, 1)
        rho = np.asarray(state["mb_data"]["dens"][block], dtype=np.float64)
        assert np.ptp(rho) > 1.0e-5
        kinetic_scale = sum(np.asarray(state["mb_data"][f"mom{i}"][block],
                                       dtype=np.float64)**2 for i in (1, 2, 3)) / rho
        for name in names:
            tolerance = (5.0e-7 * np.max(kinetic_scale) if name.startswith("tau_")
                         else 5.0e-7 * np.max(np.abs(fine_scale[name])) + 1.0e-14)
            np.testing.assert_allclose(
                sgs["mb_data"][name][block],
                expected[name],
                rtol=5.0e-6,
                atol=tolerance,
            )
    return state, sgs


@pytest.mark.parametrize("block_width,factor", [(8, 1), (8, 2), (8, 8),
                                               (12, 12), (16, 16)])
@pytest.mark.parametrize("mhd", [False, True])
def test_3d_sgs_output_matches_direct_filter(tmp_path, block_width, factor, mhd):
    """Cover single chunks, full blocks, multiple chunks, and an incomplete tail."""
    overrides = [f"mesh/nx1={2 * block_width}", f"output2/coarsen_factor={factor}",
                 "time/nlim=8", "time/tlim=1", "turb_driving/sol_fraction=0.3",
                 "turb_driving/accel_rms=2", "output1/dcycle=8", "output2/dcycle=8"]
    overrides += [f"mesh/nx{axis}={block_width}" for axis in (2, 3)]
    overrides += [f"meshblock/nx{axis}={block_width}" for axis in (1, 2, 3)]
    require_success(run_athena(tmp_path, *overrides,
                               input_file=MHD_INPUT if mhd else INPUT))
    check_sgs(tmp_path, factor, mhd)


def test_ideal_hydro_3d_sgs_matches_direct_filter(tmp_path):
    input_file = tmp_path / "ideal.athinput"
    input_file.write_text(INPUT.read_text() +
                         "\n<hydro>\neos = ideal\ngamma = 1.4\nrsolver = hllc\n")
    require_success(run_athena(
        tmp_path, "time/nlim=8", "time/tlim=1", "turb_driving/sol_fraction=0.3",
        "turb_driving/accel_rms=2", "output1/dcycle=8", "output2/dcycle=8",
        input_file=input_file,
    ))
    check_sgs(tmp_path, 2)


def test_2d_mhd_sgs_output_matches_direct_filter(tmp_path):
    input_file = tmp_path / "2d.athinput"
    input_file.write_text(MHD_INPUT.read_text() +
                         "\n<turb_driving>\nmin_kz = 0\nmax_kz = 0\n")
    require_success(run_athena(
        tmp_path, "mesh/nx1=16", "mesh/nx2=16", "mesh/nx3=1",
        "meshblock/nx3=1", "output2/coarsen_factor=8", input_file=input_file,
    ))
    check_sgs(tmp_path, 8, mhd=True)


def test_mhd_sgs_binary_and_higher_moments_preserve_raw_contract(tmp_path):
    """The ordinary bin59 and cbin236 products retain every raw-field definition."""
    input_file = tmp_path / "moments.athinput"
    input_file.write_text(MHD_INPUT.read_text() +
        "\n<output2>\ncompute_moments = true\n"
        "\n<output3>\nfile_type = bin\nvariable = mhd_sgs\nid = raw\ndcycle = 8\n")
    require_success(run_athena(tmp_path, input_file=input_file))
    state = read_binary(str(latest(tmp_path / "bin", "*.state.*.bin")))
    raw = read_binary(str(latest(tmp_path / "bin", "*.raw.*.bin")))
    moments = read_coarsened_binary(str(latest(tmp_path / "cbin_sgs_2", "*.cbin")))
    assert raw["var_names"] == [f"mhd_sgs_{i}" for i in range(1, 60)]
    assert len(moments["var_names"]) == 236
    assert moments["number_of_moments"] == 4
    for block in range(state["n_mbs"]):
        for name, values in expected_mhd_fields(state, block, 1).items():
            np.testing.assert_allclose(raw["mb_data"][name][block], values,
                                       rtol=5.0e-6, atol=5.0e-7 * np.max(np.abs(values)))
            for power, suffix in enumerate(("1st", "2nd", "3rd", "4th"), start=1):
                np.testing.assert_allclose(
                    moments["mb_data"][f"{name}_{suffix}"][block],
                    cube_mean(values**power, 2), rtol=5.0e-6,
                    atol=5.0e-7 * np.max(np.abs(values)**power),
                )


@pytest.mark.parametrize("scalars", [0, 1])
@pytest.mark.parametrize("file_type", ["bin", "cbin"])
def test_mhd_sgs_rejects_isothermal_even_with_scalar_slot(tmp_path, scalars, file_type):
    input_file = tmp_path / "isothermal.athinput"
    input_file.write_text(MHD_INPUT.read_text() +
        f"\n<mhd>\neos = isothermal\niso_sound_speed = 1\nnscalars = {scalars}\n")
    result = run_athena(tmp_path, f"output2/file_type={file_type}", input_file=input_file)
    assert result.returncode != 0
    assert "mhd_sgs" in result.stdout + result.stderr
    assert "ideal" in result.stdout + result.stderr
