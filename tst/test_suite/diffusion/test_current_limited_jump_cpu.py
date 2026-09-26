"""Experimental jump estimator: cached-cycle semantics and short dynamic checks."""

import re

import numpy as np
import pytest

from test_suite.diffusion.test_current_limited_cpu import run_model
from test_suite.diffusion.test_sts_rkl2_cpu import cells, read_binary, read_state
from test_suite.diffusion.test_sts_rkl2_mpicpu import AMR, compare_states

BREC_OUTPUT = """
<output3>
file_type=bin
variable=mhd_brec
id=brec
data_precision=real
dt=100
"""


def run_jump(directory, flags=(), **kwargs):
    defaults = ["mhd/b_rec_method=jump", "mhd/b_rec_radius=0.0625",
                "mhd/b_rec_floor=0.01", "mhd/d_i=0.1", "mesh/nghost=6",
                "mesh/nx1=64", "meshblock/nx1=32", "time/tlim=0.02"]
    return run_model(directory, defaults + list(flags), **kwargs)


def test_jump_cache_is_frozen_for_full_sts_cycle(tmp_path):
    result = run_jump(tmp_path, ["time/nlim=1", "time/tlim=10"], mode="sts",
                      extra=BREC_OUTPUT)
    counts = re.search(r"STS sweeps = (\d+) STS stages = (\d+)", result.stdout)
    assert counts and int(counts[1]) == 2 and int(counts[2]) >= 6
    snapshots = sorted((tmp_path / "bin").glob("*.brec.*.bin"))
    before, after = [cells(read_binary(path))[2]["brec"] for path in snapshots]
    np.testing.assert_array_equal(before, after)
    for first, cached in [(True, before), (False, after)]:
        _, _, state = cells(read_state(tmp_path, first=first))
        b = np.array([state[f"bcc{d}"] for d in (1, 2, 3)])
        estimate = np.maximum(np.linalg.norm(np.roll(b, -4, axis=1)
                                             - np.roll(b, 4, axis=1), axis=0) / 2, 0.01)
        if first:
            np.testing.assert_allclose(cached, estimate, rtol=2e-14, atol=2e-14)
        else:
            assert np.max(abs(cached - estimate)) > 1e-5


@pytest.mark.parametrize("amr, nghost", [(False, 4), (True, 8)])
def test_jump_rejects_insufficient_future_ghost_width(tmp_path, amr, nghost):
    # Uniform m=4 needs five cells; potential AMR m=8 needs nine, even before
    # any refinement actually occurs. A valid present mesh is not sufficient.
    result = run_jump(tmp_path, [f"mesh/nghost={nghost}"], check=False,
                      extra=AMR if amr else "")
    assert result.returncode != 0
    assert "ghost" in result.stdout.lower() + result.stderr.lower()


@pytest.mark.parametrize("setting", ["mhd/b_rec_radius=nan", "mhd/b_rec_floor=0"])
def test_jump_rejects_invalid_parameters(tmp_path, setting):
    result = run_jump(tmp_path, [setting], check=False)
    assert result.returncode != 0
    assert "finite positive b_rec_radius and b_rec_floor" in result.stdout + result.stderr


def test_jump_rejects_meshblocks_narrower_than_halo(tmp_path):
    result = run_jump(tmp_path, ["mesh/nghost=10", "meshblock/nx1=8"], check=False)
    assert result.returncode != 0
    assert "meshblock width" in result.stdout + result.stderr


@pytest.mark.parametrize("guide", [0, 1])
def test_jump_matches_fixed_harris_short_evolution(tmp_path, guide):
    flags = ["problem/test=harris", "problem/width=0.03", "problem/sheet_density=1",
             "problem/perturbation_flux=0.001", f"problem/guide_field={guide}",
             "mhd/ohmic_resistivity=0.001", "mhd/eta_max=0.1", "mhd/d_i=0.03",
             "mhd/b_rec_radius=0.125", "mesh/nghost=10", "mesh/nx2=64",
             "meshblock/nx2=32", "time/tlim=0.005"]
    flags += [f"mesh/{side}x{axis}_bc=outflow" for axis in (1, 2) for side in ("i", "o")]
    states = []
    for method in ("constant", "jump"):
        directory = tmp_path / method
        run_jump(directory, flags + [f"mhd/b_rec_method={method}"], mode="sts")
        states.append(cells(read_state(directory))[2])
    for field in states[0]:
        # R > 4 sheet widths measures the upstream field at the central sheet.
        # This is an early-time field comparison, not a steady reconnection-rate test.
        np.testing.assert_allclose(states[0][field], states[1][field], rtol=0, atol=1e-4)


def test_jump_preserves_field_amplitude_scaling(tmp_path):
    states = []
    for b0, method in [(1, "jump"), (0.5, "jump"), (0.5, "constant")]:
        directory = tmp_path / f"{b0}_{method}"
        flags = ["problem/test=sheet", "problem/width=0.03", f"problem/b0={b0}",
                 f"problem/pressure={0.5*b0*b0}", "mhd/d_i=0.03",
                 f"mhd/ohmic_resistivity={0.001*b0}", f"mhd/eta_max={0.1*b0}",
                 f"mhd/b_rec_floor={0.01*b0}", f"mhd/b_rec_method={method}",
                 "mhd/b_rec=1", "mhd/b_rec_radius=0.125", "mesh/nghost=10",
                 f"time/tlim={0.02/b0}", "mesh/ix1_bc=outflow", "mesh/ox1_bc=outflow"]
        run_jump(directory, flags, mode="sts")
        _, _, state = cells(read_state(directory))
        state["ener"] /= b0*b0
        for field in ("mom1", "mom2", "mom3", "bcc1", "bcc2", "bcc3"):
            state[field] /= b0
        states.append(state)
    for field in states[0]:
        np.testing.assert_allclose(states[0][field], states[1][field], rtol=2e-12, atol=2e-12)
    assert np.max(abs(states[2]["bcc2"] - states[1]["bcc2"])) > 1e-4


def test_jump_lag_temporal_convergence(tmp_path):
    run_jump(tmp_path / "reference", ["time/cfl_number=0.02"])
    _, _, reference = cells(read_state(tmp_path / "reference"))
    errors = []
    for ratio in (8, 4, 2):
        directory = tmp_path / str(ratio)
        run_jump(directory, [f"time/sts_max_dt_ratio={ratio}"], mode="sts")
        _, _, state = cells(read_state(directory))
        errors.append(np.max(abs(state["bcc3"] - reference["bcc3"])))
    # Beginning-of-cycle lagging is generally first order; do not require the
    # second-order convergence of the fixed-B_rec RKL2 test here.
    assert errors[1] < 0.75 * errors[0], errors
    assert errors[2] < 0.75 * errors[1], errors


def test_jump_amr_rank_agreement_and_restart(tmp_path):
    flags = ["mesh/nx2=64", "meshblock/nx2=32", "mesh/nghost=10",
             "problem/wave_n2=1", "time/tlim=10", "time/nlim=4"]
    outputs = []
    for ranks in (1, 2):
        directory = tmp_path / str(ranks)
        run_jump(directory, flags, mode="sts", ranks=ranks, extra=AMR)
        initial, final = read_state(directory, first=True), read_state(directory)
        assert final["n_mbs"] > initial["n_mbs"]
        outputs.append(final)
    compare_states(*outputs)
    split, resumed = tmp_path / "split", tmp_path / "resumed"
    restart_output = "\n<output2>\nfile_type=rst\ndcycle=2\n"
    run_jump(split, flags + ["time/nlim=2"], mode="sts", ranks=2,
             extra=AMR + restart_output)
    checkpoints = sorted((split / "rst").rglob("*.rst"))
    checkpoints = [p for p in checkpoints if p.parent.name in ("rst", "rank_00000000")]
    assert checkpoints
    # Use only restart overrides: all model and AMR settings are serialized.
    run_model(resumed, ["time/nlim=4", "output2/dcycle=0"], mode="sts", ranks=2,
              restart=checkpoints[-1])
    compare_states(outputs[1], read_state(resumed))


def test_jump_three_dimensional_cache_rank_agreement(tmp_path):
    flags = ["mesh/nx1=32", "meshblock/nx1=16", "mesh/nx2=16", "meshblock/nx2=16",
             "mesh/nx3=16", "meshblock/nx3=16", "mesh/nghost=4", "problem/wave_n2=1",
             "problem/wave_n3=1", "time/nlim=3", "time/tlim=10"]
    outputs, caches = [], []
    for ranks in (1, 2):
        directory = tmp_path / str(ranks)
        run_jump(directory, flags, mode="sts", ranks=ranks, extra=BREC_OUTPUT)
        outputs.append(read_state(directory))
        cache_file = sorted((directory / "bin").glob("*.brec.*.bin"))[-1]
        caches.append(cells(read_binary(cache_file))[2]["brec"])
        assert np.isfinite(caches[-1]).all() and caches[-1].min() >= 0.01
    compare_states(*outputs)
    np.testing.assert_allclose(*caches, rtol=2e-13, atol=2e-14)
