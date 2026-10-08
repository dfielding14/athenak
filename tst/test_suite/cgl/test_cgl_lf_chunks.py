"""Focused LF-only RKL2 chunk regressions against an already built executable.

Run directly with pytest; ATHENAK_TEST_EXE may name a built_in_pgens executable
instead of the usual ./athena. No build is performed by this module. Set
--basetemp and pytest's cache_dir to the allocated test workspace.
"""
import json
import os
from pathlib import Path
import re
import subprocess

import numpy as np
import pytest

from test_suite.cgl import cgl_lf_sts_merge as merge


@pytest.fixture(autouse=True)
def configured_binary(monkeypatch):
    executable = Path(os.environ.get("ATHENAK_TEST_EXE", "athena")).absolute()
    assert executable.is_file(), f"Build built_in_pgens first: {executable}"
    assert executable.name == "athena", "Existing merge fixture expects an athena executable"
    monkeypatch.chdir(executable.parent)


def history(directory):
    path = next(directory.glob("*.mhd.hst"))
    header = next(line for line in path.read_text().splitlines() if "[1]=" in line)
    names = re.findall(r"\[\d+\]=(\S+)", header)
    values = np.atleast_2d(np.loadtxt(path))
    assert np.all(np.isfinite(values))
    return dict(zip(names, values.T))


def chunks(directory):
    match = re.search(
        r"CGL LF chunks: max_ratio=(\S+) logical_sweeps=(\d+) "
        r"chunks=(\d+) rhs_evaluations=(\d+)",
        (directory / "stdout.log").read_text())
    assert match, "Missing LF chunk execution summary"
    ratio, sweeps, count, rhs = map(float, match.groups())
    assert ratio > 0 and count > sweeps > 0, match.group(0)
    assert rhs >= 3 * count, match.group(0)
    return {"ratio": ratio, "logical_sweeps": sweeps, "chunks": count,
            "rhs_evaluations": rhs}


def test_default_absent_and_zero_are_identical(tmp_path):
    absent, zero = tmp_path / "absent", tmp_path / "zero"
    merge.run(absent)
    merge.run(zero, {"time/cgl_lf_max_chunk_ratio": 0})
    merge.equal_fields_and_history(absent, zero)
    for directory in (absent, zero):
        assert "CGL LF chunks:" not in (directory / "stdout.log").read_text()


def test_uniform_collision_and_limiter_cadence(tmp_path):
    # Uniform state makes LF transport identically zero. Background exponential
    # decay followed by one backward-Euler mirror relaxation is independently
    # predictable; applying the latter per chunk changes the answer.
    updates = {
        "problem/test_mode": "merge_cfl_probe", "problem/amp": 0,
        "problem/by_amp": 0, "problem/ppar0": 1, "problem/pperp0": 1.8,
        "mhd/nu_coll": 2, "mhd/mirror_limiter": "true",
        "mhd/limiter_nu_coll": 10, "mhd/mirror_threshold": 1,
        "mhd/backup_limiters": "false", "mhd/cgl_lf_strict_admissibility": "true",
        "time/sts_max_dt_ratio": 16, "time/nlim": 3, "time/tlim": .05,
        "output1/dcycle": 1,
    }
    off, on = tmp_path / "off", tmp_path / "on"
    merge.run(off, updates)
    merge.run(on, updates | {"time/cgl_lf_max_chunk_ratio": .5})
    stats = chunks(on)
    left, right = history(off), history(on)
    assert len(right["time"]) >= 3
    np.testing.assert_allclose(right["time"], left["time"], rtol=0, atol=2e-15)
    np.testing.assert_allclose(merge.clean_conserved(on), merge.clean_conserved(off),
                               rtol=0, atol=3e-12)
    assert right["lf_nstage"][-1] > left["lf_nstage"][-1]
    ratio = np.exp(right["aam-D"] / right["mass"])
    piso = (1 + 2 * 1.8) / 3
    measured = 3 * piso * (ratio - 1) / (1 + 2 * ratio)
    expected = [0.8]
    for dt in np.diff(right["time"]):
        delta = expected[-1] * np.exp(-2 * dt)
        expected.append(.5 + (delta - .5) / (1 + 10 * dt) if delta > .5 else delta)
    np.testing.assert_allclose(measured, expected, rtol=2e-11, atol=2e-13)
    (tmp_path / "collision-cadence.json").write_text(json.dumps({
        "time": right["time"].tolist(), "measured_delta_p": measured.tolist(),
        "expected_delta_p": expected, "chunks": stats,
    }, indent=2) + "\n")


def test_smooth_lf_temporal_refinement_and_conservation(tmp_path):
    # The outer timestep is held fixed. Ratios <=2 all use three-stage RKL2,
    # so reducing the chunk size refines the same LF integrator. Compare on the
    # same grid to remove the spatial truncation error from the temporal test.
    states, counts, times = {}, {}, []
    for ratio in (2, 1, .5, .125):
        directory = tmp_path / f"ratio-{ratio}"
        merge.run(directory, {
            "time/cgl_lf_max_chunk_ratio": ratio, "time/sts_max_dt_ratio": 16,
            "time/tlim": .05, "mhd/cgl_lf_strict_admissibility": "true",
            "output1/dcycle": 1,
        })
        states[ratio] = merge.clean_conserved(directory)
        counts[ratio] = chunks(directory)
        times.append(history(directory)["time"])
    for timeline in times[1:]:
        np.testing.assert_allclose(timeline, times[0], rtol=0, atol=2e-15)
    reference = states[.125][:, 7:9]
    errors = np.array([np.sqrt(np.mean((states[r][:, 7:9] - reference)**2))
                       for r in (2, 1, .5)])
    orders = np.log2(errors[:-1] / errors[1:])
    (tmp_path / "temporal-refinement.json").write_text(json.dumps({
        "chunk_ratios": [2, 1, .5], "reference_ratio": .125,
        "rms_pressure_errors": errors.tolist(), "orders": orders.tolist(),
        "execution_counts": counts,
    }, indent=2) + "\n")
    assert np.all(errors > 1e-14), errors
    assert np.all((orders > 1.7) & (orders < 2.4)), (errors, orders)


def test_chunked_restart_matches_uninterrupted_run(tmp_path):
    first, resumed, direct = (tmp_path / name for name in ("first", "resumed", "direct"))
    # Stop at a full cycle boundary, not by clipping dt to an intermediate tlim.
    updates = {"time/cgl_lf_max_chunk_ratio": 1, "time/sts_max_dt_ratio": 16,
               "time/tlim": 1, "time/nlim": 3,
               "mhd/cgl_lf_strict_admissibility": "true", "output1/dcycle": 1}
    merge.run(first, updates)
    restart = sorted((first / "rst").glob("*.rst"))[-1]
    # The chunk parameter is retained by restart; deliberately do not override it.
    merge.run(resumed, {"time/nlim": 6}, restart=restart)
    merge.run(direct, updates | {"time/nlim": 6})
    chunks(resumed)
    np.testing.assert_allclose(merge.clean_conserved(resumed), merge.clean_conserved(direct),
                               rtol=0, atol=3e-13)
    resumed_h, direct_h = history(resumed), history(direct)
    assert resumed_h.keys() == direct_h.keys()
    for name in resumed_h:
        np.testing.assert_allclose(resumed_h[name][-1], direct_h[name][-1],
                                   rtol=2e-12, atol=3e-12, err_msg=name)


def test_invalid_chunk_configurations_fail_before_evolution(tmp_path):
    cases = {
        "negative": {"time/cgl_lf_max_chunk_ratio": -1},
        "nan": {"time/cgl_lf_max_chunk_ratio": "nan"},
        "infinite": {"time/cgl_lf_max_chunk_ratio": "inf"},
        "merged": {"time/sts_merge_half_sweeps": "true"},
        "explicit": {"time/sts_integrator": "none",
                     "mhd/cgl_heat_flux_integrator": "explicit"},
        "outflow": {"mesh/ix1_bc": "outflow", "mesh/ox1_bc": "outflow"},
    }
    for name, changes in cases.items():
        directory = tmp_path / name
        directory.mkdir()
        input_path = directory / "input.athinput"
        input_path.write_text(merge._input({"time/cgl_lf_max_chunk_ratio": 1} | changes))
        env = {key: value for key, value in os.environ.items()
               if not key.startswith("ATHENAK_CGL_LF_")}
        result = subprocess.run([str(Path("athena").absolute()), "-i", str(input_path),
                                 "-d", str(directory)], env=env, text=True,
                                capture_output=True, check=False)
        log = result.stdout + result.stderr
        (directory / "stdout.log").write_text(log)
        assert result.returncode != 0, name
        assert "cgl_lf_max_chunk_ratio" in log, (name, log[-3000:])
        assert not list(directory.glob("*.mhd.hst")), name
