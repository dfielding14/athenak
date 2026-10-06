"""Shared application checks for the narrow optional LF sweep-merging path."""
from pathlib import Path
import json
import os
import re
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
INPUT = ROOT / "inputs/tests/cgl_lf_sts_merge.athinput"


def _input(updates):
    text = INPUT.read_text()
    for key, value in updates.items():
        block, name = key.split("/", 1)
        match = re.search(r"(^<" + re.escape(block) + r">\s*\n)(.*?)(?=^<|\Z)",
                          text, re.M | re.S)
        assert match, key
        body = match.group(2)
        line = re.compile(r"^" + re.escape(name) + r"\s*=.*$", re.M)
        if value is None:
            body = line.sub("", body)
        else:
            body = (line.sub(name + " = " + str(value), body) if line.search(body)
                    else body + "\n" + name + " = " + str(value) + "\n")
        text = text[:match.start(2)] + body + text[match.end(2):]
    return text


def run(directory, updates=None, *, launcher=(), restart=None, wall=None):
    directory.mkdir(parents=True)
    updates = updates or {}
    command = [*launcher, str(Path("athena").absolute())]
    if restart is None:
        staged = directory / "input.athinput"
        staged.write_text(_input(updates))
        command += ["-i", str(staged)]
    else:
        command += ["-r", str(restart)]
        command += [key + "=" + str(value) for key, value in updates.items()]
    command += ["-d", str(directory)]
    if wall:
        command += ["-t", wall]
    env = os.environ.copy()
    for key in list(env):
        if key.startswith("ATHENAK_CGL_LF_"):
            env.pop(key)
    result = subprocess.run(command, env=env, capture_output=True, text=True)
    stdout = result.stdout + result.stderr
    (directory / "stdout.log").write_text(stdout)
    assert result.returncode == 0, stdout[-5000:]
    stats = {key: float(value) for key, value in re.findall(
        r"(accepted|rejected|rejected_cfl|rejected_admissibility|accepted_stages|"
        r"rejected_stages|snapshot_seconds|deferred|consumed|flushed|pending)="
        r"([+\-.0-9eE]+)", stdout)}
    if stats:
        assert stats["pending"] == 0.0
        np.testing.assert_allclose(stats["deferred"],
                                   stats["consumed"] + stats["flushed"],
                                   rtol=3e-5, atol=1e-14)
    return stats


def table(directory):
    return np.loadtxt(sorted((directory / "tab").glob("*.tab"))[-1])


def clean_conserved(directory):
    path = next(directory.glob("*.mhd.hst"))
    lines = path.read_text().splitlines()
    names = re.findall(r"\[\d+\]=(\S+)", next(x for x in lines if "[1]=" in x))
    values = np.atleast_2d(np.loadtxt(path))
    assert np.all(np.isfinite(values))
    columns = dict(zip(names, values.T))
    for name in ("lf_dfloor", "lf_pfloor", "lf_nonfin", "lf_nonpos", "lf_hardbd",
                 "lf_hwproj"):
        assert np.all(columns[name] == 0.0), name
    for name in ("mass", "tot-E"):
        np.testing.assert_allclose(columns[name], columns[name][0], rtol=0, atol=3e-12)
    state = table(directory)
    assert np.all(np.isfinite(state))
    assert np.min(state[:, [3, 7, 8]]) > 0.0
    return state


def equal_fields_and_history(first, second):
    # Parameter headers differ intentionally (option absent/false or false/true).
    # Compare physical field/history bytes; restart synchronization is exercised
    # separately by continuation and wall-clock tests.
    def products(directory):
        return {str(p.relative_to(directory)): p.read_bytes()
                for p in directory.rglob("*") if p.suffix in (".tab", ".hst")}
    left, right = products(first), products(second)
    assert left and left == right


def default_off(tmp_path):
    left, right = tmp_path / "absent", tmp_path / "false"
    assert not run(left)
    assert not run(right, {"time/sts_merge_half_sweeps": "false"})
    equal_fields_and_history(left, right)


def fallback(tmp_path, reason):
    physics = {
        "collision": {"mhd/nu_coll": 10},
        "limiter": {"mhd/mirror_limiter": "true", "mhd/limiter_nu_coll": 0},
        "strict": {"mhd/cgl_lf_strict_admissibility": "true"},
        "explicit": {"time/sts_integrator": "none", "mhd/cgl_heat_flux_integrator": "explicit"},
    }[reason]
    left, right = tmp_path / "off", tmp_path / "on"
    assert not run(left, physics | {"time/sts_merge_half_sweeps": "false"})
    assert not run(right, physics | {"time/sts_merge_half_sweeps": "true"})
    assert "ordinary half-sweep fallback" in (right / "stdout.log").read_text()
    equal_fields_and_history(left, right)


def smooth_and_output_barriers(tmp_path):
    for name, outputs in (
        ("sparse", {}),
        ("every", {"output1/dcycle": 1, "output2/dcycle": 1, "output3/dcycle": 1}),
        ("mixed", {"output1/dcycle": 3, "output2/dcycle": None, "output2/dt": .006}),
    ):
        off, on = tmp_path / (name + "off"), tmp_path / (name + "on")
        run(off, outputs | {"time/sts_merge_half_sweeps": "false"})
        stats = run(on, outputs | {"time/sts_merge_half_sweeps": "true"})
        assert (stats["accepted"] == 0) if name == "every" else (stats["accepted"] > 0)
        if name == "every":
            equal_fields_and_history(off, on)
        else:
            for path in sorted((off / "tab").glob("*.tab")):
                np.testing.assert_allclose(np.loadtxt(path), np.loadtxt(on / "tab" / path.name),
                                           rtol=0, atol=2e-8)
        clean_conserved(on)
        if name == "mixed":
            times = [float(re.search(r"time=\s*([+\-.0-9eE]+)", p.read_text()).group(1))
                     for p in sorted((on / "tab").glob("*.tab"))]
            assert any(0.0 < time < .01 for time in times), times


def restart_sync(tmp_path):
    first, resumed, reference = (tmp_path / x for x in ("first", "resumed", "reference"))
    run(first, {"time/sts_merge_half_sweeps": "true"})
    restart = sorted((first / "rst").glob("*.rst"))[-1]
    stats = run(resumed, {"time/sts_merge_half_sweeps": "true", "time/tlim": .015},
                restart=restart)
    assert stats["accepted"] > 0
    run(reference, {"time/sts_merge_half_sweeps": "false", "time/tlim": .015})
    np.testing.assert_allclose(clean_conserved(resumed), clean_conserved(reference),
                               rtol=0, atol=2e-8)
    before = np.atleast_2d(np.loadtxt(next(first.glob("*.mhd.hst"))))[-1]
    after = np.atleast_2d(np.loadtxt(next(resumed.glob("*.mhd.hst"))))[-1]
    np.testing.assert_allclose(after[[2, 6]], before[[2, 6]], rtol=0, atol=3e-12)


def dynamic_rejection(tmp_path, launchers=((),)):
    states, decisions = [], []
    for number, launcher in enumerate(launchers):
        directory = tmp_path / str(number)
        stats = run(directory, {
            "time/sts_merge_half_sweeps": "true", "meshblock/nx1": 16,
            "time/evolution": "dynamic", "mhd/rsolver": "hlle",
            "problem/test_mode": "merge_cfl_probe", "problem/b0": .1,
            "problem/by_amp": .3, "problem/ppar0": 1, "problem/pperp0": 1.3,
            "time/sts_max_dt_ratio": -1, "time/tlim": .02, "time/nlim": 12,
        }, launcher=launcher)
        assert stats["rejected_cfl"] > 0 and stats["rejected_admissibility"] == 0
        states.append(clean_conserved(directory))
        decisions.append({key: stats[key] for key in
                          ("accepted", "rejected_cfl", "rejected_stages", "pending")})
    for state, decision in zip(states[1:], decisions[1:]):
        np.testing.assert_allclose(state, states[0], rtol=2e-12, atol=2e-13)
        assert decision == decisions[0]


def wall_sync(tmp_path):
    directory = tmp_path / "wall"
    stats = run(directory, {
        "time/sts_merge_half_sweeps": "true", "time/tlim": 100, "time/nlim": -1,
        "problem/test_mode": "merge_cfl_probe", "problem/b0": 1,
        "problem/by_amp": 0, "problem/ppar0": 1, "problem/pperp0": 1,
    }, wall="00:00:01")
    assert stats["accepted"] > 0 and stats["flushed"] > 0
    assert "Terminating on wall clock limit" in (directory / "stdout.log").read_text()
    clean_conserved(directory)
    assert list((directory / "rst").glob("*.rst"))


def temporal_refinement(tmp_path):
    chi = np.sqrt(8 / np.pi) / (2 * np.pi)
    target = 256 * .4 * .5 / (64**2 * chi)
    errors = []
    for ratio in (2, 1, .5):
        states = []
        for enabled in (False, True):
            directory = tmp_path / f"{ratio}-{enabled}"
            stats = run(directory, {"time/sts_merge_half_sweeps": str(enabled).lower(),
                                   "time/sts_max_dt_ratio": ratio, "time/tlim": target})
            if enabled:
                assert stats["accepted"] > 0
            states.append(clean_conserved(directory))
        errors.append(np.sqrt(np.mean((states[0][:, 7] - states[1][:, 7])**2)))
    orders = np.log2(np.asarray(errors[:-1]) / errors[1:])
    (tmp_path / "temporal-refinement.json").write_text(json.dumps({
        "time": target, "ratios": [2, 1, .5], "rms_temperature_errors": errors,
        "orders": orders.tolist()}, indent=2) + "\n")
    assert np.all((orders > 1.7) & (orders < 2.3)), (errors, orders)
