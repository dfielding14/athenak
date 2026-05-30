"""Dust-gas stopping-time drag relaxation tests."""

import os
from pathlib import Path
from subprocess import Popen, PIPE

import numpy as np
import pytest

import athena_read
import test_suite.testutils as testutils


def _clean(paths):
    Popen(["rm", "-rf"] + paths, stdout=PIPE, stderr=PIPE).communicate()


def _assert_relaxation_history(fname):
    data = athena_read.hst(fname)
    gas_v = data["gmom1"] / data["gmass"]
    part_v = data["pmom1"] / data["pmass"]
    rel_v = part_v - gas_v
    total_mom = data["gmom1"] + data["pmom1"]

    rel_ratio = abs(rel_v[-1] / rel_v[0])
    if rel_ratio > 0.2:
        pytest.fail(f"drag relaxation too slow: final/initial = {rel_ratio:g}")

    drift = np.max(np.abs(total_mom - total_mom[0]))
    scale = max(abs(total_mom[0]), 1.0)
    if drift / scale > 1.0e-3:
        pytest.fail(f"total momentum drift too large: {drift / scale:g}")
    return data


def _run_expect_fail(flags, expected):
    command = [
        "./athena", "-i", "inputs/tests/particle_drag_relaxation.athinput"
    ] + flags
    process = Popen(command, stdout=PIPE, stderr=PIPE, text=True)
    output, errors = process.communicate()
    combined = output + errors
    if process.returncode == 0:
        pytest.fail(f"invalid input unexpectedly succeeded: {' '.join(flags)}")
    if expected not in combined:
        pytest.fail(
            f"invalid input did not report {expected!r}; output was:\n{combined}"
        )


def test_drag_relaxation():
    """Check drag damps relative velocity while conserving total momentum."""
    try:
        results = testutils.run("inputs/tests/particle_drag_relaxation.athinput")
        assert results, "particle drag relaxation run failed"

        _assert_relaxation_history("particle_drag_relaxation.user.hst")
    finally:
        _clean(["particle_drag_relaxation.user.hst"])


@pytest.mark.parametrize(
    "inputfile, history",
    [
        ("inputs/tests/particle_drag_mhd_relaxation.athinput",
         "particle_drag_mhd_relaxation.user.hst"),
        ("inputs/particles/hello_drag_particles.athinput",
         "hello_drag_particles.user.hst"),
        ("inputs/particles/user_drag_hook.athinput",
         "user_drag_hook.user.hst"),
    ],
)
def test_drag_smoke_inputs(inputfile, history):
    """Run MHD-host, hello-world, and pgen-enrolled user-drag examples."""
    try:
        results = testutils.run(inputfile)
        assert results, f"{inputfile} run failed"
        data = _assert_relaxation_history(history)
        mass_drift = np.max(np.abs(data["pmass"] - data["pmass"][0])) / data["pmass"][0]
        if mass_drift > 1.0e-12:
            pytest.fail(f"particle mass changed in {inputfile}")
    finally:
        _clean([history])


def test_drag_energy_conservation_ideal_gas():
    """Check total gas plus particle kinetic energy when energy feedback is enabled."""
    history = "particle_drag_energy.user.hst"
    try:
        args = [
            "job/basename=particle_drag_energy",
            "time/tlim=0.1",
            "hydro/eos=ideal",
            "hydro/gamma=1.6666666666666667",
            "drag_particles/include_energy=true",
        ]
        results = testutils.run("inputs/tests/particle_drag_relaxation.athinput", args)
        assert results, "ideal-gas drag energy run failed"
        data = athena_read.hst(history)
        total_energy = data["gtotE"] + data["pkinE"]
        drift = np.max(np.abs(total_energy - total_energy[0]))
        scale = max(abs(total_energy[0]), 1.0)
        if drift / scale > 1.0e-4:
            pytest.fail(f"total energy drift too large: {drift / scale:g}")
    finally:
        _clean([history])


def test_drag_particle_outputs():
    """Write particle VTK and tracked-particle files for dust particles."""
    try:
        results = testutils.run("inputs/tests/particle_drag_outputs.athinput")
        assert results, "dust particle output run failed"
        pvtk = list(Path("pvtk").glob("particle_drag_outputs.particles.*.part.vtk"))
        if not pvtk or any(path.stat().st_size == 0 for path in pvtk):
            pytest.fail("dust particle VTK output was not written")
        tracked = Path("trk/particle_drag_outputs.trk")
        if not tracked.exists() or tracked.stat().st_size == 0:
            pytest.fail("tracked dust-particle output was not written")
    finally:
        _clean(["pvtk", "trk"])


def test_drag_restart_user_hook_reenrollment():
    """Restart drag particles and ensure a pgen user drag callback is re-enrolled."""
    try:
        _clean(["rst", "particle_drag_restart.user.hst"])
        results = testutils.run("inputs/tests/particle_drag_restart.athinput")
        assert results, "pre-restart drag run failed"
        rst = "rst/particle_drag_restart.00000.rst"
        if not os.path.exists(rst):
            pytest.fail(f"restart file was not written: {rst}")

        _clean(["particle_drag_restart.user.hst"])
        restarted = testutils.run_command(["./athena", "-r", rst, "time/tlim=0.1"])
        assert restarted, "drag restart run failed"

        data = _assert_relaxation_history("particle_drag_restart.user.hst")
        if data["time"][-1] < 0.1:
            pytest.fail("drag restart did not advance to the requested final time")
        mass_drift = np.max(np.abs(data["pmass"] - data["pmass"][0])) / data["pmass"][0]
        if mass_drift > 1.0e-12:
            pytest.fail(f"particle mass changed across restart: {mass_drift:g}")
    finally:
        _clean(["rst", "particle_drag_restart.user.hst"])


@pytest.mark.parametrize(
    "flags, expected",
    [
        (["drag_particles/enabled=false"], "<particles>/pusher"),
        (["drag_particles/model=bogus"], "<drag_particles>/model"),
        (["drag_particles/model=none"], "<particles>/pusher"),
        (["drag_particles/stopping_time=0.0"], "<drag_particles>/stopping_time"),
        (["drag_particles/cfl_drag=0.0"], "<drag_particles>/cfl_drag"),
        (["particles/cfl_part=0.0"], "<particles>/cfl_part"),
        (["particles/cfl_part=1.1"], "<particles>/cfl_part"),
        (["drag_particles/particle_mass=0.0"], "<drag_particles>/particle_mass"),
        (["drag_particles/interpolation=bogus"], "<drag_particles>/interpolation"),
        (["drag_particles/deposition=bogus"], "<drag_particles>/deposition"),
        (["particles/particle_type=cosmic_ray"], "<particles>/particle_type"),
        (["particles/particle_type=bogus"], "<particles>/particle_type"),
        (["particles/pusher=drift"], "<drag_particles>/enabled"),
        (["particles/pusher=bogus"], "<particles>/pusher"),
        (["particles/ppc=-1.0"], "<particles>/ppc"),
        (["particles/assign_tag=bogus"], "<particles>/assign_tag"),
        (["mesh/ix1_bc=outflow", "mesh/ox1_bc=outflow"], "<mesh>/*_bc"),
    ],
)
def test_drag_invalid_inputs_fail_cleanly(flags, expected):
    """Exercise fatal configuration paths with block/field names in messages."""
    _run_expect_fail(flags, expected)
