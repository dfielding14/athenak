"""LF corner states and shared-face conservation on one and four MPI ranks."""

import numpy as np
import pytest

from test_suite import testutils
from test_suite.cgl import cgl_lf_outflow as outflow
from test_suite.cgl import cgl_lf_smr_conservation as conservation
from test_suite.cgl import test_cgl_amr_gpu as amr


@pytest.mark.parametrize("nlim", [0, 4])
@pytest.mark.parametrize("primitive", [False, True])
def test_cgl_lf_smr_outflow(tmp_path, nlim, primitive):
    states = []
    for ranks in (1, 4):
        directory = tmp_path / f"ranks{ranks}"
        directory.mkdir()
        flags = outflow.flags(directory, nlim) + [
            f"mesh_refinement/prolong_primitives={str(primitive).lower()}"]
        testutils.mpi_run(outflow.INPUT, flags, threads=ranks)
        states.append(outflow.check(directory, nlim))
    for variable in states[0]["var_names"]:
        np.testing.assert_array_equal(states[0]["mb_data"][variable],
                                      states[1]["mb_data"][variable])


@pytest.mark.parametrize("primitive", [False, True])
def test_cgl_lf_smr_outflow_nonuniform(tmp_path, primitive):
    states = []
    repairs = []
    for ranks in (1, 4):
        directory = tmp_path / f"ranks{ranks}"
        directory.mkdir()
        flags = outflow.flags(directory, 4) + ["problem/density_amp=0.2",
                 "problem/scalar_amp=0.1", "problem/field_amp=0.02",
                 f"mesh_refinement/prolong_primitives={str(primitive).lower()}"]
        testutils.mpi_run(outflow.INPUT, flags, threads=ranks)
        state = outflow.bin_convert.read_binary(str(sorted(
            (directory / "bin").glob(f"{outflow.BASENAME}.state.*.bin"))[-1]))
        assert state["cycle"] == 4 and state["n_mbs"] == 25
        for variable in state["var_names"]:
            assert np.all(np.isfinite(state["mb_data"][variable]))
        for variable in ("dens", "eint", "p_perp"):
            assert np.min(state["mb_data"][variable]) > 0
        amr._assert_clean_lf(testutils.athena_read.hst(
            str(directory / f"{outflow.BASENAME}.mhd.hst")))
        user = testutils.athena_read.hst(str(directory / f"{outflow.BASENAME}.user.hst"))
        amr._assert_clean_user(user, max_ndiv=1e-12)
        # Physical-edge coarse scratch is repaired before valid fine ghosts are refilled.
        repairs.append(np.array([user[name] for name in amr.AMR_REPAIR_COLUMNS]))
        states.append(state)
    amr._assert_same_state(*states)
    np.testing.assert_array_equal(*repairs)


@pytest.mark.parametrize("dimension,primitive,integrator", conservation.CASES)
def test_cgl_lf_smr_conservation(tmp_path, dimension, primitive, integrator):
    states = []
    for ranks in (1, 4):
        directory = tmp_path / f"ranks{ranks}"
        directory.mkdir()
        testutils.mpi_run(conservation.INPUT,
                         conservation.flags(directory, dimension, primitive, integrator),
                         threads=ranks)
        states.append(conservation.check(directory, dimension))
    conservation.assert_same(*states)


@pytest.mark.parametrize("primitive", [False, True])
def test_cgl_lf_smr_sts_converges_to_explicit(tmp_path, primitive):
    """Qualify the composed refined operator at a fixed mesh and common endpoint."""
    pressures, steps, work = [], [], []
    reference_state = None
    # Resolve first-order reference error below the smaller second-order STS errors.
    for integrator, budget in (("explicit", 0.005), ("explicit", 0.0025),
                               ("sts", 8.0), ("sts", 4.0), ("sts", 2.0)):
        directory = tmp_path / f"{integrator}{budget}"
        directory.mkdir()
        flags = conservation.flags(directory, 2, primitive, integrator)
        flags.append(f"time/cfl_number={budget if integrator == 'explicit' else 0.4}")
        if integrator == "sts":
            flags.append(f"time/sts_max_dt_ratio={budget}")
        testutils.mpi_run(conservation.INPUT, flags, threads=1)
        state = conservation.check(directory, 2)
        if reference_state is None:
            reference_state = state
        np.testing.assert_array_equal(state["logical"], reference_state["logical"])
        np.testing.assert_array_equal(state["u"][:, [0, 1, 2, 3, 6]],
                                      reference_state["u"][:, [0, 1, 2, 3, 6]])
        np.testing.assert_array_equal(state["bcc"], reference_state["bcc"])
        u, bcc = state["u"], state["bcc"]
        rho = u[:, 0]
        bsqr = np.sum(bcc*bcc, axis=1)
        internal = u[:, 4] - 0.5*np.sum(u[:, 1:4]**2, axis=1)/rho - 0.5*bsqr
        ratio = np.exp(u[:, 5]/rho - 2*np.log(rho) + 1.5*np.log(bsqr))
        ppar = internal/(0.5 + ratio)
        pressure = np.stack((ppar, ratio*ppar), axis=1)
        assert np.all(np.isfinite(pressure)) and np.min(pressure) > 0
        pressures.append(pressure)
        history = testutils.athena_read.hst(
            str(directory / f"{conservation.BASENAME}.mhd.hst"))
        steps.append(float(np.max(np.diff(history["time"]))))
        rhs = (history["lf_nstage"][-1] - history["lf_nstage"][0])/rho.size
        work.append((state["cycle"], rhs/state["cycle"]))

    # Weight each cell by its refinement-level volume, not by the number of blocks.
    weights = 2.0**(-2*reference_state["logical"][:, 3])
    reference_error = np.average(np.mean(abs(pressures[0] - pressures[1]),
                                         axis=(2, 3, 4)), axis=0, weights=weights)
    errors = np.array([np.average(np.mean(abs(state - pressures[1]), axis=(2, 3, 4)),
                                 axis=0, weights=weights) for state in pressures[2:]])
    scale = np.average(np.mean(pressures[1], axis=(2, 3, 4)), axis=0, weights=weights)
    evidence = (f"primitive={primitive}; ppar/pperp L1 errors={errors}; "
                f"explicit refinement={reference_error}; max dt={steps}; "
                f"(cycles,RHS/cycle)={work}")
    print(evidence)
    assert 0.45 < steps[1]/steps[0] < 0.55, evidence
    assert np.all(np.asarray(steps[3:]) < 0.6*np.asarray(steps[2:-1])), evidence
    assert all(item[1] == 2 for item in work[:2]), evidence
    assert all(item[1] > 2 for item in work[2:]), evidence
    assert np.all(errors > 1e-13), evidence
    assert np.all(errors[1:] < 0.85*errors[:-1]), evidence
    assert np.all(reference_error < 0.5*errors[-1]), evidence
    assert np.all(errors[-1]/scale < 1e-4), evidence
