#!/usr/bin/env python3
"""Compare relativistic Boris, Vay, and Higuera-Cary pushers.

This is a standalone Phase 0 experiment.  It intentionally does not import or
exercise AthenaK production code.  The normalized variables are c = q/m = 1,
with momentum-like velocity u = gamma*v and gamma = sqrt(1 + |u|^2).

The Vay implementation follows Eqs. (9)-(12) of:
  J.-L. Vay, Physics of Plasmas 15, 056701 (2008)
  https://doi.org/10.1063/1.2837054

The Higuera-Cary implementation follows Eqs. (20), (21), and (24) of:
  A. V. Higuera and J. R. Cary, Physics of Plasmas 24, 052104 (2017)
  https://doi.org/10.1063/1.4979989
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import subprocess
import sys
import time
from datetime import datetime
from datetime import timezone
from pathlib import Path
from typing import Callable
from typing import Dict
from typing import List
from typing import Optional
from typing import Sequence
from typing import Tuple

import numpy as np

try:
    import scipy
    from scipy.integrate import solve_ivp
except ImportError:
    scipy = None
    solve_ivp = None


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    ROOT / "docs" / "source" / "modules" / "figures"
    / "cr_tracer_relativistic_acceleration"
    / "cr_relativistic_pusher_spike.json"
)
DEFAULT_CRITERIA_MANIFEST = DEFAULT_OUTPUT.with_name(
    "cr_relativistic_pusher_preregistered_criteria.json"
)
DEFAULT_TRANSCRIPT = DEFAULT_OUTPUT.with_name(
    "cr_relativistic_pusher_spike_transcript.txt"
)
HC_PAPER_DOI = "10.1063/1.4979989"
VAY_PAPER_DOI = "10.1063/1.2837054"
METHODS = ("relativistic_boris", "vay", "higuera_cary")
VOLUME_PRESERVING_METHODS = ("relativistic_boris", "higuera_cary")
EQUILIBRIUM_PRESERVING_METHODS = ("vay", "higuera_cary")
PHASE_0_WORKING_SELECTION = "higuera_cary"
EXPECTED_CRITERIA_MANIFEST_SHA256 = (
    "a8659125ecd37f21eb9eaff4dab869b806e31352b72bf3bcaea087c9b3d0af0c"
)
EXPECTED_CRITERIA_CONTRACT_SHA256 = (
    "bf80dbb678aad52db43dee0b74d04ddcbdc764364fececff65ab6a1eb3e3bb53"
)

Vector = np.ndarray
Pusher = Callable[[Vector, Vector, Vector, float], Vector]


def _vec(values: Sequence[float]) -> Vector:
    return np.asarray(values, dtype=float)


def _norm(values: Vector) -> float:
    return float(np.linalg.norm(values))


def _gamma(momentum: Vector) -> float:
    return math.sqrt(1.0 + float(np.dot(momentum, momentum)))


def _velocity(momentum: Vector) -> Vector:
    return momentum / _gamma(momentum)


def _finite(values: Vector) -> bool:
    return bool(np.all(np.isfinite(values)))


def _hc_rotation_gamma(momentum_minus: Vector,
                       beta: Vector) -> Tuple[float, Dict[str, object]]:
    """Return the Higuera-Cary Eq. (20) rotation gamma.

    The paper writes the positive quadratic root directly.  When sigma is
    negative, the conjugate form below avoids subtracting nearly equal large
    values.  This algebraically equivalent form is cancellation-resistant for
    the bounded spike inputs.
    """
    gamma_minus = _gamma(momentum_minus)
    beta2 = float(np.dot(beta, beta))
    beta_dot_u = float(np.dot(beta, momentum_minus))
    sigma = gamma_minus*gamma_minus - beta2
    positive_term = beta2 + beta_dot_u*beta_dot_u
    radical = math.hypot(sigma, 2.0*math.sqrt(positive_term))
    if sigma >= 0.0:
        gamma2 = 0.5*(sigma + radical)
        evaluation = "direct_positive_root"
    else:
        gamma2 = 2.0*positive_term/(radical - sigma)
        evaluation = "conjugate_positive_root"
    return math.sqrt(gamma2), {
        "gamma_minus": gamma_minus,
        "beta2": beta2,
        "beta_dot_u": beta_dot_u,
        "sigma": sigma,
        "positive_term": positive_term,
        "radical": radical,
        "gamma_rotation": math.sqrt(gamma2),
        "evaluation": evaluation,
    }


def _positive_biquadratic_root(sigma: float, positive_term: float) -> float:
    """Evaluate the positive root cancellation-resistantly for bounded inputs."""
    radical = math.hypot(sigma, 2.0*math.sqrt(positive_term))
    if sigma >= 0.0:
        root2 = 0.5*(sigma + radical)
    else:
        root2 = 2.0*positive_term/(radical - sigma)
    return math.sqrt(root2)


def _rotate(momentum_minus: Vector, beta: Vector,
            gamma_rotation: float) -> Vector:
    """Apply the familiar Boris rotation used by both compared methods."""
    t = beta/gamma_rotation
    s = 2.0*t/(1.0 + float(np.dot(t, t)))
    momentum_prime = momentum_minus + np.cross(momentum_minus, t)
    return momentum_minus + np.cross(momentum_prime, s)


def _relativistic_boris(momentum: Vector, electric: Vector,
                        magnetic: Vector, dt: float) -> Vector:
    epsilon = 0.5*dt*electric
    beta = 0.5*dt*magnetic
    momentum_minus = momentum + epsilon
    momentum_plus = _rotate(momentum_minus, beta, _gamma(momentum_minus))
    return momentum_plus + epsilon


def _higuera_cary(momentum: Vector, electric: Vector,
                  magnetic: Vector, dt: float) -> Vector:
    epsilon = 0.5*dt*electric
    beta = 0.5*dt*magnetic
    momentum_minus = momentum + epsilon
    gamma_rotation, _ = _hc_rotation_gamma(momentum_minus, beta)
    momentum_plus = _rotate(momentum_minus, beta, gamma_rotation)
    return momentum_plus + epsilon


def _vay(momentum: Vector, electric: Vector,
         magnetic: Vector, dt: float) -> Vector:
    """Apply the explicit Vay 2008 Eqs. (9)-(12) momentum map."""
    tau = 0.5*dt*magnetic
    momentum_sharp = momentum + dt*(
        electric + 0.5*np.cross(_velocity(momentum), magnetic))
    gamma_sharp = _gamma(momentum_sharp)
    tau2 = float(np.dot(tau, tau))
    momentum_star = float(np.dot(momentum_sharp, tau))
    sigma = gamma_sharp*gamma_sharp - tau2
    gamma_final = _positive_biquadratic_root(
        sigma, tau2 + momentum_star*momentum_star)
    t = tau/gamma_final
    scale = 1.0/(1.0 + float(np.dot(t, t)))
    return scale*(
        momentum_sharp
        + float(np.dot(momentum_sharp, t))*t
        + np.cross(momentum_sharp, t)
    )


PUSHERS: Dict[str, Pusher] = {
    "relativistic_boris": _relativistic_boris,
    "vay": _vay,
    "higuera_cary": _higuera_cary,
}


def _run_pusher(method: str, momentum0: Vector, electric: Vector,
                magnetic: Vector, dt: float, steps: int,
                target_momentum: Optional[Vector] = None) -> Dict[str, object]:
    """Advance momentum and trapezoid-integrate position for diagnostics."""
    pusher = PUSHERS[method]
    momentum = momentum0.copy()
    position = np.zeros(3)
    gamma0 = _gamma(momentum)
    field_work = 0.0
    max_abs_gamma_delta = 0.0
    max_abs_momentum_delta = 0.0
    finite = True
    started = time.perf_counter()
    for _ in range(steps):
        velocity_before = _velocity(momentum)
        momentum = pusher(momentum, electric, magnetic, dt)
        velocity_after = _velocity(momentum)
        centered_velocity = 0.5*(velocity_before + velocity_after)
        position += dt*centered_velocity
        field_work += dt*float(np.dot(electric, centered_velocity))
        finite = finite and _finite(momentum) and _finite(position)
        max_abs_gamma_delta = max(
            max_abs_gamma_delta, abs(_gamma(momentum) - gamma0))
        if target_momentum is not None:
            max_abs_momentum_delta = max(
                max_abs_momentum_delta, _norm(momentum - target_momentum))
    gamma_final = _gamma(momentum)
    kinetic_energy_change = gamma_final - gamma0
    work_closure_error = kinetic_energy_change - field_work
    return {
        "momentum": momentum,
        "position": position,
        "gamma_initial": gamma0,
        "gamma_final": gamma_final,
        "kinetic_energy_change": kinetic_energy_change,
        "electric_field_work_trapezoid": field_work,
        "work_closure_error": work_closure_error,
        "abs_work_closure_error": abs(work_closure_error),
        "max_abs_gamma_delta": max_abs_gamma_delta,
        "max_rel_gamma_delta": max_abs_gamma_delta/gamma0,
        "max_abs_momentum_delta": max_abs_momentum_delta,
        "finite": finite,
        "wall_seconds": time.perf_counter() - started,
    }


def _newtonian_boris(momentum: Vector, electric: Vector,
                     magnetic: Vector, dt: float) -> Vector:
    """Apply the non-relativistic Boris map for a low-speed limit check."""
    epsilon = 0.5*dt*electric
    beta = 0.5*dt*magnetic
    momentum_minus = momentum + epsilon
    momentum_plus = _rotate(momentum_minus, beta, 1.0)
    return momentum_plus + epsilon


def _run_map(pusher: Pusher, momentum0: Vector, electric: Vector,
             magnetic: Vector, dt: float, steps: int) -> Vector:
    momentum = momentum0.copy()
    for _ in range(steps):
        momentum = pusher(momentum, electric, magnetic, dt)
    return momentum


def _numerical_jacobian(pusher: Pusher, momentum: Vector, electric: Vector,
                        magnetic: Vector, dt: float,
                        perturbation: float) -> Vector:
    jacobian = np.zeros((3, 3))
    for column in range(3):
        delta = np.zeros(3)
        delta[column] = perturbation
        jacobian[:, column] = (
            pusher(momentum + delta, electric, magnetic, dt)
            - pusher(momentum - delta, electric, magnetic, dt)
        )/(2.0*perturbation)
    return jacobian


def _lorentz_rhs(_time_value: float, state: Vector, electric: Vector,
                 magnetic: Vector) -> Vector:
    momentum = state[3:6]
    velocity = _velocity(momentum)
    return np.concatenate((velocity, electric + np.cross(velocity, magnetic)))


def _rk4_reference(momentum0: Vector, electric: Vector, magnetic: Vector,
                   duration: float, max_step: float) -> Dict[str, object]:
    """Independent fixed-step fallback for environments without SciPy."""
    steps = max(1, int(math.ceil(abs(duration)/max_step)))
    dt = duration/steps
    state = np.concatenate((np.zeros(3), momentum0.copy()))
    time_value = 0.0
    for _ in range(steps):
        k1 = _lorentz_rhs(time_value, state, electric, magnetic)
        k2 = _lorentz_rhs(
            time_value + 0.5*dt, state + 0.5*dt*k1, electric, magnetic)
        k3 = _lorentz_rhs(
            time_value + 0.5*dt, state + 0.5*dt*k2, electric, magnetic)
        k4 = _lorentz_rhs(time_value + dt, state + dt*k3, electric, magnetic)
        state += (dt/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
        time_value += dt
    return {
        "engine": "independent_fixed_step_rk4_fallback",
        "success": True,
        "steps": steps,
        "position": state[0:3],
        "momentum": state[3:6],
    }


def _reference(momentum0: Vector, electric: Vector, magnetic: Vector,
               duration: float, max_step: float) -> Dict[str, object]:
    """Integrate the prescribed-field Lorentz ODE independently."""
    if solve_ivp is None:
        return _rk4_reference(
            momentum0, electric, magnetic, duration, max_step/8.0)

    initial = np.concatenate((np.zeros(3), momentum0.copy()))
    solution = solve_ivp(
        _lorentz_rhs, (0.0, duration), initial,
        args=(electric, magnetic), method="DOP853",
        rtol=1.0e-12, atol=1.0e-14, max_step=max_step,
    )
    if not solution.success:
        raise RuntimeError(f"solve_ivp reference failed: {solution.message}")
    final = solution.y[:, -1]
    return {
        "engine": "scipy.integrate.solve_ivp_DOP853",
        "success": bool(solution.success),
        "function_evaluations": int(solution.nfev),
        "position": final[0:3],
        "momentum": final[3:6],
    }


def _momentum_error(result: Dict[str, object],
                    reference: Dict[str, object]) -> float:
    return _norm(result["momentum"] - reference["momentum"])


def _wrapped_angle_error(actual: float, expected: float) -> float:
    return abs(math.atan2(math.sin(actual - expected),
                          math.cos(actual - expected)))


def _fit_log_slope(step_sizes: Sequence[float],
                   errors: Sequence[float]) -> float:
    xs = np.log(np.asarray(step_sizes, dtype=float))
    ys = np.log(np.asarray(errors, dtype=float))
    return float(np.polyfit(xs, ys, deg=1)[0])


def _equilibrium_momentum(electric: Vector, magnetic: Vector) -> Vector:
    magnetic2 = float(np.dot(magnetic, magnetic))
    drift = np.cross(electric, magnetic)/magnetic2
    speed2 = float(np.dot(drift, drift))
    if speed2 >= 1.0:
        raise ValueError("Crossed-field equilibrium requires |E x B|/|B|^2 < 1")
    lorentz_denominator = math.sqrt(
        (1.0 - math.sqrt(speed2)) * (1.0 + math.sqrt(speed2)))
    return drift / lorentz_denominator


def _pure_e_subcase(momentum0: Vector, electric: Vector,
                    dt: float, steps: int) -> Dict[str, object]:
    duration = dt*steps
    expected_momentum = momentum0 + duration*electric
    expected_gamma = _gamma(expected_momentum)
    expected_kinetic_energy_change = expected_gamma - _gamma(momentum0)
    methods = {}
    for method in METHODS:
        result = _run_pusher(
            method, momentum0, electric, np.zeros(3), dt, steps)
        methods[method] = {
            "finite": result["finite"],
            "endpoint_momentum": result["momentum"],
            "endpoint_momentum_error_vs_closed_form": _norm(
                result["momentum"] - expected_momentum),
            "endpoint_gamma_error_vs_closed_form": abs(
                result["gamma_final"] - expected_gamma),
            "kinetic_energy_change": result["kinetic_energy_change"],
            "electric_field_work_trapezoid": result[
                "electric_field_work_trapezoid"],
            "work_closure_error": result["work_closure_error"],
        }
    return {
        "inputs": {
            "momentum0": momentum0,
            "electric": electric,
            "magnetic": np.zeros(3),
            "dt": dt,
            "steps": steps,
            "duration": duration,
        },
        "closed_form": {
            "endpoint_momentum": expected_momentum,
            "endpoint_gamma": expected_gamma,
            "kinetic_energy_change": expected_kinetic_energy_change,
        },
        "methods": methods,
    }


def _pure_e_closed_form() -> Dict[str, object]:
    return {
        "purpose": (
            "Check pure-E acceleration and deceleration against the exact "
            "momentum law u(t) = u(0) + E*t.  The field-work diagnostic is "
            "accumulated from E dot v trapezoid quadrature, not inferred from "
            "the kinetic-energy change."),
        "limit": (
            "The closed form is exact for prescribed spatially uniform E and "
            "B = 0 in normalized c = q/m = 1 units.  The independent field-work "
            "quadrature is expected to converge rather than be exact."),
        "acceleration": _pure_e_subcase(
            _vec((0.1, 0.0, 0.0)), _vec((0.25, 0.0, 0.0)),
            dt=0.05, steps=240),
        "deceleration": _pure_e_subcase(
            _vec((3.0, 0.0, 0.0)), _vec((-0.1, 0.0, 0.0)),
            dt=0.05, steps=400),
    }


def _low_speed_agreement() -> Dict[str, object]:
    momentum0 = _vec((1.0e-4, -2.0e-4, 3.0e-4))
    electric = _vec((2.0e-5, -1.0e-5, 0.0))
    magnetic = _vec((0.0, 0.0, 0.7))
    dt = 0.01
    steps = 500
    baseline = _run_map(
        _newtonian_boris, momentum0, electric, magnetic, dt, steps)
    endpoint_momenta = {}
    methods = {}
    max_speed = _norm(_velocity(momentum0))
    for method in METHODS:
        result = _run_pusher(method, momentum0, electric, magnetic, dt, steps)
        endpoint_momenta[method] = result["momentum"]
        max_speed = max(max_speed, _norm(_velocity(result["momentum"])))
        methods[method] = {
            "finite": result["finite"],
            "endpoint_momentum": result["momentum"],
            "endpoint_momentum_error_vs_newtonian_boris": _norm(
                result["momentum"] - baseline),
        }
    pairwise = {}
    max_pairwise_separation = 0.0
    for left_index, left in enumerate(METHODS):
        for right in METHODS[left_index + 1:]:
            key = f"{left}_vs_{right}"
            separation = _norm(endpoint_momenta[left] - endpoint_momenta[right])
            pairwise[key] = separation
            max_pairwise_separation = max(max_pairwise_separation, separation)
    return {
        "purpose": (
            "Check that the three relativistic maps collapse onto each other "
            "and the Newtonian Boris limit when |v| is much less than c."),
        "limit": (
            "This is a low-speed numerical agreement probe, not a derivation "
            "of the non-relativistic limit or a production regression."),
        "inputs": {
            "momentum0": momentum0,
            "electric": electric,
            "magnetic": magnetic,
            "dt": dt,
            "steps": steps,
            "duration": dt*steps,
        },
        "newtonian_boris_endpoint_momentum": baseline,
        "maximum_sampled_endpoint_speed": max_speed,
        "methods": methods,
        "pairwise_endpoint_momentum_separation": pairwise,
        "max_pairwise_endpoint_momentum_separation": max_pairwise_separation,
    }


def _work_closure() -> Dict[str, object]:
    momentum0 = _vec((1.1, -0.3, 0.2))
    electric = _vec((0.2, -0.1, 0.05))
    magnetic = _vec((0.1, 0.2, 1.0))
    duration = 12.0
    step_sizes = (0.2, 0.1, 0.05, 0.025, 0.0125)
    methods = {}
    for method in METHODS:
        errors: List[float] = []
        rows = []
        for dt in step_sizes:
            steps = int(round(duration/dt))
            result = _run_pusher(
                method, momentum0, electric, magnetic, dt, steps)
            error = result["abs_work_closure_error"]
            errors.append(error)
            rows.append({
                "dt": dt,
                "steps": steps,
                "kinetic_energy_change": result["kinetic_energy_change"],
                "electric_field_work_trapezoid": result[
                    "electric_field_work_trapezoid"],
                "work_closure_error": result["work_closure_error"],
                "abs_work_closure_error": error,
            })
        methods[method] = {
            "rows": rows,
            "fitted_abs_work_closure_slope_finest_four": _fit_log_slope(
                step_sizes[-4:], errors[-4:]),
        }
    return {
        "purpose": (
            "Accumulate prescribed electric-field work as the independent "
            "trapezoid sum of E dot v and compare it with Delta(gamma - 1)."),
        "limit": (
            "This is a convergence diagnostic for field-work quadrature in "
            "constant prescribed fields.  It is not a self-consistent "
            "particle-field energy-conservation test."),
        "inputs": {
            "momentum0": momentum0,
            "electric": electric,
            "magnetic": magnetic,
            "duration": duration,
            "step_sizes": step_sizes,
        },
        "methods": methods,
    }


def _higuera_cary_direct_work_identity() -> Dict[str, object]:
    momentum0 = _vec((1.1, -0.3, 0.2))
    electric = _vec((0.2, -0.1, 0.05))
    magnetic = _vec((0.1, 0.2, 1.0))
    duration = 12.0
    step_sizes = (0.2, 0.1, 0.05, 0.025, 0.0125)
    rows = []
    all_finite = True
    max_abs_single_step_identity_residual = 0.0
    max_abs_accumulated_identity_residual = 0.0
    for dt in step_sizes:
        steps = int(round(duration/dt))
        momentum = momentum0.copy()
        gamma_initial = _gamma(momentum)
        direct_work = 0.0
        max_abs_step_residual = 0.0
        finite = True
        for _ in range(steps):
            momentum_before = momentum.copy()
            gamma_before = _gamma(momentum_before)
            momentum = _higuera_cary(momentum, electric, magnetic, dt)
            gamma_after = _gamma(momentum)
            delta_work = dt*float(
                np.dot(electric, momentum_before + momentum)
            )/(gamma_before + gamma_after)
            residual = (gamma_after - gamma_before) - delta_work
            direct_work += delta_work
            max_abs_step_residual = max(
                max_abs_step_residual, abs(residual))
            finite = finite and _finite(momentum)
        kinetic_energy_change = _gamma(momentum) - gamma_initial
        accumulated_identity_residual = kinetic_energy_change - direct_work
        all_finite = all_finite and finite
        max_abs_single_step_identity_residual = max(
            max_abs_single_step_identity_residual, max_abs_step_residual)
        max_abs_accumulated_identity_residual = max(
            max_abs_accumulated_identity_residual,
            abs(accumulated_identity_residual))
        rows.append({
            "dt": dt,
            "steps": steps,
            "finite": finite,
            "kinetic_energy_change": kinetic_energy_change,
            "discrete_direct_electric_work": direct_work,
            "max_abs_single_step_identity_residual": max_abs_step_residual,
            "accumulated_identity_residual": accumulated_identity_residual,
            "abs_accumulated_identity_residual": abs(
                accumulated_identity_residual),
        })
    return {
        "purpose": (
            "Check the exact discrete Higuera-Cary direct-work identity "
            "separately from the trapezoid-work convergence diagnostic."),
        "identity": (
            "DeltaW = alpha*h*cE dot (w_n+w_np1)/(gamma_n+gamma_np1)"
        ),
        "normalized_variables": {
            "alpha": 1.0,
            "h": "dt",
            "speed_of_light_c": 1.0,
            "w": "u = gamma*v",
            "DeltaW": "gamma_np1 - gamma_n",
        },
        "limit": (
            "This is a bounded prescribed-field floating-point check of the "
            "discrete identity, not a self-consistent particle-field "
            "energy-conservation test."),
        "inputs": {
            "momentum0": momentum0,
            "electric": electric,
            "magnetic": magnetic,
            "duration": duration,
            "step_sizes": step_sizes,
        },
        "rows": rows,
        "all_finite": all_finite,
        "max_abs_single_step_identity_residual": (
            max_abs_single_step_identity_residual),
        "max_abs_accumulated_identity_residual": (
            max_abs_accumulated_identity_residual),
    }


def _phase_volume_jacobian() -> Dict[str, object]:
    probes = (
        {
            "momentum": _vec((1.7, -0.8, 0.45)),
            "electric": _vec((0.4, -0.25, 0.15)),
            "magnetic": _vec((0.3, 0.5, 1.2)),
            "dt": 0.7,
        },
        {
            "momentum": _vec((4.0, 1.5, -0.5)),
            "electric": _vec((-0.2, 0.35, 0.1)),
            "magnetic": _vec((0.1, -0.4, 1.8)),
            "dt": 1.1,
        },
    )
    perturbations = (1.0e-5, 3.0e-6, 1.0e-6)
    methods = {}
    for method in METHODS:
        rows = []
        max_abs_det_minus_one = 0.0
        for probe_index, probe in enumerate(probes):
            for perturbation in perturbations:
                jacobian = _numerical_jacobian(
                    PUSHERS[method], perturbation=perturbation, **probe)
                determinant = float(np.linalg.det(jacobian))
                abs_det_minus_one = abs(determinant - 1.0)
                max_abs_det_minus_one = max(
                    max_abs_det_minus_one, abs_det_minus_one)
                rows.append({
                    "probe_index": probe_index,
                    "perturbation": perturbation,
                    "determinant": determinant,
                    "abs_determinant_minus_one": abs_det_minus_one,
                })
        methods[method] = {
            "rows": rows,
            "max_abs_determinant_minus_one": max_abs_det_minus_one,
        }
    return {
        "purpose": (
            "Numerically probe the local Jacobian determinant of each "
            "constant-field momentum map.  Boris and Higuera-Cary should "
            "preserve differential volume; Vay is included as a discriminator."),
        "limit": (
            "This central-difference probe covers the isolated momentum submap "
            "at two sample points.  It is not a proof for a full staggered "
            "position-momentum integrator or for AthenaK integration."),
        "inputs": {
            "probes": probes,
            "central_difference_perturbations": perturbations,
        },
        "methods": methods,
    }


def _uniform_b_energy_phase() -> Dict[str, object]:
    electric = _vec((0.0, 0.0, 0.0))
    magnetic = _vec((0.0, 0.0, 1.0))
    momentum0 = _vec((3.0, 0.0, 0.4))
    dt = 0.15
    steps = 800
    duration = dt*steps
    gamma0 = _gamma(momentum0)
    target_angle = -duration/gamma0
    reference = _reference(
        momentum0, electric, magnetic, duration, max_step=duration/2000.0)
    methods = {}
    for method in METHODS:
        result = _run_pusher(method, momentum0, electric, magnetic, dt, steps)
        actual_angle = math.atan2(result["momentum"][1],
                                  result["momentum"][0])
        methods[method] = {
            "finite": result["finite"],
            "max_rel_gamma_delta": result["max_rel_gamma_delta"],
            "endpoint_phase_error_radians": _wrapped_angle_error(
                actual_angle, target_angle),
            "endpoint_momentum_error_vs_reference": _momentum_error(
                result, reference),
            "endpoint_relative_momentum_error_vs_reference": (
                _momentum_error(result, reference)/_norm(
                    reference["momentum"])),
        }
    return {
        "purpose": "Check no-E energy conservation and uniform-B gyro phase.",
        "inputs": {
            "momentum0": momentum0,
            "electric": electric,
            "magnetic": magnetic,
            "dt": dt,
            "steps": steps,
            "duration": duration,
        },
        "reference": reference,
        "methods": methods,
    }


def _crossed_eb_drift() -> Dict[str, object]:
    electric = _vec((0.6, 0.0, 0.0))
    magnetic = _vec((0.0, 0.0, 1.0))
    momentum0 = _vec((0.25, 0.1, 0.0))
    dt = 0.2
    steps = 4000
    duration = dt*steps
    expected_drift = np.cross(electric, magnetic)/np.dot(magnetic, magnetic)
    reference = _reference(
        momentum0, electric, magnetic, duration, max_step=duration/10000.0)
    reference_average_velocity = reference["position"]/duration
    methods = {}
    for method in METHODS:
        result = _run_pusher(method, momentum0, electric, magnetic, dt, steps)
        average_velocity = result["position"]/duration
        methods[method] = {
            "finite": result["finite"],
            "average_velocity_trapezoid_diagnostic": average_velocity,
            "average_velocity_error_vs_expected_drift": _norm(
                average_velocity - expected_drift),
            "average_velocity_error_vs_reference": _norm(
                average_velocity - reference_average_velocity),
            "endpoint_momentum_error_vs_reference": _momentum_error(
                result, reference),
        }
    return {
        "purpose": (
            "Measure long-window E x B drift with a bounded gyro transient; "
            "the pusher position is diagnostic trapezoid quadrature only."),
        "inputs": {
            "momentum0": momentum0,
            "electric": electric,
            "magnetic": magnetic,
            "expected_drift": expected_drift,
            "dt": dt,
            "steps": steps,
            "duration": duration,
        },
        "reference": {
            **reference,
            "average_velocity": reference_average_velocity,
            "average_velocity_error_vs_expected_drift": _norm(
                reference_average_velocity - expected_drift),
        },
        "methods": methods,
    }


def _crossed_field_equilibrium() -> Dict[str, object]:
    electric = _vec((0.9, 0.0, 0.0))
    magnetic = _vec((0.0, 0.0, 1.0))
    momentum0 = _equilibrium_momentum(electric, magnetic)
    dt = 0.8
    steps = 200
    duration = dt*steps
    reference = _reference(
        momentum0, electric, magnetic, duration, max_step=duration/1000.0)
    methods = {}
    for method in METHODS:
        result = _run_pusher(
            method, momentum0, electric, magnetic, dt, steps,
            target_momentum=momentum0)
        methods[method] = {
            "finite": result["finite"],
            "max_abs_momentum_departure": result["max_abs_momentum_delta"],
            "endpoint_momentum_departure": _norm(
                result["momentum"] - momentum0),
            "endpoint_momentum_error_vs_reference": _momentum_error(
                result, reference),
        }
    return {
        "purpose": (
            "Expose exact crossed-field equilibrium cancellation as an "
            "equivalent boosted-frame discriminator.  For perpendicular "
            "|E| < |B| fields, the E x B drift frame removes E; a correct "
            "stationary drift should therefore remain force-free."),
        "limit": (
            "This checks the equivalent constant-field cancellation condition "
            "without implementing a general Lorentz transformation of fields "
            "and trajectories."),
        "inputs": {
            "momentum0": momentum0,
            "electric": electric,
            "magnetic": magnetic,
            "dt": dt,
            "steps": steps,
            "duration": duration,
        },
        "reference": reference,
        "methods": methods,
    }


def _high_gamma() -> Dict[str, object]:
    electric = _vec((0.0, 0.0, 0.0))
    magnetic = _vec((0.0, 0.0, 1.0))
    momentum0 = _vec((1000.0, 0.0, 0.0))
    gamma0 = _gamma(momentum0)
    orbits = 4
    steps_per_orbit = 32
    steps = orbits*steps_per_orbit
    duration = orbits*2.0*math.pi*gamma0
    dt = duration/steps
    target_angle = -duration/gamma0
    reference = _reference(
        momentum0, electric, magnetic, duration, max_step=duration/4000.0)
    methods = {}
    for method in METHODS:
        result = _run_pusher(method, momentum0, electric, magnetic, dt, steps)
        actual_angle = math.atan2(result["momentum"][1],
                                  result["momentum"][0])
        methods[method] = {
            "finite": result["finite"],
            "max_rel_gamma_delta": result["max_rel_gamma_delta"],
            "endpoint_phase_error_radians": _wrapped_angle_error(
                actual_angle, target_angle),
            "endpoint_momentum_error_vs_reference": _momentum_error(
                result, reference),
            "endpoint_relative_momentum_error_vs_reference": (
                _momentum_error(result, reference)/_norm(
                    reference["momentum"])),
        }
    return {
        "purpose": (
            "Compare relativistic gyro phase at gamma approximately 1000; "
            "uniform-B energy remains an invariant for both methods."),
        "inputs": {
            "momentum0": momentum0,
            "gamma0": gamma0,
            "electric": electric,
            "magnetic": magnetic,
            "orbits": orbits,
            "steps_per_orbit": steps_per_orbit,
            "dt": dt,
            "steps": steps,
            "duration": duration,
        },
        "reference": reference,
        "methods": methods,
    }


def _timestep_convergence() -> Dict[str, object]:
    electric = _vec((0.35, 0.1, 0.0))
    magnetic = _vec((0.0, 0.0, 1.0))
    momentum0 = _vec((1.5, -0.4, 0.7))
    duration = 8.0
    step_sizes = (0.4, 0.2, 0.1, 0.05, 0.025)
    reference = _reference(
        momentum0, electric, magnetic, duration, max_step=duration/4000.0)
    methods = {}
    for method in METHODS:
        errors: List[float] = []
        rows = []
        for dt in step_sizes:
            steps = int(round(duration/dt))
            result = _run_pusher(
                method, momentum0, electric, magnetic, dt, steps)
            error = _momentum_error(result, reference)
            errors.append(error)
            rows.append({
                "dt": dt,
                "steps": steps,
                "endpoint_momentum_error_vs_reference": error,
            })
        methods[method] = {
            "rows": rows,
            "fitted_endpoint_error_slope_all_points": _fit_log_slope(
                step_sizes, errors),
            "fitted_endpoint_error_slope_finest_four": _fit_log_slope(
                step_sizes[-4:], errors[-4:]),
        }
    return {
        "purpose": "Check expected second-order endpoint convergence.",
        "inputs": {
            "momentum0": momentum0,
            "electric": electric,
            "magnetic": magnetic,
            "duration": duration,
            "step_sizes": step_sizes,
        },
        "reference": reference,
        "methods": methods,
    }


def _reversibility() -> Dict[str, object]:
    electric = _vec((0.3, -0.2, 0.1))
    magnetic = _vec((0.15, 0.05, 1.0))
    momentum0 = _vec((2.0, -0.7, 0.4))
    dt = 0.37
    steps = 500
    methods = {}
    for method in METHODS:
        forward = _run_pusher(
            method, momentum0, electric, magnetic, dt, steps)
        reverse = _run_pusher(
            method, forward["momentum"], electric, magnetic, -dt, steps)
        methods[method] = {
            "finite": forward["finite"] and reverse["finite"],
            "roundtrip_abs_momentum_error": _norm(
                reverse["momentum"] - momentum0),
        }
    return {
        "purpose": "Check forward-then-backward momentum-map reversibility.",
        "inputs": {
            "momentum0": momentum0,
            "electric": electric,
            "magnetic": magnetic,
            "dt": dt,
            "steps_each_direction": steps,
        },
        "methods": methods,
    }


def _finite_state_near_limit() -> Dict[str, object]:
    ratio = 1.0 - 1.0e-12
    electric = _vec((ratio, 0.0, 0.0))
    magnetic = _vec((0.0, 0.0, 1.0))
    momentum0 = _equilibrium_momentum(electric, magnetic)
    dt = 0.5
    steps = 1000
    methods = {}
    for method in METHODS:
        result = _run_pusher(
            method, momentum0, electric, magnetic, dt, steps,
            target_momentum=momentum0)
        methods[method] = {
            "finite": result["finite"],
            "max_abs_momentum_departure": result["max_abs_momentum_delta"],
            "max_rel_momentum_departure": (
                result["max_abs_momentum_delta"]/_norm(momentum0)),
        }

    strong_momentum = _vec((1000.0, 2.0, -3.0))
    strong_beta = _vec((0.0, 0.0, 5.0e7))
    strong_gamma, strong_details = _hc_rotation_gamma(
        strong_momentum, strong_beta)
    return {
        "purpose": (
            "Check finite states close to the |E|/|B| = 1 equilibrium limit "
            "and exercise the cancellation-resistant HC positive-root "
            "evaluation when sigma < 0 within the bounded spike."),
        "near_exb_limit": {
            "inputs": {
                "electric_to_magnetic_ratio": ratio,
                "equilibrium_gamma": _gamma(momentum0),
                "momentum0": momentum0,
                "electric": electric,
                "magnetic": magnetic,
                "dt": dt,
                "steps": steps,
            },
            "methods": methods,
        },
        "strong_rotation_cancellation_guard": {
            "momentum_minus": strong_momentum,
            "beta": strong_beta,
            "gamma_rotation": strong_gamma,
            "finite": math.isfinite(strong_gamma),
            "details": strong_details,
        },
    }


def _run_sheared_field_orbit(method: str, momentum0: Vector, dt: float,
                             steps: int, a: float, b: float) -> Dict[str, object]:
    """Run a symmetric drift-kick-drift orbit for the HC 2017 test fields."""
    pusher = PUSHERS[method]
    momentum = momentum0.copy()
    position = np.zeros(3)
    gamma0 = _gamma(momentum)
    h0 = gamma0 + 0.5*a*position[0]*position[0]
    az0 = 0.5*b*position[1]*position[1]
    canonical_pz0 = momentum[2] + az0
    iy0 = momentum[1]*momentum[1] + (canonical_pz0 - az0)**2
    max_abs_h_error = 0.0
    max_abs_iy_error = 0.0
    max_abs_canonical_pz_error = 0.0
    field_work = 0.0
    finite = True
    for _ in range(steps):
        velocity_before = _velocity(momentum)
        position_mid = position + 0.5*dt*velocity_before
        electric = _vec((-a*position_mid[0], 0.0, 0.0))
        magnetic = _vec((b*position_mid[1], 0.0, 0.0))
        momentum = pusher(momentum, electric, magnetic, dt)
        velocity_after = _velocity(momentum)
        position = position_mid + 0.5*dt*velocity_after
        field_work += dt*float(
            np.dot(electric, 0.5*(velocity_before + velocity_after)))
        gamma = _gamma(momentum)
        az = 0.5*b*position[1]*position[1]
        h_value = gamma + 0.5*a*position[0]*position[0]
        canonical_pz = momentum[2] + az
        iy = momentum[1]*momentum[1] + (canonical_pz0 - az)**2
        max_abs_h_error = max(max_abs_h_error, abs(h_value - h0))
        max_abs_iy_error = max(max_abs_iy_error, abs(iy - iy0))
        max_abs_canonical_pz_error = max(
            max_abs_canonical_pz_error, abs(canonical_pz - canonical_pz0))
        finite = finite and _finite(momentum) and _finite(position)
        if not finite:
            break
    kinetic_energy_change = _gamma(momentum) - gamma0
    return {
        "finite": finite,
        "endpoint_position": position,
        "endpoint_momentum": momentum,
        "max_abs_h_error": max_abs_h_error,
        "max_abs_iy_error": max_abs_iy_error,
        "max_abs_canonical_pz_error": max_abs_canonical_pz_error,
        "kinetic_energy_change": kinetic_energy_change,
        "electric_field_work_trapezoid": field_work,
        "work_closure_error": kinetic_energy_change - field_work,
    }


def _resonance_large_step_stress() -> Dict[str, object]:
    a = 1.0
    b = 2.0
    h0 = 4.0
    dt_values = (0.025, 0.1)
    steps = 12000
    initial_py_values = (1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85)
    methods = {}
    for method in METHODS:
        rows = []
        summaries = {}
        for dt in dt_values:
            dt_rows = []
            for initial_py in initial_py_values:
                initial_px = math.sqrt(h0*h0 - 1.0 - initial_py*initial_py)
                momentum0 = _vec((initial_px, initial_py, 0.0))
                result = _run_sheared_field_orbit(
                    method, momentum0, dt, steps, a, b)
                row = {
                    "dt": dt,
                    "steps": steps,
                    "initial_py": initial_py,
                    "initial_momentum": momentum0,
                    **result,
                }
                rows.append(row)
                dt_rows.append(row)
            worst = max(dt_rows, key=lambda row: row["max_abs_h_error"])
            summaries[str(dt)] = {
                "max_abs_h_error": worst["max_abs_h_error"],
                "worst_initial_py": worst["initial_py"],
                "max_abs_iy_error": max(
                    row["max_abs_iy_error"] for row in dt_rows),
                "max_abs_canonical_pz_error": max(
                    row["max_abs_canonical_pz_error"] for row in dt_rows),
                "all_finite": all(row["finite"] for row in dt_rows),
            }
        methods[method] = {
            "rows": rows,
            "summaries_by_dt": summaries,
        }
    return {
        "purpose": (
            "Run a paper-inspired resonance and large-step stress scan in the "
            "Higuera-Cary 2017 sheared fields E = -a*x*xhat and B = b*y*xhat. "
            "The scan includes the reported py approximately 1.7 neighborhood "
            "and compares dt = 1/40 with dt = 1/10."),
        "limit": (
            "This is a bounded drift-kick-drift scalar-orbit probe.  It does "
            "not reproduce the paper's Poincare plots, scan every orbit, or "
            "validate a production subcycling policy."),
        "source_parameters": {
            "a": a,
            "b": b,
            "initial_h": h0,
            "dt_values": dt_values,
            "steps_per_orbit_probe": steps,
            "initial_py_values": initial_py_values,
        },
        "methods": methods,
    }


def _long_time_stress(stress_steps: int) -> Dict[str, object]:
    uniform_inputs = {
        "momentum0": _vec((100.0, 0.3, 0.1)),
        "electric": _vec((0.0, 0.0, 0.0)),
        "magnetic": _vec((0.0, 0.0, 1.0)),
        "dt": 0.4,
        "steps": stress_steps,
    }
    equilibrium_electric = _vec((0.95, 0.0, 0.0))
    equilibrium_magnetic = _vec((0.0, 0.0, 1.0))
    equilibrium_momentum = _equilibrium_momentum(
        equilibrium_electric, equilibrium_magnetic)
    equilibrium_inputs = {
        "momentum0": equilibrium_momentum,
        "electric": equilibrium_electric,
        "magnetic": equilibrium_magnetic,
        "dt": 0.2,
        "steps": stress_steps,
    }
    methods = {}
    for method in METHODS:
        uniform = _run_pusher(method, **uniform_inputs)
        equilibrium = _run_pusher(
            method, **equilibrium_inputs, target_momentum=equilibrium_momentum)
        methods[method] = {
            "uniform_b": {
                "finite": uniform["finite"],
                "max_rel_gamma_delta": uniform["max_rel_gamma_delta"],
                "wall_seconds": uniform["wall_seconds"],
            },
            "crossed_field_equilibrium": {
                "finite": equilibrium["finite"],
                "max_abs_momentum_departure": equilibrium[
                    "max_abs_momentum_delta"],
                "wall_seconds": equilibrium["wall_seconds"],
            },
        }
    return {
        "purpose": (
            "Run many prescribed-field updates to reveal non-finite states, "
            "no-E energy drift, or loss of crossed-field equilibrium."),
        "uniform_b_inputs": uniform_inputs,
        "crossed_field_equilibrium_inputs": equilibrium_inputs,
        "methods": methods,
        "timing_note": (
            "Wall times are Python-harness observations only, not production "
            "performance measurements."),
    }


def _git_metadata() -> Dict[str, str]:
    metadata = {}
    for key, command in {
        "branch": ["git", "rev-parse", "--abbrev-ref", "HEAD"],
        "commit": ["git", "rev-parse", "HEAD"],
    }.items():
        proc = subprocess.run(
            command, cwd=ROOT, capture_output=True, text=True, check=False)
        if proc.returncode == 0:
            metadata[key] = proc.stdout.strip()
    return metadata


def _operation_and_clarity_notes() -> Dict[str, object]:
    return {
        "higuera_cary_formula_check": {
            "source": (
                "A. V. Higuera and J. R. Cary, Structure-preserving "
                "second-order integration of relativistic charged particle "
                "trajectories in electromagnetic fields, Physics of Plasmas "
                "24, 052104 (2017)."),
            "doi": HC_PAPER_DOI,
            "equation_map": [
                "Eq. (20): HC gamma_rotation^2 positive root.",
                "Eq. (21): Boris uses gamma_minus = gamma(u_minus).",
                "Eq. (24): both methods use the familiar Boris rotation.",
                (
                    "Eq. (28)-(29): HC preserves the crossed-field E x B "
                    "stationary velocity; standard relativistic Boris does not."
                ),
            ],
        },
        "vay_formula_check": {
            "source": (
                "J.-L. Vay, Simulation of beams or plasmas crossing at "
                "relativistic velocity, Physics of Plasmas 15, 056701 (2008)."
            ),
            "doi": VAY_PAPER_DOI,
            "equation_map": [
                "Eq. (9): compute u_sharp from ui, E, and vi cross B/2.",
                "Eq. (11): solve explicitly for gamma_final.",
                "Eq. (12): compute u_final from u_sharp and t.",
            ],
        },
        "shared_push_shape": [
            "Compute electric half-kick epsilon = E*dt/2.",
            "Compute u_minus = u_initial + epsilon and beta = B*dt/2.",
            "Select the rotation gamma.",
            "Apply t = beta/gamma and s = 2*t/(1 + |t|^2).",
            "Apply u_prime = u_minus + u_minus cross t.",
            "Apply u_plus = u_minus + u_prime cross s.",
            "Compute u_final = u_plus + epsilon.",
        ],
        "shared_push_shape_limit": (
            "The shared Boris-rotation shape applies to relativistic Boris and "
            "Higuera-Cary.  Vay follows its separate explicit Eqs. (9)-(12)."),
        "standard_relativistic_boris": (
            "The rotation gamma is gamma(u_minus).  This is the shorter scalar "
            "path and shares all vector rotation work with Higuera-Cary."),
        "higuera_cary": (
            "The rotation gamma uses the Eq. (20) positive root.  Relative to "
            "standard relativistic Boris, this adds beta^2, beta dot u_minus, "
            "a discriminant, and a second square root while leaving the vector "
            "rotation unchanged.  The script uses the algebraically equivalent "
            "conjugate root when sigma < 0 to reduce cancellation within the "
            "bounded spike."),
        "vay": (
            "Vay is implemented directly as an experimental third comparator. "
            "It preserves the crossed-field cancellation condition but is not "
            "volume-preserving in general.  It is not the Phase 0 working "
            "selection."),
        "cost_scope": (
            "These are operation and clarity notes, not measured AthenaK cost "
            "claims.  Python wall times are emitted only as harness context."),
    }


def _residual_limits() -> List[str]:
    return [
        (
            "This Phase 0 harness compares isolated momentum maps in prescribed "
            "constant fields; it does not validate AthenaK C++ integration."
        ),
        (
            "The harness does not exercise field gather, position staggering, "
            "subcycling policy, AMR, MPI exchange, restart, or output formats."
        ),
        (
            "Normalized units c = q/m = 1 are used.  A production integration "
            "must map AthenaK variable conventions and units explicitly."
        ),
        (
            "The crossed-field drift position is trapezoid-integrated only as "
            "a diagnostic; it is not a proposed production position update."
        ),
        (
            "The cancellation-resistant HC positive-root branch is exercised "
            "within bounded spike inputs.  Production helpers must own "
            "overflow-resistant scaled norms and broader floating-point extremes."
        ),
        (
            "The Jacobian rows are finite-difference probes of isolated "
            "constant-field momentum maps, not a proof for a complete "
            "staggered phase-space integrator."
        ),
        (
            "The sheared-field resonance scan is paper-inspired and bounded. "
            "It does not reproduce the full Poincare-section analysis."
        ),
        (
            "The Vay implementation is an experimental Phase 0 comparator. "
            "Its inclusion does not authorize a production option."
        ),
        (
            "Python timings are not representative of device kernels or "
            "production throughput."
        ),
    ]


EXPECTED_CRITERION_CHECKS = {
    "P0-01": (("max_pure_e_endpoint_momentum_error", "<="),),
    "P0-02": (
        ("max_low_speed_pairwise_endpoint_momentum_separation", "<="),
    ),
    "P0-03": (
        ("hc_finest_abs_trapezoid_work_closure_error", "<="),
        ("hc_fitted_trapezoid_work_closure_slope_finest_four", ">="),
    ),
    "P0-04": (
        ("hc_max_crossed_field_equilibrium_momentum_departure", "<="),
    ),
    "P0-05": (
        ("boris_max_crossed_field_equilibrium_momentum_departure", ">="),
    ),
    "P0-06": (
        ("max_volume_preserving_abs_determinant_minus_one", "<="),
    ),
    "P0-07": (("vay_max_abs_determinant_minus_one", ">="),),
    "P0-08": (("hc_fitted_endpoint_error_slope_finest_four", ">="),),
    "P0-09": (
        ("hc_reversibility_roundtrip_abs_momentum_error", "<="),
    ),
    "P0-10": (
        ("hc_resonance_all_finite", "=="),
        ("hc_resonance_max_abs_h_error", "<="),
    ),
    "P0-11": (("hc_minus_vay_resonance_max_abs_h_error", "<"),),
    "P0-12": (("all_long_time_stress_finite", "=="),),
    "P0-13": (
        ("hc_uniform_b_finite", "=="),
        ("hc_uniform_b_max_rel_gamma_delta", "<="),
        ("hc_uniform_b_endpoint_phase_error_radians", "<="),
    ),
    "P0-14": (
        ("hc_crossed_eb_finite", "=="),
        ("hc_crossed_eb_average_velocity_error_vs_reference", "<="),
        ("hc_crossed_eb_average_velocity_error_vs_expected_drift", "<="),
    ),
    "P0-15": (
        ("hc_high_gamma_finite", "=="),
        ("hc_high_gamma_max_rel_gamma_delta", "<="),
        ("hc_high_gamma_endpoint_phase_error_radians", "<="),
        ("hc_high_gamma_endpoint_relative_momentum_error_vs_reference", "<="),
    ),
    "P0-16": (
        ("hc_near_exb_finite", "=="),
        ("hc_near_exb_max_rel_momentum_departure", "<="),
        ("strong_rotation_cancellation_guard_finite", "=="),
        ("strong_rotation_guard_used_cancellation_resistant_branch", "=="),
    ),
    "P0-17": (
        ("hc_direct_work_all_finite", "=="),
        ("hc_direct_work_max_abs_single_step_identity_residual", "<="),
        ("hc_direct_work_max_abs_accumulated_identity_residual", "<="),
    ),
}


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256_file(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


def _criterion_contract_sha256(criteria: List[Dict[str, object]]) -> str:
    contract = [
        [
            criterion.get("id"),
            criterion.get("requirement"),
            [
                [check.get("metric"), check.get("operator"), check.get("limit")]
                for check in criterion.get("checks", [])
            ],
        ]
        for criterion in criteria
    ]
    payload = json.dumps(contract, separators=(",", ":")).encode()
    return _sha256_bytes(payload)


def _load_preregistered_manifest(path: Path) -> Tuple[Dict[str, object], str]:
    payload = path.read_bytes()
    manifest_sha256 = _sha256_bytes(payload)
    if manifest_sha256 != EXPECTED_CRITERIA_MANIFEST_SHA256:
        raise ValueError(
            "criteria manifest bytes do not match the frozen harness digest")
    manifest = json.loads(payload)
    if not isinstance(manifest, dict):
        raise ValueError("criteria manifest must contain a JSON object")
    if manifest.get("manifest_schema_version") != 1:
        raise ValueError("criteria manifest schema version mismatch")
    if manifest.get("status") != "frozen_before_results_run":
        raise ValueError("criteria manifest is not frozen before results run")
    criteria = manifest.get("criteria")
    if not isinstance(criteria, list):
        raise ValueError("criteria manifest must contain a criteria list")
    contract_sha256 = _criterion_contract_sha256(criteria)
    if contract_sha256 != EXPECTED_CRITERIA_CONTRACT_SHA256:
        raise ValueError(
            "criteria manifest requirements or check thresholds do not match "
            "the frozen harness contract")
    criterion_ids = tuple(criterion.get("id") for criterion in criteria)
    if criterion_ids != tuple(EXPECTED_CRITERION_CHECKS):
        raise ValueError(
            "criteria manifest IDs or ordering do not match the harness")
    for criterion in criteria:
        if not isinstance(criterion.get("requirement"), str):
            raise ValueError(f"{criterion['id']} requirement must be text")
        checks = criterion.get("checks")
        if not isinstance(checks, list):
            raise ValueError(f"{criterion['id']} checks must be a list")
        actual_contract = tuple(
            (check.get("metric"), check.get("operator")) for check in checks)
        expected_contract = EXPECTED_CRITERION_CHECKS[criterion["id"]]
        if actual_contract != expected_contract:
            raise ValueError(f"{criterion['id']} metric contract mismatch")
        for check in checks:
            if "limit" not in check:
                raise ValueError(
                    f"{criterion['id']} check {check['metric']} has no limit")
    return manifest, manifest_sha256


def _manifest_unchanged(path: Path, expected_sha256: str) -> None:
    actual_sha256 = _sha256_file(path)
    if actual_sha256 != expected_sha256:
        raise RuntimeError(
            "criteria manifest changed during the results run: "
            f"expected {expected_sha256}, found {actual_sha256}")


def _paths_alias(first: Path, second: Path) -> bool:
    if first.resolve() == second.resolve():
        return True
    return first.exists() and second.exists() and first.samefile(second)


def _validate_distinct_artifact_paths(
        criteria_manifest: Path, output: Path, transcript: Path) -> None:
    pairs = (
        ("criteria manifest", criteria_manifest, "results JSON", output),
        ("criteria manifest", criteria_manifest, "transcript", transcript),
        ("results JSON", output, "transcript", transcript),
    )
    for first_label, first, second_label, second in pairs:
        if _paths_alias(first, second):
            raise ValueError(f"{first_label} and {second_label} must be distinct")


def _comparison_passed(actual: object, operator: str, limit: object) -> bool:
    if operator == "==":
        return bool(actual == limit)
    if operator == "<=":
        return bool(actual <= limit)
    if operator == ">=":
        return bool(actual >= limit)
    if operator == "<":
        return bool(actual < limit)
    raise ValueError(f"unsupported criteria operator: {operator}")


def _evaluate_acceptance(cases: Dict[str, object],
                         criteria: List[Dict[str, object]]
                         ) -> Dict[str, object]:
    pure_e = cases["pure_e_closed_form"]
    pure_e_errors = [
        pure_e[subcase]["methods"][method][
            "endpoint_momentum_error_vs_closed_form"]
        for subcase in ("acceleration", "deceleration")
        for method in METHODS
    ]
    low_speed = cases["low_speed_agreement"]
    work_closure = cases["work_closure"]["methods"]["higuera_cary"]
    work_rows = work_closure["rows"]
    cancellation = cases[
        "equivalent_boosted_frame_cancellation_discriminator"]["methods"]
    jacobian = cases["phase_volume_jacobian"]["methods"]
    max_volume_preserving_jacobian_error = max(
        jacobian[method]["max_abs_determinant_minus_one"]
        for method in VOLUME_PRESERVING_METHODS
    )
    convergence = cases["timestep_convergence"]["methods"]["higuera_cary"]
    reversibility = cases["reversibility"]["methods"]["higuera_cary"]
    resonance = cases["resonance_large_step_stress"]["methods"]
    hc_resonance = resonance["higuera_cary"]["summaries_by_dt"]["0.1"]
    vay_resonance = resonance["vay"]["summaries_by_dt"]["0.1"]
    stress = cases["long_time_stress"]["methods"]
    all_stress_finite = all(
        values[subcase]["finite"]
        for values in stress.values()
        for subcase in ("uniform_b", "crossed_field_equilibrium")
    )
    uniform_b = cases["uniform_b_energy_phase"]["methods"]["higuera_cary"]
    crossed_eb = cases["crossed_eb_drift"]["methods"]["higuera_cary"]
    high_gamma = cases["high_gamma"]["methods"]["higuera_cary"]
    near_limit = cases["finite_state_near_limit"]
    near_exb = near_limit["near_exb_limit"]["methods"]["higuera_cary"]
    rotation_guard = near_limit["strong_rotation_cancellation_guard"]
    direct_work = cases["higuera_cary_direct_work_identity"]
    actuals = {
        "P0-01": {
            "max_pure_e_endpoint_momentum_error": max(pure_e_errors),
        },
        "P0-02": {
            "max_low_speed_pairwise_endpoint_momentum_separation": low_speed[
                "max_pairwise_endpoint_momentum_separation"],
        },
        "P0-03": {
            "hc_finest_abs_trapezoid_work_closure_error": work_rows[-1][
                "abs_work_closure_error"],
            "hc_fitted_trapezoid_work_closure_slope_finest_four": work_closure[
                "fitted_abs_work_closure_slope_finest_four"],
        },
        "P0-04": {
            "hc_max_crossed_field_equilibrium_momentum_departure": cancellation[
                "higuera_cary"]["max_abs_momentum_departure"],
        },
        "P0-05": {
            "boris_max_crossed_field_equilibrium_momentum_departure": (
                cancellation["relativistic_boris"][
                    "max_abs_momentum_departure"]),
        },
        "P0-06": {
            "max_volume_preserving_abs_determinant_minus_one": (
                max_volume_preserving_jacobian_error),
        },
        "P0-07": {
            "vay_max_abs_determinant_minus_one": jacobian["vay"][
                "max_abs_determinant_minus_one"],
        },
        "P0-08": {
            "hc_fitted_endpoint_error_slope_finest_four": convergence[
                "fitted_endpoint_error_slope_finest_four"],
        },
        "P0-09": {
            "hc_reversibility_roundtrip_abs_momentum_error": reversibility[
                "roundtrip_abs_momentum_error"],
        },
        "P0-10": {
            "hc_resonance_all_finite": hc_resonance["all_finite"],
            "hc_resonance_max_abs_h_error": hc_resonance["max_abs_h_error"],
        },
        "P0-11": {
            "hc_minus_vay_resonance_max_abs_h_error": (
                hc_resonance["max_abs_h_error"]
                - vay_resonance["max_abs_h_error"]),
        },
        "P0-12": {
            "all_long_time_stress_finite": all_stress_finite,
        },
        "P0-13": {
            "hc_uniform_b_finite": uniform_b["finite"],
            "hc_uniform_b_max_rel_gamma_delta": uniform_b[
                "max_rel_gamma_delta"],
            "hc_uniform_b_endpoint_phase_error_radians": uniform_b[
                "endpoint_phase_error_radians"],
        },
        "P0-14": {
            "hc_crossed_eb_finite": crossed_eb["finite"],
            "hc_crossed_eb_average_velocity_error_vs_reference": crossed_eb[
                "average_velocity_error_vs_reference"],
            "hc_crossed_eb_average_velocity_error_vs_expected_drift": crossed_eb[
                "average_velocity_error_vs_expected_drift"],
        },
        "P0-15": {
            "hc_high_gamma_finite": high_gamma["finite"],
            "hc_high_gamma_max_rel_gamma_delta": high_gamma[
                "max_rel_gamma_delta"],
            "hc_high_gamma_endpoint_phase_error_radians": high_gamma[
                "endpoint_phase_error_radians"],
            "hc_high_gamma_endpoint_relative_momentum_error_vs_reference": (
                high_gamma[
                    "endpoint_relative_momentum_error_vs_reference"]),
        },
        "P0-16": {
            "hc_near_exb_finite": near_exb["finite"],
            "hc_near_exb_max_rel_momentum_departure": near_exb[
                "max_rel_momentum_departure"],
            "strong_rotation_cancellation_guard_finite": rotation_guard[
                "finite"],
            "strong_rotation_guard_used_cancellation_resistant_branch": (
                rotation_guard["details"]["evaluation"]
                == "conjugate_positive_root"),
        },
        "P0-17": {
            "hc_direct_work_all_finite": direct_work["all_finite"],
            "hc_direct_work_max_abs_single_step_identity_residual": direct_work[
                "max_abs_single_step_identity_residual"],
            "hc_direct_work_max_abs_accumulated_identity_residual": direct_work[
                "max_abs_accumulated_identity_residual"],
        },
    }
    results = []
    for criterion in criteria:
        criterion_id = criterion["id"]
        checks = criterion["checks"]
        actual_metrics = actuals[criterion_id]
        manifest_metrics = tuple(check["metric"] for check in checks)
        if tuple(actual_metrics) != manifest_metrics:
            raise ValueError(f"{criterion_id} evaluated metric contract mismatch")
        evaluated_checks = []
        for check in checks:
            actual = actual_metrics[check["metric"]]
            passed = _comparison_passed(
                actual, check["operator"], check["limit"])
            evaluated_checks.append({
                **check,
                "actual": actual,
                "passed": passed,
            })
        results.append({
            "id": criterion_id,
            "requirement": criterion["requirement"],
            "passed": all(check["passed"] for check in evaluated_checks),
            "checks": evaluated_checks,
        })
    return {
        "all_passed": all(result["passed"] for result in results),
        "results": results,
    }


def _to_builtin(value: object) -> object:
    if isinstance(value, np.ndarray):
        return [float(item) for item in value]
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, dict):
        return {str(key): _to_builtin(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_to_builtin(item) for item in value]
    return value


def _build_report(stress_steps: int, criteria_manifest_path: Path,
                  criteria_manifest: Dict[str, object],
                  criteria_manifest_sha256: str
                  ) -> Dict[str, object]:
    cases = {
        "pure_e_closed_form": _pure_e_closed_form(),
        "low_speed_agreement": _low_speed_agreement(),
        "work_closure": _work_closure(),
        "higuera_cary_direct_work_identity": _higuera_cary_direct_work_identity(),
        "phase_volume_jacobian": _phase_volume_jacobian(),
        "uniform_b_energy_phase": _uniform_b_energy_phase(),
        "crossed_eb_drift": _crossed_eb_drift(),
        "equivalent_boosted_frame_cancellation_discriminator": (
            _crossed_field_equilibrium()),
        "high_gamma": _high_gamma(),
        "timestep_convergence": _timestep_convergence(),
        "reversibility": _reversibility(),
        "finite_state_near_limit": _finite_state_near_limit(),
        "resonance_large_step_stress": _resonance_large_step_stress(),
        "long_time_stress": _long_time_stress(stress_steps),
    }
    criteria = criteria_manifest["criteria"]
    acceptance_evaluation = _evaluate_acceptance(cases, criteria)
    observations = [
        (
            "All compared methods conserve uniform-B energy to roundoff in this "
            "prescribed-field experiment."
        ),
        (
            "The equilibrium-cancellation case directly exposes the "
            "crossed-field limiting-solution difference discussed in the paper."
        ),
        (
            "The convergence fit checks the expected second-order endpoint "
            "behavior against an independently integrated Lorentz ODE."
        ),
        (
            "The finite-state stress includes a cancellation-resistant "
            "evaluation of the HC Eq. (20) positive root within bounded spike "
            "inputs.  Production helpers must own overflow-resistant scaled "
            "norms."
        ),
        (
            "The exact discrete Higuera-Cary direct-work identity is checked "
            "separately from the trapezoid-work convergence diagnostic."
        ),
        (
            "The Vay comparator is included directly from the primary source. "
            "Its Jacobian and resonance probes retain the bounded reason not to "
            "select it for the Phase 0 production path."
        ),
    ]
    return {
        "preregistered_acceptance_criteria": {
            "manifest_path": str(criteria_manifest_path),
            "manifest_sha256": criteria_manifest_sha256,
            "loaded_before_results": True,
            "manifest_sha256_verified_unchanged_before_emission": False,
            "note": (
                "The executable loads and prints the frozen criteria manifest "
                "before calculating case results, then refuses to emit results "
                "if the manifest bytes changed during calculation."),
            "manifest": criteria_manifest,
        },
        "experiment": "CR relativistic pusher Phase 0 comparison spike",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(
            timespec="seconds"),
        "scope": (
            "Standalone prescribed-field Python comparison only; no AthenaK "
            "production integration."),
        "normalization": {
            "speed_of_light_c": 1.0,
            "charge_to_mass_q_over_m": 1.0,
            "momentum_variable": "u = gamma*v",
            "gamma": "sqrt(1 + |u|^2)",
        },
        "phase_0_selection_boundary": {
            "compared_methods": METHODS,
            "working_selection": PHASE_0_WORKING_SELECTION,
            "production_authorization": False,
            "reason": (
                "Higuera-Cary is the compared method that preserves both the "
                "crossed-field equilibrium cancellation and differential "
                "volume.  Vay is retained as an experimental comparator; "
                "standard relativistic Boris is retained as a control."),
        },
        "environment": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__ if scipy is not None else None,
            "reference_policy": (
                "scipy.integrate.solve_ivp(method=DOP853) with tight tolerances"
                if solve_ivp is not None else
                "independent fixed-step RK4 fallback because scipy is absent"
            ),
            "git": _git_metadata(),
        },
        "operation_and_clarity_notes": _operation_and_clarity_notes(),
        "cases": cases,
        "acceptance_evaluation": acceptance_evaluation,
        "observations": observations,
        "residual_limits": _residual_limits(),
    }


def _summary_lines(report: Dict[str, object], output: Path,
                   transcript: Path, output_sha256: str) -> List[str]:
    cases = report["cases"]
    equilibrium = cases[
        "equivalent_boosted_frame_cancellation_discriminator"]["methods"]
    high_gamma = cases["high_gamma"]["methods"]
    convergence = cases["timestep_convergence"]["methods"]
    reversibility = cases["reversibility"]["methods"]
    jacobian = cases["phase_volume_jacobian"]["methods"]
    resonance = cases["resonance_large_step_stress"]["methods"]
    stress = cases["long_time_stress"]["methods"]
    direct_work = cases["higuera_cary_direct_work_identity"]
    lines = [
        f"Wrote results JSON: {output}",
        f"Results JSON sha256: {output_sha256}",
        f"Transcript path: {transcript}",
        (
            "Criteria manifest sha256 verified unchanged before emission: "
            f"{report['preregistered_acceptance_criteria']['manifest_sha256']}"
        ),
        "Crossed-field equilibrium max |delta u|:",
    ]
    for method in METHODS:
        lines.append(
            f"  {method}: "
            f"{equilibrium[method]['max_abs_momentum_departure']:.6e}")
    lines.append("High-gamma endpoint phase error [rad]:")
    for method in METHODS:
        lines.append(
            f"  {method}: "
            f"{high_gamma[method]['endpoint_phase_error_radians']:.6e}")
    lines.append("Timestep-convergence fitted slope, finest four points:")
    for method in METHODS:
        lines.append(
            f"  {method}: "
            f"{convergence[method]['fitted_endpoint_error_slope_finest_four']:.6f}")
    lines.append("Reversibility roundtrip |delta u|:")
    for method in METHODS:
        lines.append(
            f"  {method}: "
            f"{reversibility[method]['roundtrip_abs_momentum_error']:.6e}")
    lines.append("Momentum-map Jacobian max |det(J) - 1|:")
    for method in METHODS:
        lines.append(
            f"  {method}: "
            f"{jacobian[method]['max_abs_determinant_minus_one']:.6e}")
    lines.append("Resonance scan dt=0.1 max |delta H|:")
    for method in METHODS:
        lines.append(
            f"  {method}: "
            f"{resonance[method]['summaries_by_dt']['0.1']['max_abs_h_error']:.6e}")
    lines.append("Long-time stress uniform-B max relative gamma drift:")
    for method in METHODS:
        lines.append(
            f"  {method}: "
            f"{stress[method]['uniform_b']['max_rel_gamma_delta']:.6e}")
    lines.extend([
        "HC exact discrete direct-work identity residuals:",
        (
            "  max abs single-step residual: "
            f"{direct_work['max_abs_single_step_identity_residual']:.6e}"
        ),
        (
            "  max abs accumulated residual: "
            f"{direct_work['max_abs_accumulated_identity_residual']:.6e}"
        ),
    ])
    acceptance = report["acceptance_evaluation"]
    lines.append(f"Acceptance criteria all passed: {acceptance['all_passed']}")
    for result in acceptance["results"]:
        status = "PASS" if result["passed"] else "FAIL"
        lines.append(f"  {result['id']}: {status}")
        for check in result["checks"]:
            lines.append(
                f"    {check['metric']}: actual={check['actual']!r} "
                f"{check['operator']} limit={check['limit']!r} -> "
                f"{check['passed']}")
    return lines


def _preregistered_manifest_lines(
        manifest_path: Path, manifest_sha256: str,
        criteria: List[Dict[str, object]]) -> List[str]:
    lines = [
        "CR relativistic pusher Phase 0 comparison spike transcript",
        f"Frozen criteria manifest: {manifest_path}",
        f"Frozen criteria manifest sha256: {manifest_sha256}",
        "Preregistered acceptance criteria loaded before results calculation:",
    ]
    for criterion in criteria:
        lines.append(f"  {criterion['id']}: {criterion['requirement']}")
        for check in criterion["checks"]:
            lines.append(
                f"    {check['metric']} {check['operator']} {check['limit']!r}")
    return lines


def _print_lines(lines: Sequence[str]) -> None:
    for line in lines:
        print(line)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=DEFAULT_OUTPUT,
        help="JSON evidence path (default: %(default)s)")
    parser.add_argument(
        "--criteria-manifest", type=Path, default=DEFAULT_CRITERIA_MANIFEST,
        help="frozen criteria manifest path (default: %(default)s)")
    parser.add_argument(
        "--transcript", type=Path, default=DEFAULT_TRANSCRIPT,
        help="run transcript path (default: %(default)s)")
    parser.add_argument(
        "--stress-steps", type=int, default=100000,
        help="steps per long-time stress subcase (default: %(default)s)")
    args = parser.parse_args()
    if args.stress_steps <= 0:
        parser.error("--stress-steps must be positive")

    output = args.output.resolve()
    criteria_manifest_path = args.criteria_manifest.resolve()
    transcript = args.transcript.resolve()
    try:
        _validate_distinct_artifact_paths(
            criteria_manifest_path, output, transcript)
    except ValueError as exc:
        parser.error(str(exc))
    try:
        criteria_manifest, criteria_manifest_sha256 = (
            _load_preregistered_manifest(criteria_manifest_path))
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        parser.error(f"cannot load frozen criteria manifest: {exc}")
    prelude_lines = _preregistered_manifest_lines(
        criteria_manifest_path, criteria_manifest_sha256,
        criteria_manifest["criteria"])
    _print_lines(prelude_lines)

    report = _to_builtin(_build_report(
        args.stress_steps, criteria_manifest_path, criteria_manifest,
        criteria_manifest_sha256))
    _manifest_unchanged(criteria_manifest_path, criteria_manifest_sha256)
    report["preregistered_acceptance_criteria"][
        "manifest_sha256_verified_unchanged_before_emission"] = True
    report["transcript_path"] = str(transcript)
    output_payload = json.dumps(report, indent=2) + "\n"
    output_sha256 = _sha256_bytes(output_payload.encode())
    summary_lines = _summary_lines(report, output, transcript, output_sha256)
    transcript_payload = "\n".join(prelude_lines + [""] + summary_lines) + "\n"
    _manifest_unchanged(criteria_manifest_path, criteria_manifest_sha256)
    try:
        _validate_distinct_artifact_paths(
            criteria_manifest_path, output, transcript)
    except ValueError as exc:
        parser.error(str(exc))
    output.parent.mkdir(parents=True, exist_ok=True)
    transcript.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(output_payload)
    transcript.write_text(transcript_payload)
    _print_lines(summary_lines)
    return 0 if report["acceptance_evaluation"]["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
