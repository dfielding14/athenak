# Regression test for scalar_mixing_perfect_powerlaw velocity initialization.
#
# This test is only active when AthenaK is configured with
#   --cmake=-DPROBLEM=scalar_mixing_perfect_powerlaw
#
# It validates the exact-shell contract directly from the centered-shell raw
# velocity spectrum used by scalar_mixing_stream_run_analysis_spectrum_fix3.py.

import glob
import logging
import os
import sys

import numpy as np

import scripts.utils.athena as athena

sys.path.insert(0, "../vis/python")
import bin_convert_new  # noqa
import scalar_mixing_stream_run_analysis_spectrum_fix3 as spectrum_fix3  # noqa

logger = logging.getLogger("athena" + __name__[7:])

_enabled = False
_SEED = 12345

_CASES = [
    {
        "name": "stream_2d_exact_shell",
        "input": "tests/scalar_mixing_stream_2d_perfect_powerlaw_regression.athinput",
        "basename": "scalar_mix_stream_2d_perfect_powerlaw",
        "nlow": 2,
        "nhigh": 6,
        "expo": 1.6666667,
        "vrms": 1.0,
        "expect_vz_zero": True,
        "shell_tol": 5.0e-7,
        "div_tol": 2.0e-7,
        "leak_tol": 5.0e-8,
    },
    {
        "name": "projection_3d_exact_shell",
        "input": "tests/scalar_mixing_projection_3d_perfect_powerlaw_regression.athinput",
        "basename": "scalar_mix_proj_3d_perfect_powerlaw",
        "nlow": 2,
        "nhigh": 6,
        "expo": 1.6666667,
        "vrms": 1.0,
        "expect_vz_zero": False,
        "shell_tol": 2.0e-6,
        "div_tol": 2.0e-7,
        "leak_tol": 5.0e-8,
    },
]


def _compiled_problem():
    cfg_path = os.path.join("build", "config.hpp")
    if not os.path.exists(cfg_path):
        return None
    with open(cfg_path, "r", encoding="utf-8") as fobj:
        for line in fobj:
            if line.startswith("#define PROBLEM_GENERATOR"):
                return line.split('"')[1]
    return None


def _run_basename(case):
    return f"{case['basename']}_{_SEED}"


def _output_path(case):
    pattern = os.path.join("build", "src", "bin", _run_basename(case) + ".hydro_w.*.bin")
    matches = sorted(glob.glob(pattern))
    if not matches:
        raise RuntimeError("No hydro_w output found for " + case["name"])
    return matches[-1]


def _clear_old_outputs(case):
    pattern = os.path.join("build", "src", "bin", _run_basename(case) + ".hydro_w.*.bin")
    for path in glob.glob(pattern):
        os.remove(path)


def _shell_bins(shape):
    nz, ny, nx = shape
    kx = np.fft.fftfreq(nx) * nx
    ky = np.fft.fftfreq(ny) * ny if ny > 1 else np.array([0.0])
    kz = np.fft.fftfreq(nz) * nz if nz > 1 else np.array([0.0])
    kz3, ky3, kx3 = np.meshgrid(kz, ky, kx, indexing="ij")
    k2 = kx3**2 + ky3**2 + kz3**2
    return kx3, ky3, kz3, k2


def _best_fit_shell_scale(kvals, shell_sum, expo):
    target = np.asarray(kvals, dtype=np.float64) ** (-expo)
    return float(np.dot(shell_sum, target) / np.dot(target, target))


def _max_shell_relative_error(kvals, shell_sum, expo):
    scale = _best_fit_shell_scale(kvals, shell_sum, expo)
    target = scale * (np.asarray(kvals, dtype=np.float64) ** (-expo))
    return float(np.max(np.abs(shell_sum - target) / target))


def _analyze_output(case):
    data = bin_convert_new.read_binary_as_athdf(
        _output_path(case), quantities=["velx", "vely", "velz"]
    )
    vx = spectrum_fix3.base._canonical_field(np.asarray(data["velx"], dtype=np.float64))
    vy = spectrum_fix3.base._canonical_field(np.asarray(data["vely"], dtype=np.float64))
    vz = spectrum_fix3.base._canonical_field(np.asarray(data["velz"], dtype=np.float64))

    means = np.array([vx.mean(), vy.mean(), vz.mean()])
    vrms = float(np.sqrt(np.mean(vx * vx + vy * vy + vz * vz)))

    if vx.ndim == 2:
        kvals, shell_sum, shell_count = spectrum_fix3._velocity_spectrum_2d(vx, vy, vz)
    else:
        kvals, shell_sum, shell_count = spectrum_fix3._velocity_spectrum_3d(vx, vy, vz)

    mask = (kvals >= case["nlow"]) & (kvals <= case["nhigh"])
    shell_rel_err = _max_shell_relative_error(kvals[mask], shell_sum[mask], case["expo"])

    leakage_num = float(np.sum(shell_sum[(kvals < case["nlow"]) | (kvals > case["nhigh"])]))
    leakage_den = float(np.sum(shell_sum))
    leakage = leakage_num / leakage_den if leakage_den > 0.0 else 0.0

    uhat_x = np.fft.fftn(vx)
    uhat_y = np.fft.fftn(vy)
    uhat_z = np.fft.fftn(vz)
    kx3, ky3, kz3, k2 = _shell_bins((1, *vx.shape) if vx.ndim == 2 else vx.shape)
    nonzero = k2 > 0.0
    kdot_u = kx3 * uhat_x + ky3 * uhat_y + kz3 * uhat_z
    energy = np.abs(uhat_x) ** 2 + np.abs(uhat_y) ** 2 + np.abs(uhat_z) ** 2
    div_num = np.sum(np.abs(kdot_u[nonzero]) ** 2)
    div_den = np.sum(k2[nonzero] * energy[nonzero])
    div_ratio = float(np.sqrt(div_num / div_den)) if div_den > 0.0 else 0.0

    return {
        "means": means,
        "vrms": vrms,
        "div_ratio": div_ratio,
        "leakage": leakage,
        "shell_rel_err": shell_rel_err,
        "shell_count": shell_count,
        "vz_max": float(np.max(np.abs(vz))),
    }


def run(**kwargs):
    del kwargs
    global _enabled
    problem = _compiled_problem()
    _enabled = problem == "scalar_mixing_perfect_powerlaw"
    if not _enabled:
        logger.info(
            "Skipping %s because build/config.hpp reports PROBLEM=%s",
            __name__,
            problem,
        )
        return

    for case in _CASES:
        _clear_old_outputs(case)
        athena.run(
            case["input"],
            [
                f"job/basename={_run_basename(case)}",
                f"problem/turb_rseed={_SEED}",
                "problem/turb_spectrum_contract=exact_shell",
            ],
        )


def analyze():
    if not _enabled:
        return True

    analyze_status = True
    for case in _CASES:
        metrics = _analyze_output(case)
        logger.info(
            "%s: shell_rel_err=%.3e max_div=%.3e leakage=%.3e",
            case["name"],
            metrics["shell_rel_err"],
            metrics["div_ratio"],
            metrics["leakage"],
        )

        max_mean_err = float(np.max(np.abs(metrics["means"])))
        if max_mean_err > 1.0e-9:
            logger.warning("%s mean velocity too large: %.3e", case["name"], max_mean_err)
            analyze_status = False

        vrms_rel_err = abs(metrics["vrms"] - case["vrms"]) / case["vrms"]
        if vrms_rel_err > 1.0e-3:
            logger.warning("%s vrms mismatch: %.6e", case["name"], vrms_rel_err)
            analyze_status = False

        if metrics["div_ratio"] > case["div_tol"]:
            logger.warning(
                "%s divergence ratio too large: %.3e", case["name"], metrics["div_ratio"]
            )
            analyze_status = False

        if metrics["shell_rel_err"] > case["shell_tol"]:
            logger.warning(
                "%s raw centered-shell spectrum mismatch: %.3e",
                case["name"],
                metrics["shell_rel_err"],
            )
            analyze_status = False

        if metrics["leakage"] > case["leak_tol"]:
            logger.warning("%s leakage too large: %.3e", case["name"], metrics["leakage"])
            analyze_status = False

        if case["expect_vz_zero"] and metrics["vz_max"] > 1.0e-12:
            logger.warning(
                "%s expected velz=0 but max(|vz|)=%.3e", case["name"], metrics["vz_max"]
            )
            analyze_status = False

    return analyze_status
