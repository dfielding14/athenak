"""Square-pulse physics validation for MC and Ito-2 flux tracers."""

from pathlib import Path
import shutil
import sys

import numpy as np

import test_suite.testutils as testutils


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
from read_prtcl_thermo_history import read_history  # noqa: E402


INPUT = str(ROOT / "tst/inputs/particles_ito_square_pulse.athinput")
RUN_ROOT = Path("run_particles_ito_square_pulse_cpu")
VELOCITY = 1.0


def _sheet_bounds(nx1):
    cell = int(0.375 * nx1)
    center = (cell + 0.5) / nx1
    half_width = 0.1 / nx1
    return center - half_width, center + half_width


def _run(method, nx1, npart, label):
    run_dir = RUN_ROOT / label
    shutil.rmtree(run_dir, ignore_errors=True)
    run_dir.mkdir(parents=True)
    basename = f"ito_square_pulse_{label}"
    slab_min, slab_max = _sheet_bounds(nx1)
    flags = [
        "-d",
        str(run_dir),
        f"job/basename={basename}",
        f"mesh/nx1={nx1}",
        f"meshblock/nx1={nx1 // 2}",
        f"tracer_seed1/count_per_event={npart}",
        f"tracer_seed1/slab_min={slab_min:.17g}",
        f"tracer_seed1/slab_max={slab_max:.17g}",
    ]
    if method == "mc":
        flags.extend(
            [
                "particles/particle_type=lagrangian_mc",
                "particles/pusher=lagrangian_mc",
            ]
        )
    else:
        flags.extend(
            [
                "particles/particle_type=lagrangian_ito",
                "particles/pusher=ito2",
            ]
        )
    assert testutils.run(INPUT, flags)
    history = read_history(
        run_dir
        / "prtcl_thermo_history"
        / f"{basename}.prtcl_thermo_history.thp"
    )
    return _displacements(history), float(np.max(history["time"]))


def _displacements(history):
    first_cycle = int(np.min(history["cycle"]))
    last_cycle = int(np.max(history["cycle"]))
    assert first_cycle == 0
    assert last_cycle > first_cycle

    initial = history["cycle"] == first_cycle
    final = history["cycle"] == last_cycle
    initial_order = np.argsort(history["tag"][initial])
    final_order = np.argsort(history["tag"][final])
    initial_tags = history["tag"][initial][initial_order]
    final_tags = history["tag"][final][final_order]
    np.testing.assert_array_equal(initial_tags, final_tags)

    displacement = (
        history["x1"][final][final_order]
        - history["x1"][initial][initial_order]
    )
    return (displacement + 0.5) % 1.0 - 0.5


def _moments(samples):
    mean = float(np.mean(samples))
    centered = samples - mean
    variance = float(np.mean(centered**2))
    width = np.sqrt(variance)
    standardized = centered / width
    return {
        "mean": mean,
        "variance": variance,
        "width": width,
        "skewness": float(np.mean(standardized**3)),
        "excess_kurtosis": float(np.mean(standardized**4) - 3.0),
        "standardized": standardized,
    }


def _standardized_wasserstein(first, second):
    assert first.size == second.size
    return float(np.mean(np.abs(np.sort(first) - np.sort(second))))


def _assert_lower_moment_agreement(mc, ito, elapsed):
    expected_mean = VELOCITY * elapsed
    mean_se = np.sqrt(
        mc["variance"] / mc["standardized"].size
        + ito["variance"] / ito["standardized"].size
    )
    expected_se_mc = np.sqrt(mc["variance"] / mc["standardized"].size)
    expected_se_ito = np.sqrt(ito["variance"] / ito["standardized"].size)
    assert abs(mc["mean"] - expected_mean) <= 5.0 * expected_se_mc
    assert abs(ito["mean"] - expected_mean) <= 5.0 * expected_se_ito
    assert abs(mc["mean"] - ito["mean"]) <= 5.0 * mean_se

    mc_squared = mc["standardized"] ** 2 * mc["variance"]
    ito_squared = ito["standardized"] ** 2 * ito["variance"]
    variance_se = np.sqrt(
        np.var(mc_squared, ddof=1) / mc_squared.size
        + np.var(ito_squared, ddof=1) / ito_squared.size
    )
    assert abs(mc["variance"] - ito["variance"]) <= 5.0 * variance_se


def _print_comparison(label, mc, ito, shape_distance):
    print(f"{label}:")
    for name in (
        "mean",
        "variance",
        "width",
        "skewness",
        "excess_kurtosis",
    ):
        print(
            f"  {name}: MC={mc[name]:.8e}, Ito-2={ito[name]:.8e}, "
            f"delta={ito[name] - mc[name]:+.3e}"
        )
    print(f"  standardized Wasserstein distance={shape_distance:.8e}")


def test_square_pulse_mc_ito2_lower_moments_and_pdf_shape_cpu():
    """MC and Ito-2 agree in sheet drift/width, but not in full PDF shape."""
    cases = ((32, 8192), (32, 32768), (64, 32768))
    results = {}
    try:
        for nx1, npart in cases:
            label = f"nx{nx1}_n{npart}"
            mc_samples, mc_time = _run("mc", nx1, npart, f"{label}_mc")
            ito_samples, ito_time = _run("ito2", nx1, npart, f"{label}_ito2")
            assert mc_time == ito_time
            assert mc_samples.size == npart
            assert ito_samples.size == npart

            mc = _moments(mc_samples)
            ito = _moments(ito_samples)
            _assert_lower_moment_agreement(mc, ito, mc_time)
            shape_distance = _standardized_wasserstein(
                mc["standardized"], ito["standardized"]
            )
            _print_comparison(label, mc, ito, shape_distance)

            assert 3.0 / np.sqrt(npart) < shape_distance < 0.25
            results[(nx1, npart)] = (mc, ito, shape_distance)

        low_shape = results[(32, 8192)][2]
        high_shape = results[(32, 32768)][2]
        assert abs(low_shape - high_shape) / high_shape < 0.15

        for method_index in (0, 1):
            coarse_variance = results[(32, 32768)][method_index]["variance"]
            fine_variance = results[(64, 32768)][method_index]["variance"]
            coarse_scaled = coarse_variance * 32
            fine_scaled = fine_variance * 64
            assert abs(coarse_scaled - fine_scaled) / coarse_scaled < 0.05

        assert results[(64, 32768)][2] < high_shape
    finally:
        shutil.rmtree(RUN_ROOT, ignore_errors=True)
