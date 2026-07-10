from __future__ import annotations

import math
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def _function_body(source: str, signature: str) -> str:
    start = source.index(signature)
    brace = source.index("{", start)
    depth = 0
    for index in range(brace, len(source)):
        if source[index] == "{":
            depth += 1
        elif source[index] == "}":
            depth -= 1
            if depth == 0:
                return source[brace:index + 1]
    raise AssertionError(f"unterminated function: {signature}")


def test_particle_timestep_uses_species_bound_for_empty_packs() -> None:
    source = (REPO_ROOT / "src/particles/particles.cpp").read_text(
        encoding="ascii"
    )
    body = _function_body(source, "void Particles::NewTimeStep()")

    assert "if (nprtcl_thispack <= 0)" not in body
    assert "const Real qom_max = pic_species_qom_max;" in body
    assert "ParticlesNewTimeStepQom" not in body


def test_particle_gyro_bound_includes_ghost_field_storage() -> None:
    source = (REPO_ROOT / "src/particles/particles.cpp").read_text(
        encoding="ascii"
    )
    body = _function_body(source, "void Particles::NewTimeStep()")

    assert "const int nx1 = bcc.extent_int(4);" in body
    assert "const int nx2 = bcc.extent_int(3);" in body
    assert "const int nx3 = bcc.extent_int(2);" in body
    assert "const int nx1 = indcs.nx1;" not in body


def test_every_pic_feedback_path_checks_post_source_admissibility() -> None:
    source = (REPO_ROOT / "src/mhd/mhd_tasks.cpp").read_text(encoding="ascii")
    expected_calls = {
        'ValidatePICFeedbackState("mhd_src_terms", stage);',
        'ValidatePICFeedbackState("efield_src", stage);',
        'ValidatePICFeedbackState("expanding_box_feedback", stage);',
    }
    for call in expected_calls:
        assert source.count(call) == 1

    guard = _function_body(source, "void MHD::ValidatePICFeedbackState(")
    for diagnostic in (
        "nonfinite=",
        "density_floor=",
        "energy_pressure_floor=",
        "temperature_floor=",
        "entropy_floor=",
    ):
        assert diagnostic in guard
    assert "restart_utils::AbortOnFatalError();" in guard


def test_parallel_shock_limits_first_cohort_before_transaction() -> None:
    source = (
        REPO_ROOT / "src/pgen/tests/pic_parallel_shock.cpp"
    ).read_text(encoding="ascii")
    callback = _function_body(source, "void ParallelShockWorkBeforeLoop(Mesh *pm)")
    limiter = _function_body(
        source, "void LimitParallelShockInjectionTimeStep(Mesh *pm)"
    )

    assert callback.index("LimitParallelShockInjectionTimeStep(pm);") < callback.index(
        "PrepareParallelShockInjectionTransaction(pm);"
    )
    assert "ppart->pic_max_cell_cross" in limiter
    assert "injection_dt *= pm->cfl_no;" in limiter
    assert "MPI_Allreduce" in limiter
    assert "BoostRelativeVelocityFromSurface(" in limiter
    assert "const Real boost_product" in limiter
    assert "vmax_x1 = ppart->pic_cr_light_speed" not in limiter


def test_relativistic_injection_component_envelope_is_sharp() -> None:
    light_speed = 10.0
    for surface_vx, vinj in ((0.4, 2.0), (4.0, 3.0), (-3.0, 7.0)):
        gamma_surface = 1.0 / math.sqrt(
            1.0 - (surface_vx / light_speed) ** 2
        )
        boost_product = surface_vx * vinj / light_speed**2
        transverse_bound = vinj / (
            gamma_surface * math.sqrt(1.0 - boost_product**2)
        )

        def boosted_x(relative_vx: float) -> float:
            denominator = 1.0 + surface_vx * relative_vx / light_speed**2
            return (surface_vx + relative_vx) / denominator

        longitudinal_bound = max(
            abs(boosted_x(vinj)), abs(boosted_x(-vinj))
        )
        sampled_transverse_max = 0.0
        for index in range(4001):
            mu = -1.0 + 2.0 * index / 4000.0
            denominator = 1.0 + boost_product * mu
            transverse = (
                vinj
                * math.sqrt(max(0.0, 1.0 - mu * mu))
                / (gamma_surface * denominator)
            )
            sampled_transverse_max = max(sampled_transverse_max, transverse)
            assert abs(boosted_x(vinj * mu)) <= longitudinal_bound * (1.0 + 1e-14)
            assert transverse <= transverse_bound * (1.0 + 1e-14)
        assert math.isclose(
            sampled_transverse_max, transverse_bound, rel_tol=2.0e-7
        )
