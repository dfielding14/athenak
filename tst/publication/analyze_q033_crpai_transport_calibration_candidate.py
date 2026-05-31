#!/usr/bin/env python3
"""Q-033 source-local synthetic CRPAI transport-contract analyzer."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from statistics import fmean
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN_ID = "Q033-EXT-CRPAI-TRANSPORT-CALIBRATION-CANDIDATE"
ARTIFACT_ROLE = "source_local_synthetic_transport_contract_only"
QUALIFICATION_EFFECT = "none_source_local_launch_blocked_contract_only"
LAUNCH_STATUS = (
    "blocked_missing_q033_generator_applicability_registered_frontier_matrix_"
    "and_external_review"
)
DECK = REPO_ROOT / "inputs/tests/pic_q033_crpai_transport_calibration_candidate.athinput"
DECK_SHA256 = "77ffeb7bdfd1f1f2299db7e4ee257d105b580a6bd554cd244e6d9e0c8214ce58"
_PGEN_DISPATCH = REPO_ROOT / "src/pgen/pgen.cpp"
_PGEN_NAME = "q033_crpai_transport_calibration_open"
ABSOLUTE_TOLERANCE = 1.0e-12
QUASI_STEADY_RELATIVE_SPAN_MAX = 1.0e-12
DRIVING_RATE_GRID = (1.0, 4.0)
NU_IN_GRID = (1.0, 9.0)
TAIL_TIMES = (1.0, 2.0, 3.0, 4.0)
SPECTRUM_K_GRID = (1.0, 2.0, 3.0)

_SOURCE_BINDINGS = {
    "src/particles/particles.cpp":
        "4b88b1b0a0d83aed824f2af6f2fb77fc5793cdeb9b0e2d65d408535f22e45c44",
    "src/particles/particles.hpp":
        "452e6091ad6f1ce50156c770d304ac36020cb53927cc1e40ba6af73e7c3bb462",
    "src/particles/particles_tasks.cpp":
        "8310e80ffdb2739ecc60b6d65ef9eae78573aeba8f55fa1981ff98b4801c8946",
    "src/particles/particles_pushers.cpp":
        "7cd3ac4b6fc95529d3097a3799b74f064a7e865c8b7be08f830a60dbee103330",
    "src/particles/particles_moments.cpp":
        "3b206b1268bc71a0d6745dcda5422da54d255ac6cb73fddd566b89b507da6c58",
    "src/mhd/mhd.hpp":
        "655197f1582463020747212294e965e19a078a59d0708ebf7a6a127e5b652d4e",
    "src/mhd/mhd_tasks.cpp":
        "228144855b0925e6405d14e33d22dd90381a01a6f317049bf3d6cd041a99530f",
    "docs/source/engineering/pic_mhd_model_contract.md":
        "de319a56b34a6be33a41e9da3719d9a06a7646a1ff26748900d746d1b167bf80",
}
_SOURCE_REQUIRED_SNIPPETS = {
    "src/particles/particles.cpp": (
        "global_bikappa_moments_experimental requires physical ",
        "reduced ion-neutral friction without an admitted ",
        "cell-centered source split",
        '"the Newtonian single-fluid MHD source task path; does not support "',
        '"<ion-neutral> alternate task lists"',
    ),
    "src/particles/particles.hpp": (
        "Real PICAdaptiveDeltaFBackgroundValue(",
        "return pic_wave_damping_mode == PICWaveDampingMode::ion_neutral_friction;",
        "pic_deltaf_adaptive_xi, pic_deltaf_adaptive_p0,",
    ),
    "src/particles/particles_tasks.cpp": (
        "TaskStatus Particles::AdaptDeltaF(Driver *pdriver, int stage)",
        "PIC adaptive delta-f fit: time=",
    ),
    "src/particles/particles_pushers.cpp": (
        "const bool adaptive_deltaf_local = UsesAdaptiveDeltaF();",
        "PICAdaptiveDeltaFBackgroundValue(",
    ),
    "src/particles/particles_moments.cpp": (
        "physical_density_scale = geom.inv_a1*geom.inv_a2*geom.inv_a3;",
        "const Real df_weight = deltaf_local ? pr(IPDFWT,p)",
    ),
    "src/mhd/mhd_tasks.cpp": (
        "Apply reduced wave damping after final expanding-box physical-frame sources.",
        "const Real factor = exp(-ppart->pic_ion_neutral_collision_rate*",
        "const Real background_rho =",
    ),
    "docs/source/engineering/pic_mhd_model_contract.md": (
        "p_ion,perp(t + dt) = p_ion,perp(t) exp(-nu_in dt)",
        "It does not establish",
        "CRPAI transport calibration, a saturated-state scattering rate, or GPU/MPI",
    ),
}
_BOUNDARY_ARTIFACT_BINDINGS = {
    "inputs/tests/pic_mhd_expanding_box_adaptive_damping_smoke.athinput":
        "52122a6100eba205191abea2c0f560f26c30f5ca97272b7fd7dd1cfc01832ba5",
    "tst/scripts/particles/pic_mhd_expanding_box_adaptive_damping_smoke.py":
        "dd049fd2e0828900f9b7152562fd032612647766bda54907c38e5d79e15d46e3",
    "tst/scripts/particles/pic_mhd_expanding_box_adaptive_damping_restart.py":
        "933d7965ae116bb4eb7fd171122e322cb63d3878b4d079a0efc12c43e69513c3",
    "tst/publication/readiness/q032_q033_extension_local_scans_2026-05-30.json":
        "8338ee08f311ce353186a514549cec1227d17de300557645655158ee73da3dc4",
    "tst/publication/readiness/"
    "q033_expanding_box_adaptive_damping_restart_resilience_bounded_local_2026-05-30.json":
        "75e96ceb330d419f4fe7d60331df2bb1fd9e791e2a482778bbe0978f541a342c",
    "tst/publication/readiness/q008_q033_expanding_coupled_source_successor_2026-05-30.json":
        "c68a5c9ff78b5652f4bb43a90e52a376a8f7a3412bcbe5627e4868cd4bc0968e",
}
_Q022_PREREQUISITE_BINDINGS = {
    "tst/publication/readiness/q022_external_reference_private_ingest_2026-05-30.json":
        "78287ea54d8350445cfad62e72f6efc2845f720580d070d075897d889ae3a214",
    "tst/publication/readiness/q022_dataset_provenance_manifest_2026-05-30.json":
        "925ab235f4ed5141468a8643e999d37e9d1a075c528df1987b577e2deb90bc27",
    "tst/publication/readiness/q022_xcmp_ext_crpai_transport_equation_map_2026-05-30.json":
        "584cbee640b15e1e67a18d5652012e3280bfcbd4bf81e9a6a9bca873b114c786",
    "tst/publication/readiness/q022_xcmp_ext_crpai_transport_tolerance_table_2026-05-30.json":
        "de47d2b7f6d75f21b92aaf0806596b8e444992ef1ae568d50366d62ee730dd81",
}
_Q022_COMPARISON_ID = "XCMP-EXT-CRPAI-TRANSPORT"
_Q022_DATASET_ID = "Q022-DATASET-XCMP-EXT-CRPAI-TRANSPORT"
_Q022_REFERENCE_ID = "sun_bai_zhao_2024_arxiv_2409.08592"
_EXPECTED_DECK_VALUES = {
    ("time", "nlim"): "0",
    ("time", "tlim"): "0.0",
    ("particles", "particle_type"): "cosmic_ray",
    ("particles", "pusher"): "boris_tsc",
    ("particles", "deposit_moments"): "true",
    ("particles", "couple_moments_to_mhd"): "true",
    ("particles", "couple_j_to_efield_representation"): "cell_centered",
    ("particles", "couple_j_deposition_mode"): "cc_convert",
    ("particles", "couple_moments_momentum_to_mhd"): "true",
    ("particles", "couple_moments_energy_to_mhd"): "true",
    ("particles", "couple_fluid_feedback_order"): "mhd_src_terms",
    ("particles", "pic_physical_mode"): "extended_mhd_pic",
    ("particles", "pic_background_mode"): "coupled",
    ("particles", "pic_feedback_mode"): "coupled",
    ("particles", "pic_cr_hall_mode"): "off",
    ("particles", "pic_wave_damping_mode"): "ion_neutral_friction",
    ("particles", "pic_ion_neutral_collision_rate"): "1.0e-4",
    ("particles", "pic_deltaf_mode"): "physical",
    ("particles", "pic_deltaf_f0"): "kappa_aniso",
    ("particles", "pic_deltaf_kappa"): "1.25",
    ("particles", "pic_deltaf_drift_x1"): "0.0",
    ("particles", "pic_deltaf_drift_x2"): "0.0",
    ("particles", "pic_deltaf_drift_x3"): "0.0",
    ("particles", "pic_deltaf_aniso_x1"): "1.0",
    ("particles", "pic_deltaf_aniso_x2"): "1.0",
    ("particles", "pic_deltaf_aniso_x3"): "1.0",
    ("particles", "pic_deltaf_background_jx"): "0.0",
    ("particles", "pic_deltaf_background_jy"): "0.0",
    ("particles", "pic_deltaf_background_jz"): "0.0",
    ("particles", "pic_deltaf_adapt_mode"): "global_bikappa_moments_experimental",
    ("particles", "pic_deltaf_adapt_interval"): "100.0",
    ("particles", "pic_expanding_box_mode"): "on",
    ("particles", "pic_expansion_law"): "exponential",
    ("particles", "track_displacement"): "true",
    ("problem", "pgen_name"): _PGEN_NAME,
    ("q033_transport_extension", "campaign_id"): CAMPAIGN_ID,
    ("q033_transport_extension", "deck_role"):
        "launch_blocked_source_local_transport_contract_only_not_evidence",
    ("q033_transport_extension", "qualification_effect"): "none",
    ("q033_transport_extension", "launch_status"): LAUNCH_STATUS,
    ("q033_transport_extension", "implementation_scope"):
        "existing_adaptive_deltaf_expanding_box_reduced_static_neutral_friction_only",
    ("q033_transport_extension", "physical_damping_scope"):
        "reduced_static_neutral_ion_neutral_friction_map_not_calibrated_astrophysical_damping",
    ("q033_transport_extension", "physical_calibration"): "open_not_claimed",
    ("q033_transport_extension", "qualification"): "open_not_claimed",
    ("q033_transport_extension", "mpi"): "open_not_claimed",
    ("q033_transport_extension", "frontier_hip"): "open_not_claimed",
}
_BUNDLE_KEYS = {
    "schema_version",
    "campaign_id",
    "artifact_role",
    "boundary_artifact_sha256",
    "q022_prerequisite_sha256",
    "series",
}
_SERIES_KEYS = {
    "series_id",
    "driving_rate",
    "nu_in",
    "reference_speed",
    "tail_samples",
}
_SAMPLE_KEYS = {"time", "xi", "parallel_displacement_variance", "spectrum"}
_SPECTRUM_KEYS = {
    "k",
    "forward_left",
    "forward_right",
    "backward_left",
    "backward_right",
}
_SPECTRUM_COMPONENTS = tuple(sorted(_SPECTRUM_KEYS - {"k"}))


class ContractError(ValueError):
    """Raised when a Q-033 source-local contract fails closed."""


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_athinput(path: Path) -> dict[str, dict[str, str]]:
    """Parse the strict Athena input subset used by the launch-blocked deck."""
    blocks: dict[str, dict[str, str]] = {}
    current = None
    for lineno, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<"):
            if not line.endswith(">"):
                raise ContractError(f"{path}:{lineno}: malformed block header")
            current = line[1:-1].strip()
            if not current:
                raise ContractError(f"{path}:{lineno}: empty block name")
            blocks.setdefault(current, {})
            continue
        if current is None or "=" not in line:
            raise ContractError(f"{path}:{lineno}: malformed parameter line")
        name, value = (item.strip() for item in line.split("=", 1))
        if not name or not value:
            raise ContractError(f"{path}:{lineno}: empty parameter name or value")
        if name in blocks[current]:
            raise ContractError(f"{path}:{lineno}: duplicate {current}/{name}")
        blocks[current][name] = value
    return blocks


def _require_exact_bindings(bindings: dict[str, str], label: str) -> dict[str, str]:
    for relative, expected in bindings.items():
        measured = _sha256(REPO_ROOT / relative)
        if measured != expected:
            raise ContractError(
                f"{relative}: {label} checksum mismatch: expected {expected}, "
                f"measured {measured}"
            )
    return dict(bindings)


def validate_source_bindings() -> dict[str, str]:
    """Require the exact implementation bytes from which the candidate was derived."""
    bindings = _require_exact_bindings(_SOURCE_BINDINGS, "source")
    for relative, snippets in _SOURCE_REQUIRED_SNIPPETS.items():
        contents = (REPO_ROOT / relative).read_text(encoding="utf-8")
        for snippet in snippets:
            if snippet not in contents:
                raise ContractError(f"{relative}: required implementation boundary is absent")
    return bindings


def validate_boundary_artifact_bindings() -> dict[str, str]:
    """Require exact prior bounded-local mechanics records without promoting them."""
    return _require_exact_bindings(_BOUNDARY_ARTIFACT_BINDINGS, "boundary artifact")


def _load_json(relative: str) -> dict[str, Any]:
    return json.loads((REPO_ROOT / relative).read_text(encoding="utf-8"))


def validate_q022_prerequisite_bindings() -> dict[str, str]:
    """Require exact fail-closed Q-022 CRPAI placeholders without promoting them."""
    bindings = _require_exact_bindings(_Q022_PREREQUISITE_BINDINGS, "Q-022 prerequisite")
    records = {Path(relative).name: _load_json(relative) for relative in bindings}
    ingest = records["q022_external_reference_private_ingest_2026-05-30.json"]
    references = [
        item for item in ingest.get("artifacts", [])
        if item.get("reference_id") == _Q022_REFERENCE_ID
    ]
    if len(references) != 1 or references[0].get("source_locator") != "arXiv:2409.08592":
        raise ContractError("Q-022 CRPAI Sun-Bai-Zhao source reference is not frozen")
    equation_map = records[
        "q022_xcmp_ext_crpai_transport_equation_map_2026-05-30.json"
    ]
    if (
        equation_map.get("comparison_id") != _Q022_COMPARISON_ID
        or equation_map.get("dataset_provenance_id") != _Q022_DATASET_ID
        or equation_map.get("reference_ids") != [_Q022_REFERENCE_ID]
        or equation_map.get("map_status")
        != "blocked_pending_reference_specific_mapping_and_external_review"
        or equation_map.get("matched_equations") != []
        or equation_map.get("unit_map") != {}
        or equation_map.get("normalization_map") != {}
        or equation_map.get("parameter_overlap") != {}
        or equation_map.get("reviewer_disposition") != "pending external review"
    ):
        raise ContractError("Q-022 CRPAI equation-map placeholder is not fail-closed")
    tolerance_table = records[
        "q022_xcmp_ext_crpai_transport_tolerance_table_2026-05-30.json"
    ]
    if (
        tolerance_table.get("comparison_id") != _Q022_COMPARISON_ID
        or tolerance_table.get("dataset_provenance_id") != _Q022_DATASET_ID
        or tolerance_table.get("freeze_status")
        != "blocked_pending_reference_dataset_extraction_and_external_review"
        or tolerance_table.get("rows") != []
        or tolerance_table.get("reviewer_disposition") != "pending external review"
    ):
        raise ContractError("Q-022 CRPAI tolerance-table placeholder is not fail-closed")
    provenance = records["q022_dataset_provenance_manifest_2026-05-30.json"]
    candidates = [
        item for item in provenance.get("dataset_candidates", [])
        if item.get("dataset_id") == _Q022_DATASET_ID
    ]
    if (
        len(candidates) != 1
        or candidates[0].get("comparison_id") != _Q022_COMPARISON_ID
        or candidates[0].get("reference_ids") != [_Q022_REFERENCE_ID]
        or candidates[0].get("extraction_status") != "blocked_extraction_input_unavailable"
        or candidates[0].get("reviewer_disposition") != "pending external review"
    ):
        raise ContractError("Q-022 CRPAI dataset-provenance placeholder is not fail-closed")
    return bindings


def validate_launch_block() -> str:
    """Require the Q-033 candidate generator to remain intentionally unavailable."""
    if _PGEN_NAME in _PGEN_DISPATCH.read_text(encoding="utf-8"):
        raise ContractError("Q-033 launch-block generator unexpectedly became available")
    return LAUNCH_STATUS


def validate_candidate_deck(path: Path = DECK) -> dict[str, Any]:
    """Validate the frozen deck while retaining its source-local launch block."""
    measured_sha256 = _sha256(path)
    if measured_sha256 != DECK_SHA256:
        raise ContractError(
            f"{path}: deck checksum mismatch: expected {DECK_SHA256}, "
            f"measured {measured_sha256}"
        )
    blocks = parse_athinput(path)
    for (block, name), expected in _EXPECTED_DECK_VALUES.items():
        measured = blocks.get(block, {}).get(name)
        if measured != expected:
            raise ContractError(
                f"{path}: {block}/{name}: expected {expected!r}, measured {measured!r}"
            )
    validate_launch_block()
    return {
        "path": str(path.relative_to(REPO_ROOT)),
        "sha256": measured_sha256,
        "launch_status": LAUNCH_STATUS,
        "qualification_effect": "none",
        "qualifying_evidence": False,
    }


def _finite_number(container: dict[str, Any], key: str) -> float:
    value = container[key]
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ContractError(f"{key} must be a finite number")
    measured = float(value)
    if not math.isfinite(measured):
        raise ContractError(f"{key} must be a finite number")
    return measured


def _require_positive(container: dict[str, Any], key: str) -> float:
    measured = _finite_number(container, key)
    if measured <= 0.0:
        raise ContractError(f"{key} must be positive")
    return measured


def _relative_span(values: list[float]) -> float:
    scale = max(abs(fmean(values)), 1.0e-30)
    return (max(values) - min(values)) / scale


def _spectrum_power(spectrum: list[dict[str, Any]]) -> float:
    if not isinstance(spectrum, list) or len(spectrum) != len(SPECTRUM_K_GRID):
        raise ContractError("Q-033 spectrum grid length mismatch")
    k_values = []
    totals = []
    for row in spectrum:
        if not isinstance(row, dict) or set(row) != _SPECTRUM_KEYS:
            raise ContractError("Q-033 spectrum keys do not match the contract")
        k_values.append(_require_positive(row, "k"))
        components = [_finite_number(row, key) for key in _SPECTRUM_COMPONENTS]
        if any(value < 0.0 for value in components):
            raise ContractError("Q-033 spectrum powers must be non-negative")
        totals.append(sum(components))
    if tuple(k_values) != SPECTRUM_K_GRID:
        raise ContractError("Q-033 spectrum k grid mismatch")
    return sum(
        0.5 * (totals[index] + totals[index + 1])
        * (k_values[index + 1] - k_values[index])
        for index in range(len(k_values) - 1)
    )


def build_synthetic_contract_bundle() -> dict[str, Any]:
    """Build normalized synthetic series used only to exercise the analyzer schema."""
    series = []
    branch_fractions = {
        "forward_left": 0.1,
        "forward_right": 0.4,
        "backward_left": 0.2,
        "backward_right": 0.3,
    }
    for driving_rate in DRIVING_RATE_GRID:
        for nu_in in NU_IN_GRID:
            nu_eff = math.sqrt(driving_rate * nu_in)
            diffusivity = 1.0 / (3.0 * nu_eff)
            xi = math.sqrt(1.0 + driving_rate / nu_eff)
            samples = []
            for time in TAIL_TIMES:
                spectrum = []
                for k_value in SPECTRUM_K_GRID:
                    amplitude = nu_eff if k_value == 2.0 else 0.0
                    spectrum.append({
                        "k": k_value,
                        **{
                            key: fraction * amplitude
                            for key, fraction in branch_fractions.items()
                        },
                    })
                samples.append({
                    "time": time,
                    "xi": xi,
                    "parallel_displacement_variance": 2.0 * diffusivity * time,
                    "spectrum": spectrum,
                })
            series.append({
                "series_id": f"drive_{driving_rate:g}_nu_in_{nu_in:g}",
                "driving_rate": driving_rate,
                "nu_in": nu_in,
                "reference_speed": 1.0,
                "tail_samples": samples,
            })
    return {
        "schema_version": 1,
        "campaign_id": CAMPAIGN_ID,
        "artifact_role": ARTIFACT_ROLE,
        "boundary_artifact_sha256": dict(_BOUNDARY_ARTIFACT_BINDINGS),
        "q022_prerequisite_sha256": dict(_Q022_PREREQUISITE_BINDINGS),
        "series": series,
    }


def _analyze_series(series: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(series, dict) or set(series) != _SERIES_KEYS:
        raise ContractError("Q-033 synthetic series keys do not match the contract")
    series_id = series["series_id"]
    if not isinstance(series_id, str) or not series_id:
        raise ContractError("Q-033 synthetic series_id must be a non-empty string")
    driving_rate = _require_positive(series, "driving_rate")
    nu_in = _require_positive(series, "nu_in")
    reference_speed = _require_positive(series, "reference_speed")
    samples = series["tail_samples"]
    if not isinstance(samples, list) or len(samples) != len(TAIL_TIMES):
        raise ContractError("Q-033 tail sample count mismatch")

    times = []
    xi_values = []
    variances = []
    spectrum_powers = []
    for sample in samples:
        if not isinstance(sample, dict) or set(sample) != _SAMPLE_KEYS:
            raise ContractError("Q-033 tail sample keys do not match the contract")
        times.append(_finite_number(sample, "time"))
        xi_values.append(_require_positive(sample, "xi"))
        variances.append(_finite_number(sample, "parallel_displacement_variance"))
        spectrum_powers.append(_spectrum_power(sample["spectrum"]))
    if tuple(times) != TAIL_TIMES:
        raise ContractError("Q-033 tail time grid mismatch")
    if any(value < 0.0 for value in variances):
        raise ContractError("Q-033 parallel displacement variances must be non-negative")

    local_nu_eff = []
    for index in range(len(times) - 1):
        dt = times[index + 1] - times[index]
        delta_variance = variances[index + 1] - variances[index]
        if dt <= 0.0 or delta_variance <= 0.0:
            raise ContractError("Q-033 displacement variance slope must be positive")
        diffusivity = delta_variance / (2.0 * dt)
        local_nu_eff.append(reference_speed * reference_speed / (3.0 * diffusivity))
    anisotropies = [abs(xi * xi - 1.0) for xi in xi_values]
    nu_eff = fmean(local_nu_eff)
    anisotropy = fmean(anisotropies)
    spectrum_power = fmean(spectrum_powers)
    expected_nu_eff = math.sqrt(driving_rate * nu_in)
    expected_anisotropy = driving_rate / expected_nu_eff
    expected_spectrum_power = expected_nu_eff
    spans = {
        "nu_eff": _relative_span(local_nu_eff),
        "anisotropy": _relative_span(anisotropies),
        "spectrum_power": _relative_span(spectrum_powers),
    }
    quasi_steady = all(
        value <= QUASI_STEADY_RELATIVE_SPAN_MAX for value in spans.values()
    )
    consistent = (
        math.isclose(nu_eff, expected_nu_eff, rel_tol=0.0, abs_tol=ABSOLUTE_TOLERANCE)
        and math.isclose(
            anisotropy, expected_anisotropy, rel_tol=0.0, abs_tol=ABSOLUTE_TOLERANCE
        )
        and math.isclose(
            spectrum_power,
            expected_spectrum_power,
            rel_tol=0.0,
            abs_tol=ABSOLUTE_TOLERANCE,
        )
        and quasi_steady
    )
    return {
        "series_id": series_id,
        "driving_rate": driving_rate,
        "nu_in": nu_in,
        "nu_eff": nu_eff,
        "anisotropy": anisotropy,
        "polarization_resolved_spectrum_power": spectrum_power,
        "quasi_steady_relative_spans": spans,
        "synthetic_quasi_steady": quasi_steady,
        "synthetic_series_consistent": consistent,
    }


def _log_slope(low_value: float, high_value: float, low_x: float, high_x: float) -> float:
    return math.log(high_value / low_value) / math.log(high_x / low_x)


def _analyze_scaling(reports: list[dict[str, Any]]) -> dict[str, Any]:
    indexed = {(report["driving_rate"], report["nu_in"]): report for report in reports}
    driving_slopes = [
        _log_slope(
            indexed[(DRIVING_RATE_GRID[0], nu_in)]["nu_eff"],
            indexed[(DRIVING_RATE_GRID[1], nu_in)]["nu_eff"],
            DRIVING_RATE_GRID[0],
            DRIVING_RATE_GRID[1],
        )
        for nu_in in NU_IN_GRID
    ]
    friction_slopes = [
        _log_slope(
            indexed[(driving_rate, NU_IN_GRID[0])]["nu_eff"],
            indexed[(driving_rate, NU_IN_GRID[1])]["nu_eff"],
            NU_IN_GRID[0],
            NU_IN_GRID[1],
        )
        for driving_rate in DRIVING_RATE_GRID
    ]
    consistent = all(
        math.isclose(value, 0.5, rel_tol=0.0, abs_tol=ABSOLUTE_TOLERANCE)
        for value in driving_slopes + friction_slopes
    )
    return {
        "synthetic_oracle_only_not_physical_scaling_law": True,
        "nu_eff_vs_driving_rate_log_slopes": driving_slopes,
        "nu_eff_vs_nu_in_log_slopes": friction_slopes,
        "expected_synthetic_log_slope": 0.5,
        "synthetic_scaling_consistent": consistent,
    }


def analyze_synthetic_contract_bundle(bundle: dict[str, Any]) -> dict[str, Any]:
    """Recompute synthetic observables while always returning nonqualifying status."""
    source_bindings = validate_source_bindings()
    boundary_bindings = validate_boundary_artifact_bindings()
    q022_prerequisite_bindings = validate_q022_prerequisite_bindings()
    deck = validate_candidate_deck()
    if not isinstance(bundle, dict) or set(bundle) != _BUNDLE_KEYS:
        raise ContractError("Q-033 synthetic bundle keys do not match the contract")
    if bundle["schema_version"] != 1 or bundle["campaign_id"] != CAMPAIGN_ID:
        raise ContractError("Q-033 synthetic bundle identity mismatch")
    if bundle["artifact_role"] != ARTIFACT_ROLE:
        raise ContractError("Q-033 bundle is not synthetic transport-contract input")
    if bundle["boundary_artifact_sha256"] != boundary_bindings:
        raise ContractError("Q-033 synthetic bundle boundary artifact checksums mismatch")
    if bundle["q022_prerequisite_sha256"] != q022_prerequisite_bindings:
        raise ContractError("Q-033 synthetic bundle Q-022 prerequisite checksums mismatch")
    if not isinstance(bundle["series"], list):
        raise ContractError("Q-033 synthetic series must be a list")

    expected_grid = {
        (driving_rate, nu_in)
        for driving_rate in DRIVING_RATE_GRID
        for nu_in in NU_IN_GRID
    }
    measured_grid = []
    reports = []
    for series in bundle["series"]:
        report = _analyze_series(series)
        key = (report["driving_rate"], report["nu_in"])
        if key in measured_grid:
            raise ContractError(f"duplicate Q-033 synthetic series {key!r}")
        measured_grid.append(key)
        reports.append(report)
    if set(measured_grid) != expected_grid:
        raise ContractError("Q-033 synthetic series grid is incomplete or contains extras")
    reports.sort(key=lambda item: (item["driving_rate"], item["nu_in"]))
    scaling = _analyze_scaling(reports)
    consistent = (
        all(report["synthetic_series_consistent"] for report in reports)
        and scaling["synthetic_scaling_consistent"]
    )
    return {
        "schema_version": 1,
        "campaign_id": CAMPAIGN_ID,
        "artifact_role": ARTIFACT_ROLE,
        "qualification_effect": QUALIFICATION_EFFECT,
        "qualifying_evidence": False,
        "physical_calibration_claimed": False,
        "launch_status": LAUNCH_STATUS,
        "source_bindings": source_bindings,
        "boundary_artifact_bindings": boundary_bindings,
        "q022_prerequisite_bindings": q022_prerequisite_bindings,
        "q022_prerequisite_status":
            "blocked_exact_placeholders_bound_nonqualifying",
        "deck_contract": deck,
        "series_count": len(reports),
        "series": reports,
        "scaling": scaling,
        "synthetic_transport_contract_consistent": consistent,
        "status": (
            "source_local_synthetic_contract_consistent_not_physical_calibration_"
            "not_qualifying_evidence"
            if consistent
            else "source_local_synthetic_contract_mismatch_not_physical_calibration_"
            "not_qualifying_evidence"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("synthetic_bundle", type=Path)
    args = parser.parse_args()
    bundle = json.loads(args.synthetic_bundle.read_text(encoding="utf-8"))
    print(json.dumps(analyze_synthetic_contract_bundle(bundle), indent=2,
                     sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
