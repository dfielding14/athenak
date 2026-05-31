#!/usr/bin/env python3
"""Q-032 reduced static-neutral friction synthetic-contract analyzer."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN_ID = "Q032-EXT-CRSI-IN-DAMPING-CANDIDATE"
COMPARISON_ID = "XCMP-EXT-CRSI-IN-DAMPING"
CLAIM_ID = "CLAIM-EXT-CRSI-IN-DAMPING-001"
REFERENCE_ID = "plotnikov_ostriker_bai_2021_arxiv_2102.11878"
DATASET_ID = "Q022-DATASET-XCMP-EXT-CRSI-IN-DAMPING"
ARTIFACT_ROLE = "source_local_synthetic_static_neutral_contract_only"
QUALIFICATION_EFFECT = "none_source_local_launch_blocked_contract_only"
LAUNCH_STATUS = "blocked_missing_plotnikov_inputs_review_authorization_and_generator"
DECK = REPO_ROOT / "inputs/tests/pic_q032_plotnikov_damped_crsi_candidate.athinput"
DECK_SHA256 = "ae471a610d8c2c50367c58df372361c1dfeec08644224789f9cba8a55df722a9"
_PGEN_DISPATCH = REPO_ROOT / "src/pgen/pgen.cpp"
_PGEN_NAME = "q032_plotnikov_damped_crsi_open"
COLLISION_RATE_GRID = (0.0, 0.7, 1.4)
DT_GRID = (0.0, 0.125, 0.5)
ABSOLUTE_TOLERANCE = 1.0e-14

_SOURCE_BINDINGS = {
    "src/mhd/mhd_tasks.cpp":
        "339326c2cc8913c321affeafc618a8fea6cbd9133ea6b26a9b58daa2baafa6d2",
    "src/particles/particles.hpp":
        "452e6091ad6f1ce50156c770d304ac36020cb53927cc1e40ba6af73e7c3bb462",
    "src/particles/particles.cpp":
        "06f80ec9a93f027ff5a24a2631754603c7ea3f83cfbcdbd8e5b9a32be711fa79",
    "docs/source/engineering/pic_mhd_model_contract.md":
        "ebad2e21725ea15eb301e03f32d8ca8fab9c6503aec90cd7090be4637402ed64",
}
_SOURCE_REQUIRED_SNIPPETS = {
    "src/mhd/mhd_tasks.cpp": (
        "const Real factor = exp(-ppart->pic_ion_neutral_collision_rate*",
        "u(m, IM2, k, j, i) *= factor;",
        "u(m, IM3, k, j, i) *= factor;",
        "u(m, IEN, k, j, i) -= transverse_ekin_loss;",
    ),
    "src/particles/particles.hpp": (
        "Real pic_ion_neutral_collision_rate = 0.0;",
        "return pic_wave_damping_mode == PICWaveDampingMode::ion_neutral_friction;",
    ),
    "src/particles/particles.cpp": (
        'pic_wave_damping_mode_str.compare("ion_neutral_friction") == 0',
        '"an active coupled <mhd> background"',
        '"the Newtonian single-fluid MHD source task path; does not support "',
        '"<ion-neutral> alternate task lists"',
    ),
    "docs/source/engineering/pic_mhd_model_contract.md": (
        "p_ion,perp(t + dt) = p_ion,perp(t) exp(-nu_in dt)",
        "static neutrals and does not establish a CRSI dispersion comparison",
    ),
}
_PREREQUISITE_BINDINGS = {
    "tst/publication/readiness/q022_external_reference_private_ingest_2026-05-30.json":
        "78287ea54d8350445cfad62e72f6efc2845f720580d070d075897d889ae3a214",
    "tst/publication/readiness/q022_dataset_provenance_manifest_2026-05-30.json":
        "dd6dd6e3187b8e5e06bfbce6f0da4cfc0623d79143a55d7ebbcc797b21b38c78",
    "tst/publication/readiness/q022_xcmp_ext_crsi_in_damping_equation_map_2026-05-30.json":
        "87ed9f5fac74cf9ce7be567a1e57cfecae6e2b3374e1bc4e4ed74d2bbe3c8c9b",
    "tst/publication/readiness/q022_xcmp_ext_crsi_in_damping_tolerance_table_2026-05-30.json":
        "3ac72c052fc465bd54c92750a96e3ee6ee1d39a79eec2f9c52eda3eb5cc5dfd7",
    "tst/publication/readiness/q032_q033_extension_local_scans_2026-05-30.json":
        "8338ee08f311ce353186a514549cec1227d17de300557645655158ee73da3dc4",
}
_EXPECTED_DECK_VALUES = {
    ("time", "nlim"): "0",
    ("time", "tlim"): "0.0",
    ("particles", "particle_type"): "cosmic_ray",
    ("particles", "pusher"): "boris_tsc",
    ("particles", "deposit_moments"): "true",
    ("particles", "couple_moments_to_mhd"): "true",
    ("particles", "pic_physical_mode"): "extended_mhd_pic",
    ("particles", "pic_background_mode"): "coupled",
    ("particles", "pic_feedback_mode"): "coupled",
    ("particles", "pic_cr_hall_mode"): "off",
    ("particles", "pic_wave_damping_mode"): "ion_neutral_friction",
    ("particles", "pic_ion_neutral_collision_rate"): "0.7",
    ("particles", "pic_expanding_box_mode"): "off",
    ("problem", "pgen_name"): _PGEN_NAME,
    ("q032_damped_crsi_extension", "campaign_id"): CAMPAIGN_ID,
    ("q032_damped_crsi_extension", "comparison_id"): COMPARISON_ID,
    ("q032_damped_crsi_extension", "claim_id"): CLAIM_ID,
    ("q032_damped_crsi_extension", "deck_role"):
        "launch_blocked_source_local_reduced_map_only_not_plotnikov_evidence",
    ("q032_damped_crsi_extension", "qualification_effect"): "none",
    ("q032_damped_crsi_extension", "launch_status"): LAUNCH_STATUS,
    ("q032_damped_crsi_extension", "implementation_scope"):
        "reduced_static_neutral_high_frequency_transverse_friction_only",
    ("q032_damped_crsi_extension", "nu_in_parameter"):
        "particles/pic_ion_neutral_collision_rate",
    ("q032_damped_crsi_extension", "nu_in_fiducial"): "0.7",
    ("q032_damped_crsi_extension", "attenuation_factor"): "exp(-nu_in*dt)",
    ("q032_damped_crsi_extension", "collision_rate_candidate_grid"):
        "0.0,0.2,0.7,1.4",
    ("q032_damped_crsi_extension", "synthetic_collision_rate_grid"):
        "0.0,0.7,1.4",
    ("q032_damped_crsi_extension", "synthetic_dt_grid"): "0.0,0.125,0.5",
    ("q032_damped_crsi_extension", "plotnikov_reference_id"): REFERENCE_ID,
    ("q032_damped_crsi_extension", "plotnikov_map"):
        "open_reference_specific_mapping_and_review",
    ("q032_damped_crsi_extension", "plotnikov_dataset"):
        "open_authorized_orion_extraction_with_uncertainty",
    ("q032_damped_crsi_extension", "plotnikov_tolerances"):
        "open_externally_reviewed_numeric_rows",
    ("q032_damped_crsi_extension", "matched_damped_crsi_qualification"):
        "open_not_claimed",
}
_BUNDLE_KEYS = {"schema_version", "campaign_id", "artifact_role", "samples"}
_SAMPLE_KEYS = {
    "nu_in",
    "dt",
    "density_before",
    "density_after",
    "mom1_before",
    "mom1_after",
    "mom2_before",
    "mom2_after",
    "mom3_before",
    "mom3_after",
    "energy_before",
    "energy_after",
}


class ContractError(ValueError):
    """Raised when a Q-032 source-local contract fails closed."""


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_json(relative: str) -> dict[str, Any]:
    return json.loads((REPO_ROOT / relative).read_text(encoding="utf-8"))


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


def _csv_floats(value: str) -> tuple[float, ...]:
    return tuple(float(item) for item in value.split(","))


def validate_source_bindings() -> dict[str, str]:
    """Require the exact implementation bytes from which this tranche was derived."""
    for relative, expected in _SOURCE_BINDINGS.items():
        path = REPO_ROOT / relative
        measured = _sha256(path)
        if measured != expected:
            raise ContractError(
                f"{relative}: source checksum mismatch: expected {expected}, "
                f"measured {measured}"
            )
        contents = path.read_text(encoding="utf-8")
        for snippet in _SOURCE_REQUIRED_SNIPPETS.get(relative, ()):
            if snippet not in contents:
                raise ContractError(f"{relative}: required source boundary is absent")
    return dict(_SOURCE_BINDINGS)


def validate_prerequisite_boundaries() -> dict[str, str]:
    """Require Plotnikov comparison inputs to remain explicitly blocked and empty."""
    for relative, expected in _PREREQUISITE_BINDINGS.items():
        measured = _sha256(REPO_ROOT / relative)
        if measured != expected:
            raise ContractError(
                f"{relative}: prerequisite checksum mismatch: expected {expected}, "
                f"measured {measured}"
            )

    ingest = _load_json(
        "tst/publication/readiness/q022_external_reference_private_ingest_2026-05-30.json"
    )
    reference = next(
        (item for item in ingest["artifacts"] if item["reference_id"] == REFERENCE_ID),
        None,
    )
    if reference is None or reference["source_locator"] != "arXiv:2102.11878":
        raise ContractError("Q-032 Plotnikov source reference is not frozen")

    provenance = _load_json(
        "tst/publication/readiness/q022_dataset_provenance_manifest_2026-05-30.json"
    )
    dataset = next(
        (item for item in provenance["dataset_candidates"]
         if item["dataset_id"] == DATASET_ID),
        None,
    )
    if (
        dataset is None
        or dataset["comparison_id"] != COMPARISON_ID
        or dataset["reference_ids"] != [REFERENCE_ID]
        or dataset["extraction_status"] != "blocked_extraction_input_unavailable"
        or dataset["reviewer_disposition"] != "pending external review"
    ):
        raise ContractError("Q-032 Plotnikov dataset boundary is not blocked")

    equation_map = _load_json(
        "tst/publication/readiness/"
        "q022_xcmp_ext_crsi_in_damping_equation_map_2026-05-30.json"
    )
    if (
        equation_map["comparison_id"] != COMPARISON_ID
        or equation_map["reference_ids"] != [REFERENCE_ID]
        or equation_map["map_status"]
        != "blocked_pending_reference_specific_mapping_and_external_review"
        or equation_map["matched_equations"]
        or equation_map["unit_map"]
        or equation_map["normalization_map"]
        or equation_map["parameter_overlap"]
        or equation_map["reviewer_disposition"] != "pending external review"
    ):
        raise ContractError("Q-032 Plotnikov equation map is not fail-closed")

    tolerances = _load_json(
        "tst/publication/readiness/"
        "q022_xcmp_ext_crsi_in_damping_tolerance_table_2026-05-30.json"
    )
    if (
        tolerances["comparison_id"] != COMPARISON_ID
        or tolerances["dataset_provenance_id"] != DATASET_ID
        or tolerances["freeze_status"]
        != "blocked_pending_reference_dataset_extraction_and_external_review"
        or tolerances["rows"]
        or tolerances["reviewer_disposition"] != "pending external review"
    ):
        raise ContractError("Q-032 Plotnikov tolerances are not fail-closed")

    scans = _load_json(
        "tst/publication/readiness/q032_q033_extension_local_scans_2026-05-30.json"
    )
    q032_scan = scans["q032_ion_neutral_alfven_envelope"]
    if (
        q032_scan["local_result"] != "pass"
        or q032_scan["status"] != "local_scan_pass_full_question_open"
        or not q032_scan["limitations"]
    ):
        raise ContractError("Q-032 local scan boundary is not explicitly limited")
    return dict(_PREREQUISITE_BINDINGS)


def validate_launch_block() -> str:
    """Require the dedicated Q-032 candidate generator to remain unavailable."""
    if _PGEN_NAME in _PGEN_DISPATCH.read_text(encoding="utf-8"):
        raise ContractError("Q-032 launch-block generator unexpectedly became available")
    return LAUNCH_STATUS


def validate_candidate_deck(path: Path = DECK) -> dict[str, Any]:
    """Validate the frozen reduced-map deck while retaining its launch block."""
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
    metadata = blocks["q032_damped_crsi_extension"]
    if _csv_floats(metadata["synthetic_collision_rate_grid"]) != COLLISION_RATE_GRID:
        raise ContractError("Q-032 synthetic collision-rate grid mismatch")
    if _csv_floats(metadata["synthetic_dt_grid"]) != DT_GRID:
        raise ContractError("Q-032 synthetic dt grid mismatch")
    if float(metadata["nu_in_fiducial"]) <= 0.0:
        raise ContractError("Q-032 fiducial collision rate must be positive")
    validate_launch_block()
    return {
        "path": str(path.relative_to(REPO_ROOT)),
        "sha256": measured_sha256,
        "launch_status": LAUNCH_STATUS,
        "qualification_effect": "none",
        "qualifying_evidence": False,
        "matched_plotnikov_qualification": False,
        "implementation_scope": metadata["implementation_scope"],
        "nu_in_fiducial": float(metadata["nu_in_fiducial"]),
        "collision_rate_candidate_grid":
            list(_csv_floats(metadata["collision_rate_candidate_grid"])),
    }


def _finite_number(sample: dict[str, Any], key: str) -> float:
    value = sample[key]
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ContractError(f"{key} must be a finite number")
    measured = float(value)
    if not math.isfinite(measured):
        raise ContractError(f"{key} must be a finite number")
    return measured


def _expected_state(nu_in: float, dt: float) -> dict[str, float]:
    density = 2.0
    mom1 = 0.25
    mom2 = 0.75
    mom3 = -0.5
    energy = 4.0
    factor = math.exp(-nu_in * dt)
    return {
        "density_before": density,
        "density_after": density,
        "mom1_before": mom1,
        "mom1_after": mom1,
        "mom2_before": mom2,
        "mom2_after": factor * mom2,
        "mom3_before": mom3,
        "mom3_after": factor * mom3,
        "energy_before": energy,
        "energy_after": energy - 0.5 * (1.0 - factor * factor)
        * (mom2 * mom2 + mom3 * mom3) / density,
    }


def build_synthetic_contract_bundle() -> dict[str, Any]:
    """Build the fixed synthetic endpoint-contract grid used by focused tests."""
    samples = []
    for nu_in in COLLISION_RATE_GRID:
        for dt in DT_GRID:
            samples.append({"nu_in": nu_in, "dt": dt, **_expected_state(nu_in, dt)})
    return {
        "schema_version": 1,
        "campaign_id": CAMPAIGN_ID,
        "artifact_role": ARTIFACT_ROLE,
        "samples": samples,
    }


def analyze_synthetic_contract_bundle(bundle: dict[str, Any]) -> dict[str, Any]:
    """Analyze only synthetic algebra and never emit Plotnikov qualification."""
    source_bindings = validate_source_bindings()
    prerequisite_bindings = validate_prerequisite_boundaries()
    deck = validate_candidate_deck()
    if set(bundle) != _BUNDLE_KEYS:
        raise ContractError("Q-032 synthetic bundle keys do not match the contract")
    if bundle["schema_version"] != 1 or bundle["campaign_id"] != CAMPAIGN_ID:
        raise ContractError("Q-032 synthetic bundle identity mismatch")
    if bundle["artifact_role"] != ARTIFACT_ROLE:
        raise ContractError("Q-032 bundle is not synthetic endpoint-contract input")
    if not isinstance(bundle["samples"], list):
        raise ContractError("Q-032 synthetic samples must be a list")

    expected_keys = {(nu_in, dt) for nu_in in COLLISION_RATE_GRID for dt in DT_GRID}
    measured_keys = []
    reports = []
    state_keys = _SAMPLE_KEYS - {"nu_in", "dt"}
    for sample in bundle["samples"]:
        if not isinstance(sample, dict) or set(sample) != _SAMPLE_KEYS:
            raise ContractError("Q-032 synthetic sample keys do not match the contract")
        measured = {key: _finite_number(sample, key) for key in _SAMPLE_KEYS}
        key = (measured["nu_in"], measured["dt"])
        if key in measured_keys:
            raise ContractError(f"duplicate Q-032 synthetic sample {key!r}")
        measured_keys.append(key)
        expected = _expected_state(*key)
        errors = {name: abs(measured[name] - expected[name]) for name in state_keys}
        reports.append({
            "nu_in": key[0],
            "dt": key[1],
            "factor": math.exp(-key[0] * key[1]),
            "maximum_absolute_error": max(errors.values()),
            "synthetic_contract_consistent":
                max(errors.values()) <= ABSOLUTE_TOLERANCE,
        })
    if set(measured_keys) != expected_keys:
        raise ContractError("Q-032 synthetic grid is incomplete or contains extra samples")

    reports.sort(key=lambda item: (item["nu_in"], item["dt"]))
    consistent = all(item["synthetic_contract_consistent"] for item in reports)
    return {
        "schema_version": 1,
        "campaign_id": CAMPAIGN_ID,
        "comparison_id": COMPARISON_ID,
        "claim_id": CLAIM_ID,
        "artifact_role": ARTIFACT_ROLE,
        "qualification_effect": QUALIFICATION_EFFECT,
        "qualifying_evidence": False,
        "matched_plotnikov_qualification": False,
        "launch_status": LAUNCH_STATUS,
        "source_bindings": source_bindings,
        "prerequisite_bindings": prerequisite_bindings,
        "deck_contract": deck,
        "sample_count": len(reports),
        "samples": reports,
        "synthetic_static_neutral_contract_consistent": consistent,
        "status": (
            "source_local_contract_consistent_not_plotnikov_qualification"
            if consistent
            else "source_local_contract_mismatch_not_plotnikov_qualification"
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
