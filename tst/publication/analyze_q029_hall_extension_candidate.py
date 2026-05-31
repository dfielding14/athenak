#!/usr/bin/env python3
"""Q-029 source-local experimental Hall normalization contract analyzer."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN_ID = "Q029-EXT-HALL-NORMALIZATION-CANDIDATE"
ARTIFACT_ROLE = "source_local_synthetic_normalization_contract_only"
QUALIFICATION_EFFECT = "none_source_local_launch_blocked_contract_only"
LAUNCH_STATUS = (
    "blocked_missing_q029_generator_and_authorization"
)
DECK = REPO_ROOT / "inputs/tests/pic_q029_extended_hall_normalization_candidate.athinput"
DECK_SHA256 = "dc159e5f9a80aec8bb94723b69f5de9edaafd624889ac7ea9a88518b46e6bb80"
_PGEN_DISPATCH = REPO_ROOT / "src/pgen/pgen.cpp"
_PGEN_NAME = "q029_extended_hall_normalization_open"
ALPHA_H_SOURCE_ORACLE_GRID = (-1.0, 0.0, 1.0)
PROJECTED_CURRENT_GRID = (-0.75, 0.25, 1.5)
ABSOLUTE_TOLERANCE = 1.0e-14

_SOURCE_BINDINGS = {
    "src/mhd/mhd_tasks.cpp":
        "339326c2cc8913c321affeafc618a8fea6cbd9133ea6b26a9b58daa2baafa6d2",
    "src/particles/particles.hpp":
        "452e6091ad6f1ce50156c770d304ac36020cb53927cc1e40ba6af73e7c3bb462",
    "src/particles/particles.cpp":
        "83bc366df553a70fc1470f9ed461dc98e8d7506fb9327d4b5d5abb97484f1637",
    "docs/source/engineering/pic_mhd_model_contract.md":
        "ebad2e21725ea15eb301e03f32d8ca8fab9c6503aec90cd7090be4637402ed64",
}
_SOURCE_REQUIRED_SNIPPETS = {
    "src/mhd/mhd_tasks.cpp": (
        "const Real jcoef = ppart->couple_j_to_efield_coeff;",
        "e2(m,ks  ,js  ,i) += jcoef*mom(m, particles::Particles::IMOM_JY, ks, js, i);",
    ),
    "src/particles/particles.hpp": (
        "bool AddsCRCurrentToCT() const {",
        "(pic_cr_hall_mode == PICCRHallMode::current_to_ct_experimental)",
    ),
    "docs/source/engineering/pic_mhd_model_contract.md": (
        "cE_CT = cE_ideal + alpha_H P_edge[J_CR]",
        "alpha_H = <particles>/couple_j_to_efield_coeff",
    ),
}
_EXPECTED_DECK_VALUES = {
    ("time", "nlim"): "0",
    ("time", "tlim"): "0.0",
    ("particles", "particle_type"): "cosmic_ray",
    ("particles", "pusher"): "boris_tsc",
    ("particles", "deposit_moments"): "true",
    ("particles", "couple_moments_to_mhd"): "true",
    ("particles", "couple_j_to_efield_coeff"): "0.5",
    ("particles", "couple_j_to_efield_representation"): "cell_centered",
    ("particles", "couple_j_deposition_mode"): "cc_convert",
    ("particles", "pic_physical_mode"): "extended_mhd_pic",
    ("particles", "pic_background_mode"): "coupled",
    ("particles", "pic_feedback_mode"): "coupled",
    ("particles", "pic_cr_hall_mode"): "current_to_ct_experimental",
    ("particles", "pic_wave_damping_mode"): "off",
    ("problem", "pgen_name"): _PGEN_NAME,
    ("q029_hall_extension", "campaign_id"): CAMPAIGN_ID,
    ("q029_hall_extension", "deck_role"):
        "launch_blocked_source_local_normalization_only_not_evidence",
    ("q029_hall_extension", "qualification_effect"): "none",
    ("q029_hall_extension", "launch_status"): LAUNCH_STATUS,
    ("q029_hall_extension", "normalization_scope"):
        "implemented_source_only_not_physical_hall_bell_map",
    ("q029_hall_extension", "alpha_h_parameter"):
        "particles/couple_j_to_efield_coeff",
    ("q029_hall_extension", "alpha_h_fiducial"): "0.5",
    ("q029_hall_extension", "chi_h_definition"): "alpha_h*j_cr0/(u_a0*b_g0)",
    ("q029_hall_extension", "chi_h_fiducial"): "0.5",
    ("q029_hall_extension", "chi_h_candidate_grid"): "0.0,0.25,0.5,1.0",
    ("q029_hall_extension", "signed_source_oracle_grid"): "-1.0,0.0,1.0",
    ("q029_hall_extension", "rho0"): "1.0",
    ("q029_hall_extension", "b_g0"): "1.0",
    ("q029_hall_extension", "u_a0"): "1.0",
    ("q029_hall_extension", "j_cr0"): "1.0",
    ("q029_hall_extension", "c_e0"): "1.0",
    ("q029_hall_extension", "literature_map"):
        "open_external_review_bai2015_mapping",
    ("q029_hall_extension", "linear_bell"): "open_not_claimed",
    ("q029_hall_extension", "nonlinear_bell"): "open_not_claimed",
    ("q029_hall_extension", "shock_front"): "open_not_claimed",
    ("q029_hall_extension", "gpu"): "open_not_claimed",
    ("q029_hall_extension", "mpi"): "open_not_claimed",
}
_BUNDLE_KEYS = {"schema_version", "campaign_id", "artifact_role", "samples"}
_SAMPLE_KEYS = {
    "alpha_h",
    "projected_j_cr_over_j_cr0",
    "ideal_c_e_over_c_e0",
    "observed_c_e_over_c_e0",
}


class ContractError(ValueError):
    """Raised when a Q-029 source-local contract fails closed."""


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


def _require_close(label: str, measured: float, expected: float) -> None:
    if not math.isclose(measured, expected, rel_tol=0.0, abs_tol=ABSOLUTE_TOLERANCE):
        raise ContractError(f"{label}: expected {expected!r}, measured {measured!r}")


def _csv_floats(value: str) -> tuple[float, ...]:
    return tuple(float(item) for item in value.split(","))


def validate_source_bindings() -> dict[str, str]:
    """Require the exact implementation bytes from which the candidate was derived."""
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
                raise ContractError(
                    f"{relative}: required source normalization is absent"
                )
    return dict(_SOURCE_BINDINGS)


def validate_launch_block() -> str:
    """Require the Q-029 candidate generator to remain intentionally unavailable."""
    if _PGEN_NAME in _PGEN_DISPATCH.read_text(encoding="utf-8"):
        raise ContractError("Q-029 launch-block generator unexpectedly became available")
    return LAUNCH_STATUS


def validate_candidate_deck(path: Path = DECK) -> dict[str, Any]:
    """Validate the frozen normalization deck while retaining its launch block."""
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

    metadata = blocks["q029_hall_extension"]
    rho0 = float(metadata["rho0"])
    b_g0 = float(metadata["b_g0"])
    u_a0 = float(metadata["u_a0"])
    j_cr0 = float(metadata["j_cr0"])
    c_e0 = float(metadata["c_e0"])
    alpha_h = float(metadata["alpha_h_fiducial"])
    chi_h = float(metadata["chi_h_fiducial"])
    if min(rho0, b_g0, u_a0, j_cr0, c_e0) <= 0.0:
        raise ContractError("Q-029 normalization scales must be positive")
    _require_close("U_A0", u_a0, b_g0 / math.sqrt(rho0))
    _require_close("cE0", c_e0, u_a0 * b_g0)
    _require_close("chi_H", chi_h, alpha_h * j_cr0 / c_e0)
    if _csv_floats(metadata["signed_source_oracle_grid"]) != ALPHA_H_SOURCE_ORACLE_GRID:
        raise ContractError("Q-029 signed source-oracle grid mismatch")
    validate_launch_block()

    return {
        "path": str(path.relative_to(REPO_ROOT)),
        "sha256": measured_sha256,
        "launch_status": LAUNCH_STATUS,
        "qualification_effect": "none",
        "qualifying_evidence": False,
        "normalization": {
            "rho0": rho0,
            "b_g0": b_g0,
            "u_a0": u_a0,
            "j_cr0": j_cr0,
            "c_e0": c_e0,
            "alpha_h_fiducial": alpha_h,
            "chi_h_fiducial": chi_h,
            "chi_h_candidate_grid": list(_csv_floats(metadata["chi_h_candidate_grid"])),
        },
    }


def _finite_number(sample: dict[str, Any], key: str) -> float:
    value = sample[key]
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ContractError(f"{key} must be a finite number")
    measured = float(value)
    if not math.isfinite(measured):
        raise ContractError(f"{key} must be a finite number")
    return measured


def build_synthetic_contract_bundle() -> dict[str, Any]:
    """Build the fixed synthetic source-contract grid used by focused tests."""
    samples = []
    for alpha_h in ALPHA_H_SOURCE_ORACLE_GRID:
        for projected_j in PROJECTED_CURRENT_GRID:
            ideal_c_e = 0.125
            samples.append({
                "alpha_h": alpha_h,
                "projected_j_cr_over_j_cr0": projected_j,
                "ideal_c_e_over_c_e0": ideal_c_e,
                "observed_c_e_over_c_e0": ideal_c_e + alpha_h * projected_j,
            })
    return {
        "schema_version": 1,
        "campaign_id": CAMPAIGN_ID,
        "artifact_role": ARTIFACT_ROLE,
        "samples": samples,
    }


def analyze_synthetic_contract_bundle(bundle: dict[str, Any]) -> dict[str, Any]:
    """Analyze only a synthetic contract bundle and never emit qualifying evidence."""
    source_bindings = validate_source_bindings()
    deck = validate_candidate_deck()
    if set(bundle) != _BUNDLE_KEYS:
        raise ContractError("Q-029 synthetic bundle keys do not match the contract")
    if bundle["schema_version"] != 1 or bundle["campaign_id"] != CAMPAIGN_ID:
        raise ContractError("Q-029 synthetic bundle identity mismatch")
    if bundle["artifact_role"] != ARTIFACT_ROLE:
        raise ContractError("Q-029 bundle is not synthetic normalization-contract input")
    if not isinstance(bundle["samples"], list):
        raise ContractError("Q-029 synthetic samples must be a list")

    expected_keys = {
        (alpha_h, projected_j)
        for alpha_h in ALPHA_H_SOURCE_ORACLE_GRID
        for projected_j in PROJECTED_CURRENT_GRID
    }
    measured_keys = []
    reports = []
    j_cr0 = deck["normalization"]["j_cr0"]
    c_e0 = deck["normalization"]["c_e0"]
    for sample in bundle["samples"]:
        if not isinstance(sample, dict) or set(sample) != _SAMPLE_KEYS:
            raise ContractError("Q-029 synthetic sample keys do not match the contract")
        alpha_h = _finite_number(sample, "alpha_h")
        projected_j = _finite_number(sample, "projected_j_cr_over_j_cr0")
        ideal_c_e = _finite_number(sample, "ideal_c_e_over_c_e0")
        observed_c_e = _finite_number(sample, "observed_c_e_over_c_e0")
        key = (alpha_h, projected_j)
        if key in measured_keys:
            raise ContractError(f"duplicate Q-029 synthetic sample {key!r}")
        measured_keys.append(key)
        chi_h = alpha_h * j_cr0 / c_e0
        expected_c_e = ideal_c_e + chi_h * projected_j
        absolute_error = abs(observed_c_e - expected_c_e)
        reports.append({
            "alpha_h": alpha_h,
            "chi_h": chi_h,
            "projected_j_cr_over_j_cr0": projected_j,
            "ideal_c_e_over_c_e0": ideal_c_e,
            "observed_c_e_over_c_e0": observed_c_e,
            "expected_c_e_over_c_e0": expected_c_e,
            "absolute_error": absolute_error,
            "synthetic_contract_consistent": absolute_error <= ABSOLUTE_TOLERANCE,
        })
    if set(measured_keys) != expected_keys:
        raise ContractError(
            "Q-029 synthetic grid is incomplete or contains extra samples"
        )

    reports.sort(key=lambda item: (item["alpha_h"], item["projected_j_cr_over_j_cr0"]))
    consistent = all(item["synthetic_contract_consistent"] for item in reports)
    return {
        "schema_version": 1,
        "campaign_id": CAMPAIGN_ID,
        "artifact_role": ARTIFACT_ROLE,
        "qualification_effect": QUALIFICATION_EFFECT,
        "qualifying_evidence": False,
        "launch_status": LAUNCH_STATUS,
        "source_bindings": source_bindings,
        "deck_contract": deck,
        "sample_count": len(reports),
        "samples": reports,
        "synthetic_normalization_contract_consistent": consistent,
        "status": (
            "source_local_contract_consistent_not_qualifying_evidence"
            if consistent
            else "source_local_contract_mismatch_not_qualifying_evidence"
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
