"""Bounded Q-011 Section 5.4 paper-shock preparation contract."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path

logger = logging.getLogger("athena" + __name__[7:])

_REPO_ROOT = Path(__file__).resolve().parents[3]
_INPUT_DECK = (
    _REPO_ROOT
    / "inputs"
    / "publication"
    / "pic_parallel_shock_section54_paper.athinput"
)
_RESULTS = {}

_EXPECTED_VALUES = {
    ("mesh", "nx1"): "4000",
    ("mesh", "x1min"): "0.0",
    ("mesh", "x1max"): "48000.0",
    ("mesh", "ix1_bc"): "reflect",
    ("mesh", "ox1_bc"): "inflow",
    ("mesh", "nx2"): "260",
    ("mesh", "x2min"): "0.0",
    ("mesh", "x2max"): "3120.0",
    ("mesh", "ix2_bc"): "periodic",
    ("mesh", "ox2_bc"): "periodic",
    ("mesh", "nx3"): "1",
    ("meshblock", "nx1"): "20",
    ("meshblock", "nx2"): "20",
    ("mesh_refinement", "refinement"): "adaptive",
    ("mesh_refinement", "num_levels"): "3",
    ("time", "tlim"): "1200.0",
    ("mhd", "eos"): "ideal",
    ("mhd", "gamma"): "1.66666666667",
    ("particles", "particle_type"): "cosmic_ray",
    ("particles", "pusher"): "boris_tsc",
    ("particles", "pic_enable_2d3v"): "true",
    ("particles", "deposit_moments"): "true",
    ("particles", "deposit_qscale"): "9.0e-4",
    ("particles", "couple_moments_to_mhd"): "true",
    ("particles", "couple_j_to_efield_representation"): "cell_centered",
    ("particles", "couple_j_deposition_mode"): "cc_convert",
    ("particles", "couple_fluid_feedback_order"): "mhd_src_terms",
    ("particles", "couple_moments_momentum_to_mhd"): "true",
    ("particles", "couple_moments_energy_to_mhd"): "true",
    ("particles", "pic_physical_mode"): "paper_mhd_pic",
    ("particles", "pic_background_mode"): "coupled",
    ("particles", "pic_feedback_mode"): "coupled",
    ("particles", "pic_cr_light_speed"): "10000.0",
    ("particles", "pic_cr_initial_state"): "momentum",
    ("particles", "pic_cr_hall_mode"): "off",
    ("particles", "pic_wave_damping_mode"): "off",
    ("particles", "pic_theta_max"): "0.3",
    ("particles", "pic_load_balance_cost_per_particle"): "0.0",
    ("problem", "pgen_name"): "pic_parallel_shock",
    ("problem", "ps_rho0"): "1.0",
    ("problem", "ps_p0"): "0.10",
    ("problem", "ps_u0"): "30.0",
    ("problem", "ps_b0"): "1.0",
    ("problem", "ps_eta"): "1.0e-3",
    ("problem", "ps_vinj_over_u0"): "3.16227766017",
    ("problem", "ps_inject_t_start"): "0.0",
    ("problem", "ps_enable_injection"): "true",
    ("problem", "ps_enable_gas_subtraction"): "true",
    ("problem", "ps_enable_curvature_amr"): "true",
    ("problem", "ps_refine_curv"): "1.0",
    ("problem", "ps_derefine_curv"): "0.1",
    ("problem", "ps_inject_seed"): "23050101",
    ("problem", "ps_enable_frame_tracking"): "false",
    ("output1", "file_type"): "bin",
    ("output1", "variable"): "mhd_w_d",
    ("output1", "id"): "rho",
    ("output1", "dt"): "100.0",
    ("output2", "file_type"): "bin",
    ("output2", "variable"): "mhd_bmag",
    ("output2", "id"): "bmag",
    ("output2", "dt"): "100.0",
    ("output3", "file_type"): "bin",
    ("output3", "variable"): "mhd_j2",
    ("output3", "id"): "j2",
    ("output3", "dt"): "100.0",
    ("output4", "file_type"): "pvtk",
    ("output4", "variable"): "prtcl_all",
    ("output4", "id"): "prtcl_all",
    ("output4", "dt"): "100.0",
    ("output5", "file_type"): "rst",
    ("output5", "dt"): "100.0",
}

_OPEN_ITEMS = [
    "exact_shock_surface_injection_distribution_audit",
    "gas_pressure_and_unit_normalization_audit",
    "downstream_40_640_ppc_macro_mass_calibration",
    "frontier_load_balance_cost_tuning",
    "snapshot_time_selection_tolerance",
    "downstream_spectrum_fit_energy_interval",
    "amr_vs_fine_uniform_residual_tolerances",
    "clean_candidate_frontier_executable_orion_root_and_campaign_execution",
    "independent_raw_artifact_recompute_and_external_review",
]


class ContractError(ValueError):
    """Raised when the frozen Q-011 preparation contract is not satisfied."""


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_athinput(path: Path = _INPUT_DECK) -> dict[str, dict[str, str]]:
    """Parse the strict subset of Athena input syntax used by the frozen deck."""
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


def _require_exact_values(blocks: dict[str, dict[str, str]]) -> None:
    mismatches = []
    for (block, name), expected in _EXPECTED_VALUES.items():
        measured = blocks.get(block, {}).get(name)
        if measured != expected:
            mismatches.append(
                f"{block}/{name}: expected {expected!r}, measured {measured!r}"
            )
    if mismatches:
        raise ContractError("Q-011 deck contract mismatch:\n" + "\n".join(mismatches))


def validate_deck(path: Path = _INPUT_DECK) -> dict[str, object]:
    """Validate sourced values and explicit AthenaK preparation choices."""
    blocks = parse_athinput(path)
    _require_exact_values(blocks)
    dx_root = (
        (float(blocks["mesh"]["x1max"]) - float(blocks["mesh"]["x1min"]))
        / float(blocks["mesh"]["nx1"])
    )
    dy_root = (
        (float(blocks["mesh"]["x2max"]) - float(blocks["mesh"]["x2min"]))
        / float(blocks["mesh"]["nx2"])
    )
    levels = int(blocks["mesh_refinement"]["num_levels"])
    cell_sizes = [dx_root / (2**level) for level in range(levels)]
    ua0 = float(blocks["problem"]["ps_b0"]) / float(
        blocks["problem"]["ps_rho0"]
    ) ** 0.5
    mach_alfven = float(blocks["problem"]["ps_u0"]) / ua0
    light_speed_over_ua0 = float(blocks["particles"]["pic_cr_light_speed"]) / ua0
    if dx_root != 12.0 or dy_root != 12.0 or cell_sizes != [12.0, 6.0, 3.0]:
        raise ContractError("Q-011 root/fine cell-size contract mismatch")
    if mach_alfven != 30.0 or light_speed_over_ua0 != 10000.0:
        raise ContractError("Q-011 dimensionless normalization contract mismatch")
    return {
        "deck": str(path.relative_to(_REPO_ROOT)),
        "deck_sha256": _sha256(path),
        "root_cell_size_c_over_omega_pi": dx_root,
        "amr_cell_sizes_c_over_omega_pi": cell_sizes,
        "mach_alfven": mach_alfven,
        "light_speed_over_ua0": light_speed_over_ua0,
    }


def build_preparation_contract() -> dict[str, object]:
    """Return the frozen contract without inspecting campaign output."""
    return {
        "schema_version": 1,
        "gate": "Q-011",
        "claim_id": "CLAIM-PAPER-SHOCK-001",
        "qualification_effect": "source_controlled_preparation_only",
        "deck": validate_deck(),
        "paper_text_values": {
            "geometry": "reflecting_left_wall_parallel_Bx_periodic_y",
            "domain_c_over_omega_pi": [48000.0, 3120.0],
            "dimensions": "2D3V",
            "mach_alfven": 30.0,
            "gas_gamma": 5.0 / 3.0,
            "injection_efficiency": 1.0e-3,
            "injection_p_over_m_over_u0": 10.0**0.5,
            "light_speed_over_ua0": 10000.0,
            "exclude_birth_time_before_omega0_inverse": 45.0,
            "amr_cell_sizes_c_over_omega_pi": [12.0, 6.0, 3.0],
            "amr_refine_curvature_threshold": 1.0,
            "amr_derefine_curvature_threshold": 0.1,
            "snapshot_times_omega0_inverse": [500.0, 1200.0],
        },
        "artifact_manifest_contract": {
            "required_candidate_bindings": [
                "clean_candidate_git_commit",
                "clean_frontier_executable_sha256",
                "input_deck_sha256",
                "analyzer_sha256",
                "authorized_orion_artifact_root",
            ],
            "qualifying_seed_list": [
                23050101,
                23050102,
                23050103,
                23050104,
                23050105,
                23050106,
                23050107,
                23050108,
            ],
            "required_grid_variants": [
                "coarse_uniform",
                "three_level_amr_root_12_finest_3",
                "fine_uniform",
            ],
            "required_snapshot_times_omega0_inverse": [500.0, 1200.0],
            "required_raw_artifacts_per_snapshot": [
                "rho_bin",
                "bmag_bin",
                "j2_bin",
                "prtcl_all_pvtk",
            ],
            "required_run_artifacts": [
                "stdout_with_q017_telemetry",
                "restart_checkpoint",
                "attempt_status_and_failure_artifacts",
            ],
            "particle_filters": {
                "source": "shock_injected",
                "birth_time_min_omega0_inverse": 45.0,
                "region": "downstream",
            },
            "primary_observables": [
                "shock_position",
                "upstream_magnetic_amplification_at_t500",
                "downstream_weighted_energy_spectra_at_t500_and_t1200",
                "late_energy_spectrum_power_law_slope_at_t1200",
                "morphology",
            ],
        },
        "campaign_execution": {
            "platform": "Frontier",
            "bulk_artifact_root": "/lustre/orion/ast207/proj-shared/dfielding/PIC",
            "status": "open_not_executed_by_preparation_tranche",
            "forbidden_storage_systems": ["Kronos"],
        },
        "open_items": list(_OPEN_ITEMS),
        "result_metrics": [],
    }


def analyze_campaign_artifacts(_manifest: dict[str, object]) -> None:
    """Reject output inspection until the remaining preregistration fields close."""
    raise ContractError(
        "Q-011 campaign analysis is blocked: Frontier execution and pre-run "
        "contract fields remain open"
    )


def run(**kwargs) -> None:
    """Run the local static deck and fail-closed analyzer-contract regression."""
    logger.debug("Running test " + __name__)
    _RESULTS.clear()
    contract = build_preparation_contract()
    _RESULTS["deck_contract"] = bool(contract["deck"]["deck_sha256"])
    _RESULTS["frontier_campaign_execution_open"] = (
        contract["campaign_execution"]["status"]
        == "open_not_executed_by_preparation_tranche"
    )
    try:
        analyze_campaign_artifacts({})
    except ContractError:
        _RESULTS["artifact_analysis_fail_closed"] = True
    else:
        _RESULTS["artifact_analysis_fail_closed"] = False


def analyze() -> bool:
    """Report only preparation-regression status, never fabricated science results."""
    logger.info("Q-011 Section 5.4 preparation contract: %s", _RESULTS)
    return _RESULTS == {
        "deck_contract": True,
        "frontier_campaign_execution_open": True,
        "artifact_analysis_fail_closed": True,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-contract", action="store_true")
    parser.add_argument("--artifact-manifest", type=Path)
    args = parser.parse_args()
    contract = build_preparation_contract()
    if args.artifact_manifest is not None:
        manifest = json.loads(args.artifact_manifest.read_text(encoding="utf-8"))
        analyze_campaign_artifacts(manifest)
    print(json.dumps(contract, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
