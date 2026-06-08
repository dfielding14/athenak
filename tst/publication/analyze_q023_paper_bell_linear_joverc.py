#!/usr/bin/env python3
"""Strict non-authorizing corrected Q023 Bell linear predecessor contract."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
from typing import Any, Mapping, Sequence

import numpy as np

from tst.publication import analyze_q011_section54_outputs as binary
from tst.publication import analyze_q023_paper_bell_linear as legacy


REPO_ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN_ID = "Q023-PAPER-BELL-LINEAR-JOVERC"
SUPERSEDES_CAMPAIGN_ID = "Q023-PAPER-BELL-LINEAR"
PGEN_NAME = "q023_paper_bell_linear_joverc"
SCHEMA_VERSION = 1
QUALIFICATION_EFFECT = "none_non_authorizing_predecessor_preparation_only"
Q043_CURRENT_CAMPAIGN_ID = "Q043-BELL-CURRENT-VOLUME-AWARE"
Q043_RAW_ORACLE_ID = "Q043-BELL-DEPOSITED-J-OVER-C-VOLUME-AWARE"
Q043_GENERATOR_NAME = "q043_bell_current_volume_aware"
Q043_FOUNDATIONAL_BINDING_STATUS = (
    "source_deck_raw_oracle_sequence_bound_registered_admission_hardening_pending"
)
Q043_FOUNDATIONAL_LINEAGE_COMMIT = (
    "883e4679ae797392f25f86e46b952feef753e6aa"
)
Q043_CURRENT_SOURCE_SHA256 = (
    "dd44fff4fe39bf1640dad6019ded29ce7d434b50a0c9831a8749383c6241f3e3"
)
Q043_RAW_ORACLE_ANALYZER_SHA256 = (
    "f0a8db7fb47798d4bdc9ec742cef09e5ae0c69ff3ccb2f25d065602b23a26bed"
)
Q043_RAW_ORACLE_DECK_MANIFEST_SHA256 = (
    "77d6ec862314d308c41490dde1137f6c384d1da24b6b72ad43e357545cff98b0"
)
Q043_SUPERSESSION_SHA256 = (
    "2f83afc03e83038db168fed7d40bdc85bc7bca1eccbba52d18ede0d319dce342"
)
Q043_REGISTERED_ADMISSION_BINDING_STATUS = (
    "pending_hardened_successor_digest_and_schema"
)
Q043_REQUIRED_CASE_COUNT = 132
SOURCE_PATH = Path("src/pgen/tests/q023_paper_bell_linear_joverc.cpp")
CORRECTED_EIGENMODE_HEADER_PATH = Path(
    "src/pgen/tests/q023_paper_bell_linear_joverc.hpp"
)
DECK_ROOT = REPO_ROOT / "inputs/tests/q023_paper_bell_linear_joverc_predecessor"
DECK_MANIFEST = DECK_ROOT / "deck_manifest.json"
DIMENSIONS = (1, 2, 3)
EPSILON_VALUES = legacy.EPSILON_VALUES
PPC = 1
SPECIES_MASS = 2.0
SPECIES_CHARGE = 4.0 * math.pi * 1.0e-6
SPECIES_Q_OVER_MC = SPECIES_CHARGE / SPECIES_MASS
STREAM_SPEED = 2.5
B_G = 1.0
K0 = 2.0 * math.pi
EXPECTED_J_OVER_C = 2.0 * B_G * K0
CURRENT_NORMALIZATION = "deposited_j_over_c_equals_2_b_g_k0"
DEPOSITION_MEASURE = (
    "ppc_times_deposit_qscale_times_species_charge_times_species0_actual_v_cr_"
    "over_root_cell_volume"
)
PARTICLE_VELOCITY_SEMANTICS = (
    "species0_vx0_vy0_vz0_are_actual_and_must_match_particles_cr_velocity"
)
ROOT_CELL_VOLUME_DEFINITION = "global_root_mesh_extents_over_global_root_mesh_counts"
DEPOSIT_QSCALE_SEMANTICS = "root_cell_macro_charge_volume_aware"
FIXED_QSCALE_POLICY = "forbidden_across_dimensions_and_resolutions"
Q_OVER_MC_REPRESENTATION = "species_charge_over_species_mass_equals_omega_over_b_g"
DECOMPOSITION_CONTRACT = "dimension_appropriate_multidirectional_matrix_required"
RAW_OUTPUT_LAYOUT = "shared_mpi_io"
FOUNDATIONAL_ADMISSION_REQUIREMENT = (
    "required_before_any_q023_execution_or_qualification"
)
DEPENDENCY_EFFECT = "non_authorizing_downstream_predecessor_dependency_only"
DECK_ROLE = "non_authorizing_q023_joverc_linear_predecessor"
GROWTH_CONVERGENCE_ABSOLUTE_TOLERANCE = 0.02
PHASE_CONVERGENCE_ABSOLUTE_TOLERANCE = 0.02
DECOMPOSITION_GROWTH_ABSOLUTE_TOLERANCE = 0.01
DECOMPOSITION_PHASE_ABSOLUTE_TOLERANCE = 0.01
DECOMPOSITION_POLARIZATION_RELATIVE_TOLERANCE = 0.05
MIN_GROWTH_FIT_R2 = 0.995
MIN_PHASE_FIT_R2 = 0.995
LINEAR_RUNTIME_TLIM = 1.75
LINEAR_OUTPUT_DT = 0.02
LINEAR_NLIM = 2_000_000
DEPOSITION_TOLERANCE_ULPS = 16.0
FROZEN_CURRENT_GYROPERIOD_MARGIN = 4.0

_BOUNDS = {
    1: ((0.0, 1.0), (0.0, 1.0), (0.0, 1.0)),
    2: ((0.0, math.sqrt(5.0)), (0.0, math.sqrt(1.25)), (0.0, 1.0)),
    3: (
        (0.0, math.sqrt(21.0)),
        (0.0, math.sqrt(5.25)),
        (0.0, math.sqrt(1.3125)),
    ),
}
_ROOT_NX = {
    1: {"coarse": (16, 4, 1), "fine": (32, 4, 1)},
    2: {"coarse": (32, 16, 1), "fine": (64, 32, 1)},
    3: {"coarse": (64, 32, 16), "fine": (128, 64, 32)},
}
DECOMPOSITIONS_BY_DIMENSION = {
    1: {"reference_x1": (2, 1, 1)},
    2: {
        "reference_x1": (2, 1, 1),
        "split_x2": (1, 2, 1),
        "split_x1x2": (2, 2, 1),
    },
    3: {
        "reference_x1": (2, 1, 1),
        "split_x2": (1, 2, 1),
        "split_x3": (1, 1, 2),
        "split_xyz": (2, 2, 2),
    },
}
PGEN_REQUIRED_STRINGS = {
    "campaign_id": CAMPAIGN_ID,
    "supersedes_campaign_id": SUPERSEDES_CAMPAIGN_ID,
    "foundational_current_campaign_id": Q043_CURRENT_CAMPAIGN_ID,
    "foundational_raw_oracle_id": Q043_RAW_ORACLE_ID,
    "foundational_binding_status": Q043_FOUNDATIONAL_BINDING_STATUS,
    "foundational_current_lineage_commit": Q043_FOUNDATIONAL_LINEAGE_COMMIT,
    "foundational_current_source_sha256": Q043_CURRENT_SOURCE_SHA256,
    "foundational_raw_oracle_analyzer_sha256": Q043_RAW_ORACLE_ANALYZER_SHA256,
    "foundational_raw_oracle_deck_manifest_sha256": (
        Q043_RAW_ORACLE_DECK_MANIFEST_SHA256
    ),
    "foundational_supersession_sha256": Q043_SUPERSESSION_SHA256,
    "foundational_registered_admission_binding_status": (
        Q043_REGISTERED_ADMISSION_BINDING_STATUS
    ),
    "foundational_registered_admission_binding_permitted": (
        "false_until_hardened_successor_lands"
    ),
    "foundational_raw_oracle_case_count": str(Q043_REQUIRED_CASE_COUNT),
    "foundational_raw_oracle_registered_pass": (
        "required_before_linear_qualification"
    ),
    "foundational_registered_execution_admission": FOUNDATIONAL_ADMISSION_REQUIREMENT,
    "dependency_effect": DEPENDENCY_EFFECT,
    "deck_role": DECK_ROLE,
    "current_normalization": CURRENT_NORMALIZATION,
    "deposition_measure": DEPOSITION_MEASURE,
    "particle_velocity_semantics": PARTICLE_VELOCITY_SEMANTICS,
    "root_cell_volume": ROOT_CELL_VOLUME_DEFINITION,
    "deposit_qscale_semantics": DEPOSIT_QSCALE_SEMANTICS,
    "fixed_qscale": FIXED_QSCALE_POLICY,
    "decomposition_contract": DECOMPOSITION_CONTRACT,
    "raw_output_layout": RAW_OUTPUT_LAYOUT,
    "epsilon_grid": "0.1,0.2,0.4,0.6,0.8",
    "timestep": "open_clean_candidate_timestep_freeze",
    "source_mode": "corrected_linear_eigenmode",
    "q_over_mc_representation": Q_OVER_MC_REPRESENTATION,
    "initial_eigenmode": "section52_right_polarized_eigenmode",
    "qualification_effect": QUALIFICATION_EFFECT,
}
_SHA256 = re.compile(r"[0-9a-f]{64}")
_COMMIT = re.compile(r"[0-9a-f]{40}")
_Q043_DEPENDENCY_KEYS = {
    "schema_version",
    "record_type",
    "binding_kind",
    "current_campaign_id",
    "raw_oracle_id",
    "generator_name",
    "integration_checkpoint_commit",
    "current_source_sha256",
    "raw_oracle_analyzer_sha256",
    "raw_oracle_deck_manifest_sha256",
    "supersession_sha256",
    "registered_admission_binding_status",
    "registered_admission_digest_bound",
    "registered_admission_schema_bound",
    "registered_execution_qualification_check_pass",
    "registered_raw_oracle_pass",
    "complete_foundational_raw_oracle_matrix_pass",
    "required_case_count",
    "measured_case_count",
    "qualification_effect",
    "launch_authorized",
    "scientific_claim_authorized",
    "publication_authorized",
}
_MATRIX_RECORD_KEYS = {
    "member_id",
    "dimension",
    "epsilon",
    "resolution",
    "decomposition",
    "decomposition_splits",
    "physics_trace",
    "provenance",
}
_PHYSICS_TRACE_KEYS = legacy._TRACE_KEYS - {"dimension", "epsilon", "raw_provenance"}
_PROVENANCE_KEYS = {
    "kind",
    "deck_path",
    "deck_sha256",
    "source_path",
    "source_sha256",
    "corrected_eigenmode_header_path",
    "corrected_eigenmode_header_sha256",
    "q043_registered_raw_oracle_dependency_sha256",
    "authorized_artifact_root",
    "candidate_clean",
    "source_clean",
    "executable_clean",
    "executable_path",
    "executable_sha256",
    "registered_execution_receipt_path",
    "registered_execution_receipt_sha256",
    "raw_artifacts",
}
_ARTIFACT_KEYS = {"path", "sha256"}
_RAW_ARTIFACT_KEYS = {"path", "sha256", "variable", "cycle", "time"}
_MHD_FIELDS = frozenset(
    ("dens", "eint", "velx", "vely", "velz", "bcc1", "bcc2", "bcc3")
)
_RAW_OUTPUT_VARIABLES = (
    "mhd_w_bcc",
    "prtcl_rho",
    "prtcl_jx",
    "prtcl_jy",
    "prtcl_jz",
)
_PARTICLE_FIELDS = _RAW_OUTPUT_VARIABLES[1:]
_RAW_OUTPUT_INDEX = {
    variable: index for index, variable in enumerate(_RAW_OUTPUT_VARIABLES, 1)
}
_OUTPUT_BOOKKEEPING_KEYS = frozenset(("file_number", "last_time"))
_EXECUTION_RECEIPT_KEYS = {
    "schema_version",
    "record_type",
    "campaign_id",
    "member_id",
    "deck_path",
    "deck_sha256",
    "source_path",
    "source_sha256",
    "corrected_eigenmode_header_path",
    "corrected_eigenmode_header_sha256",
    "executable_path",
    "executable_sha256",
    "candidate_clean",
    "source_clean",
    "executable_clean",
    "command",
    "mpi",
    "termination",
    "stdout",
    "raw_artifacts",
    "publication_authorized",
}
_COMMAND_KEYS = {"argv", "working_directory"}
_MPI_KEYS = {
    "launcher",
    "rank_count",
    "decomposition_splits",
    "rank_topology_path",
    "rank_topology_sha256",
}
_TERMINATION_KEYS = {"status", "return_code", "final_cycle", "final_time"}
_STDOUT_KEYS = {"path", "sha256", "completion_marker"}
_RANK_TOPOLOGY_KEYS = {
    "schema_version",
    "record_type",
    "campaign_id",
    "member_id",
    "rank_count",
    "decomposition_splits",
    "ranks",
    "publication_authorized",
}
_RANK_EVIDENCE_KEYS = {"rank", "meshblock_logical_locations"}


class ContractError(ValueError):
    """Raised when the corrected Q023 predecessor contract fails closed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def _canonical_json_bytes(value: object) -> bytes:
    return (
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
        + "\n"
    ).encode("utf-8")


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256_file(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


def _float_token(value: float) -> str:
    return format(value, ".17g")


def _positive_integral_ppc(value: object) -> bool:
    try:
        measured = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(measured) and measured > 0.0 and measured.is_integer()


def _mode_basis(dimension: int) -> tuple[float, float, float]:
    raw = (1.0, 2.0 if dimension >= 2 else 0.0, 4.0 if dimension >= 3 else 0.0)
    norm = math.sqrt(sum(value * value for value in raw))
    return tuple(value / norm for value in raw)


def _root_cell_volume(
    bounds: Sequence[Sequence[float]], nx: Sequence[int]
) -> float:
    return math.prod(
        (float(upper) - float(lower)) / int(count)
        for (lower, upper), count in zip(bounds, nx)
    )


def required_deposit_qscale(
    *,
    root_cell_volume: float,
    ppc: int,
    species_charge: float = SPECIES_CHARGE,
    stream_speed: float = STREAM_SPEED,
) -> float:
    """Return qscale for deposited J_CR/c = 2 B_g k0."""
    _require(
        root_cell_volume > 0.0
        and _positive_integral_ppc(ppc)
        and species_charge > 0.0
        and stream_speed > 0.0,
        "normalization factors must be positive and PPC must be a positive integer",
    )
    return (
        EXPECTED_J_OVER_C
        * root_cell_volume
        / (int(ppc) * species_charge * stream_speed)
    )


def _epsilon_token(epsilon: float) -> str:
    return f"{epsilon:.1f}".replace(".", "p")


def expected_deck_members() -> tuple[dict[str, object], ...]:
    """Return the exact two-resolution, multidirectional decomposition matrix."""
    members: list[dict[str, object]] = []
    for dimension in DIMENSIONS:
        for resolution in ("coarse", "fine"):
            decompositions = (
                {"reference_x1": DECOMPOSITIONS_BY_DIMENSION[dimension]["reference_x1"]}
                if resolution == "coarse"
                else DECOMPOSITIONS_BY_DIMENSION[dimension]
            )
            nx = _ROOT_NX[dimension][resolution]
            bounds = _BOUNDS[dimension]
            volume = _root_cell_volume(bounds, nx)
            for decomposition, splits in decompositions.items():
                for epsilon in EPSILON_VALUES:
                    stream_speed = 1.0 / epsilon
                    member_id = (
                        f"d{dimension}-{resolution}-{decomposition}-"
                        f"eps{_epsilon_token(epsilon)}"
                    )
                    members.append(
                        {
                            "member_id": member_id,
                            "dimension": dimension,
                            "epsilon": epsilon,
                            "resolution": resolution,
                            "decomposition": decomposition,
                            "decomposition_splits": list(splits),
                            "global_nx": list(nx),
                            "meshblock_nx": [
                                count // split for count, split in zip(nx, splits)
                            ],
                            "bounds": [list(axis) for axis in bounds],
                            "root_cell_volume": volume,
                            "ppc": PPC,
                            "species_mass": SPECIES_MASS,
                            "species_charge": SPECIES_CHARGE,
                            "stream_speed": stream_speed,
                            "artificial_light_speed": 1000.0 * stream_speed,
                            "deposit_qscale": required_deposit_qscale(
                                root_cell_volume=volume,
                                ppc=PPC,
                                stream_speed=stream_speed,
                            ),
                        }
                    )
    _require(len(members) == 55, "Q023 predecessor deck matrix size drifted")
    return tuple(members)


def render_deck(member: Mapping[str, object]) -> str:
    """Render one exact non-authorizing corrected linear predecessor deck."""
    dimension = int(member["dimension"])
    nx = tuple(int(value) for value in member["global_nx"])
    mb = tuple(int(value) for value in member["meshblock_nx"])
    bounds = tuple(tuple(float(value) for value in axis) for axis in member["bounds"])
    basis = _mode_basis(dimension)
    stream_speed = float(member["stream_speed"])
    stream = tuple(stream_speed * value for value in basis)
    lines = [
        "# Corrected Q023 Bell linear predecessor preparation only.",
        "# Q043 registered raw-current foundation is required before execution.",
        "# This deck grants no launch, qualification, science, or publication authority.",
        "",
        "<comment>",
        f"problem = {CAMPAIGN_ID} {member['member_id']}",
        "",
        "<job>",
        f"basename = q023_joverc_{str(member['member_id']).replace('-', '_')}",
        "",
        "<mesh>",
        "nghost = 2",
    ]
    for axis, ((lower, upper), count) in enumerate(zip(bounds, nx), 1):
        lines.extend(
            [
                f"nx{axis} = {count}",
                f"x{axis}min = {_float_token(lower)}",
                f"x{axis}max = {_float_token(upper)}",
                f"ix{axis}_bc = periodic",
                f"ox{axis}_bc = periodic",
                "",
            ]
        )
    lines.extend(
        [
            "<meshblock>",
            f"nx1 = {mb[0]}",
            f"nx2 = {mb[1]}",
            f"nx3 = {mb[2]}",
            "",
            "<time>",
            "evolution = dynamic",
            "integrator = rk2",
            "cfl_number = 0.1",
            f"nlim = {LINEAR_NLIM}",
            f"tlim = {_float_token(LINEAR_RUNTIME_TLIM)}",
            "ndiag = 1",
            "",
            "<mhd>",
            "eos = ideal",
            "reconstruct = plm",
            "rsolver = llf",
            "gamma = 1.66666666667",
            "",
            "<particles>",
            "particle_type = cosmic_ray",
            f"ppc = {PPC}.0",
            "pusher = boris_tsc",
            "nspecies = 1",
            "cr_distribution = center",
            "deposit_moments = true",
            "deposit_order = 2",
            f"deposit_qscale = {_float_token(float(member['deposit_qscale']))}",
            "couple_moments_to_mhd = true",
            "couple_j_to_efield_coeff = 1.0",
            "couple_j_to_efield_representation = cell_centered",
            "couple_j_deposition_mode = cc_convert",
            "couple_moments_momentum_to_mhd = true",
            "couple_moments_energy_to_mhd = true",
            "couple_fluid_feedback_order = mhd_src_terms",
            "couple_moments_momentum_coeff = 1.0",
            "couple_moments_energy_coeff = 1.0",
            f"cr_vx0 = {_float_token(stream[0])}",
            f"cr_vy0 = {_float_token(stream[1])}",
            f"cr_vz0 = {_float_token(stream[2])}",
            "pic_physical_mode = paper_mhd_pic_vl2_tsc",
            "pic_background_mode = coupled",
            "pic_feedback_mode = coupled",
            "pic_interp_scheme = tsc",
            "pic_enable_2d3v = true",
            f"pic_cr_light_speed = {_float_token(float(member['artificial_light_speed']))}",
            "pic_cr_initial_state = velocity",
            "pic_cr_hall_mode = off",
            "pic_wave_damping_mode = off",
            "pic_max_cell_cross = 1",
            "pic_theta_max = 0.3",
            "pic_deltaf_mode = off",
            "pic_load_balance_cost_per_particle = 0.0",
            "pic_sort_interval = 0",
            "pic_intermediate_arrays = auto",
            "pic_expanding_box_mode = off",
            "",
            "<species0>",
            f"mass = {_float_token(SPECIES_MASS)}",
            f"charge = {_float_token(SPECIES_CHARGE)}",
            f"vx0 = {_float_token(stream[0])}",
            f"vy0 = {_float_token(stream[1])}",
            f"vz0 = {_float_token(stream[2])}",
            "",
            "<problem>",
            f"pgen_name = {PGEN_NAME}",
            "",
            f"<{PGEN_NAME}>",
        ]
    )
    for name, value in PGEN_REQUIRED_STRINGS.items():
        lines.append(f"{name} = {value}")
    lines.extend(
        [
            f"dimension = {dimension}",
            f"resolution_id = {member['resolution']}",
            f"decomposition_id = {member['decomposition']}",
            (
                "prepared_meshblock_splits = "
                + ",".join(str(value) for value in member["decomposition_splits"])
            ),
            f"epsilon_default = {_float_token(float(member['epsilon']))}",
            f"epsilon = {_float_token(float(member['epsilon']))}",
            "rho = 1.0",
            "pressure = 1.0",
            "amplitude = 1.0e-6",
            "b_g = 1.0",
            "u_a = 1.0",
            "wavelength = 1.0",
            f"k0 = {_float_token(K0)}",
            f"omega = {_float_token(SPECIES_Q_OVER_MC)}",
            "c_over_v_cr = 1000.0",
            "supersession_effect = historical_c_multiplied_normalization_invalid_for_qualification",
            "launch_authorized = false",
            "qualification_eligible = false",
            "",
        ]
    )
    for index, field in enumerate(("mhd_w_bcc", "prtcl_rho", "prtcl_jx", "prtcl_jy", "prtcl_jz"), 1):
        lines.extend(
            [
                f"<output{index}>",
                "file_type = bin",
                f"variable = {field}",
                f"id = {field}",
                f"dt = {_float_token(LINEAR_OUTPUT_DT)}",
                "ghost_zones = false",
                "single_file_per_rank = false",
                "",
            ]
        )
    lines.extend(
        [
            "<output6>",
            "file_type = rst",
            "id = rst",
            f"dt = {_float_token(LINEAR_RUNTIME_TLIM)}",
            "single_file_per_rank = false",
            "",
        ]
    )
    return "\n".join(lines)


def parse_athinput_text(text: str) -> dict[str, dict[str, str]]:
    blocks: dict[str, dict[str, str]] = {}
    block: str | None = None
    for line_number, raw_line in enumerate(text.splitlines(), 1):
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            block = line[1:-1].strip()
            _require(block and block not in blocks, f"line {line_number}: invalid block")
            blocks[block] = {}
            continue
        _require(block is not None and "=" in line, f"line {line_number}: malformed")
        name, value = (part.strip() for part in line.split("=", 1))
        _require(name and value and name not in blocks[block], f"line {line_number}: duplicate")
        blocks[block][name] = value
    return blocks


def validate_rendered_deck(member: Mapping[str, object], text: str) -> dict[str, object]:
    blocks = parse_athinput_text(text)
    _require(blocks["problem"]["pgen_name"] == PGEN_NAME, "corrected pgen drifted")
    for name, expected in PGEN_REQUIRED_STRINGS.items():
        _require(blocks[PGEN_NAME].get(name) == expected, f"pgen string {name} drifted")
    _require(blocks[PGEN_NAME]["launch_authorized"] == "false", "deck authorized launch")
    _require(blocks[PGEN_NAME]["qualification_eligible"] == "false", "deck became qualifying")
    _require(
        blocks["time"]["nlim"] == str(LINEAR_NLIM)
        and float(blocks["time"]["tlim"]) == LINEAR_RUNTIME_TLIM,
        "linear runtime horizon drifted",
    )
    for name in (
        "couple_j_to_efield_coeff",
        "couple_moments_momentum_coeff",
        "couple_moments_energy_coeff",
    ):
        _require(
            math.isclose(
                float(blocks["particles"][name]), 1.0, rel_tol=0.0, abs_tol=0.0
            ),
            f"{name} drifted",
        )
    _require(blocks[PGEN_NAME]["resolution_id"] == member["resolution"], "resolution drifted")
    _require(blocks[PGEN_NAME]["decomposition_id"] == member["decomposition"], "decomposition drifted")
    _require(
        float(blocks[PGEN_NAME]["epsilon"]) == member["epsilon"]
        and float(blocks[PGEN_NAME]["epsilon_default"]) == member["epsilon"],
        "epsilon deck binding drifted",
    )
    expected_splits = ",".join(str(value) for value in member["decomposition_splits"])
    _require(blocks[PGEN_NAME]["prepared_meshblock_splits"] == expected_splits, "split metadata drifted")
    nx = tuple(int(blocks["mesh"][f"nx{axis}"]) for axis in (1, 2, 3))
    mb = tuple(int(blocks["meshblock"][f"nx{axis}"]) for axis in (1, 2, 3))
    _require(nx == tuple(member["global_nx"]), "root-grid drifted")
    _require(mb == tuple(member["meshblock_nx"]), "MeshBlock decomposition drifted")
    measured_splits = tuple(count // block for count, block in zip(nx, mb))
    _require(measured_splits == tuple(member["decomposition_splits"]), "measured split drifted")
    ppc = float(blocks["particles"]["ppc"])
    _require(_positive_integral_ppc(ppc), "PPC must be a positive integer")
    mass = float(blocks["species0"]["mass"])
    charge = float(blocks["species0"]["charge"])
    _require(math.isclose(charge / mass, SPECIES_Q_OVER_MC, rel_tol=1.0e-13), "q/m drifted")
    bounds = tuple(
        (
            float(blocks["mesh"][f"x{axis}min"]),
            float(blocks["mesh"][f"x{axis}max"]),
        )
        for axis in (1, 2, 3)
    )
    volume = _root_cell_volume(bounds, nx)
    qscale = float(blocks["particles"]["deposit_qscale"])
    global_stream = tuple(
        float(blocks["particles"][f"cr_v{axis}0"]) for axis in ("x", "y", "z")
    )
    species_stream = tuple(
        float(blocks["species0"][f"v{axis}0"]) for axis in ("x", "y", "z")
    )
    expected_stream = tuple(
        float(member["stream_speed"]) * component for component in _mode_basis(int(member["dimension"]))
    )
    _require(
        all(
            math.isclose(measured, expected, rel_tol=1.0e-13, abs_tol=1.0e-13)
            for measured, expected in zip(species_stream, expected_stream)
        ),
        "species0 actual initialization velocity drifted",
    )
    _require(
        all(
            math.isclose(global_value, species_value, rel_tol=1.0e-13, abs_tol=1.0e-13)
            for global_value, species_value in zip(global_stream, species_stream)
        ),
        "global CR velocity and species0 actual initialization velocity differ",
    )
    speed = math.sqrt(sum(value**2 for value in species_stream))
    _require(
        math.isclose(speed, float(member["stream_speed"]), rel_tol=1.0e-13)
        and math.isclose(speed, 1.0 / float(member["epsilon"]), rel_tol=1.0e-13)
        and math.isclose(
            float(blocks["particles"]["pic_cr_light_speed"]),
            float(member["artificial_light_speed"]),
            rel_tol=1.0e-13,
        ),
        "epsilon, stream speed, or artificial light speed drifted",
    )
    deposited = int(ppc) * qscale * charge * speed / volume
    _require(math.isclose(deposited, EXPECTED_J_OVER_C, rel_tol=1.0e-13, abs_tol=1.0e-13), "deposited J/c closure drifted")
    for index in range(1, 6):
        _require(
            float(blocks[f"output{index}"]["dt"]) == LINEAR_OUTPUT_DT
            and "dcycle" not in blocks[f"output{index}"],
            "linear output cadence drifted",
        )
        _require(
            blocks[f"output{index}"]["single_file_per_rank"] == "false",
            "linear raw-output shared MPI layout drifted",
        )
    _require(
        float(blocks["output6"]["dt"]) == LINEAR_RUNTIME_TLIM,
        "restart output cadence drifted",
    )
    return {
        "member_id": member["member_id"],
        "deck_sha256": _sha256_bytes(text.encode("utf-8")),
        "actual_species_stream_speed": speed,
        "configured_deposited_j_over_c": deposited,
    }


def build_deck_manifest() -> tuple[dict[str, object], dict[str, str]]:
    decks: dict[str, str] = {}
    records = []
    for member in expected_deck_members():
        relative = f"q023-joverc-{member['member_id']}.athinput"
        text = render_deck(member)
        decks[relative] = text
        records.append({**member, "deck_path": relative, **validate_rendered_deck(member, text)})
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "record_type": "q023_paper_bell_linear_joverc_predecessor_deck_manifest",
        "campaign_id": CAMPAIGN_ID,
        "foundational_binding_status": Q043_FOUNDATIONAL_BINDING_STATUS,
        "foundational_current_lineage_commit": Q043_FOUNDATIONAL_LINEAGE_COMMIT,
        "foundational_current_source_sha256": Q043_CURRENT_SOURCE_SHA256,
        "foundational_raw_oracle_analyzer_sha256": Q043_RAW_ORACLE_ANALYZER_SHA256,
        "foundational_raw_oracle_deck_manifest_sha256": (
            Q043_RAW_ORACLE_DECK_MANIFEST_SHA256
        ),
        "foundational_supersession_sha256": Q043_SUPERSESSION_SHA256,
        "foundational_registered_admission_binding_status": (
            Q043_REGISTERED_ADMISSION_BINDING_STATUS
        ),
        "foundational_registered_admission_digest_bound": False,
        "foundational_registered_admission_schema_bound": False,
        "foundational_raw_oracle_id": Q043_RAW_ORACLE_ID,
        "foundational_raw_oracle_case_count": Q043_REQUIRED_CASE_COUNT,
        "foundational_registered_execution_admission": FOUNDATIONAL_ADMISSION_REQUIREMENT,
        "qualification_effect": QUALIFICATION_EFFECT,
        "launch_authorized": False,
        "qualification_eligible": False,
        "case_count": len(records),
        "decompositions_by_dimension": {
            str(dimension): {
                name: list(splits)
                for name, splits in DECOMPOSITIONS_BY_DIMENSION[dimension].items()
            }
            for dimension in DIMENSIONS
        },
        "cases": records,
    }
    return manifest, decks


def materialize_checked_in_decks(*, replace: bool = False) -> dict[str, object]:
    manifest, decks = build_deck_manifest()
    if DECK_ROOT.exists():
        _require(replace, "checked-in predecessor deck root already exists")
        shutil.rmtree(DECK_ROOT)
    DECK_ROOT.mkdir(parents=True)
    for relative, text in decks.items():
        (DECK_ROOT / relative).write_text(text, encoding="utf-8")
    DECK_MANIFEST.write_bytes(_canonical_json_bytes(manifest))
    return manifest


def validate_checked_in_decks() -> dict[str, object]:
    expected, decks = build_deck_manifest()
    _require(DECK_MANIFEST.is_file(), "checked-in predecessor deck manifest missing")
    _require(json.loads(DECK_MANIFEST.read_text(encoding="utf-8")) == expected, "deck manifest drifted")
    _require(
        {path.name for path in DECK_ROOT.iterdir() if path.is_file()}
        == set(decks) | {DECK_MANIFEST.name},
        "predecessor deck inventory drifted",
    )
    for relative, text in decks.items():
        _require((DECK_ROOT / relative).read_text(encoding="utf-8") == text, f"deck drifted: {relative}")
    return expected


def synthetic_q043_registered_raw_oracle_dependency() -> dict[str, object]:
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": "q043_registered_raw_oracle_dependency",
        "binding_kind": "synthetic_contract_fixture",
        "current_campaign_id": Q043_CURRENT_CAMPAIGN_ID,
        "raw_oracle_id": Q043_RAW_ORACLE_ID,
        "generator_name": Q043_GENERATOR_NAME,
        "integration_checkpoint_commit": Q043_FOUNDATIONAL_LINEAGE_COMMIT,
        "current_source_sha256": Q043_CURRENT_SOURCE_SHA256,
        "raw_oracle_analyzer_sha256": Q043_RAW_ORACLE_ANALYZER_SHA256,
        "raw_oracle_deck_manifest_sha256": Q043_RAW_ORACLE_DECK_MANIFEST_SHA256,
        "supersession_sha256": Q043_SUPERSESSION_SHA256,
        "registered_admission_binding_status": Q043_REGISTERED_ADMISSION_BINDING_STATUS,
        "registered_admission_digest_bound": False,
        "registered_admission_schema_bound": False,
        "registered_execution_qualification_check_pass": False,
        "registered_raw_oracle_pass": False,
        "complete_foundational_raw_oracle_matrix_pass": False,
        "required_case_count": Q043_REQUIRED_CASE_COUNT,
        "measured_case_count": 0,
        "qualification_effect": DEPENDENCY_EFFECT,
        "launch_authorized": False,
        "scientific_claim_authorized": False,
        "publication_authorized": False,
    }


def validate_q043_dependency(
    value: object, *, artifact_root: Path | None = None
) -> dict[str, object]:
    _require(isinstance(value, Mapping) and set(value) == _Q043_DEPENDENCY_KEYS, "Q043 dependency keys drifted")
    dependency = dict(value)
    _require(dependency["schema_version"] == SCHEMA_VERSION, "Q043 dependency schema drifted")
    _require(dependency["record_type"] == "q043_registered_raw_oracle_dependency", "Q043 dependency type drifted")
    _require(
        dependency["binding_kind"] == "synthetic_contract_fixture",
        "hardened Q043 registered-admission digest/schema binding is pending",
    )
    _require(dependency["current_campaign_id"] == Q043_CURRENT_CAMPAIGN_ID, "Q043 current campaign drifted")
    _require(dependency["raw_oracle_id"] == Q043_RAW_ORACLE_ID, "Q043 raw oracle drifted")
    _require(dependency["generator_name"] == Q043_GENERATOR_NAME, "Q043 generator drifted")
    _require(
        isinstance(dependency["integration_checkpoint_commit"], str)
        and _COMMIT.fullmatch(dependency["integration_checkpoint_commit"]) is not None,
        "Q043 integration checkpoint commit malformed",
    )
    exact_foundational_bindings = {
        "integration_checkpoint_commit": Q043_FOUNDATIONAL_LINEAGE_COMMIT,
        "current_source_sha256": Q043_CURRENT_SOURCE_SHA256,
        "raw_oracle_analyzer_sha256": Q043_RAW_ORACLE_ANALYZER_SHA256,
        "raw_oracle_deck_manifest_sha256": Q043_RAW_ORACLE_DECK_MANIFEST_SHA256,
        "supersession_sha256": Q043_SUPERSESSION_SHA256,
    }
    for name, expected in exact_foundational_bindings.items():
        _require(
            isinstance(dependency[name], str)
            and (
                _COMMIT.fullmatch(dependency[name])
                if name == "integration_checkpoint_commit"
                else _SHA256.fullmatch(dependency[name])
            )
            is not None,
            f"Q043 foundational {name} malformed",
        )
        _require(
            dependency[name] == expected,
            f"Q043 foundational {name} does not bind final integration checkpoint "
            f"{Q043_FOUNDATIONAL_LINEAGE_COMMIT}",
        )
    _require(
        dependency["registered_admission_binding_status"]
        == Q043_REGISTERED_ADMISSION_BINDING_STATUS
        and dependency["registered_admission_digest_bound"] is False
        and dependency["registered_admission_schema_bound"] is False,
        "provisional Q043 registered-admission binding is forbidden",
    )
    _require(
        dependency["registered_execution_qualification_check_pass"] is False
        and dependency["registered_raw_oracle_pass"] is False
        and dependency["complete_foundational_raw_oracle_matrix_pass"] is False,
        "pending Q043 registered admission must remain unclaimed",
    )
    _require(
        dependency["required_case_count"] == Q043_REQUIRED_CASE_COUNT
        and dependency["measured_case_count"] == 0,
        "pending Q043 foundational case count status drifted",
    )
    _require(dependency["qualification_effect"] == DEPENDENCY_EFFECT, "Q043 dependency authority drifted")
    _require(dependency["launch_authorized"] is False and dependency["scientific_claim_authorized"] is False and dependency["publication_authorized"] is False, "Q043 dependency improperly grants authority")
    return dependency


def _dependency_digest(dependency: Mapping[str, object]) -> str:
    return _sha256_bytes(_canonical_json_bytes(dependency))


def theoretical_dispersion(epsilon: float) -> tuple[float, float]:
    """Return the signed frequency and growth rate for the retained exp(-ikx) mode."""
    _require(epsilon in EPSILON_VALUES, "epsilon is outside the corrected Q023 grid")
    return -epsilon, math.sqrt(1.0 - epsilon * epsilon)


def _trace_array(trace: Mapping[str, object], name: str) -> np.ndarray:
    values = np.asarray(trace[name], dtype=float)
    _require(
        values.ndim == 1 and values.size > 0 and np.all(np.isfinite(values)),
        f"{name} must be a finite nonempty one-dimensional array",
    )
    return values


def _fit_linear_quality(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    coefficients = np.polyfit(x, y, 1)
    fitted = np.polyval(coefficients, x)
    residual = y - fitted
    ss_res = float(np.sum(residual * residual))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return float(coefficients[0]), 1.0 if ss_tot == 0.0 else 1.0 - ss_res / ss_tot


def _fit_growth_quality(
    time: np.ndarray, amplitude: np.ndarray, expected_growth: float, *, label: str
) -> tuple[float, float, np.ndarray]:
    fit_mask = (expected_growth * time >= 1.0) & (expected_growth * time <= 5.0)
    _require(
        np.count_nonzero(fit_mask) >= legacy.MIN_GROWTH_SNAPSHOTS,
        "the fixed theory growth window requires at least eight snapshots",
    )
    _require(
        np.all(amplitude[fit_mask] > 0.0),
        f"{label} right-polarized amplitude must remain positive",
    )
    measured, quality = _fit_linear_quality(time[fit_mask], np.log(amplitude[fit_mask]))
    return measured, quality, fit_mask


def _validate_phase_trace(
    trace: Mapping[str, object],
    *,
    time: np.ndarray,
    mode: np.ndarray,
    interval_name: str,
    change_name: str,
    label: str,
) -> float:
    interval = _trace_array(trace, interval_name)
    change = _trace_array(trace, change_name)
    _require(interval.size == change.size, f"{label} phase arrays must have equal length")
    _require(
        np.allclose(interval, legacy.PHASE_INTERVAL, rtol=0.0, atol=1.0e-14),
        f"{label} phase intervals drifted",
    )
    extracted_interval, extracted_change = legacy._fixed_interval_phase_trace(time, mode)
    _require(
        interval.size == extracted_interval.size
        and np.allclose(interval, extracted_interval, rtol=0.0, atol=1.0e-12)
        and np.allclose(change, extracted_change, rtol=0.0, atol=1.0e-10),
        f"{label} phase trace does not match the retained mode",
    )
    return float(np.mean(extracted_change / extracted_interval))


def _within_theory_tolerance(measured: float, expected: float) -> bool:
    error = abs(measured - expected)
    return (
        error <= legacy.ABSOLUTE_TOLERANCE
        and error / abs(expected) <= legacy.RELATIVE_TOLERANCE
    )


def _analyze_physics_trace(trace: Mapping[str, object], epsilon: float) -> dict[str, object]:
    """Apply corrected signed-frequency and robust fit-quality gates."""
    _require(set(trace) == _PHYSICS_TRACE_KEYS, "Q023 physics trace keys drifted")
    expected_phase, expected_growth = theoretical_dispersion(epsilon)
    time = _trace_array(trace, "normalized_time")
    _require(
        time.size >= legacy.MIN_GROWTH_SNAPSHOTS and np.all(np.diff(time) > 0.0),
        "normalized_time must be strictly increasing with at least eight rows",
    )

    def complex_trace(real_name: str, imag_name: str) -> np.ndarray:
        real = _trace_array(trace, real_name)
        imag = _trace_array(trace, imag_name)
        _require(
            real.size == time.size and imag.size == time.size,
            f"{real_name} and {imag_name} must match normalized_time",
        )
        return real + 1.0j * imag

    right = complex_trace("right_mode_real", "right_mode_imag")
    left = complex_trace("left_mode_real", "left_mode_imag")
    velocity_right = complex_trace(
        "velocity_right_mode_real", "velocity_right_mode_imag"
    )
    velocity_left = complex_trace(
        "velocity_left_mode_real", "velocity_left_mode_imag"
    )
    paper_mode = complex_trace(
        "paper_delta_u_y_sine_fit_real", "paper_delta_u_y_sine_fit_imag"
    )
    paper_abs = _trace_array(trace, "paper_volume_averaged_abs_delta_u")
    _require(paper_abs.size == time.size, "paper |delta u| must match normalized_time")

    measured_growth, growth_r2, fit_mask = _fit_growth_quality(
        time, np.abs(right), expected_growth, label="magnetic"
    )
    velocity_growth, velocity_growth_r2, velocity_fit_mask = _fit_growth_quality(
        time, np.abs(velocity_right), expected_growth, label="fluid velocity"
    )
    paper_growth, paper_growth_r2, paper_fit_mask = _fit_growth_quality(
        time, paper_abs, expected_growth, label="paper-literal fluid velocity"
    )
    _require(
        np.array_equal(fit_mask, velocity_fit_mask)
        and np.array_equal(fit_mask, paper_fit_mask),
        "growth fit windows differ",
    )

    interval_phase = _validate_phase_trace(
        trace,
        time=time,
        mode=right,
        interval_name="phase_interval",
        change_name="phase_change",
        label="magnetic",
    )
    velocity_interval_phase = _validate_phase_trace(
        trace,
        time=time,
        mode=velocity_right,
        interval_name="velocity_phase_interval",
        change_name="velocity_phase_change",
        label="fluid velocity",
    )
    paper_interval_phase = _validate_phase_trace(
        trace,
        time=time,
        mode=paper_mode,
        interval_name="paper_delta_u_y_phase_interval",
        change_name="paper_delta_u_y_phase_change",
        label="paper-literal delta_u_y",
    )
    measured_phase, phase_r2 = _fit_linear_quality(
        time[fit_mask], np.unwrap(np.angle(right))[fit_mask]
    )
    velocity_phase, velocity_phase_r2 = _fit_linear_quality(
        time[fit_mask], np.unwrap(np.angle(velocity_right))[fit_mask]
    )
    paper_phase, paper_phase_r2 = _fit_linear_quality(
        time[fit_mask], np.unwrap(np.angle(paper_mode))[fit_mask]
    )

    tiny = np.finfo(float).tiny
    polarization_ratio = float(
        np.min(np.abs(right[fit_mask]) / np.maximum(np.abs(left[fit_mask]), tiny))
    )
    velocity_polarization_ratio = float(
        np.min(
            np.abs(velocity_right[fit_mask])
            / np.maximum(np.abs(velocity_left[fit_mask]), tiny)
        )
    )
    expected_velocity_ratio = complex(-epsilon, -expected_growth)
    velocity_ratio_error = float(
        np.max(np.abs(velocity_right[fit_mask] / right[fit_mask] - expected_velocity_ratio))
    )

    growth_pass = (
        _within_theory_tolerance(measured_growth, expected_growth)
        and growth_r2 >= MIN_GROWTH_FIT_R2
    )
    velocity_growth_pass = (
        _within_theory_tolerance(velocity_growth, expected_growth)
        and velocity_growth_r2 >= MIN_GROWTH_FIT_R2
    )
    paper_growth_pass = (
        _within_theory_tolerance(paper_growth, expected_growth)
        and paper_growth_r2 >= MIN_GROWTH_FIT_R2
    )
    phase_pass = (
        _within_theory_tolerance(measured_phase, expected_phase)
        and _within_theory_tolerance(interval_phase, expected_phase)
        and phase_r2 >= MIN_PHASE_FIT_R2
    )
    velocity_phase_pass = (
        _within_theory_tolerance(velocity_phase, expected_phase)
        and _within_theory_tolerance(velocity_interval_phase, expected_phase)
        and velocity_phase_r2 >= MIN_PHASE_FIT_R2
    )
    paper_phase_pass = (
        _within_theory_tolerance(paper_phase, expected_phase)
        and _within_theory_tolerance(paper_interval_phase, expected_phase)
        and paper_phase_r2 >= MIN_PHASE_FIT_R2
    )
    polarization_pass = polarization_ratio >= legacy.MIN_RIGHT_TO_LEFT_RATIO
    velocity_polarization_pass = (
        velocity_polarization_ratio >= legacy.MIN_RIGHT_TO_LEFT_RATIO
    )
    velocity_ratio_pass = (
        velocity_ratio_error
        <= legacy.SOURCE_LOCAL_VELOCITY_MAGNETIC_RATIO_ABSOLUTE_TOLERANCE
    )
    return {
        "expected_growth_rate_over_k0_ua": expected_growth,
        "measured_growth_rate_over_k0_ua": measured_growth,
        "growth_fit_r2": growth_r2,
        "velocity_measured_growth_rate_over_k0_ua": velocity_growth,
        "velocity_growth_fit_r2": velocity_growth_r2,
        "paper_literal_measured_growth_rate_over_k0_ua": paper_growth,
        "paper_literal_growth_fit_r2": paper_growth_r2,
        "expected_phase_frequency_over_k0_ua": expected_phase,
        "measured_phase_frequency_over_k0_ua": measured_phase,
        "phase_fit_r2": phase_r2,
        "velocity_measured_phase_frequency_over_k0_ua": velocity_phase,
        "velocity_phase_fit_r2": velocity_phase_r2,
        "paper_literal_measured_phase_frequency_over_k0_ua": paper_phase,
        "paper_literal_phase_fit_r2": paper_phase_r2,
        "right_to_left_amplitude_ratio": polarization_ratio,
        "velocity_right_to_left_amplitude_ratio": velocity_polarization_ratio,
        "velocity_magnetic_ratio_max_absolute_error": velocity_ratio_error,
        "growth_pass": growth_pass,
        "velocity_growth_pass": velocity_growth_pass,
        "paper_literal_growth_pass": paper_growth_pass,
        "phase_pass": phase_pass,
        "velocity_phase_pass": velocity_phase_pass,
        "paper_literal_phase_pass": paper_phase_pass,
        "polarization_pass": polarization_pass,
        "velocity_polarization_pass": velocity_polarization_pass,
        "velocity_magnetic_ratio_pass": velocity_ratio_pass,
        "passed": (
            growth_pass
            and velocity_growth_pass
            and paper_growth_pass
            and phase_pass
            and velocity_phase_pass
            and paper_phase_pass
            and polarization_pass
            and velocity_polarization_pass
            and velocity_ratio_pass
        ),
    }


def synthetic_physics_trace(epsilon: float) -> dict[str, object]:
    expected_phase, expected_growth = theoretical_dispersion(epsilon)
    time = np.linspace(0.0, 6.0 / expected_growth, 49)
    right = np.exp(expected_growth * time) * np.exp(1.0j * expected_phase * time)
    left = right / 100.0
    velocity_right = (-epsilon - 1.0j * expected_growth) * right
    velocity_left = velocity_right / 100.0
    interval, change = legacy._fixed_interval_phase_trace(time, right)
    velocity_interval, velocity_change = legacy._fixed_interval_phase_trace(time, velocity_right)
    paper_interval, paper_change = legacy._fixed_interval_phase_trace(time, right)
    return {
        "normalized_time": time.tolist(),
        "right_mode_real": right.real.tolist(),
        "right_mode_imag": right.imag.tolist(),
        "left_mode_real": left.real.tolist(),
        "left_mode_imag": left.imag.tolist(),
        "velocity_right_mode_real": velocity_right.real.tolist(),
        "velocity_right_mode_imag": velocity_right.imag.tolist(),
        "velocity_left_mode_real": velocity_left.real.tolist(),
        "velocity_left_mode_imag": velocity_left.imag.tolist(),
        "phase_interval": interval.tolist(),
        "phase_change": change.tolist(),
        "velocity_phase_interval": velocity_interval.tolist(),
        "velocity_phase_change": velocity_change.tolist(),
        "paper_delta_u_y_sine_fit_real": right.real.tolist(),
        "paper_delta_u_y_sine_fit_imag": right.imag.tolist(),
        "paper_delta_u_y_phase_interval": paper_interval.tolist(),
        "paper_delta_u_y_phase_change": paper_change.tolist(),
        "paper_volume_averaged_abs_delta_u": np.abs(velocity_right).tolist(),
    }


def _manifest_members() -> dict[str, dict[str, object]]:
    return {str(member["member_id"]): member for member in validate_checked_in_decks()["cases"]}


def synthetic_provenance(
    member: Mapping[str, object], dependency: Mapping[str, object]
) -> dict[str, object]:
    deck_relative = (DECK_ROOT.relative_to(REPO_ROOT) / str(member["deck_path"])).as_posix()
    return {
        "kind": "synthetic_contract_fixture",
        "deck_path": deck_relative,
        "deck_sha256": member["deck_sha256"],
        "source_path": SOURCE_PATH.as_posix(),
        "source_sha256": _sha256_file(REPO_ROOT / SOURCE_PATH),
        "corrected_eigenmode_header_path": CORRECTED_EIGENMODE_HEADER_PATH.as_posix(),
        "corrected_eigenmode_header_sha256": _sha256_file(
            REPO_ROOT / CORRECTED_EIGENMODE_HEADER_PATH
        ),
        "q043_registered_raw_oracle_dependency_sha256": _dependency_digest(dependency),
        "authorized_artifact_root": "",
        "candidate_clean": False,
        "source_clean": False,
        "executable_clean": False,
        "executable_path": "",
        "executable_sha256": "",
        "registered_execution_receipt_path": "",
        "registered_execution_receipt_sha256": "",
        "raw_artifacts": [],
    }


def synthetic_predecessor_bundle() -> dict[str, object]:
    dependency = synthetic_q043_registered_raw_oracle_dependency()
    records = []
    for member in _manifest_members().values():
        records.append(
            {
                "member_id": member["member_id"],
                "dimension": member["dimension"],
                "epsilon": member["epsilon"],
                "resolution": member["resolution"],
                "decomposition": member["decomposition"],
                "decomposition_splits": member["decomposition_splits"],
                "physics_trace": synthetic_physics_trace(float(member["epsilon"])),
                "provenance": synthetic_provenance(member, dependency),
            }
        )
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": "q023_paper_bell_linear_joverc_predecessor_bundle",
        "campaign_id": CAMPAIGN_ID,
        "q043_registered_raw_oracle_dependency": dependency,
        "records": records,
    }


def _normalized_artifact_file(path: object, root: Path, label: str) -> Path:
    _require(isinstance(path, str) and bool(path), f"{label} must be a path")
    relative = Path(path)
    _require(not relative.is_absolute(), f"{label} must be normalized root-relative")
    try:
        resolved = (root / relative).resolve(strict=True)
        normalized = resolved.relative_to(root).as_posix()
    except (OSError, RuntimeError, ValueError) as error:
        raise ContractError(f"{label} is missing or outside the authorized root") from error
    _require(normalized == path and resolved.is_file(), f"{label} must be a normalized retained file")
    return resolved


def _validate_materialized_binding(value: object, root: Path, label: str) -> dict[str, str]:
    _require(isinstance(value, Mapping) and set(value) == _ARTIFACT_KEYS, f"{label} keys drifted")
    path = _normalized_artifact_file(value["path"], root, f"{label}/path")
    digest = value["sha256"]
    _require(isinstance(digest, str) and _SHA256.fullmatch(digest) is not None, f"{label} digest malformed")
    _require(_sha256_file(path) == digest, f"{label} digest drifted")
    return {"path": str(value["path"]), "sha256": str(digest)}


def _validate_raw_artifact_binding(
    value: object, root: Path, label: str
) -> tuple[dict[str, object], Path]:
    _require(
        isinstance(value, Mapping) and set(value) == _RAW_ARTIFACT_KEYS,
        f"{label} keys drifted",
    )
    bound = _validate_materialized_binding(
        {"path": value["path"], "sha256": value["sha256"]}, root, label
    )
    _require(value["variable"] in _RAW_OUTPUT_VARIABLES, f"{label} variable drifted")
    _require(
        type(value["cycle"]) is int and value["cycle"] >= 0,
        f"{label} cycle must be a nonnegative integer",
    )
    _require(
        type(value["time"]) in (int, float)
        and math.isfinite(float(value["time"]))
        and float(value["time"]) >= 0.0,
        f"{label} time must be finite and nonnegative",
    )
    return (
        {
            **bound,
            "variable": str(value["variable"]),
            "cycle": int(value["cycle"]),
            "time": float(value["time"]),
        },
        root / bound["path"],
    )


def _validate_rank_topology(
    value: object, *, member: Mapping[str, object], rank_count: int
) -> None:
    _require(
        isinstance(value, Mapping) and set(value) == _RANK_TOPOLOGY_KEYS,
        "Q023 MPI rank-topology evidence keys drifted",
    )
    _require(
        value["schema_version"] == SCHEMA_VERSION
        and value["record_type"]
        == "q023_paper_bell_linear_joverc_mpi_rank_topology_evidence"
        and value["campaign_id"] == CAMPAIGN_ID
        and value["member_id"] == member["member_id"]
        and value["rank_count"] == rank_count
        and value["decomposition_splits"] == member["decomposition_splits"]
        and value["publication_authorized"] is False,
        "Q023 MPI rank-topology evidence identity drifted",
    )
    ranks = value["ranks"]
    _require(
        isinstance(ranks, list) and len(ranks) == rank_count,
        "Q023 MPI rank-topology evidence rank inventory drifted",
    )
    measured_ranks: set[int] = set()
    measured_locations: set[tuple[int, int, int]] = set()
    for rank in ranks:
        _require(
            isinstance(rank, Mapping) and set(rank) == _RANK_EVIDENCE_KEYS,
            "Q023 MPI rank evidence keys drifted",
        )
        rank_id = rank["rank"]
        locations = rank["meshblock_logical_locations"]
        _require(
            type(rank_id) is int
            and 0 <= rank_id < rank_count
            and rank_id not in measured_ranks,
            "Q023 MPI rank identity drifted",
        )
        _require(
            isinstance(locations, list)
            and len(locations) == 1
            and isinstance(locations[0], list)
            and len(locations[0]) == 3
            and all(type(value) is int for value in locations[0]),
            "Q023 MPI rank must evidence exactly one MeshBlock logical location",
        )
        measured_ranks.add(rank_id)
        measured_locations.add(tuple(locations[0]))
    splits = tuple(int(value) for value in member["decomposition_splits"])
    expected_locations = {
        (x1, x2, x3)
        for x3 in range(splits[2])
        for x2 in range(splits[1])
        for x1 in range(splits[0])
    }
    _require(
        measured_ranks == set(range(rank_count))
        and measured_locations == expected_locations
        and len(measured_locations) == rank_count,
        "Q023 MPI rank topology does not evidence the prepared decomposition",
    )


def _validate_execution_receipt(
    receipt: object,
    *,
    root: Path,
    member: Mapping[str, object],
    provenance: Mapping[str, object],
    raw_artifacts: Sequence[Mapping[str, object]],
) -> None:
    _require(
        isinstance(receipt, Mapping) and set(receipt) == _EXECUTION_RECEIPT_KEYS,
        "Q023 registered execution receipt keys drifted",
    )
    expected_deck = (DECK_ROOT.relative_to(REPO_ROOT) / str(member["deck_path"])).as_posix()
    _require(
        receipt["schema_version"] == SCHEMA_VERSION
        and receipt["record_type"]
        == "q023_paper_bell_linear_joverc_registered_execution_receipt"
        and receipt["campaign_id"] == CAMPAIGN_ID
        and receipt["member_id"] == member["member_id"]
        and receipt["deck_path"] == expected_deck
        and receipt["deck_sha256"] == member["deck_sha256"]
        and receipt["source_path"] == SOURCE_PATH.as_posix()
        and receipt["source_sha256"] == provenance["source_sha256"]
        and receipt["corrected_eigenmode_header_path"]
        == CORRECTED_EIGENMODE_HEADER_PATH.as_posix()
        and receipt["corrected_eigenmode_header_sha256"]
        == provenance["corrected_eigenmode_header_sha256"]
        and receipt["executable_path"] == provenance["executable_path"]
        and receipt["executable_sha256"] == provenance["executable_sha256"]
        and receipt["candidate_clean"] is True
        and receipt["source_clean"] is True
        and receipt["executable_clean"] is True
        and receipt["publication_authorized"] is False,
        "Q023 registered execution receipt identity or clean-state binding drifted",
    )
    _require(
        receipt["raw_artifacts"] == list(raw_artifacts),
        "Q023 receipt raw inventory differs from provenance",
    )

    command = receipt["command"]
    _require(
        isinstance(command, Mapping) and set(command) == _COMMAND_KEYS,
        "Q023 execution command keys drifted",
    )
    argv = command["argv"]
    input_index = argv.index("-i") if isinstance(argv, list) and "-i" in argv else -1
    _require(
        isinstance(argv, list)
        and all(isinstance(value, str) and value for value in argv)
        and argv.count(str(receipt["executable_path"])) == 1
        and 0 <= input_index < len(argv) - 1
        and argv[input_index + 1] == expected_deck
        and isinstance(command["working_directory"], str)
        and Path(command["working_directory"]).is_absolute(),
        "Q023 registered execution command drifted",
    )

    mpi = receipt["mpi"]
    _require(
        isinstance(mpi, Mapping) and set(mpi) == _MPI_KEYS,
        "Q023 execution MPI receipt keys drifted",
    )
    rank_count = math.prod(int(value) for value in member["decomposition_splits"])
    _require(
        isinstance(mpi["launcher"], str)
        and mpi["launcher"]
        and argv[0] == mpi["launcher"]
        and mpi["rank_count"] == rank_count
        and mpi["decomposition_splits"] == member["decomposition_splits"],
        "Q023 execution MPI rank count or decomposition drifted",
    )
    rank_flag_bound = any(
        argv[index] in {"-n", "-np", "--ntasks"}
        and index + 1 < len(argv)
        and argv[index + 1] == str(rank_count)
        for index in range(len(argv))
    )
    _require(rank_flag_bound, "Q023 execution command does not bind the MPI rank count")
    topology_binding = _validate_materialized_binding(
        {
            "path": mpi["rank_topology_path"],
            "sha256": mpi["rank_topology_sha256"],
        },
        root,
        "MPI rank topology evidence",
    )
    try:
        topology = json.loads((root / topology_binding["path"]).read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ContractError("Q023 MPI rank topology evidence is not valid JSON") from error
    _validate_rank_topology(topology, member=member, rank_count=rank_count)

    termination = receipt["termination"]
    _require(
        isinstance(termination, Mapping)
        and set(termination) == _TERMINATION_KEYS
        and termination["status"] == "completed"
        and termination["return_code"] == 0
        and type(termination["final_cycle"]) is int
        and termination["final_cycle"] > 0
        and type(termination["final_time"]) in (int, float)
        and math.isfinite(float(termination["final_time"]))
        and float(termination["final_time"]) >= LINEAR_RUNTIME_TLIM,
        "Q023 registered execution termination evidence drifted",
    )
    stdout = receipt["stdout"]
    _require(
        isinstance(stdout, Mapping)
        and set(stdout) == _STDOUT_KEYS
        and stdout["completion_marker"] == "athenak_driver_completed_successfully",
        "Q023 execution stdout receipt drifted",
    )
    stdout_binding = _validate_materialized_binding(
        {"path": stdout["path"], "sha256": stdout["sha256"]}, root, "execution stdout"
    )
    try:
        stdout_text = (root / stdout_binding["path"]).read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as error:
        raise ContractError("Q023 execution stdout is unavailable") from error
    _require(
        stdout["completion_marker"] in stdout_text,
        "Q023 execution stdout completion marker is absent",
    )


def _validate_raw_runtime_parameters(
    dataset: binary.AthenaBinaryDataset, *, member: Mapping[str, object]
) -> tuple[dict[str, dict[str, str]], dict[str, tuple[int, float]]]:
    expected = parse_athinput_text(render_deck(member))
    parameters = dataset.input_parameters
    output_blocks = {block for block in parameters if block.startswith("output")}
    _require(
        output_blocks == {f"output{index}" for index in range(1, 7)},
        "Q023 raw runtime output-block inventory drifted",
    )
    normalized: dict[str, dict[str, str]] = {}
    output_state: dict[str, tuple[int, float]] = {}
    for block, expected_values in expected.items():
        _require(block in parameters, f"Q023 raw runtime block {block} is missing")
        values = parameters[block]
        for name, expected_value in expected_values.items():
            _require(
                values.get(name) == expected_value,
                f"Q023 raw runtime {block}/{name} drifted from the checked-in deck",
            )
        if block.startswith("output"):
            for key in _OUTPUT_BOOKKEEPING_KEYS:
                _require(key in values, f"Q023 raw runtime {block}/{key} is missing")
            try:
                file_number = int(values["file_number"])
                last_time = float(values["last_time"])
            except ValueError as error:
                raise ContractError(
                    f"Q023 raw runtime {block} output bookkeeping is malformed"
                ) from error
            _require(
                values["file_number"] == str(file_number)
                and file_number >= 0
                and math.isfinite(last_time)
                and last_time >= -1.0,
                f"Q023 raw runtime {block} output bookkeeping drifted",
            )
            output_state[block] = (file_number, last_time)
        normalized[block] = {
            name: value
            for name, value in values.items()
            if not (block.startswith("output") and name in _OUTPUT_BOOKKEEPING_KEYS)
        }
    return normalized, output_state


def _validate_snapshot_output_state(
    output_state: Mapping[str, tuple[int, float]],
    *,
    snapshot_index: int,
    variable: str,
) -> None:
    _require(
        set(output_state) == {f"output{index}" for index in range(1, 7)},
        "Q023 raw runtime output-state inventory drifted",
    )
    _require(variable in _RAW_OUTPUT_INDEX, "Q023 raw output variable is unknown")
    variable_output_index = _RAW_OUTPUT_INDEX[variable]
    expected = {
        f"output{index}": (
            snapshot_index if index <= variable_output_index else snapshot_index + 1,
            (
                -1.0
                if snapshot_index == 0 and index <= variable_output_index
                else (
                    snapshot_index * LINEAR_OUTPUT_DT
                    if index > variable_output_index
                    else (snapshot_index - 1) * LINEAR_OUTPUT_DT
                )
            ),
        )
        for index in range(1, 6)
    }
    expected.update(
        {
            "output6": (0, -1.0) if snapshot_index == 0 else (1, 0.0),
        }
    )
    for block, (expected_file_number, expected_last_time) in expected.items():
        file_number, last_time = output_state[block]
        _require(
            file_number == expected_file_number
            and math.isclose(
                last_time, expected_last_time, rel_tol=0.0, abs_tol=1.0e-12
            ),
            f"Q023 raw runtime {block} violates the observed sequential "
            f"{variable} publication counter progression",
        )


def _volume_mean(values: np.ndarray, grid: binary.CompositeGrid) -> float:
    volumes = (
        np.diff(grid.x3_faces)[:, None, None]
        * np.diff(grid.x2_faces)[None, :, None]
        * np.diff(grid.x1_faces)[None, None, :]
    )
    return float(np.sum(values * volumes) / np.sum(volumes))


def _representation_tolerance(variable_size: int, scale: float) -> float:
    dtype = np.float32 if variable_size == 4 else np.float64
    return DEPOSITION_TOLERANCE_ULPS * np.finfo(dtype).eps * max(abs(scale), 1.0)


def _frozen_current_relative_tolerance() -> float:
    return (
        FROZEN_CURRENT_GYROPERIOD_MARGIN
        * SPECIES_Q_OVER_MC
        * B_G
        * LINEAR_RUNTIME_TLIM
    )


def _validate_particle_current_snapshot(
    datasets: Mapping[str, binary.AthenaBinaryDataset],
    *,
    member: Mapping[str, object],
    snapshot_index: int,
) -> dict[str, object] | None:
    _require(
        set(datasets) == set(_RAW_OUTPUT_VARIABLES),
        "Q023 raw snapshot field inventory drifted",
    )
    reference = datasets["mhd_w_bcc"]
    grids: dict[str, binary.CompositeGrid] = {}
    for variable, dataset in datasets.items():
        expected_variables = _MHD_FIELDS if variable == "mhd_w_bcc" else {variable}
        _require(
            set(dataset.variable_names) == set(expected_variables),
            f"Q023 raw {variable} variable inventory drifted",
        )
        _require(
            (
                dataset.time,
                dataset.cycle,
                dataset.location_size,
                dataset.variable_size,
                dataset.root_grid_shape,
                dataset.meshblock_shape,
                dataset.domain_bounds,
            )
            == (
                reference.time,
                reference.cycle,
                reference.location_size,
                reference.variable_size,
                reference.root_grid_shape,
                reference.meshblock_shape,
                reference.domain_bounds,
            ),
            "Q023 raw cross-field runtime metadata disagrees",
        )
        if variable != "mhd_w_bcc":
            grids[variable] = binary.compose_leaf_field(dataset, variable)
    reference_grid = grids["prtcl_rho"]
    for variable in _PARTICLE_FIELDS[1:]:
        grid = grids[variable]
        _require(
            np.array_equal(grid.x1_faces, reference_grid.x1_faces)
            and np.array_equal(grid.x2_faces, reference_grid.x2_faces)
            and np.array_equal(grid.x3_faces, reference_grid.x3_faces),
            "Q023 raw cross-field particle geometry disagrees",
        )

    particle_values = [grids[variable].values for variable in _PARTICLE_FIELDS]
    if snapshot_index == 0:
        _require(
            all(np.count_nonzero(values) == 0 for values in particle_values),
            "Q023 cycle-zero particle moments must be the pre-deposition zero state",
        )
        return None

    rho = grids["prtcl_rho"].values
    current = np.stack(
        [grids[variable].values for variable in _PARTICLE_FIELDS[1:]], axis=0
    )
    basis = np.asarray(_mode_basis(int(member["dimension"])), dtype=np.float64)
    expected_vector = EXPECTED_J_OVER_C * basis
    expected_rho = (
        int(member["ppc"])
        * float(member["deposit_qscale"])
        * float(member["species_charge"])
        / float(member["root_cell_volume"])
    )
    _require(
        math.isclose(
            expected_rho * float(member["stream_speed"]),
            EXPECTED_J_OVER_C,
            rel_tol=1.0e-13,
            abs_tol=1.0e-13,
        ),
        "Q023 configured deposited charge/current closure drifted",
    )
    mean_vector = np.asarray(
        [
            _volume_mean(current[index], grids[_PARTICLE_FIELDS[index + 1]])
            for index in range(3)
        ]
    )
    parallel = np.tensordot(basis, current, axes=1)
    transverse = current - basis[:, None, None, None] * parallel[None, ...]
    residual = current - expected_vector[:, None, None, None]
    parallel_mean = float(np.dot(basis, mean_vector))
    rho_mean = _volume_mean(rho, grids["prtcl_rho"])
    max_transverse = float(np.max(np.sqrt(np.sum(transverse * transverse, axis=0))))
    max_current_nonuniformity = float(np.max(np.sqrt(np.sum(residual * residual, axis=0))))
    max_rho_nonuniformity = float(np.max(np.abs(rho - expected_rho)))
    current_tolerance = max(
        _representation_tolerance(reference.variable_size, EXPECTED_J_OVER_C),
        _frozen_current_relative_tolerance() * EXPECTED_J_OVER_C,
    )
    rho_tolerance = max(
        _representation_tolerance(reference.variable_size, expected_rho),
        _frozen_current_relative_tolerance() * expected_rho,
    )
    _require(
        abs(parallel_mean - EXPECTED_J_OVER_C) <= current_tolerance,
        "Q023 guide-projected deposited current is not effectively frozen",
    )
    _require(
        max_transverse <= current_tolerance,
        "Q023 transverse deposited current is not effectively frozen",
    )
    _require(
        max_current_nonuniformity <= current_tolerance,
        "Q023 deposited current spatial nonuniformity exceeds the frozen-current bound",
    )
    _require(
        abs(rho_mean - expected_rho) <= rho_tolerance
        and max_rho_nonuniformity <= rho_tolerance,
        "Q023 deposited charge density violates the volume-aware closure",
    )
    return {
        "cycle": reference.cycle,
        "time": reference.time,
        "mean_vector": mean_vector,
        "rho_mean": rho_mean,
        "current_tolerance": current_tolerance,
        "rho_tolerance": rho_tolerance,
    }


def _validate_frozen_current_trace(
    reports: Sequence[Mapping[str, object]],
) -> None:
    _require(bool(reports), "Q023 raw current trace requires post-step measurements")
    baseline_vector = np.asarray(reports[0]["mean_vector"], dtype=float)
    baseline_rho = float(reports[0]["rho_mean"])
    current_tolerance = float(reports[0]["current_tolerance"])
    rho_tolerance = float(reports[0]["rho_tolerance"])
    for report in reports[1:]:
        _require(
            np.linalg.norm(np.asarray(report["mean_vector"], dtype=float) - baseline_vector)
            <= current_tolerance,
            "Q023 deposited-current drift from the first post-step state exceeds the frozen-current bound",
        )
        _require(
            abs(float(report["rho_mean"]) - baseline_rho) <= rho_tolerance,
            "Q023 deposited-charge drift from the first post-step state exceeds the frozen-current bound",
        )


def _dataset_observable_mapping(
    dataset: binary.AthenaBinaryDataset,
) -> dict[str, object]:
    _require(
        set(dataset.variable_names) == _MHD_FIELDS,
        "Q023 raw mhd_w_bcc variable inventory drifted",
    )
    composites = {
        name: binary.compose_leaf_field(dataset, name) for name in _MHD_FIELDS
    }
    reference = composites["dens"]

    def centers(faces: np.ndarray) -> np.ndarray:
        return 0.5 * (faces[:-1] + faces[1:])

    return {
        "Time": dataset.time,
        "x1v": centers(reference.x1_faces),
        "x2v": centers(reference.x2_faces),
        "x3v": centers(reference.x3_faces),
        **{name: composite.values for name, composite in composites.items()},
    }


def _physics_trace_from_raw_datasets(
    datasets: Sequence[binary.AthenaBinaryDataset],
    *,
    member: Mapping[str, object],
) -> dict[str, object]:
    geometry = {
        "nx": list(member["global_nx"]),
        "xmin": [float(axis[0]) for axis in member["bounds"]],
        "extent": [float(axis[1]) - float(axis[0]) for axis in member["bounds"]],
    }
    rows = []
    for dataset in datasets:
        mapped = _dataset_observable_mapping(dataset)
        legacy_right, legacy_left, legacy_velocity_right, legacy_velocity_left = (
            legacy._spatial_modes_from_dataset(
                mapped, int(member["dimension"]), geometry
            )
        )
        # The historical extractor named the Bai et al. right-handed unstable
        # polarization "left". Keep historical code byte-preserved and bind the
        # corrected predecessor to the physical polarization here.
        right, left = legacy_left, legacy_right
        velocity_right, velocity_left = (
            legacy_velocity_left,
            legacy_velocity_right,
        )
        paper_mode, paper_abs = legacy._paper_literal_velocity_observables_from_dataset(
            mapped, int(member["dimension"]), geometry
        )
        rows.append(
            (
                float(dataset.time) * K0,
                right,
                left,
                velocity_right,
                velocity_left,
                paper_mode,
                paper_abs,
            )
        )
    rows.sort(key=lambda row: row[0])
    _require(bool(rows), "Q023 raw trace extraction requires snapshots")
    time = np.asarray([row[0] for row in rows], dtype=float)
    right = np.asarray([row[1] for row in rows], dtype=complex)
    left = np.asarray([row[2] for row in rows], dtype=complex)
    velocity_right = np.asarray([row[3] for row in rows], dtype=complex)
    velocity_left = np.asarray([row[4] for row in rows], dtype=complex)
    paper_mode = np.asarray([row[5] for row in rows], dtype=complex)
    paper_abs = np.asarray([row[6] for row in rows], dtype=float)
    phase_interval, phase_change = legacy._fixed_interval_phase_trace(time, right)
    velocity_interval, velocity_change = legacy._fixed_interval_phase_trace(
        time, velocity_right
    )
    paper_interval, paper_change = legacy._fixed_interval_phase_trace(time, paper_mode)
    return {
        "normalized_time": time.tolist(),
        "right_mode_real": right.real.tolist(),
        "right_mode_imag": right.imag.tolist(),
        "left_mode_real": left.real.tolist(),
        "left_mode_imag": left.imag.tolist(),
        "velocity_right_mode_real": velocity_right.real.tolist(),
        "velocity_right_mode_imag": velocity_right.imag.tolist(),
        "velocity_left_mode_real": velocity_left.real.tolist(),
        "velocity_left_mode_imag": velocity_left.imag.tolist(),
        "phase_interval": phase_interval.tolist(),
        "phase_change": phase_change.tolist(),
        "velocity_phase_interval": velocity_interval.tolist(),
        "velocity_phase_change": velocity_change.tolist(),
        "paper_delta_u_y_sine_fit_real": paper_mode.real.tolist(),
        "paper_delta_u_y_sine_fit_imag": paper_mode.imag.tolist(),
        "paper_delta_u_y_phase_interval": paper_interval.tolist(),
        "paper_delta_u_y_phase_change": paper_change.tolist(),
        "paper_volume_averaged_abs_delta_u": paper_abs.tolist(),
    }


def _validate_provenance(
    value: object,
    *,
    member: Mapping[str, object],
    dependency: Mapping[str, object],
    artifact_root: Path | None,
    physics_trace: object,
) -> str:
    _require(isinstance(value, Mapping) and set(value) == _PROVENANCE_KEYS, "Q023 trace provenance keys drifted")
    provenance = dict(value)
    expected_deck = (DECK_ROOT.relative_to(REPO_ROOT) / str(member["deck_path"])).as_posix()
    _require(provenance["deck_path"] == expected_deck and provenance["deck_sha256"] == member["deck_sha256"], "Q023 trace deck binding drifted")
    _require(provenance["source_path"] == SOURCE_PATH.as_posix() and provenance["source_sha256"] == _sha256_file(REPO_ROOT / SOURCE_PATH), "Q023 trace source binding drifted")
    _require(
        provenance["corrected_eigenmode_header_path"]
        == CORRECTED_EIGENMODE_HEADER_PATH.as_posix()
        and provenance["corrected_eigenmode_header_sha256"]
        == _sha256_file(REPO_ROOT / CORRECTED_EIGENMODE_HEADER_PATH),
        "Q023 corrected eigenmode header binding drifted",
    )
    _require(provenance["q043_registered_raw_oracle_dependency_sha256"] == _dependency_digest(dependency), "Q043 dependency digest binding drifted")
    if provenance["kind"] == "synthetic_contract_fixture":
        _require(provenance == synthetic_provenance(member, dependency), "synthetic Q023 provenance drifted")
        return "synthetic_contract_fixture"
    _require(provenance["kind"] == "registered_execution_trace", "Q023 trace provenance kind drifted")
    _require(
        dependency["registered_admission_digest_bound"] is True
        and dependency["registered_admission_schema_bound"] is True
        and dependency["registered_execution_qualification_check_pass"] is True
        and dependency["registered_raw_oracle_pass"] is True
        and dependency["complete_foundational_raw_oracle_matrix_pass"] is True,
        "materialized Q023 trace requires the final hardened Q043 registered admission",
    )
    _require(artifact_root is not None and Path(artifact_root).is_absolute(), "materialized Q023 provenance requires absolute artifact root")
    try:
        root = Path(artifact_root).resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise ContractError("materialized Q023 artifact root is missing") from error
    _require(root.is_dir() and provenance["authorized_artifact_root"] == str(root), "Q023 authorized artifact root drifted")
    _require(provenance["candidate_clean"] is True and provenance["source_clean"] is True and provenance["executable_clean"] is True, "Q023 materialized candidate/source/executable must be clean")
    for prefix in ("executable", "registered_execution_receipt"):
        _validate_materialized_binding(
            {"path": provenance[f"{prefix}_path"], "sha256": provenance[f"{prefix}_sha256"]},
            root,
            prefix,
        )
    artifacts = provenance["raw_artifacts"]
    _require(isinstance(artifacts, list) and bool(artifacts), "Q023 materialized trace requires raw artifacts")
    normalized_and_paths = [
        _validate_raw_artifact_binding(item, root, "raw_artifact") for item in artifacts
    ]
    normalized = [item[0] for item in normalized_and_paths]
    _require(len({item["path"] for item in normalized}) == len(normalized), "Q023 raw artifact reused")
    raw_sort_key = lambda item: (
        item["cycle"],
        item["time"],
        _RAW_OUTPUT_INDEX[item["variable"]],
        item["path"],
    )
    _require(
        normalized == sorted(normalized, key=raw_sort_key)
        and normalized[0]["cycle"] == 0
        and normalized[0]["time"] == 0.0,
        "Q023 raw artifact inventory must begin at cycle zero and be ordered",
    )
    snapshot_groups: list[
        tuple[tuple[int, float], list[tuple[dict[str, object], Path]]]
    ] = []
    for artifact_and_path in normalized_and_paths:
        artifact = artifact_and_path[0]
        key = (int(artifact["cycle"]), float(artifact["time"]))
        if not snapshot_groups or snapshot_groups[-1][0] != key:
            snapshot_groups.append((key, []))
        snapshot_groups[-1][1].append(artifact_and_path)
    _require(
        all(
            len(group) == len(_RAW_OUTPUT_VARIABLES)
            and {str(item[0]["variable"]) for item in group}
            == set(_RAW_OUTPUT_VARIABLES)
            for _, group in snapshot_groups
        ),
        "Q023 raw snapshot field inventory drifted",
    )

    mhd_datasets = []
    current_reports = []
    normalized_parameters = None
    previous_cycle = -1
    previous_time = -1.0
    for snapshot_index, ((cycle, time), group) in enumerate(snapshot_groups):
        _require(
            cycle > previous_cycle and time > previous_time,
            "Q023 raw artifact metadata or ordering drifted",
        )
        snapshot_datasets = {}
        for artifact, path in group:
            variable = str(artifact["variable"])
            try:
                dataset = binary.parse_athenak_binary_bytes(
                    path.read_bytes(), source=str(path)
                )
            except (OSError, binary.AnalysisError) as error:
                raise ContractError("Q023 raw AthenaK output failed strict parsing") from error
            _require(
                dataset.cycle == cycle and dataset.time == time,
                "Q023 raw cross-field artifact metadata disagrees",
            )
            _require(
                dataset.root_grid_shape == tuple(member["global_nx"])
                and dataset.meshblock_shape == tuple(member["meshblock_nx"]),
                "Q023 raw grid or decomposition geometry drifted",
            )
            runtime_parameters, output_state = _validate_raw_runtime_parameters(
                dataset, member=member
            )
            _validate_snapshot_output_state(
                output_state, snapshot_index=snapshot_index, variable=variable
            )
            if normalized_parameters is None:
                normalized_parameters = runtime_parameters
            else:
                _require(
                    runtime_parameters == normalized_parameters,
                    "Q023 raw immutable runtime parameters differ across snapshots and fields",
                )
            snapshot_datasets[variable] = dataset
        try:
            current_report = _validate_particle_current_snapshot(
                snapshot_datasets, member=member, snapshot_index=snapshot_index
            )
        except binary.AnalysisError as error:
            raise ContractError(
                "Q023 raw particle-current fields failed strict composition"
            ) from error
        if current_report is not None:
            current_reports.append(current_report)
        previous_cycle = cycle
        previous_time = time
        mhd_datasets.append(snapshot_datasets["mhd_w_bcc"])
    _validate_frozen_current_trace(current_reports)
    _, expected_growth = theoretical_dispersion(float(member["epsilon"]))
    _require(
        mhd_datasets[-1].time * K0 >= 5.0 / expected_growth,
        "Q023 raw trace does not cover the fixed growth-fit window",
    )
    receipt_path = _normalized_artifact_file(
        provenance["registered_execution_receipt_path"],
        root,
        "registered_execution_receipt/path",
    )
    try:
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ContractError("Q023 registered execution receipt is not valid JSON") from error
    _validate_execution_receipt(
        receipt,
        root=root,
        member=member,
        provenance=provenance,
        raw_artifacts=normalized,
    )
    _require(
        isinstance(physics_trace, Mapping)
        and _canonical_json_bytes(physics_trace)
        == _canonical_json_bytes(
            _physics_trace_from_raw_datasets(mhd_datasets, member=member)
        ),
        "Q023 materialized physics trace was not derived from the exact retained raw outputs",
    )
    return "registered_execution_trace"


def _analyze_matrix_record(
    record: Mapping[str, object],
    *,
    members: Mapping[str, Mapping[str, object]],
    dependency: Mapping[str, object],
    artifact_root: Path | None,
) -> dict[str, object]:
    _require(set(record) == _MATRIX_RECORD_KEYS, "Q023 matrix record keys drifted")
    member_id = str(record["member_id"])
    _require(member_id in members, "Q023 matrix member is unknown")
    member = members[member_id]
    for name in ("dimension", "resolution", "decomposition", "decomposition_splits"):
        _require(record[name] == member[name], f"Q023 matrix {name} drifted")
    epsilon = record["epsilon"]
    _require(
        type(epsilon) in (int, float)
        and float(epsilon) == float(member["epsilon"])
        and float(epsilon) in EPSILON_VALUES,
        "Q023 epsilon drifted",
    )
    provenance_kind = _validate_provenance(
        record["provenance"],
        member=member,
        dependency=dependency,
        artifact_root=artifact_root,
        physics_trace=record["physics_trace"],
    )
    trace = record["physics_trace"]
    _require(isinstance(trace, Mapping), "Q023 physics trace must be an object")
    physics = _analyze_physics_trace(trace, float(epsilon))
    return {
        "member_id": member_id,
        "dimension": member["dimension"],
        "epsilon": float(epsilon),
        "resolution": member["resolution"],
        "decomposition": member["decomposition"],
        "decomposition_splits": member["decomposition_splits"],
        "provenance_kind": provenance_kind,
        "growth_gate_pass": physics["growth_pass"] and physics["velocity_growth_pass"] and physics["paper_literal_growth_pass"],
        "signed_phase_gate_pass": physics["phase_pass"] and physics["velocity_phase_pass"] and physics["paper_literal_phase_pass"],
        "polarization_gate_pass": physics["polarization_pass"] and physics["velocity_polarization_pass"],
        "growth_fit_quality_gate_pass": (
            physics["growth_fit_r2"] >= MIN_GROWTH_FIT_R2
            and physics["velocity_growth_fit_r2"] >= MIN_GROWTH_FIT_R2
            and physics["paper_literal_growth_fit_r2"] >= MIN_GROWTH_FIT_R2
        ),
        "phase_fit_quality_gate_pass": (
            physics["phase_fit_r2"] >= MIN_PHASE_FIT_R2
            and physics["velocity_phase_fit_r2"] >= MIN_PHASE_FIT_R2
            and physics["paper_literal_phase_fit_r2"] >= MIN_PHASE_FIT_R2
        ),
        "velocity_magnetic_ratio_gate_pass": physics["velocity_magnetic_ratio_pass"],
        "provenance_gate_pass": True,
        "measured_growth_rate_over_k0_ua": physics["measured_growth_rate_over_k0_ua"],
        "expected_growth_rate_over_k0_ua": physics["expected_growth_rate_over_k0_ua"],
        "measured_signed_phase_frequency_over_k0_ua": physics["measured_phase_frequency_over_k0_ua"],
        "expected_signed_phase_frequency_over_k0_ua": physics["expected_phase_frequency_over_k0_ua"],
        "right_to_left_amplitude_ratio": physics["right_to_left_amplitude_ratio"],
        "growth_fit_r2": physics["growth_fit_r2"],
        "velocity_growth_fit_r2": physics["velocity_growth_fit_r2"],
        "paper_literal_growth_fit_r2": physics["paper_literal_growth_fit_r2"],
        "phase_fit_r2": physics["phase_fit_r2"],
        "velocity_phase_fit_r2": physics["velocity_phase_fit_r2"],
        "paper_literal_phase_fit_r2": physics["paper_literal_phase_fit_r2"],
        "physics_gates_pass": physics["passed"],
    }


def _convergence_and_decomposition_reports(
    reports: Sequence[Mapping[str, object]]
) -> tuple[list[dict[str, object]], list[dict[str, object]], bool]:
    by_key = {
        (int(report["dimension"]), float(report["epsilon"]), str(report["resolution"]), str(report["decomposition"])): report
        for report in reports
    }
    convergence = []
    decomposition = []
    all_pass = True
    for dimension in DIMENSIONS:
        for epsilon in EPSILON_VALUES:
            coarse = by_key[(dimension, epsilon, "coarse", "reference_x1")]
            fine = by_key[(dimension, epsilon, "fine", "reference_x1")]
            growth_difference = abs(float(fine["measured_growth_rate_over_k0_ua"]) - float(coarse["measured_growth_rate_over_k0_ua"]))
            phase_difference = abs(float(fine["measured_signed_phase_frequency_over_k0_ua"]) - float(coarse["measured_signed_phase_frequency_over_k0_ua"]))
            fine_growth_error = abs(float(fine["measured_growth_rate_over_k0_ua"]) - float(fine["expected_growth_rate_over_k0_ua"]))
            coarse_growth_error = abs(float(coarse["measured_growth_rate_over_k0_ua"]) - float(coarse["expected_growth_rate_over_k0_ua"]))
            fine_phase_error = abs(float(fine["measured_signed_phase_frequency_over_k0_ua"]) - float(fine["expected_signed_phase_frequency_over_k0_ua"]))
            coarse_phase_error = abs(float(coarse["measured_signed_phase_frequency_over_k0_ua"]) - float(coarse["expected_signed_phase_frequency_over_k0_ua"]))
            passed = (
                growth_difference <= GROWTH_CONVERGENCE_ABSOLUTE_TOLERANCE
                and phase_difference <= PHASE_CONVERGENCE_ABSOLUTE_TOLERANCE
                and fine_growth_error <= coarse_growth_error + 1.0e-12
                and fine_phase_error <= coarse_phase_error + 1.0e-12
            )
            convergence.append(
                {
                    "dimension": dimension,
                    "epsilon": epsilon,
                    "growth_difference": growth_difference,
                    "signed_phase_difference": phase_difference,
                    "fine_not_worse_than_coarse": fine_growth_error <= coarse_growth_error + 1.0e-12 and fine_phase_error <= coarse_phase_error + 1.0e-12,
                    "convergence_gate_pass": passed,
                }
            )
            all_pass = all_pass and passed
            for name in DECOMPOSITIONS_BY_DIMENSION[dimension]:
                if name == "reference_x1":
                    continue
                split = by_key[(dimension, epsilon, "fine", name)]
                growth_delta = abs(float(split["measured_growth_rate_over_k0_ua"]) - float(fine["measured_growth_rate_over_k0_ua"]))
                phase_delta = abs(float(split["measured_signed_phase_frequency_over_k0_ua"]) - float(fine["measured_signed_phase_frequency_over_k0_ua"]))
                polarization_relative_delta = abs(float(split["right_to_left_amplitude_ratio"]) / float(fine["right_to_left_amplitude_ratio"]) - 1.0)
                passed = (
                    growth_delta <= DECOMPOSITION_GROWTH_ABSOLUTE_TOLERANCE
                    and phase_delta <= DECOMPOSITION_PHASE_ABSOLUTE_TOLERANCE
                    and polarization_relative_delta <= DECOMPOSITION_POLARIZATION_RELATIVE_TOLERANCE
                )
                decomposition.append(
                    {
                        "dimension": dimension,
                        "epsilon": epsilon,
                        "decomposition": name,
                        "decomposition_splits": list(DECOMPOSITIONS_BY_DIMENSION[dimension][name]),
                        "growth_delta_from_reference": growth_delta,
                        "signed_phase_delta_from_reference": phase_delta,
                        "polarization_relative_delta_from_reference": polarization_relative_delta,
                        "decomposition_invariance_gate_pass": passed,
                    }
                )
                all_pass = all_pass and passed
    return convergence, decomposition, all_pass


def analyze_predecessor_bundle(
    bundle: Mapping[str, object], *, artifact_root: Path | None = None
) -> dict[str, object]:
    """Analyze the complete corrected linear predecessor without granting authority."""
    _require(
        set(bundle)
        == {
            "schema_version",
            "record_type",
            "campaign_id",
            "q043_registered_raw_oracle_dependency",
            "records",
        },
        "Q023 predecessor bundle keys drifted",
    )
    _require(bundle["schema_version"] == SCHEMA_VERSION, "Q023 predecessor schema drifted")
    _require(bundle["record_type"] == "q023_paper_bell_linear_joverc_predecessor_bundle", "Q023 predecessor type drifted")
    _require(bundle["campaign_id"] == CAMPAIGN_ID, "Q023 predecessor campaign drifted")
    dependency = validate_q043_dependency(
        bundle["q043_registered_raw_oracle_dependency"], artifact_root=artifact_root
    )
    members = _manifest_members()
    records = bundle["records"]
    _require(isinstance(records, list), "Q023 predecessor records must be a list")
    expected = {
        (member_id, float(member["epsilon"]))
        for member_id, member in members.items()
    }
    measured = []
    reports = []
    for record in records:
        _require(isinstance(record, Mapping), "Q023 predecessor record must be an object")
        report = _analyze_matrix_record(
            record,
            members=members,
            dependency=dependency,
            artifact_root=artifact_root,
        )
        key = (str(report["member_id"]), float(report["epsilon"]))
        _require(key not in measured, "duplicate Q023 predecessor matrix record")
        measured.append(key)
        reports.append(report)
    _require(set(measured) == expected, "Q023 predecessor matrix is incomplete or contains extra records")
    provenance_kinds = {str(report["provenance_kind"]) for report in reports}
    _require(len(provenance_kinds) == 1, "Q023 predecessor matrix cannot mix provenance kinds")
    provenance_kind = next(iter(provenance_kinds))
    registered_q043_foundation_bound = (
        dependency["registered_admission_digest_bound"] is True
        and dependency["registered_admission_schema_bound"] is True
        and dependency["registered_execution_qualification_check_pass"] is True
        and dependency["registered_raw_oracle_pass"] is True
        and dependency["complete_foundational_raw_oracle_matrix_pass"] is True
    )
    _require(
        provenance_kind == "synthetic_contract_fixture"
        or registered_q043_foundation_bound,
        "materialized Q023 predecessor evidence is not downstream of the final hardened Q043 admission",
    )
    reports.sort(key=lambda item: (item["dimension"], item["epsilon"], item["member_id"]))
    convergence, decomposition, comparison_pass = _convergence_and_decomposition_reports(reports)
    physics_pass = all(bool(report["physics_gates_pass"]) for report in reports)
    predecessor_contract_pass = physics_pass and comparison_pass
    materialized_q023_trace_matrix_bound = provenance_kind == "registered_execution_trace"
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": "q023_paper_bell_linear_joverc_predecessor_analysis",
        "campaign_id": CAMPAIGN_ID,
        "foundational_current_campaign_id": Q043_CURRENT_CAMPAIGN_ID,
        "foundational_current_lineage_commit": dependency[
            "integration_checkpoint_commit"
        ],
        "foundational_raw_oracle_id": Q043_RAW_ORACLE_ID,
        "q043_dependency_contract_pass": True,
        "q043_registered_raw_oracle_dependency_pass": registered_q043_foundation_bound,
        "registered_q043_foundation_bound": registered_q043_foundation_bound,
        "materialized_q023_trace_matrix_bound": materialized_q023_trace_matrix_bound,
        "analysis_input_kind": provenance_kind,
        "record_count": len(reports),
        "growth_gates_pass": all(bool(report["growth_gate_pass"]) for report in reports),
        "signed_phase_gates_pass": all(bool(report["signed_phase_gate_pass"]) for report in reports),
        "polarization_gates_pass": all(bool(report["polarization_gate_pass"]) for report in reports),
        "growth_fit_quality_gates_pass": all(
            bool(report["growth_fit_quality_gate_pass"]) for report in reports
        ),
        "phase_fit_quality_gates_pass": all(
            bool(report["phase_fit_quality_gate_pass"]) for report in reports
        ),
        "provenance_gates_pass": all(bool(report["provenance_gate_pass"]) for report in reports),
        "resolution_convergence_gate_pass": all(item["convergence_gate_pass"] for item in convergence),
        "multidirectional_decomposition_gate_pass": all(item["decomposition_invariance_gate_pass"] for item in decomposition),
        "predecessor_contract_pass": predecessor_contract_pass,
        "linear_qualification_prerequisites_pass": (
            predecessor_contract_pass
            and registered_q043_foundation_bound
            and materialized_q023_trace_matrix_bound
        ),
        "qualification_effect": QUALIFICATION_EFFECT,
        "launch_authorized": False,
        "qualification_eligible": False,
        "scientific_claim_authorized": False,
        "publication_authorized": False,
        "complete_authoritative_q043_foundation_claimed": False,
        "passed": False,
        "records": reports,
        "convergence": convergence,
        "decomposition_invariance": decomposition,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--materialize-checked-in-decks", action="store_true")
    parser.add_argument("--validate-checked-in-decks", action="store_true")
    parser.add_argument("--replace", action="store_true")
    parser.add_argument("--synthetic-bundle", action="store_true")
    parser.add_argument("--artifact-root", type=Path)
    parser.add_argument("bundle", type=Path, nargs="?")
    args = parser.parse_args()
    selected = sum(
        (
            args.materialize_checked_in_decks,
            args.validate_checked_in_decks,
            args.synthetic_bundle,
            args.bundle is not None,
        )
    )
    _require(selected == 1, "select exactly one Q023 predecessor operation")
    if args.materialize_checked_in_decks:
        result = materialize_checked_in_decks(replace=args.replace)
    elif args.validate_checked_in_decks:
        result = validate_checked_in_decks()
    elif args.synthetic_bundle:
        result = synthetic_predecessor_bundle()
    else:
        result = analyze_predecessor_bundle(
            json.loads(args.bundle.read_text(encoding="utf-8")),
            artifact_root=args.artifact_root,
        )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
