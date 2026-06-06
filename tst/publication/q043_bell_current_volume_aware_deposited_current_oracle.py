#!/usr/bin/env python3
"""Source-local raw-output oracle for the corrected Q-043 Bell current."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import shutil
from typing import Any, Mapping, Sequence

import numpy as np

from tst.publication import analyze_q011_section54_outputs as binary


REPO_ROOT = Path(__file__).resolve().parents[2]
CHECKED_IN_DECK_ROOT = (
    REPO_ROOT / "inputs/tests/q043_bell_current_volume_aware_deposited_current_oracle"
)
CHECKED_IN_MANIFEST = CHECKED_IN_DECK_ROOT / "deck_manifest.json"
CAMPAIGN_ID = "Q043-BELL-DEPOSITED-J-OVER-C-VOLUME-AWARE"
PGEN_CAMPAIGN_ID = "Q043-BELL-CURRENT-VOLUME-AWARE"
PGEN_NAME = "q043_bell_current_volume_aware"
SCHEMA_VERSION = 1
FIELDS = ("prtcl_rho", "prtcl_jx", "prtcl_jy", "prtcl_jz")
DIMENSIONS = (1, 2, 3)
RESOLUTIONS = ("coarse", "fine")
PPC_VALUES = (1, 4)
DECOMPOSITIONS_BY_DIMENSION = {
    1: ("single", "split_x1"),
    2: ("single", "split_x1", "split_x2", "split_x1x2"),
    3: ("single", "split_x1", "split_x2", "split_x3", "split_xyz"),
}
DECOMPOSITIONS = tuple(
    dict.fromkeys(
        decomposition
        for dimension in DIMENSIONS
        for decomposition in DECOMPOSITIONS_BY_DIMENSION[dimension]
    )
)
_DECOMPOSITION_AXES = {
    "single": (),
    "split_x1": (0,),
    "split_x2": (1,),
    "split_x1x2": (0, 1),
    "split_x3": (2,),
    "split_xyz": (0, 1, 2),
}
ARTIFICIAL_C_OVER_V_CR_VALUES = (100, 1000, 10000)
EXPECTED_CASE_COUNT = (
    sum(len(DECOMPOSITIONS_BY_DIMENSION[dimension]) for dimension in DIMENSIONS)
    * len(RESOLUTIONS)
    * len(PPC_VALUES)
    * len(ARTIFICIAL_C_OVER_V_CR_VALUES)
)
SPECIES_MASS = 2.0
SPECIES_CHARGE = 4.0 * math.pi * 1.0e-6
SPECIES_Q_OVER_MC = SPECIES_CHARGE / SPECIES_MASS
STREAM_SPEED = 2.5
B_G = 1.0
K0 = 2.0 * math.pi
EXPECTED_J_OVER_C = 2.0 * B_G * K0
EXPECTED_RHO = EXPECTED_J_OVER_C / STREAM_SPEED
RAW_ORACLE_CYCLE = 1
RAW_ORACLE_DCYCLE = 2
RUNTIME_OUTPUT_BOOKKEEPING_KEYS = frozenset(("file_number", "last_time"))
RUNTIME_OUTPUT_BLOCKS = tuple(f"output{index}" for index in range(1, len(FIELDS) + 1))
RUNTIME_METADATA_CONTRACT = (
    "strict_after_validating_observed_cycle_one_sequential_output_bookkeeping_"
    "and_normalizing_only_file_number_and_last_time"
)
_RUNTIME_DEFAULTS = {
    "coord": {"special_rel": "0", "general_rel": "0"},
    "mesh_refinement": {"refinement": "none"},
    "mhd": {
        "dfloor": "1.17549e-38",
        "pfloor": "1.17549e-38",
        "tfloor": "1.17549e-38",
        "sfloor": "1.17549e-38",
        "nscalars": "0",
        "const_accel": "0",
        "cooling_dt_factor": "1",
        "t_start_ism_cooling": "0",
        "ism_cooling": "0",
        "cgm_cooling": "0",
        "beam_source": "0",
        "rel_cooling": "0",
        "fofc": "0",
    },
    "particles": {
        "track_displacement": "0",
        "pic_ion_neutral_collision_rate": "0",
        "pic_deltaf_f0": "",
        "pic_deltaf_p0": "1",
        "pic_deltaf_kappa": "1.25",
        "pic_deltaf_drift_x1": "0",
        "pic_deltaf_drift_x2": "0",
        "pic_deltaf_drift_x3": "0",
        "pic_deltaf_aniso_x1": "1",
        "pic_deltaf_aniso_x2": "1",
        "pic_deltaf_aniso_x3": "1",
        "pic_deltaf_background_rho": "0",
        "pic_deltaf_background_jx": "0",
        "pic_deltaf_background_jy": "0",
        "pic_deltaf_background_jz": "0",
        "pic_deltaf_adapt_mode": "off",
        "pic_deltaf_adapt_interval": "0",
        "pic_q017_sync_kernel_timers": "0",
        "pic_random_seed": "0",
        "pic_expansion_law": "linear",
        "pic_expansion_rate_x1": "0",
        "pic_expansion_rate_x2": "0",
        "pic_expansion_rate_x3": "0",
        "pic_no_mhd_bx": "0",
        "pic_no_mhd_by": "0",
        "pic_no_mhd_bz": "0",
        "assign_tag": "index_order",
    },
    "problem": {"user_srcs": "0", "user_hist": "0", "user_work_in_loop": "0"},
    "time": {"start_time": "0"},
    "par_end": {},
}
DEPOSITION_TOLERANCE_ULPS = 128.0
QUALIFICATION_EFFECT = "none_source_local_output_oracle_only"

_GEOMETRY = {
    1: {
        "bounds": ((0.0, 1.0), (0.0, 1.0), (0.0, 1.0)),
        "coarse": (16, 4, 1),
        "fine": (32, 4, 1),
    },
    2: {
        "bounds": ((0.0, math.sqrt(5.0)), (0.0, math.sqrt(1.25)), (0.0, 1.0)),
        "coarse": (16, 8, 1),
        "fine": (32, 16, 1),
    },
    3: {
        "bounds": (
            (0.0, math.sqrt(21.0)),
            (0.0, math.sqrt(5.25)),
            (0.0, math.sqrt(1.3125)),
        ),
        "coarse": (16, 8, 8),
        "fine": (32, 16, 16),
    },
}


class ContractError(ValueError):
    """Raised when the source-local deposited-current oracle fails closed."""


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
    """Return qscale for volume-averaged deposited J_CR/c = 2 B_g k0."""
    _require(
        root_cell_volume > 0.0
        and ppc > 0
        and species_charge > 0.0
        and stream_speed > 0.0,
        "deposited-current normalization factors must be positive",
    )
    return (
        EXPECTED_J_OVER_C
        * root_cell_volume
        / (ppc * species_charge * stream_speed)
    )


def configured_volume_mean_j_over_c(case: Mapping[str, object]) -> float:
    """Reconstruct the configured volume-mean deposited current."""
    return (
        int(case["ppc"])
        * float(case["deposit_qscale"])
        * float(case["species_charge"])
        * STREAM_SPEED
        / float(case["root_cell_volume"])
    )


def _case_id(
    dimension: int,
    resolution: str,
    ppc: int,
    decomposition: str,
    artificial_c_over_v_cr: int,
) -> str:
    return (
        f"q043-current-oracle-d{dimension}-{resolution}-"
        f"ppc{ppc}-{decomposition}-cvr{artificial_c_over_v_cr}"
    )


def _decomposition(
    dimension: int, nx: Sequence[int], name: str
) -> tuple[tuple[int, int, int], tuple[int, int, int], int]:
    """Return a dimension-valid MeshBlock grid, MeshBlock extent, and rank count."""
    _require(
        name in DECOMPOSITIONS_BY_DIMENSION[dimension],
        f"d{dimension}: decomposition is not allowed",
    )
    partitioned_axes = _DECOMPOSITION_AXES[name]
    meshblock_grid = tuple(2 if axis in partitioned_axes else 1 for axis in range(3))
    _require(
        all(count % blocks == 0 for count, blocks in zip(nx, meshblock_grid)),
        f"d{dimension} {name}: root grid is not divisible by MeshBlock grid",
    )
    meshblock_nx = tuple(
        count // blocks for count, blocks in zip(nx, meshblock_grid)
    )
    _require(
        all(count == 1 or count >= 4 for count in meshblock_nx),
        f"d{dimension} {name}: active MeshBlock extent is below AthenaK minimum",
    )
    return meshblock_grid, meshblock_nx, math.prod(meshblock_grid)


def expected_cases() -> tuple[dict[str, object], ...]:
    """Return the exact dimension-valid source-local oracle matrix."""
    cases: list[dict[str, object]] = []
    for dimension in DIMENSIONS:
        geometry = _GEOMETRY[dimension]
        bounds = geometry["bounds"]
        for resolution in RESOLUTIONS:
            nx = geometry[resolution]
            root_cell_volume = _root_cell_volume(bounds, nx)
            for ppc in PPC_VALUES:
                qscale = required_deposit_qscale(
                    root_cell_volume=root_cell_volume, ppc=ppc
                )
                for decomposition in DECOMPOSITIONS_BY_DIMENSION[dimension]:
                    meshblock_grid, meshblock_nx, mpi_ranks = _decomposition(
                        dimension, nx, decomposition
                    )
                    for artificial_c_over_v_cr in ARTIFICIAL_C_OVER_V_CR_VALUES:
                        cases.append(
                            {
                                "case_id": _case_id(
                                    dimension,
                                    resolution,
                                    ppc,
                                    decomposition,
                                    artificial_c_over_v_cr,
                                ),
                                "dimension": dimension,
                                "resolution": resolution,
                                "ppc": ppc,
                                "decomposition": decomposition,
                                "partitioned_axes": [
                                    axis + 1
                                    for axis in _DECOMPOSITION_AXES[decomposition]
                                ],
                                "mpi_ranks": mpi_ranks,
                                "global_nx": list(nx),
                                "meshblock_grid": list(meshblock_grid),
                                "meshblock_nx": list(meshblock_nx),
                                "bounds": [list(item) for item in bounds],
                                "root_cell_volume": root_cell_volume,
                                "deposit_qscale": qscale,
                                "species_mass": SPECIES_MASS,
                                "species_charge": SPECIES_CHARGE,
                                "species_charge_over_mass": SPECIES_Q_OVER_MC,
                                "artificial_c_over_v_cr": artificial_c_over_v_cr,
                                "artificial_light_speed": (
                                    artificial_c_over_v_cr * STREAM_SPEED
                                ),
                            }
                        )
    _require(len(cases) == EXPECTED_CASE_COUNT, "oracle matrix size drifted")
    _require(
        len({case["case_id"] for case in cases}) == len(cases),
        "oracle case IDs are not unique",
    )
    for case in cases:
        _require(
            math.isclose(
                configured_volume_mean_j_over_c(case),
                EXPECTED_J_OVER_C,
                rel_tol=1.0e-13,
                abs_tol=1.0e-13,
            ),
            f"{case['case_id']}: configured deposited-current closure drifted",
        )
    return tuple(cases)


def _deck_role(dimension: int) -> str:
    return (
        "source_local_thin_2d3v_deposited_current_oracle_not_authorized"
        if dimension == 1
        else "source_local_deposited_current_oracle_not_authorized"
    )


def render_oracle_deck(case: Mapping[str, object]) -> str:
    """Render one exact runnable one-cycle output-oracle deck."""
    dimension = int(case["dimension"])
    nx = tuple(int(value) for value in case["global_nx"])
    mb = tuple(int(value) for value in case["meshblock_nx"])
    bounds = tuple(tuple(float(value) for value in item) for item in case["bounds"])
    basis = _mode_basis(dimension)
    stream = tuple(STREAM_SPEED * value for value in basis)
    per_rank = "true" if int(case["mpi_ranks"]) > 1 else "false"
    basename = str(case["case_id"]).replace("-", "_")

    lines = [
        "# Corrected Q-043 Bell deposited-current source-local oracle only.",
        "# No launch, qualification, science, or publication authority.",
        "",
        "<comment>",
        f"problem = {case['case_id']}",
        "",
        "<job>",
        f"basename = {basename}",
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
            "nlim = 1",
            "tlim = 1.0e30",
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
            f"ppc = {case['ppc']}.0",
            "pusher = boris_tsc",
            "nspecies = 1",
            "cr_distribution = center",
            "deposit_moments = true",
            "deposit_order = 2",
            f"deposit_qscale = {_float_token(float(case['deposit_qscale']))}",
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
            (
                "pic_cr_light_speed = "
                f"{_float_token(float(case['artificial_light_speed']))}"
            ),
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
            f"mass = {_float_token(float(case['species_mass']))}",
            f"charge = {_float_token(float(case['species_charge']))}",
            "",
            "<problem>",
            f"pgen_name = {PGEN_NAME}",
            "",
            f"<{PGEN_NAME}>",
            f"campaign_id = {PGEN_CAMPAIGN_ID}",
            "supersedes_campaign_id = Q023-PAPER-BELL-LINEAR",
            f"deck_role = {_deck_role(dimension)}",
            f"dimension = {dimension}",
            "source_mode = uniform_current_oracle",
            "epsilon_default = 0.4",
            "epsilon = 0.4",
            "epsilon_grid = 0.1,0.2,0.4,0.6,0.8",
            "rho = 1.0",
            "pressure = 1.0",
            "amplitude = 0.0",
            "b_g = 1.0",
            "u_a = 1.0",
            "wavelength = 1.0",
            f"k0 = {_float_token(K0)}",
            f"omega = {_float_token(SPECIES_Q_OVER_MC)}",
            f"c_over_v_cr = {case['artificial_c_over_v_cr']}.0",
            "initial_eigenmode = uniform_zero_perturbation_parallel_stream",
            "current_normalization = deposited_j_over_c_equals_2_b_g_k0",
            (
                "deposition_measure = ppc_times_deposit_qscale_times_species_"
                "charge_times_v_cr_over_root_cell_volume"
            ),
            (
                "root_cell_volume = global_root_mesh_extents_over_global_"
                "root_mesh_counts"
            ),
            "deposit_qscale_semantics = root_cell_macro_charge_volume_aware",
            "fixed_qscale = forbidden_across_dimensions_and_resolutions",
            (
                "q_over_mc_representation = "
                "species_charge_over_species_mass_equals_omega_over_b_g"
            ),
            "supersession_effect = historical_c_multiplied_normalization_invalid_for_qualification",
            f"qualification_effect = {QUALIFICATION_EFFECT}",
            "launch_authorized = false",
            "timestep = open_clean_candidate_timestep_freeze",
            "",
            "<q043_bell_current_volume_aware_deposited_current_oracle>",
            f"campaign_id = {CAMPAIGN_ID}",
            f"case_id = {case['case_id']}",
            f"resolution = {case['resolution']}",
            f"decomposition = {case['decomposition']}",
            (
                "partitioned_axes = "
                + (
                    ",".join(f"x{axis}" for axis in case["partitioned_axes"])
                    if case["partitioned_axes"]
                    else "none"
                )
            ),
            "meshblock_grid = " + ",".join(str(value) for value in case["meshblock_grid"]),
            f"mpi_ranks = {case['mpi_ranks']}",
            f"artificial_c_over_v_cr = {case['artificial_c_over_v_cr']}",
            f"species_mass = {_float_token(float(case['species_mass']))}",
            f"species_charge = {_float_token(float(case['species_charge']))}",
            (
                "species_charge_over_mass = "
                f"{_float_token(float(case['species_charge_over_mass']))}"
            ),
            f"required_output_cycle = {RAW_ORACLE_CYCLE}",
            "initial_state = uniform_zero_perturbation_parallel_stream",
            f"root_cell_volume = {_float_token(float(case['root_cell_volume']))}",
            (
                "configured_volume_mean_j_over_c = "
                f"{_float_token(configured_volume_mean_j_over_c(case))}"
            ),
            f"qualification_effect = {QUALIFICATION_EFFECT}",
            "launch_authorized = false",
            "",
        ]
    )
    for index, field in enumerate(FIELDS, 1):
        lines.extend(
            [
                f"<output{index}>",
                "file_type = bin",
                f"variable = {field}",
                f"id = {field}",
                f"dcycle = {RAW_ORACLE_DCYCLE}",
                "ghost_zones = false",
                f"single_file_per_rank = {per_rank}",
                "",
            ]
        )
    return "\n".join(lines)


def parse_athinput_text(text: str) -> dict[str, dict[str, str]]:
    """Parse the strict Athena input subset used by the generated decks."""
    blocks: dict[str, dict[str, str]] = {}
    block: str | None = None
    for line_number, raw_line in enumerate(text.splitlines(), 1):
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            block = line[1:-1].strip()
            _require(block != "" and block not in blocks, "invalid or duplicate block")
            blocks[block] = {}
            continue
        _require(block is not None and "=" in line, f"line {line_number}: malformed")
        name, value = (part.strip() for part in line.split("=", 1))
        _require(
            bool(name) and bool(value) and name not in blocks[block],
            f"line {line_number}: invalid or duplicate parameter",
        )
        blocks[block][name] = value
    return blocks


def validate_rendered_deck(case: Mapping[str, object], text: str) -> dict[str, object]:
    """Validate a rendered deck and its volume-aware configured closure."""
    blocks = parse_athinput_text(text)
    nx = tuple(int(blocks["mesh"][f"nx{axis}"]) for axis in (1, 2, 3))
    mb = tuple(int(blocks["meshblock"][f"nx{axis}"]) for axis in (1, 2, 3))
    bounds = tuple(
        (
            float(blocks["mesh"][f"x{axis}min"]),
            float(blocks["mesh"][f"x{axis}max"]),
        )
        for axis in (1, 2, 3)
    )
    root_cell_volume = _root_cell_volume(bounds, nx)
    ppc = int(float(blocks["particles"]["ppc"]))
    qscale = float(blocks["particles"]["deposit_qscale"])
    stream = tuple(
        float(blocks["particles"][f"cr_v{axis}0"]) for axis in ("x", "y", "z")
    )
    species_mass = float(blocks["species0"]["mass"])
    charge = float(blocks["species0"]["charge"])
    charge_q_over_mc = charge / species_mass
    measured = (
        ppc
        * qscale
        * charge
        * math.sqrt(sum(value * value for value in stream))
        / root_cell_volume
    )
    _require(nx == tuple(case["global_nx"]), f"{case['case_id']}: global nx drifted")
    _require(mb == tuple(case["meshblock_nx"]), f"{case['case_id']}: MeshBlock drifted")
    _require(
        math.isclose(
            root_cell_volume,
            float(case["root_cell_volume"]),
            rel_tol=1.0e-13,
            abs_tol=1.0e-13,
        ),
        f"{case['case_id']}: root-cell volume drifted",
    )
    _require(ppc == case["ppc"], f"{case['case_id']}: PPC drifted")
    _require(
        species_mass > 0.0
        and math.isclose(species_mass, float(case["species_mass"]), rel_tol=1.0e-13),
        f"{case['case_id']}: species mass drifted",
    )
    _require(
        math.isclose(charge, float(case["species_charge"]), rel_tol=1.0e-13),
        f"{case['case_id']}: species charge drifted",
    )
    _require(
        math.isclose(
            charge_q_over_mc, SPECIES_Q_OVER_MC, rel_tol=1.0e-13, abs_tol=1.0e-13
        ),
        f"{case['case_id']}: species_charge/species_mass q/(mc) drifted",
    )
    _require(
        math.isclose(qscale, float(case["deposit_qscale"]), rel_tol=1.0e-13),
        f"{case['case_id']}: qscale drifted",
    )
    _require(
        math.isclose(measured, EXPECTED_J_OVER_C, rel_tol=1.0e-13, abs_tol=1.0e-13),
        f"{case['case_id']}: volume-aware configured current drifted",
    )
    _require(
        blocks["problem"]["pgen_name"] == PGEN_NAME,
        f"{case['case_id']}: corrected pgen drifted",
    )
    _require(
        blocks["time"]["nlim"] == "1"
        and blocks[PGEN_NAME]["source_mode"] == "uniform_current_oracle"
        and blocks[PGEN_NAME]["amplitude"] == "0.0"
        and blocks[PGEN_NAME]["initial_eigenmode"]
        == "uniform_zero_perturbation_parallel_stream",
        f"{case['case_id']}: one-cycle uniform-oracle contract drifted",
    )
    _require(
        math.isclose(
            float(blocks["particles"]["pic_cr_light_speed"]),
            float(case["artificial_light_speed"]),
            rel_tol=1.0e-13,
        ),
        f"{case['case_id']}: artificial light speed drifted",
    )
    _require(
        blocks["q043_bell_current_volume_aware_deposited_current_oracle"]["case_id"]
        == case["case_id"],
        f"{case['case_id']}: oracle identity drifted",
    )
    oracle_block = blocks["q043_bell_current_volume_aware_deposited_current_oracle"]
    expected_partitioned_axes = (
        ",".join(f"x{axis}" for axis in case["partitioned_axes"])
        if case["partitioned_axes"]
        else "none"
    )
    _require(
        oracle_block["decomposition"] == case["decomposition"]
        and oracle_block["partitioned_axes"] == expected_partitioned_axes
        and oracle_block["meshblock_grid"]
        == ",".join(str(value) for value in case["meshblock_grid"])
        and int(oracle_block["mpi_ranks"]) == int(case["mpi_ranks"]),
        f"{case['case_id']}: multidirectional decomposition contract drifted",
    )
    for index, field in enumerate(FIELDS, 1):
        output = blocks[f"output{index}"]
        _require(
            output["file_type"] == "bin"
            and output["variable"] == field
            and output["id"] == field
            and output["dcycle"] == str(RAW_ORACLE_DCYCLE)
            and output["single_file_per_rank"]
            == ("true" if int(case["mpi_ranks"]) > 1 else "false"),
            f"{case['case_id']}: output contract drifted",
        )
    return {
        "case_id": case["case_id"],
        "root_cell_volume": root_cell_volume,
        "deposit_qscale": qscale,
        "species_mass": species_mass,
        "species_charge": charge,
        "species_charge_over_mass": charge_q_over_mc,
        "configured_volume_mean_j_over_c": measured,
        "deck_sha256": _sha256_bytes(text.encode("utf-8")),
    }


def build_deck_manifest() -> tuple[dict[str, object], dict[str, str]]:
    """Build and validate the exact dimension-valid source-local oracle matrix."""
    decks: dict[str, str] = {}
    records = []
    for case in expected_cases():
        relative = f"{case['case_id']}.athinput"
        text = render_oracle_deck(case)
        validation = validate_rendered_deck(case, text)
        decks[relative] = text
        records.append({**case, "deck_path": relative, **validation})
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "record_type": "q043_bell_current_volume_aware_deposited_current_oracle_deck_manifest",
        "campaign_id": CAMPAIGN_ID,
        "qualification_effect": QUALIFICATION_EFFECT,
        "launch_authorized": False,
        "scientific_claim_authorized": False,
        "publication_authorized": False,
        "current_formula": (
            "volume_mean_prtcl_j_equals_ppc_times_deposit_qscale_times_"
            "species_charge_times_v_cr_over_root_cell_volume_equals_2_b_g_k0"
        ),
        "required_output_cycle": RAW_ORACLE_CYCLE,
        "output_dcycle": RAW_ORACLE_DCYCLE,
        "output_timing": "cycle_zero_initialization_and_cycle_one_finalize_only",
        "initial_state": "uniform_zero_perturbation_parallel_stream",
        "artificial_c_in_formula": False,
        "matrix_axes": {
            "dimensions": list(DIMENSIONS),
            "resolutions": list(RESOLUTIONS),
            "ppc": list(PPC_VALUES),
            "decompositions": list(DECOMPOSITIONS),
            "decompositions_by_dimension": {
                str(dimension): list(DECOMPOSITIONS_BY_DIMENSION[dimension])
                for dimension in DIMENSIONS
            },
            "artificial_c_over_v_cr": list(ARTIFICIAL_C_OVER_V_CR_VALUES),
        },
        "case_count": len(records),
        "cases": records,
    }
    return manifest, decks


def materialize_checked_in_decks(*, replace: bool = False) -> dict[str, object]:
    """Materialize the exact additive deck matrix and its digest manifest."""
    manifest, decks = build_deck_manifest()
    if CHECKED_IN_DECK_ROOT.exists():
        _require(replace, "checked-in oracle deck root already exists")
        shutil.rmtree(CHECKED_IN_DECK_ROOT)
    CHECKED_IN_DECK_ROOT.mkdir(parents=True)
    for relative, text in decks.items():
        (CHECKED_IN_DECK_ROOT / relative).write_text(text, encoding="utf-8")
    CHECKED_IN_MANIFEST.write_bytes(_canonical_json_bytes(manifest))
    return manifest


def validate_checked_in_decks() -> dict[str, object]:
    """Validate every checked-in deck byte-for-byte against the renderer."""
    expected_manifest, decks = build_deck_manifest()
    _require(CHECKED_IN_MANIFEST.is_file(), "checked-in deck manifest is missing")
    measured_manifest = json.loads(CHECKED_IN_MANIFEST.read_text(encoding="utf-8"))
    _require(measured_manifest == expected_manifest, "checked-in deck manifest drifted")
    expected_names = set(decks) | {CHECKED_IN_MANIFEST.name}
    measured_names = {
        path.name for path in CHECKED_IN_DECK_ROOT.iterdir() if path.is_file()
    }
    _require(measured_names == expected_names, "checked-in oracle deck inventory drifted")
    for relative, text in decks.items():
        _require(
            (CHECKED_IN_DECK_ROOT / relative).read_text(encoding="utf-8") == text,
            f"checked-in oracle deck drifted: {relative}",
        )
    return expected_manifest


def _case_map() -> dict[str, dict[str, object]]:
    return {str(case["case_id"]): case for case in expected_cases()}


def expected_runtime_parameters(case: Mapping[str, object]) -> dict[str, dict[str, str]]:
    """Return the exact immutable runtime header expected from the frozen deck."""
    expected = {
        block: dict(values)
        for block, values in parse_athinput_text(render_oracle_deck(case)).items()
    }
    for block, values in _RUNTIME_DEFAULTS.items():
        expected.setdefault(block, {}).update(values)
    for block in RUNTIME_OUTPUT_BLOCKS:
        expected[block].update({"gid": "-1", "data_format": "%12.5e"})
    expected["species0"].update(
        {
            "vx0": format(float(expected["particles"]["cr_vx0"]), ".6g"),
            "vy0": format(float(expected["particles"]["cr_vy0"]), ".6g"),
            "vz0": format(float(expected["particles"]["cr_vz0"]), ".6g"),
        }
    )
    return expected


def _runtime_parameter(parameters: Mapping[str, Mapping[str, str]], block: str, name: str) -> str:
    _require(block in parameters and name in parameters[block], f"runtime {block}/{name} missing")
    return parameters[block][name]


def _normalized_runtime_parameters(
    parameters: Mapping[str, Mapping[str, str]],
    *,
    case: Mapping[str, object],
    field: str,
) -> dict[str, dict[str, str]]:
    """Validate and remove only the exact cycle-one output bookkeeping state."""
    _require(field in FIELDS, "runtime output field is unknown")
    output_blocks = {block for block in parameters if block.startswith("output")}
    _require(
        output_blocks == set(RUNTIME_OUTPUT_BLOCKS),
        f"{field}: runtime output block inventory drifted",
    )
    field_index = FIELDS.index(field) + 1
    normalized = {}
    for block, values in parameters.items():
        _require(
            isinstance(block, str) and isinstance(values, Mapping),
            f"{field}: runtime parameter block is malformed",
        )
        numbered_output = block in RUNTIME_OUTPUT_BLOCKS
        if numbered_output:
            output_index = int(block[6:])
            _require(
                values.get("file_type") == "bin"
                and values.get("variable") == FIELDS[output_index - 1]
                and values.get("id") == FIELDS[output_index - 1]
                and values.get("dcycle") == str(RAW_ORACLE_DCYCLE)
                and values.get("ghost_zones") == "false"
                and values.get("single_file_per_rank")
                == ("true" if int(case["mpi_ranks"]) > 1 else "false"),
                f"{field}: runtime {block} immutable contract drifted",
            )
            _require(
                RUNTIME_OUTPUT_BOOKKEEPING_KEYS <= set(values),
                f"{field}: runtime {block} output bookkeeping is incomplete",
            )
            file_number = values["file_number"]
            try:
                parsed_file_number = int(file_number)
            except ValueError as error:
                raise ContractError(
                    f"{field}: runtime {block} file_number must be an integer"
                ) from error
            _require(
                file_number == str(parsed_file_number) and parsed_file_number >= 0,
                f"{field}: runtime {block} file_number must be a canonical "
                "non-negative integer",
            )
            expected_file_number = 1 if output_index <= field_index else 2
            _require(
                parsed_file_number == expected_file_number,
                f"{field}: runtime {block} file_number violates the observed "
                "cycle-one sequential publication contract",
            )
            try:
                last_time = float(values["last_time"])
            except ValueError as error:
                raise ContractError(
                    f"{field}: runtime {block} last_time must be numeric"
                ) from error
            _require(
                math.isfinite(last_time) and last_time == 0.0,
                f"{field}: runtime {block} last_time violates the cycle-cadence "
                "publication contract",
            )
        normalized[block] = {
            key: value
            for key, value in values.items()
            if not (numbered_output and key in RUNTIME_OUTPUT_BOOKKEEPING_KEYS)
        }
    return normalized


def _validate_runtime_dataset(
    dataset: binary.AthenaBinaryDataset,
    *,
    case: Mapping[str, object],
    field: str,
) -> None:
    _require(
        dataset.time > 0.0 and dataset.cycle == RAW_ORACLE_CYCLE,
        "oracle output must follow one complete cycle",
    )
    _require(dataset.variable_names == (field,), f"{field}: raw variable schema drifted")
    _require(dataset.root_grid_shape == tuple(case["global_nx"]), f"{field}: root grid drifted")
    _require(dataset.meshblock_shape == tuple(case["meshblock_nx"]), f"{field}: MeshBlock drifted")
    parameters = dataset.input_parameters
    _require(
        _normalized_runtime_parameters(parameters, case=case, field=field)
        == expected_runtime_parameters(case),
        f"{field}: runtime parameters drifted from the authoritative deck and "
        "frozen default contract",
    )
    _require(_runtime_parameter(parameters, "problem", "pgen_name") == PGEN_NAME, "runtime pgen drifted")
    _require(
        _runtime_parameter(parameters, "q043_bell_current_volume_aware_deposited_current_oracle", "case_id")
        == case["case_id"],
        "runtime oracle case identity drifted",
    )
    _require(
        int(float(_runtime_parameter(parameters, "particles", "ppc"))) == case["ppc"],
        "runtime PPC drifted",
    )
    _require(
        math.isclose(
            float(_runtime_parameter(parameters, "particles", "deposit_qscale")),
            float(case["deposit_qscale"]),
            rel_tol=1.0e-13,
        ),
        "runtime qscale drifted",
    )
    _require(
        math.isclose(
            float(_runtime_parameter(parameters, "species0", "mass")),
            float(case["species_mass"]),
            rel_tol=1.0e-13,
        )
        and math.isclose(
            float(_runtime_parameter(parameters, "species0", "charge")),
            float(case["species_charge"]),
            rel_tol=1.0e-13,
        ),
        "runtime species mass or charge drifted",
    )
    _require(
        math.isclose(
            float(_runtime_parameter(parameters, "species0", "charge"))
            / float(_runtime_parameter(parameters, "species0", "mass")),
            SPECIES_Q_OVER_MC,
            rel_tol=1.0e-13,
            abs_tol=1.0e-13,
        ),
        "runtime species_charge/species_mass q/(mc) drifted",
    )
    _require(
        math.isclose(
            float(_runtime_parameter(parameters, "particles", "pic_cr_light_speed")),
            float(case["artificial_light_speed"]),
            rel_tol=1.0e-13,
        ),
        "runtime artificial light speed drifted",
    )
    _require(
        _runtime_parameter(parameters, PGEN_NAME, "amplitude") == "0.0"
        and _runtime_parameter(parameters, PGEN_NAME, "source_mode")
        == "uniform_current_oracle"
        and _runtime_parameter(parameters, PGEN_NAME, "initial_eigenmode")
        == "uniform_zero_perturbation_parallel_stream",
        "runtime zero-perturbation parallel-stream contract drifted",
    )


def _raw_binding(path: Path, *, field: str, shard_index: int) -> dict[str, object]:
    resolved = path.resolve(strict=True)
    _require(resolved.is_file(), f"raw {field} shard is not a regular file")
    payload = resolved.read_bytes()
    return {
        "field": field,
        "shard_index": shard_index,
        "path": str(resolved),
        "size": len(payload),
        "sha256": _sha256_bytes(payload),
    }


def _volume_mean(values: np.ndarray, grid: binary.CompositeGrid) -> float:
    volumes = (
        np.diff(grid.x3_faces)[:, None, None]
        * np.diff(grid.x2_faces)[None, :, None]
        * np.diff(grid.x1_faces)[None, None, :]
    )
    return float(np.sum(values * volumes) / np.sum(volumes))


def _representation_tolerance(variable_size: int) -> float:
    dtype = np.float32 if variable_size == 4 else np.float64
    return DEPOSITION_TOLERANCE_ULPS * np.finfo(dtype).eps * EXPECTED_J_OVER_C


def analyze_raw_case(
    case_id: str, field_paths: Mapping[str, Sequence[Path]]
) -> dict[str, object]:
    """Analyze one post-step raw-output case and bind all raw provenance."""
    cases = _case_map()
    _require(case_id in cases, "unknown deposited-current oracle case")
    case = cases[case_id]
    _require(set(field_paths) == set(FIELDS), f"{case_id}: raw field set drifted")
    shard_count = int(case["mpi_ranks"])
    merged: dict[str, binary.AthenaBinaryDataset] = {}
    grids: dict[str, binary.CompositeGrid] = {}
    provenance = []
    for field in FIELDS:
        paths = list(field_paths[field])
        _require(len(paths) == shard_count, f"{case_id}: {field} shard count drifted")
        _require(len({str(path) for path in paths}) == len(paths), f"{case_id}: duplicate shard path")
        datasets = []
        for shard_index, path in enumerate(paths):
            provenance.append(_raw_binding(path, field=field, shard_index=shard_index))
            try:
                dataset = binary.read_athenak_binary(path)
            except binary.AnalysisError as error:
                raise ContractError(f"{case_id}: malformed raw {field} output") from error
            _validate_runtime_dataset(dataset, case=case, field=field)
            datasets.append(dataset)
        try:
            merged[field] = binary.merge_athenak_binary_datasets(datasets)
            grids[field] = binary.compose_leaf_field(merged[field], field)
        except binary.AnalysisError as error:
            raise ContractError(f"{case_id}: incomplete or inconsistent {field} shards") from error

    reference = merged[FIELDS[0]]
    for field in FIELDS[1:]:
        dataset = merged[field]
        _require(
            (
                dataset.time,
                dataset.cycle,
                dataset.location_size,
                dataset.variable_size,
                dataset.root_grid_shape,
                dataset.meshblock_shape,
                dataset.domain_bounds,
                _normalized_runtime_parameters(
                    dataset.input_parameters, case=case, field=field
                ),
            )
            == (
                reference.time,
                reference.cycle,
                reference.location_size,
                reference.variable_size,
                reference.root_grid_shape,
                reference.meshblock_shape,
                reference.domain_bounds,
                _normalized_runtime_parameters(
                    reference.input_parameters, case=case, field=FIELDS[0]
                ),
            ),
            f"{case_id}: raw field metadata disagrees",
        )
        _require(
            np.array_equal(grids[field].x1_faces, grids[FIELDS[0]].x1_faces)
            and np.array_equal(grids[field].x2_faces, grids[FIELDS[0]].x2_faces)
            and np.array_equal(grids[field].x3_faces, grids[FIELDS[0]].x3_faces),
            f"{case_id}: raw field geometry disagrees",
        )

    current = np.stack([grids[field].values for field in FIELDS[1:]], axis=0)
    rho = grids["prtcl_rho"].values
    basis = np.asarray(_mode_basis(int(case["dimension"])), dtype=np.float64)
    expected_vector = EXPECTED_J_OVER_C * basis
    parallel = np.tensordot(basis, current, axes=1)
    transverse = current - basis[:, None, None, None] * parallel[None, ...]
    residual = current - expected_vector[:, None, None, None]
    mean_vector = np.asarray(
        [_volume_mean(current[index], grids[FIELDS[index + 1]]) for index in range(3)]
    )
    parallel_mean = float(np.dot(basis, mean_vector))
    rho_mean = _volume_mean(rho, grids["prtcl_rho"])
    max_transverse = float(np.max(np.sqrt(np.sum(transverse * transverse, axis=0))))
    max_current_nonuniformity = float(np.max(np.sqrt(np.sum(residual * residual, axis=0))))
    max_rho_nonuniformity = float(np.max(np.abs(rho - EXPECTED_RHO)))
    tolerance = _representation_tolerance(reference.variable_size)
    _require(
        abs(parallel_mean - EXPECTED_J_OVER_C) <= tolerance,
        f"{case_id}: guide-field projected deposited current failed",
    )
    _require(max_transverse <= tolerance, f"{case_id}: transverse deposited current failed")
    _require(
        max_current_nonuniformity <= tolerance,
        f"{case_id}: deposited-current spatial nonuniformity failed",
    )
    _require(abs(rho_mean - EXPECTED_RHO) <= tolerance, f"{case_id}: deposited rho mean failed")
    _require(
        max_rho_nonuniformity <= tolerance,
        f"{case_id}: deposited rho spatial nonuniformity failed",
    )
    _require(
        math.isclose(
            configured_volume_mean_j_over_c(case),
            EXPECTED_J_OVER_C,
            rel_tol=1.0e-13,
            abs_tol=1.0e-13,
        ),
        f"{case_id}: configured volume-aware closure failed",
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": "q043_bell_current_volume_aware_deposited_current_raw_case_oracle",
        "campaign_id": CAMPAIGN_ID,
        "case_id": case_id,
        "qualification_effect": QUALIFICATION_EFFECT,
        "launch_authorized": False,
        "scientific_claim_authorized": False,
        "publication_authorized": False,
        "passed": False,
        "source_local_oracle_check_pass": True,
        "case_contract": case,
        "raw_provenance": provenance,
        "cross_field_runtime_metadata": RUNTIME_METADATA_CONTRACT,
        "representation_derived_absolute_tolerance": tolerance,
        "configured_volume_mean_j_over_c": configured_volume_mean_j_over_c(case),
        "measured_volume_mean_j_over_c_vector": mean_vector.tolist(),
        "measured_guide_projected_volume_mean_j_over_c": parallel_mean,
        "measured_volume_mean_prtcl_rho": rho_mean,
        "maximum_transverse_current": max_transverse,
        "maximum_current_spatial_nonuniformity": max_current_nonuniformity,
        "maximum_rho_spatial_nonuniformity": max_rho_nonuniformity,
        "tolerance_scope": (
            "source_local_exact_deposition_oracle_representation_roundoff_only_"
            "not_science_acceptance"
        ),
    }


def analyze_raw_matrix(
    raw_cases: Mapping[str, Mapping[str, Sequence[Path]]]
) -> dict[str, object]:
    """Require and analyze the complete Q043 foundational current matrix."""
    expected = _case_map()
    _require(set(raw_cases) == set(expected), "raw oracle matrix is incomplete or unknown")
    results = [analyze_raw_case(case_id, raw_cases[case_id]) for case_id in expected]
    all_paths = [
        binding["path"] for result in results for binding in result["raw_provenance"]
    ]
    _require(len(all_paths) == len(set(all_paths)), "raw artifact reused across oracle cases")
    projected = [
        float(result["measured_guide_projected_volume_mean_j_over_c"])
        for result in results
    ]
    maximum_tolerance = max(
        float(result["representation_derived_absolute_tolerance"]) for result in results
    )
    spread = max(projected) - min(projected)
    _require(spread <= maximum_tolerance, "cross-matrix deposited-current invariance failed")
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": "q043_bell_current_volume_aware_deposited_current_raw_matrix_oracle",
        "campaign_id": CAMPAIGN_ID,
        "qualification_effect": QUALIFICATION_EFFECT,
        "launch_authorized": False,
        "scientific_claim_authorized": False,
        "publication_authorized": False,
        "passed": False,
        "source_local_oracle_check_pass": True,
        "case_count": len(results),
        "dimensions_verified": list(DIMENSIONS),
        "resolutions_verified": list(RESOLUTIONS),
        "ppc_verified": list(PPC_VALUES),
        "decompositions_verified": list(DECOMPOSITIONS),
        "artificial_c_over_v_cr_verified": list(ARTIFICIAL_C_OVER_V_CR_VALUES),
        "required_output_cycle": RAW_ORACLE_CYCLE,
        "output_dcycle": RAW_ORACLE_DCYCLE,
        "output_timing": "cycle_zero_initialization_and_cycle_one_finalize_only",
        "cross_field_runtime_metadata": RUNTIME_METADATA_CONTRACT,
        "initial_state_verified": "uniform_zero_perturbation_parallel_stream",
        "artificial_c_in_formula": False,
        "configured_current_formula": (
            "PPC*deposit_qscale*species_charge*v_CR/V_root_cell=2*B_g*k0"
        ),
        "maximum_cross_matrix_projected_current_spread": spread,
        "cases": results,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--materialize-checked-in-decks", action="store_true")
    parser.add_argument("--replace", action="store_true")
    parser.add_argument("--validate-checked-in-decks", action="store_true")
    args = parser.parse_args()
    _require(
        args.materialize_checked_in_decks != args.validate_checked_in_decks,
        "select exactly one deck operation",
    )
    result = (
        materialize_checked_in_decks(replace=args.replace)
        if args.materialize_checked_in_decks
        else validate_checked_in_decks()
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
