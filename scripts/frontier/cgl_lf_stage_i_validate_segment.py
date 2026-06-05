#!/usr/bin/env python3
"""Independently validate one schema-4 CGL-LF Stage I segment inspection."""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import stat
import struct
import subprocess
import sys
import tempfile
from typing import BinaryIO, Callable, Iterator


EXECUTION_EPOCH = "E03-forcing-policy"
EXECUTION_EPOCH_SLUG = "E03_forcing_policy"
CANONICAL_PROJECT_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/CGL")
MANIFEST_SCHEMA_VERSION = 3
INSPECTION_SCHEMA_VERSION = 4
MAX_RESTART_PARAMETER_DUMP_BYTES = 11 * 4096 + 1
MAX_TEXT_LINE_BYTES = 64 * 1024
NORMALIZATION_FLOOR = 1.0
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
GIT_REVISION_PATTERN = re.compile(r"[0-9a-f]{40}")
JOB_ID_PATTERN = re.compile(r"[1-9][0-9]*")
SEGMENT_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]{0,28}")
RANK_DIRECTORY_PATTERN = re.compile(r"rank_[0-9]{8}")
HISTORY_LABEL_PATTERN = re.compile(r"\[(\d+)\]=(\S+)")
BATCH_SCRIPT_DIGEST_PATTERN = re.compile(
    r"(?m)^BATCH_SCRIPT_SHA256=([0-9a-f]{64})$"
)
BATCH_SCRIPT_DIGEST_PLACEHOLDER = "0" * 64
ACCOUNT = "AST207"
PARTITION = "batch"
EXPECTED_RANKS_PER_NODE = 8
EXPECTED_CPUS_PER_TASK = 7
MAX_SEGMENT_SECONDS = 2 * 60 * 60
MATRIX_REPOSITORY_PATH = "inputs/cgl_lf_paper/mks24_stage_i_manifest.json"
PRODUCTION_UTILITY_REPOSITORY_PATH = "scripts/frontier/cgl_lf_stage_i.py"
VALIDATOR_REPOSITORY_PATH = "scripts/frontier/cgl_lf_stage_i_validate_segment.py"
GIT_EXECUTABLE = Path("/usr/lib/git/git")
EXPECTED_SNAPSHOT_VARIABLES = (
    "dens",
    "velx",
    "vely",
    "velz",
    "eint",
    "p_perp",
    "bcc1",
    "bcc2",
    "bcc3",
)

APPROVED_CASE_IDS = frozenset(f"R{number:02d}" for number in range(2, 18))
PASSIVE_CASE_IDS = frozenset({"R06", "R07", "R08", "R09"})
FINITE_LIMITER_CASE_IDS = frozenset({"R14", "R15"})
CASE_NODE_PROFILES = {
    "R02": frozenset({1}),
    "R03": frozenset({1}),
    **{
        f"R{number:02d}": frozenset({1, 2, 4})
        for number in range(4, 16)
    },
    "R16": frozenset({1, 2}),
    "R17": frozenset({8}),
}
CASE_NAMES = {
    "R02": "paper_standard_active_alfvenic_beta10",
    "R03": "paper_standard_active_alfvenic_beta100",
    "R04": "paper_standard_active_random_beta10",
    "R05": "paper_standard_active_random_beta100",
    "R06": "paper_standard_passive_alfvenic_beta10",
    "R07": "paper_standard_passive_alfvenic_beta100",
    "R08": "paper_standard_passive_random_beta10",
    "R09": "paper_standard_passive_random_beta100",
    "R10": "paper_compressive_active_random_beta1",
    "R11": "paper_compressive_active_random_beta100_sonic",
    "R12": "paper_heat_flux_beta10_strong",
    "R13": "paper_heat_flux_beta10_weak",
    "R14": "paper_nulim_beta100_20",
    "R15": "paper_nulim_beta100_200",
    "R16": "paper_scale_separation_active_alfvenic_beta10_nperp96",
    "R17": "paper_scale_separation_active_alfvenic_beta10_nperp384",
}
CASE_INPUT_PATHS = {
    case_id: f"inputs/cgl_lf_paper/cgl_lf_{name}.athinput"
    for case_id, name in CASE_NAMES.items()
}
CASE_INPUT_PATHS.update({
    "R16": "inputs/cgl_lf_paper/cgl_lf_paper_scale_separation_beta10_nperp96.athinput",
    "R17": "inputs/cgl_lf_paper/cgl_lf_paper_scale_separation_beta10_nperp384.athinput",
})
CASE_INPUT_CONTRACTS = {
    "R02": ("192x192x384", 10.0, False, "alfvenic", 6.283185307179586, 1.0e10),
    "R03": ("192x192x384", 100.0, False, "alfvenic", 6.283185307179586, 1.0e10),
    "R04": ("192x192x384", 10.0, False, "random", 6.283185307179586, 1.0e10),
    "R05": ("192x192x384", 100.0, False, "random", 6.283185307179586, 1.0e10),
    "R06": ("192x192x384", 10.0, True, "alfvenic", 6.283185307179586, 1.0e10),
    "R07": ("192x192x384", 100.0, True, "alfvenic", 6.283185307179586, 1.0e10),
    "R08": ("192x192x384", 10.0, True, "random", 6.283185307179586, 1.0e10),
    "R09": ("192x192x384", 100.0, True, "random", 6.283185307179586, 1.0e10),
    "R10": ("192x192x384", 1.0, False, "random", 6.283185307179586, 1.0e10),
    "R11": ("192x192x384", 100.0, False, "random_sonic", 6.283185307179586, 1.0e10),
    "R12": ("192x192x384", 10.0, False, "alfvenic", 0.06283185307179586, 1.0e10),
    "R13": ("192x192x384", 10.0, False, "alfvenic", 628.3185307179587, 1.0e10),
    "R14": ("192x192x384", 100.0, False, "alfvenic", 6.283185307179586, 20.0),
    "R15": ("192x192x384", 100.0, False, "alfvenic", 6.283185307179586, 200.0),
    "R16": ("96x96x192", 10.0, False, "alfvenic", 6.283185307179586, 1.0e10),
    "R17": ("384x384x768", 10.0, False, "alfvenic", 6.283185307179586, 1.0e10),
}
COMMON_FROZEN_PARAMETERS = {
    ("mesh", "nghost"): "2",
    ("mesh", "x1min"): "0.0",
    ("mesh", "x1max"): "1.0",
    ("mesh", "ix1_bc"): "periodic",
    ("mesh", "ox1_bc"): "periodic",
    ("mesh", "x2min"): "0.0",
    ("mesh", "x2max"): "1.0",
    ("mesh", "ix2_bc"): "periodic",
    ("mesh", "ox2_bc"): "periodic",
    ("mesh", "x3min"): "0.0",
    ("mesh", "x3max"): "2.0",
    ("mesh", "ix3_bc"): "periodic",
    ("mesh", "ox3_bc"): "periodic",
    ("meshblock", "nx1"): "32",
    ("meshblock", "nx2"): "32",
    ("meshblock", "nx3"): "64",
    ("time", "evolution"): "dynamic",
    ("time", "integrator"): "rk2",
    ("time", "sts_integrator"): "rkl2",
    ("time", "sts_max_dt_ratio"): "-1.0",
    ("time", "cfl_number"): "0.3",
    ("time", "nlim"): "-1",
    ("time", "ndiag"): "10",
    ("mhd", "eos"): "cgl",
    ("mhd", "cgl_heat_flux"): "landau_fluid",
    ("mhd", "cgl_heat_flux_integrator"): "sts",
    ("mhd", "cgl_lf_strict_admissibility"): "true",
    ("mhd", "cgl_lf_record_pressure_work"): "true",
    ("mhd", "lf_coefficient_mode"): "local",
    ("mhd", "nu_coll"): "0.0",
    ("mhd", "mirror_limiter"): "true",
    ("mhd", "firehose_limiter"): "true",
    ("mhd", "cgl_firehose_threshold"): "parallel",
    ("mhd", "backup_limiters"): "false",
    ("mhd", "reconstruct"): "plm",
    ("mhd", "rsolver"): "hlle",
    ("mhd", "gamma"): "1.666666666666667",
    ("mhd", "dfloor"): "1.0e-12",
    ("mhd", "pfloor"): "1.0e-12",
    ("mhd", "tfloor"): "1.0e-12",
    ("mhd", "sfloor"): "1.0e-12",
    ("mhd", "bfloor"): "1.0e-10",
    ("problem", "pgen_name"): "cgl_lf_paper",
    ("problem", "user_hist"): "true",
    ("problem", "paper_mode"): "turbulence",
    ("problem", "rho0"): "1.0",
    ("problem", "b0"): "1.0",
    ("problem", "analysis_t_start"): "8.0",
    ("problem", "analysis_t_end"): "10.0",
    ("turb_driving", "spectrum"): "power_law",
    ("turb_driving", "nlow"): "1",
    ("turb_driving", "nhigh"): "3",
    ("turb_driving", "physical_k_shell"): "true",
    ("turb_driving", "k_shell_unit"): "3.141592653589793",
    ("turb_driving", "isotropic_power_spectrum"): "true",
    ("turb_driving", "expo"): "2.0",
    ("turb_driving", "dedt"): "0.32",
    ("turb_driving", "record_injected_work"): "true",
    ("output1", "file_type"): "hst",
    ("output1", "data_format"): "%24.16e",
    ("output1", "dt"): "0.02",
    ("output2", "file_type"): "bin",
    ("output2", "variable"): "mhd_w_bcc",
    ("output2", "dt"): "0.25",
    ("output2", "single_file_per_rank"): "true",
    ("output3", "file_type"): "rst",
    ("output3", "dt"): "1.0",
    ("output3", "single_file_per_rank"): "true",
}
STRICT_LF_FAILURE_COLUMNS = (
    "lf_dfloor",
    "lf_pfloor",
    "lf_nonfin",
    "lf_nonpos",
    "lf_hardbd",
)
LF_REQUIRED_COLUMNS = frozenset({
    "lf_nstage",
    "lf_qface",
    "lf_qprcap",
    "lf_qpr10",
    "lf_qpecap",
    "lf_qpe10",
    "lf_qprwrk",
    "lf_qpewrk",
    "lf_hwproj",
    "lf_cpwrk",
    "lf_cawrk",
})
LF_MONOTONIC_COUNT_COLUMNS = (
    "lf_nstage",
    "lf_mirror",
    "lf_firehs",
    "lf_qface",
    "lf_qprcap",
    "lf_qpr10",
    "lf_qpecap",
    "lf_qpe10",
)
LF_WORK_COLUMNS = ("lf_qprwrk", "lf_qpewrk", "lf_cpwrk", "lf_cawrk")
LF_RESTART_HISTORY_COLUMNS = (
    "lf_nstage",
    "lf_dfloor",
    "lf_pfloor",
    "lf_nonfin",
    "lf_nonpos",
    "lf_mirror",
    "lf_firehs",
    "lf_hardbd",
    "lf_qface",
    "lf_qprcap",
    "lf_qpr10",
    "lf_qpecap",
    "lf_qpe10",
    "lf_qprwrk",
    "lf_qpewrk",
    "lf_cpwrk",
    "lf_cawrk",
    "lf_hwproj",
)
COUNT_DIAGNOSTIC_INDICES = tuple(range(13)) + (17,)
HISTORICAL_MHD_HISTORY_COLUMNS = (
    "time",
    "dt",
    "mass",
    "1-mom",
    "2-mom",
    "3-mom",
    "tot-E",
    "aam-D",
    "1-KE",
    "2-KE",
    "3-KE",
    "1-ME",
    "2-ME",
    "3-ME",
    "lf_nstage",
    "lf_dfloor",
    "lf_pfloor",
    "lf_nonfin",
    "lf_nonpos",
    "lf_mirror",
    "lf_firehs",
    "lf_hardbd",
    "lf_qface",
    "lf_qprcap",
    "lf_qpr10",
    "lf_qpecap",
    "lf_qpe10",
    "lf_qprwrk",
    "lf_qpewrk",
    "lf_hwproj",
    "lf_cpwrk",
    "lf_cawrk",
)
HISTORICAL_USER_HISTORY_COLUMNS = (
    "time",
    "dt",
    "volume",
    "p_parallel",
    "p_perp",
    "force_prp2",
    "force_prl2",
    "vel_prp2",
    "mass",
    "kinetic",
    "magnetic",
    "therm_cgl",
    "b2",
    "b4",
    "delta_p",
    "abs_dp",
    "beta",
    "mirror_vol",
    "fire_vol",
    "hard_vol",
    "nu_eff",
    "force_pwr",
    "vel_prl2",
    "force_work",
)
ARCHIVED_INPUT_DYNAMIC_PARAMETERS = frozenset({
    ("time", "tlim"),
})
ARCHIVED_INPUT_OPTIONAL_DEFAULTS = {
    ("mhd", "nscalars"): "0",
}
QUALIFIED_PRODUCT_PARAMETER_ADDITIONS = {
    ("coord", "general_rel"): "0",
    ("coord", "special_rel"): "0",
    ("mesh_refinement", "refinement"): "none",
    ("mhd", "fofc"): "0",
    ("mhd", "nscalars"): "0",
    ("mhd", "sigma_max"): "3.40282e+38",
    ("output1", "file_number"): None,
    ("output1", "ghost_zones"): "0",
    ("output1", "gid"): "-1",
    ("output1", "last_time"): None,
    ("output1", "user_hist_only"): "0",
    ("output2", "data_format"): "%12.5e",
    ("output2", "file_number"): None,
    ("output2", "ghost_zones"): "0",
    ("output2", "gid"): "-1",
    ("output2", "id"): "mhd_w_bcc",
    ("output2", "last_time"): None,
    ("output3", "data_format"): "%12.5e",
    ("output3", "file_number"): None,
    ("output3", "ghost_zones"): "0",
    ("output3", "gid"): "-1",
    ("output3", "last_time"): None,
    ("problem", "user_srcs"): "0",
    ("time", "restart_time"): None,
    ("time", "start_time"): "0",
    ("turb_driving", "center_x1"): "0",
    ("turb_driving", "center_x2"): "0",
    ("turb_driving", "center_x3"): "0",
    ("turb_driving", "dt_update"): "0.01",
    ("turb_driving", "exp_prl"): "0",
    ("turb_driving", "exp_prp"): "1.66667",
    ("turb_driving", "kpeak"): "12.5664",
    ("turb_driving", "localization"): "none",
    ("turb_driving", "max_kx"): "3",
    ("turb_driving", "max_ky"): "3",
    ("turb_driving", "max_kz"): "3",
    ("turb_driving", "min_kx"): "0",
    ("turb_driving", "min_ky"): "0",
    ("turb_driving", "min_kz"): "0",
    ("turb_driving", "normalization"): "edot",
    ("turb_driving", "sigma_x1"): "-1",
    ("turb_driving", "sigma_x2"): "-1",
    ("turb_driving", "sigma_x3"): "-1",
    ("turb_driving", "sol_fraction"): "1.0",
    ("turb_driving", "tdriv_start"): "0",
    ("turb_driving", "tile_nx"): "1",
    ("turb_driving", "tile_ny"): "1",
    ("turb_driving", "tile_nz"): "1",
    ("turb_driving", "turb_flag"): "2",
}
PRODUCT_NONNEGATIVE_INTEGER_PARAMETERS = frozenset({
    ("output1", "file_number"),
    ("output2", "file_number"),
    ("output3", "file_number"),
})
PRODUCT_NONNEGATIVE_FINITE_PARAMETERS = frozenset({
    ("time", "restart_time"),
})
PRODUCT_SCHEDULE_TIME_PARAMETERS = frozenset({
    ("output1", "last_time"),
    ("output2", "last_time"),
    ("output3", "last_time"),
})
QUALIFIED_SNAPSHOT_OPTIONAL_PARAMETERS = frozenset({
    ("time", "restart_time"),
})
TURBULENCE_METADATA_INT_FIELDS = (
    "version",
    "mode_count",
    "n_updates",
    "nlow",
    "nhigh",
    "driving_type",
    "min_kx",
    "max_kx",
    "min_ky",
    "max_ky",
    "min_kz",
    "max_kz",
    "use_npeak",
    "turb_flag",
    "tile_nx",
    "tile_ny",
    "tile_nz",
    "normalization",
    "localization",
    "spectrum",
    "projection_policy",
    "physical_k_shell",
    "isotropic_power_spectrum",
    "record_injected_work",
)
TURBULENCE_METADATA_REAL_FIELDS = (
    "tcorr",
    "dt_update",
    "dedt",
    "accel_rms",
    "sol_fraction",
    "kpeak",
    "npeak",
    "expo",
    "exp_prp",
    "exp_prl",
    "tdriv_duration",
    "tdriv_start",
    "sigma_x1",
    "sigma_x2",
    "sigma_x3",
    "center_x1",
    "center_x2",
    "center_x3",
    "k_shell_unit",
)
TURBULENCE_METADATA_CONFIGURATION_FIELDS = (
    tuple(name for name in TURBULENCE_METADATA_INT_FIELDS if name != "n_updates")
    + TURBULENCE_METADATA_REAL_FIELDS
)
QUALIFIED_RESTART_BINARY_ABIS = {
    (
        "9e07542281e4e6d125582f253df3ad2e3b8b154d",
        "68f243f9204df388b24365ae65a567f6f567dbe422a6d7a43b9fb4a499ef118c",
    ): {
        "mesh_header_size": 252,
        "mesh_time_offset_after_parameter_dump": 232,
        "mesh_time_format": "<d",
        "logical_location_size": 16,
        "cost_size": 4,
        "turbulence_metadata_size": 248,
        "rng_state_size": 296,
        "rng_state_format": "<35qi4xd",
        "rng_state_generator_format": "<35q",
        "rng_state_generator_field_count": 35,
        "rng_state_iset_format": "<i",
        "rng_state_conditional_gset_format": "<d",
        "rng_state_padding_offset": 284,
        "rng_state_padding_size": 4,
        "lf_diagnostic_count": 18,
        "real_size": 8,
        "mhd_history_columns": HISTORICAL_MHD_HISTORY_COLUMNS,
        "user_history_columns": HISTORICAL_USER_HISTORY_COLUMNS,
        "normalized_ct_divb_history": "not_retained",
        # The qualified executable writes uninitialized non-semantic coarse
        # fields in mesh_indcs for non-AMR runs. MeshBlock RegionIndcs and all
        # active mesh RegionIndcs fields remain exact and authoritatively
        # validated.
        "mesh_coarse_region_fields": "legacy_opaque_uninitialized",
        "allowed_marker_modes": frozenset({
            "full_precision",
            "legacy_default_precision",
        }),
    },
}
APPROVED_THRESHOLD_POLICIES = {
    "stage-i-standard-v1": {
        "forcing": {
            "normalized_residual_lt": 1.0e-8,
        },
        "activity": {
            "forcing_absolute_gt": 1.0e-6,
            "forcing_state_normalized_gt": 1.0e-8,
            "active_pressure_work_absolute_gt": 1.0e-6,
            "active_pressure_work_state_normalized_gt": 1.0e-8,
        },
        "restart_history": {
            "work_absolute_tolerance": 1.0e-12,
            "work_relative_tolerance": 1.0e-12,
        },
        "mass_relative_drift_le": 1.0e-12,
        "mass_relative_mismatch_le": 1.0e-12,
    },
    "r03-4762472-clean-partial-v1": {
        "forcing": {
            "normalized_residual_lt": 1.0e-8,
            "absolute_residual_lt": 1.0e-10,
            "increment_normalized_residual_lt": 1.0e-9,
            "state_energy_normalized_residual_lt": 1.0e-12,
        },
        "activity": {
            "forcing_absolute_gt": 1.0e-6,
            "forcing_state_normalized_gt": 1.0e-8,
            "active_pressure_work_absolute_gt": 1.0e-6,
            "active_pressure_work_state_normalized_gt": 1.0e-8,
        },
        "restart_history": {
            "work_absolute_tolerance": 1.0e-12,
            "work_relative_tolerance": 1.0e-12,
        },
        "mass_relative_drift_le": 1.0e-12,
        "mass_relative_mismatch_le": 1.0e-12,
        "binding": {
            "case_id": "R03",
            "job_id": "4762472",
            "segment": "s00_rankio_t0_t0p5",
            "result": "clean_partial",
            "required_time": 0.5,
            "final_time": 0.31282347945569927,
        },
    },
}


class ValidationError(ValueError):
    """Report one fail-closed independent-validation rejection."""


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse the independent validator CLI."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument(
        "--inspection",
        type=Path,
        help="schema-4 inspection; defaults to MANIFEST's segment_inspection.json",
    )
    parser.add_argument(
        "--result",
        choices=("accepted", "clean_partial"),
        default="accepted",
        help="expected segment result (default: accepted)",
    )
    parser.add_argument(
        "--policy",
        choices=tuple(sorted(APPROVED_THRESHOLD_POLICIES)),
        required=True,
        help="named, immutable independent-validation threshold policy",
    )
    return parser.parse_args(argv)


def cli_absolute_path(path: Path) -> Path:
    """Return one absolute CLI path without resolving symlinks."""

    return Path(os.path.abspath(path.expanduser()))


def declared_path(value: object, label: str) -> Path:
    """Require one normalized absolute retained path without traversal."""

    if not isinstance(value, str) or not value:
        raise ValidationError(f"{label} path is missing")
    if any(ord(character) < 32 or ord(character) == 127 for character in value):
        raise ValidationError(f"{label} path contains a control character")
    path = Path(value)
    if not path.is_absolute() or ".." in path.parts:
        raise ValidationError(f"{label} path is not normalized and absolute: {path}")
    if Path(os.path.normpath(value)) != path:
        raise ValidationError(f"{label} path is not normalized and absolute: {path}")
    return path


def require_exact_int(value: object, label: str) -> int:
    """Require one JSON integer while rejecting booleans and coercible strings."""

    if type(value) is not int:
        raise ValidationError(f"{label} is not an exact integer")
    return value


def parse_walltime(value: object, label: str) -> int:
    """Parse one exact HH:MM:SS production walltime."""

    if not isinstance(value, str):
        raise ValidationError(f"{label} is not an HH:MM:SS string")
    match = re.fullmatch(r"([0-9]{2}):([0-5][0-9]):([0-5][0-9])", value)
    if match is None:
        raise ValidationError(f"{label} is not an HH:MM:SS string")
    hours, minutes, seconds = (int(item) for item in match.groups())
    return hours * 3600 + minutes * 60 + seconds


def stat_profile(profile: os.stat_result) -> tuple[int, ...]:
    """Return fields that bind one named file version and ownership profile."""

    return (
        profile.st_dev,
        profile.st_ino,
        profile.st_mode,
        profile.st_nlink,
        profile.st_uid,
        profile.st_gid,
        profile.st_size,
        profile.st_mtime_ns,
        profile.st_ctime_ns,
    )


@contextmanager
def open_regular_binary(path: Path, label: str) -> Iterator[tuple[BinaryIO, tuple[int, ...]]]:
    """Open a stable regular file descriptor while rejecting parent symlinks."""

    if not path.is_absolute() or ".." in path.parts or len(path.parts) < 2:
        raise ValidationError(f"{label} path is not normalized and absolute: {path}")
    directory_flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
    file_flags = os.O_RDONLY
    if hasattr(os, "O_CLOEXEC"):
        directory_flags |= os.O_CLOEXEC
        file_flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        directory_flags |= os.O_NOFOLLOW
        file_flags |= os.O_NOFOLLOW
    directory_descriptor = -1
    descriptor = -1
    try:
        directory_descriptor = os.open("/", directory_flags)
        for component in path.parts[1:-1]:
            next_descriptor = os.open(
                component, directory_flags, dir_fd=directory_descriptor
            )
            os.close(directory_descriptor)
            directory_descriptor = next_descriptor
        descriptor = os.open(path.name, file_flags, dir_fd=directory_descriptor)
        opened = os.fstat(descriptor)
        if not stat.S_ISREG(opened.st_mode):
            raise ValidationError(f"{label} is not a regular file: {path}")
        opened_profile = stat_profile(opened)
        with os.fdopen(descriptor, "rb", closefd=False) as stream:
            yield stream, opened_profile
        after_profile = stat_profile(os.fstat(descriptor))
        try:
            named = os.stat(path.name, dir_fd=directory_descriptor, follow_symlinks=False)
        except OSError as error:
            raise ValidationError(f"{label} changed while open: {path}") from error
        if (
            after_profile != opened_profile
            or stat_profile(named) != opened_profile
            or not stat.S_ISREG(named.st_mode)
        ):
            raise ValidationError(f"{label} changed while open: {path}")
    except OSError as error:
        raise ValidationError(
            f"cannot securely open {label}; a component may be a symlink: {path}"
        ) from error
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        if directory_descriptor >= 0:
            os.close(directory_descriptor)


def require_directory_path(path: Path, label: str) -> None:
    """Require an existing directory path with no symlink component."""

    if not path.is_absolute() or ".." in path.parts:
        raise ValidationError(f"{label} path is not normalized and absolute: {path}")
    current = Path("/")
    for component in path.parts[1:]:
        current /= component
        try:
            profile = current.lstat()
        except OSError as error:
            raise ValidationError(f"{label} is unavailable: {path}") from error
        if stat.S_ISLNK(profile.st_mode):
            raise ValidationError(f"{label} has a symlink component: {current}")
    if not stat.S_ISDIR(path.lstat().st_mode):
        raise ValidationError(f"{label} is not a directory: {path}")


def require_tree_without_symlinks(root: Path, label: str) -> None:
    """Reject every symlink in or above one retained output tree."""

    require_directory_path(root, label)
    for directory, names, files in os.walk(root, followlinks=False):
        directory_path = Path(directory)
        for name in [*names, *files]:
            path = directory_path / name
            try:
                profile = path.lstat()
            except OSError as error:
                raise ValidationError(f"{label} entry is unavailable: {path}") from error
            if stat.S_ISLNK(profile.st_mode):
                raise ValidationError(f"{label} contains a symlink: {path}")


def read_regular_bytes(
    path: Path,
    label: str,
    profiles: dict[Path, tuple[int, ...]] | None = None,
) -> bytes:
    """Read one stable regular file and optionally remember its version."""

    with open_regular_binary(path, label) as (stream, profile):
        payload = stream.read()
    if profiles is not None:
        remember_profile(profiles, path, profile, label)
    return payload


def sha256_regular_file(
    path: Path,
    label: str,
    profiles: dict[Path, tuple[int, ...]] | None = None,
) -> tuple[str, int]:
    """Hash one stable regular file and optionally remember its version."""

    digest = hashlib.sha256()
    with open_regular_binary(path, label) as (stream, profile):
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    if profiles is not None:
        remember_profile(profiles, path, profile, label)
    return digest.hexdigest(), profile[6]


def load_json(
    path: Path,
    label: str,
    profiles: dict[Path, tuple[int, ...]] | None = None,
) -> tuple[dict[str, object], str]:
    """Load one stable JSON object and return its byte digest."""

    payload = read_regular_bytes(path, label, profiles)

    def reject_constant(value: str) -> object:
        raise ValidationError(f"{label} contains non-finite JSON value {value}: {path}")

    def unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
        value: dict[str, object] = {}
        for key, item in pairs:
            if key in value:
                raise ValidationError(f"{label} repeats JSON key {key}: {path}")
            value[key] = item
        return value

    try:
        value = json.loads(
            payload,
            object_pairs_hook=unique_object,
            parse_constant=reject_constant,
        )
    except ValidationError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValidationError(f"{label} is not valid JSON: {path}") from error
    if not isinstance(value, dict):
        raise ValidationError(f"{label} must contain a JSON object: {path}")

    pending: list[object] = [value]
    while pending:
        item = pending.pop()
        if isinstance(item, dict):
            pending.extend(item.values())
        elif isinstance(item, list):
            pending.extend(item)
        elif isinstance(item, float) and not math.isfinite(item):
            raise ValidationError(f"{label} contains a non-finite JSON number: {path}")
    return value, hashlib.sha256(payload).hexdigest()


def recheck_profiles(profiles: dict[Path, tuple[int, ...]]) -> None:
    """Require every authenticated file to remain the same named version."""

    for path, expected in profiles.items():
        with open_regular_binary(path, "authenticated retained file") as (_, actual):
            if actual != expected:
                raise ValidationError(
                    f"authenticated retained file changed during validation: {path}"
                )


def remember_profile(
    profiles: dict[Path, tuple[int, ...]],
    path: Path,
    profile: tuple[int, ...],
    label: str,
) -> None:
    """Remember the first authenticated version and reject cross-read races."""

    retained = profiles.get(path)
    if retained is not None and retained != profile:
        raise ValidationError(f"{label} changed between authenticated reads: {path}")
    profiles[path] = profile


def require_sha256(value: object, label: str) -> str:
    """Require one lowercase SHA-256 digest."""

    if not isinstance(value, str) or SHA256_PATTERN.fullmatch(value) is None:
        raise ValidationError(f"{label} is not a SHA-256 digest")
    return value


def require_revision(value: object, label: str) -> str:
    """Require one full lowercase Git revision."""

    if not isinstance(value, str) or GIT_REVISION_PATTERN.fullmatch(value) is None:
        raise ValidationError(f"{label} is not a full Git revision")
    return value


def require_job_id(value: object, label: str) -> str:
    """Require one numeric scheduler job identifier."""

    if not isinstance(value, str) or JOB_ID_PATTERN.fullmatch(value) is None:
        raise ValidationError(f"{label} is not a numeric job ID")
    return value


def require_segment_name(value: object, label: str) -> str:
    """Require one bounded path-safe segment name."""

    if not isinstance(value, str) or SEGMENT_PATTERN.fullmatch(value) is None:
        raise ValidationError(f"{label} is invalid")
    return value


def require_qualification_approval(value: object) -> dict[str, object]:
    """Validate the retained E03 token using the controller's approval contract."""

    if not isinstance(value, dict):
        raise ValidationError("qualification approval token is invalid")
    if value.get("schema_version") != 1:
        raise ValidationError("qualification approval token has wrong schema")
    if value.get("execution_epoch") != EXECUTION_EPOCH:
        raise ValidationError("qualification approval token has wrong execution epoch")
    require_revision(
        value.get("approved_executable_revision"),
        "qualification approval executable revision",
    )
    require_sha256(
        value.get("approved_executable_sha256"),
        "qualification approval executable digest",
    )
    for key in ("approved_utc", "approved_by", "review_notes"):
        if not isinstance(value.get(key), str) or not str(value[key]).strip():
            raise ValidationError(f"qualification approval token lacks nonempty {key}")
    return value


def git_run(
    arguments: list[str],
    *,
    capture_output: bool = True,
    timeout: int = 120,
) -> subprocess.CompletedProcess:
    """Run the authenticated absolute Git binary with isolated configuration."""

    environment = {
        key: value for key, value in os.environ.items() if not key.startswith("GIT_")
    }
    environment.update({
        "GIT_CONFIG_NOSYSTEM": "1",
        "GIT_CONFIG_GLOBAL": os.devnull,
        "GIT_CONFIG_SYSTEM": os.devnull,
        "GIT_TERMINAL_PROMPT": "0",
        "LC_ALL": "C",
    })
    with open_regular_binary(
        GIT_EXECUTABLE, "authenticated absolute Git executable"
    ) as (stream, profile):
        if profile[2] & 0o111 == 0:
            raise ValidationError("authenticated absolute Git executable is not executable")
        descriptor = stream.fileno()
        try:
            return subprocess.run(
                [str(GIT_EXECUTABLE), "--no-replace-objects", *arguments],
                executable=f"/proc/self/fd/{descriptor}",
                pass_fds=(descriptor,),
                check=False,
                capture_output=capture_output,
                env=environment,
                timeout=timeout,
            )
        except OSError as error:
            raise ValidationError(
                "authenticated absolute Git executable could not be run"
            ) from error


def require_within(path: Path, root: Path, label: str) -> None:
    """Require one normalized path to be lexically inside a retained root."""

    try:
        path.relative_to(root)
    except ValueError as error:
        raise ValidationError(f"{label} is outside its retained root: {path}") from error


def authenticated_file(
    path: Path,
    expected_sha256: object,
    label: str,
    profiles: dict[Path, tuple[int, ...]],
    expected_size: object | None = None,
) -> dict[str, object]:
    """Authenticate one retained file against its digest and optional size."""

    expected_digest = require_sha256(expected_sha256, f"{label} digest")
    digest, size = sha256_regular_file(path, label, profiles)
    if digest != expected_digest:
        raise ValidationError(f"{label} SHA-256 differs from retained provenance")
    if expected_size is not None:
        retained_size = require_exact_int(expected_size, f"{label} retained size")
        if retained_size != size:
            raise ValidationError(f"{label} size differs from retained provenance")
    return {"path": str(path), "size_bytes": size, "sha256": digest}


def validator_implementation_provenance(
    path: Path,
    profiles: dict[Path, tuple[int, ...]],
) -> dict[str, object]:
    """Bind the running validator bytes to one stable committed Git revision."""

    digest, size = sha256_regular_file(path, "validator source", profiles)
    try:
        source = path.resolve(strict=True)
        resolved_digest, resolved_size = sha256_regular_file(
            source, "resolved validator source", profiles
        )
        if resolved_digest != digest or resolved_size != size:
            raise ValidationError("resolved validator source differs from invoked bytes")
        root_query = git_run([
            "-C",
            str(source.parent),
            "rev-parse",
            "--show-toplevel",
        ])
    except (OSError, subprocess.TimeoutExpired) as error:
        raise ValidationError("cannot determine validator source repository") from error
    if root_query.returncode != 0:
        raise ValidationError("cannot determine validator source repository")
    try:
        repository = declared_path(
            root_query.stdout.decode("utf-8").strip(), "validator source repository"
        )
        relative = source.relative_to(repository).as_posix()
    except (UnicodeDecodeError, ValueError) as error:
        raise ValidationError("validator source is outside its Git repository") from error
    require_directory_path(repository, "validator source repository")
    if relative != VALIDATOR_REPOSITORY_PATH:
        raise ValidationError("validator source repository path is inconsistent")
    try:
        head_before = git_run(["-C", str(repository), "rev-parse", "--verify", "HEAD"])
        tracked = git_run([
            "-C",
            str(repository),
            "ls-files",
            "--error-unmatch",
            "--",
            relative,
        ])
        unstaged = git_run([
            "-C",
            str(repository),
            "diff",
            "--quiet",
            "--",
            relative,
        ])
        staged = git_run([
            "-C",
            str(repository),
            "diff",
            "--cached",
            "--quiet",
            "--",
            relative,
        ])
        retained = git_run([
            "-C",
            str(repository),
            "show",
            f"HEAD:{relative}",
        ])
        head_after = git_run(["-C", str(repository), "rev-parse", "--verify", "HEAD"])
    except (OSError, subprocess.TimeoutExpired) as error:
        raise ValidationError("cannot authenticate committed validator source") from error
    if tracked.returncode != 0:
        raise ValidationError("validator source must be tracked by Git")
    if unstaged.returncode != 0 or staged.returncode != 0:
        raise ValidationError("validator source must be committed before validation")
    if (
        head_before.returncode != 0
        or head_after.returncode != 0
        or head_before.stdout != head_after.stdout
    ):
        raise ValidationError("validator source revision changed during authentication")
    revision = require_revision(
        head_before.stdout.decode("utf-8").strip(), "validator source revision"
    )
    if retained.returncode != 0 or hashlib.sha256(retained.stdout).hexdigest() != digest:
        raise ValidationError("validator committed revision bytes differ from live bytes")
    final_digest, final_size = sha256_regular_file(
        source, "resolved validator source", profiles
    )
    if final_digest != digest or final_size != size:
        raise ValidationError("validator source changed during authentication")
    return {
        "path": str(source),
        "repository_path": relative,
        "revision": revision,
        "sha256": digest,
        "size_bytes": size,
        "committed": True,
    }


def captured_file(
    path: Path,
    output_root: Path,
    label: str,
    profiles: dict[Path, tuple[int, ...]],
) -> dict[str, object]:
    """Capture one retained output file through a stable streaming hash."""

    require_within(path, output_root, label)
    digest, size = sha256_regular_file(path, label, profiles)
    return {"path": str(path), "size_bytes": size, "sha256": digest}


def parse_athinput(payload: bytes, label: str) -> dict[str, dict[str, str]]:
    """Parse the block/key structure needed from one Athena input deck."""

    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValidationError(f"{label} is not UTF-8 text") from error
    blocks: dict[str, dict[str, str]] = {}
    block: str | None = None
    for original in text.splitlines():
        line = original.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            name = line[1:-1].strip()
            if not name or name == "par_end":
                block = None
                continue
            if name in blocks:
                raise ValidationError(f"{label} repeats block <{name}>")
            blocks[name] = {}
            block = name
            continue
        if block is None or "=" not in line:
            raise ValidationError(f"{label} contains malformed parameter text")
        key, value = (token.strip() for token in line.split("=", 1))
        if not key or not value or key in blocks[block]:
            raise ValidationError(f"{label} repeats or omits {block}/{key}")
        blocks[block][key] = value
    return blocks


def parameter(
    blocks: dict[str, dict[str, str]],
    block: str,
    key: str,
    label: str,
    default: str | None = None,
) -> str:
    """Read one required or explicitly defaulted Athena parameter."""

    value = blocks.get(block, {}).get(key, default)
    if value is None:
        raise ValidationError(f"{label} lacks required parameter {block}/{key}")
    return value


def bool_parameter(
    blocks: dict[str, dict[str, str]],
    block: str,
    key: str,
    label: str,
    default: str | None = None,
) -> bool:
    """Read one exact Athena boolean."""

    value = parameter(blocks, block, key, label, default)
    if value not in {"true", "false"}:
        raise ValidationError(f"{label} parameter {block}/{key} is not boolean")
    return value == "true"


def positive_int_parameter(
    blocks: dict[str, dict[str, str]], block: str, key: str, label: str
) -> int:
    """Read one positive Athena integer."""

    try:
        value = int(parameter(blocks, block, key, label))
    except ValueError as error:
        raise ValidationError(f"{label} parameter {block}/{key} is not an integer") from error
    if value <= 0:
        raise ValidationError(f"{label} parameter {block}/{key} is not positive")
    return value


def equivalent_parameter_value(actual: str, expected: str) -> bool:
    """Compare one frozen Athena parameter while tolerating numeric formatting."""

    if actual == expected:
        return True
    try:
        left = float(actual)
        right = float(expected)
    except ValueError:
        return False
    return math.isfinite(left) and math.isfinite(right) and left == right


def frozen_parameters_for_case(
    case_id: str, case: dict[str, object]
) -> tuple[dict[tuple[str, str], str], set[tuple[str, str]]]:
    """Return the closed semantic input contract for one frozen Stage I case."""

    try:
        resolution, beta, passive, forcing, lf_k, limiter_rate = CASE_INPUT_CONTRACTS[
            case_id
        ]
    except KeyError as error:
        raise ValidationError(f"no frozen input contract exists for {case_id}") from error
    try:
        mesh = tuple(int(value) for value in str(resolution).split("x"))
    except ValueError as error:
        raise ValidationError(f"frozen input contract resolution is invalid for {case_id}") from error
    if len(mesh) != 3 or any(value <= 0 for value in mesh):
        raise ValidationError(f"frozen input contract resolution is invalid for {case_id}")
    expected = dict(COMMON_FROZEN_PARAMETERS)
    expected.update({
        ("job", "basename"): Path(CASE_INPUT_PATHS[case_id]).stem,
        ("mesh", "nx1"): str(mesh[0]),
        ("mesh", "nx2"): str(mesh[1]),
        ("mesh", "nx3"): str(mesh[2]),
        ("mhd", "passive"): str(passive).lower(),
        ("mhd", "lf_k_parallel"): format(lf_k, ".17g"),
        ("mhd", "limiter_nu_coll"): format(limiter_rate, ".17g"),
        ("problem", "passive_delta"): str(passive).lower(),
        ("problem", "beta0"): format(beta, ".17g"),
        ("problem", "p_parallel0"): format(beta / 2.0, ".17g"),
        ("problem", "p_perp0"): format(beta / 2.0, ".17g"),
    })
    forbidden: set[tuple[str, str]] = set()
    if passive:
        expected[("mhd", "iso_sound_speed")] = format(math.sqrt(beta / 2.0), ".17g")
    else:
        forbidden.add(("mhd", "iso_sound_speed"))
    if forcing == "alfvenic":
        expected.update({
            ("turb_driving", "driving_type"): "1",
            ("turb_driving", "projection_policy"): "mks24_alfvenic_perpendicular",
            ("turb_driving", "sol_fraction"): "1.0",
            ("turb_driving", "rseed"): "271828",
            ("turb_driving", "tcorr"): "2.0",
        })
    else:
        expected.update({
            ("turb_driving", "driving_type"): "0",
            ("turb_driving", "projection_policy"): "mks24_random_unprojected",
            ("turb_driving", "rseed"): "314159",
            ("turb_driving", "tcorr"): "0.2" if forcing == "random_sonic" else "2.0",
        })
        forbidden.add(("turb_driving", "sol_fraction"))
    if case_id in FINITE_LIMITER_CASE_IDS:
        forbidden.add(("mhd", "limiter_hardwall"))
    else:
        expected[("mhd", "limiter_hardwall")] = "true"
    expected_input = CASE_INPUT_PATHS[case_id]
    if (
        case.get("name") != CASE_NAMES[case_id]
        or case.get("input") != expected_input
        or case.get("resolution") != resolution
    ):
        raise ValidationError(f"retained matrix does not match frozen {case_id} mapping")
    return expected, forbidden


def require_frozen_parameters(
    blocks: dict[str, dict[str, str]],
    expected: dict[tuple[str, str], str],
    forbidden: set[tuple[str, str]],
    label: str,
) -> None:
    """Require one parameter dump to preserve the frozen scientific contract."""

    for (block, key), value in expected.items():
        actual = parameter(blocks, block, key, label)
        if not equivalent_parameter_value(actual, value):
            raise ValidationError(
                f"{label} violates frozen {block}/{key}={value}"
            )
    for block, key in forbidden:
        if key in blocks.get(block, {}):
            raise ValidationError(f"{label} contains forbidden parameter {block}/{key}")


def flattened_parameters(
    blocks: dict[str, dict[str, str]],
) -> dict[tuple[str, str], str]:
    """Return one parameter tree as an exact block/key mapping."""

    return {
        (block, key): value
        for block, parameters in blocks.items()
        for key, value in parameters.items()
    }


def require_exact_archived_input_parameters(
    blocks: dict[str, dict[str, str]],
    expected: dict[tuple[str, str], str],
    label: str,
) -> None:
    """Reject every unqualified archived-input parameter and missing dynamic field."""

    actual = flattened_parameters(blocks)
    required = set(expected) | set(ARCHIVED_INPUT_DYNAMIC_PARAMETERS)
    allowed = required | set(ARCHIVED_INPUT_OPTIONAL_DEFAULTS)
    if set(actual) - allowed:
        raise ValidationError(
            f"{label} contains unqualified parameters: {sorted(set(actual) - allowed)}"
        )
    if required - set(actual):
        raise ValidationError(
            f"{label} lacks required parameters: {sorted(required - set(actual))}"
        )
    for key, value in ARCHIVED_INPUT_OPTIONAL_DEFAULTS.items():
        if key in actual and not equivalent_parameter_value(actual[key], value):
            raise ValidationError(
                f"{label} violates qualified optional default {key[0]}/{key[1]}={value}"
            )
    try:
        target = float(actual[("time", "tlim")])
    except ValueError as error:
        raise ValidationError(f"{label} time/tlim is not numeric") from error
    if not math.isfinite(target) or target <= 0.0:
        raise ValidationError(f"{label} time/tlim is not positive and finite")


def require_input_contract(
    blocks: dict[str, dict[str, str]], case_id: str, case: dict[str, object]
) -> dict[str, object]:
    """Bind the archived deck to the closed Stage I CGL-LF physics contract."""

    label = "archived input"
    frozen_parameters, forbidden_parameters = frozen_parameters_for_case(case_id, case)
    require_frozen_parameters(blocks, frozen_parameters, forbidden_parameters, label)
    require_exact_archived_input_parameters(blocks, frozen_parameters, label)
    passive = bool_parameter(blocks, "mhd", "passive", label)
    passive_delta = bool_parameter(blocks, "problem", "passive_delta", label)
    if passive != passive_delta or passive != (case_id in PASSIVE_CASE_IDS):
        raise ValidationError("archived input passive-Delta mode differs from mapped case")
    required_values = {
        ("mhd", "eos"): "cgl",
        ("mhd", "cgl_heat_flux"): "landau_fluid",
        ("problem", "pgen_name"): "cgl_lf_paper",
    }
    for (block, key), expected in required_values.items():
        if parameter(blocks, block, key, label) != expected:
            raise ValidationError(f"archived input violates {block}/{key}={expected}")
    for block, key in (
        ("mhd", "cgl_lf_strict_admissibility"),
        ("mhd", "cgl_lf_record_pressure_work"),
        ("turb_driving", "record_injected_work"),
    ):
        if not bool_parameter(blocks, block, key, label):
            raise ValidationError(f"archived input requires {block}/{key}=true")
    for unsupported in ("hydro", "radiation", "z4c", "adm", "particles"):
        if unsupported in blocks:
            raise ValidationError(
                f"archived input contains unsupported Stage I block <{unsupported}>"
            )
    try:
        nscalars = int(parameter(blocks, "mhd", "nscalars", label, "0"))
    except ValueError as error:
        raise ValidationError("archived input mhd/nscalars is invalid") from error
    if nscalars != 0:
        raise ValidationError("qualified Stage I restart parser requires mhd/nscalars=0")
    mesh = tuple(
        positive_int_parameter(blocks, "mesh", f"nx{axis}", label)
        for axis in range(1, 4)
    )
    meshblock = tuple(
        positive_int_parameter(blocks, "meshblock", f"nx{axis}", label)
        for axis in range(1, 4)
    )
    if any(total % block != 0 for total, block in zip(mesh, meshblock)):
        raise ValidationError("archived input mesh is not exactly tiled by meshblocks")
    divisions = tuple(total // block for total, block in zip(mesh, meshblock))
    expected_meshblocks = math.prod(divisions)
    nghost = positive_int_parameter(blocks, "mesh", "nghost", label)
    hardwall = bool_parameter(blocks, "mhd", "limiter_hardwall", label, "false")
    finite_limiter = case_id in FINITE_LIMITER_CASE_IDS
    if hardwall == finite_limiter:
        raise ValidationError("archived input limiter mode differs from mapped case")
    storage = {}
    for output, kind in (("output2", "snapshot"), ("output3", "restart")):
        expected_type = "bin" if kind == "snapshot" else "rst"
        if parameter(blocks, output, "file_type", label) != expected_type:
            raise ValidationError(f"archived input {output} is not a {expected_type} output")
        storage[kind] = (
            "per_rank"
            if bool_parameter(blocks, output, "single_file_per_rank", label)
            else "shared_mpiio"
        )
    nout = tuple(value + 2 * nghost for value in meshblock)
    restart_data_size = (
        math.prod(nout) * 6 * 8
        + (nout[0] + 1) * nout[1] * nout[2] * 8
        + nout[0] * (nout[1] + 1) * nout[2] * 8
        + nout[0] * nout[1] * (nout[2] + 1) * 8
    )
    try:
        domain = tuple(
            (
                float(parameter(blocks, "mesh", f"x{axis}min", label)),
                float(parameter(blocks, "mesh", f"x{axis}max", label)),
            )
            for axis in range(1, 4)
        )
    except ValueError as error:
        raise ValidationError("archived input mesh geometry is invalid") from error
    if any(
        not math.isfinite(lower) or not math.isfinite(upper) or upper <= lower
        for lower, upper in domain
    ):
        raise ValidationError("archived input mesh geometry is invalid")
    return {
        "passive": passive,
        "hardwall": hardwall,
        "finite_limiter": finite_limiter,
        "mesh": mesh,
        "meshblock": meshblock,
        "domain": domain,
        "divisions": divisions,
        "nghost": nghost,
        "expected_meshblocks": expected_meshblocks,
        "restart_data_size": restart_data_size,
        "storage": storage,
        "history_cadence": float(parameter(blocks, "output1", "dt", label)),
        "snapshot_cadence": float(parameter(blocks, "output2", "dt", label)),
        "restart_cadence": float(parameter(blocks, "output3", "dt", label)),
        "frozen_parameters": frozen_parameters,
        "forbidden_parameters": forbidden_parameters,
        "archived_parameter_blocks": {
            block: dict(parameters) for block, parameters in blocks.items()
        },
    }


def expected_turbulence_restart_configuration(
    input_contract: dict[str, object],
) -> dict[str, int | float]:
    """Derive the exact qualified TurbulenceDriver restart configuration."""

    blocks = input_contract.get("archived_parameter_blocks")
    domain = input_contract.get("domain")
    if not isinstance(blocks, dict) or not isinstance(domain, tuple) or len(domain) != 3:
        raise ValidationError("qualified turbulence restart contract is invalid")
    label = "qualified turbulence restart contract"

    def integer(key: str, default: str) -> int:
        try:
            return int(parameter(blocks, "turb_driving", key, label, default))
        except ValueError as error:
            raise ValidationError(
                f"qualified turbulence restart integer {key} is invalid"
            ) from error

    def real(key: str, default: str) -> float:
        try:
            value = float(parameter(blocks, "turb_driving", key, label, default))
        except ValueError as error:
            raise ValidationError(
                f"qualified turbulence restart real {key} is invalid"
            ) from error
        if not math.isfinite(value):
            raise ValidationError(
                f"qualified turbulence restart real {key} is invalid"
            )
        return value

    def boolean(key: str, default: str) -> int:
        return int(bool_parameter(blocks, "turb_driving", key, label, default))

    def choice(key: str, default: str, choices: dict[str, int]) -> int:
        value = parameter(blocks, "turb_driving", key, label, default)
        try:
            return choices[value]
        except KeyError as error:
            raise ValidationError(
                f"qualified turbulence restart choice {key} is invalid"
            ) from error

    nlow = integer("nlow", "1")
    nhigh = integer("nhigh", "3")
    driving_type = integer("driving_type", "0")
    min_kx = integer("min_kx", "0")
    max_kx = integer("max_kx", str(nhigh))
    min_ky = integer("min_ky", "0")
    max_ky = integer("max_ky", str(nhigh))
    min_kz = integer("min_kz", "0")
    max_kz = integer("max_kz", str(nhigh))
    tile_nx = integer("tile_nx", "1")
    tile_ny = integer("tile_ny", "1")
    tile_nz = integer("tile_nz", "1")
    physical_k_shell = boolean("physical_k_shell", "false")
    k_shell_unit = real("k_shell_unit", "0")
    if min(nlow, tile_nx, tile_ny, tile_nz) <= 0 or nhigh < nlow:
        raise ValidationError("qualified turbulence restart integer contract is invalid")
    try:
        lengths = tuple(float(pair[1]) - float(pair[0]) for pair in domain)
    except (IndexError, TypeError, ValueError) as error:
        raise ValidationError("qualified turbulence restart domain is invalid") from error
    if any(not math.isfinite(length) or length <= 0.0 for length in lengths):
        raise ValidationError("qualified turbulence restart domain is invalid")
    tile_lengths = tuple(
        length / count
        for length, count in zip(lengths, (tile_nx, tile_ny, tile_nz))
    )
    wave_units = tuple(2.0 * math.pi / length for length in tile_lengths)

    def driven_mode(nkx: int, nky: int, nkz: int) -> bool:
        if physical_k_shell:
            normalized_k2 = math.fsum((
                (wave_units[0] * nkx) ** 2,
                (wave_units[1] * nky) ** 2,
                (wave_units[2] * nkz) ** 2,
            )) / (k_shell_unit ** 2)
            return nlow ** 2 <= normalized_k2 <= nhigh ** 2
        if driving_type == 0:
            nsqr = nkx ** 2 + nky ** 2 + nkz ** 2
            return nlow ** 2 <= nsqr <= nhigh ** 2
        nperp_sqr = nkx ** 2 + nky ** 2
        return (
            nlow ** 2 <= nperp_sqr <= nhigh ** 2
            and nlow ** 2 <= nkz ** 2 <= nhigh ** 2
        )

    mode_count = sum(
        1
        for nkx in range(min_kx, max_kx + 1)
        for nky in range(min_ky, max_ky + 1)
        for nkz in range(min_kz, max_kz + 1)
        if (nkx != 0 or nky != 0 or nkz != 0) and driven_mode(nkx, nky, nkz)
    )
    if mode_count <= 0:
        raise ValidationError("qualified turbulence restart mode count is invalid")
    use_npeak = "npeak" in blocks.get("turb_driving", {})
    npeak = real("npeak", "0") if use_npeak else 0.0
    kpeak = (
        npeak * 2.0 * math.pi / tile_lengths[0]
        if use_npeak
        else real("kpeak", str(QUALIFIED_PRODUCT_PARAMETER_ADDITIONS[
            ("turb_driving", "kpeak")
        ]))
    )
    normalization = choice("normalization", "edot", {"edot": 0, "accel_rms": 1})
    dedt = real("dedt", "0") if normalization == 0 else 0.0
    accel_rms = real("accel_rms", "0") if normalization == 1 else 0.0
    tcorr = real("tcorr", "0")
    turb_flag = integer("turb_flag", "2")
    tdriv_duration = (
        real("tdriv_duration", format(tcorr, ".17g"))
        if turb_flag == 1
        else struct.unpack("<f", struct.pack("<I", 0x7F7FFFFF))[0]
    )
    configuration: dict[str, int | float] = {
        "version": 3,
        "mode_count": mode_count,
        "nlow": nlow,
        "nhigh": nhigh,
        "driving_type": driving_type,
        "min_kx": min_kx,
        "max_kx": max_kx,
        "min_ky": min_ky,
        "max_ky": max_ky,
        "min_kz": min_kz,
        "max_kz": max_kz,
        "use_npeak": int(use_npeak),
        "turb_flag": turb_flag,
        "tile_nx": tile_nx,
        "tile_ny": tile_ny,
        "tile_nz": tile_nz,
        "normalization": normalization,
        "localization": choice(
            "localization", "none", {"none": 0, "include": 1, "exclude": 2}
        ),
        "spectrum": choice("spectrum", "parabolic", {"parabolic": 0, "power_law": 1}),
        "projection_policy": choice(
            "projection_policy",
            "solenoidal_compressive",
            {
                "solenoidal_compressive": 0,
                "mks24_random_unprojected": 1,
                "mks24_alfvenic_perpendicular": 2,
            },
        ),
        "physical_k_shell": physical_k_shell,
        "isotropic_power_spectrum": boolean("isotropic_power_spectrum", "false"),
        "record_injected_work": boolean("record_injected_work", "false"),
        "tcorr": tcorr,
        "dt_update": real("dt_update", "0.01"),
        "dedt": dedt,
        "accel_rms": accel_rms,
        "sol_fraction": real("sol_fraction", "1.0"),
        "kpeak": kpeak,
        "npeak": npeak,
        "expo": real("expo", format(5.0 / 3.0, ".17g")),
        "exp_prp": real("exp_prp", "1.66667"),
        "exp_prl": real("exp_prl", "0"),
        "tdriv_duration": tdriv_duration,
        "tdriv_start": real("tdriv_start", "0"),
        "sigma_x1": real("sigma_x1", "-1"),
        "sigma_x2": real("sigma_x2", "-1"),
        "sigma_x3": real("sigma_x3", "-1"),
        "center_x1": real("center_x1", "0"),
        "center_x2": real("center_x2", "0"),
        "center_x3": real("center_x3", "0"),
        "k_shell_unit": k_shell_unit,
    }
    if set(configuration) != set(TURBULENCE_METADATA_CONFIGURATION_FIELDS):
        raise ValidationError("qualified turbulence restart configuration is incomplete")
    return configuration


def normalized_batch_script_sha256(payload: bytes) -> str:
    """Return the normalized self-digest of one retained batch script."""

    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValidationError("prepared batch script is not UTF-8 text") from error
    matches = BATCH_SCRIPT_DIGEST_PATTERN.findall(text)
    if len(matches) != 1:
        raise ValidationError("prepared batch script lacks one self-digest")
    normalized = BATCH_SCRIPT_DIGEST_PATTERN.sub(
        f"BATCH_SCRIPT_SHA256={BATCH_SCRIPT_DIGEST_PLACEHOLDER}", text
    )
    return hashlib.sha256(normalized.encode("utf-8")).hexdigest()


def quote(value: str | Path) -> str:
    """Quote one literal exactly as the production batch generator does."""

    return shlex.quote(str(value))


def expected_job_name(manifest: dict[str, object]) -> str:
    """Return the exact safe production Slurm job name."""

    run = manifest["run"]
    segment = str(run["segment"])
    if SEGMENT_PATTERN.fullmatch(segment) is None:
        raise ValidationError("manifest segment name is invalid")
    job_name = f"cgl_mks24_{EXECUTION_EPOCH_SLUG}_{run['case_id']}_{segment}"
    job_name = re.sub(r"[^A-Za-z0-9_]+", "_", job_name)[:60]
    if re.fullmatch(r"[A-Za-z0-9_]+", job_name) is None:
        raise ValidationError("generated Stage I job name is unsafe")
    return job_name


def validate_restart_launch_intent(
    command: dict[str, object],
    manifest_path: Path,
    storage: str,
    rank_count: int,
) -> None:
    """Validate every continuation field before launch-script regeneration."""

    records = command.get("restart_files")
    if not isinstance(records, list):
        raise ValidationError("manifest restart launch inventory is invalid")
    archive_root = manifest_path.parent / "submitted_restart"
    validated: list[tuple[Path, str]] = []
    for record in records:
        if not isinstance(record, dict):
            raise ValidationError("manifest restart launch inventory is invalid")
        path = declared_path(record.get("path"), "restart launch sibling")
        require_within(path, archive_root, "restart launch sibling")
        digest = require_sha256(record.get("sha256"), "restart launch sibling digest")
        size = require_exact_int(record.get("size_bytes"), "restart launch sibling size")
        if size <= 0 or any(path == item[0] for item in validated):
            raise ValidationError("manifest restart launch inventory is invalid")
        validated.append((path, digest))
    if validated:
        paths = [item[0] for item in validated]
        if storage == "per_rank":
            expected_parents = [
                archive_root / f"rank_{rank:08d}" for rank in range(rank_count)
            ]
            if (
                len(paths) != rank_count
                or [path.parent for path in paths] != expected_parents
                or any(path.name != paths[0].name for path in paths[1:])
            ):
                raise ValidationError("manifest restart launch inventory is invalid")
        elif storage == "shared_mpiio":
            if len(paths) != 1 or paths[0].parent != archive_root:
                raise ValidationError("manifest restart launch inventory is invalid")
        else:
            raise ValidationError("manifest restart launch storage is invalid")
        restart_file = declared_path(
            command.get("restart_file"), "restart launch representative"
        )
        source_restart = declared_path(
            command.get("source_restart_file"), "source restart representative"
        )
        require_sha256(
            command.get("restart_sha256"), "restart launch representative digest"
        )
        if (
            restart_file != validated[0][0]
            or command.get("restart_sha256") != validated[0][1]
            or not isinstance(command.get("parent_segment"), dict)
            or source_restart == restart_file
        ):
            raise ValidationError("manifest continuation launch intent is inconsistent")
    elif any(
        command.get(key) is not None
        for key in (
            "restart_file",
            "restart_sha256",
            "source_restart_file",
            "parent_segment",
        )
    ):
        raise ValidationError("manifest fresh launch intent contains restart state")


def validate_launch_intent_fields(
    manifest: dict[str, object],
    manifest_path: Path,
    root: Path,
    case: dict[str, object],
    restart_storage: str,
) -> None:
    """Validate every manifest field interpolated into the production script."""

    run = manifest.get("run")
    allocation = manifest.get("allocation")
    command = manifest.get("command")
    paths = manifest.get("paths")
    if not all(isinstance(item, dict) for item in (run, allocation, command, paths)):
        raise ValidationError("manifest lacks exact production launch dictionaries")
    case_id = str(run.get("case_id"))
    segment = str(run.get("segment"))
    expected_basename = f"{EXECUTION_EPOCH_SLUG}_{CASE_NAMES[case_id]}_{segment}"
    if (
        run.get("case_name") != CASE_NAMES[case_id]
        or run.get("resolution") != case.get("resolution")
        or run.get("run_basename") != expected_basename
        or re.fullmatch(r"[A-Za-z0-9_]+", expected_basename) is None
    ):
        raise ValidationError("manifest run identity is not exact production launch intent")

    nodes = require_exact_int(allocation.get("nodes"), "allocation nodes")
    ranks_per_node = require_exact_int(
        allocation.get("ranks_per_node"), "allocation ranks_per_node"
    )
    cpus_per_task = require_exact_int(
        allocation.get("cpus_per_task"), "allocation cpus_per_task"
    )
    requested_seconds = parse_walltime(
        allocation.get("requested_walltime"), "allocation requested_walltime"
    )
    retained_seconds = require_exact_int(
        allocation.get("requested_seconds"), "allocation requested_seconds"
    )
    reserved_node_hours = allocation.get("reserved_node_hours")
    if (
        nodes not in CASE_NODE_PROFILES[case_id]
        or ranks_per_node != EXPECTED_RANKS_PER_NODE
        or cpus_per_task != EXPECTED_CPUS_PER_TASK
        or requested_seconds <= 0
        or requested_seconds > MAX_SEGMENT_SECONDS
        or retained_seconds != requested_seconds
        or isinstance(reserved_node_hours, bool)
        or not isinstance(reserved_node_hours, (int, float))
        or not math.isfinite(float(reserved_node_hours))
        or abs(float(reserved_node_hours) - nodes * requested_seconds / 3600.0)
        > 5.0e-12
    ):
        raise ValidationError("manifest allocation is not exact production launch intent")
    athena_seconds = parse_walltime(
        command.get("athena_walltime"), "command athena_walltime"
    )
    if (
        athena_seconds <= 0
        or athena_seconds > requested_seconds - 600
    ):
        raise ValidationError("manifest Athena walltime violates production launch policy")

    expected_paths = {
        "run_dir": manifest_path.parents[1],
        "output_dir": manifest_path.parents[1] / "output",
        "environment_log": manifest_path.parent / "run_environment.txt",
        "batch_script": manifest_path.parent / "cgl_lf_stage_i.sbatch",
        "slurm_log": root / "logs" / "slurm" / "%x.%j.log",
    }
    for key, expected in expected_paths.items():
        if declared_path(paths.get(key), f"launch path {key}") != expected:
            raise ValidationError(f"manifest launch path {key} is not exact production intent")
    validate_restart_launch_intent(
        command, manifest_path, restart_storage, nodes * ranks_per_node
    )


def regenerated_batch_script(
    manifest: dict[str, object], manifest_path: Path
) -> str:
    """Regenerate the exact production launch script from retained intent."""

    run = manifest["run"]
    allocation = manifest["allocation"]
    command = manifest["command"]
    paths = manifest["paths"]
    overrides = " ".join(quote(value) for value in command["overrides"])
    if overrides:
        overrides = " " + overrides
    restart = command.get("restart_file")
    restart_literal = quote(str(restart)) if restart else "''"
    bundle = command.get("source_bundle")
    runtime_checks = [
        f"require_sha {quote(command['production_utility']['sha256'])} "
        f"{quote(command['production_utility']['path'])} production_utility",
        f"require_sha {quote(command['input_sha256'])} \"${{INPUT}}\" input",
        f"require_sha {quote(command['matrix_sha256'])} {quote(command['matrix_file'])} matrix",
        f"require_sha {quote(command['executable_sha256'])} \"${{ATHENA}}\" executable",
    ]
    if isinstance(bundle, dict):
        runtime_checks.append(
            f"require_sha {quote(bundle['sha256'])} {quote(bundle['path'])} source_bundle"
        )
    qualification = command.get("qualification_approval")
    if isinstance(qualification, dict):
        runtime_checks.append(
            f"require_sha {quote(qualification['sha256'])} "
            f"{quote(qualification['path'])} qualification_approval"
        )
    restart_records = command.get("restart_files", [])
    if not isinstance(restart_records, list):
        raise ValidationError("prepared restart sibling metadata is invalid")
    for index, record in enumerate(restart_records):
        if not isinstance(record, dict):
            raise ValidationError("prepared restart sibling metadata is invalid")
        runtime_checks.append(
            f"require_sha {quote(record['sha256'])} {quote(record['path'])} "
            f"restart_{index:04d}"
        )
    archive_root = manifest_path.parent / "submitted_restart"
    if restart_records:
        require_directory_path(archive_root, "submitted restart archive")
        runtime_checks.append(
            f'test "$(find {quote(archive_root)} -type f -print | wc -l)" '
            f'-eq {len(restart_records)} || {{ echo "restart inventory changed" >&2; '
            "exit 1; }"
        )
    runtime_checks_text = "\n".join(runtime_checks)
    return f"""#!/bin/bash
#SBATCH -J {expected_job_name(manifest)}
#SBATCH -A {ACCOUNT}
#SBATCH -o {paths["slurm_log"]}
#SBATCH -p {PARTITION}
# Frontier default normal QOS is required for production; do not add -q debug.
#SBATCH -t {allocation["requested_walltime"]}
#SBATCH -N {allocation["nodes"]}
#SBATCH --gpus-per-node=8
#SBATCH --threads-per-core=1

set -euo pipefail

RUN_MANIFEST={quote(manifest_path)}
BATCH_SCRIPT_SHA256={BATCH_SCRIPT_DIGEST_PLACEHOLDER}
ATHENA={quote(command["executable"])}
INPUT={quote(command["input_file"])}
RESTART={restart_literal}
OUT_DIR={quote(paths["output_dir"])}
ENV_LOG={quote(paths["environment_log"])}
RANKS_PER_NODE={allocation["ranks_per_node"]}
CPUS_PER_TASK={allocation["cpus_per_task"]}
NNODES="${{SLURM_NNODES:?Missing SLURM_NNODES}}"
NRANKS="$((NNODES * RANKS_PER_NODE))"

require_sha() {{
  local expected="$1"
  local path="$2"
  local label="$3"
  local actual
  test -f "${{path}}" || {{ echo "missing ${{label}}: ${{path}}" >&2; exit 1; }}
  actual="$(sha256sum "${{path}}" | awk '{{print $1}}')"
  test "${{actual}}" = "${{expected}}" || {{
    echo "checksum mismatch for ${{label}}: ${{path}}" >&2
    exit 1
  }}
}}

normalized_script_sha256() {{
  sed -E 's/^BATCH_SCRIPT_SHA256=[0-9a-f]{{64}}$/BATCH_SCRIPT_SHA256={BATCH_SCRIPT_DIGEST_PLACEHOLDER}/' "$0" |
    sha256sum | awk '{{print $1}}'
}}

test -f "${{RUN_MANIFEST}}"
test -x "${{ATHENA}}"
test -f "${{INPUT}}"
if [[ -n "${{RESTART}}" ]]; then
  test -f "${{RESTART}}"
fi
test "$(normalized_script_sha256)" = "${{BATCH_SCRIPT_SHA256}}" || {{
  echo "batch script checksum mismatch" >&2
  exit 1
}}
{runtime_checks_text}
mkdir -p "${{OUT_DIR}}"

module restore
module load PrgEnv-cray
module load craype-accel-amd-gfx90a
module load cpe/25.09 cray-mpich/9.0.1 rocm/6.4.2
module load cce/20.0.0
module unload darshan-runtime

export LD_LIBRARY_PATH=${{CRAY_LD_LIBRARY_PATH}}:${{LD_LIBRARY_PATH:-}}
export MPICH_ENV_DISPLAY=1
export MPICH_VERSION_DISPLAY=1
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=0
export MPICH_OFI_NIC_POLICY=GPU
export MPICH_GPU_IPC_CACHE_MAX_SIZE=1000
export MPICH_MPIIO_HINTS="*:romio_cb_write=disable"
export MPICH_OFI_NUM_CQ_ENTRIES=131072
export FI_MR_CACHE_MONITOR=kdreg2
export FI_CXI_RX_MATCH_MODE=software
export HSA_XNACK=1
export OMP_NUM_THREADS=1

{{
  date -u +"started_utc=%Y-%m-%dT%H:%M:%SZ"
  echo "slurm_job_id=${{SLURM_JOB_ID:?Missing SLURM_JOB_ID}}"
  echo "prepared_manifest=${{RUN_MANIFEST}}"
  echo "nodes=${{NNODES}}"
  echo "ranks=${{NRANKS}}"
  module -t list 2>&1
  env | LC_ALL=C sort | grep -E '^(MPICH_|FI_|HSA_|OMP_|ROCR_|HIP_|CRAY_)' || true
}} > "${{ENV_LOG}}"

RUN_ARGS=(-i "${{INPUT}}")
if [[ -n "${{RESTART}}" ]]; then
  RUN_ARGS=(-r "${{RESTART}}")
fi

srun -N "${{NNODES}}" -n "${{NRANKS}}" --ntasks-per-node="${{RANKS_PER_NODE}}" \
  -c "${{CPUS_PER_TASK}}" --threads-per-core=1 --cpu-bind=threads \
  --gpus-per-task=1 --gpu-bind=closest \
  "${{ATHENA}}" "${{RUN_ARGS[@]}}" -d "${{OUT_DIR}}" \
  -t {quote(command["athena_walltime"])} \
  job/basename={quote(run["run_basename"])}{overrides}

date -u +"finished_utc=%Y-%m-%dT%H:%M:%SZ" >> "${{ENV_LOG}}"
"""


def matrix_case_record(
    matrix: dict[str, object], run: dict[str, object]
) -> dict[str, object]:
    """Require the retained matrix to contain the complete unique R02-R17 mapping."""

    cases = matrix.get("cases")
    if not isinstance(cases, list) or not all(isinstance(item, dict) for item in cases):
        raise ValidationError("retained matrix lacks a case list")
    by_id: dict[str, dict[str, object]] = {}
    for item in cases:
        case_id = item.get("id")
        if not isinstance(case_id, str) or case_id in by_id:
            raise ValidationError("retained matrix case IDs are invalid or duplicated")
        by_id[case_id] = item
    if list(by_id) != sorted(APPROVED_CASE_IDS):
        raise ValidationError("retained matrix is not the complete R02-R17 mapping")
    case_id = str(run["case_id"])
    case = by_id[case_id]
    for matrix_key, run_key in (
        ("name", "case_name"),
        ("resolution", "resolution"),
        ("figure_roles", "figure_roles"),
    ):
        if case.get(matrix_key) != run.get(run_key):
            raise ValidationError(f"retained matrix {matrix_key} differs from manifest run")
    frozen_parameters_for_case(case_id, case)
    return case


def authenticated_bundle_blobs(
    path: Path,
    requests: dict[str, tuple[str, str]],
    expected_sha256: str | None = None,
) -> dict[str, bytes]:
    """Read exact committed objects through a stable private bare clone."""

    for revision, repository_path in requests.values():
        require_revision(revision, "bundle object revision")
        candidate = Path(repository_path)
        if candidate.is_absolute() or ".." in candidate.parts or not candidate.parts:
            raise ValidationError("bundle object repository path is invalid")
    try:
        with tempfile.TemporaryDirectory(prefix="cgl-lf-validator-bundle-") as temporary:
            private_bundle = Path(temporary) / "source.bundle"
            digest = hashlib.sha256()
            with open_regular_binary(path, "source bundle") as (source, _profile):
                with private_bundle.open("wb") as target:
                    while True:
                        chunk = source.read(1024 * 1024)
                        if not chunk:
                            break
                        digest.update(chunk)
                        target.write(chunk)
            if (
                expected_sha256 is not None
                and digest.hexdigest()
                != require_sha256(expected_sha256, "source bundle digest")
            ):
                raise ValidationError("source bundle changed before object authentication")
            advertised = git_run(
                ["bundle", "list-heads", str(private_bundle)],
                timeout=120,
            )
            if advertised.returncode != 0:
                raise ValidationError(f"source bundle advertised refs are invalid: {path}")
            try:
                heads = {
                    require_revision(
                        line.split(maxsplit=1)[0], "source bundle advertised revision"
                    )
                    for line in advertised.stdout.decode("utf-8").splitlines()
                    if line.strip()
                }
            except (UnicodeDecodeError, IndexError) as error:
                raise ValidationError(
                    f"source bundle advertised refs are invalid: {path}"
                ) from error
            if not heads:
                raise ValidationError(f"source bundle advertises no retained revision: {path}")
            repository = Path(temporary) / "repository.git"
            clone = git_run(
                ["clone", "--bare", "--quiet", str(private_bundle), str(repository)],
                timeout=300,
            )
            if clone.returncode != 0:
                raise ValidationError(f"source bundle is not independently cloneable: {path}")
            payloads: dict[str, bytes] = {}
            for label, (revision, repository_path) in requests.items():
                commit = git_run(
                    ["-C", str(repository), "cat-file", "-e", f"{revision}^{{commit}}"],
                    timeout=120,
                )
                if commit.returncode != 0:
                    raise ValidationError(f"source bundle contents omit {label} revision")
                covered = any(
                    git_run([
                        "-C",
                        str(repository),
                        "merge-base",
                        "--is-ancestor",
                        revision,
                        head,
                    ]).returncode
                    == 0
                    for head in heads
                )
                if not covered:
                    raise ValidationError(
                        f"source bundle advertised history omits {label} revision"
                    )
                show = git_run(
                    ["-C", str(repository), "show", f"{revision}:{repository_path}"],
                    timeout=120,
                )
                if show.returncode != 0:
                    raise ValidationError(f"source bundle contents omit {label} object")
                payloads[label] = show.stdout
            return payloads
    except (OSError, subprocess.TimeoutExpired) as error:
        raise ValidationError(f"cannot independently inspect source bundle: {path}") from error


def require_exact_path(actual: object, expected: Path, label: str) -> Path:
    """Require one retained path to equal an exact production location."""

    path = declared_path(actual, label)
    if path != expected:
        raise ValidationError(f"{label} path differs from retained production layout")
    return path


def retained_product_members(
    record: object, label: str
) -> list[dict[str, object]]:
    """Return one exact shared or rank-local retained-product member list."""

    if not isinstance(record, dict):
        raise ValidationError(f"{label} record is invalid")
    rank_files = record.get("rank_files")
    members = rank_files if isinstance(rank_files, list) else [record]
    if not members or not all(isinstance(item, dict) for item in members):
        raise ValidationError(f"{label} member inventory is invalid")
    validated = []
    for item in members:
        path = declared_path(item.get("path"), f"{label} member")
        digest = require_sha256(item.get("sha256"), f"{label} member digest")
        size = require_exact_int(item.get("size_bytes"), f"{label} member size")
        if size <= 0:
            raise ValidationError(f"{label} member size is not positive")
        validated.append({"path": path, "sha256": digest, "size_bytes": size})
    if len({item["path"] for item in validated}) != len(validated):
        raise ValidationError(f"{label} member inventory is duplicated")
    storage = record.get("storage")
    if storage == "per_rank":
        expected = [f"rank_{rank:08d}" for rank in range(len(validated))]
        if [item["path"].parent.name for item in validated] != expected:
            raise ValidationError(f"{label} rank-member inventory is not contiguous")
    elif storage == "shared_mpiio":
        if len(validated) != 1:
            raise ValidationError(f"{label} shared inventory has multiple members")
    else:
        raise ValidationError(f"{label} storage class is invalid")
    return validated


def authenticate_parent_continuation(
    manifest: dict[str, object],
    manifest_path: Path,
    input_contract: dict[str, object],
    current_restart_records: list[dict[str, object]],
    profiles: dict[Path, tuple[int, ...]],
) -> dict[str, object] | None:
    """Bind a copied continuation restart exactly to its recorded parent state."""

    command = manifest["command"]
    parent = command.get("parent_segment")
    if not current_restart_records:
        if parent is not None:
            raise ValidationError("fresh segment contains parent-continuation metadata")
        return None
    required_parent_fields = {
        "execution_epoch",
        "manifest",
        "case_id",
        "segment",
        "result",
        "restart_sha256",
        "restart_files",
        "final_time",
        "restart_time",
        "input_sha256",
        "executable_sha256",
    }
    if not isinstance(parent, dict) or set(parent) != required_parent_fields:
        raise ValidationError("prepared continuation parent metadata is incomplete or ambiguous")
    require_segment_name(parent.get("segment"), "continuation parent segment name")
    root = declared_path(manifest.get("project_root"), "project root")
    parent_manifest_path = declared_path(parent.get("manifest"), "continuation parent manifest")
    require_within(
        parent_manifest_path,
        root / "runs" / "mks24-stage-i" / EXECUTION_EPOCH,
        "continuation parent manifest",
    )
    if parent_manifest_path == manifest_path:
        raise ValidationError("prepared continuation cannot name itself as parent")
    parent_manifest, parent_manifest_sha256 = load_json(
        parent_manifest_path, "continuation parent manifest", profiles
    )
    parent_run = parent_manifest.get("run")
    parent_command = parent_manifest.get("command")
    accounting = parent_manifest.get("accounting")
    parent_inspection = parent_manifest.get("scientific_inspection")
    if not all(
        isinstance(item, dict)
        for item in (parent_run, parent_command, accounting, parent_inspection)
    ):
        raise ValidationError("continuation parent lacks recorded execution metadata")
    require_job_id(parent_manifest.get("job_id"), "continuation parent job ID")
    require_segment_name(parent_run.get("segment"), "continuation parent run segment name")
    expected_parent_path = (
        root
        / "runs"
        / "mks24-stage-i"
        / EXECUTION_EPOCH
        / str(parent_run.get("case_id"))
        / str(parent_run.get("segment"))
        / "manifest"
        / "prepared_run.json"
    )
    if (
        parent_manifest_path != expected_parent_path
        or parent_manifest.get("schema_version") != MANIFEST_SCHEMA_VERSION
        or parent_manifest.get("execution_epoch") != EXECUTION_EPOCH
        or parent_manifest.get("state") != "recorded"
        or parent_manifest.get("project_root") != str(root)
        or parent_run.get("case_id") != manifest["run"].get("case_id")
        or parent.get("execution_epoch") != EXECUTION_EPOCH
        or parent.get("case_id") != parent_run.get("case_id")
        or parent.get("segment") != parent_run.get("segment")
        or parent.get("result") != accounting.get("result")
        or accounting.get("result") not in {"accepted", "clean_partial"}
        or accounting.get("execution_epoch") != EXECUTION_EPOCH
        or accounting.get("job_id") != parent_manifest.get("job_id")
        or accounting.get("case_id") != parent_run.get("case_id")
        or accounting.get("segment") != parent_run.get("segment")
    ):
        raise ValidationError("prepared continuation parent identity is inconsistent")
    parent_inspection_path = parent_manifest_path.parent / "segment_inspection.json"
    retained_inspection, parent_inspection_sha256 = load_json(
        parent_inspection_path, "continuation parent inspection", profiles
    )
    if retained_inspection != parent_inspection:
        raise ValidationError("continuation parent recorded inspection differs from retained file")
    if (
        parent_inspection.get("schema_version") != INSPECTION_SCHEMA_VERSION
        or parent_inspection.get("execution_epoch") != EXECUTION_EPOCH
        or parent_inspection.get("manifest") != str(parent_manifest_path)
        or parent_inspection.get("job_id") != parent_manifest.get("job_id")
        or parent_inspection.get("case_id") != parent_run.get("case_id")
        or parent_inspection.get("segment") != parent_run.get("segment")
    ):
        raise ValidationError("continuation parent inspection identity is inconsistent")
    try:
        final_time = float(parent["final_time"])
        restart_time = float(parent["restart_time"])
        inspection_final = float(parent_inspection["final_time"])
        terminal_time = float(parent_inspection["terminal_restart_time"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValidationError("continuation parent physical-time metadata is invalid") from error
    if (
        not all(math.isfinite(value) for value in (
            final_time, restart_time, inspection_final, terminal_time
        ))
        or final_time != restart_time
        or final_time != inspection_final
        or final_time != terminal_time
        or parent_inspection.get("clean_for_continuation") is not True
    ):
        raise ValidationError("continuation parent physical state is inconsistent")
    if (
        parent.get("input_sha256") != parent_command.get("input_sha256")
        or parent.get("input_sha256") != command.get("input_sha256")
        or parent.get("executable_sha256") != parent_command.get("executable_sha256")
        or parent.get("executable_sha256") != command.get("executable_sha256")
    ):
        raise ValidationError("continuation parent launch provenance differs from child")
    parent_members = retained_product_members(
        parent_inspection.get("terminal_restart"), "continuation parent terminal restart"
    )
    for item in parent_members:
        require_within(
            item["path"],
            parent_manifest_path.parents[1] / "output" / "rst",
            "continuation parent terminal restart",
        )
    parent_paths = parent.get("restart_files")
    if (
        not isinstance(parent_paths, list)
        or parent_paths != [str(item["path"]) for item in parent_members]
        or parent.get("restart_sha256") != parent_members[0]["sha256"]
        or command.get("source_restart_file") != str(parent_members[0]["path"])
        or len(current_restart_records) != len(parent_members)
    ):
        raise ValidationError("continuation restart inventory differs from recorded parent")
    for source, copied in zip(parent_members, current_restart_records):
        authenticated_file(
            source["path"],
            source["sha256"],
            "continuation parent terminal restart",
            profiles,
            source["size_bytes"],
        )
        if (
            copied["sha256"] != source["sha256"]
            or copied["size_bytes"] != source["size_bytes"]
        ):
            raise ValidationError("copied continuation restart bytes differ from recorded parent")
    abi = qualified_restart_abi(command)
    try:
        parent_target = float(parent_command.get("time_tlim_target"))
    except (TypeError, ValueError) as error:
        raise ValidationError("continuation parent prepared target is invalid") from error
    if not math.isfinite(parent_target):
        raise ValidationError("continuation parent prepared target is invalid")
    runtime_contract = {
        "run_basename": parent_run.get("run_basename"),
        "target_time": parent_target,
    }
    copied_paths = [Path(str(item["path"])) for item in current_restart_records]
    initial_state = grouped_restart_evidence(
        copied_paths, abi, input_contract, runtime_contract, profiles
    )
    if initial_state["time"] != final_time:
        raise ValidationError("copied continuation restart physical time differs from parent")
    return {
        "parent_manifest": str(parent_manifest_path),
        "parent_manifest_sha256": parent_manifest_sha256,
        "parent_inspection_sha256": parent_inspection_sha256,
        "parent_result": parent["result"],
        "start_time": final_time,
        "initial_restart_cycle": initial_state["cycle"],
        "initial_restart_dt": initial_state["dt"],
        "initial_history_last_time": initial_state["history_last_time"],
        "initial_restart_ordered_location_inventory": initial_state[
            "ordered_location_inventory"
        ],
        "initial_restart_ordered_location_cost_inventory": initial_state[
            "ordered_location_cost_inventory"
        ],
        "initial_restart_sha256": [
            str(item["sha256"]) for item in current_restart_records
        ],
        "parent_acceptance_validation": {
            "status": "trusted_authenticated_controller_record",
            "authorizing": False,
            "reason": (
                "the parent's accepted or clean-partial result and clean-for-continuation "
                "decision are trusted from authenticated retained controller records and "
                "are not independently replayed or authorized"
            ),
        },
    }


def authenticate_provenance(
    manifest: dict[str, object],
    manifest_path: Path,
    profiles: dict[Path, tuple[int, ...]],
) -> tuple[dict[str, object], dict[str, object], list[dict[str, object]]]:
    """Authenticate the prepared artifacts and mapped scientific input."""

    if manifest.get("schema_version") != MANIFEST_SCHEMA_VERSION:
        raise ValidationError("manifest is not a schema-3 prepared-run record")
    if manifest.get("execution_epoch") != EXECUTION_EPOCH:
        raise ValidationError("manifest execution epoch is not approved")
    root = declared_path(manifest.get("project_root"), "project root")
    if root != CANONICAL_PROJECT_ROOT:
        raise ValidationError("manifest project root is not the canonical Stage I root")
    require_directory_path(root, "project root")
    run = manifest.get("run")
    command = manifest.get("command")
    paths = manifest.get("paths")
    if not isinstance(run, dict) or not isinstance(command, dict) or not isinstance(paths, dict):
        raise ValidationError("manifest lacks run, command, or path provenance")
    case_id = run.get("case_id")
    segment = run.get("segment")
    if case_id not in APPROVED_CASE_IDS:
        raise ValidationError("manifest case is outside approved R02-R17")
    if not isinstance(segment, str) or SEGMENT_PATTERN.fullmatch(segment) is None:
        raise ValidationError("manifest segment name is invalid")
    expected_run = (
        root / "runs" / "mks24-stage-i" / EXECUTION_EPOCH / str(case_id) / segment
    )
    expected_manifest = expected_run / "manifest" / "prepared_run.json"
    if manifest_path != expected_manifest:
        raise ValidationError("manifest path differs from retained production layout")
    run_dir = require_exact_path(paths.get("run_dir"), expected_run, "run directory")
    output_dir = require_exact_path(
        paths.get("output_dir"), expected_run / "output", "output directory"
    )
    batch_path = require_exact_path(
        paths.get("batch_script"),
        expected_run / "manifest" / "cgl_lf_stage_i.sbatch",
        "batch script",
    )
    input_path = require_exact_path(
        command.get("input_file"),
        expected_run / "manifest" / "submitted_input.athinput",
        "archived input",
    )
    matrix_path = require_exact_path(
        command.get("matrix_file"),
        expected_run / "manifest" / "mks24_stage_i_manifest.json",
        "archived matrix",
    )
    require_directory_path(run_dir, "run directory")
    require_directory_path(output_dir, "output directory")
    artifacts: list[dict[str, object]] = []
    input_record = authenticated_file(
        input_path, command.get("input_sha256"), "archived input", profiles
    )
    matrix_record = authenticated_file(
        matrix_path, command.get("matrix_sha256"), "archived matrix", profiles
    )
    artifacts.extend((input_record, matrix_record))
    input_payload = read_regular_bytes(input_path, "archived input", profiles)
    matrix, _ = load_json(matrix_path, "archived matrix", profiles)
    case = matrix_case_record(matrix, run)
    input_blocks = parse_athinput(input_payload, "archived input")
    input_contract = require_input_contract(input_blocks, str(case_id), case)
    source_input = declared_path(command.get("source_input_file"), "source input")
    matrix_input = str(case["input"])
    matrix_input_parts = Path(matrix_input).parts
    if source_input.parts[-len(matrix_input_parts):] != matrix_input_parts:
        raise ValidationError("source input path differs from retained matrix mapping")
    if command.get("input_revision") != command.get("executable_revision"):
        raise ValidationError("input revision differs from executable revision")
    executable_revision = require_revision(
        command.get("executable_revision"), "executable revision"
    )
    executable = declared_path(command.get("executable"), "prepared executable")
    require_within(executable, root / "build", "prepared executable")
    executable_record = authenticated_file(
        executable, command.get("executable_sha256"), "prepared executable", profiles
    )
    if not os.access(executable, os.X_OK):
        raise ValidationError("prepared executable is not executable")
    artifacts.append(executable_record)
    bundle = command.get("source_bundle")
    utility = command.get("production_utility")
    if not isinstance(bundle, dict) or not isinstance(utility, dict):
        raise ValidationError("manifest lacks source-bundle or production-utility provenance")
    if utility.get("committed") is not True:
        raise ValidationError("production utility was not committed at preparation")
    utility_revision = require_revision(
        utility.get("revision"), "production utility revision"
    )
    utility_sha256 = require_sha256(utility.get("sha256"), "production utility digest")
    utility_path = declared_path(utility.get("path"), "production utility")
    if utility_path.parts[-3:] != Path(PRODUCTION_UTILITY_REPOSITORY_PATH).parts:
        raise ValidationError("production utility path differs from production convention")
    bundle_path = declared_path(bundle.get("path"), "source bundle")
    require_within(bundle_path, root / "source-archives", "source bundle")
    artifacts.append(
        authenticated_file(
            bundle_path, bundle.get("sha256"), "source bundle", profiles
        )
    )
    revisions = bundle.get("verified_revisions")
    required_revisions = {
        executable_revision,
        utility_revision,
        require_revision(command.get("input_revision"), "input revision"),
    }
    if (
        not isinstance(revisions, list)
        or not all(isinstance(item, str) for item in revisions)
        or len(revisions) != len(set(revisions))
        or set(revisions) != required_revisions
        or any(GIT_REVISION_PATTERN.fullmatch(item) is None for item in revisions)
    ):
        raise ValidationError("source bundle does not bind every launch revision")
    bundle_blobs = authenticated_bundle_blobs(
        bundle_path,
        {
            "input": (executable_revision, matrix_input),
            "matrix": (executable_revision, MATRIX_REPOSITORY_PATH),
            "production utility": (
                utility_revision,
                PRODUCTION_UTILITY_REPOSITORY_PATH,
            ),
        },
        str(bundle.get("sha256")),
    )
    if bundle_blobs["input"] != input_payload:
        raise ValidationError("archived input differs from authenticated bundle object")
    matrix_payload = read_regular_bytes(matrix_path, "archived matrix", profiles)
    if bundle_blobs["matrix"] != matrix_payload:
        raise ValidationError("archived matrix differs from authenticated bundle object")
    if hashlib.sha256(bundle_blobs["production utility"]).hexdigest() != utility_sha256:
        raise ValidationError(
            "production utility digest differs from authenticated bundle object"
        )
    qualification = command.get("qualification_approval")
    if not isinstance(qualification, dict):
        raise ValidationError("manifest lacks E03 qualification approval")
    qualification_path = require_exact_path(
        qualification.get("path"),
        root
        / "accounting"
        / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_qualification_approval.json",
        "qualification approval",
    )
    qualification_record = authenticated_file(
        qualification_path,
        qualification.get("sha256"),
        "qualification approval",
        profiles,
    )
    approval, _ = load_json(qualification_path, "qualification approval", profiles)
    require_qualification_approval(approval)
    artifacts.append(qualification_record)
    if qualification.get("token") != approval:
        raise ValidationError("qualification approval token differs from retained file")
    for record in (qualification, approval):
        if (
            record.get("execution_epoch") != EXECUTION_EPOCH
            or record.get("approved_executable_revision") != executable_revision
            or record.get("approved_executable_sha256") != executable_record["sha256"]
        ):
            raise ValidationError("qualification approval does not bind the executable")
    try:
        target = float(command["time_tlim_target"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValidationError("manifest lacks a finite time/tlim target") from error
    overrides = command.get("overrides")
    if (
        not math.isfinite(target)
        or not isinstance(overrides, list)
        or len(overrides) != 1
        or overrides[0] != f"time/tlim={format(target, '.17g')}"
    ):
        raise ValidationError("manifest does not bind exactly one time/tlim override")
    if command.get("allow_missing_restart_time_marker") is not False:
        raise ValidationError("prepared command permits a forbidden restart-marker bypass")
    validate_launch_intent_fields(
        manifest,
        manifest_path,
        root,
        case,
        str(input_contract["storage"]["restart"]),
    )
    batch_payload = read_regular_bytes(batch_path, "prepared batch script", profiles)
    expected_batch_digest = require_sha256(
        command.get("batch_script_sha256"), "batch script digest"
    )
    try:
        embedded = BATCH_SCRIPT_DIGEST_PATTERN.findall(batch_payload.decode("utf-8"))
    except UnicodeDecodeError as error:
        raise ValidationError("prepared batch script is not UTF-8 text") from error
    if embedded != [expected_batch_digest]:
        raise ValidationError("prepared batch script embedded digest is inconsistent")
    if normalized_batch_script_sha256(batch_payload) != expected_batch_digest:
        raise ValidationError("prepared batch script normalized digest differs")
    try:
        expected_batch = regenerated_batch_script(manifest, manifest_path).encode("utf-8")
    except (KeyError, TypeError, ValueError) as error:
        raise ValidationError("manifest cannot regenerate production launch intent") from error
    if BATCH_SCRIPT_DIGEST_PATTERN.sub(
        f"BATCH_SCRIPT_SHA256={BATCH_SCRIPT_DIGEST_PLACEHOLDER}",
        batch_payload.decode("utf-8"),
    ) != BATCH_SCRIPT_DIGEST_PATTERN.sub(
        f"BATCH_SCRIPT_SHA256={BATCH_SCRIPT_DIGEST_PLACEHOLDER}",
        expected_batch.decode("utf-8"),
    ):
        raise ValidationError("prepared batch script differs from regenerated launch intent")
    batch_digest = hashlib.sha256(batch_payload).hexdigest()
    artifacts.append({
        "path": str(batch_path),
        "size_bytes": len(batch_payload),
        "sha256": batch_digest,
    })
    restart_records = command.get("restart_files")
    if not isinstance(restart_records, list):
        raise ValidationError("manifest lacks prepared restart sibling provenance")
    retained_restart_paths: list[Path] = []
    retained_restart_records: list[dict[str, object]] = []
    for record in restart_records:
        if not isinstance(record, dict):
            raise ValidationError("prepared restart sibling provenance is invalid")
        path = declared_path(record.get("path"), "prepared restart sibling")
        require_within(path, manifest_path.parent, "prepared restart sibling")
        retained = authenticated_file(
            path,
            record.get("sha256"),
            "prepared restart sibling",
            profiles,
            record.get("size_bytes"),
        )
        artifacts.append(retained)
        retained_restart_records.append(retained)
        retained_restart_paths.append(path)
    if retained_restart_paths:
        restart_file = declared_path(command.get("restart_file"), "prepared restart")
        if restart_file not in retained_restart_paths:
            raise ValidationError("prepared restart representative is not retained")
        representative = next(
            record for record in artifacts if record["path"] == str(restart_file)
        )
        if command.get("restart_sha256") != representative["sha256"]:
            raise ValidationError("prepared restart representative digest differs")
        if not isinstance(command.get("parent_segment"), dict):
            raise ValidationError("prepared continuation lacks parent-segment provenance")
    elif any(
        command.get(key) is not None
        for key in ("restart_file", "restart_sha256", "source_restart_file", "parent_segment")
    ):
        raise ValidationError("fresh segment has inconsistent restart provenance")
    continuation = authenticate_parent_continuation(
        manifest,
        manifest_path,
        input_contract,
        retained_restart_records,
        profiles,
    )
    return (
        {
            "root": root,
            "run": run,
            "command": command,
            "output_dir": output_dir,
            "target": target,
            "case": case,
            "bundle_object_sha256": {
                label: hashlib.sha256(payload).hexdigest()
                for label, payload in bundle_blobs.items()
            },
            "continuation": continuation,
        },
        input_contract,
        artifacts,
    )


def expected_ranks(manifest: dict[str, object]) -> int:
    """Derive the exact rank count from retained allocation metadata."""

    allocation = manifest.get("allocation")
    if not isinstance(allocation, dict):
        raise ValidationError("manifest lacks allocation metadata")
    try:
        nodes = int(allocation["nodes"])
        ranks_per_node = int(allocation["ranks_per_node"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValidationError("manifest allocation does not define integer ranks") from error
    if nodes <= 0 or ranks_per_node <= 0:
        raise ValidationError("manifest allocation rank dimensions must be positive")
    return nodes * ranks_per_node


def product_groups(
    directory: Path, suffix: str, rank_count: int, label: str
) -> list[list[Path]]:
    """Rediscover shared and complete rank-local output products."""

    require_directory_path(directory, f"{label} directory")
    direct_files = sorted(
        path
        for path in directory.iterdir()
        if path.name.endswith(suffix) and path.is_file()
    )
    rank_directories = sorted(
        path
        for path in directory.iterdir()
        if path.is_dir() and RANK_DIRECTORY_PATTERN.fullmatch(path.name)
    )
    if rank_directories:
        expected_names = [f"rank_{rank:08d}" for rank in range(rank_count)]
        if [path.name for path in rank_directories] != expected_names:
            raise ValidationError(f"{label} rank directories are not complete and contiguous")
    groups = [[path] for path in direct_files]
    if rank_directories:
        names = sorted(
            path.name
            for path in rank_directories[0].iterdir()
            if path.name.endswith(suffix) and path.is_file()
        )
        for name in names:
            members = [directory / rank_name / name for rank_name in expected_names]
            for member in members:
                try:
                    profile = member.lstat()
                except OSError as error:
                    raise ValidationError(
                        f"{label} rank member is unavailable: {member}"
                    ) from error
                if not stat.S_ISREG(profile.st_mode):
                    raise ValidationError(f"{label} rank member is not regular: {member}")
            groups.append(members)
    grouped = sorted(path for group in groups for path in group)
    discovered = sorted(
        path
        for path in directory.rglob(f"*{suffix}")
        if path.is_file()
    )
    if grouped != discovered:
        raise ValidationError(f"{label} inventory contains an incomplete product")
    return sorted(groups, key=lambda group: str(group[0]))


def product_storage(group: list[Path]) -> str:
    """Return the schema-4 storage class of one discovered product group."""

    return "per_rank" if group[0].parent.name.startswith("rank_") else "shared_mpiio"


def captured_product(
    group: list[Path],
    output_root: Path,
    rank_count: int,
    label: str,
    profiles: dict[Path, tuple[int, ...]],
) -> dict[str, object]:
    """Capture one schema-4 shared or rank-local retained product."""

    if not group:
        raise ValidationError(f"{label} product group is empty")
    records = [captured_file(path, output_root, label, profiles) for path in group]
    representative = records[0].copy()
    storage = product_storage(group)
    if storage == "per_rank":
        expected_names = [f"rank_{rank:08d}" for rank in range(rank_count)]
        if (
            len(group) != rank_count
            or [path.parent.name for path in group] != expected_names
        ):
            raise ValidationError(f"{label} rank members are not complete and contiguous")
        representative["rank_files"] = records
    elif len(group) != 1:
        raise ValidationError(f"{label} shared product has multiple representatives")
    representative["storage"] = storage
    return representative


def flatten_product_records(records: list[dict[str, object]]) -> list[dict[str, object]]:
    """Return every unique retained file represented by product records."""

    flattened: dict[str, dict[str, object]] = {}
    for record in records:
        rank_files = record.get("rank_files")
        members = rank_files if isinstance(rank_files, list) else [record]
        for member in members:
            if not isinstance(member, dict) or not isinstance(member.get("path"), str):
                raise ValidationError("captured retained product record is invalid")
            flattened[str(member["path"])] = member
    return [flattened[path] for path in sorted(flattened)]


def parse_history(
    path: Path,
    label: str,
    profiles: dict[Path, tuple[int, ...]],
) -> dict[str, list[float]]:
    """Read a finite Athena history with unique contiguous indexed labels."""

    try:
        text = read_regular_bytes(path, label, profiles).decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValidationError(f"{label} is not UTF-8 text: {path}") from error
    labels: list[str] | None = None
    rows: list[list[float]] = []
    for original in text.splitlines():
        if original.startswith("#"):
            found = HISTORY_LABEL_PATTERN.findall(original)
            if found:
                indexed = sorted((int(index), name) for index, name in found)
                indices = [index for index, _ in indexed]
                names = [name for _, name in indexed]
                if (
                    not indices
                    or indices[0] not in (0, 1)
                    or indices != list(range(indices[0], indices[0] + len(indices)))
                    or len(set(names)) != len(names)
                ):
                    raise ValidationError(f"{label} labels are not unique and contiguous")
                if labels is not None and labels != names:
                    raise ValidationError(f"{label} contains conflicting label headers")
                labels = names
            continue
        if not original.strip():
            continue
        try:
            row = [float(value) for value in original.split()]
        except ValueError as error:
            raise ValidationError(f"{label} contains a non-numeric row") from error
        if not all(math.isfinite(value) for value in row):
            raise ValidationError(f"{label} contains a non-finite row")
        rows.append(row)
    if labels is None or len(rows) < 2:
        raise ValidationError(f"{label} lacks labels or sufficient data")
    if any(len(row) != len(labels) for row in rows):
        raise ValidationError(f"{label} row width differs from its labels")
    return {name: [row[index] for row in rows] for index, name in enumerate(labels)}


def require_columns(history: dict[str, list[float]], names: set[str], label: str) -> None:
    """Require a set of history columns."""

    missing = sorted(names - set(history))
    if missing:
        raise ValidationError(f"{label} lacks required columns: {missing}")


def require_exact_historical_history_columns(
    history: dict[str, list[float]],
    expected: object,
    label: str,
) -> None:
    """Bind one history to the frozen qualified executable's exact schema."""

    if (
        not isinstance(expected, tuple)
        or not expected
        or not all(isinstance(name, str) and name for name in expected)
        or len(set(expected)) != len(expected)
    ):
        raise ValidationError(f"{label} qualified historical schema is invalid")
    if tuple(history) != expected:
        raise ValidationError(
            f"{label} columns differ from qualified historical executable schema"
        )


def require_strictly_increasing(values: list[float], label: str) -> None:
    """Require an exactly ordered physical-time sequence."""

    if any(right <= left for left, right in zip(values, values[1:])):
        raise ValidationError(f"{label} is not strictly increasing")


def require_time_coverage(
    values: list[float],
    start: float,
    final: float,
    cadence: float,
    label: str,
    *,
    require_start: bool = True,
) -> None:
    """Require exact segment endpoints and complete expected-cadence coverage."""

    if not values:
        raise ValidationError(f"{label} is empty")
    if not math.isfinite(cadence) or cadence <= 0.0:
        raise ValidationError(f"{label} cadence is invalid")
    require_strictly_increasing(values, label)
    if require_start and values[0] != start:
        raise ValidationError(f"{label} does not begin at the prepared segment start")
    if not require_start and values[0] < start:
        raise ValidationError(
            f"{label} continuation products precede the authenticated start"
        )
    if values[-1] != final:
        raise ValidationError(f"{label} does not end at the exact segment endpoint")
    maximum_gap = cadence * 1.01 + 1.0e-12
    if (
        (not require_start and values[0] - start > maximum_gap)
        or any(right - left > maximum_gap for left, right in zip(values, values[1:]))
    ):
        raise ValidationError(f"{label} does not provide expected cadence coverage")


def require_history_time_coverage(
    values: list[float],
    dts: list[float],
    start: float,
    final: float,
    target: float,
    segment_result: str,
    cadence: float,
    label: str,
    continuation_history_last_time: float | None,
) -> dict[str, object]:
    """Replay AthenaK's history schedule, failing closed on terminal ambiguity."""

    def as_float32(value: float) -> float:
        try:
            return struct.unpack("<f", struct.pack("<f", value))[0]
        except (OverflowError, struct.error) as error:
            raise ValidationError(f"{label} schedule is not representable") from error

    def latest_forward_addition_predecessor(
        current: float, dt: float, predecessor_floor: float
    ) -> float:
        """Find the latest binary64 state that can advance exactly to current."""

        upper = current
        if predecessor_floor < 0.0 or upper < predecessor_floor:
            raise ValidationError(f"{label} does not provide expected cadence coverage")

        def bits(value: float) -> int:
            return struct.unpack("<Q", struct.pack("<d", value))[0]

        def value(encoded: int) -> float:
            return struct.unpack("<d", struct.pack("<Q", encoded))[0]

        lower_encoded = bits(predecessor_floor)
        upper_exclusive = bits(upper) + 1

        def first_true(predicate: Callable[[float], bool]) -> int:
            lower = lower_encoded
            upper_bound = upper_exclusive
            while lower < upper_bound:
                middle = lower + (upper_bound - lower) // 2
                if predicate(value(middle)):
                    upper_bound = middle
                else:
                    lower = middle + 1
            return lower

        first_equal = first_true(lambda predecessor: predecessor + dt >= current)
        first_greater = first_true(lambda predecessor: predecessor + dt > current)
        if (
            first_equal == upper_exclusive
            or first_equal == first_greater
            or value(first_equal) + dt != current
        ):
            raise ValidationError(f"{label} does not provide expected cadence coverage")
        return value(first_greater - 1)

    if not values:
        raise ValidationError(f"{label} is empty")
    if len(values) != len(dts) or any(
        not math.isfinite(dt) or dt <= 0.0 for dt in dts
    ):
        raise ValidationError(f"{label} timestep evidence is invalid")
    if not math.isfinite(cadence) or cadence <= 0.0:
        raise ValidationError(f"{label} cadence is invalid")
    if segment_result not in {"accepted", "clean_partial"}:
        raise ValidationError(f"{label} segment result is invalid")
    if values[-1] != final:
        raise ValidationError(f"{label} does not end at the exact segment endpoint")
    if any(value < start for value in values):
        raise ValidationError(
            f"{label} continuation products precede the authenticated start"
        )
    duplicate_terminal = len(values) >= 2 and values[-2] == final
    if duplicate_terminal and segment_result != "clean_partial":
        raise ValidationError(f"{label} contains an impossible duplicate terminal output")
    normal_rows_end = len(values) - 1
    if any(
        right <= left
        for left, right in zip(values[:normal_rows_end], values[1:normal_rows_end])
    ) or (
        not duplicate_terminal
        and normal_rows_end > 0
        and values[normal_rows_end - 1] >= final
    ):
        raise ValidationError(f"{label} is not a valid retained output timeline")

    first_scheduled_index = 0
    if continuation_history_last_time is None:
        if values[0] != start:
            raise ValidationError(f"{label} does not begin at the prepared segment start")
        schedule_last_time = start
        first_scheduled_index = 1
        origin = "fresh_initial_output"
    else:
        if (
            not math.isfinite(continuation_history_last_time)
            or continuation_history_last_time < 0.0
            or as_float32(continuation_history_last_time - cadence)
            > as_float32(start)
        ):
            raise ValidationError(
                f"{label} authenticated parent schedule is inconsistent"
            )
        schedule_last_time = continuation_history_last_time
        if values[0] == start:
            raise ValidationError(
                f"{label} continuation contains impossible exact-start output"
            )
        origin = "authenticated_parent_restart_schedule"

    scheduled_outputs = 0
    backlog = (
        continuation_history_last_time is not None
        and as_float32(start) >= as_float32(schedule_last_time + cadence)
    )
    target_32 = as_float32(target)
    for index in range(first_scheduled_index, normal_rows_end):
        value = values[index]
        expected = schedule_last_time + cadence
        predecessor_floor = values[index - 1] if index > 0 else start
        earliest_current = predecessor_floor + dts[index]
        value_32 = as_float32(value)
        expected_32 = as_float32(expected)
        latest_previous = latest_forward_addition_predecessor(
            value, dts[index], predecessor_floor
        )
        if (
            value < earliest_current
            or (
                backlog
                and (
                    value != earliest_current
                    or latest_previous != predecessor_floor
                )
            )
            or (
                not backlog
                and as_float32(latest_previous) >= expected_32
            )
            or value_32 < expected_32
            or value_32 >= target_32
        ):
            raise ValidationError(f"{label} does not provide expected cadence coverage")
        schedule_last_time = expected
        scheduled_outputs += 1
        backlog = value_32 >= as_float32(schedule_last_time + cadence)

    next_scheduled_time = schedule_last_time + cadence
    final_32 = as_float32(final)
    terminal_normal_due = duplicate_terminal or (
        final_32 >= as_float32(next_scheduled_time) and final_32 < target_32
    )
    if segment_result == "clean_partial":
        if not duplicate_terminal and terminal_normal_due:
            raise ValidationError(f"{label} does not provide expected cadence coverage")
        terminal_form = (
            "scheduled_normal_then_finalize_duplicate"
            if duplicate_terminal
            else "unscheduled_finalize_only"
        )
    else:
        if as_float32(next_scheduled_time) < final_32:
            raise ValidationError(f"{label} does not provide expected cadence coverage")
        terminal_form = "exact_target_finalize_only"

    return {
        "status": "athenak_timestep_triggered_schedule_complete",
        "origin": origin,
        "cadence": cadence,
        "authenticated_parent_history_last_time": continuation_history_last_time,
        "retained_rows": len(values),
        "scheduled_nonterminal_outputs": scheduled_outputs,
        "next_scheduled_time": next_scheduled_time,
        "terminal_form": terminal_form,
        "terminal_normal_output_due": terminal_normal_due,
        "terminal_completeness_policy": (
            "fail_closed_no_unconsumed_nominal_schedule_before_terminal"
            if segment_result == "accepted"
            else "exact_clean_partial_normal_output_plus_finalize_form"
        ),
        "exact_terminal_output": final,
    }


def require_cycle_timeline(values: list[int], label: str) -> None:
    """Require exact nonnegative, strictly increasing output-cycle metadata."""

    if not values or any(type(value) is not int or value < 0 for value in values):
        raise ValidationError(f"{label} contains an invalid cycle")
    if any(right <= left for left, right in zip(values, values[1:])):
        raise ValidationError(f"{label} is not strictly increasing")


def require_matching_output_cycles(
    snapshots: list[dict[str, object]],
    restarts: list[dict[str, object]],
) -> None:
    """Require snapshot/restart cycles to agree whenever physical times coincide."""

    snapshot_cycles = {float(item["time"]): int(item["cycle"]) for item in snapshots}
    for item in restarts:
        time = float(item["time"])
        if time in snapshot_cycles and snapshot_cycles[time] != int(item["cycle"]):
            raise ValidationError(
                "snapshot and restart cycles disagree at a shared physical time"
            )


def require_matching_shared_output_orders(
    snapshots: list[dict[str, object]],
    restarts: list[dict[str, object]],
) -> None:
    """Bind coincident shared snapshots and restarts to one exact logical order."""

    snapshot_orders = {
        float(item["time"]): item["ordered_location_inventory"][0]
        for item in snapshots
    }
    for item in restarts:
        time = float(item["time"])
        if (
            time in snapshot_orders
            and snapshot_orders[time] != item["ordered_location_inventory"]
        ):
            raise ValidationError(
                "shared snapshot and restart ordered logical-location inventories disagree"
            )


def require_continuation_initial_state(
    continuation: dict[str, object] | None,
    snapshots: list[dict[str, object]],
    restarts: list[dict[str, object]],
) -> None:
    """Reject child products before start and bind any exact-start products."""

    if continuation is None:
        return
    start_time = float(continuation["start_time"])
    initial_cycle = int(continuation["initial_restart_cycle"])
    initial_dt = float(continuation["initial_restart_dt"])
    initial_order = continuation.get("initial_restart_ordered_location_inventory")
    initial_location_cost_pairs = continuation.get(
        "initial_restart_ordered_location_cost_inventory"
    )
    if any(
        float(item["time"]) < start_time for item in [*snapshots, *restarts]
    ):
        raise ValidationError(
            "continuation retains output before authenticated start time"
        )
    start_snapshots = [
        item for item in snapshots if float(item["time"]) == start_time
    ]
    start_restarts = [
        item for item in restarts if float(item["time"]) == start_time
    ]
    if any(int(item["cycle"]) != initial_cycle for item in start_snapshots + start_restarts):
        raise ValidationError(
            "continuation exact-start output cycles differ from authenticated initial state"
        )
    if any(float(item["dt"]) != initial_dt for item in start_restarts):
        raise ValidationError(
            "continuation exact-start restart dt differs from authenticated initial state"
        )
    if initial_order is not None and any(
        item.get("ordered_location_inventory") != initial_order for item in restarts
    ):
        raise ValidationError(
            "continuation restart logical-location order differs from authenticated parent"
        )
    if initial_location_cost_pairs is not None and any(
        item.get("ordered_location_cost_inventory") != initial_location_cost_pairs
        for item in restarts
    ):
        raise ValidationError(
            "continuation restart ordered logical-location/cost pairs differ from "
            "authenticated parent"
        )


def read_exact(stream: BinaryIO, count: int, label: str) -> bytes:
    """Read an exact positive byte count from one product."""

    payload = stream.read(count)
    if len(payload) != count:
        raise ValidationError(f"{label} is truncated")
    return payload


def require_finite_binary_values(
    stream: BinaryIO, count: int, value_size: int, label: str
) -> None:
    """Stream-decode one complete little-endian floating payload."""

    if count < 0 or value_size not in (4, 8) or count % value_size != 0:
        raise ValidationError(f"{label} has an invalid floating payload size")
    format_code = "f" if value_size == 4 else "d"
    remaining = count
    chunk_limit = (1024 * 1024 // value_size) * value_size
    while remaining:
        payload = read_exact(stream, min(remaining, chunk_limit), label)
        if any(
            not math.isfinite(value)
            for (value,) in struct.iter_unpack(f"<{format_code}", payload)
        ):
            raise ValidationError(f"{label} contains a non-finite value")
        remaining -= len(payload)


def require_product_parameter_contract(
    blocks: dict[str, dict[str, str]],
    input_contract: dict[str, object],
    runtime_contract: dict[str, object],
    label: str,
    *,
    restart_product: bool,
) -> None:
    """Require the exact qualified runtime parameter inventory and values."""

    archived = flattened_parameters(dict(input_contract["archived_parameter_blocks"]))
    expected = dict(archived)
    expected.update(QUALIFIED_PRODUCT_PARAMETER_ADDITIONS)
    optional = (
        frozenset() if restart_product else QUALIFIED_SNAPSHOT_OPTIONAL_PARAMETERS
    )
    actual = flattened_parameters(blocks)
    missing = set(expected) - set(actual)
    extra = set(actual) - set(expected)
    if extra or missing - optional:
        raise ValidationError(
            f"{label} parameter inventory differs from qualified contract; "
            f"missing={sorted(missing)}, extra={sorted(extra)}"
        )
    expected[("job", "basename")] = str(runtime_contract["run_basename"])
    expected[("time", "tlim")] = format(float(runtime_contract["target_time"]), ".17g")
    for key, value in expected.items():
        if key not in actual:
            continue
        if value is None:
            continue
        if not equivalent_parameter_value(actual[key], value):
            raise ValidationError(
                f"{label} violates qualified {key[0]}/{key[1]}={value}"
            )
    for key in PRODUCT_NONNEGATIVE_INTEGER_PARAMETERS:
        if key not in actual:
            continue
        try:
            value = int(actual[key])
        except ValueError as error:
            raise ValidationError(
                f"{label} qualified {key[0]}/{key[1]} is not an integer"
            ) from error
        if value < 0 or str(value) != actual[key]:
            raise ValidationError(
                f"{label} qualified {key[0]}/{key[1]} is not a canonical nonnegative integer"
            )
    for key in PRODUCT_NONNEGATIVE_FINITE_PARAMETERS:
        if key not in actual:
            continue
        try:
            value = float(actual[key])
        except ValueError as error:
            raise ValidationError(
                f"{label} qualified {key[0]}/{key[1]} is not numeric"
            ) from error
        if not math.isfinite(value) or value < 0.0:
            raise ValidationError(
                f"{label} qualified {key[0]}/{key[1]} is not nonnegative and finite"
            )
    for key in PRODUCT_SCHEDULE_TIME_PARAMETERS:
        try:
            value = float(actual[key])
        except ValueError as error:
            raise ValidationError(
                f"{label} qualified {key[0]}/{key[1]} is not numeric"
            ) from error
        if not math.isfinite(value) or (value != -1.0 and value < 0.0):
            raise ValidationError(
                f"{label} qualified {key[0]}/{key[1]} is not a scheduling time or sentinel"
            )


def expected_logical_locations(
    input_contract: dict[str, object],
) -> set[tuple[int, int, int, int]]:
    """Return every valid level-zero logical meshblock coordinate."""

    divisions = tuple(int(value) for value in input_contract["divisions"])
    return {
        (i, j, k, 0)
        for k in range(divisions[2])
        for j in range(divisions[1])
        for i in range(divisions[0])
    }


def expected_region_indices(
    input_contract: dict[str, object], shape: tuple[int, int, int]
) -> tuple[int, ...]:
    """Return every initialized RegionIndcs field for one Stage I region."""

    nghost = int(input_contract["nghost"])
    mesh = tuple(int(value) for value in input_contract["mesh"])
    starts = tuple(nghost if axis == 0 or mesh[axis] > 1 else 0 for axis in range(3))
    ends = tuple(start + size - 1 for start, size in zip(starts, shape))
    coarse = tuple(max(1, size // 2) for size in shape)
    coarse_starts = starts
    coarse_ends = tuple(
        start + size - 1 for start, size in zip(coarse_starts, coarse)
    )
    return (
        nghost,
        *shape,
        starts[0],
        ends[0],
        starts[1],
        ends[1],
        starts[2],
        ends[2],
        *coarse,
        coarse_starts[0],
        coarse_ends[0],
        coarse_starts[1],
        coarse_ends[1],
        coarse_starts[2],
        coarse_ends[2],
    )


def expected_mesh_geometry(input_contract: dict[str, object]) -> tuple[float, ...]:
    """Return the exact root RegionSize values implied by the archived input."""

    domain = tuple(tuple(pair) for pair in input_contract["domain"])
    mesh = tuple(int(value) for value in input_contract["mesh"])
    return (
        *(pair[0] for pair in domain),
        *(pair[1] for pair in domain),
        *((upper - lower) / cells for (lower, upper), cells in zip(domain, mesh)),
    )


def expected_block_geometry(
    input_contract: dict[str, object], location: tuple[int, int, int, int]
) -> tuple[float, ...]:
    """Return exact binary-output bounds for one level-zero MeshBlock."""

    domain = tuple(tuple(pair) for pair in input_contract["domain"])
    divisions = tuple(int(value) for value in input_contract["divisions"])
    if location[3] != 0:
        raise ValidationError("Stage I binary snapshot contains a refined MeshBlock")
    bounds = []
    for logical, count, (lower, upper) in zip(location[:3], divisions, domain):
        width = (upper - lower) / count
        bounds.extend((lower + logical * width, lower + (logical + 1) * width))
    return tuple(bounds)


def require_close_sequence(
    actual: tuple[float, ...], expected: tuple[float, ...], label: str
) -> None:
    """Require one finite geometry sequence to match its archived-input contract."""

    if len(actual) != len(expected) or any(
        not math.isfinite(value)
        or not math.isclose(value, target, rel_tol=1.0e-13, abs_tol=1.0e-14)
        for value, target in zip(actual, expected)
    ):
        raise ValidationError(f"{label} differs from archived input")


def read_bounded_line(stream: BinaryIO, label: str) -> bytes:
    """Read one bounded newline-terminated product header line."""

    line = stream.readline(MAX_TEXT_LINE_BYTES + 1)
    if not line.endswith(b"\n") or len(line) > MAX_TEXT_LINE_BYTES:
        raise ValidationError(f"{label} contains an invalid header line")
    return line


def parse_header_assignment(line: bytes, key: str, label: str) -> str:
    """Parse one exact key=value binary-output header line."""

    try:
        text = line.decode("utf-8").strip()
    except UnicodeDecodeError as error:
        raise ValidationError(f"{label} header is not UTF-8") from error
    prefix = f"{key}="
    if not text.startswith(prefix):
        raise ValidationError(f"{label} lacks header field {key}")
    return text[len(prefix):].strip()


def binary_snapshot_evidence(
    path: Path,
    input_contract: dict[str, object],
    runtime_contract: dict[str, object],
    profiles: dict[Path, tuple[int, ...]],
) -> dict[str, object]:
    """Structurally parse one complete Athena binary snapshot."""

    with open_regular_binary(path, "binary snapshot") as (stream, profile):
        remember_profile(profiles, path, profile, "binary snapshot")
        if read_bounded_line(stream, "binary snapshot") != b"Athena binary output version=1.1\n":
            raise ValidationError(f"binary snapshot has an unqualified format: {path}")
        try:
            preheader = int(
                parse_header_assignment(
                    read_bounded_line(stream, "binary snapshot"),
                    "size of preheader",
                    "binary snapshot",
                )
            )
            time = float(
                parse_header_assignment(
                    read_bounded_line(stream, "binary snapshot"), "time", "binary snapshot"
                )
            )
            cycle = int(
                parse_header_assignment(
                    read_bounded_line(stream, "binary snapshot"), "cycle", "binary snapshot"
                )
            )
            location_size = int(
                parse_header_assignment(
                    read_bounded_line(stream, "binary snapshot"),
                    "size of location",
                    "binary snapshot",
                )
            )
            variable_size = int(
                parse_header_assignment(
                    read_bounded_line(stream, "binary snapshot"),
                    "size of variable",
                    "binary snapshot",
                )
            )
            variable_count = int(
                parse_header_assignment(
                    read_bounded_line(stream, "binary snapshot"),
                    "number of variables",
                    "binary snapshot",
                )
            )
        except ValueError as error:
            raise ValidationError(f"binary snapshot header is not numeric: {path}") from error
        if (
            preheader != 5
            or not math.isfinite(time)
            or cycle < 0
            or location_size not in (4, 8)
            or variable_size not in (4, 8)
            or variable_count <= 0
        ):
            raise ValidationError(f"binary snapshot header is not loadable: {path}")
        variables_line = read_bounded_line(stream, "binary snapshot")
        try:
            variables_text = variables_line.decode("utf-8").strip()
        except UnicodeDecodeError as error:
            raise ValidationError(f"binary snapshot variables are not UTF-8: {path}") from error
        if not variables_text.startswith("variables:"):
            raise ValidationError(f"binary snapshot lacks a variable inventory: {path}")
        variables = tuple(variables_text.split(":", 1)[1].split())
        if len(variables) != variable_count or len(set(variables)) != variable_count:
            raise ValidationError(f"binary snapshot variable inventory is invalid: {path}")
        try:
            header_size = int(
                parse_header_assignment(
                    read_bounded_line(stream, "binary snapshot"),
                    "header offset",
                    "binary snapshot",
                )
            )
        except ValueError as error:
            raise ValidationError(f"binary snapshot header offset is invalid: {path}") from error
        if header_size <= 0:
            raise ValidationError(f"binary snapshot parameter header is empty: {path}")
        parameter_dump = read_exact(stream, header_size, "binary snapshot parameter header")
        if b"<par_end>\n" not in parameter_dump:
            raise ValidationError(f"binary snapshot parameter header is not loadable: {path}")
        require_product_parameter_contract(
            parse_athinput(parameter_dump, "binary snapshot parameter header"),
            input_contract,
            runtime_contract,
            "binary snapshot parameter header",
            restart_product=False,
        )
        if variables != EXPECTED_SNAPSHOT_VARIABLES:
            raise ValidationError(
                f"binary snapshot variables differ from mhd_w_bcc contract: {path}"
        )
        locations: list[tuple[int, int, int, int]] = []
        seen_locations: set[tuple[int, int, int, int]] = set()
        valid_locations = expected_logical_locations(input_contract)
        block_shape = tuple(int(value) for value in input_contract["meshblock"])
        expected_indices = expected_region_indices(input_contract, block_shape)[:10]
        geometry_format = "<6f" if location_size == 4 else "<6d"
        geometry_size = struct.calcsize(geometry_format)
        while stream.tell() < profile[6]:
            indices = struct.unpack("<10i", read_exact(stream, 40, "binary snapshot block"))
            if indices[:6] != expected_indices[4:10]:
                raise ValidationError(
                    f"binary snapshot active-region indices differ from input: {path}"
                )
            dimensions = block_shape
            location = tuple(indices[6:10])
            if location in seen_locations or location not in valid_locations:
                raise ValidationError(f"binary snapshot logical inventory is invalid: {path}")
            locations.append(location)
            seen_locations.add(location)
            geometry = struct.unpack(
                geometry_format,
                read_exact(stream, geometry_size, "binary snapshot geometry"),
            )
            require_close_sequence(
                geometry,
                expected_block_geometry(input_contract, location),
                f"binary snapshot geometry: {path}",
            )
            payload_size = math.prod(dimensions) * variable_count * variable_size
            if stream.tell() + payload_size > profile[6]:
                raise ValidationError(f"binary snapshot field payload is truncated: {path}")
            require_finite_binary_values(
                stream, payload_size, variable_size, "binary snapshot field payload"
            )
        if stream.tell() != profile[6] or not locations:
            raise ValidationError(f"binary snapshot lacks a complete block payload: {path}")
    return {
        "time": time,
        "cycle": cycle,
        "variables": variables,
        "locations": frozenset(locations),
        "ordered_locations": tuple(locations),
    }


def grouped_snapshot_evidence(
    group: list[Path],
    input_contract: dict[str, object],
    runtime_contract: dict[str, object],
    profiles: dict[Path, tuple[int, ...]],
) -> dict[str, object]:
    """Require one shared/rank-local snapshot group to cover the exact mesh."""

    evidence = [
        binary_snapshot_evidence(path, input_contract, runtime_contract, profiles)
        for path in group
    ]
    for key in ("time", "cycle", "variables"):
        if any(item[key] != evidence[0][key] for item in evidence[1:]):
            raise ValidationError(f"snapshot sibling {key} metadata disagree")
    expected_locations = expected_logical_locations(input_contract)
    union: set[tuple[int, int, int, int]] = set()
    for item in evidence:
        locations = item["locations"]
        if union.intersection(locations):
            raise ValidationError("snapshot siblings contain duplicate meshblocks")
        union.update(locations)
    if union != expected_locations:
        raise ValidationError("snapshot group does not cover the exact meshblock inventory")
    return {
        "time": evidence[0]["time"],
        "cycle": evidence[0]["cycle"],
        "locations": union,
        "ordered_location_inventory": tuple(
            item["ordered_locations"] for item in evidence
        ),
    }


def restart_parameter_dump(stream: BinaryIO, path: Path) -> tuple[str, int]:
    """Read restart parameter text and return its exact binary payload offset."""

    marker = b"<par_end>\n"
    payload = b""
    while len(payload) < MAX_RESTART_PARAMETER_DUMP_BYTES:
        block = stream.read(min(4096, MAX_RESTART_PARAMETER_DUMP_BYTES - len(payload)))
        if not block:
            break
        payload += block
        end = payload.find(marker)
        if end >= 0:
            stream.seek(end + len(marker))
            try:
                return payload[:end].decode("utf-8"), end + len(marker)
            except UnicodeDecodeError as error:
                raise ValidationError(
                    f"restart parameter dump is not UTF-8 text: {path}"
                ) from error
    raise ValidationError(f"restart parameter dump lacks loadable <par_end>: {path}")


def restart_marker_text(text: str, path: Path) -> str:
    """Read exactly one time/restart_time marker from parameter text."""

    block = ""
    markers = []
    for original in text.splitlines():
        line = original.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            block = line[1:-1].strip()
            continue
        if block == "time" and "=" in line:
            key, value = line.split("=", 1)
            if key.strip() == "restart_time":
                markers.append(value.strip())
    if len(markers) != 1:
        raise ValidationError(
            f"restart parameter dump must contain one time/restart_time marker: {path}"
        )
    return markers[0]


def qualified_restart_abi(command: dict[str, object]) -> dict[str, object]:
    """Return the closed qualified binary-restart ABI for a manifest command."""

    key = (command.get("executable_revision"), command.get("executable_sha256"))
    abi = QUALIFIED_RESTART_BINARY_ABIS.get(key)
    if abi is None:
        raise ValidationError("prepared executable has no qualified binary-restart ABI")
    return abi


def require_restart_region_indices(
    encoded: bytes,
    expected: tuple[int, ...],
    label: str,
    semantic_fields: int = 19,
) -> tuple[int, ...]:
    """Validate every semantically initialized native RegionIndcs field."""

    values = struct.unpack("<19i", encoded)
    if semantic_fields not in (10, 19) or values[:semantic_fields] != expected[:semantic_fields]:
        raise ValidationError(f"{label} indices differ from archived input")
    return values


def restart_rng_evidence(
    payload: bytes, abi: dict[str, object], path: Path, n_updates: int
) -> dict[str, object]:
    """Authenticate lifecycle-valid RNG continuation state and record ignored bytes."""

    native_format = str(abi.get("rng_state_format", ""))
    generator_format = str(abi.get("rng_state_generator_format", ""))
    generator_field_count = require_exact_int(
        abi.get("rng_state_generator_field_count"),
        "qualified RNG generator field count",
    )
    iset_format = str(abi.get("rng_state_iset_format", ""))
    gset_format = str(abi.get("rng_state_conditional_gset_format", ""))
    padding_offset = require_exact_int(
        abi.get("rng_state_padding_offset"), "qualified RNG padding offset"
    )
    padding_size = require_exact_int(
        abi.get("rng_state_padding_size"), "qualified RNG padding size"
    )
    try:
        native_size = struct.calcsize(native_format)
        generator_size = struct.calcsize(generator_format)
        iset_size = struct.calcsize(iset_format)
        gset_size = struct.calcsize(gset_format)
        native_state = struct.unpack(native_format, payload)
    except (struct.error, TypeError, ValueError) as error:
        raise ValidationError("qualified RNG_State ABI is invalid") from error
    iset = native_state[-2]
    gset = native_state[-1]
    idum_payload = payload[:struct.calcsize("<q")]
    gset_payload = payload[padding_offset + padding_size:]
    if (
        type(n_updates) is not int
        or n_updates < 0
        or native_size != len(payload)
        or generator_size + iset_size != padding_offset
        or len(native_state) != generator_field_count + 2
        or padding_offset <= 0
        or padding_size <= 0
        or padding_offset + padding_size + gset_size != len(payload)
        or len(gset_payload) != gset_size
        or iset not in (0, 1)
    ):
        raise ValidationError(f"restart RNG semantic state is invalid: {path}")
    idum = native_state[0]
    if idum < 0:
        if n_updates != 0 or iset != 0:
            raise ValidationError(f"restart RNG lifecycle state is invalid: {path}")
        lifecycle = "pre_update_seed"
        canonical_continuation_payload = idum_payload
        authenticated_field_count = 1
        gset_authentication = "ignored_pre_update_uninitialized"
        ignored_pre_update_state = {
            "idum2_raw_hex": payload[8:16].hex(),
            "iy_raw_hex": payload[16:24].hex(),
            "iv_raw_hex": payload[24:generator_size].hex(),
            "gset_raw_hex": gset_payload.hex(),
            "native_padding_raw_hex": payload[
                padding_offset:padding_offset + padding_size
            ].hex(),
        }
    elif idum > 0:
        if n_updates == 0:
            raise ValidationError(f"restart RNG lifecycle state is invalid: {path}")
        if iset == 1 and not math.isfinite(float(gset)):
            raise ValidationError(f"restart RNG semantic state is invalid: {path}")
        lifecycle = "initialized_continuation"
        canonical_continuation_payload = payload[:padding_offset]
        if iset == 1:
            canonical_continuation_payload += gset_payload
        authenticated_field_count = generator_field_count + 1 + (1 if iset == 1 else 0)
        gset_authentication = (
            "authenticated_exact_finite"
            if iset == 1
            else "ignored_inactive_uninitialized"
        )
        ignored_pre_update_state = None
    else:
        raise ValidationError(f"restart RNG lifecycle state is invalid: {path}")
    return {
        "rng_lifecycle": lifecycle,
        "rng_n_updates": n_updates,
        "rng_canonical_continuation_payload": canonical_continuation_payload,
        "rng_canonical_continuation_state_sha256": hashlib.sha256(
            canonical_continuation_payload
        ).hexdigest(),
        "rng_authenticated_field_count": authenticated_field_count,
        "rng_gset_active": iset == 1,
        "rng_gset_raw_hex": gset_payload.hex(),
        "rng_gset_authentication": gset_authentication,
        "rng_pre_update_ignored_state": ignored_pre_update_state,
        "rng_native_padding_offset": padding_offset,
        "rng_native_padding_size": padding_size,
        "rng_native_padding_hex": payload[
            padding_offset:padding_offset + padding_size
        ].hex(),
    }


def restart_binary_evidence(
    path: Path,
    abi: dict[str, object],
    input_contract: dict[str, object],
    runtime_contract: dict[str, object],
    profiles: dict[Path, tuple[int, ...]],
) -> dict[str, object]:
    """Structurally parse one qualified native restart product."""

    with open_regular_binary(path, "restart product") as (stream, profile):
        remember_profile(profiles, path, profile, "restart product")
        text, parameter_size = restart_parameter_dump(stream, path)
        marker_text = restart_marker_text(text, path)
        parameter_blocks = parse_athinput(
            text.encode("utf-8"), "restart parameter dump"
        )
        require_product_parameter_contract(
            parameter_blocks,
            input_contract,
            runtime_contract,
            "restart parameter dump",
            restart_product=True,
        )
        history_last_time_text = parameter_blocks["output1"]["last_time"]
        history_last_time = float(history_last_time_text)
        header = read_exact(stream, int(abi["mesh_header_size"]), "restart mesh header")
        nmb_total, root_level = struct.unpack_from("<ii", header, 0)
        mesh_geometry = struct.unpack_from("<9d", header, 8)
        mesh_region_contract = str(abi.get("mesh_coarse_region_fields", "initialized"))
        mesh_semantic_fields = (
            10 if mesh_region_contract == "legacy_opaque_uninitialized" else 19
        )
        mesh_region = require_restart_region_indices(
            header[80:156],
            expected_region_indices(
                input_contract, tuple(int(value) for value in input_contract["mesh"])
            ),
            "restart mesh",
            mesh_semantic_fields,
        )
        meshblock_region = require_restart_region_indices(
            header[156:232],
            expected_region_indices(
                input_contract, tuple(int(value) for value in input_contract["meshblock"])
            ),
            "restart meshblock",
        )
        require_close_sequence(
            mesh_geometry,
            expected_mesh_geometry(input_contract),
            f"restart mesh geometry: {path}",
        )
        binary_time = struct.unpack_from(str(abi["mesh_time_format"]), header, 232)[0]
        dt = struct.unpack_from("<d", header, 240)[0]
        ncycle = struct.unpack_from("<i", header, 248)[0]
        if (
            nmb_total != input_contract["expected_meshblocks"]
            or root_level < 0
            or not all(math.isfinite(value) for value in mesh_geometry)
            or not math.isfinite(binary_time)
            or not math.isfinite(dt)
            or dt <= 0.0
            or ncycle < 0
        ):
            raise ValidationError(f"restart mesh header is not loadable: {path}")
        try:
            marker_time = float(marker_text)
        except ValueError as error:
            raise ValidationError(f"restart time marker is not numeric: {path}") from error
        if marker_text == format(binary_time, ".17g"):
            marker_mode = "full_precision"
        elif marker_text == format(binary_time, ".6g"):
            marker_mode = "legacy_default_precision"
        else:
            raise ValidationError(
                f"restart marker does not authenticate binary physical time: {path}"
            )
        if marker_mode not in abi["allowed_marker_modes"] or not math.isfinite(marker_time):
            raise ValidationError(f"restart marker mode is not qualified: {path}")
        location_payload = read_exact(
            stream,
            nmb_total * int(abi["logical_location_size"]),
            "restart logical-location inventory",
        )
        raw_ordered_locations = tuple(
            struct.unpack_from("<4i", location_payload, offset)
            for offset in range(0, len(location_payload), 16)
        )
        if len(set(raw_ordered_locations)) != nmb_total or any(
            value < 0 for location in raw_ordered_locations for value in location
        ):
            raise ValidationError(f"restart logical-location inventory is invalid: {path}")
        if any(location[3] < root_level for location in raw_ordered_locations):
            raise ValidationError(f"restart refinement level is below root level: {path}")
        # Binary outputs store refinement level relative to root_level; native
        # restarts store the absolute level.
        ordered_locations = tuple(
            (location[0], location[1], location[2], location[3] - root_level)
            for location in raw_ordered_locations
        )
        if set(ordered_locations) != expected_logical_locations(input_contract):
            raise ValidationError(f"restart logical-location inventory is invalid: {path}")
        locations = frozenset(ordered_locations)
        cost_payload = read_exact(
            stream, nmb_total * int(abi["cost_size"]), "restart cost inventory"
        )
        costs = struct.unpack(f"<{nmb_total}f", cost_payload)
        if not all(math.isfinite(value) and value >= 0.0 for value in costs):
            raise ValidationError(f"restart cost inventory is invalid: {path}")
        ordered_location_cost_inventory = tuple(zip(ordered_locations, costs))
        metadata_payload = read_exact(
            stream, int(abi["turbulence_metadata_size"]), "restart turbulence metadata"
        )
        metadata_ints = struct.unpack_from("<24i", metadata_payload)
        metadata_reals = struct.unpack_from("<19d", metadata_payload, 96)
        metadata = dict(zip(TURBULENCE_METADATA_INT_FIELDS, metadata_ints))
        metadata.update(zip(TURBULENCE_METADATA_REAL_FIELDS, metadata_reals))
        expected_metadata = expected_turbulence_restart_configuration(input_contract)
        for name in TURBULENCE_METADATA_CONFIGURATION_FIELDS:
            if metadata[name] != expected_metadata[name]:
                raise ValidationError(
                    "restart turbulence configuration differs from frozen input "
                    f"for '{name}': {path}"
                )
        mode_count = metadata_ints[1]
        n_updates = metadata_ints[2]
        if (
            n_updates < 0
            or not all(math.isfinite(value) for value in metadata_reals)
        ):
            raise ValidationError(f"restart turbulence metadata is not qualified: {path}")
        rng_payload = read_exact(stream, int(abi["rng_state_size"]), "restart RNG state")
        rng_evidence = restart_rng_evidence(rng_payload, abi, path, n_updates)
        amplitude_payload = read_exact(
            stream, 6 * mode_count * 8, "restart turbulence amplitudes"
        )
        amplitudes = struct.unpack(f"<{6 * mode_count}d", amplitude_payload)
        if not all(math.isfinite(value) for value in amplitudes):
            raise ValidationError(f"restart turbulence amplitudes are non-finite: {path}")
        injected_work = struct.unpack(
            "<d", read_exact(stream, 8, "restart injected work")
        )[0]
        if not math.isfinite(injected_work):
            raise ValidationError(f"restart injected work is non-finite: {path}")
        diagnostic_count = int(abi["lf_diagnostic_count"])
        lf_diag = struct.unpack(
            f"<{diagnostic_count}d",
            read_exact(stream, diagnostic_count * 8, "restart LF diagnostics"),
        )
        if not all(math.isfinite(value) for value in lf_diag):
            raise ValidationError(f"restart LF diagnostics are non-finite: {path}")
        if any(lf_diag[index] != 0.0 for index in (1, 2, 3, 4, 7)):
            raise ValidationError(f"restart strict LF failure counter is nonzero: {path}")
        if any(
            lf_diag[index] < 0.0 or not lf_diag[index].is_integer()
            for index in COUNT_DIAGNOSTIC_INDICES
        ):
            raise ValidationError(f"restart LF count diagnostics are invalid: {path}")
        if any(lf_diag[index] > lf_diag[8] for index in (9, 10, 11, 12)):
            raise ValidationError(f"restart LF cap diagnostics exceed qface: {path}")
        if input_contract["passive"] and (lf_diag[15] != 0.0 or lf_diag[16] != 0.0):
            raise ValidationError(f"passive restart contains pressure-work feedback: {path}")
        if not input_contract["hardwall"] and lf_diag[17] != 0.0:
            raise ValidationError(f"non-hardwall restart contains hardwall projections: {path}")
        data_size = struct.unpack(
            "<Q", read_exact(stream, 8, "restart variable-data size")
        )[0]
        if data_size != input_contract["restart_data_size"]:
            raise ValidationError(f"restart variable-data size differs from archived input: {path}")
        remaining = profile[6] - stream.tell()
        if remaining <= 0 or remaining % data_size != 0:
            raise ValidationError(f"restart meshblock payload is truncated: {path}")
        local_blocks = remaining // data_size
        require_finite_binary_values(
            stream,
            remaining,
            int(abi["real_size"]),
            "restart meshblock payload",
        )
        if stream.tell() != profile[6]:
            raise ValidationError(f"restart meshblock payload is incomplete: {path}")
    return {
        "parameter_size": parameter_size,
        "time": binary_time,
        "dt": dt,
        "cycle": ncycle,
        "history_last_time": history_last_time,
        "history_last_time_text": history_last_time_text,
        "marker_mode": marker_mode,
        "locations": locations,
        "ordered_location_inventory": ordered_locations,
        "ordered_cost_inventory": costs,
        "ordered_location_cost_inventory": ordered_location_cost_inventory,
        "logical_location_inventory_sha256": hashlib.sha256(
            location_payload
        ).hexdigest(),
        "cost_inventory_sha256": hashlib.sha256(cost_payload).hexdigest(),
        "turbulence_metadata_sha256": hashlib.sha256(metadata_payload).hexdigest(),
        **rng_evidence,
        "turbulence_amplitudes_sha256": hashlib.sha256(
            amplitude_payload
        ).hexdigest(),
        "local_blocks": local_blocks,
        "mode_count": mode_count,
        "n_updates": n_updates,
        "injected_work": injected_work,
        "lf_diagnostics": lf_diag,
        "mesh_region": mesh_region,
        "meshblock_region": meshblock_region,
        "mesh_region_contract": mesh_region_contract,
    }


def grouped_restart_evidence(
    group: list[Path],
    abi: dict[str, object],
    input_contract: dict[str, object],
    runtime_contract: dict[str, object],
    profiles: dict[Path, tuple[int, ...]],
) -> dict[str, object]:
    """Require one shared/rank-local restart group to be loadable and complete."""

    evidence = [
        restart_binary_evidence(
            path, abi, input_contract, runtime_contract, profiles
        )
        for path in group
    ]
    for key in (
        "time",
        "dt",
        "cycle",
        "history_last_time",
        "history_last_time_text",
        "locations",
        "ordered_location_inventory",
        "ordered_cost_inventory",
        "ordered_location_cost_inventory",
        "logical_location_inventory_sha256",
        "cost_inventory_sha256",
        "turbulence_metadata_sha256",
        "rng_lifecycle",
        "rng_n_updates",
        "rng_canonical_continuation_payload",
        "rng_canonical_continuation_state_sha256",
        "rng_authenticated_field_count",
        "rng_gset_active",
        "rng_gset_authentication",
        "turbulence_amplitudes_sha256",
        "mode_count",
        "n_updates",
        "meshblock_region",
        "mesh_region_contract",
    ):
        if any(item[key] != evidence[0][key] for item in evidence[1:]):
            raise ValidationError(f"restart sibling {key} metadata disagree")
    storage = product_storage(group)
    local_blocks = [int(item["local_blocks"]) for item in evidence]
    expected = int(input_contract["expected_meshblocks"])
    if storage == "shared_mpiio" and local_blocks != [expected]:
        raise ValidationError("shared restart does not contain every meshblock payload")
    if storage == "per_rank" and sum(local_blocks) != expected:
        raise ValidationError("rank-local restart payload counts do not cover the mesh")
    injected_work = float(evidence[0]["injected_work"])
    if any(float(item["injected_work"]) != injected_work for item in evidence[1:]):
        raise ValidationError("restart sibling injected-work values disagree")
    if storage == "shared_mpiio":
        lf_diagnostics = tuple(float(value) for value in evidence[0]["lf_diagnostics"])
    else:
        lf_diagnostics = tuple(
            math.fsum(float(item["lf_diagnostics"][index]) for item in evidence)
            for index in range(len(evidence[0]["lf_diagnostics"]))
        )
    return {
        "time": evidence[0]["time"],
        "dt": evidence[0]["dt"],
        "cycle": evidence[0]["cycle"],
        "history_last_time": evidence[0]["history_last_time"],
        "history_last_time_text": evidence[0]["history_last_time_text"],
        "marker_modes": [str(item["marker_mode"]) for item in evidence],
        "locations": evidence[0]["locations"],
        "ordered_location_inventory": evidence[0]["ordered_location_inventory"],
        "ordered_cost_inventory": evidence[0]["ordered_cost_inventory"],
        "ordered_location_cost_inventory": evidence[0][
            "ordered_location_cost_inventory"
        ],
        "replicated_restart_state_sha256": {
            key: evidence[0][key]
            for key in (
                "logical_location_inventory_sha256",
                "cost_inventory_sha256",
                "turbulence_metadata_sha256",
                "rng_canonical_continuation_state_sha256",
                "turbulence_amplitudes_sha256",
            )
        },
        "rng_gset_evidence": [
            {
                "path": str(path),
                "active": item["rng_gset_active"],
                "authentication": item["rng_gset_authentication"],
                "raw_hex": item["rng_gset_raw_hex"],
                "ignored_for_continuation_authentication": not item["rng_gset_active"],
            }
            for path, item in zip(group, evidence)
        ],
        "rng_native_padding_evidence": [
            {
                "path": str(path),
                "offset": item["rng_native_padding_offset"],
                "size_bytes": item["rng_native_padding_size"],
                "hex": item["rng_native_padding_hex"],
                "ignored_for_continuation_authentication": True,
            }
            for path, item in zip(group, evidence)
        ],
        "rng_pre_update_ignored_state_evidence": [
            {
                "path": str(path),
                **item["rng_pre_update_ignored_state"],
                "ignored_for_continuation_authentication": True,
            }
            for path, item in zip(group, evidence)
            if item["rng_pre_update_ignored_state"] is not None
        ],
        "rng_lifecycle": evidence[0]["rng_lifecycle"],
        "rng_n_updates": evidence[0]["rng_n_updates"],
        "rng_authenticated_field_count": evidence[0]["rng_authenticated_field_count"],
        "injected_work": injected_work,
        "lf_diagnostics": lf_diagnostics,
        "mesh_region_contract": evidence[0]["mesh_region_contract"],
    }


def require_exact_identity(
    manifest: dict[str, object],
    inspection: dict[str, object],
    manifest_path: Path,
    result: str,
) -> tuple[dict[str, object], dict[str, object], float]:
    """Bind one schema-4 inspection and recorded result to its prepared manifest."""

    run = manifest.get("run")
    command = manifest.get("command")
    if not isinstance(run, dict) or not isinstance(command, dict):
        raise ValidationError("manifest lacks run or command metadata")
    state = manifest.get("state")
    if state not in {"submitted", "recorded"}:
        raise ValidationError("manifest is not a submitted or recorded segment")
    require_job_id(manifest.get("job_id"), "manifest job ID")
    require_segment_name(run.get("segment"), "manifest segment name")
    if not isinstance(run.get("case_id"), str) or not run["case_id"]:
        raise ValidationError("manifest lacks a nonempty case_id")
    try:
        required_time = float(command["time_tlim_target"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValidationError("manifest lacks a finite prepared target time") from error
    if not math.isfinite(required_time):
        raise ValidationError("manifest prepared target time is not finite")
    expected = {
        "execution_epoch": EXECUTION_EPOCH,
        "job_id": manifest.get("job_id"),
        "case_id": run.get("case_id"),
        "segment": run.get("segment"),
        "manifest": str(manifest_path),
        "required_time": required_time,
    }
    if inspection.get("schema_version") != INSPECTION_SCHEMA_VERSION:
        raise ValidationError("independent validation requires a schema-4 inspection")
    for key, value in expected.items():
        if inspection.get(key) != value:
            raise ValidationError(f"inspection {key} differs from its manifest")
    if inspection.get("restart_time_marker_bypass") is not False:
        raise ValidationError("schema-4 validation forbids restart marker bypass")
    if state == "recorded":
        if manifest.get("scientific_inspection") != inspection:
            raise ValidationError(
                "recorded manifest scientific_inspection differs from retained inspection"
            )
        accounting = manifest.get("accounting")
        if not isinstance(accounting, dict) or accounting.get("result") != result:
            raise ValidationError("recorded accounting result differs from validation result")
        for key in ("job_id", "case_id", "segment", "execution_epoch"):
            if accounting.get(key) != expected[key]:
                raise ValidationError(f"recorded accounting {key} differs from manifest")
    return run, command, required_time


def require_inspection_record(
    inspection: dict[str, object], key: str, expected: object
) -> None:
    """Require one formal inspection record to equal independent capture."""

    if inspection.get(key) != expected:
        raise ValidationError(f"inspection retained {key} differs from disk")


def bind_threshold_policy(
    policy_name: str,
    manifest: dict[str, object],
    result: str,
    required_time: float,
    final_time: float,
) -> dict[str, object]:
    """Bind an immutable named threshold policy to this exact segment."""

    policy = APPROVED_THRESHOLD_POLICIES[policy_name]
    run = manifest["run"]
    actual = {
        "case_id": run["case_id"],
        "job_id": manifest["job_id"],
        "segment": run["segment"],
        "result": result,
        "required_time": required_time,
        "final_time": final_time,
    }
    binding = policy.get("binding")
    if binding is not None and actual != binding:
        raise ValidationError("named threshold policy is not approved for this segment")
    historical_r03 = (
        actual["case_id"] == "R03"
        and actual["job_id"] == "4762472"
        and actual["result"] == "clean_partial"
    )
    if historical_r03 and policy_name != "r03-4762472-clean-partial-v1":
        raise ValidationError("historical R03 clean partial requires its approved policy")
    return policy


def relative_drift(values: list[float]) -> float:
    """Return maximum drift relative to the initial finite state."""

    return max(abs(value - values[0]) for value in values) / max(
        abs(values[0]), NORMALIZATION_FLOOR
    )


def relative_mismatch(left: list[float], right: list[float]) -> float:
    """Return maximum pointwise mismatch relative to either finite state."""

    if len(left) != len(right):
        raise ValidationError("mass histories have different row counts")
    return max(
        abs(left_value - right_value)
        / max(abs(left_value), abs(right_value), NORMALIZATION_FLOOR)
        for left_value, right_value in zip(left, right)
    )


def forcing_closure(
    mhd: dict[str, list[float]],
    user: dict[str, list[float]],
    thresholds: dict[str, float],
) -> dict[str, object]:
    """Require sampled total-energy change to close against forcing work."""

    delta_energy = mhd["tot-E"][-1] - mhd["tot-E"][0]
    delta_force_work = user["force_work"][-1] - user["force_work"][0]
    if delta_force_work == 0.0:
        raise ValidationError("active forcing work did not advance")
    absolute = abs(delta_energy - delta_force_work)
    increment = absolute / max(
        abs(delta_energy), abs(delta_force_work), NORMALIZATION_FLOOR
    )
    state_energy = absolute / max(
        abs(mhd["tot-E"][0]), abs(mhd["tot-E"][-1]), NORMALIZATION_FLOOR
    )
    residuals = {
        "absolute_residual": absolute,
        "increment_normalized_residual": increment,
        "state_energy_normalized_residual": state_energy,
    }
    if not all(math.isfinite(value) for value in residuals.values()):
        raise ValidationError("forcing-work closure produced a non-finite residual")
    if increment >= thresholds["normalized_residual_lt"]:
        raise ValidationError("active forcing-work closure normalized residual is too large")
    failures = [
        name
        for name, value in residuals.items()
        if f"{name}_lt" in thresholds and value >= thresholds[f"{name}_lt"]
    ]
    if failures:
        raise ValidationError(f"forcing-work closure exceeds named policy: {failures}")
    return {
        "claim": "active_delta_energy_closure",
        "delta_energy": delta_energy,
        "delta_force_work": delta_force_work,
        "normalized_residual": increment,
        **residuals,
        "thresholds": thresholds,
    }


def require_meaningful_activity(
    delta: float,
    state_scale: float,
    absolute_gt: float,
    normalized_gt: float,
    label: str,
) -> dict[str, float]:
    """Require an activity increment above named physical and numerical floors."""

    absolute = abs(delta)
    normalized = absolute / max(abs(state_scale), NORMALIZATION_FLOOR)
    if (
        not math.isfinite(absolute)
        or not math.isfinite(normalized)
        or absolute <= absolute_gt
        or normalized <= normalized_gt
    ):
        raise ValidationError(f"{label} did not exceed the named activity policy")
    return {
        "delta": delta,
        "absolute": absolute,
        "state_normalized": normalized,
        "absolute_gt": absolute_gt,
        "state_normalized_gt": normalized_gt,
    }


def values_close(
    actual: float, expected: float, absolute_tolerance: float, relative_tolerance: float
) -> bool:
    """Return whether two finite diagnostics agree within a named tolerance."""

    return (
        math.isfinite(actual)
        and math.isfinite(expected)
        and math.isclose(
            actual,
            expected,
            rel_tol=relative_tolerance,
            abs_tol=absolute_tolerance,
        )
    )


def validate_terminal_restart_history(
    restart: dict[str, object],
    mhd: dict[str, list[float]],
    user: dict[str, list[float]],
    policy: dict[str, object],
) -> dict[str, object]:
    """Bind terminal restart forcing/LF state to synchronized terminal histories."""

    tolerances = dict(policy["restart_history"])
    absolute = float(tolerances["work_absolute_tolerance"])
    relative = float(tolerances["work_relative_tolerance"])
    injected_work = float(restart["injected_work"])
    if not values_close(injected_work, user["force_work"][-1], absolute, relative):
        raise ValidationError("terminal restart injected work differs from user history")
    restart_dt = float(restart["dt"])
    history_dt = mhd["dt"][-1]
    if restart_dt != history_dt:
        raise ValidationError("terminal restart dt differs from terminal histories")
    diagnostics = tuple(float(value) for value in restart["lf_diagnostics"])
    if len(diagnostics) != len(LF_RESTART_HISTORY_COLUMNS):
        raise ValidationError("terminal restart LF diagnostic inventory is invalid")
    residuals = {}
    for index, name in enumerate(LF_RESTART_HISTORY_COLUMNS):
        expected = mhd[name][-1]
        actual = diagnostics[index]
        if index in COUNT_DIAGNOSTIC_INDICES:
            if actual != expected:
                raise ValidationError(
                    f"terminal restart {name} differs from MHD history"
                )
        elif not values_close(actual, expected, absolute, relative):
            raise ValidationError(f"terminal restart {name} differs from MHD history")
        residuals[name] = actual - expected
    return {
        "restart_dt": restart_dt,
        "history_dt": history_dt,
        "dt_residual": restart_dt - history_dt,
        "injected_work_residual": injected_work - user["force_work"][-1],
        "lf_diagnostic_residuals": residuals,
        "tolerances": tolerances,
    }


def validate_terminal_restart_history_schedule(
    restart: dict[str, object],
    history_schedule: dict[str, object],
) -> dict[str, object]:
    """Bind the continuation schedule serialized by the terminal restart."""

    actual = float(restart["history_last_time"])
    expected = float(history_schedule["next_scheduled_time"])
    actual_text = str(restart["history_last_time_text"])
    expected_text = format(expected, ".6g")
    if actual_text != expected_text:
        raise ValidationError("terminal restart history schedule differs from histories")
    return {
        "status": "terminal_restart_continuation_schedule_bound",
        "serialized_history_last_time": actual,
        "serialized_history_last_time_text": actual_text,
        "replayed_history_last_time": expected,
        "expected_serialized_history_last_time_text": expected_text,
        "serialization_contract": "cxx_defaultfloat_six_significant_digits",
    }


def validate_physics(
    mhd: dict[str, list[float]],
    user: dict[str, list[float]],
    abi: dict[str, object],
    input_contract: dict[str, object],
    policy: dict[str, object],
    segment_start: float,
    final_time: float,
    required_time: float,
    segment_result: str,
    continuation: dict[str, object] | None,
) -> dict[str, object]:
    """Apply case-appropriate plasma-physics and conservation gates."""

    require_exact_historical_history_columns(
        mhd, abi.get("mhd_history_columns"), "MHD history"
    )
    require_exact_historical_history_columns(
        user, abi.get("user_history_columns"), "user history"
    )
    if abi.get("normalized_ct_divb_history") != "not_retained":
        raise ValidationError("qualified historical CT-divB history contract is invalid")
    require_columns(
        mhd,
        set(STRICT_LF_FAILURE_COLUMNS)
        | set(LF_REQUIRED_COLUMNS)
        | set(LF_RESTART_HISTORY_COLUMNS)
        | {"time", "tot-E", "mass"},
        "MHD history",
    )
    require_columns(user, {"time", "mass", "hard_vol", "force_work"}, "user history")
    continuation_history_last_time = (
        float(continuation["initial_history_last_time"])
        if isinstance(continuation, dict)
        else None
    )
    mhd_history_schedule = require_history_time_coverage(
        mhd["time"],
        mhd["dt"],
        segment_start,
        final_time,
        required_time,
        segment_result,
        float(input_contract["history_cadence"]),
        "MHD history time",
        continuation_history_last_time,
    )
    user_history_schedule = require_history_time_coverage(
        user["time"],
        user["dt"],
        segment_start,
        final_time,
        required_time,
        segment_result,
        float(input_contract["history_cadence"]),
        "user history time",
        continuation_history_last_time,
    )
    if mhd["time"] != user["time"] or mhd["dt"] != user["dt"]:
        raise ValidationError("MHD and user histories are not exactly synchronized")
    mass_drift = {
        "mhd": relative_drift(mhd["mass"]),
        "user": relative_drift(user["mass"]),
    }
    if any(
        value > float(policy["mass_relative_drift_le"])
        for value in mass_drift.values()
    ):
        raise ValidationError("mass conservation exceeds named policy")
    mass_mismatch = relative_mismatch(mhd["mass"], user["mass"])
    if mass_mismatch > float(policy["mass_relative_mismatch_le"]):
        raise ValidationError("mass-history mismatch exceeds named policy")
    if any(value != 0.0 for value in user["hard_vol"]):
        raise ValidationError("hard_vol is nonzero")
    strict_maxima = {name: max(mhd[name]) for name in STRICT_LF_FAILURE_COLUMNS}
    if any(value != 0.0 for name in STRICT_LF_FAILURE_COLUMNS for value in mhd[name]):
        raise ValidationError("strict LF failure counter is nonzero")
    for name in LF_MONOTONIC_COUNT_COLUMNS:
        values = mhd[name]
        if any(value < 0.0 or not value.is_integer() for value in values):
            raise ValidationError(f"{name} contains an invalid count")
        if any(right < left for left, right in zip(values, values[1:])):
            raise ValidationError(f"{name} is not monotonic")
    for name in ("lf_qprcap", "lf_qpr10", "lf_qpecap", "lf_qpe10"):
        if any(value > qface for value, qface in zip(mhd[name], mhd["lf_qface"])):
            raise ValidationError(f"{name} exceeds lf_qface")
        increments = zip(mhd[name], mhd[name][1:], mhd["lf_qface"], mhd["lf_qface"][1:])
        if any(right - left > qright - qleft for left, right, qleft, qright in increments):
            raise ValidationError(f"{name} increment exceeds lf_qface increment")
    if (
        mhd["lf_nstage"][-1] <= mhd["lf_nstage"][0]
        or mhd["lf_qface"][-1] <= mhd["lf_qface"][0]
    ):
        raise ValidationError("LF stage or qface diagnostics did not advance")
    lf_work_finite = all(
        math.isfinite(value) for name in LF_WORK_COLUMNS for value in mhd[name]
    )
    if not lf_work_finite:
        raise ValidationError("LF work diagnostic is non-finite")
    activity_policy = dict(policy["activity"])
    state_scale = max(abs(mhd["tot-E"][0]), abs(mhd["tot-E"][-1]), NORMALIZATION_FLOOR)
    forcing_activity = require_meaningful_activity(
        user["force_work"][-1] - user["force_work"][0],
        state_scale,
        float(activity_policy["forcing_absolute_gt"]),
        float(activity_policy["forcing_state_normalized_gt"]),
        "forcing work",
    )
    pressure_activity = None
    if input_contract["passive"]:
        if any(value != 0.0 for name in ("lf_cpwrk", "lf_cawrk") for value in mhd[name]):
            raise ValidationError("passive case contains pressure-work feedback")
        closure: dict[str, object] = {
            "claim": "passive_delta_no_energy_closure_claim",
            "delta_force_work": user["force_work"][-1] - user["force_work"][0],
        }
    else:
        pressure_delta = max(
            (mhd[name][-1] - mhd[name][0] for name in ("lf_cpwrk", "lf_cawrk")),
            key=abs,
        )
        pressure_activity = require_meaningful_activity(
            pressure_delta,
            state_scale,
            float(activity_policy["active_pressure_work_absolute_gt"]),
            float(activity_policy["active_pressure_work_state_normalized_gt"]),
            "active pressure work",
        )
        closure = forcing_closure(mhd, user, dict(policy["forcing"]))
    hardwall = mhd["lf_hwproj"]
    if input_contract["finite_limiter"]:
        if any(value != 0.0 for value in hardwall):
            raise ValidationError("finite-limiter case contains hardwall projections")
    else:
        if any(value < 0.0 or not value.is_integer() for value in hardwall):
            raise ValidationError("hardwall projection diagnostic contains an invalid count")
        if any(right < left for left, right in zip(hardwall, hardwall[1:])):
            raise ValidationError("hardwall projection diagnostic is not monotonic")
    return {
        "strict_maxima": strict_maxima,
        "mass_relative_drift": mass_drift,
        "mass_relative_mismatch": mass_mismatch,
        "history_time_schedule_validation": {
            "mhd": mhd_history_schedule,
            "user": user_history_schedule,
        },
        "lf_work_diagnostics_finite": lf_work_finite,
        "activity": {
            "forcing_work": forcing_activity,
            "active_pressure_work": pressure_activity,
        },
        "forcing_work": closure,
        "final_hardwall_projection_count": hardwall[-1],
    }


def require_unchanged_output_inventory(
    output_dir: Path,
    rank_count: int,
    mhd_paths: list[Path],
    user_paths: list[Path],
    snapshot_groups: list[list[Path]],
    restart_groups: list[list[Path]],
) -> None:
    """Rediscover the complete output tree after validation to close inventory races."""

    require_tree_without_symlinks(output_dir, "output tree")
    if sorted(output_dir.glob("*.mhd.hst")) != mhd_paths:
        raise ValidationError("MHD history inventory changed during validation")
    if sorted(output_dir.glob("*.user.hst")) != user_paths:
        raise ValidationError("user history inventory changed during validation")
    if product_groups(output_dir / "bin", ".bin", rank_count, "snapshot") != snapshot_groups:
        raise ValidationError("snapshot inventory changed during validation")
    if product_groups(output_dir / "rst", ".rst", rank_count, "restart") != restart_groups:
        raise ValidationError("restart inventory changed during validation")


def final_recheck_validation_state(
    profiles: dict[Path, tuple[int, ...]],
    output_dir: Path,
    rank_count: int,
    mhd_paths: list[Path],
    user_paths: list[Path],
    snapshot_groups: list[list[Path]],
    restart_groups: list[list[Path]],
) -> None:
    """Close content and inventory races immediately before returning evidence."""

    recheck_profiles(profiles)
    require_unchanged_output_inventory(
        output_dir,
        rank_count,
        mhd_paths,
        user_paths,
        snapshot_groups,
        restart_groups,
    )
    recheck_profiles(profiles)


def validate_segment(args: argparse.Namespace) -> dict[str, object]:
    """Independently validate one retained Stage I segment."""

    profiles: dict[Path, tuple[int, ...]] = {}
    manifest_path = cli_absolute_path(args.manifest)
    inspection_path = cli_absolute_path(
        args.inspection or manifest_path.parent / "segment_inspection.json"
    )
    manifest, manifest_sha256 = load_json(manifest_path, "prepared manifest", profiles)
    inspection, inspection_sha256 = load_json(
        inspection_path, "segment inspection", profiles
    )
    if inspection_path != manifest_path.parent / "segment_inspection.json":
        raise ValidationError("inspection must be the manifest's segment_inspection.json")
    run, command, required_time = require_exact_identity(
        manifest, inspection, manifest_path, args.result
    )
    provenance, input_contract, provenance_artifacts = authenticate_provenance(
        manifest, manifest_path, profiles
    )
    abi = qualified_restart_abi(command)
    continuation = provenance["continuation"]
    segment_start = (
        float(continuation["start_time"]) if isinstance(continuation, dict) else 0.0
    )
    if required_time <= segment_start + 1.0e-12:
        raise ValidationError("prepared target does not advance beyond segment start")
    rank_count = expected_ranks(manifest)
    output_dir = Path(str(provenance["output_dir"]))
    require_tree_without_symlinks(output_dir, "output tree")
    bin_dir = output_dir / "bin"
    restart_dir = output_dir / "rst"
    require_directory_path(bin_dir, "snapshot directory")
    require_directory_path(restart_dir, "restart directory")

    expected_mhd = output_dir / f"{run['run_basename']}.mhd.hst"
    expected_user = output_dir / f"{run['run_basename']}.user.hst"
    mhd_paths = sorted(output_dir.glob("*.mhd.hst"))
    user_paths = sorted(output_dir.glob("*.user.hst"))
    if mhd_paths != [expected_mhd] or user_paths != [expected_user]:
        raise ValidationError(
            "segment histories do not match the exact prepared run_basename"
        )
    mhd_record = captured_file(mhd_paths[0], output_dir, "MHD history", profiles)
    user_record = captured_file(user_paths[0], output_dir, "user history", profiles)
    require_inspection_record(inspection, "mhd_history", mhd_record)
    require_inspection_record(inspection, "user_history", user_record)

    snapshot_groups = product_groups(bin_dir, ".bin", rank_count, "snapshot")
    restart_groups = product_groups(restart_dir, ".rst", rank_count, "restart")
    if not snapshot_groups or not restart_groups:
        raise ValidationError("segment lacks retained snapshots or restarts")
    expected_snapshot_storage = input_contract["storage"]["snapshot"]
    expected_restart_storage = input_contract["storage"]["restart"]
    if any(product_storage(group) != expected_snapshot_storage for group in snapshot_groups):
        raise ValidationError("snapshot storage layout differs from archived input")
    if any(product_storage(group) != expected_restart_storage for group in restart_groups):
        raise ValidationError("restart storage layout differs from archived input")
    snapshot_records = [
        captured_product(group, output_dir, rank_count, "snapshot", profiles)
        for group in snapshot_groups
    ]
    restart_records = [
        captured_product(group, output_dir, rank_count, "restart", profiles)
        for group in restart_groups
    ]
    require_inspection_record(inspection, "snapshots", snapshot_records)
    require_inspection_record(inspection, "restarts", restart_records)

    mhd = parse_history(mhd_paths[0], "MHD history", profiles)
    user = parse_history(user_paths[0], "user history", profiles)
    final_time = mhd.get("time", [None])[-1]
    if not isinstance(final_time, float) or inspection.get("final_time") != final_time:
        raise ValidationError("inspection final time differs from synchronized histories")
    policy = bind_threshold_policy(
        args.policy, manifest, args.result, required_time, final_time
    )
    physics = validate_physics(
        mhd,
        user,
        abi,
        input_contract,
        policy,
        segment_start,
        final_time,
        required_time,
        args.result,
        continuation,
    )
    require_inspection_record(
        inspection, "maximum_strict_failure_counts", physics["strict_maxima"]
    )
    require_inspection_record(
        inspection,
        "final_hardwall_projection_count",
        physics["final_hardwall_projection_count"],
    )

    output_runtime_contract = {
        "run_basename": run["run_basename"],
        "target_time": required_time,
    }
    snapshot_evidence = [
        grouped_snapshot_evidence(
            group, input_contract, output_runtime_contract, profiles
        )
        for group in snapshot_groups
    ]
    snapshot_times = [item["time"] for item in snapshot_evidence]
    snapshot_cycles = [int(item["cycle"]) for item in snapshot_evidence]
    require_inspection_record(inspection, "snapshot_times", snapshot_times)
    if any(value > final_time for value in snapshot_times):
        raise ValidationError("segment retains a future snapshot")
    terminal_snapshots = [
        index for index, value in enumerate(snapshot_times) if value == final_time
    ]
    if len(terminal_snapshots) != 1:
        raise ValidationError("segment does not retain one unique exact terminal snapshot")
    require_time_coverage(
        snapshot_times,
        segment_start,
        final_time,
        float(input_contract["snapshot_cadence"]),
        "snapshot times",
        require_start=continuation is None,
    )
    require_cycle_timeline(snapshot_cycles, "snapshot cycles")
    canonical_locations = snapshot_evidence[terminal_snapshots[0]]["locations"]
    if any(item["locations"] != canonical_locations for item in snapshot_evidence):
        raise ValidationError("snapshot meshblock inventories differ across outputs")
    canonical_snapshot_order = snapshot_evidence[terminal_snapshots[0]][
        "ordered_location_inventory"
    ]
    if any(
        item["ordered_location_inventory"] != canonical_snapshot_order
        for item in snapshot_evidence
    ):
        raise ValidationError(
            "snapshot ordered logical-location inventories differ across outputs"
        )

    restart_evidence = [
        grouped_restart_evidence(
            group, abi, input_contract, output_runtime_contract, profiles
        )
        for group in restart_groups
    ]
    if any(
        item["mesh_region_contract"] != restart_evidence[0]["mesh_region_contract"]
        for item in restart_evidence[1:]
    ):
        raise ValidationError("restart mesh RegionIndcs contracts disagree")
    restart_times = [item["time"] for item in restart_evidence]
    restart_cycles = [int(item["cycle"]) for item in restart_evidence]
    restart_dts = [float(item["dt"]) for item in restart_evidence]
    restart_marker_modes = [item["marker_modes"] for item in restart_evidence]
    require_inspection_record(inspection, "restart_times", restart_times)
    require_inspection_record(
        inspection, "restart_time_marker_modes", restart_marker_modes
    )
    if any(value > final_time for value in restart_times):
        raise ValidationError("segment retains a future restart")
    terminal_restarts = [
        index for index, value in enumerate(restart_times) if value == final_time
    ]
    if len(terminal_restarts) != 1:
        raise ValidationError("segment does not retain one unique exact terminal restart")
    require_time_coverage(
        restart_times,
        segment_start,
        final_time,
        float(input_contract["restart_cadence"]),
        "restart times",
        require_start=continuation is None,
    )
    require_cycle_timeline(restart_cycles, "restart cycles")
    require_matching_output_cycles(snapshot_evidence, restart_evidence)
    if (
        expected_snapshot_storage == "shared_mpiio"
        and expected_restart_storage == "shared_mpiio"
    ):
        require_matching_shared_output_orders(snapshot_evidence, restart_evidence)
    require_continuation_initial_state(
        continuation, snapshot_evidence, restart_evidence
    )
    if any(item["locations"] != canonical_locations for item in restart_evidence):
        raise ValidationError("restart inventory differs from snapshot meshblock inventory")
    if any(
        item["ordered_location_inventory"]
        != restart_evidence[0]["ordered_location_inventory"]
        for item in restart_evidence[1:]
    ):
        raise ValidationError(
            "restart ordered logical-location inventories differ across outputs"
        )
    if any(
        item["ordered_location_cost_inventory"]
        != restart_evidence[0]["ordered_location_cost_inventory"]
        for item in restart_evidence[1:]
    ):
        raise ValidationError(
            "restart ordered logical-location/cost pairs differ across outputs"
        )
    terminal_restart_index = terminal_restarts[0]
    require_inspection_record(
        inspection, "terminal_restart", restart_records[terminal_restart_index]
    )
    require_inspection_record(
        inspection, "terminal_restart_time", restart_times[terminal_restart_index]
    )
    restart_history = validate_terminal_restart_history(
        restart_evidence[terminal_restart_index], mhd, user, policy
    )
    terminal_restart_history_schedule = validate_terminal_restart_history_schedule(
        restart_evidence[terminal_restart_index],
        physics["history_time_schedule_validation"]["mhd"],
    )

    checks = {
        "required_time_reached": final_time >= required_time - 1.0e-10,
        "strict_lf_failure_counters_zero": True,
        "snapshots_retained": True,
        "terminal_snapshot_retained": True,
        "restart_retained": True,
        "terminal_restart_physical_time_matches_final": True,
    }
    require_inspection_record(inspection, "checks", checks)
    if inspection.get("clean_for_continuation") is not True:
        raise ValidationError("segment is not clean for continuation")
    if args.result == "accepted":
        if final_time != required_time:
            raise ValidationError("accepted validation requires the exact prepared endpoint")
        if inspection.get("accepted") is not True:
            raise ValidationError("formal inspection does not accept the exact endpoint")
    else:
        if final_time >= required_time or checks["required_time_reached"]:
            raise ValidationError("clean_partial requires an endpoint short of its target")
        if inspection.get("accepted") is not False:
            raise ValidationError("formal inspection does not describe a clean partial")

    product_files = [
        mhd_record,
        user_record,
        *flatten_product_records(snapshot_records),
        *flatten_product_records(restart_records),
    ]
    validator_path = cli_absolute_path(Path(__file__))
    validator_implementation = validator_implementation_provenance(
        validator_path, profiles
    )
    result = {
        "schema_version": 3,
        "record_type": "stage-i-independent-segment-validation",
        "record_role": "historical_read_only_non_authorizing",
        "authorization_effect": "none",
        "validation_accepted": True,
        "segment_result": args.result,
        "threshold_policy": args.policy,
        "execution_epoch": EXECUTION_EPOCH,
        "job_id": manifest.get("job_id"),
        "case_id": run["case_id"],
        "segment": run["segment"],
        "case_mode": "passive" if input_contract["passive"] else "active",
        "limiter_mode": (
            "finite_collisional" if input_contract["finite_limiter"] else "hardwall"
        ),
        "prepared_target_time": required_time,
        "segment_start_time": segment_start,
        "final_time": final_time,
        "expected_ranks": rank_count,
        "expected_meshblocks": input_contract["expected_meshblocks"],
        "snapshot_times": snapshot_times,
        "snapshot_cycles": snapshot_cycles,
        "snapshot_ordered_location_inventories": [
            item["ordered_location_inventory"] for item in snapshot_evidence
        ],
        "restart_binary_times": restart_times,
        "restart_cycles": restart_cycles,
        "restart_dts": restart_dts,
        "restart_marker_modes": restart_marker_modes,
        "restart_ordered_location_inventories": [
            item["ordered_location_inventory"] for item in restart_evidence
        ],
        "restart_ordered_cost_inventories": [
            item["ordered_cost_inventory"] for item in restart_evidence
        ],
        "restart_ordered_location_cost_inventories": [
            item["ordered_location_cost_inventory"] for item in restart_evidence
        ],
        "restart_replicated_state_sha256": [
            item["replicated_restart_state_sha256"] for item in restart_evidence
        ],
        "restart_rng_native_padding_evidence": [
            item["rng_native_padding_evidence"] for item in restart_evidence
        ],
        "restart_rng_gset_evidence": [
            item["rng_gset_evidence"] for item in restart_evidence
        ],
        "restart_rng_pre_update_ignored_state_evidence": [
            item["rng_pre_update_ignored_state_evidence"] for item in restart_evidence
        ],
        "restart_rng_lifecycle_evidence": [
            {
                "lifecycle": item["rng_lifecycle"],
                "n_updates": item["rng_n_updates"],
                "authenticated_field_count": item["rng_authenticated_field_count"],
                "canonical_continuation_state_sha256": item[
                    "replicated_restart_state_sha256"
                ]["rng_canonical_continuation_state_sha256"],
            }
            for item in restart_evidence
        ],
        "strict_failure_maxima": physics["strict_maxima"],
        "mass_relative_drift": physics["mass_relative_drift"],
        "mass_relative_mismatch": physics["mass_relative_mismatch"],
        "history_time_schedule_validation": physics[
            "history_time_schedule_validation"
        ],
        "terminal_restart_history_schedule_validation": (
            terminal_restart_history_schedule
        ),
        "lf_work_diagnostics_finite": physics["lf_work_diagnostics_finite"],
        "named_activity_validation": physics["activity"],
        "final_hardwall_projection_count": physics[
            "final_hardwall_projection_count"
        ],
        "sampled_history_forcing_work": physics["forcing_work"],
        "terminal_restart_history_validation": restart_history,
        "continuation_origin_validation": continuation,
        "parent_acceptance_validation": (
            continuation["parent_acceptance_validation"]
            if isinstance(continuation, dict)
            else {
                "status": "not_applicable_fresh_segment",
                "authorizing": False,
                "reason": "this fresh segment has no parent acceptance decision",
            }
        ),
        "filesystem_race_validation": {
            "status": "descriptor_profile_bound_with_final_recheck",
            "authorizing": False,
            "reason": (
                "every hashed and parsed retained file is bound to one named inode "
                "profile and rechecked around a final inventory scan, but a read-only "
                "validator cannot freeze later filesystem mutation"
            ),
        },
        "parameter_contract_validation": {
            "status": "closed_qualified_inventory",
            "authorizing": False,
            "reason": (
                "archived input and embedded product parameter inventories are closed "
                "against the qualified Stage I executable contract"
            ),
        },
        "cycle_dt_consistency_validation": {
            "status": "validated",
            "authorizing": False,
            "snapshot_cycles": snapshot_cycles,
            "restart_cycles": restart_cycles,
            "restart_dts": restart_dts,
        },
        "historical_scope_validation": {
            "status": "frozen_historical_executable_validated",
            "authorizing": False,
            "execution_epoch": EXECUTION_EPOCH,
            "executable_revision": command["executable_revision"],
            "executable_sha256": command["executable_sha256"],
            "reason": (
                "this validator accepts only the qualified historical E03 executable "
                "identity and its exact retained output schemas"
            ),
        },
        "independent_authorization": {
            "authorizing": False,
            "status": "non_authorizing",
            "reasons": [
                "the qualified historical E03 executable did not retain normalized CT divB in user histories",
                "face-centered restart fields are not independently evaluated for CT divergence",
                "no independent executable restart-load smoke test is performed",
                "scheduler, accounting, reservation, ledger, and reconciliation authorization is out of scope",
                "this historical read-only evidence cannot authorize launch, continuation, result acceptance, qualification, or publication",
                "historical execution of batch runtime checksum checks is not independently replayed",
                "finite payload decoding does not establish physical bounds or publication-quality statistics",
                "the qualified legacy restart ABI contains opaque uninitialized mesh-level coarse RegionIndcs fields",
                "restart/history binding covers retained forcing and LF diagnostics but not every internal runtime state",
                "restart RNG lifecycle-valid canonical continuation state and forcing-mode state are sibling-consistent, but pre-update/inactive bytes are ignored and stochastic continuation is not replayed",
                "rank-local restart payloads are count-complete but are not independently mapped to logical MeshBlocks",
                "the read-only validator cannot freeze the filesystem after its final state recheck",
            ],
        },
        "normalized_ct_divb_validation": {
            "status": "historically_unavailable",
            "authorizing": False,
            "executable_revision": command["executable_revision"],
            "executable_sha256": command["executable_sha256"],
            "reason": (
                "the qualified historical E03 executable did not retain normalized "
                "CT divB in its exact user-history schema"
            ),
        },
        "restart_face_field_validation": {
            "status": "unavailable",
            "authorizing": False,
            "reason": (
                "native restart values are fully finite-decoded, but face-centered "
                "magnetic fields are not independently interpreted for CT divergence"
            ),
        },
        "restart_load_validation": {
            "status": "unavailable",
            "authorizing": False,
            "reason": "the retained executable is not independently run to load each restart",
        },
        "restart_stochastic_state_validation": {
            "status": "lifecycle_valid_canonical_continuation_state_consistent_ignored_bytes_recorded_not_replayed",
            "authorizing": False,
            "authenticated_field_counts": [
                item["rng_authenticated_field_count"] for item in restart_evidence
            ],
            "lifecycle_rule": (
                "idum<0 requires n_updates=0 and iset=0 and authenticates only idum; "
                "idum>0 requires n_updates>0 and authenticates idum/idum2/iy/iv/iset "
                "plus exact finite gset only when iset=1; idum=0 is invalid"
            ),
            "native_padding_offset": int(abi["rng_state_padding_offset"]),
            "native_padding_size_bytes": int(abi["rng_state_padding_size"]),
            "reason": (
                "lifecycle-valid canonical RNG continuation state and forcing-mode amplitudes "
                "are exact across siblings; pre-update idum2/iy/iv/gset, inactive gset, and "
                "native padding are separately recorded and ignored where the source does not "
                "read them, and stochastic continuation is not replayed"
            ),
        },
        "restart_turbulence_configuration_validation": {
            "status": "exact_frozen_input_source_contract",
            "authorizing": False,
            "configuration_fields": list(TURBULENCE_METADATA_CONFIGURATION_FIELDS),
            "evolving_lifecycle_field": "n_updates",
        },
        "rank_local_restart_mapping_validation": {
            "status": (
                "unavailable"
                if expected_restart_storage == "per_rank"
                else "not_applicable"
            ),
            "authorizing": False,
            "reason": (
                "rank-local payload counts are complete, but the qualified retained "
                "ABI does not independently expose their logical-MeshBlock mapping"
            ),
        },
        "operational_authorization_validation": {
            "status": "out_of_scope",
            "authorizing": False,
            "reason": (
                "scheduler, accounting, reservation, ledger, and reconciliation "
                "records are not independently authorized by this validator; this "
                "record cannot authorize launch, continuation, acceptance, "
                "qualification, or publication"
            ),
        },
        "historical_runtime_checksum_validation": {
            "status": "unavailable",
            "authorizing": False,
            "reason": (
                "the exact batch checksum checks are regenerated and compared, but "
                "their historical execution is not independently replayed"
            ),
        },
        "scientific_interpretation_validation": {
            "status": "limited",
            "authorizing": False,
            "reason": (
                "finite complete payloads and selected histories are checked, but "
                "physical bounds and publication-quality statistics are not"
            ),
        },
        "product_payload_validation": {
            "snapshots": (
                "exact active-region indices, ordered logical inventory, block geometry, "
                "closed qualified parameter contract, expected mhd_w_bcc variables, and "
                "full finite field-payload decode"
            ),
            "restarts": (
                "qualified native-ABI active mesh and complete meshblock RegionIndcs, "
                "exact root geometry and ordered logical/cost inventories, exact replicated "
                "frozen-input-derived turbulence configuration with lifecycle-valid evolving "
                "n_updates, canonical RNG continuation state, and amplitudes, recorded ignored "
                "pre-update RNG fields, inactive RNG gset, and native RNG padding, closed "
                "qualified parameter contract, forcing/LF history binding, and full finite "
                "variable-payload decode"
            ),
            "legacy_mesh_coarse_region_fields": restart_evidence[0][
                "mesh_region_contract"
            ],
        },
        "authenticated_bundle_object_sha256": provenance["bundle_object_sha256"],
        "product_count": len(product_files),
        "product_sha256": {
            str(record["path"]): record["sha256"] for record in product_files
        },
        "product_sizes": {
            str(record["path"]): record["size_bytes"] for record in product_files
        },
        "provenance_sha256": {
            str(record["path"]): record["sha256"] for record in provenance_artifacts
        },
        "manifest_sha256": manifest_sha256,
        "inspection_sha256": inspection_sha256,
        "validator_sha256": validator_implementation["sha256"],
        "validator_revision": validator_implementation["revision"],
        "validator_implementation": {
            **validator_implementation,
            "authorizing": False,
            "launch_source_bundle_coverage": "not_claimed",
            "reason": (
                "the validator's exact committed bytes are revision-bound separately "
                "from the historical launch source bundle"
            ),
        },
    }
    final_recheck_validation_state(
        profiles,
        output_dir,
        rank_count,
        mhd_paths,
        user_paths,
        snapshot_groups,
        restart_groups,
    )
    return result


def main(argv: list[str] | None = None) -> int:
    """Run the independent validator and emit deterministic JSON evidence."""

    try:
        result = validate_segment(parse_args(argv))
    except (OSError, ValidationError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
