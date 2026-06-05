"""CPU regressions for the independent CGL-LF Stage I segment validator."""

from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import struct
import subprocess

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
VALIDATOR_PATH = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_validate_segment.py"
SPEC = importlib.util.spec_from_file_location("stage_i_validator", VALIDATOR_PATH)
assert SPEC is not None and SPEC.loader is not None
VALIDATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VALIDATOR)

EXECUTION_EPOCH = "E03-forcing-policy"
EXECUTABLE_REVISION = "9e07542281e4e6d125582f253df3ad2e3b8b154d"
UTILITY_REVISION = "dbe7e50045bfe0c3a99ce0ba0b21e8735bbf1103"
TEST_EXECUTABLE = b"qualified test executable\n"
TEST_EXECUTABLE_SHA256 = hashlib.sha256(TEST_EXECUTABLE).hexdigest()
TEST_ABI = copy.deepcopy(next(iter(VALIDATOR.QUALIFIED_RESTART_BINARY_ABIS.values())))
ORIGINAL_VALIDATOR_IMPLEMENTATION_PROVENANCE = (
    VALIDATOR.validator_implementation_provenance
)
FAKE_VALIDATOR_REVISION = "f" * 40
FINAL_TIME = 0.31282347945569927
CASE_NAMES = copy.deepcopy(VALIDATOR.CASE_NAMES)
PRODUCTION_CASE_INPUT_CONTRACTS = copy.deepcopy(VALIDATOR.CASE_INPUT_CONTRACTS)
PRODUCTION_COMMON_FROZEN_PARAMETERS = copy.deepcopy(VALIDATOR.COMMON_FROZEN_PARAMETERS)
PRODUCTION_CASE_NODE_PROFILES = copy.deepcopy(VALIDATOR.CASE_NODE_PROFILES)
PRODUCTION_EXPECTED_RANKS_PER_NODE = VALIDATOR.EXPECTED_RANKS_PER_NODE
RETAINED_R03_MANIFEST = Path(
    "/lustre/orion/ast207/proj-shared/dfielding/CGL/runs/mks24-stage-i/"
    "E03-forcing-policy/R03/s00_rankio_t0_t0p5/manifest/prepared_run.json"
)
RETAINED_R03_EXECUTABLE = Path(
    "/lustre/orion/ast207/proj-shared/dfielding/CGL/build/"
    "frontier-hip-9e07542281e4-cpe25.09-cce20-rocm6.4.2/src/athena"
)
RETAINED_R03_EXECUTABLE_SHA256 = (
    "68f243f9204df388b24365ae65a567f6f567dbe422a6d7a43b9fb4a499ef118c"
)
RETAINED_R03_QUALIFICATION = Path(
    "/lustre/orion/ast207/proj-shared/dfielding/CGL/accounting/"
    "mks24_stage_i_E03_forcing_policy_qualification_approval.json"
)
RETAINED_R03_QUALIFICATION_SHA256 = (
    "e3fec9f35da42121b902f41ef752021f375aab38bc34c5ff75ae8780bfbca635"
)
FROZEN_MATRIX_SHA256 = "bf31b88b985d1ad4ffe823108dd7c1132bdfa4d5e4a6abde51f66bb7778415c9"
FROZEN_INPUT_SOURCES = {
    "R02": ("inputs/cgl_lf_paper/cgl_lf_paper_standard_active_alfvenic_beta10.athinput", "c310509aa1638418bab7117427e8cca06e5a651c60e15212e4763baad5101206"),
    "R03": ("inputs/cgl_lf_paper/cgl_lf_paper_standard_active_alfvenic_beta100.athinput", "997f449abb3c2e4d1de50509efffa127c6230e3ecb2632612200010b8fa6c0b0"),
    "R04": ("inputs/cgl_lf_paper/cgl_lf_paper_standard_active_random_beta10.athinput", "571eea2ccec5d069ccb1b49d132ba4c5b8bdac7ce15ee8d8a04c92327c453137"),
    "R05": ("inputs/cgl_lf_paper/cgl_lf_paper_standard_active_random_beta100.athinput", "6527a2072105d515287701c904c3af2cebb07e10ab0b45c3fb41f29d3467622c"),
    "R06": ("inputs/cgl_lf_paper/cgl_lf_paper_standard_passive_alfvenic_beta10.athinput", "c6e038c198b23bf20a7cf1e2a90fa82e83e5ec544492ddc64e1cf2bb535a38ae"),
    "R07": ("inputs/cgl_lf_paper/cgl_lf_paper_standard_passive_alfvenic_beta100.athinput", "72ed6a3b38342d00855f34a7c7eb67102b44cd22b6c12ff2448ddb6f141e5145"),
    "R08": ("inputs/cgl_lf_paper/cgl_lf_paper_standard_passive_random_beta10.athinput", "3535aca5e3d47dc383e353cff262d79c56cf3c083a3ac6024d04ddaa5d353780"),
    "R09": ("inputs/cgl_lf_paper/cgl_lf_paper_standard_passive_random_beta100.athinput", "2ec05b247d917259c8306911cfd31d6be4e1c2772d6cdea544db883edfc35f2c"),
    "R10": ("inputs/cgl_lf_paper/cgl_lf_paper_compressive_active_random_beta1.athinput", "93d7ed9b0846019f59bce34feeb9b22b098afffdfc2dbaa59aaf5913655cd705"),
    "R11": ("inputs/cgl_lf_paper/cgl_lf_paper_compressive_active_random_beta100_sonic.athinput", "ab0dba23f80ea2d6c173ecbb724dbcb332b13d3b1751b8776ad05313f9ccaf6c"),
    "R12": ("inputs/cgl_lf_paper/cgl_lf_paper_heat_flux_beta10_strong.athinput", "98ddea4b4f7fec18cc40abdbf5f7c8ba5b583a91f84f4dae00e1411f23d42e7c"),
    "R13": ("inputs/cgl_lf_paper/cgl_lf_paper_heat_flux_beta10_weak.athinput", "a190129ee34c46a4fe83a392a19120a58b9a9a4064742ac58511441370c59c35"),
    "R14": ("inputs/cgl_lf_paper/cgl_lf_paper_nulim_beta100_20.athinput", "2b8d5837f8a7f3070ca2ef56f8b44e53d048839918185a8eef155cb084736578"),
    "R15": ("inputs/cgl_lf_paper/cgl_lf_paper_nulim_beta100_200.athinput", "9a698a60bef4c558ccee4635d69c3acbf3fea478401bd633d09bbc0526d943d8"),
    "R16": ("inputs/cgl_lf_paper/cgl_lf_paper_scale_separation_beta10_nperp96.athinput", "c0ac4b54248e8f8dfb0f5fd34c0cfb4414b5330529cbf2836961c5277af3f2d1"),
    "R17": ("inputs/cgl_lf_paper/cgl_lf_paper_scale_separation_beta10_nperp384.athinput", "cc1092404b82129f807308a64f7585a6da31f45f41f1d2263acad0c8d30a7e04"),
}


def retained_r03_executable_is_intact() -> bool:
    """Return whether canonical R03 retains its exact prepared executable."""

    return (
        RETAINED_R03_EXECUTABLE.is_file()
        and hashlib.sha256(RETAINED_R03_EXECUTABLE.read_bytes()).hexdigest()
        == RETAINED_R03_EXECUTABLE_SHA256
    )


def retained_r03_provenance_is_intact() -> bool:
    """Return whether mutable canonical R03 launch artifacts remain exact."""

    return (
        retained_r03_executable_is_intact()
        and RETAINED_R03_QUALIFICATION.is_file()
        and hashlib.sha256(RETAINED_R03_QUALIFICATION.read_bytes()).hexdigest()
        == RETAINED_R03_QUALIFICATION_SHA256
    )


def use_production_contracts(monkeypatch) -> None:
    """Restore production contracts inside retained canonical compatibility tests."""

    monkeypatch.setattr(
        VALIDATOR, "CASE_INPUT_CONTRACTS", copy.deepcopy(PRODUCTION_CASE_INPUT_CONTRACTS)
    )
    monkeypatch.setattr(
        VALIDATOR,
        "COMMON_FROZEN_PARAMETERS",
        copy.deepcopy(PRODUCTION_COMMON_FROZEN_PARAMETERS),
    )
    monkeypatch.setattr(
        VALIDATOR, "CASE_NODE_PROFILES", copy.deepcopy(PRODUCTION_CASE_NODE_PROFILES)
    )
    monkeypatch.setattr(
        VALIDATOR, "EXPECTED_RANKS_PER_NODE", PRODUCTION_EXPECTED_RANKS_PER_NODE
    )


def qualify_fixture_layout(monkeypatch, layout: str) -> None:
    """Qualify one fixture-only shared or rank-local retained-product layout."""

    common = copy.deepcopy(VALIDATOR.COMMON_FROZEN_PARAMETERS)
    per_rank = str(layout == "per_rank").lower()
    common[("output2", "single_file_per_rank")] = per_rank
    common[("output3", "single_file_per_rank")] = per_rank
    monkeypatch.setattr(VALIDATOR, "COMMON_FROZEN_PARAMETERS", common)


MHD_COLUMNS = VALIDATOR.HISTORICAL_MHD_HISTORY_COLUMNS
USER_COLUMNS = VALIDATOR.HISTORICAL_USER_HISTORY_COLUMNS


def fake_validator_implementation_provenance(path: Path, profiles) -> dict[str, object]:
    """Return revision-shaped fixture provenance while retaining live-byte races."""

    digest, size = VALIDATOR.sha256_regular_file(path, "validator source", profiles)
    return {
        "path": str(path),
        "repository_path": VALIDATOR.VALIDATOR_REPOSITORY_PATH,
        "revision": FAKE_VALIDATOR_REVISION,
        "sha256": digest,
        "size_bytes": size,
        "committed": True,
    }


@pytest.fixture(autouse=True)
def qualify_tiny_closed_contract(monkeypatch):
    """Use production parsers with one tiny exact-contract fixture mesh."""

    monkeypatch.setattr(
        VALIDATOR,
        "validator_implementation_provenance",
        fake_validator_implementation_provenance,
    )
    monkeypatch.setitem(
        VALIDATOR.QUALIFIED_RESTART_BINARY_ABIS,
        (EXECUTABLE_REVISION, TEST_EXECUTABLE_SHA256),
        TEST_ABI,
    )
    contracts = copy.deepcopy(VALIDATOR.CASE_INPUT_CONTRACTS)
    for case_id, contract in contracts.items():
        contracts[case_id] = ("4x1x1", *contract[1:])
    monkeypatch.setattr(VALIDATOR, "CASE_INPUT_CONTRACTS", contracts)
    monkeypatch.setattr(VALIDATOR, "EXPECTED_RANKS_PER_NODE", 2)
    monkeypatch.setattr(
        VALIDATOR,
        "CASE_NODE_PROFILES",
        {case_id: frozenset({1, 2}) for case_id in VALIDATOR.APPROVED_CASE_IDS},
    )
    common = copy.deepcopy(VALIDATOR.COMMON_FROZEN_PARAMETERS)
    common.update({
        ("mesh", "nghost"): "1",
        ("meshblock", "nx1"): "1",
        ("meshblock", "nx2"): "1",
        ("meshblock", "nx3"): "1",
    })
    monkeypatch.setattr(VALIDATOR, "COMMON_FROZEN_PARAMETERS", common)


def write_json(path: Path, value: object) -> None:
    """Write one deterministic JSON fixture."""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def duplicate_json_key(path: Path, key: str) -> None:
    """Inject one nested-parser-visible duplicate authority key."""

    lines = path.read_text().splitlines()
    index = next(
        index
        for index, line in enumerate(lines)
        if line.lstrip().startswith(f'"{key}":')
    )
    original = lines[index]
    lines[index] = original.rstrip().rstrip(",") + ","
    lines.insert(index + 1, original)
    path.write_text("\n".join(lines) + "\n")


def sha256(path: Path) -> str:
    """Hash one fixture file."""

    return hashlib.sha256(path.read_bytes()).hexdigest()


def frozen_source(path: str) -> bytes:
    """Read one exact object from the frozen qualified executable revision."""

    result = VALIDATOR.git_run([
        "-C",
        str(REPOSITORY),
        "show",
        f"{EXECUTABLE_REVISION}:{path}",
    ])
    assert result.returncode == 0, result.stderr.decode(errors="replace")
    return result.stdout


def file_record(path: Path) -> dict[str, object]:
    """Describe one fixture file independently of the validator."""

    payload = path.read_bytes()
    return {
        "path": str(path),
        "size_bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def product_record(paths: list[Path], storage: str) -> dict[str, object]:
    """Build one schema-4 retained product record."""

    members = [file_record(path) for path in paths]
    representative = members[0].copy()
    representative["storage"] = storage
    if storage == "per_rank":
        representative["rank_files"] = members
    return representative


def history_header(columns: tuple[str, ...]) -> str:
    """Return one contiguous Athena history label header."""

    return "# " + " ".join(
        f"[{index}]={name}" for index, name in enumerate(columns)
    ) + "\n"


def write_history(path: Path, columns: tuple[str, ...], rows: list[tuple[float, ...]]) -> None:
    """Write one labeled finite history fixture."""

    path.write_text(
        history_header(columns)
        + "\n".join(" ".join(str(value) for value in row) for row in rows)
        + "\n"
    )


def append_history_column(path: Path, name: str, value: float) -> None:
    """Append one well-formed but historically unqualified history column."""

    lines = path.read_text().splitlines()
    header_index = next(
        index for index, line in enumerate(lines) if VALIDATOR.HISTORY_LABEL_PATTERN.search(line)
    )
    count = len(VALIDATOR.HISTORY_LABEL_PATTERN.findall(lines[header_index]))
    lines[header_index] += f" [{count}]={name}"
    for index, line in enumerate(lines):
        if line and not line.startswith("#"):
            lines[index] += f" {value}"
    path.write_text("\n".join(lines) + "\n")


def write_histories(
    fixture: dict[str, object],
    *,
    mass_end: float = 2.0,
    user_mass_end: float | None = None,
    energy_start: float = 10.0,
    energy_end: float = 10.25,
    hard_vol_end: float = 0.0,
    passive_pressure_feedback: bool = False,
    cpwork_end: float | None = None,
    cawork_end: float | None = None,
    force_work_end: float | None = None,
    nstage_start: float = 0.0,
    nstage_end: float = 2.0,
    strict_failure_end: float = 0.0,
    qface_start: float = 0.0,
    qface_end: float = 6.0,
    mirror_start: float = 0.0,
    mirror_end: float = 0.0,
    firehose_start: float = 0.0,
    firehose_end: float = 0.0,
    cap_start: float = 0.0,
    cap_end: float = 2.0,
    qprwrk_end: float = 0.1,
    qpewrk_end: float = 0.2,
    hardwall_start: float = 0.0,
    hardwall_count: float | None = None,
) -> None:
    """Write case-appropriate MHD and user histories."""

    passive = bool(fixture["passive"])
    hardwall = bool(fixture["hardwall"])
    final_time = float(fixture["final_time"])
    hw = float(2.0 if hardwall else 0.0) if hardwall_count is None else hardwall_count
    if cpwork_end is None:
        cpwork_end = -0.05 if (not passive or passive_pressure_feedback) else 0.0
    if cawork_end is None:
        cawork_end = -0.04 if (not passive or passive_pressure_feedback) else 0.0
    times = [0.0]
    while times[-1] + 0.02 < final_time:
        times.append(round(times[-1] + 0.02, 14))
    if times[-1] != final_time:
        times.append(final_time)
    mhd_rows = []
    user_rows = []
    force_work = (
        (0.123 if passive else 0.25) if force_work_end is None else force_work_end
    )
    synchronized_user_mass = mass_end if user_mass_end is None else user_mass_end
    for index, time in enumerate(times):
        fraction = time / final_time if final_time else 0.0
        terminal = index == len(times) - 1
        strict_failure = strict_failure_end if terminal else 0.0
        mhd_values = {name: 0.0 for name in MHD_COLUMNS}
        mhd_values.update({
            "time": time,
            "dt": 1.0e-3,
            "mass": 2.0 + fraction * (mass_end - 2.0),
            "tot-E": energy_start + fraction * (energy_end - energy_start),
            "lf_nstage": nstage_end if terminal else nstage_start,
            "lf_dfloor": strict_failure,
            "lf_pfloor": strict_failure,
            "lf_nonfin": strict_failure,
            "lf_nonpos": strict_failure,
            "lf_mirror": mirror_end if terminal else mirror_start,
            "lf_firehs": firehose_end if terminal else firehose_start,
            "lf_hardbd": strict_failure,
            "lf_qface": qface_end if terminal else qface_start,
            "lf_qprcap": cap_end if terminal else cap_start,
            "lf_qpr10": cap_end if terminal else cap_start,
            "lf_qpecap": cap_end if terminal else cap_start,
            "lf_qpe10": cap_end if terminal else cap_start,
            "lf_qprwrk": fraction * qprwrk_end,
            "lf_qpewrk": fraction * qpewrk_end,
            "lf_hwproj": hw if terminal else hardwall_start,
            "lf_cpwrk": fraction * cpwork_end,
            "lf_cawrk": fraction * cawork_end,
        })
        user_values = {name: 0.0 for name in USER_COLUMNS}
        user_values.update({
            "time": time,
            "dt": 1.0e-3,
            "volume": 2.0,
            "mass": 2.0 + fraction * (synchronized_user_mass - 2.0),
            "hard_vol": fraction * hard_vol_end,
            "force_work": fraction * force_work,
        })
        mhd_rows.append(tuple(mhd_values[name] for name in MHD_COLUMNS))
        user_rows.append(tuple(user_values[name] for name in USER_COLUMNS))
    write_history(Path(fixture["mhd"]), MHD_COLUMNS, mhd_rows)
    fixture["hardwall_history_end"] = hw
    write_history(Path(fixture["user"]), USER_COLUMNS, user_rows)


def input_text(case_id: str, layout: str, hardwall: bool) -> str:
    """Return one closed-contract archived Stage I input."""

    case = {
        "name": CASE_NAMES[case_id],
        "input": VALIDATOR.CASE_INPUT_PATHS[case_id],
        "resolution": "4x1x1",
    }
    expected, _ = VALIDATOR.frozen_parameters_for_case(case_id, case)
    expected = dict(expected)
    expected[("mhd", "nscalars")] = "0"
    expected[("time", "tlim")] = "10.0"
    if hardwall:
        expected[("mhd", "limiter_hardwall")] = "true"
    else:
        expected.pop(("mhd", "limiter_hardwall"), None)
    per_rank = layout == "per_rank"
    expected[("output2", "single_file_per_rank")] = str(per_rank).lower()
    expected[("output3", "single_file_per_rank")] = str(per_rank).lower()
    order = (
        "job", "mesh", "meshblock", "time", "mhd", "problem", "turb_driving",
        "output1", "output2", "output3",
    )
    lines = []
    for block in order:
        lines.append(f"<{block}>")
        for (candidate, key), value in expected.items():
            if candidate == block:
                lines.append(f"{key} = {value}")
        lines.append("")
    return "\n".join(lines)


@pytest.mark.parametrize(
    ("case_id", "driving_type", "projection_policy", "tcorr"),
    [
        ("R03", 1, 2, 2.0),
        ("R11", 0, 1, 0.2),
    ],
)
def test_turbulence_restart_configuration_matches_frozen_source_contract(
    case_id, driving_type, projection_policy, tcorr
):
    blocks = VALIDATOR.parse_athinput(
        input_text(case_id, "per_rank", True).encode(), "fixture archived input"
    )
    domain = tuple(
        (
            float(blocks["mesh"][f"x{axis}min"]),
            float(blocks["mesh"][f"x{axis}max"]),
        )
        for axis in range(1, 4)
    )
    expected = {
        "version": 3,
        "mode_count": 11,
        "nlow": 1,
        "nhigh": 3,
        "driving_type": driving_type,
        "min_kx": 0,
        "max_kx": 3,
        "min_ky": 0,
        "max_ky": 3,
        "min_kz": 0,
        "max_kz": 3,
        "use_npeak": 0,
        "turb_flag": 2,
        "tile_nx": 1,
        "tile_ny": 1,
        "tile_nz": 1,
        "normalization": 0,
        "localization": 0,
        "spectrum": 1,
        "projection_policy": projection_policy,
        "physical_k_shell": 1,
        "isotropic_power_spectrum": 1,
        "record_injected_work": 1,
        "tcorr": tcorr,
        "dt_update": 0.01,
        "dedt": 0.32,
        "accel_rms": 0.0,
        "sol_fraction": 1.0,
        "kpeak": 12.5664,
        "npeak": 0.0,
        "expo": 2.0,
        "exp_prp": 1.66667,
        "exp_prl": 0.0,
        "tdriv_duration": 3.4028234663852886e38,
        "tdriv_start": 0.0,
        "sigma_x1": -1.0,
        "sigma_x2": -1.0,
        "sigma_x3": -1.0,
        "center_x1": 0.0,
        "center_x2": 0.0,
        "center_x3": 0.0,
        "k_shell_unit": 3.141592653589793,
    }
    assert VALIDATOR.expected_turbulence_restart_configuration({
        "archived_parameter_blocks": blocks,
        "domain": domain,
    }) == expected


def product_parameter_text(
    parameters: str,
    run_basename: str,
    target_time: float,
    restart_time: float | None,
) -> str:
    """Return one production-shaped closed runtime parameter dump."""

    blocks = VALIDATOR.parse_athinput(parameters.encode(), "fixture archived input")
    for (block, key), value in VALIDATOR.QUALIFIED_PRODUCT_PARAMETER_ADDITIONS.items():
        if (block, key) == ("time", "restart_time") and restart_time is None:
            continue
        blocks.setdefault(block, {})[key] = "0" if value is None else value
    blocks["job"]["basename"] = run_basename
    blocks["time"]["tlim"] = format(target_time, ".17g")
    if restart_time is not None:
        blocks["time"]["restart_time"] = format(restart_time, ".17g")
    for block, key in VALIDATOR.PRODUCT_NONNEGATIVE_INTEGER_PARAMETERS:
        blocks[block][key] = "0"
    for block, key in VALIDATOR.PRODUCT_NONNEGATIVE_FINITE_PARAMETERS:
        if key not in blocks.get(block, {}):
            continue
        if (block, key) != ("time", "restart_time"):
            blocks[block][key] = "0"
    order = (
        "job", "mesh", "meshblock", "time", "mhd", "problem", "turb_driving",
        "output1", "output2", "output3", "coord", "mesh_refinement",
    )
    lines = []
    for block in order:
        lines.append(f"<{block}>")
        lines.extend(f"{key} = {value}" for key, value in blocks[block].items())
        lines.append("")
    return "\n".join(lines)


def logical_locations() -> list[tuple[int, int, int, int]]:
    """Return the exact four-block fixture mesh inventory."""

    return [(index, 0, 0, 0) for index in range(4)]


def fixture_cycle(time: float) -> int:
    """Return one deterministic cycle shared by coincident fixture products."""

    cycle = int(round(time * 1.0e6))
    assert 0 <= cycle <= 2_147_483_647
    return cycle


def rank_partitions(rank_count: int) -> list[list[tuple[int, int, int, int]]]:
    """Distribute the fixture mesh inventory over a complete rank set."""

    locations = logical_locations()
    return [locations[rank::rank_count] for rank in range(rank_count)]


def write_snapshot(
    path: Path,
    time: float,
    cycle: int,
    locations: list[tuple[int, int, int, int]],
    parameters: str,
    run_basename: str,
    target_time: float,
) -> None:
    """Write one structurally complete Athena binary snapshot."""

    path.parent.mkdir(parents=True, exist_ok=True)
    runtime_parameters = product_parameter_text(
        parameters, run_basename, target_time, None
    )
    parameter_dump = (runtime_parameters + "\n<par_end>\n").encode()
    header = (
        "Athena binary output version=1.1\n"
        "  size of preheader=5\n"
        f"  time={format(time, '.17g')}\n"
        f"  cycle={cycle}\n"
        "  size of location=8\n"
        "  size of variable=4\n"
        f"  number of variables={len(VALIDATOR.EXPECTED_SNAPSHOT_VARIABLES)}\n"
        f"  variables:  {' '.join(VALIDATOR.EXPECTED_SNAPSHOT_VARIABLES)}  \n"
        f"  header offset={len(parameter_dump)}\n"
    ).encode()
    blocks = bytearray()
    for location in locations:
        blocks.extend(struct.pack("<10i", 1, 1, 0, 0, 0, 0, *location))
        x1 = float(location[0]) / 4.0
        blocks.extend(struct.pack("<6d", x1, x1 + 0.25, 0.0, 1.0, 0.0, 2.0))
        blocks.extend(
            struct.pack(
                f"<{len(VALIDATOR.EXPECTED_SNAPSHOT_VARIABLES)}f",
                *([1.0] * len(VALIDATOR.EXPECTED_SNAPSHOT_VARIABLES)),
            )
        )
    path.write_bytes(header + parameter_dump + blocks)


def region_indices(nghost: int, shape: tuple[int, int, int]) -> bytes:
    """Encode every initialized native RegionIndcs field."""

    starts = (nghost, 0, 0)
    ends = tuple(start + size - 1 for start, size in zip(starts, shape))
    coarse = tuple(max(1, size // 2) for size in shape)
    coarse_ends = tuple(start + size - 1 for start, size in zip(starts, coarse))
    return struct.pack(
        "<19i",
        nghost,
        *shape,
        starts[0],
        ends[0],
        starts[1],
        ends[1],
        starts[2],
        ends[2],
        *coarse,
        starts[0],
        coarse_ends[0],
        starts[1],
        coarse_ends[1],
        starts[2],
        coarse_ends[2],
    )


def restart_data_size() -> int:
    """Return the qualified fixture MHD restart block size."""

    nout = (3, 3, 3)
    return (
        6 * 3 * 3 * 3 * 8
        + 4 * 3 * 3 * 8
        + 3 * 4 * 3 * 8
        + 3 * 3 * 4 * 8
    )


def write_restart(
    path: Path,
    time: float,
    cycle: int,
    marker_mode: str,
    local_blocks: int,
    parameters: str,
    run_basename: str,
    target_time: float,
    *,
    passive: bool,
    hardwall: bool,
    diagnostic_divisor: int = 1,
) -> None:
    """Write one structurally complete qualified native restart."""

    if marker_mode == "full_precision":
        marker = format(time, ".17g")
    elif marker_mode == "legacy_default_precision":
        marker = format(time, ".6g")
    else:
        marker = marker_mode
    path.parent.mkdir(parents=True, exist_ok=True)
    restart_parameters = product_parameter_text(
        parameters, run_basename, target_time, time
    )
    restart_parameters = restart_parameters.replace(
        f"restart_time = {format(time, '.17g')}",
        f"restart_time = {marker}",
        1,
    )
    parameter_dump = (restart_parameters + "\n<par_end>\n").encode()
    header = bytearray()
    header.extend(struct.pack("<ii", 4, 0))
    header.extend(
        struct.pack("<9d", 0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 0.25, 1.0, 2.0)
    )
    header.extend(region_indices(1, (4, 1, 1)))
    header.extend(region_indices(1, (1, 1, 1)))
    header.extend(struct.pack("<ddi", time, 0.01, cycle))
    assert len(header) == 252
    locations = b"".join(struct.pack("<4i", *location) for location in logical_locations())
    costs = struct.pack("<4f", 1.0, 2.0, 3.0, 4.0)
    parameter_blocks = VALIDATOR.parse_athinput(
        parameters.encode(), "fixture archived input"
    )
    domain = tuple(
        (
            float(parameter_blocks["mesh"][f"x{axis}min"]),
            float(parameter_blocks["mesh"][f"x{axis}max"]),
        )
        for axis in range(1, 4)
    )
    metadata_configuration = VALIDATOR.expected_turbulence_restart_configuration({
        "archived_parameter_blocks": parameter_blocks,
        "domain": domain,
    })
    n_updates = 0 if time == 0.0 else max(1, cycle)
    metadata_ints = [
        n_updates if name == "n_updates" else int(metadata_configuration[name])
        for name in VALIDATOR.TURBULENCE_METADATA_INT_FIELDS
    ]
    metadata_reals = [
        float(metadata_configuration[name])
        for name in VALIDATOR.TURBULENCE_METADATA_REAL_FIELDS
    ]
    metadata = struct.pack("<24i19d", *metadata_ints, *metadata_reals)
    assert len(metadata) == 248
    if n_updates == 0:
        rng_state = struct.pack("<35qi4xd", -1, *([0] * 34), 0, float("nan"))
    else:
        rng_state = struct.pack(
            "<35qi4xd", *tuple(range(101, 136)), 0, float("nan")
        )
    assert len(rng_state) == 296
    lf_diag = [0.0] * 18
    lf_diag[0] = 2.0 / diagnostic_divisor if time else 0.0
    lf_diag[8] = 6.0 / diagnostic_divisor if time else 0.0
    for index in (9, 10, 11, 12):
        lf_diag[index] = 2.0 / diagnostic_divisor if time else 0.0
    lf_diag[13] = 0.1 / diagnostic_divisor if time else 0.0
    lf_diag[14] = 0.2 / diagnostic_divisor if time else 0.0
    lf_diag[15] = 0.0 if passive else (-0.05 / diagnostic_divisor if time else 0.0)
    lf_diag[16] = 0.0 if passive else (-0.04 / diagnostic_divisor if time else 0.0)
    lf_diag[17] = 2.0 / diagnostic_divisor if hardwall and time else 0.0
    data_size = restart_data_size()
    payload = (
        parameter_dump
        + header
        + locations
        + costs
        + metadata
        + rng_state
        + struct.pack(f"<{6 * metadata_ints[1]}d", *([0.0] * (6 * metadata_ints[1])))
        + struct.pack("<d", (0.123 if passive else 0.25) if time else 0.0)
        + struct.pack("<18d", *lf_diag)
        + struct.pack("<Q", data_size)
        + b"\0" * (data_size * local_blocks)
    )
    path.write_bytes(payload)


def effective_marker_mode(time: float, requested: str) -> str:
    """Return the mode authenticated by the exact marker text."""

    marker = (
        format(time, ".6g")
        if requested == "legacy_default_precision"
        else format(time, ".17g")
    )
    return "full_precision" if marker == format(time, ".17g") else requested


def finalized_batch_script(
    manifest: dict[str, object], manifest_path: Path
) -> tuple[str, str]:
    """Return one exact regenerated production launch script."""

    text = VALIDATOR.regenerated_batch_script(manifest, manifest_path)
    digest = hashlib.sha256(text.encode()).hexdigest()
    return text.replace(
        VALIDATOR.BATCH_SCRIPT_DIGEST_PLACEHOLDER, digest, 1
    ), digest


def refresh_batch_script(fixture: dict[str, object]) -> None:
    """Regenerate the retained launch script after an intentional manifest edit."""

    manifest_path = Path(fixture["manifest"])
    manifest = json.loads(manifest_path.read_text())
    batch_text, batch_digest = finalized_batch_script(manifest, manifest_path)
    Path(manifest["paths"]["batch_script"]).write_text(batch_text)
    manifest["command"]["batch_script_sha256"] = batch_digest
    write_json(manifest_path, manifest)


def mutate_manifest_value(
    fixture: dict[str, object], keys: tuple[str, ...], value: object
) -> None:
    """Replace one nested manifest value without regenerating launch evidence."""

    manifest_path = Path(fixture["manifest"])
    manifest = json.loads(manifest_path.read_text())
    target = manifest
    for key in keys[:-1]:
        target = target[key]
    target[keys[-1]] = value
    write_json(manifest_path, manifest)


def bind_fixture_artifact(
    fixture: dict[str, object], label: str, path: Path, command_key: str
) -> None:
    """Update retained manifest and fake-bundle authority for one fixture artifact."""

    manifest_path = Path(fixture["manifest"])
    manifest = json.loads(manifest_path.read_text())
    manifest["command"][command_key] = sha256(path)
    write_json(manifest_path, manifest)
    fixture["bundle_blobs"][label] = path.read_bytes()
    refresh_batch_script(fixture)


def replace_bound_input(
    fixture: dict[str, object], original: str, replacement: str
) -> None:
    """Mutate one archived-input parameter while preserving fake provenance."""

    path = Path(fixture["run_dir"]) / "manifest" / "submitted_input.athinput"
    text = path.read_text()
    assert text.count(original) == 1
    path.write_text(text.replace(original, replacement, 1))
    bind_fixture_artifact(fixture, "input", path, "input_sha256")


def write_product_group(
    *,
    directory: Path,
    basename: str,
    index: int,
    suffix: str,
    layout: str,
    rank_count: int,
    writer,
) -> list[Path]:
    """Write one shared or complete rank-local product group."""

    filename = f"{basename}.{index:05d}.{suffix}"
    if layout == "shared_mpiio":
        path = directory / filename
        writer(path, logical_locations(), 4)
        return [path]
    paths = []
    for rank, locations in enumerate(rank_partitions(rank_count)):
        path = directory / f"rank_{rank:08d}" / filename
        writer(path, locations, len(locations))
        paths.append(path)
    return paths


def refresh_inspection(fixture: dict[str, object]) -> None:
    """Refresh formal product records after an intentional fixture mutation."""

    inspection_path = Path(fixture["inspection"])
    inspection = json.loads(inspection_path.read_text())
    inspection["mhd_history"] = file_record(Path(fixture["mhd"]))
    inspection["user_history"] = file_record(Path(fixture["user"]))
    snapshots = sorted(fixture["snapshot_groups"], key=lambda item: str(item["paths"][0]))
    restarts = sorted(fixture["restart_groups"], key=lambda item: str(item["paths"][0]))
    inspection["snapshots"] = [
        product_record(item["paths"], str(fixture["layout"])) for item in snapshots
    ]
    inspection["snapshot_times"] = [item["time"] for item in snapshots]
    inspection["restarts"] = [
        product_record(item["paths"], str(fixture["layout"])) for item in restarts
    ]
    inspection["restart_times"] = [item["time"] for item in restarts]
    inspection["restart_time_marker_modes"] = [
        [
            effective_marker_mode(item["time"], item["marker_mode"])
            for _ in item["paths"]
        ]
        for item in restarts
    ]
    terminal = next(item for item in restarts if item["time"] == fixture["final_time"])
    inspection["terminal_restart"] = product_record(
        terminal["paths"], str(fixture["layout"])
    )
    inspection["final_hardwall_projection_count"] = fixture["hardwall_history_end"]
    write_json(inspection_path, inspection)
    if fixture["state"] == "recorded":
        manifest_path = Path(fixture["manifest"])
        manifest = json.loads(manifest_path.read_text())
        manifest["scientific_inspection"] = inspection
        write_json(manifest_path, manifest)


def mutate_terminal_restart_state(
    fixture: dict[str, object],
    *,
    injected_work: float | None = None,
    diagnostics: dict[int, float] | None = None,
) -> None:
    """Mutate terminal restart work/diagnostics while preserving sibling semantics."""

    terminal = next(
        item
        for item in fixture["restart_groups"]
        if item["time"] == fixture["final_time"]
    )
    divisor = (
        int(fixture["rank_count"]) if str(fixture["layout"]) == "per_rank" else 1
    )
    for path in terminal["paths"]:
        path = Path(path)
        payload = bytearray(path.read_bytes())
        parameter_end = payload.index(b"<par_end>\n") + len(b"<par_end>\n")
        turbulence_metadata_offset = (
            parameter_end
            + int(TEST_ABI["mesh_header_size"])
            + len(logical_locations()) * 16
            + len(logical_locations()) * 4
        )
        mode_count = struct.unpack_from(
            "<i", payload, turbulence_metadata_offset + 4
        )[0]
        injected_work_offset = (
            turbulence_metadata_offset
            + 248
            + 296
            + 6 * mode_count * 8
        )
        if injected_work is not None:
            struct.pack_into("<d", payload, injected_work_offset, injected_work)
        for index, value in (diagnostics or {}).items():
            struct.pack_into(
                "<d",
                payload,
                injected_work_offset + 8 + index * 8,
                value / divisor,
            )
        path.write_bytes(payload)
    refresh_inspection(fixture)


def mutate_restart_header(
    fixture: dict[str, object],
    relative_offset: int,
    encoded: bytes,
    *,
    every_terminal_sibling: bool = False,
) -> None:
    """Mutate one qualified restart-header field and refresh formal records."""

    terminal = next(
        item
        for item in fixture["restart_groups"]
        if item["time"] == fixture["final_time"]
    )
    paths = terminal["paths"] if every_terminal_sibling else terminal["paths"][-1:]
    for path in paths:
        path = Path(path)
        payload = bytearray(path.read_bytes())
        parameter_end = payload.index(b"<par_end>\n") + len(b"<par_end>\n")
        start = parameter_end + relative_offset
        payload[start:start + len(encoded)] = encoded
        path.write_bytes(payload)
    refresh_inspection(fixture)


def swap_first_two_snapshot_blocks(path: Path) -> None:
    """Apply a valid logical-MeshBlock permutation to one snapshot."""

    payload = bytearray(path.read_bytes())
    block_start = payload.index(b"<par_end>\n") + len(b"<par_end>\n")
    block_size = (
        40
        + 48
        + len(VALIDATOR.EXPECTED_SNAPSHOT_VARIABLES) * 4
    )
    first = bytes(payload[block_start:block_start + block_size])
    second = bytes(payload[block_start + block_size:block_start + 2 * block_size])
    assert len(first) == block_size and len(second) == block_size
    payload[block_start:block_start + 2 * block_size] = second + first
    path.write_bytes(payload)


def mutate_terminal_replicated_restart_state(
    fixture: dict[str, object],
    field: str,
    *,
    every_sibling: bool = False,
    group_index: int = -1,
) -> None:
    """Mutate one valid replicated restart inventory and refresh formal records."""

    group = fixture["restart_groups"][group_index]
    paths = group["paths"] if every_sibling else group["paths"][-1:]
    for source in paths:
        path = Path(source)
        payload = bytearray(path.read_bytes())
        parameter_end = payload.index(b"<par_end>\n") + len(b"<par_end>\n")
        locations = parameter_end + int(TEST_ABI["mesh_header_size"])
        costs = locations + len(logical_locations()) * 16
        metadata = costs + len(logical_locations()) * 4
        rng = metadata + int(TEST_ABI["turbulence_metadata_size"])
        amplitudes = rng + int(TEST_ABI["rng_state_size"])
        if field == "logical_location_inventory":
            first_location = bytes(payload[locations:locations + 16])
            second_location = bytes(payload[locations + 16:locations + 32])
            payload[locations:locations + 32] = second_location + first_location
            first_cost = bytes(payload[costs:costs + 4])
            second_cost = bytes(payload[costs + 4:costs + 8])
            payload[costs:costs + 8] = second_cost + first_cost
        elif field == "cost_inventory":
            struct.pack_into("<f", payload, costs, 9.0)
        elif field == "cost_inventory_permutation":
            first_cost = bytes(payload[costs:costs + 4])
            second_cost = bytes(payload[costs + 4:costs + 8])
            payload[costs:costs + 8] = second_cost + first_cost
        elif field == "turbulence_n_updates":
            struct.pack_into("<i", payload, metadata + 8, fixture_cycle(float(group["time"])) + 1)
        elif field == "turbulence_nlow":
            struct.pack_into("<i", payload, metadata + 12, 99)
        elif field == "turbulence_projection_policy":
            struct.pack_into("<i", payload, metadata + 80, 0)
        elif field == "turbulence_tcorr":
            struct.pack_into("<d", payload, metadata + 96, 9.0)
        elif field == "rng_idum":
            payload[rng] ^= 1
        elif field == "rng_positive_idum":
            struct.pack_into("<q", payload, rng, 101)
        elif field == "rng_negative_idum":
            struct.pack_into("<q", payload, rng, -1)
        elif field == "rng_zero_idum":
            struct.pack_into("<q", payload, rng, 0)
        elif field == "rng_idum2":
            payload[rng + 8] ^= 1
        elif field == "rng_iy":
            payload[rng + 16] ^= 1
        elif field == "rng_iv":
            payload[rng + 24] ^= 1
        elif field == "rng_iset":
            struct.pack_into("<i", payload, rng + 280, 1)
            struct.pack_into("<d", payload, rng + 288, 0.5)
        elif field == "rng_padding":
            payload[rng + 284] ^= 1
        elif field == "rng_gset":
            struct.pack_into("<d", payload, rng + 288, 1.0)
        elif field == "rng_inactive_nonfinite_gset":
            struct.pack_into("<Q", payload, rng + 288, 0x7FF8000000000001)
        elif field == "rng_active_nonfinite_gset":
            struct.pack_into("<i", payload, rng + 280, 1)
            struct.pack_into("<d", payload, rng + 288, float("nan"))
        elif field == "turbulence_amplitudes":
            struct.pack_into("<d", payload, amplitudes, 1.0)
        else:
            raise AssertionError(f"unknown replicated restart field: {field}")
        path.write_bytes(payload)
    refresh_inspection(fixture)


@pytest.fixture
def segment_factory(tmp_path):
    """Create production-shaped segment fixtures with real payload inventories."""

    created = 0

    def make(
        *,
        case_id: str = "R04",
        final_time: float = FINAL_TIME,
        required_time: float | None = None,
        layout: str = "per_rank",
        nodes: int = 1,
        ranks_per_node: int = 2,
        hardwall: bool | None = None,
        state: str = "submitted",
        job_id: str = "12345",
        segment: str = "s00",
    ) -> dict[str, object]:
        nonlocal created
        created += 1
        if hardwall is None:
            hardwall = case_id not in VALIDATOR.FINITE_LIMITER_CASE_IDS
        required = final_time if required_time is None else required_time
        root = tmp_path / f"root-{created}"
        run_dir = (
            root
            / "runs"
            / "mks24-stage-i"
            / EXECUTION_EPOCH
            / case_id
            / segment
        )
        output = run_dir / "output"
        manifest_dir = run_dir / "manifest"
        bin_dir = output / "bin"
        restart_dir = output / "rst"
        manifest_dir.mkdir(parents=True)
        bin_dir.mkdir(parents=True)
        restart_dir.mkdir()
        rank_count = nodes * ranks_per_node
        parameters = input_text(case_id, layout, hardwall)
        archived_input = manifest_dir / "submitted_input.athinput"
        archived_input.write_text(parameters)
        matrix = {
            "campaign": "fixture",
            "cases": [
                {
                    "id": mapped_id,
                    "name": CASE_NAMES[mapped_id],
                    "input": VALIDATOR.CASE_INPUT_PATHS[mapped_id],
                    "resolution": "4x1x1",
                    "figure_roles": ["fixture"],
                }
                for mapped_id in sorted(CASE_NAMES)
            ],
        }
        matrix_path = manifest_dir / "mks24_stage_i_manifest.json"
        write_json(matrix_path, matrix)
        executable = root / "build" / "fixture" / "src" / "athena"
        executable.parent.mkdir(parents=True)
        executable.write_bytes(TEST_EXECUTABLE)
        executable.chmod(0o755)
        utility = root / "source" / "scripts" / "frontier" / "cgl_lf_stage_i.py"
        utility.parent.mkdir(parents=True)
        utility.write_text("fixture utility\n")
        bundle = root / "source-archives" / "athenak-fixture.bundle"
        bundle.parent.mkdir()
        bundle.write_bytes(b"fixture source bundle\n")
        approval = {
            "schema_version": 1,
            "execution_epoch": EXECUTION_EPOCH,
            "approved_executable_revision": EXECUTABLE_REVISION,
            "approved_executable_sha256": TEST_EXECUTABLE_SHA256,
            "approved_utc": "2026-06-05T00:00:00+00:00",
            "approved_by": "independent validator fixture",
            "review_notes": "fixture-only historical executable approval",
        }
        approval_path = (
            root
            / "accounting"
            / "mks24_stage_i_E03_forcing_policy_qualification_approval.json"
        )
        write_json(approval_path, approval)
        batch_path = manifest_dir / "cgl_lf_stage_i.sbatch"

        basename = f"{VALIDATOR.EXECUTION_EPOCH_SLUG}_{CASE_NAMES[case_id]}_{segment}"
        passive = case_id in VALIDATOR.PASSIVE_CASE_IDS
        snapshot_groups = []
        restart_groups = []
        snapshot_times = [0.0]
        while snapshot_times[-1] + 0.25 < final_time:
            snapshot_times.append(snapshot_times[-1] + 0.25)
        if snapshot_times[-1] != final_time:
            snapshot_times.append(final_time)
        restart_times = [0.0]
        while restart_times[-1] + 1.0 < final_time:
            restart_times.append(restart_times[-1] + 1.0)
        if restart_times[-1] != final_time:
            restart_times.append(final_time)
        for index, time in enumerate(snapshot_times):
            snapshot_paths = write_product_group(
                directory=bin_dir,
                basename=basename,
                index=index,
                suffix="bin",
                layout=layout,
                rank_count=rank_count,
                writer=lambda path, locations, _count, time=time: write_snapshot(
                    path, time, fixture_cycle(time), locations, parameters, basename, required
                ),
            )
            snapshot_groups.append({"paths": snapshot_paths, "time": time})
        for index, time in enumerate(restart_times):
            marker_mode = "full_precision" if index == 0 else "legacy_default_precision"
            restart_paths = write_product_group(
                directory=restart_dir,
                basename=basename,
                index=index,
                suffix="rst",
                layout=layout,
                rank_count=rank_count,
                writer=lambda path, _locations, count, time=time, marker_mode=marker_mode: (
                    write_restart(
                        path,
                        time,
                        fixture_cycle(time),
                        marker_mode,
                        count,
                        parameters,
                        basename,
                        required,
                        passive=passive,
                        hardwall=hardwall,
                        diagnostic_divisor=(
                            rank_count if layout == "per_rank" else 1
                        ),
                    )
                ),
            )
            restart_groups.append({
                "paths": restart_paths,
                "time": time,
                "marker_mode": marker_mode,
            })

        mhd = output / f"{basename}.mhd.hst"
        user = output / f"{basename}.user.hst"
        fixture: dict[str, object] = {
            "root": root,
            "run_dir": run_dir,
            "output": output,
            "manifest": manifest_dir / "prepared_run.json",
            "inspection": manifest_dir / "segment_inspection.json",
            "mhd": mhd,
            "user": user,
            "matrix": matrix_path,
            "bundle": bundle,
            "case_id": case_id,
            "segment": segment,
            "basename": basename,
            "job_id": job_id,
            "layout": layout,
            "rank_count": rank_count,
            "parameters": parameters,
            "passive": passive,
            "hardwall": hardwall,
            "state": state,
            "final_time": final_time,
            "required_time": required,
            "snapshot_groups": snapshot_groups,
            "restart_groups": restart_groups,
        }
        write_histories(fixture)
        manifest_path = Path(fixture["manifest"])
        manifest = {
            "schema_version": 3,
            "execution_epoch": EXECUTION_EPOCH,
            "state": state,
            "job_id": job_id,
            "project_root": str(root),
            "run": {
                "case_id": case_id,
                "case_name": CASE_NAMES[case_id],
                "segment": segment,
                "resolution": "4x1x1",
                "figure_roles": ["fixture"],
            },
            "allocation": {
                "nodes": nodes,
                "ranks_per_node": ranks_per_node,
            },
            "command": {
                "production_utility": {
                    "committed": True,
                    "path": str(utility),
                    "revision": UTILITY_REVISION,
                    "sha256": sha256(utility),
                },
                "qualification_approval": {
                    "path": str(approval_path),
                    "sha256": sha256(approval_path),
                    "execution_epoch": EXECUTION_EPOCH,
                    "approved_executable_revision": EXECUTABLE_REVISION,
                    "approved_executable_sha256": TEST_EXECUTABLE_SHA256,
                    "token": approval,
                },
                "source_dir": str(root / "source"),
                "source_bundle": {
                    "path": str(bundle),
                    "sha256": sha256(bundle),
                    "verified_revisions": [EXECUTABLE_REVISION, UTILITY_REVISION],
                },
                "input_revision": EXECUTABLE_REVISION,
                "source_input_file": str(root / "source" / VALIDATOR.CASE_INPUT_PATHS[case_id]),
                "input_file": str(archived_input),
                "input_sha256": sha256(archived_input),
                "matrix_file": str(matrix_path),
                "matrix_sha256": sha256(matrix_path),
                "executable": str(executable),
                "executable_revision": EXECUTABLE_REVISION,
                "executable_sha256": TEST_EXECUTABLE_SHA256,
                "restart_file": None,
                "restart_sha256": None,
                "restart_files": [],
                "source_restart_file": None,
                "parent_segment": None,
                "allow_missing_restart_time_marker": False,
                "overrides": [f"time/tlim={format(required, '.17g')}"],
                "time_tlim_target": required,
                "athena_walltime": "00:10:00",
            },
            "paths": {
                "run_dir": str(run_dir),
                "output_dir": str(output),
                "batch_script": str(batch_path),
                "environment_log": str(manifest_dir / "run_environment.txt"),
                "slurm_log": str(root / "logs" / "slurm" / "%x.%j.log"),
            },
        }
        manifest["run"]["run_basename"] = basename
        manifest["allocation"].update({
            "requested_walltime": "00:20:00",
            "requested_seconds": 1200,
            "reserved_node_hours": nodes / 3.0,
            "cpus_per_task": 7,
        })
        batch_text, batch_digest = finalized_batch_script(manifest, manifest_path)
        batch_path.write_text(batch_text)
        manifest["command"]["batch_script_sha256"] = batch_digest
        accepted = final_time >= required - 1.0e-10
        terminal_restart = restart_groups[-1]
        inspection = {
            "schema_version": 4,
            "execution_epoch": EXECUTION_EPOCH,
            "manifest": str(manifest_path),
            "job_id": job_id,
            "case_id": case_id,
            "segment": segment,
            "required_time": required,
            "final_time": final_time,
            "maximum_strict_failure_counts": {
                name: 0.0 for name in VALIDATOR.STRICT_LF_FAILURE_COLUMNS
            },
            "checks": {
                "required_time_reached": accepted,
                "strict_lf_failure_counters_zero": True,
                "snapshots_retained": True,
                "terminal_snapshot_retained": True,
                "restart_retained": True,
                "terminal_restart_physical_time_matches_final": True,
            },
            "accepted": accepted,
            "clean_for_continuation": True,
            "mhd_history": file_record(mhd),
            "user_history": file_record(user),
            "snapshots": [
                product_record(item["paths"], layout) for item in snapshot_groups
            ],
            "snapshot_times": [item["time"] for item in snapshot_groups],
            "restarts": [
                product_record(item["paths"], layout) for item in restart_groups
            ],
            "restart_times": [item["time"] for item in restart_groups],
            "restart_time_marker_modes": [
                [
                    effective_marker_mode(item["time"], item["marker_mode"])
                    for _ in item["paths"]
                ]
                for item in restart_groups
            ],
            "terminal_restart": product_record(terminal_restart["paths"], layout),
            "terminal_restart_time": final_time,
            "restart_time_marker_bypass": False,
            "final_hardwall_projection_count": 2.0 if hardwall else 0.0,
        }
        write_json(Path(fixture["inspection"]), inspection)
        if state == "recorded":
            manifest["accounting"] = {
                "execution_epoch": EXECUTION_EPOCH,
                "job_id": job_id,
                "case_id": case_id,
                "segment": segment,
                "result": "accepted" if accepted else "clean_partial",
            }
            manifest["scientific_inspection"] = inspection
        write_json(manifest_path, manifest)
        fixture["utility"] = utility
        fixture["bundle_blobs"] = {
            "input": archived_input.read_bytes(),
            "matrix": matrix_path.read_bytes(),
            "production utility": utility.read_bytes(),
        }
        return fixture

    return make


def validate(
    fixture: dict[str, object],
    *,
    policy: str = "stage-i-standard-v1",
    result: str = "accepted",
) -> dict[str, object]:
    """Run the validator directly against one local fixture."""

    canonical_root = VALIDATOR.CANONICAL_PROJECT_ROOT
    bundle_reader = VALIDATOR.authenticated_bundle_blobs
    VALIDATOR.CANONICAL_PROJECT_ROOT = Path(fixture["root"])
    VALIDATOR.authenticated_bundle_blobs = lambda _path, _requests, _sha=None: dict(
        fixture["bundle_blobs"]
    )
    try:
        return VALIDATOR.validate_segment(
            argparse.Namespace(
                manifest=Path(fixture["manifest"]),
                inspection=None,
                policy=policy,
                result=result,
            )
        )
    finally:
        VALIDATOR.CANONICAL_PROJECT_ROOT = canonical_root
        VALIDATOR.authenticated_bundle_blobs = bundle_reader


def require_rejected(
    fixture: dict[str, object],
    message: str,
    *,
    policy: str = "stage-i-standard-v1",
    result: str = "accepted",
) -> None:
    """Require one fail-closed validator rejection."""

    with pytest.raises(VALIDATOR.ValidationError, match=re.escape(message)):
        validate(fixture, policy=policy, result=result)


def add_future_snapshot(fixture: dict[str, object], time: float) -> None:
    """Add a formally recorded but scientifically invalid future snapshot."""

    index = len(fixture["snapshot_groups"])
    paths = write_product_group(
        directory=Path(fixture["output"]) / "bin",
        basename=str(fixture["basename"]),
        index=index,
        suffix="bin",
        layout=str(fixture["layout"]),
        rank_count=int(fixture["rank_count"]),
        writer=lambda path, locations, _count: write_snapshot(
            path,
            time,
            fixture_cycle(time),
            locations,
            str(fixture["parameters"]),
            str(fixture["basename"]),
            float(fixture["required_time"]),
        ),
    )
    fixture["snapshot_groups"].append({"paths": paths, "time": time})
    refresh_inspection(fixture)


def add_restart(fixture: dict[str, object], time: float) -> None:
    """Add one formally recorded restart product at an arbitrary physical time."""

    index = len(fixture["restart_groups"])
    marker_mode = "legacy_default_precision"
    paths = write_product_group(
        directory=Path(fixture["output"]) / "rst",
        basename=str(fixture["basename"]),
        index=index,
        suffix="rst",
        layout=str(fixture["layout"]),
        rank_count=int(fixture["rank_count"]),
        writer=lambda path, _locations, count: write_restart(
            path,
            time,
            fixture_cycle(time),
            marker_mode,
            count,
            str(fixture["parameters"]),
            str(fixture["basename"]),
            float(fixture["required_time"]),
            passive=bool(fixture["passive"]),
            hardwall=bool(fixture["hardwall"]),
            diagnostic_divisor=(
                int(fixture["rank_count"])
                if str(fixture["layout"]) == "per_rank"
                else 1
            ),
        ),
    )
    fixture["restart_groups"].append({
        "paths": paths,
        "time": time,
        "marker_mode": marker_mode,
    })
    refresh_inspection(fixture)


def configure_continuation(fixture: dict[str, object], start: float = 0.1) -> None:
    """Configure one production-shaped rank-local continuation launch."""

    manifest_path = Path(fixture["manifest"])
    manifest = json.loads(manifest_path.read_text())
    root = Path(fixture["root"])
    parent_segment_name = "parent_s00"
    parent_run = (
        root
        / "runs"
        / "mks24-stage-i"
        / EXECUTION_EPOCH
        / str(fixture["case_id"])
        / parent_segment_name
    )
    parent_manifest_path = parent_run / "manifest" / "prepared_run.json"
    parent_output = parent_run / "output" / "rst"
    parent_basename = (
        f"{VALIDATOR.EXECUTION_EPOCH_SLUG}_{CASE_NAMES[str(fixture['case_id'])]}_"
        f"{parent_segment_name}"
    )
    parent_paths = write_product_group(
        directory=parent_output,
        basename=parent_basename,
        index=0,
        suffix="rst",
        layout=str(fixture["layout"]),
        rank_count=int(fixture["rank_count"]),
        writer=lambda path, _locations, count: write_restart(
            path,
            start,
            fixture_cycle(start),
            "full_precision",
            count,
            str(fixture["parameters"]),
            parent_basename,
            start,
            passive=bool(fixture["passive"]),
            hardwall=bool(fixture["hardwall"]),
            diagnostic_divisor=(
                int(fixture["rank_count"])
                if str(fixture["layout"]) == "per_rank"
                else 1
            ),
        ),
    )
    parent_terminal = product_record(parent_paths, str(fixture["layout"]))
    parent_inspection = {
        "schema_version": 4,
        "execution_epoch": EXECUTION_EPOCH,
        "manifest": str(parent_manifest_path),
        "job_id": "54321",
        "case_id": fixture["case_id"],
        "segment": parent_segment_name,
        "final_time": start,
        "terminal_restart_time": start,
        "clean_for_continuation": True,
        "terminal_restart": parent_terminal,
    }
    parent_manifest = {
        "schema_version": 3,
        "execution_epoch": EXECUTION_EPOCH,
        "state": "recorded",
        "job_id": "54321",
        "project_root": str(root),
        "run": {
            "case_id": fixture["case_id"],
            "segment": parent_segment_name,
            "run_basename": parent_basename,
        },
        "command": {
            "input_sha256": manifest["command"]["input_sha256"],
            "executable_sha256": manifest["command"]["executable_sha256"],
            "time_tlim_target": start,
        },
        "accounting": {
            "execution_epoch": EXECUTION_EPOCH,
            "job_id": "54321",
            "case_id": fixture["case_id"],
            "segment": parent_segment_name,
            "result": "accepted",
        },
        "scientific_inspection": parent_inspection,
    }
    write_json(parent_manifest_path.parent / "segment_inspection.json", parent_inspection)
    write_json(parent_manifest_path, parent_manifest)

    archive = manifest_path.parent / "submitted_restart"
    records = []
    for source in parent_paths:
        relative = (
            Path(source.parent.name) / source.name
            if str(fixture["layout"]) == "per_rank"
            else Path(source.name)
        )
        archived = archive / relative
        archived.parent.mkdir(parents=True, exist_ok=True)
        archived.write_bytes(source.read_bytes())
        records.append(file_record(archived))
    manifest["command"].update({
        "restart_file": records[0]["path"],
        "restart_sha256": records[0]["sha256"],
        "restart_files": records,
        "source_restart_file": str(parent_paths[0]),
        "parent_segment": {
            "execution_epoch": EXECUTION_EPOCH,
            "manifest": str(parent_manifest_path),
            "case_id": fixture["case_id"],
            "segment": parent_segment_name,
            "result": "accepted",
            "restart_sha256": parent_terminal["sha256"],
            "restart_files": [str(path) for path in parent_paths],
            "final_time": start,
            "restart_time": start,
            "input_sha256": manifest["command"]["input_sha256"],
            "executable_sha256": manifest["command"]["executable_sha256"],
        },
    })
    write_json(manifest_path, manifest)
    refresh_batch_script(fixture)


def remove_history_data_row(path: Path, row_index: int) -> None:
    """Remove one indexed non-comment history row."""

    lines = path.read_text().splitlines()
    data_indices = [
        index for index, line in enumerate(lines) if line and not line.startswith("#")
    ]
    del lines[data_indices[row_index]]
    path.write_text("\n".join(lines) + "\n")


def trim_history_before(path: Path, start: float) -> None:
    """Retain the exact history header and rows at or after one boundary."""

    lines = path.read_text().splitlines()
    retained = [
        line
        for line in lines
        if line.startswith("#")
        or not line
        or float(line.split()[0]) >= start
    ]
    path.write_text("\n".join(retained) + "\n")


def test_validator_accepts_real_payload_rank_local_inventory(segment_factory):
    fixture = segment_factory()
    result = validate(fixture)
    assert result["validation_accepted"]
    assert result["case_mode"] == "active"
    assert result["expected_ranks"] == 2
    assert result["expected_meshblocks"] == 4
    assert result["product_count"] == 12
    assert result["sampled_history_forcing_work"]["claim"] == "active_delta_energy_closure"
    assert result["sampled_history_forcing_work"]["normalized_residual"] < 1.0e-8
    assert result["sampled_history_forcing_work"]["thresholds"] == {
        "normalized_residual_lt": 1.0e-8
    }
    assert result["named_activity_validation"]["forcing_work"]["absolute_gt"] == 1.0e-6
    assert (
        result["named_activity_validation"]["forcing_work"]["state_normalized_gt"]
        == 1.0e-8
    )
    assert result["lf_work_diagnostics_finite"]
    assert result["schema_version"] == 3
    assert result["record_role"] == "historical_read_only_non_authorizing"
    assert result["authorization_effect"] == "none"
    assert not result["independent_authorization"]["authorizing"]
    assert result["independent_authorization"]["status"] == "non_authorizing"
    assert not result["normalized_ct_divb_validation"]["authorizing"]
    assert not result["restart_face_field_validation"]["authorizing"]
    assert "max_normalized_divb" not in result
    assert result["historical_scope_validation"]["status"] == (
        "frozen_historical_executable_validated"
    )
    assert not result["historical_scope_validation"]["authorizing"]
    assert result["validator_revision"] == FAKE_VALIDATOR_REVISION
    assert result["validator_implementation"]["committed"]
    assert not result["validator_implementation"]["authorizing"]
    assert result["snapshot_ordered_location_inventories"]
    assert result["restart_ordered_location_inventories"]
    assert result["restart_ordered_cost_inventories"][0] == (1.0, 2.0, 3.0, 4.0)
    assert result["restart_replicated_state_sha256"]
    assert result["parent_acceptance_validation"]["status"] == (
        "not_applicable_fresh_segment"
    )
    assert not result["parent_acceptance_validation"]["authorizing"]


@pytest.mark.parametrize("case_id", sorted(VALIDATOR.APPROVED_CASE_IDS))
def test_validator_accepts_each_exact_frozen_r02_r17_contract(
    segment_factory, case_id
):
    result = validate(segment_factory(case_id=case_id))
    assert result["case_id"] == case_id
    assert result["case_mode"] == (
        "passive" if case_id in VALIDATOR.PASSIVE_CASE_IDS else "active"
    )


def test_validator_accepts_exact_frozen_source_r02_r17_contracts(monkeypatch):
    use_production_contracts(monkeypatch)
    matrix_payload = frozen_source(VALIDATOR.MATRIX_REPOSITORY_PATH)
    assert hashlib.sha256(matrix_payload).hexdigest() == FROZEN_MATRIX_SHA256
    matrix = json.loads(matrix_payload)
    cases = {item["id"]: item for item in matrix["cases"]}
    assert tuple(cases) == tuple(f"R{number:02d}" for number in range(2, 18))
    assert set(FROZEN_INPUT_SOURCES) == set(cases)
    for case_id, (path, digest) in FROZEN_INPUT_SOURCES.items():
        assert cases[case_id]["input"] == path
        payload = frozen_source(path)
        assert hashlib.sha256(payload).hexdigest() == digest
        contract = VALIDATOR.require_input_contract(
            VALIDATOR.parse_athinput(payload, f"frozen source {case_id}"),
            case_id,
            cases[case_id],
        )
        assert contract["expected_meshblocks"] > 0


def test_validator_git_run_ignores_hostile_path(tmp_path, monkeypatch):
    hostile = tmp_path / "hostile"
    hostile.mkdir()
    fake_git = hostile / "git"
    fake_git.write_text("#!/bin/sh\necho hostile-path-git\nexit 73\n")
    fake_git.chmod(0o755)
    monkeypatch.setenv("PATH", str(hostile))
    result = VALIDATOR.git_run(["--version"])
    assert result.returncode == 0
    assert result.stdout.startswith(b"git version ")
    assert b"hostile-path-git" not in result.stdout


def test_validator_git_run_uses_authenticated_descriptor_and_isolated_config(
    monkeypatch,
):
    captured = {}

    def capture_run(arguments, **kwargs):
        captured["arguments"] = arguments
        captured.update(kwargs)
        return subprocess.CompletedProcess(arguments, 0, stdout=b"", stderr=b"")

    monkeypatch.setenv("GIT_CONFIG_GLOBAL", "/hostile/global")
    monkeypatch.setenv("GIT_CONFIG_SYSTEM", "/hostile/system")
    monkeypatch.setattr(VALIDATOR.subprocess, "run", capture_run)
    VALIDATOR.git_run(["status"])
    assert captured["arguments"][0] == str(VALIDATOR.GIT_EXECUTABLE)
    assert captured["executable"].startswith("/proc/self/fd/")
    assert captured["pass_fds"]
    assert captured["env"]["GIT_CONFIG_NOSYSTEM"] == "1"
    assert captured["env"]["GIT_CONFIG_GLOBAL"] == os.devnull
    assert captured["env"]["GIT_CONFIG_SYSTEM"] == os.devnull
    assert captured["env"]["GIT_TERMINAL_PROMPT"] == "0"


@pytest.mark.parametrize(
    ("payload", "message"),
    [
        ('{"outer": {"schema_version": 1, "schema_version": 1}}', "repeats JSON key schema_version"),
        ('{"value": NaN}', "contains non-finite JSON value NaN"),
        ('{"value": Infinity}', "contains non-finite JSON value Infinity"),
        ('{"value": -Infinity}', "contains non-finite JSON value -Infinity"),
        ('{"value": 1e999}', "contains a non-finite JSON number"),
    ],
)
def test_authority_json_rejects_duplicate_keys_and_nonfinite_numbers(
    tmp_path, payload, message
):
    path = tmp_path / "authority.json"
    path.write_text(payload)
    with pytest.raises(VALIDATOR.ValidationError, match=re.escape(message)):
        VALIDATOR.load_json(path, "authority record")


@pytest.mark.parametrize(
    ("authority", "key", "message"),
    [
        ("manifest", "schema_version", "prepared manifest repeats JSON key schema_version"),
        ("inspection", "schema_version", "segment inspection repeats JSON key schema_version"),
        ("matrix", "campaign", "archived matrix repeats JSON key campaign"),
        (
            "qualification",
            "schema_version",
            "qualification approval repeats JSON key schema_version",
        ),
    ],
)
def test_validator_strictly_parses_each_authority_json(
    segment_factory, authority, key, message
):
    fixture = segment_factory()
    if authority == "manifest":
        path = Path(fixture["manifest"])
        duplicate_json_key(path, key)
    elif authority == "inspection":
        path = Path(fixture["inspection"])
        duplicate_json_key(path, key)
    elif authority == "matrix":
        path = Path(fixture["matrix"])
        duplicate_json_key(path, key)
        bind_fixture_artifact(fixture, "matrix", path, "matrix_sha256")
    else:
        manifest_path = Path(fixture["manifest"])
        manifest = json.loads(manifest_path.read_text())
        path = Path(manifest["command"]["qualification_approval"]["path"])
        duplicate_json_key(path, key)
        manifest["command"]["qualification_approval"]["sha256"] = sha256(path)
        write_json(manifest_path, manifest)
        refresh_batch_script(fixture)
    require_rejected(fixture, message)


def test_validator_rejects_executable_without_historical_abi(
    segment_factory, monkeypatch
):
    fixture = segment_factory()
    monkeypatch.delitem(
        VALIDATOR.QUALIFIED_RESTART_BINARY_ABIS,
        (EXECUTABLE_REVISION, TEST_EXECUTABLE_SHA256),
    )
    require_rejected(
        fixture, "prepared executable has no qualified binary-restart ABI"
    )


def test_validator_discloses_every_remaining_non_authorizing_limitation(
    segment_factory,
):
    result = validate(segment_factory())
    reasons = set(result["independent_authorization"]["reasons"])
    assert reasons == {
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
    }
    assert not result["restart_load_validation"]["authorizing"]
    assert not result["restart_stochastic_state_validation"]["authorizing"]
    assert result["restart_stochastic_state_validation"]["status"] == (
        "lifecycle_valid_canonical_continuation_state_consistent_ignored_bytes_recorded_not_replayed"
    )
    assert set(
        result["restart_stochastic_state_validation"]["authenticated_field_counts"]
    ) == {1, 36}
    assert result["restart_stochastic_state_validation"]["native_padding_offset"] == 284
    assert result["restart_stochastic_state_validation"]["native_padding_size_bytes"] == 4
    assert result["restart_rng_native_padding_evidence"]
    assert result["restart_rng_gset_evidence"]
    assert result["restart_rng_pre_update_ignored_state_evidence"][0]
    assert result["restart_rng_lifecycle_evidence"][0]["lifecycle"] == "pre_update_seed"
    assert result["restart_rng_lifecycle_evidence"][-1]["lifecycle"] == (
        "initialized_continuation"
    )
    assert not result["rank_local_restart_mapping_validation"]["authorizing"]
    assert not result["operational_authorization_validation"]["authorizing"]
    assert not result["historical_runtime_checksum_validation"]["authorizing"]
    assert not result["scientific_interpretation_validation"]["authorizing"]
    assert result["continuation_origin_validation"] is None
    assert result["filesystem_race_validation"]["status"] == (
        "descriptor_profile_bound_with_final_recheck"
    )
    assert result["parameter_contract_validation"]["status"] == (
        "closed_qualified_inventory"
    )
    assert result["cycle_dt_consistency_validation"]["status"] == "validated"
    assert result["cycle_dt_consistency_validation"]["snapshot_cycles"] == (
        result["snapshot_cycles"]
    )
    assert result["cycle_dt_consistency_validation"]["restart_cycles"] == (
        result["restart_cycles"]
    )
    assert result["cycle_dt_consistency_validation"]["restart_dts"] == (
        result["restart_dts"]
    )
    assert "closed qualified parameter contract" in (
        result["product_payload_validation"]["snapshots"]
    )
    assert "closed qualified parameter contract" in (
        result["product_payload_validation"]["restarts"]
    )


def test_validator_rejects_noncanonical_project_root(segment_factory):
    fixture = segment_factory()
    with pytest.raises(
        VALIDATOR.ValidationError,
        match="manifest project root is not the canonical Stage I root",
    ):
        VALIDATOR.validate_segment(
            argparse.Namespace(
                manifest=Path(fixture["manifest"]),
                inspection=None,
                policy="stage-i-standard-v1",
                result="accepted",
            )
        )


def test_validator_accepts_multi_node_shared_mpiio_layout(segment_factory, monkeypatch):
    monkeypatch.setitem(
        VALIDATOR.COMMON_FROZEN_PARAMETERS,
        ("output2", "single_file_per_rank"),
        "false",
    )
    monkeypatch.setitem(
        VALIDATOR.COMMON_FROZEN_PARAMETERS,
        ("output3", "single_file_per_rank"),
        "false",
    )
    fixture = segment_factory(layout="shared_mpiio", nodes=2, ranks_per_node=2)
    result = validate(fixture)
    assert result["expected_ranks"] == 4
    assert result["product_count"] == 7
    assert all(
        record == "shared_mpiio"
        for record in (
            json.loads(Path(fixture["inspection"]).read_text())["snapshots"][0]["storage"],
            json.loads(Path(fixture["inspection"]).read_text())["restarts"][0]["storage"],
        )
    )


@pytest.mark.parametrize("case_id", sorted(VALIDATOR.PASSIVE_CASE_IDS))
def test_validator_classifies_passive_and_rejects_pressure_feedback(
    segment_factory, case_id
):
    fixture = segment_factory(case_id=case_id)
    result = validate(fixture)
    assert result["case_mode"] == "passive"
    assert (
        result["sampled_history_forcing_work"]["claim"]
        == "passive_delta_no_energy_closure_claim"
    )
    assert "normalized_residual" not in result["sampled_history_forcing_work"]
    write_histories(fixture, passive_pressure_feedback=True)
    refresh_inspection(fixture)
    require_rejected(fixture, "passive case contains pressure-work feedback")


def test_validator_binds_named_policies_and_rejects_relaxation(segment_factory):
    fixture = segment_factory()
    require_rejected(
        fixture,
        "named threshold policy is not approved for this segment",
        policy="r03-4762472-clean-partial-v1",
    )
    with pytest.raises(SystemExit):
        VALIDATOR.parse_args([
            "--manifest",
            str(fixture["manifest"]),
            "--policy",
            "relaxed-local-policy",
        ])
    with pytest.raises(SystemExit):
        VALIDATOR.parse_args([
            "--manifest",
            str(fixture["manifest"]),
            "--policy",
            "stage-i-standard-v1",
            "--forcing-absolute-residual-lt",
            "1",
        ])


def test_validator_binds_recorded_scientific_inspection(segment_factory):
    fixture = segment_factory(state="recorded")
    assert validate(fixture)["validation_accepted"]
    inspection = json.loads(Path(fixture["inspection"]).read_text())
    inspection["accepted"] = False
    write_json(Path(fixture["inspection"]), inspection)
    require_rejected(
        fixture,
        "recorded manifest scientific_inspection differs from retained inspection",
    )


def test_validator_rejects_parent_symlink_in_provenance(segment_factory):
    fixture = segment_factory()
    bundle = Path(fixture["bundle"])
    real_parent = bundle.parent.with_name("source-archives-real")
    bundle.parent.rename(real_parent)
    bundle.parent.symlink_to(real_parent, target_is_directory=True)
    require_rejected(fixture, "component may be a symlink")


def test_validator_rejects_authenticated_file_mutation_race(
    segment_factory, monkeypatch
):
    fixture = segment_factory()
    original = VALIDATOR.require_unchanged_output_inventory

    def mutate_after_inventory_check(*args, **kwargs):
        original(*args, **kwargs)
        path = Path(fixture["mhd"])
        path.write_text(path.read_text() + "# raced after authenticated read\n")

    monkeypatch.setattr(
        VALIDATOR, "require_unchanged_output_inventory", mutate_after_inventory_check
    )
    require_rejected(fixture, "authenticated retained file changed during validation")


def test_validator_rejects_final_output_inventory_addition_race(
    segment_factory, monkeypatch
):
    fixture = segment_factory()
    original = VALIDATOR.recheck_profiles
    calls = 0

    def add_product_after_first_final_profile_recheck(profiles):
        nonlocal calls
        original(profiles)
        calls += 1
        if calls == 1:
            (Path(fixture["output"]) / "bin" / "raced.bin").write_bytes(b"raced\n")

    monkeypatch.setattr(
        VALIDATOR, "recheck_profiles", add_product_after_first_final_profile_recheck
    )
    require_rejected(fixture, "snapshot inventory changed during validation")


def test_validator_rejects_semantic_matrix_provenance_mismatch(segment_factory):
    fixture = segment_factory()
    matrix_path = Path(fixture["matrix"])
    matrix = json.loads(matrix_path.read_text())
    next(item for item in matrix["cases"] if item["id"] == fixture["case_id"])[
        "name"
    ] = "wrong_case"
    write_json(matrix_path, matrix)
    bind_fixture_artifact(fixture, "matrix", matrix_path, "matrix_sha256")
    require_rejected(fixture, "retained matrix name differs from manifest run")


def test_validator_rejects_source_bundle_revision_mismatch(segment_factory):
    fixture = segment_factory()
    manifest_path = Path(fixture["manifest"])
    manifest = json.loads(manifest_path.read_text())
    manifest["command"]["source_bundle"]["verified_revisions"] = [EXECUTABLE_REVISION]
    write_json(manifest_path, manifest)
    require_rejected(fixture, "source bundle does not bind every launch revision")


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ({"schema_version": 2}, "qualification approval token has wrong schema"),
        ({"review_notes": ""}, "qualification approval token lacks nonempty review_notes"),
    ],
)
def test_validator_rejects_incomplete_qualification_approval(
    segment_factory, mutation, message
):
    fixture = segment_factory()
    manifest_path = Path(fixture["manifest"])
    manifest = json.loads(manifest_path.read_text())
    approval_path = Path(manifest["command"]["qualification_approval"]["path"])
    approval = json.loads(approval_path.read_text())
    approval.update(mutation)
    write_json(approval_path, approval)
    manifest["command"]["qualification_approval"].update({
        "sha256": sha256(approval_path),
        "token": approval,
    })
    write_json(manifest_path, manifest)
    refresh_batch_script(fixture)
    require_rejected(fixture, message)


def test_validator_rejects_archived_input_not_in_bundle(segment_factory):
    fixture = segment_factory()
    input_path = Path(fixture["run_dir"]) / "manifest" / "submitted_input.athinput"
    input_path.write_text(input_path.read_text() + "# fabricated local input\n")
    manifest_path = Path(fixture["manifest"])
    manifest = json.loads(manifest_path.read_text())
    manifest["command"]["input_sha256"] = sha256(input_path)
    write_json(manifest_path, manifest)
    refresh_batch_script(fixture)
    require_rejected(fixture, "archived input differs from authenticated bundle object")


def test_validator_rejects_archived_matrix_not_in_bundle(segment_factory):
    fixture = segment_factory()
    matrix_path = Path(fixture["matrix"])
    matrix = json.loads(matrix_path.read_text())
    matrix["fabricated_local_record"] = True
    write_json(matrix_path, matrix)
    manifest_path = Path(fixture["manifest"])
    manifest = json.loads(manifest_path.read_text())
    manifest["command"]["matrix_sha256"] = sha256(matrix_path)
    write_json(manifest_path, manifest)
    refresh_batch_script(fixture)
    require_rejected(fixture, "archived matrix differs from authenticated bundle object")


def test_validator_rejects_production_utility_not_in_bundle(segment_factory):
    fixture = segment_factory()
    fixture["bundle_blobs"]["production utility"] = b"fabricated utility object\n"
    require_rejected(
        fixture, "production utility digest differs from authenticated bundle object"
    )


def test_validator_rejects_self_digesting_altered_launch_intent(segment_factory):
    fixture = segment_factory()
    batch_path = Path(fixture["run_dir"]) / "manifest" / "cgl_lf_stage_i.sbatch"
    altered = batch_path.read_text().replace("srun -N", "echo srun -N", 1)
    normalized = VALIDATOR.BATCH_SCRIPT_DIGEST_PATTERN.sub(
        f"BATCH_SCRIPT_SHA256={VALIDATOR.BATCH_SCRIPT_DIGEST_PLACEHOLDER}", altered
    )
    digest = hashlib.sha256(normalized.encode()).hexdigest()
    altered = VALIDATOR.BATCH_SCRIPT_DIGEST_PATTERN.sub(
        f"BATCH_SCRIPT_SHA256={digest}", altered
    )
    batch_path.write_text(altered)
    manifest_path = Path(fixture["manifest"])
    manifest = json.loads(manifest_path.read_text())
    manifest["command"]["batch_script_sha256"] = digest
    write_json(manifest_path, manifest)
    require_rejected(fixture, "prepared batch script differs from regenerated launch intent")


@pytest.mark.parametrize(
    ("keys", "value", "message"),
    [
        (
            ("run", "run_basename"),
            "valid\n#SBATCH -q debug",
            "manifest run identity is not exact production launch intent",
        ),
        (
            ("allocation", "nodes"),
            "1\n#SBATCH -q debug",
            "allocation nodes is not an exact integer",
        ),
        (
            ("allocation", "ranks_per_node"),
            True,
            "allocation ranks_per_node is not an exact integer",
        ),
        (
            ("allocation", "cpus_per_task"),
            "7; id",
            "allocation cpus_per_task is not an exact integer",
        ),
        (
            ("allocation", "requested_walltime"),
            "00:20:00\n#SBATCH -q debug",
            "allocation requested_walltime is not an HH:MM:SS string",
        ),
        (
            ("allocation", "requested_walltime"),
            "٠٠:٢٠:٠٠",
            "allocation requested_walltime is not an HH:MM:SS string",
        ),
        (
            ("allocation", "requested_seconds"),
            1200.0,
            "allocation requested_seconds is not an exact integer",
        ),
        (
            ("allocation", "reserved_node_hours"),
            "0.3333333333333333",
            "manifest allocation is not exact production launch intent",
        ),
        (
            ("command", "athena_walltime"),
            "00:10:00; id",
            "command athena_walltime is not an HH:MM:SS string",
        ),
        (
            ("command", "overrides"),
            ["time/tlim=0.31282347945569927; id"],
            "manifest does not bind exactly one time/tlim override",
        ),
        (
            ("paths", "slurm_log"),
            "/tmp/log\n#SBATCH -q debug",
            "launch path slurm_log path contains a control character",
        ),
    ],
)
def test_validator_rejects_typed_launch_directive_and_shell_injection_before_regeneration(
    segment_factory, monkeypatch, keys, value, message
):
    fixture = segment_factory()
    mutate_manifest_value(fixture, keys, value)

    def regeneration_must_not_run(*_args, **_kwargs):
        pytest.fail("invalid launch fields reached script regeneration")

    monkeypatch.setattr(VALIDATOR, "regenerated_batch_script", regeneration_must_not_run)
    require_rejected(fixture, message)


def test_validator_rejects_restart_launch_injection_before_regeneration(
    segment_factory, monkeypatch
):
    fixture = segment_factory()
    mutate_manifest_value(
        fixture,
        ("command", "restart_files"),
        [{
            "path": f"{fixture['run_dir']}/manifest/submitted_restart/x\nid",
            "sha256": "0" * 64,
            "size_bytes": 1,
        }],
    )

    def regeneration_must_not_run(*_args, **_kwargs):
        pytest.fail("invalid restart launch fields reached script regeneration")

    monkeypatch.setattr(VALIDATOR, "regenerated_batch_script", regeneration_must_not_run)
    require_rejected(fixture, "restart launch sibling path contains a control character")


def test_validator_rejects_incomplete_rank_local_restart_launch_before_regeneration(
    segment_factory, monkeypatch
):
    fixture = segment_factory()
    configure_continuation(fixture)
    manifest_path = Path(fixture["manifest"])
    manifest = json.loads(manifest_path.read_text())
    manifest["command"]["restart_files"].pop()
    write_json(manifest_path, manifest)

    def regeneration_must_not_run(*_args, **_kwargs):
        pytest.fail("incomplete restart inventory reached script regeneration")

    monkeypatch.setattr(VALIDATOR, "regenerated_batch_script", regeneration_must_not_run)
    require_rejected(fixture, "manifest restart launch inventory is invalid")


def test_authenticated_bundle_blobs_reads_exact_committed_objects(tmp_path, monkeypatch):
    repository = tmp_path / "repository"
    repository.mkdir()
    subprocess.run(["git", "init", "-q", str(repository)], check=True)
    subprocess.run(
        ["git", "-C", str(repository), "config", "user.email", "validator@test"],
        check=True,
    )
    subprocess.run(
        ["git", "-C", str(repository), "config", "user.name", "Validator Test"],
        check=True,
    )
    matrix_repository_path = VALIDATOR.MATRIX_REPOSITORY_PATH
    objects = {
        matrix_repository_path: b'{"cases": []}\n',
        "inputs/cgl_lf_paper/example.athinput": b"<job>\nbasename = example\n",
        VALIDATOR.PRODUCTION_UTILITY_REPOSITORY_PATH: b"print('utility')\n",
    }
    for relative, payload in objects.items():
        path = repository / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)
    subprocess.run(["git", "-C", str(repository), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repository), "commit", "-q", "-m", "objects"], check=True)
    revision = subprocess.run(
        ["git", "-C", str(repository), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    bundle = tmp_path / "objects.bundle"
    subprocess.run(
        ["git", "-C", str(repository), "bundle", "create", str(bundle), "HEAD"],
        check=True,
    )
    monkeypatch.setenv("GIT_OBJECT_DIRECTORY", str(tmp_path / "hostile-objects"))
    actual = VALIDATOR.authenticated_bundle_blobs(
        bundle,
        {
            "matrix": (revision, matrix_repository_path),
            "input": (revision, "inputs/cgl_lf_paper/example.athinput"),
            "production utility": (
                revision,
                VALIDATOR.PRODUCTION_UTILITY_REPOSITORY_PATH,
            ),
        },
    )
    assert actual == {
        "matrix": objects[matrix_repository_path],
        "input": objects["inputs/cgl_lf_paper/example.athinput"],
        "production utility": objects[VALIDATOR.PRODUCTION_UTILITY_REPOSITORY_PATH],
    }
    with pytest.raises(
        VALIDATOR.ValidationError,
        match="source bundle changed before object authentication",
    ):
        VALIDATOR.authenticated_bundle_blobs(
            bundle, {"matrix": (revision, matrix_repository_path)}, "0" * 64
        )
    original_git_run = VALIDATOR.git_run

    def report_revision_as_uncovered(arguments, **kwargs):
        if "merge-base" in arguments:
            return subprocess.CompletedProcess(arguments, 1, stdout=b"", stderr=b"")
        return original_git_run(arguments, **kwargs)

    monkeypatch.setattr(VALIDATOR, "git_run", report_revision_as_uncovered)
    with pytest.raises(
        VALIDATOR.ValidationError,
        match="source bundle advertised history omits matrix revision",
    ):
        VALIDATOR.authenticated_bundle_blobs(
            bundle, {"matrix": (revision, matrix_repository_path)}
        )


def test_validator_implementation_provenance_requires_committed_bytes(tmp_path):
    repository = tmp_path / "repository"
    path = repository / VALIDATOR.VALIDATOR_REPOSITORY_PATH
    path.parent.mkdir(parents=True)
    path.write_text("historical validator fixture\n")
    subprocess.run(["git", "init", "-q", str(repository)], check=True)
    subprocess.run(
        ["git", "-C", str(repository), "config", "user.email", "validator@test"],
        check=True,
    )
    subprocess.run(
        ["git", "-C", str(repository), "config", "user.name", "Validator Test"],
        check=True,
    )
    subprocess.run(["git", "-C", str(repository), "add", "."], check=True)
    subprocess.run(
        ["git", "-C", str(repository), "commit", "-q", "-m", "validator"],
        check=True,
    )
    revision = subprocess.run(
        ["git", "-C", str(repository), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    record = ORIGINAL_VALIDATOR_IMPLEMENTATION_PROVENANCE(path, {})
    assert record == {
        "path": str(path),
        "repository_path": VALIDATOR.VALIDATOR_REPOSITORY_PATH,
        "revision": revision,
        "sha256": sha256(path),
        "size_bytes": path.stat().st_size,
        "committed": True,
    }
    path.write_text("dirty historical validator fixture\n")
    with pytest.raises(
        VALIDATOR.ValidationError,
        match="validator source must be committed before validation",
    ):
        ORIGINAL_VALIDATOR_IMPLEMENTATION_PROVENANCE(path, {})


@pytest.mark.parametrize("retained_size", [True, "3"])
def test_authenticated_file_rejects_coercible_retained_size(tmp_path, retained_size):
    path = tmp_path / "artifact"
    path.write_bytes(b"abc")
    with pytest.raises(
        VALIDATOR.ValidationError,
        match="artifact retained size is not an exact integer",
    ):
        VALIDATOR.authenticated_file(path, sha256(path), "artifact", {}, retained_size)


def test_validator_rejects_future_snapshot_even_when_formally_recorded(segment_factory):
    fixture = segment_factory()
    add_future_snapshot(fixture, float(fixture["final_time"]) + 0.01)
    require_rejected(fixture, "segment retains a future snapshot")


def test_validator_rejects_duplicate_exact_terminal_snapshot(segment_factory):
    fixture = segment_factory()
    add_future_snapshot(fixture, float(fixture["final_time"]))
    require_rejected(
        fixture, "segment does not retain one unique exact terminal snapshot"
    )


def test_validator_rejects_duplicate_exact_terminal_restart(segment_factory):
    fixture = segment_factory()
    add_restart(fixture, float(fixture["final_time"]))
    require_rejected(
        fixture, "segment does not retain one unique exact terminal restart"
    )


def test_validator_rejects_tail_only_histories(segment_factory):
    fixture = segment_factory()
    remove_history_data_row(Path(fixture["mhd"]), 0)
    remove_history_data_row(Path(fixture["user"]), 0)
    refresh_inspection(fixture)
    require_rejected(fixture, "MHD history time does not begin at the prepared segment start")


def test_validator_requires_history_names_to_match_run_basename(segment_factory):
    fixture = segment_factory()
    path = Path(fixture["mhd"])
    path.rename(path.with_name("forged.mhd.hst"))
    require_rejected(
        fixture, "segment histories do not match the exact prepared run_basename"
    )


def test_validator_rejects_cross_read_history_rename_swap(
    segment_factory, monkeypatch
):
    fixture = segment_factory()
    original_parse_history = VALIDATOR.parse_history
    swapped = False

    def rename_swap(path, label, profiles):
        nonlocal swapped
        if label != "MHD history" or swapped:
            return original_parse_history(path, label, profiles)
        swapped = True
        backup = path.with_name(f"{path.name}.authenticated-version")
        payload = path.read_bytes()
        path.rename(backup)
        path.write_bytes(payload)
        try:
            return original_parse_history(path, label, profiles)
        finally:
            path.unlink()
            backup.rename(path)

    monkeypatch.setattr(VALIDATOR, "parse_history", rename_swap)
    require_rejected(fixture, "MHD history changed between authenticated reads")


@pytest.mark.parametrize(
    ("kind", "parser_name", "label"),
    [
        ("snapshot", "binary_snapshot_evidence", "binary snapshot"),
        ("restart", "restart_binary_evidence", "restart product"),
    ],
)
def test_validator_rejects_cross_read_product_rename_swap(
    segment_factory, monkeypatch, kind, parser_name, label
):
    fixture = segment_factory()
    target = Path(fixture[f"{kind}_groups"][0]["paths"][0])
    original_parser = getattr(VALIDATOR, parser_name)
    swapped = False

    def rename_swap(path, *args, **kwargs):
        nonlocal swapped
        if Path(path) != target or swapped:
            return original_parser(path, *args, **kwargs)
        swapped = True
        backup = path.with_name(f"{path.name}.authenticated-version")
        payload = path.read_bytes()
        path.rename(backup)
        path.write_bytes(payload)
        try:
            return original_parser(path, *args, **kwargs)
        finally:
            path.unlink()
            backup.rename(path)

    monkeypatch.setattr(VALIDATOR, parser_name, rename_swap)
    require_rejected(fixture, f"{label} changed between authenticated reads")


def test_validator_rejects_missing_history_cadence(segment_factory):
    fixture = segment_factory()
    remove_history_data_row(Path(fixture["mhd"]), 1)
    remove_history_data_row(Path(fixture["user"]), 1)
    refresh_inspection(fixture)
    require_rejected(fixture, "MHD history time does not provide expected cadence coverage")


def test_validator_rejects_missing_snapshot_cadence(segment_factory):
    fixture = segment_factory()
    missing = fixture["snapshot_groups"].pop(1)
    for path in missing["paths"]:
        Path(path).unlink()
    refresh_inspection(fixture)
    require_rejected(fixture, "snapshot times does not provide expected cadence coverage")


def test_validator_rejects_missing_restart_cadence(segment_factory):
    fixture = segment_factory(final_time=2.2)
    missing = fixture["restart_groups"].pop(1)
    for path in missing["paths"]:
        Path(path).unlink()
    refresh_inspection(fixture)
    require_rejected(fixture, "restart times does not provide expected cadence coverage")


def test_validator_rejects_nonincreasing_snapshot_cycles(segment_factory):
    fixture = segment_factory()
    terminal = fixture["snapshot_groups"][-1]
    encoded = f"  cycle={fixture_cycle(float(fixture['final_time']))}\n".encode()
    for path in terminal["paths"]:
        path = Path(path)
        payload = path.read_bytes()
        assert payload.count(encoded) == 1
        path.write_bytes(payload.replace(encoded, b"  cycle=0\n", 1))
    refresh_inspection(fixture)
    require_rejected(fixture, "snapshot cycles is not strictly increasing")


def test_validator_rejects_shared_time_cycle_disagreement(segment_factory):
    fixture = segment_factory()
    mutate_restart_header(
        fixture,
        248,
        struct.pack("<i", fixture_cycle(float(fixture["final_time"])) + 1),
        every_terminal_sibling=True,
    )
    require_rejected(
        fixture, "snapshot and restart cycles disagree at a shared physical time"
    )


def test_validator_rejects_restart_sibling_dt_disagreement(segment_factory):
    fixture = segment_factory()
    mutate_restart_header(fixture, 240, struct.pack("<d", 0.02))
    require_rejected(fixture, "restart sibling dt metadata disagree")


def test_validator_requires_continuation_histories_to_begin_at_parent_time(
    segment_factory,
):
    fixture = segment_factory()
    configure_continuation(fixture)
    for kind in ("snapshot", "restart"):
        for group in fixture[f"{kind}_groups"][:-1]:
            for path in group["paths"]:
                Path(path).unlink()
        fixture[f"{kind}_groups"] = fixture[f"{kind}_groups"][-1:]
    refresh_inspection(fixture)
    require_rejected(fixture, "MHD history time does not begin at the prepared segment start")


def test_validator_accepts_continuation_with_only_terminal_child_products(
    segment_factory,
):
    fixture = segment_factory()
    configure_continuation(fixture, start=0.3)
    for kind in ("snapshot", "restart"):
        for group in fixture[f"{kind}_groups"][:-1]:
            for path in group["paths"]:
                Path(path).unlink()
        fixture[f"{kind}_groups"] = fixture[f"{kind}_groups"][-1:]
    trim_history_before(Path(fixture["mhd"]), 0.3)
    trim_history_before(Path(fixture["user"]), 0.3)
    refresh_inspection(fixture)
    result = validate(fixture)
    assert result["snapshot_times"] == [fixture["final_time"]]
    assert result["restart_binary_times"] == [fixture["final_time"]]
    assert result["continuation_origin_validation"]["start_time"] == 0.3
    assert result["parent_acceptance_validation"]["status"] == (
        "trusted_authenticated_controller_record"
    )
    assert not result["parent_acceptance_validation"]["authorizing"]


def test_validator_rejects_terminal_only_continuation_snapshot_cadence_gap(
    segment_factory,
):
    fixture = segment_factory(final_time=0.75)
    configure_continuation(fixture, start=0.1)
    for kind in ("snapshot", "restart"):
        for group in fixture[f"{kind}_groups"][:-1]:
            for path in group["paths"]:
                Path(path).unlink()
        fixture[f"{kind}_groups"] = fixture[f"{kind}_groups"][-1:]
    trim_history_before(Path(fixture["mhd"]), 0.1)
    trim_history_before(Path(fixture["user"]), 0.1)
    refresh_inspection(fixture)
    require_rejected(fixture, "snapshot times does not provide expected cadence coverage")


def test_validator_rejects_terminal_only_continuation_restart_cadence_gap(
    segment_factory,
):
    fixture = segment_factory(final_time=1.25)
    configure_continuation(fixture, start=0.1)
    for group in fixture["snapshot_groups"][:1]:
        for path in group["paths"]:
            Path(path).unlink()
    fixture["snapshot_groups"] = fixture["snapshot_groups"][1:]
    for group in fixture["restart_groups"][:-1]:
        for path in group["paths"]:
            Path(path).unlink()
    fixture["restart_groups"] = fixture["restart_groups"][-1:]
    trim_history_before(Path(fixture["mhd"]), 0.1)
    trim_history_before(Path(fixture["user"]), 0.1)
    refresh_inspection(fixture)
    require_rejected(fixture, "restart times does not provide expected cadence coverage")


def test_validator_rejects_continuation_history_gap_without_start_products(
    segment_factory,
):
    fixture = segment_factory()
    configure_continuation(fixture, start=0.1)
    for kind in ("snapshot", "restart"):
        for group in fixture[f"{kind}_groups"][:-1]:
            for path in group["paths"]:
                Path(path).unlink()
        fixture[f"{kind}_groups"] = fixture[f"{kind}_groups"][-1:]
    trim_history_before(Path(fixture["mhd"]), 0.1)
    trim_history_before(Path(fixture["user"]), 0.1)
    remove_history_data_row(Path(fixture["mhd"]), 1)
    remove_history_data_row(Path(fixture["user"]), 1)
    refresh_inspection(fixture)
    require_rejected(fixture, "MHD history time does not provide expected cadence coverage")


def test_validator_binds_terminal_only_continuation_restart_order_to_parent(
    segment_factory,
):
    fixture = segment_factory()
    configure_continuation(fixture, start=0.1)
    for kind in ("snapshot", "restart"):
        for group in fixture[f"{kind}_groups"][:-1]:
            for path in group["paths"]:
                Path(path).unlink()
        fixture[f"{kind}_groups"] = fixture[f"{kind}_groups"][-1:]
    trim_history_before(Path(fixture["mhd"]), 0.1)
    trim_history_before(Path(fixture["user"]), 0.1)
    mutate_terminal_replicated_restart_state(
        fixture, "logical_location_inventory", every_sibling=True
    )
    require_rejected(
        fixture,
        "continuation restart logical-location order differs from authenticated parent",
    )


def test_validator_binds_terminal_only_continuation_location_cost_pairs_to_parent(
    segment_factory,
):
    fixture = segment_factory()
    configure_continuation(fixture, start=0.1)
    for kind in ("snapshot", "restart"):
        for group in fixture[f"{kind}_groups"][:-1]:
            for path in group["paths"]:
                Path(path).unlink()
        fixture[f"{kind}_groups"] = fixture[f"{kind}_groups"][-1:]
    trim_history_before(Path(fixture["mhd"]), 0.1)
    trim_history_before(Path(fixture["user"]), 0.1)
    mutate_terminal_replicated_restart_state(
        fixture, "cost_inventory_permutation", every_sibling=True
    )
    require_rejected(
        fixture,
        "continuation restart ordered logical-location/cost pairs differ from "
        "authenticated parent",
    )


@pytest.mark.parametrize("job_id", ("0", "0123", "job-123"))
def test_validator_requires_production_top_level_job_id(segment_factory, job_id):
    fixture = segment_factory(job_id=job_id)
    require_rejected(fixture, "manifest job ID is not a numeric job ID")


def test_validator_requires_valid_parent_segment_name(segment_factory):
    fixture = segment_factory()
    configure_continuation(fixture)
    mutate_manifest_value(
        fixture, ("command", "parent_segment", "segment"), "../parent"
    )
    require_rejected(fixture, "continuation parent segment name is invalid")


@pytest.mark.parametrize("job_id", ("0", "054321", "parent-job"))
def test_validator_requires_production_parent_job_id(segment_factory, job_id):
    fixture = segment_factory()
    configure_continuation(fixture)
    child = json.loads(Path(fixture["manifest"]).read_text())
    parent_path = Path(child["command"]["parent_segment"]["manifest"])
    parent = json.loads(parent_path.read_text())
    parent["job_id"] = job_id
    write_json(parent_path, parent)
    require_rejected(fixture, "continuation parent job ID is not a numeric job ID")


def test_validator_rejects_parent_continuation_identity_rebind(segment_factory):
    fixture = segment_factory()
    configure_continuation(fixture)
    mutate_manifest_value(
        fixture, ("command", "parent_segment", "case_id"), "R05"
    )
    require_rejected(fixture, "prepared continuation parent identity is inconsistent")


def test_validator_rejects_copied_continuation_restart_not_from_parent(
    segment_factory,
):
    fixture = segment_factory()
    configure_continuation(fixture)
    manifest_path = Path(fixture["manifest"])
    manifest = json.loads(manifest_path.read_text())
    record = manifest["command"]["restart_files"][-1]
    copied = Path(record["path"])
    copied.write_bytes(copied.read_bytes() + b"not the parent restart")
    manifest["command"]["restart_files"][-1] = file_record(copied)
    write_json(manifest_path, manifest)
    refresh_batch_script(fixture)
    require_rejected(
        fixture, "copied continuation restart bytes differ from recorded parent"
    )


@pytest.mark.parametrize(
    ("snapshots", "restarts", "message"),
    [
        (
            [{"time": 0.1, "cycle": 99}],
            [{"time": 0.1, "cycle": 100, "dt": 0.01}],
            "continuation exact-start output cycles differ from authenticated initial state",
        ),
        (
            [{"time": 0.1, "cycle": 100}],
            [{"time": 0.1, "cycle": 100, "dt": 0.02}],
            "continuation exact-start restart dt differs from authenticated initial state",
        ),
        (
            [{"time": 0.0, "cycle": 100}],
            [{"time": 0.1, "cycle": 100, "dt": 0.01}],
            "continuation retains output before authenticated start time",
        ),
    ],
)
def test_validator_binds_continuation_first_retained_state(
    snapshots, restarts, message
):
    continuation = {
        "start_time": 0.1,
        "initial_restart_cycle": 100,
        "initial_restart_dt": 0.01,
    }
    with pytest.raises(VALIDATOR.ValidationError, match=re.escape(message)):
        VALIDATOR.require_continuation_initial_state(
            continuation, snapshots, restarts
        )


@pytest.mark.parametrize(("kind", "message"), [
    ("snapshot", "binary snapshot field payload is truncated"),
    ("restart", "restart meshblock payload is truncated"),
])
def test_validator_rejects_truncated_real_payload_products(
    segment_factory, kind, message
):
    fixture = segment_factory()
    groups = fixture[f"{kind}_groups"]
    terminal_path = Path(groups[-1]["paths"][-1])
    terminal_path.write_bytes(terminal_path.read_bytes()[:-1])
    refresh_inspection(fixture)
    require_rejected(fixture, message)


@pytest.mark.parametrize(
    ("kind", "message"),
    [
        ("snapshot", "binary snapshot field payload contains a non-finite value"),
        ("restart", "restart meshblock payload contains a non-finite value"),
    ],
)
def test_validator_rejects_nonfinite_complete_payloads(segment_factory, kind, message):
    fixture = segment_factory()
    path = Path(fixture[f"{kind}_groups"][-1]["paths"][-1])
    payload = bytearray(path.read_bytes())
    encoded = struct.pack("<f", float("nan")) if kind == "snapshot" else struct.pack(
        "<d", float("nan")
    )
    payload[-len(encoded):] = encoded
    path.write_bytes(payload)
    refresh_inspection(fixture)
    require_rejected(fixture, message)


@pytest.mark.parametrize("kind", ("snapshot", "restart"))
def test_validator_rejects_impossible_logical_meshblock_coordinates(
    segment_factory, kind
):
    fixture = segment_factory()
    path = Path(fixture[f"{kind}_groups"][-1]["paths"][-1])
    payload = bytearray(path.read_bytes())
    if kind == "snapshot":
        original = struct.pack("<10i", 1, 1, 0, 0, 0, 0, 3, 0, 0, 0)
        replacement = struct.pack("<10i", 1, 1, 0, 0, 0, 0, 4, 0, 0, 0)
        assert payload.count(original) == 1
        payload = payload.replace(original, replacement, 1)
        message = "binary snapshot logical inventory is invalid"
    else:
        parameter_end = payload.index(b"<par_end>\n") + len(b"<par_end>\n")
        location_offset = parameter_end + int(TEST_ABI["mesh_header_size"])
        payload[location_offset:location_offset + 16] = struct.pack("<4i", 4, 0, 0, 0)
        message = "restart logical-location inventory is invalid"
    path.write_bytes(payload)
    refresh_inspection(fixture)
    require_rejected(fixture, message)


@pytest.mark.parametrize("layout", ("per_rank", "shared_mpiio"))
def test_validator_rejects_valid_snapshot_logical_location_permutation(
    segment_factory, monkeypatch, layout
):
    qualify_fixture_layout(monkeypatch, layout)
    fixture = segment_factory(layout=layout)
    swap_first_two_snapshot_blocks(Path(fixture["snapshot_groups"][-1]["paths"][-1]))
    refresh_inspection(fixture)
    require_rejected(
        fixture,
        "snapshot ordered logical-location inventories differ across outputs",
    )


def test_validator_rejects_shared_terminal_only_cross_product_order_permutation():
    snapshots = [{
        "time": 0.5,
        "ordered_location_inventory": (((1, 0, 0, 0), (0, 0, 0, 0)),),
    }]
    restarts = [{
        "time": 0.5,
        "ordered_location_inventory": ((0, 0, 0, 0), (1, 0, 0, 0)),
    }]
    with pytest.raises(
        VALIDATOR.ValidationError,
        match="shared snapshot and restart ordered logical-location inventories disagree",
    ):
        VALIDATOR.require_matching_shared_output_orders(snapshots, restarts)


@pytest.mark.parametrize("layout", ("per_rank", "shared_mpiio"))
def test_validator_rejects_valid_restart_logical_location_and_cost_permutation(
    segment_factory, monkeypatch, layout
):
    qualify_fixture_layout(monkeypatch, layout)
    fixture = segment_factory(layout=layout)
    mutate_terminal_replicated_restart_state(
        fixture, "logical_location_inventory", every_sibling=True
    )
    require_rejected(
        fixture,
        (
            "restart ordered logical-location inventories differ across outputs"
            if layout == "per_rank"
            else "shared snapshot and restart ordered logical-location inventories disagree"
        ),
    )


@pytest.mark.parametrize("layout", ("per_rank", "shared_mpiio"))
def test_validator_rejects_cost_only_permutation_across_restart_times(
    segment_factory, monkeypatch, layout
):
    qualify_fixture_layout(monkeypatch, layout)
    fixture = segment_factory(layout=layout)
    mutate_terminal_replicated_restart_state(
        fixture, "cost_inventory_permutation", every_sibling=True
    )
    require_rejected(
        fixture,
        "restart ordered logical-location/cost pairs differ across outputs",
    )


@pytest.mark.parametrize(
    ("field", "message"),
    [
        (
            "logical_location_inventory",
            "restart sibling ordered_location_inventory metadata disagree",
        ),
        ("cost_inventory", "restart sibling ordered_cost_inventory metadata disagree"),
        (
            "turbulence_n_updates",
            "restart sibling turbulence_metadata_sha256 metadata disagree",
        ),
        (
            "rng_idum",
            "restart sibling rng_canonical_continuation_payload metadata disagree",
        ),
        (
            "rng_idum2",
            "restart sibling rng_canonical_continuation_payload metadata disagree",
        ),
        (
            "rng_iy",
            "restart sibling rng_canonical_continuation_payload metadata disagree",
        ),
        (
            "rng_iv",
            "restart sibling rng_canonical_continuation_payload metadata disagree",
        ),
        (
            "rng_iset",
            "restart sibling rng_canonical_continuation_payload metadata disagree",
        ),
        (
            "turbulence_amplitudes",
            "restart sibling turbulence_amplitudes_sha256 metadata disagree",
        ),
    ],
)
def test_validator_requires_exact_replicated_restart_sibling_state(
    segment_factory, field, message
):
    fixture = segment_factory()
    mutate_terminal_replicated_restart_state(fixture, field)
    require_rejected(fixture, message)


@pytest.mark.parametrize(
    ("field", "name"),
    [
        ("turbulence_nlow", "nlow"),
        ("turbulence_projection_policy", "projection_policy"),
        ("turbulence_tcorr", "tcorr"),
    ],
)
def test_validator_binds_turbulence_configuration_to_frozen_input(
    segment_factory, field, name
):
    fixture = segment_factory()
    mutate_terminal_replicated_restart_state(fixture, field, every_sibling=True)
    require_rejected(
        fixture,
        f"restart turbulence configuration differs from frozen input for '{name}'",
    )


def test_validator_accepts_and_records_rng_native_padding_divergence(segment_factory):
    fixture = segment_factory()
    mutate_terminal_replicated_restart_state(fixture, "rng_padding")
    result = validate(fixture)
    terminal = result["restart_rng_native_padding_evidence"][-1]
    assert {item["offset"] for item in terminal} == {284}
    assert {item["size_bytes"] for item in terminal} == {4}
    assert all(item["ignored_for_continuation_authentication"] for item in terminal)
    assert len({item["hex"] for item in terminal}) == 2
    continuation = result["restart_replicated_state_sha256"][-1]
    assert "rng_canonical_continuation_state_sha256" in continuation
    assert "rng_state_sha256" not in continuation


def test_validator_authenticates_only_seed_and_records_pre_update_ignored_rng_state(
    segment_factory,
):
    fixture = segment_factory()
    for field in ("rng_idum2", "rng_iy", "rng_iv", "rng_inactive_nonfinite_gset", "rng_padding"):
        mutate_terminal_replicated_restart_state(fixture, field, group_index=0)
    result = validate(fixture)
    lifecycle = result["restart_rng_lifecycle_evidence"][0]
    assert lifecycle["lifecycle"] == "pre_update_seed"
    assert lifecycle["n_updates"] == 0
    assert lifecycle["authenticated_field_count"] == 1
    ignored = result["restart_rng_pre_update_ignored_state_evidence"][0]
    assert len(ignored) == len(fixture["restart_groups"][0]["paths"])
    for key in (
        "idum2_raw_hex",
        "iy_raw_hex",
        "iv_raw_hex",
        "gset_raw_hex",
        "native_padding_raw_hex",
    ):
        assert len({item[key] for item in ignored}) == 2
    assert all(item["ignored_for_continuation_authentication"] for item in ignored)


def test_validator_authenticates_pre_update_rng_seed(segment_factory):
    fixture = segment_factory()
    mutate_terminal_replicated_restart_state(fixture, "rng_idum", group_index=0)
    require_rejected(
        fixture,
        "restart sibling rng_canonical_continuation_payload metadata disagree",
    )


def test_restart_rng_evidence_authenticates_initialized_continuation_state():
    semantic_state = tuple(range(-17, 18)) + (1, 0.5)
    semantic_state = (17,) + semantic_state[1:]
    payload = bytearray(struct.pack("<35qi4xd", *semantic_state))
    evidence = VALIDATOR.restart_rng_evidence(
        bytes(payload), TEST_ABI, Path("/fixture/a"), 1
    )
    payload[284:288] = b"\x01\x23\x45\x67"
    padding_variant = VALIDATOR.restart_rng_evidence(
        bytes(payload), TEST_ABI, Path("/fixture/b"), 1
    )
    assert evidence["rng_lifecycle"] == "initialized_continuation"
    assert evidence["rng_authenticated_field_count"] == len(semantic_state) == 37
    assert evidence["rng_gset_active"]
    assert evidence["rng_gset_authentication"] == "authenticated_exact_finite"
    assert evidence["rng_canonical_continuation_state_sha256"] == padding_variant[
        "rng_canonical_continuation_state_sha256"
    ]
    assert evidence["rng_native_padding_hex"] != padding_variant[
        "rng_native_padding_hex"
    ]


@pytest.mark.parametrize(
    ("field", "group_index"),
    [
        ("rng_positive_idum", 0),
        ("rng_iset", 0),
        ("rng_negative_idum", -1),
        ("rng_zero_idum", -1),
    ],
)
def test_validator_rejects_inconsistent_rng_lifecycle(
    segment_factory, field, group_index
):
    fixture = segment_factory()
    mutate_terminal_replicated_restart_state(
        fixture, field, group_index=group_index
    )
    require_rejected(fixture, "restart RNG lifecycle state is invalid")


def test_validator_accepts_and_records_divergent_inactive_rng_gset(segment_factory):
    fixture = segment_factory()
    mutate_terminal_replicated_restart_state(fixture, "rng_inactive_nonfinite_gset")
    result = validate(fixture)
    terminal = result["restart_rng_gset_evidence"][-1]
    assert {item["active"] for item in terminal} == {False}
    assert {item["authentication"] for item in terminal} == {
        "ignored_inactive_uninitialized"
    }
    assert all(item["ignored_for_continuation_authentication"] for item in terminal)
    assert len({item["raw_hex"] for item in terminal}) == 2


def test_validator_rejects_changed_active_rng_gset(segment_factory):
    fixture = segment_factory()
    mutate_terminal_replicated_restart_state(fixture, "rng_iset", every_sibling=True)
    mutate_terminal_replicated_restart_state(fixture, "rng_gset")
    require_rejected(
        fixture,
        "restart sibling rng_canonical_continuation_payload metadata disagree",
    )


@pytest.mark.parametrize(
    ("offset", "encoded", "message"),
    [
        (280, struct.pack("<i", 2), "restart RNG semantic state is invalid"),
    ],
)
def test_validator_rejects_invalid_rng_iset(
    segment_factory, offset, encoded, message
):
    fixture = segment_factory()
    terminal = fixture["restart_groups"][-1]
    path = Path(terminal["paths"][-1])
    payload = bytearray(path.read_bytes())
    parameter_end = payload.index(b"<par_end>\n") + len(b"<par_end>\n")
    rng = (
        parameter_end
        + int(TEST_ABI["mesh_header_size"])
        + len(logical_locations()) * 16
        + len(logical_locations()) * 4
        + int(TEST_ABI["turbulence_metadata_size"])
    )
    payload[rng + offset:rng + offset + len(encoded)] = encoded
    path.write_bytes(payload)
    refresh_inspection(fixture)
    require_rejected(fixture, message)


def test_validator_rejects_nonfinite_active_rng_gset(segment_factory):
    fixture = segment_factory()
    mutate_terminal_replicated_restart_state(fixture, "rng_active_nonfinite_gset")
    require_rejected(fixture, "restart RNG semantic state is invalid")


def test_validator_rejects_snapshot_active_region_index_mutation(segment_factory):
    fixture = segment_factory()
    path = Path(fixture["snapshot_groups"][-1]["paths"][-1])
    payload = bytearray(path.read_bytes())
    original = struct.pack("<10i", 1, 1, 0, 0, 0, 0, 3, 0, 0, 0)
    replacement = struct.pack("<10i", 1, 2, 0, 0, 0, 0, 3, 0, 0, 0)
    assert payload.count(original) == 1
    path.write_bytes(payload.replace(original, replacement, 1))
    refresh_inspection(fixture)
    require_rejected(fixture, "binary snapshot active-region indices differ from input")


def test_validator_rejects_snapshot_finite_geometry_mutation(segment_factory):
    fixture = segment_factory()
    path = Path(fixture["snapshot_groups"][-1]["paths"][-1])
    payload = bytearray(path.read_bytes())
    indices = struct.pack("<10i", 1, 1, 0, 0, 0, 0, 3, 0, 0, 0)
    start = payload.index(indices) + len(indices)
    payload[start:start + 48] = struct.pack("<6d", 0.75, 0.99, 0.0, 1.0, 0.0, 2.0)
    path.write_bytes(payload)
    refresh_inspection(fixture)
    require_rejected(fixture, "binary snapshot geometry:")


def test_validator_rejects_restart_active_mesh_region_index_mutation(segment_factory):
    fixture = segment_factory()
    mutate_restart_header(fixture, 80 + 5 * 4, struct.pack("<i", 5))
    require_rejected(fixture, "restart mesh indices differ from archived input")


def test_validator_rejects_restart_meshblock_coarse_region_index_mutation(
    segment_factory,
):
    fixture = segment_factory()
    mutate_restart_header(fixture, 156 + 10 * 4, struct.pack("<i", 2))
    require_rejected(fixture, "restart meshblock indices differ from archived input")


def test_validator_rejects_restart_finite_root_geometry_mutation(segment_factory):
    fixture = segment_factory()
    mutate_restart_header(fixture, 8 + 3 * 8, struct.pack("<d", 1.25))
    require_rejected(fixture, "restart mesh geometry:")


def test_validator_discloses_qualified_opaque_mesh_coarse_region_fields(
    segment_factory,
):
    fixture = segment_factory()
    mutate_restart_header(fixture, 80 + 10 * 4, struct.pack("<i", 123456))
    result = validate(fixture)
    assert (
        result["product_payload_validation"]["legacy_mesh_coarse_region_fields"]
        == "legacy_opaque_uninitialized"
    )
    assert not result["independent_authorization"]["authorizing"]


@pytest.mark.parametrize("layout", ("per_rank", "shared_mpiio"))
def test_validator_binds_terminal_restart_injected_work_to_history(
    segment_factory, monkeypatch, layout
):
    if layout == "shared_mpiio":
        monkeypatch.setitem(
            VALIDATOR.COMMON_FROZEN_PARAMETERS,
            ("output2", "single_file_per_rank"),
            "false",
        )
        monkeypatch.setitem(
            VALIDATOR.COMMON_FROZEN_PARAMETERS,
            ("output3", "single_file_per_rank"),
            "false",
        )
    fixture = segment_factory(layout=layout)
    mutate_terminal_restart_state(fixture, injected_work=0.251)
    require_rejected(fixture, "terminal restart injected work differs from user history")


@pytest.mark.parametrize(
    ("index", "value", "name"),
    [
        (0, 4.0, "lf_nstage"),
        (13, 0.11, "lf_qprwrk"),
    ],
)
def test_validator_binds_terminal_restart_lf_diagnostics_to_history(
    segment_factory, index, value, name
):
    fixture = segment_factory()
    mutate_terminal_restart_state(fixture, diagnostics={index: value})
    require_rejected(fixture, f"terminal restart {name} differs from MHD history")


def test_validator_rejects_restart_parameter_contract_mutation(segment_factory):
    fixture = segment_factory()
    path = Path(fixture["restart_groups"][-1]["paths"][-1])
    payload = path.read_bytes()
    assert payload.count(b"beta0 = 10") == 1
    path.write_bytes(payload.replace(b"beta0 = 10", b"beta0 = 11", 1))
    refresh_inspection(fixture)
    require_rejected(fixture, "restart parameter dump violates qualified problem/beta0=10")


def test_validator_rejects_unqualified_archived_input_parameter(segment_factory):
    fixture = segment_factory()
    path = Path(fixture["run_dir"]) / "manifest" / "submitted_input.athinput"
    payload = path.read_text()
    assert payload.count("<mhd>\n") == 1
    path.write_text(payload.replace("<mhd>\n", "<mhd>\nfofc = 1\n", 1))
    bind_fixture_artifact(fixture, "input", path, "input_sha256")
    require_rejected(fixture, "archived input contains unqualified parameters")


@pytest.mark.parametrize(
    ("replacement", "message"),
    [
        (b"fofc = 1", "restart parameter dump violates qualified mhd/fofc=0"),
        (
            b"evil = 0",
            "restart parameter dump parameter inventory differs from qualified contract",
        ),
    ],
)
def test_validator_rejects_product_parameter_contract_escape(
    segment_factory, replacement, message
):
    fixture = segment_factory()
    path = Path(fixture["restart_groups"][-1]["paths"][-1])
    payload = path.read_bytes()
    assert payload.count(b"fofc = 0") == 1
    path.write_bytes(payload.replace(b"fofc = 0", replacement, 1))
    refresh_inspection(fixture)
    require_rejected(fixture, message)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ({"mass_end": 2.01}, "mass conservation exceeds named policy"),
        ({"hard_vol_end": 1.0}, "hard_vol is nonzero"),
    ],
)
def test_validator_enforces_mass_drift_and_hard_volume(segment_factory, mutation, message):
    fixture = segment_factory()
    write_histories(fixture, **mutation)
    refresh_inspection(fixture)
    require_rejected(fixture, message)


def test_validator_enforces_mass_history_mismatch(segment_factory):
    fixture = segment_factory()
    write_histories(fixture, mass_end=2.0 + 8.0e-13, user_mass_end=2.0 - 8.0e-13)
    refresh_inspection(fixture)
    assert validate(fixture)["mass_relative_mismatch"] <= 1.0e-12
    write_histories(fixture, mass_end=2.0 + 1.5e-12, user_mass_end=2.0 - 1.5e-12)
    refresh_inspection(fixture)
    require_rejected(fixture, "mass-history mismatch exceeds named policy")


def test_validator_scopes_absent_normalized_ct_divb_to_historical_executable(
    segment_factory,
):
    fixture = segment_factory()
    result = validate(fixture)
    record = result["normalized_ct_divb_validation"]
    assert record["status"] == "historically_unavailable"
    assert record["executable_revision"] == EXECUTABLE_REVISION
    assert record["executable_sha256"] == TEST_EXECUTABLE_SHA256
    assert "max_normalized_divb" not in result


@pytest.mark.parametrize(
    ("fixture_key", "column", "label"),
    [
        ("user", "max_ndiv", "user history"),
        ("mhd", "future_mhd_diagnostic", "MHD history"),
    ],
)
def test_validator_rejects_future_build_history_schema(
    segment_factory, fixture_key, column, label
):
    fixture = segment_factory()
    append_history_column(Path(fixture[fixture_key]), column, 0.0)
    refresh_inspection(fixture)
    require_rejected(
        fixture,
        f"{label} columns differ from qualified historical executable schema",
    )


@pytest.mark.parametrize("case_id", sorted(VALIDATOR.FINITE_LIMITER_CASE_IDS))
def test_validator_enforces_finite_limiter_contract(segment_factory, case_id):
    fixture = segment_factory(case_id=case_id)
    result = validate(fixture)
    assert result["limiter_mode"] == "finite_collisional"
    assert result["final_hardwall_projection_count"] == 0.0
    write_histories(fixture, hardwall_count=1.0)
    refresh_inspection(fixture)
    require_rejected(fixture, "finite-limiter case contains hardwall projections")


@pytest.mark.parametrize(
    ("case_id", "hardwall", "message"),
    [
        ("R04", False, "archived input lacks required parameter mhd/limiter_hardwall"),
        ("R14", True, "archived input contains forbidden parameter mhd/limiter_hardwall"),
    ],
)
def test_validator_rejects_limiter_mode_mapping_mismatch(
    segment_factory, case_id, hardwall, message
):
    fixture = segment_factory(case_id=case_id, hardwall=hardwall)
    require_rejected(fixture, message)


@pytest.mark.parametrize(
    ("case_id", "original", "replacement", "message"),
    [
        ("R17", "nx1 = 4", "nx1 = 8", "archived input violates frozen mesh/nx1=4"),
        ("R04", "beta0 = 10", "beta0 = 11", "archived input violates frozen problem/beta0=10"),
        (
            "R04",
            "projection_policy = mks24_random_unprojected",
            "projection_policy = mks24_alfvenic_perpendicular",
            "archived input violates frozen turb_driving/projection_policy=mks24_random_unprojected",
        ),
        (
            "R12",
            "lf_k_parallel = 0.062831853071795854",
            "lf_k_parallel = 6.283185307179586",
            "archived input violates frozen mhd/lf_k_parallel=0.062831853071795854",
        ),
        (
            "R14",
            "limiter_nu_coll = 20",
            "limiter_nu_coll = 200",
            "archived input violates frozen mhd/limiter_nu_coll=20",
        ),
        (
            "R04",
            "dt = 0.02",
            "dt = 0.03",
            "archived input violates frozen output1/dt=0.02",
        ),
    ],
)
def test_validator_enforces_exact_frozen_case_contracts(
    segment_factory, case_id, original, replacement, message
):
    fixture = segment_factory(case_id=case_id)
    replace_bound_input(fixture, original, replacement)
    require_rejected(fixture, message)


def test_validator_enforces_exact_frozen_matrix_mapping(segment_factory):
    fixture = segment_factory(case_id="R17")
    matrix_path = Path(fixture["matrix"])
    matrix = json.loads(matrix_path.read_text())
    next(item for item in matrix["cases"] if item["id"] == "R17")["resolution"] = "8x1x1"
    write_json(matrix_path, matrix)
    manifest_path = Path(fixture["manifest"])
    manifest = json.loads(manifest_path.read_text())
    manifest["run"]["resolution"] = "8x1x1"
    write_json(manifest_path, manifest)
    bind_fixture_artifact(fixture, "matrix", matrix_path, "matrix_sha256")
    require_rejected(fixture, "retained matrix does not match frozen R17 mapping")


def test_validator_allows_zero_but_rejects_invalid_hardwall_counts(segment_factory):
    fixture = segment_factory()
    write_histories(fixture, hardwall_count=0.0)
    mutate_terminal_restart_state(fixture, diagnostics={17: 0.0})
    assert validate(fixture)["final_hardwall_projection_count"] == 0.0
    write_histories(fixture, hardwall_count=-1.0)
    refresh_inspection(fixture)
    require_rejected(fixture, "hardwall projection diagnostic contains an invalid count")
    write_histories(fixture, hardwall_start=2.0, hardwall_count=1.0)
    refresh_inspection(fixture)
    require_rejected(fixture, "hardwall projection diagnostic is not monotonic")


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ({"mirror_end": -1.0}, "lf_mirror contains an invalid count"),
        (
            {"firehose_start": 2.0, "firehose_end": 1.0},
            "lf_firehs is not monotonic",
        ),
    ],
)
def test_validator_enforces_cumulative_limiter_counters(
    segment_factory, mutation, message
):
    fixture = segment_factory()
    write_histories(fixture, **mutation)
    refresh_inspection(fixture)
    require_rejected(fixture, message)


def test_validator_enforces_cap_increments(segment_factory):
    fixture = segment_factory()
    write_histories(
        fixture,
        qface_start=5.0,
        qface_end=6.0,
        cap_start=0.0,
        cap_end=2.0,
    )
    refresh_inspection(fixture)
    require_rejected(fixture, "lf_qprcap increment exceeds lf_qface increment")


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ({"strict_failure_end": 1.0}, "strict LF failure counter is nonzero"),
        ({"nstage_end": 0.0}, "LF stage or qface diagnostics did not advance"),
        (
            {"nstage_start": 2.0, "nstage_end": 2.0},
            "LF stage or qface diagnostics did not advance",
        ),
        (
            {"qface_end": 0.0, "cap_end": 0.0},
            "LF stage or qface diagnostics did not advance",
        ),
        (
            {
                "qface_start": 6.0,
                "qface_end": 6.0,
                "cap_start": 1.0,
                "cap_end": 1.0,
            },
            "LF stage or qface diagnostics did not advance",
        ),
    ],
)
def test_validator_enforces_universal_lf_gates(segment_factory, mutation, message):
    fixture = segment_factory()
    write_histories(fixture, **mutation)
    refresh_inspection(fixture)
    require_rejected(fixture, message)


def test_validator_rejects_zero_active_pressure_work(segment_factory):
    fixture = segment_factory()
    write_histories(
        fixture,
        qprwrk_end=0.0,
        qpewrk_end=0.0,
        cpwork_end=0.0,
        cawork_end=0.0,
    )
    refresh_inspection(fixture)
    require_rejected(
        fixture, "active pressure work did not exceed the named activity policy"
    )


def test_validator_rejects_zero_active_forcing_work(segment_factory):
    fixture = segment_factory()
    write_histories(fixture, force_work_end=0.0)
    refresh_inspection(fixture)
    require_rejected(
        fixture, "forcing work did not exceed the named activity policy"
    )


@pytest.mark.parametrize(
    ("case_id", "mutation", "message"),
    [
        (
            "R04",
            {"force_work_end": 5.0e-7},
            "forcing work did not exceed the named activity policy",
        ),
        (
            "R06",
            {"force_work_end": 5.0e-7},
            "forcing work did not exceed the named activity policy",
        ),
        (
            "R04",
            {
                "energy_start": 1.0e8,
                "energy_end": 1.0e8 + 5.0e-5,
                "force_work_end": 5.0e-5,
            },
            "forcing work did not exceed the named activity policy",
        ),
        (
            "R04",
            {"cpwork_end": 5.0e-7, "cawork_end": 0.0},
            "active pressure work did not exceed the named activity policy",
        ),
        (
            "R04",
            {
                "energy_start": 1.0e8,
                "energy_end": 1.0e8 + 2.0,
                "force_work_end": 2.0,
                "cpwork_end": 5.0e-5,
                "cawork_end": 0.0,
            },
            "active pressure work did not exceed the named activity policy",
        ),
    ],
)
def test_validator_requires_named_absolute_and_state_normalized_activity(
    segment_factory, case_id, mutation, message
):
    fixture = segment_factory(case_id=case_id)
    write_histories(fixture, **mutation)
    refresh_inspection(fixture)
    require_rejected(fixture, message)


@pytest.mark.parametrize(
    "column", ("qprwrk_end", "qpewrk_end", "cpwork_end", "cawork_end")
)
def test_validator_requires_all_lf_work_to_be_finite(segment_factory, column):
    fixture = segment_factory()
    write_histories(fixture, **{column: float("nan")})
    refresh_inspection(fixture)
    require_rejected(fixture, "MHD history contains a non-finite row")


def test_validator_enforces_only_approved_standard_active_closure(segment_factory):
    fixture = segment_factory()
    write_histories(fixture, force_work_end=0.249999995)
    mutate_terminal_restart_state(fixture, injected_work=0.249999995)
    assert validate(fixture)["sampled_history_forcing_work"]["normalized_residual"] < 1.0e-8
    write_histories(fixture, force_work_end=0.24999998)
    refresh_inspection(fixture)
    require_rejected(
        fixture, "active forcing-work closure normalized residual is too large"
    )


def test_validator_requires_exact_partial_result(segment_factory):
    fixture = segment_factory(final_time=0.5, required_time=1.0)
    require_rejected(
        fixture, "accepted validation requires the exact prepared endpoint"
    )
    result = validate(fixture, result="clean_partial")
    assert result["segment_result"] == "clean_partial"


@pytest.mark.skipif(
    not RETAINED_R03_MANIFEST.is_file() or not retained_r03_provenance_is_intact(),
    reason="retained canonical R03 or its mutable launch provenance has drifted",
)
def test_validator_accepts_retained_real_r03_rng_lifecycle_evidence(monkeypatch):
    use_production_contracts(monkeypatch)
    result = VALIDATOR.validate_segment(
        argparse.Namespace(
            manifest=RETAINED_R03_MANIFEST,
            inspection=None,
            policy="r03-4762472-clean-partial-v1",
            result="clean_partial",
        )
    )
    assert result["validation_accepted"]
    assert result["case_id"] == "R03"
    assert result["segment_result"] == "clean_partial"
    assert all(
        item["offset"] == 284
        and item["size_bytes"] == 4
        and item["ignored_for_continuation_authentication"]
        for group in result["restart_rng_native_padding_evidence"]
        for item in group
    )
    assert any(
        item["authentication"] == "ignored_inactive_uninitialized"
        and item["ignored_for_continuation_authentication"]
        for group in result["restart_rng_gset_evidence"]
        for item in group
    )
    assert result["restart_rng_lifecycle_evidence"][0]["lifecycle"] == "pre_update_seed"
    assert result["restart_rng_lifecycle_evidence"][0]["authenticated_field_count"] == 1
    assert result["restart_rng_pre_update_ignored_state_evidence"][0]
    assert all(
        item["lifecycle"] == "initialized_continuation"
        for item in result["restart_rng_lifecycle_evidence"][1:]
    )


@pytest.mark.skipif(
    not RETAINED_R03_MANIFEST.is_file() or retained_r03_provenance_is_intact(),
    reason="retained canonical R03 mutable launch provenance has not drifted",
)
def test_validator_rejects_retained_r03_launch_provenance_drift(monkeypatch):
    use_production_contracts(monkeypatch)
    with pytest.raises(
        VALIDATOR.ValidationError,
        match="(prepared executable|qualification approval) SHA-256 differs from retained provenance",
    ):
        VALIDATOR.validate_segment(
            argparse.Namespace(
                manifest=RETAINED_R03_MANIFEST,
                inspection=None,
                policy="r03-4762472-clean-partial-v1",
                result="clean_partial",
            )
        )


@pytest.mark.skipif(
    not RETAINED_R03_MANIFEST.is_file() or retained_r03_provenance_is_intact(),
    reason="non-authorizing compatibility diagnostic is unnecessary",
)
def test_retained_r03_products_remain_compatible_beyond_executable_drift(monkeypatch):
    use_production_contracts(monkeypatch)
    retained_manifest = json.loads(RETAINED_R03_MANIFEST.read_text())
    retained_approval = retained_manifest["command"]["qualification_approval"]["token"]
    original = VALIDATOR.authenticated_file
    original_load_json = VALIDATOR.load_json

    def bypass_only_missing_historical_executable(
        path, expected_sha256, label, profiles, expected_size=None
    ):
        if label in {"prepared executable", "qualification approval"}:
            return {
                "path": str(path),
                "size_bytes": path.stat().st_size,
                "sha256": expected_sha256,
            }
        return original(path, expected_sha256, label, profiles, expected_size)

    def use_retained_historical_approval(path, label, profiles):
        if label == "qualification approval":
            return retained_approval, retained_manifest["command"][
                "qualification_approval"
            ]["sha256"]
        return original_load_json(path, label, profiles)

    monkeypatch.setattr(
        VALIDATOR, "authenticated_file", bypass_only_missing_historical_executable
    )
    monkeypatch.setattr(VALIDATOR, "load_json", use_retained_historical_approval)
    result = VALIDATOR.validate_segment(
        argparse.Namespace(
            manifest=RETAINED_R03_MANIFEST,
            inspection=None,
            policy="r03-4762472-clean-partial-v1",
            result="clean_partial",
        )
    )
    assert result["validation_accepted"]
    assert not result["independent_authorization"]["authorizing"]
