#!/usr/bin/env python3
"""Assemble and report the corrected Stage I campaign with legacy passive controls.

Corrected active cases are admitted only from the isolated corrected campaign
identity.  R06-R09 are imported only as authenticated legacy passive controls:
their source and executable identities remain explicit and are never relabeled
as corrected production.  The adapter writes ordinary direct-fast report case
products plus a composite inventory that records the execution authority and
evidence class of every case.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from typing import Iterable


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[2]
REPORTER_PATH = SCRIPT_PATH.with_name("cgl_lf_stage_i_fast_report.py")
IDENTITY_TOOL_PATH = SCRIPT_PATH.with_name("cgl_lf_stage_i_corrected_identity.py")

ACTIVE_CASES = (
    "R02",
    "R03",
    "R04",
    "R05",
    "R10",
    "R11",
    "R12",
    "R13",
    "R14",
    "R15",
    "R16",
    "R17",
)
PASSIVE_CONTROLS = ("R06", "R07", "R08", "R09")
ALL_CASES = tuple(f"R{number:02d}" for number in range(2, 18))
R14_FAILURE_CASE = "R14"
R14_FAILURE_TIME_TOLERANCE = 1.0e-12
R14_FATAL_SIGNATURE = (
    "### FATAL ERROR in turbulence driver: "
    "cannot inject non-zero dedt with a zero forcing field"
)
R14_STRICT_OVERRIDE = "mhd/cgl_lf_strict_admissibility=false"
R14_APPROVED_VARIANT = "finite_limiter_hard_bound_diagnostic_nonfatal"
R14_APPROVED_OVERRIDES = [R14_STRICT_OVERRIDE]
R14_APPROVED_MODEL_CHOICES = {
    "backup_limiters": "false",
    "cgl_lf_strict_admissibility": "false",
    "firehose_limiter": "true",
    "limiter_nu_coll": "20.0",
    "mirror_limiter": "true",
}
MATRIX_RELATIVE = Path("inputs/cgl_lf_paper/mks24_stage_i_manifest.json")
RUNS_RELATIVE = Path("runs/mks24-stage-i-fast/E03-forcing-policy")
PASSIVE_RSOLVER_RELATIVE = Path("src/mhd/rsolvers/hlle_cgl.hpp")
EOS_RELATIVE = Path("src/eos/eos.hpp")
COMPATIBILITY_COMMIT = "0fad9bc0282f9858d5c22e3ce7c87227b404ecde"
COMPATIBILITY_TEST_RELATIVE = Path(
    "tst/test_suite/cgl/test_cgl_passive_fast_path_cpu.py"
)
COMPATIBILITY_TEST_SHA256 = (
    "d6a61d09d02057e145343bb47e3adc383049afb4ef6914339d11cd7a9531ff46"
)
COMPATIBILITY_FIXTURE_RELATIVE = Path(
    "src/pgen/unit_tests/cgl_passive_fast_path_test.cpp"
)
COMPATIBILITY_FIXTURE_SHA256 = (
    "9223b5cdceb4e383efcc5d35ed08468ac4e1790f5bb9388ccb76efd194985931"
)
LEGACY_EOS_TERM = b"+ 12.0*ppar*pperp*bhatx2*bhatx2"
CORRECTED_EOS_TERM = b"+ 12.0*ppar*ppar*bhatx2*bhatx2"
PASSIVE_BRANCH = b"""if (eos.passive) {
      cl = eos.IdealMHDFastSpeed(wl_idn, bxi, wl_iby, wl_ibz);
      cr = eos.IdealMHDFastSpeed(wr_idn, bxi, wr_iby, wr_ibz);
    } else {
      cl = eos.IdealMHDFastSpeed(wl_idn, wl_ipr, wl_ipp, bxi, wl_iby, wl_ibz, bfloor);
      cr = eos.IdealMHDFastSpeed(wr_idn, wr_ipr, wr_ipp, bxi, wr_iby, wr_ibz, bfloor);
    }"""


def load_module(path: Path, name: str) -> object:
    """Import one local utility by exact path."""

    existing = sys.modules.get(name)
    if existing is not None:
        return existing
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import utility: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


report = load_module(REPORTER_PATH, "_cgl_lf_stage_i_corrected_composite_reporter")
identity_tool = load_module(
    IDENTITY_TOOL_PATH, "_cgl_lf_stage_i_corrected_composite_identity"
)


class CompositeReportError(RuntimeError):
    """A corrected/legacy composite evidence contract failed."""


@dataclass(frozen=True)
class ExecutionAuthority:
    """Immutable execution authority for one evidence class."""

    campaign_id: str
    evidence_class: str
    root: Path
    source: Path
    source_revision: str
    executable: Path
    executable_sha256: str
    matrix: Path
    matrix_sha256: str
    run_relatives: tuple[Path, ...] = (RUNS_RELATIVE,)

    @property
    def run_root(self) -> Path:
        return self.run_roots[0]

    @property
    def run_roots(self) -> tuple[Path, ...]:
        return tuple(self.root / relative for relative in self.run_relatives)


@dataclass(frozen=True)
class PassiveCompatibility:
    """Source and input identities that justify legacy passive control reuse."""

    rsolver_sha256: str
    legacy_eos_sha256: str
    corrected_eos_sha256: str
    input_sha256: tuple[tuple[str, str], ...]

    def input_hashes(self) -> dict[str, str]:
        return dict(self.input_sha256)


@dataclass(frozen=True)
class CompositeConfig:
    """Complete authority configuration for one composite assembly."""

    corrected: ExecutionAuthority
    legacy: ExecutionAuthority
    corrected_identity: Path
    corrected_identity_sha256: str
    compatibility: PassiveCompatibility


DEFAULT_CORRECTED_ROOT = Path(
    "/lustre/orion/ast207/proj-shared/dfielding/CGL/campaigns/"
    "mks24-stage-i-eos-fastdisc-ppar2-corrected-v1"
)
DEFAULT_LEGACY_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/CGL")
DEFAULT_CORRECTED_SOURCE = Path(
    "/autofs/nccs-svm1_home2/dfielding/athenak-cgl-e03-corrected-0c406312f"
)
DEFAULT_LEGACY_SOURCE = Path(
    "/autofs/nccs-svm1_home2/dfielding/athenak-cgl-e03-9e075422"
)
DEFAULT_CORRECTED_EXECUTABLE = Path(
    "/lustre/orion/ast207/proj-shared/dfielding/CGL/build/"
    "frontier-hip-0c406312fa35-cpe25.09-cce20-rocm6.4.2/src/athena"
)
DEFAULT_LEGACY_EXECUTABLE = Path(
    "/lustre/orion/ast207/proj-shared/dfielding/CGL/build/"
    "frontier-hip-9e07542281e4-cpe25.09-cce20-rocm6.4.2/src/athena"
)
DEFAULT_OUTPUT = DEFAULT_CORRECTED_ROOT / (
    "analysis/mks24-stage-i-composite/E03-forcing-policy/R02-R17-final"
)
DEFAULT_CONFIG = CompositeConfig(
    corrected=ExecutionAuthority(
        campaign_id="mks24-stage-i-eos-fastdisc-ppar2-corrected-v1",
        evidence_class="corrected_active_production",
        root=DEFAULT_CORRECTED_ROOT,
        source=DEFAULT_CORRECTED_SOURCE,
        source_revision="0c406312fa35d5e1c7041d80333b0ea24b0127ae",
        executable=DEFAULT_CORRECTED_EXECUTABLE,
        executable_sha256=(
            "0f032379c4d829fc86353cc734a716f6ef6d6eda8951a5b79a3a97ad671e9b4b"
        ),
        matrix=DEFAULT_CORRECTED_SOURCE / MATRIX_RELATIVE,
        matrix_sha256=(
            "bf31b88b985d1ad4ffe823108dd7c1132bdfa4d5e4a6abde51f66bb7778415c9"
        ),
    ),
    legacy=ExecutionAuthority(
        campaign_id="mks24-stage-i-qualified-legacy-e03",
        evidence_class="authenticated_legacy_passive_control",
        root=DEFAULT_LEGACY_ROOT,
        source=DEFAULT_LEGACY_SOURCE,
        source_revision="9e07542281e4e6d125582f253df3ad2e3b8b154d",
        executable=DEFAULT_LEGACY_EXECUTABLE,
        executable_sha256=(
            "68f243f9204df388b24365ae65a567f6f567dbe422a6d7a43b9fb4a499ef118c"
        ),
        matrix=DEFAULT_LEGACY_SOURCE / MATRIX_RELATIVE,
        matrix_sha256=(
            "bf31b88b985d1ad4ffe823108dd7c1132bdfa4d5e4a6abde51f66bb7778415c9"
        ),
        run_relatives=(RUNS_RELATIVE, report.RACE_RUNS_RELATIVE),
    ),
    corrected_identity=DEFAULT_CORRECTED_ROOT / "campaign-identity.json",
    corrected_identity_sha256=(
        "fad99d9651e6da155bb5dfae1b81a8675cce232fb17e31a36e9be36d62810643"
    ),
    compatibility=PassiveCompatibility(
        rsolver_sha256=(
            "07245ddc3a318cec061ad3fd096c0e9092a15bd3a709b8c9a7193dd6565a9d82"
        ),
        legacy_eos_sha256=(
            "3925911d48ae8d59217977648ddd5314c49ae8d94f8869d5489d53355735aea7"
        ),
        corrected_eos_sha256=(
            "4774042fdb5b8e4835954a077e1d5ed38dc8c9a65b3425b8b09527f5e2909ac6"
        ),
        input_sha256=(
            (
                "R06",
                "c6e038c198b23bf20a7cf1e2a90fa82e83e5ec544492ddc64e1cf2bb535a38ae",
            ),
            (
                "R07",
                "72ed6a3b38342d00855f34a7c7eb67102b44cd22b6c12ff2448ddb6f141e5145",
            ),
            (
                "R08",
                "3535aca5e3d47dc383e353cff262d79c56cf3c083a3ac6024d04ddaa5d353780",
            ),
            (
                "R09",
                "2ec05b247d917259c8306911cfd31d6be4e1c2772d6cdea544db883edfc35f2c",
            ),
        ),
    ),
)


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of one regular file."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_json(value: object) -> bytes:
    """Return deterministic readable JSON bytes."""

    return (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def write_json(path: Path, value: object) -> None:
    """Atomically write deterministic JSON."""

    path.parent.mkdir(parents=True, exist_ok=True)
    staged: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb", dir=path.parent, prefix=f".{path.name}.", delete=False
        ) as stream:
            staged = Path(stream.name)
            stream.write(stable_json(value))
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(staged, path)
        staged = None
    finally:
        if staged is not None:
            staged.unlink(missing_ok=True)


def load_json(path: Path, label: str) -> dict[str, object]:
    """Load one JSON object."""

    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise CompositeReportError(f"cannot load {label}: {path}") from error
    if not isinstance(value, dict):
        raise CompositeReportError(f"{label} must be a JSON object: {path}")
    return value


def canonical_existing(path: Path, label: str, *, directory: bool) -> Path:
    """Resolve and type-check one required path."""

    try:
        resolved = path.expanduser().resolve(strict=True)
    except OSError as error:
        raise CompositeReportError(f"{label} does not resolve: {path}") from error
    if directory and not resolved.is_dir():
        raise CompositeReportError(f"{label} is not a directory: {resolved}")
    if not directory and not resolved.is_file():
        raise CompositeReportError(f"{label} is not a regular file: {resolved}")
    return resolved


def require_beneath(path: Path | str, parent: Path, label: str) -> Path:
    """Require one path to resolve below an exact authority root."""

    try:
        resolved = Path(path).expanduser().resolve(strict=True)
        root = parent.expanduser().resolve(strict=True)
        relative = resolved.relative_to(root)
    except (OSError, ValueError) as error:
        raise CompositeReportError(
            f"{label} escapes or does not resolve below {parent}: {path}"
        ) from error
    if not relative.parts:
        raise CompositeReportError(f"{label} must be below {root}: {resolved}")
    return resolved


def require_hash(path: Path, expected: str, label: str) -> dict[str, object]:
    """Authenticate one exact regular-file digest."""

    resolved = canonical_existing(path, label, directory=False)
    observed = sha256(resolved)
    if observed != expected:
        raise CompositeReportError(
            f"{label} SHA-256 mismatch: observed {observed}, expected {expected}"
        )
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": observed,
    }


def git_revision(source: Path) -> str:
    """Return the exact Git revision for one source tree."""

    try:
        completed = subprocess.run(
            ["/usr/bin/git", "-C", str(source), "rev-parse", "--verify", "HEAD"],
            check=True,
            text=True,
            capture_output=True,
        )
        dirty = subprocess.run(
            [
                "/usr/bin/git",
                "-C",
                str(source),
                "status",
                "--porcelain",
                "--untracked-files=no",
            ],
            check=True,
            text=True,
            capture_output=True,
        )
    except subprocess.SubprocessError as error:
        raise CompositeReportError(f"cannot authenticate source revision: {source}") from error
    if dirty.stdout.strip():
        raise CompositeReportError(f"source has tracked modifications: {source}")
    return completed.stdout.strip()


def committed_payload(commit: str, relative: Path, expected_sha256: str) -> dict[str, object]:
    """Authenticate one exact file payload from a Git commit object."""

    try:
        payload = subprocess.run(
            ["/usr/bin/git", "-C", str(REPO_ROOT), "show", f"{commit}:{relative}"],
            check=True,
            capture_output=True,
        ).stdout
    except subprocess.SubprocessError as error:
        raise CompositeReportError(
            f"cannot read compatibility regression artifact {relative} from {commit}"
        ) from error
    observed = hashlib.sha256(payload).hexdigest()
    if observed != expected_sha256:
        raise CompositeReportError(
            f"compatibility regression artifact differs at {commit}:{relative}"
        )
    return {
        "commit": commit,
        "git_path": relative.as_posix(),
        "size_bytes": len(payload),
        "sha256": observed,
    }


def validate_compatibility_regression() -> dict[str, object]:
    """Authenticate the reviewed passive-path compatibility disposition."""

    try:
        resolved = subprocess.run(
            [
                "/usr/bin/git",
                "-C",
                str(REPO_ROOT),
                "rev-parse",
                "--verify",
                f"{COMPATIBILITY_COMMIT}^{{commit}}",
            ],
            check=True,
            text=True,
            capture_output=True,
        ).stdout.strip()
        ancestry = subprocess.run(
            [
                "/usr/bin/git",
                "-C",
                str(REPO_ROOT),
                "merge-base",
                "--is-ancestor",
                COMPATIBILITY_COMMIT,
                "HEAD",
            ],
            check=False,
        )
    except subprocess.SubprocessError as error:
        raise CompositeReportError(
            f"cannot authenticate compatibility regression commit {COMPATIBILITY_COMMIT}"
        ) from error
    if resolved != COMPATIBILITY_COMMIT:
        raise CompositeReportError(
            f"compatibility regression commit differs: {resolved}"
        )
    if ancestry.returncode != 0:
        raise CompositeReportError(
            f"compatibility regression commit is not an ancestor of HEAD: {resolved}"
        )
    return {
        "disposition": "compatible_for_legacy_passive_control_reuse",
        "basis": "reviewed_standalone_and_integration_regression",
        "commit": resolved,
        "test": committed_payload(
            resolved, COMPATIBILITY_TEST_RELATIVE, COMPATIBILITY_TEST_SHA256
        ),
        "fixture": committed_payload(
            resolved, COMPATIBILITY_FIXTURE_RELATIVE, COMPATIBILITY_FIXTURE_SHA256
        ),
        "command": (
            "python3 -m pytest -q "
            "tst/test_suite/cgl/test_cgl_passive_fast_path_cpu.py"
        ),
        "reviewed_result": {
            "result": "pass",
            "summary": "1 passed",
            "elapsed_seconds": 233.45,
        },
        "established_behavior": [
            "both pressure controls are explicitly inside the hard-bound admissible region",
            "passive signal speed uses the isothermal-MHD overload",
            "passive timesteps and flow fields are pressure independent",
            "passive physical and LLF flow fluxes are pressure independent",
            "active CGL control responds to pressure changes",
        ],
    }


def authority_record(authority: ExecutionAuthority) -> dict[str, object]:
    """Authenticate and describe one execution authority."""

    root = canonical_existing(authority.root, "campaign root", directory=True)
    source = canonical_existing(authority.source, "source root", directory=True)
    run_roots = [
        canonical_existing(path, "run root", directory=True)
        for path in authority.run_roots
    ]
    for run_root in run_roots:
        if require_beneath(run_root, root, "run root") != run_root:
            raise AssertionError("unreachable run-root validation state")
    if len(run_roots) != len(set(run_roots)):
        raise CompositeReportError(f"{authority.evidence_class} run roots are not unique")
    observed_revision = git_revision(source)
    if observed_revision != authority.source_revision:
        raise CompositeReportError(
            f"{authority.evidence_class} source revision mismatch: "
            f"{observed_revision} != {authority.source_revision}"
        )
    return {
        "campaign_id": authority.campaign_id,
        "evidence_class": authority.evidence_class,
        "campaign_root": str(root),
        "run_roots": [str(run_root) for run_root in run_roots],
        "source": {
            "path": str(source),
            "revision": observed_revision,
        },
        "executable": require_hash(
            authority.executable,
            authority.executable_sha256,
            f"{authority.evidence_class} executable",
        ),
        "matrix": require_hash(
            authority.matrix,
            authority.matrix_sha256,
            f"{authority.evidence_class} matrix",
        ),
    }


def matrix_cases(path: Path) -> dict[str, dict[str, object]]:
    """Return exact Stage I matrix cases."""

    matrix = load_json(path, "Stage I matrix")
    values = matrix.get("cases")
    if not isinstance(values, list):
        raise CompositeReportError("Stage I matrix lacks a cases list")
    result = {
        str(value["id"]): value
        for value in values
        if isinstance(value, dict) and isinstance(value.get("id"), str)
    }
    if set(result) != set(ALL_CASES):
        raise CompositeReportError(
            f"Stage I matrix case coverage differs: {sorted(result)}"
        )
    return result


def expected_input(
    authority: ExecutionAuthority, case: dict[str, object], case_id: str
) -> Path:
    """Return an authenticated matrix-selected input below an authority source."""

    relative = case.get("input")
    if not isinstance(relative, str) or not relative:
        raise CompositeReportError(f"{case_id} matrix input is malformed")
    path = require_beneath(
        authority.source / relative,
        authority.source,
        f"{case_id} {authority.evidence_class} input",
    )
    if not path.is_file():
        raise CompositeReportError(f"{case_id} input is not a regular file: {path}")
    return path


def validate_passive_compatibility(
    config: CompositeConfig, cases: dict[str, dict[str, object]]
) -> dict[str, object]:
    """Authenticate the source-level contract for reusing legacy passive controls."""

    corrected = config.corrected
    legacy = config.legacy
    contract = config.compatibility
    if corrected.matrix_sha256 != legacy.matrix_sha256:
        raise CompositeReportError("corrected and legacy matrix identities differ")

    rsolver_bindings = {}
    rsolver_payloads = []
    for label, authority in (("corrected", corrected), ("legacy", legacy)):
        path = authority.source / PASSIVE_RSOLVER_RELATIVE
        rsolver_bindings[label] = require_hash(
            path, contract.rsolver_sha256, f"{label} passive Riemann solver"
        )
        rsolver_payloads.append(path.read_bytes())
    if rsolver_payloads[0] != rsolver_payloads[1]:
        raise CompositeReportError("corrected and legacy passive Riemann solvers differ")
    if rsolver_payloads[0].count(PASSIVE_BRANCH) != 1:
        raise CompositeReportError("authenticated passive signal-speed branch differs")

    legacy_eos = legacy.source / EOS_RELATIVE
    corrected_eos = corrected.source / EOS_RELATIVE
    eos_bindings = {
        "legacy": require_hash(
            legacy_eos, contract.legacy_eos_sha256, "legacy EOS source"
        ),
        "corrected": require_hash(
            corrected_eos, contract.corrected_eos_sha256, "corrected EOS source"
        ),
    }
    legacy_payload = legacy_eos.read_bytes()
    corrected_payload = corrected_eos.read_bytes()
    if legacy_payload.count(LEGACY_EOS_TERM) != 1:
        raise CompositeReportError("legacy EOS does not contain the one qualified term")
    if corrected_payload.count(CORRECTED_EOS_TERM) != 1:
        raise CompositeReportError("corrected EOS does not contain the one corrected term")
    if legacy_payload.replace(LEGACY_EOS_TERM, CORRECTED_EOS_TERM, 1) != corrected_payload:
        raise CompositeReportError(
            "legacy/corrected EOS sources differ beyond the active fast-speed term"
        )

    expected_hashes = contract.input_hashes()
    if set(expected_hashes) != set(PASSIVE_CONTROLS):
        raise CompositeReportError("passive compatibility input inventory differs")
    inputs: dict[str, object] = {}
    for case_id in PASSIVE_CONTROLS:
        case = cases[case_id]
        legacy_input = expected_input(legacy, case, case_id)
        corrected_input = expected_input(corrected, case, case_id)
        expected_sha = expected_hashes[case_id]
        legacy_binding = require_hash(
            legacy_input, expected_sha, f"{case_id} legacy passive input"
        )
        corrected_binding = require_hash(
            corrected_input, expected_sha, f"{case_id} corrected passive input"
        )
        if legacy_input.read_bytes() != corrected_input.read_bytes():
            raise CompositeReportError(f"{case_id} corrected and legacy inputs differ")
        choices = report.model_choices_for_input(legacy_input)
        if choices.get("passive_delta") != "true":
            raise CompositeReportError(f"{case_id} is not an explicit passive-Delta input")
        inputs[case_id] = {
            "legacy": legacy_binding,
            "corrected": corrected_binding,
            "passive_delta": choices["passive_delta"],
        }
    return {
        "contract_id": "legacy-passive-controls-regression-0fad9bc02-v1",
        "result": "pass",
        "scope": list(PASSIVE_CONTROLS),
        "claim": (
            "Legacy passive controls retain their legacy execution identity and are "
            "compatible for reuse under the reviewed passive-path regression."
        ),
        "disposition_basis": validate_compatibility_regression(),
        "matrix_sha256": corrected.matrix_sha256,
        "supporting_evidence": {
            "passive_rsolver": rsolver_bindings,
            "active_only_eos_change": eos_bindings,
            "inputs": inputs,
        },
    }


def selected_corrected_paths(config: CompositeConfig) -> list[Path]:
    """Return all corrected active case roots selected for composite assembly."""

    return [config.corrected.run_root / case_id for case_id in ACTIVE_CASES]


def selected_authority_run_root(
    authority: ExecutionAuthority, path: Path | str, label: str
) -> tuple[Path, Path]:
    """Return an exact selected path and its one authorized run root."""

    try:
        resolved = Path(path).expanduser().resolve(strict=True)
    except OSError as error:
        raise CompositeReportError(f"{label} does not resolve: {path}") from error
    matches = []
    for run_root in authority.run_roots:
        root = run_root.resolve(strict=True)
        try:
            relative = resolved.relative_to(root)
        except ValueError:
            continue
        if relative.parts:
            matches.append(root)
    if len(matches) != 1:
        raise CompositeReportError(
            f"{label} matched {len(matches)} authorized run roots: {resolved}"
        )
    return resolved, matches[0]


def validate_authorities(config: CompositeConfig) -> dict[str, object]:
    """Authenticate both authorities, corrected identity, and passive compatibility."""

    corrected = authority_record(config.corrected)
    legacy = authority_record(config.legacy)
    if config.corrected.root.resolve(strict=True) == config.legacy.root.resolve(strict=True):
        raise CompositeReportError("corrected and legacy campaign roots must differ")
    identity_binding = require_hash(
        config.corrected_identity,
        config.corrected_identity_sha256,
        "corrected campaign identity",
    )
    try:
        identity_result = identity_tool.validate_identity(
            config.corrected_identity, selected_corrected_paths(config)
        )
    except identity_tool.IdentityError as error:
        raise CompositeReportError(f"corrected identity validation failed: {error}") from error
    if identity_result.get("campaign_root") != str(config.corrected.root.resolve(strict=True)):
        raise CompositeReportError("corrected identity campaign root differs")
    if identity_result.get("source_revision") != config.corrected.source_revision:
        raise CompositeReportError("corrected identity source revision differs")
    artifacts = identity_result.get("artifacts")
    if not isinstance(artifacts, dict):
        raise CompositeReportError("corrected identity artifacts are malformed")
    expected_identity_artifacts = {
        "executable": config.corrected.executable_sha256,
        "matrix": config.corrected.matrix_sha256,
        "eos": config.compatibility.corrected_eos_sha256,
    }
    for name, expected_sha in expected_identity_artifacts.items():
        binding = artifacts.get(name)
        if not isinstance(binding, dict) or binding.get("sha256") != expected_sha:
            raise CompositeReportError(f"corrected identity {name} artifact differs")

    cases = matrix_cases(config.corrected.matrix)
    compatibility = validate_passive_compatibility(config, cases)
    return {
        "authorities": {
            config.corrected.evidence_class: corrected,
            config.legacy.evidence_class: legacy,
        },
        "corrected_identity": identity_binding,
        "corrected_identity_validation": identity_result,
        "passive_compatibility": compatibility,
        "cases": cases,
    }


def authority_for_case(config: CompositeConfig, case_id: str) -> ExecutionAuthority:
    """Return the only permitted execution authority for one Stage I case."""

    if case_id in ACTIVE_CASES:
        return config.corrected
    if case_id in PASSIVE_CONTROLS:
        return config.legacy
    raise CompositeReportError(f"unsupported Stage I case: {case_id}")


def require_manifest_value(
    manifest: dict[str, object], key: str, expected: object, label: str
) -> None:
    """Require one exact manifest field."""

    if manifest.get(key) != expected:
        raise CompositeReportError(
            f"{label} {key} mismatch: observed {manifest.get(key)!r}, expected {expected!r}"
        )


def validate_segment_manifest(
    config: CompositeConfig,
    authority: ExecutionAuthority,
    case_id: str,
    case: dict[str, object],
    lineage: dict[str, object],
) -> dict[str, object]:
    """Authenticate one selected segment without relabeling its executable or source."""

    if lineage.get("kind") != "fast":
        raise CompositeReportError(f"{case_id} selected lineage contains non-fast evidence")
    manifest_binding = lineage.get("manifest")
    if not isinstance(manifest_binding, dict):
        raise CompositeReportError(f"{case_id} selected lineage lacks a manifest binding")
    manifest_path_value = manifest_binding.get("path")
    if not isinstance(manifest_path_value, str):
        raise CompositeReportError(f"{case_id} selected manifest path is malformed")
    manifest_path, selected_run_root = selected_authority_run_root(
        authority,
        manifest_path_value,
        f"{case_id} selected manifest",
    )
    if sha256(manifest_path) != manifest_binding.get("sha256"):
        raise CompositeReportError(f"{case_id} selected manifest binding is stale")
    manifest = load_json(manifest_path, f"{case_id} selected manifest")

    input_path = expected_input(authority, case, case_id)
    exact = {
        "case_id": case_id,
        "case_name": case["name"],
        "root": str(authority.root.resolve(strict=True)),
        "input": str(input_path),
        "input_sha256": sha256(input_path),
        "matrix_sha256": authority.matrix_sha256,
        "executable": str(authority.executable.resolve(strict=True)),
        "executable_sha256": authority.executable_sha256,
    }
    for key, expected in exact.items():
        require_manifest_value(manifest, key, expected, case_id)

    run_dir, run_root = selected_authority_run_root(
        authority, str(manifest.get("run_dir", "")), f"{case_id} run directory"
    )
    if run_root != selected_run_root:
        raise CompositeReportError(f"{case_id} manifest and run directory roots differ")
    if run_dir.parent.name != case_id:
        raise CompositeReportError(f"{case_id} run directory is outside its case directory")
    output_dir = require_beneath(
        str(manifest.get("output_dir", "")),
        run_dir,
        f"{case_id} output directory",
    )
    if output_dir != run_dir / "output":
        raise CompositeReportError(f"{case_id} output is not the selected run output")
    if manifest_path != run_dir / "manifest/fast_run.json":
        raise CompositeReportError(f"{case_id} manifest is outside the selected run")
    lineage_exact = {
        "segment_dir": str(run_dir),
        "source_root": str(selected_run_root),
        "input": str(input_path),
        "input_sha256": sha256(input_path),
        "matrix_sha256": authority.matrix_sha256,
        "executable": str(authority.executable.resolve(strict=True)),
        "executable_sha256": authority.executable_sha256,
    }
    for key, expected in lineage_exact.items():
        if lineage.get(key) != expected:
            raise CompositeReportError(
                f"{case_id} retained lineage {key} differs from its manifest"
            )

    if case_id in ACTIVE_CASES:
        active_exact = {
            "campaign_id": config.corrected.campaign_id,
            "campaign_root": str(config.corrected.root.resolve(strict=True)),
            "campaign_identity": str(config.corrected_identity.resolve(strict=True)),
            "campaign_identity_sha256": config.corrected_identity_sha256,
            "source": str(config.corrected.source.resolve(strict=True)),
            "source_revision": config.corrected.source_revision,
            "legacy_restart_permitted": False,
        }
        for key, expected in active_exact.items():
            require_manifest_value(manifest, key, expected, case_id)
    else:
        if manifest.get("source") not in (None, str(config.legacy.source.resolve(strict=True))):
            raise CompositeReportError(f"{case_id} legacy source was relabeled")
        if manifest.get("source_revision") not in (None, config.legacy.source_revision):
            raise CompositeReportError(f"{case_id} legacy source revision was relabeled")
    return manifest


def case_execution_record(
    authority: ExecutionAuthority,
    case_id: str,
    case: dict[str, object],
    compatibility: dict[str, object],
) -> dict[str, object]:
    """Build the explicit per-case execution/evidence authority record."""

    record = {
        "case_id": case_id,
        "evidence_class": authority.evidence_class,
        "campaign_id": authority.campaign_id,
        "campaign_root": str(authority.root.resolve(strict=True)),
        "authorized_run_roots": [
            str(run_root.resolve(strict=True)) for run_root in authority.run_roots
        ],
        "source": {
            "path": str(authority.source.resolve(strict=True)),
            "revision": authority.source_revision,
        },
        "executable": require_hash(
            authority.executable,
            authority.executable_sha256,
            f"{case_id} execution executable",
        ),
        "matrix": require_hash(
            authority.matrix, authority.matrix_sha256, f"{case_id} execution matrix"
        ),
        "input": report.artifact_binding(expected_input(authority, case, case_id)),
    }
    if case_id in PASSIVE_CONTROLS:
        record["compatibility_contract_id"] = compatibility["contract_id"]
    else:
        record["compatibility_contract_id"] = None
    return record


def validate_selected_attempt_sequences(
    manifests: list[dict[str, object]], case_id: str
) -> list[int]:
    """Validate ordered attempt IDs without treating failed attempts as lineage."""

    sequences = [manifest.get("sequence") for manifest in manifests]
    if (
        not sequences
        or sequences[0] != 0
        or any(
            not isinstance(sequence, int) or isinstance(sequence, bool)
            for sequence in sequences
        )
        or any(
            current <= previous
            for previous, current in zip(sequences, sequences[1:])
        )
    ):
        raise CompositeReportError(
            f"{case_id} selected attempt sequence is not strictly increasing: "
            f"{sequences}"
        )
    return [int(sequence) for sequence in sequences]


def finite_float(value: object, label: str) -> float:
    """Return one finite floating-point value."""

    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise CompositeReportError(f"{label} is not numeric: {value!r}") from error
    if not math.isfinite(result):
        raise CompositeReportError(f"{label} is not finite: {value!r}")
    return result


def history_final_time(path: Path, label: str) -> float:
    """Authenticate one retained history and return its terminal time."""

    try:
        history, warnings = report.read_history_source(path)
    except (OSError, report.ReportError, ValueError) as error:
        raise CompositeReportError(f"{label} is unavailable: {path}") from error
    if warnings:
        raise CompositeReportError(f"{label} is not stable: {warnings}")
    rows = history.get("rows")
    if not isinstance(rows, list) or not rows or not isinstance(rows[-1], list):
        raise CompositeReportError(f"{label} has no retained rows: {path}")
    return finite_float(rows[-1][0], f"{label} terminal time")


def slurm_log_for_manifest(
    authority: ExecutionAuthority,
    case_id: str,
    manifest: dict[str, object],
) -> Path:
    """Resolve the exact immutable Slurm log declared by one failed attempt."""

    sequence = manifest.get("sequence")
    job_id = manifest.get("job_id")
    expected_pattern = authority.root / "logs/slurm-fast/%x.%j.log"
    if (
        not isinstance(sequence, int)
        or isinstance(sequence, bool)
        or not isinstance(job_id, str)
        or re.fullmatch(r"[1-9][0-9]*", job_id) is None
        or manifest.get("slurm_log") != str(expected_pattern)
    ):
        raise CompositeReportError(
            f"{case_id} failed attempt lacks exact Slurm log metadata"
        )
    job_name = f"cglc_{case_id}_s{sequence:03d}"
    path = Path(
        str(expected_pattern).replace("%x", job_name).replace("%j", job_id)
    )
    resolved = canonical_existing(path, f"{case_id} failed-attempt Slurm log", directory=False)
    require_beneath(
        resolved,
        authority.root / "logs/slurm-fast",
        f"{case_id} failed-attempt Slurm log",
    )
    return resolved


def synthetic_lineage_for_manifest(
    authority: ExecutionAuthority, manifest_path: Path
) -> dict[str, object]:
    """Build the minimal lineage projection used to authenticate a source attempt."""

    manifest_path, run_root = selected_authority_run_root(
        authority, manifest_path, "source attempt manifest"
    )
    manifest = load_json(manifest_path, "source attempt manifest")
    run_dir = manifest_path.parent.parent
    return {
        "kind": "fast",
        "segment_dir": str(run_dir),
        "source_root": str(run_root),
        "input": manifest.get("input"),
        "input_sha256": manifest.get("input_sha256"),
        "matrix_sha256": manifest.get("matrix_sha256"),
        "executable": manifest.get("executable"),
        "executable_sha256": manifest.get("executable_sha256"),
        "manifest": report.artifact_binding(manifest_path),
    }


def authenticate_r14_failed_attempt(
    config: CompositeConfig,
    case: dict[str, object],
    manifest_path: Path,
    selected_variant: object,
    selected_overrides: list[str],
    parent_restart: Path,
    parent_restart_sha256: str,
) -> dict[str, object] | None:
    """Authenticate one source attempt replaying the selected R14 checkpoint."""

    authority = config.corrected
    try:
        preview = load_json(manifest_path, "R14 source attempt manifest")
        if (
            preview.get("restart") != str(parent_restart)
            or preview.get("restart_sha256") != parent_restart_sha256
        ):
            return None
        lineage = synthetic_lineage_for_manifest(authority, manifest_path)
        manifest = validate_segment_manifest(
            config, authority, R14_FAILURE_CASE, case, lineage
        )
        sequence = manifest.get("sequence")
        run_dir = Path(str(manifest["run_dir"])).resolve(strict=True)
        if (
            not isinstance(sequence, int)
            or isinstance(sequence, bool)
            or sequence <= 0
            or manifest.get("variant") != R14_APPROVED_VARIANT
            or selected_variant != R14_APPROVED_VARIANT
            or manifest.get("command_line_overrides") != R14_APPROVED_OVERRIDES
            or selected_overrides != R14_APPROVED_OVERRIDES
            or manifest.get("strict_admissibility") is not False
            or manifest.get("continuation_policy")
            != "finite_progress_complete_terminal_products"
            or manifest.get("runtime_segmentation_changes_physics") is not False
        ):
            return None

        exit_path = run_dir / "manifest/run_exit_code"
        exit_binding = report.artifact_binding(
            canonical_existing(
                exit_path, "R14 failed-attempt exit artifact", directory=False
            )
        )
        exit_code = int(exit_path.read_text(encoding="utf-8").strip())
        if exit_code == 0:
            return None

        output = run_dir / "output"
        mhd_path, user_path = report.history_paths(output)
        if mhd_path is None or user_path is None:
            return None
        mhd_final = history_final_time(mhd_path, "R14 failed-attempt MHD history")
        user_final = history_final_time(user_path, "R14 failed-attempt user history")
        if not math.isclose(
            mhd_final,
            user_final,
            rel_tol=0.0,
            abs_tol=R14_FAILURE_TIME_TOLERANCE,
        ):
            return None

        log_path = slurm_log_for_manifest(authority, R14_FAILURE_CASE, manifest)
        log_text = log_path.read_text(encoding="utf-8", errors="replace")
        signature_count = log_text.count(R14_FATAL_SIGNATURE)
        if signature_count <= 0:
            return None
        return {
            "sequence": sequence,
            "segment": str(run_dir),
            "job_id": manifest["job_id"],
            "run_exit_code": exit_code,
            "failure_time": mhd_final,
            "fatal_signature_count": signature_count,
            "mhd_history_sha256": report.artifact_binding(mhd_path)["sha256"],
            "user_history_sha256": report.artifact_binding(user_path)["sha256"],
            "provenance": {
                "manifest": report.artifact_binding(manifest_path),
                "run_exit_code": exit_binding,
                "mhd_history": report.artifact_binding(mhd_path),
                "user_history": report.artifact_binding(user_path),
                "slurm_log": report.artifact_binding(log_path),
            },
        }
    except (
        CompositeReportError,
        OSError,
        UnicodeError,
        ValueError,
        report.ReportError,
    ):
        return None


def retained_r14_history(
    record: dict[str, object], failure_time: float
) -> dict[str, object]:
    """Authenticate the merged R14 history retained by the composite report."""

    histories = record.get("histories")
    if not isinstance(histories, dict):
        raise CompositeReportError("R14 retained history inventory is missing")
    retained: dict[str, object] = {}
    for kind in ("mhd", "user"):
        item = histories.get(kind)
        if (
            not isinstance(item, dict)
            or item.get("available") is not True
            or item.get("errors") not in (None, [])
            or not isinstance(item.get("path"), str)
            or not isinstance(item.get("binding"), dict)
        ):
            raise CompositeReportError(f"R14 retained {kind} history is unavailable")
        path = canonical_existing(
            Path(str(item["path"])), f"R14 retained {kind} history", directory=False
        )
        if item["binding"] != report.artifact_binding(path):
            raise CompositeReportError(
                f"R14 retained {kind} history binding is stale"
            )
        final = history_final_time(path, f"R14 retained {kind} history")
        if not math.isclose(
            final,
            failure_time,
            rel_tol=0.0,
            abs_tol=R14_FAILURE_TIME_TOLERANCE,
        ):
            raise CompositeReportError(
                f"R14 retained {kind} history does not reach the failure time"
            )
        retained[kind] = report.artifact_binding(path)
    return retained


def derive_r14_terminal_disposition(
    config: CompositeConfig,
    case: dict[str, object],
    record: dict[str, object],
    manifests: list[dict[str, object]],
) -> dict[str, object]:
    """Derive the only admissible R14 partial disposition from source attempts."""

    lineage = record.get("lineage")
    if (
        record.get("status") != "failed_partial"
        or not isinstance(lineage, list)
        or len(lineage) < 2
        or len(manifests) != len(lineage)
    ):
        raise CompositeReportError("R14 is not an authenticated failed partial lineage")
    terminal = lineage[-1]
    parent = lineage[-2]
    if (
        not isinstance(terminal, dict)
        or not isinstance(parent, dict)
        or terminal.get("state") != "failed"
        or terminal.get("kind") != "fast"
        or parent.get("kind") != "fast"
    ):
        raise CompositeReportError("R14 selected terminal is not a failed fast attempt")

    selected_manifest = manifests[-1]
    selected_overrides = selected_manifest.get("command_line_overrides")
    if not isinstance(selected_overrides, list) or not all(
        isinstance(value, str) for value in selected_overrides
    ):
        raise CompositeReportError("R14 selected overrides are malformed")
    expected_model = report.model_choices_for_input(
        expected_input(config.corrected, case, R14_FAILURE_CASE),
        selected_overrides,
    )
    if (
        expected_model.get("cgl_lf_strict_admissibility") != "false"
        or record.get("model_choices") != expected_model
        or any(
            expected_model.get(key) != expected
            for key, expected in R14_APPROVED_MODEL_CHOICES.items()
        )
        or selected_manifest.get("variant") != R14_APPROVED_VARIANT
        or selected_overrides != R14_APPROVED_OVERRIDES
        or selected_manifest.get("strict_admissibility") is not False
        or selected_manifest.get("continuation_policy")
        != "finite_progress_complete_terminal_products"
        or selected_manifest.get("runtime_segmentation_changes_physics") is not False
    ):
        raise CompositeReportError(
            "R14 failed-partial disposition differs from the approved physics contract"
        )

    parent_segment = Path(str(parent.get("segment_dir", ""))).resolve(strict=True)
    parent_restart_value = selected_manifest.get("restart")
    parent_restart_sha256 = selected_manifest.get("restart_sha256")
    if (
        not isinstance(parent_restart_value, str)
        or not isinstance(parent_restart_sha256, str)
        or re.fullmatch(r"[0-9a-f]{64}", parent_restart_sha256) is None
    ):
        raise CompositeReportError("R14 selected terminal lacks a parent restart")
    parent_restart = require_beneath(
        parent_restart_value,
        parent_segment / "output/rst",
        "R14 selected parent restart",
    )
    if sha256(parent_restart) != parent_restart_sha256:
        raise CompositeReportError("R14 selected parent restart binding is stale")
    restart_time = finite_float(
        report.fast_restart_time(parent_restart), "R14 selected parent restart time"
    )
    start_time = finite_float(
        selected_manifest.get("start_time"), "R14 selected terminal start time"
    )
    parent_final = finite_float(
        parent.get("observed_final_time"), "R14 selected parent final time"
    )
    if not (
        math.isclose(
            restart_time,
            start_time,
            rel_tol=0.0,
            abs_tol=report.RESTART_TIME_TOLERANCE,
        )
        and math.isclose(
            restart_time,
            parent_final,
            rel_tol=0.0,
            abs_tol=report.RESTART_TIME_TOLERANCE,
        )
    ):
        raise CompositeReportError(
            "R14 replay attempts do not start from the selected parent checkpoint"
        )

    failure_time = finite_float(record.get("final_time"), "R14 retained final time")
    attempts: list[dict[str, object]] = []
    for run_root in config.corrected.run_roots:
        case_root = run_root / R14_FAILURE_CASE
        if not case_root.is_dir():
            continue
        for manifest_path in sorted(case_root.glob("fast_s*/manifest/fast_run.json")):
            attempt = authenticate_r14_failed_attempt(
                config,
                case,
                manifest_path,
                selected_manifest.get("variant"),
                selected_overrides,
                parent_restart,
                parent_restart_sha256,
            )
            if attempt is not None:
                attempts.append(attempt)
    attempts.sort(key=lambda value: int(value["sequence"]))
    if (
        len(attempts) < 2
        or len({str(value["segment"]) for value in attempts}) != len(attempts)
        or len({int(value["sequence"]) for value in attempts}) != len(attempts)
        or len({str(value["job_id"]) for value in attempts}) != len(attempts)
    ):
        raise CompositeReportError(
            "R14 requires at least two distinct authenticated failed replay attempts"
        )
    selected_segment = str(Path(str(terminal["segment_dir"])).resolve(strict=True))
    if selected_segment not in {str(value["segment"]) for value in attempts}:
        raise CompositeReportError(
            "R14 selected terminal is not one of the authenticated failed replays"
        )
    if any(
        not math.isclose(
            finite_float(value["failure_time"], "R14 replay failure time"),
            failure_time,
            rel_tol=0.0,
            abs_tol=R14_FAILURE_TIME_TOLERANCE,
        )
        for value in attempts
    ):
        raise CompositeReportError(
            "R14 failed replay attempts do not share one physical failure time"
        )
    mhd_hashes = {str(value["mhd_history_sha256"]) for value in attempts}
    user_hashes = {str(value["user_history_sha256"]) for value in attempts}
    if len(mhd_hashes) != 1 or len(user_hashes) != 1:
        raise CompositeReportError(
            "R14 failed replay histories are not byte-identical"
        )

    retained_history = retained_r14_history(record, failure_time)
    return {
        "schema_version": 1,
        "record_type": "cgl_lf_stage_i_terminal_disposition",
        "case_id": R14_FAILURE_CASE,
        "status": "failed_partial",
        "disposition": "reproducible_finite_time_model_runtime_failure",
        "claim_scope": "authenticated_history_through_failure_time",
        "failure_time": failure_time,
        "failure_time_absolute_tolerance": R14_FAILURE_TIME_TOLERANCE,
        "fatal_signature": R14_FATAL_SIGNATURE,
        "physics": {
            "variant": R14_APPROVED_VARIANT,
            "command_line_overrides": R14_APPROVED_OVERRIDES,
            "cgl_lf_strict_admissibility": False,
            "backup_limiters": False,
            "mirror_limiter": True,
            "firehose_limiter": True,
            "limiter_nu_coll": 20.0,
        },
        "replay_history_consensus": {
            "mhd_sha256": next(iter(mhd_hashes)),
            "user_sha256": next(iter(user_hashes)),
        },
        "selected_parent": {
            "segment": str(parent_segment),
            "sequence": manifests[-2].get("sequence"),
            "restart": {
                **report.artifact_binding(parent_restart),
                "physical_time": restart_time,
            },
        },
        "selected_terminal_segment": selected_segment,
        "attempt_count": len(attempts),
        "attempts": attempts,
        "retained_history": retained_history,
    }


def validate_and_classify_case(
    config: CompositeConfig,
    case_id: str,
    case: dict[str, object],
    record: dict[str, object],
    compatibility: dict[str, object],
    *,
    require_complete: bool,
) -> dict[str, object]:
    """Validate one assembled case and stamp its real execution evidence class."""

    authority = authority_for_case(config, case_id)
    if record.get("case_id") != case_id or record.get("case_name") != case.get("name"):
        raise CompositeReportError(f"{case_id} assembled case identity differs")
    errors = record.get("errors")
    if not isinstance(errors, list) or errors:
        raise CompositeReportError(f"{case_id} assembly retained errors: {errors}")
    if record.get("status") == "assembly_error":
        raise CompositeReportError(f"{case_id} assembly failed")
    status = record.get("status")
    if require_complete and status not in ("complete", "failed_partial"):
        raise CompositeReportError(
            f"{case_id} is not complete: {status}"
        )
    if status == "failed_partial" and case_id != R14_FAILURE_CASE:
        raise CompositeReportError(f"{case_id} is not complete: {status}")
    lineage = record.get("lineage")
    if not isinstance(lineage, list) or not lineage:
        raise CompositeReportError(f"{case_id} selected lineage is missing")
    manifests = [
        validate_segment_manifest(config, authority, case_id, case, item)
        for item in lineage
        if isinstance(item, dict)
    ]
    if len(manifests) != len(lineage):
        raise CompositeReportError(f"{case_id} selected lineage is malformed")
    validate_selected_attempt_sequences(manifests, case_id)
    first = manifests[0]
    if (
        first.get("start_time") != 0.0
        or first.get("restart") is not None
        or first.get("restart_sha256") is not None
    ):
        raise CompositeReportError(f"{case_id} selected lineage does not start fresh at t=0")

    identities = record.get("lineage_identities")
    expected_identities = {
        "input_sha256": [sha256(expected_input(authority, case, case_id))],
        "matrix_sha256": [authority.matrix_sha256],
        "executable_sha256": [authority.executable_sha256],
    }
    if identities != expected_identities:
        raise CompositeReportError(
            f"{case_id} assembled lineage identities differ: {identities}"
        )

    if status == "failed_partial":
        expected_disposition = derive_r14_terminal_disposition(
            config, case, record, manifests
        )
        retained_disposition = record.get("terminal_disposition")
        if retained_disposition is None:
            record["terminal_disposition"] = expected_disposition
        elif retained_disposition != expected_disposition:
            raise CompositeReportError(
                "R14 retained terminal disposition differs from source attempts"
            )
    elif "terminal_disposition" in record:
        raise CompositeReportError(
            f"{case_id} complete lineage retains a terminal disposition"
        )

    execution = case_execution_record(authority, case_id, case, compatibility)
    record["evidence_class"] = authority.evidence_class
    record["execution_source"] = execution["source"]
    record["execution_executable"] = execution["executable"]
    record["execution_authority"] = execution
    return record


def reject_output_symlinks(output: Path) -> None:
    """Reject any symlink in the composite report tree."""

    if output.is_symlink():
        raise CompositeReportError(f"composite output may not be a symlink: {output}")
    if not output.exists():
        return
    for parent, directories, files in os.walk(output, followlinks=False):
        for name in [*directories, *files]:
            path = Path(parent) / name
            if path.is_symlink():
                raise CompositeReportError(
                    f"composite output contains a forbidden symlink: {path}"
                )


def output_is_separate(config: CompositeConfig, output: Path) -> None:
    """Require a fresh report output outside both source run trees."""

    resolved = output.expanduser().absolute().resolve(strict=False)
    for source in (
        *config.corrected.run_roots,
        *config.legacy.run_roots,
        config.corrected.source,
        config.legacy.source,
        REPO_ROOT,
    ):
        try:
            resolved.relative_to(source.resolve(strict=True))
        except ValueError:
            continue
        raise CompositeReportError(f"composite output is inside source tree: {source}")
    if resolved.exists():
        reject_output_symlinks(resolved)
        raise CompositeReportError(f"refusing to replace existing composite output: {resolved}")


def replace_output_root(value: object, old: Path, new: Path) -> object:
    """Rewrite generated staging paths to their immutable final output root."""

    if isinstance(value, dict):
        return {key: replace_output_root(item, old, new) for key, item in value.items()}
    if isinstance(value, list):
        return [replace_output_root(item, old, new) for item in value]
    if isinstance(value, str):
        old_text = str(old)
        if value == old_text or value.startswith(old_text + os.sep):
            return str(new) + value[len(old_text):]
    return value


def build_inventory(
    config: CompositeConfig,
    output: Path,
    validation: dict[str, object],
    records: dict[str, dict[str, object]],
) -> dict[str, object]:
    """Build the composite inventory while retaining base reporter compatibility."""

    return {
        "schema_version": 1,
        "record_type": "cgl_lf_stage_i_corrected_composite_report",
        "assembled_utc": report.utc_now(),
        "root": str(config.corrected.root.resolve(strict=True)),
        "legacy_control_root": str(config.legacy.root.resolve(strict=True)),
        "output": str(output),
        "matrix": report.artifact_binding(config.corrected.matrix),
        # Existing downstream report/acceptance helpers authenticate this engine.
        "adapter": report.artifact_binding(REPORTER_PATH),
        "composite_adapter": report.artifact_binding(SCRIPT_PATH),
        "corrected_identity": validation["corrected_identity"],
        "corrected_identity_validation": validation["corrected_identity_validation"],
        "authorities": validation["authorities"],
        "passive_compatibility": validation["passive_compatibility"],
        "evidence_classes": {
            config.corrected.evidence_class: list(ACTIVE_CASES),
            config.legacy.evidence_class: list(PASSIVE_CONTROLS),
        },
        "searched_fast_run_roots": [
            *(
                str(run_root.resolve(strict=True))
                for run_root in config.corrected.run_roots
            ),
            *(str(run_root.resolve(strict=True)) for run_root in config.legacy.run_roots),
        ],
        "cases": records,
    }


def assemble_composite(
    config: CompositeConfig,
    output: Path,
) -> Path:
    """Assemble all corrected active cases and authenticated legacy passive controls."""

    validation = validate_authorities(config)
    output = output.expanduser().absolute().resolve(strict=False)
    output_is_separate(config, output)
    cases = validation["cases"]
    compatibility = validation["passive_compatibility"]
    if not isinstance(cases, dict) or not isinstance(compatibility, dict):
        raise CompositeReportError("authority validation result is malformed")

    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=f".{output.name}.assemble-", dir=output.parent)
    ).resolve(strict=True)
    try:
        records: dict[str, dict[str, object]] = {}
        for case_id in ALL_CASES:
            case = cases[case_id]
            if not isinstance(case, dict):
                raise CompositeReportError(f"{case_id} matrix case is malformed")
            authority = authority_for_case(config, case_id)
            record = report.assemble_fast_case(
                authority.root, authority.source, staging, case_id, case
            )
            records[case_id] = validate_and_classify_case(
                config,
                case_id,
                case,
                record,
                compatibility,
                require_complete=True,
            )
            write_json(staging / "cases" / case_id / "lineage.json", records[case_id])

        rewritten = replace_output_root(records, staging, output)
        if not isinstance(rewritten, dict):
            raise AssertionError("rewritten records must remain a dictionary")
        records = rewritten
        for case_id, record in records.items():
            write_json(staging / "cases" / case_id / "lineage.json", record)
        inventory = build_inventory(config, output, validation, records)
        write_json(staging / "inventory.json", inventory)
        write_json(
            staging / "manifest.json",
            {
                "schema_version": 1,
                "record_type": "cgl_lf_stage_i_corrected_composite_report",
                "created_utc": report.utc_now(),
                "output": str(output),
                "matrix": inventory["matrix"],
                "adapter": inventory["adapter"],
                "composite_adapter": inventory["composite_adapter"],
                "corrected_identity": inventory["corrected_identity"],
                "commands": [
                    "preflight",
                    "assemble",
                    "validate",
                    "analyze-case",
                    "aggregate",
                    "render",
                    "manuscript",
                    "verify",
                ],
            },
        )
        reject_output_symlinks(staging)
        os.replace(staging, output)
        staging = Path()
    finally:
        if staging != Path() and staging.exists():
            shutil.rmtree(staging)
    validate_composite_inventory(config, output)
    return output / "inventory.json"


def validate_case_inventory_record(
    config: CompositeConfig,
    output: Path,
    cases: dict[str, dict[str, object]],
    compatibility: dict[str, object],
    case_id: str,
    inline: object,
) -> None:
    """Authenticate one retained composite inventory case record."""

    if not isinstance(inline, dict):
        raise CompositeReportError(f"{case_id} composite inventory record is malformed")
    path = output / "cases" / case_id / "lineage.json"
    retained = load_json(path, f"{case_id} retained lineage")
    if retained != inline:
        raise CompositeReportError(f"{case_id} retained lineage differs from inventory")
    expected_class = authority_for_case(config, case_id).evidence_class
    if inline.get("evidence_class") != expected_class:
        raise CompositeReportError(f"{case_id} evidence class differs")
    expected_execution = case_execution_record(
        authority_for_case(config, case_id), case_id, cases[case_id], compatibility
    )
    if inline.get("execution_authority") != expected_execution:
        raise CompositeReportError(f"{case_id} execution authority differs")
    if inline.get("execution_source") != expected_execution["source"]:
        raise CompositeReportError(f"{case_id} execution source differs")
    if inline.get("execution_executable") != expected_execution["executable"]:
        raise CompositeReportError(f"{case_id} execution executable differs")
    validate_and_classify_case(
        config,
        case_id,
        cases[case_id],
        inline,
        compatibility,
        require_complete=True,
    )


def validate_composite_inventory(
    config: CompositeConfig, output: Path
) -> dict[str, object]:
    """Revalidate a retained composite inventory and all case execution authorities."""

    output = canonical_existing(output, "composite output", directory=True)
    reject_output_symlinks(output)
    validation = validate_authorities(config)
    inventory_path = output / "inventory.json"
    inventory = load_json(inventory_path, "composite inventory")
    if inventory.get("record_type") != "cgl_lf_stage_i_corrected_composite_report":
        raise CompositeReportError("composite inventory record type differs")
    if inventory.get("output") != str(output):
        raise CompositeReportError("composite inventory output root differs")
    if inventory.get("root") != str(config.corrected.root.resolve(strict=True)):
        raise CompositeReportError("composite inventory corrected root differs")
    if inventory.get("legacy_control_root") != str(config.legacy.root.resolve(strict=True)):
        raise CompositeReportError("composite inventory legacy root differs")
    if inventory.get("adapter") != report.artifact_binding(REPORTER_PATH):
        raise CompositeReportError("composite inventory reporter binding differs")
    if inventory.get("composite_adapter") != report.artifact_binding(SCRIPT_PATH):
        raise CompositeReportError("composite inventory adapter binding differs")
    if inventory.get("matrix") != report.artifact_binding(config.corrected.matrix):
        raise CompositeReportError("composite inventory matrix binding differs")
    expected_classes = {
        config.corrected.evidence_class: list(ACTIVE_CASES),
        config.legacy.evidence_class: list(PASSIVE_CONTROLS),
    }
    if inventory.get("evidence_classes") != expected_classes:
        raise CompositeReportError("composite inventory evidence classes differ")
    expected_search_roots = [
        *(str(run_root.resolve(strict=True)) for run_root in config.corrected.run_roots),
        *(str(run_root.resolve(strict=True)) for run_root in config.legacy.run_roots),
    ]
    if inventory.get("searched_fast_run_roots") != expected_search_roots:
        raise CompositeReportError("composite inventory searched run roots differ")
    for key in (
        "corrected_identity",
        "corrected_identity_validation",
        "authorities",
        "passive_compatibility",
    ):
        if inventory.get(key) != validation[key]:
            raise CompositeReportError(f"composite inventory {key} differs")
    cases = validation["cases"]
    compatibility = validation["passive_compatibility"]
    inline_cases = inventory.get("cases")
    if not isinstance(cases, dict) or not isinstance(compatibility, dict):
        raise CompositeReportError("authority validation result is malformed")
    if not isinstance(inline_cases, dict) or set(inline_cases) != set(ALL_CASES):
        raise CompositeReportError("composite inventory case coverage differs")
    for case_id in ALL_CASES:
        validate_case_inventory_record(
            config, output, cases, compatibility, case_id, inline_cases[case_id]
        )
    return {
        "result": "pass",
        "inventory": report.artifact_binding(inventory_path),
        "output": str(output),
        "cases": list(ALL_CASES),
        "evidence_classes": inventory["evidence_classes"],
        "terminal_dispositions": {
            case_id: inline_cases[case_id]["terminal_disposition"]
            for case_id in ALL_CASES
            if isinstance(inline_cases[case_id], dict)
            and inline_cases[case_id].get("terminal_disposition") is not None
        },
    }


def corrected_verification(
    output: Path,
    validation: dict[str, object],
) -> dict[str, object]:
    """Remove only the base-verifier errors authorized by the R14 disposition."""

    base_path = output / "verify.base.json"
    base = load_json(base_path, "base verification")
    inventory = load_json(output / "inventory.json", "composite inventory")
    cases = inventory.get("cases")
    dispositions = validation.get("terminal_dispositions")
    if not isinstance(cases, dict) or not isinstance(dispositions, dict):
        raise CompositeReportError(
            "corrected verification evidence is malformed"
        )
    if set(dispositions) not in (set(), {R14_FAILURE_CASE}):
        raise CompositeReportError(
            "corrected verification has unsupported terminal dispositions"
        )

    case_results = base.get("cases")
    base_errors = base.get("errors")
    base_warnings = base.get("warnings")
    if (
        not isinstance(case_results, dict)
        or not isinstance(base_errors, list)
        or not isinstance(base_warnings, list)
    ):
        raise CompositeReportError("base verification record is malformed")
    flattened = [
        f"{case_id}: {error}"
        for case_id, value in case_results.items()
        if isinstance(value, dict)
        for error in value.get("errors", [])
    ]
    if base_errors != flattened:
        raise CompositeReportError(
            "base verification top-level errors differ from case errors"
        )

    adjusted = json.loads(json.dumps(base))
    accepted_top: list[str] = []
    if dispositions:
        r14 = cases.get(R14_FAILURE_CASE)
        disposition = dispositions[R14_FAILURE_CASE]
        if not isinstance(r14, dict) or not isinstance(disposition, dict):
            raise CompositeReportError("R14 corrected verification evidence is malformed")
        lineage = r14.get("lineage")
        attempts = disposition.get("attempts")
        if (
            not isinstance(lineage, list)
            or not lineage
            or not isinstance(attempts, list)
        ):
            raise CompositeReportError("R14 corrected verification lineage is malformed")
        terminal = lineage[-1]
        if not isinstance(terminal, dict):
            raise CompositeReportError("R14 corrected verification terminal is malformed")
        terminal_segment = str(Path(str(terminal.get("segment_dir"))).resolve())
        exit_code = terminal.get("run_exit_code")
        matching = [
            value
            for value in attempts
            if isinstance(value, dict)
            and str(Path(str(value.get("segment"))).resolve()) == terminal_segment
        ]
        if (
            terminal.get("state") != "failed"
            or terminal_segment != disposition.get("selected_terminal_segment")
            or not isinstance(exit_code, int)
            or isinstance(exit_code, bool)
            or exit_code == 0
            or len(matching) != 1
            or matching[0].get("run_exit_code") != exit_code
        ):
            raise CompositeReportError(
                "R14 corrected verification terminal differs from authenticated attempts"
            )
        accepted = [
            (
                f"selected lineage segment {len(lineage) - 1} has nonzero "
                f"run exit code: {exit_code}"
            ),
        ]
        if base.get("require_complete") is True:
            accepted.append("case is not complete: failed_partial")
        r14_result = case_results.get(R14_FAILURE_CASE)
        if not isinstance(r14_result, dict) or not isinstance(
            r14_result.get("errors"), list
        ):
            raise CompositeReportError("base R14 verification result is malformed")
        if any(r14_result["errors"].count(error) != 1 for error in accepted):
            raise CompositeReportError(
                "base verification lacks the exact authenticated R14 errors"
            )
        adjusted_r14 = adjusted["cases"][R14_FAILURE_CASE]
        adjusted_r14["errors"] = [
            error for error in adjusted_r14["errors"] if error not in accepted
        ]
        accepted_top = [f"R14: {error}" for error in accepted]
        adjusted["errors"] = [
            error for error in adjusted["errors"] if error not in accepted_top
        ]
    adjusted["result"] = (
        "fail"
        if adjusted["errors"]
        else "warnings"
        if adjusted["warnings"]
        else "pass"
    )
    adjusted.update(
        {
            "record_type": "cgl_lf_stage_i_corrected_composite_verification",
            "adapter": report.artifact_binding(SCRIPT_PATH),
            "base_adapter": base.get("adapter"),
            "base_verification": report.artifact_binding(base_path),
            "terminal_dispositions": dispositions,
            "accepted_base_errors": accepted_top,
        }
    )
    return adjusted


def command_verify(output: Path, values: list[str]) -> int:
    """Run the strict base verifier and adapt only authenticated R14 errors."""

    validation = validate_composite_inventory(DEFAULT_CONFIG, output)
    report.main(reporter_arguments(output, "verify", values))
    verify_path = output / "verify.json"
    if not verify_path.is_file():
        raise CompositeReportError("base verifier did not write verify.json")
    base_path = output / "verify.base.json"
    os.replace(verify_path, base_path)
    adjusted = corrected_verification(output, validation)
    write_json(verify_path, adjusted)
    return 1 if adjusted["errors"] else 0


def reporter_arguments(output: Path, command: str, values: Iterable[str]) -> list[str]:
    """Build a delegated base-reporter command over the composite output."""

    return [
        "--root",
        str(DEFAULT_CONFIG.corrected.root),
        "--output",
        str(output),
        "--frozen-source",
        str(DEFAULT_CONFIG.corrected.source),
        "--matrix",
        str(DEFAULT_CONFIG.corrected.matrix),
        command,
        *values,
    ]


def build_parser() -> argparse.ArgumentParser:
    """Build the composite report command-line interface."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    commands = parser.add_subparsers(dest="command", required=True)

    commands.add_parser("preflight", help="authenticate both execution authorities")
    assemble = commands.add_parser("assemble", help="assemble the composite inventory")
    commands.add_parser("validate", help="revalidate the composite inventory")

    analyze = commands.add_parser("analyze-case", help="delegate one case analysis")
    analyze.add_argument("arguments", nargs=argparse.REMAINDER)
    aggregate = commands.add_parser("aggregate", help="delegate campaign aggregation")
    aggregate.add_argument("arguments", nargs=argparse.REMAINDER)
    commands.add_parser("render", help="delegate figure rendering")
    commands.add_parser("manuscript", help="delegate manuscript-ready product rendering")
    verify = commands.add_parser("verify", help="delegate base report verification")
    verify.add_argument("arguments", nargs=argparse.REMAINDER)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the corrected composite report adapter without submitting jobs."""

    args = build_parser().parse_args(argv)
    try:
        if args.command == "preflight":
            validation = validate_authorities(DEFAULT_CONFIG)
            print(
                json.dumps(
                    {
                        "result": "pass",
                        "case_count": len(validation["cases"]),
                        "authorities": validation["authorities"],
                        "corrected_identity": validation["corrected_identity"],
                        "passive_compatibility": validation["passive_compatibility"],
                    },
                    indent=2,
                    sort_keys=True,
                )
            )
            return 0
        if args.command == "assemble":
            inventory = assemble_composite(DEFAULT_CONFIG, args.output)
            print(inventory)
            return 0
        if args.command == "validate":
            print(
                json.dumps(
                    validate_composite_inventory(DEFAULT_CONFIG, args.output),
                    indent=2,
                    sort_keys=True,
                )
            )
            return 0
        validate_composite_inventory(DEFAULT_CONFIG, args.output)
        values = getattr(args, "arguments", [])
        result = (
            command_verify(args.output, values)
            if args.command == "verify"
            else report.main(reporter_arguments(args.output, args.command, values))
        )
        validate_composite_inventory(DEFAULT_CONFIG, args.output)
        return result
    except (
        CompositeReportError,
        report.ReportError,
        OSError,
        ValueError,
        KeyError,
        json.JSONDecodeError,
        subprocess.SubprocessError,
    ) as error:
        print(f"corrected composite report error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
