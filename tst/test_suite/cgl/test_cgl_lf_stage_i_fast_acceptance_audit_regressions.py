"""Regression tests for audited direct-fast Stage I acceptance blockers."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import re
import sys

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
ADAPTER = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_fast_acceptance.py"
PUBLICATION = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_fast_publication.py"
WORKFLOW = REPOSITORY / "scripts/cgl_lf_workflow.py"
MATRIX = REPOSITORY / "inputs/cgl_lf_paper/mks24_stage_i_manifest.json"
HISTORY_LABEL = re.compile(r"\[(\d+)\]=(\S+)")


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def fast_acceptance():
    return load_module("cgl_lf_stage_i_fast_acceptance_audit_regressions", ADAPTER)


@pytest.fixture(scope="module")
def publication():
    return load_module("cgl_lf_stage_i_fast_publication_audit_regressions", PUBLICATION)


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def plain_binding(path: Path) -> dict[str, object]:
    resolved = path.resolve()
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256(resolved),
    }


def matrix_cases() -> dict[str, dict[str, object]]:
    manifest = json.loads(MATRIX.read_text(encoding="utf-8"))
    return {str(case["id"]): case for case in manifest["cases"]}


def expected_model_choices(input_path: Path, overrides: list[str]) -> dict[str, str]:
    workflow = load_module("cgl_lf_stage_i_fast_acceptance_test_workflow", WORKFLOW)
    source = input_path.read_text(encoding="utf-8")
    choices = {
        str(key): str(value)
        for key, value in workflow.model_choices(source, overrides).items()
    }
    strict = workflow.input_block_value(
        source, "mhd", "cgl_lf_strict_admissibility"
    ) or "false"
    for override in overrides:
        prefix = "mhd/cgl_lf_strict_admissibility="
        if override.startswith(prefix):
            strict = override[len(prefix):]
    choices["cgl_lf_strict_admissibility"] = strict
    return choices


def write_history(path: Path, columns: dict[str, list[float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    labels = list(columns)
    rows = zip(*(columns[label] for label in labels))
    text = "# " + " ".join(
        f"[{index}]={label}" for index, label in enumerate(labels, start=1)
    ) + "\n"
    text += "\n".join(
        " ".join(format(value, ".17g") for value in row) for row in rows
    )
    path.write_text(text + "\n", encoding="utf-8")


def history_fixture(
    root: Path,
    *,
    counters: bool = True,
    final_time: float = 10.0,
) -> tuple[Path, Path]:
    times = [0.0, 4.0, 6.0, 8.0, final_time]
    mhd: dict[str, list[float]] = {"time": times}
    if counters:
        for name in ("lf_dfloor", "lf_pfloor", "lf_nonfin", "lf_nonpos", "lf_hardbd"):
            mhd[name] = [0.0] * len(times)
    user = {
        "time": times,
        "kinetic": [1.0, 1.1, 1.2, 1.1, 1.0],
    }
    mhd_path = root / "fixture.mhd.hst"
    user_path = root / "fixture.user.hst"
    write_history(mhd_path, mhd)
    write_history(user_path, user)
    return mhd_path, user_path


def history_record(path: Path, *, declared_sha256: str | None = None) -> dict[str, object]:
    return {
        "available": True,
        "path": str(path.absolute()),
        "binding": {
            "path": str(path.absolute()),
            "size_bytes": path.stat().st_size,
            "sha256": declared_sha256 or sha256(path),
        },
    }


def lineage_fixture(
    root: Path,
    *,
    case_id: str = "R02",
    status: str = "complete",
    counters: bool = True,
    final_time: float = 10.0,
    strict: bool = True,
    variants: list[str] | None = None,
    overrides: list[str] | None = None,
) -> dict[str, object]:
    mhd, user = history_fixture(root, counters=counters, final_time=final_time)
    executable = root / "qualified-athena"
    executable.write_bytes(b"qualified executable fixture\n")
    selected_overrides = overrides or []
    matrix_case = matrix_cases()[case_id]
    input_path = REPOSITORY / str(matrix_case["input"])
    model_choices = expected_model_choices(input_path, selected_overrides)
    model_choices["cgl_lf_strict_admissibility"] = str(strict).lower()
    return {
        "case_id": case_id,
        "case_name": matrix_case["name"],
        "matrix_case": matrix_case,
        "status": status,
        "target_time": 10.0,
        "errors": [],
        "lineage_variants": variants or [],
        "lineage_command_line_overrides": selected_overrides,
        "lineage_claim_scopes": ["standard"],
        "model_choices": model_choices,
        "input": {
            "path": str(input_path.resolve()),
            "size_bytes": input_path.stat().st_size,
            "sha256": sha256(input_path),
        },
        "lineage_identities": {
            "input_sha256": [sha256(input_path)],
            "matrix_sha256": [sha256(MATRIX)],
            "executable_sha256": [sha256(executable)],
        },
        "test_qualified_executable": plain_binding(executable),
        "histories": {
            "mhd": history_record(mhd),
            "user": history_record(user),
        },
    }


def policy_fixture(
    required: list[str] | None = None,
    lineage: dict[str, object] | None = None,
) -> dict[str, object]:
    policy = {
        "criteria_binding": {"path": "criteria.json", "sha256": "c" * 64},
        "review_binding": {"path": "review.json", "sha256": "d" * 64},
        "manifest": json.loads(MATRIX.read_text(encoding="utf-8")),
        "verified_sources": {
            "stage_i_manifest": {
                "path": str(MATRIX.resolve()),
                "size_bytes": MATRIX.stat().st_size,
                "sha256": sha256(MATRIX),
            },
        },
        "criteria": {
            "required_cases": required or ["R02"],
            "analysis_windows": {
                "full": [4.0, 10.0],
                "early": [4.0, 8.0],
                "late": [6.0, 10.0],
            },
            "statistics_policy": {
                "bootstrap_replicates": 8,
                "gap_policy": {
                    "expected_history_cadence": 0.25,
                    "maximum_gap_expected_cadence_multiplier": 8.0,
                    "maximum_gap_forcing_tcorr_fraction": 0.25,
                },
            },
            "case_metrics": {
                "kinetic": {
                    "history": "user",
                    "column": "kinetic",
                    "stationarity_kind": "energy",
                },
            },
            "family_gates": {
                "active_passive": {"pairs": []},
            },
        },
    }
    if lineage is not None:
        executable = lineage["test_qualified_executable"]
        approval_path = Path(str(executable["path"])).parent / "qualification.json"
        write_json(
            approval_path,
            {
                "schema_version": 1,
                "execution_epoch": "E03-forcing-policy",
                "approved_executable": executable["path"],
                "approved_executable_sha256": executable["sha256"],
                "approved_executable_revision": "fixture-revision",
            },
        )
        policy["verified_sources"]["qualification_approval"] = plain_binding(
            approval_path
        )
    return policy


class FakeAcceptanceError(RuntimeError):
    """Fixture equivalent of the reviewed utility's input error."""


class FakeAcceptance:
    """Small reviewed-acceptance API stub for adapter contract tests."""

    AcceptanceError = FakeAcceptanceError

    def __init__(self, policy: dict[str, object]):
        self.policy = policy

    def load_validated_policy(self, _criteria: Path, _review: Path) -> dict[str, object]:
        return self.policy

    @staticmethod
    def case_name(policy: dict[str, object], case_id: str) -> str:
        for case in policy["manifest"]["cases"]:
            if case["id"] == case_id:
                return str(case["name"])
        raise FakeAcceptanceError(f"missing fixture case: {case_id}")

    @staticmethod
    def canonical_bundle_path(case_id: str) -> Path:
        return Path(f"/nonexistent/canonical/{case_id}/manifest.json")

    @staticmethod
    def authenticate_case_bundle(
        _policy: dict[str, object],
        case_id: str,
        bundle_path: Path,
        _mhd_binding: dict[str, object],
        _user_binding: dict[str, object],
    ) -> tuple[dict[str, object], list[dict[str, object]], dict[str, object]]:
        bundle = {
            "bundle_manifest": plain_binding(bundle_path),
            "case_id": case_id,
            "canonical_campaign_authority_eligible": False,
        }
        return bundle, [bundle["bundle_manifest"]], {
            "name": "accepted_case_bundle_lineage",
            "result": "pass",
            "reason": "fixture accepted bundle authenticated",
            "observations": bundle,
        }

    def load_history(
        self, path: Path, _label: str
    ) -> tuple[dict[str, list[float]], dict[str, object]]:
        payload = path.read_text(encoding="utf-8")
        labels: list[str] | None = None
        rows: list[list[float]] = []
        for line in payload.splitlines():
            if line.startswith("#"):
                found = HISTORY_LABEL.findall(line)
                if found:
                    labels = [
                        name for _, name in sorted(found, key=lambda item: int(item[0]))
                    ]
            elif line.strip():
                rows.append([float(value) for value in line.split()])
        assert labels is not None
        return (
            {
                label: [row[index] for row in rows]
                for index, label in enumerate(labels)
            },
            {
                "path": str(path.resolve()),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            },
        )

    @staticmethod
    def _forcing_tcorr_from_bundle(value: object) -> float:
        if isinstance(value, dict):
            candidate = value.get("forcing_tcorr")
            return float(candidate) if candidate is not None else 0.0
        if isinstance(value, (str, Path)) and Path(value).is_file():
            record = json.loads(Path(value).read_text(encoding="utf-8"))
            for key in ("forcing_tcorr", "model_choices"):
                candidate = record.get(key)
                if isinstance(candidate, dict):
                    candidate = candidate.get("forcing_tcorr")
                if candidate is not None:
                    return float(candidate)
        return 0.0

    def evaluate_case(self, *args, **kwargs) -> dict[str, object]:
        case_id = str(args[1])
        bundle = kwargs.get("bundle_path", args[6] if len(args) > 6 else None)
        tcorr = self._forcing_tcorr_from_bundle(bundle)
        statistics = self._statistics(tcorr)
        return {
            "record_type": "stage-i-scientific-case-evidence",
            "case_id": case_id,
            "result": "pass",
            "metrics": {
                "kinetic": {
                    "full": statistics,
                    "early": statistics,
                    "late": statistics,
                    "sampling_adequacy": "pass",
                    "stationarity": {"result": "pass"},
                },
            },
        }

    @staticmethod
    def _statistics(minimum_block_duration: float) -> dict[str, object]:
        return {
            "mean": 1.0,
            "standard_deviation": 0.1,
            "standard_error": 0.05,
            "confidence_interval_95": [0.9, 1.1],
            "effective_sample_count": 4.0,
            "independent_time_block_count": 2,
            "gap_adequacy": "pass",
            "sample_count": 5,
            "method": {
                "required_block_duration": minimum_block_duration,
                "minimum_block_duration": minimum_block_duration,
            },
        }

    def window_statistics(
        self,
        _times: list[float],
        _values: list[float],
        _start: float,
        _end: float,
        *,
        minimum_block_duration: float,
        **_kwargs,
    ) -> dict[str, object]:
        return self._statistics(minimum_block_duration)

    def metric_statistics(
        self,
        _history: dict[str, list[float]],
        _values: list[float],
        _metric: str,
        _policy: dict[str, object],
        *,
        kind: str,
        minimum_block_duration: float,
    ) -> dict[str, object]:
        statistics = self._statistics(minimum_block_duration)
        return {
            "full": statistics,
            "early": statistics,
            "late": statistics,
            "sampling_adequacy": "pass",
            "stationarity": {"result": "pass", "kind": kind},
        }

    @staticmethod
    def analyzer_metrics(
        _diagnostics: dict[str, object],
        _name: str,
        _policy: dict[str, object],
        _minimum_block_duration: float,
    ) -> dict[str, object]:
        return {}

    @staticmethod
    def convergence_products(
        _diagnostics: dict[str, object],
        _name: str,
    ) -> dict[str, object]:
        return {
            "velocity_spectrum_shape": {
                "x": [1.0, 2.0, 3.0],
                "y": [3.0, 2.0, 1.0],
            },
        }

    @staticmethod
    def active_energy_closure_gate(
        _policy: dict[str, object],
        case_id: str,
        _mhd: dict[str, list[float]],
        _user: dict[str, list[float]],
    ) -> dict[str, object]:
        return {
            "name": "active_energy_closure",
            "result": "pass",
            "reason": f"fixture active-energy closure passed for {case_id}",
            "observations": {},
        }

    @staticmethod
    def history_delta(
        _history: dict[str, list[float]],
        _column: str,
        _start: float,
        _end: float,
    ) -> float:
        return 1.0

    @staticmethod
    def gate(
        name: str,
        result: str,
        *,
        reason: str,
        observations: object = None,
        limits: object = None,
    ) -> dict[str, object]:
        return {
            "name": name,
            "result": result,
            "reason": reason,
            "observations": observations,
            "limits": limits,
        }

    @staticmethod
    def aggregate_gate_result(gates: list[dict[str, object]]) -> str:
        results = [str(gate["result"]) for gate in gates]
        if "fail" in results:
            return "fail"
        if "inconclusive" in results:
            return "inconclusive"
        return "pass"

    @staticmethod
    def seal_evidence(value: dict[str, object]) -> dict[str, object]:
        return value

    @staticmethod
    def pair_contrast(
        _active: dict[str, object],
        _passive: dict[str, object],
        _policy: dict[str, object],
    ) -> dict[str, object]:
        return {"result": "pass", "reason": "fixture pair passed", "metrics": {}}

    @staticmethod
    def scalar_from_case(
        _case: dict[str, object], _metric: str, _window: str
    ) -> dict[str, float]:
        return {"mean": 2.0, "standard_error": 0.01, "standard_deviation": 0.1}


class FakeSnapshotAnalyzer:
    """Replay fixture for one retained snapshot ensemble."""

    def __init__(self, ensemble: dict[str, object]):
        self.ensemble = ensemble

    def average_snapshot_records(
        self, _records: dict[str, object]
    ) -> dict[str, object]:
        return dict(self.ensemble)

    @staticmethod
    def validate_snapshot_provenance_record(
        provenance: object, _context: str
    ) -> dict[str, object]:
        assert isinstance(provenance, dict)
        assert provenance["layout"] == "single_file"
        assert provenance["expected_rank_count"] == 1
        assert provenance["rank_directory_names"] == []
        assert provenance["representative_path"] == provenance["files"][0]["path"]
        return provenance

    @staticmethod
    def revalidate_snapshot_provenance_record(
        provenance: dict[str, object], expected_ranks: int, _context: str
    ) -> None:
        assert expected_ranks == 1
        assert provenance["files"] == [
            {
                **plain_binding(Path(str(provenance["representative_path"]))),
                "symlink_target": None,
            }
        ]


class FakeReport:
    """Replay fixture for direct-fast history and snapshot diagnostics."""

    def __init__(
        self,
        windows: dict[str, object],
        ensemble: dict[str, object],
        records: dict[str, object] | None = None,
    ):
        self.windows = windows
        self.ensemble = json.loads(json.dumps(ensemble))
        self.records = json.loads(json.dumps(records or {}))
        self.analyzer = FakeSnapshotAnalyzer(self.ensemble)

    def load_pure_analyzer(self) -> FakeSnapshotAnalyzer:
        return self.analyzer

    def window_summaries(self, *_args) -> dict[str, object]:
        return self.windows

    def analyze_snapshot_paths_bounded(self, *_args) -> tuple[
        dict[str, object], dict[str, object]
    ]:
        return (
            json.loads(json.dumps(self.records)),
            json.loads(json.dumps(self.ensemble)),
        )

    @staticmethod
    def merge_histories(
        sources: list[tuple[str, Path]], destination: Path, _label: str
    ) -> dict[str, object]:
        assert len(sources) == 1
        destination.write_bytes(sources[0][1].read_bytes())
        return {"available": True, "errors": []}


def add_authenticated_execution_fixture(
    fast_acceptance,
    root: Path,
    policy: dict[str, object],
    lineage: dict[str, object],
) -> None:
    """Add one exact accepted-bundle or qualified-fast execution lineage."""

    case_id = str(lineage["case_id"])
    histories = lineage["histories"]
    source_bindings: dict[str, dict[str, object]] = {}
    if case_id == "R02":
        bundle_manifest = root / "accepted-bundle/manifest.json"
        write_json(bundle_manifest, {"fixture": "accepted R02 bundle"})
        for kind in ("mhd", "user"):
            source = bundle_manifest.parent / "history" / f"accepted.{kind}.hst"
            source.parent.mkdir(parents=True, exist_ok=True)
            source.write_bytes(Path(str(histories[kind]["path"])).read_bytes())
            source_bindings[kind] = plain_binding(source)
            histories[kind]["sources"] = [source_bindings[kind]]
        lineage["lineage"] = [{
            "kind": "accepted_r02_bundle",
            "order": 0,
            "segment": "R02 accepted bundle",
            "output": str(bundle_manifest.parent),
            "manifest": plain_binding(bundle_manifest),
            "mhd_history": source_bindings["mhd"]["path"],
            "user_history": source_bindings["user"]["path"],
        }]
        return

    executable = lineage["test_qualified_executable"]
    approval_path = Path(str(executable["path"])).parent / "qualification.json"
    write_json(
        approval_path,
        {
            "schema_version": 1,
            "execution_epoch": "E03-forcing-policy",
            "approved_executable": executable["path"],
            "approved_executable_sha256": executable["sha256"],
            "approved_executable_revision": "fixture-revision",
        },
    )
    policy["verified_sources"]["qualification_approval"] = plain_binding(approval_path)

    segment = root / "selected-fast" / case_id / "fast_s000_t0_to_t10"
    output = segment / "output"
    output.mkdir(parents=True)
    for kind in ("mhd", "user"):
        source = output / f"selected.{kind}.hst"
        source.write_bytes(Path(str(histories[kind]["path"])).read_bytes())
        source_bindings[kind] = plain_binding(source)
        histories[kind]["sources"] = [source_bindings[kind]]
    variants = lineage.get("lineage_variants", [])
    variant = variants[0] if variants else None
    overrides = list(lineage.get("lineage_command_line_overrides", []))
    manifest_path = segment / "manifest/fast_run.json"
    write_json(
        manifest_path,
        {
            "schema_version": 1,
            "case_id": case_id,
            "case_name": lineage["case_name"],
            "input": lineage["input"]["path"],
            "input_sha256": lineage["lineage_identities"]["input_sha256"][0],
            "matrix_sha256": lineage["lineage_identities"]["matrix_sha256"][0],
            "executable": executable["path"],
            "executable_sha256": executable["sha256"],
            "run_dir": str(segment),
            "output_dir": str(output),
            "sequence": 0,
            "start_time": 0.0,
            "target_time": 10.0,
            "restart": None,
            "restart_sha256": None,
            "variant": variant,
            "command_line_overrides": overrides,
        },
    )
    exit_path = segment / "manifest/run_exit_code"
    exit_path.write_text("0\n", encoding="utf-8")
    lineage["lineage"] = [{
        "kind": "fast",
        "order": 0,
        "segment": "fast_s000_t0_to_t10",
        "segment_dir": str(segment),
        "output": str(output),
        "observed_final_time": 10.0,
        "manifest": plain_binding(manifest_path),
        "input_sha256": lineage["lineage_identities"]["input_sha256"][0],
        "matrix_sha256": lineage["lineage_identities"]["matrix_sha256"][0],
        "executable_sha256": executable["sha256"],
        "variant": variant,
        "command_line_overrides": overrides,
        "run_exit_code": 0,
        "run_exit_code_artifact": plain_binding(exit_path),
        "mhd_history": source_bindings["mhd"]["path"],
        "user_history": source_bindings["user"]["path"],
    }]

def write_complete_diagnostics(
    fast_acceptance,
    case_dir: Path,
    lineage: dict[str, object],
) -> FakeReport:
    """Write one fully bound direct-fast diagnostics fixture."""

    lineage_path = case_dir / "lineage.json"
    snapshot_index = case_dir / "snapshots.json"
    snapshot_source = case_dir / "snapshot.athdf"
    snapshot_source.write_bytes(b"snapshot fixture\n")
    write_json(
        snapshot_index,
        {
            "schema_version": 1,
            "snapshots": [{
                "complete": True,
                "time": 9.0,
                "expected_ranks": 1,
                "representative": str(snapshot_source.resolve()),
            }],
        },
    )
    source_binding = {
        **fast_acceptance.binding(snapshot_source),
        "symlink_target": None,
    }
    files = [source_binding]
    provenance = {
        "files": files,
        "layout": "single_file",
        "expected_rank_count": 1,
        "rank_directory_names": [],
        "representative_path": str(snapshot_source.resolve()),
        "snapshot_time": 9.0,
        "aggregate_sha256": hashlib.sha256(
            fast_acceptance.canonical_json_bytes(files)
        ).hexdigest(),
    }
    windows = {"fixture": "replayed history diagnostics"}
    ensemble = {
        "snapshot_count": 1,
        "time_start": 8.0,
        "time_end": 10.0,
        "spectra": {
            "velocity": {
                "k": [1.0, 2.0, 3.0],
                "power_per_dk": [3.0, 2.0, 1.0],
            },
        },
    }
    records = {
        str(snapshot_source.resolve()): {
            "time": 9.0,
            "snapshot_provenance": provenance,
        },
    }
    write_json(
        case_dir / "diagnostics.json",
        {
            "schema_version": 1,
            "case_id": lineage["case_id"],
            "case_name": lineage["case_name"],
            "assembly_status": "complete",
            "analysis_status": "complete",
            "model_choices": lineage["model_choices"],
            "health": {"result": "clean", "structural_errors": []},
            "selected_snapshot_count": 1,
            "snapshot_analysis_status": "complete",
            "analysis_errors": [],
            "windows": windows,
            "snapshots": records,
            "snapshot_ensemble": ensemble,
            "compat": {
                "analysis_window": {"time_start": 8.0, "time_end": 10.0},
                "snapshot_ensemble": ensemble,
            },
            "provenance": {
                "lineage": fast_acceptance.binding(lineage_path),
                "snapshot_index": fast_acceptance.binding(snapshot_index),
                "merged_mhd_history": lineage["histories"]["mhd"]["binding"],
                "merged_user_history": lineage["histories"]["user"]["binding"],
                "analyzer": fast_acceptance.binding(fast_acceptance.PAPER_ANALYZER),
                "adapter": fast_acceptance.binding(
                    fast_acceptance.FAST_REPORT_UTILITY
                ),
            },
        },
    )
    return FakeReport(windows, ensemble, records)


def add_strict_failure_fixture(
    root: Path,
    lineage: dict[str, object],
    *,
    job_id: str,
    failure_time: float,
    hard_bound: int,
    variant: str | None = "standard",
) -> None:
    """Retain one exact strict hard-bound failure beside a variant lineage."""

    case_id = str(lineage["case_id"])
    segment = root / f"runs/strict/{case_id}/fast_s000_t0_to_t10"
    manifest = segment / "manifest/fast_run.json"
    executable = lineage["test_qualified_executable"]
    write_json(
        manifest,
        {
            "schema_version": 1,
            "root": str(root.absolute()),
            "case_id": case_id,
            "case_name": matrix_cases()[case_id]["name"],
            "job_id": job_id,
            "run_dir": str(segment.absolute()),
            "output_dir": str((segment / "output").absolute()),
            "sequence": 0,
            "start_time": 0.0,
            "target_time": 10.0,
            "slurm_log": str(root.absolute() / "logs/slurm-fast/%x.%j.log"),
            "variant": variant,
            "command_line_overrides": [],
            "input": lineage["input"]["path"],
            "executable": executable["path"],
            "input_sha256": lineage["lineage_identities"]["input_sha256"][0],
            "matrix_sha256": lineage["lineage_identities"]["matrix_sha256"][0],
            "executable_sha256": executable["sha256"],
        },
    )
    (segment / "manifest/run_exit_code").write_text("143\n", encoding="utf-8")
    log = root / f"logs/slurm-fast/cglf_{case_id}_s000.{job_id}.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text(
        f"elapsed=1 cycle=10030 time={failure_time:.17g} dt=1e-4\n"
        "CGL Landau-fluid strict admissibility failed after a split stage: "
        "sweep=post stage=1/19 dfloor=0 pfloor=0 nonfinite=0 "
        f"nonpositive=0 hard_bound={hard_bound}\n",
        encoding="utf-8",
    )
    lineage["unselected_lineages"] = [{
        "segments": [{
            "segment": str(segment.absolute()),
            "variant": variant,
            "command_line_overrides": [],
            "job_id": job_id,
            "state": "failed",
        }],
    }]


def add_r15_strict_failure_fixture(
    root: Path,
    lineage: dict[str, object],
    *,
    variant: str | None = "standard",
) -> None:
    """Retain the observed strict R15 failure beside a variant lineage."""

    add_strict_failure_fixture(
        root,
        lineage,
        job_id="4771183",
        failure_time=1.275643,
        hard_bound=2,
        variant=variant,
    )


def summarize(
    fast_acceptance,
    tmp_path: Path,
    case_id: str,
    lineage: dict[str, object],
    *,
    policy: dict[str, object] | None = None,
    complete_diagnostics: bool = False,
    monkeypatch=None,
) -> tuple[dict[str, object], dict[str, object] | None]:
    selected_policy = policy or policy_fixture([case_id], lineage)
    case_dir = tmp_path / "report" / "cases" / case_id
    lineage_path = case_dir / "lineage.json"
    if complete_diagnostics:
        add_authenticated_execution_fixture(
            fast_acceptance,
            tmp_path / "execution-fixture",
            selected_policy,
            lineage,
        )
    write_json(lineage_path, lineage)
    if complete_diagnostics:
        assert monkeypatch is not None
        report = write_complete_diagnostics(fast_acceptance, case_dir, lineage)
        monkeypatch.setattr(fast_acceptance, "load_report_module", lambda: report)
    return fast_acceptance.summarize_case(
        FakeAcceptance(selected_policy),
        selected_policy,
        case_id,
        lineage,
        fast_acceptance.binding(lineage_path),
        "fixture",
        tmp_path / "acceptance" / "cases" / case_id,
    )


def test_wrong_case_identity_is_rejected(fast_acceptance, tmp_path):
    lineage = lineage_fixture(tmp_path / "lineage", case_id="R02")

    summary, reviewed = summarize(fast_acceptance, tmp_path, "R17", lineage)

    assert summary["result"] != "pass"
    assert summary["health"]["result"] == "fail"
    assert any(
        "case_id differs" in error
        for error in summary["health"]["structural_errors"]
    )
    assert reviewed is None


@pytest.mark.parametrize(
    ("case_id", "strict", "variants", "overrides"),
    [
        ("R03", False, ["standard"], ["mhd/cgl_lf_strict_admissibility=false"]),
        (
            "R14",
            True,
            ["finite_limiter_hard_bound_diagnostic_nonfatal"],
            [
                "mhd/cgl_lf_strict_admissibility=false",
                "mhd/cgl_lf_strict_admissibility=true",
            ],
        ),
    ],
)
def test_strict_policy_and_configuration_mismatches_are_ineligible(
    fast_acceptance,
    tmp_path,
    case_id,
    strict,
    variants,
    overrides,
):
    lineage = lineage_fixture(
        tmp_path / case_id,
        case_id=case_id,
        strict=strict,
        variants=variants,
        overrides=overrides,
    )

    errors = fast_acceptance.lineage_identity_errors(
        FakeAcceptance(policy_fixture([case_id])),
        policy_fixture([case_id]),
        case_id,
        lineage,
        {"path": "fixture-lineage.json", "sha256": "a" * 64},
    )

    assert errors
    assert any(
        "strict-admissibility" in error or "override" in error
        for error in errors
    )
    if case_id == "R14":
        scope = fast_acceptance.scientific_scope(case_id, lineage)
        assert scope["campaign_interpretation_eligible"] is False
        assert scope["uniform_strict_diagnostics_eligible"] is False


def test_missing_strict_model_choice_legacy_exception_is_r02_only(
    fast_acceptance, tmp_path
):
    policy = policy_fixture(["R02", "R03"])
    acceptance = FakeAcceptance(policy)
    r02 = lineage_fixture(tmp_path / "R02", case_id="R02")
    r03 = lineage_fixture(tmp_path / "R03", case_id="R03")
    del r02["model_choices"]["cgl_lf_strict_admissibility"]
    del r03["model_choices"]["cgl_lf_strict_admissibility"]

    r02_errors = fast_acceptance.lineage_identity_errors(
        acceptance, policy, "R02", r02, {"path": "R02-lineage.json"}
    )
    r03_errors = fast_acceptance.lineage_identity_errors(
        acceptance, policy, "R03", r03, {"path": "R03-lineage.json"}
    )

    assert not r02_errors
    assert any("model choices" in error for error in r03_errors)


def test_spoofed_model_choices_and_forcing_tcorr_are_rejected(
    fast_acceptance, tmp_path
):
    lineage = lineage_fixture(tmp_path / "lineage", case_id="R03")
    lineage["model_choices"]["forcing_tcorr"] = "0.000001"
    lineage["model_choices"]["beta0"] = "999"

    summary, comparison = summarize(fast_acceptance, tmp_path, "R03", lineage)

    assert summary["health"]["result"] == "fail"
    assert comparison is None
    assert any(
        "model choices differ" in error
        for error in summary["health"]["structural_errors"]
    )


def test_bound_input_must_match_matrix_authoritative_input(
    fast_acceptance, tmp_path
):
    lineage = lineage_fixture(tmp_path / "lineage", case_id="R03")
    original = Path(str(lineage["input"]["path"]))
    attacker_input = tmp_path / "attacker.athinput"
    attacker_input.write_text(
        original.read_text(encoding="utf-8").replace(
            "tcorr = 2.0", "tcorr = 0.000001"
        ),
        encoding="utf-8",
    )
    lineage["input"] = {
        "path": str(attacker_input.resolve()),
        "size_bytes": attacker_input.stat().st_size,
        "sha256": sha256(attacker_input),
    }
    lineage["lineage_identities"]["input_sha256"] = [sha256(attacker_input)]
    lineage["model_choices"] = expected_model_choices(attacker_input, [])

    summary, comparison = summarize(fast_acceptance, tmp_path, "R03", lineage)

    assert summary["health"]["result"] == "fail"
    assert comparison is None
    assert any(
        "matrix-authoritative input" in error
        for error in summary["health"]["structural_errors"]
    )


def test_declared_history_binding_mismatch_is_rejected(fast_acceptance, tmp_path):
    lineage = lineage_fixture(tmp_path / "lineage")
    lineage["histories"]["mhd"]["binding"]["sha256"] = "0" * 64
    reviewed = FakeAcceptance(policy_fixture())

    histories, _bindings, errors = fast_acceptance.load_histories(reviewed, lineage)

    assert "mhd" not in histories
    assert any("sha256" in error or "binding" in error for error in errors)


def test_standard_r15_is_not_admitted_to_campaign_claim_scope(
    fast_acceptance, tmp_path
):
    lineage = lineage_fixture(tmp_path / "lineage", case_id="R15")
    policy = policy_fixture(["R15"])

    errors = fast_acceptance.lineage_identity_errors(
        FakeAcceptance(policy),
        policy,
        "R15",
        lineage,
        {"path": "fixture-lineage.json"},
    )
    scope = fast_acceptance.scientific_scope("R15", lineage)

    assert any("exact admitted variant" in error for error in errors)
    assert scope["classification"] == "r15_variant_not_authenticated"
    assert scope["campaign_interpretation_eligible"] is False


def test_qualified_executable_mismatch_blocks_execution_lineage(
    fast_acceptance, tmp_path
):
    lineage = lineage_fixture(tmp_path / "lineage", case_id="R03")
    policy = policy_fixture(["R03"], lineage)
    add_authenticated_execution_fixture(
        fast_acceptance, tmp_path / "execution-fixture", policy, lineage
    )
    manifest_path = Path(str(lineage["lineage"][0]["manifest"]["path"]))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["executable_sha256"] = "0" * 64
    write_json(manifest_path, manifest)
    lineage["lineage"][0]["manifest"] = plain_binding(manifest_path)
    acceptance = FakeAcceptance(policy)
    _histories, bindings, errors = fast_acceptance.load_histories(acceptance, lineage)
    assert not errors

    with pytest.raises(fast_acceptance.FastAcceptanceError, match="manifest differs"):
        fast_acceptance.direct_fast_execution_lineage_evidence(
            acceptance, policy, "R03", lineage, bindings
        )


def test_merged_history_must_match_exact_source_replay(fast_acceptance, tmp_path):
    lineage = lineage_fixture(tmp_path / "lineage", case_id="R03")
    policy = policy_fixture(["R03"], lineage)
    add_authenticated_execution_fixture(
        fast_acceptance, tmp_path / "execution-fixture", policy, lineage
    )
    merged = Path(str(lineage["histories"]["user"]["path"]))
    merged.write_text(merged.read_text(encoding="utf-8") + "# attacker row\n")
    lineage["histories"]["user"]["binding"] = plain_binding(merged)
    acceptance = FakeAcceptance(policy)
    _histories, bindings, errors = fast_acceptance.load_histories(acceptance, lineage)
    assert not errors
    report = FakeReport({}, {})

    original = fast_acceptance.load_report_module
    fast_acceptance.load_report_module = lambda: report
    try:
        with pytest.raises(fast_acceptance.FastAcceptanceError, match="exact replay"):
            fast_acceptance.direct_fast_execution_lineage_evidence(
                acceptance, policy, "R03", lineage, bindings
            )
    finally:
        fast_acceptance.load_report_module = original


def test_failed_complete_cases_are_excluded_from_reviewed_comparisons(
    fast_acceptance, tmp_path
):
    policy = policy_fixture(["R14", "R15"])
    acceptance = FakeAcceptance(policy)
    lineage = lineage_fixture(
        tmp_path / "R14",
        case_id="R14",
        strict=False,
        variants=["finite_limiter_hard_bound_diagnostic_nonfatal"],
        overrides=["mhd/cgl_lf_strict_admissibility=false"],
    )
    histories, _bindings, errors = fast_acceptance.load_histories(acceptance, lineage)
    scope = fast_acceptance.scientific_scope("R14", lineage)
    failed_health = fast_acceptance.numerical_health(
        "R14",
        lineage,
        scope,
        histories,
        errors,
        ["fixture provenance failure"],
    )
    reviewed, _reviewed_reason = fast_acceptance.reviewed_complete_case_evidence(
        acceptance, policy, "R14", lineage, failed_health
    )
    comparison, _comparison_reason = fast_acceptance.comparison_case_evidence(
        acceptance,
        policy,
        "R14",
        lineage,
        scope,
        failed_health,
        histories,
    )
    failed_summary = {
        "result": "fail",
        "health": failed_health,
        "scope": scope,
    }
    passed_summary = {
        "result": "pass",
        "health": {"result": "pass", "complete_to_target": True},
        "scope": {
            "classification": "standard_claim_scope",
            "campaign_interpretation_eligible": True,
        },
    }
    assert reviewed is None
    assert comparison is None

    campaign = fast_acceptance.build_campaign_evidence(
        acceptance,
        policy,
        {"R14": failed_summary, "R15": passed_summary},
        {"R15": {"result": "pass"}},
    )
    gate_by_name = {gate["name"]: gate for gate in campaign["gates"]}

    reviewed_gate = gate_by_name["scoped_exact_window_comparison_evidence"]
    limiter_gate = gate_by_name["finite_limiter_ordering:R15_gt_R14"]
    assert "R14" not in reviewed_gate["observations"]["available_cases"]
    assert limiter_gate["result"] == "inconclusive"
    assert "R14" not in campaign["claim_scope"]["claim_grade_cases"]
    assert gate_by_name["explicit_claim_scope"]["result"] == "fail"


def test_comparison_gates_require_pass_case_comparison_evidence(fast_acceptance):
    policy = policy_fixture(["R02", "R06", "R14", "R15"])
    policy["criteria"]["family_gates"]["active_passive"]["pairs"] = [["R02", "R06"]]
    acceptance = FakeAcceptance(policy)
    standard_scope = {
        "classification": "standard_claim_scope",
        "campaign_interpretation_eligible": True,
    }
    r14_scope = {
        "classification": "scoped_nonfatal_hard_bound_variant",
        "campaign_interpretation_eligible": True,
    }
    summaries = {
        case_id: {
            "result": "inconclusive",
            "health": {"result": "pass", "complete_to_target": True},
            "scope": r14_scope if case_id == "R14" else standard_scope,
        }
        for case_id in ("R02", "R06", "R14", "R15")
    }
    comparison_cases = {
        case_id: {"case_id": case_id, "result": "inconclusive"}
        for case_id in summaries
    }

    campaign = fast_acceptance.build_campaign_evidence(
        acceptance, policy, summaries, comparison_cases
    )
    gates = {gate["name"]: gate for gate in campaign["gates"]}

    assert gates["active_passive_pair:R02:R06"]["result"] == "inconclusive"
    assert gates["finite_limiter_ordering:R15_gt_R14"]["result"] == "inconclusive"


def test_complete_case_statistics_use_authenticated_forcing_tcorr(
    fast_acceptance, tmp_path
):
    lineage = lineage_fixture(tmp_path / "lineage")

    summary, comparison = summarize(fast_acceptance, tmp_path, "R02", lineage)

    assert comparison is None
    assert summary["comparison_evidence"]["result"] == "inconclusive"
    statistics = summary["history_statistics"]["kinetic"]["windows"]["full"][
        "statistics"
    ]
    assert statistics["method"]["required_block_duration"] == pytest.approx(2.0)


def test_complete_direct_fast_science_is_claim_grade_without_canonical_package(
    fast_acceptance, tmp_path, monkeypatch
):
    lineage = lineage_fixture(tmp_path / "lineage", case_id="R02")

    summary, comparison = summarize(
        fast_acceptance,
        tmp_path,
        "R02",
        lineage,
        complete_diagnostics=True,
        monkeypatch=monkeypatch,
    )

    evidence = json.loads(
        (
            tmp_path
            / "acceptance/cases/R02/reviewed_case_evidence.json"
        ).read_text(encoding="utf-8")
    )
    assert summary["result"] == "pass"
    assert comparison is not None and comparison["result"] == "pass"
    assert evidence["result"] == "pass"
    assert evidence["authority"] == "non-authorizing-direct-fast-scientific-assessment"
    assert evidence["release_authorizing"] is False
    assert evidence["campaign_authority_eligible"] is False
    assert evidence["evaluation_inputs"]["accepted_bundle_manifest"]["sha256"] == sha256(
        tmp_path / "execution-fixture/accepted-bundle/manifest.json"
    )
    assert evidence["evaluation_inputs"]["diagnostics"]["direct_fast_complete"][
        "sha256"
    ] == sha256(tmp_path / "report/cases/R02/diagnostics.json")
    gate_results = {gate["name"]: gate["result"] for gate in evidence["gates"]}
    assert gate_results["direct_fast_execution_lineage"] == "pass"
    assert gate_results["direct_fast_complete_diagnostics"] == "pass"
    assert gate_results["sampled_restart_ct_divb"] == "blocked_out_of_scope"
    assert gate_results["canonical_comparison_panel_products"] == "blocked_out_of_scope"
    assert "sampled-restart CT divergence evidence" in evidence[
        "scientific_kernel_scope"
    ]["excluded_from_direct_fast_assessment"]


def test_direct_fast_claim_grade_uses_reviewed_active_energy_gate(
    fast_acceptance, tmp_path, monkeypatch
):
    lineage = lineage_fixture(tmp_path / "lineage", case_id="R02")
    policy = policy_fixture(["R02"])
    policy["criteria"]["active_energy_policy"] = {"revision_id": "fixture"}
    policy["criteria"]["family_gates"]["active_passive"].update({
        "active_cases": ["R02"],
        "passive_cases": [],
        "activity_absolute_gt": 0.0,
    })

    summary, _comparison = summarize(
        fast_acceptance,
        tmp_path,
        "R02",
        lineage,
        policy=policy,
        complete_diagnostics=True,
        monkeypatch=monkeypatch,
    )
    evidence = json.loads(
        (
            tmp_path
            / "acceptance/cases/R02/reviewed_case_evidence.json"
        ).read_text(encoding="utf-8")
    )

    assert summary["result"] == "pass"
    assert {
        gate["name"]: gate["result"] for gate in evidence["gates"]
    }["active_energy_closure"] == "pass"


def test_stale_direct_fast_diagnostics_cannot_become_claim_grade(
    fast_acceptance, tmp_path, monkeypatch
):
    lineage = lineage_fixture(tmp_path / "lineage", case_id="R02")
    policy = policy_fixture(["R02"], lineage)
    case_dir = tmp_path / "report/cases/R02"
    add_authenticated_execution_fixture(
        fast_acceptance, tmp_path / "execution-fixture", policy, lineage
    )
    write_json(case_dir / "lineage.json", lineage)
    report = write_complete_diagnostics(fast_acceptance, case_dir, lineage)
    monkeypatch.setattr(fast_acceptance, "load_report_module", lambda: report)
    diagnostics_path = case_dir / "diagnostics.json"
    diagnostics = json.loads(diagnostics_path.read_text(encoding="utf-8"))
    diagnostics["provenance"]["merged_user_history"]["sha256"] = "0" * 64
    write_json(diagnostics_path, diagnostics)

    summary, comparison = fast_acceptance.summarize_case(
        FakeAcceptance(policy),
        policy,
        "R02",
        lineage,
        fast_acceptance.binding(case_dir / "lineage.json"),
        "fixture",
        tmp_path / "acceptance/cases/R02",
    )

    assert summary["result"] == "inconclusive"
    assert summary["reviewed_case_evidence"]["binding"] is None
    assert "declared binding differs" in summary["reviewed_case_evidence"]["reason"]
    assert comparison is None


def test_snapshot_index_mismatch_cannot_become_claim_grade(
    fast_acceptance, tmp_path, monkeypatch
):
    lineage = lineage_fixture(tmp_path / "lineage", case_id="R02")
    policy = policy_fixture(["R02"], lineage)
    case_dir = tmp_path / "report/cases/R02"
    add_authenticated_execution_fixture(
        fast_acceptance, tmp_path / "execution-fixture", policy, lineage
    )
    write_json(case_dir / "lineage.json", lineage)
    report = write_complete_diagnostics(fast_acceptance, case_dir, lineage)
    monkeypatch.setattr(fast_acceptance, "load_report_module", lambda: report)
    snapshot_index = case_dir / "snapshots.json"
    index = json.loads(snapshot_index.read_text(encoding="utf-8"))
    index["snapshots"][0]["expected_ranks"] = 2
    write_json(snapshot_index, index)
    diagnostics_path = case_dir / "diagnostics.json"
    diagnostics = json.loads(diagnostics_path.read_text(encoding="utf-8"))
    diagnostics["provenance"]["snapshot_index"] = plain_binding(snapshot_index)
    write_json(diagnostics_path, diagnostics)

    summary, comparison = fast_acceptance.summarize_case(
        FakeAcceptance(policy),
        policy,
        "R02",
        lineage,
        plain_binding(case_dir / "lineage.json"),
        "fixture",
        tmp_path / "acceptance/cases/R02",
    )

    assert summary["result"] == "inconclusive"
    assert summary["reviewed_case_evidence"]["binding"] is None
    assert "snapshot provenance differs" in summary["reviewed_case_evidence"]["reason"]
    assert comparison is None


def test_stored_snapshot_science_must_match_byte_level_replay(
    fast_acceptance, tmp_path, monkeypatch
):
    lineage = lineage_fixture(tmp_path / "lineage", case_id="R02")
    policy = policy_fixture(["R02"], lineage)
    case_dir = tmp_path / "report/cases/R02"
    add_authenticated_execution_fixture(
        fast_acceptance, tmp_path / "execution-fixture", policy, lineage
    )
    write_json(case_dir / "lineage.json", lineage)
    report = write_complete_diagnostics(fast_acceptance, case_dir, lineage)
    monkeypatch.setattr(fast_acceptance, "load_report_module", lambda: report)
    diagnostics_path = case_dir / "diagnostics.json"
    diagnostics = json.loads(diagnostics_path.read_text(encoding="utf-8"))
    source = next(iter(diagnostics["snapshots"]))
    diagnostics["snapshots"][source]["forged_science"] = {"mean": 1.0e99}
    write_json(diagnostics_path, diagnostics)

    summary, comparison = fast_acceptance.summarize_case(
        FakeAcceptance(policy),
        policy,
        "R02",
        lineage,
        plain_binding(case_dir / "lineage.json"),
        "fixture",
        tmp_path / "acceptance/cases/R02",
    )

    assert summary["result"] == "inconclusive"
    assert summary["reviewed_case_evidence"]["binding"] is None
    assert "byte-level science replay" in summary["reviewed_case_evidence"]["reason"]
    assert comparison is None


def test_r10_remains_exploratory_with_complete_direct_fast_science(
    fast_acceptance, tmp_path, monkeypatch
):
    lineage = lineage_fixture(tmp_path / "lineage", case_id="R10")

    summary, comparison = summarize(
        fast_acceptance,
        tmp_path,
        "R10",
        lineage,
        complete_diagnostics=True,
        monkeypatch=monkeypatch,
    )
    evidence = json.loads(
        (
            tmp_path
            / "acceptance/cases/R10/reviewed_case_evidence.json"
        ).read_text(encoding="utf-8")
    )

    assert evidence["result"] == "inconclusive"
    assert summary["result"] == "inconclusive"
    assert summary["scope"]["classification"] == "exploratory_only"
    assert comparison is None


def test_exact_r15_nonfatal_variant_is_scoped_and_preserves_strict_failure(
    fast_acceptance, tmp_path, monkeypatch
):
    lineage = lineage_fixture(
        tmp_path / "lineage",
        case_id="R15",
        strict=False,
        variants=["finite_limiter_hard_bound_diagnostic_nonfatal"],
        overrides=["mhd/cgl_lf_strict_admissibility=false"],
    )
    add_r15_strict_failure_fixture(tmp_path / "campaign", lineage)

    summary, comparison = summarize(
        fast_acceptance,
        tmp_path,
        "R15",
        lineage,
        complete_diagnostics=True,
        monkeypatch=monkeypatch,
    )
    evidence = json.loads(
        (
            tmp_path
            / "acceptance/cases/R15/reviewed_case_evidence.json"
        ).read_text(encoding="utf-8")
    )

    assert summary["result"] == "pass"
    assert comparison is not None and comparison["result"] == "pass"
    assert summary["scope"]["classification"] == "scoped_nonfatal_hard_bound_variant"
    assert summary["scope"]["uniform_strict_diagnostics_eligible"] is False
    assert "uniform strict-admissibility compliance" in summary["scope"][
        "excluded_claims"
    ]
    retained = evidence["retained_strict_failure_evidence"]
    assert retained[0]["job_id"] == "4771183"
    assert retained[0]["failure_time"] == pytest.approx(1.275643)
    assert retained[0]["failure_counters"]["lf_hardbd"] == 2
    assert {
        gate["name"]: gate["result"] for gate in evidence["gates"]
    }["r15_prior_strict_failure_retained"] == "pass"


def test_exact_r14_nonfatal_variant_is_scoped_and_preserves_strict_failure(
    fast_acceptance, tmp_path, monkeypatch
):
    lineage = lineage_fixture(
        tmp_path / "lineage",
        case_id="R14",
        strict=False,
        variants=["finite_limiter_hard_bound_diagnostic_nonfatal"],
        overrides=["mhd/cgl_lf_strict_admissibility=false"],
    )
    add_strict_failure_fixture(
        tmp_path / "campaign",
        lineage,
        job_id="4770854",
        failure_time=0.1326939,
        hard_bound=167,
    )

    summary, comparison = summarize(
        fast_acceptance,
        tmp_path,
        "R14",
        lineage,
        complete_diagnostics=True,
        monkeypatch=monkeypatch,
    )
    evidence = json.loads(
        (
            tmp_path
            / "acceptance/cases/R14/reviewed_case_evidence.json"
        ).read_text(encoding="utf-8")
    )

    assert summary["result"] == "pass"
    assert comparison is not None and comparison["result"] == "pass"
    retained = evidence["retained_strict_failure_evidence"]
    assert retained[0]["job_id"] == "4770854"
    assert retained[0]["failure_time"] == pytest.approx(0.1326939)
    assert retained[0]["failure_counters"]["lf_hardbd"] == 167
    assert {
        gate["name"]: gate["result"] for gate in evidence["gates"]
    }["r14_prior_strict_failure_retained"] == "pass"


@pytest.mark.parametrize("variant", [None, "standard"])
def test_r15_strict_failure_accepts_only_exact_standard_variant_representations(
    fast_acceptance, tmp_path, variant
):
    lineage = lineage_fixture(tmp_path / "lineage", case_id="R15")
    add_r15_strict_failure_fixture(tmp_path / "campaign", lineage, variant=variant)

    records = fast_acceptance.retained_strict_failure_evidence("R15", lineage)

    assert len(records) == 1
    assert records[0]["job_id"] == "4771183"
    assert records[0]["failure_time"] == pytest.approx(1.275643)
    assert records[0]["failure_counters"]["lf_hardbd"] == 2


def test_r15_strict_failure_rejects_nonstandard_variant(
    fast_acceptance, tmp_path
):
    lineage = lineage_fixture(tmp_path / "lineage", case_id="R15")
    add_r15_strict_failure_fixture(
        tmp_path / "campaign", lineage, variant="attacker_variant"
    )

    assert fast_acceptance.retained_strict_failure_evidence("R15", lineage) == []


def test_r15_strict_failure_rejects_spoofed_executable_identity(
    fast_acceptance, tmp_path
):
    lineage = lineage_fixture(tmp_path / "lineage", case_id="R15")
    root = tmp_path / "campaign"
    add_r15_strict_failure_fixture(root, lineage)
    manifest_path = (
        root / "runs/strict/R15/fast_s000_t0_to_t10/manifest/fast_run.json"
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["executable"] = lineage["input"]["path"]
    write_json(manifest_path, manifest)

    assert fast_acceptance.retained_strict_failure_evidence("R15", lineage) == []


def test_r15_strict_failure_rejects_other_fatal_counters(
    fast_acceptance, tmp_path
):
    lineage = lineage_fixture(tmp_path / "lineage", case_id="R15")
    root = tmp_path / "campaign"
    add_r15_strict_failure_fixture(root, lineage)
    log = root / "logs/slurm-fast/cglf_R15_s000.4771183.log"
    log.write_text(
        log.read_text(encoding="utf-8").replace("dfloor=0", "dfloor=1"),
        encoding="utf-8",
    )

    assert fast_acceptance.retained_strict_failure_evidence("R15", lineage) == []


def test_r15_nonfatal_variant_without_retained_strict_failure_is_inconclusive(
    fast_acceptance, tmp_path, monkeypatch
):
    lineage = lineage_fixture(
        tmp_path / "lineage",
        case_id="R15",
        strict=False,
        variants=["finite_limiter_hard_bound_diagnostic_nonfatal"],
        overrides=["mhd/cgl_lf_strict_admissibility=false"],
    )

    summary, comparison = summarize(
        fast_acceptance,
        tmp_path,
        "R15",
        lineage,
        complete_diagnostics=True,
        monkeypatch=monkeypatch,
    )

    assert summary["result"] == "inconclusive"
    assert summary["reviewed_case_evidence"]["result"] == "inconclusive"
    assert comparison is not None and comparison["result"] == "inconclusive"


def test_rerun_removes_stale_reviewed_case_evidence(fast_acceptance, tmp_path):
    output = tmp_path / "acceptance" / "cases" / "R02"
    stale = output / "reviewed_case_evidence.json"
    write_json(stale, {"record_type": "stage-i-scientific-case-evidence"})
    assert stale.is_file()
    policy = policy_fixture()
    reviewed = FakeAcceptance(policy)
    partial = lineage_fixture(
        tmp_path / "partial",
        status="partial",
        final_time=2.0,
    )

    fast_acceptance.summarize_case(
        reviewed,
        policy,
        "R02",
        partial,
        {"path": "fixture-lineage.json", "sha256": "a" * 64},
        "fixture",
        output,
    )

    assert not stale.exists()


def test_missing_histories_are_inconclusive_for_in_progress_case(
    fast_acceptance,
):
    scope = {
        "hard_bound_is_fatal": True,
        "campaign_interpretation_eligible": True,
    }
    lineage = {"status": "partial", "target_time": 10.0, "errors": []}

    health = fast_acceptance.numerical_health(
        "R02",
        lineage,
        scope,
        {},
        [
            "evidence gap: assembled mhd history is unavailable",
            "evidence gap: assembled user history is unavailable",
        ],
        [],
    )

    assert health["result"] == "inconclusive"


def test_missing_fatal_counters_are_inconclusive_for_complete_case(
    fast_acceptance, tmp_path
):
    lineage = lineage_fixture(tmp_path / "lineage", counters=False)
    reviewed = FakeAcceptance(policy_fixture())
    histories, _bindings, errors = fast_acceptance.load_histories(reviewed, lineage)
    scope = fast_acceptance.scientific_scope("R02", lineage)

    health = fast_acceptance.numerical_health(
        "R02", lineage, scope, histories, errors, []
    )

    assert health["result"] == "inconclusive"
    assert all(
        health["fatal_counter_maxima"][name] is None
        for name in fast_acceptance.FATAL_COUNTERS
    )


def test_symlink_output_cannot_write_inside_fast_report_source(
    fast_acceptance, tmp_path, monkeypatch
):
    source = tmp_path / "fast-report"
    source.mkdir()
    inventory = source / "inventory.json"
    lineage = lineage_fixture(source / "histories")
    write_json(
        inventory,
        {
            "output": str(source),
            "adapter": fast_acceptance.binding(fast_acceptance.FAST_REPORT_UTILITY),
            "cases": {"R02": lineage},
        },
    )
    output_link = tmp_path / "acceptance-output"
    output_link.symlink_to(source / "nested-acceptance", target_is_directory=True)
    policy = policy_fixture()
    monkeypatch.setattr(
        fast_acceptance,
        "load_acceptance_module",
        lambda: FakeAcceptance(policy),
    )

    with pytest.raises(SystemExit):
        fast_acceptance.main(
            [
                "--inventory",
                str(inventory),
                "--output",
                str(output_link),
                "--cases",
                "R02",
            ]
        )

    assert not (source / "nested-acceptance").exists()


def test_inventory_reporter_binding_must_match_current_adapter(
    fast_acceptance,
):
    current = fast_acceptance.binding(fast_acceptance.FAST_REPORT_UTILITY)

    assert fast_acceptance.authenticate_inventory_reporter(
        {"adapter": current}
    ) == current
    with pytest.raises(
        fast_acceptance.FastAcceptanceError,
        match="inventory reporter adapter declared binding differs",
    ):
        fast_acceptance.authenticate_inventory_reporter({
            "adapter": {**current, "sha256": "0" * 64},
        })


def test_stale_inventory_reporter_fails_before_claim_outputs(
    fast_acceptance, tmp_path
):
    source = tmp_path / "fast-report"
    source.mkdir()
    current = fast_acceptance.binding(fast_acceptance.FAST_REPORT_UTILITY)
    inventory = source / "inventory.json"
    write_json(inventory, {
        "output": str(source),
        "adapter": {**current, "sha256": "0" * 64},
        "cases": {},
    })
    output = tmp_path / "acceptance-output"

    with pytest.raises(SystemExit):
        fast_acceptance.main([
            "--inventory",
            str(inventory),
            "--output",
            str(output),
            "--cases",
            "R02",
        ])

    assert not output.exists()


def test_case_lineage_must_match_authenticated_inventory_record(
    fast_acceptance, tmp_path
):
    output = tmp_path / "fast-report"
    lineage = lineage_fixture(tmp_path / "lineage")
    write_json(output / "cases/R02/lineage.json", {**lineage, "status": "partial"})
    inventory = {"output": str(output), "cases": {"R02": lineage}}

    with pytest.raises(
        fast_acceptance.FastAcceptanceError,
        match="differs from its authenticated inventory record",
    ):
        fast_acceptance.case_record("R02", output, inventory)


def test_case_lineage_parse_and_binding_do_not_use_split_reads(
    fast_acceptance, tmp_path, monkeypatch
):
    output = tmp_path / "fast-report"
    lineage = lineage_fixture(tmp_path / "lineage")
    lineage_path = output / "cases/R02/lineage.json"
    write_json(lineage_path, lineage)
    inventory = {"output": str(output), "cases": {"R02": lineage}}
    original_binding = fast_acceptance.binding

    def mutate_on_split_binding(path: Path) -> dict[str, object]:
        if path.resolve() == lineage_path.resolve():
            write_json(lineage_path, {**lineage, "status": "mutated_after_parse"})
        return original_binding(path)

    monkeypatch.setattr(fast_acceptance, "binding", mutate_on_split_binding)
    record, record_binding, source = fast_acceptance.case_record(
        "R02", output, inventory
    )

    assert record == lineage
    assert source == "case_directory"
    assert record_binding == plain_binding(lineage_path)
    assert json.loads(lineage_path.read_text(encoding="utf-8")) == lineage


def test_nested_symlink_output_escape_is_rejected_before_writing(
    fast_acceptance, tmp_path, monkeypatch
):
    source = tmp_path / "fast-report"
    source.mkdir()
    inventory = source / "inventory.json"
    lineage = lineage_fixture(source / "histories")
    write_json(inventory, {
        "output": str(source),
        "adapter": fast_acceptance.binding(fast_acceptance.FAST_REPORT_UTILITY),
        "cases": {"R02": lineage},
    })
    escaped = source / "escaped-acceptance"
    escaped.mkdir()
    output = tmp_path / "acceptance-output"
    (output / "cases").mkdir(parents=True)
    (output / "cases" / "R02").symlink_to(escaped, target_is_directory=True)
    policy = policy_fixture()
    monkeypatch.setattr(
        fast_acceptance,
        "load_acceptance_module",
        lambda: FakeAcceptance(policy),
    )

    with pytest.raises(SystemExit):
        fast_acceptance.main(
            [
                "--inventory",
                str(inventory),
                "--output",
                str(output),
                "--cases",
                "R02",
            ]
        )

    assert not (escaped / "case_acceptance.json").exists()


def test_publication_discovers_scoped_direct_fast_campaign_record(
    publication, tmp_path
):
    analysis = tmp_path / "analysis"
    analysis.mkdir()
    acceptance = tmp_path / "acceptance"
    write_json(
        acceptance / "campaign_evidence.json",
        {
            "record_type": "cgl-lf-stage-i-direct-fast-campaign-evidence",
            "result": "inconclusive",
            "case_results": {"R10": "inconclusive", "R14": "pass"},
            "claim_scope": {
                "R10": {"classification": "exploratory_only"},
                "R14": {
                    "classification": "scoped_nonfatal_hard_bound_variant"
                },
            },
            "gates": [
                {
                    "name": "explicit_claim_scope",
                    "result": "pass",
                    "reason": "fixture",
                }
            ],
        },
    )

    data = publication.discover_data(analysis, [acceptance])

    assert data.campaign_acceptance is not None
    assert (
        data.campaign_acceptance["record_type"]
        == "cgl-lf-stage-i-direct-fast-campaign-evidence"
    )
    assert any(
        record.get("record_type")
        == "cgl-lf-stage-i-direct-fast-campaign-evidence"
        for record in data.acceptance_records
    )
