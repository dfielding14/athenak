"""Adversarial regressions for isolated Stage I qualification."""

from __future__ import annotations

from copy import deepcopy
from contextlib import contextmanager
from datetime import datetime, timedelta, timezone
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shlex
import struct
import subprocess

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
UTILITY = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_qualification.py"


def load_utility():
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_stage_i_qualification", UTILITY
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


qualification = load_utility()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def stat_mode(path: Path) -> int:
    return path.stat().st_mode & 0o777


def run_git(repository: Path, *arguments: str) -> str:
    return subprocess.run(
        ["git", "-C", str(repository), *arguments],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def assert_hardened_scheduler_environment(environment: dict[str, str], *,
                                          exact_timestamps: bool = False) -> None:
    assert environment == qualification.scheduler_environment(
        exact_timestamps=exact_timestamps
    )


def render_parameters(values: dict[str, str]) -> str:
    blocks: dict[str, list[tuple[str, str]]] = {}
    for qualified, value in values.items():
        block, key = qualified.split("/", 1)
        blocks.setdefault(block, []).append((key, value))
    lines = []
    for block in sorted(blocks):
        lines.append(f"<{block}>")
        lines.extend(
            f"{key} = {value}" for key, value in sorted(blocks[block])
        )
        lines.append("")
    return "\n".join(lines)


@pytest.fixture
def qualification_fixture(tmp_path, monkeypatch):
    root = tmp_path / "project"
    runs = root / "runs"
    source_archives = root / "source-archives"
    build = root / "build"
    build_manifest = root / "build-manifest"
    accounting = root / "accounting"
    for directory in (runs, source_archives, build, build_manifest, accounting):
        directory.mkdir(parents=True)
    (root / qualification.STAGE_I_LOCK_NAME).write_bytes(b"")
    (root / qualification.STAGE_I_LOCK_NAME).chmod(0o644)
    write_json(accounting / qualification.STAGE_I_RESERVATIONS_NAME, [])
    for name in qualification.STAGE_I_TRANSACTION_NAMES:
        (accounting / name).mkdir()

    policies = deepcopy(qualification.CASE_POLICIES)
    for case_id, policy in policies.items():
        if case_id == qualification.R17_CASE_ID:
            policy["mesh_shape"] = (12, 12, 12)
            policy["meshblock_shape"] = (1, 1, 1)
        else:
            policy["mesh_shape"] = (4, 1, 1)
            policy["meshblock_shape"] = (1, 1, 1)
        policy["minimum_restart_bytes"] = 512
    r12_profile = deepcopy(qualification.R12_FRESH_RERUN_PROFILE)
    r12_profile["ranks_per_node"] = 1
    r12_profile["total_ranks"] = r12_profile["nodes"]
    policies[qualification.R12_CASE_ID]["fixed_fresh_profile"] = r12_profile
    monkeypatch.setattr(qualification, "CASE_POLICIES", policies)
    monkeypatch.setattr(qualification, "RANKS_PER_NODE", 1)
    monkeypatch.setattr(qualification, "R12_FRESH_RERUN_PROFILE", r12_profile)
    monkeypatch.setattr(qualification, "CPUS_PER_TASK", 1)

    source = tmp_path / "source"
    helper = source / "scripts/frontier/cgl_lf_stage_i_qualification.py"
    helper.parent.mkdir(parents=True)
    helper.write_bytes(UTILITY.read_bytes())
    helper.chmod(0o755)
    cases = []
    selected_paths = []
    for number in range(2, 18):
        case_id = f"R{number:02d}"
        if case_id in policies:
            policy = policies[case_id]
            relative = Path(str(policy["input_relative_path"]))
            input_path = source / relative
            input_path.parent.mkdir(parents=True, exist_ok=True)
            values = qualification.expected_case_configuration(
                case_id,
                basename=str(policy["input_basename"]),
                tlim=qualification.COMMON_INPUT_CONTRACT["time/tlim"],
            )
            input_path.write_text(render_parameters(values))
            policy["input_sha256"] = hashlib.sha256(input_path.read_bytes()).hexdigest()
            name = policy["case_name"]
            resolution = policy["resolution"]
            selected_paths.append(input_path)
        else:
            relative = Path(f"inputs/{case_id}.athinput")
            name = f"fixture_{case_id}"
            resolution = "fixture"
        cases.append({
            "id": case_id,
            "name": name,
            "input": str(relative),
            "resolution": resolution,
            "figure_roles": [],
        })
    matrix_path = source / "inputs/cgl_lf_paper/mks24_stage_i_manifest.json"
    write_json(matrix_path, {"cases": cases})
    monkeypatch.setattr(
        qualification,
        "FROZEN_MATRIX_SHA256",
        hashlib.sha256(matrix_path.read_bytes()).hexdigest(),
    )

    run_git(source.parent, "init", str(source))
    run_git(source, "config", "user.email", "qualification@example.invalid")
    run_git(source, "config", "user.name", "Qualification Tests")
    run_git(source, "add", ".")
    run_git(source, "commit", "-m", "fixture source")
    revision = run_git(source, "rev-parse", "HEAD")

    monkeypatch.setattr(qualification, "DEFAULT_ROOT", root.resolve())
    monkeypatch.setattr(qualification, "REPOSITORY_ROOT", source.resolve())
    monkeypatch.setattr(
        qualification,
        "UTILITY_RELATIVE_PATH",
        Path("scripts/frontier/cgl_lf_stage_i_qualification.py"),
    )
    monkeypatch.setattr(qualification, "__file__", str(helper.resolve()))
    monkeypatch.setattr(qualification, "FROZEN_SOURCE_REVISION", revision)

    executable = build / "athena"
    executable.write_bytes(b"fixture executable\n")
    executable.chmod(0o755)
    executable_sha256 = hashlib.sha256(executable.read_bytes()).hexdigest()
    (build_manifest / "athena.sha256").write_text(
        f"{executable_sha256}  athena\n"
    )
    (build_manifest / "environment.txt").write_text(
        f"git_revision={revision}\n"
    )
    bundle = source_archives / "source.bundle"
    run_git(source, "bundle", "create", str(bundle), "--all")

    source_revision = qualification.committed_source_revision(
        source, [matrix_path, *selected_paths]
    )
    helper_provenance = qualification.utility_provenance()
    build_provenance = qualification.build_provenance(executable, build_manifest)
    test_abis = dict(qualification.QUALIFIED_RESTART_BINARY_ABIS)
    test_abis[(
        build_provenance["executable"]["revision"],
        build_provenance["executable"]["sha256"],
    )] = {
        "mesh_time_offset_after_parameter_dump": 232,
        "mesh_time_format": "<d",
        "allowed_marker_modes": frozenset({
            "full_precision", "legacy_default_precision",
        }),
    }
    monkeypatch.setattr(qualification, "QUALIFIED_RESTART_BINARY_ABIS", test_abis)
    bundle_provenance = qualification.source_bundle_provenance(
        bundle,
        root,
        [
            source_revision,
            helper_provenance["revision"],
            build_provenance["executable"]["revision"],
        ],
    )
    provenance = {
        "source": {
            "directory": str(source.resolve()),
            "revision": source_revision,
        },
        "source_bundle": bundle_provenance,
        "matrix": qualification.retained_file(matrix_path, "matrix"),
        **build_provenance,
        "qualification_helper": helper_provenance,
    }

    def build_wave(case_ids=("R04",), max_nodes=10):
        _, selected = qualification.load_case_records(
            matrix_path, source, list(case_ids)
        )
        return qualification.build_prepared_wave(
            root.resolve(),
            (runs / "qualification-focused-tests").resolve(),
            selected,
            provenance,
            max_nodes,
        )

    return {
        "root": root.resolve(),
        "source": source.resolve(),
        "matrix": matrix_path.resolve(),
        "bundle": bundle.resolve(),
        "helper": helper.resolve(),
        "executable": executable.resolve(),
        "build_manifest": build_manifest.resolve(),
        "provenance": provenance,
        "build_wave": build_wave,
    }


def packet_map(wave: dict[str, object]) -> dict[tuple[str, int], dict[str, object]]:
    return {
        (
            packet["execution_intent"]["case_id"],
            packet["execution_intent"]["allocation"]["nodes"],
        ): packet
        for wave_record in wave["waves"]
        for packet in wave_record["packets"]
    }


class EmptyQueue:
    def __init__(self, stdout=""):
        self.initial_stdout = stdout
        self.submitted = {}
        self.calls = []

    def add_submission(self, job_id, script):
        text = script.decode()
        name = re.search(r"(?m)^#SBATCH -J (.+)$", text).group(1)
        nodes = int(re.search(r"(?m)^#SBATCH -N ([0-9]+)$", text).group(1))
        self.submitted[job_id] = (
            f"{job_id}|{job_id}|N/A|{name}|PENDING|{nodes}|producer\n"
        )

    def remove_submission(self, job_id):
        self.submitted.pop(job_id, None)

    def __call__(self, argv, **kwargs):
        self.calls.append((list(argv), kwargs))
        assert argv[:3] == [str(qualification.SQUEUE), "-A", qualification.ACCOUNT]
        assert "--array" in argv
        assert {
            key: value for key, value in kwargs.items() if key != "env"
        } == {"check": True, "capture_output": True, "text": True}
        assert_hardened_scheduler_environment(kwargs["env"])
        return subprocess.CompletedProcess(
            argv, 0,
            stdout=self.initial_stdout + "".join(self.submitted.values()),
            stderr="",
        )


class PostSubmissionRaceQueue(EmptyQueue):
    def __call__(self, argv, **kwargs):
        completed = super().__call__(argv, **kwargs)
        if self.submitted:
            completed.stdout += (
                "990|990|N/A|unrelated|RUNNING|10|producer\n"
            )
        return completed


class SbatchSubmitter:
    def __init__(self, start=1000, hook=None, queue=None):
        self.next_job_id = start
        self.calls = []
        self.hook = hook
        self.queue = queue
        self.scripts = {}

    def __call__(self, argv, **kwargs):
        script = kwargs["input"]
        assert isinstance(script, bytes)
        assert kwargs["check"] is True and kwargs["capture_output"] is True
        assert "text" not in kwargs
        assert_hardened_scheduler_environment(kwargs["env"])
        assert all(not argument.endswith(".sbatch") for argument in argv)
        self.calls.append((list(argv), script))
        if self.hook is not None:
            self.hook(len(self.calls), argv, script)
        job_id = self.next_job_id
        self.next_job_id += 1
        self.scripts[str(job_id)] = script
        if self.queue is not None:
            self.queue.add_submission(str(job_id), script)
        return subprocess.CompletedProcess(
            argv, 0, stdout=f"{job_id};frontier\n".encode(), stderr=b""
        )

    def read_batch_script(self, argv, **kwargs):
        assert argv[:3] == [
            str(qualification.SCONTROL), "write", "batch_script",
        ]
        assert argv[4] == "-"
        assert {
            key: value for key, value in kwargs.items() if key != "env"
        } == {"check": True, "capture_output": True}
        assert_hardened_scheduler_environment(kwargs["env"])
        return subprocess.CompletedProcess(
            argv, 0, stdout=self.scripts[argv[3]], stderr=b""
        )


class AmbiguousSubmitter(SbatchSubmitter):
    def __call__(self, argv, **kwargs):
        super().__call__(argv, **kwargs)
        raise subprocess.TimeoutExpired(argv, 30)


class ErrorAfterIdSubmitter(SbatchSubmitter):
    def __call__(self, argv, **kwargs):
        completed = super().__call__(argv, **kwargs)
        raise subprocess.TimeoutExpired(
            argv, 30, output=completed.stdout, stderr=completed.stderr
        )


class CancelRunner:
    def __init__(self, queue):
        self.queue = queue
        self.calls = []

    def __call__(self, argv, **kwargs):
        assert argv[0] == str(qualification.SCANCEL)
        assert {
            key: value for key, value in kwargs.items() if key != "env"
        } == {"check": False, "capture_output": True, "text": True}
        assert_hardened_scheduler_environment(kwargs["env"])
        self.calls.append(list(argv))
        self.queue.remove_submission(argv[1])
        return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")


def submit_wave(wave, wave_path, *, queue=None, submitter=None, canceler=None):
    queue = queue or EmptyQueue()
    submitter = submitter or SbatchSubmitter(queue=queue)
    submitter.queue = queue
    canceler = canceler or CancelRunner(queue)
    result = qualification.submit_all_waves(
        wave_path,
        runner=submitter,
        queue_runner=queue,
        batch_script_runner=submitter.read_batch_script,
        cancel_runner=canceler,
    )
    return result, submitter, queue


def history_text(labels: list[str], rows: list[list[float]]) -> str:
    header = "# " + " ".join(
        f"[{index}]={name}" for index, name in enumerate(labels, start=1)
    )
    body = [" ".join(f"{value:.16e}" for value in row) for row in rows]
    return "\n".join(["# Athena++ history data", header, *body]) + "\n"


def snapshot_payload(time_value: float, logical_locations, *, intent,
                     variables=qualification.REQUIRED_SNAPSHOT_VARIABLES,
                     nonfinite=False, nonpositive=None, bad_nghost=False,
                     bad_geometry=False, bad_preheader=False,
                     bad_location_size=False, bad_variable_size=False,
                     bad_active_indices=False,
                     hard_bound_violation=False) -> bytes:
    policy = qualification.CASE_POLICIES[str(intent["case_id"])]
    parameters = qualification.expected_case_configuration(
        str(intent["case_id"]),
        basename=str(intent["run_basename"]),
        tlim=str(intent["target_time"]),
    )
    if bad_nghost:
        parameters["mesh/nghost"] = "0"
    nghost = int(parameters["mesh/nghost"])
    location_size = 4 if bad_location_size else 8
    variable_size = 8 if bad_variable_size else 4
    parameter_header = render_parameters(parameters).encode()
    preheader = [
        f"time={time_value:.17g}\n".encode(),
        b"cycle=1\n",
        f"size of location={location_size}\n".encode(),
        f"size of variable={variable_size}\n".encode(),
    ]
    if bad_preheader:
        preheader.append(b"unreviewed metadata=yes\n")
    header = (
        b"Athena binary output version=1.1\n"
        + f"size of preheader={len(preheader) + 1}\n".encode()
        + b"".join(preheader)
        + f"number of variables={len(variables)}\n".encode()
        + ("variables: " + " ".join(variables) + "\n").encode()
        + f"header offset={len(parameter_header)}\n".encode()
        + parameter_header
    )
    blocks = []
    for index, logical in enumerate(logical_locations):
        cells = math.prod(policy["meshblock_shape"])
        values = [
            float(number + 1)
            for number in range(len(variables))
            for _ in range(cells)
        ]
        if nonfinite and index == 0:
            values[0] = math.nan
        if nonpositive is not None and index == 0:
            values[list(variables).index(nonpositive) * cells] = -1.0
        if hard_bound_violation and index == 0:
            values[list(variables).index("eint") * cells] = 1.0
            values[list(variables).index("p_perp") * cells] = 3.0
            for name in ("bcc1", "bcc2", "bcc3"):
                values[list(variables).index(name) * cells] = 0.1
        geometry = []
        for axis in range(3):
            lower, upper = policy["mesh_bounds"][axis]
            count = policy["mesh_shape"][axis] // policy["meshblock_shape"][axis]
            width = (upper - lower) / count
            geometry.extend([
                lower + logical[axis] * width,
                lower + (logical[axis] + 1) * width,
            ])
        if bad_geometry and index == 0:
            geometry[1] += 0.125
        indices = [
            value
            for size in policy["meshblock_shape"]
            for value in (nghost, nghost + size - 1)
        ]
        if bad_active_indices:
            indices = [value + 1 for value in indices]
        blocks.append(
            struct.pack("<6i", *indices)
            + struct.pack("<4i", *logical)
            + struct.pack("<6" + ("f" if location_size == 4 else "d"), *geometry)
            + struct.pack(
                "<" + ("f" if variable_size == 4 else "d") * len(values), *values
            )
        )
    return header + b"".join(blocks)


def restart_payload(intent, marker_time: float, *, binary_time=None,
                    config_drift=False, tiny=False) -> bytes:
    if tiny:
        return b"tiny restart\n"
    values = qualification.expected_case_configuration(
        str(intent["case_id"]),
        basename=str(intent["run_basename"]),
        tlim=str(intent["target_time"]),
    )
    values["time/restart_time"] = format(marker_time, ".17g")
    if config_drift:
        values["mhd/cgl_heat_flux"] = "none"
    parameter_dump = (render_parameters(values) + "<par_end>\n").encode()
    payload = bytearray(640)
    struct.pack_into(
        "<d", payload, 232,
        marker_time if binary_time is None else binary_time,
    )
    return parameter_dump + payload


def retain_outputs(packet: dict[str, object], provenance: dict[str, object], *,
                   bad_start=False, bad_cadence=False, nonmonotonic=False,
                   no_work=False, tail_noop=False, nonfinite_snapshot=False,
                   nonpositive_snapshot=None, bad_snapshot_nghost=False,
                   bad_snapshot_geometry=False, missing_snapshot_cadence=False,
                   bad_divb=False, cap_increment_violation=False,
                   near_noop=False, missing_max_ndiv=False,
                   fractional_count=False,
                   bad_snapshot_preheader=False,
                   bad_snapshot_location_size=False,
                   bad_snapshot_variable_size=False,
                   bad_snapshot_active_indices=False,
                   wrong_variables=False, duplicate_logical=False,
                   missing_logical=False, tiny_restart=False,
                   restart_config_drift=False, restart_binary_time=None,
                   hard_bound_violation=False, work_scale=1.0,
                   imbalanced_meshblocks=False,
                   smoke="ok") -> None:
    intent = packet["execution_intent"]
    output = Path(intent["paths"]["output_dir"])
    target = float(intent["target_time"])
    cadence = float(qualification.CASE_POLICIES[intent["case_id"]]["history_dt"])
    times = qualification.expected_history_times(target, cadence)
    if bad_start:
        times[0] = cadence / 4.0
    if bad_cadence:
        times[1] += cadence / 4.0
    if nonmonotonic:
        times[2] = times[1] - cadence / 4.0
    mhd_labels = [
        "time", "mass", "tot-E", "lf_nstage", "lf_dfloor", "lf_pfloor",
        "lf_nonfin", "lf_nonpos", "lf_hardbd", "lf_qface", "lf_qprcap",
        "lf_qpecap", "lf_qprwrk", "lf_qpewrk", "lf_hwproj", "lf_cpwrk",
        "lf_cawrk",
    ]
    user_labels = [
        "time", "mass", "hard_vol", "force_pwr", "force_work", "max_ndiv",
    ]
    mhd_rows = []
    user_rows = []
    for index, time_value in enumerate(times):
        fraction = max(0.0, time_value / target)
        work = 0.0 if no_work else fraction
        work *= work_scale
        if near_noop:
            work *= 1.0e-9
        if tail_noop and index == len(times) - 1:
            work = max(0.0, times[-2] / target)
        mhd_rows.append([
            time_value, 2.0, 10.0 + work, float(index),
            0.0, 0.0, 0.0, 0.0, 0.0,
            float(index * 10), float(index), float(index),
            0.1 * work, 0.2 * work, float(index),
            0.3 * work, 0.4 * work,
        ])
        user_rows.append([
            time_value,
            2.0,
            0.0,
            0.0 if no_work else (1.0e-9 if near_noop else 0.32 * work_scale),
            work,
            1.0e-3 if bad_divb else 1.0e-14,
        ])
    if missing_max_ndiv:
        user_labels.pop()
        user_rows = [row[:-1] for row in user_rows]
    if cap_increment_violation:
        cap_index = mhd_labels.index("lf_qprcap")
        qface_index = mhd_labels.index("lf_qface")
        mhd_rows[-1][cap_index] = (
            mhd_rows[-2][cap_index]
            + mhd_rows[-1][qface_index]
            - mhd_rows[-2][qface_index]
            + 1.0
        )
    if fractional_count:
        mhd_rows[-1][mhd_labels.index("lf_qface")] += 0.5
    basename = intent["run_basename"]
    (output / f"{basename}.mhd.hst").write_text(history_text(mhd_labels, mhd_rows))
    (output / f"{basename}.user.hst").write_text(history_text(user_labels, user_rows))
    binding = qualification.load_json(
        Path(intent["paths"]["job_binding"]), "fixture job binding"
    )
    Path(intent["paths"]["environment_log"]).write_text(
        "\n".join([
            "started_utc=2026-01-01T00:00:00Z",
            f"slurm_job_id={binding['job_id']}",
            (
                "qualification_execution_contract_sha256="
                f"{intent['execution_contract_sha256']}"
            ),
            f"prepared_manifest={intent['paths']['prepared_manifest']}",
            f"nodes={intent['allocation']['nodes']}",
            (
                "ranks="
                f"{int(intent['allocation']['nodes']) * int(intent['allocation']['ranks_per_node'])}"
            ),
            "finished_utc=2026-01-01T00:00:02Z",
            "",
        ])
    )

    ranks = int(intent["allocation"]["nodes"]) * int(
        intent["allocation"]["ranks_per_node"]
    )
    policy = qualification.CASE_POLICIES[intent["case_id"]]
    all_locations = [
        (lx1, lx2, lx3, 0)
        for lx1 in range(policy["mesh_shape"][0] // policy["meshblock_shape"][0])
        for lx2 in range(policy["mesh_shape"][1] // policy["meshblock_shape"][1])
        for lx3 in range(policy["mesh_shape"][2] // policy["meshblock_shape"][2])
    ]
    locations_by_rank = [list(all_locations[rank::ranks]) for rank in range(ranks)]
    if imbalanced_meshblocks:
        locations_by_rank[0].append(locations_by_rank[1].pop())
    variables = list(qualification.REQUIRED_SNAPSHOT_VARIABLES)
    if wrong_variables:
        variables[-1] = "wrong"
    for rank in range(ranks):
        rank_locations = locations_by_rank[rank]
        terminal_locations = list(rank_locations)
        if duplicate_logical and rank == 0:
            terminal_locations[1] = terminal_locations[0]
        if missing_logical and rank == 0:
            terminal_locations.pop()
        bin_dir = output / "bin" / f"rank_{rank:08d}"
        rst_dir = output / "rst" / f"rank_{rank:08d}"
        bin_dir.mkdir(parents=True)
        rst_dir.mkdir(parents=True)
        snapshot_times = qualification.expected_history_times(
            target,
            float(qualification.CASE_POLICIES[intent["case_id"]]["snapshot_dt"]),
        )
        if missing_snapshot_cadence and len(snapshot_times) > 2:
            snapshot_times.pop(1)
        for number, snapshot_time in enumerate(snapshot_times):
            snapshot_locations = (
                terminal_locations
                if abs(snapshot_time - target) <= qualification.ENDPOINT_TOLERANCE
                else rank_locations
            )
            (bin_dir / f"snapshot.{number:05d}.bin").write_bytes(
                snapshot_payload(
                    snapshot_time,
                    snapshot_locations,
                    intent=intent,
                    variables=variables,
                    nonfinite=(
                        nonfinite_snapshot
                        and rank == 0
                        and abs(snapshot_time - target)
                        <= qualification.ENDPOINT_TOLERANCE
                    ),
                    nonpositive=(
                        nonpositive_snapshot
                        if rank == 0
                        and abs(snapshot_time - target)
                        <= qualification.ENDPOINT_TOLERANCE
                        else None
                    ),
                    bad_nghost=bad_snapshot_nghost,
                    bad_geometry=bad_snapshot_geometry,
                    bad_preheader=bad_snapshot_preheader,
                    bad_location_size=bad_snapshot_location_size,
                    bad_variable_size=bad_snapshot_variable_size,
                    bad_active_indices=bad_snapshot_active_indices,
                    hard_bound_violation=(
                        hard_bound_violation
                        and rank == 0
                        and abs(snapshot_time - target)
                        <= qualification.ENDPOINT_TOLERANCE
                    ),
                )
            )
        (rst_dir / "restart.00000.rst").write_bytes(
            restart_payload(
                intent,
                target,
                binary_time=restart_binary_time,
                config_drift=restart_config_drift,
                tiny=tiny_restart,
            )
        )

    if smoke != "missing":
        smoke_dir = Path(intent["paths"]["restart_smoke_dir"])
        smoke_dir.mkdir()
        Path(intent["paths"]["restart_smoke_log"]).write_text("restart load completed\n")
        terminal = output / "rst/rank_00000000/restart.00000.rst"
        digest = qualification.sha256(terminal)
        if smoke == "wrong":
            digest = "0" * 64
        write_json(
            Path(intent["paths"]["restart_smoke_result"]),
            {
                "schema_version": 1,
                "record_type": "cgl_lf_stage_i_restart_load_smoke",
                "execution_contract_sha256": intent["execution_contract_sha256"],
                "executable_sha256": provenance["executable"]["sha256"],
                "job_id": binding["job_id"],
                "restart_path": str(terminal),
                "restart_sha256": digest,
                "restart_size_bytes": terminal.stat().st_size,
                "exit_code": 1 if smoke == "failed" else 0,
                "completed_utc": "2026-01-01T00:00:01Z",
            },
        )


class LiveScheduler:
    def __init__(self, wave, runtimes=None, overlap=False, stale=False):
        self.calls = []
        self.records = {}
        self.account_records = {}
        self.scripts = {}
        self.private_data = "none"
        runtimes = runtimes or {}
        base = datetime.now(timezone.utc).replace(microsecond=0) - timedelta(
            days=3 if stale else 0, hours=12 if not stale else 0
        )
        previous_end = base
        for wave_record in wave["waves"]:
            wave_start = base if overlap and wave_record["wave"] > 1 else previous_end
            wave_end = wave_start
            for packet in wave_record["packets"]:
                intent = packet["execution_intent"]
                binding = qualification.load_json(
                    Path(intent["paths"]["job_binding"]), "fixture binding"
                )
                self.scripts[binding["job_id"]] = Path(
                    intent["paths"]["batch_script"]
                ).read_bytes()
                key = (intent["case_id"], intent["allocation"]["nodes"])
                elapsed = int(runtimes.get(key, 1200 / intent["allocation"]["nodes"]))
                end = wave_start + timedelta(seconds=elapsed)
                submit = wave_start - timedelta(seconds=60)
                nodes = int(intent["allocation"]["nodes"])
                fields = [
                    binding["job_id"], intent["job_name"], "COMPLETED", "0:0",
                    str(nodes), str(elapsed),
                    submit.strftime("%Y-%m-%dT%H:%M:%S%z"),
                    wave_start.strftime("%Y-%m-%dT%H:%M:%S%z"),
                    end.strftime("%Y-%m-%dT%H:%M:%S%z"),
                    qualification.PARTITION, qualification.ACCOUNT.casefold(),
                    f"billing={nodes},node={nodes}",
                    f"billing={nodes},node={nodes}",
                    str(intent["allocation"]["walltime_seconds"] // 60),
                    intent["execution_contract_sha256"],
                    shlex.join(binding["submission_argv"]),
                ]
                self.records[binding["job_id"]] = "|".join(fields) + "\n"
                self.account_records[binding["job_id"]] = "|".join([
                    *fields[:11], "producer",
                ]) + "\n"
                timestamp = lambda value: value.strftime("%Y-%m-%dT%H:%M:%SZ")
                Path(intent["paths"]["environment_log"]).write_text(
                    "\n".join([
                        f"started_utc={timestamp(wave_start)}",
                        f"slurm_job_id={binding['job_id']}",
                        (
                            "qualification_execution_contract_sha256="
                            f"{intent['execution_contract_sha256']}"
                        ),
                        f"prepared_manifest={intent['paths']['prepared_manifest']}",
                        f"nodes={nodes}",
                        (
                            "ranks="
                            f"{nodes * int(intent['allocation']['ranks_per_node'])}"
                        ),
                        f"finished_utc={timestamp(end)}",
                        "",
                    ])
                )
                smoke_result_path = Path(intent["paths"]["restart_smoke_result"])
                smoke_result = json.loads(smoke_result_path.read_text())
                smoke_result["completed_utc"] = timestamp(
                    end - timedelta(seconds=1)
                )
                write_json(smoke_result_path, smoke_result)
                output_binding_path = Path(intent["paths"]["output_binding"])
                output_binding = qualification.output_binding_record(
                    packet,
                    Path(wave["qualification_root"]),
                    wave["provenance"],
                    job_id=binding["job_id"],
                    completed_utc=timestamp(end),
                )
                write_json(output_binding_path, output_binding)
                os.chmod(output_binding_path, 0o440)
                wave_end = max(wave_end, end)
            previous_end = wave_end + timedelta(seconds=60)

    def __call__(self, argv, **kwargs):
        self.calls.append((list(argv), kwargs))
        if argv[0] == str(qualification.SCONTROL):
            if argv[1:] == ["show", "config"]:
                assert {
                    key: value for key, value in kwargs.items() if key != "env"
                } == {
                    "check": True, "capture_output": True, "text": True,
                }
                assert_hardened_scheduler_environment(kwargs["env"])
                return subprocess.CompletedProcess(
                    argv, 0,
                    stdout=f"PrivateData              = {self.private_data}\n",
                    stderr="",
                )
            assert argv[:3] == [
                str(qualification.SCONTROL), "write", "batch_script",
            ]
            assert argv[4] == "-"
            assert {
                key: value for key, value in kwargs.items() if key != "env"
            } == {"check": True, "capture_output": True}
            assert_hardened_scheduler_environment(kwargs["env"])
            return subprocess.CompletedProcess(
                argv, 0, stdout=self.scripts[argv[3]], stderr=b""
            )
        assert argv[0] == str(qualification.SACCT)
        assert_hardened_scheduler_environment(
            kwargs["env"], exact_timestamps=True
        )
        if argv[1] == "-a":
            assert argv[1:8] == [
                "-a", "-A", qualification.ACCOUNT, "-X", "--array", "-n", "-P",
            ]
            assert "-S" in argv and "-E" in argv
            assert argv[-2:] == [
                "-o", ",".join(qualification.ACCOUNT_SACCT_FIELDS),
            ]
            return subprocess.CompletedProcess(
                argv, 0, stdout="".join(self.account_records.values()), stderr=""
            )
        assert argv[1:6] == ["-j", argv[2], "-X", "-n", "-P"]
        job_id = argv[2]
        return subprocess.CompletedProcess(
            argv, 0, stdout=self.records[job_id], stderr=""
        )


def retain_complete_wave(fixture, case_ids=("R04",), *, runtimes=None,
                         overlap=False, output_options=None, max_nodes=10):
    wave = fixture["build_wave"](case_ids, max_nodes=max_nodes)
    wave_path = qualification.materialize_prepared_wave(wave)
    submit_wave(wave, wave_path)
    for packet in packet_map(wave).values():
        intent = packet["execution_intent"]
        key = (intent["case_id"], intent["allocation"]["nodes"])
        retain_outputs(
            packet,
            wave["provenance"],
            **(output_options or {}).get(key, {}),
        )
    scheduler = LiveScheduler(wave, runtimes=runtimes, overlap=overlap)
    evidence_paths = qualification.retain_wave_audit(wave_path, scheduler)
    return wave, wave_path, {
        path.name.split(".", 1)[0]: path for path in evidence_paths
    }, scheduler


def prepared_packet(fixture, case_id="R16", nodes=1):
    wave = fixture["build_wave"]((case_id,))
    wave_path = qualification.materialize_prepared_wave(wave)
    submit_wave(wave, wave_path)
    return wave, packet_map(wave)[(case_id, nodes)]


def build_r17_wave(fixture, monkeypatch):
    monkeypatch.setattr(qualification, "RANKS_PER_NODE", 8)
    return fixture["build_wave"]((qualification.R17_CASE_ID,), max_nodes=8)


def test_exact_frozen_matrix_input_and_physics_are_required(qualification_fixture):
    qualification_fixture["build_wave"](("R04", "R12", "R16"))
    matrix = qualification_fixture["matrix"]
    original = matrix.read_bytes()
    matrix.write_bytes(original + b"\n")
    with pytest.raises(ValueError, match="reviewed frozen matrix"):
        qualification_fixture["build_wave"](("R04",))
    matrix.write_bytes(original)

    policy = qualification.CASE_POLICIES["R04"]
    input_path = qualification_fixture["source"] / policy["input_relative_path"]
    input_path.write_text(input_path.read_text().replace(
        "cgl_heat_flux = landau_fluid", "cgl_heat_flux = none"
    ))
    policy["input_sha256"] = hashlib.sha256(input_path.read_bytes()).hexdigest()
    with pytest.raises(ValueError, match="physics contract"):
        qualification.validate_case_input_contract("R04", input_path)


@pytest.mark.parametrize("case_id", ["R04", "R12", "R16", "R17"])
def test_each_reviewed_qualification_case_contract_is_exact(
    qualification_fixture, case_id,
):
    policy = qualification.CASE_POLICIES[case_id]
    input_path = qualification_fixture["source"] / policy["input_relative_path"]
    input_path.write_text(input_path.read_text().replace(
        "beta0 = 10.0", "beta0 = 20.0"
    ))
    policy["input_sha256"] = hashlib.sha256(input_path.read_bytes()).hexdigest()
    with pytest.raises(ValueError, match="physics contract"):
        qualification.validate_case_input_contract(case_id, input_path)


def test_r12_measured_fresh_profile_policy_is_exact_parentless_and_keeps_cadence():
    profile = qualification.require_r12_fresh_profile_policy()
    assert profile == {
        "case_id": "R12",
        "segment": "s01_rankio_t0_t0p12",
        "start_time": 0.0,
        "target_time": 0.12,
        "nodes": 4,
        "ranks_per_node": 8,
        "total_ranks": 32,
        "walltime": "02:00:00",
        "athena_walltime": "01:50:00",
        "parent_job_id": None,
        "parent_result": None,
        "parent_segment": None,
        "restart_file": None,
        "restart_file_sha256": None,
        "restart_time": None,
    }
    policy = qualification.CASE_POLICIES["R12"]
    assert policy["node_profiles"] == (4,)
    assert policy["target_time"] == 0.12
    assert policy["history_dt"] == 0.02
    assert policy["snapshot_dt"] == 0.25
    assert qualification.COMMON_INPUT_CONTRACT["output1/dt"] == "0.02"
    assert qualification.COMMON_INPUT_CONTRACT["output2/dt"] == "0.25"
    assert qualification.COMMON_INPUT_CONTRACT["output3/dt"] == "1.0"
    assert "s00_rankio_t0_t0p25" not in json.dumps(profile)
    assert "4766856" not in json.dumps(profile)


@pytest.mark.parametrize(
    ("section", "key", "value"),
    [
        ("profile", "segment", "s01_rankio_t0_t0p25"),
        ("profile", "segment", "s00_rankio_t0_t0p25"),
        ("profile", "target_time", 0.25),
        ("profile", "nodes", 2),
        ("profile", "total_ranks", 16),
        ("profile", "walltime", "01:00:00"),
        ("profile", "athena_walltime", "00:50:00"),
        ("profile", "parent_job_id", "4766856"),
        ("profile", "parent_segment", "s00_rankio_t0_t0p25"),
        ("profile", "restart_file", "/historical/restart.rst"),
        ("policy", "history_dt", 0.01),
        ("policy", "snapshot_dt", 0.12),
    ],
)
def test_r12_measured_fresh_profile_policy_rejects_drift(
    qualification_fixture, monkeypatch, section, key, value,
):
    policy = qualification.CASE_POLICIES["R12"]
    target = policy["fixed_fresh_profile"] if section == "profile" else policy
    monkeypatch.setitem(target, key, value)
    with pytest.raises(ValueError, match="R12 qualification policy differs"):
        qualification.require_r12_fresh_profile_policy()


def test_r12_preparation_is_one_exact_fresh_four_node_packet(
    qualification_fixture, monkeypatch,
):
    profile = deepcopy(qualification.R12_FRESH_RERUN_PROFILE)
    profile["ranks_per_node"] = 8
    profile["total_ranks"] = 32
    policy = qualification.CASE_POLICIES["R12"]
    policy["fixed_fresh_profile"] = profile
    monkeypatch.setattr(qualification, "R12_FRESH_RERUN_PROFILE", profile)
    monkeypatch.setattr(qualification, "RANKS_PER_NODE", 8)

    wave = qualification_fixture["build_wave"](("R12",))
    packets = packet_map(wave)
    assert set(packets) == {("R12", 4)}
    intent = packets[("R12", 4)]["execution_intent"]
    assert intent["packet_id"] == "R12_n04_t0p12"
    assert intent["target_time"] == 0.12
    assert intent["allocation"] == {
        "nodes": 4,
        "walltime": "02:00:00",
        "walltime_seconds": 7200,
        "athena_walltime": "01:50:00",
        "athena_walltime_seconds": 6600,
        "ranks_per_node": 8,
        "cpus_per_task": 1,
    }
    script = intent["authenticated_commands"]["batch_script_text"]
    main_launch, restart_smoke = script.split("TERMINAL_RESTART=", 1)
    assert "time/tlim=0.12" in main_launch
    assert " -r " not in main_launch
    assert " -r " in restart_smoke


def test_r17_preparation_is_exactly_one_exclusive_eight_node_packet(
    qualification_fixture, monkeypatch,
):
    wave = build_r17_wave(qualification_fixture, monkeypatch)
    packets = packet_map(wave)
    assert set(packets) == {("R17", 8)}
    assert wave["policy"]["operational_only"] is True
    assert wave["policy"]["max_concurrent_nodes"] == 8
    assert len(wave["waves"]) == 1
    assert wave["waves"][0]["total_nodes"] == 8
    intent = packets[("R17", 8)]["execution_intent"]
    assert intent["allocation"]["nodes"] == 8
    assert intent["allocation"]["ranks_per_node"] == 8

    with pytest.raises(ValueError, match="one exclusive eight-node wave"):
        qualification_fixture["build_wave"](("R17",), max_nodes=10)
    with pytest.raises(ValueError, match="one exclusive eight-node wave"):
        qualification_fixture["build_wave"](("R16", "R17"), max_nodes=8)


def test_full_build_manifest_inventory_is_authenticated(
    qualification_fixture,
):
    inventory = qualification_fixture["provenance"]["build_manifest"]["inventory"]
    assert [record["name"] for record in inventory] == [
        "athena.sha256", "environment.txt",
    ]
    assert qualification.retained_inventory_sha256(inventory) == (
        qualification_fixture["provenance"]["build_manifest"]["inventory_sha256"]
    )
    (qualification_fixture["build_manifest"] / "unreviewed.txt").write_text("drift\n")
    selected = [
        qualification_fixture["source"]
        / qualification.CASE_POLICIES["R04"]["input_relative_path"]
    ]
    with pytest.raises(ValueError, match="inventory differs"):
        qualification.validate_prepared_provenance(
            qualification_fixture["provenance"],
            qualification_fixture["root"],
            selected,
        )


def test_provenance_hardlinks_and_mismatched_build_revision_are_rejected(
    qualification_fixture,
):
    policy = qualification.CASE_POLICIES["R04"]
    input_path = qualification_fixture["source"] / policy["input_relative_path"]
    external = input_path.with_suffix(".external")
    input_path.rename(external)
    input_path.hardlink_to(external)
    with pytest.raises(ValueError, match="changed while it was opened"):
        qualification.validate_case_input_contract("R04", input_path)

    provenance = deepcopy(qualification_fixture["provenance"])
    provenance["executable"]["revision"] = "0" * 40
    selected = [
        qualification_fixture["source"]
        / qualification.CASE_POLICIES[case_id]["input_relative_path"]
        for case_id in ("R04", "R12", "R16")
    ]
    with pytest.raises(ValueError, match="differs from frozen source revision"):
        qualification.validate_prepared_provenance(
            provenance, qualification_fixture["root"], selected
        )


def test_source_bundle_rebinding_race_is_rejected(
    qualification_fixture, monkeypatch,
):
    bundle = qualification_fixture["bundle"]
    original = qualification.open_regular_fd
    raced = False

    @contextmanager
    def racing_open(path, label, **kwargs):
        nonlocal raced
        with original(path, label, **kwargs) as descriptor:
            if label == "source bundle" and not raced:
                raced = True
                bundle.rename(bundle.with_suffix(".moved"))
                bundle.write_bytes(b"fabricated bundle\n")
            yield descriptor

    monkeypatch.setattr(qualification, "open_regular_fd", racing_open)
    with pytest.raises(ValueError, match="changed during inspection"):
        qualification.source_bundle_provenance(
            bundle,
            qualification_fixture["root"],
            [qualification.FROZEN_SOURCE_REVISION],
        )


def test_scheduler_and_reexec_environments_strip_caller_execution_controls(
    tmp_path, monkeypatch,
):
    poisoned = {
        "HOME": str(tmp_path / "hostile-home"),
        "PATH": str(tmp_path / "hostile-bin"),
        "SBATCH_ACCOUNT": "attacker",
        "SBATCH_EXPORT": "ALL,POISON=1",
        "SLURM_CONF": str(tmp_path / "hostile-slurm.conf"),
        "SLURM_TIME_FORMAT": "attacker-format",
        "LD_PRELOAD": str(tmp_path / "hostile.so"),
        "PYTHONPATH": str(tmp_path / "hostile-python"),
        "PERL5OPT": "-Mhostile",
        "RUBYOPT": "-rhostile",
        "BASH_ENV": str(tmp_path / "hostile-bash-env"),
        "GIT_CONFIG_GLOBAL": str(tmp_path / "hostile.gitconfig"),
        qualification.SELF_DESCRIPTOR_ENV: "99",
        qualification.SELF_SOURCE_ENV: "/hostile/source",
        qualification.REPOSITORY_ROOT_ENV: "/hostile/repository",
        "_CGL_LF_OTHER_PRIVATE_REEXEC": "hostile",
        "KEEP_ME": "retained",
    }
    for key, value in poisoned.items():
        monkeypatch.setenv(key, value)

    environment = qualification.scheduler_environment()
    assert environment["PATH"] == qualification.TRUSTED_SYSTEM_PATH
    assert environment["LC_ALL"] == "C"
    assert environment["KEEP_ME"] == "retained"
    assert not any(
        key in environment
        for key in poisoned
        if key not in {"HOME", "PATH", "KEEP_ME"}
    )
    exact = qualification.scheduler_environment(exact_timestamps=True)
    assert exact["SLURM_TIME_FORMAT"] == qualification.SLURM_TIME_FORMAT
    assert not any(
        key.startswith(("SBATCH_", "SLURM_"))
        for key in exact
        if key != "SLURM_TIME_FORMAT"
    )

    reexec = qualification.reexec_environment(7, UTILITY, REPOSITORY)
    assert reexec[qualification.SELF_DESCRIPTOR_ENV] == "7"
    assert reexec[qualification.SELF_SOURCE_ENV] == str(UTILITY)
    assert reexec[qualification.REPOSITORY_ROOT_ENV] == str(REPOSITORY)
    assert reexec["HOME"] == "/nonexistent"
    assert reexec["PATH"] == qualification.TRUSTED_SYSTEM_PATH
    assert "LD_PRELOAD" not in reexec
    assert "PYTHONPATH" not in reexec
    assert "SBATCH_ACCOUNT" not in reexec
    assert "SLURM_CONF" not in reexec


def test_reexec_uses_authenticated_system_python_descriptor_and_isolated_startup(
    monkeypatch,
):
    module = load_utility()
    calls = []

    class ReexecIntercept(Exception):
        pass

    def fake_execve(executable, argv, environment):
        python_descriptor = int(executable.rsplit("/", 1)[1])
        profile = os.fstat(python_descriptor)
        assert profile.st_uid == 0
        assert profile.st_nlink == 1
        assert profile.st_mode & 0o111
        assert os.get_inheritable(python_descriptor)
        calls.append((executable, argv, environment))
        raise ReexecIntercept

    monkeypatch.setenv("HOME", "/attacker/home")
    monkeypatch.setenv("PYTHONPATH", "/attacker/python")
    monkeypatch.setattr(module.os, "execve", fake_execve)
    descriptor = os.open(UTILITY, os.O_RDONLY)
    try:
        with pytest.raises(ReexecIntercept):
            module.reexec_authenticated_self(
                descriptor, UTILITY.resolve(), REPOSITORY.resolve(), ["--help"]
            )
    finally:
        os.close(descriptor)

    assert len(calls) == 1
    executable, argv, environment = calls[0]
    assert executable.startswith("/proc/self/fd/")
    assert argv[:4] == [str(module.SYSTEM_PYTHON), "-I", "-S", "-B"]
    assert argv[4].startswith("/proc/self/fd/")
    assert argv[5:] == ["--help"]
    assert environment["HOME"] == "/nonexistent"
    assert environment["PATH"] == module.TRUSTED_SYSTEM_PATH
    assert "PYTHONPATH" not in environment


def test_authenticated_system_python_rejects_user_controlled_interpreter(
    tmp_path, monkeypatch,
):
    module = load_utility()
    hostile = tmp_path / "python3.11"
    hostile.write_text("#!/bin/sh\nexit 0\n")
    hostile.chmod(0o755)
    monkeypatch.setattr(module, "SYSTEM_PYTHON", hostile)
    descriptor = os.open(UTILITY, os.O_RDONLY)
    try:
        with pytest.raises(ValueError, match="trusted system directory profile"):
            module.reexec_authenticated_self(
                descriptor, UTILITY.resolve(), REPOSITORY.resolve(), []
            )
    finally:
        os.close(descriptor)


def test_exclusive_write_rejects_parent_rebinding_before_creation(
    tmp_path, monkeypatch,
):
    parent = tmp_path / "retained"
    parent.mkdir()
    destination = parent / "artifact.json"
    displaced = tmp_path / "retained.displaced"
    original = qualification.open_directory_fd
    raced = False

    @contextmanager
    def racing_open(path, label):
        nonlocal raced
        with original(path, label) as retained:
            if Path(path) == parent and not raced:
                raced = True
                parent.rename(displaced)
                parent.mkdir()
            yield retained

    monkeypatch.setattr(qualification, "open_directory_fd", racing_open)
    with pytest.raises(ValueError, match="path changed before mutation"):
        qualification.write_exclusive(destination, b"authenticated\n")
    assert not destination.exists()
    assert not (displaced / destination.name).exists()


def test_git_uses_authenticated_root_owned_descriptor_and_isolated_environment(
    tmp_path, monkeypatch,
):
    for key, value in {
        "PATH": str(tmp_path / "hostile-bin"),
        "LD_PRELOAD": str(tmp_path / "hostile.so"),
        "PYTHONPATH": str(tmp_path / "hostile-python"),
        "GIT_CONFIG_GLOBAL": str(tmp_path / "hostile-global"),
        "GIT_CONFIG_SYSTEM": str(tmp_path / "hostile-system"),
        "GIT_CONFIG_PARAMETERS": "'core.repositoryformatversion=99'",
        "GIT_EXEC_PATH": str(tmp_path / "hostile-git-exec"),
    }.items():
        monkeypatch.setenv(key, value)
    calls = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        return subprocess.CompletedProcess(command, 0, stdout="git version fixture\n", stderr="")

    monkeypatch.setattr(qualification.subprocess, "run", fake_run)
    completed = qualification.git_run(["--version"], capture_output=True, text=True)
    assert completed.returncode == 0
    assert len(calls) == 1
    command, kwargs = calls[0]
    assert command[:2] == [str(qualification.GIT), "--no-replace-objects"]
    assert kwargs["executable"].startswith("/proc/self/fd/")
    assert kwargs["pass_fds"]
    assert kwargs["stdin"] == subprocess.DEVNULL
    assert kwargs["env"] == qualification.hardened_git_environment()
    assert kwargs["env"]["PATH"] == qualification.TRUSTED_SYSTEM_PATH
    assert kwargs["env"]["GIT_CONFIG_GLOBAL"] == os.devnull
    assert kwargs["env"]["GIT_CONFIG_SYSTEM"] == os.devnull
    assert kwargs["env"]["GIT_CONFIG_NOSYSTEM"] == "1"
    assert "LD_PRELOAD" not in kwargs["env"]
    assert "GIT_CONFIG_PARAMETERS" not in kwargs["env"]

    fake_git = tmp_path / "git"
    fake_git.write_text("#!/bin/sh\nexit 0\n")
    fake_git.chmod(0o755)
    monkeypatch.setattr(qualification, "GIT", fake_git)
    with pytest.raises(ValueError, match="trusted system (directory|executable) profile"):
        qualification.git_run(["--version"])


def test_committed_source_revision_ignores_hostile_path_and_git_configuration(
    qualification_fixture, tmp_path, monkeypatch,
):
    hostile_bin = tmp_path / "bin"
    hostile_bin.mkdir()
    marker = tmp_path / "caller-git-ran"
    hostile_git = hostile_bin / "git"
    hostile_git.write_text(f"#!/bin/sh\n: > {marker}\nexit 99\n")
    hostile_git.chmod(0o755)
    hostile_config = tmp_path / "hostile.gitconfig"
    hostile_config.write_text("[core]\n\trepositoryformatversion = 99\n")
    for key, value in {
        "PATH": str(hostile_bin),
        "GIT_CONFIG_GLOBAL": str(hostile_config),
        "GIT_CONFIG_SYSTEM": str(hostile_config),
        "GIT_CONFIG_NOSYSTEM": "0",
        "GIT_CONFIG_PARAMETERS": "'core.repositoryformatversion=99'",
        "GIT_EXEC_PATH": str(hostile_bin),
    }.items():
        monkeypatch.setenv(key, value)
    revision = qualification.committed_source_revision(
        qualification_fixture["source"], [qualification_fixture["matrix"]]
    )
    assert revision == qualification.FROZEN_SOURCE_REVISION
    assert not marker.exists()


def test_qualification_rejects_orphaned_private_reexec_metadata(monkeypatch):
    module = load_utility()
    monkeypatch.setenv(module.SELF_SOURCE_ENV, "/attacker/source.py")
    monkeypatch.setenv(module.REPOSITORY_ROOT_ENV, "/attacker/repository")
    with pytest.raises(
        ValueError,
        match="private qualification reexecution path is forbidden",
    ):
        module.authenticate_self([])
    monkeypatch.setenv(module.SELF_DESCRIPTOR_ENV, "0")
    with pytest.raises(ValueError, match="descriptor marker is not attached"):
        module.authenticate_self([])


def test_authenticated_qualification_descriptor_controls_private_source_metadata(
    qualification_fixture, monkeypatch,
):
    module = load_utility()
    monkeypatch.setattr(module, "require_authenticated_reexec_runtime", lambda: None)
    source = qualification_fixture["helper"]
    repository = qualification_fixture["source"]
    descriptor = os.open(source, os.O_RDONLY)
    try:
        monkeypatch.setattr(module, "__file__", f"/proc/self/fd/{descriptor}")
        monkeypatch.setenv(module.SELF_DESCRIPTOR_ENV, str(descriptor))
        monkeypatch.setenv(module.SELF_SOURCE_ENV, str(source))
        monkeypatch.setenv(module.REPOSITORY_ROOT_ENV, str(repository))
        authenticated = module.authenticate_self([])
        assert authenticated[:2] == (source, repository)
        assert authenticated[2] == qualification.sha256(source)

        monkeypatch.setenv(module.SELF_SOURCE_ENV, str(qualification_fixture["matrix"]))
        with pytest.raises(ValueError, match="source/repository relationship differs"):
            module.authenticate_self([])
    finally:
        os.close(descriptor)


def test_forged_authenticated_descriptor_rejects_nonisolated_runtime(
    qualification_fixture, monkeypatch,
):
    module = load_utility()
    source = qualification_fixture["helper"]
    repository = qualification_fixture["source"]
    descriptor = os.open(source, os.O_RDONLY)
    try:
        monkeypatch.setattr(module, "__file__", f"/proc/self/fd/{descriptor}")
        monkeypatch.setenv(module.SELF_DESCRIPTOR_ENV, str(descriptor))
        monkeypatch.setenv(module.SELF_SOURCE_ENV, str(source))
        monkeypatch.setenv(module.REPOSITORY_ROOT_ENV, str(repository))
        monkeypatch.setenv("HOME", "/attacker/home")
        with pytest.raises(ValueError, match="HOME is not sanitized"):
            module.authenticate_self([])
    finally:
        os.close(descriptor)


def test_reexec_runtime_requires_isolated_python_flags(monkeypatch):
    module = load_utility()
    monkeypatch.setenv("HOME", "/nonexistent")
    with pytest.raises(ValueError, match="reexecution is not isolated"):
        module.require_authenticated_reexec_runtime()


def test_materialized_script_is_exact_valid_and_cli_rejects_scheduler_sources(
    qualification_fixture,
):
    wave = qualification_fixture["build_wave"](("R16",))
    qualification.materialize_prepared_wave(wave)
    packet = packet_map(wave)[("R16", 1)]
    script = Path(packet["execution_intent"]["paths"]["batch_script"])
    subprocess.run(["/bin/bash", "-n", str(script)], check=True)
    assert "--wrap" not in packet["execution_intent"]["authenticated_commands"][
        "base_submission_argv"
    ]
    with pytest.raises(SystemExit):
        qualification.parser().parse_args([
            "audit-wave", "--prepared-wave", "/tmp/wave.json",
            "--scheduler-raw", "1=/tmp/fabricated",
        ])


def test_submission_uses_exact_stdin_dependency_journal_and_response_binding(
    qualification_fixture,
):
    wave = qualification_fixture["build_wave"](("R04", "R12"))
    wave_path = qualification.materialize_prepared_wave(wave)
    result, submitter, queue = submit_wave(wave, wave_path)
    assert len(result["job_bindings"]) == 4
    assert len(queue.calls) == len(wave["waves"]) + len(submitter.calls)
    packets = [
        packet for wave_record in wave["waves"] for packet in wave_record["packets"]
    ]
    for job_id, packet, (argv, payload) in zip(
        range(1000, 1000 + len(packets)), packets, submitter.calls
    ):
        commands = packet["execution_intent"]["authenticated_commands"]
        assert payload == commands["batch_script_text"].encode()
        assert hashlib.sha256(payload).hexdigest() == commands["batch_script_sha256"]
        assert argv[:3] == commands["base_submission_argv"]
        assert submitter.scripts[str(job_id)] == payload
        scheduler_script = Path(
            packet["execution_intent"]["paths"]["scheduler_batch_script"]
        )
        assert scheduler_script.read_bytes() == payload
        assert scheduler_script.stat().st_mode & 0o222 == 0
    journal = Path(wave["qualification_root"]) / "submission-journal"
    assert journal.stat().st_mode & 0o222 == 0
    assert all(path.stat().st_mode & 0o222 == 0 for path in journal.iterdir())
    first_binding = qualification.load_json(
        Path(packets[0]["execution_intent"]["paths"]["job_binding"]), "binding"
    )
    assert first_binding["sbatch_stdout"] == "1000;frontier\n"
    assert first_binding["submission_stdin_size_bytes"] == len(submitter.calls[0][1])


def test_submission_path_swap_cannot_change_authenticated_stdin(
    qualification_fixture,
):
    wave = qualification_fixture["build_wave"](("R16",))
    wave_path = qualification.materialize_prepared_wave(wave)
    packet = packet_map(wave)[("R16", 1)]
    script_path = Path(packet["execution_intent"]["paths"]["batch_script"])
    expected = script_path.read_bytes()

    def swap_first(call, _argv, payload):
        if call == 1:
            assert payload == expected
            moved = script_path.with_suffix(".moved")
            script_path.rename(moved)
            script_path.write_bytes(b"#!/bin/bash\nexit 99\n")

    submitter = SbatchSubmitter(hook=swap_first)
    submit_wave(wave, wave_path, submitter=submitter)
    assert submitter.calls[0][1] == expected


def test_submission_rejects_scheduler_stored_script_byte_mismatch(
    qualification_fixture,
):
    wave = qualification_fixture["build_wave"](("R16",))
    wave_path = qualification.materialize_prepared_wave(wave)
    queue = EmptyQueue()
    submitter = SbatchSubmitter(queue=queue)

    def mismatched_batch_script(argv, **kwargs):
        completed = submitter.read_batch_script(argv, **kwargs)
        completed.stdout += b"\n# scheduler mutation\n"
        return completed

    canceler = CancelRunner(queue)
    with pytest.raises(ValueError, match="Slurm-stored batch script"):
        qualification.submit_all_waves(
            wave_path,
            runner=submitter,
            queue_runner=queue,
            batch_script_runner=mismatched_batch_script,
            cancel_runner=canceler,
        )
    packet = packet_map(wave)[("R16", 1)]
    paths = packet["execution_intent"]["paths"]
    recovery = qualification.load_json(
        Path(paths["submission_recovery"]), "submission recovery"
    )
    assert canceler.calls == [[str(qualification.SCANCEL), "1000"]]
    assert queue.submitted == {}
    assert recovery["state"] == "cancelled_confirmed"
    assert recovery["cancellation"]["confirmed_absent"] is True
    assert not Path(paths["job_binding"]).exists()


def test_ambiguous_submit_is_durably_recovered_before_any_retry(
    qualification_fixture,
):
    wave = qualification_fixture["build_wave"](("R16",))
    wave_path = qualification.materialize_prepared_wave(wave)
    queue = EmptyQueue()
    submitter = AmbiguousSubmitter(queue=queue)
    canceler = CancelRunner(queue)
    packet = packet_map(wave)[("R16", 1)]
    paths = packet["execution_intent"]["paths"]

    with pytest.raises(ValueError, match="outcome is ambiguous"):
        qualification.submit_all_waves(
            wave_path,
            runner=submitter,
            queue_runner=queue,
            batch_script_runner=submitter.read_batch_script,
            cancel_runner=canceler,
        )
    receipt = qualification.load_json(
        Path(paths["submission_receipt"]), "submission receipt"
    )
    assert receipt["state"] == "ambiguous_sbatch_outcome"
    assert receipt["job_id"] is None
    assert "1000" in queue.submitted
    assert not canceler.calls

    result = qualification.recover_ambiguous_submission(
        wave_path,
        packet["execution_intent"]["packet_id"],
        "1000",
        "operator identified the candidate from the account queue",
        batch_script_runner=submitter.read_batch_script,
        cancel_runner=canceler,
        queue_runner=queue,
    )
    assert result["cancellation"]["confirmed_absent"] is True
    assert queue.submitted == {}
    recovery = qualification.load_json(
        Path(paths["submission_recovery"]), "submission recovery"
    )
    assert recovery["state"] == "cancelled_confirmed"
    assert recovery["job_id"] == "1000"
    with pytest.raises(ValueError, match="journal already exists"):
        qualification.submit_all_waves(
            wave_path,
            runner=submitter,
            queue_runner=queue,
            batch_script_runner=submitter.read_batch_script,
            cancel_runner=canceler,
        )


def test_sbatch_error_with_exact_job_id_is_automatically_cancelled(
    qualification_fixture,
):
    wave = qualification_fixture["build_wave"](("R16",))
    wave_path = qualification.materialize_prepared_wave(wave)
    queue = EmptyQueue()
    submitter = ErrorAfterIdSubmitter(queue=queue)
    canceler = CancelRunner(queue)
    packet = packet_map(wave)[("R16", 1)]
    paths = packet["execution_intent"]["paths"]

    with pytest.raises(ValueError, match="failed after exposing an exact job ID"):
        qualification.submit_all_waves(
            wave_path,
            runner=submitter,
            queue_runner=queue,
            batch_script_runner=submitter.read_batch_script,
            cancel_runner=canceler,
        )
    receipt = qualification.load_json(
        Path(paths["submission_receipt"]), "submission receipt"
    )
    recovery = qualification.load_json(
        Path(paths["submission_recovery"]), "submission recovery"
    )
    assert receipt["state"] == "scheduler_job_id_recovered_from_error"
    assert receipt["job_id"] == "1000"
    assert recovery["state"] == "cancelled_confirmed"
    assert queue.submitted == {}


def test_known_job_receipt_failure_cancels_and_retains_boundary_recovery(
    qualification_fixture, monkeypatch,
):
    wave = qualification_fixture["build_wave"](("R16",))
    wave_path = qualification.materialize_prepared_wave(wave)
    queue = EmptyQueue()
    submitter = SbatchSubmitter(queue=queue)
    canceler = CancelRunner(queue)
    packet = packet_map(wave)[("R16", 1)]
    paths = packet["execution_intent"]["paths"]

    def fail_receipt(*_args, **_kwargs):
        raise OSError("receipt storage unavailable")

    monkeypatch.setattr(qualification, "write_submission_receipt", fail_receipt)
    with pytest.raises(ValueError, match="receipt storage unavailable"):
        qualification.submit_all_waves(
            wave_path,
            runner=submitter,
            queue_runner=queue,
            batch_script_runner=submitter.read_batch_script,
            cancel_runner=canceler,
        )
    recovery = qualification.load_json(
        Path(paths["submission_recovery"]), "submission recovery"
    )
    assert queue.submitted == {}
    assert recovery["state"] == "cancelled_confirmed"
    assert recovery["submission_receipt"] is None
    assert recovery["submission_boundary"]["path"].endswith("-boundary.json")


@pytest.mark.parametrize(
    ("queue_text", "message"),
    [
        ("900|900|N/A|production|RUNNING|4|producer\n", "ten-node ceiling"),
        (
            "901|901|N/A|cglq_other|PENDING|1|producer\n",
            "another qualification root",
        ),
    ],
)
def test_live_project_queue_enforces_global_ceiling_and_root_exclusion(
    qualification_fixture, queue_text, message,
):
    wave = qualification_fixture["build_wave"](("R04",))
    wave_path = qualification.materialize_prepared_wave(wave)
    with pytest.raises(ValueError, match=message):
        submit_wave(wave, wave_path, queue=EmptyQueue(queue_text))


def test_r17_submission_requires_empty_account_queue(
    qualification_fixture, monkeypatch,
):
    wave = build_r17_wave(qualification_fixture, monkeypatch)
    wave_path = qualification.materialize_prepared_wave(wave)
    queue = EmptyQueue("900|900|N/A|production|RUNNING|1|producer\n")
    with pytest.raises(ValueError, match="requires an empty account queue"):
        submit_wave(wave, wave_path, queue=queue)


@pytest.mark.parametrize("state_kind", ["reservation", "transaction"])
def test_r17_submission_requires_empty_canonical_stage_i_state(
    qualification_fixture, monkeypatch, state_kind,
):
    accounting = qualification_fixture["root"] / "accounting"
    if state_kind == "reservation":
        write_json(
            accounting / qualification.STAGE_I_RESERVATIONS_NAME,
            [{"state": "prepared"}],
        )
        message = "active reservations"
    else:
        (
            accounting / qualification.STAGE_I_TRANSACTION_NAMES[0] / "pending.json"
        ).write_text("{}\n")
        message = "transaction store is not empty"
    wave = build_r17_wave(qualification_fixture, monkeypatch)
    wave_path = qualification.materialize_prepared_wave(wave)
    with pytest.raises(ValueError, match=message):
        submit_wave(wave, wave_path)


def test_r17_submission_takes_existing_stage_i_lock(
    qualification_fixture, monkeypatch,
):
    wave = build_r17_wave(qualification_fixture, monkeypatch)
    wave_path = qualification.materialize_prepared_wave(wave)
    with qualification.stage_i_exclusivity_lock(qualification_fixture["root"]):
        with pytest.raises(ValueError, match="another Stage I mutation"):
            submit_wave(wave, wave_path)


def test_r17_lock_rebinding_after_admission_fails_before_scheduler_mutation(
    qualification_fixture, monkeypatch,
):
    wave = build_r17_wave(qualification_fixture, monkeypatch)
    wave_path = qualification.materialize_prepared_wave(wave)
    root = qualification_fixture["root"]
    lock_path = root / qualification.STAGE_I_LOCK_NAME
    displaced = root / f"{qualification.STAGE_I_LOCK_NAME}.displaced"

    class RebindingQueue(EmptyQueue):
        def __init__(self):
            super().__init__()
            self.raced = False

        def __call__(self, argv, **kwargs):
            completed = super().__call__(argv, **kwargs)
            if not self.raced:
                self.raced = True
                lock_path.rename(displaced)
                lock_path.write_bytes(b"")
                lock_path.chmod(0o644)
            return completed

    queue = RebindingQueue()
    submitter = SbatchSubmitter(queue=queue)
    with pytest.raises(ValueError, match="mutation lock.*changed"):
        submit_wave(wave, wave_path, queue=queue, submitter=submitter)
    assert submitter.calls == []


def test_r17_post_submission_foreign_job_race_is_cancelled(
    qualification_fixture, monkeypatch,
):
    wave = build_r17_wave(qualification_fixture, monkeypatch)
    wave_path = qualification.materialize_prepared_wave(wave)
    queue = PostSubmissionRaceQueue()
    canceler = CancelRunner(queue)
    with pytest.raises(ValueError, match="lost exclusive account ownership"):
        submit_wave(wave, wave_path, queue=queue, canceler=canceler)
    assert canceler.calls == [[str(qualification.SCANCEL), "1000"]]


def test_r17_retained_admission_fails_closed_on_fabricated_empty_state(
    qualification_fixture, monkeypatch,
):
    wave = build_r17_wave(qualification_fixture, monkeypatch)
    wave_path = qualification.materialize_prepared_wave(wave)
    submit_wave(wave, wave_path)
    packet = packet_map(wave)[("R17", 8)]
    binding = qualification.load_json(
        Path(packet["execution_intent"]["paths"]["job_binding"]), "R17 binding"
    )
    boundary = qualification.load_json(
        Path(binding["submission_boundary"]["path"]), "R17 boundary"
    )
    admission = deepcopy(boundary["queue_admission"])
    admission["canonical_stage_i_state"]["active_reservations"] = 1
    with pytest.raises(ValueError, match="admission was not empty"):
        qualification.validate_r17_submission_admission(
            admission, qualification_fixture["root"]
        )


def test_live_queue_accepts_frontier_array_job_identifiers():
    queue = EmptyQueue(
        "".join(
            f"123_{task}|123|{task}|production|PENDING|2|producer\n"
            for task in range(1, 5)
        )
    )
    jobs = qualification.live_queue_jobs(queue)
    assert len(jobs) == 4
    assert sum(job["nodes"] for job in jobs) == 8


def test_live_queue_rejects_condensed_array_identity():
    queue = EmptyQueue(
        "123_[1-4]|123|N/A|production|PENDING|2|producer\n"
    )
    with pytest.raises(ValueError, match="invalid job identity"):
        qualification.live_queue_jobs(queue)


def test_expanded_arrays_enforce_account_ceiling(
    qualification_fixture,
):
    wave = qualification_fixture["build_wave"](("R16",))
    wave_path = qualification.materialize_prepared_wave(wave)
    array_queue = EmptyQueue(
        "".join(
            f"123_{task}|123|{task}|production|PENDING|2|producer\n"
            for task in range(1, 5)
        )
    )
    with pytest.raises(ValueError, match="ten-node ceiling"):
        submit_wave(wave, wave_path, queue=array_queue)


def test_post_submission_queue_race_is_reaudited(qualification_fixture):
    wave = qualification_fixture["build_wave"](("R16",))
    wave_path = qualification.materialize_prepared_wave(wave)
    queue = PostSubmissionRaceQueue()
    canceler = CancelRunner(queue)
    with pytest.raises(ValueError, match="post-submission account queue"):
        submit_wave(wave, wave_path, queue=queue, canceler=canceler)
    packet = packet_map(wave)[("R16", 1)]
    recovery = qualification.load_json(
        Path(packet["execution_intent"]["paths"]["submission_recovery"]),
        "submission recovery",
    )
    assert canceler.calls == [[str(qualification.SCANCEL), "1000"]]
    assert recovery["cancellation"]["confirmed_absent"] is True


def test_global_qualification_lock_excludes_parallel_roots(qualification_fixture):
    wave = qualification_fixture["build_wave"](("R16",))
    wave_path = qualification.materialize_prepared_wave(wave)
    with qualification.global_qualification_lock(qualification_fixture["root"]):
        with pytest.raises(ValueError, match="holds the lock"):
            submit_wave(wave, wave_path)


def test_journal_mutation_is_rejected(qualification_fixture):
    wave, _, evidence, scheduler = retain_complete_wave(
        qualification_fixture, ("R16",)
    )
    record = sorted(
        (Path(wave["qualification_root"]) / "submission-journal").iterdir()
    )[0]
    os.chmod(record, 0o640)
    record.write_bytes(record.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="journal records must be immutable"):
        qualification.select_profile(evidence["R16"], scheduler)


def test_live_scheduler_timestamp_command_and_fabrication_rejection(
    qualification_fixture,
):
    wave, _, evidence, scheduler = retain_complete_wave(
        qualification_fixture, ("R16",)
    )
    assert scheduler.calls
    assert all(
        call[1]["env"]["SLURM_TIME_FORMAT"] == qualification.SLURM_TIME_FORMAT
        for call in scheduler.calls
        if call[0][0] == str(qualification.SACCT)
    )
    packet = packet_map(wave)[("R16", 1)]
    binding = qualification.load_json(
        Path(packet["execution_intent"]["paths"]["job_binding"]), "binding"
    )
    scheduler.records[binding["job_id"]] = scheduler.records[
        binding["job_id"]
    ].replace("|COMPLETED|", "|FAILED|")
    with pytest.raises(ValueError, match="differs from live Slurm"):
        qualification.select_profile(evidence["R16"], scheduler)


def test_scheduler_submit_line_fabrication_is_rejected(qualification_fixture):
    wave = qualification_fixture["build_wave"](("R16",))
    wave_path = qualification.materialize_prepared_wave(wave)
    submit_wave(wave, wave_path)
    for packet in packet_map(wave).values():
        retain_outputs(packet, wave["provenance"])
    scheduler = LiveScheduler(wave)
    first = next(iter(scheduler.records))
    scheduler.records[first] = scheduler.records[first].replace(
        "/usr/bin/sbatch", "/bin/true"
    )
    with pytest.raises(ValueError, match="authenticated job binding"):
        qualification.retain_wave_audit(wave_path, scheduler)


def test_selection_reauthenticates_live_scheduler_stored_script(
    qualification_fixture,
):
    _, _, evidence, scheduler = retain_complete_wave(
        qualification_fixture, ("R16",)
    )
    job_id = next(iter(scheduler.scripts))
    scheduler.scripts[job_id] += b"\n# live scheduler mutation\n"
    with pytest.raises(ValueError, match="Slurm-stored batch script"):
        qualification.select_profile(evidence["R16"], scheduler)


def test_retained_post_submission_queue_rejects_array_rebinding(
    qualification_fixture,
):
    wave = qualification_fixture["build_wave"](("R16",))
    wave_path = qualification.materialize_prepared_wave(wave)
    submit_wave(wave, wave_path)
    packet = packet_map(wave)[("R16", 1)]
    binding = qualification.load_json(
        Path(packet["execution_intent"]["paths"]["job_binding"]), "binding"
    )
    record = deepcopy(binding["post_submission_queue"])
    record["live_jobs"][0]["array_task_id"] = "1"
    with pytest.raises(ValueError, match="condensed or ambiguous array"):
        qualification.validate_post_submission_queue(
            record,
            ["1000"],
            expected_wave_nodes=3,
            expected_current_job=("1000", packet["execution_intent"]["job_name"], 1),
        )


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"nonfinite_snapshot": True}, "non-finite values"),
        ({"wrong_variables": True}, "required variable inventory"),
        ({"duplicate_logical": True}, "logical meshblock is invalid"),
        ({"missing_logical": True}, "logical meshblock inventory"),
    ],
)
def test_snapshot_payload_inventory_and_logical_mesh_fail_closed(
    qualification_fixture, kwargs, message,
):
    wave, packet = prepared_packet(qualification_fixture)
    retain_outputs(packet, wave["provenance"], **kwargs)
    with pytest.raises(ValueError, match=message):
        qualification.inspect_pilot_output(
            packet, Path(wave["qualification_root"]), wave["provenance"]
        )


def test_r17_requires_exactly_27_meshblocks_per_rank(
    qualification_fixture, monkeypatch,
):
    wave = build_r17_wave(qualification_fixture, monkeypatch)
    wave_path = qualification.materialize_prepared_wave(wave)
    submit_wave(wave, wave_path)
    packet = packet_map(wave)[("R17", 8)]
    retain_outputs(
        packet, wave["provenance"], imbalanced_meshblocks=True
    )
    with pytest.raises(ValueError, match="exactly 27 meshblocks per rank"):
        qualification.inspect_pilot_output(
            packet, Path(wave["qualification_root"]), wave["provenance"]
        )


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"nonpositive_snapshot": "dens"}, "nonpositive dens"),
        ({"nonpositive_snapshot": "eint"}, "nonpositive eint"),
        ({"nonpositive_snapshot": "p_perp"}, "nonpositive p_perp"),
        ({"bad_snapshot_nghost": True}, "reviewed contract"),
        ({"bad_snapshot_geometry": True}, "metadata or geometry"),
        ({"bad_snapshot_preheader": True}, "preheader metadata differs"),
        ({"bad_snapshot_location_size": True}, "snapshot metadata differs"),
        ({"bad_snapshot_variable_size": True}, "snapshot metadata differs"),
        ({"bad_snapshot_active_indices": True}, "metadata or geometry"),
        ({"missing_snapshot_cadence": True}, "exact output cadence"),
        (
            {"hard_bound_violation": True},
            "independently violates CGL hard bounds",
        ),
    ],
)
def test_snapshot_positivity_exact_metadata_geometry_and_cadence_are_required(
    qualification_fixture, kwargs, message,
):
    wave, packet = prepared_packet(qualification_fixture)
    retain_outputs(packet, wave["provenance"], **kwargs)
    with pytest.raises(ValueError, match=message):
        qualification.inspect_pilot_output(
            packet, Path(wave["qualification_root"]), wave["provenance"]
        )


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"tiny_restart": True}, "completeness floor"),
        ({"restart_config_drift": True}, "reviewed contract"),
        ({"restart_binary_time": 0.5}, "marker does not authenticate"),
        ({"smoke": "missing"}, "restart smoke"),
        ({"smoke": "failed"}, "smoke result differs"),
        ({"smoke": "wrong"}, "smoke result differs"),
    ],
)
def test_restart_completeness_abi_configuration_and_load_smoke_are_required(
    qualification_fixture, kwargs, message,
):
    wave, packet = prepared_packet(qualification_fixture)
    retain_outputs(packet, wave["provenance"], **kwargs)
    with pytest.raises(ValueError, match=message):
        qualification.inspect_pilot_output(
            packet, Path(wave["qualification_root"]), wave["provenance"]
        )


def test_unqualified_restart_abi_is_rejected(qualification_fixture, monkeypatch):
    wave, packet = prepared_packet(qualification_fixture)
    retain_outputs(packet, wave["provenance"])
    monkeypatch.setattr(qualification, "QUALIFIED_RESTART_BINARY_ABIS", {})
    with pytest.raises(ValueError, match="no qualified binary-restart ABI"):
        qualification.inspect_pilot_output(
            packet, Path(wave["qualification_root"]), wave["provenance"]
        )


def test_environment_log_must_bind_job_contract_and_smoke_interval(
    qualification_fixture,
):
    wave, packet = prepared_packet(qualification_fixture)
    retain_outputs(packet, wave["provenance"])
    environment = Path(packet["execution_intent"]["paths"]["environment_log"])
    environment.write_text(environment.read_text().replace(
        "slurm_job_id=1000", "slurm_job_id=9999"
    ))
    with pytest.raises(ValueError, match="smoke result differs"):
        qualification.inspect_pilot_output(
            packet, Path(wave["qualification_root"]), wave["provenance"]
        )


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"bad_start": True}, "exact start/cadence"),
        ({"bad_cadence": True}, "exact start/cadence"),
        ({"nonmonotonic": True}, "exact start/cadence"),
        ({"no_work": True}, "nontrivial forcing or pressure work"),
        ({"tail_noop": True}, "nontrivial forcing or pressure work"),
    ],
)
def test_history_start_cadence_monotonicity_and_nontrivial_work_are_required(
    qualification_fixture, kwargs, message,
):
    wave, packet = prepared_packet(qualification_fixture)
    retain_outputs(packet, wave["provenance"], **kwargs)
    with pytest.raises(ValueError, match=message):
        qualification.inspect_pilot_output(
            packet, Path(wave["qualification_root"]), wave["provenance"]
        )


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"cap_increment_violation": True}, "cap increment exceeds interval"),
        ({"near_noop": True}, "nontrivial forcing or pressure work"),
        ({"bad_divb": True}, "normalized CT-divB"),
        ({"missing_max_ndiv": True}, "lacks required columns"),
        ({"fractional_count": True}, "count history is not integral"),
    ],
)
def test_interval_cap_work_floor_and_normalized_divb_gates_are_required(
    qualification_fixture, kwargs, message,
):
    wave, packet = prepared_packet(qualification_fixture)
    retain_outputs(packet, wave["provenance"], **kwargs)
    with pytest.raises(ValueError, match=message):
        qualification.inspect_pilot_output(
            packet, Path(wave["qualification_root"]), wave["provenance"]
        )


@pytest.mark.parametrize("link_kind", ["symlink", "hardlink"])
def test_rank_product_links_are_rejected(qualification_fixture, link_kind):
    wave, packet = prepared_packet(qualification_fixture)
    retain_outputs(packet, wave["provenance"])
    product = (
        Path(packet["execution_intent"]["paths"]["output_dir"])
        / "bin/rank_00000000/snapshot.00001.bin"
    )
    external = Path(packet["execution_intent"]["paths"]["run_dir"]) / "external.bin"
    product.rename(external)
    if link_kind == "symlink":
        product.symlink_to(external)
        message = "traverses a symlink"
    else:
        product.hardlink_to(external)
        message = "exactly one filesystem link"
    with pytest.raises(ValueError, match=message):
        qualification.inspect_pilot_output(
            packet, Path(wave["qualification_root"]), wave["provenance"]
        )


def test_descriptor_rebinding_race_is_rejected(
    qualification_fixture, monkeypatch,
):
    wave, packet = prepared_packet(qualification_fixture)
    retain_outputs(packet, wave["provenance"])
    product = (
        Path(packet["execution_intent"]["paths"]["output_dir"])
        / "bin/rank_00000000/snapshot.00000.bin"
    )
    original = qualification.structural_snapshot_profile
    raced = False

    def race(path, intent=None):
        nonlocal raced
        result = original(path, intent)
        if not raced:
            raced = True
            payload = product.read_bytes()
            product.rename(product.with_suffix(".old"))
            product.write_bytes(payload)
        return result

    monkeypatch.setattr(qualification, "structural_snapshot_profile", race)
    with pytest.raises(ValueError, match="changed during inspection"):
        qualification.inspect_pilot_output(
            packet, Path(wave["qualification_root"]), wave["provenance"]
        )


def test_dependency_wave_overlap_is_rejected(qualification_fixture):
    wave = qualification_fixture["build_wave"](("R04", "R12", "R16"))
    wave_path = qualification.materialize_prepared_wave(wave)
    submit_wave(wave, wave_path)
    for packet in packet_map(wave).values():
        retain_outputs(packet, wave["provenance"])
    with pytest.raises(ValueError, match="waves overlap"):
        qualification.retain_wave_audit(wave_path, LiveScheduler(wave, overlap=True))


def test_r17_operational_evidence_retains_decomposition_and_interval_exclusivity(
    qualification_fixture, monkeypatch,
):
    monkeypatch.setattr(qualification, "RANKS_PER_NODE", 8)
    wave, _, evidence, scheduler = retain_complete_wave(
        qualification_fixture, ("R17",), max_nodes=8
    )
    packet = packet_map(wave)[("R17", 8)]
    scientific = qualification.load_json(
        Path(packet["execution_intent"]["paths"]["scientific_evidence"]),
        "R17 scientific evidence",
    )
    assert scientific["accepted_for_profile_selection"] is False
    assert scientific["accepted_for_operational_qualification"] is True
    assert len(scientific["terminal_rank_local_outputs"]) == 64
    assert len(scientific["terminal_rank_local_restarts"]) == 64
    assert scientific["physics_measurements"]["normalized_ct_divb_below_threshold"] is True
    assert qualification.retained_inventory_sha256(
        scientific["terminal_rank_local_outputs"]
    ) == scientific["terminal_rank_local_output_inventory_sha256"]
    assert qualification.retained_inventory_sha256(
        scientific["terminal_rank_local_restarts"]
    ) == scientific["terminal_rank_local_restart_inventory_sha256"]
    decomposition = scientific["r17_decomposition"]
    assert decomposition["resolution"] == "384x384x768"
    assert decomposition["mesh_shape"] == [384, 384, 768]
    assert decomposition["meshblock_shape"] == [32, 32, 64]
    assert decomposition["logical_meshblocks"] == 1728
    assert decomposition["ranks"] == 64
    assert decomposition["meshblocks_per_rank"] == 27
    assert len(decomposition["complete_block_rank_inventory"]) == 64
    assert all(
        len(record["logical_meshblocks"]) == 27
        for record in decomposition["complete_block_rank_inventory"]
    )
    assert qualification.stable_json_sha256(
        decomposition["complete_block_rank_inventory"]
    ) == decomposition["complete_block_rank_inventory_sha256"]
    assert qualification.validate_r17_decomposition_evidence(
        decomposition, scientific["terminal_rank_local_output_inventory_sha256"]
    ) == decomposition

    retained = qualification.retain_r17_operational_qualification(
        evidence["R17"], "qualification measurement agent", scheduler
    )
    assert retained["schema_version"] == 2
    assert retained["independent_review_created"] is False
    assert retained["self_approved"] is False
    assert retained["launch_authorized"] is False
    root = qualification_fixture["root"]
    qualification_path = root / retained["operational_qualification"]["path"]
    operational = qualification.load_json(
        qualification_path, "R17 operational qualification"
    )
    assert operational["record_type"] == "stage-i-r17-operational-qualification"
    assert operational["schema_version"] == 2
    assert set(operational) == {
        "schema_version",
        "record_type",
        "execution_epoch",
        "completed_utc",
        "measured_utc",
        "measured_by",
        "job_id",
        "state",
        "exit_code",
        "nodes",
        "ranks",
        "prepared_wave",
        "qualification_evidence",
        "scientific_evidence",
        "scheduler_evidence",
        "account_scheduler_evidence",
        "account_exclusivity_evidence",
        "executable_sha256",
        "build_manifest_inventory",
        "build_manifest_inventory_sha256",
        "rank_local_outputs",
        "rank_local_output_inventory_sha256",
        "rank_local_restarts",
        "rank_local_restart_inventory_sha256",
        "decomposition_evidence",
        "restart_load_evidence",
        "physics_validation_evidence",
        "frozen_science_build_contract",
        "independent_review_contract",
        "authority",
    }
    assert operational["nodes"] == 8
    assert operational["ranks"] == 64
    assert operational["measured_by"] == "qualification measurement agent"
    assert operational["decomposition_evidence"] == decomposition
    assert operational["build_manifest_inventory"] == (
        wave["provenance"]["build_manifest"]["inventory"]
    )
    assert operational["build_manifest_inventory_sha256"] == (
        wave["provenance"]["build_manifest"]["inventory_sha256"]
    )
    assert len(operational["rank_local_outputs"]) == 64
    assert len(operational["rank_local_restarts"]) == 64
    assert operational["rank_local_output_inventory_sha256"] == (
        scientific["terminal_rank_local_output_inventory_sha256"]
    )
    assert operational["rank_local_restart_inventory_sha256"] == (
        scientific["terminal_rank_local_restart_inventory_sha256"]
    )
    frozen = operational["frozen_science_build_contract"]
    assert frozen["case_id"] == "R17"
    assert frozen["resolution"] == "384x384x768"
    assert frozen["mesh_shape"] == [384, 384, 768]
    assert frozen["meshblock_shape"] == [32, 32, 64]
    assert frozen["source_revision"] == qualification.FROZEN_SOURCE_REVISION
    assert frozen["matrix_sha256"] == qualification.FROZEN_MATRIX_SHA256
    assert frozen["input_sha256"] == qualification.CASE_POLICIES["R17"]["input_sha256"]
    assert frozen["execution_intent_sha256"] == packet["execution_intent_sha256"]
    assert frozen["execution_contract_sha256"] == (
        packet["execution_intent"]["execution_contract_sha256"]
    )
    assert frozen["executable_sha256"] == operational["executable_sha256"]
    assert frozen["build_manifest_inventory_sha256"] == (
        operational["build_manifest_inventory_sha256"]
    )
    assert qualification.stable_json_sha256(frozen["parameter_contract"]) == (
        frozen["parameter_contract_sha256"]
    )
    review_contract = operational["independent_review_contract"]
    assert review_contract["required"] is True
    assert review_contract["schema_version"] == 1
    assert review_contract["record_type"] == (
        "stage-i-r17-operational-qualification-independent-review"
    )
    assert review_contract["decision"] == "approved"
    assert review_contract["candidate_path"] == str(qualification_path)
    assert review_contract["candidate_sha256_required"] is True
    assert review_contract["reviewer_must_differ_from"] == [
        "qualification measurement agent"
    ]
    assert review_contract["reviewed_after_utc"] == operational["measured_utc"]
    assert retained["required_independent_review"] == review_contract
    assert retained["authority"] == operational["authority"] == {
        "r17_launch_authorized": False,
        "scheduler_mutation_authorized": False,
        "canonical_mutation_authorized": False,
    }
    review_path = root / review_contract["path"]
    assert stat_mode(qualification_path) == 0o444
    assert not review_path.exists()
    restart_path = root / operational["restart_load_evidence"]["path"]
    physics_path = root / operational["physics_validation_evidence"]["path"]
    scheduler_path = root / operational["scheduler_evidence"]["path"]
    account_scheduler_path = root / operational["account_scheduler_evidence"]["path"]
    account_exclusivity_path = (
        root / operational["account_exclusivity_evidence"]["path"]
    )
    assert stat_mode(restart_path) == 0o444
    assert stat_mode(physics_path) == 0o444
    assert stat_mode(scheduler_path) == 0o444
    assert stat_mode(account_scheduler_path) == 0o444
    assert stat_mode(account_exclusivity_path) == 0o444
    account_exclusivity = qualification.load_json(
        account_exclusivity_path, "R17 account exclusivity"
    )
    assert account_exclusivity["exclusive_entire_execution_interval"] is True
    assert account_exclusivity["overlapping_job_ids"] == [operational["job_id"]]
    assert account_exclusivity["query_contract"]["all_users"] is True
    assert account_exclusivity["query_contract"]["allocations_only"] is True
    assert account_exclusivity["query_contract"]["expanded_arrays"] is True
    assert account_exclusivity["visibility_contract"] == {
        "private_data": "none",
        "all_users_job_visibility": True,
    }
    assert qualification.sha256(account_scheduler_path) == (
        account_exclusivity["raw_account_scheduler_sha256"]
    )
    assert qualification.stable_json_sha256(account_exclusivity["account_jobs"]) == (
        account_exclusivity["account_jobs_sha256"]
    )
    scheduler_evidence = qualification.load_json(
        Path(packet["execution_intent"]["paths"]["scheduler_evidence"]),
        "R17 scheduler evidence",
    )
    assert qualification.validate_r17_account_exclusivity_evidence(
        account_exclusivity,
        scheduler_evidence,
        account_scheduler_path.read_text(),
    ) == account_exclusivity
    assert qualification.validate_r17_operational_qualification_contract(
        operational, root, qualification_path
    ) == operational
    fabricated_account_digest = deepcopy(account_exclusivity)
    fabricated_account_digest["account_jobs_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="exclusivity evidence differs"):
        qualification.validate_r17_account_exclusivity_evidence(
            fabricated_account_digest,
            scheduler_evidence,
            account_scheduler_path.read_text(),
        )
    physics = qualification.load_json(physics_path, "R17 physics validation")
    assert set(physics["measurements"]) == {
        "rank_local_output_inventory_sha256",
        "finite_rank_outputs",
        "mass_relative_drift_max",
        "mhd_user_mass_mismatch_max",
        "lf_bad_counts_total",
        "normalized_ct_divb_max",
        "normalized_ct_divb_threshold",
        "normalized_ct_divb_below_threshold",
    }
    with pytest.raises(ValueError, match="not profile-selection evidence"):
        qualification.select_profile(evidence["R17"], scheduler)


def test_r17_operational_contract_rejects_noncanonical_schema_and_binding_drift(
    qualification_fixture, monkeypatch,
):
    monkeypatch.setattr(qualification, "RANKS_PER_NODE", 8)
    _, _, evidence, scheduler = retain_complete_wave(
        qualification_fixture, ("R17",), max_nodes=8
    )
    retained = qualification.retain_r17_operational_qualification(
        evidence["R17"], "qualification measurement agent", scheduler
    )
    root = qualification_fixture["root"]
    qualification_path = root / retained["operational_qualification"]["path"]
    canonical = qualification.load_json(
        qualification_path, "R17 operational qualification"
    )

    def duplicate_output_rank(value):
        value["rank_local_outputs"][0] = deepcopy(value["rank_local_outputs"][1])
        value["rank_local_outputs"].sort(key=lambda item: item["path"])
        value["rank_local_output_inventory_sha256"] = (
            qualification.retained_inventory_sha256(value["rank_local_outputs"])
        )

    def drop_restart(value):
        value["rank_local_restarts"].pop()
        value["rank_local_restart_inventory_sha256"] = (
            qualification.retained_inventory_sha256(value["rank_local_restarts"])
        )

    mutations = (
        lambda value: value.__setitem__("schema_version", 1),
        lambda value: value.__setitem__("legacy_schema_1_fallback", {}),
        lambda value: value.__setitem__(
            "rank_local_output_inventory_sha256", "0" * 64
        ),
        duplicate_output_rank,
        drop_restart,
        lambda value: value["frozen_science_build_contract"].__setitem__(
            "input_sha256", "0" * 64
        ),
        lambda value: value["independent_review_contract"].__setitem__(
            "candidate_path", str(qualification_path) + ".different"
        ),
        lambda value: value["authority"].__setitem__(
            "r17_launch_authorized", True
        ),
        lambda value: value.__setitem__(
            "account_exclusivity_evidence", deepcopy(value["scheduler_evidence"])
        ),
    )
    for mutate in mutations:
        adversarial = deepcopy(canonical)
        mutate(adversarial)
        with pytest.raises(ValueError):
            qualification.validate_r17_operational_qualification_contract(
                adversarial, root, qualification_path
            )


def test_r17_canonical_publication_rejects_rebound_stage_i_lock_before_write(
    qualification_fixture, monkeypatch,
):
    monkeypatch.setattr(qualification, "RANKS_PER_NODE", 8)
    _, _, evidence, scheduler = retain_complete_wave(
        qualification_fixture, ("R17",), max_nodes=8
    )
    root = qualification_fixture["root"]
    lock_path = root / qualification.STAGE_I_LOCK_NAME
    displaced = root / f"{qualification.STAGE_I_LOCK_NAME}.publication-displaced"
    original = qualification.require_empty_canonical_stage_i_state

    def rebind_after_admission(project_root):
        state = original(project_root)
        lock_path.rename(displaced)
        lock_path.write_bytes(b"")
        lock_path.chmod(0o644)
        return state

    monkeypatch.setattr(
        qualification,
        "require_empty_canonical_stage_i_state",
        rebind_after_admission,
    )
    with pytest.raises(ValueError, match="mutation lock.*changed"):
        qualification.retain_r17_operational_qualification(
            evidence["R17"], "qualification measurement agent", scheduler
        )
    assert not any(
        "R17_operational_qualification" in path.name
        for path in (root / "accounting").iterdir()
    )


def test_r17_decomposition_evidence_rejects_digest_and_rank_map_fabrication(
    qualification_fixture, monkeypatch,
):
    monkeypatch.setattr(qualification, "RANKS_PER_NODE", 8)
    wave, _, _, _ = retain_complete_wave(
        qualification_fixture, ("R17",), max_nodes=8
    )
    packet = packet_map(wave)[("R17", 8)]
    scientific = qualification.load_json(
        Path(packet["execution_intent"]["paths"]["scientific_evidence"]),
        "R17 scientific evidence",
    )
    decomposition = scientific["r17_decomposition"]
    digest_drift = deepcopy(decomposition)
    digest_drift["complete_block_rank_inventory_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="reviewed contract"):
        qualification.validate_r17_decomposition_evidence(
            digest_drift, scientific["terminal_rank_local_output_inventory_sha256"]
        )
    duplicate = deepcopy(decomposition)
    duplicate["complete_block_rank_inventory"][0]["logical_meshblocks"][0] = (
        duplicate["complete_block_rank_inventory"][0]["logical_meshblocks"][1]
    )
    duplicate["complete_block_rank_inventory"][0]["logical_meshblocks"].sort()
    duplicate["complete_block_rank_inventory_sha256"] = qualification.stable_json_sha256(
        duplicate["complete_block_rank_inventory"]
    )
    with pytest.raises(ValueError, match="27 unique blocks per rank"):
        qualification.validate_r17_decomposition_evidence(
            duplicate, scientific["terminal_rank_local_output_inventory_sha256"]
        )


def test_r17_operational_retention_rejects_full_interval_account_overlap(
    qualification_fixture, monkeypatch,
):
    monkeypatch.setattr(qualification, "RANKS_PER_NODE", 8)
    wave, _, evidence, scheduler = retain_complete_wave(
        qualification_fixture, ("R17",), max_nodes=8
    )
    job_id = next(iter(scheduler.account_records))
    target = scheduler.account_records[job_id].strip().split("|")

    scheduler.account_records.pop(job_id)
    with pytest.raises(ValueError, match="account scheduler raw evidence is empty"):
        qualification.retain_r17_operational_qualification(
            evidence["R17"], "qualification measurement agent", scheduler
        )

    scheduler.account_records[job_id] = "|".join(target) + "\n"
    overlap = list(target)
    overlap[0] = "900"
    overlap[1] = "foreign_account_job"
    overlap[4] = "1"
    overlap[-1] = "another-account-user"
    scheduler.account_records["900"] = "|".join(overlap) + "\n"
    with pytest.raises(ValueError, match="overlapping account job"):
        qualification.retain_r17_operational_qualification(
            evidence["R17"], "qualification measurement agent", scheduler
        )

    overlap[0] = "901"
    overlap[2] = "RUNNING"
    overlap[3] = "0:0"
    overlap[8] = "Unknown"
    scheduler.account_records = {
        job_id: "|".join(target) + "\n",
        "901": "|".join(overlap) + "\n",
    }
    with pytest.raises(ValueError, match="overlapping account job"):
        qualification.retain_r17_operational_qualification(
            evidence["R17"], "qualification measurement agent", scheduler
        )

    scheduler.account_records = {job_id: "|".join(target) + "\n"}
    scheduler.private_data = "jobs"
    with pytest.raises(ValueError, match="prevents account-wide job evidence"):
        qualification.retain_r17_operational_qualification(
            evidence["R17"], "qualification measurement agent", scheduler
        )
    assert not any(
        "R17_operational_qualification" in path.name
        for path in (qualification_fixture["root"] / "accounting").iterdir()
    )


def test_r17_account_interval_query_preserves_scheduler_timezone():
    scheduler = {
        "job_id": "123",
        "job_name": "cglq_R17_n08",
        "state": "COMPLETED",
        "exit_code": "0:0",
        "nodes": 8,
        "elapsed_seconds": 60,
        "submit_utc": "2026-06-05T03:00:00-0400",
        "start_utc": "2026-06-05T03:10:00-0400",
        "end_utc": "2026-06-05T03:11:00-0400",
        "partition": qualification.PARTITION,
        "account": qualification.ACCOUNT.casefold(),
    }
    contract = qualification.account_scheduler_query_contract(
        qualification.parse_scheduler_time(scheduler["start_utc"], "start"),
        qualification.parse_scheduler_time(scheduler["end_utc"], "end"),
    )
    assert contract["start_argument"] == "2026-06-05T03:09:59"
    assert contract["end_argument"] == "2026-06-05T03:11:01"
    assert contract["start_utc"] == "2026-06-05T07:09:59+00:00"
    raw = qualification.ACCOUNT_SCHEDULER_HEADER + "\n" + "|".join([
        scheduler["job_id"], scheduler["job_name"], scheduler["state"],
        scheduler["exit_code"], str(scheduler["nodes"]),
        str(scheduler["elapsed_seconds"]), scheduler["submit_utc"],
        scheduler["start_utc"], scheduler["end_utc"], scheduler["partition"],
        scheduler["account"], "producer",
    ]) + "\n"
    evidence = qualification.build_r17_account_exclusivity_evidence(
        scheduler,
        raw,
        contract,
        {"private_data": "none", "all_users_job_visibility": True},
        "2026-06-05T07:12:00+00:00",
    )
    assert evidence["exclusive_entire_execution_interval"] is True


def test_audit_requires_job_bound_immutable_output_tree(qualification_fixture):
    wave = qualification_fixture["build_wave"](("R16",))
    wave_path = qualification.materialize_prepared_wave(wave)
    submit_wave(wave, wave_path)
    for packet in packet_map(wave).values():
        retain_outputs(packet, wave["provenance"])
    scheduler = LiveScheduler(wave)
    packet = packet_map(wave)[("R16", 1)]
    binding_path = Path(packet["execution_intent"]["paths"]["output_binding"])
    binding_path.unlink()
    with pytest.raises(ValueError, match="pilot output binding is missing"):
        qualification.retain_wave_audit(wave_path, scheduler)


def test_audit_rejects_output_mutated_after_completed_job_binding(
    qualification_fixture,
):
    wave = qualification_fixture["build_wave"](("R16",))
    wave_path = qualification.materialize_prepared_wave(wave)
    submit_wave(wave, wave_path)
    for packet in packet_map(wave).values():
        retain_outputs(packet, wave["provenance"])
    scheduler = LiveScheduler(wave)
    packet = packet_map(wave)[("R16", 1)]
    output = Path(packet["execution_intent"]["paths"]["output_dir"])
    history = next(output.glob("*.user.hst"))
    history.write_bytes(history.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="differs from live output tree"):
        qualification.retain_wave_audit(wave_path, scheduler)


def test_audit_rejects_output_binding_outside_completed_job_interval(
    qualification_fixture,
):
    wave = qualification_fixture["build_wave"](("R16",))
    wave_path = qualification.materialize_prepared_wave(wave)
    submit_wave(wave, wave_path)
    for packet in packet_map(wave).values():
        retain_outputs(packet, wave["provenance"])
    scheduler = LiveScheduler(wave)
    packet = packet_map(wave)[("R16", 1)]
    binding_path = Path(packet["execution_intent"]["paths"]["output_binding"])
    binding = json.loads(binding_path.read_text())
    binding["completed_utc"] = "2000-01-01T00:00:00Z"
    os.chmod(binding_path, 0o640)
    write_json(binding_path, binding)
    os.chmod(binding_path, 0o440)
    with pytest.raises(ValueError, match="outside Slurm job"):
        qualification.retain_wave_audit(wave_path, scheduler)


@pytest.mark.parametrize(
    ("runtimes", "selected"),
    [
        ({1: 3000, 2: 2450, 4: 1200}, 1),
        ({1: 3000, 2: 1300, 4: 1200}, 2),
        ({1: 4000, 2: 1800, 4: 1000}, 4),
        ({1: 4000, 2: 1000, 4: 950}, 2),
    ],
)
def test_selection_math_requires_complete_live_reproduced_audit(
    qualification_fixture, runtimes, selected,
):
    mapped = {("R04", nodes): runtime for nodes, runtime in runtimes.items()}
    _, _, evidence, scheduler = retain_complete_wave(
        qualification_fixture, ("R04",), runtimes=mapped
    )
    result = qualification.select_profile(evidence["R04"], scheduler)
    assert result["selected_nodes"] == selected


def test_selection_requires_cross_profile_scientific_agreement(
    qualification_fixture,
):
    _, _, evidence, scheduler = retain_complete_wave(
        qualification_fixture,
        ("R16",),
        output_options={("R16", 2): {"work_scale": 1.01}},
    )
    with pytest.raises(ValueError, match="cross-profile scientific agreement failed"):
        qualification.select_profile(evidence["R16"], scheduler)


def test_r12_selection_retains_exact_fixed_fresh_profile(
    qualification_fixture,
):
    _, _, evidence, scheduler = retain_complete_wave(
        qualification_fixture, ("R12",)
    )
    selection = qualification.select_profile(evidence["R12"], scheduler)
    assert selection["selected_nodes"] == 4
    assert selection["target_time"] == 0.12
    assert selection["selection_reason"] == (
        "selected exact measured parentless fresh R12 profile"
    )
    assert selection["cross_profile_scientific_agreement"] == {
        "schema_version": 1,
        "record_type": "cgl_lf_stage_i_cross_profile_scientific_agreement",
        "baseline_nodes": 4,
        "relative_tolerance": qualification.CROSS_PROFILE_RELATIVE_TOLERANCE,
        "absolute_tolerance": qualification.CROSS_PROFILE_ABSOLUTE_TOLERANCE,
        "comparisons": [],
        "all_profiles_agree": True,
    }
    assert selection["consumption_contract"]["fixed_fresh_profile"] == (
        qualification.R12_FRESH_RERUN_PROFILE
    )
    assert selection["consumption_contract"]["selected_profile"]["nodes"] == 4


@pytest.mark.parametrize("drift", ["source", "bundle", "build", "helper"])
def test_selection_revalidates_all_live_provenance(qualification_fixture, drift):
    _, _, evidence, scheduler = retain_complete_wave(
        qualification_fixture, ("R16",)
    )
    if drift == "source":
        path = (
            qualification_fixture["source"]
            / qualification.CASE_POLICIES["R16"]["input_relative_path"]
        )
        path.write_bytes(path.read_bytes() + b"\n")
    elif drift == "bundle":
        qualification_fixture["bundle"].write_bytes(b"fabricated bundle\n")
    elif drift == "build":
        (qualification_fixture["build_manifest"] / "environment.txt").write_text(
            "git_revision=" + "0" * 40 + "\n"
        )
    else:
        qualification_fixture["helper"].write_bytes(
            qualification_fixture["helper"].read_bytes() + b"\n"
        )
    with pytest.raises(ValueError):
        qualification.select_profile(evidence["R16"], scheduler)


def test_selection_rejects_ambiguous_profile_results(qualification_fixture):
    _, _, evidence, scheduler = retain_complete_wave(
        qualification_fixture, ("R16",)
    )
    value = json.loads(evidence["R16"].read_text())
    value["results"][1] = deepcopy(value["results"][0])
    write_json(evidence["R16"], value)
    with pytest.raises(ValueError, match="differs from reproduced wave audit"):
        qualification.select_profile(evidence["R16"], scheduler)


def test_selection_rejects_stale_scheduler_evidence(qualification_fixture):
    _, _, evidence, scheduler = retain_complete_wave(
        qualification_fixture, ("R16",)
    )
    future = datetime.now(timezone.utc) + timedelta(hours=13)
    with pytest.raises(ValueError, match="scheduler evidence is stale"):
        qualification.select_profile(
            evidence["R16"], scheduler, now=future
        )


def test_durable_selection_is_authenticated_expiring_and_non_authorizing(
    qualification_fixture,
):
    _, _, evidence, scheduler = retain_complete_wave(
        qualification_fixture, ("R16",)
    )
    generated = datetime.now(timezone.utc).replace(microsecond=0)
    retained = qualification.retain_profile_selection(
        evidence["R16"], scheduler, now=generated
    )
    selection = retained["selection"]
    selection_path = Path(retained["retained_selection"]["path"])
    contract = selection["consumption_contract"]
    assert selection_path.stat().st_mode & 0o222 == 0
    assert contract["advisory_only"] is True
    assert contract["canonical_acceptance_eligible"] is False
    assert contract["production_authorization"] is False
    assert contract["production_controller_consumption_implemented"] is False
    assert contract["reviewed_promotion_required"] is True
    assert contract["cross_profile_scientific_agreement_required"] is True
    assert selection["cross_profile_scientific_agreement"]["all_profiles_agree"] is True
    assert qualification.validate_profile_selection(
        selection_path, scheduler, now=generated + timedelta(minutes=1)
    ) == selection
    with pytest.raises(FileExistsError):
        qualification.retain_profile_selection(
            evidence["R16"], scheduler, now=generated
        )
    with pytest.raises(ValueError, match="expired or future-dated"):
        qualification.validate_profile_selection(
            selection_path, scheduler, now=generated + timedelta(hours=7)
        )
    os.chmod(selection_path, 0o640)
    with pytest.raises(ValueError, match="selection path differs"):
        qualification.validate_profile_selection(
            selection_path, scheduler, now=generated + timedelta(minutes=1)
        )


def test_durable_selection_tampering_is_rejected(qualification_fixture):
    _, _, evidence, scheduler = retain_complete_wave(
        qualification_fixture, ("R16",)
    )
    retained = qualification.retain_profile_selection(evidence["R16"], scheduler)
    path = Path(retained["retained_selection"]["path"])
    value = json.loads(path.read_text())
    value["selected_nodes"] = 99
    os.chmod(path, 0o640)
    write_json(path, value)
    os.chmod(path, 0o440)
    with pytest.raises(ValueError, match="selection digest differs"):
        qualification.validate_profile_selection(path, scheduler)
