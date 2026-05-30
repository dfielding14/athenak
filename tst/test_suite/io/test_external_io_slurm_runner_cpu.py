"""Local mock coverage for the external Slurm IO qualification runner."""

from pathlib import Path
import os
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[3]
RUNNER = ROOT / "scripts" / "run_external_io_qualification_slurm.sh"


def _write_executable(path: Path, payload: str):
    path.write_text(payload)
    path.chmod(0o755)


def _mock_environment(tmp_path: Path):
    bin_dir = tmp_path / "bin"
    repo = tmp_path / "repo"
    build = tmp_path / "build"
    bin_dir.mkdir()
    (repo / "tst" / "inputs").mkdir(parents=True)
    (build / "src").mkdir(parents=True)
    (repo / "tst" / "inputs" / "io_node_sharding.athinput").write_text("mock deck\n")
    _write_executable(build / "src" / "athena", "#!/usr/bin/env bash\nexit 0\n")
    _write_executable(
        bin_dir / "timeout",
        """#!/usr/bin/env bash
shift
if [[ "${MOCK_TIMEOUT_RANK_MAP:-0}" == 1 && " $* " == *" hostname "* ]]; then
  exit 124
fi
exec "$@"
""",
    )
    _write_executable(
        bin_dir / "srun",
        """#!/usr/bin/env bash
output_dir=
manifest=
previous=
for argument in "$@"; do
  if [[ "$previous" == -d ]]; then output_dir="$argument"; fi
  if [[ "$previous" == -r ]]; then manifest="$argument"; fi
  previous="$argument"
done
printf '%s\n' "${ATHENAK_RESTART_MANIFEST_TIMING:-unset}" >> "$MOCK_SRUN_LOG"
if [[ -n "$output_dir" ]]; then
  mkdir -p "$output_dir"
  if [[ "${MOCK_OUTPUT_ASSEMBLED:-0}" == 1 ]]; then
    : > "$output_dir/forbidden.assembled"
  fi
fi
if [[ -n "$manifest" && -n "${MOCK_MANIFEST_SIDECAR_SUFFIX:-}" ]]; then
  : > "$manifest${MOCK_MANIFEST_SIDECAR_SUFFIX}"
fi
exit "${MOCK_SRUN_EXIT:-0}"
""",
    )
    _write_executable(
        bin_dir / "hook",
        """#!/usr/bin/env bash
if [[ "${MOCK_HOOK_FAIL_PHASE:-}" == "$1" ]]; then
  exit 9
fi
printf 'hook phase=%s row=%s sample=%s output=%s\n' "$1" "$2" "$3" "$4"
""",
    )
    manifest = tmp_path / "manifest.rst"
    manifest.write_text("mock manifest\n")
    env = os.environ.copy()
    env.update(
        {
            "PATH": f"{bin_dir}:{env['PATH']}",
            "REPO": str(repo),
            "BUILD": str(build),
            "PYTHON": "/usr/bin/python3",
            "ROW_TIMEOUT": "10",
            "MOCK_SRUN_LOG": str(tmp_path / "srun.log"),
        }
    )
    return env, manifest, bin_dir / "hook"


def _run_runner(
    tmp_path: Path,
    env,
    row: str,
    action: str,
    *arguments: str,
    attempt: str = "attempt-1",
    check=True,
):
    packet = tmp_path / f"packet-{attempt}"
    run_env = env.copy()
    run_env.update({"PACKET": str(packet), "ATTEMPT_ID": attempt})
    completed = subprocess.run(
        ["/bin/bash", str(RUNNER), row, action, *arguments],
        check=check,
        capture_output=True,
        text=True,
        env=run_env,
        timeout=30,
    )
    return completed, packet


def _index_rows(packet: Path):
    lines = (packet / "packet-index.tsv").read_text().splitlines()
    keys = lines[0].split("\t")
    return [dict(zip(keys, line.split("\t"))) for line in lines[1:]]


def test_slurm_runner_is_bash3_portable_and_rejects_duplicate_launches(tmp_path):
    env, _, _ = _mock_environment(tmp_path)
    _run_runner(tmp_path, env, "ED-1", "generate", attempt="first")
    duplicate, packet = _run_runner(
        tmp_path, env, "ED-1", "generate", attempt="first", check=False
    )

    assert duplicate.returncode == 2
    rows = _index_rows(packet)
    assert [row["sample"] for row in rows] == ["rank-map", "run"]
    assert all(int(row["monotonic_elapsed_ns"]) >= 0 for row in rows)

    hostfile = tmp_path / "mr2.hostfile"
    hostfile.write_text("host-a\nhost-b\nhost-b\n")
    env["MR2_HOSTFILE"] = str(hostfile)
    _run_runner(tmp_path, env, "MR-2", "generate", attempt="portable")


def test_slurm_runner_resume_has_warmup_five_measured_samples_and_timers(tmp_path):
    env, manifest, hook = _mock_environment(tmp_path)
    env["FS_ACCOUNTING_HOOK"] = str(hook)
    _, packet = _run_runner(
        tmp_path, env, "MR-1", "resume", str(manifest), attempt="resume"
    )

    rows = _index_rows(packet)
    assert [row["sample"] for row in rows] == [
        "rank-map",
        "warmup",
        "sample-1",
        "sample-2",
        "sample-3",
        "sample-4",
        "sample-5",
    ]
    assert [row["measured"] for row in rows] == ["0", "0", "1", "1", "1", "1", "1"]
    assert all(int(row["monotonic_elapsed_ns"]) >= 0 for row in rows)
    assert (tmp_path / "srun.log").read_text().splitlines() == [
        "unset",
        "0",
        "1",
        "1",
        "1",
        "1",
        "1",
    ]


@pytest.mark.parametrize(
    ("suffix", "expected"),
    (
        (".assembled", "manifest.rst.assembled"),
        (".assembled.tmp", "manifest.rst.assembled.tmp"),
    ),
)
def test_slurm_runner_rejects_manifest_sidecar_staging(tmp_path, suffix, expected):
    env, manifest, hook = _mock_environment(tmp_path)
    env.update(
        {
            "FS_ACCOUNTING_HOOK": str(hook),
            "MOCK_MANIFEST_SIDECAR_SUFFIX": suffix,
        }
    )
    completed, packet = _run_runner(
        tmp_path, env, "MR-1", "resume", str(manifest),
        attempt="sidecar", check=False
    )

    assert completed.returncode == 2
    rows = _index_rows(packet)
    assert rows[-1]["sample"] == "warmup"
    assert rows[-1]["disposition"] == "failed-forbidden-assembled"
    scan = packet / "inventory" / "MR-1.resume.sidecar.warmup.staging-scan.txt"
    assert expected in scan.read_text()


def test_slurm_runner_indexes_output_staging_and_rank_map_timeout(tmp_path):
    env, _, _ = _mock_environment(tmp_path)
    env["MOCK_OUTPUT_ASSEMBLED"] = "1"
    completed, packet = _run_runner(
        tmp_path, env, "ED-1", "generate", attempt="assembled", check=False
    )
    assert completed.returncode == 2
    assert _index_rows(packet)[-1]["disposition"] == "failed-forbidden-assembled"

    env.pop("MOCK_OUTPUT_ASSEMBLED")
    env["MOCK_TIMEOUT_RANK_MAP"] = "1"
    completed, packet = _run_runner(
        tmp_path, env, "ED-1", "generate", attempt="timeout", check=False
    )
    assert completed.returncode == 124
    assert _index_rows(packet)[-1]["disposition"] == "timeout-rank-map"


@pytest.mark.parametrize(
    ("corrupt_after", "disposition"),
    (
        (0, "failed-invalid-rank-map-timing"),
        (1, "failed-invalid-launch-timing"),
    ),
)
def test_slurm_runner_rejects_negative_timer_artifacts(
    tmp_path, corrupt_after, disposition
):
    env, _, _ = _mock_environment(tmp_path)
    count = tmp_path / "timer-count"
    count.write_text("0\n")
    corrupt_python = tmp_path / "corrupt-python"
    _write_executable(
        corrupt_python,
        f"""#!/usr/bin/env bash
count="$(cat "{count}")"
printf '%s\n' "$((count + 1))" > "{count}"
if [[ "$count" -ge {corrupt_after} ]]; then
  printf 'monotonic_elapsed_ns\\t-1\\n' > "$2"
  exit 0
fi
exec /usr/bin/python3 "$@"
""",
    )
    env["PYTHON"] = str(corrupt_python)
    completed, packet = _run_runner(
        tmp_path, env, "ED-1", "generate", attempt="negative", check=False
    )

    assert completed.returncode == 2
    assert _index_rows(packet)[-1]["disposition"] == disposition


@pytest.mark.parametrize(
    ("hook_phase", "disposition"),
    (
        ("before", "incomplete-accounting-before"),
        ("after", "incomplete-accounting-after"),
    ),
)
def test_slurm_runner_indexes_accounting_hook_failures(
    tmp_path, hook_phase, disposition
):
    env, manifest, hook = _mock_environment(tmp_path)
    env.update(
        {
            "FS_ACCOUNTING_HOOK": str(hook),
            "MOCK_HOOK_FAIL_PHASE": hook_phase,
        }
    )
    completed, packet = _run_runner(
        tmp_path, env, "MR-1", "resume", str(manifest),
        attempt=hook_phase, check=False
    )

    assert completed.returncode == 9
    rows = _index_rows(packet)
    assert rows[-1]["sample"] == "sample-1"
    assert rows[-1]["disposition"] == disposition
