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
if [[ "${MOCK_TIMEOUT_RANK_MAP:-0}" == 1 && "$*" == *hostname* ]]; then
  exit 124
fi
exec "$@"
""",
    )
    _write_executable(
        bin_dir / "find",
        """#!/usr/bin/env bash
case "${MOCK_FIND_FAIL_MODE:-}:$*" in
  inventory:*outputs*) exit 9 ;;
  staging:*assembled*) exit 9 ;;
  packet-tree:*packet-*) exit 9 ;;
esac
exec /usr/bin/find "$@"
""",
    )
    _write_executable(
        bin_dir / "srun",
        """#!/usr/bin/env bash
output_dir=
manifest=
ntasks=4
rank_map_dir=
previous=
for argument in "$@"; do
  if [[ "$previous" == -d ]]; then output_dir="$argument"; fi
  if [[ "$previous" == -r ]]; then manifest="$argument"; fi
  if [[ "$previous" == rank-map-task ]]; then rank_map_dir="$argument"; fi
  case "$argument" in --ntasks=*) ntasks="${argument#*=}";; esac
  previous="$argument"
done
printf '%s\n' "${ATHENAK_RESTART_MANIFEST_TIMING:-unset}" >> "$MOCK_SRUN_LOG"
if [[ -n "$rank_map_dir" && "${MOCK_NO_RANK_MAP:-0}" != 1 ]]; then
  mkdir -p "$rank_map_dir"
  hosts=()
  if [[ -n "${MOCK_RANK_MAP_HOSTS:-}" ]]; then
    IFS=, read -r -a hosts <<< "$MOCK_RANK_MAP_HOSTS"
  elif [[ -n "${SLURM_HOSTFILE:-}" ]]; then
    while IFS= read -r host || [[ -n "$host" ]]; do hosts+=("$host"); done \
      < "$SLURM_HOSTFILE"
  else
    for ((rank = 0; rank < ntasks; rank++)); do
      if [[ "$rank" -lt "$((ntasks / 2))" ]]; then
        hosts+=(host-a)
      else
        hosts+=(host-b)
      fi
    done
  fi
  for ((rank = 0; rank < ntasks; rank++)); do
    printf '%s\t%s\n' "$rank" "${hosts[$rank]:-missing}" > "$rank_map_dir/rank-$rank.tsv"
  done
  if [[ -n "${MOCK_RANK_MAP_ALIAS_MODE:-}" ]]; then
    target="${MOCK_RANK_MAP_ALIAS_TARGET:?}"
    cp "$rank_map_dir/rank-0.tsv" "$target"
    rm "$rank_map_dir/rank-0.tsv"
    if [[ "$MOCK_RANK_MAP_ALIAS_MODE" == symlink ]]; then
      ln -s "$target" "$rank_map_dir/rank-0.tsv"
    else
      ln "$target" "$rank_map_dir/rank-0.tsv"
    fi
  fi
  if [[ -n "${MOCK_RANK_MAP_DIR_ALIAS_TARGET:-}" ]]; then
    target="$MOCK_RANK_MAP_DIR_ALIAS_TARGET"
    mkdir -p "$target"
    cp "$rank_map_dir"/rank-*.tsv "$target/"
    rm -rf "$rank_map_dir"
    ln -s "$target" "$rank_map_dir"
  fi
  if [[ -n "${MOCK_RANK_MAP_DIR_REPLACE_TARGET:-}" ]]; then
    target="$MOCK_RANK_MAP_DIR_REPLACE_TARGET"
    mv "$rank_map_dir" "$target"
    cp -R "$target" "$rank_map_dir"
  fi
  if [[ -n "${MOCK_RANK_MAP_INVENTORY_ALIAS_PATH:-}" ]]; then
    ln -s "${MOCK_RANK_MAP_INVENTORY_ALIAS_TARGET:?}" \
      "$MOCK_RANK_MAP_INVENTORY_ALIAS_PATH"
  fi
fi
if [[ -n "${MOCK_TIMING_ALIAS_PATH:-}" ]]; then
  ln -s "${MOCK_TIMING_ALIAS_TARGET:?}" "$MOCK_TIMING_ALIAS_PATH"
fi
if [[ -n "${MOCK_TIMING_PARENT_ALIAS_PATH:-}" ]]; then
  mkdir -p "${MOCK_TIMING_PARENT_ALIAS_TARGET:?}"
  rm -rf "$MOCK_TIMING_PARENT_ALIAS_PATH"
  ln -s "$MOCK_TIMING_PARENT_ALIAS_TARGET" "$MOCK_TIMING_PARENT_ALIAS_PATH"
fi
if [[ -n "${MOCK_PACKET_ROOT_ALIAS_PATH:-}" ]]; then
  mkdir -p "${MOCK_PACKET_ROOT_ALIAS_TARGET:?}"
  rm -rf "$MOCK_PACKET_ROOT_ALIAS_PATH"
  ln -s "$MOCK_PACKET_ROOT_ALIAS_TARGET" "$MOCK_PACKET_ROOT_ALIAS_PATH"
fi
if [[ -n "$output_dir" ]]; then
  mkdir -p "$output_dir"
  if [[ -n "${MOCK_OUTPUT_DIR_REPLACE_TARGET:-}" ]]; then
    target="$MOCK_OUTPUT_DIR_REPLACE_TARGET"
    mv "$output_dir" "$target"
    cp -R "$target" "$output_dir"
  fi
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
if [[ "${MOCK_HOOK_ALIAS_PHASE:-}" == "$1" ]]; then
  ln -s "${MOCK_HOOK_ALIAS_TARGET:?}" "${MOCK_HOOK_ALIAS_PATH:?}"
fi
if [[ "${MOCK_HOOK_INDEX_ALIAS_PHASE:-}" == "$1" ]]; then
  mv "${MOCK_HOOK_INDEX_ALIAS_PATH:?}" "${MOCK_HOOK_INDEX_ALIAS_TARGET:?}"
  ln -s "$MOCK_HOOK_INDEX_ALIAS_TARGET" "$MOCK_HOOK_INDEX_ALIAS_PATH"
fi
if [[ "${MOCK_HOOK_REPLACE_PACKET_PHASE:-}" == "$1" ]]; then
  mv "${MOCK_HOOK_REPLACE_PACKET_PATH:?}" "${MOCK_HOOK_REPLACE_PACKET_TARGET:?}"
  cp -R "$MOCK_HOOK_REPLACE_PACKET_TARGET" "$MOCK_HOOK_REPLACE_PACKET_PATH"
fi
if [[ "${MOCK_HOOK_REPLACE_INDEX_PHASE:-}" == "$1" ]]; then
  mv "${MOCK_HOOK_REPLACE_INDEX_PATH:?}" "${MOCK_HOOK_REPLACE_INDEX_TARGET:?}"
  cp "$MOCK_HOOK_REPLACE_INDEX_TARGET" "$MOCK_HOOK_REPLACE_INDEX_PATH"
fi
if [[ "${MOCK_HOOK_REPLACE_LOGS_PHASE:-}" == "$1" ]]; then
  mv "${MOCK_HOOK_REPLACE_LOGS_PATH:?}" "${MOCK_HOOK_REPLACE_LOGS_TARGET:?}"
  cp -R "$MOCK_HOOK_REPLACE_LOGS_TARGET" "$MOCK_HOOK_REPLACE_LOGS_PATH"
fi
if [[ "${MOCK_HOOK_REPLACE_INVENTORY_PHASE:-}" == "$1" ]]; then
  mv "${MOCK_HOOK_REPLACE_INVENTORY_PATH:?}" \
    "${MOCK_HOOK_REPLACE_INVENTORY_TARGET:?}"
  cp -R "$MOCK_HOOK_REPLACE_INVENTORY_TARGET" \
    "$MOCK_HOOK_REPLACE_INVENTORY_PATH"
fi
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


def _corrupt_timer_python(tmp_path: Path, payload: str, corrupt_after: int):
    count = tmp_path / "timer-count"
    count.write_text("0\n")
    corrupt_python = tmp_path / "corrupt-python"
    _write_executable(
        corrupt_python,
        f"""#!/usr/bin/env bash
count="$(cat "{count}")"
printf '%s\n' "$((count + 1))" > "{count}"
if [[ "$count" -ge {corrupt_after} ]]; then
  printf '%b' '{payload}' > "$2"
  exit 0
fi
exec /usr/bin/python3 "$@"
""",
    )
    return corrupt_python


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


@pytest.mark.parametrize(
    "hosts",
    (
        "host-a\n\n\n",
        "host-a\n \n \n",
        "host-a\nhost-b\n",
        "host-a\nhost-a\nhost-a\n",
        "host-a\nhost-b\nhost-c\n",
        " host-a\nhost-b\nhost-b\n",
        "host-a \nhost-b\nhost-b\n",
        "host a\nhost-b\nhost-b\n",
        "host-a\r\nhost-b\r\nhost-b\r\n",
    ),
)
def test_slurm_runner_rejects_noncanonical_mr2_hosts(tmp_path, hosts):
    env, _, _ = _mock_environment(tmp_path)
    hostfile = tmp_path / "mr2.hostfile"
    hostfile.write_text(hosts)
    env["MR2_HOSTFILE"] = str(hostfile)
    completed, packet = _run_runner(
        tmp_path, env, "MR-2", "generate", attempt="blank-host", check=False
    )

    assert completed.returncode == 2
    assert "MR2_HOSTFILE" in completed.stderr
    assert _index_rows(packet) == []


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
    ("env_update", "disposition"),
    (
        ({"MOCK_NO_RANK_MAP": "1"}, "failed-invalid-rank-map"),
        (
            {"MOCK_RANK_MAP_HOSTS": "host-a,host-b,host-a,host-b"},
            "failed-invalid-rank-map",
        ),
    ),
)
def test_slurm_runner_rejects_missing_or_incorrect_rank_map(
    tmp_path, env_update, disposition
):
    env, _, _ = _mock_environment(tmp_path)
    env.update(env_update)
    completed, packet = _run_runner(
        tmp_path, env, "ED-1", "generate", attempt="rank-map", check=False
    )

    assert completed.returncode == 2
    assert _index_rows(packet)[-1]["disposition"] == disposition


@pytest.mark.parametrize(
    ("fail_mode", "disposition"),
    (
        ("inventory", "failed-inventory-scan"),
        ("staging", "failed-staging-scan"),
    ),
)
def test_slurm_runner_indexes_scan_failures(tmp_path, fail_mode, disposition):
    env, _, _ = _mock_environment(tmp_path)
    env["MOCK_FIND_FAIL_MODE"] = fail_mode
    completed, packet = _run_runner(
        tmp_path, env, "ED-1", "generate", attempt="scan", check=False
    )

    assert completed.returncode == 2
    assert _index_rows(packet)[-1]["disposition"] == disposition


def test_slurm_runner_rejects_preseeded_malformed_packet_index(tmp_path):
    env, _, _ = _mock_environment(tmp_path)
    packet = tmp_path / "packet-preseed"
    packet.mkdir()
    (packet / "packet-index.tsv").write_text("wrong\theader\n")

    completed, _ = _run_runner(
        tmp_path, env, "ED-1", "generate", attempt="preseed", check=False
    )

    assert completed.returncode != 0
    assert "invalid packet index" in completed.stderr


def test_slurm_runner_rejects_packet_index_symlink_before_scheduler_work(tmp_path):
    env, _, _ = _mock_environment(tmp_path)
    packet = tmp_path / "packet-symlink"
    external_index = tmp_path / "external-index.tsv"
    packet.mkdir()
    external_index.write_text("external\n")
    (packet / "packet-index.tsv").symlink_to(external_index)

    completed, _ = _run_runner(
        tmp_path, env, "ED-1", "generate", attempt="symlink", check=False
    )

    assert completed.returncode == 2
    assert "must not contain symlinks" in completed.stderr
    assert external_index.read_text() == "external\n"


def test_slurm_runner_rejects_preexisting_deck_hard_link_before_copy(tmp_path):
    env, _, _ = _mock_environment(tmp_path)
    packet = tmp_path / "packet-hardlink"
    external_deck = tmp_path / "external-deck"
    (packet / "decks").mkdir(parents=True)
    external_deck.write_text("external\n")
    os.link(external_deck, packet / "decks" / "io_node_sharding.athinput")

    completed, _ = _run_runner(
        tmp_path, env, "ED-1", "generate", attempt="hardlink", check=False
    )

    assert completed.returncode == 2
    assert "must not contain hard-linked files" in completed.stderr
    assert external_deck.read_text() == "external\n"
    assert not (tmp_path / "srun.log").exists()


def test_slurm_runner_rejects_ambiguous_packet_index_atomic_temporary(tmp_path):
    env, _, _ = _mock_environment(tmp_path)
    packet = tmp_path / "packet-append-temporary"
    packet.mkdir()
    (packet / ".packet-index.tsv.tmp.operator").write_text("ambiguous\n")

    completed, _ = _run_runner(
        tmp_path, env, "ED-1", "generate", attempt="append-temporary", check=False
    )

    assert completed.returncode == 2
    assert "reserved helper-temporary paths" in completed.stderr
    assert not (tmp_path / "srun.log").exists()


def test_slurm_runner_rejects_preexisting_control_character_path(tmp_path):
    env, _, _ = _mock_environment(tmp_path)
    packet = tmp_path / "packet-control-path"
    packet.mkdir()
    (packet / "bad\tartifact").write_text("bad\n")

    completed, _ = _run_runner(
        tmp_path, env, "ED-1", "generate", attempt="control-path", check=False
    )

    assert completed.returncode == 2
    assert "artifact paths must not contain control characters" in completed.stderr
    assert not (tmp_path / "srun.log").exists()


def test_slurm_runner_rejects_packet_tree_scan_failure_before_scheduler_work(tmp_path):
    env, _, _ = _mock_environment(tmp_path)
    env["MOCK_FIND_FAIL_MODE"] = "packet-tree"

    completed, _ = _run_runner(
        tmp_path, env, "ED-1", "generate", attempt="scan-failure", check=False
    )

    assert completed.returncode == 2
    assert "packet symlink scan failed" in completed.stderr
    assert not (tmp_path / "srun.log").exists()


@pytest.mark.parametrize("alias_mode", ("symlink", "hardlink"))
def test_slurm_runner_rejects_aliased_rank_map_evidence(tmp_path, alias_mode):
    env, _, _ = _mock_environment(tmp_path)
    env.update(
        {
            "MOCK_RANK_MAP_ALIAS_MODE": alias_mode,
            "MOCK_RANK_MAP_ALIAS_TARGET": str(tmp_path / f"rank-map-{alias_mode}.tsv"),
        }
    )

    completed, packet = _run_runner(
        tmp_path, env, "ED-1", "generate", attempt=alias_mode, check=False
    )

    assert completed.returncode == 2
    assert _index_rows(packet) == []


def test_slurm_runner_rejects_scheduler_created_rank_map_directory_symlink(tmp_path):
    env, _, _ = _mock_environment(tmp_path)
    env["MOCK_RANK_MAP_DIR_ALIAS_TARGET"] = str(tmp_path / "external-rank-map")

    completed, packet = _run_runner(
        tmp_path, env, "ED-1", "generate", attempt="dir-symlink", check=False
    )

    assert completed.returncode == 2
    assert "must not contain symlinks" in completed.stderr
    assert _index_rows(packet) == []


def test_slurm_runner_rejects_scheduler_created_rank_map_inventory_symlink(tmp_path):
    env, _, _ = _mock_environment(tmp_path)
    packet = tmp_path / "packet-inventory-symlink"
    inventory = packet / "inventory" / "ED-1.generate.inventory-symlink.rank-map.tsv"
    external_target = tmp_path / "outside-rank-map-inventory.tsv"
    env.update(
        {
            "MOCK_RANK_MAP_INVENTORY_ALIAS_PATH": str(inventory),
            "MOCK_RANK_MAP_INVENTORY_ALIAS_TARGET": str(external_target),
        }
    )

    completed, packet = _run_runner(
        tmp_path, env, "ED-1", "generate", attempt="inventory-symlink", check=False
    )

    assert completed.returncode == 2
    assert "must not contain symlinks" in completed.stderr
    assert _index_rows(packet) == []
    assert not external_target.exists()


def test_slurm_runner_rejects_scheduler_created_timing_symlink_without_write_through(
    tmp_path,
):
    env, _, _ = _mock_environment(tmp_path)
    packet = tmp_path / "packet-timing-symlink"
    timing = packet / "logs" / "ED-1.generate.timing-symlink.rank-map.launch-timing.tsv"
    external_target = tmp_path / "outside-timing.tsv"
    env.update(
        {
            "MOCK_TIMING_ALIAS_PATH": str(timing),
            "MOCK_TIMING_ALIAS_TARGET": str(external_target),
        }
    )

    completed, packet = _run_runner(
        tmp_path, env, "ED-1", "generate", attempt="timing-symlink", check=False
    )

    assert completed.returncode == 2
    assert "must not contain symlinks" in completed.stderr
    assert _index_rows(packet) == []
    assert not external_target.exists()


def test_slurm_runner_rejects_scheduler_created_timing_parent_symlink(tmp_path):
    env, _, _ = _mock_environment(tmp_path)
    packet = tmp_path / "packet-timing-parent-symlink"
    external_logs = tmp_path / "outside-logs"
    env.update(
        {
            "MOCK_TIMING_PARENT_ALIAS_PATH": str(packet / "logs"),
            "MOCK_TIMING_PARENT_ALIAS_TARGET": str(external_logs),
        }
    )

    completed, packet = _run_runner(
        tmp_path, env, "ED-1", "generate", attempt="timing-parent-symlink",
        check=False
    )

    assert completed.returncode == 2
    assert "packet logs directory must remain" in completed.stderr
    assert _index_rows(packet) == []
    assert list(external_logs.iterdir()) == []


def test_slurm_runner_rejects_scheduler_created_packet_root_symlink(tmp_path):
    env, _, _ = _mock_environment(tmp_path)
    packet = tmp_path / "packet-root-symlink"
    external_packet = tmp_path / "outside-packet"
    env.update(
        {
            "MOCK_PACKET_ROOT_ALIAS_PATH": str(packet),
            "MOCK_PACKET_ROOT_ALIAS_TARGET": str(external_packet),
        }
    )

    completed, packet = _run_runner(
        tmp_path, env, "ED-1", "generate", attempt="root-symlink", check=False
    )

    assert completed.returncode == 2
    assert "must remain a regular non-symlink directory" in completed.stderr
    assert list(external_packet.iterdir()) == []


def test_slurm_runner_rejects_accounting_hook_alias_before_measured_launch(tmp_path):
    env, manifest, hook = _mock_environment(tmp_path)
    packet = tmp_path / "packet-hook-alias"
    alias_path = packet / "inventory" / "hook-created-alias"
    external_target = tmp_path / "outside-hook-target"
    env.update(
        {
            "FS_ACCOUNTING_HOOK": str(hook),
            "MOCK_HOOK_ALIAS_PHASE": "before",
            "MOCK_HOOK_ALIAS_PATH": str(alias_path),
            "MOCK_HOOK_ALIAS_TARGET": str(external_target),
        }
    )

    completed, packet = _run_runner(
        tmp_path, env, "MR-1", "resume", str(manifest),
        attempt="hook-alias", check=False
    )

    assert completed.returncode == 2
    assert "must not contain symlinks" in completed.stderr
    assert [row["sample"] for row in _index_rows(packet)] == ["rank-map", "warmup"]
    assert not external_target.exists()


def test_slurm_runner_rejects_failing_accounting_hook_index_alias_before_append(tmp_path):
    env, manifest, hook = _mock_environment(tmp_path)
    packet = tmp_path / "packet-hook-index-alias"
    external_index = tmp_path / "outside-packet-index.tsv"
    env.update(
        {
            "FS_ACCOUNTING_HOOK": str(hook),
            "MOCK_HOOK_INDEX_ALIAS_PHASE": "before",
            "MOCK_HOOK_INDEX_ALIAS_PATH": str(packet / "packet-index.tsv"),
            "MOCK_HOOK_INDEX_ALIAS_TARGET": str(external_index),
            "MOCK_HOOK_FAIL_PHASE": "before",
        }
    )

    completed, packet = _run_runner(
        tmp_path, env, "MR-1", "resume", str(manifest),
        attempt="hook-index-alias", check=False
    )

    assert completed.returncode == 2
    assert "must not contain symlinks" in completed.stderr
    assert (packet / "packet-index.tsv").is_symlink()
    assert [row["sample"] for row in _index_rows(packet)] == ["rank-map", "warmup"]


def test_slurm_runner_rejects_clean_packet_root_replacement_by_hook(tmp_path):
    env, manifest, hook = _mock_environment(tmp_path)
    packet = tmp_path / "packet-replaced-root"
    outside_original = tmp_path / "outside-original-packet"
    env.update(
        {
            "FS_ACCOUNTING_HOOK": str(hook),
            "MOCK_HOOK_REPLACE_PACKET_PHASE": "before",
            "MOCK_HOOK_REPLACE_PACKET_PATH": str(packet),
            "MOCK_HOOK_REPLACE_PACKET_TARGET": str(outside_original),
        }
    )

    completed, packet = _run_runner(
        tmp_path, env, "MR-1", "resume", str(manifest),
        attempt="replaced-root", check=False
    )

    assert completed.returncode == 2
    assert "identity changed" in completed.stderr
    assert outside_original.is_dir()
    assert [row["sample"] for row in _index_rows(outside_original)] == [
        "rank-map",
        "warmup",
    ]


def test_slurm_runner_rejects_clean_packet_index_replacement_by_hook(tmp_path):
    env, manifest, hook = _mock_environment(tmp_path)
    packet = tmp_path / "packet-replaced-index"
    outside_original = tmp_path / "outside-original-index.tsv"
    env.update(
        {
            "FS_ACCOUNTING_HOOK": str(hook),
            "MOCK_HOOK_REPLACE_INDEX_PHASE": "before",
            "MOCK_HOOK_REPLACE_INDEX_PATH": str(packet / "packet-index.tsv"),
            "MOCK_HOOK_REPLACE_INDEX_TARGET": str(outside_original),
        }
    )

    completed, packet = _run_runner(
        tmp_path, env, "MR-1", "resume", str(manifest),
        attempt="replaced-index", check=False
    )

    assert completed.returncode == 2
    assert "packet index identity changed" in completed.stderr
    lines = outside_original.read_text().splitlines()
    keys = lines[0].split("\t")
    rows = [dict(zip(keys, line.split("\t"))) for line in lines[1:]]
    assert [row["sample"] for row in rows] == [
        "rank-map",
        "warmup",
    ]


@pytest.mark.parametrize("directory", ("logs", "inventory"))
def test_slurm_runner_rejects_clean_fixed_directory_replacement_by_hook(
    tmp_path, directory
):
    env, manifest, hook = _mock_environment(tmp_path)
    packet = tmp_path / f"packet-replaced-{directory}"
    outside_original = tmp_path / f"outside-original-{directory}"
    upper = directory.upper()
    env.update(
        {
            "FS_ACCOUNTING_HOOK": str(hook),
            f"MOCK_HOOK_REPLACE_{upper}_PHASE": "before",
            f"MOCK_HOOK_REPLACE_{upper}_PATH": str(packet / directory),
            f"MOCK_HOOK_REPLACE_{upper}_TARGET": str(outside_original),
        }
    )

    completed, _ = _run_runner(
        tmp_path, env, "MR-1", "resume", str(manifest),
        attempt=f"replaced-{directory}", check=False
    )

    assert completed.returncode == 2
    assert f"packet {directory} directory identity changed" in completed.stderr
    assert outside_original.is_dir()


def test_slurm_runner_rejects_clean_rank_map_directory_replacement(tmp_path):
    env, _, _ = _mock_environment(tmp_path)
    outside_original = tmp_path / "outside-original-rank-map"
    env["MOCK_RANK_MAP_DIR_REPLACE_TARGET"] = str(outside_original)

    completed, packet = _run_runner(
        tmp_path, env, "ED-1", "generate", attempt="replaced-rank-map", check=False
    )

    assert completed.returncode == 2
    assert "rank-map directory identity changed" in completed.stderr
    assert outside_original.is_dir()
    assert _index_rows(packet) == []


def test_slurm_runner_rejects_clean_output_directory_replacement(tmp_path):
    env, _, _ = _mock_environment(tmp_path)
    outside_original = tmp_path / "outside-original-output"
    env["MOCK_OUTPUT_DIR_REPLACE_TARGET"] = str(outside_original)

    completed, packet = _run_runner(
        tmp_path, env, "ED-1", "generate", attempt="replaced-output", check=False
    )

    assert completed.returncode == 2
    assert "launch output directory identity changed" in completed.stderr
    assert outside_original.is_dir()
    assert [row["sample"] for row in _index_rows(packet)] == ["rank-map"]


def test_slurm_runner_rejects_control_characters_in_index_fields(tmp_path):
    env, _, _ = _mock_environment(tmp_path)
    env["FS_ACCOUNTING_HOOK"] = "hook\tbreak"
    completed, packet = _run_runner(
        tmp_path, env, "ED-1", "generate", attempt="hook-tab", check=False
    )

    assert completed.returncode == 2
    assert "hooks must contain no control characters" in completed.stderr
    assert not (packet / "packet-index.tsv").exists()


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
    env["PYTHON"] = str(
        _corrupt_timer_python(tmp_path, r"monotonic_elapsed_ns\t-1\n", corrupt_after)
    )
    completed, packet = _run_runner(
        tmp_path, env, "ED-1", "generate", attempt="negative", check=False
    )

    assert completed.returncode == 2
    assert _index_rows(packet)[-1]["disposition"] == disposition


@pytest.mark.parametrize(
    "payload",
    (
        "",
        r"monotonic_elapsed_ns\t7\nextra\tbad\n",
        r"monotonic_elapsed_ns\t7\textra\n",
        r"monotonic_elapsed_ns\tbad\n",
        r"wrong_key\t7\n",
    ),
)
@pytest.mark.parametrize(
    ("corrupt_after", "disposition"),
    (
        (0, "failed-invalid-rank-map-timing"),
        (1, "failed-invalid-launch-timing"),
    ),
)
def test_slurm_runner_rejects_malformed_timer_artifacts(
    tmp_path, payload, corrupt_after, disposition
):
    env, _, _ = _mock_environment(tmp_path)
    env["PYTHON"] = str(_corrupt_timer_python(tmp_path, payload, corrupt_after))
    completed, packet = _run_runner(
        tmp_path, env, "ED-1", "generate", attempt="malformed", check=False
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
