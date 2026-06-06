#!/bin/bash
#SBATCH --account=AST207
#SBATCH --partition=batch
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --job-name=pic-q011-pressure-gate-validate
#SBATCH --output=/lustre/orion/ast207/proj-shared/dfielding/PIC/logs/slurm/%x.%j.log

set -euo pipefail

export PATH=/opt/cray/pe/python/3.11.7/bin:/usr/bin:/bin
export PYTHONDONTWRITEBYTECODE=1
unset BASH_ENV ENV GIT_ALTERNATE_OBJECT_DIRECTORIES GIT_CONFIG_COUNT
unset GIT_CONFIG_GLOBAL GIT_CONFIG_KEY_0 GIT_CONFIG_NOSYSTEM GIT_CONFIG_SYSTEM
unset GIT_CONFIG_VALUE_0 GIT_DIR GIT_INDEX_FILE GIT_OBJECT_DIRECTORY GIT_WORK_TREE
unset PYTHONHOME PYTHONINSPECT PYTHONPATH PYTHONSTARTUP PYTHONUSERBASE

REPO_ROOT=/autofs/nccs-svm1_home2/dfielding/athenak-pic
PYTHON=/opt/cray/pe/python/3.11.7/bin/python3
TRUSTED_SYSTEM_SITE_PACKAGES=/opt/cray/pe/python/3.11.7/lib/python3.11/site-packages
PUBLICATION_USER_SITE_PACKAGES=/autofs/nccs-svm1_home2/dfielding/.local/lib/python3.11/site-packages
PUBLICATION_RUNTIME_VIEW_RECORD_COUNT=20677
PUBLICATION_RUNTIME_VIEW_SHA256=4ebe4e9218aba12bdb551421b9e3ad1e8ef66fb790c99f744dc7795c7314caa9
TRUSTED_CHECKOUT_TEST_MODULE=tst.publication.test_publish_q011_section54_pressure_pilot_bundle
EXPECTED_GIT_COMMIT="${1:?usage: sbatch $0 FULL_GIT_COMMIT}"
[[ "$EXPECTED_GIT_COMMIT" =~ ^[0-9a-f]{40}$ ]]
PYTHON_ENV=(
  /usr/bin/env -i
  HOME=/
  LANG=C
  LC_ALL=C
  PATH=/usr/bin:/bin
  TMPDIR=/tmp
)
PUBLICATION_TEST_PYTHON_ENV=(
  /usr/bin/env -i
  HOME=/
  LANG=C
  LC_ALL=C
  TMPDIR=/tmp
  PYTHONDONTWRITEBYTECODE=1
  PYTHONNOUSERSITE=1
)

set_authenticated_python_roots() {
  AUTHENTICATED_SOURCE_ROOT="${1:?usage: set_authenticated_python_roots SOURCE_ROOT}"
  AUTHENTICATED_CONTROL_PLANE_ROOT="$AUTHENTICATED_SOURCE_ROOT/tst/publication/frontier_control_plane"
  test -d "$AUTHENTICATED_SOURCE_ROOT"
  test -d "$AUTHENTICATED_CONTROL_PLANE_ROOT"
  test -d "$TRUSTED_SYSTEM_SITE_PACKAGES"
}

run_authenticated_stdin_python() {
  "${PYTHON_ENV[@]}" "$PYTHON" -I -B -S -c \
    'import sys; source_root, control_plane_root, site_packages, *args = sys.argv[1:]; sys.path[:0] = [source_root, control_plane_root, site_packages]; sys.argv = ["-", *args]; source = sys.stdin.buffer.read(); exec(compile(source, "<stdin>", "exec"), {"__name__": "__main__"})' \
    "$AUTHENTICATED_SOURCE_ROOT" "$AUTHENTICATED_CONTROL_PLANE_ROOT" \
    "$TRUSTED_SYSTEM_SITE_PACKAGES" "$@"
}

run_authenticated_module() {
  local module="${1:?usage: run_authenticated_module MODULE [ARG ...]}"
  shift
  "${PYTHON_ENV[@]}" "$PYTHON" -I -B -S -c \
    'import runpy, sys; source_root, control_plane_root, site_packages, module, *args = sys.argv[1:]; sys.path[:0] = [source_root, control_plane_root, site_packages]; sys.argv = [module, *args]; runpy.run_module(module, run_name="__main__", alter_sys=True)' \
    "$AUTHENTICATED_SOURCE_ROOT" "$AUTHENTICATED_CONTROL_PLANE_ROOT" \
    "$TRUSTED_SYSTEM_SITE_PACKAGES" "$module" "$@"
}

set_publication_test_python_roots() {
  PUBLICATION_TEST_SOURCE_ROOT="${1:?usage: set_publication_test_python_roots SOURCE_ROOT}"
  PUBLICATION_TEST_CONTROL_PLANE_ROOT="$PUBLICATION_TEST_SOURCE_ROOT/tst/publication/frontier_control_plane"
  test -d "$PUBLICATION_TEST_SOURCE_ROOT"
  test -d "$PUBLICATION_TEST_CONTROL_PLANE_ROOT"
  test -x "$PUBLICATION_TEST_PYTHON"
  test -d "$PUBLICATION_TEST_USER_RUNTIME_ROOT"
  test -d "$PUBLICATION_TEST_SYSTEM_RUNTIME_ROOT"
}

manage_publication_test_dependencies() {
  local mode="${1:?usage: manage_publication_test_dependencies materialize|verify}"
  "${PYTHON_ENV[@]}" "$PYTHON" -I -B -S - \
    "$mode" "$PUBLICATION_USER_SITE_PACKAGES" "$TRUSTED_SYSTEM_SITE_PACKAGES" \
    "$PUBLICATION_TEST_SITE_PACKAGES_ROOT" "$PYTHON" \
    "$PUBLICATION_RUNTIME_VIEW_RECORD_COUNT" "$PUBLICATION_RUNTIME_VIEW_SHA256" <<'PY'
import hashlib
import json
import os
from pathlib import Path
import stat
import sys

mode, user_site, system_site, destination, base_python, count, digest = sys.argv[1:]
source_roots = [
    ("user", Path(user_site).resolve(strict=True)),
    ("system", Path(system_site).resolve(strict=True)),
]
site_packages_root = Path(destination).resolve(strict=True)
runtime_root = site_packages_root / "runtime"
venv_root = site_packages_root.parents[2]
pth_path = site_packages_root / "publication-runtime.pth"
pth_payload = (
    'import os,site,sys; p=os.path.join(sys.prefix,"lib/python3.11/site-packages/runtime"); '
    'site.addsitedir(os.path.join(p,"user")); site.addsitedir(os.path.join(p,"system"))\n'
)
expected_count = int(count)
expected_digest = digest


def require_venv_scaffold():
    root_info = venv_root.lstat()
    if not stat.S_ISDIR(root_info.st_mode) or root_info.st_mode & 0o222:
        raise SystemExit("publication-test venv root drifted")
    expected_links = {
        venv_root / "bin/python": "python3",
        venv_root / "bin/python3": base_python,
        venv_root / "bin/python3.11": "python3",
        venv_root / "lib64": "lib",
    }
    observed_links = {}
    for directory, dirnames, filenames in os.walk(venv_root):
        for name in dirnames:
            path = Path(directory) / name
            info = path.lstat()
            if stat.S_ISLNK(info.st_mode):
                observed_links[path] = os.readlink(path)
                continue
            if not stat.S_ISDIR(info.st_mode) or info.st_mode & 0o222:
                raise SystemExit(f"publication-test venv directory drifted: {path}")
        for name in filenames:
            path = Path(directory) / name
            info = path.lstat()
            if stat.S_ISLNK(info.st_mode):
                observed_links[path] = os.readlink(path)
                continue
            if not stat.S_ISREG(info.st_mode) or info.st_mode & 0o222:
                raise SystemExit(f"publication-test venv file drifted: {path}")
    if observed_links != expected_links:
        raise SystemExit("publication-test venv link drifted")


def normalized_pyvenv_config():
    config_path = venv_root / "pyvenv.cfg"
    config = {}
    for line in config_path.read_text(encoding="utf-8").splitlines():
        key, separator, value = line.partition(" = ")
        if separator != " = " or not key or key in config:
            raise SystemExit("publication-test pyvenv.cfg drifted")
        config[key] = value
    expected = {
        "home": str(Path(base_python).parent),
        "include-system-site-packages": "false",
        "version": "3.11.7",
        "executable": str(Path(base_python).resolve(strict=True)),
        "command": f"{base_python} -m venv --without-pip {venv_root}",
    }
    if config != expected:
        raise SystemExit("publication-test pyvenv.cfg drifted")
    config["command"] = f"{base_python} -m venv --without-pip {{VENV_ROOT}}"
    return json.dumps(config, separators=(",", ":"), sort_keys=True).encode("utf-8")


def view_records():
    records = [{"kind": "directory", "path": "."}]
    for directory, dirnames, filenames in os.walk(venv_root):
        for name in dirnames:
            path = Path(directory) / name
            relative = path.relative_to(venv_root).as_posix()
            info = path.lstat()
            if stat.S_ISLNK(info.st_mode):
                records.append(
                    {"kind": "symlink", "path": relative, "target": os.readlink(path)}
                )
            elif stat.S_ISDIR(info.st_mode):
                records.append({"kind": "directory", "path": relative})
            else:
                raise SystemExit(f"publication-test venv namespace drifted: {path}")
        for name in filenames:
            path = Path(directory) / name
            relative = path.relative_to(venv_root).as_posix()
            info = path.lstat()
            if stat.S_ISLNK(info.st_mode):
                records.append(
                    {"kind": "symlink", "path": relative, "target": os.readlink(path)}
                )
            elif stat.S_ISREG(info.st_mode):
                payload = (
                    normalized_pyvenv_config()
                    if path == venv_root / "pyvenv.cfg"
                    else path.read_bytes()
                )
                records.append(
                    {
                        "kind": "regular",
                        "path": relative,
                        "sha256": hashlib.sha256(payload).hexdigest(),
                    }
                )
            else:
                raise SystemExit(f"publication-test venv namespace drifted: {path}")
    return sorted(records, key=lambda record: record["path"])


if mode == "materialize":
    if any(site_packages_root.iterdir()):
        raise SystemExit("publication-test runtime destination is not empty")
    for label, source_root in source_roots:
        destination_root = runtime_root / label
        for directory, dirnames, filenames in os.walk(source_root):
            directory_path = Path(directory)
            relative_directory = directory_path.relative_to(source_root)
            for name in dirnames:
                info = (directory_path / name).lstat()
                if not stat.S_ISDIR(info.st_mode):
                    raise SystemExit(f"publication-test source directory drifted: {name}")
            for name in filenames:
                source = directory_path / name
                relative = relative_directory / name
                before = source.lstat()
                if not stat.S_ISREG(before.st_mode):
                    raise SystemExit(f"publication-test source file drifted: {source}")
                if source.suffix == ".pyc" or "__pycache__" in relative.parts:
                    continue
                payload = source.read_bytes()
                after = source.lstat()
                if (
                    before.st_dev,
                    before.st_ino,
                    before.st_size,
                    before.st_mtime_ns,
                ) != (
                    after.st_dev,
                    after.st_ino,
                    after.st_size,
                    after.st_mtime_ns,
                ):
                    raise SystemExit(f"publication-test source changed while reading: {source}")
                target = destination_root / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(payload)
                target.chmod(0o400)
    pth_path.write_text(pth_payload, encoding="utf-8")
    pth_path.chmod(0o400)
    for directory, dirnames, _filenames in os.walk(runtime_root, topdown=False):
        for name in dirnames:
            (Path(directory) / name).chmod(0o500)
    runtime_root.chmod(0o500)
    site_packages_root.chmod(0o500)
    for directory, dirnames, filenames in os.walk(venv_root, topdown=False):
        for name in filenames:
            path = Path(directory) / name
            if not path.is_symlink():
                path.chmod(0o400)
        for name in dirnames:
            path = Path(directory) / name
            if not path.is_symlink():
                path.chmod(0o500)
    venv_root.chmod(0o500)
elif mode != "verify":
    raise SystemExit(f"unsupported publication-test runtime operation: {mode}")

require_venv_scaffold()
records = view_records()
actual_digest = hashlib.sha256(
    json.dumps(records, separators=(",", ":"), sort_keys=True).encode("utf-8")
).hexdigest()
if len(records) != expected_count or actual_digest != expected_digest:
    raise SystemExit("publication-test runtime view digest drifted")
PY
}

verify_publication_test_child_runtime() {
  "${PUBLICATION_TEST_PYTHON_ENV[@]}" \
    PATH="$PUBLICATION_TEST_VENV_ROOT/bin:/usr/bin:/bin" \
    MPLBACKEND=Agg \
    MPLCONFIGDIR="$PUBLICATION_TEST_MPLCONFIG_ROOT" \
    XDG_CACHE_HOME="$PUBLICATION_TEST_MPLCONFIG_ROOT" \
    "$PUBLICATION_TEST_PYTHON" -I -B -c \
    'import matplotlib, numpy, pathlib, scipy, sys; venv, user_root, system_root, original_user, original_system = map(pathlib.Path, sys.argv[1:]); roots = [pathlib.Path(entry).resolve() for entry in sys.path]; user_root = user_root.resolve(); system_root = system_root.resolve(); venv.resolve() == pathlib.Path(sys.prefix).resolve() or sys.exit("publication-test child venv drifted"); user_root in roots and system_root in roots or sys.exit("publication-test child runtime roots absent"); pathlib.Path(original_user).resolve() not in roots and pathlib.Path(original_system).resolve() not in roots or sys.exit("publication-test child admitted mutable source roots"); expected = [(numpy, user_root, "1.26.4"), (scipy, system_root, "1.10.1"), (matplotlib, user_root, "3.10.7")]; [(pathlib.Path(module.__file__).resolve().is_relative_to(root) and module.__version__ == version) or sys.exit(f"publication-test child module drifted: {module.__name__}") for module, root, version in expected]' \
    "$PUBLICATION_TEST_VENV_ROOT" "$PUBLICATION_TEST_USER_RUNTIME_ROOT" \
    "$PUBLICATION_TEST_SYSTEM_RUNTIME_ROOT" "$PUBLICATION_USER_SITE_PACKAGES" \
    "$TRUSTED_SYSTEM_SITE_PACKAGES"
}

run_publication_test_module() {
  local module="${1:?usage: run_publication_test_module MODULE [ARG ...]}"
  shift
  "${PUBLICATION_TEST_PYTHON_ENV[@]}" \
    PATH="$PUBLICATION_TEST_VENV_ROOT/bin:/usr/bin:/bin" \
    MPLBACKEND=Agg \
    MPLCONFIGDIR="$PUBLICATION_TEST_MPLCONFIG_ROOT" \
    XDG_CACHE_HOME="$PUBLICATION_TEST_MPLCONFIG_ROOT" \
    "$PUBLICATION_TEST_PYTHON" -I -B -S -c \
    'import runpy, sys; source_root, control_plane_root, user_root, system_root, module, *args = sys.argv[1:]; sys.path[:0] = [source_root, control_plane_root, user_root, system_root]; sys.version.split()[0] == "3.11.7" or sys.exit("publication-test Python drifted"); sys.argv = [module, *args]; runpy.run_module(module, run_name="__main__", alter_sys=True)' \
    "$PUBLICATION_TEST_SOURCE_ROOT" "$PUBLICATION_TEST_CONTROL_PLANE_ROOT" \
    "$PUBLICATION_TEST_USER_RUNTIME_ROOT" "$PUBLICATION_TEST_SYSTEM_RUNTIME_ROOT" \
    "$module" "$@"
}

cd "$REPO_ROOT"
SOURCE_STATUS=$(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
    /usr/bin/git --no-replace-objects -c core.fsmonitor=false -c core.hooksPath=/dev/null \
    status --porcelain --untracked-files=all
)
test -z "$SOURCE_STATUS"
SOURCE_COMMIT=$(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
    /usr/bin/git --no-replace-objects -c core.fsmonitor=false -c core.hooksPath=/dev/null rev-parse HEAD
)
test "$SOURCE_COMMIT" = "$EXPECTED_GIT_COMMIT"
test "$(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
    /usr/bin/git --no-replace-objects -c core.fsmonitor=false -c core.hooksPath=/dev/null \
    rev-parse origin/PIC
)" = "$EXPECTED_GIT_COMMIT"
echo "source_commit=$SOURCE_COMMIT"

SNAPSHOT_ROOT=$(/usr/bin/mktemp -d "${TMPDIR:-/tmp}/pic-q011-pressure-gate-validate.XXXXXXXX")
FULL_SUITE_ROOT=
PUBLICATION_TEST_MPLCONFIG_ROOT=$(
  /usr/bin/mktemp -d "${TMPDIR:-/tmp}/pic-q011-pressure-gate-matplotlib.XXXXXXXX"
)
PUBLICATION_TEST_VENV_ROOT=$(
  /usr/bin/mktemp -d "${TMPDIR:-/tmp}/pic-q011-pressure-gate-venv.XXXXXXXX"
)
cleanup() {
  /usr/bin/rm -rf "$SNAPSHOT_ROOT"
  if test -n "$FULL_SUITE_ROOT"; then
    /usr/bin/rm -rf "$FULL_SUITE_ROOT"
  fi
  /usr/bin/rm -rf "$PUBLICATION_TEST_MPLCONFIG_ROOT"
  /usr/bin/find "$PUBLICATION_TEST_VENV_ROOT" -type d -exec /usr/bin/chmod u+w {} + \
    2>/dev/null || true
  /usr/bin/rm -rf "$PUBLICATION_TEST_VENV_ROOT"
}
trap cleanup EXIT
/usr/bin/rmdir "$PUBLICATION_TEST_VENV_ROOT"
"${PYTHON_ENV[@]}" "$PYTHON" -I -B -m venv --without-pip "$PUBLICATION_TEST_VENV_ROOT"
/usr/bin/rm -- \
  "$PUBLICATION_TEST_VENV_ROOT/bin/Activate.ps1" \
  "$PUBLICATION_TEST_VENV_ROOT/bin/activate" \
  "$PUBLICATION_TEST_VENV_ROOT/bin/activate.csh" \
  "$PUBLICATION_TEST_VENV_ROOT/bin/activate.fish"
PUBLICATION_TEST_PYTHON="$PUBLICATION_TEST_VENV_ROOT/bin/python3"
PUBLICATION_TEST_SITE_PACKAGES_ROOT="$PUBLICATION_TEST_VENV_ROOT/lib/python3.11/site-packages"
PUBLICATION_TEST_USER_RUNTIME_ROOT="$PUBLICATION_TEST_SITE_PACKAGES_ROOT/runtime/user"
PUBLICATION_TEST_SYSTEM_RUNTIME_ROOT="$PUBLICATION_TEST_SITE_PACKAGES_ROOT/runtime/system"
test "$(/usr/bin/readlink -f "$PUBLICATION_TEST_PYTHON")" = "$(/usr/bin/readlink -f "$PYTHON")"
/usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
  /usr/bin/git --no-replace-objects -c core.fsmonitor=false -c core.hooksPath=/dev/null \
  archive "$SOURCE_COMMIT" | /usr/bin/tar -xf - -C "$SNAPSHOT_ROOT"
cd "$SNAPSHOT_ROOT"
set_authenticated_python_roots "$SNAPSHOT_ROOT"
manage_publication_test_dependencies materialize
verify_publication_test_child_runtime

publication_python_output=$(/usr/bin/find tst/publication -type f -name '*.py' | /usr/bin/sort)
publication_shell_output=$(/usr/bin/find tst/publication -type f -name '*.sh' | /usr/bin/sort)
publication_json_output=$(/usr/bin/find tst/publication -type f -name '*.json' | /usr/bin/sort)
test -n "$publication_python_output"
test -n "$publication_shell_output"
test -n "$publication_json_output"
mapfile -t publication_python <<< "$publication_python_output"
mapfile -t publication_shell <<< "$publication_shell_output"
mapfile -t publication_json <<< "$publication_json_output"
test "${#publication_python[@]}" -eq 147
test "${#publication_shell[@]}" -eq 17
test "${#publication_json[@]}" -eq 292

run_authenticated_stdin_python "${publication_python[@]}" <<'PY'
import sys
from pathlib import Path

for relative in sys.argv[1:]:
    path = Path(relative)
    compile(path.read_bytes(), str(path), "exec")
PY

/usr/bin/bash -n "${publication_shell[@]}"

run_authenticated_stdin_python "${publication_json[@]}" <<'PY'
import json
import sys
from pathlib import Path

for relative in sys.argv[1:]:
    path = Path(relative)
    with path.open(encoding="utf-8") as stream:
        json.load(stream)
PY

run_authenticated_stdin_python <<'PY'
import json
from pathlib import Path

from tst.publication.frontier_control_plane import generate_prepared_pic_artifact_inventory

path = Path("tst/publication/frontier_control_plane/prepared_pic_artifact_inventory.json")
with path.open(encoding="utf-8") as stream:
    committed = json.load(stream)
generated = generate_prepared_pic_artifact_inventory.prepared_artifact_inventory()
if committed != generated:
    raise SystemExit("Committed prepared PIC artifact inventory is stale")
PY

run_authenticated_stdin_python <<'PY'
import fcntl
import os

from tst.publication import publish_q011_section54_pressure_pilot_bundle as publisher

pic_root, _publication_root = publisher._publication_root(publisher.AUTHORIZED_PIC_ROOT)
acceptance_root = publisher._publication_acceptance_root(pic_root)
anchor = publisher._publication_transaction_anchor(pic_root)
anchor_descriptor = publisher._open_absolute_directory(anchor)
acceptance_descriptor = publisher._open_absolute_directory(acceptance_root)
try:
    publisher._lock_publication_transaction(anchor_descriptor, acceptance_descriptor)
    publisher._require_same_directory(
        anchor, anchor_descriptor, "stable publication transaction anchor"
    )
    publisher._require_same_directory(
        acceptance_root,
        acceptance_descriptor,
        "authorized PIC publication acceptance root",
    )
    for path in (anchor, acceptance_root):
        competitor = os.open(path, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
        try:
            try:
                fcntl.flock(competitor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                pass
            else:
                raise SystemExit(f"publication transaction lock did not serialize: {path}")
        finally:
            os.close(competitor)
finally:
    os.close(acceptance_descriptor)
    os.close(anchor_descriptor)
PY

run_authenticated_module unittest \
  tst.publication.frontier_control_plane.test_q011_pressure_review_packet_verifier \
  tst.publication.test_publish_q011_section54_pressure_selection \
  tst.publication.test_q011_section54_pressure_selection \
  tst.publication.test_q011_section54_qualifying_campaign_execution \
  tst.publication.test_q011_section54_attempt_manifest_materializer \
  tst.publication.test_analyze_q011_section54_campaign \
  tst.publication.test_analyze_q011_section54_numerical_qualification \
  tst.publication.test_pic_qualification_manifest \
  tst.publication.frontier_control_plane.test_control_plane.SnapshotTests.test_q011_helper_source_order_matches_execution_and_analyzer \
  tst.publication.frontier_control_plane.test_control_plane.SnapshotTests.test_planner_retention_rejects_missing_pressure_packet_binding \
  tst.publication.frontier_control_plane.test_control_plane.SnapshotTests.test_planner_retention_rejects_pressure_packet_hash_drift \
  tst.publication.frontier_control_plane.test_control_plane.SnapshotTests.test_planner_retention_rejects_cross_bound_pressure_packet_aggregate \
  tst.publication.frontier_control_plane.test_control_plane.SnapshotTests.test_planner_retention_rejects_pressure_bundle_manifest_hash_drift \
  tst.publication.frontier_control_plane.test_control_plane.SnapshotTests.test_planner_retention_rejects_pressure_aggregate_analysis_hash_drift \
  tst.publication.frontier_control_plane.test_control_plane.SnapshotTests.test_planner_retention_rejects_pressure_aggregate_descriptor_drift

FULL_SUITE_ROOT=$(/usr/bin/mktemp -d "${TMPDIR:-/tmp}/pic-q011-pressure-gate-full-suite.XXXXXXXX")
/usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
  /usr/bin/git --no-replace-objects -c core.fsmonitor=false -c core.hooksPath=/dev/null \
  clone --no-checkout --no-local "$REPO_ROOT" "$FULL_SUITE_ROOT"
cd "$FULL_SUITE_ROOT"
/usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
  /usr/bin/git --no-replace-objects -c core.fsmonitor=false -c core.hooksPath=/dev/null \
  checkout --detach "$SOURCE_COMMIT"
test "$(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
    /usr/bin/git --no-replace-objects -c core.fsmonitor=false -c core.hooksPath=/dev/null rev-parse HEAD
)" = "$SOURCE_COMMIT"
test "$(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
    /usr/bin/git --no-replace-objects -c core.fsmonitor=false -c core.hooksPath=/dev/null \
    rev-parse origin/PIC
)" = "$SOURCE_COMMIT"
SOURCE_STATUS=$(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
    /usr/bin/git --no-replace-objects -c core.fsmonitor=false -c core.hooksPath=/dev/null \
    status --porcelain --untracked-files=all
)
test -z "$SOURCE_STATUS"
TAXONOMY_DRY_RUN_SEED="$FULL_SUITE_ROOT/tst/build/src/bin/pic_bell_pub_taxonomy.taxonomy.bin"
test ! -e "$TAXONOMY_DRY_RUN_SEED"
set_publication_test_python_roots "$FULL_SUITE_ROOT"

modules_output=$(
  /usr/bin/find tst/publication -type f -name 'test_*.py' |
    /usr/bin/sed -e 's#/#.#g' -e 's#\.py$##' |
    /usr/bin/grep -vx "$TRUSTED_CHECKOUT_TEST_MODULE" |
    /usr/bin/sort
)
test -n "$modules_output"
mapfile -t modules <<< "$modules_output"
test "${#modules[@]}" -eq 65

manage_publication_test_dependencies verify
verify_publication_test_child_runtime
run_publication_test_module unittest "${modules[@]}"
manage_publication_test_dependencies verify
verify_publication_test_child_runtime
test "$(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
    /usr/bin/git --no-replace-objects -c core.fsmonitor=false -c core.hooksPath=/dev/null rev-parse HEAD
)" = "$SOURCE_COMMIT"
SOURCE_STATUS=$(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
    /usr/bin/git --no-replace-objects -c core.fsmonitor=false -c core.hooksPath=/dev/null \
    status --porcelain --untracked-files=all
)
test -z "$SOURCE_STATUS"
test ! -e "$TAXONOMY_DRY_RUN_SEED"

cd "$REPO_ROOT"
test "$(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
    /usr/bin/git --no-replace-objects -c core.fsmonitor=false -c core.hooksPath=/dev/null \
    rev-parse HEAD
)" = "$SOURCE_COMMIT"
SOURCE_STATUS=$(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
    /usr/bin/git --no-replace-objects -c core.fsmonitor=false -c core.hooksPath=/dev/null \
    status --porcelain --untracked-files=all
)
test -z "$SOURCE_STATUS"
TRUSTED_TAXONOMY_DRY_RUN_SEED="$REPO_ROOT/tst/build/src/bin/pic_bell_pub_taxonomy.taxonomy.bin"
test ! -e "$TRUSTED_TAXONOMY_DRY_RUN_SEED"
set_publication_test_python_roots "$REPO_ROOT"
manage_publication_test_dependencies verify
verify_publication_test_child_runtime
run_publication_test_module unittest "$TRUSTED_CHECKOUT_TEST_MODULE"
manage_publication_test_dependencies verify
verify_publication_test_child_runtime
test "$(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
    /usr/bin/git --no-replace-objects -c core.fsmonitor=false -c core.hooksPath=/dev/null \
    rev-parse HEAD
)" = "$SOURCE_COMMIT"
SOURCE_STATUS=$(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
    /usr/bin/git --no-replace-objects -c core.fsmonitor=false -c core.hooksPath=/dev/null \
    status --porcelain --untracked-files=all
)
test -z "$SOURCE_STATUS"
test ! -e "$TRUSTED_TAXONOMY_DRY_RUN_SEED"
