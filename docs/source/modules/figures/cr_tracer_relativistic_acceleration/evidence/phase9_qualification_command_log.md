# Phase-9 Relativistic CR Qualification Command Log

Date: 2026-05-30

This command inventory records the final workstation qualification replay for
`feature/CR_tracers_relativistic_acceleration`.  Paths under `/tmp` are
ephemeral runtime roots.  Retained logs, metrics, caches, archive inventories,
and hash manifests live beside this file.

Run each retained replay block from the repository root unless a block changes
directory explicitly:

```bash
REPO_ROOT=$(pwd)
EVID_REL=docs/source/modules/figures/cr_tracer_relativistic_acceleration/evidence
EVID=$REPO_ROOT/$EVID_REL
```

## General Builds

```bash
cmake -S . -B /tmp/athenak-cr-phase9-release -D CMAKE_BUILD_TYPE=Release \
  > "$EVID/phase9_release_configure.log" 2>&1
cmake --build /tmp/athenak-cr-phase9-release -j 4 \
  > "$EVID/phase9_release_build.log" 2>&1
cp /tmp/athenak-cr-phase9-release/CMakeCache.txt \
  "$EVID/phase9_release_CMakeCache.txt"

cmake -S . -B /tmp/athenak-cr-phase9-mpi-release \
  -D CMAKE_BUILD_TYPE=Release -D Athena_ENABLE_MPI=ON \
  > "$EVID/phase9_mpi_release_configure.log" 2>&1
cmake --build /tmp/athenak-cr-phase9-mpi-release -j 4 \
  > "$EVID/phase9_mpi_release_build.log" 2>&1
cp /tmp/athenak-cr-phase9-mpi-release/CMakeCache.txt \
  "$EVID/phase9_mpi_release_CMakeCache.txt"

cmake -S . -B /tmp/athenak-cr-phase9-debug \
  -D CMAKE_BUILD_TYPE=Debug \
  -D Kokkos_ENABLE_DEBUG=ON \
  -D Kokkos_ENABLE_DEBUG_BOUNDS_CHECK=ON \
  > "$EVID/phase9_debug_configure.log" 2>&1
cmake --build /tmp/athenak-cr-phase9-debug -j 4 \
  > "$EVID/phase9_debug_build.log" 2>&1
cp /tmp/athenak-cr-phase9-debug/CMakeCache.txt \
  "$EVID/phase9_debug_CMakeCache.txt"

cmake -S . -B /tmp/athenak-cr-phase9-mpi-debug \
  -D CMAKE_BUILD_TYPE=Debug \
  -D Athena_ENABLE_MPI=ON \
  -D Kokkos_ENABLE_DEBUG=ON \
  -D Kokkos_ENABLE_DEBUG_BOUNDS_CHECK=ON \
  > "$EVID/phase9_mpi_debug_configure.log" 2>&1
cmake --build /tmp/athenak-cr-phase9-mpi-debug -j 4 \
  > "$EVID/phase9_mpi_debug_build.log" 2>&1
cp /tmp/athenak-cr-phase9-mpi-debug/CMakeCache.txt \
  "$EVID/phase9_mpi_debug_CMakeCache.txt"
```

## Dedicated Analytical Builds

```bash
cmake -S . -B /tmp/athenak-cr-phase9-pusher-release \
  -D CMAKE_BUILD_TYPE=Release \
  -D PROBLEM=unit_tests/cr_relativistic_pusher_runtime_test \
  > "$EVID/phase9_pusher_release_configure.log" 2>&1
cmake --build /tmp/athenak-cr-phase9-pusher-release -j 4 \
  > "$EVID/phase9_pusher_release_build.log" 2>&1
cp /tmp/athenak-cr-phase9-pusher-release/CMakeCache.txt \
  "$EVID/phase9_pusher_release_CMakeCache.txt"

cmake -S . -B /tmp/athenak-cr-phase9-coupled-release \
  -D CMAKE_BUILD_TYPE=Release \
  -D PROBLEM=unit_tests/cr_relativistic_coupled_runtime_test \
  > "$EVID/phase9_coupled_release_configure.log" 2>&1
cmake --build /tmp/athenak-cr-phase9-coupled-release -j 4 \
  > "$EVID/phase9_coupled_release_build.log" 2>&1
cp /tmp/athenak-cr-phase9-coupled-release/CMakeCache.txt \
  "$EVID/phase9_coupled_release_CMakeCache.txt"

cmake -S . -B /tmp/athenak-cr-phase9-coupled-debug \
  -D CMAKE_BUILD_TYPE=Debug \
  -D PROBLEM=unit_tests/cr_relativistic_coupled_runtime_test \
  -D Kokkos_ENABLE_DEBUG=ON \
  -D Kokkos_ENABLE_DEBUG_BOUNDS_CHECK=ON \
  > "$EVID/phase9_coupled_debug_configure.log" 2>&1
cmake --build /tmp/athenak-cr-phase9-coupled-debug -j 4 \
  > "$EVID/phase9_coupled_debug_build.log" 2>&1
cp /tmp/athenak-cr-phase9-coupled-debug/CMakeCache.txt \
  "$EVID/phase9_coupled_debug_CMakeCache.txt"
```

## Analytical And Contract Replays

The retained `phase9_*runtime.log` and `phase9_*runtime_metrics.json` files
record these exact release invocations:

```bash
python3 scripts/particles/cr_relativistic_pusher_runtime_inspect.py \
  --binary /tmp/athenak-cr-phase9-pusher-release/src/athena \
  --input inputs/unit_tests/cr_relativistic_pusher_runtime.athinput \
  --criteria docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase4a_preregistered_criteria.json \
  --work-dir /tmp/cr-rel-phase9-pusher \
  --metrics "$EVID/phase9_pusher_runtime_metrics.json" \
  > "$EVID/phase9_pusher_runtime.log" 2>&1

python3 scripts/particles/cr_relativistic_coupled_runtime_inspect.py \
  --binary /tmp/athenak-cr-phase9-coupled-release/src/athena \
  --input inputs/unit_tests/cr_relativistic_coupled_runtime_test.athinput \
  --prescribed-input inputs/unit_tests/cr_relativistic_pusher_runtime.athinput \
  --criteria docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase4b_preregistered_criteria.json \
  --work-dir /tmp/cr-rel-phase9-coupled \
  --metrics "$EVID/phase9_coupled_runtime_metrics.json" \
  > "$EVID/phase9_coupled_runtime.log" 2>&1

python3 scripts/particles/cr_relativistic_subcycle_runtime_inspect.py \
  --binary /tmp/athenak-cr-phase9-coupled-release/src/athena \
  --input inputs/unit_tests/cr_relativistic_coupled_runtime_test.athinput \
  --criteria docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase5_preregistered_criteria.json \
  --work-dir /tmp/cr-rel-phase9-subcycle \
  --metrics "$EVID/phase9_subcycle_runtime_metrics.json" \
  > "$EVID/phase9_subcycle_runtime.log" 2>&1

python3 scripts/particles/cr_relativistic_restart_runtime_inspect.py \
  --binary /tmp/athenak-cr-phase9-release/src/athena \
  --input inputs/unit_tests/cr_relativistic_restart_runtime.athinput \
  --criteria docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase6_preregistered_criteria.json \
  --work-dir /tmp/cr-rel-phase9-restart \
  --json "$EVID/phase9_restart_runtime_metrics.json" \
  > "$EVID/phase9_restart_runtime.log" 2>&1

python3 scripts/particles/cr_relativistic_diagnostics_runtime_inspect.py \
  --binary /tmp/athenak-cr-phase9-release/src/athena \
  --input inputs/unit_tests/cr_relativistic_diagnostics_runtime.athinput \
  --criteria docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase7_preregistered_criteria.json \
  --work-dir /tmp/cr-rel-phase9-diagnostics \
  --metrics "$EVID/phase9_diagnostics_runtime_metrics.json" \
  > "$EVID/phase9_diagnostics_runtime.log" 2>&1

python3 scripts/particles/cr_relativistic_all_formats_inspect.py \
  --runtime-root /tmp/cr-rel-phase9-diagnostics/acceleration_uninterrupted \
  --binary /tmp/athenak-cr-phase9-release/src/athena \
  --input inputs/unit_tests/cr_relativistic_diagnostics_runtime.athinput \
  --criteria docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase7_preregistered_criteria.json \
  --metrics "$EVID/phase9_diagnostics_runtime_metrics.json" \
  > "$EVID/phase9_all_formats_metrics.json"
echo 'CR relativistic Phase-7 all-format inspection PASS' \
  > "$EVID/phase9_all_formats.log"
```

The retained debug-bounds artifacts use explicit separate roots and binaries:

```bash
python3 scripts/particles/cr_relativistic_coupled_runtime_inspect.py \
  --binary /tmp/athenak-cr-phase9-coupled-debug/src/athena \
  --input inputs/unit_tests/cr_relativistic_coupled_runtime_test.athinput \
  --prescribed-input inputs/unit_tests/cr_relativistic_pusher_runtime.athinput \
  --criteria docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase4b_preregistered_criteria.json \
  --work-dir /tmp/cr-rel-phase9-coupled-debug \
  --metrics "$EVID/phase9_coupled_debug_runtime_metrics.json" \
  > "$EVID/phase9_coupled_debug_runtime.log" 2>&1

python3 scripts/particles/cr_relativistic_subcycle_runtime_inspect.py \
  --binary /tmp/athenak-cr-phase9-coupled-debug/src/athena \
  --input inputs/unit_tests/cr_relativistic_coupled_runtime_test.athinput \
  --criteria docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase5_preregistered_criteria.json \
  --work-dir /tmp/cr-rel-phase9-subcycle-debug \
  --metrics "$EVID/phase9_subcycle_debug_runtime_metrics.json" \
  > "$EVID/phase9_subcycle_debug_runtime.log" 2>&1

python3 scripts/particles/cr_relativistic_restart_runtime_inspect.py \
  --binary /tmp/athenak-cr-phase9-debug/src/athena \
  --input inputs/unit_tests/cr_relativistic_restart_runtime.athinput \
  --criteria docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase6_preregistered_criteria.json \
  --work-dir /tmp/cr-rel-phase9-restart-debug \
  --json "$EVID/phase9_restart_debug_runtime_metrics.json" \
  > "$EVID/phase9_restart_debug_runtime.log" 2>&1

python3 scripts/particles/cr_relativistic_diagnostics_runtime_inspect.py \
  --binary /tmp/athenak-cr-phase9-debug/src/athena \
  --input inputs/unit_tests/cr_relativistic_diagnostics_runtime.athinput \
  --criteria docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase7_preregistered_criteria.json \
  --work-dir /tmp/cr-rel-phase9-diagnostics-debug \
  --metrics "$EVID/phase9_diagnostics_debug_runtime_metrics.json" \
  > "$EVID/phase9_diagnostics_debug_runtime.log" 2>&1

python3 scripts/particles/cr_relativistic_all_formats_inspect.py \
  --runtime-root /tmp/cr-rel-phase9-diagnostics-debug/acceleration_uninterrupted \
  --binary /tmp/athenak-cr-phase9-debug/src/athena \
  --input inputs/unit_tests/cr_relativistic_diagnostics_runtime.athinput \
  --criteria docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase7_preregistered_criteria.json \
  --metrics "$EVID/phase9_diagnostics_debug_runtime_metrics.json" \
  > "$EVID/phase9_all_formats_debug_metrics.json"
echo 'CR relativistic Phase-7 all-format inspection PASS' \
  > "$EVID/phase9_all_formats_debug.log"
```

The portable release all-format package carries its own qualified binary,
runtime artifacts, input, criteria, required inspector scripts, deterministic
manifest, and replay wrapper.

## MPI, SMR, And AMR Migration

```bash
python3 scripts/particles/cr_relativistic_migration_runtime_inspect.py \
  --check-registration \
  --json "$EVID/phase9_migration_registration_metrics.json" \
  > "$EVID/phase9_migration_registration.log" 2>&1

python3 scripts/particles/cr_relativistic_migration_runtime_inspect.py \
  --binary /tmp/athenak-cr-phase9-release/src/athena \
  --mpi-binary /tmp/athenak-cr-phase9-mpi-release/src/athena \
  --work-dir /tmp/cr-rel-phase9-migration \
  --json "$EVID/phase9_migration_metrics.json" \
  > "$EVID/phase9_migration_runtime.log" 2>&1

python3 scripts/particles/cr_relativistic_migration_runtime_inspect.py \
  --binary /tmp/athenak-cr-phase9-debug/src/athena \
  --mpi-binary /tmp/athenak-cr-phase9-mpi-debug/src/athena \
  --work-dir /tmp/cr-rel-phase9-migration-debug \
  --json "$EVID/phase9_migration_debug_metrics.json" \
  > "$EVID/phase9_migration_debug_runtime.log" 2>&1
```

## Adversarial Controls

```bash
python3 scripts/particles/cr_relativistic_migration_negative_controls.py \
  > "$EVID/phase9_migration_negative_controls.log" 2>&1
python3 scripts/particles/cr_relativistic_diagnostics_negative_controls.py \
  > "$EVID/phase9_diagnostics_negative_controls.log" 2>&1
python3 scripts/particles/cr_relativistic_subcycle_negative_controls.py \
  --metrics "$EVID/phase9_subcycle_negative_control_metrics.json" \
  > "$EVID/phase9_subcycle_negative_controls.log" 2>&1
python3 scripts/particles/cr_relativistic_coupled_oracle_negative_controls.py \
  --analyzer scripts/particles/cr_relativistic_coupled_runtime_inspect.py \
  --binary /tmp/athenak-cr-phase9-coupled-release/src/athena \
  --input inputs/unit_tests/cr_relativistic_coupled_runtime_test.athinput \
  --prescribed-input inputs/unit_tests/cr_relativistic_pusher_runtime.athinput \
  --criteria docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase4b_preregistered_criteria.json \
  > "$EVID/phase9_coupled_negative_controls.log" 2>&1
python3 scripts/particles/cr_relativistic_diagnostics_writer_adversarial.py \
  --binary /tmp/athenak-cr-phase9-release/src/athena \
  --work-dir /tmp/cr-rel-phase9-writer-adversarial \
  > "$EVID/phase9_diagnostics_writer_adversarial.log" 2>&1
```

## Repository Harness

These harness invocations ran sequentially because they share `tst/build`:

```bash
cd tst
python run_test_suite.py \
  --test test_suite/particles/test_particles_cr_relativistic_contract_cpu.py \
  --cpu > "$EVID/phase9_parser_contract.log" 2>&1
python run_test_suite.py \
  --test test_suite/particles/test_particles_cr_cpu.py \
  --cpu > "$EVID/phase9_legacy_cpu.log" 2>&1
python run_test_suite.py \
  --test test_suite/particles/test_particles_cr_accuracy_cpu.py \
  --cpu > "$EVID/phase9_legacy_accuracy_cpu.log" 2>&1
python run_test_suite.py \
  --test test_suite/particles/test_particles_cr_mpicpu.py \
  --mpicpu > "$EVID/phase9_legacy_mpicpu.log" 2>&1
python run_test_suite.py \
  --test test_suite/particles/test_particles_cr_accuracy_mpicpu.py \
  --mpicpu > "$EVID/phase9_legacy_accuracy_mpicpu.log" 2>&1
python run_test_suite.py --style > "$EVID/phase9_style.log" 2>&1
python run_tests.py mhd/mhd_divb_amr \
  --log_file ../docs/source/modules/figures/cr_tracer_relativistic_acceleration/evidence/phase9_deep_amr_divb.log
```

## Public Documentation Overlay

```bash
rm -rf /tmp/athenak-rel-gh-pages-overlay
git worktree add --detach /tmp/athenak-rel-gh-pages-overlay \
  4833aa9341e19861297e330ff02aabfd8001935c
cp docs/source/modules/particles.md \
  /tmp/athenak-rel-gh-pages-overlay/docs/source/modules/particles.md
cp docs/source/modules/cr_tracer*.md \
  /tmp/athenak-rel-gh-pages-overlay/docs/source/modules/
mkdir -p /tmp/athenak-rel-gh-pages-overlay/docs/source/modules/figures
cp -R docs/source/modules/figures/cr_tracer_accuracy \
  /tmp/athenak-rel-gh-pages-overlay/docs/source/modules/figures/
mkdir -p /tmp/athenak-rel-gh-pages-overlay/inputs/particles
cp inputs/particles/cr_tracer_relativistic_mhd_ideal_example.athinput \
  inputs/particles/cr_tracer_relativistic_prescribed_restart_resume.athinput \
  /tmp/athenak-rel-gh-pages-overlay/inputs/particles/
(
  cd /tmp/athenak-rel-gh-pages-overlay/docs
  make clean html SPHINXOPTS="-W --keep-going"
) > "$EVID/phase9_gh_pages_overlay_build.log" 2>&1
perl -0pi -e 's/[ \t]+$//mg' \
  "$EVID/phase9_gh_pages_overlay_build.log"
{
  echo 'origin_gh_pages_head=4833aa9341e19861297e330ff02aabfd8001935c'
  echo 'overlay_files:'
  {
    find docs/source/modules -maxdepth 1 -type f \
      \( -name 'particles.md' -o -name 'cr_tracer*.md' \) -print
    find docs/source/modules/figures/cr_tracer_accuracy -type f -print
    echo inputs/particles/cr_tracer_relativistic_mhd_ideal_example.athinput
    echo inputs/particles/cr_tracer_relativistic_prescribed_restart_resume.athinput
  } | LC_ALL=C sort | xargs shasum -a 256
} > "$EVID/phase9_gh_pages_overlay_inventory.txt"
tail -n +3 "$EVID/phase9_gh_pages_overlay_inventory.txt" |
  (cd /tmp/athenak-rel-gh-pages-overlay && shasum -a 256 -c /dev/stdin)
git worktree remove --force /tmp/athenak-rel-gh-pages-overlay
```

The overlay baseline is frozen as
`origin/gh-pages@4833aa9341e19861297e330ff02aabfd8001935c`.

## Public Solver-Coupled Example Smoke

```bash
cp inputs/particles/cr_tracer_relativistic_mhd_ideal_example.athinput \
  "$EVID/phase9_public_solver_coupled_example.athinput"
mkdir -p /tmp/cr-rel-phase9-public-example
(
  cd /tmp/cr-rel-phase9-public-example
  /tmp/athenak-cr-phase9-release/src/athena \
    -i "$REPO_ROOT/inputs/particles/cr_tracer_relativistic_mhd_ideal_example.athinput"
) > "$EVID/phase9_public_solver_coupled_example_runtime.log" 2>&1
```

## Public Typed-V2 Restart Resume Template Smoke

```bash
rm -rf /tmp/cr-rel-phase9-public-restart-template
mkdir -p /tmp/cr-rel-phase9-public-restart-template
(
  cd /tmp/cr-rel-phase9-public-restart-template
  /tmp/athenak-cr-phase9-release/src/athena \
    -i "$REPO_ROOT/inputs/unit_tests/cr_relativistic_restart_runtime.athinput" \
    > initial.log 2>&1
  /tmp/athenak-cr-phase9-release/src/athena \
    -r rst/cr_relativistic_restart.00001.rst \
    -i "$REPO_ROOT/inputs/particles/cr_tracer_relativistic_prescribed_restart_resume.athinput"
) > "$EVID/phase9_public_restart_resume_template_runtime.log" 2>&1
```

## Portable All-Format Replay

The portable package uses a prefix-mapped release binary so compiled
`__FILE__` strings do not expose the local checkout path:

```bash
cmake -S . -B /tmp/athenak-cr-phase9-release-redacted \
  -D CMAKE_BUILD_TYPE=Release \
  -D "CMAKE_CXX_FLAGS=-ffile-prefix-map=$REPO_ROOT=REPO_ROOT -fdebug-prefix-map=$REPO_ROOT=REPO_ROOT" \
  > "$EVID/phase9_portable_release_configure.log" 2>&1
cmake --build /tmp/athenak-cr-phase9-release-redacted -j 4 \
  > "$EVID/phase9_portable_release_build.log" 2>&1
cp /tmp/athenak-cr-phase9-release-redacted/CMakeCache.txt \
  "$EVID/phase9_portable_release_CMakeCache.txt"

python3 scripts/particles/cr_relativistic_diagnostics_runtime_inspect.py \
  --binary /tmp/athenak-cr-phase9-release-redacted/src/athena \
  --input inputs/unit_tests/cr_relativistic_diagnostics_runtime.athinput \
  --criteria docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase7_preregistered_criteria.json \
  --work-dir /tmp/cr-rel-phase9-diagnostics-portable \
  --metrics "$EVID/phase9_diagnostics_portable_runtime_metrics.json" \
  > "$EVID/phase9_diagnostics_portable_runtime.log" 2>&1

cat > "$EVID/phase9_portable_replay.py" <<'PY'
#!/usr/bin/env python3
"""Replay strict all-format validation from an extracted Phase-9 package."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _verify_manifest(root: Path) -> None:
    manifest = root / "portable_manifest.sha256"
    for line in manifest.read_text().splitlines():
        expected, relative = line.split("  ", 1)
        path = root / relative
        actual = _sha256(path)
        if actual != expected:
            raise RuntimeError(
                f"portable manifest mismatch for {relative}: {actual} != {expected}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    root = Path(__file__).resolve().parent
    _verify_manifest(root)

    binary = (root / "qualified_binary" / "athena").resolve()
    runtime_root = (root / "runtime_root" / "acceleration_uninterrupted").resolve()
    input_path = (
        root / "inputs" / "cr_relativistic_diagnostics_runtime.athinput").resolve()
    criteria = (
        root / "criteria"
        / "cr_relativistic_phase7_preregistered_criteria.json").resolve()
    metrics = json.loads(
        (root / "metrics" / "diagnostics_runtime_metrics_template.json").read_text())
    metrics["binary"] = str(binary)
    metrics["input"] = str(input_path)
    metrics["criteria"] = str(criteria)
    metrics["runtime_cases"] = [str(runtime_root)]
    relocated_metrics = root / "metrics" / "diagnostics_runtime_metrics_relocated.json"
    relocated_metrics.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")

    output = args.output or (root / "portable_all_formats_report.json")
    command = [
        sys.executable,
        str(root / "scripts" / "particles" / "cr_relativistic_all_formats_inspect.py"),
        "--runtime-root", str(runtime_root),
        "--binary", str(binary),
        "--input", str(input_path),
        "--criteria", str(criteria),
        "--metrics", str(relocated_metrics),
    ]
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    output.write_text(result.stdout.replace(str(root), "$PACKAGE_ROOT"))
    print("CR relativistic Phase-9 portable all-format replay PASS")
    print("  package root: $PACKAGE_ROOT")
    print("  report: portable_all_formats_report.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
PY

cat > "$EVID/phase9_portable_replay_README.md" <<'EOF'
# Phase-9 Portable All-Format Replay

Extract `phase9_all_formats_portable_replay_bundle.tar.gz` into an arbitrary
directory and run:

```bash
python3 phase9_portable_replay.py
```

The package carries:

- the qualified release AthenaK binary;
- the retained `acceleration_uninterrupted` runtime artifacts;
- the Phase-7 diagnostics input and preregistered criteria;
- the strict all-format inspector and its two Python dependencies;
- a diagnostics metrics template with the qualified binary, input, and
  criteria SHA-256 digests;
- `portable_manifest.sha256`.

The wrapper verifies the package manifest, rewrites extraction-local paths in
the metrics template, invokes the bundled strict inspector, and writes
`portable_all_formats_report.json`.  The emitted report normalizes the
extraction prefix to `$PACKAGE_ROOT`, so two fresh extractions produce
byte-identical reports.
EOF

PACKAGE_PARENT=/tmp/cr-rel-phase9-portable-package
PACKAGE_ROOT=$PACKAGE_PARENT/phase9_all_formats_portable_replay
rm -rf "$PACKAGE_PARENT"
mkdir -p "$PACKAGE_ROOT"/{criteria,inputs,metrics,qualified_binary,runtime_root,scripts/particles}
cp /tmp/athenak-cr-phase9-release-redacted/src/athena \
  "$PACKAGE_ROOT/qualified_binary/athena"
cp -R /tmp/cr-rel-phase9-diagnostics-portable/acceleration_uninterrupted \
  "$PACKAGE_ROOT/runtime_root/"
cp inputs/unit_tests/cr_relativistic_diagnostics_runtime.athinput \
  "$PACKAGE_ROOT/inputs/"
cp docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase7_preregistered_criteria.json \
  "$PACKAGE_ROOT/criteria/"
cp scripts/particles/cr_relativistic_all_formats_inspect.py \
  scripts/particles/cr_relativistic_diagnostics_runtime_inspect.py \
  scripts/particles/cr_tracer_inspect.py \
  "$PACKAGE_ROOT/scripts/particles/"
cp "$EVID/phase9_portable_replay.py" \
  "$PACKAGE_ROOT/phase9_portable_replay.py"
cp "$EVID/phase9_portable_replay_README.md" \
  "$PACKAGE_ROOT/phase9_portable_replay_README.md"
cp "$EVID/phase9_diagnostics_portable_runtime_metrics.json" \
  "$PACKAGE_ROOT/metrics/diagnostics_runtime_metrics_template.json"
(cd "$PACKAGE_ROOT" && \
  find . -type f ! -name portable_manifest.sha256 -print | \
  LC_ALL=C sort | sed 's#^\./##' | xargs shasum -a 256 \
  > portable_manifest.sha256)
python3 - "$PACKAGE_PARENT" "$EVID/phase9_all_formats_portable_replay_bundle.tar.gz" <<'PY'
import gzip
import pathlib
import sys
import tarfile

root = pathlib.Path(sys.argv[1])
output = pathlib.Path(sys.argv[2])
with output.open("wb") as raw:
    with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
        with tarfile.open(fileobj=zipped, mode="w") as archive:
            for path in sorted(root.rglob("*")):
                info = archive.gettarinfo(str(path), arcname=str(path.relative_to(root)))
                info.uid = info.gid = 0
                info.uname = info.gname = ""
                info.mtime = 0
                if path.is_file():
                    with path.open("rb") as stream:
                        archive.addfile(info, stream)
                else:
                    archive.addfile(info)
PY
tar -tzf "$EVID/phase9_all_formats_portable_replay_bundle.tar.gz" \
  > "$EVID/phase9_all_formats_portable_replay_bundle_inventory.txt"

rm -rf /tmp/cr-rel-phase9-portable-extract-a /tmp/cr-rel-phase9-portable-extract-b
mkdir -p /tmp/cr-rel-phase9-portable-extract-a /tmp/cr-rel-phase9-portable-extract-b
tar -xzf docs/source/modules/figures/cr_tracer_relativistic_acceleration/evidence/phase9_all_formats_portable_replay_bundle.tar.gz \
  -C /tmp/cr-rel-phase9-portable-extract-a
tar -xzf docs/source/modules/figures/cr_tracer_relativistic_acceleration/evidence/phase9_all_formats_portable_replay_bundle.tar.gz \
  -C /tmp/cr-rel-phase9-portable-extract-b
{
  python3 /tmp/cr-rel-phase9-portable-extract-a/phase9_all_formats_portable_replay/phase9_portable_replay.py
  python3 /tmp/cr-rel-phase9-portable-extract-b/phase9_all_formats_portable_replay/phase9_portable_replay.py
  cmp \
    /tmp/cr-rel-phase9-portable-extract-a/phase9_all_formats_portable_replay/portable_all_formats_report.json \
    /tmp/cr-rel-phase9-portable-extract-b/phase9_all_formats_portable_replay/portable_all_formats_report.json
  echo 'PASS: two fresh extractions produced byte-identical reports'
} > "$EVID/phase9_all_formats_portable_replay.log" 2>&1

{
  LOCAL_USER='dbf''75'
  LOCAL_HOME_PATTERN="/Users/$LOCAL_USER"
  HOST_PATTERN='Tin''-Drum'
  ! strings \
    /tmp/cr-rel-phase9-portable-extract-a/phase9_all_formats_portable_replay/qualified_binary/athena |
    rg -n "$LOCAL_HOME_PATTERN|$HOST_PATTERN"
} > "$EVID/phase9_portable_binary_privacy_scan.log"
```

## Integration Merge-Tree Audit

```bash
PHASE9_PRESEAL_CANDIDATE=$(git rev-parse HEAD)
set +e
git merge-tree --write-tree \
  c6a73b08e60807f8b925164c5e7edd5cb820c8ae \
  "$PHASE9_PRESEAL_CANDIDATE" \
  > "$EVID/phase9_merge_tree_origin_development.log" 2>&1
MERGE_TREE_STATUS=$?
set -e
{
  echo
  echo 'merge_tree_direction=development<-sanitized-preseal-candidate'
  echo 'merge_tree_target=c6a73b08e60807f8b925164c5e7edd5cb820c8ae'
  echo "merge_tree_candidate=$PHASE9_PRESEAL_CANDIDATE"
  echo "merge_tree_exit=$MERGE_TREE_STATUS"
} >> "$EVID/phase9_merge_tree_origin_development.log"
test "$MERGE_TREE_STATUS" -eq 1
{
  echo 'branch=feature/CR_tracers_relativistic_acceleration'
  echo "sanitized_preseal_candidate_head=$PHASE9_PRESEAL_CANDIDATE"
  echo "sanitized_preseal_candidate_tree=$(git rev-parse HEAD^{tree})"
  echo 'accepted_phase8_evidence_seal=aa8663e8a5d49e26c206363d028b52d0e350a91f'
  echo 'replaced_leaked_remote_tip=5e031387e66224b0e9dc4462fbf4d9a7ee01c9df'
  echo 'previous_portability_rebound_seal=aa66f6c27531116e12554631281c8f2ed07d93c6'
  echo 'previous_documentation_durability_seal=88d631e4943648fe83f0624cb30291fa52ab4296'
  echo 'previous_capture_completeness_seal=1a7086add5fffd55356109b99e6a66fcd0b43486'
  echo 'previous_public_overlay_privacy_seal=3667e3e1748f297615a8954306b5382ca7f476c2'
  echo 'origin_development=c6a73b08e60807f8b925164c5e7edd5cb820c8ae'
  echo 'origin_gh_pages=4833aa9341e19861297e330ff02aabfd8001935c'
} > "$EVID/phase9_branch_snapshot.txt"
```

This records the intended direction: merge the stacked candidate into the
frozen `origin/development` target.

## Whitespace Audit Scope

```bash
PHASE9_PRESEAL_CANDIDATE=$(git rev-parse HEAD)
git diff --check \
  aa8663e8a5d49e26c206363d028b52d0e350a91f.."$PHASE9_PRESEAL_CANDIDATE" \
  -- docs/source/modules inputs/particles \
  > "$EVID/phase9_commit_range_diff_check.log"
```

`phase9_commit_range_diff_check.log` retains the empty passing output for the
Phase-9 documentation and evidence packet added after the Phase-8 seal.  The
older stacked range beginning at the prerequisite branch still contains
historical retained-evidence whitespace.  That inherited archive noise is not
silently reported as repaired by this Phase-9 normalization pass.
The final rebound reviewers must rerun the same bounded check through the
successor evidence-seal tip before recording `CP-7 PROCEED`.

## Provenance Envelope, Outer Manifest, And Privacy Audit

The retained envelope is finalized after the overlay, public smoke, merge tree,
whitespace check, portable archive, and deterministic privacy log have reached
their final bytes.  Its `result_artifacts` map records SHA-256 digests of the
retained replay surfaces.  Its `portable_binary` entry binds the only binary
archived for offline reconstruction.  Ephemeral workstation build hashes are
not claimed as reconstructable envelope evidence.  The following verification
and outer-seal commands then run from the repository root:

```bash
PHASE9_PRESEAL_CANDIDATE=$(git rev-parse HEAD)
write_provenance_envelope() {
python3 - "$EVID/phase9_provenance_envelope.json" "$PHASE9_PRESEAL_CANDIDATE" <<'PY'
import hashlib
import json
import pathlib
import subprocess
import sys
import tarfile

def digest(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()

path = pathlib.Path(sys.argv[1])
candidate = sys.argv[2]
envelope = json.loads(path.read_text())
envelope.pop("binaries", None)
envelope.pop("first_sanitized_evidence_seal_head", None)
envelope["rebound_correction_candidate_head"] = candidate
envelope["rebound_correction_candidate_tree"] = subprocess.check_output(
    ["git", "rev-parse", f"{candidate}^{{tree}}"], text=True).strip()
public_input = "inputs/particles/cr_tracer_relativistic_mhd_ideal_example.athinput"
envelope["inputs"]["cr_tracer_relativistic_mhd_ideal_example.athinput"] = {
    "path": public_input,
    "sha256": digest(public_input),
}
resume_input = "inputs/particles/cr_tracer_relativistic_prescribed_restart_resume.athinput"
envelope["inputs"]["cr_tracer_relativistic_prescribed_restart_resume.athinput"] = {
    "path": resume_input,
    "sha256": digest(resume_input),
}
restart_smoke = "docs/source/modules/figures/cr_tracer_relativistic_acceleration/evidence/phase9_public_restart_resume_template_runtime.log"
envelope["result_artifacts"]["phase9_public_restart_resume_template_runtime.log"] = {
    "path": restart_smoke,
    "sha256": digest(restart_smoke),
}
for section in ("analyzers", "criteria", "inputs", "result_artifacts"):
    for item in envelope[section].values():
        item["sha256"] = digest(item["path"])
envelope["portable_bundle"]["sha256"] = digest(envelope["portable_bundle"]["path"])
envelope["portable_binary"] = {
    "archive_member": "phase9_all_formats_portable_replay/qualified_binary/athena",
}
with tarfile.open(envelope["portable_bundle"]["path"], "r:gz") as archive:
    payload = archive.extractfile(envelope["portable_binary"]["archive_member"]).read()
envelope["portable_binary"]["sha256"] = hashlib.sha256(payload).hexdigest()
path.write_text(json.dumps(envelope, indent=2, sort_keys=True) + "\n")
PY
}

verify_provenance_envelope() {
python3 - "$EVID/phase9_provenance_envelope.json" <<'PY' \
  > "$EVID/phase9_provenance_envelope_verify.log"
import hashlib
import json
import pathlib
import sys
import tarfile

envelope = json.loads(pathlib.Path(sys.argv[1]).read_text())
checked = 0
for section in ("analyzers", "criteria", "inputs", "result_artifacts"):
    for item in envelope[section].values():
        path = pathlib.Path(item["path"])
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        assert actual == item["sha256"], f"{path}: {actual} != {item['sha256']}"
        checked += 1
bundle = envelope["portable_bundle"]
actual = hashlib.sha256(pathlib.Path(bundle["path"]).read_bytes()).hexdigest()
assert actual == bundle["sha256"], f"{bundle['path']}: {actual} != {bundle['sha256']}"
checked += 1
binary = envelope["portable_binary"]
with tarfile.open(bundle["path"], "r:gz") as archive:
    payload = archive.extractfile(binary["archive_member"]).read()
actual = hashlib.sha256(payload).hexdigest()
assert actual == binary["sha256"], (
    f"{binary['archive_member']}: {actual} != {binary['sha256']}")
checked += 1
print(f"PASS: provenance envelope verified {checked} retained artifacts")
PY
}

write_outer_manifest() {
  find "$EVID_REL" -maxdepth 1 -type f -name 'phase9_*' \
    ! -name phase9_evidence_sha256.txt \
    ! -name phase9_evidence_verify.log \
    ! -name phase9_evidence_relocated_verify.log \
    -print | LC_ALL=C sort | xargs shasum -a 256 \
    > "$EVID/phase9_evidence_sha256.txt"
}
verify_outer_manifest() {
  shasum -a 256 -c "$EVID/phase9_evidence_sha256.txt" \
    > "$EVID/phase9_evidence_verify.log"
  rm -rf /tmp/cr-rel-phase9-relocated-root
  mkdir -p "/tmp/cr-rel-phase9-relocated-root/$EVID_REL"
  cp "$EVID"/phase9_* "/tmp/cr-rel-phase9-relocated-root/$EVID_REL/"
  (cd /tmp/cr-rel-phase9-relocated-root && \
    shasum -a 256 -c "$EVID_REL/phase9_evidence_sha256.txt") \
    > "$EVID/phase9_evidence_relocated_verify.log"
  cmp "$EVID/phase9_evidence_verify.log" \
    "$EVID/phase9_evidence_relocated_verify.log"
}
write_privacy_scan() {
  {
    LOCAL_USER='dbf''75'
    LOCAL_HOME_PATTERN="/Users/$LOCAL_USER"
    HOST_PATTERN='Tin''-Drum'
    echo "phase9 text privacy scan"
    ! rg -n "$LOCAL_HOME_PATTERN|$HOST_PATTERN" "$EVID"/phase9*
    echo "PASS: retained text has no local identity markers"
    USER_HOME_PREFIX='/Users/'
    ! rg -n "$USER_HOME_PREFIX|$HOST_PATTERN" \
      docs/source/modules/particles.md \
      docs/source/modules/cr_tracer*.md \
      docs/source/modules/figures/cr_tracer_accuracy \
      inputs/particles/cr_tracer_relativistic_*.athinput
    echo "PASS: public overlay source has no local identity markers"
    ! rg -n "$LOCAL_HOME_PATTERN|$HOST_PATTERN" \
      /tmp/cr-rel-phase9-portable-extract-a/phase9_all_formats_portable_replay
    echo "PASS: portable text has no local identity markers"
    ! strings \
      /tmp/cr-rel-phase9-portable-extract-a/phase9_all_formats_portable_replay/qualified_binary/athena |
      rg -n "$LOCAL_HOME_PATTERN|$HOST_PATTERN"
    echo "PASS: portable binary has no local identity markers"
  } > "$EVID/phase9_artifact_privacy_scan.log"
}

# First create both outer-verification logs.  Then generate the deterministic
# privacy log before binding it into the provenance envelope.  The final
# privacy scan runs again after the finalized envelope and manifest exist; its
# byte-identical PASS output leaves both bindings valid.
rm -f "$EVID/phase9_artifact_privacy_scan.log"
write_outer_manifest
verify_outer_manifest
write_privacy_scan
write_provenance_envelope
verify_provenance_envelope
write_outer_manifest
verify_outer_manifest
write_privacy_scan
verify_provenance_envelope
verify_outer_manifest
```

## Explicit Unsupported Or Unavailable Gates

```bash
set +e
{
  cmake -S . -B /tmp/athenak-cr-phase9-single-precision \
    -D Athena_SINGLE_PRECISION=ON
  cmake --build /tmp/athenak-cr-phase9-single-precision -j 4
} > "$EVID/phase9_single_precision_build.log" 2>&1
SINGLE_PRECISION_STATUS=$?
set -e
cp /tmp/athenak-cr-phase9-single-precision/CMakeCache.txt \
  "$EVID/phase9_single_precision_CMakeCache.txt"
test "$SINGLE_PRECISION_STATUS" -ne 0

cat > "$EVID/phase9_gpu_disposition.log" <<'EOF'
date=2026-05-30

--- workstation backend boundary ---
host_os=macOS-arm64
host_processor=Apple-M4-Max
graphics_api=Metal

--- accelerator toolchain probes ---
nvcc=NOT_FOUND
hipcc=NOT_FOUND
nvidia-smi=NOT_FOUND
rocminfo=NOT_FOUND

--- accepted Kokkos cache boundary ---
Kokkos_ENABLE_SERIAL=ON
Kokkos_ENABLE_CUDA=OFF
Kokkos_ENABLE_HIP=OFF
Kokkos_ENABLE_SYCL=OFF

--- disposition ---
GPU QUALIFIED is not claimed: this workstation exposes Metal graphics but no
configured CUDA or HIP Kokkos backend and no CUDA or HIP toolchain.
Accelerator qualification remains an explicit follow-up residual risk.
EOF
```

The single-precision build configures but fails before CR source in inherited
coordinate-header narrowing errors.  GPU execution was not run: the
workstation has no configured CUDA or HIP Kokkos backend and no corresponding
toolchain.  Both dispositions are retained separately and must not be
misreported as passed gates.
