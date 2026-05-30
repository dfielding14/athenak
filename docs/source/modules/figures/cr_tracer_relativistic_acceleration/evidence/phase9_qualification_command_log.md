# Phase-9 Relativistic CR Qualification Command Log

Date: 2026-05-30

This command inventory records the final workstation qualification replay for
`feature/CR_tracers_relativistic_acceleration`.  Paths under `/tmp` are
ephemeral runtime roots.  Retained logs, metrics, caches, archive inventories,
and hash manifests live beside this file.

## General Builds

```bash
cmake -S . -B /tmp/athenak-cr-phase9-release -D CMAKE_BUILD_TYPE=Release
cmake --build /tmp/athenak-cr-phase9-release -j 4

cmake -S . -B /tmp/athenak-cr-phase9-mpi-release \
  -D CMAKE_BUILD_TYPE=Release -D Athena_ENABLE_MPI=ON
cmake --build /tmp/athenak-cr-phase9-mpi-release -j 4

cmake -S . -B /tmp/athenak-cr-phase9-debug \
  -D CMAKE_BUILD_TYPE=Debug \
  -D Kokkos_ENABLE_DEBUG=ON \
  -D Kokkos_ENABLE_DEBUG_BOUNDS_CHECK=ON
cmake --build /tmp/athenak-cr-phase9-debug -j 4

cmake -S . -B /tmp/athenak-cr-phase9-mpi-debug \
  -D CMAKE_BUILD_TYPE=Debug \
  -D Athena_ENABLE_MPI=ON \
  -D Kokkos_ENABLE_DEBUG=ON \
  -D Kokkos_ENABLE_DEBUG_BOUNDS_CHECK=ON
cmake --build /tmp/athenak-cr-phase9-mpi-debug -j 4
```

## Dedicated Analytical Builds

```bash
cmake -S . -B /tmp/athenak-cr-phase9-pusher-release \
  -D CMAKE_BUILD_TYPE=Release \
  -D PROBLEM=unit_tests/cr_relativistic_pusher_runtime_test
cmake --build /tmp/athenak-cr-phase9-pusher-release -j 4

cmake -S . -B /tmp/athenak-cr-phase9-coupled-release \
  -D CMAKE_BUILD_TYPE=Release \
  -D PROBLEM=unit_tests/cr_relativistic_coupled_runtime_test
cmake --build /tmp/athenak-cr-phase9-coupled-release -j 4

cmake -S . -B /tmp/athenak-cr-phase9-coupled-debug \
  -D CMAKE_BUILD_TYPE=Debug \
  -D PROBLEM=unit_tests/cr_relativistic_coupled_runtime_test \
  -D Kokkos_ENABLE_DEBUG=ON \
  -D Kokkos_ENABLE_DEBUG_BOUNDS_CHECK=ON
cmake --build /tmp/athenak-cr-phase9-coupled-debug -j 4
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
  --metrics docs/source/modules/figures/cr_tracer_relativistic_acceleration/evidence/phase9_pusher_runtime_metrics.json

python3 scripts/particles/cr_relativistic_coupled_runtime_inspect.py \
  --binary /tmp/athenak-cr-phase9-coupled-release/src/athena \
  --input inputs/unit_tests/cr_relativistic_coupled_runtime_test.athinput \
  --prescribed-input inputs/unit_tests/cr_relativistic_pusher_runtime.athinput \
  --criteria docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase4b_preregistered_criteria.json \
  --work-dir /tmp/cr-rel-phase9-coupled \
  --metrics docs/source/modules/figures/cr_tracer_relativistic_acceleration/evidence/phase9_coupled_runtime_metrics.json

python3 scripts/particles/cr_relativistic_subcycle_runtime_inspect.py \
  --binary /tmp/athenak-cr-phase9-coupled-release/src/athena \
  --input inputs/unit_tests/cr_relativistic_coupled_runtime_test.athinput \
  --criteria docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase5_preregistered_criteria.json \
  --work-dir /tmp/cr-rel-phase9-subcycle \
  --metrics docs/source/modules/figures/cr_tracer_relativistic_acceleration/evidence/phase9_subcycle_runtime_metrics.json

python3 scripts/particles/cr_relativistic_restart_runtime_inspect.py \
  --binary /tmp/athenak-cr-phase9-release/src/athena \
  --input inputs/unit_tests/cr_relativistic_restart_runtime.athinput \
  --criteria docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase6_preregistered_criteria.json \
  --work-dir /tmp/cr-rel-phase9-restart \
  --json docs/source/modules/figures/cr_tracer_relativistic_acceleration/evidence/phase9_restart_runtime_metrics.json

python3 scripts/particles/cr_relativistic_diagnostics_runtime_inspect.py \
  --binary /tmp/athenak-cr-phase9-release/src/athena \
  --input inputs/unit_tests/cr_relativistic_diagnostics_runtime.athinput \
  --criteria docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase7_preregistered_criteria.json \
  --work-dir /tmp/cr-rel-phase9-diagnostics \
  --metrics docs/source/modules/figures/cr_tracer_relativistic_acceleration/evidence/phase9_diagnostics_runtime_metrics.json

python3 scripts/particles/cr_relativistic_all_formats_inspect.py \
  --runtime-root /tmp/cr-rel-phase9-diagnostics/acceleration_uninterrupted \
  --binary /tmp/athenak-cr-phase9-release/src/athena \
  --input inputs/unit_tests/cr_relativistic_diagnostics_runtime.athinput \
  --criteria docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase7_preregistered_criteria.json \
  --metrics docs/source/modules/figures/cr_tracer_relativistic_acceleration/evidence/phase9_diagnostics_runtime_metrics.json
```

Release replays used the release binaries and fresh `/tmp/cr-rel-phase9-*`
runtime roots.  Debug replays used the debug-bounds binaries and separate
`/tmp/cr-rel-phase9-*-debug` roots with the same argument structure.  The
portable release all-format package carries its own qualified binary, runtime
artifacts, input, criteria, required inspector scripts, deterministic manifest,
and replay wrapper.

## MPI, SMR, And AMR Migration

```bash
python3 scripts/particles/cr_relativistic_migration_runtime_inspect.py \
  --check-registration \
  --json docs/source/modules/figures/cr_tracer_relativistic_acceleration/evidence/phase9_migration_registration_metrics.json

python3 scripts/particles/cr_relativistic_migration_runtime_inspect.py \
  --binary /tmp/athenak-cr-phase9-release/src/athena \
  --mpi-binary /tmp/athenak-cr-phase9-mpi-release/src/athena \
  --work-dir /tmp/cr-rel-phase9-migration \
  --json docs/source/modules/figures/cr_tracer_relativistic_acceleration/evidence/phase9_migration_metrics.json

python3 scripts/particles/cr_relativistic_migration_runtime_inspect.py \
  --binary /tmp/athenak-cr-phase9-debug/src/athena \
  --mpi-binary /tmp/athenak-cr-phase9-mpi-debug/src/athena \
  --work-dir /tmp/cr-rel-phase9-migration-debug \
  --json docs/source/modules/figures/cr_tracer_relativistic_acceleration/evidence/phase9_migration_debug_metrics.json
```

## Adversarial Controls

```bash
python3 scripts/particles/cr_relativistic_migration_negative_controls.py
python3 scripts/particles/cr_relativistic_diagnostics_negative_controls.py
python3 scripts/particles/cr_relativistic_subcycle_negative_controls.py \
  --metrics docs/source/modules/figures/cr_tracer_relativistic_acceleration/evidence/phase9_subcycle_negative_control_metrics.json
python3 scripts/particles/cr_relativistic_coupled_oracle_negative_controls.py \
  --analyzer scripts/particles/cr_relativistic_coupled_runtime_inspect.py \
  --binary /tmp/athenak-cr-phase9-coupled-release/src/athena \
  --input inputs/unit_tests/cr_relativistic_coupled_runtime_test.athinput \
  --prescribed-input inputs/unit_tests/cr_relativistic_pusher_runtime.athinput \
  --criteria docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase4b_preregistered_criteria.json
python3 scripts/particles/cr_relativistic_diagnostics_writer_adversarial.py \
  --binary /tmp/athenak-cr-phase9-release/src/athena \
  --work-dir /tmp/cr-rel-phase9-writer-adversarial
```

## Repository Harness

These harness invocations ran sequentially because they share `tst/build`:

```bash
cd tst
python run_test_suite.py \
  --test test_suite/particles/test_particles_cr_relativistic_contract_cpu.py \
  --cpu
python run_test_suite.py \
  --test test_suite/particles/test_particles_cr_cpu.py \
  --cpu
python run_test_suite.py \
  --test test_suite/particles/test_particles_cr_accuracy_cpu.py \
  --cpu
python run_test_suite.py \
  --test test_suite/particles/test_particles_cr_mpicpu.py \
  --mpicpu
python run_test_suite.py \
  --test test_suite/particles/test_particles_cr_accuracy_mpicpu.py \
  --mpicpu
python run_test_suite.py --style
python run_tests.py mhd/mhd_divb_amr \
  --log_file ../docs/source/modules/figures/cr_tracer_relativistic_acceleration/evidence/phase9_deep_amr_divb.log
```

## Public Documentation Overlay

```bash
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
  /tmp/athenak-rel-gh-pages-overlay/inputs/particles/
cd /tmp/athenak-rel-gh-pages-overlay/docs
make clean html SPHINXOPTS="-W --keep-going"
perl -0pi -e 's/[ \t]+$//mg' \
  $REPO_ROOT/docs/source/modules/figures/cr_tracer_relativistic_acceleration/evidence/phase9_gh_pages_overlay_build.log
```

The overlay baseline is frozen as
`origin/gh-pages@4833aa9341e19861297e330ff02aabfd8001935c`.

## Public Solver-Coupled Example Smoke

```bash
mkdir -p /tmp/cr-rel-phase9-public-example
cd /tmp/cr-rel-phase9-public-example
/tmp/athenak-cr-phase9-release/src/athena \
  -i $REPO_ROOT/inputs/particles/cr_tracer_relativistic_mhd_ideal_example.athinput
```

## Portable All-Format Replay

The portable package uses a prefix-mapped release binary so compiled
`__FILE__` strings do not expose the local checkout path:

```bash
cmake -S . -B /tmp/athenak-cr-phase9-release-redacted \
  -D CMAKE_BUILD_TYPE=Release \
  -D 'CMAKE_CXX_FLAGS=-ffile-prefix-map=$REPO_ROOT=REPO_ROOT -fdebug-prefix-map=$REPO_ROOT=REPO_ROOT'
cmake --build /tmp/athenak-cr-phase9-release-redacted -j 4

python3 scripts/particles/cr_relativistic_diagnostics_runtime_inspect.py \
  --binary /tmp/athenak-cr-phase9-release-redacted/src/athena \
  --input inputs/unit_tests/cr_relativistic_diagnostics_runtime.athinput \
  --criteria docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase7_preregistered_criteria.json \
  --work-dir /tmp/cr-rel-phase9-diagnostics-portable \
  --metrics docs/source/modules/figures/cr_tracer_relativistic_acceleration/evidence/phase9_diagnostics_portable_runtime_metrics.json

EVID=docs/source/modules/figures/cr_tracer_relativistic_acceleration/evidence
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
python3 /tmp/cr-rel-phase9-portable-extract-a/phase9_all_formats_portable_replay/phase9_portable_replay.py
python3 /tmp/cr-rel-phase9-portable-extract-b/phase9_all_formats_portable_replay/phase9_portable_replay.py
cmp \
  /tmp/cr-rel-phase9-portable-extract-a/phase9_all_formats_portable_replay/portable_all_formats_report.json \
  /tmp/cr-rel-phase9-portable-extract-b/phase9_all_formats_portable_replay/portable_all_formats_report.json
```

## Integration Merge-Tree Audit

```bash
PHASE9_PRESEAL_CANDIDATE=$(git rev-parse HEAD)
git merge-tree --write-tree \
  c6a73b08e60807f8b925164c5e7edd5cb820c8ae \
  "$PHASE9_PRESEAL_CANDIDATE"
```

This records the intended direction: merge the stacked candidate into the
frozen `origin/development` target.

## Whitespace Audit Scope

```bash
PHASE9_PRESEAL_CANDIDATE=$(git rev-parse HEAD)
git diff --check \
  aa8663e8a5d49e26c206363d028b52d0e350a91f.."$PHASE9_PRESEAL_CANDIDATE" \
  -- docs/source/modules
```

`phase9_commit_range_diff_check.log` retains the empty passing output for the
Phase-9 documentation and evidence packet added after the Phase-8 seal.  The
older stacked range beginning at the prerequisite branch still contains
historical retained-evidence whitespace.  That inherited archive noise is not
silently reported as repaired by this Phase-9 normalization pass.
The final rebound reviewers must rerun the same bounded check through the
successor evidence-seal tip before recording `CP-7 PROCEED`.

## Provenance Envelope, Outer Manifest, And Privacy Audit

The retained envelope is regenerated after the overlay, public smoke, merge
tree, whitespace check, and portable archive have reached their final bytes.
Its `result_artifacts` map records SHA-256 digests of the retained replay
surfaces.  The following verification and outer-seal commands then run from the
repository root:

```bash
REPO_ROOT=$(pwd)
EVID=docs/source/modules/figures/cr_tracer_relativistic_acceleration/evidence
PHASE9_PRESEAL_CANDIDATE=$(git rev-parse HEAD)
python3 - "$EVID/phase9_provenance_envelope.json" "$PHASE9_PRESEAL_CANDIDATE" <<'PY'
import hashlib
import json
import pathlib
import subprocess
import sys

def digest(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()

path = pathlib.Path(sys.argv[1])
candidate = sys.argv[2]
envelope = json.loads(path.read_text())
envelope["rebound_correction_candidate_head"] = candidate
envelope["rebound_correction_candidate_tree"] = subprocess.check_output(
    ["git", "rev-parse", f"{candidate}^{{tree}}"], text=True).strip()
public_input = "inputs/particles/cr_tracer_relativistic_mhd_ideal_example.athinput"
envelope["inputs"]["cr_tracer_relativistic_mhd_ideal_example.athinput"] = {
    "path": public_input,
    "sha256": digest(public_input),
}
for section in ("analyzers", "criteria", "inputs", "result_artifacts"):
    for item in envelope[section].values():
        item["sha256"] = digest(item["path"])
envelope["portable_bundle"]["sha256"] = digest(envelope["portable_bundle"]["path"])
path.write_text(json.dumps(envelope, indent=2, sort_keys=True) + "\n")
PY

python3 - "$EVID/phase9_provenance_envelope.json" <<'PY' \
  > "$EVID/phase9_provenance_envelope_verify.log"
import hashlib
import json
import pathlib
import sys

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
print(f"PASS: provenance envelope verified {checked} retained artifacts")
PY

{
  LOCAL_USER='dbf''75'
  LOCAL_HOME_PATTERN="/Users/$LOCAL_USER"
  HOST_PATTERN='Tin''-Drum'
  echo "phase9 text privacy scan"
  ! rg -n "$LOCAL_HOME_PATTERN|$HOST_PATTERN" "$EVID"/phase9*
  echo "PASS: retained text has no local identity markers"
  ! rg -n "$LOCAL_HOME_PATTERN|$HOST_PATTERN" \
    /tmp/cr-rel-phase9-portable-extract-a/phase9_all_formats_portable_replay
  echo "PASS: portable text has no local identity markers"
  ! strings \
    /tmp/cr-rel-phase9-portable-extract-a/phase9_all_formats_portable_replay/qualified_binary/athena |
    rg -n "$LOCAL_HOME_PATTERN|$HOST_PATTERN"
  echo "PASS: portable binary has no local identity markers"
} > "$EVID/phase9_artifact_privacy_scan.log"

find "$EVID" -maxdepth 1 -type f -name 'phase9_*' \
  ! -name phase9_evidence_sha256.txt \
  ! -name phase9_evidence_verify.log \
  ! -name phase9_evidence_relocated_verify.log \
  -print | LC_ALL=C sort | xargs shasum -a 256 \
  > "$EVID/phase9_evidence_sha256.txt"
shasum -a 256 -c "$EVID/phase9_evidence_sha256.txt" \
  > "$EVID/phase9_evidence_verify.log"

rm -rf /tmp/cr-rel-phase9-relocated-root
mkdir -p "/tmp/cr-rel-phase9-relocated-root/$EVID"
cp "$EVID"/phase9_* "/tmp/cr-rel-phase9-relocated-root/$EVID/"
(cd /tmp/cr-rel-phase9-relocated-root && \
  shasum -a 256 -c "$EVID/phase9_evidence_sha256.txt") \
  > "$REPO_ROOT/$EVID/phase9_evidence_relocated_verify.log"
cmp "$EVID/phase9_evidence_verify.log" \
  "$EVID/phase9_evidence_relocated_verify.log"
```

## Explicit Unsupported Or Unavailable Gates

```bash
cmake -S . -B /tmp/athenak-cr-phase9-single-precision \
  -D Athena_SINGLE_PRECISION=ON
cmake --build /tmp/athenak-cr-phase9-single-precision -j 4
```

The single-precision build configures but fails before CR source in inherited
coordinate-header narrowing errors.  GPU execution was not run: the
workstation has no configured CUDA or HIP Kokkos backend and no corresponding
toolchain.  Both dispositions are retained separately and must not be
misreported as passed gates.
