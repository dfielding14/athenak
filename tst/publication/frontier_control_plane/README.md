# Frontier PIC Control Plane

These scripts implement the source-controlled accounting and submission
foundation required by `../PIC_PRODUCTION_READINESS_PLAN.md`. The mutable
repository copy may install a version, but it may not create manifests,
initialize the real ledger, reserve node-hours, submit jobs or reconcile jobs.
Run those operations only from a checksummed snapshot under
`$PIC_ROOT/control_plane/<version>/`.

The ledger is an append-only canonical-JSON event chain. Every primary event is
written locally, copied directly to the mounted Project Home mirror with the
exact reviewed `filesystem_copy` transport, and followed by exactly one local
non-recursive mirror receipt. Normal reads and appends reject a missing,
duplicate or divergent receipt and reject any Project Home drift.
`node_hours.csv` is a derived RFC-4180 index only.
Ledger mutation and policy promotion additionally serialize through a
descriptor lock on the site-owned Orion project directory
`/lustre/orion/ast207`, outside the user-replaceable PIC tree. Inner pinned
directory locks and compatibility lock-file identity checks remain mandatory.

Project Home is the operational ledger and control-plane mirror only. The user
selected the authorized Orion root as the sole bulk-evidence location for
simulation output, immutable sign-off bundles, private reference-artifact
staging and restore-drill evidence. This intentionally removes Kronos, DTN and
Globus from the execution boundary. Orion-only retention carries a documented
durability risk and is not represented as an institutional archive.

No real ledger genesis or reservation is permitted until a reviewed
`../readiness/storage_policy.json` is atomically promoted to the canonical
`${PIC_ROOT}/policy/storage_policy.json`, mirrored to Project Home, and bound by
matching read-only `policy/active_promotion.json` records. The active promotion
record binds the installed control-plane version and exact policy SHA-256.
Genesis and reservation read only that anchored policy and fail closed. The
exact authorized paths, `AST207` account, `batch` partition, serial-submission
policy, 10000-node-hour cap, Project Home usage, `filesystem_copy` ledger
transport and Orion bulk-evidence policy are also checked in code.
For a fresh ledger only, the promoted pre-genesis policy must explicitly open
`ledger_genesis_allowed=true` with no initialized genesis fields. That state
authorizes one-time genesis initialization only; it is not a reservation or
submission unlock. After initialization, the active policy must close
`ledger_genesis_allowed=false` and bind the exact genesis event and mirror
receipt digests before any reservation. Every non-empty ledger read also
requires matching read-only `ledger/genesis_anchor.json` files in Orion and
Project Home.

Registered science reservations have an additional fail-closed gate. They are
blocked while `science_submission_freeze.status` is
`pending_clean_candidate_freeze`, and are allowed only when the reviewed policy
authorizes one exact Orion clean-candidate manifest path and SHA-256 digest.

## Install

Install the same frozen control-plane version in Orion and Project Home only
from a reviewed clean Git commit. The production installer rejects untracked or
modified control-plane source files:

```bash
export PIC_ROOT=/lustre/orion/ast207/proj-shared/dfielding/PIC
export PROJECT_HOME_MIRROR_ROOT=/ccs/proj/ast207/proj-shared/PIC
PYTHON=/opt/cray/pe/python/3.11.7/bin/python3

"$PYTHON" -I tst/publication/frontier_control_plane/run_control_plane.py \
  install_control_plane.py \
  --pic-root /lustre/orion/ast207/proj-shared/dfielding/PIC
"$PYTHON" -I tst/publication/frontier_control_plane/run_control_plane.py \
  install_control_plane.py \
  --pic-root /ccs/proj/ast207/proj-shared/PIC
```

Both commands must print the same digest. Installation retains a pinned
descriptor for the authorized output parent, creates and populates its staging
directory relative to that descriptor, syncs every staged file and directory,
publishes the completed version with a descriptor-relative atomic rename,
syncs its parent directory, verifies that the lexical parent still names the
pinned directory, and marks the installed snapshot read-only. The same pinned
descriptor-relative publication rule applies to clean-candidate freezes and
pre-submit manifest snapshots. After the authenticated storage preflight is
complete, update the reviewed storage policy to record that digest and the
successful probes. Freeze `project_home_ledger_mirror_transport` as
`filesystem_copy`. Keep the Orion-only bulk-evidence selection and documented
durability risk explicit. Do not weaken the policy merely to create genesis.

Promote the reviewed policy through the installed control plane:

```bash
export VERSION=<reviewed-control-plane-digest>
export CONTROL_PLANE_DIR="${PIC_ROOT}/control_plane/${VERSION}"
CONTROL_PLANE=("$PYTHON" -I "${CONTROL_PLANE_DIR}/run_control_plane.py")

"${CONTROL_PLANE[@]}" promote_active_policy.py \
  --reviewed-policy /ccs/home/dfielding/athenak-pic/tst/publication/readiness/storage_policy.json
```

Promotion validates the reviewed policy, serializes promoters, durably
publishes each read-only Orion and Project Home file, and performs a final
coherent reread of both policy copies and both promotion records. Publication
cannot be collectively atomic across the two filesystems, so readers fail
closed while any partial generation is visible. A direct edit to the
repository policy or either active copy is not an authorization.

When promoting a successor over the historical pre-anchor ledger, install the
successor in both roots and migrate its already-audited genesis before changing
the active policy:

```bash
"${CONTROL_PLANE[@]}" initialize_frontier_ledger.py \
  --ledger-jsonl "${PIC_ROOT}/ledger/node_hours.jsonl" \
  --ledger-csv "${PIC_ROOT}/ledger/node_hours.csv" \
  --mirror-receipts "${PIC_ROOT}/ledger/mirror_receipts.jsonl" \
  --mirror-jsonl "${PROJECT_HOME_MIRROR_ROOT}/ledger/node_hours.jsonl" \
  --migrate-existing-anchor
```

This transition reads the mirrored predecessor policy anchor, requires its
audited genesis and receipt digests to match the existing chains, and creates
the paired read-only anchors exactly once. It never appends a second genesis.
If interruption leaves only one anchor visible, rerun the same audited
`--migrate-existing-anchor` command. It validates the surviving read-only
anchor against both chains and the active predecessor policy before publishing
only the missing mate.

## Genesis

Initialize the ledger exactly once through the installed version:

```bash
export VERSION=<reviewed-control-plane-digest>
export CONTROL_PLANE_DIR="${PIC_ROOT}/control_plane/${VERSION}"

"${CONTROL_PLANE[@]}" initialize_frontier_ledger.py \
  --ledger-jsonl "${PIC_ROOT}/ledger/node_hours.jsonl" \
  --ledger-csv "${PIC_ROOT}/ledger/node_hours.csv" \
  --mirror-receipts "${PIC_ROOT}/ledger/mirror_receipts.jsonl" \
  --mirror-jsonl "${PROJECT_HOME_MIRROR_ROOT}/ledger/node_hours.jsonl" \
  --mirror-transport filesystem_copy \
  --notes "Reviewed Frontier PIC ledger genesis after Orion-only storage selection"
```

The initializer verifies both installed inventories, the reviewed policy,
authorized paths, exact Project Home `filesystem_copy` ledger transport and
exact Orion bulk-evidence policy. It creates a new ledger only when all ledger
state is absent. Re-running the same initializer may complete only an exact
interrupted fresh-bootstrap prefix; any divergent or later ledger state fails
closed. Every later append requires the explicit mirrored genesis event and
paired read-only anchors. After this initial command, update the reviewed policy with
`ledger_genesis_allowed=false` plus the exact initialized genesis fields,
promote that closed policy, and only then reserve work.

## Clean Candidate

Build and profile the clean source commit through the installed Orion control
plane. `write_orion_build_profile.py` is the build authority: it accepts only
the authorized source root, one profile ID and the expected full Git commit.
The Orion PIC root is fixed in the installed control plane; there is no
caller-selected build root, executable path, output path, log path or Project
Home build location:

```bash
ENV_FILE="${CONTROL_PLANE_DIR}/frontier_pic_environment.sh"
source "$ENV_FILE" || exit $?

"${CONTROL_PLANE[@]}" write_orion_build_profile.py \
  --source-root /ccs/home/dfielding/athenak-pic \
  --profile-id hip-mpi-release-paper-pic \
  --expected-git-commit <full-lowercase-git-commit>
```

The writer measures and requires the exact reviewed module stack before it
creates any build path. Source the installed profile in the build shell first;
an unactivated clean shell is intentionally rejected.

For `<commit12>` equal to the first twelve characters of that commit and
`<profile>` equal to the profile ID, the writer derives this Orion-only layout:

```text
${PIC_ROOT}/build/<commit12>/<profile>/source/
${PIC_ROOT}/build/<commit12>/<profile>/cmake/
${PIC_ROOT}/bin/<commit12>/<profile>/athena
${PIC_ROOT}/bin/<commit12>/<profile>/build_profile.json
${PIC_ROOT}/bin/<commit12>/<profile>/profile_receipt.json
${PIC_ROOT}/bin/<commit12>/<profile>/CMakeCache.txt
${PIC_ROOT}/bin/<commit12>/<profile>/modules.txt
${PIC_ROOT}/bin/<commit12>/<profile>/toolchain.txt
${PIC_ROOT}/bin/<commit12>/<profile>/build-invocations.json
${PIC_ROOT}/bin/<commit12>/<profile>/git_status.preconfigure.txt
${PIC_ROOT}/bin/<commit12>/<profile>/git_status.txt
${PIC_ROOT}/bin/<commit12>/<profile>/submodule_status.txt
${PIC_ROOT}/bin/<commit12>/<profile>/environment.allowlist.txt
${PIC_ROOT}/bin/<commit12>/<profile>/build-environment.json
${PIC_ROOT}/logs/build/<commit12>.<profile>.configure.log
${PIC_ROOT}/logs/build/<commit12>.<profile>.build.log
```

The writer rejects a dirty or drifting authorized source closure and refuses
to reuse an existing build directory, bin directory or configure/build log. It
clones the local authorized source into the derived Orion `source/` directory
with `--no-hardlinks --no-checkout`, checks out the exact commit detached, and
reconstructs every recursive submodule from the corresponding validated local
source submodule with its pinned detached commit. No network fetch or mutable
source-tree build is part of this path.

The production profile fixes `ROCM_PATH=/opt/rocm-6.2.4`, the Cray `CC` wrapper
and a minimal build subprocess environment. It unloads the inactive default
`darshan-runtime` module so its site-Spack pkg-config path cannot drift into
the closed compiler-wrapper environment. The writer rejects caller
overrides, records the exact retained `build-environment.json`, generates the
exact `build-invocations.json` artifact below and executes both argv arrays
directly without a caller-provided shell recipe, command file or `tee`
pipeline:

```json
{
  "build": [
    "/usr/bin/cmake",
    "--build",
    "${PIC_ROOT}/build/<commit12>/<profile>/cmake",
    "--parallel",
    "32"
  ],
  "configure": [
    "/usr/bin/cmake",
    "-S",
    "${PIC_ROOT}/build/<commit12>/<profile>/source",
    "-B",
    "${PIC_ROOT}/build/<commit12>/<profile>/cmake",
    "-DCMAKE_BUILD_TYPE=Release",
    "-DAthena_ENABLE_MPI=ON",
    "-DKokkos_ENABLE_HIP=ON",
    "-DKokkos_ARCH_ZEN3=ON",
    "-DKokkos_ARCH_AMD_GFX90A=ON",
    "-DCMAKE_CXX_COMPILER=/opt/cray/pe/craype/2.7.33/bin/CC",
    "-DCMAKE_CXX_FLAGS=-I/opt/rocm-6.2.4/include",
    "-DCMAKE_EXE_LINKER_FLAGS=-L/opt/rocm-6.2.4/lib -lamdhip64",
    "-DPROBLEM=built_in_pgens"
  ]
}
```

Immediately before configuration and immediately after the build, the writer
captures and requires empty
`git status --ignore-submodules=none --porcelain --untracked-files=all`
output from the fresh detached checkout. It retains those exact bytes as
`git_status.preconfigure.txt` and `git_status.txt`, captures recursive
`submodule_status.txt`, copies the built executable and CMake cache into the
derived bin directory, and records modules, the fixed toolchain description,
the redacted runtime environment allowlist and the minimal build subprocess
environment. All artifact writes are exclusive.

The writer then publishes schema-v3 `build_profile.json` and its adjacent
`profile_receipt.json` without replacing either file. The profile binds the
authorized and fresh source roots, source archive, raw commit object,
recursive-submodule closure, toolchain, exact `build-invocations.json` digest,
executable digest and hashes of all eleven retained provenance inputs. The receipt
additionally binds the installed control-plane version, exact profile path and
digest, source-bundle digest, fresh source root, invocation digest, pre/post
Git-status digests, configure/build log digests and executable path and digest.

Create one candidate freeze from those derived Orion artifacts:

```bash
"${CONTROL_PLANE[@]}" create_clean_candidate_freeze.py \
  --source-root /ccs/home/dfielding/athenak-pic \
  --executable "${PIC_ROOT}/bin/<git-commit12>/<profile>/athena" \
  --build-profile "${PIC_ROOT}/bin/<git-commit12>/<profile>/build_profile.json" \
  --build-profile-id hip-mpi-release-paper-pic
```

The generated input build profile is structured JSON and binds the generated
`git archive`, raw commit object, recursive-submodule source closure, exact
build invocation artifact, toolchain and executable:

```json
{
  "schema_version": 3,
  "profile_id": "hip-mpi-release-paper-pic",
  "authorized_source_root": "/ccs/home/dfielding/athenak-pic",
  "fresh_source_root": "${PIC_ROOT}/build/<git-commit12>/<profile>/source",
  "git_commit": "<full-git-commit>",
  "git_tree": "<git-tree>",
  "source_archive_sha256": "<sha256 of git archive HEAD>",
  "source_commit_sha256": "<sha256 of raw git cat-file commit HEAD>",
  "source_bundle_sha256": "<sha256 of canonical parent archive, raw commit object and recursive-submodule binding>",
  "toolchain": "<reviewed toolchain description>",
  "build_invocations_sha256": "<sha256 of exact build-invocations.json>",
  "executable_sha256": "<sha256 of athena>",
  "provenance_inputs": {
    "configure_log": {"path": "<documented Orion path>", "sha256": "<sha256>"},
    "build_log": {"path": "<documented Orion path>", "sha256": "<sha256>"},
    "cmake_cache": {"path": "<documented Orion path>", "sha256": "<sha256>"},
    "module_list": {"path": "<documented Orion path>", "sha256": "<sha256>"},
    "toolchain": {"path": "<documented Orion path>", "sha256": "<sha256>"},
    "build_invocations": {"path": "<documented Orion path>", "sha256": "<sha256>"},
    "git_status_preconfigure": {"path": "<documented Orion path>", "sha256": "<sha256>"},
    "git_status": {"path": "<documented Orion path>", "sha256": "<sha256>"},
    "submodule_status": {"path": "<documented Orion path>", "sha256": "<sha256>"},
    "environment_allowlist": {"path": "<documented Orion path>", "sha256": "<sha256>"},
    "build_environment": {"path": "<documented Orion path>", "sha256": "<sha256>"}
  },
  "submodules": [
    {
      "path": "kokkos",
      "archive_sha256": "<sha256 of git -C kokkos archive HEAD>",
      "commit_sha256": "<sha256 of raw git -C kokkos cat-file commit HEAD>",
      "git_commit": "<pinned kokkos commit>",
      "git_tree": "<kokkos HEAD tree>"
    }
  ]
}
```

The freeze creator revalidates the documented Orion layout, exact schema-v3
profile, adjacent receipt, installed control-plane generation, eleven provenance
hashes, empty pre/post Git-status captures, executable digest and recursively
clean pinned source closure. It rejects tracked symlink payloads, regenerates
the source archive and raw commit object plus deterministic archive and raw
commit object pairs for every clean pinned submodule, verifies each archive and
commit-object identity against the reconstructed Git tree identity, and
requires canonical author and committer headers in each raw commit object. It
stages the archives, profile, receipt, executable, exact eleven provenance bytes
and manifest under a temporary directory, marks the completed candidate
read-only, syncs every staged file and directory, atomically renames it beneath
`${PIC_ROOT}/clean_candidates/<freeze-id>/`, and syncs the parent directory.
Review that manifest and then
update `science_submission_freeze` to `status=authorized` with its exact
`manifest_path` and `manifest_sha256`, followed by another active-policy
promotion. The freeze creator never authorizes itself.

The fresh detached checkout, direct argv execution, exclusive paths, dual
clean-source captures, receipt and retained logs reduce stale, mixed-build and
ordinary caller-controlled mistakes. They do not cryptographically prove that
the compiler honored the recorded command or establish compiler semantics.
They also do not prevent a malicious process running as the same Unix UID from
modifying writable source or Orion state between checks. Reproducibility,
qualification and a stronger site-level trust boundary remain separate gates.

## Submit

Create each immutable pre-submit manifest with the installed
`create_pre_submit_manifest.py`, then invoke only the installed wrapper:

```bash
squeue -u "$USER" -h -o '%i|%P|%q|%T|%j|%k' \
  > "${PIC_ROOT}/jobs/<campaign>/queue_snapshot.txt"

MANIFEST="$(
  "${CONTROL_PLANE[@]}" create_pre_submit_manifest.py \
    --config "${PIC_ROOT}/jobs/<campaign>/pre_submit_config.json"
)"

"${CONTROL_PLANE_DIR}/submit_frontier_job.sh" "$MANIFEST"
```

The config must reference the exact authorized Orion root, an immutable Slurm
directive template containing `#SBATCH -A AST207`, `#SBATCH -p batch`, one
allowed QOS, node count, walltime and exactly
`#SBATCH -o ${PIC_ROOT}/logs/slurm/%x.%j.log`. Slurm stdout is not permitted
outside that dedicated PIC log path. The template body is never executed. It
must also reference the executable, input deck, environment profile, queue
snapshot, analysis scripts and timeout-margin artifact. Every config must
declare exactly one `submission_scope`. `registered_science` configs must also
reference the policy-authorized `clean_candidate_manifest`. Every config must declare
`job_script_executable_env=PIC_EXECUTABLE`. The creator promotes a read-only
submission snapshot through a pinned descriptor for its authorized Orion
parent. Every registered artifact directory is exactly
`${PIC_ROOT}/runs/<campaign>/<submission-id>`; reruns use new submission IDs
and never reuse an existing directory. Reservation rechecks the original
candidate manifest, archived source, clean Git commit and tree attestation, structured
build profile, exact candidate executable digest, and snapshotted executable
digest.

Every config must also define a closed `launch_contract`. The installed
trampoline converts each `athena` action into a fixed `/usr/bin/srun` argv whose
executable token is always the verified immutable Athena snapshot. Config data
may supply only the closed Athena argument grammar: `-i` immediately followed
by the snapshotted input-deck token, `-d` immediately followed by an
Orion-artifact-directory token, bounded `-n`/`-c` flags and path-free Athena
`section/name=value` overrides. Restart `-r` remains rejected until trusted
restart snapshots are implemented. It cannot select a shell, executable,
interpreter or arbitrary path-bearing literal. Before `srun`, an
inventory-bound helper applies the installed `frontier_pic_environment.sh`
module profile, writes the action-specific redacted runtime allowlist artifact
`<action_id>.environment.allowlist.txt` through a trampoline-opened inherited
file descriptor, and then `exec`s the fixed argv. Qualification manifests bind
that exact runtime artifact rather than an unrestricted environment dump. The
trampoline strips Bash startup hooks, exported functions and shell-option
startup state before launching that wrapper. The wrapper syncs the
runtime-allowlist descriptor, changes it to mode `0400`, syncs it again, syncs
its pinned parent-directory descriptor, closes both descriptors and removes
both environment bindings before executing the workload. Athena and its
descendants therefore cannot rewrite the captured evidence through the
inherited launch descriptors.
The snapshotted environment profile must match that installed trusted profile:

```json
{
  "schema_version": 1,
  "executor": "trusted_trampoline_athena_argv_v1",
  "pre_actions": [],
  "actions": [
    {
      "action_id": "f0-parser",
      "kind": "athena",
      "resources": {
        "nodes": 1,
        "tasks": 1,
        "cpus_per_task": 7,
        "gpus_per_task": 1,
        "gpu_bind": "closest"
      },
      "arguments": [
        {"literal": "-i"},
        {"snapshot_role": "input-deck"},
        {"literal": "-n"}
      ],
      "stdout_artifact": "athena_stdout.txt",
      "stderr_artifact": "athena_stderr.txt"
    }
  ],
  "post_actions": []
}
```

Optional `pre_actions` and `post_actions` are declarative built-ins only:
`snapshot_sha256`, `artifact_sha256` and `artifact_nonempty`. They cannot run
user code. Snapshotted Python analysis remains an offline evidence step after
the registered job, not a trusted launch action.

The existing F0/F1 shell files remain usable as reviewed Slurm-directive and
argv references, but their bodies are not executable launch contracts. Before
submitting them through this control plane, translate each Athena `srun` into a
structured `actions` entry. This covers the F0 parser run, the F1 gyro run and
both F1 coupling-coefficient runs. Legacy ROCm probes, GPU-mapping probes,
`ldd`, log collation and Python reduction do not run implicitly; capture them
as separately reviewed evidence steps. The trusted launch intentionally fails
closed instead of interpreting those shell bodies.

`frontier_admission_smoke` is a narrow pre-freeze exception for the registered
F0 parser-contract admission smoke only: campaign `f0_hipmpi_smoke`, test
`pic_parser_contract_guards`, evidence class
`frontier_f0_admission_smoke_candidate`, physical mode
`extended_mhd_pic_parser_contract`, debug QoS and a registered short
non-production job. The policy additionally binds the exact F0 Slurm script,
input deck, environment profile, single analysis-script digest and executable
snapshot digest, exact canonical launch-contract digest, with a one-node and
15-minute ceiling. Until that exact executable and contract are reviewed and
promoted, F0 remains
`pending_exact_executable_binding`. The exception cannot be used for registered
science evidence.

The timeout-margin artifact must include
`athena_walltime_seconds`, `scheduler_walltime_seconds`,
`environment_profile_sha256`, `measured_utc` and `expires_utc`. Its profile
checksum must match the snapshotted environment script and it must be valid at
reservation time. `site_policy_checked_utc` must be no more than 24 hours old.

The installed `frontier_pic_environment.sh` defaults to
`frontier_minimum_supported`. It also exposes the controlled
`frontier_xnack1_experimental` and `frontier_ofi_tuned_experimental` profiles
for matched Q-038 A/B evidence. `frontier_selected_production` is deliberately
undefined until those comparisons support an explicit reviewed selection.
The profile validates its selector before changing modules, fails closed if a
module operation fails, and exports `SLURM_EXPORT_ENV=ALL`. Qualification
accepts only the exact ordered redacted allowlist emitted immediately before
each trusted `srun`; extra, missing or duplicate keys are rejected.

The wrapper verifies its own installed inventory, reserves worst-case
node-hours under one lock before `sbatch`, records a recovery marker before the
reservation append, and submits only the installed trusted trampoline with
`sbatch --hold --export=NIL` through a closed environment pinned to the
`frontier` cluster. The scheduled runner receives its immutable bindings as
argv and the trampoline rechecks them against the mirrored ledger; it does not
reconstruct the submitter's login environment. The reservation ledger event
binds the active policy and active-promotion SHA-256 digests. Every wrapper
metadata read used to construct the scheduler request requires that live
mirrored reservation to remain in `reserved` state. The validator rechecks the
bound policy generation before dispatch, before scheduler-ID attachment and in
the trampoline before workload execution.

Immediately before the held `sbatch`, the wrapper durably advances the recovery
marker to `scheduler_dispatch_started`. Only a failure known to occur before
that marker may cancel the unused reservation automatically. At or after that
marker, an ambiguous failure retains the reserved accounting and pending marker
and requires reviewed Slurm inspection and reconciliation. After `sbatch`
returns a scheduler ID, the validator durably advances the marker through
`scheduler_job_id_received` and `submitted_not_attached`, verifies the held
scheduler job through `scontrol`, mirrors the scheduler-ID attachment and only
then runs `scontrol release`. A terminal job discovered before attachment is
reconciled from the durable marker.

Wrapper lookups and trampoline checks use the canonical mirrored ledger
reservation as their authority, rather than mutable snapshot attachments. The
trampoline rechecks the scheduled manifest checksum, reservation attachments,
submission ID, immutable snapshots, exact Slurm-template digest and exact
executable digest. It sets `PIC_EXECUTABLE` to the verified executable snapshot
and runs only structured Athena actions through its fixed `/usr/bin/srun` argv
builder. It never opens the snapshotted Slurm template as a shell script.

The registered path rejects symlink aliases beneath writable Orion and Project
Home roots for installed control planes, policy anchors, ledgers, manifests,
clean candidates, snapshotted dependencies and launch artifacts. Provenance and
scheduler subprocesses are pinned to reviewed absolute paths:
`/usr/bin/git`, `/usr/bin/squeue`, `/usr/bin/scontrol`, `/usr/bin/sacct`,
`/usr/bin/sbatch`, `/usr/bin/scancel`, `/usr/bin/srun` and
`/opt/cray/pe/python/3.11.7/bin/python3`. Reconciliation accepts only an explicit
allowlist of terminal Slurm states, queries `sacct --clusters=frontier`, and
fails closed on an unknown state.

Pinned descriptor-relative staging closes pathname-redirection races caused by
ancestor swaps during publication. This repository control plane still cannot
prevent a process running as the same Unix UID from modifying writable inputs
or published Orion files, invoking `sbatch` directly, or running an unrelated
executable. Site-level enforcement would require scheduler-side controls and a
stronger storage trust boundary. Such a direct submission is outside the
registered evidence path and must not be accepted as readiness evidence. The
trusted path reduces accidental and ordinary caller-controlled bypasses; it is
not a site security boundary against a malicious account holder.

Scientific qualification is a separate freeze step. A Frontier qualification
manifest is accepted only when its candidate matches the live policy and its
resource identity matches a completed `registered_science` reconciliation in
the mirrored ledger: reservation ID, submission ID, scheduler job ID,
control-plane version, clean-candidate digest, immutable pre-submit-manifest
path and checksum, and run-artifact directory must all agree.

## Reconcile

After the job reaches a terminal Slurm state, run:

```bash
"${CONTROL_PLANE[@]}" reconcile_frontier_job.py \
  --job-id <job-id> \
  --ledger-jsonl "${PIC_ROOT}/ledger/node_hours.jsonl" \
  --ledger-csv "${PIC_ROOT}/ledger/node_hours.csv" \
  --receipts-jsonl "${PIC_ROOT}/ledger/mirror_receipts.jsonl" \
  --mirror-jsonl "${PROJECT_HOME_MIRROR_ROOT}/ledger/node_hours.jsonl"
```

Reconciliation queries `sacct` itself and accepts only the allocation row whose
job ID, `AST207` account and `pic-reservation=<reservation-id>` comment match
the mirrored reservation. Caller-supplied terminal state, time and nodes are
not accepted.

If `pending_submission.json` remains, stop new submissions. Its durable states
include `reserved_not_submitted`, `scheduler_dispatch_started`,
`scheduler_job_id_received` and `submitted_not_attached`. If a reservation
intent was recorded but its primary append did not finish, or if a reservation
event was appended but its immutable attachments were interrupted before
dispatch, repair the marker through the installed validator:

```bash
"${CONTROL_PLANE[@]}" validate_and_reserve_frontier_job.py \
  repair-reservation-attachments \
  --ledger-jsonl "${PIC_ROOT}/ledger/node_hours.jsonl" \
  --ledger-csv "${PIC_ROOT}/ledger/node_hours.csv" \
  --receipts-jsonl "${PIC_ROOT}/ledger/mirror_receipts.jsonl" \
  --mirror-jsonl "${PROJECT_HOME_MIRROR_ROOT}/ledger/node_hours.jsonl" \
  --reservation-id <reservation-id>
```

If the primary Orion append completed but Project Home mirroring or the local
receipt append was interrupted, repair only the missing suffix through the
installed validator:

```bash
"${CONTROL_PLANE[@]}" validate_and_reserve_frontier_job.py \
  repair-ledger-mirror \
  --ledger-jsonl "${PIC_ROOT}/ledger/node_hours.jsonl" \
  --ledger-csv "${PIC_ROOT}/ledger/node_hours.csv" \
  --receipts-jsonl "${PIC_ROOT}/ledger/mirror_receipts.jsonl" \
  --mirror-jsonl "${PROJECT_HOME_MIRROR_ROOT}/ledger/node_hours.jsonl"
```

Repair accepts canonical Orion as authoritative only when Project Home and
receipt files are exact prefixes. It rejects divergence, corruption and
duplicate receipts. If a terminal scheduler job is present before attachment,
run the installed reconciler so it can durably attach from the marker and
record terminal accounting in one reviewed recovery path. If scheduler
attachment, cancellation or accounting remains ambiguous, stop submissions and
perform reviewed manual recovery. Recovery transitions intentionally remain
available when the storage policy is relocked. Terminal reconciliation cleanup
is idempotent: rerunning reconciliation for an already-recorded terminal event
clears a matching stale pending marker and returns the existing event without
duplicating accounting.
