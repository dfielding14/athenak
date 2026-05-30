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

Registered science reservations have an additional fail-closed gate. They are
blocked while `science_submission_freeze.status` is
`pending_clean_candidate_freeze`, and are allowed only when the reviewed policy
authorizes one exact Orion clean-candidate manifest path and SHA-256 digest.

## Install

Install the same frozen control-plane version in Orion and Project Home from
the mutable repository copy:

```bash
export PIC_ROOT=/lustre/orion/ast207/proj-shared/dfielding/PIC
export PROJECT_HOME_MIRROR_ROOT=/ccs/proj/ast207/proj-shared/PIC

python3 tst/publication/frontier_control_plane/install_control_plane.py \
  --pic-root /lustre/orion/ast207/proj-shared/dfielding/PIC
python3 tst/publication/frontier_control_plane/install_control_plane.py \
  --pic-root /ccs/proj/ast207/proj-shared/PIC
```

Both commands must print the same digest. Installation uses a temporary
directory, atomically renames the completed version into place, and marks the
installed snapshot read-only. After the authenticated storage preflight is
complete, update the reviewed storage policy to record that digest and the
successful probes. Freeze `project_home_ledger_mirror_transport` as
`filesystem_copy`. Keep the Orion-only bulk-evidence selection and documented
durability risk explicit. Do not weaken the policy merely to create genesis.

Promote the reviewed policy through the installed control plane:

```bash
export VERSION=<reviewed-control-plane-digest>
export CONTROL_PLANE_DIR="${PIC_ROOT}/control_plane/${VERSION}"

python3 "${CONTROL_PLANE_DIR}/promote_active_policy.py" \
  --reviewed-policy /ccs/home/dfielding/athenak-pic/tst/publication/readiness/storage_policy.json
```

Promotion validates the reviewed policy, atomically installs read-only Orion
and Project Home copies, and atomically writes matching read-only promotion
records. A direct edit to the repository policy or either active copy is not an
authorization.

## Genesis

Initialize the ledger exactly once through the installed version:

```bash
export VERSION=<reviewed-control-plane-digest>
export CONTROL_PLANE_DIR="${PIC_ROOT}/control_plane/${VERSION}"

python3 "${CONTROL_PLANE_DIR}/initialize_frontier_ledger.py" \
  --ledger-jsonl "${PIC_ROOT}/ledger/node_hours.jsonl" \
  --ledger-csv "${PIC_ROOT}/ledger/node_hours.csv" \
  --mirror-receipts "${PIC_ROOT}/ledger/mirror_receipts.jsonl" \
  --mirror-jsonl "${PROJECT_HOME_MIRROR_ROOT}/ledger/node_hours.jsonl" \
  --mirror-transport filesystem_copy \
  --notes "Reviewed Frontier PIC ledger genesis after Orion-only storage selection"
```

The initializer verifies both installed inventories, the reviewed policy,
authorized paths, exact Project Home `filesystem_copy` ledger transport, exact
Orion bulk-evidence policy and the absence of any previous ledger state. Every
later append requires the explicit mirrored genesis event.

## Clean Candidate

After a clean source commit has been built under Orion, create one candidate
freeze through the installed control plane:

```bash
python3 "${CONTROL_PLANE_DIR}/create_clean_candidate_freeze.py" \
  --source-root /ccs/home/dfielding/athenak-pic \
  --executable "${PIC_ROOT}/build/<profile>/athena" \
  --build-profile "${PIC_ROOT}/build/<profile>/build_profile.json" \
  --build-profile-id hip-mpi-release-paper-pic
```

The input build profile is structured JSON and must already bind the generated
`git archive`, toolchain, build command and executable:

```json
{
  "schema_version": 1,
  "profile_id": "hip-mpi-release-paper-pic",
  "source_archive_sha256": "<sha256 of git archive HEAD>",
  "toolchain": "<reviewed toolchain description>",
  "build_command": "<reviewed build command>",
  "executable_sha256": "<sha256 of athena>"
}
```

The creator forces `git status --ignore-submodules=none`, conservatively rejects
any submodule, generates a `git archive` for the declared HEAD, verifies the
archive commit identity and reconstructed Git tree identity, and requires an
exact structured-profile match. It stages the archive, profile, executable and
manifest under a temporary directory, marks the completed candidate read-only,
and atomically renames it beneath
`${PIC_ROOT}/clean_candidates/<freeze-id>/`. Review that manifest and then
update `science_submission_freeze` to `status=authorized` with its exact
`manifest_path` and `manifest_sha256`, followed by another active-policy
promotion. The freeze creator never authorizes itself.

## Submit

Create each immutable pre-submit manifest with the installed
`create_pre_submit_manifest.py`, then invoke only the installed wrapper:

```bash
squeue -u "$USER" -h -o '%i|%P|%q|%T|%j|%k' \
  > "${PIC_ROOT}/jobs/<campaign>/queue_snapshot.txt"

MANIFEST="$(
  python3 "${CONTROL_PLANE_DIR}/create_pre_submit_manifest.py" \
    --config "${PIC_ROOT}/jobs/<campaign>/pre_submit_config.json"
)"

"${CONTROL_PLANE_DIR}/submit_frontier_job.sh" "$MANIFEST"
```

The config must reference the exact authorized Orion root, an immutable Slurm
directive template containing `#SBATCH -A AST207`, `#SBATCH -p batch`, one
allowed QOS, node count, walltime and an Orion-rooted output path. The template
body is never executed. It must also reference the executable, input deck,
environment profile, queue snapshot, analysis scripts and timeout-margin
artifact. Every config must declare exactly one
`submission_scope`. `registered_science` configs must also reference the
policy-authorized `clean_candidate_manifest`. Every config must declare
`job_script_executable_env=PIC_EXECUTABLE`. The creator atomically promotes a
read-only submission snapshot. Reservation rechecks the original candidate
manifest, archived source, clean Git commit and tree attestation, structured
build profile, exact candidate executable digest, and snapshotted executable
digest.

Every config must also define a closed `launch_contract`. The installed
trampoline converts each `athena` action into a fixed `/usr/bin/srun` argv whose
executable token is always the verified immutable Athena snapshot. Config data
may supply only argument literals, the snapshotted input-deck token and
Orion-artifact-directory tokens. It cannot select a shell, executable or
interpreter. Before `srun`, an inventory-bound helper applies the installed
`frontier_pic_environment.sh` module profile and then `exec`s the fixed argv.
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
snapshot digest, with a one-node and 15-minute ceiling. Until that exact
executable is reviewed and promoted, F0 remains
`pending_exact_executable_binding`. The exception cannot be used for registered
science evidence.

The timeout-margin artifact must include
`athena_walltime_seconds`, `scheduler_walltime_seconds`,
`environment_profile_sha256`, `measured_utc` and `expires_utc`. Its profile
checksum must match the snapshotted environment script and it must be valid at
reservation time. `site_policy_checked_utc` must be no more than 24 hours old.

The wrapper verifies its own installed inventory, reserves worst-case
node-hours under one lock before `sbatch`, records a recovery marker before the
reservation append, and submits only the installed trusted trampoline. Wrapper
lookups and trampoline checks use the canonical mirrored ledger reservation as
their authority, rather than mutable snapshot attachments. The trampoline
rechecks the scheduled manifest checksum, reservation attachments, submission
ID, immutable snapshots, exact Slurm-template digest and exact executable
digest. It sets `PIC_EXECUTABLE` to the verified executable snapshot and runs
only structured Athena actions through its fixed `/usr/bin/srun` argv builder.
It never opens the snapshotted Slurm template as a shell script. After `sbatch`,
the wrapper queries `scontrol` and attaches a scheduler ID only when its live
account, reservation comment and state match. It leaves a persistent recovery
marker if attachment cannot be confirmed.

The registered path rejects symlink aliases beneath writable Orion and Project
Home roots for installed control planes, policy anchors, ledgers, manifests,
clean candidates, snapshotted dependencies and launch artifacts. Provenance and
scheduler subprocesses are pinned to reviewed absolute paths:
`/usr/bin/git`, `/usr/bin/squeue`, `/usr/bin/scontrol`, `/usr/bin/sacct`,
`/usr/bin/sbatch`, `/usr/bin/scancel`, `/usr/bin/srun` and
`/opt/cray/pe/python/3.11.7/bin/python3`. Reconciliation accepts only an explicit
allowlist of terminal Slurm states and fails closed on an unknown state.

This repository control plane cannot prevent the account holder from invoking
`sbatch` directly, running an unrelated executable, or modifying writable Orion
files between checks. Site-level enforcement would require scheduler-side
controls and a stronger storage trust boundary. Such a direct submission is
outside the registered evidence path and must not be accepted as readiness
evidence. The trusted path reduces accidental and ordinary caller-controlled
bypasses; it is not a site security boundary against a malicious account
holder.

## Reconcile

After the job reaches a terminal Slurm state, run:

```bash
python3 "${CONTROL_PLANE_DIR}/reconcile_frontier_job.py" \
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

If `pending_submission.json` remains, stop new submissions. If a reservation
intent was recorded but its primary append did not finish, or if a reservation
event was appended but its immutable attachments were interrupted, repair the
marker through the installed validator:

```bash
python3 "${CONTROL_PLANE_DIR}/validate_and_reserve_frontier_job.py" \
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
python3 "${CONTROL_PLANE_DIR}/validate_and_reserve_frontier_job.py" \
  repair-ledger-mirror \
  --ledger-jsonl "${PIC_ROOT}/ledger/node_hours.jsonl" \
  --ledger-csv "${PIC_ROOT}/ledger/node_hours.csv" \
  --receipts-jsonl "${PIC_ROOT}/ledger/mirror_receipts.jsonl" \
  --mirror-jsonl "${PROJECT_HOME_MIRROR_ROOT}/ledger/node_hours.jsonl"
```

Repair accepts canonical Orion as authoritative only when Project Home and
receipt files are exact prefixes. It rejects divergence, corruption and
duplicate receipts. If scheduler attachment, cancellation or accounting
remains ambiguous, stop submissions and perform reviewed manual recovery.
Recovery transitions intentionally remain available when the storage policy is
relocked.
