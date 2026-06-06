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
record binds the installed control-plane version, exact policy SHA-256, and one
unique canonical promotion UUID. Re-promoting byte-identical policy therefore
produces a different promotion digest, so exact policy/promotion compare-and-
swap checkpoints cannot be replayed after an A-to-B-to-A sequence. Policy
publication is one durable mirrored rollback transaction: the promoter creates
rollback anchors and mirrored prepared `.active_promotion_transaction.json`
markers under both policy-parent locks, publishes and validates all four active
anchors, then durably upgrades both markers to committed before cleanup. Each
marker binds the exact predecessor and successor digest for all four anchors.
Recovery finalizes any complete valid successor forward after validating all
four successor digests, the coherent active generation, its paired installed
controller, and any authorized clean candidate, including when the surviving
marker state is prepared or mixed prepared/committed. A complete successor that
fails semantic, controller, or authorized-candidate validation retains locked
recovery evidence and is never rolled back. Any committed marker over a partial
successor retains all transaction evidence for reviewed manual recovery. Only
a prepared-only strict partial reachable publication prefix is eligible for
rollback. Partial prepared recovery validates the
exact rollback predecessor and controller before mutation, restores in reverse
publication order, and revalidates the restored predecessor before deleting
evidence. Committed completion removes rollback anchors before transaction
markers, so an anchor-cleanup failure leaves readers closed behind a marker.
Active readers fail closed while either marker or any entry in the reserved
rollback-anchor namespace exists. The next locked promoter finalizes a complete
valid successor or rolls back a strict partial prepared transaction before
validating a predecessor. Markerless or unexpected rollback-anchor evidence is
never cleaned automatically and requires reviewed manual recovery. An
unreadable or changed visible marker after commit publication begins retains
locked recovery evidence. For a complete-predecessor transaction, absent
post-commit markers leave rollback-anchor evidence and readers closed. A
partial prepared first-ever promotion with no predecessor has no rollback
anchors and requires explicit reviewed manual recovery. Never remove a
transaction marker or rollback anchor manually.

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
from a reviewed clean Git commit that exactly matches the true remote `PIC`
tip. The source runner and production installer independently require that
exact commit, and the clean true-remote gate is repeated immediately before
each install:

```bash
export PIC_ROOT=/lustre/orion/ast207/proj-shared/dfielding/PIC
export PROJECT_HOME_POLICY_ROOT=/autofs/nccs-svm1_proj/ast207/proj-shared/PIC
export PROJECT_HOME_MIRROR_ROOT=/ccs/proj/ast207/proj-shared/PIC
SOURCE_REPO=/autofs/nccs-svm1_home2/dfielding/athenak-pic
PYTHON=/opt/cray/pe/python/3.11.7/bin/python3
GIT=(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin
  /usr/bin/git -c core.fsmonitor=false -c core.hooksPath=/dev/null
)
FULL_GIT_COMMIT="$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse HEAD)"
require_exact_reviewed_source() {
  local remote_pic_tip source_status
  source_status="$(
    "${GIT[@]}" -C "$SOURCE_REPO" status --porcelain --untracked-files=all
  )"
  test -z "$source_status"
  test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse HEAD)"
  test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse origin/PIC)"
  remote_pic_tip="$(
    "${GIT[@]}" -C "$SOURCE_REPO" ls-remote --exit-code origin refs/heads/PIC |
      /usr/bin/awk '$2 == "refs/heads/PIC" {count += 1; tip = $1} END {if (count != 1) exit 1; print tip}'
  )"
  [[ "$remote_pic_tip" =~ ^[0-9a-f]{40}$ ]]
  test "$FULL_GIT_COMMIT" = "$remote_pic_tip"
}
SOURCE_CONTROL_PLANE=(
  "$PYTHON" -I -B
  "${SOURCE_REPO}/tst/publication/frontier_control_plane/run_control_plane.py"
  --expected-git-commit "$FULL_GIT_COMMIT"
)

require_exact_reviewed_source
"${SOURCE_CONTROL_PLANE[@]}" \
  install_control_plane.py \
  --pic-root "$PIC_ROOT" \
  --expected-git-commit "$FULL_GIT_COMMIT"
require_exact_reviewed_source
"${SOURCE_CONTROL_PLANE[@]}" \
  install_control_plane.py \
  --pic-root "$PROJECT_HOME_POLICY_ROOT" \
  --expected-git-commit "$FULL_GIT_COMMIT"
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

Project Home has two intentionally distinct lexical roles. Fresh storage
evidence, policy mirrors and installed control planes use the physical
`${PROJECT_HOME_POLICY_ROOT}` path. The append-only ledger remains bound to
`${PROJECT_HOME_MIRROR_ROOT}` because retained receipts, genesis anchors and
sealed attestations authenticate that historical `/ccs/proj` spelling.

Promote the reviewed policy through the installed control plane:

```bash
export VERSION='<reviewed-control-plane-digest>'
export CONTROL_PLANE_DIR="${PIC_ROOT}/control_plane/${VERSION}"
CONTROL_PLANE=("$PYTHON" -I "${CONTROL_PLANE_DIR}/run_control_plane.py")

"${CONTROL_PLANE[@]}" promote_active_policy.py \
  --reviewed-policy '<reviewed-successor-policy.json>'
```

Promotion validates the reviewed policy, serializes promoters, creates
mirrored prepared transaction markers and four predecessor rollback anchors,
publishes each read-only Orion and Project Home file, and performs a final
coherent reread of both policy copies and both promotion records. It commits
only after both transaction markers are durably upgraded to committed.
Publication cannot be collectively atomic across the two filesystems, so
readers fail closed while either transaction marker or any reserved rollback-
anchor entry exists. A surviving marker binds the exact successor digest of all
four active anchors, the coherent active policy/promotion generation, and
paired installed controller. Those bindings are revalidated before evidence is
deleted, including recovery from one surviving prepared, mixed, or committed
marker. Markerless or unexpected rollback-anchor evidence is preserved for
reviewed manual recovery. A prepared or mixed marker over the complete valid
successor is finalized forward. Any committed marker over a partial successor
retains all transaction evidence for reviewed manual recovery; only a
prepared-only strict partial reachable publication prefix may roll back. Before
rollback the next locked promoter validates the exact rollback predecessor
semantically. It restores in reverse publication order and revalidates the
canonical predecessor while evidence remains. An unreadable or changed visible
marker after commit publication begins retains locked recovery evidence. For a
complete-predecessor transaction, absent post-commit markers leave rollback-
anchor evidence and readers closed. If the successor contains an authorized
clean-candidate freeze, promotion fully revalidates that digest-bound bundle
before any successor anchor is published for an exact migration or freeze
replacement, before committed-marker publication, after both markers are
committed on the normal path, and during complete-successor recovery; exact
successor-anchor checks bracket the post-publication validations. The same visible transaction
generation must remain visible for each validation. The exact launch-prohibited
generation verifier authenticates its executing installed controller and holds
the stable serialization anchor across its active-state, ledger and candidate
proof. Committed completion removes rollback anchors before marker cleanup. The
successor controller must be the verified
current pair, while the candidate receipt must match the freeze's explicit
`build_profile_control_plane_version`, which may be a paired historical
controller when an exact migration preserves the existing freeze. Exact freeze
replacement additionally requires a candidate built by the successor
controller. Candidate or active-anchor drift at the complete commit boundary
retains locked recovery evidence with readers closed; only an already-strict-
partial publication may roll back. Later direct candidate mutation cannot alter
the authorized digest and remains rejected by downstream candidate validation.
Normal active promotion records use schema v2 and one
unique canonical promotion UUID. The one-use exact reviewed predecessor
migration is the sole path that may consume the exact literal live predecessor
bound by the executing successor controller. The currently authorized
replacement predecessor uses schema v2 and binds its historical executing
common-module digest. A direct edit to the repository policy or either active
copy is not an authorization.
Do not promote the retained historical `readiness/storage_policy.json`
directly. The live historical migration must use the generated one-use
retirement successor sequence below.

## Manual Allocation Accounting

Direct scheduler allocations outside the registered submission path are
deviations. They are never scientific evidence. If a reviewed prerequisite
replay used direct `srun` or direct `sbatch`, publish one exact read-only
authorization JSON file below
`${PIC_ROOT}/policy/manual_accounting_authorizations/` and its byte-identical
Project Home mirror below
`${PROJECT_HOME_POLICY_ROOT}/policy/manual_accounting_authorizations/`, bind
both paths and the SHA-256 digest in the reviewed
storage policy, promote that policy, and reconcile the reviewed allocation IDs
through the paired installed control plane:

```bash
export AUTHORIZATION_ID='<reviewed-authorization-id>'

"${CONTROL_PLANE[@]}" reconcile_manual_frontier_allocations.py \
  --authorization "${PIC_ROOT}/policy/manual_accounting_authorizations/${AUTHORIZATION_ID}.json" \
  --ledger-jsonl "${PIC_ROOT}/ledger/node_hours.jsonl" \
  --ledger-csv "${PIC_ROOT}/ledger/node_hours.csv" \
  --receipts-jsonl "${PIC_ROOT}/ledger/mirror_receipts.jsonl" \
  --mirror-jsonl "${PROJECT_HOME_MIRROR_ROOT}/ledger/node_hours.jsonl"
```

The helper accepts only terminal Slurm allocation records with the authorized
account, partition and QoS, requires an empty trusted Frontier queue, no PIC
pending marker and no active reservation, and appends one explicitly
nonqualifying event per allocation. Those node-hours count against the tracked
budget. Before the first append, it publishes matching read-only
`ledger/pending_manual_accounting.json` markers in Orion and Project Home. An
incomplete marker blocks unrelated ledger writers and policy promotion.
Interrupted retries may repair and complete only the exact authorized terminal
suffix after the marker-bound pre-tranche sequence number, authorized-prefix
length and chain head. The marker also binds the installed control-plane
version and promoted policy hashes. Recovery validates fresh reviewed Slurm
accounting against every scheduler-bound suffix field before publishing missing
mirror bytes. Retries reject older mirror or receipt truncation, always
regenerate the derived CSV index after mirrored publication is coherent, and
clear only the matching marker pair after the full authorized tranche is
durable. A predecessor markerless partial prefix retains its historical
control-plane and policy bindings; only its newly appended suffix is bound to
the active successor.

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
export VERSION='<reviewed-control-plane-digest>'
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
  --expected-git-commit '<full-lowercase-git-commit>'
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

The production profile fixes `ROCM_PATH=/opt/rocm-6.2.4`, the Cray `CC` wrapper,
double precision through `Athena_SINGLE_PRECISION=OFF`, and a minimal build
subprocess environment. It unloads the inactive default
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
    "-DAthena_SINGLE_PRECISION=OFF",
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
  --build-profile-id hip-mpi-release-paper-pic \
  --prepared-artifact-inventory \
    tst/publication/frontier_control_plane/prepared_pic_artifact_inventory.json
```

Regenerate and review `prepared_pic_artifact_inventory.json` before the source
commit is frozen. Freeze creation reads that committed JSON from `source.tar`
and revalidates every listed PIC deck and publication analyzer checksum from
archived bytes.

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
`manifest_path`, `manifest_sha256`, and reviewed
`build_profile_control_plane_version`, followed by another active-policy
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

"${CONTROL_PLANE_DIR}/submit_frontier_job.sh" \
  "$MANIFEST" "${PRE_SUBMIT_WRAPPER_ATTESTATION}"
```

The validator-bound queue snapshot above is intentionally the exact six-field
`%i|%P|%q|%T|%j|%k` artifact. The separately archived operator-isolation
attestation captures the seven-field `%i|%a|%P|%q|%T|%j|%k` snapshot including
the account. These artifacts have different roles and are not interchangeable.

The config must reference the exact authorized Orion root, an immutable Slurm
directive template containing `#SBATCH -A AST207`, `#SBATCH -p batch`, one
allowed QOS, node count, walltime and exactly
`#SBATCH -o ${PIC_ROOT}/logs/slurm/%x.%j.log`. Slurm stdout is not permitted
outside that dedicated PIC log path. The template body is never executed. It
must also reference the executable, input deck, environment profile, queue
snapshot, analysis scripts and timeout-margin artifact. Every config must
declare exactly one `submission_scope`. `registered_science` configs must also
reference the policy-authorized `clean_candidate_manifest`, one exact
`registered_science_authorization_id`, and the freshly sealed `pre_manifest`
operator attestation. The installed wrapper requires a separately sealed
`pre_submit_wrapper` attestation and binds both attestation digests into the
reservation ledger before dispatch. The active policy allowlists each
registered slice separately, binding its campaign, test identity, evidence class,
physical mode, minimum-supported runtime profile, QoS, short-job classification,
node, walltime and attempt ceilings, template, deck, environment profile,
analysis-script set, executable, launch contract and clean-candidate digests.
Reservation creation consumes one attempt, including a later cancellation or
failure. Unknown, exhausted and drifting slices fail before reservation. Every config must declare
`job_script_executable_env=PIC_EXECUTABLE`. The creator promotes a read-only
submission snapshot through a pinned descriptor for its authorized Orion
parent. Every registered artifact directory is exactly
`${PIC_ROOT}/runs/<campaign>/<submission-id>`; reruns use new submission IDs
and new reviewed authorization IDs, and never reuse an existing directory.
Reservation rechecks the original
candidate manifest, archived source, clean Git commit and tree attestation, structured
build profile, exact candidate executable digest, and snapshotted executable
digest. A promoted controller successor may continue to use an already authorized
clean candidate: it validates the frozen build-profile receipt against the exact
policy-bound historical controller version and its immutable inventory on both
Orion and Project Home instead of rewriting that receipt to name the successor.

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
The task-local helper also requires a numeric Slurm rank and numeric
`ROCR_VISIBLE_DEVICES` binding, verifies `libamdhip64`, `libmpi_amd` and
`libmpi_gtl_hsa` through `/usr/bin/ldd` on the pinned executable descriptor,
and emits one trusted GPU-launch preflight line before `exec`. After all actions
and declarative post-actions complete, the trampoline freezes every launch
artifact read-only, publishes a checksummed `artifact_inventory.json`, freezes
the artifact root read-only, rejects empty artifact subtrees, and leaves only
the original empty owner-only `analysis/` directory writable. Artifact-root,
declared launch-directory and analysis-staging creation start below a stable
account serialization anchor outside the replaceable `${PIC_ROOT}` name. Each
new directory is created exclusively, inspected without following aliases and
bound to the opened descriptor before use. The trampoline creates `analysis/`
under an unpredictable staging name, pins its descriptor, renames it into place
and verifies that the published entry still names the pinned directory. It
retains the artifact-root descriptor throughout launch, rechecks the lexical
run path around execution and inventory publication, retains every registered
launch-directory identity through workload execution and final publication,
captures workload-created descendant identities after each action and retains
those observed identities through final publication. It rejects replacement
while freezing, including substitution between pre-open metadata inspection
and descriptor acquisition.
Regular-file freeze applies the same pre-open identity binding before hashing.
The freeze retains every workload-payload descriptor through a closing
namespace and byte sweep, then rechecks the retained generated-inventory
descriptor. These checks reject byte-identical replacement.
The generated immutable inventory retains its exclusive-creation descriptor
through read-only freeze and final publication verification.
POSIX `mkdirat` does not return the created directory descriptor. Registered
launch therefore requires same-account process isolation from directory
creation until the first no-follow open and retained-ancestry binding complete.
Workload-created descendants require the same isolation from their creation
through the post-action recursive-capture handoff. After those observation
boundaries, retained ancestry and identity checks reject later substitution.
Offline F1 analysis retains one no-follow
artifact-root descriptor for inventory load, exact tree-closure validation,
every checksummed read and result publication. It rejects duplicate inventory
keys, noncanonical relative paths, unlisted artifacts, root substitution,
redirected `analysis/` and directory or regular-file replacement between
pre-open metadata inspection and descriptor acquisition. It retains the exact
observed inventory bytes and payload-file identity map through analysis, then
publishes `analysis/analysis.json` once without replacement and retains each
published result descriptor through receipt publication. Invoke the snapshotted analyzer only through
`/opt/cray/pe/python/3.11.7/bin/python3 -I -B`; it publishes a second immutable
`analysis/offline_analysis_receipt.json` that binds that runner, the snapshotted
analyzer and support-module digests, the frozen inventory and the passing
analysis result. Frontier qualification requires all three evidence digests
beneath the ledger-bound run directory, verifies the snapshotted analyzer and
support-module bytes, traverses the lexical PIC root plus `runs/`, `manifests/`
and `snapshot/` below the stable account serialization anchor without following
aliases, pins the pre-submit-manifest ancestry before reading it, retains that
manifest descriptor plus descriptors for the run tree, the inventory, result
and receipt files, their directory ancestry and both source files, binds the
exact inventory digest into the child, executes the analyzer through its
inherited `/proc/self/fd` descriptor, and rejects pre-submit-manifest, run-root,
evidence-file, ancestry or source replacement around no-write recomputation
before accepting the qualification manifest.
The registered v2 analyzers reject Athena stderr unless it is byte-identical to
the reviewed Cray MPICH informational transcript. Additional diagnostics,
warnings or changed settings require a new reviewed authorization ID.
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
the registered job, not a trusted launch action. The first authorized analysis
script is snapshotted as `000-<basename>` for explicit operator invocation.
Additional authorized analysis support modules retain their original basenames
inside the same immutable `snapshot/analysis/` directory so an isolated
analyzer can explicitly load only the adjacent snapshotted hashed helper bytes.

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
After the F0 parser-contract smoke passes, promote `frontier_admission_smoke` to
`closed_after_pass`; that terminal state contains no reusable executable fields
and rejects later admission-smoke reservations. The policy rejects every
registered-science allowlist until this terminal closure has been promoted.

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
compute-node trampoline obtains that authority from a descriptor-pinned,
byte-stable, read-only snapshot because the Frontier compute mount does not
provide the mirrored writer-lock operation. Login-side reads and every mutation
retain the mirrored ledger lock. The trampoline rechecks the scheduled manifest
checksum, reservation attachments, submission ID, immutable snapshots, exact
Slurm-template digest and exact executable digest. It sets `PIC_EXECUTABLE` to
the verified executable snapshot and runs only structured Athena actions
through its fixed `/usr/bin/srun` argv builder. It never opens the snapshotted
Slurm template as a shell script.

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

Before each `registered_science` pre-submit manifest is created and again
immediately before its submission wrapper is invoked, archive and review the
same-account process-isolation attestation defined by
`../readiness/q027_frontier_registered_science_same_account_isolation_attestation_template_2026-05-30.json`.
Use the source-controlled helper to capture the exact `ps`, `squeue`,
pending-marker and mirrored-ledger snapshots listed by that template, review
the hidden staging tree, and seal it only after that review:

```bash
ATTESTATION_HELPER=/ccs/home/dfielding/athenak-pic/tst/publication/capture_frontier_pre_policy_promotion_attestation.py
STAGING="$(
  "$PYTHON" "$ATTESTATION_HELPER" capture \
    --authorization-id <registered-science-authorization-id> \
    --control-plane-version "$VERSION" \
    --phase pre_manifest
)"
PRE_MANIFEST_ATTESTATION_ROOT="$(
  "$PYTHON" "$ATTESTATION_HELPER" seal \
    --staging-dir "$STAGING" \
    --attest-reviewed
)"
PRE_MANIFEST_ATTESTATION="${PRE_MANIFEST_ATTESTATION_ROOT}/attestation.json"
```

Repeat that capture-review-seal sequence with `--phase pre_submit_wrapper`
immediately before invoking `submit_frontier_job.sh`. The helper publishes
each immutable tree beneath
`${PIC_ROOT}/operator_attestations/<timestamp>-<authorization-id>-<phase>/`.
Use `--phase pre_policy_promotion` for the same reviewed boundary immediately
before a registered policy promotion, then bind that sealed artifact into the
installed promoter:

```bash
"${CONTROL_PLANE[@]}" promote_active_policy.py \
  --reviewed-policy "${PIC_ROOT}/policy/reviewed_registered_successor.json" \
  --pre-policy-promotion-attestation "${PRE_POLICY_PROMOTION_ATTESTATION}" \
  --pre-policy-promotion-authorization-id "${POLICY_PROMOTION_AUTHORIZATION_ID}"
```

The helper retains the reviewed capture-time process snapshot and refreshes the
process snapshot during sealing. It binds every retained capture snapshot by
digest and rejects capture-to-seal intervals over fifteen minutes. Do not
submit if another same-account process is authorized to mutate either PIC root
throughout registered launch and publication until the frozen artifact
inventory is durably published.

Keep immutable rejected-manifest chronology separate from current live
preflight. Immediately before policy promotion and each registered-science
submission boundary, verify coherent Orion, mirror-receipt and Project Home
chains, an absent Orion `pending_submission.json`, zero active reservations and
zero live ledger hits for any rejected pre-reservation submission UUID. Archive
the current result separately; later valid reservations may advance live
ledger counts without changing the historical incident record.

## Q011 Committed Repair-Validation Boundary

Live mutation and pressure publication require a successful clean committed
repair-validation worker first. Submit it only after the reviewed repair is
committed and pushed. This is a non-science validation worker. Before its
`sbatch`, inspect all same-user processes and Frontier jobs and confirm that no
same-user command can launch science or mutate either PIC root concurrently:

```bash
(
set -euo pipefail
test "${PIC_ROOT:-}" = /lustre/orion/ast207/proj-shared/dfielding/PIC
SOURCE_REPO=/autofs/nccs-svm1_home2/dfielding/athenak-pic
SLURM_ENV=(/usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin SLURM_CLUSTERS=frontier)
GIT=(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin
  /usr/bin/git -c core.fsmonitor=false -c core.hooksPath=/dev/null
)
SOURCE_STATUS="$("${GIT[@]}" -C "$SOURCE_REPO" status --porcelain --untracked-files=all)"
test -z "$SOURCE_STATUS"
FULL_GIT_COMMIT="$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse HEAD)"
test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse origin/PIC)"
REMOTE_PIC_TIP="$(
  "${GIT[@]}" -C "$SOURCE_REPO" ls-remote --exit-code origin refs/heads/PIC |
    /usr/bin/awk '$2 == "refs/heads/PIC" {count += 1; tip = $1} END {if (count != 1) exit 1; print tip}'
)"
[[ "$REMOTE_PIC_TIP" =~ ^[0-9a-f]{40}$ ]]
test "$FULL_GIT_COMMIT" = "$REMOTE_PIC_TIP"
/usr/bin/ps -u "$(/usr/bin/id -u)" -ww -o pid=,ppid=,lstart=,args=
"${SLURM_ENV[@]}" /usr/bin/squeue --clusters=frontier -u "$(/usr/bin/id -un)" \
  -h -o '%i|%j|%a|%q|%T|%Z|%o'
REVIEWED_NO_SCIENCE_OR_SAME_USER_PIC_MUTATOR='<yes only after reviewing both snapshots>'
test "$REVIEWED_NO_SCIENCE_OR_SAME_USER_PIC_MUTATOR" = yes
cd "$SOURCE_REPO"
REPAIR_VALIDATION_JOB_TOKEN="$("${SLURM_ENV[@]}" /usr/bin/sbatch --parsable --export=NIL \
  "${SOURCE_REPO}/tst/publication/frontier_q011_section54_pressure_gate_validation_job.sh" \
  "$FULL_GIT_COMMIT")"
printf 'repair_validation_job_token=%q\n' "$REPAIR_VALIDATION_JOB_TOKEN"
REPAIR_VALIDATION_JOB_ID="${REPAIR_VALIDATION_JOB_TOKEN%;frontier}"
[[ "$REPAIR_VALIDATION_JOB_ID" =~ ^[0-9]+$ ]]
test "$REPAIR_VALIDATION_JOB_TOKEN" = "$REPAIR_VALIDATION_JOB_ID" ||
  test "$REPAIR_VALIDATION_JOB_TOKEN" = "${REPAIR_VALIDATION_JOB_ID};frontier"
printf 'repair_validation_commit=%s\n' "$FULL_GIT_COMMIT"
printf 'repair_validation_job_id=%s\n' "$REPAIR_VALIDATION_JOB_ID"
)
```

After that worker terminates, fill the printed values into this read-only
checkpoint verifier. Preserve both values for the migration and publication
procedures below:

```bash
(
set -euo pipefail
test "${PIC_ROOT:-}" = /lustre/orion/ast207/proj-shared/dfielding/PIC
SOURCE_REPO=/autofs/nccs-svm1_home2/dfielding/athenak-pic
SLURM_ENV=(/usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin SLURM_CLUSTERS=frontier)
GIT=(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin
  /usr/bin/git -c core.fsmonitor=false -c core.hooksPath=/dev/null
)
REPAIR_VALIDATION_COMMIT='<repair_validation_commit printed by submission>'
REPAIR_VALIDATION_JOB_ID='<repair_validation_job_id printed by submission>'
[[ "$REPAIR_VALIDATION_COMMIT" =~ ^[0-9a-f]{40}$ ]]
[[ "$REPAIR_VALIDATION_JOB_ID" =~ ^[0-9]+$ ]]
SOURCE_STATUS="$("${GIT[@]}" -C "$SOURCE_REPO" status --porcelain --untracked-files=all)"
test -z "$SOURCE_STATUS"
test "$REPAIR_VALIDATION_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse HEAD)"
test "$REPAIR_VALIDATION_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse origin/PIC)"
REMOTE_PIC_TIP="$(
  "${GIT[@]}" -C "$SOURCE_REPO" ls-remote --exit-code origin refs/heads/PIC |
    /usr/bin/awk '$2 == "refs/heads/PIC" {count += 1; tip = $1} END {if (count != 1) exit 1; print tip}'
)"
test "$REPAIR_VALIDATION_COMMIT" = "$REMOTE_PIC_TIP"
REPAIR_VALIDATION_RECORD="$(
  "${SLURM_ENV[@]}" /usr/bin/sacct -X --clusters=frontier -n \
    -j "$REPAIR_VALIDATION_JOB_ID" \
    --format=JobIDRaw,JobName,Account,Partition,QOS,State,ExitCode,WorkDir --parsable2 |
    /usr/bin/awk -F'|' -v job="$REPAIR_VALIDATION_JOB_ID" '
      $1 == job {count += 1; record = $1 FS $2 FS $3 FS $4 FS $5 FS $6 FS $7 FS $8}
      END {if (count != 1) exit 1; print record}
    '
)"
test "$REPAIR_VALIDATION_RECORD" = \
  "${REPAIR_VALIDATION_JOB_ID}|pic-q011-pressure-gate-validate|ast207|batch|debug|COMPLETED|0:0|${SOURCE_REPO}"
REPAIR_VALIDATION_LOG="${PIC_ROOT}/logs/slurm/pic-q011-pressure-gate-validate.${REPAIR_VALIDATION_JOB_ID}.log"
test -r "$REPAIR_VALIDATION_LOG"
REPAIR_VALIDATION_LOG_COMMIT="$(
  /usr/bin/sed -n 's/^source_commit=//p' "$REPAIR_VALIDATION_LOG"
)"
test "$REPAIR_VALIDATION_LOG_COMMIT" = "$REPAIR_VALIDATION_COMMIT"
test "$(printf '%s\n' "$REPAIR_VALIDATION_LOG_COMMIT" | /usr/bin/wc -l)" -eq 1
)
```

Treat submission of this worker as a durable checkpoint. If a session stops
after `sbatch`, recover the one printed job ID and verify that exact job; do not
silently submit a replacement.

Frontier may return either `<job-id>` or `<job-id>;frontier` for
`sbatch --parsable`. Every live snippet below prints the raw token with shell
escaping immediately after submission, strips only that exact optional cluster
suffix, then requires the normalized job ID to be numeric. Reject any other
scheduler token, but retain its escaped diagnostic for recovery.

## Q011 Pressure-Pilot Serial Boundary

The four Q011 Section 5.4 pressure-sensitivity pilots are one engineering
calibration tranche, but they were four separate registered-science slices.
The v1 `ps_p0_1p00` slice failed closed as scheduler job `4754211` when the
registered HIP executable attempted a strided host-to-device particle subview
copy. Its terminal ledger event is retained as immutable chronology. The four
rebuilt v2 slices completed and are consumed historical evidence. Do not
relaunch them.

The current live policy contains authenticated mirrored schema-v2
storage-preflight evidence produced by the prior repaired controller. Capture
one successor preflight binding from clean tracked source, persist the emitted
fragment as a read-only review artifact under `${PIC_ROOT}/policy`, then
install the paired successor controller. Materialize the exact
predecessor-migration successor only after both installs print the same digest.
The one-use migration flag accepts only the exact reviewed live policy,
promotion, controller, preflight binding, and complete schema-v2 evidence
source-authentication tuple. The successor must preserve every policy field
except the newer controller version and a different, strictly newer
strict-current preflight binding; this preserves the authorized clean-candidate
freeze and empty registered-slice allowlist. This sequence authorizes no
science. The only permitted `sbatch` below is the exact non-science
build/freeze worker. Never invoke `pilot-policy-successor`, `policy-slices`,
`pre-submit-config`, `submit_frontier_job.sh`,
`launch_with_frontier_profile.sh`, or the historical launch block while
executing this sequence. Every reviewed `ps` and `squeue` gate must confirm
that no same-user process or job can launch science or mutate either PIC root.
Every source-only control-plane invocation below binds the same
`FULL_GIT_COMMIT` independently at both the source runner and source-only
entrypoint. Each live Q011 policy materialization repeats the clean true-remote
gate, creates a fresh detached worktree snapshot at that exact commit, and
passes the materializer's own required commit binding. The final read-only
publication verification instead executes from a fresh `git archive` snapshot
of that reviewed commit, never from mutable checkout files.
The active predecessor and successor generation verifiers both accept schema
v2. Independently bracket migration with the exact mirrored anchors, paired
active predecessor controller, and explicitly historical preserved-candidate
build controller as shown below:

```bash
(
set -euo pipefail
test "${PIC_ROOT:-}" = /lustre/orion/ast207/proj-shared/dfielding/PIC
test "${PYTHON:-}" = /opt/cray/pe/python/3.11.7/bin/python3
PROJECT_HOME_POLICY_ROOT=/autofs/nccs-svm1_proj/ast207/proj-shared/PIC
PROJECT_HOME_MIRROR_ROOT=/ccs/proj/ast207/proj-shared/PIC
SOURCE_REPO=/autofs/nccs-svm1_home2/dfielding/athenak-pic
SLURM_ENV=(/usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin SLURM_CLUSTERS=frontier)
GIT=(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin
  /usr/bin/git -c core.fsmonitor=false -c core.hooksPath=/dev/null
)
ACTIVE_PREDECESSOR_CONTROL_PLANE_VERSION=b56d96b40f2c666b9fa5b421fac589d6a4d6500a716f20b479c354a3d239cb48
ACTIVE_PREDECESSOR_POLICY_SHA256=48e74f3151b51ba84b72e3214ef4a66c03441c745912b833395f98f6f60a0bd1
ACTIVE_PREDECESSOR_PROMOTION_SHA256=4ecb45bb399cee750ce69198286afc9fb903dde92cffc7507600b0a53f1c547f
PRESERVED_CLEAN_CANDIDATE_MANIFEST="${PIC_ROOT}/clean_candidates/98a372c9-2ea0-47e6-ad34-e66343e7eea1/clean_candidate_manifest.json"
PRESERVED_CLEAN_CANDIDATE_MANIFEST_SHA256=dee6be45657e99ec477eec513c45be4ba6ac43b4fd5deb5655750f51c668e42f
PRESERVED_CLEAN_CANDIDATE_BUILD_PROFILE_CONTROL_PLANE_VERSION=821d185856722bd0178acb9427f78ac82671a4b6670779ec8400fbac54c6d721
ACTIVE_PREDECESSOR_ORION_CONTROL_PLANE="${PIC_ROOT}/control_plane/${ACTIVE_PREDECESSOR_CONTROL_PLANE_VERSION}"
ACTIVE_PREDECESSOR_PROJECT_HOME_CONTROL_PLANE="${PROJECT_HOME_POLICY_ROOT}/control_plane/${ACTIVE_PREDECESSOR_CONTROL_PLANE_VERSION}"
ACTIVE_PREDECESSOR_CONTROL_PLANE=("$PYTHON" -I -B "${ACTIVE_PREDECESSOR_ORION_CONTROL_PLANE}/run_control_plane.py")
require_no_policy_recovery_entries() {
  local policy_parent recovery_entry
  for policy_parent in "${PIC_ROOT}/policy" "${PROJECT_HOME_POLICY_ROOT}/policy"; do
    test -d "$policy_parent"
    test ! -L "$policy_parent"
    recovery_entry="$(
      /usr/bin/find "$policy_parent" -mindepth 1 -maxdepth 1 \
        \( -name '.active_promotion_transaction.json' \
        -o -name '.storage_policy.json.transaction-rollback-*' \
        -o -name '.active_promotion.json.transaction-rollback-*' \
        -o -name '.storage_policy.json.tmp-*' \
        -o -name '.active_promotion.json.tmp-*' \
        -o -name '.storage_policy.json.rollback-*' \
        -o -name '.active_promotion.json.rollback-*' \) \
        -print -quit
    )"
    test -z "$recovery_entry"
  done
}
require_no_staging_entries() {
  local staging_parent staging_entry
  for staging_parent in \
    "${PIC_ROOT}/policy/storage_preflight_bindings" \
    "${PIC_ROOT}/policy/storage_preflight_evidence" \
    "${PROJECT_HOME_POLICY_ROOT}/policy/storage_preflight_evidence" \
    "${PIC_ROOT}/control_plane" \
    "${PROJECT_HOME_POLICY_ROOT}/control_plane" \
    "${PIC_ROOT}/clean_candidates"; do
    test -d "$staging_parent"
    test ! -L "$staging_parent"
    staging_entry="$(
      /usr/bin/find "$staging_parent" -mindepth 1 -maxdepth 1 \
        \( -name '.staging.*' -o -name '.tmp-*' -o -name '*.recovery-staging' \) \
        -print -quit
    )"
    test -z "$staging_entry"
  done
  for staging_parent in \
    "${PIC_ROOT}/policy/storage_preflight_evidence" \
    "${PROJECT_HOME_POLICY_ROOT}/policy/storage_preflight_evidence"; do
    staging_entry="$(
      /usr/bin/find "$staging_parent" -mindepth 1 -maxdepth 1 \
        \( ! -type f -o ! -name '*.json' -o -perm /222 -o ! -links 1 \) \
        -print -quit
    )"
    test -z "$staging_entry"
  done
  for staging_parent in "$PIC_ROOT" "$PROJECT_HOME_POLICY_ROOT"; do
    staging_entry="$(
      /usr/bin/find "$staging_parent" -mindepth 1 -maxdepth 1 \
        -name '.pic-storage-preflight-*' -print -quit
    )"
    test -z "$staging_entry"
  done
  local orion_evidence project_home_evidence evidence_name
  orion_evidence="$(
    /usr/bin/find "${PIC_ROOT}/policy/storage_preflight_evidence" \
      -mindepth 1 -maxdepth 1 -type f -name '*.json' -printf '%f\n' |
      /usr/bin/sort
  )"
  project_home_evidence="$(
    /usr/bin/find "${PROJECT_HOME_POLICY_ROOT}/policy/storage_preflight_evidence" \
      -mindepth 1 -maxdepth 1 -type f -name '*.json' -printf '%f\n' |
      /usr/bin/sort
  )"
  test -n "$orion_evidence"
  test "$orion_evidence" = "$project_home_evidence"
  while IFS= read -r evidence_name; do
    /usr/bin/cmp -s \
      "${PIC_ROOT}/policy/storage_preflight_evidence/${evidence_name}" \
      "${PROJECT_HOME_POLICY_ROOT}/policy/storage_preflight_evidence/${evidence_name}"
  done <<< "$orion_evidence"
}
require_exact_reviewed_source() {
  local remote_pic_tip source_status
  source_status="$(
    "${GIT[@]}" -C "$SOURCE_REPO" status --porcelain --untracked-files=all
  )"
  test -z "$source_status"
  test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse HEAD)"
  test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse origin/PIC)"
  remote_pic_tip="$(
    "${GIT[@]}" -C "$SOURCE_REPO" ls-remote --exit-code origin refs/heads/PIC |
      /usr/bin/awk '$2 == "refs/heads/PIC" {count += 1; tip = $1} END {if (count != 1) exit 1; print tip}'
  )"
  [[ "$remote_pic_tip" =~ ^[0-9a-f]{40}$ ]]
  test "$FULL_GIT_COMMIT" = "$remote_pic_tip"
}
run_reviewed_q011_materializer() {
  local snapshot snapshot_root status
  require_exact_reviewed_source
  snapshot_root="$(/usr/bin/mktemp -d /tmp/athenak-pic-q011-reviewed.XXXXXX)"
  snapshot="${snapshot_root}/source"
  status=0
  "${GIT[@]}" -C "$SOURCE_REPO" worktree add --detach "$snapshot" \
    "$FULL_GIT_COMMIT" || status=$?
  if test "$status" -eq 0; then
    if test -f "$snapshot/tst/publication/q011_section54_pressure_pilot_execution.py" &&
      test ! -L "$snapshot/tst/publication/q011_section54_pressure_pilot_execution.py"; then
      "$PYTHON" -I -B \
        "$snapshot/tst/publication/q011_section54_pressure_pilot_execution.py" \
        "$@" || status=$?
    else
      status=1
    fi
  fi
  if test -e "$snapshot"; then
    "${GIT[@]}" -C "$SOURCE_REPO" worktree remove --force "$snapshot" || status=1
  fi
  /usr/bin/rm -rf -- "$snapshot_root" || status=1
  return "$status"
}
verify_exact_active_predecessor_launch_prohibited_generation() {
  local remote_pic_tip
  remote_pic_tip="$(
    "${GIT[@]}" -C "$SOURCE_REPO" ls-remote --exit-code origin refs/heads/PIC |
      /usr/bin/awk '$2 == "refs/heads/PIC" {count += 1; tip = $1} END {if (count != 1) exit 1; print tip}'
  )"
  test "$FULL_GIT_COMMIT" = "$remote_pic_tip"
  test -d "$ACTIVE_PREDECESSOR_ORION_CONTROL_PLANE"
  test -d "$ACTIVE_PREDECESSOR_PROJECT_HOME_CONTROL_PLANE"
  /usr/bin/cmp -s \
    "${ACTIVE_PREDECESSOR_ORION_CONTROL_PLANE}/inventory.json" \
    "${ACTIVE_PREDECESSOR_PROJECT_HOME_CONTROL_PLANE}/inventory.json"
  test "$(
    "${ACTIVE_PREDECESSOR_CONTROL_PLANE[@]}" \
      validate_and_reserve_frontier_job.py verify-control-plane
  )" = "$ACTIVE_PREDECESSOR_CONTROL_PLANE_VERSION"
  /usr/bin/cmp -s \
    "${PIC_ROOT}/policy/storage_policy.json" \
    "${PROJECT_HOME_POLICY_ROOT}/policy/storage_policy.json"
  /usr/bin/cmp -s \
    "${PIC_ROOT}/policy/active_promotion.json" \
    "${PROJECT_HOME_POLICY_ROOT}/policy/active_promotion.json"
  test "$ACTIVE_PREDECESSOR_POLICY_SHA256" = "$(
    /usr/bin/sha256sum "${PIC_ROOT}/policy/storage_policy.json" |
      /usr/bin/awk '{print $1}'
  )"
  test "$ACTIVE_PREDECESSOR_PROMOTION_SHA256" = "$(
    /usr/bin/sha256sum "${PIC_ROOT}/policy/active_promotion.json" |
      /usr/bin/awk '{print $1}'
  )"
  "$PYTHON" -I -B - \
    "$ACTIVE_PREDECESSOR_ORION_CONTROL_PLANE" \
    "$ACTIVE_PREDECESSOR_PROJECT_HOME_CONTROL_PLANE" \
    "$PIC_ROOT" \
    "$PROJECT_HOME_POLICY_ROOT" \
    "$ACTIVE_PREDECESSOR_CONTROL_PLANE_VERSION" \
    "$ACTIVE_PREDECESSOR_POLICY_SHA256" \
    "$ACTIVE_PREDECESSOR_PROMOTION_SHA256" \
    "$PRESERVED_CLEAN_CANDIDATE_MANIFEST" \
    "$PRESERVED_CLEAN_CANDIDATE_MANIFEST_SHA256" \
    "$PRESERVED_CLEAN_CANDIDATE_BUILD_PROFILE_CONTROL_PLANE_VERSION" <<'PY'
import sys
from pathlib import Path

old_orion = Path(sys.argv[1])
old_project_home = Path(sys.argv[2])
pic_root = Path(sys.argv[3])
project_home_root = Path(sys.argv[4])
version, policy_sha256, promotion_sha256 = sys.argv[5:8]
preserved_manifest, preserved_manifest_sha256, preserved_build_controller = sys.argv[8:11]
sys.path.insert(0, str(old_orion))

from control_plane_common import require_storage_policy_unlock_snapshot
from control_plane_common import verify_installed_control_plane
from control_plane_common import project_home_ledger_root
from ledger import require_no_incomplete_manual_accounting_marker
from promote_active_policy import _require_no_outstanding_submissions

orion_inventory = verify_installed_control_plane(
    old_orion,
    authorized_pic_root=pic_root,
)
project_home_inventory = verify_installed_control_plane(
    old_project_home,
    authorized_pic_root=project_home_root,
)
assert orion_inventory == project_home_inventory
assert orion_inventory["version"] == version

policy, snapshot = require_storage_policy_unlock_snapshot(
    control_plane_version=version,
    authorized_pic_root=pic_root,
    authorized_project_home_root=project_home_root,
    allow_pending_genesis=True,
)
assert snapshot == {
    "active_policy_sha256": policy_sha256,
    "active_promotion_sha256": promotion_sha256,
}
assert policy["registered_science_slices"] == []
assert policy["science_submission_freeze"] == {
    "build_profile_control_plane_version": preserved_build_controller,
    "manifest_path": preserved_manifest,
    "manifest_sha256": preserved_manifest_sha256,
    "status": "authorized",
}
require_no_incomplete_manual_accounting_marker(
    pic_root / "ledger" / "node_hours.jsonl",
    project_home_ledger_root(project_home_root) / "ledger" / "node_hours.jsonl",
)
_require_no_outstanding_submissions(policy, pic_root, project_home_root)
PY
  "${ACTIVE_PREDECESSOR_CONTROL_PLANE[@]}" revalidate_clean_candidate.py \
    --manifest "$PRESERVED_CLEAN_CANDIDATE_MANIFEST" \
    --expected-manifest-sha256 "$PRESERVED_CLEAN_CANDIDATE_MANIFEST_SHA256"
}
SOURCE_STATUS="$(
  "${GIT[@]}" -C "$SOURCE_REPO" status --porcelain --untracked-files=all
)"
test -z "$SOURCE_STATUS"
FULL_GIT_COMMIT="$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse HEAD)"
test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse origin/PIC)"
REMOTE_PIC_TIP="$(
  "${GIT[@]}" -C "$SOURCE_REPO" ls-remote --exit-code origin refs/heads/PIC |
    /usr/bin/awk '$2 == "refs/heads/PIC" {count += 1; tip = $1} END {if (count != 1) exit 1; print tip}'
)"
[[ "$REMOTE_PIC_TIP" =~ ^[0-9a-f]{40}$ ]]
test "$FULL_GIT_COMMIT" = "$REMOTE_PIC_TIP"
SOURCE_CONTROL_PLANE=(
  "$PYTHON" -I -B
  "${SOURCE_REPO}/tst/publication/frontier_control_plane/run_control_plane.py"
  --expected-git-commit "$FULL_GIT_COMMIT"
)
REPAIR_VALIDATION_COMMIT='<repair_validation_commit from the completed checkpoint>'
REPAIR_VALIDATION_JOB_ID='<repair_validation_job_id from the completed checkpoint>'
[[ "$REPAIR_VALIDATION_COMMIT" =~ ^[0-9a-f]{40}$ ]]
[[ "$REPAIR_VALIDATION_JOB_ID" =~ ^[0-9]+$ ]]
test "$REPAIR_VALIDATION_COMMIT" = "$FULL_GIT_COMMIT"
REPAIR_VALIDATION_RECORD="$(
  "${SLURM_ENV[@]}" /usr/bin/sacct -X --clusters=frontier -n \
    -j "$REPAIR_VALIDATION_JOB_ID" \
    --format=JobIDRaw,JobName,Account,Partition,QOS,State,ExitCode,WorkDir --parsable2 |
    /usr/bin/awk -F'|' -v job="$REPAIR_VALIDATION_JOB_ID" '
      $1 == job {count += 1; record = $1 FS $2 FS $3 FS $4 FS $5 FS $6 FS $7 FS $8}
      END {if (count != 1) exit 1; print record}
    '
)"
test "$REPAIR_VALIDATION_RECORD" = \
  "${REPAIR_VALIDATION_JOB_ID}|pic-q011-pressure-gate-validate|ast207|batch|debug|COMPLETED|0:0|${SOURCE_REPO}"
REPAIR_VALIDATION_LOG="${PIC_ROOT}/logs/slurm/pic-q011-pressure-gate-validate.${REPAIR_VALIDATION_JOB_ID}.log"
test -r "$REPAIR_VALIDATION_LOG"
REPAIR_VALIDATION_LOG_COMMIT="$(
  /usr/bin/sed -n 's/^source_commit=//p' "$REPAIR_VALIDATION_LOG"
)"
test "$REPAIR_VALIDATION_LOG_COMMIT" = "$FULL_GIT_COMMIT"
test "$(printf '%s\n' "$REPAIR_VALIDATION_LOG_COMMIT" | /usr/bin/wc -l)" -eq 1
require_no_policy_recovery_entries
require_no_staging_entries
/usr/bin/ps -u "$(/usr/bin/id -u)" -ww -o pid=,ppid=,lstart=,args=
"${SLURM_ENV[@]}" /usr/bin/squeue --clusters=frontier -u "$(/usr/bin/id -un)" \
  -h -o '%i|%j|%a|%q|%T|%Z|%o'
REVIEWED_NO_SCIENCE_OR_SAME_USER_PIC_MUTATOR='<yes only after reviewing both snapshots>'
test "$REVIEWED_NO_SCIENCE_OR_SAME_USER_PIC_MUTATOR" = yes
verify_exact_active_predecessor_launch_prohibited_generation
STORAGE_PREFLIGHT_BINDING_DIR="${PIC_ROOT}/policy/storage_preflight_bindings"
test -d "$STORAGE_PREFLIGHT_BINDING_DIR"
test ! -L "$STORAGE_PREFLIGHT_BINDING_DIR"
test "$(
  /usr/bin/stat -Lc '%F|%a|%u' "$STORAGE_PREFLIGHT_BINDING_DIR"
)" = "directory|2700|$(/usr/bin/id -u)"
STORAGE_PREFLIGHT_BINDING_DIR_ID="$(
  /usr/bin/stat -Lc '%d:%i' "$STORAGE_PREFLIGHT_BINDING_DIR"
)"
require_exact_reviewed_source
STORAGE_PREFLIGHT_FRAGMENT="$(
  "${SOURCE_CONTROL_PLANE[@]}" capture_storage_preflight_evidence.py \
    --expected-git-commit "$FULL_GIT_COMMIT"
)"
test "$(
  /usr/bin/stat -Lc '%d:%i' "$STORAGE_PREFLIGHT_BINDING_DIR"
)" = "$STORAGE_PREFLIGHT_BINDING_DIR_ID"
PROBE_ID="$(
  printf '%s\n' "$STORAGE_PREFLIGHT_FRAGMENT" |
    "$PYTHON" -I -B -c \
      'import json, sys; print(json.load(sys.stdin)["storage_preflight_evidence"]["probe_id"])'
)"
[[ "$PROBE_ID" =~ ^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$ ]]
STORAGE_PREFLIGHT_BINDING="${STORAGE_PREFLIGHT_BINDING_DIR}/${PROBE_ID}.json"
test ! -e "$STORAGE_PREFLIGHT_BINDING"
STORAGE_PREFLIGHT_STAGING="$(
  /usr/bin/mktemp -p "$STORAGE_PREFLIGHT_BINDING_DIR" .staging.XXXXXX
)"
printf '%s\n' "$STORAGE_PREFLIGHT_FRAGMENT" > "$STORAGE_PREFLIGHT_STAGING"
/usr/bin/chmod 0400 "$STORAGE_PREFLIGHT_STAGING"
/usr/bin/ln -T -- "$STORAGE_PREFLIGHT_STAGING" "$STORAGE_PREFLIGHT_BINDING"
/usr/bin/rm "$STORAGE_PREFLIGHT_STAGING"
test "$(
  /usr/bin/stat -Lc '%d:%i' "$STORAGE_PREFLIGHT_BINDING_DIR"
)" = "$STORAGE_PREFLIGHT_BINDING_DIR_ID"
"$PYTHON" -I -B -c \
  'import os, sys; fd = os.open(sys.argv[1], os.O_RDONLY | os.O_NOFOLLOW); os.fsync(fd); os.close(fd); fd = os.open(os.path.dirname(sys.argv[1]), os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW); os.fsync(fd); os.close(fd)' \
  "$STORAGE_PREFLIGHT_BINDING"

require_no_policy_recovery_entries
require_no_staging_entries
verify_exact_active_predecessor_launch_prohibited_generation
require_exact_reviewed_source
ORION_CONTROL_PLANE="$(
  "${SOURCE_CONTROL_PLANE[@]}" install_control_plane.py \
    --pic-root "$PIC_ROOT" \
    --expected-git-commit "$FULL_GIT_COMMIT"
)"
printf 'one_sided_orion_control_plane=%s\n' "$ORION_CONTROL_PLANE"
require_no_policy_recovery_entries
require_no_staging_entries
verify_exact_active_predecessor_launch_prohibited_generation
require_exact_reviewed_source
PROJECT_HOME_CONTROL_PLANE="$(
  "${SOURCE_CONTROL_PLANE[@]}" install_control_plane.py \
    --pic-root "$PROJECT_HOME_POLICY_ROOT" \
    --expected-git-commit "$FULL_GIT_COMMIT"
)"
printf 'one_sided_project_home_control_plane=%s\n' "$PROJECT_HOME_CONTROL_PLANE"
VERSION="${ORION_CONTROL_PLANE##*/}"
[[ "$VERSION" =~ ^[0-9a-f]{64}$ ]]
test "$ORION_CONTROL_PLANE" = "${PIC_ROOT}/control_plane/${VERSION}"
test "$PROJECT_HOME_CONTROL_PLANE" = \
  "${PROJECT_HOME_POLICY_ROOT}/control_plane/${VERSION}"
test "${PROJECT_HOME_CONTROL_PLANE##*/}" = "$VERSION"
/usr/bin/cmp -s \
  "${ORION_CONTROL_PLANE}/inventory.json" \
  "${PROJECT_HOME_CONTROL_PLANE}/inventory.json"
require_no_staging_entries
CONTROL_PLANE=("$PYTHON" -I -B "${ORION_CONTROL_PLANE}/run_control_plane.py")
POLICY_SUFFIX="${VERSION:0:8}_${PROBE_ID}"
MIGRATION_POLICY="${PIC_ROOT}/policy/reviewed_exact_preflight_predecessor_successor_${POLICY_SUFFIX}.json"
CANDIDATE_ONLY_POLICY="${PIC_ROOT}/policy/reviewed_candidate_only_successor_${POLICY_SUFFIX}.json"
test ! -e "$MIGRATION_POLICY"
test ! -e "$CANDIDATE_ONLY_POLICY"
require_no_policy_recovery_entries
require_no_staging_entries
verify_exact_active_predecessor_launch_prohibited_generation
/usr/bin/ps -u "$(/usr/bin/id -u)" -ww -o pid=,ppid=,lstart=,args=
"${SLURM_ENV[@]}" /usr/bin/squeue --clusters=frontier -u "$(/usr/bin/id -un)" \
  -h -o '%i|%j|%a|%q|%T|%Z|%o'
REVIEWED_NO_SCIENCE_OR_SAME_USER_PIC_MUTATOR_BEFORE_MIGRATION='<yes only after reviewing both snapshots>'
test "$REVIEWED_NO_SCIENCE_OR_SAME_USER_PIC_MUTATOR_BEFORE_MIGRATION" = yes

run_reviewed_q011_materializer \
  exact-reviewed-preflight-predecessor-policy-successor \
  --expected-git-commit "$FULL_GIT_COMMIT" \
  --baseline-policy "${PIC_ROOT}/policy/storage_policy.json" \
  --control-plane-version "$VERSION" \
  --storage-preflight-binding "$STORAGE_PREFLIGHT_BINDING" \
  --output "$MIGRATION_POLICY"

"${CONTROL_PLANE[@]}" promote_active_policy.py \
  --reviewed-policy "$MIGRATION_POLICY" \
  --migrate-exact-reviewed-storage-preflight-predecessor

EXPECTED_ACTIVE_POLICY_SHA256="$(
  /usr/bin/sha256sum "${PIC_ROOT}/policy/storage_policy.json" |
    /usr/bin/awk '{print $1}'
)"
EXPECTED_ACTIVE_PROMOTION_SHA256="$(
  /usr/bin/sha256sum "${PIC_ROOT}/policy/active_promotion.json" |
    /usr/bin/awk '{print $1}'
)"
[[ "$EXPECTED_ACTIVE_POLICY_SHA256" =~ ^[0-9a-f]{64}$ ]]
[[ "$EXPECTED_ACTIVE_PROMOTION_SHA256" =~ ^[0-9a-f]{64}$ ]]
require_no_policy_recovery_entries
require_no_staging_entries
"${CONTROL_PLANE[@]}" promote_active_policy.py \
  --verify-active-launch-prohibited-generation \
  --expected-control-plane-version "$VERSION" \
  --expected-active-policy-sha256 "$EXPECTED_ACTIVE_POLICY_SHA256" \
  --expected-active-promotion-sha256 "$EXPECTED_ACTIVE_PROMOTION_SHA256" \
  --expected-authorized-freeze-manifest "$PRESERVED_CLEAN_CANDIDATE_MANIFEST" \
  --expected-authorized-freeze-manifest-sha256 "$PRESERVED_CLEAN_CANDIDATE_MANIFEST_SHA256" \
  --expected-authorized-freeze-build-controller \
    "$PRESERVED_CLEAN_CANDIDATE_BUILD_PROFILE_CONTROL_PLANE_VERSION"
printf 'full_git_commit=%s\n' "$FULL_GIT_COMMIT"
printf 'probe_id=%s\n' "$PROBE_ID"
printf 'control_plane_version=%s\n' "$VERSION"
printf 'expected_active_policy_sha256=%s\n' "$EXPECTED_ACTIVE_POLICY_SHA256"
printf 'expected_active_promotion_sha256=%s\n' "$EXPECTED_ACTIVE_PROMOTION_SHA256"
printf 'preserved_clean_candidate_manifest=%s\n' "$PRESERVED_CLEAN_CANDIDATE_MANIFEST"
printf 'preserved_clean_candidate_manifest_sha256=%s\n' "$PRESERVED_CLEAN_CANDIDATE_MANIFEST_SHA256"
printf 'preserved_clean_candidate_build_profile_control_plane_version=%s\n' \
  "$PRESERVED_CLEAN_CANDIDATE_BUILD_PROFILE_CONTROL_PLANE_VERSION"

/usr/bin/ps -u "$(/usr/bin/id -u)" -ww -o pid=,ppid=,lstart=,args=
"${SLURM_ENV[@]}" /usr/bin/squeue --clusters=frontier -u "$(/usr/bin/id -un)" \
  -h -o '%i|%j|%a|%q|%T|%Z|%o'
REVIEWED_NO_SCIENCE_OR_SAME_USER_PIC_MUTATOR_BEFORE_BUILD_FREEZE='<yes only after reviewing both snapshots>'
test "$REVIEWED_NO_SCIENCE_OR_SAME_USER_PIC_MUTATOR_BEFORE_BUILD_FREEZE" = yes
cd "$SOURCE_REPO"
BUILD_FREEZE_JOB_TOKEN="$("${SLURM_ENV[@]}" /usr/bin/sbatch --parsable --export=NIL \
  "${SOURCE_REPO}/tst/publication/frontier_q011_clean_candidate_build_freeze_job.sh" \
  "$FULL_GIT_COMMIT" "$VERSION" \
  "$EXPECTED_ACTIVE_POLICY_SHA256" "$EXPECTED_ACTIVE_PROMOTION_SHA256" \
  "$PRESERVED_CLEAN_CANDIDATE_MANIFEST" \
  "$PRESERVED_CLEAN_CANDIDATE_MANIFEST_SHA256" \
  "$PRESERVED_CLEAN_CANDIDATE_BUILD_PROFILE_CONTROL_PLANE_VERSION")"
printf 'build_freeze_job_token=%q\n' "$BUILD_FREEZE_JOB_TOKEN"
BUILD_FREEZE_JOB_ID="${BUILD_FREEZE_JOB_TOKEN%;frontier}"
[[ "$BUILD_FREEZE_JOB_ID" =~ ^[0-9]+$ ]]
test "$BUILD_FREEZE_JOB_TOKEN" = "$BUILD_FREEZE_JOB_ID" ||
  test "$BUILD_FREEZE_JOB_TOKEN" = "${BUILD_FREEZE_JOB_ID};frontier"
printf 'build_freeze_job_id=%s\n' "$BUILD_FREEZE_JOB_ID"
)
```

If authenticated preflight capture stops after publishing only one evidence
mirror, or after publishing the complete pair but before its Orion binding,
preserve every evidence and staging artifact. Do not recapture and do not
manually remove anything. First identify one exact reviewed `PROBE_ID`, its
SHA-256, and the surviving role through read-only inspection. Then fill those
bindings into this standalone commit-forward recovery block. The authenticated
audit rejects ambiguous residue, malformed evidence, divergent mirrors and any
other one-sided pair before recovery mutates only the missing evidence mirror
and publishes the derived Orion binding:

```bash
(
set -euo pipefail
test "${PIC_ROOT:-}" = /lustre/orion/ast207/proj-shared/dfielding/PIC
test "${PYTHON:-}" = /opt/cray/pe/python/3.11.7/bin/python3
PROJECT_HOME_POLICY_ROOT=/autofs/nccs-svm1_proj/ast207/proj-shared/PIC
SOURCE_REPO=/autofs/nccs-svm1_home2/dfielding/athenak-pic
SLURM_ENV=(/usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin SLURM_CLUSTERS=frontier)
GIT=(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin
  /usr/bin/git -c core.fsmonitor=false -c core.hooksPath=/dev/null
)
FULL_GIT_COMMIT='<exact committed source identity used for interrupted capture>'
PROBE_ID='<exact reviewed interrupted-capture probe UUID>'
EXPECTED_EVIDENCE_SHA256='<sha256 of the exact surviving or complete evidence bytes>'
EXPECTED_EXISTING_ROLE='<orion_simulation_root or project_home_mirror_root>'
ACTIVE_PREDECESSOR_POLICY_SHA256=48e74f3151b51ba84b72e3214ef4a66c03441c745912b833395f98f6f60a0bd1
ACTIVE_PREDECESSOR_PROMOTION_SHA256=4ecb45bb399cee750ce69198286afc9fb903dde92cffc7507600b0a53f1c547f
[[ "$FULL_GIT_COMMIT" =~ ^[0-9a-f]{40}$ ]]
[[ "$PROBE_ID" =~ ^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$ ]]
[[ "$EXPECTED_EVIDENCE_SHA256" =~ ^[0-9a-f]{64}$ ]]
case "$EXPECTED_EXISTING_ROLE" in
  orion_simulation_root|project_home_mirror_root) ;;
  *) exit 1 ;;
esac
require_exact_reviewed_source() {
  local remote_pic_tip source_status
  source_status="$(
    "${GIT[@]}" -C "$SOURCE_REPO" status --porcelain --untracked-files=all
  )"
  test -z "$source_status"
  test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse HEAD)"
  test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse origin/PIC)"
  remote_pic_tip="$(
    "${GIT[@]}" -C "$SOURCE_REPO" ls-remote --exit-code origin refs/heads/PIC |
      /usr/bin/awk '$2 == "refs/heads/PIC" {count += 1; tip = $1} END {if (count != 1) exit 1; print tip}'
  )"
  [[ "$remote_pic_tip" =~ ^[0-9a-f]{40}$ ]]
  test "$FULL_GIT_COMMIT" = "$remote_pic_tip"
}
SOURCE_CONTROL_PLANE=(
  "$PYTHON" -I -B
  "${SOURCE_REPO}/tst/publication/frontier_control_plane/run_control_plane.py"
  --expected-git-commit "$FULL_GIT_COMMIT"
)
require_exact_reviewed_source
/usr/bin/cmp -s \
  "${PIC_ROOT}/policy/storage_policy.json" \
  "${PROJECT_HOME_POLICY_ROOT}/policy/storage_policy.json"
/usr/bin/cmp -s \
  "${PIC_ROOT}/policy/active_promotion.json" \
  "${PROJECT_HOME_POLICY_ROOT}/policy/active_promotion.json"
test "$ACTIVE_PREDECESSOR_POLICY_SHA256" = "$(
  /usr/bin/sha256sum "${PIC_ROOT}/policy/storage_policy.json" |
    /usr/bin/awk '{print $1}'
)"
test "$ACTIVE_PREDECESSOR_PROMOTION_SHA256" = "$(
  /usr/bin/sha256sum "${PIC_ROOT}/policy/active_promotion.json" |
    /usr/bin/awk '{print $1}'
)"
/usr/bin/ps -u "$(/usr/bin/id -u)" -ww -o pid=,ppid=,lstart=,args=
"${SLURM_ENV[@]}" /usr/bin/squeue --clusters=frontier -u "$(/usr/bin/id -un)" \
  -h -o '%i|%j|%a|%q|%T|%Z|%o'
REVIEWED_NO_SCIENCE_OR_SAME_USER_PIC_MUTATOR='<yes only after reviewing both snapshots>'
test "$REVIEWED_NO_SCIENCE_OR_SAME_USER_PIC_MUTATOR" = yes
AUDIT="$(
  "${SOURCE_CONTROL_PLANE[@]}" capture_storage_preflight_evidence.py \
    --expected-git-commit "$FULL_GIT_COMMIT" \
    --audit-exact-pair \
    --probe-id "$PROBE_ID" \
    --expected-evidence-sha256 "$EXPECTED_EVIDENCE_SHA256"
)"
AUDIT_STATE="$(
  printf '%s\n' "$AUDIT" |
    "$PYTHON" -I -B -c 'import json, sys; print(json.load(sys.stdin)["state"])'
)"
case "$AUDIT_STATE:$EXPECTED_EXISTING_ROLE" in
  valid_identical_pair:orion_simulation_root|\
  valid_identical_pair:project_home_mirror_root|\
  valid_orion_only:orion_simulation_root|\
  valid_project_home_only:project_home_mirror_root) ;;
  *) exit 1 ;;
esac
require_exact_reviewed_source
STORAGE_PREFLIGHT_FRAGMENT="$(
  "${SOURCE_CONTROL_PLANE[@]}" capture_storage_preflight_evidence.py \
    --expected-git-commit "$FULL_GIT_COMMIT" \
    --recover-exact-pair \
    --probe-id "$PROBE_ID" \
    --expected-evidence-sha256 "$EXPECTED_EVIDENCE_SHA256" \
    --expected-existing-role "$EXPECTED_EXISTING_ROLE"
)"
test "$(
  printf '%s\n' "$STORAGE_PREFLIGHT_FRAGMENT" |
    "$PYTHON" -I -B -c \
      'import json, sys; value = json.load(sys.stdin)["storage_preflight_evidence"]; print(value["probe_id"] + "|" + value["sha256"])'
)" = "${PROBE_ID}|${EXPECTED_EVIDENCE_SHA256}"
test "$(
  "${SOURCE_CONTROL_PLANE[@]}" capture_storage_preflight_evidence.py \
    --expected-git-commit "$FULL_GIT_COMMIT" \
    --audit-exact-pair \
    --probe-id "$PROBE_ID" \
    --expected-evidence-sha256 "$EXPECTED_EVIDENCE_SHA256" |
    "$PYTHON" -I -B -c 'import json, sys; print(json.load(sys.stdin)["state"])'
)" = valid_identical_pair
STORAGE_PREFLIGHT_BINDING_DIR="${PIC_ROOT}/policy/storage_preflight_bindings"
test "$(
  /usr/bin/stat -Lc '%F|%a|%u' "$STORAGE_PREFLIGHT_BINDING_DIR"
)" = "directory|2700|$(/usr/bin/id -u)"
STORAGE_PREFLIGHT_BINDING_DIR_ID="$(
  /usr/bin/stat -Lc '%d:%i' "$STORAGE_PREFLIGHT_BINDING_DIR"
)"
STORAGE_PREFLIGHT_BINDING="${STORAGE_PREFLIGHT_BINDING_DIR}/${PROBE_ID}.json"
test ! -e "$STORAGE_PREFLIGHT_BINDING"
STORAGE_PREFLIGHT_STAGING="$(
  /usr/bin/mktemp -p "$STORAGE_PREFLIGHT_BINDING_DIR" .staging.XXXXXX
)"
printf '%s\n' "$STORAGE_PREFLIGHT_FRAGMENT" > "$STORAGE_PREFLIGHT_STAGING"
/usr/bin/chmod 0400 "$STORAGE_PREFLIGHT_STAGING"
/usr/bin/ln -T -- "$STORAGE_PREFLIGHT_STAGING" "$STORAGE_PREFLIGHT_BINDING"
/usr/bin/rm "$STORAGE_PREFLIGHT_STAGING"
test "$(
  /usr/bin/stat -Lc '%d:%i' "$STORAGE_PREFLIGHT_BINDING_DIR"
)" = "$STORAGE_PREFLIGHT_BINDING_DIR_ID"
"$PYTHON" -I -B -c \
  'import os, sys; fd = os.open(sys.argv[1], os.O_RDONLY | os.O_NOFOLLOW); os.fsync(fd); os.close(fd); fd = os.open(os.path.dirname(sys.argv[1]), os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW); os.fsync(fd); os.close(fd)' \
  "$STORAGE_PREFLIGHT_BINDING"
printf 'recovered_probe_id=%s\n' "$PROBE_ID"
printf 'recovered_storage_preflight_binding=%s\n' "$STORAGE_PREFLIGHT_BINDING"
)
```

If the session stops after exactly one `one_sided_*_control_plane` checkpoint,
do not rerun the primary block and do not remove any `.tmp-*` entry. A
`.tmp-*` entry requires reviewed manual recovery. Otherwise, fill the exact
printed installed path, its opposite missing root, the retained preflight
binding, and the committed source identity into this one-sided recovery block.
It re-verifies the unchanged exact active-predecessor launch-prohibited
generation and mutates
only the missing half of the pair:

```bash
(
set -euo pipefail
test "${PIC_ROOT:-}" = /lustre/orion/ast207/proj-shared/dfielding/PIC
test "${PYTHON:-}" = /opt/cray/pe/python/3.11.7/bin/python3
PROJECT_HOME_POLICY_ROOT=/autofs/nccs-svm1_proj/ast207/proj-shared/PIC
SOURCE_REPO=/autofs/nccs-svm1_home2/dfielding/athenak-pic
SLURM_ENV=(/usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin SLURM_CLUSTERS=frontier)
GIT=(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin
  /usr/bin/git -c core.fsmonitor=false -c core.hooksPath=/dev/null
)
FULL_GIT_COMMIT='<exact committed source identity used for preflight and first install>'
PROBE_ID='<probe_id from the retained published preflight binding>'
VERSION='<digest basename from the one_sided_*_control_plane checkpoint>'
EXISTING_CONTROL_PLANE='<exact one_sided_*_control_plane checkpoint path>'
MISSING_CONTROL_PLANE_ROOT='<the exact opposite root>'
ACTIVE_PREDECESSOR_CONTROL_PLANE_VERSION=b56d96b40f2c666b9fa5b421fac589d6a4d6500a716f20b479c354a3d239cb48
ACTIVE_PREDECESSOR_POLICY_SHA256=48e74f3151b51ba84b72e3214ef4a66c03441c745912b833395f98f6f60a0bd1
ACTIVE_PREDECESSOR_PROMOTION_SHA256=4ecb45bb399cee750ce69198286afc9fb903dde92cffc7507600b0a53f1c547f
PRESERVED_CLEAN_CANDIDATE_MANIFEST="${PIC_ROOT}/clean_candidates/98a372c9-2ea0-47e6-ad34-e66343e7eea1/clean_candidate_manifest.json"
PRESERVED_CLEAN_CANDIDATE_MANIFEST_SHA256=dee6be45657e99ec477eec513c45be4ba6ac43b4fd5deb5655750f51c668e42f
PRESERVED_CLEAN_CANDIDATE_BUILD_PROFILE_CONTROL_PLANE_VERSION=821d185856722bd0178acb9427f78ac82671a4b6670779ec8400fbac54c6d721
ACTIVE_PREDECESSOR_ORION_CONTROL_PLANE="${PIC_ROOT}/control_plane/${ACTIVE_PREDECESSOR_CONTROL_PLANE_VERSION}"
ACTIVE_PREDECESSOR_PROJECT_HOME_CONTROL_PLANE="${PROJECT_HOME_POLICY_ROOT}/control_plane/${ACTIVE_PREDECESSOR_CONTROL_PLANE_VERSION}"
ACTIVE_PREDECESSOR_CONTROL_PLANE=("$PYTHON" -I -B "${ACTIVE_PREDECESSOR_ORION_CONTROL_PLANE}/run_control_plane.py")
[[ "$FULL_GIT_COMMIT" =~ ^[0-9a-f]{40}$ ]]
[[ "$PROBE_ID" =~ ^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$ ]]
[[ "$VERSION" =~ ^[0-9a-f]{64}$ ]]
test -r "${PIC_ROOT}/policy/storage_preflight_bindings/${PROBE_ID}.json"
if test "$EXISTING_CONTROL_PLANE" = "${PIC_ROOT}/control_plane/${VERSION}"; then
  test "$MISSING_CONTROL_PLANE_ROOT" = "$PROJECT_HOME_POLICY_ROOT"
else
  test "$EXISTING_CONTROL_PLANE" = "${PROJECT_HOME_POLICY_ROOT}/control_plane/${VERSION}"
  test "$MISSING_CONTROL_PLANE_ROOT" = "$PIC_ROOT"
fi
test -d "$EXISTING_CONTROL_PLANE"
test ! -e "${MISSING_CONTROL_PLANE_ROOT}/control_plane/${VERSION}"
EXISTING_CONTROL_PLANE_COMMAND=(
  "$PYTHON" -I -B "${EXISTING_CONTROL_PLANE}/run_control_plane.py"
)
test "$(
  "${EXISTING_CONTROL_PLANE_COMMAND[@]}" \
    validate_and_reserve_frontier_job.py verify-control-plane
)" = "$VERSION"
SOURCE_STATUS="$(
  "${GIT[@]}" -C "$SOURCE_REPO" status --porcelain --untracked-files=all
)"
test -z "$SOURCE_STATUS"
test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse HEAD)"
test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse origin/PIC)"
REMOTE_PIC_TIP="$(
  "${GIT[@]}" -C "$SOURCE_REPO" ls-remote --exit-code origin refs/heads/PIC |
    /usr/bin/awk '$2 == "refs/heads/PIC" {count += 1; tip = $1} END {if (count != 1) exit 1; print tip}'
)"
[[ "$REMOTE_PIC_TIP" =~ ^[0-9a-f]{40}$ ]]
test "$FULL_GIT_COMMIT" = "$REMOTE_PIC_TIP"
require_exact_reviewed_source() {
  local remote_pic_tip source_status
  source_status="$(
    "${GIT[@]}" -C "$SOURCE_REPO" status --porcelain --untracked-files=all
  )"
  test -z "$source_status"
  test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse HEAD)"
  test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse origin/PIC)"
  remote_pic_tip="$(
    "${GIT[@]}" -C "$SOURCE_REPO" ls-remote --exit-code origin refs/heads/PIC |
      /usr/bin/awk '$2 == "refs/heads/PIC" {count += 1; tip = $1} END {if (count != 1) exit 1; print tip}'
  )"
  [[ "$remote_pic_tip" =~ ^[0-9a-f]{40}$ ]]
  test "$FULL_GIT_COMMIT" = "$remote_pic_tip"
}
SOURCE_CONTROL_PLANE=(
  "$PYTHON" -I -B
  "${SOURCE_REPO}/tst/publication/frontier_control_plane/run_control_plane.py"
  --expected-git-commit "$FULL_GIT_COMMIT"
)
for POLICY_PARENT in "${PIC_ROOT}/policy" "${PROJECT_HOME_POLICY_ROOT}/policy"; do
  test -d "$POLICY_PARENT"
  test ! -L "$POLICY_PARENT"
  RECOVERY_ENTRY="$(
    /usr/bin/find "$POLICY_PARENT" -mindepth 1 -maxdepth 1 \
      \( -name '.active_promotion_transaction.json' \
      -o -name '.storage_policy.json.transaction-rollback-*' \
      -o -name '.active_promotion.json.transaction-rollback-*' \
      -o -name '.storage_policy.json.tmp-*' \
      -o -name '.active_promotion.json.tmp-*' \
      -o -name '.storage_policy.json.rollback-*' \
      -o -name '.active_promotion.json.rollback-*' \) -print -quit
  )"
  test -z "$RECOVERY_ENTRY"
done
for STAGING_PARENT in \
  "${PIC_ROOT}/policy/storage_preflight_bindings" \
  "${PIC_ROOT}/policy/storage_preflight_evidence" \
  "${PROJECT_HOME_POLICY_ROOT}/policy/storage_preflight_evidence" \
  "${PIC_ROOT}/control_plane" \
  "${PROJECT_HOME_POLICY_ROOT}/control_plane" \
  "${PIC_ROOT}/clean_candidates"; do
  test -d "$STAGING_PARENT"
  test ! -L "$STAGING_PARENT"
  STAGING_ENTRY="$(
    /usr/bin/find "$STAGING_PARENT" -mindepth 1 -maxdepth 1 \
      \( -name '.staging.*' -o -name '.tmp-*' -o -name '*.recovery-staging' \) \
      -print -quit
  )"
  test -z "$STAGING_ENTRY"
done
for STORAGE_ROOT in "$PIC_ROOT" "$PROJECT_HOME_POLICY_ROOT"; do
  STORAGE_PROBE_ENTRY="$(
    /usr/bin/find "$STORAGE_ROOT" -mindepth 1 -maxdepth 1 \
      -name '.pic-storage-preflight-*' -print -quit
  )"
  test -z "$STORAGE_PROBE_ENTRY"
done
ORION_EVIDENCE="$(
  /usr/bin/find "${PIC_ROOT}/policy/storage_preflight_evidence" \
    -mindepth 1 -maxdepth 1 -type f -name '*.json' -printf '%f\n' |
    /usr/bin/sort
)"
PROJECT_HOME_EVIDENCE="$(
  /usr/bin/find "${PROJECT_HOME_POLICY_ROOT}/policy/storage_preflight_evidence" \
    -mindepth 1 -maxdepth 1 -type f -name '*.json' -printf '%f\n' |
    /usr/bin/sort
)"
test -n "$ORION_EVIDENCE"
test "$ORION_EVIDENCE" = "$PROJECT_HOME_EVIDENCE"
while IFS= read -r EVIDENCE_NAME; do
  /usr/bin/cmp -s \
    "${PIC_ROOT}/policy/storage_preflight_evidence/${EVIDENCE_NAME}" \
    "${PROJECT_HOME_POLICY_ROOT}/policy/storage_preflight_evidence/${EVIDENCE_NAME}"
done <<< "$ORION_EVIDENCE"
/usr/bin/cmp -s \
  "${ACTIVE_PREDECESSOR_ORION_CONTROL_PLANE}/inventory.json" \
  "${ACTIVE_PREDECESSOR_PROJECT_HOME_CONTROL_PLANE}/inventory.json"
test "$(
  "${ACTIVE_PREDECESSOR_CONTROL_PLANE[@]}" \
    validate_and_reserve_frontier_job.py verify-control-plane
)" = "$ACTIVE_PREDECESSOR_CONTROL_PLANE_VERSION"
/usr/bin/cmp -s \
  "${PIC_ROOT}/policy/storage_policy.json" \
  "${PROJECT_HOME_POLICY_ROOT}/policy/storage_policy.json"
/usr/bin/cmp -s \
  "${PIC_ROOT}/policy/active_promotion.json" \
  "${PROJECT_HOME_POLICY_ROOT}/policy/active_promotion.json"
test "$ACTIVE_PREDECESSOR_POLICY_SHA256" = "$(
  /usr/bin/sha256sum "${PIC_ROOT}/policy/storage_policy.json" |
    /usr/bin/awk '{print $1}'
)"
test "$ACTIVE_PREDECESSOR_PROMOTION_SHA256" = "$(
  /usr/bin/sha256sum "${PIC_ROOT}/policy/active_promotion.json" |
    /usr/bin/awk '{print $1}'
)"
"$PYTHON" -I -B - \
  "$ACTIVE_PREDECESSOR_ORION_CONTROL_PLANE" \
  "$ACTIVE_PREDECESSOR_PROJECT_HOME_CONTROL_PLANE" \
  "$PIC_ROOT" \
  "$PROJECT_HOME_POLICY_ROOT" \
  "$ACTIVE_PREDECESSOR_CONTROL_PLANE_VERSION" \
  "$ACTIVE_PREDECESSOR_POLICY_SHA256" \
  "$ACTIVE_PREDECESSOR_PROMOTION_SHA256" \
  "$PRESERVED_CLEAN_CANDIDATE_MANIFEST" \
  "$PRESERVED_CLEAN_CANDIDATE_MANIFEST_SHA256" \
  "$PRESERVED_CLEAN_CANDIDATE_BUILD_PROFILE_CONTROL_PLANE_VERSION" <<'PY'
import sys
from pathlib import Path

old_orion = Path(sys.argv[1])
old_project_home = Path(sys.argv[2])
pic_root = Path(sys.argv[3])
project_home_root = Path(sys.argv[4])
version, policy_sha256, promotion_sha256 = sys.argv[5:8]
preserved_manifest, preserved_manifest_sha256, preserved_build_controller = sys.argv[8:11]
sys.path.insert(0, str(old_orion))

from control_plane_common import require_storage_policy_unlock_snapshot
from control_plane_common import verify_installed_control_plane
from control_plane_common import project_home_ledger_root
from ledger import require_no_incomplete_manual_accounting_marker
from promote_active_policy import _require_no_outstanding_submissions

orion_inventory = verify_installed_control_plane(
    old_orion,
    authorized_pic_root=pic_root,
)
project_home_inventory = verify_installed_control_plane(
    old_project_home,
    authorized_pic_root=project_home_root,
)
assert orion_inventory == project_home_inventory
assert orion_inventory["version"] == version

policy, snapshot = require_storage_policy_unlock_snapshot(
    control_plane_version=version,
    authorized_pic_root=pic_root,
    authorized_project_home_root=project_home_root,
    allow_pending_genesis=True,
)
assert snapshot == {
    "active_policy_sha256": policy_sha256,
    "active_promotion_sha256": promotion_sha256,
}
assert policy["registered_science_slices"] == []
assert policy["science_submission_freeze"] == {
    "build_profile_control_plane_version": preserved_build_controller,
    "manifest_path": preserved_manifest,
    "manifest_sha256": preserved_manifest_sha256,
    "status": "authorized",
}
require_no_incomplete_manual_accounting_marker(
    pic_root / "ledger" / "node_hours.jsonl",
    project_home_ledger_root(project_home_root) / "ledger" / "node_hours.jsonl",
)
_require_no_outstanding_submissions(policy, pic_root, project_home_root)
PY
"${ACTIVE_PREDECESSOR_CONTROL_PLANE[@]}" revalidate_clean_candidate.py \
  --manifest "$PRESERVED_CLEAN_CANDIDATE_MANIFEST" \
  --expected-manifest-sha256 "$PRESERVED_CLEAN_CANDIDATE_MANIFEST_SHA256"
/usr/bin/ps -u "$(/usr/bin/id -u)" -ww -o pid=,ppid=,lstart=,args=
"${SLURM_ENV[@]}" /usr/bin/squeue --clusters=frontier -u "$(/usr/bin/id -un)" \
  -h -o '%i|%j|%a|%q|%T|%Z|%o'
REVIEWED_NO_SCIENCE_OR_SAME_USER_PIC_MUTATOR='<yes only after reviewing both snapshots>'
test "$REVIEWED_NO_SCIENCE_OR_SAME_USER_PIC_MUTATOR" = yes
require_exact_reviewed_source
RECOVERED_CONTROL_PLANE="$(
  "${SOURCE_CONTROL_PLANE[@]}" install_control_plane.py \
    --pic-root "$MISSING_CONTROL_PLANE_ROOT" \
    --expected-git-commit "$FULL_GIT_COMMIT"
)"
test "$RECOVERED_CONTROL_PLANE" = \
  "${MISSING_CONTROL_PLANE_ROOT}/control_plane/${VERSION}"
/usr/bin/cmp -s \
  "${EXISTING_CONTROL_PLANE}/inventory.json" \
  "${RECOVERED_CONTROL_PLANE}/inventory.json"
test -z "$(
  /usr/bin/find "${MISSING_CONTROL_PLANE_ROOT}/control_plane" \
    -mindepth 1 -maxdepth 1 -name '.tmp-*' -print -quit
)"
printf 'recovered_paired_control_plane_version=%s\n' "$VERSION"
)
```

After this recovery succeeds, or after an interruption that left both
controller installs complete but did not start migration, fill the exact
retained checkpoints into this standalone post-install resume block. It
revalidates the completed repair worker, clean true-remote source, paired
controller, active-predecessor launch-prohibited generation, complete mirrored evidence pairs
and empty recovery namespaces. It executes the Q011 migration materializer from
a fresh detached worktree snapshot at the reviewed commit, performs only the
one-use migration, verifies the resulting launch-prohibited generation, and
stops before any `sbatch`:

```bash
(
set -euo pipefail
test "${PIC_ROOT:-}" = /lustre/orion/ast207/proj-shared/dfielding/PIC
test "${PYTHON:-}" = /opt/cray/pe/python/3.11.7/bin/python3
PROJECT_HOME_POLICY_ROOT=/autofs/nccs-svm1_proj/ast207/proj-shared/PIC
SOURCE_REPO=/autofs/nccs-svm1_home2/dfielding/athenak-pic
SLURM_ENV=(/usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin SLURM_CLUSTERS=frontier)
GIT=(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin
  /usr/bin/git -c core.fsmonitor=false -c core.hooksPath=/dev/null
)
FULL_GIT_COMMIT='<exact committed source identity used for preflight and installs>'
REPAIR_VALIDATION_JOB_ID='<repair_validation_job_id from the completed checkpoint>'
PROBE_ID='<probe_id from the retained published preflight binding>'
VERSION='<shared digest basename from the completed paired installs>'
ACTIVE_PREDECESSOR_CONTROL_PLANE_VERSION=b56d96b40f2c666b9fa5b421fac589d6a4d6500a716f20b479c354a3d239cb48
ACTIVE_PREDECESSOR_POLICY_SHA256=48e74f3151b51ba84b72e3214ef4a66c03441c745912b833395f98f6f60a0bd1
ACTIVE_PREDECESSOR_PROMOTION_SHA256=4ecb45bb399cee750ce69198286afc9fb903dde92cffc7507600b0a53f1c547f
PRESERVED_CLEAN_CANDIDATE_MANIFEST="${PIC_ROOT}/clean_candidates/98a372c9-2ea0-47e6-ad34-e66343e7eea1/clean_candidate_manifest.json"
PRESERVED_CLEAN_CANDIDATE_MANIFEST_SHA256=dee6be45657e99ec477eec513c45be4ba6ac43b4fd5deb5655750f51c668e42f
PRESERVED_CLEAN_CANDIDATE_BUILD_PROFILE_CONTROL_PLANE_VERSION=821d185856722bd0178acb9427f78ac82671a4b6670779ec8400fbac54c6d721
[[ "$FULL_GIT_COMMIT" =~ ^[0-9a-f]{40}$ ]]
[[ "$REPAIR_VALIDATION_JOB_ID" =~ ^[0-9]+$ ]]
[[ "$PROBE_ID" =~ ^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$ ]]
[[ "$VERSION" =~ ^[0-9a-f]{64}$ ]]
require_exact_reviewed_source() {
  local remote_pic_tip source_status
  source_status="$(
    "${GIT[@]}" -C "$SOURCE_REPO" status --porcelain --untracked-files=all
  )"
  test -z "$source_status"
  test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse HEAD)"
  test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse origin/PIC)"
  remote_pic_tip="$(
    "${GIT[@]}" -C "$SOURCE_REPO" ls-remote --exit-code origin refs/heads/PIC |
      /usr/bin/awk '$2 == "refs/heads/PIC" {count += 1; tip = $1} END {if (count != 1) exit 1; print tip}'
  )"
  [[ "$remote_pic_tip" =~ ^[0-9a-f]{40}$ ]]
  test "$FULL_GIT_COMMIT" = "$remote_pic_tip"
}
require_no_policy_recovery_entries() {
  local policy_parent recovery_entry
  for policy_parent in "${PIC_ROOT}/policy" "${PROJECT_HOME_POLICY_ROOT}/policy"; do
    test -d "$policy_parent"
    test ! -L "$policy_parent"
    recovery_entry="$(
      /usr/bin/find "$policy_parent" -mindepth 1 -maxdepth 1 \
        \( -name '.active_promotion_transaction.json' \
        -o -name '.storage_policy.json.transaction-rollback-*' \
        -o -name '.active_promotion.json.transaction-rollback-*' \
        -o -name '.storage_policy.json.tmp-*' \
        -o -name '.active_promotion.json.tmp-*' \
        -o -name '.storage_policy.json.rollback-*' \
        -o -name '.active_promotion.json.rollback-*' \) \
        -print -quit
    )"
    test -z "$recovery_entry"
  done
}
require_no_staging_entries() {
  local staging_parent staging_entry
  for staging_parent in \
    "${PIC_ROOT}/policy/storage_preflight_bindings" \
    "${PIC_ROOT}/policy/storage_preflight_evidence" \
    "${PROJECT_HOME_POLICY_ROOT}/policy/storage_preflight_evidence" \
    "${PIC_ROOT}/control_plane" \
    "${PROJECT_HOME_POLICY_ROOT}/control_plane" \
    "${PIC_ROOT}/clean_candidates"; do
    test -d "$staging_parent"
    test ! -L "$staging_parent"
    staging_entry="$(
      /usr/bin/find "$staging_parent" -mindepth 1 -maxdepth 1 \
        \( -name '.staging.*' -o -name '.tmp-*' -o -name '*.recovery-staging' \) \
        -print -quit
    )"
    test -z "$staging_entry"
  done
  for staging_parent in "$PIC_ROOT" "$PROJECT_HOME_POLICY_ROOT"; do
    staging_entry="$(
      /usr/bin/find "$staging_parent" -mindepth 1 -maxdepth 1 \
        -name '.pic-storage-preflight-*' -print -quit
    )"
    test -z "$staging_entry"
  done
  local orion_evidence project_home_evidence evidence_name
  orion_evidence="$(
    /usr/bin/find "${PIC_ROOT}/policy/storage_preflight_evidence" \
      -mindepth 1 -maxdepth 1 -type f -name '*.json' -printf '%f\n' |
      /usr/bin/sort
  )"
  project_home_evidence="$(
    /usr/bin/find "${PROJECT_HOME_POLICY_ROOT}/policy/storage_preflight_evidence" \
      -mindepth 1 -maxdepth 1 -type f -name '*.json' -printf '%f\n' |
      /usr/bin/sort
  )"
  test -n "$orion_evidence"
  test "$orion_evidence" = "$project_home_evidence"
  while IFS= read -r evidence_name; do
    /usr/bin/cmp -s \
      "${PIC_ROOT}/policy/storage_preflight_evidence/${evidence_name}" \
      "${PROJECT_HOME_POLICY_ROOT}/policy/storage_preflight_evidence/${evidence_name}"
  done <<< "$orion_evidence"
}
run_reviewed_q011_materializer() {
  local snapshot snapshot_root status
  require_exact_reviewed_source
  snapshot_root="$(/usr/bin/mktemp -d /tmp/athenak-pic-q011-reviewed.XXXXXX)"
  snapshot="${snapshot_root}/source"
  status=0
  "${GIT[@]}" -C "$SOURCE_REPO" worktree add --detach "$snapshot" \
    "$FULL_GIT_COMMIT" || status=$?
  if test "$status" -eq 0; then
    if test -f "$snapshot/tst/publication/q011_section54_pressure_pilot_execution.py" &&
      test ! -L "$snapshot/tst/publication/q011_section54_pressure_pilot_execution.py"; then
      "$PYTHON" -I -B \
        "$snapshot/tst/publication/q011_section54_pressure_pilot_execution.py" \
        "$@" || status=$?
    else
      status=1
    fi
  fi
  if test -e "$snapshot"; then
    "${GIT[@]}" -C "$SOURCE_REPO" worktree remove --force "$snapshot" || status=1
  fi
  /usr/bin/rm -rf -- "$snapshot_root" || status=1
  return "$status"
}
ACTIVE_PREDECESSOR_ORION_CONTROL_PLANE="${PIC_ROOT}/control_plane/${ACTIVE_PREDECESSOR_CONTROL_PLANE_VERSION}"
ACTIVE_PREDECESSOR_PROJECT_HOME_CONTROL_PLANE="${PROJECT_HOME_POLICY_ROOT}/control_plane/${ACTIVE_PREDECESSOR_CONTROL_PLANE_VERSION}"
ACTIVE_PREDECESSOR_CONTROL_PLANE=("$PYTHON" -I -B "${ACTIVE_PREDECESSOR_ORION_CONTROL_PLANE}/run_control_plane.py")
ORION_CONTROL_PLANE="${PIC_ROOT}/control_plane/${VERSION}"
PROJECT_HOME_CONTROL_PLANE="${PROJECT_HOME_POLICY_ROOT}/control_plane/${VERSION}"
CONTROL_PLANE=("$PYTHON" -I -B "${ORION_CONTROL_PLANE}/run_control_plane.py")
STORAGE_PREFLIGHT_BINDING="${PIC_ROOT}/policy/storage_preflight_bindings/${PROBE_ID}.json"
POLICY_SUFFIX="${VERSION:0:8}_${PROBE_ID}"
MIGRATION_POLICY="${PIC_ROOT}/policy/reviewed_exact_preflight_predecessor_successor_${POLICY_SUFFIX}.json"
CANDIDATE_ONLY_POLICY="${PIC_ROOT}/policy/reviewed_candidate_only_successor_${POLICY_SUFFIX}.json"
require_exact_reviewed_source
REPAIR_VALIDATION_RECORD="$(
  "${SLURM_ENV[@]}" /usr/bin/sacct -X --clusters=frontier -n \
    -j "$REPAIR_VALIDATION_JOB_ID" \
    --format=JobIDRaw,JobName,Account,Partition,QOS,State,ExitCode,WorkDir --parsable2 |
    /usr/bin/awk -F'|' -v job="$REPAIR_VALIDATION_JOB_ID" '
      $1 == job {count += 1; record = $1 FS $2 FS $3 FS $4 FS $5 FS $6 FS $7 FS $8}
      END {if (count != 1) exit 1; print record}
    '
)"
test "$REPAIR_VALIDATION_RECORD" = \
  "${REPAIR_VALIDATION_JOB_ID}|pic-q011-pressure-gate-validate|ast207|batch|debug|COMPLETED|0:0|${SOURCE_REPO}"
REPAIR_VALIDATION_LOG="${PIC_ROOT}/logs/slurm/pic-q011-pressure-gate-validate.${REPAIR_VALIDATION_JOB_ID}.log"
test -r "$REPAIR_VALIDATION_LOG"
REPAIR_VALIDATION_LOG_COMMIT="$(
  /usr/bin/sed -n 's/^source_commit=//p' "$REPAIR_VALIDATION_LOG"
)"
test "$REPAIR_VALIDATION_LOG_COMMIT" = "$FULL_GIT_COMMIT"
test "$(printf '%s\n' "$REPAIR_VALIDATION_LOG_COMMIT" | /usr/bin/wc -l)" -eq 1
test -d "$ORION_CONTROL_PLANE"
test -d "$PROJECT_HOME_CONTROL_PLANE"
/usr/bin/cmp -s \
  "${ORION_CONTROL_PLANE}/inventory.json" \
  "${PROJECT_HOME_CONTROL_PLANE}/inventory.json"
test "$(
  "${CONTROL_PLANE[@]}" validate_and_reserve_frontier_job.py verify-control-plane
)" = "$VERSION"
test -r "$STORAGE_PREFLIGHT_BINDING"
test ! -e "$MIGRATION_POLICY"
test ! -e "$CANDIDATE_ONLY_POLICY"
require_no_policy_recovery_entries
require_no_staging_entries
/usr/bin/cmp -s \
  "${ACTIVE_PREDECESSOR_ORION_CONTROL_PLANE}/inventory.json" \
  "${ACTIVE_PREDECESSOR_PROJECT_HOME_CONTROL_PLANE}/inventory.json"
test "$(
  "${ACTIVE_PREDECESSOR_CONTROL_PLANE[@]}" validate_and_reserve_frontier_job.py verify-control-plane
)" = "$ACTIVE_PREDECESSOR_CONTROL_PLANE_VERSION"
/usr/bin/cmp -s \
  "${PIC_ROOT}/policy/storage_policy.json" \
  "${PROJECT_HOME_POLICY_ROOT}/policy/storage_policy.json"
/usr/bin/cmp -s \
  "${PIC_ROOT}/policy/active_promotion.json" \
  "${PROJECT_HOME_POLICY_ROOT}/policy/active_promotion.json"
test "$ACTIVE_PREDECESSOR_POLICY_SHA256" = "$(
  /usr/bin/sha256sum "${PIC_ROOT}/policy/storage_policy.json" |
    /usr/bin/awk '{print $1}'
)"
test "$ACTIVE_PREDECESSOR_PROMOTION_SHA256" = "$(
  /usr/bin/sha256sum "${PIC_ROOT}/policy/active_promotion.json" |
    /usr/bin/awk '{print $1}'
)"
"$PYTHON" -I -B - \
  "$ACTIVE_PREDECESSOR_ORION_CONTROL_PLANE" \
  "$ACTIVE_PREDECESSOR_PROJECT_HOME_CONTROL_PLANE" \
  "$PIC_ROOT" \
  "$PROJECT_HOME_POLICY_ROOT" \
  "$ACTIVE_PREDECESSOR_CONTROL_PLANE_VERSION" \
  "$ACTIVE_PREDECESSOR_POLICY_SHA256" \
  "$ACTIVE_PREDECESSOR_PROMOTION_SHA256" \
  "$PRESERVED_CLEAN_CANDIDATE_MANIFEST" \
  "$PRESERVED_CLEAN_CANDIDATE_MANIFEST_SHA256" \
  "$PRESERVED_CLEAN_CANDIDATE_BUILD_PROFILE_CONTROL_PLANE_VERSION" <<'PY'
import sys
from pathlib import Path

old_orion = Path(sys.argv[1])
old_project_home = Path(sys.argv[2])
pic_root = Path(sys.argv[3])
project_home_root = Path(sys.argv[4])
version, policy_sha256, promotion_sha256 = sys.argv[5:8]
preserved_manifest, preserved_manifest_sha256, preserved_build_controller = sys.argv[8:11]
sys.path.insert(0, str(old_orion))

from control_plane_common import require_storage_policy_unlock_snapshot
from control_plane_common import verify_installed_control_plane
from control_plane_common import project_home_ledger_root
from ledger import require_no_incomplete_manual_accounting_marker
from promote_active_policy import _require_no_outstanding_submissions

orion_inventory = verify_installed_control_plane(old_orion, authorized_pic_root=pic_root)
project_home_inventory = verify_installed_control_plane(
    old_project_home,
    authorized_pic_root=project_home_root,
)
assert orion_inventory == project_home_inventory
assert orion_inventory["version"] == version
policy, snapshot = require_storage_policy_unlock_snapshot(
    control_plane_version=version,
    authorized_pic_root=pic_root,
    authorized_project_home_root=project_home_root,
    allow_pending_genesis=True,
)
assert snapshot == {
    "active_policy_sha256": policy_sha256,
    "active_promotion_sha256": promotion_sha256,
}
assert policy["registered_science_slices"] == []
assert policy["science_submission_freeze"] == {
    "build_profile_control_plane_version": preserved_build_controller,
    "manifest_path": preserved_manifest,
    "manifest_sha256": preserved_manifest_sha256,
    "status": "authorized",
}
require_no_incomplete_manual_accounting_marker(
    pic_root / "ledger" / "node_hours.jsonl",
    project_home_ledger_root(project_home_root) / "ledger" / "node_hours.jsonl",
)
_require_no_outstanding_submissions(policy, pic_root, project_home_root)
PY
"${ACTIVE_PREDECESSOR_CONTROL_PLANE[@]}" revalidate_clean_candidate.py \
  --manifest "$PRESERVED_CLEAN_CANDIDATE_MANIFEST" \
  --expected-manifest-sha256 "$PRESERVED_CLEAN_CANDIDATE_MANIFEST_SHA256"
/usr/bin/ps -u "$(/usr/bin/id -u)" -ww -o pid=,ppid=,lstart=,args=
"${SLURM_ENV[@]}" /usr/bin/squeue --clusters=frontier -u "$(/usr/bin/id -un)" \
  -h -o '%i|%j|%a|%q|%T|%Z|%o'
REVIEWED_NO_SCIENCE_OR_SAME_USER_PIC_MUTATOR='<yes only after reviewing both snapshots>'
test "$REVIEWED_NO_SCIENCE_OR_SAME_USER_PIC_MUTATOR" = yes
run_reviewed_q011_materializer \
  exact-reviewed-preflight-predecessor-policy-successor \
  --expected-git-commit "$FULL_GIT_COMMIT" \
  --baseline-policy "${PIC_ROOT}/policy/storage_policy.json" \
  --control-plane-version "$VERSION" \
  --storage-preflight-binding "$STORAGE_PREFLIGHT_BINDING" \
  --output "$MIGRATION_POLICY"
require_no_policy_recovery_entries
require_no_staging_entries
"${CONTROL_PLANE[@]}" promote_active_policy.py \
  --reviewed-policy "$MIGRATION_POLICY" \
  --migrate-exact-reviewed-storage-preflight-predecessor
EXPECTED_ACTIVE_POLICY_SHA256="$(
  /usr/bin/sha256sum "${PIC_ROOT}/policy/storage_policy.json" |
    /usr/bin/awk '{print $1}'
)"
EXPECTED_ACTIVE_PROMOTION_SHA256="$(
  /usr/bin/sha256sum "${PIC_ROOT}/policy/active_promotion.json" |
    /usr/bin/awk '{print $1}'
)"
[[ "$EXPECTED_ACTIVE_POLICY_SHA256" =~ ^[0-9a-f]{64}$ ]]
[[ "$EXPECTED_ACTIVE_PROMOTION_SHA256" =~ ^[0-9a-f]{64}$ ]]
require_no_policy_recovery_entries
require_no_staging_entries
"${CONTROL_PLANE[@]}" promote_active_policy.py \
  --verify-active-launch-prohibited-generation \
  --expected-control-plane-version "$VERSION" \
  --expected-active-policy-sha256 "$EXPECTED_ACTIVE_POLICY_SHA256" \
  --expected-active-promotion-sha256 "$EXPECTED_ACTIVE_PROMOTION_SHA256" \
  --expected-authorized-freeze-manifest "$PRESERVED_CLEAN_CANDIDATE_MANIFEST" \
  --expected-authorized-freeze-manifest-sha256 "$PRESERVED_CLEAN_CANDIDATE_MANIFEST_SHA256" \
  --expected-authorized-freeze-build-controller \
    "$PRESERVED_CLEAN_CANDIDATE_BUILD_PROFILE_CONTROL_PLANE_VERSION"
printf 'full_git_commit=%s\n' "$FULL_GIT_COMMIT"
printf 'probe_id=%s\n' "$PROBE_ID"
printf 'control_plane_version=%s\n' "$VERSION"
printf 'expected_active_policy_sha256=%s\n' "$EXPECTED_ACTIVE_POLICY_SHA256"
printf 'expected_active_promotion_sha256=%s\n' "$EXPECTED_ACTIVE_PROMOTION_SHA256"
printf 'preserved_clean_candidate_manifest=%s\n' "$PRESERVED_CLEAN_CANDIDATE_MANIFEST"
printf 'preserved_clean_candidate_manifest_sha256=%s\n' "$PRESERVED_CLEAN_CANDIDATE_MANIFEST_SHA256"
printf 'preserved_clean_candidate_build_profile_control_plane_version=%s\n' \
  "$PRESERVED_CLEAN_CANDIDATE_BUILD_PROFILE_CONTROL_PLANE_VERSION"
)
```

After this standalone resume succeeds, use only the worker-only resume block
below with the eight printed post-migration values. Do not recapture preflight,
reinstall either side, or replay the one-use policy migration.

Preserve the eight values printed immediately after exact predecessor migration
as the post-migration checkpoint. They are emitted before `sbatch`, so an
interruption cannot require replaying the one-use migration to recover the
candidate-replacement CAS anchors. Each primary or resumed build/freeze
submission must pass exactly seven positional worker arguments, in this order:
commit, successor controller, active policy SHA-256, active promotion SHA-256,
preserved manifest path, preserved manifest SHA-256, and preserved build
controller. Do not replace them with named flags.

If the session stops after migration but before build/freeze submission, run
the worker-only resume block below. If it stops during `sbatch`, first perform
reviewed scheduler inspection and recover the accepted job ID when one exists;
run this block only when that inspection proves no build/freeze job was
accepted. It rechecks the exact post-migration anchors and submits only the
build/freeze worker. It does not rerun preflight capture, installation, policy
migration, reservation, or any science launch:

```bash
(
set -euo pipefail
test "${PIC_ROOT:-}" = /lustre/orion/ast207/proj-shared/dfielding/PIC
test "${PYTHON:-}" = /opt/cray/pe/python/3.11.7/bin/python3
PROJECT_HOME_POLICY_ROOT=/autofs/nccs-svm1_proj/ast207/proj-shared/PIC
SOURCE_REPO=/autofs/nccs-svm1_home2/dfielding/athenak-pic
SLURM_ENV=(/usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin SLURM_CLUSTERS=frontier)
GIT=(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin
  /usr/bin/git -c core.fsmonitor=false -c core.hooksPath=/dev/null
)
FULL_GIT_COMMIT='<full_git_commit from the post-migration checkpoint>'
PROBE_ID='<probe_id from the post-migration checkpoint>'
VERSION='<control_plane_version from the post-migration checkpoint>'
EXPECTED_ACTIVE_POLICY_SHA256='<expected_active_policy_sha256 from the post-migration checkpoint>'
EXPECTED_ACTIVE_PROMOTION_SHA256='<expected_active_promotion_sha256 from the post-migration checkpoint>'
PRESERVED_CLEAN_CANDIDATE_MANIFEST='<preserved_clean_candidate_manifest from the post-migration checkpoint>'
PRESERVED_CLEAN_CANDIDATE_MANIFEST_SHA256='<preserved_clean_candidate_manifest_sha256 from the post-migration checkpoint>'
PRESERVED_CLEAN_CANDIDATE_BUILD_PROFILE_CONTROL_PLANE_VERSION='<preserved_clean_candidate_build_profile_control_plane_version from the post-migration checkpoint>'
REVIEWED_NO_ACCEPTED_BUILD_FREEZE_JOB='<yes only after reviewed scheduler inspection>'
[[ "$FULL_GIT_COMMIT" =~ ^[0-9a-f]{40}$ ]]
[[ "$PROBE_ID" =~ ^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$ ]]
[[ "$VERSION" =~ ^[0-9a-f]{64}$ ]]
[[ "$EXPECTED_ACTIVE_POLICY_SHA256" =~ ^[0-9a-f]{64}$ ]]
[[ "$EXPECTED_ACTIVE_PROMOTION_SHA256" =~ ^[0-9a-f]{64}$ ]]
test "$PRESERVED_CLEAN_CANDIDATE_MANIFEST" = \
  "${PIC_ROOT}/clean_candidates/98a372c9-2ea0-47e6-ad34-e66343e7eea1/clean_candidate_manifest.json"
[[ "$PRESERVED_CLEAN_CANDIDATE_MANIFEST_SHA256" =~ ^[0-9a-f]{64}$ ]]
[[ "$PRESERVED_CLEAN_CANDIDATE_BUILD_PROFILE_CONTROL_PLANE_VERSION" =~ ^[0-9a-f]{64}$ ]]
test "$REVIEWED_NO_ACCEPTED_BUILD_FREEZE_JOB" = yes
require_no_policy_recovery_entries() {
  local policy_parent recovery_entry
  for policy_parent in "${PIC_ROOT}/policy" "${PROJECT_HOME_POLICY_ROOT}/policy"; do
    test -d "$policy_parent"
    test ! -L "$policy_parent"
    recovery_entry="$(
      /usr/bin/find "$policy_parent" -mindepth 1 -maxdepth 1 \
        \( -name '.active_promotion_transaction.json' \
        -o -name '.storage_policy.json.transaction-rollback-*' \
        -o -name '.active_promotion.json.transaction-rollback-*' \
        -o -name '.storage_policy.json.tmp-*' \
        -o -name '.active_promotion.json.tmp-*' \
        -o -name '.storage_policy.json.rollback-*' \
        -o -name '.active_promotion.json.rollback-*' \) \
        -print -quit
    )"
    test -z "$recovery_entry"
  done
}
require_no_staging_entries() {
  local staging_parent staging_entry
  for staging_parent in \
    "${PIC_ROOT}/policy/storage_preflight_bindings" \
    "${PIC_ROOT}/policy/storage_preflight_evidence" \
    "${PROJECT_HOME_POLICY_ROOT}/policy/storage_preflight_evidence" \
    "${PIC_ROOT}/control_plane" \
    "${PROJECT_HOME_POLICY_ROOT}/control_plane" \
    "${PIC_ROOT}/clean_candidates"; do
    test -d "$staging_parent"
    test ! -L "$staging_parent"
    staging_entry="$(
      /usr/bin/find "$staging_parent" -mindepth 1 -maxdepth 1 \
        \( -name '.staging.*' -o -name '.tmp-*' -o -name '*.recovery-staging' \) \
        -print -quit
    )"
    test -z "$staging_entry"
  done
  for staging_parent in "$PIC_ROOT" "$PROJECT_HOME_POLICY_ROOT"; do
    staging_entry="$(
      /usr/bin/find "$staging_parent" -mindepth 1 -maxdepth 1 \
        -name '.pic-storage-preflight-*' -print -quit
    )"
    test -z "$staging_entry"
  done
  local orion_evidence project_home_evidence evidence_name
  orion_evidence="$(
    /usr/bin/find "${PIC_ROOT}/policy/storage_preflight_evidence" \
      -mindepth 1 -maxdepth 1 -type f -name '*.json' -printf '%f\n' |
      /usr/bin/sort
  )"
  project_home_evidence="$(
    /usr/bin/find "${PROJECT_HOME_POLICY_ROOT}/policy/storage_preflight_evidence" \
      -mindepth 1 -maxdepth 1 -type f -name '*.json' -printf '%f\n' |
      /usr/bin/sort
  )"
  test -n "$orion_evidence"
  test "$orion_evidence" = "$project_home_evidence"
  while IFS= read -r evidence_name; do
    /usr/bin/cmp -s \
      "${PIC_ROOT}/policy/storage_preflight_evidence/${evidence_name}" \
      "${PROJECT_HOME_POLICY_ROOT}/policy/storage_preflight_evidence/${evidence_name}"
  done <<< "$orion_evidence"
}
SOURCE_STATUS="$(
  "${GIT[@]}" -C "$SOURCE_REPO" status --porcelain --untracked-files=all
)"
test -z "$SOURCE_STATUS"
test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse HEAD)"
test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse origin/PIC)"
REMOTE_PIC_TIP="$(
  "${GIT[@]}" -C "$SOURCE_REPO" ls-remote --exit-code origin refs/heads/PIC |
    /usr/bin/awk '$2 == "refs/heads/PIC" {count += 1; tip = $1} END {if (count != 1) exit 1; print tip}'
)"
test "$FULL_GIT_COMMIT" = "$REMOTE_PIC_TIP"
ORION_CONTROL_PLANE="${PIC_ROOT}/control_plane/${VERSION}"
PROJECT_HOME_CONTROL_PLANE="${PROJECT_HOME_POLICY_ROOT}/control_plane/${VERSION}"
test -d "$ORION_CONTROL_PLANE"
test -d "$PROJECT_HOME_CONTROL_PLANE"
/usr/bin/cmp -s \
  "${ORION_CONTROL_PLANE}/inventory.json" \
  "${PROJECT_HOME_CONTROL_PLANE}/inventory.json"
CONTROL_PLANE=("$PYTHON" -I -B "${ORION_CONTROL_PLANE}/run_control_plane.py")
POLICY_SUFFIX="${VERSION:0:8}_${PROBE_ID}"
STORAGE_PREFLIGHT_BINDING="${PIC_ROOT}/policy/storage_preflight_bindings/${PROBE_ID}.json"
MIGRATION_POLICY="${PIC_ROOT}/policy/reviewed_exact_preflight_predecessor_successor_${POLICY_SUFFIX}.json"
CANDIDATE_ONLY_POLICY="${PIC_ROOT}/policy/reviewed_candidate_only_successor_${POLICY_SUFFIX}.json"
test -r "$STORAGE_PREFLIGHT_BINDING"
test -r "$MIGRATION_POLICY"
test ! -e "$CANDIDATE_ONLY_POLICY"
test "$EXPECTED_ACTIVE_POLICY_SHA256" = "$(
  /usr/bin/sha256sum "${PIC_ROOT}/policy/storage_policy.json" |
    /usr/bin/awk '{print $1}'
)"
test "$EXPECTED_ACTIVE_PROMOTION_SHA256" = "$(
  /usr/bin/sha256sum "${PIC_ROOT}/policy/active_promotion.json" |
    /usr/bin/awk '{print $1}'
)"
require_no_policy_recovery_entries
require_no_staging_entries
"${CONTROL_PLANE[@]}" promote_active_policy.py \
  --verify-active-launch-prohibited-generation \
  --expected-control-plane-version "$VERSION" \
  --expected-active-policy-sha256 "$EXPECTED_ACTIVE_POLICY_SHA256" \
  --expected-active-promotion-sha256 "$EXPECTED_ACTIVE_PROMOTION_SHA256" \
  --expected-authorized-freeze-manifest "$PRESERVED_CLEAN_CANDIDATE_MANIFEST" \
  --expected-authorized-freeze-manifest-sha256 "$PRESERVED_CLEAN_CANDIDATE_MANIFEST_SHA256" \
  --expected-authorized-freeze-build-controller \
    "$PRESERVED_CLEAN_CANDIDATE_BUILD_PROFILE_CONTROL_PLANE_VERSION"

/usr/bin/ps -u "$(/usr/bin/id -u)" -ww -o pid=,ppid=,lstart=,args=
"${SLURM_ENV[@]}" /usr/bin/squeue --clusters=frontier -u "$(/usr/bin/id -un)" \
  -h -o '%i|%j|%a|%q|%T|%Z|%o'
REVIEWED_NO_SCIENCE_OR_SAME_USER_PIC_MUTATOR='<yes only after reviewing both snapshots>'
test "$REVIEWED_NO_SCIENCE_OR_SAME_USER_PIC_MUTATOR" = yes
cd "$SOURCE_REPO"
BUILD_FREEZE_JOB_TOKEN="$("${SLURM_ENV[@]}" /usr/bin/sbatch --parsable --export=NIL \
  "${SOURCE_REPO}/tst/publication/frontier_q011_clean_candidate_build_freeze_job.sh" \
  "$FULL_GIT_COMMIT" "$VERSION" \
  "$EXPECTED_ACTIVE_POLICY_SHA256" "$EXPECTED_ACTIVE_PROMOTION_SHA256" \
  "$PRESERVED_CLEAN_CANDIDATE_MANIFEST" \
  "$PRESERVED_CLEAN_CANDIDATE_MANIFEST_SHA256" \
  "$PRESERVED_CLEAN_CANDIDATE_BUILD_PROFILE_CONTROL_PLANE_VERSION")"
printf 'build_freeze_job_token=%q\n' "$BUILD_FREEZE_JOB_TOKEN"
BUILD_FREEZE_JOB_ID="${BUILD_FREEZE_JOB_TOKEN%;frontier}"
[[ "$BUILD_FREEZE_JOB_ID" =~ ^[0-9]+$ ]]
test "$BUILD_FREEZE_JOB_TOKEN" = "$BUILD_FREEZE_JOB_ID" ||
  test "$BUILD_FREEZE_JOB_TOKEN" = "${BUILD_FREEZE_JOB_ID};frontier"
printf 'build_freeze_job_id=%s\n' "$BUILD_FREEZE_JOB_ID"
)
```

### Accepted build/freeze timeout successor recovery

Build/freeze job `4769975` is retained as an accepted failed attempt, not as a
worker-only-resume candidate. Its exact top-level scheduler state is
`TIMEOUT|0:0`, its log SHA-256 is
`d07a92eb89c5646bbb25977b9076b3e05e7a1b3935abdb1e1f4448741446cffa`,
and the log binds the seven post-migration worker inputs before ending during
the first source clone. It contains no `clean_candidate_manifest` binding.
No configure, build, freeze, revalidation, candidate-only policy, or science
launch occurred.

Preserve these exact incomplete paths in place as non-authoritative chronology:

```text
/lustre/orion/ast207/proj-shared/dfielding/PIC/build/f6471610a116/hip-mpi-release-paper-pic
/lustre/orion/ast207/proj-shared/dfielding/PIC/bin/f6471610a116/hip-mpi-release-paper-pic
```

Do not delete, rename, reuse, complete, or treat either path as build evidence.
Do not run the post-worker phase with job `4769975`, and do not use the
worker-only resume block with its source/controller checkpoint. The accepted
timeout exposed that local `git clone --no-hardlinks` copied the linked
worktree's complete shared object database into Lustre.

The reviewed successor repair uses a new exact source commit and new paired
controller. Its build authority transfers the already-authenticated top-level
revision through a standalone, detached, shallow `--no-local`
`--revision=<commit>` clone. Recursive submodules use standalone full transport
clones because their exact raw `git submodule status --recursive` provenance
includes history/tag-derived descriptions. The trusted Git capability gate
runs before any build or artifact path is created. Both paths reject
alternates, promisor/partial-clone configuration, symlinked or hard-linked
object storage, a non-detached or wrong head, and any identity drift; the
top-level clone additionally rejects extra reachable history, while submodules
reject shallow source history. Exact recursive submodule-status bytes are
compared immediately after materialization, before configure or build. Because
the canonical build path is commit-derived, the repaired successor uses a
different path while retaining the failed attempt untouched.

Before submitting the repaired build/freeze worker, run the full primary
no-science preflight, paired-install, and exact-predecessor migration phase
against the exact current active predecessor:

```text
policy SHA-256:    48e74f3151b51ba84b72e3214ef4a66c03441c745912b833395f98f6f60a0bd1
promotion SHA-256: 4ecb45bb399cee750ce69198286afc9fb903dde92cffc7507600b0a53f1c547f
controller:        b56d96b40f2c666b9fa5b421fac589d6a4d6500a716f20b479c354a3d239cb48
preflight probe:   bc399b56-8fbb-4b67-b1dc-5df1dbff62b6
preflight SHA-256: d4a289ce9f4cd7c406dbea8d457864790f3112488e47ea24766cb192beaafab2
```

This replacement migration must preserve the authorized historical clean
candidate, empty registered-science allowlist, launch prohibition, failed
attempt residue, and all prior evidence. It must capture a fresh strictly newer
preflight, pair-install the repaired controller, and print a new post-migration
checkpoint before any replacement submission. Only the replacement job ID from
that new checkpoint may enter the post-worker phase.

After the worker job completes successfully, fill the values printed by the
post-migration checkpoint and the accepted job ID into this independent
post-worker phase. It requires exactly one successful top-level Slurm record
and exact single-valued log bindings for all seven worker inputs before it
reads candidate output. It then independently revalidates the immutable
candidate through the installed controller before candidate-only promotion:

```bash
(
set -euo pipefail
test "${PIC_ROOT:-}" = /lustre/orion/ast207/proj-shared/dfielding/PIC
test "${PYTHON:-}" = /opt/cray/pe/python/3.11.7/bin/python3
PROJECT_HOME_POLICY_ROOT=/autofs/nccs-svm1_proj/ast207/proj-shared/PIC
SOURCE_REPO=/autofs/nccs-svm1_home2/dfielding/athenak-pic
SLURM_ENV=(/usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin SLURM_CLUSTERS=frontier)
GIT=(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin
  /usr/bin/git -c core.fsmonitor=false -c core.hooksPath=/dev/null
)
FULL_GIT_COMMIT='<full_git_commit from the post-migration checkpoint>'
PROBE_ID='<probe_id from the post-migration checkpoint>'
VERSION='<control_plane_version from the post-migration checkpoint>'
EXPECTED_ACTIVE_POLICY_SHA256='<expected_active_policy_sha256 from the post-migration checkpoint>'
EXPECTED_ACTIVE_PROMOTION_SHA256='<expected_active_promotion_sha256 from the post-migration checkpoint>'
PRESERVED_CLEAN_CANDIDATE_MANIFEST='<preserved_clean_candidate_manifest from the post-migration checkpoint>'
PRESERVED_CLEAN_CANDIDATE_MANIFEST_SHA256='<preserved_clean_candidate_manifest_sha256 from the post-migration checkpoint>'
PRESERVED_CLEAN_CANDIDATE_BUILD_PROFILE_CONTROL_PLANE_VERSION='<preserved_clean_candidate_build_profile_control_plane_version from the post-migration checkpoint>'
BUILD_FREEZE_JOB_ID='<build_freeze_job_id from the accepted build/freeze submission>'
[[ "$FULL_GIT_COMMIT" =~ ^[0-9a-f]{40}$ ]]
[[ "$PROBE_ID" =~ ^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$ ]]
[[ "$VERSION" =~ ^[0-9a-f]{64}$ ]]
[[ "$EXPECTED_ACTIVE_POLICY_SHA256" =~ ^[0-9a-f]{64}$ ]]
[[ "$EXPECTED_ACTIVE_PROMOTION_SHA256" =~ ^[0-9a-f]{64}$ ]]
test "$PRESERVED_CLEAN_CANDIDATE_MANIFEST" = \
  "${PIC_ROOT}/clean_candidates/98a372c9-2ea0-47e6-ad34-e66343e7eea1/clean_candidate_manifest.json"
[[ "$PRESERVED_CLEAN_CANDIDATE_MANIFEST_SHA256" =~ ^[0-9a-f]{64}$ ]]
[[ "$PRESERVED_CLEAN_CANDIDATE_BUILD_PROFILE_CONTROL_PLANE_VERSION" =~ ^[0-9a-f]{64}$ ]]
[[ "$BUILD_FREEZE_JOB_ID" =~ ^[0-9]+$ ]]
require_no_policy_recovery_entries() {
  local policy_parent recovery_entry
  for policy_parent in "${PIC_ROOT}/policy" "${PROJECT_HOME_POLICY_ROOT}/policy"; do
    test -d "$policy_parent"
    test ! -L "$policy_parent"
    recovery_entry="$(
      /usr/bin/find "$policy_parent" -mindepth 1 -maxdepth 1 \
        \( -name '.active_promotion_transaction.json' \
        -o -name '.storage_policy.json.transaction-rollback-*' \
        -o -name '.active_promotion.json.transaction-rollback-*' \
        -o -name '.storage_policy.json.tmp-*' \
        -o -name '.active_promotion.json.tmp-*' \
        -o -name '.storage_policy.json.rollback-*' \
        -o -name '.active_promotion.json.rollback-*' \) \
        -print -quit
    )"
    test -z "$recovery_entry"
  done
}
require_no_staging_entries() {
  local staging_parent staging_entry
  for staging_parent in \
    "${PIC_ROOT}/policy/storage_preflight_bindings" \
    "${PIC_ROOT}/policy/storage_preflight_evidence" \
    "${PROJECT_HOME_POLICY_ROOT}/policy/storage_preflight_evidence" \
    "${PIC_ROOT}/control_plane" \
    "${PROJECT_HOME_POLICY_ROOT}/control_plane" \
    "${PIC_ROOT}/clean_candidates"; do
    test -d "$staging_parent"
    test ! -L "$staging_parent"
    staging_entry="$(
      /usr/bin/find "$staging_parent" -mindepth 1 -maxdepth 1 \
        \( -name '.staging.*' -o -name '.tmp-*' -o -name '*.recovery-staging' \) \
        -print -quit
    )"
    test -z "$staging_entry"
  done
  for staging_parent in "$PIC_ROOT" "$PROJECT_HOME_POLICY_ROOT"; do
    staging_entry="$(
      /usr/bin/find "$staging_parent" -mindepth 1 -maxdepth 1 \
        -name '.pic-storage-preflight-*' -print -quit
    )"
    test -z "$staging_entry"
  done
  local orion_evidence project_home_evidence evidence_name
  orion_evidence="$(
    /usr/bin/find "${PIC_ROOT}/policy/storage_preflight_evidence" \
      -mindepth 1 -maxdepth 1 -type f -name '*.json' -printf '%f\n' |
      /usr/bin/sort
  )"
  project_home_evidence="$(
    /usr/bin/find "${PROJECT_HOME_POLICY_ROOT}/policy/storage_preflight_evidence" \
      -mindepth 1 -maxdepth 1 -type f -name '*.json' -printf '%f\n' |
      /usr/bin/sort
  )"
  test -n "$orion_evidence"
  test "$orion_evidence" = "$project_home_evidence"
  while IFS= read -r evidence_name; do
    /usr/bin/cmp -s \
      "${PIC_ROOT}/policy/storage_preflight_evidence/${evidence_name}" \
      "${PROJECT_HOME_POLICY_ROOT}/policy/storage_preflight_evidence/${evidence_name}"
  done <<< "$orion_evidence"
}
SOURCE_STATUS="$(
  "${GIT[@]}" -C "$SOURCE_REPO" status --porcelain --untracked-files=all
)"
test -z "$SOURCE_STATUS"
test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse HEAD)"
test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse origin/PIC)"
REMOTE_PIC_TIP="$(
  "${GIT[@]}" -C "$SOURCE_REPO" ls-remote --exit-code origin refs/heads/PIC |
    /usr/bin/awk '$2 == "refs/heads/PIC" {count += 1; tip = $1} END {if (count != 1) exit 1; print tip}'
)"
[[ "$REMOTE_PIC_TIP" =~ ^[0-9a-f]{40}$ ]]
test "$FULL_GIT_COMMIT" = "$REMOTE_PIC_TIP"
require_exact_reviewed_source() {
  local remote_pic_tip source_status
  source_status="$(
    "${GIT[@]}" -C "$SOURCE_REPO" status --porcelain --untracked-files=all
  )"
  test -z "$source_status"
  test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse HEAD)"
  test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse origin/PIC)"
  remote_pic_tip="$(
    "${GIT[@]}" -C "$SOURCE_REPO" ls-remote --exit-code origin refs/heads/PIC |
      /usr/bin/awk '$2 == "refs/heads/PIC" {count += 1; tip = $1} END {if (count != 1) exit 1; print tip}'
  )"
  [[ "$remote_pic_tip" =~ ^[0-9a-f]{40}$ ]]
  test "$FULL_GIT_COMMIT" = "$remote_pic_tip"
}
run_reviewed_q011_materializer() {
  local snapshot snapshot_root status
  require_exact_reviewed_source
  snapshot_root="$(/usr/bin/mktemp -d /tmp/athenak-pic-q011-reviewed.XXXXXX)"
  snapshot="${snapshot_root}/source"
  status=0
  "${GIT[@]}" -C "$SOURCE_REPO" worktree add --detach "$snapshot" \
    "$FULL_GIT_COMMIT" || status=$?
  if test "$status" -eq 0; then
    if test -f "$snapshot/tst/publication/q011_section54_pressure_pilot_execution.py" &&
      test ! -L "$snapshot/tst/publication/q011_section54_pressure_pilot_execution.py"; then
      "$PYTHON" -I -B \
        "$snapshot/tst/publication/q011_section54_pressure_pilot_execution.py" \
        "$@" || status=$?
    else
      status=1
    fi
  fi
  if test -e "$snapshot"; then
    "${GIT[@]}" -C "$SOURCE_REPO" worktree remove --force "$snapshot" || status=1
  fi
  /usr/bin/rm -rf -- "$snapshot_root" || status=1
  return "$status"
}
ORION_CONTROL_PLANE="${PIC_ROOT}/control_plane/${VERSION}"
CONTROL_PLANE=("$PYTHON" -I -B "${ORION_CONTROL_PLANE}/run_control_plane.py")
POLICY_SUFFIX="${VERSION:0:8}_${PROBE_ID}"
STORAGE_PREFLIGHT_BINDING="${PIC_ROOT}/policy/storage_preflight_bindings/${PROBE_ID}.json"
MIGRATION_POLICY="${PIC_ROOT}/policy/reviewed_exact_preflight_predecessor_successor_${POLICY_SUFFIX}.json"
CANDIDATE_ONLY_POLICY="${PIC_ROOT}/policy/reviewed_candidate_only_successor_${POLICY_SUFFIX}.json"
test -r "$STORAGE_PREFLIGHT_BINDING"
test -r "$MIGRATION_POLICY"
test ! -e "$CANDIDATE_ONLY_POLICY"
test "$EXPECTED_ACTIVE_POLICY_SHA256" = "$(
  /usr/bin/sha256sum "${PIC_ROOT}/policy/storage_policy.json" |
    /usr/bin/awk '{print $1}'
)"
test "$EXPECTED_ACTIVE_PROMOTION_SHA256" = "$(
  /usr/bin/sha256sum "${PIC_ROOT}/policy/active_promotion.json" |
    /usr/bin/awk '{print $1}'
)"
require_no_policy_recovery_entries
require_no_staging_entries
"${CONTROL_PLANE[@]}" promote_active_policy.py \
  --verify-active-launch-prohibited-generation \
  --expected-control-plane-version "$VERSION" \
  --expected-active-policy-sha256 "$EXPECTED_ACTIVE_POLICY_SHA256" \
  --expected-active-promotion-sha256 "$EXPECTED_ACTIVE_PROMOTION_SHA256" \
  --expected-authorized-freeze-manifest "$PRESERVED_CLEAN_CANDIDATE_MANIFEST" \
  --expected-authorized-freeze-manifest-sha256 "$PRESERVED_CLEAN_CANDIDATE_MANIFEST_SHA256" \
  --expected-authorized-freeze-build-controller \
    "$PRESERVED_CLEAN_CANDIDATE_BUILD_PROFILE_CONTROL_PLANE_VERSION"
BUILD_FREEZE_RECORD="$(
  "${SLURM_ENV[@]}" /usr/bin/sacct -X --clusters=frontier -n -j "$BUILD_FREEZE_JOB_ID" \
    --format=JobIDRaw,JobName,Account,Partition,QOS,State,ExitCode,WorkDir --parsable2 |
    /usr/bin/awk -F'|' -v job="$BUILD_FREEZE_JOB_ID" '
      $1 == job {count += 1; record = $1 FS $2 FS $3 FS $4 FS $5 FS $6 FS $7 FS $8}
      END {if (count != 1) exit 1; print record}
    '
)"
test "$BUILD_FREEZE_RECORD" = \
  "${BUILD_FREEZE_JOB_ID}|pic-q011-build-freeze|ast207|batch|debug|COMPLETED|0:0|${SOURCE_REPO}"
BUILD_FREEZE_LOG="${PIC_ROOT}/logs/slurm/pic-q011-build-freeze.${BUILD_FREEZE_JOB_ID}.log"
test -r "$BUILD_FREEZE_LOG"
read_exact_build_freeze_log_binding() {
  local key="$1" value
  value="$(
    /usr/bin/awk -F= -v key="$key" '
      $1 == key {count += 1; value = substr($0, length(key) + 2)}
      END {if (count != 1 || value == "") exit 1; print value}
    ' "$BUILD_FREEZE_LOG"
  )"
  printf '%s\n' "$value"
}
BUILD_FREEZE_LOG_COMMIT="$(
  read_exact_build_freeze_log_binding source_commit
)"
BUILD_FREEZE_LOG_CONTROL_PLANE_VERSION="$(
  read_exact_build_freeze_log_binding control_plane_version
)"
BUILD_FREEZE_LOG_EXPECTED_ACTIVE_POLICY_SHA256="$(
  read_exact_build_freeze_log_binding expected_active_policy_sha256
)"
BUILD_FREEZE_LOG_EXPECTED_ACTIVE_PROMOTION_SHA256="$(
  read_exact_build_freeze_log_binding expected_active_promotion_sha256
)"
BUILD_FREEZE_LOG_EXPECTED_AUTHORIZED_FREEZE_MANIFEST="$(
  read_exact_build_freeze_log_binding expected_authorized_freeze_manifest
)"
BUILD_FREEZE_LOG_EXPECTED_AUTHORIZED_FREEZE_MANIFEST_SHA256="$(
  read_exact_build_freeze_log_binding expected_authorized_freeze_manifest_sha256
)"
BUILD_FREEZE_LOG_EXPECTED_AUTHORIZED_FREEZE_BUILD_CONTROLLER="$(
  read_exact_build_freeze_log_binding expected_authorized_freeze_build_controller
)"
test "$BUILD_FREEZE_LOG_COMMIT" = "$FULL_GIT_COMMIT"
test "$BUILD_FREEZE_LOG_CONTROL_PLANE_VERSION" = "$VERSION"
test "$BUILD_FREEZE_LOG_EXPECTED_ACTIVE_POLICY_SHA256" = \
  "$EXPECTED_ACTIVE_POLICY_SHA256"
test "$BUILD_FREEZE_LOG_EXPECTED_ACTIVE_PROMOTION_SHA256" = \
  "$EXPECTED_ACTIVE_PROMOTION_SHA256"
test "$BUILD_FREEZE_LOG_EXPECTED_AUTHORIZED_FREEZE_MANIFEST" = \
  "$PRESERVED_CLEAN_CANDIDATE_MANIFEST"
test "$BUILD_FREEZE_LOG_EXPECTED_AUTHORIZED_FREEZE_MANIFEST_SHA256" = \
  "$PRESERVED_CLEAN_CANDIDATE_MANIFEST_SHA256"
test "$BUILD_FREEZE_LOG_EXPECTED_AUTHORIZED_FREEZE_BUILD_CONTROLLER" = \
  "$PRESERVED_CLEAN_CANDIDATE_BUILD_PROFILE_CONTROL_PLANE_VERSION"
CLEAN_CANDIDATE_MANIFEST="$(
  read_exact_build_freeze_log_binding clean_candidate_manifest
)"
CLEAN_CANDIDATE_MANIFEST_SHA256="$(
  read_exact_build_freeze_log_binding clean_candidate_manifest_sha256
)"
[[ "$CLEAN_CANDIDATE_MANIFEST_SHA256" =~ ^[0-9a-f]{64}$ ]]
EXECUTABLE="${CLEAN_CANDIDATE_MANIFEST%/*}/athena"
ENVIRONMENT_PROFILE="${ORION_CONTROL_PLANE}/frontier_pic_environment.sh"
test -x "$EXECUTABLE"
test -r "$ENVIRONMENT_PROFILE"
"${CONTROL_PLANE[@]}" revalidate_clean_candidate.py \
  --manifest "$CLEAN_CANDIDATE_MANIFEST" \
  --expected-manifest-sha256 "$CLEAN_CANDIDATE_MANIFEST_SHA256" \
  --expected-git-commit "$FULL_GIT_COMMIT" \
  --expected-receipt-control-plane-version "$VERSION"

run_reviewed_q011_materializer \
  candidate-only-policy-successor \
  --expected-git-commit "$FULL_GIT_COMMIT" \
  --baseline-policy "${PIC_ROOT}/policy/storage_policy.json" \
  --control-plane-version "$VERSION" \
  --storage-preflight-binding "$STORAGE_PREFLIGHT_BINDING" \
  --clean-candidate-manifest "$CLEAN_CANDIDATE_MANIFEST" \
  --executable "$EXECUTABLE" \
  --environment-profile "$ENVIRONMENT_PROFILE" \
  --output "$CANDIDATE_ONLY_POLICY"

EXPECTED_FINAL_ACTIVE_POLICY_SHA256="$(
  /usr/bin/sha256sum "$CANDIDATE_ONLY_POLICY" | /usr/bin/awk '{print $1}'
)"
[[ "$EXPECTED_FINAL_ACTIVE_POLICY_SHA256" =~ ^[0-9a-f]{64}$ ]]
require_no_policy_recovery_entries
require_no_staging_entries
/usr/bin/ps -u "$(/usr/bin/id -u)" -ww -o pid=,ppid=,lstart=,args=
"${SLURM_ENV[@]}" /usr/bin/squeue --clusters=frontier -u "$(/usr/bin/id -un)" \
  -h -o '%i|%j|%a|%q|%T|%Z|%o'
REVIEWED_NO_SCIENCE_OR_SAME_USER_PIC_MUTATOR='<yes only after reviewing both snapshots>'
test "$REVIEWED_NO_SCIENCE_OR_SAME_USER_PIC_MUTATOR" = yes
"${CONTROL_PLANE[@]}" promote_active_policy.py \
  --reviewed-policy "$CANDIDATE_ONLY_POLICY" \
  --replace-exact-authorized-clean-candidate-freeze \
  --expected-active-policy-sha256 "$EXPECTED_ACTIVE_POLICY_SHA256" \
  --expected-active-promotion-sha256 "$EXPECTED_ACTIVE_PROMOTION_SHA256"

FINAL_ACTIVE_POLICY_SHA256="$(
  /usr/bin/sha256sum "${PIC_ROOT}/policy/storage_policy.json" |
    /usr/bin/awk '{print $1}'
)"
FINAL_ACTIVE_PROMOTION_SHA256="$(
  /usr/bin/sha256sum "${PIC_ROOT}/policy/active_promotion.json" |
    /usr/bin/awk '{print $1}'
)"
test "$FINAL_ACTIVE_POLICY_SHA256" = "$EXPECTED_FINAL_ACTIVE_POLICY_SHA256"
[[ "$FINAL_ACTIVE_PROMOTION_SHA256" =~ ^[0-9a-f]{64}$ ]]
/usr/bin/cmp -s \
  "${PIC_ROOT}/policy/storage_policy.json" \
  "${PROJECT_HOME_POLICY_ROOT}/policy/storage_policy.json"
/usr/bin/cmp -s \
  "${PIC_ROOT}/policy/active_promotion.json" \
  "${PROJECT_HOME_POLICY_ROOT}/policy/active_promotion.json"
"${CONTROL_PLANE[@]}" revalidate_clean_candidate.py \
  --manifest "$CLEAN_CANDIDATE_MANIFEST" \
  --expected-manifest-sha256 "$CLEAN_CANDIDATE_MANIFEST_SHA256" \
  --expected-git-commit "$FULL_GIT_COMMIT" \
  --expected-receipt-control-plane-version "$VERSION"
require_no_policy_recovery_entries
require_no_staging_entries
"${CONTROL_PLANE[@]}" promote_active_policy.py \
  --verify-active-launch-prohibited-generation \
  --expected-control-plane-version "$VERSION" \
  --expected-active-policy-sha256 "$FINAL_ACTIVE_POLICY_SHA256" \
  --expected-active-promotion-sha256 "$FINAL_ACTIVE_PROMOTION_SHA256" \
  --expected-authorized-freeze-manifest "$CLEAN_CANDIDATE_MANIFEST" \
  --expected-authorized-freeze-manifest-sha256 "$CLEAN_CANDIDATE_MANIFEST_SHA256" \
  --expected-authorized-freeze-build-controller "$VERSION"
printf 'final_active_policy_sha256=%s\n' "$FINAL_ACTIVE_POLICY_SHA256"
printf 'final_active_promotion_sha256=%s\n' "$FINAL_ACTIVE_PROMOTION_SHA256"
printf 'final_authorized_clean_candidate_manifest=%s\n' "$CLEAN_CANDIDATE_MANIFEST"
printf 'final_authorized_clean_candidate_manifest_sha256=%s\n' \
  "$CLEAN_CANDIDATE_MANIFEST_SHA256"
)
```

Treat each mutation in the phases above as a durable checkpoint. If the session
stops during authenticated preflight capture, preserve every published evidence
or staging artifact and use only the exact commit-forward evidence recovery
block above. Do not delete evidence, recapture, or run a fresh primary block.
If the Orion binding was published, retain that immutable binding and resume
with its exact `PROBE_ID`; do not silently capture a replacement.

If the session stops after only one paired controller installation, use only
the exact one-sided recovery block above. If both installs completed, verify
both exact directories and their byte-identical inventories before resuming;
do not blindly rerun either install. If it stops after migration-policy
promotion, preserve the printed post-migration checkpoint, inspect the exact
active policy and promotion record, then use the worker-only resume block; do
not retry the one-use predecessor migration. If it stops during or after
`sbatch`, recover the accepted `BUILD_FREEZE_JOB_ID` through reviewed scheduler
inspection and its Slurm log. Do not resubmit unless that inspection proves no
job was accepted. Run the post-worker phase only after the exact eight-field
top-level `sacct` identity is `COMPLETED` with `ExitCode=0:0` and its log binds
the exact source commit.

A replacement submission requires a reviewed failure record. Either promotion
command may first recover an interrupted durable policy transaction while
holding both policy-parent locks. Every operator checkpoint must reject both
mirrored transaction markers, every reserved rollback anchor, and every
unexpected `.staging.*` or `.tmp-*` entry; never remove them manually. Do not
infer a committed state from only one active anchor. Candidate promotion is not
complete until the mandatory final verifier confirms the exact controller,
mirrored final policy/promotion hashes, exact newly authorized freeze, empty
science allowlist, expected source commit, and an empty recovery/staging
namespace.

The historical v2 launch procedure below is retained for audit only. It must
not be replayed. For one selected case, it required the operator to
capture and seal its `pre_manifest` attestation, write the six-field queue
snapshot, materialize one seed-timeout artifact and one selected-case config,
create its immutable manifest, capture and seal a fresh `pre_submit_wrapper`
attestation, and invoke only `submit_frontier_job.sh`. Do not materialize the
next selected-case config yet. The config materializer validates the repaired
generator bytes inside the frozen `source.tar`, rejects the failed v1 carrier,
and rejects every non-first case unless all exact preregistered predecessors
have completed reconciliation events and immutable descriptor checksums. Use
one UUID-specific handoff tree per case.

The installed reservation boundary repeats the authoritative checks under the
mirrored-ledger lock. It rejects Q011-equivalent authorization aliases unless
the authorization ID, campaign, test ID, job script, input deck and launch
contract all match one exact preregistered case. For predecessor descriptors,
it opens the immutable manifest, analyzer and helper snapshots with
`O_NOFOLLOW`, hashes stable descriptor-pinned bytes, rejects extra analysis
snapshot entries such as `__pycache__`, and executes verified source through
`compile()` without consulting cached bytecode. Source-local materialization is
only a review-artifact pre-screen; reservation through the installed controller
is the launch authority.

Set `PRIOR_CASE_CLOSURES` to the exact ordered mapping below. Omit it for
`ps_p0_1p00`; append one `--prior-case-closure` pair after each completed raw
analysis:

| selected case | required ordered `PRIOR_CASE_CLOSURES` values |
| --- | --- |
| `ps_p0_1p00` | none |
| `ps_p0_0p05` | `ps_p0_1p00=<submission-id>=<descriptor-sha256>` |
| `ps_p0_0p10` | `ps_p0_1p00=<submission-id>=<descriptor-sha256>`, then `ps_p0_0p05=<submission-id>=<descriptor-sha256>` |
| `ps_p0_0p20` | `ps_p0_1p00=<submission-id>=<descriptor-sha256>`, then `ps_p0_0p05=<submission-id>=<descriptor-sha256>`, then `ps_p0_0p10=<submission-id>=<descriptor-sha256>` |

Use this exact selected-case mapping:

| case | authorization ID | campaign |
| --- | --- | --- |
| `ps_p0_1p00` | `q011-section54-pressure-ps-p0-1p00-v2` | `q011_section54_pressure_ps_p0_1p00` |
| `ps_p0_0p05` | `q011-section54-pressure-ps-p0-0p05-v2` | `q011_section54_pressure_ps_p0_0p05` |
| `ps_p0_0p10` | `q011-section54-pressure-ps-p0-0p10-v2` | `q011_section54_pressure_ps_p0_0p10` |
| `ps_p0_0p20` | `q011-section54-pressure-ps-p0-0p20-v2` | `q011_section54_pressure_ps_p0_0p20` |

```bash
CASE='<one-preregistered-case-id>'
AUTHORIZATION_ID='<matching-preregistered-authorization-id>'
CAMPAIGN='<matching-preregistered-campaign>'
# Example for the first case: PRIOR_CASE_CLOSURES=()
# Example for the second case:
# PRIOR_CASE_CLOSURES=(--prior-case-closure "ps_p0_1p00=<submission-id>=<descriptor-sha256>")
SUBMISSION_ID="$("$PYTHON" -I -c 'import uuid; print(uuid.uuid4())')"
mkdir -p "${PIC_ROOT}/jobs/${CAMPAIGN}"
"$PYTHON" -I -c \
  'import os, sys; fd = os.open(sys.argv[1], os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW); os.fsync(fd); os.close(fd)' \
  "${PIC_ROOT}/jobs/${CAMPAIGN}"
HANDOFF_ROOT="${PIC_ROOT}/jobs/${CAMPAIGN}/${SUBMISSION_ID}"
mkdir "$HANDOFF_ROOT"

STAGING="$(
  "$PYTHON" "$ATTESTATION_HELPER" capture \
    --authorization-id "$AUTHORIZATION_ID" \
    --control-plane-version "$VERSION" \
    --phase pre_manifest
)"
PRE_MANIFEST_ATTESTATION_ROOT="$(
  "$PYTHON" "$ATTESTATION_HELPER" seal \
    --staging-dir "$STAGING" \
    --attest-reviewed
)"
PRE_MANIFEST_ATTESTATION="${PRE_MANIFEST_ATTESTATION_ROOT}/attestation.json"

(
  umask 077
  set -o noclobber
  squeue -u "$USER" -h -o '%i|%P|%q|%T|%j|%k' \
    > "${HANDOFF_ROOT}/queue_snapshot.txt"
)
chmod 0444 "${HANDOFF_ROOT}/queue_snapshot.txt"

"$PYTHON" -I -B /ccs/home/dfielding/athenak-pic/tst/publication/q011_section54_pressure_pilot_execution.py \
  seed-timeout-margin \
  --case-id "$CASE" \
  --environment-profile "$ENVIRONMENT_PROFILE" \
  --materialized-utc "$MATERIALIZED_UTC" \
  --expires-utc "$EXPIRES_UTC" \
  --output-root "${HANDOFF_ROOT}/timeout"

"$PYTHON" -I -B /ccs/home/dfielding/athenak-pic/tst/publication/q011_section54_pressure_pilot_execution.py \
  pre-submit-config \
  --case-id "$CASE" \
  --submission-id "$SUBMISSION_ID" \
  --control-plane-version "$VERSION" \
  --clean-candidate-manifest "$CLEAN_CANDIDATE_MANIFEST" \
  --executable "$EXECUTABLE" \
  --environment-profile "$ENVIRONMENT_PROFILE" \
  --pre-manifest-attestation "$PRE_MANIFEST_ATTESTATION" \
  --timeout-margin-artifact "${HANDOFF_ROOT}/timeout/timeout_margin.json" \
  --queue-snapshot "${HANDOFF_ROOT}/queue_snapshot.txt" \
  --site-policy-checked-utc "$LAST_PREFLIGHT_UTC" \
  "${PRIOR_CASE_CLOSURES[@]}" \
  --output-root "${HANDOFF_ROOT}/config"

MANIFEST="$(
  "${CONTROL_PLANE[@]}" create_pre_submit_manifest.py \
    --config "${HANDOFF_ROOT}/config/pre_submit_config.json"
)"

STAGING="$(
  "$PYTHON" "$ATTESTATION_HELPER" capture \
    --authorization-id "$AUTHORIZATION_ID" \
    --control-plane-version "$VERSION" \
    --phase pre_submit_wrapper
)"
PRE_SUBMIT_WRAPPER_ATTESTATION_ROOT="$(
  "$PYTHON" "$ATTESTATION_HELPER" seal \
    --staging-dir "$STAGING" \
    --attest-reviewed
)"
PRE_SUBMIT_WRAPPER_ATTESTATION="${PRE_SUBMIT_WRAPPER_ATTESTATION_ROOT}/attestation.json"

"${CONTROL_PLANE_DIR}/submit_frontier_job.sh" \
  "$MANIFEST" "$PRE_SUBMIT_WRAPPER_ATTESTATION"
```

After that job reaches a terminal scheduler state, run the installed
`reconcile_frontier_job.py` command below. Then invoke the snapshotted first
analysis script through the trusted isolated runner and retain the printed raw
descriptor SHA-256:

```bash
env -u PIC_F1_ANALYSIS_HELPER_FD /opt/cray/pe/python/3.11.7/bin/python3 -I -B \
  "${PIC_ROOT}/manifests/<campaign>/<submission-id>/snapshot/analysis/000-analyze_q011_section54_pressure_pilot_case.py" \
  --artifact-dir "${PIC_ROOT}/runs/<campaign>/<submission-id>" \
  --case-id '<case-id>'
```

Before starting the next case, require the reconciled terminal ledger event,
an absent `pending_submission.json`, zero active reservations, an empty
same-account PIC queue, and the immutable raw-case descriptor. This quiescent
terminal boundary is mandatory between every pair of cases.

After all four descriptors exist, publish the descriptor-verified aggregate
bundle once. The publisher accepts exactly four case directories and four
descriptor checksums. Create and sync the canonical publication root and its
fixed sibling acceptance authority once before publication; do not use a
Project Home bulk-artifact directory:

The first live acceptance-root provisioning attempt on 2026-06-03 failed
closed before `sbatch`: Orion inherited the parent setgid bit and left the exact
empty sibling root at mode `02700`, rather than the required `0700`. Preserve
that checkpoint. After the repaired helper is committed, pushed, validated by
the clean committed repair worker, and independently rereviewed from the exact
latest patch, run this incident-specific recovery once. It accepts only the
exact empty inherited-setgid checkpoint, requires an empty receipt-staging
namespace, mutates the retained no-follow descriptor to `0700`, syncs the child
and parent, and publishes one durable recovery receipt:

```bash
(
set -euo pipefail
test "${PIC_ROOT:-}" = /lustre/orion/ast207/proj-shared/dfielding/PIC
test "${PYTHON:-}" = /opt/cray/pe/python/3.11.7/bin/python3
SOURCE_REPO=/autofs/nccs-svm1_home2/dfielding/athenak-pic
SLURM_ENV=(/usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin SLURM_CLUSTERS=frontier)
GIT=(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin
  /usr/bin/git -c core.fsmonitor=false -c core.hooksPath=/dev/null
)
SOURCE_STATUS="$("${GIT[@]}" -C "$SOURCE_REPO" status --porcelain --untracked-files=all)"
test -z "$SOURCE_STATUS"
FULL_GIT_COMMIT="$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse HEAD)"
test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse origin/PIC)"
REMOTE_PIC_TIP="$(
  "${GIT[@]}" -C "$SOURCE_REPO" ls-remote --exit-code origin refs/heads/PIC |
    /usr/bin/awk '$2 == "refs/heads/PIC" {count += 1; tip = $1} END {if (count != 1) exit 1; print tip}'
)"
[[ "$REMOTE_PIC_TIP" =~ ^[0-9a-f]{40}$ ]]
test "$FULL_GIT_COMMIT" = "$REMOTE_PIC_TIP"
REPAIR_VALIDATION_COMMIT='<repair_validation_commit from the completed checkpoint>'
REPAIR_VALIDATION_JOB_ID='<repair_validation_job_id from the completed checkpoint>'
[[ "$REPAIR_VALIDATION_COMMIT" =~ ^[0-9a-f]{40}$ ]]
[[ "$REPAIR_VALIDATION_JOB_ID" =~ ^[0-9]+$ ]]
test "$REPAIR_VALIDATION_COMMIT" = "$FULL_GIT_COMMIT"
REPAIR_VALIDATION_RECORD="$(
  "${SLURM_ENV[@]}" /usr/bin/sacct -X --clusters=frontier -n \
    -j "$REPAIR_VALIDATION_JOB_ID" \
    --format=JobIDRaw,JobName,Account,Partition,QOS,State,ExitCode,WorkDir --parsable2 |
    /usr/bin/awk -F'|' -v job="$REPAIR_VALIDATION_JOB_ID" '
      $1 == job {count += 1; record = $1 FS $2 FS $3 FS $4 FS $5 FS $6 FS $7 FS $8}
      END {if (count != 1) exit 1; print record}
    '
)"
test "$REPAIR_VALIDATION_RECORD" = \
  "${REPAIR_VALIDATION_JOB_ID}|pic-q011-pressure-gate-validate|ast207|batch|debug|COMPLETED|0:0|${SOURCE_REPO}"
REPAIR_VALIDATION_LOG="${PIC_ROOT}/logs/slurm/pic-q011-pressure-gate-validate.${REPAIR_VALIDATION_JOB_ID}.log"
test -r "$REPAIR_VALIDATION_LOG"
REPAIR_VALIDATION_LOG_COMMIT="$(
  /usr/bin/sed -n 's/^source_commit=//p' "$REPAIR_VALIDATION_LOG"
)"
test "$REPAIR_VALIDATION_LOG_COMMIT" = "$FULL_GIT_COMMIT"
test "$(printf '%s\n' "$REPAIR_VALIDATION_LOG_COMMIT" | /usr/bin/wc -l)" -eq 1
/usr/bin/ps -u "$(/usr/bin/id -u)" -ww -o pid=,ppid=,lstart=,args=
"${SLURM_ENV[@]}" /usr/bin/squeue --clusters=frontier -u "$(/usr/bin/id -un)" \
  -h -o '%i|%j|%a|%q|%T|%Z|%o'
REVIEWED_NO_SCIENCE_OR_SAME_USER_PIC_MUTATOR='<yes only after reviewing both snapshots>'
test "$REVIEWED_NO_SCIENCE_OR_SAME_USER_PIC_MUTATOR" = yes
ACCEPTANCE_HELPER_RELATIVE=tst/publication/provision_q011_pressure_publication_acceptance_root.py
EXPECTED_ACCEPTANCE_HELPER_SHA256=431787450a6fe2a37a6ee1e1a37e1626444e2f6af3d7516117f820bc17963920
ACCEPTANCE_HELPER_SNAPSHOT="$(
  /usr/bin/mktemp -p /tmp q011-acceptance-helper.XXXXXX.py
)"
ACCEPTANCE_RECOVERY_RECEIPT="${PIC_ROOT}/policy/q011_pressure_publication_acceptance_root_recovery_1554766c-21e2-48b1-8cfe-b1e7e4e75aa2.json"
test ! -e "$ACCEPTANCE_RECOVERY_RECEIPT"
trap '/usr/bin/rm -f "$ACCEPTANCE_HELPER_SNAPSHOT"' EXIT
"${GIT[@]}" -C "$SOURCE_REPO" show \
  "${FULL_GIT_COMMIT}:${ACCEPTANCE_HELPER_RELATIVE}" \
  > "$ACCEPTANCE_HELPER_SNAPSHOT"
/usr/bin/chmod 0444 "$ACCEPTANCE_HELPER_SNAPSHOT"
ACCEPTANCE_HELPER_SHA256="$(
  /usr/bin/sha256sum "$ACCEPTANCE_HELPER_SNAPSHOT" | /usr/bin/awk '{print $1}'
)"
test "$ACCEPTANCE_HELPER_SHA256" = "$EXPECTED_ACCEPTANCE_HELPER_SHA256"
ACCEPTANCE_RECOVERY_PAYLOAD="$("$PYTHON" -I -B "$ACCEPTANCE_HELPER_SNAPSHOT" \
  --validated-source-commit "$FULL_GIT_COMMIT" \
  --expected-helper-sha256 "$ACCEPTANCE_HELPER_SHA256" \
  --recover-exact-empty-inherited-setgid-root
)"
ACCEPTANCE_RECOVERY_STAGING="$(
  /usr/bin/mktemp -p "${PIC_ROOT}/policy" .q011-acceptance-recovery.XXXXXX
)"
trap '/usr/bin/rm -f "$ACCEPTANCE_RECOVERY_STAGING" "$ACCEPTANCE_HELPER_SNAPSHOT"' EXIT
printf '%s\n' "$ACCEPTANCE_RECOVERY_PAYLOAD" > "$ACCEPTANCE_RECOVERY_STAGING"
"$PYTHON" -I -B "$ACCEPTANCE_HELPER_SNAPSHOT" \
  --validated-source-commit "$FULL_GIT_COMMIT" \
  --expected-helper-sha256 "$ACCEPTANCE_HELPER_SHA256" \
  --publish-recovery-receipt "$ACCEPTANCE_RECOVERY_STAGING"
"$PYTHON" -I -B "$ACCEPTANCE_HELPER_SNAPSHOT" \
  --validated-source-commit "$FULL_GIT_COMMIT" \
  --expected-helper-sha256 "$ACCEPTANCE_HELPER_SHA256" \
  --verify-recovery-receipt "$ACCEPTANCE_RECOVERY_RECEIPT"
/usr/bin/sha256sum "$ACCEPTANCE_RECOVERY_RECEIPT"
)
```

If the explicit recovery is interrupted after the retained inode reaches
`0700` but before its receipt is durably published, inspect that checkpoint and
rerun the same block with
`--reconcile-exact-empty-normalized-root` in place of
`--recover-exact-empty-inherited-setgid-root`. That mode accepts only the same
reviewed empty inode already at exact `0700`; it does not normalize a new path.
Before rerunning, require that `policy/` contains no retained
`.q011-acceptance-recovery.*` alias. A killed shell can leave an empty, partial,
or complete pre-link staging file after normalization. Stop for reviewed
inspection and removal of that exact orphan before retrying; do not let an
ordinary retry silently consume or overwrite it.
If the fixed recovery receipt already exists, verify it and do not republish it.
If that verification reports a receipt link-count drift after an interruption
inside the hard-link publication interval, inspect the fixed `policy/`
directory and invoke the same authenticated helper snapshot once with
`--reconcile-linked-recovery-receipt`. That explicit mode accepts only one
canonical read-only receipt with one matching retained
`.q011-acceptance-recovery.*` alias beneath the reviewed policy inode, removes
that alias, syncs the receipt and policy directory, and revalidates single-link
closure. Then run `--verify-recovery-receipt` again.

The reviewed recovery completed on 2026-06-04 without replacing the retained
inode. Aggregate worker job `4764415` then failed closed before exposing a
public bundle, receipt or analysis artifact because the offline analyzer
assumed one tuple order for the exact `mhd_w_bcc` variable inventory. Read-only
forensics found one uniform eight-field inventory across all twenty retained
snapshots, with no missing or extra variables. The committed fourteenth-pass
repair and worker validation passed. Replacement aggregate worker job `4764674`
then failed closed before exposing a public artifact because AthenaK emits
interior `dt`-scheduled products on the first committed step after each nominal
cadence while the analyzer still required exact nominal times. Read-only
forensics verified all `100` retained snapshot products and froze the exact
observed particle-header time, six-significant-digit mesh-header time and common
cycle tuple for every immutable case slot. Do not replay job `4764415` or
`4764674`. The fifteenth-pass repair was committed and clean worker job
`4766444` passed. Replacement aggregate job `4766456` then failed closed
before exposing a public artifact because Orion Lustre rejects nonzero
`renameat2` flags, including `RENAME_NOREPLACE`. Do not replay job `4766456`.
Submit one reviewed replacement aggregate worker only after the sixteenth-pass
preflighted guarded direct-final-name, hard-link no-replace, nondestructive
post-exposure reconciliation, and seal-before-final-guard-removal compatibility
repair is committed, pushed, clean-worker validated and independently
rereviewed from the exact latest patch. The replacement block below verifies
the existing durable recovery receipt; it must not republish it.

```bash
(
set -euo pipefail
test "${PIC_ROOT:-}" = /lustre/orion/ast207/proj-shared/dfielding/PIC
test "${PYTHON:-}" = /opt/cray/pe/python/3.11.7/bin/python3
SOURCE_REPO=/autofs/nccs-svm1_home2/dfielding/athenak-pic
SLURM_ENV=(/usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin SLURM_CLUSTERS=frontier)
GIT=(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin
  /usr/bin/git -c core.fsmonitor=false -c core.hooksPath=/dev/null
)
SOURCE_STATUS="$("${GIT[@]}" -C "$SOURCE_REPO" status --porcelain --untracked-files=all)"
test -z "$SOURCE_STATUS"
FULL_GIT_COMMIT="$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse HEAD)"
test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse origin/PIC)"
REMOTE_PIC_TIP="$(
  "${GIT[@]}" -C "$SOURCE_REPO" ls-remote --exit-code origin refs/heads/PIC |
    /usr/bin/awk '$2 == "refs/heads/PIC" {count += 1; tip = $1} END {if (count != 1) exit 1; print tip}'
)"
[[ "$REMOTE_PIC_TIP" =~ ^[0-9a-f]{40}$ ]]
test "$FULL_GIT_COMMIT" = "$REMOTE_PIC_TIP"
REPAIR_VALIDATION_COMMIT='<repair_validation_commit from the completed checkpoint>'
REPAIR_VALIDATION_JOB_ID='<repair_validation_job_id from the completed checkpoint>'
[[ "$REPAIR_VALIDATION_COMMIT" =~ ^[0-9a-f]{40}$ ]]
[[ "$REPAIR_VALIDATION_JOB_ID" =~ ^[0-9]+$ ]]
test "$REPAIR_VALIDATION_COMMIT" = "$FULL_GIT_COMMIT"
REPAIR_VALIDATION_RECORD="$(
  "${SLURM_ENV[@]}" /usr/bin/sacct -X --clusters=frontier -n \
    -j "$REPAIR_VALIDATION_JOB_ID" \
    --format=JobIDRaw,JobName,Account,Partition,QOS,State,ExitCode,WorkDir --parsable2 |
    /usr/bin/awk -F'|' -v job="$REPAIR_VALIDATION_JOB_ID" '
      $1 == job {count += 1; record = $1 FS $2 FS $3 FS $4 FS $5 FS $6 FS $7 FS $8}
      END {if (count != 1) exit 1; print record}
    '
)"
test "$REPAIR_VALIDATION_RECORD" = \
  "${REPAIR_VALIDATION_JOB_ID}|pic-q011-pressure-gate-validate|ast207|batch|debug|COMPLETED|0:0|${SOURCE_REPO}"
REPAIR_VALIDATION_LOG="${PIC_ROOT}/logs/slurm/pic-q011-pressure-gate-validate.${REPAIR_VALIDATION_JOB_ID}.log"
test -r "$REPAIR_VALIDATION_LOG"
REPAIR_VALIDATION_LOG_COMMIT="$(
  /usr/bin/sed -n 's/^source_commit=//p' "$REPAIR_VALIDATION_LOG"
)"
test "$REPAIR_VALIDATION_LOG_COMMIT" = "$FULL_GIT_COMMIT"
test "$(printf '%s\n' "$REPAIR_VALIDATION_LOG_COMMIT" | /usr/bin/wc -l)" -eq 1
/usr/bin/ps -u "$(/usr/bin/id -u)" -ww -o pid=,ppid=,lstart=,args=
"${SLURM_ENV[@]}" /usr/bin/squeue --clusters=frontier -u "$(/usr/bin/id -un)" \
  -h -o '%i|%j|%a|%q|%T|%Z|%o'
REVIEWED_NO_SCIENCE_OR_SAME_USER_PIC_MUTATOR='<yes only after reviewing both snapshots>'
test "$REVIEWED_NO_SCIENCE_OR_SAME_USER_PIC_MUTATOR" = yes
test -d "${PIC_ROOT}/publication" || /usr/bin/mkdir "${PIC_ROOT}/publication"
ACCEPTANCE_HELPER_RELATIVE=tst/publication/provision_q011_pressure_publication_acceptance_root.py
EXPECTED_ACCEPTANCE_HELPER_SHA256=431787450a6fe2a37a6ee1e1a37e1626444e2f6af3d7516117f820bc17963920
ACCEPTANCE_RECOVERY_VALIDATED_SOURCE_COMMIT=10bb501df0fa66d70f95f8983494a2156dd9ebd6
ACCEPTANCE_HELPER_SNAPSHOT="$(
  /usr/bin/mktemp -p /tmp q011-acceptance-helper.XXXXXX.py
)"
trap '/usr/bin/rm -f "$ACCEPTANCE_HELPER_SNAPSHOT"' EXIT
"${GIT[@]}" -C "$SOURCE_REPO" show \
  "${FULL_GIT_COMMIT}:${ACCEPTANCE_HELPER_RELATIVE}" \
  > "$ACCEPTANCE_HELPER_SNAPSHOT"
/usr/bin/chmod 0444 "$ACCEPTANCE_HELPER_SNAPSHOT"
ACCEPTANCE_HELPER_SHA256="$(
  /usr/bin/sha256sum "$ACCEPTANCE_HELPER_SNAPSHOT" | /usr/bin/awk '{print $1}'
)"
test "$ACCEPTANCE_HELPER_SHA256" = "$EXPECTED_ACCEPTANCE_HELPER_SHA256"
ACCEPTANCE_RECOVERY_RECEIPT="${PIC_ROOT}/policy/q011_pressure_publication_acceptance_root_recovery_1554766c-21e2-48b1-8cfe-b1e7e4e75aa2.json"
"$PYTHON" -I -B "$ACCEPTANCE_HELPER_SNAPSHOT" \
  --validated-source-commit "$ACCEPTANCE_RECOVERY_VALIDATED_SOURCE_COMMIT" \
  --expected-helper-sha256 "$ACCEPTANCE_HELPER_SHA256"
"$PYTHON" -I -B "$ACCEPTANCE_HELPER_SNAPSHOT" \
  --validated-source-commit "$ACCEPTANCE_RECOVERY_VALIDATED_SOURCE_COMMIT" \
  --expected-helper-sha256 "$ACCEPTANCE_HELPER_SHA256" \
  --verify-recovery-receipt "$ACCEPTANCE_RECOVERY_RECEIPT"
for directory in \
  "${PIC_ROOT}/publication" \
  "${PIC_ROOT}/publication_acceptance" \
  "${PIC_ROOT}"
do
  "$PYTHON" -I -c \
    'import os, sys; fd = os.open(sys.argv[1], os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW); os.fsync(fd); os.close(fd)' \
    "$directory"
done
cd "$SOURCE_REPO"
AGGREGATE_JOB_TOKEN="$("${SLURM_ENV[@]}" /usr/bin/sbatch --parsable --export=NIL \
  "${SOURCE_REPO}/tst/publication/frontier_q011_section54_pressure_pilot_publish_job.sh" \
  "$FULL_GIT_COMMIT")"
printf 'aggregate_job_token=%q\n' "$AGGREGATE_JOB_TOKEN"
AGGREGATE_JOB_ID="${AGGREGATE_JOB_TOKEN%;frontier}"
[[ "$AGGREGATE_JOB_ID" =~ ^[0-9]+$ ]]
test "$AGGREGATE_JOB_TOKEN" = "$AGGREGATE_JOB_ID" ||
  test "$AGGREGATE_JOB_TOKEN" = "${AGGREGATE_JOB_ID};frontier"
printf 'full_git_commit=%s\n' "$FULL_GIT_COMMIT"
printf 'aggregate_job_id=%s\n' "$AGGREGATE_JOB_ID"
)
```

After the aggregate worker completes, fill the two values printed above into
this packet-submission phase:

```bash
(
set -euo pipefail
test "${PIC_ROOT:-}" = /lustre/orion/ast207/proj-shared/dfielding/PIC
SOURCE_REPO=/autofs/nccs-svm1_home2/dfielding/athenak-pic
SLURM_ENV=(/usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin SLURM_CLUSTERS=frontier)
GIT=(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin
  /usr/bin/git -c core.fsmonitor=false -c core.hooksPath=/dev/null
)
FULL_GIT_COMMIT='<full_git_commit printed by aggregate submission>'
AGGREGATE_JOB_ID='<aggregate_job_id printed by aggregate submission>'
[[ "$FULL_GIT_COMMIT" =~ ^[0-9a-f]{40}$ ]]
[[ "$AGGREGATE_JOB_ID" =~ ^[0-9]+$ ]]
SOURCE_STATUS="$("${GIT[@]}" -C "$SOURCE_REPO" status --porcelain --untracked-files=all)"
test -z "$SOURCE_STATUS"
test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse HEAD)"
test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse origin/PIC)"
REMOTE_PIC_TIP="$(
  "${GIT[@]}" -C "$SOURCE_REPO" ls-remote --exit-code origin refs/heads/PIC |
    /usr/bin/awk '$2 == "refs/heads/PIC" {count += 1; tip = $1} END {if (count != 1) exit 1; print tip}'
)"
[[ "$REMOTE_PIC_TIP" =~ ^[0-9a-f]{40}$ ]]
test "$FULL_GIT_COMMIT" = "$REMOTE_PIC_TIP"
AGGREGATE_RECORD="$(
  "${SLURM_ENV[@]}" /usr/bin/sacct -X --clusters=frontier -n \
    -j "$AGGREGATE_JOB_ID" \
    --format=JobIDRaw,JobName,Account,Partition,QOS,State,ExitCode,WorkDir --parsable2 |
    /usr/bin/awk -F'|' -v job="$AGGREGATE_JOB_ID" '
      $1 == job {count += 1; record = $1 FS $2 FS $3 FS $4 FS $5 FS $6 FS $7 FS $8}
      END {if (count != 1) exit 1; print record}
    '
)"
test "$AGGREGATE_RECORD" = \
  "${AGGREGATE_JOB_ID}|pic-q011-pressure-publish|ast207|batch|debug|COMPLETED|0:0|${SOURCE_REPO}"
AGGREGATE_LOG="${PIC_ROOT}/logs/slurm/pic-q011-pressure-publish.${AGGREGATE_JOB_ID}.log"
test -r "$AGGREGATE_LOG"
AGGREGATE_LOG_COMMIT="$(
  /usr/bin/sed -n 's/^source_commit=//p' "$AGGREGATE_LOG"
)"
test "$AGGREGATE_LOG_COMMIT" = "$FULL_GIT_COMMIT"
test "$(printf '%s\n' "$AGGREGATE_LOG_COMMIT" | /usr/bin/wc -l)" -eq 1
/usr/bin/ps -u "$(/usr/bin/id -u)" -ww -o pid=,ppid=,lstart=,args=
"${SLURM_ENV[@]}" /usr/bin/squeue --clusters=frontier -u "$(/usr/bin/id -un)" \
  -h -o '%i|%j|%a|%q|%T|%Z|%o'
REVIEWED_NO_SCIENCE_OR_SAME_USER_PIC_MUTATOR='<yes only after reviewing both snapshots>'
test "$REVIEWED_NO_SCIENCE_OR_SAME_USER_PIC_MUTATOR" = yes
cd "$SOURCE_REPO"
REVIEW_PACKET_JOB_TOKEN="$("${SLURM_ENV[@]}" /usr/bin/sbatch --parsable --export=NIL \
  "${SOURCE_REPO}/tst/publication/frontier_q011_section54_pressure_pilot_review_packet_job.sh" \
  "$FULL_GIT_COMMIT")"
printf 'review_packet_job_token=%q\n' "$REVIEW_PACKET_JOB_TOKEN"
REVIEW_PACKET_JOB_ID="${REVIEW_PACKET_JOB_TOKEN%;frontier}"
[[ "$REVIEW_PACKET_JOB_ID" =~ ^[0-9]+$ ]]
test "$REVIEW_PACKET_JOB_TOKEN" = "$REVIEW_PACKET_JOB_ID" ||
  test "$REVIEW_PACKET_JOB_TOKEN" = "${REVIEW_PACKET_JOB_ID};frontier"
printf 'full_git_commit=%s\n' "$FULL_GIT_COMMIT"
printf 'review_packet_job_id=%s\n' "$REVIEW_PACKET_JOB_ID"
)
```

After the review-packet worker completes, fill its printed values into this
read-only verification phase:

```bash
(
set -euo pipefail
test "${PIC_ROOT:-}" = /lustre/orion/ast207/proj-shared/dfielding/PIC
test "${PYTHON:-}" = /opt/cray/pe/python/3.11.7/bin/python3
SOURCE_REPO=/autofs/nccs-svm1_home2/dfielding/athenak-pic
SLURM_ENV=(/usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin SLURM_CLUSTERS=frontier)
GIT=(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin
  /usr/bin/git -c core.fsmonitor=false -c core.hooksPath=/dev/null
)
FULL_GIT_COMMIT='<full_git_commit printed by packet submission>'
REVIEW_PACKET_JOB_ID='<review_packet_job_id printed by packet submission>'
[[ "$FULL_GIT_COMMIT" =~ ^[0-9a-f]{40}$ ]]
[[ "$REVIEW_PACKET_JOB_ID" =~ ^[0-9]+$ ]]
SOURCE_STATUS="$("${GIT[@]}" -C "$SOURCE_REPO" status --porcelain --untracked-files=all)"
test -z "$SOURCE_STATUS"
test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse HEAD)"
test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse origin/PIC)"
REMOTE_PIC_TIP="$(
  "${GIT[@]}" -C "$SOURCE_REPO" ls-remote --exit-code origin refs/heads/PIC |
    /usr/bin/awk '$2 == "refs/heads/PIC" {count += 1; tip = $1} END {if (count != 1) exit 1; print tip}'
)"
[[ "$REMOTE_PIC_TIP" =~ ^[0-9a-f]{40}$ ]]
test "$FULL_GIT_COMMIT" = "$REMOTE_PIC_TIP"
REVIEW_PACKET_RECORD="$(
  "${SLURM_ENV[@]}" /usr/bin/sacct -X --clusters=frontier -n \
    -j "$REVIEW_PACKET_JOB_ID" \
    --format=JobIDRaw,JobName,Account,Partition,QOS,State,ExitCode,WorkDir --parsable2 |
    /usr/bin/awk -F'|' -v job="$REVIEW_PACKET_JOB_ID" '
      $1 == job {count += 1; record = $1 FS $2 FS $3 FS $4 FS $5 FS $6 FS $7 FS $8}
      END {if (count != 1) exit 1; print record}
    '
)"
test "$REVIEW_PACKET_RECORD" = \
  "${REVIEW_PACKET_JOB_ID}|pic-q011-pressure-review|ast207|batch|debug|COMPLETED|0:0|${SOURCE_REPO}"
REVIEW_PACKET_LOG="${PIC_ROOT}/logs/slurm/pic-q011-pressure-review.${REVIEW_PACKET_JOB_ID}.log"
test -r "$REVIEW_PACKET_LOG"
REVIEW_PACKET_LOG_COMMIT="$(
  /usr/bin/sed -n 's/^source_commit=//p' "$REVIEW_PACKET_LOG"
)"
test "$REVIEW_PACKET_LOG_COMMIT" = "$FULL_GIT_COMMIT"
test "$(printf '%s\n' "$REVIEW_PACKET_LOG_COMMIT" | /usr/bin/wc -l)" -eq 1
SOURCE_STATUS="$("${GIT[@]}" -C "$SOURCE_REPO" status --porcelain --untracked-files=all)"
test -z "$SOURCE_STATUS"
test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse HEAD)"
test "$FULL_GIT_COMMIT" = "$("${GIT[@]}" -C "$SOURCE_REPO" rev-parse origin/PIC)"
REMOTE_PIC_TIP="$(
  "${GIT[@]}" -C "$SOURCE_REPO" ls-remote --exit-code origin refs/heads/PIC |
    /usr/bin/awk '$2 == "refs/heads/PIC" {count += 1; tip = $1} END {if (count != 1) exit 1; print tip}'
)"
[[ "$REMOTE_PIC_TIP" =~ ^[0-9a-f]{40}$ ]]
test "$FULL_GIT_COMMIT" = "$REMOTE_PIC_TIP"
SOURCE_SNAPSHOT="$(/usr/bin/mktemp -d /tmp/athenak-pic-final-review.XXXXXX)"
trap '/usr/bin/rm -rf -- "${SOURCE_SNAPSHOT:-}"' EXIT
"${GIT[@]}" -C "$SOURCE_REPO" archive --format=tar "$FULL_GIT_COMMIT" |
  /usr/bin/tar -xf - -C "$SOURCE_SNAPSHOT"
for SNAPSHOT_SOURCE in \
  publish_q011_section54_pressure_pilot_bundle.py \
  render_q011_section54_pressure_pilot_review_packet.py \
  readiness/plotting_environment_lock_candidate_2026-05-30.json; do
  test -f "${SOURCE_SNAPSHOT}/tst/publication/${SNAPSHOT_SOURCE}"
  test ! -L "${SOURCE_SNAPSHOT}/tst/publication/${SNAPSHOT_SOURCE}"
done
run_snapshot_python() {
  local script="${1:?usage: run_snapshot_python SCRIPT [ARG ...]}"
  shift
  "$PYTHON" -I -B -c \
    'import runpy, sys; root, script, *args = sys.argv[1:]; sys.path.insert(0, root); sys.argv = [script, *args]; runpy.run_path(script, run_name="__main__")' \
    "${SOURCE_SNAPSHOT}/tst/publication" \
    "${SOURCE_SNAPSHOT}/tst/publication/${script}" "$@"
}
run_snapshot_python publish_q011_section54_pressure_pilot_bundle.py \
  --verify-published-receipt \
  "${PIC_ROOT}/publication/q011_section54_pressure_pilot_bundle_receipt.json"
run_snapshot_plot_python() {
  local script="${1:?usage: run_snapshot_plot_python SCRIPT [ARG ...]}"
  shift
  "$PYTHON" -I -B -c \
    'import importlib.metadata, json, runpy, sys; root, script, package_root, lock_path, *args = sys.argv[1:]; sys.path.insert(0, package_root); lock = json.load(open(lock_path, encoding="utf-8")); actual_python = sys.version.split()[0]; expected_python = lock["python"]; actual_python == expected_python or sys.exit(f"plot Python drifted: {actual_python} != {expected_python}"); [(importlib.metadata.version(name) == version) or sys.exit(f"plot dependency drifted: {name}") for name, version in lock["dependencies"].items()]; sys.path.insert(0, root); sys.argv = [script, *args]; runpy.run_path(script, run_name="__main__")' \
    "${SOURCE_SNAPSHOT}/tst/publication" \
    "${SOURCE_SNAPSHOT}/tst/publication/${script}" \
    /autofs/nccs-svm1_home2/dfielding/.local/lib/python3.11/site-packages \
    "${SOURCE_SNAPSHOT}/tst/publication/readiness/plotting_environment_lock_candidate_2026-05-30.json" \
    "$@"
}
run_snapshot_plot_python render_q011_section54_pressure_pilot_review_packet.py \
  --verify-published-receipt \
  "${PIC_ROOT}/publication/q011_section54_pressure_pilot_review_packet_receipt.json"
/usr/bin/rm -rf -- "$SOURCE_SNAPSHOT"
trap - EXIT
)
```

The publisher reruns raw-case verification while retaining each case-root
descriptor, emits a separate immutable aggregate-analysis result and receipt,
and recomputes the engineering overlay analysis before and after publication.
The receipt binds the registered-execution preregistration, publisher,
aggregate analyzer and persisted result digests.

The review packet is an engineering-review artifact. Its worker and read-only
verification paths load the fixed user-site package directory explicitly and
reject any Python or dependency version that differs from
`plotting_environment_lock_candidate_2026-05-30.json`. That candidate lock does
not close the separate qualification plotting and external-export review gate.

Treat publication provisioning and each `sbatch` as durable checkpoints. If a
session stops after the acceptance root is provisioned, resume through the
helper's verification-only existing-root branch above; do not delete or
recreate it. The explicit inherited-setgid recovery option is reserved only for
the retained empty `02700` checkpoint described above. If the session stops
after aggregate submission, recover the one printed
`AGGREGATE_JOB_ID` through `sacct` and inspect its log and visible artifacts
before considering a reviewed replacement. Once the aggregate receipt verifies,
do not republish it. Apply the same rule to the one printed
`REVIEW_PACKET_JOB_ID`: recover and verify that exact packet job and receipt
before considering any reviewed replacement.

These four pilots remain engineering calibration only; they do not qualify
Section 5.4 physics.

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
job ID and scheduler-reported account match the mirrored reservation. Frontier
may report the configured `AST207` account canonically as `ast207`. A non-empty
`sacct` comment must exactly match `pic-reservation=<reservation-id>`. If
`sacct` omits the comment, reconciliation fails closed unless a fresh `scontrol`
query proves the exact job ID, account and reservation comment binding.
Caller-supplied terminal state, time and nodes are not accepted.

An immutable successor controller may reconcile a stranded terminal
`scheduler_job_id_received` marker from a prior installed generation only
through an immutable mirrored handoff created by the successor's installed
`terminal_recovery_handoff.py` tool. Supply the resulting Orion path as
`--terminal-recovery-handoff <path>`. The prior Orion and Project Home
controller copies and active policy snapshot must remain immutable and valid,
and the marker must exactly match the reserved ledger event. Paired reviewed
installation is the explicit successor-authorization boundary: do not invoke a
stale installed generation to create a recovery handoff.

If `scontrol` already purged a stranded held cancellation, stop for reviewed
recovery. The successor handoff tool has one explicit exceptional mode:

```bash
"${CONTROL_PLANE[@]}" terminal_recovery_handoff.py \
  --job-id <job-id> \
  --ledger-jsonl "${PIC_ROOT}/ledger/node_hours.jsonl" \
  --ledger-csv "${PIC_ROOT}/ledger/node_hours.csv" \
  --receipts-jsonl "${PIC_ROOT}/ledger/mirror_receipts.jsonl" \
  --mirror-jsonl "${PROJECT_HOME_MIRROR_ROOT}/ledger/node_hours.jsonl" \
  --authorize-purged-cancelled-zero-execution \
  --attest-reviewed-purged-reservation-job-binding
```

That mode fails closed unless fresh `sacct`, `scontrol` and `squeue` queries
prove the trusted wrapper job name, exact scheduler-record purge response,
empty live queue, `CANCELLED` state, zero elapsed seconds, zero allocated
nodes, empty accounting comment, authorized account, no start time, equal
submit and end times, and `0:0` exit code. The immutable mirrored handoff
freezes that accounting snapshot. It does not authorize broader
missing-`scontrol` recovery. A purged record cannot independently prove which
reservation produced the scheduler job ID, so the second flag records a
required reviewed operator attestation for that exact reservation-to-job
binding in the immutable mirrored handoff. Do not use the exceptional mode
without preserving and reviewing the submission incident evidence.

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
  --reservation-id '<reservation-id>'
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

## Q011 Stage-4 pressure-selection publication

This is a no-science publication boundary. Machine preparation seals only the
exact deterministic reanalysis and then stops at a human pressure-selection
gate. It does not register a science slice, authorize admission smoke, submit a
job, or grant a science claim. A candidate receipt can exist only after the
sealer consumes an independently supplied decision record that explicitly
binds the sealed reanalysis. The software authenticates that record's exact
reviewed content, private location, read-only mode, and strict post-reanalysis
chronology; it does not cryptographically authenticate human authorship. The
operational human stop and independent review are mandatory authority
controls.

Run this sequence only after the Stage-4 source commit has passed the exact
clean-worker wrapper. The active candidate-only generation must remain:

- controller `930a04d1d39c873ea49abfcf500069011f6d5759240a8f5c5b3341a6d243b246`;
- policy SHA-256 `aeab7e4ef92c7cbbd5b84fcd139c046f2a96f21b4f280fa6dd89c80deca981d1`;
- promotion SHA-256 `ef11cb301ec4917cd32aaca8af56e4f2d753367682c6ba28fc048c004904613e`;
- clean-candidate manifest SHA-256
  `ea5f295096b04d7e5f338873c2a29213f173677568fd33220e1e34ea239739e4`;
- empty `registered_science_slices`, closed admission smoke, no pending
  submission, and no incomplete manual-accounting marker.

Machine reanalysis, human sealing, and publication are separate mutations. All
must execute from the same read-only extracted Git archive for the exact
committed Stage-4 source. Run the machine-reanalysis command first, retain its
output, and stop for human review:

```bash
set -euo pipefail
umask 077

REPO_ROOT=/ccs/home/dfielding/athenak-pic
PIC_ROOT=/lustre/orion/ast207/proj-shared/dfielding/PIC
PROJECT_HOME_ROOT=/autofs/nccs-svm1_proj/ast207/proj-shared/PIC
PROJECT_LEDGER_ROOT=/ccs/proj/ast207/proj-shared/PIC
PYTHON=/opt/cray/pe/python/3.11.7/bin/python3
SLURM_ENV=(/usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin SLURM_CLUSTERS=frontier)

trusted_git() {
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
    /usr/bin/git --no-replace-objects \
    -c core.fsmonitor=false -c core.hooksPath=/dev/null "$@"
}

cd "$REPO_ROOT"
FULL_GIT_COMMIT=$(trusted_git rev-parse HEAD)
test "$FULL_GIT_COMMIT" = "$(trusted_git rev-parse origin/PIC)"
test "$FULL_GIT_COMMIT" = "$(
  trusted_git ls-remote origin refs/heads/PIC | /usr/bin/awk '{print $1}'
)"
test -z "$(trusted_git status --porcelain --untracked-files=all)"

test "$(/usr/bin/sha256sum "$PIC_ROOT/policy/storage_policy.json" | /usr/bin/awk '{print $1}')" \
  = aeab7e4ef92c7cbbd5b84fcd139c046f2a96f21b4f280fa6dd89c80deca981d1
test "$(/usr/bin/sha256sum "$PIC_ROOT/policy/active_promotion.json" | /usr/bin/awk '{print $1}')" \
  = ef11cb301ec4917cd32aaca8af56e4f2d753367682c6ba28fc048c004904613e
test "$(/usr/bin/sha256sum "$PIC_ROOT/clean_candidates/83dce7b7-0b03-4be2-b6da-17bb211d1fd4/clean_candidate_manifest.json" | /usr/bin/awk '{print $1}')" \
  = ea5f295096b04d7e5f338873c2a29213f173677568fd33220e1e34ea239739e4
PRE_RECOVERY_QUEUE_SNAPSHOT=$(
  "${SLURM_ENV[@]}" /usr/bin/squeue --clusters=frontier -u "$(/usr/bin/id -un)" \
    -h -o '%i|%j|%a|%q|%T|%Z|%o'
)
printf '%s\n' "$PRE_RECOVERY_QUEUE_SNAPSHOT"
test -z "$(printf '%s\n' "$PRE_RECOVERY_QUEUE_SNAPSHOT" |
  /usr/bin/grep -Ei 'pic|q011|athena' || true)"
test ! -e "$PIC_ROOT/ledger/pending_submission.json"
test ! -e "$PIC_ROOT/ledger/pending_manual_accounting.json"
test ! -e "$PROJECT_LEDGER_ROOT/ledger/pending_manual_accounting.json"

SOURCE_ARCHIVE=$(/usr/bin/mktemp "${TMPDIR:-/tmp}/pic-q011-stage4-source.XXXXXXXX.tar")
SOURCE_SNAPSHOT=$(/usr/bin/mktemp -d "${TMPDIR:-/tmp}/pic-q011-stage4-snapshot.XXXXXXXX")
trusted_git archive --format=tar "$FULL_GIT_COMMIT" > "$SOURCE_ARCHIVE"
/usr/bin/tar -xf "$SOURCE_ARCHIVE" -C "$SOURCE_SNAPSHOT"
/usr/bin/find "$SOURCE_SNAPSHOT" -type f -exec /usr/bin/chmod 0400 {} +
/usr/bin/find "$SOURCE_SNAPSHOT" -type d -exec /usr/bin/chmod 0500 {} +
/usr/bin/chmod 0400 "$SOURCE_ARCHIVE"

cd "$SOURCE_SNAPSHOT"

run_stage4() {
  /usr/bin/env -i \
    HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin TMPDIR=/tmp \
    PIC_PRESSURE_PUBLICATION_SOURCE_ARCHIVE_PATH="$SOURCE_ARCHIVE" \
    PIC_PRESSURE_PUBLICATION_SOURCE_SNAPSHOT_ROOT="$SOURCE_SNAPSHOT" \
    "$PYTHON" -I -B -S -c \
    'import runpy, sys; root, site, *args = sys.argv[1:]; sys.path[:0] = [root, root + "/tst/publication/frontier_control_plane", site]; sys.argv = ["publish_q011_section54_pressure_selection.py", *args]; runpy.run_path(root + "/tst/publication/publish_q011_section54_pressure_selection.py", run_name="__main__")' \
    "$SOURCE_SNAPSHOT" \
    /opt/cray/pe/python/3.11.7/lib/python3.11/site-packages \
    "$@"
}

PRIVATE_ROOT_RECOVERY_JSON=$(
  /usr/bin/mktemp "${TMPDIR:-/tmp}/pic-q011-stage4-private-root-recovery.XXXXXXXX.json"
)
REANALYSIS_JSON=$(
  /usr/bin/mktemp "${TMPDIR:-/tmp}/pic-q011-stage4-reanalysis.XXXXXXXX.json"
)
/usr/bin/chmod 0600 "$PRIVATE_ROOT_RECOVERY_JSON" "$REANALYSIS_JSON"

/usr/bin/ps -u "$(/usr/bin/id -u)" -ww -o pid=,ppid=,lstart=,args=
"${SLURM_ENV[@]}" /usr/bin/squeue --clusters=frontier -u "$(/usr/bin/id -un)" \
  -h -o '%i|%j|%a|%q|%T|%Z|%o'
REVIEWED_NO_SCIENCE_OR_SAME_USER_PIC_MUTATOR='<yes only after reviewing both snapshots>'
test "$REVIEWED_NO_SCIENCE_OR_SAME_USER_PIC_MUTATOR" = yes
PRESSURE_GATE_ATTESTATION_ROOT="$PIC_ROOT/pressure_gate_attestations"
test -d "$PRESSURE_GATE_ATTESTATION_ROOT"
test "$(/usr/bin/stat -c '%d|%i|%u|%g' "$PRESSURE_GATE_ATTESTATION_ROOT")" \
  = '135357496|720587399016566207|18664|31114'
test -z "$(/usr/bin/find "$PRESSURE_GATE_ATTESTATION_ROOT" -mindepth 1 -print -quit)"
test ! -e "$PIC_ROOT/pressure_gate_human_decisions"
test ! -e "$PIC_ROOT/pressure_gate_candidates"
test ! -e "$PIC_ROOT/publication/q011_section54_pressure_selection_receipt.json"
case "$(/usr/bin/stat -c '%a' "$PRESSURE_GATE_ATTESTATION_ROOT")" in
  2700)
    PRIVATE_ROOT_COMMAND=recover-pressure-gate-attestation-root
    PRIVATE_ROOT_ACTION=recovered_exact_empty_inherited_setgid_root
    ;;
  700)
    PRIVATE_ROOT_COMMAND=reconcile-pressure-gate-attestation-root
    PRIVATE_ROOT_ACTION=reconciled_exact_empty_normalized_root
    ;;
  *)
    exit 1
    ;;
esac
run_stage4 \
  "$PRIVATE_ROOT_COMMAND" \
  --expected-git-commit "$FULL_GIT_COMMIT" \
  --authorized-pic-root "$PIC_ROOT" \
  --authorized-project-home-root "$PROJECT_HOME_ROOT" \
  | /usr/bin/tee "$PRIVATE_ROOT_RECOVERY_JSON"
test "$(/usr/bin/jq -r '.action' "$PRIVATE_ROOT_RECOVERY_JSON")" \
  = "$PRIVATE_ROOT_ACTION"
test "$(/usr/bin/jq -r '.authority' "$PRIVATE_ROOT_RECOVERY_JSON")" \
  = none
test "$(/usr/bin/jq -r '.qualification_effect' "$PRIVATE_ROOT_RECOVERY_JSON")" \
  = none_no_reanalysis_no_human_selection_no_science_or_launch_authority
test "$(/usr/bin/stat -c '%d|%i|%u|%g|%a' "$PRESSURE_GATE_ATTESTATION_ROOT")" \
  = '135357496|720587399016566207|18664|31114|700'
test -z "$(/usr/bin/find "$PRESSURE_GATE_ATTESTATION_ROOT" -mindepth 1 -print -quit)"
test ! -e "$PIC_ROOT/pressure_gate_human_decisions"
test ! -e "$PIC_ROOT/pressure_gate_candidates"
test ! -e "$PIC_ROOT/publication/q011_section54_pressure_selection_receipt.json"

run_stage4 \
  prepare-reanalysis \
  --reanalysis-operator-id codex \
  --expected-git-commit "$FULL_GIT_COMMIT" \
  --authorized-pic-root "$PIC_ROOT" \
  --authorized-project-home-root "$PROJECT_HOME_ROOT" \
  | /usr/bin/tee "$REANALYSIS_JSON"

POST_REANALYSIS_QUEUE_SNAPSHOT=$(
  "${SLURM_ENV[@]}" /usr/bin/squeue --clusters=frontier -u "$(/usr/bin/id -un)" \
    -h -o '%i|%j|%a|%q|%T|%Z|%o'
)
printf '%s\n' "$POST_REANALYSIS_QUEUE_SNAPSHOT"
test -z "$(printf '%s\n' "$POST_REANALYSIS_QUEUE_SNAPSHOT" |
  /usr/bin/grep -Ei 'pic|q011|athena' || true)"
test ! -e "$PIC_ROOT/ledger/pending_submission.json"
test ! -e "$PIC_ROOT/ledger/pending_manual_accounting.json"
test ! -e "$PROJECT_LEDGER_ROOT/ledger/pending_manual_accounting.json"
test "$(/usr/bin/jq -r '.candidate_pressure_selection_receipt // empty' "$REANALYSIS_JSON")" = ""
test "$(/usr/bin/jq -r '.reviewer_attestation // empty' "$REANALYSIS_JSON")" = ""
printf 'private_root_recovery_json=%s\n' "$PRIVATE_ROOT_RECOVERY_JSON"
printf 'reanalysis_json=%s\n' "$REANALYSIS_JSON"
```

The private-root command immediately before `prepare-reanalysis` is the
incident-specific recovery boundary for the retained empty
`pressure_gate_attestations` inode that inherited exact mode `02700` from the
setgid PIC root. It changes only that retained inode to `0700`, emits
`authority=none`, and cannot prepare reanalysis, create a human decision or
candidate, publish a receipt, or grant science or launch authority. Ordinary
Stage-4 commands continue to reject a pre-existing `02700` private root.
Fresh private roots created by those commands remove inherited ACLs and
normalize only the inode created by that same invocation. If recovery is
interrupted after normalization, inspect the exact empty retained inode and
rerun only the explicit reconciliation command shown by the mode branch; never
use shell `chmod`, delete/recreate the root, or pass recovery behavior into
`prepare-reanalysis`. The publication locks coordinate cooperating Stage-4
tools; they cannot exclude a non-cooperating process running as the same Unix
identity. The reviewed process and scheduler snapshots immediately before
recovery are therefore a mandatory operational authority boundary.

**Stop here for the human pressure-selection gate.** The reviewer must inspect
the immutable pressure-review packet and the exact sealed reanalysis, then
independently create one canonical read-only decision JSON as a direct child
of the private `human_decision_root` printed above. It must bind the printed
`authoritative_reanalysis_attestation`, include `"schema_version": 1` and
`"record_type": "q011_section54_pressure_selection_human_decision"`, state
reviewer `dfielding`, select `ps_p0_1p00` with `problem_ps_p0=1.0`, preserve
the reviewed rationale exactly, include a `reviewed_utc` strictly later than
the sealed reanalysis, and include the exact `reviewer_statement` printed
under `required_human_confirmation`. It must also bind the printed
`stage4_preparation_attestation`, which durably binds the reanalysis to the
authenticated Stage-4 publisher archive. The reanalysis operator must not
create or pre-populate this decision. The sealer authenticates the decision
record, not the human author's identity.

Only after that independently supplied decision exists:

```bash
HUMAN_DECISION="${HUMAN_DECISION:?set to the reviewed direct-child decision JSON}"
test "$(stat -c '%a' "$HUMAN_DECISION")" = 400

run_stage4 \
  seal-human-selection \
  --human-decision "$HUMAN_DECISION" \
  --expected-git-commit "$FULL_GIT_COMMIT" \
  --authorized-pic-root "$PIC_ROOT" \
  --authorized-project-home-root "$PROJECT_HOME_ROOT" \
  | /usr/bin/tee /tmp/pic-q011-stage4-human-seal.json

CANDIDATE_RECEIPT=$(
  /usr/bin/jq -r '.candidate_pressure_selection_receipt.path' \
    /tmp/pic-q011-stage4-human-seal.json
)
CANDIDATE_AUTHORIZATION=$(
  /usr/bin/jq -r '.candidate_publication_authorization.path' \
    /tmp/pic-q011-stage4-human-seal.json
)
test -f "$CANDIDATE_RECEIPT"
test -f "$CANDIDATE_AUTHORIZATION"

if run_stage4 \
  publish-selection \
  --receipt "$CANDIDATE_RECEIPT" \
  --candidate-authorization "$CANDIDATE_AUTHORIZATION" \
  --controller-operator-id codex \
  --expected-git-commit "$FULL_GIT_COMMIT" \
  --authorized-pic-root "$PIC_ROOT" \
  --authorized-project-home-root "$PROJECT_HOME_ROOT" \
  > /tmp/pic-q011-stage4-publish.json \
  2> /tmp/pic-q011-stage4-publish-failure.json; then
  :
else
  /usr/bin/jq -e '.status == "failed_closed"' \
    /tmp/pic-q011-stage4-publish-failure.json
  /usr/bin/jq -e '.recovery != null' \
    /tmp/pic-q011-stage4-publish-failure.json
  CONTROLLER_PATH=$(/usr/bin/jq -r '.recovery.controller_state_attestation.path' /tmp/pic-q011-stage4-publish-failure.json)
  CONTROLLER_SHA256=$(/usr/bin/jq -r '.recovery.controller_state_attestation.sha256' /tmp/pic-q011-stage4-publish-failure.json)
  run_stage4 \
    reconcile-published \
    --controller-state-attestation "$CONTROLLER_PATH" \
    --expected-controller-state-sha256 "$CONTROLLER_SHA256" \
    --expected-git-commit "$FULL_GIT_COMMIT" \
    --authorized-pic-root "$PIC_ROOT" \
    --authorized-project-home-root "$PROJECT_HOME_ROOT" \
    | /usr/bin/tee /tmp/pic-q011-stage4-reconcile.json
fi

run_stage4 \
  verify-published \
  --receipt "$PIC_ROOT/publication/q011_section54_pressure_selection_receipt.json" \
  --authorized-pic-root "$PIC_ROOT" \
  --authorized-project-home-root "$PROJECT_HOME_ROOT" \
  | /usr/bin/tee /tmp/pic-q011-stage4-verify.json
```

Keep `SOURCE_ARCHIVE`, `SOURCE_SNAPSHOT`, the reanalysis output, human decision,
human-seal output, publication output, and any reconciliation output until the
receipt, success seal, recovery-guard state, controller-state attestation,
paired live active state, and absent publication guard have been independently
reviewed. If publication fails, inspect its canonical failure JSON. Do not
republish or delete artifacts. A non-null `.recovery` must be reconciled in the
same shell from the same authenticated source snapshot, as the command block
above does; a null `.recovery` is a failed-closed stop requiring review.
The sealed Stage-4 preparation attestation, human decision, candidate
publication authorization, and controller-state attestation form the durable
chain requiring preparation, human sealing, and publication to use the same
replacement-ref-disabled authenticated publisher archive.
