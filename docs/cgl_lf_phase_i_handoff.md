# CGL-LF Phase I Production Handoff

Checkpoint refreshed: `2026-06-02T23:58:06Z`

## Read This First

This is the durable resume point for the corrected-E03 MKS24 Stage I
production campaign. Use
`docs/cgl_lf_weak_guide_manuscript_plan.md` as the primary scientific guiding
document, `docs/cgl_lf_mks24_reproduction_implementation_plan.md` as the
detailed implementation and evidence record, and
`docs/cgl_lf_mks24_stage_i_protocol_review.md` as the stop-the-line production
protocol. This file records the newer operational boundary that had not yet
been promoted into those longer historical records.

Do not submit a new production job immediately. F-101 appeared canonical-only
outside the retained-companion workflow and has now been adopted exactly once
with the retained companion's explicit `adopt-legacy-canonical` action. The
adoption truthfully attests the present canonical bytes without claiming
observation of the original publication transition. Do not rerun adoption.
The separate Stage I helper hardening checkpoint is also reviewed, pushed,
archived, and cataloged. The next action is a fresh empty-queue and free-lock
boundary, authenticated hardened reconciliation, and an independently audited
bounded `R02/s19_rankio_t7p25_t7p5` readiness packet.

## Repository Boundary

Work in:

```text
/autofs/nccs-svm1_home2/dfielding/athenak-df
```

The branch at this checkpoint is `feature/cgl-landau-fluid`. Its durable
retained-utility checkpoint is:

```text
89ba4143c448c26fd8a66111d0409b4bcd3eea89
Add retained Stage I recost lifecycle utility
```

It is pushed and archived as
`source-archives/athenak-feature-cgl-through-89ba4143c.bundle` with SHA-256
`16f5b14610b47e1e3ca0f695bf9e5a94fba6d450e4dfe25357f0cde54fcbf597`.

The separate post-adoption Stage I helper hardening checkpoint is:

```text
f675bd677cb582a46bbda8b5f55c335e6f259dd3
Harden Stage I lock and scheduler routing
```

It is pushed and archived as
`source-archives/athenak-feature-cgl-through-f675bd677.bundle` with SHA-256
`eb32cd98a97ff57311d1a783ca6b80d3c22192370c09e85fc72843fc59ae7c42`.
Its Stage I helper SHA-256 is
`54ec671bb45aa27735a174d40b4b2e6009070716346ea09699bbe62421bbfada`.

The shared production root is:

```text
/lustre/orion/ast207/proj-shared/dfielding/CGL
```

Do not modify or stage unrelated working-tree state. At checkpoint capture,
the repository also contained:

```text
 M kokkos
?? scripts/frontier/rollback_codex_node_local_runtime.sh
```

Those paths are outside this campaign checkpoint. Treat them as user or
concurrent-process work unless the user explicitly assigns them.

## Completed Through Accepted R02 s18

Corrected-E03 `R02/s18_rankio_t7_t7p25` completed as Slurm job `4754394` and
is formally inspected, independently validated, recorded `accepted`, and
reconciled. The exact scheduler row is:

```text
4754394|cgl_mks24_E03_forcing_policy_R02_s18_rankio_t7_t7p25|COMPLETED|0:0|1|2340|2026-06-02T00:26:43|2026-06-02T01:05:59
```

The segment reached exact `t = 7.25`, used `2340` elapsed seconds
(`0.650000` node-hours), and brings the accepted corrected-E03 prefix to
`60420 / 3600 = 16.783333333333333` exact node-hours
(`16.783334` displayed). It retained one complete eight-sibling snapshot
group and terminal eight-sibling restart group, zero strict LF failure
counters, finite synchronized histories, terminal `lf_hwproj = 269865939181`,
segment forcing-work residual `6.449514054893265e-12`, and accepted-prefix
residual `3.377143242247066e-12`.

Retained evidence:

| Artifact | SHA-256 |
| --- | --- |
| `runs/mks24-stage-i/E03-forcing-policy/R02/s18_rankio_t7_t7p25/manifest/segment_inspection.json` | `400cec346d2fe9be2fd15c8de0b44c9f1edf24e7c2903e22dd6949960ca6d636` |
| `runs/mks24-stage-i/E03-forcing-policy/R02/s18_rankio_t7_t7p25/manifest/prepared_run.json` after record | `0e049354e66b4a2e4d0b7ae992ad8284e37b0687952bdf041219116e8e419a27` |
| `accounting/4754394.stage_i.independent_validation.json` | `5760966f1ccaf2c88f4c14307de281864e1aa93cdee01a33ff513feac3662905` |
| `accounting/4754394.stage_i.sacct.txt` | `1478076dc09c898a3ea7ddb8280f36d8a3c781434317bf08fe3272dff92c4971` |

Recorded reconciliation closed with `19/19/19`
ledger rows/manifests/reservations, no active reservation, no transaction,
and `issues = []`.

## Current F-101 Boundary

The retained F-101 recost generator is:

```text
accounting/utilities/generate_cgl_r02_t7p25_recost.py
SHA-256 6f5db0b32085a90f644d810e6ec7bdd00d838617b9cb22711c010aa4837fad72
mode 0755, one link
```

It generated a reviewed staged-only F-101 recost artifact:

```text
accounting/mks24_stage_i_E03_forcing_policy_R02_t7p25_recost_evidence.json.staged
SHA-256 df208829aec85c88fcc2caa757a21ad3ff3372aeb82628cbfe7bd41851faed9a
mode 0644, one link
```

At the original checkpoint capture, the canonical counterpart was absent:

```text
accounting/mks24_stage_i_E03_forcing_policy_R02_t7p25_recost_evidence.json
```

The user queue was empty and
`accounting/mks24_stage_i_E03_forcing_policy_transactions` contained no
entries at checkpoint capture. A later read-only documentation-audit probe
found the canonical Stage I lock persistently busy while the queue,
transaction directory, and staged-only recost namespace remained unchanged.
Treat lock availability as a dynamic hard gate: identify or wait out the
holder and require a fresh free-lock boundary before any promotion attempt.
The lock is the retained root-level advisory-flock file
`.mks24_stage_i_E03_forcing_policy.lock`; its existence alone is not a stale
lock signal. Read-only probes from the resumed session found no visible local
holder, so a distributed Lustre client remains possible. Unrelated PIC jobs
have also cycled through the user queue. Never bypass the flock or overlap a
queue-visible user job.

A resumed read-only probe later found an externally changed canonical-only
namespace:

```text
accounting/mks24_stage_i_E03_forcing_policy_R02_t7p25_recost_evidence.json
SHA-256 df208829aec85c88fcc2caa757a21ad3ff3372aeb82628cbfe7bd41851faed9a
mode 0644, one link
ctime 2026-06-02T22:57:18Z
```

The `.staged` name is now absent. The canonical bytes exactly match the
reviewed staged SHA-256. No retained-companion recost transaction directory,
forensic directory, or publication audit existed when this state was
discovered. The original publication transition was not observed by the
retained companion and must not be described retroactively as companion
publication. Treat this as a legacy canonical-only adoption boundary.

That legacy boundary has now been adopted exactly once under the retained
companion. A fresh pre-mutation probe found an empty live user queue and a
free root flock with profile `mode0644|links1|uid18664|regularTrue`. The
retained companion then ran `adopt-legacy-canonical`, and a separate
`verify-promoted-recost` invocation passed against the live root. The durable
adoption audit is:

```text
accounting/mks24_stage_i_E03_forcing_policy_R02_t7p25_recost_evidence.json.publication_audit.json
SHA-256 9bbfc6ee1007b0e17a4408444acb1bcdb2f59cbde9183569d506fe2606c1837d
mode 0644, one link
record_type legacy-canonical-adoption
adopted_utc 2026-06-02T23:35:28+00:00
transaction_id 2026-06-02T233522+0000-1b99386ad14d4eec84302c5005ea98f3
original_publication_transition_observed false
original_publication_method unknown
original_publisher unknown
```

Its independent retained forensic copy is:

```text
accounting/mks24_stage_i_E03_forcing_policy_recost_forensics/
  mks24_stage_i_E03_forcing_policy_R02_t7p25_recost_evidence.json/
  2026-06-02T233522+0000-1b99386ad14d4eec84302c5005ea98f3.mks24_stage_i_E03_forcing_policy_R02_t7p25_recost_evidence.json.forensic
SHA-256 df208829aec85c88fcc2caa757a21ad3ff3372aeb82628cbfe7bd41851faed9a
mode 0444, one link
```

Post-adoption read-only inspection found
`accounting/mks24_stage_i_E03_forcing_policy_recost_transactions` present
and empty, the staged twin absent, the canonical artifact unchanged at mode
`0644` and one link, the live user queue empty, and the root flock free with
the strict profile above.

The artifact has been independently audited against the retained
generator, accepted ledger, `s18` evidence, and arithmetic. Its reviewed
authorization recommendation is only:

```text
R02/s19_rankio_t7p25_t7p5
one node
Slurm walltime 01:05:00
Athena timeout 00:55:00
one-segment threshold 2700 seconds
600-second Athena and Slurm guards
```

The recost uses local `2340`, selected `2559`, projection
`829.1011116666666`, and margin `70.89888833333339` node-hours. This is a
recommendation inside retained evidence, not permission to prepare or submit
`s19` before the fresh live queue/free-lock boundary, authenticated hardened
reconciliation, and independently audited readiness packet close.

## Interrupted Promotion Review

The first reconstructed transient promoter had SHA-256:

```text
fdce32b3179f3b9128333e8330123a73be55aa1ffacdd8f68b102ccfbae1a20b
```

An independent read-only review correctly rejected it before execution:

1. It authenticated the live helper with mode `0755`, but the live helper is
   mode `0644`.
2. Its final queue check was not immediately before `os.link()`, leaving a
   publication race after expensive reconciliation and stale-state scans.

A transient corrected reconstruction was compiled with SHA-256:

```text
9c28d485e80896d8df9e5d8c3d0ca99a356410ca2f9af87e56c14b4bffd7d3a4
```

That revision changed the helper profile to `0644` and inserted a fresh
`require_empty_queue()` directly before the publication link. It had not
completed independent re-review or executed when checkpointing began. Its
`/tmp` file subsequently disappeared. Do not infer executable bytes from the
hash and do not publish with an unaudited reconstruction.

## Retained Routine-Lifecycle Utility

The transient promoter reconstruction has been superseded by a retained,
reviewed companion:

```text
scripts/frontier/cgl_lf_stage_i_checkpoint.py
SHA-256 10156515c4bcbfdcf57a2fe54220c2a80f0477f1a7bddbde322b9433946615c2
mode 0755, one link
```

Those exact bytes were independently reviewed with no findings, committed as
`89ba4143c448c26fd8a66111d0409b4bcd3eea89`, pushed, bundled as
`source-archives/athenak-feature-cgl-through-89ba4143c.bundle`, verified, and
cataloged. The bundle SHA-256 is
`16f5b14610b47e1e3ca0f695bf9e5a94fba6d450e4dfe25357f0cde54fcbf597`.

Its focused isolated fixture suite is:

```text
PYTHONDONTWRITEBYTECODE=1 python3 -B -m pytest -q -p no:cacheprovider \
  tst/test_suite/cgl/test_cgl_lf_stage_i_checkpoint.py
67 passed
```

The suite copies the historical Stage I helper into a temporary Git
repository with a local sentinel default root, so offline tests do not query
the live scheduler or resolve the shared production root. The retained
companion provides `audit-recost`, `verify-staged-recost`, `promote-recost`,
`verify-promoted-recost`, `finalize-linked-pair`, `retire-preparing`, and
`adopt-legacy-canonical`.
Canonical use authenticates the tracked clean companion and tracked clean
Stage I helper from the pinned repository root, executes reconcile from an
authenticated helper descriptor, requires the root flock and an empty live
user queue, authenticates the complete artifact authorization context, and
uses a same-directory `link/fsync/unlink/fsync` publication transition.
Separate durable recovery journals and per-artifact forensic directories
preserve fail-closed recovery after pre-link, linked-pair, canonical-only, or
post-audit interruption states. The adoption action is separate from normal
publication: it requires a fresh free lock, empty queue, exact canonical-only
namespace, authenticated existing bytes, clean reconciliation, an
adoption-specific durable journal, a mode-`0444` forensic copy, and an audit
record with `record_type = "legacy-canonical-adoption"`,
`original_publication_transition_observed = false`, and unknown original
publication method and publisher. Its three durable phases are resumable.

## Fail-Closed Resume Checklist

Resume sequentially. Shared-root mutations and queue submissions must not
overlap with another agent. Steps 1 through 5 are closed at this checkpoint;
do not rerun adoption.

1. Closed: commit, push, independently review, archive, catalog, and use the
   retained `scripts/frontier/cgl_lf_stage_i_checkpoint.py` companion.
2. Closed: adopt F-101 exactly once with `adopt-legacy-canonical` at a fresh
   free-lock and empty-queue boundary, without claiming observation of the
   original publication transition.
3. Closed: invoke `verify-promoted-recost` separately and inspect the exact
   audit, canonical artifact, retained mode-`0444` forensic copy, empty recost
   transaction directory, empty queue, and free strict root lock recorded
   above.
4. Closed: commit and push the documentation-only post-adoption checkpoint, create
   and catalog a complete-history source bundle, verify the full checksum
   ledger, and independently audit the archive.
5. Closed: before any `s19` lifecycle mutation, harden the retained Stage I helper's
   canonical flock opening to the companion's `O_NOFOLLOW`, regular-file,
   owner, mode-`0644`, one-link profile and post-lock pathname identity;
   harden all Slurm control-plane routing; and reject production fixtures and
   `--skip-slurm-test` in a separate reviewed, committed, pushed, archived,
   and cataloged transition. F-101 binds the prior helper SHA-256, so this
   helper transition remains after F-101 adoption rather than before it.
6. Require a fresh empty queue and free root-lock boundary before each
   production mutation. Run authenticated hardened reconciliation and
   require `19/19/19`, no active reservation, no Stage I transaction, and
   `issues = []` before preparing `s19`.
7. Prepare and independently audit the bounded `s19` readiness packet. Only
   then run `check-submit` and `submit`, preserving explicit acknowledgement
   of the reviewed stale shared-root campaign
   `beta25-accel05-gamma10001-purecgl-256`. Never mutate that stale campaign.

## Overall Phase I Status

The corrected forcing-policy epoch is qualified. `R02` is accepted from
fresh `t = 0` through exact `t = 7.25`. Retained companion implementation,
review, archival, first production use, and F-101 legacy adoption are
complete. The separate Stage I helper hardening checkpoint is also complete.
The immediate blocker is the fresh queue/free-lock reconciliation and bounded
readiness-packet gate before `s19`. After that closes, finish `R02` through
exact `t = 10`, then execute `R03` through `R16` sequentially and `R17` last
under the existing protocol. Stage II and manuscript-result claims remain
out of scope until Phase I production and analysis gates are complete.
