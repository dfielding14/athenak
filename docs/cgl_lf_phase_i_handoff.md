# CGL-LF Phase I Production Handoff

Checkpoint captured: `2026-06-02T21:55:26Z`

## Read This First

This is the durable resume point for the corrected-E03 MKS24 Stage I
production campaign. Use
`docs/cgl_lf_weak_guide_manuscript_plan.md` as the primary scientific guiding
document, `docs/cgl_lf_mks24_reproduction_implementation_plan.md` as the
detailed implementation and evidence record, and
`docs/cgl_lf_mks24_stage_i_protocol_review.md` as the stop-the-line production
protocol. This file records the newer operational boundary that had not yet
been promoted into those longer historical records.

Do not submit a new production job immediately. The next action is to finish
the fail-closed F-101 recost publication for accepted `R02/s18`, verify it,
refresh the three longer documents, and archive that committed documentation
state. Only then may the bounded `R02/s19_rankio_t7p25_t7p5` readiness packet
be prepared and audited.

## Repository Boundary

Work in:

```text
/autofs/nccs-svm1_home2/dfielding/athenak-df
```

The branch at this checkpoint is `feature/cgl-landau-fluid`. Its last
committed controller/documentation state before this handoff was:

```text
61c2ee98f745e08e8e6042361ebee365da381240
Record E03 R02 t7 quarter recost
```

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

At checkpoint capture, the canonical counterpart was absent:

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

The staged artifact has been independently audited against the retained
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
recommendation inside staged evidence, not permission to prepare or submit
`s19` before F-101 is canonically promoted, verified, documented, committed,
and archived.

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

## Fail-Closed Resume Checklist

Resume sequentially. Shared-root mutations and queue submissions must not
overlap with another agent.

1. Confirm the repository branch and read this handoff plus the three guiding
   documents.
2. Confirm that the user queue is empty with `squeue -u "$USER"`.
3. Confirm the recost namespace contains only the staged F-101 file above,
   with the expected checksum, mode, and one link. Confirm the canonical path
   and hidden forensic variants are absent.
4. Confirm
   `accounting/mks24_stage_i_E03_forcing_policy_transactions` is empty and
   the canonical lock is available.
5. Run authenticated hardened reconciliation:

   ```bash
   scripts/frontier/cgl_lf_stage_i.py \
     --root /lustre/orion/ast207/proj-shared/dfielding/CGL reconcile
   ```

   Require `19/19/19`, no active reservation, no transaction, and
   `issues = []`.
6. Reconstruct a one-shot byte-preserving promoter from reviewed requirements
   or recover equivalent reviewed source. Retain or otherwise durably record
   its exact bytes before execution. It must authenticate the committed live
   helper SHA-256
   `1c633ebb58294938a0a0609742ae8f8d1f88cd15242ffe649578796edeb39375`
   as mode `0644`, authenticate the generator and staged artifact above,
   hold the Stage I lock, fail closed on namespace ambiguity, preserve a
   forensic copy after an ambiguous link attempt, and perform a fresh queue
   check immediately before `os.link()`.
7. Obtain independent static review of that exact promoter and repeat live
   queue, transaction, namespace, and lock checks.
8. Promote F-101 exactly once, then verify canonical-only namespace, canonical
   SHA-256
   `df208829aec85c88fcc2caa757a21ad3ff3372aeb82628cbfe7bd41851faed9a`,
   mode `0644`, one link, empty queue, empty transaction directory, and clean
   hardened reconciliation.
9. Update the three guiding documents with the accepted `s18`/F-101 record,
   commit documentation only, create and catalog a complete-history source
   bundle, and independently audit the archive.
10. Prepare and independently audit the bounded `s19` readiness packet. Only
    then run `check-submit` and `submit`, preserving explicit acknowledgement
    of the reviewed stale shared-root campaign
    `beta25-accel05-gamma10001-purecgl-256`. Never mutate that stale campaign.

## Overall Phase I Status

The corrected forcing-policy epoch is qualified. `R02` is accepted from
fresh `t = 0` through exact `t = 7.25`. F-101 publication is the immediate
blocker. After it closes, finish `R02` through exact `t = 10`, then execute
`R03` through `R16` sequentially and `R17` last under the existing protocol.
Stage II and manuscript-result claims remain out of scope until Phase I
production and analysis gates are complete.
