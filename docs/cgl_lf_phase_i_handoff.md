# CGL-LF Phase I Production Handoff

Checkpoint refreshed: `2026-06-03T19:49:18Z`

## Read This First

This is the durable resume point for the corrected-E03 MKS24 Stage I
production campaign. Use
`docs/cgl_lf_weak_guide_manuscript_plan.md` as the primary scientific guiding
document, `docs/cgl_lf_mks24_reproduction_implementation_plan.md` as the
detailed implementation and evidence record, and
`docs/cgl_lf_mks24_stage_i_protocol_review.md` as the stop-the-line production
protocol.

Do not submit a new production job immediately. Corrected-E03
`R02/s27_rankio_t9p25_t9p75` is accepted and its F-111 recost artifact has
completed the retained companion's normal observed-publication lifecycle.
Before any `s28` mutation, commit and push these documentation bytes, create
and catalog the resulting complete-history source bundle, verify the checksum
ledger, run authenticated hardened reconciliation, require no queued `cgl_*`
workflow job and a strict free root lock, and independently audit the bounded
`R02/s28_rankio_t9p75_t10` readiness packet. Actual AthenaK submission
retains the Stage I helper's separate all-user queue preflight.

## Repository Boundary

Work in:

```text
/autofs/nccs-svm1_home2/dfielding/athenak-df
```

The branch is `feature/cgl-landau-fluid`. The retained lifecycle companion is:

```text
scripts/frontier/cgl_lf_stage_i_checkpoint.py
SHA-256 9b54c840ea1d5df6a54d9eec2cecdd239067bd13f60bfeb38f9495fc189b8ac6
mode 0755, one link
commit 7866f8cbb50bbb9ab799ff8a3cc2140f5457945b
```

The hardened Stage I helper is:

```text
scripts/frontier/cgl_lf_stage_i.py
SHA-256 54ec671bb45aa27735a174d40b4b2e6009070716346ea09699bbe62421bbfada
mode 0644, one link
commit f675bd677cb582a46bbda8b5f55c335e6f259dd3
```

The last documentation checkpoint is:

```text
b9fcf7a5eaa21025270e0a79d404dfda658e778f
Update final F-110 launch pointer
```

It is pushed and archived as:

```text
source-archives/athenak-feature-cgl-through-b9fcf7a5e.bundle
SHA-256 e85c66511cc399cbe893b9abcbf7c64160907b2167a35745d2135d2a9cec761e
```

The shared production root is:

```text
/lustre/orion/ast207/proj-shared/dfielding/CGL
```

Do not modify or stage unrelated working-tree state:

```text
 M kokkos
?? scripts/frontier/rollback_codex_node_local_runtime.sh
```

## Accepted R02 s27

Corrected-E03 `R02/s27_rankio_t9p25_t9p75` completed as Slurm job `4759856`
and is formally inspected, independently validated, recorded `accepted`, and
reconciled. The exact retained scheduler row is:

```text
4759856|cgl_mks24_E03_forcing_policy_R02_s27_rankio_t9p25_t9p75|COMPLETED|0:0|1|4889|2026-06-03T13:20:13|2026-06-03T14:50:49
```

The allocation reached exact `t = 9.75`, used `4889` elapsed seconds
(`1.358056` displayed node-hours), and brings the accepted corrected-E03
prefix to `84051 / 3600 = 23.3475` exact node-hours
(`23.347502` displayed). It retained two complete eight-sibling snapshot
groups and one terminal eight-sibling restart group, zero strict LF failure
counters, finite synchronized histories, no restart-marker bypass, terminal
`lf_hwproj = 308629878188`, and sampled-history forcing-work residual
`4.667144504856836e-13`.

Retained evidence:

| Artifact | SHA-256 |
| --- | --- |
| `runs/mks24-stage-i/E03-forcing-policy/R02/s27_rankio_t9p25_t9p75/manifest/segment_inspection.json` | `16c840e2ccc62d9a9648ca00fc7bc0240157fc99ef675dbf7f290a2ed827f004` |
| `runs/mks24-stage-i/E03-forcing-policy/R02/s27_rankio_t9p25_t9p75/manifest/prepared_run.json` after record | `4ec34805bc682adff0e9ea7a5e55e6fbe45b699b823d04c0534326ade1810c2e` |
| `accounting/4759856.stage_i.independent_validation.json` | `32cf08cbc44fde55a3c4adfffb95262fb1e10e285204cec9ff752dc7fa10b513` |
| `accounting/4759856.stage_i.sacct.txt` | `ff88a519cb7e8f357a9c252a04803381d76452451d155cf7648fc5d0c33af723` |
| `accounting/utilities/validate_cgl_stage_i_segment.py` | `340288e71183bf80b6a8122d200daa136d08362a465c6e7a46dd20448ca00896` |

Recorded reconciliation closed with `28/28/28`
ledger rows/manifests/reservations, no active reservation, no Stage I
transaction, and `issues = []`.

## Current F-111 Boundary

The retained and independently reviewed F-111 generator is:

```text
accounting/utilities/generate_cgl_r02_t9p75_recost.py
SHA-256 bf0cbee8aca8431d3669ecaf8a3682cab04e3282d59001653c365ff90e996199
mode 0755, one link
```

Pre-install review verified the fully bound `s27` endpoint, strict
existing-file-only root-lock opening, fixed scheduler paths with stripped
Slurm environment, fixed source-git routing with stripped Git environment,
derived threshold prose, distinct predecessor and successor walltime profiles,
the generalized multi-snapshot validator, and the strict scoped `%i|%j|%T`
queue parser. The parser blocks queued
`cgl_*` workflow jobs and malformed rows while permitting unrelated account
jobs. The corrected exact bytes above passed independent review before
installation, regenerated the retained independent validation byte-for-byte,
and emitted the staged recost artifact atomically.

The retained companion separately ran `verify-staged-recost`, staged
`audit-recost`, `promote-recost`, `verify-promoted-recost`, and promoted
`audit-recost`. Normal publication used its
same-directory `link/fsync/unlink/fsync` transition. The canonical artifact
is:

```text
accounting/mks24_stage_i_E03_forcing_policy_R02_t9p75_recost_evidence.json
SHA-256 354dde3f93a6ebe73532e2b7f851a8068675b7798964c707146c37428cb9b195
mode 0644, one link
```

The durable observed-publication audit is:

```text
accounting/mks24_stage_i_E03_forcing_policy_R02_t9p75_recost_evidence.json.publication_audit.json
SHA-256 bd20e3d7e1e9d1df8cb77ea1b443601a9586f55cdb413f3b6f205cc2283152e6
mode 0644, one link
record_type observed-publication
publication same-directory-link-fsync-unlink-fsync
published_utc 2026-06-03T19:35:41+00:00
transaction_id 2026-06-03T193540+0000-d18a65ff58164015bb42b30c0de86642
```

Its independent forensic copy is:

```text
accounting/mks24_stage_i_E03_forcing_policy_recost_forensics/
  mks24_stage_i_E03_forcing_policy_R02_t9p75_recost_evidence.json/
  2026-06-03T193540+0000-d18a65ff58164015bb42b30c0de86642.mks24_stage_i_E03_forcing_policy_R02_t9p75_recost_evidence.json.forensic
SHA-256 354dde3f93a6ebe73532e2b7f851a8068675b7798964c707146c37428cb9b195
mode 0444, one link
```

Post-publication inspection found the staged twin absent, both Stage I
transaction directories empty, no queued `cgl_*` workflow job, and the root
lock free with strict profile
`mode0644|links1|uid18664|regularTrue|same_inodeTrue`.

The artifact's sole recommendation is:

```text
R02/s28_rankio_t9p75_t10
one node
Slurm walltime 01:10:00
Athena timeout 01:00:00
one-segment threshold 3000 seconds
600-second Athena and Slurm guards
```

The recost uses observed half-unit elapsed time `4889`, rounded twenty-percent
final-quarter estimate `3000`, selected final-quarter estimate `3000`,
projection `897.7011116666667`, and margin `2.2988883333332524` node-hours.
The threshold remains a reviewed authorization
and post-run recost threshold, not a controller-enforced kill deadline, and
does not ratchet automatically.

## Aggressive Successor Policy

Finish R02 through exact `t = 10` and analyze the complete standard-layout
case once. For R03--R16 use the largest evidence-backed standard-layout
segments that fit the hard two-hour Slurm limit while retaining both
`600`-second guards. The current reviewed aggressive profile is half-unit
segmentation, not the R02 quarter-unit calibration cadence. Execute R17 last
and derive its largest bounded segments separately from retained eight-node
high-resolution timing evidence.

Supporting correspondence is retained as
`Majeski_AthenaK_turbulence_driver.pdf`, SHA-256
`8f30cd53f630bdc2e4af181fcdb64c092e02279bf7d952e46d1dc0ce764f674b`.
Stephen Majeski reports expected active/passive turbulence-driver behavior:
matching magnetic-field and anisotropy distributions, instability-threshold
volume fractions, Kolmogorov spectra, active perpendicular pressure balance,
its passive absence, and selective suppression of `bb:grad(u)`. This supports
the scientific expectations but does not replace retained campaign validation
or the F-111 timing evidence that authorizes the aggressive schedule.

## Historical F-101 And F-102

F-101 remains truthful historical evidence. Its externally canonical bytes
were adopted exactly once under the retained companion without claiming
observation of the original publication transition:

```text
accounting/mks24_stage_i_E03_forcing_policy_R02_t7p25_recost_evidence.json
SHA-256 df208829aec85c88fcc2caa757a21ad3ff3372aeb82628cbfe7bd41851faed9a
publication audit SHA-256 9bbfc6ee1007b0e17a4408444acb1bcdb2f59cbde9183569d506fe2606c1837d
record_type legacy-canonical-adoption
adopted_utc 2026-06-02T23:35:28+00:00
```

Do not rerun F-101 adoption. F-102 separately hardened canonical lock opening
and Slurm control-plane routing before any `s19` mutation. F-103 is the first
normal observed publication under those hardened helper bytes.

## Fail-Closed Resume Checklist

Resume sequentially. Shared-root mutations and queue submissions must not
overlap with another agent.

1. Closed: retain, review, commit, push, archive, catalog, and use the
   routine lifecycle companion.
2. Closed: adopt and separately verify historical F-101 exactly once.
3. Closed: harden the Stage I helper lock profile and scheduler routing as
   F-102 before any `s19` mutation.
4. Closed: prepare, submit, inspect, independently validate, record, and
   reconcile `R02/s19_rankio_t7p25_t7p5`.
5. Closed: retain and independently review the F-103 generator; generate,
   verify staged, independently audit, normally promote, separately verify
   promoted, and inspect canonical-only F-103 evidence.
6. Closed: commit, push, archive, catalog, and checksum-validate the post-F-103
   documentation checkpoint and complete-history source bundle.
7. Closed: prepare, submit, inspect, independently validate, record, and
   reconcile `R02/s20_rankio_t7p5_t7p75`.
8. Closed: retain and independently review the F-104 generator; generate,
   verify staged, independently audit, normally promote, separately verify
   promoted, and inspect canonical-only F-104 evidence.
9. Closed: commit, push, archive, catalog, and checksum-validate the post-F-104
   documentation checkpoint and complete-history source bundle.
10. Closed: prepare, submit, inspect, independently validate, record, and
   reconcile `R02/s21_rankio_t7p75_t8`.
11. Closed: retain and independently review the F-105 generator; generate,
   verify staged, independently audit, normally promote, separately verify
   promoted, and inspect canonical-only F-105 evidence.
12. Closed: commit, push, archive, catalog, and checksum-validate the post-F-105
   documentation checkpoint and complete-history source bundle.
13. Closed: prepare, submit, inspect, independently validate, record, and
   reconcile `R02/s22_rankio_t8_t8p25`.
14. Closed: retain and independently review the F-106 generator; generate,
   verify staged, independently audit, normally promote, separately verify
   promoted, and inspect canonical-only F-106 evidence.
15. Closed: commit and push the post-F-106 documentation checkpoint,
   create and catalog its complete-history source bundle, and validate the
   full source-archive checksum ledger.
16. Closed for `s23`: require an empty user queue,
   a free strict root lock, authenticated hardened reconciliation with
   `23/23/23`, no active reservation, no transaction, and `issues = []`.
17. Closed: prepare and independently audit only the bounded
   `R02/s23_rankio_t8p25_t8p5` packet. Only then run `check-submit` and
   `submit`, preserving explicit acknowledgement of reviewed stale shared-root
   campaign `beta25-accel05-gamma10001-purecgl-256`. Never mutate that stale
   campaign.
18. Closed: inspect, independently validate, record, reconcile, retain and
   independently review the F-107 generator; generate, verify staged,
   independently audit, normally promote, separately verify promoted, and
   inspect canonical-only F-107 evidence.
19. Closed: commit and push the post-F-107 documentation checkpoint,
   create and catalog its complete-history source bundle, and validate the
   full source-archive checksum ledger.
20. Closed for `s24`: require an empty user queue,
   a free strict root lock, authenticated hardened reconciliation with
   `24/24/24`, no active reservation, no transaction, and `issues = []`.
21. Closed: prepare and independently audit only the bounded
   `R02/s24_rankio_t8p5_t8p75` packet. Only then run `check-submit` and
   `submit`, preserving explicit acknowledgement of reviewed stale shared-root
   campaign `beta25-accel05-gamma10001-purecgl-256`. Never mutate that stale
   campaign.
22. Closed: commit and push the post-F-108 documentation checkpoint,
   create and catalog its complete-history source bundle, and validate the
   full source-archive checksum ledger.
23. Closed for `s25`: require an empty user queue,
   a free strict root lock, authenticated hardened reconciliation with
   `25/25/25`, no active reservation, no transaction, and `issues = []`.
24. Closed: prepare and independently audit only the bounded
   `R02/s25_rankio_t8p75_t9` packet. Only then run `check-submit` and `submit`,
   preserving explicit acknowledgement of reviewed stale shared-root campaign
   `beta25-accel05-gamma10001-purecgl-256`. Never mutate that stale campaign.
25. Closed: commit and push the post-F-109 documentation checkpoint,
   create and catalog its complete-history source bundle, and validate the
   full source-archive checksum ledger.
26. Closed for `s26`: require an empty user queue,
   a free strict root lock, authenticated hardened reconciliation with
   `26/26/26`, no active reservation, no transaction, and `issues = []`.
27. Closed: prepare and independently audit only the bounded
   `R02/s26_rankio_t9_t9p25` packet. Only then run `check-submit` and `submit`,
   preserving explicit acknowledgement of reviewed stale shared-root campaign
   `beta25-accel05-gamma10001-purecgl-256`. Never mutate that stale campaign.
28. Closed: inspect, independently validate, record, reconcile, retain and
   independently review the F-110 generator; generate, verify staged,
   independently audit, normally promote, separately verify promoted, and
   inspect canonical-only F-110 evidence.
29. Closed: commit and push the post-F-110 documentation checkpoint,
   create and catalog its complete-history source bundle, and validate the
   full source-archive checksum ledger.
30. Closed: prepare, submit, inspect, independently validate, record,
   reconcile, retain and independently review the F-111 generator; generate,
   verify staged, independently audit, normally promote, separately verify
   promoted, and inspect canonical-only F-111 evidence for the aggressive bounded
   `R02/s27_rankio_t9p25_t9p75` packet. Require no queued `cgl_*` workflow job,
   a free strict root lock, authenticated hardened reconciliation with
   `27/27/27`, no active reservation, no transaction, and `issues = []`.
   Before actual AthenaK submission, retain the Stage I helper's all-user
   queue preflight and explicit stale beta-25 acknowledgement.
31. Required before `s28`: commit and push this documentation checkpoint,
   create and catalog its complete-history source bundle, and validate the
   full source-archive checksum ledger.
32. Prepare and independently audit only the bounded final-quarter
   `R02/s28_rankio_t9p75_t10` packet. Require no queued `cgl_*` workflow job,
   a free strict root lock, authenticated hardened reconciliation with
   `28/28/28`, no active reservation, no transaction, and `issues = []`.
   Before actual AthenaK submission, retain the Stage I helper's all-user
   queue preflight and explicit stale beta-25 acknowledgement.

## Overall Phase I Status

The corrected forcing-policy epoch is qualified. `R02` is accepted from
fresh `t = 0` through exact `t = 9.75`. Retained companion implementation,
review, archival, F-101 legacy adoption, F-102 helper hardening, and hardened
normal observed publications at F-103 through F-111 are complete. The
immediate blocker is the documentation commit, complete-history archive, and
fresh scoped-queue/free-lock reconciliation and readiness-packet gate before
`s28`. After that closes, finish and analyze `R02` through exact `t = 10`,
execute `R03` through `R16` with aggressive largest bounded standard-layout
segments, and execute `R17` last under separately recosted eight-node bounds.
Stage II and manuscript-result claims remain out of scope until Phase I
production and analysis gates are complete.
