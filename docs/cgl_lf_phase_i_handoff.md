# CGL-LF Phase I Production Handoff

Checkpoint refreshed: `2026-06-03T11:24:40Z`

## Read This First

This is the durable resume point for the corrected-E03 MKS24 Stage I
production campaign. Use
`docs/cgl_lf_weak_guide_manuscript_plan.md` as the primary scientific guiding
document, `docs/cgl_lf_mks24_reproduction_implementation_plan.md` as the
detailed implementation and evidence record, and
`docs/cgl_lf_mks24_stage_i_protocol_review.md` as the stop-the-line production
protocol.

Do not submit a new production job immediately. Corrected-E03
`R02/s24_rankio_t8p5_t8p75` is accepted and its F-108 recost artifact has
completed the retained companion's normal observed-publication lifecycle.
Before any `s25` mutation, commit and push these documentation bytes, create
and catalog the resulting complete-history source bundle, verify the checksum
ledger, run authenticated hardened reconciliation, require a fresh empty user
queue and strict free root lock, and independently audit the bounded
`R02/s25_rankio_t8p75_t9` readiness packet.

## Repository Boundary

Work in:

```text
/autofs/nccs-svm1_home2/dfielding/athenak-df
```

The branch is `feature/cgl-landau-fluid`. The retained lifecycle companion is:

```text
scripts/frontier/cgl_lf_stage_i_checkpoint.py
SHA-256 10156515c4bcbfdcf57a2fe54220c2a80f0477f1a7bddbde322b9433946615c2
mode 0755, one link
commit 89ba4143c448c26fd8a66111d0409b4bcd3eea89
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
6d1f032a20b4dc069e454cfe5c98c9b9f984ae28
Record accepted R02 s23 F-107 checkpoint
```

It is pushed and archived as:

```text
source-archives/athenak-feature-cgl-through-6d1f032a2.bundle
SHA-256 62587378e8480cbe94caa174734f6d3cf6e4c56c2755b2b96ffe0e2bba525ef4
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

## Accepted R02 s24

Corrected-E03 `R02/s24_rankio_t8p5_t8p75` completed as Slurm job `4758576`
and is formally inspected, independently validated, recorded `accepted`, and
reconciled. The exact retained scheduler row is:

```text
4758576|cgl_mks24_E03_forcing_policy_R02_s24_rankio_t8p5_t8p75|COMPLETED|0:0|1|2346|2026-06-03T06:00:56|2026-06-03T06:40:27
```

The allocation reached exact `t = 8.75`, used `2346` elapsed seconds
(`0.651667` displayed node-hours), and brings the accepted corrected-E03
prefix to `74585 / 3600 = 20.718055555555555` exact node-hours
(`20.718057` displayed). It retained one complete eight-sibling snapshot
group and one terminal eight-sibling restart group, zero strict LF failure
counters, finite synchronized histories, no restart-marker bypass, terminal
`lf_hwproj = 299627446834`, and sampled-history forcing-work residual
`1.4737883591844796e-12`.

Retained evidence:

| Artifact | SHA-256 |
| --- | --- |
| `runs/mks24-stage-i/E03-forcing-policy/R02/s24_rankio_t8p5_t8p75/manifest/segment_inspection.json` | `24a7013abca332a65196552f0ce1a6031314deb56e8ddbebb130b75be4b887de` |
| `runs/mks24-stage-i/E03-forcing-policy/R02/s24_rankio_t8p5_t8p75/manifest/prepared_run.json` after record | `24f407dcc74011479f783ccb16d08536ce87a35314d1637f93b78135adb27ed8` |
| `accounting/4758576.stage_i.independent_validation.json` | `0a271b1f53f3363c2598a76356a7ca3063a38e5d8c0ab63a64e15b4ebcfa2229` |
| `accounting/4758576.stage_i.sacct.txt` | `06b86448d2c4eb8c3765d9c3912ee1cb9d946f0d07ad606dfae6456ab1075f7f` |

Recorded reconciliation closed with `25/25/25`
ledger rows/manifests/reservations, no active reservation, no Stage I
transaction, and `issues = []`.

## Current F-108 Boundary

The retained and independently reviewed F-108 generator is:

```text
accounting/utilities/generate_cgl_r02_t8p75_recost.py
SHA-256 c488c4d6f51a2d73e00e212f0e485d920876eb2ed79e8e18988788969e344d05
mode 0755, one link
```

Pre-install review verified the fully bound `s24` endpoint, strict
existing-file-only root-lock opening, fixed scheduler paths with stripped
Slurm environment, fixed source-git routing with stripped Git environment,
derived threshold prose, and an explicit final-review marker. The corrected
exact bytes above passed independent review before installation, regenerated
the already retained independent validation byte-for-byte, and emitted the
staged recost artifact atomically.

The retained companion separately ran `verify-staged-recost`, staged
`audit-recost`, `promote-recost`, `verify-promoted-recost`, and promoted
`audit-recost`. Normal publication used its
same-directory `link/fsync/unlink/fsync` transition. The canonical artifact
is:

```text
accounting/mks24_stage_i_E03_forcing_policy_R02_t8p75_recost_evidence.json
SHA-256 dbd3f5e3b9f21532da546caca4bdc66263d69c3b60b35938972edd04b2ab3e19
mode 0644, one link
```

The durable observed-publication audit is:

```text
accounting/mks24_stage_i_E03_forcing_policy_R02_t8p75_recost_evidence.json.publication_audit.json
SHA-256 503dad9d617456a3154f77931726b8d1edf3ccdc1602802e0db553670d0ca79c
mode 0644, one link
record_type observed-publication
publication same-directory-link-fsync-unlink-fsync
published_utc 2026-06-03T11:24:40+00:00
transaction_id 2026-06-03T112440+0000-33f82f43f34f42eb9b061cabcaa3debf
```

Its independent forensic copy is:

```text
accounting/mks24_stage_i_E03_forcing_policy_recost_forensics/
  mks24_stage_i_E03_forcing_policy_R02_t8p75_recost_evidence.json/
  2026-06-03T112440+0000-33f82f43f34f42eb9b061cabcaa3debf.mks24_stage_i_E03_forcing_policy_R02_t8p75_recost_evidence.json.forensic
SHA-256 dbd3f5e3b9f21532da546caca4bdc66263d69c3b60b35938972edd04b2ab3e19
mode 0444, one link
```

Post-publication inspection found the staged twin absent, both Stage I
transaction directories empty, the live user queue empty, and the root lock
free with strict profile
`mode0644|links1|uid18664|regularTrue|same_inodeTrue`.

The artifact's sole recommendation is:

```text
R02/s25_rankio_t8p75_t9
one node
Slurm walltime 01:05:00
Athena timeout 00:55:00
one-segment threshold 2700 seconds
600-second Athena and Slurm guards
```

The recost uses observed `2346`, local acceleration estimate `2407`,
selected estimate `2560`, projection `829.2566672222222`, and margin
`70.74333277777782` node-hours. The threshold remains a reviewed authorization
and post-run recost threshold, not a controller-enforced kill deadline, and
does not ratchet automatically.

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
22. Required before `s25`: commit and push this documentation checkpoint,
   create and catalog its complete-history source bundle, and validate the
   full source-archive checksum ledger.
23. Required before each production mutation: require an empty user queue,
   a free strict root lock, authenticated hardened reconciliation with
   `25/25/25`, no active reservation, no transaction, and `issues = []`.
24. Prepare and independently audit only the bounded
   `R02/s25_rankio_t8p75_t9` packet. Only then run `check-submit` and `submit`,
   preserving explicit acknowledgement of reviewed stale shared-root campaign
   `beta25-accel05-gamma10001-purecgl-256`. Never mutate that stale campaign.

## Overall Phase I Status

The corrected forcing-policy epoch is qualified. `R02` is accepted from
fresh `t = 0` through exact `t = 8.75`. Retained companion implementation,
review, archival, F-101 legacy adoption, F-102 helper hardening, and hardened
normal observed publications at F-103 through F-108 are complete. The
immediate blocker is the documentation commit, complete-history archive, and
fresh queue/free-lock reconciliation and readiness-packet gate before `s25`.
After that closes, finish `R02` through exact `t = 10`, execute `R03` through
`R16` sequentially, and execute `R17` last under the same fail-closed
protocol. Stage II and manuscript-result claims remain out of scope until
Phase I production and analysis gates are complete.
