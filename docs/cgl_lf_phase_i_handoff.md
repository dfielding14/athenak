# CGL-LF Phase I Production Handoff

Checkpoint refreshed: `2026-06-03T07:36:16Z`

## Read This First

This is the durable resume point for the corrected-E03 MKS24 Stage I
production campaign. Use
`docs/cgl_lf_weak_guide_manuscript_plan.md` as the primary scientific guiding
document, `docs/cgl_lf_mks24_reproduction_implementation_plan.md` as the
detailed implementation and evidence record, and
`docs/cgl_lf_mks24_stage_i_protocol_review.md` as the stop-the-line production
protocol.

Do not submit a new production job immediately. Corrected-E03
`R02/s22_rankio_t8_t8p25` is accepted and its F-106 recost artifact has
completed the retained companion's normal observed-publication lifecycle.
Before any `s23` mutation, commit and push these documentation bytes, create
and catalog the resulting complete-history source bundle, verify the checksum
ledger, run authenticated hardened reconciliation, require a fresh empty user
queue and strict free root lock, and independently audit the bounded
`R02/s23_rankio_t8p25_t8p5` readiness packet.

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
4bb580e228d190756ff040c5ce1f0ff6c7ba198b
Record accepted R02 s21 F-105 checkpoint
```

It is pushed and archived as:

```text
source-archives/athenak-feature-cgl-through-4bb580e22.bundle
SHA-256 8e4f6dbc5a7318b3e17db33348ca8fdd346c5ba8d161faaf20c23a0fdeb49411
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

## Accepted R02 s22

Corrected-E03 `R02/s22_rankio_t8_t8p25` completed as Slurm job `4758398`
and is formally inspected, independently validated, recorded `accepted`, and
reconciled. The exact retained scheduler row is:

```text
4758398|cgl_mks24_E03_forcing_policy_R02_s22_rankio_t8_t8p25|COMPLETED|0:0|1|2345|2026-06-03T02:25:21|2026-06-03T03:06:00
```

The allocation reached exact `t = 8.25`, used `2345` elapsed seconds
(`0.651389` displayed node-hours), and brings the accepted corrected-E03
prefix to `69954 / 3600 = 19.43166666666667` exact node-hours
(`19.431668` displayed). It retained one complete eight-sibling snapshot
group and one terminal eight-sibling restart group, zero strict LF failure
counters, finite synchronized histories, no restart-marker bypass, terminal
`lf_hwproj = 290890469421`, and sampled-history forcing-work residual
`3.021717626579328e-12`.

Retained evidence:

| Artifact | SHA-256 |
| --- | --- |
| `runs/mks24-stage-i/E03-forcing-policy/R02/s22_rankio_t8_t8p25/manifest/segment_inspection.json` | `77479a51cda10b1a399a9060c5a7153dec115389cd4f0ba5a6fec167c462be49` |
| `runs/mks24-stage-i/E03-forcing-policy/R02/s22_rankio_t8_t8p25/manifest/prepared_run.json` after record | `b02e754cf521a0b0e2205bcd13947c0c12b6a443adb268f9ea91a2eeeff5a803` |
| `accounting/4758398.stage_i.independent_validation.json` | `e8854e4b233d7b1b64fcc53e69975041e052e98477b85e74b8e0dc1d5a332d08` |
| `accounting/4758398.stage_i.sacct.txt` | `e6af51d0cbbd1aef3e83df0fa99f64da77d7aa991a59c0375b21b4fe07c44ea1` |

Recorded reconciliation closed with `23/23/23`
ledger rows/manifests/reservations, no active reservation, no Stage I
transaction, and `issues = []`.

## Current F-106 Boundary

The retained and independently reviewed F-106 generator is:

```text
accounting/utilities/generate_cgl_r02_t8p25_recost.py
SHA-256 341fda3812143288806285b6f8ec3969d43d0e0bd22f150190e5e4f44c57b3b6
mode 0755, one link
```

Pre-install reviews rejected stale endpoint overrides and then required strict
existing-file-only root-lock opening, fixed scheduler paths with stripped
Slurm environment, fixed source-git routing with stripped Git environment,
derived threshold prose, and an explicit final-review marker. The corrected
exact bytes above passed independent review before installation and generated
the retained independent validation and staged recost artifact atomically.

The retained companion separately ran `verify-staged-recost`,
`promote-recost`, and `verify-promoted-recost`. Normal publication used its
same-directory `link/fsync/unlink/fsync` transition. The canonical artifact
is:

```text
accounting/mks24_stage_i_E03_forcing_policy_R02_t8p25_recost_evidence.json
SHA-256 fe78232b5396d9d3ee6a5ef0e39f7aa445cf2c17125691a7d3e63d4a67978bc2
mode 0644, one link
```

The durable observed-publication audit is:

```text
accounting/mks24_stage_i_E03_forcing_policy_R02_t8p25_recost_evidence.json.publication_audit.json
SHA-256 a8289b06194f93492f20fdc965ef6e709d3b1179fc251de29efa9af01c524b99
mode 0644, one link
record_type observed-publication
publication same-directory-link-fsync-unlink-fsync
published_utc 2026-06-03T07:36:16+00:00
transaction_id 2026-06-03T073616+0000-4f48500c5fe540069eb08da1527c1aa3
```

Its independent forensic copy is:

```text
accounting/mks24_stage_i_E03_forcing_policy_recost_forensics/
  mks24_stage_i_E03_forcing_policy_R02_t8p25_recost_evidence.json/
  2026-06-03T073616+0000-4f48500c5fe540069eb08da1527c1aa3.mks24_stage_i_E03_forcing_policy_R02_t8p25_recost_evidence.json.forensic
SHA-256 fe78232b5396d9d3ee6a5ef0e39f7aa445cf2c17125691a7d3e63d4a67978bc2
mode 0444, one link
```

Post-publication inspection found the staged twin absent, both Stage I
transaction directories empty, the live user queue empty, and the root lock
free with strict profile
`mode0644|links1|uid18664|regularTrue|same_inodeTrue`.

The artifact's sole recommendation is:

```text
R02/s23_rankio_t8p25_t8p5
one node
Slurm walltime 01:05:00
Athena timeout 00:55:00
one-segment threshold 2700 seconds
600-second Athena and Slurm guards
```

The recost uses observed `2345`, local acceleration estimate `2345`,
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
15. Required before `s23`: commit and push this documentation checkpoint,
   create and catalog its complete-history source bundle, and validate the
   full source-archive checksum ledger.
16. Required before each production mutation: require an empty user queue,
   a free strict root lock, authenticated hardened reconciliation with
   `23/23/23`, no active reservation, no transaction, and `issues = []`.
17. Prepare and independently audit only the bounded
   `R02/s23_rankio_t8p25_t8p5` packet. Only then run `check-submit` and
   `submit`, preserving explicit acknowledgement of reviewed stale shared-root
   campaign `beta25-accel05-gamma10001-purecgl-256`. Never mutate that stale
   campaign.

## Overall Phase I Status

The corrected forcing-policy epoch is qualified. `R02` is accepted from
fresh `t = 0` through exact `t = 8.25`. Retained companion implementation,
review, archival, F-101 legacy adoption, F-102 helper hardening, and hardened
normal observed publications at F-103 through F-106 are complete. The
immediate blocker is the documentation commit, complete-history archive, and
fresh queue/free-lock reconciliation and readiness-packet gate before `s23`.
After that closes, finish `R02` through exact `t = 10`, execute `R03` through
`R16` sequentially, and execute `R17` last under the same fail-closed
protocol. Stage II and manuscript-result claims remain out of scope until
Phase I production and analysis gates are complete.
