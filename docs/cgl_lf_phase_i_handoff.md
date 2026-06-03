# CGL-LF Phase I Production Handoff

Checkpoint refreshed: `2026-06-03T03:35:12Z`

## Read This First

This is the durable resume point for the corrected-E03 MKS24 Stage I
production campaign. Use
`docs/cgl_lf_weak_guide_manuscript_plan.md` as the primary scientific guiding
document, `docs/cgl_lf_mks24_reproduction_implementation_plan.md` as the
detailed implementation and evidence record, and
`docs/cgl_lf_mks24_stage_i_protocol_review.md` as the stop-the-line production
protocol.

Do not submit a new production job immediately. Corrected-E03
`R02/s20_rankio_t7p5_t7p75` is accepted and its F-104 recost artifact has
completed the retained companion's normal observed-publication lifecycle.
Before any `s21` mutation, commit and push these documentation bytes, create
and catalog the resulting complete-history source bundle, verify the checksum
ledger, run authenticated hardened reconciliation, require a fresh empty user
queue and strict free root lock, and independently audit the bounded
`R02/s21_rankio_t7p75_t8` readiness packet.

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
07e044440a8a35ff8f958802d1ecc1a2e4a126f2
Record accepted R02 s19 F-103 checkpoint
```

It is pushed and archived as:

```text
source-archives/athenak-feature-cgl-through-07e044440.bundle
SHA-256 4dc5c3820b385401e0d60bfaecde9bf91fbb16467ae4ff5c6604c35652beec15
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

## Accepted R02 s20

Corrected-E03 `R02/s20_rankio_t7p5_t7p75` completed as Slurm job `4757761`
and is formally inspected, independently validated, recorded `accepted`, and
reconciled. The exact retained scheduler row is:

```text
4757761|cgl_mks24_E03_forcing_policy_R02_s20_rankio_t7p5_t7p75|COMPLETED|0:0|1|2370|2026-06-02T22:01:48|2026-06-02T22:41:35
```

The allocation reached exact `t = 7.75`, used `2370` elapsed seconds
(`0.658333` displayed node-hours), and brings the accepted corrected-E03
prefix to `65240 / 3600 = 18.122222222222224` exact node-hours
(`18.122223` displayed). It retained one complete eight-sibling snapshot
group and one terminal eight-sibling restart group, zero strict LF failure
counters, finite synchronized histories, no restart-marker bypass, terminal
`lf_hwproj = 282606449055`, and sampled-history forcing-work residual
`3.7666708324783834e-12`.

Retained evidence:

| Artifact | SHA-256 |
| --- | --- |
| `runs/mks24-stage-i/E03-forcing-policy/R02/s20_rankio_t7p5_t7p75/manifest/segment_inspection.json` | `4094abc01fb9fa2f4e690f17c549270d519a1fbbeab63c49ff8856bb4cfba791` |
| `runs/mks24-stage-i/E03-forcing-policy/R02/s20_rankio_t7p5_t7p75/manifest/prepared_run.json` after record | `06e70f554196b0da40cb343bf0c5a7bd8e9b7ee2fafdf0e88e4019a9980a6ea0` |
| `accounting/4757761.stage_i.independent_validation.json` | `eec64e5b065414912385549c2b42406683ec6fce4787bd525c357a4226ef958f` |
| `accounting/4757761.stage_i.sacct.txt` | `0d726cffe32b00f14ccad5627bf4823eb9092da491c2c5aba878a52b25040573` |

Recorded reconciliation closed with `21/21/21`
ledger rows/manifests/reservations, no active reservation, no Stage I
transaction, and `issues = []`.

## Current F-104 Boundary

The retained and independently reviewed F-104 generator is:

```text
accounting/utilities/generate_cgl_r02_t7p75_recost.py
SHA-256 ca3eebaaf0d04dd82a3daa9ee40006a888e0ac82ca74872e8bc63afb8a54e971
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
accounting/mks24_stage_i_E03_forcing_policy_R02_t7p75_recost_evidence.json
SHA-256 f8d05a4483a902f21f49c9ebb8935dba3e21d4186c40eda050ae57d169f6bb41
mode 0644, one link
```

The durable observed-publication audit is:

```text
accounting/mks24_stage_i_E03_forcing_policy_R02_t7p75_recost_evidence.json.publication_audit.json
SHA-256 625f41d9822a31ca5c9e5905a3d411885b56a03605d0c27f3107f11cc740bb5b
mode 0644, one link
record_type observed-publication
publication same-directory-link-fsync-unlink-fsync
published_utc 2026-06-03T03:21:44+00:00
transaction_id 2026-06-03T032137+0000-573af3d0ed544bb1a5e2b7b8f8f49c13
```

Its independent forensic copy is:

```text
accounting/mks24_stage_i_E03_forcing_policy_recost_forensics/
  mks24_stage_i_E03_forcing_policy_R02_t7p75_recost_evidence.json/
  2026-06-03T032137+0000-573af3d0ed544bb1a5e2b7b8f8f49c13.mks24_stage_i_E03_forcing_policy_R02_t7p75_recost_evidence.json.forensic
SHA-256 f8d05a4483a902f21f49c9ebb8935dba3e21d4186c40eda050ae57d169f6bb41
mode 0444, one link
```

Post-publication inspection found the staged twin absent, both Stage I
transaction directories empty, the live user queue empty, and the root lock
free with strict profile
`mode0644|links1|uid18664|regularTrue|same_inodeTrue`.

The artifact's sole recommendation is:

```text
R02/s21_rankio_t7p75_t8
one node
Slurm walltime 01:05:00
Athena timeout 00:55:00
one-segment threshold 2700 seconds
600-second Athena and Slurm guards
```

The recost uses observed `2370`, local acceleration estimate `2370`,
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
9. Required before `s21`: commit and push this documentation checkpoint,
   create and catalog its complete-history source bundle, and validate the
   full source-archive checksum ledger.
10. Required before each production mutation: require an empty user queue,
   a free strict root lock, authenticated hardened reconciliation with
   `21/21/21`, no active reservation, no transaction, and `issues = []`.
11. Prepare and independently audit only the bounded
   `R02/s21_rankio_t7p75_t8` packet. Only then run `check-submit` and
   `submit`, preserving explicit acknowledgement of reviewed stale shared-root
   campaign `beta25-accel05-gamma10001-purecgl-256`. Never mutate that stale
   campaign.

## Overall Phase I Status

The corrected forcing-policy epoch is qualified. `R02` is accepted from
fresh `t = 0` through exact `t = 7.75`. Retained companion implementation,
review, archival, F-101 legacy adoption, F-102 helper hardening, and hardened
normal observed publications at F-103 and F-104 are complete. The
immediate blocker is the documentation commit, complete-history archive, and
fresh queue/free-lock reconciliation and readiness-packet gate before `s21`.
After that closes, finish `R02` through exact `t = 10`, execute `R03` through
`R16` sequentially, and execute `R17` last under the same fail-closed
protocol. Stage II and manuscript-result claims remain out of scope until
Phase I production and analysis gates are complete.
