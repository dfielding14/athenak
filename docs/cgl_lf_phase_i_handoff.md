# CGL-LF Phase I Production Handoff

Checkpoint refreshed: `2026-06-03T09:28:03Z`

## Read This First

This is the durable resume point for the corrected-E03 MKS24 Stage I
production campaign. Use
`docs/cgl_lf_weak_guide_manuscript_plan.md` as the primary scientific guiding
document, `docs/cgl_lf_mks24_reproduction_implementation_plan.md` as the
detailed implementation and evidence record, and
`docs/cgl_lf_mks24_stage_i_protocol_review.md` as the stop-the-line production
protocol.

Do not submit a new production job immediately. Corrected-E03
`R02/s23_rankio_t8p25_t8p5` is accepted and its F-107 recost artifact has
completed the retained companion's normal observed-publication lifecycle.
Before any `s24` mutation, commit and push these documentation bytes, create
and catalog the resulting complete-history source bundle, verify the checksum
ledger, run authenticated hardened reconciliation, require a fresh empty user
queue and strict free root lock, and independently audit the bounded
`R02/s24_rankio_t8p5_t8p75` readiness packet.

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
849b052dd0211bbab0148ef2eeec0b84823c9927
Record accepted R02 s22 F-106 checkpoint
```

It is pushed and archived as:

```text
source-archives/athenak-feature-cgl-through-849b052dd.bundle
SHA-256 adf897b8f60149bdcf89bc0a42f69f224ec6486155341411df1ccd8bd2d7b8ed
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

## Accepted R02 s23

Corrected-E03 `R02/s23_rankio_t8p25_t8p5` completed as Slurm job `4758475`
and is formally inspected, independently validated, recorded `accepted`, and
reconciled. The exact retained scheduler row is:

```text
4758475|cgl_mks24_E03_forcing_policy_R02_s23_rankio_t8p25_t8p5|COMPLETED|0:0|1|2285|2026-06-03T04:11:20|2026-06-03T04:49:40
```

The allocation reached exact `t = 8.5`, used `2285` elapsed seconds
(`0.634722` displayed node-hours), and brings the accepted corrected-E03
prefix to `72239 / 3600 = 20.066388888888888` exact node-hours
(`20.066390` displayed). It retained one complete eight-sibling snapshot
group and one terminal eight-sibling restart group, zero strict LF failure
counters, finite synchronized histories, no restart-marker bypass, terminal
`lf_hwproj = 295614667845`, and sampled-history forcing-work residual
`2.3508932051548193e-12`.

Retained evidence:

| Artifact | SHA-256 |
| --- | --- |
| `runs/mks24-stage-i/E03-forcing-policy/R02/s23_rankio_t8p25_t8p5/manifest/segment_inspection.json` | `318786f1be49ef894606f579709942bdb38f0de37c4ab23b357569b0deb36b2b` |
| `runs/mks24-stage-i/E03-forcing-policy/R02/s23_rankio_t8p25_t8p5/manifest/prepared_run.json` after record | `aa01c14794007564e708125c4be3ce1f683c3e33a92084b0f65e67944d454149` |
| `accounting/4758475.stage_i.independent_validation.json` | `08cd3c894106f295eda79571faca85dc23c1aaf621ba1a5c1eb09c253bb152fe` |
| `accounting/4758475.stage_i.sacct.txt` | `86bc48723ad2365729c415623ab4fc7baf47a69ccc83cf10410f84b01440e4fe` |

Recorded reconciliation closed with `24/24/24`
ledger rows/manifests/reservations, no active reservation, no Stage I
transaction, and `issues = []`.

## Current F-107 Boundary

The retained and independently reviewed F-107 generator is:

```text
accounting/utilities/generate_cgl_r02_t8p5_recost.py
SHA-256 94548cf12706bc4abe9449bcb0d917d09ab58fa241cc3f3ec48f730ed701c479
mode 0755, one link
```

Pre-install review verified the fully bound `s23` endpoint, strict
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
accounting/mks24_stage_i_E03_forcing_policy_R02_t8p5_recost_evidence.json
SHA-256 1e4db459d7c0b5211108dee07fe1efa643ac7a655adadb907deb57f36822d699
mode 0644, one link
```

The durable observed-publication audit is:

```text
accounting/mks24_stage_i_E03_forcing_policy_R02_t8p5_recost_evidence.json.publication_audit.json
SHA-256 7dfca5ec22cc1bdd265a75e9aec275c81848c193562f5d7661bae1a79751e7f4
mode 0644, one link
record_type observed-publication
publication same-directory-link-fsync-unlink-fsync
published_utc 2026-06-03T09:28:03+00:00
transaction_id 2026-06-03T092803+0000-c25223679853420ebb957da333440f9b
```

Its independent forensic copy is:

```text
accounting/mks24_stage_i_E03_forcing_policy_recost_forensics/
  mks24_stage_i_E03_forcing_policy_R02_t8p5_recost_evidence.json/
  2026-06-03T092803+0000-c25223679853420ebb957da333440f9b.mks24_stage_i_E03_forcing_policy_R02_t8p5_recost_evidence.json.forensic
SHA-256 1e4db459d7c0b5211108dee07fe1efa643ac7a655adadb907deb57f36822d699
mode 0444, one link
```

Post-publication inspection found the staged twin absent, both Stage I
transaction directories empty, the live user queue empty, and the root lock
free with strict profile
`mode0644|links1|uid18664|regularTrue|same_inodeTrue`.

The artifact's sole recommendation is:

```text
R02/s24_rankio_t8p5_t8p75
one node
Slurm walltime 01:05:00
Athena timeout 00:55:00
one-segment threshold 2700 seconds
600-second Athena and Slurm guards
```

The recost uses observed `2285`, local acceleration estimate `2285`,
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
19. Required before `s24`: commit and push this documentation checkpoint,
   create and catalog its complete-history source bundle, and validate the
   full source-archive checksum ledger.
20. Required before each production mutation: require an empty user queue,
   a free strict root lock, authenticated hardened reconciliation with
   `24/24/24`, no active reservation, no transaction, and `issues = []`.
21. Prepare and independently audit only the bounded
   `R02/s24_rankio_t8p5_t8p75` packet. Only then run `check-submit` and
   `submit`, preserving explicit acknowledgement of reviewed stale shared-root
   campaign `beta25-accel05-gamma10001-purecgl-256`. Never mutate that stale
   campaign.

## Overall Phase I Status

The corrected forcing-policy epoch is qualified. `R02` is accepted from
fresh `t = 0` through exact `t = 8.5`. Retained companion implementation,
review, archival, F-101 legacy adoption, F-102 helper hardening, and hardened
normal observed publications at F-103 through F-107 are complete. The
immediate blocker is the documentation commit, complete-history archive, and
fresh queue/free-lock reconciliation and readiness-packet gate before `s24`.
After that closes, finish `R02` through exact `t = 10`, execute `R03` through
`R16` sequentially, and execute `R17` last under the same fail-closed
protocol. Stage II and manuscript-result claims remain out of scope until
Phase I production and analysis gates are complete.
