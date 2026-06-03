# CGL-LF Phase I Production Handoff

Checkpoint refreshed: `2026-06-03T01:32:34Z`

## Read This First

This is the durable resume point for the corrected-E03 MKS24 Stage I
production campaign. Use
`docs/cgl_lf_weak_guide_manuscript_plan.md` as the primary scientific guiding
document, `docs/cgl_lf_mks24_reproduction_implementation_plan.md` as the
detailed implementation and evidence record, and
`docs/cgl_lf_mks24_stage_i_protocol_review.md` as the stop-the-line production
protocol.

Do not submit a new production job immediately. Corrected-E03
`R02/s19_rankio_t7p25_t7p5` is accepted and its F-103 recost artifact has
completed the retained companion's normal observed-publication lifecycle.
Before any `s20` mutation, commit and push these documentation bytes, create
and catalog the resulting complete-history source bundle, verify the checksum
ledger, run authenticated hardened reconciliation, require a fresh empty user
queue and strict free root lock, and independently audit the bounded
`R02/s20_rankio_t7p5_t7p75` readiness packet.

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
28113fd06c8da05c5458bb4f54c7b3477751d3a6
Record Stage I helper hardening checkpoint
```

It is pushed and archived as:

```text
source-archives/athenak-feature-cgl-through-28113fd06.bundle
SHA-256 e1910db1dcb008a6b39888536e7b6fabf527fac525fb3550b466a7801ffafaf0
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

## Accepted R02 s19

Corrected-E03 `R02/s19_rankio_t7p25_t7p5` completed as Slurm job `4757300`
and is formally inspected, independently validated, recorded `accepted`, and
reconciled. The exact retained scheduler row is:

```text
4757300|cgl_mks24_E03_forcing_policy_R02_s19_rankio_t7p25_t7p5|COMPLETED|0:0|1|2450|2026-06-02T20:18:10|2026-06-02T20:59:12
```

The allocation reached exact `t = 7.5`, used `2450` elapsed seconds
(`0.680556` displayed node-hours), and brings the accepted corrected-E03
prefix to `62870 / 3600 = 17.46388888888889` exact node-hours
(`17.463890` displayed). It retained one complete eight-sibling snapshot
group and one terminal eight-sibling restart group, zero strict LF failure
counters, finite synchronized histories, no restart-marker bypass, terminal
`lf_hwproj = 279167641195`, and sampled-history forcing-work residual
`8.520423303056633e-12`.

Retained evidence:

| Artifact | SHA-256 |
| --- | --- |
| `runs/mks24-stage-i/E03-forcing-policy/R02/s19_rankio_t7p25_t7p5/manifest/segment_inspection.json` | `8b9cfb6f3d0ec215234d9481e031299da1286117424d756f654270930881f744` |
| `runs/mks24-stage-i/E03-forcing-policy/R02/s19_rankio_t7p25_t7p5/manifest/prepared_run.json` after record | `30445e1c64dc35c64db38f119ff6838c904ee4d71f78ec75b219c4d476bb7b30` |
| `accounting/4757300.stage_i.independent_validation.json` | `6ebd5054e8d8082d75343b9cd176d7faf8c3b773774d8e145302bf3bfd859dcf` |
| `accounting/4757300.stage_i.sacct.txt` | `9a45c253369f87f8f7a31ce1189ac2df5b312b293c058a6f5bd9aa2c9e953f48` |

Recorded reconciliation closed with `20/20/20`
ledger rows/manifests/reservations, no active reservation, no Stage I
transaction, and `issues = []`.

## Current F-103 Boundary

The retained and independently reviewed F-103 generator is:

```text
accounting/utilities/generate_cgl_r02_t7p5_recost.py
SHA-256 131ab294888fd70bad7140a781ebc43ea6cb05bc4eddcb15dbee91c2a445ba28
mode 0755, one link
```

An initial reviewed candidate was rejected before installation because one
diagnostic string still named F-101. The corrected exact bytes above passed
independent review before installation. A first generator invocation with
relative CLI paths failed closed before any retained or temporary output.
The subsequent absolute-path invocation regenerated the retained independent
validation and staged recost artifact atomically.

The retained companion separately ran `verify-staged-recost`,
`promote-recost`, and `verify-promoted-recost`. Normal publication used its
same-directory `link/fsync/unlink/fsync` transition. The canonical artifact
is:

```text
accounting/mks24_stage_i_E03_forcing_policy_R02_t7p5_recost_evidence.json
SHA-256 ab7fa80746fa4d2c65f515b3ace71d5668586e27fde91bbdeb578fc3a4509d98
mode 0644, one link
```

The durable observed-publication audit is:

```text
accounting/mks24_stage_i_E03_forcing_policy_R02_t7p5_recost_evidence.json.publication_audit.json
SHA-256 28e8526eb3c225cde028476dea95657b170264ea7f77156b8a839e2fe03bbf6b
mode 0644, one link
record_type observed-publication
publication same-directory-link-fsync-unlink-fsync
published_utc 2026-06-03T01:28:52+00:00
transaction_id 2026-06-03T012845+0000-1871451eb479441faa3d8a29f17d3e3c
```

Its independent forensic copy is:

```text
accounting/mks24_stage_i_E03_forcing_policy_recost_forensics/
  mks24_stage_i_E03_forcing_policy_R02_t7p5_recost_evidence.json/
  2026-06-03T012845+0000-1871451eb479441faa3d8a29f17d3e3c.mks24_stage_i_E03_forcing_policy_R02_t7p5_recost_evidence.json.forensic
SHA-256 ab7fa80746fa4d2c65f515b3ace71d5668586e27fde91bbdeb578fc3a4509d98
mode 0444, one link
```

Post-publication inspection found the staged twin absent, both Stage I
transaction directories empty, the live user queue empty, and the root lock
free with strict profile
`mode0644|links1|uid18664|regularTrue|same_inodeTrue`.

The artifact's sole recommendation is:

```text
R02/s20_rankio_t7p5_t7p75
one node
Slurm walltime 01:05:00
Athena timeout 00:55:00
one-segment threshold 2700 seconds
600-second Athena and Slurm guards
```

The recost uses observed `2450`, local acceleration estimate `2560`,
selected estimate `2560`, projection `829.2566672222222`, and margin
`70.74333277777782` node-hours. The threshold is a reviewed authorization
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
6. Required before `s20`: commit and push this documentation checkpoint,
   create and catalog its complete-history source bundle, and validate the
   full source-archive checksum ledger.
7. Required before each production mutation: require an empty user queue,
   a free strict root lock, authenticated hardened reconciliation with
   `20/20/20`, no active reservation, no transaction, and `issues = []`.
8. Prepare and independently audit only the bounded
   `R02/s20_rankio_t7p5_t7p75` packet. Only then run `check-submit` and
   `submit`, preserving explicit acknowledgement of reviewed stale shared-root
   campaign `beta25-accel05-gamma10001-purecgl-256`. Never mutate that stale
   campaign.

## Overall Phase I Status

The corrected forcing-policy epoch is qualified. `R02` is accepted from
fresh `t = 0` through exact `t = 7.5`. Retained companion implementation,
review, archival, F-101 legacy adoption, F-102 helper hardening, and the
first hardened normal observed publication at F-103 are complete. The
immediate blocker is the documentation commit, complete-history archive, and
fresh queue/free-lock reconciliation and readiness-packet gate before `s20`.
After that closes, finish `R02` through exact `t = 10`, execute `R03` through
`R16` sequentially, and execute `R17` last under the same fail-closed
protocol. Stage II and manuscript-result claims remain out of scope until
Phase I production and analysis gates are complete.
