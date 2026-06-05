# PIC Sun and Bai release qualification handoff

Last updated: 2026-06-04T17:27:44Z

## Purpose

This document is the durable restart point for the AthenaK MHD-PIC release
qualification effort. Read it before changing code, launching jobs, or removing
artifacts. The release goal is a rock-solid MHD-PIC implementation that passes
the available test matrix and reliably reproduces the non-relativistic shock
acceleration experiment in Section 5.4 and Figure 8 of Sun and Bai,
*"The Magnetohydrodynamic-Particle-In-Cell Module in Athena++:
Implementation and Code Tests."*

The later driven turbulent-box science campaign is intentionally out of scope
for this qualification effort. External terminal review is also deferred until
the local qualification work is complete.

## Executive status

Status: **partially complete; do not launch the qualifying Section 5.4 campaign
yet.**

Completed:

1. The dirty development tree was curated into reviewable commits and pushed
   to `origin/PIC`.
2. A clean executable candidate was frozen under the Orion PIC artifact root.
   Its source commit, executable checksum, deck checksums, analyzer checksums,
   policy chain, attestations, manifests, ledger entries, and run artifacts are
   retained.
3. The controller boundary was hardened to require the Orion artifact root,
   source binding, checksums, serialized execution, registered science scope,
   reconciliation, and a 10,000 node-hour cap.
4. Four small pressure-engineering calibration slices completed successfully,
   were reconciled, and have immutable case descriptors.
5. Qualifying-campaign planning, attempt publication, and numerical
   qualification scaffolding were implemented and committed as explicitly
   launch-prohibited code.
6. Before the current hardening pass, the focused Q011 source tranche passed
   129 tests and the explicit full publication test suite passed all 58 test
   modules.
7. The fifth-pass pressure publication repair passed its worker-routed
   48-test suite in job `4759380`. No aggregate artifact is public yet.
8. The sixth-pass control-plane, planner, and campaign-lifetime repair passed
   its worker-routed 441-test suite in job `4759417`. A final clean committed
   snapshot validation remains required.
9. The seventh-pass pressure terminal-seal and caller-level reconciliation
   repair passed its exact amended worker-routed 45-test suite in job
   `4759724`. Independent read-only rereview found no remaining
   pressure-publication code blocker in the modeled path-scoped writer
   boundary.
10. The eighth-pass authenticated storage-preflight, historical-policy
    migration and candidate-revalidation repair passed its worker-routed
    499-test suite in job `4759732`.
11. The ninth-pass exact historical-anchor and retained candidate-lifecycle
    repair passed its final worker-routed 507-test suite in job `4759760`.
    Independent read-only adversarial rereview found no remaining retained
    clean-candidate lifecycle blocker.
12. The integrated committed snapshot passed worker job `4759862`: 200 focused
    tests and all 1244 explicit publication tests passed, with two intentional
    skips.
13. The first live strict-storage preflight attempt exposed a canonical-path
    defect before controller installation or policy mutation. The tenth-pass
    repair passed a bounded read-only replay of the actual live predecessor.
    Preliminary worker job `4760218` passed its then-current 449-test suite
    before independent rereview identified the required canonical-policy and
    frozen-ledger-root split.
14. Worker job `4760511` passed the amended 611-test split-root suite. Worker
    job `4760559` passed the superseding exact-final-controller-byte 611-test
    suite after a formatting-only closure.
15. The reviewed tenth-pass operational repair was committed and pushed as
    `ee18be0f04dab3c7110c3372ccd904288d13ec7c`. Its first closed worker
    validation submission became Frontier job `4761480`.
16. The eleventh-pass scheduler-token repair was committed and pushed as
    `7e3f30a38d201ae438c434f9ce2ee2049f1f6455`. Its reviewed replacement
    worker validation submission became Frontier job `4761489`.
17. The twelfth-pass fixed Cray-Python worker path was committed and pushed as
    `749c95bb7492d67bd3765eadd5287f54663c4052`. Worker job `4761549` passed
    `200` archived focused tests and all `1246` explicit publication tests with
    two intentional skips.
18. Authenticated storage probe
    `1554766c-21e2-48b1-8cfe-b1e7e4e75aa2`, paired controller install
    `821d185856722bd0178acb9427f78ac82671a4b6670779ec8400fbac54c6d721`,
    historical-slice retirement, worker build-freeze job `4761634`, immutable
    candidate freeze `98a372c9-2ea0-47e6-ad34-e66343e7eea1`, independent
    candidate revalidation, and candidate-only policy promotion all passed.
19. The first live acceptance-root provision failed closed before `sbatch`:
    Orion inherited the parent setgid bit and created the exact empty sibling
    directory with mode `02700`, not the required `0700`.
20. The thirteenth-pass helper was committed and pushed as
    `10bb501df0fa66d70f95f8983494a2156dd9ebd6`. Worker job `4764398` passed
    `224` archived focused tests and all `1270` explicit publication tests with
    two intentional skips. Two exact-latest-patch rereviews passed.
21. The reviewed one-time descriptor-relative recovery normalized the retained
    acceptance-root inode from `02700` to `0700` without replacement and
    published durable receipt SHA-256
    `784bc2ba4d3ee3cc198d1b34bdab393f1edd92149e4304a99c0b165423a0c9fc`.
22. Aggregate worker job `4764415` failed closed before exposing any public
    artifact. All twenty retained `mhd_w_bcc` snapshots have one exact
    eight-field inventory, but AthenaK writes the fields in producer order while
    the analyzer assumed one different tuple order.
23. The fourteenth-pass exact unordered inventory repair was committed and
    pushed as `0ac4bb1bf8d0fe10867db2b7d776e692acea9b0b`. Worker job `4764639`
    passed `240` archived focused tests and all `1273` explicit publication
    tests with two intentional skips. An exact-latest-patch rereview passed.
24. Replacement aggregate worker job `4764674` failed closed before exposing
    any public artifact. AthenaK writes interior `dt=15` products on the first
    committed step after each nominal cadence, while the analyzer still
    required exact nominal times.
25. Read-only forensics verified all `100` retained snapshot products and
    froze the exact observed particle-header time, six-significant-digit
    mesh-header time, and common cycle tuple for every immutable case slot.
26. Keep a separate qualification-launch blocker open: before any qualifying
    materialization, audit the attempt-manifest binary/particle time equality
    and the numerical analyzer's nominal-time tolerance against AthenaK's
    committed-step output scheduling. This is not part of the immutable pilot
    publication recovery and must not broaden its accepted contract.
27. Independent software and scientific reviews found and repaired a stale
    review-packet binary-reader call contract, a compatibility-sidecar
    hash-check/use reopen window, a predecessor chronology fork, and missing
    coherent unauthorized-cycle coverage.

In progress:

1. Pressure-pilot aggregate publication remains prohibited. No public bundle,
   receipt, aggregate analysis artifact, or review packet exists.
2. Preserve the exact empty recovered `publication_acceptance/` inode and its
   durable receipt. Do not delete or recreate the directory and do not
   republish the receipt.
3. Commit, push, validate, and independently rereview the fifteenth-pass exact
   retained snapshot-metadata repair before submitting one reviewed
   replacement aggregate publication worker.

Remaining:

1. Commit and push the fifteenth-pass aggregate snapshot-metadata repair,
   pass one clean committed worker validation, and independently rereview the
   exact latest patch.
2. Publish and verify the pressure-pilot
   aggregate receipt and review packet on workers.
3. Perform a human pressure-choice review from the immutable four-slice
   calibration bundle. Do not invent or silently auto-select this science
   decision.
4. Materialize and validate a qualifying campaign plan only after the repaired
   launch boundary passes adversarial review.
5. Run the prerequisite slices, cheaper paper-test suites, the full Section 5.4
   reproduction campaign, and final publication analysis in staged order.
6. Finish optional extensions, final hardening, and the deferred external
   review before treating the implementation as production-ready.

## Mandatory operating constraints

1. Do not use Kronos. Put retained artifacts only under:

   ```text
   /lustre/orion/ast207/proj-shared/dfielding/PIC
   ```

2. Do not launch the qualifying Section 5.4 campaign until the known
   adversarial findings are repaired, reviewed, and committed.
3. Do not submit production simulations directly with `sbatch`. Use the
   registered control-plane wrappers and preserve policy, attestation,
   reservation, manifest, and ledger bindings. Direct `/usr/bin/sbatch` is
   allowed only for documented non-simulation worker validation and
   publication wrappers.
4. Keep Frontier jobs serialized until the qualification boundary is reviewed:
   one active job at a time.
5. Do not discard dirty worktree edits or hidden publication staging trees
   until their ownership and value have been checked.
6. Do not commit generated `__pycache__` directories.
7. Bind every production-facing run to a source commit, executable checksum,
   deck checksum, analyzer checksum, and the Orion PIC artifact root.
8. Submit every simulation and every materially long publication or analysis
   task to Slurm worker nodes. Login-node work is limited to bounded
   orchestration, inspection, and short focused checks.

## Repository state

Primary worktree:

```text
/autofs/nccs-svm1_home2/dfielding/athenak-pic
```

The retained clean-candidate builder still uses the reviewed lexical spelling
`/ccs/home/dfielding/athenak-pic`. Do not rewrite that build provenance
spelling during the Project Home migration.

Repair base before the staged tenth-pass commit:

```text
branch:     PIC
HEAD:       7cb4a20eb626113d112a0512726d2754fb81569e
origin/PIC: 7cb4a20eb626113d112a0512726d2754fb81569e
subject:    Harden Q011 validation checkout boundary
```

After restart, require `git rev-parse HEAD origin/PIC` to report the same
current pushed tip. Treat that tip, not the historical repair base above, as
the source checkpoint for the next operational step.

The integrated repair commit curates the previously dirty fourth through ninth
repair tranches. They span pressure publication, planner retention, raw-attempt
materialization, retained publication rollback, numerical restart provenance,
immutable-tree helpers, strict storage migration, retained candidate
revalidation, readiness chronology, and their focused tests. Inspect any new
diff before editing or committing. The expected repair ownership is:

| Repair area | Primary files |
| --- | --- |
| Numerical evidence binding | `analyze_q011_section54_numerical_qualification.py`, its test |
| Campaign admission graph | `analyze_q011_section54_campaign.py`, its test |
| Attempt publication race | `publish_q011_section54_campaign_attempt.py`, its test |
| Planner publication race and helper closure | `q011_section54_qualifying_campaign_execution.py`, its test |
| Generic pre-submit planner provenance | `frontier_control_plane/`, its tests |
| Pressure aggregate publication | `publish_q011_section54_pressure_pilot_bundle.py`, its test |
| Pressure review packet | `render_q011_section54_pressure_pilot_review_packet.py`, its test |
| Retained-tree descriptor pinning | `immutable_orion_tree.py`, its test |
| Canonical Project Home policy and frozen lexical ledger roots | `frontier_control_plane/`, `q011_section54_pressure_pilot_execution.py`, registry and controller tests |
| Hermetic Q011 workers and post-run source authorization | `frontier_q011_*job.sh`, `analyze_q011_section54_pressure_pilot.py`, readiness successor, prepared inventory |
| Publication-acceptance retained-inode recovery | `provision_q011_pressure_publication_acceptance_root.py`, its focused test, controller runbook, readiness transition |

Recent pushed history:

```text
749c95bb7 Pin Q011 validation worker Python path
7e3f30a38 Normalize Frontier parsable job tokens
ee18be0f0 Harden Q011 canonical storage operational boundary
7cb4a20eb Harden Q011 validation checkout boundary
7d157d52b Split Q011 archived and trusted checkout validation
e93f7bfd5 Run Q011 validation from read-only Git clone
d87a8e434 Refresh Q011 clean-candidate fixtures
f37155643 Clarify Q011 packet verification ordering
442242516 Fix Q011 operational migration runbook
ec47056ee Harden Q011 publication and lifecycle boundary
7324031b1 Run Q011 repair validation from a clean worker snapshot
f3128598b Record Q011 third adversarial repair transition
04959b44b Bind Q011 aggregate analysis parser compatibility successor
01141bd8a Accept empty optional Athena binary header values
7a813e0bd Expand worker Q011 repair validation coverage
eefbd60d6 Add worker-routed Q011 pressure review packet renderer
dd0f66053 Record worker-routed Q011 second hardening pass
08a235f35 Add durable PIC release qualification handoff
ecc90185c Add launch-prohibited Q011 qualifying campaign scaffolding
29a180d4b Record Q011 four-slice pilot policy promotion
0d0caf468 Close local PIC shock planning reviews
9c0524a38 Record repaired Q011 clean candidate freeze
672818de7 Record strict Q011 launch-prohibited controller boundary
3ba8c9577 Harden Q011 retry reservation and promotion boundary
597e0e371 Repair Q011 HIP particle copy and register retry
5a3ed6a3c Record hardened PIC baseline promotion chronology
```

The governing plans are:

```text
tst/publication/PIC_PRODUCTION_READINESS_PLAN.md
tst/publication/PIC_SUN_BAI_RELEASE_PHASED_WORK_PLAN.md
```

## Frozen candidate and controller bindings

Artifact root:

```text
PIC_ROOT=/lustre/orion/ast207/proj-shared/dfielding/PIC
```

Current fresh candidate-only checkpoint:

```text
FREEZE_ROOT=$PIC_ROOT/clean_candidates/98a372c9-2ea0-47e6-ad34-e66343e7eea1
CLEAN_CANDIDATE_MANIFEST=$FREEZE_ROOT/clean_candidate_manifest.json
EXECUTABLE=$FREEZE_ROOT/athena
ENVIRONMENT_PROFILE=$PIC_ROOT/control_plane/821d185856722bd0178acb9427f78ac82671a4b6670779ec8400fbac54c6d721/frontier_pic_environment.sh
```

| Item | Value |
| --- | --- |
| Frozen source commit | `749c95bb7492d67bd3765eadd5287f54663c4052` |
| Controller version | `821d185856722bd0178acb9427f78ac82671a4b6670779ec8400fbac54c6d721` |
| Storage probe ID | `1554766c-21e2-48b1-8cfe-b1e7e4e75aa2` |
| Clean-candidate manifest SHA-256 | `dee6be45657e99ec477eec513c45be4ba6ac43b4fd5deb5655750f51c668e42f` |
| Executable SHA-256 | `53e72bf57818a451e219cac950a891b48f23dc333fe8d2648a7907e51950c3ce` |
| Candidate-only active policy SHA-256 | `23a73b868146f63d4b363713f988d55e9dadffa15b26f2b2f1d07da66331f5c3` |
| Candidate-only active promotion SHA-256 | `4824ea825e7b9e42becdca4b9a8b72c0454bd1a2e02d5e365b1878ed94c53243` |

Historical consumed-slice clean candidate:

```text
FREEZE_ROOT=$PIC_ROOT/clean_candidates/0d1cabe6-a282-408b-953a-c2db8f0e1310
CLEAN_CANDIDATE_MANIFEST=$FREEZE_ROOT/clean_candidate_manifest.json
EXECUTABLE=$FREEZE_ROOT/athena
ENVIRONMENT_PROFILE=$PIC_ROOT/control_plane/6dc456e5858b2cbdbd7ef80cbc324b748e2cce5075023b49f7115eab4e84dfdd/frontier_pic_environment.sh
```

Bindings:

| Item | Value |
| --- | --- |
| Frozen source commit | `672818de70dcff51543a2fa76c88aa5af1a59a65` |
| Controller version | `6dc456e5858b2cbdbd7ef80cbc324b748e2cce5075023b49f7115eab4e84dfdd` |
| Clean-candidate manifest SHA-256 | `ef527ed467995bd60fda07b5a3b09b56ea871595ace12fd64a948e246720dbe3` |
| Executable SHA-256 | `0c249d339f5131a72b45ce244d84d26421fea08cab2493f0db1e3ff6cdf2c420` |
| Installed runtime profile SHA-256 | `7ff6cff3e263e8b4552bcbac1d5b55a85264a2a7af5565bbc9856584e35402c1` |
| Frozen build profile JSON SHA-256 | `a2b1ef2266a2b1097e4453d05f8838fcae67ecb819c759773210c04e22d0ad05` |
| Active policy SHA-256 | `647573e109852da0e343cd1c80033dbbe4ad4d2d672588cdea4e852340bad3d8` |
| Active promotion SHA-256 | `073da1d4fe7f2eb054da9f3ec2bf2860f2643c3f44ca88d6183f4592eb0681bf` |

Historical consumed-slice pressure-pilot policy:

```text
$PIC_ROOT/policy/reviewed_q011_pressure_pilot_successor_v2_0d1cabe6.json
```

The cap is 10,000 node-hours. The ledger tail after the fourth engineering
calibration reconciliation reported cumulative consumption of
`1.3019444444444446` node-hours.

The older controller, candidate and policy bindings are retained historical
chronology. Do not relaunch them. The current fresh candidate-only checkpoint
above supersedes them for the next staged qualification steps.

## Completed pressure-engineering calibration slices

These are calibration runs, not a Section 5.4 reproduction claim and not
authorization for a qualifying launch.

| Slice | Slurm job | Submission ID | Consumed node-hours | Descriptor SHA-256 |
| --- | --- | --- | --- | --- |
| `ps_p0_1p00` | `4754520` | `a9db403f-a174-45ae-956d-f96185a94db1` | `0.06444444444444444` | `da322586da4831c7747e3a2f5774adf6b999835dc8636b9d971411e1dc9471b2` |
| `ps_p0_0p05` | `4754533` | `3eb84dfa-0be7-4a62-b77b-64131c08dc6f` | `0.06277777777777778` | `978cfc63efc9503b9dd0a043793d759e930f390340b843a39d621c2b859ff85c` |
| `ps_p0_0p10` | `4754558` | `07bd4d96-ffc8-4446-98ef-562d4d83da0c` | `0.04083333333333333` | `dba35e7b157c1646444792cfda68963dbfd2062b4887961bfb380e51a113f5d4` |
| `ps_p0_0p20` | `4754581` | `8b4065a3-b02a-4312-8238-20d72d958bbb` | `0.03916666666666667` | `ef842e1440cbe1fd51e3ee79af5c838704e782708e2a827ec0bfcc89ca1ae312` |

Each immutable descriptor is retained as:

```text
$PIC_ROOT/runs/<campaign>/<submission-id>/analysis/analysis.json
```

The fourth reconciliation event is:

```text
event SHA-256: 50ed1c9010cd2c837ea5a83eba9d96db7b3474edd3db17cc26d830ae7129bd31
ledger file:   $PIC_ROOT/ledger/node_hours.jsonl
```

## Pressure-pilot aggregate publication

An initial login-node publisher attempt proved materially long-running and did
not expose a final receipt. Five Frontier Slurm worker publication attempts
then ended without exposing a public bundle:

| Slurm job | Terminal state | Elapsed | Reason |
| --- | --- | --- | --- |
| `4756951` | `FAILED` | `00:09:39` | Strict offline binary parser rejected runtime-added empty optional `particles/pic_deltaf_f0` |
| `4757026` | `CANCELLED` | `00:02:16` | Cancelled when the checkout mutated during publication; superseded by worker-local `git archive HEAD` snapshots |
| `4757047` | `CANCELLED` | `00:04:36` | Cancelled after fresh review found late-failure rollback, descriptor-anchor, and source-authorization gaps |
| `4764415` | `FAILED` | `00:05:37` | Strict offline `mhd_w_bcc` analyzer assumed one tuple order for an exact name-addressed AthenaK variable inventory |
| `4764674` | `FAILED` | `00:05:34` | Strict offline analyzer assumed exact nominal times for AthenaK interior `dt`-scheduled products |

The integrated fourteenth-pass tree and failed job `4764674` are historical.
Aggregate publication remains prohibited until the fifteenth-pass exact
retained snapshot-metadata repair is committed, pushed, clean-worker validated
and independently rereviewed from the exact latest patch. Do not invoke the
publisher directly on a login node. Use only the exact hermetic, commit-bound,
`--export=NIL` publication sequence in
[`frontier_control_plane/README.md`](frontier_control_plane/README.md).

Target outputs:

```text
$PIC_ROOT/publication/q011_section54_pressure_pilot_bundle
$PIC_ROOT/publication/q011_section54_pressure_pilot_bundle_receipt.json
$PIC_ROOT/publication/q011_section54_pressure_pilot_analysis.json
```

The two hidden pre-repair staging trees were audited and removed on
`2026-06-02T22:59:11Z`. One was empty. The other contained only three partial
`ps_p0_1p00` magnetic-field binary copies and no receipt. Neither was accepted
evidence:

```text
$PIC_ROOT/publication/.q011_section54_pressure_pilot_bundle.staging-a38200f5-06a0-47c8-8f24-046e0fdb556a
$PIC_ROOT/publication/.q011_section54_pressure_pilot_bundle.staging-6b50871f-7318-41f5-9c36-c5cb4ce1dfa8
```

Do not assume that either hidden tree still exists. The publication root was
empty immediately after the reviewed cleanup.

Check the publisher and visible outputs:

```bash
export PIC_ROOT=/lustre/orion/ast207/proj-shared/dfielding/PIC
pgrep -af publish_q011_section54_pressure_pilot_bundle || true
find "$PIC_ROOT/publication" -maxdepth 2 -printf '%M %s %p\n' | sort
```

After the corrected publisher has been reviewed, committed, pushed, and
worker-validated, inspect the hidden staging trees and rerun from the reviewed
Slurm wrapper. The wrapper performs final receipt verification on the worker
node. Do not use a login-node direct invocation or a raw wrapper submission
without the reviewed commit argument.

## Known launch blockers under repair

An independent adversarial review found five release-blocking defects in the
launch-prohibited qualifying scaffolding. None invalidates the completed
pressure-engineering calibration slices. All must be repaired before
qualifying launch:

1. `analyze_q011_section54_numerical_qualification.py` accepted caller-supplied
   evidence dictionaries without reopening retained artifacts and
   independently recomputing the admitted evidence.
2. `analyze_q011_section54_campaign.py` could admit a self-authored local
   pressure-selection and campaign-plan graph instead of proving provenance
   from the immutable published pilot receipt and the authorized materializer.
3. `publish_q011_section54_campaign_attempt.py` validated a staging tree before
   rename but could expose an invalid retained destination if a same-owner
   mutation occurred between validation and publication.
4. `q011_section54_qualifying_campaign_execution.py` created the deterministic
   final output root before population and lacked a strict exact-member
   whitelist, leaving a publication injection window.
5. The qualifying planner helper closure omitted qualification, particle,
   spatial, artifact, and pressure-publication verification dependencies from
   the retained source set.

Treat the code as launch-prohibited until the repaired implementations have
focused tests, the full publication suite, diff checks, bytecode checks, and a
fresh independent adversarial review.

### Second adversarial pass

The first repair tranche closed the original five defects structurally, but a
fresh independent pass found additional blockers:

1. Paired AMR-versus-fine residual records must be recomputed from retained raw
   attempt trees, not trusted from caller-authored derived JSON.
2. Restart parity observations must be reconstructed from retained restart and
   post-checkpoint artifacts, not trusted from caller-authored dictionaries.
3. Retained numerical recompute bundles must be canonically contained below
   the authorized Orion PIC root.
4. Campaign-plan admission must validate the exact production planner graph
   and every referenced child, not accept a minimal tree with dangling
   bindings.
5. Campaign-plan admission must cross-link the environment profile and restart
   preregistration to retained authoritative bytes.
6. The generic retained derived-bundle publisher must use atomic no-replace
   publication, descriptor-pinned postverification, and fail-closed rollback.
7. Add a real planner-to-attempt manifest materializer and production-path
   integration test; synthetic fixtures alone are insufficient.

The second repair tranche is active. Do not materialize or launch a qualifying
campaign until it passes another adversarial review.

### Third adversarial pass

The second repair tranche passed a focused `145`-test integration matrix, but
a fresh independent pass found additional blockers before commit:

1. Pressure aggregate publication needs rollback after every late failure,
   descriptor-anchored freeze and cleanup, and a complete reviewed source
   successor closure.
2. Pressure review-packet publication needs the equivalent rollback, anchor,
   source-snapshot, and receipt hardening.
3. Planner and attempt publishers must guarantee zero public exposure when a
   mutable staging pathname is substituted immediately before rename.
4. Attempt final admission self-checks must stay inside the rollback envelope.
5. Raw attempt materialization must bind a reconciled registered-execution
   receipt and exact deterministic planner-authorized roots transactionally
   before mutation.
6. Restart qualification must reject aliased branches and bind the planned
   checksum-cross-linked registered continuation.
7. Numerical recompute bundles need post-open root binding, exact directory
   closure, and exact analyzer-binding closure.

The durable transition record is
[`q011_section54_third_adversarial_repair_transition_2026-06-02.json`](readiness/q011_section54_third_adversarial_repair_transition_2026-06-02.json).
The third repair tranche is active. Do not publish the aggregate bundle,
materialize a qualifying campaign, or launch qualifying work until it passes
another independent adversarial review.

### Fourth adversarial pass

The first third-tranche implementations were independently rereviewed before
worker publication. The rereview rejected them and opened a fourth repair
tranche:

1. Consumed historical pressure-v2 slices must be structurally incapable of
   authorizing fresh launches, even if historical source bytes are restored.
2. Pressure review packets must verify the exact aggregate-receipt-bound
   analysis digest before rendering metrics or figures.
3. Worker archive provenance must be enforced rather than accepted from
   caller-supplied environment strings.
4. Operational registered-execution artifacts below `$PIC_ROOT/runs` and
   deterministic retained publication below `$PIC_ROOT/campaigns` must remain
   distinct namespaces with an immutable cross-link.
5. Reconciliation must emit the read-only registered-execution receipt only
   after the mirrored ledger event is durably appended, deriving receipt
   fields from the immutable pre-submit manifest and reconciliation event.
6. The qualifying planner needs a reviewed materializer for the exact
   `planner_retention` object injected into pre-submit configuration. Do not
   hand-author that binding during launch orchestration.
7. The generic pre-submit boundary must independently prove that the
   `planner_retention` overlay came from immutable qualifying-planner bytes.
   Shape-only validation of caller-authored root and `argv` fields is
   insufficient.
8. Completed-attempt materialization must rebind the opened `/runs/.../raw`
   descriptor to its public pathname before and after mutation.
9. Retained-publication rollback must not delete a substituted destination
   while leaving the originally published inode exposed under a moved public
   name.
10. Pressure publication archive proof must reject synthetic read-only tar
    files that merely carry a forged 40-hex PAX comment. Aggregate and packet
    publishers must recheck retained publication-root identity after their
    final pathname verifier and reject boolean schema-version aliases.
11. Restart qualification must prove authoritative execution receipts,
   admitted source-checkpoint lineage, and the admitted source-plan carrier.
12. Restart post-checkpoint outputs must be checksum-cross-linked to immutable
    registered raw-tree inventories rather than accepted as bundle-authored
    copies. Fixed-path receipt and mirrored-ledger consumption must retain
    pinned ancestry for the complete validation lifetime.
13. Numerical recompute closure needs the complete direct analyzer dependency
    set and a second stability scan.

The durable transition record is
[`q011_section54_fourth_adversarial_repair_transition_2026-06-02.json`](readiness/q011_section54_fourth_adversarial_repair_transition_2026-06-02.json).
The fourth repair tranche is active. Do not publish, freeze a replacement
candidate, materialize a qualifying campaign, or launch qualifying work until
the repaired boundary passes worker validation and another independent
adversarial review.

### Fifth pressure-publication pass

The fourth-pass pressure implementation passed its scoped worker tests, but an
independent rereview found that a coordinated byte-identical replacement could
leave a canonical receipt verifiable after rollback failed. It also found that
possessing genuine `git archive HEAD` bytes did not prove execution from the
worker-extracted source snapshot. The fifth pressure pass closes both gaps:

1. Aggregate and packet publishers arm a receipt-specific fail-closed guard
   before canonical receipt rename.
2. Public receipt consumers reject guarded receipts before and after retained
   closure verification.
3. Successful publication removes the guard only as its final action.
4. Failed rollback leaves the guard retained whenever canonical withdrawal
   cannot be proven, including coordinated byte-identical replacement.
5. Production source binding proves that the running Python module and reviewed
   source members come from the read-only worker-extracted `git archive HEAD`
   snapshot. Direct API execution from the trusted checkout is rejected.

Worker validation job `4759380` passed all `48` scoped pressure tests in
`143.870s`. The durable transition record is
[`q011_section54_fifth_adversarial_repair_transition_2026-06-03.json`](readiness/q011_section54_fifth_adversarial_repair_transition_2026-06-03.json).
Pressure aggregate publication remains prohibited until the integrated
clean-worker suite and a fresh independent rereview pass.

### Sixth control-plane and campaign-lifetime pass

The concurrent sixth pass closes the remaining retained-lifetime and planner
provenance gaps found during integration:

1. Campaign admission keeps the mirrored execution-ledger snapshot pinned
   through retained raw-product, restart, and telemetry consumption.
2. Ledger snapshot retention rejects transient receipt hide-and-restore
   mutation even when final bytes and inode match.
3. Hidden staged campaign publication uses the same retained ledger-snapshot
   lifetime as direct aggregate admission.
4. Qualifying planner retention reconstructs the complete reviewed graph,
   helper closure, frozen candidate, exact matrix, contracts, restart carrier,
   recompute plan, policy fragment, and pressure-selection review binding.
5. Pre-submit creation, reservation, ledger retention, and reconciliation
   cross-bind planner retention to the authorized submission clean-candidate
   manifest digest.

Worker validation job `4759417` passed all `441` control-plane tests in
`140.865s`. The durable transition record is
[`q011_section54_sixth_adversarial_repair_transition_2026-06-03.json`](readiness/q011_section54_sixth_adversarial_repair_transition_2026-06-03.json).
The prepared artifact inventory now binds `110` decks and `20` analyzers with
SHA-256 `72a88a2b8f71d6d1fa7f7ed6195767b91b475367081dca397459f644c7a8b251`.
Launch and aggregate publication remain prohibited until the integrated tree
is committed, pushed, validated from a clean worker snapshot, and rereviewed.

### Seventh pressure terminal-seal pass

The next pressure rereview found a rollback-after-commit edge: a seal helper
wrapper could call the real terminal rename, prove the result consumable and
then raise before its caller assigned `seal_committed = True`. Aggregate and
review-packet publishers now reconcile that broader caller-level exception
path by reopening the exact receipt-bound seal and requiring the public
fail-closed guard to remain absent before returning committed success.

Worker job `4759724` passed the exact amended `45`-test pressure suite in
`248.894s`. Independent read-only rereview found no remaining
pressure-publication code blocker under the modeled path-scoped ordinary
concurrent-writer threat. The durable transition record is
[`q011_section54_seventh_adversarial_repair_transition_2026-06-03.json`](readiness/q011_section54_seventh_adversarial_repair_transition_2026-06-03.json).

The live fixed sibling authority directory was initially absent. Its first
reviewed provisioning attempt failed closed before `sbatch`: Orion inherited
the parent setgid bit and created the exact empty
`$PIC_ROOT/publication_acceptance/` inode with mode `02700`, not `0700`.
The thirteenth pass preserved and recovered that inode in place to exact mode
`0700`, then published its durable receipt. Preserve both. The same-UID PIC-root
mutation caveat and Lustre power-loss durability caveat remain explicit
operational risks.

### Eighth storage migration and candidate-revalidation pass

The historical live policy predates authenticated mirrored storage evidence.
The eighth pass adds a strict migration boundary:

1. Source-only storage capture probes the fixed Orion artifact root and
   Project Home ledger-mirror root, then publishes byte-identical read-only
   evidence below both roots.
2. Normal policy promotion, unlock and reservation require that evidence.
3. One narrow promotion-only flag permits the exact historical predecessor to
   be retired into a newer-controller, empty-allowlist, pending-freeze policy.
   Runtime unlock and reservation have no legacy mode.
4. Clean-candidate reads watch retained PIC-root ancestry and reject parent
   substitution or rename-away-and-restore.
5. Worker freeze creation immediately invokes the installed read-only
   candidate verifier against the exact emitted manifest digest.

Worker job `4759732` passed the integrated `499`-test lifecycle suite in
`133.923s`. The staged controller digest is
`dab112c006506b99415fd4885efa1dc778632ae446120dd576ff860dab5e4be9`.
The durable transition record is
[`q011_section54_eighth_adversarial_repair_transition_2026-06-03.json`](readiness/q011_section54_eighth_adversarial_repair_transition_2026-06-03.json).
The ninth pass below supersedes this intermediate lifecycle checkpoint.

### Ninth exact-anchor and retained-candidate lifecycle pass

Adversarial rereview of the eighth pass found avoidable watcher-lifetime gaps
around ancestry recheck and teardown. It also identified that the generic
historical migration shape could be narrower. The ninth pass closes those
boundaries:

1. The one-use migration branch binds the exact mirrored historical live policy
   SHA-256 `647573e109852da0e343cd1c80033dbbe4ad4d2d672588cdea4e852340bad3d8`
   and promotion SHA-256
   `073da1d4fe7f2eb054da9f3ec2bf2860f2643c3f44ca88d6183f4592eb0681bf`.
2. Candidate reads drain events after ancestry recheck and close retained
   members, candidate directories and PIC-root ancestry descriptors before the
   final watcher drain.
3. Missing candidate-root construction failures close the already-created
   watcher before returning.
4. Six adversarial late-window substitutions and the missing-root fd-count
   cleanup now have focused regressions.

Worker job `4759760` passed the exact final `507`-test lifecycle suite in
`133.031s`. Independent read-only rereview replayed every identified interval
and found no remaining retained clean-candidate lifecycle blocker. The staged
controller digest is
`4bb093662299911c870bf68d57de9f22897e197d95baf8d303e34b0f5e083c8c`.
The durable transition record is
[`q011_section54_ninth_adversarial_repair_transition_2026-06-03.json`](readiness/q011_section54_ninth_adversarial_repair_transition_2026-06-03.json).

### Tenth canonical Project Home path pass

The first reviewed live preflight attempt found that strict no-symlink storage
probing and the historical `/ccs/proj` Project Home alias were incompatible by
construction. Both attempted captures failed closed before evidence
publication, controller installation, or policy promotion. The active policy
and promotion hashes remained unchanged.

The tenth pass moves fresh successor policy state to the physical Project Home
root `/autofs/nccs-svm1_proj/ast207/proj-shared/PIC`. The one-use retirement
branch still authenticates the historical `/ccs/proj` predecessor spelling
after checking that it resolves to the authorized physical root. That narrow
compatibility is required for the retained predecessor manual-accounting
mirrors and sealed predecessor attestation; it is not allowed for fresh
evidence.

Independent rereview after worker job `4760218` found a second-order
requirement: the append-only ledger mirror must retain the historical lexical
root `/ccs/proj/ast207/proj-shared/PIC`. Existing ledger receipts, genesis
anchors and sealed operator attestations bind those exact path bytes. Fresh
storage evidence, installed policy controllers, successor policy mirrors and
successor manual-accounting authorization mirrors use the canonical `/autofs`
root. Submission wrappers continue to use `/ccs/proj` only for the frozen
ledger mirror.

A bounded read-only replay passed against the actual live predecessor and all
`123` mirrored ledger records after the split-root repair. The focused
storage-preflight and Q011 successor suite passed all `35` tests. Worker job
`4760218` remains the preliminary canonical-path validation checkpoint: it
passed the then-current 449-test affected suite in `141.748s` before the
second-order ledger split was identified. Worker job `4760511` passed the final
amended 611-test affected split-root suite in `157.366s`. The staged controller
received a formatting-only closure, then superseding worker job `4760559`
passed the exact-final-byte 611-test affected split-root suite in `155.181s`.
Independent read-only code rereview found no remaining code-path blocker; its
two runbook findings were repaired before commit. The staged controller digest
is
`821d185856722bd0178acb9427f78ac82671a4b6670779ec8400fbac54c6d721`.
The durable transition record is
[`q011_section54_tenth_canonical_project_home_repair_transition_2026-06-03.json`](readiness/q011_section54_tenth_canonical_project_home_repair_transition_2026-06-03.json).

### Tenth-pass operational wrapper closure

Operational rereview of the split-root migration path found additional
fail-closed requirements before live mutation. The staged closure now:

1. submits build-freeze, pressure aggregate, and pressure review-packet workers
   with `--export=NIL`, fixed paths, sanitized environment variables, and an
   exact reviewed Git commit;
2. authenticates the installed controller before the build worker sources the
   installed environment;
3. runs committed repair validation from a read-only Git archive and rejects
   source drift before and after validation;
4. carries an append-only post-run pressure-source authorization successor for
   the modified publication wrappers without reauthorizing historical slices;
5. renders the engineering review packet only after proving the explicit
   plotting-lock candidate dependency versions; and
6. documents checkpoint-specific recovery for preflight evidence publication,
   paired controller installation, one-use retirement promotion, worker
   submission, and receipt publication; and
7. makes one clean committed repair-validation worker a hard prerequisite for
   live migration and pressure publication, using only pinned core discovery
   tools inside the sanitized worker.

The plotting lock remains an engineering-review candidate. It does not close
the later qualification plotting or external-export review gate. Independent
read-only operational rereview passed the exact latest dirty tranche with no
remaining blocker after the pinned-tool and hard-prerequisite repairs. Live
mutation remains prohibited until clean commit and push and committed-tree
worker validation pass.

The regenerated prepared-artifact inventory binds `110` decks and `20`
analyzers with SHA-256
`4cb878a431271d85ff64b330bffc4db1abaa7124a1a5669755da12ddb4677871`.

### Eleventh Frontier parsable-token pass

The first committed worker-validation submission after the tenth-pass repair
was accepted by Frontier as job `4761480`, but `sbatch --parsable` returned a
cluster-qualified token of the form `<job-id>;frontier`. The operator shell
required digits only and therefore exited after submission without printing the
normalized ID. No live policy mutation began.

Job `4761480` remained useful validation evidence. It passed all `200` archived
focused tests in `46.003s`, then failed closed before the full sweep because the
trusted checkout became dirty while the README token parser was repaired. Its
top-level state is `FAILED`, elapsed time is `00:00:55`, and exit code is `1:0`.

The eleventh pass normalizes only the exact optional `;frontier` suffix for all
four documented `sbatch --parsable` boundaries: repair validation, build
freeze, pressure aggregate publication, and review-packet publication. Each
snippet prints the raw token with `%q` shell escaping immediately after
submission, before strict normalization and rejection, so an accepted job
remains recoverable even if Frontier returns an unexpected token. Any token
other than `<job-id>` or `<job-id>;frontier` still fails closed.

Do not mutate live policy state from commit `ee18be0f0`. Freeze and push the
eleventh-pass documentation closure, then run the exact committed worker
validation again without editing the trusted checkout while it executes.
The failed-closed `4761480` record is the reviewed reason for exactly one
replacement validation submission after that push.

The durable append-only transition record is
[`q011_section54_eleventh_frontier_sbatch_token_transition_2026-06-03.json`](readiness/q011_section54_eleventh_frontier_sbatch_token_transition_2026-06-03.json).

### Twelfth repair-validation worker Python-path pass

The reviewed eleventh-pass replacement validation submission was accepted as
Frontier job `4761489` with raw token `4761489;frontier`; token normalization
worked as intended. The archived focused tranche passed all `200` tests in
`45.363s`. The full trusted-checkout sweep then ran `1246` tests in `457.364s`
and failed one standalone-oracle test with two intentional skips.

The isolated failure was operational, not scientific. The hardened worker set
`PATH=/usr/bin:/bin`, while the executable oracle shebang is
`#!/usr/bin/env python3`. On the compute image `/usr/bin/python3` is Python
`3.6.15`, which rejects `from __future__ import annotations`. The reviewed
runtime is Cray Python `3.11.7`.

The twelfth pass sets the worker's fixed trusted path to
`/opt/cray/pe/python/3.11.7/bin:/usr/bin:/bin`. Its separate hermetic Git probes
remain constrained to `/usr/bin:/bin`. A closed-environment standalone oracle
replay and its focused four-test module pass, and independent read-only rereview
found no remaining worker-path blocker.

Do not mutate live policy state from commit `7e3f30a38`. Freeze and push this
twelfth-pass closure, then run one reviewed replacement committed worker
validation without editing the trusted checkout while it executes. The
failed-closed `4761489` record is the reviewed reason for that replacement.

The durable append-only transition record is
[`q011_section54_twelfth_worker_python_path_transition_2026-06-03.json`](readiness/q011_section54_twelfth_worker_python_path_transition_2026-06-03.json).

### Thirteenth acceptance-root inherited-setgid pass

The first post-freeze pressure-publication provisioning attempt failed closed
before aggregate `sbatch`. The reviewed shell path used
`install -d -m 0700`, but Orion inherited the parent setgid bit. The retained
sibling authority is empty, owned by `dfielding:ast207`, has no ACL xattrs, and
has exact mode `02700`.

Two independent read-only audits agreed that the directory must not be deleted
or recreated. The staged source-side helper opens the fixed Orion PIC root
component-by-component with no-follow descriptors, opens the fixed sibling
relative to the retained parent, verifies identity, owner, group, xattrs and
contents, and permits only an explicit exact-empty `02700 -> 0700` recovery.
Fresh creation also normalizes inherited mode through `fchmod`, while ordinary
production resume is verification-only for the exact reviewed existing `0700`
root. Receipt publication remains beneath the reviewed `policy/` inode and
requires a fixed-digest authenticated helper, canonical read-only receipt, and
single-link closure. A separate explicit mode reconciles the narrowly bounded
post-link staging-alias interruption state.

Commit and push the helper, pass one clean committed repair-validation worker,
independently rereview the exact latest patch, then run the incident-specific
durable recovery-receipt block in the controller README. That block requires an
empty policy-side receipt-staging namespace before normalization; any retained
pre-link orphan is a reviewed stop, not an implicit retry input. The durable
append-only transition record is
[`q011_section54_thirteenth_acceptance_root_setgid_repair_transition_2026-06-03.json`](readiness/q011_section54_thirteenth_acceptance_root_setgid_repair_transition_2026-06-03.json).

### Fourteenth aggregate MHD-inventory-order pass

The committed thirteenth-pass worker and exact-latest rereviews passed, and the
one-time recovery retained the reviewed sibling inode while normalizing it to
`0700`. Aggregate worker job `4764415` then failed closed before exposing any
public artifact. A read-only forensic replay of all twenty retained
`mhd_w_bcc` snapshots found one exact eight-field inventory with no missing or
extra variables. AthenaK emits those variables in producer order:

```text
dens velx vely velz eint bcc1 bcc2 bcc3
```

The analyzer preregistration lists the same required names with `eint` earlier.
The binary parser already maps payload arrays by header name and rejects
duplicates. The staged compatibility repair therefore requires exact unordered
inventory equality, still rejecting missing or extra variables, and composes
every required field by name. It does not broaden the accepted physics product,
change estimators or thresholds, or authorize launch. The append-only transition
record is
[`q011_section54_fourteenth_aggregate_mhd_inventory_order_repair_transition_2026-06-04.json`](readiness/q011_section54_fourteenth_aggregate_mhd_inventory_order_repair_transition_2026-06-04.json).

### Fifteenth aggregate snapshot-metadata pass

The committed fourteenth-pass worker and exact-latest rereview passed.
Replacement aggregate worker job `4764674` then failed closed before exposing
any public artifact. AthenaK increments time before testing whether `dt`-based
outputs are due, so each interior product is written on the first committed
step after its nominal cadence. Mesh binary headers serialize that time with
default six-significant-digit stream precision; particle VTK headers retain
`max_digits10`.

Read-only forensics verified all `100` retained snapshot products. Within every
immutable case slot, the four mesh products share one binary-header time and all
five products share one cycle. The particle-header lateness above nominal is at
most `0.054575889532856081`, and every mesh-header time is exactly the
six-significant-digit rendering of its particle-header time. The staged repair
does not add a general tolerance. It freezes the exact observed tuple for every
retained case slot, checks those tuples against a strict `< 0.1` engineering
compatibility cap, and requires exact tuple matches during publication. It does
not change the scientific contract or authorize launch.

## Validation baseline

The committed launch-prohibited checkpoint passed these validations before the
current repair pass:

```bash
cd /ccs/home/dfielding/athenak-pic
PYTHONDONTWRITEBYTECODE=1 \
PYTHONPATH="$PWD:$PWD/tst/publication/frontier_control_plane" \
python3 -B -m unittest \
  tst.publication.test_analyze_q011_section54_campaign \
  tst.publication.test_q011_section54_model \
  tst.publication.test_q011_section54_particles \
  tst.publication.test_q011_section54_spatial \
  tst.publication.test_q011_section54_restart \
  tst.publication.test_q011_section54_artifacts \
  tst.publication.test_q011_section54_qualifying_campaign_preregistration \
  tst.publication.test_q011_section54_pressure_selection \
  tst.publication.test_q011_section54_qualifying_campaign_execution \
  tst.publication.test_publish_q011_section54_campaign_attempt \
  tst.publication.test_analyze_q011_section54_numerical_qualification
```

Result: `129` focused tests passed.

The full publication suite also passed when run as an explicit module list:

```bash
cd /ccs/home/dfielding/athenak-pic
mapfile -t modules < <(
  rg --files tst/publication -g 'test_*.py' |
    sed -e 's#/#.#g' -e 's#\.py$##' |
    sort
)
PYTHONDONTWRITEBYTECODE=1 \
PYTHONPATH="$PWD:$PWD/tst/publication/frontier_control_plane" \
python3 -B -m unittest "${modules[@]}"
```

Result: all `58` explicit test modules passed. A plain
`unittest discover -s tst/publication -t .` invocation is not the correct
driver because the start directory is not importable in that form.

Committed clean-snapshot worker job `4759862` superseded that older baseline:
it passed 200 focused tests and all 1244 explicit publication tests, with two
intentional skips.

After integrating the repairs, rerun focused tests, the explicit full suite,
`git diff --check`, and `python3 -B -m py_compile` for every changed Python
module before committing.

## New-agent restart checklist

Run these first:

```bash
cd /autofs/nccs-svm1_home2/dfielding/athenak-pic
export PIC_ROOT=/lustre/orion/ast207/proj-shared/dfielding/PIC

git status --short --branch
git rev-parse HEAD origin/PIC
git log --oneline --decorate -12

squeue -u "$USER" -o '%i|%j|%T|%M|%L|%N|%k'
tail -n 5 "$PIC_ROOT/ledger/node_hours.jsonl"
find "$PIC_ROOT" -path '*/markers/*' -type f -print | sort

pgrep -af publish_q011_section54_pressure_pilot_bundle || true
find "$PIC_ROOT/publication" -maxdepth 2 -printf '%M %s %p\n' | sort

find tst/publication -type d -name __pycache__ -print | sort
git diff --stat
git diff -- tst/publication
```

Then:

1. Read this file and both governing plans.
2. Verify the exact empty recovered `0700` acceptance-root inode and its durable
   receipt. Preserve both.
3. Inspect, commit, and push the fifteenth-pass exact retained
   snapshot-metadata repair.
4. Run focused and full validation suites from the clean committed worker
   snapshot.
5. Request a fresh independent adversarial review.
6. Verify `$PIC_ROOT/publication_acceptance/` at exact mode `0700` and its
   existing durable receipt. Finish the pressure-pilot worker aggregate publication and
   review packet, then verify both receipts.
7. Ask the human collaborator to review the immutable four-slice pressure
   calibration evidence and select the qualifying pressure option.
8. Only then materialize the qualifying campaign plan and proceed through the
   staged release plan.

## Phased route to the overall goal

The remaining work should stay staged:

| Phase | Purpose | Launch condition |
| --- | --- | --- |
| Repair and aggregate publication | Close known software blockers and publish the four-slice calibration evidence | No qualifying launch |
| Human pressure selection | Choose the Section 5.4 pressure parameter from immutable evidence | Human-reviewed selection receipt |
| Qualifying campaign planning | Materialize retained plan, helper closure, manifests, and receipts | Repaired boundary plus adversarial review |
| Prerequisite slices | Exercise restart, spatial, particle, artifact, and analysis paths cheaply | Each slice reconciled and reviewed |
| Paper-test reproduction suite | Reproduce the cheaper Sun and Bai validation tests | Quantitative acceptance criteria met |
| Section 5.4 reproduction | Run non-relativistic shock acceleration reproduction through retained analysis | Figure-8-style diagnostics reviewed |
| Release hardening | Optional extensions, regression closure, documentation, and clean freeze | Full local matrix passes |
| External review | Independent terminal review | Deferred until local work is complete |

Do not skip directly to the expensive Section 5.4 campaign. The qualification
boundary exists to make later shock and turbulence science trustworthy.
