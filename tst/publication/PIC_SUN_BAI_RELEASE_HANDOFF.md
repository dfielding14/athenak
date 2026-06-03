# PIC Sun and Bai release qualification handoff

Last updated: 2026-06-03T17:06:00Z

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

In progress:

1. Pressure-pilot aggregate publication remains prohibited. No public bundle,
   receipt, aggregate analysis artifact, or review packet exists.
2. Preserve all integrated dirty worktree edits while the final clean-worker
   validation, commit and push finish.
3. The live controller and policy intentionally remain at the historical
   consumed-slice generation until the reviewed one-time strict-storage
   migration is performed from the committed source.

Remaining:

1. Commit and push the integrated source and pass the clean committed worker
   validation.
2. Capture authenticated mirrored storage evidence, install the paired
   controller, retire the consumed historical slices, freeze and independently
   revalidate a fresh worker-built candidate, and promote the candidate-only
   policy.
3. Provision the fixed sibling publication-acceptance authority through the
   reviewed one-time transition, then publish and verify the pressure-pilot
   aggregate receipt and review packet on workers.
4. Perform a human pressure-choice review from the immutable four-slice
   calibration bundle. Do not invent or silently auto-select this science
   decision.
5. Materialize and validate a qualifying campaign plan only after the repaired
   launch boundary passes adversarial review.
6. Run the prerequisite slices, cheaper paper-test suites, the full Section 5.4
   reproduction campaign, and final publication analysis in staged order.
7. Finish optional extensions, final hardening, and the deferred external
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
/ccs/home/dfielding/athenak-pic
```

Repair base before the integrated ninth-pass commit:

```text
branch:     PIC
HEAD:       7324031b1cc50395c7903e73dbb153d4fdf71290
origin/PIC: 7324031b1cc50395c7903e73dbb153d4fdf71290
subject:    Run Q011 repair validation from a clean worker snapshot
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

Recent pushed history:

```text
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

These controller, candidate and policy bindings are retained historical
chronology. Do not relaunch them. Replace them only through the ninth-pass
strict-storage retirement, fresh worker freeze and candidate-only promotion
sequence documented below.

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
not expose a final receipt. Three Frontier Slurm worker publication attempts
then ended without exposing a public bundle:

| Slurm job | Terminal state | Elapsed | Reason |
| --- | --- | --- | --- |
| `4756951` | `FAILED` | `00:09:39` | Strict offline binary parser rejected runtime-added empty optional `particles/pic_deltaf_f0` |
| `4757026` | `CANCELLED` | `00:02:16` | Cancelled when the checkout mutated during publication; superseded by worker-local `git archive HEAD` snapshots |
| `4757047` | `CANCELLED` | `00:04:36` | Cancelled after fresh review found late-failure rollback, descriptor-anchor, and source-authorization gaps |

Aggregate publication remains prohibited until the integrated ninth-pass tree
is committed, pushed, clean-worker validated, installed and used to freeze a
fresh candidate. Provision the reviewed sibling acceptance authority before
rerunning the worker publisher. Do not invoke the publisher directly on a
login node.

```text
/usr/bin/sbatch tst/publication/frontier_q011_section54_pressure_pilot_publish_job.sh
```

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
node. Do not use a login-node direct invocation.

```bash
cd /ccs/home/dfielding/athenak-pic
/usr/bin/sbatch tst/publication/frontier_q011_section54_pressure_pilot_publish_job.sh
```

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

The live fixed sibling authority directory
`$PIC_ROOT/publication_acceptance/` is intentionally absent. Provision it once
with reviewed owner, group, mode and ACL immediately before worker
publication. The same-UID PIC-root mutation caveat and Lustre power-loss
durability caveat remain explicit operational risks.

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

After integrating the repairs, rerun focused tests, the explicit full suite,
`git diff --check`, and `python3 -B -m py_compile` for every changed Python
module before committing.

## New-agent restart checklist

Run these first:

```bash
cd /ccs/home/dfielding/athenak-pic
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
2. Preserve and inspect the in-flight repair edits.
3. Commit and push the integrated repaired boundary.
4. Run focused and full validation suites from the clean committed worker
   snapshot.
5. Request a fresh independent adversarial review.
6. Capture authenticated mirrored storage evidence, install the paired
   controller, promote the historical-slice retirement successor, freeze and
   independently revalidate a fresh worker-built candidate, and promote the
   candidate-only successor.
7. Provision `$PIC_ROOT/publication_acceptance/` through the reviewed one-time
   transition. Finish the pressure-pilot worker aggregate publication and
   review packet, then verify both receipts.
8. Ask the human collaborator to review the immutable four-slice pressure
   calibration evidence and select the qualifying pressure option.
9. Only then materialize the qualifying campaign plan and proceed through the
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
