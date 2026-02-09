# AthenaK PIC Publication Readiness Plan (Step-by-Step)

This plan upgrades the current publication-priority figure set from internal
quick-look quality to publication-ready quality, with strict gates between each
step. It is aligned with the existing test inventory in:
`/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENT_PIC_TEST_PROBLEM_CATALOG.md`.

Primary constraints:
1. Preserve current regression/safety tests.
2. Add publication-physics runs as separate inputs/scripts where needed.
3. Enforce test and review gates between every implementation step.
4. Commit each completed step so rollback is trivial.

## 1) Alignment With Existing Catalog

This plan uses the catalog sections directly:
1. Section 3: Entity-derived PIC tests.
2. Section 4: Athena++ paper-inspired tests.
3. Section 6: publication-oriented priority set.
4. Section 7: immediate plotting batch.
5. Section 8: publication tooling and campaign assets.

Figure IDs used in this plan:
1. F01: entity deposit maps.
2. F02: EM vacuum convergence.
3. F03: Langmuir oscillation.
4. F04: two-stream and Weibel growth.
5. F05: Bell growth.
6. F06: multispecies oscillation and energy drift.
7. F07: CRSI polarization.
8. F08: CRPAI polarization.
9. F09: expanding/compressing anisotropy.

Current acceptance status baseline (updated from run
`tst/.codex/pic_publication_runs/20260208_115449`):
1. Publishable now: F01, F03, F09.
2. Publishable after styling only: F02, F04, F06.
3. Needs physics-run redesign (new input deck): F05, F07, F08.

## 2) Mandatory Gate Template (Used After Every Step)

Every step in this plan must end with the same gate workflow before moving on.

### Gate A: Correctness and Reproducibility
1. Re-run only affected tests first.
2. Re-run neighboring tests in the same physics family.
3. Confirm deterministic output under fixed seed and rank count.

### Gate B: MPI and Decomposition Safety
1. Run at least one serial (`np1`) and one MPI (`np2` or `np4`) case.
2. Check serial/MPI parity metric or bounded deviation metric.
3. Confirm no new warnings/fatals in run logs.

### Gate C: Tooling and Style
1. `python3 -m py_compile` on changed Python files.
2. `python3 -m flake8` on changed Python files.
3. Run changed-style checks as applicable.

### Gate D: Figure Review
1. Verify axis units, titles, legends, and annotations are physically correct.
2. Verify log/linear scales are intentional.
3. Verify fitted-window and uncertainty notes are visible where needed.

### Gate E: Change Review and Commit
1. Inspect diff file-by-file for unintended behavior changes.
2. Stage only files for the completed step.
3. Commit with a step-tagged message:
   `PIC-PUB Step <N>: <short summary>`.
4. Record run IDs and figure paths in step notes.

No step may start until the previous step passed Gates A-E.

## 3) Global Preflight (Step 0)

### Step 0.1: Baseline Snapshot
1. Record git branch/head and dirty-state summary.
2. Record current acceptance decisions for F01-F09.
3. Record baseline run ID:
   `tst/.codex/pic_publication_runs/20260208_102612`.

### Step 0.2: Environment Snapshot
1. Record `./athena -c` output.
2. Record MPI launcher path and version.
3. Record Python environment and plotting package versions.

### Step 0.3: Baseline Repeatability Check
1. Re-run `entity_core,benchmarks` once with MPI on.
2. Confirm identical pass/fail statuses and comparable metrics.
3. Archive as `baseline_repeatability` run ID.

### Step 0 Gate Tests
1. `python3 tst/publication/run_pic_publication_suite.py --groups entity_core,benchmarks --tier local --enable-mpi --mpiexec /opt/homebrew/bin/mpiexec --keep-going`
2. `python3 tst/publication/plot_pic_publication_figures.py --run-dir <run_id> --bundle all`
3. Confirm all cases `status=ok` in `summary.json`.

### Step 0 Deliverables
1. Baseline snapshot markdown note under `tst/.codex/pic_publication_runs/<run_id>/`.
2. Commit metadata-only notes if needed.

## 4) Styling-Only Upgrade Track (F01, F03, F09)

These steps keep physics inputs unchanged and improve publication presentation,
traceability, and figure semantics.

## 4.1) Step 1: F01 Entity Deposit Publication Panel

### Objective
Upgrade F01 into a clean methods-quality panel with explicit decomposition
parity evidence and neutral-species context.

### Implementation
1. Keep primary panel as `pdens` support map.
2. Add compact parity inset/table:
   `max|npN-np1|` for available ranks.
3. Add charge-neutrality text:
   `sigma(rho_prtcl)` and total charge integral.
4. Use publication colormap and fixed colorbar units.
5. Add panel labels `(a)` Minkowski and `(b)` Reflect.

### Required Tests
1. `particles.pic_entity_deposit_mink` in `np1,np2,np3,np4`.
2. `particles.pic_entity_deposit_reflect` in `np1,np2,np4`.
3. Verify parity metrics remain exact or machine-level.

### Step 1 Review Checklist
1. Does panel communicate support and parity, not false charge structure?
2. Are all annotations physically meaningful and concise?
3. Is color normalization consistent across reruns?

### Step 1 Exit Criteria
1. Figure visually clear and methods-appendix ready.
2. No regression test failures.
3. Commit completed.

## 4.2) Step 2: F03 Langmuir Publication Panel

### Objective
Keep current data but make the figure publishable with analytic context.

### Implementation
1. Add analytic sinusoid overlays using expected frequency.
2. Add residual panel or inset for phase error vs time.
3. Report extracted frequencies for `np1` and `np2`.
4. Use normalized time axis if required by target paper format.
5. Keep raw integrated-current panel available in supplementary output.

### Required Tests
1. `particles.pic_langmuir_frequency_proxy` in serial and MPI.
2. Verify extracted frequencies within existing thresholds.
3. Verify serial/MPI phase drift remains bounded.

### Step 2 Review Checklist
1. Do measured and analytic curves separate clearly if any error exists?
2. Is the frequency-fit method documented in caption text?
3. Is visual clutter minimized?

### Step 2 Exit Criteria
1. Main panel + residual panel are publication quality.
2. Existing proxy test still passes unchanged.
3. Commit completed.

## 4.3) Step 3: F09 Anisotropy Publication Panel

### Objective
Promote anisotropy trend plot to publishable by adding theory context and
uncertainty.

### Implementation
1. Keep startup-point filtering (`Jy` near zero excluded).
2. Add fitted trend lines and confidence intervals.
3. Add theoretical reference trend from chosen expansion model.
4. Add MPI overlay and slope parity annotation.
5. Produce two variants:
   one for main text, one with full diagnostic annotations for supplement.

### Required Tests
1. `particles.pic_expanding_box_anisotropy_proxy` in `np1,np2,np4`.
2. Confirm slope sign and separation checks still pass.
3. Confirm MPI slope differences stay under thresholds.

### Step 3 Review Checklist
1. Are sign and magnitude claims explicitly justified in-figure?
2. Is model assumption noted in caption metadata?
3. Are fitted windows clearly indicated?

### Step 3 Exit Criteria
1. F09 moves to "publishable now" quality.
2. No proxy regressions introduced.
3. Commit completed.

## 5) Physics-Run Redesign Track (F02, F04, F05, F06, F07, F08)

These figures need new dedicated publication input decks and analysis logic.
Existing proxy decks stay as regression gates and are not deleted.

## 5.1) Step 4: Create Publication Deck Naming and Manifest Pattern

### Objective
Introduce a stable naming convention so publication physics runs do not alter
proxy regression behavior.

### Implementation
1. New input decks use suffix `_publication`.
2. New analysis scripts use suffix `_publication`.
3. Manifest groups split into:
   `entity_core_proxy`, `entity_core_publication`,
   `benchmarks_proxy`, `benchmarks_publication`.
4. Plotting script consumes only publication groups for final pack.

### Required Tests
1. Existing proxy groups still runnable and unchanged.
2. New publication groups list correctly with `--list`.

### Step 4 Exit Criteria
1. Tooling split is stable.
2. No renamed/deleted proxy assets.
3. Commit completed.

## 5.2) Step 5: F02 EM Vacuum Convergence Redesign

### Objective
Produce a true convergence figure with enough points to support measured order.

### Implementation
1. Add publication deck with resolution sweep:
   `N = 16, 24, 32, 48, 64` (local) and `N = 96` (HPC optional).
2. Add script to fit observed order `p` from `L1 ~ N^{-p}`.
3. Run both serial and MPI at each resolution.
4. Add parity plot panel:
   `|L1_np1 - L1_np2| / L1_np1` vs `N`.

### Required Tests
1. Existing `particles.pic_em_vacuum_wave` passes unchanged.
2. New publication script must pass:
   `p` within target band and serial/MPI mismatch below bound.

### Review Gate Criteria
1. Minimum 5 resolution points used in fit.
2. Fit excludes under-resolved outlier only with logged rationale.
3. Figure caption records fit window and confidence interval.

### Step 5 Exit Criteria
1. F02 classification improves to "publishable after styling only".
2. Commit completed.

## 5.3) Step 6: F04 Two-Stream and Weibel Growth Redesign

### Objective
Replace stability-guard curves with physically growing instability runs.

### Implementation
1. Keep current proxy decks/scripts unchanged for regression.
2. Add new publication decks:
   one for two-stream growth regime, one for Weibel growth regime.
3. Add seeded perturbation controls to make linear phase measurable.
4. Increase output cadence in linear window.
5. Add script for automatic linear-window detection and growth fit.
6. Add comparison against expected growth trend from selected theory model.

### Required Tests
1. Proxy scripts still pass unchanged.
2. Publication scripts pass growth criteria:
   positive gamma with robust fit (`R2` threshold), serial/MPI consistency.
3. Rank sweep at least `np1,np2` and optional `np4`.

### Review Gate Criteria
1. Linear window is not hand-picked silently.
2. Saturation region excluded from main growth fit.
3. Error bars or fit uncertainty shown.

### Step 6 Exit Criteria
1. F04 becomes scientifically interpretable growth evidence.
2. Commit completed.

## 5.4) Step 7: F05 Bell Growth Redesign

### Objective
Produce a clear Bell-growth panel with a robust linear-growth fit.

### Implementation
1. Create `pic_bell_growth_publication.athinput`:
   longer `tlim`, higher dynamic range, controlled CR loading.
2. Add parameter sweep over CR loading and resolution.
3. Enforce coupled vs uncoupled pair for each sweep point.
4. Fit growth in linear phase and report uncertainty.
5. Add MPI parity overlays (`np1,np2` at minimum).

### Required Tests
1. Existing `particles.pic_bell_growth_proxy` remains unchanged and passing.
2. Publication Bell script checks:
   uncoupled near-flat, coupled positive gamma,
   fit quality above threshold,
   serial/MPI gamma agreement within tolerance.

### Review Gate Criteria
1. `R2` sufficiently high for claimed window.
2. Amplification factor high enough for publication signal quality.
3. Sweep picks one canonical "main text" run plus supplementary scan panel.

### Step 7 Exit Criteria
1. F05 leaves redesign class.
2. Commit completed.

## 5.5) Step 8: F06 Multispecies Oscillation Redesign

### Objective
Obtain physically consistent oscillation and controlled drift across mesh modes.

### Implementation
1. Add publication variants for uniform/SMR/AMR cases with harmonized setup.
2. Increase particles-per-cell and tune timestep/CFL for drift control.
3. Add phase and frequency error extraction relative to analytic oscillator.
4. Separate amplitude panel and drift panel with confidence bands.
5. Add rank parity checks for each mesh mode.

### Required Tests
1. Existing proxy script still passes.
2. Publication script verifies:
   frequency accuracy bounds,
   phase error bounds,
   energy drift bounds,
   serial/MPI consistency per mesh mode.

### Review Gate Criteria
1. Drift curves do not dominate physics interpretation.
2. AMR/SMR differences are explained and quantified.
3. No hidden normalization tricks in presentation.

### Step 8 Exit Criteria
1. F06 suitable for publication benchmarking panel.
2. Commit completed.

## 5.6) Step 9: F07 CRSI Redesign

### Objective
Produce rank-robust CRSI branch-growth evidence with clear dominant mode.

### Implementation
1. Add publication CRSI deck with tuned background and distribution params.
2. Run both delta-f and full-f publication variants if needed.
3. Add ensemble runs over random seeds for noise characterization.
4. Fit growth rates for left/right branches with uncertainty.
5. Add branch dominance ratio evolution panel.

### Required Tests
1. Existing `pic_crsi_deltaf_proxy` remains a fast gate.
2. Publication CRSI script checks:
   branch dominance stability across ranks and seeds,
   gamma agreement tolerances,
   acceptable fit quality.

### Review Gate Criteria
1. Dominant branch should not flip unpredictably across rank/seed.
2. If flips persist, parameter regime is rejected.
3. Final chosen regime has documented physical rationale.

### Step 9 Exit Criteria
1. F07 moves to "publishable after styling only" or better.
2. Commit completed.

## 5.7) Step 10: F08 CRPAI Redesign

### Objective
Get robust prolate/oblate branch-selection behavior with rank stability.

### Implementation
1. Add publication prolate and oblate decks with stronger branch separation.
2. Run rank and seed matrices matching CRSI methodology.
3. Quantify branch dominance confidence intervals.
4. Add side-by-side panel with consistent axis scales and fit windows.

### Required Tests
1. Existing `pic_crpai_polarization_proxy` remains unchanged and passing.
2. Publication CRPAI script checks:
   expected branch preference for prolate/oblate,
   rank/seed robustness thresholds,
   fit quality requirements.

### Review Gate Criteria
1. Opposite branch preference must persist across accepted runs.
2. Outlier seeds/ranks must be reported, not hidden.
3. Final panel includes uncertainty or spread.

### Step 10 Exit Criteria
1. F08 exits redesign class.
2. Commit completed.

## 6) Campaign Integration and Final Figure Production

## 6.1) Step 11: Manifest and Automation Integration

### Objective
Integrate all publication variants into `tst/publication` workflow while
retaining fast proxy suite.

### Implementation
1. Add publication case entries to manifest.
2. Add group-level switches:
   `--groups entity_core_publication,benchmarks_publication`.
3. Keep proxy and publication runs independently executable.
4. Ensure output archive structure clearly separates both classes.

### Required Tests
1. Run proxy-only suite.
2. Run publication-only suite.
3. Run mixed suite and confirm no collisions in output filenames.

### Step 11 Exit Criteria
1. Stable run orchestration for both modes.
2. Commit completed.

## 6.2) Step 12: Figure Styling Finalization Pass

### Objective
Apply publication styling uniformly once physics content is accepted.

### Implementation
1. Standardize typography, line weights, and panel lettering.
2. Standardize color palette by figure family.
3. Add consistent annotation templates:
   fit window, gamma, confidence interval, rank parity note.
4. Export at required DPI and vector formats when possible.

### Required Tests
1. Regenerate full publication figure pack.
2. Verify no data logic changed during styling-only edits.
3. Confirm script-level hash or checksum consistency for data arrays.

### Step 12 Exit Criteria
1. All accepted figures pass visual quality checklist.
2. Commit completed.

## 6.3) Step 13: Final Acceptance and Lock

### Objective
Produce the final acceptance report and lock figure versions.

### Implementation
1. Run full publication suite matrix with MPI enabled.
2. Generate final figure pack path and immutable run ID.
3. Write `FINAL_PUBLICATION_FIGURE_ACCEPTANCE.md` with:
   per-figure pass/fail, run IDs, fit metrics, and residual risks.
4. Tag figure pack manifest with commit hash.

### Required Tests
1. Full publication suite local matrix.
2. Optional HPC reruns for heavy or high-resolution cases.
3. Spot-check reproducibility by rerunning one figure per family.

### Step 13 Exit Criteria
1. Every figure classified either:
   publishable now,
   publishable after styling only,
   or deferred with explicit reason.
2. Final lock commit completed.

## 7) Detailed Test Matrix by Figure

Use this matrix during Steps 1-10.

1. F01 Entity deposit:
   `np1,np2,np3,np4` for Minkowski, `np1,np2,np4` for Reflect.
2. F02 EM vacuum:
   `N=16,24,32,48,64` for `np1,np2`; optional `N=96` HPC.
3. F03 Langmuir:
   baseline resolution in `np1,np2`, plus one higher resolution in `np1`.
4. F04 two-stream/Weibel publication:
   at least `np1,np2`, plus one resolution scan.
5. F05 Bell publication:
   coupled/uncoupled for `np1,np2`, CR-loading scan.
6. F06 multispecies publication:
   uniform/SMR/AMR each in `np1,np2`.
7. F07 CRSI publication:
   `np1,np2,np4` with at least 3 seeds.
8. F08 CRPAI publication:
   prolate/oblate in `np1,np2,np4` with at least 3 seeds.
9. F09 anisotropy:
   expanding/compressing in `np1,np2,np4`.

## 8) Required Review Artifacts Per Step

For each completed step, archive:
1. Run command and run ID.
2. Case status table from `summary.json`.
3. Key numeric metrics table (fit slopes/gamma/error/drift).
4. Figure path list.
5. Reviewer notes:
   what changed, what improved, what remains.

Store artifacts under:
`tst/.codex/pic_publication_runs/<run_id>/reviews/step_<N>/`.

## 9) Commit and Branch Discipline

1. One step per commit.
2. Commit message template:
   `PIC-PUB Step <N>: <scope>`.
3. Include run ID(s) in commit body.
4. Do not combine physics redesign with styling in one commit.
5. Do not modify proxy thresholds unless redesign explicitly requires it.

## 10) Stop Conditions and Escalation Rules

Stop and escalate before next step if any of these occur:
1. Existing regression proxy begins failing after unrelated edits.
2. Rank-parity behavior regresses without clear physical cause.
3. Fit quality requires manual cherry-picking of windows.
4. Publication claim cannot be made from current parameter regime.

Escalation output must include:
1. exact failing metric,
2. run ID,
3. suspected root cause,
4. minimal proposed corrective action.

## 11) Final Mapping Back to Catalog

This plan is catalog-compatible because it:
1. keeps Section 3 and Section 4 proxy tests as regression anchors,
2. adds publication variants without replacing existing IDs,
3. reuses Section 8 tooling (`manifest`, `runner`, `plotter`),
4. directly executes the Section 6 publication priority set in staged form.

## 12) Execution Status Snapshot (2026-02-08)

### 12.1 Completed Infrastructure Steps
1. Added publication-specific decks under `inputs/tests/*_publication.athinput`
   for two-stream, Weibel, Bell, multispecies, CRSI, and CRPAI.
2. Added publication-specific analysis modules under
   `tst/scripts/particles/*_publication.py`.
3. Added windowed exponential fitting utility in
   `tst/scripts/particles/pic_analysis_utils.py`.
4. Split manifest execution groups into proxy and publication families in
   `tst/publication/pic_publication_manifest.py`.
5. Updated run harness defaults and publication plotting to prefer publication
   artifacts with proxy fallback:
   - `tst/publication/run_pic_publication_suite.py`
   - `tst/publication/plot_pic_publication_figures.py`
6. Added metrics extraction utility:
   `tst/publication/compute_pic_publication_metrics.py`.

### 12.2 Latest Validated Run IDs
1. Full local proxy+publication matrix:
   `tst/.codex/pic_publication_runs/20260208_115449`
   - `overall_status = ok` for all 22 selected cases.
2. Focused publication-threshold regression rerun:
   `tst/.codex/pic_publication_runs/20260208_115405`
   - `overall_status = ok` for Bell/MSO/CRPAI publication modules.

### 12.3 Strict Figure Acceptance (from run `20260208_115449`)
1. F01 `01_entity_deposit_maps.png`:
   `publishable now`.
   - Rationale: exact decomposition parity annotations and clean support maps.
2. F02 `02_em_vacuum_convergence.png`:
   `publishable after styling only`.
   - Rationale: 5-point convergence and MPI overlap are good; title still says
     "proxy" and should be publication-labelled.
3. F03 `03_langmuir_currents.png`:
   `publishable now`.
   - Rationale: analytic overlays plus residual panel and rank overlays.
4. F04 `04_two_stream_weibel_modes.png`:
   `publishable after styling only`.
   - Rationale: growth is visible and fitted; panel still needs cleaner fit
     window annotation and axis-label polish.
5. F05 `05_bell_growth.png`:
   `needs physics-run redesign (new input deck)`.
   - Rationale: current deck shows early-time impulse then decay; final
     amplification is ~1.01 and fitted gamma sign depends on window treatment.
6. F06 `06_multispecies_oscillation.png`:
   `publishable after styling only`.
   - Rationale: frequency and drift are controlled and rank-consistent; visual
     normalization/caption text needs publication cleanup.
7. F07 `07_crsi_polarization.png`:
   `needs physics-run redesign (new input deck)`.
   - Rationale: branch dominance flips across rank (`np1 right`, `np2/np4 left`)
     indicating regime is not rank-robust for publication claim.
8. F08 `08_crpai_polarization.png`:
   `needs physics-run redesign (new input deck)`.
   - Rationale: prolate/oblate branch preference is not consistently separated
     across rank; asymmetry remains regime-sensitive.
9. F09 `09_expanding_box_anisotropy.png`:
   `publishable now`.
   - Rationale: opposite-signed slopes with uncertainties and clean trend fits.

### 12.4 Metrics Anchors (run `20260208_115449`)
1. F02 order: `np1_order=np2_order=1.7036`.
2. F04 growth fits:
   - two-stream: `gamma_np1=6.217`, `R2_np1=0.880`,
     `gamma_np2=4.422`, `R2_np2=0.917`.
   - Weibel: `gamma_np1=6.965`, `R2_np1=0.924`,
     `gamma_np2=3.916`, `R2_np2=0.967`.
3. F05 Bell:
   - uncoupled: ratio `1.000`, gamma ~`0`.
   - coupled: ratio `1.010`, fitted early-window gamma `33.15`
     (not publication-stable due to impulse/decay profile).
4. F06 drifts:
   - uniform `1.09e-2`, SMR `2.03e-2`, AMR `2.86e-2`.
5. F09 slopes:
   - expanding `-0.2314`, compressing `+0.4519`, rank overlays consistent.

### 12.5 Remaining Work to Reach Full "Publishable Now"
1. F05 Bell:
   redesign deck toward sustained linear growth (reduce startup impulse,
   increase linear-phase duration, tune CR loading/background).
2. F07 CRSI:
   lock a rank-stable branch-dominant parameter regime; likely requires seed
   control and/or stronger branch separation.
3. F08 CRPAI:
   tune prolate/oblate setups for robust opposite branch preference across
   rank and seeds.
4. Styling-only polish for F02/F04/F06 after physics locks above.
