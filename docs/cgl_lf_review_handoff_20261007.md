# CGL-LF benchmark and stepping review handoff

Status checked at **2026-10-08 00:34 UTC (October 7, 20:34 Eastern)**. This is a review handoff, not a claim of numerical or physical validation. Development is paused pending review; no new simulations were launched for this handoff.

## Workspace and evidence locations

Development checkout: `/autofs/nccs-svm1_home2/dfielding/athenak-cgl-wo2`, branch `c/cgl-lf-physics-benchmark`, committed HEAD `359e814ac`. The original `/autofs/nccs-svm1_home2/dfielding/athenak-cgl` checkout is separate. All builds, tests, simulation outputs and scratch must remain under `/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2`.

In the paths below, define:

```text
M=/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/physics-benchmark/matched
D=$M/floor-diagnosis
```

Start with these retained reports:

- [Interrupted original experiment: reviewed report and figures](validation/cgl_lf_matched_physics_benchmark/interrupted_20261007/report.md). This directory also contains compact metrics, provenance and PNG/PDF figures.
- [Failure reproduction, timestep controls and costs](validation/cgl_lf_matched_physics_benchmark/forensics_20261007/report.md), with detailed evidence paths in adjacent `report.json`.
- [Performance logging and timing](validation/cgl_lf_matched_physics_benchmark/timing_20261007/report.md).
- [Matched benchmark guide](cgl_lf_matched_physics_benchmark.md).

Large original outputs: `$M/active-to14/` and `$M/passive-to14/`. Full original analysis: `$M/interrupted-analysis/`. Do not put large outputs in Git or overwrite retained simulation provenance.

## Intended experiment and current status

The user requested a matched active/passive turbulence comparison, with the earlier inexpensive run retained as a low-resolution reference. Canonical settings are 192×192×384, periodic box (1,1,2), mean B along z, initial rho=B0=1 and both pressures=5 (beta=10), PPM4/HLLE, RK2/RKL2, CFL=0.3, STS safety=0.9. Conservative A and the accepted B4 sharp-contact limitation must be preserved.

LF uses k_L=2π, background collision rate zero, soft anisotropy thresholds X=−2,+1 with X=2(p_perp−p_parallel)/B², finite limiter rate 10¹⁰, backups off. Forcing is OU with correlation time 2, seed 271828, physical shell [π,3π], k⁻² power and mixed projection. The amplitude blend `sol_fraction=1/(1+sqrt(2))` prescribes equal expected solenoidal/compressive innovation power. `dedt=0.32` is power per volume, hence total power 0.64; the user explicitly accepted this normalization difference from the paper. Do not retune forcing to prescribed Mach or magnetic fluctuation amplitudes.

Correct build route: `PROBLEM=built_in_pgens` with `problem/pgen_name=cgl_lf_paper`; the legacy `PROBLEM=cgl_lf_paper` is a different pgen route. Planned production target is t=14 and candidate averaging window [6,14], subject to actual stationarity and sampling review. Main comparison figures use dbfplot, mathematical labels, power-law eye guides, log PDF ordinates, and a **linear** magnetic-strength PDF abscissa per the latest user instruction. Pressure balance includes signed, scale-dependent statistics.

The original passive member failed at t≈2.5077. The original active member was then canceled at the user's request; last history t≈6.16005, last snapshot t≈6.00019. Reviewed startup comparisons [0.5,2.5] and active-only [2,6] show narrower active magnetic-strength distribution, suppressed field-aligned strain with appreciable parallel velocity, and pressure compensation at smaller resolved scales. These are encouraging directional signatures, but startup sampling, the passive failure and limited spectral range make the physical comparison **inconclusive**. Block ranges in figures are extrema of time-block means, not confidence intervals or independent-sample errors.

## Latest simulations: now finished, not currently running

At the status check, `squeue -u dfielding` was empty. Allocation **5631477** completed normally; the two large-block preflights stopped at their configured wall limits:

| Member | Directory under M | Final time | Final cycle | Solver wall duration |
| --- | --- | ---: | ---: | ---: |
| Active | `large-block-evolution-active/` | 1.632665 | 4441 | ≈3300 s |
| Passive | `large-block-evolution-passive/` | 2.139086 | 4653 | ≈4200 s |

Each used one node/eight GPUs, **one 96×96×192 block per GPU**, retaining the full 192×192×384 domain. Both returned zero but did **not** reach their requested t=3. They used the original uncapped RKL2 integration, not the proposed repair. The passive preflight did not reach the original failure time. Their clean stops therefore do not qualify a stepping fix.

Automatic analysis job **5631755** completed in 12 min 30 s on one node. Read:

```text
$M/large-block-preflight-analysis/attempts/20261007T170643.398424Z-5631755/report.md
$M/large-block-preflight-analysis/attempts/20261007T170643.398424Z-5631755/comparison/report.md
```

This analysis covers [0.5,1.5] in two 0.5-unit blocks. Its gate reports retained checks consistent, but **the new figures still require scientific review**. Analysis success is not physical validation. Simulation and analysis source revisions are recorded separately.

## What failed and what remains uncertain

The original passive failure is reproducible from a native checkpoint at cycle 5800, t=2.507714514270586. In pre-LF RKL2 stage 13/15, one interior cell reaches raw p_perp=−0.0087581473, while p_parallel=0.1827353911. At entry to that sweep, its pressures were p_parallel=0.1572801904 and p_perp=0.04455490248. An independent reconstruction of the multistage recurrence reproduces the negative pressure. This is finite positivity loss, not merely a roundoff-sized wall crossing, a plotting issue, or a stale magnetic-field cache.

The fatal `pfloor=1` matters; the accompanying `hard_bound=2254` counts intermediate wall crossings, which are allowed before scheduled wall treatment. Do not equate those counters with physical threshold occupancy or use floors to conceal a failed stage.

Changing rank/block ownership produced initially roundoff-sized pressure differences that later became large in a localized cold cell. A proper one-step 8-GPU versus 64-GPU comparison did not reveal interface-concentrated errors. This does not exclude a later implementation defect. The LF face collision rate includes threshold scattering, switching sharply between background zero and 10¹⁰; this is a plausible source of nonlinear sensitivity, not a demonstrated complete explanation. Formation of the depleted cell before the fatal stage remains an important concern.

Exact cold checkpoint: `$D/cold-state64/checkpoint5800/`, indexed by the forensic JSON. Restart SHA256: `004a60b4a94a7f5ac86cbf16cd91cd940af4d741db43aa45838421a8ba3160b5`. Original integrator fails from this state on both tested rank layouts. Controls and full-precision comparisons are under `$D/cold-probes-node1/`.

## How much the timestep was reduced

All tests retained **CFL 0.3**. RKL2 retained its separate **STS safety 0.9**. The original failing outer timestep was 3.08486590479×10⁻⁴. Identical-checkpoint controls advanced one original-step duration:

| Method | Initial outer timestep | Initial reduction factor | Outer steps to endpoint | Result |
| --- | ---: | ---: | ---: | --- |
| Original RKL2 | 3.08487×10⁻⁴ | 1 | — | Fails |
| RKL2, ratio cap 14 | 3.766734×10⁻⁵ | 8.19 | 8 | Positive, but large local discrepancy |
| RKL2, ratio cap 4 | 1.076210×10⁻⁵ | 28.66 | 25 | Positive, closer to explicit |
| Explicit LF | 8.968413×10⁻⁷ | 343.97 | 289 | Positive; refinement still required |

Timesteps evolve, so counts differ from initial reduction factors. `sts_max_dt_ratio=N` imposes `dt_outer <= N * (0.9 * dt_LF_raw)`; it is neither CFL=N nor a reduction by N. The explicit LF case instead uses CFL×its parabolic bound. These controls reduce the entire outer step, affecting hydro, forcing and collision cadence as well as LF. They are not isolated LF-only comparisons.

At the formerly failing cell, final (p_parallel,p_perp) are (1.46652,1.35381) for explicit, (1.44434,1.33163) for cap 4, and (1.02855,0.915847) for cap 14. Cap 4 differs locally by 1.5–1.6%; cap 14 by 30–32%. Whole-box relative pressure L2 differences are only about 2×10⁻⁴ and conceal the local discrepancy. The explicit result is a comparison, not yet a temporally converged truth.

Measured local costs on one node with the old 48³ block layout, excluding initialization/final output, are 6.296 node-hours per physical-time unit for cap 4, 2.916 for cap 14, and 50.572 for explicit. These short depleted-state measurements are not whole-run forecasts. No cap is qualified for production. The user chose **develop efficient, robust LF stepping first**, rather than spend many hours on a capped production pair.

## Precisely what passive means

The passive flow is **isothermal MHD**. Its dynamical pressure is

$$P_{\rm dyn}=c_{\rm iso}^{2}\rho=5\rho,\qquad c_{\rm iso}=\sqrt5.$$

Density, momentum and magnetic induction use native isothermal MHD fluxes, magnetic pressure/tension and forcing. The evolved CGL pressures do not supply momentum stress or determine MHD wave speeds. There is no active CGL total-energy budget for this isothermal flow.

Alongside it, two thermodynamic fields evolve passively. Their persistent conservative coordinates are

$$J=\rho\ln\!\left(\frac{p_\parallel B^2}{\rho^3}\right),\qquad
A=\rho\ln\!\left(\frac{p_\perp\rho^2}{p_\parallel B^3}\right).$$

They are advected using the native isothermal mass flux. Without LF or scattering, the continuum material invariants imply p_parallel B²/rho³ and p_perp/(rho B) remain constant along fluid trajectories. LF temporarily evolves actual thermal energy density U=p_parallel/2+p_perp and magnetic moment p_perp/|B|, then converts back. Conservative passive A remains part of the formulation.

Thus **p_perp is an evolved thermodynamic field, not merely a plotted quantity**. It enters its material invariant, perpendicular LF heat flux (temperature and magnetic-strength-gradient terms), anisotropy thresholds, scattering/pressure relaxation, threshold-dependent LF conductivity, LF parabolic timestep estimates, admissibility checks and diagnostics. Those operations modify the passive pressures but do not apply their stresses to momentum or induction. Their stiffness can still change numerical timing or abort the run. In particular, adding an outer LF timestep cap changes the numerical cadence of the passive flow despite the absence of physical pressure feedback.

Passive pressure balance must distinguish evolved tracer p_perp from actual momentum pressure 5rho; the analysis labels the latter separately. Passive pressure-work histories are diagnostic stress contractions, not work by a force applied to the isothermal fluid. KE+ME+passive U is not the active conserved total-energy budget.

Source entry points: `src/eos/cgl_passive.hpp`, passive branch of `src/mhd/rsolvers/hlle_cgl.hpp`, native isothermal prepass in `src/mhd/mhd_tasks.cpp`, passive wavespeeds in `src/mhd/mhd_newdt.cpp`, and face flux/timestep functions in `src/diffusion/cgl_landau_fluid.cpp`.

## Proposed fix and incomplete implementation

Candidate approach: adaptively subdivide LF half-sweeps internally while preserving the outer hydro/forcing step and original collision/wall cadence. Reject invalid trial states globally across MPI ranks and retry smaller chunks. Control local pressure discrepancies between a full chunk and two half chunks; positivity alone is inadequate. With discontinuous conductivity switches, this discrepancy is an empirical error indicator requiring convergence tests, not an automatically justified smooth second-order estimator. Do not assume division by three or use Richardson extrapolation without evidence.

Only an initial trial-handling layer is currently edited, **uncommitted and unbuilt**:

- `src/mhd/mhd_sts.cpp`
- `src/diffusion/cgl_landau_fluid.hpp`
- `src/diffusion/cgl_landau_fluid.cpp`

It introduces raw-pressure checks and retained trial-failure flags, suppresses immediate trial-time aborts so communication can drain, and bypasses repeated representation/collision hooks inside trial chunks. **Driver integration is missing**; these edits are not a runnable or validated fix. Root ownership was reserved for `src/driver/driver.cpp/.hpp`. Review the full diff before continuing.

Rollback must restore conserved/primitive fields including ghosts, representation, physical diagnostics/EOS counters, recurrence registers, caches and any modified timestep state. Rejected work must remain visible in computational cost but not be counted as accepted physical work. Intermediate anisotropy wall crossings alone must not trigger rejection. Wall/rate treatment must retain its original logical timing.

Design and validation notes:

```text
$D/integration-controller-options.md
$D/robust-stepping-validation-plan.md
$D/decomposition-source-audit.md
```

## Review priorities before a new production pair

1. Independently review the reproduced failure and whether the proposed controller addresses accuracy as well as positivity; investigate preceding cold-cell depletion.
2. Validate isolated LF evolution and coupled checkpoint replay with explicit timestep refinement, inspecting local pressure minima/errors as well as integrals. Compare cost at matched accuracy against LF-only explicit/SSPRK2; adaptive RKL overhead may erase its advantage.
3. Verify MPI rejection/rollback, conservation and diagnostic ledgers, restart behavior, unchanged outer collision cadence and passive flow identity at matched outer timesteps. Preserve existing damping checks and B4 expected limitations.
4. Review the automatically generated large-block preflight figures against actual data. They supply startup/layout evidence, not an accepted physical reference.
5. Only after a qualified numerical fix, start a fresh matched pair from t=0, retain exact input/overrides/source/build/seed/commands, and assess late-time sampling. Do not silently promote or change the lineage of old preflights.

Performance logging already reports interval wall seconds per completed cycle and global zone-cycles/s; it includes intervening communication/I/O and excludes STS stages from the cycle count. It was ported after examining the `scaling-tests` branch. An early matched-time comparison found the one-node large-block layout used about 47% fewer node-hours than the old eight-node layout, while taking longer wall time. This is useful efficiency evidence, not a controlled late-time scaling study.

Primary physical references: [Squire 2023](https://arxiv.org/html/2303.00468v2) and [MKS24](https://arxiv.org/html/2405.02418v2). Preserve distinctions between their forcing/limiter conventions and ours; Steve's apparent firehose cutoff near X=−1.4 is not our X=−2. An active/passive contrast does not establish resolution convergence or independently validate LF coefficients. Classify scientific evidence as consistent, concerning or inconclusive with reasons; do not invent tolerances from paper images.
