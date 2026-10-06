# Task 7: fused LF face kernels

Triage: **implemented**. Three directional face launches become one padded-cell
launch while each copied directional body retains its exact arithmetic.

The candidate combines directional LF face kernels on a padded cell lattice.
Each thread evaluates its lower x/y/z face when that face exists; guards include
the high caps and precede every stencil load. Separate lexical scopes retain
each direction's original arithmetic. The generator copies and hashes all six
full/none directional bodies without changing expressions, face-field choice,
safety arithmetic, limiter decisions, or diagnostic face ownership.

The independently timed candidate uses Task1+Task0-refill source. It adds no class/header layout or
runtime switch. Existing detailed profiling retains the original directional
path. Ordinary fused execution has a named Kokkos region and aggregate
heat_flux_total timing; the old directional buckets are unused in that path.
CCE reported occupancy-target warnings. Hardware resource behavior must be
judged from measurements rather than assumed load-traffic savings.

## Standalone correctness

Seventeen actual-HIP configurations cover 1D/2D/3D, safe/fast arithmetic,
full/none diagnostics, local/background coefficients, unequal spacing, reversal,
checkerboard, and multilevel stencils, including four-rank 2D/3D runs. Two
below-floor cases execute zero LF stages and agree in all state outputs. The
other fifteen compare 42 snapshots with 1,845,516 evolved face values having
exact uint64-bit identity, including signed zero.

All binary fields and every physical restart byte agree. All history columns
except the floating qpar_work/qperp_work sums agree exactly. Fusion changes the
diagnostic reduction order: five final restart files differ only in those two
named doubles, with maximum relative difference 3.60189268286223e-16. The
comparison bound is 256 double epsilons relative and zero absolute allowance.
Every byte outside those fields and the existing 36-byte unused-root-index
normalization must match exactly. Raw hashes, exact differing offsets, named
field values, serialization source lines, and complete audits are retained.
The diagnostic bytes are inspected, not normalized.

Earlier harness failures remain available: entirely below-floor cases produce
no LF traces because they have zero LF stages, and reversal generators require
at least 20 cycles rather than the generic four-cycle check. Corrected checks
retain those requirements; no scientific threshold was loosened.

## Standalone unprofiled GPU timing

Allocation 5628672 was reserved for these measurements. Each case used one
warmup per binary, three alternating paired repeats, 40 cycles, identical saved
inputs, LF profiling disabled, and initial/final field outputs. Solver-reported
elapsed time starts after initialization and initial outputs; it includes the
evolution loop, evolved/final outputs, and final problem-generator validation.
Launch-inclusive elapsed is recorded separately. Every pair passes the same state/counter checks.

| Case | Arithmetic / diagnostics | Ranks | Unfused ms/cycle | Fused ms/cycle | Speedup |
| --- | --- | ---: | ---: | ---: | ---: |
| lf2d | safe-full | 1 | 42.318 | 35.247 | 1.2006 |
| lf2d | safe-full | 4 | 42.848 | 35.632 | 1.2025 |
| lf2d | fast-none | 1 | 27.727 | 26.592 | 1.0427 |
| lf2d | fast-none | 4 | 26.326 | 26.202 | 1.0047 |
| lf3d | safe-full | 1 | 90.338 | 59.947 | 1.5070 |
| lf3d | safe-full | 4 | 63.832 | 43.240 | 1.4762 |
| lf3d | fast-none | 1 | 34.558 | 32.603 | 1.0600 |
| lf3d | fast-none | 4 | 33.319 | 31.545 | 1.0562 |

Default safe/full configurations show about 20% gain in 2D and 48–51% in 3D
here. Fast/none improvements are smaller and partly overlap run-to-run spread.
The safe/full and fast/none measurements change both arithmetic and diagnostic
mode; they do not isolate the cause of the speedup. No memory-bandwidth or
scaling mechanism is claimed from these measurements.

## Final integrated acceptance

The fusion was regenerated from the immutable combined unfused source after
Tasks 0–6. All 1,282 source files agree except the intended insertion in the LF
translation unit. The complete LF TU and all six extracted directional bodies
also agree with the earlier measured composition; final compilation uses the
integrated headers and other TUs. Both CPU and HIP links verify that every reused
full-build object/archive stayed unchanged. Scratch trace instrumentation adds
only a header and a post-flux snapshot to the final STS TU, preserving passive
canonical pressure refreshes and all other lifecycle statements.

All 23 actual-GPU paired functional checks passed in allocation 5629018: the
17 active configurations plus six passive checks covering 2D/3D, safe/fast,
full/none diagnostics, collisions, and one/four ranks. Across 66 snapshots,
2,018,316 evolved face values agree bit for bit. Every physical restart byte,
binary field, and non-work history column agrees, and all six recorded safety
counters remain zero. Seven final restart files differ only in the two permitted
q-work diagnostic sums; the largest relative discrepancy is 5.606513348085635e-16,
within the unchanged 256-epsilon relative / zero-absolute allowance. The same
36 unused root-index bytes are the only restart normalization; diagnostic bytes
are inspected individually, not normalized.

These final-composition runs are functional checks with tracing, not timing
measurements. The raw runner JSON retains an unused generic timing-template
description; its phase and one-pair-per-case records are explicit. The eight
speedups above apply only to the separately measured Task0+Task1+Task7 composition.
No new final-composition speedup or whole-work-order speedup is inferred.

The ordinary profiling regression now expects `heat_flux_total`, since the
fused path reports aggregate timing. The detailed directional regression is
unchanged. Both tests pass on final CPU and HIP binaries. Their first CPU launch
failed before solver initialization because the scratch runner requested GPU
MPI support without GTL; corrected backend-specific environment and all original
failure logs are retained. The immutable numerical source archives were kept,
and a separate test overlay validates this one-token expectation change.

## Durable evidence

- [Final functional summary](final-correctness-summary.json) contains every case,
  exact face counts, binary/input hashes, and named diagnostic-only differences.
- [Standalone correctness](standalone-correctness-summary.json) preserves the
  initial acceptance; [timing summary](timing-summary.json) retains all eight
  configurations and all 64 warmup/measured application records.
- [Source/build provenance](source-build-provenance.json) records extraction,
  face ownership, source comparison, compile/link commands and binary hashes.
- [Profiling regression](profiling-regression-summary.json) records both CPU/HIP
  checks and the one-token test expectation correction.
- [Artifact index](artifact-index.json) hashes the complete original results,
  logs, source patches and tools under WO2. Raw failures and outputs remain intact.
- `analysis/` preserves the source generator, ownership checker and runtime audit.
