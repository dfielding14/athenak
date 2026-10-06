# Final WO2 integration and release checks

This record covers the combined implementation on Frontier allocation `5629018`
(`frontier06527` and `frontier09344`). It uses the build and runtime contract in
the [main report](../README.md). All raw work remains under
`/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2`.
The final allocation was released after all applications finished and reports
`COMPLETED` with exit status `0:0` in the retained scheduler accounting.

All selected implementation and validation work is complete. The final release
source matches the implementation worktree byte for byte across all 761
manifested source, test, input and build-configuration files. Documentation and
evidence added after the build are recorded separately.

## Combined acceptance

| Suite | Completed result | Qualification |
| --- | --- | --- |
| CPU CGL regressions | 235 passed, 2 strict B4 expected failures | Ordered union of the initial invocation and explicit corrective reruns. Six passive groups are recorded separately. |
| Broad HIP-launch CGL regressions | 235 passed, 2 strict B4 expected failures | Ordered union; self-building Serial helpers remain CPU coverage. |
| Dedicated GPU regressions | 52 passed | Actual HIP execution, including one/four-rank boundary, conservation, merge and AMR checks. |
| CPU MPI regressions | 20 passed | One/four-rank execution; passive CPU-MPI groups are separate. |
| HIP MPI regressions | 20 passed | Exact accepted GPU runtime, one/four-rank execution. |
| MPI GPU AMR | 2 passed | Four-rank restart through regridding and two-node 16-rank 3D churn. |
| GPU full workflow | 29 passed | Exact supplied HIP executable with `--no-build`. |
| GPU paper smoke | 3 passed | Both active cases and the newly released passive case. |
| Final passive package | 18 groups / 240 applications passed | CPU1, HIP1 and HIP4, including exact native-isothermal flow and restart comparisons. |
| Final release passive CPU-MPI | 6 groups / 80 applications / 65 checks passed | Four ranks, all five reconstructors, forced 3D flow and actual resumed full-state restart identity. |
| Final fusion comparison | 23 pairs passed | 66 snapshots and 2,018,316 face values agree bitwise; only two diagnostic reduction sums differ at rounding level. |
| Clean release equivalence | 25 CPU + 25 HIP pairs passed | 100 successful applications; exact fields, histories and normalized restart state, including all diagnostics and live RNG. |

Counts overlap where a focused rerun exercises an existing test. They must not
be summed into a new number of independent physical experiments. The aggregate
JSON files retain every observation and identify which final result replaces a
failed initial attempt. Skipped tests are not counted as passes. The unchanged
19 Stage I utility checks passed in the frozen baseline and were not rerun as
part of the final application suite.

The broader CPU-named suite also ran through a HIP application launcher. Tests
which configure their own Serial helpers remain CPU coverage. The
[CPU](evidence/final-research/cpu-aggregate.json),
[HIP-launch](evidence/final-research/hip-launch-aggregate.json),
[dedicated GPU](evidence/final-research/gpu-aggregate.json),
[CPU MPI](evidence/final-research/mpi-cpu-aggregate.json) and
[HIP MPI](evidence/final-research/mpi-hip-aggregate.json) aggregates retain the
original failed attempts and every corrective observation. The additional
[passive CPU-MPI report](evidence/task4-research/release-mpicpu-delivery/README.md)
records final-release coverage separately from the earlier passive package.

## Diagnosed failures and retained corrections

The first frozen-source test copy carried a relative Kokkos Git pointer that
became invalid after copying. Helper CMake configuration failed before any
numerical check. A separate test tree uses the canonical, same-revision Kokkos
symlink; the immutable application source and binaries remain retained. A
workflow unit test also required real Git metadata, so its unchanged assertions
were rerun from the implementation worktree with outputs under WO2.

The new timestep bound legitimately changes step sizes and counts. The
[Task1 followup](../task1-final-followup/README.md) preserves all original input
decks and scientific assertions while using existing STS ratio caps to retain
the intended refine/derefine sequence and the decay tests' minimum 100 cycles.
It also records the derived 13-stage heated post-sweep and zero stages for a
zero LF operator. The profiling test now expects fusion's aggregate timer;
the detailed directional test is unchanged.

The passive header added two list-initialization expressions that narrowed to
float. Explicit casts around the complete expressions restore the focused
float checks. The [compiler proof](../passive-cast-equivalence/README.md) shows
identical double CPU objects and HIP device code/constants. Moving the literal
pressure-traction helper to a shared header avoids importing unrelated GR
float errors into the previously supported CGL fast-speed check. The pressure
floor helper now links the actual Kokkos abort implementation. All double/float
focused checks pass, without a full-application float-portability claim.

Task4's first final four-GPU wrapper omitted parts of the accepted runtime
environment. Both fused and unfused executables intermittently faulted. Restoring
the complete contract produced six consecutive passing probes and a passing
80-application HIP4 suite. The [Task4 evidence](../task4/README.md) retains the
faults and exact environments. This establishes the verified runtime contract;
it does not identify the precise origin of the faulting host address.

## Clean release build

The clean release snapshot follows Task7 commit
`230b6e192cc53518b40795b54ddb4d0798def182`, with the archived formatting and test
followups. Both complete CPU and HIP builds pass. A read-only audit accounts
for all 340 existing numerical source files and one added helper header:
relative to the accepted fused snapshot,
changes are the two identity casts, literal helper relocation, and formatting
or comments. All changed C++ files pass the retained style checks; unrelated
pre-existing whole-tree diagnostics are separately recorded.

The initial CPU link inadvertently inherited a ROCm library dependency. Its
failed comparisons exited in the loader before solver startup. Relinking the
same objects under CPU-only modules removes that dependency; every object,
archive and the original failed binary remains unchanged and hashed. The
original failed executable is retained and is not the released CPU binary.

| Backend | Released executable under WO2 | SHA-256 |
| --- | --- | --- |
| CPU | `final-research/release-bin/athena-cpu` | `671441a4bbfb62f7192b7b6aaeaf8cc522603a84ea0482cab7f4d0408ab12d46` |
| HIP | `final-research/release-bin/athena-hip` | `6530c1b31a7d99077f3ff5932d4b1c790e33be38c1fa222cfda097f87b339f9e` |

Each backend's release matrix compares the accepted fused binary with this
clean build on 25 paired configurations, including forced passive 3D cases on
one/four ranks. It requires exact fields, histories and live restart state;
only the previously proven unused restart bytes are normalized. These are
ordinary evolved-state comparisons, not new face-trace coverage or performance
measurements. All 50 pairs pass: each backend has 50 exact field, 27 exact
history and 50 normalized restart file comparisons. Six recorded repair
counters are zero throughout. The [release report](evidence/final-research/RELEASE_EQUIVALENCE_REPORT.md)
and [summary](evidence/final-research/release-equivalence-summary.json) give
the exact restart normalization ranges, source provenance and retained loader
failure. The independent final-release CPU-MPI passive package additionally
validates all five reconstructors and actual resumed runs.

The final documentation build passes with warnings treated as errors, excluding
only the same two pre-existing orphan-page warnings. The user's dirty engineering
guide in the original checkout is preserved; implementation and evidence are in
the separate `c/cgl-lf-wo2` worktree.
An independent read-only review found no implementation blockers and confirmed
that current source matches the tested snapshot.

## Limits

The accepted B4 sharp-contact limitation remains explicit. The LF timestep
derivation bounds the local frozen-coefficient face Jacobian; it is not a global
nonlinear RKL2 or complete AMR stability theorem. Passive mode retains its stated
uniform-periodic and source restrictions, and has no irreversible shock heating.
Original paper forcing is planar with zero LF work; the separately labeled 3D
case supplies finite-time nonzero-transport evidence, not turbulence resolution
convergence. Individual task timings are not multiplied into a combined speedup.
Optional Tasks8–9 remain unselected.

The [artifact manifest](artifact-manifest.json) hashes the durable evidence
copied into this report. Raw field/restart files, failed attempts and complete
build trees remain under WO2 at the paths recorded in those artifacts.
