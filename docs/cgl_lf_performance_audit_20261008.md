# CGL-LF turbulence cost audit, 2026-10-08

The exact archived scaling workload now measures **203.17 million
zone-cycles/s/node**, versus **103.22 million** for the original executable
repeated on the same nodes. This is whole-step throughput, not a kernel speedup:
the corrected LF timestep bound reduces work from 42 to 14 RHS evaluations per
outer step on that uniform, aligned-field workload. The turbulent checkpoint
requires about 146 evaluations with chunk cap 4 and is a different workload.

Detailed LF heat-flux diagnostics explain the largest measured avoidable part
of the turbulence run's cost. Disabling only those diagnostics gave **4.65×
faster outer timesteps**, preserving strict checks, safe arithmetic, weighted
fluxes, the chunk cap, and all saved endpoint fields at float32 precision.

The archived scaling result is 103.50 million zone-cycles/s/node, compared
with 3.730 million for the completed passive continuation: a 27.74× gap.
Both count outer cell updates, not internal LF stages. The relevant historical
source is `f2d0a259`, job 4925508, under `CGL/scaling/current_consistent_primary_f2d0a259`;
the older top-level scaling report describes a slower, earlier implementation.

## Correction: compiler settings differ

The first audit omitted a material build difference. The current application
and Kokkos libraries add **`-fno-cray -mno-daz-ftz`**, and linking adds
`-no-pie`; the historical scaling build used none of these. `-fno-cray` disables
Cray compiler enhancements, including GPU optimizations. A subsequent matched
full-build comparison measures only a **4.9% throughput gain** from removing
the corrective flags on the uniform workload with identical stage counts.
The turbulent checkpoint gains **9.1%**, but also changes its pressure
trajectory and takes 3,489 rather than 3,653 RHS evaluations over 25 steps;
that comparison does not isolate compiler kernel cost. Those flags do not
explain the original factor of thirty.
Runtime `cgl_lf_arithmetic=fast` is a separate switch and does not recreate the
historical build. Release/HIP/`-O3` matching alone does not establish compiler
equivalence.

WO2 introduced these flags for reproduced correctness failures, documented in
[Task 0](validation/wo2/README.md#task-0-establish-the-actual-gpu-baseline).
A CCE20 HIP minimum reduction returned a timestep 100.5 times too large on a
nonuniform density test; `-fno-cray` repaired it where `precise`/`strict` modes
did not. A separate CPU extreme-value test required the startup subnormal fix.
The [upstream Kokkos issue](https://github.com/kokkos/kokkos/issues/9476)
independently reports the minimum-reduction problem. This does not negate the
historical uniform scaling measurement.

Job **5638165** rebuilt the complete application and Kokkos independently for
both flag sets, from application revision `4af7c9a55` and Kokkos
`08ceff92bcf3a828844480bc1e6137eb74028517`. The historical flags reproduced the
100.5-fold timestep error; the corrective flags returned the analytic timestep.
Keep the corrective flags. Neither build explicitly enables unsafe floating
point atomics, and the diagnostic reduction is not an FP atomic-add loop.

Job **5638232** measures the archived executable and both fresh builds. Each
has one warmup invocation and two measured invocations in alternating order.
The uniform test uses the byte-identical archived 20-step input, two nodes,
16 GPUs and one 256³ block/GPU; timing spans cycles 5–20. The checkpoint test
uses one node, eight GPUs, 25 steps, no output writers and cycles 12070–12090.
Separate profile and endpoint-output runs are excluded from those timings.

| Workload / executable | Million zone-cycles/s/node, repeat mean | Repeat range |
| --- | ---: | ---: |
| Archived uniform / original executable | 103.217 | 103.155–103.279 |
| Archived uniform / current, corrective flags | 203.167 | 202.502–203.833 |
| Archived uniform / current, historical flags (incorrect reduction) | 213.189 | 213.179–213.199 |
| Turbulent checkpoint / current, corrective flags | 15.441 | 15.306–15.576 |
| Turbulent checkpoint / current, historical flags (incorrect reduction) | 16.850 | 16.848–16.853 |

The 1.968× current/original uniform throughput ratio is real for the stated
input but includes the change in internal work. The corrected build's uniform
heat-flux profile actually takes longer per RHS than the old executable. Do
not present these measurements as a twofold kernel optimization or infer
turbulence cost from a uniform box.

The stage reduction can be checked analytically. The archived built-in
turbulence pgen initializes `rho=1`, `p=1/gamma=0.6`, and `Bz=sqrt(0.12)`, with
uniform spacing `h=1/256`. Its parallel diffusivity is 0.19672783565. The old
isotropic bound, including hydro CFL 0.3, was `0.3 h²/(6 chi)`. The current
aligned face bound uses independent STS safety 0.9: `0.9 h²/(2 chi)`, nine
times larger. The same outer timestep therefore selects 21 versus 7 stages
per half-sweep. These counts were measured in separate profiles. See the
[LF bound definition](source/modules/cgl_lf_timestep.md) and
[WO2 bound qualification](validation/wo2/task1/README.md) for its scope and
nonlinear limitations.

## Direct same-checkpoint comparison

Job 5637965 used one node/eight GPUs, one 96×96×192 block per GPU, the same
time-4.976 checkpoint, binary `4af7c9a55`, and 25 outer steps per case. All
cases retained `dedt=0.16`, CFL 0.3, strict admissibility, finite limiters,
pressure-work recording, and LF chunk cap 4. The timing excludes five warmup
steps and the forced final-output interval, leaving 19 measured intervals.

| LF arithmetic / diagnostics / flux | Seconds/step | Million zone-cycles/s/node | Speedup |
| --- | ---: | ---: | ---: |
| Safe / full / weighted | 4.272 | 3.314 | 1.00× |
| Safe / none / weighted | 0.920 | 15.395 | 4.65× |
| Fast / none / weighted | 0.787 | 17.986 | 5.43× |
| Fast / none / physical | 0.771 | 18.353 | 5.54× |

Safe/full and safe/none have identical outer timestep schedules, internal work
counts, and saved fields. Fast modes keep the same outer schedule but change
internal work counts slightly and cause pressure differences: about 0.026% RMS
and up to 3.7% locally. These are short comparisons, not full-precision identity
or long-run accuracy claims. The minimal recommended setting is
**safe / none / weighted**, with strict checks and cap 4 retained. No production
settings were changed and no new production continuation was launched by this
audit. `none` disables detailed q-face/cap/q-work accounting, not admissibility
or separately enabled pressure/forcing-work ledgers. Suppressed q diagnostics
must be treated as unavailable.

A separate profile attributed 90.62 rank-mean seconds of 109.16 elapsed seconds
to LF heat-flux calculation, about 83%. The expensive q diagnostics run inside
that bucket. The timestep-bound reduction took 3.96 seconds; admissibility took
0.58 seconds. Profiling fences perturb timing, and nested buckets must not be
summed twice. This does not support disabling strictness for speed.

## Why the remaining comparison differs

The historical run used uniform unforced active CGL, PLM, no limiters/output,
fast/no-diagnostic/physical LF, and one 256³ block/GPU. Its retained timesteps
and stage-count rule imply 42 LF RHS evaluations per outer step. The current
long continuation measured 137.60, and the short audit about 146: **3.3–3.5×
more LF evaluations per outer step** from bounded LF subcycling. This is a work
count, not a measured independent wall-time factor.

Historical cells/GPU were 16.78 million versus 1.77 million now, a 9.48× ratio.
That ratio is not itself a measured speed penalty. Current LF timestep bounds
also inspect faces, transverse derivatives and drift terms and are recomputed
between chunks. Remaining effects include that work, layout/communication,
turbulent state, PPM, forcing, source changes and the compiler differences
above. Both builds use HIP/VEGA90A Release `-O3`; existing LF fast paths remain
present. These shared settings do not rule out a compiler performance effect.

## Evidence and reproduction

The audit completed successfully in **340 node-seconds (0.0944 node-hours)**.
All artifacts remain under the required WO2 test workspace:

`/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/physics-benchmark/matched/power032/lf-subcycling/performance-audit`

The [full report](/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/physics-benchmark/matched/power032/lf-subcycling/performance-audit/report.md)
contains exact repetition commands and primary historical links.
[Metrics](/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/physics-benchmark/matched/power032/lf-subcycling/performance-audit/metrics.json)
retain timing definitions, internal work counts, profile buckets, endpoint
differences, executable/restart hashes, and separate analysis provenance.

The passive continuation itself stopped cleanly at its wall limit at
time 4.976016897678366, cycle 12065, with a passing final-snapshot check and
zero recorded failure counters. The t=5 target and matched active/passive
experiment remain incomplete. The audit does not change that scientific status.

After review, the canonical matched input was changed to
`cgl_lf_diagnostics=none`, retaining safe arithmetic, weighted fluxes, and strict
checks. The measurements above precede that input-only change; their retained
inputs and artifacts remain unchanged. No new production run was launched.
The full-build and matched timing evidence is retained at
`performance-audit/../compiler-comparison`, including `comparison.json`, both
build manifests, and `timing/job-5638232/summary.json`. Each manifest records
actual compile/link commands, application/Kokkos identities, library and
executable hashes. The historical flags remain disqualified regardless of
their timing. Subsequent chunk experiments completed, but no larger setting
was adopted: short strict runs survived while local pressure comparisons were
not converged. The user stopped chunk tuning and speculative kernel work.
No experimental kernel patch was applied and no simulation is running.

The user explicitly rejected further full-diagnostic optimization. A partial
seven-value reducer experiment was stopped in job 5638265; its patch was never
applied to the repository and is not a proposed benchmark improvement. Retain
`cgl_lf_diagnostics=none`. The earlier full/none measurements identify the
configuration mistake behind the original slow run; they are not a matched
comparison to the user's diagnostic-free scaling workload.

The small [matched timing metrics](validation/cgl_lf_performance_20261008/matched-timing.json)
are retained in Git. Full builds/reduction controls cost 614 node-seconds;
matched timing and endpoints cost 469 seconds on two nodes (938 node-seconds).
Reproduce them with the retained `compiler-comparison/build-and-check.sbatch`
and `compiler-comparison/timing/run.sbatch` launchers, respectively. Their
bounded diagnostic runs do not automatically resume the physical benchmark.

An additional 127-second isolated compiler-probe matrix in job 5638265 tested
9,576 finite, exactly represented float/double extrema cases per compiler
variant plus the original density/timestep reproducer. `-fno-cray` passed all
of them. Default CCE20 and each narrower GPU/math/vector/unroll-disable flag
failed 6,080 extrema cases and returned the density timestep 100.5 times too
large. The defects affect Min, Max and MinMax. hipcc passed the generic extrema
checks, but its density probe did not compile because the standalone setup
lacked the MPI include path; it remains unqualified. No compiler-policy change
follows from this experiment. Evidence is under
`compiler-comparison/workaround-probes/runs-5638265-1791482860`.

## Current investigation: source changes against the pre-WO1 baseline

Compare actual code changes before proposing further settings or kernel work.
The historical scaling revision is older than the exact pre-WO1 baseline;
keep these source boundaries distinct:

- Historical scaling: `f2d0a25978e4fad9beae6874f279148c0460533b`.
- Pre-WO1: `8222de3aae4aaebd886653e7c61846e5f17987b4`.
- WO1 final: `543a1d7ec967a592f155fdd0d2213a5dfbe20b5c`.
- Pre-WO2: `9a4030b09307bc43865d5e597638f8645a6388f8`.
- WO2 release: `7a37710f6c224e24e7c7f364e7e0b812b3a9494c`.
- Latest solver change: `4af7c9a55cd5c5621f7d8a337c1ab132e4440968`.

The source audit is in progress. Existing throughput measurements above are
evidence for their specified workloads, not attribution of the passive
slowdown to a particular WO1/WO2 change. The LF chunk loop was introduced
after WO2 in `4af7c9a55`; it was not present in the historical scaling code.
