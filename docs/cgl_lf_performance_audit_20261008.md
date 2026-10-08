# CGL-LF turbulence cost audit, 2026-10-08

Detailed LF heat-flux diagnostics explain the largest avoidable part of the
turbulence run's cost. Disabling only those diagnostics gave **4.65× faster
outer timesteps**, preserving strict checks, safe arithmetic, weighted fluxes,
the chunk cap, and all saved endpoint fields at float32 precision.

The archived scaling result is 103.50 million zone-cycles/s/node, compared
with 3.730 million for the completed passive continuation: a 27.74× gap.
Both count outer cell updates, not internal LF stages. The relevant historical
source is `f2d0a259`, job 4925508, under `CGL/scaling/current_consistent_primary_f2d0a259`;
the older top-level scaling report describes a slower, earlier implementation.

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
turbulent state, PPM, forcing and source changes. Both builds use optimized
HIP/VEGA90A Release `-O3`; existing LF fast paths remain present.

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
