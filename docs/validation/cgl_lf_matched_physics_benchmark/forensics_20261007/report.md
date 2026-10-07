# Passive LF pressure-floor forensics, 2026-10-07

**Numerical behavior is concerning; the root cause and physical validation
remain inconclusive.** The original passive member failed strict checking
because an intermediate LF update produced a finite negative perpendicular
pressure. A replay reproduced it. A different MPI layout followed a different
thermal trajectory and did not fail through $t=2.6$; that is not a qualified
fix. Subsequent controls from the actual cold checkpoint traverse the failed
step with capped STS and explicit LF. Cap 4 agrees more closely with explicit
integration than cap 14, but is substantially more expensive than the original
scheme. No guard is qualified for a full turbulent run.

This report concerns the original $192\times192\times384$, PPM4/HLLE,
RK2/RKL2 calculation with $48^3$ blocks, STS safety 0.9 and no outer
STS-ratio cap. The retained production revision is
`ab9b543e7e6a972ebfc026d0c98ea3dfee9cb55b`; exact executable and analysis
hashes, measurements and raw-evidence paths are in [report.json](report.json).
The [benchmark guide](../../../cgl_lf_matched_physics_benchmark.md) defines the
physical experiment. Large outputs remain outside Git.

## Reproduced failure

The original 64-rank/eight-node diagnostic replay, starting from the retained
$t\simeq2$ checkpoint, fails at cycle 5800,
$t=2.5077145142705861$, pre-sweep stage 13/15. The violating cell is an
interior cell at $(0.3671875,0.8255208333,0.5807291667)$, global block 51.

| Raw stage quantity | Value |
| --- | ---: |
| Thermal energy density $U$ | 0.0826095483 |
| Magnetic moment $μ=p_\perp/\lvert B\rvert$ | -0.0260856442 |
| $p_\parallel$ | 0.1827353911 |
| $p_\perp$ | **-0.0087581473** |
| Cached and independently computed $\lvert B\rvert$ | 0.3357458681, identical |

An independent reconstruction of the five RKL2 contributions reproduces $U$
exactly at retained precision and $\mu$ within $2.1\times10^{-17}$.
The initial, stage-11 and stage-12 perpendicular pressures are positive:
0.04455, 0.03651 and 0.01538. The current weighted flux contribution increases
$p_\perp$ by 0.003361, but the complete multistage combination makes it
negative. This establishes finite intermediate-stage positivity loss, not
subtraction roundoff, a ghost-cell artifact or a cached-field discrepancy.
It does not establish why the preceding stage states and spatial fluxes arose.

The log's `pfloor=1` triggers the abort. Its `hard_bound=2254` records
intermediate wall crossings, which are permitted until the scheduled wall
projection. In general floor counters include refreshed ghosts and the other
admissibility counters inspect interior cells after repair; the instrumented
raw capture is what identifies this particular event as interior and negative.
Strict checks, pressure floors, closure, forcing and limiter thresholds were
not relaxed. The accepted B4 sharp-contact limitation remains separate.

## Bounded controls

| Control | Result |
| --- | --- |
| Same 64-rank/eight-node replay from $t\simeq2$ | Reproduces the original stage-13 pressure floor |
| Eight ranks/one node, same $48^3$ block grid | Reaches cycle 5800 with a different thermal state; the next original outer step passes |
| Direct cycle 5700→5801 versus an extra restart at cycle 5800 | All retained full-precision interior flow, J/A, B and pressure fields are exactly equal |
| Eight-rank continuation from cycle 5700 to $t=2.6$ | Completes cycle 6103 with zero accumulated LF floors, nonfinite or nonpositive events |
| Initial planned cap/explicit sequence | Stopped when the recreated one-node control passed; later controls use the actual cold checkpoint described below |
| Eight-rank/eight-GPU and 64-rank/64-GPU HIP one-step controls at $t\simeq2$ | Both pass to cycle 4261 at exactly $t=2.0005274584221415$; differences summarized below |
| 64 ranks sharing eight GPUs, intended one-step comparison | HIP initialization GPU hang before evolution; no final state to compare |

The ordinary eight-rank/eight-GPU two-cycle health smoke subsequently passed.
The shared-GPU launch failure is infrastructure evidence, not a CGL result.
No CPU solver control was run; full-precision restart decoding and comparison
were performed independently on CPU.

## First-step differences across MPI layouts

The proper one-GPU-per-rank controls use the same binary and identical
$t\simeq2$ checkpoint. Full-precision comparison of all active cells gives:

| Field | RMS absolute difference | Maximum absolute difference | Maximum's distance from a block face (cells) |
| --- | ---: | ---: | ---: |
| $p_\parallel$ | $2.65\times10^{-15}$ | $1.92\times10^{-13}$ | 10 |
| $p_\perp$ | $8.25\times10^{-16}$ | $6.04\times10^{-14}$ | 6 |
| Stored conservative $J$ | $4.88\times10^{-16}$ | $3.69\times10^{-14}$ | 4 |
| Stored conservative $A$ | $5.37\times10^{-16}$ | $3.72\times10^{-14}$ | 4 |

Density, momentum and cell-centered magnetic components differ by at most
$6.66\times10^{-16}$. These first-step differences are at floating-point
scale, far smaller than the later thermal divergence. Pressure is independently
decoded from the retained J/A variables; the raw J/A comparison does not depend
on that decoding.

The unchanged $48^3$ block grid has 384 of its 768 directed face-neighbor
relations change from MPI to local communication. Rank ownership is reconstructed
from the verified uniform restart block costs and the source's contiguous
load-balancing rule. Distance is the number of active cell layers from a face,
with boundary cells at distance zero. For cells within two layers of a changed
face, pressure RMS differences are $2.72\times10^{-15}$ and
$8.39\times10^{-16}$; for cells at least four layers from those faces, they are
$2.64\times10^{-15}$ and $8.23\times10^{-16}$.
There is no strong interface concentration in this first-step comparison.
The JSON also retains four-layer strata, maximum locations and face-neighbor
rank mappings. Edge/corner-only neighbor relations are not a separate stratum.

The cumulative intermediate hard-bound counter differs by one out of roughly
$1.38\times10^{10}$ counted cell-stage events; floor, nonfinite and nonpositive
counters agree. This counter is not a physical threshold occupancy. A
roundoff-sized initial difference followed by later thermal amplification is
consistent with nonlinear threshold sensitivity, but does not exclude a defect
that becomes important later. No new acceptance tolerance is asserted.

## Thermal divergence precedes the fatal step

The retained 64-rank and eight-rank snapshots both have cycle 5776,
$t\simeq2.5003003035$; their times differ by $8.9\times10^{-16}$.
Their flow and B differences are at float32 quantization scale: at most
$1.2\times10^{-7}$ in any stored value and $5.5\times10^{-11}$ RMS.
Only one density cell differs, with at most 104 differing cells per velocity
component and 52 per B component, out of 14,155,776 cells.

| Absolute pressure difference | Median | 99th percentile | RMS | Maximum |
| --- | ---: | ---: | ---: | ---: |
| $p_\parallel$ | 0.000266 | 0.03341 | 0.008129 | 4.00334 |
| $p_\perp$ | 0.000200 | 0.02549 | 0.006191 | 4.01946 |

Exactly **one cell** differs by more than one pressure unit in each component:
the cell that later fails. There its original pressures are 4.04852/3.93480,
versus 8.05186/7.95426 in the one-node continuation; stored flow and B are
identical at that cell. The thermal outlier therefore exists about 24 cycles
before the failure. Whole-box pressure distributions remain very similar and
would conceal this localized problem. These are descriptive differences,
not invented acceptance tolerances.

## Controls from the actual cold checkpoint

A proper 64-rank replay retained a native cycle-5800 checkpoint at
$t=2.507714514270586$, before the fatal step. Its failing-cell state is
$p_\parallel=0.1572801904$, $p_\perp=0.04455490248$, and $U=0.1231949977$,
matching the instrumented sweep's initial state. Restarting this checkpoint on
64 ranks reproduces the original stage-13/15 pressure floor. An eight-rank
restart reproduces the same stage and pressure-floor count, establishing a
usable negative control. Its rank-local intermediate wall count is 18,709
rather than 2,254 because each rank now owns 16 blocks rather than two.

From this identical state, one-node/eight-GPU controls advance to exactly
$t=2.508023000861065$, one original outer-step duration
$\Delta t=0.0003084865904790016$. Only the time-integration settings change.
Closure, limiter rate and thresholds, forcing, conservative variables and strict
checks remain fixed. All three alternatives finish with zero accumulated LF
floor, nonfinite or nonpositive events.

| Integration | Outer cycles | Final cold-cell $p_\parallel$ | Final cold-cell $p_\perp$ | Pressure RMS difference from explicit ($\parallel$, $\perp$) |
| --- | ---: | ---: | ---: | ---: |
| Explicit LF | 289 | 1.46652 | 1.35381 | — |
| RKL2, outer ratio cap 4 | 25 | 1.44434 | 1.33163 | 0.001079, 0.001032 |
| RKL2, outer ratio cap 14 | 8 | 1.02855 | 0.915847 | 0.001345, 0.001230 |

The cold cell is also the final pressure minimum in each case. Cap 4 differs
there from explicit by 1.5–1.6%; cap 14 differs by 30–32%. The maximum absolute
pressure differences anywhere are 0.1903 for cap 4 and 0.4380 for cap 14.
Cap 4's maximum occurs in a different, warmer cell; cap 14's maximum occurs at
the original cold cell. Whole-box relative $L^2$ pressure differences are only
about $2\times10^{-4}$ and conceal the important local discrepancy.
These are measured comparisons with a smaller-step explicit integration,
not a demonstration of time convergence or an exact reference solution.
Cap 14's successful exit does not qualify it as an accurate guard.

The volume-integrated passive thermal energies differ from explicit by
$4.09\times10^{-8}$ for cap 4 and $1.89\times10^{-7}$ for cap 14, out of
17.17078535. The small integral differences also conceal local disagreement.
The JSON retains kinetic, face-based magnetic and passive thermal integrals;
these compare the full coupled integrations and are not an isolated LF energy
budget or a physical conservation claim.

Source inspection identifies a plausible source of nonlinear sensitivity:
`BuildCGLLFFaceState` adds `LimiterCollisionRate` to the face collision rate.
With the prescribed $10^{10}$ rate, crossing a soft threshold sharply changes
the LF conductivity between stages. The timestep bound deliberately uses only
background collisions, conservatively ignoring this suppression. That bound
does not prove positivity for a nonlinear multistage trajectory. These
controls retain the original conductivity convention; they do not establish
that this switch alone caused the earlier thermal depletion.

## Measured cost of the bounded guards

Timing below subtracts the first progress record from the penultimate one,
excluding initialization, final output and the last clipped step. Physical
duration is the sum of the corresponding printed timesteps, rather than a
difference of rounded absolute times. Measurements use the old $48^3$ block
layout on one node, and are local to this short, depleted-state interval.

| Integration | Timed cycles | Physical time advanced | Measured wall seconds | Node-hours per unit physical time |
| --- | ---: | ---: | ---: | ---: |
| Cap 4 | 24 | 0.00030522112 | 6.91761 | 6.296 |
| Explicit LF | 288 | 0.00030822585 | 56.1153 | 50.572 |
| Cap 14 | 7 | 0.00030011069 | 3.15077 | 2.916 |

All cases retain configured CFL 0.3. Their untruncated outer timesteps span
$[1.07621,1.35400]\times10^{-5}$ for cap 4,
$[8.73298,11.5566]\times10^{-7}$ for explicit, and
$[3.76673,4.74485]\times10^{-5}$ for cap 14. Explicit starts with
$\Delta t=8.968413\times10^{-7}$, confirming that restart initialization
recomputed the timestep after the override. STS safety remains 0.9; the outer
ratio cap and explicit CFL use different multipliers on the LF parabolic limit.

These measurements confirm a substantial cost for cap 4. They are not a
forecast for an entire run with evolving stiffness, different block sizes and
output overhead. The four solver launches plus retained-field analysis used a
131-second one-node allocation, or 0.03639 node-hours. No longer capped run was
launched on the strength of this single-interval result.

## Provenance, cost and remaining work

Raw evidence is under
`/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/physics-benchmark/matched/floor-diagnosis/`.
Each solver case retains its exact command, restart hash, effective input,
executable manifest, log and timing. The compact JSON indexes those files.
`bounded_stability_probes.py` retains the reproducible diagnostic procedure;
`reconstruct_capture.py` independently evaluates the captured recurrence.

The six follow-up solver invocations listed in the JSON used 625.14 launch-wall
seconds on one node/eight GPUs, including a 1.08-second administrative parser
rejection before evolution. This excludes allocation idle time and the earlier
$t\simeq2\to5700$ recreation. It is not the cost of a completed benchmark.
The subsequent proper 64-GPU one-step control used a 27-second eight-node
allocation, or 0.060 node-hours; its solver-launch wall time was 6.63 seconds.
Capturing the actual cold state and repeating the original failure used a
342-second eight-node allocation, or 0.760 node-hours.
Output cadence changes retained a matching $t\simeq2.5$ snapshot without
editing the checkpoint state. The independent restart reader matches native
kinetic/thermal history integrals to summation roundoff and native face-based
magnetic energy exactly; cell-centered magnetic energy is labeled separately.

The amplification of initially roundoff-sized thermal differences across
numerical layouts remains to be localized. Nonlinear threshold sensitivity is plausible, but a spatial,
communication or timestep-control defect has not been excluded. The reproducible
cold state now provides a retained test for candidate integrators. Cap 4 is a
more promising guard than cap 14 in this test, with a material computational
cost and no demonstrated long-run guarantee. Fresh strict matched runs and
adequate time sampling are still required before any active/passive physical
validation claim.
