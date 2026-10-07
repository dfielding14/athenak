# Passive LF pressure-floor forensics, 2026-10-07

**Numerical behavior is concerning; the root cause and physical validation
remain inconclusive.** The original passive member failed strict checking
because an intermediate LF update produced a finite negative perpendicular
pressure. A replay reproduced it. A different MPI layout followed a different
thermal trajectory and did not fail through $t=2.6$; that is not a qualified
fix. No capped-STS or explicit-LF mitigation was tested successfully or accepted.

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
| Cap-4, cap-14 and explicit-LF variants | **Not executed:** the required failing one-node control did not reproduce |
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
Output cadence changes retained a matching $t\simeq2.5$ snapshot without
editing the checkpoint state. The independent restart reader matches native
kinetic/thermal history integrals to summation roundoff and native face-based
magnetic energy exactly; cell-centered magnetic energy is labeled separately.

The amplification of initially roundoff-sized thermal differences across
numerical layouts remains to be localized. Nonlinear threshold sensitivity is plausible, but a spatial,
communication or timestep-control defect has not been excluded. A reproducible
failing synchronized state is needed to qualify a numerical mitigation. Fresh
strict matched runs and adequate time sampling are still required before any
active/passive physical validation claim.
