# Q011 Resource Scaling Audit and Preproduction Plan

Date: 2026-06-06

This is a resource-planning record only. It does not alter a historical
preregistration, authorize execution, inspect qualifying output, or close a
Section 5.4 claim.

## Audit Result

The current qualifying matrix is three grid variants by eight paired seeds:
24 baseline attempts, plus an AMR restart-continuation carrier, failed-attempt
reserve, preproduction work, and output/publication overhead.

The active Frontier cap is 10,000 node-hours. The validated ledger head is
sequence 122 with 1.3019444444444446 node-hours consumed, zero reserved, and
9998.698055555555 node-hours currently unreserved. This cap is project-wide,
not Q011-only; a protected reserve for the nonlinear Bell campaign and other
required methods work remains unbound.

The production geometry is:

| Variant | Grid or bound | Cells | Projected cycles |
|---|---:|---:|---:|
| coarse uniform | 4000 x 260, dx=12 | 1,040,000 | 17,480 |
| three-level AMR | root 4000 x 260, finest dx=3 | 1,040,000 to 16,640,000 active-cell equivalent | about 69,920 after finest-level activation |
| fine uniform | 16000 x 1040, dx=3 | 16,640,000 | 69,920 |

The selected `p0=1.0` pressure pilot used one GPU, five MeshBlocks, 2,000
active cells, 874 cycles, and 142,536 terminal particles. It recorded
54.266364554 driver seconds, 146.669497363 output-publication seconds, and
232 scheduler seconds. This is a valid engineering point, but it is severely
underfilled and cannot determine production multi-node throughput or AMR
occupancy.

## Work Equations

With startup-particle removal at `t=45`, the retained particle-time integral is

```text
I(T) = 0.5 * 45^2 + 0.5 * (T - 45)^2, T >= 45.
```

Thus `I(60)=1125` and `I(1200)=668025`. Using the production-to-pilot
transverse-width ratio of 13 and the four-times-smaller fine timestep:

```text
R_particle(coarse)   = I(1200)/I(60) * 13     = 7719.4
R_particle(amr/fine) = I(1200)/I(60) * 13 * 4 = 30877.6
```

The zone-work ratios are:

```text
R_zone(coarse) = (1,040,000 / 2,000) * (17,480 / 874) = 10,400
R_zone(fine)   = (16,640,000 / 2,000) * (69,920 / 874) = 665,600
R_zone(amr)    = 665,600 * f_amr
```

where `f_amr` is the unknown time-weighted active-cell fraction relative to the
full fine domain.

If the pilot's mixed work rates were unchanged and all eight GPUs on a Frontier
node were used ideally, the driver-only brackets would be:

| Variant | Driver node-hours per run |
|---|---:|
| coarse | 14.55 to 19.60 |
| AMR | particle floor 58.18; 78.38 at root-only active cells; 1254.16 at full-fine occupancy |
| fine | 58.18 to 1254.16 |

The interval is intentionally wide. The pressure pilot supplies only one mixed
field-plus-particle timing equation, so it cannot identify separate work
coefficients. It also does not measure parallel output bandwidth.

For the current eight paired triads, the constant-rate driver-only bracket is
1047.26 to 20223.27 node-hours. Adding only a 25 percent margin raises the
upper value to 25279.08 node-hours. A scenario-only likely envelope, assuming
4-8 times better large-grid GPU throughput, 0.65-0.8 multi-node efficiency,
AMR occupancy of 0.2-0.4, and a 1.35-1.75 total uplift, is about 3098-9768
node-hours. That scenario is not an authorization estimate; it shows why the
current matrix might fit only under favorable scaling and would leave no
robust Bell or failure reserve.

## Recommended Matrix

The recommended minimum is three complete paired triads:

```text
variants: coarse_uniform_dx12, three_level_amr_root_dx12_finest_dx3,
          fine_uniform_dx3
seeds:    23050101, 23050102, 23050103
attempts: 9
```

Three paired seeds preserve direct grid comparisons and provide the minimum
ensemble that can expose seed variability and support a sample variance. Two
paired seeds are not adequate for a publication-grade stochastic-method
comparison. An asymmetric matrix is not recommended because the current
AMR-versus-fine release criteria require a matching fine run for each compared
AMR seed.

Additional complete paired triads may be added, in existing seed order, only
from resource telemetry gathered before qualifying-output inspection. The
historical 24-run preregistration remains unchanged; adopting this recommendation
requires a separate versioned qualifying-preregistration successor.

## Scaling Pilots

All pilots use `problem/ps_p0=1.0`, engineering seeds disjoint from the
qualifying seeds, registered control-plane execution, and a hard total pilot
cap of 500 consumed node-hours. Physical fields may not be inspected. Only
scheduler accounting, runtime telemetry, active-cell and MeshBlock histories,
particle counts and updates, memory high-water marks, output timers, and
artifact byte counts may be used.

### Phase 1: Full-Geometry Short-Step Scaling

Use exact production geometry and physics with a separately tested
engineering-only output-suppressed deck and `nlim=512`.

| Variant | Node ladder |
|---|---|
| coarse | 4, 8, 16 |
| AMR | 8, 16, 32 |
| fine | 32, 64, 128 |

Run each point once, then repeat the selected node count three times. Select
the smallest node count whose median node-hours per completed cycle is within
10 percent of the observed minimum and that passes memory and walltime gates.

### Phase 2: Reduced-Transverse Full-Time Calibration

Run all three variants to `t=1200` with full production output cadence but a
transverse extent of 240 instead of 3120:

```text
coarse/AMR root: 4000 x 20
fine:            16000 x 80
```

This preserves the production x extent, terminal time, injection history, and
cell sizes while reducing transverse work by a factor of 13. It measures late
particle loading, AMR occupancy history, output/checkpoint cost, and stability
without creating a qualifying dataset.

### Phase 3: Held-Out Validation

Run one additional full-geometry short-step case per variant at the selected
node count. Fit a nonnegative additive driver model from active-cell cycles and
particle updates, and model I/O separately from bytes and output timers.

For authorization, use the larger of the model projection and the direct
reduced-transverse extrapolation, multiplied by:

```text
max(1.25, 1 + maximum held-out absolute fractional residual).
```

Reject the model if held-out driver-plus-output time differs by more than
20 percent.

## Decision Gates

Production remains blocked until all of the following pass:

1. Exact source/evidence bindings and no qualifying-output inspection.
2. Every variant completes full-geometry scaling without OOM or runtime failure.
3. Reduced-transverse runs reach exact terminal time; production-scaled peak
   device memory is at most 70 percent; no projected output event consumes more
   than 20 percent of requested walltime.
4. The held-out cost model is within 20 percent.
5. The project-wide equation passes:

```text
C_consumed + C_reserved + C_preproduction + C_baseline + C_restart
+ max(C_fine, 0.15 * C_baseline) + R_nonQ011 <= 10000.
```

Here `C_fine` is the conservative cost of one selected fine-uniform baseline
attempt, so the failed-attempt replacement reserve is at least one expensive
run and at least 15 percent of the selected baseline matrix. `R_nonQ011` is a
separately bound protected reserve for the nonlinear Bell campaign and other
required non-Q011 methods work.

6. At least three complete paired triads fit. If they do not, do not fall back
   to two; redesign the campaign or obtain more resources.

The exact machine-readable record is
`q011_resource_scaling_preproduction_plan_2026-06-06.json`.
