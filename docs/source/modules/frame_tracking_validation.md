# Frame Tracking Medium-Resolution Validation

This page reports the first medium-resolution scientific validation runs for
the frame tracker. These runs are intended to test whether tracked and
untracked calculations agree after transformation to lab coordinates and
whether restart continuation preserves the evolving solution. They do not
currently satisfy the production-candidate acceptance gate.

## Status

The required uniform-grid serial tracking-on versus tracking-off comparison
failed for both shipped frame-aware examples:

| Problem | Medium run | Selected-mass relative difference | Limit events | Gate result |
| --- | --- | ---: | ---: | --- |
| Cloud crushing | `96 x 32 x 32`, `tlim=0.04` | `1.4178e-2` | `0` | Fail: exceeds `1.0e-2`. |
| TRML | `32 x 32 x 64`, `tlim=0.25` | `6.0872e-2` | `3` | Fail: exceeds `1.0e-2` and is slew limited. |

Centroid comparisons pass the one-cell acceptance criterion in both cases, and
all inspected snapshots and history data are finite. The failed selected-mass
criteria are sufficient to withhold production guidance. Longer cloud
restart/MPI/AMR variants are therefore not presented as completed validation;
they should be run after resolving or explicitly accepting the physical
discrepancy.

## Executed Comparisons

The calculations use a Release CPU build of the feature branch at executable
revision `505f2df4`. The example inputs include an explicit canonical
`enabled = true` entry so the disabled reference can be selected with
`frame_tracking/enabled=false` on the command line. Enabled behavior is
unchanged from the prior default-enabled input.

| Comparison | Variants | Outcome |
| --- | --- | --- |
| Cloud physical comparison | Serial uniform grid, tracking enabled versus disabled | Fails selected-mass tolerance; no miss or limit events. |
| TRML physical comparison | Serial uniform grid, tracking enabled versus disabled | Fails selected-mass tolerance; three limit events. |
| TRML restart split | Enabled continuous versus restarts at approximately 25, 50, and 75 percent of final time | Strict field/controller continuation criterion not met; see raw CSV. |
| TRML MPI comparison | Enabled serial versus 4-rank MPI, uniform grid | Field and controller differences satisfy `1.0e-10`; the same three limit events retain a failed health row. |

The restart comparisons apply the specified `100 * machine_epsilon` threshold
for the double-precision build. Controller differences are small, but some
fluid-field norms and final frame velocities exceed that deliberately strict
threshold. This is a blocking result, not a tolerance adjustment request.

## Commands

Cloud medium-resolution serial runs:

```bash
./build_cloud_crushing/src/athena \
  -i inputs/hydro/cloud_crushing_snr.athinput \
  -d run_cloud_medium_on \
  mesh/nx1=96 mesh/nx2=32 mesh/nx3=32 \
  meshblock/nx1=32 meshblock/nx2=16 meshblock/nx3=16 \
  time/tlim=0.04 time/nlim=-1 output1/dt=0.004 output2/dt=0.01

./build_cloud_crushing/src/athena \
  -i inputs/hydro/cloud_crushing_snr.athinput \
  -d run_cloud_medium_off \
  mesh/nx1=96 mesh/nx2=32 mesh/nx3=32 \
  meshblock/nx1=32 meshblock/nx2=16 meshblock/nx3=16 \
  time/tlim=0.04 time/nlim=-1 output1/dt=0.004 output2/dt=0.01 \
  frame_tracking/enabled=false
```

TRML medium-resolution serial runs:

```bash
./build_trml_frame_tracking/src/athena \
  -i inputs/hydro/TRML/TRML_frame_tracking.athinput \
  -d run_trml_medium_on \
  mesh/nx1=32 mesh/nx2=32 mesh/nx3=64 \
  meshblock/nx1=16 meshblock/nx2=16 meshblock/nx3=32 \
  time/tlim=0.25 time/nlim=-1 \
  output1/dt=0.025 output2/dt=0.0625 output3/dt=0.0625

./build_trml_frame_tracking/src/athena \
  -i inputs/hydro/TRML/TRML_frame_tracking.athinput \
  -d run_trml_medium_off \
  mesh/nx1=32 mesh/nx2=32 mesh/nx3=64 \
  meshblock/nx1=16 meshblock/nx2=16 meshblock/nx3=32 \
  time/tlim=0.25 time/nlim=-1 \
  output1/dt=0.025 output2/dt=0.0625 output3/dt=0.0625 \
  frame_tracking/enabled=false
```

The reusable comparator reads native binary snapshots and
`.frame_tracker.hst`, transforms enabled snapshots back to lab coordinates,
and emits a common CSV schema:

```bash
python scripts/compare_frame_tracking_validation.py \
  --reference-dir run_cloud_medium_off \
  --candidate-dir run_cloud_medium_on \
  --output docs/source/_static/frame_tracking_validation_summary.csv \
  --problem cloud_crushing --resolution 96x32x32 \
  --tracking-mode enabled_vs_disabled \
  --comparison-reference serial_tracking_disabled \
  --comparison-kind physical --axis x1 --selection density \
  --target-min 5 --target-max 200
```

Download the measured rows:
[frame_tracking_validation_summary.csv](../_static/frame_tracking_validation_summary.csv).

## Required Next Work

1. Diagnose the selected-mass discrepancy in the transformed tracked versus
   untracked solutions, beginning with selection sensitivity and frame-aware
   boundary/source consistency.
2. Tune or justify the TRML conservative slew limits so the documented
   configuration does not silently saturate during its validation interval.
3. Investigate the strict restart-continuity mismatch before relying on long
   split production runs.
4. After those blockers are resolved, rerun the full serial, restart,
   four-rank MPI, and AMR matrix and publish the updated CSV and field-norm
   comparisons. The current uniform-grid TRML MPI comparison is a useful
   numerical diagnostic, not a passed production validation.

Until these steps pass their stated criteria, the feature is supported for
wiring tests and controlled method development only.
