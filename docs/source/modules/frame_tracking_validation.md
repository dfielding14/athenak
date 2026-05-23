# Frame Tracking Medium-Resolution Validation

This page reports diagnostic validation of the frame tracker on the shipped
cloud-crushing and TRML frame-aware examples. The results diagnose and fix a
restart discontinuity, but they do not satisfy the physical production gate.
Production guidance remains withheld.

## Status

Serial tracking-on versus tracking-off comparisons remain blocking:

| Problem and metric | Resolution and duration | Relative difference | Health events | Result |
| --- | --- | ---: | ---: | --- |
| Cloud density-selected mass | `48 x 16 x 16`, `tlim=0.04` | `2.0914e-2` | `0` | Fail: exceeds `1.0e-2`. |
| Cloud density-selected mass | `96 x 32 x 32`, `tlim=0.04` | `1.4178e-2` | `0` | Fail: improves with resolution but remains above tolerance. |
| TRML temperature-window mass | `16 x 16 x 32`, `tlim=0.25` | `4.0290e-2` | `5` limit events | Fail. |
| TRML temperature-window mass | `32 x 32 x 64`, `tlim=0.25` | `6.0872e-2` | `3` limit events | Fail. |
| TRML passive-tracer mass diagnostic | `16 x 16 x 32`, `tlim=0.25` | `2.4241e-3` | Same controller events | Measured pass for mass only. |
| TRML passive-tracer mass diagnostic | `32 x 32 x 64`, `tlim=0.25` | `4.0503e-3` | Same controller events | Measured pass for mass only. |

Centroid checks pass in these runs, and all inspected histories and snapshots
are finite. The tracer result is a diagnostic of selection sensitivity, not a
replacement acceptance criterion: the documented controller still selects a
temperature interval and remains slew limited in the TRML case.

Download the current on/off diagnostic rows:
[frame_tracking_resolution_sensitivity.csv](../_static/frame_tracking_resolution_sensitivity.csv).

## Boundary And Source Audit

The shipped frame-aware boundaries use the required transforms:

```text
x_lab = x_grid + FrameDisplacement(axis)
v_grid = v_lab - FrameVelocity(axis)
```

The cloud cooling and TRML cooling sources operate on local thermodynamic
state; their cooling rate does not require a position or velocity-frame
transformation. The applied Galilean boost preserves internal energy.

Two boundary defects were identified during this audit and corrected:

1. After a frame displacement or boost, the tracker previously left
   non-periodic physical and user ghost states at the prior frame until the
   normal boundary task ran. A continuous run therefore consumed different
   first-post-update ghost states than a restarted run, which refills
   boundaries during initialization.
2. The TRML user boundary initialized density, momentum, and energy but did
   not populate its passive cold-fraction scalar in ghost cells. Tracer-mass
   diagnostics therefore lacked a complete frame-aware inflow contract.

The tracker now reapplies physical/user boundary conditions after a changed
frame state and converts only ghost-zone primitives, avoiding an unnecessary
full active-mesh conversion. The TRML boundary now fills the passive scalar
with the same evaluated cold fraction used for initialization.

## Restart Diagnosis

The earlier strict TRML restart failure combined a runtime defect with an
analysis-contract defect:

- The stale boundary-state path described above made uninterrupted and
  restarted moving-boundary runs follow different post-checkpoint fluxes.
- Native AthenaK binary snapshots store fluid fields as `float`, so they
  cannot certify a `100 * double epsilon` full-field requirement.
- Restart and MPI continuations on an identical grid were unnecessarily
  interpolated in displaced lab coordinates rather than compared directly in
  their common grid frame.

The comparison script now compares same-grid restart/MPI native snapshots
directly, labels them as binary float32 diagnostics, and prevents a field
tolerance below their representable precision. Strict controller quantities
remain tested with double-precision history output. A high-precision
TRML boundary restart regression additionally compares formatted conserved
output at `100 * machine_epsilon` scaled tolerance.

On the corrected medium TRML case, 25, 50, and 75 percent restart splits have
passing controller-continuity rows and passing binary field diagnostics. The
run itself still contains controller limit events, so this is evidence for the
restart fix, not a passed production configuration.

Download the corrected restart diagnostic rows:
[frame_tracking_restart_diagnostic.csv](../_static/frame_tracking_restart_diagnostic.csv).

## Commands

Cloud medium enabled run for the corrected implementation:

```bash
./build_cloud_crushing/src/athena \
  -i inputs/hydro/cloud_crushing_snr.athinput \
  -d run_cloud_medium_on \
  mesh/nx1=96 mesh/nx2=32 mesh/nx3=32 \
  meshblock/nx1=32 meshblock/nx2=16 meshblock/nx3=16 \
  time/tlim=0.04 time/nlim=-1 output1/dt=0.004 output2/dt=0.04
```

The disabled reference uses the same command plus
`frame_tracking/enabled=false`.

TRML medium serial run:

```bash
./build_trml_frame_tracking/src/athena \
  -i inputs/hydro/TRML/TRML_frame_tracking.athinput \
  -d run_trml_medium_on \
  mesh/nx1=32 mesh/nx2=32 mesh/nx3=64 \
  meshblock/nx1=16 meshblock/nx2=16 meshblock/nx3=32 \
  time/tlim=0.25 time/nlim=-1 \
  output1/dt=0.025 output2/dt=0.25 output3/dt=0.0625
```

Temperature-window and passive-tracer diagnostics use the same transformed
snapshots:

```bash
python scripts/compare_frame_tracking_validation.py \
  --reference-dir run_trml_medium_off \
  --candidate-dir run_trml_medium_on \
  --output diagnostics.csv \
  --problem TRML --resolution 32x32x64 \
  --tracking-mode enabled_vs_disabled_temperature \
  --comparison-reference serial_tracking_disabled \
  --comparison-kind physical --axis x3 --selection temperature \
  --target-min 0.015 --target-max 0.08

python scripts/compare_frame_tracking_validation.py \
  --reference-dir run_trml_medium_off \
  --candidate-dir run_trml_medium_on \
  --output diagnostics.csv --append \
  --problem TRML --resolution 32x32x64 \
  --tracking-mode enabled_vs_disabled_tracer \
  --comparison-reference serial_tracking_disabled \
  --comparison-kind physical --axis x3 --selection scalar \
  --target-min 0 --target-max 1
```

## Required Next Work

1. Define and justify a production selected-material observable for TRML.
   The hard instantaneous temperature window fails while conserved tracer mass
   passes, so changing the gate requires a documented scientific decision.
2. Develop a TRML conservative controller recipe that eliminates limit events
   without worsening the accepted physical observable.
3. Add a high-precision full conserved-state path for medium production
   restart certification; native binary snapshots are diagnostic only.
4. Only after serial physical gates pass, run cloud/TRML restart, four-rank
   MPI, and AMR production comparisons and publish those results.

Until these steps pass their stated criteria, the feature is supported for
wiring tests and controlled method development only.
