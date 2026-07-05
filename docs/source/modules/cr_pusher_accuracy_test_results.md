# CR Pusher Accuracy Test Results

Run location:

```text
/lustre/orion/ast207/proj-shared/dfielding/AMR/particles/testing/cr_pusher_accuracy_20260704_153322
/lustre/orion/ast207/proj-shared/dfielding/AMR/particles/testing/cr_pusher_accuracy_latest
```

Source binding:

```text
branch: feature/CR_tracers_followup_architecture
commit: e7dd583f485d97cfd662ad6ca129b55e6b538230
executable: build-frontier/src/athena
executable sha256: 80696df02b07a07b961f248793fccfe21e9fb7b7c819a94e137bf980fbec7e93
```

The run was intentionally made from an uncommitted test harness state. The
dirty status recorded in `source.status` contains only the new accuracy-test
plan, three input files, the analyzer, and the runner script.

## Matrix

The executed sweep covers the first three tests in
`docs/source/modules/cr_pusher_accuracy_test_plan.md`:

```text
uniform_b:
  global gyro_fraction=0.05
  per-particle gyro_fraction=0.05, 0.10, 0.20, 0.30

smooth_divb0:
  strict reference per-particle gyro_fraction=0.01
  per-particle gyro_fraction=0.02, 0.05, 0.10, 0.20, 0.30
  global gyro_fraction=0.05

boundary_crossing:
  single-block reference
  multi-block global gyro_fraction=0.05
  multi-block per-particle gyro_fraction=0.05, 0.10
```

The harness wrote 16 tracked-particle files. The analyzer's track parser reads
427 frames total and finds zero nonfinite particle records.

## Result

The full analyzer pass succeeds:

```text
Overall hard-check status: PASS
hard_errors: []
rows: 56
max_dE_E: 1.83980563e-07
smooth_error_increases_with_gyro_fraction: true
boundary_max_step_dx: 1.47872981e-01
```

Key accuracy results:

```text
uniform_b global g=0.05 species0 rms_dv:       1.38117828e-04
uniform_b per-particle g=0.05 species0 rms_dv: 1.59660898e-04
uniform_b per-particle g=0.10 species0 rms_dv: 5.51838005e-04
uniform_b per-particle g=0.20 species0 rms_dv: 1.24012584e-03
uniform_b per-particle g=0.30 species0 rms_dv: 4.92266216e-03

smooth_divb0 global g=0.05 species0 rms_dv:       5.37495234e-05
smooth_divb0 per-particle g=0.05 species0 rms_dv: 5.67909386e-05

boundary_crossing tracked-particle dx/dv vs single-block reference: 0
```

Tiny-test wall times are written in each case's `run.time`. These runs are too
small to benchmark production throughput, but they do show the expected
direction from removing global over-subcycling:

```text
uniform_b global g=0.05:       2.12 s
uniform_b per-particle g=0.05: 0.88 s
uniform_b per-particle g=0.10: 0.79 s
smooth_divb0 global g=0.05:    0.59 s
smooth_divb0 per-particle g=0.05: 0.77 s
```

The smooth-field magnetic moment changes are not treated as hard failures in
this test, because the field is intentionally nonuniform. The convergence
check is the trajectory error relative to the strict `gyro_fraction=0.01`
reference, and that error increases monotonically with relaxed gyro fraction.

## Conclusion

`subcycle_per_particle_gyro=true` with `gyro_fraction=0.05` passes the clean
pusher tests and agrees closely with the conservative global-subcycle reference
while avoiding the production over-subcycling cost. `gyro_fraction=0.10` also
looks plausible in these clean tests, but it should remain a candidate setting
until it is compared on production-like frozen-field restarts. `0.20` and
`0.30` show visibly larger phase/trajectory error and should be treated as
experimental.

One parser issue was found and fixed in the new analyzer: the binary `trk`
reader must skip exactly the text separator line after each frame header. An
older whitespace-skipping parser can accidentally consume valid binary payload
bytes and report false nonfinite records.
