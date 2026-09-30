WO1 scratch performance and byte-comparison harness
=================================================

Status: the definitive post-G baseline is complete at runs/post-g. Its source
commit is 9655659fa7958724ad1c9f127b9989e01675443d. All 24 configurations passed
three-repeat byte comparison, including full restart files. Binary hashes and
provenance are recorded there. Do not restage inputs after this freeze.
Earlier harness-check-* runs only validated harness repeatability.

Python: /Users/dbf75/.uv/envs/interactive/.venv/bin/python3
Script: /tmp/cgl-wo1-perf/harness.py

1. Run `harness.py stage` once with the final post-G inputs. The input manifest
   records source and staged SHA-256 hashes. Do not stage again after baseline.
2. Run:
   harness.py run --label post-g --serial /path/to/serial/athena \
     --mpi /path/to/mpi/athena --repeats 3 --mpi-ranks 1 4
3. After each candidate, use the identical invocation with a fresh label and
   `--compare post-g`. Both executables are copied and hashed at run start.
   Optional `--cases NAME ...` restricts the matrix. Existing labels are never
   overwritten. Runs use the same fixed absolute input path, without overrides.

Required cases: lf1d, lf2d, lf3d, smr, pure_cgl.
Risk cases: shear, smr_varied, smr_outflow, outflow, density, low_b, velocity.
MPI 1 and 4 ranks: lf2d, lf3d, smr, shear, smr_varied, smr_outflow. All cases also run serial.

The required 1D/2D/SMR cases run for about one decay time; 3D runs ten cycles.
The pure-CGL case uses the existing wave regression input. Risk cases exercise
shearing boundaries, a two-region stepped SMR interface with primitive
prolongation/nonuniform density/nonzero velocity/a passive scalar, its variant
with refinement touching an outflow boundary, uniform-grid outflow,
a density jump of 200, sub-floor B, and finite-velocity field waves. The shear
case uses the existing cgl_lf_sbox generator with amp=1e-4, bx0=1, and an STS
ratio of 20 so strict admissibility passes on the control binary.

Every case writes conserved+B and primitive+B binary snapshots initially,
every ten cycles, and finally; high-precision history every cycle; and Real
precision restart snapshots. Both varied SMR cases have binary ghost zones. Files
are compared byte-for-byte with matching inventories across repeats and
against the baseline, including full restart metadata/payload. Binary field
payloads are float32 in AthenaK; full restart equality additionally checks
full-precision conserved and face-field values including ghost zones. If a
restart metadata field becomes nondeterministic, diagnose it explicitly;
do not silently discard differing bytes. Nothing is normalized currently.

Results are in runs/LABEL/results.json, with per-run stdout.log and output
files retained. SHA-256 hashes, cycle/stage counts, all timing samples, medians,
and candidate/baseline ratios are recorded. The default is three repeats.
A byte mismatch produces a nonzero exit. Runtime or geometry/count failures
also abort with their run directory retained.

Timing definitions:
- solver_s_per_cycle: AthenaK's `cpu time used` divided by completed cycles.
  Despite the historical label, this is the Kokkos driver wall-clock timer.
- wall_s_per_cycle: external process wall time divided by completed cycles;
  includes startup, MPI launch, shutdown, and outputs.
- sts_profile_s_per_stage: sum of existing exclusive LF/STS compute buckets,
  using rank-mean seconds for MPI, divided by actual global RKL stage count.
  heat_flux_total includes its precompute/directional/work subregions, so these
  nested buckets are not added again. The total stage count is cross-checked
  against lf_nstage / actual physical cell count, not treated as cell visits.
- shared_transport_s_per_stage: existing shared transport bucket sum divided
  by the same stage count. These buckets include RK and initialization calls;
  this is NOT a measurement of STS-only communication time.

Whole-STS wall time is not separately available from existing timers. The
profiled LF/STS compute time per stage is a partial measured cost, not a pure
whole-STS wall time. Profiling fences and mandatory outputs affect absolute
costs; only same-input, same-output, same-rank comparisons are meaningful.

After each completed candidate, `report.py LABEL` writes comparison.md with
the five required serial cases and 2D/3D/SMR MPI4 cases, before/after medians
and ratios. Rejected candidates and their exact failure lists remain intact.
For incremental performance, `report.py p3 --reference p2` writes the P2-to-P3
table and adds cumulative post-G ratios; byte comparison still uses post-G.
