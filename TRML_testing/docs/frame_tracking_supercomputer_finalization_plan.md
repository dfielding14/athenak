# TRML Frame-Tracking Finalization Plan for Supercomputer Campaign

## 1. Current Implementation Status (from latest `TRML_testing/`)

### 1.1 What is working
- Modern-only parser and validation are active.
- New diagnostics (`z_centroid_ft`, `z_midpoint_ft`, `z_ctrl_ft`, `v_p`, `v_d`, `v_i`) are present in runtime logs.
- Smoke and screening runs complete without NaN/Inf or command-cap chatter:
  - `TRML_testing/nextstep_impl_runs_20260221_172836/smoke_metrics.json`
  - `TRML_testing/nextstep_impl_runs_20260221_172836/screen_metrics.json`

### 1.2 What remains unresolved
- Long-horizon robustness is not yet proven for extended-z, high-`xi` setups.
- In `xi=100` (`TRML_testing/n1_zext_xi_suite_20260221_221837/xi_100/run.log`):
  - run reached `t=100`,
  - late-time `FrameTracking skipped` appears repeatedly (`165` skips),
  - final `mean_z_ft` is large positive (`~0.476`), indicating interface escape from the intended tracked regime.
- `xi=300` campaign was manually stopped at `t~46.3`:
  - `TRML_testing/n1_zext_xi_suite_20260221_221837/suite_run.log`

### 1.3 Interpretation
- Controller numerics are stable; failures are dominated by **observable/band robustness** over long times, not immediate instability.
- Shear bounds are useful diagnostics but **not a hard failure criterion** for finalization.

---

## 2. Finalization Objectives

The supercomputer campaign should close four risks:
1. **Tracking continuity**: avoid prolonged loss of tracked gas in the chosen thermal band.
2. **Centering quality**: suppress slow bias drift in `mean_z_ft` over long windows.
3. **Operational robustness**: maintain stable behavior across rank counts and restart boundaries.
4. **Scalability**: quantify frame-tracking overhead in production-scale MPI runs.

---

## 3. Required Test Matrix

## Gate A: Build + Schema (must pass first)
| Test ID | Purpose | Run/Command | Pass Criteria |
|---|---|---|---|
| A1 | Hydro build | `cmake --build <build-hydro-mpi> -j` | Clean build |
| A2 | MHD build | `cmake --build <build-mhd-mpi> -j` | Clean build |
| A3 | Schema-positive | modern-only input boots | Starts and advances |
| A4 | Schema-negative (legacy keys) | inject `boost_every` in `<problem>` | Fails fast with migration message |
| A5 | Schema-negative (new fields) | bad `ft_position_signal`, bad blend, bad `ft_int_leak_tau` | Correct fatal validation message |

Reference logs: `TRML_testing/nextstep_impl_schema_20260221_172736/*.log`.

## Gate B: Single-node regression (fast)
Use `128^3`, `64^3` meshblocks, `8` ranks, `tlim=12`.

| Test ID | Parameters | Purpose | Pass Criteria |
|---|---|---|---|
| B1 | `R0` (centroid, no integral) | Regression baseline | No NaN/Inf; diagnostics emitted |
| B2 | `N1` candidate | Current best centering candidate | Better or equal centering vs B1 |
| B3 | `N2` candidate | Bound-strong alternative | Stable + no pathological skips |
| B4 | `N3` candidate | Conservation-lean alternative | Stable + acceptable centering |

Ranking signal (late window `t>=6`):
- Primary: lowest `abs(mean(mean_z_ft))`.
- Secondary: lowest `std(mean_z_ft)`.
- Tertiary: lowest skip fraction and momentum drift.

---

## Gate C: Extended-z robustness (primary unresolved risk)
Use:
- domain `x3 in [-1.0, 0.5]`,
- resolution `128 x 128 x 192`,
- meshblocks `64 x 64 x 96`,
- `8` ranks minimum for functional parity,
- `tlim=100`.

Run exactly:
1. `xi=100`
2. `xi=300`
3. `xi=1000`

### Required outputs per run
- `run.log`
- `TRML.user.hst`
- `TRML.hydro.hst`
- binary snapshots/slices
- plot set from `TRML_testing/plot_run_diagnostics.py`

### Gate C pass criteria (evaluate on `t>=20`)
- `mean_z_ft` remains bounded with no monotonic runaway trend.
- Skip behavior:
  - no persistent skip streaks after warmup (`t>5`),
  - skip fraction remains low enough that controller remains active most of the time.
- Command behavior:
  - no prolonged cap-locking (`vel_z_shift` at hard cap for long contiguous intervals).
- Conservation diagnostics remain physically acceptable for chosen forcing/cooling regime.

Note: shear extrema (`zmin_shear`, `zmax_shear`) are tracked and reported, but do not hard-fail the run.

---

## Gate D: Restart continuity
For the Gate C winner:
1. Continuous run `0 -> 100`.
2. Restarted run `0 -> 50` + restart `50 -> 100`.

Pass criteria:
- no persistent discontinuity in `z_filt`, `v_i`, `vel_z_shift` beyond one actuation interval,
- late-time statistics agree within normal turbulent variance.

---

## Gate E: MPI scale-out validation on supercomputer

## E1. Strong scaling
- Fixed problem: at least `256^3` (or larger if queue policy allows).
- Rank ladder: `64, 128, 256, 512` (and higher if available).
- Record:
  - total walltime,
  - time spent in frame-tracking cadence section (if instrumented),
  - effective zone-cycles/s.

## E2. Weak scaling
- Keep zones-per-rank approximately fixed (same meshblock shape policy).
- Increase ranks/nodes to production target.
- Confirm control metrics (skip fraction, centering, command behavior) stay statistically comparable.

## E3. Reduction-overhead check
- Compare with `use_frame_tracking=false` at identical decomposition.
- Quantify cost increase from frame-tracking reductions and control kernel.

---

## 4. Standard Analysis Workflow

For each run directory:

```bash
mkdir -p history
ln -sfn ../TRML.user.hst history/TRML.user.hst
ln -sfn ../TRML.hydro.hst history/TRML.hydro.hst
python3 TRML_testing/plot_run_diagnostics.py --run-dir . --log-file run.log --out-dir plots
```

Minimum plots to archive:
- `frame_tracking_timeseries.png`
- `interface_z_ranges.png`
- `hydro_history.png`
- `temperature_slices_3panel.png`
- `velx_slices_3panel.png`
- `vely_slices_3panel.png`
- `velz_slices_3panel.png`
- `density_projection_y_3panel.png`

---

## 5. Required Artifacts for Final Sign-off

For every candidate run:
1. input file (`TRML.athinput`),
2. raw logs/histories,
3. plot directory,
4. machine metadata (node count, ranks/node, compiler, build flags),
5. concise metrics summary JSON/CSV.

For final selected configuration:
1. full `t=100` continuous run package,
2. matching restart package,
3. strong/weak scaling summary,
4. short decision memo documenting why the selected parameterization is production-ready.

