# Matched active/passive CGL–Landau-fluid turbulence benchmark

This experiment compares active CGL with a passive CGL control on the same
`192×192×384` periodic grid, using PPM (`ppm4`) and the same retained forcing
parameters. It asks whether anisotropic-pressure feedback changes magnetic
strength fluctuations, pressure-anisotropy occupancy, velocity gradients and
pressure compensation in this finite-time realization.

**Physical validation status: inconclusive; the lower-power passive run also failed.**
The canonical injection is now `dedt=0.16` per unit volume, giving nominal
total power `0.32` in the volume-two box. Start both members from the prescribed
initial state using the existing committed solver. First run passive to `t=3`;
only if it completes healthily continue that same lineage to `t=14`, then run
the fresh active member. Adaptive LF stepping is not part of this experiment.
Lower forcing is a changed experiment, not an established cure for the LF
failure or a claim of physical acceptance.

The lower-power passive member passed `t=3` in job 5635675, then job 5636171
failed at cycle 10316, `t=4.398059` (located by replay): pre-LF RKL2 stage 7/17 triggered one pressure-floor
event. The two allocations used 2.630833 node-hours. The automatic continuation
stopped; active and the paired analysis have not started. Evidence is retained
under `matched/power032`, with the failed log in `passive-continue-to14/run.log`.
The original executable was replayed to cycle 10300 to retain a checkpoint
immediately before the failure. A small optional LF-only chunk cap
(`time/cgl_lf_max_chunk_ratio`, numerical revision `4af7c9a55`) preserves
the outer CFL, forcing, and logical collision cadence. Five focused GPU
regressions passed: unchanged default, collision/limiter cadence, smooth LF
refinement and conservation, restart consistency, and configuration rejection.
Qualification job 5637361 then stopped after 55 seconds because its launcher
tried to override a key absent from the old restart. The corrected launcher
uses a minimal input overlay and verifies merged parameters before evolving.

Corrected qualification 5637416 completed in 625 seconds on one node. Its
disabled-cap control reproduced the original stage failure. Chunk ratios
4, 2, and 1 all passed the same 25 outer steps with zero floor/nonfinite/
nonpositive counters and identical mechanical fields at snapshot precision.
Relative to ratio 1, ratio 4 differed by 0.029% RMS in parallel pressure and
0.027% in perpendicular pressure, but worst local differences remained 14.3%
and 13.9%. Ratio 1 is not a converged reference: local accuracy remains
unresolved despite second-order refinement in the separate smooth test.
The ratio-4 probe cost 124.58 seconds for 0.00838 simulation time, equivalent
to 4.13 node-hours per time unit locally, including startup and output. This
short stiff interval does not measure the cost of a complete run.

Job 5637477 completed a bounded ratio-4 passive continuation from `t=4.40107`
to `t=4.9760169`, stopping cleanly at its 1h50 wall-clock limit before `t=5`.
The final snapshot passed its health check and all recorded failure counters
were zero. It used one node, eight large blocks, `dedt=0.16`, and CFL `0.3`.
This is a robustness diagnostic, not a fresh matched reference or an automatic
launch of the active member. Long-run stability is unproven. Its measured rate
was 3.730 million zone-cycles/s/node. The subsequent
[performance audit](cgl_lf_performance_audit_20261008.md) found a 4.65× speedup
by disabling only detailed LF q diagnostics, preserving strict checks and saved
endpoint fields in a short comparison. Following that audit, the canonical input
now disables those detailed q sums while retaining safe arithmetic, weighted
fluxes, strict admissibility, and pressure/forcing-work recording. Historical
inputs remain unchanged. A resumed run with a changed diagnostic mode requires
an explicitly recorded override; do not relabel the old segment's settings.
The full-build compiler comparison and matched timing have now completed
(jobs 5638165/5638232). The current corrected build reaches 203.2 million
zone-cycles/s/node on the exact archived uniform scaling input versus 103.2
million for the archived executable, partly because the corrected timestep
bound reduces LF stages from 42 to 14 per outer step. The turbulent checkpoint
still uses about 146 stages and measures 15.44 million with diagnostics off.
Removing correctness flags gives only 4.9% on the equal-work uniform test and
reproduces an incorrect minimum reduction, so retain those flags. The fresh
corrected build matches all retained endpoint fields at snapshot precision.
The user explicitly requires detailed LF diagnostics off; do not spend further
work on full-mode diagnostic optimization or enable it for this benchmark.
Chunk tuning and speculative kernel experiments have also been stopped at the
user's request. No larger chunk setting or experimental kernel patch was
adopted. The current task is a source and commit comparison against the exact
pre-WO1 baseline, separating WO1, WO2, and later changes. The matched
scientific run is still incomplete, and no simulation is currently running.
See `matched/power032/lf-subcycling/README.md`, `comparison.json`,
`continuation-submission.json`, and the retained qualification records.
The commands below describe the original
lower-power experiment and must not be read as evidence that it completed.

The first 192-grid pair used 48-cubed blocks and eight nodes per member.
It and the subsequent higher-power preflights used `dedt=0.32`, nominal total
power `0.64`; their inputs, measurements and reports remain historical evidence.
Its passive member failed strict admissibility at `t=2.5077145`, cycle 5800;
the active member was stopped on request after `t=6.16`. Preserve both as
failed/interrupted evidence, not a completed stationary comparison.
An unchanged-settings replay reproduced an interior perpendicular pressure
of `-0.00875815` at RKL2 pre-sweep stage 13/15. This is a finite stage
positivity failure, distinct from the earlier contraction issue and the
accepted B4 sharp-contact limitation. Strict checks remain enabled.
The short shared interval `[0.5,2.5]` is dominated by startup, and active-only
`[2,6]` cannot establish the requested late-time active/passive contrast.
The [reviewed interrupted-run report](validation/cgl_lf_matched_physics_benchmark/interrupted_20261007/report.md)
retains the figures, measurements, uncertainty definitions, provenance and
17.449 allocated node-hour cost. It is separate from any replacement run.

The initial passive plumbing run
exposed a one-ULP pressure-decode contraction difference across GPU kernels.
The explicit contraction-order fix preserves the strict firehose wall and
passes its CPU/GPU regression, the failed-checkpoint replay, and fresh active
and passive 3D plumbing through `t=0.5`. The unpatched GPU negative control
fails the independent wall roundtrip. These are numerical checks, not a
completed turbulent comparison. Keep the failing preflight and its provenance.

The existing [96-grid PLM experiment](cgl_lf_physics_benchmark.md) and its
retained reports remain a frozen reference. Changing both resolution and
reconstruction means the new pair is **not a resolution-convergence test** of
that reference. Neither experiment independently validates LF coefficients.

## Models and matched parameters

Use the single [matched input](../inputs/cgl_lf_paper/cgl_lf_physics_benchmark_matched_beta10.athinput).
The [launcher](../scripts/run_cgl_lf_matched_benchmark.py) selects the member by
setting both `mhd/passive` and `problem/passive_delta` to the same boolean.
Do not override either flag separately or switch modes across a restart.

| Setting | Shared value |
| --- | --- |
| Box and grid | `(1,1,2)`, `192×192×384`, spacing `1/192`, periodic throughout |
| Blocks and halos | 8 blocks of `96×96×192`, one per GPU on one eight-GPU node; `nghost=3` |
| Initial state | `rho=1`, `B=(0,0,1)`, `u=0`, `p_parallel=p_perp=5`, initial beta ten |
| Numerics | `ppm4`, HLLE, RK2/RKL2, CFL `0.3`, STS safety `0.9`, merged sweeps off, FOFC off |
| LF | Local coefficients, `lf_k_parallel=2*pi`, safe arithmetic, weighted fluxes, detailed q diagnostics off; strict admissibility and pressure-work recording on |
| Collisions and limiters | `nu_coll=0`, soft thresholds `X=-2,+1`, finite rate `1e10`, `limiter_hardwall=false`, backups off |
| Sound-speed parameter | Explicit `iso_sound_speed=sqrt(5)` in both inputs |
| Forcing | OU, seed `271828`, `tcorr=2`, `dt_update=0.01`, continuous driving |
| Modes and mixture | Type zero, physical shell `pi<=|k|<=3*pi`, full signed bounds `-3..3`, power-law `expo=2`, solenoidal/compressive amplitude blend `1/(1+sqrt(2))` |
| Injection | `normalization=edot`, `dedt=0.16` per volume per time; nominal total-box power `0.32` |
| Initial target | Passive numerical gate at `t=3`, then the same passive lineage to `14` if healthy, followed by fresh active to `14`; no cycle cap; common analysis interval `[6,14]` |

The physical firehose wall `p_perp-p_parallel>=-B²` remains active even with
backups off. Strict checks and thresholds must not be relaxed to complete a
run. Both strict B4 expected failures remain separate sharp-contact limitations.

In the **active** member the CGL pressure tensor acts on momentum; total energy
and conservative anisotropy A are evolved. In the **passive** member the flow
is isothermal MHD with dynamical pressure `p_dyn=c_iso²*rho=5*rho`.
The separately evolved CGL thermal pressures use J/A material invariants,
including LF transport, collisions and limiters, but do not act on momentum.
See the [passive model contract](source/modules/cgl_passive.md).
The active member retains the same `iso_sound_speed` input for explicit
matching; this parameter does not replace its CGL characteristic speeds.

This is a comparison of the two complete dynamical models, including their
different thermal feedback. Equal input parameters do not imply identical
adaptive timesteps, realized density, Mach number or evolving thermal beta.
There is no thermal thermostat for the CGL pressures; initial beta ten is not
a fixed-beta ensemble.

The seed and fixed OU update interval match the modal stochastic construction.
They do **not** make the applied acceleration identical in the two flows:
normalization to `dedt` depends on density, velocity and the actual timestep.
Retain each member's forcing snapshots and accumulated `force_work`, and measure
the realized Helmholtz mixture. The blend gives equal expected innovation
powers, not an exact instantaneous half-and-half partition. The
[original guide](cgl_lf_physics_benchmark.md#setup-units-and-references) explains
the modal spectrum and normalization convention. Its `dedt=0.32` experiment
is historical; the new pair uses `0.16`. The total-power normalization now
matches MKS24's stated `0.32`, while the forcing projection and other model
differences remain.

## Build and retain the executable

Use `-DPROBLEM=built_in_pgens` with the input's
`problem/pgen_name=cgl_lf_paper`. The relevant implementation is
[`src/pgen/tests/cgl_lf_paper.cpp`](../src/pgen/tests/cgl_lf_paper.cpp).
The older `-DPROBLEM=cgl_lf_paper` custom generator has a different interface
and history and is not this experiment.

Use the complete [Frontier compiler/runtime contract](validation/wo2/README.md#task-0-establish-the-actual-gpu-baseline)
and the [fresh-build procedure](cgl_lf_physics_benchmark.md#build-the-correct-problem-generator),
with the build directory under the new matched artifact root:

```bash
export BENCH_SOURCE=/autofs/nccs-svm1_home2/dfielding/athenak-cgl-wo2
export BENCH_ROOT=/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/physics-benchmark
export PAIR_ROOT="$BENCH_ROOT/matched/power032"
mkdir -p "$PAIR_ROOT/tmp" "$PAIR_ROOT/cache" "$PAIR_ROOT/mpl-cache"
cd "$PAIR_ROOT"
export TMPDIR="$PAIR_ROOT/tmp" XDG_CACHE_HOME="$PAIR_ROOT/cache"
export MPLCONFIGDIR="$PAIR_ROOT/mpl-cache" PYTHONDONTWRITEBYTECODE=1
```

For the fresh lower-power pair, reuse the retained committed timing build below;
do not compile the parked adaptive-controller work. For an independent rebuild
of that same numerical revision, load the modules and clear inherited compiler
include overrides as described in that procedure. The configuration is:

```bash
cmake -S "$BENCH_SOURCE" -B "$PAIR_ROOT/build-hip" \
  -DPROBLEM=built_in_pgens -DCMAKE_BUILD_TYPE=Release \
  -DAthena_SINGLE_PRECISION=OFF -DAthena_ENABLE_MPI=ON \
  -DKokkos_ENABLE_HIP=ON -DKokkos_ARCH_ZEN3=ON -DKokkos_ARCH_VEGA90A=ON \
  -DCMAKE_CXX_COMPILER=CC -DCMAKE_CXX_FLAGS="-fno-cray -mno-daz-ftz" \
  -DCMAKE_EXE_LINKER_FLAGS="-no-pie"
cmake --build "$PAIR_ROOT/build-hip" --parallel 8
```

Retain the exact committed numerical revision, source status, Kokkos revision,
commands, build logs, cache, compiler/modules and executable SHA256 in a build
manifest. The launcher requires at least `revision` and `binary_sha256`, checks
the executable against them, and checks that the current numerical source
matches the recorded revision. An audited replacement-object build must also
retain its unchanged-input verification and original cache. A CMake cache alone
does not identify such an executable. Use the same retained executable and
build provenance for both members, including every continuation.

The retained executable has revision
`635b42562e5e23fc9966c9faa01c9da62134f2ec` and SHA256
`fb11d6355486dce126be1d4564f99dbcc0e9a13e7227a05f3a2f262ae466d47a`.
It contains the original committed solver and the verified interval timing
logging, not adaptive LF trial stepping. Its numerical plumbing checks do not
establish success of this fresh scientific pair. Use its retained provenance:

```bash
export PAIR_BINARY="$BENCH_ROOT/matched/build-timing-hip/athena"
export PAIR_BUILD_MANIFEST="$BENCH_ROOT/matched/build-timing-hip/manifest.json"
export PAIR_BUILD_CACHE="$BENCH_ROOT/matched/build-timing-hip/CMakeCache.txt"
```

## Launch one member per new directory

Obtain an appropriate GPU allocation and preserve its scheduler record. Load
the modules and exports from the
[self-contained runtime block](cgl_lf_physics_benchmark.md#run-resume-and-retain-provenance),
including `HSA_XNACK=1`, GPU-aware MPICH and `FI_CXI_ATS=0`. Restore the matched
`TMPDIR`, `MPLCONFIGDIR` and `XDG_CACHE_HOME` above afterward if needed.
The launcher records the environment; it does not load modules or establish it.

The following commands use one node and eight ranks/GPUs, or one block per
rank. They do not reserve an allocation themselves. Preserve the same rank
layout for both members and their continuations. The prescribed order is the
passive gate, passive continuation if healthy, then fresh active. These are
sequential launches; any performance comparison requires its own controlled
timing protocol. All new segments live under `matched/power032`, leaving the
older `matched/active-to14` and `matched/passive-to14` evidence untouched.

Start passive from the initial state:

```bash
python3 "$BENCH_SOURCE/scripts/run_cgl_lf_matched_benchmark.py" \
  "$PAIR_ROOT/passive-initial-to3" --mode passive \
  --executable "$PAIR_BINARY" --build-manifest "$PAIR_BUILD_MANIFEST" \
  --build-cache "$PAIR_BUILD_CACHE" --nodes 1 --ranks-per-node 8 \
  --job-id "$SLURM_JOB_ID" \
  --classification "lower-power passive numerical gate; not physical validation" \
  --set turb_driving/dedt=0.16 --set time/tlim=3 --set time/nlim=-1
```

Before continuing, inspect the exit status, final log/checkpoint time and strict
admissibility counters. Reaching `t=3` without pressure/density-floor, nonfinite
or nonpositive events is a numerical gate; inspect pressure minima and timestep
behavior as well. Intermediate LF hard-bound crossings are separate activity,
not a substitute for those failure checks. A clean wall-clock stop short of
three is incomplete. Stop and report a failure without relaxing the checks or
changing the physics. Passing this gate does not establish temporal accuracy,
stationarity or agreement with the papers.

After a healthy gate, resume its actual complete checkpoint with the same
forcing state and solver:

```bash
python3 "$BENCH_SOURCE/scripts/run_cgl_lf_matched_benchmark.py" \
  "$PAIR_ROOT/passive-continue-to14" --mode passive \
  --executable "$PAIR_BINARY" --build-manifest "$PAIR_BUILD_MANIFEST" \
  --build-cache "$PAIR_BUILD_CACHE" --nodes 1 --ranks-per-node 8 \
  --job-id "$SLURM_JOB_ID" \
  --restart "$PAIR_ROOT/passive-initial-to3/rst/ACTUAL_CHECKPOINT.rst" \
  --set turb_driving/dedt=0.16 --set time/tlim=14 --set time/nlim=-1
```

Then launch active from its prescribed initial state:

```bash
python3 "$BENCH_SOURCE/scripts/run_cgl_lf_matched_benchmark.py" \
  "$PAIR_ROOT/active-to14" --mode active \
  --executable "$PAIR_BINARY" --build-manifest "$PAIR_BUILD_MANIFEST" \
  --build-cache "$PAIR_BUILD_CACHE" --nodes 1 --ranks-per-node 8 \
  --job-id "$SLURM_JOB_ID" \
  --set turb_driving/dedt=0.16 --set time/tlim=14 --set time/nlim=-1
```

Retain the two passive segments in an aggregate manifest for analysis:

```bash
mkdir "$PAIR_ROOT/passive-union14"
printf '%s\n' '{"schema_version":1,"segments":["../passive-initial-to3","../passive-continue-to14"]}' \
  > "$PAIR_ROOT/passive-union14/benchmark_metadata.json"
```

Use new segment names if either directory contains previous data. Add
`--wall-time HH:MM:SS` when a clean wall-clock checkpoint is needed before an
allocation ends; choose a value shorter than the remaining allocation time.
This does not change the physical target. A zero exit code following a
wall-clock stop does not establish that `t=14` was reached: inspect the final
log and checkpoint time. Do not analyze missing endpoint coverage as a full
window.

Each launch retains canonical/effective input, copied launcher, command,
executable/cache/manifest hashes, mode, resource layout, optional restart hash,
environment, completion status and output inventory. Keep the raw files:

| Output | Cadence and meaning |
| --- | --- |
| `*.user.hst`, `*.mhd.hst` | `0.02`; full printed precision, physical integrals, forcing work and LF health counters |
| `bin/*.mhd_w_bcc.*.bin` | `0.25`; float32 primitives and cell-centered B, no ghosts; `eint` denotes parallel pressure |
| `bin/*.turb_force.*.bin` | `0.25`; acceleration components for actual Helmholtz decomposition |
| `rst/*.rst` | `1.0` plus termination outputs; double state and restart-persistent forcing/counters |

Passive `thermal-U` is physical thermal energy, whereas `cgl-J` is a conserved
material invariant, not total energy. Passive `force_work` records kinetic
energy injected by forcing kicks. The active `Delta tot-E - Delta force_work`
budget is not applicable to an isothermal passive flow; do not manufacture that
closure from J or from passive thermal pressure.

### Read the performance log

The current driver ports interval throughput from `scaling-tests`, commit
`46a6f704b563587025d0faa87fdd7a1623458d57`, and also reports wall seconds per
cycle. A progress line retains `elapsed`, `cycle`, `time`, and `dt`, then adds
`interval_cycles`, `wall_seconds_per_cycle`, and `zone-cycles/s`.
Set `time/ndiag=1` for individual cycles; the canonical cadence of 100 gives
the mean over each reporting interval. Initialization and restart reset the
baseline; a no-work interval is labeled `performance_interval=warmup`.

With interval wall time \(\Delta t_w\), completed cycles \(\Delta n\), and
global active-zone updates \(\Delta N_z\), these quantities are
\(\Delta t_w/\Delta n\) and \(\Delta N_z/\Delta t_w\). The uniform grid has
14,155,776 active cells, excluding ghosts. Throughput is global, not per GPU,
and counts complete cycles rather than RK or RKL stages. Divide by the
retained GPU or node count when comparing resource efficiency. Aggregate
unequal intervals using total work divided by total wall time.

Timing uses rank zero's elapsed wall clock and adds no fence, barrier, or MPI
reduction. It includes communication and output between reports. Intervals
containing snapshots/restarts can be slower; setup before execution is
excluded. The historical final labels `cpu time used` and
`zone-cycles/cpu_second` also refer to elapsed wall time, not CPU core seconds.
Compare physical time advanced per wall hour as well as zone cycles: LF
stiffness and changing stage counts alter work per cycle. Report scheduler
allocation node-hours separately, including idle time after a member fails.

The [retained timing check](validation/cgl_lf_matched_physics_benchmark/timing_20261007/report.md)
measured startup medians of 0.415 s/cycle active and 0.528 s/cycle passive on
the one-node layout. It also verified restart timing and zero floor/invalid
pressure counters. These short checks establish memory use and logging, not
late-time stability or a controlled scaling comparison. The lower-power pair
still requires its passive numerical gate; reaching the old failure time with
different forcing does not reproduce the old thermal state or establish a fix.
When launching concurrent members in a larger allocation, optional
`--nodelist NODE` pins a launcher to its intended node and is retained in the
run command and metadata.

## Resume each member's own lineage

Select a complete checkpoint from that member with the same mesh, PPM/halo,
physics and forcing configuration. Verify its header and time; do not use the
old 96/PLM reference, a higher-power `dedt=0.32` checkpoint, or a reduced-grid
plumbing checkpoint. Changing injection requires the fresh initial state, not
resuming the old higher-power trajectory. The checkpoint is
authoritative for its state and serialized parameters: specifying the new
canonical input cannot convert an incompatible checkpoint into this experiment.
The launcher rejects conflicting mode overrides and compares the checkpoint's
retained physical/numerical settings against the intended member before launch.
On resume, `effective.athinput` records the checkpoint header plus explicit
command-line overrides, including serialized defaults and the passive encoding.
Preserve its forcing RNG and counters and use a new directory.

For an interrupted passive gate, retain target 3 until the gate is complete.
For later interrupted production segments, retain target 14; for a planned
extension set target 18. This active example extends a completed first segment:

```bash
python3 "$BENCH_SOURCE/scripts/run_cgl_lf_matched_benchmark.py" \
  "$PAIR_ROOT/active-to18" --mode active \
  --executable "$PAIR_BINARY" --build-manifest "$PAIR_BUILD_MANIFEST" \
  --build-cache "$PAIR_BUILD_CACHE" --nodes 1 --ranks-per-node 8 \
  --job-id "$SLURM_JOB_ID" \
  --restart "$PAIR_ROOT/active-to14/rst/ACTUAL_CHECKPOINT.rst" \
  --set turb_driving/dedt=0.16 --set time/tlim=18 --set time/nlim=-1
```

Repeat for passive in `passive-to18`, using its own checkpoint from
`passive-continue-to14` and `--mode passive`. Replace the
checkpoint placeholder with the actual completed file; its numeric filename
is an output counter, not sufficient evidence of physical time. If additional
interrupted segments were required, retain every segment in lineage order.

Create a separate aggregate directory for each member without modifying child
metadata. For example:

```bash
mkdir "$PAIR_ROOT/active-union18" "$PAIR_ROOT/passive-union18"
printf '%s\n' '{"schema_version":1,"segments":["../active-to14","../active-to18"]}' \
  > "$PAIR_ROOT/active-union18/benchmark_metadata.json"
printf '%s\n' '{"schema_version":1,"segments":["../passive-initial-to3","../passive-continue-to14","../passive-to18"]}' \
  > "$PAIR_ROOT/passive-union18/benchmark_metadata.json"
```

Keep earlier aggregate manifests and analyses unchanged. Restart branching,
duplicate times and actual endpoint coverage are audited by the analyzer.

## Compare the same interval and bins

The predeclared first comparison is **[6,14]**, four contiguous two-unit blocks.
The blocks describe variability; they are not assumed independent. Examine
achieved Mach/beta, forcing mixture, energy, occupancy trends and correlation
times. If support remains inadequate, extend **both** members to 18 and analyze
[6,18], preserving [6,14]. Do not move the start time or select a favorable
interval to match a figure.

Install the local dbfplot package into an analysis environment under the
artifact root, as in the
[analysis environment instructions](cgl_lf_physics_benchmark.md#averaging-analysis-and-interpretation).
Use that environment for the comparison:

```bash
python3 "$BENCH_SOURCE/scripts/compare_cgl_lf_physics_benchmark.py" \
  "$PAIR_ROOT/active-to14" "$PAIR_ROOT/passive-union14" \
  --time-start 6 --time-end 14 --block-duration 2 --near-width 0.05 \
  --output-dir "$PAIR_ROOT/comparison-6-14"
```

For an extension, substitute the two union directories, end 18 and a new output
directory. The comparison reads the retained data; it never launches simulations.
It checks shared effective mesh, numerical, closure and forcing parameters,
allowing the explicit mode/encoding differences and listed output/stop controls.
Both members must cover the same requested interval. A parameter match alone
does not certify solver health or statistical convergence.

The comparison scans both primitive streams, including endpoint-bracketing
snapshots, to choose common **true histogram edges**. It passes those edges to
the single-run analyzer through `--pdf-edges`; clipped tails are rejected rather
than silently renormalized. B-strength PDFs use a **linear x-axis** and a
logarithmic density axis. Solid/dashed curves distinguish active/passive;
quantities retain consistent colors and physical normalization in dbfplot.

Main outputs are `report.md`, `metrics.json`, `shared_pdf_edges.json`,
`figure-audit.json` and PNG/PDF pairs `pdfs`, `occupancy`, `spectra`, and
`pressure_balance_scale`. The `active/` and `passive/` subdirectories retain
their individual metrics, reports and supplementary figures. Shared bins,
source/input/output hashes, matching exceptions, signed contrasts and block
support remain machine-readable.

## What the scale-dependent pressure comparison measures

Set `a=p_perp-<p_perp>` and `b=p_B-<p_B>`, where `p_B=B²/2`. For each physical
perpendicular shell compute:

```text
Paa = sum_shell |FFT(a)/N|² / dk
Pbb = sum_shell |FFT(b)/N|² / dk
Pab = sum_shell Re[(FFT(a)/N) * conjugate(FFT(b)/N)] / dk
C = Pab / sqrt(Paa*Pbb)
R = (Paa + Pbb + 2*Pab) / (Paa+Pbb)
```

Average Paa, Pbb and signed Pab in physical time **before** forming C and R,
separately for the full interval and each complete block. Thus C is signed
shell correlation, not squared coherence or a mean of snapshot ratios. Exact
equal-amplitude compensation gives C=-1,R=0. C=-1 with R>0 reveals an amplitude
mismatch. Zero denominators are unavailable, not assigned an artificial value.

Both members use common `dk=2*pi`, shell sums divided by dk, and midpoint
coordinates. Retain the first perpendicular shell in metrics and Parseval
sums, but omit it from logarithmic-wavenumber plots: it contains pure-parallel
modes, despite its positive bin midpoint. Main perpendicular spectra and C/R
curves sum all parallel wavenumbers. The dotted transverse eight-cell scale
and gray region above 0.75 times the transverse Nyquist wavenumber are plotting
guides, not measured boundaries of an inertial range. A separate counterpart
retained in metrics restricts
full `0<|k|<=min(Nyquist_xyz)/4`; use the same mask on all three powers.
Temporal block ranges are descriptive spread, not confidence intervals.

This C/R comparison uses physical **thermal** p_perp in both members.
In passive dynamics that pressure does not push the fluid. The separately
labeled `c_iso²*rho` spectrum and dynamic-pressure balance diagnose its actual
isothermal momentum pressure. Likewise report the passive isothermal Mach
number separately from the common CGL thermal-pressure sound-speed proxy.
With volume averages, their definitions are

```text
u_rms² = <|u-<u>|²>, p_iso=(p_parallel+2*p_perp)/3
M_proxy = u_rms / sqrt[(5/3)*<p_iso>/<rho>]
M_iso = u_rms / c_iso                    (passive dynamics only)
```

The thermal proxy is not a CGL characteristic-wave Mach number. Neither
definition is silently retuned to a target Mach number as the thermal state
drifts.

The velocity-gradient spectra differentiate u before projecting the tensor
onto the local magnetic direction. They are not spectra of the magnitudes of
derivatives of a previously projected velocity. Centered differences attenuate
short wavelengths, perpendicular shells near the square-grid corner have
incomplete annuli, and a short spectrum is not a reliable inertial-range fit.
The normalization audit found correct shell sums, component powers and
Parseval identities; it did not establish that all measured steepening was
physical or quantify reconstruction error.
Specifically, with `G_ij=partial_j u_i` and `P=I-b*b`, the four spectra use
`b.G.b`, `P.G.b`, `b.G.P` and `P.G.P`, summing vector or tensor component
powers. This follows the local-gradient convention of Squire Eq. 24. MKS24
Figure 6b uses line styles to compare forcing types; our line styles instead
identify active/passive dynamics. The positive one-third gradient and
negative five-thirds energy/pressure slopes are illustrative guides, not
universal targets or fitted acceptance criteria.

## Interpretation and retained limits

[Squire et al. (2023)](https://arxiv.org/html/2303.00468v2) motivate active/passive
comparisons of magnetic-strength fluctuations, anisotropy and flow gradients.
[Majeski, Kunz & Squire (2024)](https://arxiv.org/html/2405.02418v2) provide
pressure-balance, gradient and spectral comparisons. These are scientific
reference points, not pixel targets or numerical acceptance bands. In
particular, the forcing geometry still differs. The fresh pair prescribes total
power `0.32` (`dedt=0.16` per volume), matching MKS24's stated total-power
normalization. The earlier total-power `0.64` runs remain separately labeled
historical experiments; neither their measurements nor their validation status
transfer automatically to the fresh pair.

Evaluate the diagnostic groups together, including achieved regime, numerical
health, sampling and model differences. Label a comparison consistent,
concerning or inconclusive with its reason; do not declare a closure valid,
an asymptotic slope converged, or a solver defective from visual agreement or
disagreement alone. Passing the repaired passive preflight does not establish
any of these physical conclusions.
