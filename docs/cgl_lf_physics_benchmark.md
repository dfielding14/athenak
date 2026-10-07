# CGL–Landau-fluid turbulence physics benchmark

This benchmark asks whether one forced, active-CGL box exhibits a mutually
consistent set of pressure-anisotropy, pressure-balance, velocity-strain and
energy statistics. It complements the independent wave, decay, conservation
and restart tests in [WO2 validation](validation/wo2/README.md). One box cannot
establish resolution convergence, active-versus-passive causality, or the
accuracy of an LF closure coefficient. Agreement with a turbulence figure is
not a substitute for those tests.

Changing field strength and density drives different parallel and perpendicular
pressures. Their anisotropic stress then changes the flow that produced them,
while parallel heat transport and threshold scattering modify the pressures.
The diagnostics examine these linked effects: residence near instability
thresholds, perpendicular pressure compensation, field-strength-changing
strain, and the achieved forcing/energy regime.

Use the [canonical input](../inputs/cgl_lf_paper/cgl_lf_physics_benchmark_beta10.athinput).
Extend its sampling interval without retuning its physical parameters. This is
a finite-time, initially beta-ten experiment with no cooling; its thermal state
is not held at beta ten.

## Setup, units and references

The code absorbs the electromagnetic factor into B: magnetic pressure is
`p_B=B²/2` and `v_A0=B0/sqrt(rho0)`. Units are `L_perp=1`, `v_A0=1`, `rho0=1`,
pressure `rho0*v_A0²=1`, and time `L_perp/v_A0=1`. The mean field points along z;
`L_parallel=2` and the box volume is two.

| Setting | Canonical choice |
| --- | --- |
| Domain/grid | Periodic `(1,1,2)`, `96×96×192`, spacing `1/96`; 64 blocks of `24×24×48` |
| Initial state | `rho=1`, `B=(0,0,1)`, `u=0`, `p_parallel=p_perp=5` |
| Evolution | Active CGL, PLM/HLLE, RK2/RKL2, CFL `0.3`, STS safety `0.9`, merged sweeps off |
| LF closure | Local coefficients, `lf_k_parallel=2*pi`, safe arithmetic, weighted fluxes, full diagnostics |
| Collisions/limiters | `nu_coll=0`, thresholds `X=-2,+1`, finite `limiter_nu_coll=1e10`, `limiter_hardwall=false`, backups off |
| OU forcing | `tcorr=2`, `dt_update=0.01`, seed `271828`, continuous driving |
| Force geometry | Type zero, `solenoidal_compressive`, physical shell `pi<=|k|<=3*pi`, `power_law`, `expo=2` |
| Force mixture | `sol_fraction=1/(1+sqrt(2))=0.4142135623730951`, an amplitude blend |
| Injection | `normalization=edot`, `dedt=0.32` per volume per time |
| Initial duration | `tlim=10`, no cycle limit; predeclared analysis `[6,10]`, blocks of duration two |

Strict LF admissibility and pressure-work recording remain on; `fofc=false`
preserves the input's stated reconstruction behavior. Finite limiter relaxation
preserves `U=p_perp+p_parallel/2` and relaxes anisotropy toward the activated
threshold. Backups off still retains the unconditional physical firehose wall
`p_perp-p_parallel>=-B²`. Thus this is not unconstrained finite-rate firehose
transport. See [the CGL helpers](../src/eos/cgl_physics.hpp).

The driver blends acceleration amplitudes as `s*f_sol+(1-s)*f_comp`.
Two transverse and one longitudinal degrees of freedom give equal expected
innovation powers when `2*s²=(1-s)²`. The parameter is not a solenoidal energy
fraction, and the realized finite-time mixture must be measured. Full signed
mode bounds `-3..3` give 38 enumerated entries, or 19 opposite-wavevector pairs,
in this physical shell. This differs from the historical positive-octant
selection and from unprojected random forcing.

[NormalizeForce](../src/srcterms/turb_driver.cpp) divides both work moments by
volume before solving its normalization quadratic. Therefore `dedt=0.32`
corresponds to nominal total-box power `0.64`. Measure actual injection using
increments of `force_work`. Keeping `dedt=0.32` is a deliberate benchmark choice;
it does not match MKS24's stated total-box power of `0.32`.

| Primary reference | Comparison and qualification |
| --- | --- |
| [Squire et al. (2023), Fig. 3](https://arxiv.org/html/2303.00468v2) | Magnetic-strength and anisotropy PDFs, and threshold occupancy. Its anisotropy axis is `4*pi*Delta p/B_phys²=X/2`; its active/passive comparison and forcing selection are separate experiments. |
| [Majeski, Kunz & Squire (2024), Figs. 2b, 5a, 6a–b](https://arxiv.org/html/2405.02418v2) | Unstable volume, velocity spectra, thermal/magnetic-pressure spectra, and local-field gradients. Their random forcing is unprojected and their stated injection normalization is total-box power. |
| [MKS24, Fig. 13b](https://arxiv.org/html/2405.02418v2) | A beta-100 limiter-rate scan, useful context rather than a target PDF for this beta-ten realization. |

The supplied “Steve” comparator with a cutoff near `X=-1.4` differs from this
input's `-2`. Record that mismatch rather than moving the canonical threshold.
Paper images support qualitative comparisons; they do not supply tight
numerical tolerances. The startup exclusion `t<6` is predeclared, consistent
with the developed-turbulence discussion in
[MKS24 §3.3](https://arxiv.org/html/2405.02418v2).

## Build the correct problem generator

Configure `-DPROBLEM=built_in_pgens` and select `problem/pgen_name=cgl_lf_paper`
in the input. The implementation is
[`src/pgen/tests/cgl_lf_paper.cpp`](../src/pgen/tests/cgl_lf_paper.cpp), using
`paper_mode=turbulence`, lowercase `b0`, and the 22-column user history below.
The older `-DPROBLEM=cgl_lf_paper` route selects `src/pgen/cgl_lf_paper.cpp`, with
a different interface and history. It is not this benchmark build.

The unchanged numerical implementation can reuse the released Frontier binary:

```text
/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/final-research/release-bin/athena-hip
SHA256 6530c1b31a7d99077f3ff5932d4b1c790e33be38c1fa222cfda097f87b339f9e
Numerical source revision 7a37710f6, with retained release-source audit
```

For a fresh build, use the committed Kokkos revision and the complete
[Frontier compiler/runtime contract](validation/wo2/README.md). After loading
that environment:

```bash
BENCH_SOURCE=/autofs/nccs-svm1_home2/dfielding/athenak-cgl-wo2
BENCH_ROOT=/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/physics-benchmark
mkdir -p "$BENCH_ROOT/build-hip/tmp"
export TMPDIR="$BENCH_ROOT/build-hip/tmp"
cmake -S "$BENCH_SOURCE" -B "$BENCH_ROOT/build-hip" \
  -DPROBLEM=built_in_pgens -DCMAKE_BUILD_TYPE=Release \
  -DAthena_SINGLE_PRECISION=OFF -DAthena_ENABLE_MPI=ON \
  -DKokkos_ENABLE_HIP=ON -DKokkos_ARCH_ZEN3=ON -DKokkos_ARCH_VEGA90A=ON \
  -DCMAKE_CXX_COMPILER=CC -DCMAKE_CXX_FLAGS="-fno-cray -mno-daz-ftz" \
  -DCMAKE_EXE_LINKER_FLAGS="-no-pie"
cmake --build "$BENCH_ROOT/build-hip" --parallel 8
```

## Run, resume and retain provenance

All builds, temporary files, logs and analysis belong under `BENCH_ROOT`.
The following self-contained procedure uses the released binary. For a new
build, set the binary, expected hash, numerical source revision, dirty-state
record and build cache to that build's retained provenance. Do not substitute
the current analysis checkout's revision for the simulation binary's revision.
Acquire a site-approved eight-GPU allocation first, and retain its allocation
record. The application step below assigns eight blocks per GPU.

```bash
BENCH_SOURCE=/autofs/nccs-svm1_home2/dfielding/athenak-cgl-wo2
BENCH_ROOT=/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/physics-benchmark
mkdir -p "$BENCH_ROOT"
cd "$BENCH_ROOT"
module restore
module load PrgEnv-cray craype-accel-amd-gfx90a cpe/25.09 cray-mpich/9.0.1
module load rocm/6.4.2 cce/20.0.0 cray-python/3.11.7
module unload darshan-runtime
unset CPATH C_INCLUDE_PATH CPLUS_INCLUDE_PATH LIBRARY_PATH
# Ensure constructor environment overrides cannot silently change this input.
unset ATHENAK_CGL_LF_DIAGNOSTICS ATHENAK_CGL_LF_ARITHMETIC ATHENAK_CGL_LF_STS_FLUX
unset ATHENAK_CGL_LF_PROFILE ATHENAK_CGL_LF_PROFILE_DETAIL
export INCLUDE_PATH_X86_64="$(CC --print-resource-dir)/include:${CC_X86_64}/include/craylibs"
export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH:-}"
export MPICH_GPU_SUPPORT_ENABLED=1 MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=0
export MPICH_OFI_NIC_POLICY=GPU MPICH_GPU_IPC_CACHE_MAX_SIZE=1000 HSA_XNACK=1
export MPICH_MPIIO_HINTS='*:romio_cb_write=disable'
export MPICH_OFI_NUM_CQ_ENTRIES=131072 FI_MR_CACHE_MONITOR=kdreg2
export FI_CXI_RX_MATCH_MODE=software OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export TMPDIR="$BENCH_ROOT/tmp" MPLCONFIGDIR="$BENCH_ROOT/mpl-cache"
export XDG_CACHE_HOME="$BENCH_ROOT/cache"
mkdir -p "$TMPDIR" "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
export BENCH_SOURCE BENCH_ROOT
export BENCH_BINARY=/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/final-research/release-bin/athena-hip
export BENCH_CACHE=/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/final-research/release-build-hip/CMakeCache.txt
export BENCH_BINARY_SHA=6530c1b31a7d99077f3ff5932d4b1c790e33be38c1fa222cfda097f87b339f9e
export BENCH_SIM_REVISION=7a37710f6c224e24e7c7f364e7e0b812b3a9494c
export BENCH_SEGMENT=canonical BENCH_END=10 BENCH_RESTART=
# Use a new segment name if this directory already exists; never overwrite a run.
python3 - <<'PYRUN'
import hashlib, json, os, re, shutil, subprocess, time
from datetime import datetime, timezone
from pathlib import Path
e = os.environ
assert e.get("SLURM_JOB_ID"), "An allocated GPU job is required"
source, root = Path(e["BENCH_SOURCE"]), Path(e["BENCH_ROOT"])
binary = Path(e["BENCH_BINARY"])
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(binary) == e["BENCH_BINARY_SHA"]
run = root / e["BENCH_SEGMENT"]
run.mkdir(parents=True, exist_ok=False)
original = source / "inputs/cgl_lf_paper/cgl_lf_physics_benchmark_beta10.athinput"
shutil.copy2(original, run / "input.athinput")
shutil.copy2(e["BENCH_CACHE"], run / "CMakeCache.txt")
overrides = ["time/tlim=" + e["BENCH_END"], "time/nlim=-1"]
text = original.read_text()
for key, value in [("tlim", e["BENCH_END"]), ("nlim", "-1")]:
    text, count = re.subn(r"(?m)^"+key+r"\s*=.*$", key+" = "+value, text)
    assert count == 1
(run / "effective.athinput").write_text(text)
restart = Path(e["BENCH_RESTART"]).resolve() if e["BENCH_RESTART"] else None
arguments = ["-r", str(restart)] if restart else ["-i", "input.athinput"]
command = ["srun", "--exact", "-N", "1", "-n", "8", "--ntasks-per-node", "8",
           "--threads-per-core=1", "--cpu-bind=threads", "-c", "7",
           "--gpus-per-task=1", "--gpu-bind=closest", str(binary),
           *arguments, *overrides]
git = lambda *args: subprocess.check_output(["git", "-C", str(source), *args], text=True).strip()
metadata = {"schema_version": 1, "simulation": {
    "revision": e["BENCH_SIM_REVISION"], "dirty": False,
    "revision_basis": "Released binary and numerical source audit",
    "executable": str(binary), "executable_sha256": sha(binary),
    "input_path": "effective.athinput", "input_sha256": sha(run/"effective.athinput"),
    "canonical_input_path": "input.athinput", "canonical_input_sha256": sha(original),
    "build_cache_path": "CMakeCache.txt", "build_cache_sha256": sha(run/"CMakeCache.txt")},
    "outputs": {"snapshot_glob": "bin/*.mhd_w_bcc.*.bin",
                "forcing_glob": "bin/*.turb_force.*.bin"},
    "launch": {"command": command, "overrides": overrides, "cwd": str(run),
        "ranks": 8, "gpus": 8, "slurm_job_id": e["SLURM_JOB_ID"],
        "restart": str(restart) if restart else None,
        "restart_sha256": sha(restart) if restart else None,
        "checkout_revision": git("rev-parse", "HEAD"),
        "checkout_status": git("status", "--short", "--ignore-submodules=all"),
        "environment": {k:v for k,v in e.items() if k.startswith(("MPICH_", "FI_", "HSA_", "OMP_"))
                        or k in ("LOADEDMODULES", "LD_LIBRARY_PATH")},
        "started_utc": datetime.now(timezone.utc).isoformat()},
    "planned_analysis": {"initial_window": [6,10], "block_duration": 2}}
def save():
    (run/"benchmark_metadata.json").write_text(json.dumps(metadata, indent=2)+"\n")
save()
start = time.monotonic()
with (run/"run.log").open("w") as log:
    result = subprocess.run(command, cwd=run, stdout=log, stderr=subprocess.STDOUT)
metadata["launch"].update(returncode=result.returncode,
    ended_utc=datetime.now(timezone.utc).isoformat(), wall_seconds=time.monotonic()-start,
    binary_hash_unchanged=(sha(binary) == e["BENCH_BINARY_SHA"]))
for key, suffix in [("user_history", "*.user.hst"), ("mhd_history", "*.mhd.hst")]:
    metadata["outputs"][key] = [str(p.relative_to(run)) for p in sorted(run.glob(suffix))]
save()
raise SystemExit(result.returncode)
PYRUN
```

Retain `module -t list`, `scontrol show job "$SLURM_JOB_ID"`, compiler/runtime
versions, GPU/driver information and the canonical forcing-mode audit with the
manifest. The executed experiment's scratch `scripts/stage_run.py` and
`runtime.sh` preserve equivalent launch evidence; they are not required by the
self-contained procedure above.

To extend, select a complete checkpoint whose log/header confirms the intended
physical time, then repeat the Python capture-and-launch block with:

```bash
export BENCH_SEGMENT=canonical-to14 BENCH_END=14
export BENCH_RESTART=/absolute/path/to/canonical/rst/cgl_lf_physics_benchmark_beta10.XXXXX.rst
```

Replace `XXXXX` with the actual checkpoint counter. Keep the same binary, rank
layout, forcing parameters and RNG state. The `-r` launch preserves forcing
phase and output counters; it does not regenerate initial conditions. Preserve
both segments and their original metadata. The effective input records physical
parameters and explicit overrides; the checkpoint remains authoritative for
its state and serialized counters.

## Outputs and weighting

Histories are written every `0.02`, fields/acceleration every `0.25`, and
restarts every `1.0` in physical time. Binary field snapshots are float32;
restart state is double. The canonical files are common MPI files and the
analyzed snapshots contain no ghost zones.

| File | Meaning |
| --- | --- |
| `bin/*.mhd_w_bcc.*.bin` | `dens`, `velx/vely/velz`, legacy `eint` meaning **parallel pressure**, `p_perp`, `bcc1/bcc2/bcc3` |
| `bin/*.turb_force.*.bin` | Applied acceleration `force1/force2/force3`; use the matching physical time |
| `*.user.hst` | Built-in problem's volume integrals and accumulated forcing work |
| `*.mhd.hst` | Conserved quantities including `tot-E`, and LF health/work diagnostics |
| `rst/*.rst` | State, forcing RNG/phase, cumulative forcing work and output counters |

Spatial means are volume weighted; cells have equal volume here, ranks and
blocks do not. User columns `volume`, `mass`, `kinetic`, `magnetic`, `therm_cgl`
contain `V`, `integral rho`, `integral rho*u²/2`, `integral B²/2`, and
`integral(p_perp+p_parallel/2)`. Divide extensive quantities by V for densities.
`p_parallel`, `p_perp`, `delta_p`, `abs_dp`, `beta`, `b2`, `b4`, and `nu_eff`
are also volume integrals.

`vel_prp2/vel_prl2` integrate `ux²+uy²` and `uz²`; `force_prp2/force_prl2`
integrate analogous acceleration components. These use the fixed z direction,
not the local field or a Helmholtz projection. `mirror_vol/fire_vol` include
threshold equality; `hard_vol` counts violated physical/active backup bounds.
These are instantaneous post-operator samples, not every limiter activation
within a cycle. `force_pwr` is instantaneous `integral rho*u dot f`;
`force_work` is accumulated actually applied forcing energy and must not be
integrated in time again.

For this input (`nu_coll=0`, both soft limiters on, backups off),
`nu_eff/(1e10*V)` measures the strictly activated soft-rate fraction at the
history checkpoint. It differs from inclusive `mirror_vol/fire_vol` and from
near-threshold residence; none measures the fraction of time spent unstable
between outputs.

Check numerical health (`lf_dfloor`, `lf_pfloor`, `lf_nonfin`, `lf_nonpos`)
separately from intended physical scattering. `lf_hardbd` counts unprojected
LF-stage bound crossings, and `lf_hwproj` counts physical wall projections;
their activity alone is not a numerical failure. Strict mode checks every LF
stage for floor, nonfinite and nonpositive states, and checks hard bounds after
the physical wall. Post-operator `hard_vol` should be zero. Report these
different checkpoints explicitly rather than treating every limiter counter
as a failed admissibility check.
`lf_nstage` is a cell-stage total. `lf_qprwrk/lf_qpewrk` record face-flux work;
`lf_cpwrk/lf_cawrk` record configured pressure work. Offline gradients are not
an exact reconstruction of the limited face stencil.

## Four diagnostic groups

The analyzer produces `metrics.json`, `report.md`, and four PNG/PDF pairs:
`marginality`, `pressure_balance`, `gradients`, and `spectra_energy`. Each result
must retain definitions, selected times and sampling support. Undefined
correlations or ratios remain explicitly unavailable. `metrics.json` retains
input/output hashes and the analysis, reader and binary-reader script hashes;
the analysis checkout revision is separate from the simulation revision in the
launch metadata.

**Marginality.** Define `Delta p=p_perp-p_parallel`,
`p_iso=(p_parallel+2*p_perp)/3`, `beta=2*p_iso/B²` and
`X=beta*Delta p/p_iso=2*Delta p/B²`, with local cell-centered B. Plot volume PDFs
of `|B|/B0` and X. Density-normalized histograms integrate to one with bin widths;
report clipped tails. Separate strict exceedance `X<-2` or `X>1`, inclusive
threshold occupancy, and symmetric near-threshold bands
`|X-X_threshold|<=0.05`. The band width is a diagnostic convention, not a limiter
or acceptance tolerance. Report its interior-only part separately. The analyzer
also propagates half-ULP float32 storage intervals for pressure and B into an
X rounding envelope. Crossings outside that envelope withstand the stated
storage rounding assumption; the envelope is not a dynamical-error bound.
Use double history/restart evidence when interpreting small crossings.

**Pressure balance.** Subtract each snapshot's volume mean and compare spectra
of `p_parallel`, `p_perp`, and `p_B=B²/2`. Pressure-spectrum units are
pressure squared times length. The spectrum of magnetic pressure is distinct
from magnetic energy. Matching auto-spectra does not establish compensation;
also report the signed correlation and normalized residual:

```text
C = <delta p_perp * delta p_B> / sqrt(<delta p_perp²> * <delta p_B²>)
R = <(delta p_perp + delta p_B)²> / (<delta p_perp²> + <delta p_B²>)
```

Exact equal-amplitude cancellation gives `C=-1`, `R=0`. Do not independently
rescale pressure spectra before comparing their amplitudes. These signed
correlation and residual diagnostics are real-space statistics, not additional
spectral cross-power outputs.

**Gradients.** Use local `b=B/|B|` and second-order periodic centered derivatives
with actual mesh spacings. Compute signed strain
`S=sum_ij b_i*b_j*partial_j u_i`, divergence `D=div u`, and `S-D`, the
ideal-induction material derivative of `ln|B|`. Retain `u_parallel=u dot b`
and `b dot grad(u dot b)` separately:

```text
b dot grad(u dot b) = S + u dot [(b dot grad)b].
```

This is a continuum identity; the centered-difference product rule has a
discretization residual. The curvature term distinguishes the two gradients.
Small parallel velocity alone does not show suppressed field-strength-changing
strain. With `G_ij=partial_j u_i` and `P_ij=delta_ij-b_i*b_j`, the four local
gradient components are `b.G.b`, `P.G.b`, `b.G.P` and `P.G.P`. Apply the
projectors to the already differentiated velocity, as in
[Squire et al. Eq. 24](https://arxiv.org/html/2303.00468v2#S3.E24); these
components contain no derivatives of b. Differentiating a projected velocity
would add field-direction terms and measure something different. The analyzer
reports these spectra and RMS diagnostics, with separate `b.grad(u.b)`, `S-D`
and local parallel velocity; it does not produce signed gradient PDFs.
Near-grid-scale gradients depend strongly on the derivative stencil and output
precision.

**Spectra and energy.** Let `U=<rho*u>/<rho>` be the mass-weighted bulk velocity.
The kinetic spectral field `sqrt(rho)*(u-U)/sqrt(2)` represents turbulent kinetic
energy density; the magnetic field `(B-<B>)/sqrt(2)` represents magnetic
fluctuation energy density. Raw velocity variance is a separate comparison.
The kinetic transform retains its own mean, which need not vanish despite
subtracting mass-weighted bulk velocity; no second mean removal is applied.
The FFT is normalized by cell count. Shell powers are sums divided by `dk`,
not averages over mode count. Use physical `ki=2*pi*ni/Li` and perpendicular
shells `k_perp=sqrt(kx²+ky²)` relative to the initial z guide. Their width is
`dk=min(2*pi/Lx,2*pi/Ly)=2*pi` here, with bins `[n*dk,(n+1)*dk)`. Keep the
`k_perp=0` shell when checking integrated power; it may contain parallel modes.
Account explicitly for any removed mean of a transformed spectral field in its
Parseval energy check. Pressure and gradient spectra have no energy factor 1/2.
The acceleration Helmholtz decomposition uses the full physical wavevector.

Plot the forcing support and a conservative `k_Nyquist/4` guide, corresponding
to eight cells per wavelength; distinguish the perpendicular cutoff band
`k_perp>=0.75*k_Nyquist`. The scalar resolved strain ratio additionally selects
full `|k|<=min(k_Nyquist,x,k_Nyquist,y,k_Nyquist,z)/4` and excludes the first
perpendicular shell. Perpendicular shells otherwise sum over all parallel
wavenumbers, so small `k_perp` alone does not guarantee a resolved gradient.
Centered derivatives attenuate short wavelengths. These marks organize
interpretation and do not establish convergence or prescribe a fitted slope.

Report `delta B_rms/B0`, evolving beta, energies, and the declared isotropic
sound-speed proxy
`M_delta=sqrt(<|u-<u>|²>)/sqrt(gamma*<p_iso>/<rho>)`, with `gamma=5/3`.
The Mach numerator subtracts volume-mean velocity; the energy spectrum above
uses mass-weighted bulk velocity, so these are explicitly different weightings.
This is not a CGL wave speed. The beta history divided by V is
`<2*p_iso/B²>`, which differs from `2*<p_iso>/<B²>`. Distinguish vector-field
fluctuations from fluctuations of `|B|`.

Over synchronized history samples compare `Delta tot-E` to `Delta force_work`
and report their residual and normalization. Actual total power is
`Delta force_work/Delta t`, and volumetric power divides this by V. Energy
bookkeeping alone does not yield a viscous-heating fraction or cascade flux.

## Averaging, analysis and interpretation

The first window is **[6,10]**, declared before simulation, with contiguous
physical-time blocks of duration two. It supplies only two nominal blocks,
not two guaranteed independent samples. Weight snapshots and histories in
physical time with trapezoidal integration. Output times can overshoot their
nominal cadence: use actual stored times and require bracketing samples for
both requested endpoints. Interpolate diagnostic values, histogram bins and
spectral bins at those endpoints, never fields or unsupported extrapolations.
Anchor complete blocks at the requested start time. Check complete output inventories
and remove repeated restart-boundary samples with an explicit provenance rule.

```bash
python3 "$BENCH_SOURCE/scripts/analyze_cgl_lf_physics_benchmark.py" \
  "$BENCH_ROOT/canonical" --time-start 6 --time-end 10 --block-duration 2 \
  --output-dir "$BENCH_ROOT/analysis-6-10"
```

Inspect energy, Mach, beta, occupancy and pressure-balance trends, and estimate
temporal autocorrelation. Report block means, their range and sample standard
deviation as descriptive variability, not independent-sample confidence
intervals. At least four complete blocks are desirable; increase block duration
if correlations persist longer than two. Thermal heating may coexist with
settled velocity statistics: identify which quantities are approximately
stationary.

If support is inadequate, extend the **same realization** to 14, retaining the
original [6,10] report, and analyze [6,14] for four blocks. Extend to 18 if needed.
Create a separate aggregate directory without changing child metadata:

```bash
mkdir "$BENCH_ROOT/analysis-run"
cat > "$BENCH_ROOT/analysis-run/benchmark_metadata.json" <<'JSON'
{"schema_version": 1, "segments": ["../canonical", "../canonical-to14"]}
JSON
python3 "$BENCH_SOURCE/scripts/analyze_cgl_lf_physics_benchmark.py" \
  "$BENCH_ROOT/analysis-run" --time-start 6 --time-end 14 --block-duration 2 \
  --output-dir "$BENCH_ROOT/analysis-6-14"
```

The analyzer retains child provenance, trims superseded restart branches and
handles duplicate times. Check its retained-time inventory before interpreting
the expanded interval. An extension to 18 appends its segment to this explicit
chronological list. Do not retune forcing, lower resolution, reset initial
conditions, or select only
intervals that resemble a figure. If adequate sampling remains unavailable,
label the affected comparisons inconclusive.

| Classification | Required explanation |
| --- | --- |
| **Consistent** | Aligned definitions and achieved regime, adequate sampling, and qualitatively compatible signatures across retained blocks; list remaining differences. |
| **Concerning** | A reproducible inconsistency survives definition, provenance and sampling checks; identify its candidate cause without declaring a solver defect from a visual mismatch alone. |
| **Inconclusive** | Sampling, scale separation, rounding, forcing differences or trends prevent the intended inference; identify the missing evidence. |

Troubleshoot in this order: definitions/units/weighting and restart sampling;
measured forcing mixture/power, Mach/beta and fluctuation amplitudes; resolution
and forcing-to-dissipation scale separation; then implementation with an
independent small test. A single `96×96×192` run need not offer a long inertial
range. State what each diagnostic establishes instead of assigning tight
acceptance bands to paper images.
