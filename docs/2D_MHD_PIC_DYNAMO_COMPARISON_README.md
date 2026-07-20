# 2D MHD-PIC Turbulent Dynamo Comparison Setup

This note is for comparing the 2D MHD-PIC turbulent-box runs in this repo
against another code.

## Branch And Scope

Use branch:

```bash
git fetch origin
git checkout PIC
git pull --ff-only origin PIC
```

The relevant upstream branch is `origin/PIC`. The 2D comparison setup is in the
repo, not just in local scratch directories.

Main files:

- `inputs/particles/pic_turbulent_dynamo_512_2d_mhd_gamma1p00001.athinput`
- `inputs/particles/pic_turbulent_dynamo_512_2d_cr_fid_gamma1p00001.athinput`
- `inputs/tests/initial_perturbations_2d_inplane.athinput`
- `src/srcterms/initial_perturbations.cpp`
- `src/srcterms/initial_perturbations.hpp`

The initializer supports `magnetic_geometry = in_plane_2d`, which seeds a weak,
zero-net-flux, in-plane magnetic field on a 2D mesh. The CR deck uses 2D spatial
geometry but keeps all three particle velocity components and all three magnetic
field components active: 2D3V3B.

## Physics Setup

Common mesh and MHD settings:

- Domain: periodic box, `x1,x2 in [-0.5, 0.5]`, `nx1 = nx2 = 512`, `nx3 = 1`.
- Meshblocks: `256 x 128 x 1`, giving 8 meshblocks for one Frontier node with
  8 GPUs.
- Integrator/reconstruction/Riemann solver: `rk2`, `plm`, `hlld` for the
  MHD-only deck and `vl2`, `plm`, `hlld` for the coupled CR deck.
- Gas: nearly isothermal ideal gas with `gamma = 1.00001` and
  `iso_sound_speed = 1.0`.
- Driving: solenoidal turbulence, `constant_edot = true`, `dedt = 0.15`,
  forcing on modes `n = 1..3`, with target sonic Mach number around 0.5.
- Initial magnetic field: tangled in-plane field, zero net flux, `Bz = 0`,
  `magnetic_rms = 1.0e-3`, modes `n = 1..4`.

CR settings in the PIC deck:

- `particle_type = cosmic_ray`
- `ppc = 6.0`
- `pusher = boris_tsc`
- `pic_enable_2d3v = true`
- `pic_cr_light_speed = 100.0`
- `pic_cr_initial_state = velocity`
- Six species with velocities `(+/-5,0,0)`, `(0,+/-5,0)`, `(0,0,+/-5)`.
- `mass = 1.0`, `charge = 100.0` for each species.
- Coupled MHD-PIC controls: `time/integrator = vl2`, TSC push and deposition,
  `pic_background_mode = coupled`, `pic_feedback_mode = coupled`, and
  `pic_cr_hall_mode = off`. Set only the last control to `full` for a matched
  full-CR-Hall run.

The six-beam loading is intended as an initially nearly isotropic CR population.
It is not a sampled thermal distribution. The code converts the input velocities
to the internal momentum/Lorentz-factor state using the configured reduced speed
of light.

## Two Science Variants

Run both variants for a clean comparison:

1. `cr-from-ic`: start the CR deck from `t = 0`.
2. `cr-injected-at-saturation`: run the MHD-only deck first, then restart with
   the CR deck from the saturated or late-time MHD state.

The MHD-only deck currently has `time/tlim = 25.0`. The CR deck has
`time/tlim = 40.0`. Adjust those if the comparison needs a longer saturated
phase.

Important caveat: this is a controlled 2D turbulent MHD-PIC comparison, not a
true 3D turbulent dynamo. Treat "saturation" operationally as a late-time
statistical state of this 2D driven setup.

## Build On Frontier

From a Frontier login node, build from the pulled `PIC` branch:

```bash
export REPO=/path/to/athenak-pic
export BUILD=/path/to/build-frontier-2d-inplane

export PIC_FRONTIER_PROFILE=frontier_minimum_supported
source "$REPO/tst/publication/frontier_control_plane/frontier_pic_environment.sh"
export OMP_NUM_THREADS=1

cmake -S "$REPO" -B "$BUILD" \
  -DCMAKE_BUILD_TYPE=Release \
  -DAthena_SINGLE_PRECISION=OFF \
  -DAthena_ENABLE_MPI=ON \
  -DKokkos_ENABLE_HIP=ON \
  -DKokkos_ARCH_ZEN3=ON \
  -DKokkos_ARCH_AMD_GFX90A=ON \
  -DCMAKE_CXX_COMPILER=/opt/cray/pe/craype/2.7.33/bin/CC \
  -DCMAKE_CXX_FLAGS=-I/opt/rocm-6.2.4/include \
  '-DCMAKE_EXE_LINKER_FLAGS=-L/opt/rocm-6.2.4/lib -lamdhip64' \
  -DPROBLEM=turb

cmake --build "$BUILD" --parallel 32
```

The executable is `$BUILD/src/athena`.

## Run Commands

Use one Frontier node with 8 MPI ranks and one GPU per rank.

MHD-only saturation run:

```bash
export REPO=/path/to/athenak-pic
export BUILD=/path/to/build-frontier-2d-inplane
export RUN=/path/to/run/mhd-saturation

srun -N1 -n8 --ntasks-per-node=8 --cpus-per-task=1 \
  --gpus-per-task=1 --gpu-bind=closest --kill-on-bad-exit=1 \
  "$BUILD/src/athena" \
  -i "$REPO/inputs/particles/pic_turbulent_dynamo_512_2d_mhd_gamma1p00001.athinput" \
  -d "$RUN" \
  -t 01:55:00 time/tlim=25.0 time/nlim=-1
```

CRs from the initial condition:

```bash
export RUN=/path/to/run/cr-from-ic

srun -N1 -n8 --ntasks-per-node=8 --cpus-per-task=1 \
  --gpus-per-task=1 --gpu-bind=closest --kill-on-bad-exit=1 \
  "$BUILD/src/athena" \
  -i "$REPO/inputs/particles/pic_turbulent_dynamo_512_2d_cr_fid_gamma1p00001.athinput" \
  -d "$RUN" \
  -t 01:55:00 time/tlim=40.0 time/nlim=-1
```

CRs injected into a late MHD state:

```bash
export MHD_RUN=/path/to/run/mhd-saturation
export RUN=/path/to/run/cr-injected-at-saturation
export RST=$(find "$MHD_RUN/rst" -maxdepth 1 -name '*.rst' | sort | tail -1)

srun -N1 -n8 --ntasks-per-node=8 --cpus-per-task=1 \
  --gpus-per-task=1 --gpu-bind=closest --kill-on-bad-exit=1 \
  "$BUILD/src/athena" \
  -r "$RST" \
  -i "$REPO/inputs/particles/pic_turbulent_dynamo_512_2d_cr_fid_gamma1p00001.athinput" \
  -d "$RUN" \
  -t 01:55:00 time/tlim=40.0 time/nlim=-1
```

For longer walltime-limited campaigns, restart from the newest file in
`$RUN/rst`. The `-t 01:55:00` option prints a restart before the 2 hour queue
limit.

## Outputs To Compare

MHD outputs:

- History: `*.hst`.
- Cell fields: `*.mhd_w_bcc.*.bin`, every `dt = 0.5`.
- Divergence monitor: `*.mhd_divb.*.bin`, every `dt = 0.5`.
- Restarts: `rst/*.rst`, every `dt = 5.0`.

CR outputs:

- Particle VTK: `*.prtcl_all.*.pvtk`, every `dt = 5.0`.
- Gridded CR diagnostics every `dt = 1.0`:
  `prtcl_ecr`, `prtcl_pcr`, `prtcl_wcr`,
  `prtcl_ekin_flux_x`, `prtcl_ekin_flux_y`, `prtcl_ekin_flux_z`,
  `prtcl_pcr_aniso`.

Good comparison quantities:

- Time histories of kinetic, thermal, magnetic, and CR energy.
- RMS Mach number. With `cs = 1`, use the RMS velocity directly.
- `B_rms`, component energies, and growth of `Bz` from initially zero.
- 2D magnetic and kinetic power spectra.
- Slices of density, pressure, velocity, `Bx`, `By`, `Bz`, `|B|`, and CR energy
  density.
- CR pressure, CR energy flux components, pressure anisotropy, and particle
  pitch-angle/gyro-radius distributions if using particle dumps.
- `max(abs(divB))` and RMS `divB`.

For code-to-code comparison, do not require cell-by-cell agreement unless the
forcing phases and random streams are intentionally matched. The useful targets
are statistical agreement, conservation behavior, spectra, morphology, and
transport diagnostics.

## Validation Artifacts

Existing validation checks on ORNL filesystems:

- Host 2D in-plane IC smoke:
  `/lustre/orion/ast207/proj-shared/dfielding/PIC/turb_box/validation/2d-inplane-smoke-host-20260706T140247/summary.json`
- Frontier GPU smoke:
  `/lustre/orion/ast207/proj-shared/dfielding/PIC/turb_box/validation/2d-frontier-gpu-smoke-20260706T141713/cr-from-ic/gpu_smoke_summary.json`

Expected IC sanity checks from those runs:

- `B_rms = 1.000000047e-3`
- `mean(Bx)`, `mean(By)`, `mean(Bz)` are near zero.
- Initial `Bz` is zero to roundoff.
- `max(abs(divB))` is at roundoff scale.
- The GPU smoke produced `prtcl_ecr`, confirming the gridded CR derived output
  path works.

ORNL-local helper scripts for the current campaign live outside the repo at:

```text
/lustre/orion/ast207/proj-shared/dfielding/PIC/turb_box/scripts/
```

The main helper is `frontier_submit_2d_inplane_campaign.sh`. It submits the
build, both CR variants, the MHD precursor, and the dependent analysis job.
