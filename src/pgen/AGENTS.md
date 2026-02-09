# AGENTS.md

## Purpose
This directory hosts AthenaK problem generators. These are responsible for
initial conditions and for enrolling optional user callbacks (BCs, source
terms, refinement criteria, history outputs, and finalization hooks).
It includes both built-in regression tests (`tests/`) and a large set of
compile-time user problems.
See `../../AGENTS.md` for repository-wide conventions and workflow.

---

## Key Entry Points

### Core types and selection
- `pgen.hpp`: `ProblemGenerator` definition and callback pointers:
  `user_bcs_func`, `user_srcs_func`, `user_ref_func`, `user_hist_func`,
  `pgen_final_func`, plus built-in test declarations.
- `pgen.cpp`: `ProblemGenerator` constructors, restart loading, and selection
  logic for built-in tests vs `UserProblem`.

### Selection flow (verified in `pgen.cpp`)
- New runs:
  - Sets `user_bcs=true` if any mesh face uses `BoundaryFlag::user`.
  - Reads `problem/user_srcs` and `problem/user_hist`.
  - If `USER_PROBLEM_ENABLED`, calls `UserProblem(pin,false)`.
  - Otherwise dispatches on `problem/pgen_name` to a built-in test.
  - After initialization, enforces only requested BC/src/history callbacks were
    enrolled (no checks for `user_ref_func` or `pgen_final_func`).
- Restart runs:
  - Reads restart data into physics arrays.
  - Re-calls `UserProblem(pin,true)` or the built-in test with `restart=true`
    to re-enroll callbacks and re-init any derived state.
  - Performs the same callback checks as the new-run path.

### Built-in test suite (`src/pgen/tests`)
These are the only problems selectable with `problem/pgen_name` when
`USER_PROBLEM_ENABLED` is not set:
- `advection` -> `tests/advection.cpp` (`ProblemGenerator::Advection`)
- `cpaw` -> `tests/cpaw.cpp` (`ProblemGenerator::AlfvenWave`)
- `gr_bondi` -> `tests/gr_bondi.cpp` (`ProblemGenerator::BondiAccretion`)
- `tetrad` -> `tests/rad_check_tetrad.cpp` (`ProblemGenerator::CheckOrthonormalTetrad`)
- `hohlraum` -> `tests/rad_hohlraum.cpp` (`ProblemGenerator::Hohlraum`)
- `linear_wave` -> `tests/linear_wave.cpp` (`ProblemGenerator::LinearWave`)
- `implode` -> `tests/lw_implode.cpp` (`ProblemGenerator::LWImplode`)
- `gr_monopole` -> `tests/gr_monopole.cpp` (`ProblemGenerator::Monopole`)
- `orszag_tang` -> `tests/orszag_tang.cpp` (`ProblemGenerator::OrszagTang`)
- `pic_parallel_shock` -> `tests/pic_parallel_shock.cpp`
  (`ProblemGenerator::PICParallelShock`)
- `rad_linear_wave` -> `tests/rad_linear_wave.cpp`
  (`ProblemGenerator::RadiationLinearWave`)
- `shock_tube` -> `tests/shock_tube.cpp` (`ProblemGenerator::ShockTube`)
- `z4c_linear_wave` -> `tests/z4c_linear_wave.cpp`
  (`ProblemGenerator::Z4cLinearWave`)
- `spherical_collapse` -> `tests/collapse.cpp`
  (`ProblemGenerator::SphericalCollapse`)
- `diffusion` -> `tests/diffusion.cpp` (`ProblemGenerator::Diffusion`)

---

## Problem Generators in This Directory (compile-time `UserProblem` candidates)

### General hydro/MHD tests
- `blast.cpp`: Spherical blast wave; initializes hydro and/or MHD (neutral vs ion
  fluids) and supports coordinate warps for DynGRMHD.
- `current_sheet.cpp`: Double Harris current sheet with optional perturbations.
- `field_loop.cpp`: Field-loop and density-cylinder advection tests with shearing-box
  support and multiple loop orientations.
- `kh.cpp`: Kelvin-Helmholtz instability with multiple profiles; requires
  `nscalars > 0` for a tracer.
- `mri2d.cpp`: 2D shearing-sheet MRI; hydro or MHD; enforces shearing-box sources.
- `mri3d.cpp`: 3D shearing-box MRI (MHD only); seeds pressure perturbations.
- `rt.cpp`: Rayleigh-Taylor instability (2D/3D); uses `const_accel_val` gravity and
  optional smooth interface.
- `shwave.cpp`: Linear shearing-wave tests in a shearing box; hydro (ipert=1/2)
  or MHD (ipert=3).
- `shock_cloud.cpp`: Planar shock impacting a spherical cloud; sets inflow BCs
  via Rankine-Hugoniot relations.
- `shu_osher.cpp`: Shu-Osher shock tube (hydro only).
- `slotted_cyl.cpp`: Slotted-cylinder passive scalar advection (hydro only).

### Turbulence / mass removal / thermal instability
- `turb.cpp`: Turbulent box ICs for hydro, MHD, and ion-neutral two-fluid;
  enrolls `TurbulentHistory`.
- `turbulent_box.cpp`: Uniform hydro box; enrolls `user_srcs_func` (stub) and
  `TurbulentHistory`. (Header comment is stale.)
- `thermal_instability.cpp`: Uniform hydro ICs with `user_srcs_func` (stub) and
  `TurbulentHistory`. (Header comment is stale.)
- `mass_removal_test.cpp`: Uniform hydro ICs with `UserSource` removing mass
  within a radius and adding cooling.
- `turb_amr_test.cpp`: Uniform box for AMR/SMR turbulence driving; defines a
  standalone `RefinementCondition` helper but does not enroll `user_ref_func`.
- `turb_timed_amr.cpp`: Turbulence ICs plus time-triggered refinement
  (`user_ref_func`).
- `turb_timed_amr_bstd.cpp`: Turbulence ICs plus AMR based on B-field
  std/mean threshold (`user_ref_func`).

### Hydrostatic / CGM setups
- `hydrostatic_1d.cpp`: Stratified 1D hydrostatic atmosphere with constant
  gravity; enrolls `UserSource` and `Static` BCs. (Header comment is stale.)
- `hydrostatic_3d.cpp`: 3D hydrostatic atmosphere in a spherical potential;
  enrolls `UserSource` and `Static` BCs. (Header comment is stale.)
- `cgm_static.cpp`: Hydrostatic CGM in a composite potential; gravity source
  and static BCs.
- `cgm_cooling_flow.cpp`: Cooling-flow CGM initialized from profile files; adds
  gravity/mass-loss sources, custom BCs, and `pgen_final_func=FreeProfile`.
- `cgm_cooling_flow_magnetized.cpp`: Magnetized cooling-flow variant with an
  initial B field; same profile/BC/source structure.
- `cgm_cooling_flow_zoom.cpp`: Cooling-flow plus disk-profile initialization
  without AMR.
- `cgm_cooling_flow_amr.cpp`: Cooling-flow plus disk-profile initialization
  with refinement driven by density gradients.
- `cgm_cooling_flow_amr_metals.cpp`: Adds metallicity and SN injection
  (uses `sn_scheduler` and particles) plus AMR refinement.

### Radiation tests
- `rad_beam.cpp`: GR beam test; sets a beam mask and zero-intensity BCs.
- `rad_diffusion.cpp`: Static/dynamic diffusion test; initializes hydro and
  angular intensities.
- `rad_relax.cpp`: Thermal relaxation test; uniform hydro and radiation fields.
- `rad_shadow.cpp`: Shadowing test; requires geodesic mesh (`nlevel=2`,
  `rotate_geo=false`) and sets inflow radiation BCs.
- `rad_snake.cpp`: "Snake" beam test; overrides metric/tetrad and recomputes
  `na` using existing angular mesh data (does not change neighbor lists).

### Particles
- `part_random.cpp`: Random particle positions/velocities; sets particle
  timestep (`ppart->dtnew`).

### Two-fluid / ion-neutral
- `cshock.cpp`: C-shock profile in ion-neutral MHD; integrates ODE on host and
  sets inflow states.
- `twofluid.cpp`: Simple two-fluid (hydro + MHD) initialization with B field.

### GR / numerical relativity
- `gr_torus.cpp`: GR equilibrium torus (Fishbone-Moncrief or Chakrabarti) in
  Kerr-Schild coordinates; enrolls `NoInflowTorus` BCs and `TorusFluxes`
  history; populates `spherical_grids`.
- `dyngr_tov.cpp`: Dynamical-GR TOV star; selects EOS policy; enrolls
  `VacuumBC` and `TOVHistory`.
- `elliptica_bns.cpp`: Binary neutron star initial data via Elliptica; enrolls
  `EllipticaHistory`.
- `lorene_bns.cpp`: Binary neutron star initial data via LORENE; enrolls
  `BNSHistory`.
- `sgrid_bns.cpp`: Binary neutron star initial data via SGRID; enrolls history
  and AMR refinement.
- `z4c_one_puncture.cpp`: Single puncture initial data; calls Z4c ADM
  conversions and refinement.
- `z4c_two_puncture.cpp`: Two punctures via TwoPunctures; sets Z4c state and
  refinement.
- `z4c_spectre_bbh.cpp`: SpECTRE BBH initial data loader; sets Z4c state and
  refinement.

---

## Implementation Notes
- Kokkos-first initialization: ICs are set with `par_for` on device. Keep new
  ICs device-resident unless external readers require host loops.
- Orszag-Tang scan controls: `tests/orszag_tang.cpp` now accepts optional
  `<problem>` knobs `ot_B0`, `ot_d0`, `ot_p0`, `ot_v0`, and `ot_mach`.
  When `ot_mach > 0`, velocity amplitude is derived from
  `Mach * sqrt(gamma * p0 / d0)`; defaults preserve historical ICs.
- Callback enrollment: if `problem/user_srcs` or `problem/user_hist` is set,
  `UserProblem` must set `user_srcs_func` / `user_hist_func`, or the constructor
  will abort. User BCs must be set when any mesh face uses `BoundaryFlag::user`.
- MHD face fields: when filling `b0`, remember to set the extra face at block
  edges (`i==ie`, `j==je`, `k==ke`).
- Restart behavior: `UserProblem(..., true)` is called after restart read; most
  pgens early-return to avoid reinitializing state but still re-enroll callbacks.
- Coupled PIC restart section: when `couple_moments_to_mhd=true` and
  `deposit_moments=true`, `pgen.cpp` restores particle real/int arrays, moments,
  and optional edge-current state from the PR3a section for both single-file and
  `single_file_per_rank` restarts. The per-rank path uses restart metadata
  (`rank_eachmb`, `gids_eachrank`, `nmb_eachrank`) to remap each local
  MeshBlock's source file/offset before loading.
- `tests/pic_parallel_shock.cpp` is the physics-benchmark pgen for Section 5.4
  style parallel-shock studies (reflecting wall, `B0 || x`, eta injection,
  conservative gas subtraction, and curvature AMR via `g_rho`/`g_P`).
  Keep this benchmark path separate from Orszag-Tang smoke/pipeline workflows.

---

## Cautions
- Several files have stale header comments (notably `hydrostatic_1d.cpp`,
  `hydrostatic_3d.cpp`, `thermal_instability.cpp`, `turbulent_box.cpp`,
  `twofluid.cpp`). Trust the code, not the banner, when documenting behavior.
- Many GR/BNS pgens depend on optional external libraries (Elliptica, LORENE,
  SGRID, TwoPunctures, SpECTRE). Keep their includes and data-loading paths
  intact.
- Coupled PIC restart load now hard-validates section marker/version, array
  dimensions, and per-MeshBlock particle-count tables before copying data into
  Kokkos views; any mismatch aborts with a fatal error.
