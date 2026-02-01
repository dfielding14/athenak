# AGENTS.md

## Purpose
This directory implements the general-relativistic radiation module. It provides
angular transport on a geodesic grid, explicit RK updates for specific intensity,
implicit radiation-matter coupling, and tetrad construction in Cartesian
Kerr-Schild coordinates. Radiation can run standalone or coupled to Hydro/MHD.
See `../../AGENTS.md` for repository-wide conventions and workflow.

---

## Key Files and Responsibilities

### Core module
- `radiation.hpp`: `radiation::Radiation` class, task IDs, and key data members.
- `radiation.cpp`: constructor/destructor, input parsing, array allocation, and
  geodesic-grid + tetrad setup.

### Time integration and tasks
- `radiation_tasks.cpp`: task-list assembly and boundary/AMR wrappers.
- `radiation_update.cpp`: explicit RK update for `i0` (spatial + angular fluxes).
- `radiation_newdt.cpp`: radiation timestep estimate (computed once at init).

### Fluxes, sources, and geometry
- `radiation_fluxes.cpp`: spatial and angular flux calculations on the geodesic grid.
- `radiation_source.cpp`: implicit radiation-matter coupling and Compton update.
- `radiation_tetrad.cpp` / `radiation_tetrad.hpp`: tetrad construction in CKS.
- `radiation_opacities.hpp`: opacity model used by `radiation_source.cpp`.

---

## Activation and Preconditions
- Radiation **requires** GR coordinates: constructor exits if
  `pcoord->is_general_relativistic` is false.
- AMR is **not** supported: constructor exits if `pmesh->adaptive` is true.
- Two-fluid runs are **not** supported: exits if both `<hydro>` and `<mhd>` exist.
- Angular grid is provided by `GeodesicGrid` constructed from `<radiation>` params.

---

## Configuration Inputs (from `<radiation>` block)
Parsed in `radiation.cpp`:
- `nlevel` (required): geodesic grid refinement level.
- `rotate_geo` (optional, default `true`): rotate geodesic grid.
- `angular_fluxes` (optional, default `true`): enable angular advection.
- `n_0_floor` (optional, default `0.1`): floor for `n_0` in excision handling.
- `rad_source` (optional, default `true` if Hydro/MHD enabled): enable coupling.
- `fixed_fluid` (optional, default `false`): skip Hydro/MHD time integration tasks.
- `affect_fluid` (optional, default `true`): feedback radiation on fluid in coupling.
- `kappa_s` (required if `rad_source`): scattering opacity.
- `power_opacity` (optional, default `false`): use power-law opacities.
- `kappa_a`, `kappa_p` (required if `rad_source` and `power_opacity=false`).
- `compton` (optional, default `false`): Compton update (requires `<units>`).
- `arad` (required if `rad_source` and `<units>` not enabled).
- `beam_source` (optional, default `false`): enable beam source term.
- `reconstruct` (optional, default `plm`): `dc`, `plm`, `ppm4`, `ppmx`, `wenoz`.
  High-order methods require `nghost >= 3`.

Notes:
- If no Hydro/MHD, `rad_source` is forced to `false` in the constructor.
- `beam_source` is handled by `SourceTerms::BeamSource` and requires `beam_mask`.

---

## Data Layout (major arrays)
All arrays are Kokkos views with pack-major layout.

### Angular grid and tetrad
- `prgeo`: `GeodesicGrid` (angles, neighbors, arc lengths, unit flux vectors).
- `nh_c[nangles][4]`: normal-frame components at angle centers.
- `nh_f[nangles][6][4]`: normal-frame components at angle edges.
- `tet_c[nmb][4][4][k][j][i]`: tetrad components at cell centers.
- `tetcov_c[nmb][4][4][k][j][i]`: covariant tetrad components at cell centers.
- `tet_d1_x1f`, `tet_d2_x2f`, `tet_d3_x3f`: tetrad subsets at faces.
- `na[nmb][nangles][k][j][i][6]`: angular advection coefficients (if enabled).
- `norm_to_tet[nmb][4][4][k][j][i]`: normal-to-tetrad transform.

### Radiation state
- `i0[nmb][nangles][k][j][i]`: specific intensity (conserved).
- `i1`: intermediate register for RK stages.
- `iflx`: face-centered spatial fluxes (`x1f/x2f/x3f`).
- `divfa`: angular flux divergence (only if `angular_fluxes`).
- `coarse_i0`: coarse-grid buffer for SMR/AMR (when `multilevel`).
- `beam_mask`: boolean mask for beam source (only if `beam_source`).

### Boundary handling
- `pbval_i`: `MeshBoundaryValuesCC` for `i0` (and fluxes if multilevel).

---

## Task List Flow
Assembled in `radiation_tasks.cpp` with three modes.

### Radiation + MHD (when `<mhd>` and `fixed_fluid=false`)
**before_stagen**
- `Radiation::InitRecv` and `MHD::InitRecv`

**stagen**
- `Radiation::CopyCons` -> `Radiation::CalculateFluxes` -> `SendFlux` -> `RecvFlux`
- `Radiation::RKUpdate`
- `MHD::Fluxes` -> `SendFlux` -> `RecvFlux` -> `RKUpdate`
- `MHD::CornerE` -> `SendE` -> `RecvE` -> `CT`
- `Radiation::AddRadiationSourceTerm`
- `Radiation::RestrictI` -> `SendI` -> `RecvI`
- `MHD::RestrictU` -> `SendU` -> `RecvU`
- `MHD::RestrictB` -> `SendB` -> `RecvB`
- `Radiation::ApplyPhysicalBCs`
- `Radiation::Prolongate` -> `MHD::Prolongate` -> `MHD::ConToPrim`

**after_stagen**
- `Radiation::ClearSend` -> `Radiation::ClearRecv`
- `MHD::ClearSend` -> `MHD::ClearRecv`

### Radiation + Hydro (when `<hydro>` and `fixed_fluid=false`)
Same structure as MHD case, but with Hydro tasks instead of MHD/CT.

### Radiation transport only
**before_stagen**
- `Radiation::InitRecv`

**stagen**
- `CopyCons` -> `CalculateFluxes` -> `SendFlux` -> `RecvFlux` -> `RKUpdate`
- `AddRadiationSourceTerm` (no-op unless `rad_source` true)
- `RestrictI` -> `SendI` -> `RecvI`
- `ApplyPhysicalBCs` -> `Prolongate`

**after_stagen**
- `ClearSend` -> `ClearRecv`

Notes:
- `Radiation::NewTimeStep` is **not** scheduled in task lists; it is called
  from `Driver::Initialize` once at startup.

---

## Core Algorithms

### Spatial fluxes (`radiation_fluxes.cpp`)
- Upwinded using `n^1`, `n^2`, `n^3` computed from tetrads and geodesic normals.
- Reconstruction applies to the primitive `n_0 I` (using `tet_c` for conversion).

### Angular fluxes (optional)
- Uses geodesic neighbor graph (`num_neighbors`, `ind_neighbors`).
- Edge lengths from `arc_lengths` and weights from `solid_angles`.
- Uses `na` computed in `radiation_tetrad.cpp` and `unit_flux` from `GeodesicGrid`.

### Explicit update (`radiation_update.cpp`)
- RK update for `i0` with spatial and angular flux divergence.
- Enforces non-negativity on `n_0 I` and applies excision if enabled.
- Beam source term is added via `SourceTerms::BeamSource` when `beam_source`.

### Implicit coupling (`radiation_source.cpp`)
- Solves for gas temperature via a quartic (`FourthPolyRoot`).
- Updates intensities in the comoving frame and optionally feeds back to fluid
  (`affect_fluid`).
- Compton update is optional and requires `<units>`.
- Uses `OpacityFunction` for constant or power-law opacities.

### Tetrads (`radiation_tetrad.cpp`)
- Tetrads are computed in **Cartesian Kerr-Schild** using
  `ComputeMetricAndInverse`, `ComputeMetricDerivatives`, and `ComputeTetrad`.
- `nh_c` and `nh_f` are derived from the geodesic grid positions.
- `na` is computed only when `angular_fluxes` is enabled.

### Timestep (`radiation_newdt.cpp`)
- Computes min spatial `dx` and angular `dangle / |na/n0|` (if enabled).
- Applied once at initialization; later dt control is handled by other modules.

---

## Boundary and AMR Handling
- Boundary communication uses `MeshBoundaryValuesCC` (`pbval_i`).
- `SendFlux` / `RecvFlux` are only active when `pmesh->multilevel` is true.
- `RestrictI` / `Prolongate` operate on `i0` and `coarse_i0` for SMR.

---

## Cautions
- Radiation currently assumes GR in Cartesian Kerr-Schild; it does not run in
  non-GR coordinate systems.
- AMR is disabled (`pmesh->adaptive`); SMR buffers exist but AMR is not supported.
- `fixed_fluid=true` skips Hydro/MHD time integration tasks, but
  `AddRadiationSourceTerm` can still update fluid if `affect_fluid=true`.
- `beam_source` requires `beam_mask` to be initialized by the problem generator.

---

## Extension Points
- Add new opacity laws in `radiation_opacities.hpp` and thread through
  `radiation_source.cpp`.
- Add new angular transport schemes by extending `na` computation or the angular
  flux loop in `radiation_fluxes.cpp`.
- Add new radiation diagnostics by integrating over `solid_angles` and `i0`.
