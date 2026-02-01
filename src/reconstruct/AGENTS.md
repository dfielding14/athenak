# AGENTS.md

## Purpose
This directory provides header-only reconstruction routines used by Hydro, MHD,
DynGRMHD, and Radiation to build left/right interface states for fluxes. All
routines are `KOKKOS_INLINE_FUNCTION`s; the directional wrappers operate on
`DvceArray5D<Real>` inputs with `ScrArray2D<Real>` scratch outputs, while the
scalar helpers (`PLM`, `PPM4`, `PPMX`, `WENOZ`) operate on scalar values.
See `../../AGENTS.md` for repository-wide conventions and workflow.

---

## Key Files and Responsibilities

### `dc.hpp`
- Piecewise-constant (donor cell) reconstruction.
- `DonorCellX1/X2/X3` fill `ql`/`qr` with cell-centered values for each variable.
- Commented call ranges specify how to cover both L/R states on the active zone.

### `plm.hpp`
- Piecewise-linear reconstruction with a slope limiter for uniform spacing.
- `PLM(q_im1, q_i, q_ip1, ql_ip1, qr_i)` computes a limited slope:
  - Uses the harmonic-mean form `dqm = (dql*dqr)/(dql+dqr)` when slopes agree.
  - Sets slope to zero when slopes change sign.
- `PiecewiseLinearX1/X2/X3` apply `PLM` for all variables.

### `ppm.hpp`
- Piecewise-parabolic reconstruction for uniform, Cartesian-like grids.
- `PPM4`: Colella–Woodward (CW) limiters.
- `PPMX`: Colella–Sekora (CS) extremum-preserving limiters.
  - CS extensions from McCorquodale et al. are explicitly not included (per file
    header comment).
- `PiecewiseParabolicX1/X2/X3` wrappers:
  - `extremum_preserving=true` selects `PPMX`, otherwise `PPM4`.
  - `apply_floors=true` clamps `IDN`/`IEN` only in the `PPMX` branch; the `PPM4`
    path ignores `apply_floors` entirely.
  - Energy floor uses `efloor = pfloor/(gamma-1)`; code comments note ideal-gas
    assumptions only.

### `wenoz.hpp`
- Fifth-order WENO-Z reconstruction for uniform, Cartesian-like grids.
- `WENOZ`: uses WENO-Z+ style weights (`tau_5 = |beta0-beta2|`, `eps=1e-42`) and
  fixed ideal weights `{0.1, 0.6, 0.3}`.
- `WENOZX1/X2/X3` wrappers mirror the PPM wrappers and apply optional floors to
  `IDN`/`IEN` only.

---

## Interfaces and Calling Conventions

- All directional wrappers take:
  - `TeamMember_t` for Kokkos team-level parallelism.
  - `DvceArray5D<Real>` input `q` with `nvar = q.extent_int(1)`.
  - `ScrArray2D<Real>` scratch buffers for `ql`/`qr` outputs.
- Indexing semantics (from file comments and usage):
  - X1 wrappers write `ql(n,i+1)` and `qr(n,i)` for interface states.
  - X2 wrappers write `ql_jp1(n,i)` and `qr_j(n,i)` for j-interfaces.
  - X3 wrappers write `ql_kp1(n,i)` and `qr_k(n,i)` for k-interfaces.
- Stencils:
  - DC: single cell (`q_i`).
  - PLM: 3-point stencil (`q_{i-1}, q_i, q_{i+1}`).
  - PPM/WENOZ: 5-point stencil (`q_{i-2}..q_{i+2}`).
- The wrappers are dimension-agnostic; callers supply `(m,k,j, il, iu)` and the
  appropriate index permutation for the direction being reconstructed.

---

## Usage in the Codebase

### Hydro (`src/hydro/hydro_fluxes.cpp`)
- Chooses the reconstruction method from `<hydro>/reconstruct`.
- Calls `DonorCellX*`, `PiecewiseLinearX*`, `PiecewiseParabolicX*`, or `WENOZX*`.
- For WENOZ, `apply_floors=true` on primitive variables (`w0`).
- For PPM, `apply_floors=true` is passed, but only `PPMX` enforces floors; `PPM4`
  does not apply floors.
- `recon_method == ppmx` sets `extremum_preserving=true` for the PPM wrapper.

### MHD (`src/mhd/mhd_fluxes.cpp`)
- Uses the same wrappers for both primitive variables (`w0`) and cell-centered
  magnetic fields (`bcc0`).
- PPM/WENOZ floors are applied only to `w0` (`apply_floors=true`) and **not** to
  magnetic fields (`apply_floors=false`).

### DynGRMHD (`src/dyn_grmhd/dyn_grmhd_fluxes.cpp`)
- Uses the same wrappers for `w0` and cell-centered magnetic fields (`bcc0`).
- PPM/WENOZ floors are **disabled** (`apply_floors=false`) for all fields; code
  comments note that EOS floors are not used in this module.

### Radiation (`src/radiation/radiation_fluxes.cpp`)
- Uses `PLM`, `PPM4`, `PPMX`, and `WENOZ` directly (not the X1/X2/X3 wrappers).
- Performs directional upwinding based on the sign of the streaming direction
  component `n^i`.
- Reconstruction method is selected from `<radiation>/reconstruct`.

---

## Configuration and Constraints

- Reconstruction methods are enumerated in `ReconstructionMethod` in
  `src/athena.hpp`: `dc`, `plm`, `ppm4`, `ppmx`, `wenoz`.
- Selection occurs in module constructors:
  - Hydro: `src/hydro/hydro.cpp`
  - MHD: `src/mhd/mhd.cpp`
  - Radiation: `src/radiation/radiation.cpp`
- Ghost-zone requirements enforced by modules:
  - `ppm4`/`ppmx`/`wenoz` require `nghost >= 3` in Hydro/MHD/Radiation.
  - `plm` with FOFC requires `nghost >= 3` in Hydro/MHD.
  - `ppm4`/`ppmx`/`wenoz` with FOFC require `nghost >= 4` in Hydro/MHD.
- Uniform-spacing assumptions: PLM, PPM, and WENOZ are documented as valid for
  Cartesian-like coordinates with uniform mesh spacing.

---

## Extension Points

- Add a new reconstruction method by:
  - Implementing the kernel in `src/reconstruct/`.
  - Extending `ReconstructionMethod` in `src/athena.hpp`.
  - Wiring selection logic in `hydro.cpp`, `mhd.cpp`, and `radiation.cpp`.
  - Calling the new method in the corresponding `*_fluxes.cpp` routines.
- If the new method changes stencil width, update ghost-zone checks and any FOFC
  constraints accordingly.

---

## References (as documented in headers)
- Colella & Woodward (1984), Colella & Sekora (2008) for PPM.
- McCorquodale & Colella (2011) and Peterson & Hammett (2013) referenced in PPM
  comments (not fully implemented).
- Borges et al. (2008) for WENO-Z.
- Acker et al. (2016) for the WENO-Z+ weighting variant used.
- Del Zanna et al. (2007, A.18) for WENO-Z smoothness coefficients.
