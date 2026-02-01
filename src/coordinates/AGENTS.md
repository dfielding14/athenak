# AGENTS.md

## Purpose
This directory provides coordinate-system utilities, GR metric helpers, ADM storage,
and black-hole excision masks used across hydro, MHD, radiation, and dynamical GR.
See `../../AGENTS.md` for repository-wide conventions and workflow.

---

## Key Files and Responsibilities

### Coordinates core
- `coordinates.hpp` / `coordinates.cpp`: `Coordinates` class, SR/GR flags, GR source
  terms (`CoordSrcTerms`), and excision mask management.

### GR metric helpers
- `cartesian_ks.hpp`: Cartesian Kerr-Schild metric utilities:
  `ComputeMetricAndInverse`, `ComputeMetricDerivatives`, and
  `ComputeADMDecomposition` (lapse, shift, spatial metric, extrinsic curvature).

### ADM storage and helpers
- `adm.hpp` / `adm.cpp`: `adm::ADM` class (lapse, shift, g_ij, K_ij, psi4) plus
  helpers for spatial determinant/inverse and face-centered metric averages
  (`Face1Metric/Face2Metric/Face3Metric`). When Z4c is present, lapse/shift are
  shallow slices into Z4c storage.

### Excision masks
- `excision.cpp`: `SetExcisionMasks` (fixed CKS radius) and `UpdateExcisionMasks`
  (lapse-based scheme) to build `excision_floor` and `excision_flux` masks.

### Cell location utilities
- `cell_locations.hpp`: `LeftEdgeX`, `CellCenterX`, and `CellCenterIndex` for
  uniform Cartesian grids, used by mesh and coordinate logic.

---

## Coordinate Flow (high level)
1. `Coordinates` constructor determines SR/GR mode. If `<adm>` or `<z4c>` exists
   alongside hydro/MHD, it sets `is_dynamical_relativistic` and disables SR/GR
   coordinate flags.
2. For GR, it reads the `<coord>` block and sets up excision parameters. If
   excision is enabled with a fixed scheme, it allocates and computes masks once.
3. Hydro/MHD task lists call `CoordSrcTerms` for GR geometric source terms; MHD
   uses the magnetic-field overload.
4. Excision masks are consumed by EOS/primitive solvers, FOFC, and radiation
   sources. Lapse-based excision is refreshed via Z4c task hooks.

---

## Configuration Inputs (coord block)
Typical fields parsed in `Coordinates` when GR is active:
- `special_rel` (bool), `general_rel` (bool)
- `minkowski` (bool) to force flat metric
- `a` (BH spin), `excise` (bool), `dexcise`, `pexcise`
- `flux_excise_r` (override for FOFC region)
- `excision_scheme` (`fixed` or `lapse`) and `excise_lapse` threshold

Defaults for excision radius and flux excision are set to the Kerr horizon when
radiation is enabled.

---

## Extension Points and Cautions
- Adding a new coordinate system or metric requires new metric helpers and a
  selection path in `Coordinates` (and likely EOS/primitive solver support).
- If you change ADM variable layout, update `ADM_names` and any consumers in
  `dyn_grmhd` and `z4c` that assume current indices.
- Excision masks assume Cartesian Kerr-Schild geometry; if you add a new metric,
  revisit the excision logic and its radius definition.
