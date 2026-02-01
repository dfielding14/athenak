# AGENTS.md

## Purpose
This directory implements the geodesic angular mesh used by radiation and a
spherical interpolation helper used by diagnostics and wave extraction. It
provides neighbor connectivity, solid angles, arc lengths, and optional angular
unit vectors for fluxes on a geodesic sphere.
See `../../AGENTS.md` for repository-wide conventions and workflow.

---

## Key Files and Responsibilities

### Angular mesh
- `geodesic_grid.hpp` / `geodesic_grid.cpp`: `GeodesicGrid` class. Builds the
  geodesic sphere, neighbor graph, solid angles, arc lengths, and (optionally)
  unit vectors along angular edges.

### Spherical interpolation
- `spherical_grid.hpp` / `spherical_grid.cpp`: `SphericalGrid` class (subclass of
  `GeodesicGrid`) that interpolates 3D mesh data onto a fixed-radius sphere.

---

## GeodesicGrid Overview

### Construction parameters
`GeodesicGrid(int nlev, bool rotate, bool fluxes)`
- `nlevel` (`nlev`): geodesic refinement level.
- `rotate_geo` (`rotate`): rotate mesh to an "optimal" orientation.
- `geo_fluxes` (`fluxes`): compute angular unit vectors for edge fluxes.

### Number of angles
- For `nlevel > 0`: `nangles = 5 * 2 * nlevel^2 + 2`.
- For `nlevel == 0`: special 8-angle (octant) mesh for testing.
- For `nlevel < 0`: fatal error.

### Stored arrays (DualArray)
Allocated and filled when `nlevel > 0`:
- `num_neighbors[nangles]`: number of neighbors (5 or 6).
- `ind_neighbors[nangles][6]`: neighbor indices; `INT_MAX` marks missing neighbor.
- `ind_neighbors_edges[nangles][6]`: index of this angle in each neighbor's list;
  `INT_MAX` for missing.
- `solid_angles[nangles]`: face solid angles.
- `arc_lengths[nangles][6]`: edge arc lengths; `FLT_MAX` for missing.
- `cart_pos[nangles][3]`: unit cartesian position of face centers.
- `cart_pos_mid[nangles][6][3]`: unit cartesian midpoints of face edges.
- `polar_pos[nangles][2]`: `(theta, phi)` at face centers.
- `polar_pos_mid[nangles][6][2]`: `(theta, phi)` at face edges.
- `unit_flux[nangles][6][2]`: `(d zeta, d psi)` along edges (only if `geo_fluxes`).

Internal host arrays used to build the mesh:
- `amesh_normals[5][2+nlevel][2+2*nlevel][3]`: normals for 5 patches.
- `ameshp_normals[2][3]`: normals at the poles.
- `amesh_indices[5][2+nlevel][2+2*nlevel]`: 2D-to-1D map for patches.
- `ameshp_indices[2]`: indices for the poles.

### Construction sequence (nlevel > 0)
1. Build one patch of face-center normals from six base points and normalize.
2. Rotate that patch into the other four patches about the z-axis.
3. Fill ghost cells in `amesh_normals` and `amesh_indices` to provide wrap-around
   neighbor connectivity.
4. Create the 2D-to-1D index map and assign pole indices at the end of the list.
5. For each angle:
   - Compute neighbors with `Neighbors`.
   - Compute solid angle and arc lengths with `SolidAngleAndArcLengths`.
6. Build `ind_neighbors_edges` by cross-referencing neighbor lists.
7. Symmetrize shared arc lengths to remove round-off differences.
8. Optional rotation: `OptimalAngles` -> `RotateGrid`.
9. Compute `cart_pos`, `cart_pos_mid`, `polar_pos`, `polar_pos_mid`.
10. If `geo_fluxes` is true, compute `unit_flux` for edges with `UnitFluxDir`,
    then symmetrize across shared edges.
11. Sync all DualArrays from host to device.

### Special case: `nlevel == 0`
- Emits a warning; intended for testing only.
- Requires `rotate_geo == false` and `geo_fluxes == false` (otherwise fatal).
- Initializes `nangles = 8`, `solid_angles`, and `cart_pos` only.
- Neighbor/edge arrays are not populated.

---

## GeodesicGrid Helper Functions
- `GridCartPosition`: cartesian position at face center (unit vector).
- `GridCartPositionMid`: normalized midpoint between two face centers.
- `Neighbors`: returns neighbor list (5 or 6) and handles poles.
- `CircumcenterNormalized`: spherical circumcenter of a triangle on the unit
  sphere.
- `SolidAngleAndArcLengths`: polygon area (solid angle) and edge arc lengths.
- `ArcLength`: great-circle distance between two face centers.
- `OptimalAngles`: brute-force scan over `(zeta, psi)` to pick a rotation that
  maximizes the minimum component magnitude of rotated normals.
- `RotateGrid`: apply rotation to all normals and refill ghost cells.
- `UnitFluxDir`: compute unit direction along an edge in angular coordinates.
- `GreatCircleParam`: helper for great-circle parameterization.

---

## SphericalGrid Overview

### Purpose
`SphericalGrid` is a `GeodesicGrid` specialized for sampling/interpolating mesh
fields on a sphere of fixed radius.

### Construction
`SphericalGrid(MeshBlockPack *ppack, int nlev, Real rad)`:
- Calls `GeodesicGrid(nlev, true, false)` (always rotated; no angular flux vectors).
- Allocates interpolation arrays sized by `nangles` and `2*ng` (ghost zones).

### Interpolation data
- `interp_coord[nangles][3]`: target cartesian coordinates on the sphere.
- `interp_indcs[nangles][4]`: (MeshBlock index, i, j, k) of the nearest cell
  center for each angle.
- `interp_wghts[nangles][2*ng][3]`: Lagrange weights for x/y/z dimensions.
- `interp_vals[nangles][nvars]`: interpolated values (host-synced).

### Interpolation workflow
1. `SetInterpolationCoordinates`:
   - For GR (`is_general_relativistic` or `is_dynamical_relativistic`), converts
     `(theta, phi)` to Cartesian Kerr-Schild coordinates using `bh_spin`.
   - Otherwise uses standard spherical-to-Cartesian mapping.
2. `SetInterpolationIndices`:
   - Finds the MeshBlock containing each target point and computes the nearest
     cell indices.
   - If not on this pack, indices remain `-1`.
3. `SetInterpolationWeights`:
   - Builds Lagrange weights using a `2*ng` stencil in each dimension.
   - If `ii0 == -1`, weights are zeroed.
4. `InterpolateToSphere(nvars, val)`:
   - Recomputes indices/weights if `pmesh->adaptive` is true.
   - Evaluates a `2*ng` x `2*ng` x `2*ng` stencil per angle and variable.
   - If `ii0 == -1`, output is `0` on that rank.
   - Syncs `interp_vals` to host for diagnostics.

---

## Integration Points
- **Radiation** (`radiation/radiation.cpp`):
  - Builds `GeodesicGrid` using `<radiation>/nlevel`, `rotate_geo`,
    `angular_fluxes`.
  - Uses `num_neighbors`, `ind_neighbors`, `arc_lengths`, and `solid_angles` to
    compute angular flux divergence.
  - Uses `cart_pos`, `cart_pos_mid`, and `unit_flux` to construct tetrads and
    angular advection terms.
- **Derived outputs** (`outputs/derived_variables.cpp`):
  - Radiation moments integrate over `solid_angles`.
- **Problem generators / tests**:
  - `SphericalGrid` is used for diagnostics in GR problems (e.g., `gr_torus`,
    `gr_monopole`) and in Z4c wave extraction.
  - Radiation tests reference geodesic neighbor data and `unit_flux`.

---

## Cautions
- `nlevel == 0` only populates `solid_angles` and `cart_pos`; neighbor and edge
  arrays are not initialized.
- `rotate_geo` and `geo_fluxes` are incompatible with `nlevel == 0`.
- `SphericalGrid` interpolation assumes valid ghost zones and uses a `2*ng`
  stencil; ensure ghost zones are populated before interpolation.
- If an angle is not inside any MeshBlock of a pack, interpolation yields `0` for
  that rank; callers may need MPI aggregation.

---

## Extension Points
- Add new angular metrics (e.g., curvature, face areas) by extending
  `GeodesicGrid` construction and syncing new DualArrays.
- Add alternative angular flux schemes by extending `unit_flux` computation or
  by adding new edge-based data products.
- Modify interpolation stencil size or scheme in `SphericalGrid` by changing
  `SetInterpolationWeights` and `InterpolateToSphere`.
