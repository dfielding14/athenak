# AGENTS.md

## Purpose
This directory implements AthenaK's output system. It parses `<output[n]>` blocks,
selects variables, assembles host-side buffers, and writes mesh/particle/diagnostic
files in multiple formats.
See `../../AGENTS.md` for repository-wide conventions and workflow.

---

## Key Entry Points

### Core types
- `outputs.hpp`: `OutputParameters`, `BaseTypeOutput`, and `Outputs` definitions plus
  `var_choice` (the canonical list of selectable output variables).
- `outputs.cpp`: parses `<output[n]>` blocks and registers each output type.
- `basetype_output.cpp`: validates requested variables against enabled physics and
  builds the per-output variable list.
- `derived_variables.cpp`: computes derived variables into `BaseTypeOutput::derived_var`.

### I/O abstraction
- `io_wrapper.hpp` / `io_wrapper.cpp`: `IOWrapper` provides MPI-IO vs stdio wrappers
  used by several output writers.

---

## Output Types (file_type → implementation)

Mesh data:
- `vtk` → `vtk_mesh.cpp` (`MeshVTKOutput`)
- `bin` → `binary.cpp` (`MeshBinaryOutput`)
- `cbin` → `coarsened_binary.cpp` (`CoarsenedBinaryOutput`)

Particles:
- `pvtk` → `vtk_prtcl.cpp` (`ParticleVTKOutput`)
- `trk` → `track_prtcl.cpp` (`TrackedParticleOutput`); payload rows are big-endian
  float32 `(x, y, z, vx, vy, vz)` records written at offsets determined by
  particle tags `0 <= PTAG < nparticles`

Diagnostics:
- `hst` → `history.cpp` (`HistoryOutput`)
- `log` → `eventlog.cpp` (`EventLogOutput`)
- `tab` → `formatted_table.cpp` (`FormattedTableOutput`)
- `pdf` → `pdf.cpp` (`PDFOutput`)

Restart:
- `rst` → `restart.cpp` (`RestartOutput`)

The list above matches the registration in `Outputs::Outputs` (`outputs.cpp`).

---

## Configuration Notes
- `<output[n]>` blocks are discovered by name in `Outputs::Outputs`.
- `file_type` is required; `dt` or `dcycle` is required to schedule outputs.
- `variable` is required for all file types except `hst`, `log`, and `rst`.
  `trk` still parses `variable` in `Outputs::Outputs` (even though it is unused in
  `TrackedParticleOutput`), so include it to avoid parser errors.
- Only one `hst`, `log`, and `rst` block is allowed.
- Output-specific parameters (e.g., `single_file_per_rank`, `coarsen_factor`, PDF
  bin settings) are parsed in `outputs.cpp` and individual writers.
- `pdf` parameters are 1D/2D only: `bin_min`, `bin_max`, `nbin`, and `logscale` always
  apply; `variable_2` plus `bin2_*`, `nbin2`, `logscale2` enable 2D (dimension is
  `nbin2 == 0 ? 1 : 2` in `PDFOutput`).
- `pdf` rejects variable groups (`hydro_u`, `hydro_w`, `mhd_u`, `mhd_w`) in
  `Outputs::Outputs`; only scalar names from `var_choice` are accepted.
- `pdf` supports `mass_weighted` (cell volume times density from hydro/MHD/Z4c).
- `cbin` supports `compute_moments`; when true, each variable outputs 4 moments
  (mean, mean^2, mean^3, mean^4) over the coarsened cell volume.

---

## Extension Points

### PIC diagnostic derived outputs (current)
- Particle moment outputs include:
  - `prtcl_rho`, `prtcl_jx`, `prtcl_jy`, `prtcl_jz`
  - `prtcl_dpxdt`, `prtcl_dpydt`, `prtcl_dpzdt`, `prtcl_dedt`
  - `prtcl_ebdot`
  - `prtcl_jx_edge`, `prtcl_jy_edge`, `prtcl_jz_edge`
- The `prtcl_j*_edge` outputs are cell-centered projections of the edge-current
  arrays used by `MHD::EFieldSrc`; they are diagnostics, not raw staggered dumps.
- All of the above require `<particles>/deposit_moments=true`; validation is
  enforced in `BaseTypeOutput` before output object construction.

### Add a new output variable
- Update `NOUTPUT_CHOICES` and `var_choice` in `outputs.hpp`.
- Map the new variable to data in `BaseTypeOutput` (constructor in
  `basetype_output.cpp`).
- Add derived computation in `derived_variables.cpp` if needed.

### Add a new output format
- Implement a `BaseTypeOutput` subclass with `WriteOutputFile` (and `LoadOutputData`
  if needed).
- Register the new type in `Outputs::Outputs` in `outputs.cpp`.

---

## Observed Output Directories
Constructors create format-specific directories, e.g.:
- `vtk/`, `pvtk/`, `bin/`, `rst/`, `trk/`, `tab/`
- `cbin_<file_id>_<coarsen_factor>/`
- `pdf_<file_id>` (with a second variable suffix for 2D PDFs); writes
  `<basename>.bins.pdf` once and `<basename>.#####.pdf` per output time
Rank-specific subdirectories are created when `single_file_per_rank` is enabled.

---

## Particle VTK Content Notes
- `pvtk` outputs include:
  - `POINTS` coordinates for all particles,
  - big-endian 32-bit integer scalar attributes with explicit names (`gid`,
    `ptag`, plus `species` for CR particles or `sn_id` for star particles;
    fallback names use `idata_<n>`),
  - velocity vectors (`VECTORS vel`).
- This is the preferred path for publication-style particle phase-space
  diagnostics without parsing restart files.

---

## Cautions
- `FormattedTableOutput` enforces 1D slices (fails on 2D/3D without slice planes).
- Derived variables in `derived_variables.cpp` operate on active zones only.
- Keep Kokkos data layout conventions intact when adding outputs or variables.
