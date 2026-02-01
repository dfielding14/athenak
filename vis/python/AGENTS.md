# AGENTS.md

## Purpose
This directory contains Python utilities for post-processing AthenaK outputs:
readers for AthenaK file formats, binary-to-athdf/XDMF conversion helpers, and
standalone plotting/analysis scripts. These tools are not used by the runtime
solver. See `../../AGENTS.md` for repository-wide conventions and workflow.

---

## Key Files and Responsibilities

### `__init__.py`
- Empty; marks `vis/python` as a package.

### `athena_read.py`
- `check_nan_flag` and `check_nan(data)`: optional NaN validation (off by
  default).
- `error_dat(filename, **kwargs)`: `numpy.loadtxt` wrapper used in regression
  tests; keeps 2D shape and can NaN-check.
- `tab(filename)`: reads `.tab` outputs, parses the `time=` / `cycle=` header,
  drops the first data column (cell index), and returns a dict of arrays keyed
  by column headings.
- `hst(filename, raw=False)`: reads `.hst` history files, uses the most recent
  header if multiple exist, and (when `raw=False`) prunes restart branches by
  enforcing monotonic time.
- `athdf(...)`: reads `.athdf` output using `h5py`.
  - `raw=True` returns file attributes plus `Levels`, `LogicalLocations`,
    `x1f/x2f/x3f`, `x1v/x2v/x3v`, and each variable dataset.
  - `raw=False` merges MeshBlocks into a single grid, supports level selection,
    optional ghost zones, and spatial subsetting via `x1/2/3_min` and
    `x1/2/3_max`.
  - Uses `face_func_1/2/3` for user-defined face positions (when `xrat_root` is
    negative) and `vol_params[0]` for Kerr-Schild volume integrals when
    restriction is required.
- `AthenaError`: custom `RuntimeError` subclass raised by readers.

### `bin_convert.py` (legacy binary reader)
- `read_binary(filename)`: parses Athena binary output v1.1 and returns a
  `filedata` dict with:
  `header`, `time`, `cycle`, `var_names`, `nvars`, `Nx1/2/3`, `x1/2/3min/max`,
  `n_mbs`, `nx1/2/3_mb`, `nx1/2/3_out_mb`, `mb_index`, `mb_logical`,
  `mb_geometry`, `mb_data`.
- `read_coarsened_binary(filename)`: same format plus `coarsening factor` and
  `number of moments`, and scaled `Nx*`/`nx*_mb`.
- `read_all_ranks_binary(rank0_filename)`: globs `rank_*` directories, filters
  mismatched file sizes, and concatenates MeshBlock arrays.
- `read_all_ranks_coarsened_binary(rank0_filename)`: same for coarsened dumps
  (no size filter).
- `read_*_as_athdf(...)`: builds athdf-like arrays in memory with cartesian
  coordinates and simple volume/center functions. Prolongation uses
  `np.repeat`; restriction is a stub (`pass`).
- `read_single_rank_binary_as_athdf(...)`: returns one meshblock (index 0) as an
  athdf-like dict; uses `mb_geometry` columns `[0/1]`, `[2/3]`, `[4/5]` as
  coordinate bounds when building `x1f/x2f/x3f`.
- `read_coarsened_binary_as_athdf(...)`: coarsened version of the above.
- `write_xdmf_for(xdmfname, dumpname, fdata)`: writes XDMF referencing HDF5
  datasets `/uov` (hydro) and `/B` (magnetic fields).
- `convert_file(binary_fname)`: calls `read_binary`, then `write_athdf`, then
  `write_xdmf_for`; CLI converts one or more `.bin` files.

### `bin_convert_new.py` (updated binary reader)
- Module constants: `CODE_HEADER_EXPECTED`, `SUPPORTED_VERSION`,
  `COORD_CARTESIAN`.
- `_get_from_header(...)`: internal helper for header parsing.
- `read_binary` / `read_coarsened_binary`: wrappers around `_read_binary_impl`
  (BytesIO reader).
- `read_all_ranks_binary` / `read_all_ranks_coarsened_binary`: combine `rank_*`
  files; binary variant filters mismatched file sizes.
- `read_binary_as_athdf`, `read_all_ranks_*_as_athdf`,
  `read_coarsened_binary_as_athdf`: mirror `bin_convert.py` behavior (cartesian
  coordinates, prolongation via `np.repeat`, restriction stub).
- `read_single_rank_binary_as_athdf(filename, meshblock_index_in_file, ...)`:
  extracts a specific meshblock index.
- `read_rank_binary_as_athdf(filename, ...)`: merges all meshblocks in a rank
  into a full grid using `mb_index` placement.
- `write_xdmf_for`: XDMF writer with the same dataset conventions as
  `bin_convert.py`.
- `convert_file`: calls `write_athdf` then `write_xdmf_for`.
- `athinput(filename)`: parses athinput files into nested dicts with typecasts
  for ints/floats/complex (strings containing `_` remain strings).
- `__all__` exports the public API list at the bottom of the file.

### `make_athdf.py`
- CLI for converting a stem of `.bin` files to `.athdf`/`.xdmf` using
  `bin_convert.read_binary` and `bin_convert.write_athdf`.

---

## Plotting and Analysis Scripts

### `calculate_tori_equil.py`
- Generates a 2x2 figure of Chakrabarti and Fishbone-Moncrief torus density and
  pressure using default parameters at the top of the file.
- Writes `gr_equilibria.png` in the current working directory.
- Defines helper functions: `geometry`, `metric`, `l_k`, `vz`, `c_cn`, `c_l`,
  `c_h`, `c_u_t`, `fm_ls`, `fm_f`, `fm_h`, `fm_h_vec`, `calculate_rho`,
  `calculate_tt`.

### `calculate_tori_magnetization.py`
- CLI: computes volume- and mass-weighted plasma `sigma` and `beta_inv` from a
  GRMHD `.bin` dump.
- Inputs: file name, optional `--r_max` and `--rho_min`.
- Reads binary dumps directly; expects `dens`, `eint`, `velx/vely/velz`,
  `bcc1/2/3`.

### `calculate_tori_rpeak.py`
- CLI: computes `r_peak` for Fishbone-Moncrief (`fm`) or Chakrabarti (`c`) tori
  using `brentq`.
- Inputs: `torus_type`, `spin`, `r_in`, `r_out`, optional `--n`.

### `plot_hst.py`
- CLI: plots a single variable from `.hst` using `athena_read.hst`.
- Prints the data dict and calls `plt.show()`; the output argument only sets a
  non-interactive backend.

### `plot_mesh.py`
- CLI: plots `mesh_structure.dat` (generated by `athena -m`) as 3D line
  segments.
- Blank lines delimit meshblock edge segments.

### `plot_slice.py`
- CLI: plots 2D slices from 2D/3D `.bin` outputs; supports grid, horizon,
  horizon mask, and ergosphere overlays.
- Reads binary dumps directly, computes per-block extents, and loads requested
  variables.
- Derived quantities are requested via the `derived:` prefix and computed
  in-code using:
  - `set_derived_dependencies()` (required base variables)
  - `set_labels()` (colorbar labels)
  - coordinate helpers for CKS/SKS metrics and transforms

#### Derived quantity names (from `set_derived_dependencies`)
Gas and thermodynamics:
- `pgas`, `pgas_rho`, `T`, `prad_pgas`

Velocities and 4-velocities:
- `vr_nr`, `vth_nr`, `vph_nr`
- `uut`, `ut`, `ux`, `uy`, `uz`, `ur`, `uth`, `uph`
- `u_t`, `u_x`, `u_y`, `u_z`, `u_r`, `u_th`, `u_ph`
- `vx`, `vy`, `vz`, `vr_rel`, `vth_rel`, `vph_rel`

Magnetic fields and pressures:
- `Br_nr`, `Bth_nr`, `Bph_nr`, `pmag_nr`, `beta_inv_nr`, `sigma_nr`
- `bt`, `bx`, `by`, `bz`, `br`, `bth`, `bph`
- `b_t`, `b_x`, `b_y`, `b_z`, `b_r`, `b_th`, `b_ph`
- `Br_rel`, `Bth_rel`, `Bph_rel`, `pmag_rel`, `beta_inv_rel`, `sigma_rel`,
  `sigmah_rel`, `va_rel`
- `pmag_prad`

Radiation tensor and Eddington quantities:
- `prad`
- `Rtr`, `Rtth`, `Rtph`, `Rrr`, `Rthth`, `Rphph`, `Rrth`, `Rrph`, `Rthph`
- `Rtx_Rtt`, `Rty_Rtt`, `Rtz_Rtt`, `Rxx_Rtt`, `Ryy_Rtt`, `Rzz_Rtt`,
  `Rxy_Rtt`, `Rxz_Rtt`, `Ryz_Rtt`
- `Rtr_Rtt`, `Rtth_Rtt`, `Rtph_Rtt`, `Rrr_Rtt`, `Rthth_Rtt`, `Rphph_Rtt`,
  `Rrth_Rtt`, `Rrph_Rtt`, `Rthph_Rtt`
- `R01_R00_ff`, `R02_R00_ff`, `R03_R00_ff`, `R11_R00_ff`, `R22_R00_ff`,
  `R33_R00_ff`, `R12_R00_ff`, `R13_R00_ff`, `R23_R00_ff`

Opacities and optical depths:
- `kappa_a`, `kappa_s`, `kappa_t`
- `alpha_a`, `alpha_s`, `alpha_t`
- `tau_a`, `tau_s`, `tau_t`

Enthalpy and Bernoulli:
- `wgas`, `wmhd`, `wgasrad`, `wmhdrad`
- `Begas`, `Bemhd`, `Begasrad`, `Bemhdrad`

Conserved quantities:
- `cons_hydro_nr_t`, `cons_hydro_nr_x`, `cons_hydro_nr_y`, `cons_hydro_nr_z`
- `cons_em_nr_t`
- `cons_mhd_nr_t`, `cons_mhd_nr_x`, `cons_mhd_nr_y`, `cons_mhd_nr_z`
- `cons_hydro_rel_t`, `cons_hydro_rel_x`, `cons_hydro_rel_y`, `cons_hydro_rel_z`
- `cons_em_rel_t`, `cons_em_rel_x`, `cons_em_rel_y`, `cons_em_rel_z`
- `cons_mhd_rel_t`, `cons_mhd_rel_x`, `cons_mhd_rel_y`, `cons_mhd_rel_z`

### `plot_tab.py`
- CLI: plots or animates `.tab` outputs using `athena_read.tab`.
- `Player` (subclass of `FuncAnimation`) adds a slider and playback controls.
- Assumes a filename pattern with 5-digit indices (prefix + `NNNNN` + `.tab`).
- Chooses the x-axis variable by checking `x1v`, then `x2v`, then `x3v` (the
  last one found wins).

### `visualize_turb_field.py`
- CLI: visualizes turbulence driving outputs (`force1`, `force2`, `force3`) as
  a 2D slice.
- Uses `bin_convert_new.read_binary_as_athdf` and `read_binary` for MeshBlock
  metadata.
- Produces a 6-panel figure (magnitude, components, quiver, line cut, FFT
  power) and saves a PNG.
- If no filename is given, searches `bin/` for the latest `turb_force*.bin`.

---

## Cautions and Quirks

- `write_athdf` is referenced by `bin_convert.py`, `bin_convert_new.py`, and
  `make_athdf.py` but is not defined in this directory. Calls to
  `convert_file` or `make_athdf.py` will fail unless `write_athdf` is provided
  elsewhere.
- The `read_*_as_athdf` helpers in both `bin_convert` modules leave restriction
  unimplemented (`pass`); `subsample` and `fast_restrict` do not perform actual
  restriction.
- `bin_convert*_as_athdf` hardcodes `coord = 'cartesian'` and does not read
  coordinate metadata from files.
- `athena_read.athdf` requires `vol_params[0]` (spin `a`) for Kerr-Schild volume
  restriction when `vol_func` is not provided.
- `plot_hst.py` does not write output files; the output argument only selects a
  non-interactive backend.
- `plot_tab.py` assumes the `prefixNNNNN.tab` naming convention and will pick
  `x3v` over `x2v`/`x1v` if present.
- `visualize_turb_field.py` uses `mb_logical[mb, 0]` as the refinement level
  and `raw_data.get('MaxLevel', 0)` for the legend; `read_binary` does not
  populate `MaxLevel`, so the legend defaults to level 0 unless updated.
