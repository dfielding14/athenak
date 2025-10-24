# Configuring AthenaK

AthenaK reads runtime parameters from `.athinput` files. Each file consists of labelled blocks (`<blockname>`) followed by `name = value` entries and optional trailing comments (`# comment`). Blocks can appear more than once; later assignments override earlier ones.

## Example Layout

```ini
# Minimal hydrodynamic setup
<job>
basename = simulation_name   # Prefix for output files

<mesh>
nx1 = 256                    # Cells in x1
x1min = 0.0                  # Domain bounds
x1max = 1.0
ix1_bc = periodic            # Inner boundary condition
ox1_bc = periodic            # Outer boundary condition

<time>
cfl_number = 0.4             # Courant number
tlim = 10.0                  # Stop time
integrator = rk2             # Time integrator (rk1/rk2/rk3/rk4/imex*)

<hydro>
gamma = 1.4                  # Equation of state
reconstruct = plm            # Reconstruction scheme
rsolver = hllc               # Riemann solver
```

Parameters can be overridden from the command line using the `block/name=value` syntax, e.g. `./athena -i inputs/hydro/sod.athinput mesh/nx1=512`.

## Essential Blocks

| Block | Purpose | Notes |
|-------|---------|-------|
| `<job>`  | Naming, metadata | `basename` prefixes every output file |
| `<mesh>` | Grid setup | `nx[1-3]`, physical extents, `nghost`, boundary flags; combine with `<meshblock>` to control per-block resolution |
| `<time>` | Time integration | `cfl_number`, `tlim`, `dt`, `integrator`, `nlim`; CLI flag `-t hh:mm:ss` sets a wall-clock limit |
| `<output#>` | Diagnostics/output | `file_type`, cadence (`dt`/`ncycle`/`ndump`), `variable`, `single_file_per_rank`, `data_format` |
| `<hydro>` / `<mhd>` | Enable fluid/MHD modules | Presence of the block activates the module |
| `<radiation>`, `<z4c>`, ... | Additional physics | Include only when the corresponding physics is needed |

### Mesh Parameters

- `nx1`, `nx2`, `nx3`: Active cells along each direction. Set unused directions to 1.
- `nghost`: Number of ghost cells per side (defaults to 2).
- `x1min`, `x1max`, etc.: Physical domain limits.
- Boundary flags (`ix1_bc`, `ox1_bc`, ...): `periodic`, `outflow`, `reflect`, `inflow`, `diode`, `vacuum`, `shear_periodic`, or `user`.

- `<meshblock>/nx{1,2,3}` controls per-block resolution. When omitted, each mesh block defaults to the global `nx` values (`docs/source/reference/input_parameters.md`).

### Time Integration

- `cfl_number`: 0.1–0.4 for stability (module dependent).
- `tlim` / `nlim`: Stop conditions by time or cycle count.
- `dt`: Fixed timestep (omit for adaptive CFL stepping).
- `integrator`: `rk1`, `rk2`, `rk3`, `rk4`, or ImEx variants (`imex2`, `imex3`).
- Wall-clock limits are supplied via `./athena -t hh:mm:ss`.

### Hydrodynamics / MHD

Presence of `<hydro>` enables the non-relativistic solver, `<mhd>` enables MHD (adds magnetic field arrays). Key options include:

| Parameter | Values | Description |
|-----------|--------|-------------|
| `eos` | `ideal`, `isothermal` | Equation of state |
| `gamma` | float | Adiabatic index (ideal EOS) |
| `reconstruct` | `dc`, `plm`, `ppm4`, `ppmx`, `wenoz` | Spatial reconstruction |
| `rsolver` (hydro) | `llf`, `hlle`, `hllc`, `roe`, `advect` | Riemann solver |
| `rsolver` (MHD) | `llf`, `hlle`, `hlld`, `advect` | MHD flux solver |
| `viscosity`, `conductivity` | float | Optional hydro diffusion coefficients |
| `ohmic_resistivity`, `viscosity`, `conductivity` | float | Optional MHD diffusion coefficients |

Relativistic or GR modules (`<dyn_grmhd>`, `<z4c>`, `<radiation>`) follow the same pattern: adding the block activates the module, and each block documents its specific parameters.

### Output Blocks

Each `<output#>` block defines one stream:

```ini
<output1>
file_type = vtk           # vtk, tab, hst, bin, athdf, ...
dt = 0.1                  # Output cadence
variable = cons           # cons, prim, diag, user-defined
single_file_per_rank = true
```

Additional fields include `id` (label), `include_ghost_zones`, `data_format` (ascii/binary for tables), and `location` (directory override).

Multiple outputs can be defined: `<output2>`, `<output_restart>`, etc. Restart dumps typically use `file_type = rst` and are controlled by `ndump`, `ncycle`, or `dt` in the block.

## Physics-Specific Blocks

Problem generators expose their own `<problem>` block. For example:

```ini
<problem>
shock_dir = 1
xshock = 0.5
dl = 1.0
pl = 1.0
dr = 0.125
pr = 0.1
```

Other modules follow similar patterns (`<shearing_box>`, `<turb_driving>`, `<particles>`, etc.). Consult the module documentation under `docs/source/modules/` for detailed tables.

## Tips

- Anything after `#` on a line is ignored.
- Blocks may be repeated; later entries replace previous values.
- Use `./athena -i input.athinput -n` to validate and exit without running.
- Runtime overrides (command-line assignments) apply after the file is parsed.
