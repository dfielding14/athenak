# Configuring AthenaK

AthenaK reads runtime parameters from `.athinput` files. Each file consists of labelled blocks (`<blockname>`) followed by `name = value` entries and optional trailing comments (`# comment`). Blocks can appear more than once; later assignments override earlier ones.

## Example Layout

The following is an abbreviated orientation extract. Use the complete shipped
`inputs/hydro/sod.athinput` deck for a runnable starting point.

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
nx2 = 1
x2min = 0.0
x2max = 1.0
nx3 = 1
x3min = 0.0
x3max = 1.0

<meshblock>
nx1 = 256
nx2 = 1
nx3 = 1

<time>
evolution = dynamic
cfl_number = 0.4             # Courant number
tlim = 10.0                  # Stop time
integrator = rk2             # Time integrator; see supported list below

<hydro>
eos = ideal                  # Equation of state
gamma = 1.4                  # Adiabatic index
reconstruct = plm            # Reconstruction scheme
rsolver = hllc               # Riemann solver
```

Existing parameters can be overridden from the command line using the
`block/name=value` syntax, e.g. `./athena -i inputs/hydro/sod.athinput
mesh/nx1=512`. An override for a missing block or parameter is rejected.

## Essential Blocks

| Block | Purpose | Notes |
|-------|---------|-------|
| `<job>`  | Naming, metadata | `basename` prefixes every output file |
| `<mesh>` | Grid setup | `nx[1-3]`, physical extents, `nghost`, boundary flags; combine with `<meshblock>` to control per-block resolution |
| `<time>` | Time integration | `evolution`, `cfl_number`, `tlim`, `integrator`, `nlim`; CLI flag `-t hh:mm:ss` sets a wall-clock limit |
| `<output#>` | Diagnostics/output | `file_type`, cadence (`dt` or `dcycle`), `variable` where required, and format-specific options |
| `<hydro>` / `<mhd>` | Enable fluid/MHD modules | Presence of the block activates the module |
| `<radiation>`, `<z4c>`, ... | Additional physics | Include only when the corresponding physics is needed |

### Mesh Parameters

- `nx1`, `nx2`, `nx3`: Active cells along each direction. Set unused directions to 1.
- `nghost`: Number of ghost cells per side (defaults to 2).
- `x1min`, `x1max`, etc.: Physical domain limits.
- Boundary flags (`ix1_bc`, `ox1_bc`, ...): `periodic`, `outflow`, `reflect`,
  `inflow`, `diode`, `vacuum`, or `user`; `shear_periodic` is accepted only
  for the x1 faces of an appropriately configured shearing-box case.
- `inflow` requires the selected problem generator to initialize the constant
  inflow state for that face; setting the boundary flag alone is insufficient.

- `<meshblock>/nx{1,2,3}` controls each block's cell dimensions and decomposition. When omitted, each mesh block defaults to the global `nx` values (`docs/source/reference/input_parameters.md`).

### Time Integration

- `cfl_number`: CFL multiplier used when the driver reduces module timestep estimates.
- `tlim` / `nlim`: Stop conditions by time or cycle count.
- `integrator`: `rk1`, `rk2`, `rk3`, `rk4`, `imex2`, `imex2+`, or `imex3`.
- Wall-clock limits are supplied via `./athena -t hh:mm:ss`.

### Hydrodynamics / MHD

Presence of `<hydro>` enables fluid state and `<mhd>` enables magnetic fluid
state. The coordinate configuration selects non-relativistic, special
relativistic, or general relativistic behavior where implemented. The solver
rows below describe the non-relativistic/kinematic configurations:

| Parameter | Values | Description |
|-----------|--------|-------------|
| `eos` | `ideal`, `isothermal` | Equation of state for non-relativistic cases; SR and GR hydro/MHD reject `isothermal` |
| `gamma` | float | Adiabatic index (ideal EOS) |
| `reconstruct` | `dc`, `plm`, `ppm4`, `ppmx`, `wenoz` | Spatial reconstruction |
| `rsolver` (dynamic hydro) | `llf`, `hlle`, `hllc`, `roe` | Non-relativistic dynamic solver; `hllc` requires ideal EOS |
| `rsolver` (dynamic MHD) | `llf`, `hlle`, `hlld` | Non-relativistic dynamic MHD solver |
| `rsolver` (kinematic hydro/MHD) | `advect` | Advection-only kinematic solver |
| `viscosity`, `conductivity` | float | Optional hydro diffusion coefficients |
| `ohmic_resistivity`, `viscosity`, `conductivity` | float | Optional MHD diffusion coefficients |
| `fofc` | boolean | First-order flux correction; requires additional ghost zones as described below |

The default `nghost = 2` supports `dc` and `plm` without FOFC. The
higher-order choices `ppm4`, `ppmx`, and `wenoz` require `nghost >= 3`.
Enabling `fofc = true` requires `nghost >= 3` with `plm`, or `nghost >= 4`
with `ppm4`, `ppmx`, or `wenoz`, for both `<hydro>` and `<mhd>`.
When SMR or AMR is enabled, the mesh constructor also requires an even
`nghost`; consequently refined cases needing at least three ghost zones must
use `nghost >= 4`.

Additional physics blocks have compatibility constraints; check the relevant
module page and a shipped input deck before combining them.

### Output Blocks

Each `<output#>` block defines one stream:

```ini
<output1>
file_type = vtk           # vtk, tab, hst, log, pvtk, trk, cbin, pdf, bin, cart, sph, rst
dt = 0.1                  # Output cadence
variable = hydro_w        # Required for field-bearing formats
```

Common optional fields include `id`, `ghost_zones`, `gid`, and slice positions.
`data_format` controls formatted table/history values.
`single_file_per_rank` is read for `bin`, `cbin`, and `rst` output streams.

Multiple outputs can be defined: `<output2>`, `<output3>`, and so on.
Restart dumps use `file_type = rst` and are controlled by `dt` or `dcycle`.

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
- Use `./athena -i input.athinput -n` to inspect parsed parameters and
  command-line overrides without running. It exits before mesh and physics
  construction, so it does not validate module compatibility or mesh rules.
- Runtime overrides (command-line assignments) apply after the file is parsed.
