# Module: Outputs

## Overview
The Outputs module manages all data output from AthenaK simulations, supporting 14 different formats for mesh data, particle data, and diagnostics.

## Source Location
`src/outputs/`

## Output Formats (VERIFIED from outputs.cpp)

### Mesh Data Outputs

| Format | File Type | Extension | Description | Source File |
|--------|-----------|-----------|-------------|-------------|
| `vtk` | MeshVTKOutput | `.vtk` | VTK legacy format for ParaView/VisIt | `vtk_mesh.cpp` |
| `bin` | MeshBinaryOutput | `.bin` | Native binary format (fastest I/O) | `binary.cpp` |
| `rst` | RestartOutput | `.rst` | Restart files for exact state preservation | `restart.cpp` |
| `tab` | FormattedTableOutput | `.tab` | ASCII table format | `formatted_table.cpp` |
| `cbin` | CoarsenedBinaryOutput | `.cbin` | Coarsened binary (reduced resolution) | `coarsened_binary.cpp` |

### Particle Outputs (if particles enabled)

| Format | File Type | Extension | Description | Source File |
|--------|-----------|-----------|-------------|-------------|
| `pvtk` | ParticleVTKOutput | `.pvtk` | Particle VTK format | `vtk_prtcl.cpp` |
| `trk` | TrackedParticleOutput | `.trk` | Tracked particle trajectories | `track_prtcl.cpp` |
| `df` | ParticleDFOutput | `.df` | Particle distribution functions | `df_prtcl.cpp` |
| `dxh` | ParticleDxHistOutput | `.dxh` | Particle displacement histograms | `dxhist_prtcl.cpp` |
| `ppd` | ParticlePositionsOutput | `.ppd` | Particle position dumps | `posvel_prtcl.cpp` |
| `prst` | ParticleRestartOutput | `.prst` | Particle restart data | `restart_prtcl.cpp` |

### Diagnostic Outputs

| Format | File Type | Extension | Description | Source File |
|--------|-----------|-----------|-------------|-------------|
| `hst` | HistoryOutput | `.hst` | Time series data | `history.cpp` |
| `log` | EventLogOutput | `.log` | Event counter logs | `event_log.cpp` |
| `pdf` | PDFOutput | `.pdf` | Probability density functions | `pdf.cpp` |

## Configuration Parameters

### Basic Output Block
```ini
<output1>
file_type = vtk              # Required: output format
dt = 0.1                     # Time between outputs
variable = cons              # cons (conserved) or prim (primitive)
dirname = data              # Output directory (optional)
```

### Advanced Parameters

| Parameter | Type | Default | Description | Formats |
|-----------|------|---------|-------------|---------|
| `single_file_per_rank` | bool | false | Each MPI rank writes separately | bin, rst, cbin |
| `ghost_zones` | bool | false | Include ghost cells in output | vtk, bin |
| `coarsen_factor` | int | 2 | Coarsening factor (power of 2) | cbin only |
| `tlim` | Real | - | Stop output after this time | all |
| `nlim` | int | - | Maximum number of outputs | all |
| `slice_x1` | Real | - | Output x1 slice at this position | vtk, bin |
| `slice_x2` | Real | - | Output x2 slice at this position | vtk, bin |
| `slice_x3` | Real | - | Output x3 slice at this position | vtk, bin |

## File Naming Convention

```
basename.block[BLOCKID].out[OUTPUTID].[CYCLE].[FORMAT]
```

Examples:
- `shock.block0.out1.00000.vtk` - VTK output, block 0, cycle 0
- `turb.00100.rst` - Restart file at cycle 100
- `run.00050.hst` - History file at cycle 50

## Implementation Details

### Output Manager (`outputs.cpp`)
```cpp
// Line 95-200: Output type registration
if (file_type == "vtk") {
  pnew_out = new MeshVTKOutput(pout_params, pin, pmesh);
} else if (file_type == "bin") {
  pnew_out = new MeshBinaryOutput(pout_params, pin, pmesh);
} else if (file_type == "rst") {
  pnew_out = new RestartOutput(pout_params, pin, pmesh);
}
// ... etc for all 14 formats
```

### Key Classes
- `BaseOutput`: Abstract base class for all outputs
- `OutputParameters`: Configuration for each output
- `OutputManager`: Coordinates all outputs

## Usage Examples

### Basic VTK Output
```ini
<output1>
file_type = vtk
dt = 0.1
variable = prim
```

## PDF Output (N-D)

PDF outputs are N-dimensional histograms (1D–4D) written as binary arrays plus a
one-time ASCII header. Each PDF is stored in a directory named `pdf_<id>` (and
includes all PDF variables in the directory name for uniqueness).

### PDF Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `variable_1..variable_4` | string | - | PDF variables (1–4 dimensions). Each must be a **scalar** output variable. |
| `binN_min` / `binN_max` | Real | - | Bin range for dimension N. |
| `nbinN` | int | - | Number of bins for dimension N. |
| `logscaleN` | bool | false | Logarithmic binning for dimension N. Requires `binN_min > 0`. |
| `weight` | string | `volume` | `volume`, `mass`, or `variable`. |
| `weight_variable` | string | - | Required when `weight=variable`. Must be a **scalar** output variable. |

Notes:
- `mass_weighted` is **not supported**. Use `weight=mass`.
- The weight is **always multiplied by the cell volume** (`dx1*dx2*dx3`), so AMR
  cells contribute correctly by default.
- PDFs are computed over **active zones only** (ghost zones excluded).
- If a variable expands to multiple components (e.g., a group like `mhd_u` or a
  multi-moment derived variable), the run will error; PDFs require scalar variables.

### Output Layout

Each output time writes:
- `pdf_<id>/<basename>.header.pdf` (ASCII, written once)
- `pdf_<id>/<basename>.NNNNN.pdf` (binary per output)

Binary file format:
1. `double` time
2. `double` histogram values, flattened row-major across dimensions with
   **underflow/overflow bins** included in each dimension.

### Coordinate Variables for PDF Binning

These built-in derived variables can be used as PDF dimensions:
- Cartesian: `coord_x`, `coord_y`, `coord_z`
- Spherical: `coord_r`, `coord_theta`, `coord_phi`, `coord_costheta`
- Cylindrical: `coord_cyl_R`, `coord_cyl_phi`, `coord_cyl_z`

Important caveats:
- These coordinates are computed from `x1/x2/x3` using `CellCenterX`, which assumes
  **uniform spacing** in each coordinate direction.
- The spherical/cylindrical variants are derived from Cartesian `x/y/z`. They are
  correct for **Cartesian meshes**, but not for runs whose native mesh coordinates
  are spherical or cylindrical.
- Use `logscaleN` only for positive coordinates (e.g., `coord_r`, `coord_cyl_R`).

### Velocity and Flux Derived Variables (for PDFs)

These derived variables are useful for flux PDFs and multi-dimensional binning:
- Spherical velocity components: `vel_sph_r`, `vel_sph_theta`, `vel_sph_phi`
- Cylindrical velocity components: `vel_cyl_R`, `vel_cyl_phi`
- Radial mass flux (spherical): `mdot_sph` (and `mdot_sph_out`, `mdot_sph_in`)
- Radial energy flux (spherical): `edot_sph` (and `edot_sph_out`, `edot_sph_in`)
- Radial energy flux components (spherical):
  - `edot_sph_kin` (kinetic)
  - `edot_sph_th` (thermal/enthalpy)
  - `edot_sph_mag` (magnetic, MHD only)

#### Coordinate-based velocity definitions

Let `(x, y, z)` be the Cartesian cell-center position and `R = sqrt(x^2 + y^2)`,
`r = sqrt(x^2 + y^2 + z^2)`. For velocity `(vx, vy, vz)`:

- Spherical radial:
  - `v_r = (vx*x + vy*y + vz*z) / r`
- Spherical polar (theta, measured from +z):
  - `v_theta = (z * (vx*x + vy*y) / (r*R)) - vz*(R/r)` (set to 0 if `r=0` or `R=0`)
- Spherical azimuthal:
  - `v_phi = (-vx*y + vy*x) / R` (set to 0 if `R=0`)
- Cylindrical radial:
  - `v_R = (vx*x + vy*y) / R` (set to 0 if `R=0`)
- Cylindrical azimuthal:
  - `v_phi = (-vx*y + vy*x) / R` (set to 0 if `R=0`)

These formulas assume the mesh coordinates are Cartesian (`x1/x2/x3 = x/y/z`).

#### Flux derivations (ideal gas)

For ideal gas with internal energy density `e_int` and pressure
`P = (gamma - 1) * e_int`, the conservative energy fluxes are:

- **Hydro** energy flux:
  - `F_E = (E + P) * v`, where `E = e_int + 0.5*rho*v^2`
- **MHD** energy flux:
  - `F_E = (E + P + 0.5*B^2) * v - (v · B) * B`,
    where `E = e_int + 0.5*rho*v^2 + 0.5*B^2`

Taking the radial component (`v_r`, `B_r`) gives:

- Total radial energy flux:
  - `F_E,r = (0.5*rho*v^2 + gamma*e_int + B^2) * v_r - (v · B) * B_r`
  - Hydro case: set `B = 0`

The component fluxes used by the new derived variables are:

- Kinetic energy flux:
  - `F_kin,r = 0.5 * rho * v^2 * v_r`  (`edot_sph_kin`)
- Thermal (enthalpy) energy flux:
  - `F_th,r = gamma * e_int * v_r`      (`edot_sph_th`)
- Magnetic energy flux (MHD only):
  - `F_mag,r = B^2 * v_r - (v · B) * B_r` (`edot_sph_mag`)

By construction:
`F_E,r = F_kin,r + F_th,r + F_mag,r`.

### Selected Derived Variables (Mesh Outputs)

These derived variables are available to mesh outputs (vtk/bin/rst/tab/cbin) and
can also be used as PDF dimensions or weights when they are scalar.

- `mhd_divb`: Divergence of the face-centered magnetic field. Computed over
  active + ghost zones in all available dimensions (x1/x2/x3); ghost-zone values
  are only written if `ghost_zones=true` for that output.
- `mhd_curv`: Magnetic curvature magnitude from cell-centered B and its
  gradients. Cells with `|B| = 0` return `0` (avoids divide-by-zero).
- `mhd_dynamo_ks`: Eight diagnostic scalars written as
  `mhd_dynamo_B^2`, `mhd_dynamo_B^4`, `mhd_dynamo_dB^2`,
  `mhd_dynamo_BdB^2`, `mhd_dynamo_|BxJ|^2`, `mhd_dynamo_|B.J|^2`,
  `mhd_dynamo_U^2`, `mhd_dynamo_dU`. Uses finite differences in all available
  dimensions (includes x3 when 3D).
- `prtcl_d`: Particle number per cell (nearest-cell binning of particle
  positions). Convert to number density by dividing by cell volume; multiply by
  particle mass for a mass density (if masses are uniform).

Notes:
- `mhd_*` variables require MHD; `prtcl_d` requires particles enabled.

#### Using PDFs for Mdot in bins of (r, theta, T, v_r)

For mass flux PDFs, set the PDF weight to the radial mass flux and bin by
`coord_r`, `coord_theta`, `temperature`, and `vel_sph_r`. Because PDF weights
always include cell volume, the binned sum is effectively
`sum(volume * rho * v_r)`. You can divide by the radial bin width in
post-processing to form a flux density in `r`.

Example:
```ini
<output_mdot_pdf>
file_type = pdf
dt = 0.1
variable_1 = coord_r
bin1_min = 1e-2
bin1_max = 10
nbin1 = 128
logscale1 = true
variable_2 = coord_theta
bin2_min = 0
bin2_max = 3.14159
nbin2 = 64
variable_3 = temperature
bin3_min = 1e-3
bin3_max = 1e2
nbin3 = 64
logscale3 = true
variable_4 = vel_sph_r
bin4_min = -1
bin4_max = 1
nbin4 = 128
weight = variable
weight_variable = mdot_sph
```

### PDF Examples

**1D mass-weighted PDF**
```ini
<output1>
file_type = pdf
dt = 0.1
variable_1 = mhd_u_d
bin1_min = 1e-6
bin1_max = 1e-1
nbin1 = 128
logscale1 = true
weight = mass
```

**2D PDF weighted by a derived scalar (current density)**
```ini
<output2>
file_type = pdf
dt = 0.2
variable_1 = mhd_bcc1
bin1_min = -1
bin1_max = 1
nbin1 = 200
variable_2 = mhd_j2
bin2_min = 0
bin2_max = 10
nbin2 = 200
weight = variable
weight_variable = mhd_j2
```

**3D PDF including coordinates**
```ini
<output3>
file_type = pdf
dt = 0.5
variable_1 = coord_x
bin1_min = -1
bin1_max = 1
nbin1 = 128
variable_2 = coord_y
bin2_min = -1
bin2_max = 1
nbin2 = 128
variable_3 = coord_z
bin3_min = -1
bin3_max = 1
nbin3 = 128
weight = volume
```

### Parallel Binary Output
```ini
<output2>
file_type = bin
dt = 1.0
single_file_per_rank = true  # Important for large parallel runs
ghost_zones = false
```

### Restart Configuration
```ini
<output3>
file_type = rst
dt = 10.0  # Checkpoint every 10 time units
```

### Coarsened Output for Large Runs
```ini
<output4>
file_type = cbin
dt = 0.5
coarsen_factor = 4  # Reduce resolution by 4x
```

### History Output
```ini
<output5>
file_type = hst
dt = 0.01  # High frequency for time series
```

## Performance Considerations

### I/O Performance
| Format | Speed | Size | Use Case |
|--------|-------|------|----------|
| bin | Fastest | Medium | Production runs |
| vtk | Slow | Large | Visualization |
| rst | Fast | Large | Checkpointing |
| cbin | Fast | Small | Large datasets |
| tab | Slowest | Large | Small test problems |

### Parallel I/O
- Use `single_file_per_rank=true` for bin/rst on parallel filesystems
- VTK always writes one file per MeshBlock
- History files are written by rank 0 only

## Common Issues and Solutions

### Issue: Output files too large
**Solution**: Use `cbin` with `coarsen_factor` or output less frequently

### Issue: Too many files in parallel runs
**Solution**: Use `single_file_per_rank=false` (default) for bin format

### Issue: Can't visualize binary files
**Solution**: Use VTK format for visualization, binary for analysis

## Binary Format Recommended
**Note**: AthenaK uses native binary format for efficient parallel I/O. No external libraries required.

## Testing
```bash
# Test output format
./athena -i test.athinput
ls *.vtk *.bin *.hst  # Check generated files
```

## See Also
- [Configuration Guide](../configuration.md)
- [Running Simulations](../running.md)
- Source: `src/outputs/outputs.cpp` L95-200 for format registration
