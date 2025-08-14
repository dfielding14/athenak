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