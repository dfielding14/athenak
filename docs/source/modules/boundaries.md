# Module: Boundary Values

## Overview
The Boundary Values (bvals) module handles boundary conditions, ghost zone exchanges, and MPI communication for parallel simulations, supporting various physics-specific boundary types.

## Source Location
`src/bvals/`

## Key Components

| File | Purpose | Key Functions |
|------|---------|---------------|
| `bvals.hpp/cpp` | Core boundary class | BC application, MPI setup |
| `bvals_cc.cpp` | Cell-centered boundaries | Scalar/vector fields |
| `bvals_fc.cpp` | Face-centered boundaries | Magnetic fields |
| `buffs_cc.cpp` | Cell-centered buffers | MPI communication |
| `buffs_fc.cpp` | Face-centered buffers | B-field communication |
| `flux_correct_cc.cpp` | Flux correction | AMR conservation |
| `prolongation.cpp` | AMR prolongation | Fine→coarse |
| `prolong_prims.cpp` | Primitive prolongation | Better accuracy |

### Physics-Specific BCs
| File | Purpose |
|------|---------|
| `physics/hydro_bcs.cpp` | Hydro boundary conditions |
| `physics/bfield_bcs.cpp` | Magnetic field BCs |
| `physics/radiation_bcs.cpp` | Radiation BCs |
| `physics/z4c_bcs.cpp` | GR metric BCs |

## Boundary Types

### Standard Boundaries

| Type | Description | Implementation |
|------|-------------|----------------|
| `periodic` | Periodic wrap | Copy from opposite boundary |
| `outflow` | Zero gradient | Copy from last active zone |
| `reflect` | Mirror symmetry | Reflect normal velocity |
| `inflow` | Fixed inflow | Set to specified values |
| `diode` | One-way flow | Allows outflow, prevents inflow |
| `vacuum` | Vacuum BC | Zero density/pressure |
| `shear_periodic` | Shearing box | Periodic with orbital shear |
| `user` | User-defined | Call problem generator |

### Configuration
```ini
<mesh>
ix1_bc = outflow    # inner-x1 boundary
ox1_bc = outflow    # outer-x1 boundary
ix2_bc = periodic   # inner-x2 boundary
ox2_bc = periodic   # outer-x2 boundary
ix3_bc = reflect    # inner-x3 boundary
ox3_bc = reflect    # outer-x3 boundary
```

## Implementation Details

### Ghost Zones
Default: 2 ghost zones
```ini
<mesh>
nghost = 2  # Can increase for high-order methods
```

### Reflect BC Example
```cpp
// Reflect normal velocity component
if (reflect_x1) {
  prim(IVX, k, j, -i) = -prim(IVX, k, j, i);  // Flip vx
  prim(IVY, k, j, -i) = prim(IVY, k, j, i);   // Keep vy
  prim(IVZ, k, j, -i) = prim(IVZ, k, j, i);   // Keep vz
}
```

## MPI Communication

### Exchange Pattern
```{mermaid}
flowchart LR
    MB1[MeshBlock 1] -->|Send Buffer| MB2[MeshBlock 2]
    MB2 -->|Send Buffer| MB1
    MB3[MeshBlock 3] -->|Send Buffer| MB4[MeshBlock 4]
    MB4 -->|Send Buffer| MB3
```

### Communication Steps
1. **Pack**: Copy boundary data to send buffers
2. **Send**: MPI_Isend to neighbors
3. **Receive**: MPI_Irecv from neighbors
4. **Unpack**: Copy from receive buffers to ghost zones

### Buffer Management
```cpp
// Non-blocking communication
MPI_Isend(send_buffer, size, MPI_DOUBLE, 
          neighbor_rank, tag, comm, &request);
          
MPI_Irecv(recv_buffer, size, MPI_DOUBLE,
          neighbor_rank, tag, comm, &request);
          
// Wait for completion
MPI_Wait(&request, &status);
```

## AMR Boundaries

### Prolongation (Coarse→Fine)
```cpp
// Linear interpolation for refinement
fine[2i]   = coarse[i]
fine[2i+1] = 0.5*(coarse[i] + coarse[i+1])
```

### Restriction (Fine→Coarse)
```cpp
// Volume-weighted average
coarse[i] = 0.125*sum(fine[2i:2i+1, 2j:2j+1, 2k:2k+1])
```

### Flux Correction
```cpp
// Ensure conservation at coarse-fine boundaries
F_coarse = sum(F_fine * area_fine) / area_coarse
```

## Magnetic Field BCs

### Divergence Preservation
Face-centered fields require special treatment:
```cpp
// Ensure div(B) = 0 at boundaries
Bx_ghost = Bx_active  // Continuous
By_ghost = By_active + correction
Bz_ghost = Bz_active + correction
```

## User-Defined Boundaries

### Implementation
```cpp
// In problem generator
void UserBoundary(MeshBlockPack *pmbp) {
  // Custom boundary conditions
  if (time > t_inject) {
    // Inject material
    prim(IDN, k, j, 0) = rho_inject;
    prim(IVX, k, j, 0) = v_inject;
  }
}
```

### Usage
User boundaries are called automatically when boundary type is set to "user" in input file:
```ini
<mesh>
ix1_bc = user   # Will call UserBoundary function
```

## Shearing Box Boundaries

### Shearing Periodic
```cpp
// Special boundaries for orbital shear
y_shift = q * Omega * x * t
Copy with shifted y-coordinate
```

## Performance Optimization

### Communication Overlap
```cpp
// Overlap computation with communication
StartBoundaryExchange();
ComputeInterior();
WaitForBoundaries();
ComputeBoundaries();
```

### Buffer Packing
- Contiguous memory for MPI efficiency
- Vectorized packing/unpacking
- Minimize message count

## Common Issues

### Boundary Artifacts
- Check BC implementation
- Verify ghost zone filling
- Test with uniform flow

### MPI Deadlock
- Use non-blocking communication
- Check send/receive pairing
- Verify tag uniqueness

### AMR Conservation
- Enable flux correction
- Check prolongation order
- Monitor conservation errors

## Testing

Boundary condition tests:
- Uniform advection (periodic)
- Shock reflection (reflect)
- Wave propagation (outflow)
- Inflow injection (inflow)
- One-way valve (diode)

## See Also
- [Mesh Module](mesh.md)
- Source: `src/bvals/bvals.cpp`