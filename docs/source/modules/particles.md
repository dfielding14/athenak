# Module: Particles

## Overview
The Particles module implements Lagrangian particle tracking with various pusher algorithms, interpolation schemes, and particle-mesh coupling for applications including cosmic rays, dust, and test particles.

## Source Location
`src/particles/`

## Key Components

| File | Purpose | Key Functions |
|------|---------|---------------|
| `particles.hpp/cpp` | Core particle class | Initialization, management |
| `particles_pushers.cpp` | Integration algorithms | `Boris()`, `VanLeer()` |
| `particles_tasks.cpp` | Task registration | Particle evolution tasks |

## Configuration Parameters

From `<particles>` block:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `nspecies` | int | 1 | Number of species |
| `particle_type` | string | required | Type (cosmic_ray only currently) |
| `pusher` | string | required | Integration (drift, boris_lin, boris_tsc) |
| `interpolation` | string | tsc | Interpolation (ngp, cic, tsc) |
| `ppc` | Real | 1.0 | Particles per cell |
| `assign_tag` | string | index_order | Tagging method |
| `mass_log_spacing` | Real | 1.0 | Mass spectrum spacing |
| `min_mass` | Real | 1.0 | Minimum particle mass |

## Particle Data Structure

### Per-Particle Properties
```cpp
struct Particle {
  Real x, y, z;       // Position
  Real vx, vy, vz;    // Velocity
  Real mass;          // Mass
  Real charge;        // Charge
  int tag;            // Unique ID
  int property;       // User property
};
```

## Interpolation Schemes

### Nearest Grid Point (NGP)
```cpp
// 0th order - fastest
weight = 1.0 at nearest cell
```

### Cloud-in-Cell (CIC)
1st order - linear:
$$\text{weight} = (1-dx)(1-dy)(1-dz)$$

### Triangular Shaped Cloud (TSC)
2nd order - smooth:
$$\text{weight} = W(dx) \cdot W(dy) \cdot W(dz)$$
$$W(x) = \begin{cases} \frac{3}{4} - x^2 & \text{for } |x| < 0.5 \end{cases}$$

## Pusher Algorithms

### Available Pushers

#### Drift Pusher
Simple drift in fields:
$$\mathbf{x}_{n+1} = \mathbf{x}_n + \mathbf{v}_{\text{drift}} \cdot \Delta t$$

#### Boris Pusher (boris_lin, boris_tsc)
Velocity update with E and B fields:
- boris_lin: Linear interpolation
- boris_tsc: TSC interpolation

$$\mathbf{v}^- = \mathbf{v}_n + \frac{q}{m}\mathbf{E} \cdot \frac{\Delta t}{2}$$
$$\mathbf{v}' = \mathbf{v}^- + \mathbf{v}^- \times \boldsymbol{\Omega}$$
$$\mathbf{v}^+ = \mathbf{v}^- + \mathbf{v}' \times \mathbf{S}$$
$$\mathbf{v}_{n+1} = \mathbf{v}^+ + \frac{q}{m}\mathbf{E} \cdot \frac{\Delta t}{2}$$

## Execution Flow

```{mermaid}
flowchart TD
    Start[Particle Tasks] --> Interp[Interpolate Fields]
    Interp --> Push[Push Particles]
    Push --> BC[Apply Boundaries]
    BC --> Deposit[Deposit to Mesh]
    Deposit --> Comm[MPI Communication]
    Comm --> Sort[Sort/Rebalance]
    Sort --> Done[Complete]
```

## Particle-Mesh Coupling

### Field Interpolation
Get fields at particle position:
$$\mathbf{E}_p = \text{Interpolate}(\mathbf{E}_{\text{mesh}}, x_p, y_p, z_p)$$
$$\mathbf{B}_p = \text{Interpolate}(\mathbf{B}_{\text{mesh}}, x_p, y_p, z_p)$$

### Charge/Current Deposition
Deposit particle properties to mesh:
$$\rho_{\text{mesh}} \mathrel{+}= q_p \cdot W(x_p, y_p, z_p)$$
$$\mathbf{J}_{\text{mesh}} \mathrel{+}= q_p \cdot \mathbf{v}_p \cdot W(x_p, y_p, z_p)$$

## Boundary Conditions

### Particle Boundaries
- `periodic`: Wrap position
- `outflow`: Remove particles
- Particles use same BC as mesh boundaries

## MPI Communication

### Particle Exchange
```cpp
// Particles crossing MeshBlock boundaries
if (particle.x > x_max) {
  SendToNeighbor(particle, RIGHT);
}
```

### Load Balancing
```cpp
// Redistribute particles
if (n_particles > threshold) {
  RedistributeParticles();
}
```

## Output Formats

### Particle Outputs
- `pvtk`: VTK format for visualization
- `trk`: Tracked particle trajectories
- `df`: Distribution functions
- `ppd`: Position dumps
- `prst`: Restart data

## Common Applications

### Test Particles
```ini
<particles>
particle_type = lagrangian
pusher = boris
interpolation = tsc
ppc = 10
```

### Cosmic Rays
```ini
<particles>
particle_type = cosmic_ray
pusher = boris
mass_log_spacing = 2.0
min_mass = 1e-3
```

### Dust Grains
```ini
<particles>
particle_type = dust
pusher = van_leer
interpolation = cic
```

## Performance Considerations

### GPU Optimization
```cpp
// Particle operations on GPU
Kokkos::parallel_for("push", n_particles,
KOKKOS_LAMBDA(int p) {
  // Push particle p
});
```

### Memory Layout
- Structure of Arrays (SoA) for GPU
- Minimize random access
- Coalesced memory access

## Common Issues

### Particle Clustering
- Use adaptive particle splitting/merging
- Implement super-particles
- Load balance regularly

### Numerical Heating
- Use higher-order interpolation
- Reduce timestep
- Check pusher accuracy

### Conservation
- Momentum not exactly conserved
- Use conservative schemes if critical
- Monitor conservation errors

## Testing

Test problems in `src/pgen/`:
- `part_random.cpp`: Random particles
- `part_static_turb.cpp`: Particles in turbulence

## See Also
- [Outputs Module](outputs.md) - Particle output formats
- [MHD Module](mhd.md) - Field coupling
- Source: `src/particles/particles.cpp`