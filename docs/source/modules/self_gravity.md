# Self-Gravity and Multigrid

Self-gravity solves the Poisson equation

```text
del^2 phi = 4 pi G rho
```

with a multigrid solver, stores the gravitational potential in
`MeshBlockPack::pgrav->phi`, and applies the corresponding acceleration through
the hydro or MHD source-term path.

## Enabling Self-Gravity

Self-gravity requires three pieces in the input file:

```ini
<gravity>
four_pi_G       = 1.0
threshold       = 1.0e-8
full_multigrid  = true

<hydro_srcterms>
self_gravity = true
```

For MHD, use the same `<gravity>` block and enable the source term under
`<mhd_srcterms>`:

```ini
<mhd_srcterms>
self_gravity = true
```

The `<gravity>` block allocates the potential and multigrid driver. The
`self_gravity` source-term flag controls whether the solved potential is used to
update momentum and, for ideal equations of state, energy.

## Solver Controls

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `four_pi_G` | required | Gravitational coupling in code units. |
| `threshold` | required unless `niteration` is set | Relative defect norm used for convergence. |
| `niteration` | required unless `threshold` is set | Fixed multigrid iteration count. Ignored when `threshold >= 0`. |
| `omega` | `1.15` | Weighted Jacobi smoothing parameter. |
| `npresmooth` | `1` | Number of pre-smoothing sweeps per V-cycle. |
| `npostsmooth` | `1` | Number of post-smoothing sweeps per V-cycle. |
| `full_multigrid` | `false` | Start from a full multigrid solve instead of only iterating from the current potential. |
| `fmg_ncycle` | `1` | V-cycles per level in the FMG setup. |
| `mg_nghost` | `1` | Ghost-zone depth used by the multigrid hierarchy. |
| `root_on_host` | `false` | Store the root multigrid solve on host memory. |
| `show_defect` | `false` | Print the initial, iterative, and final defect norms. |
| `subtract_average` | `true` for periodic domains | Remove the mean source term before solving the periodic Poisson problem. |

For periodic domains the mean density cannot source a periodic potential, so
`subtract_average = true` should usually remain enabled. It is automatically
disabled for non-periodic or multipole-boundary solves.

## Boundary Conditions

By default, the multigrid boundary conditions follow the mesh boundary flags.
The `<gravity>/mg_bc` override can replace non-periodic multigrid boundaries:

| `mg_bc` | Meaning |
|---------|---------|
| `zerofixed` | Zero-value Dirichlet boundary. |
| `zerograd` | Zero-gradient Neumann boundary. |
| `multipole` | Isolated-boundary multipole expansion. |

Multipole boundaries support:

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `mporder` | `4` | Expansion order, either `2` or `4`. |
| `auto_mporigin` | `true` | Use the density center of mass as the expansion origin. |
| `nodipole` | `false` | Suppress the dipole term. Mutually exclusive with `auto_mporigin`. |
| `mporigin_x1`, `mporigin_x2`, `mporigin_x3` | required when `auto_mporigin = false` | User-set expansion origin. |

## Source-Term Update

The solver loads density from hydro or MHD conserved variables and solves for
`phi`. `SourceTerms::SelfGravity` then updates momentum using centered
differences of `phi`.

For ideal equations of state, the energy update uses the Godunov density fluxes
from the Riemann solve. The update is designed to be compatible with the
finite-volume flux path, but total energy conservation is still sensitive to the
finite multigrid residual.

## Output

The potential can be written with:

```ini
<output1>
file_type = bin
variable  = grav_phi
dt        = 1.0
```

`grav_phi` requires an active `<gravity>` block.

## Test Problems

The feature adds three built-in test problem generators:

| `pgen_name` | Source | Purpose |
|-------------|--------|---------|
| `gravity` | `src/pgen/tests/jeans_wave.cpp` | Periodic Jeans-wave setup for hydro or MHD. |
| `binary_gravity` | `src/pgen/tests/binary_gravity.cpp` | Two dense spheres for isolated-potential validation. |
| `be_collapse` | `src/pgen/tests/be_collapse.cpp` | Bonnor-Ebert-like collapse with Jeans refinement. |

The regression inputs are:

| Input | Coverage |
|-------|----------|
| `inputs/tests/selfgravity.athinput` | Hydro Jeans wave, multiple MeshBlocks, `grav_phi` output. |
| `inputs/tests/selfgravity_mhd.athinput` | MHD Jeans wave with the same density source. |

Run the CPU regression tests with:

```bash
cd tst
python run_test_suite.py --cpu --test test_suite/selfgravity/test_selfgravity_cpu.py
```

## Current Limits

- MeshBlocks must be logically cubic for the multigrid solve.
- The Poisson source is mass density only; magnetic fields do not source the
  potential.
- Periodic gravity solves remove the mean source term.
- The currently tested path is serial CPU with multiple MeshBlocks. MPI, GPU,
  moving AMR, and isolated multipole accuracy need larger follow-up tests before
  production use.
