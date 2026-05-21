# Self-Gravity and Multigrid

Self-gravity solves

```text
del^2 phi = 4 pi G rho
```

with a multigrid Poisson solve, stores the potential in
`MeshBlockPack::pgrav->phi`, and applies `-grad(phi)` through the hydro or MHD
source-term path.

## Input Contract

Use one fluid block and one matching source-term block:

```ini
<hydro>
eos = isothermal

<gravity>
four_pi_G = 1.0
threshold = 1.0e-8
full_multigrid = true
mg_bc = none

<hydro_srcterms>
self_gravity = true
```

For MHD, use `<mhd>` and `<mhd_srcterms>` instead. The current implementation
rejects ambiguous two-fluid self-gravity inputs that include both `<hydro>` and
`<mhd>`.

Startup validation is intentionally strict. Self-gravity currently requires:

- a `<gravity>` block whenever `self_gravity = true` is set;
- a 3D mesh;
- logically cubic MeshBlocks: `meshblock/nx1 == meshblock/nx2 == meshblock/nx3`;
- MeshBlock sizes that are powers of two;
- `gravity/four_pi_G >= 0`;
- either `gravity/threshold` or `gravity/niteration`;
- `gravity/niteration > 0` when `gravity/threshold < 0`;
- valid `gravity/mg_bc`: `none`, `zerofixed`, `zerograd`, or `multipole`;
- valid multipole options when `mg_bc = multipole`: `mporder = 2` or `4`,
  and not both `auto_mporigin = true` and `nodipole = true`.

## Minimal Periodic Self-Gravity Input

Periodic gravity removes the mean density before solving because a nonzero mean
source has no periodic potential. Use `mg_bc = none` and leave
`subtract_average = true`.

```ini
<mesh>
nghost = 3
nx1 = 32
nx2 = 32
nx3 = 32
ix1_bc = periodic
ox1_bc = periodic
ix2_bc = periodic
ox2_bc = periodic
ix3_bc = periodic
ox3_bc = periodic

<meshblock>
nx1 = 16
nx2 = 16
nx3 = 16

<gravity>
four_pi_G = 1.0
threshold = 1.0e-8
full_multigrid = true
mg_nghost = 1
mg_bc = none
subtract_average = true
```

Complete examples:

- `inputs/selfgravity/periodic_jeans_hydro.athinput`
- `inputs/selfgravity/periodic_jeans_mhd.athinput`

## Isolated Multipole Input

For isolated problems, use non-periodic mesh boundaries and override the
multigrid boundary with `mg_bc = multipole`.

```ini
<gravity>
four_pi_G = 1.0
threshold = 1.0e-6
full_multigrid = true
mg_bc = multipole
mporder = 4
auto_mporigin = true
nodipole = false
```

If `auto_mporigin = false`, set all three origin coordinates:
`mporigin_x1`, `mporigin_x2`, and `mporigin_x3`.

Complete examples:

- `inputs/selfgravity/isolated_binary_multipole.athinput`
- `inputs/selfgravity/be_collapse_hydro_smr.athinput`
- `inputs/selfgravity/be_collapse_mhd.athinput`

## Choosing Convergence Controls

Use threshold mode for production and correctness tests:

```ini
threshold = 1.0e-8
niteration = 1
```

When `threshold >= 0`, `niteration` is only a fallback/default parameter in the
input and the solve iterates to the threshold. For fixed work per solve, set:

```ini
threshold = -1.0
niteration = 4
```

Recommended starting points:

| Problem | Suggested control |
|---------|-------------------|
| Periodic Jeans tests | `threshold = 1.0e-8` |
| Isolated binary or BE smoke | `threshold = 1.0e-6` |
| Cheap performance scans | `threshold = -1`, `niteration = 4` |
| Regression comparisons | `show_defect = true` |

## Performance Knobs

The main knobs are:

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `npresmooth` | `1` | Pre-smoothing sweeps per V-cycle. |
| `npostsmooth` | `1` | Post-smoothing sweeps per V-cycle. |
| `full_multigrid` | `false` | Start from FMG instead of the existing potential. |
| `fmg_ncycle` | `1` | V-cycles per FMG level. |
| `mg_nghost` | `1` | Multigrid ghost depth. |
| `root_on_host` | `false` | Store root-grid solve arrays on host memory. |
| `profile` | `false` | Print phase timing diagnostics per gravity solve. |

With `profile = true`, each solve prints source loading, setup, root transfer,
smoothing, boundary exchange, restriction/prolongation, result retrieval, and
total solve timings. Use this before changing smoother counts or root placement.

## Output

Write the potential with:

```ini
<output1>
file_type = bin
variable = grav_phi
dt = 1.0
```

For source-term checks, write conserved fluid variables:

```ini
<output2>
file_type = bin
variable = hydro_u
dt = 1.0
```

## Tested Configurations

The regression suite includes:

- periodic Jeans hydro and MHD convergence;
- analytic periodic Jeans potential comparison;
- hydro/MHD zero-field gravity equivalence;
- binary uniform-sphere potential and force comparison with multipole BCs;
- source-term momentum update comparison against the analytic Jeans force;
- restart equivalence for two cycles run straight versus one cycle plus restart;
- `root_on_host = false` and `root_on_host = true` comparisons;
- fixed-iteration mode;
- startup failure tests for invalid self-gravity inputs;
- static refined-hierarchy BE collapse smoke;
- MPI 2-rank and 4-rank Jeans tests;
- MPI multipole binary reduction test;
- GPU smoke tests for Jeans and root-host/device equivalence when a CUDA build is
  available.

Run the serial CPU suite with:

```bash
cd tst
python run_test_suite.py --cpu --test test_suite/selfgravity/test_selfgravity_cpu.py
```

Run the MPI suite with:

```bash
cd tst
python run_test_suite.py --mpicpu --test test_suite/selfgravity/test_selfgravity_mpicpu.py
```

Run the GPU suite on a CUDA-capable node with:

```bash
cd tst
python run_test_suite.py --gpu --test test_suite/selfgravity/test_selfgravity_gpu.py
```

## Debugging Convergence Failures

Use this order:

1. Set `show_defect = true` and confirm the defect decreases every cycle.
2. Set `profile = true` and check whether time is going into smoothing,
   boundary exchange, root transfer, or setup.
3. For periodic problems, keep `subtract_average = true`.
4. For isolated problems, use `mg_bc = multipole` before tuning smoother counts.
5. If fixed-iteration mode is used, compare against threshold mode before
   trusting the result.
6. If a restart diverges from a straight run, compare `grav_phi` and the fluid
   conserved variables at the same cycle.

## Current Limits

- Self-gravity is 3D only.
- MeshBlocks must be logically cubic and power-of-two sized.
- The active source is mass density only; magnetic fields do not directly source
  the potential.
- Exactly one fluid module is supported per self-gravity solve.
- Periodic solves remove the mean source term.
- The tested multilevel path is a static refined hierarchy. Dynamically moving
  AMR with self-gravity should be treated as experimental until it has its own
  refinement/regrid regression.
