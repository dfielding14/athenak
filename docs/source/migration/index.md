# Athena++ To AthenaK Migration

AthenaK is not an API-compatible replacement for Athena++. A port should begin
from one implemented AthenaK problem or module path and be validated on the
public source baseline rather than inferred from an Athena++ input file.

## Key Changes

| Concern | AthenaK public pattern |
| --- | --- |
| Storage | Kokkos views such as `DvceArray5D<Real>` defined in `src/athena.hpp` |
| Cell loops | `par_for(..., KOKKOS_LAMBDA(...))` on `DevExeSpace()` |
| Work unit | Arrays span a `MeshBlockPack`; field access includes block index `m` |
| Initial conditions | `ProblemGenerator::UserProblem(ParameterInput *, const bool restart)` |
| Custom generator selection | Configure with `-DPROBLEM=<source_stem>` unless the built-in dispatcher lists its `pgen_name` |
| Refinement configuration | `<mesh_refinement>` and `<amr_criterion*>` blocks |

## Port A Problem Generator

1. Choose a nearby public example and copy its build/run validation path.
2. Create or edit a generator under `src/pgen/` implementing
   `ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart)`.
3. Obtain pack and mesh metadata as production generators do:

```cpp
MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
auto &indcs = pmy_mesh_->mb_indcs;
int nmb = pmbp->nmb_thispack;
auto &u0 = pmbp->phydro->u0;

par_for("my_problem", DevExeSpace(), 0, nmb - 1,
        indcs.ks, indcs.ke, indcs.js, indcs.je, indcs.is, indcs.ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    u0(m, IDN, k, j, i) = 1.0;
  });
```

4. Configure a problem-specific executable and run a minimal deck:

```bash
cmake -S . -B build-my-problem -DPROBLEM=my_problem
cmake --build build-my-problem
./build-my-problem/src/athena -i inputs/path/to/my_problem.athinput -d run-my-problem time/nlim=1
```

Use [Problem Generators](../modules/pgen.md) for built-in versus custom
selection rules and [Kokkos Guide](../kokkos_guide.md) for device-accessible
data patterns.

## Port Inputs Deliberately

Do not assume an Athena++ block or parameter is accepted. Begin with a shipped
AthenaK input matching the active physics, then apply source-backed changes:

- Set mesh and time controls using [Configuration](../configuration.md).
- Select fluid reconstruction and solvers using [Hydro](../modules/hydro.md)
  or [MHD](../modules/mhd.md).
- Configure refinement through `<mesh_refinement>` and
  `<amr_criterion*>` as documented in [Mesh](../modules/mesh.md).
- Configure output streams using [Outputs](../modules/outputs.md).

## Capability Boundaries

Only interfaces documented in stable module pages should be treated as public
porting targets. For example, the public particle module currently exposes a
cosmic-ray `drift` pusher; pages under [Developer Notes](../engineering/index.md)
that discuss alternative particle work are not stable migration instructions.

## Validation

For each port:

1. Build the selected generator and backend.
2. Run one or a few cycles with output enabled.
3. Check output names/content and diagnostic behavior.
4. Compare against an expected solution or a predecessor run when making a
   scientific-equivalence claim.

See [Common Gotchas](common_gotchas.md) for recurring implementation traps.
