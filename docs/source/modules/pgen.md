# Module: Problem Generators

Problem generators initialize simulation state and may enroll
problem-specific boundary, source, refinement, history, or finalization
callbacks. The public interface is implemented in `src/pgen/pgen.hpp` and
`src/pgen/pgen.cpp`.

## Selection Model

AthenaK supports two public selection paths:

| Path | How it is selected | Implementation |
| --- | --- | --- |
| Built-in verification generator | Set `<problem>/pgen_name` in an input deck | Methods compiled from `src/pgen/tests/*.cpp` and dispatched in `ProblemGenerator::CallProblemGenerator()` |
| Custom/source generator | Configure with `-DPROBLEM=<file-stem>` | The chosen `src/pgen/<file-stem>.cpp` implements `ProblemGenerator::UserProblem()` |

The default build accepts these `pgen_name` values:

| `pgen_name` | Built-in method |
| --- | --- |
| `advection` | `Advection` |
| `cpaw` | `AlfvenWave` |
| `gr_bondi` | `BondiAccretion` |
| `cshock` | `CShock` |
| `linear_wave` | `LinearWave` |
| `implode` | `LWImplode` |
| `gr_monopole` | `Monopole` |
| `mri3d` | `MRI3d` |
| `orszag_tang` | `OrszagTang` |
| `rad_linear_wave` | `RadiationLinearWave` |
| `rad_beam` | `RadiationBeam` |
| `shock_tube` | `ShockTube` |
| `shwave` | `Shwave` |
| `z4c_boosted_puncture` | `Z4cBoostedPuncture` |
| `z4c_linear_wave` | `Z4cLinearWave` |
| `spherical_collapse` | `SphericalCollapse` |
| `diffusion` | `Diffusion` |

For example, `inputs/hydro/sod.athinput` selects `pgen_name = shock_tube`;
`inputs/mhd/orszag_tang.athinput` selects `pgen_name = orszag_tang`.

## Custom Problem Generator

Custom files implement the class member declared in `src/pgen/pgen.hpp`:

```cpp
void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (restart) return;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  // Initialize the active module data on pmbp.
}
```

Build and run a shipped custom problem with:

```bash
cmake -S . -B build-turb -DPROBLEM=turb
cmake --build build-turb
./build-turb/src/athena -i inputs/hydro/turb.athinput -d run-turb time/nlim=1
```

Other public custom source files include `blast.cpp`, `field_loop.cpp`,
`gr_torus.cpp`, `kh.cpp`, `rt.cpp`, and radiation or numerical-relativity
setups under `src/pgen/`. A matching input deck is necessary; the source file
is authoritative for its `<problem>` parameters.

## Optional Callbacks

The `ProblemGenerator` object exposes function pointers for optional hooks:

| Hook member | Enabled by |
| --- | --- |
| `user_bcs_func` | A mesh boundary configured as `user` |
| `user_srcs_func` | `<problem>/user_srcs = true` |
| `user_hist_func` | `<problem>/user_hist = true` |
| `user_ref_func` | Mandatory when an adaptive `<amr_criterion*>` selects `method = user` |
| `pgen_final_func` | Set by a generator needing final analysis |

When `user_srcs` or `user_hist` is requested but its callback was not set by
the problem generator, construction aborts. The adaptive `method = user` path
invokes `user_ref_func`, so its generator must enroll the callback before
refinement checks run. A generator using `inflow` boundaries must initialize
the associated constant inflow states.

## Development Rule

Do not infer public support from a development-only problem-generator file or
draft input deck. A public worked example requires the source file, its input
deck, and a verified build/run command to exist on the public implementation
baseline.

## See Also

- [Configuration](../configuration.md)
- [Worked Examples](../examples/index.md)
- [Input Parameters](../reference/input_parameters.md)
