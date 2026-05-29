# Developer Interface Reference

AthenaK does not publish a separately versioned external C++ API. This page
maps the public headers and supported extension points in the documented
source tree; the headers remain authoritative when implementation details
change.

## Core Types

| Type | Declared in | Confirmed public role |
| --- | --- | --- |
| `Mesh` | `src/mesh/mesh.hpp` | Root mesh, MeshBlock metadata, time state, refinement, physics attachment |
| `MeshBlock` | `src/mesh/meshblock.hpp` | Per-pack block IDs, refinement levels, physical extents, boundary/neighbor data |
| `MeshBlockPack` | `src/mesh/meshblock_pack.hpp` | Pack ownership of coordinates, optional physics objects, and task lists |
| `Driver` | `src/driver/driver.hpp` | Initialization, task-list execution, time evolution, finalization |
| `ProblemGenerator` | `src/pgen/pgen.hpp` | Initial conditions and optional user callbacks |
| `hydro::Hydro` | `src/hydro/hydro.hpp` | Hydro state, solver selection, stage tasks |
| `mhd::MHD` | `src/mhd/mhd.hpp` | MHD/CT state, solver selection, stage tasks |
| `radiation::Radiation` | `src/radiation/radiation.hpp` | Angle-resolved radiation state and tasks |

## Header-Level Interfaces

`Mesh` exposes construction and runtime operations including:

```cpp
explicit Mesh(ParameterInput *pin);
void BuildTreeFromScratch(ParameterInput *pin);
void BuildTreeFromRestart(ParameterInput *pin, IOWrapper &resfile,
                          bool single_file_per_rank=false);
void AddCoordinatesAndPhysics(ParameterInput *pinput);
void NewTimeStep(const Real tlim);
BoundaryFlag GetBoundaryFlag(const std::string& input_string);
```

`Driver` uses the following public lifecycle calls:

```cpp
Driver(ParameterInput *pin, Mesh *pmesh, Real wtlim, Kokkos::Timer* ptimer);
void Initialize(Mesh *pmesh, ParameterInput *pin, Outputs *pout, bool rflag);
void Execute(Mesh *pmesh, ParameterInput *pin, Outputs *pout);
void Finalize(Mesh *pmesh, ParameterInput *pin, Outputs *pout);
void ExecuteTaskList(Mesh *pm, std::string tl, int stage);
```

`Hydro` and `MHD` expose their task assembly and per-stage task methods in
their headers. For example, both provide `Fluxes`, `RKUpdate`,
`ApplyPhysicalBCs`, `Prolongate`, `ConToPrim`, and `NewTimeStep`; the MHD
class additionally provides `CornerE` and `CT`.

## Kokkos Storage Aliases

`src/athena.hpp` defines array aliases as templates over the element type:

```cpp
using DevExeSpace = Kokkos::DefaultExecutionSpace;
using DevMemSpace = Kokkos::DefaultExecutionSpace::memory_space;
using LayoutWrapper = Kokkos::LayoutRight;

template <typename T>
using DvceArray4D = Kokkos::View<T ****, LayoutWrapper, DevMemSpace>;
template <typename T>
using DvceArray5D = Kokkos::View<T *****, LayoutWrapper, DevMemSpace>;
template <typename T>
using DvceArray6D = Kokkos::View<T ******, LayoutWrapper, DevMemSpace>;
```

The last index is fastest varying under `LayoutRight`. Consult
[Kokkos Guide](../kokkos_guide.md) before adding device kernels.

## Problem-Generator Extension Points

A custom generator selected at configure time with `-DPROBLEM=<stem>`
implements:

```cpp
void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart);
```

`src/pgen/pgen.hpp` exposes these optional callback members:

| Callback member | Use |
| --- | --- |
| `pgen_final_func` | Work executed during driver finalization |
| `user_bcs_func` | Required by any boundary configured with `user` |
| `user_srcs_func` | Used when `<problem>/user_srcs = true` |
| `user_ref_func` | Required when an adaptive criterion has `method = user` |
| `user_hist_func` | Used when `<problem>/user_hist = true` |

See [Problem Generators](../modules/pgen.md) for the built-in/custom selection
rules and shipped examples. `method = user` is a strict callback contract:
the adaptive-refinement path invokes the enrolled refinement function.

## Boundaries And Parameters

`Mesh::GetBoundaryFlag()` accepts `periodic`, `reflect`, `inflow`,
`outflow`, `diode`, `user`, `shear_periodic`, `vacuum`, and `undef`.
`undef` is also used internally for inactive directions; ordinary input decks
specify physical boundary flags only for active directions. `shear_periodic`
has additional mesh constraints described in [Mesh](../modules/mesh.md).

`ParameterInput` provides the typed `Get*` and `GetOrAdd*` interfaces used by
module constructors. For runnable public parameters and defaults, use
[Public Input Parameters](input_parameters.md) and inspect the selected
module source or shipped input deck for specialized controls.

## Validation Guidance

When adding or changing an interface:

1. Match declarations in the owning header and parser logic in the
   implementation.
2. Add or update a shipped input/test when the behavior is intended to be
   public.
3. Build and execute a minimal input that exercises the new path.
4. Update the relevant module and parameter-reference pages with only the
   validated behavior.
