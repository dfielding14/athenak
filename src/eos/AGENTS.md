<!-- BEGIN build-memory-table -->
# Equation of State Navigation

## Scope

This subtree owns EOS parameters, characteristic-speed helpers, floors and ceilings, and
conserved/primitive conversion for hydro and MHD. Fluid constructors select the concrete
EOS; flux, boundary, source-term, and timestep code consume `EOS_Data`. Dynamical GRMHD
uses the separate policy-based primitive solver under `primitive-solver/`.

## Task lookup

| Task | Start in | Also inspect or validate |
| --- | --- | --- |
| Change the common interface, floors, or wave-speed helpers | `src/eos/eos.hpp`, `src/eos/eos.cpp` | Hydro/MHD fluxes, FOFC, event counters |
| Change non-relativistic ideal/isothermal conversion | `src/eos/ideal_hyd.cpp`, `ideal_mhd.cpp`, `isothermal_*.cpp` | `src/hydro/hydro.cpp`, `src/mhd/mhd.cpp` |
| Change special- or general-relativistic conversion | `src/eos/ideal_sr*.cpp`, `ideal_gr*.cpp`, `ideal_c2p_*.hpp` | SR/GR linear-wave and shock-tube tests |
| Change EOS selection | `src/hydro/hydro.cpp`, `src/mhd/mhd.cpp` | Coordinate mode flags and input examples |
| Change dynamical-GRMHD primitive recovery | `src/eos/primitive-solver/`, `src/dyn_grmhd/dyn_grmhd.cpp` | `BuildDynGRMHD`, `src/z4c/AGENTS.md` |
| Change tabulated, hybrid, or Compose policies | `src/eos/primitive-solver/` | EOS policy headers/sources, unit systems, error policies, dynamical-GRMHD tests |

## Important flow

- `Hydro` or `MHD` reads its input block and coordinate mode, constructs one
  `EquationOfState` subclass, and exposes its compact `EOS_Data` to device kernels.
- Reconstruction and Riemann solvers use EOS characteristic speeds. After each update and
  boundary fill, `ConsToPrim` enforces physical limits and refreshes primitive state.
- First-order flux correction calls conversion in a test-only-floor mode to detect bad
  candidate states before replacing fluxes.
- Dynamical GRMHD builds a policy-based primitive solver and error response in
  `primitive-solver/`; it is not interchangeable with the legacy virtual overloads.

## Local validation

- `cd tst && python3 run_test_suite.py --test test_suite/sr/test_sr_lwave1d_cpu.py`:
  SR hydro/MHD conversion and wave-speed convergence.
- `cd tst && python3 run_test_suite.py --test test_suite/gr/test_gr_lwave1d_cpu.py`:
  stationary-spacetime GR conversion.
- `cd tst && python3 run_test_suite.py --test test_suite/dyngrmhd/test_dyngrmhd_nqt_shocktube_cpu.py`:
  policy-based dynamical-GRMHD primitive recovery.

## Local constraints

- Hydro and MHD conversions are distinct virtual overloads; a derived class normally
  implements only the overload matching its fluid.
- If floors modify a primitive state, audit whether conserved variables and
  `Mesh::ecounter` must also change. Relativistic failure responses are part of solver
  behavior, not merely diagnostics.
- A new compiled EOS implementation must be added to `src/CMakeLists.txt`, wired into the
  appropriate fluid constructor, and covered by both conversion and evolution tests.
<!-- END build-memory-table -->
