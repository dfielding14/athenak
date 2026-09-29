# Module: DynGRMHD (Dynamical GR MHD)

## Role in AthenaK
DynGRMHD augments the standard MHD update with curvature source terms and
stress–energy feedback so that a dynamical spacetime (Z4c/ADM) and the fluid evolve
consistently. The module is instantiated only when a simulation enables both `<mhd>` and a
general-relativistic metric block; otherwise the standard MHD pipeline is used
([src/mesh/meshblock_pack.cpp:190]).

## Prerequisites
- `<z4c>` or `<adm>` **and** `<mhd>` must be present in the input deck; attempting to pair
  dynamical metrics with `<hydro>` is rejected ([src/mesh/meshblock_pack.cpp:206]).
- When DynGRMHD is enabled the mesh pack also owns a `Tmunu` accumulator that receives the
  fluid stress–energy tensor prior to the Z4c update ([src/mesh/meshblock_pack.cpp:214]).

## Source Layout
| File | Responsibility | Key APIs |
|------|----------------|----------|
| `dyn_grmhd.hpp` | Class interfaces, task IDs, policy enums | `DynGRMHD`, `DynGRMHDPS`, `BuildDynGRMHD` |
| `dyn_grmhd.cpp` | Object construction, task registration, primitive recovery, stress–energy | `QueueDynGRMHDTasks`, `ConToPrim`, `SetTmunu` |
| `dyn_grmhd_fluxes.cpp` | Reconstruction, GR Riemann solvers, FOFC trigger | `CalcFluxes<T>` |
| `dyn_grmhd_fofc.cpp` | First-order flux correction and excision handling | `FOFC<T>` |
| `dyn_grmhd_util.hpp` | Helper routines used by flux/FOFC paths | `ExtractPrimitives`, `InsertFluxes` |
| `dyn_grmhd/rsolvers/*` | Characteristic estimates for LLF/HLLE in curved spacetime | `LLF_DYNGR`, `HLLE_DYNGR` |

## Task Graph Integration
DynGRMHD does not create its own task lists; instead it extends the Numerical Relativity
queue so the fluid and spacetime stay synchronized. `QueueDynGRMHDTasks` wires the
standard MHD sequence while substituting GR-aware pieces where needed
([src/dyn_grmhd/dyn_grmhd.cpp:130]). The stage ordering is:

- Post/post-process communications with the existing `mhd::MHD` methods.
- Insert the GR-aware flux task (`CalcFluxes`) before Z4c updates the metric so fluxes see
  the latest spacetime ([src/dyn_grmhd/dyn_grmhd.cpp:145]).
- If a Z4c object exists, queue `SetTmunu` before the Runge–Kutta update so the metric
  driver has access to the stress–energy source computed from the current primitives
  ([src/dyn_grmhd/dyn_grmhd.cpp:160]).
- Run the usual MHD stages (RK update, CT, prolongation) and finally convert the updated
  conserved variables back to primitives using the GR primitive solver; the task is marked
  optional for excised cells via `Z4c_Excise` dependency
  ([src/dyn_grmhd/dyn_grmhd.cpp:191]).

Coordinate source terms are added from inside `MHD::MHDSrcTerms`. For dynamical metrics
the MHD source routine calls back into DynGRMHD to accumulate the metric-derivative terms
([src/mhd/mhd_tasks.cpp:251], [src/mhd/mhd_tasks.cpp:262]).

## Equation of State & Primitive Recovery
`BuildDynGRMHD` instantiates a `DynGRMHDPS` template that wraps
`PrimitiveSolverHydro<EOSPolicy, ErrorPolicy>` using the EOS named by `<mhd>/dyn_eos`
(`ideal`, `piecewise_poly`, or `compose`). The only supported error policy today is
`reset_floor`, meaning the primitive solver falls back to floor values when recovery
fails ([src/dyn_grmhd/dyn_grmhd.cpp:56]).

`DynGRMHDPS` implements the full conversion pipeline:
- `ConToPrim` runs the Valencia recovery on the entire pack each stage
  ([src/dyn_grmhd/dyn_grmhd.cpp:247]).
- `ConToPrimBC`/`PrimToConInit` are used to convert ghost layers before/after boundary
  conditions are applied ([src/dyn_grmhd/dyn_grmhd.cpp:267],
  [src/dyn_grmhd/dyn_grmhd.cpp:300]).
- `ConvertInternalEnergyToPressure` recomputes pressure from the temperature supplied by
  the EOS (important for CompOSE tables) ([src/dyn_grmhd/dyn_grmhd.cpp:220]).

Setting `<mhd>/fixed = true` skips the GR-specific updates—fluxes, C2P, coordinate source
terms, and stress–energy all short-circuit—so the run can operate in the Cowling
approximation without touching the spacetime ([src/dyn_grmhd/dyn_grmhd.cpp:123],
[src/dyn_grmhd/dyn_grmhd.cpp:248], [src/dyn_grmhd/dyn_grmhd.cpp:358]).

## Fluxes, Riemann Solvers, and FOFC
`CalcFluxes<T>` reuses the MHD reconstruction choice (DC/PLM/PPM/WENOZ) and swaps in GR
characteristic speeds when invoking either the LLF or HLLE solver selected by
`<mhd>/rsolver` ([src/dyn_grmhd/dyn_grmhd_fluxes.cpp:20]). The `dyn_scratch` integer tunes
the amount of shared memory reserved for these kernels
([src/dyn_grmhd/dyn_grmhd.cpp:120], [src/dyn_grmhd/dyn_grmhd_fluxes.cpp:58]).

After the high-order fluxes are computed, DynGRMHD optionally re-evaluates the step with a
first-order flux correction. FOFC is triggered when either `<mhd>/fofc = true` or the
coordinate system marks faces for black-hole excision
([src/dyn_grmhd/dyn_grmhd_fluxes.cpp:386]). The detection pass predicts the updated
conserved state, enforces a discrete maximum principle using `enforce_maximum` and `dmp_M`,
and attempts a trial primitive recovery to flag troubled cells
([src/dyn_grmhd/dyn_grmhd_fofc.cpp:63]). Flagged faces (or excised regions) are then
recomputed with a first-order LLF/HLLE flux using the same solver family as the primary
update ([src/dyn_grmhd/dyn_grmhd_fofc.cpp:167]).

`<mhd>/fofc_method` is parsed and stored for compatibility, but the current implementation
always reuses the main solver choice when rebuilding fluxes.

## Metric Coupling and Stress–Energy Feedback
- `AddCoordTermsEOS` evaluates ADM Christoffel symbols, lapse/shift gradients, and injects
  the curvature source terms into the conservative update. It is dispatched with the
  appropriate ghost-zone width from `MHD::MHDSrcTerms`
  ([src/dyn_grmhd/dyn_grmhd.cpp:422]).
- `SetTmunu` populates the energy density, momentum, and stress tensor that Z4c consumes in
  the subsequent stage. The routine is only scheduled when a Z4c object is present and is
  skipped entirely when `fixed` is enabled ([src/dyn_grmhd/dyn_grmhd.cpp:160],
  [src/dyn_grmhd/dyn_grmhd.cpp:353]).

## Configuration
All parameters live in the `<mhd>` block unless stated otherwise.

| Parameter | Default | Notes |
|-----------|---------|-------|
| `dyn_eos` | – | Required; choose `ideal`, `piecewise_poly`, or `compose` ([src/dyn_grmhd/dyn_grmhd.cpp:62]). |
| `dyn_error` | – | Must be `reset_floor`; other values abort at startup ([src/dyn_grmhd/dyn_grmhd.cpp:74]). |
| `rsolver` | – | `llf` or `hlle` for both the primary flux update and FOFC fallback ([src/dyn_grmhd/dyn_grmhd.cpp:97]). |
| `fofc_method` | `llf` | Parsed for future use; flux correction currently follows `rsolver` ([src/dyn_grmhd/dyn_grmhd.cpp:109]). |
| `dyn_scratch` | `0` | Controls the KokkoS scratch level used by `CalcFluxes` ([src/dyn_grmhd/dyn_grmhd.cpp:120]). |
| `enforce_maximum` | `true` | Enables the discrete maximum principle inside FOFC ([src/dyn_grmhd/dyn_grmhd.cpp:121]). |
| `dmp_M` | `1.2` | Multiplier used when evaluating the maximum-principle bounds ([src/dyn_grmhd/dyn_grmhd.cpp:122]). |
| `fixed` | `false` | Freezes MHD updates for Cowling runs ([src/dyn_grmhd/dyn_grmhd.cpp:124]). |
| `fofc` | `false` | Global FOFC toggle inherited from the base MHD module ([src/mhd/mhd.cpp:176]). |

## Notes and Limitations
- Only the `reset_floor` primitive-recovery policy is wired in at present; alternative
  policies would require new template instantiations.
- `ApplyPhysicalBCs` exists for completeness but the default task wiring continues to use
  the standard MHD boundary routines ([src/dyn_grmhd/dyn_grmhd.cpp:187]).
- When `fixed = true`, no stress–energy tensor is written back to Z4c, so the spacetime
  evolves independently of the MHD state ([src/dyn_grmhd/dyn_grmhd.cpp:358]).

## See Also
- [MHD Module](mhd.md) – array ownership, reconstruction choices, and the base FOFC flag.
- [Z4c Module](z4c.md) – evolution of the spacetime that consumes `T^{\mu\nu}`.
- [TaskList Module](tasklist.md) – details on how the Numerical Relativity queue executes.
