# Module: DynGRMHD (Dynamical GR MHD)

## Role in AthenaK
DynGRMHD augments the standard MHD update with curvature source terms and
primitive recovery when a simulation enables both `<mhd>` and an `<adm>` or
`<z4c>` metric block; otherwise the standard MHD pipeline is used. With Z4c,
the module also supplies staged stress-energy feedback to the evolving
spacetime. An ADM-only configuration uses the DynGRMHD fluid path on an ADM
metric, but does not use the Z4c feedback task route
([src/mesh/meshblock_pack.cpp:190], [src/dyn_grmhd/dyn_grmhd.cpp:196]).

## Prerequisites
- `<z4c>` or `<adm>` **and** `<mhd>` must be present in the input deck; attempting to pair
  dynamical metrics with `<hydro>` is rejected ([src/mesh/meshblock_pack.cpp:206]).
- With Z4c, the mesh pack's `Tmunu` accumulator receives the fluid
  stress-energy tensor before the Z4c update ([src/mesh/meshblock_pack.cpp:214],
  [src/dyn_grmhd/dyn_grmhd.cpp:196]).

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
(`ideal`, `piecewise_poly`, `compose`, or `hybrid`). The only supported error policy today is
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

Setting `<mhd>/fixed = true` skips staged fluxes, C2P, coordinate source terms,
and scheduled `SetTmunu` refreshes. `PrimToConInit` still initializes `Tmunu`
with the fixed guard temporarily disabled, so a Z4c run can retain that
initial matter source. This is not conventional Cowling matter evolution
([src/dyn_grmhd/dyn_grmhd.cpp:249], [src/dyn_grmhd/dyn_grmhd.cpp:411]).

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
recomputed with a first-order fallback flux ([src/dyn_grmhd/dyn_grmhd_fofc.cpp:167]).

`<mhd>/fofc_method` is parsed and stored for compatibility, but the current implementation
does not uniformly follow either requested setting. In particular, the current `hlle`
path includes an x1 right-face LLF fallback call, so it must not be assumed to reuse the
primary solver family uniformly.

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
| `eos` | – | Must be `ideal` so the current MHD state layout includes the energy storage read by DynGRMHD. |
| `dyn_eos` | – | Required; choose `ideal`, `piecewise_poly`, `compose`, or `hybrid` ([src/dyn_grmhd/dyn_grmhd.cpp:78]). |
| `dyn_error` | – | Must be `reset_floor`; other values abort at startup ([src/dyn_grmhd/dyn_grmhd.cpp:74]). |
| `rsolver` | – | `llf` or `hlle` for the primary flux update; fallback behavior is qualified above ([src/dyn_grmhd/dyn_grmhd.cpp:97]). |
| `fofc_method` | `llf` | Parsed/stored; current fallback implementation does not uniformly follow either requested setting ([src/dyn_grmhd/dyn_grmhd.cpp:109]). |
| `dyn_scratch` | `0` | Controls the KokkoS scratch level used by `CalcFluxes` ([src/dyn_grmhd/dyn_grmhd.cpp:120]). |
| `enforce_maximum` | `true` | Enables the discrete maximum principle inside FOFC ([src/dyn_grmhd/dyn_grmhd.cpp:121]). |
| `dmp_M` | `1.2` | Multiplier used when evaluating the maximum-principle bounds ([src/dyn_grmhd/dyn_grmhd.cpp:122]). |
| `fixed` | `false` | Disables staged fluid-side dynamic updates and staged `Tmunu` refreshes after initialization ([src/dyn_grmhd/dyn_grmhd.cpp:124]). |
| `fofc` | `false` | Global FOFC toggle inherited from the base MHD module ([src/mhd/mhd.cpp:176]). |

## Notes and Limitations
- Only the `reset_floor` primitive-recovery policy is wired in at present; alternative
  policies would require new template instantiations.
- `ApplyPhysicalBCs` exists for completeness but the default task wiring continues to use
  the standard MHD boundary routines ([src/dyn_grmhd/dyn_grmhd.cpp:187]).
- When `fixed = true`, initialization can populate `Tmunu`, but later staged
  fluid updates and `SetTmunu` refreshes are suppressed; a Z4c evolution may
  therefore retain the initial matter source ([src/dyn_grmhd/dyn_grmhd.cpp:249],
  [src/dyn_grmhd/dyn_grmhd.cpp:411]).

## See Also
- [MHD Module](mhd.md) – array ownership, reconstruction choices, and the base FOFC flag.
- [Z4c Module](z4c.md) – evolution of the spacetime that consumes `T^{\mu\nu}`.
- [TaskList Module](tasklist.md) – details on how the Numerical Relativity queue executes.
