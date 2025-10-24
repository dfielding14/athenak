# Module: Ion-Neutral

## Role in AthenaK
The Ion-Neutral module extends AthenaK’s MHD solver with a second, neutral fluid and
couples the two through implicit drag/chemistry terms. It assumes an existing single-fluid
MHD state for the ion component, augments it with a hydro module for the neutrals, and
adds the stiff coupling pieces needed for partially-ionised plasmas
([src/ion-neutral/ion-neutral.hpp:32]).

## Activation Requirements
- `<ion-neutral>` block must appear **alongside** both `<mhd>` and `<hydro>`; the mesh
  constructor checks for this combination and aborts if either fluid block is missing
  ([src/mesh/meshblock_pack.cpp:206]).
- The driver must run an ImEx integrator. At start-up the driver verifies `nimp_stages >
  0` when ion-neutral is active and terminates otherwise
  ([src/driver/driver.cpp:340]).

## Source Layout
| File | Responsibility | Key APIs |
|------|----------------|----------|
| `ion-neutral.hpp` | Class definition, task IDs, coupling coefficients | `IonNeutral`, `FirstTwoImpRK`, `ImpRKUpdate` |
| `ion-neutral.cpp` | Constructor that reads user coefficients | `IonNeutral::IonNeutral` |
| `ion-neutral_tasks.cpp` | Task graph wiring and implicit update kernels | `AssembleIonNeutralTasks`, `ImpRKUpdate` |

## Task Graph Integration
`AssembleIonNeutralTasks` interleaves the neutral hydro tasks with the ion MHD ones inside
the standard stage lists ([src/ion-neutral/ion-neutral_tasks.cpp:20]). The ordering is:

1. `FirstTwoImpRK` runs before any explicit work and performs two implicit solves for the
   stiff drag term, while copying both fluids’ conserved states into the ImEx scratch
   registers ([src/ion-neutral/ion-neutral_tasks.cpp:42],
   [src/ion-neutral/ion-neutral_tasks.cpp:58]).
2. The ion (MHD) stage advances first, reusing the existing flux/communication tasks.
3. The neutral (hydro) stage follows immediately, using the updated ion momentum in the
   coupling source.
4. `ImpRKUpdate` executes after both explicit legs to finish the implicit solve for the
   current stage and populate the driver’s `impl_src` array for subsequent stages
   ([src/ion-neutral/ion-neutral_tasks.cpp:68],
   [src/ion-neutral/ion-neutral_tasks.cpp:128]).

Both fluids continue through the usual boundary-condition, prolongation, and C2P steps,
and each reports a timestep limit via `NewTimeStep` so the driver can take the minimum.

## Implicit Drag and Chemistry
The module applies collisional drag, ionisation, and recombination via ImEx updates:

- `drag_coeff` (`γ`) multiplies the symmetric ion/neutral momentum exchange term
  ([src/ion-neutral/ion-neutral.hpp:45], [src/ion-neutral/ion-neutral_tasks.cpp:206]).
- `ionization_coeff` (`ξ`) adds mass/momentum to the ion fluid and removes it from the
  neutrals ([src/ion-neutral/ion-neutral_tasks.cpp:220]).
- `recombination_coeff` (`α`) removes ions and returns mass to the neutrals
  ([src/ion-neutral/ion-neutral_tasks.cpp:224]).

The analytic solve in `ImpRKUpdate` blends the two fluids’ momentum using these
coefficients and updates the densities accordingly before recomputing the stiff source
terms for later stages ([src/ion-neutral/ion-neutral_tasks.cpp:172],
[src/ion-neutral/ion-neutral_tasks.cpp:232]).

## Input Parameters (`<ion-neutral>` block)
| Parameter | Default | Description |
|-----------|---------|-------------|
| `drag_coeff` | – (required) | Collisional coupling strength between ions and neutrals ([src/ion-neutral/ion-neutral.cpp:18]). |
| `ionization_coeff` | `0.0` | Converts neutral mass/momentum into ion mass/momentum ([src/ion-neutral/ion-neutral.cpp:19]). |
| `recombination_coeff` | `0.0` | Converts ion mass/momentum back to neutrals ([src/ion-neutral/ion-neutral.cpp:20]). |

## Notes
- The module manipulates the existing `hydro::Hydro` and `mhd::MHD` arrays in-place and
  therefore inherits their reconstruction, Riemann solver, and boundary settings.
- Implicit solves run only on stages `< nexp_stages`; the final stage reuses the stored
  source terms without recomputing them ([src/ion-neutral/ion-neutral_tasks.cpp:159]).
- Because the stiff source relies on the driver’s ImEx scratch (`Driver::impl_src`), runs
  without an implicit leg will abort during initialization.

## See Also
- [MHD Module](mhd.md) – ion fluid equations and EM coupling.
- [Hydrodynamics Module](hydro.md) – neutral fluid infrastructure.
- [Driver Module](driver.md) – details on ImEx integrator requirements.
