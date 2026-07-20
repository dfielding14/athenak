# Mignone-R2 Full-Hall Shock: Simulation and Injection Setup

This document specifies the current uniform-grid 2D3V shock calculation. It is
a close reproduction of the non-relativistic R2 experiment in Mignone et al.
(2018), with AthenaK's full CR-Hall MHD-PIC coupling. Two minimal safeguards
are added to the injection prescription: only swept tracer in gas with
\(p>30p_0\) is eligible, and a carrier injects only when its complete local
particle bundle can be subtracted conservatively while retaining a thermal
reserve.

The authoritative runtime deck is
[`inputs/publication/pic_parallel_shock_mignone_r2_full_hall.athinput`](inputs/publication/pic_parallel_shock_mignone_r2_full_hall.athinput).
The shock and injection implementation is
[`src/pgen/tests/pic_parallel_shock.cpp`](src/pgen/tests/pic_parallel_shock.cpp).

## Simulation setup

The normalization follows Mignone et al. (2018): velocity is measured in the
upstream Alfvén speed \(v_{A,0}\), length in the ion skin depth
\(c/\omega_{pi}=v_{A,0}/\Omega_0\), and time in \(\Omega_0^{-1}\). AthenaK
uses \(\rho_0=p_0=B_0=1\), so \(v_{A,0}=1\), with background-ion
\(q/(mc)=1\).

| Quantity | Current value |
| --- | --- |
| Geometry | Uniform 2D in \((x_1,x_2)\), with three particle velocity components |
| Domain | \(0\leq x_1\leq120000\), \(0\leq x_2\leq3000\), unit-depth \(x_3\) |
| Grid | \(11520\times288\times1\); \(\Delta x_1=\Delta x_2=10.4167\,c/\omega_{pi}\) |
| MeshBlocks | \(64\times32\times1\); no mesh refinement |
| Boundaries | Reflecting wall at inner \(x_1\), steady inflow at outer \(x_1\), periodic \(x_2\) |
| Upstream gas | \(\rho_0=1\), \(p_0=1\), \(\boldsymbol{u}=(-30,0,0)\) |
| Magnetic field | \(\boldsymbol{B}_0=(1,0,0)\), parallel to the shock normal |
| Equation of state | Ideal gas, \(\gamma=5/3\) |
| MHD numerics | PLM reconstruction, LLF Riemann solver, RK2, CFL 0.3 |
| Modeled shock speed | \(v_{sh}=10\) in the reflecting-wall frame |
| Duration | \(t=0\) to \(3000\,\Omega_0^{-1}\) |
| Initial perturbations | None |

The left reflecting wall converts the leftward inflow into a right-moving
parallel shock. Frame tracking and shock recentering are disabled; all output
coordinates are absolute simulation coordinates.

## MHD-PIC model

- The calculation begins with no CR particles. All particles are created by
  shock injection.
- There is one positively charged species with unit mass and unit charge.
- Particle pushing and field interpolation use the 2D3V Boris-TSC scheme.
  Charge and current use second-order TSC deposition.
- Gas-particle momentum and energy feedback are fully enabled with unit
  coefficients.
- The full CR-Hall electric field, induction, force, and energy terms are
  enabled with `pic_cr_hall_mode = full`. The physical closure follows Bai
  et al. (2015); it is not the Hall-neglecting approximation used by Sun & Bai
  (2023).
- The artificial light speed is \(\mathbb{C}=10^4v_{A,0}\). Particle state is
  stored as momentum per unit mass, and the maximum gyro-angle per particle
  step is 0.3.
- The intended Frontier build is Release, double precision, MPI+HIP for
  `ZEN3/GFX90A`, with `-munsafe-fp-atomics`.

## Shock detection and swept-mass tracer

The detector follows Mignone et al. (2018), Section 4.5.1, Equations 74--75.
A cell is flagged as part of the primary shock when all three conditions hold:

1. the in-plane velocity divergence is negative;
2. the sum of the separately normalized \(x_1\)- and \(x_2\)-pressure
   curvatures exceeds 0.2; and
3. the surrounding \(3\times3\) stencil has
   \(p_{\min}<15p_0\) and \(p_{\max}>250p_0\).

The values 15 and 250 are shock-detection bounds only. They are not injection
carrier thresholds.

One passive scalar stores the tagged mass density \(\rho T\). At the start of
an accumulation cycle it is zero. During every step, currently shock-flagged
cells are set to \(T=1\), and the scalar is then advected conservatively. Define

\[
M_T = \sum_i (\rho T)_i\,\Delta V_i,
\]

and the pressure-qualified swept mass

\[
M_{\rm hot} = \sum_{i:f_i=0,\ p_i>30p_0}
               (\rho T)_i\,\Delta V_i,
\]

where \(f_i\) is the current shock flag. An injection event is considered when
a shock is still detected and

\[
M_{\rm hot} > 0.8M_T.
\]

Thus \(Q=0.8\) remains the global swept-fraction trigger and \(p>30p_0\)
remains a cheap carrier prefilter. Tracer belonging to a carrier that fails the
exact affordability test below is not cleared or partially consumed: it stays
fully pending for a later event.

## Particle injection

For every injection event:

1. **Mass budget.** Before any deferral, the candidate CR mass is
   \(\Delta M_{\rm cr}=\eta M_{\rm hot}\), with \(\eta=2\times10^{-3}\).
   If a carrier is deferred, its mass is removed from the current active
   budget and remains in its tracer. Any active-budget remainder smaller than
   one macro-particle is retained for the next event.
2. **Macro-particle mass.** The target is \(N_{\rho_0}=4\) particles per
   unit-density swept cell. With the unit-depth 2D cell volume,
   \[
   m_{\rm macro}=\frac{\eta\rho_0\Delta V}{N_{\rho_0}}
                =0.0542534722222.
   \]
   All injected particles have unit numerical weight.
3. **Candidate bundles.** Carrier cells are sampled in proportion to their
   local contribution to \(M_{\rm hot}\). A systematic mass-weighted draw
   avoids a second particle-weight type and keeps each cell's realized count
   close to its proportional expectation. The resulting particles are grouped
   into one candidate bundle per carrier, with positions uniform within that
   cell.
4. **Injected momentum.** Directions are isotropic in three dimensions in the
   modeled shock frame. The momentum magnitude is
   \(p_{\rm inj}/m=\sqrt{10}\,u_0=94.8683\), corresponding to approximately
   \(10E_{sh}\), where \(E_{sh}=u_0^2/2\). The state is boosted by
   \(+v_{sh}=+10\) into the reflecting-wall frame.
5. **Carrier-local affordability.** Before changing particles, tracer, or gas,
   AthenaK applies the complete candidate bundle to a copy of the carrier's
   conserved state. The trial removes the bundle's full mass, three momentum
   components, and kinetic energy from that same cell. In addition to producing
   a finite state above the density and pressure floors, the trial must leave
   at least 50% of the cell's pre-subtraction thermal energy.
6. **Stable active set.** A carrier that fails the trial is excluded from the
   current event, all of its tracer remains pending, and none of its proposed
   particles are created. The active swept mass, CR budget, and systematic
   candidate bundles are then recomputed. This repeats until every carrier in
   the active set is affordable. Excluded mass is not redistributed to other
   carriers.
7. **Conservative commit.** Only the stable bundles are appended. Their exact
   mass, momentum, and kinetic energy are then subtracted from their respective
   carrier cells. The accumulated tracer cohort is reset, explicitly deferred
   hot carriers are restored in place, and current shock cells seed the next
   cohort. The one-cell subtraction stencil is mandatory.

Cold swept tracer with \(p\leq30p_0\) contributes neither particles nor an
injection budget and is discarded when an event commits; its mass is not
reassigned to hotter cells. The pressure
criterion does not explicitly require a carrier to lie geometrically
downstream: the prefilter is "tagged, no longer shock flagged, and
\(p>30p_0\)," followed by the exact affordability test.

## Floor behavior and early-particle removal

The gas density and pressure floors are \(10^{-6}\rho_0\) and \(10^{-8}p_0\).
The carrier-local trial is the normal control for an unaffordable subtraction;
it uses the full conserved-state update rather than pressure alone because
removing momentum also changes the remaining gas kinetic energy. There is no
energy clipping, source redistribution, or event-wide/global throttle. A hard
stage-level floor and non-finite check remains after source application as a
fatal backstop for implementation errors, not as part of the injection
algorithm.

Following the start-up treatment in Mignone et al. (2018), the code performs
one removal at \(t\geq960\): shock-injected particles with birth time
\(t_{birth}<480\) are deleted. Their mass, momentum, and energy are recorded as
an explicit particle sink and are not returned to the gas. Injection otherwise
continues for the full calculation.

## Outputs and restart

- History: every \(1\,\Omega_0^{-1}\).
- MHD primitives plus cell-centered magnetic field: every
  \(30\,\Omega_0^{-1}\).
- Deposited CR charge density and three current components: every
  \(30\,\Omega_0^{-1}\).
- Raw particle output, including tags and birth times: every
  \(60\,\Omega_0^{-1}\).
- Restart checkpoint: every \(120\,\Omega_0^{-1}\).

Restart files preserve the tracer accumulation state, sub-particle mass
reservoir, particle tags, injection/removal ledgers, and early-removal flag.
The \(Q=0.8\) trigger, \(p>30p_0\) prefilter, and carrier-affordability rule are
restart-fingerprinted. The affordability implementation changes the restart
schema, so checkpoints made before it are invalid and this calculation must
begin at \(t=0\).

## Relation to the published method

Mignone et al. use all swept tracer
\(\rho_{sh}=\rho T(1-f)\), trigger when
\(\sum\rho_{sh}>0.8\sum\rho T\), and distribute an \(\eta\)-fraction of that
full swept mass locally. AthenaK restricts that mass to the \(p>30p_0\) subset
and defers any carrier whose exact local transaction cannot retain the 50%
thermal reserve. The deferral does not clip a particle or move its source to a
different cell; it leaves that carrier's tagged mass pending. AthenaK's
numerical integrator, particle implementation, and full CR-Hall coupling are
also its own implementations, so this is a controlled reproduction target
rather than a bitwise reproduction of PLUTO.

## References

1. A. Mignone, G. Bodo, B. Vaidya, and G. Mattia (2018),
   "[A Particle Module for the PLUTO Code. I. An Implementation of the MHD-PIC
   Equations](https://arxiv.org/abs/1804.01946)," *ApJ* **859**, 13.
   See Sections 4.5 and 4.5.1, especially Equations 74--77.
2. X.-N. Bai, D. Caprioli, L. Sironi, and A. Spitkovsky (2015),
   "[Magnetohydrodynamic-Particle-in-Cell Method for Coupling Cosmic Rays with
   a Thermal Plasma: Application to Non-relativistic
   Shocks](https://arxiv.org/abs/1412.1087)," *ApJ* **809**, 55.
   See Section 2 for the full MHD-PIC and CR-Hall closure.
3. X. Sun and X.-N. Bai (2023),
   "[The Magnetohydrodynamic-Particle-In-Cell Module in Athena++:
   Implementation and Code Tests](https://arxiv.org/abs/2304.10568)," *MNRAS*
   **523**, 3328--3347. See Section 5.4 for the parallel-shock comparison.
