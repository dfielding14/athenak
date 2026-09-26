# MHD-PIC Questions and Answers

This worksheet records the answers obtained from the implementation on the
`PIC_development` branch at commit `340e74ba75c55b1e0b19f24a6b1a8aec479df199`.
Unless stated otherwise, the answers describe the standard full-$f$ configuration of
the coherent current algorithm selected by
`pic_physical_mode=paper_mhd_pic_vl2_tsc`; the historical `paper_mhd_pic`,
`engineering`, and extension paths are not interchangeable with it. Each answer has
two tiers:

1. A plain-English answer of no more than one paragraph and ideally one or two
   sentences.
2. A detailed answer, when needed, with implementation specifics, exact code
   references, tests, and qualifications.

The answers distinguish intended algorithm design, implemented behavior, and behavior
demonstrated by tests.

## Introduction: runtime model and distribution representation

### What `pic_physical_mode` means

`pic_physical_mode` selects a coherent physical and numerical contract for a run. It
is not a coefficient or a switch for one isolated feature: it fixes the broad meaning
of the particle state and determines which coupling, induction, feedback, and
integration choices are allowed together. Its purpose is to prevent a science run
from being assembled from independent legacy switches that may work individually but
do not form one validated physical model (`src/particles/particles.cpp:378-403`;
`docs/source/engineering/pic_mhd_model_contract.md:5-23`).

The current selector includes several deliberately different lanes. `engineering`
preserves the historical development interface; `paper_test_particle` isolates paper
particle mechanics without backreaction; `paper_mhd_pic` preserves the historical
pre-VL2 chronology; `paper_mhd_pic_vl2_tsc` is the active coupled science model; and
`extended_mhd_pic` contains separately named extensions that require their own
qualification. These names are model identities, not increasing levels of accuracy,
and results from one should not automatically be interpreted as results from another
(`src/particles/particles.hpp:49-51`;
`docs/source/engineering/pic_mhd_model_contract.md:12-23`).

### Why this document focuses on `paper_mhd_pic_vl2_tsc`

`paper_mhd_pic_vl2_tsc` is the current coherent implementation intended for new
coupled paper/science runs. The name summarizes its main ingredients:
the base name denotes a paper-oriented coupled MHD-PIC path, `VL2` denotes the
two-stage van Leer predictor-corrector chronology, and `TSC` denotes the
triangular-shaped-cloud particle interpolation and deposition. At a broad level it
uses relativistic CR momentum, a midpoint Boris particle update, time-centered MHD
coupling, conservative particle--gas momentum exchange, and, when total energy is
evolved, conservative energy exchange. Hall-off and the derived full CR-Hall closure
are two explicitly selected induction policies within this same algorithmic identity
(`docs/source/engineering/pic_mhd_model_contract.md:12-23,101-117`;
`src/particles/particles.cpp:1295-1380`).

### Will there eventually be only one mode?

The checked-in contract does not promise that every `pic_physical_mode` value will be
collapsed into one. What it states explicitly is that new coupled publication runs
should use `paper_mhd_pic_vl2_tsc` or a separately named extension, while
`engineering` and `paper_mhd_pic` remain for migration and historical compatibility.
The current organization is therefore best described as **one main coupled production
path, with distinct test and extension identities where the physical contract really
differs**, rather than one universal mode for every use case. Removing the remaining
identities would be a future design decision, not a commitment in the current source
or model contract (`docs/source/engineering/pic_mhd_model_contract.md:20-23`).

### What a full-$f$ configuration means

Full-$f$ describes how the CR distribution is represented; it is separate from
`pic_physical_mode` and is not the same as the `full` value of
`pic_cr_hall_mode`. In a full-$f$ calculation, the macroparticle ensemble represents
the entire distribution $f(\boldsymbol x,\boldsymbol p,t)$. When moment deposition and
feedback are enabled, each particle contributes its full macroparticle weight, without
a delta-$f$ perturbation multiplier. In a physical delta-$f$ calculation, the
distribution is written as $f=f_0+\delta f$, and particles carry signed weights
representing $\delta f$; deposited particle moments use those weights while configured
analytic background moments are supplied separately. Delta-$f$ can reduce particle
noise when the perturbation is small, but it has distinct loading, feedback, and
validation requirements (`docs/source/modules/particles.md:164-175`;
`src/particles/particles_moments.cpp:1262-1265`).

In this document, “standard full-$f$” means that physical delta-$f$ weighting is not
enabled--normally `pic_deltaf_mode=off`. With `cr_distribution=random`,
`pic_deltaf_mode=quiet_start` changes how particles are initially sampled but still
deposits full moments, so it is not physical delta-$f$. The distinction matters here
because several detailed staging statements differ for physical delta-$f$, and the
current full CR-Hall implementation explicitly requires uniform-grid full-$f$ MHD-PIC
(`src/particles/particles.cpp:655-665,1400-1414`;
`src/particles/particles.hpp:446-451`).

## 1. Time centering of the electromagnetic fields

Are the electromagnetic fields from the MHD state taken at the half MHD time step
($n+1/2$) to integrate the particles? In other words, is the particle integration
second-order accurate in time?

### Plain-English answer

Yes, for `paper_mhd_pic_vl2_tsc`: the real particle kick uses the predicted MHD
state at $n+1/2$ and the particle position at $x^{n+1/2}$. This is the implemented
second-order midpoint/VL2 particle update, although the repository does not yet have
an end-to-end temporal-convergence test of the fully coupled CR--MHD system.

### Detailed answer and evidence

The code does not interpolate an edge-centered MHD electric field. It TSC-interpolates
the midpoint cell-centered fluid velocity and magnetic field (and the Hall drift when
enabled) to the midpoint particle position, constructs

$$
\boldsymbol{cE}_p=-\left(\boldsymbol{u}_p+\boldsymbol{v}_{H,p}\right)
\times\boldsymbol{B}_p,
$$

and applies one full-$\Delta t$ Boris momentum kick.

- The VL2 coefficients are $\beta=1/2$ in stage 1 and $\beta=1$ in stage 2, with
  the corrector rebuilt from the saved $n$ state
  (`src/driver/driver.cpp:155-167`; `src/mhd/mhd_update.cpp:40-42,99-101`).
- The end of stage 1 refreshes the primitive and cell-centered magnetic states
  (`src/mhd/mhd_tasks.cpp:339-374,1754-1772`). The stage-2 pusher reads those
  `w0` and `bcc0` arrays, interpolates them at the already half-drifted particle
  position, and performs the real kick
  (`src/particles/particles_pushers.cpp:470-521,531-554`).
- Position is advanced as two half drifts, one in each VL2 stage
  (`src/particles/particles_pushers.cpp:579-649`). In full-Hall mode only, stage 1
  also makes a scratch half-kick used to predict the midpoint Hall current; it does
  not alter the true particle momentum
  (`src/particles/particles_pushers.cpp:429-432,518-528`).

The focused `pic_relativistic_gyro_paper` regression passed in this review and gave
phase-error convergence orders tending to 2 for three reduced-light-speed choices.
That validates the isolated relativistic Boris update, not the temporal order of the
entire coupled algorithm (`tst/scripts/particles/pic_relativistic_gyro_paper.py:378-386`).

## 2. Predictor-corrector structure of the coupled update

Does the MHD-PIC algorithm employ a second-order, two-stage predictor-corrector
scheme in which midpoint states are first predicted and subsequently used to perform
a time-centered update of the coupled CR-MHD system using $n+1/2$ quantities?

### Plain-English answer

Yes, for `paper_mhd_pic_vl2_tsc`: stage 1 predicts the MHD half-time state and moves
the particles to their half-step positions; stage 2 uses those midpoint quantities
for a time-centered particle kick and rebuilds the final MHD state from the saved
$n$ state.

### Detailed answer and evidence

For the full-$f$ path, the operative sequence is:

1. **Stage-1 predictor:** save the $n$-level MHD state, deposit predictor particle
   moments, drift particles from $x^n$ to $x^{n+1/2}$, and advance MHD by
   $\Delta t/2$. In full-Hall mode, a scratch half-kick predicts the momentum used
   in the midpoint Hall current; Hall-off needs no such scratch kick.
2. **Stage-2 corrector:** use $x^{n+1/2}$ and the predicted MHD
   $\boldsymbol{u}^{n+1/2},\boldsymbol{B}^{n+1/2}$ (plus
   $\boldsymbol{v}_H^{n+1/2}$ in full Hall) for the real full-step Boris kick.
   Deposit the realized momentum and energy changes at the midpoint, complete the
   second half-position drift, and recompute the endpoint MHD state from the saved
   $n$ state with the full $\Delta t$ update.

The full-$f$ particle tasks are deliberately inserted after `MHD::CopyCons` and before
`MHD::Fluxes` in each stage (`src/particles/particles_tasks.cpp:228-401`), and an
abridged MHD causal order is
`CopyCons -> Fluxes -> RKUpdate -> MHDSrcTerms -> CornerE -> CT -> ConToPrim`
(`src/mhd/mhd_tasks.cpp:329-380`). The stage-1 analytic feedback is only a predictor;
the stage-2 update is rebuilt from the saved base state and uses the realized exchange
(`src/mhd/mhd_tasks.cpp:811-857`). Full Hall therefore has a predictor/corrector
nuance: its particle kick uses predicted midpoint $\boldsymbol v_H$, while its final
grid Hall correction uses the deposited realized particle impulse.

Physical delta-$f$ retains the two-stage particle integration but inserts its moment
chain after `RKUpdate` and before `MHDSrcTerms`; it uses analytic
background-plus-perturbation $\rho\boldsymbol E+\boldsymbol J\times\boldsymbol B$
feedback in both stages instead of the full-$f$ realized-impulse corrector
(`src/particles/particles_tasks.cpp:233-239`;
`src/mhd/mhd_tasks.cpp:924-973`).

This answer should not be generalized to the historical or engineering paths, whose
generic particle push is scheduled differently
(`src/particles/particles_tasks.cpp:175-205`;
`src/particles/particles_pushers.cpp:729-737,816-868`).
The dedicated static stage-graph checker also passed during this review
(`tst/scripts/particles/pic_paper_task_stage_trace_vl2_tsc.py:140-255`).

## 3. Final CR energy-momentum feedback and conservation

Is the final energy-momentum feedback from the cosmic rays computed after all
particles have moved through a full time step, with the change in CR momentum and
energy deposited on the grid and the opposite change added as source terms to the gas
momentum and energy to ensure energy-momentum conservation?

### Plain-English answer

Yes for the full-$f$, energy-evolving VL2 path, with one timing correction: the full
particle momentum kick and the deposition of its realized momentum/energy change
occur at $x^{n+1/2}$, before the particle completes its second half-position drift.
The stage-2 gas update receives the exact opposite deposited change.

### Detailed answer and evidence

During stage 2, the code saves each particle's old momentum and kinetic energy,
performs $\boldsymbol p^n\rightarrow\boldsymbol p^{n+1}$ over the full $\Delta t$,
and stores

$$
\frac{\Delta\boldsymbol p_{\rm CR}}{\Delta t},\qquad
\frac{\Delta E_{\rm CR}}{\Delta t}
$$

with the macroparticle normalization included
(`src/particles/particles_pushers.cpp:518-554`). These rates are TSC-deposited while
the particle is still at $x^{n+1/2}$
(`src/particles/particles_moments.cpp:1293-1304`); the second half drift follows
(`src/particles/particles_tasks.cpp:391-393`;
`src/particles/particles_pushers.cpp:647-715`). The gas source then subtracts the
deposited `DPDT` from gas momentum and `DEDT` from gas total energy
(`src/mhd/mhd_tasks.cpp:844-855`). Unit feedback coefficients are enforced in this
mode, so this exchange cannot be silently rescaled
(`src/particles/particles.cpp:1371-1380`).

Thus the conservation statement is exact at the discrete exchange level in a closed
periodic domain, up to floating-point error. With nonperiodic boundaries, the budget
must also include MHD boundary fluxes and particle escape/reflection impulses. In
full-Hall mode the matched Hall Poynting flux is also part of the gas-energy balance,
so cell-by-cell gas energy need not equal only $-\Delta E_{\rm CR}$ even when the
global budget closes.
An isothermal MHD run has no evolved gas-energy equation, so only the momentum part of
this statement applies there.

The focused `pic_paper_coupling_conservation_vl2_tsc` regression passed during this
review for Hall-off and full-Hall configurations; its conservation and coefficient
invariance checks are defined at
`tst/scripts/particles/pic_paper_coupling_conservation_vl2_tsc.py:169-211`. In this
run, the maximum total-momentum residuals were $6.57\times10^{-7}$ and
$9.95\times10^{-7}$, and the total-energy residuals were $2.81\times10^{-6}$ and
$-9.67\times10^{-6}$, respectively.

## 4. Interpolated electric field and its parallel component

Is the interpolated electric field modified to eliminate a spurious component
parallel to the magnetic field introduced by interpolation? Alternatively, is it
calculated as $\boldsymbol{E}=-\boldsymbol{v}\times\boldsymbol{B}$ from separately
interpolated $\boldsymbol{v}$ and $\boldsymbol{B}$ at each particle location, so that
$E_{\parallel}=0$ automatically? If the latter is used, is that treatment an issue
when the CR Hall term is present?

### Plain-English answer

The particle pusher uses the latter construction: it interpolates
$\boldsymbol u$, $\boldsymbol B$, and, in full Hall, $\boldsymbol v_H$, then forms
$\boldsymbol{cE}_p=-(\boldsymbol u_p+\boldsymbol v_{H,p})\times\boldsymbol B_p$;
there is no separate projection. Consequently the particle-sampled field has
$E_\parallel=0$ by construction, including with the CR Hall term, apart from roundoff.

### Detailed answer and evidence

The same TSC stencil separately interpolates cell-centered `bcc0` and primitive
velocity `w0` (`src/particles/field_interpolation.hpp:90-107`). Full Hall also
interpolates `cr_hall_drift`, after which the pusher explicitly evaluates the cross
product (`src/particles/particles_pushers.cpp:493-514`). There is no cleanup step;
the code records $\boldsymbol{cE}_p\!\cdot\!\boldsymbol B_p$ as a diagnostic
(`src/particles/particles_pushers.cpp:540-550`). The Hall-off engineering/passive-MHD
`pic_boris_midpoint_eb` regression passed in this review; in this run the constructed
$\boldsymbol{cE}\cdot\boldsymbol B$ was exactly zero in the measured field and the
output diagnostic was $-4.1\times10^{-16}$. This checks the ordinary particle
cross-product construction, not the full-Hall grid corrector discussed next
(`tst/scripts/particles/pic_boris_midpoint_eb.py:299-302`).

There is an important grid-level qualification. The full-Hall stage-2 corrector does
not reconstruct its Hall field as a cross product; it uses the realized deposited
particle impulse,

$$
\boldsymbol{cE}_{H,\mathrm{corr}}
=-\frac{(\Delta\boldsymbol p_{\rm CR}/\Delta t)_{\rm dep}}
        {\alpha_i\rho},
$$

and does not project this vector against the grid magnetic field
(`src/mhd/mhd_tasks.cpp:88-95,141-159`;
`src/mhd/mhd_corner_e.cpp:339-354`). This is analytically the conservative counterpart
of $-\boldsymbol v_H\times\boldsymbol B$, but interpolation, deposition, and face/edge
reconstruction do not guarantee exact pointwise orthogonality. I found no focused
regression bounding a stage-2 **grid** $E_\parallel$ residual, so whether such a
residual is practically important is not established by the present evidence.

## 5. Feedback source terms and the CR Hall contribution to CT

Are all CR feedback terms implemented as source terms that leave the core AthenaK MHD
integration unchanged? Does the CR Hall term directly enter constrained transport;
that is, is it included when the edge-centered electric field is calculated for the
CT update?

### Plain-English answer

No. Gas momentum--and, for ideal MHD, total-energy--feedback is applied as a source,
and Hall-off leaves the induction update as ordinary ideal-MHD CT, but full CR Hall
directly modifies the magnetic-induction EMFs and the matched total-energy flux used
by the MHD/CT update.

### Detailed answer and evidence

- The conservative gas backreaction is applied in `MHD::MHDSrcTerms`
  (`src/mhd/mhd_tasks.cpp:867-1017`). This does not require changing the MHD Riemann
  solver itself.
- With `pic_cr_hall_mode=off`, the paper mode does not inject deposited CR current
  into CT through the generic `EFieldSrc` path
  (`src/particles/particles.hpp:585-590`). The conservation regression also verifies
  that changing the legacy `couple_j_to_efield_coeff` leaves paper-mode induction
  unchanged (`tst/scripts/particles/pic_paper_coupling_conservation_vl2_tsc.py:193-203`).
- With `pic_cr_hall_mode=full`, `AddCRHallFluxes` is called inside `MHD::Fluxes`,
  before first-order flux correction and the RK update
  (`src/mhd/mhd_tasks.cpp:691-697`). It adds Hall contributions to the face EMFs and
  the matched $(\boldsymbol{cE}_H\times\boldsymbol B)_n$ total-energy flux
  (`src/mhd/mhd_tasks.cpp:507-555`). In 2D/3D, `CornerE` includes the
  stage-appropriate Hall contribution in the edge-centered electric field
  (`src/mhd/mhd_corner_e.cpp:27-43,158-173,328-355`); in 1D the corrected face EMF is
  copied directly to the edge (`src/mhd/mhd_corner_e.cpp:48-65`). CT then advances
  face-centered $\boldsymbol B$ with its curl
  (`src/mhd/mhd_ct.cpp:20-22,38-84`).

Therefore full Hall is part of the core conservative flux/CT sequence, even though it
does not use the separate legacy direct-current `MHD::EFieldSrc` injection route.

## 6. Order of operations with CR subcycling

What is the order of operations when CR subcycling is used? In particular, if the
code uses a two-stage predictor-corrector scheme, how does that scheme operate when
subcycling is enabled?

### Plain-English answer

CR subcycling is not implemented on this branch, so there is no subcycled
predictor-corrector ordering to describe. Particle cell-crossing and gyro-angle limits
instead reduce the one shared MHD/PIC timestep; the two VL2 stages are integrator
stages, not particle subcycles.

### Detailed answer and evidence

`Particles::NewTimeStep` computes a cell-crossing limit and a Boris gyro-angle limit
(`src/particles/particles.cpp:2043-2149`). `Mesh::NewTimeStep` includes the resulting
particle limit in the global timestep and performs the MPI global minimum
(`src/mesh/mesh.cpp:639-648`). The driver then executes one two-stage VL2 loop for
that shared $\Delta t$, with no nested particle-substep loop
(`src/driver/driver.cpp:502-524`). A repository-wide search found no CR-subcycle
parameter or implementation.

The actual non-subcycled order is the stage-1 predictor/first half drift followed by
the stage-2 midpoint kick, realized-impulse deposition, and second half drift described
in question 2. Full Hall adds one scratch half-kick for its midpoint current prediction,
but that scratch operation is not a subcycle and does not update the true particle
momentum.

## 7. AMR compatibility

Is MHD-PIC incompatible with adaptive mesh refinement (AMR)?

### Plain-English answer

No: the blanket statement is false on `PIC_development`. Hall-off
`paper_mhd_pic_vl2_tsc` has explicit static/adaptive mesh-refinement support, but full
CR Hall is explicitly restricted to a uniform grid, and the current dynamic-AMR tests
are bounded regressions rather than broad production qualification.

### Detailed answer and evidence

On a multilevel Hall-off mesh, the paper path allocates a dedicated receiver-record
transport (`paper_smooth`) and deposits each particle's TSC moments at the resolution
of every receiving leaf block (`src/particles/particles.cpp:1479-1502`;
`src/particles/particles_moments.cpp:934-1120,1136-1145`). After adaptive
reconstruction, the mesh refreshes the retained particle module and geometrically
remaps particle ownership (`src/mesh/mesh_refinement.cpp:946-986,1765-1836`). Static
SMR and adaptive-refinement Hall-off test decks exercise these paths; the standalone
receiver-resolution oracle passed all four checks during this review
(`tst/scripts/particles/pic_paper_smooth_tsc_interface.py:330-372`;
`tst/scripts/particles/pic_q009_coupled_boundary_lifetime_vl2_tsc.py:279-365`;
`inputs/tests/pic_paper_smooth_tsc_interface.athinput:34-42,62-100`;
`inputs/tests/pic_q009_coupled_boundary_lifetime_vl2_tsc.athinput:36-42,58-99`).

The hard limitations are:

- `pic_cr_hall_mode=full` rejects every multilevel mesh
  (`src/particles/particles.cpp:1400-1414`).
- Expanding-box mode with active MHD also rejects SMR/AMR
  (`src/particles/particles.cpp:874-888`).
- Particles are remapped after refinement but are not split or merged, so refinement
  does not restore the requested `ppc` in every newly refined cell
  (`src/mesh/mesh_refinement.cpp:1789-1828`).

One checked-in engineering note still says receiver-resolution deposition is missing
(`docs/source/engineering/pic_amr_lifetime_and_interface_policy.md:47-51`); that
statement appears stale relative to the current source and tests above.

## 8. Key CR-module parameters and macroparticle normalization

What are the key parameters for the CR module? After choosing a reduced speed of
light, is a velocity magnitude specified for each species and used to determine its
Lorentz factor? Is a number density or a rest-mass density specified? Is that density
defined per macroparticle, such that the total density in a cell is the number of
particles per cell multiplied by the mass or number density assigned to one
macroparticle? Is $q/(mc)$ specified separately?

### Plain-English answer

For the standard full-$f$ loader, there is no separate generic CR number-density or
rest-mass-density input. You specify the total particles per cell, species mass/charge
and initial state, while `deposit_qscale` together with the particle weight fixes each
macroparticle's normalization; $\gamma$ is derived using the reduced light speed, and
in paper/extension modes `species_charge/species_mass` is already the code-normalized
CR $q/(mc)$.

### Detailed answer and evidence

The main controls are:

| Purpose | Key parameters | Meaning |
|---|---|---|
| Algorithm | `particle_type=cosmic_ray`, `pusher=boris_tsc`, `pic_physical_mode`, `pic_background_mode`, `pic_feedback_mode`, `pic_cr_hall_mode` | Selects the particle model and coupling closure. |
| Relativistic state | `pic_cr_light_speed`, `pic_cr_initial_state=velocity\|momentum` | Sets the reduced $C$ and interpretation of the initial components. |
| Loading/species | `ppc`, `nspecies`, `<species0>/mass`, `charge`, `vx0`, `vy0`, `vz0` (and `species1`, etc.), `cr_distribution` | Sets the aggregate macro-particle count and species properties. |
| Normalization/deposition | `deposit_qscale`, `deposit_order`, `deposit_moments` | Sets macroparticle normalization and grid deposition. The paper VL2 mode requires TSC (`deposit_order=2`). |
| Timestep | `pic_max_cell_cross`, `pic_theta_max` | Limits cell crossing and gyro angle by reducing the shared timestep. |
| Full Hall only | `pic_background_ion_q_over_mc` | Sets the background-ion $q/(mc)$ used by the Hall closure; it is distinct from the CR species ratio. |

These parameters are parsed primarily in `src/particles/particles.cpp:378-835,1690-1773`.

**Initial state and Lorentz factor.** The generic initializer takes Cartesian
components, not one velocity magnitude. Per-species components fall back to the global
`cr_vx0`, `cr_vy0`, and `cr_vz0` values
(`src/particles/particles.cpp:1714-1721`). If
`pic_cr_initial_state=velocity`, the components are physical velocity and are converted
using

$$
\gamma=\frac{1}{\sqrt{1-v^2/C^2}},\qquad
\frac{\boldsymbol p}{m}=\gamma\boldsymbol v.
$$

If `pic_cr_initial_state=momentum`--the default for the paper mode--the components are
already $\boldsymbol p/m$, and runtime kinematics use

$$
\gamma=\sqrt{1+\frac{|\boldsymbol p/m|^2}{C^2}},\qquad
\boldsymbol v=\frac{\boldsymbol p/m}{\gamma}.
$$

See `src/particles/particles.cpp:579-617,1734-1745,1930-1941` and
`src/particles/particles.hpp:137-162`. Specialized problem generators may replace the
generic vector with shells, beams, or other distributions.

**Density and macroparticle weight.** `ppc` is the aggregate number of macroparticles
per cell across all species, realized as a total count per MeshBlock; it is not a
per-species count. Species are assigned round-robin, so, when evenly divisible,
`ppc=128` with two species gives 64 particles of each species per cell **on average
over the MeshBlock**. Random placement does not guarantee exactly 64 of each species
in every cell (`src/particles/particles.cpp:1782-1792,1843-1851`). For a particle
created by the standard initializer in a leaf cell,

$$
w_p=\mathrm{IPWT}=\frac{V_{\rm leaf}}{V_{\rm root}},
$$

and the code uses

$$
M_p=\mathrm{deposit\_qscale}\,w_p m_s,\qquad
Q_p=\mathrm{deposit\_qscale}\,w_p q_s
$$

as its macro rest mass and code-normalized macro charge contribution
(`src/particles/particles.cpp:1959-1961`;
`src/particles/particles_pushers.cpp:538-543`;
`src/particles/particles_moments.cpp:446-493`). On a uniform root grid, for species
$s$ with $N_{s,\mathrm{cell}}$ particles in a cell of root-level volume $V_{\rm root}$,

$$
n_s=\frac{N_{s,\mathrm{cell}}\,\mathrm{deposit\_qscale}}{V_{\rm root}},\qquad
\rho_s=m_s n_s.
$$

Thus “density = `ppc` times mass” is true only after accounting for species allocation,
cell volume, and `deposit_qscale`. These formulas describe the standard `ppc`
initializer; specialized/manual injection paths may assign `IPWT` differently. The
deposited `prtcl_rho` is the signed code moment $Q_{\rm CR}=q_{\rm CR}/c$, not
rest-mass density or ordinary physical charge density
(`MHD_PIC_CR_HALL_CODE_MAP.md:77-95,149-169`).

**Charge-to-mass ratio.** The initializer stores

$$
\mathrm{IPM}=\frac{\mathrm{species\_charge}}{\mathrm{species\_mass}}
$$

(`src/particles/particles.cpp:1942-1945`). In the paper/extension normalization this
ratio already represents $q/(mc)$, so it is not divided by the reduced $C$ again
(`MHD_PIC_CR_HALL_CODE_MAP.md:300-316`). The separately named
`pic_background_ion_q_over_mc` is required only for the background ions in full-Hall
mode (`src/particles/particles.cpp:499-500,1400-1407`).

## Appendix A. Is the full-Hall grid-level parallel electric field a problem?

### Current assessment

Potentially, but it has not been demonstrated to be a bug. The concern is confined to
the full-Hall stage-2 **grid** corrector: if the part of its EMF parallel to a
consistently collocated grid magnetic field remains finite under convergence and makes
a non-negligible contribution to $\nabla\times\boldsymbol E$, it could produce an
unphysical magnetic-field update. The particle pusher itself is not affected by this
specific concern because its electric field is constructed explicitly as a cross
product and is perpendicular to its interpolated $\boldsymbol B$.

In the continuum Hall closure,

$$
\boldsymbol{cE}_H=-\boldsymbol v_H\times\boldsymbol B,
$$

so $\boldsymbol{E}_H\cdot\boldsymbol B=0$. The stage-2 grid corrector instead derives
the Hall EMF from the deposited, realized particle impulse
(`src/mhd/mhd_tasks.cpp:88-95,141-159`) and subsequently reconstructs it to faces and
edges (`src/mhd/mhd_corner_e.cpp:328-355`). Deposition, averaging, reconstruction, and
division by the grid density do not commute with taking
$\boldsymbol v_H\times\boldsymbol B$. Consequently, the resulting grid vector is not
guaranteed to be perpendicular to a separately reconstructed grid magnetic field.

This observation alone does **not** prove that the method is wrong. A deposited or
cell-averaged $\langle\boldsymbol v_H\times\boldsymbol B\rangle$ need not be
perpendicular to $\langle\boldsymbol B\rangle$, even when the underlying pointwise
field is perpendicular everywhere. Moreover, constrained transport stores
$\boldsymbol E$ and $\boldsymbol B$ at different staggered locations, so a discrete
statement about $E_\parallel$ is meaningful only after specifying how the magnetic
field is collocated with each edge EMF. A naive projection against a convenient
cell-centered or edge-reconstructed $\boldsymbol B$ is therefore not automatically
more accurate and could make the induction corrector inconsistent with the deposited
impulse or with the matched face/edge reconstruction.

The practical issue is the curl, not merely a nonzero dot product. A parallel residual
that converges away as grid resolution and particle sampling increase is ordinary
truncation or sampling error. A component that does not converge away and contributes
materially to the CT curl could alter magnetic evolution or topology and would indicate
that the corrector needs revision. The existing total energy-momentum conservation
regression does not directly answer this question because it does not measure the
parallel EMF or isolate its contribution to the CT update.

### Focused test needed

Use a smooth, deterministic full-Hall problem and construct a magnetic field at the
same edge locations as the stage-2 Hall EMF using a documented collocation rule. Then
measure, with safe handling of vanishing fields,

$$
\epsilon_\parallel =
\frac{|\boldsymbol E_H\cdot\boldsymbol B_*|}
     {|\boldsymbol E_H|\,|\boldsymbol B_*|+\epsilon},
$$

and the dynamically relevant ratio

$$
R_{\rm CT}=
\frac{\|\nabla\times\boldsymbol E_{H,\parallel}\|}
     {\|\nabla\times\boldsymbol E_H\|+\epsilon},
\qquad
\boldsymbol E_{H,\parallel}=
\frac{\boldsymbol E_H\cdot\boldsymbol B_*}{|\boldsymbol B_*|^2+\epsilon}
\boldsymbol B_*.
$$

The study should refine grid resolution and particles per cell separately, reporting
at least RMS and high-percentile values of $\epsilon_\parallel$ and $R_{\rm CT}$. If
both decrease at the expected discretization/noise rate, no algorithmic change is
indicated. If either approaches a nonzero value or changes the converged magnetic
solution, the next step is to compare the existing impulse-derived corrector with a
consistently projected or cross-product-based corrector; any replacement must preserve
the face/edge CT construction and the conservative energy-momentum budget.
