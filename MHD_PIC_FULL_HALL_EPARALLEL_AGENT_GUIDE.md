# Agent Guide: Diagnose the Full-Hall Grid-Level Parallel EMF

## Mission

Determine whether the stage-2, impulse-derived CR-Hall grid electromotive force (EMF)
contains a parallel component that materially corrupts AthenaK's constrained-transport
(CT) magnetic update. Diagnose the current method before proposing a fix. Do **not**
add a projection to the production update unless the investigation demonstrates a
non-convergent, dynamically significant error.

Work on the current `PIC_development` branch and record the exact starting commit. The
implementation mapped when this guide was written was commit
`340e74ba75c55b1e0b19f24a6b1a8aec479df199`. Read the repository `AGENTS.md` files
that govern every directory you touch.

Restrict the primary claim to the coherent, uniform-grid, full-$f$ path:

```text
pic_physical_mode = paper_mhd_pic_vl2_tsc
pic_cr_hall_mode  = full
```

Do not mix this investigation with the legacy `engineering` current-to-CT path,
physical delta-$f$, expanding-box physics, AMR, resistivity, or other nonideal EMFs.
Full Hall is already restricted to uniform-grid, full-$f$, ideal MHD by
`src/particles/particles.cpp:1400-1422`.

## Required outcome

Finish with one of these evidence-based classifications:

1. **No actionable concern:** the residual is roundoff in uniform fields and converges
   away under controlled **joint** refinement of grid spacing, particle sampling, and
   timestep; its CT contribution also converges to zero. One-factor sweeps may stop at
   floors set by the other two errors.
2. **Collocation-dependent residual with no demonstrated CT effect:** a cell/face dot
   product remains nonzero for a stated collocation, but the operational CT difference
   converges away under joint refinement and the magnetic solution is unchanged at
   convergence.
3. **Methodological concern:** the residual's CT contribution or magnetic-solution
   effect approaches a nonzero value after particle, grid, and timestep errors are
   controlled.
4. **Inconclusive:** particle noise, arbitrary collocation, limiter activation, or an
   incomplete diagnostic prevents a defensible conclusion.

Do not call the current method correct or incorrect before reaching one of these
outcomes. If the result is inconclusive, state exactly what remains unresolved rather
than guessing.

## Physics question

AthenaK stores $\boldsymbol{cE}=c\boldsymbol E$. In the continuum large-scale CR-Hall
closure,

$$
Q_e=\alpha_i\rho+Q_{\rm CR},
\qquad
\boldsymbol v_H=
\frac{\boldsymbol K_{\rm CR}-Q_{\rm CR}\boldsymbol u_g}{Q_e},
$$

and

$$
\boldsymbol{cE}_H=-\boldsymbol v_H\times\boldsymbol B,
\qquad
\boldsymbol{cE}_H\cdot\boldsymbol B=0.
$$

Here $Q_{\rm CR}=\rho_{q,{\rm CR}}/c$, where $\rho_{q,{\rm CR}}$ is the signed CR
charge density, and $\boldsymbol K_{\rm CR}=\boldsymbol J_{\rm CR}/c$. The input
$\alpha_i=q_i/(m_i c)$ uses the background-ion charge and mass, so
$\alpha_i\rho$ is the background-ion charge density divided by $c$. Equivalently, if
$\rho_{q,i}$ denotes that charge density, then
$\alpha_i=\rho_{q,i}/(\rho c)$. The CR Lorentz-force density satisfies

$$
\boldsymbol F_{\rm CR}
=Q_{\rm CR}\boldsymbol{cE}
 +\boldsymbol K_{\rm CR}\times\boldsymbol B
=\alpha_i\rho\,\boldsymbol v_H\times\boldsymbol B,
$$

so the equivalent continuum corrector is

$$
\boldsymbol{cE}_H=-\frac{\boldsymbol F_{\rm CR}}{\alpha_i\rho}.
$$

The stage-2 code uses the **realized deposited particle impulse** to form the
exact-discrete Hall corrector. Particle--gas exchange conservation has a separate,
more direct cause: the code deposits `DPDT/DEDT` and applies their opposites to the gas
with the same RK weight in `MHD::AddPaperVL2FeedbackSource`
(`src/mhd/mhd_tasks.cpp:844-856`). In a periodic domain this balances global particle
and gas momentum and energy to numerical tolerance. Using `DPDT` in the Hall EMF ties
the field closure to that realized impulse, but is not by itself what conserves the
exchange.

The concern is introduced by coarse graining, not by the continuum identity. Each
particle impulse may be perpendicular to the magnetic field sampled by that particle,
while

$$
\sum_p W_p\,\Delta\boldsymbol p_p
$$

need not be perpendicular to a separately averaged cell or reconstructed face
magnetic field when $\boldsymbol B$ varies across the TSC stencil. In general,

$$
\left\langle\boldsymbol v_H\times\boldsymbol B\right\rangle
\cdot\left\langle\boldsymbol B\right\rangle\ne0
$$

does not imply that the underlying pointwise field had a parallel component.
Deposition, division by grid density, component-wise PLM reconstruction, and GS07
face-to-edge assembly also need not preserve perpendicularity.

The physical danger is therefore not merely a nonzero grid dot product. It is a
non-convergent contribution to Faraday's law,

$$
\delta(\partial_t\boldsymbol B)
=-\nabla_h\times\boldsymbol{cE}_{H,\parallel}.
$$

A component with zero **discrete CT curl** makes no contribution to the discrete CT
magnetic increment. A persistent curlful residual could behave like an unmodeled
nonideal term, altering the magnetic evolution or topology. Convergence of that
discrete effect under controlled refinement is the relevant continuum evidence; a raw
dot product alone establishes neither continuum correctness nor separately collocated
Hall-energy-flux consistency. Global energy--momentum conservation and CT's
preservation of $\nabla\cdot\boldsymbol B$ do not by themselves exclude a local
induction error.

There is also no unique raw "edge $E\cdot B$" in the present staggering. `efld.x1e`,
`x2e`, and `x3e` are different vector components stored on differently oriented
edges. Any edge-level dot product must first define and document how all three EMF
components and all three magnetic components are brought to one common location.

The implementation equations and normalization are documented in
`docs/source/engineering/pic_mhd_model_contract.md:59-104,119-158` and
`MHD_PIC_CR_HALL_CODE_MAP.md:44-138,149-175,206-296`.

## Code path to trace

Follow the quantity from the particle midpoint field to CT. Do not diagnose only one
end of this chain.

| Stage | Relevant code | What to verify |
|---|---|---|
| Particle midpoint field | `src/particles/particles_pushers.cpp:420-554` | TSC interpolation of $\boldsymbol u_g$, $\boldsymbol v_H$, and $\boldsymbol B$ at `493-506`; cross-product field at `508-514`; full stage-2 Boris kick at `518-521`; realized `IPDPX:IPDPZ`, `IPDE`, sampled $B$, and `IPEBDOT` at `531-550`. |
| Stage selection | `src/particles/particles_moments.cpp:28-39` | In full Hall, the generic realized-impulse moments run only in stage 2. |
| TSC impulse deposition | `src/particles/particles_moments.cpp:446-505,1253-1304` | `IPDPX:IPDPZ` are divided by cell volume and accumulated as `IMOM_DPXDT:IMOM_DPZDT` with the TSC shape. |
| Deposition ordering | `src/particles/particles_tasks.cpp:228-401` | Deposition and synchronization finish before `MHD::Fluxes`; the final half drift follows deposition. |
| Cell corrector | `src/mhd/mhd_tasks.cpp:88-95` | Stage-2 $\boldsymbol{cE}_{H,c}=-\boldsymbol{\mathrm{DPDT}}/(\alpha_i\rho)$; there is no projection. |
| Face/cell closure construction | `src/mhd/mhd_tasks.cpp:45-195` | The predictor reconstructs $\boldsymbol v_H$ and $\boldsymbol B$ separately and then crosses them at `45-87`; the corrector reconstructs the components of $\boldsymbol{cE}_I$ at `97-160`; cell states are built at `162-195`. The same face $B$ is used for the matched $(\boldsymbol{cE}_H\times\boldsymbol B)_n$ energy flux. |
| Face flux insertion | `src/mhd/mhd_tasks.cpp:471-558,655-704` | Tangential Hall EMFs and the matched total-energy flux are added before FOFC and RK update. |
| FOFC | `src/mhd/mhd_tasks.cpp:563-648`; `src/mhd/mhd_fofc.cpp:489-495` | Donor-cell Hall induction and energy terms are restored together on fallback faces. Exclude FOFC initially, then test it separately. |
| Cell-to-edge construction | `src/mhd/mhd_corner_e.cpp:27-43,158-220,328-452` | Stage 2 uses the impulse-derived cell corrector in the 2D/3D GS07 edge construction. In 1D, corrected face EMFs are copied directly at `48-66`. |
| CT update | `src/mhd/mhd_ct.cpp:20-86` | The actual magnetic update uses the native staggered curl of `efld`. |
| Opposite gas impulse | `src/mhd/mhd_tasks.cpp:783-857`, especially `844-856` | The gas receives the opposite stage-2 `DPDT/DEDT`; this is the direct action--reaction mechanism. Do not change it during diagnosis. |
| Existing Hall diagnostics | `src/particles/particles_moments.cpp:781-823`; `src/outputs/history.cpp:241-262,336-363` | Current history exposes only `hall_Rmax` and `hall_Lmax`, not a grid $E_\parallel$ measure. |
| Existing particle diagnostic | `src/outputs/derived_variables.cpp:2718-2726` | `prtcl_ebdot` is the pusher's sampled $\boldsymbol{cE}_p\cdot\boldsymbol B_p$; it is not the grid-corrector diagnostic. |

No existing output exposes the Hall-only `MHD::efld` consumed by CT
(`src/mhd/mhd.hpp:133-141`). Existing `prtcl_dpxdt`, `prtcl_dpydt`, and
`prtcl_dpzdt` outputs (`src/outputs/derived_variables.cpp:2670-2703`) can provide a
cheap offline cell-centered screen when combined with density and `mhd_w_bcc`, but
those endpoint outputs are not a substitute for a stage-local midpoint or CT-edge
diagnostic.

## Investigation rules

- Begin with read-only tracing and reproduce the current baseline.
- Add diagnostics without changing the production EMF, face energy flux, gas source,
  task order, or particle state.
- Make any diagnostic path default-off or test-only and avoid a production performance
  cost when it is disabled.
- Preserve the existing meanings of `hall_Rmax`, `hall_Lmax`, and `prtcl_ebdot`.
- Measure stage 1 and stage 2 separately. Stage 1 constructs the Hall field as an
  explicit cross product and is an important internal control; stage 2 is the target.
- Use one periodic MeshBlock, `fofc=false`, no resistivity, and no extra electric-field
  sources for the first controlled tests. Add decomposition, FOFC, and other complexity
  only after the core result is understood.
- Do not infer a Hall-only edge field by subtracting independent Hall-off and full-Hall
  simulations: those runs evolve different particle and gas states. Build a shadow
  diagnostic from the same stage state.
- Keep double- and single-precision tolerances separate. Scale roundoff checks to
  machine precision and the local field magnitudes rather than using an unexplained
  absolute constant.

## Diagnostics to implement

### 1. Particle-level impulse control

The existing `IPEBDOT` verifies
$\boldsymbol{cE}_p\cdot\boldsymbol B_p$. The corrector, however, uses the realized
impulse, so also measure

$$
\epsilon_{p,\parallel,2}=
\left[
\frac{\displaystyle\sum_p
      \frac{(\Delta\boldsymbol p_p\cdot\boldsymbol B_p)^2}
           {|\boldsymbol B_p|^2+B_{\rm floor}^2}}
     {\displaystyle\sum_p |\Delta\boldsymbol p_p|^2+\epsilon}
\right]^{1/2}.
$$

All required particle quantities are present after
`src/particles/particles_pushers.cpp:538-550`. For a frozen midpoint field this
should be at roundoff. A failure here identifies the pusher or diagnostic itself,
before deposition is involved. The stored `IPDPX:IPDPZ` are rates rather than raw
increments, but the common timestep factor cancels in this ratio.

### 2. Cell-level corrector

At the beginning of stage-2 `MHD::AddCRHallFluxes`, use the exact midpoint arrays that
feed the live update:

$$
\boldsymbol{cE}_I
=-\frac{\boldsymbol{\mathrm{DPDT}}}{\alpha_i\rho},
\qquad
\boldsymbol{cE}_X=-\boldsymbol v_H\times\boldsymbol B_{cc}.
$$

Do not compute this later from endpoint output $B$, which is at the wrong time level.
Record:

$$
\epsilon_{\parallel,2}=
\left[
\frac{\displaystyle\sum_c V_c
      \frac{(\boldsymbol{cE}_I\cdot\boldsymbol B_{cc})^2}
           {|\boldsymbol B_{cc}|^2+B_{\rm floor}^2}}
     {\displaystyle\sum_c V_c|\boldsymbol{cE}_I|^2+\epsilon}
\right]^{1/2},
$$

$$
R_{\rm closure}=
\frac{\|\boldsymbol{cE}_I-\boldsymbol{cE}_X\|_2}
     {\|\boldsymbol{cE}_X\|_2+\epsilon}.
$$

$R_{\rm closure}$ is the complete cell-closure mismatch, not a pure parallel metric;
use $\epsilon_{\parallel,2}$ to isolate the collocated parallel projection.

Also report an absolute RMS of $\boldsymbol{cE}_I\cdot\boldsymbol B_{cc}$ and a
masked 95th or 99th percentile of the local normalized dot product. Define the mask
from $|\boldsymbol{cE}_I||\boldsymbol B|$ relative to its global maximum and report
the threshold; otherwise near-zero fields will dominate normalized maxima. For both
particle and grid projection norms, set
$B_{\rm floor}^2=\epsilon_{\rm mach}B_{\max}^2$, record $B_{\max}$ and the resulting
floor, and configure the primary tests without magnetic nulls so this guard is
numerically inactive.

A low-intrusion instrumentation point is immediately after the aliases are captured
in `src/mhd/mhd_tasks.cpp:471-500`.

### 3. Face-level corrector

Measure the same quantities using the **exact selected face states** returned by
`CRHallGridFaceState`. Expose the reconstructed $B_x,B_y,B_z$ already computed at
`src/mhd/mhd_tasks.cpp:147-152`, then accumulate orientation-specific diagnostics
immediately after the calls at `513-515`, `532-534`, and `550-552`.

Do not naively reduce over those live loop bounds. They contain tangential ghost and
extended faces, with still wider ranges when FOFC is enabled, so doing so can duplicate
faces and make the result decomposition-dependent. Either write diagnostic scratch
arrays and reduce over uniquely owned physical faces afterward, or define and test an
explicit canonical face-ownership rule.

Report stage-1 predictor and stage-2 corrector values separately. Stage 1 performs the
cross product after face reconstruction and should be perpendicular at this location;
it is the best check that the face metric and its $B$ collocation are correct.

### 4. CT-level dynamical effect

Cell and face dot products are screens, not the final decision. Define two operational
maps that return Hall-only magnetic rates without RK timestep factors:

- $\mathcal C_I[\boldsymbol e]$ follows the actual **impulse-corrector** route:
  component-wise PLM reconstruction of the cell-centered closure $\boldsymbol e$.
- $\mathcal C_X[\boldsymbol v_H,\boldsymbol B]$ follows the actual
  **cross-product** route: reconstruct $\boldsymbol v_H$ and transverse
  $\boldsymbol B$ separately, use the staggered face-normal $B$, and only then take
  the face cross product, while using $-\boldsymbol v_H\times\boldsymbol B_{cc}$ for
  the cell term.

Both maps must use the frozen live upwind and FOFC choices, matched face Hall-energy
terms, the appropriate cell and face contributions consumed by GS07 `CornerE`,
edge-boundary exchange, and finally the CT curl stencil at
`src/mhd/mhd_ct.cpp:52-83`.

At a documented cell collocation, define

$$
P_\perp\boldsymbol{cE}_I=
\boldsymbol{cE}_I-
\frac{\boldsymbol{cE}_I\cdot\boldsymbol B_{cc}}
     {|\boldsymbol B_{cc}|^2+B_{\rm floor}^2}\boldsymbol B_{cc},
$$

and evaluate

$$
\boldsymbol C_I=\mathcal C_I[\boldsymbol{cE}_I],\qquad
\boldsymbol C_\perp=\mathcal C_I[P_\perp\boldsymbol{cE}_I],\qquad
\boldsymbol C_{X|I}=\mathcal C_I[\boldsymbol{cE}_X],\qquad
\boldsymbol C_X=\mathcal C_X[\boldsymbol v_H,\boldsymbol B].
$$

$\boldsymbol C_{X|I}$ is a useful same-operator diagnostic: comparing it with
$\boldsymbol C_I$ isolates the difference between the two **cell closures** after the
same component-wise reconstruction. It does not reproduce AthenaK's cross-product
algorithm. $\boldsymbol C_X$ does reproduce that algorithm; comparing it with
$\boldsymbol C_I$ includes both the cell-closure difference and the fact that
reconstruction and the cross product do not commute.

The operational CT contribution associated with removing the collocated parallel
component is

$$
\boldsymbol C_\parallel^{\rm op}=\boldsymbol C_I-\boldsymbol C_\perp,
\qquad
R_{\parallel,{\rm CT}}=
\frac{\|\boldsymbol C_\parallel^{\rm op}\|_2}
     {\|\boldsymbol C_I\|_2+\epsilon}.
$$

This definition avoids pretending that AthenaK stores one collocated edge vector. It
also accounts for the fact that 2D/3D `CornerE` consumes both cell and face Hall terms;
projecting only face terms or one final edge component is not a consistent alternate
operator. Because limiting is nonlinear, label this an **operational projection
sensitivity**, not an exact decomposition of the live staggered edge field.

If the face-level screen dominates while the cell-projected sensitivity is negligible,
add a second, explicitly labeled comparator that projects both the selected face terms
against their face $B$ and the cell terms against $B_{cc}$. Route those two projected
sets together through GS07. Never project only the faces while leaving `CornerE`'s cell
Hall term unmodified.

Also report both comparisons explicitly,

$$
R_{{\rm closure}|I,{\rm CT}}=
\frac{\|\boldsymbol C_I-\boldsymbol C_{X|I}\|_2}
     {\|\boldsymbol C_{X|I}\|_2+\epsilon},
\qquad
R_{{\rm full-route},{\rm CT}}=
\frac{\|\boldsymbol C_I-\boldsymbol C_X\|_2}
     {\|\boldsymbol C_X\|_2+\epsilon}.
$$

The practical per-stage measure is

$$
r_{\Delta B,\parallel}=
\frac{|\beta_{\rm stage}|\,\Delta t\,
      \|\boldsymbol C_\parallel^{\rm op}\|_2}
     {\|\boldsymbol B\|_2+\epsilon}.
$$

Use the actual CT stage weight `beta_dt = beta[stage-1]*dt`, not a hard-coded
$\Delta t$; it differs in the stage-1 control. The dimensionless
$r_{\Delta B,\parallel}$ measures practical per-stage impact only. It is **not** a
timestep-convergence proof because it vanishes trivially as $\Delta t\to0$ even if the
instantaneous rate error does not. Base correctness on
$\boldsymbol C_\parallel^{\rm op}$, $R_{\parallel,{\rm CT}}$, or a rate normalized by
a fixed physical Hall timescale.

Retain every numerator and denominator separately so a nearly zero total Hall curl
cannot make a ratio misleading. Even the same-operator closure ratio can include
perpendicular or finite-timestep differences, and the full-route ratio additionally
includes reconstruction-order effects; do not mislabel either as pure $E_\parallel$.
Define every CT norm as an RMS over uniquely owned physical staggered magnetic faces,
with one documented dual-volume weighting, so MeshBlock interfaces are not double
counted and decomposition changes do not change the metric.

Construct separate diagnostic-only shadow edge arrays and give them synchronization
equivalent to the live `SendE/RecvE` sequence. Evaluate their curls only after the
shadow receive completes and before the live `CT` task; a curl taken immediately after
`CornerE` is invalid at periodic and MeshBlock boundaries. Do not alter the production
field.

As secondary evidence in a periodic smooth problem, a consistently collocated volume
integral of $\boldsymbol{cE}_H\cdot\boldsymbol B$ may be reported as a **continuum
helicity-source proxy**. It is not a discrete magnetic-helicity diagnostic without a
staggered $\boldsymbol A\cdot\boldsymbol B$ definition and the appropriate boundary
terms, and it cannot replace the native CT comparison.

## Test campaign

Run the tests in this order. Do not begin with a shock or nonlinear Bell run.

### Test 0: Baseline and instrumentation self-check

Before adding diagnostics, run:

```bash
cd tst
python3 run_tests.py particles/pic_paper_coupling_conservation_vl2_tsc
python3 run_tests.py particles/pic_boris_midpoint_eb
```

The first is the current Hall-off/full-Hall global exchange regression
(`tst/scripts/particles/pic_paper_coupling_conservation_vl2_tsc.py:1-212` and
`inputs/tests/pic_paper_coupling_conservation_vl2_tsc.athinput`). The second validates
the existing particle cross-product diagnostic but is not a full-Hall grid test
(`tst/scripts/particles/pic_boris_midpoint_eb.py:1-349`; its final invariant check is
at `298-303`). Record the commands and measured residuals. Rerun them after every
diagnostic or method experiment.

Diagnostic self-checks must show that:

- Hall-off does not populate Hall-only diagnostics;
- the full-Hall particle $\boldsymbol{cE}_p\cdot\boldsymbol B_p$ and
  $\Delta\boldsymbol p_p\cdot\boldsymbol B_p$ controls are at roundoff;
- the stage-1 face cross-product metric is at roundoff;
- all reported denominators and masks are finite and nontrivial.

### Test 1: Uniform oblique-field control

Create a small 3D periodic problem with one MeshBlock, uniform gas, a uniform magnetic
field with all three components nonzero, and a uniform CR drift not parallel to
$\boldsymbol B$. Use one particle at each cell center for the exactly uniform control,
or implement an explicit test-only subcell lattice, and choose parameters so the Hall
EMF is nonzero. Repeating particles at the same center is not increased spatial
sampling. Run only one or a few steps.

Because every deposited and reconstructed state is uniform, particle-, cell-, and
face-level parallel residuals should be roundoff, while the Hall CT curl should vanish.
A residual or Hall magnetic update above scaled roundoff is an immediate implementation
or diagnostic error. Since the exact Hall curl is zero here, judge absolute curls and
increments rather than ratios with vanishing denominators. Do not proceed to
interpreting smooth-field convergence until this control passes.

### Test 2: Smooth, varying, divergence-free field

Use a periodic 3D field with variation in at least two directions, nonzero components,
and no magnetic nulls. One suitable analytically divergence-free form is

$$
\boldsymbol B=
\left(B_{x0}-a\sin ky,\;
      B_{y0}+a\sin kx,\;
      B_{z0}\right),
$$

with $|a|$ small compared with the background magnitude. Initialize face-centered
$B$ analytically so the discrete divergence is at roundoff. Use smooth gas fields and
a nonzero Hall drift. For deterministic spatial refinement, implement an explicitly
test-only quiet subcell lattice; for sampling refinement, use random placement with
multiple recorded seeds. The existing generic full-$f$ choices are center and random,
and duplicated cell-center particles do not provide high spatial sampling. A one-step
or few-step test is preferred so nonlinear evolution and limiter activation do not
obscure the operator being tested.

If the existing `linear_wave` setup cannot create this controlled state, add one
focused test problem generator; do not distort a production problem generator.

First perform one-factor-at-a-time studies. Each can converge only to the floor set by
the axes held fixed:

| Axis | Minimum useful sweep | What it isolates |
|---|---|---|
| Grid spacing | At least three resolutions, with a quiet lattice or particle noise demonstrated subdominant | TSC/PLM/GS07 truncation and collocation error. |
| Particles per cell | At least three values at fixed grid; rescale `deposit_qscale` to hold physical CR density/current fixed | Particle sampling error. Use several seeds if loading is random. |
| Timestep | At least three **recorded actual timesteps** at fixed grid and particle load | Finite-step difference between predicted cross-product and realized-impulse correctors. Make the case CFL-limited or document whether cell-crossing, gyro-angle, or Hall-wave bounds control each run. |
| Decomposition | One MeshBlock baseline, then a split MeshBlock and serial/MPI comparison | Moment/face/edge synchronization artifacts. |

Then perform at least three **joint-refinement** levels. Reduce $h$ and the actual
$\Delta t$ together, and raise particle sampling enough that the one-factor studies
show sampling noise below the spatial/temporal signal. A claim that the total residual
converges away requires this joint sequence; three unrelated one-factor plateaus are
not sufficient.

When changing resolution or `ppc`, preserve the intended physical CR density and
current. For the standard initializer, use

$$
\mathrm{deposit\_qscale}
=\mathrm{deposit\_qscale}_{\rm ref}
 \frac{V_{\rm cell}}{V_{{\rm cell},{\rm ref}}}
 \frac{\mathrm{ppc}_{\rm ref}}{\mathrm{ppc}}.
$$

For sampling studies use random placement with an explicit `pic_random_seed`
(`src/particles/particles.cpp:777-782,1911-1915`) and repeat seeds when the first trend
is noisy. Do not use the Q043 Bell deck for a `ppc` noise sweep: it requires
`cr_distribution=center`, and multiple particles can duplicate the same cell-center
location (`src/pgen/tests/q043_bell_current_volume_aware.cpp:262-300`;
`src/particles/particles.cpp:1853-1877`).

For a smooth deterministic case, second-order spatial behavior is the natural target,
although PLM limiting can reduce the observed order. Stochastic sampling noise should
typically decrease approximately as `ppc^-1/2`. Fit and report the observed rates;
do not force these expectations onto data that are not in the asymptotic regime, and
expect a one-factor sweep to plateau once another error source dominates.

### Test 3: Representative full-Hall physics

Only after Tests 1 and 2 are understood, reuse the compact full-Hall Bell candidate at
`inputs/tests/pic_cr_hall_bell_linear.athinput`. Add the new diagnostics without
changing its normalization. Along the **live trajectory**, compare instantaneous and
accumulated one-step CT increments from the impulse, projected, same-operator
cross-product, and actual cross-product-route shadows. This is a same-state sensitivity
diagnostic, not an independently evolved solution; it cannot supply alternate growth
rates, phases, or nonlinear spectra.
Include at least one resolution repeat and one decomposition check. Do not interpret a
larger centered `ppc` as a sampling-noise study unless the loader is changed to place
those particles at distinct, controlled subcell locations.

Reuse the existing standalone analyzer at
`tst/publication/analyze_pic_cr_hall_bell.py:1-255`; it already compares Hall-off,
base full-Hall, and higher-resolution full-Hall results and measures growth, frequency,
polarization, deposited charge/current, $R$, and $\Lambda$. Its higher-resolution
normalization must continue to satisfy the volume-aware relation enforced at
`src/pgen/tests/q043_bell_current_volume_aware.cpp:359-410`.

This Bell case is a representative endpoint, not the primary oracle: its effective
dimensionality and growing mode can hide collocation or sampling errors that the smooth
3D manufactured problem exposes directly.

### Test 4: Evolved counterfactual, only if warranted

The same-state projected and cross-product shadows in Diagnostics 4 are mandatory.
Only evolve alternate algorithms if $R_{\parallel,{\rm CT}}$,
$R_{{\rm closure}|I,{\rm CT}}$, $R_{{\rm full-route},{\rm CT}}$, or an accumulated
rate difference approaches a nonzero plateau after controlled joint refinement.

Run separate, opt-in evolutions of:

1. the current impulse-derived corrector;
2. the actual cross-product-route corrector;
3. optionally, a consistently projected complete corrector.

Route each through the same face reconstruction, Hall energy flux, FOFC handling,
CornerE construction, communication, and CT stencil. An edge-only projection is not a
valid production candidate. Once evolved, these cases no longer share the same state;
only this separate experiment can compare alternate Bell growth rate, phase, or
spectrum. If a production change is eventually proposed, it must retain the opposite
`DPDT/DEDT` particle--gas exchange and use a Hall energy flux matched to the alternate
induction closure.

## Decision criteria

Use convergence, not a single arbitrary threshold.

- **Immediate concern:** the uniform control has a residual larger than scaled
  roundoff, or an absolute Hall CT update above scaled roundoff, despite a nonzero
  well-conditioned Hall EMF.
- **Likely truncation/sampling:** each one-factor study is consistent with a documented
  floor and decreases when its targeted component is resolved, while controlled joint
  refinement sends all of $\|\boldsymbol C_I-\boldsymbol C_\perp\|$,
  $\|\boldsymbol C_I-\boldsymbol C_{X|I}\|$, and
  $\|\boldsymbol C_I-\boldsymbol C_X\|$ toward zero.
- **Collocation-dependent residual with no demonstrated CT effect:** a chosen cell/face
  dot product plateaus or varies with collocation while the native complete-operator CT
  difference converges to zero under joint refinement.
- **Methodological concern:** an instantaneous CT-rate metric, or a separately evolved
  magnetic growth/phase/spectrum difference, approaches a nonzero plateau after grid,
  particle, and timestep errors are jointly controlled.
- **Not evidence:** passing global conservation or maintaining
  $\nabla\cdot\boldsymbol B=0$ alone. Both can hold while the local curl is inaccurate.

A full-route plateau alone does not diagnose $E_\parallel$: it can include the
noncommutativity of reconstruction and the cross product. If the full-route difference
persists while $\boldsymbol C_I-\boldsymbol C_\perp$ converges away, report a broader
closure/reconstruction discrepancy rather than a parallel-EMF defect.

Do not use $r_{\Delta B,\parallel}$ versus timestep as a correctness-convergence metric:
its explicit RK timestep factor forces it downward. Use it only to state the practical
size of one stage at a fixed run configuration. When the reference Hall curl is nearly
zero, report the absolute rate and, if useful, normalize it by one fixed physical Hall
timescale rather than by a vanishing quantity.

Compare any apparent parallel contribution with the ordinary magnetic truncation
error, for example the difference between successive grid resolutions. A residual far
below that error is not presently a limiting defect even if it is algebraically
nonzero.

## Why an immediate projection is not justified

- The CT staggering does not provide one unique three-component edge $\boldsymbol B$
  for projection.
- Projection against an averaged field can remove a legitimate covariance contained
  in $\langle\boldsymbol v_H\times\boldsymbol B\rangle$.
- The current corrector is tied exactly to the deposited particle impulse; projection
  changes that impulse-derived Hall closure even though the opposite gas source can
  remain action--reaction conservative.
- Projecting only the final edge EMF would bypass the shared face upwinding and FOFC
  construction and would not define a consistent face/cell/edge closure. A component
  parallel to a given face $B$ drops out of that face's
  $\boldsymbol{cE}\times\boldsymbol B$ energy flux, but staggered edge projection is
  not automatically projection against that same $B$.
- A smaller $E\cdot B$ under one arbitrary collocation is not sufficient evidence that
  the magnetic update is more accurate.

If a persistent CT error is demonstrated, design the alternate closure as a complete
face/energy/edge operator and re-run every conservation and staging regression.

## Deliverables

Produce all of the following:

1. A concise implementation note listing the exact starting SHA and every source file
   changed for diagnostics.
2. A focused input deck and automated analyzer for the uniform and smooth tests. Use a
   name such as `pic_full_hall_eparallel_vl2_tsc`; follow the existing
   `inputs/tests/` and `tst/scripts/particles/` conventions.
3. Tables of absolute and normalized particle, cell, face, and CT metrics for every
   refinement axis, including fitted convergence rates and mask definitions.
4. Plots only where they add information: convergence versus $h$, `ppc`, and
   $\Delta t`, plus one magnetic-solution comparison if a plateau is found.
5. Baseline and post-change results for
   `pic_paper_coupling_conservation_vl2_tsc` and `pic_boris_midpoint_eb`, plus the new
   focused regression.
6. A final classification using the four outcomes at the top of this guide, with the
   evidence that supports it and any remaining uncertainty.

Do not merge a projection or other algorithm change as part of the diagnostic phase.
Leave any such comparison clearly labeled as a shadow/counterfactual experiment for
review.
