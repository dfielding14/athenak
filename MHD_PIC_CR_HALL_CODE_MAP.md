# CR-Hall closure: paper-to-code map

This is the implementation contract for the uniform-grid, full-f
`paper_mhd_pic_vl2_tsc` CR-Hall path. It is deliberately narrow: its purpose is to
make the normalization, signs, staggering, and time centering reviewable before the
physics is trusted. It is not a general plasma-closure document.

## References and model scope

- [Bai et al. (2015)](https://arxiv.org/abs/1412.1087) is the primary equation
  source, especially equations 10--20 and Appendix A.
- [Mignone et al. (2018)](https://arxiv.org/abs/1804.01946) is the independent
  conservative implementation cross-check. Section 3.2.1 is the important
  midpoint CR-Hall predictor reference.
- [Sun & Bai (2023)](https://arxiv.org/abs/2304.10568) is the reference for the
  Athena++-style VL2/Boris/TSC integration pattern. That implementation uses
  ideal-MHD induction and does not supply the missing CR-Hall closure.

The supported first implementation is non-relativistic ideal MHD plus relativistic
full-orbit CRs on a uniform Cartesian mesh. Thermal electrons are massless, thermal-ion
mass supplies the gas density, CR mass is dynamically negligible, and resolved scales
are much larger than the thermal-ion inertial length. Conventional Hall MHD, electron
pressure/inertia, resistivity, delta-f, expanding boxes, and AMR are outside this first
scope. The formal ordering requires `|R| << 1`; the code reports `R` and `Lambda` but
does not invent an abort threshold.

## Signed physical definitions

Use signed charge and current throughout:

```math
q_i + q_e + q_{\rm cr}=0, \qquad
\boldsymbol J_{\rm cr}=\sum_s q_s n_s\boldsymbol u_s,
```

where `q_i > 0` for the background ions and `q_cr` may have either sign. Define

```math
R \equiv \frac{q_{\rm cr}}{|q_e|}
  =\frac{q_{\rm cr}}{q_i+q_{\rm cr}}, \qquad
\boldsymbol E_0=-\frac{\boldsymbol u_g\times\boldsymbol B}{c}.
```

After dropping the conventional thermal-plasma Hall and electron-pressure terms, the
large-scale CR-Hall closure is

```math
\boldsymbol E=\boldsymbol E_0
-\frac{\boldsymbol J_{\rm cr}-q_{\rm cr}\boldsymbol u_g}
       {|q_e|c}\times\boldsymbol B.
```

The direct current-difference form is the code form. Do not construct
`u_cr = J_cr/q_cr`: that quotient is noisy or undefined in empty and
charge-cancelling cells.

The force and power transferred *to the particles* are

```math
\boldsymbol F_{\rm cr}
=q_{\rm cr}\boldsymbol E
 +\frac{\boldsymbol J_{\rm cr}}{c}\times\boldsymbol B,
\qquad
P_{\rm cr}=\boldsymbol J_{\rm cr}\cdot\boldsymbol E.
```

The gas receives `-F_cr` and `-P_cr`. The gas total-energy flux also contains

```math
\boldsymbol F_{E,H}
=\frac{c}{4\pi}(\boldsymbol E-\boldsymbol E_0)\times\boldsymbol B.
```

The source and flux correction are both required. Applying one without the other is not
the Bai closure.

## AthenaK normalization

AthenaK stores the electric field used by CT and the pusher as `cE`, not `E`. Magnetic
fields use the usual Athena normalization in which magnetic energy is `B^2/2`, so the
code has absorbed `4*pi`.

Define the actual stored CR moments

```math
Q_{\rm cr}\equiv q_{\rm cr}/c, \qquad
\boldsymbol K_{\rm cr}\equiv\boldsymbol J_{\rm cr}/c,
```

and the background-ion parameter and electron denominator

```math
\alpha_i\equiv\frac{q_i}{\rho c}, \qquad
Q_e\equiv\frac{|q_e|}{c}=\alpha_i\rho+Q_{\rm cr}.
```

`<particles>/pic_background_ion_q_over_mc` supplies the positive physical `alpha_i`.
It is not a tunable Hall-strength coefficient. A non-positive or non-finite `Q_e` makes
this closure undefined and is a real input/state error.

The code equations are therefore

```math
\boldsymbol D = \boldsymbol K_{\rm cr}-Q_{\rm cr}\boldsymbol u_g,
\qquad
\boldsymbol v_H=\frac{\boldsymbol D}{Q_e},
```

```math
\boldsymbol{cE}_0=-\boldsymbol u_g\times\boldsymbol B,
\qquad
\boldsymbol{cE}_H=-\boldsymbol v_H\times\boldsymbol B,
\qquad
\boldsymbol{cE}=\boldsymbol{cE}_0+\boldsymbol{cE}_H.
```

The equivalent force identities provide a compact sign check:

```math
R=Q_{\rm cr}/Q_e,
```

```math
\boldsymbol F_{\rm cr}
=Q_{\rm cr}\boldsymbol{cE}+\boldsymbol K_{\rm cr}\times\boldsymbol B
=(1-R)\left(
 Q_{\rm cr}\boldsymbol{cE}_0+
 \boldsymbol K_{\rm cr}\times\boldsymbol B\right).
```

Particle power is

```math
P_{\rm cr}=\boldsymbol K_{\rm cr}\cdot\boldsymbol{cE}
=(1-R)\boldsymbol K_{\rm cr}\cdot\boldsymbol{cE}_0.
```

Because `4*pi` is absorbed, the Hall energy flux added in AthenaK units is simply

```math
\boldsymbol F_{E,H}=\boldsymbol{cE}_H\times\boldsymbol B.
```

There is no extra factor of the physical light speed or `pic_cr_light_speed` in any of
these grid equations.

## Stored quantities and locations

| Physical quantity | AthenaK representation | Location and role |
| --- | --- | --- |
| Particle position | `prtcl_rdata[IPX:IPZ]` | Particle; true position |
| Particle `p/m` | `prtcl_rdata[IPVX:IPVZ]` | Particle; relativistic state despite legacy names |
| Particle `q/(mc)` | `prtcl_rdata[IPM]` | Particle; `species_charge/species_mass` |
| Macro charge | `deposit_qscale*IPWT*species_charge` | Particle contribution before division by cell volume |
| `Q_cr` | `moments[IMOM_RHO]` or `cr_hall_moments[IMOM_RHO]` | Cell-centered signed TSC deposit |
| `K_cr` | `IMOM_JX:IMOM_JZ` | Cell-centered signed TSC deposit of `Q_cr*v` |
| Particle impulse rate | `IMOM_DPXDT:IMOM_DPZDT`, `IMOM_DEDT` | Cell-centered TSC deposit at `x^(n+1/2)` |
| Gas state | `MHD::w0`, `MHD::u0` | Cell-centered primitive/conserved state |
| Magnetic field | `MHD::b0`, `MHD::bcc0` | Face-centered CT field and cell-centered interpolation field |
| Hall drift | `Particles::cr_hall_drift` | Cell-centered `v_H=D/Q_e`, retained for a stage |
| MHD face EMFs | `e2x1/e3x1`, `e1x2/e3x2`, `e1x3/e2x3` | Riemann face EMFs plus reconstructed Hall correction before `CornerE` |
| MHD edge EMF | `MHD::efld` | GS07 edge-centered total `cE`; `CT` applies `-curl(cE)` |
| Gas Hall energy flux | direct face-flux divergence in `MHDSrcTerms` | Uses the identical reconstructed Hall face state as CT; not stored in `uflx(IEN)` |

`DepositCellCenteredMomentsShape<2>` in
`src/particles/particles_moments.cpp` divides each macro contribution by the local cell
volume and uses the same TSC shape for charge, current, and impulse. Moment boundary
exchange accumulates overlapping contributions before the closure is built.

The relevant code path starts in:

- mode parsing, storage, and restart fingerprinting:
  `src/particles/particles.hpp` and `src/particles/particles.cpp`;
- TSC moments and `v_H` construction: `src/particles/particles_moments.cpp`;
- staged ordering: `src/particles/particles_tasks.cpp`;
- midpoint Boris push: `src/particles/particles_pushers.cpp`;
- gas exchange and matched face induction/energy fluxes: `src/mhd/mhd_tasks.cpp`;
- GS07 construction of the total edge EMF: `src/mhd/mhd_corner_e.cpp`;
- constrained transport sign convention: `src/mhd/mhd_ct.cpp`.

## Exact VL2 time-centering contract

`paper_mhd_pic_vl2_tsc` selects a two-stage VL2 predictor-corrector in
`src/driver/driver.cpp`: stage 1 has `beta=1/2`, and stage 2 rebuilds the final state
from the saved base state with `beta=1`. `w0` and `bcc0` remain at the input time level
until `ConToPrim` closes each stage.

The existing Hall-off sequence is:

| Stage | Particle/MHD state and operation |
| --- | --- |
| 1 | MHD advances provisionally from `n` to `n+1/2`; the particle kick is a no-op; `Q_cr^n,K_cr^n` are deposited at `x^n`; the true particle position drifts to `x^(n+1/2)`; analytic base-time feedback advances the gas predictor. |
| 2 | MHD fluxes use the stage-1 half-time state; the full Boris kick advances `p^n` to `p^(n+1)` at fixed `x^(n+1/2)` using half-time fields; actual `Delta p/Delta t` and `Delta E/Delta t` are deposited there; the second half drift advances `x` to `n+1`; the gas receives the opposite deposited impulse once. |

Full Hall adds one necessary scratch predictor. The stage-2 pusher currently precedes
the ordinary stage-2 moment deposit, so reading ordinary `moments` in that pusher would
use stale base-time current. The correct sequence is:

1. At stage 1, deposit and synchronize signed `Q_cr^n,K_cr^n` at `x^n`. Build
   `v_H^n` from the base gas state.
2. Use the full base field to make a first-order half-step prediction
   `p* = p^n + (dt/2) dp/dt|_n + O(dt^2)` in scratch storage. The Mignone
   semi-implicit predictor or an algebraically equivalent half Boris step is sufficient.
   Do not change the true particle momentum.
3. Use `v_H^n` consistently for the stage-1 analytic gas feedback and matched Hall
   face induction/energy fluxes. Drift the true position to `x^(n+1/2)`.
4. At stage 2, deposit and synchronize `Q_cr^(n+1/2),K_cr^(n+1/2)` from
   `x^(n+1/2),p*`. Rebuild and retain `v_H^(n+1/2)` from the half-time gas state.
5. Advance the true momentum from `p^n` to `p^(n+1)` with the full Boris kick at
   `x^(n+1/2)`, using the same retained Hall drift as CT and the energy-flux update.
6. Deposit the actual particle momentum and kinetic-energy changes at
   `x^(n+1/2)`. Apply their opposites to the gas exactly once. Do not also apply the
   analytic force in the stage-2 corrector.
7. Reconstruct `v_H` and transverse `B` to faces with AthenaK's limited PLM,
   select one interface state by the ordinary density-flux sign, and combine it with
   the unique staggered face-normal `B`. Add its `cE_H` components to the face EMFs
   and apply the matching Hall energy flux with the stage-2 weight. `CornerE` then
   builds the edge EMF, CT updates `B`, and the second particle half drift completes.

The scratch predictor exists only to time-center the Hall field. It is not a second
particle evolution and contributes no gas impulse. If the existing `IPEX:IPEZ` slots are
reused for predicted `p/m`, their temporary meaning must end before the real push rewrites
them with their normal sampled-electric-field diagnostics.

For the particle interpolation, construct the full field from one interpolated `u_g`,
`v_H`, and `B`:

```math
\boldsymbol{cE}_p
=-[\boldsymbol u_{g,p}+\boldsymbol v_{H,p}]\times\boldsymbol B_p.
```

This keeps `cE_p dot B_p=0` without a separate cleanup. The grid discretization uses
the same retained stage drift but not the particle interpolation stencil. At each face,
limited PLM reconstructs all components of `v_H` and the two transverse components of
`bcc0`; the normal magnetic component is the exact staggered face value. The density
flux chooses the lower or upper interface state. Its `cE_H=-v_H x B` augments the face
EMF, while `CornerE` uses the total cell field `-(u_g+v_H) x bcc0` in the GS07 formula.
In one dimension the face component is copied directly. This is the Appendix-B ordering:
interface reconstruction precedes upwind selection.

## Momentum and energy bookkeeping

The sign contract is simple and should stay visible in the kernels:

| Update | Stage-1 predictor | Stage-2 corrector |
| --- | --- | --- |
| Particle momentum | scratch prediction only | `+Delta p` from full-field Boris push |
| Gas momentum | `-beta*dt*F_cr^n` | `-Delta p` from deposited actual impulse |
| Particle kinetic energy | scratch prediction only | `+Delta E` from full-field Boris push |
| Gas total energy source | `-beta*dt*K_cr dot cE` | `-Delta E` from deposited actual impulse |
| Gas Hall energy flux | `-beta*dt*div(cE_H x B)` | same at the retained half-time state |
| Magnetic field | `-beta*dt*curl(cE_0+cE_H)` | same at the retained half-time state |

The Hall energy term is a conservative flux divergence implemented as a separate update
in `MHDSrcTerms`, not by mutating `uflx(IEN)`. The same helper, interface orientation,
limited state, and staggered normal field used for the induction correction supply
`(cE_H x B)_n`. On a uniform periodic mesh, each face contribution enters its two
neighboring cells with opposite signs. It is not another particle-work source.

In a periodic domain, the volume integral of the Hall flux divergence vanishes. Gas plus
particle momentum and energy then change only by the normal MHD truncation/roundoff
budget. This is the decisive check against a sign error or double counting.

## Artificial-`C` semantics

`<particles>/pic_cr_light_speed = C` controls only relativistic particle kinematics:

```math
\gamma=\sqrt{1+|\boldsymbol p/m|^2/C^2}, \qquad
\boldsymbol v=(\boldsymbol p/m)/\gamma.
```

`CRLorentzFactor`, `CRVelocityFromState`, and `CRKineticEnergy` in
`src/particles/particles.hpp` implement these relations. `species_charge/species_mass`
and `IPM` already represent `q/(mc)` in the simulated system. Never divide them, the
deposited current, `alpha_i`, or the Lorentz force by artificial `C` again.

For the present Bai-style runs, `C` is the simulated relativistic transition speed and is
chosen comfortably above all MHD and intended CR drift speeds. It is not yet the
Ji--Hopkins physical reduced-speed-of-light formulation.

## Runtime and diagnostic contract

- `pic_cr_hall_mode=full` enables the pusher, CT, momentum/energy exchange, and Hall
  energy flux as one closure.
- `pic_cr_hall_mode=off` removes the whole Hall correction while retaining the coupled
  Sun--Bai momentum/energy path.
- `current_to_ct_experimental` and `couple_j_to_efield_coeff` are historical engineering
  controls, not alternate strengths of the physical closure.
- Historical decks remain explicitly `off`; new coupled science decks explicitly select
  `full` and provide `pic_background_ion_q_over_mc`.

Only two new global history quantities are required:

```math
\max |R|=\max\left|\frac{Q_{\rm cr}}{Q_e}\right|,
\qquad
\max\Lambda=\max\frac{|\boldsymbol D|}{Q_e v_A}
=\max\frac{|\boldsymbol v_H|}{v_A},
\quad v_A=|B|/\sqrt{\rho}.
```

These report model regime and Hall importance; they are not warning classes or automatic
termination policies. Bell dispersion analysis may additionally use a signed parallel
`Lambda`, but the runtime history maximum is the non-negative magnitude above.

The retained drift and these diagnostics are stage-centered quantities. After a cycle,
the history and Hall contribution to the next CFL estimate use the most recent stage-2
midpoint deposit; the initialization record is zero before the first deposit. This avoids
an otherwise redundant endpoint particle deposit and MPI exchange. In the intended
`|R| << 1` regime, the existing particle crossing bound is tighter than Hall advection
for the Bell, shock, and turbulent-box inputs used here.

## Compact Bell sign and dispersion check

For the linear validation with `B0=rho=1`, parallel positive CR current `K`,
uniform charge `Q`, `a=1-R`, and parallel Hall drift `V_H=K/Q_e`, define
`b_+=B_y+iB_z` and use `b_+ proportional to exp(s*t+i*k*x)`. The full closure
reduces to

```math
(s+i aQ)(s+i kV_H)=a k(aK-k).
```

The Bell-unstable branch therefore has positive `k` and

```math
\omega_r=\frac{aQ+kV_H}{2}, \qquad
\gamma^2=a k(aK-k)-\frac{(aQ-kV_H)^2}{4}.
```

This convention matters: `B_y=A cos(kx), B_z=+A sin(kx)` is the unstable
winding for `K>0`, and the corresponding analyzer coefficient is
`mean[(B_y+iB_z) exp(-ikx)]`. The historical opposite winding is a stable
oscillatory branch and must not be compared with the unstable dispersion.

The immutable-source limited-face qualification in worker job `4973137` passed all 12
checks. At `R=0.01`, `Lambda=1`, it measured `gamma=5.5185` and
`omega_r=2.5860`, versus `5.5576` and `2.5755` from the equation above.
Doubling the parallel resolution gave `gamma=5.5484` and `omega_r=2.5827`,
changes of `0.54%` and `0.13%` from the base case. The Hall-off control measured
`gamma=6.2102` and `omega_r=0.06279`, versus `6.2829` and `0.06283`.
Raw moment outputs independently gave `J_cr/c=(4*pi,0,0)` and
`Q_cr/c=0.04*pi`; the dominant wavelength and unstable polarization also passed.
These results qualify the signed algebra and the final limited-PLM spatial route.

## Short review checklist

Before accepting the implementation, verify directly that:

1. changing a CR species charge sign changes both `Q_cr` and `K_cr` signs;
2. `Q_e=alpha_i*rho+Q_cr`, not `alpha_i*rho-Q_cr`;
3. the pusher and CT both use `cE_H=-v_H x B` from the same stage state;
4. stage 2 uses predicted midpoint current, not stale stage-1 or post-kick current;
5. the actual particle impulse is returned to the gas once and with the opposite sign;
6. the Hall Poynting term is `+cE_H x B` in the flux and therefore
   `-dt*div(cE_H x B)` in the update; and
7. no artificial-`C`, free Hall coefficient, or extra `4*pi` enters the closure.
