# Q-019 Scientifically Qualifying Nonlinear Bell Saturation Campaign Design

Date: 2026-06-06

Status: **design frozen; execution and scientific qualification remain blocked**

This record designs the nonlinear Bell campaign that follows the existing
Q-023 nonlinear fixed-current-like source-local foundation, itself derived
from the Q-023 Sun--Bai Section 5.2 linear carrier. It does not authorize a
run, bind a clean executable, inspect qualifying output, set a
literature-comparison tolerance, or close
`CLAIM-PROD-BELL-NONLINEAR-NOHALL-001`.

## Decision

The recommended first qualifying campaign is a **high-rigidity,
approximately-fixed-current, no-Hall nonlinear continuation of Q-023** in
`paper_mhd_pic_vl2_tsc`.

The primary scientific claim is deliberately narrow:

> AthenaK resolves nonlinear growth and a non-box-limited saturated state, if
> one is reached, for the exact tested high-rigidity, no-Hall,
> ideal-MHD-carrier Bell setup, with the saturation mechanism classified from
> current, drift, wavelength, and energy-transfer diagnostics.

The campaign must not claim CR-deflection or anisotropic-pressure saturation
unless a separate finite-rigidity campaign closes its own source, mapping,
convergence, and acceptance gates. It must not claim an exact locked-current
experiment because AthenaK currently advances the Q-023 CR particles rather
than locking their velocities. It must not claim a direct reproduction of the
driven-boundary survey of Zacharegkas et al. because the current source-local
foundation does not continuously inject fresh CRs through CR-specific
boundaries.

## Binding to the Existing Q-023 Nonlinear Foundation

The direct engineering predecessor is
`inputs/tests/pic_q023_prod_bell_nonlinear_foundation_vl2_tsc.athinput`, bound
by
`tst/publication/readiness/q023_prod_bell_nonlinear_source_local_foundation_2026-06-06.json`.
Q-019 preserves its runtime, no-Hall/full-feedback physics, Q-023
normalization, `P0=1`, `1e-4` right-polarized carrier, current product
`PPC*deposit_qscale=2e9`, native cell size, and normalized-time convention
`tau=k0 U_A t`.

The qualifying baseline deliberately enlarges the foundation's
eight-wavelength periodic 2D domain to 16 wavelengths, adds a deterministic
seeded divergence-free symmetry-breaking perturbation, requests particle
diagnostics, and replaces its source-local single-run heuristics with the
fail-closed ensemble and convergence gates below. The foundation's
`tau=[14,22]` diagnostic window is only a pilot prior; it is not inherited as
a science acceptance window.

## Primary-Literature Interpretation

Only primary literature is used for the scientific distinctions in this
design.

### Riquelme & Spitkovsky: fixed current and plasma acceleration

Riquelme & Spitkovsky (2009) explicitly locked the CR velocities in their
constant-current experiments, thereby removing CR backreaction. They found
that plasma acceleration along the CR drift reduces the current in the plasma
frame and gives an intrinsic saturation condition near the point where the
Alfven speed in the amplified field becomes comparable to the CR drift speed.
In multiple dimensions, turbulence transfers magnetic power to longer
wavelengths; they also identify earlier saturation by CR deflection when the
CR gyroradius approaches the dominant wavelength.

AthenaK can test the **plasma-acceleration branch approximately**, using the
extreme Q-023 rigidity and a diagnostic requirement that the lab-frame CR
current and CR momentum remain effectively unchanged. AthenaK cannot currently
reproduce the exact locked-particle experiment without a separately named and
qualified locked-current source.

Primary source:
[Riquelme & Spitkovsky 2009, ApJ 694, 626](https://doi.org/10.1088/0004-637X/694/1/626);
[arXiv:0810.4565](https://arxiv.org/abs/0810.4565).

### Gargate et al.: self-consistent hybrid evolution

Gargate et al. (2010) used hybrid simulations with kinetic CRs and kinetic
background ions. Their nonlinear evolution includes background-plasma
acceleration and energy gain, CR-current evolution, and eventual CR
isotropization. They attribute saturation to the loss of scale separation and
background-plasma demagnetization while the driving current is maintained,
not to a universal anisotropic-pressure law.

AthenaK can measure the corresponding fluid acceleration, gas energization,
CR current, and CR isotropization diagnostics. It cannot test kinetic
background-ion demagnetization, ion-cyclotron-scale physics, or a
Bell-to-Weibel transition because the carrier is ideal MHD.

Primary source:
[Gargate et al. 2010, ApJL 711, L127](https://doi.org/10.1088/2041-8205/711/2/L127);
[arXiv:1002.1701](https://arxiv.org/abs/1002.1701).

### Zacharegkas et al.: self-consistent anisotropic momentum-flux scaling

Zacharegkas et al. (2024; arXiv first posted in 2022) used driven hybrid
simulations in which CRs enter and leave the domain while the thermal plasma
and fields are periodic. They report that the saturated transverse magnetic
pressure scales with the initial net CR momentum flux or anisotropic pressure,
rather than with CR energy flux alone. Their simulations also show coupling of
CRs and gas near saturation and migration of magnetic power toward scales
comparable to the CR gyroradius.

AthenaK can compute the required CR momentum tensor from particles and test an
**undriven periodic-box analogue** after a finite-rigidity distribution
generator is implemented. The current Q-023 source cannot directly reproduce
the driven-boundary experiment, and the exact AthenaK-to-reference
normalization map remains open.

Primary source:
[Zacharegkas et al. 2024, ApJ 967, 71](https://doi.org/10.3847/1538-4357/ad3960);
[arXiv:2210.08072](https://arxiv.org/abs/2210.08072).

### Q-023 and Sun--Bai foundation

The existing Q-023 source implements the Sun & Bai Section 5.2 normalization,
right-polarized eigenmode, and analytical linear dispersion relation. It uses
kinetic CR particles coupled conservatively to an ideal-MHD carrier. The
nonlinear campaign must inherit that source-local normalization and must close
the Q-023 linear predecessor before using the nonlinear output as evidence.

Primary source:
[Sun & Bai 2023, arXiv:2304.10568](https://arxiv.org/abs/2304.10568).

## What AthenaK Can Actually Support

### Supported by `paper_mhd_pic_vl2_tsc`

- Relativistic CR momentum state and Boris--TSC particle push.
- Ideal-MHD induction, with no direct CR-current Hall term.
- Self-consistent CR momentum change and conservative opposite momentum and
  kinetic-energy feedback to ideal MHD.
- Periodic 1D-carrier, 2D3V, and 3D geometries.
- Deterministic particle loading and deterministic seeded source generation.
- Built-in MHD histories, raw `mhd_w_bcc` fields, deposited particle-current
  fields, particle PVTK positions/velocities/weights, and restarts.
- Reconstruction of CR momentum from output velocity and the recorded
  artificial light speed for the nonrelativistic Q-023 baseline.

### Missing or requiring new Q-019 source-local work

- No exact locked-CR-current mode.
- No current Q-023 nonlinear ensemble generator with seeded
  symmetry-breaking perturbations, enlarged matrix boxes, and exact Q-019
  deck validation.
- No current hot drifting isotropic-shell generator for a
  Zacharegkas-like anisotropic-pressure survey.
- No CR-specific driven inflow/outflow boundary that leaves MHD periodic.
- No built-in CR pressure-tensor or CR kinetic-energy history.
- No nonlinear Bell production analyzer, independent recompute, resource
  model, or reviewed saturation window.
- No kinetic background ions/electrons, ion inertial scale, WICE physics,
  background-ion demagnetization, or electron-scale competing modes.

### Consequence for claims

The first production result can be a bounded high-rigidity no-Hall result. A
pressure-law or CR-deflection result requires a separate finite-rigidity
campaign. A driven Zacharegkas reproduction and exact Riquelme--Spitkovsky
locked-current experiment are unsupported until new, separately named source
mechanics are implemented and qualified.

## Governing Dimensionless Parameters

All parameters below must be emitted by the future deck validator and analysis
record. Exact reference-specific unit maps remain a Q-022 gate.

| Parameter | Definition in this campaign | Scientific role |
|---|---|---|
| `epsilon` | `U_A / v_CR` | Q-023 streaming-speed parameter and linear growth control |
| `M_A,CR` | `v_CR / U_A = 1/epsilon` | CR drift Mach number |
| `gamma0` | `k0 U_A sqrt(1-epsilon^2)` | Q-023 fastest-mode linear growth rate |
| `gamma0 / Omega0` | growth rate divided by initial CR gyrofrequency | High-rigidity and time-scale-separation control |
| `lambda0` | `2 pi / k0` | Q-023 fastest-mode wavelength |
| `k0 r_g0` | `k0 (p_CR/m) / Omega0` for the initial beam | Separates high-rigidity fixed-current-like behavior from CR deflection |
| `C / v_CR` | artificial light speed divided by CR speed | Reduced-light-speed validity |
| `beta_code` | `P_gas / (B0^2/2)` in AthenaK units | Exact thermal-pressure scope of the MHD carrier |
| `J_CR / (2 B0 C k0)` | Q-023 current-normalization closure | Confirms the inherited linear source |
| `rho_CR / rho0`, `E_CR / P_B0` | macro-particle mass and energy loading | Backreaction and fixed-current-like inertia scope |
| `Pi_CR,ij / P_B0` | initial CR momentum-flux tensor normalized by initial magnetic pressure | Required for any Zacharegkas-style pressure claim |
| `r_g(B) / lambda_d` | instantaneous CR gyroradius divided by dominant magnetic wavelength | Deflection/resonance discriminant |
| `L / lambda0` | box size in fastest-mode wavelengths | Nonlinear wavelength-migration and box-limit control |
| `Delta x / lambda0`, PPC, `Delta t gamma0` | numerical controls | Resolution, sampling, and time-step qualification |

For the recommended Q-023 continuation,
`k0 r_g0` is approximately `2.5e6`. This is a derived AthenaK design value,
not a literature tolerance. It makes the baseline a strong high-rigidity
test and makes a CR-deflection or pressure-saturation interpretation
inappropriate unless the diagnostics demonstrate otherwise.

## Recommended Baseline

The recommended production baseline is:

| Field | Recommended value |
|---|---|
| campaign branch | `Q019-HR-FIXED-CURRENT-LIKE-NOHALL` |
| runtime | `paper_mhd_pic_vl2_tsc` |
| dimension | 2D3V |
| induction | ideal MHD, `pic_cr_hall_mode=off` |
| feedback | coupled conservative momentum and energy feedback |
| normalization | inherit Q-023 exactly: `rho=B0=U_A=lambda0=1`, `k0=2 pi` |
| streaming | `epsilon=0.4`, hence `v_CR/U_A=2.5` |
| CR gyrofrequency | inherit Q-023: `Omega0=1e-6 k0 U_A` |
| artificial light speed | inherit Q-023: `C=1000 v_CR=2500` |
| gas | ideal, `gamma=5/3`, `P=1`; no broader beta claim |
| initial Q-023 eigenmode amplitude | `1e-4`, inherited from the nonlinear foundation |
| source perturbation | new deterministic seeded, divergence-free symmetry-breaking perturbation; exact spectrum and relative energy freeze after excluded pilots |
| box | 16 times the Q-023 2D periodic domain |
| mesh | `1024 x 512 x 1`, Q-023 nonlinear-foundation native cell size, `64 x 64 x 1` MeshBlocks |
| particles | 16 PPC, with `deposit_qscale=1.25e8`, inherited from the nonlinear foundation |
| seeds | Q-023 qualifying seeds `23050101` through `23050108` |
| pilot seeds | Q-023 excluded pilot seeds `23050091`, `23050092` |
| provisional pilot horizon | `tau=k0 U_A t=30`; exact production horizon and saturation window freeze after excluded pilots |

The enlarged box, pilot horizon, and qualifying sensitivities are
campaign-design choices, not values or tolerances quoted from the literature.
The 16-wavelength baseline is paired with a 32-wavelength box test because
nonlinear Bell power migrates to longer wavelengths and a plateau at the
fundamental box mode is not scientific saturation. Every PPC row must set
`deposit_qscale=2e9/PPC`; changing PPC without preserving this product changes
the physical current and is rejected.

## Qualifying Matrix

The exact deck checksums, node counts, walltimes, output cadence, terminal
time, saturation window, and numeric acceptance tolerances must be frozen in a
versioned preregistration successor before qualifying output is inspected.

### Core 2D3V matrix

Use all eight qualifying seeds at every row and pair the same seed across
rows.

| Row | Box scale | Resolution scale | PPC | Time-step scale | Role |
|---|---:|---:|---:|---:|---|
| `F2D` | 16 | 1 | 16 | 1 | recommended fiducial |
| `B2D` | 32 | 1 | 16 | 1 | box-size gate |
| `R2D-L` | 16 | 0.5 | 16 | 1 | low-resolution trend |
| `R2D-H` | 16 | 2 | 16 | 1 | high-resolution gate |
| `P2D-L` | 16 | 1 | 8 | 1 | low-PPC trend |
| `P2D-H` | 16 | 1 | 32 | 1 | high-PPC gate |
| `T2D-L` | 16 | 1 | 16 | 0.5 | smaller-step gate |
| `T2D-H` | 16 | 1 | 16 | 2 | larger-step trend and fail-closed stability check |

This is 64 paired qualifying attempts before any registered replacement.

### Dimensionality matrix

| Row | Geometry | Box scale | Resolution | PPC | Seeds | Claim role |
|---|---|---:|---:|---:|---:|---|
| `F1D` | Q-023 thin 2D3V carrier | 16 | native | 16 | 8 | one-dimensional control only |
| `M2D` | 2D3V | 4 | native | 2 | 8 | matched control for 3D |
| `F3D` | 3D | 4 | native (`512 x 256 x 128`) | 2 | 8 | dimensionality and morphology confirmation |

`F3D` may support a 3D saturation statement only if it independently passes
the box-limit and convergence gates. Otherwise it is a dimensionality control
and the accepted saturation claim remains explicitly 2D3V.

### Finite-rigidity self-consistent branch

This branch is viable only after a new distribution generator and a Q-022
equation/normalization map are reviewed. The exact grid in `k0 r_g0`,
initial CR momentum flux, thermal beta, and distribution shape is
**open**. It must use a separate campaign identifier and cannot reuse the
high-rigidity fixed-current-like interpretation.

The branch must distinguish:

1. CR-current reduction and transverse-momentum growth from CR deflection.
2. Gas acceleration and gas heating.
3. Migration of magnetic power toward `r_g(B)`.
4. Saturated transverse magnetic pressure versus the initial CR momentum-flux
   tensor.

An undriven periodic-box result may be compared only as an undriven analogue.
It may not be described as reproducing the driven-boundary setup of
Zacharegkas et al.

## Required Diagnostics

### Primary time series

- `sqrt(<B_perp^2>)/B0`, `<B_perp^2>/2`, and component-resolved magnetic
  energy.
- Volume-integrated gas momentum, kinetic energy, thermal energy, magnetic
  energy, and total MHD energy.
- CR lab-frame current vector and current parallel to the initial guide field.
- CR bulk velocity relative to the volume-averaged gas.
- CR kinetic energy, momentum vector, and full momentum-flux/pressure tensor.
- Conservative total gas-plus-CR momentum and energy residuals.
- Dominant magnetic wavelength and shell-integrated magnetic spectrum.

### Mechanism discriminants

| Discriminant | Fixed-current-like plasma-acceleration interpretation | CR-deflection / pressure interpretation |
|---|---|---|
| lab-frame `J_CR,parallel` | remains within a preregistered numerical invariance tolerance | decreases materially |
| CR momentum tensor | nearly unchanged | develops transverse momentum/pressure and reduced anisotropy |
| gas parallel velocity | rises and reduces relative drift | may rise, but CR evolution is also material |
| amplified Alfven speed | approaches the relative drift near saturation | not sufficient by itself |
| `r_g(B)/lambda_d` | remains well above unity | approaches order unity as deflection becomes important |
| `P_B,perp,sat / Pi_CR,0` | descriptive only | primary pressure-scaling endpoint after Q-022 mapping |

No numeric threshold in this table is yet a literature-comparison tolerance.
All such thresholds remain open until reference extraction, excluded pilots,
and external review are complete.

### Spectral and morphology diagnostics

- Three-dimensional or two-dimensional FFTs of magnetic-field components,
  with the guide-field mean removed.
- Right/left circular polarization during the inherited linear interval.
- `k_peak(t)`, integral scale, and power in the fundamental box shell.
- Density and magnetic-pressure PDFs; density--magnetic-pressure correlation.
- Quantified filament/cavity measures, not image-only judgments.
- Per-seed morphology panels at preregistered normalized times.

### Required raw evidence

- Frequent built-in MHD history output.
- Raw `mhd_w_bcc` snapshots.
- Deposited `prtcl_rho`, `prtcl_jx`, `prtcl_jy`, `prtcl_jz`, and feedback
  channels at the analysis cadence.
- Particle PVTK output at preregistered sparse checkpoints, including
  positions, velocity, species, source cohort, and macro weight.
- Restarts bracketing the selected nonlinear and saturation windows.
- Stdout runtime identity, terminal proof, Q-017 telemetry, scheduler
  accounting, inventories, and hashes.

For any finite-rigidity relativistic branch, velocity-only float PVTK output
is not automatically an adequate momentum-tensor oracle. That branch must add
or independently validate a momentum-preserving raw diagnostic.

## Fail-Closed Acceptance Logic

The analyzer must emit one of `accepted`, `limited`, `rejected`, or
`not_saturated`. It must never infer success from job completion alone.

### Universal gates

1. Exact clean candidate, executable, deck, analyzer, policy, seed, artifact
   inventory, and observed-time bindings pass.
2. Q-023 linear growth, signed phase, and polarization predecessor closes for
   the selected runtime.
3. The initial Q-023 normalization, current, guide field, pressure, particle
   state, and dimensionless parameters match the frozen deck.
4. No nonfinite output, corrupt artifact, undeclared floor/safety intervention,
   or terminal-time failure occurs.
5. Gas-plus-CR momentum and energy conservation pass preregistered numerical
   tolerances derived from excluded convergence pilots.
6. Saturation is assessed only in a fixed preregistered normalized-time
   window. If no accepted plateau exists by the frozen terminal time, the
   result is `not_saturated`; the window may not be moved after inspection.
7. A magnetic plateau is not accepted as physical saturation if the dominant
   wavelength reaches the fundamental box shell or if the box-doubling row
   fails the frozen agreement criterion.
8. Resolution, PPC, and time-step high-side comparisons pass their frozen
   paired-seed criteria. Finite completed outliers remain in the estimate.
9. The fixed eight-seed median, per-seed values, arithmetic mean, sample
   standard deviation, BCa interval, and inherited Holm family-wise control
   are reported.
10. An independent implementation recomputes the primary table from retained
    raw artifacts.

### Fixed-current-like classification gate

The Riquelme--Spitkovsky plasma-acceleration interpretation is admissible only
if the lab-frame CR current and CR momentum remain invariant within a
preregistered numerical tolerance while gas acceleration reduces the relative
drift. If this gate fails, the result may still be a valid self-consistent
high-rigidity calculation, but it must not be called a fixed-current
saturation result.

### CR-deflection / pressure classification gate

A CR-deflection or anisotropic-pressure interpretation is admissible only if:

- the finite-rigidity source and reference map are frozen;
- CR momentum/current evolution is resolved and converged;
- the pressure tensor is independently recomputed;
- the `r_g(B)/lambda_d` and pressure-scaling endpoints pass their pre-run
  criteria; and
- the claim is scoped to the actual undriven or driven boundary condition.

### Literature-comparison gate

Every numeric comparison tolerance for Riquelme & Spitkovsky, Gargate et al.,
or Zacharegkas et al. remains **open**. Q-022 must first freeze the
equation map, normalization map, parameter overlap, extracted reference data,
extraction uncertainty, and reviewer-approved tolerance rows. Qualitative
mechanism agreement cannot be silently promoted into quantitative
reproduction.

## Resource-Aware Staged Plan

The Frontier cap is shared project-wide:

```text
maximum project total = 10000 node-hours
validated consumed at design time = 1.3019444444444446 node-hours
validated reserved at design time = 0 node-hours
```

Protect a **hard 1500-node-hour Q-019 envelope** before authorizing Q-011
production. This is a planning ceiling, not a scientific tolerance.

| Stage | Work | Hard Q-019 ceiling |
|---|---|---:|
| 0 | source-local generator, analyzer, host/unit tests, Q-023 predecessor closure | no production allocation |
| 1 | excluded window/resource pilots and strong-scaling measurements | 250 node-hours |
| 2 | 64-run paired 2D3V qualifying matrix | 400 node-hours |
| 3 | 1D/matched-2D/3D dimensionality matrix | 600 node-hours |
| 4 | restart carrier, registered replacements, independent recompute, and failure reserve | 250 node-hours |
| total | all mandatory Q-019 work | 1500 node-hours |

### Stage 0: source and analysis closure

Create a new Q-019 problem generator derived from the Q-023 helper math without
mutating the historical Q-023 source or nonlinear foundation. It must validate
the inherited nonlinear-foundation normalization, amplitude, current product,
native cell size, and time convention; create the enlarged periodic boxes;
apply the frozen seeded symmetry-breaking perturbation; preserve restart
state; and emit no qualification claim. Extend the source-local nonlinear
analyzer into the qualifying analyzer and implement its independent recompute
contract before Frontier pilots.

### Stage 1: excluded pilots

Use only pilot seeds `23050091` and `23050092`. These pilots are archived,
explicitly excluded from every qualifying estimate, and may be inspected to
freeze:

- perturbation spectrum and amplitude;
- exact terminal time and saturation window;
- output cadence and storage envelope;
- numerical conservation and convergence thresholds;
- node ladder, memory envelope, and walltime;
- the final resource estimate for every qualifying row.

Run short full-geometry scaling points before long pilots. Include a 3D
The existing source-local `tau=[14,22]` window may guide pilot placement but
must not become the qualifying window without the excluded-pilot freeze.

### Stage 2: core 2D3V production

Run the matrix serially through the registered control plane. Do not inspect
or tune based on partial qualifying output. Reconcile and inventory each
attempt before the next submission.

### Stage 3: dimensionality production

Run `F1D`, `M2D`, and `F3D` only after Stage 2 and the measured project budget
leave the full Stage 3 envelope. If `F3D` is box-limited, retain it as a
dimensionality result and do not claim 3D saturation.

### Stage 4: closure and reserve

Use the reserve only for preregistered restart continuation, scheduler or
artifact failures, and independent recomputation. A scientifically finite
completed outlier is not a failed attempt and cannot be replaced.

### Project-wide authorization equation

Before every Q-019 stage:

```text
C_consumed + C_reserved + C_Q019_remaining
+ C_Q011_authorized + R_other_required <= 10000.
```

If the measured Q-019 core projection exceeds 1500 node-hours, stop and create
a reviewed versioned redesign before any qualifying output is generated. Do
not silently reduce seeds, drop sensitivity rows, or spend the Q-011 reserve.
The finite-rigidity extension is deferred until the high-rigidity core and
other mandatory project claims fit under the cap.

## Open Mappings That Block Execution

1. Exact Q-019 seeded perturbation spectrum and energy.
2. Exact production terminal time, saturation window, output cadence, and
   numerical acceptance thresholds from excluded pilots.
3. Exact AthenaK particle macro-mass, `rho_CR/rho0`, CR energy loading,
   current, and momentum-flux mapping to each independent reference; Q-019
   still requires `PPC*deposit_qscale=2e9`.
4. Exact code-unit magnetic-pressure and beta mapping in every Q-022
   comparison record.
5. Exact finite-rigidity distribution, `k0 r_g0`, CR momentum flux, and
   thermal-beta grid.
6. Whether a new exact locked-current control will be implemented; it is not
   required for the recommended first campaign.
7. Whether a driven CR boundary will be implemented; without it, no direct
   Zacharegkas driven-boundary reproduction is admissible.
8. All literature-derived numeric tolerances and extracted reference values.
9. Clean-candidate, executable, deck, analyzer, policy, Orion root, resource,
   and retention bindings.

## Concrete Next Implementation

The next implementation should create, in new files, a Q-019 high-rigidity
nonlinear problem generator and analyzer that:

1. imports or reuses the Q-023 basis/eigenmode math and preserves the tracked
   nonlinear-foundation contract;
2. rejects any change to the inherited Q-023 current normalization;
3. supports the 16- and 32-wavelength 2D boxes plus the stated 1D/3D controls;
4. freezes deterministic pilot and qualifying seeds;
5. records all mechanism discriminants and conservation budgets; and
6. refuses to emit a fixed-current, CR-deflection, pressure-saturation, or 3D
   saturation claim unless the corresponding gate passes.
