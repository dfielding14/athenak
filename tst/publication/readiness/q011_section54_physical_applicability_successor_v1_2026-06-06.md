# Q011 Section 5.4 physical-applicability successor v1

Status: **source-local diagnostic/gate design only; no launch authority; no
policy authority; no qualifying-output inspection authority; no claim
closure; no production result**

## Purpose

The existing Q011 Section 5.4 campaign can reproduce the specified no-Hall
AthenaK equations. That numerical reproduction is not sufficient to establish
that Hall omission, the MHD treatment, or the finite transverse box remain
physically applicable throughout shock acceleration.

This additive successor supplies the missing fail-closed analysis boundary. It
does not modify historical Q011 reducers, campaign contracts, decks, source, or
execution policy. It consumes already decoded matched production-science
snapshots and a future runtime time/escape record.

## Exact normalization boundary

`Q011-APP-NORM` requires an exact supplied mapping. The reducer does not infer
normalization from a deck or generator:

- one selected positive-ion species, index zero;
- `|q/m|=1`, `rho0=1`, `B0=1`, and particle light speed `C=10000`;
- gas density is `rho_g`;
- `prtcl_rho` is the code-normalized deposited `rho_CR/c`, named `rho_q`;
- `prtcl_j*` is the code-normalized deposited `J_CR/c`, named `J_q`;
- `v_A=|B|/sqrt(rho_g)`;
- `J_q,gas=J_q-rho_q v_g`;
- particle momentum is `p=gamma(v)v`;
- the conservative containment gyroradius is `r_g*=|p|/(|q/m||B_local|)`.

Missing, additional, type-drifted, or value-drifted fields fail before any
applicability record is emitted.

## Snapshot diagnostics

The reducer reuses the existing Q011 production-science matched-grid decoder.
For every cell it computes

```text
R       = rho_q / (rho_g + rho_q)
Lambda  = |J_q - rho_q v_g| / ((rho_g + rho_q) v_A)
d_i     = rho_g^(-1/2)
S_delta = min(actual leaf dx1, actual leaf dx2) / d_i
```

The immutable in-memory result contains per-cell `R`, `Lambda`, `d_i`,
`S_delta`, actual leaf spacings, and gas-frame current magnitude maps.
Physical-applicability evidence must retain those maps with the bound grid
faces; summary JSON alone is insufficient. Per-cell maxima are mandatory and
there is no rare-cell waiver.

Detected-front-centered summaries are produced for the full domain,
downstream `[-1200,-120]`, precursor `[120,1200]`, and far-upstream
`[1200,2400]` offsets in `c/omega_pi`. Each region reports area-weighted and
`|J_q,gas|`-weighted statistics. A region with exactly zero current reports an
explicit unavailable current-weighted statistic rather than inventing one.

## Ion-scale separation

The `S_delta` test uses each composite cell's physical source level to recover
the actual leaf spacing. The intentionally unresolved shock-transition layer,
defined here as the detected front plus or minus `120 c/omega_pi`, is excluded
from ion-scale qualification. Microscopic shock-structure and self-consistent
injection claims remain prohibited regardless of gate outcome.

The precursor magnetic characteristic scale uses a two-dimensional,
Hann-windowed vector `delta B` spectrum. Its characteristic wavelength is
`2*pi` divided by the magnetic-power-weighted mean wavenumber. The sub-`10d_i`
fraction is the fraction of precursor magnetic fluctuation power at
wavelengths shorter than ten times the largest local precursor `d_i`. Absent
resolvable precursor fluctuation power fails closed.

## Particle containment

Active shock-injected positive-weight particles are selected using
`cr_source==1` and `birth_time>=45`. The reducer reconstructs
`p=gamma(v)v`, samples the matched finest-composite magnetic field with TSC,
and reports macro-weighted and CR-energy-weighted gyroradius quantiles,
maximum gyroradius, and the energy fraction above `Ly/4`.

The x2 TSC stencil is periodic. The nonperiodic x1 stencil must remain inside
the retained state or the reducer fails. At least 1000 positive-weight active
particles are required to report q999.

The composite-field TSC sample is a snapshot applicability diagnostic. Exact
runtime pusher-field exposure remains a future runtime-evidence obligation,
especially for AMR particles near refinement interfaces.

The same TSC stencil samples `R` and `Lambda` at every selected particle. The
successor reports macro-weighted and CR-kinetic-energy-weighted means,
quantiles, and threshold-exceedance fractions for all active particles,
detected-front upstream particles, and a fixed high-energy tail defined by the
macro-weighted q990 specific kinetic energy.

## Conservative gates

The numerical thresholds below are **AthenaK-selected**, not numeric
tolerances prescribed by the cited literature.

| Gate | Acceptance |
| --- | --- |
| `Q011-APP-NORM` | exact normalization and representation mapping |
| `Q011-APP-R` | full-domain per-cell `Rmax <= 0.01` |
| `Q011-APP-LAMBDA` | full-domain per-cell `Lambdamax <= 0.1` |
| `Q011-APP-DI` | outside the shock transition, `S_delta,min >= 1`; precursor `lambda_B,char/d_i,max >= 10`; precursor sub-`10d_i` power fraction `<= 0.05` |
| `Q011-APP-RG` | macro and energy q999 `<= Ly/8`; maximum `<= Ly/2`; CR-energy fraction with `r_g*>Ly/4 <= 1e-3` |
| `Q011-APP-TIME` | complete contiguous every-cycle post-startup coverage through `t=1200`, complete particle exposure including pre-destruction events, and a closed boundary-escape ledger with no unaccounted particles, weight, or energy |

## Literature relationship

- Bai et al. (2015) requires a dilute CR population, identifies the CR-Hall
  parameter `Lambda`, states that no-Hall behavior requires `Lambda << 1`, and
  requires MHD-PIC interpretation on scales much larger than the ion inertial
  length. It does not prescribe `0.01`, `0.1`, or a decade as universal
  thresholds. The selected `Rmax` bound limits the corresponding density
  correction to approximately one percent. At `Lambda=0.1`, the linear
  correction factor `1+(Lambda/2)^2` is `1.0025`.
- Sun & Bai Section 5.4 states that the transverse domain is large enough to
  contain several high-energy CR gyroradii. It does not prescribe a numeric
  containment tolerance. The q999, maximum, and high-gyroradius energy-fraction
  limits are conservative AthenaK interpretations.

Primary references:

- Bai, Caprioli, Sironi & Spitkovsky 2015, ApJ 809, 55,
  DOI `10.1088/0004-637X/809/1/55`, arXiv `1412.1087`.
- Sun & Bai 2023, MNRAS, DOI `10.1093/mnras/stad1548`,
  arXiv `2304.10568`.

## Runtime time/escape interface

The history reducer requires a future
`q011_section54_runtime_time_escape_applicability_v1` record with exact keys:

- exact normalization and non-authorizing fields;
- cycle coverage from post-startup removal at `t=45` through `t=1200`,
  including first/last cycle, exact contiguous count, gap count, and complete
  restart-segment coverage;
- all-cycle extrema for every `R`, `Lambda`, ion-scale, and gyroradius gate
  observable; those extrema must conservatively dominate every supplied
  snapshot observable;
- complete particle exposure from every active update and every
  pre-destruction boundary event, including escaped particles, with sampled
  `R`/`Lambda` maxima and cumulative macro/energy-weighted exceedance
  fractions consistent with the all-cycle cell extrema;
- complete `ix1` and `ox1` escape ledgers, explicit periodic x2 faces, exact
  face-to-total count/weight/energy closure, zero unaccounted quantities, and
  zero source-census residuals.

The interface permits accounted boundary escape. It does not permit silently
destroyed particles or an exposure record that omits escaped particles. No
runtime telemetry conforming to this interface exists yet.

## Claim rejection rules

- `APP-NORM` or `APP-R` failure rejects all physical MHD-PIC Bell, shock, and
  DSA claims.
- `APP-LAMBDA` failure rejects Hall-negligible, Bell-mechanism, physical
  magnetic-amplification, and physical DSA-scattering claims.
- Any observed `Lambda>=1` additionally rejects the claim that the no-Hall run
  approximates the target plasma.
- `APP-DI` failure rejects physical precursor-turbulence and scattering
  interpretation.
- `APP-RG` failure rejects `Emax`, high-energy slope/cutoff, acceleration-rate,
  and acceleration-efficiency claims.
- Missing or incomplete `APP-TIME` evidence rejects all history/global
  applicability claims.
- Microscopic shock-structure and self-consistent injection claims are always
  excluded.

Passing every gate still grants no launch, policy, inspection, or claim
authority.

## Open production obligations

1. Implement and independently review the runtime all-cycle exposure and
   boundary-escape ledger in the separate C++ workstream.
2. Bind retained raw snapshots and returned cell maps to immutable production
   attempt provenance.
3. Confirm the composite TSC gyroradius diagnostic against runtime sampled
   particle fields, especially across AMR interfaces.
4. Execute the registered production campaign and evaluate this successor
   without changing thresholds after inspecting results.
5. Independently review physical interpretation and claim scope.
