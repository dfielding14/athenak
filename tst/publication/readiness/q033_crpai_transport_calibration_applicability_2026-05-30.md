# Q-033 CRPAI Transport-Calibration Applicability Boundary

## Scope

This source-local note freezes the already implemented extension boundary for a
future adaptive-delta-f physical-damping CRPAI transport-calibration campaign.
It is not physical calibration, qualification evidence, or authorization to
launch a campaign.

The existing `extended_mhd_pic` path can compose:

- physical `kappa_aniso` delta-f with the experimental global bi-kappa refit;
- expanding-box momentum, magnetic-flux, and endpoint feedback maps;
- the reduced high-frequency static-neutral ion-neutral friction map; and
- cell-centered `cc_convert` conservative MHD feedback with Hall mode off.

The exact code boundary is frozen by the paired analyzer's source checksums.
Existing bounded host records cover parser guards, one-cycle manufactured-source
ordering, cadence sensitivity, and restart continuity. They do not establish a
physical transport coefficient.

## Reduced Damping Boundary

The implemented damping map is

```text
p_ion,perp(t + dt) = p_ion,perp(t) exp(-nu_in dt)
nu_in = <particles>/pic_ion_neutral_collision_rate
```

and removes the corresponding transverse ion kinetic energy from ideal MHD.
For expanding-box runs it executes after the endpoint physical-frame feedback
map. This is the bounded static-neutral reduction, not the repository's general
two-fluid ion-neutral module and not an astrophysical damping calibration.

## Candidate Observable Contract

The launch-blocked candidate reserves an analysis contract for future reviewed
runtime artifacts:

```text
anisotropy = abs(xi^2 - 1)
kappa_parallel = slope(var(delta x_parallel)) / 2
nu_eff = v_ref^2 / (3 kappa_parallel)
```

It also requires polarization-resolved wave-power spectra, reviewed tail-window
quasi-steady criteria, and a registered matrix for `nu_eff` scaling versus the
driving rate and `nu_in`.

The paired analyzer exercises those calculations only with normalized synthetic
series. Its synthetic relation and thresholds are parser-oracle fixtures. They
are not predicted CRPAI scaling laws, extracted reference values, accepted
runtime tolerances, or qualifying evidence.

## Launch Block

The Q-033 deck names the intentionally unavailable
`q033_crpai_transport_calibration_open` problem generator, sets `nlim=0` and
`tlim=0.0`, and records `qualification_effect=none`. No local or synthetic
result may claim physical calibration or qualification.

The following remain open:

- Sun-Bai-Zhao applicability review and frozen reference extraction;
- reviewed physical normalization, runtime decks, windows, tolerances, and seed
  policy;
- registered Frontier HIP and MPI campaign matrices;
- immutable runtime artifact inventory and independent recomputation;
- decomposition, resolution, timestep, momentum-bin, pitch-angle-bin, and
  box-size trends;
- external reviewer disposition.
