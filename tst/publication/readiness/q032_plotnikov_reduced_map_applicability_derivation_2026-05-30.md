# Q-032 Plotnikov Reduced-Map Applicability Derivation

## Disposition

This is a bounded source-local applicability derivation for the existing
`pic_wave_damping_mode=ion_neutral_friction` map. It is not a Plotnikov-matched
equation map, extracted reference dataset, tolerance table, campaign
authorization, or qualification artifact.

```text
artifact_role = source_local_applicability_derivation_only
qualification_effect = none
claim_closure = false
plotnikov_matched_qualification = false
```

The Q-022 comparison records remain authoritative for any future matched
comparison. This note deliberately does not populate their empty equation,
normalization, parameter-overlap, excluded-regime, dataset-extraction, or
numeric-tolerance fields.

## Source Anchors

The implementation boundary is frozen separately by
`q032_reduced_static_neutral_plotnikov_boundary_2026-05-30.md` and the Q-032
readiness sidecars. In `src/mhd/mhd_tasks.cpp`,
`MHD::ApplyPICWaveDampingMap` applies one final-stage source map with

```text
D_src = exp(-nu_in dt)
nu_in = <particles>/pic_ion_neutral_collision_rate
```

to the active-MHD conserved state.

The private Orion archive is registered by
`q022_external_reference_private_ingest_2026-05-30.json`. Its
`plotnikov_ostriker_bai_2021_arxiv_2102.11878` entry identifies
`arxiv_2102.11878.pdf` with source locator `arXiv:2102.11878`. That PDF was read
as a source reference only. No figure, table, or runtime value has been
extracted from it for this artifact.

The archived paper's Section 2.2.1 and Section 3.1 are reading anchors, not a
registered AthenaK-to-literature equation map. They distinguish the
high-frequency one-fluid wave asymptote from the per-step transverse ion
momentum update. This note preserves that distinction without claiming that a
future AthenaK campaign satisfies the paper's applicability conditions.

## Implemented Source Substep

For one non-expanding endpoint source application, AthenaK leaves density,
longitudinal momentum, and magnetic fields unchanged by this map and updates

```text
M1' = M1
M2' = D_src M2
M3' = D_src M3
rho' = rho
E' = E - (1/2) (1 - D_src^2) (M2^2 + M3^2) / rho
```

for ideal MHD. Because density is unchanged during the substep,

```text
u_ion,perp' = exp(-nu_in dt) u_ion,perp
K_ion,perp' = exp(-2 nu_in dt) K_ion,perp
```

for that source operation alone. The direct factor `D_src` is a transverse ion
momentum or velocity sink. It is not, by itself, an Alfven-wave amplitude
attenuation factor or a magnetic wave-energy attenuation factor.

The local Q-032 runtime carrier is intentionally narrower still: it uses no
simulation particles, no particle feedback, and no Plotnikov comparison. Its
control-versus-damped result is a source-map mechanics regression only.

## Conditional Wave-Mode Asymptote

The archived reference describes a static-neutral one-fluid limit in which
neutral dynamics are omitted and transverse ion friction enters the linear
Alfven-wave problem. In that conditional limit, the source-local linear
oscillator has the form

```text
omega^2 + i nu_in omega = k^2 V_A,i^2
```

and, only when the decoupled high-frequency condition is satisfied,

```text
k V_A,i >> max(nu_in, nu_ni)
```

the wave-mode envelope has

```text
-Im(omega) = Gamma_d = nu_in / 2
D_wave_amp(t) = exp(-nu_in t / 2)
D_wave_energy(t) = exp(-nu_in t)
```

where `D_wave_energy` applies to wave energy or amplitude-squared. The
observable roles must remain distinct even where the factors are algebraically
related:

```text
for nu_in dt > 0:
D_src(dt) = exp(-nu_in dt)
D_wave_amp(dt) = sqrt(D_src(dt))
D_wave_energy(dt) = D_src(dt)
D_ion_perp_kinetic_energy(dt) = D_src(dt)^2
```

The equality `D_wave_energy(dt) = D_src(dt)` is a conditional algebraic
coincidence between different observables, not permission to identify the
source sink with magnetic wave energy. The direct source substep leaves the
magnetic field unchanged. A matched campaign needs a reviewed eigenmode,
normalization, and observable mapping before it may compare an AthenaK
measurement to a literature damping asymptote.

## Unresolved Applicability

The following items remain open and fail closed:

- Static-neutral reduction: AthenaK applies a single-fluid active-MHD
  transverse source sink and does not evolve a neutral fluid in this path.
  Review whether the intended comparison stays inside the decoupled
  static-neutral regime.
- High-frequency condition: the conditional asymptote requires
  `k V_A,i >> max(nu_in, nu_ni)`. Freeze the candidate wavenumber range, units,
  collision-frequency normalization, and excluded low-frequency regimes.
- Longitudinal drag: the AthenaK reduced map leaves `M1` unchanged. Review the
  intentional reduction against the selected observables and comparison scope.
- Wave observable: the source map directly updates transverse ion momentum;
  its same-substep magnetic field is unchanged. Freeze whether amplitude,
  amplitude-squared, magnetic wave energy, growth rate, bandwidth, saturation,
  or isotropization is measured and how it is normalized.
- CRSI setup: the bounded runtime carrier has no simulation particles and no
  CRSI feedback. Freeze a separately reviewed true-CRSI generator,
  distribution, seed-wave, resolution, timestep, and sensitivity matrix.
- External reference data: the archived PDF is a source reference only.
  Perform authorized extraction with provenance and uncertainty before creating
  reviewed tolerance rows.
- Review and execution: Q-032 remains open. Obtain external review and separate
  registered MPI and Frontier authorization before any qualifying execution.

No extracted reference values, matched equation mappings, or numeric
tolerances are supplied by this artifact.

## Particle Phase Scrambling

The archived reference's numerical-method summary states that particles are
phase-scrambled when they cross the periodic system boundary and re-enter, to
mimic a larger numerical box. The current Q-032 source-local records do not
bind or demonstrate an AthenaK counterpart and do not establish that omission
is acceptable for a matched campaign.

```text
particle_phase_scrambling_for_source_local_sink_derivation = not_required
particle_phase_scrambling_for_plotnikov_matched_campaign = unresolved_fail_closed
```

Particle phase scrambling is not required to derive or regress the
particle-free active-MHD source sink. Before any Plotnikov-matched CRSI claim,
the campaign must either bind and validate the required phase-scrambling
behavior or obtain a reviewed, scoped justification for omitting it.

## Fail-Closed Boundary

The existing Q-022 placeholders intentionally remain blocked:

```text
q022_dataset_status = blocked_extraction_input_unavailable
q022_equation_map_status = blocked_pending_reference_specific_mapping_and_external_review
q022_numeric_tolerance_rows = intentionally_empty
```

Nothing in this derivation promotes the bounded local source-map regression,
the conditional algebra above, or the private source-reference archive into
Plotnikov-matched damped-CRSI evidence.
