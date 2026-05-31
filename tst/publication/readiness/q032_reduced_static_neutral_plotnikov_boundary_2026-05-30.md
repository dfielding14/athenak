# Q-032 Reduced Static-Neutral Friction Boundary

## Scope

This page freezes the code-level boundary of the existing opt-in
`pic_wave_damping_mode=ion_neutral_friction` map. It is source-local
preparation for a separately named Plotnikov-comparison candidate. It is not
ion-neutral-damped CRSI qualification and does not supply a literature map.

The selected PIC extension is separate from the repository's general
two-fluid `src/ion-neutral/` module. It does not evolve neutral density,
neutral momentum, ionization, recombination, or backreaction onto a neutral
fluid.

## Implemented Map

For a non-expanding run, the implemented endpoint map is applied once after
the final explicit RK source update. Define

```text
D = exp(-nu_in dt)
nu_in = <particles>/pic_ion_neutral_collision_rate > 0.
```

The code updates the active-MHD conserved state as

```text
M1' = M1
M2' = D M2
M3' = D M3
rho' = rho
E' = E - (1/2) (1 - D^2) (M2^2 + M3^2) / rho
```

for ideal MHD. Magnetic fields are not modified by this map. Expanding-box
runs route the same friction map after the endpoint physical-frame feedback
map. The source-local Q-032 candidate intentionally keeps expanding-box mode
off so it freezes only the direct non-expanding route.

## Applicability Boundary

The implementation is a bounded high-frequency, static-neutral reduction. It
does not by itself establish:

- a reference-specific equation or normalization map to Plotnikov, Ostriker
  and Bai (2021), arXiv:2102.11878;
- an extracted reference dataset with uncertainty;
- reviewed numerical tolerance rows;
- unstable-bandwidth, growth-rate, saturation, or isotropization agreement;
- MPI decomposition parity or Frontier HIP/GPU parity;
- external-review closure.

The existing Q022 comparison records intentionally retain these items as open.
The Q-032 synthetic analyzer checks only the exact endpoint algebra above and
requires those Plotnikov-comparison records to remain fail-closed.

## Launch Block

The Q-032 deck names the intentionally unavailable
`q032_plotnikov_damped_crsi_open` problem generator, sets `nlim=0`, and records
`qualification_effect=none`. The paired analyzer accepts only a synthetic
endpoint-contract bundle and always emits both
`qualifying_evidence=false` and `matched_plotnikov_qualification=false`.

No synthetic or bounded-local result may be promoted into a matched
ion-neutral-damped CRSI qualification claim.
