# Q-029 Experimental Hall Source Normalization

## Scope

This page freezes the code-level normalization of the existing opt-in
`current_to_ct_experimental` source. It is a source-local preparation artifact,
not Hall-Bell qualification evidence and not a physical literature map.

The implementation in `src/mhd/mhd_tasks.cpp` applies

```text
cE_CT = cE_ideal + alpha_H P_edge[J_CR]
alpha_H = <particles>/couple_j_to_efield_coeff
```

when `Particles::AddsCRCurrentToCT()` admits the separately named
`extended_mhd_pic` mode and coupled particle moments. `P_edge` is the selected
cell-centered-to-edge conversion or the explicit staggered-current
representation.

## Dimensionless Form

Choose reference values

```text
rho0 > 0
B_g0 > 0
U_A0 = B_g0 / sqrt(rho0)
cE0 = U_A0 B_g0
J_CR0 > 0
```

and normalized fields

```text
hat(cE) = cE / cE0
hat(J_CR) = J_CR / J_CR0
chi_H = alpha_H J_CR0 / (U_A0 B_g0).
```

The exact implemented source becomes

```text
hat(cE_CT) = hat(cE_ideal) + chi_H P_edge[hat(J_CR)].
```

The source-local candidate deck freezes

```text
rho0 = 1
B_g0 = 1
U_A0 = 1
cE0 = 1
J_CR0 = 1
alpha_H,fiducial = 0.5
chi_H,fiducial = 0.5
chi_H candidate grid = {0, 0.25, 0.5, 1}
signed source-oracle alpha_H grid = {-1, 0, 1}.
```

The signed grid exists only to verify the additive source contract and its
odd-in-coefficient increment. Negative `alpha_H` is not registered as a
physical extension campaign point.

## Fail-Closed Boundary

The Q-029 candidate deck names the intentionally unavailable
`q029_extended_hall_normalization_open` problem generator, sets `nlim=0`, and
records `qualification_effect=none`. The paired analyzer accepts only a
synthetic normalization-contract bundle and always emits
`qualifying_evidence=false`.

The following remain open and cannot be inferred from this source-local
preparation:

- the reviewed physical CR-Hall mapping to Bai et al. (2015);
- linear or nonlinear Bell behavior;
- shock-front applicability;
- decomposition or MPI qualification;
- Frontier HIP or GPU qualification;
- external-review closure.
