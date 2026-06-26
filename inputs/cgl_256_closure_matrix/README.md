# CGL 256^3 Closure Matrix

This matrix isolates two fixed-input effects in active beta-10 CGL turbulence:

1. uniform background pressure-isotropization at zero heat flux; and
2. collisionless Landau-fluid heat flux relative to zero heat flux.

| Case | Uniform `nu_coll` | LF heat flux | `lf_k_parallel` |
|---|---:|---|---:|
| `cgl256_b10_nolf_nu0` | 0 | off | N/A |
| `cgl256_b10_nolf_nu10` | 10 | off | N/A |
| `cgl256_b10_lf_k2pi_nu0` | 0 | on | `2*pi` |
| `cgl256_b10_lf_kpi_nu0` | 0 | on | `pi` |

All cases use a periodic `256^3`, `L=1` box, active CGL pressure feedback,
`B0=1`, initial ion beta 10, seed `271828`, and the same perpendicular
Alfvenic forcing realization. The total injection target is `dedt=0.16` and
the OU correlation time is 1, obtained from the R17 values by accounting for
the halved box volume and parallel length.

The mirror/parallel-firehose hard-wall limiter is held fixed in every case.
Thus `nu_coll=0` means no uniform background collisions; it does not mean that
threshold-local microinstability regularization has been removed. This fixed
limiter policy prevents a closure or collision comparison from being
confounded by different admissibility handling.
