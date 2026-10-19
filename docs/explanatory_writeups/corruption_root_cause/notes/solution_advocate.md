# Solution Advocate

## Strongest Fair Case

The shared-buffer explanation is unusually strong because it is predictive,
not merely compatible with the archive. The destructive
`Kokkos::realloc(derived_var, ...)` appears exactly inside
`coord_abscostheta`, after some coordinate rows have already been computed.
That predicts order-dependent survival:

- radius-first streams lose radius and collapse into radius underflow;
- `output29` computes the destructive angle first, shrinking the first extent
  from `nmb_max` to `nmb`; the next radius call triggers the common capacity
  realloc back to `nmb_max`, erases angle, and then writes radius;
- `output31` never requests `coord_abscostheta`, so it remains an independent
  control;
- totals can remain nearly correct because the weights are still accumulated,
  just into the wrong bins.

Traceable source history:

- regression introduced by `ec009a27b87029db6f35a95045ae624ff5c0d1e4`
- destructive allocation removed by `db953863b6c0e06d893b4a9738d4c111279e42fe`

The real-data comparison has the same mixed signature. Full rebuilt
`output29` appropriately disagrees with the broken archive
(`E_L1 = 1.9176`), while its radial marginal agrees
(`E_L1 = 1.60e-11`). The intact `output31` control agrees closely
(`E_L1 = 3.87e-7` full 4D), and signed totals agree near `1e-12`.

## Claim Boundary

This strongly supports the root cause of the archived PDF pattern. It does not
replace a same-state old-line/new-line AthenaK A/B reproduction, prove that no
other PDF-output bug exists, or certify complete-snapshot AMR closure.
