# Q011 Section 5.4 pressure-selection options

Status: **human review memo only; no pressure selected; no execution
authorization; no science claim**

## Immutable evidence basis

- Source commit:
  `ce41d4b29bc646b4f0740468e1026f7308b29026`
- Aggregate receipt:
  `/lustre/orion/ast207/proj-shared/dfielding/PIC/publication/q011_section54_pressure_pilot_bundle_receipt.json`
  with SHA-256
  `9117b3dbc7573187b2d080568e69bdbbee0642f2a965aa543273ab3ea3d67be9`
- Review-packet receipt:
  `/lustre/orion/ast207/proj-shared/dfielding/PIC/publication/q011_section54_pressure_pilot_review_packet_receipt.json`
  with SHA-256
  `3f20d3d26a479aa508439f9d038ec6510643bf407aa081fae22959a57571de5d`
- Review-packet inventory SHA-256:
  `ba38e86575baee720871e1f624e20edb5482f9a4a00afb4286ee136f821112c1`

Both publication receipts passed their documented read-only consumer
verification after worker publication. These four short pilots remain
engineering calibration only, not Section 5.4 reproduction evidence.

## Physical interpretation

All four cases hold fixed the prescribed parameters that define the intended
parallel-shock and current-driven setup: upstream density and magnetic field,
inflow speed `u0=30`, nominal inflow-frame Alfvenic Mach number `M_A=30`,
injection efficiency, injected-particle momentum, numerical light speed,
geometry, and resolution. The only explicit case override is upstream gas
pressure. Because the code defines the fail-closed injection-subtraction guard
as a fixed fraction of `p0`, its absolute threshold also changes
proportionally; the guard fraction remains fixed and the retained minima stay
roughly eight orders of magnitude above it. This guard aborts before clipping
and is distinct from AthenaK's separately configured EOS pressure floor.

For `rho0=B0=1` and `gamma=5/3`:

| Rank | Case | `p0` | `beta0=2p0/B0^2` | `cs` | nominal inflow-frame sonic Mach | nominal inflow-frame parallel fast Mach | `p0/(rho0 u0^2)` |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | `ps_p0_1p00` | 1.00 | 2.0 | 1.29099 | 23.2379 | 23.2379 | 1/900 |
| 2 | `ps_p0_0p20` | 0.20 | 0.4 | 0.57735 | 51.9615 | 30.0 | 1/4500 |
| 3 | `ps_p0_0p10` | 0.10 | 0.2 | 0.40825 | 73.4847 | 30.0 | 1/9000 |
| 4 | `ps_p0_0p05` | 0.05 | 0.1 | 0.28868 | 103.9230 | 30.0 | 1/18000 |

These Mach values use the prescribed inflow speed `u0`, not the upstream speed
relative to the moving shock. They therefore should not be read as measured
shock-frame Mach numbers or exact matches to a literature shock-frame Mach
convention.

Every option is a strong-shock candidate by thermal-to-inflow-ram-pressure
ratio. Lowering `p0` changes plasma beta and nominal sonic Mach number, but
does not change the fixed nominal Alfvenic Mach number or the injection
prescription. The realized CR current that drives Bell growth is an evolved
quantity that depends on particle transport, scattering, and backreaction;
these pilots do not establish that it is pressure independent.

## Literature relationship

1. Sun and Bai (2023), Section 5.4, uses a high-Mach parallel shock with an
   initial `M_A=30` and fixed artificial injection. It says the setup is
   similar to Bai et al. (2015), but does not specify the upstream gas
   pressure.
2. Bai et al. (2015) explicitly sets `T0=P0=1` and states that this choice is
   irrelevant while `T0 << v0^2`. Among the directly linked setup papers
   checked here, it provides the explicit upstream pressure value and therefore
   gives `p0=1` the strongest provenance.
3. Bell (2004) identifies CR current, magnetic field, and density as primary
   controls in its idealized fixed-current linear analysis. Later
   hybrid-shock work also emphasizes shock speed, Alfvenic Mach number, and
   the self-consistently evolved accelerated-particle population. The
   prescribed upstream, inflow, and injection parameters are fixed across
   these pressure pilots, but the realized CR current is not prescribed or
   demonstrated to be identical.
4. Caprioli and Spitkovsky (2014) shows that strong-shock magnetic
   amplification and the resonant-to-non-resonant transition depend strongly
   on shock strength and the accelerated-particle population. The present
   MHD-PIC setup prescribes injection, so a colder upstream does not make the
   injection microphysics more self-consistent.

Primary references:

- Sun & Bai 2023, MNRAS, DOI `10.1093/mnras/stad1548`,
  arXiv `2304.10568`.
- Bai, Caprioli, Sironi & Spitkovsky 2015, ApJ 809, 55,
  DOI `10.1088/0004-637X/809/1/55`, arXiv `1412.1087`.
- Bell 2004, MNRAS 353, 550,
  DOI `10.1111/j.1365-2966.2004.08097.x`.
- Caprioli & Spitkovsky 2014, ApJ 794, 46,
  DOI `10.1088/0004-637X/794/1/46`, arXiv `1401.7679`.

## Common read-only findings

- All cases reach the terminal pilot time at cycle `874`.
- The prescribed ideal injection surface is at `x=600`; the retained
  grid-scale diagnostics do not establish the physical shock position exactly
  at that coordinate.
- The maximum of each terminal `t=60` full-domain transverse-averaged density
  profile, divided by `rho0`, is approximately `4.073-4.086`, close to the
  strong-shock value of four.
- The terminal density, magnetic-field, and PIC-current panels are
  qualitatively similar across all four cases. This visual comparison is not
  a quantitative demonstration that the realized CR currents are identical.
- Terminal particle counts differ by only eight particles:
  `142536` for `p0=1` and `142528` for each colder case.
- The load-balancing diagnostic `load.particle_efficiency` is `1.0` for every
  case; this is not an injection or particle-acceleration efficiency.
- Tracked GPU memory differs by less than `0.004%`.
- Minimum normalized pressure remains positive in every case and more than
  `9.8957e7` times the configured proportional injection-subtraction guard.
  The minimum `p/p0` values are `0.99596`, `0.99377`, `0.99226`, and `0.98957`
  in ranked order.
- Invalid particle records, retained terminal early-cohort particles, and the
  printed source-transaction residual are exactly zero in every case.
- The largest terminal global magnetic amplification is only
  `Bmax/B0=1.00479`; these short pilots do not contain developed Bell
  amplification.
- Maximum normalized profile RMS differences between adjacent ranked cases
  are:

  | Pair | Maximum normalized RMS difference |
  | --- | ---: |
  | `p0=1.0` versus `p0=0.2` | `0.0319738` |
  | `p0=0.2` versus `p0=0.1` | `0.00406863` |
  | `p0=0.1` versus `p0=0.05` | `0.00204834` |

For this memo-side recomputation, each comparison uses the terminal `t=60`
full-domain transverse-averaged one-dimensional profiles. It forms the
pointwise difference between two cases, computes its RMS over all retained
`x` cells, normalizes density by `rho0=1`, gas pressure by
`rho0*u0^2=900`, velocity by `u0=30`, and magnetic field by `B0=1`, then
reports the maximum of those normalized component RMS values. The terminal
transverse-profile-maximum density range and these RMS summaries are derived
from immutable retained profiles but are not separately bound outputs of the
aggregate analyzer or review packet.

These are unaligned wall-frame profile comparisons. Changing `p0` changes the
finite-Mach shock propagation, so an accumulating front displacement can
increase pointwise density RMS even when the local shock structures remain
similar. The RMS values therefore cannot trigger a move-down decision by
themselves; any pressure-linked material-difference claim needs a shock-aligned
or region-based recomputation.

The normalized-profile comparisons and the `0.05` material-difference band
below are post-hoc advisory engineering heuristics, not literature-derived
tolerances or preregistered scientific acceptance criteria. Each pressure has
one retained deterministic realization; there are no repeated seeds,
uncertainty estimates, or statistical evidence that the observed differences
are pressure-linked. The comparisons therefore support only a cautious
early-time operational ranking.

These pilots do not run long enough to choose between cases using late Bell
amplification, diffusive-shock-acceleration spectra, or the target `t=500` and
`t=1200` Section 5.4 diagnostics. They also provide no AMR or
restart-continuation evidence: AMR was explicitly disabled and no restart
continuation was exercised.
The retained `prtcl_jx` panels are qualitatively similar, but no bound
quantitative current-profile metric is included in this ranking. A future
Bell-driving diagnostic must define particle current relative to the upstream
gas, include the coupled particle charge density, and distinguish precursor
and far-upstream populations. The present evidence therefore supports a
provenance-first engineering recommendation, not a physics-preferred Bell or
DSA baseline.

## Ranked option exploration

### 1. `p0=1.0`: provenance-first baseline

Pros:

- It is the explicit upstream pressure used by Bai et al. (2015).
- It has the largest absolute thermal-pressure reserve against subtraction,
  cancellation, and roundoff.
- It is still strongly ram-pressure dominated and has high nominal
  inflow-frame Mach numbers.
- It is the fastest retained pilot by both zone cycles/s and particle
  updates/s.
- No retained early engineering metric shows a defect relative to the colder
  cases under the post-hoc advisory heuristic.

Cons:

- It is the warmest option, with `beta0=2` and nominal inflow-frame parallel
  fast Mach `23.24` instead of `30`.
- Its terminal density profile has the largest difference from the colder
  cluster, although it remains below the post-hoc advisory band.

Exploration result:

- The observed `0.03197` maximum normalized profile difference is below the
  post-hoc advisory material-difference band of `0.05`; this is not a
  scientific pass/fail threshold.
- The largest component of that difference is density RMS; pressure,
  velocity, and magnetic RMS differences are respectively `0.01479`,
  `0.01119`, and approximately `1.0e-7`.
- Given the evidence limitations above, the retained packet provides no
  operational reason to leave the literature-anchored baseline.

### 2. `p0=0.2`: first cold-shock alternative

Pros:

- It is the least aggressive option with `cs < vA`, giving nominal
  inflow-frame parallel fast Mach `30`.
- It retains substantially more pressure reserve than `p0=0.1` or `0.05`.
- It lies close to the colder cases in the retained profiles.

Cons:

- It has no direct pressure-value support in the principal literature.
- No retained early engineering metric demonstrates an advantage over
  `p0=1`; the pilots cannot establish scientific superiority.
- Moving to it would trade provenance and pressure reserve for a nominal
  inflow-frame fast-Mach match that Sun and Bai do not explicitly require.

Exploration result:

- This remains a viable fallback if later evidence demonstrates a
  pressure-linked defect in `p0=1`.
- The current immutable packet does not demonstrate such a defect.

### 3. `p0=0.1`: aggressive cold-shock alternative

Pros:

- It lies in the observed early-time colder-case profile cluster: its maximum
  normalized difference from `p0=0.05` is `0.00205`. This does not establish
  a cold-case limit or convergence.
- It gives high nominal inflow-frame sonic Mach while retaining more pressure
  than `p0=0.05`.

Cons:

- It has half the pressure reserve of `p0=0.2` with no demonstrated advantage
  in the retained early engineering diagnostics.
- Higher nominal sonic Mach alone is not evidence of a better Section 5.4
  reproduction. The nominal Alfvenic Mach number and injection prescription
  are unchanged, but equality of the realized CR current is not established.
- It is the slowest retained pilot.

Exploration result:

- It is a viable diagnostic fallback, not a preferred qualifying baseline.
- There is no retained evidence supporting a move from `p0=0.2` to `0.1`.

### 4. `p0=0.05`: numerical stress endpoint

Pros:

- It provides the coldest available stress test.
- It provides a third colder early-time profile for comparison; it does not
  establish physical or numerical convergence.

Cons:

- It has the least absolute pressure reserve and greatest sensitivity to
  subtraction, cancellation, and roundoff. The pilots do not show EOS
  pressure-floor activity.
- No colder case exists to corroborate a pressure-linked anomaly.
- It provides no distinct literature or retained early engineering-evidence
  advantage.

Exploration result:

- Treat it as a useful stress endpoint, not as a primary qualifying choice.
- Any defect unique to this endpoint should block or motivate a separately
  authorized investigation, not force selection of this case.

## Decision rule

Advisory provenance-first ranking: `p0=1.0`, then `0.2`, then `0.1`, then
`0.05`.

The move-down rule, material-difference band, and factor-of-two criterion are
post-hoc operational heuristics. They are not literature-derived scientific
thresholds, and the single retained realization per pressure supplies no
uncertainty estimate.

Move down only when the current option has a repeated, pressure-linked failure
or material difference confirmed with shock-aligned or region-based metrics,
the next option satisfies the post-hoc advisory pass band or improves the
failure by at least a factor of two, and no new lower-pressure safety
escalation appears. A larger unaligned wall-frame RMS, a larger early magnetic
fluctuation, a single throughput measurement, or a higher nominal sonic Mach
number is not sufficient move-down evidence. A physics-motivated move also
requires quantitative diagnostics of the realized accelerated-particle
current relative to the upstream gas, the coupled particle charge density,
precursor and far-upstream populations, and late-time Bell/DSA behavior.

Advisory recommendation: retain `p0=1.0` as the qualifying baseline. This
recommendation is not a human pressure-selection receipt and must not be
consumed as one.
