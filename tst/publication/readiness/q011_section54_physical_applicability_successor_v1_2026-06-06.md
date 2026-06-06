# Q011 Section 5.4 physical-applicability successor v1

Status: **clean repair successor to independently rejected commits
`4fd567f5a4db0662c765668d9cb56aa63be22a3c` and `72792e172`; additive
source-local diagnostic/gate only; source base remains rejected; no production
telemetry or result; no authorization; fresh independent review required**

## Purpose and boundary

This Python successor determines whether a future Q011 non-relativistic shock
run remains inside the physical applicability domain required for no-Hall
MHD-PIC and finite-box acceleration claims. It does not launch a run, mutate
policy, inspect protected qualifying output, or close a claim. Passing every
gate has no authorization effect.

The reducer preserves the Bai et al. normalization and formulas:

```text
R       = rho_q / (rho_g + rho_q)
v_A     = |B| / sqrt(rho_g)
Lambda  = |J_q - rho_q v_g| / ((rho_g + rho_q) v_A)
d_i     = rho_g^(-1/2)
S_delta = min(actual leaf dx1, actual leaf dx2) / d_i
```

## Independent-review repairs

The rejected predecessor's seven findings are repaired as follows.

1. **Actual complete time inventory.** History accepts only at least 1000
   unique immutable per-cycle canonical telemetry artifacts opened and checked
   against their byte counts and SHA-256 digests. Every artifact binds the
   exact attempt, executable, normalization, trusted Q011 source commit, and
   sealed escape implementation identity. Inline invented records, a summary
   interval, an arbitrary source, or a one-cycle assertion cannot pass.
2. **Evidence-bound normalization and snapshots.** Normalization is read from a
   bound canonical runtime record and cross-bound to decoded runtime input
   parameters, deck, canonical source manifest, source archive, executable, and
   exact trusted source commit. Each snapshot manifest binds all seven raw
   products. The reducer opens and strictly decodes mesh/current products from
   actual Athena binary bytes and particles from actual AthenaK binary PVTK
   bytes, then proves equality with the supplied decoded objects and manifest
   digests. Arbitrary ASCII raw products fail. History reopens and reparses
   every retained raw, manifest, normalization, deck, source, and executable
   artifact before accepting a snapshot. History rejects stripped dictionaries
   and requires complete `ApplicabilitySnapshot` objects with retained evidence
   roots, immutable digest-bound maps, all regions, particle populations,
   spectrum, and provenance.
3. **Escape applicability and closure.** The runtime ledger closes exact
   particle counts and bounded macro-mass, kinetic-energy, and vector-momentum
   face residuals. Escaped particles must be included in all-cycle exposure,
   high-energy, maximum-energy, and gyroradius extrema. Escaped count,
   macro-mass, or kinetic-energy fraction above the claim-specific
   preregistered `1e-3` bounds rejects acceleration claims.
4. **Finite scattering-amplitude floor.** `Q011-APP-DI` requires precursor
   `delta B_rms/B0 >= 0.1`; nonzero machine-scale fluctuation power is
   insufficient.
5. **Claim scope.** `Q011-APP-LAMBDA`, `Q011-APP-DI`, or incomplete
   time/escape evidence rejects `Emax`, high-energy slope/cutoff,
   acceleration-rate, and acceleration-efficiency claims.
6. **Strong high-energy containment.** Active-plus-escaped global maximum and
   macro-q990-specific-energy high-tail maximum must each satisfy
   `r_g/Ly <= 1/8`. The original macro/energy `q999 <= Ly/8` and energy
   fraction above `Ly/4 <= 1e-3` checks remain explicit, but cannot dilute a
   failing high-energy or global maximum.
7. **Actual-leaf-aware magnetic estimator.** The reducer rejects a
   finest-composite FFT interpretation. It restricts to the coarsest source
   level present in the precursor, computes resolved power there, and
   conservatively assigns all subgrid residual power to the high-wavenumber
   upper bound. A mixed-level `[0,1]` adversarial state exercises this path.

The latest independent rejection of `72792e172` is additionally repaired:

1. every admitted cycle record is decoded only from a separately opened bound
   artifact; 1155 inline invented records and arbitrary source identities are
   explicitly non-admitting;
2. raw snapshot products are decoded from retained Athena binary/PVTK bytes,
   rather than accepted from caller-fabricable decoded digests;
3. each nominal restart checkpoint must be the first committed cycle whose end
   crosses the slot and must bind its previous committed cycle and time;
4. reason-coded physical escape is accepted only on `outer_x1`; any nonzero
   `inner_x1` escape rejects;
5. high-energy slope/cutoff applicability requires a separately bound complete
   active-plus-escaped binwise/tail reduction whose bins cover the all-cycle
   maximum energy; absent or biased evidence rejects that claim;
6. runtime coverage starts at the first cycle whose start is at or above
   `t=45` and whose bound previous committed time is below `t=45`; exact
   equality of the first start to `t=45` is not required.

## Snapshot evidence

Every retained snapshot includes:

- exact canonical normalization/source/deck/executable bindings;
- all retained raw-product bindings and decoded-payload digests, reopened and
  reparsed during history reduction;
- immutable digest-bound maps of `R`, `Lambda`, `d_i`, `S_delta`, actual leaf
  spacings, and gas-frame current magnitude;
- full-domain and detected-front-centered downstream, precursor, and
  far-upstream area/current-weighted statistics;
- a shock-transition-excluded actual-leaf ion-scale result and full precursor
  magnetic estimator record;
- macro/energy-weighted active, upstream, and high-energy-tail particle
  exposure statistics;
- local-B TSC particle gyroradius statistics and strict global/high-energy
  maxima.

Per-cell maxima are mandatory. No rare-cell waiver exists. The shock transition
within `120 c/omega_pi` of the detected front is excluded only from the
`S_delta` test. Microscopic shock-structure and self-consistent-injection claims
remain permanently excluded.

## Sealed escape-accounting contract

The successor binds the sealed physical-boundary escape implementation at
`d614e5c84aad3a541dc96af68ef1178dabc66f71`.

At every production restart checkpoint with nominal slot
`t=100,200,...,1200`, the reducer requires:

- the first committed cycle crossing that nominal slot, with exact previous
  committed cycle/time chronology bound to the all-cycle telemetry;
- one digest-bound restart payload and one digest-bound full `prtcl_all` PVTK
  payload with the exact observed committed cycle and time;
- complete schema-1 `ps_escape` and complete schema-3 CR source ledgers from
  the actual restart header;
- `ps_escape_audit_calls == 2 * cycle` for paper VL2 and
  `ps_escape_last_audit_time == observed_committed_time`;
- completed startup-cohort removal and
  `ps_escaped_initial_cr_count_global == 0`;
- exact injected, startup-removed, escaped-injected, and active count closure,
  plus mass closure at relative residual `<=1e-12`;
- active injected count and mass reproduced from the bound PVTK payload;
- active kinetic energy inside the exact float32 PVTK reconstruction interval;
- cumulative escaped maximum energy and gyroradius equal to the bound per-cycle
  inventory through the checkpoint cycle.

The terminal sealed ledger is also cross-bound to the independent runtime
boundary ledger for escaped count, mass, kinetic energy, and all three momentum
components. The boundary ledger requires exact reason codes and permits
physical escape only through `outer_x1`; any nonzero `inner_x1` escape rejects.
A positive escaped census must carry positive mass, energy, and global
gyroradius evidence. Counts, masses, energies, audit counts, audit times, and
escaped maxima must be monotonic. Signed momentum components are quantified and
residual-bounded but are not required to be monotonic because physically valid
cumulative signed momentum may change direction.

The conservative kinetic-energy escape fraction uses the lower bound of active
energy reconstructed from the terminal float32 PVTK payload. Thus active
particles cannot disappear from snapshots and make an acceleration claim pass.
Absent, incomplete, forged, unbound, cadence-incomplete, or internally
inconsistent restart/PVTK/telemetry evidence fails closed.

## AthenaK-selected gates

The cited literature gives qualitative applicability conditions, not these
numeric cutoffs. Every numeric threshold below is explicitly
**AthenaK-selected**.

| Gate | Acceptance |
| --- | --- |
| `Q011-APP-NORM` | exact bound normalization, source, executable, raw and decoded snapshot provenance |
| `Q011-APP-R` | full-domain per-cell `Rmax <= 0.01` |
| `Q011-APP-LAMBDA` | full-domain per-cell `Lambdamax <= 0.1` |
| `Q011-APP-DI` | outside shock transition `S_delta,min >= 1`; precursor `lambda_B,char/d_i,max >= 10`; sub-`10d_i` power upper bound `<= 0.05`; `delta B_rms/B0 >= 0.1` |
| `Q011-APP-RG` | snapshot macro/energy `q999/Ly <= 1/8` and energy fraction above `Ly/4 <= 1e-3`; all-cycle active-plus-escaped global and high-energy-tail `r_g,max/Ly <= 1/8`; claim-specific escaped count, mass, and kinetic-energy fractions each `<= 1e-3` |
| `Q011-APP-TIME` | separately opened byte-bound contiguous per-cycle telemetry from the first start crossing `t=45` through exact `t=1200`; first-crossing checkpoint cadence; complete active/boundary/escaped exposure; outer-`x1`-only escape; closed count/mass/energy/momentum ledgers |

## Claim-specific escape rules

All four acceleration claims use separately named, preregistered count,
macro-mass, and kinetic-energy fraction bounds of `1e-3`:

- `Emax_claim`: active-plus-escaped all-cycle maximum and fraction bounds;
- `high_energy_slope_or_cutoff_claim`: separately bound complete active plus
  escaped energy-bin and high-energy-tail evidence, with each high-energy-bin
  and aggregate-tail escaped count/mass/energy fraction at most `1e-3`; absent,
  incomplete, non-covering, or biased evidence rejects the claim;
- `acceleration_rate_claim`: active-plus-escaped all-cycle maximum and fraction
  bounds;
- `acceleration_efficiency_claim`: escaped energy explicitly bounded, with
  count and mass bounds retained conservatively.

Incomplete exposure or escape accounting rejects all four claims. No claim may
pass because energetic particles disappeared from the active population.

## Literature relationship and no-Hall boundary

Bai et al. (2015, ApJ 809, 55; arXiv:1412.1087) requires dilute CR loading,
identifies the CR-Hall parameter, requires `Lambda << 1` for no-Hall behavior,
and limits MHD-PIC interpretation to scales much larger than `d_i`. It does not
prescribe `0.01`, `0.1`, a decade, `5%`, or `delta B/B0=0.1`.

Sun & Bai (2023, MNRAS; arXiv:2304.10568), Section 5.4, states that the
transverse domain contains several high-energy CR gyroradii. It does not
prescribe a numeric containment cutoff. Requiring both the global and
high-energy-tail maxima below `Ly/8` is the controlling AthenaK-selected
interpretation: the transverse box contains at least eight of the largest
qualifying gyroradii. The q999 and energy-fraction checks are retained as
additional preregistered diagnostics.

This successor does not add Hall physics. `Lambda > 0.1` rejects the no-Hall
physical interpretation and related acceleration claims; `Lambda >= 1`
additionally rejects any claim that the no-Hall run approximates the target
plasma.

## Residual production obligations

No conforming production runtime artifact exists in this worktree. Tests use
synthetic adversarial evidence only. Remaining production obligations are:

1. implement and independently review the corresponding C++ separately
   retained per-cycle telemetry, particle-exposure, first-crossing checkpoint,
   reason-coded outer-`x1` escape, and slope/cutoff binwise-tail reductions;
2. execute the registered production campaign without changing thresholds
   after observing results;
3. retain complete bound raw, decoded, map, runtime, restart, PVTK, and escape
   artifacts;
4. independently review the resulting physical interpretation and claim scope.
