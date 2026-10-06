# Task 1: face-stencil LF timestep and independent STS safety

Triage: **adapted and implemented**. The requested cell-centered estimate was
insufficient for the staggered stencil. The implementation bounds the absolute
rows of the actual local, frozen temperature Jacobian, including unequal cell
widths, staggered normal B, neighboring density, both temperatures, exact VL4
limiter derivatives and grad-B reverse coupling. `time/sts_safety` defaults to
0.9, must be finite in (0,1], and replaces advective CFL scaling in both STS
budget selection paths. Explicit processes retain their CFL factor. RKL2
coefficients, odd-stage selection and collision scheduling are unchanged.

The [derivation](../../../source/modules/cgl_lf_timestep.md) gives the formulas
and assumptions. This is a frozen uniform-grid Jacobian magnitude bound, not a
proof of nonlinear RKL2 stability or a bound on the composite AMR operator.
Conductivity/cap derivatives, discontinuous scattering switches and interleaved
projections require separate validation. The work order's proposed 45-degree
coefficient of five does not bound the implemented VL4 stencil: its smooth
nonzero-slope row is six, and its zero-slope kink envelope is twenty, in units
of chi/h^2. A bounded secant coefficient is not a limiter derivative.

## Independent mathematical and runtime evidence

The committed JSON records include 41 matrix/envelope cases, 13 grad-B
finite-difference Jacobians (maximum relative row error 8.02e-9), production face
formulas compared with independent stencils (approximately 5e-16 agreement),
and compiled exceptional-arithmetic checks against 80-digit Decimal references.
The latter cover seven coefficient/geometry/density cases, six limiter cases
and four local sound-speed cases. Cases include directional fields, random
fields, narrow reversals, density jumps, unequal spacing and 3D geometry.

A divergence-free checkerboard has unit normal face fields but cell-centered
B=(0,0,epsilon). Its true parallel stiffness grows as epsilon^-2. The new
2D/3D regression checks this exact row and catches the old estimate. Another
independent limiter case has radius 72.9081, exceeding a proposed secant bound
54.8585; the implemented derivative envelope is 99.8419. Its nonlinear RKL
trajectory stayed bounded, so it is a mathematical counterexample to that
proposed bound, not an observed production instability.

Both complete CPU/HIP builds passed. The HIP full workflow passed 29 cases,
and both active paper smoke cases passed. The focused original run on each
backend gave 33 passes and six failures: five undeclared input overrides and
one stale heated-step expectation. After correcting only those test fixtures,
all 12 affected/adjacent cases passed on each backend. Original failures remain
in scratch. No growth, accuracy, admissibility or repair tolerance was relaxed.
Wider final integration acceptance is recorded separately in the WO2 report.

Changed expected values follow independent formulas:

- Density contact C=200: initial dt is 1.704732834526797e-4, from
  `40*0.9/[64^2*chi*(C+3)]`; seven stages per half-sweep remain expected.
- The 1D/2D reversal reference budgets become 2.4229069212191616 and
  0.006978246033076956. Their inputs retain 20 cycles and use tlim=1000 so the
  larger budget cannot end the 1D regression early.
- The uniformly heated step needs 13 post stages plus seven pre stages,
  replacing the stale expectation of 11 post stages. The test independently
  computes both the initial budget and the heated pressure/stage requirement.

## GPU performance and result changes

All measurements used MI250X/gfx90a, CCE20 with the verified Task0 compiler
contract, Kokkos 08ceff92, double precision and Cray MPICH9.0.1 in allocation
5628672. Each configuration used one warmup and three unprofiled repeats in
an exclusive application window. Wave/SMR repeats alternated reference and
candidate. Input/output schedules and rank layouts were held fixed. The
reference is the Task0 physical-corner correction; periodic cases are
unchanged by that correction. Separate synchronized profiles are diagnostic
only and excluded from the timing medians.

The printed solver timer starts after initialization and initial outputs. It
includes the evolution loop, in-run/final outputs and final problem diagnostics;
it excludes process launch and initial setup. Launch-inclusive wall times are
retained separately. These are application timings, not isolated kernel costs.

| Input | Ranks | Reference/candidate s per cycle | Mean stages per half-sweep | Measured improvement |
| --- | ---: | --- | --- | --- |
| Original paper turbulence, 64 cycles | 1 | 0.183986 / 0.082471 | 9 / 3 | 2.23x per cycle |
| Supplemental 3D turbulence, 64 cycles | 1 | 0.184772 / 0.103999 | 9 / 4.1875 | 1.78x per cycle |
| Propagating slow eigenwave, t=0.05 | 1 | 0.026953 / 0.011358 | 12.9718 / 3 | 2.37x per cycle and common time |
| Sinusoidal parallel decay | 1 | 0.014841 / 0.019234 | 5 / 4.9785 | 1.73x to common time; slower per cycle |
| SMR decay | 1 | 0.052342 / 0.055759 | 6.9429 / 7 | 1.26x to common time; slower per cycle |
| SMR decay | 4 | 0.044181 / 0.047911 | 6.9429 / 7 | 1.24x to common time; slower per cycle |

The decay run takes 208 versus 93 cycles; SMR takes 70 versus 52. Their faster
completion is due to fewer cycles, despite a costlier bound per cycle. All
three have identical final float32 field outputs. An accidental four-rank
launch of the one-block 1D input was rejected as expected: eight failed launches
are retained and explicitly excluded; they are not performance data.

The propagating slow wave retains the original 256 cells, amplitude 1e-5 and
strict 1e-3 analytic tolerances, extending only tlim to 0.05. All eight runs
passed and repeat exactly within each binary. Both controllers take 284 cycles;
median total solver time is 7.65465 versus 3.22578 seconds. The final maximum
velocity output difference is 2.73e-11; density, pressures and B agree at output
precision. This is the propagating-wave benchmark; the sinusoidal decay fixture
is reported separately.

### Turbulence qualification

The original paper input is exactly planar with initial B along z. Its generic
forcing policy and default non-isotropic spectrum have no z phase. Full Real
restart inspection confirms zero transverse B, zero z variation and exactly
zero applied LF q-work. Its 2.23x result measures LF overhead reduction and flow
throughput; it does not validate nonzero nonlinear heat transport. Both original
full-duration runs reach t=1 in 620 cycles with zero repairs and mass drift
1.11e-16. Maximum final density/E/A differences are 5.42e-14/8.86e-13/5.48e-13.

A separately labeled supplemental input changes only the forcing projection to
`mks24_alfvenic_perpendicular`. It keeps the original 48x48x96 grid, 24x24x48
blocks, seed, forcing strength, correlation time, LF/limiter settings and tlim=1.
The original deck is preserved. The original strict-admissibility default is
false; these runs preserve it and explicitly require all repair counters zero.
All twelve supplemental applications pass: warmup, three 64-cycle repeats,
separate profile and full-duration run for each executable.

Both full 3D runs reach t=1 in 619 cycles, with zero density/pressure floor,
nonfinite, nonpositive, hard-bound and wall counters and mass drift 1.11e-16.
Actual parallel/perpendicular q-work magnitudes are approximately 1.38e-3 and
4.59e-4, far above the activity threshold 3.41e-12. Transverse B reaches 0.3114
and within-block z variation 0.3037. The final maximum Real density/E/A
changes are 4.7573e-7/7.3080e-6/9.0803e-6; relative L1 changes are
4.2103e-8/5.8129e-8/2.5934e-5. Full-run solver times are 107.1858 versus
65.37339 seconds, single samples distinct from the three-repeat timing result.
This supports finite-time nonlinear activity/stability and controller agreement;
it is not a turbulence convergence campaign or a closed forced-energy ledger.

## Provenance and restart comparison

All original commands, binaries, inputs, logs and outputs are under
`/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2`. The
[artifact index](artifact-index.json) hashes the copied mathematical/comparison
records and identifies their original paths. The integration candidate source/docs
manifest is [task1-documented-candidate.json](task1-documented-candidate.json).
One documentation sentence was corrected after independent review: nonfinite
VL4 slopes fail closed; they do not use the finite kink envelope. All eleven
tested source/input/test files remain byte-identical. The CPU/HIP source LF hash is
`1631a9e4aa308a019249179594e26b2073365aea161c10e1c4de6646a2a1331d`.

The full CPU/HIP executables are `build/task1-cpu-refill/src/athena`
(SHA256 `a0ef5a65c4d4f1f3962d7ee1823b8aaee3b150d38c58ccf2d926ac12473ca93a`)
and `build/task1-hip-refill/src/athena`
(`5f7fdae897632d6b56cbef37220791394cbbbb6c3256dfe83ed1e082f8e74efe`).
The custom paper HIP executable is
`task1-research/paper-build/bin/athena-cgl-lf-paper-hip`
(`5ee92eacc28f2b52f900cd068e2e6ef1f694c344482362a785f8443579c81c7f`).

The [restart comparator](../restart_compare.py) defaults to excluding only the
36 unused root coarse-index bytes diagnosed in Task0. Turbulence comparisons
explicitly opt into two further layout-validated exclusions: dormant startup
RNG storage, and four bytes of native RNG struct alignment padding. The latter
occurs between `iset` and `gset`, including in evolved records. The
[compiled ABI probe](rng_restart_abi.cpp) establishes offsets and sizes. All live
RNG members, seed/cache flags, diagnostics and physical fields remain compared.
Raw hashes and exact excluded ranges are retained; outputs are never rewritten.
Initial raw repeatability failures remain archived beside these narrow audits.
