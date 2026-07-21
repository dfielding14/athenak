# Strong guide-field MHD ringing reproducer

## Summary

This reproducer evolves continuously driven, nearly isothermal ideal MHD in a
periodic cube with a strong magnetic guide field aligned exactly with the grid.
The Riemann solver and time integrator are held fixed at HLLD and RK2.  The only
comparison variable is reconstruction: PLM, PPMX, or WENO-Z.

At 256 cubed resolution, all three runs develop sharp, nearly grid-aligned
magnetic structures.  PPMX and especially WENO-Z develop a coherent
cell-to-cell alternating wake beside some of these structures.  PLM is much
smoother but more diffusive and is not a satisfactory replacement for the
higher-order methods.  The artifact is most obvious in magnetic-field
curvature and current, but it is directly visible in the unsmoothed magnetic
magnitude.

This is a reconstruction-sensitive HLLD/RK2 problem, not a comparison of time
integrators or Riemann solvers.

## Code and problem generator

- Branch: <https://github.com/dfielding14/athenak/tree/ringing-fix>
- Problem generator:
  <https://github.com/dfielding14/athenak/blob/ringing-fix/src/pgen/turb.cpp>
- Canonical input:
  <https://github.com/dfielding14/athenak/blob/ringing-fix/inputs/mhd/turb_ringing_reproducer.athinput>
- Frontier build helper:
  <https://github.com/dfielding14/athenak/blob/ringing-fix/configure_build_turb.sh>

The executable must be configured with `-DPROBLEM=turb`.  In this pgen,
`ifield=2` selects a uniform face-centered B3 field.  With the ideal-gas setup
below, `beta=0.04` gives an initial guide-field magnitude of approximately
7.071, five times the corresponding beta-one field.

## Fixed numerical and physical setup

| Item | Value |
| --- | --- |
| Domain | Periodic `[-0.5, 0.5]` cubed |
| Resolution | 256 cubed; cubic cells with delta-x = 1/256 |
| MeshBlocks | 64 cubed cells; 4 cubed blocks total |
| Initial state | density 1, zero velocity, uniform grid-aligned B3 |
| Equation of state | Ideal gas, gamma = 1.00001 |
| Plasma beta | 0.04 |
| Riemann solver | HLLD |
| Time integrator | RK2 |
| CFL number | 0.33 |
| Driving | Solenoidal OU forcing, seed 1 |
| Driving power | `dedt=0.25`, `tcorr=1.0` |
| Driven modes | k = 1 through 3, parabolic spectrum peaked at k = 2 |
| End time | 1.5 |

The input writes every 0.25 time unit, so output 00006 is the t=1.5 state.
Three ghost cells are retained for an identical mesh across all reconstruction
choices.

## Build

On Frontier, from the repository root:

```bash
git submodule update --init --recursive
./configure_build_turb.sh
```

This produces
`build-frontier-hip-cpe25.09-cce20-rocm6.4.2-turb/src/athena` using a Release
MPI+HIP build for MI250X GPUs.  On another platform, the essential CMake
selection is:

```bash
cmake -S . -B build-turb \
  -DCMAKE_BUILD_TYPE=Release \
  -DPROBLEM=turb \
  -DAthena_ENABLE_MPI=ON
cmake --build build-turb --parallel
```

Add the platform-appropriate Kokkos backend and architecture options.

## Run exactly three variants

Keep every parameter fixed and override only the basename and reconstruction
method.  For example, with 64 MPI ranks:

```bash
mkdir -p run/plm run/ppmx run/wenoz

mpiexec -n 64 build-turb/src/athena \
  -i inputs/mhd/turb_ringing_reproducer.athinput -d run/plm \
  job/basename=turb_ringing_plm mhd/reconstruct=plm

mpiexec -n 64 build-turb/src/athena \
  -i inputs/mhd/turb_ringing_reproducer.athinput -d run/ppmx \
  job/basename=turb_ringing_ppmx mhd/reconstruct=ppmx

mpiexec -n 64 build-turb/src/athena \
  -i inputs/mhd/turb_ringing_reproducer.athinput -d run/wenoz \
  job/basename=turb_ringing_wenoz mhd/reconstruct=wenoz
```

The original Frontier realization used eight nodes, one MPI rank per GPU, and
64 ranks total.  Each of its three final dumps contained all 64 rank shards.

## What to inspect

Read the binary shards with
`vis/python/bin_convert.py::read_binary_as_athdf` or the equivalent function in
`bin_convert_new.py`.  Display exact cell values: do not interpolate, resample,
filter, or smooth the data.

At t=1.5, two compact discriminator regions are:

1. The x1 cell-centered plane at `x1=0.001953125`, around `x2=-0.1` and
   `x3=0.3`.
2. The x2 cell-centered plane at `x2=0.017578125`, around `x1=-0.1` and
   `x3=0.3`.  This is the cell selected by the requested input slice coordinate
   `x2=0.01953125`.

The characteristic symptom is a narrow, grid-aligned sequence of alternating
cell values immediately adjacent to or trailing a sharp magnetic structure.
It is coherent across multiple transverse rows rather than a single bad cell.
Curvature makes the alternating pattern visually strongest; current magnitude
(the square root of the `mhd_j2` output) also highlights it, and magnetic
magnitude confirms that it exists in the evolved field itself.

The observed amplitude ordering is:

```text
WENO-Z ringing > PPMX ringing > PLM ringing.
```

For the first region in the recorded t=1.5 realization, exact-cell magnetic
magnitude diagnostics were:

| Reconstruction | Second-difference RMS | RMS above 0.5 Nyquist |
| --- | ---: | ---: |
| PLM | 0.00271 | 0.000206 |
| PPMX | 0.00494 | 0.000520 |
| WENO-Z | 0.00589 | 0.000578 |

These values document the known realization; they are not proposed as strict
regression tolerances.  The qualitative cell-alternating morphology and the
ordering across the three reconstructions are the intended reproduction
criteria.

The three-way baseline was produced before the independent HLLD correction
described below.  Repeating the high-order controls with corrected HLLD left
the visible bands and ringing metrics nearly unchanged, so exact cycle counts
or cell values should not be expected to match bit for bit.

## Separate HLLD correctness problem found during the audit

This is a real solver defect, but it is not the cause of the reconstruction
ordering above and fixing it does not cure the turbulence ringing.

The ideal-gas HLLD implementation contained a legacy weak-normal-field
shortcut:

```text
0.5 Bn^2 < 1e-4 pt_star  =>  U**L = U*L and U**R = U*R.
```

That shortcut is invalid for a small but finite face-normal field.  The two
double-star states carry the rotational/Alfven-wave jumps.  Collapsing them to
the single-star states zeros the associated Rankine-Hugoniot flux corrections,
allowing the solver to return an F* flux in a wave-fan sector that requires
F**.  A weak normal field makes the rotational speeds approach the contact
speed; it does not make transverse velocity and magnetic-field jumps vanish.

This condition is particularly easy to satisfy when a strong transverse guide
field dominates the total star pressure while the turbulent face-normal
component is weak.  In this reproducer, `pt_star` is approximately 26, so the
old shortcut applied for roughly `abs(Bn) < 0.072`.

The branch ports the corresponding
[Athena++ correction](https://github.com/PrincetonUniversity/athena/commit/0e3ca3dde5a6faffcfd72d8e41e5df0db1ce5111),
also tracked as [AthenaK issue
#601](https://github.com/IAS-Astrophysics/athenak/issues/601).  The corrected
solver always constructs the Miyoshi-Kusano double-star states for finite Bn;
the separate denominator guards remain intact.

The focused regression is a 16-cell moving rotational discontinuity with
RK2, donor-cell reconstruction, HLLD, `Bn=0.01`, and opposing transverse
velocity and field.  Recorded behavior is:

| Solver | Maximum transverse velocity | L-infinity errors |
| --- | ---: | --- |
| Old weak-Bn shortcut | 1.0006068 | approximately 6.1e-4 |
| Corrected double-star states | 1.0 | velocity 3.9e-6; field 1.9e-6 |

The regression requires both errors to remain below 2e-5 and rejects velocity
overshoot.  The correction affects Newtonian ideal-EOS HLLD with finite weak
face-normal field and a nonsmooth transverse jump.  The isothermal HLLD path
is separate.
