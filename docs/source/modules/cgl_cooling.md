# History-controlled CGL cooling

The `cgl_lf_cooling` problem generator reuses `cgl_lf_paper` initialization and
diagnostics, adding a local cooling source that preserves
$\Delta p=p_\perp-p_\parallel$. It supports active CGL and passive CGL on an
adiabatic MHD background, on a uniform periodic mesh with or without
Landau-fluid conduction and collisions.

Build the production cooling executable with:

```sh
cmake -S . -B build-cgl-cooling -DPROBLEM=cgl_lf_cooling
cmake --build build-cgl-cooling --parallel
```

CMake resolves this name to `src/pgen/cgl_lf_cooling.cpp`. This custom build
always selects cooling, independent of `problem/pgen_name`. Keep the following
name in input files so they also work with a `built_in_pgens` executable.

In an existing `cgl_lf_paper` input, change/add:

```ini
<problem>
pgen_name = cgl_lf_cooling
cooling_piso = 5.0  # target volume-averaged isotropic pressure, in code units
cooling_tau = 0.2   # controller response time, in code time units

<output1>
file_type = hst
dt = 0.02
```

These numbers illustrate a target matching the existing beta-10, unit-field
input; choose the target and response time for the actual experiment. The pgen
enables user history automatically. There must be exactly one history output,
with a positive cadence. A positive `dcycle` can replace `dt`. Do not enable `problem/user_srcs`:
this pgen enrolls its own split user source, outside the RK stage-source hook.

At scheduled history samples it reuses the local thermal-energy and volume
sums, then performs one two-scalar MPI all-reduction. It does not read the history
file or traverse the mesh again. The usual history reduction remains unchanged.
For active CGL, the coefficient held until the next scheduled sample is

$$
\alpha=\frac{1}{\tau}\max\left(1-\frac{(3/2)p_{\rm target}V}{E_{\rm th}},0\right),
\qquad E_{\rm th}=\int(p_\perp+p_\parallel/2)\,dV.
$$

The user history adds `cool_alpha`, the coefficient for subsequent steps.
Cooling stops below the target. At approximately steady heating, the mean energy
excess is roughly $\overline P_{\rm heat}\tau$. Choose the response time longer
than the interval between history samples, but short compared with the thermal
heating time. Time-based history requires `0 < dt <= cooling_tau`; for cycle-based
history, check the resulting physical interval yourself.

Each full step removes the same pressure decrement from both components:

$$
\delta p=\min\left[
 p_{\rm iso}(1-e^{-\alpha\Delta t}),
 \max(\min(p_\perp,p_\parallel)-p_{\rm stop},0)
\right],\qquad p_{\rm iso}=\frac{2p_\perp+p_\parallel}{3}.
$$

Here $p_{\rm stop}=\max(p_{\rm floor},|\Delta p|/9)$ stops cooling before the
larger pressure exceeds ten times the smaller one. This fixed numerical guard
preserves $\Delta p$ and leaves the mirror/firehose limiters unchanged. Cells
already beyond 10:1 receive no cooling; the guard neither repairs them nor
prevents transport or fluid evolution from creating them. Restart before the
extreme ratios develop to test this prevention. The cutoff can reduce the
cooling power and should be revisited if many cells reach it.

User history reports `cool_hits`, the number of owned cell-updates whose cooling
was reduced by the 10:1 guard since the previous history output, and `cool_frac`,
that count divided by all owned cell-updates over the interval. These include
complete suppression in already extreme cells and exclude pressure-floor-only
clipping and ghost cells. A cell can contribute on successive steps. Restarting
resets these diagnostics. Repeated output at the same cycle repeats the same
interval values. The counters use the cooling source pass and are read only at
history output. Use a fresh history file when adopting the additional columns.

Active CGL total energy decreases by $3\delta p/2$, and the stored anisotropy is updated
consistently. Density, momentum, and magnetic field are unchanged by cooling.
The common cap protects both pressure floors. If encoding in total energy loses
the retained heat to roundoff, the existing EOS floor repair supplies a
representable state. This can change $\Delta p$ at roundoff level; a repair that
would add net energy is discarded instead. Floating-point storage limits still
apply when thermal energy is tiny compared with kinetic or magnetic energy.

The source runs once before the existing collision task, after the hyperbolic
update. For split LF, the post-LF collision update remains in its existing place.
The exponential is exact for isolated, uncapped cooling with a held coefficient;
the source has first-order split
coupling to the fluid evolution. A finite exponential update at every RK source
stage would instead corrupt the RK combination of the conserved anisotropy.
Uniform periodic ghost cells receive the same source as owned cells, avoiding
another boundary exchange. AMR and isothermal passive CGL are rejected.

The same executable can load its own checkpoint with the
`inputs/particles/frozen_cgl.athinput` overlay for frozen-field CR tracing.
During a frozen restart, the saved controller and cooling history columns are
retained while cooling and controller updates stop. Use a new output directory
for the particle run. The ordinary fluid restart restores cooling normally.

The held coefficient and nominal history schedule are saved at full precision in
the normal restart parameter header. Extra finalization history rows do not
update the controller. On restart, the pgen restores the scheduled history clock,
including when the run stopped between samples. Keep the cooling parameters and
history cadence unchanged for identical continuation. An older checkpoint without
controller state starts with cooling off until the next history sample.

## Adiabatic passive comparison

To carry passive CGL pressures and LF transport on an adiabatic MHD flow, add:

```ini
<mhd>
eos = cgl
passive = true
passive_eos = adiabatic
gamma = 1.666666666666667

<problem>
pgen_name = cgl_lf_cooling
passive_delta = true
p_gas0 = 5.0
cooling_piso = 5.0
cooling_tau = 0.2
```

`p_gas0` sets the independent MHD gas pressure; its default is the initial
CGL isotropic pressure. Use $\gamma=5/3$ and the same initial pressure, mesh,
forcing, reconstruction, cooling target, response time, and history cadence as
the active run. The example
`inputs/cgl_lf_paper/cgl_lf_paper_cooling_passive_adiabatic.athinput`
is a small runnable case. Copy the production active run's other settings when
making its passive companion.

The controller uses the adiabatic gas's `therm_mhd` history sum,
$E_{\rm gas}=\int p_{\rm gas}/(\gamma-1)\,dV$, with target
$p_{\rm target}V/(\gamma-1)$. Thus its held coefficient is
$\alpha=\max(1-p_{\rm target}V/[(\gamma-1)E_{\rm gas}],0)/\tau$.
The gas internal energy receives the same exponential cooling law, limited by
its pressure, temperature, and entropy floors. The same coefficient cools the
passive CGL pair using the fixed-$\Delta p$ update above. Their 10:1 cap applies
only to that pair; it cannot reduce the gas cooling or feed passive
thermodynamics back into the flow.
The passive pressures receive no dissipated MHD heat and are not separately
thermostatted, so their mean need not remain at the gas pressure target.

The MHD total energy is stored separately from the passive J/A invariants.
LF, collisions, limiters, and passive pressure repairs do not change it.
User history reports `therm_cgl`, `therm_mhd`, and `p_gas` separately. Primitive
output adds `gas_eint`; the legacy `eint` field still means $p_\parallel$.
Matching controls does not imply identical turbulent trajectories: the active
pressure stress changes the flow, and different timesteps sample the forcing
process differently.

New adiabatic passive checkpoints contain the extra gas energy and use their
own passive encoding. Start this comparison from initial conditions; changing
`passive_eos` on an isothermal passive checkpoint is rejected.

The cooling implementation stays in `src/pgen/cgl_lf_cooling.cpp`. A
callback after the generic history reduction would avoid the small extra
collective but require changing the output framework;
reading the `.hst` file would add I/O and restart dependencies. Neither is needed.

Run the focused CPU/MPI tests with a double-precision `built_in_pgens` executable:

```sh
ATHENAK_TEST_EXE=/absolute/path/to/athena PYTHONPATH=tst \
  python3 -m pytest -q tst/test_suite/cgl/test_cgl_lf_cooling.py \
  tst/test_suite/cgl/test_cgl_adiabatic_passive.py
```
