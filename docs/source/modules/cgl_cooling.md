# History-controlled CGL cooling

The `cgl_lf_cooling` problem generator reuses `cgl_lf_paper` initialization and
diagnostics, adding a local cooling source that preserves
$\Delta p=p_\perp-p_\parallel$. It supports active CGL on a uniform periodic mesh,
with or without Landau-fluid conduction and collisions. Existing problem
generators are unchanged.

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

At scheduled history samples it reuses the existing local `therm_cgl` and volume
sums, then performs one two-scalar MPI all-reduction. It does not read the history
file or traverse the mesh again. The usual history reduction remains unchanged.
The coefficient held until the next scheduled sample is

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
 \max(\min(p_\perp,p_\parallel)-p_{\rm floor},0)
\right],\qquad p_{\rm iso}=\frac{2p_\perp+p_\parallel}{3}.
$$

Total energy decreases by $3\delta p/2$, and the stored CGL anisotropy is updated
consistently. Density, momentum, and magnetic field are unchanged by cooling.
The common cap protects both pressure floors. If encoding in total energy loses
the retained heat to roundoff, the existing EOS floor repair supplies a
representable state. This can change $\Delta p$ at roundoff level; a repair that
would add net energy is discarded instead. Floating-point storage limits still
apply when thermal energy is tiny compared with kinetic or magnetic energy.

The source runs once before the existing collision task, after the hyperbolic
update. For split LF, the post-LF collision update remains in its existing place.
This is exact for isolated cooling with a held coefficient and first-order split
coupling to the fluid evolution. A finite exponential update at every RK source
stage would instead corrupt the RK combination of the conserved anisotropy.
Uniform periodic ghost cells receive the same source as owned cells, avoiding
another boundary exchange. AMR, passive CGL, and frozen backgrounds are rejected.

The held coefficient and nominal history schedule are saved at full precision in
the normal restart parameter header. Extra finalization history rows do not
update the controller. On restart, the pgen restores the scheduled history clock,
including when the run stopped between samples. Keep the cooling parameters and
history cadence unchanged for identical continuation. An older checkpoint without
controller state starts with cooling off until the next history sample.

The implementation stays in `src/pgen/tests/cgl_lf_cooling.cpp`, with only normal
pgen registration elsewhere. A callback after the generic history reduction would
avoid the small extra collective but require changing the output framework;
reading the `.hst` file would add I/O and restart dependencies. Neither is needed.

Run the focused CPU/MPI tests with a double-precision `built_in_pgens` executable:

```sh
ATHENAK_TEST_EXE=/absolute/path/to/athena PYTHONPATH=tst \
  python3 -m pytest -q tst/test_suite/cgl/test_cgl_lf_cooling.py
```
