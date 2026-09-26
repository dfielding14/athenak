# Current-limited Ohmic resistivity

This experimental Newtonian closure increases magnetic diffusivity when a current
sheet approaches an ion inertial length. It does not impose a reconnection rate
or resolve Hall/kinetic physics. The species convention is singly ionized ions
with fixed composition; magnetic pressure is $B^2/2$ and $\mathbf J=\nabla\times\mathbf B$.

```ini
<mhd>
resistivity_model = current_limited
ohmic_resistivity = 0.000001
eta_max = 0.005
d_i = 0.005
b_rec_method = constant
b_rec = 1.0
resistivity_integrator = explicit
```

`ohmic_resistivity` is the background diffusivity $\eta_0>0$, `eta_max` must be
at least $\eta_0$, and `d_i` is the positive ion length at unit density in code
units. `b_rec` is a positive, prescribed reconnecting field, independent of the
local magnetic-field magnitude or guide field. All parameters must be finite.
The default model is `constant`; its existing update is preserved. Setting
`eta_max = ohmic_resistivity` also uses those constant kernels.

With individually floored neighboring densities averaged to the edge,

$$
q=\frac{|\mathbf J|d_i}{B_{\rm rec}\sqrt{\rho}},\qquad
q_* = 1-\sqrt{\eta_0/\eta_{\max}},\qquad
\eta_* = \sqrt{\eta_0\eta_{\max}},
$$

$$
\eta(q)=\begin{cases}
\eta_0/(1-q),&q\le q_*,\\
\eta_{\max}-(q_*/q)(\eta_{\max}-\eta_*),&q>q_*.
\end{cases}
$$

This is the total diffusivity. The code neither adds another background term nor
clips the current. A sheet with upstream field $B_{\rm up}$ activates at thickness
approximately $d_i(\rho)B_{\rm up}/B_{\rm rec}$; it need not activate at $d_i(\rho)$.

Each native curl component lives on a different constrained-transport edge. The
two transverse components are averaged onto the edge before evaluating $q$.
The same edge coefficient multiplies the EMF and its contribution to the
conservative resistive Poynting flux. There is no separate Joule-heating source.
Two ghost cells suffice for this fixed-field stencil.

The coefficient is recomputed at every explicit or STS stage. The timestep bound
uses `eta_max` in the existing $\Delta x^2/(2N\eta)$ estimate, where $N$ is the
number of active dimensions. To use RKL2 also set `resistivity_integrator = sts`
and `<time>/sts_integrator = rkl2`. The differential-diffusivity bound motivates
this estimate; numerical tests are still required for the interpolated nonlinear
operator. Relativity and nonzero `eta_ad` are unsupported.

## Diagnostics

`mhd_eta` and `mhd_q` are cell-centered reconstructed diagnostics: average each
native current component to the cell center, then evaluate the closure using the
cell density. They are proxies, not averages of the nonlinear edge coefficients.
`mhd_eta` also works with constant resistivity; `mhd_q` requires `current_limited`.
Use the resistive test problem's edge maxima to monitor actual activation.
The cell-centered $\eta|\mathbf J|^2$ diagnostic is not an exact discrete heating
budget. Check total, magnetic, kinetic, and inferred internal energy separately.

The built-in `problem/pgen_name = resistive_tests` supports `test = forcefree`,
`gaussian`, `sheet`, and `harris`. Example decks are in
`inputs/mhd/resistive_{forcefree,gaussian,sheet,harris,gem}.athinput`. The force-free
wave accepts integer `wave_n1`, `wave_n2`, `wave_n3`; its analytic face averages
preserve discrete divergence for oblique wavevectors. Harris uses a discrete
curl of a vector potential, pressure balance, and complete total energy.

With `problem/user_hist = true`, the user history records actual-edge maxima
`q_max_edge`, `eta_max_ed`; cell-proxy volume fractions `frac_qstar`, `frac_q1`;
integrated `etaJ2_cc`; and its high-$q$ fraction `heat_frac`. Fractions use global
integrals and maxima use MPI maxima. Harris additionally records `x_q`, `x_etaJ`,
`x_etaJz`, `x_Ez` at the symmetry-centered X edge, with electric fields normalized
by $B_0v_{A0}$. `x_Ez` uses interpolated velocity and field plus the native-edge
resistive EMF; it is not the Riemann solver's numerical EMF. These diagnostics
assume the principal X point remains at the domain center.
Harris histories also record `ref_etaJz` and `ref_Ez` at the positive-x boundary
on the sheet midplane. The open single-sheet deck initially has no O point;
the boundary is a fixed flux reference, and its generally nonzero electric
field must be subtracted. For $\psi=A_z(x_{\rm ref},0)-A_z(0,0)$,
$d\psi/dt=E_z(0,0)-E_z(x_{\rm ref},0)$. Uniform sheet diffusion therefore gives
zero reference-subtracted rate even though the individual electric fields are
nonzero.

For an independent flux-based rate, supply verified X/O positions on the same
midplane to the postprocessor; it does not track changing topology:

```sh
python vis/python/reconnection_flux.py run/bin/*.state.*.bin \
  --xpoint 0 --opoint 12.8 --midplane 0 --b0 1 --rho-up 0.2 --di 1 > flux.csv
```

The example coordinates are for the GEM deck. Use full, unsharded, unrefined 2-D
`mhd_u_bcc` dumps. The signed flux is $A_z(O)-A_z(X)$; the script differentiates
it in time and normalizes by $B_0v_{A0}$. Energy integrals are per unit depth.
Output cadence, cell-centered quadrature, and point placement affect the rate.
The script also reports slope thickness $\delta=B_0/|\partial_yB_x|$ at X, its
ratio to local $d_i(\rho_X)$, and the connected half-maximum $|J_z|$ half-length
along the midplane. Their ratio is a geometric proxy for aspect ratio. If the
current layer reaches the boundary without falling below half maximum, its
length/aspect are unmeasured (`nan`); these proxies assume an axis-aligned sheet.

## Acceptance and scope

The fast numerical suite tests the scalar closure, unchanged constant runs,
the degenerate limit, weak Gaussian diffusion, nonlinear force-free decay,
sheet spreading, decomposition, refinement, and restart behavior. Run it with
the same CPU/MPI test environment as the STS suite; set `ATHENA_BASELINE` to a
pre-change executable to enable the baseline comparison.

Reconnection input decks are experimental benchmarks. Resolution, $d_i/L$,
background diffusivity, cap, and guide-field scans are needed before making
claims about a steady reconnection rate or its independence of system size.
Small GEM runs provide a qualitative comparison only. Fast regression success
does not substitute for these scientific acceptance studies.

For the $L=1$, $\rho_{\rm up}=B_0=1$ Harris deck, uniform-$\eta$ control at
$S=2\times10^5$ uses `resistivity_model=constant`, `ohmic_resistivity=5e-6`.
The ideal control uses `ohmic_resistivity=0`, `resistivity_integrator=explicit`,
and `time/sts_integrator=none`. For both controls, change output 5 from `mhd_q`
to `mhd_bmag` (and its ID to `bmag`); the user-history $q$ proxy remains available.
The fiducial sheet has width $0.02L$, whereas $d_i=0.005L$. A jump radius chosen
to measure its initial upstream field must span several sheet widths; $3d_i$
does not yet reach that upstream field. Keep this distinction explicit when
comparing the two methods during sheet formation.

## Experimental jump estimate

Set `b_rec_method = jump`, `b_rec_radius = R`, and positive `b_rec_floor` to use

$$
B_{\rm rec}(\mathbf x)=\max\left[B_{\rm floor},\max_{k\ {\rm active}}
\frac{|\mathbf B_c(\mathbf x+m_k\Delta x_k\mathbf e_k)
-\mathbf B_c(\mathbf x-m_k\Delta x_k\mathbf e_k)|}{2}\right],\quad
m_k=\lceil R/\Delta x_k\rceil.
$$

This measures a field jump, so a uniform guide field cancels. Adjacent cell
estimates are averaged to the operator edge. `b_rec` is used only by the constant
method. The magnetic halo must satisfy `nghost >= max(m_k)+1`, including the
finest permitted AMR level. Meshblocks must accommodate that halo; multilevel
meshes also require even `nghost` and enough coarse cells. A large physical
radius can make the communication and memory cost substantial.

The cache is computed after initialization/refinement boundary repair and before
each complete timestep. It stays fixed through the pre-STS, RK, and post-STS
updates; current and density still change at each stage. Restart reconstructs
the cache, rather than serializing it. `mhd_brec` reports the cache used by the
last completed cycle (freshly reconstructed at initialization or refinement),
not an estimate recomputed from the output magnetic field. The associated
`mhd_eta`/`mhd_q` use that same cache and the output state.

Freezing the estimate introduces a temporal lag: the jump model must not be
assumed second order even though the fixed-field integrator is second order.
Its estimator, timestep sensitivity, decomposition, restart, and cost require
separate validation. Agreement of short runs is not evidence for a converged
steady reconnection rate or acceptance for turbulent production simulations.

## Local validation

The MPI CPU build uses Kokkos bounds checking: 102 diffusion/estimator/diagnostic
tests and 97 broader regressions pass, with no skipped baseline checks.
Constant-model outputs and
timesteps match the saved pre-change STS executable bit for bit; equal-coefficient
current-limited runs also reproduce the constant path with both integrators.
The numerical checks cover oblique three-dimensional currents, nonlinear
force-free decay, high-beta sheet spreading through $q<1$, conservative energy,
mesh refinement, and restart. All five example decks have also been initialized
and evolved in small smoke runs; these are not the resolution studies described
above. GPU execution has not been tested locally.

Fixed-field force-free spatial orders are 2.00--2.03; the finer temporal pair
gives order 2.01. Relative total-energy errors are at most $3.20\times10^{-14}$.
The sheet test reaches peak $q=0.658$ with both explicit and STS integration.

For the jump-model temporal test, halving the STS timestep cap gives field errors
$6.9997\times10^{-4}$, $3.4736\times10^{-4}$, and $1.7151\times10^{-4}$ against a
small-step reference: observed orders 1.01 and 1.02. In short Harris comparisons,
the maximum $B_x$ difference from fixed $B_{\rm rec}$ is $3.65\times10^{-5}$
without a guide field and $3.69\times10^{-5}$ with $B_g=B_0$. These are early-time
field comparisons, not measurements of a converged steady reconnection rate.

Small local cost measurements use one MPI rank, Debug/O1 with bounds checks,
and medians of three runs including process startup. Increasing ghost width
from 2 to 16 for a single $32^3$ block and 30 constant-diffusivity cycles changes
wall time from 2.75 to 4.62 s (1.68 times) and cell-array volume by 5.62 times,
with identical physical states. This measures local halo work, not inter-node
communication. On a $16^3$ block with 20 cycles and the same diffusion timestep
bound, nonlinear fixed-field STS costs 1.95 s versus 0.73 s for constant
diffusivity (2.66 times); both execute 40 sweeps and 200 stages. These timings
measure overhead, not a production speedup or equivalent physical evolution.

For a configured executable, run the diffusion tests from the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="$PWD/tst:$PWD/vis/python" \
  ATHENA=/absolute/path/to/athena \
  ATHENA_BASELINE=/absolute/path/to/prechange/athena \
  python -m pytest -q tst/test_suite/diffusion
```

An MPI launcher, NumPy, h5py, pytest, and a C++ compiler are required. The
`ATHENA_BASELINE` comparison is skipped if that executable is not supplied.
Final run logs include STS sweep/stage counts in addition to cycles and runtime;
these work counters restart from zero when a checkpoint is resumed.
