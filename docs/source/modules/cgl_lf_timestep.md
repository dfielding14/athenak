# LF Discrete Timestep Bound

The LF timestep calculation follows the face stencil used by the heat flux.
It bounds the absolute row sums of a local temperature Jacobian, with the
coefficients specified below held fixed. The resulting parabolic reference
step is `2 / Lambda`. The RKL2 controller multiplies it by
`time/sts_safety` (default `0.9`, finite and in `(0, 1]`). Advective and explicit
parabolic limits continue to use `time/cfl_number`; see
[Super Time Stepping](super_time_stepping.md#timestep-control).

## Face Geometry And Temperature Variables

During an LF sweep, density and magnetic field are fixed. Write
$T_\parallel=p_\parallel/\rho$ and $T_\perp=p_\perp/\rho$.
For a face $f$ between cells $l,r$, let

$$
\overline B_f=\frac{|B_l|+|B_r|}{2},\qquad
\rho_f=\frac{\rho_l+\rho_r}{2},\qquad
\beta_{if}=\frac{|B_i|}{\overline B_f}.
$$

The normal field component comes from the staggered face array; transverse
components are the averages of neighboring cell-centered fields. Define
$b_{nf}=B_{nf}/\overline B_f$ and the corresponding transverse $b_{tf}$.
The quantity $b_{nf}$ can exceed one in magnitude. In particular, face fields
can oscillate while their cell-centered averages cancel. A cell-centered
unit-direction estimate does not bound that stencil.

Let $h_n$ and $h_t$ denote the receiving block's directional cell widths.
Faces rejected by the closure's weak-field admission gate contribute zero.
The same density and pressure floors as the face closure are used. Where a
pressure clamp is active, its derivative is bounded by one; at its kink the
adjacent branches give the same conservative bound. The derivation below uses
the unclamped notation for clarity.

## Transverse Limiter Derivative

For four strictly same-sign slopes, the VL4 limiter is

$$
H=\frac{4}{1/a+1/b+1/c+1/d},\qquad
h_r=\frac{\partial H}{\partial a_r}=\frac{H^2}{4a_r^2}.
$$

The four slopes are the forward and backward temperature differences in each
of the two cells adjoining the face. Their six distinct temperature
coefficients have absolute sum

$$
\frac{\Gamma}{h_t}
=\frac{2\,[\max(h_a,h_b)+\max(h_c,h_d)]}{h_t},
\qquad 0\leq\Gamma\leq8.
$$

The implementation evaluates $h_r=4(w_r/\sum_s w_s)^2$, with
$w_r=\min_s|a_s|/|a_r|$, to avoid squaring extreme slopes. Equal nonzero
slopes give $\Gamma=1$. Strictly mixed signs give a locally zero limiter and
$\Gamma=0$, including when another slope is zero. Sign-consistent slopes
containing zero select the envelope $\Gamma=8$, covering neighboring
differentiable branches. A non-finite slope returns infinity and fails closed;
a non-finite primitive state remains an admissibility failure.

Representing the flux value with a bounded secant coefficient is insufficient
here: that coefficient's temperature derivative also contributes to the
Jacobian. The timestep uses limiter partial derivatives instead.

## Two Coupled Row Envelopes

Compute separate transverse factors for each temperature and define

$$
S_\parallel=\frac{2|b_n|}{h_n}
 +\sum_{t\ne n}\frac{|b_t|\Gamma_{\parallel,t}}{h_t},
\qquad
S_\perp=\frac{2|b_n|}{h_n}
 +\sum_{t\ne n}\frac{|b_t|\Gamma_{\perp,t}}{h_t}.
$$

Sums include only active dimensions. The normal two-point temperature
difference has row norm $2/h_n$. Let $g=(b\cdot\nabla B)/\overline B$ use the
same centered magnetic-magnitude gradients as the face closure, and set
$r=p_{\perp f}/p_{\parallel f}$. The grad-B response is

$$
f=T_{\perp f}\left(1-\frac{T_{\perp f}}{T_{\parallel f}}\right),
\qquad
\frac{\partial f}{\partial T_{\perp f}}=1-2r,
\qquad
\frac{\partial f}{\partial T_{\parallel f}}=r^2.
$$

At isotropy, $f=0$ but these derivatives do not vanish. Face temperatures
are density-weighted means with nonnegative weights summing to one; pressure
clamping can only reduce their derivative norms. Define

$$
P_f=\chi_{\parallel f}|b_{nf}|S_{\parallel f},
\qquad
Q_f=\chi_{\perp f}|b_{nf}|
\left[S_{\perp f}+\left(|1-2r_f|+r_f^2\right)|g_f|\right].
$$

The timestep conductivities include background collisions and omit additional
nonnegative limiter scattering, giving an upper envelope for a frozen face
coefficient. Flux caps are omitted from the bound: with their scale frozen,
the scalar capped-flux response has derivative magnitude at most one.

The energy and magnetic-moment fluxes give the receiving-cell row bounds

$$
R_{\parallel i}=\sum_{f\in\partial i}
 \frac{\rho_f}{\rho_i h_n}
 \left[P_f+2|1-\beta_{if}|Q_f\right],
\qquad
R_{\perp i}=\sum_{f\in\partial i}
 \frac{\rho_f}{\rho_i h_n}\beta_{if}Q_f.
$$

The factor $2(1-\beta)$ follows from recovering parallel temperature from
thermal energy and perpendicular temperature. Summing each face norm is
conservative even when different faces share temperature stencil entries.
Set

$$
\Lambda=\max_i\{R_{\parallel i},R_{\perp i}\},
\qquad \Delta t_{\rm ref}=\frac{2}{\Lambda}.
$$

A zero row adds no timestep restriction. The global timestep selection takes
the minimum reference budget over blocks, MPI ranks, and registered processes.
The refreshed post-sweep budget uses the current state and the same safety
factor as the initial selection.

## What The Bound Establishes

For the local uniform-grid temperature operator with density, magnetic field,
face conductivity, scattering factors and flux-cap scales frozen, these rows
bound the Jacobian infinity norm and hence its spectral radius. The grad-B
pressure-ratio derivative is included. The limiter kink envelope covers
neighboring differentiable branches.

A spectral-magnitude bound alone does not locate eigenvalues on the negative
real axis. Therefore `2 / Lambda` does not prove forward-Euler or RKL2 stability
for an arbitrary complex or positive spectrum. It is the conventional reference
step used by the parabolic controller, subject to numerical validation.

The bound does not include temperature derivatives of local conductivity or
cap scales, discontinuous hard-scattering switches, stage-varying coefficients,
the complete density/energy coordinate map, or interleaved floor and wall
projections. It also does not bound the complete AMR ghost/interpolation and
flux-correction operator. Nonlinear, refined-mesh, and full-physics acceptance
must be checked independently; this derivation is not a global nonlinear RKL2
theorem.

## Exceptional Arithmetic And Regression Cases

The normal path evaluates finite positive face rows directly. If an
intermediate conductivity, normalization, density ratio or complete row
underflows or overflows, a logarithmic path combines those factors before
taking the reciprocal. Thus a finite stiffness is retained when an intermediate
conductivity alone is not representable. That path uses the conservative
factor $(1+r)^2\geq|1-2r|+r^2$ and absolute magnetic-gradient sums. Invalid
positive face sound speeds yield a zero timestep and fail the existing mesh
timestep check. Configured `lf_k_parallel` and background `lf_c_parallel0`
must be finite and positive.

The focused regressions include:

- A divergence-free staggered checkerboard in two and three dimensions. For
  equal widths $h$, face amplitudes one and uniform cell-centered guide field
  $\epsilon$, its parallel row is
  $\Lambda=4\chi_\parallel(2/\epsilon^2+[3\mathrm{D}])/h^2$.
- A field-aligned density contact of contrast $C$. The low-density interface
  cell has one ordinary face and one high-density face, giving
  $\Lambda=\chi_\parallel(C+3)/h^2$, rather than using the maximum face twice.
- Magnetic reversals, post-sweep refresh after uniform heating, finite safety
  parameter validation, STS-only safety scaling, restart loading, and extreme
  rows whose complete stiffness remains representable.

These tests retain their scientific growth, floor and admissibility criteria.
Changed timestep expectations follow the face-row formulas and the independent
STS safety factor.
