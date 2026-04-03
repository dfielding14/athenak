Below is a “spec + implementation playbook” you can hand directly to an AI coding agent. It’s written so the agent can work *inside this one pgen file* (`scalar_mixing.cpp`) and only touch other files if absolutely necessary.

I’ll assume:

* You evolve **conserved** scalar variables (q_s \equiv \rho,\theta_s) stored in `u0(m, nhydro + s, ...)`.
* The corresponding **primitive** scalar (\theta_s) is available in `w0(m, nhydro + s, ...)` (as your history routine already assumes).
* BCs are triply periodic.
* (\rho) and (\mathbf u) are frozen when `scalar_only=true`.

The four forcing methods to implement (per-scalar selectable):

1. Mean gradient forcing
2. Mass-conserving Allen–Cahn (bistable, (\theta_1>0), (\theta_2) unstable, (\theta_3) stable, and **mass conservation**)
3. Two-reservoir fringe/sponge forcing (no global reduction)
4. Reservoir-forced reaction–diffusion (sponge + local KPP or bistable reaction; no global reduction)

---

# A. High-level design decision (important)

AthenaK only lets you enroll **one** `user_srcs_func`. So implement a **single** source term function that:

* Loops over all scalars `ns=0..nscalars-1`
* Applies the chosen method for each scalar based on per-scalar runtime config

You will replace:

```cpp
user_srcs_func = MeanGradientForcing;
```

with:

```cpp
user_srcs_func = ScalarForcingSource;
```

(where `ScalarForcingSource` implements all four methods).

Keep the old behavior as the default for scalar 0 so existing input files still work.

---

# B. Add per-scalar configuration (new struct + enums)

## B1) Add enums and a config struct in the anonymous namespace

Near your existing:

```cpp
namespace {
  Real mean_gradient_G = 1.0;
}
```

Replace/extend with something like:

* `ScalarForcingMode` enum (0..4)
* `ReactionType` enum (for method 4)
* `ScalarForcingConfig` struct holding **all** parameters needed by any method
* A global `std::vector<ScalarForcingConfig>` sized to `nscalars`

Suggested fields (keep POD/simple types so it can be captured into Kokkos lambdas):

* `int mode;`  // 0 none, 1 mean-grad, 2 mc-AllenCahn, 3 sponge, 4 sponge+reaction
* Mean gradient: `Real G;`
* Reaction (AC/bistable): `Real tau, theta1, theta2, theta3; Real a;` (optional `a` and derive `theta2`)
* Sponge: `Real chiL0, chiR0; Real thetaL, thetaR; Real xL0, xL1, xR0, xR1; Real mask_w;`
* Method 4 only: `int reaction_type;` // 1 KPP, 2 bistable
* Optional safety: `Real theta_floor; bool clip_reaction_interval;` (recommended for KPP)

Also store a global “have we initialized configs” flag if useful.

**Naming convention for inputs** (recommended):

* `problem/scalar0_mode`, `problem/scalar1_mode`, ...
* And similarly `scalar0_tau`, `scalar0_theta1`, `scalar0_chiL0`, ...

This is much easier than trying to parse arrays/lists from the input file.

---

# C. Parse inputs in `UserProblem` for **all scalars** (including restart)

## C1) Always read per-scalar parameters before `if (restart) return;`

Right now you read `mean_gradient_G` and enroll sources before you `return` on restart — good. Extend that logic so the per-scalar config vector is filled even on restart.

Pseudo steps:

1. Get `nscalars = pmbp->phydro->nscalars;`
2. Resize global config vector: `cfg.resize(nscalars);`
3. For each scalar `ns`:

   * read mode (default: `ns==0 ? 1 : 0` so scalar0 remains mean-grad by default)
   * read parameters with defaults

### Defaults (reasonable)

These defaults aim to not break positivity:

* `theta1 = 1e-6` (must be >0)
* `theta3 = 1.0`
* `a = 0.5` and set `theta2 = theta1 + a*(theta3-theta1)` unless user supplies `theta2`
* `tau = 1.0`
* Sponge:

  * `chiL0 = chiR0 = 1.0` (relaxation rate)
  * `thetaL = theta1`, `thetaR = theta3`
  * slabs default somewhere like `[0.20Lx, 0.30Lx]` and `[0.70Lx, 0.80Lx]`
  * smoothing `mask_w = 0.02Lx`
* Method 4:

  * `reaction_type = 2` (bistable default; KPP is more fragile for positivity)

### Important: slab coordinates

Use physical coordinates in the same units as the mesh. Compute `Lx = x1max - x1min` from `pm->mesh_size` in `UserProblem` for default placement.

---

# D. Initialize multiple scalars

You already do:

```cpp
for (int ns = 0; ns < nscalars; ns++) {
  u0(m, nhydro+ns, k, j, i) = rho * scalar_init;
}
```

Modify so each scalar can have its own init:

* Input: `problem/scalar0_init`, `problem/scalar1_init`, ...
* Default: fall back to existing `scalar_init`

Also enforce `scalar_init > 0` (at least the configured `theta_floor`) because you said the code cannot handle `theta <= 0`.

---

# E. Implement the single source term function for all 4 methods

## E1) New function signature

Add at top:

```cpp
void ScalarForcingSource(Mesh* pm, const Real bdt);
```

Enroll it:

```cpp
user_srcs_func = ScalarForcingSource;
```

## E2) Indexing rules (critical)

For scalar `ns`:

* Primitive value: `theta = w0(m, nhydro + ns, k, j, i)`
* Conserved variable update: `u0(m, nhydro + ns, k, j, i) += rho * source(theta, x) * bdt`

This matches your current mean-gradient forcing style.

---

# F. Method-specific formulas & implementation details

Below are the exact source terms the agent should implement.

---

## Method 1: Mean gradient forcing (existing, generalized to multiple scalars)

**Per scalar** parameter: `G_s`

**Source:**
[
S(\mathbf x,\theta_s) = G_s,v_x(\mathbf x)
]
and in conservative form:
[
\partial_t(\rho \theta_s) ;+=; \rho,G_s,v_x
]

**Kernel skeleton:**

* Get `density = w0(IDN)`
* Get `vx = w0(IVX)`
* `u0(nhydro+ns) += density * G_s * vx * bdt`

---

## Method 2: Mass-conserving Allen–Cahn (bistable + global reduction)

**Reaction (bistable cubic with (\theta_1>0)):**
[
R(\theta)= -\frac{1}{\tau}(\theta-\theta_1)(\theta-\theta_2)(\theta-\theta_3)
]
with (\theta_1<\theta_2<\theta_3), (\theta_1>0).

**Mass-weighted mean (because you evolve (\rho\theta)):**
[
\langle R\rangle_\rho = \frac{\int \rho,R(\theta),dV}{\int \rho,dV}.
]

**Mass-conserving source:**
[
S(\mathbf x,\theta)=R(\theta)-\langle R\rangle_\rho.
]

**Implementation steps per scalar:**

1. Compute numerator and denominator by reduction over all interior cells on this rank:

   * `num_local = Σ rho * R(theta) * vol`
   * `den_local = Σ rho * vol`
2. MPI_Allreduce to get global sums, compute:

   * `Rbar = num_global / den_global`
3. Apply update:

   * `u0 += rho * (R(theta) - Rbar) * bdt`

### Notes / pitfalls to tell the agent

* This requires **one global reduction per scalar per source call** (unless you batch scalars into one reduction).
* Because of the `-Rbar` shift, the *exact* fixed points are not strictly (\theta_1) and (\theta_3) everywhere; bulk phases settle near values solving (R(\theta)=Rbar). That’s expected for mass-conserving Allen–Cahn. (But equilibria will still be near (\theta_1,\theta_3) if the mixture fractions aren’t extreme.)
* Positivity: choose (\theta_1) not too tiny (e.g. (10^{-6}) or larger). If you clip/floor after update, you will break strict mass conservation.

---

## Method 3: Two-reservoir fringe/sponge forcing (no global reduction)

You want two thin slabs in x with smooth masks (\chi_L(x)), (\chi_R(x)).

**Source:**
[
S(\mathbf x,\theta)= -\chi_L(x)(\theta-\theta_L) - \chi_R(x)(\theta-\theta_R)
]
with (\theta_L,\theta_R>0).

### Smooth mask implementation (recommended)

Define a smooth top-hat window on ([x_0,x_1]):

[
W(x;x_0,x_1,w)=\tfrac12\left[\tanh!\left(\frac{x-x_0}{w}\right)-\tanh!\left(\frac{x-x_1}{w}\right)\right]
]

Then
[
\chi_L(x)=\chi_{L0},W(x;x_{L0},x_{L1},w),\quad
\chi_R(x)=\chi_{R0},W(x;x_{R0},x_{R1},w).
]

### Kokkos implementation detail

You’ll need to compute cell-center `x` in the source kernel, exactly like you do in `InitTurbulentVelocity`:

```cpp
Real &x1min = size.d_view(m).x1min;
Real &x1max = size.d_view(m).x1max;
Real x = CellCenterX(i-is, nx1, x1min, x1max);
```

Then compute `chiL, chiR` and apply update.

**Positivity benefit:** This forcing always relaxes toward positive reservoir values.

---

## Method 4: Reservoir-forced reaction–diffusion (sponge + local reaction; no global reduction)

**Source:**
[
S(\mathbf x,\theta)= -\chi_L(x)(\theta-\theta_L) - \chi_R(x)(\theta-\theta_R) + R(\theta).
]

Two reaction choices:

### 4a) KPP/logistic between (\theta_1) and (\theta_3) (monostable)

Use a shifted logistic in (\theta\in[\theta_1,\theta_3]):
[
R_{\text{KPP}}(\theta)=\frac{1}{\tau(\theta_3-\theta_1)}(\theta-\theta_1)(\theta_3-\theta).
]

**Important:** for (\theta<\theta_1), this becomes negative and would push (\theta) downward (bad for positivity).
So implement KPP with **clamping** in the source evaluation:

* `theta_eff = min(max(theta, theta1), theta3)`
* evaluate `R_KPP(theta_eff)`

This avoids the reaction driving (\theta) below (\theta_1).

### 4b) Bistable (Allen–Cahn-type), shifted so (\theta_1>0) (recommended default)

Same cubic as Method 2, but **without** the mean subtraction:

[
R_{\text{bi}}(\theta)= -\frac{1}{\tau}(\theta-\theta_1)(\theta-\theta_2)(\theta-\theta_3).
]

**Why this is attractive:** it is self-bounding: for (\theta<\theta_1), the cubic is positive and pushes upward.

---

# G. Minimal code organization the agent should follow

Implement helper functions usable inside Kokkos lambdas:

1. `SmoothTopHat(x, x0, x1, w)` — marked `KOKKOS_INLINE_FUNCTION`
2. `BistableReaction(theta, tau, t1, t2, t3)` — `KOKKOS_INLINE_FUNCTION`
3. `KppReaction(theta, tau, t1, t3, clip=true)` — `KOKKOS_INLINE_FUNCTION`

Then in `ScalarForcingSource`:

* Get references: `u0`, `w0`, `size`, `indcs`, `nhydro`, `nscalars`
* Loop `ns`
* Switch on `cfg[ns].mode` and run the appropriate kernel(s)

For Method 2, do:

* reduction (Kokkos + MPI) → `Rbar`
* then update kernel

---

# H. Concrete “diff-style” guidance (safe edits only)

These are edits the agent can apply confidently inside **this file** without touching the solver/integrator.

### H1) Enroll the new source function

In `UserProblem`:

```diff
-  user_srcs_func = MeanGradientForcing;
+  user_srcs_func = ScalarForcingSource;
```

### H2) Generalize scalar init

In `pgen_scalar`, change to per-scalar init values (agent will implement reading `scalarN_init`).

### H3) Add new function prototype

At top:

```diff
-void MeanGradientForcing(Mesh* pm, const Real bdt);
+void ScalarForcingSource(Mesh* pm, const Real bdt);
```

(You can keep `MeanGradientForcing` if you want, but it becomes unused.)

---

# I. Input file examples (give these to the agent as acceptance tests)

### Example 1: scalar0 mean gradient; scalar1 mass-conserving AC

```ini
<hydro>
scalar_only = true
nscalars = 2

<problem>
# Scalar 0: mean gradient
scalar0_mode = 1
scalar0_mean_gradient = 1.0
scalar0_init = 0.5

# Scalar 1: mass-conserving Allen–Cahn
scalar1_mode = 2
scalar1_init = 0.5
scalar1_tau = 1.0
scalar1_theta1 = 1e-3
scalar1_theta3 = 1.0
scalar1_a = 0.5     # theta2 = theta1 + a*(theta3-theta1)
```

### Example 2: sponge reservoirs only

```ini
<hydro>
scalar_only = true
nscalars = 1

<problem>
scalar0_mode = 3
scalar0_init = 0.2
scalar0_thetaL = 0.1
scalar0_thetaR = 1.0
scalar0_chiL0 = 2.0
scalar0_chiR0 = 2.0
scalar0_xL0 = 0.2
scalar0_xL1 = 0.3
scalar0_xR0 = 0.7
scalar0_xR1 = 0.8
scalar0_mask_w = 0.02
```

### Example 3: sponge + bistable reaction

```ini
<problem>
scalar0_mode = 4
scalar0_reaction_type = 2   # 2=bistable, 1=KPP
scalar0_tau = 0.5
scalar0_theta1 = 1e-3
scalar0_theta3 = 1.0
scalar0_a = 0.3
# plus sponge params as above
```

---

# J. Validation checklist (tell the agent to run these)

1. **Regression:** With `nscalars=1` and `scalar0_mode=1`, results should match the current mean-gradient case.
2. **Mass conservation test (Method 2):**

   * Turn off diffusion and set `u=0` (or very small).
   * Verify (\int q,dV) is constant to roundoff over time.
3. **No global takeover:** For Method 2, start near the unstable point and verify the final state is a two-phase mixture consistent with the conserved mean.
4. **Positivity:** Stress test with stronger forcing; ensure no scalar hits (\theta\le 0). If it does, increase `theta1`, decrease `bdt`, decrease forcing rates, or add reaction clamping for KPP.
