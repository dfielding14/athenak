# GPU Harris-sheet comparison plan

Date: 2026-09-26. Branch: `fast_reconnection`.
Status: planned; the GPU campaign and movies have not been run or produced.

## Question and scope

Measure how the current-limited resistivity changes the evolution and
reconnection rate of an otherwise identical Harris sheet. Run separate copies
of the same single-sheet initial condition, with exactly the same mesh and
background Ohmic diffusivity. Produce synchronized current-density movies and
reconnection-rate curves for the control and a small parameter scan.

Here “same Ohmic resistivity” means the same background
$\eta_0=10^{-6}$. The control has $\eta=\eta_0$ everywhere. The new model uses
the total coefficient $\eta(q)$, which already contains the background; do not
add $\eta_0$ again. The equations and implementation are documented in
[current-limited resistivity](current_limited_resistivity.md).

The [completed local pilot](../reports/status-update-2026-09-26/STATUS_UPDATE.md)
establishes CPU checks and diagnostic limitations. Its uniform-resistivity
control used $5\times10^{-6}$, so it is not the matched-background control
requested here. Its coarse continuation reached only $q_{\max}=0.2266$ by
$t=5$, without entering the upper branch. Neither GPU parity nor resolved fast
reconnection has been demonstrated.

## Four primary runs

| Case | Model | `ohmic_resistivity` | `d_i` | `eta_max` | Comparison |
|---|---|---:|---:|---:|---|
| H0 | `constant` | $10^{-6}$ | diagnostic only | unused | Matched Ohmic control |
| H1 | `current_limited` | $10^{-6}$ | 0.005 | 0.005 | Existing fiducial physics |
| H2 | `current_limited` | $10^{-6}$ | 0.020 | 0.005 | Larger ion scale at fixed cap |
| H3 | `current_limited` | $10^{-6}$ | 0.020 | 0.001 | Lower cap at fixed ion scale |

Use `b_rec_method=constant` and `b_rec=1` throughout. H1 versus H2 tests the
ion scale; H2 versus H3 tests the cap. All three share H0 because these closure
parameters do not change the initialized density, pressure, velocity, or field.
For H0, label any reconstructed $q$ by the chosen diagnostic ion scale; it has
no dynamical ion scale.

The larger-ion-scale cases start below the upper branch, but may reach it more
readily. This is a hypothesis to test. Changing only `eta_max` leaves the
constitutive law unchanged while both solutions remain below their respective
$q_*=1-\sqrt{\eta_0/\eta_{\max}}$, though timestep/stage counts and numerical
errors can still change. The model already enhances diffusivity for
$0<q<q_*$; “upper-branch onset” must not be described as an on/off switch.

Defer guide-field and jump-estimator scans. Changing the mesh, background
resistivity, seed, sheet width, or other initial conditions requires a newly
matched control. A later ion-scale scan can reuse H0 on this same grid.

## Shared setup

Start from [the Harris input](../inputs/mhd/resistive_harris.athinput), with
these common settings and only the table's model parameters changed:

| Setting | Value |
|---|---|
| Mesh | $1600\times800\times1$, uniform, no AMR |
| Domain | $x\in[-0.5,0.5]$, $y\in[-0.25,0.25]$, unit depth |
| Boundaries | Outflow in $x,y$; periodic in the inactive $z$ direction |
| Initial sheet | $B_0=1$, width 0.02, guide field 0 |
| Gas | $\rho_{\rm up}=1$, sheet density enhancement 1, $p_{\rm up}=0.1$, $\gamma=5/3$ |
| Seed | Gaussian flux perturbation $10^{-4}$, width 0.1 |
| Numerics | Double precision, RK2, PLM, HLLD, CFL 0.4 |
| Diffusion integration | `resistivity_integrator=sts`, `sts_integrator=rkl2`, `sts_max_dt_ratio=32` for all four |
| Initial meshblock choice | $40\times40$; tune once on the GPU and keep the selected layout common |

Here $v_{A0}=B_0/\sqrt{\rho_{\rm up}}=1$ and time is in $L/v_{A0}$ with
$L=1$. The spacing is $0.000625$: H1 has eight cells per upstream ion length
and about 5.7 per initial central ion length. The relevant local length is
$d_i/\sqrt{\rho}$, so count cells in the evolving active sheet, not just upstream.

The same integration settings do not imply identical timesteps or STS stage
counts: the current-limited timestep bound uses `eta_max`. Record actual
timesteps and stages. Repeat H0 and an informative model case at CFL 0.2 over
the same physical interval before attributing rate differences to the closure.

## GPU execution sequence

1. Fetch `origin/fast_reconnection`; build Release with the GPU backend and
   architecture appropriate to the machine, with `Athena_SINGLE_PRECISION=OFF`.
   Save the source revision, full configuration, compiler/runtime versions,
   executable hash, GPU model, rank/device mapping, and complete inputs.
2. Run the existing diffusion/nonlinear STS checks, including the ratio-32
   cases, and a short CPU/GPU Harris parity comparison. Use numerical
   tolerances appropriate to backend reduction order, not bitwise identity.
   Verify that all four new cases initialize identical physical state arrays.
3. Run a short GPU timing and output pilot. Measure memory, disk use, wall time,
   cycles, and STS stages before committing to the full series. Do not infer
   an STS speedup from stage counts; that requires matched explicit timings.
4. Evolve every case to the same physical horizon, initially $t=20$ in
   checkpointed chunks of five. Inspect thinning, topology, local resolution,
   and upper-branch occupancy between chunks. Extend the common horizon if
   onset or a useful post-onset interval has not yet been captured and measured
   cost permits it. There is no guarantee of onset by $t=20$.
5. Refine H0 and an informative model case together to $3200\times1600$ over
   an overlapping physical interval, starting from the same analytic initial
   condition. Assess onset, flux histories, and rate statistics against this
   refinement and the timestep/cadence checks.

Keep the original seed and initial sheet for this series. If the whole family
fails to reach the upper branch within an affordable horizon, report that
bounded result. A thinner sheet or stronger seed would start a separate,
fully matched series. Do not tune only the model-on run to obtain a desired rate.

## Movies and saved outputs

Make a control-versus-H1 movie first, then H0-versus-H2 and H0-versus-H3 using
the same rendered control frames. A four-panel overview is useful if legible.
The primary field is signed $J_z$, with magnetic-flux contours; mark X/O points
only when their classification has been verified. Add $q$ and $\eta$ views to
explain where the model changes the sheet. Label case, parameters, physical
time, and resolution.

Use identical axes, aspect ratio, and a common fixed symmetric color scale
for $J_z$ across cases and frames. If its dynamic range requires a symmetric
log scale, state the linear threshold. Use common scales for corresponding
$q$ and $\eta$ panels too. Never use independent per-frame autoscaling.

Proposed starting cadences are field dumps every 0.05, histories every 0.001,
and restarts every 1. Save full, unsharded 2-D `mhd_u_bcc` state dumps. Reuse
the existing binary reader and reconstruct the same $J_z$ diagnostic in every
case; `mhd_jz` is also available as a derived output. Diagnostic cell currents
are proxies for the native CT-edge currents used by the closure. Retain edge
maxima from histories and selected `mhd_eta`/`mhd_q`/divergence outputs for QA.
`mhd_q` is unavailable for H0; reconstruct its diagnostic $q$ if needed.

Output intervals are crossed by timesteps; dumps need not land at exactly
the requested times. Synchronize by actual timestamps, never by frame index.
For nearest-frame pairing, show the time offsets and require offsets below
5% of the movie cadence; increase output cadence if this fails. If visual
interpolation is used, label it and keep rate measurements on the original
snapshots. Increase field cadence around rapid events and test rate sensitivity
to cadence; dense histories alone cannot validate a sparse flux derivative.

At $t=20$, the proposed cadence gives about 401 frames, or a 20-second movie at
20 fps. Eight double-precision state variables require approximately 82 MB
per frame: roughly 33 GB per case and 132 GB for four, before headers, extra
scalar fields, and restarts. Check available storage in the GPU pilot; adjust
the common cadence or optional outputs together. Keep raw data outside Git.

## Reconnection-rate measurements

Deliver signed flux and rate time series for every case, using $B_0L$ for
flux normalization and $B_0v_{A0}$ for rate normalization. Show the following
two diagnostics separately:

1. **Fixed boundary-reference flux.** While the central point remains an X
   point, use the existing diagnostic
   $$\psi_{\rm ref}(t)=A_z(x_{\rm ref},0)-A_z(0,0)
   =-\int_0^{x_{\rm ref}}B_y(x,0)\,dx,$$
   $$R_{\rm ref}=\frac{d\psi_{\rm ref}/dt}{B_0v_{A0}}.$$
   Independently compare with
   $[E_z(X)-E_z(x_{\rm ref},0)]/(B_0v_{A0})$, including both ideal and resistive
   terms. The outflow boundary is a fixed flux reference, not an O point, and
   its electric field is generally nonzero. Retain the flux series if central
   topology changes, but flag its loss of a central-X interpretation.
2. **Verified X–O flux.** Once an island exists, reconstruct $A_z$, identify
   saddle/extremum pairs, and follow a consistently identified reconnecting
   sheet and island. Use $\psi_{XO}=A_z(O)-A_z(X)$ and
   $R_{XO}=(d\psi_{XO}/dt)/(B_0v_{A0})$, with an independent
   $[E_z(X)-E_z(O)]/(B_0v_{A0})$ check at the actual nulls. Save point positions
   and identities. Segment curves at mergers, exits, or pair reassignment;
   never differentiate across an identity-change jump. Mark intervals with
   no valid pair as unmeasured, rather than inventing an O point.

Keep signs in the CSV and state the convention in plots. If showing $|R|$,
label it explicitly. Save unsmoothed fluxes and document derivative/smoothing
choices. Do not infer reconnection from $\eta J_z$ alone. Existing histories
reconstruct physical EMFs, not the full numerical CT EMF; their disagreement
with flux changes is a numerical-error diagnostic. The local startup runs
showed substantial spatial contamination despite small timestep sensitivity.

Alongside rates, plot edge $q_{\max}$, diffusivity, sheet thickness, local ion
length in cells, upper-branch occupancy, and energy/heating diagnostics.
Control upper-branch fractions are not directly meaningful; recompute with
the corresponding model threshold if a comparison is needed. Check finite
states, positive density/pressure, and CT divergence. Open-boundary energy
changes are expected; the existing approximate boundary-flux ledger is not
an exact conservation proof.

## Work to do on the GPU machine and deliverables

The implementation and CPU evidence are ready to transfer. The following
campaign-specific work remains:

- Generate four audited input decks from one common template. Do not use
  `benchmarks/reconnection/run_local.py --model constant` unchanged: it sets
  $\eta=5\times10^{-6}$ and uses explicit diffusion. Its mesh argument is also
  tied to the fiducial ion scale. Use explicit fixed mesh dimensions for H0–H3.
- Reuse the existing analyzer and comparison code, but give every case a
  unique ID containing its model parameters. `compare.py` currently keys
  comparisons by model, resolution, and CFL; several cap/ion-scale variants
  would collide. Preserve the historical CPU report and its input files.
- Add synchronized movie rendering using the existing binary reader,
  Matplotlib, and FFmpeg. The requested movie renderer does not yet exist.
- Extend the current central/midplane topology checks to track actual X/O
  pairs and sample their physical EMFs when they move off those locations.
  The current tools do not provide full 2-D null tracking.

Deliver MP4 comparisons, rate/flux CSVs and figures, diagnostic plots, exact
inputs and provenance, and a short interpretation over a common time window.
Report onset and transient behavior first; quote a steady mean only if a
stationary interval actually appears. A physical enhancement claim requires
stability under spatial, timestep, and output-cadence refinement, with at least
4–5 cells per local ion length in the active layer as an initial resolution
screen, not a substitute for convergence. There is no prescribed target rate.

All code and retained small CPU results travel with this branch. Ignored local
builds, raw dumps, and restart files do not. Start the matched GPU series from
its common initial condition; transfer old raw data separately only if needed.
