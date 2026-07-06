# Auto-settled top-Mdot frame tracking implementation guide

## Scope

Implement a new frame-tracking mode that uses the upper-boundary mass flux as
both the settling detector and the physical velocity estimator:

```text
v_frame_target = Mdot_top / (rho_cold * A_top)
```

The existing fixed-window `mode=top_mdot` path should remain available for
comparison. The new path should be a separate mode, tentatively
`mode=top_mdot_auto`.

The first validation run should use the TRML setup at:

- resolution `128 x 128 x 256`;
- `xi = 1e3`, named as `xi1e3` in files and job names;
- `Mach_rel = 0.5`, keeping the current `velocity = 0.6454972243679028`;
- `chi = 10^1.75`, keeping `rho_cold = 56.23413251903491`;
- no labels such as `10x cooling`.

Assuming the current input convention `xi = t_shear / t_cool_min` and
`t_shear = 1.5491933384829668`, the xi=1e3 input should use:

```text
t_cool_min = 0.0015491933384829668
```

## Sign convention

The sign convention must be explicit in comments, startup logs, and
machine-readable diagnostics.

Definitions:

- `Mdot_top = integral_top rho * v3 * dA`, sampled in the upper active-cell
  layer.
- Positive `Mdot_top` means material moving in the positive x3 direction,
  outward through the top boundary.
- Negative `Mdot_top` means material moving in the negative x3 direction,
  inward/downward through the top boundary.
- `v_frame_target = Mdot_top / (rho_cold * A_top)` is the target lab-frame
  velocity of the moving coordinate frame.
- The stored fluid velocities are frame-relative. When the frame velocity
  changes by `Delta v_frame`, the code must add
  `fluid_boost = -Delta v_frame` to the stored fluid velocity.

Therefore:

```text
Delta v_frame = v_frame_target_now - v_frame_current
fluid_boost   = -Delta v_frame
```

Sanity checks:

- If `Mdot_top > 0`, then `v_frame_target > 0`; the frame moves upward and the
  stored fluid receives a negative boost during the ramp.
- If `Mdot_top < 0`, then `v_frame_target < 0`; the frame moves downward and the
  stored fluid receives a positive boost during the ramp.

This relationship is the main sign guardrail. Diagnostics must expose both
`v_frame_target_now` and `fluid_boost` so a sign error is visible immediately.

## Controller design

Use a simple state machine:

```text
observing -> locked/ramping -> locked/full-speed
```

### Observing phase

At every frame-tracker application, sample instantaneous `Mdot_top`.

Maintain two exponentially weighted moving averages:

```text
Mdot_fast
Mdot_slow
```

Use time-aware EMA updates:

```text
alpha_fast = 1 - exp(-dt / tau_fast)
alpha_slow = 1 - exp(-dt / tau_slow)
Mdot_fast += alpha_fast * (Mdot_top - Mdot_fast)
Mdot_slow += alpha_slow * (Mdot_top - Mdot_slow)
```

Initialize both filters to the first valid sampled value after
`mdot_auto_detect_start_time`.

Define the relative settling metric:

```text
settle_rel = abs(Mdot_fast - Mdot_slow)
             / max(abs(Mdot_slow), mdot_auto_abs_floor)
```

Declare the flux settled only if:

```text
time >= mdot_auto_detect_start_time
settle_rel <= mdot_auto_rel_tol
```

for a continuous duration:

```text
mdot_auto_hold_time
```

Suggested first defaults:

```text
mdot_auto_detect_start_time = 0.5
mdot_auto_fast_tau          = 0.5
mdot_auto_slow_tau          = 1.5
mdot_auto_rel_tol           = 0.075
mdot_auto_abs_floor         = 1.0e-12
mdot_auto_hold_time         = 0.5
```

These are intentionally configurable. The previous xi=100-style diagnostic run
showed `Mdot_top` approaching a near-plateau by roughly `t = 3--5`; the xi=1e3
case may settle earlier or later, so hard-coded times should be avoided.

### Lock and ramp phase

At the first settled time:

```text
t_settle = current time
Mdot_locked = Mdot_slow
v_frame_locked = Mdot_locked / (rho_cold * A_top)
v_frame_at_lock = current frame velocity
ramp_duration = clamp(mdot_auto_ramp_fraction * t_settle,
                      mdot_auto_ramp_min_time,
                      mdot_auto_ramp_max_time)
```

Use the user's intended scaling:

```text
mdot_auto_ramp_fraction = 0.25
```

So if the detector settles at `t = 5`, ramp over `1.25`; if it settles at
`t = 10`, ramp over `2.5`.

Suggested bounds:

```text
mdot_auto_ramp_min_time = 0.25
mdot_auto_ramp_max_time = 5.0
```

During the ramp:

```text
f = smoothstep((time - t_settle) / ramp_duration)
v_frame_target_now = v_frame_at_lock
                   + f * (v_frame_locked - v_frame_at_lock)
```

Then apply the existing Galilean boost machinery using:

```text
Delta v_frame = v_frame_target_now - frame_velocity_x3
fluid_boost   = -Delta v_frame
```

The existing slew limiter may still apply, but if it limits the ramp, the
history must say so through the existing `ft_limit`/`ft_skip` style fields.

## Proposed input interface

Add the following parameters under `<frame_tracking>` for `mode=top_mdot_auto`:

```text
mode = top_mdot_auto
mdot_density = 56.23413251903491
mdot_area = -1.0                    # default: infer Lx * Ly

mdot_auto_detect_start_time = 0.5
mdot_auto_fast_tau = 0.5
mdot_auto_slow_tau = 1.5
mdot_auto_rel_tol = 0.075
mdot_auto_abs_floor = 1.0e-12
mdot_auto_hold_time = 0.5

mdot_auto_ramp_fraction = 0.25
mdot_auto_ramp_min_time = 0.25
mdot_auto_ramp_max_time = 5.0
```

Existing generic controls still apply:

```text
apply_every = 1
diagnostic_every = 2000
axes = x3
max_abs_boost = 0.2
max_boost_change_mode = per_time
max_boost_change_rate = 0.02
```

For this auto mode, `start_time` should not be the settling trigger. It can
remain available for compatibility, but the auto mode should use
`mdot_auto_detect_start_time` and the measured settling condition.

## Diagnostics

Retain the existing per-axis fields:

```text
ft_vf_x3     # current lab-frame velocity of the moving frame
ft_dx_x3     # integrated frame displacement
ft_dv_x3     # instantaneous fluid boost added to stored v3
ft_limit     # slew/max-boost limiter active
ft_skip      # no boost applied
```

Add auto-top-Mdot fields. Keep labels short enough to remain readable in the
history header:

```text
ft_mtop      # instantaneous signed Mdot_top
ft_mfast     # fast EMA of Mdot_top
ft_mslow     # slow EMA of Mdot_top; used for lock
ft_mrel      # abs(Mfast-Mslow)/max(abs(Mslow), floor)
ft_hold      # current continuous settled-hold time
ft_settle    # 1 after the target has locked, otherwise 0
ft_tset      # lock time; -1 before lock
ft_mlock     # locked Mdot value; 0 before lock
ft_vlock     # locked target frame velocity; Mlock/(rho_cold*A)
ft_vtgt      # ramped target frame velocity for this update
ft_ramp      # smoothstep ramp factor, 0..1
```

Startup and diagnostic stdout should include a one-line sign contract:

```text
Mdot_top = integral rho*v3*dA at upper x3; positive is +x3/outward.
v_frame_target = Mdot_top/(rho_cold*A); fluid_boost = -Delta v_frame.
```

The validation script should explicitly verify the sign relationship:

```text
sign(ft_vlock) == sign(ft_mlock), unless ft_mlock == 0
sign(ft_dv_x3) == -sign(delta ft_vf_x3), unless delta ft_vf_x3 == 0
```

## Code-change checklist

1. Add enum and parser support:
   - add `kFTTopMdotAuto`;
   - accept `top_mdot_auto`, `auto_top_mdot`, and possibly
     `settled_top_mdot` as aliases;
   - keep existing `top_mdot` behavior unchanged.

2. Add state members:
   - auto detector configuration;
   - EMA state;
   - settled hold time;
   - locked time, locked Mdot, locked frame velocity;
   - frame velocity at lock;
   - current ramp factor and current target frame velocity.

3. Add restart persistence:
   - bump `kFrameTrackingStateVersion`;
   - write/read the auto-mode state;
   - initialize cleanly from older restart states.

4. Reuse `SampleTopMassFlux()`:
   - do not duplicate the flux sampling kernel;
   - confirm it samples only blocks adjacent to the upper x3 boundary;
   - keep the sign as `rho * v3 * dA`.

5. Factor the x3 Galilean boost application if useful:
   - both `top_mdot` and `top_mdot_auto` should use the same final sign path:
     `fluid_boost = -Delta v_frame`;
   - update momentum and ideal-gas total energy using the existing formula.

6. Extend `FillHistoryData()`:
   - include auto-mode diagnostics;
   - ensure `NHISTORY_VARIABLES` capacity is sufficient;
   - verify non-root ranks still return zeros for frame-tracker history.

7. Extend validation and documentation comments:
   - add sign comments near `SampleTopMassFlux()` and the boost application;
   - add a short recipe note to the relevant frame-tracking docs after the
     code path is validated.

## Test input and run naming

Create a dedicated test input under the Lustre TRML simple input area, not by
overwriting the canonical input:

```text
/lustre/orion/ast207/proj-shared/dfielding/TRML/simple/inputs/
  TRML_xi1e3_M0p5_chi10p1p75_topmdotauto_128x128x256.athinput
```

Suggested basename:

```text
TRML_xi1e3_M0p5_chi10p1p75_topmdotauto_128x128x256
```

Required input changes from the canonical guide input:

```text
<mesh>
nx1 = 128
nx2 = 128
nx3 = 256

<problem>
rho_cold = 56.23413251903491
velocity = 0.6454972243679028
t_cool_min = 0.0015491933384829668

<frame_tracking>
mode = top_mdot_auto
mdot_density = 56.23413251903491
...
```

Use meshblocks appropriate for the requested Frontier allocation. For a
one-node/eight-GPU smoke run, an eight-block layout is desirable if it is stable
and memory-safe; otherwise record the actual decomposition in the run script
and final status.

Output, script, and log names should include `xi1e3`, not `10x`.

## Verification plan

### Build verification

After implementation:

```text
/ccs/home/dfielding/athenak-trml-tracers-tracking/configure.sh
```

The configure script should build, not only configure. Save the executable with
a unique name such as:

```text
/lustre/orion/ast207/proj-shared/dfielding/TRML/simple/athena_topmdotauto
```

Record a checksum before submitting jobs:

```text
sha256sum athena_topmdotauto
```

### Runtime verification

Use the debug queue for the first smoke run if available. The target smoke run
should be short enough to verify settling and ramping, for example through
`t = 12--15`, unless the detector has not locked by then.

Check the frame-tracker history immediately after the run starts:

- `ft_mtop`, `ft_mfast`, and `ft_mslow` finite;
- `ft_mrel` decreases when the flux plateaus;
- `ft_tset` appears at a physically sensible time;
- `ft_vlock = ft_mlock / (rho_cold * A_top)`;
- `ft_vtgt` ramps smoothly from `v_frame_at_lock` to `ft_vlock`;
- `ft_dv_x3 = -Delta ft_vf_x3` up to output cadence and limiter effects.

### Physics sanity checks

Compare against the previous fixed-window top-Mdot run and the no-frame lower
BC comparison:

- hot reservoir temperature remains close to the intended hot phase;
- no perfectly horizontal, frame-induced ridge appears;
- the frame-tracker locks before the mixed layer reaches the top of the box;
- lower boundary behavior matches the selected `lower_x3_user_bc`;
- mass and energy histories show no abrupt jump at lock time or ramp start.

## Validated smoke-run recipe

The first validated smoke run used:

```text
executable:
  /lustre/orion/ast207/proj-shared/dfielding/TRML/simple/athena_topmdotauto
input:
  /lustre/orion/ast207/proj-shared/dfielding/TRML/simple/inputs/
    TRML_xi1e3_M0p5_chi10p1p75_topmdotauto_128x128x256.athinput
script:
  /lustre/orion/ast207/proj-shared/dfielding/TRML/simple/run_scripts/
    run_TRML_xi1e3_topmdotauto_128_batch_smoke.sh
output:
  /lustre/orion/ast207/proj-shared/dfielding/TRML/simple/data/
    TRML_xi1e3_M0p5_chi10p1p75_topmdotauto_128x128x256_batch_smoke
```

The debug queue was attempted first, but Slurm rejected debug access for this
account, so the smoke run used `batch` on one node with eight ranks/eight GPUs.
The test input used `128 x 128 x 256` zones, `64 x 64 x 128` meshblocks,
`xi = 1e3`, `Mach_rel = 0.5`, and `chi = 10^1.75`. The lower user-x3 boundary
in this smoke input was `outflow`; use that fact when comparing against
reflecting-lower-boundary runs.

Job `4942156` completed to `tlim = 15.0` with Slurm state `COMPLETED` and exit
code `0:0`. The frame tracker locked at:

```text
t_settle = 5.9033561356249935
Mdot_locked = -0.27008543772526483
v_frame_locked = -0.0048028737285818815
ramp_end_first_history_time = 7.500264628933981
```

The sign checks passed:

```text
v_frame_locked = Mdot_locked/(rho_cold*A_top)
sign(Mdot_locked) == sign(v_frame_locked)
same-sign count for Delta ft_vf_x3 and ft_dv_x3 = 0
```

The final frame-tracker history row at `t = 15.0` had `ft_ramp = 1`,
`ft_skip = 1`, `ft_vf_x3 = ft_vtgt = ft_vlock`, and no further boost applied.

## Judicious subagent use

Use subagents for bounded, independently checkable tasks. Avoid parallel edits
to the same source file. The main agent should remain the integrator that
decides what gets patched, built, submitted, and reported.

Recommended subagent assignments:

| Subagent role | When to use | Scope | Should edit files? |
| --- | --- | --- | --- |
| Sign-audit reviewer | Before or immediately after implementation | Trace `Mdot_top`, `v_frame_target`, `Delta v_frame`, `fluid_boost`, momentum, and energy signs. Produce a written sign report. | No |
| Patch reviewer | After the main implementation patch exists | Review `frame_tracker.hpp/cpp` for state initialization, restart persistence, history capacity, and mode separation. | No, unless explicitly assigned a non-overlapping fix |
| Build runner | After patch review | Run configure/build, capture compiler errors, and report exact failing lines. | No |
| Input/script reviewer | After the test input and Slurm script are drafted | Check names use `xi1e3`, resolution is `128x128x256`, `Mach_rel=0.5`, `chi=10^1.75`, and `t_cool_min` follows the xi convention. | Prefer no; if editing, only input/script files |
| Runtime monitor | After submission | Watch `squeue`, log output, and early `.frame_tracker.hst`; report trigger/ramp status. | No |
| Verification analyst | After run completes or reaches enough time | Parse history data and compare sign relationships and physics sanity checks. | No |

Subagent prompt templates:

```text
Sign-audit reviewer:
Inspect src/srcterms/frame_tracker.cpp and the proposed top_mdot_auto patch.
Do not edit files. Verify the sign convention from Mdot_top sampling through
frame velocity update, stored fluid boost, momentum update, energy update, and
history labels. Report any sign ambiguity or mismatch.
```

```text
Patch reviewer:
Review the top_mdot_auto implementation only. Do not edit files. Check parser
aliases, defaults, restart state, history capacity, old top_mdot preservation,
and behavior for older restart files. Report concrete issues with file/line
references where possible.
```

```text
Runtime monitor:
Monitor the submitted xi1e3 128x128x256 top_mdot_auto job. Do not edit files.
Report queue state, log errors, and the first available frame-tracker history
rows showing Mdot, lock time, ramp factor, target frame velocity, actual frame
velocity, and fluid boost.
```

Subagent outputs should be treated as evidence to integrate, not as automatic
authority. The main agent must still inspect critical code paths and verify the
final run.

## Stop conditions

Stop and report before launching a production-length run if any of the
following occurs:

- the sign audit finds a mismatch between target frame velocity and fluid
  boost;
- the code changes the old `top_mdot` behavior unintentionally;
- build succeeds but the history fields do not expose enough information to
  verify the sign convention;
- the smoke run locks at an obviously pathological time, such as immediately at
  the first sample or only after the mixed layer has reached the top boundary;
- the smoke run shows a discontinuous mass/energy jump at lock or ramp start.
