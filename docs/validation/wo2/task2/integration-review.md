# Independent Task2 integration review

Read-only review of the uncommitted Task2 driver change over Task1 commit
822cad5e8a49f5f510bece9dcb4b601cacc6f205 on 2026-10-06. No applications were run,
and no production source was changed during this review.

No missing rollback, output, or MPI dependency was found within the explicit
eligibility gate: uniform strictly periodic active CGL LF with one RKL2 process,
zero collision rates, no configured limiters or strict-admissibility mode, and
no additional field updates, user source/boundary callbacks, or excluded physics.
This review does not extend that scope or establish nonlinear stability.

- Trial state includes the complete conserved/primitive arrays, including ghosts,
  plus all mesh EOS events, LF diagnostics, A/mu representation and both local
  timestep estimates. Eligible LF task lists do not change face/cell magnetic
  fields, coarse maps, forcing state, or output scheduling.
- Stage 1 resets the private heat-work recurrence; AddHeatFluxes replaces its
  stage work. Hyperbolic pressure-work bookkeeping is not touched by the trial.
  Invalidating the fused primitive cache and zeroing u_sts1/u_sts_rhs prevents
  rejected nonfinite values from entering zero-coefficient recurrence terms.
  Other STS registers are repopulated from restored u0 or those cleared values.
- Each stage completes its ordinary send/receive clearing tasks before returning.
  The MPI_BOR rejection decision makes the accepted/rejected path collective.
  Timing/profile counters intentionally retain attempted work; physical/EOS
  diagnostics are restored on rejection.
- The deferral predicate duplicates the actual float32 time and cycle output
  tests. No half-sweep debt is serialized: output, final-time/cycle exits and
  unexpected wall-clock exits flush first. Wall time is already broadcast across
  ranks. Finalization repeats the flush defensively.
- A rejected trial restores local budgets, completes the old pending half and
  selects the next timestep from that synchronized state using the old completed
  cycle for the growth cap, before the ordinary new pre sweep. Accepted merged
  sweeps refresh local advective/LF budgets and enforce the fresh advective limit.
  Forcing still runs only afterward in its original before-timeintegrator task.

Source reviewed: src/driver/driver.cpp and driver.hpp, the enrolled MHD STS task
sequence, mhd_sts.cpp state/register lifecycle, LF private diagnostic recurrence,
CGL collision hooks, and mesh timestep selection. Existing runtime acceptance and
performance evidence are owned by the Task2 implementation report.
