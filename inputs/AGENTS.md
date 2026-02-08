# AGENTS.md

## Purpose
This directory holds AthenaK runtime input decks (plain-text `.athinput` files)
used for examples, problem setups, and regression tests. Files here are consumed
by `ParameterInput` and passed to the executable via `-i` in `src/main.cpp`.

---

## Directory Map (Current)
- `dyngr/`: dynamical GRMHD input decks.
- `grhydro/`: GR hydrodynamics input decks.
- `grmhd/`: GRMHD input decks.
- `hydro/`: Newtonian hydro input decks.
- `ion-neutral/`: two-fluid ion-neutral input decks.
- `mhd/`: Newtonian MHD input decks.
- `particles/`: particle-module input decks.
- `radiation/`: radiation input decks.
- `shearing_box/`: shearing-box input decks.
- `srhydro/`: special-relativistic hydro input decks.
- `srmhd/`: special-relativistic MHD input decks.
- `tests/`: regression-test input decks referenced by the test harness.
- `z4c/`: Z4c numerical-relativity input decks.
- `turb_timed_amr_stage1.athinput`: top-level standalone input deck.

---

## Input File Format (ParameterInput)
Format rules are enforced by `src/parameter_input.cpp`:

- **Blocks**: `<blockname>` must appear on its own line (after leading spaces).
  The parser reads until the first `>`; missing `>` is fatal.
- **Parameters**: `name = value` (spaces around `=` are optional). Each parameter
  must be on a single line. The parser assumes `=` is present.
- **Comments**: `#` starts a comment. Lines whose first non-space character is
  `#` are ignored. Inline comments after `#` are stripped from the value.
- **Whitespace**: tabs are stripped; blank/whitespace-only lines are ignored.
- **Ordering**: any `name=value` line before the first block is a fatal error.
- **Duplicates**: repeated block names are merged; repeated parameter names in a
  block overwrite earlier values.
- **Terminator**: `<par_end>` stops parsing early (used when reading restarts).
- **Case**: block/parameter names are case-sensitive string matches.

---

## Parameter Access Semantics
- All values are stored as strings and converted on demand by `Get*` functions.
- `GetInteger/GetReal/GetBoolean/GetString` require the block and parameter to
  exist; missing entries are fatal.
- `GetOrAdd*` inserts a default if the parameter is missing.
- Boolean parsing accepts `0/1` or `true/false` (case-insensitive).
- There is no global schema validation; extra parameters are stored and ignored
  unless a module explicitly reads them.

---

## Command-Line Overrides and Restarts
- `src/main.cpp` loads input files via `-i <file>` and restarts via `-r <file>`.
- If both `-r` and `-i` are provided, the input file is loaded *after* the
  restart header and overrides its parameters.
- `ParameterInput::ModifyFromCmdline` accepts overrides of the form
  `block/param=value` and applies them after file loading. The block and
  parameter must already exist or the code aborts.

---

## Tests and Tooling
- The test harness (`tst/scripts/utils/athena.py`) constructs paths as
  `inputs/<filename>`, so tests assume files live under this directory. Renaming
  or moving input decks requires updating the test scripts.

### PR2 PIC Coupling Regression Decks
- Main coupled deck:
  - `inputs/tests/pic_mhd_current_coupling.athinput`
- Passive-background isolation deck:
  - `inputs/tests/pic_mhd_passive_mode.athinput`
- No-MHD Boris isolation deck:
  - `inputs/tests/pic_no_mhd_boris.athinput`
- Midpoint E+B Boris verification deck:
  - `inputs/tests/pic_boris_midpoint_eb.athinput`
  - includes diagnostics outputs for `prtcl_dpxdt/dpydt/dpzdt`, `prtcl_dedt`,
    and `prtcl_ebdot` to validate midpoint-frozen-in and feedback bookkeeping.
- Guard decks used for fatal-path validation:
  - `inputs/tests/pic_mhd_coupling_guard_radiation.athinput`
  - `inputs/tests/pic_mhd_coupling_guard_nr.athinput`
  - `inputs/tests/pic_mhd_coupling_guard_hydro.athinput`
  - `inputs/tests/pic_mhd_coupling_guard_isothermal.athinput`
  - `inputs/tests/pic_mhd_coupling_guard_edge_relativistic.athinput`
- PR2 coupling-sensitive overrides commonly used by tests:
  - `particles/couple_moments_to_mhd`
  - `particles/couple_j_to_efield_coeff`
  - `particles/couple_j_to_efield_representation`
  - `particles/couple_moments_momentum_to_mhd`
  - `particles/couple_moments_energy_to_mhd`
  - `particles/couple_fluid_feedback_order`
  - deterministic CR init knobs: `particles/cr_vx0`, `particles/cr_vy0`,
    `particles/cr_vz0`.
  - optional per-species drift overrides:
    `speciesN/vx0`, `speciesN/vy0`, `speciesN/vz0` (fallback to global
    `particles/cr_v*0` when omitted).
- Staged PIC runtime-control overrides (parse/guard stage):
  - `particles/pic_background_mode = coupled|passive_mhd|no_mhd`
  - `particles/pic_feedback_mode = coupled|test_particle`
  - `particles/pic_interp_scheme = tsc`
  - `particles/pic_cr_light_speed`, `particles/pic_max_cell_cross`,
    `particles/pic_theta_max`
  - `particles/pic_deltaf_mode = off|on` (with `particles/pic_deltaf_f0`)
  - `particles/pic_sort_interval`
  - `particles/pic_intermediate_arrays = auto|off`
  - `particles/pic_expanding_box_mode = off|on`
  - `particles/pic_expansion_rate_x1/x2/x3` (require expanding-box mode on)
- Passive-mode guard overrides used by regressions:
  - `particles/pic_background_mode = coupled|passive_mhd|no_mhd`
  - `particles/pic_feedback_mode = test_particle`
  - `particles/couple_moments_to_mhd` (rejected in passive mode)
- No-MHD guard overrides used by regressions:
  - `particles/pic_background_mode = no_mhd|coupled`
  - `particles/pic_feedback_mode = test_particle|coupled`
  - `particles/pic_no_mhd_bx`, `particles/pic_no_mhd_by`,
    `particles/pic_no_mhd_bz`

### Entity-Mirroring PIC Regression Decks
- `inputs/tests/pic_entity_deposit_mink.athinput`
  - Pair-species (`q=-1,+1`) periodic deposit case used for decomposition
    invariance checks across MPI partitionings.
- `inputs/tests/pic_entity_deposit_reflect.athinput`
  - Pair-species reflect-boundary deposit case used for non-periodic
    decomposition/neutrality checks.
- `inputs/tests/pic_em_vacuum_wave.athinput`
  - EM-vacuum-style convergence anchor using `linear_wave` in strict
    test-particle mode for current AthenaK architecture.
- `inputs/tests/pic_langmuir_frequency_proxy.athinput`
  - Frequency-accuracy proxy in `pic_background_mode=no_mhd` using a uniform
    `Bz` Boris orbit; regression extracts dominant frequency from deposited
    `prtcl_jx/prtcl_jy` time series.
- `inputs/tests/pic_two_stream_growth_proxy.athinput`
  - Counter-streaming no-MHD proxy with two species and opposite `vx0` drifts;
    regression fits the dominant density-mode growth rate and gates out
    positive exponential growth.

---

## Extension Points
- Add new input decks for problem generators or regression tests in the
  appropriate subdirectory.
- Keep new parameter names aligned with the modules that read them (use
  `GetOrAdd*` for optional values and `Get*` for required ones).
