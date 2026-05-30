> **HISTORICAL ONLY.** Do not execute commands or infer qualification,
> publication, branch, or Frontier-submission authority from this file. Use
> `tst/publication/PIC_PRODUCTION_READINESS_PLAN.md` as the sole controlling plan.

# AthenaK MHD and PIC Test Problem Catalog (Internal)

This is a concise working catalog of currently implemented regression test
problems for:
- AthenaK baseline MHD/hydro physics tests,
- Entity-derived (or Entity-mirroring) PIC tests, and
- Athena++ MHD-PIC paper-inspired PIC/MHD-PIC tests.

Primary use:
- track what already exists,
- identify which cases are strongest for publication-quality plots,
- guide expansion toward richer CR-shock MHD-PIC results.

## 1) Inventory Summary

- Baseline hydro/MHD tests in harness: 3 scripts
  - `tst/scripts/hydro/hydro_linwave.py`
  - `tst/scripts/mhd/mhd_linwave.py`
  - `tst/scripts/mhd/divb_amr.py`
- Particle/PIC tests in harness: 25 scripts under `tst/scripts/particles/`
  (including utility self-test `pic_analysis_utils.py`).

## 2) AthenaK Baseline MHD/Hydro Tests

| Test ID | Script | Input deck(s) | Core check | Plot value |
| --- | --- | --- | --- | --- |
| `hydro.hydro_linwave` | `tst/scripts/hydro/hydro_linwave.py` | `inputs/tests/linear_wave_hydro.athinput` | L1 error thresholds + convergence across integrator/reconstruction/Riemann choices | Convergence plots (L1 vs resolution) for sound/entropy waves |
| `mhd.mhd_linwave` | `tst/scripts/mhd/mhd_linwave.py` | `inputs/tests/linear_wave_mhd.athinput` | L1 error thresholds + convergence for fast/Alfven/slow/entropy waves | Convergence plots for MHD wave families |
| `mhd.divb_amr` | `tst/scripts/mhd/divb_amr.py` | `inputs/mhd/ffc_divb.athinput` | `max(|divB|) < 1e-10` on AMR output | `divB` map / histogram for AMR cleanliness |

## 3) Entity-Derived / Entity-Mirroring PIC Tests

| Test ID | Script | Input deck(s) | Physical target | Current gate | Plot value |
| --- | --- | --- | --- | --- | --- |
| `particles.pic_entity_deposit_mink` | `tst/scripts/particles/pic_entity_deposit_mink.py` | `inputs/tests/pic_entity_deposit_mink.athinput` | Minkowski deposit parity | Exact/near-exact `Q,J,npart` and field parity vs baseline | Deposit stencil parity figures |
| `particles.pic_entity_deposit_reflect` | `tst/scripts/particles/pic_entity_deposit_reflect.py` | `inputs/tests/pic_entity_deposit_reflect.athinput` | Reflecting-boundary deposit parity | Exact/near-exact `Q,J` checks and decomposition consistency | Boundary deposit parity figures |
| `particles.pic_em_vacuum_wave` | `tst/scripts/particles/pic_em_vacuum_wave.py` | `inputs/tests/pic_em_vacuum_wave.athinput` | EM vacuum propagation accuracy | L1 error convergence and serial/MPI agreement | Error-vs-resolution EM wave panel |
| `particles.pic_langmuir_frequency_proxy` | `tst/scripts/particles/pic_langmuir_frequency_proxy.py` | `inputs/tests/pic_langmuir_frequency_proxy.athinput` | Langmuir oscillation frequency | Dominant frequency close to expected value | Frequency spectrum + fitted dominant mode |
| `particles.pic_two_stream_growth_proxy` | `tst/scripts/particles/pic_two_stream_growth_proxy.py` | `inputs/tests/pic_two_stream_growth_proxy.athinput` | Two-stream mode growth behavior | Growth-rate envelope and serial/MPI consistency | Mode amplitude vs time (log scale) |
| `particles.pic_weibel_growth_proxy` | `tst/scripts/particles/pic_weibel_growth_proxy.py` | `inputs/tests/pic_weibel_growth_proxy.athinput` | Weibel mode growth behavior | Growth-rate envelope and serial/MPI consistency | Current/field mode growth panel |

## 4) Athena++ MHD-PIC Paper-Inspired Tests

| Test ID | Script | Input deck(s) | Paper-aligned target | Current gate | Plot value |
| --- | --- | --- | --- | --- | --- |
| `particles.pic_no_mhd_boris` | `tst/scripts/particles/pic_no_mhd_boris.py` | `inputs/tests/pic_no_mhd_boris.athinput` | Test-particle/no-MHD orbit path | Invariants + guard checks | Orbit/momentum phase portrait |
| `particles.pic_boris_midpoint_eb` | `tst/scripts/particles/pic_boris_midpoint_eb.py` | `inputs/tests/pic_boris_midpoint_eb.athinput` | Section-2 midpoint E+B Boris behavior | Frozen-in and momentum/energy exchange consistency | Midpoint-Boris validation panel |
| `particles.pic_bell_growth_proxy` | `tst/scripts/particles/pic_bell_growth_proxy.py` | `inputs/tests/pic_bell_growth_proxy.athinput` | Bell instability growth | Coupled growth lower bound, uncoupled suppression | `B` growth-rate comparison |
| `particles.pic_multispecies_backreaction_oscillation` | `tst/scripts/particles/pic_multispecies_backreaction_oscillation.py` | `inputs/tests/pic_multispecies_osc_uniform.athinput`, `inputs/tests/pic_multispecies_osc_smr.athinput`, `inputs/tests/pic_multispecies_osc_amr_proxy.athinput` | Electron/positron-gas oscillation parity | Frequency window, energy drift, SMR/AMR consistency | Frequency + energy drift vs mesh strategy |
| `particles.pic_crsi_deltaf_proxy` | `tst/scripts/particles/pic_crsi_deltaf_proxy.py` | `inputs/tests/pic_crsi_deltaf_proxy.athinput` | CRSI branch growth with `delta f` | Dominant branch growth/r2/noise checks | Polarization branch growth and noise reduction |
| `particles.pic_crpai_polarization_proxy` | `tst/scripts/particles/pic_crpai_polarization_proxy.py` | `inputs/tests/pic_crpai_prolate_proxy.athinput`, `inputs/tests/pic_crpai_oblate_proxy.athinput` | CRPAI polarization behavior | Dominant branch/r2 and split checks | Left/right polarization spectra |
| `particles.pic_expanding_box_anisotropy_proxy` | `tst/scripts/particles/pic_expanding_box_anisotropy_proxy.py` | `inputs/tests/pic_expanding_box_proxy.athinput`, `inputs/tests/pic_compressing_box_proxy.athinput` | Expanding/compressing anisotropy trend | Opposite slope sign and separation checks | Anisotropy ratio trend vs time |
| `particles.pic_refinement_boundary_characterization` | `tst/scripts/particles/pic_refinement_boundary_characterization.py` | `inputs/tests/pic_refinement_boundary_smr_proxy.athinput`, `inputs/tests/pic_refinement_boundary_amr_proxy.athinput` | Coarse/fine boundary behavior | Drift/CV/jump envelopes + cross-decomp checks | Boundary homogeneity metrics |
| `particles.pic_amr_shock_lb_smoke` | `tst/scripts/particles/pic_amr_shock_lb_smoke.py` | `inputs/tests/pic_amr_shock_lb_smoke.athinput` | Reduced CR shock + AMR/LB smoke | Shock signatures, tail growth, AMR creation, MPI parity | Shock profile + AMR/LB diagnostics |

## 5) AthenaK MHD-PIC Integration and Safety Tests

| Test ID | Script | Input deck(s) | Focus |
| --- | --- | --- | --- |
| `particles.pic_analysis_utils` | `tst/scripts/particles/pic_analysis_utils.py` | n/a | Utility fit-function self-check |
| `particles.pic_deposit_conservation` | `tst/scripts/particles/pic_deposit_conservation.py` | `inputs/tests/pic_deposit_conservation.athinput` | Deposit conservation + parser/guard negatives |
| `particles.pic_decomp_invariance` | `tst/scripts/particles/pic_decomp_invariance.py` | `inputs/tests/pic_deposit_conservation.athinput` | Decomposition invariance for deposited moments |
| `particles.pic_mhd_passive_mode` | `tst/scripts/particles/pic_mhd_passive_mode.py` | `inputs/tests/pic_mhd_passive_mode.athinput` | Passive MHD mode behavior + guards |
| `particles.pic_mhd_current_coupling` | `tst/scripts/particles/pic_mhd_current_coupling.py` | `inputs/tests/pic_mhd_current_coupling.athinput`, `inputs/tests/pic_mhd_current_coupling_default_mode.athinput` | Main coupling correctness, direct-vs-convert parity, continuity, guards |
| `particles.pic_mhd_coupling_decomp` | `tst/scripts/particles/pic_mhd_coupling_decomp.py` | `inputs/tests/pic_mhd_current_coupling.athinput` | Coupled decomposition invariance + continuity |
| `particles.pic_mhd_coupling_multilevel` | `tst/scripts/particles/pic_mhd_coupling_multilevel.py` | `inputs/tests/pic_mhd_coupling_multilevel.athinput` | Multilevel coupling parity and response checks |
| `particles.pic_mhd_coupling_nonperiodic` | `tst/scripts/particles/pic_mhd_coupling_nonperiodic.py` | `inputs/tests/pic_mhd_coupling_nonperiodic.athinput`, `inputs/tests/pic_mhd_coupling_nonperiodic_default_mode.athinput` | Nonperiodic coupling behavior + BC guard checks |
| `particles.pic_mhd_restart_fidelity` | `tst/scripts/particles/pic_mhd_restart_fidelity.py` | `inputs/tests/pic_mhd_restart_fidelity.athinput` | Restart full-vs-segment consistency, direct-vs-convert ratios |
| `particles.pic_restart_safety_guards` | `tst/scripts/particles/pic_restart_safety_guards.py` | `inputs/tests/pic_restart_safety_guards.athinput` | Runtime safety guards + per-rank restart watch |

## 6) Publication-Oriented Priority Set

Recommended first-pass figure candidates (highest value):

1. Entity PIC core validation
- `pic_entity_deposit_mink`, `pic_entity_deposit_reflect`,
  `pic_em_vacuum_wave`, `pic_langmuir_frequency_proxy`,
  `pic_two_stream_growth_proxy`, `pic_weibel_growth_proxy`.
- Goal: establish solver correctness and dispersion/growth credibility.

2. Athena++ benchmark reproductions
- `pic_bell_growth_proxy`, `pic_multispecies_backreaction_oscillation`,
  `pic_crsi_deltaf_proxy`, `pic_crpai_polarization_proxy`,
  `pic_expanding_box_anisotropy_proxy`.
- Goal: show known CR-driven instability and anisotropy trends with
  decomposition robustness.

3. CR-shock MHD-PIC storyline (highest-impact target)
- Start from `pic_amr_shock_lb_smoke`.
- Current status: strong smoke/regression gate, not yet a full science case.
- Needed to become publication-grade:
  - longer run duration and larger dynamic range in shock downstream,
  - richer outputs for CR distribution diagnostics (`f(x,p)`, spectra),
  - controlled parameter scans (Mach number, CR loading, resolution),
  - multi-panel plots: density/B/current structure, CR phase-space, spectra,
    and AMR/load-balance context.

## 7) Immediate Next Plotting Batch (Practical)

For near-term internal figure generation, prioritize:

1. `pic_bell_growth_proxy`
- Plot: `log(|B|)` vs time for coupled vs uncoupled and MPI parity.

2. `pic_crsi_deltaf_proxy` and `pic_crpai_polarization_proxy`
- Plot: left/right polarization branch amplitudes with fitted growth rates.

3. `pic_multispecies_backreaction_oscillation`
- Plot: oscillation frequency extraction and energy-drift comparison
  (uniform vs SMR vs AMR-proxy).

4. `pic_amr_shock_lb_smoke`
- Plot: shock profile snapshots (`rho`, `|B|`, `|J|`), tail metric evolution,
  and AMR-created block history as a bridge to the full CR-shock science run.

## 8) Publication Tooling and Campaign Assets

Publication workflow utilities now live in:
- `tst/publication/pic_publication_manifest.py`
- `tst/publication/run_pic_publication_suite.py`
- `tst/publication/plot_pic_publication_figures.py`
- `tst/publication/run_pic_shock_scan.py`
- `tst/publication/pvtk_particles.py`

New shock publication decks:
- `inputs/tests/pic_amr_shock_lb_publication_local.athinput`
- `inputs/tests/pic_amr_shock_lb_publication_hpc.athinput`

These decks are scan-ready with explicit Orszag-Tang controls in `<problem>`:
- `ot_mach` (Mach-proxy, maps to velocity amplitude),
- `ot_v0`, `ot_p0`, `ot_d0`, `ot_B0` (direct IC controls when needed).

Particle phase-space support for publication plots:
- `pvtk` particle outputs now include velocity vectors (`VECTORS vel`) and
  explicit integer scalar headers (`gid`, `ptag`, `species`/`sn_id`), enabling
  direct speed-spectrum and phase-space post-processing from `.part.vtk`.

## 9) Publication Upgrade Execution Plan

Detailed step-by-step execution plan for moving all F01-F09 figures to
publication readiness is documented in:
- `AGENT_PIC_PUBLICATION_READINESS_PLAN.md`

How it fits this catalog:
1. Uses Section 3 and Section 4 tests as the unchanged regression baseline.
2. Uses Section 6 and Section 7 as the publication figure target set.
3. Uses Section 8 tooling (`manifest`, run harness, plotting scripts) for all
   staged runs and figure exports.
4. Adds publication-physics variants without replacing current proxy tests.
