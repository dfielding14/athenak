<!-- BEGIN build-memory-table -->
# Particle Science Deck Navigation

## Scope

This directory contains hand-maintained particle example and science decks,
centered on turbulent MHD-PIC boxes. Particle implementation lives in `src/`;
small regression carriers live in `inputs/tests/`; Bell and shock campaign
freezes live in `inputs/publication/`.

## Task lookup

| Task | Start in | Also inspect or validate |
| --- | --- | --- |
| Change the matched three-dimensional turbulent-box matrix | `pic_turbulent_dynamo_128_*` | Keep MHD, test-particle, weak-feedback, and fiducial coupled controls synchronized |
| Change the near-isothermal restart-injection path | `pic_turbulent_dynamo_128_*_gamma1p00001.athinput` | Restart injection guards, output/restart code, and the precursor/CR pair |
| Change the controlled two-dimensional comparison | `pic_turbulent_dynamo_512_2d_*` | `docs/2D_MHD_PIC_DYNAMO_COMPARISON_README.md` and both MHD/CR variants |
| Change forcing calibration | `pic_turbulent_dynamo_calibration_64.athinput` | `src/srcterms/turb_driver.*` and production forcing parameters |
| Change the legacy drift example | `random_particle_drift.athinput` | Particle drift initialization and pusher behavior |

## Important flow

The turbulent-box decks require a build with `-DPROBLEM=turb`. Their initial
conditions and history callback come from `src/pgen/turb.cpp`; magnetic seeds
come from `src/srcterms/initial_perturbations.*`; forcing comes from
`src/srcterms/turb_driver.*`; CR modes, six-beam loading, feedback, and restart
injection come from `src/particles/` and restart/output code.

The shared initialization path has a focused regression:

- `cd tst && python3 run_tests.py particles/pic_turbulent_dynamo_smoke --cmake=-DPROBLEM=turb`

That smoke test uses `inputs/tests/pic_turbulent_dynamo_smoke.athinput`; it
checks initialization, beam cancellation, one driven MHD/PIC energy budget,
and 2D3V/3-D forcing spectra. It does not establish production runtime scale
or scientific convergence.

## Local constraints

- Keep matched MHD/CR decks aligned in mesh, EOS, forcing, perturbation
  spectrum/seeds, and output cadence unless a difference is intentional and
  documented.
- Restart injection requires an explicitly enabled CR deck and a compatible
  particle-free MHD precursor. Changes must preserve restart fingerprints and
  avoid reapplying new-run initialization.
- Preserve opposite six-beam momentum/current cancellation when changing
  species definitions or loading.
- Describe the 512-square setup as a controlled two-dimensional turbulent
  MHD-PIC comparison, not a true three-dimensional dynamo.
<!-- END build-memory-table -->
