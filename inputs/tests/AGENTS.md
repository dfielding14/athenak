<!-- BEGIN build-memory-table -->
# Regression Input Navigation

## Scope

This directory contains hand-maintained runtime fixtures for AthenaK's compiled
regression and standalone validation tests. The decks do not implement physics;
they select code paths in `src/` and are normally owned together with the test
module that references them.

## Task lookup

| Task | Start with | Also inspect or validate |
| --- | --- | --- |
| Change a classical hydro, MHD, or reconstruction fixture | The matching deck and `tst/scripts/hydro/` or `tst/scripts/mhd/` module | Problem generator, EOS/reconstruction path, and focused harness run |
| Change mesh, AMR, boundary, or load-balance coverage | `*amr*`, `*smr*`, boundary, and restart decks | `src/mesh/`, `src/bvals/`, top-level `tst/test_*load_balance.py`, and decomposition variants |
| Change radiation, GR, or Z4c coverage | The named deck | Matching suite under `tst/scripts/` and source module guide |
| Change a particle pusher or timestep fixture | `pic_*gyro*`, `pic_boris_*`, or `pic_no_mhd_*` | `tst/scripts/particles/` and `src/particles/particles_pushers.cpp` |
| Change deposition or MHD-PIC coupling | `pic_*deposit*`, `pic_mhd_*`, or `pic_paper_*` | Particle test guide, particle moments/tasks, MHD tasks, and conservation checks |
| Change particle restart, AMR, or physical-boundary behavior | `pic_*restart*`, `pic_*refinement*`, or boundary decks | Particle migration/boundary code and serial/MPI variants |
| Change parallel-shock behavior | `pic_parallel_shock_*` | `src/pgen/tests/pic_parallel_shock.cpp`, matching particle tests, and any publication deck sibling |
| Change a bounded PIC proxy | Bell, CRSI, CRPAI, multispecies, expanding-box, or turbulence deck | Its analyzer/test and evidence classification |

## Test-to-deck flow

- `tst/scripts/utils/athena.py` resolves test input names below `inputs/`.
  Trace `_INPUT*` constants and command-line overrides in the owning test before
  renaming or changing a deck.
- Pairing is not one-to-one: one test may use several decks, and several tests
  may share one deck. Update every consumer when a fixture contract changes.
- Run focused compiled tests from `tst/`, for example
  `python3 run_tests.py particles/pic_relativistic_gyro_timestep`.
- Preparation, proxy, candidate, and publication-like names do not establish
  scientific qualification or campaign authority.

Current focused anchors include
`pic_relativistic_gyro_timestep.athinput`,
`pic_parallel_shock_split_removal_restart.athinput`, and
`pic_parallel_shock_mignone_tracer_smoke.athinput`; keep each synchronized with
the same-named module under `tst/scripts/particles/`.

## Generated Q043 oracle decks

`q043_bell_current_volume_aware_deposited_current_oracle/` is a generated,
manifest-bound matrix. Do not hand-edit individual `.athinput` files. Regenerate
and validate the complete matrix with
`tst/publication/q043_bell_current_volume_aware_deposited_current_oracle.py`,
then update `deck_manifest.json` and its host test together.
<!-- END build-memory-table -->
