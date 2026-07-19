<!-- BEGIN build-memory-table -->
# Outputs and Restart Navigation

## Scope

This subtree owns parsing `<outputN>` blocks, mapping runtime fields to output variables,
derived diagnostics, all file-format writers, and low-level serial/MPI I/O. Restart reads
are completed in `src/main.cpp`, `src/mesh/build_tree.cpp`, and `src/pgen/pgen.cpp`, so a
restart format change crosses those boundaries.

## Task lookup

| Task | Start in | Also inspect or validate |
| --- | --- | --- |
| Add an output format or change `<outputN>` parsing | `src/outputs/outputs.cpp`, `outputs.hpp` | `Driver` output scheduling, `src/CMakeLists.txt` |
| Add or remap an output variable | `src/outputs/outputs.hpp`, `basetype_output.cpp` | `derived_variables.cpp`, owning physics arrays |
| Change a derived diagnostic | `src/outputs/derived_variables.cpp` | Ghost-zone/stencil requirements and matching tests |
| Change history or event reductions | `src/outputs/history.cpp`, `eventlog.cpp` | Pgen history hooks, source-term/frame-tracker data |
| Change mesh, particle, PDF, or projection formats | Matching writer in `src/outputs/` | Readers in `vis/python/` and analysis scripts |
| Change restart serialization | `src/outputs/restart.cpp` | `src/pgen/pgen.cpp`, `src/mesh/build_tree.cpp`, `src/main.cpp` |
| Change shared/per-rank/per-node I/O | `src/file_sharding.hpp`, `src/outputs/io_wrapper.*` | Format writers, `docs/file_sharding_explained.md` |

## Important flow

1. `Outputs` scans input blocks and constructs a `BaseTypeOutput` subclass for each
   configured format; restart output is kept last because it serializes updated output
   counters and parameters.
2. `Driver` checks time/cycle cadence, calls `LoadOutputData`, then calls the format's
   `WriteOutputFile`.
3. `BaseTypeOutput` selects MeshBlocks and variables, computes requested derived fields,
   and stages device data to host arrays. Each writer owns its final layout and MPI
   coordination.
4. Restart output serializes global metadata plus physics payloads. The restart path first
   reconstructs the mesh, then `ProblemGenerator` restores module arrays and reenrolls
   problem-specific callbacks.

## Local validation

- `cd tst && python3 run_test_suite.py --test test_suite/nr/test_nr_frame_tracking_restart_cpu.py`:
  restart continuation and history/output state on CPU.
- `cd tst && python3 run_tests.py restart/per_node_restart --cmake=-DAthena_ENABLE_MPI=ON`:
  legacy MPI regression for shared, per-rank, and per-node restart layouts.
- When changing a data format, run a consuming test or script using the corresponding
  reader under `vis/python/`; file existence alone does not validate layout compatibility.

## Local constraints

- Adding a selectable variable requires keeping `NOUTPUT_CHOICES`, `var_choice`, variable
  validation/loading, labels, and any derived computation in sync.
- Adding a format requires a `BaseTypeOutput` subclass, constructor dispatch in
  `outputs.cpp`, and source registration in `src/CMakeLists.txt`.
- Treat restart layout changes as schema changes: update writer, every read/shard path,
  size/offset accounting, and compatibility tests together.
- Shared, per-rank, and per-node modes change coordination and paths, not ownership of
  physics data. Preserve zero-count participation where a collective MPI-IO path requires
  it.
<!-- END build-memory-table -->
