<!-- BEGIN build-memory-table -->
# Boundary Values Navigation

## Scope

This subtree owns ghost-zone exchange for cell- and face-centered mesh fields, flux
correction across refinement boundaries, physical boundary kernels, and particle
migration. Mesh topology and neighbor metadata come from `src/mesh/`; physics task lists
decide when these operations run.

## Task lookup

| Task | Start in | Also inspect or validate |
| --- | --- | --- |
| Change boundary flags, buffer types, or the public interface | `src/bvals/bvals.hpp`, `src/bvals/bvals.cpp` | `src/mesh/mesh.hpp`, all physics callers |
| Change cell-centered exchange | `src/bvals/bvals_cc.cpp`, `src/bvals/buffs_cc.cpp` | Hydro, radiation, Z4c task paths |
| Change face-centered exchange | `src/bvals/bvals_fc.cpp`, `src/bvals/buffs_fc.cpp` | `src/mhd/mhd_tasks.cpp`, constrained transport |
| Change AMR flux correction | `src/bvals/flux_correct_cc.cpp`, `src/bvals/flux_correct_fc.cpp` | Fluid update ordering and AMR convergence tests |
| Change boundary prolongation | `src/bvals/prolongation.cpp`, `src/bvals/prolong_prims.cpp` | `src/mesh/prolongation.hpp`, `src/eos/AGENTS.md` |
| Change a physical boundary condition | `src/bvals/physics/` | Matching fluid state layout, pgen user-boundary hooks |
| Change particle crossing/migration | `src/bvals/bvals_part.cpp` | `src/particles/particles_tasks.cpp`, particle tests |
| Change task adapters or cleanup ordering | `src/bvals/bvals_tasks.cpp` and module `*_tasks.cpp` | `src/tasklist/task_list.hpp` |

## Important flow

- `MeshBoundaryValues::InitializeBuffers` derives fixed neighbor-buffer ranges from mesh
  dimensionality and AMR level relationships. The index ordering must match
  `MeshBlock::nghbr` and `src/mesh/nghbr_index.hpp`.
- A fluid stage normally initializes receives, computes and corrects fluxes, updates the
  state, restricts and exchanges evolved variables, applies physical boundaries,
  prolongates coarse/fine boundaries, then converts boundary conserved state to
  primitives.
- `MeshBoundaryValuesCC` handles conserved/radiation/Z4c arrays;
  `MeshBoundaryValuesFC` separately handles magnetic face fields and edge-field flux
  correction. Particle boundary values use message lists rather than mesh-field buffers.

## Local validation

- `cd tst && python3 run_test_suite.py --test test_suite/rad/test_rad_lwave1d_amr_cpu.py`:
  cell-centered AMR boundary exchange.
- `cd tst && python3 run_test_suite.py --test test_suite/sr/test_sr_lwave2d_amr_mpicpu.py`:
  MPI coarse/fine exchange for hydro and face-centered MHD state.
- `cd tst && python3 run_test_suite.py --test test_suite/z4c/test_z4c_lwave2d_amr_mpicpu.py`:
  Z4c's higher-order boundary path.

## Local constraints

- Do not reorder the fixed neighbor/buffer indices without updating both mesh neighbor
  construction and every send/receive mapping.
- Preserve the split between variable buffers and flux buffers and between same-, coarse-,
  and fine-level ranges.
- Physical boundary kernels must honor the exact variable parity and face staggering of
  their physics module. User boundaries are enrolled through `ProblemGenerator`.
- Z4c communicates extra same-level coarse data for higher-order interpolation; do not
  assume its CC layout is identical to fluid CC exchange.
<!-- END build-memory-table -->
