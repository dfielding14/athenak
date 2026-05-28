# `single_file_per_node` Review Notes And Follow-Up Plan

This branch now implements `single_file_per_node` as a general file-sharding
mode for restart, binary, coarsened-binary, spherical-slice, and PDF outputs.
The design is sound: it reduces file count from one file per MPI rank to about
one file per physical node, while keeping the heaviest coordination within a
node-local MPI communicator.

## Current Status

- Restart output uses a shared manifest plus one payload shard per node under
  `rst/node_XXXXXXXX/`.
- Restart input accepts either the manifest path or a path inside a rank/node
  shard directory and normalizes the request to the correct layout.
- Large MPI-IO transfers use byte-wise chunking so code paths are not limited
  by `int` MPI counts.
- Collective MPI-IO paths keep ranks with zero local payload in the collective
  call, which avoids deadlocks and incomplete collectives.
- Binary, coarsened-binary, spherical-slice, and PDF readers can reassemble
  rank- or node-sharded files from any shard path.
- Regression tests cover per-node restart round trips, restart path variants,
  stale payload collisions, forced small MPI-IO chunks, binary/coarsened-binary
  empty-shard reassembly, spherical-slice equivalence, and N-dimensional PDF
  weighting compatibility.

## Remaining Performance Work

These are useful follow-up optimizations, but they are separable from the
correctness and usability fixes in this review.

1. Stream binary and coarsened-binary MeshBlock payloads directly to their final
   offsets instead of staging `write_mbs * data_size` bytes in one host buffer.
   This will reduce peak host memory for very large shards.
2. Replace the spherical-slice per-node gather with a sparse ownership exchange
   if `ntheta*nphi*nvars` becomes large enough that root-side node staging is a
   measurable bottleneck.
3. Add an optional lightweight benchmark script that records shard counts,
   payload bytes, open counts, and wall time for shared, per-rank, and per-node
   modes on the same input.
4. Add a compatibility fixture for a committed legacy per-node restart file if
   the project decides to preserve older experimental files permanently.

## Documentation Integration

The publishable documentation page is
[`docs/source/modules/output_file_sharding.md`](source/modules/output_file_sharding.md).
It is intentionally independent of this planning note and can be copied into the
GitHub Pages branch without bringing along review-only TODOs.

Suggested `gh-pages` integration:

1. Add `modules/output_file_sharding` to the hidden modules toctree in
   `docs/source/index.md`.
2. Add an entry in `docs/source/modules/index.md` under Support Systems.
3. Add a short link from `docs/source/modules/outputs.md`, preferably near the
   `single_file_per_rank` output-parameter discussion.
