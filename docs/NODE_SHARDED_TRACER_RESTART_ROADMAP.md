# Node-Sharded Tracer Restart Roadmap

## Current boundary

Lagrangian Monte Carlo tracer restarts are rank-local: each rank appends its own
particle header, schedule state, real arrays, and integer arrays to its restart
file. MPI runs therefore require `single_file_per_rank=true`. Node-sharded
tracer restarts, plus shared files in MPI, fail before restart output or input
begins; serial shared files remain supported.

## Smallest viable extension

Keep the existing node payload and transaction machinery. Append one tracer
section to each node payload instead of adding another file family. Version the
node manifest and add only the data needed to locate and validate those sections:

- global tracer count, particle dimensions, schedule state, and `next_tracer_tag`;
- per-node tracer count, byte offset, and byte length;
- checked totals proving that the node inventories cover the global tracer set.

Continue reading version-1 manifests as mesh-only restarts. Rank-sharded restart
files and their particle section stay unchanged.

## Implementation path

1. **Write node-local sections.** Gather particle counts within each node, compute
   rank prefix offsets, and have ranks write their existing host mirrors directly
   into the node payload with collective MPI-IO. Do not gather all particles onto
   the node leader. Store schedule state once and verify that all ranks agree.
2. **Publish atomically.** Close and validate every payload before publishing the
   manifest, using the existing temporary-file and reservation flow. Any particle
   write failure must discard the incomplete generation.
3. **Read independently of the saved rank count.** Validate all tracer inventory
   records before allocation. Divide the saved particle records across the current
   ranks, read bounded contiguous slices, then reuse the existing particle
   ownership/redistribution path to move each tracer to the rank owning its mesh
   location.
4. **Restore global state.** Broadcast the saved schedules and tag counter, reject
   inconsistent dimensions or duplicate schedule identifiers, and update local and
   global particle counts after redistribution.
5. **Remove the guard.** Enable `single_file_per_node=true` for Lagrangian Monte
   Carlo restarts only after both write and read paths pass the checks below.

## Required checks

- zero tracers globally, empty ranks, and nodes with no local tracers;
- one and multiple ranks per node, plus restart with a different MPI rank count;
- AMR ownership changes across restart and tracers on block boundaries;
- exact schedule, tag, real-array, and integer-array round trips;
- truncated sections, bad offsets/counts, missing payloads, and interrupted writes;
- CPU MPI tests followed by multi-node HIP runs on the target AMD system, including
  restart time, peak host memory, and bytes written per node.

The first implementation should target the existing particle record layout and
same executable precision. A portable self-describing particle format can wait
until cross-build restart compatibility is an actual requirement.
