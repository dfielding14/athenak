# CR Rich TRK Recovery Implementation Guide

This guide fixes the current production-blocking `trk` regression: the active
CR branch writes only position and velocity, while the original PK production
tracker wrote the magnetic-field and geometry diagnostics needed for analysis.

Ground-truth source for the old format is
`pkempski/athenak-PK:particlesPK`, specifically:

```text
src/outputs/outputs.hpp
src/outputs/track_prtcl.cpp
```

The old buffered production record writes 18 floats per tracked particle:

```text
tag,time,x,y,z,vx,vy,vz,Bx,By,Bz,K1,K2,K3,dB1,dB2,dB3,jmag
```

Do not infer this list from memory.  Preserve this order unless every consumer
is updated deliberately.

## Objectives

1. Restore the rich PK `trk` payload.
2. Support three file layouts for `file_type = trk`:

```text
single_file_per_rank = true,  single_file_per_node = false  -> trk/rank_XXXXXXXX/<basename>.trk
single_file_per_rank = false, single_file_per_node = true   -> trk/node_XXXXXXXX/<basename>.trk
single_file_per_rank = false, single_file_per_node = false  -> trk/<basename>.trk
```

`single_file_per_rank` and `single_file_per_node` must be mutually exclusive.
There should be no `track_single_file_per_rank` compatibility alias.
Any input deck that will set these through command-line overrides must still
declare both keys in the relevant `output` block, because AthenaK rejects
overrides for parameters that are absent from the parsed deck.

3. Reuse the sharding design from
`origin/feature/single-file-per-node-outputs` rather than inventing a parallel
API.  In that branch, `src/file_sharding.hpp` defines
`FileShardMode::{shared, per_node, per_rank}` and helpers such as
`ShardDirectoryName`, `GatherShardCounts`, `PrefixCountBeforeMe`, and
`ShardCommunicator`.

4. Make the shared single-file mode use real MPI-IO.  This mode is useful for
small tests and compatibility, but it should not be the production default at
large scale.

5. Run a complete replacement PM suite after validation, because the existing
CR evolution runs did not write the required tracking fields.

## Implementation Work Order

### 1. Port output sharding support

Bring the relevant infrastructure from
`origin/feature/single-file-per-node-outputs` onto
`feature/CR_tracers_followup_architecture`:

```text
src/file_sharding.hpp
globals: node_id, rank_in_node, ranks_per_node, nnodes, rank_to_node, node_comm
main.cpp: MPI_Comm_split_type(MPI_COMM_TYPE_SHARED) initialization and cleanup
outputs.cpp: parse single_file_per_rank/single_file_per_node into FileShardMode
outputs.hpp: OutputParameters::file_shard_mode
IOWrapper support for setting the shard communicator if needed
```

Keep this port surgical.  Do not wholesale replace unrelated output code from
the node-output branch.

The parser behavior should match the node-output branch:

```text
single_file_per_rank=true and single_file_per_node=true -> fatal error
single_file_per_rank=true -> FileShardMode::per_rank
single_file_per_node=true -> FileShardMode::per_node
neither set -> FileShardMode::shared
```

Wire `file_type = trk` through this parser.

### 2. Restore the rich tracked-particle data structure

Extend `TrackedParticleData` in `src/outputs/outputs.hpp` from:

```text
tag,x,y,z,vx,vy,vz
```

to the PK fields:

```text
tag,x,y,z,vx,vy,vz,Bx,By,Bz,K1,K2,K3,dB1,dB2,dB3,jmag
```

Do not add `time` to the in-memory struct unless it simplifies buffering.  The
time value is known at write time and is part of the on-disk record.

### 3. Restore field population in `track_prtcl.cpp`

For each tracked particle:

```text
Bx,By,Bz = pr(IPBX), pr(IPBY), pr(IPBZ)
K1,K2,K3 = (bhat . grad) bhat
dB1,dB2,dB3 = grad |B|
jmag = |curl B|
```

Use the PK branch as the starting formula source.  Then harden it:

- validate that MHD fields exist before using `pm->pmb_pack->pmhd`;
- clamp or safely reject derivative indices that would leave the available
  ghost-zone range;
- guard against `|B| = 0` before normalizing `bhat`;
- preserve current `track_per_species` tag mapping and duplicate/missing-tag
  checks.

The current branch's tracked-tag ordering and fail-fast checks are useful; keep
them.

### 4. Define a self-describing rich `trk` frame

Each frame should have a text header followed by binary `float` payload:

```text
# AthenaK tracked particle data at time= ... nranks= ... cycle= ...
# trk_format=rich_v1 nfields=18 record_count=...
# fields=tag,time,x,y,z,vx,vy,vz,bx,by,bz,k1,k2,k3,db1,db2,db3,jmag
```

Then write `record_count * 18` floats.

Every file layout should write the same record format.  This is essential:
rank/node/shared readers should differ only in file discovery, not in record
schema.

### 5. Implement the three `trk` layouts

#### Per-rank

Use for production by default.

```text
output901/file_type = trk
output901/single_file_per_rank = true
output901/single_file_per_node = false
```

Each rank writes only local tracked records to:

```text
trk/rank_XXXXXXXX/<basename>.trk
```

No global gather is needed.  The record contains explicit `tag` and `time`, so
files can be merged or read independently.

#### Per-node

Use when file count is a problem but shared MPI-IO is too risky.

```text
output901/single_file_per_rank = false
output901/single_file_per_node = true
```

Within each `node_comm`, gather local rich records to the node leader and write:

```text
trk/node_XXXXXXXX/<basename>.trk
```

This should mirror the node-output branch's `FileShardMode::per_node` model.
It avoids one file per rank while avoiding a full-run rank-0 funnel.

#### Shared single file

Use only for small tests or compatibility.

```text
output901/single_file_per_rank = false
output901/single_file_per_node = false
```

All ranks write to:

```text
trk/<basename>.trk
```

This must use MPI-IO, not rank-0 gather.  For each frame:

1. Compute local record count.
2. Gather shard counts or compute an `MPI_Exscan` prefix.
3. Have rank 0 append the frame header and determine payload offset.
4. Use collective MPI-IO writes so each rank writes its local packed records at
   `payload_offset + 18*sizeof(float)*prefix`.

This mode will likely fail or perform poorly at full Frontier scale, so do not
make it the production default.

### 6. Update readers and analyzers

Update `scripts/analyze_cr_pusher_accuracy.py` and any campaign track readers
to parse:

```text
legacy 6-float current files: x,y,z,vx,vy,vz
PK buffered 18-float files: tag,time,x,y,z,vx,vy,vz,Bx,By,Bz,K1,K2,K3,dB1,dB2,dB3,jmag
new rich_v1 files with explicit nfields and fields headers
```

The accuracy tests only need `x,y,z,vx,vy,vz`, but the parser must not silently
misread rich frames.  Tests should assert that rich fields are finite for CR
production-like runs.

### 7. Use subagents deliberately

Use subagents for parallel checking, not for making unsupervised design changes.
Recommended split:

```text
agent-1: compare PK particlesPK tracker against current tracker and produce an exact field/order diff
agent-2: inspect origin/feature/single-file-per-node-outputs and summarize minimal sharding code to port
agent-3: after implementation, review trk frame parser/writer consistency and binary offsets
agent-4: after build, inspect smoke outputs and validate headers, field counts, file counts, finite values
agent-5: before production restart, audit submit scripts and athinputs for correct trk settings
```

The main agent should own final edits, resolve conflicts, run builds/tests, and
make the go/no-go decision.

## Validation Ladder

### Local/static checks

Run:

```bash
git diff -- src/outputs/track_prtcl.cpp src/outputs/outputs.hpp src/outputs/outputs.cpp src/file_sharding.hpp src/main.cpp
cmake --build build-frontier --target athena -j 8
```

Check:

```text
no remaining output901/track_single_file_per_rank
trk uses single_file_per_rank/single_file_per_node
tracked record field count is 18 in rich mode
shared mode does not use rank-0 gather for payload
```

### Tiny GPU smoke

Run tiny tracked-particle cases for all three layouts:

```text
single_file_per_rank=true,  single_file_per_node=false
single_file_per_rank=false, single_file_per_node=true
single_file_per_rank=false, single_file_per_node=false
```

For each, verify:

```text
expected file layout exists
headers contain trk_format=rich_v1 and nfields=18
record_count matches payload size
tag/time columns are present
Bx/By/Bz/K1/K2/K3/dB1/dB2/dB3/jmag are finite
existing position/velocity accuracy checks still pass
```

### Production-like smoke

Use the current frozen restart machinery in:

```text
/lustre/orion/ast207/proj-shared/dfielding/AMR/particles
```

Submit a short run from the same restart style as the PM suite with:

```text
output901/file_type = trk
output901/dt = 0.0005
output901/nparticles = 2048
output901/ncycle = 2000
output901/buffer_size = 40000000
output901/single_file_per_rank = true
output901/single_file_per_node = false
```

Verify file volume and cadence:

```text
2048 tracks/species * 19 species = 38912 records/frame
38912 records * 18 floats * 4 bytes ~= 2.80 MB/frame
dt=0.0005 -> ~= 5.6 GB per simulation time unit
100 L/c -> ~= 560 GB total trk payload before headers/filesystem overhead
```

This is much larger than the current broken 6-field output, but it is the
required analysis product.

## Production Rerun

Once the production-like smoke passes, the existing CR evolution products
should be treated as incomplete for track-based science.  Launch a replacement
PM suite from the appropriate latest particle restarts or fluid restarts.

Required production output settings:

```text
trk rich fields enabled by default
single_file_per_rank = true
single_file_per_node = false
shared single-file mode disabled for production
psamp disabled unless explicitly needed
pmom dt = 0.5
pspec dt = 0.5
prst dt = 5.0
time/ndiag = 1000
particles/subcycle_per_particle_gyro = true
particles/subcycle_gyro_fraction = 0.05
particles/exchange_gyro_only_substeps = false
particles/update_global_counts_each_exchange = false
```

Rerun all PM values in the suite, monitor early throughput, and cancel/revise
if rich `trk` output dominates runtime more than expected.

Minimum acceptance before leaving jobs unattended:

```text
one fresh rich trk frame exists per PM run
reader reports nfields=18
sampled records have finite B/K/dB/jmag
restart cadence is correct
performance is within the new rich-output expectation
no old 6-field trk-only production jobs remain queued/running as if valid
```
