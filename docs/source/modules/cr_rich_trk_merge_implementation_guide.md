# CR Rich TRK Merge Implementation Guide

This guide specifies the post-processing tool for converting high-cadence
rank-sharded rich `trk` output into one particle-major HDF5 file per PM run.

The production `trk` output is intentionally sharded by rank to avoid large-scale
MPI-IO failures:

```text
trk/rank_XXXXXXXX/<basename>.trk
```

Those files are ownership shards, not trajectories.  A particle can move between
MPI ranks, so following one particle requires reading all rank shards, grouping
records by particle identity, and sorting each particle by time/cycle.

## Current Production Scale

The `t=6` to `t=106` rich-track PM suite uses:

```text
ntrack_per_species = 2048
nspecies = 19
track_per_species = 1
ntracked_prtcls = 38912
dt_trk = 0.0005
duration = 100 L/c
expected samples ~= 200000
nfields = 18
```

The completed PM1 and PM3 `trk/` directories are about `1.1 TB` each as raw
rank-sharded output.  This is larger than the scientific payload because every
rank writes frame headers at high cadence, including many frames with few or no
local tracked particles.

The full science-preserving merged payload is not the 2048-particle estimate.
It is closer to:

```text
38912 particles * 200000 samples * 18 float32 fields ~= 560 GB
```

The planned dense HDF5 layout factors out redundant `tag` and `time`, so the
main values dataset is:

```text
38912 particles * 200000 samples * 16 float32 values ~= 498 GB
```

This is still large, but it is far easier to analyze and should be about a
factor of two smaller than the rank-sharded directories.

## Tool Location

Implement the merger as:

```text
scripts/merge_rich_trk_hdf5.py
```

The implementation language is Python.  Use `mpi4py` for production-scale
parallelism, `numpy` for parsing/sorting/packing, and `h5py` for the final HDF5
file.  Do not require `numba` in the first implementation; the main bottlenecks
are filesystem I/O, record redistribution, sorting, and large contiguous writes.
Add optional `numba` only if profiling shows a specific loop dominates.

On the current Frontier environment, `mpi4py` is available but `h5py` is not
MPI-enabled:

```text
mpi4py: available
h5py.get_config().mpi: false
```

Therefore the first implementation must not depend on parallel HDF5.  Use MPI
for input parsing and redistribution, then do a serial final HDF5 assembly from
sorted temporary shards.

## Output HDF5 Format

Write one HDF5 file per PM run, for example:

```text
Pm1_4096_eta3e-6_static19_richtrk_ppc1e-4_t6_to_t106.tracks.h5
```

Use this layout:

```text
/
  attrs:
    format = "athenak_rich_trk_merged_v1"
    source_trk_format = "rich_v1"
    source_run_dir
    source_commit
    source_layout = "rank"
    source_nfiles
    source_bytes
    nfields_source = 18
    nvalue_fields = 16
    ntracked_prtcls
    ntrack_per_species
    nspecies
    track_per_species
    merge_command
    merge_mpi_size

/particles
  compound table with one row per output particle:
    output_tag int64
    species int32
    track_tag int64
    row int64
    count int64
    first_time float64
    last_time float64

/times
  float64 [ntimes]

/cycles
  int64 [ntimes]

/values
  float32 [nparticles, ntimes, 16]
  fields:
    x,y,z,vx,vy,vz,bx,by,bz,k1,k2,k3,db1,db2,db3,jmag
```

The original rich record fields are recoverable as:

```text
tag,time,x,y,z,vx,vy,vz,bx,by,bz,k1,k2,k3,db1,db2,db3,jmag
```

where `tag` comes from `/particles/output_tag`, `time` comes from `/times`, and
the remaining 16 fields come from `/values[row, :, :]`.

`/times` must be derived from the payload `time` field, not from the printed
header time.  The writer stores payload time as `float32` and prints header time
with limited stream precision, so they can differ slightly.  Use header `cycle`
as the authoritative ordering/duplicate key.

Default HDF5 choices:

```text
compression = none
shuffle = false
values chunks = (1, 32768, 16)
times chunks = contiguous or one chunk
particles chunks = contiguous or one chunk
```

Compression should be an optional command-line setting, not the default.  The
first priority is predictable throughput and simple random access to one
particle trajectory.

## Command-Line Interface

Initial CLI:

```bash
python scripts/merge_rich_trk_hdf5.py \
  --run-dir /path/to/Pm1_..._t6_to_t106 \
  --output /path/to/Pm1_...tracks.h5 \
  --tmp-dir /path/to/scratch/merge_pm1 \
  --layout rank \
  --require-complete \
  --delete-temp
```

MPI production invocation:

```bash
srun -N <nodes> -n <ranks> python scripts/merge_rich_trk_hdf5.py ...
```

Useful options:

```text
--layout auto|rank|node|shared
--require-complete
--allow-missing
--keep-cycle
--compression none|lzf|gzip
--max-files-per-rank N
--exchange-rows N
--dry-run
--validate-only
--delete-temp
```

`--require-complete` is the default for production.  If every particle does not
have one record for every global frame, the tool should fail and report the
first missing/duplicate `(output_tag, cycle)` pairs.  `--allow-missing` can be
added later to write NaN-filled gaps plus a `/present` mask, but that is not the
first production path.

## Parallel Algorithm

### 1. Build a source manifest

Rank zero enumerates source files:

```text
rank layout:   trk/rank_*/*.trk
node layout:   trk/node_*/*.trk
shared layout: trk/*.trk
```

For each file, store path and byte size.  Broadcast the manifest to all MPI
ranks.  Assign files round-robin or by greedy byte balance so each rank gets a
similar number of input bytes.

### 2. Parse frames in parallel

Each MPI rank reads its assigned files sequentially.  The parser must:

- find `# AthenaK tracked particle data at time= ...`,
- parse `cycle`, `ntracked_prtcls`, `ntrack_per_species`,
  `track_per_species`, and `record_count`,
- require `trk_format=rich_v1`,
- require `nfields=18`,
- require the exact field list
  `tag,time,x,y,z,vx,vy,vz,bx,by,bz,k1,k2,k3,db1,db2,db3,jmag`,
- read `record_count * 18` little-endian `float32` values,
- reject malformed/truncated payloads.

The parser should use `numpy.frombuffer` or `numpy.fromfile` and reshape to
`(record_count, 18)`.  Avoid Python objects per record.

### 3. Redistribute records by particle identity

After the first valid header determines `ntracked_prtcls`, assign contiguous tag
ranges to MPI ranks:

```text
owner_rank = floor(output_tag * mpi_size / ntracked_prtcls)
```

Contiguous ownership makes the final particle table naturally sorted by
`output_tag`.

Process input in synchronized batches.  `MPI_Alltoallv` is collective, so ranks
must not flush opportunistically at different row counts.  For the production
rank layout, the natural synchronized batch is one source rank file per MPI
rank:

1. Convert `records[:, 0]` to integer `output_tag`.
2. Compute `owner_rank`.
3. Group rows by owner.
4. Use `MPI_Alltoallv` on packed `MPI.BYTE` buffers containing the structured
   temp records.
5. Each owner appends received rows to one local unsorted temp file with a
   fixed little-endian structured schema:

```text
tag:int64, cycle:int64, time:float64, values:float32[16]
```

The source payload stores `time` but not `cycle`; `cycle` must be carried from
the frame header because it is the authoritative ordering key.

```text
tmp/owner_000123.unsorted.bin
```

This avoids writing one temporary file per bucket per rank and avoids keeping
the full run in memory.  The first implementation may ignore `--exchange-rows`
inside a single source file; rank-layout production files are small enough that
file-boundary exchange is the safer simple implementation.

### 4. Sort owner shards

After all records are redistributed, each owner rank reads or memory maps its
local structured unsorted file.

Sort by:

```text
output_tag, cycle
```

The dense output does not need to store per-record cycle because `/cycles`
stores the global frame cycle array.  It also does not need to store per-record
`tag` or `time`; those are reconstructed from `/particles` and `/times`.

For each owned particle:

1. Check tags are contiguous and equal to the rank's assigned tag interval.
2. Check no duplicate times/cycles for the same tag.
3. Check the count equals `ntimes` under `--require-complete`.
4. Strip `tag,cycle,time` and write the 16 value fields to:

```text
tmp/owner_000123.values.f32
tmp/owner_000123.index.npy
```

The local index should include `output_tag`, `species`, `track_tag`, local row
offset, count, first time, and last time.

### 5. Assemble the final HDF5 file

Because current `h5py` is serial, rank zero assembles the final file after a
barrier:

1. Read all owner indexes.
2. Validate that particle rows cover `output_tag = 0..ntracked_prtcls-1`.
3. Create `/particles`, `/times`, `/cycles`, and `/values`.
4. For each owner shard, copy its particle-major values into the correct
   `/values[row0:row1, :, :]` hyperslab.
5. Flush and close the HDF5 file.
6. If `--delete-temp` is set, remove temp shards only after the HDF5 file has
   passed read-back validation.

Serial assembly is acceptable for the first version because the expensive part
is the 8192-file, roughly terabyte-scale parse.  The final write is a large
contiguous HDF5 write of roughly 500 GB.

## Validation Requirements

The merger must fail fast on format mismatches:

```text
missing rich_v1 header
nfields != 18
field list differs from the rich PK-compatible order
inconsistent ntracked_prtcls across files
inconsistent ntrack_per_species across files
truncated payload
non-finite values in required diagnostics
```

The merger must report, but not necessarily fail, if the source contains empty
rank frames.  Empty local frames are normal.  The global frame count should still
sum to `ntracked_prtcls`.

Production acceptance checks:

```text
all cycles/times are globally consistent
every global frame has exactly ntracked_prtcls records
every output_tag appears exactly once per frame
every particle has count == ntimes
particle times are strictly increasing
HDF5 read-back of several random particles matches source records bitwise for the 16 value fields
```

Read-back spot checks should include:

```text
first particle
last particle
one particle per species
several random output_tag values
particles known to migrate between ranks, if easy to identify
```

## Test Plan

### Unit/parser tests

Use existing tiny rich `trk` layout smoke outputs.  Test all supported source
layouts:

```text
shared
rank
node
```

Each test should parse headers, payload sizes, field names, and frame counts.

### Small MPI integration test

Run with 4 to 8 MPI ranks on a small rich-track directory:

```bash
srun -n 4 python scripts/merge_rich_trk_hdf5.py \
  --run-dir <tiny-rich-trk-run> \
  --output <tmp>/tracks.h5 \
  --tmp-dir <tmp>/merge \
  --require-complete
```

Compare merged HDF5 trajectories against the existing analyzer's merged frames.

### Production dry run

Before writing the full HDF5 file for a PM run:

```bash
srun -N 1 -n 8 python scripts/merge_rich_trk_hdf5.py \
  --run-dir <PM-run> \
  --tmp-dir <scratch>/merge_pm1 \
  --dry-run \
  --validate-only
```

This should report:

```text
source_nfiles = 8192
source_bytes ~= 1.1 TB for PM1/PM3
ntracked_prtcls = 38912
ntrack_per_species = 2048
nspecies = 19
ntimes ~= 200000
estimated_values_bytes ~= 498 GB
```

### Full production merge

Run one PM case first.  Use enough MPI ranks that each owner shard is only a few
GB:

```text
64 ranks  -> owner payload roughly 8 to 9 GB
128 ranks -> owner payload roughly 4 to 5 GB
256 ranks -> owner payload roughly 2 to 3 GB
```

Prefer 128 ranks for the first full merge unless queue constraints make that
awkward.  After one PM case passes validation, run the rest of the suite.

## Reader Convenience API

Add a small helper module or examples to read trajectories:

```python
import h5py

with h5py.File("Pm1.tracks.h5", "r") as f:
    particles = f["particles"]
    times = f["times"][:]
    values = f["values"]

    row = particles["row"][particles["output_tag"][:] == output_tag][0]
    trajectory = values[row, :, :]
```

Field lookup should come from `f["values"].attrs["fields"]`, not hard-coded in
analysis notebooks.

## Work Order

1. Implement the strict rich `trk` frame parser.
2. Implement manifest generation and byte-balanced MPI file assignment.
3. Implement MPI `Alltoallv` redistribution by contiguous output-tag ownership.
4. Implement owner-rank sorting, per-particle validation, and compact temp shard
   writing.
5. Implement serial HDF5 assembly and read-back validation.
6. Test on tiny rich layout outputs.
7. Test on a short PM directory subset.
8. Run one full PM merge.
9. Run the full PM suite merge.
10. Add a short usage note to `docs/source/modules/particles.md` after the tool
    is working.
