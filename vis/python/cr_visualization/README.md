# AthenaK MHD and Particle Visualization

These readers pair AthenaK binary MHD output with particle-major merged track
files. They are path- and resolution-independent. The Pm1 production data below
is a worked example, not a hard-coded dataset. Use the
`feature/CR_tracers_followup_architecture` branch.

## Environment

On Frontier, use the same Python stack used to create the merged tracks. No
virtual environment is required.

```bash
module restore
module load cpe/24.07 cray-mpich/8.1.30 cray-python/3.11.7

REPO=/ccs/home/dfielding/athenak-cr-tracers-followup-architecture
cd "$REPO"
export PYTHONPATH="$REPO/vis/python${PYTHONPATH:+:$PYTHONPATH}"
```

The code uses `numpy`, `h5py`, and, for the distributed example, `mpi4py`.

## One MeshBlock and Its Particle Tracks

Choose one MHD rank file and a MeshBlock index within it. The script reads that
MeshBlock with `bin_convert.read_single_rank_binary_as_athdf()`, scans the
corresponding raw track shard for particle identities recorded inside its
physical bounds, and reads only those particle rows from the merged HDF5 file.
It does not scan the full merged track payload.

To inspect the MeshBlocks present in an unfamiliar rank file:

```python
import bin_convert

raw = bin_convert.read_binary_as_athdf(mhd_bin, raw=True)
print(raw["n_mbs"])
print(raw["mb_logical"])   # logical i, j, k, refinement level
print(raw["mb_geometry"])  # x1, x2, and x3 lower/upper bounds
```

For the Pm1 run, the particles were seeded from restart `00024` at `t=6`. The
matching frozen MHD snapshot is `full_mhd_w_bcc.00024.bin`:

```bash
MHD_ROOT=/lustre/orion/ast207/proj-shared/dfielding/AMR/data/Pm1_4096_eta3e-6/bin
RUN=/lustre/orion/ast207/proj-shared/dfielding/AMR/particles/campaign-pm1-4096-eta3e-6-static19-richtrk-v1/Pm1_4096_eta3e-6_static19_richtrk_ppc1e-4_t6_to_t106
BASE=Pm1_4096_eta3e-6_static19_richtrk_ppc1e-4_t6_to_t106

python3 vis/python/cr_visualization/read_meshblock_and_tracks.py \
  --mhd-bin "$MHD_ROOT/rank_00000000/Pm1_S4_eta3e-6.full_mhd_w_bcc.00024.bin" \
  --raw-trk "$RUN/trk/rank_00000000/$BASE.trk" \
  --merged-tracks "$RUN/$BASE.tracks.h5" \
  --meshblock 0 \
  --max-particles 32 \
  --time-stride 20 \
  --output "$RUN/Pm1_rank0_mb0_tracks32_stride20.h5"
```

The optional output is a self-contained HDF5 bundle:

```text
/mhd/{x1f,x1v,x2f,x2v,x3f,x3v}
/mhd/{MHD field names}
/mhd/bounds
/tracks/particles
/tracks/{times,cycles,values}
/tracks/inside_meshblock
```

The bundled track array is ordered `[particle, time, field]`. Its `fields`
attribute names the final dimension. `inside_meshblock` is a Boolean array with
shape `[particle, time]`. Add `--inside-only` to replace samples outside the
chosen MeshBlock with `NaN`, which breaks plotted lines at the block boundary.

The same operation can be used directly from Python:

```python
from cr_visualization.cr_data import (
    find_track_visits,
    read_meshblock,
    read_merged_track_subset,
)

mhd = read_meshblock(mhd_bin, meshblock_index=0)
visits = find_track_visits(raw_track_files, mhd["Bounds"])
tracks = read_merged_track_subset(
    merged_tracks,
    visits.output_tags,
    bounds=mhd["Bounds"],
    max_particles=32,
    time_stride=20,
)
```

`mhd["velx"]`, `mhd["vely"]`, and `mhd["velz"]` are the fluid velocity.
`mhd["bcc1"]`, `mhd["bcc2"]`, and `mhd["bcc3"]` are the cell-centered magnetic
field. Every MHD field is ordered `[z, y, x]`; the matching coordinates are
`x3v`, `x2v`, and `x1v`. Particle `vx`, `vy`, and `vz` are particle velocity,
not fluid velocity.

The merged `particles` table supplies `output_tag`, `track_tag`, and `species`.
For the logarithmic CR species used by the Pm suite, mass is reconstructed from
`particles/min_mass` and `particles/mass_log_spacing` in the run's
runtime-overrides file as `min_mass * mass_log_spacing**species`. A configured
field-line tracer species is massless and should be handled separately.

The raw shard provides spatial ownership records, not complete histories. A
particle can cross many ranks, so the script uses the shard only to discover
`output_tag` values and obtains complete trajectories from the merged HDF5.
Pass multiple files after `--raw-trk` when tracks were written per node, through
a shared file, or with a decomposition that does not map one shard to the
selected MHD rank. Positions are always checked against the actual MeshBlock
bounds.

Useful selection options are:

```text
--quantities dens velx vely velz bcc1 bcc2 bcc3
--track-fields x y z bx by bz k1 k2 k3
--species 0 1 2
--time-min 6 --time-max 16
--max-particles 32 --seed 1
--time-stride 20
```

If `--quantities` or `--track-fields` is omitted, all available fields are
returned. Legacy ASCII-framed and compact-header track shards are both
supported. A raw track format may contain additional fields; the spatial query
requires only `tag` (or `output_tag`), `time`, `x`, `y`, and `z`.

## Entire MHD Snapshot and All Particle Tracks

For a small simulation, the complete MHD snapshot can be assembled serially:

```python
import bin_convert

mhd = bin_convert.read_all_ranks_binary_as_athdf(
    rank0_bin,
    quantities=["velx", "vely", "velz", "bcc1", "bcc2", "bcc3"],
)
```

Do not use that serial path for a 4096 cubed snapshot. The Pm1 MHD output is
about 2 TiB, and the merged rich tracks are 533 GiB. For production data, each
MPI process reads a share of the MHD rank files with
`bin_convert.read_binary_as_athdf(raw=True)` and a contiguous share of particle
rows from the merged HDF5. No data is gathered onto MPI rank zero.

```bash
srun -N 16 -n 128 \
  python3 vis/python/cr_visualization/read_full_dataset_mpi.py \
  --mhd-rank0 "$MHD_ROOT/rank_00000000/Pm1_S4_eta3e-6.full_mhd_w_bcc.00024.bin" \
  --merged-tracks "$RUN/$BASE.tracks.h5" \
  --quantities velx vely velz bcc1 bcc2 bcc3 \
  --time-stride 20 \
  --particle-batch 4
```

Collectively, `local_meshblocks` covers every rank file and `local_tracks`
covers every particle row. The distributed example deliberately leaves these
objects local:

```python
meshblocks, tracks = read_local_data(args, MPI.COMM_WORLD)
```

Each entry of `meshblocks` contains one MeshBlock, its coordinates, logical
location, bounds, time, cycle, and requested fields. `tracks` contains the
local contiguous particle range, times, cycles, field names, and values. Dave
and Ken can pass these local objects directly into their distributed
visualization pipeline.

At full fidelity the Pm1 track array has shape `(38912, 200000, 16)`. Use
`--time-stride 1` when every saved point is required. For interactive rendering,
a stride of 10 or 20 retains all particles while reducing trajectory geometry
substantially. Long time ranges can also be processed as consecutive windows
with `--time-min` and `--time-max` instead of residing in memory at once.
`--particle-batch` bounds temporary HDF5 read memory independently of the local
particle partition size.

For an evolving MHD calculation, run the MHD reader once per desired binary
snapshot and select the corresponding particle time interval. The Pm pushing
runs use `frozen_mhd=true`, so their single restart-matched MHD snapshot applies
to the full particle trajectory.
