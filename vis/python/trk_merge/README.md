# Rich Track Merge

`merge_rich_trk.py` converts AthenaK `trk_format=rich_v1` tracked-particle shards into one particle-major HDF5 file. It is intentionally specific to the rich format with magnetic-field, curvature, gradient, and current-magnitude fields; older or leaner track formats need a separate schema-aware merger. The raw `trk` files are written by the rank or node that owns each particle at each output time, so one particle's history can be spread across many files. The merger parses each frame, redistributes records by `output_tag` with MPI, sorts each owned trajectory by cycle, and assembles `/values[particle, time, field]` along with `/particles`, `/times`, and `/cycles`.

Example workflow:

```bash
module load cray-python/3.11.7 cray-mpich

srun -N 16 -n 128 python3 vis/python/trk_merge/merge_rich_trk.py \
  --run-dir /path/to/Pm1_4096_eta3e-6_static19_richtrk_ppc1e-4_t6_to_t106 \
  --layout rank \
  --tmp-dir /path/to/scratch/merge_tmp/Pm1 \
  --output /path/to/Pm1_4096_eta3e-6_static19_richtrk_ppc1e-4_t6_to_t106.tracks.h5
```

For large production runs, use `--layout rank`. `--tmp-dir` is scratch space for the MPI merge: each rank writes temporary owner-sorted binary shards there before rank 0 assembles the final HDF5. Put it on Lustre, not a login-node local filesystem, and budget roughly one extra copy of the final `/values` payload while the merge is running. The directory is deleted after a successful merge unless `--keep-temp` is set. The output HDF5 stores rich fields in this order:

```text
x,y,z,vx,vy,vz,bx,by,bz,k1,k2,k3,db1,db2,db3,jmag
```
