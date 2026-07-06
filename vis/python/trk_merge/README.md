# Rich Track Merge

`merge_rich_trk.py` converts AthenaK rich tracked-particle shards into one particle-major HDF5 file. It supports both legacy ASCII-framed `trk_format=rich_v1` inputs and compact-header `trk_format=rich_v2` inputs; both carry the same rich magnetic-field, curvature, gradient, and current-magnitude payload. Older or leaner track formats need a separate schema-aware merger. The raw `trk` files are written by the rank or node that owns each particle at each output time, so one particle's history can be spread across many files. The merger parses each frame, redistributes records by `output_tag` with MPI, sorts each owned trajectory by cycle, and assembles `/values[particle, time, field]` along with `/particles`, `/times`, and `/cycles`.

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

## Plot One Particle

`plot_random_particle.py` is a small example script for inspecting one trajectory from a merged `.tracks.h5` file. It picks `--row` if supplied, otherwise a reproducible random particle from `--seed`, then plots the magnetic moment proxy `mu_M = v_perp^2 / 2B` and the curvature scaled by `2 pi c / Omega`.

```bash
python3 vis/python/trk_merge/plot_random_particle.py \
  /path/to/Pm1_4096_eta3e-6_static19_richtrk_ppc1e-4_t6_to_t106.tracks.h5 \
  --row 1024 \
  --tmax 5 \
  --history-file /lustre/orion/ast207/proj-shared/dfielding/AMR/data/Pm1_4096_eta1e-6/Pm1_S4_eta1e-6.user.hst \
  --out particle_track.png
```

The script needs `B_rms` to normalize `2 pi c / Omega`. Prefer `--history-file`, pointing at the matching MHD `.hst` file; the script reads the nearest-time `B^2` history column and uses `B_rms = sqrt(B^2)`. You can override this with `--b-rms`. If neither is supplied, the script loudly warns and uses `B_rms = 1`, which is only useful for a quick shape check.

The merged HDF5 stores particle `species`, but not the actual particle mass. For these production runs the mass is inferred from the run directory recorded in the HDF5 metadata, specifically the nearby `*.runtime_overrides.txt` file with `particles/min_mass` and `particles/mass_log_spacing`. If that file is absent, pass those two values explicitly with `--min-mass` and `--mass-log-spacing`.
