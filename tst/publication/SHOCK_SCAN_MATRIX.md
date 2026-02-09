# Default CR-Shock Scan Matrix (Publication)

This is the default grid emitted by `run_pic_shock_scan.py` when no custom
lists are provided.

## Axes

- Mach proxy (`problem/ot_mach`): `0.8, 1.0, 1.2`
- CR loading (`particles/ppc`): `2, 4, 8`
- Resolution (`mesh/nx1 = mesh/nx2`): `64, 96, 128`
  - with `mesh/nx3 = max(8, nx1/4)`.

Total combinations: `3 x 3 x 3 = 27`.

## Intended Use

1. Local staging:
- run a reduced subset (for example `--max-cases 3`) to confirm stability,
  diagnostics, and figure pipeline plumbing.

2. Cluster campaign:
- emit full command script with `--deck hpc --nproc <cluster_ranks>`,
- execute on GPU-capable resources,
- archive outputs per basename for post-processing.

## Output Naming

Each run uses a deterministic basename:
- `pic_shock_scan_<deck_tag>_m<MACH>_ppc<PPC>_n<NX>`

This enables straightforward grouping by scan axis during figure generation.

