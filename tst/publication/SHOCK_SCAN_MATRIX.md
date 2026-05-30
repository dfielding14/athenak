# Default Orszag-Tang Shock-Rich Engineering Stress Scan Matrix

This is the default grid emitted by `run_pic_shock_scan.py` when no custom
lists are provided.

This is an engineering proxy and stress workflow. It is not a Sun and Bai
(2023) Section 5.4 parallel-shock reproduction matrix, a production
qualification campaign, or Frontier submission authorization. Use
`PIC_PRODUCTION_READINESS_PLAN.md` for qualification requirements and execution
policy.

## Local Defaults

- Mach proxy (`problem/ot_mach`): `0.8, 1.0, 1.2`
- CR loading (`particles/ppc`): `2, 4, 8`
- Resolution (`mesh/nx1 = mesh/nx2`): `32, 48, 64`
  - with `mesh/nx3 = max(8, nx1/4)`.

Total combinations: `3 x 3 x 3 = 27`.

## HPC Emit-Only Defaults

- Mach proxy (`problem/ot_mach`): `1.0, 1.2, 1.4`
- CR loading (`particles/ppc`): `8, 16, 24`
- Resolution (`mesh/nx1 = mesh/nx2`): `256, 384, 512`
  - with `mesh/nx3 = max(8, nx1/4)`.

Total combinations: `3 x 3 x 3 = 27`.

## Intended Use

1. Local staging:
- run a reduced subset (for example `--max-cases 3`) to confirm stability,
  diagnostics, and figure pipeline plumbing.

2. Approved external engineering-stress campaign:
- emit an intentionally inert planning script with
  `--deck hpc --nproc <cluster_ranks>`,
- select the matching matrix with `--tier hpc`,
- translate selected commands into immutable per-submission snapshots and
  execute only through the canonical validated wrapper under an explicitly
  approved resource policy,
- archive outputs per basename for post-processing.

## Output Naming

Each run uses a deterministic basename:
- `pic_shock_scan_<tier>_m<MACH>_ppc<PPC>_n<NX>`

This enables straightforward grouping by scan axis during figure generation.
Generated shell files are intentionally inert planning artifacts at both tiers.
Use the Python helper with explicit local acknowledgement for workstation runs.
