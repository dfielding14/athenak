# Final release four-rank CPU passive acceptance

The six permanent `test_cgl_passive_mpicpu.py` groups pass on four CPU MPI ranks:
80 application launches and 65 acceptance checks. This is additional CPU-MPI
coverage; the earlier CPU1/HIP1/HIP4 result counts remain unchanged.

The exact final relinked CPU binary is recorded by path and SHA-256 in
`manifest.json`. The runner uses CPU-only Cray modules,
`MPICH_GPU_SUPPORT_ENABLED=0`, four MPI ranks and zero GPUs. The original CPU
release binary with unintended ROCm link dependencies was not used. Every
application's working directory, temporary directory and outputs reside under
WO2. These concurrent functional runs support no timing claims.

Coverage includes all five released reconstructors, exact native-isothermal
flow and timestep identity, forcing in three dimensions, floor/weak-field
conditions, independent periodic heating, linear response, thermal advection,
full-state resumed restart identity, and explicit unsupported-mode fences.
The package ran its unchanged assertions on the frozen `source-release` tests.

`pytest.log` records all six passing groups. `groups/` contains the complete
per-group result records and commands; the full-precision restart artifacts
remain at the paths in those records. `manifest.json` includes source and
runner hashes, binary identity, loaded runtime contract and completed counts.
