# Final release equivalence runner

Preparation only; no release comparison applications have been launched yet.

After root confirms the rebuilt binaries are immutable, from WO2 run:

```bash
SLURM_JOB_ID=5629018 bash final-research/run_release_equivalence.sh hip \
  --output /lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/final-research/release-equivalence-hip-5629018
SLURM_JOB_ID=5629018 bash final-research/run_release_equivalence.sh cpu \
  --output /lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/final-research/release-equivalence-cpu-5629018
```

Defaults compare the accepted `final-research/bin/athena-{cpu,hip}` (verified
against the existing immutable manifest) to
`final-research/release-build-{cpu,hip}/src/athena`. `--new-binary` can select a
different completed release binary. The two paths must differ and stay under
WO2. All binary hashes are checked before launches and after the complete run.
No existing binary, input or result is overwritten.

Each backend runs the exact 23 saved Task7 configurations and two additional
forced passive 3D safe/full cases at one/four ranks. The latter use the accepted
Task7 passive fixture with 12 cycles and the tested Task4 forcing block, changing seed 519 to 1
so the existing seed-1-only dormant startup audit remains applicable;
fixture provenance is in `release-equivalence-fixtures/manifest.json`.
`--list-cases` prepares no applications and lists names, ranks and hashes.
`--cases NAME ...` chooses exact case names; unknown/empty selections fail.
`--node NODE` can constrain steps to an allocated node if root requests it.

This is **final-state equivalence**, with one untraced accepted/release pair per
configuration. It provides no new face-array coverage and no performance claim.
Every binary field and history file must match byte for byte. All normalized
restart bytes must match, including diagnostics and live forcing/RNG state;
there is no q-work or other roundoff allowance. Raw hashes are retained.

The existing 36-byte unused root-index normalization applies to all restarts.
Only forced fixtures explicitly opt in to the separately proven startup-dormant
RNG fields and four alignment-padding bytes, after the helper validates the
native v3 metadata/signatures. Each file records exact normalized ranges and
raw hashes. No state bytes or live RNG members are ignored.

The HIP wrapper enforces the complete accepted runtime contract after loading
modules: GPU-aware MPI, managed-memory flag, HSA_XNACK, GPU NIC policy, IPC cache,
and FI settings. The Python runner checks these values before launching.
CPU explicitly disables GPU MPI support. Both paths clear LF profiling/tracing
and Kokkos tool overrides before each run, and place temporary files under WO2.
Slurm uses exact one-node allocation steps with one GPU per HIP rank; CPU tests
reserve no GPUs and may overlap other untimed work. Fresh per-case logs/results
are saved incrementally, including failures.
