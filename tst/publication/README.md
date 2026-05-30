# PIC Publication Toolkit (Internal)

This folder provides a manifest-driven workflow for the publication-priority
PIC/MHD-PIC set:

1. Entity PIC core validation
2. Athena++ benchmark reproductions
3. CR-shock MHD-PIC storyline

## Files

- `pic_publication_manifest.py`
  - Canonical case list, grouping, run mode (`module` or direct deck), and
    output artifact globs.
- `run_pic_publication_suite.py`
  - Runs selected groups/cases, archives outputs and logs under
    `tst/.codex/pic_publication_runs/<run_id>/`.
- `plot_pic_publication_figures.py`
  - Builds quick-look PNG figures from archived artifacts.
- `run_pic_shock_scan.py`
  - Emits (and optionally executes) Mach/CR-loading/resolution scan commands
    for publication shock decks.
- `pvtk_particles.py`
  - Lightweight parser for AthenaK `pvtk/*.part.vtk` particle files.
- `PIC_LARGE_MACHINE_VALIDATION.md`
  - Required larger-machine validation tiers, commands, pass criteria, and
    archive checklist for publication/HPC campaigns that are too large for a
    laptop-class run.

## Typical Workflow

1. Run local publication groups:
```bash
python3 tst/publication/run_pic_publication_suite.py \
  --groups entity_core_publication,benchmarks_publication,shock_story \
  --tier local
```

2. Plot figures from one run directory:
```bash
python3 tst/publication/plot_pic_publication_figures.py \
  --run-dir tst/.codex/pic_publication_runs/<run_id>
```

3. Generate only one figure bundle when needed:
```bash
python3 tst/publication/plot_pic_publication_figures.py \
  --run-dir tst/.codex/pic_publication_runs/<run_id> \
  --bundle proxy_regression

python3 tst/publication/plot_pic_publication_figures.py \
  --run-dir tst/.codex/pic_publication_runs/<run_id> \
  --bundle publication_physics
```

4. Emit shock scan commands for cluster execution:
```bash
python3 tst/publication/run_pic_shock_scan.py \
  --deck hpc \
  --nproc 64
```

By default, the scan tool only writes command artifacts (`emit-only`). Pass
`--run` to execute commands immediately.

## New Shock Decks

- `inputs/tests/pic_amr_shock_lb_publication_local.athinput`
  - Richer diagnostics and moderate runtime for local/GPU workstations.
- `inputs/tests/pic_amr_shock_lb_publication_hpc.athinput`
  - Heavier profile intended for HPC cluster campaigns.

Both decks include:
- grid diagnostics (`rho`, `|B|`, `J^2`, particle moment fields),
- coupled PIC controls,
- `pvtk` particle dumps for phase-space/spectrum post-processing,
- restart cadence for long campaigns.

## Notes

- For MPI validation evidence, configure with:
```bash
python3 tst/publication/run_pic_publication_suite.py \
  --groups entity_core_publication,benchmarks_publication \
  --tier local \
  --enable-mpi \
  --mpiexec /opt/homebrew/bin/mpiexec
```
- Standard `tst/run_tests.py` suite discovery skips `*_publication.py`
  modules by default. Use this manifest runner for publication campaigns, pass
  an explicit test such as `particles/pic_bell_growth_publication`, or set
  `ATHENA_INCLUDE_PUBLICATION_TESTS=1` to include publication modules in broad
  suite discovery.
- Heavy cases should be run on GPU-capable systems.
- The plotting script is robust to missing cases; it will generate only figures
  with available artifacts.
- `proxy_regression` includes guardrail diagnostics (including stability checks).
- `publication_physics` excludes guard-only two-stream/Weibel proxy traces.
- For Mach scans, use `problem/ot_mach`; for CR loading, vary `particles/ppc`;
  for resolution scans, vary `mesh/nx*` and `meshblock/nx*`.
