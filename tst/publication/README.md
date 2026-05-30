# PIC Exploratory Artifact Toolkit (Internal)

This folder provides a manifest-driven workflow for exploratory PIC/MHD-PIC
artifact generation:

1. Entity-inspired particle-numerics regressions
2. extended engineering proxies and artifact candidates
3. CR-shock AMR/output/load-balance stress workflows

The checked-in manifest is **not** a Sun and Bai (2023) paper-reproduction
manifest. Its remaining publication-oriented case IDs and filenames are
historical compatibility identifiers. Generated summaries and figures classify
them as `engineering_proxy` and mark them
`not_sun_bai_reproduction=true`. Use
`PIC_PRODUCTION_READINESS_PLAN.md` as the only controlling release,
qualification and Frontier-authorization plan.

## Files

- `pic_publication_manifest.py`
  - Historical exploratory case list, grouping, run mode (`module` or direct
    deck), and output artifact globs. It is not a release-gate registry.
- `run_pic_publication_suite.py`
  - Runs selected groups/cases, archives outputs and logs under
    `tst/.codex/pic_publication_runs/<run_id>/`.
- `plot_pic_publication_figures.py`
  - Builds visibly labeled quick-look PNG/PDF figures and checksummed bundle
    manifests from archived artifacts.
- `run_pic_shock_scan.py`
  - Emits Mach/CR-loading/resolution scan commands for Orszag-Tang-style
    shock-rich engineering stress decks. Direct execution is local-tier only.
- `artifact_lineage.py`
  - Writes checksummed lineage companions for retained standalone quicklooks.
- `pvtk_particles.py`
  - Lightweight parser for AthenaK `pvtk/*.part.vtk` particle files.
- `test_pic_artifact_taxonomy.py`
  - Regression checks for engineering-proxy labels, legacy aliases, checksummed
    archives, HPC emit-only safety and the retired-runbook boundary.
- `PIC_LARGE_MACHINE_VALIDATION.md`
  - Retired historical note. Do not execute commands or infer sign-off authority
    from that file.
- `PIC_PRODUCTION_READINESS_PLAN.md`
  - Canonical release, qualification, archival and Frontier-authorization plan.

## Typical Workflow

1. Run local exploratory artifact groups:
```bash
python3 tst/publication/run_pic_publication_suite.py \
  --groups extended_entity_engineering,extended_benchmark_engineering,shock_engineering \
  --tier local \
  --local-exploratory-run
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
  --bundle engineering_proxy

python3 tst/publication/plot_pic_publication_figures.py \
  --run-dir tst/.codex/pic_publication_runs/<run_id> \
  --bundle extended_engineering
```

4. Emit Orszag-Tang-style engineering stress scan commands:
```bash
python3 tst/publication/run_pic_shock_scan.py \
  --tier hpc \
  --deck hpc \
  --nproc 64
```

Every generated scan shell file is intentionally inert planning material. The
HPC tier also rejects `--run`. Translate selected planned commands into immutable
per-submission snapshots and execute them only through the canonical plan's
approved pre-submit manifest, validator, ledger and scheduler-wrapper workflow.
The local tier retains `--run --local-exploratory-run` for bounded workstation
engineering checks. Direct helper execution is rejected on Frontier and inside
Slurm allocations.

## Engineering Stress Decks

- `inputs/tests/pic_amr_shock_lb_publication_local.athinput`
  - Richer diagnostics and moderate runtime for local/GPU workstations.
- `inputs/tests/pic_amr_shock_lb_publication_hpc.athinput`
  - Heavier engineering-stress profile. It is not a paper Section 5.4
    reproduction deck.

Both decks include:
- grid diagnostics (`rho`, `|B|`, `J^2`, particle moment fields),
- coupled PIC controls,
- `pvtk` particle dumps for phase-space/spectrum post-processing,
- restart cadence for bounded engineering recovery checks.

## Taxonomy Regression

Run the artifact-safety regression after changing this toolkit:

```bash
python3 tst/publication/test_pic_artifact_taxonomy.py
```

## Notes

- For exploratory MPI parity artifacts, configure with:
```bash
python3 tst/publication/run_pic_publication_suite.py \
  --groups extended_entity_engineering,extended_benchmark_engineering \
  --tier local \
  --local-exploratory-run \
  --enable-mpi \
  --mpiexec /opt/homebrew/bin/mpiexec
```
- Standard `tst/run_tests.py` suite discovery skips `*_publication.py`
  modules by default. Use this manifest runner for exploratory artifact runs, pass
  an explicit test such as `particles/pic_bell_growth_publication`, or set
  `ATHENA_INCLUDE_PUBLICATION_TESTS=1` to include publication modules in broad
  suite discovery.
- Heavy cases should be run only under an explicitly approved execution policy.
- The plotting script is robust to missing cases; it will generate only figures
  with available artifacts.
- `proxy_regression` includes guardrail diagnostics (including stability checks).
- `proxy_regression` and `publication_physics` remain accepted as deprecated CLI
  aliases for `engineering_proxy` and `extended_engineering`. Neither implies
  physical qualification.
- The plotting tool rejects `--bundle sun_bai_2023_reproduction` until a
  qualified reproduction bundle is implemented through the canonical plan.
- For Mach scans, use `problem/ot_mach`; for CR loading, vary `particles/ppc`;
  for resolution scans, vary `mesh/nx*` and `meshblock/nx*`.
