# PIC Large-Machine Validation Plan

This document is the handoff for PIC/MHD-PIC validation that is too large to
finish responsibly on a laptop-class CPU run. It assumes the current branch is
`c/pic-review` and that the local regression/proxy evidence has already been
collected. Re-run the commands below on the target machine and keep the archived
run directories with the branch or PR record.

## Local Evidence Already Collected

The following checks were run successfully on the local workstation:

- `python3 run_tests.py particles --cmake=-DCMAKE_BUILD_TYPE=Debug --cmake=-DAthena_ENABLE_MPI=ON`
  - Result: 25/25 default particle tests passed.
  - Normal discovery excludes `*_publication.py` and helper `*_utils.py` modules.
- `python3 tst/publication/run_pic_publication_suite.py --groups entity_core_proxy,benchmarks_proxy --tier local --enable-mpi --cmake-arg=-DCMAKE_BUILD_TYPE=Debug`
  - Result after running outside the sandbox for MPI sockets:
    `overall_status: ok`.
  - Archived proxy case artifacts under
    `tst/.codex/pic_publication_runs/pic_proxy_archive_retry`.
- `python3 tst/publication/compute_pic_publication_metrics.py --run-dir <proxy_run_dir>`
  - Result: nonempty `metrics_summary.json` for F01-F09 proxy metrics.
- `env MPLBACKEND=Agg MPLCONFIGDIR=/private/tmp/athenak_pic_mplconfig XDG_CACHE_HOME=/private/tmp/athenak_pic_xdg python3 tst/publication/plot_pic_publication_figures.py --run-dir <proxy_run_dir> --bundle all --dpi 160`
  - Result: proxy figures generated; `publication_physics` correctly excludes
    guard-only two-stream/Weibel proxy panels when publication artifacts are
    absent.

## Environment Setup

Use a clean worktree for the branch and record the exact revision:

```bash
git status -sb
git rev-parse HEAD
git diff --check
```

Use a site-supported compiler/MPI/Kokkos stack. The exact GPU CMake flags are
site-specific, but record `./athena -c` in every archived run. For CPU/MPI
validation, this build is sufficient:

```bash
python3 tst/publication/run_pic_publication_suite.py \
  --groups entity_core_proxy,benchmarks_proxy \
  --tier local \
  --enable-mpi \
  --cmake-arg=-DCMAKE_BUILD_TYPE=Release \
  --run-id proxy_release_mpi
```

For GPU-capable systems, use the local site equivalent of the Kokkos backend
flags, for example CUDA/HIP/SYCL settings as appropriate, and keep the full
configure log. Do not silently mix CPU and GPU evidence; name the run IDs with
the backend.

## Required Validation Tiers

### 1. Default Particle Regression

Run the normal regression suite first. It should remain publication-free by
default.

```bash
cd tst
python3 run_tests.py particles \
  --cmake=-DCMAKE_BUILD_TYPE=Release \
  --cmake=-DAthena_ENABLE_MPI=ON
```

Pass criteria:

- Summary reports all default particle tests passed.
- The discovered set does not include `*_publication.py`.
- No `pic_analysis_utils` or other helper modules appear as tests.

### 2. Publication Proxy Archive

Run the manifest proxy archive. This verifies the manifest runner, output
copying, F01-F09 metrics, and figure generation on a complete archived dataset.

```bash
python3 tst/publication/run_pic_publication_suite.py \
  --groups entity_core_proxy,benchmarks_proxy \
  --tier local \
  --enable-mpi \
  --cmake-arg=-DCMAKE_BUILD_TYPE=Release \
  --run-id proxy_release_mpi
```

Then post-process it:

```bash
python3 tst/publication/compute_pic_publication_metrics.py \
  --run-dir tst/.codex/pic_publication_runs/proxy_release_mpi

env MPLBACKEND=Agg MPLCONFIGDIR=/tmp/athenak_pic_mpl XDG_CACHE_HOME=/tmp/athenak_pic_xdg \
python3 tst/publication/plot_pic_publication_figures.py \
  --run-dir tst/.codex/pic_publication_runs/proxy_release_mpi \
  --bundle all \
  --dpi 240
```

Pass criteria:

- `summary.json` has `overall_status: ok`.
- Every selected case has `status: ok` and nonzero `copied_outputs`.
- `metrics_summary.json` has nonempty `F01_mink`, `F01_reflect`, `F02_em`,
  `F03_langmuir`, `F04_twostream`, `F04_weibel`, `F05_bell`, `F06_mso`,
  `F07_crsi`, `F08_crpai`, and `F09_aniso`.
- Proxy two-stream and Weibel metrics use the stability fit, not a positive
  windowed growth search. Their `np*_gamma` values should be at or below the
  proxy test tolerance near zero.
- `figures/proxy_regression` contains figures `01` through `09`.
- `figures/publication_physics` does not contain
  `04_two_stream_weibel_modes.png` unless publication two-stream/Weibel
  artifacts are present.

### 3. Local Publication Campaign

Run the publication-oriented local tier. This is larger than the proxy archive
and may be slow on CPU-only machines.

```bash
python3 tst/publication/run_pic_publication_suite.py \
  --groups entity_core_publication,benchmarks_publication,shock_story \
  --tier local \
  --enable-mpi \
  --cmake-arg=-DCMAKE_BUILD_TYPE=Release \
  --run-id publication_local_release_mpi \
  --keep-going
```

Post-process:

```bash
python3 tst/publication/compute_pic_publication_metrics.py \
  --run-dir tst/.codex/pic_publication_runs/publication_local_release_mpi

env MPLBACKEND=Agg MPLCONFIGDIR=/tmp/athenak_pic_mpl XDG_CACHE_HOME=/tmp/athenak_pic_xdg \
python3 tst/publication/plot_pic_publication_figures.py \
  --run-dir tst/.codex/pic_publication_runs/publication_local_release_mpi \
  --bundle all \
  --dpi 320
```

Pass criteria:

- Manifest `overall_status: ok`.
- EM vacuum publication convergence is positive and consistent across MPI
  ranks.
- Two-stream and Weibel publication metrics report positive growth,
  sufficient `R2`, and MPI-identical time/amplitude series.
- Bell publication reports quiet uncoupled behavior and coupled magnetic
  growth above the publication threshold.
- Multispecies publication reports finite amplitudes, sufficient turn counts
  for uniform/SMR, bounded energy drift, and exact serial/MPI agreement for the
  discrete turn metric.
- CRSI/CRPAI publication cases report separated polarization branches and
  stable MPI overlays.
- Shock-story local deck archives grid, particle, restart, and PVTK outputs.

### 4. Full EM Publication Sweep

The EM publication script has a gated full-resolution mode. Run this on the
larger machine if the paper-quality convergence figure needs the full sweep:

```bash
cd tst
ATHENA_PIC_PUBLICATION_FULL=1 python3 run_tests.py \
  particles/pic_em_vacuum_wave_publication \
  --cmake=-DCMAKE_BUILD_TYPE=Release \
  --cmake=-DAthena_ENABLE_MPI=ON
```

Pass criteria:

- Serial and MPI error tables are produced through the full requested
  resolution list.
- Convergence order is positive and consistent with the default local
  publication run.
- Serial and MPI final errors agree to the script tolerance.

### 5. HPC Shock Campaign

The heavy shock deck is intended for cluster-scale runs:

```bash
python3 tst/publication/run_pic_publication_suite.py \
  --cases amr_shock_publication_hpc \
  --tier hpc \
  --enable-mpi \
  --cmake-arg=-DCMAKE_BUILD_TYPE=Release \
  --mpiexec srun \
  --run-id amr_shock_hpc_release \
  --keep-going
```

If the site uses `mpiexec`, replace `--mpiexec srun` with the site launcher.
The runner will call `<launcher> -n <rank_count> ...`; verify that this matches
the scheduler policy before launching a production allocation.

Pass criteria:

- `amr_shock_publication_hpc` exits cleanly.
- Archived outputs include mesh fields, particle moment fields, restarts, and
  PVTK particle files.
- Shock quick-look figure generation succeeds.
- Particle counts, positions, velocities, gids, species IDs, and tags are
  finite and valid in PVTK outputs.
- Restart cadence is sufficient to resume a failed production campaign.

### 6. Shock Parameter Scans

Use the scan emitter to generate reproducible command lists before running
expensive sweeps:

```bash
python3 tst/publication/run_pic_shock_scan.py \
  --deck hpc \
  --nproc 64
```

Review the emitted commands, then run only the approved grid. Vary one axis at
a time:

- Mach number: `problem/ot_mach`
- CR loading: `particles/ppc`
- Resolution: `mesh/nx*` and `meshblock/nx*`

Archive every scan run with:

- input deck and command line,
- `./athena -c`,
- stdout/stderr log,
- final restart,
- figure bundle,
- metrics summary when available.

## Failure Triage

If any case fails:

1. Keep `summary.json`, each `case.json`, and every `run.log`.
2. Record `git rev-parse HEAD`, `git diff --check`, and `./athena -c`.
3. Identify whether the failure is build/configuration, MPI launcher, runtime
   instability, analysis threshold, or archive/post-processing.
4. For analysis-threshold failures, inspect the raw time series before changing
   thresholds. Do not loosen a threshold unless the diagnostic is still proving
   the intended physical or numerical property.
5. For AMR timestep collapse or NaNs, archive the last good restart and the
   first failed output cycle.

## Final Sign-Off Criteria

Treat the branch as large-machine validated only when all of the following are
true:

- Default particle regressions pass in Release+MPI on the target machine.
- Proxy manifest archive passes and produces coherent F01-F09 metrics/figures.
- Local publication manifest passes and produces coherent publication figures.
- Full EM publication sweep passes if full-resolution convergence is required
  for the paper figure.
- HPC shock deck passes, or any remaining shock limitations are explicitly
  scoped as future campaign work rather than branch correctness work.
- All generated `metrics_summary.json` files and figure bundles are archived
  with the final PR/review record.
- `git diff --check` is clean after any fixes made during the campaign.
