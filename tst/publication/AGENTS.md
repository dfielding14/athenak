# AGENTS.md

## Purpose
This directory contains non-regression tooling for publication-oriented PIC and
MHD-PIC campaigns. Scripts here are for orchestrating selected runs, collecting
artifacts, and generating internal quick-look figures.

## Contents
- `pic_publication_manifest.py`
  - Case inventory for the publication priority set.
- `run_pic_publication_suite.py`
  - Executes selected manifest cases and archives logs/artifacts.
- `plot_pic_publication_figures.py`
  - Produces figure candidates from archived artifacts in separate
    `proxy_regression` and `publication_physics` bundles.
- `run_pic_shock_scan.py`
  - Emits/executes scan command grids for shock campaigns.
- `pvtk_particles.py`
  - Reader for AthenaK particle VTK files.

## Usage Constraints
- These scripts are not part of `tst/run_tests.py` regression discovery.
- Heavy shock scans are expected to run on external HPC resources.
- Keep output paths deterministic and rooted under `tst/.codex/` unless
  explicitly overridden by CLI options.
- Treat proxy guardrail plots as regression diagnostics, not publication physics
  evidence.
