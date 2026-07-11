# AGENTS.md

## Purpose
This directory contains non-regression tooling for exploratory PIC and MHD-PIC
artifact campaigns. Scripts here orchestrate selected runs, collect artifacts,
and generate internal quick-look figures. The checked-in manifest is not a
Sun and Bai (2023) reproduction manifest.

## Contents
- `pic_publication_manifest.py`
  - Historical exploratory artifact-case inventory. Publication-oriented names
    must not be interpreted as qualification.
- `run_pic_publication_suite.py`
  - Executes selected manifest cases and archives logs/artifacts.
- `plot_pic_publication_figures.py`
  - Produces figure candidates from archived artifacts in separate
    `engineering_proxy` and `extended_engineering` bundles. Historical bundle
    names remain deprecated CLI aliases only.
- `run_pic_shock_scan.py`
  - Emits Orszag-Tang-style shock-rich engineering stress grids. Direct
    execution is acknowledged local-only use and is rejected on Frontier.
- `artifact_lineage.py`
  - Writes checksummed companions for retained standalone exploratory artifacts.
- `test_pic_artifact_taxonomy.py`
  - Enforces artifact classifications, lineage hashes, aliases and HPC
    emit-only safety.
- `readiness/`
  - Source-controlled release claims, findings, external-artifact inventory,
    storage policy and validation-manifest schema. These records are compact
    control metadata, not qualification evidence.
- `test_pic_readiness_registry.py`
  - Enforces the authorized Frontier root and budget, selected extension claims,
    initial findings, paper-source checksum and minimum validation schema.
- `pvtk_particles.py`
  - Reader for AthenaK particle VTK files.
- `analyze_pic_cr_hall_bell.py` and `analyze_pic_cr_hall_bell_nonlinear.py`
  - Compact linear and nonlinear full-CR-Hall Bell qualification reducers.
- `analyze_pic_parallel_shock_cr_hall_controlled.py`
  - Controlled delayed-injection no-CR/Hall-off/full-Hall shock reducer.
- `analyze_pic_turbulent_dynamo_cr_hall.py`
  - Current-HEAD 128-cubed full-Hall turbulence reducer and Hall-off decision aid.

## Usage Constraints
- These scripts are not part of `tst/run_tests.py` regression discovery.
- `PIC_PRODUCTION_READINESS_PLAN.md` is the only controlling release,
  qualification and Frontier-authorization document.
- `PIC_LARGE_MACHINE_VALIDATION.md` is retired historical evidence only. Do not
  execute commands or infer authority from it.
- Heavy shock scans require an explicitly approved execution policy.
- Keep output paths deterministic and rooted under `tst/.codex/` unless
  explicitly overridden by CLI options.
- Treat proxy guardrail plots as regression diagnostics, not publication physics
  evidence.
- Require generated exploratory summaries and figures to state
  `evidence_class=engineering_proxy` and `not_sun_bai_reproduction=true`.
- Run `python3 tst/publication/test_pic_artifact_taxonomy.py` after changing
  this toolkit.
- Run `python3 tst/publication/test_pic_readiness_registry.py` after changing
  release registries, schema, budget or storage policy.
- Do not cite a checked-in exploratory case as `sun_bai_2023_reproduction`
  evidence until it is promoted through the canonical plan's claims registry.
