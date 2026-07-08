<!-- BEGIN build-memory-table -->
# Test Navigation

## Scope

`tst/` owns AthenaK's regression harness, validation scripts, style checks, and
standalone source-level tests. Runtime decks live under `inputs/`; implementation
changes belong under `src/`. Exploratory artifact and campaign tooling under
`tst/publication/` follows its own `AGENTS.md` and is outside normal regression
discovery.

## Task lookup

| Task | Start in | Also inspect or validate |
| --- | --- | --- |
| Add or change a compiled regression | `tst/scripts/<suite>/` | Referenced `inputs/tests/*.athinput` and a focused `run_tests.py` invocation |
| Change PIC/MHD-PIC validation | `tst/scripts/particles/` | `tst/scripts/particles/AGENTS.md` and the affected source modules |
| Change test build or execution behavior | `tst/run_tests.py`, `tst/scripts/utils/athena.py` | `.github/workflows/main.yml` and `.gitlab-ci.yml` |
| Change load-balancing source checks | `tst/test_mesh_load_balance.py`, `tst/test_pic_static_load_balance.py` | `src/mesh/`, `src/driver/`, and paired restart decks |
| Change C++ or Python style checks | `tst/scripts/style/`, `setup.cfg` | CI lint jobs |
| Change artifact, qualification, or campaign tooling | `tst/publication/` | Read `tst/publication/AGENTS.md` first |

## Harness flow

- Run `run_tests.py` from `tst/`; test names are relative to `tst/scripts/`.
- The runner imports selected modules, builds AthenaK once, then calls each
  module's `run()` and `analyze()`.
- Suite discovery skips `utils`, `style`, `*_utils.py`, and, by default,
  `*_publication.py`. Explicit publication module names still run;
  `ATHENA_INCLUDE_PUBLICATION_TESTS=1` includes them in broad discovery.
- `tst/scripts/utils/athena.py` configures with `cmake3`, builds with `make -j8`,
  and resolves test input names below `inputs/`.
- Every harness invocation removes `tst/build` before and after execution. Do not
  use that directory for retained evidence.
- Top-level pytest files and `tst/publication/` tests are separate from
  `run_tests.py`.

## Local validation

- `cd tst && python3 run_tests.py particles/pic_relativistic_gyro_paper`:
  focused compiled particle regression.
- `cd tst && python3 run_tests.py hydro mhd radiation`: CPU suite selection used
  by CI.
- `python3 -m pytest -q tst/test_mesh_load_balance.py`: standalone mesh
  load-balancer checks.
- `python3 -m flake8`: repository Python lint configured by `setup.cfg`.
- `bash tst/scripts/style/check_athena_cpp_style.sh`: full C++ style check.

Pass multiple build options by repeating `--cmake`, for example
`--cmake=-DCMAKE_BUILD_TYPE=Debug --cmake=-DAthena_ENABLE_MPI=ON`.
<!-- END build-memory-table -->
