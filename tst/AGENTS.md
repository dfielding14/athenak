<!-- BEGIN build-memory-table -->
# Regression Test Navigation

## Scope

This subtree owns the current pytest-based regression harness, regression input files,
style checks, and a legacy `run()`/`analyze()` test framework. Run the current harness
from `tst/`; its relative paths, output readers, and temporary build layout assume that
working directory.

## Task lookup

| Task | Start in | Also inspect or validate |
| --- | --- | --- |
| Run all tests for one device class | `tst/run_test_suite.py` | CPU/MPI/GPU commands below and CI flags |
| Add a current regression test | `tst/test_suite/<area>/test_*_{cpu,mpicpu,gpu}.py` | Closest test in that area, device suffix contract |
| Add or change regression input data | `tst/inputs/` | The test's command-line overrides and expected outputs |
| Change build/run/convergence helpers | `tst/test_suite/testutils.py` | Every test area and `vis/python/athena_read.py` |
| Change source style checks | `tst/test_suite/style/` | `CPPLINT.cfg`, `setup.cfg`, GitHub CI |
| Maintain a legacy test | `tst/run_tests.py`, `tst/scripts/` | Its module's `run()` and `analyze()` functions |
| Change an output reader used by tests | `vis/python/athena_read.py` | Producing writer in `src/outputs/` and representative tests |

## Important flow

1. `run_test_suite.py` chooses tests by `_cpu`, `_mpicpu`, or `_gpu` in the filename; a
   `--test` path infers the device class from that suffix.
2. `testutils.clean_make` configures the repository into `tst/build`, builds
   `tst/build/src/athena`, and links `tst/inputs` into the executable directory.
3. Pytest cases call `testutils.run` or `mpi_run`, pass `block/parameter=value` overrides,
   read generated files with `athena_read`, assert errors/convergence/invariants, and clean
   their outputs.
4. The legacy runner imports modules under `tst/scripts/`, builds once per configuration,
   then calls each module's `run()` and `analyze()`.

## Local validation

- `python3 run_test_suite.py --style`: C++ style and repository Python flake8 checks.
- `python3 run_test_suite.py --test test_suite/rad/test_rad_lwave1d_amr_cpu.py`:
  build and run one focused current test.
- `python3 run_test_suite.py --cpu`: full CPU-selected pytest suite.
- `python3 run_test_suite.py --mpicpu`: full MPI CPU suite; requires MPI toolchain/runtime.
- `python3 run_test_suite.py --gpu "-DKokkos_ARCH_<TARGET>=On -DCMAKE_CXX_COMPILER=<nvcc_wrapper>"`:
  GPU suite with site-specific Kokkos flags.
- `python3 run_tests.py <legacy-test-path>`: run one or more modules relative to
  `tst/scripts/`.

## Local constraints

- Use the device suffix exactly; the runner uses it both for collection and build-mode
  selection.
- The current harness recreates `tst/build`; do not store hand-built artifacts there.
- Keep numerical thresholds explicit and justify updates with an algorithmic change or
  reproducible platform effect. Test both absolute error and convergence when the test is
  designed to do so.
- Clean generated files in `finally` blocks so a failed parameterization does not
  contaminate later cases.
<!-- END build-memory-table -->
