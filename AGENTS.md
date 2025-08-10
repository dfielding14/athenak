# Repository Guidelines

## Project Structure & Module Organization
- `src/`: C++ sources; entrypoint `main.cpp`, modular subdirs (`hydro/`, `mhd/`, `z4c/`, `outputs/`, etc.).
- `kokkos/`: vendored Kokkos dependency used by the build.
- `inputs/`: sample `.athinput` problem files.
- `tst/`: regression harness (`run_tests.py`) and test suites under `tst/scripts/`.
- `docs/`: focused notes (e.g., SFB turbulence).
- `scripts/`: utility runners (e.g., SLURM example).
- Build artifacts belong in a separate `build/` directory (in‑source builds are blocked).

## Build, Test, and Development Commands
- Configure + build (CPU):
  - `mkdir -p build && cd build && cmake .. && make -j`
- Configure for CUDA (example from CI):
  - `cmake .. -DKokkos_ENABLE_CUDA=On -DKokkos_ARCH_VOLTA70=On -DCMAKE_CXX_COMPILER=../kokkos/bin/nvcc_wrapper`
- Run regression tests (CPU):
  - `python3 tst/run_tests.py hydro mhd radiation`
- Run regression tests (CUDA):
  - `python3 tst/run_tests.py --cmake=-DKokkos_ENABLE_CUDA=On --cmake=-DKokkos_ARCH_VOLTA70=On --cmake=-DCMAKE_CXX_COMPILER=kokkos/bin/nvcc_wrapper`
- Cluster example: see `configure.sh` (HIP/Cray modules, MPI).

## Coding Style & Naming Conventions
- C++: Google‑style lint via `tst/scripts/style/check_athena_cpp_style.sh`; format with `clang-format` using `tst/scripts/style/clang_format_v2.cfg`.
- Python: `flake8` with `setup.cfg` (max line length 90).
- Naming: files `lower_snake_case.cpp|hpp`; types `CamelCase`; functions/variables `lower_snake_case`.

## Testing Guidelines
- Framework: Python regression harness discovers tests under `tst/scripts/<suite>/*.py`.
- Add tests near related domain area (e.g., `tst/scripts/hydro/your_test.py`). Expose `run()` and `analyze()`.
- Local smoke run: `python3 tst/run_tests.py <suite>` or a single test `python3 tst/run_tests.py suite/test_name`.
- CI runs style checks and CPU/GPU regression subsets; ensure new tests are stable and time‑bounded.

## Commit & Pull Request Guidelines
- Commits: short, imperative subject (<72 chars), scope first when helpful (e.g., `mhd: fix flux limiter`).
- PRs: include problem context, what/why, minimal repro or input file, and notes on performance impact. Link issues. Ensure:
  - CI green (C++/Python lint + regressions)
  - Added/updated tests and docs where applicable
  - Screenshots or output snippets for new diagnostics helpful

## Security & Configuration Tips
- Prefer out‑of‑source builds; avoid committing `build_*` outputs.
- GPU/MPI flags vary by machine; mirror CI examples or adapt `configure.sh` for your environment.
