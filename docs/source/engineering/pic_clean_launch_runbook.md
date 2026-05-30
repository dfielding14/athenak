# MHD-PIC Clean-Launch Runbook

## Purpose

Use this runbook to establish a clean local launch for the explicit
`paper_mhd_pic` and `extended_mhd_pic` runtime identities. These checks are
host-side development evidence only. They do not qualify MPI runtime, HIP/GPU
runtime, Frontier portability, scientific reproduction, or publication claims.

The governing model boundary is documented in
[MHD-PIC Runtime Model Contract](pic_mhd_model_contract.md). Supported build
scope is documented in
[MHD-PIC Supported Toolchains](pic_supported_toolchains.md).

## Host Serial Build

Build a Debug serial executable in the location expected by the direct parser
harness:

```bash
cd /ccs/home/dfielding/athenak-pic
cmake -S . -B tst/build \
  -DCMAKE_BUILD_TYPE=Debug \
  -DAthena_ENABLE_MPI=OFF \
  -DAthena_ENABLE_OPENMP=OFF
cmake --build tst/build -j 8
./tst/build/src/athena -c
```

The compiled-configuration output must report `built_in_pgens`, double
precision, MPI `OFF`, and OpenMP `OFF`. The `-c` output identifies compile-time
configuration only; it does not select a PIC runtime model.

Use a separate directory for the Release serial launch check:

```bash
cd /ccs/home/dfielding/athenak-pic
cmake -S . -B /tmp/athenak-pic-host-release \
  -DCMAKE_BUILD_TYPE=Release \
  -DAthena_ENABLE_MPI=OFF \
  -DAthena_ENABLE_OPENMP=OFF
cmake --build /tmp/athenak-pic-host-release -j 8
/tmp/athenak-pic-host-release/src/athena -c
```

## Parser Contract Suite

Run the parser contract suite directly after building `tst/build/src/athena`.
The explicit `PYTHONPATH` and working directory are required because the
harness imports `scripts.utils.athena` and resolves `tst/build/src/athena`
relative to the current directory.

```bash
cd /ccs/home/dfielding/athenak-pic/tst
PYTHONPATH="$PWD" python3 - <<'PY'
from scripts.particles import pic_parser_contract_guards

pic_parser_contract_guards.run()
if not pic_parser_contract_guards.analyze():
    raise SystemExit("PIC parser contract guards failed")
PY
```

This suite must exercise positive launches for `paper_mhd_pic` and
`extended_mhd_pic`, then reject unsupported modes and invalid cross-mode
compositions. In particular, `paper_mhd_pic` must reject direct-current CT,
Hall, reduced ion-neutral, and adaptive-delta-f extension selections.

## Runtime Identity

For rank zero cosmic-ray launches, AthenaK prints one identity line beginning
with:

```text
PIC runtime model: physical_mode=<mode> state=<state> C=<value> ...
```

The identity line must be retained with the run log. For both `paper_mhd_pic`
and `extended_mhd_pic`, expect:

```text
state=momentum_p_over_m
deposition=tsc
restart_schema=7
```

For `paper_mhd_pic`, also expect:

```text
physical_mode=paper_mhd_pic
induction=ideal_mhd_only
```

An `extended_mhd_pic` launch without the Hall source also reports
`induction=ideal_mhd_only`. A separately selected
`pic_cr_hall_mode=current_to_ct_experimental` extension reports
`induction=cr_current_to_ct`. An extension identity is not paper-mode evidence.

## Optional Host Build Checks

Use separate build directories for OpenMP and MPI compile checks:

```bash
cd /ccs/home/dfielding/athenak-pic
cmake -S . -B /tmp/athenak-pic-host-openmp \
  -DCMAKE_BUILD_TYPE=Debug \
  -DAthena_ENABLE_MPI=OFF \
  -DAthena_ENABLE_OPENMP=ON \
  -DKokkos_ENABLE_OPENMP=ON
cmake --build /tmp/athenak-pic-host-openmp -j 8
OMP_PROC_BIND=false /tmp/athenak-pic-host-openmp/src/athena -c

cmake -S . -B /tmp/athenak-pic-host-mpi \
  -DCMAKE_BUILD_TYPE=Debug \
  -DAthena_ENABLE_MPI=ON \
  -DKokkos_ENABLE_MPI=ON
cmake --build /tmp/athenak-pic-host-mpi -j 8
/tmp/athenak-pic-host-mpi/src/athena -c
```

An MPI compile check is not an MPI rank-launch check. Qualifying MPI and
HIP/GPU runtime checks remain controlled Frontier work.

## Sanitizer Diagnostic

On the current Cray/AMD login toolchain, an AddressSanitizer executable aborts
at startup with an AddressSanitizer ODR-violation diagnostic for duplicated
`.str` globals. Do not record that as a passing ASan run and do not suppress
the diagnostic in qualifying evidence. Use the UBSan-only fallback command in
[MHD-PIC Supported Toolchains](pic_supported_toolchains.md) for bounded host
diagnostics while the ASan startup issue remains unresolved.

## Frontier Stop Boundary

Do not submit a Frontier job, initialize a real node-hour ledger, or treat a
hand-submitted job as qualifying evidence until all of the following are true:

1. The reviewed storage policy records the user-selected Orion-only bulk
   evidence root, the Project Home ledger/control-plane mirror, the explicit
   durability risk, and `ledger_genesis_allowed=true`.
2. The authorized Project Home mirror and initialized mirrored ledger exist.
3. The exact immutable control-plane version authorized by the storage policy
   is installed under the Frontier PIC root.
4. The immutable pre-submit manifest is created and the installed
   `submit_frontier_job.sh` wrapper is used.

## Verified Against

- `CMakeLists.txt:35`
- `CMakeLists.txt:50`
- `CMakeLists.txt:66`
- `src/particles/particles.cpp:1129`
- `tst/scripts/particles/pic_parser_contract_guards.py:29`
- `tst/publication/frontier_control_plane/control_plane_common.py:187`
- `tst/publication/frontier_control_plane/submit_frontier_job.sh:2`
