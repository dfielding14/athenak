# MHD-PIC Clean-Launch Runbook

## Purpose

Use this runbook to establish a clean local launch for the explicit PIC runtime
interface. These checks are host-side development evidence only. They do not
qualify MPI runtime, HIP/GPU runtime, Frontier portability, scientific
reproduction, or publication claims.

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
from scripts.particles import pic_parser_contract_guards_vl2_tsc

pic_parser_contract_guards_vl2_tsc.run()
if not pic_parser_contract_guards_vl2_tsc.analyze():
    raise SystemExit("PIC parser contract guards failed")
PY
```

This suite exercises coupled `rk1`, coupled Sun--Bai `vl2`, test-particle, and
bounded adaptive/damped expanding-box launches, then rejects unsupported
values and inconsistent compositions. In particular, `vl2` requires the
explicit coupled MHD, TSC deposition, and conservative-feedback controls;
`pic_cr_hall_mode=full` requires `vl2`. Hall accepts only `off` or `full`.

## Runtime Identity

For rank zero cosmic-ray launches, AthenaK prints one identity line beginning
with:

```text
PIC runtime model: integrator=<integrator> state=<state> C=<value> ...
```

The identity line must be retained with the run log. A coupled VL2 cosmic-ray
Boris launch reports:

```text
state=momentum_p_over_m
deposition=tsc
restart_schema=8
```

For a coupled VL2 launch with `pic_cr_hall_mode=off`, also expect:

```text
integrator=vl2
induction=ideal_mhd_only
```

A full-closure launch instead reports `induction=cr_hall_full` and then prints
the positive `PIC CR-Hall: background_ion_q_over_mc=<value>` normalization on
the next line. The same identity line reports the explicit background,
feedback, delta-f, damping, and expanding-box selections. `rk2` denotes Heun;
it does not select the Sun--Bai staged chronology.

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

A fresh ledger bootstrap and a registered-science submission are separate
unlocks. A hand-submitted job is never qualifying evidence.
The sole bulk-evidence root is
`/lustre/orion/ast207/proj-shared/dfielding/PIC`; Kronos is not an execution,
continuation, or bulk-publication path. Project Home is limited to the small
ledger and control-plane mirror.

For the one-time pre-genesis bootstrap only:

1. Install the exact immutable control-plane version in both the authorized
   Orion root and Project Home mirror from a reviewed clean Git commit. The
   production installer rejects modified or untracked control-plane source.
2. Promote a reviewed policy that records the user-selected Orion-only bulk
   evidence root, the Project Home ledger/control-plane mirror, the explicit
   durability risk, `ledger_genesis_allowed=true`, and no initialized
   `ledger_genesis` fields.
3. Run the installed `initialize_frontier_ledger.py` procedure exactly once to
   create the mirrored genesis event and paired read-only genesis anchors.

The open `ledger_genesis_allowed=true` policy authorizes initialization only.
It is not a reservation or submission unlock. Before any registered-science
reservation or submission:

1. Confirm that the initialized mirrored ledger and paired read-only genesis
   anchors exist.
2. Update and promote the reviewed policy with
   `ledger_genesis_allowed=false` and the exact initialized genesis event and
   mirror-receipt digests.
3. Confirm that the promoted policy authorizes the installed immutable
   control-plane version.
4. Create a clean-candidate freeze and promote
   `science_submission_freeze.status=authorized` with its exact Orion manifest
   path and SHA-256 digest.
5. Create the immutable pre-submit manifest and use only the installed
   `submit_frontier_job.sh` wrapper.

The installed wrapper reserves and durably mirrors worst-case node-hours. The
reservation binds the active policy and active-promotion digests; those digests
are rechecked before dispatch, before scheduler-ID attachment, and in the
trampoline before workload execution. Wrapper metadata reads used to construct
the scheduler request require a live mirrored reservation in `reserved` state.
Immediately before submitting its trusted trampoline with
`sbatch --hold --export=NIL` through a closed environment pinned to the
`frontier` cluster, the wrapper durably records `scheduler_dispatch_started`.
The scheduled runner receives immutable bindings through argv and rechecks
them against the mirrored ledger without reconstructing the submitter's login
environment. It then records
`scheduler_job_id_received` and `submitted_not_attached`, attaches the verified
scheduler ID, and only then runs `scontrol release`.

Only a failure known to precede `scheduler_dispatch_started` may automatically
cancel an unused reservation. Any ambiguity at or after that marker retains
accounting and requires reviewed Slurm inspection and reconciliation. Reconcile
a terminal job from the durable marker if it reaches a terminal state before
attachment. Terminal reconciliation cleanup is idempotent: rerunning it clears
a matching stale pending marker without duplicating accounting. Slurm stdout is
restricted to `${PIC_ROOT}/logs/slurm/%x.%j.log`. Each run uses exactly
`/lustre/orion/ast207/proj-shared/dfielding/PIC/runs/<campaign>/<submission-id>`;
reruns use a new submission ID.

Immediately before workload execution, the launch wrapper syncs the inherited
runtime-allowlist descriptor, changes it to mode `0400`, syncs it again, syncs
its pinned directory descriptor, closes both descriptors, and removes both
environment bindings. Installed snapshots, clean-candidate freezes and
pre-submit snapshots are staged and published relative to pinned authorized
parent descriptors, with a fail-closed lexical-parent identity check after
publication. This closes ancestor-swap pathname redirection during publication;
it does not establish an integrity boundary against a malicious process running
as the same Unix UID.

The separately policy-bound Frontier admission smoke is a narrow
non-production exception documented in
`tst/publication/frontier_control_plane/README.md`; it is not a
registered-science unlock.

## Verified Against

- `CMakeLists.txt:35`
- `CMakeLists.txt:50`
- `CMakeLists.txt:66`
- `src/particles/particles.cpp:1129`
- `tst/scripts/particles/pic_parser_contract_guards.py:29`
- `tst/publication/frontier_control_plane/control_plane_common.py:187`
- `tst/publication/frontier_control_plane/submit_frontier_job.sh:2`
