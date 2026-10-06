# WO2 Task 2: optional adjacent LF half-sweep merging

Triage: **adapted and implemented**. The default-off implementation is suitable for the tested, narrow eligibility
scope. With Task 1's timestep bound and the final-refill fix, paired GPU timings
show 1.77× and 1.49× median speedups in a separate eligible thermal-wave fixture.
The final build and regression pass containing all accepted WO2 tasks is still
pending; this report describes the independently composed Task 1 + refill + Task 2
CPU/HIP binaries.

## Scope and scheduling

`<time>/sts_merge_half_sweeps=true` enables merging only for active CGL, LF-only
RKL2 on a uniform periodic mesh, with zero collision/scattering rates, disabled
collision limiters and `strict_admissibility=false`. AMR/SMR, passive CGL,
nonperiodic/user boundaries, additional diffusion, and the other coupled/source
features excluded by the driver predicate use the ordinary schedule. The absent
parameter and explicit `false` preserve the ordinary path without snapshot
allocation. The original strict wave inputs remain unchanged. The original
paper turbulence deck has strictness false, but its limiter settings make it
ineligible; it receives no claimed Task 2 speedup.

An eligible cycle may defer its final LF half until the next cycle's first half.
Output, restart, final-cycle and mesh-change barriers require completed states.
A late wall-clock stop flushes the old pending half and recalculates the next dt
before the final checkpoint. Collision-rate updates are never moved across a
merged boundary because nonzero rates are ineligible.

Each attempted merge snapshots conserved and primitive state plus diagnostics
and timestep state. A fresh advective-CFL check and admissibility checks make one
collective MPI accept/reject decision. A rejected attempt restores the state,
clears STS registers that could contain nonfinite trial values, completes the old
half, and resumes the ordinary schedule with a refreshed dt. Strict mode uses the
ordinary path because a fatal strict check cannot be transactionally recovered.

## Correctness evidence

Permanent tests cover absent/false identity, collision/limiter/strict/explicit
fallback, sparse and every-cycle output, mixed time/cycle output, restart,
wall-clock flush, smooth temporal refinement, and a deterministic physical CFL
rejection fixture. CPU and HIP core coverage passed; their MPI wrappers passed
with one and four ranks. The standalone HIP checks also retain byte comparisons
against the corrected original baseline. The MPI rejection fixture accepts three
trials and rejects two for CFL on both rank counts; final fields agree exactly,
with no pending duration.

The first permanent suite had nine passes and one test-reader failure on each
backend: a restarted one-row history needs `np.atleast_2d` around `np.loadtxt`.
The corrected restart, mixed-output and temporal tests all passed again, as did
both MPI tests. The mixed-output test also removes an existing `dcycle` key so
that the requested time cadence is actually exercised. Original logs are retained;
no C++ fix or scientific tolerance relaxation was needed for these harness fixes.

At the common time 0.04921753108038256, the integrated smooth merged-versus-ordinary
parallel-temperature RMS differences on three successively halved cycle timesteps
were 1.1841361447e-10, 2.9782991616e-11 and 7.4725521222e-12, exactly matching between
CPU and HIP. The observed orders were 1.9913 and 1.9948. The original analytic
amplitude/phase checks and independent conservation/health checks passed.

The stationary wall fixture exercises a pending-half flush after an unpredicted
wall stop; the resumed-run test checks synchronization and conserved quantities.
The actual CFL rejection branch is exercised. Corrupted-state admissibility
rollback was reviewed but is not established by a deliberate nonfinite runtime
injection. AMR is fenced out and was not separately exercised by these Task 2
runtime tests. These tests do not establish a global nonlinear stability theorem.

## Paired GPU timing

All 16 apps ran serially in one exclusive window in allocation 5628672: one warmup
per mode, followed by three alternating ordinary/merged repeats for each ratio.
The fixture is independently specified: 64³ cells, eight 32³ MeshBlocks, 96 cycles,
periodic collisionless nonstrict kinematic LF with an x-directed thermal wave.
It retains amplitude 1e-4, decay tolerance 0.03 and phase tolerance 0.02 from the
rotated-decay input. Safe kernels, full diagnostics and weighted registers are
explicit; profiling is disabled. The solver timer starts after initialization
and initial output; it includes evolution, in-run/final output and final problem
diagnostics. Process launch and initial setup are excluded and recorded separately.
The hardware/toolchain is the same MI250X/gfx90a, CCE20, ROCm6.4.2, Cray MPICH9.0.1
and Kokkos08ceff92 configuration described in the parent report. This measures an eligible smooth LF workload,
not the unchanged ineligible paper configurations or turbulent active CGL.

| STS ratio cap | Ordinary seconds/cycle: median [min, max] | Merged seconds/cycle: median [min, max] | Median speedup | Stages/cycle, ordinary → merged |
|---|---|---|---|---|
| 2 | 0.170389 [0.168536, 0.223745] | 0.096444 [0.096180, 0.097465] | 1.7667× | 6 → 3.03125 |
| 20 | 0.358103 [0.356936, 0.362271] | 0.239928 [0.239246, 0.243127] | 1.4925× | 14 → 9.05208 |

The relatively wide ratio-2 ordinary timing range is retained, with no discarded
sample. Every enabled run accepted 95 merges, rejected zero, and ended with zero
pending duration. Snapshot time was 0.008031–0.008077 seconds per 96-cycle run.
All 16 apps passed their unchanged analytic checks, conserved mass exactly in
reported history, and ended with all six health counters zero, including
`lf_hardbd` in the separate posthoc audit.

All eight repeat-to-repeat physical comparisons passed. In fact, all 40 compared
files were raw-byte identical, including all 16 restart comparisons. The audit
also independently applied only the documented 36 unused root `RegionIndcs`
coarse-index bytes for restart comparison; it permitted no other normalization.
Raw hashes and raw differing-byte offsets are retained. The original timing
results were left unchanged by the audit.

For six double-precision variables, two ghosted snapshots cost 96 bytes per stored
cell and 192 bytes of copy traffic per accepted attempt. This fixture stores
8 × 36³ cells: 35,831,808 bytes (34.17 MiB) of snapshots and 71,663,616 bytes
(68.34 MiB) of copy traffic per attempt. Rejection adds restoration and STS-register
clearing work. These costs are included in the timing above; they can outweigh
savings for smaller workloads or frequent output barriers.

## Provenance and retained artifacts

The [manifest](evidence-manifest.json) records hashes and original paths for every
archived artifact. The [source patch](evidence/task2-permanent-candidate.patch)
is against Task 1 commit `822cad5e8a49f5f510bece9dcb4b601cacc6f205`; its SHA256 is
`91c394625528eb074746aa34081dd1d1f9474b6faa9e5f8c48327db0fc6a9fe1`.
The production application differs from the tested driver only in indentation
and line wrapping at report time; final combined-build verification remains due.

Tested CPU binary SHA256:
`4d0b368c12e766ce021c58c2a0a34bc3c90f811225e26713864a7de8600c346f`.
Tested HIP binary SHA256:
`c729a817c7508016ea73a75ac0a579b5b28045a1f8c71f071423a8870106f3b0`.
Compatible full Task 1 objects were reused after recompiling the audited Driver
layout consumers and diagnostic pgen; original-baseline Mesh-layout objects were
not mixed into these links. [CPU](evidence/cpu-build.json) and
[HIP](evidence/hip-build.json) manifests retain exact compiler/link commands.

The archive includes [timing results](evidence/timing-results.json), the separate
[repeatability/health audit](evidence/timing-repeatability.json), all 16 run
commands/logs, staged inputs, runners, source/test package, test XML/logs,
[CPU temporal refinement](evidence/cpu-temporal-refinement.json) and
[HIP temporal refinement](evidence/hip-temporal-refinement.json). Large original
field/restart files remain under the recorded WO2 paths. The earlier chronological
[audit](evidence/audit.txt) retains prior unresolved states followed by their
resolution; this README is the final Task 2 summary for the tested composition.
