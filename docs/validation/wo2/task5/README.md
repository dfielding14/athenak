# WO2 Task 5: representation-aware LF physical boundaries

Triage: **adapted and implemented**. Representation-aware boundary fills replace
the active-CGL constructor fence. Task 4 passive mode has a separate periodic-only
restriction because its energy slot has a different meaning.

The selected design gives boundary conversion an explicit representation flag.
During an LF sweep, built-in magnetic ghosts are filled before hydrodynamic ghosts.
Each fixed-inflow face decodes its immutable stored total-E/A state using those
final face fields and writes mu exactly once. It never reads concurrently written
cell variables. Copy/reflection faces propagate mu, and later fixed faces overwrite
corners from their own prescribed states. Outside LF the original fill order and
ordinary A encoding remain. `CGLMHD::PrimToCons` emits the MHD module's current
representation for all destination arrays, including temporary buffers used by
user callbacks. Scalar `SingleP2C_CGLMHD` retains its ordinary-A meaning.

The Task 0 post-prolongation refill now calls `ApplyPhysicalBCs`, including user
callbacks, so transverse corner donors are fresh. The new callback contract requires
idempotent ghost fills, fresh face-to-cell magnetic fields before `PrimToCons`, and
no active-cell updates or one-time side effects. Direct conserved writes must respect
`cgl_slot_representation`. The new `cgl_lf_boundary` problem supplies independently
encoded stored inflow and a complete primitive-based callback.

Validation already executed on Frontier job 5628672:

- 72 successful application runs across CCE20 Serial and HIP/gfx90a, one/four ranks,
  explicit LF and RKL2. Uniform analytic states include nonzero velocity, oblique B,
  pressure anisotropy and a passive scalar, on 1D, 2D physical corners, 25-block SMR,
  and 3D faces/corners. All active cells and the full one-cell LF halo meet the
  independent constant-state reference. Binary field tolerance is 2e-7 because
  the output format is float32.
- Nonuniform 1D stored inflow and the independent primitive callback give identical
  field outputs in all eight backend/rank/integrator pairs. Full-precision restart
  parsing independently compares every conserved value and face field, including
  ghosts: maximum difference is exactly zero in every pair. Parameter strings and
  diagnostics differ by construction between boundary types and are not described
  as raw-file equality. Source data and raw hashes remain retained.
- The first CPU harness stopped on a missing override parameter before one test
  could launch. The input was corrected to declare `amp`; its failed artifact and
  the successful rerun are retained. No numerical tolerance was changed.
- Independent source review found no blocking write race or representation gap in
  the active-CGL explicit/RKL2 paths.

Artifacts: `runs/{cpu-prototype,cpu-prototype-recheck,hip-prototype,
cpu-mpi-prototype,hip-mpi-prototype}/results.json`, `inflow-user-full-precision.json`,
`build/{cpu,hip}/manifest.json`, and `task5-complete-draft.patch`, all beneath the
user's WO2 test directory. The early `task5-core-draft.patch` is explicitly marked
superseded: it contained an obsolete restriction to the primary conserved array.

The repository regression suite passed 14 CPU and 14 HIP checks, including
one/four-rank tests and an isotropic nonuniform inflow reference. The same 28
checks passed against the combined Tasks0/1/3/6/5 build; their JUnit results and
binary/source hashes are retained in the acceptance manifest. Its analytic
uniform state includes the complete one-cell LF halo at initialization and cycle4.

An independent negative control restores the baseline representation-blind boundary
and primitive conversion while bypassing only its constructor fence. All four
combinations of explicit/RKL2 and stored-inflow/user boundary fail the physical
oracle or actual CGL admissibility check after entering evolution. Input-parser
failures do not count. This establishes that the regression detects the defect.

No timing benefit is claimed for this correctness change. The temporary allocations in the test callback are a regression fixture,
not a proposed application boundary implementation.

All changed numerical source matches the tested combined build byte-for-byte;
the new test initializer has only token-verified whitespace cleanup. See
[source-integration.json](source-integration.json).
