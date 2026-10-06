# Task 3: compact LF communication

Triage: **adapted**. Once-per-sweep A/mu conversion and light primitive refresh
already existed; this task implements compact communication while retaining all
required boundary, magnetic, and projection work.

The candidate sends only IEN/IAN during each multilevel LF stage. It preserves
all boundary fills, magnetic communication, coarse/fine projection, restriction,
Bcc refresh, and the existing A/mu state lifecycle. The baseline already converts
once per sweep and uses a light primitive refresh; the remaining change is the
communicated variable range, not removal of required conversion work.

Generic cell-centered pack/unpack and coarse/fine flux-correction APIs accept a
contiguous variable range. Existing callers default to the full range. During
LF stages the MHD task path selects two variables, including receive sizes;
other tasks retain full exchanges. Shearing keeps its existing full exchange.
Allocated buffer capacity remains unchanged. With one scalar, the communicated
payload is 2/7 of the original, a 71.43% reduction in those cell-centered and flux
messages. This is a byte-count result, not a claimed wall-time or memory-capacity
reduction.

The reference is the immutable corrected-compiler baseline plus the independent
Task0 physical-corner refill. The Task6 same-level flux synchronization is excluded
from both sides of this proof and has its own patch.

## Correctness

Actual CPU and HIP runs on one and four MPI ranks preserve every recorded state:
110 snapshots on one rank and 440 snapshots across four ranks, per backend, with
zero differing snapshots. The 12/48 evolved E/mu flux snapshots also agree exactly.
Every history, binary field output, and physical restart state is bitwise equal.
Only the separately documented 36 bytes of uninitialized root coarse-index header
storage are normalized; raw hashes and exact ranges are retained.

Task3's original diagnosis also demonstrates why magnetic boundary fills remain:
skipping them changes corner magnetic ghosts, then reconstructed pressure, then
subsequent LF fluxes. Compact transport alone does not produce that change.

## GPU throughput

All six case/rank comparisons passed every output comparison. The measurements
used allocation 5628672, one warmup and three alternating paired repeats, fixed
12 cycles, LF profiling disabled, initial/final state output, and identical inputs
for both binaries. Values below are reference median solver-reported seconds per
cycle divided by candidate median; values above 1 favor the candidate. This
solver timer starts after initialization and initial outputs and includes the
evolution loop, evolved/final outputs, and final problem-generator validation.

| Case | 1 rank | 4 ranks |
| --- | ---: | ---: |
| smr | 0.9926 | 1.0335 |
| smr_varied | 0.9921 | 1.0154 |
| smr_outflow | 1.0369 | 0.9999 |

The changes range from about -0.8% to +3.7% and are mixed against the observed
run-to-run spread. These cases do not establish a consistent throughput benefit.
The supported conclusions are reduced message volume and preserved physical
results. No large-grid or scaling improvement is claimed.

## Provenance and archive

The compact [artifact index](artifact-index.json) preserves original absolute
paths and SHA256 hashes; raw outputs remain under
`/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/p1-research`.
The archive includes [trace summaries](trace-summary.json), all original
CPU/HIP output comparison hashes, [timing repeats](timing-summary.json), and
the rejected magnetic-boundary [negative control](negative-control.json).
The timing harness allowed diagnostic summation roundoff for shared Task7 use;
a separate offline audit confirms every Task3 history is actually byte-identical.
Restart state compares exactly after only the documented 36 header bytes.

Isolated builds added forwarding overloads to link unchanged baseline callers;
production uses the equivalent default arguments and requires recompiling callers.
The combined Task0/1/3/6/5 full CPU/HIP builds exercise that production API.
No scratch forwarding overload is included in the production patch.

The following paths are relative to the retained raw evidence root, unless
provided as files in this archive.

## Review artifacts

- Production patch after Task0 refill: `task3.after-refill.patch`.
- Full scratch source: `source_task3/`.
- Stage and output evidence: `runs/task3-{cpu,hip}/`.
- GPU timing, inputs, commands, raw hashes, and all repeats:
  `runs/optimization-task3-time-5628672/results.json`.
- Independent Task6 integration: `task6-sync.after-task3.patch`.

Patch SHA256: `accaf5fd08602b4905ccede591f647a4f8f0d5054c5abbac4229d8f249754b56`.
The production patch contains no scratch tracing or new runtime switch.
