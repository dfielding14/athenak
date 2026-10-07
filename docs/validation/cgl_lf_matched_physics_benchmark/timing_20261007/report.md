# Larger-block timing and restart plumbing, 2026-10-07

**Plumbing checks passed. Physical evidence remains inconclusive.** These are
20 startup cycles per member plus three active restart cycles, far before the
previous passive failure at approximately t=2.51. They establish memory/layout,
logging and short restart operation; they do not establish late-time stability,
physical agreement or a scaling result.

Both startup checks used the full 192 × 192 × 384 mesh on **one Frontier node,
eight GPU devices and eight MPI ranks**, with one **96 × 96 × 192 MeshBlock per
rank/device**. All three steps ran on frontier03053 in allocation 5631477.

| Check | Completed cycles | Final physical time | Selected timing intervals | Median seconds/cycle | Median global zone-cycles/s |
| --- | ---: | ---: | --- | ---: | ---: |
| Active startup | 0–20 | 0.007953347 | 2–19, 18 intervals | 0.4153661 | 3.408024 × 10⁷ |
| Passive startup | 0–20 | 0.01225940 | 2–19, 18 intervals | 0.5278932 | 2.681561 × 10⁷ |
| Active restart | 20–23 | 0.009136146 | Cycle 22 only | 0.4393325 | 3.222110 × 10⁷ |

The timing selections exclude the initial no-work report, the first completed
cycle and finalization. Initial input/output setup precedes the no-work report.
Finalization writes full outputs and is included in the last printed interval:
6.295364 s active, 2.723286 s passive and 2.739233 s restart. Those values are
not representative solver-step costs. No scheduled outputs fell inside the
selected startup intervals. Active/passive physical timesteps and LF stage
counts can differ; these measurements must not be interpreted as an isolated
comparison of kernel speeds or extrapolated to mature turbulence.

All three launchers returned zero and confirmed the unchanged binary hash. The
logs have no fatal, strict-admissibility or out-of-memory error. Retained final
histories have zero density-floor, pressure-floor, nonfinite and nonpositive LF
counters. Each rank reports exactly one MeshBlock. The restart begins with a
cycle-20 `warmup` report, then advances through cycles 21, 22 and 23 without
counting the previous segment's cycles as work in the new process.

## Existing branch implementation and log meanings

The audit inspected **84 retained local and origin refs** using Git reads,
without switching branches or fetching. The selected implementation was
`46a6f704b563587025d0faa87fdd7a1623458d57`, “Report interval zone-cycle
throughput,” from `scaling-tests`. Unlike the alternative TRML tracer version,
its numerator uses actual cumulative global MeshBlock updates and remains
correct if AMR changes the block count between reports. The broader TRML
synchronized-region profiler was not imported: it deliberately adds fences
and includes many unrelated subsystem callsites.

The surgical port is commit **635b42562e5e23fc9966c9faa01c9da62134f2ec**.
Existing `elapsed`, `cycle`, `time` and `dt` fields are retained. New fields are:

- `interval_cycles`: completed top-level cycles since the preceding report.
- `wall_seconds_per_cycle`: rank-zero wall-time difference divided by
  `interval_cycles`; an interval mean unless `time/ndiag=1`.
- `zone-cycles/s`: global completed MeshBlock-update difference times active
  cells per block, divided by the same wall-time difference. Ghost zones and
  RK/STS stage counts are excluded. This is a global rate, not a per-GPU rate.

The initial no-work interval is labeled `performance_interval=warmup` and omits
undefined rates. No timing barriers or fences were added. The timer is rank-zero
elapsed wall time, not an MPI maximum or accumulated CPU/GPU time. Communication
and outputs between reports contribute to it.

Every GPU fixture interval passed an independent check against successive
printed timestamps and cycle counts, allowing exactly the rounding intervals
implied by their printed decimal precision. For this uniform mesh the numerator
is **14,155,776 × interval_cycles**. This is a logging-arithmetic check, not a
physical acceptance tolerance.

## Build and retained evidence

The executable revision is 635b42562e5e23fc9966c9faa01c9da62134f2ec; its SHA256 is
`fb11d6355486dce126be1d4564f99dbcc0e9a13e7227a05f3a2f262ae466d47a`.
The production build rebuilt **154** translation units with transitive driver
header dependencies and verified **21** reused linked inputs against the prior
manifest. The source snapshot and reused input hashes were unchanged after
linking. Compilation/linking took **558.57 s** with four concurrent compiles on
frontier03053. Compile warnings and exact commands are retained with the build.

The separate launcher fix, **6fabbc7d0**, replaces initially empty explicit
history inventories with supported glob metadata. An interrupted launcher can
therefore discover histories already on disk; successful finalization still
retains exact inventories. Eighteen matched tests passed, including actual
metadata-writing and analyzer-discovery checks for interrupted and completed
launches. This Python-only change needs no numerical rebuild. Simulation
revision, checkout revision at launch and the exact retained launcher hash are
recorded separately in [metrics.json](metrics.json).

All raw outputs, exact input/override records, commands, logs, restarts, build
manifest/cache and verification script are retained beneath:

```text
/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/physics-benchmark/matched/
  build-timing-hip/
  large-block-plumbing-active/
  large-block-plumbing-passive/
  large-block-plumbing-active-restart/
  scripts/large-block-restart-plumbing.sh
  scripts/review-large-block-plumbing.py
  audits/branch-performance-logging.md
  audits/branch-performance-logging-refs.json
```

The compact [metrics.json](metrics.json) retains exact run commands and hashes.
The verification script reads only retained logs and histories; it launches no
simulation. Re-run it from the matched scratch directory with:

```bash
python3 scripts/review-large-block-plumbing.py
```
