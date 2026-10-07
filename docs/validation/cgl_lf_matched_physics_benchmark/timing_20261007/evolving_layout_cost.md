# Evolving-layout cost comparison, 2026-10-07

For the same **900 early passive cycles (100–1000)**, the one-node layout used
**46.8% fewer node-hours**, or **1.88 times as much work per node-hour**, while
taking **4.26 times longer in elapsed wall time**. This is one observed startup
interval, not a controlled scaling study or a prediction of late-time cost.

Both runs use the full 192 × 192 × 384 grid. Their endpoint physical times agree
at the printed precision: **0.05853484 to 0.5214464**. The measurements below use
differences of the original `elapsed` values; they do not substitute medians
from the separate [20-cycle plumbing checks](report.md).

| Layout | Nodes / GPU devices | MeshBlocks | Blocks per device | Elapsed seconds | Node-hours |
| --- | ---: | --- | ---: | ---: | ---: |
| Original `passive-to14` | 8 / 64 | 128 × (48 × 48 × 48) | 2 | 155.53340 | 0.34562978 |
| New `large-block-evolution-passive` | 1 / 8 | 8 × (96 × 96 × 192) | 1 | 661.81403 | 0.18383723 |

`Node-hours = nodes × elapsed seconds / 3600`. These are resources occupied by
each member during the selected work interval. They exclude reservation idle
time, startup before cycle 100, finalization, later work and other steps sharing
the allocation. They therefore must not be added to the whole-allocation costs
in the [interrupted-run review](../interrupted_20261007/report.md).

The wall-clock difference includes communication, solver work, progress logging
and scheduled I/O between the two diagnostic timestamps. Both retained inputs
request histories every 0.02, primitive-plus-B and forcing snapshots every
0.25, and restarts every 1.0. The selected interval crosses the two snapshot
times 0.25 and 0.5, but no restart time. This is not a kernel-only measurement.
The timer is rank zero's elapsed wall clock; it is neither an MPI maximum nor
accumulated GPU time. No measurement fences or barriers were introduced.

The endpoints were printed with six digits after the scientific-format
mantissa decimal. Propagating their rounding gives an elapsed-difference
half-width of 0.000055 s for each member. Agreement of printed physical times
does not establish bitwise identity of their fields. Repeated runs are absent,
so these precision bounds are not performance uncertainty estimates.

The original executable revision is
`ab9b543e7e6a972ebfc026d0c98ea3dfee9cb55b`; the new executable revision is
`635b42562e5e23fc9966c9faa01c9da62134f2ec`. The production source difference adds
interval progress logging and counters in `driver.cpp` and `driver.hpp`; the
other changed source is a unit test. No physics update changed. Retained
manifests and metadata agree on each revision and executable SHA256. The
retained effective inputs differ only in MeshBlock dimensions, final time
(14 versus 3) and progress-report cadence (100 versus 25 cycles). Both preserve
the same seed, forcing, closure, limiters, reconstruction, integration and
output cadence. Exact metadata, input and manifest hashes, commands, layout,
endpoint log lines and definitions are in
[evolving_layout_cost.json](evolving_layout_cost.json).

Node count, block dimensions, block count and blocks per device change together.
Hardware placement, concurrent workloads, output costs and logging cadence may
also matter. This observation cannot isolate the effect of one versus two
blocks per GPU. The cost rises as turbulence develops and LF stage counts
change; do not project this short interval directly to the scientific window.
The original passive run subsequently failed strict admissibility. These cost
measurements establish neither completion of the new preflight nor physical
validation.

The lightweight extraction script reads only retained logs, metadata, input
files, build manifests and Git source differences. It submits no jobs and loads
no field arrays. Reproduce its external JSON with:

```bash
M=/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/physics-benchmark/matched
python3 "$M/scripts/review-evolving-layout-cost.py"
```

The result is `$M/audits/evolving-layout-cost.json`. Its script SHA256 and
analysis checkout revision are recorded separately from both simulation
revisions. Fixed log prefixes through cycle 1000 are hashed, so a later append
to the new run's log does not alter the selected evidence.
