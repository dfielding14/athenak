# Synchronized Performance Timing

AthenaK has an opt-in set of coarse timing regions for short performance
diagnostics. The timers attribute asynchronous accelerator work to named parts
of the driver, TRML cooling, frame tracking, and Lagrangian Monte-Carlo tracer
implementation.

The diagnostics are disabled by default. They deliberately synchronize the
accelerator at every region boundary, so they are intended for controlled
benchmarks rather than production runs.

## Enabling the Timers

Add `performance_timing = true` to the existing `<time>` block:

```ini
<time>
integrator = rk2
cfl_number = 0.3
nlim = 200
tlim = 1.0e9
ndiag = 20
performance_timing = true
```

No additional problem-generator or particle setting is required. If the
parameter is absent, it defaults to `false`.

Run AthenaK normally, for example:

```bash
srun -N 2 -n 16 --gpus-per-task=1 ./athena -i benchmark.athinput
```

The report is written to standard output by rank zero when the simulation
finalizes cleanly. A crash or forced cancellation may prevent the report from
being printed.

## How It Works

Each instrumented region does the following:

1. Fence all outstanding Kokkos work.
2. Start a rank-local `Kokkos::Timer`.
3. Execute the instrumented operation.
4. Fence again and accumulate the elapsed time and call count.

The final report uses MPI reductions to compute the maximum and mean time over
all ranks. The maximum is the most relevant value for elapsed-time bottlenecks,
because all ranks ultimately wait for the slowest rank.

The fences are essential for GPU attribution. Without them, a timer could
measure only the time required to launch a kernel while charging its execution
to a later operation. The tradeoff is that the fences suppress some normal CPU,
GPU, and communication overlap and therefore perturb total runtime.

When timing is disabled, `StartPerformanceRegion()` and
`StopPerformanceRegion()` return immediately without fencing or reading a
timer.

## Report Format

A report begins with:

```text
PERFORMANCE_TIMING synchronized=true nested_regions=true run_seconds_max=1.619837e+01
PERFORMANCE_REGION name max_rank_seconds mean_rank_seconds fraction_of_run calls_max
```

It is followed by one line for each region:

```text
PERFORMANCE_REGION stage_tasks 1.503801e+01 1.488858e+01 9.283659e-01 240
PERFORMANCE_REGION trml_cooling 2.823469e-01 2.812812e-01 1.743058e-02 240
PERFORMANCE_REGION particle_seed 1.795832e-03 1.576349e-03 1.108650e-04 120
```

The fields are:

| Field | Meaning |
| --- | --- |
| `run_seconds_max` | Maximum total runtime over MPI ranks after initialization. It includes evolution and final outputs written before the report. |
| `max_rank_seconds` | Largest accumulated region time on any MPI rank. Use this first when locating wall-clock bottlenecks. |
| `mean_rank_seconds` | Region time averaged over all MPI ranks. A large gap from the maximum suggests rank imbalance. |
| `fraction_of_run` | `max_rank_seconds / run_seconds_max`. This is a convenient scale estimate, not an additive profile. |
| `calls_max` | Largest region call count on any rank. Task polling and retries can make this larger than the number of cycles or stages. |

The regions are nested. For example, `trml_cooling` is inside `stage_tasks`, and
the three detailed frame regions are inside `frame_total`. Consequently, region
fractions must not be added together. The maximum for a region and maximum
total runtime can also occur on different ranks.

## Current Regions

| Region | Scope |
| --- | --- |
| `before_integrator` | The complete `before_timeintegrator` task list once per cycle. |
| `before_stage` | The `before_stagen` task list for every explicit integrator stage. |
| `stage_tasks` | The main `stagen` task list for every explicit integrator stage. |
| `after_stage` | The `after_stagen` task list for every explicit integrator stage. |
| `after_integrator` | The complete `after_timeintegrator` task list once per cycle. |
| `outputs` | Scheduled outputs during evolution and all final outputs. Initial outputs occur before the timing counters are reset and are excluded. |
| `amr` | Adaptive mesh refinement and load-balancing work. |
| `timestep` | The driver-level new-timestep calculation. |
| `trml_cooling` | The `simple_TRML_extended` cooling kernel and reduction. This region is nested inside `stage_tasks`. |
| `tracer_save_flux` | Saving RK-accumulated density fluxes required by Lagrangian Monte-Carlo tracers. This region is nested inside `stage_tasks`. |
| `frame_total` | Complete frame-tracker application. This region is nested inside `after_integrator`. |
| `frame_control` | Target measurement and controller update within `frame_total`. |
| `frame_boundary` | Frame-dependent physical-boundary and ghost-primitive refresh within `frame_total`. |
| `frame_timestep` | Timestep refresh following a changed boost within `frame_total`. |
| `particle_push` | Particle pusher work. |
| `particle_comm` | Particle ownership, count exchange, send, receive, and completion tasks. |
| `particle_adjust` | Particle position adjustment after mesh refinement. |
| `particle_seed` | Scheduled tracer-seeding check and any due seeding event. |

Regions that do not apply to a run are printed with zero seconds and zero
calls. The TRML cooling region is currently specific to
`simple_TRML_extended`; the driver and particle regions are code-wide.

## Recommended Benchmark Workflow

Use the timers to compare controlled, short runs:

1. Build the same optimized executable for every case.
2. Use the same nodes, MPI ranks, GPUs, mesh, MeshBlock layout, integrator, and
   output schedule.
3. Change one feature at a time, such as tracers enabled versus disabled.
4. Prefer a cycle limit such as `nlim = 100` to `500`. A cycle limit makes the
   amount of work directly comparable when physical timesteps differ.
5. Keep large data outputs disabled unless output performance is the quantity
   being tested.
6. Repeat short cases when run-to-run variability matters and compare medians.
7. After locating the bottleneck, disable synchronized timing and measure normal
   end-to-end throughput using AthenaK's `zone-cycles/cpu_second` report.

Extract the complete timing report with:

```bash
rg '^PERFORMANCE_(TIMING|REGION)' slurm-job.out
```

Print regions in descending order of maximum-rank time with:

```bash
awk '$1 == "PERFORMANCE_REGION" && $2 != "name" {
  printf "%-24s %12.6f s %8.3f%% %8d calls\n", $2, $3, 100*$5, $6
}' slurm-job.out | sort -k2,2nr
```

For MPI imbalance, compare `max_rank_seconds` with `mean_rank_seconds`. A region
whose maximum greatly exceeds its mean should be investigated by rank before
optimizing its local kernels.

## Tracer Slowdown Case Study

The diagnostics were introduced to investigate a TRML run with a
`512 x 512 x 1024` mesh, sixteen `256^3` MeshBlocks, sixteen GPUs on two
Frontier nodes, and 822 Lagrangian Monte-Carlo tracers. A matched no-tracer run
used the same hydro, cooling, frame tracking, decomposition, and cycle count.

Before the fix, `particle_seed` consumed 57.6 seconds of a 73.7-second,
120-cycle benchmark. `SeedTracersAtTime()` allocated host mirrors and copied the
complete fluid state from every GPU to the CPU on every cycle before checking
whether a seeding event was due. Almost every copy was unnecessary because the
configured particles had already been seeded.

The corrected implementation checks the schedule metadata first. When no event
is due, it refreshes the global particle counts and returns without allocating
host fluid arrays or copying the mesh state. Fluid data are copied only when an
actual seeding event requires cell eligibility and sampling weights.

The final validation produced:

| Case | Total zone-cycles/s | Zone-cycles/s/node | Approximate wall time per simulation-time unit |
| --- | ---: | ---: | ---: |
| Tracers, old seeding path | `4.37e8` | `2.19e8` | `29.6 min` |
| Tracers, corrected seeding path | `1.99e9` | `9.94e8` | `6.5 min` |
| No tracers | `2.10e9` | `1.05e9` | `6.2 min` |

After the fix, `particle_seed` used 0.0018 seconds over 120 cycles. The corrected
tracer run was 5.7 percent slower than the no-tracer run in the complete
benchmark, including an additional final particle-history output. The
evolution-only tracer overhead was approximately 2 to 3 percent and was mainly
the saved density flux and particle communication, not movement of the 822
particles themselves.

## Adding a Region

New regions are declared in `PerformanceRegion` in `src/driver/driver.hpp`. A
matching name must be added in the same order to the report-name array in
`Driver::ReportPerformanceTiming()` in `src/driver/driver.cpp`.

Instrument an operation with:

```cpp
pdrive->StartPerformanceRegion(PerformanceRegion::my_region);
DoWork();
pdrive->StopPerformanceRegion(PerformanceRegion::my_region);
```

Every control-flow path after `StartPerformanceRegion()` must reach the matching
stop call. Do not recursively or concurrently enter the same region on one
rank: each region currently owns one timer rather than a timer stack. Prefer a
small number of scientifically meaningful coarse regions over fencing every
kernel.

After adding a region:

1. Build at least one Release configuration with the relevant backend.
2. Run a short case with timing disabled to check the normal path.
3. Run with timing enabled and verify the region call count and nesting.
4. Compare the sum only among mutually exclusive sibling regions; never sum a
   parent with its children.
