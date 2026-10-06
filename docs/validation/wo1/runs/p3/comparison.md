# p3 versus post-g

24 configurations, 3 repeats; 0 byte-comparison failures.

| Case | Solver ms/cycle before | After | After/before | Partial compute µs/stage before | After | After/before |
|---|---:|---:|---:|---:|---:|---:|
| lf1d-serial-0 | 0.51253 | 0.474336 | 0.9255 | 30.2 | 25.8468 | 0.8559 |
| lf2d-serial-0 | 17.9485 | 16.4182 | 0.9147 | 1175.4 | 1061.74 | 0.9033 |
| lf3d-serial-0 | 222.032 | 204.298 | 0.9201 | 14492 | 13244.3 | 0.9139 |
| smr-serial-0 | 17.6533 | 16.1661 | 0.9158 | 927.783 | 819.099 | 0.8829 |
| pure_cgl-serial-0 | 0.414306 | 0.435514 | 1.0512 | — | — | — |
| lf2d-mpi-4 | 5.31476 | 4.82113 | 0.9071 | 311.004 | 277.251 | 0.8915 |
| lf3d-mpi-4 | 58.974 | 55.0435 | 0.9334 | 3739.98 | 3438.47 | 0.9194 |
| smr-mpi-4 | 5.45213 | 5.11799 | 0.9387 | 245.654 | 217.781 | 0.8865 |

Partial compute per stage sums existing exclusive LF/STS buckets (rank mean for MPI). Whole-STS wall time is not separately available. Shared transport includes RK/initialization and is excluded here; all samples and timing categories remain in results.json.

Baseline JSON: /private/tmp/cgl-wo1-perf/runs/post-g/results.json
Candidate JSON: /private/tmp/cgl-wo1-perf/runs/p3/results.json
