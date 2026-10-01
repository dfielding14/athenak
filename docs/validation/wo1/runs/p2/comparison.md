# p2 versus post-g

24 configurations, 3 repeats; 0 byte-comparison failures.

| Case | Solver ms/cycle before | After | After/before | Partial compute µs/stage before | After | After/before |
|---|---:|---:|---:|---:|---:|---:|
| lf1d-serial-0 | 0.51253 | 0.489406 | 0.9549 | 30.2 | 27.232 | 0.9017 |
| lf2d-serial-0 | 17.9485 | 16.8225 | 0.9373 | 1175.4 | 1093.94 | 0.9307 |
| lf3d-serial-0 | 222.032 | 211.354 | 0.9519 | 14492 | 13710.5 | 0.9461 |
| smr-serial-0 | 17.6533 | 16.2939 | 0.9230 | 927.783 | 825.396 | 0.8896 |
| pure_cgl-serial-0 | 0.414306 | 0.448139 | 1.0817 | — | — | — |
| lf2d-mpi-4 | 5.31476 | 4.92957 | 0.9275 | 311.004 | 283.743 | 0.9123 |
| lf3d-mpi-4 | 58.974 | 55.7028 | 0.9445 | 3739.98 | 3499.44 | 0.9357 |
| smr-mpi-4 | 5.45213 | 5.04405 | 0.9252 | 245.654 | 215.159 | 0.8759 |

Partial compute per stage sums existing exclusive LF/STS buckets (rank mean for MPI). Whole-STS wall time is not separately available. Shared transport includes RK/initialization and is excluded here; all samples and timing categories remain in results.json.

Baseline JSON: /private/tmp/cgl-wo1-perf/runs/post-g/results.json
Candidate JSON: /private/tmp/cgl-wo1-perf/runs/p2/results.json
