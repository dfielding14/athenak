# p6 versus post-g

24 configurations, 3 repeats; 0 byte-comparison failures.

| Case | Solver ms/cycle before | After | After/before | Partial compute µs/stage before | After | After/before |
|---|---:|---:|---:|---:|---:|---:|
| lf1d-serial-0 | 0.51253 | 0.46879 | 0.9147 | 30.2 | 25.4872 | 0.8439 |
| lf2d-serial-0 | 17.9485 | 16.2852 | 0.9073 | 1175.4 | 1052.44 | 0.8954 |
| lf3d-serial-0 | 222.032 | 202.263 | 0.9110 | 14492 | 13095.8 | 0.9037 |
| smr-serial-0 | 17.6533 | 16.0113 | 0.9070 | 927.783 | 809.41 | 0.8724 |
| pure_cgl-serial-0 | 0.414306 | 0.421236 | 1.0167 | — | — | — |
| lf2d-mpi-4 | 5.31476 | 4.82474 | 0.9078 | 311.004 | 276.791 | 0.8900 |
| lf3d-mpi-4 | 58.974 | 54.86 | 0.9302 | 3739.98 | 3427.13 | 0.9164 |
| smr-mpi-4 | 5.45213 | 5.10856 | 0.9370 | 245.654 | 215.664 | 0.8779 |

Partial compute per stage sums existing exclusive LF/STS buckets (rank mean for MPI). Whole-STS wall time is not separately available. Shared transport includes RK/initialization and is excluded here; all samples and timing categories remain in results.json.

Baseline JSON: /private/tmp/cgl-wo1-perf/runs/post-g/results.json
Candidate JSON: /private/tmp/cgl-wo1-perf/runs/p6/results.json
