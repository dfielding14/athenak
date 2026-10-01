# p1 versus post-g

24 configurations, 3 repeats; 9 byte-comparison failures.

| Case | Solver ms/cycle before | After | After/before | Partial compute µs/stage before | After | After/before |
|---|---:|---:|---:|---:|---:|---:|
| lf1d-serial-0 | 0.51253 | 0.522001 | 1.0185 | 30.2 | 31.167 | 1.0320 |
| lf2d-serial-0 | 17.9485 | 18.2275 | 1.0155 | 1175.4 | 1199.71 | 1.0207 |
| lf3d-serial-0 | 222.032 | 223.315 | 1.0058 | 14492 | 14702 | 1.0145 |
| smr-serial-0 | 17.6533 | 17.1706 | 0.9727 | 927.783 | 939.888 | 1.0130 |
| pure_cgl-serial-0 | 0.414306 | 0.438819 | 1.0592 | — | — | — |
| lf2d-mpi-4 | 5.31476 | 5.29176 | 0.9957 | 311.004 | 311.314 | 1.0010 |
| lf3d-mpi-4 | 58.974 | 58.7645 | 0.9964 | 3739.98 | 3763.08 | 1.0062 |
| smr-mpi-4 | 5.45213 | 5.26446 | 0.9656 | 245.654 | 246.226 | 1.0023 |

Partial compute per stage sums existing exclusive LF/STS buckets (rank mean for MPI). Whole-STS wall time is not separately available. Shared transport includes RK/initialization and is excluded here; all samples and timing categories remain in results.json.

Baseline JSON: /private/tmp/cgl-wo1-perf/runs/post-g/results.json
Candidate JSON: /private/tmp/cgl-wo1-perf/runs/p1/results.json

Exact failed files:

```json
[
  {
    "case": "smr_outflow-serial-0",
    "repeat": 0,
    "type": "baseline mismatch",
    "files": [
      "bin/smr_outflow.mhd_u_bcc.00001.bin",
      "bin/smr_outflow.mhd_w_bcc.00001.bin",
      "rst/smr_outflow.00001.rst",
      "smr_outflow.mhd.hst",
      "smr_outflow.user.hst"
    ]
  },
  {
    "case": "smr_outflow-serial-0",
    "repeat": 1,
    "type": "baseline mismatch",
    "files": [
      "bin/smr_outflow.mhd_u_bcc.00001.bin",
      "bin/smr_outflow.mhd_w_bcc.00001.bin",
      "rst/smr_outflow.00001.rst",
      "smr_outflow.mhd.hst",
      "smr_outflow.user.hst"
    ]
  },
  {
    "case": "smr_outflow-serial-0",
    "repeat": 2,
    "type": "baseline mismatch",
    "files": [
      "bin/smr_outflow.mhd_u_bcc.00001.bin",
      "bin/smr_outflow.mhd_w_bcc.00001.bin",
      "rst/smr_outflow.00001.rst",
      "smr_outflow.mhd.hst",
      "smr_outflow.user.hst"
    ]
  },
  {
    "case": "smr_outflow-mpi-1",
    "repeat": 0,
    "type": "baseline mismatch",
    "files": [
      "bin/smr_outflow.mhd_u_bcc.00001.bin",
      "bin/smr_outflow.mhd_w_bcc.00001.bin",
      "rst/smr_outflow.00001.rst",
      "smr_outflow.mhd.hst",
      "smr_outflow.user.hst"
    ]
  },
  {
    "case": "smr_outflow-mpi-1",
    "repeat": 1,
    "type": "baseline mismatch",
    "files": [
      "bin/smr_outflow.mhd_u_bcc.00001.bin",
      "bin/smr_outflow.mhd_w_bcc.00001.bin",
      "rst/smr_outflow.00001.rst",
      "smr_outflow.mhd.hst",
      "smr_outflow.user.hst"
    ]
  },
  {
    "case": "smr_outflow-mpi-1",
    "repeat": 2,
    "type": "baseline mismatch",
    "files": [
      "bin/smr_outflow.mhd_u_bcc.00001.bin",
      "bin/smr_outflow.mhd_w_bcc.00001.bin",
      "rst/smr_outflow.00001.rst",
      "smr_outflow.mhd.hst",
      "smr_outflow.user.hst"
    ]
  },
  {
    "case": "smr_outflow-mpi-4",
    "repeat": 0,
    "type": "baseline mismatch",
    "files": [
      "bin/smr_outflow.mhd_u_bcc.00001.bin",
      "bin/smr_outflow.mhd_w_bcc.00001.bin",
      "rst/smr_outflow.00001.rst",
      "smr_outflow.mhd.hst",
      "smr_outflow.user.hst"
    ]
  },
  {
    "case": "smr_outflow-mpi-4",
    "repeat": 1,
    "type": "baseline mismatch",
    "files": [
      "bin/smr_outflow.mhd_u_bcc.00001.bin",
      "bin/smr_outflow.mhd_w_bcc.00001.bin",
      "rst/smr_outflow.00001.rst",
      "smr_outflow.mhd.hst",
      "smr_outflow.user.hst"
    ]
  },
  {
    "case": "smr_outflow-mpi-4",
    "repeat": 2,
    "type": "baseline mismatch",
    "files": [
      "bin/smr_outflow.mhd_u_bcc.00001.bin",
      "bin/smr_outflow.mhd_w_bcc.00001.bin",
      "rst/smr_outflow.00001.rst",
      "smr_outflow.mhd.hst",
      "smr_outflow.user.hst"
    ]
  }
]
```
