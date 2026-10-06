# Final release equivalence

All 25 CPU and 25 HIP accepted/release pairs passed in allocation 5629018: 100 successful applications, no numerical mismatches. Each backend includes 20 one-rank pairs and five four-rank pairs. The accepted fused binaries were preserved. The release binaries use the final source snapshot, including formatting, two complete-expression `Real` casts, and the literal pressure-traction helper relocation.

| Exact file comparisons per backend | Pairs |
| --- | ---: |
| Field binary files | 50 |
| History files, including forcing history | 27 |
| Normalized restart files | 50 |

Field and history files match byte for byte. Restart comparisons include every diagnostic and live RNG member, without a floating-point tolerance. The existing explicit normalization covers only 36 unused root-index bytes for ordinary cases; forced cases additionally enable the validated dormant seed1 startup RNG members and four native RNG padding bytes. Initial forced restart normalization totals 320 bytes, evolved forced normalization totals 40 bytes. Exact ranges and both raw hashes remain in the result files. There are 24 CPU and 26 HIP raw restart differences; all disappear under those existing rules. No additional bytes were excluded.

The saved 23 Task7 inputs cover active and passive physics, safe/fast arithmetic, full/no diagnostics, dimensionality and refinement. Two extra forced passive 3D safe/full cases use one and four ranks for 12 cycles, each accumulating 491520 LF cell-stages. Both binaries see the identical staged seed1 input; the fixture manifest records its derivation from the earlier seed519 forcing setup. All six recorded repair counters are zero in every accepted and release application. All passive fixtures use PLM. These are ordinary-run state and emitted-restart comparisons, not resumed-run validation, new face-array coverage, or performance measurements.

The HIP runner enforces the complete previously accepted GPU-aware MPI/HSA/FI environment, recorded in each result. The CPU runner uses CPU-only modules and disables GPU-aware MPI. Binary SHA-256 values were checked before each launch and after completion.

The initial CPU comparison attempt is retained separately. Thirteen accepted applications succeeded while the initial release processes exited at dynamic loading because that CPU link inherited a ROCm dependency. No release solver initialized in that attempt. Root relinked unchanged CPU objects under CPU-only modules, preserved the original binary and logs, and recorded the command/environment/hashes. All 25 CPU pairs were then rerun in a fresh directory and passed.

Evidence under this directory:

- `release-equivalence-summary.json`: compact results, per-case input hashes, source/build manifest hashes, runtime environment, and raw-different restart metadata.
- `release-equivalence-{cpu,hip}-ready-5629018/results.json`: all 50 pair records, logs, commands, file hashes, and normalization ranges.
- `release-bin/manifest-{cpu,hip}.json` and `release-source-manifest.json`: immutable release provenance.
- `run_release_equivalence.py`, `run_release_equivalence.sh`, and `summarize_release_equivalence.py`: retained runner and independent completion summary.
- `release-equivalence-cpu-5629018/`: retained initial loader failure, including `loader-failure-audit.json`.

Release CPU SHA-256: `671441a4bbfb62f7192b7b6aaeaf8cc522603a84ea0482cab7f4d0408ab12d46`.

Release HIP SHA-256: `6530c1b31a7d99077f3ff5932d4b1c790e33be38c1fa222cfda097f87b339f9e`.
