# Final GPU workflow runners

Prepared without launching applications. Run only after the final source/binary
and a coordinated allocation are assigned. The shell wrapper establishes the
same corrected Frontier environment used by accepted WO2 HIP runs.

```
SLURM_JOB_ID=ALLOCATION bash run_final_gpu_workflows.sh \
  --source /lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/final-research/source \
  --binary /lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/final-research/build-hip/src/athena \
  --build-dir /lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/final-research/build-hip \
  --expected-binary-sha256 FINAL_SHA256 \
  --git-metadata-source /autofs/nccs-svm1_home2/dfielding/athenak-cgl-wo2 \
  --run-name acceptance-hip --execute
```

All three paths are explicit; source, build and binary must be below WO2.
The workflow script must reside in the supplied source tree (not resolve through
a symlink to an external source). Outputs, cwd, temporary files and caches remain
under final-research. Omit `--execute` for a preparation-only manifest, and use a
separate `--run-name` for the subsequent executed run; retained outputs are never
overwritten.

Both commands pass `--no-build`. The launcher hashes the exact provided binary
before every direct or MPI launch and uses `srun --exact` without overlap. It
records actual launch commands and allocation. Full acceptance requires all 29
cases; paper-smoke requires all declared cases (three after passive release),
with no disabled case. No scientific overrides are supplied. The runner clears
inherited LF/profiling environment controls, retains source/input/binary/compiler
cache hashes, and checks for changes after the workflows complete.

Python syntax and shell syntax checks passed. The immutable snapshot has no usable
Git worktree metadata. `--git-metadata-source` supplies an existing metadata
directory with `GIT_WORK_TREE` set to the frozen source and optional locks disabled.
The relocated Kokkos `.git` pointer is ignored only for Git status; this exception
is explicit in the runner manifest and does not change the binary, source, or
input hashes. Two initial metadata failures are retained and launched no apps.

The current executed run is `acceptance-hip-final2` on allocation 5629018, against
`source-fused` and the exact HIP binary SHA256
`ce444cd1c34992fbb9cb8e1b1ebb1d9ec569b4d0a93f1f969f7470ca395ca027`.
The active acceptance is functional only, concurrent with other final regression
suites; elapsed times are not performance claims.
