# Final fused GPU workflows

The exact final HIP binary passed all **29 full-workflow cases** and all **three
paper-smoke cases**, including passive CGL. No case was disabled or skipped.
These are functional acceptance runs, concurrent with other final suites in
allocation 5629018; elapsed times are not performance measurements.

Binary SHA256:
`ce444cd1c34992fbb9cb8e1b1ebb1d9ec569b4d0a93f1f969f7470ca395ca027`.
The immutable source was `WO2/final-research/source-fused`. Both workflow commands
used `--no-build` and the retained exact-binary launcher. Each launch checked the
binary hash, and the final runner confirmed unchanged binary, workflow, source
files and input decks. No scientific overrides were added.

The [full manifest](evidence/full/manifest.json) retains all scientific diagnostics
and case commands. The [paper-smoke manifest](evidence/paper-smoke/manifest.json)
contains the active Alfvénic, active random and passive Alfvénic cases. The
[runner manifest](evidence/runner-manifest.json) records source/input/binary hashes,
compiler-cache provenance and the expected-versus-executed case lists.

Two initial attempts failed before any application launch because the immutable
snapshot has no Git worktree metadata and its copied Kokkos submodule pointer is
relocated. Their logs are retained. The final runner supplies read-only production
Git metadata with the frozen source as its worktree and optional locks disabled;
only the unusable Kokkos submodule Git-status lookup is ignored. This affects
metadata collection, not source files, solver code or scientific acceptance.

The [evidence manifest](evidence-manifest.json) hashes every archived command,
input, log, history and summary. Larger original case outputs remain under
`WO2/final-research/acceptance-hip-final2`. This report covers these workflows;
the broader CPU/HIP/MPI suites and their separately documented test adaptations
are independent final acceptance evidence.
