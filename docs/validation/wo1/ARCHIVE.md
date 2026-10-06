# WO1 evidence archive

The top-level run directories contain the original September 29, 2026 measurements, before the subsequent
review fixes. They are evidence for that code state, not new measurements of
later commits. Historical absolute paths and commit IDs are preserved verbatim.

- `P1-tested-rejected.patch` is the complete rejected candidate.
- `inputs/` and `inputs.json` preserve the fixed input bytes and hashes.
- `runs/post-g/provenance.json` records the baseline source, compiler and executable
  hashes. Each later run directory preserves its own provenance and source patch.
- `runs/*/results.json.gz` contains the original full JSON results compressed
  losslessly: output hashes, repeatability, comparison results and timing samples.
- `runs/p1/difference-analysis.json` and `failure-summary.txt` record the active-cell
  SMR/outflow discrepancy, including the first differing energy-history cycle.
- `harness.py`, `report.py`, and `README.txt` preserve the original scratch harness
  and its timing/comparison definitions. They contain historical local paths;
  reconstruct that layout or adapt those paths in a working copy before rerunning.
  Decompress the result JSON files if using the original report script.
- `pre-rebase-log.txt`, `final-commits.txt`, and `rebase-verification.txt` record
  the original and cleaned histories and exact final-tree equality.
- `SHA256SUMS.json` hashes every copied evidence file. It intentionally excludes
  this explanatory page and the checksum manifest itself.

Raw simulation outputs and executable copies remain in `/tmp/cgl-wo1-perf/`.
This compact committed archive preserves the requested key evidence independently
of that temporary directory; it does not contain all raw binary outputs.

`review/` contains the subsequent shear/AMR diagnoses, local test logs, final
binary hashes and a fresh 24-configuration, three-repeat output comparison.
Its copy of the comparison harness recognizes the added strict boundary checks
without counting them as LF stages. Timing samples from this concurrent
validation are not a revised performance claim. Weak-field failure and full
float-build blocker evidence are retained separately in that directory.

`review/weak-field-transition/` preserves the October 1 B4 stage trace, failed
bounded-reference experiment, and scratch physical-log-ratio prototype. Its
[report](review/weak-field-transition/STATUS_UPDATE.md) preserves the historical
investigation. The October 6 [reference comparison and acceptance decision](review/weak-field/README.md)
supersedes its pending-redesign status: retain A and accept the shared extreme
sharp-contact limitation, with the strict expected failure still visible.

`review/weak-field-reference/` contains the October 6 reference comparison and
compact reproduction evidence. `review/closeout/` records final local B4 checks
and the successful Linux CPU/MPI/CUDA compile run at production commit `7b3345fd`.

P1 is an unresolved correctness investigation for WO2. Compact exchange and
skipped magnetic boundary work must be isolated to determine which frozen ghost
values affect the coarse/fine outflow update. The existing comparison rejects
the candidate; it does not establish which implementation caused the difference.
