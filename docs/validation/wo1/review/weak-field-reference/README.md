# Compact B4 reference-comparison evidence

Decision: retain conservative A and the current implementation. On 2026-10-06 the user accepted the shared extreme sharp-contact limitation. The strict expected failure and numerical bounds remain unchanged. See [the report](STATUS_UPDATE.md).

This package copies the October 6 comparison's reports, figures, numerical summaries, commands, build logs, patches, source/binary hashes, and scripts. The dated decision update supersedes the original report's pending-decision prose; historical measurements are unchanged. Family manifests preserve the original scratch inventories and are not inventories of this compact package.

`final-smooth-states.tar.gz` contains only the 12 final smooth-state tables required by the existing analysis script: six current-code cases and three cases for each reference. [The subset inventory](final-smooth-states.json) records their original bytes and hashes. It contains no complete time series or full stage traces. The current summary supplies the sharp-contact histories and energy-budget records used by the script. It is stored losslessly as `evidence/current/summary.json.gz`; decompress it before plotting or replaying the analysis.

From this directory, reproduce the figures and derived metrics using Python with NumPy and Matplotlib:

```sh
gzip -dc evidence/current/summary.json.gz > evidence/current/summary.json
tar -xzf final-smooth-states.tar.gz -C evidence
python3 analyze_comparison.py --runs-root evidence --output figures
```

The full local bundle remains at `/Users/dbf75/Documents/Codex/2026-10-06/cgl-b4-comparison`. Its three `evidence/{current,majeski,squire}/runs.tar.gz` archives (108.5 MB combined) are deliberately omitted here. [FULL_BUNDLE_SHA256SUMS.json](FULL_BUNDLE_SHA256SUMS.json) is the unchanged inventory of that original bundle, including its original report and omitted raw archive hashes. `SHA256SUMS.json` instead inventories this compact package. Absolute paths within original provenance, commands, and scripts are historical records; adapt working copies when replaying elsewhere. The relative links and final-state subset make the present findings reviewable without those paths.

Source rebuilding requires the pinned source revisions and dependencies recorded by each family. The Squire development snapshot is not a verified publication revision, and its disabled RMS accumulation makes both nominal floor settings duplicate all-magnetized controls. No numerical behavior from that archive was silently repaired.
