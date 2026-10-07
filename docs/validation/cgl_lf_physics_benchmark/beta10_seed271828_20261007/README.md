# Beta-ten turbulence reference, seed 271828

This is the first scientifically reviewed reference for the reusable
[CGL-LF physics benchmark](../../../cgl_lf_physics_benchmark.md), completed
2026-10-07 UTC on Frontier. The corrected realization reaches t=18 at cycle
32604. All physics parameters retain the canonical input, including
`dedt=0.32` per volume (nominal total power 0.64), with no cooling.

Read the [report and four figure groups](report.md), the
[scientific review](scientific-review.md), and the
[numerical correction/checks](numerical-checks.md). Numerical integrity is
consistent in the corrected run. Several physical signatures are consistent
with the intended behavior, but quantitative paper agreement remains
inconclusive. This is a descriptive reference, not an acceptance oracle.

## Measurements and sampling

The fixed averaging interval is [6,18]: 50 bracketing snapshots and six
duration-two blocks. Values below are time means ± sample SD of block means,
not independent-sample confidence intervals.

| Measurement | Mean ± block SD |
| --- | ---: |
| Isotropic sound-speed Mach proxy | 0.21366 ± 0.01626 |
| Vector delta-B RMS / B0 | 0.80664 ± 0.07136 |
| Volume mean of local beta | 9.2424 ± 1.4971 |
| Signed perpendicular-pressure correlation C | −0.15199 ± 0.13213 |
| Normalized pressure residual variance R | 0.87297 ± 0.10487 |
| Mirror near-band fraction, halfwidth 0.05 | 0.03288 ± 0.01943 |
| Firehose near-band fraction, halfwidth 0.05 | 0.01140 ± 0.00956 |

Resolved parallel-strain/perpendicular-gradient power is 0.004826; substantial
local parallel velocity fluctuations remain. The spectrum near the numerical
cutoff has a different ratio and is not convergence evidence. Global pressure
compensation is modest and time dependent. A bounded
[signed pressure-band audit](pressure-episode/README.md) explains why the
positive global correlation around t=14 does not imply loss of compensation
at every resolved scale. That audit contains three instantaneous snapshots,
not another temporal ensemble.

The [predeclared sampling plan](experiment-plan.json) and
[second extension decision](fixed-extension-to18-decision.json) retain the
reason for extending the same run from 10 to 14 to 18, keeping the averaging
start at six. Earlier reports remain under the raw-data root. Two-unit blocks
are correlated; the review also pairs them into three four-unit blocks.
Heating and beta drift prevent a stationary thermal interpretation. No further
extension was chosen merely to seek agreement with a paper.

## Reproduction and provenance

The [guide](../../../cgl_lf_physics_benchmark.md) gives the exact fresh build,
run, resume and analysis procedure. Build with `PROBLEM=built_in_pgens` and
select `problem/pgen_name=cgl_lf_paper`; the legacy `PROBLEM=cgl_lf_paper`
route is a different pgen. The built-in pgen itself is unchanged.

For the retained data, after loading the guide's Python environment:

```bash
BENCH_SOURCE=/autofs/nccs-svm1_home2/dfielding/athenak-cgl-wo2
BENCH_ROOT=/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/physics-benchmark
export PYTHONDONTWRITEBYTECODE=1
export MPLCONFIGDIR="$BENCH_ROOT/mpl-cache" TMPDIR="$BENCH_ROOT/tmp"
python3 "$BENCH_SOURCE/scripts/analyze_cgl_lf_physics_benchmark.py" \
  "$BENCH_ROOT/analysis-run-fixed-to18" \
  --time-start 6 --time-end 18 --block-duration 2 \
  --output-dir "$BENCH_ROOT/analysis-fixed-6-18-reproduced"
```

The aggregate lists `canonical-fixed`, `canonical-fixed-to14` and
`canonical-fixed-to18`, in that order. Their exact launch commands, restart
hashes, overrides, input, source/build identity, backend and module environment
are embedded in [metrics.json](metrics.json). Simulation revision is
`71ad25ebce73d33db048defd8f585a7dce0528c2`; the analysis checkout revision is
`64d652317bec173d1d33e23906d0e93f21463409`. The corrected executable SHA256 is
`bd699cf711d0c2a21037690d28e8639b140f0ed96ad19c2d0406636c331527a3`.
The original numerical-analysis script SHA256 is
`a8241542d0d3cc5e4ecc83dd17b4cd7819467f9c2079eb707ebf231a5df812aa`.

[Publication provenance](provenance.json) pins the original generated
artifacts. Metrics remain unchanged. The figures were subsequently rendered
from those retained metrics using `dbfplot`, with logarithmic PDFs,
mathematical definitions and illustrative power-law guides. The PDF bands
show the pointwise min–max of six duration-two averaged curves, not confidence
intervals. [figure-rendering.json](figure-rendering.json) records the separate
plotting revision, source hash, exact command and prior rendering record;
[figure-audit.json](figure-audit.json) records the strict export audits.
All four groups use PDF and 240-dpi PNG; only the multi-panel canvas size is
exempted from the default single-panel audit profile.
The published report adds a labeled scientific-review preface and trims
trailing whitespace.
[Artifact hashes](artifact-hashes.json) cover the final retained files.
Large simulation binaries/restarts remain outside Git under `BENCH_ROOT`
(about 9.4 GiB for the three corrected run directories). Original failed
restart evidence is preserved separately there.

## Measured cost and limits

The three corrected segments used **4465.65 s = 74.43 min** of launcher wall
time on eight GPU compute devices, **9.9237 device-hours**, on one Frontier
node (four dual-device MI250X packages). The final CPU analysis took
**224.86 s**, with peak RSS **1,151,828 KiB**. The original t=10 discovery run
and failed continuation used another **5.2428 device-hours**; they are not
part of the corrected reference cost. [costs.json](costs.json) retains
per-segment costs, plumbing/analysis timings and separate Slurm allocation
wall accounting. Application device-hours exclude compilation, allocation
idle time and CPU analysis; no unmeasured build cost is invented.
Both allocations completed and were released: their combined allocated wall
time is **2 h 28 min 20 s (2.4722 node-hours)**, including the discovery run,
investigations, analysis and idle intervals.

Actual applied work over [6,18] is 7.68; conserved-energy change minus work is
−5.72e−12. Floor/nonfinite/nonpositive counters and post-operator hard-volume
history remain zero. The numerical-checks note distinguishes the reserved
`lf_hwproj` counter from measured activity and records the corrected real
restart tests. A successful analysis does not establish physical validation.

One active box cannot establish an active/passive causal contrast, spatial
convergence, a precise inertial-range slope or LF coefficient accuracy. Keep
the independent physical damping tests and the accepted B4 sharp-contact
limitations. For future changes, first check definitions/sampling, then
forcing and achieved regime, then resolution, then implementation.
