# Task 1 final regression followup

Final broad tests exposed four expectations tied to the old timestep calculation.
The adaptations below change only six test files. Both physical AMR input decks
are byte-identical to the frozen final source; their hashes are retained in the
[input audit](evidence/task1-input-immutability.json). No production C++ or
scientific acceptance tolerance changed.

## AMR event sampling

The unchanged 3D fixture switches its refinement threshold at time 2.55e-4.
Its frozen initial cycle dt was 1.7221175563e-4, with smaller refined-grid steps,
so derefinement occurred after cycle three. With the new LF bound, the ordinary
advective dt of 7.4992842625e-4 steps past that threshold before the first AMR
check; the original final test therefore exercises no refinement at all. The 2D
case likewise steps past its unchanged 3.5e-4 switch.

The existing `sts_max_dt_ratio` controls only a cycle cap. Initial budget probes
show a 3D STS budget of 1.3127838454e-3. A cap of 0.13 restores the intended event
sequence with times 0, 1.7066189991e-4, 2.0958526452e-4 and 2.5769884706e-4.
The 2D cap is 0.1, giving first dt 2.7262521112e-4. These are explicit flags in
2D/3D GPU churn, the CPU wall regression and the two-node 16-rank churn test;
the fixture fields, amplitudes, physical thresholds, refinement rules, mesh and
strictness remain unchanged.

CPU probes and the original assertions confirm 2D cell counts 576→2304→576 and
3D counts 13824→110592→13824. The wall test again has 216 blocks immediately
before cycle-three derefinement and 27 afterward. Its original initial margin
exceeds 0.03; the pre-derefinement minimum is 0.06035 (>0.05), the final minimum is
5.3154e-7, and exactly eight cells meet the unchanged 2e-6 wall tolerance. All
six LF health counters are zero. Both complete churn tests retain and pass their
original conservation, repair, div-B, sampling-count and topology assertions.

## Heated post-sweep stage count

For the uniform aligned fixture,
`chi=sqrt(8/pi)/(2*pi)` and `dt_FE=h²/(2 chi)`. Its cycle cap gives
`dt=20*0.9*dt_FE=0.008651519135223495`. Heating at 3000 changes the pressure to
`1+(2/3)*3000*dt=18.30303827044699`, so the post-half ratio is
`10*sqrt(p)=42.78205029033437`. The RKL2 capacity `(s²+s-2)/4` requires 13 stages.
The unheated case needs seven. The CPU helper already derives these values; the
MPI test's stale 11 was corrected to 13, retaining one/four-rank state/timestep
agreement and the stage total 64×(7+13)=1280.

## Zero-field operator

When every face is below `bfloor`, the discrete LF row is zero and the controller
correctly schedules zero LF stages. The stale positive-stage assertion now
requires zero. The pgen's independent pressure-state check at absolute tolerance
1e-13, zero face/work assertions and all repair checks are unchanged. The retained
initial failure confirms this was solely `0 > 0`, not a physical-state failure.

## Resolved quantitative decay

The larger STS budget reduced both parallel quantitative fixtures to 93 cycles,
below their unchanged requirement of at least 100. Their original application
decay validations passed. At fixed aligned background coefficients,
`dt_FE=dx²/(2 chi_parallel)`, so the old test cap of ten gives only
`2*128²/(0.9*10*(2*pi)²)=92.22` full-step equivalents in one parallel damping
time. A test-only cap of eight raises this to 115.28. All four quantitative cases
now pass, including the original `cycles >= 100` assertion, one-damping-time
check, 0.003 relative amplitude tolerance, and deliberately incorrect closure
negative references. Their physical input decks remain byte-identical; the
[audit](evidence/decay-cap-audit.json) records input hashes and measured errors.
The four relative errors range from 2.029e-4 to 2.761e-4.

## Validation status

The exact final CPU binary passed five focused tests (wall, low field, two
single-rank heating cases, and one/four-rank heating agreement), plus both original
churn tests and all four quantitative decay cases. The corrected MPI heating check also passed independently before
that combined run. Evidence includes original failures, commands, probe histories,
[test patch](evidence/task1-final-regression.patch), logs and JUnit XML in the
[manifest](evidence-manifest.json).

The exact final HIP binary also passed all 11 focused checks: both capped
2D/3D churn tests, the wall test, low field, two single-rank heating cases,
one/four-rank heating agreement, and four quantitative decay cases. The separate
two-node 16-GPU churn test passed its original topology, conservation, div-B and
repair assertions. These runs used the final `run_pytest.sh` runtime in allocation
5629018, with binary SHA `ce444cd1c34992fbb9cb8e1b1ebb1d9ec569b4d0a93f1f969f7470ca395ca027`.
The [GPU audit](evidence/task1-final-gpu-followup.json) lists each passed case and
its JUnit evidence. GPU execution is confirmed independently of the earlier CPU
runs of GPU-named tests. All these checks are functional; concurrent tests make
their elapsed times unsuitable for performance claims.
