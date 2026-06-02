# MKS24 Reproduction-First CGL-LF Manuscript and Later Extension Plan

Status: execution plan, revised after accepted corrected-E03 R02 job `4754008`
reached exact `t = 7.0`, the corrected F-100 recost was retained, and F-083/F-084
controller hardening remained promoted. A new committed post-F-100 documentation
archive remains required before successor preparation. The preceding read-only
Frontier-root audit was completed on 2026-05-29 EDT (2026-05-30 UTC). The
first paper-production
segment, mapped case `R16/s00_t0_t2`, completed on Frontier as job `4674731`
on 2026-05-25 and was recorded only as a clean partial diagnostic result
(`t = 1.66204848945`, `1.253333` node-hours), not an accepted scientific
result. After that legacy shared-MPI-I/O segment incurred an approximately
twenty-minute pause at its `t = 1` snapshot/checkpoint boundary, commit
`1f52784d` changed retained production snapshots and restarts to rank-local
files with grouped checksum/analysis support. A fresh ranked-output
`R16/s01_rankio_t0_t2` segment began running as Frontier job `4675322` on
2026-05-25; its complete initial eight-rank snapshot/restart groups were
successfully read by the rank-aware analyzer, and its complete `t = 1`
eight-rank snapshot/restart groups were emitted without the legacy-scale
boundary pause. It completed as an inspected clean partial at
`t = 1.87769834013` using `1.253056` node-hours. Same-digest continuation
`R16/s02_rankio_t1p877698_t2` (job `4676557`) reached `t = 2.0`, passed
inspection with zero strict LF safety counters and complete terminal
rank-local products, and used `0.106389` node-hours. It is an accepted prefix
gate, not a completed `t = 10` case. Same-digest continuation
`R16/s03_rankio_t2_t3p5` (job `4676696`) completed on 2026-05-25 at
16:50:22 EDT, reached exactly `t = 3.5`, passed formal inspection with zero
strict LF safety counters and complete terminal rank-local snapshot/restart
products, and used `1.383611` node-hours. It extends the accepted ranked
prefix only. Continuation `R16/s04_rankio_t3p5_t5` (job `4679803`)
completed on 2026-05-25 at 21:53:04 EDT, reached exactly `t = 5.0`, passed
formal inspection with zero strict LF safety counters and complete terminal
rank-local snapshot/restart products, and used `1.567778` node-hours.
Continuation `R16/s05_rankio_t5_t6p5` (job `4681476`) terminated cleanly
on its Athena wall-clock limit on 2026-05-26 at 00:11:57 EDT, reached
`t = 6.43828031751 < 6.5`, passed partial-continuation inspection with zero
strict LF safety counters and complete terminal ranked products, and used
`1.753333` node-hours. Short continuation `R16/s06_rankio_t6p438280_t6p5`
(job `4683291`) completed on 2026-05-26 at 00:34:30 EDT, reached exactly
`t = 6.5`, passed formal inspection with zero strict LF safety counters and
complete terminal ranked products, and used `0.080556` node-hours. It extends
the accepted ranked prefix to `t = 6.5`. Continuation
`R16/s07_rankio_t6p5_t7p5` (job `4683918`) completed on 2026-05-26 at
01:58:33 EDT, reached exactly `t = 7.5`, passed formal inspection with zero
strict LF safety counters and complete terminal ranked products, and used
`1.281944` node-hours. It extends the accepted ranked prefix to `t = 7.5`.
Continuation `R16/s08_rankio_t7p5_t8p5` (job `4686032`) completed on
2026-05-27 at exact `t = 8.5`, passed formal inspection with zero strict LF
safety counters and complete terminal ranked products, and used `1.282778`
node-hours. On 2026-05-29 it was recorded `rejected` for current reproduction
admission because the merged turbulent-driver replacement changes forcing
evolution and restart state/layout. The pre-replacement restart lineage
cannot continue under the replacement executable. At that E01 stop line,
Stage I actual use was `9.962778` node-hours with no active reservation. The
recovery decision was a new execution epoch, not an attempt to salvage the
pre-replacement lineage. Preserve the existing `R16` tree, inspections, and
ledger as immutable historical epoch `E01-pre-modal-driver`. The separate
`E02-modal-driver` namespace and ledger were created before the next MKS24
submission. That
namespace, its initial pilot-only reservation, and its shared-root submission
checks are committed through `7fae0bcf`. Fresh replacement-driver `E02`
jobs `4743735`, `4743933`, `4743988`, `4744019`, `4744056`, `4744120`,
and `4744158` are now accepted through exact `R16 t = 10.0`, using
`6.145556` node-hours. The completed lineage bundle passes the `t = 8`--`10`
production-window analyzer. Standard-layout `R02` jobs `4744198` and
`4744205` are accepted through native snapshot boundary `t = 0.25`, use
`0.473333` node-hours, and project `19.000000` node-hours per standard
`t = 10` case from the authenticated continuation rate. High-resolution `R17`
jobs `4744210` and `4744230` are accepted through `t = 0.10`, use `4.235556`
node-hours, and project `42.488889` node-hours per simulated time unit from
the authenticated continuation rate. F-071 closes the residual runtime and
storage bracket and authorizes only the frozen sixteen-case `R02`--`R17`
mapped matrix under sequential inspection inside a `900.000000` node-hour
`E02` Stage I envelope. R02 job `4744249` then advances the first mapped
production lineage through exact `t = 1.0` in `1.452778` node-hours. F-072
updates the continuation-aware matrix projection to `702.195370` node-hours,
leaving `197.804630` node-hours of envelope margin. R02 job `4744518`
continues the same inspected lineage through exact `t = 1.5` in `1.052222`
node-hours. F-073 updates the projection to `725.464938` node-hours, leaving
`174.535062` node-hours of envelope margin. R02 job `4744913` continues the
same inspected lineage through exact `t = 2.0` in `1.068889` node-hours.
F-074 updates the projection to `730.081697` node-hours, leaving `169.918303`
node-hours of envelope margin. R02 job `4745305` continues the same inspected
lineage through exact `t = 2.5` in `1.061944` node-hours. F-075 updates the
projection to `728.164877` node-hours, leaving `171.835123` node-hours of
envelope margin. F-075 historically authorized continuation of the frozen
matrix one inspected segment at a time with `R17` last. Continuation
`R02/s07_rankio_t2p5_t3` was prepared
from the inspected `s06` terminal restart siblings and submitted as job
`4745498` on 2026-05-30. An independent Phase A source-to-code audit then
found that the retained `E02` turbulent driver does not exactly implement the
published MKS24 forcing semantics: random roles inherit a generic projection
and implicit parabolic-spectrum default, while planar roles do not generally
enforce `grad_perp dot u_perp = 0` for retained `k_z != 0` modes. Job
`4745498` was cancelled after `498` seconds and recorded `aborted`, consuming
`0.138333` node-hours and bringing E02 cumulative use to `15.628610`.
Stop-line evidence JSON SHA-256 is
`de85c443c3bc9bab795f7f1903e5de224d7d86b553ac40e9dd7249986c3df211`.
Preserve all E02 products as pipeline and cost evidence, but prohibit every
new E02 preparation or submission. F-078 closes explicit-policy
implementation and corrected-build qualification. The reviewed E03 token is
retained, reconciliation passed, and fresh `R02/s00_rankio_t0_t0p1` job
`4745922` was submitted from `t = 0`, formally inspected, and recorded
`accepted` at exact `t = 0.1` for `0.188889` node-hours.
Authenticated continuation `R02/s01_rankio_t0p1_t0p25` job `4746154` is
also formally inspected and recorded `accepted` through exact `t = 0.25` for
`0.289167` node-hours. Corrected E03 Stage I use is now `0.478056`
node-hours. Retained E03 R02 recost evidence JSON SHA-256 is
`eb6071cf53d453b2f75707003616ac9cbdfe025d62969d00e34e7948bee3310c`;
its provisional corrected-matrix projection is `700.868889` node-hours
inside the `900.000000` envelope.
Authenticated `R02/s02_rankio_t0p25_t1` job `4746182` is also formally
inspected and recorded `accepted` through exact `t = 1.0` for `1.465556`
node-hours. Corrected E03 Stage I use is now `1.943612` node-hours. Updated
R02 recost evidence JSON SHA-256 is
`7d11e8a0004e24417167aeb1f186f8b1f68c8eb9d6206467a708a7cad16cf252`;
preserve it as superseded arithmetic history. The corrected conservative
projection at exact `t = 1.0` is `704.604815` node-hours inside the
`900.000000` envelope. Authenticated `R02/s03_rankio_t1_t1p5` job `4746356`
is also formally and independently inspected and recorded `accepted` through
exact `t = 1.5` for `1.048611` node-hours. Corrected E03 Stage I use is now
`2.992223` node-hours. Retained corrected R02 recost evidence JSON SHA-256 is
`b6fd6e8fcd939ca1baa06411a3cd6f04359b3bc2ff5b2a0097e2cc87b1763fc7`;
its authorization-facing corrected-matrix projection is `724.645556`
node-hours inside the `900.000000` envelope, leaving `175.354444`
node-hours of margin.
Authenticated `R02/s04_rankio_t1p5_t2` job `4746435` is also formally and
independently inspected and recorded `accepted` through exact `t = 2.0` for
`1.063611` node-hours. Corrected E03 Stage I use is now `4.055834`
node-hours. Retained corrected R02 recost evidence JSON SHA-256 is
`c382b78ab9e21466648adf9d7ea38d2b407f0e986579f128bdfe6bc44bef79dd`;
its authorization-facing corrected-matrix projection is `728.845556`
node-hours inside the `900.000000` envelope, leaving `171.154444`
node-hours of margin.
Authenticated `R02/s05_rankio_t2_t2p5` job `4746663` is also formally and
independently inspected and recorded `accepted` through exact `t = 2.5` for
`1.062500` node-hours. Corrected E03 Stage I use is now `5.118334`
node-hours. Retained corrected R02 recost evidence JSON SHA-256 is
`39115585cf7ca866e3aade1e53c69aab6ee9217361d75db7946748905e1c0e5c`;
its authorization-facing corrected-matrix projection is `728.534445`
node-hours inside the `900.000000` envelope, leaving `171.465555`
node-hours of margin.
Authenticated `R02/s06_rankio_t2p5_t3` job `4746773` is also formally and
independently inspected and recorded `accepted` through exact `t = 3.0` for
`1.077500` node-hours. Corrected E03 Stage I use is now `6.195834`
node-hours. Retained corrected R02 recost evidence JSON SHA-256 is
`d2c5e27553426beb03c8c1882bd6473b682029c6ba67acd32bb06c5fb9a3ff01`;
its authorization-facing corrected-matrix projection is `732.734445`
node-hours inside the `900.000000` envelope, leaving `167.265555`
node-hours of margin.
Authenticated `R02/s07_rankio_t3_t3p5` job `4746953` is also formally and
independently inspected and recorded `accepted` through exact `t = 3.5` for
`1.151111` node-hours. Corrected E03 Stage I use is now `7.346945`
node-hours. Retained corrected R02 recost evidence JSON SHA-256 is
`468e787db48d98661cd3ffe125829a78f135bb7ca4c9104c766a00b693f586a3`;
its authorization-facing corrected-matrix projection is `753.345556`
node-hours inside the `900.000000` envelope, leaving `146.654444`
node-hours of margin. The clean `4144`-second interval exceeded the prior
`4000`-second operational guard. The reviewed next-segment threshold is
`4500` seconds for `R02/s08_rankio_t3p5_t4` only.
Authenticated `R02/s08_rankio_t3p5_t4` job `4747015` is also formally and
independently inspected and recorded `accepted` through exact `t = 4.0` for
`1.200000` node-hours. Corrected E03 Stage I use is now `8.546945`
node-hours. Retained corrected R02 recost evidence JSON SHA-256 is
`dd5b0f5b82cf81005c1a481ee83d1127aa43861dcdd8c22765170e28b0ad6c99`;
its authorization-facing corrected-matrix projection is `767.034445`
node-hours inside the `900.000000` envelope, leaving `132.965555`
node-hours of margin. The clean `4320`-second interval satisfies the reviewed
`s08` cap. A stricter independently reviewed next profile authorizes only
`R02/s09_rankio_t4_t4p5` with Slurm walltime `01:40:00`, Athena timeout
`01:30:00`, and a one-segment `4800`-second threshold. It does not ratchet.
Authenticated `R02/s09_rankio_t4_t4p5` job `4747087` is also formally and
independently inspected and recorded `accepted` through exact `t = 4.5` for
`1.208056` node-hours. Corrected E03 Stage I use is now `9.755001`
node-hours. Retained corrected R02 recost evidence JSON SHA-256 is
`8b10b1786db520740b687d989207f3cda6e0123f9bf2b62e75a5df3849511561`;
its authorization-facing corrected-matrix projection is `769.290001`
node-hours inside the `900.000000` envelope, leaving `130.710000`
node-hours of margin. The clean `4349`-second interval satisfies the reviewed
`s09` cap. Independent review renews the same profile only for
`R02/s10_rankio_t4p5_t5`: Slurm walltime `01:40:00`, Athena timeout
`01:30:00`, and a one-segment `4800`-second threshold retaining `600`
seconds on both timeout margins. It does not ratchet automatically.
Authenticated `R02/s10_rankio_t4p5_t5` job `4747146` is also formally and
independently inspected and recorded `accepted` through exact `t = 5.0` for
`1.226389` node-hours. Corrected E03 Stage I use is now `10.981390`
node-hours. Retained corrected R02 recost evidence JSON SHA-256 is
`5256fb7fd9b75f175a16c2c16ee8a995b8d6d9281a92b70a979ea62171835bf2`;
its authorization-facing corrected-matrix projection is `774.423334`
node-hours inside the `900.000000` envelope, leaving `125.576666`
node-hours of margin. The clean `4415`-second interval satisfies the reviewed
`s10` cap. Independent review renews the same profile only for
`R02/s11_rankio_t5_t5p5`: Slurm walltime `01:40:00`, Athena timeout
`01:30:00`, and a one-segment `4800`-second threshold retaining `600`
seconds on both timeout margins. It does not ratchet automatically.
After clean `s11` inspection, accounting, reconciliation, and recost, propose quarter units if elapsed time exceeds `4800` seconds, remaining cap headroom falls below `300` seconds, or the reviewed trend estimator projects the next half-unit at or above `4800` seconds. Any scientific, provenance, scheduler, storage, budget, or reconciliation failure blocks successor preparation.
Authenticated `R02/s11_rankio_t5_t5p5` job `4747202` is formally and
independently inspected and recorded `accepted` through exact `t = 5.5` for
`1.268889` node-hours. Corrected E03 Stage I use is now `12.250279`
node-hours. Retained corrected R02 recost evidence JSON SHA-256 is
`323d7cb3cbe108eee944045c189cff75ffc07d06833c49bcfa8315f07fad8d51`;
its authorization-facing corrected-matrix projection is `786.323334`
node-hours inside the `900.000000` envelope, leaving `113.676666`
node-hours of margin. The clean `4568`-second interval leaves `232` seconds
of headroom beneath the reviewed `s11` cap, below the reviewed `300`-second
floor. At that boundary, independent review historically authorized only
`R02/s12_rankio_t5p5_t5p75` on one node from authenticated `s11` terminal
siblings: Slurm walltime `01:05:00`, Athena timeout
`00:55:00`, and a one-segment `2700`-second threshold retaining `600`
seconds on both timeout margins. It does not ratchet automatically. The
reviewed timing estimator selects `4721` seconds for the next half-unit and
`2361` seconds for the proposed quarter-unit segment.
Authenticated `R02/s12_rankio_t5p5_t5p75` job `4747500` is formally and
independently inspected and recorded `accepted` through exact `t = 5.75` for
`0.660000` node-hours. Corrected E03 Stage I use is now `12.910279`
node-hours. Retained corrected R02 recost evidence JSON SHA-256 is
`d0a60e222138971e1bc9aae978bb3d62aee460f09000ee62f7f9ffdfe534165f`;
its authorization-facing corrected-matrix projection is `814.945556`
node-hours inside the `900.000000` envelope, leaving `85.054444`
node-hours of margin. The clean `2376`-second interval leaves `324` seconds
beneath the reviewed `s12` threshold. Independent review historically authorized
only `R02/s13_rankio_t5p75_t6` on one node from authenticated `s12` terminal
siblings: Slurm walltime `01:05:00`, Athena timeout `00:55:00`, and a
one-segment `2700`-second threshold retaining `600` seconds on both timeout
margins. It does not ratchet automatically. The reviewed timing estimator
selects `2468` seconds for the proposed `s13` quarter-unit segment.
Authenticated `R02/s13_rankio_t5p75_t6` job `4747834` is formally and
independently inspected and recorded `accepted` through exact `t = 6.0` for
`0.599722` node-hours. Corrected E03 Stage I use is now `13.510001`
node-hours. Retained corrected R02 recost evidence JSON SHA-256 is
`4bc19f5b44587e73869982425b1cf0ad459547c3e0c5de1ab1b2b5cb366c6322`;
its authorization-facing corrected-matrix projection is `814.945556`
node-hours inside the `900.000000` envelope, leaving `85.054444`
node-hours of margin. The clean `2159`-second interval leaves `541` seconds
beneath the reviewed `s13` threshold. Independent review historically authorized
only `R02/s14_rankio_t6_t6p25` on one node from authenticated `s13` terminal
siblings: Slurm walltime `01:05:00`, Athena timeout `00:55:00`, and a
one-segment `2700`-second threshold retaining `600` seconds on both timeout
margins. It does not ratchet automatically. The reviewed timing estimator
selects `2468` seconds for the proposed `s14` quarter-unit segment.
Authenticated `R02/s14_rankio_t6_t6p25` job `4748138` is formally and
independently inspected and recorded `accepted` through exact `t = 6.25` for
`0.655278` node-hours. Corrected E03 Stage I use is now `14.165279`
node-hours. Retained corrected R02 recost evidence JSON SHA-256 is
`5ba7aec52d3e7824cf4ccf52f95ce1e46368c84ed59f0525e2bbbcfc167fde91`;
its authorization-facing corrected-matrix projection is `829.101112`
node-hours inside the `900.000000` envelope, leaving `70.898888`
node-hours of margin. The clean `2359`-second interval leaves `341` seconds
beneath the reviewed `s14` threshold. Independent review historically authorized
only `R02/s15_rankio_t6p25_t6p5` on one node from authenticated `s14` terminal
siblings: Slurm walltime `01:05:00`, Athena timeout `00:55:00`, and a
one-segment `2700`-second threshold retaining `600` seconds on both timeout
margins. It did not ratchet automatically. The reviewed timing estimator
selected `2559` seconds for the proposed `s15` quarter-unit segment.
Authenticated `R02/s15_rankio_t6p25_t6p5` job `4752118` is formally and
independently inspected and recorded `accepted` through exact `t = 6.5` for
`0.653611` node-hours. Corrected E03 Stage I use was then `14.818890`
node-hours. Retained corrected R02 recost evidence JSON SHA-256 is
`3603077400bd4530001155b83bfc4ac1053ec18505c56a6c29ffadc795bd1255`;
its authorization-facing corrected-matrix projection is `829.101112`
node-hours inside the `900.000000` envelope, leaving `70.898888`
node-hours of margin. The clean `2353`-second interval leaves `347` seconds
beneath the reviewed `s15` threshold. Independent review historically authorized
only `R02/s16_rankio_t6p5_t6p75` on one node from authenticated `s15` terminal
siblings: Slurm walltime `01:05:00`, Athena timeout `00:55:00`, and a
one-segment `2700`-second threshold retaining `600` seconds on both timeout
margins. It did not ratchet automatically. The reviewed timing estimator
selected `2559` seconds for the proposed `s16` quarter-unit segment.
Authenticated `R02/s16_rankio_t6p5_t6p75` job `4753294` is formally and
independently inspected and recorded `accepted` through exact `t = 6.75` for
`0.647222` node-hours. Corrected E03 Stage I use is now `15.466112`
node-hours. Retained corrected R02 recost evidence JSON SHA-256 is
`5b997623a3f8d83c200034f42e0c9f1b9a09c811bcb46e2fc1a0f7c1eb58815e`;
its authorization-facing corrected-matrix projection is `829.101112`
node-hours inside the `900.000000` envelope, leaving `70.898888`
node-hours of margin. The clean `2330`-second interval leaves `370` seconds
beneath the reviewed `s16` threshold. Independent review historically authorized
only `R02/s17_rankio_t6p75_t7` on one node from authenticated `s16` terminal
siblings after the post-F-099 documentation archive, reconciliation, and
queue/shared-root audit gates closed. It used Slurm walltime `01:05:00`, Athena
timeout `00:55:00`, and a one-segment `2700`-second threshold retaining `600`
seconds on both timeout margins. It did not ratchet automatically. The reviewed
timing estimator selected `2559` seconds for that `s17` quarter-unit segment.
Authenticated `R02/s17_rankio_t6p75_t7` job `4754008` is formally and
independently inspected and recorded `accepted` through exact `t = 7.0` for
`0.667222` displayed node-hours. Corrected E03 Stage I
use is now `16.133334` displayed node-hours
(`58080 / 3600 = 16.133333333333333` exact
node-hours). Retained corrected R02 F-100 recost evidence JSON SHA-256 is
`63438d8843b86815a2e25d4d6de6ef78ce9756c561fe9db3e762965ada73202f`; its authorization-facing corrected-matrix projection is
`829.101112` node-hours inside the `900.000000`
envelope, leaving `70.898888` node-hours of
margin. The clean `2402`-second interval leaves
`298` seconds beneath the reviewed `s17`
threshold. Independent review narrowly and conditionally authorizes only
`R02/s18_rankio_t7_t7p25` on one node from authenticated `s17` terminal siblings
after these documentation bytes are committed, the
resulting actual post-F-100 source bundle is
archived and cataloged with SHA-256, hardened
reconciliation is rerun, and the queue/shared-root audit repeats the explicit
stale beta-25 acknowledgement. Use Slurm walltime `01:05:00`, Athena timeout
`00:55:00`, and a one-segment `2700`-second threshold retaining `600` seconds on
both timeout margins. It does not ratchet automatically. The reviewed timing
estimator selects `2559` seconds for the proposed
`s18` quarter-unit segment.
Reset the planning ceiling for that new epoch to an incremental `4000`
node-hours while continuing to report the historical `0.851670` debug and
`9.962778` Stage I node-hours. A reset is an accounting boundary for future
authorization; it does not erase consumed allocation.
The corrected worktree now selects explicit restart-retained
`mks24_random_unprojected` and `mks24_alfvenic_perpendicular` forcing
policies, sets `spectrum = power_law` in every paper deck, records an
authoritative `time/restart_time` marker in new restart dumps, and isolates
future work beneath `E03-forcing-policy`. The E03 helper was configured to
fail closed until a reviewed qualification token bound the immutable corrected
executable and revision. No E03 preparation was authorized until that
executable completed Frontier qualification and the token was retained. The
schema-v2 panel-status
inventory separately retains eleven admitted comparison families, eleven
reference-blocked rows, and the external Figure 10 model boundary; admitted
comparisons remain `not_run` until reviewed numeric criteria are recorded.
The independent E03 helper audit is closed after targeted replay and recovery
probes. Corrected immutable Frontier qualification is also closed for revision
`9e07542281e4e6d125582f253df3ad2e3b8b154d`, executable SHA-256
`68f243f9204df388b24365ae65a567f6f567dbe422a6d7a43b9fb4a499ef118c`,
and retained source-bundle SHA-256
`c39d55809989d20aa5438711803f4fd43237fa9284c7e15183e84fdd693d3687`.
Jobs `g024` and `g025` qualify the explicit planar-Alfvenic and unprojected
random policies. `g026` supplies an uninterrupted eight-rank reference;
`g027` resumes authenticated restart siblings at
`time/restart_time = 0.00968183`, crosses an OU refresh, and matches the
reference within retained-format tolerances; `g028` closes one-rank/eight-rank
decomposition identity. After two pre-submit wording corrections, `g029c`
qualifies passive-Delta semantics. `g030b` reaches exact `t = 2.0` with zero
strict counters, terminal `lf_hwproj = 10084905222`, forcing-work relative
residual `1.8681072370071585e-12`, nine shared snapshots, five shared
checkpoints, and fourteen MPI-I/O records. `g031` launches the standard
`192 x 192 x 384` layout with 27 meshblocks per GPU rank, reaches exact
`t = 0.01`, retains complete initial/final eight-sibling ranked products, and
closes startup, one-node memory-fit, and retained-output sizing only.
Canonical `g031` evidence JSON SHA-256 is
`a16415c5f9c557a0dad35936c0bb0e89658da9e3b78a4d86852ead4e006f0d98`;
the supplemental independently generated record remains retained with
SHA-256
`e000a789d10e749b476d09a1be5e6526f0f0dc2dbe0e46c47b0ef09bc0f56ff6`.
Reviewed approval token SHA-256
`e3fec9f35da42121b902f41ef752021f375aab38bc34c5ff75ae8780bfbca635`
was retained atomically at `2026-05-30T21:33:31+00:00`; it binds the corrected
revision and executable. Reconciliation then passed, and fresh E03
`R02/s00_rankio_t0_t0p1` job `4745922` was submitted from `t = 0`, formally
inspected, and recorded `accepted` at exact `t = 0.1` for `0.188889`
node-hours. Its record-time `segment_inspection.json` and embedded
`scientific_inspection` payload SHA-256 are
`f4b0a31304e31610ec9ec83133f79b747fb0d8cf0f9b3b013febc1cddc2ae784`;
the forcing-work relative residual is `1.920820308941694e-11`, strict LF
safety counters remain zero, and complete initial/final eight-sibling ranked
snapshot and restart groups are retained.
Authenticated continuation `R02/s01_rankio_t0p1_t0p25` job `4746154`
reaches exact `t = 0.25` in `0.289167` node-hours. Its record-time
`segment_inspection.json` SHA-256 is
`f0c4229d9a2cd86bed76e87efd32204cc504fbd2d8585d0f5d4468bd25328d6f`;
strict LF failure counters remain zero, terminal `lf_hwproj = 28701760`,
complete eight-rank terminal snapshot/restart groups are retained, and the
sampled-history forcing-work relative residual is
`1.7741014899016423e-11`.
Authenticated continuation `R02/s02_rankio_t0p25_t1` job `4746182` reaches
exact `t = 1.0` in `1.465556` node-hours. Its record-time
`segment_inspection.json` SHA-256 is
`482c7f0f9552e319cc50e92f6d84f17565a90f9562b06f1d67ef1c57d240a2dd`;
strict LF failure counters remain zero, terminal `lf_hwproj = 65252911379`,
complete terminal eight-rank restart siblings are retained, and the
sampled-history forcing-work relative residual is
`8.616579960442532e-12`.
Authenticated continuation `R02/s03_rankio_t1_t1p5` job `4746356` reaches
exact `t = 1.5` in `1.048611` node-hours. Its record-time
`segment_inspection.json` SHA-256 is
`5a66f4e037fe08f176b51c7f949db0377bc9270e666ae1e1c69dda85d7f96633`;
strict LF failure counters remain zero, terminal `lf_hwproj = 106136472818`,
complete two-group eight-rank snapshots and terminal eight-rank restart
siblings are retained, and the sampled-history forcing-work relative
residual is `4.978727845741857e-13`.
Authenticated continuation `R02/s04_rankio_t1p5_t2` job `4746435` reaches
exact `t = 2.0` in `1.063611` node-hours. Its record-time
`segment_inspection.json` SHA-256 is
`fee92e6d8b484d1050266f6d5600577d19099df02335299114da59dd42050606`;
strict LF failure counters remain zero, terminal `lf_hwproj = 123445535090`,
complete two-group eight-rank snapshots and terminal eight-rank restart
siblings are retained, and the sampled-history forcing-work relative
residual is `3.034159691118515e-12`.
Authenticated continuation `R02/s05_rankio_t2_t2p5` job `4746663` reaches
exact `t = 2.5` in `1.062500` node-hours. Its record-time
`segment_inspection.json` SHA-256 is
`57fb0b9f1a19c0125aec0619802f2c41bd91ac97c51f32e70ceae9a434355da4`;
strict LF failure counters remain zero, terminal `lf_hwproj = 171943591926`,
complete two-group eight-rank snapshots and terminal eight-rank restart
siblings are retained, and the sampled-history forcing-work relative
residual is `3.0678457367645077e-12`.
Authenticated continuation `R02/s06_rankio_t2p5_t3` job `4746773` reaches
exact `t = 3.0` in `1.077500` node-hours. Its record-time
`segment_inspection.json` SHA-256 is
`4b32798174222c9b9ca165285397ed9a02a411bb3e3dbf09e1f10611eebc70db`;
strict LF failure counters remain zero, terminal `lf_hwproj = 184131904300`,
complete two-group eight-rank snapshots and terminal eight-rank restart
siblings are retained, and the sampled-history forcing-work relative
residual is `2.8133353495278692e-12`.
Authenticated continuation `R02/s07_rankio_t3_t3p5` job `4746953` reaches
exact `t = 3.5` in `1.151111` node-hours. Its record-time
`segment_inspection.json` SHA-256 is
`4714af847c1ce10517258ea34f6151eb2a60d7fbd1bf1717d2de38d2e6822933`;
strict LF failure counters remain zero, terminal `lf_hwproj = 192268745176`,
complete two-group eight-rank snapshots and terminal eight-rank restart
siblings are retained, and the sampled-history forcing-work relative
residual is `2.6602098521945525e-13`.
Authenticated continuation `R02/s08_rankio_t3p5_t4` job `4747015` reaches
exact `t = 4.0` in `1.200000` node-hours. Its record-time
`segment_inspection.json` SHA-256 is
`be705544f6c8f4a71f5501130b1e317555cdf1cac0663fba267acec82d967649`;
strict LF failure counters remain zero, terminal `lf_hwproj = 203117871743`,
complete two-group eight-rank snapshots and terminal eight-rank restart
siblings are retained, and the sampled-history forcing-work relative
residual is `1.3875872397103518e-13`.
Authenticated continuation `R02/s09_rankio_t4_t4p5` job `4747087` reaches
exact `t = 4.5` in `1.208056` node-hours. Its record-time
`segment_inspection.json` SHA-256 is
`0373fa3505beaadbd818db5271dae13bf69a41b1683bb561d558d8d71eff7581`;
strict LF failure counters remain zero, terminal `lf_hwproj = 207989256307`,
complete two-group eight-rank snapshots and terminal eight-rank restart
siblings are retained, and the sampled-history forcing-work relative
residual is `1.9837412357887464e-12`.
Authenticated continuation `R02/s10_rankio_t4p5_t5` job `4747146` reaches
exact `t = 5.0` in `1.226389` node-hours. Its record-time
`segment_inspection.json` SHA-256 is
`8fb8069dc1f75448a8f79a38de98141c5d6dfed83bc3530eb633f8c009655014`;
strict LF failure counters remain zero, terminal `lf_hwproj = 220146873111`,
complete two-group eight-rank snapshots and terminal eight-rank restart
siblings are retained, and the sampled-history forcing-work relative
residual is `5.300795515586916e-12`.
Authenticated continuation `R02/s11_rankio_t5_t5p5` job `4747202` reaches
exact `t = 5.5` in `1.268889` node-hours. Its record-time
`segment_inspection.json` SHA-256 is
`928c281788b8ba359eb499aa607529b6c8eba2f8236e374f808d57fe239d9c1e`;
strict LF failure counters remain zero, terminal `lf_hwproj = 224321212896`,
complete two-group eight-rank snapshots and terminal eight-rank restart
siblings are retained, and the sampled-history forcing-work relative
residual is `5.093197346120908e-12`.
Commits `d210cdd5`, `eab6e12b`, and `7fae0bcf` retain bundle-backed source
provenance, complete rank-local debug-restart archival, an all-user-job debug
preflight, reviewed top-level shared-root acknowledgement, and prepare-time
rejection of absent override targets. The bundle used for `g019`,
`source-archives/athenak-feature-cgl-through-ec98c25e.bundle`, has SHA-256
`938170bd0e4757befe614ee9c471875d539de92fde93670d32a96e6d1771cd58`.
Replacement-driver debug job `4743463` (`g014`) then exposed and retained a
deck/override preflight gap before physics startup, consuming `0.002778`
node-hours. Corrected startup job `4743465` (`g015`) used explicit hard-wall
and rank-local settings in retained reduced input SHA-256
`f41e8a3fc98bc2d9881ad9281fa9af76c2ecb2043bec95a19f2dfda015d39ec5`,
reached `t = 0.01` with zero strict LF safety counters, and consumed
`0.005556` node-hours. Its retained evidence JSON has SHA-256
`e2f4a047f601e6d1539d77a85a09d65e4df8b338c22eb6638441cdd6af75cd70`.
Follow-up `g016` (job `4743470`) loaded the authenticated terminal eight-rank
restart set and reached `t = 0.02` cleanly. The first uninterrupted comparator,
`g017` (job `4743472`), was retained as `inconclusive` because stopping `g015`
exactly at `t = 0.01` clipped its final timestep and changed the subsequent
integration sequence. Natural-cycle continuation `g018` (job `4743473`) then
resumed the complete `g017` `.00013.rst` set and matched the uninterrupted
`g017` terminal histories, rank-local snapshot fields, and turbulence-force
slice within declared retained-format tolerances. Its evidence JSON has
SHA-256
`4d898cf4aa4427b81164e4eb688f011adbbb277874b1ca7453a7a2c1e60b580f`.
This closes reduced modal restart identity. One-rank comparator `g019` (job
`4743474`) then matched the uninterrupted eight-rank `g017` trajectory within
declared retained-format tolerances, closing reduced GPU/MPI decomposition.
Its evidence JSON has SHA-256
`509ec34d4f1f27dfdbb9b507d2ca663d7233d168a0b83775b0ab532dc283c752`.
Reduced passive-Alfvénic `g020` (job `4743480`) and active-random `g021` (job
`4743483`) smoke then passed with zero strict counters; their evidence JSON
SHA-256 values are
`4abeddeb25db489b3334be240b837e703d035585d0814fc79d004e273cc34b3f`
and `aaee51f186ac0f223fa4a1c35b238bb1c2e22f2a1ac9075ad1e28ffea53e8498`.
Reduced nonlinear hard-wall `g022` (job `4743484`) then reached `t = 2.0`
with zero strict counters, active algebraic hard-wall projection
(`lf_hwproj = 10055979902`), forcing-work relative residual
`1.863e-12`, nine readable shared snapshots, five retained shared
checkpoints, and fourteen retained MPI-I/O net-write records. It used
`0.431667` node-hours. Its evidence JSON has SHA-256
`339a3309276a93be7bb55c196f367ead9754f7c2ebdf94afdb9aeedd37f026cf`.
This closes reduced replacement-driver nonlinear correctness only. The
standard-layout rank-local startup and retained-output sizing gate was then
closed by `g023` (job `4743662`). That one-node/eight-rank run launched the
guarded `192 x 192 x 384` layout, advanced through `t = 0.01` in `79`
allocated seconds with zero strict counters, forcing-work relative residual
`1.053e-11`, two complete eight-rank binary groups, and two complete
eight-rank modal restart groups. Each terminal snapshot shard is `63710281`
bytes (`509682248` bytes aggregated), and each modal restart shard is
`172670767` bytes (`1381366136` bytes aggregated). The assembled endpoint
reader passes at shape `384 x 192 x 192`. `g023` used `0.021944` node-hours;
its evidence JSON has SHA-256
`effba10b616246e1c9f6863c6a068f674551fa9ac4b4a02bdb92ead0f963c14f`.
This closes replacement-driver Frontier debug startup and sizing
qualification only, not production runtime or long-time statistics.

The focused post-merge local gate is also retained. The serial CGL suites
passed `39` tests in `23.89` seconds; evidence JSON SHA-256 is
`3fc869c0d7455f1221d7838de900c0db7426673b2730cd954f264b8934fc0c29`.
Scheduled one-node MPI CPU job `4743666` passed both one-rank/four-rank
decomposition regressions in `2.81` pytest seconds and consumed `0.002222`
node-hours separately from the debug-helper ledger; evidence JSON SHA-256 is
`397aced4e90aa7e77ab17da05a522ad4c99c9d0e15bd48d86ce20aa195fb36d3`.
Utility commits `476f9dbd`, `9c5c1b40`, and `bad8ba05` select only the
explicit accepted restart lineage during bundling and admit normal
sampled-history separation across authenticated restart boundaries while
rejecting gaps inconsistent with the recorded cadence. They assembled the
ranked prefix through `t = 2.0` into diagnostic bundle
`runs/bundles/diagnostic-prefixes/R16_rankio_t2_prefix_20260525`, excluding
the abandoned shared-output branch. The default production-window
`paper-analyze` action selected no prefix snapshots because its declared
window is `t = 8`--`10` and wrote one history-only diagnostic figure; a
separately labeled transient-window (`t = 0`--`2`) analyzer check read all ten
retained snapshots, selected 101 history rows, passed its synthetic check,
and wrote nine diagnostic figures. Neither action is reproduction evidence.
The first scientific objective is an independent AthenaK reproduction
of the attainable numerical results in Majeski, Kunz, and Squire (2024,
hereafter MKS24), before any new weak-guide or background-collisionality
experiment is executed.

Primary manuscript source: `docs/cgl_lf_validation.tex`

Governing prose guide: `docs/writing_style_guide.md`

Detailed implementation and evidence record:
`docs/cgl_lf_mks24_reproduction_implementation_plan.md`

Pinned reference source: official arXiv `2405.02418v2`, staged by
`scripts/stage_cgl_lf_mks24_reference.py` with checksums retained outside the
tracked source tree.

Completion boundary: this document defines what must be done. It is not
evidence that MKS24 has been reproduced, that a new scientific result exists,
or that a manuscript is complete.

## 1. Revision Decision

The previous version of this plan did **not** include everything needed to
reproduce MKS24. It defined a new cubic-box experiment with weak and stronger
guide-field anchors, a scan in background collision frequency `nu_coll`, and
new target values for magnetic fluctuation amplitude and sonic Mach number.
Those are extension experiments, not the published MKS24 experiment.

MKS24 instead studies collisionless high-beta CGL-Landau-fluid (CGL-LF)
turbulence in an elongated mean-field box. Its collisionality study is a
variation of the anomalous microinstability scattering limiter
`nu_lim`, applied only after mirror or firehose thresholds are crossed, not a
general uniform background-collision scan. The manuscript must reproduce this
distinction exactly before using the extension to ask new questions.

The project is reorganized into two non-interchangeable stages:

| Stage | Objective | Permission boundary |
| --- | --- | --- |
| I. MKS24 reproduction | Reproduce the published CGL-LF simulation setup, numerical comparisons, and all quantitatively attainable figure claims; account honestly for the external hybrid-kinetic comparator. | This is the only production-science campaign authorized for planning under the current `4000` node-hour ceiling. Execution still requires the gates below. |
| II. New extension | Study weak versus stronger guide fields and background collisionality after Stage I establishes a credible baseline. | Retained as a later candidate protocol; it is not presently authorized and must be recosted and approved after Stage I. |

The manuscript claim discipline is:

| Claim level | Permitted statement | Required evidence |
| --- | --- | --- |
| Implementation | AthenaK implements the MKS24 CGL-LF case definitions and diagnostics used here. | Input/workflow tests, physical-equation audit, forcing and limiter audit, restart/GPU qualification, archived provenance. |
| CGL-LF numerical reproduction | AthenaK independently reproduces a specified MKS24 CGL-LF result within stated comparison uncertainty. | Matched production cases, explicit normalizations, reference data or admissible digitization, convergence/statistical checks, and a per-panel pass/fail table. |
| Full-paper reproduction | Every claimed MKS24 result, including non-AthenaK external comparisons, has been reproduced or independently re-evaluated. | In addition to CGL-LF reproduction, a qualified path for MKS24 Figure 10 and any other externally generated result. |
| Extension inference | A new guide-field or background-collisionality result extends the reproduced baseline. | Stage I acceptance followed by separately approved, funded, executed, and analyzed Stage II cases. |
| Not permitted | The method is perfect or reproduces results that have not passed the stated gates. | No numerical campaign establishes perfection; absent data are not evidence. |

## 2. MKS24 Source Contract

The source to reproduce is:

```text
S. Majeski, M. W. Kunz, and J. Squire (2024)
Self-organization in collisionless, high-beta turbulence
Journal of Plasma Physics 90, 535900601
arXiv:2405.02418v2
```

Before modifying reproduction inputs or interpreting a panel, future agents
must restage or verify the pinned source and reread the governing passages:

| Source portion | Reproduction relevance |
| --- | --- |
| `MKS24.tex:204-239` | CGL-MHD equations, 3+1 Landau-fluid heat fluxes, cold-electron and microinstability-model assumptions. |
| `MKS24.tex:463-469` | Domain, resolution, beta, standard closure parameters, forcing, and duration. |
| `MKS24.tex:478-496` | Fourier spectra, pressure-stress transfer, and scale-dependent alignment diagnostics. |
| `MKS24.tex:503-655` | Main numerical results through heat-flux sensitivity. |
| `MKS24.tex:662-698` | Microinstability limiter-rate scan and alternative firehose-threshold discussion. |
| `MKS24.tex:760-918` | Appendices relating the reduced theory to RKMHD, Braginskii stress, and larger density fluctuations. |

The implementation record is already substantially advanced: it defines
guarded production inputs, diagnostic products, reference extractors, known
normalization gaps, reduced tests, and GPU qualification evidence. This
manuscript plan does not supersede those technical records. It determines
which of them must be completed and accepted before writing scientific claims.

### 2.1 What "completely reproduce" can mean

MKS24 contains theory, CGL-LF simulations, and a hybrid-kinetic comparison
whose data originated outside the CGL-LF solver. An AthenaK CGL-LF campaign
can independently reproduce:

1. The stated CGL-LF closure, numerical regime, active/passive construction,
   forcing family, limiter-rate variation, and CGL-LF diagnostics.
2. Every CGL-LF simulation panel for which matching observables and a
   qualified normalization/reference comparison are available or obtained.
3. The paper's analytic scalings and their connection to the reproduced
   numerical observables, by a documented derivation audit.

It cannot silently claim independent reproduction of MKS24 Figure 10, which
uses a `Pegasus++` hybrid-kinetic simulation reported from prior work. A
full-paper reproduction claim therefore requires one of the following:

1. Obtain the underlying hybrid-kinetic dataset and reproduce the plotted
   alignment analysis with provenance.
2. Obtain and execute a scientifically equivalent hybrid-kinetic workflow,
   including its model assumptions and convergence evidence.
3. Restrict the claim explicitly to complete reproduction of MKS24's CGL-LF
   simulation results, treating Figure 10 as external comparison context.

Until that decision is resolved, the manuscript must use the third wording
and must not say that every result in MKS24 was independently reproduced.

## 3. Stage I Physical and Numerical Specification

### 3.1 CGL-LF model

Use the CGL pressure equations with the MKS24 3+1 Landau-fluid closure:

```math
q_{\perp}=-{v_{\rm th,\parallel}\over\sqrt{\pi}|k_\parallel|}
\left[\rho\nabla_\parallel\left({p_\perp\over\rho}\right)
-p_\perp\left(1-{p_\perp\over p_\parallel}\right)
{\nabla_\parallel B\over B}\right],
```

```math
q_{\parallel}=-{2v_{\rm th,\parallel}\over\sqrt{\pi}|k_\parallel|}
\rho\nabla_\parallel\left({p_\parallel\over\rho}\right),\qquad
v_{\rm th,\parallel}=\sqrt{2p_\parallel/\rho}.
```

Record AthenaK's unit conversion and coefficient mapping in every production
bundle. Preserve the paper assumptions of initially isotropic ion pressure,
cold electrons for the directly reproduced CGL-LF calculation, and
unresolved ion-Larmor-scale instabilities represented through the declared
limiter closure.

### 3.2 Published setup to match

| Quantity | MKS24 reproduction value |
| --- | --- |
| Domain | Fully periodic `[L_x, L_y, L_z] = [1, 1, 2]`, with `L_z = L_parallel` and `L_x = L_y = L_perp`. |
| Guide field | Uniform `B_0` along `z`; AthenaK decks state the selected unit mapping explicitly. |
| Standard mesh | `n_perp = 192`, `n_parallel = 384`. |
| Scale-separation meshes | Active Alfvenic beta-10 cases at `n_perp = 96, 192, 384`, with `n_parallel = 2 n_perp`. |
| Initial beta values | `beta0 = 1, 10, 100`. |
| Standard LF scale | `|k_parallel| = 4*pi/L_parallel = 2*pi` in the stated box. |
| Standard limiter | Hard-wall equivalent of `nu_lim = 10^10 v_A/L_perp`. |
| Instability thresholds | Mirror `beta Delta > 1`; parallel/fluid firehose `beta Delta < -2`. |
| Background collisions | `nu_coll = 0.0` for direct MKS24 CGL-LF reproduction. |
| Duration | At least `t_f = 10 L_perp/v_A`; statistical windows beyond `t v_A/L_perp = 6` and no shorter than `2 L_perp/v_A`. |

The paper reports a beta-100 pressure normalization in its Figure 13
discussion that does not by itself fix all AthenaK dimensional spectral
ordinate transforms. The input normalization must remain explicit, and
dimensional panel comparisons must be blocked until their mapping is
qualified rather than normalized by guesswork.

### 3.3 Published forcing to match

| Quantity | MKS24 reproduction value |
| --- | --- |
| Process | Ornstein-Uhlenbeck correlated velocity forcing. |
| Fiducial injection | `d_t E_K = 0.32 rho0 v_A^2 L_perp^3`. |
| Fiducial correlation time | `t_corr = L_parallel/v_A`. |
| Forced shell | `k in (2*pi/L_parallel) * [1, 3]`. |
| Power distribution | Proportional to `k^-2`. |
| Alfvenic forcing | Perturb only velocity perpendicular to `z`, incompressibly at the outer scale. |
| Random forcing | Velocity forcing without the Alfvenic directional constraint, while retaining declared time correlation. |
| Sonic-correlation case | The beta-100 random-forcing compressive test uses `t_corr` corresponding to `L_parallel/v_th`. |

If the MKS24 source does not publish a forcing random seed or an output needed
for pointwise identity, record that as unavailable source information. The
reproduction target then becomes statistical agreement of declared
observables, not equality of turbulent fields.

### 3.4 Collisionality terminology

For Stage I, distinguish two different inputs:

| Parameter | Meaning | Stage I role |
| --- | --- | --- |
| `nu_coll` | Uniform background collisional relaxation in AthenaK. | Fixed to `0.0`; a nonzero scan is not an MKS24 reproduction case. |
| `nu_lim` / AthenaK limiter-scattering representation | Anomalous scattering activated only in cells beyond mirror/firehose thresholds and driving pressure back toward marginality. | Reproduce MKS24's limiter sensitivity at beta 100 using `20`, `200`, and hard-wall-equivalent `10^10 v_A/L_perp`. |

The later uniform-`nu_coll` experiment remains useful, but it answers a new
question and must not be described as the reproduction of MKS24's
microinstability-limiter scan.

## 4. Stage I Run Matrix

The repository already defines guarded inputs for the required CGL-LF case
families. Reuse one accepted run wherever a plotted role is physically and
numerically identical; do not execute aliases twice merely because a figure
uses them twice.

### 4.1 Sixteen mapped CGL-LF simulations and one unresolved definition

| ID | Existing input role | Configuration | MKS24 purpose |
| --- | --- | --- | --- |
| R02 | `standard_active_alfvenic_beta10` | Active, Alfvenic, beta `10`, standard mesh/closure | Central backbone: Figures 1, 2, 4-9, 11, and 12 nominal role. |
| R03 | `standard_active_alfvenic_beta100` | Active, Alfvenic, beta `100`, hard wall | Figures 3-4 beta results and Figure 13 hard-wall role. |
| R04 | `standard_active_random_beta10` | Active, random, beta `10`, standard closure | Forcing comparison. |
| R05 | `standard_active_random_beta100` | Active, random, beta `100`, Alfvenic correlation time | Forcing/density comparison. |
| R06 | `standard_passive_alfvenic_beta10` | Passive, Alfvenic, beta `10` | Active/passive comparison. |
| R07 | `standard_passive_alfvenic_beta100` | Passive, Alfvenic, beta `100` | Unstable-volume/beta comparison and Figure 13 passive context. |
| R08 | `standard_passive_random_beta10` | Passive, random, beta `10` | Forcing/passive comparison. |
| R09 | `standard_passive_random_beta100` | Passive, random, beta `100` | Forcing/passive and Figure 3 beta comparison. |
| R10 | `compressive_active_random_beta1` | Active, random, beta `1` | Figure 3 beta-dependent compressive spectrum. |
| R11 | `compressive_active_random_beta100_sonic` | Active, random, beta `100`, sonic `t_corr` | Figure 3 sonic-correlation experiment. |
| R12 | `heat_flux_beta10_strong` | Active, Alfvenic, beta `10`, `|k_parallel|/100` | Figure 12 stronger-heat-flux limit. |
| R13 | `heat_flux_beta10_weak` | Active, Alfvenic, beta `10`, `100 |k_parallel|` | Figure 12 weak/double-adiabatic-like limit. |
| R14 | `nulim_beta100_20` | Active, Alfvenic, beta `100`, `nu_lim = 20` | Figure 13 limiter sensitivity. |
| R15 | `nulim_beta100_200` | Active, Alfvenic, beta `100`, `nu_lim = 200` | Figure 13 limiter sensitivity. |
| R16 | `scale_separation_beta10_nperp96` | Active, Alfvenic, beta `10`, `96 x 96 x 192` | Figure 11 low-resolution case. |
| R17 | `scale_separation_beta10_nperp384` | Active, Alfvenic, beta `10`, `384 x 384 x 768` | Figure 11 high-resolution case. |

The repository also defines `standard_active_alfvenic_beta1`. The current
checksum-qualified extractor mappings do not map that case to a displayed
MKS24 result: the source explicitly identifies Figure 3's beta-1 curve as
randomly driven (`R10`), while the admitted Figure 4 and Figure 9 mappings use
`R02`-`R09`. Treat the separate Alfvenic beta-1 deck as unresolved inventory,
not as an authorized production case, unless the Phase A source audit
identifies and records a specific published claim requiring it.

`R03` supplies the otherwise duplicated beta-100 hard-wall `nu_lim` role, and
`R02` supplies the `n_perp = 192` scale-separation and nominal heat-flux roles.
The currently mapped run matrix is therefore sixteen unique calculations
after role and input-alias reuse.

### 4.2 Alias and execution contract

The workflow exposes some published roles through more than one campaign
entry. The production manifest must resolve these aliases before reserving
compute:

| Physical case | Existing workflow role(s) | Required handling |
| --- | --- | --- |
| `R02` | `paper_standard_active_alfvenic_beta10`; reused by the nominal Figure 12 comparison and the `n_perp = 192` Figure 11 comparison. | Execute once and link all three panel roles to the same accepted run bundle. |
| `R03` | `paper_standard_active_alfvenic_beta100` and `paper_nulim_beta100_hardwall`. | Compare parsed physics, numerical, forcing, output, and analysis parameters after ignoring the basename/comments; if they remain identical, execute once and link the hard-wall role to that output. |

If an intended alias has a material parameter difference, it is not an alias:
assign a new `Rxx` identifier, revise the run count and budget, and document
why the additional calculation is necessary before submission. In particular,
the five paper workflow commands must not be executed independently without
deduplicating their overlapping case roles.

### 4.3 Reproduction is not achieved by input files alone

The existing input definitions are necessary but not sufficient. Before a
case contributes to a claim, it must have:

1. A retained submitted input, executable hash, code revision, module/build
   provenance, forcing metadata, and case identity.
2. Strict numerical safety and restart evidence appropriate to its limiter
   model and resolution.
3. A validated late-time analysis interval satisfying the published duration
   and averaging requirements.
4. Retained source snapshots/histories, derived products, analysis
   configuration, checksums, and cost/storage accounting.
5. A comparison outcome for every MKS24 panel to which it contributes.

## 5. Figure-by-Figure Reproduction Gate

The paper shall not state "MKS24 is reproduced" until the following table has
been completed with pass/fail/blocked outcomes and justified uncertainties.

| MKS24 figure or result | Required cases/products | Reproduction gate |
| --- | --- | --- |
| Figure 1: beta-10 active/passive `beta Delta` slices | `R02`, `R06`; common-time or statistically declared slice visual; instability mask. | Recreate the qualitative active/passive spatial contrast with source-provenance caption; do not treat a hand-selected attractive snapshot as quantitative proof. |
| Figure 2(a): pressure-density joint PDFs | `R02`; joint `delta p_parallel`, `delta p_perp`, density surface. | Use the existing source-qualified sampled-surface transform or a superior author dataset; report residuals and uncertainty. |
| Figure 2(b): unstable-volume histories | Active/passive Alfvenic/random beta `10`, `100` roles (`R02`-`R09`). | Compare all eight curves against the retained checksum-qualified reference histories and the reported active/passive volume-fraction conclusion. |
| Figure 3: compressive velocity spectra | Panel (a): `R03`, `R05`, `R11`; panel (b): `R10`, `R04`, `R11`, `R09`. | Qualify absolute spectral normalization or obtain reference data; then test forcing-correlation and beta-dependence conclusions. |
| Figure 4: density PDF and density spectra | The eight Alfvenically correlated active/passive beta-10/beta-100 roles, `R02`-`R09`. | Compare admitted normalized density PDFs; unblock spectral ordinate mapping before claiming reproduction of panel (b). |
| Figure 5: energy spectra and local-field eddy anisotropy | `R02`, `R04`; local-field structure/eddy products. | Compare admitted normalized eddy-scale curves; qualify panel (a) spectral ordinate mapping. |
| Figure 6: pressure/magnetic and strain spectra | `R02`, `R04`. | Establish comparison normalization for dimensional spectra and test pressure-balance/parallel-strain-suppression claims. |
| Figure 7: pressure-anisotropy gradient and transfer | `R02`, `R04`, `R06`; normalized transfer product. | Compare admitted dimensionless transfer curves; unblock upper-panel dimensional spectra before a complete panel claim. |
| Figure 8: scale-dependent alignment PDFs | `R02`, `R06`. | Compare existing calibrated selected-shell alignment PDFs and preserve color-map uncertainty. |
| Figure 9: alignment peaks across cases | The eight standard active/passive beta-10/beta-100 roles, `R02`-`R09`. | Compare the existing admitted peak-alignment curves and test active/passive separation. |
| Figure 10: hybrid-kinetic comparison | External `Pegasus++` beta-16 result, not generated by CGL-LF AthenaK. | Block a full-paper reproduction claim until qualified source data or a compatible kinetic rerun is obtained; otherwise label it external context. |
| Figure 11: resolution/scale separation | `R16`, `R02`, `R17`. | Compare dimensionless alignment peaks; unblock energy-spectrum ordinate mapping; verify the scale at which alignment ends trends with dissipation scale. |
| Figure 12: heat-flux sensitivity | `R12`, `R02`, `R13`, plus nominal passive `R06`. | Compare dimensionless alignment peaks; unblock lower-panel gradient-spectrum mapping; test claimed insensitivity across factor-100 `|k_parallel|` variations. |
| Figure 13: limiter-induced collisionality | `R14`, `R15`, `R03`, with passive hard-wall context `R07` where plotted. | Compare admitted `beta Delta` PDF and normalized transfer curves; unblock energy/strain dimensional panels; test nonmonotonic response to `nu_lim`. |
| Figure 13 associated quantitative text | `R15` postprocessing with parallel and oblique firehose masks; outer-scale rate estimate; external kinetic context. | Test the reported `10.4%` versus `18.5%` unstable-volume reclassification over the same late-time interval and the estimated rate near `235 v_A/L_perp`; treat the reported hybrid-kinetic `17.9%` comparison as external unless its data are obtained. |
| Analytic/reduced-model conclusions and appendices | Equation audit plus reproduced diagnostics above. | Trace each claimed analytic prediction to an equation and a reproduced numerical test; explicitly separate derivation verification from independent simulation evidence. |

### 5.1 Theory and appendix audit

The numerical panel gate alone does not reproduce all of MKS24. The
derivations and stated analytic limits must be checked independently:

| MKS24 theory component | Required reproduction work | Acceptance record |
| --- | --- | --- |
| Main-text high-beta reduced CGL-MHD ordering and invariants | Re-derive the retained ordering, pressure-stress suppression, heat-flux scaling, and conserved quantities from the stated closure and conventions. | Equation-by-equation derivation notebook or manuscript appendix, with assumptions and any discrepancies listed. |
| Appendix: comparison with reduced kinetic MHD | Reproduce the reduction and identify exactly which magneto-immutability signatures survive and which non-local pressure effects do not. | Checked derivation linked to the main-text claim table. |
| Appendix: Braginskii viscous stress | Reproduce the collisional-limit argument and the criterion `beta k_parallel v_A/nu >> 1`, keeping it distinct from the threshold-activated `nu_lim` simulations. | Checked derivation plus terminology review against Stage I Section 3.4. |
| Appendix: larger density fluctuations | Reproduce the amended ordering and identify which conclusions depend on the small-density-fluctuation assumption. | Checked limitation statement connected to reproduced density diagnostics. |

The staged source audit identifies thirteen displayed result figures and
these three theoretical appendices. A complete claim must account for both
the panel gate above and this derivation audit.

### 5.2 Blocking reference-data work

The implementation record currently documents a real comparison boundary:
dimensionless products and a pressure-coordinate surface transform are
admitted, while several absolute spectral and strain ordinates are not yet
qualified because the source does not specify the discrete transform
normalization required for a faithful mapping to AthenaK products.

For a complete CGL-LF reproduction, those panels cannot remain silently
excluded. Before scientific manuscript completion, do at least one of:

1. Obtain machine-readable panel data or analysis code from the authors or an
   archival deposit.
2. Obtain an explicit normalization convention sufficient to implement and
   validate the transform.
3. Reanalyze source simulation dumps if made available with their analysis
   definitions.

If none is available, the manuscript may still report a carefully delimited
partial reproduction, but it must list the blocked panels and may not use
"complete reproduction."

## 6. Diagnostics, Visuals, and Acceptance Rules

### 6.1 Required MKS24 observables

| Observable family | Minimum retained products |
| --- | --- |
| Basic state | `rho`, `u`, `B`, `p_parallel`, `p_perp`, beta, `Delta p`, forcing work, and energy accounting. |
| Instability exposure | Mirror/firehose active fractions, limiter action, `nu_lim` policy, hard-wall projection/counters, and any invalid-state counters. |
| Pressure/density | Joint pressure-density PDFs, normalized density PDFs, and density fluctuation spectra. |
| Compressibility | Fourier-projected compressive velocity spectrum and its fraction of flow energy. |
| Cascade/pressure structure | Kinetic, magnetic, `p_parallel`, `p_perp`, magnetic-pressure, `Delta p`, and local parallel/perpendicular gradient spectra. |
| Energetics | `Delta p`-stress transfer with the MKS24 stated normalization, applied pressure/LF work, forcing work, and global accounting quality. |
| Organization | Scale-dependent rate-of-strain/local-field alignment PDFs and peak curves; local-field structure-function eddy anisotropy. |
| Reproduction comparison | Per-panel reference data, normalization metadata, residuals, uncertainty, and pass/fail/blocked status. |

### 6.2 Additional credibility diagnostics

These diagnostics strengthen the reproduction without changing the target
experiment:

| Addition | Reason |
| --- | --- |
| Low-field fractions `f_low(q)` and field-direction regularization counters | Protect local-field operations from hidden near-zero-field failures. |
| Block/bootstrap uncertainty over the declared late-time interval | Prevent one realization or burst from determining a result. |
| Exact forcing/pressure/LF work ledgers and restart identity | Distinguish physical stress transfer from implementation/accounting error. |
| Resolution and output-cadence sensitivity for selected products | Demonstrate that the plotted conclusion is not created by numerical or sampling choices. |

### 6.3 Visual-selection rule

For reproduction figures, first use the panel definition in MKS24. For any
additional AthenaK visual:

1. Select snapshots by a predeclared statistical rule, such as the saved
   output closest to the midpoint of the accepted averaging window.
2. Use fixed or pooled-percentile color limits documented in provenance.
3. Label supplementary event-focused snapshots as selected by an explicit
   criterion; do not substitute them for statistics.

### 6.4 Acceptance categories

Each panel/product must be assigned one of:

| Status | Meaning |
| --- | --- |
| `passed` | Setup and observable are matched, comparison uncertainty is declared, and AthenaK is quantitatively consistent with the reference criterion. |
| `failed` | A qualified comparison is inconsistent; retain it and revise the claim. |
| `blocked_reference` | A matching AthenaK product exists, but a defensible published-data transform/reference dataset is unavailable. |
| `external_model` | The published panel uses physics/code outside the AthenaK CGL-LF reproduction scope, such as Figure 10. |
| `not_run` | The required production simulation has not been executed. |

## 7. Stage I Implementation Work Before Production

The detailed record should remain the source of truth for completed
implementation evidence. For manuscript execution, the immediate required
work is:

| Work item | Existing starting point | Required completion evidence |
| --- | --- | --- |
| Freeze MKS24 case inventory | `inputs/cgl_lf_paper/mks24_stage_i_manifest.json` and guarded workflow `paper-mks24-stage-i` define sixteen source-mapped executions, reuse aliases, and exclude the unmapped active-Alfvenic beta-1 deck. | Preserve the committed manifest and validate it against every submitted input bundle. |
| Close reference-data boundary | Many dimensionless curves/surfaces are extracted; dimensional panels remain blocked. | Author/archive data, explicit transform proof, or manuscript-scoped blocked-panel decision. |
| Verify numerical model identity | Corrected collisions, limiter threshold selection, forcing work, passive feedback, and reduced/GPU gates are recorded. `docs/cgl_lf_mks24_stage_i_protocol_review.md` preserves the E02 stop line and records the corrected E03 disposition. | Preserve strict inspection for fresh mapped production from `t = 0`. |
| Replacement-driver qualification | Immutable E02 build revision `462b9dbd53e085dea46c2478b567781576d7d03e` remains pipeline and cost evidence only. Corrected E03 revision `9e07542281e4e6d125582f253df3ad2e3b8b154d` and executable SHA-256 `68f243f9204df388b24365ae65a567f6f567dbe422a6d7a43b9fb4a499ef118c` pass `g024`--`g031`: explicit planar/random policies, post-refresh restart identity, decomposition identity, passive semantics, nonlinear hard wall, and standard-layout sizing. Fresh R02 jobs `4745922`, `4746154`, `4746182`, `4746356`, `4746435`, `4746663`, `4746773`, `4746953`, `4747015`, `4747087`, `4747146`, `4747202`, `4747500`, `4747834`, `4748138`, `4752118`, `4753294`, and `4754008` are accepted through exact `t = 7.0`. | Preserve the retained token and continue only from the authenticated R02 checkpoint. |
| Execution-epoch isolation | `scripts/frontier/cgl_lf_stage_i.py` records new work beneath `runs/mks24-stage-i/E03-forcing-policy/`, uses `mks24_stage_i_E03_forcing_policy_*` accounting files, requires the reviewed token, tags manifests and jobs with the epoch, and rejects cross-epoch restart ancestry and bundle discovery. E01 and E02 trees remain outside E03 discovery. | Preserve passing offline isolation/recovery regressions and clean `reconcile` before every E03 submission. |
| Shared-root coordination | `/lustre/orion/ast207/proj-shared/dfielding/CGL` also contains exploratory non-MKS24 runs. A read-only audit found cancelled job `4743020` and two-node debug job `4743106` beneath `runs/beta25-accel05-*`; neither is Stage I evidence. Slurm now records `4743106` as `COMPLETED 0:0`, but its top-level campaign manifest remains stale at `state = running`. The helper therefore continues to require explicit acknowledgement after a read-only isolation review. | Before every Stage I prepare and submit action, inspect all user jobs and all active CGL-root campaign records. Do not modify the stale exploratory manifest. Do not overlap Stage I with another root-writing CGL campaign unless an explicit isolation and concurrency review authorizes it. Account exploratory runs separately from the MKS24 ledgers. |
| Retained source provenance | Both launchers require a checksummed in-root `--source-bundle` containing the prepared revision. The corrected qualified bundle is `athenak-feature-cgl-through-9e075422.bundle` with SHA-256 `c39d55809989d20aa5438711803f4fd43237fa9284c7e15183e84fdd693d3687`. F-080 controller transition `05cb4c324bfd8feec72ebdeb33b1961c9fde70bf` is archived as `athenak-feature-cgl-through-05cb4c324.bundle` with SHA-256 `1381918e471730d8c9639014475566f93fc87bac66b62738dd8316fc69c03570`; retained evidence JSON SHA-256 is `46fc1c4054e75f4224302be1895c1b9eaae88097544512ac067c53376e542e34`. Post-record lifecycle commits `9480e62764528a3f40066d22a192f0e99b369891` and `ef1e42fa088203ac9ef6ec8e47e668db4fb95a3c` require the live helper to remain committed during historical authentication, reserve case-level `analysis/`, and reject `--segment analysis`; focused Stage I tests pass (`10 passed`). Job `4746182` used verified launch bundle `athenak-feature-cgl-through-ef1e42fa.bundle` with SHA-256 `80c211c41a8de32688c3be577d8c727ec3268c347c2a7c3b2d5143e1c34593ec`. Job `4746356` used verified launch bundle `athenak-feature-cgl-through-377937639.bundle` with SHA-256 `91f8241b68247972f626c2210856a4e2c3cfe9ac9fed335b3c03f1131d2b2d68`. Job `4746435` used verified launch bundle `athenak-feature-cgl-through-a4ca8fb35.bundle` with SHA-256 `cda4f5009b63164fac390c8ae2abfabc8200fbda396e98add171a9c1e0c2dcaa`. Job `4746663` used verified launch bundle `athenak-feature-cgl-through-3fdbe6152.bundle` with SHA-256 `18eee2385721fd15cbbe9da8886cca6e0f4900901b6a86dcd87d3ca6659ffd85`. Job `4746773` used verified launch bundle `athenak-feature-cgl-through-bff57e4e6.bundle` with SHA-256 `91f22fa494412de766e83afcaa187d7d45168602ab8c498c98bcf4a8077c989b`. Job `4746953` used verified launch bundle `athenak-feature-cgl-through-6b16f3351.bundle` with SHA-256 `dc3ac356af07c0848825d456f151d11761b48cbc03ce1a1e9b30eeaeb5fb2644`. Job `4747015` used verified launch bundle `athenak-feature-cgl-through-cd52b8c55.bundle` with SHA-256 `b8dcf41f8714676e58e09af903364cf4e33408492d50c601cfdefd4b39ac1844`. Job `4747087` used verified launch bundle `athenak-feature-cgl-through-787424596.bundle` with SHA-256 `ddc1aa10443a08fb699b255772e777e2b2cb1782038eaf0320a36d37277369bc`. Job `4747146` used verified launch bundle `athenak-feature-cgl-through-eacb52678.bundle` with SHA-256 `a95d74aefbe1648c4633c05706e81baec7ae8bc75a807b81742e82e64d46e072`. Job `4747202` used verified launch bundle `athenak-feature-cgl-through-7669f2655.bundle` with SHA-256 `ea04bb9a4dacbc614b77070e97debfa763dea0cf0b2182e4c2222f2bd8a87bc2`. Job `4747500` used verified launch bundle `athenak-feature-cgl-through-8a0a406b6.bundle` with SHA-256 `b6e3188aba94f0f48e3d129b455bc28ea67e1cfac63eb55ba2101fc19e0c1fb5`. Job `4747834` used verified launch bundle `athenak-feature-cgl-through-b7d1828cf.bundle` with SHA-256 `e8d83fe47a3c72365dbbefbed8eceb10a326e85926f2033158a6fd7a0e4c3203`. F-083/F-084 hardening is promoted at canonical revision `9689c269bf329542815a1b2b137881126964b05c`, helper SHA-256 `1c633ebb58294938a0a0609742ae8f8d1f88cd15242ffe649578796edeb39375`: rendered lifecycle commands use the canonical absolute helper path; R02--R16 require one node; R17 requires eight nodes, accepted `t = 10` R02--R16 predecessors, and permanent lower-case lockout after R17 starts; replay validates reservation metadata, scheduler evidence, ledger arithmetic, cumulative ceilings, and positive finite thresholds. The helpers retain complete sibling restart sets, reject absent override targets before reservation, and reject queued user jobs plus unacknowledged active top-level records during submission preflight. `Prepared` and `submitted` manifests remain strict against live helper bytes; each `recorded` historical helper blob authenticates from its checksum-bound retained bundle. | Require helper authentication of retained corrected artifacts at every E03 lifecycle transition. |
| Production accounting path | `scripts/frontier/cgl_lf_stage_i.py` validates mapped cases and aliases, uses `batch` with default production `normal` QOS, requires the E03 approval token, archives executable/input/bundle/utility provenance, submits atomically under lock, enforces sequential reservations, rejects requests above two hours, authenticates restart markers, and requires formal inspection before record. E03 started empty by design and now records accepted R02 jobs `4745922`, `4746154`, `4746182`, `4746356`, `4746435`, `4746663`, `4746773`, `4746953`, `4747015`, `4747087`, `4747146`, `4747202`, `4747500`, `4747834`, `4748138`, `4752118`, `4753294`, and `4754008` through exact `t = 7.0` at `16.133334` cumulative node-hours. Hardened reconciliation reports `18/18/18` ledger rows/manifests/reservations, no active reservation, no transaction, and `issues = []`. | Before preparing a successor, commit these documentation bytes, archive and catalog the resulting actual post-F-100 source bundle and its SHA-256, rerun hardened reconciliation, and repeat the queue/shared-root audit with the explicit stale beta-25 acknowledgement. Then prepare only `R02/s18_rankio_t7_t7p25` on one node from the authenticated `s17` terminal restart with Slurm walltime `01:05:00`, Athena timeout `00:55:00`, and its reviewed one-segment `2700`-second threshold retaining `600` seconds on both timeout margins. The threshold is scoped only to `s18` and cannot ratchet automatically. Any scientific, provenance, scheduler, storage, budget, or reconciliation failure blocks successor preparation. Inspect, account, reconcile, and recost, then finish R02, execute R03--R16 sequentially, and execute `R17` last. |
| Production analysis/plot orchestration | Existing analyzer and MKS24 extractors implement many products; `paper-analyze` composes repeated checksum-qualified split reference manifests while rejecting duplicate product identifiers, and commit `449e297b` adds explicit partial-case comparison recording plus header-only snapshot-window selection; `cgl_lf_stage_i.py bundle-case` follows the explicit accepted restart lineage, merges sampled histories that need not repeat authenticated restart boundary rows while bounding gaps by their retained cadence and boundary timesteps, and links time-deduplicated shared or rank-local snapshot sets only after a case reaches its required final time, while `bundle-campaign` requires all sixteen completed mapped cases for cross-case comparisons. | Workflow command regenerating the per-panel gate table and retained figures from accepted output bundles. |
| Manuscript conversion | Illustrated validation note exists. | TeX manuscript structured around reproduction claims, blocked boundaries, and only later an extension. |

Job `4748138` used verified launch bundle
`athenak-feature-cgl-through-493fd3059.bundle` with SHA-256
`9e79a393a76baf15a9252bf96920e7c0f2a5067a956f28279c232a2f7e4b3623`.
Job `4752118` used verified launch bundle
`athenak-feature-cgl-through-368bc86e4.bundle` with SHA-256
`49a26705e9e34a203271e5a9e1b330d42bf97a8d64d1ad557d6edbc47ee76c68`.
Job `4753294` used verified launch bundle
`athenak-feature-cgl-through-6c0739806.bundle` with SHA-256
`a837ebe60d62a922dde3fc3ddd87350eec8a9e3f439f671f783346ee6290617f`.
Job `4754008` used verified launch bundle
`athenak-feature-cgl-through-e119e2dcf.bundle` with SHA-256
`25ffc9cd279a02debf8b55e318f0674ac3c14bf9a05cf05827218ae2056a5f14`.
The sole candidate `s18` launch requires a new committed post-F-100
documentation archive. Its actual bundle name
and SHA-256 must be recorded before
preparation; they are intentionally not guessed here.

No weak/strong-guide case family or uniform-`nu_coll` extension case should be
implemented as a production priority while these reproduction blockers remain.

## 8. Stage I Compute and Storage Budget

The user has authorized a clean planning reset if it improves the campaign.
The earlier revision adopted that reset for `E02-modal-driver` because the
modal turbulent-driver replacement required a new lineage from `t = 0`.
Historical allocation remains visible. The F-076 forcing-contract audit now
requires another fresh epoch from `t = 0`; no new `E02` submission is
authorized.

The live debug-helper ceiling remains `1000` node-hours and currently records
all qualification epochs. The earlier `1068.888889` Stage I reservation is
archived pre-replacement planning history; it does not authorize new
submissions and must not be silently copied into E03.

| Allocation record | Calculation basis | Node-hours |
| --- | --- | ---: |
| Historical debug qualification | Retained pre-replacement CGL-LF debug ledger | `0.851670` |
| Fresh replacement-driver debug qualification used so far | Retained `g014` through `g023`, including the failed preflight execution and inconclusive clipped-boundary comparator | `0.510280` |
| Corrected E03 qualification | Retained `g024` through `g031`, including explicit forcing policies, post-refresh restart identity, decomposition, passive semantics, nonlinear hard wall, and standard-layout sizing | `0.485834` |
| Current recorded debug total | Historical debug qualification plus E02 replacement-driver and corrected E03 debug use | `1.847784` |
| Accepted corrected E03 R02 pilot | Job `4745922`; one node; exact `t = 0.1`; formally inspected and recorded with no active reservation | `0.188889` |
| Accepted corrected E03 R02 native-boundary continuation | Job `4746154`; one node; exact `t = 0.25`; formally inspected and recorded with no active reservation | `0.289167` |
| Corrected E03 R02 prefix through exact `t = 0.25` | Jobs `4745922` and `4746154`; retained recost evidence JSON SHA-256 `eb6071cf53d453b2f75707003616ac9cbdfe025d62969d00e34e7948bee3310c` | `0.478056` |
| Accepted corrected E03 R02 continuation through exact `t = 1.0` | Job `4746182`; retained inspection SHA-256 `482c7f0f9552e319cc50e92f6d84f17565a90f9562b06f1d67ef1c57d240a2dd` | `1.465556` |
| Corrected E03 R02 prefix through exact `t = 1.0` | Jobs `4745922`, `4746154`, and `4746182`; retained superseded arithmetic-history recost evidence JSON SHA-256 `7d11e8a0004e24417167aeb1f186f8b1f68c8eb9d6206467a708a7cad16cf252` | `1.943612` |
| Accepted corrected E03 R02 continuation through exact `t = 1.5` | Job `4746356`; retained inspection SHA-256 `5a66f4e037fe08f176b51c7f949db0377bc9270e666ae1e1c69dda85d7f96633` | `1.048611` |
| Corrected E03 R02 prefix through exact `t = 1.5` | Jobs `4745922`, `4746154`, `4746182`, and `4746356`; retained corrected recost evidence JSON SHA-256 `b6fd6e8fcd939ca1baa06411a3cd6f04359b3bc2ff5b2a0097e2cc87b1763fc7` | `2.992223` |
| Accepted corrected E03 R02 continuation through exact `t = 2.0` | Job `4746435`; retained inspection SHA-256 `fee92e6d8b484d1050266f6d5600577d19099df02335299114da59dd42050606` | `1.063611` |
| Corrected E03 R02 prefix through exact `t = 2.0` | Jobs `4745922`, `4746154`, `4746182`, `4746356`, and `4746435`; retained corrected recost evidence JSON SHA-256 `c382b78ab9e21466648adf9d7ea38d2b407f0e986579f128bdfe6bc44bef79dd` | `4.055834` |
| Accepted corrected E03 R02 continuation through exact `t = 2.5` | Job `4746663`; retained inspection SHA-256 `57fb0b9f1a19c0125aec0619802f2c41bd91ac97c51f32e70ceae9a434355da4` | `1.062500` |
| Corrected E03 R02 prefix through exact `t = 2.5` | Jobs `4745922`, `4746154`, `4746182`, `4746356`, `4746435`, and `4746663`; retained corrected recost evidence JSON SHA-256 `39115585cf7ca866e3aade1e53c69aab6ee9217361d75db7946748905e1c0e5c` | `5.118334` |
| Accepted corrected E03 R02 continuation through exact `t = 3.0` | Job `4746773`; retained inspection SHA-256 `4b32798174222c9b9ca165285397ed9a02a411bb3e3dbf09e1f10611eebc70db` | `1.077500` |
| Corrected E03 R02 prefix through exact `t = 3.0` | Jobs `4745922`, `4746154`, `4746182`, `4746356`, `4746435`, `4746663`, and `4746773`; retained corrected recost evidence JSON SHA-256 `d2c5e27553426beb03c8c1882bd6473b682029c6ba67acd32bb06c5fb9a3ff01` | `6.195834` |
| Accepted corrected E03 R02 continuation through exact `t = 3.5` | Job `4746953`; retained inspection SHA-256 `4714af847c1ce10517258ea34f6151eb2a60d7fbd1bf1717d2de38d2e6822933` | `1.151111` |
| Corrected E03 R02 prefix through exact `t = 3.5` | Jobs `4745922`, `4746154`, `4746182`, `4746356`, `4746435`, `4746663`, `4746773`, and `4746953`; retained corrected recost evidence JSON SHA-256 `468e787db48d98661cd3ffe125829a78f135bb7ca4c9104c766a00b693f586a3` | `7.346945` |
| Accepted corrected E03 R02 continuation through exact `t = 4.0` | Job `4747015`; retained inspection SHA-256 `be705544f6c8f4a71f5501130b1e317555cdf1cac0663fba267acec82d967649` | `1.200000` |
| Corrected E03 R02 prefix through exact `t = 4.0` | Jobs `4745922`, `4746154`, `4746182`, `4746356`, `4746435`, `4746663`, `4746773`, `4746953`, and `4747015`; retained corrected recost evidence JSON SHA-256 `dd5b0f5b82cf81005c1a481ee83d1127aa43861dcdd8c22765170e28b0ad6c99` | `8.546945` |
| Accepted corrected E03 R02 continuation through exact `t = 4.5` | Job `4747087`; retained inspection SHA-256 `0373fa3505beaadbd818db5271dae13bf69a41b1683bb561d558d8d71eff7581` | `1.208056` |
| Corrected E03 R02 prefix through exact `t = 4.5` | Jobs `4745922`, `4746154`, `4746182`, `4746356`, `4746435`, `4746663`, `4746773`, `4746953`, `4747015`, and `4747087`; retained corrected recost evidence JSON SHA-256 `8b10b1786db520740b687d989207f3cda6e0123f9bf2b62e75a5df3849511561` | `9.755001` |
| Accepted corrected E03 R02 continuation through exact `t = 5.0` | Job `4747146`; retained inspection SHA-256 `8fb8069dc1f75448a8f79a38de98141c5d6dfed83bc3530eb633f8c009655014` | `1.226389` |
| Corrected E03 R02 prefix through exact `t = 5.0` | Jobs `4745922`, `4746154`, `4746182`, `4746356`, `4746435`, `4746663`, `4746773`, `4746953`, `4747015`, `4747087`, and `4747146`; retained corrected recost evidence JSON SHA-256 `5256fb7fd9b75f175a16c2c16ee8a995b8d6d9281a92b70a979ea62171835bf2` | `10.981390` |
| Accepted corrected E03 R02 continuation through exact `t = 5.5` | Job `4747202`; retained inspection SHA-256 `928c281788b8ba359eb499aa607529b6c8eba2f8236e374f808d57fe239d9c1e` | `1.268889` |
| Corrected E03 R02 prefix through exact `t = 5.5` | Jobs `4745922`, `4746154`, `4746182`, `4746356`, `4746435`, `4746663`, `4746773`, `4746953`, `4747015`, `4747087`, `4747146`, and `4747202`; retained corrected recost evidence JSON SHA-256 `323d7cb3cbe108eee944045c189cff75ffc07d06833c49bcfa8315f07fad8d51` | `12.250279` |
| Accepted corrected E03 R02 quarter-unit continuation through exact `t = 5.75` | Job `4747500`; retained inspection SHA-256 `24dbc74639d09222b321f805af0d47adf5b917f2cd83d765bb75f69ddadb5680`; retained independent validation SHA-256 `7d27cd3a67eb06fdde1b96cb0c5809fef524eeb08d71a745fca2927f0d203a0e` | `0.660000` |
| Corrected E03 R02 prefix through exact `t = 5.75` | Jobs `4745922`, `4746154`, `4746182`, `4746356`, `4746435`, `4746663`, `4746773`, `4746953`, `4747015`, `4747087`, `4747146`, `4747202`, and `4747500`; retained corrected recost evidence JSON SHA-256 `d0a60e222138971e1bc9aae978bb3d62aee460f09000ee62f7f9ffdfe534165f` | `12.910279` |
| Accepted corrected E03 R02 quarter-unit continuation through exact `t = 6.0` | Job `4747834`; retained inspection SHA-256 `91cf8bbd7f9e1592ae987f92bd296a5d2c66e7fcc941ddecd6c05b3e046c9b29`; retained independent validation SHA-256 `fdc05bac65e860a1ac5dc3c6f5a9b8f045b54a9fd421c6ef8030fc902008a448` | `0.599722` |
| Corrected E03 R02 prefix through exact `t = 6.0` | Jobs `4745922`, `4746154`, `4746182`, `4746356`, `4746435`, `4746663`, `4746773`, `4746953`, `4747015`, `4747087`, `4747146`, `4747202`, `4747500`, and `4747834`; retained corrected recost evidence JSON SHA-256 `4bc19f5b44587e73869982425b1cf0ad459547c3e0c5de1ab1b2b5cb366c6322` | `13.510001` |
| Accepted corrected E03 R02 quarter-unit continuation through exact `t = 6.25` | Job `4748138`; retained inspection SHA-256 `e9979a42e3e7572b233ce788b4c79a8512932746fc6c9cd46ffacc72e233854a`; retained independent validation SHA-256 `f7068e64ba870a93f6b525e7168e7fa815d9364ade01b870cc3cf39b532ffb13` | `0.655278` |
| Corrected E03 R02 prefix through exact `t = 6.25` | Jobs `4745922`, `4746154`, `4746182`, `4746356`, `4746435`, `4746663`, `4746773`, `4746953`, `4747015`, `4747087`, `4747146`, `4747202`, `4747500`, `4747834`, and `4748138`; retained corrected recost evidence JSON SHA-256 `5ba7aec52d3e7824cf4ccf52f95ce1e46368c84ed59f0525e2bbbcfc167fde91` | `14.165279` |
| Accepted corrected E03 R02 quarter-unit continuation through exact `t = 6.5` | Job `4752118`; retained inspection SHA-256 `fd1fa6f3e4218b611814b71a9081140d91f2c1d866a2f6f565de33ec040c7492`; retained independent validation SHA-256 `cba6a5ff33e7c5920874e78e55c8288e7863de43b80a1accfa018d45e8547968` | `0.653611` |
| Corrected E03 R02 prefix through exact `t = 6.5` | Jobs `4745922`, `4746154`, `4746182`, `4746356`, `4746435`, `4746663`, `4746773`, `4746953`, `4747015`, `4747087`, `4747146`, `4747202`, `4747500`, `4747834`, `4748138`, and `4752118`; retained corrected recost evidence JSON SHA-256 `3603077400bd4530001155b83bfc4ac1053ec18505c56a6c29ffadc795bd1255` | `14.818890` |
| Accepted corrected E03 R02 quarter-unit continuation through exact `t = 6.75` | Job `4753294`; retained inspection SHA-256 `1149de945447eae3b274052985363c1c821720a117254662657e33d8b21467de`; retained independent validation SHA-256 `4889d08083004d3d89f6cc8beede3521ce7f8143983c7b4bb233209d59844b05` | `0.647222` |
| Corrected E03 R02 prefix through exact `t = 6.75` | Jobs `4745922`, `4746154`, `4746182`, `4746356`, `4746435`, `4746663`, `4746773`, `4746953`, `4747015`, `4747087`, `4747146`, `4747202`, `4747500`, `4747834`, `4748138`, `4752118`, and `4753294`; retained corrected recost evidence JSON SHA-256 `5b997623a3f8d83c200034f42e0c9f1b9a09c811bcb46e2fc1a0f7c1eb58815e` | `15.466112` |
| Accepted corrected E03 R02 quarter-unit continuation through exact `t = 7.0` | Job `4754008`; retained inspection SHA-256 `f5028f8d6c6157e837b34e06712d68360040bfebdfb5bb1bfdc0221cd8431db8`; retained independent validation SHA-256 `af4e188ae27904e41b6bf7d517784688ac2b46be03de70d7285c96a4b1e8ccdd` | `0.667222` |
| Corrected E03 R02 prefix through exact `t = 7.0` | Jobs `4745922`, `4746154`, `4746182`, `4746356`, `4746435`, `4746663`, `4746773`, `4746953`, `4747015`, `4747087`, `4747146`, `4747202`, `4747500`, `4747834`, `4748138`, `4752118`, `4753294`, and `4754008`; retained corrected recost evidence JSON SHA-256 `63438d8843b86815a2e25d4d6de6ef78ce9756c561fe9db3e762965ada73202f` | `16.133334` |
| Separately retained scheduled MPI CPU local-gate allocation | Job `4743666`; not charged through the debug-helper ledger | `0.002222` |
| Historical `E01-pre-modal-driver` Stage I | Retained production segments through rejected job `4686032` | `9.962778` |
| Pre-E02 historical subtotal | Historical debug qualification plus `E01-pre-modal-driver` Stage I; reported for provenance, not charged against the fresh planning reset | `10.814448` |
| Historical `E02-modal-driver` project ceiling | Superseded incremental ceiling for archived replacement-driver qualification and Stage I work | `4000.000000` |
| Historical initial `E02` replacement-driver qualification envelope | Smallest-discriminating local and Frontier debug reruns before archived production recosting | `50.000000` |
| Archived `E02` `R16` pipeline run | Seven accepted jobs through exact `t = 10.0`; zero strict counters; completed analyzer pipeline check | `6.145556` |
| Historical `E02` `R16` completion reservation | F-068 envelope; completed R16 left `0.854444` node-hours unused | `7.000000` |
| Archived `E02` `R02` standard-layout timing pilot | F-069 bounded measurement; two accepted jobs through native snapshot boundary `t = 0.25` | `0.473333` |
| Archived `E02` `R17` high-resolution timing pilot | F-071 bounded eight-node measurement; two accepted jobs through authenticated continuation `t = 0.10` | `4.235556` |
| Historical `E02` Stage I mapped-matrix envelope | Superseded F-071 measured projection `697.019444` plus `202.980556` node-hours of margin; only frozen `R02`--`R17` cases | `900.000000` |
| Archived `E02` `R02` mapped continuation through `t = 1.0` | F-072 longer developed standard-layout interval; updated matrix projection `702.195370` leaves `197.804630` margin | `1.452778` |
| Archived `E02` `R02` mapped continuation from `t = 1.0` through `t = 1.5` | F-073 late-time standard-layout interval; updated matrix projection `725.464938` leaves `174.535062` margin | `1.052222` |
| Archived `E02` `R02` mapped continuation from `t = 1.5` through `t = 2.0` | F-074 late-time standard-layout interval; updated matrix projection `730.081697` leaves `169.918303` margin | `1.068889` |
| Archived `E02` `R02` mapped continuation from `t = 2.0` through `t = 2.5` | F-075 late-time standard-layout interval; updated matrix projection `728.164877` leaves `171.835123` margin | `1.061944` |

The archived estimate is not a measured replacement-driver benchmark. Preserve
the `E02` archive and carry its operational lessons into the corrected fresh
epoch:

1. Preserve the committed bundle-backed source-provenance, shared-root, and
   absent-override-key controls through `7fae0bcf`.
2. Preserve completed `R16` and retained completion-recost evidence JSON
   SHA-256 `d42a0432858127a94ab4aee5aab0c0b79572b70730085352aec53147670adad3`.
3. Preserve completed R02 timing-pilot evidence JSON SHA-256
   `81c9f6e56cca9eb5a69449a1166bf5d3e3cc8a0a4102b6b4fd0c3f73f3c2123b`.
4. Preserve completed R17 timing-pilot evidence JSON SHA-256
   `79147715729a2c9017400e170de2191a3f7d6b0bf411acee09a948ada464b822`.
   Use its authenticated developed continuation rate and measured rank-local
   groups rather than carrying forward a cell/CFL bracket.
5. Preserve accepted R02 `t = 1.0` continuation evidence JSON SHA-256
   `7ba577644cd65eefb4c9e3a782f284c36be9e86894af1ec17e6242eba1c17b06`.
   Reject normal-QOS requests above two hours during preparation.
6. Preserve accepted R02 `t = 1.5` continuation evidence JSON SHA-256
   `0319284b72628f2456087f378588515726f0fa8586b20660d37da2a20a1ba6ba`.
   Use the measured late-time rate for conservative standard-layout recosting.
7. Preserve accepted R02 `t = 2.0` continuation evidence JSON SHA-256
   `d0c7e0f9cdd8031a9d88c195adee9822943a358eb55193b0e1146bf3b95a08e4`.
   Continue to use the most recent measured late-time rate conservatively.
8. Preserve accepted R02 `t = 2.5` continuation evidence JSON SHA-256
   `f823ef6993045d5d90f016f7052518a353f348f7f83b9efdc19d32ee097f0c59`.
   Continue to use the most recent measured late-time rate conservatively.
9. Prepare or submit no additional `E02` segment.
10. Preserve the completed corrected-policy qualification and retained
    approval token. Fresh mapped work began from `R02/s00_rankio_t0_t0p1`
    job `4745922` at `t = 0`; authenticated continuation job `4746154` is
    accepted through exact `t = 0.25`. The corrected prefix uses `0.478056`
    node-hours. Retained recost evidence JSON SHA-256 is
    `eb6071cf53d453b2f75707003616ac9cbdfe025d62969d00e34e7948bee3310c`.
    Authenticated continuation job `4746182` is also accepted through exact
    `t = 1.0`; the corrected prefix now uses `1.943612` node-hours. Updated
    recost evidence JSON SHA-256
    `7d11e8a0004e24417167aeb1f186f8b1f68c8eb9d6206467a708a7cad16cf252`
    is preserved as superseded arithmetic history. Authenticated continuation
    job `4746356` is accepted through exact `t = 1.5`; the corrected prefix
    now uses `2.992223` node-hours. Retained corrected recost evidence JSON
    SHA-256 is
    `b6fd6e8fcd939ca1baa06411a3cd6f04359b3bc2ff5b2a0097e2cc87b1763fc7`.
    Authenticated continuation job `4746435` is accepted through exact
    `t = 2.0`; the corrected prefix now uses `4.055834` node-hours. Retained
    corrected recost evidence JSON SHA-256 is
    `c382b78ab9e21466648adf9d7ea38d2b407f0e986579f128bdfe6bc44bef79dd`.
    Authenticated continuation job `4746663` is accepted through exact
    `t = 2.5`; the corrected prefix now uses `5.118334` node-hours. Retained
    corrected recost evidence JSON SHA-256 is
    `39115585cf7ca866e3aade1e53c69aab6ee9217361d75db7946748905e1c0e5c`.
    Authenticated continuation job `4746773` is accepted through exact
    `t = 3.0`; the corrected prefix now uses `6.195834` node-hours. Retained
    corrected recost evidence JSON SHA-256 is
    `d2c5e27553426beb03c8c1882bd6473b682029c6ba67acd32bb06c5fb9a3ff01`.
    Authenticated continuation job `4746953` is accepted through exact
    `t = 3.5`; the corrected prefix now uses `7.346945` node-hours. Retained
    corrected recost evidence JSON SHA-256 is
    `468e787db48d98661cd3ffe125829a78f135bb7ca4c9104c766a00b693f586a3`.
    Authenticated continuation job `4747015` is accepted through exact
    `t = 4.0`; the corrected prefix now uses `8.546945` node-hours. Retained
    corrected recost evidence JSON SHA-256 is
    `dd5b0f5b82cf81005c1a481ee83d1127aa43861dcdd8c22765170e28b0ad6c99`.
    Authenticated continuation job `4747087` is accepted through exact
    `t = 4.5`; the corrected prefix now uses `9.755001` node-hours. Retained
    corrected recost evidence JSON SHA-256 is
    `8b10b1786db520740b687d989207f3cda6e0123f9bf2b62e75a5df3849511561`.
    Authenticated continuation job `4747146` is accepted through exact
    `t = 5.0`; the corrected prefix now uses `10.981390` node-hours. Retained
    corrected recost evidence JSON SHA-256 is
    `5256fb7fd9b75f175a16c2c16ee8a995b8d6d9281a92b70a979ea62171835bf2`.
    Authenticated continuation job `4747202` is accepted through exact
    `t = 5.5`; the corrected prefix now uses `12.250279` node-hours. Retained
    corrected recost evidence JSON SHA-256 is
    `323d7cb3cbe108eee944045c189cff75ffc07d06833c49bcfa8315f07fad8d51`.
    Authenticated quarter-unit continuation job `4747500` is accepted through
    exact `t = 5.75`; the corrected prefix now uses `12.910279` node-hours.
    Retained corrected recost evidence JSON SHA-256 is
    `d0a60e222138971e1bc9aae978bb3d62aee460f09000ee62f7f9ffdfe534165f`.
    Authenticated quarter-unit continuation job `4747834` is accepted through
    exact `t = 6.0`; the corrected prefix now uses `13.510001` node-hours.
    Retained corrected recost evidence JSON SHA-256 is
    `4bc19f5b44587e73869982425b1cf0ad459547c3e0c5de1ab1b2b5cb366c6322`.
    Authenticated quarter-unit continuation job `4748138` is accepted through
    exact `t = 6.25`; the corrected prefix now uses `14.165279` node-hours.
    Retained corrected recost evidence JSON SHA-256 is
    `5ba7aec52d3e7824cf4ccf52f95ce1e46368c84ed59f0525e2bbbcfc167fde91`.
    Authenticated quarter-unit continuation job `4752118` is accepted through
    exact `t = 6.5`; the corrected prefix now uses `14.818890` node-hours.
    Retained corrected recost evidence JSON SHA-256 is
    `3603077400bd4530001155b83bfc4ac1053ec18505c56a6c29ffadc795bd1255`.
    Authenticated quarter-unit continuation job `4753294` is accepted through
    exact `t = 6.75`; the corrected prefix now uses `15.466112` displayed
    node-hours (`55678 / 3600 = 15.466111111111111` exact node-hours).
    Retained corrected recost evidence JSON SHA-256 is
    `5b997623a3f8d83c200034f42e0c9f1b9a09c811bcb46e2fc1a0f7c1eb58815e`.
    Authenticated quarter-unit continuation job `4754008` is accepted through
    exact `t = 7.0`; the corrected prefix now uses
    `16.133334` displayed node-hours
    (`58080 / 3600 = 16.133333333333333`
    exact node-hours). Retained corrected recost evidence JSON SHA-256 is
    `63438d8843b86815a2e25d4d6de6ef78ce9756c561fe9db3e762965ada73202f`.
11. In that fresh epoch, execute no more than one production segment at a time,
    review strict counters, work accounting, storage, and measured cost after
    every segment, finish R02, execute R03--R16 sequentially, complete R17
    last, and stop to recost if projected Stage I consumption exceeds its
    reserved ceiling.

The archived estimated sequential-retention storage allocation, conservatively
including the excluded beta-1 inventory definition only as an unexecuted contingency, is approximately
`880.368 GB` with margin; the no-pruning alternative is approximately
`1140.924 GB`. Those values are planning history, not an `E02` reservation:
they predate rank-local production retention and the modal restart payload.
Completed R16 retains `4686854064` raw segment-tree snapshot/restart bytes
(`4.365` GiB); its deduplicated 41-snapshot/11-restart cadence is
`4514126544` bytes. R02 terminal groups measure `509682424` snapshot bytes
and `1381366304` restart bytes. They project `36092008728` bytes (`33.613`
GiB) for one standard-layout case at the same cadence. R17 measures maximum
rank-local snapshot and restart groups of `4077462784` and `11052866880`
bytes. The measured sixteen-case logical raw estimate is `798559758560`
bytes (`743.717` GiB). Running R17 last, retaining accepted snapshots and
only the final two restart groups, yields a sequential peak of `622953051824`
bytes and an F-071 storage reservation of `747543662189` bytes after a
20-percent margin. The no-pruning alternative with the same margin is
`958271710272` bytes.
The contingency does not authorize a beta-1 production run absent a
source-mapped role.

### 8.1 Archived `E01-pre-modal-driver` Production Record

The following rows are immutable historical execution provenance, not
replacement-driver reproduction evidence:

| Case/segment | Frontier job | Requested allocation | Input/build provenance | Acceptance gate |
| --- | ---: | --- | --- | --- |
| `R16/s00_t0_t2` | `4674731`, `COMPLETED`, `0:0` | `1` node, `01:15:12` charged elapsed; `1.253333` actual node-hours | Legacy shared-MPI-I/O input at committed revision `e129a7ff8be60eb177c8268da0a73e7700d8b543`; executable SHA-256 `de99066f3546976093f6ff656660952032e30f31c939736c0212629cfb015157`; archived submitted input and matrix in the run manifest directory. | `clean_partial`: terminal `t = 1.66204848945 < 2.0`, zero strict LF safety counters and retained terminal restart; retained as diagnostic evidence only and not continued into the changed output-protocol digest. |
| `R16/s01_rankio_t0_t2` | `4675322`, `COMPLETED`, `0:0` | `1` node, `01:15:11` charged elapsed; `1.253056` actual node-hours | Rank-local retained-output input/utility revision `1f52784d291718fff3bf1dc984a682ada52780ab`; executable revision `e129a7ff8be60eb177c8268da0a73e7700d8b543` and SHA-256 `de99066f3546976093f6ff656660952032e30f31c939736c0212629cfb015157`; complete `t = 1` ranked snapshot/restart boundary and terminal ranked products retained, initial snapshot read through the rank-aware analyzer, and no legacy-scale output pause observed. | `clean_partial`: terminal `t = 1.87769834013 < 2.0`, zero strict LF safety counters and complete terminal rank-local snapshot/restart products; eligible for same-digest continuation only. |
| `R16/s02_rankio_t1p877698_t2` | `4676557`, `COMPLETED`, `0:0` | `1` node, `00:06:23` charged elapsed; `0.106389` actual node-hours | Continued from inspected terminal eight-rank restart set of `s01` with unchanged input/executable digests. | `accepted`: reached `t = 2.0`, zero strict LF safety counters, and complete terminal rank-local snapshot/restart products; establishes the ranked prefix gate only. |
| `R16/s03_rankio_t2_t3p5` | `4676696`, `COMPLETED`, `0:0` | `1` node, `01:23:01` charged elapsed; `1.383611` actual node-hours | Continued from inspected accepted terminal eight-rank restart set of `s02` with unchanged input/executable digests; terminal outputs contain five complete eight-rank snapshot groups at `t = 2.50014742786`, `2.75033411668`, `3.00023812685`, `3.25023774201`, and `3.5`, plus the terminal eight-rank restart group. | `accepted`: reached `t = 3.5`, formal inspection records zero strict LF safety counters and complete terminal rank-local products; establishes a longer ranked prefix gate only. |
| `R16/s04_rankio_t3p5_t5` | `4679803`, `COMPLETED`, `0:0` | `1` node, `01:34:04` charged elapsed; `1.567778` actual node-hours | Continued from the inspected accepted terminal eight-rank `.00004.rst` set of `s03` with unchanged input/executable digests. Terminal outputs contain six complete eight-rank snapshot groups at `t = 3.75032903725`, `4.00000477772`, `4.25028947444`, `4.49999984569`, `4.75028920409`, and `5.0`, plus the terminal eight-rank restart group. | `accepted`: reached `t = 5.0`, formal inspection records zero strict LF safety counters and complete terminal rank-local products; establishes a longer ranked prefix gate only. |
| `R16/s05_rankio_t5_t6p5` | `4681476`, `COMPLETED`, `0:0` | `1` node, `01:45:12` charged elapsed; `1.753333` actual node-hours | Continued from the inspected accepted terminal eight-rank `.00005.rst` set of `s04` using unchanged input/executable digests. Terminal outputs contain six complete eight-rank snapshot groups at the `t = 5.25`, `5.5`, `5.75`, `6.0`, `6.25`, and terminal `6.43828031751` states, plus complete scheduled `t = 6.0` and terminal restart groups. | `clean_partial`: clean application-walltime termination at `t = 6.43828031751 < 6.5`; formal inspection records zero strict LF safety counters and complete terminal ranked products, allowing same-digest continuation only. |
| `R16/s06_rankio_t6p438280_t6p5` | `4683291`, `COMPLETED`, `0:0` | `1` node, `00:04:50` charged elapsed; `0.080556` actual node-hours | Continued from the inspected clean-partial terminal eight-rank `.00007.rst` set of `s05` using unchanged input/executable digests and a `time/tlim=6.5` override; terminal outputs include complete eight-rank snapshot and restart groups at `t = 6.5`. | `accepted`: reached `t = 6.5`, formal inspection records zero strict LF safety counters and complete terminal rank-local products; promotes the accepted ranked prefix through the clean-partial parent. |
| `R16/s07_rankio_t6p5_t7p5` | `4683918`, `COMPLETED`, `0:0` | `1` node, `01:16:55` charged elapsed; `1.281944` actual node-hours | Continued from the inspected accepted terminal eight-rank `.00008.rst` set of `s06` using unchanged input/executable digests and a `time/tlim=7.5` override; terminal outputs contain three complete eight-rank snapshot groups at `t = 7.00026203265`, `7.25019568336`, and `7.5`, plus the terminal eight-rank restart group. | `accepted`: reached `t = 7.5`, formal inspection records zero strict LF safety counters and complete terminal rank-local products; extends the accepted ranked prefix but remains before the analysis window. |
| `R16/s08_rankio_t7p5_t8p5` | `4686032`, `COMPLETED`, `0:0` | `1` node, `01:16:58` charged elapsed; `1.282778` actual node-hours | Ran from the inspected accepted terminal eight-rank `.00009.rst` set of `s07` under pre-replacement executable revision `e129a7ff8be60eb177c8268da0a73e7700d8b543` and SHA-256 `de99066f3546976093f6ff656660952032e30f31c939736c0212629cfb015157`; terminal products reached exact `t = 8.5`. | `rejected`: formal inspection records zero strict LF safety counters and complete terminal ranked products, but the merged turbulent-driver replacement changes forcing/restart state, so this old-executable continuation is archival evidence only. |

The ranked `t = 0`--`2` prefix was assembled as the diagnostic bundle
`runs/bundles/diagnostic-prefixes/R16_rankio_t2_prefix_20260525` after the
bundler was made lineage-aware, then reassembled successfully under the
cadence-bounded implementation as
`runs/bundles/diagnostic-prefixes/R16_rankio_t2_prefix_cadence_20260525`.
Both bundles contain only `s01` and `s02` and link ten retained rank-local
snapshots. The first bundle's default `paper-analyze`
execution produced only a history summary because a `t = 2` prefix has no
samples in the production `t = 8`--`10` analysis window. A distinct,
infrastructure-only transient-window check over `t = 0`--`2` then read all
ten snapshots, selected 101 history rows, passed the synthetic analyzer
check, and rendered nine diagnostic figures. It is not a final-time `R16`
bundle and is not panel reproduction evidence.

After `s03` acceptance, the lineage-checked bundle
`runs/bundles/diagnostic-prefixes/R16_rankio_t3p5_prefix_cadence_20260525`
was assembled from only `s01`, `s02`, and `s03`; its manifest records
`accepted_final_time = 3.5`, `status = accepted_for_analysis`, and fifteen
deduplicated rank-local snapshots. It verifies accepted-prefix assembly at
the new gate, but it still has no samples in the required production window
`t = 8`--`10`: a default-window analyzer check found `176` merged history
rows overall but selected zero history rows and zero snapshots in that
window, while passing its synthetic diagnostic test. It is not panel
reproduction evidence.

After `s04` acceptance, the further lineage-checked bundle
`runs/bundles/diagnostic-prefixes/R16_rankio_t5_prefix_cadence_20260525`
was assembled from only `s01` through `s04`; its manifest records
`accepted_final_time = 5.0`, `status = accepted_for_analysis`, and twenty-one
deduplicated rank-local snapshots. A default-window analyzer check found
`251` merged history rows overall but selected zero history rows and zero
snapshots in `t = 8`--`10`, while passing its synthetic diagnostic test.
This is accepted continuation evidence, not panel reproduction evidence.

After `s06` acceptance, including the inspected `s05` clean partial in its
continued lineage, the accepted-prefix bundle
`runs/bundles/diagnostic-prefixes/R16_rankio_t6p5_prefix_cadence_20260526`
was assembled from `s01` through `s06`; its manifest records
`accepted_final_time = 6.5`, `status = accepted_for_analysis`, and
twenty-eight deduplicated ranked snapshots. A default-window analyzer check
found `326` merged history rows overall but selected zero history rows and
zero snapshots in `t = 8`--`10`, while passing its synthetic diagnostic test.
This is accepted continuation evidence, not panel reproduction evidence.

After `s07` acceptance, the accepted-prefix bundle
`runs/bundles/diagnostic-prefixes/R16_rankio_t7p5_prefix_cadence_20260526`
was assembled from `s01` through `s07`; its manifest records
`accepted_final_time = 7.5`, `status = accepted_for_analysis`, and
thirty-one deduplicated ranked snapshots. A generic default-window analysis
found `376` merged history rows overall but selected zero history rows and
zero snapshots in `t = 8`--`10`, while passing its synthetic diagnostic
test and writing one history-only figure. A requested Figure 11 reference
comparison fails closed because the present case selects an empty analysis;
this prefix is not panel reproduction evidence.

Commit `449e297b` adds an explicitly scoped
`--allow-partial-reference-cases` comparison mode, which records products
omitted because their cases are not present in a single-case bundle, and
selects snapshot windows from binary headers before loading full fields. A
local terminal-snapshot plumbing probe at `t = 2.0` evaluated only the
available Figure 11 `R16` product and recorded the absent `R02` and `R17`
products as omissions. That transient, out-of-window probe is not production
comparison evidence.

Job `4683918` is inspected and recorded as an accepted continuation through
`t = 7.5` after the clean partial `4681476` at `t = 6.43828031751` and
accepted terminal-gate segment `4683291`. Job `4686032` reached `t = 8.5`
and passed inspection, but was recorded `rejected` after the replacement
turbulent driver changed its forcing/restart-state contract. Do not continue
the old restart lineage or describe `R16` as a completed accepted case;
replacement-driver production must begin from `t = 0` after qualification
and immutable-build archival.

### 8.2 Archived `E02-modal-driver` Pipeline and Cost Record

The following rows are archived replacement-driver pipeline and cost
evidence. F-076 prohibits their use as direct MKS24 reproduction evidence and
prohibits any new E02 preparation or submission:

| Case/segment | Frontier job | Requested allocation | Input/build provenance | Acceptance gate |
| --- | ---: | --- | --- | --- |
| `R16/s00_rankio_t0_t2` | `4743735`, `COMPLETED`, `0:0` | `1` node, `00:59:43` charged elapsed; `0.995278` actual node-hours | Fresh `E02-modal-driver` rank-local run from `t = 0`; source bundle `athenak-feature-cgl-through-d93abefb.bundle` SHA-256 `4a3653afdf8394daad5b4c2f7c0542c5017437d89d60bcbb91e51abde0e6240a`; immutable executable SHA-256 `df87684e9d2b7af33b36c2757779d15de87b84ef051ca9c4489f2efa358f5c48`. | `accepted`: exact `t = 2.0`, zero strict LF failure counters, nine complete eight-rank snapshot groups, and three complete modal-restart groups. |
| `R16/s01_rankio_t2_t3p5` | `4743933`, `COMPLETED`, `0:0` | `1` node, `00:52:15` charged elapsed; `0.870833` actual node-hours | Continued from the inspected terminal eight-rank modal restart set of `s00` with matching input and executable digests. | `accepted`: exact `t = 3.5`, zero strict LF failure counters, six complete eight-rank snapshot groups, and two complete modal-restart groups. |
| `R16/s02_rankio_t3p5_t5` | `4743988`, `COMPLETED`, `0:0` | `1` node, `00:57:09` charged elapsed; `0.952500` actual node-hours | Continued from inspected `s01` modal restart siblings with source bundle `athenak-feature-cgl-through-80fa53cc.bundle` SHA-256 `ecd5c02d1e21a712625b45419b8c2d87f0542c283dadb38848bf5baa82fa143d`. | `accepted`: exact `t = 5.0`, zero strict counters, and complete terminal ranked products. |
| `R16/s03_rankio_t5_t6p5` | `4744019`, `COMPLETED`, `0:0` | `1` node, `00:59:20` charged elapsed; `0.988889` actual node-hours | Continued from inspected `s02` modal restart siblings with unchanged input/executable digests. | `accepted`: exact `t = 6.5`, zero strict counters, and complete terminal ranked products. |
| `R16/s04_rankio_t6p5_t8` | `4744056`, `COMPLETED`, `0:0` | `1` node, `00:59:13` charged elapsed; `0.986944` actual node-hours | Continued from inspected `s03` modal restart siblings with unchanged input/executable digests. | `accepted`: exact `t = 8.0`, zero strict counters, and complete terminal ranked products; reaches the lower analysis-window boundary only. |
| `R16/s05_rankio_t8_t9p5` | `4744120`, `COMPLETED`, `0:0` | `1` node, `00:59:44` charged elapsed; `0.995556` actual node-hours | Continued from inspected `s04` modal restart siblings with unchanged input/executable digests. | `accepted`: exact `t = 9.5`, zero strict counters, and complete terminal ranked products; enters but does not complete the production window. |
| `R16/s06_rankio_t9p5_t10` | `4744158`, `COMPLETED`, `0:0` | `1` node, `00:21:20` charged elapsed; `0.355556` actual node-hours | Continued from inspected `s05` modal restart siblings with unchanged input/executable digests. | `accepted`: exact `t = 10.0`, zero strict counters, and complete terminal ranked products; closes the production window. |
| `R02/s00_rankio_t0_t0p1` | `4744198`, `COMPLETED`, `0:0` | `1` node, `00:11:18` charged elapsed; `0.188333` actual node-hours | Fresh standard-layout timing pilot from `t = 0`; source bundle `athenak-feature-cgl-through-c5bd782b.bundle` SHA-256 `e25e86ec61f743fdfce76ec8f942cc3093620aa99b0576851211b3b7cc7a74c1`; immutable executable SHA-256 `df87684e9d2b7af33b36c2757779d15de87b84ef051ca9c4489f2efa358f5c48`. | `accepted`: exact `t = 0.1`, zero strict counters, complete initial/terminal ranked products, and authenticated modal restart siblings. |
| `R02/s01_rankio_t0p1_t0p25` | `4744205`, `COMPLETED`, `0:0` | `1` node, `00:17:06` charged elapsed; `0.285000` actual node-hours | Continued from inspected `s00` modal restart siblings with matching input and executable digests. | `accepted`: exact native snapshot boundary `t = 0.25`, zero strict counters, complete terminal ranked products, and measured developed continuation cost. |
| `R02/s03_rankio_t0p25_t1` | `4744249`, `COMPLETED`, `0:0` | `1` node, `01:27:10` charged elapsed; `1.452778` actual node-hours | Continued from inspected `s01` modal restart siblings with matching input and executable digests. The preceding `s02` reservation was released before submission after Slurm normal-QOS preflight rejected a `02:15:00` request above the 120-minute limit. | `accepted`: exact native restart boundary `t = 1.0`, zero strict counters, expected active hard-wall projection, complete terminal ranked products, and measured longer developed cost. |
| `R02/s04_rankio_t1_t1p5` | `4744518`, `COMPLETED`, `0:0` | `1` node, `01:03:08` charged elapsed; `1.052222` actual node-hours | Continued from inspected `s03` modal restart siblings with matching input and executable digests; source bundle `athenak-feature-cgl-through-5089701f.bundle` SHA-256 `c3e62cb50fd89c823b5d2610d13c024a9d293639b7f1404af01dd801c37066f0`. | `accepted`: exact native restart boundary `t = 1.5`, zero strict counters, expected active hard-wall projection, complete terminal ranked products, and measured late-time cost. |
| `R02/s05_rankio_t1p5_t2` | `4744913`, `COMPLETED`, `0:0` | `1` node, `01:04:08` charged elapsed; `1.068889` actual node-hours | Continued from inspected `s04` modal restart siblings with matching input and executable digests; source bundle `athenak-feature-cgl-through-0814d835.bundle` SHA-256 `a629b8da254cfc000dd4ceb75de57a40113ed7070c261c3ebddf907bb94ea4ad`. | `accepted`: exact native restart boundary `t = 2.0`, zero strict counters, expected active hard-wall projection, complete terminal ranked products, and measured late-time cost. |
| `R02/s06_rankio_t2_t2p5` | `4745305`, `COMPLETED`, `0:0` | `1` node, `01:03:43` charged elapsed; `1.061944` actual node-hours | Continued from inspected `s05` modal restart siblings with matching input and executable digests; source bundle `athenak-feature-cgl-through-9d9ee2f7.bundle` SHA-256 `c3053e7c312cfe67035cceae7b0ec407afa3cfc8b76ee27a27f71f07b01d2bac`. | `accepted`: exact native restart boundary `t = 2.5`, zero strict counters, expected active hard-wall projection, complete terminal ranked products, and measured late-time cost. |
| `R17/s00_rankio_t0_t0p05` | `4744210`, `COMPLETED`, `0:0` | `8` nodes, `00:15:50` charged elapsed; `2.111111` actual node-hours | Fresh high-resolution timing pilot from `t = 0`; source bundle `athenak-feature-cgl-through-cd63f81e.bundle` SHA-256 `361b3d45a0ecaca29cd84d17c231db965bab5c90fd18e652b6b5fedff4fa2b67`; immutable executable SHA-256 `df87684e9d2b7af33b36c2757779d15de87b84ef051ca9c4489f2efa358f5c48`. | `accepted`: exact `t = 0.05`, zero strict counters, complete 64-rank products, and balanced 27-meshblock-per-rank decomposition. |
| `R17/s01_rankio_t0p05_t0p1` | `4744230`, `COMPLETED`, `0:0` | `8` nodes, `00:15:56` charged elapsed; `2.124444` actual node-hours | Continued from inspected `s00` 64-rank modal restart siblings with matching input and executable digests. | `accepted`: exact `t = 0.10`, zero strict counters, complete terminal ranked products, and measured developed high-resolution cost. |

Completed lineage bundle
`runs/mks24-stage-i/E02-modal-driver/bundles/completed/R16_rankio_t10_20260530`
contains all seven accepted segments and 41 deduplicated ranked snapshots.
Its manifest SHA-256 is
`9955085f3557a15c9775722bd1574dbd3636adbbb8df4c4ba5e4e3b9b37d6aa8`.
The `t = 8`--`10` analyzer selects nine snapshots and 101 history rows,
passes its synthetic test, and closes applied forcing work to relative
residual `8.806e-13`; diagnostics JSON SHA-256 is
`b3547ff697391dc7d56588028da2d7dd595e9b1ad33fbf7887e0e98188f9e6ab`.
Completion-recost evidence JSON SHA-256
`d42a0432858127a94ab4aee5aab0c0b79572b70730085352aec53147670adad3`
records `6.145556` actual node-hours and a broad mapped-matrix runtime bracket
from `1087.763412` cell-count-only node-hours to `2956.012436` CFL-aware
node-hours. R16 is retained pipeline and cost evidence for its mapped Figure 11
role; the F-076 forcing-contract audit prevents using it as direct MKS24
reproduction evidence or as authorization for the unresolved standard-layout
matrix.

Standard-layout timing pilot
`runs/mks24-stage-i/E02-modal-driver/R02` contains two accepted segments
through native snapshot boundary `t = 0.25`. The authenticated
`t = 0.1`--`0.25` continuation uses `0.285000` node-hours, or `1.900000`
node-hours per simulated time unit, projecting `19.000000` node-hours for one
standard-layout `t = 10` case. R02 recost evidence JSON SHA-256
`81c9f6e56cca9eb5a69449a1166bf5d3e3cc8a0a4102b6b4fd0c3f73f3c2123b`
initially narrowed the mapped-matrix runtime bracket to
`424.145556`--`576.145556` node-hours. At the F-070 review, the remaining
spread was the `R17` high-resolution factor. R02 remains timing evidence only,
not an accepted paper-statistics bundle.

High-resolution timing pilot
`runs/mks24-stage-i/E02-modal-driver/R17` contains two accepted eight-node,
64-rank segments through authenticated continuation `t = 0.10`. Each rank
owns the intended standard-layout load of 27 meshblocks. The developed
`t = 0.05`--`0.10` continuation uses `2.124444` node-hours, or `42.488889`
node-hours per simulated time unit. Together with the reusable R02 and R17
prefixes, that projects `697.019444` node-hours for the full mapped E02
matrix. R17 recost evidence JSON SHA-256
`79147715729a2c9017400e170de2191a3f7d6b0bf411acee09a948ada464b822`
records the runtime, strict counters, product sizes, and F-071 `900.000000`
node-hour measured envelope.

First mapped-matrix continuation
`R02/s03_rankio_t0p25_t1` reaches exact native restart boundary `t = 1.0`
in `1.452778` node-hours with zero strict counters and expected active
hard-wall projection. The developed `t = 0.25`--`1.0` interval measures
`1.937037` node-hours per simulated time unit, updating the continuation-aware
matrix projection to `702.195370` node-hours while leaving `197.804630`
node-hours of F-071 envelope margin. F-072 evidence JSON SHA-256
`7ba577644cd65eefb4c9e3a782f284c36be9e86894af1ec17e6242eba1c17b06`
also records the unsubmitted `s02` preflight attempt and requires preparation
to reject normal-QOS requests above two hours before restart archival.

Second mapped-matrix continuation
`R02/s04_rankio_t1_t1p5` reaches exact native restart boundary `t = 1.5`
in `1.052222` node-hours with zero strict counters and expected active
hard-wall projection. The late-time `t = 1.0`--`1.5` interval measures
`2.104444` node-hours per simulated time unit. Applying that rate
conservatively to the unfinished standard-layout work updates the
continuation-aware matrix projection to `725.464938` node-hours while leaving
`174.535062` node-hours of F-071 envelope margin. F-073 evidence JSON SHA-256
`0319284b72628f2456087f378588515726f0fa8586b20660d37da2a20a1ba6ba`
records its inspection, terminal products, accounting, and recost.

Third mapped-matrix continuation
`R02/s05_rankio_t1p5_t2` reaches exact native restart boundary `t = 2.0`
in `1.068889` node-hours with zero strict counters and expected active
hard-wall projection. The late-time `t = 1.5`--`2.0` interval measures
`2.137778` node-hours per simulated time unit. Applying that rate
conservatively to the unfinished standard-layout work updates the
continuation-aware matrix projection to `730.081697` node-hours while leaving
`169.918303` node-hours of F-071 envelope margin. F-074 evidence JSON SHA-256
`d0c7e0f9cdd8031a9d88c195adee9822943a358eb55193b0e1146bf3b95a08e4`
records its inspection, terminal products, accounting, and recost.

Fourth mapped-matrix continuation
`R02/s06_rankio_t2_t2p5` reaches exact native restart boundary `t = 2.5`
in `1.061944` node-hours with zero strict counters and expected active
hard-wall projection. The late-time `t = 2.0`--`2.5` interval measures
`2.123888` node-hours per simulated time unit. Applying that rate
conservatively to the unfinished standard-layout work updates the
continuation-aware matrix projection to `728.164877` node-hours while leaving
`171.835123` node-hours of F-071 envelope margin. F-075 evidence JSON SHA-256
`f823ef6993045d5d90f016f7052518a353f348f7f83b9efdc19d32ee097f0c59`
records its inspection, terminal products, accounting, and recost.

### 8.3 Archived `E02-modal-driver` Entry Gate and Corrected-Epoch Reset

The archived `E02-modal-driver` campaign required all of the following before
its first production submission. Retain the evidence because it defines the
minimum qualification categories for the corrected epoch, but do not treat
the old build or any old gate result as authorization to prepare or submit a
new `E02` segment:

1. Focused local CPU and MPI tests for modal-driver restart identity,
   MKS24 forcing orientation/cadence/work accounting, CGL-LF AMR interaction,
   and F-056 layout-independent cap fractions.
2. One immutable HIP/MPI build whose source revision, Kokkos revision,
   toolchain, executable checksum, and build manifest are archived for Stage I.
   The retained replacement-driver build uses revision
   `462b9dbd53e085dea46c2478b567781576d7d03e` and executable SHA-256
   `df87684e9d2b7af33b36c2757779d15de87b84ef051ca9c4489f2efa358f5c48`;
   retain its Git bundle before submitting qualification.
3. Frontier debug evidence for rank-local restart identity, one-rank versus
   eight-rank GPU/MPI behavior, reduced active/passive/random paper smoke, and
   a corrected hard-wall active-beta-10 run through at least `t = 2`.
4. A standard-layout startup probe using the replacement executable and the
   intended rank-local snapshot/restart layout, so modal restart size and I/O
   behavior are measured rather than inferred from shared-file history.
5. The committed epoch-aware `cgl_lf_stage_i.py` path that stored E02 runs
   beneath `runs/mks24-stage-i/E02-modal-driver/`, keeps separate accounting
   files, refuses cross-epoch continuation or bundle discovery, preserves
   completed R16, R02, and R17 accounting, and limits current authorization to
   the frozen mapped `R02`--`R17` matrix under sequential inspection.
6. An all-user-job preflight confirming that no unrelated live CGL-root
   campaign will be modified or unintentionally overlapped.

Those historical gates passed for modal restart mechanics, I/O, and cost
measurement.
Fresh `E02` jobs `4743735` through `4744158` execute R16
from `t = 0` through exact `t = 10.0`, pass segment inspections with zero
strict counters, authenticate restart continuation, assemble the completed
bundle, and pass production-window analysis. R02 jobs `4744198` and `4744205`
complete the standard-layout timing pilot through native cadence `t = 0.25`.
R17 jobs `4744210` and `4744230` close high-resolution timing through
authenticated continuation `t = 0.10`. R02 job `4744249` then advances the
first mapped lineage through exact `t = 1.0`, and job `4744518` continues it
through exact `t = 1.5`. Job `4744913` continues it through exact `t = 2.0`.
Job `4745305` continues it through exact `t = 2.5`. F-076 then supersedes the
E02 continuation authorization: preserve these products as pipeline and cost
evidence. The F-078-qualified corrected executable and retained token now bind
fresh `R02/s00_rankio_t0_t0p1` job `4745922`, accepted from `t = 0` through
exact `t = 0.1`, and authenticated continuation `R02/s01_rankio_t0p1_t0p25`
job `4746154`, accepted through exact `t = 0.25`. Retained recost evidence
authorizes only `R02/s02_rankio_t0p25_t1` next. That continuation is now
accepted as job `4746182` through exact `t = 1.0`; updated retained recost
evidence historically authorized only `R02/s03_rankio_t1_t1p5` next. That
continuation is now accepted as job `4746356` through exact `t = 1.5`;
retained corrected recost evidence historically authorized only
`R02/s04_rankio_t1p5_t2` next. That continuation is now accepted as job
`4746435` through exact `t = 2.0`; retained corrected recost evidence
historically authorized only `R02/s05_rankio_t2_t2p5` next. That continuation
is now accepted as job `4746663` through exact `t = 2.5`; retained corrected
recost evidence historically authorized only `R02/s06_rankio_t2p5_t3` next.
That continuation is now accepted as job `4746773` through exact `t = 3.0`;
retained corrected recost evidence historically authorized only
`R02/s07_rankio_t3_t3p5` next. That continuation is now accepted as job
`4746953` through exact `t = 3.5`; retained corrected recost evidence
historically authorized only `R02/s08_rankio_t3p5_t4` next under its reviewed
one-segment `4500`-second threshold. That continuation is now accepted as job
`4747015` through exact `t = 4.0`; retained corrected recost evidence
historically authorized only `R02/s09_rankio_t4_t4p5` next with Slurm
walltime `01:40:00`, Athena timeout `01:30:00`, and a reviewed one-segment
`4800`-second threshold that could not ratchet automatically. That
continuation is now accepted as job `4747087` through exact `t = 4.5`;
retained corrected recost evidence historically authorized only
`R02/s10_rankio_t4p5_t5` next with Slurm walltime `01:40:00`, Athena timeout
`01:30:00`, and a reviewed one-segment `4800`-second threshold retaining
`600` seconds on both timeout margins. That continuation is now accepted as
job `4747146` through exact `t = 5.0`; retained corrected recost evidence
historically authorized only `R02/s11_rankio_t5_t5p5` next with Slurm walltime `01:40:00`,
Athena timeout `01:30:00`, and a reviewed one-segment `4800`-second threshold
retaining `600` seconds on both timeout margins. The threshold is scoped only
to `s11` and cannot ratchet automatically.
After clean `s11` inspection, accounting, reconciliation, and recost, propose quarter units if elapsed time exceeds `4800` seconds, remaining cap headroom falls below `300` seconds, or the reviewed trend estimator projects the next half-unit at or above `4800` seconds. Any scientific, provenance, scheduler, storage, budget, or reconciliation failure blocks successor preparation.
That continuation is now accepted as job `4747202` through exact `t = 5.5`;
retained corrected recost evidence historically authorized only
`R02/s12_rankio_t5p5_t5p75` on one node from authenticated `s11` terminal
siblings next, with Slurm walltime `01:05:00`, Athena timeout
`00:55:00`, and a reviewed one-segment `2700`-second threshold retaining
`600` seconds on both timeout margins. The threshold is scoped only to `s12`
and cannot ratchet automatically.
That quarter-unit continuation is now accepted as job `4747500` through exact
`t = 5.75`; retained corrected recost evidence historically authorized only
`R02/s13_rankio_t5p75_t6` on one node from authenticated `s12` terminal
siblings next, with Slurm walltime `01:05:00`, Athena timeout `00:55:00`, and a
reviewed one-segment `2700`-second threshold retaining `600` seconds on both
timeout margins. The threshold is scoped only to `s13` and cannot ratchet
automatically.
That quarter-unit continuation is now accepted as job `4747834` through exact
`t = 6.0`; retained corrected recost evidence historically authorized only
`R02/s14_rankio_t6_t6p25` on one node from authenticated `s13` terminal
siblings next, with Slurm walltime `01:05:00`, Athena timeout `00:55:00`, and a
reviewed one-segment `2700`-second threshold retaining `600` seconds on both
timeout margins. The threshold is scoped only to `s14` and cannot ratchet
automatically.
That quarter-unit continuation is now accepted as job `4748138` through exact
`t = 6.25`; retained corrected recost evidence historically authorized only
`R02/s15_rankio_t6p25_t6p5` on one node from authenticated `s14` terminal
siblings next, with Slurm walltime `01:05:00`, Athena timeout `00:55:00`, and a
reviewed one-segment `2700`-second threshold retaining `600` seconds on both
timeout margins. The threshold was scoped only to `s15` and could not ratchet
automatically.
That quarter-unit continuation is now accepted as job `4752118` through exact
`t = 6.5`; retained corrected recost evidence historically authorized only
`R02/s16_rankio_t6p5_t6p75` on one node from authenticated `s15` terminal
siblings next, with Slurm walltime `01:05:00`, Athena timeout `00:55:00`, and a
reviewed one-segment `2700`-second threshold retaining `600` seconds on both
timeout margins. The threshold was scoped only to `s16` and could not ratchet
automatically.
That quarter-unit continuation is now accepted as job `4753294` through exact
`t = 6.75`; promoted F-099 recost evidence historically authorized only
`R02/s17_rankio_t6p75_t7` on one node from authenticated `s16` terminal siblings
under the reviewed one-segment non-ratcheting profile.
That quarter-unit continuation is now accepted as job `4754008` through exact
`t = 7.0`; promoted F-100 recost evidence narrowly and conditionally authorizes
only `R02/s18_rankio_t7_t7p25` on one node from authenticated `s17` terminal
siblings next after these documentation bytes are committed, the resulting
actual post-F-100 source bundle and its SHA-256 are archived and cataloged,
hardened reconciliation is rerun, and the
queue/shared-root audit repeats the explicit stale beta-25 acknowledgement. Use
Slurm walltime `01:05:00`, Athena timeout `00:55:00`, and a reviewed one-segment
`2700`-second threshold retaining `600` seconds on both timeout margins. The
threshold is scoped only to `s18` and cannot ratchet automatically.

Progress on 2026-05-30: `g014` retained an input-parse failure because its
archived sizing deck lacked the `mhd/limiter_hardwall` key targeted by a
command-line override. The helper now rejects absent override keys during
`prepare`. Corrected `g015` passed the reduced one-node/eight-GPU startup gate
with explicit deck settings, two complete eight-rank binary groups, fifteen
complete eight-rank restart groups, and `28,403,848` bytes per rank-local
restart file. This does not satisfy the restart-identity, nonlinear, or
standard-layout requirements above.

Further progress on 2026-05-30: rank-local restart-load smoke `g016` passed.
The intended `g016` versus uninterrupted `g017` comparison was retained as
inconclusive because the exact `t = 0.01` stop clipped the integration
timestep. Natural-cycle comparison `g018` instead resumed `g017` checkpoint
`.00013.rst` and passed with maximum terminal MHD-history, user-history,
float32 snapshot-field, and turbulence-force-slice differences
`6.10e-17`, `1.82e-13`, `3.73e-9`, and `2.04e-14`. Histories/tabs use
`rtol = atol = 1e-12`; float32 snapshots use `rtol = 1e-6`,
`atol = 1e-7`. Future identity checks must branch from a natural-cycle
checkpoint unless both comparator legs intentionally share the same clipped
segment boundary.

One-rank comparator `g019` then followed the same 25-cycle uninterrupted
trajectory as eight-rank `g017` through `t = 0.02`. Maximum terminal
MHD-history, user-history, assembled float32 snapshot-field, and
turbulence-force-slice differences were `4.44e-17`, `7.11e-15`, `3.73e-9`,
and `8.88e-16`. This closes the reduced GPU/MPI decomposition gate; future
decomposition checks must compare assembled physical grids rather than
rank-local shard boundaries.

Reduced passive-Alfvénic smoke `g020` passed with perpendicular-only forcing
(`force_prl2 = 0`, `force_prp2 = 128`) and exactly zero passive
`lf_cpwrk`/`lf_cawrk`. Reduced active-random smoke `g021` retained
`driving_type = 0`, `rseed = 314159`, activated all three force components,
and produced finite nonzero pressure-work ledgers. Both reached `t = 0.01`
with zero strict counters and retained restart checkpoints.

Reduced nonlinear hard-wall `g022` then reached `t = 2.0` with active
algebraic projection, zero strict counters, finite retained ledgers, and
forcing-work relative residual `1.863e-12`. Standard-layout rank-local
startup and sizing run `g023` launched the guarded `192 x 192 x 384` layout,
reached `t = 0.01` cleanly, retained complete initial/final eight-rank
snapshot and modal-restart groups, and passed assembled endpoint reads. Its
startup-only products close the remaining Frontier debug gate without
establishing production runtime.

Focused post-merge local retention is complete: the serial CGL suites passed
`39` tests, and scheduled one-node MPI CPU job `4743666` passed both retained
decomposition regressions. With no queued user job after that allocation, the
six `E02-modal-driver` operational entry conditions above were satisfied. R16 is retained
through exact `t = 10.0`. R02 then closes the standard-layout timing gate
through native snapshot boundary `t = 0.25`. R17 closes the high-resolution
timing gate through authenticated continuation `t = 0.10`, holding the active
Alfvenic beta-10 physics fixed while changing from `192 x 192 x 384` to
`384 x 384 x 768`. R02 subsequently reaches exact `t = 1.0`. The next
inspected continuation reaches exact `t = 1.5` in job `4744518`. The next
inspected continuation reaches exact `t = 2.0` in job `4744913`. The next
inspected continuation reaches exact `t = 2.5` in job `4745305`. F-076
subsequently prohibits every new E02 preparation or submission because the
forcing contract is not exact. Do not continue this lineage.

### 8.4 Why Stage II Is Not Currently Reserved

The former full extension plan reserved `3735.627090` node-hours including
the earlier ledger. Before the `E02` reset, combining that extension matrix
with the archived Stage I planning reservation would have required:

```text
0.851670 + 1068.888889 + 3734.775420 = 4804.515979 node-hours,
```

which exceeds the `4000` node-hour ceiling by `804.515979` node-hours. The
fresh `E02` ceiling does not make that stale extension matrix executable.
Stage II remains unauthorized until accepted Stage I results and measured
replacement-driver costs exist. After Stage I, use the actual `E02` ledger to
design a reduced extension within the remaining ceiling or request a larger
allocation.

## 9. Stage I Manuscript Program

`docs/cgl_lf_validation.tex` presently documents method validation. Convert
it into the reproduction manuscript only after the production inputs,
reference-data decision, and execution authorization have been reviewed.

| Manuscript section | Required argument and evidence |
| --- | --- |
| Abstract | State the independent reproduction question, the exact CGL-LF scope achieved, principal quantitative outcomes once available, and any external/normalization limitation. |
| Introduction | Motivate magneto-immutability and the need to reproduce MKS24 before testing new guide-field/collisionality regimes. |
| MKS24 target contract | State the closure, box, forcing, beta/limiter/LF-scale cases, and what direct reproduction excludes unless externally supplied. |
| AthenaK method validation | Condense existing operator, forcing, limiter, restart, GPU, and work-accounting evidence needed to trust the independent reproduction. |
| Reproduction campaign | Present the sixteen mapped cases, any audited beta-1 disposition, output/analysis protocol, cost/storage provenance, and per-panel comparison method. |
| Reproduced equation of state and compressive behavior | Figures 1-4 comparisons, with blocked status if required normalization remains unresolved. |
| Reproduced cascade and self-organization | Figures 5-9 and 11 comparisons: spectra, transfer, alignment, structure, and resolution. |
| Reproduced closure sensitivities | Figures 12-13: heat-flux scale and limiter-induced scattering rate. |
| External comparison boundary | Treat Figure 10 and any unregenerated hybrid-kinetic content explicitly. |
| Discussion and conclusions | State which MKS24 CGL-LF results passed, failed, or remain blocked; introduce Stage II only as subsequent work unless it has later been executed. |
| Reproducibility appendix/data statement | Inputs, hashes, source checksums, extraction/analysis configuration, comparison table, ledger, and retained-data locations. |

The main reproduction visuals should follow the published panel logic rather
than replacing it with a new narrative. Additional numerical-control and
provenance figures may be placed in appendices or supplemental material.

### 9.1 Manuscript acceptance checklist

A manuscript may claim completion of the Stage I reproduction only when:

1. It follows `docs/writing_style_guide.md` and makes no perfection claim.
2. Every Stage I run used for a claim is accepted, archived, and linked to its
   case role, analysis interval, provenance, and budget entry.
3. Every MKS24 simulation figure/result is assigned a status from Section 6.4.
4. All status `passed` statements have quantitative comparison data and
   uncertainty; all blocked or external panels are disclosed.
5. A full-paper reproduction claim is made only if the Figure 10 external
   model gate and any other external-result gates have been resolved.
6. TeX, bibliography, compact plotted source data, rendered figures, and
   provenance are tracked; large production data are retained externally with
   checksums.
7. The PDF builds through the declared TeX/BibTeX workflow without missing
   figures or references.

## 10. Stage II: Later Guide-Field and Background-Collisionality Extension

Stage II preserves the new experiment requested before the MKS24 priority was
identified. It is not part of Stage I reproduction and must not begin until
the Stage I acceptance decision is recorded.

### 10.1 New questions

After reproducing MKS24, ask:

1. How do magneto-immutability, anisotropic pressure work, and structural
   diagnostics change when field wandering is much larger than in MKS24?
2. How does a uniform background collisional relaxation rate `nu_coll`
   change the turbulent state, separately from MKS24's threshold-activated
   anomalous limiter scattering?
3. Are guide-field-dependent responses robust across beta and LF closure
   scale?

### 10.2 Extension anchors

Use a periodic cubic box and isotropic solenoidal driving only as an
extension design, not a reproduction configuration. Define:

| Extension guide regime | Collisionless anchor target |
| --- | --- |
| Weak guide `wg` | Late-time `R_B = B_rms/B_mean = 2.0 +/- 0.2`, `M_s = 0.50 +/- 0.05`, and primary `beta_rms approximately 10`. |
| Stronger guide `sg` | Late-time `R_deltaB = delta_B_rms/B_mean = 0.50 +/- 0.05`, `M_s = 0.50 +/- 0.05`, and primary `beta_rms approximately 10`. |

The retained target definitions are:

```math
\boldsymbol{B}_{\rm mean}=\langle\boldsymbol{B}\rangle_V,\qquad
B_{\rm mean}=|\boldsymbol{B}_{\rm mean}|,\qquad
B_{\rm rms}=\langle|\boldsymbol{B}|^2\rangle_V^{1/2},
```

```math
\delta B_{\rm rms}=
\langle|\boldsymbol{B}-\boldsymbol{B}_{\rm mean}|^2\rangle_V^{1/2},
\qquad
p_{\rm iso}={2p_\perp+p_\parallel\over3},
```

```math
M_s={u_{\rm rms}\over
\sqrt{(5/3)\langle p_{\rm iso}\rangle_V/\langle\rho\rangle_V}},
\qquad
\beta_{\rm rms}={2\langle p_{\rm iso}\rangle_V\over B_{\rm rms}^2}.
```

The weak-guide definition deliberately uses total `B_rms`; the stronger-guide
definition deliberately uses fluctuating `delta_B_rms`. Report both magnetic
ratios for both anchors and do not interchange these targets during tuning.

Calibrate anchors with a calibration-only seed and hold science seeds out of
tuning. Measure each collisionless anchor turnover time
`tau_eddy,0,G` and define the background-collision scan:

```math
C_G = \nu_{\rm coll}\tau_{{\rm eddy},0,G}
    \in \{0,\ 0.1,\ 1,\ 10,\ 100\}.
```

For nonzero `C_G`, keep each accepted anchor's imposed field, forcing,
thermodynamics, LF scale, and seed fixed while changing only `nu_coll` and
case metadata. Drift in magnetic amplitude, Mach number, or beta is then a
physical response rather than a calibration failure. Separately tuned
state-matched runs may be secondary comparisons only.

### 10.3 Candidate Stage II matrix

This is the retained scientific design to be recosted after Stage I, not an
execution authorization. A post-reproduction approval may reduce it only by
recording which resulting inference is relinquished.

| Tier | Guide regimes and `C_G` | Resolutions | Models and seeds | Purpose |
| --- | --- | --- | --- | --- |
| Calibration | `wg`, `sg`; `0` | `32^3`, then `64^3` | Active finalists and matched passive confirmation; calibration seed `161803` only. | Freeze anchor inputs and `tau_eddy,0,G`; never enter scientific averages. |
| Collisionless anchor ladder | `wg`, `sg`; `0` | `64^3`, `96^3`, `192^3`, `384^3` | Active/passive; seed `271828`. | Resolve guide-field and feedback baselines. |
| Collisionless seed repeat | `wg`, `sg`; `0` | `96^3`, `192^3` | Active/passive; seed `314159`. | Estimate realization uncertainty. |
| Background-collision sweep | `wg`, `sg`; `0.1`, `1`, `10`, `100` | `96^3`, `192^3` | Active/passive; seed `271828`. | Measure fixed-input `nu_coll` response. |
| Strongly collisional endpoint | `wg`, `sg`; `100` | `384^3` | Active/passive; seed `271828`. | Test resolution of the largest collision contrast. |
| Collisional seed repeat | `wg`, `sg`; `1`, `100` | `192^3` | Active/passive; seed `314159`. | Check sampling at the transition and endpoint. |
| Selected sensitivities | Both guide regimes; collisionless and `C_G = 1` roles selected before execution. | `96^3`, with `192^3` beta checks only if budgeted. | Active/passive beta checks at `beta_rms approximately 1,100`; active-only LF scale factors `1/2,2`; seed `271828`. | Bound beta and closure-scale dependence of principal inferences. |

Stable production identities should encode guide regime, collision coordinate,
feedback choice, resolution, and seed, for example
`g_wg_c1_active_N192_s271828`. For every active/passive pair, reuse the
identical archived forcing realization. For every fixed-input collision
series, change only `nu_coll` and identifying metadata after its collisionless
anchor has been frozen.

### 10.4 Extension safeguards and acceptance contract

| Safeguard or secondary check | Requirement |
| --- | --- |
| Held-out calibration | Use tuning-only seed `161803`; science seeds `271828` and `314159` are not opened before anchor freeze. |
| Low-field validity | Record `min(|B|/B_mean)`, `f_low(q)` for `q = 0.01, 0.05, 0.10`, and field-direction regularization counters; no unqualified claim if local-field validity fails. |
| Quantitative structure | Use isotropic shell spectra as primary weak-guide scale-space measures and local-field structure functions for anisotropy; guide-relative spectra are labeled secondary. |
| Visual selection | Pre-register midpoint snapshots, central planes, and pooled percentile color limits before reviewing visuals. |
| Beta sensitivity | After primary `beta_rms approximately 10`, assess selected cases at `beta_rms approximately 1` and `100`. |
| LF-scale sensitivity | At minimum vary `lf_k_parallel` by factors of two around the baseline in selected anchor/transitional cases. |
| Estimands | Separate active/passive feedback effects, fixed-input background-collisionality effects, and weak/strong-guide contrasts. |
| Late-time sampling | Require a stationary accepted window spanning at least four measured turnover times, with block/bootstrap uncertainty retained for reported statistics. |
| Safety and accounting | Require zero invalid-state/undefined-direction failures and retain forcing, applied pressure/LF work, limiter, and energy-residual evidence. |
| Resolution claims | Describe a guide/collision effect as resolution-robust only after its sign survives the declared seed checks and the applicable `384^3` endpoint comparison. |

The minimum Stage II analysis bundle must contain target/stationarity
decisions, low-field and limiter exposure, active/passive and collision
contrasts with uncertainty, isotropic spectra and local-field structure
products, selected common-scale slices, exact submitted-input/build hashes,
and compute/storage accounting. No Stage II result enters the manuscript
unless that bundle records its accepted case and comparison role.

### 10.5 Stage II recosting rule

Do not reuse the former `3735.627090` node-hour extension reservation after
Stage I is inserted ahead of it. After Stage I:

1. Use measured MKS24 performance and remaining ledger rather than the prior
   speculative scaling.
2. Decide whether a reduced extension can answer the main new question within
   the measured remainder of the corrected fresh epoch.
3. If the full two-guide, five-collisionality, beta/LF-sensitivity, paired
   `384^3` program remains scientifically necessary, request an amended total
   budget before submission.

## 11. Execution Order and Stop Conditions

| Phase | Work | Deliverable | Stop condition |
| --- | --- | --- | --- |
| A. Reproduction specification audit | Verify pinned MKS24 source, sixteen-run mapped alias map, disposition of the unmapped active-Alfvenic beta-1 definition, closure/forcing/limiter identity, and figure/status table. | Reviewed Stage I protocol. | A published simulation role or observable remains unidentified. |
| B. Reference-data closure | Resolve dimensional panel normalization/data boundary and the Figure 10 external-model decision. | Qualified reference manifests or explicit scoped limitation. | A "complete" claim would depend on unqualified data. |
| C. Production readiness | Retain and commit the bundle-backed source-provenance follow-up; validate required local/GPU/restart/work-accounting tests for the replacement driver; complete debug qualification of the corrected immutable build. | Passing readiness review, canonical corrected-epoch qualification token, and fail-closed corrected-epoch ledger. | Numerical gates, retained provenance, epoch isolation, qualification approval, or submission controls fail. |
| D. MKS24 production | Preserve completed E02 R16 and R02/R17 timing products only as pipeline and cost evidence. The reviewed E03 token and clean reconciliation now bind fresh R02 jobs `4745922`, `4746154`, `4746182`, `4746356`, `4746435`, `4746663`, `4746773`, `4746953`, `4747015`, `4747087`, `4747146`, `4747202`, `4747500`, `4747834`, `4748138`, `4752118`, `4753294`, and `4754008`, accepted from `t = 0` through exact `t = 7.0` for `16.133334` cumulative node-hours. F-079--F-100 controller provenance, hardening, and accepted-segment recosting are archived or promoted. Before preparing a successor, commit these documentation bytes, archive and catalog the resulting actual post-F-100 source bundle and its SHA-256, rerun hardened reconciliation, and repeat the queue/shared-root audit with the explicit stale beta-25 acknowledgement. Then prepare only `R02/s18_rankio_t7_t7p25` on one node from the authenticated `s17` terminal siblings with Slurm walltime `01:05:00`, Athena timeout `00:55:00`, and its reviewed one-segment `2700`-second threshold retaining `600` seconds on both timeout margins, inspect, account, reconcile, and recost, then finish R02 and execute R03--R16 sequentially. The threshold is scoped only to `s18` and cannot ratchet automatically. Any scientific, provenance, scheduler, storage, budget, or reconciliation failure blocks successor preparation. Execute R17 last. The beta-1 inventory definition remains excluded unless Phase A assigns it a published role. | Accepted corrected-epoch run bundles and measured ledger/storage. | Forcing fidelity, safety, stationarity, accounting, storage, or budget gate fails. |
| E. MKS24 analysis | Generate each required product and complete the figure-by-figure status table. | Quantitative reproduction report. | Any claimed result lacks qualified comparison evidence. |
| F. Reproduction manuscript | Transform the TeX note into a buildable Stage I manuscript. | Manuscript and reproducibility package. | Claims exceed the accepted status table. |
| G. Extension decision | Recalculate the Stage II design using measured Stage I results/costs. | Approved extension matrix or deferred-work record. | Insufficient budget or unresolved baseline reproduction. |
| H. Extension execution/manuscript integration | Execute only an approved Stage II matrix and extend the manuscript if scientifically supported. | New-result supplement or expanded manuscript. | Any extension gate fails or funding is absent. |

## 12. Immediate Handoff

As of the accepted job `4754008` record, corrected F-100 recost, and inspected
timestamp `2026-06-02T03:10:59+00:00`:

1. The repository contains the TeX validation note, the writing style guide,
   a detailed MKS24 implementation/evidence plan, substantial MKS24 analysis
   infrastructure, guarded inputs for sixteen source-mapped CGL-LF production
   roles plus one unresolved active-Alfvenic beta-1 definition, a committed
   sixteen-run mapped execution manifest, a production-QOS accounting utility
   with terminal-time inspection and accepted-bundle assembly gates, and
   reduced/GPU qualification evidence. The stop-line protocol review is
   recorded in `docs/cgl_lf_mks24_stage_i_protocol_review.md`; it is not a
   corrected-build scientific signoff.
2. Frontier job `4674731` completed as the legacy shared-output
   `R16/s00_t0_t2` diagnostic segment and was recorded `clean_partial` at
   `t = 1.66204848945`. The ranked-output replacement is inspected and
   accepted through its `t = 6.5` prefix via jobs `4675322`, `4676557`,
   `4676696`, `4679803`, `4681476`, and `4683291`;
   its lineage-only diagnostic prefix bundle passed a separately labeled
   transient-window analyzer pipeline check without promoting that prefix to
   a reproduced result. Segment `R16/s05_rankio_t5_t6p5`
   terminated cleanly on application walltime at `t = 6.43828031751 < 6.5`,
   was recorded `clean_partial` with complete terminal ranked products and
   zero strict LF safety counters, and used `1.753333` node-hours. Short
   continuation `R16/s06_rankio_t6p438280_t6p5` reached exact `t = 6.5`,
   was recorded `accepted` after formal inspection, and used `0.080556`
   node-hours. The accepted-prefix `t = 6.5` bundle contains twenty-eight
   deduplicated ranked snapshots but remains outside the production
   comparison window. Continuation `R16/s07_rankio_t6p5_t7p5` reached
   exact `t = 7.5`, was recorded `accepted`, and used `1.281944`
   node-hours; its generic prefix analysis still selects zero samples in the
   `t = 8`--`10` window. Continuation `R16/s08_rankio_t7p5_t8p5` reached
   exact `t = 8.5` as job `4686032`, passed formal clean inspection, and
   used `1.282778` node-hours, but is recorded `rejected` for current
   reproduction admission because the replacement turbulent driver changes
   forcing/restart state.
   That archived epoch contains no complete admitted `t = 10` `R16` run bundle
   or reproduced panel. Fresh `E02-modal-driver` R16 completion is recorded
   separately below.
3. It does **not** contain completed paper-scale MKS24 production
   simulations, completed quantitative comparisons for the full figure
   program, resolved dimensional reference mappings for all panels,
   independent reproduction of the Figure 10 hybrid-kinetic comparator, or a
   completed paper manuscript.
4. The shared Frontier root also contains exploratory, non-MKS24 comparison
   campaigns beneath `runs/beta25-accel05-*`. Job `4743020` was cancelled and
   cleaned while preserving provenance. Two-node debug job `4743106` ran from
   revision `89835fea` and Slurm now records it as `COMPLETED 0:0`; its
   top-level campaign manifest remains stale at `state = running`. It is not
   Stage I evidence and must be accounted separately. Do not modify or prune
   that tree; explicitly acknowledge the reviewed stale record during
   preflight.
5. The next task is not to launch the weak-guide extension or continue old
   restarts. The isolated `E02-modal-driver` accounting namespace,
   shared-root helper checks, immutable HIP/MPI executable, and retained
   source bundles are archived. Reduced nonlinear hard-wall `g022`,
   standard-layout rank-local sizing `g023`, focused serial CGL regression,
   and scheduled MPI CPU regression gates passed. Fresh `R16` jobs `4743735`
   through `4744158` reach exact `t = 10.0`; the completed
   bundle passes its `t = 8`--`10` analysis. R02 jobs `4744198` and `4744205`
   complete the standard-layout timing pilot through native cadence `t = 0.25`
   in `0.473333` node-hours. R17 jobs `4744210` and `4744230` complete the
   high-resolution timing pilot through authenticated continuation `t = 0.10`
   in `4.235556` node-hours. F-071 historically authorized only the frozen
   mapped matrix under a `900.000000` node-hour envelope and sequential
   inspection. R02 job
   `4744249` advances the first mapped lineage through exact `t = 1.0` in
   `1.452778` node-hours; F-072 updates the matrix projection to `702.195370`
   node-hours and hardens the two-hour normal-QOS preparation limit. R02 job
   `4744518` continues the same inspected lineage through exact `t = 1.5` in
   `1.052222` node-hours; F-073 updates the projection to `725.464938`
   node-hours. R02 job `4744913` continues through exact `t = 2.0` in
   `1.068889` node-hours; F-074 updates the projection to `730.081697`
   node-hours. R02 job `4745305` continues through exact `t = 2.5` in
   `1.061944` node-hours; F-075 updates the projection to `728.164877`
   node-hours. Continuation `R02/s07_rankio_t2p5_t3` was prepared from the
   inspected `s06` terminal restart siblings and submitted as job `4745498`
   on 2026-05-30. The independent Phase A source-to-code forcing audit then
   invalidated further direct-reproduction use of E02. Job `4745498` was
   cancelled after `498` seconds and recorded `aborted`, consuming `0.138333`
   node-hours. E02 cumulative use is `15.628610` node-hours with no active
   reservation. Preserve E02 as pipeline and cost evidence. Corrected E03
   Frontier qualification jobs `4745621`, `4745643`, `4745651`, `4745756`,
   `4745825`, `4745843`, `4745848`, and `4745890` now pass as `g024` through
   `g031` against immutable executable SHA-256
   `68f243f9204df388b24365ae65a567f6f567dbe422a6d7a43b9fb4a499ef118c`.
   The reviewed E03 approval token was retained and the initially empty
   canonical E03 ledger reconciled cleanly. Fresh
   `R02/s00_rankio_t0_t0p1` job `4745922` was submitted from `t = 0`,
   formally inspected, and recorded `accepted` at exact `t = 0.1` for
   `0.188889` node-hours. F-079/F-080 controller provenance is committed,
   archived, and reconciled. Authenticated `R02/s01_rankio_t0p1_t0p25` job
   `4746154` is also accepted through exact `t = 0.25` for `0.289167`
   node-hours. Authenticated `R02/s02_rankio_t0p25_t1` job `4746182` is
   accepted through exact `t = 1.0` for `1.465556` node-hours, bringing
   corrected E03 use to `1.943612`. Authenticated
   `R02/s03_rankio_t1_t1p5` job `4746356` is accepted through exact
   `t = 1.5` for `1.048611` node-hours, bringing corrected E03 use to
   `2.992223`. Retained corrected recost evidence JSON SHA-256 is
   `b6fd6e8fcd939ca1baa06411a3cd6f04359b3bc2ff5b2a0097e2cc87b1763fc7`;
   its authorization-facing projection is `724.645556` node-hours with
   `175.354444` margin. F-083/F-084 controller hardening is promoted at
   canonical revision `9689c269bf329542815a1b2b137881126964b05c`, helper SHA-256
   `1c633ebb58294938a0a0609742ae8f8d1f88cd15242ffe649578796edeb39375`;
   hardened reconciliation is clean. Authenticated `R02/s04_rankio_t1p5_t2`
   job `4746435` is accepted through exact `t = 2.0` for `1.063611`
   node-hours, bringing corrected E03 use to `4.055834`. Retained corrected
   recost evidence JSON SHA-256 is
   `c382b78ab9e21466648adf9d7ea38d2b407f0e986579f128bdfe6bc44bef79dd`;
   its authorization-facing projection is `728.845556` node-hours with
   `171.154444` margin. Authenticated `R02/s05_rankio_t2_t2p5` job `4746663`
   is accepted through exact `t = 2.5` for `1.062500` node-hours, bringing
   corrected E03 use to `5.118334`. Retained corrected recost evidence JSON
   SHA-256 is
   `39115585cf7ca866e3aade1e53c69aab6ee9217361d75db7946748905e1c0e5c`;
   its authorization-facing projection is `728.534445` node-hours with
   `171.465555` margin. Authenticated `R02/s06_rankio_t2p5_t3` job `4746773`
   is accepted through exact `t = 3.0` for `1.077500` node-hours, bringing
   corrected E03 use to `6.195834`. Retained corrected recost evidence JSON
   SHA-256 is
   `d2c5e27553426beb03c8c1882bd6473b682029c6ba67acd32bb06c5fb9a3ff01`;
   its authorization-facing projection is `732.734445` node-hours with
   `167.265555` margin. Authenticated `R02/s07_rankio_t3_t3p5` job `4746953`
   is accepted through exact `t = 3.5` for `1.151111` node-hours, bringing
   corrected E03 use to `7.346945`. Retained corrected recost evidence JSON
   SHA-256 is
   `468e787db48d98661cd3ffe125829a78f135bb7ca4c9104c766a00b693f586a3`;
   its authorization-facing projection is `753.345556` node-hours with
   `146.654444` margin. The accepted clean interval took `4144` seconds,
   above the prior `4000`-second operational guard. Authenticated
   `R02/s08_rankio_t3p5_t4` job `4747015` is accepted through exact `t = 4.0`
   for `1.200000` node-hours, bringing corrected E03 use to `8.546945`.
   Retained corrected recost evidence JSON SHA-256 is
   `dd5b0f5b82cf81005c1a481ee83d1127aa43861dcdd8c22765170e28b0ad6c99`;
   its authorization-facing projection is `767.034445` node-hours with
   `132.965555` margin. The accepted clean interval took `4320` seconds,
   satisfying the reviewed `s08` cap. Hardened reconciliation closes with
   `9/9/9` ledger rows/manifests/reservations, no active reservation, and no
   transaction. Authenticated `R02/s09_rankio_t4_t4p5` job `4747087` is
   accepted through exact `t = 4.5` for `1.208056` node-hours, bringing
   corrected E03 use to `9.755001`. Retained corrected recost evidence JSON
   SHA-256 is
   `8b10b1786db520740b687d989207f3cda6e0123f9bf2b62e75a5df3849511561`;
   its authorization-facing projection is `769.290001` node-hours with
   `130.710000` margin. The accepted clean interval took `4349` seconds,
   satisfying the reviewed `s09` cap. Hardened reconciliation closes with
   `10/10/10` ledger rows/manifests/reservations, no active reservation, and
   no transaction. Authenticated `R02/s10_rankio_t4p5_t5` job `4747146` is
   accepted through exact `t = 5.0` for `1.226389` node-hours, bringing
   corrected E03 use to `10.981390`. Retained corrected recost evidence JSON
   SHA-256 is
   `5256fb7fd9b75f175a16c2c16ee8a995b8d6d9281a92b70a979ea62171835bf2`;
   its authorization-facing projection is `774.423334` node-hours with
   `125.576666` margin. The accepted clean interval took `4415` seconds,
   satisfying the reviewed `s10` cap. Hardened reconciliation closes with
   `11/11/11` ledger rows/manifests/reservations, no active reservation, and
   no transaction. Authenticated `R02/s11_rankio_t5_t5p5` job `4747202` is
   accepted through exact `t = 5.5` for `1.268889` node-hours, bringing
   corrected E03 use to `12.250279`. Retained corrected recost evidence JSON
   SHA-256 is
   `323d7cb3cbe108eee944045c189cff75ffc07d06833c49bcfa8315f07fad8d51`;
   its authorization-facing projection is `786.323334` node-hours with
   `113.676666` margin. The accepted clean interval took `4568` seconds,
   leaving `232` seconds below the reviewed `s11` cap. Hardened reconciliation
   closes with `12/12/12` ledger rows/manifests/reservations, no active
   reservation, and no transaction. Authenticated
   `R02/s12_rankio_t5p5_t5p75` job `4747500` is accepted through exact
   `t = 5.75` for `0.660000` node-hours, bringing corrected E03 use to
   `12.910279`. Retained corrected recost evidence JSON SHA-256 is
   `d0a60e222138971e1bc9aae978bb3d62aee460f09000ee62f7f9ffdfe534165f`;
   its authorization-facing projection is `814.945556` node-hours with
   `85.054444` margin. The accepted clean interval took `2376` seconds,
   leaving `324` seconds below the reviewed `s12` threshold. Hardened
   reconciliation closes with `13/13/13` ledger rows/manifests/reservations, no
   active reservation, and no transaction. Authenticated
   `R02/s13_rankio_t5p75_t6` job `4747834` is accepted through exact
   `t = 6.0` for `0.599722` node-hours, bringing corrected E03 use to
   `13.510001`. Retained corrected recost evidence JSON SHA-256 is
   `4bc19f5b44587e73869982425b1cf0ad459547c3e0c5de1ab1b2b5cb366c6322`;
   its authorization-facing projection is `814.945556` node-hours with
   `85.054444` margin. The accepted clean interval took `2159` seconds,
   leaving `541` seconds below the reviewed `s13` threshold. Hardened
   reconciliation closes with `14/14/14` ledger rows/manifests/reservations, no
   active reservation, and no transaction. Authenticated
   `R02/s14_rankio_t6_t6p25` job `4748138` is accepted through exact
   `t = 6.25` for `0.655278` node-hours, bringing corrected E03 use to
   `14.165279`. Retained corrected recost evidence JSON SHA-256 is
   `5ba7aec52d3e7824cf4ccf52f95ce1e46368c84ed59f0525e2bbbcfc167fde91`;
   its authorization-facing projection is `829.101112` node-hours with
   `70.898888` margin. The accepted clean interval took `2359` seconds,
   leaving `341` seconds below the reviewed `s14` threshold. Hardened
   reconciliation closes with `15/15/15` ledger rows/manifests/reservations, no
   active reservation, and no transaction. Authenticated
   `R02/s15_rankio_t6p25_t6p5` job `4752118` is accepted through exact
   `t = 6.5` for `0.653611` node-hours, bringing corrected E03 use to
   `14.818890`. Retained corrected recost evidence JSON SHA-256 is
   `3603077400bd4530001155b83bfc4ac1053ec18505c56a6c29ffadc795bd1255`;
   its authorization-facing projection is `829.101112` node-hours with
   `70.898888` margin. The accepted clean interval took `2353` seconds,
   leaving `347` seconds below the reviewed `s15` threshold. Hardened
   reconciliation closes with `16/16/16` ledger rows/manifests/reservations, no
   active reservation, and no transaction. Authenticated
   `R02/s16_rankio_t6p5_t6p75` job `4753294` is accepted through exact
   `t = 6.75` for `0.647222` node-hours, bringing corrected E03 displayed use
   to `15.466112` (`55678 / 3600 = 15.466111111111111` exact node-hours).
   Retained corrected recost evidence JSON SHA-256 is
   `5b997623a3f8d83c200034f42e0c9f1b9a09c811bcb46e2fc1a0f7c1eb58815e`;
   its authorization-facing projection is `829.101112` node-hours with
   `70.898888` margin. The accepted clean interval took `2330` seconds,
   leaving `370` seconds below the reviewed `s16` threshold. Its
   sampled-history forcing-work relative residual is
   `6.437727059090897e-12`; accepted-prefix residual is
   `3.1385227986388927e-12`; terminal `lf_hwproj = 259077916291`.
   Hardened reconciliation closes with `17/17/17` ledger
   rows/manifests/reservations, no active reservation, no transaction, and no
   issue. That historical authorization was only for
   `R02/s17_rankio_t6p75_t7` under the full F-099 profile. Authenticated
   `R02/s17_rankio_t6p75_t7` job `4754008` is accepted through exact `t = 7.0`
   for `0.667222` displayed node-hours, bringing
   corrected E03 use to `16.133334` displayed node-hours
   (`58080 / 3600 = 16.133333333333333`
   exact node-hours). Retained corrected F-100 recost evidence JSON SHA-256 is
   `63438d8843b86815a2e25d4d6de6ef78ce9756c561fe9db3e762965ada73202f`; its authorization-facing projection is
   `829.101112` node-hours with
   `70.898888` margin. The accepted clean
   interval took `2402` seconds, leaving
   `298` seconds below the reviewed `s17`
   threshold. Its sampled-history forcing-work relative residual is
   `7.10881221809853e-12`; accepted-prefix residual is
   `3.2704293797689304e-12`; terminal `lf_hwproj =
   265802141460`. Hardened reconciliation closes with
   `18/18/18 ledger rows/manifests/reservations, no active reservation, no transaction, and issues = []`. Before preparing a successor, commit
   these documentation bytes, archive and catalog
   the resulting actual post-F-100 source bundle
   with SHA-256, rerun hardened
   reconciliation, and repeat the queue/shared-root audit with the explicit
   stale beta-25 acknowledgement. Then prepare only
   `R02/s18_rankio_t7_t7p25` on one node from the authenticated `s17` terminal
   siblings with Slurm walltime `01:05:00`, Athena timeout `00:55:00`, and its
   reviewed one-segment `2700`-second threshold retaining `600` seconds on both
   timeout margins. This threshold is scoped only to `s18` and cannot ratchet
   automatically. Any scientific, provenance, scheduler, storage, budget, or
   reconciliation failure blocks successor preparation. Inspect, account,
   reconcile, and recost. Finish R02, execute R03--R16 sequentially, and
   execute R17 last.

## 13. Reproducibility Record

For every accepted Stage I or later Stage II calculation, archive:

1. Git revision, executable/build/platform metadata, exact submitted input,
   workflow role, model/forcing/seed metadata, parent restart, and case alias
   mapping.
2. Model controls: `beta0`, domain/resolution, LF scale, thresholds,
   `nu_coll`, limiter scattering policy/rate, forcing shell/injection/time
   correlation, and output cadence.
3. Consumed/reserved node-hours, walltime, storage footprint, retention and
   pruning decisions, and allocation identifier.
4. Analysis intervals, raw histories, required snapshots, products, scripts,
   comparison reference checksums, normalization rules, uncertainties, figure
   files, and pass/fail/blocked decisions.
5. All failures and excluded results; an unsuccessful case may not be silently
   replaced by a selected successful realization.

Before drafting a result statement, classify its evidence:

| Verb | Use only when |
| --- | --- |
| `reproduces` | A matched published result passes the declared quantitative comparison and provenance gates. |
| `shows` | A result is directly measured in accepted simulations and survives its stated uncertainty/gates. |
| `suggests` | An inference is supported but limited by model, sampling, reference-data, or resolution bounds. |
| `does not establish` | The desired claim, including perfection or complete reproduction across an unresolved external-model boundary, exceeds evidence. |

## 14. References and Requirement Traceability

### 14.1 References and local sources

1. Majeski, S., Kunz, M. W., and Squire, J. (2024),
   *Self-organization in collisionless, high-beta turbulence*,
   Journal of Plasma Physics 90, 535900601, arXiv:2405.02418.
2. Squire, J., Schekochihin, A. A., Quataert, E., and Kunz, M. W. (2019),
   *Magneto-immutable turbulence in weakly collisional plasmas*,
   Journal of Plasma Physics, arXiv:1811.12421.
3. Squire, J., Kunz, M. W., Arzamasskiy, L., Johnston, Z., Quataert, E.,
   and Schekochihin, A. A. (2023), *Pressure anisotropy and viscous heating
   in weakly collisional plasma turbulence*, Journal of Plasma Physics,
   arXiv:2303.00468.
4. `docs/cgl_lf_mks24_reproduction_implementation_plan.md`: implementation,
   reference-extraction, qualification, and current evidence record.
5. `docs/writing_style_guide.md`: manuscript structure and claim-calibration
   requirements.
6. `docs/cgl_lf_validation.tex`: current manuscript source to convert.

### 14.2 Traceability

| Requirement | Where specified | Evidence required before completion |
| --- | --- | --- |
| Reproduce MKS24 before extension | Sections 1-2 and 11 | Stage I status table and accepted manuscript claim boundary. |
| Match MKS24 CGL-LF physical/numerical setup | Section 3 | Audited manifests, accepted production cases, and provenance. |
| Include MKS24 collisionality result correctly | Sections 3.4, 4, and 5 | Beta-100 `nu_lim` comparisons for `20`, `200`, and hard wall; no conflation with `nu_coll`. |
| Cover all MKS24 numerical figures/results and theory appendices | Section 5 | Per-panel pass/fail/blocked/external-model table and checked derivation record. |
| Address non-CGL Figure 10 honestly | Sections 2.1 and 5 | External dataset/kinetic rerun or explicitly limited reproduction claim. |
| Follow `writing_style_guide.md` in a complete manuscript | Section 9 | Converted TeX, figure/provenance package, and successful PDF build. |
| Use up to `4000` node-hours responsibly | Section 8 | Fail-closed ledger, measured updates, and Stage I-first reservation. |
| Retain requested weak/strong-guide and uniform-collisionality study | Section 10 | Post-reproduction, recosted and approved Stage II protocol. |
