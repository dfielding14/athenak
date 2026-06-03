# MKS24 Stage I Protocol Review: Stop-The-Line Audit

## Decision

Read `docs/cgl_lf_phase_i_handoff.md` first for the durable current production
boundary. Accepted corrected-E03 `R02/s19_rankio_t7p25_t7p5` job `4757300`
reached exact `t = 7.5` in `2450` seconds. Formal and independently
regenerated inspection pass with `18` retained products, zero strict LF
failure counters, no restart-marker bypass, terminal
`lf_hwproj = 279167641195`, and forcing-work residual
`8.520423303056633e-12`. Hardened reconciliation closes at `20/20/20`, with
no active reservation, no transaction, and `issues = []`.

A retained fail-closed recost lifecycle companion exists at
`scripts/frontier/cgl_lf_stage_i_checkpoint.py` with SHA-256
`10156515c4bcbfdcf57a2fe54220c2a80f0477f1a7bddbde322b9433946615c2`;
its isolated fixture suite passes (`67 passed`). Historical F-101 remains a
truthfully adopted legacy canonical boundary. F-102 hardened the Stage I
helper before `s19`; the helper SHA-256 is
`54ec671bb45aa27735a174d40b4b2e6009070716346ea09699bbe62421bbfada`,
and focused helper slice passes (`10 passed, 35 deselected`).

F-103 is the first hardened normal observed publication under the retained
companion. Its canonical recost artifact SHA-256 is
`ab7fa80746fa4d2c65f515b3ace71d5668586e27fde91bbdeb578fc3a4509d98`;
its observed-publication audit SHA-256 is
`28e8526eb3c225cde028476dea95657b170264ea7f77156b8a839e2fe03bbf6b`.
It authorizes only `R02/s20_rankio_t7p5_t7p75` on one node with Slurm
`01:05:00`, Athena `00:55:00`, and the reviewed one-segment `2700`-second
threshold. Before preparing `s20`, commit, push, archive, and catalog the
current documentation checkpoint, validate the full source-archive checksum
ledger, then require a fresh empty queue, free strict root lock, authenticated
hardened reconciliation, `20/20/20`, no active reservation, no transaction,
`issues = []`, and an independently audited bounded readiness packet.

The frozen `R02`-`R17` inventory in
`inputs/cgl_lf_paper/mks24_stage_i_manifest.json` is accepted as the Stage I
case inventory for the AthenaK reproduction of Majeski, Kunz, and Squire
(2024, arXiv:2405.02418v2), but the `E02-modal-driver` execution protocol is
not scientifically signed off. An independent source-to-code audit found
that the retained production driver does not exactly implement the published
forcing semantics. Prohibit every new `E02` preparation or submission,
preserve its outputs as pipeline and cost evidence, and use the now-qualified
explicit E03 forcing policies. The reviewed token is retained, reconciliation
passed. Fresh `R02/s00_rankio_t0_t0p1` job `4745922` was then submitted from
`t = 0`, formally inspected, and recorded `accepted` at exact `t = 0.1` for
`0.188889` node-hours. Authenticated continuation
`R02/s01_rankio_t0p1_t0p25` job `4746154` is also formally inspected and
recorded `accepted` through exact `t = 0.25` for `0.289167` node-hours,
bringing corrected E03 Stage I use to `0.478056` node-hours. Authenticated
`R02/s02_rankio_t0p25_t1` job `4746182` is formally inspected and recorded
`accepted` through exact `t = 1.0` for `1.465556` node-hours, bringing
corrected E03 Stage I use to `1.943612` node-hours. Authenticated
`R02/s03_rankio_t1_t1p5` job `4746356` is formally and independently
inspected and recorded `accepted` through exact `t = 1.5` for `1.048611`
node-hours, bringing corrected E03 Stage I use to `2.992223` node-hours.
Authenticated `R02/s04_rankio_t1p5_t2` job `4746435` is formally and
independently inspected and recorded `accepted` through exact `t = 2.0` for
`1.063611` node-hours, bringing corrected E03 Stage I use to `4.055834`
node-hours.
Authenticated `R02/s05_rankio_t2_t2p5` job `4746663` is formally and
independently inspected and recorded `accepted` through exact `t = 2.5` for
`1.062500` node-hours, bringing corrected E03 Stage I use to `5.118334`
node-hours.
Authenticated `R02/s06_rankio_t2p5_t3` job `4746773` is formally and
independently inspected and recorded `accepted` through exact `t = 3.0` for
`1.077500` node-hours, bringing corrected E03 Stage I use to `6.195834`
node-hours.
Authenticated `R02/s07_rankio_t3_t3p5` job `4746953` is formally and
independently inspected and recorded `accepted` through exact `t = 3.5` for
`1.151111` node-hours, bringing corrected E03 Stage I use to `7.346945`
node-hours.
Authenticated `R02/s08_rankio_t3p5_t4` job `4747015` is formally and
independently inspected and recorded `accepted` through exact `t = 4.0` for
`1.200000` node-hours, bringing corrected E03 Stage I use to `8.546945`
node-hours.
Authenticated `R02/s09_rankio_t4_t4p5` job `4747087` is formally and
independently inspected and recorded `accepted` through exact `t = 4.5` for
`1.208056` node-hours, bringing corrected E03 Stage I use to `9.755001`
node-hours.
Authenticated `R02/s10_rankio_t4p5_t5` job `4747146` is formally and
independently inspected and recorded `accepted` through exact `t = 5.0` for
`1.226389` node-hours, bringing corrected E03 Stage I use to `10.981390`
node-hours.
Authenticated `R02/s11_rankio_t5_t5p5` job `4747202` is formally and
independently inspected and recorded `accepted` through exact `t = 5.5` for
`1.268889` node-hours, bringing corrected E03 Stage I use to `12.250279`
node-hours.
Authenticated `R02/s12_rankio_t5p5_t5p75` job `4747500` is formally and
independently inspected and recorded `accepted` through exact `t = 5.75` for
`0.660000` node-hours, bringing corrected E03 Stage I use to `12.910279`
node-hours.
Authenticated `R02/s13_rankio_t5p75_t6` job `4747834` is formally and
independently inspected and recorded `accepted` through exact `t = 6.0` for
`0.599722` node-hours, bringing corrected E03 Stage I use to `13.510001`
node-hours.
Authenticated `R02/s14_rankio_t6_t6p25` job `4748138` is formally and
independently inspected and recorded `accepted` through exact `t = 6.25` for
`0.655278` node-hours, bringing corrected E03 Stage I use to `14.165279`
node-hours.
Authenticated `R02/s15_rankio_t6p25_t6p5` job `4752118` is formally and
independently inspected and recorded `accepted` through exact `t = 6.5` for
`0.653611` node-hours, bringing corrected E03 Stage I use to `14.818890`
node-hours.
Authenticated `R02/s16_rankio_t6p5_t6p75` job `4753294` is formally and
independently inspected and recorded `accepted` through exact `t = 6.75` for
`0.647222` displayed node-hours, bringing corrected E03 Stage I use to
`15.466112` displayed node-hours.
Authenticated `R02/s17_rankio_t6p75_t7` job `4754008` is formally and
independently inspected and recorded `accepted` through exact `t = 7.0` for
`0.667222` displayed node-hours, bringing corrected E03
Stage I use to `16.133334` displayed node-hours.
Authenticated `R02/s18_rankio_t7_t7p25` job `4754394` is formally and
independently inspected and recorded `accepted` through exact `t = 7.25` for
`0.650000` node-hours, bringing corrected E03 Stage I use to `16.783334`
displayed node-hours.
Authenticated `R02/s19_rankio_t7p25_t7p5` job `4757300` is formally and
independently inspected and recorded `accepted` through exact `t = 7.5` for
`0.680556` displayed node-hours, bringing corrected E03 Stage I use to
`17.463890` displayed node-hours.

The `900.000000` node-hour `E02` value is an authorized measurement-based
projected envelope, not a fully measured matrix cost. It no longer authorizes
any new `E02` preparation or submission. This review does not authorize
Stage II, the excluded active-Alfvenic beta-1 inventory deck, or a full-paper
reproduction claim.

The pinned paper source is
`build-cgl-implementation/cgl_lf_reference/arXiv-2405.02418v2-verified/source/MKS24.tex`
with SHA-256
`6d6e748fd1883c5d33167be653d67aed2f84a9b364267d9e90ea075df184af4c`.
The immutable production executable was built from revision
`462b9dbd53e085dea46c2478b567781576d7d03e` and has SHA-256
`df87684e9d2b7af33b36c2757779d15de87b84ef051ca9c4489f2efa358f5c48`.
That E02 executable remains historical pipeline and cost evidence only.
Corrected E03 qualification is closed for revision
`9e07542281e4e6d125582f253df3ad2e3b8b154d`, executable SHA-256
`68f243f9204df388b24365ae65a567f6f567dbe422a6d7a43b9fb4a499ef118c`,
and retained source-bundle SHA-256
`c39d55809989d20aa5438711803f4fd43237fa9284c7e15183e84fdd693d3687`.

## Source-To-Code Contract

| Contract item | Pinned source | AthenaK implementation and submitted controls | Review disposition |
| --- | --- | --- | --- |
| CGL pressure equations | `MKS24.tex:210-239` | `src/eos/cgl_mhd.cpp`; `src/diffusion/cgl_landau_fluid.cpp:154-197` | Matched for Stage I. |
| Landau-fluid closure | `MKS24.tex:230-239` defines the 3+1 `q_perp` and `q_parallel` closure with characteristic `abs(k_parallel)` | `src/diffusion/cgl_landau_fluid.cpp:75,101,172-197,209-285`; nominal decks set `mhd/lf_k_parallel = 6.283185307179586` | Matched for the nominal `4*pi/L_parallel = 2*pi` convention. |
| Electron treatment | `MKS24.tex:239` assumes cold electrons and isotropic background plasma pressure | The Stage I decks evolve the ion CGL-LF model without a separate electron-pressure closure | Matched scope; do not claim weakly collisional electron physics. |
| Active and passive Delta | `MKS24.tex:465` distinguishes active CGL-MHD from passive-Delta isothermal-MHD comparisons | `src/eos/cgl_mhd.cpp:68-73`; `src/pgen/tests/cgl_lf_paper.cpp:159-172`; `src/mhd/rsolvers/hlle_cgl.hpp:194`; active and passive decks set both `mhd/passive` and `problem/passive_delta` consistently | Matched. |
| Mirror and firehose thresholds | `MKS24.tex:465` uses mirror `beta*Delta > 1` and parallel-firehose `beta*Delta < -2`; `MKS24.tex:695` separately discusses the oblique alternative | `src/eos/cgl_physics.hpp:15,49`; every Stage I paper deck explicitly sets `mhd/cgl_firehose_threshold = parallel` | Matched. The code default is intentionally not relied upon. For Figure 13, `R15` supports late-time parallel-versus-oblique reclassification, not a separate oblique execution; the paper reports `10.4%`, `18.5%`, and external hybrid-kinetic `17.9%` context. |
| Hard-wall limiter | `MKS24.tex:467` uses nominal `nu_lim = 1e10 v_A/L_perp` as a hard wall | Standard decks set `mhd/limiter_nu_coll = 1.0e10` and `mhd/limiter_hardwall = true`; `src/eos/cgl_mhd.cpp:75-94` parses the controls; `src/eos/cgl_physics.hpp:49` applies the algebraic projection and `src/diffusion/cgl_landau_fluid.cpp:77` suppresses LF coefficients when the threshold policy requires it | Matched production policy. The hard wall combines algebraic projection with threshold-dependent LF-coefficient suppression; it is not a claim to resolve microinstability kinetics. |
| Periodic elongated box | `MKS24.tex:467` specifies `[L_x,L_y,L_z] = [1,1,2]` with `L_z = L_parallel` and standard `192 x 192 x 384` resolution | Canonical standard decks set mesh extents `[1,1,2]` and resolution `192 x 192 x 384` | Matched. |
| Final time | `MKS24.tex:467` requires at least `t_f = 10 L_perp/v_A` | Canonical matrix decks set `time/tlim = 10.0`; accepted continuations may shorten only the submitted segment target and must preserve the same inspected lineage | Matched. |
| OU forcing | `MKS24.tex:469` specifies velocity forcing, `d_t E_K = 0.32`, `t_corr = L_parallel/v_A = 2`, forced shell `[1,3]` in units of `2*pi/L_parallel`, and `k^-2` power | The archived E02 decks set `dedt = 0.32`, `tcorr = 2.0`, `physical_k_shell = true`, `k_shell_unit = pi`, `nlow = 1`, `nhigh = 3`, and `expo = 2.0`, but omit `spectrum = power_law`; the E02 executable therefore applies its `parabolic` default to random roles. The corrected E03 decks set `spectrum = power_law` explicitly in every paper deck. | Qualified for E03 by explicit-deck audit and corrected-policy jobs `g024` and `g025`; preserve the E02 discrepancy as historical stop-line evidence. |
| Forcing families | `MKS24.tex:469` distinguishes unconstrained random three-component forcing from planar Alfvenic forcing with `grad_perp dot u_perp = 0` | The archived E02 executable applies the generic projection to all paper roles. Random roles are therefore projected rather than unconstrained; planar roles use a full-`abs(k)^2` denominator and do not generally satisfy perpendicular incompressibility for retained `k_z != 0` modes. Corrected E03 adds restart-retained `mks24_random_unprojected` and `mks24_alfvenic_perpendicular` policies, enforces `sol_fraction = 1` for the latter, and sets the matching policy explicitly in every paper deck. | Qualified for E03: `g024` retains zero `f_z` and vanishing perpendicular divergence for retained `k_z != 0` modes; `g025` retains nonzero `f_z` and unprojected full divergence. Generic defaults remain available for non-paper users. |
| Modal forcing restart identity | Required for segmented AthenaK production | `src/srcterms/turb_driver.cpp:1630-1728`, `src/outputs/restart.cpp:294-354`, and `src/pgen/pgen.cpp:175-215` retain and validate modal forcing state. Corrected E03 also records `time/restart_time` in every restart parameter dump so the helper can authenticate the selected terminal restart against inspected physical time before continuation. | Qualified for E03 by `g026`/`g027`: the resumed eight-sibling checkpoint crosses an OU refresh and matches the uninterrupted reference within retained-format tolerances. |

## Corrected E03 Qualification Disposition

Corrected immutable Frontier qualification is closed by reviewed jobs `g024`
through `g031`. `g024` and `g025` qualify explicit planar and random forcing;
`g026`/`g027` close restart identity across an OU refresh; `g028` closes
one-rank/eight-rank GPU decomposition identity; `g029c` closes passive-Delta
semantics; `g030b` reaches exact `t = 2.0` with zero strict counters,
terminal `lf_hwproj = 10084905222`, forcing-work relative residual
`1.8681072370071585e-12`, nine shared snapshots, five checkpoints, and
fourteen MPI-I/O records; `g031` closes standard-layout startup, one-node
memory fit, and ranked retained-output sizing only. Canonical `g031` evidence
JSON SHA-256 is
`a16415c5f9c557a0dad35936c0bb0e89658da9e3b78a4d86852ead4e006f0d98`.
The full reviewed `g024`--`g031` checksum set is recorded in F-078 of
`docs/cgl_lf_mks24_reproduction_implementation_plan.md`. Reviewed approval
token SHA-256
`e3fec9f35da42121b902f41ef752021f375aab38bc34c5ff75ae8780bfbca635`
was retained atomically at `2026-05-30T21:33:31+00:00`.

## Frozen Matrix Review

The canonical matrix contains sixteen unique executions:

| Cases | Purpose | Distinguishing controls |
| --- | --- | --- |
| `R02`-`R09` | Eight paper-standard active/passive, Alfvenic/random, beta-10/beta-100 decks | Standard layout, nominal heat flux, hard-wall limiter |
| `R10`-`R11` | Figure 3 compressive cases | Random forcing; beta `1` and sonic-correlation beta `100`, with `R11 tcorr = 0.2 = L_parallel/v_th` |
| `R12`-`R13` | Figure 12 heat-flux-strength variants | `lf_k_parallel = 0.06283185307179586` and `628.3185307179587`, giving nominal uncapped closure amplitudes 100x stronger and weaker than nominal; applied fluxes still depend on caps and threshold-conditioned denominators |
| `R14`-`R15` | Figure 13 limiter-rate variants | `limiter_nu_coll = 20` and `200`, compared with hard-wall alias `R03` |
| `R16`, `R17` | Figure 11 scale separation | `96 x 96 x 192` and `384 x 384 x 768`, compared with standard-layout alias `R02` |

The explicit manifest aliases and reused mapped roles are part of the case
inventory:

- `R02` supplies standard active-Alfvenic beta-10, nominal-active Figure 12,
  and `n_perp = 192` Figure 11 roles.
- `R03` supplies standard active-Alfvenic beta-100 and nominal hard-wall
  Figure 13 roles.
- `R06` supplies the reused nominal passive Figure 12 role; it is not an
  explicit manifest alias.

The separate
`inputs/cgl_lf_paper/cgl_lf_paper_standard_active_alfvenic_beta1.athinput`
deck remains an inventory definition only. The displayed beta-1 compressive
role is random forcing and is mapped to `R10`; no published-result role has
been established for the active-Alfvenic beta-1 deck. It is therefore
`unmapped_not_authorized`.

## Observable Boundary

Stage I must report panel-level status rather than silently widening the
claim:

| Status | Meaning |
| --- | --- |
| `passed` or `failed` | A checksum-qualified reference product, accepted required cases, and reviewed comparison criterion exist. |
| `not_run` | An admitted comparison still lacks one or more required accepted bundles or products. |
| `blocked_reference` | The matching AthenaK product exists or is planned, but a defensible published-data transform or reference dataset is unavailable. |
| `external_model` | The published panel uses physics or code outside the AthenaK CGL-LF reproduction scope. |

The source defines shell-binned spectra at `MKS24.tex:480-486` but does not
state the discrete Fourier normalization needed to map absolute digitized
ordinates into AthenaK products. Absolute spectral and strain panels remain
`blocked_reference` unless author/archive numeric data, donor diagnostic code,
or an explicit normalization statement is retained. Figure 10 remains
`external_model`: it is a hybrid-kinetic `Pegasus++` comparison, not an
AthenaK CGL-LF result (`MKS24.tex:621,626`).

## Manuscript Appendix Work

The corrected E03 simulation protocol entry gate is closed through Frontier
qualification and token retention. Fresh mapped production has started with
accepted R02 jobs `4745922` and `4746154` through exact `t = 0.25`.
Authenticated job `4746182` now extends that accepted lineage through exact
`t = 1.0`. Authenticated job `4746356` now extends that accepted lineage
through exact `t = 1.5`. Authenticated job `4746435` now extends that
accepted lineage through exact `t = 2.0`. Authenticated job `4746663` now
extends that accepted lineage through exact `t = 2.5`.
Authenticated job `4746773` now extends that accepted lineage through exact
`t = 3.0`.
Authenticated job `4746953` now extends that accepted lineage through exact
`t = 3.5`.
Authenticated job `4747015` now extends that accepted lineage through exact
`t = 4.0`.
Authenticated job `4747087` now extends that accepted lineage through exact
`t = 4.5`.
Authenticated job `4747146` now extends that accepted lineage through exact
`t = 5.0`.
Authenticated job `4747202` now extends that accepted lineage through exact
`t = 5.5`.
Authenticated job `4747500` now extends that accepted lineage through exact
`t = 5.75`.
Authenticated job `4747834` now extends that accepted lineage through exact
`t = 6.0`.
Authenticated job `4748138` now extends that accepted lineage through exact
`t = 6.25`.
Authenticated job `4752118` now extends that accepted lineage through exact
`t = 6.5`.
Authenticated job `4753294` now extends that accepted lineage through exact
`t = 6.75`.
Authenticated job `4754008` now extends that accepted lineage through exact
`t = 7.0`.
Authenticated job `4754394` now extends that accepted lineage through exact
`t = 7.25`.
Authenticated job `4757300` now extends that accepted lineage through exact
`t = 7.5`.
The manuscript also still needs
compact derivation prose. That prose must:

1. State the 3+1 `q_perp` and `q_parallel` closure and the cold-electron,
   isotropic-background assumptions from `MKS24.tex:230-239`.
2. Explain why the nominal constant `abs(k_parallel) = 4*pi/L_parallel` is a
   model closure scale and why the Figure 12 decks vary it by factors of 100.
3. Distinguish active-Delta feedback from passive-Delta diagnostics.
4. Distinguish nominal algebraic hard-wall projection from the finite
   `nu_lim = 20, 200` scan and from unresolved kinetic microinstability
   dynamics.
5. Carry the published reduced-order argument into a short observable map:
   suppressed `grad_parallel Delta p`, pressure-stress transfer, alignment,
   heat-flux insensitivity, scale separation, and limiter sensitivity.
6. State the absolute-spectrum normalization boundary and the Figure 10
   external-model boundary without presenting either as a failed AthenaK
   result.
7. Reproduce the main-text high-beta reduced CGL ordering, stress and
   heat-flux suppression, invariants, and conserved quantities.
8. Reproduce the reduced-kinetic-MHD comparison, separating the signatures
   that survive from the non-local pressure effects that do not.
9. Reproduce the Braginskii collisional-limit argument and
   `beta k_parallel v_A / nu >> 1`, keeping it distinct from the
   threshold-activated `nu_lim` scan.
10. Reproduce the larger-density-fluctuation ordering and identify which
    conclusions no longer follow.

These appendix items remain writing requirements. The forcing discrepancy
separately required a corrected production executable, requalification, and a
fresh execution epoch; the corrected E03 disposition above closes that entry
gate without admitting any E02 restart or any E02 output as a paper-production
result.

## Operational Handoff

Use `docs/cgl_lf_weak_guide_manuscript_plan.md` as the primary guiding
document and `docs/cgl_lf_mks24_reproduction_implementation_plan.md` as the
detailed evidence log. The canonical Frontier root is
`/lustre/orion/ast207/proj-shared/dfielding/CGL`. Preserve the stale
exploratory `beta25-accel05-gamma10001-purecgl-256` campaign record untouched;
explicitly acknowledge it only after a read-only isolation review. Do not
prepare or submit any new `E02-modal-driver` segment. Preserve the
accepted `E02` products as pipeline and cost evidence. Corrected E03
qualification is closed. The reviewed E03 approval token is retained and
reconciliation passed. Fresh `R02/s00_rankio_t0_t0p1` job `4745922` was
formally inspected and recorded `accepted` at exact `t = 0.1`, using
`0.188889` node-hours. F-079/F-080 controller provenance is committed,
archived, and reconciled. Authenticated `R02/s01_rankio_t0p1_t0p25` job
`4746154` is also accepted through exact `t = 0.25` for `0.289167`
node-hours. Its sampled-history forcing-work relative residual is
`1.7741014899016423e-11`; strict LF failure counters remain zero; terminal
`lf_hwproj = 28701760`; complete eight-rank snapshot and restart groups are
retained. Corrected E03 use is `0.478056` node-hours. Retained recost evidence
JSON SHA-256 is
`eb6071cf53d453b2f75707003616ac9cbdfe025d62969d00e34e7948bee3310c`;
it projects `700.868889` provisional matrix node-hours inside the
`900.000000` envelope. Authenticated `R02/s02_rankio_t0p25_t1` job `4746182`
is now accepted through exact `t = 1.0` for `1.465556` node-hours. Its
sampled-history forcing-work relative residual is
`8.616579960442532e-12`; strict LF failure counters remain zero; terminal
`lf_hwproj = 65252911379`; complete terminal eight-rank restart siblings are
retained. Corrected E03 use is `1.943612` node-hours. Updated recost evidence
JSON SHA-256 is
`7d11e8a0004e24417167aeb1f186f8b1f68c8eb9d6206467a708a7cad16cf252`;
preserve it as superseded arithmetic history. The corrected conservative
projection at exact `t = 1.0` is `704.604815` matrix node-hours. Authenticated
`R02/s03_rankio_t1_t1p5` job `4746356` is now accepted through exact
`t = 1.5` for `1.048611` node-hours. Its sampled-history forcing-work
relative residual is `4.978727845741857e-13`; strict LF failure counters
remain zero; terminal `lf_hwproj = 106136472818`; complete two-group
eight-rank snapshots and terminal eight-rank restart siblings are retained.
Corrected E03 use is `2.992223` node-hours. Retained corrected recost evidence
JSON SHA-256 is
`b6fd6e8fcd939ca1baa06411a3cd6f04359b3bc2ff5b2a0097e2cc87b1763fc7`;
it projects at most `724.645556` matrix node-hours inside the `900.000000`
envelope, leaving `175.354444` node-hours of margin. F-083/F-084 controller
hardening is promoted at canonical revision
`9689c269bf329542815a1b2b137881126964b05c`, helper SHA-256
`1c633ebb58294938a0a0609742ae8f8d1f88cd15242ffe649578796edeb39375`.
Focused helper tests, Sphinx warnings-as-errors, syntax, diff, and hardened
reconciliation pass. Authenticated `R02/s04_rankio_t1p5_t2` job `4746435`
is now accepted through exact `t = 2.0` for `1.063611` node-hours. Its
sampled-history forcing-work relative residual is `3.034159691118515e-12`;
strict LF failure counters remain zero; terminal `lf_hwproj = 123445535090`;
complete two-group eight-rank snapshots and terminal eight-rank restart
siblings are retained. Corrected E03 use is `4.055834` node-hours. Retained
corrected recost evidence JSON SHA-256 is
`c382b78ab9e21466648adf9d7ea38d2b407f0e986579f128bdfe6bc44bef79dd`;
it projects at most `728.845556` matrix node-hours inside the `900.000000`
envelope, leaving `171.154444` node-hours of margin. Retained F-087 evidence
historically authorized only the bounded `t = 2.0`--`2.5` continuation.
Authenticated
`R02/s05_rankio_t2_t2p5` job `4746663` is now accepted through exact
`t = 2.5` for `1.062500` node-hours. Its sampled-history forcing-work
relative residual is `3.0678457367645077e-12`; strict LF failure counters
remain zero; terminal `lf_hwproj = 171943591926`; complete two-group
eight-rank snapshots and terminal eight-rank restart siblings are retained.
Corrected E03 use is `5.118334` node-hours. Retained corrected recost
evidence JSON SHA-256 is
`39115585cf7ca866e3aade1e53c69aab6ee9217361d75db7946748905e1c0e5c`;
it projects at most `728.534445` matrix node-hours inside the `900.000000`
envelope, leaving `171.465555` node-hours of margin. Retained F-088 evidence
historically authorized only the bounded `t = 2.5`--`3.0` continuation.
Authenticated `R02/s06_rankio_t2p5_t3` job `4746773` is now accepted through
exact `t = 3.0` for `1.077500` node-hours. Its sampled-history forcing-work
relative residual is `2.8133353495278692e-12`; strict LF failure counters
remain zero; terminal `lf_hwproj = 184131904300`; complete two-group
eight-rank snapshots and terminal eight-rank restart siblings are retained.
Corrected E03 use is `6.195834` node-hours. Retained corrected recost
evidence JSON SHA-256 is
`d2c5e27553426beb03c8c1882bd6473b682029c6ba67acd32bb06c5fb9a3ff01`;
it projects at most `732.734445` matrix node-hours inside the `900.000000`
envelope, leaving `167.265555` node-hours of margin. Retained F-089 evidence
historically authorized only the bounded `t = 3.0`--`3.5` continuation.
Authenticated `R02/s07_rankio_t3_t3p5` job `4746953` is now accepted through
exact `t = 3.5` for `1.151111` node-hours. Its sampled-history forcing-work
relative residual is `2.6602098521945525e-13`; strict LF failure counters
remain zero; terminal `lf_hwproj = 192268745176`; complete two-group
eight-rank snapshots and terminal eight-rank restart siblings are retained.
Corrected E03 use is `7.346945` node-hours. Retained corrected recost
evidence JSON SHA-256 is
`468e787db48d98661cd3ffe125829a78f135bb7ca4c9104c766a00b693f586a3`;
it projects at most `753.345556` matrix node-hours inside the `900.000000`
envelope, leaving `146.654444` node-hours of margin. The accepted clean
segment took `4144` seconds, above the prior `4000`-second operational guard.
Linnaeus independently approved one bounded `s08` continuation under an
explicit `4500`-second ceiling that retains `300` seconds before the Athena
timeout and cannot ratchet automatically. Archive and catalog the current
committed controller state, reconcile, and then prepare only
`R02/s08_rankio_t3p5_t4` from the authenticated terminal siblings were the
historical next actions. That continuation is accepted below.
Authenticated `R02/s08_rankio_t3p5_t4` job `4747015` is now accepted through
exact `t = 4.0` for `1.200000` node-hours. Its sampled-history forcing-work
relative residual is `1.3875872397103518e-13`; strict LF failure counters
remain zero; terminal `lf_hwproj = 203117871743`; complete two-group
eight-rank snapshots and terminal eight-rank restart siblings are retained.
Corrected E03 use is `8.546945` node-hours. Retained corrected recost
evidence JSON SHA-256 is
`dd5b0f5b82cf81005c1a481ee83d1127aa43861dcdd8c22765170e28b0ad6c99`;
it projects at most `767.034445` matrix node-hours inside the `900.000000`
envelope, leaving `132.965555` node-hours of margin. Linnaeus blocked the
first staged artifact because it mislabeled the immediate predecessor
threshold. The regenerated artifact distinguishes the `4000`-second baseline,
the `4500`-second `s08` cap, and the one-segment `4800`-second `s09` cap.
Byte-preserving promotion is independently approved. Archive and catalog the
current committed controller state. Hardened reconciliation closes with
`9/9/9` ledger rows/manifests/reservations, no active reservation, and no
transaction. Reconcile, then prepare only
`R02/s09_rankio_t4_t4p5` from the authenticated terminal siblings with Slurm
walltime `01:40:00`, Athena timeout `01:30:00`, and a `4800`-second threshold
retaining `600` seconds on both timeout margins. The threshold is scoped only
to `s09` and cannot ratchet automatically. Those were the historical next
actions. Authenticated `R02/s09_rankio_t4_t4p5` job `4747087` is now accepted
through exact `t = 4.5` for `1.208056` node-hours. Its sampled-history
forcing-work relative residual is `1.9837412357887464e-12`; strict LF failure
counters remain zero; terminal `lf_hwproj = 207989256307`; complete two-group
eight-rank snapshots and terminal eight-rank restart siblings are retained.
Corrected E03 use is `9.755001` node-hours. Retained corrected recost evidence
JSON SHA-256 is
`8b10b1786db520740b687d989207f3cda6e0123f9bf2b62e75a5df3849511561`;
it projects at most `769.290001` matrix node-hours inside the `900.000000`
envelope, leaving `130.710000` node-hours of margin. Hardened reconciliation
closes with `10/10/10` ledger rows/manifests/reservations, no active
reservation, and no transaction. Archive and catalog the current committed
controller state, reconcile, then prepare only `R02/s10_rankio_t4p5_t5` from
the authenticated `s09` terminal siblings with Slurm walltime `01:40:00`,
Athena timeout `01:30:00`, and a `4800`-second threshold retaining `600`
seconds on both timeout margins. The threshold is scoped only to `s10` and
cannot ratchet automatically. At that boundary, those were the historical next
actions.
Authenticated `R02/s10_rankio_t4p5_t5` job `4747146` is now accepted through
exact `t = 5.0` for `1.226389` node-hours. Its sampled-history forcing-work
relative residual is `5.300795515586916e-12`; strict LF failure counters
remain zero; terminal `lf_hwproj = 220146873111`; complete two-group
eight-rank snapshots and terminal eight-rank restart siblings are retained.
Corrected E03 use is `10.981390` node-hours. Retained corrected recost
evidence JSON SHA-256 is
`5256fb7fd9b75f175a16c2c16ee8a995b8d6d9281a92b70a979ea62171835bf2`;
it projects at most `774.423334` matrix node-hours inside the `900.000000`
envelope, leaving `125.576666` node-hours of margin. Hardened reconciliation
closes with `11/11/11` ledger rows/manifests/reservations, no active
reservation, and no transaction. At that boundary, the historical next action
was the reviewed `s11` probe.
Authenticated `R02/s11_rankio_t5_t5p5` job `4747202` is now accepted through
exact `t = 5.5` for `1.268889` node-hours. Its sampled-history forcing-work
relative residual is `5.093197346120908e-12`; strict LF failure counters
remain zero; terminal `lf_hwproj = 224321212896`; complete two-group
eight-rank snapshots and terminal eight-rank restart siblings are retained.
Corrected E03 use is `12.250279` node-hours. Retained corrected recost
evidence JSON SHA-256 is
`323d7cb3cbe108eee944045c189cff75ffc07d06833c49bcfa8315f07fad8d51`;
it projects at most `786.323334` matrix node-hours inside the `900.000000`
envelope, leaving `113.676666` node-hours of margin. Hardened reconciliation
closes with `12/12/12` ledger rows/manifests/reservations, no active
reservation, and no transaction. At that boundary, the historical next action
was to archive and catalog the committed controller state, reconcile, then
prepare only `R02/s12_rankio_t5p5_t5p75` on
one node from the authenticated `s11` terminal siblings with Slurm walltime `01:05:00`,
Athena timeout `00:55:00`, and a `2700`-second threshold retaining `600`
seconds on both timeout margins. The threshold is scoped only to `s12` and
cannot ratchet automatically. Any scientific,
provenance, scheduler, storage, budget, or reconciliation failure blocks
successor preparation. Inspect, account, reconcile, and recost before any
further extension. Finish R02, execute R03--R16 sequentially, and complete
`R17` last.
Authenticated `R02/s12_rankio_t5p5_t5p75` job `4747500` is now accepted
through exact `t = 5.75` for `0.660000` node-hours. Its sampled-history
forcing-work relative residual is `6.6019081979506825e-12`; strict LF failure
counters remain zero; terminal `lf_hwproj = 235196163710`; one complete
eight-rank snapshot group and terminal eight-rank restart siblings are
retained. Corrected E03 use is `12.910279` node-hours. Retained corrected
recost evidence JSON SHA-256 is
`d0a60e222138971e1bc9aae978bb3d62aee460f09000ee62f7f9ffdfe534165f`;
it projects at most `814.945556` matrix node-hours inside the `900.000000`
envelope, leaving `85.054444` node-hours of margin. Hardened reconciliation
closes with `13/13/13` ledger rows/manifests/reservations, no active
reservation, and no transaction. At that boundary, the historical next action
was to archive and catalog the committed controller state, reconcile, then
prepare only `R02/s13_rankio_t5p75_t6` on
one node from the authenticated `s12` terminal siblings with Slurm walltime
`01:05:00`, Athena timeout `00:55:00`, and a `2700`-second threshold retaining
`600` seconds on both timeout margins. The threshold is scoped only to `s13`
and cannot ratchet automatically. Any scientific, provenance, scheduler,
storage, budget, or reconciliation failure blocks successor preparation.
Inspect, account, reconcile, and recost before any further extension. Finish
R02, execute R03--R16 sequentially, and complete `R17` last.
Authenticated `R02/s13_rankio_t5p75_t6` job `4747834` is now accepted
through exact `t = 6.0` for `0.599722` node-hours. Its sampled-history
forcing-work relative residual is `5.579303003548508e-12`; strict LF failure
counters remain zero; terminal `lf_hwproj = 239804982453`; one complete
eight-rank snapshot group and terminal eight-rank restart siblings are
retained. Corrected E03 use is `13.510001` node-hours. Retained corrected
recost evidence JSON SHA-256 is
`4bc19f5b44587e73869982425b1cf0ad459547c3e0c5de1ab1b2b5cb366c6322`;
it projects at most `814.945556` matrix node-hours inside the `900.000000`
envelope, leaving `85.054444` node-hours of margin. Hardened reconciliation
closes with `14/14/14` ledger rows/manifests/reservations, no active
reservation, and no transaction. At that boundary, the historical next action
was to archive and catalog the committed controller state, reconcile, then
prepare only `R02/s14_rankio_t6_t6p25` on
one node from the authenticated `s13` terminal siblings with Slurm walltime
`01:05:00`, Athena timeout `00:55:00`, and a `2700`-second threshold retaining
`600` seconds on both timeout margins. The threshold is scoped only to `s14`
and cannot ratchet automatically. Any scientific, provenance, scheduler,
storage, budget, or reconciliation failure blocks successor preparation.
Inspect, account, reconcile, and recost before any further extension. Finish
R02, execute R03--R16 sequentially, and complete `R17` last.
Authenticated `R02/s14_rankio_t6_t6p25` job `4748138` is now accepted
through exact `t = 6.25` for `0.655278` node-hours. Its sampled-history
forcing-work relative residual is `8.23749299161874e-12`; strict LF failure
counters remain zero; terminal `lf_hwproj = 246541685689`; one complete
eight-rank snapshot group and terminal eight-rank restart siblings are
retained. Corrected E03 use is `14.165279` node-hours. Retained corrected
recost evidence JSON SHA-256 is
`5ba7aec52d3e7824cf4ccf52f95ce1e46368c84ed59f0525e2bbbcfc167fde91`;
it projects at most `829.101112` matrix node-hours inside the `900.000000`
envelope, leaving `70.898888` node-hours of margin. Hardened reconciliation
closes with `15/15/15` ledger rows/manifests/reservations, no active
reservation, and no transaction. At that boundary, the historical next action
was to archive and catalog the committed controller state, reconcile, then
prepare only `R02/s15_rankio_t6p25_t6p5` on one node from the authenticated
`s14` terminal siblings with Slurm walltime `01:05:00`, Athena timeout
`00:55:00`, and a `2700`-second threshold retaining `600` seconds on both
timeout margins. The threshold was scoped only to `s15` and could not ratchet
automatically. Any scientific, provenance, scheduler, storage, budget, or
reconciliation failure blocks successor preparation. Inspect, account,
reconcile, and recost before any further extension. Finish R02, execute
R03--R16 sequentially, and complete `R17` last.
Authenticated `R02/s15_rankio_t6p25_t6p5` job `4752118` is now accepted
through exact `t = 6.5` for `0.653611` node-hours. Its formal inspection
SHA-256 is
`fd1fa6f3e4218b611814b71a9081140d91f2c1d866a2f6f565de33ec040c7492`;
recorded manifest SHA-256 is
`5cc7a5b4d67cd4f3926f5c435e23139b1103b35f5ef58d9acf892e11bb35e7f8`;
retained independent validation SHA-256 is
`cba6a5ff33e7c5920874e78e55c8288e7863de43b80a1accfa018d45e8547968`;
retained scheduler evidence SHA-256 is
`f1bafbfc2bd875e159b7ce80ef1b7415d464da812060d3a605d13c90192fe40a`.
Its sampled-history forcing-work relative residual is
`7.631778891465068e-12`; accepted-prefix residual is
`3.0155808676444225e-12`; strict LF failure counters remain zero; terminal
`lf_hwproj = 253499800091`; one complete eight-rank snapshot group and terminal
eight-rank restart siblings are retained. Retained quarter validator SHA-256 is
`0ec3004247f62e2266939c94ef07cdb94513d6e8edc2ceddcf0504e3ac812789`;
retained one-use recost generator SHA-256 is
`317a9ff15b8446cdc09e5e43a18dd261a80ae69ff274c28b831ec5279dc32e78`.
Job `4752118` used source launch bundle
`athenak-feature-cgl-through-368bc86e4.bundle` with SHA-256
`49a26705e9e34a203271e5a9e1b330d42bf97a8d64d1ad557d6edbc47ee76c68`.
Corrected E03 use is `14.818890` node-hours. Retained corrected recost evidence
JSON SHA-256 is
`3603077400bd4530001155b83bfc4ac1053ec18505c56a6c29ffadc795bd1255`;
it projects at most `829.101112` matrix node-hours inside the `900.000000`
envelope, leaving `70.898888` node-hours of margin. The accepted `2353`-second
interval leaves `347` seconds beneath its reviewed threshold, and the reviewed
`s16` estimator selects `2559 <= 2700` seconds. Hardened reconciliation closes
with `16/16/16` ledger rows/manifests/reservations, no active reservation, and
no transaction. At that boundary, the historical next action was to archive
and catalog the committed controller state, reconcile, then prepare only
`R02/s16_rankio_t6p5_t6p75` on one node from the authenticated `s15` terminal
siblings with Slurm walltime `01:05:00`, Athena timeout `00:55:00`, and a
`2700`-second threshold retaining `600` seconds on both timeout margins. The
threshold was scoped only to `s16` and could not ratchet automatically. Any
scientific, provenance, scheduler, storage, budget, or reconciliation failure
blocks successor preparation. Inspect, account, reconcile, and recost before
any further extension. Finish R02, execute R03--R16 sequentially, and complete
`R17` last.
Authenticated `R02/s16_rankio_t6p5_t6p75` job `4753294` is now accepted
through exact `t = 6.75` for `0.647222` displayed node-hours. Its exact
scheduler row is
`4753294|cgl_mks24_E03_forcing_policy_R02_s16_rankio_t6p5_t6p75|COMPLETED|0:0|1|2330|2026-06-01T20:40:10|2026-06-01T21:19:07`;
retained scheduler evidence SHA-256 is
`11b3d7bbae54e4c51fa7548134ae0e7e837fe8b624af698a2e32f30c9cd82fc9`.
Its previously audited, non-retained submitted-state manifest SHA-256 is
`10c469c424159ecc3ca283a2c951aa2c3456d9b7fbd01841926ee55b973c42f0`;
formal inspection SHA-256 is
`1149de945447eae3b274052985363c1c821720a117254662657e33d8b21467de`;
recorded manifest SHA-256 is
`a87787fdaf4a0a8491bdf42e834fba7399f084ffef1cd9ff5e9c4aae84037089`;
retained independent validation SHA-256 is
`4889d08083004d3d89f6cc8beede3521ce7f8143983c7b4bb233209d59844b05`.
Its sampled-history forcing-work relative residual is
`6.437727059090897e-12`; accepted-prefix residual is
`3.1385227986388927e-12`; strict LF failure counters remain zero; terminal
`lf_hwproj = 259077916291`; one complete eight-rank snapshot group and
terminal eight-rank restart siblings are retained. The accepted prefix totals
`55678` seconds (`15.466111111111111` exact and `15.466112` displayed
node-hours). Retained quarter validator SHA-256 is
`0ec3004247f62e2266939c94ef07cdb94513d6e8edc2ceddcf0504e3ac812789`;
retained one-use recost generator SHA-256 is
`e9991ec75559113fb0da5b5fae5d117f7fe6701160316237227169dde04d05af`;
reviewed one-shot, non-retained recost promoter SHA-256 is
`08fae00aa1cc365cb39d0b935f059a732955465eb26977e7fd459791ae762883`.
Job `4753294` used source launch bundle
`athenak-feature-cgl-through-6c0739806.bundle` with SHA-256
`a837ebe60d62a922dde3fc3ddd87350eec8a9e3f439f671f783346ee6290617f`.
Retained corrected recost evidence JSON SHA-256 is
`5b997623a3f8d83c200034f42e0c9f1b9a09c811bcb46e2fc1a0f7c1eb58815e`;
it projects at most `829.1011116666666` exact (`829.101112` displayed) matrix
node-hours inside the `900.000000` envelope, leaving `70.89888833333339`
exact (`70.898888` displayed) node-hours of margin. The accepted `2330`-second
interval leaves `370` seconds beneath its reviewed threshold. With prior
selected `2559`, previous observed `2353`, observed `2330`, and local `2330`,
the reviewed `s17` estimator selects `2559 <= 2700` seconds. Hardened
reconciliation closes with `17/17/17` ledger rows/manifests/reservations, no
active reservation, no transaction, and `issues = []`. Historical
authorization was only `R02/s17_rankio_t6p75_t7` on one node from the
authenticated `s16` terminal siblings under the full F-099 profile.
Authenticated `R02/s17_rankio_t6p75_t7` job `4754008` is now accepted through
exact `t = 7.0` for `0.667222` displayed node-hours. Its
exact scheduler row is `4754008|cgl_mks24_E03_forcing_policy_R02_s17_rankio_t6p75_t7|COMPLETED|0:0|1|2402|2026-06-01T22:24:23|2026-06-01T23:04:53`; retained scheduler evidence
SHA-256 is `4ca296c0eb44863dbde4fec34e35b145eba7691df6938688f0c55f2752560fae`. Its previously audited, non-retained
submitted-state manifest SHA-256 is
`a1d986135c397f68a986accff0a4db6621caf644f57cc8c5ae720ef35dcaeea6`;
formal inspection SHA-256 is `f5028f8d6c6157e837b34e06712d68360040bfebdfb5bb1bfdc0221cd8431db8`; recorded
manifest SHA-256 is `73b970e5014ac425f805013ad2901935b7a884bbbcde54ba9a785a53222d63ea`; retained independent
validation SHA-256 is `af4e188ae27904e41b6bf7d517784688ac2b46be03de70d7285c96a4b1e8ccdd`. Its sampled-history
forcing-work relative residual is `7.10881221809853e-12`;
accepted-prefix residual is `3.2704293797689304e-12`; strict LF
failure counters remain zero; terminal `lf_hwproj = 265802141460`;
one complete eight-rank snapshot group and terminal eight-rank restart siblings
are retained. The accepted prefix totals `58080`
seconds (`16.133333333333333` exact and
`16.133334` displayed node-hours). Retained quarter
validator SHA-256 is
`0ec3004247f62e2266939c94ef07cdb94513d6e8edc2ceddcf0504e3ac812789`;
retained one-use recost generator SHA-256 is `b55c741516ccd7e5be236856599e6abf6adf859bb2ebacc801e6164d994db4aa`;
reviewed one-shot, non-retained recost promoter SHA-256 is
`b3d62c75bd175a73df8dcf42a607a37551774261c2bca1193b0a2be125e5d8e4`. Job `4754008` used
source launch bundle `athenak-feature-cgl-through-e119e2dcf.bundle` with SHA-256
`25ffc9cd279a02debf8b55e318f0674ac3c14bf9a05cf05827218ae2056a5f14`.
Retained corrected recost evidence JSON SHA-256 is `63438d8843b86815a2e25d4d6de6ef78ce9756c561fe9db3e762965ada73202f`; it
projects at most `829.1011116666666` exact
(`829.101112` displayed) matrix node-hours
inside the `900.000000` envelope, leaving
`70.89888833333339` exact
(`70.898888` displayed) node-hours of margin.
The accepted `2402`-second interval leaves
`298` seconds beneath its reviewed threshold.
With prior selected `2559`, previous observed `2330`, observed
`2402`, and local `2474`, the
reviewed `s18` estimator selects `2559 <= 2700`
seconds. Hardened reconciliation closes with `18/18/18` ledger
rows/manifests/reservations, no active reservation, no transaction, and
`issues = []`.
Commit the post-F-100 documentation state, archive
and catalog its resulting actual controller bundle and its SHA-256,
rerun reconciliation, repeat the queue/shared-root audit with explicit
acknowledgement of the stale beta-25 record, then prepare only
`R02/s18_rankio_t7_t7p25` on one node from the authenticated `s17` terminal
siblings with Slurm walltime `01:05:00`, Athena timeout `00:55:00`, and a
`2700`-second threshold retaining `600` seconds on both timeout margins. The
threshold is scoped only to `s18` and cannot ratchet automatically. Any
scientific, provenance, scheduler, storage, budget, or reconciliation failure
blocks successor preparation. Inspect, account, reconcile, and recost before
any further extension. Finish R02, execute R03--R16 sequentially, and complete
`R17` last.

That F-100 authorization is historical. Authenticated
`R02/s18_rankio_t7_t7p25` job `4754394` is accepted through exact `t = 7.25`
for `2340` seconds (`0.650000` node-hours). Historical F-101 canonical recost
evidence SHA-256
`df208829aec85c88fcc2caa757a21ad3ff3372aeb82628cbfe7bd41851faed9a`
was truthfully adopted under the retained companion; its audit SHA-256 is
`9bbfc6ee1007b0e17a4408444acb1bcdb2f59cbde9183569d506fe2606c1837d`.
F-102 then hardened lock opening and scheduler control-plane routing before any
`s19` mutation. Authenticated `R02/s19_rankio_t7p25_t7p5` job `4757300` is
accepted through exact `t = 7.5` for `2450` seconds (`0.680556` displayed
node-hours). Its exact scheduler row is
`4757300|cgl_mks24_E03_forcing_policy_R02_s19_rankio_t7p25_t7p5|COMPLETED|0:0|1|2450|2026-06-02T20:18:10|2026-06-02T20:59:12`;
formal inspection SHA-256 is
`8b9cfb6f3d0ec215234d9481e031299da1286117424d756f654270930881f744`;
recorded manifest SHA-256 is
`30445e1c64dc35c64db38f119ff6838c904ee4d71f78ec75b219c4d476bb7b30`;
retained independent validation SHA-256 is
`6ebd5054e8d8082d75343b9cd176d7faf8c3b773774d8e145302bf3bfd859dcf`;
retained scheduler evidence SHA-256 is
`9a45c253369f87f8f7a31ce1189ac2df5b312b293c058a6f5bd9aa2c9e953f48`.
It retains one complete eight-rank snapshot group and terminal eight-rank
restart siblings, zero strict LF failure counters, finite synchronized
histories, no restart-marker bypass, terminal `lf_hwproj = 279167641195`, and
sampled-history forcing-work relative residual `8.520423303056633e-12`.
The accepted prefix totals `62870` seconds (`17.46388888888889` exact and
`17.463890` displayed node-hours). Hardened reconciliation closes with
`20/20/20` ledger rows/manifests/reservations, no active reservation, no
transaction, and `issues = []`.

The retained F-103 generator SHA-256 is
`131ab294888fd70bad7140a781ebc43ea6cb05bc4eddcb15dbee91c2a445ba28`.
Canonical F-103 recost evidence SHA-256 is
`ab7fa80746fa4d2c65f515b3ace71d5668586e27fde91bbdeb578fc3a4509d98`;
observed-publication audit SHA-256 is
`28e8526eb3c225cde028476dea95657b170264ea7f77156b8a839e2fe03bbf6b`.
The normal retained-companion publication leaves canonical-only bytes, a
mode-`0444` forensic copy, empty transaction directories, an empty user queue,
and a free strict lock. With prior selected `2559`, previous observed `2340`,
observed `2450`, and local acceleration `2560`, the reviewed `s20` estimator
selects `2560 <= 2700` seconds. It projects `829.2566672222222` matrix
node-hours with `70.74333277777782` margin. Commit, push, archive, and catalog
the post-F-103 documentation checkpoint, validate the archive ledger, rerun
hardened reconciliation, and repeat the queue/shared-root audit with explicit
acknowledgement of the stale beta-25 record. Require completion of the durable
handoff checklist items 6--8 before preparation. Then prepare only
`R02/s20_rankio_t7p5_t7p75` on one node from authenticated `s19` terminal
siblings with Slurm walltime `01:05:00`, Athena timeout `00:55:00`, and a
reviewed one-segment `2700`-second threshold retaining `600` seconds on both
timeout margins. The threshold is scoped only to `s20` and cannot ratchet
automatically. Any scientific, provenance, scheduler, storage, budget, or
reconciliation failure blocks successor preparation. Inspect, account,
reconcile, and recost before any further extension. Finish R02, execute
R03--R16 sequentially, and complete `R17` last.

The first R02 preflight exposed a nonblocking preview-rendering defect:
`check-submit` enforced the reviewed shared-root acknowledgement but omitted it
from the printed follow-up command. Atomic submission retained the required
acknowledgement, so job `4745922` is valid. F-079 records the corrected
equals-style rendering and parser-round-trip regression required before any
later submission.

Applying that fix after the pilot record exposed a second lifecycle boundary:
one mutable live-helper checksum cannot authenticate both retained history and
future helper revisions. F-080 keeps `prepared` and `submitted` manifests
strict against live helper bytes, while `recorded` manifests authenticate each
historical helper blob from its checksum-bound source bundle. Final controller
transition `05cb4c324bfd8feec72ebdeb33b1961c9fde70bf` is archived in
`athenak-feature-cgl-through-05cb4c324.bundle` with SHA-256
`1381918e471730d8c9639014475566f93fc87bac66b62738dd8316fc69c03570`.
Retained F-080 evidence JSON SHA-256 is
`46fc1c4054e75f4224302be1895c1b9eaae88097544512ac067c53376e542e34`;
post-archive reconciliation passes. The first post-transition continuation
was accepted as job `4746154`; later authenticated continuations `4746182`
and `4746356` reach exact `t = 1.5`; authenticated continuation `4746435`
reaches exact `t = 2.0`; authenticated continuation `4746663` now reaches
exact `t = 2.5`; authenticated continuation `4746773` now reaches exact
`t = 3.0`; authenticated continuation `4746953` now reaches exact `t = 3.5`.
Authenticated continuation `4747015` now reaches exact `t = 4.0`.
Authenticated continuation `4747087` now reaches exact `t = 4.5`.
Authenticated continuation `4747146` now reaches exact `t = 5.0`.
Authenticated continuation `4747202` now reaches exact `t = 5.5`.
Authenticated continuation `4747500` now reaches exact `t = 5.75`.
Authenticated continuation `4747834` now reaches exact `t = 6.0`.
Authenticated continuation `4748138` now reaches exact `t = 6.25`.
Authenticated continuation `4752118` now reaches exact `t = 6.5`.
Authenticated continuation `4753294` now reaches exact `t = 6.75`.
Authenticated continuation `4754008` now reaches exact `t = 7.0`.
Authenticated continuation `4754394` now reaches exact `t = 7.25`.
Authenticated continuation `4757300` now reaches exact `t = 7.5`.
After the post-F-103 documentation checkpoint is committed, pushed, archived,
and cataloged, the archive ledger is validated, reconciliation is rerun, and
the queue/shared-root audit repeats explicit acknowledgement of the stale
beta-25 record. Require completion of the durable handoff checklist items
6--8, then use only job `4757300`'s authenticated terminal checkpoint for the
next bounded R02 continuation `R02/s20_rankio_t7p5_t7p75`.
Commits `9480e62764528a3f40066d22a192f0e99b369891`
and `ef1e42fa088203ac9ef6ec8e47e668db4fb95a3c` additionally require the live
helper to remain committed during historical authentication, reserve the
case-level `analysis/` evidence namespace, and reject `--segment analysis`.
The focused Stage I helper subset passes (`10 passed`). Job `4746182` used
the verified `athenak-feature-cgl-through-ef1e42fa.bundle` launch bundle with
SHA-256
`80c211c41a8de32688c3be577d8c727ec3268c347c2a7c3b2d5143e1c34593ec`.
