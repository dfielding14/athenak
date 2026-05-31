# MKS24 Stage I Protocol Review: Stop-The-Line Audit

## Decision

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
accepted lineage through exact `t = 2.0`.
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
envelope, leaving `171.154444` node-hours of margin. Archive and catalog the
current committed controller state, reconcile, then prepare only
`R02/s05_rankio_t2_t2p5` from the
authenticated terminal siblings. Inspect, account, reconcile, and recost
before any further extension. Finish R02, execute R03--R16 sequentially, and
complete `R17` last.

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
now reaches exact `t = 2.0`. Use only job `4746435`'s
authenticated terminal checkpoint for the next bounded R02 continuation.
Commits `9480e62764528a3f40066d22a192f0e99b369891`
and `ef1e42fa088203ac9ef6ec8e47e668db4fb95a3c` additionally require the live
helper to remain committed during historical authentication, reserve the
case-level `analysis/` evidence namespace, and reject `--segment analysis`.
The focused Stage I helper subset passes (`10 passed`). Job `4746182` used
the verified `athenak-feature-cgl-through-ef1e42fa.bundle` launch bundle with
SHA-256
`80c211c41a8de32688c3be577d8c727ec3268c347c2a7c3b2d5143e1c34593ec`.
