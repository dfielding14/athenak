# MKS24 Stage I Protocol Review: Stop-The-Line Audit

## Decision

The frozen `R02`-`R17` inventory in
`inputs/cgl_lf_paper/mks24_stage_i_manifest.json` is accepted as the Stage I
case inventory for the AthenaK reproduction of Majeski, Kunz, and Squire
(2024, arXiv:2405.02418v2), but the `E02-modal-driver` execution protocol is
not scientifically signed off. An independent source-to-code audit found
that the retained production driver does not exactly implement the published
forcing semantics. Prohibit every new `E02` preparation or submission,
preserve its outputs as pipeline and cost evidence, implement and qualify
explicit MKS24 forcing policies, then start a fresh execution epoch from
`t = 0`.

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
| OU forcing | `MKS24.tex:469` specifies velocity forcing, `d_t E_K = 0.32`, `t_corr = L_parallel/v_A = 2`, forced shell `[1,3]` in units of `2*pi/L_parallel`, and `k^-2` power | The archived E02 decks set `dedt = 0.32`, `tcorr = 2.0`, `physical_k_shell = true`, `k_shell_unit = pi`, `nlow = 1`, `nhigh = 3`, and `expo = 2.0`, but omit `spectrum = power_law`; the E02 executable therefore applies its `parabolic` default to random roles. The corrected worktree sets `spectrum = power_law` explicitly in every paper deck. | **Blocked pending corrected-build qualification:** preserve the E02 finding as historical audit evidence, qualify the corrected worktree, and start a fresh epoch. |
| Forcing families | `MKS24.tex:469` distinguishes unconstrained random three-component forcing from planar Alfvenic forcing with `grad_perp dot u_perp = 0` | The archived E02 executable applies the generic projection to all paper roles. Random roles are therefore projected rather than unconstrained; planar roles use a full-`abs(k)^2` denominator and do not generally satisfy perpendicular incompressibility for retained `k_z != 0` modes. The corrected worktree adds restart-retained `mks24_random_unprojected` and `mks24_alfvenic_perpendicular` policies, enforces `sol_fraction = 1` for the latter, and sets the matching policy explicitly in every paper deck. | **Blocked pending corrected-build qualification:** preserve generic defaults for non-paper users, qualify the corrected policies, and start a fresh epoch. |
| Modal forcing restart identity | Required for segmented AthenaK production | `src/srcterms/turb_driver.cpp:1630-1728`, `src/outputs/restart.cpp:294-354`, and `src/pgen/pgen.cpp:175-215` retain and validate modal forcing state. The corrected worktree also records `time/restart_time` in every restart parameter dump so the E03 helper can authenticate the selected terminal restart against the inspected physical time before continuation. | Modal-state mechanics were qualified in the archived `E02-modal-driver` epoch. The forcing-correct E03 executable and its explicit restart-time gate still require Frontier qualification. |

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

The simulation protocol is blocked on forcing correction. The manuscript also
still needs compact derivation prose. That prose must:

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

These appendix items are writing requirements. The forcing discrepancy
separately requires a corrected production executable, requalification, and a
fresh execution epoch.

## Operational Handoff

Use `docs/cgl_lf_weak_guide_manuscript_plan.md` as the primary guiding
document and `docs/cgl_lf_mks24_reproduction_implementation_plan.md` as the
detailed evidence log. The canonical Frontier root is
`/lustre/orion/ast207/proj-shared/dfielding/CGL`. Preserve the stale
exploratory `beta25-accel05-gamma10001-purecgl-256` campaign record untouched;
explicitly acknowledge it only after a read-only isolation review. Do not
prepare or submit any new `E02-modal-driver` segment. Preserve the
accepted `E02` products as pipeline and cost evidence. After forcing
correction and replacement-driver qualification, start a fresh epoch from
`t = 0`, run one formally inspected segment at a time, and complete `R17`
last.
