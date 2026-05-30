# Q-034 Unsupported-Capability Local Closure Sidecar

Date: 2026-05-30

Qualification effect: `descriptive_local_closure_only`

This sidecar audits the local handling of every row in the explicit
unsupported-capability register. It records checked-in local artifacts and the
remaining closure requirement. It does not promote engineering proxies,
extension source smokes, or control-plane scaffolding into scientific or
production qualification.

## Parser-Suite Evidence

The bounded host parser suite is:

- `inputs/tests/pic_parser_contract_guards.athinput`
- `tst/scripts/particles/pic_parser_contract_guards.py`

Direct host execution on 2026-05-30 passed:

```text
pic_parser_contract_guards: PASS
{'positive_cases': 4, 'rejection_cases': 22, 'total_cases': 26}
```

The zero-step suite accepts coherent `paper_test_particle`, `paper_mhd_pic`,
`extended_mhd_pic`, and adaptive-delta-f extension configurations. It checks
precise rejection diagnostics for unsupported parser values and for
cross-model misuse, including direct-current CT in paper mode, relativistic
MHD feedback in paper mode, Hall use outside extension mode, ion-neutral
friction misuse, and adaptive-delta-f misuse. This is bounded serial host
evidence, not a Frontier, MPI-rank-launch, or GPU qualification result.

## Register Audit

### 1. Full electromagnetic PIC physics

- Register handling: `rename_proxy`, `docs_limit`
- Local disposition: `local_handling_present_residual_review`
- Local artifacts:
  - `tst/publication/pic_publication_manifest.py` classifies every checked-in
    case as `engineering_proxy`; the EM-vacuum case maps to
    `mhd_linear_wave_with_inactive_particles`.
  - `tst/publication/plot_pic_publication_figures.py` emits explicitly
    unqualified engineering figures, including a visible proxy watermark.
- Residual closure requirement: retain the explicit known-limitations wording
  in every exported bundle and external figure review. No checked-in proxy may
  be cited as a full electromagnetic-PIC qualification artifact.

### 2. Self-consistent Langmuir-wave physics

- Register handling: `rename_proxy`, `docs_limit`
- Local disposition: `local_handling_present_residual_review`
- Local artifacts:
  - `tst/publication/pic_publication_manifest.py` maps the compatibility-named
    Langmuir case to `nonrelativistic_uniform_b_gyrofrequency_anchor`.
  - `tst/publication/plot_pic_publication_figures.py` titles the output as a
    Langmuir proxy current oscillation.
- Residual closure requirement: preserve the narrow-oracle wording in exported
  captions and known limitations. The compatibility filename must not be
  interpreted as self-consistent Langmuir-wave qualification.

### 3. Physical two-stream or Weibel fidelity

- Register handling: `rename_proxy`, `docs_limit`
- Local disposition: `local_handling_present_residual_review`
- Local artifacts:
  - `tst/publication/pic_publication_manifest.py` maps the cases to
    `two_stream_like_engineering_control` and
    `transverse_current_weibel_like_engineering_control`.
  - `tst/publication/compute_pic_publication_metrics.py` labels generated
    metrics `engineering_proxy` and `unqualified_exploratory_metrics`.
  - `tst/publication/plot_pic_publication_figures.py` generates unqualified
    proxy-watermarked figures.
- Residual closure requirement: retain proxy labeling in every export and
  reviewed caption. Separate physical implementations and qualification
  campaigns are required before making two-stream or Weibel fidelity claims.

### 4. Hall-dominated Bell or shock-front behavior

- Register handling: `extension_gate`
- Local disposition: `extension_gate_open`
- Local artifacts:
  - `docs/source/engineering/pic_mhd_model_contract.md` isolates
    `pic_cr_hall_mode=current_to_ct_experimental` under `extended_mhd_pic` and
    documents its narrow source normalization.
  - `tst/scripts/particles/pic_extended_hall_ct_smoke.py` is the bounded
    odd-in-coefficient manufactured-source smoke.
  - `tst/scripts/particles/pic_parser_contract_guards.py` rejects Hall use in
    paper mode and Hall use without coupled moments.
- Residual closure requirement: derive and preregister the physical Hall
  normalization, then complete linear/nonlinear Bell, shock-front,
  decomposition, MPI, and GPU qualification. Until then, reject Hall-dominated
  Bell and shock-front claims.

### 5. Injection from the thermal pool

- Register handling: `docs_limit`
- Local disposition: `unsupported_simplified_prescription_only`
- Local artifacts:
  - `src/pgen/tests/pic_parallel_shock.cpp` identifies the generator as a
    parallel-shock benchmark and implements eta-based shock-surface CR
    injection with conservative gas subtraction.
  - `inputs/tests/pic_parallel_shock_*.athinput` expose
    `problem/ps_enable_injection` for the unqualified parallel-shock scaffold.
- Residual closure requirement: every Section 5.4 reproduction manifest must
  state that injection is a simplified supra-thermal shock-surface
  prescription. It is not thermal-pool injection microphysics.

### 6. Relativistic MHD background fluid

- Register handling: `parser_reject`, `docs_limit`
- Local disposition: `parser_reject_present_residual_review`
- Local artifacts:
  - `tst/scripts/particles/pic_parser_contract_guards.py` verifies that
    `paper_mhd_pic` with `coord/special_rel=true` fails with the bounded
    non-relativistic-MHD diagnostic.
  - `tst/scripts/particles/pic_mhd_current_coupling.py` retains broader
    relativistic-coupling rejection checks.
  - `docs/source/engineering/pic_mhd_model_contract.md` scopes the supported
    paper equations and extension boundary.
- Residual closure requirement: archive the parser-suite output with the
  release evidence and retain an explicit supported-mode limitation. A
  separately named implementation and qualification campaign are required
  before claiming a relativistic MHD background.

### 7. Physical damping and calibrated CR transport coefficients

- Register handling: `extension_gate`
- Local disposition: `extension_gate_open`
- Local artifacts:
  - `docs/source/engineering/pic_mhd_model_contract.md` documents the bounded
    static-neutral ion-friction reduction and adaptive bi-kappa delta-f fit.
  - `tst/scripts/particles/pic_ion_neutral_friction_smoke.py` checks the exact
    host momentum-decay and energy-sink map.
  - `tst/scripts/particles/pic_adaptive_deltaf_smoke.py` checks the host fit,
    parser guards, restart fingerprint, and restart parity.
  - `tst/scripts/particles/pic_parser_contract_guards.py` verifies extension
    isolation and bounded-configuration rejection paths.
- Residual closure requirement: complete physical damping comparisons,
  effective-scattering and saturation analysis, transport-coefficient
  calibration, decomposition, MPI, and GPU qualification. The local source
  smokes do not establish calibrated transport coefficients.

### 8. Oblique-shock generality

- Register handling: `parser_reject`, `docs_limit`
- Local disposition: `unsupported_documented_no_selector`
- Local artifacts:
  - `src/pgen/tests/pic_parallel_shock.cpp` is explicitly a parallel-shock
    generator with `B0` parallel to `x`.
  - `inputs/tests/pic_parallel_shock_*.athinput` are explicitly named
    parallel-shock scaffolds.
  - This sidecar records the unsupported disposition.
- Residual closure requirement: no safe obliquity selector or oblique-shock
  parser path exists locally. Do not infer oblique-shock support from the
  parallel generator. If oblique shocks are implemented later, add a
  separately named model selector, parser contract, documentation, and
  qualification campaign.

### 9. Frontier work outside the authorized scheduling and budget boundary

- Register handling: `docs_limit`, `extension_gate`
- Local disposition: `fail_closed_policy_scaffold_present`
- Local artifacts:
  - `tst/publication/readiness/storage_policy.json` records `batch`,
    `debug_preferred_normal_fallback`, and the `10000` node-hour cap.
  - `tst/publication/frontier_control_plane/control_plane.schema.json` limits
    immutable pre-submit manifests to `debug` or `normal`.
  - `tst/publication/frontier_control_plane/submit_frontier_job.sh` reserves
    against the `10000` node-hour cap before `sbatch`.
  - `tst/publication/frontier_control_plane/README.md` prohibits real ledger
    genesis until storage preflight succeeds.
- Residual closure requirement: ledger genesis is complete under the
  user-selected Orion-only bulk-evidence policy. Keep submissions inside the
  recorded scheduling policy and cap; work outside either boundary requires
  explicit revised authorization and user permission.

## Local Q-034 Disposition

This audit closes the local sidecar requirement only. The parser-rejection
evidence is present, proxy rows are visibly classified, oblique-shock
generality is explicitly unsupported, and selected extension rows remain open
by design. Scientific qualification, export review, external review, and
terminal release sign-off remain separate gates.
