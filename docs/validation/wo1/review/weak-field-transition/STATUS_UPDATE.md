# Project status update: B4 transport encoding and floor transition

> **Superseded decision, 2026-10-06:** retain conservative A and the current implementation. The user accepted the sharp-contact failure as a known limitation shared with the checked references; the strict expected failure remains. The Q prototype and continuation guidance below are historical investigations, not active proposals or required work. See the [reference comparison and accepted decision](../weak-field-reference/STATUS_UPDATE.md).

- Date: 2026-10-01
- Exact timestamp: 2026-10-01T16:34:08Z
- Design-review addendum: 2026-10-01T16:54:32Z; Q remains an experimental candidate, not a selected production replacement.
- Project or repository: `/Users/dbf75/.codex/worktrees/bf22/athenak-DF`
- Report profile: compact; focused diagnosis and formulation experiment, not a production-readiness report.
- Status: investigation complete; prototype validated in the cases below; production B4 fix remains incomplete.
- Branch: `c/cgl-lf-wo1`
- Commit: `7b3345fd262a1cf826e9476fb40ca20a28a15881`
- Worktree state: production source unchanged by this investigation; this evidence directory is new. Root work orders `WO1_CGL_LF_fixes.md` and `WO2_CGL_LF_followups.md` remain untracked.
- Agent identifier: `/root`, with independent analysis by `/root/b4_discrete_analysis` and `/root/b4_reset_audit`.
- Data or simulations analyzed: first-five-cycle full-application traces, a 50-cycle reference-field counterfactual, 24 prototype contacts, and nine prototype compression runs.
- Compute environment: local Darwin arm64; Apple Clang 21.0.0; Release, double precision, Kokkos Serial; MPI/OpenMP/CUDA off. Kokkos commit `08ceff92bcf3a828844480bc1e6137eb74028517`.

# Tier 0: What happened and why it matters

The task was to identify how the first five updates generate anisotropy and choose a consistent encoding and weak-field transition. The failing magnetized cell was traced through face fluxes, the first-order flux correction (FOFC) predictor, Runge–Kutta (RK) updates, constrained transport, EOS conversion, and the end-of-cycle wall projection.

The failure is already present before EOS repair. The first cell exceeding the test's pressure-ratio bound has $B=0.08547$ at cycle 5, far above `bfloor=1e-10`. It never undergoes an EOS reset during these five cycles. Changing only the sub-floor isotropy reference to a finite field of 1 still produces a ratio of $2.55\times10^{10}$ after 50 cycles.

The current conserved quantity contains $-3\ln B$. Mixing it linearly while mixing $B$ linearly does not preserve isotropy. The actual RK blend produces this discrepancy even when both constituent fields are magnetized. Floor resets are present upstream, but they are not necessary for the demonstrated failure.

The experimental formulation evolves $Q=\rho\ln(p_\perp/p_\parallel)$ during hyperbolic RK, including its required CGL velocity-strain source. Below the field floor, use $Q=0$ at fixed internal energy and disable the entire anisotropy source. When the field rises above the floor, retain the evolved $Q$ and resume that source. No neighbour or floor field enters the encoding. The subsequent design review below withdraws the initial recommendation to adopt this formulation: the evidence establishes a discrete defect, but not that replacing A is necessary or preferable.

A scratch implementation passes all 24 contact cases for 50 cycles, with pressure ratios between 0.99410 and 1.02520. Nine smooth-compression tests verify the source sign, coefficient, approximately second-order convergence, and isotropic sub-floor limit. This supports further investigation; it does not validate shocks, multidimensional gradients, actual FOFC fallback, AMR, LF coupling, nonideal induction, or restart compatibility. A matched turbulence comparison is a primary adoption gate. The existing strict expected failures remain unchanged.

# Tier 2: Detailed methods, implementation, and validation

## Trace and counterfactuals

The input is `inputs/unit_tests/cgl_weak_field_transport.athinput`: 128 cells, density 1, uniform velocity $v_x=10$, transverse field $B_y=10^{-12}$ on the left and 1 on the right, and isotropic gas pressures 1.5 and 1. Total pressure is initially continuous. The method is RK2, CFL 0.4, PLM/HLLE, with FOFC enabled. Optional collisions and kinetic limiters are disabled; the mandatory fluid-firehose wall remains active.

The original code stores

$$
s=\frac A\rho=\ln\frac{p_\perp}{p_\parallel}+2\ln\rho-3\ln B.
$$

Sub-floor isotropy uses $B_{\rm floor}$, giving $s=69.07755$ for this input. HLLE's anisotropy flux is $F_A=F_\rho s_{\rm upwind}$; it is not a separate standard HLLE dissipative jump in cell-centred $A$. The existing B4 face rule substitutes the other reconstructed state's field when the upwind state is sub-floor. Cell storage still uses the floor.

The first failing cell is zero-based cell 64, centred at $x=0.50390625$ (internal index 67). Values below are after the full cycle, including the mandatory wall:

| Cycle | Field magnitude | $p_\perp/p_\parallel$ | Wall change in $A$ |
| --- | ---: | ---: | ---: |
| 1 | 0.67887020 | 0.69567330 | +0.22258976 |
| 2 | 0.42472229 | 0.87475477 | +0.29741703 |
| 3 | 0.25363780 | 1.12115764 | 0 |
| 4 | 0.14788253 | 1.78092077 | 0 |
| 5 | 0.08547084 | 3.00522003 | 0 |

EOS conversion changes none of this cell's density, energy, or $A$ in these cycles. Sub-floor EOS resets occur in upstream cell 63. The wall raises $A$ in cycles 1 and 2; the later excess exists before that wall. Before/after-FOFC face fluxes are identical, so fallback is not responsible. The reversed-flow trace mirrors the result in cell 63. Physical pressures are evaluated after the magnetic update; combining post-RK fluid variables with pre-CT fields would give a misleading intermediate ratio.

Three controls discriminate the proposed causes:

- Raising `bfloor` from $10^{-10}$ to $10^{-6}$ leaves all physical table outputs through cycle 5 byte-identical.
- Lowering it to $10^{-14}$ makes the input magnetized everywhere. The ratio then reaches $5.5232\times10^9$ in cycle 1 without any EOS reset in that cycle. This uses the ordinary magnetized encoding at both face states.
- Encoding sub-floor cell isotropy with a fixed reference field of 1 changes the initial weak-cell $s$ from 69.07755 to 0. Physical outputs remain byte-identical through cycle 40. They first differ at cycle 41; at cycle 50 the maximum ratio is $2.5484\times10^{10}$ versus $2.4630\times10^{10}$ in the original. This is a failed candidate, retained in [bounded-reference.patch](bounded-reference.patch).

## Why linear mixing fails

For two isotropic states with equal density and aligned positive fields, linear mixing of $s$ and $B$ with fraction $f$ gives

$$
\frac{p_\perp}{p_\parallel}\bigg|_{\rm mix}
=\left[\frac{(1-f)B_1+fB_2}{B_1^{1-f}B_2^f}\right]^3\ge1.
$$

This is the arithmetic/geometric-mean mismatch, with no floor required. Numerical spatial mixing and the RK convex combination both need a compatible treatment.

The actual cycle-5 RK2 combination demonstrates the problem directly. Reconstruct its second Euler component algebraically as $U_E=2U_{\rm final}-U_{\rm start}$ and $B_E=2B_{\rm final}-B_{\rm start}$, before the final wall. The initial and Euler-component physical log ratios are 0.57713052 and -0.66116873. Their density-weighted mean is -0.04189505, but the code's blend yields 1.10035079, exceeding both. The excess is 1.14224584, with $B_E=0.02305915>B_{\rm floor}$. This isolates a discrete incompatibility in the existing representation, rather than assigning the entire error to a particular face flux. All five budgets are in [rk-mixing-budget.json](rk-mixing-budget.json).

## Experimental encoding and transition

Use $q=\ln(p_\perp/p_\parallel)$ and $Q=\rho q$. For smooth ideal CGL,

$$
\frac{D\ln p_\perp}{Dt}=\hat b_i\hat b_j\partial_jv_i-2\nabla\cdot\boldsymbol v,
\qquad
\frac{D\ln p_\parallel}{Dt}=-2\hat b_i\hat b_j\partial_jv_i-\nabla\cdot\boldsymbol v,
$$

so the required conservative transport equation with a strain source is

$$
\partial_tQ+\nabla\cdot(Q\boldsymbol v)
=\rho\left(3\hat b_i\hat b_j\partial_jv_i-\nabla\cdot\boldsymbol v\right).
$$

Here $\hat{\boldsymbol b}$ is the magnetic unit vector. The scratch face flux is $F_Q=F_\rho q_{\rm upwind}$. At equal $q$, density-weighted mixing preserves that $q$ independently of $B$; the strain source supplies physical anisotropy.

1. For $B\le B_{\rm floor}$, pressure conversion sets $p_\parallel=p_\perp=2e_{\rm int}/3$, $Q=0$, and the whole source is zero. This retains the intended isotropic weak-field fallback at fixed internal energy.
2. For $B>B_{\rm floor}$, decode the evolved $Q/\rho$ directly and enable the CGL source. A cell becoming magnetized retains its transported/updated $Q$; it receives no reference-field offset or additional isotropy reset.
3. Use the same representation in every RK register, reconstructed face, pressure conversion, and FOFC prediction during the operator. Floors and wall projections must re-encode in that same representation.

The scratch binary changes the IAN slot globally to $Q$, uses centred one-dimensional velocity gradients evaluated from pre-stage primitives, and adds the same stage-weighted source to the actual and predicted updates. It rejects multidimensional or LF-split RK evolution. Total-energy and magnetic updates are unchanged. This is design evidence only: $q$ has no magnetic-reference divergence, but can still diverge as a pressure tends to zero.

## Validation and limits

| Purpose and method | Expected result / tolerance | Actual result | Status | Caveat |
| --- | --- | --- | --- | --- |
| Compare instrumented and original physical outputs | Exact bytes | Identical | passed | Trace-enabled CPU build |
| Original and bounded-reference 50-cycle contacts | Positive pressures and ratio in [0.5, 2] | Maximum ratios $2.46\times10^{10}$ and $2.55\times10^{10}$ | failed | Demonstrates inadequate reference-only fix |
| Prototype: dc/plm × velocities +10, -10, +1, 0 × floors $10^{-14},10^{-10},10^{-6}$ | All 24 cases, cycles 0–50, positive pressures and ratio in [0.5, 2] | Range [0.99410, 1.02520]; supersonic PLM [0.99630, 1.00515] | passed | 1D transverse contacts |
| Open-boundary energy budget in the same 24 runs | Absolute residual below $10^{-13}$ | At most $7.11\times10^{-15}$; boundary primitives unchanged | passed | Budget includes boundary energy flux |
| Parallel/perpendicular affine compression, nx 128/256/512 | $q$ error below $10^{-4}$, decreasing with resolution | Errors below $1.12\times10^{-6}$, approximately fourfold reduction on doubling nx | passed | Smooth source, no shocks |
| Sub-floor affine compression at all three resolutions | $q=0$ | Exactly zero | passed | Does not test anisotropic threshold crossings |
| Prototype multidimensional, shock, forced FOFC, AMR, LF, restart and GPU/MPI checks | Consistent production behaviour | Unavailable | not run | Required before production acceptance |

The contact energy budget is $E_f=E_0+0.25|v_x|t$ on the unit domain. All 33 prototype runs retain positive pressures without recorded EOS floor or FOFC events. Compression uses $v_x=-0.2x$ on $[-50,50]$ through $t=0.2$, measuring $|x|<10$; magnetized fields have magnitude 10. With $\rho=1/0.96$, the exact values are $q_\parallel=2\ln(0.96)$ and $q_\perp=-\ln(0.96)$.

| Compression | Maximum $q$ error, nx 128 | nx 256 | nx 512 |
| --- | ---: | ---: | ---: |
| Parallel | 1.1190e-6 | 2.8275e-7 | 7.1566e-8 |
| Perpendicular | 5.5633e-7 | 1.4004e-7 | 3.5279e-8 |
| Sub-floor | 0 | 0 | 0 |

The source equation is equivalent to smooth ideal CGL, but a centred nonconservative strain source is not proof of equivalent shock weak solutions or heating partition. Nonideal induction also needs an explicit policy: preserving the old advected-$A$ model adds $3\rho\hat{\boldsymbol b}\cdot(\partial_t\boldsymbol B)_{\rm nonideal}/B$ to the $Q$ equation. Neither issue was settled by these experiments.

## Subsequent design review: turbulence, discrete sources, and walls

The user's review identifies the central missing test: [Squire et al. (2023), section 2 and Appendix A.3](https://doi.org/10.1017/S0022377823000727), motivates the conservative formulation by errors associated with explicit parallel strain and describes grid-scale compressive modes with $k_\perp\gg k_\parallel$. Successful smooth one-dimensional tests do not address this failure mode. A matched A/Q experiment must compare pressure-anisotropy PDFs, threshold volume fractions, and spectra, including directional high-wavenumber power. Use identical physical inputs, forcing realization, LF/limiter rules, resolutions, and time windows; inspect pre-projection anisotropy and wall activity so a wall cannot conceal noise. Agreement should be assessed against sampling uncertainty and resolution dependence. A difference is not automatically evidence against Q because A also has a measured mixing bias.

The requested field-based source is correct as a smooth identity, but fixed-cell CT differences contain advection and are not material changes. Moreover, exact discrete equivalence to the existing A update imports its bias. Define $h=3\ln B-2\ln\rho$ and $H=\rho h$, so $Q=A+H$. With matched reconstructed states and $F_H=F_\rho h_{\rm upwind}$, the Euler-stage source increment that exactly reproduces A is

$$
\Delta Q_{\rm src}=H^{n+1}-H^n+\Delta t\,\nabla_h\cdot F_H.
$$

Since $F_Q=F_A+F_H$, subtracting $H^{n+1}$ leaves exactly the original conservative A update. For RK, replace $H^n$ by $\gamma_0H^{\rm current}+\gamma_1H^{\rm saved}$ and $\Delta t$ by $\beta\Delta t$. The nonlinear H change from RK mixing is therefore restored too. This derivation concerns magnetized states with matched flux construction; the production sub-floor face substitution requires separate treatment.

A direct scalar check mixed equal-density isotropic states with aligned $B_1=0.085$, $B_2=1$, and fraction $f=0.3$. Both A and this exactly equivalent Q update give $q=2.10759514$, or $p_\perp/p_\parallel=8.22842927$. A field-based material source can differ: taking the logarithm after transporting the fields differs from transporting their logarithms. It may preserve this contact, but is not automatically the same A shock discretization. A complete staggered-CT material-remap scheme has not been specified or tested.

Physical nonideal induction and numerical mixing also require separate interpretation. Writing the extra induction term as $\boldsymbol N$ gives the algebraic contribution $3\rho\hat{\boldsymbol b}\cdot\boldsymbol N/B$. Including it reproduces the field dependence of advected A, but does not by itself prove magnetic-moment conservation for a dissipative model. At fixed density, preserving both double-adiabatic invariants requires $p_\perp\propto B$, $p_\parallel\propto B^{-2}$, and $\dot e_{\rm int}=(p_\perp-p_\parallel)\dot B/B$. An independent Joule-heating contribution requires a pressure-heating partition. Conserving total E does not alone guarantee identical thermalization or pressure partition in different evolving solutions.

The current code has no unconditional fluid-mirror wall. The pure-CGL fluid boundary in code units is $p_\perp^2\le3p_\parallel(B^2+2p_\perp)$, distinct from the kinetic mirror threshold; see [Bhoriya et al. (2024), equations 18–20](https://arxiv.org/html/2405.17487v1#S4.SS1). Adding it would protect hyperbolicity but would not fix the cycle-5 contact error: the measured ratio 3.0052 is below its fluid-mirror limit of approximately 6.0118. Both fluid walls also permit $q\to-\infty$ at low beta: $B^2=1$, $p_\parallel=0.1$, $p_\perp=\epsilon\to0$ satisfies both. Positive pressure floors remain necessary.

The 3D LF turbulence comparison is **not run**. The existing Q prototype explicitly rejects 3D and LF. Current production input definitions include `inputs/cgl_lf_paper/cgl_lf_paper_standard_active_alfvenic_beta10.athinput` and its beta-100/random-forcing counterparts; `scripts/analyze_cgl_lf_paper.py` already computes the main PDFs, occupancy, and spectra. The $8\times8\times16$, $t=0.01$ smoke case cannot establish turbulent statistics. Which Q source to extend is a scientific choice: the centred-strain prototype tests the stated strain-noise concern directly, whereas a newly designed field-based source is a different candidate.

# Tier 3: Reproducibility, audit trail, and handoff

## Repository state and retained outputs

All experiments start from the commit in the metadata. Production source and the two strict B4 expected failures remain unchanged. Independent review checked the source sign, stage weighting, EOS/reset paths, and migration risks. No production merge acceptance is implied.

| Artifact | Contents and use |
| --- | --- |
| [provenance.json](provenance.json) | Source/Kokkos commits, original scratch paths, executable SHA-256 hashes, exact-output and no-FOFC checks |
| [instrumentation.patch](instrumentation.patch), [instrument.py](instrument.py) | Full trace instrumentation, including the added header; historical construction script |
| [run_trace.py](run_trace.py), [analyze_trace.py](analyze_trace.py), [trace-commands.json](trace-commands.json), [trace-summary.json](trace-summary.json) | Exact five-case command arguments and per-cycle reset/wall summaries |
| [plm-state.csv.gz](plm-state.csv.gz), [plm-stage-log.txt.gz](plm-stage-log.txt.gz) | Compressed full active-cell phase trace and face/EOS diagnostic log for the primary PLM case |
| [rk-mixing-budget.json](rk-mixing-budget.json), [bounded-reference-result.json](bounded-reference-result.json), [bounded-reference.patch](bounded-reference.patch) | RK decomposition and failed reference-field candidate |
| [q-prototype.patch](q-prototype.patch), [q-prototype-README.md](q-prototype-README.md) | Complete scratch candidate, including compression fixture, with method and limitations |
| [run_q_contacts.py](run_q_contacts.py), [run_q_compression.py](run_q_compression.py), [q-contact.athinput](q-contact.athinput), [q-compression.athinput](q-compression.athinput) | Runnable experiment scripts and input bytes; archive names differ from scratch names |
| [q-contact-results.json.gz](q-contact-results.json.gz), [q-compression-results.json](q-compression-results.json) | All 24 contact cases' per-cycle extrema and nine compression results |
| [apply_q_prototype.py](apply_q_prototype.py), [add_q_compression_fixture.py](add_q_compression_fixture.py) | Historical patch-construction scripts; prefer the complete patch for replay |

The 22 evidence files total approximately 166 kB before this report. The repository's `docs/validation/wo1/SHA256SUMS.json` covers this archive. Raw outputs and build logs remain in `/tmp/cgl-b4-transition-20261001/` and `/tmp/cgl-b4-q-prototype-20261001/`; not every raw table is copied here. The trace CSV and summarized prototype results are retained independently of those temporary directories. No figures were needed for this narrow comparison.

## Commands, environment, and replay

Original trace commands are in `trace-commands.json`, run using `B4_TRACE=1` and the preserved `athena-trace-baseline`. The trace scratch `build/src/athena` now contains the bounded-reference candidate, so it must not be mistaken for the baseline executable. Instrumentation and counterfactual patches are separate; the prototype patch is also based directly on the recorded production commit.

Original prototype build and run commands, from `/tmp/cgl-b4-q-prototype-20261001`, were:

```sh
cmake -S source -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j8
/Users/dbf75/.uv/envs/interactive/.venv/bin/python3 run_prototype.py
/Users/dbf75/.uv/envs/interactive/.venv/bin/python3 run_affine.py
```

For replay, extract a fresh `git archive` of the recorded commit into an isolated `source/`, supply the recorded Kokkos checkout, and apply `q-prototype.patch` there. Copy `run_q_contacts.py` to the experiment root as `run_prototype.py`, `run_q_compression.py` as `run_affine.py`, and `q-contact.athinput` as `input.athinput`. Build and run the commands above using a Python environment with NumPy. The compression script generates `affine.athinput`. Historical trace/construction scripts contain absolute paths; adapt only a working copy. The archived patch application was checked against a fresh export; a fresh full rebuild from this renamed archive was not run.

Compute was on one local machine, with eight parallel build jobs and serial simulations. No scheduler, GPU allocation, or MPI run was used. Aggregate wall time and node-hours were not recorded; there is no performance claim.

## Continuation requirements

Retain production A while comparing candidate discrete updates. The contact and smooth-source experiments do not establish Q as the preferred production formulation. If Q is pursued, decide explicitly whether production keeps the existing A representation outside that operator or migrates globally; the scratch prototype only tests the latter locally. Preserving external A would reduce restart/LF format changes, but it still requires all RK registers and boundary/AMR operations inside the operator to agree on Q.

Before applying a production change, cover multidimensional strain and oblique shear, shocks, actual FOFC fallback, and initially anisotropic cells crossing the threshold in both directions. Audit EOS, RK baselines, LF $A\leftrightarrow\mu$ conversion, AMR restriction/reflux/prolongation, boundary conditions, direct problem-generator encodings, diagnostics, and nonideal terms. Global $Q$ must convert or reject old $A$ restarts; silently reinterpreting their IAN slot is invalid. Keep the existing 50-cycle acceptance bounds and strict expected failures until the real production implementation passes them.
