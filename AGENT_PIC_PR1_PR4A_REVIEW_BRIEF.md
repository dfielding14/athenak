# AthenaK PIC PR4 Review Brief (PR4a-PR4f)

## 1) Mission
You are a review-focused agent.
Your job is to audit **PR4 only** (PR4a through PR4f) in AthenaK PIC and verify
that the current implementation is numerically sane, test-covered, and aligned
with the documented AthenaK/Entity contract.

This is a review/verification task. Do not implement new PR4 features.

## 2) Start State (Current Working Tree)
- Repo: `/Users/dbf75/Work/Research/AthenaK/athenak-DF`
- Branch: `dev/PIC`
- Expected HEAD: `bb07aa58`
- Working tree is expected to be **dirty** with in-progress PR4f hardening:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_multilevel.py`
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_restart_fidelity.py`
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENT_PIC_HANDOFF.md`
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENT_PIC_IMPLEMENTATION_GUIDE.md`
  - this brief file itself

Run before review:
```bash
cd /Users/dbf75/Work/Research/AthenaK/athenak-DF
git status --short
git rev-parse --abbrev-ref HEAD
git rev-parse --short HEAD
```

If branch/HEAD differs, report that first and continue review against actual state.

## 3) PR4 Scope Under Review
- **PR4a**: runtime selector + compatibility/default guards
- **PR4b**: particle-side dataflow/wrapper plumbing for direct edge current mode
- **PR4c**: direct staggered current deposition kernel (first-order trajectory path)
- **PR4d**: coupling integration into MHD E-field source path
- **PR4e**: continuity/decomposition hardening and guard coverage
- **PR4f**: multilevel/restart direct-vs-conversion closeout checks

Out of scope for this review:
- PR5+ work
- changing compatibility policy (keep to evaluation only)
- higher-order direct deposition implementation

## 4) Mandatory Reading Order
1. `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENTS.md`
2. `/Users/dbf75/Work/Research/AthenaK/athenak-DF/inputs/AGENTS.md`
3. `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/mhd/AGENTS.md`
4. `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/AGENTS.md`
5. `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENT_PIC_HANDOFF.md` (PR4 sections)
6. `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENT_PIC_IMPLEMENTATION_GUIDE.md` (PR4 section)
7. `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.hpp`
8. `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp`
9. `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_tasks.cpp`
10. `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_moments.cpp`
11. `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/mhd/mhd_tasks.cpp`
12. `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_current_coupling.py`
13. `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_decomp.py`
14. `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_multilevel.py`
15. `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_restart_fidelity.py`
16. `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_nonperiodic.py`

## 5) Entity Parity References (Required)
1. `/Users/dbf75/Work/Research/AthenaK/entity/src/engines/srpic.hpp:114`
2. `/Users/dbf75/Work/Research/AthenaK/entity/src/engines/srpic.hpp:521`
3. `/Users/dbf75/Work/Research/AthenaK/entity/src/kernels/currents_deposit.hpp:170`
4. `/Users/dbf75/Work/Research/AthenaK/entity/src/kernels/currents_deposit.hpp:406`
5. `/Users/dbf75/Work/Research/AthenaK/entity/src/kernels/ampere_mink.hpp:142`

Use these to classify each axis as:
- aligned
- intentional divergence
- unplanned divergence

## 6) As-Built PR4 Anchors You Must Verify

### PR4a (selector + guard)
- mode enum/default:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.hpp:40`
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.hpp:111`
- mode parsing:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp:265`
- guard direct->edge representation:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp:384`
- guard direct->strict periodic:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp:394`

### PR4b (dataflow wiring)
- edge-current boundary task IDs:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.hpp:67`
- edge-current boundary object:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.hpp:145`
- pbval_jedge allocation:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp:460`
- old-position staging for direct mode:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_moments.cpp:82`
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_moments.cpp:102`

### PR4c (direct deposition kernel)
- direct kernel launch:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_moments.cpp:263`
- zig-zag style trajectory weighting / atomics:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_moments.cpp:300`
- direct-deposition path gate:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_moments.cpp:249`

### PR4d (coupling integration)
- direct edge-current comm wrappers inserted before `MHD::EFieldSrc`:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_tasks.cpp:198`
- cc-convert path still inserted before `MHD::EFieldSrc`:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_tasks.cpp:267`
- additive edge-current receive accumulation:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_moments.cpp:748`
- MHD edge representation consumer:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/mhd/mhd_tasks.cpp:437`

### PR4e/PR4f (test hardening + closeout)
- direct-mode + continuity checks in coupling test:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_current_coupling.py:392`
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_current_coupling.py:904`
- decomposition direct parity + continuity checks:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_decomp.py:526`
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_decomp.py:575`
- multilevel closeout checks (added):
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_multilevel.py:266`
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_multilevel.py:306`
- restart closeout checks (added):
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_restart_fidelity.py:272`
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_restart_fidelity.py:288`
- per-rank restart negative guard:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_restart_fidelity.py:244`

## 7) Required Review Questions
Answer with exact `file:line` evidence.

1. Is default behavior still compatibility-preserving (`cc_convert`) and opt-in for direct mode?
2. Is direct mode only active under intended runtime/geometry constraints?
3. Does runtime ordering stay coherent: push -> deposit -> sync/comm -> E source usage?
4. Is direct edge-current communication additive and deterministic across decomposition?
5. Are PR4e/PR4f tests meaningful (not only nonzero checks), and are thresholds numerically sane?
6. Any unplanned divergence from Entity on:
   - ordering,
   - trajectory deposition semantics,
   - E += coeff*J application,
   - boundary/comm policy?

## 8) Required Test Matrix
Run exactly:
```bash
cd /Users/dbf75/Work/Research/AthenaK/athenak-DF/tst && python run_tests.py particles/pic_mhd_current_coupling --cmake=-DAthena_ENABLE_MPI=ON --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx
cd /Users/dbf75/Work/Research/AthenaK/athenak-DF/tst && python run_tests.py particles/pic_mhd_coupling_decomp --cmake=-DAthena_ENABLE_MPI=ON --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx
cd /Users/dbf75/Work/Research/AthenaK/athenak-DF/tst && python run_tests.py particles/pic_mhd_restart_fidelity --cmake=-DAthena_ENABLE_MPI=ON --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx
cd /Users/dbf75/Work/Research/AthenaK/athenak-DF/tst && python run_tests.py particles/pic_mhd_coupling_multilevel --cmake=-DAthena_ENABLE_MPI=ON --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx
cd /Users/dbf75/Work/Research/AthenaK/athenak-DF/tst && python run_tests.py particles/pic_mhd_coupling_nonperiodic --cmake=-DAthena_ENABLE_MPI=ON --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx
cd /Users/dbf75/Work/Research/AthenaK/athenak-DF && bash tst/scripts/style/check_athena_cpp_style_changed.sh
cd /Users/dbf75/Work/Research/AthenaK/athenak-DF && bash tst/scripts/style/check_python_style_changed.sh
```

## 9) Deliverable Format (Required)
Return one review report with:

1. `Executive Verdict`
- Healthy / Needs action (one paragraph)

2. `Findings (By Severity)`
- P0/P1/P2 findings first, each with exact `file:line`
- If none, state “No blocking findings”

3. `PR4a→PR4f Status Table`
- Scope intent, as-built status, evidence, risk

4. `Entity Parity Matrix`
- ordering, deposition semantics, coupling semantics, comm/BC semantics

5. `Validation Evidence`
- commands run, pass/fail, and key numeric excerpts

6. `Residual Risks + Next Actions`
- concrete, ordered, scoped to PR4 completion

## 10) Guardrails
1. Do not revert unrelated local changes.
2. Do not use destructive git commands.
3. Do not commit unless explicitly asked.
4. Keep conclusions tied to source and tests, not docs alone.
