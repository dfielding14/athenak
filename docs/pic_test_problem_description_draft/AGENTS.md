# AGENTS.md -- PIC Test Problem Description Draft

## Purpose
This directory contains the publication-style methods manuscript for the AthenaK
MHD-PIC test-problem suite. Edits here must prioritize technical accuracy while
improving readability, narrative flow, and interpretive value.

## Writing Approach (Required)
Use a two-pass editorial workflow on every substantial revision.

### Pass 1: Structure, Flow, and Narrative Architecture
- Reorder sections so the reader moves from scope -> evidence -> interpretation.
- Ensure each section starts with a clear purpose sentence.
- Keep transitions explicit: why this section follows the previous one.
- Group tests by physical role (baseline, PIC core, benchmarks, integration).
- Keep a clear distinction between:
  - what is measured,
  - what is observed now,
  - what remains to be improved.

### Pass 2: Sentence-Level Technical Clarity
Review each section line by line. For every test-task description, include:
- `Logic`: what the test is designed to probe.
- `Results`: what the current run/gate/figure demonstrates.
- `Implications`: what those results mean for code credibility/capability.
- `Next steps` (optional but encouraged): what upgrade closes remaining gaps.

## Style Targets
- Use precise, concrete language; avoid vague claims.
- Keep prose lively and human, but do not sacrifice correctness.
- Prefer short declarative sentences over overloaded clauses.
- Avoid hype language and avoid dry checklist-only prose.
- Tie claims directly to gates, figures, or known diagnostics.

## Accuracy and Claim Discipline
- Do not claim publication-grade physics where diagnostics are still proxy/smoke.
- Flag known limitations (e.g., shock storyline still smoke-level unless upgraded
  diagnostics are present).
- Distinguish decomposition consistency from physical convergence.
- Keep uncertainty visible when inference is provisional.

## Figure-Caption Expectations
Each caption should state:
1. what is plotted,
2. what comparison is being made,
3. what conclusion is supported,
4. what limitation remains (if any).

## Formatting and Build
- Keep this manuscript in AASTeX two-column format.
- Keep lines reasonably wrapped for readability.
- Rebuild after edits and resolve hard TeX errors before handing off.

## Scope Boundary
This directory is for manuscript assets only. Do not modify simulation kernels,
input decks, or test harness logic here unless explicitly requested.
