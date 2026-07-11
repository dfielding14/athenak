# AthenaK MHD-PIC Project Ethos

## Purpose

We are building MHD-PIC in AthenaK to do science: turbulent plasmas with kinetic
cosmic rays, non-relativistic shocks, and self-consistent Bell amplification. The code
must be trustworthy, but trust is not measured by the number of guards, test permutations,
or process documents we create. It comes from a coherent physical model, clean code,
decisive validation, and successful use on the problems we care about.

We aim for the productive middle ground between carelessness and paralysis. We do not
accept physically inconsistent shortcuts. We also do not spend months protecting against
imagined failures in configurations we do not intend to run.

## The non-negotiable standard

The implementation must be physically correct and internally consistent within its stated
model and regime of validity.

- Start from an explicit equation set and state its assumptions.
- Carry signs, units, normalizations, time centering, and staggering consistently into the
  code.
- Treat coupled terms as one closure. A switch must not leave an inconsistent subset of
  induction, particle forcing, gas feedback, or energy exchange.
- Preserve momentum, energy, charge/current normalization, particle identity, and magnetic
  divergence to the degree promised by the numerical method.
- Never hide a physical failure with a floor, silent correction, stale diagnostic, or
  misleading label.
- State real limitations plainly. A controlled approximation is acceptable; an
  undocumented change of model is not.

Physical correctness does not mean pretending that a reduced model contains every plasma
effect. It means implementing the chosen approximation faithfully and using it where its
assumptions are defensible.

## How we write code

Changes should be surgical: the smallest coherent modification that solves the present
problem and fits AthenaK's design.

- Read nearby code first and match its naming, control flow, data layout, Kokkos patterns,
  error handling, and formatting.
- Prefer boring, direct code over clever abstractions.
- Reuse established storage and task paths unless they cause the problem.
- Avoid parallel implementations of the same physical operation.
- Keep diffs narrow enough that their physical and numerical consequences are easy to
  review.
- Do not mix opportunistic cleanup into a physics fix unless correctness requires it.

Comments should explain information the code cannot express clearly: a physical sign,
normalization, frame, time level, conservation argument, or non-obvious constraint. Do not
add comments that merely narrate the next line.

## Robustness without defensive clutter

Protect failures that could silently corrupt science, hang a run, or waste a large
allocation. Do not scatter checks through every loop for states that cannot arise in a
supported workflow.

Validate important inputs at construction or injection boundaries. Reject inadmissible
states before they are silently repaired. Enforce timestep limits associated with forces
actually present. Preserve ownership and identity across migration and restart. Put each
check at the narrowest shared boundary that covers the real failure.

When a real bug appears, fix it and add the smallest regression that would have caught it.
Robustness should grow from evidence, not imagination.

## Validation proportional to the claim

Validation should answer a question rather than fill a matrix. The normal sequence is:

1. A local algebraic, manufactured, or invariant test that catches a wrong sign, factor,
   normalization, centering choice, or exchange term.
2. One physical benchmark with a known orbit, conservation result, dispersion relation,
   jump condition, or growth rate.
3. A small production-shaped pilot on the architecture and decomposition that matter.
4. A targeted sensitivity test only when the claim depends on it or the pilot reveals a
   sensitivity.

Do not take the Cartesian product of precision, dimensionality, resolution, particle
count, MPI decomposition, restart, CPU, GPU, and every runtime option. Choose
representative coverage and expand it when a result or planned configuration gives a
reason.

Exploratory runs and plots should happen early. They become evidence for a claim after the
relevant benchmark and sensitivity checks pass.

## Full physics first, switches for comparison

The default science path should implement the complete selected MHD-PIC closure. Optional
modes are for controlled comparisons, reproduction, and term-isolation tests. They should
be few, explicit, and physically coherent.

Prefer one atomic model choice over switches for every connected subterm. An experimental
coefficient is not a physical model unless it represents a defined physical parameter.

## Efficiency

Efficiency is part of scientific capability, but optimization should be evidence-driven.
Avoid obvious repeated allocation, global communication, or host-device movement in a hot
path. Profile a production-shaped run before a substantial rewrite, optimize the measured
bottleneck, and rerun one representative physics case afterward.

Clear code with adequate performance is better than complicated code optimized for a
hypothetical workload.

## What we fix now

Before starting a safeguard, refactor, or validation campaign, ask:

1. Could the issue silently change the physics or invalidate a planned result?
2. Can it occur in a supported configuration we intend to run soon?
3. Has it been observed, demonstrated in the code, or established by a direct scaling
   argument?
4. What is the smallest implementation and test that resolves the actual risk?

A demonstrated physical or conservation defect in a planned workflow is immediate. So is
a known crash, hang, corruption path, or severe measured bottleneck. A limitation outside
the supported workflow should be documented and deferred. A plausible but unobserved issue
belongs in a short backlog. A speculative issue with no use case can be left alone.

## Tracking work

Trackers should make the next decision easier, not preserve every possibility considered
during a review.

- Separate active correctness work, known limitations, and optional backlog items.
- Collapse closed findings to the fix commit and decisive test.
- Do not make deferred AMR, every platform combination, speculative optimization, or
  unrelated whole-code hardening a gate for uniform-grid science.
- Remove obsolete acceptance policies when the physical target changes.
- Delete temporary trackers when active items are closed or transferred.

## When a change is done

A change is done when it matches the intended equations and AthenaK style, is clear and
surgical, passes the focused regression and physical benchmark, has representative
coverage for the intended execution path, and has no known issue that invalidates its next
use.

Then we use the code. If it breaks in a new way, we learn from the failure, make the
smallest correct fix, and continue.
