# Human editor pass

## Editorial verdict

This draft has the strongest opening of the three:

> Why not launch all 161 snapshots immediately?

and:

> The parts have passed tests. The assembled production machine has not.

and:

> It repeats the same mistake 161 times.

Those lines make the strategy feel necessary rather than ceremonial. Keep
them.

The middle becomes more administrative. The exact golden-snapshot table and
release checklist are necessary, but the prose around them sounds like a gate
document rather than a collaborator explaining why each condition matters.
The revision should preserve the exact requirements while making the decision
logic more visible.

## Highest-priority changes

### 1. Replace the status line

Current:

> Strongly supported release strategy; the golden Frontier gate has not yet
> run.

Suggested:

> The argument for one golden run is strong. The run itself is still pending.

This is plainer and does not confuse confidence in the strategy with evidence
from the missing test.

### 2. Introduce the golden snapshot as a deliberate choice

Current:

> The selected golden snapshot is fixed:

Suggested:

> The golden case is not an arbitrary sample. It is the exact 8 pc snapshot
> with archived controls and a manageable 128-node footprint:

This tells the reader why the table exists.

### 3. Explain “global AMR closure” once in ordinary language

The phrase appears repeatedly and is central, but a non-specialist reader may
not know what it buys.

Add a plain sentence near the first use:

> In plain language, global AMR closure means that every leaf cell is counted
> once, no region is missing, no region is counted twice, and the accepted
> volume adds back up to the full domain.

After that, the shorthand is fine.

### 4. Give the release checklist a spoken lead-in

Current:

> Release requires all of the following:

Suggested:

> A job exit code is not enough. The golden run passes only if all of the
> following are true:

This makes the checklist feel like a scientific decision rule rather than
process overhead.

### 5. Turn “How this could be wrong” into a real argument with itself

Current:

> The strategy could be too conservative if the full run adds no meaningful
> information over partial tests.

Suggested:

> Maybe this is overkill. If the complete run taught us nothing that the
> eight-shard test had not already shown, then spending 128 nodes on a gate
> would be ceremony.

Then answer directly:

> But the complete run is the first test that can catch missing shards,
> duplicate leaves, overlapping AMR levels, a wrong total domain volume, and
> failures in the 1024-rank reduction.

This is more conversational and makes the defense concrete.

## Figure guidance

### Validation ladder

The current explanation identifies the right conceptual boundary, but it can
point more directly at the picture:

> Everything to the left of the eight-shard step is evidence from a partial
> domain. The 1024-shard step crosses a new line: for the first time, the code
> can check whether the whole AMR volume is present exactly once.

If the figure visually marks completed and pending steps by color, name those
colors:

> The green steps are finished. The first amber step is the one that can still
> expose a campaign-wide failure.

The caption is honest but dry. Suggested wording:

> Completed tests reach eight real shards. The pending golden run extends that
> same path to all 1024 shards. Its cell count is projected from measured shard
> sizes, not from a completed run.

### Production inventory

The existing “look at the yellow part” paragraph is strong because it directs
the eye and explains why the feature matters. Keep it.

The last sentence can be sharper:

Current:

> That makes “run everything and inspect later” especially weak.

Suggested:

> For those 21 snapshots, there is no archived answer to rescue us after the
> fact. The pipeline has to earn trust before they run.

### Cartoon

The flow is visually clear, but “fail: revise, do not multiply” is the best
idea in the drawing and deserves a sentence in the prose:

> The important arrow is the one pointing backward. A failed golden run sends
> us back to the tested Andes path; it does not fan out into production.

## Phrases that sound bureaucratic or over-polished

Current:

> Because the remaining unknowns are common-mode risks.

Suggested:

> Because the remaining failures would hit every snapshot the same way.

Current:

> The validation ladder explains why the golden run is logically different.

Suggested:

> The ladder shows the missing test plainly: none of the completed steps sees
> the whole domain.

Current:

> The green steps strongly support the numerical design.

Suggested:

> The green steps tell us the equations and partial pipeline are behaving.

Current:

> Matching real-data controls sit orders of magnitude inside proposed
> full-validator thresholds.

Suggested:

> The matched real-data controls pass by orders of magnitude.

Current:

> The opposite risk is treating one golden success as universal proof.

This is good and should stay. It is direct, skeptical, and sets up the second
canary clearly.

Current:

> The golden-first strategy is not generic caution.

Suggested:

> The golden run is not caution for its own sake.

Current:

> Prove one exact complete reconstruction, profile it on the real production
> hardware, then multiply.

Suggested:

> First prove one whole snapshot on the hardware that will do the work. Then
> multiply.

## Rhythm and repetition

- “Golden,” “complete,” “run,” and “Frontier” necessarily recur, but several
  sentences use two or three of them together. Prefer “whole snapshot,” “the
  gate,” or “the 1024-shard case” where the reference is unambiguous.
- The exact table and eight-item checklist create a long formal block. Follow
  them with one short sentence: “That is the gate. Anything less is a partial
  success.”
- “Correctness and throughput are both acceptable” is abstract. Name the
  decision: controls pass, global AMR checks pass, and the measured rate is
  affordable for 161 snapshots.
- The “What does not work” section is excellent. Its short negative sentences
  sound candid and should not be polished away.

## Recommended pass-two shape

1. Keep the opening almost unchanged.
2. Explain global AMR closure in one plain sentence.
3. Frame the table and checklist as the reason this exact snapshot is the
   gate.
4. Point more literally at the completed/pending boundary in the ladder.
5. End with the operational rule: one whole snapshot, then one 256-node
   canary, then controlled production.
