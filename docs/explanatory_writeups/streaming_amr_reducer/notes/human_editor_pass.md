# Human editor pass

## Editorial verdict

The central explanation is very good. These lines sound exactly like a human
collaborator reducing the problem to its essential shape:

> The desired object is a weighted histogram. It is not a picture of a uniform
> grid. That matters.

and:

> Read one bounded chunk of MeshBlocks. Compute the needed coordinates and
> weights on the GPU. Add them into rank-local histograms. Throw the chunk
> away.

Keep both.

The draft becomes less human when it shifts from explaining the calculation to
defending an “architecture.” Words such as “operational path,” “correctness
case,” “preferred architecture,” and “internal consistency” make the middle
sound like a design review. The underlying argument is simpler: this computes
the quantity we actually want, it fits in memory, and the tests say it is doing
the accounting correctly.

## Highest-priority changes

### 1. Replace the status line

Current:

> Strongly supported preferred architecture; full Frontier performance remains
> untested.

Suggested:

> This is the right calculation. We have not yet shown that the current kernel
> is fast enough on Frontier.

This separates scientific correctness from performance without release-note
language.

### 2. Guide the reader through the memory equation

The memory equation is useful, but the prose immediately becomes compressed.
Tell the reader what part matters:

> There is one line worth noticing in this estimate: memory follows the chunk
> size and histogram catalog, not the size of the full snapshot. That is the
> whole reason the method remains practical.

Then state the measured histogram footprints.

The 4 PB dense-cube estimate is an excellent scale argument. Give it a spoken
lead-in:

> Now compare that with the tempting alternative. A uniform 8 pc cube across
> 400 kpc contains about \(1.25\times10^{14}\) cells...

### 3. Make the predictions sound falsifiable

The prediction bullets are clear, but they read like requirements. Introduce
them with:

> If direct AMR reduction is really the right abstraction, it has to survive
> six checks.

Then retain the list. This turns the section into a test of the idea rather
than a specification.

### 4. Replace the “What works” laundry list with a confidence ladder

The bullets currently mix oracle tests, format compatibility, GPU execution,
real-data runs, and publication behavior without telling the reader how those
pieces build confidence.

Suggested opening:

> The evidence comes in layers. The oracle checks the mathematics. The
> one-rank/two-rank tests check decomposition. The real-shard runs check the
> actual I/O, GPU, MPI, and publication path. The archived controls then check
> the result against surviving external evidence.

After that, retain only bullets that add a concrete fact.

### 5. State the bottleneck without sounding defensive

Current:

> The current all-product kernel is deliberately simple.

Suggested:

> The current kernel is simple enough to audit, and simple enough to
> bottleneck.

Then keep the concrete list of serialized reads, global FP64 atomics, and
host-staged reduction. This sounds candid rather than preemptively defensive.

## Figure guidance

### Streaming data-flow cartoon

The cartoon is clear. The caption should emphasize the one irreversible design
choice:

> Each bounded AMR chunk contributes to persistent histograms and is then
> discarded. At no point is a global dense cube assembled. Schematic, not
> evidence.

After the figure, add a one-sentence reading instruction:

> Follow the blue path left to right. The histograms survive each loop; the
> cell chunk does not.

That sentence would make the bounded-memory mechanism immediately intuitive.

### Andes runtime/throughput figure

The current paragraph makes the right caveat but does not first tell the reader
what visual feature to inspect.

Suggested replacement:

> Look at the move from one K80 to four K80s. The job processes eight times as
> much real data, but aggregate throughput improves by only about
> \(2.09\times\). So the full path works, but it is not scaling cleanly yet.

Then explain what “the full path” includes. This puts observation before
interpretation.

Current:

> The useful result is not that K80 is fast. It is not.

This is lively and worth keeping, but place it after the visual observation.

The source paths make the caption cumbersome. Keep a short human-readable
caption and move exact log paths to the evidence-audit note or a source line
below the figure.

### Same-weight closure figure

“Look at the scale” is a good start, but it needs a visible object:

> Look at how tightly the three same-weight families cluster near zero spread.
> Their relative differences are around \(10^{-13}\).

Then explain why that is an accounting check. “Far below the full-validator
limit” is implementation language and can move to a caveat or audit note.

### Key-figures section

Current:

> Together they answer three different questions: does it run, does it conserve
> its own accounting, and does it reproduce independent surviving evidence?

The three questions are useful. The lead-in is polished and generic. A more
spoken version:

> These figures do different jobs. One says the pipeline runs. One says the
> bookkeeping closes. The archived controls ask the harder question: did we
> reconstruct the right answer?

## Phrases that sound bureaucratic or over-polished

Current:

> This demonstrates the architecture, not production efficiency.

Suggested:

> This tells us the design works end to end. It does not tell us the Frontier
> run will be fast.

Current:

> They are optimization questions inside the correct reduction architecture.

Suggested:

> Those are reasons to tune this reducer, not reasons to build a dense cube.

Current:

> The independent oracle fixes those definitions explicitly.

Suggested:

> The independent oracle pins those definitions down.

Current:

> One complete ... reconstruction on Frontier is the decisive correctness gate.

Suggested:

> One complete ... Frontier reconstruction is the next test that can still
> change our mind.

Current:

> The correctness case is strong.

Suggested:

> The tests say the reducer is doing the right calculation.

Current:

> without scientifically unnecessary resampling

Suggested:

> without inventing a dense grid that the PDF calculation never needed

## Rhythm and repetition

- “Architecture,” “path,” “correctness,” and “complete” recur often. Replace
  them with concrete actions whenever possible: reads, bins, closes, rejects,
  publishes.
- The formal section has several long paragraphs in a row. Break the dense-cube
  comparison into its own short paragraph so it lands as the obvious
  consequence.
- “One complete real shard” and “eight complete real shards” appear in several
  sections. State them fully once; later references can say “the Andes
  real-data runs.”
- Keep the blunt sentence “They are not arguments for a dense cube.” It is one
  of the draft's best pivots.

## Recommended pass-two shape

1. Preserve the simple-picture section.
2. Add one sentence that teaches the reader how to follow the cartoon.
3. Explain the memory equation before giving catalog numbers.
4. Recast the evidence as a ladder from oracle to real archived controls.
5. Make the Frontier gap concrete: right calculation, current kernel speed
   still unknown.
