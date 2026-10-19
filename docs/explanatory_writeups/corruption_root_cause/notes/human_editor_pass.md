# Human editor pass

## Editorial verdict

This draft mostly sounds like a sharp collaborator at a blackboard. The best
parts are short, specific, and slightly impatient with vague explanations:

> The simulation state was fine. The histogram formulas were mostly fine. The
> temporary coordinate workspace was not.

and:

> The annoying part is that the file writer can still count every cell.

Keep that voice. It is doing real explanatory work.

The weaker passages drift into incident-report language. They summarize the
evidence correctly, but they stop guiding the reader through the mechanism.
The largest problem is the cartoon discussion: the prose says that axis order
matters, but it does not actually walk the reader through how an erased row can
later be repopulated. That missing sentence is important because it explains
the otherwise odd survival of the `output29` radial marginal.

## Highest-priority changes

### 1. Make the cartoon explain axis-order dependence

The current caption says:

> The precise rows that matter depend on PDF axis order.

That is true, but too compressed. A reader can reasonably look at the cartoon,
see radius erased, and then wonder why the text later says the `output29`
radial marginal survives.

Add a paragraph immediately after the cartoon along these lines:

> Read this left to right. Before the angle routine runs, the shared workspace
> already contains coordinates written by earlier routines. `realloc` replaces
> that storage, and the angle routine fills back only its own row. A coordinate
> computed later can still repopulate its row. That is why axis order decides
> which marginals survive.

This is the single most useful prose revision in the document.

### 2. Shorten the source-history paragraph

The two full commit hashes interrupt the explanation just when the reader has
the mechanism in hand. Keep the exact hashes for traceability, but do not make
the reader carry them through a sentence.

Current:

> The source history is unusually clean. Commit
> `ec009a27b87029db6f35a95045ae624ff5c0d1e4` changed...

Suggested body wording:

> The source history gives us an unusually clean before-and-after. Commit
> `ec009a27` introduced the unconditional allocation inside
> `coord_abscostheta`; commit `db953863` removed it. The full hashes are listed
> in the traceability notes.

### 3. Turn the prediction table into a challenge to the explanation

The table is useful, but it arrives in report voice. Introduce it with a more
human sentence:

> A root-cause story is only useful if it gets the exceptions right. This one
> has four places where it could have failed.

That gives the table a purpose beyond inventory.

### 4. Rewrite “What works” around the exceptions

The current bullet list repeats material already covered in the prediction
table. The real rhetorical point is that the explanation gets the mixed
outcome right.

Suggested opening:

> This story earns confidence because it explains the exceptions, not just the
> failures. It predicts the twelve radius collapses, the surviving `output29`
> radial marginal, and the intact `output31` control with the same mechanism.

The existing bullets can follow, but several will then be redundant.

## Figure guidance

### Shared-buffer cartoon

The drawing is useful, but the reader needs explicit instructions about what
the gray “erased / zero” boxes mean and why a later coordinate may survive.
The added axis-order paragraph above should fix this.

The caption can also be more direct:

> Computing `coord_abscostheta` replaces the shared coordinate buffer. Earlier
> rows disappear; only rows computed again afterward can recover. Schematic,
> not evidence.

### `output29` corruption signature

The existing “look at” paragraph is strong. It tells the reader where to look
and why the comparison matters. Tighten the first two sentences so the visual
action comes first:

> Start with the lower-left panel. The archived and rebuilt radial curves sit
> on top of each other. Now move to the lower-right panel: the angular
> marginals separate sharply.

Then give the interpretation. This ordering lets the reader see before being
told what the result means.

The phrase “the source hydro data, radial coordinate, energy-flux weight, and
reduction agree” is accurate but dense. A more spoken version is:

> That agreement tells us the source state and radial accounting survived.
> The failure appears only when angle matters.

### Archived-control errors

This is already the best figure explanation in the draft:

> The red bar is not a failure of the reconstruction. It is the required
> negative control.

Keep that line. Add one concrete visual cue before it:

> The red bar is isolated by roughly eleven orders of magnitude from the green
> controls.

Only use that wording if the plotted scale visibly supports it; otherwise say
“far above the green controls.”

The long generator path in the first figure caption is traceable but visually
heavy. Prefer a short caption and put the full source path in the surrounding
text or evidence-audit note.

## Phrases that sound bureaucratic or over-polished

Current:

> The explanation is specific enough to predict the mixed archive.

Suggested:

> The important thing is that this explanation predicts the mixed result.

Current:

> Those alternatives struggle with the control pattern.

Suggested:

> Those alternatives do not explain why one marginal survives while the joint
> distribution fails.

Current:

> For this archive, the issue is decided well enough to act.

Suggested:

> For the existing archive, we know enough to act.

Current:

> The remaining uncertainty is not the root cause. It is how broadly future
> derived-variable code could repeat the same class of shared-buffer error.

Suggested:

> The root cause is no longer the open question. The open question is whether
> another derived-coordinate routine can make the same mistake.

Current:

> For intuition, the useful schematic is...

Suggested:

> The cartoon explains the mechanism. The paired radial-and-angular comparison
> is the evidence that makes it convincing.

## Rhythm and repetition

- “Remains useful” appears several times. Vary with “survives,” “is still
  trustworthy,” or “is an intact control.”
- “Strong” and “strongest” are used as verdict words. Prefer showing why the
  case is strong, especially near figures.
- Several paragraphs begin with “The.” Most are fine, but the “key figures”
  and “what works” sections would benefit from a sharper spoken opening.
- Keep the short three-sentence paragraph in “The simple picture.” It gives
  the document its best rhythm.

## Recommended pass-two shape

1. Preserve the opening and simple-picture sections almost unchanged.
2. Add the missing axis-order walkthrough after the cartoon.
3. Shorten commit hashes in the body.
4. Let the figures carry more of the case; remove repeated summary language
   from “The key figures” and “What works.”
5. Keep the blunt distinction between what is broken and what remains useful.
