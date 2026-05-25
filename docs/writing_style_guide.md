Across both papers, the signature style is lucid, theory-forward, and physically anchored. The author is not trying to sound ornate. He is trying to make a complicated mechanism feel intelligible, constrained, and worth taking seriously. The prose works because it repeatedly does four things: it identifies the central uncertainty, reduces the problem to a small set of controlling parameters, derives regime-dependent consequences, and then checks or calibrates those consequences against simulations or observations.

Here is a writing style guide an AI agent could use to emulate that voice.

## 1. The governing aesthetic

Write like a patient expert doing live theory in front of the reader.

The tone is calm, exact, and slightly pedagogical. It is confident without swagger, speculative without vagueness, and technical without trying to dazzle. The prose assumes the reader is smart but busy, so it constantly answers three questions:

1. What is the actual problem?
2. What controls it?
3. What follows, physically, from those controls?

The style is not primarily “beautiful” at the sentence level. Its beauty comes from architecture: each section narrows the uncertainty, each equation earns its place, and each conclusion is tied back to either a scaling argument, a simulation, or an observational constraint.

The deepest rule is this: **imitate the reasoning structure more than the vocabulary.** The surface phrases matter, but the real signature is the sequence
**stakes → uncertainty → simplified model → key dimensionless parameters → analytic estimate → validation/calibration → implications → caveats**.

---

## 2. High-level organization

### Abstract

The abstract follows a very consistent pattern:

* Open with a broad unresolved problem or uncertainty.
* State exactly what the paper does.
* Name the method in one compact phrase.
* Give the core bifurcation or regime split.
* State the observational or physical implication.
* End with broader significance.

The rhythm is usually:
“X remains uncertain. In this paper, we do Y. We use A validated by B to derive C. For regime 1, result 1; by contrast, for regime 2, result 2. Given observational constraints, implication 3. We discuss implications for Z.”

This style favors **compact, high-density abstracts** with very little throat-clearing.

### Introduction

The introductions in both papers are exemplary. They proceed in layers:

1. Start with the large astrophysical stakes.
2. Name the persistent puzzles or tensions.
3. Introduce the candidate mechanism of interest.
4. Summarize prior work by category, not by citation dump.
5. Identify the central physical uncertainty.
6. Explain why that uncertainty matters for modeling.
7. Bring in an observational handle if one exists.
8. State the paper’s scope and roadmap.

This creates a feeling that the paper is entering an existing conversation at exactly the right point.

For AI use: never open with generic praise for the field. Open with a meaningful unresolved issue. Then narrow quickly.

### Core body structure

The first paper uses a sequence like:

* analytic approximations
* numerical validation
* observational calibration
* discussion

The second uses:

* energetic plausibility
* spatial distribution
* transport uncertainty
* observational prospects
* discussion

That tells you the preferred organizational principle: **build from simplest argument to strongest, most reality-connected argument**.

A section should usually do one of four jobs:

* isolate a mechanism
* derive a scaling
* validate that scaling
* translate it into observable or broader consequences

### Discussion and summary

The discussion does not merely repeat results. It does four distinct things:

* restates the physical picture in plain language
* clarifies where the argument is robust and where it depends on uncertain microphysics
* compares with prior work
* identifies the next observational or numerical tests

The summaries zoom out without losing technical fidelity.

---

## 3. The conceptual habits to imitate

### A. Start from uncertainty, not from novelty

The author almost always frames the paper around a real uncertainty:
transport physics, diffusion coefficient, the link between feedback and observables, the spatial distribution of CR pressure.

That is different from saying “we present a novel model.” The style is less self-promotional and more problem-oriented.

Use formulations like:

* “X remains a key uncertainty in assessing whether Y.”
* “It is not clear whether…”
* “A central question is…”
* “To assess the implications of X, we must have a handle on Y.”

### B. Reduce the problem aggressively, but transparently

The papers repeatedly say, in effect: here is the simplest model that still captures the controlling physics.

That move is never hidden. It is announced.

Use sentences like:

* “A key simplification can be made…”
* “With these approximations…”
* “We consider a simple model in which…”
* “To provide an order of magnitude estimate…”

The simplification is framed as a tool, not an apology.

### C. Organize around regimes

A huge part of the style is regime logic.

The author is always asking:

* what is the controlling inequality?
* what happens above it?
* what happens below it?

That produces clean contrasts:
rapid diffusion vs slow diffusion, group vs cluster, early feedback vs maintenance feedback, calorimetric vs escape-dominated, dynamically important vs negligible.

Make these contrasts explicit and early:
“For X larger than Y… By contrast, for X smaller than Y…”

### D. Translate math into physical meaning immediately

No equation is left uninterpreted for long. After deriving something, the prose quickly answers:
What does this mean physically?
Which terms matter?
What scaling is most important?
What does it predict for a Milky Way-like or group-mass system?

This is essential to the voice.

### E. Re-express results in intuitive or observable forms

The author frequently derives a quantity one way and then rewrites it in terms of:

* star formation rate
* hydrostatic pressure
* halo mass
* gamma-ray luminosity
* momentum flux relative to radiation
* dimensionless pressure ratios

This move is central to the style. It says: the derivation matters, but interpretation matters more.

### F. Make the strongest claim justified by the assumptions

The prose is careful, but not timid.

It uses calibrated confidence:

* “shows” for direct derivations or numerical results
* “suggests” when the inference is one step removed
* “likely” when the evidence is good but not definitive
* “may” or “we speculate” for broader implications

This is one of the author’s strongest habits.

---

## 4. Section-level playbook

### Section opening

A section usually begins by telling the reader exactly what it is for.

Examples of the move, paraphrased:

* “In this section we focus on analytic approximations…”
* “To assess the implications of these results, we must…”
* “Here we present numerical models with…”
* “We now derive simple estimates to elucidate…”

A good section opener has three ingredients:
purpose, scope, and why it matters.

### Middle of section

The middle usually alternates between:

* one conceptual paragraph
* one equation or estimate
* one interpretation paragraph
* one concrete example or figure reference

That alternation makes dense material readable.

### Section closing

Sections often end by restating the main takeaway before moving on.

Use endings like:

* “A key conclusion from this section is…”
* “Taken together, these estimates imply…”
* “This motivates the numerical calculations that follow.”
* “We return to this point below.”

---

## 5. Paragraph architecture

The standard paragraph template is:

1. **Topic sentence:** make the claim or define the task.
2. **Development:** give the estimate, equation, or reasoning.
3. **Interpretation:** explain what the result means physically.
4. **Closure:** state the takeaway or transition.

A paragraph should do one job. Avoid paragraphs that simultaneously introduce a new model, summarize literature, derive an equation, and make three implications.

The style is modular.

A very typical paragraph in this voice looks like:

“To assess whether X can drive Y, we compare timescale A to timescale B. In the limit A ≪ B, the governing equation simplifies to…. This shows that the system is controlled primarily by parameter Z. Physically, this is because…. For fiducial parameters appropriate to…, the result is…. This approximation is valid so long as….”

That is almost the whole style in miniature.

---

## 6. Sentence-level rules

### Voice

Prefer active first-person plural:
“We show,” “we find,” “we estimate,” “we argue,” “we model,” “we return to…”

This voice feels collaborative and authoritative without being self-important.

### Rhythm

Use medium-to-long analytical sentences, but break them with short, clean payoff sentences.

Pattern:
long explanatory sentence → short interpretive sentence.

For example:
“This scaling implies that the diffusion time decreases rapidly with increasing compactness of the system, so the CR pressure cannot build to as large a value.”
“Thus the wind is weak.”

That short final sentence gives clarity and force.

### Syntax

The prose often opens with logical framing rather than the grammatical subject:

* “By contrast, …”
* “Given these assumptions, …”
* “For reference, …”
* “In the limit…, …”
* “Physically, …”
* “Together with…, …”
* “More generally, …”

This is a major stylistic trait. It makes the argument feel chained together.

### Contrast structure

Use contrast relentlessly and precisely:

* “For…, … . By contrast, for…, … .”
* “The former…, while the latter…”
* “This is true in X, but not in Y.”

### Definition discipline

Define terms once, then reuse the same phrase. Do not hunt for synonyms just to avoid repetition.

This author repeats technical nouns on purpose:
diffusion coefficient, CR pressure, scale-height, sonic point, group-mass halo, outer halo, momentum flux.

Clarity beats variety.

---

## 7. Word choice and diction

### Preferred verbs

Use verbs of inference and mechanism, not branding:

* derive
* estimate
* quantify
* assess
* calibrate
* constrain
* validate
* reproduce
* imply
* suppress
* dominate
* accelerate
* regulate
* reconcile
* motivate
* interpret
* compare
* elucidate

### Preferred qualifiers

The style uses qualifiers constantly, but with discipline:

* roughly
* approximately
* modestly
* significantly
* strongly
* plausibly
* likely
* uncertain
* representative
* fiducial
* dynamically important
* of order
* factor of few
* to good accuracy
* in this limit
* under these assumptions

These qualifiers narrow the claim instead of weakening it.

### Preferred connective phrases

A phrase bank close to this style:

* “A key result is…”
* “A useful way to think about this is…”
* “This is because…”
* “This immediately shows that…”
* “It is useful to assess whether…”
* “For reference…”
* “As an example…”
* “We stress that…”
* “This does not prove…, but it does suggest…”
* “Taken together…”
* “One uncertainty in applying these results is…”
* “A natural next step is…”

### Words to avoid

Avoid too much startup-adjacent or hype language (use it but use it sparingingly):

* groundbreaking
* transformative
* revolutionary
* elegant solution
* leverages
* unlocks
* game-changing
* disruptive
* unprecedented

The author’s prose earns significance through argument, not adjectives.

---

## 8. How to handle equations in this style

### Every equation needs a job

Do not drop equations as decoration. Introduce them with a sentence that says what they are for.

Good lead-ins:

* “The key simplification is…”
* “The relevant ratio is…”
* “This yields a compact expression for…”
* “A good approximate solution is…”

### Name the controlling parameters

The author loves identifying “the two key dimensionless parameters” or “the relevant criterion.”

That is a hallmark. Always say which combinations of variables actually matter.

### Give the limiting form

After a general expression, simplify it in the physically relevant limit.

This is a core move:
general result → asymptotic limit → physical interpretation.

### Anchor with a concrete example

After a derivation, give one numeric anchor:
Milky Way-like values, dwarf galaxy values, group-mass halo values, etc.

This grounds the prose.

### State the validity domain

Always say when the estimate applies:
so long as diffusion dominates,
provided the transport time is shorter than the loss time,
assuming the CRs are proton-dominated,
neglecting pionic losses,
interior to the sonic point,
at large radii.

This is crucial to sounding like this author.

### Recast in intuitive form

Whenever possible, rewrite the result as:

* a ratio to a familiar rate
* a dimensionless pressure fraction
* a comparison of timescales
* a threshold condition

That is how the prose bridges formal derivation and physical understanding.

---

## 9. How to write about figures and tables

In these papers, captions are not labels. They are mini-arguments.

A strong caption in this style does four things:

* states what is plotted
* names what parameter is being varied
* explains what pattern matters
* tells the reader what physical conclusion to draw

In the main text, figure discussion follows a pattern:

1. announce what the figure shows
2. point out the key pattern
3. explain the pattern physically
4. connect it back to the analytic expectation

Example skeleton:
“Figure X shows the radial profiles of…. The most striking feature is that…. This occurs because…. The behavior bifurcates near…. This is consistent with the analytic estimate in equation….”

Use tables similarly: summarize parameters and outcomes, then tell the reader why those columns matter.

---

## 10. How to handle prior literature

The literature review style is highly functional.

Do not cite paper after paper without organizing them. Group prior work by role:

* foundational analytic models
* numerical simulations
* phenomenological/observational constraints
* related but different approaches

Then say what each group established and what remains unresolved.

The style is respectful but not deferential. It will say:
“Earlier work found X.”
“Our results are consistent with Y.”
“A key difference from Z is…”
“This conclusion is qualitatively similar to…”
“Our results strongly disfavor…”

That balance is important: precise, fair, and willing to differentiate.

---

## 11. How to handle uncertainty, caveats, and speculation

This is one of the author’s best habits.

### State assumptions early

Do not smuggle them in later. Say them where they first matter.

### Name the dominant uncertainty

Not every uncertainty. The dominant one.

These papers keep returning to the same core uncertainty rather than opening infinite side doors.

### Put caveats near the relevant claim

Do not banish all caveats to the end. Put them where they sharpen interpretation.

### Label speculation explicitly

Use a clear ladder:

* “shows” for direct result
* “suggests” for reasonable inference
* “likely” for strong but indirect inference
* “may” for open possibility
* “we speculate” for bolder extrapolation

### Never let uncertainty erase the argument

This author does not say “everything is uncertain.” He says, effectively:
“Given these uncertainties, here is what still seems robust.”

That is the ideal.

---

## 12. The characteristic rhetorical moves

There are recurring moves worth teaching directly to an AI.

### The simplification move

“A key simplification can be made…”

Use when reducing a complex system to its controlling balance.

### The regime split

“For X above Y…, by contrast for X below Y…”

Use when the paper’s main result depends on a threshold.

### The physical interpretation move

“Physically, this is because…”

Use after a derivation or figure.

### The calibration move

“To assess the implications, we must calibrate…”

Use when an uncertain parameter must be tied to data.

### The robustness move

“This does not prove…, but it does show…”

Use when evidence is suggestive but not definitive.

### The synthesis move

“Taken together…”

Use to join analytic, numerical, and observational threads.

---

## 13. What not to imitate too literally

Do not parody the style by repeating the same stock phrases every paragraph. The goal is not to say “physically” and “we reiterate” every few lines. The goal is to preserve the logic:

* explicit structure
* interpretable parameters
* regime-based reasoning
* calibrated confidence
* theory tied to observable consequences

Also do not become so cautious that the prose loses force. These papers are careful, but they still make claims.

---

## 14. A practical template for AI agents

Here is a compact operating template.

### Opening paragraph template

“X plays a key role in Y, but the physics of Z remains uncertain. Existing models suggest A, B, and C, yet it is not clear whether [specific mechanism] can account for [specific phenomenon]. In this work, we isolate that mechanism and ask under what conditions it is dynamically important.”

### Theory paragraph template

“To estimate this, we consider a simple model in which…. The relevant ratio is…. In the limit…, the governing equation reduces to…. This implies that the problem is controlled primarily by…. Physically, this is because….”

### Validation paragraph template

“Figure X shows that the numerical solutions follow the same overall trends. In particular, for…, the system…, whereas for…, it…. The agreement is not exact, but it is good enough to support the analytic scaling.”

### Caveated implication template

“Taken together, these results suggest that…. This conclusion depends on…, which remains uncertain, but it appears robust so long as…. An important consequence is that….”

---

## 15. Revision checklist

Before calling a draft finished, ask:

* Does the introduction identify a real uncertainty rather than just announce a topic?
* Does each section say why it exists?
* Does each equation earn its place and get interpreted in words?
* Have I identified the key dimensionless parameters or threshold conditions?
* Have I made the main regime split explicit?
* Have I given at least one fiducial numerical example?
* Have I translated at least one theoretical result into an observable or intuitive ratio?
* Are my caveats local, specific, and ranked by importance?
* Am I making the strongest claim justified by the assumptions?
* Does the discussion compare with prior work and point to a next test?

---

## 16. Condensed system prompt for an AI agent

Write in a calm, theory-driven astrophysical style. Begin from a meaningful physical uncertainty, not from hype. Narrow quickly to the specific mechanism of interest. Use simple models openly and purposefully. Identify the controlling dimensionless parameters and organize the argument around regime splits. Introduce each equation with a reason, define terms clearly, simplify in relevant limits, and immediately interpret the result physically. Re-express formal results in intuitive or observable forms. Use active first-person plural, precise transitions, and repeated technical nouns rather than decorative synonyms. Calibrate confidence carefully: show what follows directly, suggest what is inferred, and label speculation. Integrate figures and tables as parts of the argument. End by stating what is robust, what remains uncertain, how the results compare to prior work, and what observation or simulation would test the idea next.

The shortest summary is: **make the argument feel inevitable, but never sound more certain than the physics allows.**
