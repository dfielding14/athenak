# Clebsch Stream Follow-Up Plan

## Direct quote from Camillo

> Hi Drummond,
>
>            what the LLM is suggesting is not good. In 6.2 it assumes that "the dominant interactions are local in scale, so that |p|, |q|, and |ξ| are all of order K". I don't see a reason for that, if the two potentials are picked at random everything interacts with equal probability, so there will be terms in the summation where |q| is of order K=|\xi| and |p| of order 1. In fact these are the ones which will give the leading order. So, my choice for the exponents m_1 and m_2 (assuming they are equal) would be -\zeta-4 and not -\zeta/2-4 as suggested by the LLM. But then it's important that the summation in (30) takes place over ALL frequencies and you don't instead take only p and q of approximately the same order...
>
> That is, if we are not missing something and, when taking the product of two Clebsch potentials, the probability does something funky and the product ends up being more regular than the separate potentials. Maybe Svitlana and I should think about this...
>
> Best,
>
>                                 Camillo

> Hi Drummond,
>
>              I did a bit more sanity checking. And I don't think Svitlana and I are missing anything, actually. If you follow what the LLM is telling you, essentially it is that to get the correct spectrum for the vector product of the two Clebsch potential, you need a spectrum for each of the factors which is 1/2 as regular as the spectrum you are trying to get for the product. Note that the derivatives don't really play any role in its heuristic, it's just this assumption that only interactions at roughly the same frequency matter in your summation in (30). So, let's look at this rule in 1d and without derivatives, in a much better understood situation.
>
> Let's say we were trying to create a random function in 1d of regularity 1/2 by putting i.i.d Gaussian random variables on the Fourier coefficients. Then (d-1)=0 and so what you put on the Fourier coefficient of frequency k is |k|^{-1}. This gives you precisely a 1d Brownian motion, right? So, now take the product of two such random function. The LLM in your writeup is instead suggesting that you should put \alpha/2 regularity in each factor to get \alpha regularity in the product. With that heuristic the spectrum of the product of two Brownian motions should be compatible with (almost) 1, right? But the product of two Brownian motions is instead as a regular as a single Brownian motion. Trying in fact asking this specific question to the LLM: "is the product of two Brownian motions more regular than a single Brownian motion, i.e. more regular than Hoelder 1/2?". The LLM will reply that it is not, in fact it will tell you that the regularity of the product of two Brownian motions is still 1/2-\varepsilon almost surely, and almost surely not 1/2+\varepsilon ... at least that's what answers to me :-)
>
> Best,
>
>                                        Camillo

## Why this matters

Camillo's criticism targets the weakest part of the current 3D Clebsch writeup: the continuum asymptotic argument that assumes the dominant contributions to

\[
\hat{\mathbf u}(\xi)
=
-\sum_{p+q=\xi}
(p\times q)\,\hat\phi_1(p)\hat\phi_2(q)
\]

come from the scale-local sector `|p| ~ |q| ~ |\xi|`.

That assumption is not built into the exact AthenaK shell-response tensor. The code sums over all retained `p,q` pairs. But the assumption does appear in the documentation, and it also appears indirectly in the asymptotic seed used to initialize the 3D shell fit.

So there are two separate questions:

1. Is the mathematical derivation in the note wrong or at least unjustified?
2. Even if the exact fitter is more robust than the note, is the current asymptotic seed in the code biasing the 3D stream construction in the wrong direction?

## Current assessment

The present working assessment should be:

- The exact discrete shell-response tensor in the code is still the right object to fit.
- The current asymptotic explanation in the LaTeX note is not trustworthy.
- The current seed `beta_asym = (expo + 8) / 4` in the 3D stream fitter should now be treated as provisional rather than justified.
- The 3D stream documentation should be revised immediately so it does not present the `m_1 + m_2 = \zeta + 8` balance as an established continuum law.

## Files implicated

### Documentation

- `docs/fourier_coefficient_variance_note.tex`
- `docs/scalar_mixing_master.tex`
- `docs/scalar_mixing_stream_function_notes.md`
- possibly `docs/scalar_only_vs_kinematic.tex` if the same heuristic is echoed there

### Code

- `src/pgen/scalar_mixing.cpp`

Key code location:

- the exact shell-response tensor and 3D shell fit live in the `BuildClebschShellResponseTensor` and `FitClebschShellWeights` path
- the current asymptotic seed appears in
  - `beta_asym = 0.25 * (cfg.expo + 8.0)`

## Immediate response to Camillo

We should respond quickly and directly:

- acknowledge that his criticism is substantive and likely correct
- distinguish between the exact tensor fit in the code and the weaker continuum heuristic in the note
- say explicitly that the derivation and the asymptotic seed are being revisited
- commit to sending back a corrected note plus numerical evidence

A reasonable summary sentence is:

> The exact retained-mode shell tensor is still the object we fit in AthenaK, but the continuum asymptotic argument in the note is too weak as written, and we are revisiting both the derivation and the asymptotic seed used to initialize the fit.

## Next-step plan

### Phase 1: Repair the mathematical story

Goal: replace the unjustified locality heuristic with a mathematically defensible analysis.

Tasks:

1. Re-derive the 3D Clebsch asymptotics starting from the exact second-moment identity
   - separate the convolution into interaction sectors:
     - local-local
     - low-high
     - high-low
     - high-high with near-cancellation
2. Compare the result with paradifferential/Bony-product intuition
   - check whether the Clebsch product inherits the weaker factor regularity rather than gaining regularity from both factors
3. Write down the symmetric candidate law suggested by Camillo
   - if `m_1 = m_2 = m`, the proposed guide is
     - `m = \zeta + 4`
   - equivalently, the coefficient amplitude guide would be
     - `|\hat\phi_j(k)| ~ |k|^{-(\zeta+4)/2}`
4. Identify precisely what remains heuristic after this re-derivation
   - continuum asymptotics
   - finite-band effects
   - lattice shell-count effects
   - importance-sampling effects

Deliverable:

- a corrected mathematical note that states what is exact, what is heuristic, and what is still open

### Phase 2: Audit the current code against the new picture

Goal: determine whether the implementation is only poorly explained, or whether it is genuinely biased by the bad asymptotic seed.

Tasks:

1. Instrument the shell-response tensor in the 3D stream path
   - record how much of each output shell comes from:
     - low-high
     - high-low
     - local-local
     - high-high cancellation sectors
2. Measure those sector contributions for several candidate potential laws
   - current seed
   - Camillo-style symmetric seed
   - a few nearby exponents for sensitivity
3. Check whether the current fitted result depends strongly on the initial seed
4. Check whether the current potential cutoff `phi_nhigh ~ nhigh` is too restrictive
   - if nonlocal interactions dominate, the potential support may need to extend beyond the target velocity band

Deliverable:

- a short numerical audit of which interaction sectors actually dominate the discrete retained tensor used by AthenaK

### Phase 3: Decide whether the 3D fitter must change

Goal: modify the code only if the audit shows that the current seed or cutoff materially biases the result.

Possible code changes:

1. Replace the asymptotic seed
   - current:
     - `beta_asym = (expo + 8) / 4`
   - candidate replacement:
     - `beta_asym = (expo + 4) / 2`
2. Reduce or remove the imposed asymmetric split in the initial guess
   - current seed uses a fixed offset `beta_delta`
   - if the symmetric law is the correct guide, the default start should probably be symmetric
3. Extend the potential shell support
   - allow `phi_nhigh` to exceed the target velocity `nhigh`
4. Add diagnostics
   - save the seed model used
   - save sector-resolved shell-response summaries
   - save potential cutoff information

Acceptance criterion:

- the realized 3D stream velocity spectrum should become less sensitive to arbitrary seed choices
- any change should preserve or improve the current spectrum quality and eliminate hidden dependence on the old heuristic

### Phase 4: Correct the documentation

Goal: ensure the writeup no longer overclaims.

Required changes:

1. `docs/fourier_coefficient_variance_note.tex`
   - remove the scale-local asymptotic section in its current form
   - replace it with:
     - the exact variance identity
     - a sector-based discussion
     - a clear statement that the continuum asymptotics are under active revision
2. `docs/scalar_mixing_master.tex`
   - remove the statement that `m_1 + m_2 = \zeta + 8` is the governing asymptotic rule
   - replace it with a statement that the code fits the exact discrete shell-response tensor, and that the continuum guide is not yet settled
3. `docs/scalar_mixing_stream_function_notes.md`
   - update the description of the 3D shell-fit seed
4. Any other note that currently repeats the same asymptotic claim

Deliverable:

- rebuilt PDFs and Markdown notes with no unsupported continuum claim

### Phase 5: Add focused regression tests

Goal: make sure this issue does not come back in disguised form.

Tests to add:

1. A toy random-product diagnostic
   - for example a 1D Fourier-series product sanity test mirroring Camillo's Brownian example
2. A 3D stream seed-sensitivity test
   - compare realized spectra from multiple initial seeds
3. A 3D stream cutoff-sensitivity test
   - vary potential support beyond `nhigh`
4. A sector-dominance diagnostic test
   - ensure the code can report which interaction sectors dominate the fitted response

Deliverable:

- reproducible evidence that the revised 3D stream setup is not resting on an untested asymptotic shortcut

## Recommended execution order

1. Send Camillo a brief acknowledgement and say we are revisiting both the derivation and the fit seed.
2. Patch the LaTeX note immediately so it no longer presents the locality heuristic as established fact.
3. Instrument the 3D shell-response tensor and run the sector-dominance audit.
4. Test whether the current 3D fitter is sensitive to the old asymptotic seed.
5. If necessary, revise the seed and potential cutoff in `scalar_mixing.cpp`.
6. Rebuild the notes and send Camillo a corrected writeup with numerical evidence.

## Minimal immediate actions

If time is tight, the first three actions should be:

1. stop repeating the `m_1 + m_2 = \zeta + 8` claim in the notes as if it were settled
2. check whether the code's exact shell fit is robust to replacing the current seed
3. quantify whether low-high interactions actually dominate the discrete shell-response tensor

That will tell us quickly whether this is only a documentation problem or a real implementation problem.
