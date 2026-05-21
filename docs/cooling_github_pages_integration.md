# Cooling Documentation GitHub Pages Integration

This note records the documentation review for integrating the cooling feature
branch into the `origin/gh-pages` documentation framework.

## Review Findings

1. The cooling feature documentation is compatible with the current GitHub Pages
   build stack.  The live Pages branch uses Sphinx, MyST, and
   `docs/Makefile`; the page `docs/source/modules/cooling.md` builds as a MyST
   Markdown document without syntax changes.
2. The cooling figures can remain under `docs/source/modules/figures/`.  Sphinx
   copies the relative Markdown image targets into the generated HTML correctly.
3. The page must be registered in the Pages branch toctree.  If
   `modules/cooling` is not added to `docs/source/index.md`, Sphinx can still
   parse the file when explicitly linked, but it will not appear in the main
   navigation.
4. The existing `origin/gh-pages:docs/source/modules/srcterms.md` still
   documents the removed legacy `<hydro>/ism_cooling` and `<hydro>/cgm_cooling`
   modes.  That page should be updated when cooling docs are integrated so the
   live documentation does not simultaneously advertise removed legacy inputs
   and the new standalone `<cooling>` block.
5. `docs/source/modules/cooling_implementation_plan.md` is useful as an
   engineering design record, but it should not be added to the public module
   navigation unless we intentionally want design-history docs on the live site.
   The user-facing module page should be `docs/source/modules/cooling.md`.

## Files To Carry Into `gh-pages`

Copy or merge these files from `feature/cooling`:

```text
docs/source/modules/cooling.md
docs/source/modules/figures/cooling_constant_solution.png
docs/source/modules/figures/cooling_model_rate_errors.png
docs/source/modules/figures/cooling_table_interpolation_convergence.png
```

The examples and validator are not Sphinx sources, but the cooling page points
to them and they should travel with the feature branch:

```text
inputs/cooling/
tools/validate_cooling_table.py
tst/test_suite/cooling/
```

## Required Pages Branch Edits

In `docs/source/index.md`, add `Cooling` to the Module Reference bullet list:

```markdown
- Numerical methods: [Reconstruction](modules/reconstruction.md), [Riemann Solvers](modules/riemann_solvers.md), [EOS](modules/eos.md), [Diffusion](modules/diffusion.md), [Outputs](modules/outputs.md), [Boundaries](modules/boundaries.md), [Source Terms](modules/srcterms.md), [Cooling](modules/cooling.md), [Shearing Box](modules/shearing_box.md), [Problem Generators](modules/pgen.md)
```

Also add `modules/cooling` to the hidden Modules toctree:

```text
modules/srcterms
modules/cooling
modules/shearing_box
```

In `docs/source/modules/index.md`, update the support-system count and add a
row:

```markdown
## Support Systems (6+ modules)

| **Cooling Source Terms** | `src/srcterms/` | Configurable radiative cooling and heating | [cooling.md](cooling.md) |
```

In `docs/source/modules/srcterms.md`, add a short pointer near the start of the
source-term family section:

```markdown
## Cooling

The generalized cooling and heating implementation is documented separately in
[Cooling Source Terms](cooling.md). New inputs should use the standalone
`<cooling>` block; legacy `<hydro>/ism_cooling`, `<mhd>/ism_cooling`,
`<hydro_srcterms>/ism_cooling`, `<mhd_srcterms>/ism_cooling`, and corresponding
`cgm_cooling` keys are intentionally not supported by the cooling feature.
```

Then remove or rewrite the old `### ISM Cooling (<hydro>)` and
`### CGM Cooling (<hydro>)` subsections so they do not advertise legacy
configuration keys as current behavior.

## Validation Performed

I tested the integration in a detached worktree created from `origin/gh-pages`
after fetching the remote branch.  The test copied the cooling page and figures,
registered `modules/cooling` in the main toctree, added a module-index row, and
inserted the `srcterms.md` pointer above.

The Pages-branch build command was:

```bash
cd docs
make clean html SPHINXOPTS="-W --keep-going"
```

Result: the build succeeded with warnings treated as errors.  The Sphinx output
read `modules/cooling`, generated the HTML page, and copied all three cooling
figures.

## Integration Checklist

1. Merge or cherry-pick the cooling feature docs and assets into `gh-pages`.
2. Apply the toctree and module-index edits above.
3. Update `srcterms.md` so it no longer contradicts the new cooling block.
4. Run `cd docs && make clean html SPHINXOPTS="-W --keep-going"`.
5. Push `gh-pages` only after the local Sphinx build passes.
