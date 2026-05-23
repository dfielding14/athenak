# Frame Tracking GitHub Pages Integration

This feature branch is based on `origin/main`, while the editable Sphinx/MyST
documentation tree lives on `origin/gh-pages`. The files below are intended to
be copied or merged into that docs branch without replacing existing broad
module pages.

## Drop-In Files

- `docs/source/modules/frame_tracking.md`
- `docs/source/modules/frame_tracking_next_steps.md`
- `docs/source/modules/frame_tracking_recipes.md`
- `docs/source/examples/cloud_crushing_snr.md`
- `docs/source/examples/trml_frame_tracking.md`
- `docs/FRAME_TRACKING_NEXT_STEPS.md`
- `docs/source/_static/sedov_boundary_cgs.png`
- `docs/source/_static/sedov_boundary_cgs.csv`
- `docs/source/_static/cloud_crushing_lowres_validation.png`
- `docs/source/_static/cloud_crushing_lowres_midplane.png`
- `docs/source/_static/cloud_crushing_lowres_dense_mass.png`
- `docs/source/_static/cloud_crushing_lowres_dense_mass.csv`
- `docs/source/_static/cloud_crushing_lowres_frame_tracker.csv`
- `docs/source/_static/cloud_crushing_lowres_summary.csv`
- `docs/source/_static/cloud_crushing_lowres_density_slices_6_vertical_equal_aspect.png`

Do not copy `docs/source/modules/pgen.md` or
`docs/source/modules/srcterms.md` from older versions of this feature branch
into `gh-pages`; those short feature notes were replaced by the dedicated
`modules/frame_tracking.md` page so the existing `gh-pages` module pages remain
intact.

## Suggested `docs/source/index.md` Additions

Add this link to the Module Reference bullet list:

```markdown
[Frame Tracking](modules/frame_tracking.md)
[Frame Tracking Production Plan](modules/frame_tracking_next_steps.md)
[Frame Tracking Recipes And Migration](modules/frame_tracking_recipes.md)
```

Add this entry to the hidden Modules toctree:

```markdown
modules/frame_tracking
modules/frame_tracking_next_steps
modules/frame_tracking_recipes
```

Add these links to the Worked Examples bullet list:

```markdown
- [Sedov Cloud Crushing With Frame Tracking](examples/cloud_crushing_snr.md)
- [TRML Frame Tracking Example](examples/trml_frame_tracking.md)
```

Add these entries to the hidden Examples toctree:

```markdown
examples/cloud_crushing_snr
examples/trml_frame_tracking
```

## Suggested `docs/source/modules/index.md` Addition

Add a support-system row:

```markdown
| **Frame Tracking** | `src/srcterms/frame_tracker.*` | Shared moving-frame controller for hydro/MHD examples | [frame_tracking.md](frame_tracking.md) |
```

## Validation Commands

From a worktree on `origin/gh-pages` after copying the files and applying the
index edits:

```bash
python3 -m venv docs/.venv
docs/.venv/bin/pip install -r docs/requirements.txt
docs/.venv/bin/sphinx-build -b html docs/source docs/_build/html
```

The feature docs use MyST Markdown, local PNG/CSV assets under
`docs/source/_static`, and relative links that match the existing GitHub Pages
framework. `docs/FRAME_TRACKING_NEXT_STEPS.md` is the maintained roadmap
published alongside this integration handoff; it records remaining benchmark
and scientific-validation gates without overstating low-resolution runs.
