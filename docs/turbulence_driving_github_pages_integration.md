# GitHub Pages Integration: Turbulence Driving

This feature branch contains Pages-ready source files without copying the
entire publication branch into the implementation branch.  When integrating
documentation into `gh-pages`, bring across:

```text
docs/source/modules/turbulence_driving.md
docs/source/examples/turbulence_driving.md
docs/source/_static/turbulence_driving_projections.png
```

Add these entries to the corresponding toctrees on `gh-pages`:

```markdown
# docs/source/modules/index.md
- [Turbulence Driving](turbulence_driving.md)

# docs/source/index.md, under Worked Examples
- [Localized and AMR Turbulence Driving](examples/turbulence_driving.md)
```

and include the page paths in the hidden `toctree` blocks:

```text
modules/turbulence_driving
examples/turbulence_driving
```

## Validation

Validate the rendered documentation against the live Pages tree after
overlaying these files and nav entries:

```bash
python -m pip install -r docs/requirements.txt
make -C docs clean html SPHINXOPTS="-W --keep-going"
```

The module page documents the normalization and AMR invariants.  The example
page includes only projections generated from the three committed input
problems using `vis/python/plot_turbulence_driving_projection.py`.
