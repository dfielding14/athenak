# GitHub Pages Integration: Turbulence Driving

This feature branch retains only the documentation source and generated
projection asset needed for integration into the publication branch. Overlay
the following files onto `gh-pages`:

```text
docs/source/modules/turbulence_driving.md
docs/source/examples/turbulence_driving.md
docs/source/_static/turbulence_driving_projections.png
```

Wire the module page into `docs/source/modules/index.md` alongside the other
source-term modules and add the example page to `docs/source/index.md` under
worked examples. Include the page paths in any hidden toctrees:

```text
modules/turbulence_driving
examples/turbulence_driving
```

## Documentation Contract

The module page records the public API, equations, AMR invariants, restart
state, and explicit rejection of prototype option names. The example page
contains only simulation-derived checks and a projection figure regenerated
from the committed inputs with:

```bash
python vis/python/plot_turbulence_driving_projection.py \
  --panel "Constant Edot=run/turb_edot_fixed/bin/turb_edot_fixed.force.00004.bin" \
  --panel "Tiled + Gaussian include=run/turb_tiled_include/bin/turb_tiled_include.force.00004.bin" \
  --panel "AMR + Gaussian exclude=run/turb_amr_exclude/bin/turb_amr_exclude.force.00004.bin" \
  --output docs/source/_static/turbulence_driving_projections.png
```

## Pages Validation

Validate the overlay against the live Pages tree rather than treating this
feature checkout as a standalone Sphinx project:

```bash
python -m pip install -r docs/requirements.txt
make -C docs clean html SPHINXOPTS="-W --keep-going"
```

The feature is ready to publish only when that strict Sphinx build succeeds
with the module/example navigation entries present.
