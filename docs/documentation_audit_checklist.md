# AthenaK Documentation Audit Checklist

Use the checkboxes below to track progress. Whenever you update a document, record the findings in `docs/documentation_audit_log.md` **before** checking the box.

## Urgent To Dos:
- [ ] make sure that all math renders properly

## Sphinx Site (docs/build/html)

- [ ] Landing Page — `docs/build/html/index.html` (`_sources/index.rst.txt`)
- [ ] Overview — `docs/build/html/overview.html` (`_sources/overview.md.txt`)
- [ ] Quickstart — `docs/build/html/quickstart.html` (`_sources/quickstart.md.txt`)
- [ ] Configuration Guide — `docs/build/html/configuration.html` (`_sources/configuration.md.txt`)
- [ ] Building Guide — `docs/build/html/building.html` (`_sources/building.md.txt`)
- [ ] Running Guide — `docs/build/html/running.html` (`_sources/running.md.txt`)
- [ ] Troubleshooting — `docs/build/html/troubleshooting.html` (`_sources/troubleshooting.md.txt`)
- [ ] Kokkos Guide — `docs/build/html/kokkos_guide.html` (`_sources/kokkos_guide.md.txt`)
- [ ] Contributing Docs — `docs/build/html/contributing_docs.html` (`_sources/contributing_docs.md.txt`)
- [ ] Glossary — `docs/build/html/glossary.html` (`_sources/glossary.rst.txt`)
- [ ] CGM Cooling Flow (Metals) — `docs/build/html/cgm_cooling_flow_metals.html` (`_sources/cgm_cooling_flow_metals.md.txt`)

### Examples

- [ ] Example: Binary Merger — `docs/build/html/examples/binary_merger.html`
- [ ] Example: Blast Wave — `docs/build/html/examples/blast_wave.html`
- [ ] Example: MRI Turbulence — `docs/build/html/examples/mri_turbulence.html`
- [ ] Example: Shock Tube — `docs/build/html/examples/shock_tube.html`
- [ ] Example: Turbulence — `docs/build/html/examples/turbulence.html`

### Flowcharts

- [x] Runtime Flowchart — `docs/build/html/flowcharts/runtime.html`
- [x] System Architecture Flowchart — `docs/build/html/flowcharts/system_architecture.html`

### Migration Guides

- [ ] Migration Index — `docs/build/html/migration/index.html`
- [ ] Migration – Common Gotchas — `docs/build/html/migration/common_gotchas.html`

### Module Reference

- [ ] Modules Overview — `docs/build/html/modules/index.html`
- [ ] Mesh Module — `docs/build/html/modules/mesh.html`
- [ ] Coordinates Module — `docs/build/html/modules/coordinates.html`
- [ ] Task List Module — `docs/build/html/modules/tasklist.html`
- [ ] Driver Module — `docs/build/html/modules/driver.html`
- [x] Hydro Module — `docs/build/html/modules/hydro.html`
- [x] MHD Module — `docs/build/html/modules/mhd.html`
- [ ] Dyn GRMHD Module — `docs/build/html/modules/dyn_grmhd.html`
- [ ] Ion-Neutral Module — `docs/build/html/modules/ion_neutral.html`
- [ ] Shearing Box Module — `docs/build/html/modules/shearing_box.html`
- [ ] Diffusion Module — `docs/build/html/modules/diffusion.html`
- [x] Source Terms Module — `docs/build/html/modules/srcterms.html`
- [ ] Outputs Module — `docs/build/html/modules/outputs.html`
- [ ] Particles Module — `docs/build/html/modules/particles.html`
- [ ] Radiation Module — `docs/build/html/modules/radiation.html`
- [ ] Reconstruction Module — `docs/build/html/modules/reconstruction.html`
- [ ] Riemann Solvers Module — `docs/build/html/modules/riemann_solvers.html`
- [ ] Problem Generators Module — `docs/build/html/modules/pgen.html`
- [ ] Boundaries Module — `docs/build/html/modules/boundaries.html`
- [ ] Z4c Module — `docs/build/html/modules/z4c.html`
- [ ] EOS Module — `docs/build/html/modules/eos.html`

### Reference Material

- [ ] Reference Index — `docs/build/html/reference/api_reference.html`
- [ ] File Reference — `docs/build/html/reference/file_reference.html`
- [ ] Input Parameters Reference — `docs/build/html/reference/input_parameters.html`

### Miscellaneous Pages

- [ ] Search Page — `docs/build/html/search.html`
- [ ] General Index — `docs/build/html/genindex.html`

> _Reminder_: The repository stores pre-built HTML alongside `_sources/*.txt` extracts. Update the Markdown/ReST sources under `docs/source/` wherever possible so future rebuilds stay consistent.
