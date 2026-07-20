# AthenaK Documentation Audit Checklist

Use the checkboxes below to track progress. Whenever you update a document, record the findings in `docs/documentation_audit_log.md` **before** checking the box.

## Current Rebuild Status (2026-05-30)

- [x] Rebuilt and reviewed the current Sphinx HTML site with the declared
  requirements in isolated `/tmp/athenak-pic-docs-venv`. The available system
  Sphinx 2.3.1 stack is missing the required `myst_parser` extension, so use
  that isolated dependency environment for the current source candidate.
- [x] Removed the stale glossary navigation entry. No glossary page exists in
  the current 58-page site.

## Sphinx Site (docs/build/html)

- [x] Landing Page — `docs/build/html/index.html` (`_sources/index.rst.txt`)
- [x] Overview — `docs/build/html/overview.html` (`_sources/overview.md.txt`)
- [x] Quickstart — `docs/build/html/quickstart.html` (`_sources/quickstart.md.txt`)
- [x] Configuration Guide — `docs/build/html/configuration.html` (`_sources/configuration.md.txt`)
- [x] Building Guide — `docs/build/html/building.html` (`_sources/building.md.txt`)
- [x] Running Guide — `docs/build/html/running.html` (`_sources/running.md.txt`)
- [x] Troubleshooting — `docs/build/html/troubleshooting.html` (`_sources/troubleshooting.md.txt`)
- [x] Kokkos Guide — `docs/build/html/kokkos_guide.html` (`_sources/kokkos_guide.md.txt`)
- [x] Contributing Docs — `docs/build/html/contributing_docs.html` (`_sources/contributing_docs.md.txt`)
- [x] CGM Cooling Flow (Metals) — `docs/build/html/cgm_cooling_flow_metals.html` (`_sources/cgm_cooling_flow_metals.md.txt`)

### Examples

- [x] Example: Binary Merger — `docs/build/html/examples/binary_merger.html`
- [x] Example: Blast Wave — `docs/build/html/examples/blast_wave.html`
- [x] Example: MRI Turbulence — `docs/build/html/examples/mri_turbulence.html`
- [x] Example: Shock Tube — `docs/build/html/examples/shock_tube.html`
- [x] Example: Turbulence — `docs/build/html/examples/turbulence.html`

### Flowcharts

- [x] Runtime Flowchart — `docs/build/html/flowcharts/runtime.html`
- [x] System Architecture Flowchart — `docs/build/html/flowcharts/system_architecture.html`

### Migration Guides

- [x] Migration Index — `docs/build/html/migration/index.html`
- [x] Migration – Common Gotchas — `docs/build/html/migration/common_gotchas.html`

### Module Reference

- [x] Modules Overview — `docs/build/html/modules/index.html`
- [x] Mesh Module — `docs/build/html/modules/mesh.html`
- [x] Coordinates Module — `docs/build/html/modules/coordinates.html`
- [x] Task List Module — `docs/build/html/modules/tasklist.html`
- [x] Driver Module — `docs/build/html/modules/driver.html`
- [x] Hydro Module — `docs/build/html/modules/hydro.html`
- [x] MHD Module — `docs/build/html/modules/mhd.html`
- [x] Dyn GRMHD Module — `docs/build/html/modules/dyn_grmhd.html`
- [x] Ion-Neutral Module — `docs/build/html/modules/ion_neutral.html`
- [ ] Shearing Box Module — `docs/build/html/modules/shearing_box.html`
- [ ] Diffusion Module — `docs/build/html/modules/diffusion.html`
- [x] Source Terms Module — `docs/build/html/modules/srcterms.html`
- [x] Outputs Module — `docs/build/html/modules/outputs.html`
- [x] Particles Module — `docs/build/html/modules/particles.html`
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

### Engineering Notes

- [x] MHD-PIC CR-Hall Paper-to-Code Map — `docs/source/engineering/pic_cr_hall_code_map.md`
- [x] MHD-PIC Project Ethos — `docs/source/engineering/pic_project_ethos.md`
- [x] Mignone-R2 Full-Hall Shock Setup — `docs/source/engineering/pic_mignone_r2_shock_setup.md`
- [x] MHD-PIC Clean-Launch Runbook — `docs/build/html/engineering/pic_clean_launch_runbook.html`
- [x] MHD-PIC Supported Toolchains — `docs/build/html/engineering/pic_supported_toolchains.html`
- [x] MHD-PIC Runtime Model Contract — `docs/build/html/engineering/pic_mhd_model_contract.html`
- [x] MHD-PIC AMR Lifetime and Interface Policy — `docs/build/html/engineering/pic_amr_lifetime_and_interface_policy.html`
- [x] MHD-PIC Q-016 Bounded Particle Provenance and Spectra — `docs/build/html/engineering/pic_q016_particle_provenance_spectra.html`

> _Reminder_: The repository stores pre-built HTML alongside `_sources/*.txt` extracts. Update the Markdown/ReST sources under `docs/source/` wherever possible so future rebuilds stay consistent.
