# AthenaK Documentation Audit Ledger

Audit started: 2026-05-29

## Baseline

- Documentation source branch: `origin/gh-pages` at `8eb329959244d68bffc3b2347432a9c85a622394`.
- Public implementation baseline: `origin/main` at `886dd2a1437e45a3a30b3eeebf2adfa838328f73`.
  `origin/gh-pages` is a docs/deployment-only descendant of this commit and has
  no `src/`, `inputs/`, `tst/`, or `vis/` changes relative to `origin/main`.
- Development comparison: `origin/development` at
  `01f24197e850cc6396fde7408e3174de5694ef0a` includes additional code, inputs,
  tests, and documentation that are not stable public behavior for this site.
- Deployment: `.github/workflows/docs.yml` builds `docs/source/` on pushes to
  `gh-pages` and deploys `docs/build/html` through GitHub Pages.
- Baseline build: `cd docs && make clean html SPHINXBUILD=/tmp/athenak-doc-audit-venv/bin/sphinx-build SPHINXOPTS='-W --keep-going'`
  passed before edits; Sphinx read 55 source pages.
- Baseline rendered check: the deployed home page and the strict local render
  load the same navigation structure; the home-page
  `_static/athenak_fluid_sim.html` iframe is present and contains its live
  `canvas` element. It must remain present through final verification.

## Inventory Rules

Status values used below:

- `Pending`: no current-pass claim audit has been completed.
- `In progress`: the page has been read and evidence review or editing is active.
- `Checked`: material claims have been checked against public-baseline evidence.
- `Clear`: no edit is required after checking.
- `Edited`: evidence-backed edits have been made but final gates remain.
- `Reviewed`: an independent reviewer reported no unresolved actionable issue.
- `Verified`: strict build, link, rendered, and post-review checks are complete.

Sphinx's baseline environment includes all 55 Markdown documents below through
the `index` toctree and nested section toctrees. There are no unlinked
Markdown or reStructuredText source pages under `docs/source/` in the baseline.
Generated `search.html` and `genindex.html` are rendered-site checks, not
source-page claim audits.

## Published Source Pages

| Page path | Audience | Purpose | Evidence to inspect | Accuracy | Completeness | Usability | Required edits | Independent reviewer | Final verification |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `source/index.md` | mixed | Home page, routing, simulation banner | toctrees, `_static/athenak_fluid_sim.html`, render | Checked | Edited | Edited | Preserved simulation; improved routes and iframe title | Confucius: clean | Verified |
| `source/quickstart.md` | beginner | First clone/build/run/output path | build files, Sod/Orszag decks, CLI/output code | Checked | Edited | Edited | Corrected first-run, VTK and restart workflows | Poincare: clean | Verified |
| `source/building.md` | user/developer | CPU, MPI, GPU and custom builds | top-level CMake and Kokkos config | Checked | Edited | Edited | Corrected version, backend and custom-pgen build guidance | Poincare: clean | Verified |
| `source/configuration.md` | user | Input blocks and common parameters | parser and module constructors | Checked | Edited | Edited | Corrected implemented parameter and solver boundaries | Poincare: clean | Verified |
| `source/running.md` | user | Runtime CLI and execution workflows | `src/main.cpp`, driver and outputs | Checked | Edited | Edited | Corrected CLI, directories, outputs and restart use | Poincare: clean | Verified |
| `source/troubleshooting.md` | beginner/user | Diagnose failures | build, parser, restart and timestep paths | Checked | Edited | Edited | Replaced unsupported recovery advice | Poincare: clean | Verified |
| `source/overview.md` | mixed | Architecture and public capability orientation | `main.cpp`, `meshblock_pack.cpp`, constructors | Checked | Edited | Edited | Restricted overview to public-baseline behavior | Confucius: clean | Verified |
| `source/flowcharts/system_architecture.md` | developer | Component relationships | mesh/physics/task construction | Checked | Edited | Edited | Corrected module edges and coupling labels | Anscombe: clean | Verified |
| `source/flowcharts/runtime.md` | developer/expert | Startup and timestep sequence | main/driver/task lifecycle | Checked | Edited | Edited | Corrected parsing, mesh, driver and timestep ordering | Confucius: clean | Verified |
| `source/kokkos_guide.md` | developer/expert | Device-safe coding guidance | `src/athena.hpp`, representative kernels | Checked | Edited | Edited | Aligned aliases and wrapper guidance with public API | Confucius: clean | Verified |
| `source/modules/index.md` | user/developer | Module catalogue | public `src/` tree and child pages | Checked | Edited | Edited | Rebuilt supported module routing | Copernicus: clean | Verified |
| `source/modules/mesh.md` | developer/expert | Mesh packs and refinement | `src/mesh/`, parser and callbacks | Checked | Edited | Edited | Corrected static/multilevel/adaptive behavior | Heisenberg: clean | Verified |
| `source/modules/driver.md` | developer/expert | Lifecycle and integrators | `src/driver/`, `src/mesh/mesh.cpp` | Checked | Edited | Edited | Corrected CFL storage and stopping semantics | Carver: clean | Verified |
| `source/modules/tasklist.md` | developer/expert | Scheduling API | task registrations and coupled graphs | Checked | Edited | Edited | Corrected forcing placement and coupled graph ownership | Halley: clean after iterations | Verified |
| `source/modules/coordinates.md` | user/developer | Metrics and excision | `src/coordinates/`, GR inputs | Checked | Edited | Edited | Corrected DynGRMHD fragment and radiation-forced excision | Carver: clean after iterations | Verified |
| `source/modules/hydro.md` | user/developer | Hydro solver | `src/hydro/`, EOS and inputs | Checked | Edited | Edited | Corrected activation, solvers and task behavior | Heisenberg: clean | Verified |
| `source/modules/mhd.md` | user/developer | MHD and CT | `src/mhd/`, boundaries and inputs | Checked | Edited | Edited | Corrected CT, coupling and solver constraints | Heisenberg: clean | Verified |
| `source/modules/radiation.md` | user/developer | Radiation module | `src/radiation/`, shipped inputs | Checked | Edited | Edited | Corrected public activation and beam prerequisites | Anscombe: clean | Verified |
| `source/modules/z4c.md` | expert/developer | Numerical relativity | `src/z4c/`, pgen and Z4c decks | Checked | Edited | Edited | Flagged broken decks; documented AMR bridge and indexing | Halley: clean after iterations | Verified |
| `source/modules/dyn_grmhd.md` | expert/developer | Dynamical GRMHD | `src/dyn_grmhd/`, ADM/Z4c wiring | Checked | Edited | Edited | Distinguished ADM/Z4c and fixed/FOFC/EOS behavior | Carver: clean after iterations | Verified |
| `source/modules/ion_neutral.md` | user/developer | Ion-neutral coupling | `src/ion-neutral/`, construction and decks | Checked | Edited | Edited | Corrected IMEX order and metric incompatibility | Halley: clean after iterations | Verified |
| `source/modules/particles.md` | user/developer | Particle capabilities | `src/particles/`, particle input | Checked | Edited | Edited | Restricted public path to validated 3D periodic drift | Carver: clean | Verified |
| `source/modules/reconstruction.md` | expert/developer | Reconstruction methods | reconstruct code and fluid setup | Checked | Edited | Edited | Clarified DynGRMHD ownership and limits | Halley: clean | Verified |
| `source/modules/riemann_solvers.md` | expert/developer | Flux solver selection | solver constructors and DynGRMHD | Checked | Edited | Edited | Added DynGRMHD solver dispatch | Halley: clean | Verified |
| `source/modules/eos.md` | user/developer | EOS choices | EOS and DynGRMHD constructors | Checked | Edited | Edited | Added required DynGRMHD selector/state constraints | Carver: clean | Verified |
| `source/modules/diffusion.md` | user/developer | Explicit diffusion | diffusion code and task routing | Checked | Edited | Edited | Scoped saturation and DynGRMHD limitations | Carver: clean | Verified |
| `source/modules/outputs.md` | user/expert | Outputs and restarts | output sources and public readers | Checked | Edited | Edited | Corrected formats/restarts; qualified failing PDF path | Anscombe: clean | Verified |
| `source/modules/boundaries.md` | developer/expert | Boundary behavior | parser and boundary sources | Checked | Edited | Edited | Completed parser list and implementation limits | Carver: clean | Verified |
| `source/modules/srcterms.md` | user/developer | Public source terms | public `src/srcterms/`, development diff | Checked | Edited | Edited | Removed development-only support claims | Copernicus: clean | Verified |
| `source/modules/shearing_box.md` | user/developer | Shearing-box workflows | shearing tasks and MRI decks | Checked | Edited | Edited | Qualified 2D shifted/remap path and MRI2D state | Halley: clean | Verified |
| `source/modules/pgen.md` | developer | Problem generator guide | `src/pgen/`, CMake and decks | Checked | Edited | Edited | Corrected built-in/custom and callback contract | Anscombe: clean | Verified |
| `source/reference/input_parameters.md` | expert/user | Parameter catalogue | constructor/parser reads and inputs | Checked | Edited | Edited | Replaced unsupported authority claims with verified scope | Anscombe: clean | Verified |
| `source/reference/file_reference.md` | developer | File map | public source/input/visualization tree | Checked | Edited | Edited | Removed absent/development files | Anscombe: clean | Verified |
| `source/reference/api_reference.md` | developer | Interface reference | public headers and API implementations | Checked | Edited | Edited | Corrected callback and type references | Heisenberg: clean | Verified |
| `source/migration/index.md` | developer | Athena++ porting workflow | pgen/API/Kokkos patterns | Checked | Edited | Edited | Rebuilt around current public extension route | Confucius: clean | Verified |
| `source/migration/common_gotchas.md` | developer | Porting pitfalls | host/device and current API rules | Checked | Edited | Edited | Corrected current migration hazards | Confucius: clean | Verified |
| `source/examples/index.md` | beginner/user | Example selection | public shipped inputs | Checked | Edited | Edited | Removed unstable routes and clarified choices | Nash: clean | Verified |
| `source/examples/shock_tube.md` | beginner/user | Sod workflow | Sod deck, pgen and outputs | Checked | Edited | Edited | Executed and corrected output procedure | Nash: clean | Verified |
| `source/examples/blast_wave.md` | user | Blast workflow | blast pgen/input and outputs | Checked | Edited | Edited | Corrected custom-build/run path; executed smoke run | Nash: clean | Verified |
| `source/examples/turbulence.md` | user/expert | Turbulence workflow | public turbulence source/input | Checked | Edited | Edited | Restricted to public executable route; executed smoke run | Nash: clean | Verified |
| `source/examples/mri_turbulence.md` | user | MRI workflow | shearing pgens and inputs | Checked | Edited | Edited | Verified MRI3D; flagged non-runnable MRI2D | Nash: clean | Verified |
| `source/examples/binary_merger.md` | expert | BNS prerequisites | `sgrid_bns.cpp` and input | Checked | Edited | Edited | Documented external library/data gate | Nash: clean | Verified |
| `source/cgm_cooling_flow_metals.md` | user/expert | Development-only problem note | public/development source comparison | Checked | Edited | Edited | Marked absent from public runnable baseline | Copernicus, Nietzsche: clean | Verified |
| `source/tools/visualization.md` | user | Post-processing tools | public `vis/python/` utilities | Checked | Edited | Edited | Documented public reader/conversion paths only | Copernicus: clean | Verified |
| `source/contributing_docs.md` | maintainer | Documentation workflow | Sphinx and deployment workflow | Checked | Edited | Edited | Added evidence/review/build/banner requirements | Confucius: clean | Verified |
| `source/engineering/index.md` | maintainer/developer | Engineering record landing | public/development comparison | Checked | Edited | Edited | Added stable-versus-record boundary and links | Copernicus, Nietzsche: clean | Verified |
| `source/engineering/amr_turbulence_implementation.md` | maintainer | AMR turbulence record | public turbulence implementation | Checked | Edited | Edited | Added development-record warning | Nietzsche: clean | Verified |
| `source/engineering/divb_amr_completion_plan.md` | maintainer | AMR div(B) plan | public MHD/boundary tree | Checked | Edited | Edited | Added development-record warning | Nietzsche: clean | Verified |
| `source/engineering/divb_amr_final_status.md` | maintainer | AMR div(B) record | public source/tests | Checked | Edited | Edited | Added development-record warning | Nietzsche: clean | Verified |
| `source/engineering/divb_amr_implementation_plan.md` | maintainer | AMR div(B) design | public/development implementation | Checked | Edited | Edited | Added development-record warning | Nietzsche: clean | Verified |
| `source/engineering/divb_amr_implementation_summary.md` | maintainer | AMR div(B) summary | public source/tests | Checked | Edited | Edited | Added development-record warning | Nietzsche: clean | Verified |
| `source/engineering/face_field_correction_implementation.md` | maintainer | Face-field record | public MHD/boundary source | Checked | Edited | Edited | Added development-record warning | Nietzsche: clean | Verified |
| `source/engineering/particle_merge_cr_pushers.md` | maintainer | Particle/CR record | public particle implementation | Checked | Edited | Edited | Added development-record warning | Nietzsche: clean | Verified |
| `source/engineering/particle_merge_implementation.md` | maintainer | Particle merge record | public particle implementation | Checked | Edited | Edited | Added development-record warning | Nietzsche: clean | Verified |
| `source/engineering/particle_merge_unified_data.md` | maintainer | Particle data record | public particle structures | Checked | Edited | Edited | Added development-record warning | Nietzsche: clean | Verified |

## Non-Published Documentation Material

| Location | Classification | Initial action |
| --- | --- | --- |
| `docs/README.md`, `docs/documentation_audit_guide.md`, `docs/documentation_audit_checklist.md`, `docs/documentation_audit_log.md` | Maintainer working material outside Sphinx source | Reconcile with this ledger and final verification evidence. |
| `docs/OLD_possibly_delete/*.md` | Duplicate/legacy copies of several engineering records | Confirm they are not published; do not present them as public guidance. |
| Files existing only under `origin/development:docs/` or its `docs/source/` | Development-only documentation | Do not import into public docs without evidence-backed policy decision. |

## Journey Verification

| Journey | Starting page | Tested path | Findings / edits | Independent review | Final status |
| --- | --- | --- | --- | --- | --- |
| First clone, build, and run | Home | Quickstart -> Building -> Running -> Shock Tube | Corrected build/run/output/restart commands; Sod run and restart executed | Poincare | Verified in render |
| Configure a new simulation | Home | Examples -> Configuration -> Input Parameters | Corrected block parameters, solver constraints, overrides and compatibility gates | Poincare, Anscombe | Verified in render |
| Select a scientific example | Home | Worked Examples -> selected example -> Outputs | Kept runnable public paths; gated MRI2D and binary merger prerequisites | Nash, Anscombe | Verified in render |
| Understand architecture/modules | Home | Overview -> Modules -> selected module | Corrected task/module routing and stable/public scope | Confucius, Copernicus | Verified in render |
| Write/port Kokkos kernels | Home | Kokkos Guide -> Migration -> API/reference | Aligned aliases, callbacks and public porting patterns to code | Confucius, Heisenberg | Verified in render |
| Port Athena++ work | Home | Migration -> Gotchas -> Pgen/Modules | Replaced legacy assumptions with current custom-pgen and Kokkos path | Confucius, Anscombe | Verified in render |
| Find outputs, restarts, algorithms | Home | Reference/Modules -> Outputs and methods pages | Corrected restart/output semantics and deep technical limits | Anscombe, Carver, Halley | Verified in render |
| Locate engineering records safely | Home | Developer Notes -> record -> stable module guidance | Added development/historical warnings and stable links | Copernicus, Nietzsche | Verified in render |

## Reviewer Coverage

| Reviewer | Page scope | Resolution outcome |
| --- | --- | --- |
| Poincare | Beginner start/build/configure/run/troubleshoot path | Cleared after command and output corrections |
| Confucius | Home/navigation, overview/runtime, Kokkos, migration and contributor flow | Cleared after routing and lifecycle corrections |
| Copernicus | Module landing/source terms, visualization, CGM and engineering landing | Cleared after public/development separation |
| Nash | All example pages | Cleared after executable/gated-route corrections |
| Heisenberg | Mesh, Hydro, MHD and API reference | Cleared after interface and compatibility corrections |
| Anscombe | System architecture, parameters/files, outputs, pgen and radiation | Cleared after parser/output/reference corrections |
| Nietzsche | Engineering child records and development routing | Cleared after warning/routing review |
| Carver | Boundaries, coordinates, diffusion, driver, DynGRMHD, EOS and particles | Iterated through findings; final clean report received |
| Halley | TaskList, Z4c, ion-neutral, shearing, reconstruction and solvers | Iterated through findings; final clean report received |

## Confirmed Limitations

| Area | Confirmed public status and documentation disposition |
| --- | --- |
| PDF output | The registered `pdf` output path failed its smoke run with a Kokkos deep-copy extent error; documented as currently failing/unvalidated. |
| MRI2D shearing input | The custom MRI2D run fails because shearing-box source terms are not enabled for that path; documented as non-runnable. |
| Binary merger | Requires external SGRID library and runtime data unavailable in this audit environment; documented as prerequisite-gated, not claimed executed. |
| Radiation beam input | The shipped beam input needs the relevant `pgen_name` selection to run the audited beam path; documented explicitly. |
| Z4c adaptive/boosted decks | Shipped decks contain pgen, extraction-index, AMR-method, and/or missing criterion-bridge defects; documented rather than presented as working examples. |
| Particles | Current public supported guidance is restricted to 3D periodic drift because 2D/non-periodic paths are unsafe or unvalidated. |
| CPU compiler output | Public CPU build succeeds but emits pre-existing AppleClang variable-length-array warnings in `src/mesh/mesh.cpp`; no source change was made in this docs audit. |

## Verification Log

| Gate | Evidence / result | Status |
| --- | --- | --- |
| Strict baseline HTML build | Sphinx read 55 pages; `-W --keep-going` succeeded before edits | Complete |
| Baseline deployed home page visual inspection | Sidebar/search loaded; simulation iframe and canvas present | Complete |
| Baseline local home page visual inspection | Strict local render loaded; simulation iframe and canvas present | Complete |
| Complete published source inventory | 55 of 55 pages included through Sphinx toctrees; no orphaned source page | Complete |
| First-run executable check | Public-baseline build completed Sod; observed `Sod.hydro.hst` and `tab/Sod.hydro_w.*.tab` | Complete |
| VTK example executable check | One-cycle default-build Orszag-Tang run emitted both documented VTK streams | Complete |
| Restart executable check | Added temporary `rst` stream to Sod, observed checkpoints, and resumed from `Sod.00001.rst` | Complete |
| Intermediate strict HTML build after user-flow edits | Sphinx read all 55 pages with `-W --keep-going`; build succeeded | Complete |
| User-flow rendered inspection | Local Quickstart and Running pages rendered with corrected commands and readable code blocks/sidebar | Complete |
| Additional executable checks | Blast, turbulence, MRI3D and ion-neutral one-cycle runs succeeded; binary conversion and fixed radiation beam audit route checked | Complete |
| Deliberate failure/gate checks | MRI2D and PDF failures observed and documented; binary merger dependency absence and malformed Z4c decks dispositioned | Complete |
| Strict final HTML build | `git diff --check` and `make clean html ... SPHINXOPTS='-W --keep-going'` succeeded for all 55 source pages | Complete |
| Final link check | `make linkcheck ... SPHINXOPTS='-W --keep-going'` succeeded; `build/linkcheck/output.txt` has 0 lines | Complete |
| Final rendered landing/deep-page inspection | Browser loaded home, eight journey destinations, and corrected coordinates/TaskList/Z4c/DynGRMHD/engineering pages with expected markers | Complete |
| Interactive banner preservation after edits | Rendered home iframe is present with title; banner iframe contains live `1280 x 427` canvas; `_static/athenak_fluid_sim.html` is unchanged | Complete |
| Reviewer coverage and disposition audit | Every published page has reviewer coverage; all Carver/Halley iterative findings were fixed and re-reviewed clean | Complete |
| Repository scope check | `git diff --exit-code -- src inputs tst vis CMakeLists.txt .github` and banner-asset diff checks passed; only documentation files changed | Complete |
