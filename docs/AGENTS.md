<!-- BEGIN build-memory-table -->
# Documentation Navigation

## Scope

`docs/` contains the Sphinx manual, engineering contracts, audit records, and
historical/reference material. Documentation explains the repository but does
not override runtime behavior: verify technical claims in `src/`, the matching
input decks, and focused tests before recording them here.

## Task lookup

| Task | Start in | Also inspect or validate |
| --- | --- | --- |
| Change site navigation or Sphinx configuration | `source/index.md`, `source/conf.py` | The affected toctree, links, and a local HTML build |
| Change a module, example, migration, or tool page | `source/modules/`, `source/examples/`, `source/migration/`, or `source/tools/` | Matching source guide and implementation under `src/` |
| Change an MHD-PIC model or implementation contract | `source/engineering/pic_mhd_model_contract.md`, `source/engineering/pic_cr_hall_code_map.md` | `MHD_PIC_NEXT_STEPS_GUIDE.md`, `src/particles/`, `src/mhd/`, and focused regressions |
| Audit documentation against code | `AGENT_PRIMER.md` | `documentation_audit_guide.md`, `documentation_audit_checklist.md`, and `documentation_audit_log.md` |
| Change the PIC manuscript draft | `pic_test_problem_description_draft/` | Read its local `AGENTS.md`; keep manuscript and implementation evidence distinct |
| Change static assets or theme behavior | `source/_static/`, `source/_templates/`, `source/conf.py` | Render the pages that consume the asset |

## Build and validation flow

- `source/` is the authored Sphinx tree; `build/` is generated output and must
  not be hand-edited or treated as source.
- Keep new pages reachable from an appropriate toctree in `source/index.md` or
  a child index.
- Build HTML from the repository root with `make -C docs html`.
- Check external and internal links with `make -C docs linkcheck` when links or
  navigation change.
- For iterative authoring, use `make -C docs live` if the optional
  `sphinx-autobuild` dependency is installed.

## Local constraints

- Follow `AGENT_PRIMER.md` for evidence-backed audits; do not infer current
  behavior from stale prose.
- `reference_paper/` is imported reference material, and
  `OLD_possibly_delete/` is historical. Do not modernize either tree as part of
  ordinary manual edits.
- Preserve scientific notation, equations, and provenance when updating model
  documents. State explicitly whether a test is a regression, engineering
  proxy, preparation, or qualification result.
<!-- END build-memory-table -->
