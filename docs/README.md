# AthenaK Documentation Guide

Welcome to the documentation toolkit for AthenaK. This folder now contains a conventional Sphinx project that is easy to extend and rebuild locally.

## Quick Start

```bash
cd docs
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
make html
```

The generated site will live in `docs/build/html/index.html`.

### Live Preview

During authoring you can run:

```bash
make live
```

This uses `sphinx-autobuild` to rebuild the site and refresh the browser automatically (served at `http://127.0.0.1:8000` by default).

## Writing Content

- Place Markdown (`.md`) or reStructuredText (`.rst`) files under `docs/source`. The build is configured with [MyST Markdown](https://myst-parser.readthedocs.io/en/latest/) so you can use familiar Markdown syntax alongside Sphinx directives.
- Add new documents to the appropriate `toctree` in `docs/source/index.md` (or the relevant section index) so they appear in the navigation.
- Reuse the existing section structure:
  - `examples/`
  - `flowcharts/`
  - `migration/`
  - `modules/`
  - `reference/`
- Keep large diagrams or static assets in `docs/source/_static`.

## Style & Best Practices

1. Keep sentences short and prefer active voice.
2. Use fenced code blocks with language clues (e.g. ```cpp).
3. Reference source code with relative paths (`src/module/file.cpp:42`) when explaining behaviour.
4. Run `make linkcheck` in CI or locally for new pages that add external links.
5. Capture any inconsistencies you find between docs and code in `documentation_audit_log.md`.

## Troubleshooting

- **Missing theme or parser**: ensure dependencies from `requirements.txt` are installed (consider `pip install -r docs/requirements.txt --upgrade`).
- **Cross-link warnings**: the build runs in `nitpicky` mode, so reference targets must exist or be added to `nitpick_ignore` in `conf.py`.
- **Python path imports**: `conf.py` prepends the repository root to `sys.path` for autodoc; adjust if you relocate modules.

Happy documenting!
