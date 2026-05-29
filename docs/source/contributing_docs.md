# Contributing Documentation

Published documentation lives under `docs/source/` on the documentation
branch. Public user guidance must describe code, input decks, tests, or
observed behavior present on that public baseline. Branch experiments and
design histories belong in [Developer Notes](engineering/index.md) and must
be labelled as such.

## Source Layout

| Path | Purpose |
| --- | --- |
| `docs/source/index.md` | Home page and top-level navigation |
| `docs/source/modules/` | Implemented module behavior and configuration |
| `docs/source/examples/` | Runnable or explicitly prerequisite-gated cases |
| `docs/source/reference/` | Source-backed reference maps |
| `docs/source/migration/` | Developer porting guidance |
| `docs/source/engineering/` | Qualified implementation records |
| `docs/source/_static/` | Static assets, including the interactive home banner |

The interactive fluid-simulation banner on the home page is a published site
asset. Changes to its integration must keep it visible and functional and
must be checked in a rendered browser view.

## Evidence Standard

For each technical change:

1. Identify the relevant public source file, shipped deck, build option, or
   runtime result.
2. State defaults and limitations only when they are verified from that
   evidence.
3. Mark a workflow as prerequisite-gated or development-only when it is not
   runnable from the public baseline.
4. Add navigation links only after the target page has the same support
   status.

Line numbers are helpful during review, but source file and symbol names are
more durable in published prose.

## Build And Link Check

Install the documentation dependencies, then build with warnings treated as
errors:

```bash
python -m venv .venv-docs
. .venv-docs/bin/activate
python -m pip install -r docs/requirements.txt
sphinx-build -W --keep-going -b html docs/source docs/build/html
sphinx-build -W --keep-going -b linkcheck docs/source docs/build/linkcheck
```

The repository `docs/Makefile` also exposes `html`, `linkcheck`, and `live`
targets:

```bash
make -C docs html SPHINXOPTS="-W --keep-going"
make -C docs linkcheck SPHINXOPTS="-W --keep-going"
```

## Review Checklist

- The page is linked from the intended navigation path or intentionally
  recorded as not published.
- Commands use actual build flags, input paths, and executable behavior.
- Parameters and defaults match their owning source.
- Example output claims were checked with a focused run where feasible.
- Stable pages do not present branch-only work as available functionality.
- Tables, diagrams, code blocks, and the home banner render correctly.
