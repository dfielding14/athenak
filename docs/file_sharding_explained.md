# File Sharding Explained

The GitHub Pages-ready version of this documentation lives at
[`docs/source/modules/output_file_sharding.md`](source/modules/output_file_sharding.md).
That file is written in MyST Markdown and matches the `origin/gh-pages`
documentation layout.

To integrate it into the published documentation branch:

1. Copy `docs/source/modules/output_file_sharding.md` into the same path on
   `gh-pages`.
2. Add `modules/output_file_sharding` to the hidden modules toctree in
   `docs/source/index.md`.
3. Add a row or link for file sharding in `docs/source/modules/index.md`.
4. Optionally link it from the Outputs page near the output-parameter table.

The page documents the three file-layout modes, the restart manifest/payload
layout, supported output types, reader behavior, and operational guidance for
choosing `shared`, `single_file_per_rank`, or `single_file_per_node`.
