# Phase-9 Portable All-Format Replay

Extract `phase9_all_formats_portable_replay_bundle.tar.gz` into an arbitrary
directory and run:

```bash
python3 phase9_portable_replay.py
```

The package carries:

- the qualified release AthenaK binary;
- the retained `acceleration_uninterrupted` runtime artifacts;
- the Phase-7 diagnostics input and preregistered criteria;
- the strict all-format inspector and its two Python dependencies;
- a diagnostics metrics template with the qualified binary, input, and
  criteria SHA-256 digests;
- `portable_manifest.sha256`.

The wrapper verifies the package manifest, rewrites extraction-local paths in
the metrics template, invokes the bundled strict inspector, and writes
`portable_all_formats_report.json`.  The emitted report normalizes the
extraction prefix to `$PACKAGE_ROOT`, so two fresh extractions produce
byte-identical reports.
