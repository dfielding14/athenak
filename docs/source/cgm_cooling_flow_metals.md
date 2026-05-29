# Development Record: CGM Cooling Flow With Metals And Supernovae

```{warning}
This is not a runnable public AthenaK example. The public implementation
baseline used by this documentation site does not contain
`src/pgen/cgm_cooling_flow_amr_metals.cpp`, its supporting profile/SN helpers,
or a shipped CGM cooling-flow input deck.
```

## Status

An earlier version of this page described a CGM cooling-flow problem generator
with metallicity, supernova feedback, and additional cooling support as if it
were available in the published implementation. Audit comparison against the
public `main` implementation shows those files and inputs are absent.

The page is retained only as a record that such feature work has been
investigated. It must not be used for public run commands, parameter lookup, or
claims about supported physics.

## Supported Public Alternatives

- [Source Terms](modules/srcterms.md) documents the public source-term
  implementation, including the shipped ISM cooling function and relativistic
  cooling flag.
- [Worked Examples](examples/index.md) lists runnable public example paths.
- [Configuration](configuration.md) documents supported input-file structure
  and output streams.

## Publication Gate For Future CGM Guidance

Before a CGM cooling-flow workflow can appear as a stable example, its owning
branch must provide and validate all of the following on the public
implementation baseline:

1. Merged problem-generator and supporting source files.
2. Shipped input decks with documented required external data, if any.
3. Verified build and run commands with observable output and restart behavior.
4. Source-backed parameter documentation and representative validation tests.

Until those gates are met, CGM-specific configurations and physics claims
belong in development records rather than the stable user path.
