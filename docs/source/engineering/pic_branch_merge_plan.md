# MHD-PIC Breaking Merge Plan

## Goal

Eventually merge `PIC_development` into `PIC` as the single supported MHD-PIC
implementation. This will be a deliberate breaking cutover, not a compatibility
release. The merged branch will support the explicit model we intend to use for
new science runs; it will not carry parsers, aliases, restart fields, tests, or
control-plane machinery whose only purpose is to reproduce older development
interfaces.

Old results remain reproducible through immutable Git tags and their archived
executables. In particular, the completed non-relativistic shock is tied to
`pic-nonrelshock-t3000-20260720`, and the cleaned uniform-grid implementation is
tied to `pic-uniform-science-candidate-20260720`. Post-merge executables will not
be required to restart checkpoints made by those older interfaces.

## Compatibility Policy

The merge intentionally provides no backward compatibility for:

- inputs using the former `pic_physical_mode` selector or its implicit defaults;
- restart files whose particle section predates the post-merge schema;
- retired Hall experiments, raw-current-to-CT controls, or ambiguous delta-f
  aliases;
- superseded generated campaign matrices and their frozen launch machinery; or
- historical analyzers whose contracts are bound to removed source layouts.

We will not add an input translator, dual restart reader, compatibility mode, or
deprecated alias layer. Anyone needing an old checkpoint must use the tagged
source and executable that created it. Existing output files remain available
for analysis independent of restart support.

## Cutover Sequence

1. **Choose the cutover point.** Finish or freeze any active run that still
   needs an old checkpoint. Record the exact commit, deck, executable, and
   analysis entry point for every retained result.
2. **Tag both sides.** Create a final pre-merge tag on `PIC` and confirm the
   current `PIC_development` science-candidate tag. These tags are the complete
   backward-compatibility mechanism.
3. **Create a temporary integration branch.** Merge `PIC_development` into the
   latest `PIC`, resolve genuine branch conflicts, and make no physics changes
   while resolving them.
4. **Remove compatibility code in one reviewable commit.** Delete the
   restart-only obsolete-selector shim, `RestartLegacyModeTag()`,
   `MatchesRestartLegacyModeTag()`, the schema-8 legacy mode slot, and any tests
   that exist only to read the old layout. Introduce a compact new restart
   schema containing only the state needed by the explicit runtime model.
5. **Remove obsolete repository baggage.** Delete remaining frozen scripts,
   readiness bindings, and documentation that describe unsupported launch
   paths. Keep only current decks, analyzers, physical references, and concise
   release notes pointing to the archival tags.
6. **Validate the clean interface.** Run the normal-precision focused
   conservation test, one MPI single-precision turbulence smoke, one short
   full-Hall shock canary, and a new-schema uninterrupted-versus-restart
   comparison. Test only the retained interface; do not test old schemas or
   aliases.
7. **Review and merge.** Inspect the final diff for physics changes, confirm the
   retained science decks are explicit and reproducible, then merge into `PIC`
   only with project-owner approval.

## What the Final Interface Keeps

The merged code should expose the controls that directly select an algorithm or
physical term:

- `time/integrator=vl2` for the staged MHD-PIC chronology;
- explicit particle pusher, artificial light speed, and initializer semantics;
- explicit deposition, background, feedback, and gas-coupling controls; and
- `pic_cr_hall_mode=off|full` for the physical CR-Hall closure.

Defaults should be simple and safe, while every retained science deck states its
full model explicitly. Descriptive names such as `paper_smooth` may remain when
they identify a real numerical algorithm; names and code that exist solely for
legacy compatibility should not.

## Definition of Done

The breaking merge is ready when:

- active source and inputs contain no compatibility selector or legacy mode tag;
- the restart layout contains no unused compatibility fields and round-trips
  only the new schema;
- all retained decks use the explicit interface;
- old campaign machinery is recoverable only through pre-merge tags, not the
  active tree;
- the compact validation set passes in the supported Frontier configurations;
  and
- `PIC` and `PIC_development` have a documented, reviewable merge history.

