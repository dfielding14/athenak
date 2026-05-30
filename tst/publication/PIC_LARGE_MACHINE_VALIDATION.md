# RETIRED: Historical PIC Large-Machine Validation Note

> **Historical evidence only. Do not use this file as a release gate, scientific
> qualification plan, Frontier runbook, or authorization to submit jobs.**

The canonical and controlling document is
`tst/publication/PIC_PRODUCTION_READINESS_PLAN.md`. It supersedes the scope,
terminology, pass criteria, operational instructions, commands, and sign-off
language that previously appeared here.

## Historical Context

The removed runbook described an exploratory workflow built around an obsolete
`c/pic-review` branch and a historical dirty-tree snapshot. It recorded:

- a successful local default-particle regression run;
- a successful exploratory proxy-manifest run;
- F01-F09 engineering metric extraction and quick-look plotting;
- proposed larger-machine exploratory proxy, plotting, and shock-rich
  Orszag-Tang stress workflows.

Those observations remain useful provenance for why the broader readiness plan
was written. They are not current evidence for physical correctness, Sun and
Bai (2023) reproduction, production readiness, or a scoped state-of-the-art
claim.

## Historical Artifact Locations

The old workflow wrote local exploratory artifacts below:

- `tst/.codex/pic_publication_runs/`
- `tst/build/src/bin/`
- `tst/build/src/pvtk/`

Historical names containing `publication` are compatibility identifiers only.
They do not assign an evidence class. Current artifact tooling must emit an
explicit evidence class and `not_sun_bai_reproduction=true` for unqualified
engineering artifacts.

## Current Procedure

Do not recover or execute commands from prior revisions of this file. Start
with `tst/publication/PIC_PRODUCTION_READINESS_PLAN.md`, register the intended
claim, verify its prerequisites, use only the authorized Frontier phase, and
record node-hour reservations and reconciled usage in the required ledger.
