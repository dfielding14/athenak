# Visual Diagram Plan

## Visual Language

Blue means trusted or preserved, red means destructive or broken, amber means
partial, and gray means infrastructure. Every diagram is a schematic, not
evidence.

## Opening Cartoon

Show two shared-buffer timelines. In radius-first streams, the inner
`coord_abscostheta` realloc invalidates radius. In `output29`, the inner realloc
writes angle into an `nmb` buffer, then the next call's common capacity check
reallocates back to `nmb_max`, invalidates angle, and writes radius. This
two-reallocation sequence is the key to the mixed outcome.

Prompt: `diagram_prompts/shared_buffer_erasure.md`.

## Summary Diagram

Place the fourteen streams on a short timeline split by coordinate order:
twelve radius-first streams collapse, `output29` loses angle but recomputes
radius afterward, and `output31` avoids the destructive coordinate entirely.

Prompt: `diagram_prompts/mixed_survival_timeline.md`.

## Evidence Pair

Use `output29_corruption_signature.pdf` as the mechanism fingerprint and
`archived_control_errors.pdf` as the strong negative-control summary.
