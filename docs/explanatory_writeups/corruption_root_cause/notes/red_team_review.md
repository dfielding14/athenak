# Red-Team Review

## Verdict

The shared-buffer `Kokkos::realloc` is the leading explanation and is strong
enough to justify rebuilding the PDFs. The draft still overstates the result as
a demonstrated causal proof. The source diff and archive fingerprint are
highly suggestive, but the document does not show a controlled old-code versus
fixed-code reproduction, and its simplified mechanism is incomplete for
`output29`.

## Prioritized Findings

### P0: The stated mechanism does not fully explain why `output29` loses angle

The cartoon says the destructive `coord_abscostheta` call erases earlier rows
and then writes a surviving angle row. That explains radius-first streams, but
not `output29`, whose axis order is `coord_abscostheta` followed by `coord_r`.
After the inner realloc, the angle is written first and radius is written
later. On the simple story in the draft, both should survive.

The missing step is visible in the old source. The inner realloc changes the
first extent from `nmb_max` to `nmb`. A later
`ComputeDerivedVariable(coord_r)` call can then trigger the common capacity
check and realloc back to `nmb_max`, erasing the angle before writing radius.
Whether that second realloc occurs depends on the MeshBlock-pack extents. The
draft needs this two-stage sequence, including the relevant `nmb` versus
`nmb_max` condition. Without it, the central mixed-survival explanation is
internally inconsistent.

**What would make the conclusion wrong:** instrumenting the broken executable
shows that no second realloc occurs for `output29`, or shows that the angle row
survives the full derived-variable sequence while the archived angle still
collapses.

### P0: “Demonstrated root cause” is stronger than the presented causal test

The evidence is a source diff plus an archive pattern plus agreement from a new
reducer. That is a strong diagnosis, but it is not the clean causal
intervention the status label implies. There is no shown regression test that
runs the same AthenaK PDF calculation on the same state with the offending line
present and absent and reproduces/removes the exact mixed pattern.

The eight-shard rebuild establishes that saved `hydro_w` can reproduce the
preserved marginal and `output31`. It does not directly establish that the
specific old online writer operation caused every archived discrepancy.

**What would make the conclusion wrong:** the old-line/new-line A/B test fails
to create and remove the observed radius-first collapse and `output29`
angle-loss pattern. A different PDF-writer bug that survives removal of this
line would also invalidate the strong version.

### P1: The archive-wide impact claim is not traceable to a completed audit artifact

The draft says twelve of fourteen streams are unusable across all five
simulations. The available corruption summary says the latest snapshots were
checked, while the all-snapshot audit script documents an intended full scan.
No generated full-audit CSV or Markdown result is present under the documented
output directory. The draft therefore blurs three different claims:

- the stream definitions are structurally vulnerable;
- the latest snapshots in five simulations exhibit the pattern;
- every archived time sample in those streams is unusable.

Those are not equivalent. This matters because an early output written before
the regression, or under a different pack state, could still be scientifically
salvageable.

**What would make the conclusion wrong:** a full archive audit finds any
radius-first snapshot with meaningful non-underflow radial support, or finds
time ranges written before the bad commit that do not show the corruption.

### P1: `output31` is a useful control, but not an independent ground truth

`output31` avoids `coord_abscostheta`, so it is a good negative control for
this specific bug. It still comes from the same AthenaK PDF machinery, uses the
same archived sparse format, and is compared against a reducer implementing
similar coordinate and energy-flux definitions. Calling it an “independent
control” invites a stronger interpretation than the evidence supports.

The full-4D mismatch of `3.87e-7` is small and consistent with float32 saved
fields versus online precision, but that explanation has not been isolated
with a precision-controlled experiment.

**What would make the conclusion wrong:** an independently calculated
cell-level `output31` reference disagrees materially with both the archive and
the reducer, or the observed agreement is shown to result from a shared
definition or reader error.

### P1: The draft treats zeroed coordinates as universal, but `realloc` behavior and pack state matter

The text repeatedly says rows were “erased / zero” and cells were counted using
a coordinate “replaced by zeros.” The observed radius-underflow signature
supports zero-like values in the affected archive, but the source-level
guarantee is weaker: `realloc` does not preserve contents. The exact result can
depend on Kokkos initialization behavior, backend, allocation size, and whether
the requested extents changed. The explanation should distinguish the unsafe
operation from the empirically observed zero/underflow outcome.

**What would make the conclusion wrong:** reproducing the old code on the
production backend yields preserved or nonzero stale rows rather than the
observed deterministic underflow pattern.

### P2: The figure is a strong fingerprint but is not unique to this mechanism

The matched eight-shard figure shows that total and radial marginal agree while
angle does not. That rules out many broad failures. It does not uniquely select
the shared-buffer realloc: any online-writer defect that corrupts only the
angular bin map after retaining the same weights could produce the same
high-level figure. The source diff makes the realloc explanation compelling,
but the figure alone is not a causal discriminator.

**What would make the conclusion wrong:** another angular-coordinate defect,
axis-stride defect, or stale-row indexing defect reproduces the same controls
and remains after the realloc is removed.

## Required Changes Before Calling It Demonstrated

1. Show the full two-reallocation timeline for radius-first products and
   `output29`, including `nmb` and `nmb_max`.
2. Add a same-state broken-line/fixed-line AthenaK regression result, or weaken
   the status to “strongly supported root cause.”
3. State the exact temporal coverage of the completed archive audit and link
   its generated result.
4. Describe `output31` as a negative control for `coord_abscostheta`, not fully
   independent ground truth.
5. Separate “unsafe contents are not preserved” from the observed
   “affected archived values land at zero/underflow.”
