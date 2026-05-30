# Deferred Documentation Application And Validation

## Do Not Publish Before Code Merge

This procedure supports a detached pre-merge preview and a post-code-merge
restaging pass. Building pages in a detached worktree is permitted before
publication. Pushing or merging Pages content before the code exists on the
primary branch is not.

## Strict Preview

Create a detached worktree from the inspected Pages baseline:

```bash
git worktree add --detach /tmp/athenak-gh-pages-io-docs origin/gh-pages
python3 scripts/stage_gh_pages_io_docs.py /tmp/athenak-gh-pages-io-docs
python3 scripts/stage_gh_pages_io_docs.py \
  --verify-staged /tmp/athenak-gh-pages-io-docs
```

The helper refuses a symbolic branch target, a dirty strict target, baseline
or protected-blob drift, source-fragment drift, missing or ambiguous anchors,
unexpected add-target presence, non-idempotent transformations, contradiction
search failures, and any modified-file set other than the manifest allowlist.

It prints a complete review diff and these strict documentation commands:

```bash
cd /tmp/athenak-gh-pages-io-docs/docs
make clean html SPHINXOPTS="-W --keep-going"
make linkcheck SPHINXOPTS="-W --keep-going"
```

Pass `--run-builds` to either strict staging or `--verify-staged` when the
helper should run those two commands itself.

## Drift Reconciliation

If the Pages baseline has moved, do not edit around the failure manually.
Create a fresh detached worktree from the new baseline and generate an
external packet:

```bash
python3 scripts/stage_gh_pages_io_docs.py \
  --reviewed-drift /tmp/athenak-gh-pages-io-docs-drift \
  /tmp/athenak-gh-pages-io-docs
```

Reviewed-drift mode writes `base/`, `current/`, `proposed/`, `diffs/`, source
payload copies, and `report.json` outside the target worktree. It fingerprints
the target before and after packet creation and fails if target bytes or Git
status change. An auditor must review the packet and approve a reconciled
manifest baseline before a later strict staging pass writes anything.

## Post-Code-Merge Restaging

After the IO code branch merges:

1. Update local `origin/gh-pages` explicitly outside the helper.
2. Create a fresh detached worktree from that current remote-tracking ref.
3. Run strict staging. If the baseline moved, use reviewed-drift mode and
   reconcile deliberately.
4. Run `--verify-staged`.
5. Run warnings-as-errors HTML and link-check builds.
6. Review the exact nine-file diff and the contradiction-search report.
7. Inspect the rendered Modules Support Systems table. Confirm that Outputs and
   Boundary Values remain rows in one table and that both links resolve.
8. Follow the rendered route from Examples to IO Outputs And Sharding, inspect
   the PDF implementation table, verify that its PDF entry remains a table row,
   verify the `payload_rank` fail-closed prose, and confirm that the home-page
   iframe remains present.
9. Assign an independent code-to-doc audit.
10. Open a separate `gh-pages` review only after discrepancies are resolved or
   explicitly scoped out.

The helper never fetches, commits, pushes, or edits a checked-out `gh-pages`
branch.

## Expected Validation Record

Record:

- merged code commit and Pages baseline commit;
- helper command and exact allowlist diff;
- successful `--verify-staged` result;
- warnings-as-errors HTML and link-check results;
- rendered Modules-table, neighboring-link, example-route, PDF-table-entry,
  `payload_rank`, and home-iframe checks;
- contradiction-search report;
- independent code-to-doc and navigation audit dispositions;
- real CUDA-capable GPU execution status, distinguishing it from CPU smoke;
- intended deployment backend and separate HIP packet status or explicit `N/A`
  justification;
- real multi-node qualification status for per-node output and restart paths;
  and
- retained boundaries, especially sliced or adaptive `cbin`.
