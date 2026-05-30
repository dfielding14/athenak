# IO Feature Branch Process Artifact Disposition

## Purpose

This record applies the `RCP-10` process-packaging rule without deleting useful
history prematurely. The branch still has open external CUDA and scheduler-backed
multi-node gates, so the large execution records remain available for review until
those gates close.

Do not remove, relocate, or condense a listed artifact silently. Any later cleanup
commit must update this table and identify where the durable content moved.

## Current Disposition

| Artifact | Current disposition | Final merge disposition | Reason and migration rule |
| --- | --- | --- | --- |
| `IO_FORMAT_COMPATIBILITY.md` | Retain | Retain in the final code branch | This is the durable public compatibility contract for legacy and new layouts. Update it when qualification evidence changes. |
| `IO_FEATURE_BRANCH_DECISION_LOG.md` | Retain during review | Summarize durable decisions into the compatibility contract and user-facing docs; preserve the full log in the pull-request record unless maintainers explicitly request an in-repository development archive | The alternatives and reversal paths are valuable during review, but the full implementation diary is larger than the long-term runtime contract. |
| `IO_FEATURE_AUDIT_LEDGER.md` | Retain during review | Preserve the final evidence summary in the pull-request record; archive in repository only if maintainers want a durable development record | The ledger is the authoritative checkpoint and residual-risk record while the feature is being qualified. |
| `IO_FEATURE_BRANCH_GUIDE.md` | Retain during review | Preserve in the pull-request record or an explicitly approved development archive | It records the original feature-isolation plan and should not become ordinary user documentation. |
| `IO_FEATURE_BRANCH_INTEGRATION_PLAN.md` | Retain during review | Preserve in the pull-request record or an explicitly approved development archive | It explains selective integration from the historical branches and remains useful for reviewer provenance. |
| `IO_FEATURE_BRANCH_ROBUSTIFICATION_GUIDE.md` | Retain at its frozen root path until robustification and external gates close | Preserve in the pull-request record or an explicitly approved development archive after final checksum verification | The frozen guide is still an active execution contract. Moving it before closure would invalidate the recorded path and complicate audit comparisons. |
| `IO_RESTART_MANIFEST_SCALING_QUALIFICATION.md` | Retain | Retain until `RCP-08B`; then keep the completed evidence summary and move raw scheduler logs outside the repository | This is the preregistered scaling contract needed to interpret external timing evidence. |
| `IO_EXTERNAL_IO_QUALIFICATION_PLAN.md` | Retain | Retain until external CUDA and multi-node qualification close, plus a HIP packet or explicit `N/A` disposition for the intended deployment backend; then summarize results in the pull request and compatibility contract | This prevents local workstation evidence from being overstated as deployment qualification. |
| `deferred_docs/gh-pages/io-output-formats-and-sharding/` | Retain | Retain until the separate post-code-merge `gh-pages` change merges | Public documentation must remain staged and unpublished on the code branch. |
| `scripts/stage_gh_pages_io_docs.py` | Retain | Retain with the deferred Pages bundle while it remains useful for deterministic restaging | The helper is part of the safe two-phase Pages lifecycle, not a publication side effect. |
| `scripts/run_external_io_qualification_slurm.sh` | Retain | Retain until external CUDA, multi-node, filesystem, and `RCP-08B` evidence closes; then preserve or archive with the qualification record | The executable Slurm runner keeps preregistered topology, repetition, timeout, log, and filesystem-accounting requirements reproducible. |
| `scripts/finalize_external_io_qualification_packet.sh` | Retain | Retain with the external qualification runner while immutable packet evidence remains required | The finalizer makes the acyclic inner-manifest and outer-index checksum contract executable, publishes metadata through synced packet-local atomic renames, re-syncs admitted retained metadata, recovers writable inconsistent metadata pairs, removes packet write permissions before archival recording, and admits the archive sink used by its final substantive append. Retained packet-local reserved paths and archive-adjacent destination-scoped replacement candidates fail closed for operator adjudication. Use it within the documented trusted lifecycle from the first runner invocation through finalization without inter-invocation packet-child replacement. |
| `scripts/validate_external_io_packet_index.py` | Retain | Retain with the external qualification runner and packet finalizer | The shared validator enforces strict runner-index lifecycle and canonical outer-archive TSV grammars, plus descriptor-safe metadata-pair reset, synced local metadata publication, checked complete sibling-temporary atomic replacement for packet-index and admitted outer-archive publication, and fail-closed archive-directory admission when a destination-scoped replacement candidate remains. |

## Prohibited Cleanup

Do not:

1. delete the frozen robustification guide while its checksum is still a checkpoint
   gate;
2. remove the ledger before external qualification findings and residual risks have
   a durable replacement;
3. publish deferred Pages content from the code feature branch;
4. hide qualification gaps by condensing them into a generic follow-up sentence; or
5. move raw scheduler logs into the source repository when a retained evidence path
   and summarized result are sufficient.

## Final Packaging Reflection

The branch remains one coherent IO feature: output distribution, format readers,
runtime policy, tests, deferred Pages content, and qualification contracts all serve
the same user-facing workflows. The process records are intentionally verbose
because the branch is still under qualification. Their later reduction is a
review-stage packaging step, not a reason to discard evidence during implementation.
