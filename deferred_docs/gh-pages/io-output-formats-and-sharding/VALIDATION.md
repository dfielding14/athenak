# Deferred Documentation Application And Validation

## Do Not Publish Before Code Merge

This procedure is for applying the staged pages to a separate `gh-pages`
change after the IO code branch has merged. CP-03 native direct restart loading
has replaced the earlier transient rank-0 `<manifest>.assembled` staging
design. Building pages in a temporary worktree is permitted before
publication; pushing or merging Pages content before the code exists on the
primary branch is not.

## Integration Procedure

1. Re-check the merged code rather than trusting this staged package alone:

   ```bash
   rg -n "single_file_per_node|final_output_policy|output_timing|sphslice|scale[1-4]|weight_variable" \
     src inputs/io tst/test_suite/io vis/python
   rg -n "bin_convert_new" vis/python inputs tst
   rg -n "assembled|CanonicalPayloadPath|LoadLocalBlocks|Read_bytes_at_all|number of nodes|sparse_angles" \
     src tst/test_suite/io vis/python
   ```

   The second search should show only intentional rejection/tests or migration
   discussion, not an imported public module or promoted example. For the
   third search, `.assembled` references should be no-production-staging tests
   or clearly marked historical documentation, not a runtime assembly path.

2. Create a temporary worktree for the current Pages baseline:

   ```bash
   git fetch origin gh-pages
   git worktree add --detach /tmp/athenak-gh-pages-io-docs origin/gh-pages
   ```

3. Copy the complete staged pages from `overlay/docs/source/` to their matching
   locations under `/tmp/athenak-gh-pages-io-docs/docs/source/`.

4. Merge the two reviewed insertion fragments into the corresponding live
   reference pages:

   - add `input_parameters.io_outputs.md` material to
     `docs/source/reference/input_parameters.md`;
   - add `file_reference.io_outputs.md` material to
     `docs/source/reference/file_reference.md`.

   Retain all unrelated live reference material. These are insertions, not
   whole-page replacements.

   For a validation-only build before editing the live page bodies, copy each
   fragment beside its target using an `.inc` suffix and temporarily append a
   MyST ``{include}`` directive in the parent page. Do not copy a section-level
   fragment into `docs/source/` as a standalone `.md` page because Sphinx will
   correctly warn that it is not a complete top-level document.

5. Verify navigation:

   - `docs/source/modules/index.md` still routes to `modules/outputs`;
   - `docs/source/index.md` still routes to the tools and examples catalogues;
   - `docs/source/examples/index.md` now routes to
     `examples/io_outputs_and_sharding`.

6. Audit contradictory language in the temporary Pages worktree:

   ```bash
   rg -n "bin_convert_new|logscale[1-4]|mass_weighted is not|unversioned|single_file_per_node|sphslice|final_output_policy|output_timing" \
     /tmp/athenak-gh-pages-io-docs/docs/source
   ```

   Resolve stale converter/PDF descriptions before the Pages change is
   proposed. References to new names should agree with the merged source and
   executable examples.

7. Build the candidate site with warnings treated as errors:

   ```bash
   cd /tmp/athenak-gh-pages-io-docs/docs
   make clean html SPHINXOPTS="-W --keep-going"
   ```

8. Re-run the IO tests or cite the accepted code-branch qualification result,
   with particular attention to public example input decks and reader
   commands.

9. Assign a read-only independent audit comparing:

   - parser keys and defaults in `src/outputs/outputs.cpp` and
     `src/driver/driver.cpp`;
   - layouts implemented by output writers and native direct restart loading;
   - reader behavior in `vis/python/`;
   - test/example inputs; and
   - candidate Pages content and toctrees.

10. Only after discrepancies are resolved or explicitly scoped out should a
    separate `gh-pages` review be opened.

## Expected Validation Record

The Pages change description should include:

- code merge commit and Pages base commit;
- replacement pages and inserted sections;
- successful `make clean html SPHINXOPTS="-W --keep-going"` result;
- confirmation that public documentation no longer promotes
  `bin_convert_new.py`;
- confirmation that node restart uses the public manifest only, performs native
  direct payload reads, and has no production `.assembled` staging path;
- frozen `origin/main` shared and per-rank restart resume qualification status,
  in addition to immutable fixture checksum status;
- CUDA-capable GPU execution status for the IO regression, distinguishing a
  real device run from the CPU smoke path;
- real multi-node qualification status for per-node output/restart paths; and
- any retained boundary, especially sliced `cbin`.
