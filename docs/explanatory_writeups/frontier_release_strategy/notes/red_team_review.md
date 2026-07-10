# Red-Team Review

## Verdict

Golden-snapshot-first is the right release strategy. The current draft
describes a stronger gate than the job system actually enforces. As written,
an operator can submit production without a saved golden-validator pass, the
golden job does not run the full validator automatically, and routine
production jobs do not execute the archived controls the document says remain
required.

## Prioritized Findings

### P0: The “golden gate” is advisory, not an enforced production dependency

`frontier_golden_8pc.sbatch` runs the reducer and the lightweight manifest
check, then only prints an instruction to run `validate_real_8pc.py`.
`submit_frontier_production.sh` requires a generic confirmation token but does
not require, locate, or verify a passing golden JSON report. Nothing prevents
production submission before the golden run or after a failed golden
validator.

The document predicts that a failed golden run “should block every production
array.” The current tooling does not implement that property.

**What would make the conclusion wrong:** successfully submit an executing
production array with the confirmation token when no passing golden report
exists. The current scripts appear to allow exactly that.

### P0: The release artifact is not pinned to an executable or source revision

The jobs accept any executable at the configured path, and the rebuild
manifest records source data metadata but no reducer commit, source-tree
status, executable checksum, compiler/toolchain identity, or product-catalog
version. The current worktree is also uncommitted/untracked. A golden pass can
therefore validate one binary while production silently runs another.

For a one-off campaign, this is a more plausible common-mode failure than many
of the risks in the narrative.

**What would make the conclusion wrong:** replace or rebuild the executable
after the golden pass and show that production submission proceeds without
detecting the change.

### P0: The global AMR “proof” is weaker than claimed

The complete run activates unique logical-key, ancestor/descendant, and summed
volume checks. Those checks do not prove exact spatial coverage because the
reducer does not verify logical-key bounds, geometry-to-key consistency,
same-level physical overlap, or holes offset by extra equal-volume blocks.
The draft says the complete run can “prove” global AMR uniqueness and closure;
it can only pass the implemented invariants.

This weakness matters more on the 21 snapshots without archived PDF controls,
where internal invariants carry most of the release argument.

**What would make the conclusion wrong:** a synthetic malformed full snapshot
with compensating hole and overlap passes the golden global checks.

### P1: Per-snapshot archived controls are promised but not implemented in production

The document says per-snapshot manifest validation and archived controls remain
required. The production worker only runs `validate_rebuild_manifest`, which
checks counts, finite sums, and file existence. It logs
`control_pdf_sequence` but does not compare against it. The full
`validate_real_8pc.py` is hard-coded to the one golden identity.

Consequently, the 140 mapped controls are inventory metadata, not an active
production safety net.

**What would make the conclusion wrong:** a production output with wrong bin
contents but finite totals and complete files passes the worker’s post-run
validation. The current lightweight validator would permit that.

### P1: One 8 pc golden snapshot does not cover the campaign’s main variants

The selected golden snapshot covers one simulation, one phase, 1024 ranks, and
128 nodes. It does not cover:

- 4 pc snapshots at 256 nodes;
- uniform-phase snapshots with no archived controls;
- high- and low-accretion variants;
- maximum block counts or worst filesystem pressure;
- every source-header/unit variant present in the five manifests.

The draft mentions a 256-node canary, but it is not part of the formal release
requirements or the submission guard. A single success should release only the
validated class, not automatically all five simulations.

**What would make the conclusion wrong:** the 8 pc golden passes while the
first 4 pc, uniform-phase, or high/low-accretion snapshot fails due to a
variant-specific header, scale, AMR, or resource issue.

### P1: “Acceptable throughput” is undefined and therefore cannot gate release

The formal list requires acceptable `original`, `science`, and `all`
throughput, but gives no numeric walltime, cost, memory, I/O, or scaling
threshold. The golden job requests two hours; production workers request six.
Without a predeclared threshold, a slow result can be rationalized after the
fact, and the release decision is not reproducible.

The validation thresholds are also proposed after partial results are known.
They may be reasonable, but the draft does not justify them from precision
loss, reduction order, or a held-out comparison.

**What would make the conclusion wrong:** reviewers given the same golden logs
reach different release decisions because no quantitative performance
criterion exists.

### P1: The golden external controls cover a narrow scientific slice

The full validator compares only `output29` and `output31`, both using
`edot_sph`. That is a good repair-path check, but it does not independently
validate the mass, volume, scalar, metal, dust, and other science products.
The golden strategy currently treats correctness and throughput of
`science`/`all` as if the same gate covers both; in practice, only throughput
would be newly observed for many products.

**What would make the conclusion wrong:** `original` controls pass and
`science`/`all` run quickly, but an independently checked science product uses
the wrong physical conversion or weight.

### P2: “Atomic publication” and output refusal do not prevent incomplete visible products

Product files are individually renamed before the final manifest exists.
Refusing an existing output directory prevents accidental overwrite, but a
failed job can leave final-looking PDFs in a claimed directory. The gate is
safe only if every consumer and release script requires the final manifest and
validator report.

**What would make the conclusion wrong:** a downstream workflow discovers and
uses final-named PDFs from a failed golden or production directory without
checking the manifest/report.

## Required Changes Before Calling This a Release Gate

1. Make production submission require a passing, immutable golden report tied
   to the exact executable checksum and source revision.
2. Run `validate_real_8pc.py` inside the golden job and fail the job on any
   validator failure.
3. Add geometry-to-logical-key/domain checks before saying global coverage is
   proved.
4. Implement generic per-snapshot archived-control validation for the 140
   mapped snapshots.
5. Make the 256-node 4 pc canary and at least one no-control uniform snapshot
   explicit release stages.
6. Define numeric throughput, walltime, memory, and validation thresholds
   before observing the Frontier result.
