# Gotham IO Feature Branch Guide

## Status And Purpose

This document defines the proposed scope, behavior, implementation approach, and
validation plan for a new AthenaK IO feature branch. The branch will isolate reusable
IO functionality developed on `origin/gotham-1.0` while avoiding Gotham-specific
physics, production inputs, and incidental runtime changes.

This is a porting and design guide, not a record of an already completed branch. The
new implementation will be rebuilt cleanly on the current code baseline rather than
created by merging or cherry-picking the Gotham branch wholesale. The objective is
not merely to recover Gotham functionality: it is to turn that functionality into a
well-designed, maintainable, rigorously tested, documented IO subsystem.

### CP-07 Implementation Refresh

The planning narrative below is retained as decision history. The current feature
implementation has now landed the following contracts, which supersede earlier
forward-looking statements where they differ:

- node restart loading is native and direct: the public `.rst` manifest is parsed
  strictly, each rank routes its local MeshBlock spans from node payloads through
  CP-02 chunked `IOWrapper::Read_bytes_at_all` reads, and production loading never
  stages `<manifest>.assembled`;
- the public node restart entry point is the manifest only; generated
  `node_########/*.g<generation>.payload.rst` paths are rejected as direct `-r`
  arguments, and a payload content marker also rejects hard-link or copied-byte
  aliases that cannot be identified reliably from their pathname;
- manifest loading validates relative generated paths, canonical containment after
  symlink resolution, ordered dense node inventory, exact segment coverage and
  node-local offsets, byte counts, consistent generations, completion, and
  byte-identical replicated payload headers; payload and segment inventories are
  bounded, and segment counts must be positive;
- node-sharded `.bin` and full-volume `.cbin` add preheader inventory metadata
  (`distribution`, `node`, `number of nodes`, and `number of meshblocks`) without
  changing legacy shared or rank files; readers accept explicit empty shards and
  reject incomplete or duplicate node inventories, bound cumulative metadata and
  payload reads, reject non-positive variable counts, and preflight aggregate
  reconstruction plus athdf-like conversion allocations incrementally;
- spherical slices publish through `<file>.tmp` followed by rename, declare
  `layout=dense` or `layout=sparse_angles`, record rank/node inventory metadata for
  sharded files, preserve explicit empty shards, validate complete ownership, bound
  whole-file reads plus coordinate allocations and embedded-header offsets, and
  normalize reconstructed metadata while combining sparse siblings incrementally;
- modern PDF headers and payloads publish independently through checked `<file>.tmp`
  files followed by atomic rename; sparse rank/node headers declare sibling-inventory
  metadata, sibling headers must agree on V2 declaration, declared V2 headers
  require V2 payload preambles with matching cycles, and shipped readers reject
  malformed shard aliases and unreasonable metadata, payload, or aggregate
  reconstruction sizes while retaining only inventory summaries and normalizing
  reconstructed metadata; readers preserve historical transitional unversioned
  dense/sparse payload compatibility;
- frozen `origin/main` shared and per-rank restart fixtures are resumed by regression
  tests in addition to checksum verification;
- `vis/python/bin_convert.py` is the only supported binary conversion API; and
- sliced `.cbin` remains deliberately excluded in every shard mode and is rejected
  during construction.

The requested feature set includes:

1. `single_file_per_node` output and restart sharding.
2. N-dimensional PDF output.
3. Spherical-slice output.
4. Generic diagnostic output variables required by those formats.
5. Output timing diagnostics for measuring IO scaling.
6. A controlled final-output policy suitable for large production jobs.
7. One authoritative binary conversion/reading module, with `bin_convert_new.py`
   folded into `bin_convert.py`.
8. Runnable examples and complete reader/helper tooling for all new formats.
9. A documentation integration bundle designed for later application to
   `origin/gh-pages`, after the code feature branch is complete and merged.

## Non-Negotiable Branch Decisions

### Branch Base

The new feature branch should be created from `origin/main`, currently represented in
the investigation worktree by commit `886dd2a1`.

`origin/gotham-1.0` is not based on current `origin/main`; it has a long, heavily
divergent history. A branch-level merge would import unrelated changes across problem
generators, particles, source terms, mesh behavior, documentation, build settings,
and other modules.

### Source Of Truth For Feature Behavior

The original commits on `origin/gotham-1.0` are the historical evidence for the
features. They will be inspected and selectively reimplemented against the clean base.

The branch `origin/feature/single-file-per-node-outputs` must not be used as an
implementation source. It represents a later attempted extraction and is considered
unreliable for this work. It may be ignored entirely during implementation.

### No Wholesale Cherry-Picks

Even Gotham commits with accurate feature names contain coupled application-specific
edits or intermediate mistakes. The new branch will port behavior and tests, not
commit objects.

### Improve The Design, Do Not Merely Reproduce It

The Gotham implementation is an evidence source, not a specification that overrides
engineering judgment. For every feature family, the implementing agent must first
answer:

1. What user-visible capability is required?
2. What file-format or API compatibility is valuable to preserve?
3. What Gotham design is brittle, incomplete, redundant, untested, or coupled to
   application-specific code?
4. What cleaner implementation can provide the same capability on current
   `origin/main`?
5. What tests, examples, and documentation prove that improvement?

The agent is expected to improve implementation quality where needed. Improvements
may include cleaner abstractions, stronger validation, clearer input names, better
error messages, consolidated Python APIs, safer restart semantics, better timing
measurement, and more robust file metadata. Any intentional change to a Gotham
on-disk layout or user-facing API must be documented, justified, and tested against
the compatibility decision made for that format.

### Production Behavior Must Be Explicit

Two Gotham runtime behaviors are valuable but need cleaner integration:

- Output timing diagnostics should be available without unconditionally changing
  ordinary-run stdout behavior or timing synchronization.
- Checkpoint-only finalization should be configurable, not silently replace the
  established final-output behavior for every user.

### Work Slowly And Require Independent Audit

This is a broad, high-risk IO feature touching file formats, MPI behavior, restart
correctness, Python analysis interfaces, and documentation. The implementing agent
must work deliberately and should not collapse discovery, implementation, and
approval into a single hurried pass.

The required execution discipline is:

1. Start with an explicit implementation plan and update it as each feature family
   is understood, implemented, and verified.
2. Inspect current `origin/main` interfaces and the relevant original Gotham diffs
   before editing each subsystem.
3. Implement one coherent capability at a time and run its focused checks before
   moving on.
4. Spawn subagents for independent, well-scoped audits throughout the work. At a
   minimum, independent audits are required for:
   - Gotham-history scope separation and exclusions;
   - C++ output/file-format and MPI/restart behavior;
   - Python reader/converter API consolidation and backward compatibility;
   - regression tests, negative tests, and runnable examples;
   - the deferred `gh-pages` documentation integration package.
5. Use subagents for bounded code or test work only with disjoint file ownership;
   tell them they are not alone in the codebase and must not overwrite other edits.
6. After each delegated result, the main implementing agent must review the evidence,
   inspect the resulting diff when applicable, and independently run the relevant
   verification. A subagent opinion is not by itself a pass condition.
7. Before merge readiness is claimed, run a final independent audit that looks for
   untested behavior, accidental Gotham-specific imports, format incompatibility,
   stale documentation, and missing helper tooling.

Audit results should be recorded in the branch notes or pull request description with
the files reviewed, commands/tests run, findings addressed, and any unresolved
qualification work.

### Independent Audit Ledger And Stop Gates

The branch must maintain an audit ledger rather than relying on informal statements
that review occurred. Each independent audit entry must record:

- the commit SHA or exact working-tree snapshot reviewed and the feature scope;
- files, fixtures, documentation pages, and historical commits inspected;
- commands, test outputs, comparison scripts, or builds used as evidence;
- findings classified as blocking, non-blocking follow-up, or rejected with reason;
- the resolution commit or change set for each blocking finding; and
- the result of a re-audit after blocking findings are resolved.

The following signoffs are separate stop gates. No one signoff substitutes for
another:

| Audit Gate | Required Reviewer Focus | Merge Blocker |
| --- | --- | --- |
| File-format compatibility | Headers, schema, versioning, filenames, fixtures, old-reader/new-reader boundaries | Unexplained layout or compatibility change |
| Python API consolidation | `bin_convert.py` surface, migrated consumers, readers, CLI, error handling | Residual supported `bin_convert_new` dependency or lost workflow |
| MPI and restart correctness | Node ownership, empty shards, collectives, manifests, resumed numbering | Data loss, deadlock risk, or checkpoint overwrite risk |
| Tests and examples | Harness placement, oracle quality, negative cases, executable workflows | A promoted behavior lacks automated evidence |
| Deferred Pages integration | Final API names, parameters, navigation, Sphinx build instructions | Publishable docs contradict implemented behavior |

An implementing agent must pause at each relevant stop gate, resolve blocking
findings, and request re-audit before proceeding to final merge preparation.

## Proposed Branch Name And Commit Shape

A suitable branch name is:

```text
feature/io-output-formats-and-sharding
```

The anticipated commit sequence is:

1. Add generic derived output variables needed by diagnostic IO.
2. Add N-dimensional PDF representation, input parsing, and reader support.
3. Add spherical-slice output and reader support.
4. Add output-shard mode infrastructure and `single_file_per_node` support.
5. Add per-node restart manifest and payload support.
6. Add output timing diagnostics.
7. Add explicit final-output policy with restart-safe numbering.
8. Consolidate binary conversion into the sole supported `bin_convert.py` API.
9. Add reader/helper scripts, runnable examples, and rigorous regression tests.
10. Add the deferred `gh-pages` documentation integration bundle and validation
    instructions.

The implementation may combine or split commits where compilation dependencies make
that necessary, but each committed change should remain reviewable by feature family.

## Scope Summary

| Capability | Include | Reason |
| --- | --- | --- |
| N-dimensional PDFs | Yes | Core requested IO capability |
| PDF linear/log/symlog axes | Yes | Reusable capability already present in final Gotham implementation |
| PDF volume/mass/variable weighting | Yes | Required for useful multi-field PDFs |
| Sparse COO PDF shard encoding | Yes | Necessary for efficient rank/node-sharded PDFs |
| Spherical-slice output | Yes | Core requested IO capability |
| Spherical-slice Python reader | Yes | Output is not usable without post-processing support |
| `single_file_per_node` mesh outputs | Yes | Core requested IO scaling capability |
| `single_file_per_node` restarts | Yes | Required for large production checkpointing |
| Generic coordinate/flux/passive-scalar output variables | Yes | They are general IO fields and PDF axes |
| Output timing diagnostics | Yes, cleaned design | Directly measures sharding benefit |
| Configurable final-output policy | Yes, cleaned design | Avoids costly unrequested terminal dumps at scale |
| Consolidated `bin_convert.py` reader/converter | Yes | Eliminate duplicated tools and make the newer API canonical |
| Complete new-format reader/helper scripts | Yes | Outputs must be usable in analysis workflows |
| Runnable example inputs and analysis examples | Yes | Demonstrate features and provide regression fixtures |
| Deferred `gh-pages` documentation integration bundle | Yes | Candidate documentation without premature docs-branch changes |
| Independent subagent audits and recorded evidence | Yes | IO/restart/file-format scope warrants separate review passes |
| Automatic terminal diagnostic dump removal by default | No | Would be a backwards-incompatible policy change |
| Cooling-flow `cooling_time` derived variable | No | Coupled to CGM cooling source terms |
| Gotham problem generators and input decks | No | Simulation-specific rather than IO infrastructure |
| Dust and supernova changes | No | Physics-specific |
| New particle restart implementation and MPI hang fixes | No by default | Separate particle-IO feature unless explicitly added |

## Historical Source Map

The following Gotham commits contain the behavior to be reconstructed. Inclusion in
this table does not mean the commit should be cherry-picked. Gotham is historical
evidence and a comparison target, not the specification.

Before implementing each feature family, record the observed Gotham behavior, any
defect or design debt found, the behavior retained, the behavior redesigned, its
compatibility consequence, and tests proving the chosen outcome. Explicit review
targets include file metadata/versioning, parser validation, shard discovery and
empty-shard behavior, MPI collective safety, restart crash/recovery and numbering,
Python API quality, and diagnostic error messages.

| Commit | Subject | Proposed Disposition |
| --- | --- | --- |
| `6774d871` | Add new derived output variables for hydro/MHD moments | Port only generic output-variable portions needed or desirable for IO |
| `82d924d0` | Add N-D PDF output and coordinate/flux derived variables | Port the core N-D PDF and generic diagnostic changes |
| `dcd7c2d1` | Rename H_hydro to enthalpy_plus_ke for clarity | Adopt corrected generic terminology/calculation where flux fields are included |
| `abefea9f` | Fix derived variable output bugs causing GPU memory faults | Port required correctness fixes |
| `d0625a33` | Extend PDF outputs, derived vars, and CC85 test | Port generic PDF/derived-output work; exclude CC85 application content |
| `5efcd425` | finish merging in ND-pdfs | Port only necessary parser/interface fixes |
| `1acd1e8e` | add single_file_per_rank for pdfs | Integrate PDF sharding compatibility into clean design |
| `965cf7ce` | add pdf support for individual passive scalars and update gotham inputs | Port scalar fields; exclude Gotham inputs/pgen edits |
| `d6887323` | add spherical slices output | Port spherical-slice format implementation |
| `d8080a40` | add a python reader for the spherical slice | Port reader functionality |
| `139c094b` | use COO for sparse sfpr pdfs | Port sparse shard representation and reader behavior |
| `97fd9e85` | add abs(costheta) as a pdf coordinate | Port generic coordinate field |
| `edbcd4f3` | Add single_file_per_node output sharding | Port sharding behavior and independent PDF scale-mode behavior |
| `8fae12b1` | Fix reading binary slices when sfpn is used | Port applicable binary-reader correctness fix |
| `868c588e` | Improve single_file_per_node restart and sphslice handling | Port applicable robustness behavior and reconstruct tests |
| `95f5bf56` | outputs: rename sphslice radius tag | Adopt the final stable filename radius token `r_<value>` |
| `59f4b3fe` | Add timing diagnostic for outputs | Port cleaned timing feature only; exclude turb-driver cleanup |
| `e4d5fc2c` | remove final output except restart | Port as explicit selectable final-output policy |
| `b313c5a1` | final outputs don't increment counters | Do not reproduce its numbering semantics; see finalization design below |
| `4b01b047` | hardware flag for fast atomicAdd | Use only as evidence that Gotham corrected its final-output call; do not port unrelated changes |

## Baseline And Compatibility Requirements

The branch must preserve existing output behavior unless a user opts into new
functionality.

### Existing Runs

An existing input file that does not mention new options must continue to:

- parse successfully;
- produce existing output formats with their current names and binary layout;
- produce terminal outputs according to the historical/default policy;
- read historical restart files;
- read historical PDF data using the updated reader.

### New Runs

New functionality is activated through new output options or explicit runtime policy:

- an output block selecting `file_type = sphslice`;
- an output block selecting the N-dimensional PDF form;
- an output block setting `single_file_per_node = true`;
- a run-level setting enabling output timing;
- a run-level setting selecting checkpoint-only finalization.

### Binary-Format Compatibility Contract And Golden Fixtures

The implementation must not claim broad backward compatibility without an explicit
release contract. Before changing any writer or reader, create a matrix for each
supported representation:

| Output Family | Layouts Requiring A Decision |
| --- | --- |
| Binary (`.bin`) | Existing shared/per-rank; new shared/per-node; sliced variants |
| Coarsened binary (`.cbin`) | Existing shared/partitioned; new per-node if supported |
| PDF (`.pdf`) | Existing supported layout; new N-D dense/shared and sparse shards |
| Spherical slice (`.sph.bin`) | New shared/per-rank/per-node representations |
| Restart (`.rst` and manifest/payload family) | Existing shared/per-rank; new per-node |

For every applicable matrix cell, the implementation record must state exactly one
compatibility level:

| Compatibility Level | Meaning |
| --- | --- |
| Byte-identical | Files generated for the legacy mode are byte-for-byte unchanged |
| Header/schema-compatible | Payload/header semantics remain compatible although nonsemantic bytes may differ |
| New-reader compatible | The canonical new reader reads both; older readers need not read the new representation |
| Deliberately incompatible | The change requires a format-version bump and explicit migration documentation |

The compatibility work must include frozen baseline fixtures produced from a clean
`origin/main` checkout before implementation changes are applied. Baseline fixtures
must never be silently regenerated from the feature branch. A fixture manifest must
record the producing commit, input deck, command, build configuration, file size, and
checksum.

Tests against the contract must examine magic/version headers, preheader fields,
numeric type sizes, variable ordering, dimensional metadata, embedded parameter
header behavior, filenames, and shard discovery. They must prove that canonical new
readers consume baseline fixtures under the declared contract, and that legacy output
layouts remain readable wherever the contract promises compatibility.

### Sharding Exclusivity

For any output block supporting shard modes:

```ini
single_file_per_rank = true
single_file_per_node = true
```

must be rejected as invalid. An output instance must use exactly one of:

| Mode | Configuration |
| --- | --- |
| Shared | Both sharding flags false or absent |
| Per-rank | `single_file_per_rank = true` |
| Per-node | `single_file_per_node = true` |

## Feature Family 1: Generic Diagnostic Output Variables

### Motivation

N-dimensional PDFs become substantially more useful when coordinates and fluxes can
be binned without writing simulation-specific derived-data machinery. Gotham includes
a set of broadly applicable derived variables. These belong with IO functionality
because they define reusable output fields and PDF coordinate axes, not a specific
physical setup.

### Included Coordinate Variables

The following output variables should be available for PDFs and any other compatible
cell-centered diagnostic output:

| Variable | Meaning |
| --- | --- |
| `coord_x` | Cartesian `x` coordinate |
| `coord_y` | Cartesian `y` coordinate |
| `coord_z` | Cartesian `z` coordinate |
| `coord_r` | Spherical radius |
| `coord_theta` | Spherical polar angle |
| `coord_phi` | Spherical azimuthal angle |
| `coord_cyl_R` | Cylindrical radius |
| `coord_cyl_phi` | Cylindrical azimuthal angle |
| `coord_cyl_z` | Cylindrical vertical coordinate |
| `coord_costheta` | `cos(theta)` |
| `coord_abscostheta` | `abs(cos(theta))` |

### Included Flux And Velocity Variables

The generic spherical and vertical diagnostics should include:

| Variable Family | Purpose |
| --- | --- |
| `mdot_sph`, `mdot_sph_out`, `mdot_sph_in` | Radial mass flux and directional partitions |
| `edot_sph`, `edot_sph_out`, `edot_sph_in` | Radial energy flux and directional partitions |
| `mdot_vert`, `mdot_vert_out`, `mdot_vert_in` | Vertical mass flux relative to the midplane |
| `edot_vert`, `edot_vert_out`, `edot_vert_in` | Vertical energy flux relative to the midplane |
| `vel_sph_r`, `vel_sph_theta`, `vel_sph_phi` | Velocity components in spherical coordinates |
| `vel_cyl_R`, `vel_cyl_phi` | Velocity components in cylindrical coordinates |
| `edot_sph_kin`, `edot_sph_th`, `edot_sph_mag` | Decomposed radial energy flux components |

When energy-flux expressions are ported, the cleaned implementation should retain
Gotham's corrected `enthalpy_plus_ke` terminology rather than the misleading
intermediate name `H_hydro`.

### Included Passive Scalar Variables

The branch should expose individual hydro and MHD passive scalars as output fields,
following the generic Gotham capability, so that PDFs can be built from individual
scalar components rather than only aggregate quantities.

The implementation must:

- validate scalar indices against the enabled scalar count;
- avoid advertising variables that cannot exist for the active physics package;
- retain established naming if it is unambiguous on the base branch;
- add tests using at least one scalar as either a PDF axis or a PDF weight.

### Excluded Derived Fields

`cooling_time`, introduced on Gotham in `5ad60ec4`, is intentionally excluded from
this branch. It calls into CGM-cooling-specific source-term behavior and is not
required for generic IO formats. It should be proposed separately if needed.

### Correctness Requirements

The port must incorporate the derived-variable allocation/indexing fixes motivated by
`abefea9f`, which addressed GPU memory faults. Validation must include a device-capable
build or, where hardware is unavailable, at minimum compilation and targeted CPU/MPI
coverage plus a recorded GPU validation requirement before merge.

## Feature Family 2: N-Dimensional PDF Output

### User-Facing Capability

The new PDF system should support histograms of one through four variables. The
historical one- or two-variable input form must remain readable for compatibility.

A proposed modern three-dimensional PDF configuration is:

```ini
<output2>
file_type       = pdf
dt              = 1.0
variable_1      = hydro_w_d
bin1_min        = 1.0e-3
bin1_max        = 1.0e2
nbin1           = 128
scale1          = log
variable_2      = hydro_w_vx
bin2_min        = -1.0
bin2_max        = 1.0
nbin2           = 128
scale2          = symlog
linthresh2      = 1.0e-4
variable_3      = coord_abscostheta
bin3_min        = 0.0
bin3_max        = 1.0
nbin3           = 32
scale3          = linear
weight          = mass
```

### Axis Definition

Each axis `N` supports:

| Parameter | Required | Meaning |
| --- | --- | --- |
| `variable_N` | Yes for active axes | Output variable binned along this axis |
| `binN_min` | Yes | Physical lower bin edge |
| `binN_max` | Yes | Physical upper bin edge |
| `nbinN` | Yes | Number of bins |
| `scaleN` | Optional | `linear`, `log`, or `symlog` |
| `linthreshN` | Required for `symlog` | Half-width of the linear region |

Validation must reject:

- gaps in active dimensions, such as `variable_1` and `variable_3` without
  `variable_2`;
- more than four dimensions;
- non-positive `nbinN`;
- `binN_min >= binN_max`;
- non-positive bounds for logarithmic axes;
- non-positive or missing `linthreshN` for symlog axes;
- `linthreshN` supplied for a non-symlog axis unless the implementation clearly
  documents that it is ignored.

### Axis Scales

The branch should include the PDF scale handling found in the final Gotham tree:

| Scale | Intended Use |
| --- | --- |
| `linear` | Coordinates or signed/unsigned quantities with uniform resolution |
| `log` | Positive variables spanning orders of magnitude, such as density |
| `symlog` | Signed variables spanning orders of magnitude, such as velocity or flux |

The C++ writer and Python reader must use matching forward and inverse transforms.
Round-trip tests must verify reconstructed bin edges or centers for all three modes.

### Weights

The modern PDF form should support:

| Setting | Meaning |
| --- | --- |
| `weight = volume` | Accumulate cell volume |
| `weight = mass` | Accumulate cell mass |
| `weight = variable` | Accumulate a specified output field |
| `weight_variable = <name>` | Field used when `weight = variable` |

Legacy `mass_weighted` behavior must either:

- remain supported and be translated into `weight = mass`; or
- be rejected only when mixed inconsistently with the modern options.

The compatibility rule must be documented and tested explicitly.

### Shared And Sharded PDF Storage

Shared PDF output should retain a dense representation. Partitioned PDF output
should use sparse coordinate-list records so empty bins do not multiply file size by
the number of shards.

The planned representation is:

| Layout | Stored Data |
| --- | --- |
| Shared | Dense histogram over all bins |
| Per-rank | Sparse COO histogram for each rank shard |
| Per-node | Sparse COO histogram for each node shard |

The PDF metadata/header must record enough information for readers to reconstruct:

- dimensionality;
- axis variable names;
- bin definitions and scale modes;
- weighting mode;
- dense versus sparse representation;
- shared, rank, or node distribution;
- total bin count.

Each modern PDF header and payload file must be written to `<file>.tmp`, checked,
closed, and atomically renamed independently. This prevents a reader from observing
partial individual files. It is not a filesystem transaction covering the header and
payload as one indivisible family.

### Reader Requirements

`vis/python/read_pdf.py` should:

- read legacy dense PDFs;
- read modern shared dense PDFs;
- read sparse per-rank PDFs by locating matching sibling shards;
- read sparse per-node PDFs by locating matching sibling shards;
- sum sparse records into a dense histogram;
- validate compatible metadata across shards;
- validate canonical sibling paths, shard identifiers, and complete declared rank or
  node inventories before reconstruction;
- reject unreasonable file-controlled dense allocations before construction;
- report missing, duplicated, or inconsistent shard data clearly;
- expose physical bin coordinates for linear, log, and symlog axes.

## Feature Family 3: Spherical-Slice Output

### User-Facing Capability

The branch should introduce `file_type = sphslice` for evaluating selected output
variables on an origin-centered spherical surface.

This is a new binary analysis-oriented format. It must coexist with the existing
`file_type = sph` VTK surface output from `origin/main`; the branch must not silently
rename, remove, or reinterpret legacy `sph`, its `radius` parameter, or its existing
filename convention.

A representative configuration is:

```ini
<output3>
file_type            = sphslice
variable             = hydro_w
dt                   = 1.0
slice_r              = 10.0
ntheta               = 128
nphi                 = 256
single_file_per_node = true
```

Gotham's final input parameter is `slice_r`. Commit `95f5bf56` does not rename that
input; it corrects the radius token in the output filename from `r=<value>` to
`r_<value>`. The clean implementation should use `slice_r` consistently and retain
the corrected filename token.

### Data Model

The output represents a two-dimensional spherical sampling grid:

| Axis | Resolution |
| --- | --- |
| Polar/angular coordinate | `ntheta` |
| Azimuthal coordinate | `nphi` |

The writer should interpolate requested cell-centered fields onto the requested
surface using the implementation's established interpolation behavior.

### Shared Versus Partitioned Files

| Mode | Expected Representation |
| --- | --- |
| Shared | Complete dense spherical surface |
| Per-rank | Sparse angular ownership written by rank shards |
| Per-node | Sparse angular ownership written by node shards |

For partitioned output, each shard must contain enough angular-index information to
reassemble a complete surface without assuming every shard owns data.

The implemented writer publishes every spherical-slice file through
`<file>.tmp` followed by rename. Shared files declare `layout=dense`; rank and node
files declare `layout=sparse_angles`, include shard-inventory metadata, and retain a
valid explicit empty shard when the selected owner has no angles.

### Reader Requirements

`vis/python/read_sphslice.py` should:

- read shared spherical slices;
- accept a path to any rank or node shard;
- discover matching sibling shards;
- validate common metadata;
- validate non-overlapping and complete angular coverage where a dense surface is
  requested;
- reassemble a dense array with a documented variable/angle ordering.

### Spherical-Slice Tests

Tests should compare shared, per-rank, and per-node representations of the same
small problem after reader reassembly. The comparison must cover:

- at least one scalar output;
- at least one multi-component output if supported by the format;
- empty-shard behavior where the decomposition permits it;
- corrected radius-parameter naming.

## Feature Family 4: Output Sharding And `single_file_per_node`

### Goal

The output sharding feature exists to avoid two poor large-run extremes:

- all MPI ranks contending on a single output file; and
- every MPI rank creating its own file for every output event.

`single_file_per_node` writes approximately one data file per compute node while
preserving normal rank ownership of simulation data.

### Shard Mode Abstraction

The implementation should introduce one internal shard-mode abstraction, equivalent
in purpose to Gotham's `FileShardMode`, rather than threading unrelated Boolean tests
through every writer.

Conceptually:

```cpp
enum class FileShardMode {
  shared,
  per_rank,
  per_node
};
```

All output types that support sharding should receive or derive this mode in one
consistent way.

### Node Identification And Communication

Per-node output requires:

- construction of a node-local MPI communicator;
- stable dense node identifiers for directory/file names;
- designation of shard writers where necessary;
- collective behavior that includes ranks with zero owned output data;
- cleanup of communicator resources.

Node sharding affects only IO coordination and on-disk layout. It must not change:

- MeshBlock ownership;
- load-balancing behavior;
- physics data distribution;
- AMR decisions.

### Supported Output Types

The new branch should support the following modes:

| Output Type | Shared | Per-Rank | Per-Node | Notes |
| --- | --- | --- | --- | --- |
| Binary mesh output | Yes | Existing/retained | Add | Include slice-reader correctness |
| Coarsened binary output | Yes for full volume | Existing/retained for full volume | Add for full volume | Sliced `cbin` remains excluded and is rejected explicitly |
| PDF | Yes | Retain/complete | Add | Sparse COO for partitioned files |
| Spherical slice | Yes | Add/retain as implemented | Add | Sparse angular ownership |
| Restart | Yes | Existing/retained | Add | Requires manifest/payload design |

Other formats should remain unchanged unless a small interface adjustment is necessary
for final-output policy handling.

### On-Disk Directory Layout

Gotham's final tree provides the following deterministic layout to preserve or
consciously migrate from:

```text
bin/rank_00000000/<basename>.<id>.00000.bin
bin/node_00000000/<basename>.<id>.00000.bin

pdf_<id>_<axis-names>/rank_00000000/<basename>.<id>.00000.pdf
pdf_<id>_<axis-names>/node_00000000/<basename>.<id>.00000.pdf

bin/rank_00000000/<basename>.<id>.r_<radius>.00000.sph.bin
bin/node_00000000/<basename>.<id>.r_<radius>.00000.sph.bin

cbin_<id>_<factor>/rank_00000000/<basename>.<id>.<number>.cbin
cbin_<id>_<factor>/node_00000000/<basename>.<id>.<number>.cbin

rst/rank_00000000/<basename>.00000.rst
rst/<basename>.00000.rst
rst/node_00000000/<basename>.00000.g<generation>.payload.rst
```

Exact prefixes, axis-name encoding, and number formatting must be checked while
porting onto `origin/main`. Any deliberate layout change requires corresponding
reader updates, compatibility tests, and documentation; otherwise the final Gotham
layout above is the behavioral reference.

### Empty-Shard Handling

An output event may have a rank or node with no data for a particular slice or sparse
diagnostic. Writers must not deadlock or repeatedly retry because a shard owns zero
payload.

Required behavior:

- collective operations remain collective across the relevant communicator;
- node-sharded writers publish a valid explicit empty representation when a node
  owns no records for the selected product;
- node-sharded binary and full-volume coarsened-binary shards carry additive
  `distribution`, `node`, `number of nodes`, and `number of meshblocks` preheader
  fields, and spherical-slice shards carry layout plus rank/node inventory metadata;
- output schedule counters advance once for an ordinary scheduled output event;
- the Python reader accepts legitimate empty shards but rejects missing, duplicate,
  malformed, or non-dense node inventories when the additive metadata is present.

## Feature Family 5: Per-Node Restart Files

### Purpose

At large scale, restart output is often the most important and largest IO operation.
Per-node restart support is therefore a first-class part of `single_file_per_node`,
not an optional follow-up.

### Proposed Manifest/Payload Layout

Per-node restart output should write:

```text
rst/<basename>.00000.rst
rst/node_00000000/<basename>.00000.g<generation>.payload.rst
rst/node_00000001/<basename>.00000.g<generation>.payload.rst
...
```

The root file is a manifest containing global information and sufficient layout
metadata to locate and interpret payload shards. Each node file contains the heavy
restart payload for data assigned to that node.

### Restart Entry Points

The supported node restart entry point is the public manifest path only. Direct
restart from a generated node payload is rejected with an instruction to use the
manifest. Ordinary shared restart files remain valid even if an unrelated filename
ends in `.payload.rst`. Each generated node payload also carries an explicit marker
after the replicated parameter dump. Manifest loading requires and consumes the
marker, while ordinary restart loading rejects marked bytes. This closes direct
restart through hard links or copied payload bytes without changing legacy shared
or per-rank restart files.

### Native Direct Loading

The earlier implementation temporarily assembled a rank-0 shared restart file before
calling legacy restart loading. That prior state is retained as useful decision
history, but it is not the production path. The current loader:

1. Parses and validates the public manifest.
2. Opens one canonical payload for the replicated restart header.
3. Checks byte-identical replicated headers across payloads.
4. Routes each rank's local MeshBlock spans to the declared node payloads.
5. Reads those spans directly through CP-02 chunked `IOWrapper` positioned reads.

No production `<manifest>.assembled` or `<manifest>.assembled.tmp` artifact is
created. Tests reserve those names as directories during resume to prove that native
loading does not depend on them.

### Manifest Validation

The manifest loader rejects absolute or traversing paths, malformed node-directory
or generation components, payload symlink escapes after canonicalization, incomplete
or reordered fixed records, non-contiguous node inventories, mixed generations,
inconsistent byte counts, segment gaps or overlaps, inconsistent node-local offsets,
missing or truncated payloads, replicated-header mismatches, oversized declared
payload or segment inventories, non-positive segment counts, and absent or corrupt
payload markers.

### MPI-IO Robustness

The cleaned implementation should incorporate the restart robustness intent from
Gotham while avoiding unrelated particle-output behavior:

- no MPI collective may exclude zero-payload participants in its communicator;
- large byte transfers should be chunkable so they do not depend on an `int`-sized
  MPI count;
- open/close operations and header/payload roles should be unambiguous;
- stale or conflicting node-shard directories must be detected or safely replaced
  according to the final file-overwrite policy.

### Particle Restart Boundary

The codebase may already have `rst_prtcl` output and finalization should classify it
as a restart-style output when it exists. However, the new branch does not
automatically include Gotham's particle-restart rewrite or its MPI append-mode hang
fix. That work should be a separately reviewed extension unless explicitly requested.

## Feature Family 6: Output Timing Diagnostics

### Why It Belongs In This Branch

`single_file_per_node` is motivated by output scalability. It is not sufficient to
provide new layouts without a practical mechanism to measure whether they improve
runtime on target systems.

Gotham commit `59f4b3fe` demonstrated the basic instrumentation:

- fence before output;
- start a `Kokkos::Timer`;
- load and write output;
- fence after output;
- print elapsed time on rank 0.

That behavior is useful, but the clean port should address its limitations.

### Proposed Interface

Timing should be opt-in so ordinary simulations retain their established log volume
and avoid synchronization solely for diagnostics.

A proposed run-level parameter is:

```ini
<time>
output_timing = true
```

If a different existing diagnostics block or naming convention is more appropriate on
the base branch, use that convention while preserving the opt-in semantics.

### Events To Time

When enabled, output timing should cover:

| Event | Reporting |
| --- | --- |
| Initial output set | Aggregate timing for startup outputs |
| Each scheduled output block | Timing identified by output block and file type |
| Final output set | Aggregate timing for outputs actually selected by final policy |

### MPI Reporting

The clean implementation should measure local elapsed time and report the MPI maximum
over all ranks involved in the simulation. Rank 0 should print that maximum.

The maximum is preferable to Gotham's rank-0-only report because overall progress is
limited by the slowest participant in collective or synchronized IO.

A representative log format is:

```text
[output-io] block=output2 type=pdf mode=per_node elapsed_max_s=1.284
```

### Timing Semantics

The guide intentionally accepts synchronization overhead when timing is enabled:

- a pre-output fence excludes outstanding compute work from the output timer;
- a post-output fence ensures device copies and writer-related kernels are included;
- MPI reduction adds small diagnostic overhead;
- timing-disabled runs do not incur these extra measurement fences/reductions solely
  for reporting.

### Timing Exclusions

Do not port unrelated changes from `59f4b3fe`, specifically deletion of debug prints
from `src/srcterms/turb_driver.cpp`. That cleanup is outside this IO scope.

### Timing Validation

Tests should verify:

- the parameter enables timing messages;
- disabled timing does not print the new timing records;
- messages identify shared, per-rank, and per-node modes correctly;
- MPI timing reporting functions with multiple ranks;
- timing does not alter output data or output scheduling.

Performance benchmarking on a real multi-node system is a release qualification task,
not a deterministic regression test.

## Feature Family 7: Final-Output Policy

### Existing Baseline Behavior

On current `origin/main`, `Driver::Finalize()` writes every configured output object
after the execution loop completes. This creates a terminal diagnostic snapshot even
if the diagnostic was not scheduled exactly at the final simulation time.

This behavior is convenient for small runs but potentially expensive at scale,
especially once PDFs and spherical slices are large or node-sharded.

### Gotham Behavior

Gotham changed finalization in stages:

1. `f86b6eb0` disabled all final output as an intermediate application-driven change.
2. `e4d5fc2c` restored final output only for `rst` and `rst_prtcl`.
3. `b313c5a1` added no-counter-advance logic for final output writers, but mistakenly
   called final restart writing with `is_final = false`.
4. `4b01b047` corrected the Gotham finalization call to `is_final = true`.

The final Gotham branch therefore writes checkpoint-style files at shutdown but not
ordinary terminal diagnostic files.

### Cleaned Design: Explicit Policy

The branch should provide Gotham's useful checkpoint-only behavior as a selectable
runtime policy rather than silently changing every existing run.

A proposed run-level setting is:

```ini
<time>
final_output_policy = restart_only
```

Supported policies should be:

| Policy | Finalization Behavior | Recommended Use |
| --- | --- | --- |
| `all` | Write all configured outputs at shutdown | Backwards-compatible default and small tests |
| `restart_only` | Write `rst` and existing `rst_prtcl` outputs only | Large production runs |
| `none` | Write no automatic shutdown outputs | Explicit expert use or external checkpoint control |

Recommended default:

```ini
final_output_policy = all
```

This retains current behavior for existing input files while giving Gotham-style
production jobs an explicit `restart_only` path.

### Why The Gotham Counter Semantics Should Not Be Copied

Gotham's final implementation writes a numbered terminal restart while suppressing
advancement of `file_number` and `last_time`. On resume in the same output directory,
the next scheduled restart can reuse the same output number and overwrite the terminal
checkpoint.

That is not acceptable as a default checkpointing behavior. A final checkpoint is an
actual emitted output and should consume its output sequence number.

### Proposed Counter Rule

Every output file that is actually written, including a finalization output, advances
its scheduling metadata exactly as an ordinary output event does.

For a final restart:

1. Construct the filename using the current `file_number`.
2. Advance `file_number` and relevant schedule metadata before embedding the
   parameter dump in the restart file, matching established restart behavior.
3. Write the restart file.
4. On resume, read the advanced output metadata from the restart file.
5. Any later restart write uses a new filename and cannot overwrite the terminal
   checkpoint by repeating its number.

This keeps the benefit of `restart_only` finalization without inheriting Gotham's
numbering risk.

### API Design

It is not necessary to duplicate Gotham's `is_final` behavior in every output writer
solely to suppress counters. The clean implementation may choose one of these
approaches:

| Approach | Assessment |
| --- | --- |
| Driver selects output types at finalization; writers always advance when writing | Preferred: simplest semantics |
| Writers receive a finalization context but still advance for emitted files | Acceptable if needed for headers/logging |
| Writers receive `is_final` and suppress counters | Reject: recreates numbering risk |

The final implementation should favor policy selection in `Driver::Finalize()` and
keep writer scheduling behavior uniform.

### Finalization Tests

Regression coverage must verify:

| Test | Expected Result |
| --- | --- |
| Default policy with diagnostic output | Existing final diagnostic dump remains present |
| `restart_only` with PDF configured | No unscheduled terminal PDF dump |
| `restart_only` with spherical slice configured | No unscheduled terminal slice dump |
| `restart_only` with binary configured | No unscheduled terminal binary dump |
| `restart_only` with restart configured | Final restart exists |
| `none` with restart configured | No automatic final restart exists |
| Resume from final restart and evolve to another restart | New restart number is used; final checkpoint is not overwritten |
| Per-node final restart then resume | Manifest and payload shards remain valid and numbering advances |

## Reader And Analysis Tooling

### Tooling Is Part Of The Feature, Not A Follow-Up

No new output format is complete merely because C++ can write it. The branch must
ship the Python interfaces needed to read, validate, combine, inspect, and convert
every output layout introduced here. The reader implementations, CLI behavior,
examples, and tests are first-class acceptance criteria.

### One Canonical Binary Converter: `bin_convert.py`

The existing Gotham history exposes a tooling problem: both
`vis/python/bin_convert.py` and `vis/python/bin_convert_new.py` exist, and downstream
users can end up depending on different implementations. This branch must remove
that ambiguity.

The user's supported workflow is based on `bin_convert_new.py`. Therefore the clean
port must:

1. Treat the design and functionality of `bin_convert_new.py` as the primary source
   for the consolidated implementation.
2. Produce one canonical supported module named `vis/python/bin_convert.py`.
3. Fold into that canonical module any still-needed behavior currently provided only
   by the older `bin_convert.py`, including HDF5/XDMF writing helpers or public API
   functions on which scripts/tests depend.
4. Add all new binary/per-rank/per-node/sliced-output reading behavior to the
   canonical `bin_convert.py`.
5. Remove `vis/python/bin_convert_new.py` from the final intended API and update all
   in-repository imports, documentation, and examples to import `bin_convert`.
6. Avoid maintaining a compatibility shim named `bin_convert_new.py` unless a
   concrete external compatibility requirement is identified, documented, and
   approved. A permanent two-module alias would reintroduce the redundancy this
   branch is intended to remove.

The implementing agent must inventory both modules before editing. At minimum, the
inventory must identify:

- public callable functions present in each module;
- code paths importing either module, such as visualization or conversion scripts;
- whether the newer implementation calls back into the older one for helpers such as
  `write_athdf`;
- CLI usage or scripts outside the module that would break when the duplicate file is
  removed;
- behavior required for old shared binary and coarsened-binary files.

The consolidated module should have a clearly declared public API, shared internal
parsing/reassembly primitives, type annotations where practical, and no duplicated
alternate implementation hidden elsewhere in `vis/python`.

### Canonical Converter API Requirements

The exact final API should be confirmed against current code and test consumers, but
the unified `bin_convert.py` must provide the equivalent of the useful modern
surface:

| API Family | Required Role |
| --- | --- |
| Raw binary readers | Read shared, per-rank, and per-node binary data |
| Coarsened binary readers | Read supported shared and partitioned coarsened data |
| `*_as_athdf` views | Assemble data into the existing in-memory analysis convention |
| HDF5/XDMF output helpers | Preserve necessary conversion workflows from the older module |
| File conversion CLI/function | Convert supported binary output layouts through one entry point |
| Input/header helpers | Preserve only documented, used helpers or provide a migration path |

Per-node discovery and reassembly must be integrated into the same public entry
points wherever possible. Users should not need a separate converter just because an
output was generated with `single_file_per_node = true`.

### Canonical Converter Migration Ledger

Before deleting or renaming converter code, create a migration ledger containing each
public function, signature/return-shape decision, known repository consumer,
preserve/migrate/remove decision, and its required test. The first pass must include
the known existing workflow in `vis/python/make_athdf.py`, which depends on
`read_binary`, `write_athdf`, and `write_xdmf_for`. Because the modern Gotham
`bin_convert_new.py` implementation does not itself provide every legacy writer
helper, simply replacing the older file is insufficient.

The minimum ledger is:

| Public Capability | Required Decision | Required Evidence |
| --- | --- | --- |
| `read_binary` and related raw readers | Retain modern behavior in `bin_convert.py`; extend for node shards | Import test and shared/per-rank/per-node fixture readback |
| `read_*_as_athdf` assembly interfaces | Retain or document changed signatures/shapes | Key/shape/value assertions on old and new layouts |
| `write_athdf` | Preserve in canonical module or deliberately migrate `make_athdf.py` | HDF5 generation test through supported workflow |
| `write_xdmf_for` | Preserve if the ATHDF workflow remains supported | XDMF output smoke/contents test |
| CLI conversion entry point | Preserve or announce replacement explicitly | CLI conversion tests for supported layouts |
| Module import name | Public users and examples import `bin_convert` only | Repository import sweep and documentation check |

Any API removal must identify the in-repository consumers inspected, the replacement
path, and the tests updated. The merged documentation and all promoted examples must
expose only the canonical `bin_convert.py` interface.

### Included Python Files

The feature branch should contain user-facing tooling alongside writers:

| Tool | Required Capability |
| --- | --- |
| `vis/python/read_pdf.py` | Dense and sparse N-D PDF reading; rank/node shard reconstruction; scale metadata |
| `vis/python/read_sphslice.py` | Shared and partitioned spherical-slice reconstruction |
| `vis/python/bin_convert.py` | Sole supported binary/coarsened-binary reader and converter, based primarily on the newer implementation and extended for new IO |

The writer-to-tool contract is:

| Emitted Output | Required User-Facing Read Path |
| --- | --- |
| Binary, sliced binary, and coarsened binary | Unified `vis/python/bin_convert.py` |
| N-dimensional PDF | `vis/python/read_pdf.py` |
| Spherical slice | `vis/python/read_sphslice.py` |
| Per-node restart manifest and payloads | Strict public-manifest validation followed by native routed payload reads; generated payload paths are not public restart entry points and no production `.assembled` staging is used |

No additional reader helper is required merely to split logic by sharding mode. If a
shared utility module is warranted to remove meaningful duplication between readers,
it must have a narrow documented role and dedicated tests.

### Tooling Compatibility

Python tools should be able to read:

- output generated by the clean feature branch;
- historical compatible shared output from the base branch;
- Gotham-style files where the clean on-disk representation intentionally matches it.

The new implementation is not required to support malformed files or every
intermediate Gotham layout. Any intentionally unsupported intermediate format should
be documented.

### Reader Validation And Error Handling

Every reader must validate sufficient metadata before silently combining shards.
Required failure cases include:

- sibling shards with incompatible output time, variable set, dimensions, scale
  metadata, or file-format version;
- duplicate or missing spherical-slice angle ownership;
- duplicate or out-of-range sparse PDF entries where not permitted;
- a missing node/rank shard when completeness is required;
- a requested file that belongs to a sharded family but cannot be resolved safely;
- binary or coarsened-binary MeshBlocks with nonuniform emitted extents inside
  one file;
- malformed headers or unsupported layout versions.

Error messages should identify the file family and the incompatible field rather than
surface an opaque NumPy reshape or file-open error.

### Converter And Reader Migration Tests

Testing for the tool consolidation must include:

1. An API inventory test or explicit review checklist demonstrating that required
   `bin_convert_new.py` functionality exists in the final `bin_convert.py`.
2. An import sweep proving repository code and examples no longer import
   `bin_convert_new`.
3. Conversion of existing shared binary and coarsened-binary fixtures through the
   unified module.
4. Conversion or in-memory read of per-rank and per-node binary outputs.
5. Binary slice reconstruction, including an empty node shard when applicable.
6. Preservation of required HDF5/XDMF workflows formerly supplied by the older
   module.
7. Command-line smoke tests if the converter retains a CLI entry point.
8. Reader tests for `read_pdf.py` and `read_sphslice.py` covering both valid
   reassembly and corrupted/inconsistent shard metadata.

The branch must not delete `bin_convert_new.py` until equivalent supported behavior is
present in `bin_convert.py` and all in-repository dependencies are updated and
verified.

## Runnable Examples And Reference Workflows

### Purpose

Examples are not decorative documentation. They must serve as small, executable
workflows demonstrating how a user produces and reads each new format, and they
should be reused by regression testing wherever feasible.

### Required Example Inputs

The feature branch should add compact example or test input decks covering:

| Example | Must Demonstrate |
| --- | --- |
| N-D PDF example | At least three axes, mixed scale modes, a useful weight, and Python readback |
| Spherical-slice example | `slice_r`, resolution controls, shared output, and Python readback |
| Per-node binary example | `single_file_per_node`, conversion through canonical `bin_convert.py` |
| Per-node PDF/spherical-slice example | Sharded files and reassembly into analysis arrays |
| Per-node restart example | Manifest/payload output and restart invocation |
| Production-output policy example | `output_timing` plus `final_output_policy = restart_only` |

Where one carefully designed small input can exercise several capabilities without
becoming obscure, it may be reused. Otherwise, keep examples focused and readable.

### Required Python Example Workflows

The branch should include executable snippets or scripts that demonstrate:

- loading a modern N-D PDF and inspecting bin axes/weights;
- loading a shared and sharded spherical slice;
- reading/converting a shared and per-node binary file with `bin_convert.py`;
- locating and using a per-node restart manifest;
- interpreting output timing records for a scaling comparison.

Examples must import the final public module names. They must never teach users to
import `bin_convert_new`.

### Example Verification

At least the small examples intended for users must be exercised in automated tests
or a repeatable smoke-test command. The test should verify both successful output
generation and successful analysis readback, not just process exit status.

## Anticipated Code Impact

The exact diff will depend on the current base implementation, but the intended write
scope is approximately:

| Area | Expected Files | Purpose |
| --- | --- | --- |
| Output registration and parameters | `src/outputs/outputs.cpp`, `src/outputs/outputs.hpp` | New formats, PDF options, sharding/policy context |
| Derived variables | `src/outputs/basetype_output.cpp`, `src/outputs/derived_variables.cpp` | Generic diagnostic fields |
| PDF writer | `src/outputs/pdf.cpp` | N-D accumulation, scale handling, sparse sharding |
| Spherical-slice writer | `src/outputs/spherical_slice.cpp` plus build registration | New output format |
| Binary/coarsened writers | `src/outputs/binary.cpp`, `src/outputs/coarsened_binary.cpp` | Node shard writes and empty-shard behavior |
| Restart writer/reader plumbing | `src/outputs/restart.cpp`, `src/main.cpp`, `src/mesh/*` as required | Per-node manifest and payload loading |
| IO wrapper/sharding utility | `src/outputs/io_wrapper.*`, new sharding helper if warranted, globals only if unavoidable | Node communicators and robust IO operations |
| Driver | `src/driver/driver.cpp` | Timing and final-output policy |
| Input parsing | `src/parameter_input.*` only if necessary | New option parsing/support |
| Python readers/converters | `vis/python/read_pdf.py`, `vis/python/read_sphslice.py`, `vis/python/bin_convert.py`; remove `vis/python/bin_convert_new.py` after migration | Canonical reassembly and analysis support |
| Importing Python utilities | Current scripts importing either converter | Migrate to the sole `bin_convert` API |
| Regression tests | `tst/inputs/*`, new `tst/test_suite/io/*` tests, and a chosen immutable-fixture directory | Focused, negative, and compatibility coverage using the actual harness |
| User examples | `inputs/io/*` decks and `vis/python/examples/*` readback scripts | Executable user workflows following repository input conventions |
| Deferred docs integration package | Feature-branch staging artifacts plus an integration manifest | Ready later application to `gh-pages` |

Changes outside this scope require a specific justification in the implementation
commit or review description.

## Excluded Code And Features

The following must not enter the IO feature branch merely because they are adjacent
in Gotham history:

| Excluded Item | Reason |
| --- | --- |
| `inputs/gotham/*` and Gotham job scripts/tables | Production setup, not reusable IO infrastructure |
| `src/pgen/cgm_cooling_flow_full.cpp` and related cooling-flow problem generators | Application-specific |
| Dust-model work | Physics feature unrelated to IO |
| Supernova momentum changes | Physics behavior unrelated to IO |
| CGM `cooling_time` derived output and its source-term edits | Requires separate physics-backed scope |
| Turbulence-driver debug-print cleanup in `59f4b3fe` | Incidental cleanup |
| Global hardware/atomic changes from `4b01b047` | Unrelated to final-output correction |
| Particle restart append-mode rewrite from `3f77e59a` | Separate particle IO concern unless later approved |
| Any implementation taken from `origin/feature/single-file-per-node-outputs` | Explicitly rejected source branch |
| A continuing `bin_convert_new.py` / `bin_convert.py` dual public API | Explicitly rejected tooling redundancy |

## Regression Test Plan

### Fixture, Example, And Harness Layout

The current test harness discovers Python tests below `tst/test_suite/`, with target
suffixes such as `_cpu`, `_mpicpu`, and `_gpu`, and uses input decks below
`tst/inputs/`. The feature implementation must follow those conventions rather than
inventing an unconnected test layout.

Required placement and separation:

| Artifact | Required Location Or Rule |
| --- | --- |
| Runtime IO tests | Add a coherent `tst/test_suite/io/` family, using `_cpu`, `_mpicpu`, or `_gpu` only as the needed build/run target dictates |
| Producer input decks | Small deterministic decks under `tst/inputs/` |
| Golden baseline files | A clearly named committed fixture directory chosen consistently with test conventions, plus a manifest of source commit, generation command, sizes, and checksums |
| Generated example outputs | Test scratch/output directories, never confused with immutable compatibility fixtures |
| Reader-only malformed fixtures | Small corrupt, truncated, or metadata-inconsistent artifacts runnable without rebuilding or launching AthenaK |

All IO-producing tests must use isolated run directories or equivalent comprehensive
cleanup. The existing general-purpose cleanup does not remove every new IO directory
or artifact and therefore cannot by itself prevent cross-test contamination.

Tests should be added beside each coherent feature implementation, not postponed
until all formats and policies have been ported.

### Build Matrix

At minimum:

| Build | Purpose |
| --- | --- |
| Serial CPU build | Fast format and compatibility tests |
| MPI CPU build | Per-rank/per-node semantics and restart tests |
| GPU-capable build or documented CI target | Derived-variable and device-memory safety |

### Test Design Standards

The test suite must be designed as part of implementation, not appended after the
formats appear to work. For each output family, require:

- a small deterministic producer input;
- a Python readback/validation path using the shipped public reader;
- cross-mode comparisons where shared, per-rank, and per-node representations should
  represent identical data;
- invalid-input and corrupted-file negative tests where validation logic is added;
- restart or repeat-write tests where scheduling state or overwrite behavior matters;
- assertions on metadata and file discovery, not just numerical payload values;
- clear tolerances and precision rationale for interpolated or reduced data.

Where an existing test framework cannot support a needed scenario, add the smallest
reusable harness improvement required rather than substituting manual verification.

### PDF Tests

Required automated coverage includes:

1. Legacy one-dimensional PDF compatibility.
2. Legacy two-dimensional PDF compatibility if supported by the base.
3. Three-dimensional PDF with mixed `linear`, `log`, and `symlog` axes.
4. Four-dimensional parsing and allocation on a small dataset.
5. `weight = volume`.
6. `weight = mass`.
7. `weight = variable`.
8. Generic coordinate axes, including `coord_abscostheta`.
9. A passive scalar axis or weight.
10. Shared versus per-rank reconstructed equality.
11. Shared versus per-node reconstructed equality.
12. Invalid-input rejection cases for axis scales and dimensions.

### Spherical-Slice Tests

Required automated coverage includes:

1. Shared spherical-slice creation and reader load.
2. Per-rank reassembly.
3. Per-node reassembly.
4. Equality of reconstructed values between modes within tolerance.
5. Radius parameter parsing through `slice_r` and filename emission using
   `r_<value>`.
6. Invalid or incomplete shard detection.

### Binary And Coarsened-Binary Tests

Required automated coverage includes:

1. Shared binary output.
2. Per-rank binary reconstruction if already available.
3. Per-node binary reconstruction.
4. Sliced binary output with a node that owns no selected slice data.
5. Equivalent tests for coarsened binary where supported.

### Canonical `bin_convert.py` Tests

Required automated coverage includes:

1. Legacy shared binary file reading and conversion.
2. Legacy shared coarsened-binary reading and conversion.
3. Modern API behavior migrated from `bin_convert_new.py`.
4. Per-rank shard discovery and reconstruction.
5. Per-node shard discovery and reconstruction.
6. Sliced binary data in all supported layouts.
7. HDF5/XDMF output generation if retained as supported API.
8. Repository import scanning to ensure no supported script/example imports
   `bin_convert_new`.
9. Absence of duplicated converter logic or a documented temporary migration
   exception with an explicit removal gate.

### Restart Tests

Required automated coverage includes:

1. Shared restart round trip.
2. Existing per-rank restart round trip.
3. Per-node restart manifest and payload creation.
4. Restart by manifest path.
5. Restart by normalized node-shard path if supported.
6. A decomposition containing zero-payload ranks/nodes where feasible.
7. Forced small transfer chunks to exercise chunked IO logic.
8. Terminal checkpoint numbering followed by resumed evolution.
9. `restart_only` final policy with per-node restarts.

### Timing Tests

Required automated coverage includes:

1. No new timing record when timing is disabled.
2. Timing records for scheduled shared output when enabled.
3. Timing records identifying `per_rank` and `per_node` modes.
4. Final timing report includes only output selected under the final policy.

### Final-Output Policy Tests

Required automated coverage spans all policies:

| Policy | Required Coverage |
| --- | --- |
| `all` | Retains historical terminal diagnostic output |
| `restart_only` | Writes terminal checkpoint, suppresses terminal diagnostics |
| `none` | Writes no automatic final output |

### Example And Documentation Tests

Required automated tests or repeatable validation commands include:

1. Every promoted example input produces the expected output family.
2. Each associated Python snippet/script reads its generated output successfully
   using only public shipped tooling.
3. Examples for node-sharded output verify reassembly rather than merely listing
   output files.
4. Parameter names and defaults quoted in documentation match the implemented parser.
5. The deferred `gh-pages` content builds without warnings in a temporary
   `origin/gh-pages` worktree after it is staged for integration.

## Performance Qualification Plan

Functional tests cannot establish the value of output sharding. Before merge or
production adoption, run a controlled scaling comparison:

| Comparison | Measure |
| --- | --- |
| Shared vs per-rank vs per-node binary output | Output elapsed maximum, file count, total bytes |
| Shared vs per-node PDF output | Output elapsed maximum and reader reconstruction cost |
| Shared vs per-node spherical slice | Output elapsed maximum and file-count reduction |
| Shared vs per-node restart | Checkpoint elapsed maximum and restart read elapsed time |

The benchmark should record:

- machine/filesystem;
- nodes, MPI ranks, and ranks per node;
- MeshBlock count and representative payload size;
- enabled output blocks;
- output timing records;
- file counts and total output size;
- any failures, stalls, or unexpected filesystem behavior.

## Deferred `gh-pages` Documentation Integration

### Publication Timing Rule

The feature branch must include extensive, merge-ready documentation content, but it
must not modify, merge into, or publish the live `gh-pages` branch before the code
feature branch has been accepted and merged. CP-03 native direct restart loading has
replaced the earlier transient rank-0 `<manifest>.assembled` staging design. Do not
publish or open a Pages integration review until the staged pages are reapplied,
rebuilt, and re-audited against the then-current `origin/gh-pages`. Documentation
may be validated against a temporary worktree of `origin/gh-pages`; validation is
not publication.

### Live Documentation Framework To Target

The authoritative documentation framework is the separate `origin/gh-pages` branch.
At the time this guide was drafted, its relevant routing is:

| Destination | Role |
| --- | --- |
| `docs/source/modules/outputs.md` | Existing Outputs module page to extend or supersede carefully |
| `docs/source/modules/index.md` | Modules catalogue and hidden module toctree |
| `docs/source/index.md` | Global navigation and tools/examples routing |
| `docs/source/tools/visualization.md` | Current binary-converter/visualization user guidance |
| `docs/source/reference/input_parameters.md` | Input-parameter reference to update with finalized options |
| `docs/source/reference/file_reference.md` | File/read-path reference where new output artifacts warrant additions |
| `docs/source/examples/index.md` | Entry point for runnable output examples |

The live visualization page currently describes `bin_convert_new.py` as the preferred
modern interface and `bin_convert.py` as a legacy dependency. The eventual
documentation integration must deliberately replace that guidance with the single
canonical `bin_convert.py` interface.

### Known Live Documentation Reconciliation Tasks

The eventual documentation integration is an update to existing guidance, not an
assumption that the live Pages branch is empty. At minimum it must reconcile:

| Live Page | Required Reconciliation |
| --- | --- |
| `docs/source/modules/outputs.md` | Replace or update existing N-D PDF parameter guidance, including existing `logscaleN` material versus the finalized scale/symlog design; add sharding, spherical slices, restarts, timing, and final-output policy where appropriate |
| `docs/source/tools/visualization.md` | Replace the two-converter narrative with canonical `bin_convert.py` and document all supported readers |
| `docs/source/configuration.md` and `docs/source/running.md` | Reconcile user-facing sharding, timing, and final-output policy configuration guidance where applicable |
| `docs/source/reference/input_parameters.md` | Verify finalized option names, defaults, accepted values, and validation failures against parser code |
| `docs/source/reference/file_reference.md` | Add new file families, naming/layout conventions, and readers if warranted by its final role |
| `docs/source/examples/index.md` | Register the executable IO example workflows and route users to their readback steps |
| `docs/source/modules/index.md` and `docs/source/index.md` | Verify existing routing; change them only if a new page target requires registration |

### Feature-Branch Documentation Package

Once names, defaults, layouts, and readers are implemented and tested, the feature
branch should add a deferred documentation package suitable for copying or applying
to `origin/gh-pages` after merge. If `origin/main` still does not carry the live
Sphinx/MyST documentation tree, do not import that whole tree merely for this feature.
Instead, keep a clearly named staged integration bundle and manifest that records the
target `gh-pages` paths and navigation edits.

The package must include publishable MyST content for:

| Proposed Page Or Update | Required Content |
| --- | --- |
| Outputs module update or focused IO extension page | Formats, sharding modes, data layout, compatibility, final-output policy |
| Visualization/tools update | Canonical `bin_convert.py`, `read_pdf.py`, `read_sphslice.py`, Python examples |
| Input-parameter reference update | All final parser keys, defaults, validation rules, and allowed values |
| File-reference update where appropriate | Final output extensions, shard naming/layout, manifest/payload relationship, reader paths |
| Worked-example page/update | Runnable input decks and readback commands/scripts |
| Integration manifest | Exact target files, navigation changes, stale-page cleanup, build command and expected checks |

The staging location may be selected during implementation to match repo conventions,
but it must be unmistakably deferred content rather than an attempted live
`gh-pages` update.

### Required User Documentation Content

The eventual Pages content must document:

1. N-dimensional `pdf` axes, scales, weights, dense/shared layout, and sparse
   sharded layout.
2. Spherical-slice configuration, `slice_r`, output filenames, data ordering, and
   reader usage.
3. Shared, per-rank, and per-node file layouts for supported output types.
4. Per-node restart manifests, payload shards, restart invocation, and numbering
   semantics under finalization.
5. `output_timing` behavior, overhead, and interpretation of reported MPI maximums.
6. `final_output_policy` values, default behavior, and production guidance.
7. The single supported `bin_convert.py` API, migration away from
   `bin_convert_new.py`, and all new-format readers.
8. Complete runnable examples and expected output/readback checks.
9. Known compatibility boundaries, including any intentionally unsupported
   intermediate Gotham layouts.

### Merge-Time Documentation Workflow

After the code branch is merged, the integrating agent must:

1. Re-read the final implemented parameters, file layouts, public Python APIs, and
   test examples from code rather than trusting this planning guide.
2. Create or update the candidate MyST pages in a temporary worktree based on
   `origin/gh-pages`.
3. Apply page additions/updates to the appropriate locations under
   `docs/source/modules/`, `docs/source/tools/`, `docs/source/reference/`, and
   `docs/source/examples/`.
4. Add new navigation targets to both `docs/source/index.md` and the appropriate
   catalogue/toctree in `docs/source/modules/index.md` or `docs/source/examples/index.md`.
5. Search existing live pages for contradictions, especially stale references to
   `bin_convert_new.py`, outdated output defaults, or absent sharding modes.
6. Run:

   ```bash
   cd docs
   make clean html SPHINXOPTS="-W --keep-going"
   ```

7. Record the target files changed, navigation edits, stale guidance removed, and
   successful warnings-as-errors build result.
8. Only then submit the separate `gh-pages` integration change.

### Documentation Independent Audit

Before the separate `gh-pages` change is submitted, a subagent must independently
compare:

- final code parser/defaults;
- final file-format and reader behavior;
- test/example inputs and readback scripts; and
- proposed Pages content and navigation wiring.

The main integrating agent must resolve or explicitly record every discrepancy before
claiming the documentation is publishable.

## Review Gates

The feature branch is ready for review only when all of the following are true:

- The implementation is based on clean `origin/main`.
- No changes are taken from `origin/feature/single-file-per-node-outputs`.
- Gotham-specific pgen, dust, cooling-flow, and SN edits are absent.
- Existing inputs continue to run without enabling new behavior.
- N-dimensional PDF output and readers pass focused tests.
- Modern PDF headers and payloads publish atomically per file, sparse shards declare
  complete inventory metadata, and shipped readers reject malformed aliases and
  unreasonable file-controlled allocations.
- Spherical-slice output and reader pass focused tests.
- Shared/per-rank/per-node comparisons pass for supported output formats.
- Per-node restart round trips pass.
- Node-payload hard links, byte copies, and corrupt content markers cannot bypass the
  public-manifest-only restart contract.
- `bin_convert.py` is the only supported binary conversion module, incorporates the
  required modern functionality, and all repository imports/examples use it.
- All shipped reader/helper scripts have valid-file and corrupt/inconsistent-file
  tests appropriate to their format.
- User-facing example workflows execute successfully and validate readback data.
- Output timing is opt-in and reports a useful MPI elapsed maximum.
- Final-output policy is explicit and defaults compatibly.
- Final restart numbering does not overwrite a terminal checkpoint after resume.
- A format compatibility matrix and frozen baseline fixtures establish the promised
  legacy-reader/new-reader and on-disk compatibility boundaries.
- A deferred, target-path-specific `gh-pages` documentation package exists, and no
  premature live `gh-pages` branch modification has been made.
- The staged Pages content has been checked against final code and is accompanied by
  instructions for warnings-as-errors validation after the code branch merges.
- Independent subagent audits of C++ IO/restart behavior, Python tooling, tests and
  examples, and deferred docs integration have been reviewed and addressed.
- The audit ledger identifies blocking findings, their resolutions, and completed
  re-audits for every required stop gate.
- Any GPU-sensitive derived-variable behavior has an explicit validation result or
  merge-blocking validation requirement.

## Open Decisions For Implementation Review

The following details are proposed here but should be confirmed during implementation
review:

| Decision | Proposed Choice | Reason |
| --- | --- | --- |
| Branch name | `feature/io-output-formats-and-sharding` | User-selected feature namespace and scoped name |
| Timing parameter location | `<time>/output_timing` | Runtime-wide diagnostic setting |
| Timing default | `false` | Avoid new overhead/log noise unless requested |
| Timing reduction | MPI maximum | Captures slowest-output cost |
| Final-output parameter location | `<time>/final_output_policy` | Runtime behavior rather than one format |
| Final-output default | `all` | Preserves existing behavior |
| Production policy | `restart_only` | Avoids costly terminal diagnostics |
| Final restart numbering | Consume next number normally | Prevent overwrite after resume |
| Canonical converter module | `vis/python/bin_convert.py` | User uses modern implementation; avoid duplicated APIs |
| Converter migration basis | Fold `bin_convert_new.py` into `bin_convert.py`, preserving required old helpers | Modern functionality with one maintained entry point |
| Example status | Executable and tested, not prose-only | Demonstrates production-facing workflows |
| Documentation timing | Stage in feature branch; integrate into `gh-pages` only after code merge, reapplication, rebuild, and re-audit | Avoid publishing behavior before it exists |
| Documentation validation | Temporary `origin/gh-pages` worktree with Sphinx warnings as errors | Fits live Pages framework safely |
| Review method | Multiple bounded subagent audits plus main-agent verification | High-risk IO work requires independent checks |
| Particle restart enhancements | Separate follow-up | Keeps IO branch focused |

## Definition Of Done

The IO feature branch is complete when it provides a clean, test-backed, documented
implementation of:

- generic diagnostic output fields needed for analysis;
- N-dimensional PDFs with useful scaling and weighting support;
- spherical-slice output and Python reconstruction;
- `single_file_per_node` sharding for supported mesh diagnostics and restarts;
- robust per-node restart loading;
- opt-in output timing appropriate for assessing IO scaling; and
- explicit, restart-safe final-output control;
- one authoritative `bin_convert.py` implementation with complete new-format reader
  support and no unresolved `bin_convert_new.py` duplication;
- executable, validated examples for producing and analyzing the new outputs; and
- extensive deferred documentation content prepared for post-merge integration into
  the `gh-pages` framework.

The branch must deliver these features without importing unrelated Gotham simulation
content, without depending on the rejected attempted extraction branch, and without
introducing silent default behavior changes that risk user data or restart history.
The implementing agent must proceed deliberately, use independent subagents to audit
the critical feature families and integration artifacts, and resolve their findings
before describing the branch as complete.
