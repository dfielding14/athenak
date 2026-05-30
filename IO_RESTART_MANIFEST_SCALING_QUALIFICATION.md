# IO Restart Manifest Scaling Qualification

## Purpose

This record preregisters measurements collected during external `RCP-09`
multi-node qualification for the `RCP-08B` node-restart manifest-validation
disposition. Run it before deciding whether replicated validation should remain
unchanged or move to a separate scaling branch.

The current implementation deliberately validates the complete manifest and
replicated payload headers on every rank. That path is correct and locally
tested. Do not refactor it before collecting scheduler-backed evidence.

## Opt-In Instrumentation

Set:

```bash
export ATHENAK_RESTART_MANIFEST_TIMING=1
```

Node-manifest resume then emits one structured record per rank for each phase:

```text
[restart-manifest] phase=validate rank=... estimate_scope=manifest_parse_plus_replicated_header_compare payload_count=... manifest_bytes=... header_bytes=... logical_header_validation_bytes=... logical_validation_pressure_bytes=... logical_validation_pressure_capped=... elapsed_s=...
[restart-manifest] phase=startup_parse ...
[restart-manifest] phase=load_local_blocks ...
```

The flag is default-off. It changes diagnostics only. It does not alter restart
metadata, read routing, collectives, or output files.

The `phase=validate` byte fields are a scoped logical validation-pressure
estimate, not observed filesystem IO. For a valid manifest with `P` payloads,
manifest size `M`, and replicated header size `H`, each rank reports:

```text
logical_header_validation_bytes = 2 * (P - 1) * H
logical_validation_pressure_bytes = M + logical_header_validation_bytes
```

This scope counts one logical manifest parse plus both sides of each replicated
payload-header comparison in `NodeRestartManifest::Load()`. It excludes the
signature probe before `Load()`, metadata-only payload-size checks, page-cache
effects, filesystem read-ahead, storage-server behavior, and subsequent local
MeshBlock reads. `logical_validation_pressure_capped=1` means diagnostic-only
arithmetic exceeded `UINT64_MAX` and the byte fields were capped at that value;
restart correctness is not allowed to fail because an estimate overflowed.

The logical estimate is mandatory telemetry, but it does not close the
filesystem-read-amplification row. External qualification must also collect an
attributable observed filesystem-read measurement from the target platform.
Use supported scheduler, process, or filesystem accounting and retain the raw
commands and logs. If the platform cannot expose a reliable attributable
measurement, mark that row incomplete and record the limitation. Do not
substitute logical pressure for observed filesystem reads.

## Preregistered Topologies

Run every required topology on distinct physical hosts. Capture the scheduler
job ID and `hostname` for every rank.

| Tier | Nodes | Ranks per node | Total ranks | Required |
| --- | ---: | ---: | ---: | --- |
| correctness-small | 2 | 2 | 4 | Yes |
| scaling-low | 4 | 2 | 8 | Yes |
| scaling-dense | 4 | 8 | 32 | Yes |
| production-representative | Operator must record | Operator must record | Operator must record | Yes, when it differs materially from the rows above |

Run five measured resumes per topology after one unmeasured warm-up resume.
Record all five samples. Do not discard an outlier without preserving it and
explaining the exclusion.

## Frozen Correctness Decks

Run these decks before the representative scaling case. Preserve the exact
launcher command, rank-to-host mapping, manifest, stdout, stderr, and recursive
output inventory.

### ED-1: Empty sliced-diagnostic node

Use `tst/inputs/io_node_sharding.athinput` with two physical nodes and four ranks.
Place world ranks `0,1` on host A and ranks `2,3` on host B. Apply:

```text
mesh/nx1=32
output1/single_file_per_node=true
output2/single_file_per_node=true
output3/single_file_per_node=true
output4/single_file_per_node=true
output5/single_file_per_node=true
output6/single_file_per_node=true
time/final_output_policy=none
```

The frozen mesh has four x1 MeshBlocks. With contiguous world-rank placement,
each physical node owns two MeshBlocks. For sliced `.bin` `output2` at
`slice_x1=0.25`, require these exact shard files:

```text
bin/node_00000000/io_node.slice.00000.bin
bin/node_00000001/io_node.slice.00000.bin
```

The node-0 shard must declare zero selected MeshBlocks, the node-1 shard must
declare one selected MeshBlock, and canonical reader assembly must equal a
shared-output run with the same `mesh/nx1=32` override. Preserve the full
incidental output inventory, but use this bounded sliced-stream inventory as
the ED-1 pass/fail check.

### MR-1: Balanced minimal restart resume

Generate a checkpoint from `tst/inputs/io_node_sharding.athinput` with the ED-1
rank placement and:

```text
mesh/nx1=32
output1/dt=-1
output2/dt=-1
output3/dt=-1
output4/dt=-1
output5/dt=-1
output6/single_file_per_node=true
time/final_output_policy=none
```

Require one public manifest, two node payloads, `payload_count=2`, payload block
inventory `[2, 2]`, and four one-block segment records. The exact bounded
inventory is:

```text
rst/io_node.00000.rst
rst/node_00000000/io_node.00000.g<G>.payload.rst
rst/node_00000001/io_node.00000.g<G>.payload.rst
```

`<G>` must be one identical unsigned generation token in both payload paths.
Resume only through the public manifest with:

```text
time/tlim=0
time/nlim=0
time/final_output_policy=none
time/output_timing=false
```

Set `ATHENAK_RESTART_MANIFEST_TIMING=1` only for measured resumes. Require one
record per rank for each registered phase and no `.assembled` or
`.assembled.tmp` file.

### MR-2: Minimal imbalanced restart resume

Use two physical nodes and three ranks. Place world rank `0` on host A and ranks
`1,2` on host B. Generate and resume as in MR-1, except use `mesh/nx1=24`.
Require the same bounded manifest/payload path inventory, `payload_count=2`,
payload block inventory `[1, 2]`, and three one-block segment records. This is
the frozen smallest-practical-node-owned-payload deck.

## Preregistered Payload Targets

Measure two checkpoint families:

| Family | Target | Purpose |
| --- | --- | --- |
| correctness | `tst/inputs/io_node_sharding.athinput` at SHA-256 `92f54bb10841f428fa3942fbecbae2f9b455e458661aded78729261bf67e3ff8` | Prove layout, empty/non-owning-node behavior, direct resume, and log collection. Do not substitute the promoted example deck for this frozen qualification input. |
| representative | At least `256 MiB` payload bytes per writing node and at least `4 MiB` replicated header bytes when the production mesh can reach that header size | Measure logical validation pressure and startup cost under a meaningful payload inventory. |

If the production case cannot reach the representative targets, record the
largest practical checkpoint and the limiting resource. Preserve the observed
payload and header sizes in the evidence packet.

## Required Environment Record

Complete this table before consuming scheduler time:

| Field | Record |
| --- | --- |
| Cluster and scheduler | Pending external operator |
| Parallel filesystem and mount | Pending external operator |
| AthenaK branch SHA | Pending external operator |
| Frozen robustification-guide SHA-256 | `cdf53351104c135f2d79f9ee2f36c2908001b00a31339db18eb1cf66bcae11ff` |
| Compiler and version | Pending external operator |
| MPI implementation and version | Pending external operator |
| Kokkos configuration | Pending external operator |
| Scheduler job IDs | Pending external operator |
| Rank-to-host mapping command | `srun --label hostname` or scheduler-equivalent |
| Rank-separated log template | `--label --output=<packet>/logs/<row>.<action>.<attempt>.<sample>.%J.%t.out --error=<packet>/logs/<row>.<action>.<attempt>.<sample>.%J.%t.err` or scheduler-equivalent; Slurm `%J` includes job and step IDs |
| Scheduler launcher | Pending external operator |
| Scheduler allocation walltime | Pending external operator |
| Per-row timeout | `300 s` default; operator must record any justified override |
| Observed filesystem-read command | Pending external operator; must be attributable to the measured resume |
| Repeat count | One warm-up plus five measured resumes per topology and checkpoint family |
| Reported percentiles | median, p95, and maximum |

## Aggregation And Startup Endpoint

Retain every per-rank structured record. For each measured resume:

1. define `validate_run_max_s` as the maximum per-rank `phase=validate
   elapsed_s`;
2. define `startup_parse_run_max_s` as the maximum per-rank
   `phase=startup_parse elapsed_s`;
3. define `logical_validation_pressure_run_bytes` as the sum of per-rank
   `logical_validation_pressure_bytes`;
4. retain the per-rank values as well as those run-level aggregates; and
5. reject the evidence packet if any representative record has
   `logical_validation_pressure_capped=1`.

For each topology and checkpoint family, compute median, p95, and maximum over
the five measured run-level samples. Use the nearest-rank empirical percentile:
sort ascending and select index `ceil(0.95 * N) - 1`; with the preregistered
`N=5`, p95 is the maximum. Do not pool ranks across runs and do not interpolate.

Measure `external_minimal_resume_wall_s` separately around every measured
representative no-step resume launcher from inside an existing scheduler
allocation so queue wait is excluded. For each measured run, start a monotonic
external timer immediately before invoking that run's scheduler launcher and
stop it immediately after the launcher returns success after all ranks exit.
The frozen `time/tlim=0`, `time/nlim=0`, and
`time/final_output_policy=none` overrides make this a no-step minimal-resume
endpoint. It is an operational startup/control-plane proxy, not a claim that
the launcher interval isolates filesystem service time. For each run compute
`validate_startup_fraction = validate_run_max_s /
external_minimal_resume_wall_s`, then apply the same five-sample p95 rule.

## Required Measurements

Retain raw stdout and stderr for every run. For each rank and run, capture:

1. payload count;
2. rank count;
3. node count from the scheduler rank-to-host mapping;
4. replicated header bytes;
5. manifest bytes;
6. per-rank `phase=validate elapsed_s`;
7. aggregate scoped logical validation pressure, computed by summing
   `logical_validation_pressure_bytes` across ranks;
8. `phase=startup_parse startup_before_s`, `startup_after_s`, and `elapsed_s`;
9. `phase=load_local_blocks elapsed_s`, `local_bytes`, and `span_count`;
10. filesystem, scheduler, compiler, MPI, and Kokkos environment;
11. `external_minimal_resume_wall_s` from the operational endpoint above; and
12. confirmation that no `.assembled` or `.assembled.tmp` file appears.
13. attributable observed filesystem-read bytes for the measured resume, the
    command used to collect them, and raw accounting logs; and
14. observed read amplification, defined as attributable filesystem-read bytes
    divided by the checkpoint bytes required by one logical resume.

## Executable Scheduler Templates

Run every correctness and scaling row with rank-separated logs. The checked-in
Slurm runner is the normative executable path; preserve an equivalent script on
other schedulers:

```bash
export PACKET=/durable/path/to/packet
export REPO=/path/to/athenak
export BUILD=/path/to/mpi-build
export PYTHON=/path/to/python
export ROW_TIMEOUT="${ROW_TIMEOUT:-300}"
export ATTEMPT_ID=attempt-1
export FS_ACCOUNTING_HOOK=/site/path/to/attributable-read-accounting-hook

"$REPO/scripts/run_external_io_qualification_slurm.sh" ED-1 generate
"$REPO/scripts/run_external_io_qualification_slurm.sh" MR-1 generate
"$REPO/scripts/run_external_io_qualification_slurm.sh" MR-1 resume "$MR1_MANIFEST"
export MR2_HOSTFILE="$PACKET/decks/mr2.hostfile"
"$REPO/scripts/run_external_io_qualification_slurm.sh" MR-2 generate
"$REPO/scripts/run_external_io_qualification_slurm.sh" MR-2 resume "$MR2_MANIFEST"
export ARCHIVE_RECORD=/durable/path/to/io-qualification-archive.tsv
"$REPO/scripts/finalize_external_io_qualification_packet.sh" \
  "$PACKET" "$ARCHIVE_RECORD"
```

The runner fixes ED-1 and MR-1 at two nodes with two ranks per node. Its
three-line `MR2_HOSTFILE` must list canonical hostname-token host A once and
canonical hostname-token host B twice, forcing the MR-2 world-rank placement.
Resume rows execute one
unmeasured warm-up plus five measured samples. Only measured samples enable
`ATHENAK_RESTART_MANIFEST_TIMING=1`. Every launcher interval is measured inside
one timer process and published through revalidated packet-root and
log-directory descriptors; every rank log path includes row and sample
identity; every
exit code and timeout disposition lands in `packet-index.tsv`; packet-local
attempt IDs and duplicate-key rejection protect retries; every Athena launch
retains a forbidden-staging scan; and every measured sample invokes the
site-filled accounting hook before and after the launcher. Resume scans include
the packet and manifest directory and reject both `*.assembled` and
`*.assembled.tmp`. If the hook is unset, the observed-read row remains
explicitly incomplete. Recursive-inventory and forbidden-staging scan failures
are retained and terminally indexed. Rank-map tasks retain one record per rank;
the runner requires regular single-link rank-map records and validates the
observed ED-1, MR-1, MR-2, or scaling topology before indexing success. Before
deck copying and every scheduler launch, it fail-closed scans the packet tree
for symlinks, hard-link aliases, control-character paths, and scanner failures;
it repeats that scan immediately after every scheduler return. The retained
rank-map inventory is published without following an existing target. Measured
resume rows require the site-filled accounting hook unless they terminally
record hook unavailability. Packet-tree admission repeats after each accounting
hook execution and requires the packet-root identity pinned at launcher
startup. Candidate TSV rows are validated before append. The runner also pins
its packet-local index inode and fixed evidence directories for the invocation,
plus each rank-map and launch-output directory after creation, and revalidates
those identities after scheduler and hook transitions. The packet finalizer
is idempotent and retryable:
it validates the runner lifecycle, rejects existing or dangling archive-record
symlinks and conflicting outer claims for one canonical packet path, and
publishes metadata through packet-local exclusive temporaries, file sync,
same-directory atomic rename, and packet-directory sync. It re-syncs admitted
retained metadata and resets writable inconsistent metadata pairs on retry.
Retained packet-local reserved temporary paths and archive-adjacent
destination-scoped replacement candidates fail closed for explicit operator
adjudication rather than being deleted from filename shape alone. It verifies
checksums, removed write
permissions, and a regenerated canonical inventory before appending the
external outer-index archive record. It repeats packet-tree alias and
reserved-path admission after permission removal and before durable
publication. Packet-index and archive-row appends use checked complete writes
to sibling temporaries, file sync, same-directory atomic replacement, and
parent-directory sync. Archive-sink admission always re-syncs its parent.
The descriptor-bound primitive requires the archive-directory, archive-sink,
and packet identities admitted at finalizer startup. The append is the helper's
last substantive operation. From the first runner invocation through
finalization, use a trusted environment with no concurrent or inter-invocation
same-owner packet or archive mutator and do not replace packet-local children.
Owner write-bit removal is not a defense against an owner deliberately
restoring permissions after the helper returns.

For `ED-1`, generate outputs from `tst/inputs/io_node_sharding.athinput` with
the exact ED-1 overrides, then assemble the sliced `.bin` shard with the
canonical reader and compare it numerically against the corresponding shared
run. For `MR-1` and `MR-2`, generate the exact preregistered manifest family,
run the direct manifest resume template above, and check the declared payload
and segment inventory before accepting the row. Record all substitutions,
launcher exit codes, and timeout dispositions in the packet index.

## Empty And Non-Owning Node Arrangement

ED-1 is the required physical-node-with-no-selected-sliced-records deck. MR-2 is
the required smallest-practical-node-owned restart payload deck. Record the
observed inventories and stop qualification if they differ from the frozen
expectations.

## Decision Thresholds

Use the following preregistered thresholds for `RCP-08B`:

| Result | Decision |
| --- | --- |
| p95 `validate_run_max_s` is no more than `2 s` and p95 `validate_startup_fraction` is no more than `5%` for every representative tier | Keep replicated validation in this branch and document measured acceptability. |
| p95 `validate_run_max_s` exceeds `5 s` or p95 `validate_startup_fraction` exceeds `10%` on any representative tier | Open a separate scaling branch to evaluate rank-0 parse plus structured broadcast and node-leader header validation. Do not hide that protocol change inside this feature branch. |
| Result lies between those thresholds | Preserve the current protocol in this branch, record the evidence, and require an explicit follow-up decision before production rollout. |

Correctness failures always block rollout regardless of timing.
Unavailable or unattributable observed filesystem-read measurements leave the
filesystem-amplification row incomplete even when logical-pressure and timing
thresholds pass.

## Candidate Alternatives For RCP-08B

Record one final disposition after measurements:

1. keep replicated validation;
2. parse and validate on rank 0, then broadcast structured metadata;
3. validate payload headers on node leaders and distribute node-local metadata;
4. open a separate scaling branch because the optimization is material but too
   risky for this feature branch.

Compare correctness risk, implementation complexity, measured performance, and
required requalification. Any protocol refactor must rerun malformed-manifest,
header-mismatch, missing-payload, marker, traversal, alias, changed-rank-count,
multi-node, empty/non-owning-node, forced-chunking, collective-error, and
no-`.assembled` or `.assembled.tmp` checks.

## RCP-08A Reflection

The current direction remains measurement-first. Per-rank strict validation is
easy to reason about and already covered locally. A centralized protocol might
reduce logical validation pressure, but without real topology data it would
expand the restart blast radius based on speculation.
