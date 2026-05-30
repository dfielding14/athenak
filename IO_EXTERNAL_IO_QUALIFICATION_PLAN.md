# External IO Qualification Plan

## Status

This workstation cannot close `RCP-09`: it has no CUDA or HIP toolchain, no
scheduler launcher, and no second physical host. This record defines the
external evidence required before the branch can be called merge-ready.

## Frozen Inputs

| Field | Value |
| --- | --- |
| Feature branch | `feature/io-output-formats-and-sharding` |
| Final tested local branch SHA | Pending: record the exact committed robustification SHA before external qualification |
| Robustification guide SHA-256 | `cdf53351104c135f2d79f9ee2f36c2908001b00a31339db18eb1cf66bcae11ff` |
| Guide size | `2222` lines, `89840` bytes |

Record the exact tested branch SHA after local robustification commits land.

## Immutable Evidence Packets

Create one immutable evidence-packet directory per qualification run family.
Do not treat scheduler scratch paths or terminal summaries as durable evidence.
Each packet must retain:

1. a copied input deck and any runtime overrides;
2. exact working directory, environment variables, build commands, execution
   commands, reader commands, and failure-injection commands;
3. exit codes and scheduler job IDs;
4. host, device, compiler, MPI, Kokkos, scheduler, and filesystem inventory;
5. raw scheduler, test, timing, and reader logs;
6. generated manifests, output inventories, and recursive file listings;
7. numerical-comparison results;
8. a packet index connecting each reported result to its raw evidence; and
9. an `artifacts.sha256` manifest covering every retained packet artifact except
   itself and the packet index.

Archive the packet root somewhere durable before closing a qualification row.
Record that archive path, the `artifacts.sha256` digest, and an outer SHA-256
digest of the completed packet index in the evidence template below. The packet
index must name `artifacts.sha256`, record its digest, and state that the inner
manifest intentionally excludes itself and the index. The ledger or archive
record holds the outer index digest. Do not create a self-referential checksum
claim.

Finalize a completed packet with the checked-in helper. Keep the archive record
outside the packet so the packet can become read-only after the outer digest is
retained:

```bash
export ARCHIVE_RECORD=/durable/path/to/io-qualification-archive.tsv
"$REPO/scripts/finalize_external_io_qualification_packet.sh" \
  "$PACKET" "$ARCHIVE_RECORD"
```

The finalizer is idempotent for completed packets and retryable after admitted
publication failures. It rejects packet aliases, existing
or dangling archive-record symlinks, archive-record hard-link aliases,
control-character paths, structurally invalid runner histories, and conflicting
outer claims for one canonical packet path. It publishes metadata through a
packet-local exclusive temporary, file sync, same-directory atomic rename, and
packet-directory sync; re-syncs admitted retained metadata on retry; and resets
a writable inconsistent metadata pair on retry. Retained packet-local reserved
temporary paths and archive-adjacent destination-scoped replacement candidates
fail closed for explicit operator adjudication rather than being deleted from
filename shape alone. It
verifies inner checksums; removes and verifies packet write permissions; repeats
packet-tree alias and reserved-path admission; verifies inner checksums again;
regenerates the canonical inventory to reject late artifacts; and only then
appends the durable
outer-index record through an alias-safe descriptor-bound primitive that
requires archive-directory, archive-sink, and packet identities admitted at
finalizer startup. Packet-index and archive-row appends publish a completely
written and synced sibling temporary through same-directory atomic replacement,
then sync the parent directory. Archive-sink admission always re-syncs its
parent directory. The final append is the helper's last substantive operation.
A retry after interrupted metadata, permission, or archive-record publication
resumes from admitted retained state or fails closed on an ambiguous
packet-local or archive-adjacent reserved temporary for operator adjudication.
Run the entire first-runner-invocation-through-finalization lifecycle from a
trusted environment with no concurrent or inter-invocation same-owner packet or
archive mutator. Removing owner write bits is an operator-visible immutability
transition, not a defense against an owner deliberately restoring permissions
after the helper returns.

## CUDA Qualification

Use a CUDA-capable AthenaK build. Retain:

1. host and selected GPU inventory;
2. compiler, CUDA toolkit, MPI, and Kokkos versions;
3. exact CMake and build commands;
4. exact test commands and exit codes;
5. generated output inventory;
6. warnings; and
7. an independent evidence audit.

At minimum run:

```bash
cd "$CUDA_BUILD/src"
"$PYTHON" -m pytest -p no:cacheprovider -q \
  "$REPO/tst/test_suite/io/test_output_formats_gpu.py"
```

Also exercise representative Hydro and MHD modern PDFs, derived coordinate and
velocity axes, scalar weighting, and the adaptive-pack-growth derived-PDF
regression on the actual device build. CPU execution of that module is smoke
coverage only.

Execute the spherical-slice producer on the device backend as well:

```bash
cd "$CUDA_BUILD/src"
"$PYTHON" -m pytest -p no:cacheprovider -q \
  "$REPO/tst/test_suite/io/test_output_formats_cpu.py" -k sphslice
"$PYTHON" -m pytest -p no:cacheprovider -q \
  "$REPO/tst/test_suite/io/test_output_formats_mpicpu.py" -k sphslice
```

These commands intentionally run the producer modules from the CUDA build
directory. Retain evidence for shared, rank, and node layouts; the analytical
cross-face oracle; mixed-level AMR rebuild against the binary ghost-snapshot
oracle; and finite-double-to-serialized-float overflow rejection. A HIP
deployment must repeat the same spherical-slice rows in its separate HIP
packet.

Preregister one explicit MHD device command and retain its oracle:

```bash
cd "$CUDA_BUILD/src"
./athena -i "$REPO/inputs/mhd/blast_mhd.athinput" -d "$PACKET/run_mhd_pdf" \
  output1/file_type=pdf output1/id=device_mhd \
  output1/variable_1=coord_r output1/nbin1=16 \
  output1/bin1_min=1.0e-3 output1/bin1_max=1.0 output1/scale1=log \
  output1/variable_2=edot_sph output1/nbin2=32 \
  output1/bin2_min=-1.0e2 output1/bin2_max=1.0e2 output1/scale2=symlog \
  output1/linthresh2=1.0e-2 output1/weight=volume
"$PYTHON" "$REPO/vis/python/examples/read_io_outputs.py" pdf \
  "$PACKET/run_mhd_pdf/pdf_device_mhd_coord_r_edot_sph/"*.pdf
```

Retain the reader output, require finite histogram values, and record a
nonempty finite histogram sum. If the concrete target MHD deck changes before
execution, record the replacement deck and rationale in the packet index.

The frozen robustification guide requires CUDA qualification. HIP deployment is
not qualified by that lane. If this branch will be deployed on a HIP backend,
run and retain a separate HIP evidence packet with the same representative
device regressions, environment inventory, and independent audit before making
a HIP production-readiness claim.

## Multi-Node Qualification

Use at least two physical nodes with more than one rank per node on the target
parallel filesystem. Capture scheduler job ID and rank-to-host mapping.

Run and retain evidence for:

1. full-volume node-sharded `.bin`;
2. sliced node-sharded `.bin`;
3. uniform 3D active-zone full-volume node-sharded `.cbin`;
4. node-sharded modern PDF;
5. node-sharded `sphslice`;
6. node-sharded restart publication;
7. direct manifest restart resume;
8. changed rank count;
9. changed rank distribution across nodes where feasible;
10. a genuinely empty or non-owning node for diagnostic shard families;
11. output timing labels;
12. expected and observed shard inventories;
13. canonical-reader assembly and equality checks;
14. no `.assembled` or `.assembled.tmp` staging;
15. production-filesystem rename, unlink, and cleanup behavior;
16. injected post-payload-write and post-publication restart rollback;
17. the measurements preregistered in
    `IO_RESTART_MANIFEST_SCALING_QUALIFICATION.md`; and
18. an independent topology-evidence audit.

Explicit empty diagnostic shards are required for `.bin`, `.cbin`, modern PDF,
and `sphslice` where the selected output may leave a node without records.
Restart payload inventories remain positive by contract; use the preregistered
imbalanced-positive `MR-2` restart case instead of an empty restart payload.

## Scheduler Runbook

Before consuming an allocation, record these packet-local values:

```bash
export PACKET=/durable/path/to/packet
export REPO=/path/to/athenak
export BUILD=/path/to/mpi-build
export PYTHON=/path/to/python
export ROW_TIMEOUT=300
export ATTEMPT_ID=attempt-1
mkdir -p "$PACKET"/{decks,logs,outputs,inventory}
```

Request enough walltime for one warm-up plus five measured resumes per scaling
row, correctness decks, reader checks, and failure injection. Record the
requested walltime and actual scheduler allocation. Use rank-separated
stdout/stderr. For Slurm, use the checked-in runner:

```bash
"$REPO/scripts/run_external_io_qualification_slurm.sh" ED-1 generate
"$REPO/scripts/run_external_io_qualification_slurm.sh" MR-1 generate
"$REPO/scripts/run_external_io_qualification_slurm.sh" MR-1 resume "$MR1_MANIFEST"
export MR2_HOSTFILE="$PACKET/decks/mr2.hostfile"
"$REPO/scripts/run_external_io_qualification_slurm.sh" MR-2 generate
"$REPO/scripts/run_external_io_qualification_slurm.sh" MR-2 resume "$MR2_MANIFEST"
```

The script fixes ED-1 and MR-1 at two nodes with two ranks per node. It requires
an explicit three-line `MR2_HOSTFILE` for MR-2: canonical hostname-token host A
once, then canonical hostname-token host B twice.
It records unique row/action/attempt/sample logs with Slurm job-step IDs,
launcher exit codes, timeout dispositions, single-process monotonic launcher
intervals published through revalidated packet-root and log-directory
descriptors, retained
forbidden-staging scans, and a site-filled
`FS_ACCOUNTING_HOOK`. Duplicate packet-local launch keys are rejected. Resume
scans include the packet and manifest directory and reject both `*.assembled`
and `*.assembled.tmp`. An unset hook leaves the observed-filesystem-read row
explicitly incomplete. Before deck copying and before every scheduler launch,
the runner fail-closed scans the packet tree for symlinks, hard-link aliases,
control-character paths, and scan failures. It repeats that scan immediately
after each scheduler return. Observed rank-map files must be regular single-link
files, and the retained topology inventory is created without following an
existing publication target. Measured resume rows require an attributable-read
hook unless they terminally record that the hook was unavailable. Packet-tree
admission repeats after each accounting-hook execution and requires the
packet-root identity pinned at launcher startup. Candidate TSV rows are
validated before append. The runner also pins the packet-local index inode and
fixed `decks/`, `logs/`, `outputs/`, and `inventory/` directories for its
invocation, plus each rank-map and launch-output directory after creation. It
revalidates those identities after scheduler and accounting-hook transitions.
Do not replace packet-local children between runner invocations or between the
last runner invocation and finalization.
Use
scheduler-equivalent scripts on a non-Slurm cluster only if they retain and
validate the observed rank map and terminally index recursive-inventory or
forbidden-staging scan failures. Preserve the exact substitutions. For `ED-1`,
`MR-1`, and `MR-2`, use the exact overrides and
oracles preregistered in
`IO_RESTART_MANIFEST_SCALING_QUALIFICATION.md`. After each output row, retain
the recursive file inventory and run the applicable canonical Python reader
with `assemble_shards=True`. For restart rows, retain direct public-manifest
resume commands and prove that no `.assembled` or `.assembled.tmp` file exists.

Failure-injection rows must distinguish precommit and postcommit behavior.
Inject failures after payload publication, after public-manifest publication,
and during generation-reservation removal. Precommit failures must not leave a
published manifest. Once manifest rename succeeds, the checkpoint is committed:
postcommit failures may return a failure status but must preserve the public
manifest and its declared payloads as a resumable checkpoint.

Record a per-row timeout disposition. A timeout is a failed or incomplete row,
never a silently discarded sample.

## Evidence Template

```markdown
### External Qualification: Short Title

| Field | Record |
| --- | --- |
| Environment | Host, scheduler, filesystem, compiler, selected device backend, MPI, Kokkos |
| Deployment backend | CUDA, HIP, or explicit `N/A` justification |
| Branch snapshot | Exact SHA |
| Guide snapshot | Exact SHA-256 checksum and line count |
| Scheduler evidence | Job ID plus physical rank-to-host mapping |
| Device evidence | CUDA or HIP inventory and selected device when applicable |
| HIP packet status | Immutable packet path and digest, or explicit `N/A` justification |
| Build | Exact commands |
| Runtime topology | Nodes, ranks, ranks per node, empty/non-owning-node arrangement |
| Commands | Exact execution and reader commands with exit codes |
| Evidence packet | Immutable archive root plus packet-index path |
| Packet checksums | `artifacts.sha256` path and digest plus outer packet-index digest |
| Checksum rule | Inner manifest excludes itself and the index; packet index records the inner digest; ledger or archive record holds the outer index digest |
| Retained logs | Paths to scheduler, test, and timing logs |
| Output inventory | Expected and observed files |
| Results | Test results and numerical comparisons |
| Timing | Relevant measurements |
| Observed filesystem reads | Attributable measurement command, raw logs, per-rank or per-job values, and amplification result |
| Logical validation pressure | Structured AthenaK estimate; never substituted for observed filesystem reads |
| Timeout policy | Requested walltime, row timeout, and any timeout disposition |
| Working tree | State before and after validation |
| Failures | Any failure and disposition |
| Auditor | Independent qualification reviewer |
| Status | Passed, failed, or incomplete |
```

## Stop Rule

Do not describe the branch as merge-ready or production-qualified until the
required CUDA and multi-node lanes pass and their independent auditors accept
the evidence. Do not make a HIP production-readiness claim unless the
deployment backend is explicitly `N/A` or a separate HIP packet passes.
