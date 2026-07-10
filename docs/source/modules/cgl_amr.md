# CGL Adaptive Mesh Refinement

## Scope

AthenaK supports adaptive and static mesh refinement for nonrelativistic CGL
MHD, including the Landau-fluid (LF) heat-flux closure integrated with RKL2
super-time-stepping (STS). The implementation has been exercised with uniform
refinement, static coarse/fine interfaces, repeated adaptive refinement and
derefinement, restart across regrids, and multi-node MPI+GPU execution.

This guide explains the CGL-specific AMR algorithm, the configuration that has
been validated, the tests supporting it, and how to run the supplied
production-scale example. For the CGL equations and LF closure, first see
[CGL Method Physics Primer](cgl_mhd_method.md) and
[CGL Landau-Fluid Heat Flux](cgl_landau_fluid.md).

The strongest validation applies to Cartesian, periodic, nonrelativistic CGL
MHD using PLM reconstruction, the HLLE solver, primitive AMR prolongation, and
LF+RKL2 STS. Other combinations may share the same implementation but should
be qualified separately.

## Why CGL Needs A Specialized AMR Path

Outside an LF sweep, the sixth CGL conserved slot, `IAN`, stores conserved
pressure anisotropy. Primitive arrays instead store parallel and perpendicular
pressure. During an LF sweep, `IAN` temporarily stores magnetic moment. A
generic AMR operation cannot safely interpolate, restrict, or consume this
slot without knowing its current representation.

CGL AMR therefore separates the state into two groups:

- **Primary conserved state:** density, momentum, total energy, and passive
  scalar densities. Generic conservative AMR transfer remains authoritative
  for these fields.
- **CGL thermodynamic encoding:** pressure anisotropy or magnetic moment in
  `IAN`, depending on the active operator. CGL-aware completion reconstructs
  this slot from an admissible thermodynamic state and the correct magnetic
  field.

This separation is the central conservation invariant: ordinary CGL
projection must not rewrite a valid primary conserved state merely to make the
derived anisotropy encoding admissible.

## AMR Algorithm

### Active refinement

When a block is created at a finer level:

1. Generic `RefineCC` prolongs density, momentum, total energy, and scalar
   densities conservatively.
2. Generic face-field prolongation constructs the fine magnetic field.
3. The CGL completion decodes the parent thermodynamics, prolongs and limits
   the pressure difference, and checks each child for admissibility.
4. `IAN` is rebuilt for each child using the already-prolonged primary state
   and magnetic field.
5. Primary fields are changed only if a generic child state is itself invalid;
   such fallbacks are separately counted as AMR repairs.

For an admissible parent with no primary repair, child averages recover the
parent density, momentum, energy, and scalar densities to floating-point
roundoff.

### Restriction and derefinement

Restriction follows the complementary ordering:

1. Primary conserved fields are volume-averaged with generic conservative
   restriction.
2. Face-centered magnetic fields are restricted and finalized through the
   constrained-transport (CT) path.
3. CGL thermodynamics are projected against that finalized coarse magnetic
   field.
4. The coarse `IAN` encoding is rebuilt without altering valid primary
   fields.
5. After any face-field repair, boundary data and primitives are refreshed so
   the next reconstruction never consumes stale ghosts.

The live primitive state used by corner electric fields and CT is not mutated
early by the restriction helper. This preserves the normal MHD task ordering
and prevents CGL bookkeeping from changing the magnetic update.

### LF and STS representation

LF stages operate while `IAN` stores magnetic moment. AMR boundary operations
on that path explicitly require the magnetic-moment representation. At the end
of the split LF sweep, AthenaK converts back to conserved anisotropy before
ordinary hyperbolic evolution, output, or restart.

Representation guards fail fast if an ordinary physical or user boundary
would interpret the temporary LF representation as conserved anisotropy.
Periodic boundaries are covered by the production validation described below;
custom boundaries require their own qualification.

### Magnetic fields and diagnostics

AMR face-field transfer and repair preserve the CT representation. Nonfinite
face fields are fatal rather than silently converted into a recoverable CGL
state. Repair counters are cumulative, reduced across MPI ranks, and preserved
in restart files.

Derived output buffers are resized whenever AMR changes the number of local
MeshBlocks. Repeated outputs such as `mhd_divb` are therefore safe across
refinement, derefinement, and load balancing.

## Source Map

| Responsibility | Primary source |
| --- | --- |
| Generic refinement, restriction, redistribution, and CGL completion | `src/mesh/mesh_refinement.cpp` |
| CGL admissibility intervals and repair masks | `src/eos/cgl_amr_projection.hpp` |
| CGL primitive recovery and representation checks | `src/eos/cgl_mhd.cpp` |
| MHD AMR and CT task ordering | `src/mhd/mhd_tasks.cpp` |
| LF representation lifecycle and STS stages | `src/diffusion/cgl_landau_fluid.cpp`, `src/diffusion/sts_rkl2.cpp` |
| CGL-aware boundary behavior | `src/bvals/physics/bfield_bcs.cpp` |
| AMR-derived fields and output resizing | `src/mesh/refinement_criteria.cpp`, `src/outputs/derived_variables.cpp` |
| Focused AMR problem generator and repair histories | `src/pgen/tests/divb_amr.cpp` |

## Configuration

### Common CGL and LF settings

Production CGL+LF AMR runs should start from the following constraints:

```ini
<time>
integrator = rk2
sts_integrator = rkl2

<mhd>
eos = cgl
passive = false
reconstruct = plm
rsolver = hlle
cgl_heat_flux = landau_fluid
cgl_heat_flux_integrator = sts
cgl_lf_strict_admissibility = true
cgl_lf_record_pressure_work = false
```

`cgl_lf_record_pressure_work` must remain disabled with primitive AMR
prolongation because that diagnostic communication path has not been audited
for this representation lifecycle.

The fast production LF mode is compatible with the tested AMR path:

```ini
<mhd>
cgl_lf_diagnostics = none
cgl_lf_arithmetic = fast
cgl_lf_sts_flux = physical
cgl_lf_profile = false
cgl_lf_profile_detail = false
```

Strict safety counters remain active when detailed LF diagnostics are disabled.

### Static or uniform refinement

For a static refined region, enable primitive prolongation and specify one or
more refined regions:

```ini
<mesh_refinement>
refinement = static
prolong_primitives = true
max_nmb_per_rank = 128

<refined_region1>
level = 1
x1min = -0.5
x1max = 0.5
x2min = -0.5
x2max = 0.5
x3min = -0.6666666666666666
x3max = 0.6666666666666666
```

Using the full domain as the refined region gives uniform refinement through
the AMR transfer machinery. Add a second full-domain region at `level = 2` for
two additional levels. If no coarse/fine interface or AMR transfer is needed,
setting the root mesh directly to the target resolution is simpler and avoids
AMR overhead.

### Adaptive current-based refinement

The validated production workflow tags blocks using the maximum
cell-centered current magnitude, `|curl B|`:

```ini
<mesh_refinement>
refinement = adaptive
num_levels = 3
ncycle_check = 8
refinement_interval = 16
prolong_primitives = true
max_nmb_per_rank = 64

<amr_criterion1>
method = min_max
variable = mhd_current
value_max = 14.416394851879101
```

`num_levels = 3` means the root level plus two additional physical refinement
levels. Blocks whose maximum current exceeds `value_max` are refined; blocks
below it are candidates for derefinement, subject to the mesh nesting rules and
the refinement interval.

The numerical threshold above is specific to the paired production-scale
decks supplied with this guide. For a different base state, resolution, box,
or normalization, select a new threshold from the intended physical
criterion. Do not reuse this value merely because the mesh dimensions match.

Choose `max_nmb_per_rank` with headroom. The validated run ended at 62--63
blocks per rank with a cap of 64; a changed threshold or node count can require
a larger cap.

## Supplied Production-Scale Workflow

Two input files reproduce the accepted 3-D workflow:

- `inputs/cgl_lf_paper/cgl_lf_paper_amr_192x192x256_base.athinput`
- `inputs/cgl_lf_paper/cgl_lf_paper_amr_192x192x256_current_l2.athinput`

The physical box is `1 x 1 x 4/3`, the root grid is `192 x 192 x 256`, and
MeshBlocks contain `48 x 48 x 64` cells. The base deck evolves the unrefined
mesh from `t=0` to `t=1`. Its restart output is shared so the continuation can
change MPI rank count. The AMR overlay then evolves from `t=1` to `t=2` with
two additional current-selected levels and rank-local checkpoints.

Run the stages in separate, initially empty directories using the normal MPI
or scheduler launcher for the machine. For example, define the source and
executable paths, then run:

```bash
src=/path/to/athenak
exe=/path/to/athena
work="$PWD/run"
mkdir -p "$work/base" "$work/amr"

(cd "$work/base" && "$exe" \
  -i "$src/inputs/cgl_lf_paper/cgl_lf_paper_amr_192x192x256_base.athinput")

(cd "$work/amr" && "$exe" \
  -r "$work/base/rst/cgl_amr192_base.00004.rst" \
  -i "$src/inputs/cgl_lf_paper/cgl_lf_paper_amr_192x192x256_current_l2.athinput")
```

The second file is an overlay: mesh geometry, MHD state, turbulence state, and
parameters not explicitly replaced are loaded from the restart. Preserve the
rank layout when restarting a rank-local AMR checkpoint. A shared checkpoint
is required when changing rank count.

For long queue-limited runs, use AthenaK's internal wall-clock limit to write a
clean terminal checkpoint, then resume from the complete rank-zero path under
`rst/rank_00000000/`; AthenaK maps that path to the corresponding rank-local
file on every rank.

## Validation Evidence

### Focused regression coverage

The checked-in CGL suites exercise:

- conservative primary-state refinement for density, all momentum components,
  total energy, and passive scalar density;
- smooth primitive prolongation, strong anisotropy, low magnetic field, mirror
  and firehose constraints, and separately counted repair fallbacks;
- static coarse/fine interfaces and comparisons with uniform references;
- current-triggered refinement followed by derefinement;
- repeated 3-D LF+STS refinement churn, including a 27-to-216-to-27 block
  transition;
- restart immediately before or after a regrid, including repair-counter
  continuity;
- post-CT coarse-state encoding and face-field ghost refresh;
- repeated derived `mhd_divb` output while the local block count changes;
- MPI reductions, MPI+GPU restart, and multi-node 3-D LF AMR execution.

Focused conservation tests require primary conserved histories to remain
within an absolute `5e-12` of their initial values. Smooth tests require zero
unexpected AMR primary repairs and clean fatal LF counters.

The main regression entry points are:

- `tst/test_suite/cgl/test_cgl_amr_gpu.py`
- `tst/test_suite/cgl/test_cgl_amr_mpi_gpu.py`
- `tst/test_suite/cgl/test_cgl_landau_fluid_cpu.py`

The Frontier qualification used Cray CCE 20.0.0, Cray MPICH 9.0.1, ROCm
6.4.2, and a Release HIP build. It passed the CPU/MPI, single-GPU, and two-node
16-rank MPI+GPU phases.

### Full-scale 3-D result

The production rehearsal used the supplied decks and completed at `t=2`,
cycle 9003, on 64 Frontier GPUs. The final topology contained 4026 leaf
MeshBlocks: 10 at level 1 and 4016 at level 2, distributed at 62--63 blocks per
rank. Continuation logs after the `t~=1.5` checkpoint alone recorded 1547 block
creations and 413 deletions.

The run passed the following gates:

- no density, energy, nonfinite, or empty-interval AMR primary repairs;
- zero fatal LF density-floor, pressure-floor, nonfinite, nonpositive, and
  hard-bound counters;
- complete 64-rank restart and full-state output groups;
- populated level-1 and level-2 leaf topology within the per-rank cap;
- successful repeated restart through evolving topology;
- successful terminal `mhd_divb` and full-state derived output after the local
  block count changed.

The fixed-grid base conserved mass to `3.3e-16` relative. Across the complete
AMR interval from `t=1` to `t=2`, the history mass changed by
`-3.214915000882e-10` relative. An independent long-double integration of the
rank-local density output from `t=1.5001283551` to `t=2` measured
`-2.439244823569e-10`; the binary and in-code totals agreed within `5.4e-12`.
No primary-field repair accounted for the change. This is consistent with
finite-precision accumulation through the full-scale AMR workload, not the
nonconservative primitive-reconstruction defect guarded by the focused tests.

For this scale, the accepted terminal mass gate is `1e-9` relative. This does
not replace the stricter `5e-12` focused-regression tolerance. Conservation
tolerances should scale with workload while remaining far below changes that
would indicate a nonconservative transfer path.

## Recommended Run-Time Gates

For a new CGL AMR configuration, require all of the following before using it
for production science:

1. The requested terminal time is reached and all rank-local output groups are
   complete.
2. `lf_dfloor`, `lf_pfloor`, `lf_nonfin`, `lf_nonpos`, and `lf_hardbd` remain
   zero in strict validation.
3. Density, energy, nonfinite, and empty-interval AMR primary-repair counters
   remain zero unless the run is explicitly designed as a repair stress test.
4. Mass conservation meets a documented tolerance appropriate to the problem
   size and number of regrids.
5. Every intended refinement level is populated and no rank exceeds
   `max_nmb_per_rank`.
6. Face-field divergence remains finite and consistent with the corresponding
   uniform or static-reference problem.
7. At least one restart crosses a period of active refinement or derefinement.

Mirror, firehose, anisotropy, and hard-wall projection counters are not
automatically failures in limiter-active physics. Interpret them according to
the configured closure. Primary repairs and fatal LF counters are the stronger
numerical safety indicators.

## Validated Envelope And Remaining Gaps

Confidence is high for Cartesian, periodic CGL and CGL+LF+STS runs using PLM,
HLLE, primitive AMR prolongation, up to two additional levels, and the tested
MPI+GPU restart/output workflow. Both nearly uniform fine coverage and
aggressive refinement/derefinement have been exercised.

The following configurations require separate validation:

- custom physical or user boundaries during LF stages;
- non-Cartesian coordinates;
- SR, GR, or dynamical-GR MHD, where CGL is not currently supported;
- substantially deeper AMR hierarchies or different MeshBlock shapes;
- reconstruction, Riemann solver, or integrator combinations outside the
  tested PLM/HLLE/RK2/RKL2 path;
- extreme near-floor or near-zero-field states intended to trigger primary
  repairs;
- CGL LF pressure-work recording with primitive AMR prolongation.

These are coverage boundaries, not known failures. Qualify them with focused
conservation, topology, restart, and representation-safety tests before use.
