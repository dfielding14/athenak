# Q-016 Bounded Particle Provenance And Spectra

This note records a bounded host-regression path for particle provenance and
weighted spectra. It does not qualify a production shock campaign.

## Persistent Particle Metadata

Every cosmic-ray particle carries:

- integer `PGID`: current owning MeshBlock ID;
- integer `PTAG`: persistent tracking ID;
- integer `PSP`: species ID;
- integer `PCRSOURCE`: persistent source cohort (`0=initial`,
  `1=shock_injected`);
- real `IPWT`: macro-particle weight;
- real `IPT_BIRTH`: creation time;
- real `IPF0`: sampled initial analytic delta-f background value;
- real `IPDFWT`: evolving delta-f perturbation weight.

The initial population sets `PCRSOURCE=initial` and records the mesh start time.
The `pic_parallel_shock` injection source sets `PCRSOURCE=shock_injected` and
records the injection time. Generic migration and restart loops copy all real
and integer particle slots. Because the persistent record layout changes,
particle restart schema version `7` must be used by the writer and both loader
paths.

## Exported Fields

Particle VTK output exposes named scalar fields:

```text
gid ptag species cr_source macro_weight birth_time deltaf_f0 deltaf_weight
```

The `vel` vector contains physical velocity. The tracked-particle stream keeps
its existing six-float row payload and labels the row index as `ptag` with
columns `x,y,z,vx,vy,vz`.

## Bounded Spectrum Helper

`tst/publication/q016_particle_spectra.py` reads one particle VTK file and emits
a JSON physical-speed spectrum resolved by species, source, and bounded
birth-time cohort. Both physical-speed and birth-time bin edges are explicit.
The emitted density normalization is:

```text
weighted_sum_in_bins / diff(bin_edges)
```

The helper has two explicit weight semantics:

- `full_f`: histogram weight is `macro_weight`.
- `delta_f_perturbation`: histogram weight is
  `macro_weight * deltaf_weight`. This is the signed sampled perturbation only;
  the analytic background represented by `deltaf_f0` is intentionally not
  added.

## Regression Boundary

`tst/publication/test_pic_q016_particle_provenance_spectra.py` reconstructs
typed integer PVTK decoding and both spectrum modes independently from a
synthetic binary fixture.

`tst/scripts/particles/pic_q016_particle_provenance.py` runs a bounded serial
parallel-shock fixture, requires restart schema `7`, compares uninterrupted and
restarted particle metadata, observes MeshBlock migration, and independently
reconstructs exported spectra.

The particle VTK reader requires exactly one `POINTS` section and one matching
`POINT_DATA` section, rejects duplicate scalar or vector sections, requires
float vector payloads, and rejects trailing unknown content. The harness merges
schema-compatible same-cycle optional `gid` slices before tag-sorted
comparison. This synthetic `gid`-slice coverage is not MPI decomposition
evidence.

The serial harness invokes `./athena` directly by default. A controlled MPI
execution sets `ATHENA_Q016_NPROC` and may set `ATHENA_Q016_LAUNCHER` or
`MPIEXEC`; the harness rejects multi-rank requests when the selected executable
reports MPI support disabled.

The remaining qualification gaps are MPI decomposition coverage, repeated AMR
refinement/derefinement lifetime coverage, HIP parity, and a preregistered
qualifying shock campaign.
