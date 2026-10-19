# GOTHAM full-cube PDF reconstruction

This is a deliberately one-off MPI+Kokkos reducer for the five GOTHAM
simulations. It repairs the corrupted AthenaK N-D PDF streams from saved
full-volume `hydro_w` AMR shards without constructing a dense Cartesian cube.

The executable accepts only the GOTHAM schema:

```text
dens velx vely velz eint s_00 s_01 s_02
```

Geometry is read as float64 and fields as float32. MeshBlock dimensions are
read from each snapshot's embedded input deck; production snapshots use
`64^3` cells per block.

## What it computes

- `--products original`: exact-compatible replacements for `output23` through
  `output36`, including AthenaK underflow/overflow and C-order flattening.
- `--products science`: 69 compact science products covering volume-filling
  phases, signed opening angle, gross inflow/outflow, kinetic/enthalpy
  transport, Mach structure, vertical breakout, angular momentum, metals,
  dust, and signed net CGM cooling diagnostics, including cooling-weighted
  phase and thermokinematic PDFs.
- `--products all`: both sets in one streaming pass.

The original set uses 0.292 GiB of histogram storage per rank. The complete
83-product set uses 1.696 GiB.

The 69-science/83-complete product implementation was frozen and passed the
full Frontier release sequence on 2026-06-09 under identity
`04bbaa1e76bd25becff45d2b6e2975a64e477fdd71905def1eb70fd8ad260b9b`.
A science-product campaign then completed all 161 inventory snapshots. The
review changes made after that campaign correct the cooling-time regularizer
and make release-source verification fail closed, so this source tree requires
a new identity and release sequence before any future production run.

The cooling products use the same signed net source convention as AthenaK's
saved `cooling_time` field: positive values are cooling-dominated and negative
values are heating-dominated. `science_r_theta_edot_cool_volume` is the local
net source-rate density in `erg s^-1 cm^-3`, not a cell-integrated luminosity,
so it remains comparable across the 4 pc and 8 pc meshes. The
cooling-weighted temperature products instead use cell-integrated signed net
cooling luminosity in `erg s^-1`.

## Execution model

- One MPI rank per accelerator.
- Independent shard assignment by sorted shard index.
- POSIX `pread` of complete fixed-size MeshBlock records.
- Raw chunk copy to the device.
- One Kokkos derive-and-bin kernel per chunk.
- Product-by-product root-only MPI reductions.
- Dense output compatible with
  `/lustre/orion/ast207/proj-shared/gotham/analysis/gotham_analysis/readers/pdf.py`.
- Atomic `.partial` publication and a final reconstruction manifest.

Full-snapshot runs reject mismatched complete header digests, malformed block
records, invalid hydro state, duplicate AMR leaf keys, ancestor/descendant
overlap, non-contiguous expected shard manifests, and domain-volume closure
errors.

## Build

CPU on Andes:

```bash
cmake -S tools/gotham_pdf_rebuild -B build-gotham-pdf-cpu \
  -DCMAKE_BUILD_TYPE=Release -DKokkos_ENABLE_SERIAL=ON
cmake --build build-gotham-pdf-cpu --parallel 16
```

Andes K80 CUDA uses normal GCC9 OpenMPI because reductions stage through host
memory. Do not use `openmpi/3.1.6-gpu`; it injects incompatible RPATH and CUDA
runtime dependencies. Bundled Kokkos 4.7.02 also needs the one-line Kepler
compile fix in `patches/kokkos-4.7.02-kepler-compile.patch`.

```bash
module restore
module load cuda/11.2.2
(cd kokkos && patch -p1 --forward < \
  ../tools/gotham_pdf_rebuild/patches/kokkos-4.7.02-kepler-compile.patch)
cmake -S tools/gotham_pdf_rebuild -B build-gotham-pdf-andes-cuda \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER="$PWD/kokkos/bin/nvcc_wrapper" \
  -DMPI_CXX_COMPILER="$(command -v mpicxx)" \
  -DKokkos_ENABLE_CUDA=ON -DKokkos_ENABLE_CUDA_LAMBDA=ON \
  -DKokkos_ARCH_KEPLER37=ON
cmake --build build-gotham-pdf-andes-cuda --parallel 16
```

Frontier HIP:

```bash
env -i HOME="$HOME" USER="$USER" LOGNAME="$USER" SHELL=/bin/bash \
  PATH=/usr/local/bin:/usr/bin:/bin:/opt/cray/pe/lmod/lmod/libexec \
  bash --noprofile --norc <<'FRONTIER_BUILD'
set -euo pipefail
source /opt/cray/pe/lmod/lmod/init/bash
module use \
  /opt/cray/pe/lmod/modulefiles/core \
  /opt/cray/pe/modulefiles/Core \
  /opt/cray/pe/lmod/lmod/modulefiles/Core \
  /opt/cray/pe/lmod/modulefiles/craype-targets/default \
  /sw/frontier/spack-envs/modules/Core/25.03
module load cpe/25.09
module load craype-x86-trento craype-network-ofi craype-accel-amd-gfx90a
module load PrgEnv-cray cray-mpich/9.0.1 rocm/6.4.2 cce/20.0.0 cmake/3.30.5

cmake -S tools/gotham_pdf_rebuild -B build-gotham-pdf-frontier \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=CC -DMPI_CXX_COMPILER=CC \
  -DKokkos_ENABLE_HIP=ON -DKokkos_ARCH_ZEN3=ON \
  -DKokkos_ARCH_VEGA90A=ON \
  -DCMAKE_CXX_FLAGS="-I${ROCM_PATH}/include -munsafe-fp-atomics" \
  -DCMAKE_EXE_LINKER_FLAGS="-L${ROCM_PATH}/lib -lamdhip64"
cmake --build build-gotham-pdf-frontier --parallel 16
FRONTIER_BUILD
```

The sanitized shell is deliberate. After the July 2026 Frontier module
transition, `module restore` can retain CCE 18 include variables and mix them
with CCE 20. The recipe above was compile-tested with CCE 20.0.0, ROCm 6.4.2,
and Cray MPICH 9.0.1.

## Run

Example mapped repair:

```bash
srun ... gotham_pdf_rebuild \
  --input-dir /path/to/res_8pc/phase2/bin \
  --sequence 00028 \
  --output-number 00062 \
  --output-time 10.520002963096713 \
  --expected-shards 1024 \
  --output-dir /separate/rebuild/tree/res_8pc/phase2/00028 \
  --products original \
  --chunk-blocks 32
```

Never write into the archived simulation PDF directories. The static Frontier
manifests under `jobs/manifests/` contain exact output mappings for all 161
available full cubes: 140 map to archived PDFs and 21 uniform-phase cubes have
no archived PDF control.

## Validation

```bash
module reset
module load gcc/9.3.0 python/.3.11-anaconda3
source /lustre/orion/ast207/proj-shared/gotham/venv_gotham/bin/activate
GOTHAM_PDF_REBUILD_EXE="$PWD/build-gotham-pdf-cpu/gotham_pdf_rebuild" \
  python -m pytest -q tools/gotham_pdf_rebuild/tests
```

The independent Python oracle checks all original streams, bin boundaries,
derived variables, weights, MPI invariance, `32^3` and `64^3` layouts, and the
production PDF reader.

After the complete Frontier `res_8pc/phase2/00028` repair:

```bash
python tools/gotham_pdf_rebuild/validate_real_8pc.py \
  /path/to/rebuilt/00028 \
  /lustre/orion/ast207/proj-shared/brent/gotham/res_8pc/phase2
```

The latest pre-review Frontier release sequence completed on 2026-06-09 using
immutable identity
`04bbaa1e76bd25becff45d2b6e2975a64e477fdd71905def1eb70fd8ad260b9b`.
The original, science, and all-product 1024-shard goldens, the 256-node
2048-shard canary, and the 512-shard no-control canary all passed their full
validators and predeclared performance budgets. Those reports validate only
their frozen source and executable; they do not authorize the reviewed source
tree. See `FRONTIER_HANDOFF.md` for paths, job IDs, metrics, campaign status,
and the operator procedure for a new freeze.

## Plot every rebuilt PDF

`plot_rebuilt_pdfs.py` follows the slice-plot conventions in the shared GOTHAM
analysis tree. It makes informative pairwise 2D marginals, summing over the
remaining axes. It does not plot raw weighted volume versus radius, polar
angle, or a physical quantity because those profiles mostly show geometry or
normalization rather than physical state. Mass- and volume-weighted
radius--polar-angle heatmaps are also omitted because they reproduce weighted
geometry rather than the physical quantity carried by the product;
transport-weighted radial-angle heatmaps are retained.

Products containing radius, polar angle, and another physical variable make
three shared summaries:

- a 2x2 quantity--radius figure for all angles, the polar 30 degrees, the
  midplane 30 degrees, and the intermediate band; and
- a 3x2 quantity--polar-angle figure at the nearest PDF radial bins to 2, 4,
  8, 16, 32, and 64 kpc; and
- a 2x1 line figure showing conditional mean and median quantity versus
  radius for the four angular selections, then versus polar angle for the six
  radial shells.

The renderer uses physical labels, configured colormaps, four-sided inward
ticks, matched colorbars, phase/Mach reference lines, and fixed-width titles
showing time since the first SN. Product IDs and under/overflow warnings are
not drawn on figures; they remain available in metadata.

These panels are bin-integrated weighted histograms, not
probability-density-normalized distributions. As in the shared PDF preview
renderer, every axis is restricted to interior bins before marginalization;
metadata records underflow/overflow by axis. A recorded visual-support mask
may omit at most `1e-4` of a panel's absolute weight so empty and negligible
bins do not look like signal. Volume-weighted radial panels are divided by the
selected angular sector's analytic shell volume. Fixed-shell angular panels
are divided by each radial/angular wedge volume. Their colors therefore show
conditional volume fractions instead of raw geometric volume. Angular sector
boundaries that fall inside an existing PDF bin use fractional bin overlap. A
plot from a partial or geometrically unvalidated reducer manifest is visibly
watermarked. Conditional mean/median lines use the underlying histogram
weights directly; any common normalization cancels within each selected
radius or angle bin. Signed-weight products omit median profiles because a
signed median is undefined.

Radial transport panels retaining `coord_r` are divided by radial-bin width
and shown as shell-crossing diagnostics. Transport panels without radius
remain raw volume-integrated moments with an extra length dimension. Energy
transport includes kinetic plus enthalpy advection, not gravitational,
magnetic, or radiative terms.

```bash
module reset
module load gcc/9.3.0 python/.3.11-anaconda3
source /lustre/orion/ast207/proj-shared/gotham/venv_gotham/bin/activate

python tools/gotham_pdf_rebuild/plot_rebuilt_pdfs.py \
  /path/to/rebuilt/00028 \
  --output-dir /path/to/rebuilt_plots \
  --render-profile production \
  --workers 8 \
  --skip-existing
```

Different snapshots can append to the same shallow plot root. PNGs are under
`png/` and provenance is under `metadata/`. PNG filenames follow the short
semantic form `scientific_view.simulation.t########.###Myr.png`; the
zero-padded time is absolute simulation time so lexical glob order remains
chronological across pre- and post-SN frames. Figure titles show time relative
to the first SN. Specific transport weights and marginalized axes are included
only when needed to distinguish scientifically different views; signed
`cos(theta)` is distinct from folded `|cos(theta)|`. The renderer refuses
semantic filename collisions instead of silently overwriting a plot.
Vertical-geometry views use cylindrical radius on the x axis and absolute
height on the y axis, with filenames beginning
`cylindrical_radius_absolute_z`.
Uniform-phase branches use distinct descriptive simulation names such as
`res_8pc_uniform`, preventing same-time branch collisions and producing
separate movies.

For example, the volume-weighted pressure products produce these movie-ready
sequences:

```bash
png/pressure_radius_costheta.res_4pc_highmdot.t*.png
png/pressure_radius_theta_cuts.res_4pc_highmdot.t*.png
png/pressure_costheta_radius_shells.res_4pc_highmdot.t*.png
```

`make_rebuilt_pdf_movies.py` groups every descriptive view by simulation and
encodes H.264 High/YUV420p MP4s with a front-loaded index for QuickTime. The
completed production suite contains 26,565 PNGs and 1,155 movies:

```bash
python tools/gotham_pdf_rebuild/make_rebuilt_pdf_movies.py \
  /path/to/rebuilt_plots/png \
  --output-dir /path/to/rebuilt_plots/movies \
  --ffmpeg /path/to/ffmpeg \
  --ffprobe /path/to/ffprobe \
  --max-dimension 1440 \
  --workers 14
```

Use repeatable `--product` shell patterns for a subset. `--overwrite` replaces
only the matching snapshot/product files and leaves unrelated snapshots
alone; `--skip-existing` reuses plots only when source, plotter, shared
analysis helpers, style, library versions, and rendering options match. The
renderer requires the reducer's final `rebuild_manifest.json`, rejects
non-finite weights and negative values in positive-weight products, and locks
the output tree against concurrent writers.

See `FRONTIER_HANDOFF.md` before production.
