# Reproducing the 1024^2 scalar-mixing stream-function run

This note records the exact setup used for the 2D stream-function scalar-mixing
production run with a left/right scalar step.

## Build

Configure and build AthenaK with MPI and the `scalar_mixing` problem generator:

```bash
cmake -S . -B build-mpi -DPROBLEM=scalar_mixing -DAthena_ENABLE_MPI=ON
cmake --build build-mpi -j8
```

## Run

Launch the run on 4 MPI ranks:

```bash
mkdir -p build-mpi/runs/scalar_mix_stream_1024_step
cd build-mpi/runs/scalar_mix_stream_1024_step
mpiexec -n 4 ../../src/athena \
  -i /path/to/repo/inputs/hydro/scalar_mixing_stream_1024_step.athinput
```

The committed input file uses:

- `1024 x 1024` cells on `[-0.5, 0.5]^2`
- `2 x 2` meshblocks of size `512 x 512`
- `turb_use_stream_function = true`
- target velocity slope `-5/3`
- `turb_nlow = 2`, `turb_nhigh = 512`
- left/right scalar step: `theta = 0` for `x < 0`, `theta = 1` for `x >= 0`
- `scalar_diffusivity = 1e-3`
- `tlim = 10`
- binary `hydro_w` outputs every `dt = 1`

For practicality on a laptop-scale local MPI run, the committed production input
uses stream-mode importance sampling with:

```text
turb_k_crit = 16
```

This keeps the requested target spectrum while making the direct mode-sum
initialization tractable on local hardware.

## Post-processing

Generate the scalar-spectrum overlay and scalar-map figures:

```bash
python3 vis/python/scalar_mixing_stream_run_analysis.py
```

By default the script reads:

```text
build-mpi/runs/scalar_mix_stream_1024_step/bin
```

and writes figures to:

```text
build-mpi/runs/scalar_mix_stream_1024_step/analysis
```
