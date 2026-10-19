# Synthetic/reference validation

This directory validates the standalone reducer against an independent Python
implementation. The synthetic writer emits exact AthenaK binary output v1.1
GOTHAM `hydro_w` shards: double geometry, float32 fields, eight variables, and
runtime-configurable `nx1/nx2/nx3`.

The tests cover:

- C-order dense PDF layout, dimension strides, and emitted bin edges
- regular, underflow, and overflow bins
- all stored and derived variables used by the reducer
- volume, mass, signed flux, inflow/outflow, kinetic, and thermal weights
- invariant global results with one and two MPI ranks
- exhaustive `32^3` fixtures plus an actual processed `64^3` production-layout block
- all 14 original broken streams and the production analysis PDF reader

The quick suite writes the compact science product set. The `large` test also
writes the original product set, which is approximately 299 MiB per run.

On Andes:

```bash
module load python/3.7-anaconda3
cd /ccs/home/dfielding/athenak-gotham-pdf-rebuild/tools/gotham_pdf_rebuild/tests
pytest -v -m "not large"
pytest -v -m large
```

Override paths when testing another build or analysis checkout:

```bash
export GOTHAM_PDF_REBUILD_EXE=/path/to/gotham_pdf_rebuild
export GOTHAM_ANALYSIS_ROOT=/path/to/gotham/analysis
export GOTHAM_MPIEXEC=/path/to/mpiexec
```

Generate standalone synthetic shards with:

```bash
python synthetic_gotham.py /tmp/gotham-synthetic --blocks 4 --shards 2 --nx1 32
python synthetic_gotham.py /tmp/gotham-production-layout --blocks 1 --shards 1 --nx1 64
```
