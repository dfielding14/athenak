## Convert And Read IO Outputs

Add a binary output stream to a run-specific input deck:

```ini
<output3>
file_type = bin
variable = hydro_w
dt = 0.1
```

Run that deck into its own directory, then convert one dump to HDF5/XDMF:

```bash
./build/src/athena -i my_sod_binary.athinput -d run-bin
python vis/python/bin_convert.py run-bin/bin/Sod.hydro_w.00000.bin
```

For programmatic use:

```python
from vis.python import bin_convert

raw = bin_convert.read_binary("run/bin/simulation.prim.00000.bin")
node = bin_convert.read_binary(
    "run/bin/node_00000000/simulation.prim.00000.bin",
    assemble_shards=True,
)
```

`bin_convert.py` reads `.bin` and supported `.cbin` files and retains the
conversion helpers used by existing workflows. Passing `assemble_shards=True`
assembles canonical `rank_########` or `node_########` sibling families.

Modern PDF products use the `AKPDFV2` representation and are read with
`read_pdf.py`. Fixed-radius binary angular slices are read with
`read_sphslice.py`. Both readers automatically assemble a sharded family when
given a shard path:

```python
from vis.python import read_pdf, read_sphslice
from vis.python.io_reader_common import ReaderLimits

limits = ReaderLimits(max_live_bytes=512 * 1024 * 1024)
pdf_path = "run/pdf_rho/node_00000000/simulation.00000.pdf"
pdf_header_path = "run/pdf_rho/node_00000000/simulation.header.pdf"
pdf_header = read_pdf.read_pdf_header(pdf_header_path, limits=limits)
pdf = read_pdf.read_pdf(pdf_path, limits=limits)
surface_path = (
    "run/bin/node_00000000/"
    "simulation.shell.r_5.0000000000000000e-01.00000.sph.bin"
)
surface_header = read_sphslice.read_sphslice_header(surface_path, limits=limits)
surface = read_sphslice.read_sphslice(
    surface_path, limits=limits
)
```

`ReaderLimits` also exposes `max_header_read_bytes` and
`max_payload_read_bytes`. The header-only APIs validate one header's intrinsic
metadata and declared shard identity without materializing the full numerical
payload. Full readers discover siblings and reject incomplete or inconsistent
families. Use the header-only APIs for bounded declaration inspection before
large analysis jobs.

For runnable commands, see
[IO Outputs And Sharding](../examples/io_outputs_and_sharding.md).
