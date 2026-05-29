# Visualization Utilities (`vis/python`)

The public AthenaK tree ships lightweight Python readers and plotting helpers
under `vis/python/`. Use the output type written by your input deck as the
starting point; the default Sod deck writes tabular and history data, while
binary conversion requires an output stream with `file_type = bin`.

## Shipped Utilities

| Script | Public behavior |
| --- | --- |
| `athena_read.py` | Reads `.tab`, `.hst`, and `.athdf` data into Python objects. |
| `plot_tab.py` | Displays a selected variable from one or more formatted `.tab` dumps. |
| `plot_hst.py` | Displays a selected history variable against simulation time. |
| `bin_convert.py` | Reads Athena binary mesh output and can write `.athdf` plus `.xdmf` companions. |
| `make_athdf.py` | Batch wrapper around `bin_convert.py` for a filename stem. |
| `plot_mesh.py`, `plot_slice.py` | Additional plotting tools for mesh/slice analysis; consult their command-line help for required arguments. |
| `cartgrid.py` | Reader/helper for Cartesian-grid binary output. |

`bin_convert_new.py`, `read_pdf.py`, and `read_sphslice.py` are not present in
the public baseline and are therefore not public interfaces on this site.

## Plot The Quickstart Output

After running the Sod command in [Quickstart](../quickstart.md), inspect the
table and history outputs with:

```bash
python vis/python/plot_tab.py \
  -i run-sod/tab/Sod.hydro_w.00000.tab -v dens

python vis/python/plot_hst.py \
  -i run-sod/Sod.hydro.hst -v mass
```

These plotting scripts require Matplotlib in addition to NumPy. For scripted
analysis without plotting:

```python
from vis.python import athena_read

table = athena_read.tab("run-sod/tab/Sod.hydro_w.00000.tab")
history = athena_read.hst("run-sod/Sod.hydro.hst")
print(table["dens"][:5])
print(history["mass"][-1])
```

## Convert Binary Mesh Output

Add a binary output stream to a run-specific hydrodynamic input deck:

```ini
<output3>
file_type = bin
variable  = hydro_w
dt        = 0.1
```

Run that deck into its own directory, then convert one dump to HDF5/XDMF:

```bash
./build/src/athena -i my_sod_binary.athinput -d run-bin
python vis/python/bin_convert.py run-bin/bin/Sod.hydro_w.00000.bin
```

The script writes `run-bin/bin/Sod.hydro_w.00000.athdf` and the corresponding
`.athdf.xdmf` file. It requires NumPy and `h5py`; `tqdm` is optional for
progress display.

For a set of binary dumps sharing a stem:

```bash
python vis/python/make_athdf.py run-bin/bin/Sod.hydro_w. -v
```

For programmatic use, `bin_convert.py` exposes `read_binary`,
`read_coarsened_binary`, `read_all_ranks_binary`,
`read_all_ranks_coarsened_binary`, the `*_as_athdf` assembly functions,
`write_athdf`, `write_xdmf_for`, and `convert_file`.

## Output Guidance

- Use `vtk` output directly in ParaView or VisIt when conversion is unnecessary.
- Use `tab` and `hst` for low-overhead plots and diagnostics.
- Use `bin` plus `bin_convert.py` when you need HDF5/XDMF products or
  programmatic mesh-block access.

See [Outputs](../modules/outputs.md) for output-stream configuration and
[Running Simulations](../running.md) for filename conventions.
