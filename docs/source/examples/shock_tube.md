# Example: Shock Tube

The Sod hydrodynamic shock tube is the public first-run example. It uses the
built-in `shock_tube` problem generator selected in
`inputs/hydro/sod.athinput`.

## Run Sod

```bash
cmake -S . -B build
cmake --build build
./build/src/athena -i inputs/hydro/sod.athinput -d run-sod
```

A verified public run produced:

```text
run-sod/Sod.hydro.hst
run-sod/tab/Sod.hydro_w.00000.tab
...
run-sod/tab/Sod.hydro_w.00025.tab
```

The table header identifies columns including `x1v`, `dens`, `velx`, `vely`,
`velz`, and `eint`. Plot density from the final table with:

```python
import matplotlib.pyplot as plt
from vis.python import athena_read

data = athena_read.tab("run-sod/tab/Sod.hydro_w.00025.tab")
plt.plot(data["x1v"], data["dens"])
plt.xlabel("x1")
plt.ylabel("density")
plt.show()
```

## Related Built-In Shock Input

The Brio-Wu MHD shock tube uses the same built-in dispatcher:

```bash
./build/src/athena -i inputs/mhd/bw.athinput -d run-bw
```

The shipped Brio-Wu deck writes primitive and cell-centered magnetic table
outputs plus history data.

## Configuration Notes

The Sod deck uses `nghost = 2`, `reconstruct = plm`, and `rsolver = llf`.
Switching to `ppm4`, `ppmx`, or `wenoz` requires `nghost >= 3`; combining
higher-order reconstruction with `fofc = true` requires `nghost >= 4`.

For a small resolution experiment, edit a copy of the input deck. Runtime
overrides update existing parameters only:

```bash
cp inputs/hydro/sod.athinput my_sod.athinput
./build/src/athena -i my_sod.athinput -d run-sod-128 mesh/nx1=128 meshblock/nx1=128
```

## See Also

- [5-Minute Quickstart](../quickstart.md)
- [Configuration](../configuration.md)
- [Problem Generators](../modules/pgen.md)
