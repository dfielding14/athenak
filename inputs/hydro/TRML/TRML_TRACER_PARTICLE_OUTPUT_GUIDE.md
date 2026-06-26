# TRML tracer-particle output guide

This guide covers the tracer outputs produced by
`TRML_with_Tracers_and_Tracking.athinput`. The thermodynamic history (`.thp`) is
the primary analysis product. Particle VTK files are optional visualization
snapshots, while restart files preserve simulation state and are not intended
for analysis.

## Quick start

For the canonical basename and output ID, the main history file is

```text
<run-directory>/prtcl_thermo_history/
  TRML_with_Tracers_and_Tracking.thermo.thp
```

Inspect it or convert it to NumPy's `.npz` format with the repository reader:

```bash
REPO=/path/to/athenak-DF
RUN=/path/to/run-directory
THP="${RUN}/prtcl_thermo_history/TRML_with_Tracers_and_Tracking.thermo.thp"

python "${REPO}/vis/python/read_prtcl_thermo_history.py" \
  "${THP}" --npz "${RUN}/tracer_thermo.npz"
```

The command prints the record count and column names. The `.thp` file is a
versioned binary file; do not attempt to read it as text or with `numpy.load`.

It can also be read directly from Python:

```python
from pathlib import Path
import sys

repo = Path("/path/to/athenak-DF")
run = Path("/path/to/run-directory")
sys.path.insert(0, str(repo / "vis" / "python"))

from read_prtcl_thermo_history import read_history

path = (
    run
    / "prtcl_thermo_history"
    / "TRML_with_Tracers_and_Tracking.thermo.thp"
)
data = read_history(path)
print(data.keys())
```

`data` is a dictionary of one-dimensional NumPy arrays in long-table form:
every array has one element per particle sample. It is not a rectangular
`particle x time` array because new tracers are introduced throughout the run.

## Canonical tracer populations

The final run contains two populations:

| `seed_id` | Sampling | Schedule | Final count |
| --- | --- | --- | ---: |
| 1 | Uniform by volume throughout the initial domain | 777 particles at `t=0` | 777 |
| 2 | Uniform by volume in the uppermost active-cell slab | 45 particles at each of 51 times, `t=0, 1.5, ..., 75` | 2,295 |

Thus a complete run should contain 3,072 unique tags. Fewer particles are
present at earlier times because most of the `seed_id=2` population has not yet
been introduced.

The two populations have different sampling designs. Analyze them separately
unless there is a deliberate reason to combine them. A raw particle count is
not automatically a volume- or mass-weighted physical integral.

## Columns in the canonical `.thp` file

Every record contains the following metadata:

| Column | Meaning |
| --- | --- |
| `time` | Simulation time of this sample. |
| `cycle` | Simulation cycle of this sample. |
| `tag` | Globally unique, persistent particle identifier. Use this to join a particle across outputs and restarts. |
| `seed_id` | Population that created the particle: 1 for the initial volume or 2 for timed upper-boundary injection. |
| `x1`, `x2`, `x3` | Particle's current grid-frame cell-center coordinates. |
| `gid` | Current MeshBlock global ID. This can change after AMR, load balancing, or restart and is not a particle identifier. |

The canonical output then samples these fluid quantities at the cell containing
each particle:

| Column | Definition and interpretation |
| --- | --- |
| `density` | Cell density, `rho`. |
| `pressure` | Ideal-gas pressure, `(gamma - 1) eint`. |
| `temperature` | Code temperature, `pressure/density`. |
| `entropy` | Specific-entropy proxy, `log(pressure/density**gamma)`. |
| `v1`, `v2`, `v3` | Cell fluid velocity components in the tracked grid frame. These are not particle velocities. |
| `mach` | Grid-frame fluid speed divided by `sqrt(gamma*pressure/density)`. |
| `scalar0` | Cold-material fraction: approximately 1 for cold material and 0 for hot material. |

All quantities are in simulation code units. With the canonical values,

```text
T_cold = pres/rho_cold = 1/56.23413251903491
T_hot  = pres/rho_hot  = 1
```

The configured diagnostic phase boundaries are `T/T_cold = 10` and 30, while
radiative cooling is active only up to `T/T_cold = 20`.

## What a Monte-Carlo tracer history means

These tracers follow cell-to-cell mass flux probabilistically. Consequently:

- positions lie at cell centers and jump between neighboring cells;
- the output cadence can skip intermediate cell moves;
- finite-differencing `x1`, `x2`, or `x3` does not give a useful instantaneous
  particle velocity;
- thermodynamic columns are snapshots of the current cell's fluid state, not
  continuously carried particle properties;
- repeated samples from one tag are correlated and should not be treated as
  independent measurements.

For a time-dependent population, avoid pooling every record into one histogram:
particles introduced earlier would be counted more times. Instead, group by
simulation time for snapshots or by particle age for post-injection evolution.

## Basic integrity checks

This example reports the available records and first appearance of every tag:

```python
import numpy as np

print("records:", data["time"].size)
print("unique tags:", np.unique(data["tag"]).size)
print("time range:", data["time"].min(), data["time"].max())

tags, first_index = np.unique(data["tag"], return_index=True)
first_seed = data["seed_id"][first_index]
first_time = data["time"][first_index]
first_x3 = data["x3"][first_index]

for seed_id in (1, 2):
    selected = first_seed == seed_id
    print(
        f"seed {seed_id}: {selected.sum()} particles, "
        f"first-observation times {first_time[selected].min():g}--"
        f"{first_time[selected].max():g}, "
        f"x3 {first_x3[selected].min():g}--{first_x3[selected].max():g}"
    )

for name, values in data.items():
    if np.issubdtype(values.dtype, np.floating):
        assert np.all(np.isfinite(values)), f"non-finite values in {name}"
```

For the canonical output cadence, the first recorded time of an injected tracer
normally equals its scheduled injection time. In a run whose injection and
output cadences are not aligned, first appearance can lag injection by up to the
next history output.

Restart boundaries can occasionally produce duplicate samples with the same
`(cycle, tag)`. Check for them before calculating statistics:

```python
key = np.rec.fromarrays([data["cycle"], data["tag"]])
_, keep = np.unique(key, return_index=True)
keep = np.sort(keep)
if keep.size != data["tag"].size:
    print("duplicate (cycle, tag) rows:", data["tag"].size - keep.size)
    data = {name: values[keep] for name, values in data.items()}
```

## Plotting temperature evolution

Particle age is the most useful time coordinate for the continuously injected
population. The following produces individual temperature histories and
ensemble 16th/50th/84th percentiles for each seed population:

```python
import matplotlib.pyplot as plt
import numpy as np

T_cold = 1.0 / 56.23413251903491

# Map each record to the first recorded time of its tag.
tags, inverse = np.unique(data["tag"], return_inverse=True)
birth_time = np.full(tags.size, np.inf)
np.minimum.at(birth_time, inverse, data["time"])
age = data["time"] - birth_time[inverse]

colors = {1: "tab:blue", 2: "tab:orange"}
labels = {1: "initial volume", 2: "upper injection"}
rng = np.random.default_rng(24680)

fig, (ax_tracks, ax_stats) = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)

for seed_id in (1, 2):
    seed_tags = np.unique(data["tag"][data["seed_id"] == seed_id])
    sample = rng.choice(seed_tags, size=min(24, seed_tags.size), replace=False)
    for tag in sample:
        rows = np.flatnonzero(data["tag"] == tag)
        rows = rows[np.argsort(data["time"][rows])]
        ax_tracks.plot(
            age[rows],
            data["temperature"][rows] / T_cold,
            color=colors[seed_id],
            alpha=0.22,
            linewidth=0.8,
        )

    selected = data["seed_id"] == seed_id
    max_age = age[selected].max()
    edges = np.linspace(0.0, max_age, 61)
    centers = 0.5 * (edges[:-1] + edges[1:])
    q16 = np.full(centers.size, np.nan)
    q50 = np.full(centers.size, np.nan)
    q84 = np.full(centers.size, np.nan)
    for n, (left, right) in enumerate(zip(edges[:-1], edges[1:])):
        rows = selected & (age >= left) & (age < right)
        if np.count_nonzero(rows) >= 20:
            q16[n], q50[n], q84[n] = np.quantile(
                data["temperature"][rows] / T_cold, [0.16, 0.50, 0.84]
            )
    ax_stats.plot(centers, q50, color=colors[seed_id], label=labels[seed_id])
    ax_stats.fill_between(centers, q16, q84, color=colors[seed_id], alpha=0.2)

for axis in (ax_tracks, ax_stats):
    axis.set_yscale("log")
    axis.set_xlabel("time since first tracer observation")
    axis.axhline(1.0, color="0.4", linewidth=0.7, linestyle=":")
    axis.axhline(10.0, color="0.4", linewidth=0.7, linestyle="--")
    axis.axhline(20.0, color="0.4", linewidth=0.7, linestyle="-.")
    axis.axhline(30.0, color="0.4", linewidth=0.7, linestyle="--")

ax_tracks.set_ylabel(r"$T/T_{\rm cold}$")
ax_tracks.set_title("Example particle histories")
ax_stats.set_title("Population median and 16--84% range")
ax_stats.legend(frameon=False)
fig.tight_layout()
fig.savefig("tracer_temperature_evolution.png", dpi=180)
```

For `seed_id=1`, age is simulation time. For `seed_id=2`, it is time since the
particle first appears in the history. The latter is the appropriate coordinate
for asking how quickly newly injected hot material cools or mixes.

## Converting to the laboratory frame

The stored particle coordinates and fluid velocities are in the tracked grid
frame. Density, pressure, temperature, entropy, and `scalar0` are unchanged by a
Galilean frame transformation. For the canonical x3-only tracker,

```text
x3_lab = x3_grid + ft_dx_x3
v3_lab = v3_grid + ft_vf_x3
```

Use the frame history at matching times:

```python
from pathlib import Path
import sys
import numpy as np

repo = Path("/path/to/athenak-DF")
run = Path("/path/to/run-directory")
sys.path.insert(0, str(repo / "vis" / "python"))
import athena_read

frame = athena_read.hst(
    str(run / "TRML_with_Tracers_and_Tracking.frame_tracker.hst")
)
order = np.argsort(frame["time"])
frame_time = frame["time"][order]
frame_dx3 = frame["ft_dx_x3"][order]
frame_v3 = frame["ft_vf_x3"][order]

dx3 = np.interp(data["time"], frame_time, frame_dx3)
vf3 = np.interp(data["time"], frame_time, frame_v3)
x3_lab = data["x3"] + dx3
v3_lab = data["v3"] + vf3

gamma = 1.666666667
sound_speed = np.sqrt(gamma * data["pressure"] / data["density"])
mach_lab = np.sqrt(data["v1"]**2 + data["v2"]**2 + v3_lab**2) / sound_speed
```

Do not use the stored `mach` as a laboratory-frame Mach number: its velocity
magnitude is evaluated in the grid frame and must be recomputed as above.

## MPI, restarts, and file management

- MPI ranks gather their particle samples to rank 0, producing one `.thp` file,
  not one file per rank.
- The writer appends new time blocks to the existing file. Continue a restart in
  the same run directory with the same thermodynamic variable list.
- A changed variable list, precision, or schema cannot be appended; AthenaK
  exits with a schema-mismatch error.
- Start an unrelated rerun in a clean directory. Reusing a compatible old file
  would append a second run and create overlapping times and tags.
- Particle tags and seed schedules are preserved by restart. The canonical
  `single_file_per_rank=true` setting applies to restart files, not to `.thp`
  history output.
- The current reader loads the complete file into memory. This is modest for the
  canonical 3,072-particle run but should be considered when greatly increasing
  particle count or output frequency.

## Particle VTK snapshots

The canonical input also defines `file_type=pvtk`, but its `dt=1.0e9` effectively
disables recurring snapshots. To produce ParaView-readable particle snapshots,
override that cadence, for example:

```bash
output4/dt=5.0
```

Files then appear as

```text
<run-directory>/pvtk/
  TRML_with_Tracers_and_Tracking.tracers.00000.part.vtk
```

These files are useful for spatial visualization and contain positions plus
particle metadata such as tag, seed ID, and MeshBlock ID. Use the `.thp` stream,
not VTK, for thermodynamic trajectory analysis.

## Common mistakes

- Do not interpret `v1`, `v2`, or `v3` as particle velocities.
- Do not finite-difference the cell-center trajectory to estimate fluid velocity.
- Do not mix `seed_id=1` and `seed_id=2` without accounting for their different
  selection and injection schedules.
- Do not pool all times when estimating a snapshot distribution.
- Do not treat `gid` as a stable particle identity; use `tag`.
- Do not compare grid-frame `x3`, `v3`, or `mach` directly with lab-frame fluid
  diagnostics without applying the frame-history correction.
- Do not reuse an old `.thp` file for a fresh run.
