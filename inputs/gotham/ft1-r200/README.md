# GOTHAM: $f_T=1$ with $R_{\rm sonic}=R_{200}$

This run starts from the selected cooling-flow CGM with **no initial gas disk**.
The fixed stellar-disk gravitational potential is retained.

| Parameter | Value |
| --- | --- |
| Sonic radius $R_{\rm sonic}=R_{200}$ | 144.76736315360156 kpc |
| Matching radius $R_b=1.5R_{200}$ | 217.15104473040236 kpc |
| Temperature condition | $T(R_b)=\mu m_p v_c^2(R_b)/(\gamma k_B)$ |
| Boundary temperature | 147007.741674 K |
| Inflow rate | $4.11368408045\,M_\odot\,{\rm yr}^{-1}$ |
| CGM mass inside $R_{200}$, before mesh sampling | $6.83863\times10^9\,M_\odot$ |
| Profile coverage | 0.1 to 354.52914961347017 kpc, including corner ghost cells |
| Cluster seed / count | 42 / 1500 |

`athinput.gotham` retains the fiducial potential, production mesh and AMR settings,
CGM rotation (`r_circ=10` kpc), initial turbulent motions, ongoing turbulence
driving, local heating, cooling, dust settings, and star clusters from
`../athinput.gotham`. These additions and nonspherical gravity will cause the
spherical IC to relax. Only the gas disk is omitted, through
`problem/initialize_disk=false`; no disk table is needed.

The run ends at `time/tlim=10` (1 Gyr), with `SN/delay=10` preserved. This is the
relaxation phase and ends before the first supernova. A later feedback phase
requires restarting with a larger time limit and the matching particle restart;
keep the supernova delay at 10.

## Run

Build AthenaK from this checkout with `-DPROBLEM=gotham` and the MPI/GPU settings
appropriate to the target machine. The cooling tables are compiled into GOTHAM.
Copy this entire directory to the run location, then launch the executable from
that directory so the relative table paths resolve:

```bash
cp -R inputs/gotham/ft1-r200 /path/to/new-run
cd /path/to/new-run
/path/to/athena -i athinput.gotham
```

Use the machine's MPI launcher around the final command for a parallel run.
The supplied production mesh is intended for a parallel simulation.

- `tables/cgm_profiles.txt`: radius, $P/\rho$, density, and positive inward speed,
  all in code units. GOTHAM supplies the inward sign.
- `tables/clusters.txt`: positions, velocities, and formation times in code
  units; cluster masses in solar masses.
- `ic_parameters.json`: generation parameters. Its absolute `source_input`
  records the original local source; the portable source is `../athinput.gotham`.

## Regenerate

From the repository root, initialize the pinned submodules and install
`external/CGM_ICs/requirements.txt` into a Python environment. Generate into a
new directory; the exporter refuses to overwrite an existing one:

```bash
git submodule update --init --recursive
python external/CGM_ICs/generate_gotham_ics.py \
    inputs/gotham/athinput.gotham build/gotham-ft1-r200 \
    --r-sonic 144.76736315360156 --r-boundary 217.15104473040236 \
    --f-t 1 --no-disk --seed 42
```

## Validation

The 4644 exported CGM rows reproduce the selected comparison profile to
floating-point precision. The mass flux is constant to $10^{-12}$ relative
tolerance, and the boundary temperature, units, coverage, and retained template
parameters were checked. All 15 exporter tests pass, including no-disk generation
from both the CLI option and input setting. The compiled C++/Python gravity and
cooling comparison also passes.

A fresh serial CPU GOTHAM build loaded these tables and clusters without a disk
file and advanced two steps. Initial density matched the spherical CGM profile
to $5.36\times10^{-8}$ relative precision in the binary output. Initial and evolved
fields remained finite with positive density and internal energy. A separate
default-disk case also passed. These startup tests used a reduced mesh, no AMR,
and fewer turbulent modes; production MPI/GPU execution and the 1 Gyr relaxation
have not been run.
