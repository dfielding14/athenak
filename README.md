# AthenaK

Block-based AMR framework with fluid, particle and numerical relativity solvers in Kokkos.

## Overview

AthenaK is a complete rewrite of the AMR framework and fluid solvers in the [Athena++](https://github.com/PrincetonUniversity/athena) astrophysical MHD code using the [Kokkos](https://kokkos.org/) programming model.

Using Kokkos enables *performance-portability*.  AthenaK will run on any hardware supported by Kokkos, including CPU, GPUs from various vendors, and ARM processors.

AthenaK is targeting challenging problems that require exascale resources, and as such it does not implement all of the features of Athena++.  Current code features are:
- Block-based AMR with dynamical execution via a task list
- Non-relativistic (Newtonian) hydrodynamics and MHD
- Special relativistic (SR) hydrodynamics and MHD
- General relativistic (GR) hydrodynamics and MHD in stationary spacetimes
- Relativistic radiation transport
- Lagrangian tracer particles, and charged test particles
- Numerical relativity solver using the Z4c formalism
- GR hydrodynamics and MHD in dynamical spacetimes

The numerical algorithms implemented in AthenaK are all based on higher-order finite volume methods with a variety of reconstruction algorithms, Riemann solvers, and time integration methods.

## Getting Started

The code is designed to be user-friendly with as few external dependencies as possible.

Documention is permanently under construction on the [wiki](https://github.com/IAS-Astrophysics/athenak/wiki) pages.

In particular, see the complete list of [requirements](https://github.com/IAS-Astrophysics/athenak/wikis/Requirements), or
instructions on how to [download](https://github.com/IAS-Astrophysics/athenak/wikis/Download) and [build](https://github.com/IAS-Astrophysics/athenak/wikis/Build) the code for various devices.
Other pages give instructions for running the code.

Since AthenaK is very similar to Athena++, the [Athena++ documention](https://github.com/PrincetonUniversity/athena/wiki) may also be helpful.

## GOTHAM initial conditions

`gotham_dev` pins the matching generator in `external/CGM_ICs`. From the repository root:

```bash
git submodule update --init --recursive
python -m venv .venv-ics
.venv-ics/bin/python -m pip install -r external/CGM_ICs/requirements.txt
.venv-ics/bin/python external/CGM_ICs/generate_gotham_ics.py \
    inputs/gotham/athinput.gotham build/gotham-run --r-sonic 0.7
```

The generator creates a new run directory containing the profiles, clusters,
`athinput.gotham`, and `ic_parameters.json`. It refuses to overwrite an existing
directory. Build AthenaK with `-DPROBLEM=gotham`, then run the resulting executable
from `build/gotham-run` with `-i athinput.gotham`. Table paths are relative to that
directory. The template's run length and supernova delay are preserved.

The CGM temperature is matched at $R_b=1.5R_{200}$ by default, with
$c_s^2(R_b)=f_Tv_c^2(R_b)$ and adiabatic sound speed $c_s$. Use `--f-t` (default 1)
and `--r-boundary` (kpc) to change this condition. The matching radius is
independent of the mesh; the exported profile must still cover its corner ghost
cells. See the [1 and 20 kpc sonic-radius comparison](external/CGM_ICs/docs/review/outer_temperature_comparison.md)
for results at $f_T=0.5,1,2$, including cases that fail that coverage requirement.

Alternatively, `--cooling-to-flow` imposes $t_{\rm cool}/t_{\rm flow}=1$ at the
same radius; supply a value such as `--cooling-to-flow 2` to change the ratio.
This replaces `--f-t`. Cooling time uses the net cooling rate and
$t_{\rm flow}=r/|v_r|$. See the [comparison of both boundary conditions](external/CGM_ICs/docs/review/outer_boundary_comparison.md)
for ratios 0.5, 1, and 2. The ratio need not select a unique flow: the generator
rejects detected multiple roots instead of choosing a branch automatically.

See the generator's [physics review and export conventions](external/CGM_ICs/docs/review/gotham_integration.md)
for the corrections, retained model parameters, and remaining approximations.
The ICs use a spherical cooling flow plus an isothermal disk; they are not an exact
three-dimensional rotating equilibrium.

To check the generator and its agreement with the compiled GOTHAM physics:

```bash
.venv-ics/bin/python -m pip install -r external/CGM_ICs/requirements-dev.txt
.venv-ics/bin/python -m pytest -q external/CGM_ICs/tests
cd tst
../.venv-ics/bin/python -m scripts.hydro.gotham_ics
```

The last check needs a C++17 compiler. For reproducible runs, use the pinned
submodule commit (`git submodule update --init`), without `--remote`. When updating
the generator, publish its commit before committing the new gitlink in AthenaK.

## Code papers

For more details on the features and algorithms implemented in AthenaK, see the code papers:
- [Stone et al (2024)](https://ui.adsabs.harvard.edu/abs/2024arXiv240916053S/abstract): basic framework
- [Zhu et al. (2024)](https://ui.adsabs.harvard.edu/abs/2024arXiv240910383Z/abstract): numerical relativity solver
- [Fields at al. (2024)](https://ui.adsabs.harvard.edu/abs/2024arXiv240910384F/abstract): GR hydro and MHD solver in dynamical spacetimes

Please reference these papers as appropriate for any publications that use AthenaK.
