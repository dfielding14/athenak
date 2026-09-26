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

## Super-time-stepping for diffusion

RKL2 super-time-stepping (STS) is available for Newtonian hydro and MHD. Enable
it globally and select the diffusion terms to advance with STS, for example:

```ini
<time>
sts_integrator = rkl2
sts_max_dt_ratio = 20  # optional cap relative to the explicit diffusion timestep

<mhd>
ohmic_resistivity = 0.01
resistivity_integrator = sts
```

For either `<hydro>` or `<mhd>`, use `viscosity_integrator = sts` with positive
`viscosity`, or `conductivity_integrator = sts` with positive `conductivity`.
Each selector defaults to `explicit`, so selected and explicit terms can coexist.
The global default is `sts_integrator = none`; `sts_max_dt_ratio = -1` removes
the optional cap. Hyperbolic and source-term timestep limits still apply.

STS uses second-order RKL2 half-sweeps around the ordinary RK update, with
diffusion flux/EMF exchanges at each stage and restarts at cycle boundaries.
It supports constant isotropic viscosity, constant unsaturated isotropic
conductivity, and constant Ohmic resistivity. Temperature-dependent or saturated
conduction, relativistic physics, radiation, ion-neutral coupling, and
shearing-box/orbital advection are not supported with STS.

## Code papers

For more details on the features and algorithms implemented in AthenaK, see the code papers:
- [Stone et al (2024)](https://ui.adsabs.harvard.edu/abs/2024arXiv240916053S/abstract): basic framework
- [Zhu et al. (2024)](https://ui.adsabs.harvard.edu/abs/2024arXiv240910383Z/abstract): numerical relativity solver
- [Fields at al. (2024)](https://ui.adsabs.harvard.edu/abs/2024arXiv240910384F/abstract): GR hydro and MHD solver in dynamical spacetimes

Please reference these papers as appropriate for any publications that use AthenaK.
