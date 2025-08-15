.. _glossary:

========
Glossary
========

.. glossary::
   :sorted:

   AMR
      Adaptive Mesh Refinement - Dynamic mesh refinement technique to resolve multiple scales in simulations.
      See :doc:`modules/mesh` for implementation details.

   MeshBlock
      Fundamental unit of domain decomposition in AthenaK. Each MeshBlock contains a fixed number of cells
      and can be refined or coarsened during AMR. See :doc:`modules/mesh`.

   Kokkos
      Performance portability layer that allows AthenaK to run on CPUs and GPUs from a single codebase.
      See :doc:`kokkos_guide`.

   TaskList
      Task-based execution system that manages computational dependencies automatically.
      See :doc:`modules/tasklist`.

   Riemann Solver
      Numerical methods for solving discontinuities in hyperbolic conservation laws.
      Available solvers include LLF, HLLC, HLLD, and Roe. See :doc:`modules/riemann_solvers`.

   EOS
      Equation of State - Closure relation for fluid dynamics equations.
      Options include ideal gas, isothermal, and primitive solver. See :doc:`modules/eos`.

   MHD
      Magnetohydrodynamics - Physics module for magnetic field evolution.
      See :doc:`modules/mhd`.

   GRMHD
      General Relativistic Magnetohydrodynamics - MHD in curved spacetime.
      See :doc:`modules/dyn_grmhd`.

   Z4c
      Numerical relativity formalism for evolving Einstein's equations.
      See :doc:`modules/z4c`.

   Problem Generator
      User-defined function that sets initial conditions for simulations.
      Located in ``src/pgen/``. See :doc:`modules/pgen`.

   Boundary Conditions
      Methods for handling domain boundaries: periodic, outflow, reflecting, user-defined.
      See :doc:`modules/boundaries`.

   Ghost Cells
      Extra cells at MeshBlock boundaries for applying boundary conditions and MPI communication.
      Typically 2-4 cells deep.

   CFL Condition
      Courant-Friedrichs-Lewy stability condition that limits timestep based on wave speeds.
      Set via ``cfl_number`` parameter.

   Reconstruction
      Spatial interpolation methods: PLM (piecewise linear), PPM (piecewise parabolic), WENOZ.
      See :doc:`modules/reconstruction`.

   CT
      Constrained Transport - Method for maintaining divergence-free magnetic fields in MHD.
      
   FOFC
      First-Order Flux Correction - Fallback to first-order accuracy in problematic cells.

   Particles
      Lagrangian particle tracking module for test particles, cosmic rays, or dust.
      See :doc:`modules/particles`.

   Source Terms
      Additional physics: gravity, turbulence driving, cooling, etc.
      See :doc:`modules/srcterms`.

   Turbulence Driver
      Stochastic forcing for maintaining turbulence in simulations.
      Available in Fourier and Spherical Fourier-Bessel modes.

   Outputs
      Data output formats: VTK, HDF5, binary, restart files.
      See :doc:`modules/outputs`.

   Restart
      Checkpoint files for resuming simulations. Created with ``file_type = rst``.

   Load Balancing
      Dynamic redistribution of MeshBlocks across MPI ranks for optimal performance.

   Primitives
      Physical variables (density, velocity, pressure) vs conserved variables (momentum, energy).

   C2P
      Conserved-to-Primitive variable conversion, critical in relativistic codes.

   Coordinates
      Coordinate systems: Cartesian, spherical, cylindrical, Cartesian-Kerr-Schild.
      See :doc:`modules/coordinates`.

   Shearing Box
      Boundary conditions for local disk simulations with orbital shear.
      See :doc:`modules/shearing_box`.

   Diffusion
      Physical dissipation processes: viscosity, resistivity, thermal conduction.
      See :doc:`modules/diffusion`.

   Radiation
      Radiation transport module using M1 closure or other methods.
      See :doc:`modules/radiation`.

   Ion-Neutral
      Two-fluid physics for partially ionized plasmas.
      See :doc:`modules/ion_neutral`.

   CGM
      Circumgalactic Medium - Gas surrounding galaxies, focus of cooling flow problems.

   ISM
      Interstellar Medium - Gas and dust between stars.

   SN
      Supernova - Stellar explosion that injects energy, mass, and metals.

   Input File
      ``.athinput`` file containing simulation parameters organized in blocks.

   Parameter Input
      Runtime configuration system using hierarchical parameter blocks.

.. index::
   single: Adaptive Mesh Refinement
   single: MeshBlock
   single: Kokkos
   single: TaskList
   single: Riemann solver
   single: Equation of State
   single: Magnetohydrodynamics
   single: General Relativity
   single: Numerical Relativity
   single: Problem Generator
   single: Boundary Conditions
   single: Ghost Cells
   single: CFL Condition
   single: Reconstruction Methods
   single: Constrained Transport
   single: Particles
   single: Source Terms
   single: Turbulence
   single: Output Formats
   single: Restart Files
   single: Load Balancing
   single: Coordinates
   single: Shearing Box
   single: Diffusion
   single: Radiation Transport
   single: Two-fluid
   single: Circumgalactic Medium
   single: Interstellar Medium
   single: Supernova Feedback