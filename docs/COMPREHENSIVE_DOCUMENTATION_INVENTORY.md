# AthenaK Codebase Documentation Inventory

## Overview
This document provides a comprehensive structured inventory of all files in the AthenaK `/src` directory, organized by module with descriptions, key components, and dependencies.

## Core System Files

### Main Entry Point
- **`main.cpp`** - Main application entry point
  - Purpose: Program startup and initialization
  - Key components: `main()`, command line parsing, Kokkos initialization
  - Used by: Executable entry point

### System-Wide Headers
- **`athena.hpp`** - Core type definitions and Kokkos abstractions 
  - Purpose: Central header with fundamental types, enums, and parallel execution patterns
  - Key components: Real type alias, Kokkos View templates (DvceArray1D-6D), parallel loop wrappers (par_for), variable indices (IDN, IM1, etc.), face/edge field structs
  - Used by: Every module in the codebase

- **`athena_tensor.hpp`** - Tensor operations for general relativity
  - Purpose: 4D tensor manipulation utilities for spacetime calculations
  - Key components: 4-index tensor operations, metric manipulations
  - Used by: Z4c, coordinates, general relativity modules

- **`globals.hpp/cpp`** - Global variables
  - Purpose: MPI rank information and global state
  - Key components: `global_variable::my_rank`, `global_variable::nranks`
  - Used by: All MPI-parallel operations

- **`parameter_input.hpp/cpp`** - Input file parsing
  - Purpose: Parse and manage parameters from .athinput files
  - Key components: `ParameterInput` class, `InputBlock`, `InputLine`, parameter accessors
  - Used by: All modules for configuration

## Module-by-Module Inventory

## mesh/
Core mesh infrastructure and adaptive mesh refinement.

### Header Files
- **`mesh.hpp`** - Main mesh class and data structures
  - Purpose: Overall grid structure management, domain decomposition
  - Key components: `Mesh` class, `RegionSize`, `RegionIndcs`, `NeighborBlock`, `LogicalLocation`, `EventCounters`
  - Used by: All physics modules, driver, outputs

- **`meshblock.hpp`** - Individual mesh block definition
  - Purpose: Single mesh block (patch) data structure
  - Key components: `MeshBlock` class, local grid information
  - Used by: Mesh, physics modules

- **`meshblock_pack.hpp`** - Collections of mesh blocks for GPU performance  
  - Purpose: Group mesh blocks together for efficient GPU kernels
  - Key components: `MeshBlockPack` class, pack management
  - Used by: All physics modules

- **`meshblock_tree.hpp`** - Tree structure for AMR
  - Purpose: Hierarchical tree for adaptive mesh refinement
  - Key components: `MeshBlockTree` class, tree operations
  - Used by: Mesh, mesh refinement

- **`mesh_refinement.hpp`** - AMR logic and algorithms
  - Purpose: Adaptive mesh refinement criteria and operations
  - Key components: `MeshRefinement` class, refinement functions
  - Used by: Mesh, driver

- **`nghbr_index.hpp`** - Neighbor indexing utilities
  - Purpose: Fast neighbor finding and indexing
  - Key components: Neighbor index calculations
  - Used by: Boundary values, mesh operations

- **`prolongation.hpp`** - Interpolation between refinement levels
  - Purpose: Transfer data between coarse and fine grids
  - Key components: Prolongation operators
  - Used by: Boundary values, AMR

- **`restriction.hpp`** - Coarsening operations
  - Purpose: Transfer data from fine to coarse grids  
  - Key components: Restriction operators
  - Used by: Boundary values, AMR

### Source Files
- **`mesh.cpp`** - Mesh class implementation
- **`meshblock.cpp`** - MeshBlock implementation
- **`meshblock_pack.cpp`** - MeshBlockPack implementation
- **`meshblock_tree.cpp`** - Tree operations
- **`build_tree.cpp`** - Tree construction algorithms
- **`load_balance.cpp`** - Load balancing for parallel execution
- **`mesh_refinement.cpp`** - AMR implementation

## hydro/
Non-relativistic hydrodynamics solver for Euler equations.

### Header Files
- **`hydro.hpp`** - Main hydro class
  - Purpose: Hydrodynamics evolution and data management
  - Key components: `Hydro` class, `HydroTaskIDs`, Riemann solver enums, variable arrays (u0, w0), boundary values, task functions
  - Used by: Driver, task lists, boundary values

### Riemann Solvers (`rsolvers/`)
- **`advect_hyd.hpp`** - Simple advection solver
- **`hlle_hyd.hpp`** - HLLE approximate Riemann solver
- **`hllc_hyd.hpp`** - HLLC Riemann solver  
- **`roe_hyd.hpp`** - Roe linearized solver
- **`llf_hyd.hpp`** - Local Lax-Friedrichs solver
- **`hlle_grhyd.hpp`** - HLLE for general relativistic hydro
- **`llf_grhyd.hpp`** - LLF for general relativistic hydro
- **`hlle_srhyd.hpp`** - HLLE for special relativistic hydro
- **`hllc_srhyd.hpp`** - HLLC for special relativistic hydro  
- **`llf_srhyd.hpp`** - LLF for special relativistic hydro
- **`llf_hyd_singlestate.hpp`** - Single-state LLF solver

### Source Files
- **`hydro.cpp`** - Main hydro implementation
- **`hydro_fluxes.cpp`** - Flux calculations
- **`hydro_fofc.cpp`** - First-order flux correction
- **`hydro_newdt.cpp`** - Time step calculation
- **`hydro_tasks.cpp`** - Task list integration
- **`hydro_update.cpp`** - Conservative variable updates

## mhd/
Magnetohydrodynamics solver with constrained transport.

### Header Files  
- **`mhd.hpp`** - Main MHD class
  - Purpose: MHD evolution with magnetic field constraint transport
  - Key components: `MHD` class, `MHDTaskIDs`, magnetic field arrays (b0, bcc0), electric field arrays, CT operations, MHD boundary functions
  - Used by: Driver, task lists, boundary values

### Riemann Solvers (`rsolvers/`)
- **`advect_mhd.hpp`** - Simple advection for MHD
- **`hlle_mhd.hpp`** - HLLE MHD solver
- **`hlld_mhd.hpp`** - HLLD MHD solver (handles contact/Alfven waves)
- **`llf_mhd.hpp`** - Local Lax-Friedrichs MHD
- **`hlle_grmhd.hpp`** - HLLE for general relativistic MHD
- **`llf_grmhd.hpp`** - LLF for general relativistic MHD  
- **`hlle_srmhd.hpp`** - HLLE for special relativistic MHD
- **`llf_srmhd.hpp`** - LLF for special relativistic MHD
- **`llf_mhd_singlestate.hpp`** - Single-state LLF MHD solver

### Source Files
- **`mhd.cpp`** - Main MHD implementation
- **`mhd_fluxes.cpp`** - MHD flux calculations  
- **`mhd_ct.cpp`** - Constrained transport implementation
- **`mhd_corner_e.cpp`** - Corner electric field calculations
- **`mhd_fofc.cpp`** - First-order flux correction for MHD
- **`mhd_newdt.cpp`** - MHD time step calculation
- **`mhd_tasks.cpp`** - MHD task list integration  
- **`mhd_update.cpp`** - MHD variable updates

## bvals/
Boundary values and inter-block communication.

### Header Files
- **`bvals.hpp`** - Boundary value classes and MPI communication
  - Purpose: Handle all boundary conditions and inter-block data exchange
  - Key components: `MeshBoundaryValues` base class, `MeshBoundaryValuesCC/FC` derived classes, `MeshBoundaryBuffer`, `ParticlesBoundaryValues`, MPI communication functions, boundary condition enums
  - Used by: All physics modules, mesh operations

### Source Files  
- **`bvals.cpp`** - Base boundary value implementation
- **`bvals_cc.cpp`** - Cell-centered boundary values
- **`bvals_fc.cpp`** - Face-centered boundary values  
- **`bvals_part.cpp`** - Particle boundary values
- **`bvals_tasks.cpp`** - Boundary value task integration
- **`buffs_cc.cpp`** - Cell-centered communication buffers
- **`buffs_fc.cpp`** - Face-centered communication buffers
- **`flux_correct_cc.cpp`** - Flux correction for cell-centered variables
- **`flux_correct_fc.cpp`** - Flux correction for face-centered variables
- **`prolongation.cpp`** - Prolongation operations
- **`prolong_prims.cpp`** - Primitive variable prolongation

### Physics-Specific Boundary Conditions (`physics/`)
- **`hydro_bcs.cpp`** - Hydro boundary conditions
- **`bfield_bcs.cpp`** - Magnetic field boundary conditions  
- **`radiation_bcs.cpp`** - Radiation boundary conditions
- **`z4c_bcs.cpp`** - Z4c spacetime boundary conditions

## coordinates/
Coordinate systems and metric calculations.

### Header Files
- **`coordinates.hpp`** - Coordinate system management
  - Purpose: Handle different coordinate systems and spacetime metrics
  - Key components: `Coordinates` class, `CoordData`, relativistic flags, excision masks, coordinate source terms
  - Used by: All physics modules, metric calculations

- **`adm.hpp`** - ADM formalism for general relativity
  - Purpose: 3+1 decomposition of spacetime
  - Key components: ADM variables, lapse, shift
  - Used by: Z4c, general relativity modules

- **`cartesian_ks.hpp`** - Cartesian Kerr-Schild coordinates
- **`cell_locations.hpp`** - Cell position calculations

### Source Files
- **`coordinates.cpp`** - Main coordinate implementation
- **`adm.cpp`** - ADM implementation
- **`excision.cpp`** - Black hole excision algorithms

## driver/
Main simulation driver and time integration.

### Header Files
- **`driver.hpp`** - Main simulation driver
  - Purpose: Orchestrate entire simulation execution and time integration
  - Key components: `Driver` class, time evolution control, integrator weights, task execution, wall clock timing
  - Used by: Main function, task orchestration

### Source Files  
- **`driver.cpp`** - Driver implementation

## eos/
Equations of state for different fluids and relativity.

### Header Files
- **`eos.hpp`** - Equation of state base classes and derived implementations
  - Purpose: Convert between primitive and conserved variables for different EOS
  - Key components: `EquationOfState` base class, `EOS_Data` struct, derived classes for isothermal/ideal gas in Newtonian/SR/GR, sound speed functions, wave speed calculations
  - Used by: All fluid physics modules

### Primitive Solver (`primitive-solver/`)
- **`primitive_solver.hpp`** - Advanced primitive recovery
- **`eos_compose.hpp`** - CompOSE EOS interface
- **`piecewise_polytrope.hpp`** - Piecewise polytropic EOS
- **`unit_system.hpp`** - Physical unit systems
- **`error_policy_interface.hpp`** - Error handling for primitive recovery
- **`eos_policy_interface.hpp`** - EOS policy interface
- **`ps_error.hpp`** - Primitive solver error handling
- **`ps_types.hpp`** - Primitive solver type definitions
- **`reset_floor.hpp`** - Variable floor application
- **`geom_math.hpp`** - Geometric mathematical utilities
- **`idealgas.hpp`** - Ideal gas utilities
- **`numtools_root.hpp`** - Root finding utilities

### Source Files
- **`eos.cpp`** - Base EOS implementation
- **`ideal_hyd.cpp`** - Ideal gas non-relativistic hydro
- **`ideal_mhd.cpp`** - Ideal gas non-relativistic MHD
- **`ideal_srhyd.cpp`** - Ideal gas special relativistic hydro
- **`ideal_srmhd.cpp`** - Ideal gas special relativistic MHD  
- **`ideal_grhyd.cpp`** - Ideal gas general relativistic hydro
- **`ideal_grmhd.cpp`** - Ideal gas general relativistic MHD
- **`isothermal_hyd.cpp`** - Isothermal non-relativistic hydro
- **`isothermal_mhd.cpp`** - Isothermal non-relativistic MHD

### Additional Headers
- **`ideal_c2p_hyd.hpp`** - Ideal gas hydro conservative-to-primitive
- **`ideal_c2p_mhd.hpp`** - Ideal gas MHD conservative-to-primitive  
- **`primitive_solver_hyd.hpp`** - Hydro primitive solver interface

## outputs/
Data output in various formats.

### Header Files
- **`outputs.hpp`** - Output system management and all output types
  - Purpose: Handle all types of simulation data output
  - Key components: `Outputs` class, `BaseTypeOutput` base class, derived output classes (VTK, Binary, History, etc.), `OutputParameters`, output variable choices array, history data management
  - Used by: Driver, analysis tools

- **`io_wrapper.hpp`** - I/O abstraction layer
  - Purpose: Abstract file I/O operations
  - Key components: File I/O wrapper functions
  - Used by: All output modules

### Source Files
- **`outputs.cpp`** - Main output system
- **`vtk_mesh.cpp`** - VTK mesh output
- **`vtk_prtcl.cpp`** - VTK particle output
- **`binary.cpp`** - Binary mesh data output
- **`coarsened_binary.cpp`** - Coarsened binary output  
- **`restart.cpp`** - Restart file output
- **`history.cpp`** - History/time series output
- **`formatted_table.cpp`** - Tabular data output
- **`basetype_output.cpp`** - Base output type implementation
- **`derived_variables.cpp`** - Derived variable calculations
- **`eventlog.cpp`** - Event logging output
- **`pdf.cpp`** - Probability distribution function output
- **`track_prtcl.cpp`** - Tracked particle output
- **`pos_prtcl.cpp`** - Particle position output
- **`df_prtcl.cpp`** - Particle distribution function output
- **`dxhist_prtcl.cpp`** - Particle displacement histogram output  
- **`rst_prtcl.cpp`** - Particle restart output
- **`io_wrapper.cpp`** - I/O wrapper implementation

## particles/
Lagrangian particle tracking and cosmic ray transport.

### Header Files
- **`particles.hpp`** - Particle system management
  - Purpose: Lagrangian particle evolution and tracking
  - Key components: `Particles` class, `ParticlesTaskIDs`, particle pusher enums, particle data arrays (position, velocity, properties), boundary values
  - Used by: Driver, outputs, boundary values

### Source Files
- **`particles.cpp`** - Main particle implementation
- **`particles_pushers.cpp`** - Particle integration schemes  
- **`particles_tasks.cpp`** - Particle task integration

## reconstruct/
Spatial reconstruction methods for high-order accuracy.

### Header Files
- **`dc.hpp`** - Donor cell (first-order) reconstruction
- **`plm.hpp`** - Piecewise linear method (second-order)
- **`ppm.hpp`** - Piecewise parabolic method (third-order)  
- **`wenoz.hpp`** - WENO-Z reconstruction (high-order)

## radiation/
Radiation transport with moment methods.

### Header Files
- **`radiation.hpp`** - Radiation transport system
  - Purpose: Radiation moment transport and fluid-radiation coupling
  - Key components: `Radiation` class, `RadiationTaskIDs`, intensity arrays, tetrad calculations, angular mesh, opacity parameters, radiation source terms
  - Used by: Driver, task lists, outputs

- **`radiation_tetrad.hpp`** - Tetrad frame calculations
  - Purpose: Coordinate transformations for radiation transport
  - Key components: Tetrad construction and transformations
  - Used by: Radiation module

- **`radiation_opacities.hpp`** - Opacity calculations
  - Purpose: Radiation-matter interaction coefficients
  - Key components: Opacity functions and tables
  - Used by: Radiation module

### Source Files
- **`radiation.cpp`** - Main radiation implementation
- **`radiation_fluxes.cpp`** - Radiation flux calculations
- **`radiation_newdt.cpp`** - Radiation time step limits
- **`radiation_source.cpp`** - Radiation source terms
- **`radiation_tasks.cpp`** - Radiation task integration
- **`radiation_tetrad.cpp`** - Tetrad implementation
- **`radiation_update.cpp`** - Radiation variable updates

## srcterms/
Source terms for various physical effects.

### Header Files  
- **`srcterms.hpp`** - Source term management
  - Purpose: Physical source terms (gravity, cooling, turbulence, etc.)
  - Key components: `SourceTerms` class, constant acceleration, cooling functions, shearing box, turbulence driver interface
  - Used by: Hydro, MHD, driver

- **`turb_driver.hpp`** - Turbulence driving
  - Purpose: Stochastic turbulence forcing  
  - Key components: `TurbulenceDriver` class, Fourier space forcing
  - Used by: Source terms

- **`cooling_tables.hpp`** - Cooling rate tables
- **`ismcooling.hpp`** - Interstellar medium cooling

### Source Files
- **`srcterms.cpp`** - Main source term implementation
- **`srcterms_newdt.cpp`** - Source term time step limits
- **`turb_driver.cpp`** - Turbulence driver implementation

### Additional Files
- **`README.md`** - Source terms documentation
- **`turb_driver_diagnostic.patch`** - Diagnostic patch
- **`turb_driver_old_working.cpp`** - Legacy implementation

## diffusion/
Physical diffusion processes.

### Header Files
- **`viscosity.hpp`** - Viscous stress tensor
  - Purpose: Viscous diffusion in fluid flows
  - Key components: Viscosity coefficients, stress calculations
  - Used by: Hydro, MHD

- **`resistivity.hpp`** - Magnetic diffusion  
  - Purpose: Ohmic dissipation of magnetic fields
  - Key components: Resistivity coefficients, current calculations
  - Used by: MHD

- **`conduction.hpp`** - Thermal conduction
  - Purpose: Heat diffusion processes
  - Key components: Thermal conductivity, heat flux calculations  
  - Used by: Hydro, MHD

- **`current_density.hpp`** - Current density calculations

### Source Files
- **`viscosity.cpp`** - Viscosity implementation
- **`resistivity.cpp`** - Resistivity implementation  
- **`conduction.cpp`** - Conduction implementation

## z4c/
Z4c formulation for numerical relativity.

### Header Files
- **`z4c.hpp`** - Z4c spacetime evolution
  - Purpose: Numerical relativity using Z4c formalism
  - Key components: `Z4c` class, `Z4cTaskIDs`, evolved variable indices, constraint calculations, gauge conditions
  - Used by: Driver, outputs, coordinates

- **`z4c_amr.hpp`** - Z4c adaptive mesh refinement
- **`z4c_macros.hpp`** - Z4c utility macros  
- **`compact_object_tracker.hpp`** - Black hole tracking
- **`tmunu.hpp`** - Stress-energy tensor

### Source Files
- **`z4c.cpp`** - Main Z4c implementation
- **`z4c_adm.cpp`** - ADM variable calculations
- **`z4c_amr.cpp`** - Z4c AMR implementation  
- **`z4c_calcrhs.cpp`** - Right-hand side calculations
- **`z4c_gauge.cpp`** - Gauge condition implementation
- **`z4c_newdt.cpp`** - Z4c time step calculation
- **`z4c_tasks.cpp`** - Z4c task integration
- **`z4c_update.cpp`** - Z4c variable updates
- **`z4c_Sbc.cpp`** - Sommerfeld boundary conditions
- **`z4c_calculate_weyl_scalars.cpp`** - Weyl scalar calculations
- **`z4c_wave_extr.cpp`** - Gravitational wave extraction
- **`compact_object_tracker.cpp`** - Compact object tracking
- **`tmunu.cpp`** - Stress-energy tensor implementation

## dyn_grmhd/
General relativistic MHD in dynamical spacetimes.

### Header Files
- **`dyn_grmhd.hpp`** - Dynamic spacetime GRMHD
  - Purpose: GRMHD evolution in dynamic spacetime backgrounds  
  - Key components: `DynGRMHD` class, spacetime coupling
  - Used by: Z4c, MHD

- **`dyn_grmhd_util.hpp`** - GRMHD utilities

### Riemann Solvers (`rsolvers/`)  
- **`flux_dyn_grmhd.hpp`** - Flux calculations for dynamic GRMHD
- **`hlle_dyn_grmhd.hpp`** - HLLE solver for dynamic GRMHD
- **`llf_dyn_grmhd.hpp`** - LLF solver for dynamic GRMHD

### Source Files
- **`dyn_grmhd.cpp`** - Main dynamic GRMHD implementation
- **`dyn_grmhd_fluxes.cpp`** - Dynamic GRMHD flux calculations
- **`dyn_grmhd_fofc.cpp`** - First-order flux correction

## shearing_box/
Shearing box simulations for accretion disk physics.

### Header Files
- **`shearing_box.hpp`** - Shearing box framework
  - Purpose: Local shear flow simulations for accretion disks
  - Key components: `ShearingBox` class, orbital advection, remapping
  - Used by: Hydro, MHD

- **`remap_fluxes.hpp`** - Flux remapping utilities

### Source Files  
- **`shearing_box.cpp`** - Main shearing box implementation
- **`shearing_box_cc.cpp`** - Cell-centered shearing box  
- **`shearing_box_fc.cpp`** - Face-centered shearing box
- **`shearing_box_srcterms.cpp`** - Shearing box source terms
- **`shearing_box_tasks.cpp`** - Shearing box task integration
- **`orbital_advection.cpp`** - Orbital advection implementation
- **`orbital_advection_cc.cpp`** - Cell-centered orbital advection
- **`orbital_advection_fc.cpp`** - Face-centered orbital advection

## ion-neutral/
Two-fluid ion-neutral plasma physics.

### Header Files
- **`ion-neutral.hpp`** - Ion-neutral fluid coupling
  - Purpose: Two-fluid plasma with ion-neutral interactions
  - Key components: `IonNeutral` class, collision terms
  - Used by: Driver, two-fluid simulations

### Source Files
- **`ion-neutral.cpp`** - Ion-neutral implementation  
- **`ion-neutral_tasks.cpp`** - Ion-neutral task integration

## geodesic-grid/
Spherical grids for radiation transport and analysis.

### Header Files
- **`geodesic_grid.hpp`** - Geodesic grid construction
  - Purpose: Spherical angular grids for radiation transport
  - Key components: `GeodesicGrid` class, icosahedral subdivision
  - Used by: Radiation module

- **`spherical_grid.hpp`** - Spherical coordinate grids  
  - Purpose: Spherical coordinate analysis tools
  - Key components: `SphericalGrid` class
  - Used by: Analysis, problem generators

### Source Files
- **`geodesic_grid.cpp`** - Geodesic grid implementation
- **`spherical_grid.cpp`** - Spherical grid implementation

## tasklist/
Task list management for execution scheduling.

### Header Files
- **`task_list.hpp`** - Task scheduling framework
  - Purpose: Manage execution dependencies and scheduling
  - Key components: Task definitions, dependency management
  - Used by: Driver, all physics modules

- **`numerical_relativity.hpp`** - Numerical relativity task lists
  - Purpose: Specialized task lists for GR simulations
  - Key components: NR-specific task scheduling
  - Used by: Z4c, numerical relativity

### Source Files  
- **`numerical_relativity.cpp`** - NR task list implementation

## units/
Physical unit systems and conversions.

### Header Files
- **`units.hpp`** - Unit system management
  - Purpose: Handle physical units and conversions
  - Key components: `Units` class, unit conversions
  - Used by: EOS, problem generators

### Source Files
- **`units.cpp`** - Unit system implementation

## utils/
General utility functions and tools.

### Header Files
- **`utils.hpp`** - General utilities
  - Purpose: Common utility functions used throughout codebase
  - Key components: Mathematical utilities, helper functions
  - Used by: Many modules

- **`finite_diff.hpp`** - Finite difference operations
- **`current.hpp`** - Current density calculations  
- **`random.hpp`** - Random number generation
- **`profile_reader.hpp`** - Profile data reading
- **`lagrange_interpolator.hpp`** - Lagrange interpolation
- **`tr_table.hpp`** - Table reading utilities
- **`tr_utils.hpp`** - Table utilities

### Source Files
- **`change_rundir.cpp`** - Runtime directory changes
- **`show_config.cpp`** - Configuration display
- **`lagrange_interpolator.cpp`** - Lagrange interpolation implementation
- **`tr_table.cpp`** - Table implementation

## pgen/
Problem generators for initial conditions.

### Header Files
- **`pgen.hpp`** - Problem generator framework
  - Purpose: Set up initial conditions for various astrophysical problems
  - Key components: `ProblemGenerator` class, user function pointers, predefined problems, spherical grids for analysis
  - Used by: Driver, mesh initialization

### Test Problems (`tests/`)
- **`advection.cpp`** - Simple advection test
- **`collapse.cpp`** - Gravitational collapse
- **`cpaw.cpp`** - Circularly polarized Alfven wave
- **`diffusion.cpp`** - Diffusion test problems
- **`gr_bondi.cpp`** - General relativistic Bondi accretion
- **`gr_monopole.cpp`** - GR monopole test
- **`linear_wave.cpp`** - Linear wave tests
- **`lw_implode.cpp`** - Liska-Wendroff implosion
- **`orszag_tang.cpp`** - Orszag-Tang MHD vortex
- **`rad_check_tetrad.cpp`** - Radiation tetrad verification
- **`rad_hohlraum.cpp`** - Radiation hohlraum
- **`rad_linear_wave.cpp`** - Radiation linear waves
- **`shock_tube.cpp`** - Shock tube problems
- **`z4c_linear_wave.cpp`** - Z4c gravitational wave tests

### Astrophysical Problems
- **`blast.cpp`** - Blast wave problems
- **`field_loop.cpp`** - Magnetic field loop advection
- **`gr_torus.cpp`** - General relativistic torus
- **`kh.cpp`** - Kelvin-Helmholtz instability
- **`rt.cpp`** - Rayleigh-Taylor instability  
- **`shock_cloud.cpp`** - Shock-cloud interaction
- **`shu_osher.cpp`** - Shu-Osher shock tube
- **`slotted_cyl.cpp`** - Slotted cylinder advection
- **`turb.cpp`** - Turbulence problems
- **`thermal_instability.cpp`** - Thermal instability

### Specialized Problems
- **`cgm_cooling_flow.cpp`** - Circumgalactic medium cooling
- **`current_sheet.cpp`** - Magnetic current sheet
- **`cshock.cpp`** - C-type shock
- **`dynamo.cpp`** - Dynamo simulations  
- **`dyngr_tov.cpp`** - Dynamic GR Tolman-Oppenheimer-Volkoff
- **`elliptica_bns.cpp`** - Binary neutron star (Elliptica)
- **`hydrostatic_1d.cpp`** - 1D hydrostatic equilibrium
- **`hydrostatic_3d.cpp`** - 3D hydrostatic equilibrium
- **`lorene_bns.cpp`** - Binary neutron star (LORENE)
- **`mass_removal_test.cpp`** - Mass removal test
- **`mri2d.cpp`** - 2D magnetorotational instability
- **`mri3d.cpp`** - 3D magnetorotational instability
- **`part_random.cpp`** - Random particle distribution
- **`part_static_turb.cpp`** - Static turbulent particles
- **`sgrid_bns.cpp`** - Binary neutron star (SpEC grid)
- **`sfb_turb.cpp`** - Spherical Fourier-Bessel turbulence
- **`sfb_turb_amr_test.cpp`** - SFB turbulence AMR test
- **`turbulent_box.cpp`** - Turbulent box simulations
- **`twofluid.cpp`** - Two-fluid problems
- **`z4c_one_puncture.cpp`** - Single black hole Z4c
- **`z4c_spectre_bbh.cpp`** - Binary black hole from SpECTRE
- **`z4c_two_puncture.cpp`** - Binary black hole Z4c

### Additional Specialized
- **`rad_beam.cpp`** - Radiation beam test
- **`rad_diffusion.cpp`** - Radiation diffusion
- **`rad_relax.cpp`** - Radiation relaxation
- **`rad_shadow.cpp`** - Radiation shadow test
- **`rad_snake.cpp`** - Radiation snake test
- **`rt_mhd_amr.cpp`** - RT instability with MHD and AMR
- **`turb_amr_test.cpp`** - Turbulence AMR test
- **`turb_amr_wave_test.cpp`** - Turbulence AMR wave test
- **`turb_mhd_amr_wave.cpp`** - MHD turbulence AMR wave

## Key Dependencies and Relationships

### Core Dependencies
- **athena.hpp**: Required by all modules for basic types and Kokkos abstractions
- **parameter_input.hpp**: Used by all modules for configuration
- **mesh/**: Foundation for all physics modules
- **bvals/**: Required by all physics modules for boundary conditions

### Physics Module Dependencies
- **hydro/** → eos, coordinates, bvals, mesh, reconstruct
- **mhd/** → hydro dependencies + magnetic field handling
- **radiation/** → hydro/mhd + geodesic-grid for angular mesh
- **z4c/** → coordinates, mesh, bvals for spacetime evolution
- **particles/** → mesh, bvals for Lagrangian tracking

### Utility Dependencies  
- **outputs/** → All physics modules for data extraction
- **driver/** → All modules for orchestration
- **reconstruct/** → hydro, mhd, radiation for high-order accuracy
- **srcterms/** → hydro, mhd for additional physics

This inventory provides a comprehensive overview of the AthenaK codebase structure, with each module's purpose, key components, and usage relationships clearly documented for efficient navigation and understanding.