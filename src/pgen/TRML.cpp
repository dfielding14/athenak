//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file TRML.cpp
//! \brief Problem generator for Turbulent Radiative Mixing Layer (TRML) simulations.
//!
//! Written by Rajsekhar Mohapatra 2022-2023
//! (AthenaK version of Drummond Fielding's TRML setup)
//!
//! ================================ PHYSICS OVERVIEW ================================
//!
//! This problem generator sets up a shearing interface between hot and cold gas phases
//! to study turbulent radiative mixing layers. The setup models the interaction between:
//!   - A hot, diffuse phase (e.g., the circumgalactic medium)
//!   - A cold, dense phase (e.g., cool clouds or filaments)
//!
//! Key physical processes:
//!   1. SHEAR FLOW: A velocity gradient across the interface drives Kelvin-Helmholtz
//!      instability and turbulent mixing.
//!   2. RADIATIVE COOLING: Gas cools radiatively according to a cooling function,
//!      which can cause runaway cooling in the mixing layer.
//!   3. THERMAL INSTABILITY: The competition between cooling and mixing determines
//!      whether the cold gas grows (condensation) or shrinks (evaporation).
//!
//! ================================ KEY PARAMETERS ================================
//!
//! Temperature/Density Structure:
//!   - T_hot, T_cold: Hot and cold phase temperatures
//!   - T_peak: Temperature at the peak of the cooling curve (typically ~10^5 K)
//!   - contrast: Temperature (or density) contrast between phases
//!
//! Timescales:
//!   - t_shear = L/v: Shear timescale (eddy turnover time)
//!   - t_cool: Cooling time at T_peak
//!   - xi = t_shear/t_cool: Dimensionless ratio controlling mixing layer evolution
//!     * xi << 1: Cooling-dominated (cold gas grows)
//!     * xi >> 1: Mixing-dominated (cold gas destroyed)
//!
//! Cooling Function Parameters:
//!   - beta_lo, beta_hi: Power-law slopes of cooling function at low/high T
//!   - T_cutoff: Temperature above which cooling is turned off
//!
//! ================================ SPECIAL FEATURES ================================
//!
//! 1. FRAME TRACKING:
//!    The simulation can automatically apply a Galilean velocity shift to keep the
//!    mixing layer interface centered in the domain. This prevents the interface
//!    from drifting into the boundary regions over long simulation times.
//!    Enabled via: use_frame_tracking = true
//!
//! 2. USER-DEFINED TIMESTEP (Cooling CFL):
//!    A cooling CFL condition limits dt so that temperature changes by at most
//!    a fraction (cfl_cool) per timestep, preventing numerical instability in
//!    rapidly cooling regions.
//!
//! 3. DIAGNOSTIC OUTPUT:
//!    Periodic output of min/max temperatures, total mass, cooling rates, etc.
//!    Controlled via: ndiag parameter
//!
//! 4. MHD SUPPORT:
//!    Can be run with magnetic fields. The field geometry and strength are
//!    controlled via plasma_beta and B_direction parameters.
//!
//! ================================ INPUT FILE PARAMETERS ==========================
//!
//! See inputs/hydro/TRML.athinput or inputs/mhd/TRML.athinput for examples.
//!
//! Required in <problem> block:
//!   - user_srcs = true      (enables custom cooling/heating)
//!   - user_dt = true        (enables cooling CFL)
//!   - user_work_in_loop = true  (enables frame tracking and diagnostics)
//!
//! ================================================================================

#include <algorithm>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <cstring>
#include <random>
#include <iomanip>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/cell_locations.hpp"
#include "eos/eos.hpp"
#include "eos/ideal_c2p_hyd.hpp"
#include "eos/ideal_c2p_mhd.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "pgen.hpp"
#include "particles/particles.hpp"
#include "diffusion/conduction.hpp"
#include "srcterms/srcterms.hpp"
#include "srcterms/ismcooling.hpp"
#include "globals.hpp"
#include "units/units.hpp"
#include "utils/marching_cubes.hpp"  // Marching cubes algorithm for isosurface area calculations
#include "utils/random.hpp"
#include "srcterms/turb_driver.hpp"
// #include "turb_init.hpp"
// #include "srcterms/TurbGen.h"


// Number of slabs for z-direction decomposition in parallel reductions.
// Used for computing mean quantities across the mixing layer.
#define NREDUCTION_SLAB 24

namespace {

//----------------------------------------------------------------------------------------
//! \struct pgen_trml
//! \brief Data structure holding all problem-specific parameters for the TRML simulation.
//!
//! This structure stores configuration parameters read from the input file and computed
//! derived quantities needed throughout the simulation. It is allocated once at
//! initialization and persists for the entire simulation run.

struct pgen_trml {
  // ====================================================================================
  // EQUATION OF STATE AND FLOOR VALUES
  // ====================================================================================
  Real gamma_adi;              //!< Adiabatic index (ratio of specific heats, gamma)
  Real dfloor;                 //!< Density floor (minimum allowed density)
  Real pfloor;                 //!< Pressure floor (minimum allowed pressure)

  // ====================================================================================
  // INITIAL STATE PARAMETERS
  // ====================================================================================
  Real rho_0;                  //!< Reference density (hot phase density in code units)
  Real pgas_0;                 //!< Reference pressure (initial pressure, uniform)
  Real contrast;               //!< Temperature contrast ratio (T_hot/T_cold = rho_cold/rho_hot)
  Real velocity;               //!< Shear velocity amplitude (hot phase moves at +v, cold at -v)
  Real t_shear;                //!< Shear timescale = 1/velocity (eddy turnover time)

  // ====================================================================================
  // COOLING FUNCTION PARAMETERS
  // ====================================================================================
  Real xi;                     //!< Dimensionless ratio t_shear/t_cool (key parameter!)
                               //!< Controls whether cooling or mixing dominates
  Real t_cool_start;           //!< Simulation time at which to turn on radiative cooling
  Real beta_lo;                //!< Power-law slope of cooling function below T_peak
  Real beta_hi;                //!< Power-law slope of cooling function above T_peak
  Real alpha_heat;             //!< Derived heating coefficient exponent
  Real heat_coefficient;       //!< Derived heating amplitude coefficient
  Real dt_cutoff;              //!< Minimum allowed timestep (prevents runaway)
  Real cfl_cool;               //!< Cooling CFL number (limits dT/T per timestep)

  // ====================================================================================
  // TEMPERATURE THRESHOLDS
  // ====================================================================================
  Real T_cold;                 //!< Cold phase temperature (= pgas_0/rho_0/contrast)
  Real T_hot;                  //!< Hot phase temperature (= pgas_0/rho_0)
  Real T_peak;                 //!< Peak cooling temperature (mixing layer forms here)
  Real T_cutoff;               //!< Temperature above which cooling is turned off
  Real T_floor;                //!< Temperature floor for adjusting cold gas
  Real epsilon_T;              //!< Temperature tolerance parameter

  //!< Temperature bands for phase classification (logarithmically spaced)
  Real T_peak_lo;              //!< Lower bound of T_peak band (for diagnostics)
  Real T_peak_hi;              //!< Upper bound of T_peak band (for diagnostics)
  Real T_hot_lo;               //!< Lower temperature considered "hot"
  Real T_cold_hi;              //!< Upper temperature considered "cold"

  // ====================================================================================
  // GEOMETRY PARAMETERS
  // ====================================================================================
  Real ztop;                   //!< Upper z-boundary of the domain
  Real zbot;                   //!< Lower z-boundary of the domain
  Real smoothing_thickness;    //!< Width of the tanh transition at the interface
  Real z_interface;            //!< z-position of the initial interface

  // ====================================================================================
  // MHD PARAMETERS
  // ====================================================================================
  Real alpha_magdens;          //!< Power-law index for B-field scaling with density
                               //!< (used when B_dens_ratio = true)

  // ====================================================================================
  // ADAPTIVE MESH REFINEMENT THRESHOLDS
  // ====================================================================================
  Real density_ratio_threshold; //!< Refine if density gradient exceeds this
  Real vel2_rms_threshold;      //!< Refine if velocity RMS exceeds this
  Real T_max_threshold;         //!< Refine if T_max in block exceeds this
  Real T_min_threshold;         //!< Refine if T_min in block is below this

  // ====================================================================================
  // COOLING CONTROL FLAGS AND PARAMETERS
  // ====================================================================================
  int  ndiag;                  //!< Output diagnostics every ndiag timesteps (-1 = off)
  bool cool_subc;              //!< Enable cooling subcycling (multiple cooling steps per hydro step)
  Real cool_subfac;            //!< Subcycling factor (dt_cool = cool_subfac * t_cool)
  bool use_temp_floor_cool;    //!< Enable temperature floor for cooling
  Real T_floor_cool;           //!< Minimum temperature for cooling (cooling off below this)
  bool use_temp_ceiling;       //!< Enable maximum temperature cap
  Real T_ceiling;              //!< Maximum allowed temperature in simulation
  bool adjust_temp_floor;      //!< If T < T_floor, average with neighbors
  bool use_dens_ceiling;       //!< Enable density ceiling for cooling
  Real dens_ceiling;           //!< Disable cooling above this density

  // ====================================================================================
  // DIAGNOSTIC TRACKING VARIABLES
  // ====================================================================================
  Real tot_mass;               //!< Total mass in domain (for conservation tracking)
  Real tot_coolrate;           //!< Total cooling rate (energy loss per unit time)

  // ====================================================================================
  // FRAME TRACKING PARAMETERS
  // Frame tracking applies a Galilean velocity shift to keep the mixing layer
  // interface centered in the domain. This prevents the interface from drifting
  // into boundary regions over long simulation times.
  // ====================================================================================
  bool use_frame_tracking;     //!< Enable frame tracking
  Real t_frame_tracking_start; //!< Time at which to start frame tracking
  int  n_frame_track;          //!< Apply frame tracking every n_frame_track timesteps
  Real max_vz_frame_tracking;  //!< Maximum velocity shift allowed per application
  Real T_ft_lo;                //!< Lower temperature bound for tracking (track gas in this range)
  Real T_ft_hi;                //!< Upper temperature bound for tracking
  bool ft_use_uppzlim;         //!< Enable upper z-limit for frame tracking
  Real ft_uppzlim;             //!< Upper z-limit (shift velocities if interface above this)
  bool ft_use_lowzlim;         //!< Enable lower z-limit for frame tracking
  Real ft_lowzlim;             //!< Lower z-limit (shift velocities if interface below this)
};

// Global pointer to the TRML parameters structure.
// Allocated in ProblemGenerator::UserProblem() and persists for the simulation lifetime.
pgen_trml* ptrml = new pgen_trml();

// ====================================================================================
// FUNCTION DECLARATIONS
// ====================================================================================

//! \brief Add user-defined source terms (cooling/heating) to conserved variables.
//! Called every RK substep via user_srcs_func.
void AddUserSrcs(Mesh *pm, const Real bdt);

//! \brief Write custom history file output (mass in different phases, cooling rates, etc.).
void UserHistOutput(HistoryData *pdata, Mesh *pm);

//! \brief User-defined refinement criterion for AMR based on temperature/density gradients.
void UserRefine(MeshBlockPack* pmbp);

//! \brief Apply custom boundary conditions in the z-direction (reflecting for velocities).
void UserZBoundaryConditions(Mesh *pm);

//! \brief Compute user-defined timestep constraint (cooling CFL).
//! Limits dt so that temperature changes by at most cfl_cool per timestep.
void UserTimeStep(Mesh *pm);

//! \brief Output diagnostic information (min/max T, mass, cooling rates) every ndiag steps.
void Diagnostic(Mesh *pm);

//! \brief Apply radiative cooling and background heating to the gas.
//! Uses a power-law cooling function with parameters beta_lo, beta_hi, T_peak.
void AddCoolingHeating(Mesh *pm, const Real bdt, DvceArray5D<Real> &u0,
                   const DvceArray5D<Real> &w0, const EOS_Data &eos_data);

//! \brief Main user work-in-loop function called after each timestep.
//! Dispatches to Diagnostic() and FrameTracking() as needed.
void UserWorkInLoop(Mesh *pm);

//! \brief Apply Galilean velocity shift to keep the mixing layer interface centered.
//! Finds the mean position and velocity of gas near T_peak and shifts all velocities
//! to keep this gas near the domain center.
void FrameTracking(Mesh *pm);

//! \brief Apply a temperature ceiling to prevent runaway heating.
void ApplyTempCeiling(Mesh *pm, const Real bdt, DvceArray5D<Real> &u0,
                   const DvceArray5D<Real> &w0, const EOS_Data &eos_data);

//! \brief Adjust temperatures below T_floor by averaging with neighboring cells.
void AdjustTempTFloor(Mesh *pm, const Real bdt, DvceArray5D<Real> &u0,
                   const DvceArray5D<Real> &w0, const EOS_Data &eos_data);

//! \brief Cleanup function called at simulation end (deallocates ptrml).
void TRMLFinalWork(ParameterInput *pin, Mesh *pm);

// ====================================================================================
// MHD INITIALIZATION HELPER FUNCTIONS (Newton-Raphson solver for density profile)
// These functions solve for the density profile when B-field scales with density.
// ====================================================================================

//! \brief Density equation for Newton-Raphson solver (used when B_dens_ratio = true).
KOKKOS_INLINE_FUNCTION
Real density_equation(Real rho, Real rho_hot, Real beta_hot, Real contrast, Real alpha_magdens);

//! \brief Derivative of density equation for Newton-Raphson solver.
KOKKOS_INLINE_FUNCTION
Real density_equation_derivative(Real rho, Real rho_hot, Real beta_hot, Real contrast, Real alpha_magdens);

//! \brief Newton-Raphson iteration to solve for density given total pressure constraint.
KOKKOS_INLINE_FUNCTION
Real newton_density_solver(Real rho_guess, const Real rho_hot, const Real beta_hot, const Real contrast,
  Real alpha_magdens, Real tol, int max_iter);

//! \brief Solve for density profile at given z-position using Newton-Raphson method.
KOKKOS_INLINE_FUNCTION
Real solve_density_profile_NR(const Real z, const Real rho_hot, const Real beta_hot,
  const Real contrast, Real alpha_magdens);

} // namespace

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()
//! \brief Problem Generator for nonlinear thermal instability

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  user_srcs_func = AddUserSrcs;
  user_hist_func = UserHistOutput;
  user_bcs_func = UserZBoundaryConditions;
  user_time_step_func = UserTimeStep;
  user_work_in_loop_func = UserWorkInLoop;
  user_ref_func= UserRefine;
  pgen_final_func = TRMLFinalWork;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->phydro == nullptr && pmbp->pmhd == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
       << "TRML problem generator can only be run with Hydro and/or MHD, "
       << "but no <hydro> or <mhd> block in input file" << std::endl;
    exit(EXIT_FAILURE);
  }
  bool is_hydro = (pmbp->phydro != nullptr) ? true : false;
  bool is_mhd = (pmbp->pmhd != nullptr) ? true : false;
  EOS_Data &eos = (is_mhd) ?
                pmbp->pmhd->peos->eos_data : pmbp->phydro->peos->eos_data;
  Real lz_max = pmy_mesh_->mesh_size.x3max;
  Real lz_min = pmy_mesh_->mesh_size.x3min;

  // Read general parameters from input file
  ptrml->gamma_adi         = eos.gamma;
  ptrml->dfloor            = eos.dfloor;
  ptrml->pfloor            = eos.pfloor;
  Real number_density      = pin->GetOrAddReal("problem", "n0", 0.1);
  Real rho_0               = pin->GetReal("problem", "rho_0");
  ptrml->rho_0             = rho_0;
  Real pgas_0              = pin->GetReal("problem", "pgas_0");
  ptrml->pgas_0            = pgas_0;
  // In MHD case with scaled B-field, temperature and density contrast are different.
  // In this notation, we use contrast to denote the temperature contrast 
  // and calculate the density contrast from the temperature contrast 
  // by equating the total (magnetic+thermal) pressure in the initial state.
  Real contrast            = pin->DoesParameterExist("problem", "density_contrast") ?
                             pin->GetReal("problem", "density_contrast") : 
                             pin->GetReal("problem", "contrast");
  ptrml->contrast          = contrast; 
  Real velocity            = pin->GetReal("problem", "velocity");
  ptrml->velocity          = velocity;
  ptrml->t_shear           = 1.0/velocity;

  // Read cooling-table-related parameters from input file
  ptrml->t_cool_start      = pin->GetOrAddReal("problem", "t_cool_start", 0.0);
  ptrml->dt_cutoff         = pin->GetOrAddReal("problem", "dt_cutoff", 3.0e-5);
  ptrml->cfl_cool          = pin->GetOrAddReal("problem", "cfl_cool", 0.1);
  Real T_cold              = pgas_0/rho_0 / contrast;
  ptrml->T_cold            = T_cold;
  Real T_hot               = pgas_0/rho_0;
  ptrml->T_hot             = T_hot;
  Real T_peak              = T_cold * pin->GetReal("problem", "T_peak_over_T_cold");
  ptrml->T_peak            = T_peak;
  ptrml->T_cutoff          = T_hot * pin->GetReal("problem", "T_cutoff_over_T_hot");
  ptrml->adjust_temp_floor = pin->GetOrAddBoolean("problem", "adjust_temp_floor", false);
  ptrml->T_floor           = T_cold * pin->GetOrAddReal("problem", "T_floor", 0.1);

  ptrml->T_peak_lo         = pow(10., std::log10(T_peak) - std::log10(contrast)/20.);
  ptrml->T_peak_hi         = pow(10., std::log10(T_peak) + std::log10(contrast)/20.);
  ptrml->T_cold_hi         = pow(10., std::log10(T_cold) + 2 * std::log10(contrast)/20.);
  ptrml->T_hot_lo          = pow(10., std::log10(T_hot) - 2 * std::log10(contrast)/20.);

  Real beta_lo            = pin->GetReal("problem", "beta_lo");
  ptrml->beta_lo           = beta_lo;
  Real beta_hi            = pin->GetReal("problem", "beta_hi");
  ptrml->beta_hi           = beta_hi;
  ptrml->alpha_magdens     = pin->GetOrAddReal("problem", "alpha_magdens", 1.0/2.0);
  ptrml->xi                = pin->GetReal("problem", "xi");
  ptrml->alpha_heat        = (beta_lo-beta_hi) * (std::log(T_cold/T_peak) / std::log(contrast)) - beta_hi;
  ptrml->heat_coefficient  = pow(T_cold/T_peak, (beta_hi-beta_lo) * (std::log(T_cold/T_peak) / std::log(contrast) + 1));
  ptrml->epsilon_T         = pin->GetOrAddReal("problem", "epsilon_T", 0.05);

  // get ztop and bottom
  ptrml->ztop              = pin->GetReal("mesh", "x3max");
  ptrml->zbot              = pin->GetReal("mesh", "x3min");

  // Initial conditions and Boundary values
  Real smoothing_thickness = pin->GetReal("problem", "smoothing_thickness");
  ptrml->smoothing_thickness = smoothing_thickness;

  Real z_interface              = pin->GetOrAddReal("problem", "z_interface", 0.0);
  ptrml->z_interface       = z_interface;
  Real t_cool_0            = pin->GetOrAddReal("problem", "t_cool_0", -1.0);

  ptrml->use_temp_floor_cool = pin->GetOrAddBoolean("problem","use_temp_floor_cool", false);
  if (ptrml->use_temp_floor_cool) ptrml->T_floor_cool = pin->GetOrAddReal("problem","temp_floor_cool", 2e4)/pmbp->punit->temperature_cgs();
  else ptrml->T_floor_cool = eos.tfloor; // If user has not set a floor on the cooling, pass the eos floor

  ptrml->use_temp_ceiling = pin->GetOrAddBoolean("problem","use_temp_ceiling", false);
  if (ptrml->use_temp_ceiling) ptrml->T_ceiling = pin->GetOrAddReal("problem","temp_ceiling", 1e9)/pmbp->punit->temperature_cgs();
  ptrml->cool_subc = pin->GetOrAddBoolean("problem","cool_subc",false);
  ptrml->cool_subfac = pin->GetOrAddReal("problem","cool_subfac", 1.0);
  ptrml->use_dens_ceiling = pin->GetOrAddBoolean("problem","use_dens_ceiling_cool", false);
  if (ptrml->use_dens_ceiling){
    Real n0_ceiling = pin->GetOrAddReal("problem","n0_ceiling_cool", 20.0);
    ptrml->dens_ceiling = n0_ceiling*pmbp->punit->mu()*pmbp->punit->atomic_mass_unit_cgs/pmbp->punit->density_cgs();
  }
  else{
    ptrml->dens_ceiling = FLT_MAX;
  }
  ptrml->ndiag = pin->GetOrAddInteger("problem","ndiag",-1);

  ptrml->use_frame_tracking = pin->GetOrAddBoolean("problem","use_frame_tracking", false); // switch on/off frame tracking so that gas close to T_peak stays close to the center
  ptrml->max_vz_frame_tracking = pin->GetOrAddReal("problem","max_vz_frame_tracking", 0.01*velocity); // Maximum velocity that we impart to the domain during frame-tracking
  ptrml->t_frame_tracking_start = pin->GetOrAddReal("problem","t_frame_tracking_start", ptrml->t_cool_start); // time after which we switch on frame-tracking
  ptrml->n_frame_track = pin->GetOrAddInteger("problem","n_frame_track", 20); // perform frame-tracking every n_frame_track time steps
  // Temperature limits for frame tracking
  ptrml->T_ft_lo = T_cold * pin->GetOrAddReal("problem", "T_ft_lo_over_T_cold", 1.0); 
  ptrml->T_ft_hi = T_cold * pin->GetOrAddReal("problem", "T_ft_hi_over_T_cold", ptrml->T_peak_hi/T_cold);
  ptrml->ft_use_uppzlim   = pin->GetOrAddBoolean("problem", "ft_use_uppzlim", false);
  ptrml->ft_uppzlim       = pin->GetOrAddReal("problem", "ft_uppzlim", (z_interface+0.8*(lz_max-z_interface)));
  ptrml->ft_use_lowzlim   = pin->GetOrAddBoolean("problem", "ft_use_lowzlim", false);
  ptrml->ft_lowzlim       = pin->GetOrAddReal("problem", "ft_lowzlim", (lz_min+0.2*(z_interface-lz_min)));


  // Print info
  if (t_cool_0 > 0){
    if (global_variable::my_rank == 0) {
      std::cout << "xi was = " << ptrml->xi << "\n";
      std::cout << "setting t_cool_0 = " << t_cool_0 << "\n";
      ptrml->xi = ptrml->t_shear / t_cool_0;
      std::cout << "xi is now = " << ptrml->xi << "\n";
    }
    ptrml->xi = ptrml->t_shear / t_cool_0;
  }

  // capture variables for kernel
  auto &indcs = pmy_mesh_->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int nmb1 = (pmbp->nmb_thispack-1);

  int nx1 = indcs.nx1;
  int nx2 = indcs.nx2;
  int nx3 = indcs.nx3;

  int &ng = indcs.ng;
  int n1m1 = nx1 + 2*ng - 1;
  int n2m1 = (nx2 > 1)? (nx2 + 2*ng - 1) : 0;
  int n3m1 = (nx3 > 1)? (nx3 + 2*ng - 1) : 0;

  int64_t seed_perturb = pin->GetOrAddInteger("problem","seed_perturb",-1);
  Real sigma_perturb = pin->GetOrAddReal("problem","sigma_perturb",0.0);

  // Print info
  if (global_variable::my_rank == 0) {
    std::cout << std::setprecision(16) << "rho_0 = " << rho_0 << "\n";
    std::cout << std::setprecision(16) << "pgas_0 = " << pgas_0 << "\n";
    std::cout << std::setprecision(16) << "T_peak = " << T_peak << "\n";
    std::cout << std::setprecision(16) << "T_cold = " << T_cold << "\n";
    std::cout << std::setprecision(16) << "beta_hi = " << beta_hi << "\n";
    std::cout << std::setprecision(16) << "beta_lo = " << beta_lo << "\n";
    std::cout << std::setprecision(16) << "alpha_magdens = " << ptrml->alpha_magdens << "\n";
    std::cout << std::setprecision(16) << "xi = " << ptrml->xi << "\n";
    std::cout << std::setprecision(16) << "t_shear = " << ptrml->t_shear << "\n";
    std::cout << std::setprecision(16) << "z_interface = " << ptrml->z_interface << "\n";
    std::cout << std::setprecision(16) << "smoothing_thickness = " << ptrml->smoothing_thickness << "\n";
    std::cout << std::setprecision(16) << "contrast = " << contrast << "\n";
    std::cout << std::setprecision(16) << "velocity = " << velocity << "\n";
    std::cout << std::setprecision(16) << "T_cutoff = " << ptrml->T_cutoff << "\n";
    std::cout << std::setprecision(16) << "T_peak_lo = " << ptrml->T_peak_lo << "\n";
    std::cout << std::setprecision(16) << "T_peak_hi = " << ptrml->T_peak_hi << "\n";
    std::cout << std::setprecision(16) << "T_ft_lo = " << ptrml->T_ft_lo << "\n";
    std::cout << std::setprecision(16) << "T_ft_hi = " << ptrml->T_ft_hi << "\n";
    std::cout << std::setprecision(16) << "T_cold_hi = " << ptrml->T_cold_hi << "\n";
    std::cout << std::setprecision(16) << "T_hot_lo = " << ptrml->T_hot_lo << "\n";
    std::cout << std::setprecision(16) << "dfloor = " << ptrml->dfloor << "\n";
    std::cout << std::setprecision(16) << "pfloor = " << ptrml->pfloor << "\n";
  }
  // End print info

  auto adaptive = pmy_mesh_->adaptive;
  if (adaptive){
    ptrml->density_ratio_threshold = pin->GetOrAddReal("problem", "density_ratio_threshold", 0.0);
    ptrml->vel2_rms_threshold = pin->GetOrAddReal("problem", "vel2_rms_threshold", 0.0);
    ptrml->T_min_threshold = pin->GetOrAddReal("problem", "T_min_threshold", 0.0);
    ptrml->T_max_threshold = pin->GetOrAddReal("problem", "T_max_threshold", 0.0);
  }

  if (restart) return;

  // Initialize Hydro variables -------------------------------
  if (pmbp->phydro != nullptr || pmbp->pmhd != nullptr) {

    // Initialize Hydro/MHD variables -------------------------------
    auto &w0 = (is_mhd) ? pmbp->pmhd->w0 : pmbp->phydro->w0;
    auto &u0 = (is_mhd) ? pmbp->pmhd->u0 : pmbp->phydro->u0;
    EOS_Data &eos = (is_mhd) ?
                  pmbp->pmhd->peos->eos_data : pmbp->phydro->peos->eos_data;
    int &nscalars = (is_mhd) ?
                  pmbp->pmhd->nscalars : pmbp->phydro->nscalars;
    int &nfluid = (is_mhd) ? pmbp->pmhd->nmhd : pmbp->phydro->nhydro;
    Real gm1 = eos.gamma - 1.0;
    auto &size = pmbp->pmb->mb_size;

    // Set initial conditions
    if (global_variable::my_rank == 0) {
      std::cout << "Now initializing hydro/MHD variables." << "\n";
    }
    par_for("pgen_trml", DevExeSpace(),0,nmb1,ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real &x3min = size.d_view(m).x3min;
      Real &x3max = size.d_view(m).x3max;
      int nx3 = indcs.nx3;
      Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);

      Real rho = rho_0 * (contrast + (1.0-contrast) * 0.5 * (1.0+std::tanh((x3v-z_interface)/smoothing_thickness)));
      Real velx = velocity * 0.5 * std::tanh((x3v-z_interface)/smoothing_thickness);
      w0(m,IDN,k,j,i) = rho;
      w0(m,IVX,k,j,i) = velx;
      w0(m,IVY,k,j,i) = 0.0;
      w0(m,IVZ,k,j,i) = 0.0;
      if (eos.is_ideal) {
        w0(m,IEN,k,j,i) = pgas_0/gm1;
      }
      // add passive scalars
      if(nscalars>0){
        w0(m,nfluid,k,j,i) = 0.5 * (1.0+std::tanh((z_interface-x3v)/smoothing_thickness));
      }
    });
    // Convert primitives to conserved
    if (is_hydro) {
      pmbp->phydro->peos->PrimToCons(w0, u0, is, ie, js, je, ks, ke);
    }
    else if (is_mhd) {
      if (global_variable::my_rank == 0) {
        std::cout << "Now initializing MHD (only) variables." << "\n";
      }
      auto &b0_ = pmbp->pmhd->b0;
      auto &bcc0_ = pmbp->pmhd->bcc0;

      // First initialize fields to zero
      par_for("pgen_bondi_bfield", DevExeSpace(), 0,nmb1,0,n3m1,0,n2m1,0,n1m1,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        b0_.x1f(m,k,j,i) = 0.0;
        if (i==n1m1) b0_.x1f(m,k,j,i+1) = 0.0;
        b0_.x2f(m,k,j,i) = 0.0;
        if (j==n2m1) b0_.x2f(m,k,j+1,i) = 0.0;
        b0_.x3f(m,k,j,i) = 0.0;
        if (k==n3m1) b0_.x3f(m,k+1,j,i) = 0.0;
        bcc0_(m,IBX,k,j,i) = 0.0;
        bcc0_(m,IBY,k,j,i) = 0.0;
        bcc0_(m,IBZ,k,j,i) = 0.0;
      });

      bool tanh_mag=false, sin_mag=false;
      std::string B_z_prof = pin->GetOrAddString("problem", "B_z_prof", "uni");  // magnetic field profile
      if (global_variable::my_rank == 0) {
        std::cout << "B_z_prof = " << B_z_prof << "\n";
      }
      if(B_z_prof == "tanh") tanh_mag = true;
      else if(B_z_prof == "sin") sin_mag = true;
      Real plasma_beta = pin->GetOrAddReal("problem", "plasma_beta", 100.0);
      Real B_z_thickness = pin->GetOrAddReal("problem", "B_z_thickness", smoothing_thickness);
      bool B_use_z_boundlim = pin->GetOrAddBoolean("problem", "B_use_z_boundlim", false);
      Real B_z_uppboundary = pin->GetOrAddReal("problem", "B_z_uppboundary", lz_max-4.*smoothing_thickness);
      Real B_z_lowboundary = pin->GetOrAddReal("problem", "B_z_lowboundary", lz_min+4.*smoothing_thickness);
      int B_direction = pin->GetOrAddInteger("problem", "B_direction", 1); // 0 = x, 1 = y, 2 = z
      // parameter to scale the initial magnetic field according to the compression ratio
      bool B_dens_ratio = pin->GetOrAddBoolean("problem", "B_dens_ratio", false); 
      Real alpha_magdens = ptrml->alpha_magdens;  // power-law index for B-field scaling with density profile in cold gas
      
      // Print MHD initialization parameters
      if (global_variable::my_rank == 0) {
        std::cout << std::setprecision(16);
        std::cout << "MHD Initialization Parameters:\n";
        std::cout << "  plasma_beta = " << plasma_beta << "\n";
        std::cout << "  B_direction = " << B_direction << " (0=x, 1=y, 2=z)\n";
        std::cout << "  B_z_thickness = " << B_z_thickness << "\n";
        std::cout << "  B_use_z_boundlim = " << (B_use_z_boundlim ? "true" : "false") << "\n";
        if (B_use_z_boundlim) {
          std::cout << "  B_z_lowboundary = " << B_z_lowboundary << "\n";
          std::cout << "  B_z_uppboundary = " << B_z_uppboundary << "\n";
        }
        std::cout << "  B_dens_ratio = " << (B_dens_ratio ? "true" : "false") << "\n";
        std::cout << "  alpha_magdens = " << alpha_magdens << "\n";
      }
      
      par_for("pgen_set_bfield", DevExeSpace(),0,nmb1,ks,ke,js,je,is,ie,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        Real &x3min = size.d_view(m).x3min;
        Real &x3max = size.d_view(m).x3max;
        int nx3 = indcs.nx3;
        Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
        Real z_weight = 1.0;
        if(tanh_mag){
          z_weight = std::tanh((x3v-z_interface)/B_z_thickness);
        }
        else if(sin_mag){
          z_weight = std::sin(M_PI*(x3v-z_interface)/B_z_thickness);
        }
        else{
          z_weight = 1.0;
        }
        if(B_use_z_boundlim){
          // Smoothly turn off magnetic field at the boundaries, use an exponential function
          if(x3v > B_z_uppboundary) z_weight *= std::exp(-2.*(x3v-B_z_uppboundary)/B_z_thickness);
          if(x3v < B_z_lowboundary) z_weight *= std::exp(-2.*(B_z_lowboundary-x3v)/B_z_thickness);
        }
        Real pmag_z = pgas_0/plasma_beta;

        if(B_dens_ratio){
          // Scale the magnetic field according to the compression ratio
          Real temp_0 = pgas_0/rho_0;
          Real ptot_0 = pgas_0*(1.0+1.0/plasma_beta);
          Real contrast_z = (contrast + (1.0-contrast) * 0.5 * (1.0+std::tanh((x3v-z_interface)/smoothing_thickness)));
          Real temp_z = temp_0/contrast_z;
          Real rho_z = solve_density_profile_NR(x3v, rho_0, plasma_beta, contrast_z, alpha_magdens);
          Real pgas_z = rho_z*temp_z;
          pmag_z = fmax(ptot_0 - pgas_z, 0.0);
          // Now update the density and pressure to be consistent with the magnetic field
          w0(m,IDN,k,j,i) = rho_z;
          if (eos.is_ideal) {
            w0(m,IEN,k,j,i) = pgas_z/gm1;
          }
        }
        b0_.x1f(m,k,j,i) = (B_direction == 0) ? z_weight*sqrt(2*pmag_z): 0.0;
        b0_.x2f(m,k,j,i) = (B_direction == 1) ? z_weight*sqrt(2*pmag_z): 0.0;
        b0_.x3f(m,k,j,i) = (B_direction == 2) ? z_weight*sqrt(2*pmag_z): 0.0;

        if (i==ie) b0_.x1f(m,k,j,i+1) = (B_direction == 0) ? z_weight*sqrt(2*pmag_z): 0.0;
        if (j==je) b0_.x2f(m,k,j+1,i) = (B_direction == 1) ? z_weight*sqrt(2*pmag_z): 0.0;
        if (k==ke) b0_.x3f(m,k+1,j,i) = (B_direction == 2) ? z_weight*sqrt(2*pmag_z): 0.0;
      });
      par_for("pgen_set_bfield", DevExeSpace(),0,nmb1,ks,ke,js,je,is,ie,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        bcc0_(m,IBX,k,j,i) = 0.5*(b0_.x1f(m,k,j,i) + b0_.x1f(m,k  ,j  ,i+1));
        bcc0_(m,IBY,k,j,i) = 0.5*(b0_.x2f(m,k,j,i) + b0_.x2f(m,k  ,j+1,i  ));
        bcc0_(m,IBZ,k,j,i) = 0.5*(b0_.x3f(m,k,j,i) + b0_.x3f(m,k+1,j  ,i  ));
      });
      // Convert primitives to conserved
      pmbp->pmhd->peos->PrimToCons(w0, bcc0_, u0, 0, n1m1, 0, n2m1, 0, n3m1);
    }
    if (global_variable::my_rank == 0) {
      std::cout << "Hydro/MHD variables initialization done" << "\n";
    }

    // Initialize particles - within Hydro/MHD block check
    if (pmbp->ppart != nullptr) {
      // captures for the kernel
      auto &mblev = pmbp->pmb->mb_lev;
      auto gids = pmbp->gids;

      // Get fixed seed for reproducibility - will be offset by MeshBlock ID
      int64_t pos_init_seed = pin->GetOrAddInteger("particles","pos_init_seed",280496);
      
      // Get distribution type: false = by mass (default), true = uniform by volume
      bool uniform_by_volume = pin->GetOrAddBoolean("particles","uniform_by_volume",true);

      // count total mass/volume across the domain
      Real total_mass = 0.0;
      Real total_volume = 0.0;
      const int nmkji = (pmbp->nmb_thispack)*indcs.nx3*indcs.nx2*indcs.nx1;
      const int nkji = indcs.nx3*indcs.nx2*indcs.nx1;
      const int nji  = indcs.nx2*indcs.nx1;

      Kokkos::parallel_reduce("pgen_mass_vol", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
      KOKKOS_LAMBDA(const int &idx, Real &total_mass, Real &total_volume) {
        // compute m,k,j,i indices of thread and evaluate
        int m = (idx)/nkji;
        int k = (idx - m*nkji)/nji;
        int j = (idx - m*nkji - k*nji)/indcs.nx1;
        int i = (idx - m*nkji - k*nji - j*indcs.nx1) + is;
        k += ks;
        j += js;

        Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
        total_mass += u0(m,IDN,k,j,i) * vol;
        total_volume += vol;
      }, total_mass, total_volume);

      Real total_mass_thispack = total_mass;
      Real total_volume_thispack = total_volume;

  #if MPI_PARALLEL_ENABLED
      // get total mass/volume over all MPI ranks
      MPI_Allreduce(MPI_IN_PLACE, &total_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &total_volume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #endif

      // get number of particles for this mbpack using MC to deal with fractional particles
      Real target_nparticles = pin->GetOrAddReal("particles","target_count",100000.0);
      Real mass_per_particle = total_mass / target_nparticles;
      Real volume_per_particle = total_volume / target_nparticles;

      // create shared array to hold number of particles per zone
      DualArray2D<int> nparticles_per_zone("partperzone", nmkji,2);
      par_for("particle_count", DevExeSpace(), 0,nmkji-1,
      KOKKOS_LAMBDA(int idx) {
        int m = (idx)/nkji;
        int k = (idx - m*nkji)/nji;
        int j = (idx - m*nkji - k*nji)/indcs.nx1;
        int i = (idx - m*nkji - k*nji - j*indcs.nx1) + is;
        k += ks;
        j += js;
        
        // Create deterministic pseudo-random number using hash function
        // This ensures same particle count for same zone regardless of MPI decomposition
        int64_t mb_gid = gids + m;  // Global MeshBlock ID
        int64_t zone_seed = pos_init_seed + mb_gid * 999983 + i * 7919 + j * 104729 + k * 524287;
        
        // Simple hash-based pseudo-random number generation (device-compatible)
        // Using a variant of the splitmix64 algorithm
        uint64_t z = static_cast<uint64_t>(zone_seed);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        z = z ^ (z >> 31);
        Real rand_val = static_cast<Real>(z & 0x7FFFFFFFULL) / static_cast<Real>(0x80000000ULL);
        
        Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
        
        // Calculate number of particles per zone based on distribution type
        Real nppc;
        if (uniform_by_volume) {
          // Uniform distribution: each zone gets same number per unit volume
          nppc = vol / volume_per_particle;
        } else {
          // Mass-weighted distribution: zones get particles proportional to mass
          nppc = u0(m,IDN,k,j,i) * vol / mass_per_particle;
        }
        
        int nparticles = static_cast<int>(nppc);
        nppc = fabs(fmod(nppc, 1.0));
        if (rand_val < nppc) {
          nparticles += 1;
        }
        nparticles_per_zone.d_view(idx,0) = nparticles;
      });

      // count total number of particles in this pack
      nparticles_per_zone.template modify<DevExeSpace>();
      nparticles_per_zone.template sync<HostMemSpace>();
      int nparticles_thispack = 0;
      for (int i=0; i<nmkji; ++i) {
        nparticles_per_zone.h_view(i,1) = nparticles_thispack;
        nparticles_thispack += nparticles_per_zone.h_view(i,0);
      }
      nparticles_per_zone.template modify<HostMemSpace>();
      nparticles_per_zone.template sync<DevMemSpace>();

      // helpful debug statement
      if (uniform_by_volume) {
        if (global_variable::my_rank == 0) {
          std::cout << "Particle distribution: UNIFORM BY VOLUME" << std::endl;
        }
        std::cout << "Rank " << global_variable::my_rank
                  << ": total volume across domain: " << total_volume
                  << ", total volume in pack: " << total_volume_thispack
                  << ", target nparticles: " << target_nparticles
                  << ", nparticles in pack: " << nparticles_thispack
                  << std::endl;
      } else {
        if (global_variable::my_rank == 0) {
          std::cout << "Particle distribution: BY MASS (default)" << std::endl;
        }
        std::cout << "Rank " << global_variable::my_rank
                  << ": total mass across domain: " << total_mass
                  << ", total mass in pack: " << total_mass_thispack
                  << ", target nparticles: " << target_nparticles
                  << ", nparticles in pack: " << nparticles_thispack
                  << std::endl;
      }

      // reallocate space for particles and get relevant pointers
      pmbp->ppart->ReallocateParticles(nparticles_thispack);

      auto &pr = pmbp->ppart->prtcl_rdata;
      auto &pi = pmbp->ppart->prtcl_idata;

      // initialize particles. only intended for Lagrangian-type particles
      par_for("part_init", DevExeSpace(), 0,nmkji-1,
      KOKKOS_LAMBDA(int idx) {
        int m = (idx)/nkji;
        int k = (idx - m*nkji)/nji;
        int j = (idx - m*nkji - k*nji)/indcs.nx1;
        int i = (idx - m*nkji - k*nji - j*indcs.nx1) + is;
        k += ks;
        j += js;

        int nparticles_in_zone = nparticles_per_zone.d_view(idx,0);
        int starting_index = nparticles_per_zone.d_view(idx,1);

        for (int p=0; p<nparticles_in_zone; ++p) {
          int pidx = p + starting_index;

          pi(PGID,pidx) = gids + m;
          pi(PLASTLEVEL,pidx) = mblev.d_view(m);

          // set particle to zone center
          pr(IPX,pidx) = CellCenterX(i-is, nx1, size.d_view(m).x1min,
                                    size.d_view(m).x1max);
          pr(IPY,pidx) = CellCenterX(j-js, nx2, size.d_view(m).x2min,
                                    size.d_view(m).x2max);
          pr(IPZ,pidx) = CellCenterX(k-ks, nx3, size.d_view(m).x3min,
                                    size.d_view(m).x3max) -
                        size.d_view(m).dx3/2;
        }
      });
    }
  }
  return;
}

namespace {

//----------------------------------------------------------------------------------------
//! \fn UserZBoundaryConditions
//  \brief Sets boundary condition on surfaces of computational domain
// For this problem we shall set density and pressure to be constant along outer z boundaries
// and maintain a constant v_x
// At lower boundary we set no inflow
void UserZBoundaryConditions(Mesh *pm) {
  auto &indcs = pm->mb_indcs;
  MeshBlockPack *pmbp = pm->pmb_pack;
  int &ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng) : 1;

  auto &mb_bcs = pmbp->pmb->mb_bcs;
  bool is_hydro = (pmbp->phydro != nullptr) ? true : false;
  bool is_mhd = (pmbp->pmhd != nullptr) ? true : false;

  Real rho_0 = ptrml->rho_0;
  Real pgas_0 = ptrml->pgas_0;
  Real velocity = ptrml->velocity;
  // Initialize Hydro/MHD variables -------------------------------
  int nmb = pmbp->nmb_thispack;
  auto u0_ = (is_mhd) ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  auto w0_ = (is_mhd) ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  EOS_Data &eos = (is_mhd) ?
                  pmbp->pmhd->peos->eos_data : pmbp->phydro->peos->eos_data;
  int nvar = u0_.extent_int(1);

  int ks = indcs.ks, ke = indcs.ke;

  auto gm1 = eos.gamma - 1.0;
  Real use_e = eos.use_e;

  // x3-Boundary
  // ConsToPrim over all x3 ghost zones *and* at the innermost/outermost x3-active zones
  // of Meshblocks, even if Meshblock face is not at the edge of computational domain
  if (is_hydro) {
    pmbp->phydro->peos->ConsToPrim(u0_,w0_,false,0,(n1-1),0,(n2-1),ks-ng,ks);
    pmbp->phydro->peos->ConsToPrim(u0_,w0_,false,0,(n1-1),0,(n2-1),ke,ke+ng);
  } else if (is_mhd) {
    auto &b0 = pmbp->pmhd->b0;
    auto &bcc = pmbp->pmhd->bcc0;
    pmbp->pmhd->peos->ConsToPrim(u0_,b0,w0_,bcc,false,0,(n1-1),0,(n2-1),ks-ng,ks);
    pmbp->pmhd->peos->ConsToPrim(u0_,b0,w0_,bcc,false,0,(n1-1),0,(n2-1),ke,ke+ng);
  }
  else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
        << "Should be either hydro or MHD setup" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  par_for("userBC_x3", DevExeSpace(),0,(nmb-1),0,(nvar-1),0,(n2-1),0,(n1-1),
  KOKKOS_LAMBDA(int m, int n, int j, int i) {
    if (mb_bcs.d_view(m,BoundaryFace::inner_x3) == BoundaryFlag::user) {
      for (int k=0; k<ng; ++k) {
        // noinflow at inner x3 BC
        // if (n==(IVZ)){
        //   w0_(m,n,ks-k-1,j,i) = 0.0;
        // }
        // else {
        w0_(m,n,ks-k-1,j,i) = w0_(m,n,ks,j,i);
        // }
      }
    }
    if (mb_bcs.d_view(m,BoundaryFace::outer_x3) == BoundaryFlag::user) {
      for (int k=0; k<ng; ++k) {
        if (n == (IDN)){
          w0_(m,n,ke+k+1,j,i) = rho_0;
        }
        else if (n == (IEN) || (n==ITM)){
          if (use_e){
            w0_(m,n,ke+k+1,j,i) = pgas_0 /gm1 ;
          }
          else {
            w0_(m,n,ke+k+1,j,i) = pgas_0 /rho_0 ;
          }
        }
        else if (n == (IVX)){
          w0_(m,n,ke+k+1,j,i) = 0.5*velocity;
        }
        else {
          w0_(m,n,ke+k+1,j,i) = w0_(m,n,ke,j,i);
        }
      }
    }
  });

  // Now B-fields
  // Set x3-BCs on b0 and bcc0 if Meshblock face is at the edge of computational domain
  if (is_mhd) {
    auto &b0 = pmbp->pmhd->b0;
    auto &bcc_ = pmbp->pmhd->bcc0;
    par_for("noinflow_field_x3", DevExeSpace(),0,(nmb-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int j, int i) {
      if (mb_bcs.d_view(m,BoundaryFace::inner_x3) == BoundaryFlag::user) {
        for (int k=0; k<ng; ++k) {
          b0.x1f(m,ks-k-1,j,i) = b0.x1f(m,ks,j,i);
          if (i == n1-1) {b0.x1f(m,ks-k-1,j,i+1) = b0.x1f(m,ks,j,i+1);}
          b0.x2f(m,ks-k-1,j,i) = b0.x2f(m,ks,j,i);
          if (j == n2-1) {b0.x2f(m,ks-k-1,j+1,i) = b0.x2f(m,ks,j+1,i);}
          b0.x3f(m,ks-k-1,j,i) = b0.x3f(m,ks,j,i);
        }
      }
      if (mb_bcs.d_view(m,BoundaryFace::outer_x3) == BoundaryFlag::user) {
        for (int k=0; k<ng; ++k) {
          b0.x1f(m,ke+k+1,j,i) = b0.x1f(m,ke,j,i);
          if (i == n1-1) {b0.x1f(m,ke+k+1,j,i+1) = b0.x1f(m,ke,j,i+1);}
          b0.x2f(m,ke+k+1,j,i) = b0.x2f(m,ke,j,i);
          if (j == n2-1) {b0.x2f(m,ke+k+1,j+1,i) = b0.x2f(m,ke,j+1,i);}
          b0.x3f(m,ke+k+2,j,i) = b0.x3f(m,ke+1,j,i);
        }
      }
    });
    par_for("noinflow_field_x3", DevExeSpace(),0,(nmb-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int j, int i) {
      if (mb_bcs.d_view(m,BoundaryFace::inner_x3) == BoundaryFlag::user) {
        for (int k=0; k<ng; ++k) {
          bcc_(m,IBX,ks-k-1,j,i) = 0.5*(b0.x1f(m,ks-k-1,j,i) + b0.x1f(m,ks-k-1,j  ,i+1));
          bcc_(m,IBY,ks-k-1,j,i) = 0.5*(b0.x2f(m,ks-k-1,j,i) + b0.x2f(m,ks-k-1,j+1,i  ));
          bcc_(m,IBZ,ks-k-1,j,i) = 0.5*(b0.x3f(m,ks-k-1,j,i) + b0.x3f(m,ks-k  ,j  ,i  ));
        }
      }
      if (mb_bcs.d_view(m,BoundaryFace::outer_x3) == BoundaryFlag::user) {
        for (int k=0; k<ng; ++k) {
          bcc_(m,IBX,ke+k+1,j,i) = 0.5*(b0.x1f(m,ke+k+1,j,i) + b0.x1f(m,ke+k+1,j  ,i+1));
          bcc_(m,IBY,ke+k+1,j,i) = 0.5*(b0.x2f(m,ke+k+1,j,i) + b0.x2f(m,ke+k+1,j+1,i  ));
          bcc_(m,IBZ,ke+k+1,j,i) = 0.5*(b0.x3f(m,ke+k+1,j,i) + b0.x3f(m,ke+k+2,j  ,i  ));
        }
      }
    });
  }
  // PrimToCons on x3 ghost zones 
  if (is_hydro) {
    pmbp->phydro->peos->PrimToCons(w0_,u0_,0,(n1-1),0,(n2-1),ks-ng,ks-1);
    pmbp->phydro->peos->PrimToCons(w0_,u0_,0,(n1-1),0,(n2-1),ke+1,ke+ng);
  } else if (is_mhd) {
    auto &bcc0_ = pmbp->pmhd->bcc0;
    pmbp->pmhd->peos->PrimToCons(w0_,bcc0_,u0_,0,(n1-1),0,(n2-1),ks-ng,ks-1);
    pmbp->pmhd->peos->PrimToCons(w0_,bcc0_,u0_,0,(n1-1),0,(n2-1),ke+1,ke+ng);
  }

  // Now implement frame tracking for the mixing layer
  // Perform a Galilean shift such that the interface stays close to the center and does not interact with the boundaries.
  // int nmb1 = pmbp->nmb_thispack - 1;

  // if (ptrml->use_frame_tracking && pm->time > ptrml->t_frame_tracking_start && pm->ncycle % ptrml->n_frame_track == 0) {

  //   // if (is_hydro) {
  //   //   pmbp->phydro->peos->ConsToPrim(u0_,w0_,false,0,(n1-1),0,(n2-1),0,(n3-1));
  //   // } else if (is_mhd) {
  //   //   auto &b0 = pmbp->pmhd->b0;
  //   //   auto &bcc = pmbp->pmhd->bcc0;
  //   //   pmbp->pmhd->peos->ConsToPrim(u0_,b0,w0_,bcc,false,0,(n1-1),0,(n2-1),0,(n3-1));
  //   // }
  //   Real T_peak = ptrml->T_peak;
  //   Real T_peak_lo = ptrml->T_peak_lo;
  //   Real T_peak_hi = ptrml->T_peak_hi;

  //   Real mass_peak = 0.0;
  //   Real mean_z_peak = 0.0;
  //   Real mean_velz_peak = 0.0;
  //   Real min_z_peak = pm->mesh_size.x3max;
  //   Real max_z_peak = pm->mesh_size.x3min;

  //   Real max_z_frame_tracking = ptrml->max_z_frame_tracking; // Maximum allowed z-value of gas near T_peak
  //   Real min_z_frame_tracking = ptrml->min_z_frame_tracking; // Minimum allowed z-value of gas near T_peak
  //   Real mean_z_frame_tracking = ptrml->z_interface; // We'll try to keep the interface near the initial value
  //   Real max_vz_frame_tracking = ptrml->max_vz_frame_tracking; // Maximum velocity that we impart to the domain during frame-tracking

  //   Real x3min_mesh = pm->mesh_size.x3min;
  //   Real x3max_mesh = pm->mesh_size.x3max;

  //   Kokkos::parallel_reduce("frametracking", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  //   KOKKOS_LAMBDA(const int &idx, Real &masspeak_, Real &mean_z_peak_, Real &mean_velz_peak_, Real &min_z_peak_, Real &max_z_peak_) {
  //     int m = (idx)/nkji;
  //     int k = (idx - m*nkji)/nji;
  //     int j = (idx - m*nkji - k*nji)/nx1;
  //     int i = (idx - m*nkji - k*nji - j*nx1) + is;
  //     k += ks;
  //     j += js;
  //     Real dens=1.0, temp = 1.0, eint = 1.0;
  //     dens = w0_(m,IDN,k,j,i);
  //     Real dvol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
  //     if (use_e) {
  //       temp = w0_(m,IEN,k,j,i)/dens*gm1;
  //       eint = w0_(m,IEN,k,j,i);
  //     } else {
  //       temp = w0_(m,ITM,k,j,i);
  //       eint = w0_(m,ITM,k,j,i)*dens/gm1;
  //     }
  //     Real pres = eint*gm1;

  //     Real &x3min = size.d_view(m).x3min;
  //     Real &x3max = size.d_view(m).x3max;
  //     int nx3 = indcs.nx3;
  //     Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);

  //     if ((temp<T_peak_hi) && (temp>T_peak_lo)){
  //       masspeak_ += dens*dvol;
  //       mean_z_peak_ += x3v*dens*dvol;
  //       mean_velz_peak_ += w0_(m,IVZ,k,j,i)*dens*dvol;
  //       min_z_peak_ = fmin(x3v,min_z_peak_);
  //       max_z_peak_ = fmax(x3v,max_z_peak_);
  //     }
      
  //   }, Kokkos::Sum<Real>(mass_peak), 
  //     Kokkos::Sum<Real>(mean_z_peak),
  //     Kokkos::Sum<Real>(mean_velz_peak),
  //     Kokkos::Min<Real>(min_z_peak),
  //     Kokkos::Max<Real>(max_z_peak));
  // #if MPI_PARALLEL_ENABLED
  //   MPI_Allreduce(MPI_IN_PLACE, &mass_peak, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  //   MPI_Allreduce(MPI_IN_PLACE, &mean_z_peak, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  //   MPI_Allreduce(MPI_IN_PLACE, &mean_velz_peak, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  //   MPI_Allreduce(MPI_IN_PLACE, &min_z_peak, 1, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
  //   MPI_Allreduce(MPI_IN_PLACE, &max_z_peak, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
  // #endif

  //   if (global_variable::my_rank == 0 && ptrml->ndiag>0 && pm->ncycle % ptrml->ndiag == 0) {
  //     std::cout << " mass_peak=" << mass_peak << std::endl;
  //   }
  //   if(mass_peak > 0.){
  //     mean_z_peak = mean_z_peak/mass_peak;
  //     mean_velz_peak = mean_velz_peak/mass_peak;

  //     // First we try to maintain the mean z near the intial value
  //     // Real z_shift = mean_z_frame_tracking - mean_z_peak;

  //     // First we try to maintain the mean velz of gas at T_peak close to zero
  //     Real vel_z_shift = -mean_velz_peak;
  //     Real z_shift = vel_z_shift*pm->dt;

  //     // If some gas goes too close to the upper boundary we apply a negative shift
  //     if((max_z_peak+z_shift) > max_z_frame_tracking) z_shift = max_z_frame_tracking - max_z_peak;
  //     // If some gas goes too close to the bottom boundary we add a positive shift
  //     if((min_z_peak+z_shift) < min_z_frame_tracking) z_shift = min_z_frame_tracking - min_z_peak; 

  //     // Now update vel_z_shift after applying these conditions 
  //     vel_z_shift = z_shift/pm->dt;
  //     // We set a maximum velocity that the frame-tracking can impart
  //     if(fabs(vel_z_shift) > fabs(max_vz_frame_tracking)) vel_z_shift = vel_z_shift * fabs(max_vz_frame_tracking/vel_z_shift);

  //     if (global_variable::my_rank == 0) {
  //       std::cout << " mean_z_peak=" << mean_z_peak << std::endl;
  //       std::cout << " mean_velz_peak=" << mean_velz_peak << std::endl;
  //       std::cout << " min_z_peak=" << min_z_peak << std::endl;
  //       std::cout << " max_z_peak=" << max_z_peak << std::endl;
  //       std::cout << " vel_z_shift=" << vel_z_shift << std::endl;
  //     }
  //     // Now perform the frame shift ()
  //     par_for("velzframetracking", DevExeSpace(), 0, nmb1, ks-ng, ke+ng, js-ng, je+ng, is-ng, ie+ng,
  //     KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
  //       Real dmomx3 = w0_(m,IDN,k,j,i)*vel_z_shift;
  //       Real kin_energy_old = 0.5*w0_(m,IDN,k,j,i)*(SQR(w0_(m,IVX,k,j,i))+SQR(w0_(m,IVY,k,j,i))+SQR(w0_(m,IVZ,k,j,i)));
  //       Real kin_energy_new = 0.5*w0_(m,IDN,k,j,i)*(SQR(w0_(m,IVX,k,j,i))+SQR(w0_(m,IVY,k,j,i))+SQR(w0_(m,IVZ,k,j,i) + vel_z_shift));
  //       Real denergy = kin_energy_new - kin_energy_old;
  //       w0_(m,IVZ,k,j,i) += vel_z_shift;
  //       u0_(m,IM3,k,j,i) += dmomx3;
  //       u0_(m,IEN,k,j,i) += denergy;
  //     });
  //   }
  //   // Now PrimToCons 
  //   // if (is_hydro) {
  //   //   pmbp->phydro->peos->PrimToCons(w0_,u0_,0,(n1-1),0,(n2-1),0,(n3-1));
  //   // } else if (is_mhd) {
  //   //   auto &bcc0_ = pmbp->pmhd->bcc0;
  //   //   pmbp->pmhd->peos->PrimToCons(w0_,bcc0_,u0_,0,(n1-1),0,(n2-1),0,(n3-1));
  //   // }
  // }

}

//----------------------------------------------------------------------------------------
// ! \fn UserTimeStep
// ! \brief Sets the time step for cooling
// ----------------------------------------------------------------------------------------

void UserTimeStep(Mesh *pm){
  auto &indcs = pm->mb_indcs;
  int is = indcs.is, nx1 = indcs.nx1;
  int js = indcs.js, nx2 = indcs.nx2;
  int ks = indcs.ks, nx3 = indcs.nx3;
  MeshBlockPack *pmbp = pm->pmb_pack;
  const int nmkji = (pmbp->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  Real dtnew_loc = static_cast<Real>(std::numeric_limits<float>::max());

  bool is_mhd = (pmbp->pmhd != nullptr) ? true : false;

  auto u0 = (is_mhd) ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  auto w0 = (is_mhd) ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  EOS_Data &eos = (is_mhd) ?
                  pmbp->pmhd->peos->eos_data : pmbp->phydro->peos->eos_data;
  Real use_e = eos.use_e;
  Real gamma = eos.gamma;
  Real gm1 = gamma - 1.0;

  Real T_cutoff = ptrml->T_cutoff;
  Real T_hot = ptrml->T_hot;
  Real T_peak = ptrml->T_peak;
  Real epsilon_T = ptrml->epsilon_T;
  Real xi = ptrml->xi;
  Real pgas_0 = ptrml->pgas_0;
  Real t_shear = ptrml->t_shear;
  Real beta_lo = ptrml->beta_lo;
  Real beta_hi = ptrml->beta_hi;
  Real heat_coefficient = ptrml->heat_coefficient;
  Real alpha_heat = ptrml->alpha_heat;
  Real dt_cutoff = ptrml->dt_cutoff;
  Real cfl_cool = ptrml->cfl_cool;

  // For applying density ceiling
  bool use_dens_ceiling = ptrml->use_dens_ceiling;
  Real dens_ceiling = ptrml->dens_ceiling;

  // find smallest (e/cooling_rate) in each cell
  Kokkos::parallel_reduce("srcterms_cooling_newdt",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, Real &min_dt) {
    // compute m,k,j,i indices of thread and call function
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;

    // temperature in cgs unit
    Real temp = 1.0;
    Real eint = 1.0;
    if (use_e) {
      temp = w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1;
      eint = w0(m,IEN,k,j,i);
    } else {
      temp = w0(m,ITM,k,j,i);
      eint = w0(m,ITM,k,j,i)*w0(m,IDN,k,j,i)/gm1;
    }

    Real pres = eint*gm1;

    Real cooling_heating = FLT_MIN;

    if (temp>T_cutoff){
      cooling_heating = FLT_MIN;
    } else if (use_dens_ceiling && w0(m,IDN,k,j,i) > dens_ceiling) {
      // If density is above ceiling, we set cooling/heating to zero
      cooling_heating = FLT_MIN;
    } else {
      Real Edot_cool = 1.5 * xi * pgas_0/t_shear * SQR(pres/pgas_0);
      Edot_cool *= (temp<T_peak) ? pow(temp/T_peak, -beta_lo):pow(temp/T_peak, -beta_hi);

      Real Edot_heat = 1.5 * xi * pgas_0/t_shear * heat_coefficient * SQR(pres/pgas_0);
      Edot_heat *= (temp<(1+epsilon_T)*T_hot) ? pow((temp/T_peak),alpha_heat) : pow(((1+epsilon_T)*T_hot/T_peak),alpha_heat) * pow( (temp/((1+epsilon_T)*T_hot)),(-beta_hi-0.5));

      cooling_heating=FLT_MIN+Edot_cool-Edot_heat;
    }
    Real dt_cool = cfl_cool * (eint/fabs(cooling_heating));
    dt_cool = fmax(dt_cool,dt_cutoff);
    min_dt = fmin(min_dt,dt_cool);
  }, Kokkos::Min<Real>(dtnew_loc));
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, &dtnew_loc, 1, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
#endif
  pm->pgen->dtnew = dtnew_loc;

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void AddUserSrcs()
//! \brief Add User Source Terms
// NOTE source terms must all be computed using primitive (w0) and NOT conserved (u0) vars
void AddUserSrcs(Mesh *pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  bool is_mhd = (pmbp->pmhd != nullptr) ? true : false;
  DvceArray5D<Real> &u0 = (is_mhd) ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  const DvceArray5D<Real> &w0 = (is_mhd) ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  const EOS_Data &eos_data = (is_mhd) ?
                  pmbp->pmhd->peos->eos_data : pmbp->phydro->peos->eos_data;
  if (ptrml->adjust_temp_floor) AdjustTempTFloor(pm,bdt,u0,w0,eos_data);
  if (pm->time > ptrml->t_cool_start) {
    AddCoolingHeating(pm,bdt,u0,w0,eos_data);
  }
  // if (ptrml->use_frame_tracking && pm->time > ptrml->t_frame_tracking_start && pm->ncycle % ptrml->n_frame_track == 0) FrameTracking(pm,bdt,u0,w0,eos_data);
  if (ptrml->use_temp_ceiling) ApplyTempCeiling(pm,bdt,u0,w0,eos_data);
  return;
}


// ----------------------------------------------------------------------------------------
// ! \fn void SourceTerms::AddCoolingHeating()
// ! \brief Add explict ISM cooling and heating source terms in the energy equations.
// NOTE source terms must all be computed using primitive (w0) and NOT conserved (u0) vars
void AddCoolingHeating(Mesh *pm, const Real bdt, DvceArray5D<Real> &u0,
                   const DvceArray5D<Real> &w0, const EOS_Data &eos_data) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  int is = indcs.is, nx1 = indcs.nx1;
  int js = indcs.js, nx2 = indcs.nx2;
  int ks = indcs.ks, nx3 = indcs.nx3;
  const int nmkji = (pmbp->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  auto &size = pmbp->pmb->mb_size;
  Real bta_time = bdt/pm->dt;
  Real use_e = eos_data.use_e;
  Real gamma = eos_data.gamma;
  Real gm1 = gamma - 1.0;

  Real T_cutoff = ptrml->T_cutoff;
  Real T_hot = ptrml->T_hot;
  Real T_peak = ptrml->T_peak;
  Real epsilon_T = ptrml->epsilon_T;
  Real xi = ptrml->xi;
  Real pgas_0 = ptrml->pgas_0;
  Real t_shear = ptrml->t_shear;
  Real beta_lo = ptrml->beta_lo;
  Real beta_hi = ptrml->beta_hi;
  Real heat_coefficient = ptrml->heat_coefficient;
  Real alpha_heat = ptrml->alpha_heat;

  // For applying density ceiling
  bool use_dens_ceiling = ptrml->use_dens_ceiling;
  Real dens_ceiling = ptrml->dens_ceiling;

  Real net_cool=0.0, net_heat=0.0;
  Kokkos::parallel_reduce("coolingheating", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, Real &net_cooling, Real &net_heating) {
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;
    Real dens=1.0, temp = 1.0, eint = 1.0;
    dens = w0(m,IDN,k,j,i);
    Real dvol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
    if (use_e) {
      temp = w0(m,IEN,k,j,i)/dens*gm1;
      eint = w0(m,IEN,k,j,i);
    } else {
      temp = w0(m,ITM,k,j,i);
      eint = w0(m,ITM,k,j,i)*dens/gm1;
    }
    Real pres = eint*gm1;

    Real cooling_heating=0.0;

    if (temp>T_cutoff){
      cooling_heating=  FLT_MIN; // No cooling/heating if temperature exceeds cutoff
    } else if (use_dens_ceiling && dens > dens_ceiling) {
      cooling_heating = FLT_MIN; // No cooling/heating if density exceeds ceiling
    } else {
      Real Edot_cool = 1.5 * xi * pgas_0/t_shear * SQR(pres/pgas_0);
      Edot_cool *= (temp<T_peak) ? pow(temp/T_peak, -beta_lo):pow(temp/T_peak, -beta_hi);
      net_cooling += Edot_cool*dvol;

      Real Edot_heat = 1.5 * xi * pgas_0/t_shear * heat_coefficient * SQR(pres/pgas_0);
      Edot_heat *= (temp<(1+epsilon_T)*T_hot) ? pow((temp/T_peak),alpha_heat) : pow(((1+epsilon_T)*T_hot/T_peak),alpha_heat) * pow( (temp/((1+epsilon_T)*T_hot)),(-beta_hi-0.5));
      net_heating += Edot_heat*dvol;

      cooling_heating=Edot_cool-Edot_heat;
    }
    u0(m,IEN,k,j,i) -= bdt * cooling_heating;

  }, Kokkos::Sum<Real>(net_cool), Kokkos::Sum<Real>(net_heat));
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &net_cool, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &net_heat, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (global_variable::my_rank == 0) {
    if (ptrml->ndiag>0 && pm->ncycle % ptrml->ndiag == 0 && bta_time > 0.99) {
        std::cout << " net_cool=" << net_cool << std::endl
                  << " net_heat=" << net_heat << std::endl;
    }
  }
  ptrml->tot_coolrate = net_cool-net_heat; // Update the global variable
  return;
}


// ----------------------------------------------------------------------------------------
// ! \fn void SourceTerms::ApplyTempCeiling()
// ! \brief Set a temperature ceiling as specified in the input file.
// NOTE source terms must all be computed using primitive (w0) and NOT conserved (u0) vars
void ApplyTempCeiling(Mesh *pm, const Real bdt, DvceArray5D<Real> &u0,
                   const DvceArray5D<Real> &w0, const EOS_Data &eos_data) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  int is = indcs.is, nx1 = indcs.nx1;
  int js = indcs.js, nx2 = indcs.nx2;
  int ks = indcs.ks, nx3 = indcs.nx3;
  const int nmkji = (pmbp->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  Real bta_time = bdt/pm->dt;
  Real use_e = eos_data.use_e;
  Real gamma = eos_data.gamma;
  Real gm1 = gamma - 1.0;
  Real temp_unit = pmbp->punit->temperature_cgs();
  int nceiling_count=0;
  Real T_ceiling_ = ptrml->T_ceiling;
  Kokkos::parallel_reduce("Temp_ceiling", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, int &sum) {
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;
    Real dens=1.0, temp = 1.0, eint = 1.0;
    Real eint_ceiling = 1.0; Real deint = 0.0;
    Real dTemp = 0.0;
    dens = w0(m,IDN,k,j,i);
    if (use_e) {
      temp = w0(m,IEN,k,j,i)/dens*gm1;
      eint = w0(m,IEN,k,j,i);
      if(temp>T_ceiling_) {
        eint_ceiling = T_ceiling_ * dens/gm1;
        deint = eint_ceiling - eint;
        sum ++;
        // Now update internal energy
        u0(m,IEN,k,j,i) = u0(m,IEN,k,j,i) + bta_time * deint;
      }
    } else {
      temp = w0(m,ITM,k,j,i);
      eint = w0(m,ITM,k,j,i)*dens/gm1;
      if(temp>T_ceiling_) {
        dTemp = T_ceiling_ - temp;
        deint = dTemp*dens/gm1;
        sum ++;
        // Now update internal energy
        u0(m,IEN,k,j,i) = u0(m,IEN,k,j,i) + bta_time * deint;
      }
    }
  }, Kokkos::Sum<int>(nceiling_count));
  // Can print a statement here if temperature ceiling is applied
  #if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &nceiling_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (global_variable::my_rank == 0) {
    if (ptrml->ndiag>0 && pm->ncycle % ptrml->ndiag == 0 && bta_time > 0.99) {
      if (nceiling_count >0) {
        std::cout << " Temperature ceiling implementation is active." << std::endl;
        std::cout << " nceiling_count=" << nceiling_count << std::endl;
      }
    }
  }
  return;
}

// ----------------------------------------------------------------------------------------
// ! \fn void SourceTerms::AdjustTempTFloor()
// ! \brief Set temperature to the average of temperature in neighbouring cells
// ! whenever temperature drops below the floor value.
// NOTE source terms must all be computed using primitive (w0) and NOT conserved (u0) vars
void AdjustTempTFloor(Mesh *pm, const Real bdt, DvceArray5D<Real> &u0,
                   const DvceArray5D<Real> &w0, const EOS_Data &eos_data) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  int is = indcs.is, nx1 = indcs.nx1;
  int js = indcs.js, nx2 = indcs.nx2;
  int ks = indcs.ks, nx3 = indcs.nx3;
  const int nmkji = (pmbp->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  Real bta_time = bdt/pm->dt;
  Real use_e = eos_data.use_e;
  Real gamma = eos_data.gamma;
  Real gm1 = gamma - 1.0;
  Real temp_unit = pmbp->punit->temperature_cgs();
  int nfloor_count=0;
  Real T_floor_ = ptrml->T_floor;
  Kokkos::parallel_reduce("Temp_floor", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, int &sum) {
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;
    Real dens=1.0, temp = 1.0, eint = 1.0;
    Real eint_new = 1.0; Real deint = 0.0;
    dens = w0(m,IDN,k,j,i);
    if (use_e) {
      eint = w0(m,IEN,k,j,i);
      Real eint_km1 = w0(m,IEN,k-1,j,i);
      Real eint_kp1 = w0(m,IEN,k+1,j,i);
      Real eint_jm1 = w0(m,IEN,k,j-1,i);
      Real eint_jp1 = w0(m,IEN,k,j+1,i);
      Real eint_im1 = w0(m,IEN,k,j,i-1);
      Real eint_ip1 = w0(m,IEN,k,j,i+1);
      temp = eint/dens*gm1;
      if(temp<T_floor_) {
        eint_new = (1./6.) * (eint_km1 + eint_kp1 + eint_jm1 +
                              eint_jp1 + eint_im1 + eint_ip1);
        deint = eint_new - eint;
        // Now update internal energy
        u0(m,IEN,k,j,i) += deint;
        sum ++;
      }
    } else {
      temp = w0(m,ITM,k,j,i);
      eint = temp*dens/gm1;
      Real eint_km1 = w0(m,ITM,k-1,j,i)*w0(m,IDN,k-1,j,i)/gm1;
      Real eint_kp1 = w0(m,ITM,k+1,j,i)*w0(m,IDN,k+1,j,i)/gm1;
      Real eint_jm1 = w0(m,ITM,k,j-1,i)*w0(m,IDN,k,j-1,i)/gm1;
      Real eint_jp1 = w0(m,ITM,k,j+1,i)*w0(m,IDN,k,j+1,i)/gm1;
      Real eint_im1 = w0(m,ITM,k,j,i-1)*w0(m,IDN,k,j,i-1)/gm1;
      Real eint_ip1 = w0(m,ITM,k,j,i+1)*w0(m,IDN,k,j,i+1)/gm1;
      if(temp<T_floor_) {
        eint_new = (1./6.) * (eint_km1 + eint_kp1 + eint_jm1 +
                              eint_jp1 + eint_im1 + eint_ip1);
        deint = eint_new - eint;
        // Now update internal energy
        u0(m,IEN,k,j,i) += deint;
        sum ++;
      }
    }
  }, Kokkos::Sum<int>(nfloor_count));
  // Can print a statement here if temperature floor is applied
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &nfloor_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (global_variable::my_rank == 0) {
    if (ptrml->ndiag>0 && pm->ncycle % ptrml->ndiag == 0 && bta_time > 0.99) {
      if (nfloor_count >0) {
        std::cout << " Temperature floor implementation is active." << std::endl;
        std::cout << " nfloor_count=" << nfloor_count << std::endl;
      }
    }
  }
  return;
}

// Refinement criteria

void UserRefine(MeshBlockPack* pmbp) {
  // capture variables for kernels
  Mesh *pm = pmbp->pmesh;
  auto &indcs = pm->mb_indcs;
  int &is = indcs.is, nx1 = indcs.nx1;
  int &js = indcs.js, nx2 = indcs.nx2;
  int &ks = indcs.ks, nx3 = indcs.nx3;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;

  // check (on device) Hydro/MHD refinement conditions over all MeshBlocks
  auto refine_flag_ = pm->pmr->refine_flag;
  int nmb = pmbp->nmb_thispack;
  int mbs = pm->gids_eachrank[global_variable::my_rank];
  bool is_hydro = (pmbp->phydro != nullptr) ? true : false;
  bool is_mhd = (pmbp->pmhd != nullptr) ? true : false;
  EOS_Data &eos = (is_mhd) ?
                  pmbp->pmhd->peos->eos_data : pmbp->phydro->peos->eos_data;
  Real use_e = eos.use_e;
  Real gm1 = eos.gamma - 1.0;

  // maximum intercell density ratio
  Real density_ratio_threshold = ptrml->density_ratio_threshold;
  Real vel2_rms_threshold = ptrml->vel2_rms_threshold;
  Real T_max_threshold = ptrml->T_max_threshold;
  Real T_min_threshold = ptrml->T_min_threshold;
  if ((is_hydro) || (is_mhd)) {
    // if (global_variable::my_rank == 0) {printf("UserRefine\n");}
    auto &w0 = (is_hydro)? pmbp->phydro->w0 : pmbp->pmhd->w0;
    par_for_outer("RefineCond",DevExeSpace(), 0, 0, 0, (nmb-1),
    KOKKOS_LAMBDA(TeamMember_t tmember, const int m) {
      // density ratio threshold between neighbouring cells
      Real density_ratio_max = 1.0;
      Real vel2_rms = 0.0;
      Real counter = 0.0;
      int trigger = 0;
      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tmember, nkji),
      [=](const int idx, Real& dratio_max, Real &vel2_sq_sum, int &trigger_, Real &counter_) {
        int k = (idx)/nji;
        int j = (idx - k*nji)/nx1;
        int i = (idx - k*nji - j*nx1) + is;
        j += js;
        k += ks;
        vel2_sq_sum += SQR(w0(m,IVY,k,j,i));
        counter_ += 1.0;

        dratio_max = fmax(w0(m,IDN,k,j,i-1)/w0(m,IDN,k,j,i),1.0);
        dratio_max = fmax(w0(m,IDN,k,j,i+1)/w0(m,IDN,k,j,i),dratio_max);
        dratio_max = fmax(w0(m,IDN,k,j,i)/w0(m,IDN,k,j,i-1),dratio_max);
        dratio_max = fmax(w0(m,IDN,k,j,i)/w0(m,IDN,k,j,i+1),dratio_max);

        if (nx2>1) {
          dratio_max = fmax(w0(m,IDN,k,j-1,i)/w0(m,IDN,k,j,i),dratio_max);
          dratio_max = fmax(w0(m,IDN,k,j+1,i)/w0(m,IDN,k,j,i),dratio_max);
          dratio_max = fmax(w0(m,IDN,k,j,i)/w0(m,IDN,k,j-1,i),dratio_max);
          dratio_max = fmax(w0(m,IDN,k,j,i)/w0(m,IDN,k,j+1,i),dratio_max);
        }
        if (nx3>1) {
          dratio_max = fmax(w0(m,IDN,k-1,j,i)/w0(m,IDN,k,j,i),dratio_max);
          dratio_max = fmax(w0(m,IDN,k+1,j,i)/w0(m,IDN,k,j,i),dratio_max);
          dratio_max = fmax(w0(m,IDN,k,j,i)/w0(m,IDN,k-1,j,i),dratio_max);
          dratio_max = fmax(w0(m,IDN,k,j,i)/w0(m,IDN,k+1,j,i),dratio_max);
        }
        Real temp=1.0;
        if (use_e) {
          temp = w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1;
        } else {
          temp = w0(m,ITM,k,j,i);
        }
        if ((temp < T_max_threshold) and (temp > T_min_threshold)) {
          trigger_ += 1;
        }
        else trigger_ += 0;

      },Kokkos::Max<Real>(density_ratio_max), Kokkos::Sum<Real>(vel2_rms)
       , Kokkos::Sum<int>(trigger), Kokkos::Sum<Real>(counter));

      vel2_rms = sqrt(vel2_rms/counter);
      if ((density_ratio_max > density_ratio_threshold) || (vel2_rms>vel2_rms_threshold) || (trigger>0)) {
        refine_flag_.d_view(m+mbs) = 1;
      }
      else if ((density_ratio_max < 0.25 * density_ratio_threshold) || ((vel2_rms < 0.25*vel2_rms_threshold)) || (trigger==0)) {
        refine_flag_.d_view(m+mbs) = -1;
      }
    });
  }
}


// ----------------------------------------------------------------------------------------
// ! \fn void UserWorkInLoop()
// ! \brief Function called in hydro or mhd tasks in "after_timeintegrator" stage
void UserWorkInLoop(Mesh *pm) {
  // Frame tracking

  if (ptrml->ndiag>0 && pm->ncycle % ptrml->ndiag == 0) {
    Diagnostic(pm);
  }
  if (ptrml->use_frame_tracking && pm->time > ptrml->t_frame_tracking_start && pm->ncycle % ptrml->n_frame_track == 0) FrameTracking(pm);
}

// ----------------------------------------------------------------------------------------
// ! \fn void UserWorkInLoop::FrameTracking()
// ! \brief Add z-momentum such that the interface stays close to the center and does not interact with the boundaries.

void FrameTracking(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  int &ng = indcs.ng;
  int is = indcs.is, nx1 = indcs.nx1;
  int js = indcs.js, nx2 = indcs.nx2;
  int ks = indcs.ks, nx3 = indcs.nx3;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng) : 1;
  int nmb1 = pmbp->nmb_thispack - 1;
  const int nmkji = (pmbp->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  auto &size = pmbp->pmb->mb_size;
  bool is_hydro = (pmbp->phydro != nullptr) ? true : false;
  bool is_mhd = (pmbp->pmhd != nullptr) ? true : false;
  EOS_Data &eos = (is_mhd) ?
                  pmbp->pmhd->peos->eos_data : pmbp->phydro->peos->eos_data;
  Real use_e = eos.use_e;
  Real gamma = eos.gamma;
  Real gm1 = gamma - 1.0;

  Real T_ft_lo = ptrml->T_ft_lo;
  Real T_ft_hi = ptrml->T_ft_hi;

  Real mass_ft = 0.0;
  Real mean_z_ft = 0.0;
  Real mean_velz_ft = 0.0;
  Real min_z_ft = pm->mesh_size.x3max;
  Real max_z_ft = pm->mesh_size.x3min;

  Real ft_uppzlim = ptrml->ft_uppzlim; // Maximum allowed z-value of frame tracked gas (used when ft_use_uppzlim is true)
  Real ft_lowzlim = ptrml->ft_lowzlim; // Minimum allowed z-value of frame tracked gas (used when ft_use_lowzlim is true)
  Real mean_z_frame_tracking = ptrml->z_interface; // We'll try to keep the interface near the initial value
  Real max_vz_frame_tracking = ptrml->max_vz_frame_tracking; // Maximum velocity that we impart to the domain during frame-tracking

  auto &u0 = (is_mhd) ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  auto &w0 = (is_mhd) ? pmbp->pmhd->w0 : pmbp->phydro->w0;

  Kokkos::parallel_reduce("frametracking", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, Real &massft_, Real &mean_z_ft_, Real &mean_velz_ft_, Real &min_z_ft_, Real &max_z_ft_) {
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;
    Real dens=1.0, temp = 1.0, eint = 1.0;
    dens = w0(m,IDN,k,j,i);
    Real dvol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
    if (use_e) {
      temp = w0(m,IEN,k,j,i)/dens*gm1;
      eint = w0(m,IEN,k,j,i);
    } else {
      temp = w0(m,ITM,k,j,i);
      eint = w0(m,ITM,k,j,i)*dens/gm1;
    }

    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    int nx3 = indcs.nx3;
    Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);

    if ((temp<T_ft_hi) && (temp>T_ft_lo)){
      massft_ += dens*dvol;
      mean_z_ft_ += x3v*dens*dvol;
      mean_velz_ft_ += w0(m,IVZ,k,j,i)*dens*dvol;
      min_z_ft_ = fmin(x3v,min_z_ft_);
      max_z_ft_ = fmax(x3v,max_z_ft_);
    }
    
  }, Kokkos::Sum<Real>(mass_ft), 
    Kokkos::Sum<Real>(mean_z_ft),
    Kokkos::Sum<Real>(mean_velz_ft),
    Kokkos::Min<Real>(min_z_ft),
    Kokkos::Max<Real>(max_z_ft));
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &mass_ft, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &mean_z_ft, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &mean_velz_ft, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &min_z_ft, 1, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &max_z_ft, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
#endif

  if (global_variable::my_rank == 0 && ptrml->ndiag>0 && pm->ncycle % ptrml->ndiag == 0) {
    std::cout << " mass_ft=" << mass_ft << std::endl;
  }
  if(mass_ft > 0.){
    mean_z_ft = mean_z_ft/mass_ft;
    mean_velz_ft = mean_velz_ft/mass_ft;

    // First we try to maintain the mean velz of gas at T_ft close to zero
    Real vel_z_shift = -mean_velz_ft;

    if(ptrml->ft_use_uppzlim || ptrml->ft_use_lowzlim) {
      Real z_shift = vel_z_shift*pm->dt;
      // If some gas goes too close to the bottom boundary we add a positive shift
      if(ptrml->ft_use_lowzlim && (min_z_ft+z_shift) < ft_lowzlim) z_shift = ft_lowzlim - min_z_ft; 
      // If some gas goes too close to the upper boundary we apply a negative shift
      if(ptrml->ft_use_uppzlim && (max_z_ft+z_shift) > ft_uppzlim) z_shift = ft_uppzlim - max_z_ft;
      // Now update vel_z_shift after applying these conditions 
      vel_z_shift = z_shift/pm->dt;
    }
    // We set a maximum velocity that the frame-tracking can impart
    if(fabs(vel_z_shift) > fabs(max_vz_frame_tracking)) vel_z_shift = vel_z_shift * fabs(max_vz_frame_tracking/vel_z_shift);

    if (global_variable::my_rank == 0) {
      std::cout << " mean_z_ft=" << mean_z_ft << std::endl;
      std::cout << " mean_velz_ft=" << mean_velz_ft << std::endl;
      std::cout << " min_z_ft=" << min_z_ft << std::endl;
      std::cout << " max_z_ft=" << max_z_ft << std::endl;
      std::cout << " vel_z_shift=" << vel_z_shift << std::endl;
    }
    // Now perform the frame shift ()
    par_for("velzframetracking", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      w0(m,IVZ,k,j,i) += vel_z_shift;
    });
    // Now PrimToCons 
    if (is_hydro) {
      pmbp->phydro->peos->PrimToCons(w0,u0,0,(n1-1),0,(n2-1),0,(n3-1));
    } else if (is_mhd) {
      auto &bcc0 = pmbp->pmhd->bcc0;
      pmbp->pmhd->peos->PrimToCons(w0,bcc0,u0,0,(n1-1),0,(n2-1),0,(n3-1));
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
// ! \fn void UserWorkInLoop::Diagnostic()
// ! \brief Compute volume and mass averages and print them to the screen
void Diagnostic(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  int is = indcs.is, nx1 = indcs.nx1;
  int js = indcs.js, nx2 = indcs.nx2;
  int ks = indcs.ks, nx3 = indcs.nx3;
  const int nmkji = (pmbp->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  auto &size = pmbp->pmb->mb_size;
  bool is_mhd = (pmbp->pmhd != nullptr) ? true : false;
  EOS_Data &eos = (is_mhd) ?
                  pmbp->pmhd->peos->eos_data : pmbp->phydro->peos->eos_data;

  auto &w0 = (is_mhd) ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  DvceArray5D<Real> bcc0; // Initialize as an empty view
  if(is_mhd) bcc0 = pmbp->pmhd->bcc0;

  Real T_peak_lo = ptrml->T_peak_lo;
  Real T_peak_hi = ptrml->T_peak_hi;
  Real velocity  = ptrml->velocity;

  // For calculating cooling rate

  Real use_e = eos.use_e;
  Real gamma = eos.gamma;
  Real gm1 = gamma - 1.0;

  Real dtnew = std::numeric_limits<Real>::max();

  Real min_dens = std::numeric_limits<Real>::max();
  Real min_vtot = std::numeric_limits<Real>::max();
  Real min_temp = std::numeric_limits<Real>::max();
  Real min_eint = std::numeric_limits<Real>::max();
  Real max_dens = std::numeric_limits<Real>::min();
  Real max_vtot = std::numeric_limits<Real>::min();
  Real max_temp = std::numeric_limits<Real>::min();
  Real max_eint = std::numeric_limits<Real>::min();
  Real tot_mass = 0.;
  Real tot_eint = 0.;
  Real tot_vsq = 0.0;
  Real tot_bsq = 0.0;
  Real tot_vol = 0.0;

  Real min_z_peak = 10.0;
  Real max_z_peak = -10.0;

  Real min_z_shear = 10.0;
  Real max_z_shear = -10.0;

  Real min_z_vely = 10.0;
  Real max_z_vely = -10.0;

  auto force_ = pmbp->pturb->force;
  // find smallest (e/cooling_rate) in each cell
  Kokkos::parallel_reduce("diagnostic", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, Real &min_dt, Real &min_d, Real &min_v, Real &min_t,
  Real &min_e, Real &max_d, Real &max_v, Real &max_t, Real &max_e, Real &tot_m,
  Real &tot_e, Real &totvol_, Real &totvsq_, Real &totbsq_,
  Real &min_z_peak_, Real &max_z_peak_, Real &min_z_shear_, Real &max_z_shear_,
  Real &min_z_vely_, Real &max_z_vely_) {
    // compute m,k,j,i indices of thread and call function
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;

    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    int nx3 = indcs.nx3;
    Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);

    Real dx = fmin(fmin(size.d_view(m).dx1,size.d_view(m).dx2),size.d_view(m).dx3);
    Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;

    // temperature in cgs unit
    Real temp = 1.0, eint = 1.0;
    Real dens = w0(m,IDN,k,j,i);
    Real velx = w0(m,IVX,k,j,i);
    Real vely = w0(m,IVY,k,j,i);
    Real velz = w0(m,IVZ,k,j,i);
    if (use_e) {
      temp = w0(m,IEN,k,j,i)/dens*gm1;
      eint = w0(m,IEN,k,j,i);
    } else {
      temp = w0(m,ITM,k,j,i);
      eint = w0(m,ITM,k,j,i)*dens/gm1;
    }

    Real vtot = sqrt(SQR(w0(m,IVX,k,j,i))+SQR(w0(m,IVY,k,j,i))+SQR(w0(m,IVZ,k,j,i)));
    min_dt = fmin(dx/sqrt(gamma*temp), min_dt);
    min_d = fmin(dens, min_d);
    min_v = fmin(vtot,min_v);
    min_t = fmin(temp, min_t);
    min_e = fmin(eint, min_e);
    max_d = fmax(dens, max_d);
    max_v = fmax(vtot,max_v);
    max_t = fmax(temp, max_t);
    max_e = fmax(eint, max_e);
    tot_m += dens*vol;
    tot_e += w0(m,IEN,k,j,i)*vol;


    totvol_ += vol;
    totvsq_ += vol*SQR(vtot);
    if(is_mhd) {
      totbsq_ += vol*(SQR(bcc0(m,IBX,k,j,i))+SQR(bcc0(m,IBY,k,j,i))+SQR(bcc0(m,IBZ,k,j,i)));
    }

    if ((temp<T_peak_hi) && (temp>T_peak_lo)){
      min_z_peak_ = fmin(x3v,min_z_peak_);
      max_z_peak_ = fmax(x3v,max_z_peak_);
    }

    if (velx > -0.45 * velocity) {
      min_z_shear_ = fmin(x3v,min_z_shear_); // minimum height of gas with v_x > -0.45 velocity
    }

    if (velx < 0.45 * velocity) {
      max_z_shear_ = fmax(x3v,max_z_shear_); // maximum height of gas with v_x < 0.45 velocity
    }

    if (vely > 5e-2 * velocity){
      min_z_vely_ = fmin(x3v, min_z_vely_);
      max_z_vely_ = fmax(x3v, max_z_vely_);
    }
  }, Kokkos::Min<Real>(dtnew),
     Kokkos::Min<Real>(min_dens),
     Kokkos::Min<Real>(min_vtot),
     Kokkos::Min<Real>(min_temp),
     Kokkos::Min<Real>(min_eint),
     Kokkos::Max<Real>(max_dens),
     Kokkos::Max<Real>(max_vtot),
     Kokkos::Max<Real>(max_temp),
     Kokkos::Max<Real>(max_eint),
     Kokkos::Sum<Real>(tot_mass),
     Kokkos::Sum<Real>(tot_eint),
     Kokkos::Sum<Real>(tot_vol),
     Kokkos::Sum<Real>(tot_vsq),
     Kokkos::Sum<Real>(tot_bsq),
     Kokkos::Min<Real>(min_z_peak),
     Kokkos::Max<Real>(max_z_peak),
     Kokkos::Min<Real>(min_z_shear),
     Kokkos::Max<Real>(max_z_shear),
     Kokkos::Min<Real>(min_z_vely),
     Kokkos::Max<Real>(max_z_vely));
  Real dt_hyd  = (is_mhd) ? pmbp->pmhd->dtnew        : pmbp->phydro->dtnew;
  auto *pcond = (is_mhd) ? pmbp->pmhd->pcond : pmbp->phydro->pcond;
  Real dt_cond = FLT_MAX;
  if (pcond != nullptr) {
    Real dt_cond = (is_mhd) ? pmbp->pmhd->pcond->dtnew : pmbp->phydro->pcond->dtnew;
  }
  Real dt_src  = FLT_MAX;
  auto *psrc = (is_mhd) ? pmbp->pmhd->psrc : pmbp->phydro->psrc;
  if (psrc != nullptr) {
    dt_src  = psrc->dtnew;
  }
  Real dt_user  = FLT_MAX;
  if (pmbp->pmesh->pgen->user_dt) {
    dt_user = pmbp->pmesh->pgen->dtnew;
  }
#if MPI_PARALLEL_ENABLED
  Real m_min[12] = {dtnew,min_dens,min_vtot,min_temp,min_eint,dt_hyd,dt_cond,dt_src,dt_user,min_z_peak,min_z_shear,min_z_vely};
  Real m_max[7] = {max_dens,max_vtot,max_temp,max_eint,max_z_peak,max_z_shear,max_z_vely};
  Real gm_min[12];
  Real gm_max[7];
  Real loc_sum[5] = {tot_mass,tot_eint,tot_vol,tot_vsq,tot_bsq};
  Real glob_sum[5];
  for(int index=0; index<5; index++){
    glob_sum[index]=0.0;
  }
  //MPI_Allreduce(MPI_IN_PLACE, &dtnew, 1, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(m_min, gm_min, 12, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(m_max, gm_max, 7, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(loc_sum, glob_sum, 5, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  dtnew = gm_min[0];
  min_dens = gm_min[1];
  min_vtot = gm_min[2];
  min_temp = gm_min[3];
  min_eint = gm_min[4];
  dt_hyd   = gm_min[5];
  dt_cond  = gm_min[6];
  dt_src   = gm_min[7];
  dt_user  = gm_min[8];
  min_z_peak=gm_min[9];
  min_z_shear=gm_min[10];
  min_z_vely=gm_min[11];
  max_dens = gm_max[0];
  max_vtot = gm_max[1];
  max_temp = gm_max[2];
  max_eint = gm_max[3];
  max_z_peak=gm_max[4];
  max_z_shear=gm_max[5];
  max_z_vely=gm_max[6];

  tot_mass = glob_sum[0];
  tot_eint = glob_sum[1];
  tot_vol = glob_sum[2];
  tot_vsq = glob_sum[3];
  tot_bsq = glob_sum[4];
#endif
  Real totcoolrate = ptrml->tot_coolrate;
  ptrml->tot_mass = tot_mass;
  // ptrml->tot_coolrate = totcoolrate;
  Real vrms = sqrt(tot_vsq/tot_vol);
  Real brms = sqrt(tot_bsq/tot_vol);
  if (global_variable::my_rank == 0) {
    std::cout << std::setprecision(16)
              << " min_d=" << min_dens << " max_d=" << max_dens << std::endl
              << " min_T=" << min_temp << " max_T=" << max_temp << std::endl
              << " min_v=" << min_vtot    << " max_v=" << max_vtot    << std::endl
              << " min_e=" << min_eint << " max_e=" << max_eint << std::endl
              << " tot_m=" << tot_mass << " tot_e=" << tot_eint << std::endl
              << " dt_temp=" << dtnew   << " dt_hyd=" << dt_hyd << std::endl
              << " dt_cond=" << dt_cond << " dt_src=" << dt_src << std::endl
              << " dt_user=" << dt_user << " totcoolrate=" << totcoolrate << std::endl
              << " v_rms=" << vrms << " b_rms=" << brms << std::endl
              << " min_z_peak=" << min_z_peak << std::endl << " max_z_peak=" << max_z_peak << std::endl
              << " min_z_shear=" << min_z_shear << std::endl << " max_z_shear=" << max_z_shear << std::endl
              << " min_z_vely=" << min_z_vely << std::endl << " max_z_vely=" << max_z_vely << std::endl;
  }
  return;
}



//----------------------------------------------------------------------------------------
//! \fn UserHistOutput
//! \brief Sets user-defined history output

void UserHistOutput(HistoryData *pdata, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  bool is_mhd = (pmbp->pmhd != nullptr) ? true : false;
  int &nscalars = (is_mhd) ?
                  pmbp->pmhd->nscalars : pmbp->phydro->nscalars;
  int &nfluid = (is_mhd) ? pmbp->pmhd->nmhd : pmbp->phydro->nhydro;
  int nsum = 30;
  // Store nsum0 in the first array and rest in the second array
  int nsum0 = 18; 
  if (nscalars>0) nsum = nsum+2;
  int nsum1 = nsum-nsum0;
  int nother_hist = 7; // history variables that are not global sums
  int nuser = nsum+nother_hist;
  int nuser1 = nsum1+nother_hist; // nother_hist are stored in the second array
  pdata->nhist = nuser;
    // std::cout << "nhist = " << pdata->nhist << std::endl;
  if (pdata->nhist > NHISTORY_VARIABLES) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "User history function specified pdata->nhist larger than"
              << " NHISTORY_VARIABLES" << std::endl;
    exit(EXIT_FAILURE);
  }
  int nhi=0;
  pdata->label[nhi++] = "Mflx_top "; // 0
  pdata->label[nhi++] = "Eflx_top ";
  pdata->label[nhi++] = "Mflx_bot ";
  pdata->label[nhi++] = "vysq_pk ";
  pdata->label[nhi++] = "V_pk ";
  pdata->label[nhi++] = "Zavg_pk "; // 5
  pdata->label[nhi++] = "vysq_c ";
  pdata->label[nhi++] = "momx_c ";
  pdata->label[nhi++] = "kex_c ";
  pdata->label[nhi++] = "V_c ";
  pdata->label[nhi++] = "vysq_ci "; // 10
  pdata->label[nhi++] = "V_ci ";
  pdata->label[nhi++] = "vysq_h ";
  pdata->label[nhi++] = "momx_h ";
  pdata->label[nhi++] = "kex_h ";
  pdata->label[nhi++] = "V_h "; // 15
  pdata->label[nhi++] = "vysq_hi ";
  pdata->label[nhi++] = "V_hi ";
  pdata->label[nhi++] = "area1 ";
  pdata->label[nhi++] = "area2 ";
  pdata->label[nhi++] = "area4 "; // 20
  pdata->label[nhi++] = "vd_pk ";
  pdata->label[nhi++] = "vd_pk_dx ";
  pdata->label[nhi++] = "ell_T_pk ";
  pdata->label[nhi++] = "vd_ci ";
  pdata->label[nhi++] = "vd_ci_dx "; // 25
  pdata->label[nhi++] = "ell_T_ci ";
  pdata->label[nhi++] = "vd_hi ";
  pdata->label[nhi++] = "vd_hi_dx ";
  pdata->label[nhi++] = "ell_T_hi "; // 29
  if (nscalars>0){
    pdata->label[nhi++] = "scalar_xmom ";
    pdata->label[nhi++] = "scalar_xener ";
  }
  pdata->label[nhi++] = "tot_cool_rate ";
  pdata->label[nhi++] = "zmin_shear ";
  pdata->label[nhi++] = "zmax_shear ";
  pdata->label[nhi++] = "zmin_peak ";
  pdata->label[nhi++] = "zmax_peak ";
  pdata->label[nhi++] = "zmin_vy ";
  pdata->label[nhi++] = "zmax_vy ";

  EOS_Data &eos = (is_mhd) ?
                  pmbp->pmhd->peos->eos_data : pmbp->phydro->peos->eos_data;

  Real use_e = eos.use_e;
  Real gm1 = eos.gamma - 1.0;
  Real gamma = eos.gamma;
  // capture class variables for kernel
  //auto &u0_ = pmbp->phydro->u0;
  auto &w0 = (is_mhd) ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  auto &size = pmbp->pmb->mb_size;

  // loop over all MeshBlocks in this pack
  auto &indcs = pm->mb_indcs;
  int is = indcs.is; int nx1 = indcs.nx1;
  int js = indcs.js; int nx2 = indcs.nx2;
  int ks = indcs.ks; int nx3 = indcs.nx3;

  // Get the extent of the boxes along z
  Real lxmax3 = pm->mesh_size.x3max;
  Real lxmin3 = pm->mesh_size.x3min;

  Real T_peak = ptrml->T_peak;
  Real T_peak_lo = ptrml->T_peak_lo;
  Real T_peak_hi = ptrml->T_peak_hi;
  Real T_cold_hi = ptrml->T_cold_hi;
  Real T_hot_lo  = ptrml->T_hot_lo;
  Real velocity  = ptrml->velocity;

  Real min_z_peak = 10.0;
  Real max_z_peak = -10.0;

  Real min_z_shear = 10.0;
  Real max_z_shear = -10.0;

  Real min_z_vely = 10.0;
  Real max_z_vely = -10.0;

  auto force_ = pmbp->pturb->force;

  const int nmkji = (pmbp->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  array_sum::GlobalSum sum_this_mb0, sum_this_mb1;
  // store data into hdata array
  for (int n=0; n<NREDUCTION_VARIABLES; ++n) {
    sum_this_mb0.the_array[n] = 0.0;
    sum_this_mb1.the_array[n] = 0.0;
  }


  Kokkos::parallel_reduce("UserHistSums",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &mb_sum0, array_sum::GlobalSum &mb_sum1,
  Real &min_z_peak_, Real &max_z_peak_, Real &min_z_shear_, Real &max_z_shear_,
  Real &min_z_vely_, Real &max_z_vely_) {
    // compute n,k,j,i indices of thread
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;

    Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;

    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    int nx3 = indcs.nx3;
    Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
    Real dz = size.d_view(m).dx3;
    Real dA = size.d_view(m).dx1*size.d_view(m).dx2;

    Real dens = w0(m,IDN,k,j,i);
    Real velx = w0(m,IVX,k,j,i);
    Real vely = w0(m,IVY,k,j,i);
    Real velz = w0(m,IVZ,k,j,i);

    Real temp = 0.0;
    Real eint = 0.0;
    if (use_e) {
      temp = w0(m,IEN,k,j,i)/dens*gm1;
      eint = w0(m,IEN,k,j,i);
    } else {
      temp = w0(m,ITM,k,j,i);
      eint = w0(m,ITM,k,j,i)*dens/gm1;
    }
    Real ekin = 0.5*w0(m,IDN,k,j,i) * (SQR(w0(m,IVX,k,j,i))+SQR(w0(m,IVY,k,j,i))+SQR(w0(m,IVZ,k,j,i)));
    Real etot = eint+ekin;
    Real work = (gamma-1.0)*eint;

    Real scalar_mass_fraction = 0.0;
    if (nscalars > 0){
      scalar_mass_fraction = w0(m,nfluid,k,j,i);
    }


    array_sum::GlobalSum hvars0, hvars1;
    // Let's put first 18 variables in hvars0 and the rest in hvars1
    if (fabs(x3v-lxmax3) < dz) {
      hvars0.the_array[0] = dens*dA*velz; // top mass flux
      hvars0.the_array[1] = (etot+work)*dA*velz; // top energy flux
    }
    else if (fabs(x3v-lxmin3) < dz) {
      hvars0.the_array[2] = dens*dA*velz; // mass flux bottom
    }
    if ((temp<T_peak_hi) && (temp>T_peak_lo)){
      hvars0.the_array[3] = SQR(vely)*vol; // rms y velocity near T_peak
      hvars0.the_array[4] = vol; // volume near T_peak
      hvars0.the_array[5] = x3v*vol;// average height of gas near T_peak
      min_z_peak_ = fmin(x3v,min_z_peak_);
      max_z_peak_ = fmax(x3v,max_z_peak_);
    }

    if ((temp<T_cold_hi)){
      hvars0.the_array[6] = SQR(vely)*vol; // rms y velocity in T_cold gas
      hvars0.the_array[7] = dens*velx*vol; // x-momentum in T_cold gas
      hvars0.the_array[8] = 0.5*dens*SQR(velx)*vol; // x-kinetic energy in T_cold gas
      hvars0.the_array[9] = vol; // volume of T_cold gas
    }

    if ((temp>T_cold_hi) && (temp<T_peak_lo)) {
      hvars0.the_array[10] = SQR(vely)*vol; // rms y velocity just above T_cold
      hvars0.the_array[11] = vol; // volume just above T_cold
    }

    if ((temp>T_hot_lo)){
      hvars0.the_array[12] = SQR(vely)*vol; // rms y velocity in T_hot gas
      hvars0.the_array[13] = dens*velx*vol; // x-momentum in T_hot gas
      hvars0.the_array[14] = 0.5*dens*SQR(velx)*vol; // x-kinetic energy in T_hot gas
      hvars0.the_array[15] = vol; // volume in T_hot gas
    }

    if ((temp<T_hot_lo) && (temp>T_peak_hi)){
      hvars0.the_array[16] = SQR(vely)*vol; // rms y velocity below T_hot
      hvars0.the_array[17] = vol; // volume below T_hot
    }
    
    // Get the area contributed by the surface above and to the right
    // of the current cell
    // Initialize Cube to store relative data
    Real v0,v1,v2,v3,v4,v5,v6,v7;
    Real iso = log10(T_peak);
    int ns = 0;
    for (int s = 1; s<8; s*=2) {
      if (use_e) {
        v0 = log10(gm1*w0(m,IEN,k  ,j  ,i  )/w0(m,IDN,k  ,j  ,i  )) - iso;
        v1 = log10(gm1*w0(m,IEN,k  ,j  ,i+s)/w0(m,IDN,k  ,j  ,i+s)) - iso;
        v2 = log10(gm1*w0(m,IEN,k  ,j+s,i+s)/w0(m,IDN,k  ,j+s,i+s)) - iso;
        v3 = log10(gm1*w0(m,IEN,k  ,j+s,i  )/w0(m,IDN,k  ,j+s,i  )) - iso;
        v4 = log10(gm1*w0(m,IEN,k+s,j  ,i  )/w0(m,IDN,k+s,j  ,i  )) - iso;
        v5 = log10(gm1*w0(m,IEN,k+s,j  ,i+s)/w0(m,IDN,k+s,j  ,i+s)) - iso;
        v6 = log10(gm1*w0(m,IEN,k+s,j+s,i+s)/w0(m,IDN,k+s,j+s,i+s)) - iso;
        v7 = log10(gm1*w0(m,IEN,k+s,j+s,i  )/w0(m,IDN,k+s,j+s,i  )) - iso;
      } else {
        v0 = log10(w0(m,ITM,k  ,j  ,i  )) - iso;
        v1 = log10(w0(m,ITM,k  ,j  ,i+s)) - iso;
        v2 = log10(w0(m,ITM,k  ,j+s,i+s)) - iso;
        v3 = log10(w0(m,ITM,k  ,j+s,i  )) - iso;
        v4 = log10(w0(m,ITM,k+s,j  ,i  )) - iso;
        v5 = log10(w0(m,ITM,k+s,j  ,i+s)) - iso;
        v6 = log10(w0(m,ITM,k+s,j+s,i+s)) - iso;
        v7 = log10(w0(m,ITM,k+s,j+s,i  )) - iso;
      }
      Cube c = Cube(v0, v1, v2, v3, v4, v5, v6, v7);
      hvars1.the_array[ns] = process_cube(c)*dA;
      ns++;
    }
    // get small scale diffusion velocity for gas selected to
    // be the peak cooling region
    if ((temp<T_hot_lo) && (temp>T_cold_hi)){
      // get local temperature gradient
      Real gradTx = 0.0;
      Real gradTy = 0.0;
      Real gradTz = 0.0;
      if (use_e) {
        gradTx = gm1*(w0(m,IEN,k,j,i+1)/w0(m,IDN,k,j,i+1)-w0(m,IEN,k,j,i-1)/w0(m,IDN,k,j,i-1))/(2.0*size.d_view(m).dx1);
        gradTy = gm1*(w0(m,IEN,k,j+1,i)/w0(m,IDN,k,j+1,i)-w0(m,IEN,k,j-1,i)/w0(m,IDN,k,j-1,i))/(2.0*size.d_view(m).dx2);
        gradTz = gm1*(w0(m,IEN,k+1,j,i)/w0(m,IDN,k+1,j,i)-w0(m,IEN,k-1,j,i)/w0(m,IDN,k-1,j,i))/(2.0*size.d_view(m).dx3);
      } else {
        gradTx = (w0(m,ITM,k,j,i+1)-w0(m,ITM,k,j,i-1))/(2.0*size.d_view(m).dx1);
        gradTy = (w0(m,ITM,k,j+1,i)-w0(m,ITM,k,j-1,i))/(2.0*size.d_view(m).dx2);
        gradTz = (w0(m,ITM,k+1,j,i)-w0(m,ITM,k-1,j,i))/(2.0*size.d_view(m).dx3);
      }
      Real gradT = sqrt(SQR(gradTx)+SQR(gradTy)+SQR(gradTz));
      // define temperature scale height
      Real ell_T = w0(m,ITM,k,j,i)/gradT;
      // get normal vector defined by temperature gradient
      Real nx = gradTx/gradT;
      Real ny = gradTy/gradT;
      Real nz = gradTz/gradT;

      // get velocity components in direction of temperature gradient
      // defined in the local region
      Real vin_km1 = w0(m,IVX,k-1,j,i)*nx + w0(m,IVY,k-1,j,i)*ny + w0(m,IVZ,k-1,j,i)*nz;
      Real vin_kp1 = w0(m,IVX,k+1,j,i)*nx + w0(m,IVY,k+1,j,i)*ny + w0(m,IVZ,k+1,j,i)*nz;
      Real vin_jm1 = w0(m,IVX,k,j-1,i)*nx + w0(m,IVY,k,j-1,i)*ny + w0(m,IVZ,k,j-1,i)*nz;
      Real vin_jp1 = w0(m,IVX,k,j+1,i)*nx + w0(m,IVY,k,j+1,i)*ny + w0(m,IVZ,k,j+1,i)*nz;
      Real vin_im1 = w0(m,IVX,k,j,i-1)*nx + w0(m,IVY,k,j,i-1)*ny + w0(m,IVZ,k,j,i-1)*nz;
      Real vin_ip1 = w0(m,IVX,k,j,i+1)*nx + w0(m,IVY,k,j,i+1)*ny + w0(m,IVZ,k,j,i+1)*nz;

      Real vin_gradx = (vin_ip1-vin_im1)/(2.0*size.d_view(m).dx1);
      Real vin_grady = (vin_jp1-vin_jm1)/(2.0*size.d_view(m).dx2);
      Real vin_gradz = (vin_kp1-vin_km1)/(2.0*size.d_view(m).dx3);
      Real vin = -1*ell_T*(vin_gradx*nx + vin_grady*ny + vin_gradz*nz);
      Real vin_dx = -2*size.d_view(m).dx1*(vin_gradx*nx + vin_grady*ny + vin_gradz*nz);
      if ((temp<T_peak_hi) && (temp>T_peak_lo)){
        hvars1.the_array[3] = vin*vol; // the small scale diffusive velocity near the cooling surface
        hvars1.the_array[4] = vin_dx*vol; // different measurement of vdiff
        hvars1.the_array[5] = ell_T*vol; // the temperature scale height near the cooling surface
      } else if (temp < T_peak_lo) {
        hvars1.the_array[6] = vin*vol; // vdiff in slightly colder gas
        hvars1.the_array[7] = vin_dx*vol; // vdiff in slightly colder gas
        hvars1.the_array[8] = ell_T*vol; // ell_T in slightly colder gas
      } else {
        hvars1.the_array[9] = vin*vol; // vdiff in slightly hotter gas
        hvars1.the_array[10] = vin_dx*vol; // vdiff in slightly hotter gas
        hvars1.the_array[11] = ell_T*vol; // ell_T in slightly hotter gas
      }
    }

    if (nscalars > 0){
      hvars1.the_array[12] = scalar_mass_fraction * dens * velx * vol; // scalar weighted x-momentum
      hvars1.the_array[13] = scalar_mass_fraction * 0.5 * dens * SQR(velx)*vol; // scalar weighted x-kin energy
    }
    // fill min/max and cooling parts of the_array with zeros
    for (int n=nsum1; n<nsum1+nother_hist; ++n) {
      hvars1.the_array[n] = 0.0;
    }
    // fill redundant variables with 0
    for (int n=nsum0; n<NREDUCTION_VARIABLES; ++n) {
      hvars0.the_array[n] = 0.0;
    }
    for (int n=nuser1; n<NREDUCTION_VARIABLES; ++n) {
      hvars1.the_array[n] = 0.0;
    }

    // sum into parallel reduce
    mb_sum0 += hvars0;
    mb_sum1 += hvars1;

    if (velx > -0.45 * velocity) {
      min_z_shear_ = fmin(x3v,min_z_shear_); // minimum height of gas with v_x > -0.45 velocity
    }

    if (velx < 0.45 * velocity) {
      max_z_shear_ = fmax(x3v,max_z_shear_); // maximum height of gas with v_x < 0.45 velocity
    }

    if (vely > 5e-2 * velocity){
      min_z_vely_ = fmin(x3v, min_z_vely_);
      max_z_vely_ = fmax(x3v, max_z_vely_);
    }
  }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb0),
     Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb1),
     Kokkos::Min<Real>(min_z_peak),
     Kokkos::Max<Real>(max_z_peak),
     Kokkos::Min<Real>(min_z_shear),
     Kokkos::Max<Real>(max_z_shear),
     Kokkos::Min<Real>(min_z_vely),
     Kokkos::Max<Real>(max_z_vely));

  #if MPI_PARALLEL_ENABLED
  Real m_min[3] = {min_z_peak,min_z_shear,min_z_vely};
  Real m_max[3] = {max_z_peak,max_z_shear,max_z_vely};
  Real gm_min[3];
  Real gm_max[3];
  //MPI_Allreduce(MPI_IN_PLACE, &dtnew, 1, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(m_min, gm_min, 3, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(m_max, gm_max, 3, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
  min_z_peak  = gm_min[0];
  min_z_shear = gm_min[1];
  min_z_vely  = gm_min[2];
  max_z_peak  = gm_max[0];
  max_z_shear = gm_max[1];
  max_z_vely  = gm_max[2];
#endif

  // store data into hdata array
  for (int n=0; n<nsum; ++n) {
    pdata->hdata[n] = (n < nsum0) ? sum_this_mb0.the_array[n] : sum_this_mb1.the_array[n-nsum0];
  }
  // The history variables added undergo an MPI_Reduce before writing the output.
  // Since we have already reduced the following variables, we shall only store them in the master processor.
  // For quantities that are already calculated in other parts of the code or which are maxima/minima,
  // we store their values in the master processer, since the current history function can only perform
  // a global sum
  if (global_variable::my_rank == 0){
    pdata->hdata[nsum  ] = ptrml->tot_coolrate;
    pdata->hdata[nsum+1] = min_z_shear;
    pdata->hdata[nsum+2] = max_z_shear;
    pdata->hdata[nsum+3] = min_z_peak;
    pdata->hdata[nsum+4] = max_z_peak;
    pdata->hdata[nsum+5] = min_z_vely;
    pdata->hdata[nsum+6] = max_z_vely;
  }
  else{
    for(int n=nsum; n<nsum+nother_hist; ++n){
      pdata->hdata[n] = 0.0;
    }
  }
  return;
}



void TRMLFinalWork(ParameterInput *pin, Mesh *pm) {
  delete ptrml;
}
// Some utility functions for calculating the initial density profile in MHD problems
// Calculate density equation value
KOKKOS_INLINE_FUNCTION
Real density_equation(Real rho, Real rho_hot, Real beta_hot, Real contrast, Real alpha_magdens) {
  Real left_side = rho + (pow(rho, 2.0*alpha_magdens) * contrast)/(beta_hot * pow(rho_hot, 2.0*alpha_magdens-1.0));
  Real right_side = (1.0 + 1.0/beta_hot) * rho_hot * contrast;
  return left_side - right_side;
}

// Calculate derivative of density equation
KOKKOS_INLINE_FUNCTION
Real density_equation_derivative(Real rho, Real rho_hot, Real beta_hot, Real contrast, Real alpha_magdens) {
  return 1.0 + 2.0 * alpha_magdens * pow(rho, 2.0*alpha_magdens-1.0) * contrast / (beta_hot * pow(rho_hot, 2.0*alpha_magdens-1.0));
}

KOKKOS_INLINE_FUNCTION
// Newton-Raphson solver for density equation
Real newton_density_solver(Real rho_guess, const Real rho_hot, const Real beta_hot, const Real contrast, 
  Real alpha_magdens, Real tol=1e-8, int max_iter=100) {
  Real rho = rho_guess;

  for (int iter = 0; iter < max_iter; ++iter) {
  // Function value
  Real f = density_equation(rho, rho_hot, beta_hot, contrast, alpha_magdens);

  // Derivative of function
  Real df = density_equation_derivative(rho, rho_hot, beta_hot, contrast, alpha_magdens);

  // Newton step
  Real rho_new = rho - f/df;

  // Check convergence
  if (abs(rho_new - rho) < tol) {
  return rho_new;
  }

  rho = rho_new;
}

return rho; // Return last value if max_iter reached
}

KOKKOS_INLINE_FUNCTION 
Real solve_density_profile_NR(const Real z, const Real rho_hot, const Real beta_hot,
  const Real contrast, Real alpha_magdens) {
  return newton_density_solver(rho_hot, rho_hot, beta_hot, contrast, alpha_magdens);
}

} //namespace

// Define a cooling function and a device cooling function