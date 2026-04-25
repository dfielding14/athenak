#include <iostream>
#include <cmath>
#include <string>
#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "pgen.hpp"
#include "units/units.hpp"
#include "utils/profile_reader.hpp"
#include "utils/sn_scheduler.hpp"
#include "utils/random.hpp"
#include "particles/particles.hpp"

//===========================================================================//
//                               Globals                                     //
//===========================================================================//

KOKKOS_INLINE_FUNCTION
void SetEquilibriumState(const DvceArray5D<Real> &u0,
  int m, int k, int j, int i,
  Real x1v, Real x2v, Real x3v,
  Real x1l, Real x1r, Real x2l, Real x2r,
  Real G, Real r_s, Real rho_s, Real m_g,
  Real a_g, Real z_g, Real r_m, Real rho_m,
  Real gm1, int IZS, int IDS, int IDL,
  Real dz, Real Z, Real Zsol,
  const ProfileReader &disk_profile);

KOKKOS_INLINE_FUNCTION
void SetCoolingFlowState(const DvceArray5D<Real> &u0,
                         int m, int k, int j, int i,
                         Real x1v, Real x2v, Real x3v,
                         Real gm1, const ProfileReader &profile);

KOKKOS_INLINE_FUNCTION
void SetRotation(const DvceArray5D<Real> &u0,
                 int m, int k, int j, int i,
                 Real x1v, Real x2v, Real x3v,
                 Real r_circ, Real v_circ);

KOKKOS_INLINE_FUNCTION
void SetDustScalars(const DvceArray5D<Real> &u0,
                    int m, int k, int j, int i,
                    int IZS, int IDS, int IDL,
                    Real dz, Real Z, Real Zsol);

KOKKOS_INLINE_FUNCTION
Real GravPot(Real x1, Real x2, Real x3,
             Real G, Real r_s, Real rho_s,
             Real M_gal, Real a_gal, Real z_gal,
             Real c_outer, Real rho_mean);

KOKKOS_INLINE_FUNCTION
Real LimitDustTransfer(const Real requested,
                       const Real donor_mass,
                       const Real dust_donor_cap_frac);

KOKKOS_INLINE_FUNCTION
Real ToroidalVectorPotential(const Real x1, const Real x2, const Real b0,
                             const Real r_soft);

KOKKOS_INLINE_FUNCTION
Real TangledVectorPotential(const int dir, const Real x1, const Real x2,
                            const Real x3, const int nmodes,
                            const DvceArray2D<Real> &k_modes,
                            const DvceArray2D<Real> &aka,
                            const DvceArray2D<Real> &akb);

void UserSource(Mesh* pm, const Real bdt);
void GravitySource(Mesh* pm, const Real bdt);
void SNSource(Mesh* pm, const Real bdt);
void DustSource(Mesh* pm, const Real bdt);
void UserBoundary(Mesh* pm);
void FreeProfile(ParameterInput *pin, Mesh *pm);
void RefinementCondition(MeshBlockPack* pmbp);
Real ReferenceMagneticField(const ProfileReader &profile, const Real r_ref,
                            const Real beta);
void InitializeMagneticField(ParameterInput *pin, MeshBlockPack *pmbp,
                             const std::string &bfield_type,
                             const Real b_ref);
void ApplyUserMagneticBoundary(Mesh *pm);
void AddUserBoundaryMagneticEnergy(Mesh *pm);

namespace {
  enum class MagneticFieldType { none, vertical, toroidal, tangled };

  // Constants for gravitational potential
  Real r_scale;
  Real rho_scale;
  Real m_gal;
  Real a_gal;
  Real z_gal;
  Real r_200;
  Real rho_mean;

  // Constants for rotation
  Real r_circ;
  Real v_circ;

  // Constant for metallicity
  Real Z;
  Real Z_sol = 0.02;

  // Constants for SN
  Real r_inj;
  Real e_sn;
  Real m_ej;
  Real Z_ej;
  Real sn_delay;

  // Constants for dust
  Real d_z_init;
  Real d_z_sn;
  Real rho_gr;
  Real a_gr_s;
  Real a_gr_l;
  Real min_dust_frac;
  Real dust_donor_cap_frac;
  Real accretion_nH_min;
  Real accretion_nH_max;
  Real accretion_T_max;
  Real accretion_refractory_frac;

  // Profiles
  ProfileReader profile_reader;
  ProfileReader disk_profile_reader;

  // Refinment condition threshold
  Real ddens_threshold;

  // SN injection persistent buffer
  DvceArray2D<Real> sn_centers_buffer;
  Kokkos::View<int> sn_counter;

  // SN injection flags
  int last_sn_detect_cycle = -1;
  int num_sn_this_cycle;

  // Passive scalar index offsets
  int scalar_IZS;  // metallicity
  int scalar_IDS;  // small dust
  int scalar_IDL;  // large dust

  MagneticFieldType ParseMagneticFieldType(const std::string &name) {
    if (name == "none") return MagneticFieldType::none;
    if (name == "vertical") return MagneticFieldType::vertical;
    if (name == "toroidal") return MagneticFieldType::toroidal;
    if (name == "tangled") return MagneticFieldType::tangled;

    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<problem>/bfield_type must be none, vertical, toroidal, or tangled."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  bool UsesMagneticField(const MagneticFieldType type) {
    return type != MagneticFieldType::none;
  }
}

//===========================================================================//
//                               Initialize                                  //
//===========================================================================//

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;

  // Enroll user functions
  user_srcs_func  = UserSource;
  user_bcs_func   = UserBoundary;
  pgen_final_func = FreeProfile;
  user_ref_func   = RefinementCondition;

  if (global_variable::my_rank == 0) {
    std::cout << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Units                                         " << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Unit Length         : " << 1.0/pmbp->punit->kpc()
                                          << " kpc"   << std::endl;
    std::cout << "Unit Temperature    : " << 1.0/pmbp->punit->kelvin()
                                          << " K"     << std::endl;
    std::cout << "Unit Number Density : " << std::pow(pmbp->punit->cm(),3)
                                          << " cm^-3" << std::endl;
    std::cout << "Unit Velocity       : " << 1.0/pmbp->punit->km_s()
                                          << " km/s"  << std::endl;
    std::cout << "Unit Time           : " << 1.0/pmbp->punit->myr()
                                          << " Myr"   << std::endl;
    std::cout << std::endl;
  }

  // Read in constants
  r_scale   = pin->GetReal("potential", "r_scale");
  rho_scale = pin->GetReal("potential", "rho_scale");
  m_gal     = pin->GetReal("potential", "mass_gal");
  a_gal     = pin->GetReal("potential", "scale_gal");
  z_gal     = pin->GetReal("potential", "z_gal");
  r_200     = pin->GetReal("potential", "r_200");
  rho_mean  = pin->GetReal("potential", "rho_mean");
  r_circ    = pin->GetReal("problem", "r_circ");
  v_circ    = pin->GetReal("problem", "v_circ");
  Z         = pin->GetOrAddReal("problem", "metallicity", 1.0/3);

  // Read in SN injection radius and compute energy and mass injection densities
  r_inj = pin->GetReal("SN","r_inj"); // Input in code unis
  const Real sphere_vol = (4.0/3.0)*M_PI*std::pow(r_inj,3);
  const Real E_def = 1e51; // Default 10^51 ergs
  const Real M_def = 8.4;  // Default 8.4 solar masses
  e_sn  = pin->GetOrAddReal("SN","E_sn",E_def)*pmbp->punit->erg()/sphere_vol;
  m_ej  = pin->GetOrAddReal("SN","M_ej",M_def)*pmbp->punit->msun()/sphere_vol;
  Z_ej  = pin->GetOrAddReal("SN","Z_ej",0.1);
  sn_delay = pin->GetOrAddReal("SN","delay",0.0);

  // Set passive scalar indices
  const bool use_mhd = (pmbp->pmhd != nullptr);
  if (!use_mhd && pmbp->phydro == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Gotham requires either <hydro> or <mhd>." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const int nscalars = use_mhd ? pmbp->pmhd->nscalars : pmbp->phydro->nscalars;
  const int nfluid = use_mhd ? pmbp->pmhd->nmhd : pmbp->phydro->nhydro;
  if (nscalars != 3) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Dust requires the active gas module to have nscalars = 3, but got "
              << nscalars << std::endl;
    std::exit(EXIT_FAILURE);
  }
  scalar_IZS = nfluid;
  scalar_IDS = nfluid + 1;
  scalar_IDL = nfluid + 2;

  // Read in dust model parameters
  d_z_init = pin->GetReal("dust","d_z_init");
  d_z_sn   = pin->GetReal("dust","d_z_sn");
  rho_gr   = pin->GetReal("dust", "rho_gr");
  a_gr_s   = pin->GetReal("dust", "a_gr_s");
  a_gr_l   = pin->GetReal("dust", "a_gr_l");
  min_dust_frac = pin->GetOrAddReal("dust", "D_floor", 1.0e-20);
  dust_donor_cap_frac = pin->GetOrAddReal("dust", "dust_donor_cap_frac",
		                          1.0 - min_dust_frac);
  accretion_nH_min = pin->GetOrAddReal("dust", "accretion_nH_min", 0.1);
  accretion_nH_max = pin->GetOrAddReal("dust", "accretion_nH_max", 1.0e3);
  accretion_T_max = pin->GetOrAddReal("dust", "accretion_T_max", 1.0e4);
  accretion_refractory_frac = pin->GetOrAddReal("dust",
                                                "accretion_refractory_frac",
                                                0.1);

  auto FatalDustInput = [&](const char *msg) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << msg << std::endl;
    std::exit(EXIT_FAILURE);
  };
  if (a_gr_s <= 0.0) FatalDustInput("<dust>/a_gr_s must be > 0.");
  if (a_gr_l <= 0.0) FatalDustInput("<dust>/a_gr_l must be > 0.");
  if (accretion_refractory_frac <= 0.0 || accretion_refractory_frac > 1.0) {
    FatalDustInput("<dust>/accretion_refractory_frac must satisfy 0 < f <= 1.");
  }
  if (accretion_nH_min <= 0.0) {
    FatalDustInput("<dust>/accretion_nH_min must be > 0.");
  }
  if (accretion_nH_max < accretion_nH_min) {
    FatalDustInput("<dust>/accretion_nH_max must be >= <dust>/accretion_nH_min.");
  }
  if (accretion_T_max <= 0.0) {
    FatalDustInput("<dust>/accretion_T_max must be > 0.");
  }

  // Read the density gradient threshold for refinement
  ddens_threshold = pin->GetReal("problem", "ddens_max");

  // Output parameter information
  if (global_variable::my_rank == 0) {
    std::cout << "==============================================" << std::endl;
    std::cout << "Potential Parameters                          " << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "r_scale             : " << r_scale              << std::endl;
    std::cout << "rho_scale           : " << rho_scale            << std::endl;
    std::cout << "m_gal               : " << m_gal                << std::endl;
    std::cout << "a_gal               : " << a_gal                << std::endl;
    std::cout << "z_gal               : " << z_gal                << std::endl;
    std::cout << "r_200               : " << r_200                << std::endl;
    std::cout << "rho_mean            : " << rho_mean             << std::endl;
    std::cout << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Other Parameters                              " << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "r_circ              : " << r_circ               << std::endl;
    std::cout << "v_circ              : " << v_circ               << std::endl;
    std::cout << "metallicity         : " << Z                    << std::endl;
    std::cout << "ddens_threshold     : " << ddens_threshold      << std::endl;
    std::cout << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Supernova Parameters                          " << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "r_inj               : " << r_inj                << std::endl;
    std::cout << "e_sn                : " << e_sn                 << std::endl;
    std::cout << "m_ej                : " << m_ej                 << std::endl;
    std::cout << "Z_ej                : " << Z_ej                 << std::endl;
    std::cout << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Dust Parameters                               " << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "d_z_init            : " << d_z_init             << std::endl;
    std::cout << "d_z_sn              : " << d_z_sn               << std::endl;
    std::cout << "rho_gr              : " << rho_gr               << std::endl;
    std::cout << "a_gr_s              : " << a_gr_s               << std::endl;
    std::cout << "a_gr_l              : " << a_gr_l               << std::endl;
    std::cout << "D_floor             : " << min_dust_frac        << std::endl;
    std::cout << "dust_donor_cap_frac : " << dust_donor_cap_frac  << std::endl;
    std::cout << "accretion_nH_min    : " << accretion_nH_min     << std::endl;
    std::cout << "accretion_nH_max    : " << accretion_nH_max     << std::endl;
    std::cout << "accretion_T_max     : " << accretion_T_max      << std::endl;
    std::cout << "accretion_f_ref     : " << accretion_refractory_frac
              << std::endl;
    std::cout << std::endl;
  }

  // Read the CGM cooling flow profile file
  std::string profile_file = pin->GetString("problem", "profile_file");
  ProfileReaderHost profile_reader_host;
  profile_reader_host.ReadProfiles(profile_file);
  profile_reader = profile_reader_host.CreateDeviceReader();
  if (global_variable::my_rank==0) {
    std::cout << "Successfully loaded CGM profiles from "
              << profile_file << std::endl;
  }

  // Read in the disk profile file
  std::string disk_profile_file = pin->GetString("problem", "disk_profile_file");
  ProfileReaderHost disk_profile_reader_host;
  disk_profile_reader_host.ReadProfiles(disk_profile_file);
  disk_profile_reader = disk_profile_reader_host.CreateDeviceReader();
  if (global_variable::my_rank==0) {
    std::cout << "Successfully loaded disk profiles from "
	      << disk_profile_file << std::endl;
  }

  std::string bfield_type = pin->GetOrAddString("problem", "bfield_type", "none");
  MagneticFieldType bfield = ParseMagneticFieldType(bfield_type);
  Real b_ref = 0.0;
  if (UsesMagneticField(bfield)) {
    if (pmbp->pmhd == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<problem>/bfield_type != none requires an <mhd> block."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (!pin->DoesParameterExist("problem", "beta")) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<problem>/beta is required when bfield_type != none."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    const Real beta = pin->GetReal("problem", "beta");
    const Real beta_r_ref = pin->GetOrAddReal("problem", "beta_r_ref", 250.0);
    if (beta <= 0.0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<problem>/beta must be > 0." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    b_ref = ReferenceMagneticField(profile_reader, beta_r_ref, beta);
    if (global_variable::my_rank == 0) {
      std::cout << "==============================================" << std::endl;
      std::cout << "Magnetic Field Parameters                     " << std::endl;
      std::cout << "==============================================" << std::endl;
      std::cout << "bfield_type         : " << bfield_type          << std::endl;
      std::cout << "beta                : " << beta                 << std::endl;
      std::cout << "beta_r_ref          : " << beta_r_ref           << std::endl;
      std::cout << "B_ref               : " << b_ref                << std::endl;
      std::cout << std::endl;
    }
  }

  // Count total particles and initialize SN centers buffer
  pmy_mesh_->CountParticles();
  sn_centers_buffer = DvceArray2D<Real>("sn_centers_buffer", 3, pmy_mesh_->nprtcl_total);
  sn_counter = Kokkos::View<int>("sn_counter");
  if (global_variable::my_rank==0) {
    std::cout << "Successfully initialized " << pmy_mesh_->nprtcl_total
              << " particles!" << std::endl;
  }

  if (restart) return;

  // Generate the initial turbulent field
  int nlow   = pin->GetOrAddInteger("problem", "cgm_turb_nlow", 1);
  int nhigh  = pin->GetOrAddInteger("problem", "cgm_turb_nhigh", 8);
  Real expo  = pin->GetOrAddReal("problem", "cgm_turb_expo", 5.0/3.0);
  Real v_rms = pin->GetOrAddReal("problem", "cgm_turb_rms", 0.1);
  Real cgm_turb_xscale = pin->GetOrAddReal("problem", "cgm_turb_xscale", 0.01);
  Real cgm_turb_yscale = pin->GetOrAddReal("problem", "cgm_turb_yscale", 0.01);
  Real cgm_turb_zscale = pin->GetOrAddReal("problem", "cgm_turb_zscale", 0.01);

  // Initialize random state
  RNG_State rstate;
  rstate.idum = -1;

  // Domain size
  Real lx = pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min;
  Real ly = pmy_mesh_->mesh_size.x2max - pmy_mesh_->mesh_size.x2min;
  Real lz = pmy_mesh_->mesh_size.x3max - pmy_mesh_->mesh_size.x3min;
  Real dkx = 2.0*M_PI/lx;
  Real dky = 2.0*M_PI/ly;
  Real dkz = 2.0*M_PI/lz;

  // Count modes
  int nmodes = 0;
  for (int nkx = -nhigh; nkx <= nhigh; nkx++) {
    for (int nky = -nhigh; nky <= nhigh; nky++) {
      for (int nkz = -nhigh; nkz <= nhigh; nkz++) {
        if (nkx == 0 && nky == 0 && nkz == 0) continue;
        int nsqr = nkx*nkx + nky*nky + nkz*nkz;
        if (nsqr >= nlow*nlow && nsqr <= nhigh*nhigh) {
          nmodes++;
        }
  }}}

  // Allocate arrays
  DualArray2D<Real> k_modes, aka, akb;
  Kokkos::realloc(k_modes, 3, nmodes);
  Kokkos::realloc(aka, 3, nmodes);
  Kokkos::realloc(akb, 3, nmodes);

  // Generate modes
  int nmode = 0;
  Real total_energy = 0.0;
  for (int nkx = -nhigh; nkx <= nhigh; nkx++) {
    for (int nky = -nhigh; nky <= nhigh; nky++) {
      for (int nkz = -nhigh; nkz <= nhigh; nkz++) {
        if (nkx == 0 && nky == 0 && nkz == 0) continue;
        int nsqr = nkx*nkx + nky*nky + nkz*nkz;
        if (nsqr >= nlow*nlow && nsqr <= nhigh*nhigh) {
          Real kx = dkx*nkx, ky = dky*nky, kz = dkz*nkz;
          Real kiso = sqrt(kx*kx + ky*ky + kz*kz);

          k_modes.h_view(0, nmode) = kx;
          k_modes.h_view(1, nmode) = ky;
          k_modes.h_view(2, nmode) = kz;

          Real norm = 1.0/pow(kiso, (expo+2.0)/2.0);

          Real aval[3], bval[3];
          for (int dir = 0; dir < 3; dir++) {
            aval[dir] = norm * RanGaussianSt(&rstate);
            bval[dir] = norm * RanGaussianSt(&rstate);
          }

          Real k_dirs[3] = {kx, ky, kz};
          Real ka = kx*aval[0] + ky*aval[1] + kz*aval[2];
          Real kb = kx*bval[0] + ky*bval[1] + kz*bval[2];

          for (int dir = 0; dir < 3; dir++) {
            aval[dir] -= k_dirs[dir]*ka/(kiso*kiso);
            bval[dir] -= k_dirs[dir]*kb/(kiso*kiso);

            aka.h_view(dir,nmode) = aval[dir];
            akb.h_view(dir,nmode) = bval[dir];

            total_energy += 0.5*(aval[dir]*aval[dir] + bval[dir]*bval[dir]);
          }
          nmode++;
        }
      }
    }
  }

  Real v_norm = v_rms/sqrt(total_energy);

  k_modes.template modify<HostMemSpace>();
  k_modes.template sync<DevExeSpace>();
  aka.template modify<HostMemSpace>();
  aka.template sync<DevExeSpace>();
  akb.template modify<HostMemSpace>();
  akb.template sync<DevExeSpace>();

  // Capture variables for kernel
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  auto &size = pmbp->pmb->mb_size;

  int IZS = scalar_IZS;
  int IDS = scalar_IDS;
  int IDL = scalar_IDL;

  auto &u0 = (pmbp->pmhd != nullptr) ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  EOS_Data &eos = (pmbp->pmhd != nullptr) ? pmbp->pmhd->peos->eos_data :
                                            pmbp->phydro->peos->eos_data;
  Real gm1 = eos.gamma - 1.0;

  auto &profile = profile_reader;
  auto &disk_profile = disk_profile_reader;

  Real G = pmbp->punit->grav_constant();
  Real r_s = r_scale;
  Real rho_s = rho_scale;
  Real m_g = m_gal;
  Real a_g = a_gal;
  Real z_g = z_gal;
  Real r_m = r_200;
  Real rho_m = rho_mean;

  Real r_c = r_circ;
  Real v_c = v_circ;

  Real Zsol = Z_sol;
  Real Z_ = Z;
  Real dz_init = d_z_init;
  Real min_df = min_dust_frac;

  // Use loaded profiles
  par_for("pgen_ic", DevExeSpace(), 0, (pmbp->nmb_thispack-1), ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    int nx1 = indcs.nx1;
    Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
    Real x1l =   LeftEdgeX(i-is, nx1, x1min, x1max);
    Real x1r = LeftEdgeX(i+1-is, nx1, x1min, x1max);

    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    int nx2 = indcs.nx2;
    Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
    Real x2l =   LeftEdgeX(j-js, nx2, x2min, x2max);
    Real x2r = LeftEdgeX(j+1-js, nx2, x2min, x2max);

    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    int nx3 = indcs.nx3;
    Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);

    SetCoolingFlowState(u0, m, k, j, i, x1v, x2v, x3v, gm1, profile);
    SetRotation(u0, m, k, j, i, x1v, x2v, x3v, r_c, v_c);
    SetDustScalars(u0, m, k, j, i, IZS, IDS, IDL, min_df, Z_, Zsol);
    SetEquilibriumState(u0, m, k, j, i, x1v, x2v, x3v,
                        x1l, x1r, x2l, x2r, G, r_s, rho_s,
                        m_g, a_g, z_g, r_m, rho_m, gm1,
                        IZS, IDS, IDL, dz_init, Z_, Zsol,
                        disk_profile);

    // Compute turbulent velocities by summing Fourier modes
    Real vx = 0.0, vy = 0.0, vz = 0.0;
    for (int n = 0; n < nmodes; n++) {
      Real phase = k_modes.d_view(0,n)*x1v
	         + k_modes.d_view(1,n)*x2v
	         + k_modes.d_view(2,n)*x3v;
      Real cos_phase = cos(phase);
      Real sin_phase = sin(phase);

      vx += aka.d_view(0,n)*cos_phase - akb.d_view(0,n)*sin_phase;
      vy += aka.d_view(1,n)*cos_phase - akb.d_view(1,n)*sin_phase;
      vz += aka.d_view(2,n)*cos_phase - akb.d_view(2,n)*sin_phase;
    }

    // Attenuate in the center by 1 - Gaussian
    Real att = 1.0 - exp(-0.5 * ( SQR(x1v)/SQR(cgm_turb_xscale)
                                + SQR(x2v)/SQR(cgm_turb_yscale)
                                + SQR(x3v)/SQR(cgm_turb_zscale)));

    // Normalize to desired RMS velocity
    vx *= v_norm*att; vy *= v_norm*att; vz *= v_norm*att;

    // Add to conserved variables
    Real rho = u0(m,IDN,k,j,i);
    Real rho_v1 = u0(m,IM1,k,j,i);
    Real rho_v2 = u0(m,IM2,k,j,i);
    Real rho_v3 = u0(m,IM3,k,j,i);

    u0(m,IEN,k,j,i) += 0.5 * rho * (SQR(vx) + SQR(vy) + SQR(vz));
    u0(m,IEN,k,j,i) += rho_v1 * vx + rho_v2 * vy + rho_v3 * vz;
    u0(m,IM1,k,j,i) += rho * vx;
    u0(m,IM2,k,j,i) += rho * vy;
    u0(m,IM3,k,j,i) += rho * vz;
  });

  if (global_variable::my_rank == 0) {
    std::cout << "Successfully initialized grid!" << std::endl;
  }

  if (pmbp->pmhd != nullptr) {
    InitializeMagneticField(pin, pmbp, bfield_type, b_ref);
  }

  return;
}

Real ReferenceMagneticField(const ProfileReader &profile, const Real r_ref,
                            const Real beta) {
  DvceArray1D<Real> b_ref("gotham_b_ref", 1);
  Kokkos::parallel_for("gotham_b_ref", Kokkos::RangePolicy<>(DevExeSpace(), 0, 1),
  KOKKOS_LAMBDA(const int n) {
    const Real pressure = profile.GetDensity(r_ref)*profile.GetTemperature(r_ref);
    b_ref(n) = sqrt(2.0*pressure/beta);
  });
  DevExeSpace().fence();

  auto h_b_ref = Kokkos::create_mirror_view(b_ref);
  Kokkos::deep_copy(h_b_ref, b_ref);
  return h_b_ref(0);
}

KOKKOS_INLINE_FUNCTION
Real ToroidalVectorPotential(const Real x1, const Real x2, const Real b0,
                             const Real r_soft) {
  return -b0*sqrt(x1*x1 + x2*x2 + r_soft*r_soft);
}

KOKKOS_INLINE_FUNCTION
Real TangledVectorPotential(const int dir, const Real x1, const Real x2,
                            const Real x3, const int nmodes,
                            const DvceArray2D<Real> &k_modes,
                            const DvceArray2D<Real> &aka,
                            const DvceArray2D<Real> &akb) {
  Real avec = 0.0;
  for (int n = 0; n < nmodes; ++n) {
    const Real phase = k_modes(0,n)*x1 + k_modes(1,n)*x2 + k_modes(2,n)*x3;
    avec += aka(dir,n)*cos(phase) - akb(dir,n)*sin(phase);
  }
  return avec;
}

void InitializeMagneticField(ParameterInput *pin, MeshBlockPack *pmbp,
                             const std::string &bfield_type,
                             const Real b_ref) {
  auto bfield = ParseMagneticFieldType(bfield_type);
  auto &indcs = pmbp->pmesh->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int nmb = pmbp->nmb_thispack;
  auto &size = pmbp->pmb->mb_size;
  auto &b0 = pmbp->pmhd->b0;
  auto &bcc0 = pmbp->pmhd->bcc0;
  auto &u0 = pmbp->pmhd->u0;

  const bool use_tangled = (bfield == MagneticFieldType::tangled);
  const bool use_toroidal = (bfield == MagneticFieldType::toroidal);
  const bool use_vertical = (bfield == MagneticFieldType::vertical);

  if (bfield == MagneticFieldType::none || use_vertical) {
    par_for("gotham_b_vertical", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      b0.x1f(m,k,j,i) = 0.0;
      b0.x2f(m,k,j,i) = 0.0;
      b0.x3f(m,k,j,i) = use_vertical ? b_ref : 0.0;
      if (i == ie) b0.x1f(m,k,j,i+1) = 0.0;
      if (j == je) b0.x2f(m,k,j+1,i) = 0.0;
      if (k == ke) b0.x3f(m,k+1,j,i) = use_vertical ? b_ref : 0.0;
    });
  } else {
    int ncells1 = indcs.nx1 + 2*(indcs.ng);
    int ncells2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*(indcs.ng)) : 2;
    int ncells3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*(indcs.ng)) : 2;
    DvceArray4D<Real> a1, a2, a3;
    Kokkos::realloc(a1, nmb, ncells3, ncells2, ncells1);
    Kokkos::realloc(a2, nmb, ncells3, ncells2, ncells1);
    Kokkos::realloc(a3, nmb, ncells3, ncells2, ncells1);

    int nmodes = 0;
    DualArray2D<Real> k_modes, aka, akb;
    if (use_tangled) {
      const int nlow = pin->GetOrAddInteger("problem", "bfield_tangled_nlow", 1);
      const int nhigh = pin->GetOrAddInteger("problem", "bfield_tangled_nhigh", 8);
      const Real expo = pin->GetOrAddReal("problem", "bfield_tangled_expo", 5.0/3.0);
      const int seed = pin->GetOrAddInteger("problem", "bfield_tangled_seed", -1);
      const Real lx = pmbp->pmesh->mesh_size.x1max - pmbp->pmesh->mesh_size.x1min;
      const Real ly = pmbp->pmesh->mesh_size.x2max - pmbp->pmesh->mesh_size.x2min;
      const Real lz = pmbp->pmesh->mesh_size.x3max - pmbp->pmesh->mesh_size.x3min;
      const Real dkx = 2.0*M_PI/lx;
      const Real dky = 2.0*M_PI/ly;
      const Real dkz = 2.0*M_PI/lz;

      for (int nkx = -nhigh; nkx <= nhigh; ++nkx) {
        for (int nky = -nhigh; nky <= nhigh; ++nky) {
          for (int nkz = -nhigh; nkz <= nhigh; ++nkz) {
            if (nkx == 0 && nky == 0 && nkz == 0) continue;
            const int nsqr = nkx*nkx + nky*nky + nkz*nkz;
            if (nsqr >= nlow*nlow && nsqr <= nhigh*nhigh) ++nmodes;
          }
        }
      }
      if (nmodes == 0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Tangled magnetic field requested with zero Fourier modes."
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }

      Kokkos::realloc(k_modes, 3, nmodes);
      Kokkos::realloc(aka, 3, nmodes);
      Kokkos::realloc(akb, 3, nmodes);

      RNG_State rstate;
      rstate.idum = (seed > 0) ? -seed : seed;
      if (rstate.idum == 0) rstate.idum = -1;
      int nmode = 0;
      for (int nkx = -nhigh; nkx <= nhigh; ++nkx) {
        for (int nky = -nhigh; nky <= nhigh; ++nky) {
          for (int nkz = -nhigh; nkz <= nhigh; ++nkz) {
            if (nkx == 0 && nky == 0 && nkz == 0) continue;
            const int nsqr = nkx*nkx + nky*nky + nkz*nkz;
            if (nsqr < nlow*nlow || nsqr > nhigh*nhigh) continue;

            const Real kx = dkx*nkx;
            const Real ky = dky*nky;
            const Real kz = dkz*nkz;
            const Real kiso = sqrt(kx*kx + ky*ky + kz*kz);
            const Real norm = 1.0/pow(kiso, 0.5*(expo + 4.0));

            k_modes.h_view(0,nmode) = kx;
            k_modes.h_view(1,nmode) = ky;
            k_modes.h_view(2,nmode) = kz;
            for (int dir = 0; dir < 3; ++dir) {
              aka.h_view(dir,nmode) = norm*RanGaussianSt(&rstate);
              akb.h_view(dir,nmode) = norm*RanGaussianSt(&rstate);
            }
            ++nmode;
          }
        }
      }
      k_modes.template modify<HostMemSpace>();
      k_modes.template sync<DevExeSpace>();
      aka.template modify<HostMemSpace>();
      aka.template sync<DevExeSpace>();
      akb.template modify<HostMemSpace>();
      akb.template sync<DevExeSpace>();
    }

    auto k_modes_d = k_modes.d_view;
    auto aka_d = aka.d_view;
    auto akb_d = akb.d_view;
    auto &nghbr = pmbp->pmb->nghbr;
    auto &mblev = pmbp->pmb->mb_lev;
    const Real r_soft = pin->GetOrAddReal("problem", "bfield_toroidal_r_soft", 0.0);

    par_for("gotham_vector_potential", DevExeSpace(), 0, nmb-1, ks, ke+1, js, je+1,
            is, ie+1, KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      Real &x2min = size.d_view(m).x2min;
      Real &x2max = size.d_view(m).x2max;
      Real &x3min = size.d_view(m).x3min;
      Real &x3max = size.d_view(m).x3max;

      const int nx1 = indcs.nx1;
      const int nx2 = indcs.nx2;
      const int nx3 = indcs.nx3;
      const Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
      const Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
      const Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
      const Real x1f = LeftEdgeX(i-is, nx1, x1min, x1max);
      const Real x2f = LeftEdgeX(j-js, nx2, x2min, x2max);
      const Real x3f = LeftEdgeX(k-ks, nx3, x3min, x3max);
      const Real dx1 = size.d_view(m).dx1;
      const Real dx2 = size.d_view(m).dx2;
      const Real dx3 = size.d_view(m).dx3;

      a1(m,k,j,i) = use_tangled ?
          TangledVectorPotential(0, x1v, x2f, x3f, nmodes, k_modes_d, aka_d, akb_d) : 0.0;
      a2(m,k,j,i) = use_tangled ?
          TangledVectorPotential(1, x1f, x2v, x3f, nmodes, k_modes_d, aka_d, akb_d) : 0.0;
      a3(m,k,j,i) = use_toroidal ?
          ToroidalVectorPotential(x1f, x2f, b_ref, r_soft) :
          TangledVectorPotential(2, x1f, x2f, x3v, nmodes, k_modes_d, aka_d, akb_d);

      if (use_tangled &&
          ((nghbr.d_view(m,8 ).lev > mblev.d_view(m) && j==js) ||
           (nghbr.d_view(m,9 ).lev > mblev.d_view(m) && j==js) ||
           (nghbr.d_view(m,10).lev > mblev.d_view(m) && j==js) ||
           (nghbr.d_view(m,11).lev > mblev.d_view(m) && j==js) ||
           (nghbr.d_view(m,12).lev > mblev.d_view(m) && j==je+1) ||
           (nghbr.d_view(m,13).lev > mblev.d_view(m) && j==je+1) ||
           (nghbr.d_view(m,14).lev > mblev.d_view(m) && j==je+1) ||
           (nghbr.d_view(m,15).lev > mblev.d_view(m) && j==je+1) ||
           (nghbr.d_view(m,24).lev > mblev.d_view(m) && k==ks) ||
           (nghbr.d_view(m,25).lev > mblev.d_view(m) && k==ks) ||
           (nghbr.d_view(m,26).lev > mblev.d_view(m) && k==ks) ||
           (nghbr.d_view(m,27).lev > mblev.d_view(m) && k==ks) ||
           (nghbr.d_view(m,28).lev > mblev.d_view(m) && k==ke+1) ||
           (nghbr.d_view(m,29).lev > mblev.d_view(m) && k==ke+1) ||
           (nghbr.d_view(m,30).lev > mblev.d_view(m) && k==ke+1) ||
           (nghbr.d_view(m,31).lev > mblev.d_view(m) && k==ke+1) ||
           (nghbr.d_view(m,40).lev > mblev.d_view(m) && j==js && k==ks) ||
           (nghbr.d_view(m,41).lev > mblev.d_view(m) && j==js && k==ks) ||
           (nghbr.d_view(m,42).lev > mblev.d_view(m) && j==je+1 && k==ks) ||
           (nghbr.d_view(m,43).lev > mblev.d_view(m) && j==je+1 && k==ks) ||
           (nghbr.d_view(m,44).lev > mblev.d_view(m) && j==js && k==ke+1) ||
           (nghbr.d_view(m,45).lev > mblev.d_view(m) && j==js && k==ke+1) ||
           (nghbr.d_view(m,46).lev > mblev.d_view(m) && j==je+1 && k==ke+1) ||
           (nghbr.d_view(m,47).lev > mblev.d_view(m) && j==je+1 && k==ke+1))) {
        const Real xl = x1v + 0.25*dx1;
        const Real xr = x1v - 0.25*dx1;
        a1(m,k,j,i) = 0.5*(
            TangledVectorPotential(0, xl, x2f, x3f, nmodes, k_modes_d, aka_d, akb_d) +
            TangledVectorPotential(0, xr, x2f, x3f, nmodes, k_modes_d, aka_d, akb_d));
      }

      if (use_tangled &&
          ((nghbr.d_view(m,0 ).lev > mblev.d_view(m) && i==is) ||
           (nghbr.d_view(m,1 ).lev > mblev.d_view(m) && i==is) ||
           (nghbr.d_view(m,2 ).lev > mblev.d_view(m) && i==is) ||
           (nghbr.d_view(m,3 ).lev > mblev.d_view(m) && i==is) ||
           (nghbr.d_view(m,4 ).lev > mblev.d_view(m) && i==ie+1) ||
           (nghbr.d_view(m,5 ).lev > mblev.d_view(m) && i==ie+1) ||
           (nghbr.d_view(m,6 ).lev > mblev.d_view(m) && i==ie+1) ||
           (nghbr.d_view(m,7 ).lev > mblev.d_view(m) && i==ie+1) ||
           (nghbr.d_view(m,24).lev > mblev.d_view(m) && k==ks) ||
           (nghbr.d_view(m,25).lev > mblev.d_view(m) && k==ks) ||
           (nghbr.d_view(m,26).lev > mblev.d_view(m) && k==ks) ||
           (nghbr.d_view(m,27).lev > mblev.d_view(m) && k==ks) ||
           (nghbr.d_view(m,28).lev > mblev.d_view(m) && k==ke+1) ||
           (nghbr.d_view(m,29).lev > mblev.d_view(m) && k==ke+1) ||
           (nghbr.d_view(m,30).lev > mblev.d_view(m) && k==ke+1) ||
           (nghbr.d_view(m,31).lev > mblev.d_view(m) && k==ke+1) ||
           (nghbr.d_view(m,32).lev > mblev.d_view(m) && i==is && k==ks) ||
           (nghbr.d_view(m,33).lev > mblev.d_view(m) && i==is && k==ks) ||
           (nghbr.d_view(m,34).lev > mblev.d_view(m) && i==ie+1 && k==ks) ||
           (nghbr.d_view(m,35).lev > mblev.d_view(m) && i==ie+1 && k==ks) ||
           (nghbr.d_view(m,36).lev > mblev.d_view(m) && i==is && k==ke+1) ||
           (nghbr.d_view(m,37).lev > mblev.d_view(m) && i==is && k==ke+1) ||
           (nghbr.d_view(m,38).lev > mblev.d_view(m) && i==ie+1 && k==ke+1) ||
           (nghbr.d_view(m,39).lev > mblev.d_view(m) && i==ie+1 && k==ke+1))) {
        const Real xl = x2v + 0.25*dx2;
        const Real xr = x2v - 0.25*dx2;
        a2(m,k,j,i) = 0.5*(
            TangledVectorPotential(1, x1f, xl, x3f, nmodes, k_modes_d, aka_d, akb_d) +
            TangledVectorPotential(1, x1f, xr, x3f, nmodes, k_modes_d, aka_d, akb_d));
      }

      if ((nghbr.d_view(m,0 ).lev > mblev.d_view(m) && i==is) ||
          (nghbr.d_view(m,1 ).lev > mblev.d_view(m) && i==is) ||
          (nghbr.d_view(m,2 ).lev > mblev.d_view(m) && i==is) ||
          (nghbr.d_view(m,3 ).lev > mblev.d_view(m) && i==is) ||
          (nghbr.d_view(m,4 ).lev > mblev.d_view(m) && i==ie+1) ||
          (nghbr.d_view(m,5 ).lev > mblev.d_view(m) && i==ie+1) ||
          (nghbr.d_view(m,6 ).lev > mblev.d_view(m) && i==ie+1) ||
          (nghbr.d_view(m,7 ).lev > mblev.d_view(m) && i==ie+1) ||
          (nghbr.d_view(m,8 ).lev > mblev.d_view(m) && j==js) ||
          (nghbr.d_view(m,9 ).lev > mblev.d_view(m) && j==js) ||
          (nghbr.d_view(m,10).lev > mblev.d_view(m) && j==js) ||
          (nghbr.d_view(m,11).lev > mblev.d_view(m) && j==js) ||
          (nghbr.d_view(m,12).lev > mblev.d_view(m) && j==je+1) ||
          (nghbr.d_view(m,13).lev > mblev.d_view(m) && j==je+1) ||
          (nghbr.d_view(m,14).lev > mblev.d_view(m) && j==je+1) ||
          (nghbr.d_view(m,15).lev > mblev.d_view(m) && j==je+1) ||
          (nghbr.d_view(m,16).lev > mblev.d_view(m) && i==is && j==js) ||
          (nghbr.d_view(m,17).lev > mblev.d_view(m) && i==is && j==js) ||
          (nghbr.d_view(m,18).lev > mblev.d_view(m) && i==ie+1 && j==js) ||
          (nghbr.d_view(m,19).lev > mblev.d_view(m) && i==ie+1 && j==js) ||
          (nghbr.d_view(m,20).lev > mblev.d_view(m) && i==is && j==je+1) ||
          (nghbr.d_view(m,21).lev > mblev.d_view(m) && i==is && j==je+1) ||
          (nghbr.d_view(m,22).lev > mblev.d_view(m) && i==ie+1 && j==je+1) ||
          (nghbr.d_view(m,23).lev > mblev.d_view(m) && i==ie+1 && j==je+1)) {
        const Real xl = x3v + 0.25*dx3;
        const Real xr = x3v - 0.25*dx3;
        a3(m,k,j,i) = use_toroidal ?
            ToroidalVectorPotential(x1f, x2f, b_ref, r_soft) :
            0.5*(TangledVectorPotential(2, x1f, x2f, xl, nmodes, k_modes_d, aka_d, akb_d) +
                 TangledVectorPotential(2, x1f, x2f, xr, nmodes, k_modes_d, aka_d, akb_d));
      }
    });

    par_for("gotham_b_curl", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      const Real dx1 = size.d_view(m).dx1;
      const Real dx2 = size.d_view(m).dx2;
      const Real dx3 = size.d_view(m).dx3;

      b0.x1f(m,k,j,i) = ((a3(m,k,j+1,i) - a3(m,k,j,i))/dx2 -
                         (a2(m,k+1,j,i) - a2(m,k,j,i))/dx3);
      b0.x2f(m,k,j,i) = ((a1(m,k+1,j,i) - a1(m,k,j,i))/dx3 -
                         (a3(m,k,j,i+1) - a3(m,k,j,i))/dx1);
      b0.x3f(m,k,j,i) = ((a2(m,k,j,i+1) - a2(m,k,j,i))/dx1 -
                         (a1(m,k,j+1,i) - a1(m,k,j,i))/dx2);

      if (i == ie) {
        b0.x1f(m,k,j,i+1) = ((a3(m,k,j+1,i+1) - a3(m,k,j,i+1))/dx2 -
                             (a2(m,k+1,j,i+1) - a2(m,k,j,i+1))/dx3);
      }
      if (j == je) {
        b0.x2f(m,k,j+1,i) = ((a1(m,k+1,j+1,i) - a1(m,k,j+1,i))/dx3 -
                             (a3(m,k,j+1,i+1) - a3(m,k,j+1,i))/dx1);
      }
      if (k == ke) {
        b0.x3f(m,k+1,j,i) = ((a2(m,k+1,j,i+1) - a2(m,k+1,j,i))/dx1 -
                             (a1(m,k+1,j+1,i) - a1(m,k+1,j,i))/dx2);
      }
    });

    if (use_tangled) {
      const int nx1 = indcs.nx1;
      const int nx2 = indcs.nx2;
      const int nx3 = indcs.nx3;
      const int nmkji = nmb*nx3*nx2*nx1;
      const int nkji = nx3*nx2*nx1;
      const int nji = nx2*nx1;
      Real b2_sum = 0.0;
      Real vol_sum = 0.0;
      Kokkos::parallel_reduce("gotham_b_rms", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
      KOKKOS_LAMBDA(const int idx, Real &b2_total, Real &vol_total) {
        int m = idx/nkji;
        int k = (idx - m*nkji)/nji;
        int j = (idx - m*nkji - k*nji)/nx1;
        int i = (idx - m*nkji - k*nji - j*nx1) + is;
        k += ks;
        j += js;
        const Real bx = 0.5*(b0.x1f(m,k,j,i) + b0.x1f(m,k,j,i+1));
        const Real by = 0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i));
        const Real bz = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
        const Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
        b2_total += (bx*bx + by*by + bz*bz)*vol;
        vol_total += vol;
      }, Kokkos::Sum<Real>(b2_sum), Kokkos::Sum<Real>(vol_sum));

#if MPI_PARALLEL_ENABLED
      Real local[2] = {b2_sum, vol_sum};
      Real global[2];
      MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      b2_sum = global[0];
      vol_sum = global[1];
#endif

      const Real b_rms = sqrt(b2_sum/vol_sum);
      const Real b_scale = (b_rms > 0.0) ? b_ref/b_rms : 0.0;
      par_for("gotham_b_scale", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        b0.x1f(m,k,j,i) *= b_scale;
        b0.x2f(m,k,j,i) *= b_scale;
        b0.x3f(m,k,j,i) *= b_scale;
        if (i == ie) b0.x1f(m,k,j,i+1) *= b_scale;
        if (j == je) b0.x2f(m,k,j+1,i) *= b_scale;
        if (k == ke) b0.x3f(m,k+1,j,i) *= b_scale;
      });
    }
  }

  par_for("gotham_bcc_energy", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real bx = 0.5*(b0.x1f(m,k,j,i) + b0.x1f(m,k,j,i+1));
    const Real by = 0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i));
    const Real bz = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
    bcc0(m,IBX,k,j,i) = bx;
    bcc0(m,IBY,k,j,i) = by;
    bcc0(m,IBZ,k,j,i) = bz;
    u0(m,IEN,k,j,i) += 0.5*(bx*bx + by*by + bz*bz);
  });

  if (global_variable::my_rank == 0) {
    std::cout << "Successfully initialized magnetic field!" << std::endl;
  }
}

KOKKOS_INLINE_FUNCTION
void SetCoolingFlowState(const DvceArray5D<Real> &u0,
                         int m, int k, int j, int i,
                         Real x1v, Real x2v, Real x3v,
                         Real gm1, const ProfileReader &profile) {
    // Calculate radius
    Real r = sqrt(x1v*x1v + x2v*x2v + x3v*x3v);

    // Get values from profiles via interpolation
    Real rho  = profile.GetDensity(r);
    Real temp = profile.GetTemperature(r);
    Real vr   = profile.GetVelocity(r);
    Real rmin = profile.GetRmin();

    // set vr to go to zero if r < rmin
    if (r < rmin) {
      vr *= (r / rmin);
    }

    // Calculate pressure from temperature
    Real press = rho * temp;

    // Set radial velocity components based on position
    Real v1 = 0.0, v2 = 0.0, v3 = 0.0;
    constexpr Real tiny = 1.0e-20;
    if (r > tiny) {  // Avoid division by zero
      // Negative sign accounts for inflowing vr
      v1 = -vr * x1v / r;
      v2 = -vr * x2v / r;
      v3 = -vr * x3v / r;
    }

    // Set state variables
    u0(m, IDN, k, j, i) = rho;
    u0(m, IM1, k, j, i) = rho * v1;
    u0(m, IM2, k, j, i) = rho * v2;
    u0(m, IM3, k, j, i) = rho * v3;
    u0(m, IEN, k, j, i) = press/gm1 + 0.5*rho*(SQR(v1) + SQR(v2) + SQR(v3));
}

KOKKOS_INLINE_FUNCTION
void SetRotation(const DvceArray5D<Real> &u0,
                 int m, int k, int j, int i,
                 Real x1v, Real x2v, Real x3v,
                 Real r_circ, Real v_circ) {
  // Calculate radius
  Real r = sqrt(x1v*x1v + x2v*x2v + x3v*x3v);
  Real R = sqrt(x1v*x1v + x2v*x2v);

  // Calculate azimuthal velocity
  Real vx = 0.0, vy = 0.0, vz = 0.0;
  constexpr Real tiny = 1.0e-20;
  if (r > tiny && R > tiny) {  // Avoid division by zero
    Real v_phi = 0.0;
    Real sin_theta = R / r;

    if (r < r_circ) {
      v_phi = v_circ * sin_theta;
    }
    if (r > r_circ) {
      v_phi = v_circ * sin_theta * r_circ / r;
    }

    // Calculate azimuthal velocity components
    vx = -v_phi * x2v / R;
    vy = v_phi * x1v / R;
  }

  // Set state variables
  Real rho = u0(m,IDN,k,j,i);
  Real rho_v1 = u0(m,IM1,k,j,i);
  Real rho_v2 = u0(m,IM2,k,j,i);
  Real rho_v3 = u0(m,IM3,k,j,i);

  u0(m,IEN,k,j,i) += 0.5 * rho * (SQR(vx) + SQR(vy) + SQR(vz));
  u0(m,IEN,k,j,i) += rho_v1 * vx + rho_v2 * vy + rho_v3 * vz;
  u0(m,IM1,k,j,i) += rho * vx;
  u0(m,IM2,k,j,i) += rho * vy;
  u0(m,IM3,k,j,i) += rho * vz;
}

KOKKOS_INLINE_FUNCTION
void SetDustScalars(const DvceArray5D<Real> &u0,
                    int m, int k, int j, int i,
                    int IZS, int IDS, int IDL,
                    Real dz, Real Z, Real Zsol) {
  Real rho = u0(m, IDN, k, j, i);
  Real Z_total = Z * Zsol * rho;
  u0(m, IZS, k, j, i) = (1.0 - dz) * Z_total;
  u0(m, IDS, k, j, i) = 0.5 * dz * Z_total;
  u0(m, IDL, k, j, i) = 0.5 * dz * Z_total;
}

KOKKOS_INLINE_FUNCTION
void SetEquilibriumState(const DvceArray5D<Real> &u0,
                         int m, int k, int j, int i,
                         Real x1v, Real x2v, Real x3v,
                         Real x1l, Real x1r, Real x2l, Real x2r,
                         Real G, Real r_s, Real rho_s, Real m_g,
                         Real a_g, Real z_g, Real r_m, Real rho_m,
			 Real gm1, int IZS, int IDS, int IDL,
			 Real dz, Real Z, Real Zsol,
			 const ProfileReader &disk_profile) {
    // Calculate radius
    Real R = sqrt(x1v * x1v + x2v * x2v);
    Real R1l = sqrt(x1l * x1l + x2v * x2v);
    Real R1r = sqrt(x1r * x1r + x2v * x2v);
    Real R2l = sqrt(x1v * x1v + x2l * x2l);
    Real R2r = sqrt(x1v * x1v + x2r * x2r);

    // Don't extrapolate past the last entry in the table
    if (R > disk_profile.GetRmax()) return;

    // Calculate Gravitational Potentialsi
    Real c_out = (4.0/3.0) * pow(5 * r_m, 1.5);

    Real phi0    = GravPot(x1v, x2v, 0.0, G, r_s, rho_s, m_g, a_g, z_g, c_out, rho_m);
    Real phi0_1l = GravPot(x1l, x2v, 0.0, G, r_s, rho_s, m_g, a_g, z_g, c_out, rho_m);
    Real phi0_1r = GravPot(x1r, x2v, 0.0, G, r_s, rho_s, m_g, a_g, z_g, c_out, rho_m);
    Real phi0_2l = GravPot(x1v, x2l, 0.0, G, r_s, rho_s, m_g, a_g, z_g, c_out, rho_m);
    Real phi0_2r = GravPot(x1v, x2r, 0.0, G, r_s, rho_s, m_g, a_g, z_g, c_out, rho_m);

    Real phi   = GravPot(x1v, x2v, x3v, G, r_s, rho_s, m_g, a_g, z_g, c_out, rho_m);
    Real phi1r = GravPot(x1r, x2v, x3v, G, r_s, rho_s, m_g, a_g, z_g, c_out, rho_m);
    Real phi2l = GravPot(x1v, x2l, x3v, G, r_s, rho_s, m_g, a_g, z_g, c_out, rho_m);
    Real phi2r = GravPot(x1v, x2r, x3v, G, r_s, rho_s, m_g, a_g, z_g, c_out, rho_m);
    Real phi1l = GravPot(x1l, x2v, x3v, G, r_s, rho_s, m_g, a_g, z_g, c_out, rho_m);

    Real f_x1_ = -(phi1r - phi1l) / (x1r - x1l);
    Real f_x2_ = -(phi2r - phi2l) / (x2r - x2l);
    Real dPhi_dR = f_x1_ * x1v / R  + f_x2_ * x2v / R;

    // Compute sound speed squared
    Real temp = disk_profile.GetTemperature(R);
    Real cs2 = temp; // isothermal, no factor of gamma

    // Get densities from profiles via interpolation
    Real rho = disk_profile.GetDensity(R) * exp(-(phi - phi0) / cs2);
    Real rho_1l = disk_profile.GetDensity(R1l) * exp(-(phi1l - phi0_1l) / cs2);
    Real rho_1r = disk_profile.GetDensity(R1r) * exp(-(phi1r - phi0_1r) / cs2);
    Real rho_2l = disk_profile.GetDensity(R2l) * exp(-(phi2l - phi0_2l) / cs2);
    Real rho_2r = disk_profile.GetDensity(R2r) * exp(-(phi2r - phi0_2r) / cs2);
    Real p_x1_ = temp*(rho_1r - rho_1l) / (x1r - x1l);
    Real p_x2_ = temp*(rho_2r - rho_2l) / (x2r - x2l);

    Real dP_dR_over_rho = 0.0;
    constexpr Real tiny = 1.0e-20;
    if (rho > tiny) {
      dP_dR_over_rho = (p_x1_ * x1v  + p_x2_ * x2v) / (R * rho);
    }

    // Calculate circular velocity
    Real v_phi = sqrt(R*fmax(dP_dR_over_rho - dPhi_dR, 0.0));

    // Calculate azimuthal velocity
    Real v1 = 0.0, v2 = 0.0;
    if (R > tiny) {  // Avoid division by zero
      // Calculate azimuthal velocity components
      v1 = -v_phi * x2v / R;
      v2 =  v_phi * x1v / R;
    }

    // Combine cgm and disk material consistently
    Real rho_cgm  = u0(m, IDN, k, j, i);
    Real mom1_cgm = u0(m, IM1, k, j, i);
    Real mom2_cgm = u0(m, IM2, k, j, i);
    Real mom3_cgm = u0(m, IM3, k, j, i);
    Real E_cgm    = u0(m, IEN, k, j, i);

    Real KE_cgm  = 0.5*(SQR(mom1_cgm) + SQR(mom2_cgm) + SQR(mom3_cgm))/rho_cgm;
    Real Eth_cgm = E_cgm - KE_cgm;

    Real rho_tot  = rho_cgm + rho;
    Real mom1_tot = mom1_cgm + rho * v1;
    Real mom2_tot = mom2_cgm + rho * v2;
    Real mom3_tot = mom3_cgm;

    Real KE_tot = 0.5*(SQR(mom1_tot)+SQR(mom2_tot)+SQR(mom3_tot))/rho_tot;

    // Set state variables
    u0(m, IDN, k, j, i) = rho_tot;
    u0(m, IM1, k, j, i) = mom1_tot;
    u0(m, IM2, k, j, i) = mom2_tot;
    u0(m, IM3, k, j, i) = mom3_tot;
    u0(m, IEN, k, j, i) = Eth_cgm + (rho * temp)/gm1 + KE_tot;

    // Update dust only with disk material
    u0(m, IZS, k, j, i) += (1.0 - dz) * Z * Zsol * rho;
    u0(m, IDS, k, j, i) += 0.5 * dz * Z * Zsol * rho;
    u0(m, IDL, k, j, i) += 0.5 * dz * Z * Zsol * rho;
}

//===========================================================================//
//                              Source Terms                                 //
//===========================================================================//

void UserSource(Mesh* pm, const Real bdt) {
  GravitySource(pm, bdt);
  DustSource(pm, bdt);
  SNSource(pm, bdt);
  return;
}

void GravitySource(Mesh* pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb1 = pmbp->nmb_thispack - 1;
  auto &size = pmbp->pmb->mb_size;
  auto &u0 = (pmbp->pmhd != nullptr) ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  auto &w0 = (pmbp->pmhd != nullptr) ? pmbp->pmhd->w0 : pmbp->phydro->w0;

  Real G = pmbp->punit->grav_constant();
  Real r_s = r_scale;
  Real rho_s = rho_scale;
  Real m_g = m_gal;
  Real a_g = a_gal;
  Real z_g = z_gal;
  Real rho_m = rho_mean;
  Real c_out = (4.0/3.0) * pow(5 * r_200, 1.5);

  par_for("gravity_source", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real x1min = size.d_view(m).x1min, x1max = size.d_view(m).x1max;
    const Real x2min = size.d_view(m).x2min, x2max = size.d_view(m).x2max;
    const Real x3min = size.d_view(m).x3min, x3max = size.d_view(m).x3max;

    const Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
    const Real x1l = LeftEdgeX(i-is,   nx1, x1min, x1max);
    const Real x1r = LeftEdgeX(i+1-is, nx1, x1min, x1max);

    const Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
    const Real x2l = LeftEdgeX(j-js,   nx2, x2min, x2max);
    const Real x2r = LeftEdgeX(j+1-js, nx2, x2min, x2max);

    const Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
    const Real x3l = LeftEdgeX(k-ks,   nx3, x3min, x3max);
    const Real x3r = LeftEdgeX(k+1-ks, nx3, x3min, x3max);

    Real phi1l = GravPot(x1l,x2v,x3v,G,r_s,rho_s,m_g,a_g,z_g,c_out,rho_m);
    Real phi1r = GravPot(x1r,x2v,x3v,G,r_s,rho_s,m_g,a_g,z_g,c_out,rho_m);

    Real phi2l = GravPot(x1v,x2l,x3v,G,r_s,rho_s,m_g,a_g,z_g,c_out,rho_m);
    Real phi2r = GravPot(x1v,x2r,x3v,G,r_s,rho_s,m_g,a_g,z_g,c_out,rho_m);

    Real phi3l = GravPot(x1v,x2v,x3l,G,r_s,rho_s,m_g,a_g,z_g,c_out,rho_m);
    Real phi3r = GravPot(x1v,x2v,x3r,G,r_s,rho_s,m_g,a_g,z_g,c_out,rho_m);

    constexpr Real tiny = 1e-20;
    Real f_x1_ = -(phi1r - phi1l) / fmax(x1r - x1l, tiny);
    Real f_x2_ = -(phi2r - phi2l) / fmax(x2r - x2l, tiny);
    Real f_x3_ = -(phi3r - phi3l) / fmax(x3r - x3l, tiny);

    Real density = w0(m, IDN, k, j, i);
    Real src_x1 = bdt * density * f_x1_;
    Real src_x2 = bdt * density * f_x2_;
    Real src_x3 = bdt * density * f_x3_;

    u0(m,IM1,k,j,i) += src_x1;
    u0(m,IM2,k,j,i) += src_x2;
    u0(m,IM3,k,j,i) += src_x3;
    u0(m,IEN,k,j,i) += (src_x1 * w0(m,IVX,k,j,i) +
                        src_x2 * w0(m,IVY,k,j,i) +
                        src_x3 * w0(m,IVZ,k,j,i));

  });

  return;
}

KOKKOS_INLINE_FUNCTION
Real GravPot(Real x1, Real x2, Real x3,
             Real G, Real r_s, Real rho_s,
             Real M_gal, Real a_gal, Real z_gal,
             Real c_outer, Real rho_mean) {
  const Real R2 = fma(x1, x1 , x2*x2);
  const Real R  = sqrt(R2);
  const Real r2 = fma(x3 , x3 , R2);
  const Real r  = sqrt(fmax(r2, 1e-20));

  // NFW component
  Real x = r / r_s;
  Real phi_NFW = -4 * M_PI * G * rho_s * SQR(r_s) * log1p(x) / x;

  // Miyamoto-Nagai model
  Real phi_MN = -G * M_gal / sqrt(R2 + SQR(sqrt(fma(x3 , x3 , z_gal*z_gal)) + a_gal));

  // Outer component
  Real phi_Outer = 4 * M_PI * G * rho_mean * (c_outer * sqrt(r) + (1.0/6.0) * r2);

  // Total potential
  return phi_NFW + phi_MN + phi_Outer;
}

void SNSource(Mesh* pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb1 = pmbp->nmb_thispack - 1;
  auto &size = pmbp->pmb->mb_size;
  int IZS = scalar_IZS;
  int IDS = scalar_IDS;
  int IDL = scalar_IDL;

  auto &u0 = (pmbp->pmhd != nullptr) ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  auto &w0 = (pmbp->pmhd != nullptr) ? pmbp->pmhd->w0 : pmbp->phydro->w0;

  Real dr = r_inj;
  Real dz_sn = d_z_sn;

  if (pm->ncycle != last_sn_detect_cycle) {
    last_sn_detect_cycle = pm->ncycle;

    auto &pr = pmbp->ppart->prtcl_rdata;
    auto &pi = pmbp->ppart->prtcl_idata;
    int npart = pmbp->ppart->nprtcl_thispack;

    auto gids = pmbp->gids;

    Real time = pm->time;
    int nrdata = pmbp->ppart->nrdata;
    Real unit_time = pmbp->punit->time_cgs();

    // Array of positions where SNs go off at this timestep
    auto &sn_centers = sn_centers_buffer;
    Kokkos::deep_copy(sn_counter, 0);
    auto d_counter = sn_counter;

    par_for("sn_source", DevExeSpace(), 0, npart-1, KOKKOS_LAMBDA(const int p) {

      Real next_sn_time = pr(nrdata-1, p);

      if (time > next_sn_time) {
        // Update particle sn tracking
        pi(2, p) += 1;
        int sn_idx = pi(2, p);
        Real par_t_create = pr(nrdata-3, p);
        Real cluster_mass = pr(nrdata-2, p);
        pr(nrdata-1, p) = GetNthSNTime(cluster_mass, par_t_create, unit_time, sn_idx);

        // Register SN center location and adjust for boundaries
        // Warning : if r_inj is large this will lead to unexpected behavior at the start
        int idx = Kokkos::atomic_fetch_add(&d_counter(), 1);
        int m = pi(PGID, p) - gids;

        Real x1min = size.d_view(m).x1min;
        Real x1max = size.d_view(m).x1max;
        Real x2min = size.d_view(m).x2min;
        Real x2max = size.d_view(m).x2max;
        Real x3min = size.d_view(m).x3min;
        Real x3max = size.d_view(m).x3max;

        sn_centers(0, idx) = fmin(fmax(pr(IPX,p), x1min+dr), x1max-dr);
        sn_centers(1, idx) = fmin(fmax(pr(IPY,p), x2min+dr), x2max-dr);
        sn_centers(2, idx) = fmin(fmax(pr(IPZ,p), x3min+dr), x3max-dr);
      }
    });

    DevExeSpace().fence();

    // Get number of SNe on host
    Kokkos::deep_copy(num_sn_this_cycle, d_counter);
  }

  if (num_sn_this_cycle > 0) {
    Real beta = bdt / pm->dt;
    Real e_sn_ = e_sn * beta;
    Real m_ej_ = m_ej * beta;
    Real Z_ej_ = Z_ej;

    auto &sn_centers = sn_centers_buffer;
    int num_sn = num_sn_this_cycle;

    par_for("sn_injection", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real x1min = size.d_view(m).x1min;
      Real x1max = size.d_view(m).x1max;
      Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

      Real x2min = size.d_view(m).x2min;
      Real x2max = size.d_view(m).x2max;
      Real x2v = CellCenterX(j-js, nx2, x2min, x2max);

      Real x3min = size.d_view(m).x3min;
      Real x3max = size.d_view(m).x3max;
      Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);

      for (int sn = 0; sn < num_sn; ++sn) {
        Real sn_x = sn_centers(0, sn);
        Real sn_y = sn_centers(1, sn);
        Real sn_z = sn_centers(2, sn);

        // Calculate distance from SN center
        Real dx = x1v - sn_x;
        Real dy = x2v - sn_y;
        Real dz = x3v - sn_z;
        Real r = sqrt(dx*dx + dy*dy + dz*dz);

        // Inject if within injection radius
        if (r <= dr) {
	  u0(m,IDN,k,j,i) += m_ej_;
          u0(m,IEN,k,j,i) += e_sn_;
	  u0(m, IZS, k, j, i) += (1.-dz_sn) * Z_ej_ * m_ej_;
          u0(m, IDS, k, j, i) += 0.5 * dz_sn * Z_ej_ * m_ej_;
          u0(m, IDL, k, j, i) += 0.5 * dz_sn * Z_ej_ * m_ej_;
	}
      }

    });
  }

  return;
}

KOKKOS_INLINE_FUNCTION
Real LimitDustTransfer(const Real requested,
                       const Real donor_mass,
                       const Real max_drain_frac) {
  if (!(requested > 0.0) || !(donor_mass > 0.0)) return 0.0;
  return fmin(requested, donor_mass * max_drain_frac);
}

void DustSource(Mesh* pm, const Real bdt) {

  if (pm->time < sn_delay) return;

  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmbp->nmb_thispack - 1;

  EOS_Data &eos = (pmbp->pmhd != nullptr) ? pmbp->pmhd->peos->eos_data :
                                            pmbp->phydro->peos->eos_data;
  Real gamma = eos.gamma;
  Real gm1 = gamma - 1.0;

  int IZS = scalar_IZS;
  int IDS = scalar_IDS;
  int IDL = scalar_IDL;

  auto &u0 = (pmbp->pmhd != nullptr) ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  auto &w0 = (pmbp->pmhd != nullptr) ? pmbp->pmhd->w0 : pmbp->phydro->w0;

  // --- grain parameters (from input deck) --------------------------------
  // Representative bin radii used by all dust processes.
  Real rho_grain_cgs    = rho_gr;    // internal grain density [g cm^-3]
  Real a_small          = a_gr_s;    // PAH-like / ultrasmall carbonaceous [μm]
  Real a_large          = a_gr_l;    // representative large-grain radius [μm]

  // --- physical constants ------------------------------------------------
  const Real X_H              = 0.75;      // hydrogen mass fraction
  const Real Z_solar          = Z_sol;     // solar metallicity
  const Real f_coag           = 0.5;       // coagulation efficiency
  const Real v_coag_ref_kms   = 0.1;       // coag. reference velocity [km/s]

  // --- safety parameters -------------------------------------------------
  const Real D_floor          = min_dust_frac;
  const Real max_drain_frac   = dust_donor_cap_frac;
  const Real acc_nH_min       = accretion_nH_min;
  const Real acc_nH_max       = accretion_nH_max;
  const Real acc_T_max        = accretion_T_max;
  const Real acc_f_ref        = accretion_refractory_frac;

  // --- unit conversions --------------------------------------------------
  //  temp_to_K      : code temperature → Kelvin
  //  rho_to_nH      : code density → hydrogen number density [cm^-3]
  //  vel_to_kms     : code velocity → km/s
  //  dt_Myr         : sub-step duration in Myr
  const Real temp_to_K  = pmbp->punit->temperature_cgs();
  const Real rho_to_nH  = X_H * pmbp->punit->density_cgs()
                              / pmbp->punit->atomic_mass_unit_cgs;
  const Real vel_to_kms = 1.0 / pmbp->punit->km_s();
  const Real dt_Myr     = bdt / pmbp->punit->myr();

  par_for("dust_source", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real rho_gas = w0(m, IDN, k, j, i);                 // code units
    const Real nH      = rho_to_nH * rho_gas;                 // cm^-3
    const Real T_K     = temp_to_K * gm1 * w0(m, IEN, k, j, i) / rho_gas; // K
    const Real cs_kms  = vel_to_kms // Sound speed in km/s (used by shattering timescale)
                       * sqrt(gamma * gm1 * w0(m, IEN, k, j, i) / rho_gas);

    const Real Z_gas   = w0(m, IZS, k, j, i);   // gas-phase metal fraction
    const Real D_small = w0(m, IDS, k, j, i);   // small-grain dust fraction
    const Real D_large = w0(m, IDL, k, j, i);   // large-grain dust fraction
    const Real D_total = D_small + D_large;     // total dust fraction

    // Dust *mass densities* in code units (what the conserved scalars store).
    const Real rho_Ds = D_small * rho_gas;
    const Real rho_Dl = D_large * rho_gas;

    // Accumulators for net conserved-density changes (code units).
    Real delta_gas_metal = 0.0;   // Δ(ρ * Z_gas)
    Real delta_small     = 0.0;   // Δ(ρ * D_small)
    Real delta_large     = 0.0;   // Δ(ρ * D_large)

    // =================================================================
    //  1.  THERMAL SPUTTERING  (dust → gas-phase metals)
    //      Tsai & Mathews (1995), McKinnon et al. (2017)
    //
    //      Note : τ_sp = t_sp / a where [a] = μm
    //      τ_sp = 1657 Myr × (10^-3 / nH) × [1 + (T / 2×10^6)^{-5/2}]
    //
    //      Δρ_D = -dt × (3 / a) × ρ_D / τ_sp
    //
    //      The factor 3/a arises because grain mass ∝ a³ and the
    //      absolute erosion rate da/dt is size-independent, so
    //      d(ln M)/dt = 3 (da/dt)/a = 3 / (a τ_sp).
    // =================================================================
    {
      Real tau_sput_Myr = 1657.0 * (1.0e-3 / nH);
      const Real T_ratio = T_K / 2.0e6;
      tau_sput_Myr *= 1.0 + 1.0 / (T_ratio * T_ratio * sqrt(T_ratio));

      Real sput_small = 0.0, sput_large = 0.0;
      if (D_small > D_floor) {
        const Real req = dt_Myr * 3.0 * rho_Ds / (tau_sput_Myr * a_small);
        sput_small = LimitDustTransfer(req, rho_Ds, max_drain_frac);
      }
      if (D_large > D_floor) {
        const Real req = dt_Myr * 3.0 * rho_Dl / (tau_sput_Myr * a_large);
        sput_large = LimitDustTransfer(req, rho_Dl, max_drain_frac);
      }

      delta_gas_metal += sput_small + sput_large;
      delta_small     -= sput_small;
      delta_large     -= sput_large;
    }

    // =================================================================
    //  2.  GAS-PHASE ACCRETION  (refractory gas-phase metals → dust)
    //      Conservative revised t_grow model with Le Bourlot sticking
    //
    //      Z_acc = f_ref × Z_gas
    //      α_LB(T) = [1 + 10^-4 T^(3/2)]^-1
    //      τ_acc = 150 Myr × (100 / n_H) × √(50 K / T)
    //                      × (Z_sol / Z_acc) × (1 + D_tot / Z_acc)
    //                      / α_LB(T)
    //
    //      Δρ_D = +dt × ρ_D / (τ_acc × a)
    // =================================================================
    {
      const Real Z_acc = acc_f_ref * Z_gas;
      Real tau_acc_Myr = 150.0;
      Real accr_small = 0.0, accr_large = 0.0;
      const bool accretion_active = (T_K > 0.0) &&
                                    (T_K < acc_T_max) &&
                                    (nH >= acc_nH_min) &&
                                    (nH <= acc_nH_max) &&
                                    (Z_acc > D_floor);
      if (accretion_active) {
        const Real alpha_LB = 1.0 / (1.0 + 1.0e-4 * T_K * sqrt(T_K));
        tau_acc_Myr *= (100.0 / nH);
        tau_acc_Myr *= sqrt(50.0 / T_K);
        tau_acc_Myr *= (Z_solar / Z_acc);
        tau_acc_Myr *= (1.0 + D_total / Z_acc);
        tau_acc_Myr /= alpha_LB;

        if (D_small > D_floor)
          accr_small = dt_Myr * rho_Ds / (tau_acc_Myr * a_small);
        if (D_large > D_floor)
          accr_large = dt_Myr * rho_Dl / (tau_acc_Myr * a_large);
      }

      // Cap total accretion so it cannot drain more than max_drain_frac
      // of the available refractory gas-phase metal reservoir.
      const Real accr_total = accr_small + accr_large;
      if (accr_total > 0.0) {
        const Real max_accretion = Z_acc * rho_gas * max_drain_frac;
        if (accr_total > max_accretion) {
          const Real scale = max_accretion / accr_total;
          accr_small *= scale;
          accr_large *= scale;
        }
      }

      delta_gas_metal -= (accr_small + accr_large);   // removed from gas
      delta_small     += accr_small;
      delta_large     += accr_large;
    }

    // =================================================================
    //  3.  SHATTERING  (large grains → small grains)
    //      Dubois et al. (2024)
    //
    //      τ_shatt = 540 Myr × (1 / n_H) × (ρ_gr / 3 g cm^-3)
    //                        × (0.01 / D_large) × (10 km/s / c_s)
    //
    //      Δρ_Dl = -dt × ρ_Dl / (τ_shatt × a_large)
    // =================================================================
    {
      Real tau_shatt_Myr = 540.0;
      tau_shatt_Myr *= (1.0 / nH);
      tau_shatt_Myr *= (rho_grain_cgs / 3.0);
      tau_shatt_Myr *= (0.01 / D_large);
      tau_shatt_Myr *= (10.0 / cs_kms);

      Real shatt = 0.0;
      if (D_large > D_floor) {
        const Real req = dt_Myr * rho_Dl / (tau_shatt_Myr * a_large);
        shatt = LimitDustTransfer(req, rho_Dl, max_drain_frac);
      }

      delta_large -= shatt;
      delta_small += shatt;
    }

    // =================================================================
    //  4.  COAGULATION  (small grains → large grains)
    //      Dubois et al. (2024)
    //
    //      τ_coag = 2.7 Myr × (ρ_gr / 3 g cm^-3) × (10^3 / n_H)
    //                        × (0.01 / D_small) × (0.1 km/s / v_coag)
    //                        × f_coag
    //
    //      Δρ_Ds = -dt × ρ_Ds / (τ_coag × a_small / 0.05)
    // =================================================================
    {
      Real tau_coag_Myr = 2.7;
      tau_coag_Myr *= (rho_grain_cgs / 3.0);
      tau_coag_Myr *= (1.0e3 / nH);
      tau_coag_Myr *= (0.01 / D_small);
      tau_coag_Myr *= (0.1 / v_coag_ref_kms);
      tau_coag_Myr /= f_coag;

      Real coag = 0.0;
      if (D_small > D_floor) {
        const Real req = dt_Myr * rho_Ds / (tau_coag_Myr * (a_small / 0.005));
        coag = LimitDustTransfer(req, rho_Ds, max_drain_frac);
      }

      delta_small -= coag;
      delta_large += coag;
    }

    // =====================================================================
    //  Apply accumulated changes to conserved scalar densities
    // =====================================================================
    u0(m, IZS, k, j, i) += delta_gas_metal;
    u0(m, IDS, k, j, i) += delta_small;
    u0(m, IDL, k, j, i) += delta_large;

    // =====================================================================
    //  Safety clamps
    // =====================================================================

    // Gas-phase metallicity: clamp to [0, ρ_gas].
    u0(m, IZS, k, j, i) = fmin(fmax(u0(m, IZS, k, j, i), 0.0), rho_gas);

    // Total dust: clamp so D_total ≤ ρ_gas  (dust fraction ≤ 1).
    const Real rho_dust_total = u0(m, IDS, k, j, i) + u0(m, IDL, k, j, i);
    const Real rescale = (rho_dust_total > 0.0)
                       ? fmin(1.0, rho_gas / rho_dust_total)
                       : 0.0;
    u0(m, IDS, k, j, i) *= rescale;
    u0(m, IDL, k, j, i) *= rescale;

    // Floor: ensure dust densities never drop below D_floor × ρ_gas.
    u0(m, IDS, k, j, i) = fmax(D_floor * rho_gas, u0(m, IDS, k, j, i));
    u0(m, IDL, k, j, i) = fmax(D_floor * rho_gas, u0(m, IDL, k, j, i));
  });

  return;
}

void ApplyUserMagneticBoundary(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp->pmhd == nullptr) return;

  auto &indcs = pm->mb_indcs;
  int &ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int nmb1 = pmbp->nmb_thispack - 1;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  auto &b0 = pmbp->pmhd->b0;

  if (pm->mesh_bcs[BoundaryFace::inner_x1] != BoundaryFlag::periodic) {
    par_for("gotham_b_user_x1", DevExeSpace(), 0, nmb1, 0, (n3-1), 0, (n2-1),
    KOKKOS_LAMBDA(int m, int k, int j) {
      if (mb_bcs.d_view(m, BoundaryFace::inner_x1) == BoundaryFlag::user) {
        for (int i = 0; i < ng; ++i) {
          b0.x1f(m,k,j,is-i-1) = b0.x1f(m,k,j,is);
          b0.x2f(m,k,j,is-i-1) = b0.x2f(m,k,j,is);
          if (j == n2-1) b0.x2f(m,k,j+1,is-i-1) = b0.x2f(m,k,j+1,is);
          b0.x3f(m,k,j,is-i-1) = b0.x3f(m,k,j,is);
          if (k == n3-1) b0.x3f(m,k+1,j,is-i-1) = b0.x3f(m,k+1,j,is);
        }
      }
      if (mb_bcs.d_view(m, BoundaryFace::outer_x1) == BoundaryFlag::user) {
        for (int i = 0; i < ng; ++i) {
          b0.x1f(m,k,j,ie+i+2) = b0.x1f(m,k,j,ie+1);
          b0.x2f(m,k,j,ie+i+1) = b0.x2f(m,k,j,ie);
          if (j == n2-1) b0.x2f(m,k,j+1,ie+i+1) = b0.x2f(m,k,j+1,ie);
          b0.x3f(m,k,j,ie+i+1) = b0.x3f(m,k,j,ie);
          if (k == n3-1) b0.x3f(m,k+1,j,ie+i+1) = b0.x3f(m,k+1,j,ie);
        }
      }
    });
  }
  if (pm->one_d) return;

  if (pm->mesh_bcs[BoundaryFace::inner_x2] != BoundaryFlag::periodic) {
    par_for("gotham_b_user_x2", DevExeSpace(), 0, nmb1, 0, (n3-1), 0, (n1-1),
    KOKKOS_LAMBDA(int m, int k, int i) {
      if (mb_bcs.d_view(m, BoundaryFace::inner_x2) == BoundaryFlag::user) {
        for (int j = 0; j < ng; ++j) {
          b0.x1f(m,k,js-j-1,i) = b0.x1f(m,k,js,i);
          if (i == n1-1) b0.x1f(m,k,js-j-1,i+1) = b0.x1f(m,k,js,i+1);
          b0.x2f(m,k,js-j-1,i) = b0.x2f(m,k,js,i);
          b0.x3f(m,k,js-j-1,i) = b0.x3f(m,k,js,i);
          if (k == n3-1) b0.x3f(m,k+1,js-j-1,i) = b0.x3f(m,k+1,js,i);
        }
      }
      if (mb_bcs.d_view(m, BoundaryFace::outer_x2) == BoundaryFlag::user) {
        for (int j = 0; j < ng; ++j) {
          b0.x1f(m,k,je+j+1,i) = b0.x1f(m,k,je,i);
          if (i == n1-1) b0.x1f(m,k,je+j+1,i+1) = b0.x1f(m,k,je,i+1);
          b0.x2f(m,k,je+j+2,i) = b0.x2f(m,k,je+1,i);
          b0.x3f(m,k,je+j+1,i) = b0.x3f(m,k,je,i);
          if (k == n3-1) b0.x3f(m,k+1,je+j+1,i) = b0.x3f(m,k+1,je,i);
        }
      }
    });
  }
  if (pm->two_d) return;

  if (pm->mesh_bcs[BoundaryFace::inner_x3] != BoundaryFlag::periodic) {
    par_for("gotham_b_user_x3", DevExeSpace(), 0, nmb1, 0, (n2-1), 0, (n1-1),
    KOKKOS_LAMBDA(int m, int j, int i) {
      if (mb_bcs.d_view(m, BoundaryFace::inner_x3) == BoundaryFlag::user) {
        for (int k = 0; k < ng; ++k) {
          b0.x1f(m,ks-k-1,j,i) = b0.x1f(m,ks,j,i);
          if (i == n1-1) b0.x1f(m,ks-k-1,j,i+1) = b0.x1f(m,ks,j,i+1);
          b0.x2f(m,ks-k-1,j,i) = b0.x2f(m,ks,j,i);
          if (j == n2-1) b0.x2f(m,ks-k-1,j+1,i) = b0.x2f(m,ks,j+1,i);
          b0.x3f(m,ks-k-1,j,i) = b0.x3f(m,ks,j,i);
        }
      }
      if (mb_bcs.d_view(m, BoundaryFace::outer_x3) == BoundaryFlag::user) {
        for (int k = 0; k < ng; ++k) {
          b0.x1f(m,ke+k+1,j,i) = b0.x1f(m,ke,j,i);
          if (i == n1-1) b0.x1f(m,ke+k+1,j,i+1) = b0.x1f(m,ke,j,i+1);
          b0.x2f(m,ke+k+1,j,i) = b0.x2f(m,ke,j,i);
          if (j == n2-1) b0.x2f(m,ke+k+1,j+1,i) = b0.x2f(m,ke,j+1,i);
          b0.x3f(m,ke+k+2,j,i) = b0.x3f(m,ke+1,j,i);
        }
      }
    });
  }
}

void AddUserBoundaryMagneticEnergy(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp->pmhd == nullptr) return;

  auto &indcs = pm->mb_indcs;
  int &ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int nmb1 = pmbp->nmb_thispack - 1;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  auto &u0 = pmbp->pmhd->u0;
  auto &b0 = pmbp->pmhd->b0;

  par_for("gotham_user_b_energy", DevExeSpace(), 0, nmb1, 0, (n3-1), 0, (n2-1),
          0, (n1-1), KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const bool user_x1 = ((i < is) &&
        (mb_bcs.d_view(m, BoundaryFace::inner_x1) == BoundaryFlag::user)) ||
        ((i > ie) && (mb_bcs.d_view(m, BoundaryFace::outer_x1) == BoundaryFlag::user));
    const bool user_x2 = ((j < js) &&
        (mb_bcs.d_view(m, BoundaryFace::inner_x2) == BoundaryFlag::user)) ||
        ((j > je) && (mb_bcs.d_view(m, BoundaryFace::outer_x2) == BoundaryFlag::user));
    const bool user_x3 = ((k < ks) &&
        (mb_bcs.d_view(m, BoundaryFace::inner_x3) == BoundaryFlag::user)) ||
        ((k > ke) && (mb_bcs.d_view(m, BoundaryFace::outer_x3) == BoundaryFlag::user));

    if (user_x1 || user_x2 || user_x3) {
      const Real bx = 0.5*(b0.x1f(m,k,j,i) + b0.x1f(m,k,j,i+1));
      const Real by = 0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i));
      const Real bz = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
      u0(m,IEN,k,j,i) += 0.5*(bx*bx + by*by + bz*bz);
    }
  });
}

//===========================================================================//
//                             User Boundary                                 //
//===========================================================================//

void UserBoundary(Mesh* pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  int &ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng) : 1;
  int &is = indcs.is;  int &ie  = indcs.ie;
  int &js = indcs.js;  int &je  = indcs.je;
  int &ks = indcs.ks;  int &ke  = indcs.ke;

  int nmb1 = pmbp->nmb_thispack - 1;
  auto &size = pmbp->pmb->mb_size;
  auto &mb_bcs = pm->pmb_pack->pmb->mb_bcs;

  int IZS = scalar_IZS;
  int IDS = scalar_IDS;
  int IDL = scalar_IDL;

  auto &u0 = (pmbp->pmhd != nullptr) ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  auto &eos = (pmbp->pmhd != nullptr) ? pmbp->pmhd->peos->eos_data :
                                        pmbp->phydro->peos->eos_data;
  Real gm1 = eos.gamma - 1.0;

  auto &profile = profile_reader;

  Real r_c = r_circ;
  Real v_c = v_circ;

  Real Zsol = Z_sol;
  Real Z_ = Z;
  Real dz_init = d_z_init;
  Real min_df = min_dust_frac;

  if (pmbp->pmhd != nullptr) {
    ApplyUserMagneticBoundary(pm);
  }

  // Handle X1 boundaries
  par_for("static_x1", DevExeSpace(), 0, nmb1, 0, (n3-1), 0, (n2-1), 0, (ng-1),
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;

    Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);
    Real x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);

    // Inner X1 boundary
    if (mb_bcs.d_view(m, BoundaryFace::inner_x1) == BoundaryFlag::user) {
      Real x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);

      SetCoolingFlowState(u0, m, k, j, i, x1v, x2v, x3v, gm1, profile);
      SetRotation(u0, m, k, j, i, x1v, x2v, x3v, r_c, v_c);
      SetDustScalars(u0, m, k, j, i, IZS, IDS, IDL, min_df, Z_, Zsol);
    }

    // Outer X1 boundary
    if (mb_bcs.d_view(m, BoundaryFace::outer_x1) == BoundaryFlag::user) {
      int i_out = ie + i + 1;
      Real x1v = CellCenterX(i_out-is, indcs.nx1, x1min, x1max);

      SetCoolingFlowState(u0, m, k, j, i_out, x1v, x2v, x3v, gm1, profile);
      SetRotation(u0, m, k, j, i_out, x1v, x2v, x3v, r_c, v_c);
      SetDustScalars(u0, m, k, j, i_out, IZS, IDS, IDL, min_df, Z_, Zsol);
    }
  });

  // Handle X2 boundaries
  par_for("static_x2", DevExeSpace(), 0, nmb1, 0, (n3-1), 0, (ng-1), 0, (n1-1),
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;

    Real x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);
    Real x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);

    // Inner X2 boundary
    if (mb_bcs.d_view(m, BoundaryFace::inner_x2) == BoundaryFlag::user) {
      Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);

      SetCoolingFlowState(u0, m, k, j, i, x1v, x2v, x3v, gm1, profile);
      SetRotation(u0, m, k, j, i, x1v, x2v, x3v, r_c, v_c);
      SetDustScalars(u0, m, k, j, i, IZS, IDS, IDL, min_df, Z_, Zsol);
    }

    // Outer X2 boundary
    if (mb_bcs.d_view(m, BoundaryFace::outer_x2) == BoundaryFlag::user) {
      int j_out = je + j + 1;
      Real x2v = CellCenterX(j_out-js, indcs.nx2, x2min, x2max);

      SetCoolingFlowState(u0, m, k, j_out, i, x1v, x2v, x3v, gm1, profile);
      SetRotation(u0, m, k, j_out, i, x1v, x2v, x3v, r_c, v_c);
      SetDustScalars(u0, m, k, j_out, i, IZS, IDS, IDL, min_df, Z_, Zsol);
    }
  });

  // Handle X3 boundaries
  par_for("static_x3", DevExeSpace(), 0, nmb1, 0, (ng-1), 0, (n2-1), 0, (n1-1),
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;

    Real x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);
    Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);

    // Inner X3 boundary
    if (mb_bcs.d_view(m, BoundaryFace::inner_x3) == BoundaryFlag::user) {
      Real x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);

      SetCoolingFlowState(u0, m, k, j, i, x1v, x2v, x3v, gm1, profile);
      SetRotation(u0, m, k, j, i, x1v, x2v, x3v, r_c, v_c);
      SetDustScalars(u0, m, k, j, i, IZS, IDS, IDL, min_df, Z_, Zsol);
    }

    // Outer X3 boundary
    if (mb_bcs.d_view(m, BoundaryFace::outer_x3) == BoundaryFlag::user) {
      int k_out = ke + k + 1;
      Real x3v = CellCenterX(k_out-ks, indcs.nx3, x3min, x3max);

      SetCoolingFlowState(u0, m, k_out, j, i, x1v, x2v, x3v, gm1, profile);
      SetRotation(u0, m, k_out, j, i, x1v, x2v, x3v, r_c, v_c);
      SetDustScalars(u0, m, k_out, j, i, IZS, IDS, IDL, min_df, Z_, Zsol);
    }
  });

  if (pmbp->pmhd != nullptr) {
    AddUserBoundaryMagneticEnergy(pm);
  }

}

//===========================================================================//
//                              Refinement                                   //
//===========================================================================//

// Refine region based on density gradient threshold
void RefinementCondition(MeshBlockPack* pmbp) {
  Mesh *pmesh       = pmbp->pmesh;
  int nmb           = pmbp->nmb_thispack;
  int mbs           = pmesh->gids_eachrank[global_variable::my_rank];
  auto &refine_flag = pmesh->pmr->refine_flag;
  auto &indcs       = pmesh->mb_indcs;
  int &is = indcs.is, nx1 = indcs.nx1;
  int &js = indcs.js, nx2 = indcs.nx2;
  int &ks = indcs.ks, nx3 = indcs.nx3;
  const int nkji = nx3 * nx2 * nx1;
  const int nji  = nx2 * nx1;
  auto &u0       = (pmbp->pmhd != nullptr) ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  auto &w0       = (pmbp->pmhd != nullptr) ? pmbp->pmhd->w0 : pmbp->phydro->w0;

  auto &ddens_thresh = ddens_threshold;

  par_for_outer("UserRefineCond",DevExeSpace(), 0, 0, 0, (nmb-1),
  KOKKOS_LAMBDA(TeamMember_t tmember, const int m) {

    Real team_ddmax;
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tmember, nkji),
    [=](const int idx, Real& ddmax) {
      int k = (idx)/nji;
      int j = (idx - k*nji)/nx1;
      int i = (idx - k*nji - j*nx1) + is;
      j += js;
      k += ks;

      // Calculate density gradient
      Real d2 = (SQR(u0(m,IDN,k,j,i+1) - u0(m,IDN,k,j,i-1))
               + SQR(u0(m,IDN,k,j+1,i) - u0(m,IDN,k,j-1,i))
               + SQR(u0(m,IDN,k+1,j,i) - u0(m,IDN,k-1,j,i)));
      ddmax = fmax((sqrt(d2)/u0(m,IDN,k,j,i)), ddmax);

      // Calculate pressure gradient
      Real p2 = (SQR(w0(m,IEN,k,j,i+1) - w0(m,IEN,k,j,i-1))
               + SQR(w0(m,IEN,k,j+1,i) - w0(m,IEN,k,j-1,i))
	       + SQR(w0(m,IEN,k+1,j,i) - w0(m,IEN,k-1,j,i)));
      ddmax = fmax((sqrt(p2)/w0(m,IEN,k,j,i)), ddmax);
    },Kokkos::Max<Real>(team_ddmax));

    if (team_ddmax > ddens_thresh) {refine_flag.d_view(m+mbs) = 1;}
    if (team_ddmax < 0.25*ddens_thresh) {refine_flag.d_view(m+mbs) = -1;}

  });

  // sync host and device
  refine_flag.template modify<DevExeSpace>();
  refine_flag.template sync<HostMemSpace>();
}

//===========================================================================//
//                            Post Main Loop                                 //
//===========================================================================//

void FreeProfile(ParameterInput *pin, Mesh *pm) {
  // Free Kokkos views before Kokkos::finalize is called
  profile_reader = ProfileReader();
  disk_profile_reader = ProfileReader();
  sn_centers_buffer = DvceArray2D<Real>();
  sn_counter = Kokkos::View<int>();
}
