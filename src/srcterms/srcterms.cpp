//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file srcterms.cpp
//  Implements various (physics) source terms to be added to the Hydro or MHD eqns.
//  Source terms objects are stored in the respective fluid class, so that
//  Hydro/MHD can have different source terms

#include "srcterms.hpp"

#include <iostream>
#include <string> // string

#include "athena.hpp"
#include "coordinates/cartesian_ks.hpp"
#include "coordinates/cell_locations.hpp"
#include "eos/eos.hpp"
#include "geodesic-grid/geodesic_grid.hpp"
#include "hydro/hydro.hpp"
#include "ismcooling.hpp"
#include "mesh/mesh.hpp"
#include "mhd/mhd.hpp"
#include "parameter_input.hpp"
#include "radiation/radiation.hpp"
#include "turb_driver.hpp"
#include "units/units.hpp"
#include "cooling_tables.hpp"

//----------------------------------------------------------------------------------------
// constructor, parses input file and initializes data structures and parameters
// Only source terms specified in input file are initialized.

SourceTerms::SourceTerms(std::string block, MeshBlockPack *pp, ParameterInput *pin) :
  pmy_pack(pp),
  Tbins("Tbins",1), nHbins("nHbins",1), He_mf_bins("He_mf_bins",1),
  Metal_Cooling("Metal_Cooling",1,1), H_He_Cooling("H_He_Cooling",1,1,1),
  shearing_box_r_phi(false) {
  // (1) (constant) gravitational acceleration
  const_accel = pin->GetOrAddBoolean(block, "const_accel", false);
  if (const_accel) {
    const_accel_val = pin->GetReal(block, "const_accel_val");
    const_accel_dir = pin->GetInteger(block, "const_accel_dir");
    if (const_accel_dir < 1 || const_accel_dir > 3) {
      std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__ << std::endl
                << "const_accle_dir must be 1,2, or 3" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // (2) Optically thin (ISM) cooling
  ism_cooling = pin->GetOrAddBoolean(block, "ism_cooling", false);
  if (ism_cooling) {
    hrate = pin->GetReal(block, "hrate");
  }

  // (2b) CGM cooling
  cgm_cooling = pin->GetOrAddBoolean(block, "cgm_cooling", false);
  if (cgm_cooling) {
    // Load cooling tables from cooling tables.hpp
    Kokkos::realloc(Tbins, Tbins_TOTAL_SIZE);
    Kokkos::realloc(nHbins, nHbins_TOTAL_SIZE);
    Kokkos::realloc(He_mf_bins, He_mf_bins_TOTAL_SIZE);
    Kokkos::realloc(Metal_Cooling, Metal_Cooling_DIM_0, Metal_Cooling_DIM_1);
    Kokkos::realloc(H_He_Cooling, H_He_Cooling_DIM_0, H_He_Cooling_DIM_1, H_He_Cooling_DIM_2);

    for (int i = 0; i < Tbins_TOTAL_SIZE; ++i) Tbins(i) = Tbins_ARR[i];
    for (int i = 0; i < nHbins_TOTAL_SIZE; ++i) nHbins(i) = nHbins_ARR[i];
    for (int i = 0; i < He_mf_bins_TOTAL_SIZE; ++i) He_mf_bins(i) = He_mf_bins_ARR[i];

    for (int i = 0; i < Metal_Cooling_DIM_0; ++i) {
      for (int j = 0; j < Metal_Cooling_DIM_1; ++j) {
        Metal_Cooling(i, j) = Metal_Cooling_ARR[i * H_He_Cooling_DIM_2 + j];
    }}

    for (int i = 0; i < H_He_Cooling_DIM_0; ++i) {
      for (int j = 0; j < H_He_Cooling_DIM_1; ++j) {
        for (int k = 0; k < H_He_Cooling_DIM_2; ++k) {
          H_He_Cooling(i, j, k) = H_He_Cooling_ARR[i * (H_He_Cooling_DIM_1 * H_He_Cooling_DIM_2) 
                                                   + j * (H_He_Cooling_DIM_2) + k];
    }}}
  }

  // (3) beam source (radiation)
  beam = pin->GetOrAddBoolean(block, "beam_source", false);
  if (beam) {
    dii_dt = pin->GetReal(block, "dii_dt");
  }

  // (4) cooling (relativistic)
  rel_cooling = pin->GetOrAddBoolean(block, "rel_cooling", false);
  if (rel_cooling) {
    crate_rel = pin->GetReal(block, "crate_rel");
    cpower_rel = pin->GetOrAddReal(block, "cpower_rel", 1.);
  }

  // (5) shearing box
  if (pin->DoesBlockExist("shearing_box")) {
    shearing_box = true;
    qshear = pin->GetReal("shearing_box","qshear");
    omega0 = pin->GetReal("shearing_box","omega0");
  } else {
    shearing_box = false;
  }
}

//----------------------------------------------------------------------------------------
// destructor

SourceTerms::~SourceTerms() {
}

//----------------------------------------------------------------------------------------
//! \fn
// Add constant acceleration
// NOTE source terms must be computed using primitive (w0) and NOT conserved (u0) vars

void SourceTerms::ConstantAccel(const DvceArray5D<Real> &w0, const EOS_Data &eos_data,
                                const Real bdt, DvceArray5D<Real> &u0) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;

  Real &g = const_accel_val;
  int &dir = const_accel_dir;

  par_for("const_acc", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real src = bdt*g*w0(m,IDN,k,j,i);
    u0(m,dir,k,j,i) += src;
    if (eos_data.is_ideal) { u0(m,IEN,k,j,i) += src*w0(m,dir,k,j,i); }
  });

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void SourceTerms::ISMCooling()
//! \brief Add explict ISM cooling and heating source terms in the energy equations.
// NOTE source terms must be computed using primitive (w0) and NOT conserved (u0) vars

void SourceTerms::ISMCooling(const DvceArray5D<Real> &w0, const EOS_Data &eos_data,
                             const Real bdt, DvceArray5D<Real> &u0) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  Real use_e = eos_data.use_e;
  Real gamma = eos_data.gamma;
  Real gm1 = gamma - 1.0;
  Real heating_rate = hrate;
  Real temp_unit = pmy_pack->punit->temperature_cgs();
  Real n_unit = pmy_pack->punit->density_cgs()/pmy_pack->punit->mu()
                /pmy_pack->punit->atomic_mass_unit_cgs;
  Real cooling_unit = pmy_pack->punit->pressure_cgs()/pmy_pack->punit->time_cgs()
                      /n_unit/n_unit;
  Real heating_unit = pmy_pack->punit->pressure_cgs()/pmy_pack->punit->time_cgs()/n_unit;

  par_for("ism_cooling", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // temperature in cgs unit
    Real temp = 1.0;
    if (use_e) {
      temp = temp_unit*w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1;
    } else {
      temp = temp_unit*w0(m,ITM,k,j,i);
    }

    Real lambda_cooling = ISMCoolFn(temp)/cooling_unit;
    Real gamma_heating = heating_rate/heating_unit;

    u0(m,IEN,k,j,i) -= bdt * w0(m,IDN,k,j,i) *
                        (w0(m,IDN,k,j,i) * lambda_cooling - gamma_heating);
  });

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void SourceTerms::CGMCooling()
//! \brief Add explict CGM cooling source term in the energy equations.
// NOTE source terms must be computed using primitive (w0) and NOT conserved (u0) vars

void SourceTerms::CGMCooling(const DvceArray5D<Real> &w0, const EOS_Data &eos_data,
                             const Real bdt, DvceArray5D<Real> &u0) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  Real use_e = eos_data.use_e;
  Real gamma = eos_data.gamma;
  Real gm1 = gamma - 1.0;

  auto &units = pmy_pack->punit; 
  Real temp_unit = units->temperature_cgs();
  Real n_unit = units->density_cgs()/units->mu()/units->atomic_mass_unit_cgs;
  Real cooling_unit = units->pressure_cgs()/units->time_cgs()/n_unit/n_unit;

  // Read only cooling tables
  DvceArray1D<const Real> Tbins_ = Tbins;
  DvceArray1D<const Real> nHbins_ = nHbins;
  DvceArray1D<const Real> He_mf_bins_ = He_mf_bins;
  DvceArray2D<const Real> Metal_Cooling_ = Metal_Cooling;
  DvceArray3D<const Real> H_He_Cooling_ = H_He_Cooling;

  par_for("cgm_cooling", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real temp = 1.0; // temperature in cgs units
    if (use_e) {
      temp = temp_unit*w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1;
    } else {
      temp = temp_unit*w0(m,ITM,k,j,i);
    }
    Real nH = n_unit*w0(m,IDN,k,j,i); // density in cgs units
    Real Z = 1.0/3; // metallicity [Zsun]

    // WiersmaCooling at redshift z = 0 taken from Wiersma et al (2009)
    Real lambda_cooling = 0.0; // Ensure we are in range of cooling table
    if (temp > Tbins_(0) && temp < Tbins_(Tbins_.extent(0) - 1)
        && nH > nHbins_(0)  && nH < nHbins_(nHbins_.extent(0) -1)) {
      // Compute Helium Mass Fraction
      Real He2Habundance = pow(10, -1.07) * (0.71553 + 0.28447 * Z); // Asplund+09, Groves+04
      Real X = (1.0 - 0.014 * Z) / (1.0 + 4.0 * He2Habundance);
      Real Y = 4.0 * He2Habundance * X;

      int iHe = 0;
      // Check if val is outside the bounds of the He_mf_bins array
      if (Y > He_mf_bins_(He_mf_bins_.extent(0) - 1)) {
        iHe = He_mf_bins_.extent(0) - 1;
      }
      else if (Y > He_mf_bins_(0)) {
        Real dx = He_mf_bins_.extent(1) - He_mf_bins_.extent(0);
        // We assume that the spacing in He_mf_bins is regular
        iHe = static_cast<int>((Y + (dx / 2) - He_mf_bins_(0)) / dx);
      }
      
      // Convert input values to log space
      Real log_temp = log10(temp);
      Real log_density = log10(nH);

      // Locate indices in Tbins and nHbins
      int i = 0, j = 0;
      while (i < Tbins_.extent(0) - 1 && log10(Tbins_(i + 1)) < log_temp) ++i; 
      while (j < nHbins_.extent(0) - 1 && log10(nHbins_(j + 1)) < log_density) ++j;

      // Logarithms of the Tbins and nHbins bounding the point
      Real log_T0 = log10(Tbins_(i));
      Real log_T1 = log10(Tbins_(i + 1));
      Real log_nH0 = log10(nHbins_(j));
      Real log_nH1 = log10(nHbins_(j + 1));

      // Compute weights for bilinear interpolation
      Real t = (log_temp - log_T0) / (log_T1 - log_T0);
      Real u = (log_density - log_nH0) / (log_nH1 - log_nH0);
  
      // Corner values in log space from the H_He cooling grid
      Real log_C00 = log10(H_He_Cooling_(iHe, i, j));
      Real log_C10 = log10(H_He_Cooling_(iHe, i + 1, j));
      Real log_C01 = log10(H_He_Cooling_(iHe, i, j + 1));
      Real log_C11 = log10(H_He_Cooling_(iHe, i + 1, j + 1));
  
      // Corner values in log space from the Metal cooling grid
      Real log_M00 = log10(Metal_Cooling_(i, j));
      Real log_M10 = log10(Metal_Cooling_(i + 1, j));
      Real log_M01 = log10(Metal_Cooling_(i, j + 1));
      Real log_M11 = log10(Metal_Cooling_(i + 1, j + 1));
  
      // Bilinear interpolation in log space
      Real log_interpolated_prim_cooling =
          (1 - t) * (1 - u) * log_C00 +
          t * (1 - u) * log_C10 +
          (1 - t) * u * log_C01 +
          t * u * log_C11;
      Real log_interpolated_metal_cooling =
          (1 - t) * (1 - u) * log_M00 +
          t * (1 - u) * log_M10 +
          (1 - t) * u * log_M01 +
          t * u * log_M11;
  
      // Convert back to linear space and return
      Real prim_cooling = pow(10.0, log_interpolated_prim_cooling);
      Real metal_cooling = pow(10.0, log_interpolated_metal_cooling);

  
      lambda_cooling = prim_cooling + Z * metal_cooling;
    }

    u0(m,IEN,k,j,i) -= bdt * pow(w0(m,IDN,k,j,i),2) * lambda_cooling/cooling_unit;
  });

  return;
}

KOKKOS_INLINE_FUNCTION
int find_he_bin(Real val) {
    Real he_0 = He_mf_bins_ARR[0];
    Real he_1 = He_mf_bins_ARR[1];
    Real he_N = He_mf_bins_ARR[He_mf_bins_TOTAL_SIZE - 1];

    // Check if val is outside the bounds of the He_mf_bins array
    if (val < he_0) return 0;    
    if (val > he_N) return He_mf_bins_TOTAL_SIZE - 1;

    // Calculate the step size (dx) between consecutive elements
    Real dx = he_1 - he_0;

    // Calculate the bin index and return it
    // We assume that the spacing in He_mf_bins is regular
    return static_cast<int>((val + (dx / 2) - he_0) / dx);
}

//----------------------------------------------------------------------------------------
//! \fn void SourceTerms::RelCooling()
//! \brief Add explict relativistic cooling in the energy and momentum equations.
// NOTE source terms must be computed using primitive (w0) and NOT conserved (u0) vars

void SourceTerms::RelCooling(const DvceArray5D<Real> &w0, const EOS_Data &eos_data,
                             const Real bdt, DvceArray5D<Real> &u0) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  Real use_e = eos_data.use_e;
  Real gamma = eos_data.gamma;
  Real gm1 = gamma - 1.0;
  Real cooling_rate = crate_rel;
  Real cooling_power = cpower_rel;

  par_for("cooling", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // temperature in cgs unit
    Real temp = 1.0;
    if (use_e) {
      temp = w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1;
    } else {
      temp = w0(m,ITM,k,j,i);
    }

    auto &ux = w0(m,IVX,k,j,i);
    auto &uy = w0(m,IVY,k,j,i);
    auto &uz = w0(m,IVZ,k,j,i);

    auto ut = 1.0 + ux*ux + uy*uy + uz*uz;
    ut = sqrt(ut);

    u0(m,IEN,k,j,i) -= bdt*w0(m,IDN,k,j,i)*ut*pow((temp*cooling_rate), cooling_power);
    u0(m,IM1,k,j,i) -= bdt*w0(m,IDN,k,j,i)*ux*pow((temp*cooling_rate), cooling_power);
    u0(m,IM2,k,j,i) -= bdt*w0(m,IDN,k,j,i)*uy*pow((temp*cooling_rate), cooling_power);
    u0(m,IM3,k,j,i) -= bdt*w0(m,IDN,k,j,i)*uz*pow((temp*cooling_rate), cooling_power);
  });

  return;
}

//----------------------------------------------------------------------------------------
//! \fn SourceTerms::BeamSource()
// \brief Add beam of radiation

void SourceTerms::BeamSource(DvceArray5D<Real> &i0, const Real bdt) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = (pmy_pack->nmb_thispack-1);
  int nang1 = (pmy_pack->prad->prgeo->nangles-1);

  auto &nh_c_ = pmy_pack->prad->nh_c;
  auto &tt = pmy_pack->prad->tet_c;
  auto &tc = pmy_pack->prad->tetcov_c;

  auto &excise = pmy_pack->pcoord->coord_data.bh_excise;
  auto &rad_mask_ = pmy_pack->pcoord->excision_floor;
  Real &n_0_floor_ = pmy_pack->prad->n_0_floor;

  auto &beam_mask_ = pmy_pack->prad->beam_mask;
  Real &dii_dt_ = dii_dt;
  par_for("beam_source",DevExeSpace(),0,nmb1,0,nang1,ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(int m, int n, int k, int j, int i) {
    if (beam_mask_(m,n,k,j,i)) {
      Real n0 = tt(m,0,0,k,j,i);
      Real n_0 = tc(m,0,0,k,j,i)*nh_c_.d_view(n,0) + tc(m,1,0,k,j,i)*nh_c_.d_view(n,1)
               + tc(m,2,0,k,j,i)*nh_c_.d_view(n,2) + tc(m,3,0,k,j,i)*nh_c_.d_view(n,3);
      i0(m,n,k,j,i) += n0*n_0*dii_dt_*bdt;
      // handle excision
      // NOTE(@pdmullen): exicision criterion are not finalized.  The below zeroes all
      // intensities within rks <= 1.0 and zeroes intensities within angles where n_0
      // is about zero.  This needs future attention.
      if (excise) {
        if (rad_mask_(m,k,j,i) || fabs(n_0) < n_0_floor_) { i0(m,n,k,j,i) = 0.0; }
      }
    }
  });

  return;
}
