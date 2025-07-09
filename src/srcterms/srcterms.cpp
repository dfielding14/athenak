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
  Tbins("Tbins",1), nHbins("nHbins",1),
  Metal_Cooling("Metal_Cooling",1,1), H_He_Cooling("H_He_Cooling",1,1),
  Metal_Cooling_CIE("Metal_Cooling_CIE",1), H_He_Cooling_CIE("H_He_Cooling_CIE",1),
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
    // Ensure that ISM cooling is not enabled 
    if (ism_cooling) {
      std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__ << std::endl
                << "CGM cooling and ISM cooling are incompatible" << std::endl;
      std::exit(EXIT_FAILURE);
    };

    // User input parameters for ISM heating
    hrate = pin->GetOrAddReal(block, "hrate", 0.0); // Heating rate in cgs units
    // Normalization factor in density code units 
    hscale_norm = pin->GetOrAddReal(block, "hscale_norm", 0.0); 
    hscale_height = pin->GetOrAddReal(block, "hscale_height", 0.0); // Scale height in code units
    hscale_radius = pin->GetOrAddReal(block, "hscale_radius", 0.0); // Scale radius in code units
    hscale_alpha = pin->GetOrAddReal(block, "hscale_alpha", 0.0); // Scale coeff in code units

    // Initialize Cooling Tables to the right dimensions from cooling_tables.hpp
    Kokkos::realloc(Tbins, Tbins_TOTAL_SIZE);
    Kokkos::realloc(nHbins, nHbins_TOTAL_SIZE);
    Kokkos::realloc(Metal_Cooling, Metal_Cooling_DIM_0, Metal_Cooling_DIM_1);
    Kokkos::realloc(H_He_Cooling, H_He_Cooling_DIM_0, H_He_Cooling_DIM_1);
    Kokkos::realloc(Metal_Cooling_CIE, Metal_Cooling_CIE_DIM_0);
    Kokkos::realloc(H_He_Cooling_CIE, H_He_Cooling_CIE_DIM_0);
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

  Initialize();
}

//----------------------------------------------------------------------------------------
// destructor

SourceTerms::~SourceTerms() {
}

//----------------------------------------------------------------------------------------
//! \fn
//  \brief Function to initialize

void SourceTerms::Initialize() {
  if (cgm_cooling) {
    // Load Cooling Tables from cooling_tables.hpp into host device
    for (int i = 0; i < Tbins_DIM_0; ++i) Tbins.h_view(i) = Tbins_ARR[i];
    for (int i = 0; i < nHbins_DIM_0; ++i) nHbins.h_view(i) = nHbins_ARR[i];

    for (int i = 0; i < Metal_Cooling_DIM_0; ++i) {
      for (int j = 0; j < Metal_Cooling_DIM_1; ++j) {
        Metal_Cooling.h_view(i, j) = Metal_Cooling_ARR[i * Metal_Cooling_DIM_1 + j];
    }}

    for (int i = 0; i < H_He_Cooling_DIM_0; ++i) {
      for (int j = 0; j < H_He_Cooling_DIM_1; ++j) {
        H_He_Cooling.h_view(i, j) = H_He_Cooling_ARR[i * H_He_Cooling_DIM_1 + j];
    }}

    for (int i = 0; i < Metal_Cooling_CIE_DIM_0; ++i) {
      Metal_Cooling_CIE.h_view(i) = Metal_Cooling_CIE_ARR[i];
    }

    for (int i = 0; i < H_He_Cooling_CIE_DIM_0; ++i) {
      H_He_Cooling_CIE.h_view(i) = H_He_Cooling_CIE_ARR[i];
    }

    // Synchronize Cooling Tables to device memory
    Tbins.template modify<HostMemSpace>();
    Tbins.template sync<DevExeSpace>();
    nHbins.template modify<HostMemSpace>();
    nHbins.template sync<DevExeSpace>();
    Metal_Cooling.template modify<HostMemSpace>();
    Metal_Cooling.template sync<DevExeSpace>();
    H_He_Cooling.template modify<HostMemSpace>();
    H_He_Cooling.template sync<DevExeSpace>();
    Metal_Cooling_CIE.template modify<HostMemSpace>();
    Metal_Cooling_CIE.template sync<DevExeSpace>();
    H_He_Cooling_CIE.template modify<HostMemSpace>();
    H_He_Cooling_CIE.template sync<DevExeSpace>();
  }
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
  auto &size = pmy_pack->pmb->mb_size;
  int nscalars = pmy_pack->phydro->nscalars;
  int nhydro = pmy_pack->phydro->nhydro;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  Real use_e = eos_data.use_e;
  Real gamma = eos_data.gamma;
  Real gm1 = gamma - 1.0;

  auto &units = pmy_pack->punit; 
  Real temp_unit = units->temperature_cgs();
  Real nH_unit = units->density_cgs()/units->atomic_mass_unit_cgs;
  Real cooling_unit = units->pressure_cgs()/units->time_cgs()/nH_unit/nH_unit;
  Real heating_unit = units->pressure_cgs()/units->time_cgs()/nH_unit;
  Real length_unit = units->length_cgs();

  auto Tbins_ = Tbins.d_view;
  auto nHbins_ = nHbins.d_view;
  auto Metal_Cooling_ = Metal_Cooling.d_view;
  auto H_He_Cooling_ = H_He_Cooling.d_view;
  auto Metal_Cooling_CIE_ = Metal_Cooling_CIE.d_view;
  auto H_He_Cooling_CIE_ = H_He_Cooling_CIE.d_view;

  auto Tfloor  = Tbins_ARR[0];
  auto Tceil   = Tbins_ARR[Tbins_DIM_0 - 1];
  auto nHfloor = nHbins_ARR[0];
  auto nHceil  = nHbins_ARR[nHbins_DIM_0 - 1];

  Real X = 0.75; // Hydrogen mass fraction
  Real Zsol = 0.02; // Solar metallicity

  Real h_rate = hrate;
  Real h_norm = hscale_norm;
  Real h_height = hscale_height;
  Real h_radius = hscale_radius;
  Real h_alpha = hscale_alpha;

  par_for("cgm_cooling", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    int nx1 = indcs.nx1;
    Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    int nx2 = indcs.nx2;
    Real x2v = CellCenterX(j-js, nx2, x2min, x2max);

    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    int nx3 = indcs.nx3;
    Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);

    // Get metallicity [Zsun]
    Real Z = 1./3;
    if (nscalars > 0) {
      Z = w0(m,nhydro,k,j,i)/Zsol; // Assumes Z is the first passive scalar
    }

    // Get azimuthal radius for heating calculations
    Real R = sqrt(x1v * x1v + x2v * x2v); 

    Real temp = 1.0; // temperature in cgs units
    if (use_e) {
      temp = temp_unit * w0(m,IEN,k,j,i) / w0(m,IDN,k,j,i) * gm1;
    } else {
      temp = temp_unit * w0(m,ITM,k,j,i);
    }
    Real nH = X * nH_unit * w0(m,IDN,k,j,i); // density in cgs units

    // Caculate PIE cooling
    // WiersmaCooling at redshift z = 0 taken from Wiersma et al (2009)
    Real lambda_cooling_PIE = 0.0;     
    // Ensure we are in range of cooling table
    if (temp > Tfloor && temp < Tceil && nH > nHfloor && nH < nHceil) {
      // Convert input values to log space
      Real log_temp = log10(temp);
      Real log_density = log10(nH);

      // Locate indices in Tbins and nHbins
      int i = 0, j = 0;
      while (i < Tbins_DIM_0 - 2 && log10(Tbins_(i + 1)) < log_temp) ++i; 
      while (j < nHbins_DIM_0 - 2 && log10(nHbins_(j + 1)) < log_density) ++j;

      // Logarithms of the Tbins and nHbins bounding the point
      Real log_T0 = log10(Tbins_(i));
      Real log_T1 = log10(Tbins_(i + 1));
      Real log_nH0 = log10(nHbins_(j));
      Real log_nH1 = log10(nHbins_(j + 1));

      // Compute weights for bilinear interpolation
      Real t = (log_temp - log_T0) / (log_T1 - log_T0);
      Real u = (log_density - log_nH0) / (log_nH1 - log_nH0);
  
      // Corner values in log space from the H_He cooling grid
      Real C00 = H_He_Cooling_(i, j);
      Real C10 = H_He_Cooling_(i + 1, j);
      Real C01 = H_He_Cooling_(i, j + 1);
      Real C11 = H_He_Cooling_(i + 1, j + 1);
 
      // Corner values in log space from the Metal cooling grid
      Real M00 = Metal_Cooling_(i, j);
      Real M10 = Metal_Cooling_(i + 1, j);
      Real M01 = Metal_Cooling_(i, j + 1);
      Real M11 = Metal_Cooling_(i + 1, j + 1);

      // Bilinear interpolation in log space
      Real prim_cooling =
          (1 - t) * (1 - u) * C00 +
          t * (1 - u) * C10 +
          (1 - t) * u * C01 +
          t * u * C11;
      Real metal_cooling =
          (1 - t) * (1 - u) * M00 +
          t * (1 - u) * M10 +
          (1 - t) * u * M01 +
          t * u * M11;
 
      lambda_cooling_PIE = prim_cooling + Z * metal_cooling;
    } 
    
    // Caculate CIE cooling
    Real lambda_cooling_CIE = 0.0;  
    if (temp > Tfloor && temp < Tceil) {
      // Convert input values to log space
      Real log_temp = log10(temp);

      // Locate index in Tbins
      int i = 0;
      while (i < Tbins_DIM_0 - 2 && log10(Tbins_(i + 1)) < log_temp) ++i;

      // Logarithms of the Tbins bounding the point
      Real log_T0 = log10(Tbins_(i));
      Real log_T1 = log10(Tbins_(i + 1));

      // Compute weight for linear interpolation
      Real t = (log_temp - log_T0) / (log_T1 - log_T0);

      // Bin values in log space from the H_He cooling grid
      Real C0 = H_He_Cooling_CIE_(i);
      Real C1 = H_He_Cooling_CIE_(i + 1);

      // Bin values in log space from the Metal cooling grid
      Real M0 = Metal_Cooling_CIE_(i);
      Real M1 = Metal_Cooling_CIE_(i + 1);

      // Linear Interpolation
      Real prim_cooling = C0 + t * (C1 - C0);
      Real metal_cooling = M0 + t * (M1 - M0);
      
      lambda_cooling_CIE = prim_cooling + Z * metal_cooling;
    }
    else if (temp < Tfloor) {
      // for temperatures less than 100 K, use Koyama & Inutsuka (2002)
      lambda_cooling_CIE = Z * (2.0e-19 * exp(-1.184e5 / (temp + 1.0e3))
                                + 2.8e-28 * sqrt(temp) * exp(-92.0 / temp));
    }

    // Calculate heating
    Real horz_falloff = exp(-R / h_radius);
    Real vert_falloff = exp(-(x3v*x3v) / sqrt(h_height*h_height + h_alpha*R*R));
    Real gamma_heating = h_rate * h_norm * X * nH_unit * horz_falloff * vert_falloff;
    if (temp > 1e4) gamma_heating *= pow(temp / 1e4, -8.0);

    // Combine CIE and PIE cooling
    Real dx = size.d_view(m).dx1 * length_unit;
    Real neutral_frac = 1 - (0.5 * (1 + tanh((temp - 8e3) / 1.5e3)));
    Real tau = neutral_frac * nH * 1e-17 * dx; // optical depth of cell
    Real frac = exp(-tau);
    Real lambda_cooling = (1 - frac) * lambda_cooling_CIE + frac * lambda_cooling_PIE;
    gamma_heating *= (1 - frac);
    
    // Add cooling and heating source terms to energy equation
    u0(m,IEN,k,j,i) -= bdt * X * w0(m,IDN,k,j,i) *
                       (X * w0(m,IDN,k,j,i) * lambda_cooling / cooling_unit 
                                             - gamma_heating / heating_unit);
  });

  return;
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
