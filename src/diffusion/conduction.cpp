//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file conduction.cpp
//! \brief Implements functions for Conduction class. This includes isotropic thermal
//! conduction, in which heat flux is proportional to negative local temperature gradient.
//! Conduction may be added to Hydro and/or MHD independently.

#include <float.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <string>
#include <iostream> // cout

// Athena++ headers
#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "eos/eos.hpp"
#include "conduction.hpp"
#include "units/units.hpp"

namespace {

parabolic::ParabolicIntegratorMode ParseConductivityIntegrator(std::string block,
    ParameterInput *pin) {
  std::string integrator = pin->GetOrAddString(block, "conductivity_integrator", "explicit");
  if (integrator == "explicit") {
    return parabolic::ParabolicIntegratorMode::explicit_mode;
  } else if (integrator == "sts") {
    return parabolic::ParabolicIntegratorMode::sts;
  }

  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
            << "<" << block << ">/conductivity_integrator = '" << integrator
            << "' must be 'explicit' or 'sts'" << std::endl;
  std::exit(EXIT_FAILURE);
}

} // namespace

// VanLeer Limiter which takes 2 slopes
KOKKOS_INLINE_FUNCTION
Real VLL2State(const Real a, const Real b) {
  if (a*b > 0) {
    return 2.0*a*b/(a+b);
  } else {
    return 0.0;
  }
}

// VanLeer Limiter which takes 4 slopes
KOKKOS_INLINE_FUNCTION
Real VLL4State(const Real a, const Real b, const Real c, const Real d) {
  return VLL2State(VLL2State(a,b), VLL2State(c,d));
}

KOKKOS_INLINE_FUNCTION
Real CGLBMag(const Real bx, const Real by, const Real bz) {
  return sqrt(SQR(bx) + SQR(by) + SQR(bz));
}

KOKKOS_INLINE_FUNCTION
Real CGLLimitedHeatFlux(const Real q_unlimited, const Real q_max) {
  if (q_max <= 0.0) {
    return 0.0;
  }
  return q_unlimited*q_max/(q_max + fabs(q_unlimited));
}

KOKKOS_INLINE_FUNCTION
Real CGLLimiterCollisionRate(const Real ppar, const Real pperp,
                             const Real bx, const Real by, const Real bz,
                             const Real lim_coll, const bool mlim,
                             const bool flim, const bool backup_lim) {
  const Real paniso = pperp - ppar;
  const Real bsqr = SQR(bx) + SQR(by) + SQR(bz);
  const Real limiter_nu = fmax(lim_coll, static_cast<Real>(0.0));
  const Real backup_nu = static_cast<Real>(1.0e10);
  Real nu_eff = 0.0;

  if (flim && backup_lim) {
    if ((paniso <= static_cast<Real>(-0.7)*bsqr) && (paniso > -bsqr)) {
      nu_eff = fmax(nu_eff, limiter_nu);
    } else if (paniso <= -bsqr) {
      nu_eff = fmax(nu_eff, backup_nu);
    }
  } else if (flim) {
    if (paniso <= static_cast<Real>(-0.7)*bsqr) {
      nu_eff = fmax(nu_eff, limiter_nu);
    }
  }

  if (mlim && backup_lim) {
    if ((paniso >= static_cast<Real>(0.5)*bsqr) && (paniso < bsqr)) {
      nu_eff = fmax(nu_eff, limiter_nu);
    } else if (paniso >= static_cast<Real>(0.5)*bsqr) {
      nu_eff = fmax(nu_eff, backup_nu);
    }
  } else if (mlim) {
    if (paniso >= static_cast<Real>(0.5)*bsqr) {
      nu_eff = fmax(nu_eff, limiter_nu);
    }
  }

  return nu_eff;
}

KOKKOS_INLINE_FUNCTION
Real CGLChiPerp(const Real cpar, const Real lf_k_parallel, const Real nu_eff) {
  const Real denom = static_cast<Real>(2.5066282746310002)*cpar*lf_k_parallel
                   + nu_eff;
  if (denom <= 0.0) {
    return 0.0;
  }
  return static_cast<Real>(2.0)*SQR(cpar)/denom;
}

KOKKOS_INLINE_FUNCTION
Real CGLChiParallel(const Real cpar, const Real lf_k_parallel, const Real nu_eff) {
  const Real denom = static_cast<Real>(5.013256549262000)*cpar*lf_k_parallel
                   + static_cast<Real>(1.4247779607693793)*nu_eff;
  if (denom <= 0.0) {
    return 0.0;
  }
  return static_cast<Real>(8.0)*SQR(cpar)/denom;
}

KOKKOS_INLINE_FUNCTION
void CGLLandauFluidFaceFlux(const Real rho_l, const Real rho_r,
                            const Real ppar_l, const Real ppar_r,
                            const Real pperp_l, const Real pperp_r,
                            const Real bx, const Real by, const Real bz,
                            const Real bh_dir,
                            const Real grad_tpar_x, const Real grad_tpar_y,
                            const Real grad_tpar_z,
                            const Real grad_tperp_x, const Real grad_tperp_y,
                            const Real grad_tperp_z,
                            const Real grad_bmag_x, const Real grad_bmag_y,
                            const Real grad_bmag_z,
                            const Real lf_k_parallel,
                            const bool lf_coeff_local,
                            const Real lf_c_parallel0,
                            const Real nu_coll,
                            const Real lim_coll,
                            const bool mlim,
                            const bool flim,
                            const bool backup_lim,
                            const Real dfloor,
                            const Real pfloor,
                            const Real tfloor,
                            const Real bfloor,
                            Real &energy_flux,
                            Real &moment_flux) {
  energy_flux = 0.0;
  moment_flux = 0.0;

  const Real bmag = CGLBMag(bx, by, bz);
  if (bmag <= bfloor || lf_k_parallel <= 0.0) {
    return;
  }

  const Real bhx = bx/bmag;
  const Real bhy = by/bmag;
  const Real bhz = bz/bmag;
  const Real rho = fmax(0.5*(rho_l + rho_r), dfloor);
  const Real ppar = fmax(0.5*(ppar_l + ppar_r), pfloor);
  const Real pperp = fmax(0.5*(pperp_l + pperp_r), pfloor);
  const Real cpar = lf_coeff_local ? sqrt(fmax(ppar/rho, tfloor)) : lf_c_parallel0;
  const Real nu_eff = fmax(nu_coll, static_cast<Real>(0.0)) +
      CGLLimiterCollisionRate(ppar, pperp, bx, by, bz, lim_coll, mlim, flim, backup_lim);
  const Real chi_perp = CGLChiPerp(cpar, lf_k_parallel, nu_eff);
  const Real chi_parallel = CGLChiParallel(cpar, lf_k_parallel, nu_eff);

  const Real gradpar_tpar = bhx*grad_tpar_x + bhy*grad_tpar_y + bhz*grad_tpar_z;
  const Real gradpar_tperp = bhx*grad_tperp_x + bhy*grad_tperp_y + bhz*grad_tperp_z;
  const Real gradpar_bmag = bhx*grad_bmag_x + bhy*grad_bmag_y + bhz*grad_bmag_z;
  const Real q_parallel_l = -chi_parallel*rho*gradpar_tpar;
  const Real q_perp_l = -chi_perp*(rho*gradpar_tperp -
                                   pperp*(1.0 - pperp/ppar)*gradpar_bmag/bmag);
  const Real q_parallel_max = static_cast<Real>(1.5957691216057308)*cpar*ppar;
  const Real q_perp_max = static_cast<Real>(0.7978845608028654)*cpar*pperp;

  // Squire et al. (2023), eqs. 3.1-3.2.
  const Real q_parallel = CGLLimitedHeatFlux(q_parallel_l, q_parallel_max);
  const Real q_perp = CGLLimitedHeatFlux(q_perp_l, q_perp_max);

  energy_flux = bh_dir*(q_perp + 0.5*q_parallel);
  moment_flux = bh_dir*q_perp/bmag;
}

//----------------------------------------------------------------------------------------
//! \fn Real TempDepKappa()
//! \brief Temperature-dependent conductivity given by Parker (1953) and Spitzer (1962)

KOKKOS_INLINE_FUNCTION
Real TempDepKappa(Real temp, Real limit) {
  if (temp < 6.5e4) {
    return 2.5e3 * pow(temp, 0.5);
  } else {
    return fmin(6.0e-7*pow(temp, 2.5), limit);
  }
}

//----------------------------------------------------------------------------------------
//! \brief Conduction constructor
// Note first argument passes string ("hydro" or "mhd") denoting in wihch class this
// object is being constructed, and therefore which <block> in the input file from which
// the parameters are read.
// Note that the coefficient of thermal conduction, kappa, corresponds to conductivity,
// not diffusivity. This is different from the coefficient used in Athena++.

Conduction::Conduction(std::string block, MeshBlockPack *pp, ParameterInput *pin) :
    pmy_pack(pp) {
  const bool has_iso_conduction = pin->DoesParameterExist(block,"isotropic_conduction");
  const bool has_cgl_heat_flux = pin->DoesParameterExist(block,"cgl_heat_flux");

  if (has_iso_conduction && has_cgl_heat_flux) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<" << block << "> cannot set both isotropic_conduction and "
              << "cgl_heat_flux" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Read parameters for isotropic thermal conduction (if any)
  if (has_iso_conduction) {
    iso_cond_type = pin->GetString(block,"isotropic_conduction");
    // Check for valid type
    if ((iso_cond_type.compare("constant") != 0) &&
        (iso_cond_type.compare("spitzer") != 0) &&
        (iso_cond_type.compare("spitzer_limited") != 0)) {
      std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__ << std::endl
                << "Invalid choice for isotropic thermal conduction type" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // constant conductivity
    if (iso_cond_type.compare("constant") == 0) {
      kappa_iso = pin->GetReal(block,"kappa_iso");
    }
    kappa_iso_limit = pin->GetOrAddReal(block,"kappa_iso_limit",
                      static_cast<Real>(std::numeric_limits<float>::max()));
    mode = ParseConductivityIntegrator(block, pin);
  }

  // Read parameters for CGL Landau-fluid heat flux (if any). This mode is intentionally
  // separate from ordinary isotropic conduction because the flux is field-aligned and
  // updates total energy plus magnetic moment during the STS sweep.
  if (has_cgl_heat_flux) {
    cgl_heat_flux_type = pin->GetString(block,"cgl_heat_flux");
    if (cgl_heat_flux_type.compare("landau_fluid") != 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Invalid choice for CGL heat flux: " << cgl_heat_flux_type
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    lf_k_parallel = pin->GetReal(block,"lf_k_parallel");
    if (lf_k_parallel <= 0.0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<" << block << ">/lf_k_parallel must be positive and is interpreted "
                << "as the Landau-fluid |k_parallel|"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    const std::string coeff_mode = pin->GetOrAddString(block, "lf_coefficient_mode",
                                                       "local");
    if (coeff_mode == "local") {
      lf_coeff_local = true;
    } else if (coeff_mode == "background") {
      lf_coeff_local = false;
      lf_c_parallel0 = pin->GetReal(block, "lf_c_parallel0");
      if (lf_c_parallel0 <= 0.0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "<" << block << ">/lf_c_parallel0 must be positive when "
                  << "lf_coefficient_mode = background" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<" << block << ">/lf_coefficient_mode = '" << coeff_mode
                << "' must be 'local' or 'background'" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    mode = ParseConductivityIntegrator(block, pin);
  }
}

//----------------------------------------------------------------------------------------
//! \brief Conduction destructor

Conduction::~Conduction() {
}

//----------------------------------------------------------------------------------------
//! \fn void AddHeatFluxes()
//! \brief Wrapper function that adds heat fluxes for different types of thermal
//! conduction to face-centered fluxes of conserved variables

void Conduction::AddHeatFluxes(const DvceArray5D<Real> &w0, const EOS_Data &eos,
    DvceFaceFld5D<Real> &flx) {
  if (IsCGLLandauFluidHeatFlux()) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "CGL Landau-fluid heat flux requires cell-centered magnetic fields. "
              << "Call AddCGLLandauFluidHeatFluxes." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (iso_cond_type.compare("constant") == 0) {
    AddIsotropicHeatFluxConstCond(w0, eos, flx);
  } else if ((iso_cond_type.compare("spitzer") == 0) ||
             (iso_cond_type.compare("spitzer_limited") == 0)) {
    AddIsotropicHeatFluxSpitzerCond(w0, eos, flx);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void AddCGLLandauFluidHeatFluxes()
//! \brief Adds field-aligned CGL Landau-fluid heat fluxes to total energy and magnetic
//! moment fluxes. The IAN slot is interpreted as p_perp/|B| only during the CGL STS
//! update; the hyperbolic solve continues to store conserved anisotropy A there.

void Conduction::AddCGLLandauFluidHeatFluxes(const DvceArray5D<Real> &w0,
    const DvceArray5D<Real> &bcc, const EOS_Data &eos, DvceFaceFld5D<Real> &flx) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto size = pmy_pack->pmb->mb_size;
  bool multi_d = pmy_pack->pmesh->multi_d;
  bool three_d = pmy_pack->pmesh->three_d;
  Real lf_k_parallel_ = lf_k_parallel;
  bool lf_coeff_local_ = lf_coeff_local;
  Real lf_c_parallel0_ = lf_c_parallel0;
  Real dfloor_ = eos.dfloor;
  Real pfloor_ = eos.pfloor;
  Real tfloor_ = eos.tfloor;
  Real bfloor_ = eos.bfloor;
  Real nu_coll_ = eos.nu_coll;
  Real lim_coll_ = eos.lim_coll;
  bool mlim_ = eos.mlim;
  bool flim_ = eos.flim;
  bool backup_lim_ = eos.backup_lim;

  if (lf_k_parallel_ <= 0.0) {
    return;
  }

  // fluxes in x1-direction
  auto &flx1 = flx.x1f;
  par_for("cgl_lf_hflux1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real tpl_l = w0(m,IPR,k,j,i-1)/w0(m,IDN,k,j,i-1);
    const Real tpl_r = w0(m,IPR,k,j,i  )/w0(m,IDN,k,j,i  );
    const Real tpp_l = w0(m,IPP,k,j,i-1)/w0(m,IDN,k,j,i-1);
    const Real tpp_r = w0(m,IPP,k,j,i  )/w0(m,IDN,k,j,i  );
    const Real bmag_l = CGLBMag(bcc(m,IBX,k,j,i-1), bcc(m,IBY,k,j,i-1),
                                bcc(m,IBZ,k,j,i-1));
    const Real bmag_r = CGLBMag(bcc(m,IBX,k,j,i), bcc(m,IBY,k,j,i),
                                bcc(m,IBZ,k,j,i));
    Real grad_pl_x = (tpl_r - tpl_l)/size.d_view(m).dx1;
    Real grad_pp_x = (tpp_r - tpp_l)/size.d_view(m).dx1;
    Real grad_b_x = (bmag_r - bmag_l)/size.d_view(m).dx1;
    Real grad_pl_y = 0.0, grad_pp_y = 0.0;
    Real grad_pl_z = 0.0, grad_pp_z = 0.0;
    Real grad_b_y = 0.0, grad_b_z = 0.0;

    if (multi_d) {
      const Real tpl_ju_i   = w0(m,IPR,k,j+1,i  )/w0(m,IDN,k,j+1,i  );
      const Real tpl_jl_i   = w0(m,IPR,k,j-1,i  )/w0(m,IDN,k,j-1,i  );
      const Real tpl_ju_im1 = w0(m,IPR,k,j+1,i-1)/w0(m,IDN,k,j+1,i-1);
      const Real tpl_jl_im1 = w0(m,IPR,k,j-1,i-1)/w0(m,IDN,k,j-1,i-1);
      const Real tpp_ju_i   = w0(m,IPP,k,j+1,i  )/w0(m,IDN,k,j+1,i  );
      const Real tpp_jl_i   = w0(m,IPP,k,j-1,i  )/w0(m,IDN,k,j-1,i  );
      const Real tpp_ju_im1 = w0(m,IPP,k,j+1,i-1)/w0(m,IDN,k,j+1,i-1);
      const Real tpp_jl_im1 = w0(m,IPP,k,j-1,i-1)/w0(m,IDN,k,j-1,i-1);
      const Real b_ju_i = CGLBMag(bcc(m,IBX,k,j+1,i), bcc(m,IBY,k,j+1,i),
                                  bcc(m,IBZ,k,j+1,i));
      const Real b_jl_i = CGLBMag(bcc(m,IBX,k,j-1,i), bcc(m,IBY,k,j-1,i),
                                  bcc(m,IBZ,k,j-1,i));
      const Real b_ju_im1 = CGLBMag(bcc(m,IBX,k,j+1,i-1),
                                    bcc(m,IBY,k,j+1,i-1),
                                    bcc(m,IBZ,k,j+1,i-1));
      const Real b_jl_im1 = CGLBMag(bcc(m,IBX,k,j-1,i-1),
                                    bcc(m,IBY,k,j-1,i-1),
                                    bcc(m,IBZ,k,j-1,i-1));
      grad_pl_y = 0.25*((tpl_ju_i - tpl_jl_i) + (tpl_ju_im1 - tpl_jl_im1))/
                  size.d_view(m).dx2;
      grad_pp_y = 0.25*((tpp_ju_i - tpp_jl_i) + (tpp_ju_im1 - tpp_jl_im1))/
                  size.d_view(m).dx2;
      grad_b_y = 0.25*((b_ju_i - b_jl_i) + (b_ju_im1 - b_jl_im1))/
                 size.d_view(m).dx2;
    }

    if (three_d) {
      const Real tpl_ku_i   = w0(m,IPR,k+1,j,i  )/w0(m,IDN,k+1,j,i  );
      const Real tpl_kl_i   = w0(m,IPR,k-1,j,i  )/w0(m,IDN,k-1,j,i  );
      const Real tpl_ku_im1 = w0(m,IPR,k+1,j,i-1)/w0(m,IDN,k+1,j,i-1);
      const Real tpl_kl_im1 = w0(m,IPR,k-1,j,i-1)/w0(m,IDN,k-1,j,i-1);
      const Real tpp_ku_i   = w0(m,IPP,k+1,j,i  )/w0(m,IDN,k+1,j,i  );
      const Real tpp_kl_i   = w0(m,IPP,k-1,j,i  )/w0(m,IDN,k-1,j,i  );
      const Real tpp_ku_im1 = w0(m,IPP,k+1,j,i-1)/w0(m,IDN,k+1,j,i-1);
      const Real tpp_kl_im1 = w0(m,IPP,k-1,j,i-1)/w0(m,IDN,k-1,j,i-1);
      const Real b_ku_i = CGLBMag(bcc(m,IBX,k+1,j,i), bcc(m,IBY,k+1,j,i),
                                  bcc(m,IBZ,k+1,j,i));
      const Real b_kl_i = CGLBMag(bcc(m,IBX,k-1,j,i), bcc(m,IBY,k-1,j,i),
                                  bcc(m,IBZ,k-1,j,i));
      const Real b_ku_im1 = CGLBMag(bcc(m,IBX,k+1,j,i-1),
                                    bcc(m,IBY,k+1,j,i-1),
                                    bcc(m,IBZ,k+1,j,i-1));
      const Real b_kl_im1 = CGLBMag(bcc(m,IBX,k-1,j,i-1),
                                    bcc(m,IBY,k-1,j,i-1),
                                    bcc(m,IBZ,k-1,j,i-1));
      grad_pl_z = 0.25*((tpl_ku_i - tpl_kl_i) + (tpl_ku_im1 - tpl_kl_im1))/
                  size.d_view(m).dx3;
      grad_pp_z = 0.25*((tpp_ku_i - tpp_kl_i) + (tpp_ku_im1 - tpp_kl_im1))/
                  size.d_view(m).dx3;
      grad_b_z = 0.25*((b_ku_i - b_kl_i) + (b_ku_im1 - b_kl_im1))/
                 size.d_view(m).dx3;
    }

    const Real bx = 0.5*(bcc(m,IBX,k,j,i-1) + bcc(m,IBX,k,j,i));
    const Real by = 0.5*(bcc(m,IBY,k,j,i-1) + bcc(m,IBY,k,j,i));
    const Real bz = 0.5*(bcc(m,IBZ,k,j,i-1) + bcc(m,IBZ,k,j,i));
    const Real bmag = CGLBMag(bx, by, bz);
    const Real bhx = (bmag > bfloor_) ? bx/bmag : 0.0;
    Real energy_flux = 0.0, moment_flux = 0.0;
    CGLLandauFluidFaceFlux(w0(m,IDN,k,j,i-1), w0(m,IDN,k,j,i),
                           w0(m,IPR,k,j,i-1), w0(m,IPR,k,j,i),
                           w0(m,IPP,k,j,i-1), w0(m,IPP,k,j,i),
                           bx, by, bz, bhx,
                           grad_pl_x, grad_pl_y, grad_pl_z,
                           grad_pp_x, grad_pp_y, grad_pp_z,
                           grad_b_x, grad_b_y, grad_b_z,
                           lf_k_parallel_, lf_coeff_local_, lf_c_parallel0_,
                           nu_coll_, lim_coll_, mlim_, flim_, backup_lim_,
                           dfloor_, pfloor_, tfloor_, bfloor_,
                           energy_flux, moment_flux);
    flx1(m,IEN,k,j,i) += energy_flux;
    flx1(m,IAN,k,j,i) += moment_flux;
  });
  if (pmy_pack->pmesh->one_d) {return;}

  // fluxes in x2-direction
  auto &flx2 = flx.x2f;
  par_for("cgl_lf_hflux2", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real tpl_l = w0(m,IPR,k,j-1,i)/w0(m,IDN,k,j-1,i);
    const Real tpl_r = w0(m,IPR,k,j  ,i)/w0(m,IDN,k,j  ,i);
    const Real tpp_l = w0(m,IPP,k,j-1,i)/w0(m,IDN,k,j-1,i);
    const Real tpp_r = w0(m,IPP,k,j  ,i)/w0(m,IDN,k,j  ,i);
    const Real bmag_l = CGLBMag(bcc(m,IBX,k,j-1,i), bcc(m,IBY,k,j-1,i),
                                bcc(m,IBZ,k,j-1,i));
    const Real bmag_r = CGLBMag(bcc(m,IBX,k,j,i), bcc(m,IBY,k,j,i),
                                bcc(m,IBZ,k,j,i));
    Real grad_pl_x = 0.0, grad_pp_x = 0.0;
    Real grad_pl_y = (tpl_r - tpl_l)/size.d_view(m).dx2;
    Real grad_pp_y = (tpp_r - tpp_l)/size.d_view(m).dx2;
    Real grad_b_y = (bmag_r - bmag_l)/size.d_view(m).dx2;
    Real grad_pl_z = 0.0, grad_pp_z = 0.0;
    Real grad_b_x = 0.0, grad_b_z = 0.0;

    const Real tpl_j_iup = w0(m,IPR,k,j  ,i+1)/w0(m,IDN,k,j  ,i+1);
    const Real tpl_j_i   = w0(m,IPR,k,j  ,i-1)/w0(m,IDN,k,j  ,i-1);
    const Real tpl_jm_iup = w0(m,IPR,k,j-1,i+1)/w0(m,IDN,k,j-1,i+1);
    const Real tpl_jm_i   = w0(m,IPR,k,j-1,i-1)/w0(m,IDN,k,j-1,i-1);
    const Real tpp_j_iup = w0(m,IPP,k,j  ,i+1)/w0(m,IDN,k,j  ,i+1);
    const Real tpp_j_i   = w0(m,IPP,k,j  ,i-1)/w0(m,IDN,k,j  ,i-1);
    const Real tpp_jm_iup = w0(m,IPP,k,j-1,i+1)/w0(m,IDN,k,j-1,i+1);
    const Real tpp_jm_i   = w0(m,IPP,k,j-1,i-1)/w0(m,IDN,k,j-1,i-1);
    const Real b_j_iup = CGLBMag(bcc(m,IBX,k,j,i+1), bcc(m,IBY,k,j,i+1),
                                 bcc(m,IBZ,k,j,i+1));
    const Real b_j_i = CGLBMag(bcc(m,IBX,k,j,i-1), bcc(m,IBY,k,j,i-1),
                               bcc(m,IBZ,k,j,i-1));
    const Real b_jm_iup = CGLBMag(bcc(m,IBX,k,j-1,i+1),
                                  bcc(m,IBY,k,j-1,i+1),
                                  bcc(m,IBZ,k,j-1,i+1));
    const Real b_jm_i = CGLBMag(bcc(m,IBX,k,j-1,i-1),
                                bcc(m,IBY,k,j-1,i-1),
                                bcc(m,IBZ,k,j-1,i-1));
    grad_pl_x = 0.25*((tpl_j_iup - tpl_j_i) + (tpl_jm_iup - tpl_jm_i))/
                size.d_view(m).dx1;
    grad_pp_x = 0.25*((tpp_j_iup - tpp_j_i) + (tpp_jm_iup - tpp_jm_i))/
                size.d_view(m).dx1;
    grad_b_x = 0.25*((b_j_iup - b_j_i) + (b_jm_iup - b_jm_i))/
               size.d_view(m).dx1;

    if (three_d) {
      const Real tpl_ku_j   = w0(m,IPR,k+1,j  ,i)/w0(m,IDN,k+1,j  ,i);
      const Real tpl_kl_j   = w0(m,IPR,k-1,j  ,i)/w0(m,IDN,k-1,j  ,i);
      const Real tpl_ku_jm1 = w0(m,IPR,k+1,j-1,i)/w0(m,IDN,k+1,j-1,i);
      const Real tpl_kl_jm1 = w0(m,IPR,k-1,j-1,i)/w0(m,IDN,k-1,j-1,i);
      const Real tpp_ku_j   = w0(m,IPP,k+1,j  ,i)/w0(m,IDN,k+1,j  ,i);
      const Real tpp_kl_j   = w0(m,IPP,k-1,j  ,i)/w0(m,IDN,k-1,j  ,i);
      const Real tpp_ku_jm1 = w0(m,IPP,k+1,j-1,i)/w0(m,IDN,k+1,j-1,i);
      const Real tpp_kl_jm1 = w0(m,IPP,k-1,j-1,i)/w0(m,IDN,k-1,j-1,i);
      const Real b_ku_j = CGLBMag(bcc(m,IBX,k+1,j,i), bcc(m,IBY,k+1,j,i),
                                  bcc(m,IBZ,k+1,j,i));
      const Real b_kl_j = CGLBMag(bcc(m,IBX,k-1,j,i), bcc(m,IBY,k-1,j,i),
                                  bcc(m,IBZ,k-1,j,i));
      const Real b_ku_jm1 = CGLBMag(bcc(m,IBX,k+1,j-1,i),
                                    bcc(m,IBY,k+1,j-1,i),
                                    bcc(m,IBZ,k+1,j-1,i));
      const Real b_kl_jm1 = CGLBMag(bcc(m,IBX,k-1,j-1,i),
                                    bcc(m,IBY,k-1,j-1,i),
                                    bcc(m,IBZ,k-1,j-1,i));
      grad_pl_z = 0.25*((tpl_ku_j - tpl_kl_j) + (tpl_ku_jm1 - tpl_kl_jm1))/
                  size.d_view(m).dx3;
      grad_pp_z = 0.25*((tpp_ku_j - tpp_kl_j) + (tpp_ku_jm1 - tpp_kl_jm1))/
                  size.d_view(m).dx3;
      grad_b_z = 0.25*((b_ku_j - b_kl_j) + (b_ku_jm1 - b_kl_jm1))/
                 size.d_view(m).dx3;
    }

    const Real bx = 0.5*(bcc(m,IBX,k,j-1,i) + bcc(m,IBX,k,j,i));
    const Real by = 0.5*(bcc(m,IBY,k,j-1,i) + bcc(m,IBY,k,j,i));
    const Real bz = 0.5*(bcc(m,IBZ,k,j-1,i) + bcc(m,IBZ,k,j,i));
    const Real bmag = CGLBMag(bx, by, bz);
    const Real bhy = (bmag > bfloor_) ? by/bmag : 0.0;
    Real energy_flux = 0.0, moment_flux = 0.0;
    CGLLandauFluidFaceFlux(w0(m,IDN,k,j-1,i), w0(m,IDN,k,j,i),
                           w0(m,IPR,k,j-1,i), w0(m,IPR,k,j,i),
                           w0(m,IPP,k,j-1,i), w0(m,IPP,k,j,i),
                           bx, by, bz, bhy,
                           grad_pl_x, grad_pl_y, grad_pl_z,
                           grad_pp_x, grad_pp_y, grad_pp_z,
                           grad_b_x, grad_b_y, grad_b_z,
                           lf_k_parallel_, lf_coeff_local_, lf_c_parallel0_,
                           nu_coll_, lim_coll_, mlim_, flim_, backup_lim_,
                           dfloor_, pfloor_, tfloor_, bfloor_,
                           energy_flux, moment_flux);
    flx2(m,IEN,k,j,i) += energy_flux;
    flx2(m,IAN,k,j,i) += moment_flux;
  });
  if (pmy_pack->pmesh->two_d) {return;}

  // fluxes in x3-direction
  auto &flx3 = flx.x3f;
  par_for("cgl_lf_hflux3", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real tpl_l = w0(m,IPR,k-1,j,i)/w0(m,IDN,k-1,j,i);
    const Real tpl_r = w0(m,IPR,k  ,j,i)/w0(m,IDN,k  ,j,i);
    const Real tpp_l = w0(m,IPP,k-1,j,i)/w0(m,IDN,k-1,j,i);
    const Real tpp_r = w0(m,IPP,k  ,j,i)/w0(m,IDN,k  ,j,i);
    const Real bmag_l = CGLBMag(bcc(m,IBX,k-1,j,i), bcc(m,IBY,k-1,j,i),
                                bcc(m,IBZ,k-1,j,i));
    const Real bmag_r = CGLBMag(bcc(m,IBX,k,j,i), bcc(m,IBY,k,j,i),
                                bcc(m,IBZ,k,j,i));
    Real grad_pl_x = 0.0, grad_pp_x = 0.0;
    Real grad_pl_y = 0.0, grad_pp_y = 0.0;
    Real grad_pl_z = (tpl_r - tpl_l)/size.d_view(m).dx3;
    Real grad_pp_z = (tpp_r - tpp_l)/size.d_view(m).dx3;
    Real grad_b_z = (bmag_r - bmag_l)/size.d_view(m).dx3;
    Real grad_b_x = 0.0, grad_b_y = 0.0;

    const Real tpl_k_iup = w0(m,IPR,k  ,j,i+1)/w0(m,IDN,k  ,j,i+1);
    const Real tpl_k_i   = w0(m,IPR,k  ,j,i-1)/w0(m,IDN,k  ,j,i-1);
    const Real tpl_km_iup = w0(m,IPR,k-1,j,i+1)/w0(m,IDN,k-1,j,i+1);
    const Real tpl_km_i   = w0(m,IPR,k-1,j,i-1)/w0(m,IDN,k-1,j,i-1);
    const Real tpp_k_iup = w0(m,IPP,k  ,j,i+1)/w0(m,IDN,k  ,j,i+1);
    const Real tpp_k_i   = w0(m,IPP,k  ,j,i-1)/w0(m,IDN,k  ,j,i-1);
    const Real tpp_km_iup = w0(m,IPP,k-1,j,i+1)/w0(m,IDN,k-1,j,i+1);
    const Real tpp_km_i   = w0(m,IPP,k-1,j,i-1)/w0(m,IDN,k-1,j,i-1);
    const Real b_k_iup = CGLBMag(bcc(m,IBX,k,j,i+1), bcc(m,IBY,k,j,i+1),
                                 bcc(m,IBZ,k,j,i+1));
    const Real b_k_i = CGLBMag(bcc(m,IBX,k,j,i-1), bcc(m,IBY,k,j,i-1),
                               bcc(m,IBZ,k,j,i-1));
    const Real b_km_iup = CGLBMag(bcc(m,IBX,k-1,j,i+1),
                                  bcc(m,IBY,k-1,j,i+1),
                                  bcc(m,IBZ,k-1,j,i+1));
    const Real b_km_i = CGLBMag(bcc(m,IBX,k-1,j,i-1),
                                bcc(m,IBY,k-1,j,i-1),
                                bcc(m,IBZ,k-1,j,i-1));
    grad_pl_x = 0.25*((tpl_k_iup - tpl_k_i) + (tpl_km_iup - tpl_km_i))/
                size.d_view(m).dx1;
    grad_pp_x = 0.25*((tpp_k_iup - tpp_k_i) + (tpp_km_iup - tpp_km_i))/
                size.d_view(m).dx1;
    grad_b_x = 0.25*((b_k_iup - b_k_i) + (b_km_iup - b_km_i))/
               size.d_view(m).dx1;

    if (multi_d) {
      const Real tpl_k_ju = w0(m,IPR,k  ,j+1,i)/w0(m,IDN,k  ,j+1,i);
      const Real tpl_k_jl = w0(m,IPR,k  ,j-1,i)/w0(m,IDN,k  ,j-1,i);
      const Real tpl_km_ju = w0(m,IPR,k-1,j+1,i)/w0(m,IDN,k-1,j+1,i);
      const Real tpl_km_jl = w0(m,IPR,k-1,j-1,i)/w0(m,IDN,k-1,j-1,i);
      const Real tpp_k_ju = w0(m,IPP,k  ,j+1,i)/w0(m,IDN,k  ,j+1,i);
      const Real tpp_k_jl = w0(m,IPP,k  ,j-1,i)/w0(m,IDN,k  ,j-1,i);
      const Real tpp_km_ju = w0(m,IPP,k-1,j+1,i)/w0(m,IDN,k-1,j+1,i);
      const Real tpp_km_jl = w0(m,IPP,k-1,j-1,i)/w0(m,IDN,k-1,j-1,i);
      const Real b_k_ju = CGLBMag(bcc(m,IBX,k,j+1,i), bcc(m,IBY,k,j+1,i),
                                  bcc(m,IBZ,k,j+1,i));
      const Real b_k_jl = CGLBMag(bcc(m,IBX,k,j-1,i), bcc(m,IBY,k,j-1,i),
                                  bcc(m,IBZ,k,j-1,i));
      const Real b_km_ju = CGLBMag(bcc(m,IBX,k-1,j+1,i),
                                   bcc(m,IBY,k-1,j+1,i),
                                   bcc(m,IBZ,k-1,j+1,i));
      const Real b_km_jl = CGLBMag(bcc(m,IBX,k-1,j-1,i),
                                   bcc(m,IBY,k-1,j-1,i),
                                   bcc(m,IBZ,k-1,j-1,i));
      grad_pl_y = 0.25*((tpl_k_ju - tpl_k_jl) + (tpl_km_ju - tpl_km_jl))/
                  size.d_view(m).dx2;
      grad_pp_y = 0.25*((tpp_k_ju - tpp_k_jl) + (tpp_km_ju - tpp_km_jl))/
                  size.d_view(m).dx2;
      grad_b_y = 0.25*((b_k_ju - b_k_jl) + (b_km_ju - b_km_jl))/
                 size.d_view(m).dx2;
    }

    const Real bx = 0.5*(bcc(m,IBX,k-1,j,i) + bcc(m,IBX,k,j,i));
    const Real by = 0.5*(bcc(m,IBY,k-1,j,i) + bcc(m,IBY,k,j,i));
    const Real bz = 0.5*(bcc(m,IBZ,k-1,j,i) + bcc(m,IBZ,k,j,i));
    const Real bmag = CGLBMag(bx, by, bz);
    const Real bhz = (bmag > bfloor_) ? bz/bmag : 0.0;
    Real energy_flux = 0.0, moment_flux = 0.0;
    CGLLandauFluidFaceFlux(w0(m,IDN,k-1,j,i), w0(m,IDN,k,j,i),
                           w0(m,IPR,k-1,j,i), w0(m,IPR,k,j,i),
                           w0(m,IPP,k-1,j,i), w0(m,IPP,k,j,i),
                           bx, by, bz, bhz,
                           grad_pl_x, grad_pl_y, grad_pl_z,
                           grad_pp_x, grad_pp_y, grad_pp_z,
                           grad_b_x, grad_b_y, grad_b_z,
                           lf_k_parallel_, lf_coeff_local_, lf_c_parallel0_,
                           nu_coll_, lim_coll_, mlim_, flim_, backup_lim_,
                           dfloor_, pfloor_, tfloor_, bfloor_,
                           energy_flux, moment_flux);
    flx3(m,IEN,k,j,i) += energy_flux;
    flx3(m,IAN,k,j,i) += moment_flux;
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void AddIsotropicHeatFluxConstCond()
//! \brief Adds isotropic heat flux computed using constant conductivity to face-centered
//! fluxes of conserved variables

void Conduction::AddIsotropicHeatFluxConstCond(const DvceArray5D<Real> &w0,
    const EOS_Data &eos, DvceFaceFld5D<Real> &flx) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto size = pmy_pack->pmb->mb_size;
  bool cgl = eos.is_cgl;
  Real gm1 = eos.gamma-1.0;
  Real &kappa_ = kappa_iso;

  // fluxes in x1-direction
  auto &flx1 = flx.x1f;
  par_for("conduct1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (cgl) {
      // Legacy scalar pressure-temperature prototype. The MHD constructor fences this
      // path for CGL so it is not confused with Landau-fluid heat flux.
      Real dtempdxprl = (w0(m,IPR,k,j,i)/w0(m,IDN,k,j,i) -
                         w0(m,IPR,k,j,i-1)/w0(m,IDN,k,j,i-1)) / size.d_view(m).dx1;
      Real dtempdxprp = (w0(m,IPP,k,j,i)/w0(m,IDN,k,j,i) -
                         w0(m,IPP,k,j,i-1)/w0(m,IDN,k,j,i-1)) / size.d_view(m).dx1;
      flx1(m,IPR,k,j,i) -= kappa_ * dtempdxprl;
      flx1(m,IPP,k,j,i) -= kappa_ * dtempdxprp;
    } else {
      Real dtempdx = (w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i) -
                      w0(m,IEN,k,j,i-1)/w0(m,IDN,k,j,i-1)) * gm1 /
                     size.d_view(m).dx1;
      flx1(m,IEN,k,j,i) -= kappa_ * dtempdx;
    }
  });
  if (pmy_pack->pmesh->one_d) {return;}

  // fluxes in x2-direction
  auto &flx2 = flx.x2f;
  par_for("conduct2",DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (cgl) {
      // Legacy scalar pressure-temperature prototype; fenced for CGL at runtime.
      Real dtempdxprl = (w0(m,IPR,k,j,i)/w0(m,IDN,k,j,i) -
                         w0(m,IPR,k,j-1,i)/w0(m,IDN,k,j-1,i)) / size.d_view(m).dx2;
      Real dtempdxprp = (w0(m,IPP,k,j,i)/w0(m,IDN,k,j,i) -
                         w0(m,IPP,k,j-1,i)/w0(m,IDN,k,j-1,i)) / size.d_view(m).dx2;
      flx2(m,IPR,k,j,i) -= kappa_ * dtempdxprl;
      flx2(m,IPP,k,j,i) -= kappa_ * dtempdxprp;
    } else {
      Real dtempdx = (w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i) -
                      w0(m,IEN,k,j-1,i)/w0(m,IDN,k,j-1,i)) * gm1 /
                     size.d_view(m).dx2;
      flx2(m,IEN,k,j,i) -= kappa_ * dtempdx;
    }
  });
  if (pmy_pack->pmesh->two_d) {return;}

  // fluxes in x3-direction
  auto &flx3 = flx.x3f;
  par_for("conduct3",DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (cgl) {
      // Legacy scalar pressure-temperature prototype; fenced for CGL at runtime.
      Real dtempdxprl = (w0(m,IPR,k,j,i)/w0(m,IDN,k,j,i) -
                         w0(m,IPR,k-1,j,i)/w0(m,IDN,k-1,j,i)) / size.d_view(m).dx3;
      Real dtempdxprp = (w0(m,IPP,k,j,i)/w0(m,IDN,k,j,i) -
                         w0(m,IPP,k-1,j,i)/w0(m,IDN,k-1,j,i)) / size.d_view(m).dx3;
      flx3(m,IPR,k,j,i) -= kappa_ * dtempdxprl;
      flx3(m,IPP,k,j,i) -= kappa_ * dtempdxprp;
    } else {
      Real dtempdx = (w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i) -
                      w0(m,IEN,k-1,j,i)/w0(m,IDN,k-1,j,i)) * gm1 /
                     size.d_view(m).dx3;
      flx3(m,IEN,k,j,i) -= kappa_ * dtempdx;
    }
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void TempDependentHeatFlux()
//! \brief Adds heat flux to face-centered fluxes of conserved variables with
//! temperature-dependent conductivity

void Conduction::AddIsotropicHeatFluxSpitzerCond(const DvceArray5D<Real> &w0,
    const EOS_Data &eos, DvceFaceFld5D<Real> &flx) {
/*
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto size = pmy_pack->pmb->mb_size;
  const bool &sat_hflux_ = sat_hflux;
  bool &multi_d = pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_pack->pmesh->three_d;
  Real gm1 = eos.gamma-1.0;
  Real kappaceil = kappa_ceiling;
  Real temp_unit = pmy_pack->punit->temperature_cgs();
  Real kappa_unit = pmy_pack->punit->pressure_cgs()*pmy_pack->punit->velocity_cgs()*
                    pmy_pack->punit->length_cgs()/pmy_pack->punit->temperature_cgs();

  // fluxes in x1-direction
  auto &flx1 = flx.x1f;
  par_for("conduct1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // Add heat fluxes into fluxes of conserved variables: energy
    Real temp_l = w0(m,IEN,k,j,i-1)/w0(m,IDN,k,j,i-1)*gm1;
    Real temp_r = w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1;
    Real pres_l = w0(m,IEN,k,j,i-1)*gm1;
    Real pres_r = w0(m,IEN,k,j,i)*gm1;
    Real kappaf = 0.5*(TempDepKappa(temp_unit*temp_l,kappaceil)+
                  TempDepKappa(temp_unit*temp_r,kappaceil))/kappa_unit;
    Real dtempdx1 = (temp_r-temp_l)/size.d_view(m).dx1;
    Real hflx = kappaf*dtempdx1;
    // Saturation of thermal conduction by harmonic mean
    if (sat_hflux_) {
      Real dtempdx2 = 0.0, dtempdx3 = 0.0;
      if (multi_d) {
        temp_ll = w0(m,IEN,k,j-1,i-1)/w0(m,IDN,k,j-1,i-1)*gm1;
        temp_lr = w0(m,IEN,k,j+1,i-1)/w0(m,IDN,k,j+1,i-1)*gm1;
        temp_rl = w0(m,IEN,k,j-1,i)/w0(m,IDN,k,j-1,i)*gm1;
        temp_rr = w0(m,IEN,k,j+1,i)/w0(m,IDN,k,j+1,i)*gm1;
        dtempdx2 = VanLeerLimiter4State(temp_rr-temp_r,temp_r-temp_rl,
                                        temp_lr-temp_l,temp_l-temp_ll)/size.d_view(m).dx2;
      }
      if (three_d) {
        temp_ll = w0(m,IEN,k-1,j,i-1)/w0(m,IDN,k-1,j,i-1)*gm1;
        temp_lr = w0(m,IEN,k+1,j,i-1)/w0(m,IDN,k+1,j,i-1)*gm1;
        temp_rl = w0(m,IEN,k-1,j,i)/w0(m,IDN,k-1,j,i)*gm1;
        temp_rr = w0(m,IEN,k+1,j,i)/w0(m,IDN,k+1,j,i)*gm1;
        dtempdx3 = VL4Limiter(temp_rr-temp_r,temp_r-temp_rl,
                              temp_lr-temp_l,temp_l-temp_ll)/size.d_view(m).dx3;
      }
      Real tempgrad = sqrt(SQR(dtempdx1)+SQR(dtempdx2)+SQR(dtempdx3));
      Real pres_cs = 0.5*(pres_l*sqrt(temp_l)+pres_r*sqrt(temp_r));
      Real sat_fac = 1.0/(1.0+kappaf*tempgrad/(1.5*pres_cs));
      hflx *= sat_fac;
    }
    flx1(m,IEN,k,j,i) -= hflx;
  });
  if (pmy_pack->pmesh->one_d) {return;}

  // fluxes in x2-direction
  auto &flx2 = flx.x2f;
  par_for("conduct2",DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // Add heat fluxes into fluxes of conserved variables: energy
    Real temp_l = 0.0, temp_r = 0.0, pres_l = 0.0, pres_r = 0.0;
    Real temp_ll = 0.0, temp_lr = 0.0, temp_rl = 0.0, temp_rr = 0.0;
    temp_l = w0(m,IEN,k,j-1,i)/w0(m,IDN,k,j-1,i)*gm1;
    temp_r = w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1;
    pres_l = w0(m,IEN,k,j-1,i)*gm1;
    pres_r = w0(m,IEN,k,j,i)*gm1;
    Real kappaf = 0.5*(TempDepKappa(temp_unit*temp_l,kappaceil)+
                  TempDepKappa(temp_unit*temp_r,kappaceil))/kappa_unit;
    Real dtempdx2 = (temp_r-temp_l)/size.d_view(m).dx2;
    Real hflx = kappaf*dtempdx2;
    // Saturation of thermal conduction
    if (sat_hflux_) {
      Real dtempdx1 = 0.0, dtempdx3 = 0.0;
      temp_ll = w0(m,IEN,k,j-1,i-1)/w0(m,IDN,k,j-1,i-1)*gm1;
      temp_lr = w0(m,IEN,k,j-1,i+1)/w0(m,IDN,k,j-1,i+1)*gm1;
      temp_rl = w0(m,IEN,k,j,i-1)/w0(m,IDN,k,j,i-1)*gm1;
      temp_rr = w0(m,IEN,k,j,i+1)/w0(m,IDN,k,j,i+1)*gm1;
      dtempdx1 = VL4Limiter(temp_rr-temp_r,temp_r-temp_rl,
                            temp_lr-temp_l,temp_l-temp_ll)/size.d_view(m).dx1;
      if (three_d) {
        temp_ll = w0(m,IEN,k-1,j-1,i)/w0(m,IDN,k-1,j-1,i)*gm1;
        temp_lr = w0(m,IEN,k+1,j-1,i)/w0(m,IDN,k+1,j-1,i)*gm1;
        temp_rl = w0(m,IEN,k-1,j,i)/w0(m,IDN,k-1,j,i)*gm1;
        temp_rr = w0(m,IEN,k+1,j,i)/w0(m,IDN,k+1,j,i)*gm1;
        dtempdx3 = VL4Limiter(temp_rr-temp_r,temp_r-temp_rl,
                              temp_lr-temp_l,temp_l-temp_ll)/size.d_view(m).dx3;
      }
      Real tempgrad = sqrt(SQR(dtempdx1)+SQR(dtempdx2)+SQR(dtempdx3));
      Real pres_cs = 0.5*(pres_l*sqrt(temp_l)+pres_r*sqrt(temp_r));
      Real sat_fac = 1.0/(1.0+kappaf*tempgrad/(1.5*pres_cs));
      hflx *= sat_fac;
    }
    flx2(m,IEN,k,j,i) -= hflx;
  });
  if (pmy_pack->pmesh->two_d) {return;}

  // fluxes in x3-direction
  auto &flx3 = flx.x3f;
  par_for("conduct3",DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // Add heat fluxes into fluxes of conserved variables: energy
    Real temp_l = 0.0, temp_r = 0.0, pres_l = 0.0, pres_r = 0.0;
    Real temp_ll = 0.0, temp_lr = 0.0, temp_rl = 0.0, temp_rr = 0.0;
    temp_l = w0(m,IEN,k-1,j,i)/w0(m,IDN,k-1,j,i)*gm1;
    temp_r = w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1;
    pres_l = w0(m,IEN,k-1,j,i)*gm1;
    pres_r = w0(m,IEN,k,j,i)*gm1;
    Real kappaf = 0.5*(TempDepKappa(temp_unit*temp_l,kappaceil)+
                  TempDepKappa(temp_unit*temp_r,kappaceil))/kappa_unit;
    Real dtempdx3 = (temp_r-temp_l)/size.d_view(m).dx3;
    Real hflx = kappaf*dtempdx3;
    // Saturation of thermal conduction
    if (sat_hflux_) {
      Real dtempdx1 = 0.0, dtempdx2 = 0.0;
      temp_ll = w0(m,IEN,k-1,j,i-1)/w0(m,IDN,k-1,j,i-1)*gm1;
      temp_lr = w0(m,IEN,k-1,j,i+1)/w0(m,IDN,k-1,j,i+1)*gm1;
      temp_rl = w0(m,IEN,k,j,i-1)/w0(m,IDN,k,j,i-1)*gm1;
      temp_rr = w0(m,IEN,k,j,i+1)/w0(m,IDN,k,j,i+1)*gm1;
      dtempdx1 = VL4Limiter(temp_rr-temp_r,temp_r-temp_rl,
                            temp_lr-temp_l,temp_l-temp_ll)/size.d_view(m).dx1;
      temp_ll = w0(m,IEN,k-1,j-1,i)/w0(m,IDN,k-1,j-1,i)*gm1;
      temp_lr = w0(m,IEN,k-1,j+1,i)/w0(m,IDN,k-1,j+1,i)*gm1;
      temp_rl = w0(m,IEN,k,j-1,i)/w0(m,IDN,k,j-1,i)*gm1;
      temp_rr = w0(m,IEN,k,j+1,i)/w0(m,IDN,k,j+1,i)*gm1;
      dtempdx2 = VL4Limiter(temp_rr-temp_r,temp_r-temp_rl,
                            temp_lr-temp_l,temp_l-temp_ll)/size.d_view(m).dx2;
      Real tempgrad = sqrt(SQR(dtempdx1)+SQR(dtempdx2)+SQR(dtempdx3));
      Real pres_cs = 0.5*(pres_l*sqrt(temp_l)+pres_r*sqrt(temp_r));
      Real sat_fac = 1.0/(1.0+kappaf*tempgrad/(1.5*pres_cs));
      hflx *= sat_fac;
    }
    flx3(m,IEN,k,j,i) -= hflx;
  });

*/
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Conduction::NewTimeStep()
//! \brief Compute new time step for thermal conduction.

void Conduction::NewTimeStep(const DvceArray5D<Real> &w0, const EOS_Data &eos_data) {
  dtnew = static_cast<Real>(std::numeric_limits<float>::max());
  Real fac;
  if (pmy_pack->pmesh->three_d) {
    fac = 1.0/6.0;
  } else if (pmy_pack->pmesh->two_d) {
    fac = 0.25;
  } else {
    fac = 0.5;
  }
//  if (sat_hflux == true) {
//    dtnew = static_cast<Real>(std::numeric_limits<float>::max());
//    return;
//  }

  // set flag for Spitzer conductivity
  bool spitzer = false;
  if (!IsCGLLandauFluidHeatFlux() &&
      (iso_cond_type.compare("spitzer") == 0 ||
       iso_cond_type.compare("spitzer_limited") == 0)) {
    spitzer = true;
  }
  Real limit_ = kappa_iso_limit;
  Real temp_unit=0.0, kappa_unit=0.0;

  if (spitzer) {
    temp_unit = pmy_pack->punit->temperature_cgs();
    kappa_unit = pmy_pack->punit->pressure_cgs()*pmy_pack->punit->velocity_cgs()*
                 pmy_pack->punit->length_cgs()/pmy_pack->punit->temperature_cgs();
  }

  // capture variables for kernel
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, nx1 = indcs.nx1;
  int js = indcs.js, nx2 = indcs.nx2;
  int ks = indcs.ks, nx3 = indcs.nx3;
  const int nmkji = (pmy_pack->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  auto &w0_ = w0;
  auto &multi_d = pmy_pack->pmesh->multi_d;
  auto &three_d = pmy_pack->pmesh->three_d;
  auto &size = pmy_pack->pmb->mb_size;

  if (IsCGLLandauFluidHeatFlux()) {
    Real lf_k_parallel_ = lf_k_parallel;
    bool lf_coeff_local_ = lf_coeff_local;
    Real lf_c_parallel0_ = lf_c_parallel0;
    Real tfloor_ = eos_data.tfloor;
    if (lf_k_parallel_ <= 0.0) {
      return;
    }
    Kokkos::parallel_reduce("cgl_lf_newdt", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
    KOKKOS_LAMBDA(const int &idx, Real &min_dt) {
      int m = (idx)/nkji;
      int k = (idx - m*nkji)/nji;
      int j = (idx - m*nkji - k*nji)/nx1;
      int i = (idx - m*nkji - k*nji - j*nx1) + is;
      k += ks;
      j += js;

      const Real cpar = lf_coeff_local_ ?
                        sqrt(fmax(w0_(m,IPR,k,j,i)/w0_(m,IDN,k,j,i), tfloor_)) :
                        lf_c_parallel0_;
      const Real chi = static_cast<Real>(1.5957691216057308)*cpar/lf_k_parallel_;
      if (chi <= 0.0) {
        return;
      }
      min_dt = fmin(min_dt, SQR(size.d_view(m).dx1)/chi);
      if (multi_d) {
        min_dt = fmin(min_dt, SQR(size.d_view(m).dx2)/chi);
      }
      if (three_d) {
        min_dt = fmin(min_dt, SQR(size.d_view(m).dx3)/chi);
      }
    }, Kokkos::Min<Real>(dtnew));
    dtnew *= fac;
    return;
  }

  Real gm1 = eos_data.gamma-1.0;
  Real pressure_to_temperature = (eos_data.is_cgl) ? 1.0 : gm1;
  Real kappa0 = kappa_iso;

  if (kappa0 <= 0.0) {
    return;
  }

  // find smallest timestep for thermal conduction in each cell
  // Note loop over all cells needed even for constant conductivity
  Kokkos::parallel_reduce("cond_newdt", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, Real &min_dt) {
    // compute m,k,j,i indices of thread and call function
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;

    Real kappa_ = kappa0;
    if (spitzer) {
      Real temp = w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*pressure_to_temperature;
      kappa_ = TempDepKappa(temp*temp_unit, limit_)/kappa_unit;
    }

    min_dt = fmin(min_dt, SQR(size.d_view(m).dx1)/kappa_*w0_(m,IDN,k,j,i)/
                  pressure_to_temperature);
    if (multi_d) {
      min_dt = fmin(min_dt, SQR(size.d_view(m).dx2)/kappa_*w0_(m,IDN,k,j,i)/
                    pressure_to_temperature);
    }
    if (three_d) {
      min_dt = fmin(min_dt, SQR(size.d_view(m).dx3)/kappa_*w0_(m,IDN,k,j,i)/
                    pressure_to_temperature);
    }
  }, Kokkos::Min<Real>(dtnew));
  dtnew *= fac;

  return;
}
