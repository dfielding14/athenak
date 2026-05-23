//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_landau_fluid.cpp
//! \brief CGL Landau-fluid heat-flux closure and parabolic timestep bound.

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "diffusion/cgl_landau_fluid.hpp"

namespace {

parabolic::ParabolicIntegratorMode ParseCGLHeatFluxIntegrator(ParameterInput *pin) {
  const std::string integrator =
      pin->GetOrAddString("mhd", "cgl_heat_flux_integrator", "sts");
  if (integrator == "sts") {
    return parabolic::ParabolicIntegratorMode::sts;
  }
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
            << "<mhd>/cgl_heat_flux_integrator = '" << integrator
            << "' is not implemented; the initial CGL Landau-fluid feature is STS-only."
            << std::endl;
  std::exit(EXIT_FAILURE);
}

KOKKOS_INLINE_FUNCTION
Real CGLLimiterCollisionRate(const Real ppar, const Real pperp, const Real bsqr,
                             const Real limiter_rate, const bool mirror,
                             const bool firehose, const bool backup) {
  const Real paniso = pperp - ppar;
  const Real rate = fmax(limiter_rate, static_cast<Real>(0.0));
  const Real backup_rate = static_cast<Real>(1.0e10);
  Real nu = 0.0;
  if (firehose) {
    if (backup && paniso <= -bsqr) {
      nu = backup_rate;
    } else if (paniso <= static_cast<Real>(-0.7)*bsqr) {
      nu = rate;
    }
  }
  if (mirror) {
    if (backup && paniso >= bsqr) {
      nu = fmax(nu, backup_rate);
    } else if (paniso >= static_cast<Real>(0.5)*bsqr) {
      nu = fmax(nu, rate);
    }
  }
  return nu;
}

KOKKOS_INLINE_FUNCTION
Real CGLLimitedHeatFlux(const Real q, const Real qmax) {
  return (qmax > 0.0) ? q*qmax/(qmax + fabs(q)) : 0.0;
}

struct CGLLFFaceState {
  Real rho;
  Real ppar;
  Real pperp;
  Real bmag;
  Real bhx;
  Real bhy;
  Real bhz;
  Real bhdir;
  Real chi_parallel;
  Real chi_perp;
  Real cparallel;
};

KOKKOS_INLINE_FUNCTION
bool BuildCGLLFFaceState(const Real rho_l, const Real rho_r,
                         const Real ppar_l, const Real ppar_r,
                         const Real pperp_l, const Real pperp_r,
                         const Real bx, const Real by, const Real bz, const int dir,
                         const Real lf_k, const bool coeff_local, const Real cparallel0,
                         const EOS_Data &eos, CGLLFFaceState &face) {
  const Real bsqr = SQR(bx) + SQR(by) + SQR(bz);
  const Real bmag = sqrt(bsqr);
  if (bmag <= eos.bfloor || lf_k <= 0.0) {
    return false;
  }
  face.rho = fmax(static_cast<Real>(0.5)*(rho_l + rho_r), eos.dfloor);
  face.ppar = fmax(static_cast<Real>(0.5)*(ppar_l + ppar_r), eos.pfloor);
  face.pperp = fmax(static_cast<Real>(0.5)*(pperp_l + pperp_r), eos.pfloor);
  face.bmag = bmag;
  face.bhx = bx/bmag;
  face.bhy = by/bmag;
  face.bhz = bz/bmag;
  face.bhdir = (dir == 0) ? face.bhx : ((dir == 1) ? face.bhy : face.bhz);
  face.cparallel =
      coeff_local ? sqrt(fmax(face.ppar/face.rho, eos.tfloor)) : cparallel0;
  const Real nu = fmax(eos.nu_coll, static_cast<Real>(0.0)) +
      CGLLimiterCollisionRate(face.ppar, face.pperp, bsqr, eos.lim_coll,
                              eos.mlim, eos.flim, eos.backup_lim);
  const Real denom_perp = static_cast<Real>(2.5066282746310002)*face.cparallel*lf_k + nu;
  const Real denom_parallel = static_cast<Real>(5.013256549262000)*face.cparallel*lf_k
                            + static_cast<Real>(1.4247779607693793)*nu;
  face.chi_perp = (denom_perp > 0.0) ? 2.0*SQR(face.cparallel)/denom_perp : 0.0;
  face.chi_parallel =
      (denom_parallel > 0.0) ? 8.0*SQR(face.cparallel)/denom_parallel : 0.0;
  return true;
}

KOKKOS_INLINE_FUNCTION
void CGLLFFlux(const CGLLFFaceState &face, const Real gtpar_x, const Real gtpar_y,
               const Real gtpar_z, const Real gtperp_x, const Real gtperp_y,
               const Real gtperp_z, const Real gb_x, const Real gb_y, const Real gb_z,
               Real &eflux, Real &muflux) {
  const Real grad_tpar =
      face.bhx*gtpar_x + face.bhy*gtpar_y + face.bhz*gtpar_z;
  const Real grad_tperp =
      face.bhx*gtperp_x + face.bhy*gtperp_y + face.bhz*gtperp_z;
  const Real grad_b = face.bhx*gb_x + face.bhy*gb_y + face.bhz*gb_z;
  const Real qpar_l = -face.chi_parallel*face.rho*grad_tpar;
  const Real qperp_l = -face.chi_perp*(face.rho*grad_tperp -
      face.pperp*(1.0 - face.pperp/face.ppar)*grad_b/face.bmag);
  const Real qpar = CGLLimitedHeatFlux(
      qpar_l, static_cast<Real>(1.5957691216057308)*face.cparallel*face.ppar);
  const Real qperp = CGLLimitedHeatFlux(
      qperp_l, static_cast<Real>(0.7978845608028654)*face.cparallel*face.pperp);
  eflux = face.bhdir*(qperp + static_cast<Real>(0.5)*qpar);
  muflux = face.bhdir*qperp/face.bmag;
}

} // namespace

CGLLandauFluid::CGLLandauFluid(MeshBlockPack *pp, ParameterInput *pin) :
    dtnew(static_cast<Real>(std::numeric_limits<float>::max())),
    lf_k_parallel(0.0),
    lf_coeff_local(true),
    lf_c_parallel0(0.0),
    mode(parabolic::ParabolicIntegratorMode::sts),
    pmy_pack(pp),
    tpar_("cgl_lf_tpar", 1, 1, 1, 1),
    tperp_("cgl_lf_tperp", 1, 1, 1, 1),
    bmag_("cgl_lf_bmag", 1, 1, 1, 1) {
  const std::string model = pin->GetString("mhd", "cgl_heat_flux");
  if (model != "landau_fluid") {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<mhd>/cgl_heat_flux = '" << model << "' must be 'landau_fluid'"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  lf_k_parallel = pin->GetReal("mhd", "lf_k_parallel");
  if (lf_k_parallel <= 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<mhd>/lf_k_parallel must be positive." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const std::string coeff_mode =
      pin->GetOrAddString("mhd", "lf_coefficient_mode", "local");
  if (coeff_mode == "background") {
    lf_coeff_local = false;
    lf_c_parallel0 = pin->GetReal("mhd", "lf_c_parallel0");
    if (lf_c_parallel0 <= 0.0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<mhd>/lf_c_parallel0 must be positive." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else if (coeff_mode != "local") {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<mhd>/lf_coefficient_mode = '" << coeff_mode
              << "' must be 'local' or 'background'." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  mode = ParseCGLHeatFluxIntegrator(pin);
}

void CGLLandauFluid::AddHeatFluxes(const DvceArray5D<Real> &w,
                                   const DvceArray5D<Real> &bcc,
                                   const EOS_Data &eos, DvceFaceFld5D<Real> &f) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int ncells1 = indcs.nx1 + 2*indcs.ng;
  const int ncells2 = (indcs.nx2 > 1) ? indcs.nx2 + 2*indcs.ng : 1;
  const int ncells3 = (indcs.nx3 > 1) ? indcs.nx3 + 2*indcs.ng : 1;
  if (tpar_.extent(0) != static_cast<std::size_t>(pmy_pack->nmb_thispack) ||
      tpar_.extent(1) != static_cast<std::size_t>(ncells3) ||
      tpar_.extent(2) != static_cast<std::size_t>(ncells2) ||
      tpar_.extent(3) != static_cast<std::size_t>(ncells1)) {
    Kokkos::realloc(tpar_, pmy_pack->nmb_thispack, ncells3, ncells2, ncells1);
    Kokkos::realloc(tperp_, pmy_pack->nmb_thispack, ncells3, ncells2, ncells1);
    Kokkos::realloc(bmag_, pmy_pack->nmb_thispack, ncells3, ncells2, ncells1);
  }

  auto tpar = tpar_;
  auto tperp = tperp_;
  auto bmag = bmag_;
  par_for("cgl_lf_precompute", DevExeSpace(), 0, nmb1, 0, ncells3 - 1,
          0, ncells2 - 1, 0, ncells1 - 1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real rho = fmax(w(m,IDN,k,j,i), eos.dfloor);
    tpar(m,k,j,i) = w(m,IPR,k,j,i)/rho;
    tperp(m,k,j,i) = w(m,IPP,k,j,i)/rho;
    bmag(m,k,j,i) = sqrt(SQR(bcc(m,IBX,k,j,i)) + SQR(bcc(m,IBY,k,j,i)) +
                           SQR(bcc(m,IBZ,k,j,i)));
  });

  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  auto size = pmy_pack->pmb->mb_size;
  const Real lf_k = lf_k_parallel;
  const bool local = lf_coeff_local;
  const Real cpar0 = lf_c_parallel0;
  auto &f1 = f.x1f;
  par_for("cgl_lf_flux1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie + 1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real tx = (tpar(m,k,j,i) - tpar(m,k,j,i-1))/size.d_view(m).dx1;
    Real px = (tperp(m,k,j,i) - tperp(m,k,j,i-1))/size.d_view(m).dx1;
    Real bxg = (bmag(m,k,j,i) - bmag(m,k,j,i-1))/size.d_view(m).dx1;
    Real ty = 0.0, py = 0.0, byg = 0.0, tz = 0.0, pz = 0.0, bzg = 0.0;
    if (multi_d) {
      ty = 0.25*(tpar(m,k,j+1,i) - tpar(m,k,j-1,i) +
                 tpar(m,k,j+1,i-1) - tpar(m,k,j-1,i-1))/size.d_view(m).dx2;
      py = 0.25*(tperp(m,k,j+1,i) - tperp(m,k,j-1,i) +
                 tperp(m,k,j+1,i-1) - tperp(m,k,j-1,i-1))/size.d_view(m).dx2;
      byg = 0.25*(bmag(m,k,j+1,i) - bmag(m,k,j-1,i) +
                  bmag(m,k,j+1,i-1) - bmag(m,k,j-1,i-1))/size.d_view(m).dx2;
    }
    if (three_d) {
      tz = 0.25*(tpar(m,k+1,j,i) - tpar(m,k-1,j,i) +
                 tpar(m,k+1,j,i-1) - tpar(m,k-1,j,i-1))/size.d_view(m).dx3;
      pz = 0.25*(tperp(m,k+1,j,i) - tperp(m,k-1,j,i) +
                 tperp(m,k+1,j,i-1) - tperp(m,k-1,j,i-1))/size.d_view(m).dx3;
      bzg = 0.25*(bmag(m,k+1,j,i) - bmag(m,k-1,j,i) +
                  bmag(m,k+1,j,i-1) - bmag(m,k-1,j,i-1))/size.d_view(m).dx3;
    }
    const Real bx = 0.5*(bcc(m,IBX,k,j,i-1) + bcc(m,IBX,k,j,i));
    const Real by = 0.5*(bcc(m,IBY,k,j,i-1) + bcc(m,IBY,k,j,i));
    const Real bz = 0.5*(bcc(m,IBZ,k,j,i-1) + bcc(m,IBZ,k,j,i));
    CGLLFFaceState face;
    Real eflux = 0.0, muflux = 0.0;
    if (BuildCGLLFFaceState(w(m,IDN,k,j,i-1), w(m,IDN,k,j,i),
                            w(m,IPR,k,j,i-1), w(m,IPR,k,j,i),
                            w(m,IPP,k,j,i-1), w(m,IPP,k,j,i),
                            bx, by, bz, 0, lf_k, local, cpar0, eos, face)) {
      CGLLFFlux(face, tx, ty, tz, px, py, pz, bxg, byg, bzg, eflux, muflux);
    }
    f1(m,IEN,k,j,i) = eflux;
    f1(m,IAN,k,j,i) = muflux;
  });
  if (pmy_pack->pmesh->one_d) return;

  auto &f2 = f.x2f;
  par_for("cgl_lf_flux2", DevExeSpace(), 0, nmb1, ks, ke, js, je + 1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real tx = 0.25*(tpar(m,k,j,i+1) - tpar(m,k,j,i-1) +
                          tpar(m,k,j-1,i+1) - tpar(m,k,j-1,i-1))/size.d_view(m).dx1;
    const Real px = 0.25*(tperp(m,k,j,i+1) - tperp(m,k,j,i-1) +
                          tperp(m,k,j-1,i+1) - tperp(m,k,j-1,i-1))/size.d_view(m).dx1;
    const Real bxg = 0.25*(bmag(m,k,j,i+1) - bmag(m,k,j,i-1) +
                           bmag(m,k,j-1,i+1) - bmag(m,k,j-1,i-1))/size.d_view(m).dx1;
    const Real ty = (tpar(m,k,j,i) - tpar(m,k,j-1,i))/size.d_view(m).dx2;
    const Real py = (tperp(m,k,j,i) - tperp(m,k,j-1,i))/size.d_view(m).dx2;
    const Real byg = (bmag(m,k,j,i) - bmag(m,k,j-1,i))/size.d_view(m).dx2;
    Real tz = 0.0, pz = 0.0, bzg = 0.0;
    if (three_d) {
      tz = 0.25*(tpar(m,k+1,j,i) - tpar(m,k-1,j,i) +
                 tpar(m,k+1,j-1,i) - tpar(m,k-1,j-1,i))/size.d_view(m).dx3;
      pz = 0.25*(tperp(m,k+1,j,i) - tperp(m,k-1,j,i) +
                 tperp(m,k+1,j-1,i) - tperp(m,k-1,j-1,i))/size.d_view(m).dx3;
      bzg = 0.25*(bmag(m,k+1,j,i) - bmag(m,k-1,j,i) +
                  bmag(m,k+1,j-1,i) - bmag(m,k-1,j-1,i))/size.d_view(m).dx3;
    }
    const Real bx = 0.5*(bcc(m,IBX,k,j-1,i) + bcc(m,IBX,k,j,i));
    const Real by = 0.5*(bcc(m,IBY,k,j-1,i) + bcc(m,IBY,k,j,i));
    const Real bz = 0.5*(bcc(m,IBZ,k,j-1,i) + bcc(m,IBZ,k,j,i));
    CGLLFFaceState face;
    Real eflux = 0.0, muflux = 0.0;
    if (BuildCGLLFFaceState(w(m,IDN,k,j-1,i), w(m,IDN,k,j,i),
                            w(m,IPR,k,j-1,i), w(m,IPR,k,j,i),
                            w(m,IPP,k,j-1,i), w(m,IPP,k,j,i),
                            bx, by, bz, 1, lf_k, local, cpar0, eos, face)) {
      CGLLFFlux(face, tx, ty, tz, px, py, pz, bxg, byg, bzg, eflux, muflux);
    }
    f2(m,IEN,k,j,i) = eflux;
    f2(m,IAN,k,j,i) = muflux;
  });
  if (pmy_pack->pmesh->two_d) return;

  auto &f3 = f.x3f;
  par_for("cgl_lf_flux3", DevExeSpace(), 0, nmb1, ks, ke + 1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real tx = 0.25*(tpar(m,k,j,i+1) - tpar(m,k,j,i-1) +
                          tpar(m,k-1,j,i+1) - tpar(m,k-1,j,i-1))/size.d_view(m).dx1;
    const Real px = 0.25*(tperp(m,k,j,i+1) - tperp(m,k,j,i-1) +
                          tperp(m,k-1,j,i+1) - tperp(m,k-1,j,i-1))/size.d_view(m).dx1;
    const Real bxg = 0.25*(bmag(m,k,j,i+1) - bmag(m,k,j,i-1) +
                           bmag(m,k-1,j,i+1) - bmag(m,k-1,j,i-1))/size.d_view(m).dx1;
    const Real ty = 0.25*(tpar(m,k,j+1,i) - tpar(m,k,j-1,i) +
                          tpar(m,k-1,j+1,i) - tpar(m,k-1,j-1,i))/size.d_view(m).dx2;
    const Real py = 0.25*(tperp(m,k,j+1,i) - tperp(m,k,j-1,i) +
                          tperp(m,k-1,j+1,i) - tperp(m,k-1,j-1,i))/size.d_view(m).dx2;
    const Real byg = 0.25*(bmag(m,k,j+1,i) - bmag(m,k,j-1,i) +
                           bmag(m,k-1,j+1,i) - bmag(m,k-1,j-1,i))/size.d_view(m).dx2;
    const Real tz = (tpar(m,k,j,i) - tpar(m,k-1,j,i))/size.d_view(m).dx3;
    const Real pz = (tperp(m,k,j,i) - tperp(m,k-1,j,i))/size.d_view(m).dx3;
    const Real bzg = (bmag(m,k,j,i) - bmag(m,k-1,j,i))/size.d_view(m).dx3;
    const Real bx = 0.5*(bcc(m,IBX,k-1,j,i) + bcc(m,IBX,k,j,i));
    const Real by = 0.5*(bcc(m,IBY,k-1,j,i) + bcc(m,IBY,k,j,i));
    const Real bz = 0.5*(bcc(m,IBZ,k-1,j,i) + bcc(m,IBZ,k,j,i));
    CGLLFFaceState face;
    Real eflux = 0.0, muflux = 0.0;
    if (BuildCGLLFFaceState(w(m,IDN,k-1,j,i), w(m,IDN,k,j,i),
                            w(m,IPR,k-1,j,i), w(m,IPR,k,j,i),
                            w(m,IPP,k-1,j,i), w(m,IPP,k,j,i),
                            bx, by, bz, 2, lf_k, local, cpar0, eos, face)) {
      CGLLFFlux(face, tx, ty, tz, px, py, pz, bxg, byg, bzg, eflux, muflux);
    }
    f3(m,IEN,k,j,i) = eflux;
    f3(m,IAN,k,j,i) = muflux;
  });
}

void CGLLandauFluid::NewTimeStep(const DvceArray5D<Real> &w, const EOS_Data &eos) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, nx1 = indcs.nx1;
  const int js = indcs.js, nx2 = indcs.nx2;
  const int ks = indcs.ks, nx3 = indcs.nx3;
  const int nmkji = pmy_pack->nmb_thispack*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji = nx2*nx1;
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  const Real kpar = lf_k_parallel;
  const bool local = lf_coeff_local;
  const Real cpar0 = lf_c_parallel0;
  auto size = pmy_pack->pmb->mb_size;
  dtnew = static_cast<Real>(std::numeric_limits<float>::max());
  Kokkos::parallel_reduce("cgl_lf_newdt", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, Real &min_dt) {
    const int m = idx/nkji;
    int k = (idx - m*nkji)/nji + ks;
    int j = (idx - m*nkji - (k - ks)*nji)/nx1 + js;
    const int i = idx - m*nkji - (k - ks)*nji - (j - js)*nx1 + is;
    const Real rho = fmax(w(m,IDN,k,j,i), eos.dfloor);
    const Real cpar = local ? sqrt(fmax(w(m,IPR,k,j,i)/rho, eos.tfloor)) : cpar0;
    const Real chi = static_cast<Real>(1.5957691216057308)*cpar/kpar;
    if (chi > 0.0) {
      min_dt = fmin(min_dt, SQR(size.d_view(m).dx1)/chi);
      if (multi_d) min_dt = fmin(min_dt, SQR(size.d_view(m).dx2)/chi);
      if (three_d) min_dt = fmin(min_dt, SQR(size.d_view(m).dx3)/chi);
    }
  }, Kokkos::Min<Real>(dtnew));
  const Real fac = three_d ? static_cast<Real>(1.0/6.0) :
                   (multi_d ? static_cast<Real>(0.25) : static_cast<Real>(0.5));
  dtnew *= fac;
}
