//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_landau_fluid.cpp
//! \brief CGL Landau-fluid heat-flux closure and parabolic timestep bound.

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

#include "athena.hpp"
#include "diffusion/limiters.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "mesh/nghbr_index.hpp"
#include "eos/eos.hpp"
#include "eos/cgl_physics.hpp"
#include "eos/ideal_c2p_mhd.hpp"
#include "diffusion/cgl_landau_fluid.hpp"
#include "diffusion/cgl_landau_fluid_arithmetic.hpp"

namespace {

bool CGLProfileEnvValue(const char *name, bool fallback) {
  const char *value = std::getenv(name);
  if (value == nullptr) {
    return fallback;
  }
  const std::string text(value);
  if (text == "1" || text == "true" || text == "TRUE" ||
      text == "yes" || text == "YES" || text == "on" || text == "ON") {
    return true;
  }
  if (text == "0" || text == "false" || text == "FALSE" ||
      text == "no" || text == "NO" || text == "off" || text == "OFF") {
    return false;
  }
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
            << name << " = '" << text
            << "' is not a boolean value." << std::endl;
  std::exit(EXIT_FAILURE);
}

std::string CGLEnvStringValue(const char *name, const std::string &fallback) {
  const char *value = std::getenv(name);
  if (value == nullptr) {
    return fallback;
  }
  return std::string(value);
}

CGLLFDiagnosticsMode ParseCGLLFDiagnosticsMode(const std::string &mode) {
  if (mode == "full") {
    return CGLLFDiagnosticsMode::full;
  }
  if (mode == "none") {
    return CGLLFDiagnosticsMode::none;
  }
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
            << "<mhd>/cgl_lf_diagnostics = '" << mode
            << "' is not implemented; valid choices are [full,none]."
            << std::endl;
  std::exit(EXIT_FAILURE);
}

CGLLFArithmeticMode ParseCGLLFArithmeticMode(const std::string &mode) {
  if (mode == "safe") {
    return CGLLFArithmeticMode::safe;
  }
  if (mode == "fast") {
    return CGLLFArithmeticMode::fast;
  }
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
            << "<mhd>/cgl_lf_arithmetic = '" << mode
            << "' is not implemented; valid choices are [safe,fast]."
            << std::endl;
  std::exit(EXIT_FAILURE);
}

CGLLFSTSFluxMode ParseCGLLFSTSFluxMode(const std::string &mode) {
  if (mode == "weighted") {
    return CGLLFSTSFluxMode::weighted;
  }
  if (mode == "physical") {
    return CGLLFSTSFluxMode::physical;
  }
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
            << "<mhd>/cgl_lf_sts_flux = '" << mode
            << "' is not implemented; valid choices are [weighted,physical]."
            << std::endl;
  std::exit(EXIT_FAILURE);
}

const char *CGLLFProfileBucketName(CGLLFProfileBucket bucket) {
  switch (bucket) {
    case CGLLFProfileBucket::heat_flux_total:
      return "heat_flux_total";
    case CGLLFProfileBucket::heat_flux_precompute:
      return "heat_flux_precompute";
    case CGLLFProfileBucket::heat_flux_flux1:
      return "heat_flux_flux1";
    case CGLLFProfileBucket::heat_flux_flux1_gradients:
      return "heat_flux_flux1_gradients";
    case CGLLFProfileBucket::heat_flux_flux1_face_state:
      return "heat_flux_flux1_face_state";
    case CGLLFProfileBucket::heat_flux_flux1_closure:
      return "heat_flux_flux1_closure";
    case CGLLFProfileBucket::heat_flux_flux2:
      return "heat_flux_flux2";
    case CGLLFProfileBucket::heat_flux_flux2_gradients:
      return "heat_flux_flux2_gradients";
    case CGLLFProfileBucket::heat_flux_flux2_face_state:
      return "heat_flux_flux2_face_state";
    case CGLLFProfileBucket::heat_flux_flux2_closure:
      return "heat_flux_flux2_closure";
    case CGLLFProfileBucket::heat_flux_flux3:
      return "heat_flux_flux3";
    case CGLLFProfileBucket::heat_flux_flux3_gradients:
      return "heat_flux_flux3_gradients";
    case CGLLFProfileBucket::heat_flux_flux3_face_state:
      return "heat_flux_flux3_face_state";
    case CGLLFProfileBucket::heat_flux_flux3_closure:
      return "heat_flux_flux3_closure";
    case CGLLFProfileBucket::heat_flux_work_diagnostics:
      return "heat_flux_work_diagnostics";
    case CGLLFProfileBucket::timestep_reduction:
      return "timestep_reduction";
    case CGLLFProfileBucket::sweep_begin_conversion:
      return "sweep_begin_conversion";
    case CGLLFProfileBucket::sts_clear_flux:
      return "sts_clear_flux";
    case CGLLFProfileBucket::sts_update_copies:
      return "sts_update_copies";
    case CGLLFProfileBucket::sts_update_kernel:
      return "sts_update_kernel";
    case CGLLFProfileBucket::primitive_refresh:
      return "primitive_refresh";
    case CGLLFProfileBucket::admissibility:
      return "admissibility";
    case CGLLFProfileBucket::sweep_end_conversion:
      return "sweep_end_conversion";
    case CGLLFProfileBucket::post_sweep_collisions:
      return "post_sweep_collisions";
    case CGLLFProfileBucket::parabolic_init_recv:
      return "parabolic_init_recv";
    case CGLLFProfileBucket::parabolic_send_flux:
      return "parabolic_send_flux";
    case CGLLFProfileBucket::parabolic_recv_flux:
      return "parabolic_recv_flux";
    case CGLLFProfileBucket::parabolic_restrict_u:
      return "parabolic_restrict_u";
    case CGLLFProfileBucket::parabolic_send_u:
      return "parabolic_send_u";
    case CGLLFProfileBucket::parabolic_recv_u:
      return "parabolic_recv_u";
    case CGLLFProfileBucket::parabolic_physical_bcs:
      return "parabolic_physical_bcs";
    case CGLLFProfileBucket::parabolic_prolongate:
      return "parabolic_prolongate";
    case CGLLFProfileBucket::count:
      break;
  }
  return "unknown";
}

parabolic::ParabolicIntegratorMode ParseCGLHeatFluxIntegrator(ParameterInput *pin) {
  const std::string integrator =
      pin->GetOrAddString("mhd", "cgl_heat_flux_integrator", "sts");
  if (integrator == "sts") {
    return parabolic::ParabolicIntegratorMode::sts;
  }
  if (integrator == "explicit") {
    return parabolic::ParabolicIntegratorMode::explicit_mode;
  }
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
            << "<mhd>/cgl_heat_flux_integrator = '" << integrator
            << "' is not implemented; valid choices are [sts,explicit]."
            << std::endl;
  std::exit(EXIT_FAILURE);
}

struct CGLLFFaceState {
  Real rho;
  Real ppar;
  Real pperp;
  Real bmag_inv;
  Real bhx;
  Real bhy;
  Real bhz;
  Real bhdir;
  Real cparallel;
  Real lf_k;
  Real nu;
};

KOKKOS_INLINE_FUNCTION
Real ScaledMagneticMagnitude(const Real bx, const Real by, const Real bz) {
  const Real scale = fmax(fabs(bx), fmax(fabs(by), fabs(bz)));
  if (scale == 0.0) {
    return 0.0;
  }
  const Real sx = bx/scale;
  const Real sy = by/scale;
  const Real sz = bz/scale;
  const Real scaled_magnitude = sqrt(sx*sx + sy*sy + sz*sz);
  const Real maximum = std::numeric_limits<Real>::max();
  return (scaled_magnitude > maximum/scale)
             ? maximum
             : scale*scaled_magnitude;
}

// Coefficient one-norm of the transverse VL4 gradient derivative, before
// division by the transverse cell width. Closure coefficients are frozen.
KOKKOS_INLINE_FUNCTION
Real CGLLFVL4DerivativeNorm(const Real a, const Real b, const Real c, const Real d) {
  if (!Kokkos::isfinite(a) || !Kokkos::isfinite(b) ||
      !Kokkos::isfinite(c) || !Kokkos::isfinite(d)) {
    return std::numeric_limits<Real>::infinity();
  }
  const bool nonnegative = a >= 0.0 && b >= 0.0 && c >= 0.0 && d >= 0.0;
  const bool nonpositive = a <= 0.0 && b <= 0.0 && c <= 0.0 && d <= 0.0;
  if (!nonnegative && !nonpositive) return 0.0;
  // VL4 has no unique derivative at zero slopes; cover every adjacent branch.
  const Real minimum = fmin(fmin(fabs(a),fabs(b)),fmin(fabs(c),fabs(d)));
  if (minimum == 0.0) return 8.0;
  // h_r = dVL4/da_r = 4 (w_r/sum w)^2, w_r=min|a|/|a_r|.
  // Scaling avoids squared slopes and overflowing harmonic-mean intermediates.
  const Real wa=minimum/fabs(a), wb=minimum/fabs(b);
  const Real wc=minimum/fabs(c), wd=minimum/fabs(d);
  const Real sum=wa+wb+wc+wd;
  return 8.0*(SQR(fmax(wa,wb)/sum)+SQR(fmax(wc,wd)/sum));
}

// Positive sums represented in log2 space. Used only when direct face-row
// arithmetic loses a finite contribution through overflow or underflow.
KOKKOS_INLINE_FUNCTION
Real CGLLFLogAbs(const Real value) {
  return (value == 0.0) ? -std::numeric_limits<Real>::infinity() :
                          Kokkos::log2(fabs(value));
}

KOKKOS_INLINE_FUNCTION
Real CGLLFLogAdd(const Real a, const Real b) {
  if (Kokkos::isnan(a) || Kokkos::isnan(b)) {
    return std::numeric_limits<Real>::quiet_NaN();
  }
  const Real larger = fmax(a,b), smaller = fmin(a,b);
  if (Kokkos::isinf(larger)) return larger;
  return larger + Kokkos::log2(1.0+Kokkos::exp2(smaller-larger));
}

KOKKOS_INLINE_FUNCTION
Real CGLLFLogScale(const Real logarithm, const Real factor) {
  return (factor == 0.0) ? -std::numeric_limits<Real>::infinity() :
                           logarithm + CGLLFLogAbs(factor);
}

KOKKOS_INLINE_FUNCTION
Real CGLLFLogDiffusivity(const Real speed, const Real kpar, const Real nu,
                         const Real numerator, const Real collision) {
  const Real log_speed = CGLLFLogAbs(speed);
  return CGLLFLogAbs(numerator) + 2.0*log_speed -
      CGLLFLogAdd(CGLLFLogAbs(kpar)+log_speed,
                  CGLLFLogAbs(collision)+CGLLFLogAbs(nu));
}

KOKKOS_INLINE_FUNCTION
bool BuildCGLLFFaceState(const Real rho_l, const Real rho_r,
                         const Real ppar_l, const Real ppar_r,
                         const Real pperp_l, const Real pperp_r,
                         const Real bx, const Real by, const Real bz, const Real bbar,
                         const int dir,
                         const Real lf_k, const bool coeff_local, const Real cparallel0,
                         const bool backup, const EOS_Data &eos,
                         CGLLFFaceState &face) {
  if (bbar <= eos.bfloor || lf_k <= 0.0) {
    return false;
  }
  face.rho = fmax(static_cast<Real>(0.5)*rho_l +
                  static_cast<Real>(0.5)*rho_r, eos.dfloor);
  face.ppar = fmax(static_cast<Real>(0.5)*ppar_l +
                   static_cast<Real>(0.5)*ppar_r, eos.pfloor);
  face.pperp = fmax(static_cast<Real>(0.5)*pperp_l +
                    static_cast<Real>(0.5)*pperp_r, eos.pfloor);
  // Do not renormalize the averaged vector: conduction must weaken when
  // neighboring fields reverse, while all inverse-field factors use Bbar.
  face.bmag_inv = static_cast<Real>(1.0)/bbar;
  face.bhx = bx/bbar;
  face.bhy = by/bbar;
  face.bhz = bz/bbar;
  face.bhdir = (dir == 0) ? face.bhx : ((dir == 1) ? face.bhy : face.bhz);
  face.cparallel =
      coeff_local ? sqrt(fmax(face.ppar/face.rho, eos.tfloor)) : cparallel0;
  if (coeff_local && (!Kokkos::isfinite(face.cparallel) || face.cparallel == 0.0)) {
    // The pressure/density ratio may overflow or underflow although its square
    // root is representable. Preserve ordinary arithmetic outside this corner.
    face.cparallel = fmax(sqrt(face.ppar)/sqrt(face.rho),
                           sqrt(fmax(eos.tfloor,static_cast<Real>(0.0))));
  }
  const Real maximum = std::numeric_limits<Real>::max();
  const Real sqrt_max = sqrt(maximum);
  const Real bsqr = (bbar <= sqrt_max) ? bbar*bbar : maximum;
  const Real nu = fmax(eos.nu_coll, static_cast<Real>(0.0)) +
      cgl::LimiterCollisionRate(face.ppar, face.pperp, bsqr, eos, backup);
  face.lf_k = lf_k;
  face.nu = nu;
  return true;
}

KOKKOS_INLINE_FUNCTION
void CGLLFFlux(const CGLLFFaceState &face, const Real gtpar_x, const Real gtpar_y,
               const Real gtpar_z, const Real gtperp_x, const Real gtperp_y,
               const Real gtperp_z, const Real gb_x, const Real gb_y, const Real gb_z,
               const Real dt_sweep, const Real rkl_weight, Real &eflux,
               Real &muflux, cgl_lf::ScaledValue &weighted_qpar_flux,
               cgl_lf::ScaledValue &weighted_qperp_flux, Real &qpar_ratio,
               Real &qperp_ratio) {
  const Real grad_tpar =
      face.bhx*gtpar_x + face.bhy*gtpar_y + face.bhz*gtpar_z;
  const Real grad_tperp =
      face.bhx*gtperp_x + face.bhy*gtperp_y + face.bhz*gtperp_z;
  const Real grad_b = face.bhx*gb_x + face.bhy*gb_y + face.bhz*gb_z;
  Real signed_qpar_ratio = 0.0;
  signed_qpar_ratio = cgl::ParallelHeatFluxRatio(
      face.cparallel, face.rho, face.ppar, face.lf_k, face.nu, grad_tpar);
  qpar_ratio = fabs(signed_qpar_ratio);
  Real signed_qperp_ratio = 0.0;
  signed_qperp_ratio = cgl::PerpendicularHeatFluxRatio(
      face.cparallel, face.rho, face.ppar, face.pperp, face.bmag_inv,
      face.lf_k, face.nu, grad_tperp, grad_b);
  qperp_ratio = fabs(signed_qperp_ratio);

  weighted_qpar_flux = cgl_lf::LimitedHeatFlux(
      signed_qpar_ratio, face.cparallel, cgl::kSqrtEightOverPi, face.ppar);
  weighted_qpar_flux = cgl_lf::Multiply(weighted_qpar_flux, face.bhdir);
  weighted_qpar_flux = cgl_lf::Multiply(weighted_qpar_flux, dt_sweep);
  weighted_qpar_flux = cgl_lf::Multiply(weighted_qpar_flux, rkl_weight);

  weighted_qperp_flux = cgl_lf::LimitedHeatFlux(
      signed_qperp_ratio, face.cparallel, cgl::kSqrtTwoOverPi, face.pperp);
  weighted_qperp_flux = cgl_lf::Multiply(weighted_qperp_flux, face.bhdir);
  weighted_qperp_flux = cgl_lf::Multiply(weighted_qperp_flux, dt_sweep);
  weighted_qperp_flux = cgl_lf::Multiply(weighted_qperp_flux, rkl_weight);

  // The common stage weight commutes with AMR restriction and remains attached
  // until the final face values are representable.
  eflux = cgl_lf::Materialize(cgl_lf::Add(
      weighted_qperp_flux,
      cgl_lf::Multiply(weighted_qpar_flux, static_cast<Real>(0.5))));
  muflux = cgl_lf::Materialize(
      cgl_lf::Multiply(weighted_qperp_flux, face.bmag_inv));
}

KOKKOS_INLINE_FUNCTION
void CGLLFFluxFast(const CGLLFFaceState &face, const Real gtpar_x,
                   const Real gtpar_y, const Real gtpar_z,
                   const Real gtperp_x, const Real gtperp_y,
                   const Real gtperp_z, const Real gb_x, const Real gb_y,
                   const Real gb_z, const Real dt_sweep,
                   const Real rkl_weight, Real &eflux, Real &muflux,
                   Real &weighted_qpar_flux, Real &weighted_qperp_flux,
                   Real &qpar_ratio, Real &qperp_ratio) {
  const Real grad_tpar =
      face.bhx*gtpar_x + face.bhy*gtpar_y + face.bhz*gtpar_z;
  const Real grad_tperp =
      face.bhx*gtperp_x + face.bhy*gtperp_y + face.bhz*gtperp_z;
  const Real grad_b = face.bhx*gb_x + face.bhy*gb_y + face.bhz*gb_z;
  const Real signed_qpar_ratio = cgl::ParallelHeatFluxRatio(
      face.cparallel, face.rho, face.ppar, face.lf_k, face.nu, grad_tpar);
  const Real signed_qperp_ratio = cgl::PerpendicularHeatFluxRatio(
      face.cparallel, face.rho, face.ppar, face.pperp, face.bmag_inv,
      face.lf_k, face.nu, grad_tperp, grad_b);
  qpar_ratio = fabs(signed_qpar_ratio);
  qperp_ratio = fabs(signed_qperp_ratio);

  const Real qpar = cgl_lf::LimitedRatio(signed_qpar_ratio)*face.cparallel*
                    cgl::kSqrtEightOverPi*face.ppar;
  const Real qperp = cgl_lf::LimitedRatio(signed_qperp_ratio)*face.cparallel*
                     cgl::kSqrtTwoOverPi*face.pperp;
  const Real common = face.bhdir*dt_sweep*rkl_weight;
  weighted_qpar_flux = qpar*common;
  weighted_qperp_flux = qperp*common;
  eflux = weighted_qperp_flux + static_cast<Real>(0.5)*weighted_qpar_flux;
  muflux = weighted_qperp_flux*face.bmag_inv;
}

KOKKOS_INLINE_FUNCTION
bool OwnsHeatFluxDiagnosticFace(const int m, const int direction, const int face_index,
                                const int lower, const int upper, const bool multilevel,
                                const int my_level,
                                const DualArray2D<NeighborBlock> &nghbr) {
  if (face_index > lower && face_index <= upper) {
    return true;
  }
  if (face_index != lower && face_index != upper + 1) {
    return false;
  }
  if (!multilevel) {
    return face_index == lower;
  }

  const int side = (face_index == lower) ? -1 : 1;
  for (int n1 = 0; n1 < 2; ++n1) {
    for (int n2 = 0; n2 < 2; ++n2) {
      const int neighbor = (direction == 0) ? NeighborIndex(side, 0, 0, n1, n2)
                         : ((direction == 1) ? NeighborIndex(0, side, 0, n1, n2)
                                             : NeighborIndex(0, 0, side, n1, n2));
      if (nghbr.d_view(m, neighbor).gid < 0) {
        continue;
      }
      const int neighbor_level = nghbr.d_view(m, neighbor).lev;
      if (face_index == lower && neighbor_level > my_level) {
        return false;
      }
      if (face_index == upper + 1 && neighbor_level < my_level) {
        return true;
      }
    }
  }
  return face_index == lower;
}

KOKKOS_INLINE_FUNCTION
void AccumulateCGLLFDiagnosticFace(array_sum::GlobalSum &qstats,
                                   const Real qpar_ratio,
                                   const Real qperp_ratio, const Real area,
                                   const Real delta_tpar,
                                   const Real delta_tperp,
                                   const cgl_lf::ScaledValue &weighted_qpar_flux,
                                   const cgl_lf::ScaledValue &weighted_qperp_flux) {
  qstats.the_array[0] += 1.0;
  if (qpar_ratio > 1.0) qstats.the_array[1] += 1.0;
  if (qpar_ratio > 10.0) qstats.the_array[2] += 1.0;
  if (qperp_ratio > 1.0) qstats.the_array[3] += 1.0;
  if (qperp_ratio > 10.0) qstats.the_array[4] += 1.0;
  auto qpar_work = cgl_lf::Multiply(weighted_qpar_flux, -area);
  qpar_work = cgl_lf::Multiply(qpar_work, delta_tpar);
  auto qperp_work = cgl_lf::Multiply(weighted_qperp_flux, -area);
  qperp_work = cgl_lf::Multiply(qperp_work, delta_tperp);
  qstats.the_array[5] += cgl_lf::Materialize(qpar_work);
  qstats.the_array[6] += cgl_lf::Materialize(qperp_work);
}

KOKKOS_INLINE_FUNCTION
void AccumulateCGLLFDiagnosticFace(array_sum::GlobalSum &qstats,
                                   const Real qpar_ratio,
                                   const Real qperp_ratio, const Real area,
                                   const Real delta_tpar,
                                   const Real delta_tperp,
                                   const Real weighted_qpar_flux,
                                   const Real weighted_qperp_flux) {
  qstats.the_array[0] += 1.0;
  if (qpar_ratio > 1.0) qstats.the_array[1] += 1.0;
  if (qpar_ratio > 10.0) qstats.the_array[2] += 1.0;
  if (qperp_ratio > 1.0) qstats.the_array[3] += 1.0;
  if (qperp_ratio > 10.0) qstats.the_array[4] += 1.0;
  qstats.the_array[5] += -area*weighted_qpar_flux*delta_tpar;
  qstats.the_array[6] += -area*weighted_qperp_flux*delta_tperp;
}

} // namespace

CGLLandauFluid::CGLLandauFluid(MeshBlockPack *pp, ParameterInput *pin) :
    dtnew(static_cast<Real>(std::numeric_limits<float>::max())),
    lf_k_parallel(0.0),
    lf_coeff_local(true),
    lf_c_parallel0(0.0),
    strict_admissibility(false),
    effective_backup_limiter(false),
    diagnostics_mode(CGLLFDiagnosticsMode::full),
    arithmetic_mode(CGLLFArithmeticMode::safe),
    sts_flux_mode(CGLLFSTSFluxMode::weighted),
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
  if (!(lf_k_parallel > 0.0) || !std::isfinite(lf_k_parallel)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<mhd>/lf_k_parallel must be finite and positive." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const std::string coeff_mode =
      pin->GetOrAddString("mhd", "lf_coefficient_mode", "local");
  if (coeff_mode == "background") {
    lf_coeff_local = false;
    lf_c_parallel0 = pin->GetReal("mhd", "lf_c_parallel0");
    if (!(lf_c_parallel0 > 0.0) || !std::isfinite(lf_c_parallel0)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<mhd>/lf_c_parallel0 must be finite and positive." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else if (coeff_mode != "local") {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<mhd>/lf_coefficient_mode = '" << coeff_mode
              << "' must be 'local' or 'background'." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  mode = ParseCGLHeatFluxIntegrator(pin);
  diagnostics_mode = ParseCGLLFDiagnosticsMode(
      CGLEnvStringValue(
          "ATHENAK_CGL_LF_DIAGNOSTICS",
          pin->GetOrAddString("mhd", "cgl_lf_diagnostics", "full")));
  arithmetic_mode = ParseCGLLFArithmeticMode(
      CGLEnvStringValue(
          "ATHENAK_CGL_LF_ARITHMETIC",
          pin->GetOrAddString("mhd", "cgl_lf_arithmetic", "safe")));
  sts_flux_mode = ParseCGLLFSTSFluxMode(
      CGLEnvStringValue(
          "ATHENAK_CGL_LF_STS_FLUX",
          pin->GetOrAddString("mhd", "cgl_lf_sts_flux", "weighted")));
  if (sts_flux_mode == CGLLFSTSFluxMode::physical &&
      diagnostics_mode == CGLLFDiagnosticsMode::full) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<mhd>/cgl_lf_sts_flux = 'physical' requires "
              << "<mhd>/cgl_lf_diagnostics = 'none'; q-work diagnostics are "
              << "defined for the weighted STS RHS." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (sts_flux_mode == CGLLFSTSFluxMode::physical &&
      arithmetic_mode != CGLLFArithmeticMode::fast) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<mhd>/cgl_lf_sts_flux = 'physical' requires "
              << "<mhd>/cgl_lf_arithmetic = 'fast'; safe arithmetic keeps the "
              << "weighted STS RHS." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  strict_admissibility =
      pin->GetOrAddBoolean("mhd", "cgl_lf_strict_admissibility", false);
  const bool configured_backup =
      pin->GetOrAddBoolean("mhd", "backup_limiters", false);
  const bool instability_limiter_active =
      pin->GetOrAddBoolean("mhd", "mirror_limiter", false) ||
      pin->GetOrAddBoolean("mhd", "firehose_limiter", false);
  effective_backup_limiter = cgl::EffectiveBackupLimiter(
      configured_backup, true, instability_limiter_active,
      strict_admissibility);

  profile_enabled_ =
      CGLProfileEnvValue("ATHENAK_CGL_LF_PROFILE",
                         pin->GetOrAddBoolean("mhd", "cgl_lf_profile", false));
  profile_detail_enabled_ =
      CGLProfileEnvValue(
          "ATHENAK_CGL_LF_PROFILE_DETAIL",
          pin->GetOrAddBoolean("mhd", "cgl_lf_profile_detail", false));
  profile_detail_enabled_ = profile_enabled_ && profile_detail_enabled_;
  if (profile_enabled_ && global_variable::my_rank == 0) {
    std::cout << "CGL Landau-fluid profiling enabled; timing regions use Kokkos "
              << "fences and are intended for profiling runs only." << std::endl;
    if (profile_detail_enabled_) {
      std::cout << "CGL Landau-fluid detailed profiling enabled; additional "
                << "probe kernels replay directional heat-flux sub-work and "
                << "do not update evolved state." << std::endl;
    }
  }
}

CGLLFProfileRegion::CGLLFProfileRegion(CGLLandauFluid *profile,
                                       CGLLFProfileBucket bucket) :
    profile_(profile),
    bucket_(bucket),
    active_(profile != nullptr && profile->ProfileEnabled()) {
  if (active_) {
    Kokkos::fence();
    timer_.reset();
  }
}

CGLLFProfileRegion::~CGLLFProfileRegion() {
  if (active_) {
    Kokkos::fence();
    profile_->AddProfileTime(bucket_, static_cast<Real>(timer_.seconds()));
  }
}

void CGLLandauFluid::AddProfileTime(CGLLFProfileBucket bucket, Real seconds) {
  if (!profile_enabled_) {
    return;
  }
  const int index = static_cast<int>(bucket);
  if (index < 0 || index >= kCGLLFProfileBucketCount) {
    return;
  }
  profile_seconds_[index] += seconds;
  profile_counts_[index] += 1;
}

void CGLLandauFluid::ReportProfile(const char *context) const {
  if (!profile_enabled_) {
    return;
  }

  Real local_seconds[kCGLLFProfileBucketCount] = {};
  Real local_counts[kCGLLFProfileBucketCount] = {};
  Real sum_seconds[kCGLLFProfileBucketCount] = {};
  Real max_seconds[kCGLLFProfileBucketCount] = {};
  Real sum_counts[kCGLLFProfileBucketCount] = {};
  Real max_counts[kCGLLFProfileBucketCount] = {};
  for (int n = 0; n < kCGLLFProfileBucketCount; ++n) {
    local_seconds[n] = profile_seconds_[n];
    local_counts[n] = static_cast<Real>(profile_counts_[n]);
  }

  int max_nstages = profile_max_nstages_;
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(local_seconds, sum_seconds, kCGLLFProfileBucketCount, MPI_ATHENA_REAL,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(local_seconds, max_seconds, kCGLLFProfileBucketCount, MPI_ATHENA_REAL,
                MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(local_counts, sum_counts, kCGLLFProfileBucketCount, MPI_ATHENA_REAL,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(local_counts, max_counts, kCGLLFProfileBucketCount, MPI_ATHENA_REAL,
                MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&profile_max_nstages_, &max_nstages, 1, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);
#else
  for (int n = 0; n < kCGLLFProfileBucketCount; ++n) {
    sum_seconds[n] = local_seconds[n];
    max_seconds[n] = local_seconds[n];
    sum_counts[n] = local_counts[n];
    max_counts[n] = local_counts[n];
  }
#endif

  if (global_variable::my_rank != 0) {
    return;
  }

  const Mesh *pm = pmy_pack->pmesh;
  const auto &mesh_indcs = pm->mesh_indcs;
  const auto &mb_indcs = pm->mb_indcs;
  const char *slurm_nodes = std::getenv("SLURM_JOB_NUM_NODES");
  const Real inv_nranks =
      (global_variable::nranks > 0) ? static_cast<Real>(1.0/global_variable::nranks)
                                    : static_cast<Real>(1.0);

  std::cout << std::endl
            << "CGL Landau-fluid profile summary ("
            << (context != nullptr ? context : "final") << ")" << std::endl;
  std::cout << "  ranks=" << global_variable::nranks;
  if (slurm_nodes != nullptr) {
    std::cout << " slurm_nodes=" << slurm_nodes;
  }
  std::cout << " meshblocks_total=" << pm->nmb_total
            << " meshblocks_rank0=" << pmy_pack->nmb_thispack
            << " meshblock_cells=" << mb_indcs.nx1 << "x" << mb_indcs.nx2
            << "x" << mb_indcs.nx3
            << " global_cells=" << mesh_indcs.nx1 << "x" << mesh_indcs.nx2
            << "x" << mesh_indcs.nx3 << std::endl;
  std::cout << "  lf_k_parallel=" << lf_k_parallel
            << " coeff_mode=" << (lf_coeff_local ? "local" : "background")
            << " strict_admissibility=" << (strict_admissibility ? "true" : "false")
            << " backup_limiter=" << (effective_backup_limiter ? "true" : "false")
            << " diagnostics="
            << (diagnostics_mode == CGLLFDiagnosticsMode::full ? "full" : "none")
            << " arithmetic="
            << (arithmetic_mode == CGLLFArithmeticMode::safe ? "safe" : "fast")
            << " sts_flux="
            << (sts_flux_mode == CGLLFSTSFluxMode::weighted ? "weighted" : "physical")
            << " profile_detail=" << (profile_detail_enabled_ ? "true" : "false")
            << " last_nstages=" << profile_last_nstages_
            << " max_nstages=" << max_nstages << std::endl;
  std::cout << "  timing: rank_mean_s is averaged over ranks; rank_max_s is the "
            << "slowest rank. heat_flux_total includes the three directional flux "
            << "regions." << std::endl;
  if (profile_detail_enabled_) {
    std::cout << "  detailed heat-flux buckets are profile-only probe kernels: "
              << "gradients, face_state, and closure replay sub-work without "
              << "updating evolved flux arrays." << std::endl;
  }
  std::cout << std::left << std::setw(34) << "bucket"
            << std::right << std::setw(16) << "rank_mean_s"
            << std::setw(16) << "rank_max_s"
            << std::setw(18) << "calls_rank_mean"
            << std::setw(16) << "calls_rank_max"
            << std::setw(18) << "max_s_per_call" << std::endl;
  for (int n = 0; n < kCGLLFProfileBucketCount; ++n) {
    if (sum_counts[n] <= 0.0) {
      continue;
    }
    const Real calls_mean = sum_counts[n]*inv_nranks;
    const Real seconds_mean = sum_seconds[n]*inv_nranks;
    const Real seconds_per_call =
        (max_counts[n] > 0.0) ? max_seconds[n]/max_counts[n] : 0.0;
    std::cout << std::left << std::setw(34)
              << CGLLFProfileBucketName(static_cast<CGLLFProfileBucket>(n))
              << std::right << std::scientific << std::setprecision(6)
              << std::setw(16) << seconds_mean
              << std::setw(16) << max_seconds[n]
              << std::fixed << std::setprecision(1)
              << std::setw(18) << calls_mean
              << std::setw(16) << max_counts[n]
              << std::scientific << std::setprecision(6)
              << std::setw(18) << seconds_per_call << std::endl;
  }
}

void CGLLandauFluid::AccumulateHeatFluxDiagnostics(const array_sum::GlobalSum &stats) {
  diagnostics.qfaces += static_cast<std::uint64_t>(stats.the_array[0]);
  diagnostics.qpar_cap += static_cast<std::uint64_t>(stats.the_array[1]);
  diagnostics.qpar_cap10 += static_cast<std::uint64_t>(stats.the_array[2]);
  diagnostics.qperp_cap += static_cast<std::uint64_t>(stats.the_array[3]);
  diagnostics.qperp_cap10 += static_cast<std::uint64_t>(stats.the_array[4]);
  stage_qpar_work_ += stats.the_array[5];
  stage_qperp_work_ += stats.the_array[6];
}

void CGLLandauFluid::ResetHeatFluxDiagnostics() {
  diagnostics.qfaces = 0;
  diagnostics.qpar_cap = 0;
  diagnostics.qpar_cap10 = 0;
  diagnostics.qperp_cap = 0;
  diagnostics.qperp_cap10 = 0;
  diagnostics.qpar_work = 0.0;
  diagnostics.qperp_work = 0.0;
  stage_qpar_work_ = 0.0;
  stage_qperp_work_ = 0.0;
  sweep_qpar_work_ = 0.0;
  sweep_qperp_work_ = 0.0;
  sweep_qpar_work1_ = 0.0;
  sweep_qperp_work1_ = 0.0;
  sweep_qpar_work2_ = 0.0;
  sweep_qperp_work2_ = 0.0;
  sweep_qpar_rhs_ = 0.0;
  sweep_qperp_rhs_ = 0.0;
}

// Uniform LF-only sweeps keep density, momentum, and B fixed between stages.
void CGLLandauFluid::RefreshPrimitives(DvceArray5D<Real> &cons,
                                       const DvceArray5D<Real> &bcc,
                                       DvceArray5D<Real> &prim, const EOS_Data &eos_in,
                                       int il, int iu, int jl, int ju, int kl, int ku) {
  const EOS_Data eos = eos_in;
  const int ni = iu - il + 1;
  const int nji = (ju - jl + 1)*ni;
  const int nkji = (ku - kl + 1)*nji;
  const int nmkji = pmy_pack->nmb_thispack*nkji;
  auto tpar = tpar_, tperp = tperp_, bmag_c2p = bmag_c2p_;
  int nfloord = 0, nfloore = 0, nfloort = 0;
  Kokkos::parallel_reduce("cgl_lf_refresh_and_precompute",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int idx, int &sumd, int &sume, int &sumt) {
    const int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/ni;
    const int i = idx - m*nkji - k*nji - j*ni + il;
    j += jl;
    k += kl;
    MHDCons1D u;
    u.d = cons(m,IDN,k,j,i);
    u.mx = cons(m,IM1,k,j,i);
    u.my = cons(m,IM2,k,j,i);
    u.mz = cons(m,IM3,k,j,i);
    u.e = cons(m,IEN,k,j,i);
    u.mu = cons(m,IAN,k,j,i);
    u.bx = bcc(m,IBX,k,j,i);
    u.by = bcc(m,IBY,k,j,i);
    u.bz = bcc(m,IBZ,k,j,i);
    HydPrim1D w;
    bool dfloor_used=false, efloor_used=false, tfloor_used=false, bfloor_used=false;
    SingleC2P_CGLMHDFromMagneticMoment(u, eos, w, dfloor_used, efloor_used,
                                      tfloor_used, bfloor_used, bmag_c2p(m,k,j,i));
    if (dfloor_used) {
      cons(m,IDN,k,j,i) = u.d;
      prim(m,IDN,k,j,i) = w.d;
      prim(m,IVX,k,j,i) = w.vx;
      prim(m,IVY,k,j,i) = w.vy;
      prim(m,IVZ,k,j,i) = w.vz;
      ++sumd;
    }
    if (efloor_used) {
      cons(m,IEN,k,j,i) = u.e;
      cons(m,IAN,k,j,i) = u.mu;
      ++sume;
    }
    if (bfloor_used) cons(m,IAN,k,j,i) = u.mu;
    prim(m,IPR,k,j,i) = w.e;
    prim(m,IPP,k,j,i) = w.pp;
    const Real rho = fmax(w.d, eos.dfloor);
    tpar(m,k,j,i) = w.e/rho;
    tperp(m,k,j,i) = w.pp/rho;
    (void) sumt;
  }, Kokkos::Sum<int>(nfloord), Kokkos::Sum<int>(nfloore), Kokkos::Sum<int>(nfloort));
  pmy_pack->pmesh->ecounter.neos_dfloor += nfloord;
  pmy_pack->pmesh->ecounter.neos_efloor += nfloore;
  pmy_pack->pmesh->ecounter.neos_tfloor += nfloort;
}

void CGLLandauFluid::AddHeatFluxes(const DvceArray5D<Real> &w,
                                   const DvceArray5D<Real> &bcc,
                                   const DvceFaceFld4D<Real> &b,
                                   const EOS_Data &eos_in, Real dt_sweep,
                                   Real rkl_weight, DvceFaceFld5D<Real> &f) {
  CGLLFProfileRegion heat_flux_profile(this, CGLLFProfileBucket::heat_flux_total);
  const EOS_Data eos = eos_in;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int ncells1 = indcs.nx1 + 2*indcs.ng;
  const int ncells2 = (indcs.nx2 > 1) ? indcs.nx2 + 2*indcs.ng : 1;
  const int ncells3 = (indcs.nx3 > 1) ? indcs.nx3 + 2*indcs.ng : 1;
  stage_qpar_work_ = 0.0;
  stage_qperp_work_ = 0.0;
  if (tpar_.extent(0) != static_cast<std::size_t>(pmy_pack->nmb_thispack) ||
      tpar_.extent(1) != static_cast<std::size_t>(ncells3) ||
      tpar_.extent(2) != static_cast<std::size_t>(ncells2) ||
      tpar_.extent(3) != static_cast<std::size_t>(ncells1)) {
    Kokkos::realloc(tpar_, pmy_pack->nmb_thispack, ncells3, ncells2, ncells1);
    Kokkos::realloc(tperp_, pmy_pack->nmb_thispack, ncells3, ncells2, ncells1);
    Kokkos::realloc(bmag_, pmy_pack->nmb_thispack, ncells3, ncells2, ncells1);
  }

  if (fused_primitive_refresh_ &&
      (bmag_c2p_.extent(0) != tpar_.extent(0) ||
       bmag_c2p_.extent(1) != tpar_.extent(1) ||
       bmag_c2p_.extent(2) != tpar_.extent(2) ||
       bmag_c2p_.extent(3) != tpar_.extent(3))) {
    Kokkos::realloc(bmag_c2p_, pmy_pack->nmb_thispack, ncells3, ncells2, ncells1);
  }
  auto tpar = tpar_;
  auto tperp = tperp_;
  auto bmag = bmag_;
  auto bmag_c2p = bmag_c2p_;
  const bool fused = fused_primitive_refresh_;
  if (!fused || !precomputed_) {
    CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_precompute);
    par_for("cgl_lf_precompute", DevExeSpace(), 0, nmb1, 0, ncells3 - 1,
            0, ncells2 - 1, 0, ncells1 - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real rho = fmax(w(m,IDN,k,j,i), eos.dfloor);
      tpar(m,k,j,i) = w(m,IPR,k,j,i)/rho;
      tperp(m,k,j,i) = w(m,IPP,k,j,i)/rho;
      bmag(m,k,j,i) = ScaledMagneticMagnitude(
          bcc(m,IBX,k,j,i), bcc(m,IBY,k,j,i), bcc(m,IBZ,k,j,i));
      if (fused) {
        const Real bsqr = SQR(bcc(m,IBX,k,j,i)) + SQR(bcc(m,IBY,k,j,i))
                         + SQR(bcc(m,IBZ,k,j,i));
        bmag_c2p(m,k,j,i) = sqrt(bsqr);
      }
    });
    precomputed_ = fused;
  }

  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  const bool multilevel = pmy_pack->pmesh->multilevel;
  auto size = pmy_pack->pmb->mb_size;
  auto nghbr = pmy_pack->pmb->nghbr;
  auto mblev = pmy_pack->pmb->mb_lev;
  const Real lf_k = lf_k_parallel;
  const bool local = lf_coeff_local;
  const Real cpar0 = lf_c_parallel0;
  const bool backup = effective_backup_limiter;
  const bool collect_heat_flux_diagnostics =
      diagnostics_mode == CGLLFDiagnosticsMode::full;
  const bool fast_arithmetic = arithmetic_mode == CGLLFArithmeticMode::fast;
  const bool weighted_sts_flux = sts_flux_mode == CGLLFSTSFluxMode::weighted;
  if (!weighted_sts_flux) {
    dt_sweep = 1.0;
    rkl_weight = 1.0;
  }
  if (!collect_heat_flux_diagnostics) {
    ResetHeatFluxDiagnostics();
  }
  // Ordinary execution owns each cell's lower x/y/z face, including padded
  // high caps. Directional guards precede every stencil read. Detailed profiling
  // retains the original directional kernels and their separate replay buckets.
  if (!profile_detail_enabled_) {
    Kokkos::Profiling::pushRegion("cgl_lf_fluxes_fused");
    auto f1 = f.x1f;
    auto f2 = f.x2f;
    auto f3 = f.x3f;
    const int ni = ie - is + 2;
    const int nj = multi_d ? je - js + 2 : 1;
    const int nk = three_d ? ke - ks + 2 : 1;
    const int nji = nj*ni;
    const int nkji = nk*nji;
    const int nmkji = (nmb1 + 1)*nkji;
    if (collect_heat_flux_diagnostics) {
      array_sum::GlobalSum qstats_fused;
      Kokkos::parallel_reduce("cgl_lf_fluxes_fused",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
      KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &qstats) {
      const int m = idx/nkji;
      const int k = (idx - m*nkji)/nji + ks;
      const int j = (idx - m*nkji - (k - ks)*nji)/ni + js;
      const int i = idx - m*nkji - (k - ks)*nji - (j - js)*ni + is;
      if (j <= je && k <= ke) {
        Real tx = (tpar(m,k,j,i) - tpar(m,k,j,i-1))/size.d_view(m).dx1;
        Real px = (tperp(m,k,j,i) - tperp(m,k,j,i-1))/size.d_view(m).dx1;
        Real bxg = (bmag(m,k,j,i) - bmag(m,k,j,i-1))/size.d_view(m).dx1;
        Real ty = 0.0, py = 0.0, byg = 0.0, tz = 0.0, pz = 0.0, bzg = 0.0;
        if (multi_d) {
          ty = VL4Limiter(tpar(m,k,j+1,i) - tpar(m,k,j,i),
                          tpar(m,k,j,i) - tpar(m,k,j-1,i),
                          tpar(m,k,j+1,i-1) - tpar(m,k,j,i-1),
                          tpar(m,k,j,i-1) - tpar(m,k,j-1,i-1))/size.d_view(m).dx2;
          py = VL4Limiter(tperp(m,k,j+1,i) - tperp(m,k,j,i),
                          tperp(m,k,j,i) - tperp(m,k,j-1,i),
                          tperp(m,k,j+1,i-1) - tperp(m,k,j,i-1),
                          tperp(m,k,j,i-1) - tperp(m,k,j-1,i-1))/size.d_view(m).dx2;
          byg = 0.25*(bmag(m,k,j+1,i) - bmag(m,k,j-1,i) +
                      bmag(m,k,j+1,i-1) - bmag(m,k,j-1,i-1))/size.d_view(m).dx2;
        }
        if (three_d) {
          tz = VL4Limiter(tpar(m,k+1,j,i) - tpar(m,k,j,i),
                          tpar(m,k,j,i) - tpar(m,k-1,j,i),
                          tpar(m,k+1,j,i-1) - tpar(m,k,j,i-1),
                          tpar(m,k,j,i-1) - tpar(m,k-1,j,i-1))/size.d_view(m).dx3;
          pz = VL4Limiter(tperp(m,k+1,j,i) - tperp(m,k,j,i),
                          tperp(m,k,j,i) - tperp(m,k-1,j,i),
                          tperp(m,k+1,j,i-1) - tperp(m,k,j,i-1),
                          tperp(m,k,j,i-1) - tperp(m,k-1,j,i-1))/size.d_view(m).dx3;
          bzg = 0.25*(bmag(m,k+1,j,i) - bmag(m,k-1,j,i) +
                      bmag(m,k+1,j,i-1) - bmag(m,k-1,j,i-1))/size.d_view(m).dx3;
        }
        const Real bx = b.x1f(m,k,j,i);
        const Real by = 0.5*bcc(m,IBY,k,j,i-1) + 0.5*bcc(m,IBY,k,j,i);
        const Real bz = 0.5*bcc(m,IBZ,k,j,i-1) + 0.5*bcc(m,IBZ,k,j,i);
        CGLLFFaceState face;
        Real eflux = 0.0, muflux = 0.0;
        Real qpar_ratio = 0.0, qperp_ratio = 0.0;
        if (BuildCGLLFFaceState(w(m,IDN,k,j,i-1), w(m,IDN,k,j,i),
                                w(m,IPR,k,j,i-1), w(m,IPR,k,j,i),
                                w(m,IPP,k,j,i-1), w(m,IPP,k,j,i),
                                bx, by, bz, 0.5*bmag(m,k,j,i-1) + 0.5*bmag(m,k,j,i),
                                0, lf_k, local, cpar0, backup, eos, face)) {
          if (fast_arithmetic) {
            Real weighted_qpar_flux = 0.0, weighted_qperp_flux = 0.0;
            CGLLFFluxFast(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                          dt_sweep, rkl_weight, eflux, muflux,
                          weighted_qpar_flux, weighted_qperp_flux,
                          qpar_ratio, qperp_ratio);
            if (OwnsHeatFluxDiagnosticFace(m, 0, i, is, ie, multilevel,
                                           mblev.d_view(m), nghbr)) {
              const Real area = size.d_view(m).dx2*size.d_view(m).dx3;
              AccumulateCGLLFDiagnosticFace(qstats, qpar_ratio, qperp_ratio, area,
                                            tpar(m,k,j,i) - tpar(m,k,j,i-1),
                                            tperp(m,k,j,i) - tperp(m,k,j,i-1),
                                            weighted_qpar_flux, weighted_qperp_flux);
            }
          } else {
            cgl_lf::ScaledValue weighted_qpar_flux, weighted_qperp_flux;
            CGLLFFlux(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                      dt_sweep, rkl_weight, eflux, muflux, weighted_qpar_flux,
                      weighted_qperp_flux, qpar_ratio, qperp_ratio);
            if (OwnsHeatFluxDiagnosticFace(m, 0, i, is, ie, multilevel,
                                           mblev.d_view(m), nghbr)) {
              const Real area = size.d_view(m).dx2*size.d_view(m).dx3;
              AccumulateCGLLFDiagnosticFace(qstats, qpar_ratio, qperp_ratio, area,
                                            tpar(m,k,j,i) - tpar(m,k,j,i-1),
                                            tperp(m,k,j,i) - tperp(m,k,j,i-1),
                                            weighted_qpar_flux, weighted_qperp_flux);
            }
          }
        }
        f1(m,IEN,k,j,i) = eflux;
        f1(m,IAN,k,j,i) = muflux;
      }
      if (multi_d && i <= ie && k <= ke) {
        const Real tx = VL4Limiter(tpar(m,k,j,i+1) - tpar(m,k,j,i),
                          tpar(m,k,j,i) - tpar(m,k,j,i-1),
                          tpar(m,k,j-1,i+1) - tpar(m,k,j-1,i),
                          tpar(m,k,j-1,i) - tpar(m,k,j-1,i-1))/size.d_view(m).dx1;
        const Real px = VL4Limiter(tperp(m,k,j,i+1) - tperp(m,k,j,i),
                          tperp(m,k,j,i) - tperp(m,k,j,i-1),
                          tperp(m,k,j-1,i+1) - tperp(m,k,j-1,i),
                          tperp(m,k,j-1,i) - tperp(m,k,j-1,i-1))/size.d_view(m).dx1;
        const Real bxg = 0.25*(bmag(m,k,j,i+1) - bmag(m,k,j,i-1) +
                               bmag(m,k,j-1,i+1) - bmag(m,k,j-1,i-1))/size.d_view(m).dx1;
        const Real ty = (tpar(m,k,j,i) - tpar(m,k,j-1,i))/size.d_view(m).dx2;
        const Real py = (tperp(m,k,j,i) - tperp(m,k,j-1,i))/size.d_view(m).dx2;
        const Real byg = (bmag(m,k,j,i) - bmag(m,k,j-1,i))/size.d_view(m).dx2;
        Real tz = 0.0, pz = 0.0, bzg = 0.0;
        if (three_d) {
          tz = VL4Limiter(tpar(m,k+1,j,i) - tpar(m,k,j,i),
                          tpar(m,k,j,i) - tpar(m,k-1,j,i),
                          tpar(m,k+1,j-1,i) - tpar(m,k,j-1,i),
                          tpar(m,k,j-1,i) - tpar(m,k-1,j-1,i))/size.d_view(m).dx3;
          pz = VL4Limiter(tperp(m,k+1,j,i) - tperp(m,k,j,i),
                          tperp(m,k,j,i) - tperp(m,k-1,j,i),
                          tperp(m,k+1,j-1,i) - tperp(m,k,j-1,i),
                          tperp(m,k,j-1,i) - tperp(m,k-1,j-1,i))/size.d_view(m).dx3;
          bzg = 0.25*(bmag(m,k+1,j,i) - bmag(m,k-1,j,i) +
                      bmag(m,k+1,j-1,i) - bmag(m,k-1,j-1,i))/size.d_view(m).dx3;
        }
        const Real bx = 0.5*bcc(m,IBX,k,j-1,i) + 0.5*bcc(m,IBX,k,j,i);
        const Real by = b.x2f(m,k,j,i);
        const Real bz = 0.5*bcc(m,IBZ,k,j-1,i) + 0.5*bcc(m,IBZ,k,j,i);
        CGLLFFaceState face;
        Real eflux = 0.0, muflux = 0.0;
        Real qpar_ratio = 0.0, qperp_ratio = 0.0;
        if (BuildCGLLFFaceState(w(m,IDN,k,j-1,i), w(m,IDN,k,j,i),
                                w(m,IPR,k,j-1,i), w(m,IPR,k,j,i),
                                w(m,IPP,k,j-1,i), w(m,IPP,k,j,i),
                                bx, by, bz, 0.5*bmag(m,k,j-1,i) + 0.5*bmag(m,k,j,i),
                                1, lf_k, local, cpar0, backup, eos, face)) {
          if (fast_arithmetic) {
            Real weighted_qpar_flux = 0.0, weighted_qperp_flux = 0.0;
            CGLLFFluxFast(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                          dt_sweep, rkl_weight, eflux, muflux,
                          weighted_qpar_flux, weighted_qperp_flux,
                          qpar_ratio, qperp_ratio);
            if (OwnsHeatFluxDiagnosticFace(m, 1, j, js, je, multilevel,
                                           mblev.d_view(m), nghbr)) {
              const Real area = size.d_view(m).dx1*size.d_view(m).dx3;
              AccumulateCGLLFDiagnosticFace(qstats, qpar_ratio, qperp_ratio, area,
                                            tpar(m,k,j,i) - tpar(m,k,j-1,i),
                                            tperp(m,k,j,i) - tperp(m,k,j-1,i),
                                            weighted_qpar_flux, weighted_qperp_flux);
            }
          } else {
            cgl_lf::ScaledValue weighted_qpar_flux, weighted_qperp_flux;
            CGLLFFlux(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                      dt_sweep, rkl_weight, eflux, muflux, weighted_qpar_flux,
                      weighted_qperp_flux, qpar_ratio, qperp_ratio);
            if (OwnsHeatFluxDiagnosticFace(m, 1, j, js, je, multilevel,
                                           mblev.d_view(m), nghbr)) {
              const Real area = size.d_view(m).dx1*size.d_view(m).dx3;
              AccumulateCGLLFDiagnosticFace(qstats, qpar_ratio, qperp_ratio, area,
                                            tpar(m,k,j,i) - tpar(m,k,j-1,i),
                                            tperp(m,k,j,i) - tperp(m,k,j-1,i),
                                            weighted_qpar_flux, weighted_qperp_flux);
            }
          }
        }
        f2(m,IEN,k,j,i) = eflux;
        f2(m,IAN,k,j,i) = muflux;
      }
      if (three_d && i <= ie && j <= je) {
        const Real tx = VL4Limiter(tpar(m,k,j,i+1) - tpar(m,k,j,i),
                          tpar(m,k,j,i) - tpar(m,k,j,i-1),
                          tpar(m,k-1,j,i+1) - tpar(m,k-1,j,i),
                          tpar(m,k-1,j,i) - tpar(m,k-1,j,i-1))/size.d_view(m).dx1;
        const Real px = VL4Limiter(tperp(m,k,j,i+1) - tperp(m,k,j,i),
                          tperp(m,k,j,i) - tperp(m,k,j,i-1),
                          tperp(m,k-1,j,i+1) - tperp(m,k-1,j,i),
                          tperp(m,k-1,j,i) - tperp(m,k-1,j,i-1))/size.d_view(m).dx1;
        const Real bxg = 0.25*(bmag(m,k,j,i+1) - bmag(m,k,j,i-1) +
                               bmag(m,k-1,j,i+1) - bmag(m,k-1,j,i-1))/size.d_view(m).dx1;
        const Real ty = VL4Limiter(tpar(m,k,j+1,i) - tpar(m,k,j,i),
                          tpar(m,k,j,i) - tpar(m,k,j-1,i),
                          tpar(m,k-1,j+1,i) - tpar(m,k-1,j,i),
                          tpar(m,k-1,j,i) - tpar(m,k-1,j-1,i))/size.d_view(m).dx2;
        const Real py = VL4Limiter(tperp(m,k,j+1,i) - tperp(m,k,j,i),
                          tperp(m,k,j,i) - tperp(m,k,j-1,i),
                          tperp(m,k-1,j+1,i) - tperp(m,k-1,j,i),
                          tperp(m,k-1,j,i) - tperp(m,k-1,j-1,i))/size.d_view(m).dx2;
        const Real byg = 0.25*(bmag(m,k,j+1,i) - bmag(m,k,j-1,i) +
                               bmag(m,k-1,j+1,i) - bmag(m,k-1,j-1,i))/size.d_view(m).dx2;
        const Real tz = (tpar(m,k,j,i) - tpar(m,k-1,j,i))/size.d_view(m).dx3;
        const Real pz = (tperp(m,k,j,i) - tperp(m,k-1,j,i))/size.d_view(m).dx3;
        const Real bzg = (bmag(m,k,j,i) - bmag(m,k-1,j,i))/size.d_view(m).dx3;
        const Real bx = 0.5*bcc(m,IBX,k-1,j,i) + 0.5*bcc(m,IBX,k,j,i);
        const Real by = 0.5*bcc(m,IBY,k-1,j,i) + 0.5*bcc(m,IBY,k,j,i);
        const Real bz = b.x3f(m,k,j,i);
        CGLLFFaceState face;
        Real eflux = 0.0, muflux = 0.0;
        Real qpar_ratio = 0.0, qperp_ratio = 0.0;
        if (BuildCGLLFFaceState(w(m,IDN,k-1,j,i), w(m,IDN,k,j,i),
                                w(m,IPR,k-1,j,i), w(m,IPR,k,j,i),
                                w(m,IPP,k-1,j,i), w(m,IPP,k,j,i),
                                bx, by, bz, 0.5*bmag(m,k-1,j,i) + 0.5*bmag(m,k,j,i),
                                2, lf_k, local, cpar0, backup, eos, face)) {
          if (fast_arithmetic) {
            Real weighted_qpar_flux = 0.0, weighted_qperp_flux = 0.0;
            CGLLFFluxFast(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                          dt_sweep, rkl_weight, eflux, muflux,
                          weighted_qpar_flux, weighted_qperp_flux,
                          qpar_ratio, qperp_ratio);
            if (OwnsHeatFluxDiagnosticFace(m, 2, k, ks, ke, multilevel,
                                           mblev.d_view(m), nghbr)) {
              const Real area = size.d_view(m).dx1*size.d_view(m).dx2;
              AccumulateCGLLFDiagnosticFace(qstats, qpar_ratio, qperp_ratio, area,
                                            tpar(m,k,j,i) - tpar(m,k-1,j,i),
                                            tperp(m,k,j,i) - tperp(m,k-1,j,i),
                                            weighted_qpar_flux, weighted_qperp_flux);
            }
          } else {
            cgl_lf::ScaledValue weighted_qpar_flux, weighted_qperp_flux;
            CGLLFFlux(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                      dt_sweep, rkl_weight, eflux, muflux, weighted_qpar_flux,
                      weighted_qperp_flux, qpar_ratio, qperp_ratio);
            if (OwnsHeatFluxDiagnosticFace(m, 2, k, ks, ke, multilevel,
                                           mblev.d_view(m), nghbr)) {
              const Real area = size.d_view(m).dx1*size.d_view(m).dx2;
              AccumulateCGLLFDiagnosticFace(qstats, qpar_ratio, qperp_ratio, area,
                                            tpar(m,k,j,i) - tpar(m,k-1,j,i),
                                            tperp(m,k,j,i) - tperp(m,k-1,j,i),
                                            weighted_qpar_flux, weighted_qperp_flux);
            }
          }
        }
        f3(m,IEN,k,j,i) = eflux;
        f3(m,IAN,k,j,i) = muflux;
      }
      }, Kokkos::Sum<array_sum::GlobalSum>(qstats_fused));
      AccumulateHeatFluxDiagnostics(qstats_fused);
    } else {
      Kokkos::parallel_for("cgl_lf_fluxes_fused",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
      KOKKOS_LAMBDA(const int idx) {
      const int m = idx/nkji;
      const int k = (idx - m*nkji)/nji + ks;
      const int j = (idx - m*nkji - (k - ks)*nji)/ni + js;
      const int i = idx - m*nkji - (k - ks)*nji - (j - js)*ni + is;
      if (j <= je && k <= ke) {
        Real tx = (tpar(m,k,j,i) - tpar(m,k,j,i-1))/size.d_view(m).dx1;
        Real px = (tperp(m,k,j,i) - tperp(m,k,j,i-1))/size.d_view(m).dx1;
        Real bxg = (bmag(m,k,j,i) - bmag(m,k,j,i-1))/size.d_view(m).dx1;
        Real ty = 0.0, py = 0.0, byg = 0.0, tz = 0.0, pz = 0.0, bzg = 0.0;
        if (multi_d) {
          ty = VL4Limiter(tpar(m,k,j+1,i) - tpar(m,k,j,i),
                          tpar(m,k,j,i) - tpar(m,k,j-1,i),
                          tpar(m,k,j+1,i-1) - tpar(m,k,j,i-1),
                          tpar(m,k,j,i-1) - tpar(m,k,j-1,i-1))/size.d_view(m).dx2;
          py = VL4Limiter(tperp(m,k,j+1,i) - tperp(m,k,j,i),
                          tperp(m,k,j,i) - tperp(m,k,j-1,i),
                          tperp(m,k,j+1,i-1) - tperp(m,k,j,i-1),
                          tperp(m,k,j,i-1) - tperp(m,k,j-1,i-1))/size.d_view(m).dx2;
          byg = 0.25*(bmag(m,k,j+1,i) - bmag(m,k,j-1,i) +
                      bmag(m,k,j+1,i-1) - bmag(m,k,j-1,i-1))/size.d_view(m).dx2;
        }
        if (three_d) {
          tz = VL4Limiter(tpar(m,k+1,j,i) - tpar(m,k,j,i),
                          tpar(m,k,j,i) - tpar(m,k-1,j,i),
                          tpar(m,k+1,j,i-1) - tpar(m,k,j,i-1),
                          tpar(m,k,j,i-1) - tpar(m,k-1,j,i-1))/size.d_view(m).dx3;
          pz = VL4Limiter(tperp(m,k+1,j,i) - tperp(m,k,j,i),
                          tperp(m,k,j,i) - tperp(m,k-1,j,i),
                          tperp(m,k+1,j,i-1) - tperp(m,k,j,i-1),
                          tperp(m,k,j,i-1) - tperp(m,k-1,j,i-1))/size.d_view(m).dx3;
          bzg = 0.25*(bmag(m,k+1,j,i) - bmag(m,k-1,j,i) +
                      bmag(m,k+1,j,i-1) - bmag(m,k-1,j,i-1))/size.d_view(m).dx3;
        }
        const Real bx = b.x1f(m,k,j,i);
        const Real by = 0.5*bcc(m,IBY,k,j,i-1) + 0.5*bcc(m,IBY,k,j,i);
        const Real bz = 0.5*bcc(m,IBZ,k,j,i-1) + 0.5*bcc(m,IBZ,k,j,i);
        CGLLFFaceState face;
        Real eflux = 0.0, muflux = 0.0;
        Real qpar_ratio = 0.0, qperp_ratio = 0.0;
        if (BuildCGLLFFaceState(w(m,IDN,k,j,i-1), w(m,IDN,k,j,i),
                                w(m,IPR,k,j,i-1), w(m,IPR,k,j,i),
                                w(m,IPP,k,j,i-1), w(m,IPP,k,j,i),
                                bx, by, bz, 0.5*bmag(m,k,j,i-1) + 0.5*bmag(m,k,j,i),
                                0, lf_k, local, cpar0, backup, eos, face)) {
          if (fast_arithmetic) {
            Real weighted_qpar_flux = 0.0, weighted_qperp_flux = 0.0;
            CGLLFFluxFast(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                          dt_sweep, rkl_weight, eflux, muflux,
                          weighted_qpar_flux, weighted_qperp_flux,
                          qpar_ratio, qperp_ratio);
          } else {
            cgl_lf::ScaledValue weighted_qpar_flux, weighted_qperp_flux;
            CGLLFFlux(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                      dt_sweep, rkl_weight, eflux, muflux, weighted_qpar_flux,
                      weighted_qperp_flux, qpar_ratio, qperp_ratio);
          }
        }
        f1(m,IEN,k,j,i) = eflux;
        f1(m,IAN,k,j,i) = muflux;
      }
      if (multi_d && i <= ie && k <= ke) {
        const Real tx = VL4Limiter(tpar(m,k,j,i+1) - tpar(m,k,j,i),
                          tpar(m,k,j,i) - tpar(m,k,j,i-1),
                          tpar(m,k,j-1,i+1) - tpar(m,k,j-1,i),
                          tpar(m,k,j-1,i) - tpar(m,k,j-1,i-1))/size.d_view(m).dx1;
        const Real px = VL4Limiter(tperp(m,k,j,i+1) - tperp(m,k,j,i),
                          tperp(m,k,j,i) - tperp(m,k,j,i-1),
                          tperp(m,k,j-1,i+1) - tperp(m,k,j-1,i),
                          tperp(m,k,j-1,i) - tperp(m,k,j-1,i-1))/size.d_view(m).dx1;
        const Real bxg = 0.25*(bmag(m,k,j,i+1) - bmag(m,k,j,i-1) +
                               bmag(m,k,j-1,i+1) - bmag(m,k,j-1,i-1))/size.d_view(m).dx1;
        const Real ty = (tpar(m,k,j,i) - tpar(m,k,j-1,i))/size.d_view(m).dx2;
        const Real py = (tperp(m,k,j,i) - tperp(m,k,j-1,i))/size.d_view(m).dx2;
        const Real byg = (bmag(m,k,j,i) - bmag(m,k,j-1,i))/size.d_view(m).dx2;
        Real tz = 0.0, pz = 0.0, bzg = 0.0;
        if (three_d) {
          tz = VL4Limiter(tpar(m,k+1,j,i) - tpar(m,k,j,i),
                          tpar(m,k,j,i) - tpar(m,k-1,j,i),
                          tpar(m,k+1,j-1,i) - tpar(m,k,j-1,i),
                          tpar(m,k,j-1,i) - tpar(m,k-1,j-1,i))/size.d_view(m).dx3;
          pz = VL4Limiter(tperp(m,k+1,j,i) - tperp(m,k,j,i),
                          tperp(m,k,j,i) - tperp(m,k-1,j,i),
                          tperp(m,k+1,j-1,i) - tperp(m,k,j-1,i),
                          tperp(m,k,j-1,i) - tperp(m,k-1,j-1,i))/size.d_view(m).dx3;
          bzg = 0.25*(bmag(m,k+1,j,i) - bmag(m,k-1,j,i) +
                      bmag(m,k+1,j-1,i) - bmag(m,k-1,j-1,i))/size.d_view(m).dx3;
        }
        const Real bx = 0.5*bcc(m,IBX,k,j-1,i) + 0.5*bcc(m,IBX,k,j,i);
        const Real by = b.x2f(m,k,j,i);
        const Real bz = 0.5*bcc(m,IBZ,k,j-1,i) + 0.5*bcc(m,IBZ,k,j,i);
        CGLLFFaceState face;
        Real eflux = 0.0, muflux = 0.0;
        Real qpar_ratio = 0.0, qperp_ratio = 0.0;
        if (BuildCGLLFFaceState(w(m,IDN,k,j-1,i), w(m,IDN,k,j,i),
                                w(m,IPR,k,j-1,i), w(m,IPR,k,j,i),
                                w(m,IPP,k,j-1,i), w(m,IPP,k,j,i),
                                bx, by, bz, 0.5*bmag(m,k,j-1,i) + 0.5*bmag(m,k,j,i),
                                1, lf_k, local, cpar0, backup, eos, face)) {
          if (fast_arithmetic) {
            Real weighted_qpar_flux = 0.0, weighted_qperp_flux = 0.0;
            CGLLFFluxFast(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                          dt_sweep, rkl_weight, eflux, muflux,
                          weighted_qpar_flux, weighted_qperp_flux,
                          qpar_ratio, qperp_ratio);
          } else {
            cgl_lf::ScaledValue weighted_qpar_flux, weighted_qperp_flux;
            CGLLFFlux(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                      dt_sweep, rkl_weight, eflux, muflux, weighted_qpar_flux,
                      weighted_qperp_flux, qpar_ratio, qperp_ratio);
          }
        }
        f2(m,IEN,k,j,i) = eflux;
        f2(m,IAN,k,j,i) = muflux;
      }
      if (three_d && i <= ie && j <= je) {
        const Real tx = VL4Limiter(tpar(m,k,j,i+1) - tpar(m,k,j,i),
                          tpar(m,k,j,i) - tpar(m,k,j,i-1),
                          tpar(m,k-1,j,i+1) - tpar(m,k-1,j,i),
                          tpar(m,k-1,j,i) - tpar(m,k-1,j,i-1))/size.d_view(m).dx1;
        const Real px = VL4Limiter(tperp(m,k,j,i+1) - tperp(m,k,j,i),
                          tperp(m,k,j,i) - tperp(m,k,j,i-1),
                          tperp(m,k-1,j,i+1) - tperp(m,k-1,j,i),
                          tperp(m,k-1,j,i) - tperp(m,k-1,j,i-1))/size.d_view(m).dx1;
        const Real bxg = 0.25*(bmag(m,k,j,i+1) - bmag(m,k,j,i-1) +
                               bmag(m,k-1,j,i+1) - bmag(m,k-1,j,i-1))/size.d_view(m).dx1;
        const Real ty = VL4Limiter(tpar(m,k,j+1,i) - tpar(m,k,j,i),
                          tpar(m,k,j,i) - tpar(m,k,j-1,i),
                          tpar(m,k-1,j+1,i) - tpar(m,k-1,j,i),
                          tpar(m,k-1,j,i) - tpar(m,k-1,j-1,i))/size.d_view(m).dx2;
        const Real py = VL4Limiter(tperp(m,k,j+1,i) - tperp(m,k,j,i),
                          tperp(m,k,j,i) - tperp(m,k,j-1,i),
                          tperp(m,k-1,j+1,i) - tperp(m,k-1,j,i),
                          tperp(m,k-1,j,i) - tperp(m,k-1,j-1,i))/size.d_view(m).dx2;
        const Real byg = 0.25*(bmag(m,k,j+1,i) - bmag(m,k,j-1,i) +
                               bmag(m,k-1,j+1,i) - bmag(m,k-1,j-1,i))/size.d_view(m).dx2;
        const Real tz = (tpar(m,k,j,i) - tpar(m,k-1,j,i))/size.d_view(m).dx3;
        const Real pz = (tperp(m,k,j,i) - tperp(m,k-1,j,i))/size.d_view(m).dx3;
        const Real bzg = (bmag(m,k,j,i) - bmag(m,k-1,j,i))/size.d_view(m).dx3;
        const Real bx = 0.5*bcc(m,IBX,k-1,j,i) + 0.5*bcc(m,IBX,k,j,i);
        const Real by = 0.5*bcc(m,IBY,k-1,j,i) + 0.5*bcc(m,IBY,k,j,i);
        const Real bz = b.x3f(m,k,j,i);
        CGLLFFaceState face;
        Real eflux = 0.0, muflux = 0.0;
        Real qpar_ratio = 0.0, qperp_ratio = 0.0;
        if (BuildCGLLFFaceState(w(m,IDN,k-1,j,i), w(m,IDN,k,j,i),
                                w(m,IPR,k-1,j,i), w(m,IPR,k,j,i),
                                w(m,IPP,k-1,j,i), w(m,IPP,k,j,i),
                                bx, by, bz, 0.5*bmag(m,k-1,j,i) + 0.5*bmag(m,k,j,i),
                                2, lf_k, local, cpar0, backup, eos, face)) {
          if (fast_arithmetic) {
            Real weighted_qpar_flux = 0.0, weighted_qperp_flux = 0.0;
            CGLLFFluxFast(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                          dt_sweep, rkl_weight, eflux, muflux,
                          weighted_qpar_flux, weighted_qperp_flux,
                          qpar_ratio, qperp_ratio);
          } else {
            cgl_lf::ScaledValue weighted_qpar_flux, weighted_qperp_flux;
            CGLLFFlux(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                      dt_sweep, rkl_weight, eflux, muflux, weighted_qpar_flux,
                      weighted_qperp_flux, qpar_ratio, qperp_ratio);
          }
        }
        f3(m,IEN,k,j,i) = eflux;
        f3(m,IAN,k,j,i) = muflux;
      }
      });
    }
    Kokkos::Profiling::popRegion();
    return;
  }
  auto f1 = f.x1f;
  const int ni1 = ie - is + 2;
  const int nj1 = je - js + 1;
  const int nk1 = ke - ks + 1;
  const int nji1 = nj1*ni1;
  const int nkji1 = nk1*nji1;
  const int nmkji1 = (nmb1 + 1)*nkji1;
  if (collect_heat_flux_diagnostics) {
  array_sum::GlobalSum qstats1;
  {
    CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux1);
    Kokkos::parallel_reduce("cgl_lf_flux1",
        Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji1),
    KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &qstats) {
    const int m = idx/nkji1;
    const int k = (idx - m*nkji1)/nji1 + ks;
    const int j = (idx - m*nkji1 - (k - ks)*nji1)/ni1 + js;
    const int i = idx - m*nkji1 - (k - ks)*nji1 - (j - js)*ni1 + is;
    Real tx = (tpar(m,k,j,i) - tpar(m,k,j,i-1))/size.d_view(m).dx1;
    Real px = (tperp(m,k,j,i) - tperp(m,k,j,i-1))/size.d_view(m).dx1;
    Real bxg = (bmag(m,k,j,i) - bmag(m,k,j,i-1))/size.d_view(m).dx1;
    Real ty = 0.0, py = 0.0, byg = 0.0, tz = 0.0, pz = 0.0, bzg = 0.0;
    if (multi_d) {
      ty = VL4Limiter(tpar(m,k,j+1,i) - tpar(m,k,j,i),
                      tpar(m,k,j,i) - tpar(m,k,j-1,i),
                      tpar(m,k,j+1,i-1) - tpar(m,k,j,i-1),
                      tpar(m,k,j,i-1) - tpar(m,k,j-1,i-1))/size.d_view(m).dx2;
      py = VL4Limiter(tperp(m,k,j+1,i) - tperp(m,k,j,i),
                      tperp(m,k,j,i) - tperp(m,k,j-1,i),
                      tperp(m,k,j+1,i-1) - tperp(m,k,j,i-1),
                      tperp(m,k,j,i-1) - tperp(m,k,j-1,i-1))/size.d_view(m).dx2;
      byg = 0.25*(bmag(m,k,j+1,i) - bmag(m,k,j-1,i) +
                  bmag(m,k,j+1,i-1) - bmag(m,k,j-1,i-1))/size.d_view(m).dx2;
    }
    if (three_d) {
      tz = VL4Limiter(tpar(m,k+1,j,i) - tpar(m,k,j,i),
                      tpar(m,k,j,i) - tpar(m,k-1,j,i),
                      tpar(m,k+1,j,i-1) - tpar(m,k,j,i-1),
                      tpar(m,k,j,i-1) - tpar(m,k-1,j,i-1))/size.d_view(m).dx3;
      pz = VL4Limiter(tperp(m,k+1,j,i) - tperp(m,k,j,i),
                      tperp(m,k,j,i) - tperp(m,k-1,j,i),
                      tperp(m,k+1,j,i-1) - tperp(m,k,j,i-1),
                      tperp(m,k,j,i-1) - tperp(m,k-1,j,i-1))/size.d_view(m).dx3;
      bzg = 0.25*(bmag(m,k+1,j,i) - bmag(m,k-1,j,i) +
                  bmag(m,k+1,j,i-1) - bmag(m,k-1,j,i-1))/size.d_view(m).dx3;
    }
    const Real bx = b.x1f(m,k,j,i);
    const Real by = 0.5*bcc(m,IBY,k,j,i-1) + 0.5*bcc(m,IBY,k,j,i);
    const Real bz = 0.5*bcc(m,IBZ,k,j,i-1) + 0.5*bcc(m,IBZ,k,j,i);
    CGLLFFaceState face;
    Real eflux = 0.0, muflux = 0.0;
    Real qpar_ratio = 0.0, qperp_ratio = 0.0;
    if (BuildCGLLFFaceState(w(m,IDN,k,j,i-1), w(m,IDN,k,j,i),
                            w(m,IPR,k,j,i-1), w(m,IPR,k,j,i),
                            w(m,IPP,k,j,i-1), w(m,IPP,k,j,i),
                            bx, by, bz, 0.5*bmag(m,k,j,i-1) + 0.5*bmag(m,k,j,i),
                            0, lf_k, local, cpar0, backup, eos, face)) {
      if (fast_arithmetic) {
        Real weighted_qpar_flux = 0.0, weighted_qperp_flux = 0.0;
        CGLLFFluxFast(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                      dt_sweep, rkl_weight, eflux, muflux,
                      weighted_qpar_flux, weighted_qperp_flux,
                      qpar_ratio, qperp_ratio);
        if (OwnsHeatFluxDiagnosticFace(m, 0, i, is, ie, multilevel,
                                       mblev.d_view(m), nghbr)) {
          const Real area = size.d_view(m).dx2*size.d_view(m).dx3;
          AccumulateCGLLFDiagnosticFace(qstats, qpar_ratio, qperp_ratio, area,
                                        tpar(m,k,j,i) - tpar(m,k,j,i-1),
                                        tperp(m,k,j,i) - tperp(m,k,j,i-1),
                                        weighted_qpar_flux, weighted_qperp_flux);
        }
      } else {
        cgl_lf::ScaledValue weighted_qpar_flux, weighted_qperp_flux;
        CGLLFFlux(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                  dt_sweep, rkl_weight, eflux, muflux, weighted_qpar_flux,
                  weighted_qperp_flux, qpar_ratio, qperp_ratio);
        if (OwnsHeatFluxDiagnosticFace(m, 0, i, is, ie, multilevel,
                                       mblev.d_view(m), nghbr)) {
          const Real area = size.d_view(m).dx2*size.d_view(m).dx3;
          AccumulateCGLLFDiagnosticFace(qstats, qpar_ratio, qperp_ratio, area,
                                        tpar(m,k,j,i) - tpar(m,k,j,i-1),
                                        tperp(m,k,j,i) - tperp(m,k,j,i-1),
                                        weighted_qpar_flux, weighted_qperp_flux);
        }
      }
    }
    f1(m,IEN,k,j,i) = eflux;
    f1(m,IAN,k,j,i) = muflux;
    }, Kokkos::Sum<array_sum::GlobalSum>(qstats1));
  }
  AccumulateHeatFluxDiagnostics(qstats1);
  } else {
    CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux1);
    Kokkos::parallel_for("cgl_lf_flux1",
        Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji1),
    KOKKOS_LAMBDA(const int idx) {
    const int m = idx/nkji1;
    const int k = (idx - m*nkji1)/nji1 + ks;
    const int j = (idx - m*nkji1 - (k - ks)*nji1)/ni1 + js;
    const int i = idx - m*nkji1 - (k - ks)*nji1 - (j - js)*ni1 + is;
    Real tx = (tpar(m,k,j,i) - tpar(m,k,j,i-1))/size.d_view(m).dx1;
    Real px = (tperp(m,k,j,i) - tperp(m,k,j,i-1))/size.d_view(m).dx1;
    Real bxg = (bmag(m,k,j,i) - bmag(m,k,j,i-1))/size.d_view(m).dx1;
    Real ty = 0.0, py = 0.0, byg = 0.0, tz = 0.0, pz = 0.0, bzg = 0.0;
    if (multi_d) {
      ty = VL4Limiter(tpar(m,k,j+1,i) - tpar(m,k,j,i),
                      tpar(m,k,j,i) - tpar(m,k,j-1,i),
                      tpar(m,k,j+1,i-1) - tpar(m,k,j,i-1),
                      tpar(m,k,j,i-1) - tpar(m,k,j-1,i-1))/size.d_view(m).dx2;
      py = VL4Limiter(tperp(m,k,j+1,i) - tperp(m,k,j,i),
                      tperp(m,k,j,i) - tperp(m,k,j-1,i),
                      tperp(m,k,j+1,i-1) - tperp(m,k,j,i-1),
                      tperp(m,k,j,i-1) - tperp(m,k,j-1,i-1))/size.d_view(m).dx2;
      byg = 0.25*(bmag(m,k,j+1,i) - bmag(m,k,j-1,i) +
                  bmag(m,k,j+1,i-1) - bmag(m,k,j-1,i-1))/size.d_view(m).dx2;
    }
    if (three_d) {
      tz = VL4Limiter(tpar(m,k+1,j,i) - tpar(m,k,j,i),
                      tpar(m,k,j,i) - tpar(m,k-1,j,i),
                      tpar(m,k+1,j,i-1) - tpar(m,k,j,i-1),
                      tpar(m,k,j,i-1) - tpar(m,k-1,j,i-1))/size.d_view(m).dx3;
      pz = VL4Limiter(tperp(m,k+1,j,i) - tperp(m,k,j,i),
                      tperp(m,k,j,i) - tperp(m,k-1,j,i),
                      tperp(m,k+1,j,i-1) - tperp(m,k,j,i-1),
                      tperp(m,k,j,i-1) - tperp(m,k-1,j,i-1))/size.d_view(m).dx3;
      bzg = 0.25*(bmag(m,k+1,j,i) - bmag(m,k-1,j,i) +
                  bmag(m,k+1,j,i-1) - bmag(m,k-1,j,i-1))/size.d_view(m).dx3;
    }
    const Real bx = b.x1f(m,k,j,i);
    const Real by = 0.5*bcc(m,IBY,k,j,i-1) + 0.5*bcc(m,IBY,k,j,i);
    const Real bz = 0.5*bcc(m,IBZ,k,j,i-1) + 0.5*bcc(m,IBZ,k,j,i);
    CGLLFFaceState face;
    Real eflux = 0.0, muflux = 0.0;
    Real qpar_ratio = 0.0, qperp_ratio = 0.0;
    if (BuildCGLLFFaceState(w(m,IDN,k,j,i-1), w(m,IDN,k,j,i),
                            w(m,IPR,k,j,i-1), w(m,IPR,k,j,i),
                            w(m,IPP,k,j,i-1), w(m,IPP,k,j,i),
                            bx, by, bz, 0.5*bmag(m,k,j,i-1) + 0.5*bmag(m,k,j,i),
                            0, lf_k, local, cpar0, backup, eos, face)) {
      if (fast_arithmetic) {
        Real weighted_qpar_flux = 0.0, weighted_qperp_flux = 0.0;
        CGLLFFluxFast(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                      dt_sweep, rkl_weight, eflux, muflux,
                      weighted_qpar_flux, weighted_qperp_flux,
                      qpar_ratio, qperp_ratio);
      } else {
        cgl_lf::ScaledValue weighted_qpar_flux, weighted_qperp_flux;
        CGLLFFlux(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                  dt_sweep, rkl_weight, eflux, muflux, weighted_qpar_flux,
                  weighted_qperp_flux, qpar_ratio, qperp_ratio);
      }
    }
    f1(m,IEN,k,j,i) = eflux;
    f1(m,IAN,k,j,i) = muflux;
    });
  }
  if (profile_detail_enabled_) {
    Real detail = 0.0;
    {
      CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux1_gradients);
      Kokkos::parallel_reduce("cgl_lf_flux1_profile_gradients",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji1),
      KOKKOS_LAMBDA(const int idx, Real &sum) {
      const int m = idx/nkji1;
      const int k = (idx - m*nkji1)/nji1 + ks;
      const int j = (idx - m*nkji1 - (k - ks)*nji1)/ni1 + js;
      const int i = idx - m*nkji1 - (k - ks)*nji1 - (j - js)*ni1 + is;
      const Real tx = (tpar(m,k,j,i) - tpar(m,k,j,i-1))/size.d_view(m).dx1;
      const Real px = (tperp(m,k,j,i) - tperp(m,k,j,i-1))/size.d_view(m).dx1;
      const Real bxg = (bmag(m,k,j,i) - bmag(m,k,j,i-1))/size.d_view(m).dx1;
      Real ty = 0.0, py = 0.0, byg = 0.0, tz = 0.0, pz = 0.0, bzg = 0.0;
      if (multi_d) {
        ty = VL4Limiter(tpar(m,k,j+1,i) - tpar(m,k,j,i),
                      tpar(m,k,j,i) - tpar(m,k,j-1,i),
                      tpar(m,k,j+1,i-1) - tpar(m,k,j,i-1),
                      tpar(m,k,j,i-1) - tpar(m,k,j-1,i-1))/size.d_view(m).dx2;
        py = VL4Limiter(tperp(m,k,j+1,i) - tperp(m,k,j,i),
                      tperp(m,k,j,i) - tperp(m,k,j-1,i),
                      tperp(m,k,j+1,i-1) - tperp(m,k,j,i-1),
                      tperp(m,k,j,i-1) - tperp(m,k,j-1,i-1))/size.d_view(m).dx2;
        byg = 0.25*(bmag(m,k,j+1,i) - bmag(m,k,j-1,i) +
                    bmag(m,k,j+1,i-1) - bmag(m,k,j-1,i-1))/size.d_view(m).dx2;
      }
      if (three_d) {
        tz = VL4Limiter(tpar(m,k+1,j,i) - tpar(m,k,j,i),
                      tpar(m,k,j,i) - tpar(m,k-1,j,i),
                      tpar(m,k+1,j,i-1) - tpar(m,k,j,i-1),
                      tpar(m,k,j,i-1) - tpar(m,k-1,j,i-1))/size.d_view(m).dx3;
        pz = VL4Limiter(tperp(m,k+1,j,i) - tperp(m,k,j,i),
                      tperp(m,k,j,i) - tperp(m,k-1,j,i),
                      tperp(m,k+1,j,i-1) - tperp(m,k,j,i-1),
                      tperp(m,k,j,i-1) - tperp(m,k-1,j,i-1))/size.d_view(m).dx3;
        bzg = 0.25*(bmag(m,k+1,j,i) - bmag(m,k-1,j,i) +
                    bmag(m,k+1,j,i-1) - bmag(m,k-1,j,i-1))/size.d_view(m).dx3;
      }
      sum += fabs(tx) + fabs(px) + fabs(bxg) + fabs(ty) + fabs(py) + fabs(byg)
             + fabs(tz) + fabs(pz) + fabs(bzg);
      }, Kokkos::Sum<Real>(detail));
    }
    profile_detail_sink_ += detail;
    detail = 0.0;
    {
      CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux1_face_state);
      Kokkos::parallel_reduce("cgl_lf_flux1_profile_face_state",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji1),
      KOKKOS_LAMBDA(const int idx, Real &sum) {
      const int m = idx/nkji1;
      const int k = (idx - m*nkji1)/nji1 + ks;
      const int j = (idx - m*nkji1 - (k - ks)*nji1)/ni1 + js;
      const int i = idx - m*nkji1 - (k - ks)*nji1 - (j - js)*ni1 + is;
      const Real bx = b.x1f(m,k,j,i);
      const Real by = 0.5*bcc(m,IBY,k,j,i-1) + 0.5*bcc(m,IBY,k,j,i);
      const Real bz = 0.5*bcc(m,IBZ,k,j,i-1) + 0.5*bcc(m,IBZ,k,j,i);
      CGLLFFaceState face;
      if (BuildCGLLFFaceState(w(m,IDN,k,j,i-1), w(m,IDN,k,j,i),
                              w(m,IPR,k,j,i-1), w(m,IPR,k,j,i),
                              w(m,IPP,k,j,i-1), w(m,IPP,k,j,i),
                              bx, by, bz, 0.5*bmag(m,k,j,i-1) + 0.5*bmag(m,k,j,i),
                              0, lf_k, local, cpar0, backup, eos, face)) {
        sum += face.cparallel + face.bmag_inv + fabs(face.bhdir) + face.nu;
      }
      }, Kokkos::Sum<Real>(detail));
    }
    profile_detail_sink_ += detail;
    detail = 0.0;
    {
      CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux1_closure);
      Kokkos::parallel_reduce("cgl_lf_flux1_profile_closure",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji1),
      KOKKOS_LAMBDA(const int idx, Real &sum) {
      const int m = idx/nkji1;
      const int k = (idx - m*nkji1)/nji1 + ks;
      const int j = (idx - m*nkji1 - (k - ks)*nji1)/ni1 + js;
      const int i = idx - m*nkji1 - (k - ks)*nji1 - (j - js)*ni1 + is;
      const Real tx = (tpar(m,k,j,i) - tpar(m,k,j,i-1))/size.d_view(m).dx1;
      const Real px = (tperp(m,k,j,i) - tperp(m,k,j,i-1))/size.d_view(m).dx1;
      const Real bxg = (bmag(m,k,j,i) - bmag(m,k,j,i-1))/size.d_view(m).dx1;
      Real ty = 0.0, py = 0.0, byg = 0.0, tz = 0.0, pz = 0.0, bzg = 0.0;
      if (multi_d) {
        ty = VL4Limiter(tpar(m,k,j+1,i) - tpar(m,k,j,i),
                      tpar(m,k,j,i) - tpar(m,k,j-1,i),
                      tpar(m,k,j+1,i-1) - tpar(m,k,j,i-1),
                      tpar(m,k,j,i-1) - tpar(m,k,j-1,i-1))/size.d_view(m).dx2;
        py = VL4Limiter(tperp(m,k,j+1,i) - tperp(m,k,j,i),
                      tperp(m,k,j,i) - tperp(m,k,j-1,i),
                      tperp(m,k,j+1,i-1) - tperp(m,k,j,i-1),
                      tperp(m,k,j,i-1) - tperp(m,k,j-1,i-1))/size.d_view(m).dx2;
        byg = 0.25*(bmag(m,k,j+1,i) - bmag(m,k,j-1,i) +
                    bmag(m,k,j+1,i-1) - bmag(m,k,j-1,i-1))/size.d_view(m).dx2;
      }
      if (three_d) {
        tz = VL4Limiter(tpar(m,k+1,j,i) - tpar(m,k,j,i),
                      tpar(m,k,j,i) - tpar(m,k-1,j,i),
                      tpar(m,k+1,j,i-1) - tpar(m,k,j,i-1),
                      tpar(m,k,j,i-1) - tpar(m,k-1,j,i-1))/size.d_view(m).dx3;
        pz = VL4Limiter(tperp(m,k+1,j,i) - tperp(m,k,j,i),
                      tperp(m,k,j,i) - tperp(m,k-1,j,i),
                      tperp(m,k+1,j,i-1) - tperp(m,k,j,i-1),
                      tperp(m,k,j,i-1) - tperp(m,k-1,j,i-1))/size.d_view(m).dx3;
        bzg = 0.25*(bmag(m,k+1,j,i) - bmag(m,k-1,j,i) +
                    bmag(m,k+1,j,i-1) - bmag(m,k-1,j,i-1))/size.d_view(m).dx3;
      }
      const Real bx = b.x1f(m,k,j,i);
      const Real by = 0.5*bcc(m,IBY,k,j,i-1) + 0.5*bcc(m,IBY,k,j,i);
      const Real bz = 0.5*bcc(m,IBZ,k,j,i-1) + 0.5*bcc(m,IBZ,k,j,i);
      CGLLFFaceState face;
      Real eflux = 0.0, muflux = 0.0;
      Real qpar_ratio = 0.0, qperp_ratio = 0.0;
      if (BuildCGLLFFaceState(w(m,IDN,k,j,i-1), w(m,IDN,k,j,i),
                              w(m,IPR,k,j,i-1), w(m,IPR,k,j,i),
                              w(m,IPP,k,j,i-1), w(m,IPP,k,j,i),
                              bx, by, bz, 0.5*bmag(m,k,j,i-1) + 0.5*bmag(m,k,j,i),
                              0, lf_k, local, cpar0, backup, eos, face)) {
        if (fast_arithmetic) {
          Real weighted_qpar_flux = 0.0, weighted_qperp_flux = 0.0;
          CGLLFFluxFast(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                        dt_sweep, rkl_weight, eflux, muflux,
                        weighted_qpar_flux, weighted_qperp_flux,
                        qpar_ratio, qperp_ratio);
        } else {
          cgl_lf::ScaledValue weighted_qpar_flux, weighted_qperp_flux;
          CGLLFFlux(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                    dt_sweep, rkl_weight, eflux, muflux, weighted_qpar_flux,
                    weighted_qperp_flux, qpar_ratio, qperp_ratio);
        }
        sum += fabs(eflux) + fabs(muflux) + qpar_ratio + qperp_ratio;
      }
      }, Kokkos::Sum<Real>(detail));
    }
    profile_detail_sink_ += detail;
  }
  if (pmy_pack->pmesh->one_d) {
    return;
  }

  auto f2 = f.x2f;
  const int ni2 = ie - is + 1;
  const int nj2 = je - js + 2;
  const int nk2 = ke - ks + 1;
  const int nji2 = nj2*ni2;
  const int nkji2 = nk2*nji2;
  const int nmkji2 = (nmb1 + 1)*nkji2;
  if (collect_heat_flux_diagnostics) {
  array_sum::GlobalSum qstats2;
  {
    CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux2);
    Kokkos::parallel_reduce("cgl_lf_flux2",
        Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji2),
    KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &qstats) {
    const int m = idx/nkji2;
    const int k = (idx - m*nkji2)/nji2 + ks;
    const int j = (idx - m*nkji2 - (k - ks)*nji2)/ni2 + js;
    const int i = idx - m*nkji2 - (k - ks)*nji2 - (j - js)*ni2 + is;
    const Real tx = VL4Limiter(tpar(m,k,j,i+1) - tpar(m,k,j,i),
                      tpar(m,k,j,i) - tpar(m,k,j,i-1),
                      tpar(m,k,j-1,i+1) - tpar(m,k,j-1,i),
                      tpar(m,k,j-1,i) - tpar(m,k,j-1,i-1))/size.d_view(m).dx1;
    const Real px = VL4Limiter(tperp(m,k,j,i+1) - tperp(m,k,j,i),
                      tperp(m,k,j,i) - tperp(m,k,j,i-1),
                      tperp(m,k,j-1,i+1) - tperp(m,k,j-1,i),
                      tperp(m,k,j-1,i) - tperp(m,k,j-1,i-1))/size.d_view(m).dx1;
    const Real bxg = 0.25*(bmag(m,k,j,i+1) - bmag(m,k,j,i-1) +
                           bmag(m,k,j-1,i+1) - bmag(m,k,j-1,i-1))/size.d_view(m).dx1;
    const Real ty = (tpar(m,k,j,i) - tpar(m,k,j-1,i))/size.d_view(m).dx2;
    const Real py = (tperp(m,k,j,i) - tperp(m,k,j-1,i))/size.d_view(m).dx2;
    const Real byg = (bmag(m,k,j,i) - bmag(m,k,j-1,i))/size.d_view(m).dx2;
    Real tz = 0.0, pz = 0.0, bzg = 0.0;
    if (three_d) {
      tz = VL4Limiter(tpar(m,k+1,j,i) - tpar(m,k,j,i),
                      tpar(m,k,j,i) - tpar(m,k-1,j,i),
                      tpar(m,k+1,j-1,i) - tpar(m,k,j-1,i),
                      tpar(m,k,j-1,i) - tpar(m,k-1,j-1,i))/size.d_view(m).dx3;
      pz = VL4Limiter(tperp(m,k+1,j,i) - tperp(m,k,j,i),
                      tperp(m,k,j,i) - tperp(m,k-1,j,i),
                      tperp(m,k+1,j-1,i) - tperp(m,k,j-1,i),
                      tperp(m,k,j-1,i) - tperp(m,k-1,j-1,i))/size.d_view(m).dx3;
      bzg = 0.25*(bmag(m,k+1,j,i) - bmag(m,k-1,j,i) +
                  bmag(m,k+1,j-1,i) - bmag(m,k-1,j-1,i))/size.d_view(m).dx3;
    }
    const Real bx = 0.5*bcc(m,IBX,k,j-1,i) + 0.5*bcc(m,IBX,k,j,i);
    const Real by = b.x2f(m,k,j,i);
    const Real bz = 0.5*bcc(m,IBZ,k,j-1,i) + 0.5*bcc(m,IBZ,k,j,i);
    CGLLFFaceState face;
    Real eflux = 0.0, muflux = 0.0;
    Real qpar_ratio = 0.0, qperp_ratio = 0.0;
    if (BuildCGLLFFaceState(w(m,IDN,k,j-1,i), w(m,IDN,k,j,i),
                            w(m,IPR,k,j-1,i), w(m,IPR,k,j,i),
                            w(m,IPP,k,j-1,i), w(m,IPP,k,j,i),
                            bx, by, bz, 0.5*bmag(m,k,j-1,i) + 0.5*bmag(m,k,j,i),
                            1, lf_k, local, cpar0, backup, eos, face)) {
      if (fast_arithmetic) {
        Real weighted_qpar_flux = 0.0, weighted_qperp_flux = 0.0;
        CGLLFFluxFast(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                      dt_sweep, rkl_weight, eflux, muflux,
                      weighted_qpar_flux, weighted_qperp_flux,
                      qpar_ratio, qperp_ratio);
        if (OwnsHeatFluxDiagnosticFace(m, 1, j, js, je, multilevel,
                                       mblev.d_view(m), nghbr)) {
          const Real area = size.d_view(m).dx1*size.d_view(m).dx3;
          AccumulateCGLLFDiagnosticFace(qstats, qpar_ratio, qperp_ratio, area,
                                        tpar(m,k,j,i) - tpar(m,k,j-1,i),
                                        tperp(m,k,j,i) - tperp(m,k,j-1,i),
                                        weighted_qpar_flux, weighted_qperp_flux);
        }
      } else {
        cgl_lf::ScaledValue weighted_qpar_flux, weighted_qperp_flux;
        CGLLFFlux(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                  dt_sweep, rkl_weight, eflux, muflux, weighted_qpar_flux,
                  weighted_qperp_flux, qpar_ratio, qperp_ratio);
        if (OwnsHeatFluxDiagnosticFace(m, 1, j, js, je, multilevel,
                                       mblev.d_view(m), nghbr)) {
          const Real area = size.d_view(m).dx1*size.d_view(m).dx3;
          AccumulateCGLLFDiagnosticFace(qstats, qpar_ratio, qperp_ratio, area,
                                        tpar(m,k,j,i) - tpar(m,k,j-1,i),
                                        tperp(m,k,j,i) - tperp(m,k,j-1,i),
                                        weighted_qpar_flux, weighted_qperp_flux);
        }
      }
    }
    f2(m,IEN,k,j,i) = eflux;
    f2(m,IAN,k,j,i) = muflux;
    }, Kokkos::Sum<array_sum::GlobalSum>(qstats2));
  }
  AccumulateHeatFluxDiagnostics(qstats2);
  } else {
    CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux2);
    Kokkos::parallel_for("cgl_lf_flux2",
        Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji2),
    KOKKOS_LAMBDA(const int idx) {
    const int m = idx/nkji2;
    const int k = (idx - m*nkji2)/nji2 + ks;
    const int j = (idx - m*nkji2 - (k - ks)*nji2)/ni2 + js;
    const int i = idx - m*nkji2 - (k - ks)*nji2 - (j - js)*ni2 + is;
    const Real tx = VL4Limiter(tpar(m,k,j,i+1) - tpar(m,k,j,i),
                      tpar(m,k,j,i) - tpar(m,k,j,i-1),
                      tpar(m,k,j-1,i+1) - tpar(m,k,j-1,i),
                      tpar(m,k,j-1,i) - tpar(m,k,j-1,i-1))/size.d_view(m).dx1;
    const Real px = VL4Limiter(tperp(m,k,j,i+1) - tperp(m,k,j,i),
                      tperp(m,k,j,i) - tperp(m,k,j,i-1),
                      tperp(m,k,j-1,i+1) - tperp(m,k,j-1,i),
                      tperp(m,k,j-1,i) - tperp(m,k,j-1,i-1))/size.d_view(m).dx1;
    const Real bxg = 0.25*(bmag(m,k,j,i+1) - bmag(m,k,j,i-1) +
                           bmag(m,k,j-1,i+1) - bmag(m,k,j-1,i-1))/size.d_view(m).dx1;
    const Real ty = (tpar(m,k,j,i) - tpar(m,k,j-1,i))/size.d_view(m).dx2;
    const Real py = (tperp(m,k,j,i) - tperp(m,k,j-1,i))/size.d_view(m).dx2;
    const Real byg = (bmag(m,k,j,i) - bmag(m,k,j-1,i))/size.d_view(m).dx2;
    Real tz = 0.0, pz = 0.0, bzg = 0.0;
    if (three_d) {
      tz = VL4Limiter(tpar(m,k+1,j,i) - tpar(m,k,j,i),
                      tpar(m,k,j,i) - tpar(m,k-1,j,i),
                      tpar(m,k+1,j-1,i) - tpar(m,k,j-1,i),
                      tpar(m,k,j-1,i) - tpar(m,k-1,j-1,i))/size.d_view(m).dx3;
      pz = VL4Limiter(tperp(m,k+1,j,i) - tperp(m,k,j,i),
                      tperp(m,k,j,i) - tperp(m,k-1,j,i),
                      tperp(m,k+1,j-1,i) - tperp(m,k,j-1,i),
                      tperp(m,k,j-1,i) - tperp(m,k-1,j-1,i))/size.d_view(m).dx3;
      bzg = 0.25*(bmag(m,k+1,j,i) - bmag(m,k-1,j,i) +
                  bmag(m,k+1,j-1,i) - bmag(m,k-1,j-1,i))/size.d_view(m).dx3;
    }
    const Real bx = 0.5*bcc(m,IBX,k,j-1,i) + 0.5*bcc(m,IBX,k,j,i);
    const Real by = b.x2f(m,k,j,i);
    const Real bz = 0.5*bcc(m,IBZ,k,j-1,i) + 0.5*bcc(m,IBZ,k,j,i);
    CGLLFFaceState face;
    Real eflux = 0.0, muflux = 0.0;
    Real qpar_ratio = 0.0, qperp_ratio = 0.0;
    if (BuildCGLLFFaceState(w(m,IDN,k,j-1,i), w(m,IDN,k,j,i),
                            w(m,IPR,k,j-1,i), w(m,IPR,k,j,i),
                            w(m,IPP,k,j-1,i), w(m,IPP,k,j,i),
                            bx, by, bz, 0.5*bmag(m,k,j-1,i) + 0.5*bmag(m,k,j,i),
                            1, lf_k, local, cpar0, backup, eos, face)) {
      if (fast_arithmetic) {
        Real weighted_qpar_flux = 0.0, weighted_qperp_flux = 0.0;
        CGLLFFluxFast(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                      dt_sweep, rkl_weight, eflux, muflux,
                      weighted_qpar_flux, weighted_qperp_flux,
                      qpar_ratio, qperp_ratio);
      } else {
        cgl_lf::ScaledValue weighted_qpar_flux, weighted_qperp_flux;
        CGLLFFlux(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                  dt_sweep, rkl_weight, eflux, muflux, weighted_qpar_flux,
                  weighted_qperp_flux, qpar_ratio, qperp_ratio);
      }
    }
    f2(m,IEN,k,j,i) = eflux;
    f2(m,IAN,k,j,i) = muflux;
    });
  }
  if (profile_detail_enabled_) {
    Real detail = 0.0;
    {
      CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux2_gradients);
      Kokkos::parallel_reduce("cgl_lf_flux2_profile_gradients",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji2),
      KOKKOS_LAMBDA(const int idx, Real &sum) {
      const int m = idx/nkji2;
      const int k = (idx - m*nkji2)/nji2 + ks;
      const int j = (idx - m*nkji2 - (k - ks)*nji2)/ni2 + js;
      const int i = idx - m*nkji2 - (k - ks)*nji2 - (j - js)*ni2 + is;
      const Real tx = VL4Limiter(tpar(m,k,j,i+1) - tpar(m,k,j,i),
                      tpar(m,k,j,i) - tpar(m,k,j,i-1),
                      tpar(m,k,j-1,i+1) - tpar(m,k,j-1,i),
                      tpar(m,k,j-1,i) - tpar(m,k,j-1,i-1))
                            /size.d_view(m).dx1;
      const Real px = VL4Limiter(tperp(m,k,j,i+1) - tperp(m,k,j,i),
                      tperp(m,k,j,i) - tperp(m,k,j,i-1),
                      tperp(m,k,j-1,i+1) - tperp(m,k,j-1,i),
                      tperp(m,k,j-1,i) - tperp(m,k,j-1,i-1))
                            /size.d_view(m).dx1;
      const Real bxg = 0.25*(bmag(m,k,j,i+1) - bmag(m,k,j,i-1) +
                             bmag(m,k,j-1,i+1) - bmag(m,k,j-1,i-1))
                             /size.d_view(m).dx1;
      const Real ty = (tpar(m,k,j,i) - tpar(m,k,j-1,i))/size.d_view(m).dx2;
      const Real py = (tperp(m,k,j,i) - tperp(m,k,j-1,i))/size.d_view(m).dx2;
      const Real byg = (bmag(m,k,j,i) - bmag(m,k,j-1,i))/size.d_view(m).dx2;
      Real tz = 0.0, pz = 0.0, bzg = 0.0;
      if (three_d) {
        tz = VL4Limiter(tpar(m,k+1,j,i) - tpar(m,k,j,i),
                      tpar(m,k,j,i) - tpar(m,k-1,j,i),
                      tpar(m,k+1,j-1,i) - tpar(m,k,j-1,i),
                      tpar(m,k,j-1,i) - tpar(m,k-1,j-1,i))/size.d_view(m).dx3;
        pz = VL4Limiter(tperp(m,k+1,j,i) - tperp(m,k,j,i),
                      tperp(m,k,j,i) - tperp(m,k-1,j,i),
                      tperp(m,k+1,j-1,i) - tperp(m,k,j-1,i),
                      tperp(m,k,j-1,i) - tperp(m,k-1,j-1,i))/size.d_view(m).dx3;
        bzg = 0.25*(bmag(m,k+1,j,i) - bmag(m,k-1,j,i) +
                    bmag(m,k+1,j-1,i) - bmag(m,k-1,j-1,i))/size.d_view(m).dx3;
      }
      sum += fabs(tx) + fabs(px) + fabs(bxg) + fabs(ty) + fabs(py) + fabs(byg)
             + fabs(tz) + fabs(pz) + fabs(bzg);
      }, Kokkos::Sum<Real>(detail));
    }
    profile_detail_sink_ += detail;
    detail = 0.0;
    {
      CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux2_face_state);
      Kokkos::parallel_reduce("cgl_lf_flux2_profile_face_state",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji2),
      KOKKOS_LAMBDA(const int idx, Real &sum) {
      const int m = idx/nkji2;
      const int k = (idx - m*nkji2)/nji2 + ks;
      const int j = (idx - m*nkji2 - (k - ks)*nji2)/ni2 + js;
      const int i = idx - m*nkji2 - (k - ks)*nji2 - (j - js)*ni2 + is;
      const Real bx = 0.5*bcc(m,IBX,k,j-1,i) + 0.5*bcc(m,IBX,k,j,i);
      const Real by = b.x2f(m,k,j,i);
      const Real bz = 0.5*bcc(m,IBZ,k,j-1,i) + 0.5*bcc(m,IBZ,k,j,i);
      CGLLFFaceState face;
      if (BuildCGLLFFaceState(w(m,IDN,k,j-1,i), w(m,IDN,k,j,i),
                              w(m,IPR,k,j-1,i), w(m,IPR,k,j,i),
                              w(m,IPP,k,j-1,i), w(m,IPP,k,j,i),
                              bx, by, bz, 0.5*bmag(m,k,j-1,i) + 0.5*bmag(m,k,j,i),
                              1, lf_k, local, cpar0, backup, eos, face)) {
        sum += face.cparallel + face.bmag_inv + fabs(face.bhdir) + face.nu;
      }
      }, Kokkos::Sum<Real>(detail));
    }
    profile_detail_sink_ += detail;
    detail = 0.0;
    {
      CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux2_closure);
      Kokkos::parallel_reduce("cgl_lf_flux2_profile_closure",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji2),
      KOKKOS_LAMBDA(const int idx, Real &sum) {
      const int m = idx/nkji2;
      const int k = (idx - m*nkji2)/nji2 + ks;
      const int j = (idx - m*nkji2 - (k - ks)*nji2)/ni2 + js;
      const int i = idx - m*nkji2 - (k - ks)*nji2 - (j - js)*ni2 + is;
      const Real tx = VL4Limiter(tpar(m,k,j,i+1) - tpar(m,k,j,i),
                      tpar(m,k,j,i) - tpar(m,k,j,i-1),
                      tpar(m,k,j-1,i+1) - tpar(m,k,j-1,i),
                      tpar(m,k,j-1,i) - tpar(m,k,j-1,i-1))
                            /size.d_view(m).dx1;
      const Real px = VL4Limiter(tperp(m,k,j,i+1) - tperp(m,k,j,i),
                      tperp(m,k,j,i) - tperp(m,k,j,i-1),
                      tperp(m,k,j-1,i+1) - tperp(m,k,j-1,i),
                      tperp(m,k,j-1,i) - tperp(m,k,j-1,i-1))
                            /size.d_view(m).dx1;
      const Real bxg = 0.25*(bmag(m,k,j,i+1) - bmag(m,k,j,i-1) +
                             bmag(m,k,j-1,i+1) - bmag(m,k,j-1,i-1))
                             /size.d_view(m).dx1;
      const Real ty = (tpar(m,k,j,i) - tpar(m,k,j-1,i))/size.d_view(m).dx2;
      const Real py = (tperp(m,k,j,i) - tperp(m,k,j-1,i))/size.d_view(m).dx2;
      const Real byg = (bmag(m,k,j,i) - bmag(m,k,j-1,i))/size.d_view(m).dx2;
      Real tz = 0.0, pz = 0.0, bzg = 0.0;
      if (three_d) {
        tz = VL4Limiter(tpar(m,k+1,j,i) - tpar(m,k,j,i),
                      tpar(m,k,j,i) - tpar(m,k-1,j,i),
                      tpar(m,k+1,j-1,i) - tpar(m,k,j-1,i),
                      tpar(m,k,j-1,i) - tpar(m,k-1,j-1,i))/size.d_view(m).dx3;
        pz = VL4Limiter(tperp(m,k+1,j,i) - tperp(m,k,j,i),
                      tperp(m,k,j,i) - tperp(m,k-1,j,i),
                      tperp(m,k+1,j-1,i) - tperp(m,k,j-1,i),
                      tperp(m,k,j-1,i) - tperp(m,k-1,j-1,i))/size.d_view(m).dx3;
        bzg = 0.25*(bmag(m,k+1,j,i) - bmag(m,k-1,j,i) +
                    bmag(m,k+1,j-1,i) - bmag(m,k-1,j-1,i))/size.d_view(m).dx3;
      }
      const Real bx = 0.5*bcc(m,IBX,k,j-1,i) + 0.5*bcc(m,IBX,k,j,i);
      const Real by = b.x2f(m,k,j,i);
      const Real bz = 0.5*bcc(m,IBZ,k,j-1,i) + 0.5*bcc(m,IBZ,k,j,i);
      CGLLFFaceState face;
      Real eflux = 0.0, muflux = 0.0;
      Real qpar_ratio = 0.0, qperp_ratio = 0.0;
      if (BuildCGLLFFaceState(w(m,IDN,k,j-1,i), w(m,IDN,k,j,i),
                              w(m,IPR,k,j-1,i), w(m,IPR,k,j,i),
                              w(m,IPP,k,j-1,i), w(m,IPP,k,j,i),
                              bx, by, bz, 0.5*bmag(m,k,j-1,i) + 0.5*bmag(m,k,j,i),
                              1, lf_k, local, cpar0, backup, eos, face)) {
        if (fast_arithmetic) {
          Real weighted_qpar_flux = 0.0, weighted_qperp_flux = 0.0;
          CGLLFFluxFast(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                        dt_sweep, rkl_weight, eflux, muflux,
                        weighted_qpar_flux, weighted_qperp_flux,
                        qpar_ratio, qperp_ratio);
        } else {
          cgl_lf::ScaledValue weighted_qpar_flux, weighted_qperp_flux;
          CGLLFFlux(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                    dt_sweep, rkl_weight, eflux, muflux, weighted_qpar_flux,
                    weighted_qperp_flux, qpar_ratio, qperp_ratio);
        }
        sum += fabs(eflux) + fabs(muflux) + qpar_ratio + qperp_ratio;
      }
      }, Kokkos::Sum<Real>(detail));
    }
    profile_detail_sink_ += detail;
  }
  if (pmy_pack->pmesh->two_d) {
    return;
  }

  auto f3 = f.x3f;
  const int ni3 = ie - is + 1;
  const int nj3 = je - js + 1;
  const int nk3 = ke - ks + 2;
  const int nji3 = nj3*ni3;
  const int nkji3 = nk3*nji3;
  const int nmkji3 = (nmb1 + 1)*nkji3;
  if (collect_heat_flux_diagnostics) {
  array_sum::GlobalSum qstats3;
  {
    CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux3);
    Kokkos::parallel_reduce("cgl_lf_flux3",
        Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji3),
    KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &qstats) {
    const int m = idx/nkji3;
    const int k = (idx - m*nkji3)/nji3 + ks;
    const int j = (idx - m*nkji3 - (k - ks)*nji3)/ni3 + js;
    const int i = idx - m*nkji3 - (k - ks)*nji3 - (j - js)*ni3 + is;
    const Real tx = VL4Limiter(tpar(m,k,j,i+1) - tpar(m,k,j,i),
                      tpar(m,k,j,i) - tpar(m,k,j,i-1),
                      tpar(m,k-1,j,i+1) - tpar(m,k-1,j,i),
                      tpar(m,k-1,j,i) - tpar(m,k-1,j,i-1))/size.d_view(m).dx1;
    const Real px = VL4Limiter(tperp(m,k,j,i+1) - tperp(m,k,j,i),
                      tperp(m,k,j,i) - tperp(m,k,j,i-1),
                      tperp(m,k-1,j,i+1) - tperp(m,k-1,j,i),
                      tperp(m,k-1,j,i) - tperp(m,k-1,j,i-1))/size.d_view(m).dx1;
    const Real bxg = 0.25*(bmag(m,k,j,i+1) - bmag(m,k,j,i-1) +
                           bmag(m,k-1,j,i+1) - bmag(m,k-1,j,i-1))/size.d_view(m).dx1;
    const Real ty = VL4Limiter(tpar(m,k,j+1,i) - tpar(m,k,j,i),
                      tpar(m,k,j,i) - tpar(m,k,j-1,i),
                      tpar(m,k-1,j+1,i) - tpar(m,k-1,j,i),
                      tpar(m,k-1,j,i) - tpar(m,k-1,j-1,i))/size.d_view(m).dx2;
    const Real py = VL4Limiter(tperp(m,k,j+1,i) - tperp(m,k,j,i),
                      tperp(m,k,j,i) - tperp(m,k,j-1,i),
                      tperp(m,k-1,j+1,i) - tperp(m,k-1,j,i),
                      tperp(m,k-1,j,i) - tperp(m,k-1,j-1,i))/size.d_view(m).dx2;
    const Real byg = 0.25*(bmag(m,k,j+1,i) - bmag(m,k,j-1,i) +
                           bmag(m,k-1,j+1,i) - bmag(m,k-1,j-1,i))/size.d_view(m).dx2;
    const Real tz = (tpar(m,k,j,i) - tpar(m,k-1,j,i))/size.d_view(m).dx3;
    const Real pz = (tperp(m,k,j,i) - tperp(m,k-1,j,i))/size.d_view(m).dx3;
    const Real bzg = (bmag(m,k,j,i) - bmag(m,k-1,j,i))/size.d_view(m).dx3;
    const Real bx = 0.5*bcc(m,IBX,k-1,j,i) + 0.5*bcc(m,IBX,k,j,i);
    const Real by = 0.5*bcc(m,IBY,k-1,j,i) + 0.5*bcc(m,IBY,k,j,i);
    const Real bz = b.x3f(m,k,j,i);
    CGLLFFaceState face;
    Real eflux = 0.0, muflux = 0.0;
    Real qpar_ratio = 0.0, qperp_ratio = 0.0;
    if (BuildCGLLFFaceState(w(m,IDN,k-1,j,i), w(m,IDN,k,j,i),
                            w(m,IPR,k-1,j,i), w(m,IPR,k,j,i),
                            w(m,IPP,k-1,j,i), w(m,IPP,k,j,i),
                            bx, by, bz, 0.5*bmag(m,k-1,j,i) + 0.5*bmag(m,k,j,i),
                            2, lf_k, local, cpar0, backup, eos, face)) {
      if (fast_arithmetic) {
        Real weighted_qpar_flux = 0.0, weighted_qperp_flux = 0.0;
        CGLLFFluxFast(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                      dt_sweep, rkl_weight, eflux, muflux,
                      weighted_qpar_flux, weighted_qperp_flux,
                      qpar_ratio, qperp_ratio);
        if (OwnsHeatFluxDiagnosticFace(m, 2, k, ks, ke, multilevel,
                                       mblev.d_view(m), nghbr)) {
          const Real area = size.d_view(m).dx1*size.d_view(m).dx2;
          AccumulateCGLLFDiagnosticFace(qstats, qpar_ratio, qperp_ratio, area,
                                        tpar(m,k,j,i) - tpar(m,k-1,j,i),
                                        tperp(m,k,j,i) - tperp(m,k-1,j,i),
                                        weighted_qpar_flux, weighted_qperp_flux);
        }
      } else {
        cgl_lf::ScaledValue weighted_qpar_flux, weighted_qperp_flux;
        CGLLFFlux(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                  dt_sweep, rkl_weight, eflux, muflux, weighted_qpar_flux,
                  weighted_qperp_flux, qpar_ratio, qperp_ratio);
        if (OwnsHeatFluxDiagnosticFace(m, 2, k, ks, ke, multilevel,
                                       mblev.d_view(m), nghbr)) {
          const Real area = size.d_view(m).dx1*size.d_view(m).dx2;
          AccumulateCGLLFDiagnosticFace(qstats, qpar_ratio, qperp_ratio, area,
                                        tpar(m,k,j,i) - tpar(m,k-1,j,i),
                                        tperp(m,k,j,i) - tperp(m,k-1,j,i),
                                        weighted_qpar_flux, weighted_qperp_flux);
        }
      }
    }
    f3(m,IEN,k,j,i) = eflux;
    f3(m,IAN,k,j,i) = muflux;
    }, Kokkos::Sum<array_sum::GlobalSum>(qstats3));
  }
  AccumulateHeatFluxDiagnostics(qstats3);
  } else {
    CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux3);
    Kokkos::parallel_for("cgl_lf_flux3",
        Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji3),
    KOKKOS_LAMBDA(const int idx) {
    const int m = idx/nkji3;
    const int k = (idx - m*nkji3)/nji3 + ks;
    const int j = (idx - m*nkji3 - (k - ks)*nji3)/ni3 + js;
    const int i = idx - m*nkji3 - (k - ks)*nji3 - (j - js)*ni3 + is;
    const Real tx = VL4Limiter(tpar(m,k,j,i+1) - tpar(m,k,j,i),
                      tpar(m,k,j,i) - tpar(m,k,j,i-1),
                      tpar(m,k-1,j,i+1) - tpar(m,k-1,j,i),
                      tpar(m,k-1,j,i) - tpar(m,k-1,j,i-1))/size.d_view(m).dx1;
    const Real px = VL4Limiter(tperp(m,k,j,i+1) - tperp(m,k,j,i),
                      tperp(m,k,j,i) - tperp(m,k,j,i-1),
                      tperp(m,k-1,j,i+1) - tperp(m,k-1,j,i),
                      tperp(m,k-1,j,i) - tperp(m,k-1,j,i-1))/size.d_view(m).dx1;
    const Real bxg = 0.25*(bmag(m,k,j,i+1) - bmag(m,k,j,i-1) +
                           bmag(m,k-1,j,i+1) - bmag(m,k-1,j,i-1))/size.d_view(m).dx1;
    const Real ty = VL4Limiter(tpar(m,k,j+1,i) - tpar(m,k,j,i),
                      tpar(m,k,j,i) - tpar(m,k,j-1,i),
                      tpar(m,k-1,j+1,i) - tpar(m,k-1,j,i),
                      tpar(m,k-1,j,i) - tpar(m,k-1,j-1,i))/size.d_view(m).dx2;
    const Real py = VL4Limiter(tperp(m,k,j+1,i) - tperp(m,k,j,i),
                      tperp(m,k,j,i) - tperp(m,k,j-1,i),
                      tperp(m,k-1,j+1,i) - tperp(m,k-1,j,i),
                      tperp(m,k-1,j,i) - tperp(m,k-1,j-1,i))/size.d_view(m).dx2;
    const Real byg = 0.25*(bmag(m,k,j+1,i) - bmag(m,k,j-1,i) +
                           bmag(m,k-1,j+1,i) - bmag(m,k-1,j-1,i))/size.d_view(m).dx2;
    const Real tz = (tpar(m,k,j,i) - tpar(m,k-1,j,i))/size.d_view(m).dx3;
    const Real pz = (tperp(m,k,j,i) - tperp(m,k-1,j,i))/size.d_view(m).dx3;
    const Real bzg = (bmag(m,k,j,i) - bmag(m,k-1,j,i))/size.d_view(m).dx3;
    const Real bx = 0.5*bcc(m,IBX,k-1,j,i) + 0.5*bcc(m,IBX,k,j,i);
    const Real by = 0.5*bcc(m,IBY,k-1,j,i) + 0.5*bcc(m,IBY,k,j,i);
    const Real bz = b.x3f(m,k,j,i);
    CGLLFFaceState face;
    Real eflux = 0.0, muflux = 0.0;
    Real qpar_ratio = 0.0, qperp_ratio = 0.0;
    if (BuildCGLLFFaceState(w(m,IDN,k-1,j,i), w(m,IDN,k,j,i),
                            w(m,IPR,k-1,j,i), w(m,IPR,k,j,i),
                            w(m,IPP,k-1,j,i), w(m,IPP,k,j,i),
                            bx, by, bz, 0.5*bmag(m,k-1,j,i) + 0.5*bmag(m,k,j,i),
                            2, lf_k, local, cpar0, backup, eos, face)) {
      if (fast_arithmetic) {
        Real weighted_qpar_flux = 0.0, weighted_qperp_flux = 0.0;
        CGLLFFluxFast(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                      dt_sweep, rkl_weight, eflux, muflux,
                      weighted_qpar_flux, weighted_qperp_flux,
                      qpar_ratio, qperp_ratio);
      } else {
        cgl_lf::ScaledValue weighted_qpar_flux, weighted_qperp_flux;
        CGLLFFlux(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                  dt_sweep, rkl_weight, eflux, muflux, weighted_qpar_flux,
                  weighted_qperp_flux, qpar_ratio, qperp_ratio);
      }
    }
    f3(m,IEN,k,j,i) = eflux;
    f3(m,IAN,k,j,i) = muflux;
    });
  }
  if (profile_detail_enabled_) {
    Real detail = 0.0;
    {
      CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux3_gradients);
      Kokkos::parallel_reduce("cgl_lf_flux3_profile_gradients",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji3),
      KOKKOS_LAMBDA(const int idx, Real &sum) {
      const int m = idx/nkji3;
      const int k = (idx - m*nkji3)/nji3 + ks;
      const int j = (idx - m*nkji3 - (k - ks)*nji3)/ni3 + js;
      const int i = idx - m*nkji3 - (k - ks)*nji3 - (j - js)*ni3 + is;
      const Real tx = VL4Limiter(tpar(m,k,j,i+1) - tpar(m,k,j,i),
                      tpar(m,k,j,i) - tpar(m,k,j,i-1),
                      tpar(m,k-1,j,i+1) - tpar(m,k-1,j,i),
                      tpar(m,k-1,j,i) - tpar(m,k-1,j,i-1))
                            /size.d_view(m).dx1;
      const Real px = VL4Limiter(tperp(m,k,j,i+1) - tperp(m,k,j,i),
                      tperp(m,k,j,i) - tperp(m,k,j,i-1),
                      tperp(m,k-1,j,i+1) - tperp(m,k-1,j,i),
                      tperp(m,k-1,j,i) - tperp(m,k-1,j,i-1))
                            /size.d_view(m).dx1;
      const Real bxg = 0.25*(bmag(m,k,j,i+1) - bmag(m,k,j,i-1) +
                             bmag(m,k-1,j,i+1) - bmag(m,k-1,j,i-1))
                             /size.d_view(m).dx1;
      const Real ty = VL4Limiter(tpar(m,k,j+1,i) - tpar(m,k,j,i),
                      tpar(m,k,j,i) - tpar(m,k,j-1,i),
                      tpar(m,k-1,j+1,i) - tpar(m,k-1,j,i),
                      tpar(m,k-1,j,i) - tpar(m,k-1,j-1,i))
                            /size.d_view(m).dx2;
      const Real py = VL4Limiter(tperp(m,k,j+1,i) - tperp(m,k,j,i),
                      tperp(m,k,j,i) - tperp(m,k,j-1,i),
                      tperp(m,k-1,j+1,i) - tperp(m,k-1,j,i),
                      tperp(m,k-1,j,i) - tperp(m,k-1,j-1,i))
                            /size.d_view(m).dx2;
      const Real byg = 0.25*(bmag(m,k,j+1,i) - bmag(m,k,j-1,i) +
                             bmag(m,k-1,j+1,i) - bmag(m,k-1,j-1,i))
                             /size.d_view(m).dx2;
      const Real tz = (tpar(m,k,j,i) - tpar(m,k-1,j,i))/size.d_view(m).dx3;
      const Real pz = (tperp(m,k,j,i) - tperp(m,k-1,j,i))/size.d_view(m).dx3;
      const Real bzg = (bmag(m,k,j,i) - bmag(m,k-1,j,i))/size.d_view(m).dx3;
      sum += fabs(tx) + fabs(px) + fabs(bxg) + fabs(ty) + fabs(py) + fabs(byg)
             + fabs(tz) + fabs(pz) + fabs(bzg);
      }, Kokkos::Sum<Real>(detail));
    }
    profile_detail_sink_ += detail;
    detail = 0.0;
    {
      CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux3_face_state);
      Kokkos::parallel_reduce("cgl_lf_flux3_profile_face_state",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji3),
      KOKKOS_LAMBDA(const int idx, Real &sum) {
      const int m = idx/nkji3;
      const int k = (idx - m*nkji3)/nji3 + ks;
      const int j = (idx - m*nkji3 - (k - ks)*nji3)/ni3 + js;
      const int i = idx - m*nkji3 - (k - ks)*nji3 - (j - js)*ni3 + is;
      const Real bx = 0.5*bcc(m,IBX,k-1,j,i) + 0.5*bcc(m,IBX,k,j,i);
      const Real by = 0.5*bcc(m,IBY,k-1,j,i) + 0.5*bcc(m,IBY,k,j,i);
      const Real bz = b.x3f(m,k,j,i);
      CGLLFFaceState face;
      if (BuildCGLLFFaceState(w(m,IDN,k-1,j,i), w(m,IDN,k,j,i),
                              w(m,IPR,k-1,j,i), w(m,IPR,k,j,i),
                              w(m,IPP,k-1,j,i), w(m,IPP,k,j,i),
                              bx, by, bz, 0.5*bmag(m,k-1,j,i) + 0.5*bmag(m,k,j,i),
                              2, lf_k, local, cpar0, backup, eos, face)) {
        sum += face.cparallel + face.bmag_inv + fabs(face.bhdir) + face.nu;
      }
      }, Kokkos::Sum<Real>(detail));
    }
    profile_detail_sink_ += detail;
    detail = 0.0;
    {
      CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux3_closure);
      Kokkos::parallel_reduce("cgl_lf_flux3_profile_closure",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji3),
      KOKKOS_LAMBDA(const int idx, Real &sum) {
      const int m = idx/nkji3;
      const int k = (idx - m*nkji3)/nji3 + ks;
      const int j = (idx - m*nkji3 - (k - ks)*nji3)/ni3 + js;
      const int i = idx - m*nkji3 - (k - ks)*nji3 - (j - js)*ni3 + is;
      const Real tx = VL4Limiter(tpar(m,k,j,i+1) - tpar(m,k,j,i),
                      tpar(m,k,j,i) - tpar(m,k,j,i-1),
                      tpar(m,k-1,j,i+1) - tpar(m,k-1,j,i),
                      tpar(m,k-1,j,i) - tpar(m,k-1,j,i-1))
                            /size.d_view(m).dx1;
      const Real px = VL4Limiter(tperp(m,k,j,i+1) - tperp(m,k,j,i),
                      tperp(m,k,j,i) - tperp(m,k,j,i-1),
                      tperp(m,k-1,j,i+1) - tperp(m,k-1,j,i),
                      tperp(m,k-1,j,i) - tperp(m,k-1,j,i-1))
                            /size.d_view(m).dx1;
      const Real bxg = 0.25*(bmag(m,k,j,i+1) - bmag(m,k,j,i-1) +
                             bmag(m,k-1,j,i+1) - bmag(m,k-1,j,i-1))
                             /size.d_view(m).dx1;
      const Real ty = VL4Limiter(tpar(m,k,j+1,i) - tpar(m,k,j,i),
                      tpar(m,k,j,i) - tpar(m,k,j-1,i),
                      tpar(m,k-1,j+1,i) - tpar(m,k-1,j,i),
                      tpar(m,k-1,j,i) - tpar(m,k-1,j-1,i))
                            /size.d_view(m).dx2;
      const Real py = VL4Limiter(tperp(m,k,j+1,i) - tperp(m,k,j,i),
                      tperp(m,k,j,i) - tperp(m,k,j-1,i),
                      tperp(m,k-1,j+1,i) - tperp(m,k-1,j,i),
                      tperp(m,k-1,j,i) - tperp(m,k-1,j-1,i))
                            /size.d_view(m).dx2;
      const Real byg = 0.25*(bmag(m,k,j+1,i) - bmag(m,k,j-1,i) +
                             bmag(m,k-1,j+1,i) - bmag(m,k-1,j-1,i))
                             /size.d_view(m).dx2;
      const Real tz = (tpar(m,k,j,i) - tpar(m,k-1,j,i))/size.d_view(m).dx3;
      const Real pz = (tperp(m,k,j,i) - tperp(m,k-1,j,i))/size.d_view(m).dx3;
      const Real bzg = (bmag(m,k,j,i) - bmag(m,k-1,j,i))/size.d_view(m).dx3;
      const Real bx = 0.5*bcc(m,IBX,k-1,j,i) + 0.5*bcc(m,IBX,k,j,i);
      const Real by = 0.5*bcc(m,IBY,k-1,j,i) + 0.5*bcc(m,IBY,k,j,i);
      const Real bz = b.x3f(m,k,j,i);
      CGLLFFaceState face;
      Real eflux = 0.0, muflux = 0.0;
      Real qpar_ratio = 0.0, qperp_ratio = 0.0;
      if (BuildCGLLFFaceState(w(m,IDN,k-1,j,i), w(m,IDN,k,j,i),
                              w(m,IPR,k-1,j,i), w(m,IPR,k,j,i),
                              w(m,IPP,k-1,j,i), w(m,IPP,k,j,i),
                              bx, by, bz, 0.5*bmag(m,k-1,j,i) + 0.5*bmag(m,k,j,i),
                              2, lf_k, local, cpar0, backup, eos, face)) {
        if (fast_arithmetic) {
          Real weighted_qpar_flux = 0.0, weighted_qperp_flux = 0.0;
          CGLLFFluxFast(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                        dt_sweep, rkl_weight, eflux, muflux,
                        weighted_qpar_flux, weighted_qperp_flux,
                        qpar_ratio, qperp_ratio);
        } else {
          cgl_lf::ScaledValue weighted_qpar_flux, weighted_qperp_flux;
          CGLLFFlux(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                    dt_sweep, rkl_weight, eflux, muflux, weighted_qpar_flux,
                    weighted_qperp_flux, qpar_ratio, qperp_ratio);
        }
        sum += fabs(eflux) + fabs(muflux) + qpar_ratio + qperp_ratio;
      }
      }, Kokkos::Sum<Real>(detail));
    }
    profile_detail_sink_ += detail;
  }
}

void CGLLandauFluid::AdvanceHeatFluxWorkDiagnostics(
    const parabolic::RKL2Coefficients &coeffs, int stage, int nstages) {
  CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_work_diagnostics);
  if (profile_enabled_) {
    profile_last_nstages_ = nstages;
    if (nstages > profile_max_nstages_) {
      profile_max_nstages_ = nstages;
    }
  }
  if (stage == 1) {
    sweep_qpar_work_ = 0.0;
    sweep_qperp_work_ = 0.0;
    sweep_qpar_work1_ = 0.0;
    sweep_qperp_work1_ = 0.0;
    sweep_qpar_work2_ = 0.0;
    sweep_qperp_work2_ = 0.0;
    sweep_qpar_rhs_ = 0.0;
    sweep_qperp_rhs_ = 0.0;
  }
  if (diagnostics_mode == CGLLFDiagnosticsMode::none) {
    ResetHeatFluxDiagnostics();
    return;
  }

  const Real qpar_rhs = stage_qpar_work_;
  const Real qperp_rhs = stage_qperp_work_;
  sweep_qpar_work2_ = sweep_qpar_work1_;
  sweep_qperp_work2_ = sweep_qperp_work1_;
  sweep_qpar_work1_ = sweep_qpar_work_;
  sweep_qperp_work1_ = sweep_qperp_work_;
  // The per-sweep diagnostic starts at zero, so the RKL2 S0 term vanishes.
  const Real first_rhs_coeff =
      cgl_lf::CachedRHSCoefficient(coeffs.gammaj_tilde, nstages);
  sweep_qpar_work_ = cgl_lf::WeightedRKL2Update(
      coeffs.muj, sweep_qpar_work1_, coeffs.nuj, sweep_qpar_work2_,
      0.0, 0.0, first_rhs_coeff, sweep_qpar_rhs_, qpar_rhs);
  sweep_qperp_work_ = cgl_lf::WeightedRKL2Update(
      coeffs.muj, sweep_qperp_work1_, coeffs.nuj, sweep_qperp_work2_,
      0.0, 0.0, first_rhs_coeff, sweep_qperp_rhs_, qperp_rhs);
  if (stage == 1) {
    sweep_qpar_rhs_ = qpar_rhs;
    sweep_qperp_rhs_ = qperp_rhs;
  }
  if (stage == nstages) {
    diagnostics.qpar_work += sweep_qpar_work_;
    diagnostics.qperp_work += sweep_qperp_work_;
  }
}

void CGLLandauFluid::AdvancePressureWorkDiagnostics(Real beta_dt, Real gam0, Real gam1,
                                                     int stage, Real pressure_power,
                                                     Real anisotropic_power) {
  if (stage == 1) {
    pressure_work_cycle_start_ = diagnostics.pressure_work;
    anisotropic_pressure_work_cycle_start_ = diagnostics.anisotropic_pressure_work;
  }
  diagnostics.pressure_work =
      gam0*diagnostics.pressure_work + gam1*pressure_work_cycle_start_
      + beta_dt*pressure_power;
  diagnostics.anisotropic_pressure_work =
      gam0*diagnostics.anisotropic_pressure_work
      + gam1*anisotropic_pressure_work_cycle_start_
      + beta_dt*anisotropic_power;
}

void CGLLandauFluid::NewTimeStep(const DvceArray5D<Real> &w,
                                 const DvceArray5D<Real> &bcc,
                                 const DvceFaceFld4D<Real> &b,
                                 const EOS_Data &eos_in) {
  CGLLFProfileRegion profile(this, CGLLFProfileBucket::timestep_reduction);
  const EOS_Data eos = eos_in;
  // Limiter scattering and heat-flux caps can only reduce each face coefficient.
  // The row envelope keeps neither reduction, so switching their branches does
  // not invalidate it at fixed temperature, density and magnetic field.
  EOS_Data bound_eos = eos;
  bound_eos.mlim = bound_eos.flim = bound_eos.backup_lim = false;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, nx1 = indcs.nx1;
  const int js = indcs.js, nx2 = indcs.nx2;
  const int ks = indcs.ks, nx3 = indcs.nx3;
  const int nmb = pmy_pack->nmb_thispack;
  const int n1 = nx1 + 2*indcs.ng;
  const int n2 = (nx2 > 1) ? nx2 + 2*indcs.ng : 1;
  const int n3 = (nx3 > 1) ? nx3 + 2*indcs.ng : 1;
  if (timestep_bmag_.extent(0) != static_cast<std::size_t>(nmb) ||
      timestep_bmag_.extent(1) != static_cast<std::size_t>(n3) ||
      timestep_bmag_.extent(2) != static_cast<std::size_t>(n2) ||
      timestep_bmag_.extent(3) != static_cast<std::size_t>(n1)) {
    Kokkos::realloc(timestep_bmag_, nmb, n3, n2, n1);
    Kokkos::realloc(timestep_tpar_, nmb, n3, n2, n1);
    Kokkos::realloc(timestep_tperp_, nmb, n3, n2, n1);
  }
  auto bmag = timestep_bmag_;
  auto tpar = timestep_tpar_, tperp = timestep_tperp_;
  // This independent scratch is fresh after RK, restriction or prolongation;
  // it does not invalidate the heat-flux primitive-refresh cache.
  par_for("cgl_lf_timestep_bmag", DevExeSpace(), 0, nmb-1, 0, n3-1,
          0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    bmag(m,k,j,i) = ScaledMagneticMagnitude(
        bcc(m,IBX,k,j,i), bcc(m,IBY,k,j,i), bcc(m,IBZ,k,j,i));
    const Real rho = fmax(w(m,IDN,k,j,i), eos.dfloor);
    tpar(m,k,j,i) = w(m,IPR,k,j,i)/rho;
    tperp(m,k,j,i) = w(m,IPP,k,j,i)/rho;
  });
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  const int ndim = pmy_pack->pmesh->three_d ? 3 :
                   (pmy_pack->pmesh->multi_d ? 2 : 1);
  const Real kpar = lf_k_parallel;
  const bool local = lf_coeff_local;
  const Real cpar0 = lf_c_parallel0;
  auto size = pmy_pack->pmb->mb_size;
  dtnew = static_cast<Real>(std::numeric_limits<float>::max());
  Kokkos::parallel_reduce("cgl_lf_newdt", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int &idx, Real &min_dt) {
    const int m = idx/nkji;
    const int k = (idx - m*nkji)/nji + ks;
    const int j = (idx - m*nkji - (k-ks)*nji)/nx1 + js;
    const int i = idx - m*nkji - (k-ks)*nji - (j-js)*nx1 + is;
    const Real rho = fmax(w(m,IDN,k,j,i), eos.dfloor);
    const Real bi = bmag(m,k,j,i);
    const Real dx[3] = {size.d_view(m).dx1, size.d_view(m).dx2,
                        size.d_view(m).dx3};
    Real row_parallel = 0.0, row_perp = 0.0;
    Real log_row_parallel = -std::numeric_limits<Real>::infinity();
    Real log_row_perp = log_row_parallel;
    bool scaled_rows = false;
    bool invalid_state = false;
    for (int dir=0; dir<ndim; ++dir) {
      for (int side=-1; side<=1; side+=2) {
        // l/r always follow the positive coordinate direction.
        const int il = i + ((dir == 0 && side < 0) ? -1 : 0);
        const int jl = j + ((dir == 1 && side < 0) ? -1 : 0);
        const int kl = k + ((dir == 2 && side < 0) ? -1 : 0);
        const int ir = il + (dir == 0), jr = jl + (dir == 1);
        const int kr = kl + (dir == 2);
        Real bv[3] = {
            0.5*bcc(m,IBX,kl,jl,il) + 0.5*bcc(m,IBX,kr,jr,ir),
            0.5*bcc(m,IBY,kl,jl,il) + 0.5*bcc(m,IBY,kr,jr,ir),
            0.5*bcc(m,IBZ,kl,jl,il) + 0.5*bcc(m,IBZ,kr,jr,ir)};
        // The actual staggered normal field may exceed Bbar by any factor;
        // replacing this with a cell-centered average misses checkerboards.
        bv[dir] = (dir == 0) ? b.x1f(m,kr,jr,ir) :
                  ((dir == 1) ? b.x2f(m,kr,jr,ir) : b.x3f(m,kr,jr,ir));
        const Real bbar = 0.5*bmag(m,kl,jl,il) + 0.5*bmag(m,kr,jr,ir);
        CGLLFFaceState face;
        if (!BuildCGLLFFaceState(
            w(m,IDN,kl,jl,il), w(m,IDN,kr,jr,ir),
            w(m,IPR,kl,jl,il), w(m,IPR,kr,jr,ir),
            w(m,IPP,kl,jl,il), w(m,IPP,kr,jr,ir),
            bv[0], bv[1], bv[2], bbar, dir, kpar, local, cpar0,
            false, bound_eos, face)) continue;
        if (!(face.cparallel > 0.0) || !Kokkos::isfinite(face.cparallel)) {
          invalid_state = true;
          continue;
        }
        const Real nu = fmax(eos.nu_coll, static_cast<Real>(0.0));
        if (nu == std::numeric_limits<Real>::infinity()) continue;
        const Real bn = fabs(face.bhdir);
        if (bv[dir] == 0.0) continue;
        const Real bh[3] = {face.bhx, face.bhy, face.bhz};
        // Each transverse VL4 derivative has six coefficients, including the
        // unequal-slope central coefficients. Its row norm can approach 8/h_t;
        // the secant theta<=1 representation does not bound this Jacobian.
        Real grad_parallel = 2.0*bn/dx[dir];
        Real grad_perp = grad_parallel;
        Real derivative_parallel[3] = {0.0,0.0,0.0};
        Real derivative_perp[3] = {0.0,0.0,0.0};
        Real grad_b = face.bhdir*(bmag(m,kr,jr,ir)-bmag(m,kl,jl,il))/dx[dir];
        for (int t=0; t<ndim; ++t) {
          if (t == dir || bv[t] == 0.0) continue;
          const int di = (t == 0), dj = (t == 1), dk = (t == 2);
          derivative_parallel[t] = CGLLFVL4DerivativeNorm(
              tpar(m,kr+dk,jr+dj,ir+di)-tpar(m,kr,jr,ir),
              tpar(m,kr,jr,ir)-tpar(m,kr-dk,jr-dj,ir-di),
              tpar(m,kl+dk,jl+dj,il+di)-tpar(m,kl,jl,il),
              tpar(m,kl,jl,il)-tpar(m,kl-dk,jl-dj,il-di));
          derivative_perp[t] = CGLLFVL4DerivativeNorm(
              tperp(m,kr+dk,jr+dj,ir+di)-tperp(m,kr,jr,ir),
              tperp(m,kr,jr,ir)-tperp(m,kr-dk,jr-dj,ir-di),
              tperp(m,kl+dk,jl+dj,il+di)-tperp(m,kl,jl,il),
              tperp(m,kl,jl,il)-tperp(m,kl-dk,jl-dj,il-di));
          grad_parallel += fabs(bh[t])*derivative_parallel[t]/dx[t];
          grad_perp += fabs(bh[t])*derivative_perp[t]/dx[t];
          const Real gb = 0.25*(
              bmag(m,kr+dk,jr+dj,ir+di)-bmag(m,kr-dk,jr-dj,ir-di) +
              bmag(m,kl+dk,jl+dj,il+di)-bmag(m,kl-dk,jl-dj,il-di))/dx[t];
          grad_b += bh[t]*gb;
        }
        const Real pressure_ratio = face.pperp/face.ppar;
        // Differentiate Tperp_f*(1-Tperp_f/Tpar_f), including its
        // reverse coupling to Tpar. Both density-weighted face means
        // have nonnegative weights whose sum is one.
        const Real drift = (fabs(1.0-2.0*pressure_ratio) +
                            pressure_ratio*pressure_ratio)*fabs(grad_b)/bbar;
        const Real cp = face.cparallel;
        Real chi_parallel = cgl::kSqrtEightOverPi*cp/
            (kpar + (cgl::kThreePiMinusEight*nu/cgl::kSqrtEightPi)/cp);
        Real chi_perp = cgl::kSqrtTwoOverPi*cp/
            (kpar + (nu/(0.5*cgl::kSqrtTwoPi))/cp);
        if (!Kokkos::isfinite(chi_parallel) || (chi_parallel == 0.0 && cp > 0.0)) {
          const Real response = -cgl::ParallelHeatFluxRatio(cp,1.0,1.0,kpar,nu,1.0);
          chi_parallel = cgl::PositiveProduct4(cgl::kSqrtEightOverPi,cp,response,1.0);
        }
        if (!Kokkos::isfinite(chi_perp) || (chi_perp == 0.0 && cp > 0.0)) {
          const Real response = -cgl::PerpendicularHeatFluxRatio(
              cp,1.0,1.0,1.0,0.0,kpar,nu,1.0,0.0);
          chi_perp = cgl::PositiveProduct4(cgl::kSqrtTwoOverPi,cp,response,1.0);
        }
        const Real para_face = cgl::PositiveProduct4(chi_parallel,bn,grad_parallel,1.0);
        const Real perp_face = cgl::PositiveProduct4(chi_perp,bn,grad_perp+drift,1.0);
        const Real beta = bi/bbar;
        const Real factor = (face.rho/rho)/dx[dir];
        const Real face_parallel = factor*(para_face + 2.0*fabs(1.0-beta)*perp_face);
        const Real face_perp = factor*beta*perp_face;
        const bool direct_face = Kokkos::isfinite(para_face) && para_face > 0.0 &&
            Kokkos::isfinite(perp_face) && perp_face > 0.0 &&
            Kokkos::isfinite(factor) && factor > 0.0 &&
            Kokkos::isfinite(face_parallel) && face_parallel > 0.0 &&
            Kokkos::isfinite(face_perp) && (face_perp > 0.0 || beta == 0.0);
        if (!scaled_rows && direct_face &&
            Kokkos::isfinite(row_parallel+face_parallel) &&
            Kokkos::isfinite(row_perp+face_perp)) {
          row_parallel += face_parallel;
          row_perp += face_perp;
          continue;
        }
        if (!scaled_rows) {
          log_row_parallel = CGLLFLogAbs(row_parallel);
          log_row_perp = CGLLFLogAbs(row_perp);
          scaled_rows = true;
        }
        Real log_face_parallel, log_face_perp;
        if (direct_face) {
          log_face_parallel = CGLLFLogAbs(face_parallel);
          log_face_perp = CGLLFLogAbs(face_perp);
        } else {
          // Keep conductivity, normalization, density and spacing together;
          // a finite row need not have a representable intermediate chi or b_n.
          const Real log_bar = CGLLFLogAbs(bbar);
          const Real log_bn = CGLLFLogAbs(bv[dir])-log_bar;
          Real log_grad_parallel = 1.0+log_bn-CGLLFLogAbs(dx[dir]);
          Real log_grad_perp = log_grad_parallel;
          Real log_grad_b = log_bn +
              CGLLFLogAbs(bmag(m,kr,jr,ir)-bmag(m,kl,jl,il))-CGLLFLogAbs(dx[dir]);
          for (int t=0; t<ndim; ++t) {
            if (t == dir || bv[t] == 0.0) continue;
            const int di=(t==0), dj=(t==1), dk=(t==2);
            const Real log_bt = CGLLFLogAbs(bv[t])-log_bar;
            const Real log_h = CGLLFLogAbs(dx[t]);
            log_grad_parallel = CGLLFLogAdd(log_grad_parallel,
                log_bt+CGLLFLogAbs(derivative_parallel[t])-log_h);
            log_grad_perp = CGLLFLogAdd(log_grad_perp,
                log_bt+CGLLFLogAbs(derivative_perp[t])-log_h);
            // Absolute magnetic-gradient sums remain an upper bound when
            // individual centered differences overflow or cancel.
            const Real log_gbt = CGLLFLogAdd(
                CGLLFLogAbs(bmag(m,kr+dk,jr+dj,ir+di)-bmag(m,kr-dk,jr-dj,ir-di)),
                CGLLFLogAbs(bmag(m,kl+dk,jl+dj,il+di)-bmag(m,kl-dk,jl-dj,il-di)))
                -2.0-log_h;
            log_grad_b = CGLLFLogAdd(log_grad_b,log_bt+log_gbt);
          }
          // The cold-path factor (1+r)^2 bounds |1-2r|+r^2 without
          // materializing an overflowing pressure ratio or squaring it.
          const Real log_ratio = CGLLFLogAbs(face.pperp)-CGLLFLogAbs(face.ppar);
          const Real log_drift = 2.0*CGLLFLogAdd(0.0,log_ratio)+log_grad_b-log_bar;
          const Real log_factor = CGLLFLogAbs(face.rho)-CGLLFLogAbs(rho)-CGLLFLogAbs(dx[dir]);
          const Real log_para = CGLLFLogDiffusivity(cp,kpar,nu,cgl::kSqrtEightOverPi,
              cgl::kThreePiMinusEight/cgl::kSqrtEightPi)+log_bn+log_grad_parallel+log_factor;
          const Real log_perp = CGLLFLogDiffusivity(cp,kpar,nu,cgl::kSqrtTwoOverPi,
              2.0/cgl::kSqrtTwoPi)+log_bn+CGLLFLogAdd(log_grad_perp,log_drift)+log_factor;
          log_face_parallel = CGLLFLogAdd(log_para,
              CGLLFLogScale(log_perp,2.0*fabs(1.0-beta)));
          log_face_perp = CGLLFLogScale(log_perp,beta);
        }
        log_row_parallel = CGLLFLogAdd(log_row_parallel,log_face_parallel);
        log_row_perp = CGLLFLogAdd(log_row_perp,log_face_perp);
      }
    }
    // A Gershgorin/row-sum bound for the local temperature Jacobian with
    // density, magnetic field, closure coefficients and cap scales frozen.
    // It is not a proof of nonlinear or composite-AMR RKL2 stability.
    if (invalid_state || Kokkos::isnan(log_row_parallel) || Kokkos::isnan(log_row_perp)) {
      min_dt = 0.0;
    } else if (scaled_rows) {
      const Real log_bound = fmax(log_row_parallel,log_row_perp);
      min_dt = fmin(min_dt,Kokkos::exp2(1.0-log_bound));
    } else {
      const Real bound = fmax(row_parallel,row_perp);
      if (bound > 0.0) min_dt = fmin(min_dt,2.0/bound);
    }
  }, Kokkos::Min<Real>(dtnew));
}

void CGLLandauFluid::RecordAdmissibility(const DvceArray5D<Real> &u,
                                         const DvceArray5D<Real> &w,
                                         const DvceArray5D<Real> &bcc,
                                         const EOS_Data &eos_in,
                                         int dfloor_delta, int pfloor_delta,
                                         const char *sweep_name,
                                         int stage, int nstages, bool wall_checkpoint) {
  if (wall_checkpoint && !strict_admissibility) return;
  CGLLFProfileRegion profile(this, CGLLFProfileBucket::admissibility);
  const EOS_Data eos = eos_in;
  const bool backup = effective_backup_limiter;
  if (profile_enabled_ && !wall_checkpoint) {
    profile_last_nstages_ = nstages;
    if (nstages > profile_max_nstages_) {
      profile_max_nstages_ = nstages;
    }
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, nx1 = indcs.nx1;
  const int js = indcs.js, nx2 = indcs.nx2;
  const int ks = indcs.ks, nx3 = indcs.nx3;
  const int nmkji = pmy_pack->nmb_thispack*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji = nx2*nx1;
  int nonfinite = 0, nonpositive = 0, mirror = 0, firehose = 0, hard_bound = 0;
  Kokkos::parallel_reduce("cgl_lf_admissibility",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int idx, int &nbad, int &nnonpos, int &nmirror,
                int &nfirehose, int &nhard) {
    const int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    const int i = idx - m*nkji - k*nji - j*nx1 + is;
    k += ks;
    j += js;
    const Real rho = w(m,IDN,k,j,i);
    const Real ppar = w(m,IPR,k,j,i);
    const Real pperp = w(m,IPP,k,j,i);
    const Real energy = u(m,IEN,k,j,i);
    const Real moment = u(m,IAN,k,j,i);
    if (!Kokkos::isfinite(rho) || !Kokkos::isfinite(ppar) ||
        !Kokkos::isfinite(pperp) || !Kokkos::isfinite(energy) ||
        !Kokkos::isfinite(moment)) {
      ++nbad;
    }
    if (rho <= 0.0 || ppar <= 0.0 || pperp <= 0.0) {
      ++nnonpos;
    }
    const Real bsqr = SQR(bcc(m,IBX,k,j,i)) + SQR(bcc(m,IBY,k,j,i))
                     + SQR(bcc(m,IBZ,k,j,i));
    const Real paniso = pperp - ppar;
    if (eos.mlim && cgl::MirrorLimiterActive(paniso, bsqr, eos)) ++nmirror;
    if (eos.flim && cgl::FirehoseLimiterActive(paniso, bsqr, eos)) ++nfirehose;
    if (cgl::HardBoundViolated(paniso, bsqr, eos, backup)) ++nhard;
  }, Kokkos::Sum<int>(nonfinite), Kokkos::Sum<int>(nonpositive),
     Kokkos::Sum<int>(mirror), Kokkos::Sum<int>(firehose),
     Kokkos::Sum<int>(hard_bound));

  // Count every unprojected LF stage, including its hard-wall crossings. Entry
  // and post-wall checks validate the split state without counting extra stages.
  if (!wall_checkpoint) {
    diagnostics.nstage += static_cast<std::uint64_t>(nmkji);
    diagnostics.dfloor += static_cast<std::uint64_t>(dfloor_delta);
    diagnostics.pfloor += static_cast<std::uint64_t>(pfloor_delta);
    diagnostics.nonfinite += static_cast<std::uint64_t>(nonfinite);
    diagnostics.nonpositive += static_cast<std::uint64_t>(nonpositive);
    diagnostics.mirror += static_cast<std::uint64_t>(mirror);
    diagnostics.firehose += static_cast<std::uint64_t>(firehose);
    diagnostics.hard_bound += static_cast<std::uint64_t>(hard_bound);
  }

  // LF changes anisotropy while B is frozen and can cross a hard wall. The
  // scheduled sweep-end projection must restore it before hyperbolic evolution;
  // floors and invalid pressures remain fatal at every intermediate stage.
  if (strict_admissibility &&
      (dfloor_delta > 0 || pfloor_delta > 0 || nonfinite > 0 ||
       nonpositive > 0 || (wall_checkpoint && hard_bound > 0))) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "CGL Landau-fluid strict admissibility failed "
              << (wall_checkpoint ? "at a split-sweep wall checkpoint: "
                                  : "after a split stage: ")
              << "sweep=" << sweep_name << " stage=" << stage << "/" << nstages
              << " dfloor=" << dfloor_delta << " pfloor=" << pfloor_delta
              << " nonfinite=" << nonfinite << " nonpositive=" << nonpositive
              << " hard_bound=" << hard_bound << std::endl;
    std::exit(EXIT_FAILURE);
  }
}
