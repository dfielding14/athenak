//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file q029_hall_bell_linear.cpp
//! \brief Source-local Q-029 experimental Hall-Bell launch preparation.

#include <algorithm>
#include <cmath>

#if !defined(Q029_HALL_BELL_LINEAR_HOST_CONTRACT)
#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "pgen/pgen.hpp"
#endif

#if defined(Q029_HALL_BELL_LINEAR_HOST_CONTRACT) && \
    !defined(Q023_PAPER_BELL_LINEAR_HOST_CONTRACT)
#define Q023_PAPER_BELL_LINEAR_HOST_CONTRACT 1
#define Q029_UNDEF_Q023_PAPER_BELL_LINEAR_HOST_CONTRACT 1
#endif
#include "q023_paper_bell_linear.hpp"
#if defined(Q029_UNDEF_Q023_PAPER_BELL_LINEAR_HOST_CONTRACT)
#undef Q023_PAPER_BELL_LINEAR_HOST_CONTRACT
#undef Q029_UNDEF_Q023_PAPER_BELL_LINEAR_HOST_CONTRACT
#endif

namespace q029_hall_bell_linear {

using q023_paper_bell_linear::Add;
using q023_paper_bell_linear::Basis;
using q023_paper_bell_linear::Cross;
using q023_paper_bell_linear::Dot;
using q023_paper_bell_linear::EigenmodeAtPhase;
using q023_paper_bell_linear::ModeBasis;
using q023_paper_bell_linear::ModeParameters;
using q023_paper_bell_linear::ModeSample;
using q023_paper_bell_linear::NormalizedGrowthRate;
using q023_paper_bell_linear::Scale;
using q023_paper_bell_linear::Vector3;
using q023_paper_bell_linear::VectorPotentialAt;

inline double CurrentDensity(const double ppc, const double qscale,
                             const double charge, const double stream_speed) {
  return ppc*qscale*charge*stream_speed;
}

inline double ChiH(const double alpha_h, const double j_cr,
                   const double u_a, const double b_g) {
  return alpha_h*j_cr/(u_a*b_g);
}

inline bool IsPositivePreparedChiH(const double chi_h) {
  return chi_h > 0.0 &&
      (std::abs(chi_h - 0.25) <= 1.0e-12 ||
       std::abs(chi_h - 0.5) <= 1.0e-12 ||
       std::abs(chi_h - 1.0) <= 1.0e-12);
}

}  // namespace q029_hall_bell_linear

#if !defined(Q029_HALL_BELL_LINEAR_HOST_CONTRACT)
namespace {

[[noreturn]] void Q029Fatal(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << message << std::endl;
  std::exit(EXIT_FAILURE);
}

void Q029RequireClose(const std::string &label, const Real measured,
                      const Real expected) {
  if (!std::isfinite(measured) || !std::isfinite(expected)) {
    Q029Fatal(label + " must be finite for the Q-029 Hall-Bell preparation contract");
  }
  const Real scale = std::max(static_cast<Real>(1.0), std::abs(expected));
  if (std::abs(measured - expected) > static_cast<Real>(1.0e-12)*scale) {
    Q029Fatal(label + " does not match the Q-029 Hall-Bell preparation contract");
  }
}

void Q029RequireFinite(const std::string &label, const Real value) {
  if (!std::isfinite(value)) {
    Q029Fatal(label + " must be finite for the Q-029 Hall-Bell preparation contract");
  }
}

void Q029RequireString(ParameterInput *pin, const std::string &block,
                       const std::string &name, const std::string &expected) {
  if (pin->GetString(block, name).compare(expected) != 0) {
    Q029Fatal("<" + block + ">/" + name +
              " does not match the Q-029 Hall-Bell preparation contract");
  }
}

void Q029RequireBoolean(ParameterInput *pin, const std::string &block,
                        const std::string &name, const bool expected) {
  if (pin->GetBoolean(block, name) != expected) {
    Q029Fatal("<" + block + ">/" + name +
              " does not match the Q-029 Hall-Bell preparation contract");
  }
}

}  // namespace

void ProblemGenerator::Q029HallBellLinear(ParameterInput *pin, const bool restart) {
  using q029_hall_bell_linear::Basis;
  using q029_hall_bell_linear::ChiH;
  using q029_hall_bell_linear::CurrentDensity;
  using q029_hall_bell_linear::IsPositivePreparedChiH;
  using q029_hall_bell_linear::ModeBasis;
  using q029_hall_bell_linear::ModeParameters;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->pmhd == nullptr || pmbp->ppart == nullptr) {
    Q029Fatal("q029_hall_bell_linear requires MHD and particles");
  }
  if (!pmbp->pmhd->peos->eos_data.is_ideal) {
    Q029Fatal("q029_hall_bell_linear requires ideal-MHD EOS");
  }
  if (pmy_mesh_->multilevel) {
    Q029Fatal("q029_hall_bell_linear source-local preparation rejects AMR/SMR");
  }

  const std::string block = "q029_hall_bell_linear";
  const int dimension = pin->GetInteger(block, "dimension");
  const bool exact_topology =
      ((dimension == 1 || dimension == 2) && !pmy_mesh_->one_d &&
       pmy_mesh_->two_d && !pmy_mesh_->three_d) ||
      (dimension == 3 && !pmy_mesh_->one_d && !pmy_mesh_->two_d &&
       pmy_mesh_->three_d);
  const int expected_global_nx1 =
      (dimension == 1) ? 32 : ((dimension == 2) ? 64 : 128);
  const int expected_global_nx2 =
      (dimension == 1) ? 4 : ((dimension == 2) ? 32 : 64);
  const int expected_global_nx3 = (dimension == 3) ? 32 : 1;
  const int expected_block_nx1 = 32;
  const int expected_block_nx2 = (dimension == 1) ? 4 : 32;
  const int expected_block_nx3 = (dimension == 3) ? 32 : 1;
  if (!exact_topology || pmy_mesh_->mesh_indcs.nx1 != expected_global_nx1 ||
      pmy_mesh_->mesh_indcs.nx2 != expected_global_nx2 ||
      pmy_mesh_->mesh_indcs.nx3 != expected_global_nx3 ||
      pmy_mesh_->mb_indcs.nx1 != expected_block_nx1 ||
      pmy_mesh_->mb_indcs.nx2 != expected_block_nx2 ||
      pmy_mesh_->mb_indcs.nx3 != expected_block_nx3) {
    Q029Fatal("q029_hall_bell_linear requires its exact global and meshblock "
              "launch-preparation geometry");
  }
  const Real expected_x1max =
      (dimension == 1) ? 1.0 : ((dimension == 2) ? std::sqrt(5.0) :
                                                       std::sqrt(21.0));
  const Real expected_x2max =
      (dimension == 1) ? 1.0 : ((dimension == 2) ? std::sqrt(1.25) :
                                                       std::sqrt(5.25));
  const Real expected_x3max = (dimension == 3) ? std::sqrt(1.3125) : 1.0;
  Q029RequireClose("global x1min", pmy_mesh_->mesh_size.x1min, 0.0);
  Q029RequireClose("global x1max", pmy_mesh_->mesh_size.x1max, expected_x1max);
  Q029RequireClose("global x2min", pmy_mesh_->mesh_size.x2min, 0.0);
  Q029RequireClose("global x2max", pmy_mesh_->mesh_size.x2max, expected_x2max);
  Q029RequireClose("global x3min", pmy_mesh_->mesh_size.x3min, 0.0);
  Q029RequireClose("global x3max", pmy_mesh_->mesh_size.x3max, expected_x3max);
  Q029RequireString(pin, "particles", "pic_physical_mode", "extended_mhd_pic");
  Q029RequireString(pin, "particles", "pic_cr_initial_state", "velocity");
  Q029RequireString(
      pin, "particles", "pic_cr_hall_mode", "current_to_ct_experimental");
  Q029RequireBoolean(pin, "particles", "pic_enable_2d3v", true);
  Q029RequireString(pin, block, "campaign_id", "Q029-HALL-BELL-LINEAR");
  Q029RequireString(
      pin, block, "deck_role",
      (dimension == 1) ?
          "source_local_runnable_thin_2d3v_launch_preparation_not_authorized" :
          "source_local_runnable_launch_preparation_not_authorized");
  Q029RequireString(pin, block, "epsilon_grid", "0.1,0.2,0.4,0.6,0.8");
  Q029RequireString(pin, block, "chi_h_grid", "0.25,0.5,1.0");
  Q029RequireString(
      pin, block, "alpha_h_parameter", "particles/couple_j_to_efield_coeff");
  Q029RequireString(pin, block, "chi_h_definition", "alpha_h*j_cr/(u_a*b_g)");
  Q029RequireString(
      pin, block, "chi_h_grid_semantics",
      "positive_launch_preparation_only_no_hall_bell_qualification");
  Q029RequireString(
      pin, block, "initial_eigenmode",
      "q023_section52_seed_carrier_only_not_hall_dispersion_oracle");
  Q029RequireString(pin, block, "qualification_effect", "none");
  Q029RequireString(pin, block, "linear_hall_bell", "open_not_claimed");
  Q029RequireString(pin, block, "nonlinear_hall_bell", "open_not_claimed");
  Q029RequireString(pin, block, "timestep", "open_clean_candidate_timestep_freeze");

  const Real epsilon_default = pin->GetReal(block, "epsilon_default");
  const Real epsilon = pin->GetOrAddReal(block, "epsilon", epsilon_default);
  const Real chi_h_default = pin->GetReal(block, "chi_h_default");
  const Real chi_h = pin->GetOrAddReal(block, "chi_h", chi_h_default);
  const Real rho = pin->GetReal(block, "rho");
  const Real pressure = pin->GetReal(block, "pressure");
  const Real amplitude = pin->GetReal(block, "amplitude");
  const Real b_g = pin->GetReal(block, "b_g");
  const Real u_a = pin->GetReal(block, "u_a");
  const Real wavelength = pin->GetReal(block, "wavelength");
  const Real k0 = pin->GetReal(block, "k0");
  const Real omega = pin->GetReal(block, "omega");
  const Real c_over_v_cr = pin->GetReal(block, "c_over_v_cr");
  Q029RequireFinite("epsilon_default", epsilon_default);
  Q029RequireFinite("epsilon", epsilon);
  Q029RequireFinite("chi_h_default", chi_h_default);
  Q029RequireFinite("chi_h", chi_h);
  Q029RequireFinite("rho", rho);
  Q029RequireFinite("pressure", pressure);
  Q029RequireFinite("amplitude", amplitude);
  Q029RequireFinite("b_g", b_g);
  Q029RequireFinite("u_a", u_a);
  Q029RequireFinite("wavelength", wavelength);
  Q029RequireFinite("k0", k0);
  Q029RequireFinite("omega", omega);
  Q029RequireFinite("c_over_v_cr", c_over_v_cr);
  if (!(epsilon > 0.0 && epsilon < 1.0) || !IsPositivePreparedChiH(chi_h) ||
      pressure <= 0.0 || amplitude <= 0.0 || amplitude > 1.0e-2) {
    Q029Fatal("q029_hall_bell_linear epsilon/chi_H/pressure/amplitude contract "
              "is invalid");
  }
  Q029RequireClose("rho", rho, 1.0);
  Q029RequireClose("b_g", b_g, 1.0);
  Q029RequireClose("u_a", u_a, b_g/std::sqrt(rho));
  Q029RequireClose("wavelength", wavelength, 1.0);
  Q029RequireClose("k0", k0, 2.0*M_PI/wavelength);
  Q029RequireClose("omega", omega, 1.0e-6*k0*u_a);
  Q029RequireClose("c_over_v_cr", c_over_v_cr, 1.0e3);

  const Basis basis = ModeBasis(dimension);
  const Real cr_vx = pin->GetReal("particles", "cr_vx0");
  const Real cr_vy = pin->GetReal("particles", "cr_vy0");
  const Real cr_vz = pin->GetReal("particles", "cr_vz0");
  const Real v_cr = std::sqrt(cr_vx*cr_vx + cr_vy*cr_vy + cr_vz*cr_vz);
  const Real light_speed = pin->GetReal("particles", "pic_cr_light_speed");
  const Real ppc = pin->GetReal("particles", "ppc");
  const Real qscale = pin->GetReal("particles", "deposit_qscale");
  const Real alpha_h = pin->GetReal("particles", "couple_j_to_efield_coeff");
  const Real mass = pin->GetReal("species0", "mass");
  const Real charge = pin->GetReal("species0", "charge");
  const Real j_cr = CurrentDensity(ppc, qscale, charge, v_cr);
  Q029RequireFinite("cr_vx0", cr_vx);
  Q029RequireFinite("cr_vy0", cr_vy);
  Q029RequireFinite("cr_vz0", cr_vz);
  Q029RequireFinite("v_cr", v_cr);
  Q029RequireFinite("pic_cr_light_speed", light_speed);
  Q029RequireFinite("ppc", ppc);
  Q029RequireFinite("deposit_qscale", qscale);
  Q029RequireFinite("couple_j_to_efield_coeff", alpha_h);
  Q029RequireFinite("species mass", mass);
  Q029RequireFinite("species charge", charge);
  Q029RequireFinite("j_CR", j_cr);
  Q029RequireClose("epsilon", epsilon, u_a/v_cr);
  Q029RequireClose("cr_vx0", cr_vx, v_cr*basis.parallel.x1);
  Q029RequireClose("cr_vy0", cr_vy, v_cr*basis.parallel.x2);
  Q029RequireClose("cr_vz0", cr_vz, v_cr*basis.parallel.x3);
  Q029RequireClose("pic_cr_light_speed", light_speed, c_over_v_cr*v_cr);
  Q029RequireClose("species q/m", charge/mass, omega/b_g);
  Q029RequireClose("j_CR", j_cr, 2.0*b_g*light_speed*k0);
  Q029RequireClose("chi_H", chi_h, ChiH(alpha_h, j_cr, u_a, b_g));

  if (restart) return;

  const ModeParameters parameters = {
    dimension, epsilon, amplitude, rho, pressure, b_g, u_a, k0
  };
  auto &indcs = pmy_mesh_->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const int nmb = pmbp->nmb_thispack;
  auto &size = pmbp->pmb->mb_size;
  auto &w0 = pmbp->pmhd->w0;
  auto &u0 = pmbp->pmhd->u0;
  auto &b0 = pmbp->pmhd->b0;
  auto &bcc0 = pmbp->pmhd->bcc0;

  const int ncells1 = indcs.nx1 + 2*indcs.ng;
  const int ncells2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 2;
  const int ncells3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 2;
  DvceArray4D<Real> a1;
  DvceArray4D<Real> a2;
  DvceArray4D<Real> a3;
  Kokkos::realloc(a1, nmb, ncells3, ncells2, ncells1);
  Kokkos::realloc(a2, nmb, ncells3, ncells2, ncells1);
  Kokkos::realloc(a3, nmb, ncells3, ncells2, ncells1);

  par_for("pgen_q029_hall_bell_vector_potential", DevExeSpace(),
          0, nmb - 1, ks, ke + 1, js, je + 1, is, ie + 1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const int nx1 = indcs.nx1;
    const int nx2 = indcs.nx2;
    const int nx3 = indcs.nx3;
    const Real x1v = CellCenterX(i - is, nx1, size.d_view(m).x1min,
                                 size.d_view(m).x1max);
    const Real x1f = LeftEdgeX(i - is, nx1, size.d_view(m).x1min,
                               size.d_view(m).x1max);
    const Real x2v = CellCenterX(j - js, nx2, size.d_view(m).x2min,
                                 size.d_view(m).x2max);
    const Real x2f = LeftEdgeX(j - js, nx2, size.d_view(m).x2min,
                               size.d_view(m).x2max);
    const Real x3v = CellCenterX(k - ks, nx3, size.d_view(m).x3min,
                                 size.d_view(m).x3max);
    const Real x3f = LeftEdgeX(k - ks, nx3, size.d_view(m).x3min,
                               size.d_view(m).x3max);
    a1(m, k, j, i) = VectorPotentialAt(parameters, {x1v, x2f, x3f}).x1;
    a2(m, k, j, i) = VectorPotentialAt(parameters, {x1f, x2v, x3f}).x2;
    a3(m, k, j, i) = VectorPotentialAt(parameters, {x1f, x2f, x3v}).x3;
  });

  par_for("pgen_q029_hall_bell_ct_field", DevExeSpace(),
          0, nmb - 1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real dx1 = size.d_view(m).dx1;
    const Real dx2 = size.d_view(m).dx2;
    const Real dx3 = size.d_view(m).dx3;
    b0.x1f(m, k, j, i) = (a3(m, k, j + 1, i) - a3(m, k, j, i))/dx2 -
                          (a2(m, k + 1, j, i) - a2(m, k, j, i))/dx3;
    b0.x2f(m, k, j, i) = (a1(m, k + 1, j, i) - a1(m, k, j, i))/dx3 -
                          (a3(m, k, j, i + 1) - a3(m, k, j, i))/dx1;
    b0.x3f(m, k, j, i) = (a2(m, k, j, i + 1) - a2(m, k, j, i))/dx1 -
                          (a1(m, k, j + 1, i) - a1(m, k, j, i))/dx2;
    if (i == ie) {
      b0.x1f(m, k, j, i + 1) =
          (a3(m, k, j + 1, i + 1) - a3(m, k, j, i + 1))/dx2 -
          (a2(m, k + 1, j, i + 1) - a2(m, k, j, i + 1))/dx3;
    }
    if (j == je) {
      b0.x2f(m, k, j + 1, i) =
          (a1(m, k + 1, j + 1, i) - a1(m, k, j + 1, i))/dx3 -
          (a3(m, k, j + 1, i + 1) - a3(m, k, j + 1, i))/dx1;
    }
    if (k == ke) {
      b0.x3f(m, k + 1, j, i) =
          (a2(m, k + 1, j, i + 1) - a2(m, k + 1, j, i))/dx1 -
          (a1(m, k + 1, j + 1, i) - a1(m, k + 1, j, i))/dx2;
    }
  });

  par_for("pgen_q029_hall_bell_primitives", DevExeSpace(),
          0, nmb - 1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real x1 = CellCenterX(i - is, indcs.nx1, size.d_view(m).x1min,
                                size.d_view(m).x1max);
    const Real x2 = CellCenterX(j - js, indcs.nx2, size.d_view(m).x2min,
                                size.d_view(m).x2max);
    const Real x3 = CellCenterX(k - ks, indcs.nx3, size.d_view(m).x3min,
                                size.d_view(m).x3max);
    const Real phase = parameters.k0*
        (basis.parallel.x1*x1 + basis.parallel.x2*x2 + basis.parallel.x3*x3);
    const auto sample = EigenmodeAtPhase(parameters, phase);
    w0(m, IDN, k, j, i) = parameters.rho;
    w0(m, IVX, k, j, i) = sample.velocity.x1;
    w0(m, IVY, k, j, i) = sample.velocity.x2;
    w0(m, IVZ, k, j, i) = sample.velocity.x3;
    w0(m, IEN, k, j, i) = parameters.pressure;
    bcc0(m, IBX, k, j, i) =
        0.5*(b0.x1f(m, k, j, i) + b0.x1f(m, k, j, i + 1));
    bcc0(m, IBY, k, j, i) =
        0.5*(b0.x2f(m, k, j, i) + b0.x2f(m, k, j + 1, i));
    bcc0(m, IBZ, k, j, i) =
        0.5*(b0.x3f(m, k, j, i) + b0.x3f(m, k + 1, j, i));
  });
  pmbp->pmhd->peos->PrimToCons(w0, bcc0, u0, is, ie, js, je, ks, ke);
}
#endif
