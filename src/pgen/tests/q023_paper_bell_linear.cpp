//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file q023_paper_bell_linear.cpp
//! \brief Source-local Q-023 Sun and Bai Section 5.2 Bell eigenmode preparation.

#include <algorithm>
#include <cmath>

#if !defined(Q023_PAPER_BELL_LINEAR_HOST_CONTRACT)
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

#if defined(Q023_PAPER_BELL_LINEAR_HOST_CONTRACT)
#define Q023_INLINE inline
#else
#define Q023_INLINE KOKKOS_INLINE_FUNCTION
#endif

namespace q023_paper_bell_linear {

struct Vector3 {
  double x1;
  double x2;
  double x3;
};

struct Basis {
  Vector3 parallel;
  Vector3 transverse_a;
  Vector3 transverse_b;
};

struct ModeParameters {
  int dimension;
  double epsilon;
  double amplitude;
  double rho;
  double pressure;
  double b_g;
  double u_a;
  double k0;
};

struct ModeSample {
  Vector3 magnetic;
  Vector3 velocity;
};

Q023_INLINE Vector3 Add(const Vector3 a, const Vector3 b) {
  return {a.x1 + b.x1, a.x2 + b.x2, a.x3 + b.x3};
}

Q023_INLINE Vector3 Scale(const Vector3 value, const double factor) {
  return {factor*value.x1, factor*value.x2, factor*value.x3};
}

Q023_INLINE double Dot(const Vector3 a, const Vector3 b) {
  return a.x1*b.x1 + a.x2*b.x2 + a.x3*b.x3;
}

Q023_INLINE Vector3 Cross(const Vector3 a, const Vector3 b) {
  return {
    a.x2*b.x3 - a.x3*b.x2,
    a.x3*b.x1 - a.x1*b.x3,
    a.x1*b.x2 - a.x2*b.x1
  };
}

Q023_INLINE Vector3 Normalize(const Vector3 value) {
  const double inv_norm = 1.0/std::sqrt(Dot(value, value));
  return Scale(value, inv_norm);
}

Q023_INLINE Basis ModeBasis(const int dimension) {
  const Vector3 parallel = Normalize({
    1.0,
    (dimension >= 2) ? 2.0 : 0.0,
    (dimension >= 3) ? 4.0 : 0.0
  });
  const Vector3 transverse_a = (dimension == 1) ?
      Vector3{0.0, 1.0, 0.0} :
      Normalize(Vector3{-parallel.x2, parallel.x1, 0.0});
  return {parallel, transverse_a, Cross(parallel, transverse_a)};
}

Q023_INLINE double NormalizedGrowthRate(const double epsilon) {
  return std::sqrt(1.0 - epsilon*epsilon);
}

Q023_INLINE ModeSample EigenmodeAtPhase(const ModeParameters parameters,
                                        const double phase) {
  const Basis basis = ModeBasis(parameters.dimension);
  const double cs = std::cos(phase);
  const double sn = std::sin(phase);
  const double growth = NormalizedGrowthRate(parameters.epsilon);
  const Vector3 magnetic = Add(
      Scale(basis.parallel, parameters.b_g),
      Add(Scale(basis.transverse_a, parameters.amplitude*cs),
          Scale(basis.transverse_b, -parameters.amplitude*sn)));
  const double velocity_scale = parameters.amplitude/std::sqrt(parameters.rho);
  const Vector3 velocity = Add(
      Scale(basis.transverse_a,
            velocity_scale*(-parameters.epsilon*cs + growth*sn)),
      Scale(basis.transverse_b,
            velocity_scale*(growth*cs + parameters.epsilon*sn)));
  return {magnetic, velocity};
}

Q023_INLINE Vector3 VectorPotentialAt(const ModeParameters parameters,
                                     const Vector3 position) {
  const Basis basis = ModeBasis(parameters.dimension);
  const double phase = parameters.k0*Dot(basis.parallel, position);
  const double cs = std::cos(phase);
  const double sn = std::sin(phase);
  const Vector3 guide = Scale(Cross(basis.parallel, position),
                              0.5*parameters.b_g);
  const Vector3 perturbation = Add(
      Scale(basis.transverse_a, parameters.amplitude*cs/parameters.k0),
      Scale(basis.transverse_b, -parameters.amplitude*sn/parameters.k0));
  return Add(guide, perturbation);
}

}  // namespace q023_paper_bell_linear

#if !defined(Q023_PAPER_BELL_LINEAR_HOST_CONTRACT)
namespace {

[[noreturn]] void Q023Fatal(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << message << std::endl;
  std::exit(EXIT_FAILURE);
}

void Q023RequireClose(const std::string &label, const Real measured,
                      const Real expected) {
  const Real scale = std::max(static_cast<Real>(1.0), std::abs(expected));
  if (std::abs(measured - expected) > static_cast<Real>(1.0e-12)*scale) {
    Q023Fatal(label + " does not match the Q-023 Section 5.2 contract");
  }
}

void Q023RequireString(ParameterInput *pin, const std::string &block,
                       const std::string &name, const std::string &expected) {
  if (pin->GetString(block, name).compare(expected) != 0) {
    Q023Fatal("<" + block + ">/" + name +
              " does not match the Q-023 Section 5.2 contract");
  }
}

}  // namespace

void ProblemGenerator::Q023PaperBellLinear(ParameterInput *pin, const bool restart) {
  using q023_paper_bell_linear::Basis;
  using q023_paper_bell_linear::EigenmodeAtPhase;
  using q023_paper_bell_linear::ModeBasis;
  using q023_paper_bell_linear::ModeParameters;
  using q023_paper_bell_linear::Vector3;
  using q023_paper_bell_linear::VectorPotentialAt;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->pmhd == nullptr || pmbp->ppart == nullptr) {
    Q023Fatal("q023_paper_bell_linear requires MHD and particles");
  }
  if (!pmbp->pmhd->peos->eos_data.is_ideal) {
    Q023Fatal("q023_paper_bell_linear requires ideal-MHD EOS");
  }
  if (pmy_mesh_->multilevel) {
    Q023Fatal("q023_paper_bell_linear source-local preparation rejects AMR/SMR");
  }

  const std::string block = "q023_paper_bell_linear";
  const int dimension = pin->GetInteger(block, "dimension");
  const int mesh_dimension = pmy_mesh_->one_d ? 1 : (pmy_mesh_->two_d ? 2 : 3);
  const bool thin_2d3v_1d_carrier =
      dimension == 1 && mesh_dimension == 2 &&
      pmy_mesh_->mb_indcs.nx2 == 4 && pmy_mesh_->mb_indcs.nx3 == 1;
  const bool native_dimension =
      dimension >= 2 && dimension <= 3 && dimension == mesh_dimension;
  if (!(thin_2d3v_1d_carrier || native_dimension)) {
    Q023Fatal("<q023_paper_bell_linear>/dimension does not match the required "
              "native mesh or exact thin 2D3V carrier");
  }
  Q023RequireString(pin, block, "campaign_id", "Q023-PAPER-BELL-LINEAR");
  Q023RequireString(
      pin, block, "deck_role",
      (dimension == 1) ?
          "source_local_runnable_thin_2d3v_carrier_preparation_not_authorized" :
          "source_local_runnable_preparation_not_authorized");
  Q023RequireString(pin, block, "epsilon_grid", "0.1,0.2,0.4,0.6,0.8");
  Q023RequireString(pin, block, "initial_eigenmode",
                    "section52_right_polarized_eigenmode");
  Q023RequireString(pin, block, "timestep", "open_clean_candidate_timestep_freeze");

  const Real epsilon_default = pin->GetReal(block, "epsilon_default");
  const Real epsilon = pin->GetOrAddReal(block, "epsilon", epsilon_default);
  const Real rho = pin->GetReal(block, "rho");
  const Real pressure = pin->GetReal(block, "pressure");
  const Real amplitude = pin->GetReal(block, "amplitude");
  const Real b_g = pin->GetReal(block, "b_g");
  const Real u_a = pin->GetReal(block, "u_a");
  const Real wavelength = pin->GetReal(block, "wavelength");
  const Real k0 = pin->GetReal(block, "k0");
  const Real omega = pin->GetReal(block, "omega");
  const Real c_over_v_cr = pin->GetReal(block, "c_over_v_cr");
  if (!(epsilon > 0.0 && epsilon < 1.0) || pressure <= 0.0 ||
      amplitude <= 0.0 || amplitude > 1.0e-2) {
    Q023Fatal("q023_paper_bell_linear epsilon/pressure/amplitude contract is invalid");
  }
  Q023RequireClose("rho", rho, 1.0);
  Q023RequireClose("b_g", b_g, 1.0);
  Q023RequireClose("u_a", u_a, b_g/std::sqrt(rho));
  Q023RequireClose("wavelength", wavelength, 1.0);
  Q023RequireClose("k0", k0, 2.0*M_PI/wavelength);
  Q023RequireClose("omega", omega, 1.0e-6*k0*u_a);
  Q023RequireClose("c_over_v_cr", c_over_v_cr, 1.0e3);

  const Basis basis = ModeBasis(dimension);
  const Real cr_vx = pin->GetReal("particles", "cr_vx0");
  const Real cr_vy = pin->GetReal("particles", "cr_vy0");
  const Real cr_vz = pin->GetReal("particles", "cr_vz0");
  const Real v_cr = std::sqrt(cr_vx*cr_vx + cr_vy*cr_vy + cr_vz*cr_vz);
  const Real light_speed = pin->GetReal("particles", "pic_cr_light_speed");
  const Real ppc = pin->GetReal("particles", "ppc");
  const Real qscale = pin->GetReal("particles", "deposit_qscale");
  const Real mass = pin->GetReal("species0", "mass");
  const Real charge = pin->GetReal("species0", "charge");
  Q023RequireClose("epsilon", epsilon, u_a/v_cr);
  Q023RequireClose("cr_vx0", cr_vx, v_cr*basis.parallel.x1);
  Q023RequireClose("cr_vy0", cr_vy, v_cr*basis.parallel.x2);
  Q023RequireClose("cr_vz0", cr_vz, v_cr*basis.parallel.x3);
  Q023RequireClose("pic_cr_light_speed", light_speed, c_over_v_cr*v_cr);
  Q023RequireClose("species q/m", charge/mass, omega/b_g);
  Q023RequireClose("j_CR", ppc*qscale*charge*v_cr, 2.0*b_g*light_speed*k0);

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

  par_for("pgen_q023_bell_vector_potential", DevExeSpace(),
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

  par_for("pgen_q023_bell_ct_field", DevExeSpace(),
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

  par_for("pgen_q023_bell_primitives", DevExeSpace(),
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

#undef Q023_INLINE
