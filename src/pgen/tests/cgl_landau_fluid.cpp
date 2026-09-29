//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_landau_fluid.cpp
//! \brief Quantitative built-in tests for CGL Landau-fluid heat-flux transport.

#include <sys/stat.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "eos/cgl_physics.hpp"
#include "diffusion/cgl_landau_fluid.hpp"
#include "mhd/mhd.hpp"
#include "parameter_input.hpp"
#include "pgen/pgen.hpp"
#include "outputs/outputs.hpp"

namespace {

enum class TestMode {
  parallel_decay,
  perp_decay,
  collision_relaxation,
  grad_b,
  flux_limiter,
  limiter_heat_flux_suppression,
  limiter_stress,
  field_aligned_wave,
  paper_oblique_wave,
  paper_eigen_wave,
  rotated_decay,
  field_reversal,
  density_contact,
  timestep_refresh,
  hotspot,
  low_field
};

struct Projection {
  Real mean = 0.0;
  Real sin_amp = 0.0;
  Real cos_amp = 0.0;
};

using Complex = std::complex<Real>;

template <typename ViewType>
auto HostCopy(const ViewType &view) {
  return Kokkos::create_mirror_view_and_copy(HostMemSpace(), view);
}

[[noreturn]] void Fail(const std::string &msg) {
  std::cout << "CGL LF quantitative test failed: " << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

void Require(bool condition, const std::string &msg) {
  if (!condition) {
    Fail(msg);
  }
}

void EnsureDirectory(const std::string &dir) {
  std::string partial;
  for (const char c : dir) {
    partial.push_back(c);
    if (c == '/' && partial.size() > 1) {
      mkdir(partial.c_str(), 0775);
    }
  }
  mkdir(dir.c_str(), 0775);
}

std::ofstream OpenValidationDataFile(ParameterInput *pin, const std::string &suffix) {
  std::ofstream out;
  const bool enabled = pin->GetOrAddBoolean("problem", "validation_output", false);
  if (!enabled || global_variable::my_rank != 0) {
    return out;
  }
  const std::string dir =
      pin->GetOrAddString("problem", "validation_output_dir", "docs/figures/data");
  EnsureDirectory(dir);
  std::ostringstream fname;
  fname << dir << "/" << pin->GetString("job", "basename") << "." << suffix << ".csv";
  out.open(fname.str());
  if (!out) {
    Fail("could not open validation output file " + fname.str());
  }
  out << std::setprecision(17);
  return out;
}

Real RelativeError(const Real got, const Real expected) {
  const Real scale = std::max(std::abs(expected), static_cast<Real>(1.0e-30));
  return std::abs(got - expected)/scale;
}

Real ComplexRelativeError(const Complex got, const Complex expected, const Real scale) {
  return std::abs(got - expected)/std::max(scale, static_cast<Real>(1.0e-30));
}

void RequireRelative(const std::string &label, const Real got, const Real expected,
                     const Real rel_tol) {
  const Real scale = std::max(std::abs(expected), static_cast<Real>(1.0e-30));
  const Real rel_err = std::abs(got - expected)/scale;
  if (!std::isfinite(got) || !std::isfinite(expected) || rel_err > rel_tol) {
    std::cout << label << " got=" << got << " expected=" << expected
              << " rel_err=" << rel_err << " rel_tol=" << rel_tol << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

TestMode ParseMode(ParameterInput *pin) {
  const std::string mode = pin->GetOrAddString("problem", "test_mode", "parallel_decay");
  if (mode == "parallel_decay") return TestMode::parallel_decay;
  if (mode == "perp_decay") return TestMode::perp_decay;
  if (mode == "collision_relaxation") return TestMode::collision_relaxation;
  if (mode == "grad_b") return TestMode::grad_b;
  if (mode == "flux_limiter") return TestMode::flux_limiter;
  if (mode == "limiter_heat_flux_suppression") {
    return TestMode::limiter_heat_flux_suppression;
  }
  if (mode == "limiter_stress") return TestMode::limiter_stress;
  if (mode == "field_aligned_wave") return TestMode::field_aligned_wave;
  if (mode == "paper_oblique_wave") return TestMode::paper_oblique_wave;
  if (mode == "paper_eigen_wave") return TestMode::paper_eigen_wave;
  if (mode == "rotated_decay") return TestMode::rotated_decay;
  if (mode == "timestep_refresh") return TestMode::timestep_refresh;
  if (mode == "density_contact") return TestMode::density_contact;
  if (mode == "field_reversal") return TestMode::field_reversal;
  if (mode == "hotspot") return TestMode::hotspot;
  if (mode == "low_field") return TestMode::low_field;
  Fail("<problem>/test_mode must be parallel_decay, perp_decay, "
       "collision_relaxation, grad_b, flux_limiter, "
       "limiter_heat_flux_suppression, limiter_stress, "
       "field_aligned_wave, paper_oblique_wave, paper_eigen_wave, rotated_decay, "
       "density_contact, timestep_refresh, field_reversal, hotspot, or low_field");
}

const char *ModeName(const TestMode mode) {
  switch (mode) {
    case TestMode::parallel_decay: return "parallel_decay";
    case TestMode::perp_decay: return "perp_decay";
    case TestMode::collision_relaxation: return "collision_relaxation";
    case TestMode::grad_b: return "grad_b";
    case TestMode::flux_limiter: return "flux_limiter";
    case TestMode::limiter_heat_flux_suppression: return "limiter_heat_flux_suppression";
    case TestMode::limiter_stress: return "limiter_stress";
    case TestMode::field_aligned_wave: return "field_aligned_wave";
    case TestMode::paper_oblique_wave: return "paper_oblique_wave";
    case TestMode::paper_eigen_wave: return "paper_eigen_wave";
    case TestMode::rotated_decay: return "rotated_decay";
    case TestMode::timestep_refresh: return "timestep_refresh";
    case TestMode::density_contact: return "density_contact";
    case TestMode::field_reversal: return "field_reversal";
    case TestMode::hotspot: return "hotspot";
    case TestMode::low_field: return "low_field";
  }
  return "unknown";
}

void RequireOneDimensionalSingleBlock(Mesh *pm) {
  Require(pm->pmb_pack->nmb_thispack == 1,
          "this unit test currently expects one MeshBlock on one rank");
  Require(pm->mb_indcs.nx2 == 1 && pm->mb_indcs.nx3 == 1,
          "this unit test currently expects a 1D mesh");
}

void RequireOneDimensionalMesh(Mesh *pm) {
  Require(pm->mb_indcs.nx2 == 1 && pm->mb_indcs.nx3 == 1,
          "this unit test currently expects a 1D mesh");
}

void RequireSingleBlock(Mesh *pm) {
  Require(pm->pmb_pack->nmb_thispack == 1,
          "this unit test currently expects one MeshBlock on one rank");
}

Real XCenter(Mesh *pm, const int q) {
  const Real xmin = pm->mesh_size.x1min;
  const Real xmax = pm->mesh_size.x1max;
  return CellCenterX(q, pm->mb_indcs.nx1, xmin, xmax);
}

Real Wavenumber(ParameterInput *pin, Mesh *pm) {
  const int mode_number = pin->GetOrAddInteger("problem", "mode_number", 1);
  const Real length = pm->mesh_size.x1max - pm->mesh_size.x1min;
  return 2.0*M_PI*static_cast<Real>(mode_number)/length;
}

struct RotatedWave {
  Real kx = 0.0;
  Real ky = 0.0;
  Real kz = 0.0;
};

RotatedWave RotatedWavenumber(ParameterInput *pin, Mesh *pm) {
  const int mode_number = pin->GetOrAddInteger("problem", "mode_number", 1);
  const std::string axis = pin->GetOrAddString("problem", "rotated_axis", "x");
  const Real multiple = 2.0*M_PI*static_cast<Real>(mode_number);
  const Real kx = multiple/(pm->mesh_size.x1max - pm->mesh_size.x1min);
  const Real ky = multiple/(pm->mesh_size.x2max - pm->mesh_size.x2min);
  const Real kz = multiple/(pm->mesh_size.x3max - pm->mesh_size.x3min);
  if (axis == "x") return {kx, 0.0, 0.0};
  if (axis == "y") return {0.0, ky, 0.0};
  if (axis == "z") return {0.0, 0.0, kz};
  if (axis == "oblique") return {kx, ky, kz};
  Fail("<problem>/rotated_axis must be x, y, z, or oblique");
}

template <typename HostView>
Projection ProjectTemperature(const HostView &w, ParameterInput *pin, Mesh *pm,
                              const int pidx) {
  RequireOneDimensionalSingleBlock(pm);
  const int is = pm->mb_indcs.is;
  const int js = pm->mb_indcs.js;
  const int ks = pm->mb_indcs.ks;
  const int nx1 = pm->mb_indcs.nx1;
  const Real k_wave = Wavenumber(pin, pm);
  const Real xmin = pm->mesh_size.x1min;

  Projection p;
  for (int q = 0; q < nx1; ++q) {
    const int i = is + q;
    p.mean += w(0,pidx,ks,js,i)/w(0,IDN,ks,js,i);
  }
  p.mean /= static_cast<Real>(nx1);

  for (int q = 0; q < nx1; ++q) {
    const int i = is + q;
    const Real x = XCenter(pm, q);
    const Real phase = k_wave*(x - xmin);
    const Real value = w(0,pidx,ks,js,i)/w(0,IDN,ks,js,i) - p.mean;
    p.sin_amp += value*std::sin(phase);
    p.cos_amp += value*std::cos(phase);
  }
  p.sin_amp *= 2.0/static_cast<Real>(nx1);
  p.cos_amp *= 2.0/static_cast<Real>(nx1);
  return p;
}

template <typename HostView>
Projection ProjectRotatedTemperature(const HostView &w, ParameterInput *pin, Mesh *pm,
                                     const int pidx) {
  RequireSingleBlock(pm);
  const int is = pm->mb_indcs.is;
  const int js = pm->mb_indcs.js;
  const int ks = pm->mb_indcs.ks;
  const int nx1 = pm->mb_indcs.nx1;
  const int nx2 = pm->mb_indcs.nx2;
  const int nx3 = pm->mb_indcs.nx3;
  const RotatedWave wave = RotatedWavenumber(pin, pm);
  const Real xmin = pm->mesh_size.x1min;
  const Real ymin = pm->mesh_size.x2min;
  const Real zmin = pm->mesh_size.x3min;
  const Real ncells = static_cast<Real>(nx1*nx2*nx3);
  Projection p;
  for (int qk = 0; qk < nx3; ++qk) {
    for (int qj = 0; qj < nx2; ++qj) {
      for (int qi = 0; qi < nx1; ++qi) {
        const Real value = w(0,pidx,ks + qk,js + qj,is + qi)
                         / w(0,IDN,ks + qk,js + qj,is + qi);
        p.mean += value;
      }
    }
  }
  p.mean /= ncells;
  for (int qk = 0; qk < nx3; ++qk) {
    const Real z = CellCenterX(qk, nx3, zmin, pm->mesh_size.x3max);
    for (int qj = 0; qj < nx2; ++qj) {
      const Real y = CellCenterX(qj, nx2, ymin, pm->mesh_size.x2max);
      for (int qi = 0; qi < nx1; ++qi) {
        const Real value = w(0,pidx,ks + qk,js + qj,is + qi)
                         / w(0,IDN,ks + qk,js + qj,is + qi) - p.mean;
        const Real x = CellCenterX(qi, nx1, xmin, pm->mesh_size.x1max);
        const Real phase = wave.kx*(x - xmin) + wave.ky*(y - ymin)
                         + wave.kz*(z - zmin);
        p.sin_amp += value*std::sin(phase);
        p.cos_amp += value*std::cos(phase);
      }
    }
  }
  p.sin_amp *= 2.0/ncells;
  p.cos_amp *= 2.0/ncells;
  return p;
}

template <typename HostView>
Projection ProjectPrimitive(const HostView &w, ParameterInput *pin, Mesh *pm,
                            const int idx) {
  RequireOneDimensionalMesh(pm);
  const int is = pm->mb_indcs.is;
  const int js = pm->mb_indcs.js;
  const int ks = pm->mb_indcs.ks;
  const int nx1 = pm->mb_indcs.nx1;
  const Real k_wave = Wavenumber(pin, pm);
  const Real xmin = pm->mesh_size.x1min;
  const Real length = pm->mesh_size.x1max - pm->mesh_size.x1min;
  auto &size = pm->pmb_pack->pmb->mb_size;
  size.template sync<HostMemSpace>();

  Projection p;
  for (int m = 0; m < pm->pmb_pack->nmb_thispack; ++m) {
    for (int q = 0; q < nx1; ++q) {
      p.mean += w(m,idx,ks,js,is + q)*size.h_view(m).dx1;
    }
  }
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &p.mean, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
  p.mean /= length;

  for (int m = 0; m < pm->pmb_pack->nmb_thispack; ++m) {
    for (int q = 0; q < nx1; ++q) {
      const Real x = CellCenterX(q, nx1, size.h_view(m).x1min,
                                 size.h_view(m).x1max);
      const Real phase = k_wave*(x - xmin);
      const Real value = w(m,idx,ks,js,is + q) - p.mean;
      p.sin_amp += value*std::sin(phase)*size.h_view(m).dx1;
      p.cos_amp += value*std::cos(phase)*size.h_view(m).dx1;
    }
  }
#if MPI_PARALLEL_ENABLED
  Real amplitudes[2] = {p.sin_amp, p.cos_amp};
  MPI_Allreduce(MPI_IN_PLACE, amplitudes, 2, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  p.sin_amp = amplitudes[0];
  p.cos_amp = amplitudes[1];
#endif
  p.sin_amp *= 2.0/length;
  p.cos_amp *= 2.0/length;
  return p;
}

template <typename HostView>
Projection ProjectCellField(const HostView &bcc, ParameterInput *pin, Mesh *pm,
                            const int idx) {
  RequireOneDimensionalMesh(pm);
  const int is = pm->mb_indcs.is;
  const int js = pm->mb_indcs.js;
  const int ks = pm->mb_indcs.ks;
  const int nx1 = pm->mb_indcs.nx1;
  const Real k_wave = Wavenumber(pin, pm);
  const Real xmin = pm->mesh_size.x1min;
  const Real length = pm->mesh_size.x1max - pm->mesh_size.x1min;
  auto &size = pm->pmb_pack->pmb->mb_size;
  size.template sync<HostMemSpace>();

  Projection p;
  for (int m = 0; m < pm->pmb_pack->nmb_thispack; ++m) {
    for (int q = 0; q < nx1; ++q) {
      p.mean += bcc(m,idx,ks,js,is + q)*size.h_view(m).dx1;
    }
  }
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &p.mean, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
  p.mean /= length;

  for (int m = 0; m < pm->pmb_pack->nmb_thispack; ++m) {
    for (int q = 0; q < nx1; ++q) {
      const Real x = CellCenterX(q, nx1, size.h_view(m).x1min,
                                 size.h_view(m).x1max);
      const Real phase = k_wave*(x - xmin);
      const Real value = bcc(m,idx,ks,js,is + q) - p.mean;
      p.sin_amp += value*std::sin(phase)*size.h_view(m).dx1;
      p.cos_amp += value*std::cos(phase)*size.h_view(m).dx1;
    }
  }
#if MPI_PARALLEL_ENABLED
  Real amplitudes[2] = {p.sin_amp, p.cos_amp};
  MPI_Allreduce(MPI_IN_PLACE, amplitudes, 2, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  p.sin_amp = amplitudes[0];
  p.cos_amp = amplitudes[1];
#endif
  p.sin_amp *= 2.0/length;
  p.cos_amp *= 2.0/length;
  return p;
}

Complex ProjectionToComplex(const Projection &p) {
  return Complex(p.cos_amp, -p.sin_amp);
}

Real ByCell(ParameterInput *pin, Mesh *pm, const int q) {
  const Real by_amp = pin->GetOrAddReal("problem", "by_amp", 0.0);
  const Real k_wave = Wavenumber(pin, pm);
  const Real xmin = pm->mesh_size.x1min;
  return by_amp*std::sin(k_wave*(XCenter(pm, q) - xmin));
}

Real BMagCell(ParameterInput *pin, Mesh *pm, const int q) {
  const Real bx0 = pin->GetOrAddReal("problem", "b0", 1.0);
  const Real bz0 = pin->GetOrAddReal("problem", "bz0", 0.0);
  const Real by = ByCell(pin, pm, q);
  return std::sqrt(SQR(bx0) + SQR(by) + SQR(bz0));
}

Real LimitedHeatFlux(const Real q_unlimited, const Real q_max) {
  return cgl::LimitedHeatFlux(q_unlimited, q_max);
}

KOKKOS_INLINE_FUNCTION
Real EigenRealSpacePerturbation(const Real amp, const Real re, const Real im,
                                const Real c, const Real s) {
  return amp*(re*c - im*s);
}

bool GetBooleanOrFalse(ParameterInput *pin, const std::string &block,
                       const std::string &name) {
  return pin->DoesParameterExist(block, name) ? pin->GetBoolean(block, name) : false;
}

Real GetRealOrZero(ParameterInput *pin, const std::string &block,
                   const std::string &name) {
  return pin->DoesParameterExist(block, name) ? pin->GetReal(block, name) : 0.0;
}

Real BackgroundCollisionFrequency(ParameterInput *pin) {
  return GetRealOrZero(pin, "mhd", "nu_coll");
}

Real EffectiveCollisionFrequency(ParameterInput *pin, Mesh *pm,
                                 const Real ppar, const Real pperp,
                                 const Real bx, const Real by, const Real bz) {
  const EOS_Data &eos = pm->pmb_pack->pmhd->peos->eos_data;
  const bool backup = cgl::EffectiveBackupLimiter(
      eos.backup_lim, true, eos.mlim || eos.flim,
      GetBooleanOrFalse(pin, "mhd", "cgl_lf_strict_admissibility"));
  const Real bsqr = SQR(bx) + SQR(by) + SQR(bz);
  return std::max(eos.nu_coll, static_cast<Real>(0.0)) +
         cgl::LimiterCollisionRate(ppar, pperp, bsqr, eos, backup);
}

Real FaceCParallel(ParameterInput *pin, const Real rho, const Real ppar) {
  const std::string coeff_mode =
      pin->GetOrAddString("mhd", "lf_coefficient_mode", "local");
  if (coeff_mode == "background") {
    return pin->GetOrAddReal("mhd", "lf_c_parallel0",
                             std::sqrt(std::max(ppar/rho, static_cast<Real>(0.0))));
  }
  if (coeff_mode == "local") {
    const Real tfloor = pin->GetOrAddReal("mhd", "tfloor", 1.0e-30);
    return std::sqrt(std::max(ppar/rho, tfloor));
  }
  Fail("<mhd>/lf_coefficient_mode must be local or background");
}

Real ChiPerp(const Real cpar, const Real lf_k, const Real nu_eff) {
  const Real denom = cgl::kSqrtTwoPi*cpar*lf_k + 2.0*nu_eff;
  return (denom > 0.0) ? static_cast<Real>(2.0)*SQR(cpar)/denom : 0.0;
}

Real ChiParallel(const Real cpar, const Real lf_k, const Real nu_eff) {
  const Real denom = cgl::kSqrtEightPi*cpar*lf_k + cgl::kThreePiMinusEight*nu_eff;
  return (denom > 0.0) ? static_cast<Real>(8.0)*SQR(cpar)/denom : 0.0;
}

Real GradBMomentFlux(ParameterInput *pin, Mesh *pm, const int face) {
  const int nx1 = pm->mb_indcs.nx1;
  const Real dx = (pm->mesh_size.x1max - pm->mesh_size.x1min)/static_cast<Real>(nx1);
  const int left = (face - 1 + nx1)%nx1;
  const int right = face%nx1;
  const Real bx = pin->GetOrAddReal("problem", "b0", 1.0);
  const Real bmag_face = 0.5*BMagCell(pin, pm, left) +
                           0.5*BMagCell(pin, pm, right);
  const Real bhx = bx/bmag_face;
  const Real grad_b_x = (BMagCell(pin, pm, right) - BMagCell(pin, pm, left))/dx;
  const Real gradpar_b = bhx*grad_b_x;
  const Real rho = pin->GetOrAddReal("problem", "rho0", 1.0);
  const Real ppar = pin->GetOrAddReal("problem", "ppar0", 1.0);
  const Real pperp = pin->GetOrAddReal("problem", "pperp0", 1.2);
  const Real lf_k = pin->GetReal("mhd", "lf_k_parallel");
  const Real cpar = FaceCParallel(pin, rho, ppar);
  // The limiter depends only on field magnitude; supply the face Bbar.
  const Real nu_eff = EffectiveCollisionFrequency(
      pin, pm, ppar, pperp, bmag_face, 0.0, 0.0);
  const Real chi_perp = ChiPerp(cpar, lf_k, nu_eff);
  const Real qperp_l = -chi_perp*(-pperp*(1.0 - pperp/ppar)*gradpar_b/bmag_face);
  const Real qperp = LimitedHeatFlux(qperp_l, cgl::kSqrtTwoOverPi*cpar*pperp);
  return bhx*qperp/bmag_face;
}

struct DecayState {
  Real tpar;
  Real tperp;
};

DecayState IntegrateDecayReference(ParameterInput *pin, Mesh *pm, const TestMode mode,
                                   const Real chi_parallel, const Real chi_perp) {
  const Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  const Real ppar0 = pin->GetOrAddReal("problem", "ppar0", 1.0);
  const Real pperp0 = pin->GetOrAddReal("problem", "pperp0", 1.0);
  const Real amp = pin->GetOrAddReal("problem", "amp", 1.0e-4);
  const Real k_wave = Wavenumber(pin, pm);
  const Real nu_coll =
      std::max(BackgroundCollisionFrequency(pin), static_cast<Real>(0.0));
  DecayState s{0.0, 0.0};
  if (mode == TestMode::parallel_decay) {
    s.tpar = (ppar0/rho0)*amp;
  } else {
    s.tperp = (pperp0/rho0)*amp;
  }

  // Independent continuous linear system for coupled temperature amplitudes.
  // Background collisions conserve (T_parallel + 2*T_perp)/3 and decay their
  // difference at nu_coll. Use its exact 2x2 matrix exponential as the reference.
  const Real a = -chi_parallel*SQR(k_wave) - TWO_3RDS*nu_coll;
  const Real b = TWO_3RDS*nu_coll;
  const Real c = ONE_3RD*nu_coll;
  const Real d = -chi_perp*SQR(k_wave) - ONE_3RD*nu_coll;
  const Real half_trace = 0.5*(a + d);
  const Real gap = std::sqrt(SQR(0.5*(a - d)) + b*c);
  const Real ep = std::exp((half_trace + gap)*pm->time);
  const Real em = std::exp((half_trace - gap)*pm->time);
  const Real even = 0.5*(ep + em);
  const Real odd = (gap > 0.0) ? 0.5*(ep - em)/gap : pm->time*ep;
  return {even*s.tpar + odd*((a - half_trace)*s.tpar + b*s.tperp),
          even*s.tperp + odd*(c*s.tpar + (d - half_trace)*s.tperp)};
}

void CheckCollisionRelaxation(ParameterInput *pin, Mesh *pm) {
  auto *pmhd = pm->pmb_pack->pmhd;
  auto w = HostCopy(pmhd->w0);
  const int is = pm->mb_indcs.is;
  const int js = pm->mb_indcs.js;
  const int ks = pm->mb_indcs.ks;
  const int nx1 = pm->mb_indcs.nx1;
  const Real ppar0 = pin->GetOrAddReal("problem", "ppar0", 1.0);
  const Real pperp0 = pin->GetOrAddReal("problem", "pperp0", 1.0);
  const Real nu_coll = BackgroundCollisionFrequency(pin);
  const Real rel_tol = pin->GetOrAddReal("problem", "collision_rel_tol", 1.0e-12);
  Require(nu_coll >= 0.0, "collision_relaxation requires nonnegative <mhd>/nu_coll");

  Real measured_paniso = 0.0;
  Real measured_piso = 0.0;
  for (int q = 0; q < nx1; ++q) {
    const Real ppar = w(0,IPR,ks,js,is + q);
    const Real pperp = w(0,IPP,ks,js,is + q);
    measured_paniso += pperp - ppar;
    measured_piso += ONE_3RD*ppar + TWO_3RDS*pperp;
  }
  measured_paniso /= static_cast<Real>(nx1);
  measured_piso /= static_cast<Real>(nx1);

  const Real initial_paniso = pperp0 - ppar0;
  const Real initial_piso = ONE_3RD*ppar0 + TWO_3RDS*pperp0;
  const Real expected_paniso = initial_paniso*std::exp(-nu_coll*pm->time);
  RequireRelative("collision_relaxation anisotropy", measured_paniso,
                  expected_paniso, rel_tol);
  RequireRelative("collision_relaxation isotropic pressure", measured_piso,
                  initial_piso, rel_tol);
  std::cout << "CGL LF collision_relaxation passed: measured_paniso="
            << measured_paniso << " expected_paniso=" << expected_paniso
            << " time=" << pm->time << " nu_coll=" << nu_coll << std::endl;
}

void CheckDecay(ParameterInput *pin, Mesh *pm, const TestMode mode) {
  auto *pmhd = pm->pmb_pack->pmhd;
  auto w = HostCopy(pmhd->w0);
  const int pidx = (mode == TestMode::parallel_decay) ? IPR : IPP;
  const Projection projection = ProjectTemperature(w, pin, pm, pidx);

  const Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  const Real p0 = (mode == TestMode::parallel_decay) ?
                  pin->GetOrAddReal("problem", "ppar0", 1.0) :
                  pin->GetOrAddReal("problem", "pperp0", 1.0);
  const Real amp = pin->GetOrAddReal("problem", "amp", 1.0e-4);
  const Real rel_tol = pin->GetOrAddReal("problem", "decay_rel_tol", 5.0e-2);
  const Real phase_tol = pin->GetOrAddReal("problem", "phase_rel_tol", 2.0e-2);
  const Real lf_k = pin->GetReal("mhd", "lf_k_parallel");
  const Real k_wave = Wavenumber(pin, pm);
  const Real bx0 = pin->GetOrAddReal("problem", "b0", 1.0);
  const Real by0 = pin->GetOrAddReal("problem", "by0", 0.0);
  const Real bz0 = pin->GetOrAddReal("problem", "bz0", 0.0);
  const Real ppar0 = pin->GetOrAddReal("problem", "ppar0", 1.0);
  const Real pperp0 = pin->GetOrAddReal("problem", "pperp0", 1.0);
  const Real cpar0 = FaceCParallel(pin, rho0, ppar0);
  const Real nu_eff = EffectiveCollisionFrequency(pin, pm, ppar0, pperp0, bx0, by0, bz0);
  const Real chi_parallel = ChiParallel(cpar0, lf_k, nu_eff);
  const Real chi_perp = ChiPerp(cpar0, lf_k, nu_eff);
  const DecayState expected_state =
      IntegrateDecayReference(pin, pm, mode, chi_parallel, chi_perp);
  const Real expected_mode_amp = (mode == TestMode::parallel_decay) ?
                                 expected_state.tpar : expected_state.tperp;
  const Real initial_amp = (p0/rho0)*amp;
  const Real expected_amp = std::abs(expected_mode_amp);
  const Real measured_amp = std::abs(projection.sin_amp);
  const Real rel_err = RelativeError(measured_amp, expected_amp);
  const Real cos_limit = phase_tol*std::abs(initial_amp);

  RequireRelative(std::string(ModeName(mode)) + " Fourier amplitude",
                  measured_amp, expected_amp, rel_tol);
  Require(std::abs(projection.cos_amp) <= cos_limit,
          std::string(ModeName(mode)) + " developed an unexpected cosine component");

  auto out = OpenValidationDataFile(pin, "decay");
  if (out) {
    out << "metric,value\n"
        << "time," << pm->time << "\n"
        << "initial_amp," << initial_amp << "\n"
        << "measured_amp," << measured_amp << "\n"
        << "measured_sin_amp," << projection.sin_amp << "\n"
        << "measured_cos_amp," << projection.cos_amp << "\n"
        << "expected_amp," << expected_amp << "\n"
        << "expected_tpar_amp," << expected_state.tpar << "\n"
        << "expected_tperp_amp," << expected_state.tperp << "\n"
        << "chi_parallel," << chi_parallel << "\n"
        << "chi_perp," << chi_perp << "\n"
        << "nu_eff," << nu_eff << "\n"
        << "nu_coll," << BackgroundCollisionFrequency(pin) << "\n"
        << "k_wave," << k_wave << "\n"
        << "rel_err," << rel_err << "\n"
        << "rel_tol," << rel_tol << "\n"
        << "cos_limit," << cos_limit << "\n";
  }

  std::cout << "CGL LF " << ModeName(mode) << " passed: measured_amp="
            << measured_amp << " expected_amp=" << expected_amp
            << " rel_err=" << rel_err << " rel_tol=" << rel_tol
            << " cos_amp=" << projection.cos_amp << " cos_limit=" << cos_limit
            << std::endl;
}

void CheckRotatedDecay(ParameterInput *pin, Mesh *pm) {
  auto *pmhd = pm->pmb_pack->pmhd;
  auto w = HostCopy(pmhd->w0);
  const Projection projection = ProjectRotatedTemperature(w, pin, pm, IPR);
  const Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  const Real ppar0 = pin->GetOrAddReal("problem", "ppar0", 1.0);
  const Real pperp0 = pin->GetOrAddReal("problem", "pperp0", 1.0);
  const Real amp = pin->GetOrAddReal("problem", "amp", 1.0e-4);
  const Real rel_tol = pin->GetOrAddReal("problem", "decay_rel_tol", 5.0e-2);
  const Real phase_tol = pin->GetOrAddReal("problem", "phase_rel_tol", 2.0e-2);
  const Real bx0 = pin->GetOrAddReal("problem", "b0", 1.0);
  const Real by0 = pin->GetOrAddReal("problem", "by0", 0.0);
  const Real bz0 = pin->GetOrAddReal("problem", "bz0", 0.0);
  const Real bmag = std::sqrt(SQR(bx0) + SQR(by0) + SQR(bz0));
  Require(bmag > 0.0, "rotated_decay requires a nonzero mean field");
  const RotatedWave wave = RotatedWavenumber(pin, pm);
  const Real k_parallel = (bx0*wave.kx + by0*wave.ky + bz0*wave.kz)/bmag;
  const Real cpar0 = FaceCParallel(pin, rho0, ppar0);
  const Real nu_eff = EffectiveCollisionFrequency(pin, pm, ppar0, pperp0, bx0, by0, bz0);
  const Real chi_parallel =
      ChiParallel(cpar0, pin->GetReal("mhd", "lf_k_parallel"), nu_eff);
  const Real initial_amp = (ppar0/rho0)*amp;
  const Real expected_amp =
      initial_amp*std::exp(-chi_parallel*SQR(k_parallel)*pm->time);
  const Real measured_amp = std::abs(projection.sin_amp);
  const Real rel_err = RelativeError(measured_amp, expected_amp);
  RequireRelative("rotated_decay Fourier amplitude", measured_amp, expected_amp, rel_tol);
  Require(std::abs(projection.cos_amp) <= phase_tol*std::abs(initial_amp),
          "rotated_decay developed an unexpected cosine component");
  auto out = OpenValidationDataFile(pin, "rotated_decay");
  if (out) {
    out << "metric,value\n"
        << "time," << pm->time << "\n"
        << "initial_amp," << initial_amp << "\n"
        << "measured_amp," << measured_amp << "\n"
        << "expected_amp," << expected_amp << "\n"
        << "chi_parallel," << chi_parallel << "\n"
        << "nu_eff," << nu_eff << "\n"
        << "k_parallel," << k_parallel << "\n"
        << "rel_err," << rel_err << "\n"
        << "rel_tol," << rel_tol << "\n";
  }
  std::cout << "CGL LF rotated_decay passed: k_parallel=" << k_parallel
            << " measured_amp=" << measured_amp << " expected_amp=" << expected_amp
            << " rel_err=" << rel_err << " rel_tol=" << rel_tol << std::endl;
}

Real InitialParallelPressure(ParameterInput *pin, Mesh *pm, const int q);
Real InitialPerpPressure(ParameterInput *pin, Mesh *pm, const int q);

void CheckLowField(ParameterInput *pin, Mesh *pm) {
  auto *pmhd = pm->pmb_pack->pmhd;
  auto w = HostCopy(pmhd->w0);
  const int is = pm->mb_indcs.is;
  const int js = pm->mb_indcs.js;
  const int ks = pm->mb_indcs.ks;
  const int nx1 = pm->mb_indcs.nx1;
  const Real bx0 = pin->GetOrAddReal("problem", "b0", 1.0);
  const Real by0 = pin->GetOrAddReal("problem", "by0", 0.0);
  const Real bz0 = pin->GetOrAddReal("problem", "bz0", 0.0);
  const Real bmag = std::sqrt(SQR(bx0) + SQR(by0) + SQR(bz0));
  const Real bfloor = pin->GetReal("mhd", "bfloor");
  Require(bmag <= bfloor, "low_field requires |B| <= <mhd>/bfloor");
  Real maximum_error = 0.0;
  for (int q = 0; q < nx1; ++q) {
    const Real ppar = w(0,IPR,ks,js,is + q);
    const Real pperp = w(0,IPP,ks,js,is + q);
    Require(std::isfinite(ppar) && std::isfinite(pperp) &&
            ppar > 0.0 && pperp > 0.0,
            "low_field produced a nonfinite or nonpositive pressure");
    maximum_error = std::max(maximum_error,
        std::abs(ppar - InitialParallelPressure(pin, pm, q)));
    maximum_error = std::max(maximum_error,
        std::abs(pperp - InitialPerpPressure(pin, pm, q)));
  }
  const Real abs_tol = pin->GetOrAddReal("problem", "low_field_abs_tol", 1.0e-13);
  Require(maximum_error <= abs_tol,
          "low_field transported pressure at or below the magnetic-field floor");
  std::cout << "CGL LF low_field passed: bmag=" << bmag
            << " bfloor=" << bfloor << " maximum_pressure_change="
            << maximum_error << std::endl;
}

void CheckGradB(ParameterInput *pin, Mesh *pm) {
  auto *pmhd = pm->pmb_pack->pmhd;
  auto w = HostCopy(pmhd->w0);
  auto bcc = HostCopy(pmhd->bcc0);
  const int is = pm->mb_indcs.is;
  const int js = pm->mb_indcs.js;
  const int ks = pm->mb_indcs.ks;
  const int nx1 = pm->mb_indcs.nx1;
  const Real dx = (pm->mesh_size.x1max - pm->mesh_size.x1min)/static_cast<Real>(nx1);
  const Real pperp0 = pin->GetOrAddReal("problem", "pperp0", 1.2);
  const Real rel_tol = pin->GetOrAddReal("problem", "grad_b_rel_tol", 1.0e-1);
  auto out = OpenValidationDataFile(pin, "grad_b");
  if (out) {
    out << "x,measured_delta,expected_delta,bmag_initial,bmag_final\n";
  }

  Real err2 = 0.0;
  Real ref2 = 0.0;
  Real dot = 0.0;
  for (int q = 0; q < nx1; ++q) {
    const int i = is + q;
    const Real flux_r = GradBMomentFlux(pin, pm, q + 1);
    const Real flux_l = GradBMomentFlux(pin, pm, q);
    const Real expected_delta = pm->time*(-(flux_r - flux_l)/dx);
    const Real bmag_initial = BMagCell(pin, pm, q);
    const Real bmag_final = std::sqrt(SQR(bcc(0,IBX,ks,js,i)) +
                                      SQR(bcc(0,IBY,ks,js,i)) +
                                      SQR(bcc(0,IBZ,ks,js,i)));
    const Real measured_delta = w(0,IPP,ks,js,i)/bmag_final - pperp0/bmag_initial;
    err2 += SQR(measured_delta - expected_delta);
    ref2 += SQR(expected_delta);
    dot += measured_delta*expected_delta;
    if (out) {
      out << XCenter(pm, q) << "," << measured_delta << "," << expected_delta << ","
          << bmag_initial << "," << bmag_final << "\n";
    }
  }

  Require(ref2 > 0.0, "grad_b reference response is zero");
  const Real rel_err = std::sqrt(err2/ref2);
  Require(rel_err <= rel_tol, "grad_b response magnitude is outside tolerance");
  Require(dot > 0.0, "grad_b response has the wrong sign");
  std::cout << "CGL LF grad_b passed: rms_rel_err=" << rel_err
            << " rel_tol=" << rel_tol << std::endl;
}

Real InitialParallelPressure(ParameterInput *pin, Mesh *pm, const int q) {
  const Real ppar0 = pin->GetOrAddReal("problem", "ppar0", 1.0);
  const Real amp = pin->GetOrAddReal("problem", "amp", 0.5);
  const Real k_wave = Wavenumber(pin, pm);
  const Real xmin = pm->mesh_size.x1min;
  return ppar0*(1.0 + amp*std::sin(k_wave*(XCenter(pm, q) - xmin)));
}

Real InitialPerpPressure(ParameterInput *pin, Mesh *pm, const int q) {
  const Real pperp0 = pin->GetOrAddReal("problem", "pperp0", 1.0);
  const Real amp = pin->GetOrAddReal("problem", "amp", 0.5);
  const Real k_wave = Wavenumber(pin, pm);
  const Real xmin = pm->mesh_size.x1min;
  return pperp0*(1.0 + amp*std::sin(k_wave*(XCenter(pm, q) - xmin)));
}

Real LimitedParallelHeatFlux(ParameterInput *pin, Mesh *pm, const int face,
                             Real &q_unlimited, Real &q_max,
                             const bool force_collisionless = false) {
  const int nx1 = pm->mb_indcs.nx1;
  const Real dx = (pm->mesh_size.x1max - pm->mesh_size.x1min)/static_cast<Real>(nx1);
  const int left = (face - 1 + nx1)%nx1;
  const int right = face%nx1;
  const Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  const Real lf_k = pin->GetReal("mhd", "lf_k_parallel");
  const Real bx = pin->GetOrAddReal("problem", "b0", 1.0);
  const Real by = pin->GetOrAddReal("problem", "by0", 0.0);
  const Real bz = pin->GetOrAddReal("problem", "bz0", 0.0);
  const Real ppar_l = InitialParallelPressure(pin, pm, left);
  const Real ppar_r = InitialParallelPressure(pin, pm, right);
  const Real pperp_l = InitialPerpPressure(pin, pm, left);
  const Real pperp_r = InitialPerpPressure(pin, pm, right);
  const Real ppar_face = std::max(static_cast<Real>(0.5)*(ppar_l + ppar_r),
                                  static_cast<Real>(1.0e-30));
  const Real pperp_face = std::max(static_cast<Real>(0.5)*(pperp_l + pperp_r),
                                   static_cast<Real>(1.0e-30));
  const Real cpar = FaceCParallel(pin, rho0, ppar_face);
  const Real nu_eff = force_collisionless ? 0.0 :
      EffectiveCollisionFrequency(pin, pm, ppar_face, pperp_face, bx, by, bz);
  const Real chi_parallel = ChiParallel(cpar, lf_k, nu_eff);
  const Real grad_tpar = (ppar_r/rho0 - ppar_l/rho0)/dx;
  q_unlimited = -chi_parallel*rho0*grad_tpar;
  q_max = cgl::kSqrtEightOverPi*cpar*ppar_face;
  return LimitedHeatFlux(q_unlimited, q_max);
}

Real LimitedPerpHeatFlux(ParameterInput *pin, Mesh *pm, const int face,
                         Real &q_unlimited, Real &q_max,
                         const bool force_collisionless = false) {
  const int nx1 = pm->mb_indcs.nx1;
  const Real dx = (pm->mesh_size.x1max - pm->mesh_size.x1min)/static_cast<Real>(nx1);
  const int left = (face - 1 + nx1)%nx1;
  const int right = face%nx1;
  const Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  const Real lf_k = pin->GetReal("mhd", "lf_k_parallel");
  const Real bx = pin->GetOrAddReal("problem", "b0", 1.0);
  const Real by = pin->GetOrAddReal("problem", "by0", 0.0);
  const Real bz = pin->GetOrAddReal("problem", "bz0", 0.0);
  const Real ppar_l = InitialParallelPressure(pin, pm, left);
  const Real ppar_r = InitialParallelPressure(pin, pm, right);
  const Real pperp_l = InitialPerpPressure(pin, pm, left);
  const Real pperp_r = InitialPerpPressure(pin, pm, right);
  const Real ppar_face = std::max(static_cast<Real>(0.5)*(ppar_l + ppar_r),
                                  static_cast<Real>(1.0e-30));
  const Real pperp_face = std::max(static_cast<Real>(0.5)*(pperp_l + pperp_r),
                                   static_cast<Real>(1.0e-30));
  const Real cpar = FaceCParallel(pin, rho0, ppar_face);
  const Real nu_eff = force_collisionless ? 0.0 :
      EffectiveCollisionFrequency(pin, pm, ppar_face, pperp_face, bx, by, bz);
  const Real chi_perp = ChiPerp(cpar, lf_k, nu_eff);
  const Real grad_tperp = (pperp_r/rho0 - pperp_l/rho0)/dx;
  q_unlimited = -chi_perp*rho0*grad_tperp;
  q_max = cgl::kSqrtTwoOverPi*cpar*pperp_face;
  return LimitedHeatFlux(q_unlimited, q_max);
}

void CheckFluxLimiter(ParameterInput *pin, Mesh *pm) {
  auto *pmhd = pm->pmb_pack->pmhd;
  auto w = HostCopy(pmhd->w0);
  const int is = pm->mb_indcs.is;
  const int js = pm->mb_indcs.js;
  const int ks = pm->mb_indcs.ks;
  const int nx1 = pm->mb_indcs.nx1;
  const Real dx = (pm->mesh_size.x1max - pm->mesh_size.x1min)/static_cast<Real>(nx1);
  const Real rel_tol = pin->GetOrAddReal("problem", "flux_limiter_rel_tol", 2.0e-1);
  const Real min_unlimited_ratio =
      pin->GetOrAddReal("problem", "flux_limiter_min_unlimited_ratio", 10.0);
  auto face_out = OpenValidationDataFile(pin, "flux_limiter_faces");
  if (face_out) {
    face_out << "x_face,q_unlimited,q_limited,q_max,"
             << "q_perp_unlimited,q_perp_limited,q_perp_max\n";
  }
  auto cell_out = OpenValidationDataFile(pin, "flux_limiter_cells");
  if (cell_out) {
    cell_out << "x,measured_ppar,expected_ppar,initial_ppar,"
             << "measured_pperp,expected_pperp,initial_pperp\n";
  }

  Real err2 = 0.0;
  Real ref2 = 0.0;
  Real max_ratio_parallel = 0.0;
  Real max_ratio_perp = 0.0;
  for (int face = 0; face <= nx1; ++face) {
    Real qpar_l = 0.0, qpar_max = 0.0;
    const Real qpar = LimitedParallelHeatFlux(pin, pm, face, qpar_l, qpar_max);
    Real qperp_l = 0.0, qperp_max = 0.0;
    const Real qperp = LimitedPerpHeatFlux(pin, pm, face, qperp_l, qperp_max);
    max_ratio_parallel = std::max(max_ratio_parallel,
        std::abs(qpar_l)/std::max(qpar_max, static_cast<Real>(1.0e-30)));
    max_ratio_perp = std::max(max_ratio_perp,
        std::abs(qperp_l)/std::max(qperp_max, static_cast<Real>(1.0e-30)));
    if (face_out) {
      const Real x_face = pm->mesh_size.x1min + static_cast<Real>(face)*dx;
      face_out << x_face << "," << qpar_l << "," << qpar << "," << qpar_max
               << "," << qperp_l << "," << qperp << "," << qperp_max << "\n";
    }
    Require(std::abs(qpar) <= qpar_max*(1.0 + 1.0e-12),
            "limited parallel heat flux exceeds q_max");
    Require(std::abs(qperp) <= qperp_max*(1.0 + 1.0e-12),
            "limited perpendicular heat flux exceeds q_max");
    if (qpar_l != 0.0) {
      Require(qpar*qpar_l > 0.0, "limited parallel heat flux changed sign");
    }
    if (qperp_l != 0.0) {
      Require(qperp*qperp_l > 0.0, "limited perpendicular heat flux changed sign");
    }
  }

  for (int q = 0; q < nx1; ++q) {
    const int i = is + q;
    Real q_r_l = 0.0, q_r_max = 0.0;
    Real q_l_l = 0.0, q_l_max = 0.0;
    const Real flux_r = LimitedParallelHeatFlux(pin, pm, q + 1, q_r_l, q_r_max);
    const Real flux_l = LimitedParallelHeatFlux(pin, pm, q, q_l_l, q_l_max);
    Real qperp_r_l = 0.0, qperp_r_max = 0.0;
    Real qperp_l_l = 0.0, qperp_l_max = 0.0;
    const Real perp_flux_r = LimitedPerpHeatFlux(pin, pm, q + 1,
                                                 qperp_r_l, qperp_r_max);
    const Real perp_flux_l = LimitedPerpHeatFlux(pin, pm, q,
                                                 qperp_l_l, qperp_l_max);
    const Real expected_ppar = InitialParallelPressure(pin, pm, q)
                             + pm->time*(-(flux_r - flux_l)/dx);
    const Real expected_pperp = InitialPerpPressure(pin, pm, q)
                              + pm->time*(-(perp_flux_r - perp_flux_l)/dx);
    const Real measured_ppar = w(0,IPR,ks,js,i);
    const Real measured_pperp = w(0,IPP,ks,js,i);
    err2 += SQR(measured_ppar - expected_ppar)
          + SQR(measured_pperp - expected_pperp);
    ref2 += SQR(expected_ppar - InitialParallelPressure(pin, pm, q))
          + SQR(expected_pperp - InitialPerpPressure(pin, pm, q));
    if (cell_out) {
      cell_out << XCenter(pm, q) << "," << measured_ppar << "," << expected_ppar
               << "," << InitialParallelPressure(pin, pm, q)
               << "," << measured_pperp << "," << expected_pperp
               << "," << InitialPerpPressure(pin, pm, q) << "\n";
    }
  }

  Require(max_ratio_parallel >= min_unlimited_ratio,
          "test did not enter the strongly limited parallel heat-flux regime");
  Require(max_ratio_perp >= min_unlimited_ratio,
          "test did not enter the strongly limited perpendicular heat-flux regime");
  Require(ref2 > 0.0, "flux limiter reference response is zero");
  const Real rel_err = std::sqrt(err2/ref2);
  if (rel_err > rel_tol) {
    std::cout << "flux_limiter rel_err=" << rel_err << " rel_tol=" << rel_tol
              << std::endl;
    Fail("flux limiter update is outside tolerance");
  }
  std::cout << "CGL LF flux_limiter passed: rms_rel_err=" << rel_err
            << " rel_tol=" << rel_tol
            << " max_parallel_unlimited_over_qmax=" << max_ratio_parallel
            << " max_perp_unlimited_over_qmax=" << max_ratio_perp
            << " min_unlimited_over_qmax=" << min_unlimited_ratio << std::endl;
}

// Compare production face fluxes with a prescribed total collision rate, without
// calling the rate helper in the reference. This fixture has constant B=(1,0,0).
void CheckPrescribedCollisionRateFlux(ParameterInput *pin, Mesh *pm) {
  auto *pmhd = pm->pmb_pack->pmhd;
  auto *lf = pmhd->pcgl_lf;
  Require(pm->ncycle == 0 && lf != nullptr && !lf->lf_coeff_local,
          "prescribed-rate flux test requires nlim=0 and background LF coefficients");
  const Real nu = pin->GetReal("problem", "expected_nu_eff");
  const Real cpar = lf->lf_c_parallel0;
  const Real pi = std::acos(static_cast<Real>(-1.0));
  const Real chi_par = 8.0*SQR(cpar)/(std::sqrt(8.0*pi)*cpar*lf->lf_k_parallel +
                                      (3.0*pi - 8.0)*nu);
  const Real chi_perp = 2.0*SQR(cpar)/(std::sqrt(2.0*pi)*cpar*lf->lf_k_parallel +
                                     2.0*nu);
  const auto &indcs = pm->mb_indcs;
  DvceFaceFld5D<Real> flux("prescribed_rate_flux", 1, pmhd->nmhd, 1, 1,
                           indcs.nx1 + 2*indcs.ng);
  lf->AddHeatFluxes(pmhd->w0, pmhd->bcc0, pmhd->b0, pmhd->peos->eos_data,
                      1.0, 1.0, flux);
  auto actual = HostCopy(flux.x1f);
  auto w = HostCopy(pmhd->w0);
  auto b = HostCopy(pmhd->bcc0);
  const int js = indcs.js, ks = indcs.ks;
  const Real dx = (pm->mesh_size.x1max - pm->mesh_size.x1min)/indcs.nx1;
  const Real tolerance = 256.0*std::numeric_limits<Real>::epsilon();
  for (int i = indcs.is; i <= indcs.ie + 1; ++i) {
    Require(b(0,IBX,ks,js,i) == 1.0 && b(0,IBY,ks,js,i) == 0.0 &&
            b(0,IBZ,ks,js,i) == 0.0, "prescribed-rate test requires B=(1,0,0)");
    const Real dl = w(0,IDN,ks,js,i-1), dr = w(0,IDN,ks,js,i);
    const Real pl = w(0,IPR,ks,js,i-1), pr = w(0,IPR,ks,js,i);
    const Real tl = w(0,IPP,ks,js,i-1), tr = w(0,IPP,ks,js,i);
    const Real rho = 0.5*(dl + dr);
    const Real qpar_l = -chi_par*rho*(pr/dr - pl/dl)/dx;
    const Real qperp_l = -chi_perp*rho*(tr/dr - tl/dl)/dx;
    const Real qpar_max = std::sqrt(8.0/pi)*cpar*0.5*(pl + pr);
    const Real qperp_max = std::sqrt(2.0/pi)*cpar*0.5*(tl + tr);
    const Real qpar = qpar_l/(1.0 + std::abs(qpar_l)/qpar_max);
    const Real qperp = qperp_l/(1.0 + std::abs(qperp_l)/qperp_max);
    const Real expected_energy = qperp + 0.5*qpar;
    Require(std::abs(actual(0,IEN,ks,js,i) - expected_energy) <=
                tolerance*std::max(static_cast<Real>(1.0), std::abs(expected_energy)),
            "prescribed collision rate disagrees with production energy flux");
    Require(std::abs(actual(0,IAN,ks,js,i) - qperp) <=
                tolerance*std::max(static_cast<Real>(1.0), std::abs(qperp)),
            "prescribed collision rate disagrees with production moment flux");
  }
  std::cout << "CGL LF prescribed collision-rate flux passed: nu_eff=" << nu << std::endl;
}

void CheckLimiterHeatFluxSuppression(ParameterInput *pin, Mesh *pm) {
  if (pin->DoesParameterExist("problem", "expected_nu_eff")) {
    CheckPrescribedCollisionRateFlux(pin, pm);
  }
  auto *pmhd = pm->pmb_pack->pmhd;
  auto w = HostCopy(pmhd->w0);
  const int is = pm->mb_indcs.is;
  const int js = pm->mb_indcs.js;
  const int ks = pm->mb_indcs.ks;
  const int nx1 = pm->mb_indcs.nx1;
  const Real dx = (pm->mesh_size.x1max - pm->mesh_size.x1min)/static_cast<Real>(nx1);
  const Real max_suppression_ratio =
      pin->GetOrAddReal("problem", "limiter_suppression_max_ratio", 5.0e-2);
  const Real min_collisionless_ratio =
      pin->GetOrAddReal("problem", "limiter_suppression_min_collisionless_ratio", 10.0);
  auto face_out = OpenValidationDataFile(pin, "limiter_heat_flux_suppression_faces");
  if (face_out) {
    face_out << "x_face,qpar_unlimited,qpar_collisionless,qpar_limited,qpar_max,"
             << "qperp_unlimited,qperp_collisionless,qperp_limited,qperp_max\n";
  }

  Real max_parallel_ratio = 0.0;
  Real max_perp_ratio = 0.0;
  Real max_parallel_collisionless_over_qmax = 0.0;
  Real max_perp_collisionless_over_qmax = 0.0;
  for (int face = 0; face <= nx1; ++face) {
    Real qpar_l = 0.0, qpar_max = 0.0;
    const Real qpar = LimitedParallelHeatFlux(pin, pm, face, qpar_l, qpar_max);
    Real qpar_collisionless_l = 0.0, qpar_collisionless_max = 0.0;
    (void)LimitedParallelHeatFlux(pin, pm, face, qpar_collisionless_l,
                                  qpar_collisionless_max, true);
    Real qperp_l = 0.0, qperp_max = 0.0;
    const Real qperp = LimitedPerpHeatFlux(pin, pm, face, qperp_l, qperp_max);
    Real qperp_collisionless_l = 0.0, qperp_collisionless_max = 0.0;
    (void)LimitedPerpHeatFlux(pin, pm, face, qperp_collisionless_l,
                              qperp_collisionless_max, true);

    max_parallel_collisionless_over_qmax = std::max(
        max_parallel_collisionless_over_qmax,
        std::abs(qpar_collisionless_l)/std::max(qpar_max, static_cast<Real>(1.0e-30)));
    max_perp_collisionless_over_qmax = std::max(
        max_perp_collisionless_over_qmax,
        std::abs(qperp_collisionless_l)/std::max(qperp_max, static_cast<Real>(1.0e-30)));
    if (std::abs(qpar_collisionless_l) > 1.0e-30) {
      max_parallel_ratio = std::max(max_parallel_ratio,
                                    std::abs(qpar_l/qpar_collisionless_l));
      Require(qpar_l*qpar_collisionless_l > 0.0,
              "limiter-suppressed parallel heat flux changed sign");
    }
    if (std::abs(qperp_collisionless_l) > 1.0e-30) {
      max_perp_ratio = std::max(max_perp_ratio,
                                std::abs(qperp_l/qperp_collisionless_l));
      Require(qperp_l*qperp_collisionless_l > 0.0,
              "limiter-suppressed perpendicular heat flux changed sign");
    }
    Require(std::abs(qpar) <= qpar_max*(1.0 + 1.0e-12),
            "limiter-suppressed parallel heat flux exceeds q_max");
    Require(std::abs(qperp) <= qperp_max*(1.0 + 1.0e-12),
            "limiter-suppressed perpendicular heat flux exceeds q_max");
    if (qpar_l != 0.0) {
      Require(qpar*qpar_l > 0.0, "limited parallel heat flux changed sign");
    }
    if (qperp_l != 0.0) {
      Require(qperp*qperp_l > 0.0, "limited perpendicular heat flux changed sign");
    }

    if (face_out) {
      const Real x_face = pm->mesh_size.x1min + static_cast<Real>(face)*dx;
      face_out << x_face << "," << qpar_l << "," << qpar_collisionless_l
               << "," << qpar << "," << qpar_max << ","
               << qperp_l << "," << qperp_collisionless_l
               << "," << qperp << "," << qperp_max << "\n";
    }
  }

  for (int q = 0; q < nx1; ++q) {
    const int i = is + q;
    Require(std::isfinite(w(0,IDN,ks,js,i)) && std::isfinite(w(0,IPR,ks,js,i)) &&
            std::isfinite(w(0,IPP,ks,js,i)),
            "limiter heat-flux suppression produced a nonfinite primitive state");
  }

  Require(max_parallel_collisionless_over_qmax >= min_collisionless_ratio,
          "parallel setup did not have collisionless |q_L| >> q_max");
  Require(max_perp_collisionless_over_qmax >= min_collisionless_ratio,
          "perpendicular setup did not have collisionless |q_L| >> q_max");
  Require(max_parallel_ratio <= max_suppression_ratio,
          "parallel limiter collisionality did not strongly suppress q_L");
  Require(max_perp_ratio <= max_suppression_ratio,
          "perpendicular limiter collisionality did not strongly suppress q_L");

  std::cout << "CGL LF limiter_heat_flux_suppression passed: "
            << "max_parallel_suppressed_over_collisionless=" << max_parallel_ratio
            << " max_perp_suppressed_over_collisionless=" << max_perp_ratio
            << " limit=" << max_suppression_ratio
            << " max_parallel_collisionless_over_qmax="
            << max_parallel_collisionless_over_qmax
            << " max_perp_collisionless_over_qmax="
            << max_perp_collisionless_over_qmax << std::endl;
}

void CheckLimiterStress(ParameterInput *pin, Mesh *pm) {
  auto *pmhd = pm->pmb_pack->pmhd;
  auto w = HostCopy(pmhd->w0);
  auto u = HostCopy(pmhd->u0);
  auto bcc = HostCopy(pmhd->bcc0);
  const int is = pm->mb_indcs.is;
  const int js = pm->mb_indcs.js;
  const int ks = pm->mb_indcs.ks;
  const int nx1 = pm->mb_indcs.nx1;
  const EOS_Data &eos = pmhd->peos->eos_data;
  const std::string limiter_kind =
      pin->GetOrAddString("problem", "limiter_kind", "mirror");

  const Real amp = pin->GetOrAddReal("problem", "amp", 1.0e-4);
  const Real ppar0 = pin->GetOrAddReal("problem", "ppar0", 1.0);
  const Real pperp0 = pin->GetOrAddReal("problem", "pperp0", 1.0);
  const bool has_lf = pin->DoesParameterExist("mhd", "cgl_heat_flux");
  // This cellwise reference describes one cycle without spatial heat transport.
  // Multi-cycle relaxation is checked independently from the per-cycle history.
  const bool check_relaxation = (amp == 0.0 || !has_lf) && pm->ncycle == 1 &&
      pin->GetString("mhd", "rsolver") == "advect";
  Real max_relaxation_error = 0.0;
  for (int q = 0; q < nx1; ++q) {
    const int i = is + q;
    const Real rho = w(0,IDN,ks,js,i);
    const Real ppar = w(0,IPR,ks,js,i);
    const Real pperp = w(0,IPP,ks,js,i);
    const Real bsqr = SQR(bcc(0,IBX,ks,js,i)) + SQR(bcc(0,IBY,ks,js,i)) +
                      SQR(bcc(0,IBZ,ks,js,i));
    const Real paniso = pperp - ppar;
    Require(std::isfinite(rho) && std::isfinite(ppar) && std::isfinite(pperp) &&
            std::isfinite(u(0,IEN,ks,js,i)) && std::isfinite(u(0,IAN,ks,js,i)),
            "limiter stress produced a nonfinite state");
    Require(rho > 0.0 && ppar > 0.0 && pperp > 0.0,
            "limiter stress produced a nonpositive primitive state");
    Require(!cgl::HardBoundViolated(paniso, bsqr, eos, eos.backup_lim),
            "limiter stress exceeded a configured hard wall");
    if (check_relaxation) {
      const Real phase = Wavenumber(pin, pm)*(XCenter(pm, q) - pm->mesh_size.x1min);
      const Real seed = (limiter_kind == "mirror" ? 1.0 : -1.0)*0.25*amp*std::sin(phase);
      const Real initial_ppar = ppar0*(1.0 + seed);
      const Real initial_pperp = pperp0*(1.0 - seed);
      const Real initial_piso = (initial_ppar + 2.0*initial_pperp)/3.0;
      const Real mirror = 0.5*eos.mirror_threshold*bsqr;
      const Real firehose = -0.5*eos.firehose_threshold*bsqr;
      const Real lower = eos.backup_lim ?
          std::max(-bsqr, eos.firehose_backup_factor*firehose) : -bsqr;
      const Real upper = eos.backup_lim ? eos.mirror_backup_factor*mirror :
          std::numeric_limits<Real>::max();
      Real expected = initial_pperp - initial_ppar;
      // The LF pre sweep and RK boundary apply walls before the final rates.
      if (has_lf) expected = std::min(upper, std::max(lower, expected));
      expected *= std::exp(-eos.nu_coll*pm->time);
      if (eos.mlim && expected > mirror) {
        expected = mirror + (expected - mirror)/(1.0 + eos.lim_coll*pm->time);
      } else if (eos.flim && expected < firehose) {
        expected = firehose + (expected - firehose)/(1.0 + eos.lim_coll*pm->time);
      }
      expected = std::min(upper, std::max(lower, expected));
      const Real error = std::abs(paniso - expected)/initial_piso;
      max_relaxation_error = std::max(max_relaxation_error, error);
      Require(error <= 1.0e-12,
              "limiter stress disagrees with one full-step analytic relaxation");
      RequireRelative("limiter stress conserved isotropic pressure",
                      (ppar + 2.0*pperp)/3.0, initial_piso, 1.0e-12);
    }
  }
  std::cout << "CGL LF limiter_stress passed for " << limiter_kind
            << " analytic_relaxation=" << check_relaxation
            << " max_relaxation_error=" << max_relaxation_error << std::endl;
}

struct WaveState {
  Complex rho;
  Complex vx;
  Complex ppar;
};

WaveState WaveRHS(const WaveState &s, const Real rho0, const Real ppar0,
                  const Real k_wave, const Real chi_parallel) {
  const Complex ik(0.0, k_wave);
  WaveState rhs;
  rhs.rho = -ik*rho0*s.vx;
  rhs.vx = -ik*s.ppar/rho0;
  rhs.ppar = -static_cast<Real>(3.0)*ik*ppar0*s.vx
             - chi_parallel*SQR(k_wave)*(s.ppar - (ppar0/rho0)*s.rho);
  return rhs;
}

WaveState operator+(const WaveState &a, const WaveState &b) {
  return {a.rho + b.rho, a.vx + b.vx, a.ppar + b.ppar};
}

WaveState operator*(const Real c, const WaveState &a) {
  return {c*a.rho, c*a.vx, c*a.ppar};
}

WaveState IntegrateWaveReference(ParameterInput *pin, Mesh *pm) {
  const Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  const Real ppar0 = pin->GetOrAddReal("problem", "ppar0", 1.0);
  const Real amp = pin->GetOrAddReal("problem", "amp", 1.0e-5);
  const Real k_wave = Wavenumber(pin, pm);
  const Real pperp0 = pin->GetOrAddReal("problem", "pperp0", ppar0);
  const Real bx0 = pin->GetOrAddReal("problem", "b0", 1.0);
  const Real by0 = pin->GetOrAddReal("problem", "by0", 0.0);
  const Real bz0 = pin->GetOrAddReal("problem", "bz0", 0.0);
  const Real cpar0 = FaceCParallel(pin, rho0, ppar0);
  const Real lf_k = pin->GetReal("mhd", "lf_k_parallel");
  const Real nu_eff = EffectiveCollisionFrequency(pin, pm, ppar0, pperp0, bx0, by0, bz0);
  const Real chi_parallel = ChiParallel(cpar0, lf_k, nu_eff);
  const Real c_cgl = std::sqrt(3.0*ppar0/rho0);

  WaveState s{Complex(0.0, -rho0*amp),
              Complex(0.0, -c_cgl*amp),
              Complex(0.0, -3.0*ppar0*amp)};
  const int nsteps = pin->GetOrAddInteger("problem", "reference_steps", 20000);
  const Real dt = pm->time/static_cast<Real>(std::max(nsteps, 1));
  for (int n = 0; n < nsteps; ++n) {
    const WaveState k1 = WaveRHS(s, rho0, ppar0, k_wave, chi_parallel);
    const WaveState k2 = WaveRHS(s + 0.5*dt*k1, rho0, ppar0, k_wave, chi_parallel);
    const WaveState k3 = WaveRHS(s + 0.5*dt*k2, rho0, ppar0, k_wave, chi_parallel);
    const WaveState k4 = WaveRHS(s + dt*k3, rho0, ppar0, k_wave, chi_parallel);
    s = s + (dt/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
  }
  return s;
}

void RequireComplexRelative(const std::string &label, const Complex got,
                            const Complex expected, const Real scale,
                            const Real rel_tol) {
  const Real rel_err =
      std::abs(got - expected)/std::max(scale, static_cast<Real>(1.0e-30));
  if (!std::isfinite(std::real(got)) || !std::isfinite(std::imag(got)) ||
      rel_err > rel_tol) {
    std::cout << label << " got=(" << std::real(got) << "," << std::imag(got)
              << ") expected=(" << std::real(expected) << "," << std::imag(expected)
              << ") rel_err=" << rel_err << " rel_tol=" << rel_tol << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

void RequireWaveEvolution(const Real reference_change, const Real tolerance) {
  Require(reference_change >= 10.0*tolerance,
          "wave reference changed by less than ten tolerances; "
          "a frozen state could pass");
}

void CheckFieldAlignedWave(ParameterInput *pin, Mesh *pm) {
  auto *pmhd = pm->pmb_pack->pmhd;
  auto w = HostCopy(pmhd->w0);
  const WaveState ref = IntegrateWaveReference(pin, pm);
  const Projection rho_p = ProjectPrimitive(w, pin, pm, IDN);
  const Projection vx_p = ProjectPrimitive(w, pin, pm, IVX);
  const Projection ppar_p = ProjectPrimitive(w, pin, pm, IPR);
  const Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  const Real ppar0 = pin->GetOrAddReal("problem", "ppar0", 1.0);
  const Real amp = pin->GetOrAddReal("problem", "amp", 1.0e-5);
  const Real wave_tol = pin->GetOrAddReal("problem", "wave_rel_tol", 1.0e-3);
  const Real c_cgl = std::sqrt(3.0*ppar0/rho0);
  const Complex rho_m = ProjectionToComplex(rho_p);
  const Complex vx_m = ProjectionToComplex(vx_p);
  const Complex ppar_m = ProjectionToComplex(ppar_p);

  RequireWaveEvolution(std::max({
      ComplexRelativeError(ref.rho, Complex(0.0, -rho0*amp), rho0*amp),
      ComplexRelativeError(ref.vx, Complex(0.0, -c_cgl*amp), c_cgl*amp),
      ComplexRelativeError(ref.ppar, Complex(0.0, -3.0*ppar0*amp), 3.0*ppar0*amp)}),
      wave_tol);
  RequireComplexRelative("field_aligned_wave rho", rho_m, ref.rho, rho0*amp, wave_tol);
  RequireComplexRelative("field_aligned_wave vx", vx_m, ref.vx, c_cgl*amp, wave_tol);
  RequireComplexRelative("field_aligned_wave p_parallel", ppar_m, ref.ppar,
                         3.0*ppar0*amp, wave_tol);
  auto out = OpenValidationDataFile(pin, "field_aligned_wave");
  if (out) {
    out << "variable,measured_real,measured_imag,reference_real,reference_imag,"
        << "scale,rel_err,rel_tol\n";
    out << "rho," << std::real(rho_m) << "," << std::imag(rho_m) << ","
        << std::real(ref.rho) << "," << std::imag(ref.rho) << ","
        << rho0*amp << "," << ComplexRelativeError(rho_m, ref.rho, rho0*amp)
        << "," << wave_tol << "\n";
    out << "vx," << std::real(vx_m) << "," << std::imag(vx_m) << ","
        << std::real(ref.vx) << "," << std::imag(ref.vx) << ","
        << c_cgl*amp << "," << ComplexRelativeError(vx_m, ref.vx, c_cgl*amp)
        << "," << wave_tol << "\n";
    out << "p_parallel," << std::real(ppar_m) << "," << std::imag(ppar_m) << ","
        << std::real(ref.ppar) << "," << std::imag(ref.ppar) << ","
        << 3.0*ppar0*amp << ","
        << ComplexRelativeError(ppar_m, ref.ppar, 3.0*ppar0*amp)
        << "," << wave_tol << "\n";
  }
  std::cout << "CGL LF field_aligned_wave passed against linear Fourier reference"
            << ": rho_rel_err=" << ComplexRelativeError(rho_m, ref.rho, rho0*amp)
            << " vx_rel_err=" << ComplexRelativeError(vx_m, ref.vx, c_cgl*amp)
            << " p_parallel_rel_err="
            << ComplexRelativeError(ppar_m, ref.ppar, 3.0*ppar0*amp)
            << " rel_tol=" << wave_tol << std::endl;
}

struct PaperWaveState {
  Complex rho;
  Complex vx;
  Complex vy;
  Complex vz;
  Complex by;
  Complex bz;
  Complex ppar;
  Complex pperp;
};

PaperWaveState operator+(const PaperWaveState &a, const PaperWaveState &b) {
  return {a.rho + b.rho, a.vx + b.vx, a.vy + b.vy, a.vz + b.vz,
          a.by + b.by, a.bz + b.bz, a.ppar + b.ppar, a.pperp + b.pperp};
}

PaperWaveState operator*(const Real c, const PaperWaveState &a) {
  return {c*a.rho, c*a.vx, c*a.vy, c*a.vz, c*a.by, c*a.bz,
          c*a.ppar, c*a.pperp};
}

PaperWaveState operator*(const Complex c, const PaperWaveState &a) {
  return {c*a.rho, c*a.vx, c*a.vy, c*a.vz, c*a.by, c*a.bz,
          c*a.ppar, c*a.pperp};
}

PaperWaveState PaperWaveRHS(const PaperWaveState &s, const Real rho0,
                            const Real p0, const Real bx0, const Real by0,
                            const Real bz0, const Real k_wave,
                            const Real chi_parallel, const Real chi_perp) {
  const Complex ik(0.0, k_wave);
  const Real bmag0 = std::sqrt(SQR(bx0) + SQR(by0) + SQR(bz0));
  const Real bhx = bx0/bmag0;
  const Real bhy = by0/bmag0;
  const Real bhz = bz0/bmag0;
  const Complex delta_p = s.ppar - s.pperp;
  const Complex flux_x = s.pperp + delta_p*bhx*bhx + by0*s.by + bz0*s.bz;
  const Complex flux_y = delta_p*bhx*bhy - bx0*s.by;
  const Complex flux_z = delta_p*bhx*bhz - bx0*s.bz;

  PaperWaveState rhs;
  rhs.rho = -ik*rho0*s.vx;
  rhs.vx = -ik*flux_x/rho0;
  rhs.vy = -ik*flux_y/rho0;
  rhs.vz = -ik*flux_z/rho0;
  rhs.by = -ik*(by0*s.vx - bx0*s.vy);
  rhs.bz = -ik*(bz0*s.vx - bx0*s.vz);

  const Complex dbmag_dt = bhy*rhs.by + bhz*rhs.bz;
  rhs.ppar = p0*(static_cast<Real>(3.0)*rhs.rho/rho0 -
                 static_cast<Real>(2.0)*dbmag_dt/bmag0);
  rhs.pperp = p0*(rhs.rho/rho0 + dbmag_dt/bmag0);

  const Complex q_parallel =
      -chi_parallel*ik*bhx*(s.ppar - (p0/rho0)*s.rho);
  const Complex q_perp =
      -chi_perp*ik*bhx*(s.pperp - (p0/rho0)*s.rho);
  rhs.ppar += -ik*bhx*q_parallel;
  rhs.pperp += -ik*bhx*q_perp;
  return rhs;
}

PaperWaveState IntegratePaperWaveReference(ParameterInput *pin, Mesh *pm) {
  const Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  const Real p0 = pin->GetOrAddReal("problem", "ppar0", 5.0);
  const Real amp = pin->GetOrAddReal("problem", "amp", 1.0e-5);
  const Real bx0 = pin->GetOrAddReal("problem", "b0", 1.0);
  const Real by0 = pin->GetOrAddReal("problem", "by0", std::sqrt(2.0));
  const Real bz0 = pin->GetOrAddReal("problem", "bz0", 0.5);
  const Real k_wave = Wavenumber(pin, pm);
  Real chi_parallel = 0.0;
  Real chi_perp = 0.0;
  if (pin->DoesParameterExist("mhd", "cgl_heat_flux")) {
    const Real cpar0 = FaceCParallel(pin, rho0, p0);
    const Real lf_k = pin->GetReal("mhd", "lf_k_parallel");
    const Real nu_eff = EffectiveCollisionFrequency(pin, pm, p0, p0, bx0, by0, bz0);
    chi_parallel = ChiParallel(cpar0, lf_k, nu_eff);
    chi_perp = ChiPerp(cpar0, lf_k, nu_eff);
  }

  // This deliberately follows the pgen's transverse-velocity initial value
  // problem.  It is not an exact CGL eigenmode initialization.
  PaperWaveState s{Complex(0.0, 0.0),
                   Complex(0.0, 0.0),
                   Complex(0.0, -amp),
                   Complex(0.0, 0.0),
                   Complex(0.0, 0.0),
                   Complex(0.0, 0.0),
                   Complex(0.0, 0.0),
                   Complex(0.0, 0.0)};
  const int nsteps = pin->GetOrAddInteger("problem", "reference_steps", 20000);
  const Real dt = pm->time/static_cast<Real>(std::max(nsteps, 1));
  for (int n = 0; n < nsteps; ++n) {
    const PaperWaveState k1 = PaperWaveRHS(s, rho0, p0, bx0, by0, bz0, k_wave,
                                           chi_parallel, chi_perp);
    const PaperWaveState k2 = PaperWaveRHS(s + 0.5*dt*k1, rho0, p0, bx0, by0,
                                           bz0, k_wave, chi_parallel, chi_perp);
    const PaperWaveState k3 = PaperWaveRHS(s + 0.5*dt*k2, rho0, p0, bx0, by0,
                                           bz0, k_wave, chi_parallel, chi_perp);
    const PaperWaveState k4 = PaperWaveRHS(s + dt*k3, rho0, p0, bx0, by0, bz0,
                                           k_wave, chi_parallel, chi_perp);
    s = s + (dt/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
  }
  return s;
}

Complex GetComplexParameter(ParameterInput *pin, const std::string &base) {
  return Complex(pin->GetOrAddReal("problem", base + "_re", 0.0),
                 pin->GetOrAddReal("problem", base + "_im", 0.0));
}

PaperWaveState ReadPaperEigenVector(ParameterInput *pin) {
  return {GetComplexParameter(pin, "eigen_rho"),
          GetComplexParameter(pin, "eigen_vx"),
          GetComplexParameter(pin, "eigen_vy"),
          GetComplexParameter(pin, "eigen_vz"),
          GetComplexParameter(pin, "eigen_by"),
          GetComplexParameter(pin, "eigen_bz"),
          GetComplexParameter(pin, "eigen_ppar"),
          GetComplexParameter(pin, "eigen_pperp")};
}

void CheckPaperObliqueWave(ParameterInput *pin, Mesh *pm) {
  auto *pmhd = pm->pmb_pack->pmhd;
  auto w = HostCopy(pmhd->w0);
  auto bcc = HostCopy(pmhd->bcc0);
  const PaperWaveState ref = IntegratePaperWaveReference(pin, pm);
  const Real amp = pin->GetOrAddReal("problem", "amp", 1.0e-5);
  const Real p0 = pin->GetOrAddReal("problem", "ppar0", 5.0);
  const Real wave_tol = pin->GetOrAddReal("problem", "wave_rel_tol", 1.0e-3);
  const Complex vy_m = ProjectionToComplex(ProjectPrimitive(w, pin, pm, IVY));
  const Complex by_m = ProjectionToComplex(ProjectCellField(bcc, pin, pm, IBY));
  const Complex ppar_m = ProjectionToComplex(ProjectPrimitive(w, pin, pm, IPR));
  const Complex pperp_m = ProjectionToComplex(ProjectPrimitive(w, pin, pm, IPP));

  RequireWaveEvolution(std::max({
      ComplexRelativeError(ref.vy, Complex(0.0, -amp), amp),
      std::abs(ref.by)/amp, std::abs(ref.ppar)/(p0*amp),
      std::abs(ref.pperp)/(p0*amp)}), wave_tol);
  RequireComplexRelative("paper_oblique_wave vy", vy_m, ref.vy, amp, wave_tol);
  RequireComplexRelative("paper_oblique_wave By", by_m, ref.by, amp, wave_tol);
  RequireComplexRelative("paper_oblique_wave p_parallel", ppar_m, ref.ppar,
                         p0*amp, wave_tol);
  RequireComplexRelative("paper_oblique_wave p_perp", pperp_m, ref.pperp,
                         p0*amp, wave_tol);
  auto out = OpenValidationDataFile(pin, "paper_oblique_wave");
  if (out) {
    out << "variable,measured_real,measured_imag,reference_real,reference_imag,"
        << "scale,rel_err,rel_tol\n";
    out << "vy," << std::real(vy_m) << "," << std::imag(vy_m) << ","
        << std::real(ref.vy) << "," << std::imag(ref.vy) << ","
        << amp << "," << ComplexRelativeError(vy_m, ref.vy, amp)
        << "," << wave_tol << "\n";
    out << "By," << std::real(by_m) << "," << std::imag(by_m) << ","
        << std::real(ref.by) << "," << std::imag(ref.by) << ","
        << amp << "," << ComplexRelativeError(by_m, ref.by, amp)
        << "," << wave_tol << "\n";
    out << "p_parallel," << std::real(ppar_m) << "," << std::imag(ppar_m) << ","
        << std::real(ref.ppar) << "," << std::imag(ref.ppar) << ","
        << p0*amp << "," << ComplexRelativeError(ppar_m, ref.ppar, p0*amp)
        << "," << wave_tol << "\n";
    out << "p_perp," << std::real(pperp_m) << "," << std::imag(pperp_m) << ","
        << std::real(ref.pperp) << "," << std::imag(ref.pperp) << ","
        << p0*amp << "," << ComplexRelativeError(pperp_m, ref.pperp, p0*amp)
        << "," << wave_tol << "\n";
  }
  const char *closure = pin->DoesParameterExist("mhd", "cgl_heat_flux") ?
                        "CGL LF" : "pure CGL";
  std::cout << closure
            << " paper_oblique_wave passed against linear oblique IVP reference"
            << ": vy_rel_err=" << ComplexRelativeError(vy_m, ref.vy, amp)
            << " By_rel_err=" << ComplexRelativeError(by_m, ref.by, amp)
            << " p_parallel_rel_err=" << ComplexRelativeError(ppar_m, ref.ppar, p0*amp)
            << " p_perp_rel_err=" << ComplexRelativeError(pperp_m, ref.pperp, p0*amp)
            << " rel_tol=" << wave_tol << std::endl;
}

void CheckPaperEigenWave(ParameterInput *pin, Mesh *pm) {
  auto *pmhd = pm->pmb_pack->pmhd;
  auto w = HostCopy(pmhd->w0);
  auto bcc = HostCopy(pmhd->bcc0);
  const PaperWaveState eigen = ReadPaperEigenVector(pin);
  const Real amp = pin->GetOrAddReal("problem", "amp", 1.0e-5);
  const Real wave_tol = pin->GetOrAddReal("problem", "eigen_wave_rel_tol", 1.0e-3);
  const Real component_floor =
      pin->GetOrAddReal("problem", "eigen_component_floor", 1.0e-4);
  const Real zero_abs_tol =
      pin->GetOrAddReal("problem", "eigen_zero_abs_tol", 1.0e-8);
  const Complex lambda(pin->GetReal("problem", "eigen_lambda_re"),
                       pin->GetReal("problem", "eigen_lambda_im"));
  const Complex phase = std::exp(lambda*pm->time);
  const PaperWaveState ref = (amp*phase)*eigen;
  RequireWaveEvolution(std::abs(phase - Complex(1.0, 0.0)), wave_tol);

  const PaperWaveState measured{
      ProjectionToComplex(ProjectPrimitive(w, pin, pm, IDN)),
      ProjectionToComplex(ProjectPrimitive(w, pin, pm, IVX)),
      ProjectionToComplex(ProjectPrimitive(w, pin, pm, IVY)),
      ProjectionToComplex(ProjectPrimitive(w, pin, pm, IVZ)),
      ProjectionToComplex(ProjectCellField(bcc, pin, pm, IBY)),
      ProjectionToComplex(ProjectCellField(bcc, pin, pm, IBZ)),
      ProjectionToComplex(ProjectPrimitive(w, pin, pm, IPR)),
      ProjectionToComplex(ProjectPrimitive(w, pin, pm, IPP))};

  const std::string branch = pin->GetOrAddString("problem", "eigen_branch", "unknown");
  auto out = OpenValidationDataFile(pin, "paper_eigen_wave");
  if (out) {
    out << "variable,measured_real,measured_imag,reference_real,reference_imag,"
        << "eigen_real,eigen_imag,scale,rel_err,rel_tol,required\n";
  }

  const char *names[8] = {"rho", "vx", "vy", "vz", "By", "Bz", "p_parallel",
                          "p_perp"};
  const Complex measured_components[8] = {measured.rho, measured.vx, measured.vy,
                                          measured.vz, measured.by, measured.bz,
                                          measured.ppar, measured.pperp};
  const Complex ref_components[8] = {ref.rho, ref.vx, ref.vy, ref.vz, ref.by, ref.bz,
                                     ref.ppar, ref.pperp};
  const Complex eigen_components[8] = {eigen.rho, eigen.vx, eigen.vy, eigen.vz,
                                       eigen.by, eigen.bz, eigen.ppar, eigen.pperp};
  Real max_rel_err = 0.0;
  Real max_abs_err = 0.0;
  for (int n = 0; n < 8; ++n) {
    const Real component_amp = std::abs(eigen_components[n]);
    const bool required = component_amp >= component_floor;
    const Real scale = amp*std::max(component_amp, component_floor);
    const Real rel_err = ComplexRelativeError(measured_components[n],
                                              ref_components[n], scale);
    const Real abs_err = std::abs(measured_components[n] - ref_components[n]);
    max_rel_err = std::max(max_rel_err, rel_err);
    max_abs_err = std::max(max_abs_err, abs_err);
    if (out) {
      out << names[n] << "," << std::real(measured_components[n]) << ","
          << std::imag(measured_components[n]) << ","
          << std::real(ref_components[n]) << "," << std::imag(ref_components[n]) << ","
          << std::real(eigen_components[n]) << "," << std::imag(eigen_components[n])
          << "," << scale << "," << rel_err << "," << wave_tol << ","
          << (required ? 1 : 0) << "\n";
    }
    if (required) {
      RequireComplexRelative(std::string("paper_eigen_wave ") + names[n],
                             measured_components[n], ref_components[n], scale,
                             wave_tol);
    } else {
      Require(abs_err <= zero_abs_tol,
              std::string("paper_eigen_wave inactive component ") + names[n] +
              " exceeded zero_abs_tol");
    }
  }

  std::cout << "CGL paper_eigen_wave " << branch
            << " passed against supplied eigenmode: lambda=("
            << std::real(lambda) << "," << std::imag(lambda)
            << ") max_rel_err=" << max_rel_err
            << " rel_tol=" << wave_tol
            << " max_abs_err=" << max_abs_err
            << " zero_abs_tol=" << zero_abs_tol << std::endl;
}

Real hotspot_min_parallel, hotspot_min_perp, hotspot_initial_energy;
Real hotspot_initial_parallel, hotspot_initial_perp, hotspot_max_energy_error;

void MonitorHotSpot(Mesh *pm, const Real) {
  const auto w = HostCopy(pm->pmb_pack->pmhd->w0);
  const auto u = HostCopy(pm->pmb_pack->pmhd->u0);
  const auto &ind = pm->mb_indcs;
  Real energy = 0.0;
  for (int j=ind.js; j<=ind.je; ++j) {
    for (int i=ind.is; i<=ind.ie; ++i) {
      const Real tpar = w(0,IPR,ind.ks,j,i)/w(0,IDN,ind.ks,j,i);
      const Real tperp = w(0,IPP,ind.ks,j,i)/w(0,IDN,ind.ks,j,i);
      Require(std::isfinite(tpar) && std::isfinite(tperp), "hotspot became nonfinite");
      hotspot_min_parallel = std::min(hotspot_min_parallel, tpar);
      hotspot_min_perp = std::min(hotspot_min_perp, tperp);
      energy += u(0,IEN,ind.ks,j,i);
    }
  }
  if (hotspot_initial_energy == 0.0) hotspot_initial_energy = energy;
  hotspot_max_energy_error = std::max(hotspot_max_energy_error,
      std::abs(energy/hotspot_initial_energy - 1.0));
}

void HotSpotHistory(HistoryData *pdata, Mesh *pm) {
  MonitorHotSpot(pm, 0.0);
  pdata->nhist = 3;
  pdata->label[0] = "min_tpar";
  pdata->label[1] = "min_tperp";
  pdata->label[2] = "max_energy_error";
  pdata->hdata[0] = hotspot_min_parallel;
  pdata->hdata[1] = hotspot_min_perp;
  pdata->hdata[2] = hotspot_max_energy_error;
}

void CheckHotSpot(ParameterInput *pin, Mesh *pm) {
  MonitorHotSpot(pm, 0.0);
  const Real width = pin->GetReal("problem", "hotspot_width");
  const Real cpar = pin->GetReal("mhd", "lf_c_parallel0");
  const Real chi_perp = std::sqrt(2.0/std::acos(-1.0))*cpar/
                       pin->GetReal("mhd", "lf_k_parallel");
  const Real diffusion_times = pm->time*chi_perp/SQR(width);
  std::cout << std::setprecision(17)
            << "CGL LF hotspot: min_tpar=" << hotspot_min_parallel
            << " min_tperp=" << hotspot_min_perp
            << " energy_error=" << hotspot_max_energy_error
            << " diffusion_times=" << diffusion_times << std::endl;
  Require(diffusion_times >= 5.0, "hotspot requires at least five diffusion times");
  Require(hotspot_min_parallel >= hotspot_initial_parallel*(1.0 - 1.0e-12) &&
          hotspot_min_perp >= hotspot_initial_perp*(1.0 - 1.0e-12),
          "hotspot fell below its initial temperature minima");
  Require(hotspot_max_energy_error <= 5.0e-13, "hotspot did not conserve total energy");
}

Real reversal_peak_deviation = 0.0;
Real reversal_min_pressure = 1.0;
Real reversal_ppar0 = 1.0, reversal_pperp0 = 1.0;

// Sample both the pre-sweep state at each RK stage and the final post-sweep state.
void MonitorFieldReversal(Mesh *pm, const Real) {
  const auto w = HostCopy(pm->pmb_pack->pmhd->w0);
  const auto &indcs = pm->mb_indcs;
  for (int m=0; m<pm->pmb_pack->nmb_thispack; ++m) {
    for (int k=indcs.ks; k<=indcs.ke; ++k) {
      for (int j=indcs.js; j<=indcs.je; ++j) {
        for (int i=indcs.is; i<=indcs.ie; ++i) {
          const Real ppar = w(m,IPR,k,j,i), pperp = w(m,IPP,k,j,i);
          Require(std::isfinite(ppar) && std::isfinite(pperp),
                  "field reversal developed nonfinite pressures");
          reversal_min_pressure = std::min(reversal_min_pressure, std::min(ppar, pperp));
          reversal_peak_deviation = std::max(reversal_peak_deviation,
              std::max(std::abs(ppar/reversal_ppar0 - 1.0),
                       std::abs(pperp/reversal_pperp0 - 1.0)));
        }
      }
    }
  }
}

void CheckFieldReversal(ParameterInput *pin, Mesh *pm) {
  MonitorFieldReversal(pm, 0.0);
  const Real amp = pin->GetReal("problem", "amp");
  const Real growth = reversal_peak_deviation/std::max(amp, static_cast<Real>(1.0e-12));
  const Real sweep_ratio = 0.5*pm->dt/pm->dt_parabolic_sts;
  std::cout << "CGL LF field_reversal: cycles=" << pm->ncycle
            << " sweep_ratio=" << sweep_ratio
            << " peak_delta_p=" << reversal_peak_deviation
            << " growth=" << growth << " min_pressure=" << reversal_min_pressure
            << std::endl;
  Require(pm->ncycle >= 20, "field reversal requires at least 20 cycles");
  Require(sweep_ratio >= 9.5, "field reversal requires a sweep ratio near 10 or larger");
  Require(reversal_min_pressure > 0.0, "field reversal developed nonpositive pressures");
  Require(reversal_peak_deviation < 1.0e-2 && growth < 10.0,
          "field reversal amplified the pressure seed");
}

// A uniform RK source isolates post-sweep timestep refresh from spatial transport.
Real refresh_heating_rate = 0.0;
Real refresh_pre_dt = 0.0;
int refresh_pre_stages = 0;

void HeatTimestepRefresh(Mesh *pm, const Real beta_dt) {
  auto *pmhd = pm->pmb_pack->pmhd;
  const auto &indcs = pm->mb_indcs;
  const int nmb = pm->pmb_pack->nmb_thispack;
  const int ncell = nmb*indcs.nx1*indcs.nx2*indcs.nx3;
  if (refresh_pre_stages == 0) {
    refresh_pre_stages = pmhd->pcgl_lf->diagnostics.nstage/ncell;
    refresh_pre_dt = pmhd->pcgl_lf->dtnew;
  }
  auto u = pmhd->u0;
  const Real de = refresh_heating_rate*beta_dt;
  par_for("cgl_lf_refresh_heating", DevExeSpace(), 0, nmb - 1,
          indcs.ks, indcs.ke, indcs.js, indcs.je, indcs.is, indcs.ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    u(m,IEN,k,j,i) += de;
  });
}

void CheckTimestepRefresh(ParameterInput *pin, Mesh *pm) {
  auto *lf = pm->pmb_pack->pmhd->pcgl_lf;
  const auto &indcs = pm->mb_indcs;
  const int nmb = pm->pmb_pack->nmb_thispack;
  const int ncell = nmb*indcs.nx1*indcs.nx2*indcs.nx3;
  const int post_stages = lf->diagnostics.nstage/ncell - refresh_pre_stages;
  int stages_min[2] = {refresh_pre_stages, post_stages};
  int stages_max[2] = {refresh_pre_stages, post_stages};
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, stages_min, 2, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, stages_max, 2, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif
  Require(stages_min[0] == stages_max[0] && stages_min[1] == stages_max[1],
          "pre/post stage counts differ across MPI ranks");
  Require(pm->ncycle == 1 && refresh_pre_stages == 7,
          "timestep refresh requires one cycle with seven pre-sweep stages");
  Require(refresh_heating_rate > 0.0 ? post_stages > refresh_pre_stages :
                                      post_stages == refresh_pre_stages,
          "post-sweep stages did not respond to the RK heating");
  const Real pressure = 1.0 + (2.0/3.0)*refresh_heating_rate*pm->dt_last_completed;
  const Real dx = (pm->mesh_size.x1max - pm->mesh_size.x1min)/pm->mesh_indcs.nx1;
  const Real chi0 = std::sqrt(8.0/M_PI)/lf->lf_k_parallel;
  const Real dt0 = 0.5*dx*dx/chi0;
  RequireRelative("initial uniform LF timestep", refresh_pre_dt, dt0, 2.0e-12);
  RequireRelative("heated uniform LF timestep", lf->dtnew, dt0/std::sqrt(pressure),
                  2.0e-12);
  RequireRelative("refreshed mesh LF timestep", pm->dt_parabolic_sts,
                  pm->cfl_no*dt0/std::sqrt(pressure), 2.0e-12);
  const auto w = HostCopy(pm->pmb_pack->pmhd->w0);
  for (int m=0; m<nmb; ++m) {
    for (int i=indcs.is; i<=indcs.ie; ++i) {
      RequireRelative("heated parallel pressure", w(m,IPR,indcs.ks,indcs.js,i),
                      pressure, 2.0e-12);
      RequireRelative("heated perpendicular pressure", w(m,IPP,indcs.ks,indcs.js,i),
                      pressure, 2.0e-12);
    }
  }
  if (global_variable::my_rank == 0) {
    std::cout << "CGL LF timestep_refresh: pre_stages=" << stages_min[0]
              << " post_stages=" << stages_min[1] << " ranks_agree=true"
              << " final_pressure=" << pressure << " dt_parabolic="
              << pm->dt_parabolic_sts << std::endl;
  }
}

void CheckDensityContact(ParameterInput *pin, Mesh *pm) {
  const auto w = HostCopy(pm->pmb_pack->pmhd->w0);
  const auto &indcs = pm->mb_indcs;
  const Real tpar0 = pin->GetReal("problem", "ppar0")/pin->GetReal("problem", "rho0");
  const Real tperp0 = pin->GetReal("problem", "pperp0")/pin->GetReal("problem", "rho0");
  Real deviation = 0.0;
  for (int i=indcs.is; i<=indcs.ie; ++i) {
    const Real rho = w(0,IDN,indcs.ks,indcs.js,i);
    const Real tpar = w(0,IPR,indcs.ks,indcs.js,i)/rho;
    const Real tperp = w(0,IPP,indcs.ks,indcs.js,i)/rho;
    Require(std::isfinite(tpar) && std::isfinite(tperp) && tpar > 0.0 && tperp > 0.0,
            "density contact developed invalid temperatures");
    deviation = std::max(deviation,
        std::max(std::abs(tpar/tpar0 - 1.0), std::abs(tperp/tperp0 - 1.0)));
  }
  std::cout << "CGL LF density_contact: contrast="
            << pin->GetReal("problem", "density_contrast")
            << " cycles=" << pm->ncycle << " max_delta_T=" << deviation
            << " sweep_ratio=" << 0.5*pm->dt/pm->dt_parabolic_sts << std::endl;
  Require(pm->ncycle >= 3, "density contact requires at least three cycles");
  Require(0.5*pm->dt/pm->dt_parabolic_sts >= 10.0 - 1.0e-12,
          "density contact requires sweep ratio >= 10");
  // The regression also compares seeded/unseeded runs at every cycle.
  const Real amp = pin->GetReal("problem", "amp");
  Require(deviation < 10.0*std::max(amp, static_cast<Real>(1.0e-12)),
          "density contact amplified its temperature seed");
}

void FinalizeCGLLFQuantitative(ParameterInput *pin, Mesh *pm) {
  auto *pmhd = pm->pmb_pack->pmhd;
  Require(pmhd != nullptr && pmhd->peos->eos_data.is_cgl,
          "quantitative LF tests require <mhd>/eos = cgl");
  const TestMode mode = ParseMode(pin);
  if (mode == TestMode::rotated_decay || mode == TestMode::field_reversal ||
      mode == TestMode::hotspot) {
    RequireSingleBlock(pm);
  } else if (mode == TestMode::field_aligned_wave ||
             mode == TestMode::paper_oblique_wave ||
             mode == TestMode::paper_eigen_wave ||
             mode == TestMode::timestep_refresh) {
    RequireOneDimensionalMesh(pm);
  } else {
    RequireOneDimensionalSingleBlock(pm);
  }
  if (mode == TestMode::parallel_decay || mode == TestMode::perp_decay) {
    CheckDecay(pin, pm, mode);
  } else if (mode == TestMode::collision_relaxation) {
    CheckCollisionRelaxation(pin, pm);
  } else if (mode == TestMode::grad_b) {
    CheckGradB(pin, pm);
  } else if (mode == TestMode::flux_limiter) {
    CheckFluxLimiter(pin, pm);
  } else if (mode == TestMode::limiter_heat_flux_suppression) {
    CheckLimiterHeatFluxSuppression(pin, pm);
  } else if (mode == TestMode::limiter_stress) {
    CheckLimiterStress(pin, pm);
  } else if (mode == TestMode::field_aligned_wave) {
    CheckFieldAlignedWave(pin, pm);
  } else if (mode == TestMode::paper_oblique_wave) {
    CheckPaperObliqueWave(pin, pm);
  } else if (mode == TestMode::paper_eigen_wave) {
    CheckPaperEigenWave(pin, pm);
  } else if (mode == TestMode::rotated_decay) {
    CheckRotatedDecay(pin, pm);
  } else if (mode == TestMode::timestep_refresh) {
    CheckTimestepRefresh(pin, pm);
  } else if (mode == TestMode::density_contact) {
    CheckDensityContact(pin, pm);
  } else if (mode == TestMode::field_reversal) {
    CheckFieldReversal(pin, pm);
  } else if (mode == TestMode::hotspot) {
    CheckHotSpot(pin, pm);
  } else if (mode == TestMode::low_field) {
    CheckLowField(pin, pm);
  }
}

} // namespace

void ProblemGenerator::CGLLandauFluid(ParameterInput *pin, const bool restart) {
  pgen_final_func = FinalizeCGLLFQuantitative;
  if (restart) return;

  auto *pmbp = pmy_mesh_->pmb_pack;
  auto *pmhd = pmbp->pmhd;
  if (pmhd == nullptr || !pmhd->peos->eos_data.is_cgl) {
    Fail("quantitative LF tests require <mhd>/eos = cgl");
  }
  const TestMode mode = ParseMode(pin);
  if (mode == TestMode::rotated_decay || mode == TestMode::field_reversal ||
      mode == TestMode::hotspot) {
    RequireSingleBlock(pmy_mesh_);
  } else if (mode == TestMode::field_aligned_wave ||
             mode == TestMode::paper_oblique_wave ||
             mode == TestMode::paper_eigen_wave ||
             mode == TestMode::timestep_refresh) {
    RequireOneDimensionalMesh(pmy_mesh_);
  } else {
    RequireOneDimensionalSingleBlock(pmy_mesh_);
  }
  const Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  const Real ppar0 = pin->GetOrAddReal("problem", "ppar0", 1.0);
  const Real pperp0 = pin->GetOrAddReal("problem", "pperp0", 1.0);
  const Real amp = pin->GetOrAddReal("problem", "amp", 1.0e-4);
  const Real density_contrast = (mode == TestMode::density_contact) ?
      pin->GetOrAddReal("problem", "density_contrast", 200.0) : 1.0;
  if (mode == TestMode::field_reversal) {
    Require(user_srcs, "field_reversal requires <problem>/user_srcs = true");
    user_srcs_func = MonitorFieldReversal;
    reversal_peak_deviation = 0.0;
    reversal_min_pressure = std::min(ppar0, pperp0);
    reversal_ppar0 = ppar0;
    reversal_pperp0 = pperp0;
  }
  const Real hotspot_width = pin->GetOrAddReal("problem", "hotspot_width", 0.04);
  if (mode == TestMode::hotspot) {
    Require(pmy_mesh_->multi_d && !pmy_mesh_->three_d && user_srcs && user_hist,
            "hotspot requires 2D and user_srcs/user_hist = true");
    user_srcs_func = MonitorHotSpot;
    user_hist_func = HotSpotHistory;
    hotspot_min_parallel = hotspot_min_perp = 1.0e30;
    hotspot_initial_energy = hotspot_max_energy_error = 0.0;
  }
  if (mode == TestMode::timestep_refresh) {
    Require(user_srcs && rho0 == 1.0 && ppar0 == 1.0 && pperp0 == 1.0 &&
            pmhd->pcgl_lf != nullptr && pmhd->pcgl_lf->lf_coeff_local &&
            pmhd->peos->eos_data.nu_coll == 0.0,
            "timestep refresh requires unit uniform state, local closure, "
            "and no collisions");
    refresh_heating_rate = pin->GetOrAddReal("problem", "heating_rate", 3000.0);
    Require(refresh_heating_rate >= 0.0, "timestep refresh requires nonnegative heating");
    refresh_pre_stages = 0;
    refresh_pre_dt = 0.0;
    user_srcs_func = HeatTimestepRefresh;
  }
  const Real bx0 = pin->GetOrAddReal("problem", "b0", 1.0);
  const Real by0 = pin->GetOrAddReal("problem", "by0", 0.0);
  const Real bz0 = pin->GetOrAddReal("problem", "bz0", 0.0);
  const Real by_amp = pin->GetOrAddReal("problem", "by_amp", 0.0);
  const Real eig_rho_re = pin->GetOrAddReal("problem", "eigen_rho_re", 0.0);
  const Real eig_rho_im = pin->GetOrAddReal("problem", "eigen_rho_im", 0.0);
  const Real eig_vx_re = pin->GetOrAddReal("problem", "eigen_vx_re", 0.0);
  const Real eig_vx_im = pin->GetOrAddReal("problem", "eigen_vx_im", 0.0);
  const Real eig_vy_re = pin->GetOrAddReal("problem", "eigen_vy_re", 0.0);
  const Real eig_vy_im = pin->GetOrAddReal("problem", "eigen_vy_im", 0.0);
  const Real eig_vz_re = pin->GetOrAddReal("problem", "eigen_vz_re", 0.0);
  const Real eig_vz_im = pin->GetOrAddReal("problem", "eigen_vz_im", 0.0);
  const Real eig_by_re = pin->GetOrAddReal("problem", "eigen_by_re", 0.0);
  const Real eig_by_im = pin->GetOrAddReal("problem", "eigen_by_im", 0.0);
  const Real eig_bz_re = pin->GetOrAddReal("problem", "eigen_bz_re", 0.0);
  const Real eig_bz_im = pin->GetOrAddReal("problem", "eigen_bz_im", 0.0);
  const Real eig_ppar_re = pin->GetOrAddReal("problem", "eigen_ppar_re", 0.0);
  const Real eig_ppar_im = pin->GetOrAddReal("problem", "eigen_ppar_im", 0.0);
  const Real eig_pperp_re = pin->GetOrAddReal("problem", "eigen_pperp_re", 0.0);
  const Real eig_pperp_im = pin->GetOrAddReal("problem", "eigen_pperp_im", 0.0);
  const Real k_wave = Wavenumber(pin, pmy_mesh_);
  RotatedWave rotated_wave;
  if (mode == TestMode::rotated_decay) {
    rotated_wave = RotatedWavenumber(pin, pmy_mesh_);
  }
  const Real xmin = pmy_mesh_->mesh_size.x1min;
  const Real xlength = pmy_mesh_->mesh_size.x1max - xmin;
  const Real reversal_center = 0.5*(xmin + pmy_mesh_->mesh_size.x1max);
  const Real ymin = pmy_mesh_->mesh_size.x2min;
  const Real ymax = pmy_mesh_->mesh_size.x2max;
  const Real hotspot_x = 0.5*(xmin + pmy_mesh_->mesh_size.x1max) +
      0.5*(pmy_mesh_->mesh_size.x1max - xmin)/pmy_mesh_->mesh_indcs.nx1;
  const Real hotspot_y = 0.5*(ymin + ymax) + 0.5*(ymax - ymin)/pmy_mesh_->mesh_indcs.nx2;
  const Real zmin = pmy_mesh_->mesh_size.x3min;
  const Real zmax = pmy_mesh_->mesh_size.x3max;
  const std::string limiter_kind =
      pin->GetOrAddString("problem", "limiter_kind", "mirror");
  const int limiter_kind_id = (limiter_kind == "mirror") ? 0 :
                              (limiter_kind == "firehose") ? 1 : -1;
  if (limiter_kind_id < 0) {
    Fail("<problem>/limiter_kind must be mirror or firehose");
  }

  auto &indcs = pmy_mesh_->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmbp->nmb_thispack;
  auto w0 = pmhd->w0;
  auto bcc0 = pmhd->bcc0;
  auto b0 = pmhd->b0;
  auto size = pmbp->pmb->mb_size;

  par_for("cgl_lf_quant_init_prim", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int q = i - is;
    const RegionSize block_size = size.d_view(m);
    const Real x = CellCenterX(q, indcs.nx1, block_size.x1min, block_size.x1max);
    const Real y = CellCenterX(j - js, indcs.nx2, block_size.x2min,
                               block_size.x2max);
    const Real z = CellCenterX(k - ks, indcs.nx3, block_size.x3min,
                               block_size.x3max);
    Real phase = k_wave*(x - xmin);
    if (mode == TestMode::rotated_decay) {
      phase = rotated_wave.kx*(x - xmin) + rotated_wave.ky*(y - ymin)
            + rotated_wave.kz*(z - zmin);
    }
    const Real s = sin(phase);
    const Real c = cos(phase);
    Real rho = rho0;
    Real vx = 0.0;
    Real vy = 0.0;
    Real vz = 0.0;
    Real ppar = ppar0;
    Real pperp = pperp0;
    Real by = (mode == TestMode::grad_b) ? by_amp*s : by0;
    Real bz = bz0;

    if (mode == TestMode::density_contact) {
      // Two one-face density jumps, with uniform temperatures plus a small seed.
      if (x >= xmin + 0.25*xlength && x < xmin + 0.75*xlength) rho *= density_contrast;
      const Real seed = amp*((q%2 == 0) ? 1.0 : -1.0);
      ppar = ppar0*(rho/rho0)*(1.0 + seed);
      pperp = pperp0*(rho/rho0)*(1.0 + seed);
    } else if (mode == TestMode::hotspot) {
      const Real gaussian = exp(-(SQR(x - hotspot_x) + SQR(y - hotspot_y))/
                                  (2.0*SQR(hotspot_width)));
      ppar = ppar0*(1.0 + 99.0*gaussian);
      pperp = pperp0*(1.0 + 99.0*gaussian);
    } else if (mode == TestMode::field_reversal) {
      by = tanh((x - reversal_center)/block_size.dx1);
      const Real seed = amp*((q%2 == 0) ? 1.0 : -1.0);
      ppar = ppar0*(1.0 + seed);
      pperp = pperp0*(1.0 + seed);
    } else if (mode == TestMode::parallel_decay || mode == TestMode::rotated_decay) {
      ppar = ppar0*(1.0 + amp*s);
    } else if (mode == TestMode::perp_decay) {
      pperp = pperp0*(1.0 + amp*s);
    } else if (mode == TestMode::flux_limiter || mode == TestMode::low_field ||
               mode == TestMode::limiter_heat_flux_suppression) {
      ppar = ppar0*(1.0 + amp*s);
      pperp = pperp0*(1.0 + amp*s);
    } else if (mode == TestMode::limiter_stress) {
      if (limiter_kind_id == 0) {
        ppar = ppar0*(1.0 + 0.25*amp*s);
        pperp = pperp0*(1.0 - 0.25*amp*s);
      } else {
        ppar = ppar0*(1.0 - 0.25*amp*s);
        pperp = pperp0*(1.0 + 0.25*amp*s);
      }
    } else if (mode == TestMode::field_aligned_wave) {
      const Real c_cgl = sqrt(3.0*ppar0/rho0);
      rho = rho0*(1.0 + amp*s);
      vx = c_cgl*amp*s;
      ppar = ppar0*(1.0 + 3.0*amp*s);
      pperp = pperp0*(1.0 + amp*s);
    } else if (mode == TestMode::paper_oblique_wave) {
      vx = 0.0;
      vy = amp*s;
    } else if (mode == TestMode::paper_eigen_wave) {
      rho = rho0 + EigenRealSpacePerturbation(amp, eig_rho_re, eig_rho_im, c, s);
      vx = EigenRealSpacePerturbation(amp, eig_vx_re, eig_vx_im, c, s);
      vy = EigenRealSpacePerturbation(amp, eig_vy_re, eig_vy_im, c, s);
      vz = EigenRealSpacePerturbation(amp, eig_vz_re, eig_vz_im, c, s);
      ppar = ppar0 + EigenRealSpacePerturbation(amp, eig_ppar_re, eig_ppar_im, c, s);
      pperp = pperp0 + EigenRealSpacePerturbation(amp, eig_pperp_re, eig_pperp_im,
                                                  c, s);
      by = by0 + EigenRealSpacePerturbation(amp, eig_by_re, eig_by_im, c, s);
      bz = bz0 + EigenRealSpacePerturbation(amp, eig_bz_re, eig_bz_im, c, s);
    }

    w0(m,IDN,k,j,i) = rho;
    w0(m,IVX,k,j,i) = vx;
    w0(m,IVY,k,j,i) = vy;
    w0(m,IVZ,k,j,i) = vz;
    w0(m,IPR,k,j,i) = ppar;
    w0(m,IPP,k,j,i) = pperp;
    bcc0(m,IBX,k,j,i) = bx0;
    bcc0(m,IBY,k,j,i) = by;
    bcc0(m,IBZ,k,j,i) = bz;
  });

  par_for("cgl_lf_quant_init_b1", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie + 1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    b0.x1f(m,k,j,i) = bx0;
  });
  par_for("cgl_lf_quant_init_b2", DevExeSpace(), 0, nmb - 1, ks, ke, js, je + 1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int q = i - is;
    const RegionSize block_size = size.d_view(m);
    const Real x = CellCenterX(q, indcs.nx1, block_size.x1min, block_size.x1max);
    const Real s = sin(k_wave*(x - xmin));
    const Real c = cos(k_wave*(x - xmin));
    b0.x2f(m,k,j,i) = (mode == TestMode::field_reversal) ?
        tanh((x - reversal_center)/block_size.dx1) :
        (mode == TestMode::grad_b) ? by_amp*s :
        by0 + ((mode == TestMode::paper_eigen_wave) ?
               EigenRealSpacePerturbation(amp, eig_by_re, eig_by_im, c, s) : 0.0);
  });
  par_for("cgl_lf_quant_init_b3", DevExeSpace(), 0, nmb - 1, ks, ke + 1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int q = i - is;
    const RegionSize block_size = size.d_view(m);
    const Real x = CellCenterX(q, indcs.nx1, block_size.x1min, block_size.x1max);
    const Real s = sin(k_wave*(x - xmin));
    const Real c = cos(k_wave*(x - xmin));
    b0.x3f(m,k,j,i) = bz0 + ((mode == TestMode::paper_eigen_wave) ?
                             EigenRealSpacePerturbation(amp, eig_bz_re, eig_bz_im,
                                                        c, s) : 0.0);
  });

  pmhd->peos->PrimToCons(w0, bcc0, pmhd->u0, is, ie, js, je, ks, ke);
  if (mode == TestMode::hotspot) {
    MonitorHotSpot(pmy_mesh_, 0.0);
    hotspot_initial_parallel = hotspot_min_parallel;
    hotspot_initial_perp = hotspot_min_perp;
  }
}
