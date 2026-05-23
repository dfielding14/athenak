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
#include <sstream>
#include <string>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "eos/cgl_physics.hpp"
#include "mhd/mhd.hpp"
#include "parameter_input.hpp"
#include "pgen/pgen.hpp"

namespace {

enum class TestMode {
  parallel_decay,
  perp_decay,
  grad_b,
  flux_limiter,
  limiter_heat_flux_suppression,
  limiter_stress,
  field_aligned_wave,
  paper_oblique_wave,
  paper_eigen_wave
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
  if (mode == "grad_b") return TestMode::grad_b;
  if (mode == "flux_limiter") return TestMode::flux_limiter;
  if (mode == "limiter_heat_flux_suppression") {
    return TestMode::limiter_heat_flux_suppression;
  }
  if (mode == "limiter_stress") return TestMode::limiter_stress;
  if (mode == "field_aligned_wave") return TestMode::field_aligned_wave;
  if (mode == "paper_oblique_wave") return TestMode::paper_oblique_wave;
  if (mode == "paper_eigen_wave") return TestMode::paper_eigen_wave;
  Fail("<problem>/test_mode must be parallel_decay, perp_decay, grad_b, "
       "flux_limiter, limiter_heat_flux_suppression, limiter_stress, "
       "field_aligned_wave, paper_oblique_wave, or paper_eigen_wave");
}

const char *ModeName(const TestMode mode) {
  switch (mode) {
    case TestMode::parallel_decay: return "parallel_decay";
    case TestMode::perp_decay: return "perp_decay";
    case TestMode::grad_b: return "grad_b";
    case TestMode::flux_limiter: return "flux_limiter";
    case TestMode::limiter_heat_flux_suppression: return "limiter_heat_flux_suppression";
    case TestMode::limiter_stress: return "limiter_stress";
    case TestMode::field_aligned_wave: return "field_aligned_wave";
    case TestMode::paper_oblique_wave: return "paper_oblique_wave";
    case TestMode::paper_eigen_wave: return "paper_eigen_wave";
  }
  return "unknown";
}

void RequireOneDimensionalSingleBlock(Mesh *pm) {
  Require(pm->pmb_pack->nmb_thispack == 1,
          "this unit test currently expects one MeshBlock on one rank");
  Require(pm->mb_indcs.nx2 == 1 && pm->mb_indcs.nx3 == 1,
          "this unit test currently expects a 1D mesh");
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
Projection ProjectPrimitive(const HostView &w, ParameterInput *pin, Mesh *pm,
                            const int idx) {
  RequireOneDimensionalSingleBlock(pm);
  const int is = pm->mb_indcs.is;
  const int js = pm->mb_indcs.js;
  const int ks = pm->mb_indcs.ks;
  const int nx1 = pm->mb_indcs.nx1;
  const Real k_wave = Wavenumber(pin, pm);
  const Real xmin = pm->mesh_size.x1min;

  Projection p;
  for (int q = 0; q < nx1; ++q) {
    p.mean += w(0,idx,ks,js,is + q);
  }
  p.mean /= static_cast<Real>(nx1);

  for (int q = 0; q < nx1; ++q) {
    const Real x = XCenter(pm, q);
    const Real phase = k_wave*(x - xmin);
    const Real value = w(0,idx,ks,js,is + q) - p.mean;
    p.sin_amp += value*std::sin(phase);
    p.cos_amp += value*std::cos(phase);
  }
  p.sin_amp *= 2.0/static_cast<Real>(nx1);
  p.cos_amp *= 2.0/static_cast<Real>(nx1);
  return p;
}

template <typename HostView>
Projection ProjectCellField(const HostView &bcc, ParameterInput *pin, Mesh *pm,
                            const int idx) {
  RequireOneDimensionalSingleBlock(pm);
  const int is = pm->mb_indcs.is;
  const int js = pm->mb_indcs.js;
  const int ks = pm->mb_indcs.ks;
  const int nx1 = pm->mb_indcs.nx1;
  const Real k_wave = Wavenumber(pin, pm);
  const Real xmin = pm->mesh_size.x1min;

  Projection p;
  for (int q = 0; q < nx1; ++q) {
    p.mean += bcc(0,idx,ks,js,is + q);
  }
  p.mean /= static_cast<Real>(nx1);

  for (int q = 0; q < nx1; ++q) {
    const Real x = XCenter(pm, q);
    const Real phase = k_wave*(x - xmin);
    const Real value = bcc(0,idx,ks,js,is + q) - p.mean;
    p.sin_amp += value*std::sin(phase);
    p.cos_amp += value*std::cos(phase);
  }
  p.sin_amp *= 2.0/static_cast<Real>(nx1);
  p.cos_amp *= 2.0/static_cast<Real>(nx1);
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

Real LimiterCollisionRate(ParameterInput *pin, const Real ppar, const Real pperp,
                          const Real bx, const Real by, const Real bz) {
  const bool mlim = GetBooleanOrFalse(pin, "mhd", "mirror_limiter");
  const bool flim = GetBooleanOrFalse(pin, "mhd", "firehose_limiter");
  const bool backup_lim = GetBooleanOrFalse(pin, "mhd", "backup_limiters");
  const Real lim_coll = std::max(GetRealOrZero(pin, "mhd", "limiter_nu_coll"),
                                 static_cast<Real>(0.0));
  const Real bsqr = SQR(bx) + SQR(by) + SQR(bz);
  return cgl::LimiterCollisionRate(ppar, pperp, bsqr, lim_coll, mlim, flim,
                                   backup_lim);
}

Real EffectiveCollisionFrequency(ParameterInput *pin, const Real ppar, const Real pperp,
                                 const Real bx, const Real by, const Real bz) {
  return std::max(BackgroundCollisionFrequency(pin), static_cast<Real>(0.0)) +
         LimiterCollisionRate(pin, ppar, pperp, bx, by, bz);
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
  const Real denom = cgl::kSqrtTwoPi*cpar*lf_k + nu_eff;
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
  const Real bz = pin->GetOrAddReal("problem", "bz0", 0.0);
  const Real by = 0.5*(ByCell(pin, pm, left) + ByCell(pin, pm, right));
  const Real bmag_face = std::sqrt(SQR(bx) + SQR(by) + SQR(bz));
  const Real bhx = bx/bmag_face;
  const Real grad_b_x = (BMagCell(pin, pm, right) - BMagCell(pin, pm, left))/dx;
  const Real gradpar_b = bhx*grad_b_x;
  const Real rho = pin->GetOrAddReal("problem", "rho0", 1.0);
  const Real ppar = pin->GetOrAddReal("problem", "ppar0", 1.0);
  const Real pperp = pin->GetOrAddReal("problem", "pperp0", 1.2);
  const Real lf_k = pin->GetReal("mhd", "lf_k_parallel");
  const Real cpar = FaceCParallel(pin, rho, ppar);
  const Real nu_eff = EffectiveCollisionFrequency(pin, ppar, pperp, bx, by, bz);
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

  const auto apply_heat_flux = [=](DecayState state, const Real dt) {
    state.tpar *= std::exp(-chi_parallel*SQR(k_wave)*dt);
    state.tperp *= std::exp(-chi_perp*SQR(k_wave)*dt);
    return state;
  };
  const auto apply_collision = [=](DecayState state, const Real dt) {
    const Real piso = ONE_3RD*state.tpar + TWO_3RDS*state.tperp;
    Real paniso = state.tperp - state.tpar;
    paniso *= std::exp(-nu_coll*dt);
    state.tpar = piso - TWO_3RDS*paniso;
    state.tperp = piso + ONE_3RD*paniso;
    return state;
  };

  // Driver order with STS enabled: pre half-sweep, CGL collision after the main
  // integrator, post half-sweep, then the post-STS CGL collision.
  s = apply_heat_flux(s, 0.5*pm->time);
  s = apply_collision(s, pm->time);
  s = apply_heat_flux(s, 0.5*pm->time);
  s = apply_collision(s, pm->time);
  return s;
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
  const Real nu_eff = EffectiveCollisionFrequency(pin, ppar0, pperp0, bx0, by0, bz0);
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
      EffectiveCollisionFrequency(pin, ppar_face, pperp_face, bx, by, bz);
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
      EffectiveCollisionFrequency(pin, ppar_face, pperp_face, bx, by, bz);
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

void CheckLimiterHeatFluxSuppression(ParameterInput *pin, Mesh *pm) {
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
  const Real backup_tol = pin->GetOrAddReal("problem", "backup_bound_tol", 1.02);
  const std::string limiter_kind =
      pin->GetOrAddString("problem", "limiter_kind", "mirror");

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
    if (limiter_kind == "mirror") {
      Require(paniso <= backup_tol*bsqr, "mirror stress exceeded the backup bound");
    } else if (limiter_kind == "firehose") {
      Require(paniso >= -backup_tol*bsqr, "firehose stress exceeded the backup bound");
    } else {
      Fail("<problem>/limiter_kind must be mirror or firehose");
    }
  }
  std::cout << "CGL LF limiter_stress passed for " << limiter_kind << std::endl;
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
  const Real nu_eff = EffectiveCollisionFrequency(pin, ppar0, pperp0, bx0, by0, bz0);
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
  const Real wave_tol = pin->GetOrAddReal("problem", "wave_rel_tol", 2.5e-1);
  const Real c_cgl = std::sqrt(3.0*ppar0/rho0);
  const Complex rho_m = ProjectionToComplex(rho_p);
  const Complex vx_m = ProjectionToComplex(vx_p);
  const Complex ppar_m = ProjectionToComplex(ppar_p);

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
    const Real nu_eff = EffectiveCollisionFrequency(pin, p0, p0, bx0, by0, bz0);
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
  const Real wave_tol = pin->GetOrAddReal("problem", "wave_rel_tol", 4.0e-1);
  const Complex vy_m = ProjectionToComplex(ProjectPrimitive(w, pin, pm, IVY));
  const Complex by_m = ProjectionToComplex(ProjectCellField(bcc, pin, pm, IBY));
  const Complex ppar_m = ProjectionToComplex(ProjectPrimitive(w, pin, pm, IPR));
  const Complex pperp_m = ProjectionToComplex(ProjectPrimitive(w, pin, pm, IPP));

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
  const Real wave_tol = pin->GetOrAddReal("problem", "eigen_wave_rel_tol", 7.5e-2);
  const Real component_floor =
      pin->GetOrAddReal("problem", "eigen_component_floor", 1.0e-4);
  const Real zero_abs_tol =
      pin->GetOrAddReal("problem", "eigen_zero_abs_tol", 1.0e-8);
  const Complex lambda(pin->GetReal("problem", "eigen_lambda_re"),
                       pin->GetReal("problem", "eigen_lambda_im"));
  const Complex phase = std::exp(lambda*pm->time);
  const PaperWaveState ref = (amp*phase)*eigen;

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

void FinalizeCGLLFQuantitative(ParameterInput *pin, Mesh *pm) {
  auto *pmhd = pm->pmb_pack->pmhd;
  Require(pmhd != nullptr && pmhd->peos->eos_data.is_cgl,
          "quantitative LF tests require <mhd>/eos = cgl");
  RequireOneDimensionalSingleBlock(pm);
  const TestMode mode = ParseMode(pin);
  if (mode == TestMode::parallel_decay || mode == TestMode::perp_decay) {
    CheckDecay(pin, pm, mode);
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
  RequireOneDimensionalSingleBlock(pmy_mesh_);

  const TestMode mode = ParseMode(pin);
  const Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  const Real ppar0 = pin->GetOrAddReal("problem", "ppar0", 1.0);
  const Real pperp0 = pin->GetOrAddReal("problem", "pperp0", 1.0);
  const Real amp = pin->GetOrAddReal("problem", "amp", 1.0e-4);
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
  const Real xmin = pmy_mesh_->mesh_size.x1min;
  const Real xmax = pmy_mesh_->mesh_size.x1max;
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

  par_for("cgl_lf_quant_init_prim", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int q = i - is;
    const Real x = CellCenterX(q, indcs.nx1, xmin, xmax);
    const Real s = sin(k_wave*(x - xmin));
    const Real c = cos(k_wave*(x - xmin));
    Real rho = rho0;
    Real vx = 0.0;
    Real vy = 0.0;
    Real vz = 0.0;
    Real ppar = ppar0;
    Real pperp = pperp0;
    Real by = (mode == TestMode::grad_b) ? by_amp*s : by0;
    Real bz = bz0;

    if (mode == TestMode::parallel_decay) {
      ppar = ppar0*(1.0 + amp*s);
    } else if (mode == TestMode::perp_decay) {
      pperp = pperp0*(1.0 + amp*s);
    } else if (mode == TestMode::flux_limiter ||
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
    const Real x = CellCenterX(q, indcs.nx1, xmin, xmax);
    const Real s = sin(k_wave*(x - xmin));
    const Real c = cos(k_wave*(x - xmin));
    b0.x2f(m,k,j,i) = (mode == TestMode::grad_b) ? by_amp*s :
        by0 + ((mode == TestMode::paper_eigen_wave) ?
               EigenRealSpacePerturbation(amp, eig_by_re, eig_by_im, c, s) : 0.0);
  });
  par_for("cgl_lf_quant_init_b3", DevExeSpace(), 0, nmb - 1, ks, ke + 1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int q = i - is;
    const Real x = CellCenterX(q, indcs.nx1, xmin, xmax);
    const Real s = sin(k_wave*(x - xmin));
    const Real c = cos(k_wave*(x - xmin));
    b0.x3f(m,k,j,i) = bz0 + ((mode == TestMode::paper_eigen_wave) ?
                             EigenRealSpacePerturbation(amp, eig_bz_re, eig_bz_im,
                                                        c, s) : 0.0);
  });

  pmhd->peos->PrimToCons(w0, bcc0, pmhd->u0, is, ie, js, je, ks, ke);
}
