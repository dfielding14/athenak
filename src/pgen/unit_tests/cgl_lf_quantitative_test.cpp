//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_lf_quantitative_test.cpp
//! \brief Quantitative unit-style tests for CGL Landau-fluid heat-flux STS.

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "parameter_input.hpp"
#include "pgen/pgen.hpp"

namespace {

constexpr Real kSqrtTwoOverPi = 0.7978845608028654;
constexpr Real kSqrtEightOverPi = 1.5957691216057308;

enum class TestMode {
  parallel_decay,
  perp_decay,
  grad_b,
  flux_limiter,
  limiter_stress,
  field_aligned_wave,
  paper_oblique_wave
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
  if (mode == "limiter_stress") return TestMode::limiter_stress;
  if (mode == "field_aligned_wave") return TestMode::field_aligned_wave;
  if (mode == "paper_oblique_wave") return TestMode::paper_oblique_wave;
  Fail("<problem>/test_mode must be parallel_decay, perp_decay, grad_b, "
       "flux_limiter, limiter_stress, field_aligned_wave, or paper_oblique_wave");
}

const char *ModeName(const TestMode mode) {
  switch (mode) {
    case TestMode::parallel_decay: return "parallel_decay";
    case TestMode::perp_decay: return "perp_decay";
    case TestMode::grad_b: return "grad_b";
    case TestMode::flux_limiter: return "flux_limiter";
    case TestMode::limiter_stress: return "limiter_stress";
    case TestMode::field_aligned_wave: return "field_aligned_wave";
    case TestMode::paper_oblique_wave: return "paper_oblique_wave";
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
  const Real by = ByCell(pin, pm, q);
  return std::sqrt(SQR(bx0) + SQR(by));
}

Real LimitedHeatFlux(const Real q_unlimited, const Real q_max) {
  if (q_max <= 0.0) {
    return 0.0;
  }
  return q_unlimited*q_max/(q_max + std::abs(q_unlimited));
}

Real GradBMomentFlux(ParameterInput *pin, Mesh *pm, const int face) {
  const int nx1 = pm->mb_indcs.nx1;
  const Real dx = (pm->mesh_size.x1max - pm->mesh_size.x1min)/static_cast<Real>(nx1);
  const int left = (face - 1 + nx1)%nx1;
  const int right = face%nx1;
  const Real bx = pin->GetOrAddReal("problem", "b0", 1.0);
  const Real by = 0.5*(ByCell(pin, pm, left) + ByCell(pin, pm, right));
  const Real bmag_face = std::sqrt(SQR(bx) + SQR(by));
  const Real bhx = bx/bmag_face;
  const Real grad_b_x = (BMagCell(pin, pm, right) - BMagCell(pin, pm, left))/dx;
  const Real gradpar_b = bhx*grad_b_x;
  const Real ppar = pin->GetOrAddReal("problem", "ppar0", 1.0);
  const Real pperp = pin->GetOrAddReal("problem", "pperp0", 1.2);
  const Real cpar0 = pin->GetOrAddReal("mhd", "lf_c_parallel0", 1.0);
  const Real lf_k = pin->GetReal("mhd", "lf_k_parallel");
  const Real chi_perp = kSqrtTwoOverPi*cpar0/lf_k;
  const Real qperp_l = -chi_perp*(-pperp*(1.0 - pperp/ppar)*gradpar_b/bmag_face);
  const Real qperp = LimitedHeatFlux(qperp_l, kSqrtTwoOverPi*cpar0*pperp);
  return bhx*qperp/bmag_face;
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
  const Real cpar0 = pin->GetOrAddReal("mhd", "lf_c_parallel0", 1.0);
  const Real lf_k = pin->GetReal("mhd", "lf_k_parallel");
  const Real k_wave = Wavenumber(pin, pm);
  const Real chi = (mode == TestMode::parallel_decay ? kSqrtEightOverPi :
                                                       kSqrtTwoOverPi)*cpar0/lf_k;
  const Real initial_amp = (p0/rho0)*amp;
  const Real expected_amp = initial_amp*std::exp(-chi*SQR(k_wave)*pm->time);

  RequireRelative(std::string(ModeName(mode)) + " Fourier amplitude",
                  std::abs(projection.sin_amp), expected_amp, rel_tol);
  Require(std::abs(projection.cos_amp) <= phase_tol*std::abs(initial_amp),
          std::string(ModeName(mode)) + " developed an unexpected cosine component");
  std::cout << "CGL LF " << ModeName(mode) << " passed: measured_amp="
            << projection.sin_amp << " expected_amp=" << expected_amp << std::endl;
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
  }

  Require(ref2 > 0.0, "grad_b reference response is zero");
  const Real rel_err = std::sqrt(err2/ref2);
  Require(rel_err <= rel_tol, "grad_b response magnitude is outside tolerance");
  Require(dot > 0.0, "grad_b response has the wrong sign");
  std::cout << "CGL LF grad_b passed: rms_rel_err=" << rel_err << std::endl;
}

Real InitialParallelPressure(ParameterInput *pin, Mesh *pm, const int q) {
  const Real ppar0 = pin->GetOrAddReal("problem", "ppar0", 1.0);
  const Real amp = pin->GetOrAddReal("problem", "amp", 0.5);
  const Real k_wave = Wavenumber(pin, pm);
  const Real xmin = pm->mesh_size.x1min;
  return ppar0*(1.0 + amp*std::sin(k_wave*(XCenter(pm, q) - xmin)));
}

Real LimitedParallelHeatFlux(ParameterInput *pin, Mesh *pm, const int face,
                             Real &q_unlimited, Real &q_max) {
  const int nx1 = pm->mb_indcs.nx1;
  const Real dx = (pm->mesh_size.x1max - pm->mesh_size.x1min)/static_cast<Real>(nx1);
  const int left = (face - 1 + nx1)%nx1;
  const int right = face%nx1;
  const Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  const Real cpar0 = pin->GetOrAddReal("mhd", "lf_c_parallel0", 1.0);
  const Real lf_k = pin->GetReal("mhd", "lf_k_parallel");
  const Real ppar_l = InitialParallelPressure(pin, pm, left);
  const Real ppar_r = InitialParallelPressure(pin, pm, right);
  const Real ppar_face = std::max(static_cast<Real>(0.5)*(ppar_l + ppar_r),
                                  static_cast<Real>(1.0e-30));
  const Real grad_tpar = (ppar_r/rho0 - ppar_l/rho0)/dx;
  q_unlimited = -(kSqrtEightOverPi*cpar0/lf_k)*rho0*grad_tpar;
  q_max = kSqrtEightOverPi*cpar0*ppar_face;
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

  Real err2 = 0.0;
  Real ref2 = 0.0;
  Real max_ratio = 0.0;
  for (int face = 0; face <= nx1; ++face) {
    Real q_l = 0.0, q_max = 0.0;
    const Real q = LimitedParallelHeatFlux(pin, pm, face, q_l, q_max);
    max_ratio = std::max(max_ratio, std::abs(q_l)/std::max(q_max,
                         static_cast<Real>(1.0e-30)));
    Require(std::abs(q) <= q_max*(1.0 + 1.0e-12),
            "limited heat flux exceeds q_max");
    if (q_l != 0.0) {
      Require(q*q_l > 0.0, "limited heat flux changed sign");
    }
  }

  for (int q = 0; q < nx1; ++q) {
    const int i = is + q;
    Real q_r_l = 0.0, q_r_max = 0.0;
    Real q_l_l = 0.0, q_l_max = 0.0;
    const Real flux_r = LimitedParallelHeatFlux(pin, pm, q + 1, q_r_l, q_r_max);
    const Real flux_l = LimitedParallelHeatFlux(pin, pm, q, q_l_l, q_l_max);
    const Real expected_ppar = InitialParallelPressure(pin, pm, q)
                             + pm->time*(-(flux_r - flux_l)/dx);
    const Real measured_ppar = w(0,IPR,ks,js,i);
    err2 += SQR(measured_ppar - expected_ppar);
    ref2 += SQR(expected_ppar - InitialParallelPressure(pin, pm, q));
  }

  Require(max_ratio >= min_unlimited_ratio,
          "test did not enter the strongly limited heat-flux regime");
  Require(ref2 > 0.0, "flux limiter reference response is zero");
  const Real rel_err = std::sqrt(err2/ref2);
  if (rel_err > rel_tol) {
    std::cout << "flux_limiter rel_err=" << rel_err << " rel_tol=" << rel_tol
              << std::endl;
    Fail("flux limiter update is outside tolerance");
  }
  std::cout << "CGL LF flux_limiter passed: rms_rel_err=" << rel_err
            << " max_unlimited_over_qmax=" << max_ratio << std::endl;
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
  const Real cpar0 = pin->GetOrAddReal("mhd", "lf_c_parallel0", std::sqrt(ppar0/rho0));
  const Real lf_k = pin->GetReal("mhd", "lf_k_parallel");
  const Real chi_parallel = kSqrtEightOverPi*cpar0/lf_k;
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
  const Real rel_err = std::abs(got - expected)/std::max(scale, static_cast<Real>(1.0e-30));
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

  RequireComplexRelative("field_aligned_wave rho", ProjectionToComplex(rho_p), ref.rho,
                         rho0*amp, wave_tol);
  RequireComplexRelative("field_aligned_wave vx", ProjectionToComplex(vx_p), ref.vx,
                         c_cgl*amp, wave_tol);
  RequireComplexRelative("field_aligned_wave p_parallel", ProjectionToComplex(ppar_p),
                         ref.ppar, 3.0*ppar0*amp, wave_tol);
  std::cout << "CGL LF field_aligned_wave passed against linear Fourier reference"
            << std::endl;
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

PaperWaveState PaperWaveRHS(const PaperWaveState &s, const Real rho0,
                            const Real p0, const Real bx0, const Real by0,
                            const Real bz0, const Real k_wave,
                            const Real chi_parallel, const Real chi_perp) {
  const Complex ik(0.0, k_wave);
  const Real bmag0 = std::sqrt(SQR(bx0) + SQR(by0) + SQR(bz0));
  const Real bhx = bx0/bmag0;
  const Real bhy = by0/bmag0;
  const Real bhz = bz0/bmag0;
  const Complex delta_p = s.pperp - s.ppar;
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
    const Real cpar0 = pin->GetOrAddReal("mhd", "lf_c_parallel0", std::sqrt(p0/rho0));
    const Real lf_k = pin->GetReal("mhd", "lf_k_parallel");
    chi_parallel = kSqrtEightOverPi*cpar0/lf_k;
    chi_perp = kSqrtTwoOverPi*cpar0/lf_k;
  }

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

void CheckPaperObliqueWave(ParameterInput *pin, Mesh *pm) {
  auto *pmhd = pm->pmb_pack->pmhd;
  auto w = HostCopy(pmhd->w0);
  auto bcc = HostCopy(pmhd->bcc0);
  const PaperWaveState ref = IntegratePaperWaveReference(pin, pm);
  const Real amp = pin->GetOrAddReal("problem", "amp", 1.0e-5);
  const Real p0 = pin->GetOrAddReal("problem", "ppar0", 5.0);
  const Real wave_tol = pin->GetOrAddReal("problem", "wave_rel_tol", 4.0e-1);

  RequireComplexRelative("paper_oblique_wave vy", ProjectionToComplex(
                         ProjectPrimitive(w, pin, pm, IVY)), ref.vy, amp, wave_tol);
  RequireComplexRelative("paper_oblique_wave By", ProjectionToComplex(
                         ProjectCellField(bcc, pin, pm, IBY)), ref.by, amp, wave_tol);
  RequireComplexRelative("paper_oblique_wave p_parallel", ProjectionToComplex(
                         ProjectPrimitive(w, pin, pm, IPR)), ref.ppar, p0*amp, wave_tol);
  RequireComplexRelative("paper_oblique_wave p_perp", ProjectionToComplex(
                         ProjectPrimitive(w, pin, pm, IPP)), ref.pperp, p0*amp, wave_tol);
  std::cout << "CGL LF paper_oblique_wave passed against linear oblique reference"
            << std::endl;
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
  } else if (mode == TestMode::limiter_stress) {
    CheckLimiterStress(pin, pm);
  } else if (mode == TestMode::field_aligned_wave) {
    CheckFieldAlignedWave(pin, pm);
  } else if (mode == TestMode::paper_oblique_wave) {
    CheckPaperObliqueWave(pin, pm);
  }
}

} // namespace

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
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
    Real rho = rho0;
    Real vx = 0.0;
    Real ppar = ppar0;
    Real pperp = pperp0;

    if (mode == TestMode::parallel_decay || mode == TestMode::flux_limiter) {
      ppar = ppar0*(1.0 + amp*s);
    } else if (mode == TestMode::perp_decay) {
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
    }

    w0(m,IDN,k,j,i) = rho;
    w0(m,IVX,k,j,i) = vx;
    w0(m,IVY,k,j,i) = (mode == TestMode::paper_oblique_wave) ? amp*s : 0.0;
    w0(m,IVZ,k,j,i) = 0.0;
    w0(m,IPR,k,j,i) = ppar;
    w0(m,IPP,k,j,i) = pperp;
    bcc0(m,IBX,k,j,i) = bx0;
    bcc0(m,IBY,k,j,i) = (mode == TestMode::grad_b) ? by_amp*s : by0;
    bcc0(m,IBZ,k,j,i) = bz0;
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
    b0.x2f(m,k,j,i) = (mode == TestMode::grad_b) ? by_amp*s : by0;
  });
  par_for("cgl_lf_quant_init_b3", DevExeSpace(), 0, nmb - 1, ks, ke + 1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    b0.x3f(m,k,j,i) = bz0;
  });

  pmhd->peos->PrimToCons(w0, bcc0, pmhd->u0, is, ie, js, je, ks, ke);
}
