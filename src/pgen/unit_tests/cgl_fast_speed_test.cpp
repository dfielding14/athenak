//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_fast_speed_test.cpp
//! \brief Unit checks for the CGL fast speed, active HLLE flux, and timestep route.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/rsolvers/hlle_cgl.hpp"
#include "parameter_input.hpp"
#include "pgen/pgen.hpp"

namespace {

constexpr Real kTol = 2.0e-13;
constexpr Real kBx = 0.80622577482985496524;  // sqrt(0.65)
constexpr Real kBy = 0.59160797830996160426;  // sqrt(0.35)
constexpr Real kBz = 0.0;

struct ObliqueState {
  Real density;
  Real vx;
  Real vy;
  Real vz;
  Real ppar;
  Real pperp;
};

constexpr ObliqueState kLeft{1.0, 0.20, 0.03, -0.02, 1.0, 0.5};
constexpr ObliqueState kRight{0.8, -0.10, -0.04, 0.05, 1.0, 0.5};

void Fail(const std::string &label, const Real got, const Real expected) {
  std::cout << "CGL fast-speed test failed for " << label
            << ": got=" << got << ", expected=" << expected << std::endl;
  std::exit(EXIT_FAILURE);
}

void RequireClose(const std::string &label, const Real got, const Real expected) {
  const Real scale = std::max(static_cast<Real>(1.0), std::abs(expected));
  if (!std::isfinite(got) || std::abs(got - expected) > kTol*scale) {
    Fail(label, got, expected);
  }
}

void Require(const std::string &label, const bool condition) {
  if (!condition) {
    std::cout << "CGL fast-speed test failed for " << label << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

Real LiteratureDiscriminant(const Real ppar, const Real pperp, const Real bx,
                            const Real by, const Real bz) {
  const Real bx2 = bx*bx;
  const Real b2 = bx2 + by*by + bz*bz;
  const Real mu2 = bx2/b2;
  const Real qsq = b2 + 2.0*pperp + (2.0*ppar - pperp)*mu2;
  return qsq*qsq + 4.0*pperp*pperp*(1.0 - mu2)*mu2
       - 12.0*ppar*pperp*mu2*(2.0 - mu2)
       + 12.0*ppar*ppar*mu2*mu2 - 12.0*bx2*ppar;
}

Real LiteratureFastSpeed(const Real density, const Real ppar, const Real pperp,
                         const Real bx, const Real by, const Real bz) {
  const Real b2 = bx*bx + by*by + bz*bz;
  const Real mu2 = bx*bx/b2;
  const Real qsq = b2 + 2.0*pperp + (2.0*ppar - pperp)*mu2;
  return std::sqrt(0.5*(qsq + std::sqrt(LiteratureDiscriminant(
      ppar, pperp, bx, by, bz)))/density);
}

Real LegacyDiscriminant(const Real ppar, const Real pperp, const Real bx,
                        const Real by, const Real bz) {
  const Real bx2 = bx*bx;
  const Real b2 = bx2 + by*by + bz*bz;
  const Real mu2 = bx2/b2;
  const Real qsq = b2 + 2.0*pperp + (2.0*ppar - pperp)*mu2;
  return qsq*qsq + 4.0*pperp*pperp*(1.0 - mu2)*mu2
       - 12.0*ppar*pperp*mu2*(2.0 - mu2)
       + 12.0*ppar*pperp*mu2*mu2 - 12.0*bx2*ppar;
}

Real LegacyRegularizedFastSpeed(const Real density, const Real ppar, const Real pperp,
                                const Real bx, const Real by, const Real bz) {
  const Real b2 = bx*bx + by*by + bz*bz;
  const Real mu2 = bx*bx/b2;
  const Real qsq = b2 + 2.0*pperp + (2.0*ppar - pperp)*mu2;
  return std::sqrt(0.5*(qsq + std::sqrt(std::abs(LegacyDiscriminant(
      ppar, pperp, bx, by, bz))))/density);
}

struct CglStateAndFlux {
  MHDCons1D conserved;
  MHDCons1D flux;
  Real fast_speed;
};

CglStateAndFlux IndependentStateAndFlux(const MHDPrim1D &w, const Real bxi,
                                         const bool literature_correct) {
  CglStateAndFlux result{};
  const Real bsq = bxi*bxi + w.by*w.by + w.bz*w.bz;
  const Real bmag = std::sqrt(bsq);
  const Real pmag = 0.5*bsq;
  const Real firehose_factor = 1.0 + (w.pp - w.e)/bsq;
  const Real vsq = w.vx*w.vx + w.vy*w.vy + w.vz*w.vz;

  result.conserved.d = w.d;
  result.conserved.mx = w.d*w.vx;
  result.conserved.my = w.d*w.vy;
  result.conserved.mz = w.d*w.vz;
  result.conserved.e = 0.5*w.e + w.pp + 0.5*w.d*vsq + pmag;
  result.conserved.mu = w.d*std::log(w.pp/w.e*w.d*w.d/(bmag*bsq));
  result.conserved.by = w.by;
  result.conserved.bz = w.bz;

  result.flux.d = w.d*w.vx;
  result.flux.mx = w.d*w.vx*w.vx + pmag + w.pp - bxi*bxi*firehose_factor;
  result.flux.my = w.d*w.vy*w.vx - bxi*w.by*firehose_factor;
  result.flux.mz = w.d*w.vz*w.vx - bxi*w.bz*firehose_factor;
  result.flux.e = result.conserved.e*w.vx + w.vx*(w.pp + pmag)
                - bxi*(bxi*w.vx + w.by*w.vy + w.bz*w.vz)*firehose_factor;
  result.flux.mu = result.conserved.mu*w.vx;
  result.flux.by = w.by*w.vx - bxi*w.vy;
  result.flux.bz = w.bz*w.vx - bxi*w.vz;

  if (literature_correct) {
    result.fast_speed = LiteratureFastSpeed(w.d, w.e, w.pp, bxi, w.by, w.bz);
  } else {
    result.fast_speed = LegacyRegularizedFastSpeed(w.d, w.e, w.pp, bxi, w.by, w.bz);
  }
  return result;
}

void SubtractScaled(MHDCons1D &flux, const Real speed, const MHDCons1D &conserved) {
  flux.d -= speed*conserved.d;
  flux.mx -= speed*conserved.mx;
  flux.my -= speed*conserved.my;
  flux.mz -= speed*conserved.mz;
  flux.e -= speed*conserved.e;
  flux.mu -= speed*conserved.mu;
  flux.by -= speed*conserved.by;
  flux.bz -= speed*conserved.bz;
}

Real HlleCombine(const Real left, const Real right, const Real tmp) {
  return 0.5*(left + right) + (left - right)*tmp;
}

MHDCons1D IndependentHlleFlux(const MHDPrim1D &wl, const MHDPrim1D &wr,
                              const Real bxi, const bool literature_correct) {
  CglStateAndFlux left = IndependentStateAndFlux(wl, bxi, literature_correct);
  CglStateAndFlux right = IndependentStateAndFlux(wr, bxi, literature_correct);
  const Real al = std::min(wr.vx - right.fast_speed, wl.vx - left.fast_speed);
  const Real ar = std::max(wr.vx + right.fast_speed, wl.vx + left.fast_speed);
  const Real bp = (ar > 0.0) ? ar : 1.0e-20;
  const Real bm = (al < 0.0) ? al : -1.0e-20;
  const Real tmp = 0.5*(bp + bm)/(bp - bm);

  SubtractScaled(left.flux, bm, left.conserved);
  SubtractScaled(right.flux, bp, right.conserved);

  MHDCons1D result{};
  result.d = HlleCombine(left.flux.d, right.flux.d, tmp);
  result.mx = HlleCombine(left.flux.mx, right.flux.mx, tmp);
  result.my = HlleCombine(left.flux.my, right.flux.my, tmp);
  result.mz = HlleCombine(left.flux.mz, right.flux.mz, tmp);
  result.e = HlleCombine(left.flux.e, right.flux.e, tmp);
  result.mu = (result.d >= 0.0)
      ? result.d*left.conserved.mu/left.conserved.d
      : result.d*right.conserved.mu/right.conserved.d;
  result.by = -HlleCombine(left.flux.by, right.flux.by, tmp);
  result.bz = HlleCombine(left.flux.bz, right.flux.bz, tmp);
  return result;
}

template <typename FluxView, typename E3View, typename E2View>
void CheckXFlux(const std::string &label, const FluxView &flux, const E3View &e3x1,
                const E2View &e2x1, const MHDCons1D &expected,
                const int m, const int k, const int j, const int i) {
  RequireClose(label + ".density", flux(m,IDN,k,j,i), expected.d);
  RequireClose(label + ".momentum-x", flux(m,IM1,k,j,i), expected.mx);
  RequireClose(label + ".momentum-y", flux(m,IM2,k,j,i), expected.my);
  RequireClose(label + ".momentum-z", flux(m,IM3,k,j,i), expected.mz);
  RequireClose(label + ".energy", flux(m,IEN,k,j,i), expected.e);
  RequireClose(label + ".anisotropy", flux(m,IMU,k,j,i), expected.mu);
  RequireClose(label + ".electric-z", e3x1(m,k,j,i), expected.by);
  RequireClose(label + ".electric-y", e2x1(m,k,j,i), expected.bz);
}

void CheckObliqueRegression(const EOS_Data &eos) {
  // This state is inside the CGL hyperbolic bounds, but the legacy
  // p_parallel*p_perp coefficient makes its fast-mode discriminant negative.
  const Real density = 1.0;
  const Real ppar = 1.0;
  const Real pperp = 0.5;
  const Real bx = std::sqrt(0.65);
  const Real by = std::sqrt(0.35);
  const Real bz = 0.0;
  const Real corrected = LiteratureDiscriminant(ppar, pperp, bx, by, bz);
  const Real legacy = LegacyDiscriminant(ppar, pperp, bx, by, bz);
  Require("oblique corrected discriminant is positive", corrected > 1.0);
  Require("oblique legacy discriminant is negative", legacy < -1.0);
  RequireClose("oblique corrected discriminant", corrected, 1.083125);
  RequireClose(
      "oblique fast speed",
      eos.IdealMHDFastSpeed(density, ppar, pperp, bx, by, bz, eos.bfloor),
      1.4169920456419176);
}

void CheckDirectionalLimits(const EOS_Data &eos) {
  const Real perpendicular = eos.IdealMHDFastSpeed(
      2.0, 1.2, 0.7, 0.0, 1.0, 0.0, eos.bfloor);
  RequireClose("perpendicular fast speed", perpendicular, std::sqrt(1.2));

  const Real density = 1.3;
  const Real ppar = 0.9;
  const Real pperp = 0.6;
  const Real bx = 0.8;
  const Real parallel = eos.IdealMHDFastSpeed(
      density, ppar, pperp, bx, 0.0, 0.0, eos.bfloor);
  RequireClose(
      "parallel fast speed", parallel,
      LiteratureFastSpeed(density, ppar, pperp, bx, 0.0, 0.0));
}

MHDPrim1D MakeLocalState(const ObliqueState &state) {
  MHDPrim1D w{};
  w.d = state.density;
  w.vx = state.vx;
  w.vy = state.vy;
  w.vz = state.vz;
  w.e = state.ppar;
  w.pp = state.pperp;
  w.by = kBy;
  w.bz = kBz;
  return w;
}

void CheckActiveReconstructedHlleRoute(const EOS_Data &eos) {
  constexpr int nvars = IAN + 1;
  RegionIndcs indcs{};
  CoordData coord{};
  DualArray1D<RegionSize> size("cgl-fast-speed-size", 1);
  DvceArray4D<Real> bx("cgl-fast-speed-bx", 1, 1, 1, 1);
  DvceArray5D<Real> flux("cgl-fast-speed-flux", 1, nvars, 1, 1, 1);
  DvceArray4D<Real> e3x1("cgl-fast-speed-e3x1", 1, 1, 1, 1);
  DvceArray4D<Real> e2x1("cgl-fast-speed-e2x1", 1, 1, 1, 1);
  DvceArray5D<Real> pflux("cgl-fast-speed-pflux", 1, 6, 1, 1, 1);
  Kokkos::deep_copy(bx, kBx);

  const size_t scratch_size =
      2*ScrArray2D<Real>::shmem_size(nvars, 1)
    + 2*ScrArray2D<Real>::shmem_size(NMAG, 1);
  const ObliqueState left = kLeft;
  const ObliqueState right = kRight;
  par_for_outer("cgl_fast_speed_active_hlle", DevExeSpace(), scratch_size, 0, 0, 0,
  KOKKOS_LAMBDA(TeamMember_t member, const int) {
    ScrArray2D<Real> wl(member.team_scratch(0), nvars, 1);
    ScrArray2D<Real> wr(member.team_scratch(0), nvars, 1);
    ScrArray2D<Real> bl(member.team_scratch(0), NMAG, 1);
    ScrArray2D<Real> br(member.team_scratch(0), NMAG, 1);
    par_for_inner(member, 0, 0, [&](const int i) {
      wl(IDN,i) = left.density;
      wl(IVX,i) = left.vx;
      wl(IVY,i) = left.vy;
      wl(IVZ,i) = left.vz;
      wl(IPR,i) = left.ppar;
      wl(IPP,i) = left.pperp;
      wr(IDN,i) = right.density;
      wr(IVX,i) = right.vx;
      wr(IVY,i) = right.vy;
      wr(IVZ,i) = right.vz;
      wr(IPR,i) = right.ppar;
      wr(IPP,i) = right.pperp;
      bl(IBX,i) = kBx;
      bl(IBY,i) = kBy;
      bl(IBZ,i) = kBz;
      br(IBX,i) = kBx;
      br(IBY,i) = kBy;
      br(IBZ,i) = kBz;
    });
    member.team_barrier();
    mhd::HLLE_CGL(member, eos, indcs, size, coord, 0, 0, 0, 0, 0, IVX,
                  wl, wr, bl, br, bx, flux, e3x1, e2x1, false, pflux);
  });
  Kokkos::fence();

  auto hflux = Kokkos::create_mirror_view_and_copy(HostMemSpace(), flux);
  auto he3x1 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), e3x1);
  auto he2x1 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), e2x1);
  const MHDPrim1D wl = MakeLocalState(kLeft);
  const MHDPrim1D wr = MakeLocalState(kRight);
  const MHDCons1D corrected = IndependentHlleFlux(wl, wr, kBx, true);
  const MHDCons1D legacy = IndependentHlleFlux(wl, wr, kBx, false);
  CheckXFlux("active reconstructed HLLE literature route", hflux, he3x1, he2x1,
             corrected, 0, 0, 0, 0);
  Require("active reconstructed HLLE differs from legacy density flux",
          std::abs(corrected.d - legacy.d) > 1.0e-3);
  Require("active reconstructed HLLE differs from legacy momentum flux",
          std::abs(corrected.mx - legacy.mx) > 1.0e-3);
  std::cout << "CGL active reconstructed-HLLE route checks passed" << std::endl;
}

} // namespace

void RunCglFastSpeedChecks() {
  EOS_Data eos{};
  eos.gamma = 5.0/3.0;
  eos.bfloor = 1.0e-10;
  CheckObliqueRegression(eos);
  CheckDirectionalLimits(eos);
  CheckActiveReconstructedHlleRoute(eos);
  std::cout << "CGL fast-speed checks passed" << std::endl;
}

void RunCglFastSpeedStandaloneChecks() {
  Kokkos::initialize();
  {
    RunCglFastSpeedChecks();
  }
  Kokkos::finalize();
}

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  (void) pin;
  if (restart) return;
  RunCglFastSpeedChecks();
}
