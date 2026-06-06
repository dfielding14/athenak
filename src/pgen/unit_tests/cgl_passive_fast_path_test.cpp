//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_passive_fast_path_test.cpp
//! \brief Focused checks that passive CGL flow uses only isothermal-MHD wave speeds.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "eos/cgl_physics.hpp"
#include "mhd/rsolvers/llf_mhd_singlestate.hpp"
#include "parameter_input.hpp"
#include "pgen/pgen.hpp"

namespace {

constexpr Real kTol = 2.0e-13;
constexpr Real kBx = 0.80622577482985496524;  // sqrt(0.65)
constexpr Real kBy = 0.59160797830996160426;  // sqrt(0.35)

void Fail(const std::string &label, const Real got, const Real expected) {
  std::cout << "Passive CGL fast-path test failed for " << label
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
    std::cout << "Passive CGL fast-path test failed for " << label << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

Real IsothermalFastSpeed(const EOS_Data &eos, const Real d, const Real bx,
                         const Real by, const Real bz) {
  const Real asq = eos.iso_cs*eos.iso_cs*d;
  const Real bperp2 = by*by + bz*bz;
  const Real qsq = bx*bx + bperp2 + asq;
  const Real tmp = bx*bx + bperp2 - asq;
  return std::sqrt(0.5*(qsq + std::sqrt(tmp*tmp + 4.0*asq*bperp2))/d);
}

MHDPrim1D MakeState(const Real d, const Real vx, const Real vy, const Real vz,
                    const Real ppar, const Real pperp, const Real by, const Real bz) {
  MHDPrim1D state{};
  state.d = d;
  state.vx = vx;
  state.vy = vy;
  state.vz = vz;
  state.e = ppar;
  state.pp = pperp;
  state.by = by;
  state.bz = bz;
  return state;
}

void CheckFlowFlux(const std::string &label, const MHDCons1D &got,
                   const MHDCons1D &expected) {
  RequireClose(label + ".density", got.d, expected.d);
  RequireClose(label + ".momentum-x", got.mx, expected.mx);
  RequireClose(label + ".momentum-y", got.my, expected.my);
  RequireClose(label + ".momentum-z", got.mz, expected.mz);
  RequireClose(label + ".field-y", got.by, expected.by);
  RequireClose(label + ".field-z", got.bz, expected.bz);
}

void CheckAdmissible(const std::string &label, const MHDPrim1D &state) {
  const Real bsqr = kBx*kBx + state.by*state.by + state.bz*state.bz;
  const Real paniso = state.pp - state.e;
  Require(label + " above firehose hard bound",
          !cgl::FirehoseHardBoundViolated(paniso, bsqr));
  Require(label + " below mirror hard bound",
          !cgl::MirrorHardBoundViolated(paniso, bsqr));
}

void CheckPassiveSignalAndFluxPaths(const EOS_Data &passive) {
  Require("CGL EOS enabled", passive.is_cgl);
  Require("passive CGL enabled", passive.passive);
  RequireClose("configured isothermal sound speed", passive.iso_cs, 0.7);

  const Real expected_iso = IsothermalFastSpeed(passive, 1.0, kBx, kBy, 0.0);
  const Real passive_speed = passive.IdealMHDFastSpeed(1.0, kBx, kBy, 0.0);
  const Real active_speed =
      passive.IdealMHDFastSpeed(1.0, 1.0, 0.5, kBx, kBy, 0.0, passive.bfloor);
  RequireClose("passive isothermal overload", passive_speed, expected_iso);
  RequireClose("corrected active overload control", active_speed, 1.4169920456419176);
  Require("active and passive overloads are distinguishable",
          std::abs(active_speed - passive_speed) > 0.2);

  const MHDPrim1D left_a =
      MakeState(1.21, 0.23, -0.11, 0.07, 1.0, 0.5, kBy, 0.08);
  const MHDPrim1D right_a =
      MakeState(0.83, -0.17, 0.09, -0.04, 0.7, 0.9, kBy - 0.12, -0.05);
  const MHDPrim1D left_b =
      MakeState(1.21, 0.23, -0.11, 0.07, 1.5, 1.0, kBy, 0.08);
  const MHDPrim1D right_b =
      MakeState(0.83, -0.17, 0.09, -0.04, 1.2, 1.0, kBy - 0.12, -0.05);
  CheckAdmissible("left state A", left_a);
  CheckAdmissible("right state A", right_a);
  CheckAdmissible("left state B", left_b);
  CheckAdmissible("right state B", right_b);

  MHDCons1D ua{}, ub{}, fa{}, fb{};
  Real cfa = 0.0;
  Real cfb = 0.0;
  mhd::SingleStateLLF_CGLStateAndFlux(left_a, kBx, passive, ua, fa, cfa);
  mhd::SingleStateLLF_CGLStateAndFlux(left_b, kBx, passive, ub, fb, cfb);
  const Real expected_left =
      IsothermalFastSpeed(passive, left_a.d, kBx, left_a.by, left_a.bz);
  RequireClose("passive state signal speed A", cfa, expected_left);
  RequireClose("passive state signal speed pressure invariance", cfb, cfa);
  CheckFlowFlux("passive physical flow flux pressure invariance", fb, fa);

  MHDCons1D passive_flux_a{}, passive_flux_b{};
  mhd::SingleStateLLF_CGL(left_a, right_a, kBx, passive, passive_flux_a);
  mhd::SingleStateLLF_CGL(left_b, right_b, kBx, passive, passive_flux_b);
  CheckFlowFlux("passive LLF flow flux pressure invariance", passive_flux_b,
                passive_flux_a);

  EOS_Data active = passive;
  active.passive = false;
  MHDCons1D active_flux_a{}, active_flux_b{};
  mhd::SingleStateLLF_CGL(left_a, right_a, kBx, active, active_flux_a);
  mhd::SingleStateLLF_CGL(left_b, right_b, kBx, active, active_flux_b);
  Require("active LLF control responds to CGL pressures",
          std::abs(active_flux_b.mx - active_flux_a.mx) > 0.1);
}

} // namespace

void RunCglPassiveFastPathChecks() {
  EOS_Data passive{};
  passive.gamma = 5.0/3.0;
  passive.iso_cs = 0.7;
  passive.is_ideal = true;
  passive.is_cgl = true;
  passive.passive = true;
  passive.bfloor = 1.0e-10;
  CheckPassiveSignalAndFluxPaths(passive);
  std::cout << "Passive CGL signal-speed and flux checks passed" << std::endl;
}
