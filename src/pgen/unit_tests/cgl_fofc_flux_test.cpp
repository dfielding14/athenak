//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_fofc_flux_test.cpp
//! \brief Unit test for CGL FOFC single-state LLF flux signs and permutations.

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/rsolvers/llf_mhd_singlestate.hpp"
#include "parameter_input.hpp"
#include "pgen/pgen.hpp"

namespace {

constexpr Real kTol = 1.0e-13;

struct GlobalState {
  Real d;
  Real v[3];
  Real ppar;
  Real pperp;
  Real b[3];
};

struct GlobalFlux {
  Real d;
  Real mom[3];
  Real e;
  Real mu;
  Real emf[3];
};

void Fail(const std::string &label, Real got, Real expected) {
  std::cout << "CGL FOFC flux unit test failed: " << label
            << " got=" << got << " expected=" << expected
            << " abs_err=" << std::abs(got - expected) << std::endl;
  std::exit(EXIT_FAILURE);
}

void CheckClose(const std::string &label, Real got, Real expected) {
  Real scale = std::fmax(1.0, std::fmax(std::abs(got), std::abs(expected)));
  if (!std::isfinite(got) || !std::isfinite(expected) ||
      std::abs(got - expected) > kTol*scale) {
    Fail(label, got, expected);
  }
}

EOS_Data MakeCglEOS(bool passive=false) {
  EOS_Data eos;
  eos.gamma = 5.0/3.0;
  eos.iso_cs = 0.7;
  eos.nu_coll = 0.0;
  eos.lim_coll = 0.0;
  eos.is_ideal = true;
  eos.is_cgl = true;
  eos.passive = passive;
  eos.flim = false;
  eos.mlim = false;
  eos.coll = false;
  eos.backup_lim = false;
  eos.dfloor = 1.0e-12;
  eos.pfloor = 1.0e-12;
  eos.tfloor = 1.0e-12;
  eos.sfloor = 1.0e-12;
  eos.bfloor = 1.0e-10;
  eos.gamma_max = 1.0e12;
  eos.sigma_max = 1.0e12;
  return eos;
}

Real CglFastSpeed(const Real d, const Real pr, const Real pp, const Real bx,
                  const Real by, const Real bz, const Real bfloor) {
  // Match EOS_Data::IdealMHDFastSpeed exactly; this test isolates FOFC flux/EMF wiring.
  Real bprp2 = by*by + bz*bz;
  Real bx2 = bx*bx;
  Real b2 = bx2 + bprp2;
  if (b2 < bfloor*bfloor) {
    Real p = ONE_3RD*pr + TWO_3RDS*pp;
    Real asq = (5.0/3.0)*p;
    Real qsq = b2 + asq;
    Real tmp = b2 - asq;
    return std::sqrt(0.5*(qsq + std::sqrt(tmp*tmp + 4.0*asq*bprp2))/d);
  }
  Real bhatx2 = bx2/b2;
  Real qsq = b2 + 2.0*pp + (2.0*pr - pp)*bhatx2;
  Real radicand = qsq*qsq + 4.0*pp*pp*(1.0 - bhatx2)*bhatx2
                - 12.0*pr*pp*bhatx2*(2.0 - bhatx2)
                + 12.0*pr*pp*bhatx2*bhatx2 - 12.0*bx2*pr;
  return std::sqrt(0.5*(qsq + std::sqrt(std::abs(radicand))/d));
}

Real IsothermalFastSpeed(const Real d, const Real bx, const Real by, const Real bz,
                         const Real iso_cs) {
  Real asq = iso_cs*iso_cs*d;
  Real ct2 = by*by + bz*bz;
  Real qsq = bx*bx + ct2 + asq;
  Real tmp = bx*bx + ct2 - asq;
  return std::sqrt(0.5*(qsq + std::sqrt(tmp*tmp + 4.0*asq*ct2))/d);
}

MHDPrim1D LocalState(const GlobalState &s, int dir) {
  MHDPrim1D w;
  w.d = s.d;
  w.e = s.ppar;
  w.pp = s.pperp;
  if (dir == 0) {
    w.vx = s.v[0]; w.vy = s.v[1]; w.vz = s.v[2];
    w.by = s.b[1]; w.bz = s.b[2];
  } else if (dir == 1) {
    w.vx = s.v[1]; w.vy = s.v[2]; w.vz = s.v[0];
    w.by = s.b[2]; w.bz = s.b[0];
  } else {
    w.vx = s.v[2]; w.vy = s.v[0]; w.vz = s.v[1];
    w.by = s.b[0]; w.bz = s.b[1];
  }
  return w;
}

Real NormalField(const GlobalState &s, int dir) {
  return s.b[dir];
}

void ExpectedStateAndFlux(const MHDPrim1D &w, const Real bxi, const EOS_Data &eos,
                          MHDCons1D &u, MHDCons1D &f, Real &cf) {
  Real ppar = w.e;
  Real pperp = w.pp;
  Real bsq = bxi*bxi + w.by*w.by + w.bz*w.bz;
  Real bmag = std::sqrt(bsq);
  Real fh = 1.0;

  u.d = w.d;
  u.mx = w.d*w.vx;
  u.my = w.d*w.vy;
  u.mz = w.d*w.vz;
  u.by = w.by;
  u.bz = w.bz;

  if (bmag > eos.bfloor) {
    fh = 1.0 + (pperp - ppar)/bsq;
    u.mu = w.d*std::log(pperp/ppar*w.d*w.d/(bmag*bsq));
  } else {
    ppar = TWO_3RDS*pperp + ONE_3RD*ppar;
    pperp = ppar;
    u.mu = w.d*std::log(w.d*w.d/(eos.bfloor*eos.bfloor*eos.bfloor));
  }

  Real pB = 0.5*bsq;
  Real vsq = w.vx*w.vx + w.vy*w.vy + w.vz*w.vz;
  u.e = pperp + 0.5*ppar + 0.5*(w.d*vsq + bsq);

  f.d = w.d*w.vx;
  if (eos.passive) {
    f.mx = w.d*w.vx*w.vx + pB + eos.iso_cs*eos.iso_cs*w.d - bxi*bxi;
    f.my = w.d*w.vy*w.vx - bxi*w.by;
    f.mz = w.d*w.vz*w.vx - bxi*w.bz;
    cf = IsothermalFastSpeed(w.d, bxi, w.by, w.bz, eos.iso_cs);
  } else {
    f.mx = w.d*w.vx*w.vx + pB + pperp - bxi*bxi*fh;
    f.my = w.d*w.vy*w.vx - bxi*w.by*fh;
    f.mz = w.d*w.vz*w.vx - bxi*w.bz*fh;
    cf = CglFastSpeed(w.d, ppar, pperp, bxi, w.by, w.bz, eos.bfloor);
  }
  f.e = u.e*w.vx + w.vx*(pperp + pB)
      - bxi*(bxi*w.vx + w.by*w.vy + w.bz*w.vz)*fh;
  f.mu = u.mu*w.vx;
  f.by = w.by*w.vx - bxi*w.vy;
  f.bz = w.bz*w.vx - bxi*w.vz;
}

MHDCons1D ExpectedCglLLF(const MHDPrim1D &wl, const MHDPrim1D &wr,
                         const Real bxi, const EOS_Data &eos) {
  MHDCons1D ul, ur, fl, fr, flux;
  Real cl, cr;
  ExpectedStateAndFlux(wl, bxi, eos, ul, fl, cl);
  ExpectedStateAndFlux(wr, bxi, eos, ur, fr, cr);
  Real a = std::fmax(std::abs(wl.vx) + cl, std::abs(wr.vx) + cr);

  flux.d  = 0.5*(fl.d  + fr.d  - a*(ur.d  - ul.d ));
  flux.mx = 0.5*(fl.mx + fr.mx - a*(ur.mx - ul.mx));
  flux.my = 0.5*(fl.my + fr.my - a*(ur.my - ul.my));
  flux.mz = 0.5*(fl.mz + fr.mz - a*(ur.mz - ul.mz));
  flux.e  = 0.5*(fl.e  + fr.e  - a*(ur.e  - ul.e ));
  flux.mu = (flux.d >= 0.0) ? flux.d*(ul.mu/ul.d) : flux.d*(ur.mu/ur.d);
  flux.by = -0.5*(fl.by + fr.by - a*(ur.by - ul.by));
  flux.bz =  0.5*(fl.bz + fr.bz - a*(ur.bz - ul.bz));
  return flux;
}

GlobalFlux ToGlobalFlux(const MHDCons1D &local, int dir) {
  GlobalFlux g;
  g.d = local.d;
  g.e = local.e;
  g.mu = local.mu;
  g.mom[0] = 0.0; g.mom[1] = 0.0; g.mom[2] = 0.0;
  g.emf[0] = 0.0; g.emf[1] = 0.0; g.emf[2] = 0.0;

  if (dir == 0) {
    g.mom[0] = local.mx; g.mom[1] = local.my; g.mom[2] = local.mz;
    g.emf[2] = local.by;
    g.emf[1] = local.bz;
  } else if (dir == 1) {
    g.mom[1] = local.mx; g.mom[2] = local.my; g.mom[0] = local.mz;
    g.emf[0] = local.by;
    g.emf[2] = local.bz;
  } else {
    g.mom[2] = local.mx; g.mom[0] = local.my; g.mom[1] = local.mz;
    g.emf[1] = local.by;
    g.emf[0] = local.bz;
  }
  return g;
}

void CheckLocalFlux(const std::string &label, const MHDCons1D &got,
                    const MHDCons1D &expected) {
  CheckClose(label + ".d", got.d, expected.d);
  CheckClose(label + ".mx", got.mx, expected.mx);
  CheckClose(label + ".my", got.my, expected.my);
  CheckClose(label + ".mz", got.mz, expected.mz);
  CheckClose(label + ".e", got.e, expected.e);
  CheckClose(label + ".mu", got.mu, expected.mu);
  CheckClose(label + ".by", got.by, expected.by);
  CheckClose(label + ".bz", got.bz, expected.bz);
}

void CheckGlobalFlux(const std::string &label, const GlobalFlux &got,
                     const GlobalFlux &expected) {
  CheckClose(label + ".d", got.d, expected.d);
  CheckClose(label + ".m1", got.mom[0], expected.mom[0]);
  CheckClose(label + ".m2", got.mom[1], expected.mom[1]);
  CheckClose(label + ".m3", got.mom[2], expected.mom[2]);
  CheckClose(label + ".e", got.e, expected.e);
  CheckClose(label + ".mu", got.mu, expected.mu);
  CheckClose(label + ".e1", got.emf[0], expected.emf[0]);
  CheckClose(label + ".e2", got.emf[1], expected.emf[1]);
  CheckClose(label + ".e3", got.emf[2], expected.emf[2]);
}

void TestUnequalStates(const EOS_Data &eos, bool passive) {
  GlobalState left{1.30, {0.27, -0.19, 0.11}, 0.91, 0.57, {0.43, -0.31, 0.26}};
  GlobalState right{0.82, {-0.13, 0.23, -0.17}, 0.63, 1.04, {0.43, 0.18, -0.37}};
  if (passive) {
    left.ppar = 0.71; left.pperp = 0.88;
    right.ppar = 0.55; right.pperp = 0.92;
  }

  for (int dir=0; dir<3; ++dir) {
    MHDPrim1D wl = LocalState(left, dir);
    MHDPrim1D wr = LocalState(right, dir);
    Real bxi = NormalField(left, dir);

    MHDCons1D got;
    mhd::SingleStateLLF_CGL(wl, wr, bxi, eos, got);
    MHDCons1D expected = ExpectedCglLLF(wl, wr, bxi, eos);
    std::string label = passive ? "passive.dir" : "active.dir";
    label += std::to_string(dir + 1);
    CheckLocalFlux(label, got, expected);
    CheckGlobalFlux(label + ".global", ToGlobalFlux(got, dir), ToGlobalFlux(expected, dir));
  }
}

void TestEqualStateEMFSigns(const EOS_Data &eos) {
  GlobalState s{1.17, {0.41, -0.29, 0.23}, 0.76, 1.13, {0.37, -0.22, 0.58}};
  Real ephys[3];
  ephys[0] = s.v[2]*s.b[1] - s.v[1]*s.b[2];
  ephys[1] = s.v[0]*s.b[2] - s.v[2]*s.b[0];
  ephys[2] = s.v[1]*s.b[0] - s.v[0]*s.b[1];

  for (int dir=0; dir<3; ++dir) {
    MHDPrim1D wl = LocalState(s, dir);
    MHDPrim1D wr = LocalState(s, dir);
    Real bxi = NormalField(s, dir);

    MHDCons1D local;
    mhd::SingleStateLLF_CGL(wl, wr, bxi, eos, local);
    GlobalFlux global = ToGlobalFlux(local, dir);
    std::string label = "equal-state-emf.dir" + std::to_string(dir + 1);
    if (dir == 0) {
      CheckClose(label + ".e2", global.emf[1], ephys[1]);
      CheckClose(label + ".e3", global.emf[2], ephys[2]);
    } else if (dir == 1) {
      CheckClose(label + ".e1", global.emf[0], ephys[0]);
      CheckClose(label + ".e3", global.emf[2], ephys[2]);
    } else {
      CheckClose(label + ".e1", global.emf[0], ephys[0]);
      CheckClose(label + ".e2", global.emf[1], ephys[1]);
    }
  }
}

} // namespace

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (restart) return;

  EOS_Data active = MakeCglEOS(false);
  EOS_Data passive = MakeCglEOS(true);

  TestUnequalStates(active, false);
  TestUnequalStates(passive, true);
  TestEqualStateEMFSigns(active);

  std::cout << "CGL FOFC flux unit test passed" << std::endl;
}
