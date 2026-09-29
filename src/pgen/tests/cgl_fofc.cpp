//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_fofc.cpp
//! \brief End-to-end test for CGL FOFC mutation of flux and face-EMF Kokkos views.

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <utility>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "mhd/rsolvers/llf_mhd_singlestate.hpp"
#include "parameter_input.hpp"
#include "pgen/pgen.hpp"
#include "reconstruct/ppm.hpp"
#include "reconstruct/wenoz.hpp"

namespace {

constexpr Real kTol = 1.0e-13;
int validation_count = 0;

void Fail(const std::string &label, Real got, Real expected) {
  std::cout << "CGL FOFC end-to-end test failed: " << label
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

template <typename ViewType>
auto HostCopy(const ViewType &view) {
  return Kokkos::create_mirror_view_and_copy(HostMemSpace(), view);
}

void CheckReconstructionFloors(const EOS_Data &input_eos) {
  EOS_Data eos = input_eos;
  eos.pfloor = 1.0;
  DvceArray5D<Real> q("cgl_reconstruction_states", 1, 6, 7, 7, 7);
  DvceArray2D<Real> result("cgl_reconstruction_faces", 12, 4);
  Kokkos::deep_copy(q, 0.25);
  const size_t scratch_size = 2*ScrArray2D<Real>::shmem_size(6, 7);
  par_for_outer("cgl_reconstruction_floor_check", DevExeSpace(), scratch_size, 0,
                0, 11, KOKKOS_LAMBDA(TeamMember_t member, const int test) {
    EOS_Data local = eos;
    local.is_cgl = test < 6;
    const int method = test % 6;
    const int dir = method % 3;
    ScrArray2D<Real> ql(member.team_scratch(0), 6, 7);
    ScrArray2D<Real> qr(member.team_scratch(0), 6, 7);
    if (method == 0) {
      PiecewiseParabolicX1(member, local, true, true, 0, 3, 3, 3, 3, q, ql, qr);
    } else if (method == 1) {
      PiecewiseParabolicX2(member, local, true, true, 0, 3, 3, 3, 3, q, ql, qr);
    } else if (method == 2) {
      PiecewiseParabolicX3(member, local, true, true, 0, 3, 3, 3, 3, q, ql, qr);
    } else if (method == 3) {
      WENOZX1(member, local, true, 0, 3, 3, 3, 3, q, ql, qr);
    } else if (method == 4) {
      WENOZX2(member, local, true, 0, 3, 3, 3, 3, q, ql, qr);
    } else {
      WENOZX3(member, local, true, 0, 3, 3, 3, 3, q, ql, qr);
    }
    member.team_barrier();
    par_for_inner(member, 0, 0, [&](const int) {
      result(test,0) = ql(IPR, dir == 0 ? 4 : 3);
      result(test,1) = qr(IPR, 3);
      result(test,2) = ql(IPP, dir == 0 ? 4 : 3);
      result(test,3) = qr(IPP, 3);
    });
  });
  const auto faces = HostCopy(result);
  for (int test=0; test<12; ++test) {
    const Real expected = test < 6 ? eos.pfloor : eos.pfloor/(eos.gamma - 1.0);
    CheckClose("parallel reconstruction floor left", faces(test,0), expected);
    CheckClose("parallel reconstruction floor right", faces(test,1), expected);
    CheckClose("perpendicular reconstruction floor left", faces(test,2),
                test < 6 ? eos.pfloor : 0.25);
    CheckClose("perpendicular reconstruction floor right", faces(test,3),
                test < 6 ? eos.pfloor : 0.25);
  }
}

MHDPrim1D XState(const decltype(HostCopy(std::declval<DvceArray5D<Real>>())) &w,
                 const decltype(HostCopy(std::declval<DvceArray5D<Real>>())) &bcc,
                 int m, int k, int j, int i) {
  MHDPrim1D s;
  s.d = w(m,IDN,k,j,i);
  s.vx = w(m,IVX,k,j,i);
  s.vy = w(m,IVY,k,j,i);
  s.vz = w(m,IVZ,k,j,i);
  s.e = w(m,IPR,k,j,i);
  s.pp = w(m,IPP,k,j,i);
  s.by = bcc(m,IBY,k,j,i);
  s.bz = bcc(m,IBZ,k,j,i);
  return s;
}

MHDPrim1D YState(const decltype(HostCopy(std::declval<DvceArray5D<Real>>())) &w,
                 const decltype(HostCopy(std::declval<DvceArray5D<Real>>())) &bcc,
                 int m, int k, int j, int i) {
  MHDPrim1D s;
  s.d = w(m,IDN,k,j,i);
  s.vx = w(m,IVY,k,j,i);
  s.vy = w(m,IVZ,k,j,i);
  s.vz = w(m,IVX,k,j,i);
  s.e = w(m,IPR,k,j,i);
  s.pp = w(m,IPP,k,j,i);
  s.by = bcc(m,IBZ,k,j,i);
  s.bz = bcc(m,IBX,k,j,i);
  return s;
}

MHDPrim1D ZState(const decltype(HostCopy(std::declval<DvceArray5D<Real>>())) &w,
                 const decltype(HostCopy(std::declval<DvceArray5D<Real>>())) &bcc,
                 int m, int k, int j, int i) {
  MHDPrim1D s;
  s.d = w(m,IDN,k,j,i);
  s.vx = w(m,IVZ,k,j,i);
  s.vy = w(m,IVX,k,j,i);
  s.vz = w(m,IVY,k,j,i);
  s.e = w(m,IPR,k,j,i);
  s.pp = w(m,IPP,k,j,i);
  s.by = bcc(m,IBX,k,j,i);
  s.bz = bcc(m,IBY,k,j,i);
  return s;
}

template <typename FluxView, typename E3View, typename E2View>
void CheckXFace(const std::string &label, const FluxView &flx, const E3View &e3x1,
                const E2View &e2x1, const MHDCons1D &expected,
                int m, int k, int j, int i) {
  CheckClose(label + ".IDN", flx(m,IDN,k,j,i), expected.d);
  CheckClose(label + ".IM1", flx(m,IM1,k,j,i), expected.mx);
  CheckClose(label + ".IM2", flx(m,IM2,k,j,i), expected.my);
  CheckClose(label + ".IM3", flx(m,IM3,k,j,i), expected.mz);
  CheckClose(label + ".IEN", flx(m,IEN,k,j,i), expected.e);
  CheckClose(label + ".IMU", flx(m,IMU,k,j,i), expected.mu);
  CheckClose(label + ".E3", e3x1(m,k,j,i), expected.by);
  CheckClose(label + ".E2", e2x1(m,k,j,i), expected.bz);
}

template <typename FluxView, typename E1View, typename E3View>
void CheckYFace(const std::string &label, const FluxView &flx, const E1View &e1x2,
                const E3View &e3x2, const MHDCons1D &expected,
                int m, int k, int j, int i) {
  CheckClose(label + ".IDN", flx(m,IDN,k,j,i), expected.d);
  CheckClose(label + ".IM2", flx(m,IM2,k,j,i), expected.mx);
  CheckClose(label + ".IM3", flx(m,IM3,k,j,i), expected.my);
  CheckClose(label + ".IM1", flx(m,IM1,k,j,i), expected.mz);
  CheckClose(label + ".IEN", flx(m,IEN,k,j,i), expected.e);
  CheckClose(label + ".IMU", flx(m,IMU,k,j,i), expected.mu);
  CheckClose(label + ".E1", e1x2(m,k,j,i), expected.by);
  CheckClose(label + ".E3", e3x2(m,k,j,i), expected.bz);
}

template <typename FluxView, typename E2View, typename E1View>
void CheckZFace(const std::string &label, const FluxView &flx, const E2View &e2x3,
                const E1View &e1x3, const MHDCons1D &expected,
                int m, int k, int j, int i) {
  CheckClose(label + ".IDN", flx(m,IDN,k,j,i), expected.d);
  CheckClose(label + ".IM3", flx(m,IM3,k,j,i), expected.mx);
  CheckClose(label + ".IM1", flx(m,IM1,k,j,i), expected.my);
  CheckClose(label + ".IM2", flx(m,IM2,k,j,i), expected.mz);
  CheckClose(label + ".IEN", flx(m,IEN,k,j,i), expected.e);
  CheckClose(label + ".IMU", flx(m,IMU,k,j,i), expected.mu);
  CheckClose(label + ".E2", e2x3(m,k,j,i), expected.by);
  CheckClose(label + ".E1", e1x3(m,k,j,i), expected.bz);
}

template <typename FluxView>
void CheckXPressure(const std::string &label, const FluxView &flx,
                    const MHDCons1D &pressure, const MHDCons1D &anisotropic,
                    int m, int k, int j, int i) {
  CheckClose(label + ".PX", flx(m,mhd::ICGLPressureX,k,j,i), pressure.mx);
  CheckClose(label + ".PY", flx(m,mhd::ICGLPressureY,k,j,i), pressure.my);
  CheckClose(label + ".PZ", flx(m,mhd::ICGLPressureZ,k,j,i), pressure.mz);
  CheckClose(label + ".AX", flx(m,mhd::ICGLAnisPressureX,k,j,i), anisotropic.mx);
  CheckClose(label + ".AY", flx(m,mhd::ICGLAnisPressureY,k,j,i), anisotropic.my);
  CheckClose(label + ".AZ", flx(m,mhd::ICGLAnisPressureZ,k,j,i), anisotropic.mz);
}

template <typename FluxView>
void CheckYPressure(const std::string &label, const FluxView &flx,
                    const MHDCons1D &pressure, const MHDCons1D &anisotropic,
                    int m, int k, int j, int i) {
  CheckClose(label + ".PX", flx(m,mhd::ICGLPressureX,k,j,i), pressure.mz);
  CheckClose(label + ".PY", flx(m,mhd::ICGLPressureY,k,j,i), pressure.mx);
  CheckClose(label + ".PZ", flx(m,mhd::ICGLPressureZ,k,j,i), pressure.my);
  CheckClose(label + ".AX", flx(m,mhd::ICGLAnisPressureX,k,j,i), anisotropic.mz);
  CheckClose(label + ".AY", flx(m,mhd::ICGLAnisPressureY,k,j,i), anisotropic.mx);
  CheckClose(label + ".AZ", flx(m,mhd::ICGLAnisPressureZ,k,j,i), anisotropic.my);
}

template <typename FluxView>
void CheckZPressure(const std::string &label, const FluxView &flx,
                    const MHDCons1D &pressure, const MHDCons1D &anisotropic,
                    int m, int k, int j, int i) {
  CheckClose(label + ".PX", flx(m,mhd::ICGLPressureX,k,j,i), pressure.my);
  CheckClose(label + ".PY", flx(m,mhd::ICGLPressureY,k,j,i), pressure.mz);
  CheckClose(label + ".PZ", flx(m,mhd::ICGLPressureZ,k,j,i), pressure.mx);
  CheckClose(label + ".AX", flx(m,mhd::ICGLAnisPressureX,k,j,i), anisotropic.my);
  CheckClose(label + ".AY", flx(m,mhd::ICGLAnisPressureY,k,j,i), anisotropic.mz);
  CheckClose(label + ".AZ", flx(m,mhd::ICGLAnisPressureZ,k,j,i), anisotropic.mx);
}

void ValidateFOFCMutation(Mesh *pm, const Real /*bdt*/) {
  auto *pmhd = pm->pmb_pack->pmhd;
  if (pmhd == nullptr) {
    std::cout << "CGL FOFC end-to-end test requires MHD" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  auto w = HostCopy(pmhd->w0);
  auto bcc = HostCopy(pmhd->bcc0);
  auto b1 = HostCopy(pmhd->b0.x1f);
  auto b2 = HostCopy(pmhd->b0.x2f);
  auto b3 = HostCopy(pmhd->b0.x3f);
  auto f1 = HostCopy(pmhd->uflx.x1f);
  auto f2 = HostCopy(pmhd->uflx.x2f);
  auto f3 = HostCopy(pmhd->uflx.x3f);
  if (!pmhd->record_cgl_pressure_work) {
    std::cout << "CGL FOFC end-to-end test requires retained CGL pressure traction"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  auto pf1 = HostCopy(pmhd->cgl_pflux.x1f);
  auto pf2 = HostCopy(pmhd->cgl_pflux.x2f);
  auto pf3 = HostCopy(pmhd->cgl_pflux.x3f);
  auto e3x1 = HostCopy(pmhd->e3x1);
  auto e2x1 = HostCopy(pmhd->e2x1);
  auto e1x2 = HostCopy(pmhd->e1x2);
  auto e3x2 = HostCopy(pmhd->e3x2);
  auto e2x3 = HostCopy(pmhd->e2x3);
  auto e1x3 = HostCopy(pmhd->e1x3);
  auto fofc = HostCopy(pmhd->fofc);

  int m = 0;
  auto &indcs = pm->mb_indcs;
  int i = indcs.is + 1;
  int j = indcs.js + 1;
  int k = indcs.ks + 1;
  const EOS_Data &eos = pmhd->peos->eos_data;

  MHDCons1D expected, pressure, anisotropic;

  mhd::SingleStateLLF_CGL(XState(w, bcc, m, k, j, i-1), XState(w, bcc, m, k, j, i),
                          b1(m,k,j,i), eos, expected);
  mhd::SingleStateLLF_CGLPressureTraction(
      XState(w, bcc, m, k, j, i-1), XState(w, bcc, m, k, j, i),
      b1(m,k,j,i), eos, pressure, anisotropic);
  CheckXFace("x-lower", f1, e3x1, e2x1, expected, m, k, j, i);
  CheckXPressure("x-lower", pf1, pressure, anisotropic, m, k, j, i);

  mhd::SingleStateLLF_CGL(XState(w, bcc, m, k, j, i), XState(w, bcc, m, k, j, i+1),
                          b1(m,k,j,i+1), eos, expected);
  mhd::SingleStateLLF_CGLPressureTraction(
      XState(w, bcc, m, k, j, i), XState(w, bcc, m, k, j, i+1),
      b1(m,k,j,i+1), eos, pressure, anisotropic);
  CheckXFace("x-upper", f1, e3x1, e2x1, expected, m, k, j, i+1);
  CheckXPressure("x-upper", pf1, pressure, anisotropic, m, k, j, i+1);

  mhd::SingleStateLLF_CGL(YState(w, bcc, m, k, j-1, i), YState(w, bcc, m, k, j, i),
                          b2(m,k,j,i), eos, expected);
  mhd::SingleStateLLF_CGLPressureTraction(
      YState(w, bcc, m, k, j-1, i), YState(w, bcc, m, k, j, i),
      b2(m,k,j,i), eos, pressure, anisotropic);
  CheckYFace("y-lower", f2, e1x2, e3x2, expected, m, k, j, i);
  CheckYPressure("y-lower", pf2, pressure, anisotropic, m, k, j, i);

  mhd::SingleStateLLF_CGL(YState(w, bcc, m, k, j, i), YState(w, bcc, m, k, j+1, i),
                          b2(m,k,j+1,i), eos, expected);
  mhd::SingleStateLLF_CGLPressureTraction(
      YState(w, bcc, m, k, j, i), YState(w, bcc, m, k, j+1, i),
      b2(m,k,j+1,i), eos, pressure, anisotropic);
  CheckYFace("y-upper", f2, e1x2, e3x2, expected, m, k, j+1, i);
  CheckYPressure("y-upper", pf2, pressure, anisotropic, m, k, j+1, i);

  mhd::SingleStateLLF_CGL(ZState(w, bcc, m, k-1, j, i), ZState(w, bcc, m, k, j, i),
                          b3(m,k,j,i), eos, expected);
  mhd::SingleStateLLF_CGLPressureTraction(
      ZState(w, bcc, m, k-1, j, i), ZState(w, bcc, m, k, j, i),
      b3(m,k,j,i), eos, pressure, anisotropic);
  CheckZFace("z-lower", f3, e2x3, e1x3, expected, m, k, j, i);
  CheckZPressure("z-lower", pf3, pressure, anisotropic, m, k, j, i);

  mhd::SingleStateLLF_CGL(ZState(w, bcc, m, k, j, i), ZState(w, bcc, m, k+1, j, i),
                          b3(m,k+1,j,i), eos, expected);
  mhd::SingleStateLLF_CGLPressureTraction(
      ZState(w, bcc, m, k, j, i), ZState(w, bcc, m, k+1, j, i),
      b3(m,k+1,j,i), eos, pressure, anisotropic);
  CheckZFace("z-upper", f3, e2x3, e1x3, expected, m, k+1, j, i);
  CheckZPressure("z-upper", pf3, pressure, anisotropic, m, k+1, j, i);

  if (fofc(m,k,j,i)) {
    std::cout << "CGL FOFC end-to-end test failed: FOFC flag was not reset" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  ++validation_count;
}

// Exercise the real only_testfloors path before C2P can hide a nonfinite A.
void CheckNonfiniteDetector(Mesh *pm) {
  auto *pmhd = pm->pmb_pack->pmhd;
  const int i = pm->mb_indcs.is, j = pm->mb_indcs.js, k = pm->mb_indcs.ks;
  auto u = HostCopy(pmhd->u0);
  const int count_before = pm->ecounter.nfofc;
  for (const int n : {IEN, IAN}) {
    const Real original = u(0,n,k,j,i);
    for (const Real bad : {std::numeric_limits<Real>::quiet_NaN(),
                           std::numeric_limits<Real>::infinity(),
                           -std::numeric_limits<Real>::infinity()}) {
      u(0,n,k,j,i) = bad;
      Kokkos::deep_copy(pmhd->utest, u);
      Kokkos::deep_copy(pmhd->fofc, false);
      pmhd->peos->ConsToPrim(pmhd->utest, pmhd->b0, pmhd->w0, pmhd->bcc0,
                              true, i, i, j, j, k, k);
      const auto flags = HostCopy(pmhd->fofc);
      if (!flags(0,k,j,i)) Fail("nonfinite FOFC input", 0.0, 1.0);
    }
    u(0,n,k,j,i) = original;
  }
  Kokkos::deep_copy(pmhd->fofc, false);
  pm->ecounter.nfofc = count_before;
}

void ValidatePressureStep(Mesh *pm, const Real /*bdt*/) {
  auto *pmhd = pm->pmb_pack->pmhd;
  const auto u = HostCopy(pmhd->u0);
  const auto w = HostCopy(pmhd->w0);
  auto &indcs = pm->mb_indcs;
  for (int i=indcs.is; i<=indcs.ie; ++i) {
    for (int n=0; n<pmhd->nmhd; ++n) {
      if (!std::isfinite(u(0,n,indcs.ks,indcs.js,i)) ||
          !std::isfinite(w(0,n,indcs.ks,indcs.js,i))) {
        Fail("pressure-step finite state", 0.0, 1.0);
      }
    }
    if (!(w(0,IPR,indcs.ks,indcs.js,i) > 0.0) ||
        !(w(0,IPP,indcs.ks,indcs.js,i) > 0.0)) {
      Fail("pressure-step positive pressures", 0.0, 1.0);
    }
  }
}

void FinalizePressureStep(ParameterInput *, Mesh *pm) {
  ValidatePressureStep(pm, 0.0);
  CheckClose("pressure-step completed cycles", pm->ncycle, 50.0);
  std::cout << "CGL pressure-step test passed after 50 cycles; FOFC counts are in "
            << "the event log" << std::endl;
}

void FinalizeFOFCMutationTest(ParameterInput *, Mesh *) {
  if (validation_count != 1) {
    std::cout << "CGL FOFC end-to-end test failed: expected one validation, got "
              << validation_count << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::cout << "CGL FOFC end-to-end test passed" << std::endl;
}

} // namespace

void ProblemGenerator::CGLFOFC(ParameterInput *pin, const bool restart) {
  const bool pressure_step =
      pin->GetOrAddString("problem", "test_mode", "flux_mutation") == "pressure_step";
  pgen_final_func = pressure_step ? FinalizePressureStep : FinalizeFOFCMutationTest;
  user_srcs_func = pressure_step ? ValidatePressureStep : ValidateFOFCMutation;
  if (restart) return;

  auto *pmbp = pmy_mesh_->pmb_pack;
  auto *pmhd = pmbp->pmhd;
  if (pmhd == nullptr || !pmhd->peos->eos_data.is_cgl) {
    std::cout << "CGL FOFC end-to-end test requires <mhd>/eos = cgl" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  auto &indcs = pmy_mesh_->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb = pmbp->nmb_thispack;
  auto w0 = pmhd->w0;
  auto bcc0 = pmhd->bcc0;
  auto b0 = pmhd->b0;

  par_for("cgl_fofc_e2e_init", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real x = static_cast<Real>(i - is);
    Real y = static_cast<Real>(j - js);
    Real z = static_cast<Real>(k - ks);
    w0(m,IDN,k,j,i) = 1.0 + 0.04*x + 0.03*y + 0.02*z;
    w0(m,IVX,k,j,i) = 0.12 + 0.01*x - 0.015*y + 0.006*z;
    w0(m,IVY,k,j,i) = -0.08 + 0.012*x + 0.007*y - 0.004*z;
    w0(m,IVZ,k,j,i) = 0.05 - 0.006*x + 0.011*y + 0.009*z;
    w0(m,IPR,k,j,i) = 0.74 + 0.025*x + 0.017*y + 0.013*z;
    w0(m,IPP,k,j,i) = 1.08 + 0.019*x + 0.021*y + 0.015*z;
    bcc0(m,IBX,k,j,i) = 0.43;
    bcc0(m,IBY,k,j,i) = -0.31;
    bcc0(m,IBZ,k,j,i) = 0.26;
    if (pressure_step) {
      // A three-cell peak has 1000:1 perpendicular-pressure contrast.
      const int distance = abs(i - (is + ie)/2);
      w0(m,IDN,k,j,i) = 1.0;
      w0(m,IVX,k,j,i) = 1.0;
      w0(m,IVY,k,j,i) = 0.0;
      w0(m,IVZ,k,j,i) = 0.0;
      w0(m,IPR,k,j,i) = 1.0;
      w0(m,IPP,k,j,i) = distance == 0 ? 1.0 : (distance == 1 ? 0.25 : 0.001);
      bcc0(m,IBX,k,j,i) = 2.0;
      bcc0(m,IBY,k,j,i) = 0.0;
      bcc0(m,IBZ,k,j,i) = 0.0;
    }
  });

  par_for("cgl_fofc_e2e_b1", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie+1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    b0.x1f(m,k,j,i) = pressure_step ? 2.0 : 0.43;
  });
  par_for("cgl_fofc_e2e_b2", DevExeSpace(), 0, nmb-1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    b0.x2f(m,k,j,i) = pressure_step ? 0.0 : -0.31;
  });
  par_for("cgl_fofc_e2e_b3", DevExeSpace(), 0, nmb-1, ks, ke+1, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    b0.x3f(m,k,j,i) = pressure_step ? 0.0 : 0.26;
  });

  pmhd->peos->PrimToCons(w0, bcc0, pmhd->u0, is, ie, js, je, ks, ke);
  CheckReconstructionFloors(pmhd->peos->eos_data);
  CheckNonfiniteDetector(pmy_mesh_);
  Kokkos::deep_copy(pmhd->fofc, false);
  if (pressure_step) return;

  auto fofc = pmhd->fofc;
  int flag_i = is + 1;
  int flag_j = js + 1;
  int flag_k = ks + 1;
  par_for("cgl_fofc_e2e_flag", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(int) {
    fofc(0,flag_k,flag_j,flag_i) = true;
  });
}
