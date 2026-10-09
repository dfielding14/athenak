//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_fofc.cpp
//! \brief End-to-end test for CGL FOFC mutation of flux and face-EMF Kokkos views.

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>
#include <utility>

#include "athena.hpp"
#include "driver/driver.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "mhd/rsolvers/llf_mhd_singlestate.hpp"
#include "parameter_input.hpp"
#include "pgen/pgen.hpp"

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

template <typename ViewType>
auto Snapshot(const ViewType &view) {
  auto host = Kokkos::create_mirror(view);
  Kokkos::deep_copy(host, view);
  return host;
}

template <typename ViewType, typename HostView>
void CheckUnchanged(const std::string &label, const ViewType &view,
                    const HostView &before) {
  const auto after = HostCopy(view);
  if (std::memcmp(before.data(), after.data(),
                  before.size()*sizeof(typename ViewType::non_const_value_type)) != 0) {
    Fail(label, 0.0, 1.0);
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

void CheckCorrectedCell(Mesh *pm, const int k, const int j, const int i) {
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
}

void ValidateFOFCMutation(Mesh *pm, const Real /*bdt*/) {
  const auto &indcs = pm->mb_indcs;
  CheckCorrectedCell(pm, indcs.ks+1, indcs.js+1, indcs.is+1);
  ++validation_count;
}

void CheckNonfinitePredictorFluxes(ParameterInput *pin, Mesh *pm) {
  auto *pmhd = pm->pmb_pack->pmhd;
  const auto saved_u0 = Snapshot(pmhd->u0), saved_w0 = Snapshot(pmhd->w0);
  const auto saved_bcc = Snapshot(pmhd->bcc0), saved_u1 = Snapshot(pmhd->u1);
  const auto saved_bx = Snapshot(pmhd->b0.x1f), saved_by = Snapshot(pmhd->b0.x2f);
  const auto saved_bz = Snapshot(pmhd->b0.x3f), saved_b1x = Snapshot(pmhd->b1.x1f);
  const auto saved_b1y = Snapshot(pmhd->b1.x2f), saved_b1z = Snapshot(pmhd->b1.x3f);
  const auto saved_flags = Snapshot(pmhd->fofc);
  const auto saved_events = pm->ecounter;
  const Real saved_dt = pm->dt;
  Kokkos::Timer timer;
  Driver driver(pin, pm, -1.0, &timer);
  pm->dt = 1.0e-4;
  Kokkos::deep_copy(pmhd->u1, pmhd->u0);
  Kokkos::deep_copy(pmhd->b1.x1f, pmhd->b0.x1f);
  Kokkos::deep_copy(pmhd->b1.x2f, pmhd->b0.x2f);
  Kokkos::deep_copy(pmhd->b1.x3f, pmhd->b0.x3f);

  const DvceArray5D<Real> flux[] = {pmhd->uflx.x1f, pmhd->uflx.x2f, pmhd->uflx.x3f};
  const DvceArray5D<Real> traction[] = {pmhd->cgl_pflux.x1f, pmhd->cgl_pflux.x2f,
                                      pmhd->cgl_pflux.x3f};
  const DvceArray4D<Real> emf1[] = {pmhd->e3x1, pmhd->e1x2, pmhd->e2x3};
  const DvceArray4D<Real> emf2[] = {pmhd->e2x1, pmhd->e3x2, pmhd->e1x3};
  const Real invalid[] = {std::numeric_limits<Real>::quiet_NaN(),
                         std::numeric_limits<Real>::infinity(),
                         -std::numeric_limits<Real>::infinity()};
  const auto &indcs = pm->mb_indcs;
  const int lower[] = {indcs.is+1, indcs.js+1, indcs.ks+1};
  for (int direction=0; direction<3; ++direction) {
    int upper[] = {lower[0], lower[1], lower[2]};
    ++upper[direction];
    const int *cells[] = {lower, upper};
    for (const int field : {IEN, IAN}) {
      for (int bad=0; bad<3; ++bad) {
        const std::string label = "nonfinite predictor dir=" + std::to_string(direction)
            + " field=" + std::to_string(field) + " bad=" + std::to_string(bad);
        for (int d=0; d<3; ++d) {
          Kokkos::deep_copy(flux[d], 0.0);
          Kokkos::deep_copy(traction[d], 0.0);
          Kokkos::deep_copy(emf1[d], 0.0);
          Kokkos::deep_copy(emf2[d], 0.0);
        }
        Kokkos::deep_copy(pmhd->fofc, false);
        Kokkos::deep_copy(Kokkos::subview(flux[direction], 0, field,
            upper[2], upper[1], upper[0]), invalid[bad]);
        const int before = pm->ecounter.nfofc;
        pmhd->FOFC(&driver, 1);
        CheckClose(label + " detection count", pm->ecounter.nfofc-before, 2.0);
        CheckCorrectedCell(pm, lower[2], lower[1], lower[0]);
        CheckCorrectedCell(pm, upper[2], upper[1], upper[0]);

        // Only faces touching either detected cell may be replaced.
        for (int d=0; d<3; ++d) {
          const auto f = HostCopy(flux[d]), p = HostCopy(traction[d]);
          const auto e1 = HostCopy(emf1[d]), e2 = HostCopy(emf2[d]);
          for (int m=0; m<f.extent_int(0); ++m) {
            for (int k=0; k<f.extent_int(2); ++k) {
              for (int j=0; j<f.extent_int(3); ++j) {
                for (int i=0; i<f.extent_int(4); ++i) {
                  const int face[] = {i, j, k};
                  bool touched = false;
                  for (const int *cell : cells) {
                    bool adjacent = m == 0;
                    for (int axis=0; axis<3; ++axis) {
                      const int delta = face[axis] - cell[axis];
                      adjacent = adjacent && (axis == d ? delta == 0 || delta == 1
                                                        : delta == 0);
                    }
                    touched = touched || adjacent;
                  }
                  if (touched) continue;
                  for (int n=0; n<f.extent_int(1); ++n) {
                    if (f(m,n,k,j,i) != 0.0) Fail(label + " unrelated flux", 1.0, 0.0);
                  }
                  for (int n=0; n<p.extent_int(1); ++n) {
                    if (p(m,n,k,j,i) != 0.0) Fail(label + " unrelated traction", 1.0, 0.0);
                  }
                  // EMF scratch views omit the extra outer plane of the flux views.
                  if (k < e1.extent_int(1) && j < e1.extent_int(2) &&
                      i < e1.extent_int(3) &&
                      (e1(m,k,j,i) != 0.0 || e2(m,k,j,i) != 0.0)) {
                    Fail(label + " unrelated EMF", 1.0, 0.0);
                  }
                }
              }
            }
          }
        }
        const auto flags = HostCopy(pmhd->fofc);
        for (size_t c=0; c<flags.size(); ++c) {
          if (flags.data()[c]) Fail(label + " flag reset", 1.0, 0.0);
        }
        CheckUnchanged(label + " u0", pmhd->u0, saved_u0);
        CheckUnchanged(label + " w0", pmhd->w0, saved_w0);
        CheckUnchanged(label + " bcc0", pmhd->bcc0, saved_bcc);
        CheckUnchanged(label + " b0.x1f", pmhd->b0.x1f, saved_bx);
        CheckUnchanged(label + " b0.x2f", pmhd->b0.x2f, saved_by);
        CheckUnchanged(label + " b0.x3f", pmhd->b0.x3f, saved_bz);
      }
    }
  }
  Kokkos::deep_copy(pmhd->u1, saved_u1);
  Kokkos::deep_copy(pmhd->b1.x1f, saved_b1x);
  Kokkos::deep_copy(pmhd->b1.x2f, saved_b1y);
  Kokkos::deep_copy(pmhd->b1.x3f, saved_b1z);
  Kokkos::deep_copy(pmhd->fofc, saved_flags);
  pm->dt = saved_dt;
  pm->ecounter = saved_events;
  std::cout << "CGL nonfinite predictor-to-FOFC checks passed (18 cases)" << std::endl;
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

void FinalizeFOFCMutationTest(ParameterInput *pin, Mesh *pm) {
  if (validation_count != 1) {
    std::cout << "CGL FOFC end-to-end test failed: expected one validation, got "
              << validation_count << std::endl;
    std::exit(EXIT_FAILURE);
  }
  CheckNonfinitePredictorFluxes(pin, pm);
  std::cout << "CGL FOFC end-to-end test passed" << std::endl;
}

} // namespace

void ProblemGenerator::CGLFOFC(ParameterInput *pin, const bool restart) {
  const std::string mode = pin->GetOrAddString("problem", "test_mode", "flux_mutation");
  const bool pressure_step = mode == "pressure_step";
  if (!pressure_step && mode != "flux_mutation") {
    std::cout << "CGL FOFC test requires test_mode = flux_mutation or pressure_step"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
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
  if (pressure_step && (pmy_mesh_->nmb_total != 1 || indcs.nx2 != 1 || indcs.nx3 != 1)) {
    std::cout << "CGL pressure-step test requires a single 1D MeshBlock" << std::endl;
    std::exit(EXIT_FAILURE);
  }
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
