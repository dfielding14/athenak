//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_fofc_end_to_end_test.cpp
//! \brief End-to-end test for CGL FOFC mutation of flux and face-EMF Kokkos views.

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>

#include "athena.hpp"
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

  MHDCons1D expected;

  mhd::SingleStateLLF_CGL(XState(w, bcc, m, k, j, i-1), XState(w, bcc, m, k, j, i),
                          b1(m,k,j,i), eos, expected);
  CheckXFace("x-lower", f1, e3x1, e2x1, expected, m, k, j, i);

  mhd::SingleStateLLF_CGL(XState(w, bcc, m, k, j, i), XState(w, bcc, m, k, j, i+1),
                          b1(m,k,j,i+1), eos, expected);
  CheckXFace("x-upper", f1, e3x1, e2x1, expected, m, k, j, i+1);

  mhd::SingleStateLLF_CGL(YState(w, bcc, m, k, j-1, i), YState(w, bcc, m, k, j, i),
                          b2(m,k,j,i), eos, expected);
  CheckYFace("y-lower", f2, e1x2, e3x2, expected, m, k, j, i);

  mhd::SingleStateLLF_CGL(YState(w, bcc, m, k, j, i), YState(w, bcc, m, k, j+1, i),
                          b2(m,k,j+1,i), eos, expected);
  CheckYFace("y-upper", f2, e1x2, e3x2, expected, m, k, j+1, i);

  mhd::SingleStateLLF_CGL(ZState(w, bcc, m, k-1, j, i), ZState(w, bcc, m, k, j, i),
                          b3(m,k,j,i), eos, expected);
  CheckZFace("z-lower", f3, e2x3, e1x3, expected, m, k, j, i);

  mhd::SingleStateLLF_CGL(ZState(w, bcc, m, k, j, i), ZState(w, bcc, m, k+1, j, i),
                          b3(m,k+1,j,i), eos, expected);
  CheckZFace("z-upper", f3, e2x3, e1x3, expected, m, k+1, j, i);

  if (fofc(m,k,j,i)) {
    std::cout << "CGL FOFC end-to-end test failed: FOFC flag was not reset" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  ++validation_count;
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

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  pgen_final_func = FinalizeFOFCMutationTest;
  user_srcs_func = ValidateFOFCMutation;
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
  });

  par_for("cgl_fofc_e2e_b1", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie+1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    b0.x1f(m,k,j,i) = 0.43;
  });
  par_for("cgl_fofc_e2e_b2", DevExeSpace(), 0, nmb-1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    b0.x2f(m,k,j,i) = -0.31;
  });
  par_for("cgl_fofc_e2e_b3", DevExeSpace(), 0, nmb-1, ks, ke+1, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    b0.x3f(m,k,j,i) = 0.26;
  });

  pmhd->peos->PrimToCons(w0, bcc0, pmhd->u0, is, ie, js, je, ks, ke);
  Kokkos::deep_copy(pmhd->fofc, false);

  auto fofc = pmhd->fofc;
  int flag_i = is + 1;
  int flag_j = js + 1;
  int flag_k = ks + 1;
  par_for("cgl_fofc_e2e_flag", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(int) {
    fofc(0,flag_k,flag_j,flag_i) = true;
  });
}
