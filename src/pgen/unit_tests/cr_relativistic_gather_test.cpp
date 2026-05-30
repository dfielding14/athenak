//========================================================================================
// AthenaK astrophysical fluid dynamics & numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cr_relativistic_gather_test.cpp
//! \brief Host/device parity and analytical tests for passive relativistic CR gathers.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "parameter_input.hpp"
#include "particles/trilinear_gather.hpp"
#include "pgen/pgen.hpp"

namespace {

using particles::CellCenteredTrilinearStencil;
using particles::ConstructCellCenteredTrilinearStencil;
using particles::ConstructIdealCE;
using particles::GatherNewtonianCellCenteredFields;
using particles::NewtonianCellCenteredFields;
using particles::ParticleVector3;

struct GatherProbe {
  int m;
  Real x1, x2, x3;
};

struct GatherResult {
  NewtonianCellCenteredFields fields;
  ParticleVector3 cE;
  CellCenteredTrilinearStencil stencil;
  Real cE_dot_b;
  Real weight_sum;
};

struct HostBlockSize {
  HostArray1D<RegionSize> d_view;
};

enum class TestProfile {uniform, affine, nonlinear};

[[noreturn]] void Fail(const std::string &message) {
  std::cout << "CR relativistic gather helper test failed: " << message << std::endl;
  std::exit(EXIT_FAILURE);
}

void Require(const bool condition, const std::string &message) {
  if (!condition) {Fail(message);}
}

Real VectorNorm(const ParticleVector3 &value) {
  return std::sqrt(value.x*value.x + value.y*value.y + value.z*value.z);
}

bool NearlyEqual(const Real a, const Real b, const Real tolerance) {
  if (!std::isfinite(a) || !std::isfinite(b)) {return a == b;}
  Real scale = std::max(static_cast<Real>(1.0), std::max(std::fabs(a), std::fabs(b)));
  return std::fabs(a - b) <= tolerance*scale;
}

void RequireVectorNear(const ParticleVector3 &actual, const ParticleVector3 &expected,
                       const Real tolerance, const std::string &message) {
  Require(NearlyEqual(actual.x, expected.x, tolerance) &&
          NearlyEqual(actual.y, expected.y, tolerance) &&
          NearlyEqual(actual.z, expected.z, tolerance), message);
}

KOKKOS_INLINE_FUNCTION
ParticleVector3 UniformU() {
  return {0.25, -0.5, 0.75};
}

KOKKOS_INLINE_FUNCTION
ParticleVector3 UniformB() {
  return {-0.4, 0.6, 1.2};
}

KOKKOS_INLINE_FUNCTION
ParticleVector3 AffineU(const Real x1, const Real x2, const Real x3) {
  return {0.2 + 0.3*x1 - 0.4*x2 + 0.5*x3,
          -0.1 + 0.6*x1 + 0.2*x2 - 0.3*x3,
          0.7 - 0.2*x1 + 0.8*x2 + 0.1*x3};
}

KOKKOS_INLINE_FUNCTION
ParticleVector3 AffineB(const Real x1, const Real x2, const Real x3) {
  return {-0.3 + 0.4*x1 + 0.2*x2 - 0.1*x3,
          0.5 - 0.3*x1 + 0.6*x2 + 0.2*x3,
          0.9 + 0.1*x1 - 0.5*x2 + 0.7*x3};
}

KOKKOS_INLINE_FUNCTION
ParticleVector3 NonlinearU(const Real x1, const Real x2, const Real x3) {
  const Real twopi = 6.283185307179586476925286766559;
  return {0.15 + 0.2*Kokkos::sin(twopi*x1)*Kokkos::cos(twopi*x2),
          -0.25 + 0.15*Kokkos::cos(twopi*x2)*Kokkos::sin(twopi*x3),
          0.35 + 0.1*Kokkos::sin(twopi*x3)*Kokkos::cos(twopi*x1)};
}

KOKKOS_INLINE_FUNCTION
ParticleVector3 NonlinearB(const Real x1, const Real x2, const Real x3) {
  const Real twopi = 6.283185307179586476925286766559;
  return {-0.4 + 0.12*Kokkos::cos(twopi*x1)*Kokkos::sin(twopi*x3),
          0.6 + 0.18*Kokkos::sin(twopi*x2)*Kokkos::cos(twopi*x1),
          1.1 + 0.14*Kokkos::cos(twopi*x3)*Kokkos::sin(twopi*x2)};
}

KOKKOS_INLINE_FUNCTION
ParticleVector3 EvaluateUProfile(TestProfile profile, Real x1, Real x2, Real x3) {
  if (profile == TestProfile::affine) {return AffineU(x1, x2, x3);}
  if (profile == TestProfile::nonlinear) {return NonlinearU(x1, x2, x3);}
  return UniformU();
}

KOKKOS_INLINE_FUNCTION
ParticleVector3 EvaluateBProfile(TestProfile profile, Real x1, Real x2, Real x3) {
  if (profile == TestProfile::affine) {return AffineB(x1, x2, x3);}
  if (profile == TestProfile::nonlinear) {return NonlinearB(x1, x2, x3);}
  return UniformB();
}

KOKKOS_INLINE_FUNCTION
Real WeightSum(const CellCenteredTrilinearStencil &stencil) {
  Real total = 0.0;
  for (int k=0; k<stencil.nk; ++k) {
    for (int j=0; j<stencil.nj; ++j) {
      for (int i=0; i<2; ++i) {
        total += stencil.wei1[i]*stencil.wei2[j]*stencil.wei3[k];
      }
    }
  }
  return total;
}

template<typename PrimitiveField, typename CellField, typename BlockSize>
KOKKOS_INLINE_FUNCTION
GatherResult EvaluateGather(const PrimitiveField &w0, const CellField &bcc,
                            const BlockSize &mbsize, const GatherProbe &probe,
                            int is, int js, int ks) {
  GatherResult result;
  result.stencil = ConstructCellCenteredTrilinearStencil(
      mbsize, probe.m, probe.x1, probe.x2, probe.x3, is, js, ks, true, true);
  result.fields = GatherNewtonianCellCenteredFields(
      w0, bcc, mbsize, probe.m, probe.x1, probe.x2, probe.x3, is, js, ks, true, true);
  result.cE = ConstructIdealCE(result.fields.u, result.fields.b);
  result.cE_dot_b = result.cE.x*result.fields.b.x + result.cE.y*result.fields.b.y +
                    result.cE.z*result.fields.b.z;
  result.weight_sum = WeightSum(result.stencil);
  return result;
}

void FillFields(DvceArray5D<Real> &w0, DvceArray5D<Real> &bcc, const RegionIndcs &indcs,
                const RegionSize &size, TestProfile u_profile, TestProfile b_profile) {
  int n1 = indcs.nx1 + 2*indcs.ng;
  int n2 = indcs.nx2 + 2*indcs.ng;
  int n3 = indcs.nx3 + 2*indcs.ng;
  par_for("cr_rel_gather_fill", DevExeSpace(), 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int k, const int j, const int i) {
    Real x1 = size.x1min + (static_cast<Real>(i - indcs.is) + 0.5)*size.dx1;
    Real x2 = size.x2min + (static_cast<Real>(j - indcs.js) + 0.5)*size.dx2;
    Real x3 = size.x3min + (static_cast<Real>(k - indcs.ks) + 0.5)*size.dx3;
    ParticleVector3 u = EvaluateUProfile(u_profile, x1, x2, x3);
    ParticleVector3 b = EvaluateBProfile(b_profile, x1, x2, x3);
    w0(0,IVX,k,j,i) = u.x;
    w0(0,IVY,k,j,i) = u.y;
    w0(0,IVZ,k,j,i) = u.z;
    bcc(0,IBX,k,j,i) = b.x;
    bcc(0,IBY,k,j,i) = b.y;
    bcc(0,IBZ,k,j,i) = b.z;
  });
  Kokkos::fence();
}

KOKKOS_INLINE_FUNCTION
ParticleVector3 Difference(const ParticleVector3 &left, const ParticleVector3 &right) {
  return {left.x - right.x, left.y - right.y, left.z - right.z};
}

DualArray1D<RegionSize> MakeBlockSize(const std::string &label, const RegionSize &size) {
  DualArray1D<RegionSize> mbsize(label, 1);
  mbsize.h_view(0) = size;
  mbsize.modify<HostMemSpace>();
  mbsize.sync<DevMemSpace>();
  return mbsize;
}

HostBlockSize MakeHostBlockSize(const std::string &label, const RegionSize &size) {
  HostBlockSize mbsize{HostArray1D<RegionSize>(label, 1)};
  mbsize.d_view(0) = size;
  return mbsize;
}

void RequireStencilBounds(const CellCenteredTrilinearStencil &stencil,
                          const RegionIndcs &indcs, const std::string &label) {
  int n1 = indcs.nx1 + 2*indcs.ng;
  int n2 = indcs.nx2 + 2*indcs.ng;
  int n3 = indcs.nx3 + 2*indcs.ng;
  Require(stencil.i0 >= 0 && stencil.i0 + 1 < n1, label + " x1 stencil bounds");
  Require(stencil.j0 >= 0 && stencil.j0 + stencil.nj - 1 < n2,
          label + " x2 stencil bounds");
  Require(stencil.k0 >= 0 && stencil.k0 + stencil.nk - 1 < n3,
          label + " x3 stencil bounds");
}

template<typename BlockSize>
GatherResult EvaluateDeviceGather(const DvceArray5D<Real> &w0, const DvceArray5D<Real> &bcc,
                                  const BlockSize &mbsize, const GatherProbe &probe,
                                  const RegionIndcs &indcs) {
  DvceArray1D<GatherResult> device_result("cr_rel_gather_device_result", 1);
  par_for("cr_rel_gather_device_probe", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(const int) {
    device_result(0) = EvaluateGather(w0, bcc, mbsize, probe,
                                      indcs.is, indcs.js, indcs.ks);
  });
  auto host_result = Kokkos::create_mirror_view_and_copy(HostMemSpace(), device_result);
  return host_result(0);
}

void CheckAnalyticalCase(DvceArray5D<Real> &w0, DvceArray5D<Real> &bcc,
                         const DualArray1D<RegionSize> &mbsize, const RegionIndcs &indcs,
                         const RegionSize &size, const std::vector<GatherProbe> &probes,
                         TestProfile u_profile, TestProfile b_profile, Real tolerance,
                         const std::string &label) {
  FillFields(w0, bcc, indcs, size, u_profile, b_profile);
  int nprobe = static_cast<int>(probes.size());
  DvceArray1D<GatherProbe> device_probes("cr_rel_gather_probes", nprobe);
  DvceArray1D<GatherResult> device_results("cr_rel_gather_results", nprobe);
  auto host_probes = Kokkos::create_mirror_view(device_probes);
  for (int n=0; n<nprobe; ++n) {host_probes(n) = probes[n];}
  Kokkos::deep_copy(device_probes, host_probes);
  auto mbsize_ = mbsize;
  par_for("cr_rel_gather_probe", DevExeSpace(), 0, nprobe-1,
  KOKKOS_LAMBDA(const int n) {
    device_results(n) = EvaluateGather(w0, bcc, mbsize_, device_probes(n),
                                       indcs.is, indcs.js, indcs.ks);
  });
  auto host_results = Kokkos::create_mirror_view_and_copy(HostMemSpace(), device_results);
  auto host_w0 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), w0);
  auto host_bcc = Kokkos::create_mirror_view_and_copy(HostMemSpace(), bcc);
  HostBlockSize host_mbsize = MakeHostBlockSize("cr_rel_gather_host_size", size);

  for (int n=0; n<nprobe; ++n) {
    GatherResult host = EvaluateGather(host_w0, host_bcc, host_mbsize, probes[n],
                                       indcs.is, indcs.js, indcs.ks);
    ParticleVector3 expected_u = EvaluateUProfile(
        u_profile, probes[n].x1, probes[n].x2, probes[n].x3);
    ParticleVector3 expected_b = EvaluateBProfile(
        b_profile, probes[n].x1, probes[n].x2, probes[n].x3);
    ParticleVector3 expected_cE = ConstructIdealCE(expected_u, expected_b);
    RequireVectorNear(host_results(n).fields.u, host.fields.u, tolerance,
                      label + " host/device u parity");
    RequireVectorNear(host_results(n).fields.b, host.fields.b, tolerance,
                      label + " host/device B parity");
    RequireVectorNear(host_results(n).cE, host.cE, tolerance,
                      label + " host/device cE parity");
    RequireVectorNear(host.fields.u, expected_u, tolerance, label + " gathered u oracle");
    RequireVectorNear(host.fields.b, expected_b, tolerance, label + " gathered B oracle");
    RequireVectorNear(host.cE, expected_cE, tolerance, label + " cE oracle");
    Require(NearlyEqual(host.weight_sum, 1.0, tolerance), label + " weight sum");
    Real invariant_scale = std::max(static_cast<Real>(1.0),
                                    VectorNorm(host.cE)*VectorNorm(host.fields.b));
    Require(std::fabs(host.cE_dot_b) <= tolerance*invariant_scale,
            label + " cE dot B invariant");
    RequireStencilBounds(host_results(n).stencil, indcs, label);
  }
}

Real ComputeNonlinearRmsError(const int ncell, const bool vary_u) {
  RegionIndcs indcs{2, ncell, ncell, ncell, 2, ncell + 1, 2, ncell + 1,
                    2, ncell + 1, 0, 0, 0, 0, 0, 0};
  Real dx = 1.0/static_cast<Real>(ncell);
  RegionSize size{0.0, 0.0, 0.0, 1.0, 1.0, 1.0, dx, dx, dx};
  int ntotal = ncell + 2*indcs.ng;
  DvceArray5D<Real> w0("cr_rel_gather_convergence_w0", 1, 5, ntotal, ntotal, ntotal);
  DvceArray5D<Real> bcc("cr_rel_gather_convergence_bcc", 1, 3, ntotal, ntotal, ntotal);
  DualArray1D<RegionSize> mbsize = MakeBlockSize("cr_rel_gather_convergence_size", size);
  FillFields(w0, bcc, indcs, size,
             vary_u ? TestProfile::nonlinear : TestProfile::uniform,
             vary_u ? TestProfile::uniform : TestProfile::nonlinear);

  int nprobe = ncell*ncell*ncell;
  Real error_sum = 0.0;
  auto mbsize_ = mbsize;
  Kokkos::parallel_reduce("cr_rel_gather_convergence",
    Kokkos::RangePolicy<>(DevExeSpace(), 0, nprobe),
  KOKKOS_LAMBDA(const int idx, Real &sum) {
    int i = idx % ncell;
    int j = (idx/ncell) % ncell;
    int k = idx/(ncell*ncell);
    GatherProbe probe{0, (static_cast<Real>(i) + 0.37)*dx,
                      (static_cast<Real>(j) + 0.41)*dx,
                      (static_cast<Real>(k) + 0.58)*dx};
    GatherResult actual = EvaluateGather(w0, bcc, mbsize_, probe,
                                         indcs.is, indcs.js, indcs.ks);
    ParticleVector3 expected_u = EvaluateUProfile(
        vary_u ? TestProfile::nonlinear : TestProfile::uniform,
        probe.x1, probe.x2, probe.x3);
    ParticleVector3 expected_b = EvaluateBProfile(
        vary_u ? TestProfile::uniform : TestProfile::nonlinear,
        probe.x1, probe.x2, probe.x3);
    ParticleVector3 expected_cE = ConstructIdealCE(expected_u, expected_b);
    ParticleVector3 du = Difference(actual.fields.u, expected_u);
    ParticleVector3 db = Difference(actual.fields.b, expected_b);
    ParticleVector3 dcE = Difference(actual.cE, expected_cE);
    sum += du.x*du.x + du.y*du.y + du.z*du.z +
           db.x*db.x + db.y*db.y + db.z*db.z +
           dcE.x*dcE.x + dcE.y*dcE.y + dcE.z*dcE.z;
  }, error_sum);
  return std::sqrt(error_sum/(9.0*static_cast<Real>(nprobe)));
}

void CheckSyntheticInterfaces(const Real tolerance) {
  RegionIndcs indcs{2, 4, 4, 4, 2, 5, 2, 5, 2, 5, 0, 0, 0, 0, 0, 0};
  RegionSize left_size{0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 0.125, 0.25, 0.25};
  RegionSize right_size{0.5, 0.0, 0.0, 1.0, 1.0, 1.0, 0.125, 0.25, 0.25};
  DvceArray5D<Real> left_w0("cr_rel_gather_left_w0", 1, 5, 8, 8, 8);
  DvceArray5D<Real> left_bcc("cr_rel_gather_left_bcc", 1, 3, 8, 8, 8);
  DvceArray5D<Real> right_w0("cr_rel_gather_right_w0", 1, 5, 8, 8, 8);
  DvceArray5D<Real> right_bcc("cr_rel_gather_right_bcc", 1, 3, 8, 8, 8);
  DualArray1D<RegionSize> left_mbsize = MakeBlockSize("cr_rel_gather_left_size", left_size);
  DualArray1D<RegionSize> right_mbsize = MakeBlockSize("cr_rel_gather_right_size", right_size);
  FillFields(left_w0, left_bcc, indcs, left_size, TestProfile::affine,
             TestProfile::affine);
  FillFields(right_w0, right_bcc, indcs, right_size, TestProfile::affine,
             TestProfile::affine);
  GatherProbe interface_probe{0, 0.5, 0.47, 0.61};
  GatherResult left = EvaluateDeviceGather(
      left_w0, left_bcc, left_mbsize, interface_probe, indcs);
  GatherResult right = EvaluateDeviceGather(
      right_w0, right_bcc, right_mbsize, interface_probe, indcs);
  RequireVectorNear(left.fields.u, right.fields.u, tolerance,
                    "synthetic block interface u agreement");
  RequireVectorNear(left.fields.b, right.fields.b, tolerance,
                    "synthetic block interface B agreement");
  RequireVectorNear(left.cE, right.cE, tolerance,
                    "synthetic block interface cE agreement");

  RegionSize coarse_size{0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.25, 0.25, 0.25};
  RegionIndcs fine_indcs{2, 8, 8, 8, 2, 9, 2, 9, 2, 9, 0, 0, 0, 0, 0, 0};
  RegionSize fine_size{0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.125, 0.125, 0.125};
  DvceArray5D<Real> coarse_w0("cr_rel_gather_coarse_w0", 1, 5, 8, 8, 8);
  DvceArray5D<Real> coarse_bcc("cr_rel_gather_coarse_bcc", 1, 3, 8, 8, 8);
  DvceArray5D<Real> fine_w0("cr_rel_gather_fine_w0", 1, 5, 12, 12, 12);
  DvceArray5D<Real> fine_bcc("cr_rel_gather_fine_bcc", 1, 3, 12, 12, 12);
  DualArray1D<RegionSize> coarse_mbsize = MakeBlockSize("cr_rel_gather_coarse_size",
                                                        coarse_size);
  DualArray1D<RegionSize> fine_mbsize = MakeBlockSize("cr_rel_gather_fine_size", fine_size);
  FillFields(coarse_w0, coarse_bcc, indcs, coarse_size, TestProfile::affine,
             TestProfile::affine);
  FillFields(fine_w0, fine_bcc, fine_indcs, fine_size, TestProfile::affine,
             TestProfile::affine);
  GatherProbe refinement_probe{0, 0.37, 0.41, 0.58};
  GatherResult coarse = EvaluateDeviceGather(
      coarse_w0, coarse_bcc, coarse_mbsize, refinement_probe, indcs);
  GatherResult fine = EvaluateDeviceGather(
      fine_w0, fine_bcc, fine_mbsize, refinement_probe, fine_indcs);
  RequireVectorNear(coarse.fields.u, fine.fields.u, tolerance,
                    "synthetic coarse/fine affine u agreement");
  RequireVectorNear(coarse.fields.b, fine.fields.b, tolerance,
                    "synthetic coarse/fine affine B agreement");
  RequireVectorNear(coarse.cE, fine.cE, tolerance,
                    "synthetic coarse/fine affine cE agreement");
}

void RunGatherHelperTests() {
  const Real epsilon = std::numeric_limits<Real>::epsilon();
  const Real tolerance = 256.0*epsilon;
  RegionIndcs indcs{2, 4, 4, 4, 2, 5, 2, 5, 2, 5, 0, 0, 0, 0, 0, 0};
  RegionSize size{0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.25, 0.25, 0.25};
  DvceArray5D<Real> w0("cr_rel_gather_w0",1,5,8,8,8);
  DvceArray5D<Real> bcc("cr_rel_gather_bcc",1,3,8,8,8);
  DualArray1D<RegionSize> mbsize("cr_rel_gather_mbsize",1);
  mbsize.h_view(0) = size;
  mbsize.modify<HostMemSpace>();
  mbsize.sync<DevMemSpace>();

  std::vector<GatherProbe> probes = {
    {0, 0.37, 0.41, 0.58},
    {0, 0.0, 0.0, 0.0},
    {0, std::nextafter(0.0, -1.0), std::nextafter(1.0, 2.0), 0.5},
    {0, 1.0, 1.0, 1.0},
  };
  CheckAnalyticalCase(w0, bcc, mbsize, indcs, size, probes,
                      TestProfile::uniform, TestProfile::uniform,
                      tolerance, "uniform");
  CheckAnalyticalCase(w0, bcc, mbsize, indcs, size, probes, TestProfile::affine,
                      TestProfile::uniform,
                      512.0*epsilon, "affine-u");
  CheckAnalyticalCase(w0, bcc, mbsize, indcs, size, probes,
                      TestProfile::uniform, TestProfile::affine,
                      512.0*epsilon, "affine-B");

  ParticleVector3 zero = ConstructIdealCE({2.0, -4.0, 6.0}, {-1.0, 2.0, -3.0});
  RequireVectorNear(zero, {0.0, 0.0, 0.0}, tolerance, "parallel-field cancellation");
  RequireVectorNear(ConstructIdealCE({1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}),
                    {0.0, 0.0, -1.0}, tolerance, "x cross y sign");
  RequireVectorNear(ConstructIdealCE({0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}),
                    {-1.0, 0.0, 0.0}, tolerance, "y cross z sign");
  RequireVectorNear(ConstructIdealCE({0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}),
                    {0.0, -1.0, 0.0}, tolerance, "z cross x sign");

  FillFields(w0, bcc, indcs, size, TestProfile::affine, TestProfile::affine);
  GatherProbe trap{0, 0.39, 0.47, 0.61};
  GatherResult gathered = EvaluateDeviceGather(w0, bcc, mbsize, trap, indcs);
  ParticleVector3 pointwise = ConstructIdealCE(AffineU(trap.x1, trap.x2, trap.x3),
                                               AffineB(trap.x1, trap.x2, trap.x3));
  RequireVectorNear(gathered.cE, pointwise, 512.0*epsilon,
                    "cross gathered affine fields");
  DvceArray5D<Real> cell_cE("cr_rel_gather_cell_cE", 1, 3, 8, 8, 8);
  par_for("cr_rel_gather_fill_cell_cE", DevExeSpace(), 0, 7, 0, 7, 0, 7,
  KOKKOS_LAMBDA(const int k, const int j, const int i) {
    ParticleVector3 u{w0(0,IVX,k,j,i), w0(0,IVY,k,j,i), w0(0,IVZ,k,j,i)};
    ParticleVector3 b{bcc(0,IBX,k,j,i), bcc(0,IBY,k,j,i), bcc(0,IBZ,k,j,i)};
    ParticleVector3 cE = ConstructIdealCE(u, b);
    cell_cE(0,0,k,j,i) = cE.x;
    cell_cE(0,1,k,j,i) = cE.y;
    cell_cE(0,2,k,j,i) = cE.z;
  });
  Kokkos::fence();
  DvceArray1D<ParticleVector3> device_gathered_cell_cE(
      "cr_rel_gather_device_cell_cE", 1);
  CellCenteredTrilinearStencil trap_stencil = gathered.stencil;
  par_for("cr_rel_gather_cell_cE_probe", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(const int) {
    device_gathered_cell_cE(0) = trap_stencil.Gather(cell_cE, 0, 0, 1, 2);
  });
  auto host_gathered_cell_cE =
      Kokkos::create_mirror_view_and_copy(HostMemSpace(), device_gathered_cell_cE);
  ParticleVector3 gathered_cell_cE = host_gathered_cell_cE(0);
  Require(VectorNorm(Difference(gathered.cE, gathered_cell_cE)) > 1024.0*epsilon,
          "cross-after-gather differs from gather-after-cross");

  for (bool vary_u : {true, false}) {
    Real error8 = ComputeNonlinearRmsError(8, vary_u);
    Real error16 = ComputeNonlinearRmsError(16, vary_u);
    Real error32 = ComputeNonlinearRmsError(32, vary_u);
    Real order8to16 = std::log(error8/error16)/std::log(2.0);
    Real order16to32 = std::log(error16/error32)/std::log(2.0);
    std::string label = vary_u ? "varying-u" : "varying-B";
    std::cout << "CR relativistic gather " << label << " nonlinear RMS errors: "
              << error8 << " " << error16 << " " << error32 << "; orders: "
              << order8to16 << " " << order16to32 << std::endl;
    Require(order8to16 >= 1.8, label + " nonlinear gather 8-to-16 convergence order");
    Require(order16to32 >= 1.8, label + " nonlinear gather 16-to-32 convergence order");
  }
  CheckSyntheticInterfaces(1024.0*epsilon);

  FillFields(w0, bcc, indcs, size, TestProfile::affine, TestProfile::affine);
  GatherResult interior_before = EvaluateDeviceGather(w0, bcc, mbsize, probes[0], indcs);
  GatherResult boundary_before = EvaluateDeviceGather(w0, bcc, mbsize, probes[1], indcs);
  par_for("cr_rel_gather_poison", DevExeSpace(), 0, 7, 0, 7, 0, 7,
  KOKKOS_LAMBDA(const int k, const int j, const int i) {
    if (i < indcs.is || i > indcs.ie || j < indcs.js || j > indcs.je ||
        k < indcs.ks || k > indcs.ke) {
      w0(0,IVX,k,j,i) = 91.0;
      w0(0,IVY,k,j,i) = -73.0;
      w0(0,IVZ,k,j,i) = 57.0;
      bcc(0,IBX,k,j,i) = -41.0;
      bcc(0,IBY,k,j,i) = 37.0;
      bcc(0,IBZ,k,j,i) = -29.0;
    }
  });
  GatherResult interior_after = EvaluateDeviceGather(w0, bcc, mbsize, probes[0], indcs);
  GatherResult boundary_after = EvaluateDeviceGather(w0, bcc, mbsize, probes[1], indcs);
  RequireVectorNear(interior_after.fields.u, interior_before.fields.u, tolerance,
                    "interior gather ignores poisoned ghosts");
  Require(VectorNorm({boundary_after.fields.u.x - boundary_before.fields.u.x,
                      boundary_after.fields.u.y - boundary_before.fields.u.y,
                      boundary_after.fields.u.z - boundary_before.fields.u.z}) > 1.0,
          "boundary gather consumes ghost cells");
}

} // namespace

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  (void)pin;
  (void)restart;
  RunGatherHelperTests();
  std::cout << "CR relativistic gather helper tests passed" << std::endl;
}
