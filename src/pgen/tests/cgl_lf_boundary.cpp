//========================================================================================
// AthenaK CGL-LF fixed-state and PrimToCons user-boundary regression.
//========================================================================================
#include <cmath>
#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "eos/ideal_c2p_mhd.hpp"
#include "mhd/mhd.hpp"
#include "pgen/pgen.hpp"
#include "coordinates/cell_locations.hpp"

namespace {
MHDPrim1D boundary_state;
Real boundary_scalar;

// The callback fills ghosts without changing active cells or advancing any state.
// PrimToCons produces the current A/mu representation, including in a temporary
// buffer. Copying only owned physical ghosts avoids changing internal MB ghosts.
void FixedPrimitiveBoundary(Mesh *pm) {
  auto *pp = pm->pmb_pack;
  auto *mhd = pp->pmhd;
  const auto s = boundary_state;
  const Real scalar = boundary_scalar;
  const auto ind = pm->mb_indcs;
  const int nmb = pp->nmb_thispack;
  const int n1 = mhd->u0.extent_int(4);
  const int n2 = mhd->u0.extent_int(3);
  const int n3 = mhd->u0.extent_int(2);
  const int nv = mhd->u0.extent_int(1);
  const int nmhd = mhd->nmhd;
  const bool two = pm->multi_d, three = pm->three_d;
  auto flags = pp->pmb->mb_bcs.d_view;
  auto b1 = mhd->b0.x1f, b2 = mhd->b0.x2f, b3 = mhd->b0.x3f;
  const auto owned = KOKKOS_LAMBDA(const int m, const int k, const int j,
                                  const int i, const int component) {
    return (i < ind.is && flags(m, inner_x1) == BoundaryFlag::user) ||
        (i > ind.ie + (component == 0) &&
         flags(m, outer_x1) == BoundaryFlag::user) ||
        (two && j < ind.js && flags(m, inner_x2) == BoundaryFlag::user) ||
        (two && j > ind.je + (component == 1) &&
         flags(m, outer_x2) == BoundaryFlag::user) ||
        (three && k < ind.ks && flags(m, inner_x3) == BoundaryFlag::user) ||
        (three && k > ind.ke + (component == 2) &&
         flags(m, outer_x3) == BoundaryFlag::user);
  };
  par_for("boundary_test_b1", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2-1, 0, n1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    if (owned(m, k, j, i, 0)) b1(m, k, j, i) = s.bx;
  });
  par_for("boundary_test_b2", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2, 0, n1-1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    if (owned(m, k, j, i, 1)) b2(m, k, j, i) = s.by;
  });
  par_for("boundary_test_b3", DevExeSpace(), 0, nmb-1, 0, n3, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    if (owned(m, k, j, i, 2)) b3(m, k, j, i) = s.bz;
  });
  DvceArray5D<Real> prim("boundary_prim", nmb, nv, n3, n2, n1);
  DvceArray5D<Real> cons("boundary_cons", nmb, nv, n3, n2, n1);
  DvceArray5D<Real> bcc("boundary_bcc", nmb, 3, n3, n2, n1);
  par_for("boundary_test_prim", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    prim(m, IDN, k, j, i) = s.d;
    prim(m, IVX, k, j, i) = s.vx;
    prim(m, IVY, k, j, i) = s.vy;
    prim(m, IVZ, k, j, i) = s.vz;
    prim(m, IPR, k, j, i) = s.e;
    prim(m, IPP, k, j, i) = s.pp;
    for (int n=nmhd; n < nv; ++n) prim(m, n, k, j, i) = scalar;
    bcc(m, IBX, k, j, i) = 0.5*(b1(m, k, j, i) + b1(m, k, j, i+1));
    bcc(m, IBY, k, j, i) = 0.5*(b2(m, k, j, i) + b2(m, k, j+1, i));
    bcc(m, IBZ, k, j, i) = 0.5*(b3(m, k, j, i) + b3(m, k+1, j, i));
  });
  mhd->peos->PrimToCons(prim, bcc, cons, 0, n1-1, 0, n2-1, 0, n3-1);
  auto u = mhd->u0;
  par_for("boundary_test_copy", DevExeSpace(), 0, nmb-1, 0, nv-1,
          0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(int m, int n, int k, int j, int i) {
    if (owned(m, k, j, i, -1)) u(m, n, k, j, i) = cons(m, n, k, j, i);
  });
}
}  // namespace

void ProblemGenerator::CGLLFBoundary(ParameterInput *pin, const bool restart) {
  auto *pp = pmy_mesh_->pmb_pack;
  auto *mhd = pp->pmhd;
  boundary_state.d = pin->GetOrAddReal("problem", "rho0", 1.0);
  boundary_state.vx = pin->GetOrAddReal("problem", "vx0", 0.0);
  boundary_state.vy = pin->GetOrAddReal("problem", "vy0", 0.0);
  boundary_state.vz = pin->GetOrAddReal("problem", "vz0", 0.0);
  boundary_state.e = pin->GetOrAddReal("problem", "ppar0", 1.0);
  boundary_state.pp = pin->GetOrAddReal("problem", "pperp0", 1.01);
  boundary_state.bx = pin->GetOrAddReal("problem", "guide_b1", 0.7);
  boundary_state.by = pin->GetOrAddReal("problem", "guide_b2", 0.2);
  boundary_state.bz = pin->GetOrAddReal("problem", "guide_b3", -0.15);
  boundary_scalar = pin->GetOrAddReal("problem", "scalar0", 0.5);
  const Real amp = pin->GetOrAddReal("problem", "amp", 0.0);
  if (user_bcs) user_bcs_func = FixedPrimitiveBoundary;
  const auto s = boundary_state;
  HydCons1D fixed;
  SingleP2C_CGLMHD(s, mhd->peos->eos_data.bfloor, fixed);
  auto &uin = mhd->pbval_u->u_in;
  auto &bin = mhd->pbval_b->b_in;
  for (int f=0; f < 6; ++f) {
    uin.h_view(IDN, f) = fixed.d;
    uin.h_view(IM1, f) = fixed.mx;
    uin.h_view(IM2, f) = fixed.my;
    uin.h_view(IM3, f) = fixed.mz;
    uin.h_view(IEN, f) = fixed.e;
    uin.h_view(IAN, f) = fixed.mu;
    for (int n=mhd->nmhd; n < mhd->nmhd+mhd->nscalars; ++n) {
      uin.h_view(n, f) = s.d*boundary_scalar;
    }
    bin.h_view(IBX, f) = s.bx;
    bin.h_view(IBY, f) = s.by;
    bin.h_view(IBZ, f) = s.bz;
  }
  uin.modify_host(); uin.sync_device();
  bin.modify_host(); bin.sync_device();
  if (restart) return;
  auto w = mhd->w0, bcc = mhd->bcc0;
  const int nmb = pp->nmb_thispack;
  const int n1 = w.extent_int(4), n2 = w.extent_int(3), n3 = w.extent_int(2);
  const int nv = w.extent_int(1), nmhd = mhd->nmhd;
  const auto ind = pmy_mesh_->mb_indcs;
  const auto size = pp->pmb->mb_size.d_view;
  const Real scalar = boundary_scalar;
  par_for("boundary_test_init", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real x = CellCenterX(i-ind.is, ind.nx1, size(m).x1min, size(m).x1max);
    w(m, IDN, k, j, i) = s.d;
    w(m, IVX, k, j, i) = s.vx;
    w(m, IVY, k, j, i) = s.vy;
    w(m, IVZ, k, j, i) = s.vz;
    w(m, IPR, k, j, i) = s.e*(1.0 + amp*sin(2.0*M_PI*x));
    w(m, IPP, k, j, i) = s.pp*(1.0 + amp*sin(2.0*M_PI*x));
    for (int n=nmhd; n < nv; ++n) w(m, n, k, j, i) = scalar;
    bcc(m, IBX, k, j, i) = s.bx;
    bcc(m, IBY, k, j, i) = s.by;
    bcc(m, IBZ, k, j, i) = s.bz;
  });
  Kokkos::deep_copy(mhd->b0.x1f, s.bx);
  Kokkos::deep_copy(mhd->b0.x2f, s.by);
  Kokkos::deep_copy(mhd->b0.x3f, s.bz);
  mhd->peos->PrimToCons(w, bcc, mhd->u0, 0, n1-1, 0, n2-1, 0, n3-1);
}
