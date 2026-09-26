#ifndef DIFFUSION_CURRENT_LIMITED_RESISTIVITY_HPP_
#define DIFFUSION_CURRENT_LIMITED_RESISTIVITY_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file current_limited_resistivity.hpp
//! \brief Current-dependent Ohmic closure, edge stencils, and reconnecting-field estimate.

#include <cmath>

#include "athena.hpp"
#include "mesh/mesh.hpp"

namespace current_limited {

struct Parameters {
  Real eta0, eta_max, d_i, b_rec, q_star, eta_star;
  DvceArray4D<Real> b_rec_cells;
  int b_rec_offset = 0;  // ng-1: cache has one ghost cell in each active direction
};

struct State {
  Real j1, j2, j3, q, eta;
};

// The cache is empty for the fixed-amplitude model. For jump, it is held fixed
// through the whole cycle, including both parabolic sweeps and the RK stages.
KOKKOS_INLINE_FUNCTION
Real CellBRec(const Parameters &p, const bool multi_d, const bool three_d,
              const int m, const int k, const int j, const int i) {
  if (p.b_rec_cells.extent_int(0) == 0) return p.b_rec;
  return p.b_rec_cells(m, k-(three_d ? p.b_rec_offset : 0),
                       j-(multi_d ? p.b_rec_offset : 0), i-p.b_rec_offset);
}

KOKKOS_INLINE_FUNCTION
Real JumpEstimate(const DvceArray5D<Real> &bcc, const RegionSize &size,
                   const Real radius, const Real floor,
                   const bool multi_d, const bool three_d,
                   const int m, const int k, const int j, const int i) {
  const Real dx[3] = {size.dx1, size.dx2, size.dx3};
  const int ndim = three_d ? 3 : (multi_d ? 2 : 1);
  Real estimate = floor;
  for (int a=0; a<ndim; ++a) {
    const int reach = static_cast<int>(ceil(radius/dx[a]));
    int lo[3] = {i, j, k}, hi[3] = {i, j, k};
    lo[a] -= reach;
    hi[a] += reach;
    Real jump2 = 0.0;
    for (int c=0; c<3; ++c) {
      const Real delta = bcc(m,c,hi[2],hi[1],hi[0])-bcc(m,c,lo[2],lo[1],lo[0]);
      jump2 += delta*delta;
    }
    estimate = fmax(estimate, 0.5*sqrt(jump2));
  }
  return estimate;
}

// Total diffusivity, not an increment to eta0. Its differential diffusivity at
// fixed density and b_rec is at most eta_max. This regularizes, but does not cap J.
KOKKOS_INLINE_FUNCTION
Real Eta(const Real q, const Parameters &p) {
  if (p.eta_max == p.eta0) return p.eta0;
  if (q <= p.q_star && q < 1.0) return p.eta0/(1.0 - q);
  // Equivalent forms avoid cancellation near q_star and inf/inf for very large q.
  if (q <= 1.0) {
    return p.eta_star + (p.eta_max-p.eta_star)*((q-p.q_star)/q);
  }
  return p.eta_max - (p.q_star/q)*(p.eta_max-p.eta_star);
}

// Native CT curl: component 0 is centered in x1 and on x2/x3 faces, etc.
// Inactive dimensions are collapsed, including duplicate transverse face planes.
KOKKOS_INLINE_FUNCTION
Real NativeCurrent(const DvceFaceFld4D<Real> &b, const RegionSize &size,
                   const bool multi_d, const bool three_d, const int component,
                   const int m, int k, int j, const int i) {
  if (!multi_d) j = 0;
  if (!three_d) k = 0;
  if (component == 0) {
    Real value = 0.0;
    if (multi_d) value += (b.x3f(m,k,j,i)-b.x3f(m,k,j-1,i))/size.dx2;
    if (three_d) value -= (b.x2f(m,k,j,i)-b.x2f(m,k-1,j,i))/size.dx3;
    return value;
  }
  if (component == 1) {
    Real value = -(b.x3f(m,k,j,i)-b.x3f(m,k,j,i-1))/size.dx1;
    if (three_d) value += (b.x1f(m,k,j,i)-b.x1f(m,k-1,j,i))/size.dx3;
    return value;
  }
  Real value = (b.x2f(m,k,j,i)-b.x2f(m,k,j,i-1))/size.dx1;
  if (multi_d) value -= (b.x1f(m,k,j,i)-b.x1f(m,k,j-1,i))/size.dx2;
  return value;
}

KOKKOS_INLINE_FUNCTION
State Evaluate(const Real j1, const Real j2, const Real j3, const Real rho,
               const Parameters &p, const Real b_rec) {
  // d_i is specified at rho=1. A sheet with upstream B_up triggers at
  // thickness d_i(rho)*B_up/b_rec, not necessarily at d_i(rho).
  Real q = sqrt(j1*j1+j2*j2+j3*j3)*(p.d_i/b_rec)/sqrt(rho);
  return {j1, j2, j3, q, Eta(q, p)};
}

// Gather the other two staggered currents onto this component's edge. Each
// transverse gather spans four native edges in 3D, reduced along inactive axes.
KOKKOS_INLINE_FUNCTION
State EdgeState(const DvceFaceFld4D<Real> &b, const DvceArray5D<Real> &w,
                const RegionSize &size, const Parameters &p, const Real dfloor,
                const bool multi_d, const bool three_d, const int component,
                const int m, int k, int j, const int i) {
  if (!multi_d) j = 0;
  if (!three_d) k = 0;
  const bool active[3] = {true, multi_d, three_d};
  Real current[3];
  for (int c=0; c<3; ++c) {
    if (c == component) {
      current[c] = NativeCurrent(b, size, multi_d, three_d, c, m, k, j, i);
    } else {
      Real sum = 0.0;
      const int na = active[component] ? 2 : 1;
      const int nc = active[c] ? 2 : 1;
      for (int a=0; a<na; ++a) {
        for (int d=0; d<nc; ++d) {
          int index[3] = {i, j, k};
          index[component] += a;
          index[c] -= d;
          sum += NativeCurrent(b, size, multi_d, three_d, c, m,
                               index[2], index[1], index[0]);
        }
      }
      current[c] = sum/(na*nc);
    }
  }
  const int a = (component+1)%3, c = (component+2)%3;
  const int na = active[a] ? 2 : 1, nc = active[c] ? 2 : 1;
  Real rho = 0.0;
  Real b_rec = 0.0;
  for (int da=0; da<na; ++da) {
    for (int dc=0; dc<nc; ++dc) {
      int index[3] = {i, j, k};
      index[a] -= da;
      index[c] -= dc;
      rho += fmax(w(m,IDN,index[2],index[1],index[0]), dfloor);
      b_rec += CellBRec(p, multi_d, three_d, m, index[2], index[1], index[0]);
    }
  }
  b_rec = p.b_rec_cells.extent_int(0) == 0 ? p.b_rec : b_rec/(na*nc);
  return Evaluate(current[0], current[1], current[2], rho/(na*nc), p, b_rec);
}

// Diagnostic proxy: average each native current component to the cell center,
// then evaluate the closure there. It is not an average of the edge coefficients.
KOKKOS_INLINE_FUNCTION
State CellState(const DvceFaceFld4D<Real> &b, const DvceArray5D<Real> &w,
                const RegionSize &size, const Parameters &p, const Real dfloor,
                const bool multi_d, const bool three_d,
                const int m, const int k, const int j, const int i) {
  const bool active[3] = {true, multi_d, three_d};
  Real current[3];
  for (int c=0; c<3; ++c) {
    const int a = (c+1)%3, d = (c+2)%3;
    const int na = active[a] ? 2 : 1, nd = active[d] ? 2 : 1;
    Real sum = 0.0;
    for (int da=0; da<na; ++da) {
      for (int dd=0; dd<nd; ++dd) {
        int index[3] = {i, j, k};
        index[a] += da;
        index[d] += dd;
        sum += NativeCurrent(b, size, multi_d, three_d, c, m,
                             index[2], index[1], index[0]);
      }
    }
    current[c] = sum/(na*nd);
  }
  return Evaluate(current[0], current[1], current[2], fmax(w(m,IDN,k,j,i),dfloor), p,
                   CellBRec(p, multi_d, three_d, m, k, j, i));
}

} // namespace current_limited
#endif // DIFFUSION_CURRENT_LIMITED_RESISTIVITY_HPP_
