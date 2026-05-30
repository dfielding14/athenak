#ifndef PARTICLES_TRILINEAR_GATHER_HPP_
#define PARTICLES_TRILINEAR_GATHER_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file trilinear_gather.hpp
//  \brief device-capable cell-centered trilinear gathers for particle pushers

#include "athena.hpp"

namespace particles {

struct ParticleVector3 {
  Real x, y, z;
};

KOKKOS_INLINE_FUNCTION
int TrilinearFloorToInt(Real x) {
  int i = static_cast<int>(x);
  if (static_cast<Real>(i) > x) {--i;}
  return i;
}

struct CellCenteredTrilinearStencil {
  int i0, j0, k0;
  int nj, nk;
  Real wei1[2], wei2[2], wei3[2];

  template<typename CellField>
  KOKKOS_INLINE_FUNCTION
  ParticleVector3 Gather(const CellField &field, int m, int ivx, int ivy,
                         int ivz) const {
    ParticleVector3 value{0.0, 0.0, 0.0};
    for (int k=0; k<nk; ++k) {
      for (int j=0; j<nj; ++j) {
        for (int i=0; i<2; ++i) {
          Real w = wei1[i]*wei2[j]*wei3[k];
          value.x += w*field(m,ivx,k0+k,j0+j,i0+i);
          value.y += w*field(m,ivy,k0+k,j0+j,i0+i);
          value.z += w*field(m,ivz,k0+k,j0+j,i0+i);
        }
      }
    }
    return value;
  }
};

template<typename BlockSize>
KOKKOS_INLINE_FUNCTION
CellCenteredTrilinearStencil ConstructCellCenteredTrilinearStencil(
    const BlockSize &mbsize, int m, Real x1, Real x2, Real x3, int is, int js,
    int ks, bool multi_d, bool three_d) {
  CellCenteredTrilinearStencil stencil;

  Real a1 = (x1 - mbsize.d_view(m).x1min)/mbsize.d_view(m).dx1 +
            static_cast<Real>(is) - 0.5;
  stencil.i0 = TrilinearFloorToInt(a1);
  Real d1 = a1 - static_cast<Real>(stencil.i0);
  stencil.wei1[0] = 1.0 - d1;
  stencil.wei1[1] = d1;

  stencil.j0 = js;
  stencil.nj = 1;
  stencil.wei2[0] = 1.0;
  stencil.wei2[1] = 0.0;
  if (multi_d) {
    Real a2 = (x2 - mbsize.d_view(m).x2min)/mbsize.d_view(m).dx2 +
              static_cast<Real>(js) - 0.5;
    stencil.j0 = TrilinearFloorToInt(a2);
    Real d2 = a2 - static_cast<Real>(stencil.j0);
    stencil.nj = 2;
    stencil.wei2[0] = 1.0 - d2;
    stencil.wei2[1] = d2;
  }

  stencil.k0 = ks;
  stencil.nk = 1;
  stencil.wei3[0] = 1.0;
  stencil.wei3[1] = 0.0;
  if (three_d) {
    Real a3 = (x3 - mbsize.d_view(m).x3min)/mbsize.d_view(m).dx3 +
              static_cast<Real>(ks) - 0.5;
    stencil.k0 = TrilinearFloorToInt(a3);
    Real d3 = a3 - static_cast<Real>(stencil.k0);
    stencil.nk = 2;
    stencil.wei3[0] = 1.0 - d3;
    stencil.wei3[1] = d3;
  }

  return stencil;
}

template<typename CellField>
KOKKOS_INLINE_FUNCTION
ParticleVector3 GatherCellCenteredB(const CellCenteredTrilinearStencil &stencil,
                                    const CellField &bcc, int m) {
  return stencil.Gather(bcc, m, IBX, IBY, IBZ);
}

template<typename CellField, typename BlockSize>
KOKKOS_INLINE_FUNCTION
ParticleVector3 GatherCellCenteredB(const CellField &bcc, const BlockSize &mbsize,
                                    int m, Real x1, Real x2, Real x3, int is,
                                    int js, int ks, bool multi_d, bool three_d) {
  auto stencil = ConstructCellCenteredTrilinearStencil(
      mbsize, m, x1, x2, x3, is, js, ks, multi_d, three_d);
  return GatherCellCenteredB(stencil, bcc, m);
}

struct NewtonianCellCenteredFields {
  ParticleVector3 u;
  ParticleVector3 b;
};

template<typename PrimitiveField, typename CellField, typename BlockSize>
KOKKOS_INLINE_FUNCTION
NewtonianCellCenteredFields GatherNewtonianCellCenteredFields(
    const PrimitiveField &w0, const CellField &bcc, const BlockSize &mbsize,
    int m, Real x1, Real x2, Real x3, int is, int js, int ks, bool multi_d,
    bool three_d) {
  auto stencil = ConstructCellCenteredTrilinearStencil(
      mbsize, m, x1, x2, x3, is, js, ks, multi_d, three_d);
  NewtonianCellCenteredFields fields;
  fields.u = stencil.Gather(w0, m, IVX, IVY, IVZ);
  fields.b = GatherCellCenteredB(stencil, bcc, m);
  return fields;
}

KOKKOS_INLINE_FUNCTION
ParticleVector3 ConstructIdealCE(const ParticleVector3 &u, const ParticleVector3 &b) {
  // This is the explicitly c-normalized electric field: cE = -u x B.
  ParticleVector3 cE;
  cE.x = -u.y*b.z + u.z*b.y;
  cE.y = -u.z*b.x + u.x*b.z;
  cE.z = -u.x*b.y + u.y*b.x;
  return cE;
}

} // namespace particles
#endif // PARTICLES_TRILINEAR_GATHER_HPP_
