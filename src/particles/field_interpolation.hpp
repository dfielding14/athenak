#ifndef PARTICLES_FIELD_INTERPOLATION_HPP_
#define PARTICLES_FIELD_INTERPOLATION_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file field_interpolation.hpp
//! \brief Shared particle-field interpolation helpers.

#include "athena.hpp"

namespace particles {

template <typename SizeView>
KOKKOS_INLINE_FUNCTION
void InterpolateTSCFields(const RegionIndcs indcs, const SizeView size,
                          const DvceArray5D<Real> bcc,
                          const DvceArray5D<Real> w0,
                          const bool use_mhd_fluid_velocity, const int m,
                          Real x, Real y, Real z, Real &Bx, Real &By,
                          Real &Bz, Real &Ux, Real &Uy, Real &Uz,
                          const bool allow_2d3v) {
  const bool three_d = (indcs.nx3 > 1);
  const bool use_bz_channel = three_d || allow_2d3v;
  const Real dx1 = size.d_view(m).dx1;
  const Real dx2 = size.d_view(m).dx2;
  const Real dx3 = three_d ? size.d_view(m).dx3 : 1.0;

  const Real fx = (x - size.d_view(m).x1min) / dx1 - 0.5;
  const Real fy = (y - size.d_view(m).x2min) / dx2 - 0.5;
  const Real fz = three_d ? (z - size.d_view(m).x3min) / dx3 - 0.5 : 0.0;

  // Anchor the three-point TSC stencil on the nearest cell center. Using
  // floor(f) drops one support point when a particle moves just below a center.
  const Real fic = floor(fx + 0.5);
  const Real fjc = floor(fy + 0.5);
  const Real fkc = three_d ? floor(fz + 0.5) : 0.0;
  const int ic = static_cast<int>(fic) + indcs.is;
  const int jc = static_cast<int>(fjc) + indcs.js;
  const int kc = three_d ? static_cast<int>(fkc) + indcs.ks : indcs.ks;

  const Real di = fx - fic;
  const Real dj = fy - fjc;
  const Real dk = three_d ? (fz - fkc) : 0.0;

  auto weight = [](const Real d) {
    const Real ad = fabs(d);
    if (ad < 0.5) return 0.75 - ad * ad;
    if (ad < 1.5) {
      const Real t = 1.5 - ad;
      return 0.5 * t * t;
    }
    return 0.0;
  };

  Real wx[3] = {weight(di + 1.0), weight(di), weight(di - 1.0)};
  Real wy[3] = {weight(dj + 1.0), weight(dj), weight(dj - 1.0)};
  Real wz[3] = {1.0, 0.0, 0.0};
  int kz[3] = {kc, kc, kc};
  if (three_d) {
    wz[0] = weight(dk + 1.0);
    wz[1] = weight(dk);
    wz[2] = weight(dk - 1.0);
    kz[0] = kc - 1;
    kz[1] = kc;
    kz[2] = kc + 1;
  }

  int ix[3] = {ic - 1, ic, ic + 1};
  int iy[3] = {jc - 1, jc, jc + 1};
  // Active MHD fields carry synchronized same-level and refinement-interface
  // ghost zones. The no-MHD manufactured carrier has no boundary exchange.
  const int i_min = use_mhd_fluid_velocity ? indcs.is - indcs.ng : indcs.is;
  const int i_max = use_mhd_fluid_velocity ? indcs.ie + indcs.ng : indcs.ie;
  const int j_min = use_mhd_fluid_velocity ? indcs.js - indcs.ng : indcs.js;
  const int j_max = use_mhd_fluid_velocity ? indcs.je + indcs.ng : indcs.je;
  const int k_min = use_mhd_fluid_velocity ? indcs.ks - indcs.ng : indcs.ks;
  const int k_max = use_mhd_fluid_velocity ? indcs.ke + indcs.ng : indcs.ke;
  for (int n = 0; n < 3; ++n) {
    ix[n] = (ix[n] < i_min) ? i_min : ((ix[n] > i_max) ? i_max : ix[n]);
    iy[n] = (iy[n] < j_min) ? j_min : ((iy[n] > j_max) ? j_max : iy[n]);
    if (three_d) {
      kz[n] = (kz[n] < k_min) ? k_min : ((kz[n] > k_max) ? k_max : kz[n]);
    }
  }

  Bx = By = Bz = 0.0;
  Ux = Uy = Uz = 0.0;
  const int nk = three_d ? 3 : 1;
  for (int kk = 0; kk < nk; ++kk) {
    for (int jj = 0; jj < 3; ++jj) {
      for (int ii = 0; ii < 3; ++ii) {
        const Real w = wx[ii] * wy[jj] * wz[kk];
        Bx += w * bcc(m, IBX, kz[kk], iy[jj], ix[ii]);
        By += w * bcc(m, IBY, kz[kk], iy[jj], ix[ii]);
        if (use_mhd_fluid_velocity) {
          Ux += w*w0(m, IVX, kz[kk], iy[jj], ix[ii]);
          Uy += w*w0(m, IVY, kz[kk], iy[jj], ix[ii]);
          Uz += w*w0(m, IVZ, kz[kk], iy[jj], ix[ii]);
        }
        if (use_bz_channel) {
          Bz += w * bcc(m, IBZ, kz[kk], iy[jj], ix[ii]);
        }
      }
    }
  }
  if (!use_bz_channel) {
    Bz = 0.0;
    Uz = 0.0;
  }
}

}  // namespace particles

#endif  // PARTICLES_FIELD_INTERPOLATION_HPP_
