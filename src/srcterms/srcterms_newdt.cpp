//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file srcterms_newdt.cpp
//! \brief function to compute timestep for source terms across all MeshBlock(s) in a
//! MeshBlockPack

#include <float.h>

#include <limits>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "cooling.hpp"
#include "srcterms.hpp"

//----------------------------------------------------------------------------------------
//! \fn void SourceTerms::NewTimeStep()
//! \brief Compute new timestep for source terms.

void SourceTerms::NewTimeStep(const DvceArray5D<Real> &w0, const EOS_Data &eos_data) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, nx1 = indcs.nx1;
  int js = indcs.js, nx2 = indcs.nx2;
  int ks = indcs.ks, nx3 = indcs.nx3;
  const int nmkji = (pmy_pack->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  dtnew = static_cast<Real>(std::numeric_limits<float>::max());

  if (pcooling != nullptr) pcooling->NewTimeStep(w0, eos_data, dtnew);

  if (rel_cooling) {
    Real use_e = eos_data.use_e;
    Real gamma = eos_data.gamma;
    Real gm1 = gamma - 1.0;
    Real cooling_rate = crate_rel;
    Real cooling_power = cpower_rel;

    // find smallest (e/cooling_rate) in each cell
    Kokkos::parallel_reduce("srcterms_cooling_newdt",
                            Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
    KOKKOS_LAMBDA(const int &idx, Real &min_dt) {
      // compute m,k,j,i indices of thread and call function
      int m = (idx)/nkji;
      int k = (idx - m*nkji)/nji;
      int j = (idx - m*nkji - k*nji)/nx1;
      int i = (idx - m*nkji - k*nji - j*nx1) + is;
      k += ks;
      j += js;

      // temperature in cgs unit
      Real temp = 1.0;
      Real eint = 1.0;
      if (use_e) {
        temp = w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1;
        eint = w0(m,IEN,k,j,i);
      } else {
        temp = w0(m,ITM,k,j,i);
        eint = w0(m,ITM,k,j,i)*w0(m,IDN,k,j,i)/gm1;
      }

      auto &ux = w0(m, IVX, k, j, i);
      auto &uy = w0(m, IVY, k, j, i);
      auto &uz = w0(m, IVZ, k, j, i);

      auto ut = 1. + ux * ux + uy * uy + uz * uz;
      ut = sqrt(ut);

      // The following should be approximately correct
      // add a tiny number
      Real cooling_heating = FLT_MIN + fabs(w0(m,IDN,k,j,i) * ut *
                             pow(temp*cooling_rate, cooling_power));

      min_dt = fmin((eint/cooling_heating), min_dt);
    }, Kokkos::Min<Real>(dtnew));
  }

  return;
}
