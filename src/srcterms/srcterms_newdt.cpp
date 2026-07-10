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
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "coordinates/cell_locations.hpp"
#include "cgm_cooling.hpp"
#include "ismcooling.hpp"
#include "srcterms.hpp"
#include "units/units.hpp"
#include "cooling_tables.hpp"

//----------------------------------------------------------------------------------------
//! \fn void SourceTerms::NewTimeStep()
//! \brief Compute new timestep for source terms.

void SourceTerms::NewTimeStep(const DvceArray5D<Real> &w0, const EOS_Data &eos_data) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, nx1 = indcs.nx1;
  int js = indcs.js, nx2 = indcs.nx2;
  int ks = indcs.ks, nx3 = indcs.nx3;

  // Get active fluid variable count from the source term's physics block
  int nfluid = cgm_fluid_vars;

  const int nmkji = (pmy_pack->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  dtnew = static_cast<Real>(std::numeric_limits<float>::max());

  if (ism_cooling) {
    Real use_e = eos_data.use_e;
    Real gamma = eos_data.gamma;
    Real gm1 = gamma - 1.0;
    Real heating_rate = hrate;
    Real temp_unit = pmy_pack->punit->temperature_cgs();
    Real n_unit = pmy_pack->punit->density_cgs()/pmy_pack->punit->mu()
                  / pmy_pack->punit->atomic_mass_unit_cgs;
    Real cooling_unit = pmy_pack->punit->pressure_cgs()/pmy_pack->punit->time_cgs()
                        / n_unit/n_unit;
    Real heating_unit = pmy_pack->punit->pressure_cgs()/pmy_pack->punit->time_cgs()
                        / n_unit;

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
        temp = temp_unit*w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1;
        eint = w0(m,IEN,k,j,i);
      } else {
        temp = temp_unit*w0(m,ITM,k,j,i);
        eint = w0(m,ITM,k,j,i)*w0(m,IDN,k,j,i)/gm1;
      }

      Real lambda_cooling = ISMCoolFn(temp)/cooling_unit;
      Real gamma_heating = heating_rate/heating_unit;

      // add a tiny number
      Real cooling_heating = FLT_MIN + fabs(w0(m,IDN,k,j,i) *
                             (w0(m,IDN,k,j,i) * lambda_cooling - gamma_heating));

      min_dt = fmin((eint/cooling_heating), min_dt);
    }, Kokkos::Min<Real>(dtnew));
  }

  if (cgm_cooling && cgm_cooling_limit_dt) {
    auto &size = pmy_pack->pmb->mb_size;
    Real gm1 = eos_data.gamma - 1.0;

    auto &units = pmy_pack->punit;
    Real temp_unit = units->temperature_cgs();
    Real nH_unit = units->density_cgs()/units->atomic_mass_unit_cgs;
    Real cooling_unit = units->pressure_cgs()/units->time_cgs()/nH_unit/nH_unit;
    Real heating_unit = units->pressure_cgs()/units->time_cgs()/nH_unit;
    Real length_unit = units->length_cgs();

    auto Tbins_ = Tbins.d_view;
    auto nHbins_ = nHbins.d_view;
    auto Metal_Cooling_ = Metal_Cooling.d_view;
    auto H_He_Cooling_ = H_He_Cooling.d_view;
    auto Metal_Cooling_CIE_ = Metal_Cooling_CIE.d_view;
    auto H_He_Cooling_CIE_ = H_He_Cooling_CIE.d_view;

    auto Tfloor  = std::pow(10, Tbins_ARR[0]);
    auto Tceil   = std::pow(10, Tbins_ARR[Tbins_DIM_0 - 1]);
    auto nHfloor = std::pow(10, nHbins_ARR[0]);
    auto nHceil  = std::pow(10, nHbins_ARR[nHbins_DIM_0 - 1]);

    Real h_rate = hrate;
    Real h_norm = hscale_norm;
    Real h_height = hscale_height;
    Real h_radius = hscale_radius;
    Real T_max_ = T_max;
    Real dt_frac = cgm_cooling_dt_frac;

    // find smallest (e/cooling_rate) in each cell
    Kokkos::parallel_reduce("srcterms_cgm_cooling_newdt",
                            Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
    KOKKOS_LAMBDA(const int &idx, Real &min_dt) {
      // compute m,k,j,i indices of thread and call function
      int m = (idx)/nkji;
      int k = (idx - m*nkji)/nji;
      int j = (idx - m*nkji - k*nji)/nx1;
      int i = (idx - m*nkji - k*nji - j*nx1) + is;
      k += ks;
      j += js;

      const Real x1min = size.d_view(m).x1min, x1max = size.d_view(m).x1max;
      const Real x2min = size.d_view(m).x2min, x2max = size.d_view(m).x2max;
      const Real x3min = size.d_view(m).x3min, x3max = size.d_view(m).x3max;

      const Real x1v = CellCenterX(i - is, nx1, x1min, x1max);
      const Real x2v = CellCenterX(j - js, nx2, x2min, x2max);
      const Real x3v = CellCenterX(k - ks, nx3, x3min, x3max);

      const auto rates = cgm_cooling::CalculateRates(
          w0(m,IDN,k,j,i), w0(m,IEN,k,j,i), w0(m,nfluid,k,j,i), x1v, x2v,
          x3v, size.d_view(m).dx1, gm1, temp_unit, nH_unit, cooling_unit,
          heating_unit, length_unit, h_rate, h_norm, h_height, h_radius,
          T_max_, Tfloor, Tceil, nHfloor, nHceil, Tbins_, nHbins_,
          Metal_Cooling_, H_He_Cooling_, Metal_Cooling_CIE_, H_He_Cooling_CIE_);
      if (rates.net_cooling_rate > 0.0) {
        min_dt = fmin(dt_frac*w0(m,IEN,k,j,i)/
                      (rates.net_cooling_rate + FLT_MIN), min_dt);
      }
    }, Kokkos::Min<Real>(dtnew));

  }

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
