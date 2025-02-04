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
 
  if (cgm_cooling) {
    Real use_e = eos_data.use_e;
    Real gamma = eos_data.gamma;
    Real gm1 = gamma - 1.0;

    auto &units = pmy_pack->punit;
    Real temp_unit = units->temperature_cgs();
    Real n_unit = units->density_cgs()/units->mu()/units->atomic_mass_unit_cgs;
    Real cooling_unit = units->pressure_cgs()/units->time_cgs()/n_unit/n_unit;

    auto Tbins_ = Tbins.d_view;
    auto nHbins_ = nHbins.d_view;
    auto Metal_Cooling_ = Metal_Cooling.d_view;
    auto H_He_Cooling_ = H_He_Cooling.d_view;
    auto Metal_Cooling_CIE_ = Metal_Cooling_CIE.d_view;
    auto H_He_Cooling_CIE_ = H_He_Cooling_CIE.d_view;

    auto Tfloor  = Tbins_ARR[0];
    auto Tceil   = Tbins_ARR[Tbins_DIM_0 - 1];
    auto nHfloor = nHbins_ARR[0];
    auto nHceil  = nHbins_ARR[nHbins_DIM_0 - 1];

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

      Real temp = 1.0; // temperature in cgs units
      Real eint = 1.0;
      if (use_e) {
        temp = temp_unit*w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1;
        eint = w0(m,IEN,k,j,i);
      } else {
        temp = temp_unit*w0(m,ITM,k,j,i);
        eint = w0(m,ITM,k,j,i)*w0(m,IDN,k,j,i)/gm1;
      }
      Real nH = n_unit*w0(m,IDN,k,j,i); // density in cgs units
      Real Z = 1.0; // metallicity [Zsun]

      // WiersmaCooling at redshift z = 0 taken from Wiersma et al (2009)
      Real lambda_cooling = 0.0; // Ensure we are in range of cooling table
      if (temp > Tfloor && temp < Tceil && nH > nHfloor && nH < nHceil) {
        // Convert input values to log space
        Real log_temp = log10(temp);
        Real log_density = log10(nH);

        // Locate indices in Tbins and nHbins
        int i = 0, j = 0;
        while (i < Tbins_DIM_0 - 2 && log10(Tbins_(i + 1)) < log_temp) ++i;
        while (j < nHbins_DIM_0 - 2 && log10(nHbins_(j + 1)) < log_density) ++j;

        // Logarithms of the Tbins and nHbins bounding the point
        Real log_T0 = log10(Tbins_(i));
        Real log_T1 = log10(Tbins_(i + 1));
        Real log_nH0 = log10(nHbins_(j));
        Real log_nH1 = log10(nHbins_(j + 1));

        // Compute weights for bilinear interpolation
        Real t = (log_temp - log_T0) / (log_T1 - log_T0);
        Real u = (log_density - log_nH0) / (log_nH1 - log_nH0);

        // Corner values in log space from the H_He cooling grid
        Real C00 = H_He_Cooling_(i, j);
        Real C10 = H_He_Cooling_(i + 1, j);
        Real C01 = H_He_Cooling_(i, j + 1);
        Real C11 = H_He_Cooling_(i + 1, j + 1);

        // Corner values in log space from the Metal cooling grid
        Real M00 = Metal_Cooling_(i, j);
        Real M10 = Metal_Cooling_(i + 1, j);
        Real M01 = Metal_Cooling_(i, j + 1);
        Real M11 = Metal_Cooling_(i + 1, j + 1);

        // Bilinear interpolation in log space
        Real prim_cooling =
            (1 - t) * (1 - u) * C00 +
            t * (1 - u) * C10 +
            (1 - t) * u * C01 +
            t * u * C11;
        Real metal_cooling =
            (1 - t) * (1 - u) * M00 +
            t * (1 - u) * M10 +
            (1 - t) * u * M01 +
            t * u * M11;

        lambda_cooling = prim_cooling + Z * metal_cooling;
      } // If density is higher than ceiling, switch to CIE
      else if (temp > Tfloor && temp < Tceil && nH >= nHceil) {
        // Convert input values to log space
        Real log_temp = log10(temp);

        // Locate index in Tbins
        int i = 0;
        while (i < Tbins_DIM_0 - 2 && log10(Tbins_(i + 1)) < log_temp) ++i;

        // Logarithms of the Tbins bounding the point
        Real log_T0 = log10(Tbins_(i));
        Real log_T1 = log10(Tbins_(i + 1));

        // Compute weight for linear interpolation
        Real t = (log_temp - log_T0) / (log_T1 - log_T0);

        // Bin values in log space from the H_He cooling grid
        Real C0 = H_He_Cooling_CIE_(i);
        Real C1 = H_He_Cooling_CIE_(i + 1);

        // Bin values in log space from the Metal cooling grid
        Real M0 = Metal_Cooling_CIE_(i);
        Real M1 = Metal_Cooling_CIE_(i + 1);

        // Linear Interpolation
        Real prim_cooling = C0 + t * (C1 - C0);
        Real metal_cooling = M0 + t * (M1 - M0);

        lambda_cooling = prim_cooling + Z * metal_cooling;
      }
      else if (temp < Tfloor && nH >= nHceil) {
        // for temperatures less than 100 K, use Koyama & Inutsuka (2002)
        lambda_cooling = Z*(2.0e-19*exp(-1.184e5/(temp + 1.0e3)) +
                            2.8e-28*sqrt(temp)*exp(-92.0/temp));
      }
 
      Real cooling_heating = FLT_MIN // add a tiny number
        + fabs(pow(w0(m,IDN,k,j,i),2) * lambda_cooling/cooling_unit);

      min_dt = fmin((eint/cooling_heating), min_dt);
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
