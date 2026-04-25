#ifndef SRCTERMS_CGM_COOLING_HPP_
#define SRCTERMS_CGM_COOLING_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgm_cooling.hpp
//! \brief Shared CGM cooling/heating helpers.

#include <float.h>

#include "athena.hpp"
#include "cooling_tables.hpp"

namespace cgm_cooling {

struct CGMCoolingResult {
  Real temp;
  Real net_cooling_rate;
  Real cap_energy_change;
  Real cap_active;
};

template<typename View1D, typename View2D>
KOKKOS_INLINE_FUNCTION
CGMCoolingResult CalculateRates(const Real rho, const Real eint, const Real metallicity,
                                const Real x1v, const Real x2v, const Real x3v,
                                const Real dx1, const Real gm1, const Real temp_unit,
                                const Real nH_unit, const Real cooling_unit,
                                const Real heating_unit, const Real length_unit,
                                const Real h_rate, const Real h_norm,
                                const Real h_height, const Real h_radius,
                                const Real T_max, const Real Tfloor,
                                const Real Tceil, const Real nHfloor,
                                const Real nHceil, const View1D &Tbins,
                                const View1D &nHbins, const View2D &Metal_Cooling,
                                const View2D &H_He_Cooling,
                                const View1D &Metal_Cooling_CIE,
                                const View1D &H_He_Cooling_CIE) {
  constexpr Real X = 0.75;    // hydrogen mass fraction
  constexpr Real Zsol = 0.02; // solar metallicity

  const Real temp = temp_unit*eint/rho*gm1;
  const Real nH = X*nH_unit*rho;
  const Real Z = metallicity/Zsol;
  const Real log_temp = log10(temp);
  const Real log_nH = log10(nH);

  const Real m_cap = (temp >= T_max) ? 1.0 : 0.0;
  const Real m_lowT = (temp < Tfloor) ? 1.0 : 0.0;
  const Real m_Tin = (temp >= Tfloor && temp <= Tceil) ? 1.0 : 0.0;
  const Real m_Nin = (nH >= nHfloor && nH <= nHceil) ? 1.0 : 0.0;
  const Real m_PIE = m_Tin*m_Nin;
  const Real m_CIE = m_Tin;

  int iT = 0, jN = 0;
  while (iT < Tbins_DIM_0 - 2 && Tbins(iT + 1) < log_temp) ++iT;
  while (jN < nHbins_DIM_0 - 2 && nHbins(jN + 1) < log_nH) ++jN;

  const Real log_T0 = Tbins(iT);
  const Real log_T1 = Tbins(iT + 1);
  const Real t = (log_temp - log_T0)/(log_T1 - log_T0);
  const Real omt = 1.0 - t;

  const Real log_n0 = nHbins(jN);
  const Real log_n1 = nHbins(jN + 1);
  const Real u = (log_nH - log_n0)/(log_n1 - log_n0);
  const Real omu = 1.0 - u;

  const Real prim_PIE =
      omt*omu*H_He_Cooling(iT,   jN  ) +
      t  *omu*H_He_Cooling(iT+1, jN  ) +
      omt*u  *H_He_Cooling(iT,   jN+1) +
      t  *u  *H_He_Cooling(iT+1, jN+1);

  const Real metal_PIE =
      omt*omu*Metal_Cooling(iT,   jN  ) +
      t  *omu*Metal_Cooling(iT+1, jN  ) +
      omt*u  *Metal_Cooling(iT,   jN+1) +
      t  *u  *Metal_Cooling(iT+1, jN+1);

  const Real lambda_PIE = m_PIE*fma(Z, metal_PIE, prim_PIE);

  const Real C0 = H_He_Cooling_CIE(iT);
  const Real M0 = Metal_Cooling_CIE(iT);
  const Real prim_CIE = fma(t, H_He_Cooling_CIE(iT+1) - C0, C0);
  const Real metal_CIE = fma(t, Metal_Cooling_CIE(iT+1) - M0, M0);
  const Real lambda_CIE_tab = fma(Z, metal_CIE, prim_CIE);

  const Real e1 = exp(-1.184e5/(temp + 1.0e3));
  const Real e2 = exp(-92.0/temp);
  const Real lambda_lowT = Z*(2.0e-19*e1 + 2.8e-28*sqrt(temp)*e2);
  const Real lambda_CIE = m_CIE*lambda_CIE_tab + (1.0 - m_CIE)*(m_lowT*lambda_lowT);

  Real gamma_heating = 0.0;
  if (h_rate != 0.0 && h_norm != 0.0) {
    const Real R2 = fma(x1v, x1v, x2v*x2v);
    const Real R = sqrt(R2);
    const Real horz_falloff = exp(-R/h_radius);
    const Real vert_scale2 = h_height*h_height*(1.0 + R2/(h_radius*h_radius));
    const Real vert_falloff = exp(-(x3v*x3v)/vert_scale2);

    // Historical Gotham normalization: h_norm carries the calibrated dimensional
    // scale, while h_rate is the input heating coefficient in cgs units.
    gamma_heating = h_rate*h_norm*X*nH_unit*horz_falloff*vert_falloff;

    const Real m_hot = (temp > 1.0e4) ? 1.0 : 0.0;
    const Real inv_ratio = 1.0e4/temp;
    const Real damp_factor =
        m_hot*inv_ratio*inv_ratio*inv_ratio*inv_ratio*
              inv_ratio*inv_ratio*inv_ratio*inv_ratio + (1.0 - m_hot);
    gamma_heating *= damp_factor;
  }

  const Real dx_cgs = dx1*length_unit;
  const Real neutral_frac = 1.0 - 0.5*(1.0 + tanh((temp - 8.0e3)/1.5e3));
  const Real tau = neutral_frac*nH*1.0e-17*dx_cgs;
  const Real frac = exp(-tau);
  const Real lambda_cooling = (1.0 - frac)*lambda_CIE + frac*lambda_PIE;

  CGMCoolingResult result;
  result.temp = temp;
  result.net_cooling_rate =
      X*rho*(X*rho*(lambda_cooling/cooling_unit) - (gamma_heating/heating_unit));
  result.cap_energy_change = -((temp - T_max)/temp_unit)*rho/gm1;
  result.cap_active = m_cap;
  return result;
}

KOKKOS_INLINE_FUNCTION
Real CoolingTime(const CGMCoolingResult &result, const Real eint) {
  const Real denom = (result.net_cooling_rate >= 0.0)
                   ? result.net_cooling_rate + FLT_MIN
                   : result.net_cooling_rate - FLT_MIN;
  return (result.cap_active > 0.0) ? 0.0 : eint/denom;
}

} // namespace cgm_cooling

#endif // SRCTERMS_CGM_COOLING_HPP_
