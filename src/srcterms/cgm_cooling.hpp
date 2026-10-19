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
  // Above the table's top row (nH = 1) reuse that row: PIE has converged to CIE there.
  const Real log_nH = fmin(log10(nH), log10(nHceil));

  const Real m_cap = (temp >= T_max) ? 1.0 : 0.0;
  const Real m_lowT = (temp < Tfloor) ? 1.0 : 0.0;
  const Real m_Tin = (temp >= Tfloor && temp <= Tceil) ? 1.0 : 0.0;
  const Real m_Nin = (nH >= nHfloor) ? 1.0 : 0.0;
  const Real m_PIE = m_Tin*m_Nin;
  const Real m_CIE = m_Tin;

  // Bins are uniform (to 2e-5 dex), so the bin index is arithmetic rather than a scan.
  // Clamping before the cast also maps NaN and +-inf to a valid bin.
  const Real dlogT = (Tbins(Tbins_DIM_0 - 1) - Tbins(0))/(Tbins_DIM_0 - 1);
  const Real dlognH = (nHbins(nHbins_DIM_0 - 1) - nHbins(0))/(nHbins_DIM_0 - 1);
  const int iT = static_cast<int>(fmin(fmax(floor((log_temp - Tbins(0))/dlogT), 0.0),
                                       Tbins_DIM_0 - 2.0));
  const int jN = static_cast<int>(fmin(fmax(floor((log_nH - nHbins(0))/dlognH), 0.0),
                                       nHbins_DIM_0 - 2.0));

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

// phi(x) = (1 - e^-x)/x and psi(u) = ln(1 + u)/u, both -> 1 as their argument -> 0
KOKKOS_INLINE_FUNCTION
Real Phi(const Real x) {
  return (fabs(x) < 1.0e-8) ? 1.0 : -expm1(-x)/x;
}

KOKKOS_INLINE_FUNCTION
Real Psi(const Real u) {
  return (fabs(u) < 1.0e-8) ? 1.0 : log1p(u)/u;
}

//! \fn Real IntegrateTemperature()
//! \brief Returns T [K] after cooling/heating from T0 for dt at fixed rho, Z, position
//! and dx. The net rate R(T) (CalculateRates at e = c_V*T) is taken linear in T between
//! nodes T_k = 10^(Tbins(0) + k*dlogT), on which de/dt = -R integrates exactly: the
//! Townsend (2009) scheme generalized as in Tan & Fielding (2023, App. C), simplified.
//! CalculateRates jumps at the table edges Tfloor (node 0) and Tceil (node N-1), so these
//! are segment boundaries with one-sided rates; T rests on an edge if the rate beyond it
//! pushes back. T stays in [T_lo, T_hi]. Costs about one rate evaluation per node
//! crossed.
template<typename View1D, typename View2D>
KOKKOS_INLINE_FUNCTION
Real IntegrateTemperature(const Real T0, const Real dt, const Real rho,
                          const Real metallicity, const Real x1v, const Real x2v,
                          const Real x3v, const Real dx1, const Real gm1,
                          const Real temp_unit, const Real nH_unit,
                          const Real cooling_unit, const Real heating_unit,
                          const Real length_unit, const Real h_rate, const Real h_norm,
                          const Real h_height, const Real h_radius, const Real T_lo,
                          const Real T_hi, const Real Tfloor, const Real Tceil,
                          const Real nHfloor, const Real nHceil, const View1D &Tbins,
                          const View1D &nHbins, const View2D &Metal_Cooling,
                          const View2D &H_He_Cooling, const View1D &Metal_Cooling_CIE,
                          const View1D &H_He_Cooling_CIE) {
  const Real c_V = rho/(gm1*temp_unit);  // de/dT
  auto rate = [&](const Real T) {
    return CalculateRates(rho, c_V*T, metallicity, x1v, x2v, x3v, dx1, gm1, temp_unit,
                          nH_unit, cooling_unit, heating_unit, length_unit, h_rate,
                          h_norm, h_height, h_radius, T_hi, Tfloor, Tceil, nHfloor,
                          nHceil, Tbins, nHbins, Metal_Cooling, H_He_Cooling,
                          Metal_Cooling_CIE, H_He_Cooling_CIE).net_cooling_rate;
  };
  const Real logT_0 = Tbins(0);
  const Real dlogT = (Tbins(Tbins_DIM_0 - 1) - logT_0)/(Tbins_DIM_0 - 1);
  auto node = [&](const int k) { return pow(10.0, logT_0 + k*dlogT); };
  // segment k lies below the table (k < 0), in it, or above it (k >= N-1); rates are
  // taken on the segment's own side of the edges, with a margin against T->e->T round-off
  auto side = [](const int k) { return (k < 0) ? 0 : ((k < Tbins_DIM_0 - 1) ? 1 : 2); };
  auto side_rate = [&](const Real T, const int k) {
    constexpr Real eps = 1.0e-12;
    if (side(k) == 0) return rate(fmin(T, Tfloor*(1.0 - eps)));
    if (side(k) == 2) return rate(fmax(T, Tceil*(1.0 + eps)));
    return rate(fmin(fmax(T, Tfloor*(1.0 + eps)), Tceil*(1.0 - eps)));
  };

  Real T = fmin(fmax(T0, T_lo), T_hi);
  int k = static_cast<int>(floor((log10(T) - logT_0)/dlogT));  // T in [T_k, T_k+1]
  Real Ta = node(k), Tb = node(k + 1);
  Real Ra = side_rate(Ta, k), Rb = side_rate(Tb, k);
  Real R0 = Ra + (Rb - Ra)/(Tb - Ta)*(T - Ta);
  Real rem = dt;
  while (R0 != 0.0) {
    const Real s = (Rb - Ra)/(Tb - Ta);
    const bool down = (R0 > 0.0);  // net cooling lowers T
    const Real Te = down ? fmax(Ta, T_lo) : fmin(Tb, T_hi);
    // interpolate from the edge's own node, so Re is exact when the edge is a node
    const Real Re = down ? Ra + s*(Te - Ta) : Rb + s*(Te - Tb);
    if (Re*R0 > 0.0) {  // R keeps its sign up to the edge: T reaches it after te
      const Real te = c_V*(T - Te)/R0*Psi(s*(Te - T)/R0);
      if (te < rem) {
        rem -= te;
        T = Te;
        if (down ? (Ta <= T_lo) : (Tb >= T_hi)) break;  // at the floor or the cap
        const int kn = down ? k - 1 : k + 1;
        if (side(kn) != side(k)) {
          // crossing a table edge: re-evaluate both ends on the new side, and rest on the
          // edge if the rate there pushes T back
          k = kn;
          Ta = node(k); Tb = node(k + 1);
          Ra = side_rate(Ta, k); Rb = side_rate(Tb, k);
          R0 = down ? Rb : Ra;
          if (down ? (R0 <= 0.0) : (R0 >= 0.0)) break;
          continue;
        }
        k = kn;
        if (down) {
          Tb = Ta; Rb = Ra;
          Ta = node(k); Ra = side_rate(Ta, k);
        } else {
          Ta = Tb; Ra = Rb;
          Tb = node(k + 1); Rb = side_rate(Tb, k);
        }
        R0 = Re;  // exact node value, so round-off cannot flip the direction
        continue;
      }
    }
    T -= R0*rem/c_V*Phi(s*rem/c_V);  // ends inside this bin
    break;
  }
  return T;
}

} // namespace cgm_cooling

#endif // SRCTERMS_CGM_COOLING_HPP_
