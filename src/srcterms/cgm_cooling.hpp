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

#include <vector>

#include "athena.hpp"
#include "cooling_tables.hpp"
#include "cgm_cooling_terms.hpp"

namespace cgm_cooling {

constexpr Real X = 0.75;     // hydrogen mass fraction
constexpr Real Zsol = 0.02;  // solar metallicity

struct CGMCoolingResult {
  Real temp;
  Real net_cooling_rate;
  Real cap_active;
};

//! \fn TempTerms TempPart()
//! \brief The parts of the net rate that depend only on temperature [K]
template<typename View1D>
KOKKOS_INLINE_FUNCTION
TempTerms TempPart(const Real temp, const Real Tfloor, const Real Tceil,
                   const View1D &Tbins, const View1D &Metal_Cooling_CIE,
                   const View1D &H_He_Cooling_CIE) {
  TempTerms tt;
  tt.temp = temp;
  tt.below_table = (temp < Tfloor) ? 1.0 : 0.0;
  tt.in_table = (temp >= Tfloor && temp <= Tceil) ? 1.0 : 0.0;

  // Bins are uniform (to 2e-5 dex), so the bin index is arithmetic rather than a scan.
  // Clamping before the cast also maps NaN and +-inf to a valid bin.
  const Real log_temp = log10(temp);
  const Real dlogT = (Tbins(Tbins_DIM_0 - 1) - Tbins(0))/(Tbins_DIM_0 - 1);
  tt.iT = static_cast<int>(fmin(fmax(floor((log_temp - Tbins(0))/dlogT), 0.0),
                                Tbins_DIM_0 - 2.0));
  const Real log_T0 = Tbins(tt.iT);
  const Real log_T1 = Tbins(tt.iT + 1);
  tt.t = (log_temp - log_T0)/(log_T1 - log_T0);

  const Real C0 = H_He_Cooling_CIE(tt.iT);
  const Real M0 = Metal_Cooling_CIE(tt.iT);
  tt.prim_CIE = fma(tt.t, H_He_Cooling_CIE(tt.iT+1) - C0, C0);
  tt.metal_CIE = fma(tt.t, Metal_Cooling_CIE(tt.iT+1) - M0, M0);

  const Real e1 = exp(-1.184e5/(temp + 1.0e3));
  const Real e2 = exp(-92.0/temp);
  tt.lowT = 2.0e-19*e1 + 2.8e-28*sqrt(temp)*e2;

  const Real inv_ratio = 1.0e4/temp;
  tt.heat_damp = (temp > 1.0e4) ? inv_ratio*inv_ratio*inv_ratio*inv_ratio*
                                  inv_ratio*inv_ratio*inv_ratio*inv_ratio : 1.0;

  tt.neutral_frac = 1.0 - 0.5*(1.0 + tanh((temp - 8.0e3)/1.5e3));
  return tt;
}

//! \fn CellTerms CellPart()
//! \brief The parts of the net rate that depend only on the cell
template<typename View1D>
KOKKOS_INLINE_FUNCTION
CellTerms CellPart(const Real rho, const Real metallicity, const Real x1v,
                   const Real x2v, const Real x3v, const Real dx1, const Real nH_unit,
                   const Real length_unit, const Real h_rate, const Real h_norm,
                   const Real h_height, const Real h_radius, const Real nHfloor,
                   const Real nHceil, const View1D &nHbins) {
  CellTerms ct;
  ct.X_rho = X*rho;
  ct.nH = X*nH_unit*rho;
  ct.Z = metallicity/Zsol;
  ct.in_nH = (ct.nH >= nHfloor) ? 1.0 : 0.0;

  // Above the table's top row (nH = 1) reuse that row: PIE has converged to CIE there.
  const Real log_nH = log10(fmin(ct.nH, nHceil));
  const Real dlognH = (nHbins(nHbins_DIM_0 - 1) - nHbins(0))/(nHbins_DIM_0 - 1);
  ct.jN = static_cast<int>(fmin(fmax(floor((log_nH - nHbins(0))/dlognH), 0.0),
                                nHbins_DIM_0 - 2.0));
  const Real log_n0 = nHbins(ct.jN);
  const Real log_n1 = nHbins(ct.jN + 1);
  ct.u = (log_nH - log_n0)/(log_n1 - log_n0);

  ct.gamma0 = 0.0;
  if (h_rate != 0.0 && h_norm != 0.0) {
    const Real R2 = fma(x1v, x1v, x2v*x2v);
    const Real R = sqrt(R2);
    const Real horz_falloff = exp(-R/h_radius);
    const Real vert_scale2 = h_height*h_height*(1.0 + R2/(h_radius*h_radius));
    const Real vert_falloff = exp(-(x3v*x3v)/vert_scale2);

    // Historical Gotham normalization: h_norm carries the calibrated dimensional
    // scale, while h_rate is the input heating coefficient in cgs units.
    ct.gamma0 = h_rate*h_norm*X*nH_unit*horz_falloff*vert_falloff;
  }

  ct.dx_cgs = dx1*length_unit;
  return ct;
}

//! \fn Real NetRate()
//! \brief Net cooling rate (cooling - heating) per volume [code units] of a cell at the
//! temperature of tt
template<typename View2D>
KOKKOS_INLINE_FUNCTION
Real NetRate(const CellTerms &ct, const TempTerms &tt, const Real cooling_unit,
             const Real heating_unit, const View2D &Metal_Cooling,
             const View2D &H_He_Cooling) {
  const int iT = tt.iT;
  const int jN = ct.jN;
  const Real t = tt.t;
  const Real omt = 1.0 - t;
  const Real u = ct.u;
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

  const Real lambda_PIE = tt.in_table*ct.in_nH*fma(ct.Z, metal_PIE, prim_PIE);

  const Real lambda_CIE_tab = fma(ct.Z, tt.metal_CIE, tt.prim_CIE);
  const Real lambda_lowT = ct.Z*tt.lowT;
  const Real lambda_CIE = tt.in_table*lambda_CIE_tab +
                          (1.0 - tt.in_table)*(tt.below_table*lambda_lowT);

  const Real gamma_heating = ct.gamma0*tt.heat_damp;

  const Real tau = tt.neutral_frac*ct.nH*1.0e-17*ct.dx_cgs;
  const Real frac = exp(-tau);
  const Real lambda_cooling = (1.0 - frac)*lambda_CIE + frac*lambda_PIE;

  return ct.X_rho*(ct.X_rho*(lambda_cooling/cooling_unit) - (gamma_heating/heating_unit));
}

//! \fn CGMCoolingResult CalculateRates()
//! \brief Net cooling rate of a cell at internal energy eint: NetRate of its CellPart
//! and TempPart
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
  const Real temp = temp_unit*eint/rho*gm1;
  const CellTerms ct = CellPart(rho, metallicity, x1v, x2v, x3v, dx1, nH_unit,
                                length_unit, h_rate, h_norm, h_height, h_radius,
                                nHfloor, nHceil, nHbins);
  const TempTerms tt = TempPart(temp, Tfloor, Tceil, Tbins, Metal_Cooling_CIE,
                                H_He_Cooling_CIE);

  CGMCoolingResult result;
  result.temp = temp;
  result.net_cooling_rate = NetRate(ct, tt, cooling_unit, heating_unit, Metal_Cooling,
                                    H_He_Cooling);
  result.cap_active = (temp >= T_max) ? 1.0 : 0.0;
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

// CalculateRates is not linear in T within a table bin: chords between table nodes miss
// its slope by O(bin width), 1.2% in the isobaric TI growth rate at 3000 K. Sub-nodes
// cut slope and value errors by kNodesPerBin and its square.
constexpr int kNodesPerBin = 10;

//! Integrator nodes and the temperature terms at both ends of each segment between
//! them (cgm_cooling::TabulateNodes): segment k = [temp[k - k_lo], temp[k - k_lo + 1]],
//! with ends lo[k - k_lo] and hi[k - k_lo].
struct CoolingNodes {
  int k_lo;                       // index k of the first node and segment
  Real logT_0, dlogT;             // node k sits at log10 T = logT_0 + k*dlogT
  std::vector<Real> temp;         // node temperatures [K]
  std::vector<TempTerms> lo, hi;  // terms at the low and high end of each segment
};

//! \fn CoolingNodes TabulateNodes()
//! \brief Tabulates, on the host, the nodes T_k = 10^(Tbins(0) + k*dlogT) of
//! IntegrateTemperature, kNodesPerBin per table bin, from just below T_lo to just above
//! T_hi, with the temperature terms at both ends of every segment [T_k, T_k+1].
//! CalculateRates jumps at the table edges Tfloor (node 0) and Tceil (node kTop), so a
//! segment below, inside or above the table takes both ends on its own side of them,
//! nudged by eps against T -> e -> T round-off.
template<typename View1D>
inline CoolingNodes TabulateNodes(const Real T_lo, const Real T_hi, const Real Tfloor,
                                  const Real Tceil, const View1D &Tbins,
                                  const View1D &Metal_Cooling_CIE,
                                  const View1D &H_He_Cooling_CIE) {
  constexpr int kTop = (Tbins_DIM_0 - 1)*kNodesPerBin;  // node index of Tceil
  constexpr Real eps = 1.0e-12;
  CoolingNodes nodes;
  nodes.logT_0 = Tbins(0);
  nodes.dlogT = (Tbins(Tbins_DIM_0 - 1) - nodes.logT_0)/kTop;
  nodes.k_lo = static_cast<int>(floor((log10(T_lo) - nodes.logT_0)/nodes.dlogT)) - 1;
  const int k_hi = static_cast<int>(ceil((log10(T_hi) - nodes.logT_0)/nodes.dlogT)) + 1;

  for (int k = nodes.k_lo; k <= k_hi; ++k) {
    nodes.temp.push_back(pow(10.0, nodes.logT_0 + k*nodes.dlogT));
  }
  for (int k = nodes.k_lo; k < k_hi; ++k) {
    auto on_side = [&](const Real T) {  // T as seen from segment k
      if (k < 0) return fmin(T, Tfloor*(1.0 - eps));        // below the table
      if (k >= kTop) return fmax(T, Tceil*(1.0 + eps));     // above the table
      return fmin(fmax(T, Tfloor*(1.0 + eps)), Tceil*(1.0 - eps));
    };
    const Real Ta = nodes.temp[k - nodes.k_lo];
    const Real Tb = nodes.temp[k - nodes.k_lo + 1];
    nodes.lo.push_back(TempPart(on_side(Ta), Tfloor, Tceil, Tbins, Metal_Cooling_CIE,
                                H_He_Cooling_CIE));
    nodes.hi.push_back(TempPart(on_side(Tb), Tfloor, Tceil, Tbins, Metal_Cooling_CIE,
                                H_He_Cooling_CIE));
  }
  return nodes;
}

//! \fn Real IntegrateTemperature()
//! \brief Returns T [K] after cooling/heating from T0 for dt at fixed density,
//! metallicity, position and dx (the cell's CellTerms). The net rate is taken linear in
//! T between the nodes of TabulateNodes, on which de/dt = -R integrates exactly: the
//! Townsend (2009) scheme generalized as in Tan & Fielding (2023, App. C), simplified.
//! Each segment has its own end rates, so where the rate jumps at a table edge T rests
//! on the edge if the rate beyond it pushes back. T stays in [T_lo, T_hi]. Costs two
//! NetRate evaluations per segment visited.
template<typename RealArray, typename TermsArray, typename View2D>
KOKKOS_INLINE_FUNCTION
Real IntegrateTemperature(const Real T0, const Real dt, const Real c_V, const Real T_lo,
                          const Real T_hi, const CellTerms &cell, const int k_lo,
                          const Real logT_0, const Real dlogT, const RealArray &node_T,
                          const TermsArray &seg_lo, const TermsArray &seg_hi,
                          const Real cooling_unit, const Real heating_unit,
                          const View2D &Metal_Cooling, const View2D &H_He_Cooling) {
  // net rate at the low and high end of segment k
  auto rate_lo = [&](const int k) {
    return NetRate(cell, seg_lo(k - k_lo), cooling_unit, heating_unit, Metal_Cooling,
                   H_He_Cooling);
  };
  auto rate_hi = [&](const int k) {
    return NetRate(cell, seg_hi(k - k_lo), cooling_unit, heating_unit, Metal_Cooling,
                   H_He_Cooling);
  };

  Real T = fmin(fmax(T0, T_lo), T_hi);
  int k = static_cast<int>(floor((log10(T) - logT_0)/dlogT));  // T in [T_k, T_k+1]
  Real Ta = node_T(k - k_lo), Tb = node_T(k + 1 - k_lo);
  Real Ra = rate_lo(k), Rb = rate_hi(k);
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
        // step into the next segment. Its rate at this node equals Re except at a table
        // edge, where T rests on the edge if the rate beyond it pushes T back.
        k += down ? -1 : 1;
        Ta = node_T(k - k_lo);
        Tb = node_T(k + 1 - k_lo);
        Ra = rate_lo(k);
        Rb = rate_hi(k);
        R0 = down ? Rb : Ra;
        if (down ? (R0 <= 0.0) : (R0 >= 0.0)) break;
        continue;
      }
    }
    T -= R0*rem/c_V*Phi(s*rem/c_V);  // ends inside this segment
    break;
  }
  return T;
}

} // namespace cgm_cooling

#endif // SRCTERMS_CGM_COOLING_HPP_
