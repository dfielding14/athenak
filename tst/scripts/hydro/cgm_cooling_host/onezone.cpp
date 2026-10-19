//! \file onezone.cpp
//! \brief Host one-zone regression test of cgm_cooling::IntegrateTemperature with the
//! athinput.gotham units and heating. Reference temperatures come from an RK4
//! integration of the same rates. Exits nonzero if any case is off by more than 0.5%.

#include <cstdio>

#include "cgm_cooling.hpp"

namespace {
struct Table1D {
  const Real *a;
  Real operator()(const int i) const { return a[i]; }
};
struct Table2D {
  const Real *a;
  int n1;
  Real operator()(const int i, const int j) const { return a[i*n1 + j]; }
};
}  // namespace

int main() {
  // athinput.gotham <units>, with the constants and formulas of src/units/units.*
  const Real length_unit = 3.0856775809623245e+21, mass_cgs = 3.036951775493658e+39;
  const Real time_cgs = 3.15576e+15, mu = 0.62;
  const Real amu_cgs = 1.67262192369e-24, kb_cgs = 1.3806488e-16, yr_cgs = 3.15576e+7;
  const Real vel_cgs = length_unit/time_cgs;
  const Real density_cgs = mass_cgs/(length_unit*length_unit*length_unit);
  const Real pressure_cgs = density_cgs*vel_cgs*vel_cgs;
  const Real temp_unit = vel_cgs*vel_cgs*mu*amu_cgs/kb_cgs;
  const Real nH_unit = density_cgs/amu_cgs;
  const Real cooling_unit = pressure_cgs/time_cgs/nH_unit/nH_unit;
  const Real heating_unit = pressure_cgs/time_cgs/nH_unit;

  // athinput.gotham <hydro>; tfloor = 10 K, finest-level dx, Z = Zsol/3 mass fraction
  const Real gm1 = 1.66666667 - 1.0;
  const Real h_rate = 2e-26, h_norm = 1347.6541799847946;
  const Real h_height = 0.04950000000848683, h_radius = 0.9402361259893334;
  const Real T_lo = 10.0, T_hi = 5.0e8;
  const Real metallicity = 0.02/3.0;
  const Real dx1 = 400.0/384.0/256.0;

  const Table1D Tbins{Tbins_ARR}, nHbins{nHbins_ARR};
  const Table1D Metal_CIE{Metal_Cooling_CIE_ARR}, H_He_CIE{H_He_Cooling_CIE_ARR};
  const Table2D Metal{Metal_Cooling_ARR, Metal_Cooling_DIM_1};
  const Table2D H_He{H_He_Cooling_ARR, H_He_Cooling_DIM_1};
  const Real Tfloor = pow(10.0, Tbins_ARR[0]);
  const Real Tceil = pow(10.0, Tbins_ARR[Tbins_DIM_0 - 1]);
  const Real nHfloor = pow(10.0, nHbins_ARR[0]);
  const Real nHceil = pow(10.0, nHbins_ARR[nHbins_DIM_0 - 1]);

  struct Case { Real nH, x1, T0, dt_yr, T_ref; };  // nH [cm^-3], x1 [kpc], T [K]
  const Case cases[] = {
    {100.0, 3.0, 1.0e5, 300.0, 1.154e4},   // stiff: t_cool ~ 64 yr
    {10.0, 3.0, 1.0e6, 3.0e4, 1.061e4},
    {62.4, 0.0, 2.0e4, 1.0e4, 1.003e4},
    {1.0e-3, 50.0, 150.0, 1.0e6, 1084.0},  // photoheated from below the table
    {1.0e-3, 50.0, 150.0, 1.0e8, 1.023e4}, // ... up to equilibrium
    {1.0e-2, 3.0, 1.0e6, 300.0, 1.0e6},    // t_cool >> dt
  };

  int nfail = 0;
  for (const Case &c : cases) {
    const Real rho = c.nH/(0.75*nH_unit);  // nH = X*nH_unit*rho, as in CalculateRates
    const Real T1 = cgm_cooling::IntegrateTemperature(
        c.T0, c.dt_yr*yr_cgs/time_cgs, rho, metallicity, c.x1, 0.0, 0.0, dx1, gm1,
        temp_unit, nH_unit, cooling_unit, heating_unit, length_unit, h_rate, h_norm,
        h_height, h_radius, T_lo, T_hi, Tfloor, Tceil, nHfloor, nHceil, Tbins, nHbins,
        Metal, H_He, Metal_CIE, H_He_CIE);
    const Real err = fabs(T1/c.T_ref - 1.0);
    const bool ok = (err < 5.0e-3);
    if (!ok) ++nfail;
    std::printf("%s nH=%-6g x1=%-4g T0=%-6g dt=%-6g yr: T=%.4e K (ref %.4e, err %.1e)\n",
                ok ? "ok  " : "FAIL", c.nH, c.x1, c.T0, c.dt_yr, T1, c.T_ref, err);
  }
  return nfail;
}
