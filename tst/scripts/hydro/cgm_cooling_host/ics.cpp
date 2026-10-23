//! \file ics.cpp
//! \brief Compare the IC generator against GOTHAM's compiled gravity and PIE cooling.

#include <iomanip>
#include <iostream>

#define SQR(x) ((x)*(x))
#include "galaxy_potential.hpp"
#include "cgm_cooling.hpp"

struct Table1D {
  const Real *data;
  Real operator()(int i) const { return data[i]; }
};
struct Table2D {
  const Real *data;
  Real operator()(int i, int j) const { return data[i*nHbins_DIM_0 + j]; }
};

int main() {
  Real G, rs, rhos, mass, a, b, r200, rhom;
  if (!(std::cin >> G >> rs >> rhos >> mass >> a >> b >> r200 >> rhom)) return 1;
  const Table1D Tbins{Tbins_ARR}, nHbins{nHbins_ARR};
  const Table1D metal_cie{Metal_Cooling_CIE_ARR}, prim_cie{H_He_Cooling_CIE_ARR};
  const Table2D metal{Metal_Cooling_ARR}, prim{H_He_Cooling_ARR};
  const Real Tmin = pow(10.0, Tbins_ARR[0]);
  const Real Tmax = pow(10.0, Tbins_ARR[Tbins_DIM_0-1]);
  const Real nmin = pow(10.0, nHbins_ARR[0]);
  const Real nmax = pow(10.0, nHbins_ARR[nHbins_DIM_0-1]);
  Real x, y, z, T, nH, Z;
  std::cout << std::setprecision(17);
  while (std::cin >> x >> y >> z >> T >> nH >> Z) {
    const Real phi = GalaxyPot(x, y, z, G, rs, rhos, mass, a, b,
                               GalaxyPotOuterCoef(r200), rhom);
    // dx=0 and heating=0 select the optically thin limit used by the 1D flow.
    const auto cell = cgm_cooling::CellPart(nH/cgm_cooling::X, Z*0.02,
        x, y, z, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, nmin, nmax, nHbins);
    const auto temp = cgm_cooling::TempPart(T, Tmin, Tmax, Tbins,
                                          metal_cie, prim_cie);
    const Real lambda = cgm_cooling::NetRate(cell, temp, 1.0, 1.0, metal, prim)/(nH*nH);
    std::cout << phi << ' ' << lambda << '\n';
  }
}
