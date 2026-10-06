#ifndef MHD_RSOLVERS_CGL_PRESSURE_TRACTION_HPP_
#define MHD_RSOLVERS_CGL_PRESSURE_TRACTION_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_pressure_traction.hpp
//! \brief Shared CGL pressure traction without relativistic solver dependencies.

#include "athena.hpp"
#include "eos/eos.hpp"

namespace mhd {

//----------------------------------------------------------------------------------------
//! \brief Isolate the physical CGL pressure traction retained by an LLF face flux.

KOKKOS_INLINE_FUNCTION
void SingleStateLLF_CGLPressureTraction(const MHDPrim1D &wl, const MHDPrim1D &wr,
                                        const Real &bxi, const EOS_Data &eos,
                                        MHDCons1D &pressure, MHDCons1D &anisotropic) {
  pressure.mx = pressure.my = pressure.mz = 0.0;
  anisotropic.mx = anisotropic.my = anisotropic.mz = 0.0;
  Real pparl = wl.e;
  Real pperpl = wl.pp;
  Real pparr = wr.e;
  Real pperpr = wr.pp;
  const Real bsql = SQR(bxi) + SQR(wl.by) + SQR(wl.bz);
  const Real bsqr = SQR(bxi) + SQR(wr.by) + SQR(wr.bz);
  const Real bmagl = sqrt(bsql);
  const Real bmagr = sqrt(bsqr);
  if (bmagl <= eos.bfloor) {
    pparl = TWO_3RDS*pperpl + ONE_3RD*pparl;
    pperpl = pparl;
  }
  if (bmagr <= eos.bfloor) {
    pparr = TWO_3RDS*pperpr + ONE_3RD*pparr;
    pperpr = pparr;
  }
  const Real dl = (bmagl > eos.bfloor) ? (pperpl - pparl)/bsql : 0.0;
  const Real dr = (bmagr > eos.bfloor) ? (pperpr - pparr)/bsqr : 0.0;
  const Real alx = -SQR(bxi)*dl;
  const Real aly = -bxi*wl.by*dl;
  const Real alz = -bxi*wl.bz*dl;
  const Real arx = -SQR(bxi)*dr;
  const Real ary = -bxi*wr.by*dr;
  const Real arz = -bxi*wr.bz*dr;
  anisotropic.mx = 0.5*(alx + arx);
  anisotropic.my = 0.5*(aly + ary);
  anisotropic.mz = 0.5*(alz + arz);
  pressure.mx = 0.5*(pperpl + alx + pperpr + arx);
  pressure.my = anisotropic.my;
  pressure.mz = anisotropic.mz;
}

}  // namespace mhd

#endif  // MHD_RSOLVERS_CGL_PRESSURE_TRACTION_HPP_
