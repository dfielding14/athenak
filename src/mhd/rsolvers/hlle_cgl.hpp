#ifndef MHD_RSOLVERS_HLLE_CGL_HPP_
#define MHD_RSOLVERS_HLLE_CGL_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hlle_cgl.hpp
//! \brief HLLE Riemann solver for CGL.

#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "eos/cgl_passive.hpp"
#include "mhd/mhd.hpp"
#include "mhd/rsolvers/cgl_pressure_traction.hpp"

namespace mhd {

//----------------------------------------------------------------------------------------
//! \fn void HLLE
//! \brief The HLLE Riemann solver for hydrodynamics (both ideal gas and isothermal)

KOKKOS_INLINE_FUNCTION
void HLLE_CGL(TeamMember_t const &member, const EOS_Data &eos,
     const RegionIndcs &indcs,const DualArray1D<RegionSize> &size,const CoordData &coord,
     const int m, const int k, const int j, const int il, const int iu, const int ivx,
     const ScrArray2D<Real> &wl, const ScrArray2D<Real> &wr,
     const ScrArray2D<Real> &bl, const ScrArray2D<Real> &br, const DvceArray4D<Real> &bx,
     DvceArray5D<Real> flx, DvceArray4D<Real> ey, DvceArray4D<Real> ez,
     const bool record_pwork, DvceArray5D<Real> pflux) {
  if (eos.passive) {
    // Fluxes() has already supplied the native flow and electric fields.
    const int ivy = IVX + ((ivx-IVX)+1)%3;
    const int ivz = IVX + ((ivx-IVX)+2)%3;
    const int iby = ivy-IVX, ibz = ivz-IVX;
    par_for_inner(member, il, iu, [&](const int i) {
      const Real fd = flx(m,IDN,k,j,i);
      const bool left = fd >= 0.0;
      const Real rho = left ? wl(IDN,i) : wr(IDN,i);
      Real pp = fmax(left ? wl(IPR,i) : wr(IPR,i), eos.pfloor);
      Real pt = fmax(left ? wl(IPP,i) : wr(IPP,i), eos.pfloor);
      const Real by = left ? bl(iby,i) : br(iby,i);
      const Real bz = left ? bl(ibz,i) : br(ibz,i);
      const Real bn = bx(m,k,j,i);
      const Real bmag = sqrt(SQR(bn) + SQR(by) + SQR(bz));
      if (bmag <= eos.bfloor) {
        pp = TWO_3RDS*pt + ONE_3RD*pp;
        pt = pp;
      }
      const auto q = cgl::PassiveSpecificInvariants(rho, pp, pt,
                                                   fmax(bmag, eos.bfloor));
      flx(m,IEN,k,j,i) = fd*q.j;
      flx(m,IAN,k,j,i) = fd*q.a;
      if (record_pwork) {
        // Thermal stress diagnoses passive pressure work without accelerating the flow.
        MHDPrim1D l, r;
        l.e = wl(IPR,i); l.pp = wl(IPP,i); l.by = bl(iby,i); l.bz = bl(ibz,i);
        r.e = wr(IPR,i); r.pp = wr(IPP,i); r.by = br(iby,i); r.bz = br(ibz,i);
        MHDCons1D pressure, anisotropic;
        SingleStateLLF_CGLPressureTraction(l, r, bn, eos, pressure, anisotropic);
        pflux(m,ivx-IVX,k,j,i) = pressure.mx;
        pflux(m,ivy-IVX,k,j,i) = pressure.my;
        pflux(m,ivz-IVX,k,j,i) = pressure.mz;
        pflux(m,3+ivx-IVX,k,j,i) = anisotropic.mx;
        pflux(m,3+ivy-IVX,k,j,i) = anisotropic.my;
        pflux(m,3+ivz-IVX,k,j,i) = anisotropic.mz;
      }
    });
    return;
  }
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;
  int iby = ((ivx-IVX) + 1)%3;
  int ibz = ((ivx-IVX) + 2)%3;
  Real bfloor = eos.bfloor;

  par_for_inner(member, il, iu, [&](const int i) {
    //--- Step 1.  Create local references for L/R states (helps compiler vectorize)

    Real &wl_idn = wl(IDN,i);
    Real &wl_ivx = wl(ivx,i);
    Real &wl_ivy = wl(ivy,i);
    Real &wl_ivz = wl(ivz,i);
    Real &wl_ipr = wl(IPR,i);
    Real &wl_ipp = wl(IPP,i);
    Real &wl_iby = bl(iby,i);
    Real &wl_ibz = bl(ibz,i);

    Real &wr_idn = wr(IDN,i);
    Real &wr_ivx = wr(ivx,i);
    Real &wr_ivy = wr(ivy,i);
    Real &wr_ivz = wr(ivz,i);
    Real &wr_ipr = wr(IPR,i);
    Real &wr_ipp = wr(IPP,i);
    Real &wr_iby = br(iby,i);
    Real &wr_ibz = br(ibz,i);

    Real bxi = bx(m,k,j,i);

    // Magnetic pressure and magnitude on each side of the face.
    Real pbl = 0.5*(bxi*bxi + SQR(wl_iby) + SQR(wl_ibz));
    Real pbr = 0.5*(bxi*bxi + SQR(wr_iby) + SQR(wr_ibz));
    Real bmagl = sqrt(2.*pbl);
    Real bmagr = sqrt(2.*pbr);

    Real el = 0.5*wl_ipr + wl_ipp
              + 0.5*wl_idn*(SQR(wl_ivx)+SQR(wl_ivy)+SQR(wl_ivz)) + pbl;
    Real er = 0.5*wr_ipr + wr_ipp
              + 0.5*wr_idn*(SQR(wr_ivx)+SQR(wr_ivy)+SQR(wr_ivz)) + pbr;

    // Isotropize only weak states at fixed thermal energy; their transported A
    // uses the neighboring field to avoid injecting the arbitrary floor into it.
    Real fhl = 1.0, fhr = 1.0;
    if (bmagl > bfloor) {
      fhl += (wl_ipp - wl_ipr)/(2.0*pbl);
    } else {
      wl_ipr = TWO_3RDS*wl_ipp + ONE_3RD*wl_ipr;
      wl_ipp = wl_ipr;
    }
    if (bmagr > bfloor) {
      fhr += (wr_ipp - wr_ipr)/(2.0*pbr);
    } else {
      wr_ipr = TWO_3RDS*wr_ipp + ONE_3RD*wr_ipr;
      wr_ipp = wr_ipr;
    }

    //--- Step 3. Compute fast magnetosonic speed in the left and right states.

    const Real cl = eos.IdealMHDFastSpeed(wl_idn, wl_ipr, wl_ipp,
                                         bxi, wl_iby, wl_ibz, bfloor);
    const Real cr = eos.IdealMHDFastSpeed(wr_idn, wr_ipr, wr_ipp,
                                         bxi, wr_iby, wr_ibz, bfloor);

    //--- Step 4. Compute the max/min wave speeds from the left and right states.
    Real al = fmin((wr_ivx - cr),(wl_ivx - cl));
    Real ar = fmax((wr_ivx + cr),(wl_ivx + cl));

    // following min/max set to TINY_NUMBER to fix bug found in converging supersonic flow
    Real bp = ar > 0.0 ? ar : 1.0e-20;
    Real bm = al < 0.0 ? al : -1.0e-20;

    //--- Step 5. Compute L/R fluxes along the lines bm/bp: F_L - (S_L)U_L; F_R - (S_R)U_R

    Real vxl = wl_ivx - bm;
    Real vxr = wr_ivx - bp;

    MHDCons1D fl,fr;
    fl.d  = wl_idn*vxl;
    fr.d  = wr_idn*vxr;

    fl.mx = wl_idn*wl_ivx*vxl + pbl + wl_ipp - SQR(bxi)*fhl;
    fr.mx = wr_idn*wr_ivx*vxr + pbr + wr_ipp - SQR(bxi)*fhr;

    fl.my = wl_idn*wl_ivy*vxl - bxi*wl_iby*fhl;
    fr.my = wr_idn*wr_ivy*vxr - bxi*wr_iby*fhr;

    fl.mz = wl_idn*wl_ivz*vxl - bxi*wl_ibz*fhl;
    fr.mz = wr_idn*wr_ivz*vxr - bxi*wr_ibz*fhr;

    //Energy fluxes
    fl.e = el*vxl + wl_ivx*(wl_ipp + pbl);
    fr.e = er*vxr + wr_ivx*(wr_ipp + pbr);
    fl.e -= bxi*(bxi*wl_ivx + wl_iby*wl_ivy + wl_ibz*wl_ivz)*fhl;
    fr.e -= bxi*(bxi*wr_ivx + wr_iby*wr_ivy + wr_ibz*wr_ivz)*fhr;

    //B/E field fluxes for CT
    fl.by = wl_iby*vxl - bxi*wl_ivy;
    fr.by = wr_iby*vxr - bxi*wr_ivy;

    fl.bz = wl_ibz*vxl - bxi*wl_ivz;
    fr.bz = wr_ibz*vxr - bxi*wr_ivz;

    //--- Step 6. Compute the HLLE flux at interface.

    Real tmp=0.0;
    if (bp != bm) tmp = 0.5*(bp + bm)/(bp - bm);

    // Treat A/rho as the advected scalar, flux is mass_flx*A/rho from upwind side.
    Real fdtmp = 0.5*(fl.d  + fr.d ) + (fl.d  - fr.d )*tmp;
    const bool left = fdtmp >= 0.0;
    const Real rho = left ? wl_idn : wr_idn;
    const Real ppar = left ? wl_ipr : wr_ipr;
    const Real pperp = left ? wl_ipp : wr_ipp;
    const Real bmag = left ? bmagl : bmagr;
    Real anis;
    if (bmag > bfloor) {
      anis = cgl::SpecificAnisotropyForFlux(rho, ppar, pperp, bmag, SQR(bmag));
    } else {
      const Real bref = fmax(left ? bmagr : bmagl, bfloor);
      anis = 2.0*log(rho) - 3.0*log(bref);
    }
    Real fmutmp = fdtmp*anis;

    flx(m,IDN,k,j,i) = fdtmp;
    flx(m,ivx,k,j,i) = 0.5*(fl.mx + fr.mx) + (fl.mx - fr.mx)*tmp;
    flx(m,ivy,k,j,i) = 0.5*(fl.my + fr.my) + (fl.my - fr.my)*tmp;
    flx(m,ivz,k,j,i) = 0.5*(fl.mz + fr.mz) + (fl.mz - fr.mz)*tmp;
    flx(m,IEN,k,j,i) = 0.5*(fl.e + fr.e ) + (fl.e - fr.e)*tmp;
    flx(m,IAN,k,j,i) = fmutmp;
    ey(m,k,j,i) = -0.5*(fl.by + fr.by) - (fl.by - fr.by)*tmp;
    ez(m,k,j,i) =  0.5*(fl.bz + fr.bz) + (fl.bz - fr.bz)*tmp;

    if (record_pwork) {
      const Real danisl = (bmagl > bfloor) ? (wl_ipp - wl_ipr)/(2.0*pbl) : 0.0;
      const Real danisr = (bmagr > bfloor) ? (wr_ipp - wr_ipr)/(2.0*pbr) : 0.0;
      const Real alx = -SQR(bxi)*danisl;
      const Real aly = -bxi*wl_iby*danisl;
      const Real alz = -bxi*wl_ibz*danisl;
      const Real arx = -SQR(bxi)*danisr;
      const Real ary = -bxi*wr_iby*danisr;
      const Real arz = -bxi*wr_ibz*danisr;
      const Real plx = wl_ipp + alx;
      const Real ply = aly;
      const Real plz = alz;
      const Real prx = wr_ipp + arx;
      const Real pry = ary;
      const Real prz = arz;
      pflux(m,ivx-IVX,k,j,i) = 0.5*(plx + prx) + (plx - prx)*tmp;
      pflux(m,ivy-IVX,k,j,i) = 0.5*(ply + pry) + (ply - pry)*tmp;
      pflux(m,ivz-IVX,k,j,i) = 0.5*(plz + prz) + (plz - prz)*tmp;
      pflux(m,3+ivx-IVX,k,j,i) = 0.5*(alx + arx) + (alx - arx)*tmp;
      pflux(m,3+ivy-IVX,k,j,i) = 0.5*(aly + ary) + (aly - ary)*tmp;
      pflux(m,3+ivz-IVX,k,j,i) = 0.5*(alz + arz) + (alz - arz)*tmp;
    }
  });

  return;
}
} // namespace mhd
#endif // MHD_RSOLVERS_HLLE_CGL_HPP_
