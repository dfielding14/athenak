//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_mhd.cpp
//! \brief derived class that implements ideal gas EOS in nonrelativistic mhd

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "mhd/mhd.hpp"
#include "eos.hpp"
#include "eos/cgl_physics.hpp"
#include "eos/ideal_c2p_mhd.hpp"
#include "diffusion/cgl_landau_fluid.hpp"

namespace {

void RequireNonnegativeCGLParameter(const char *name, const Real value) {
  if (value < 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<mhd>/" << name << " must be nonnegative" << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

Real ParseCGLFirehoseThreshold(ParameterInput *pin) {
  Real legacy_threshold = 2.0;
  if (pin->DoesParameterExist("mhd", "cgl_firehose_threshold")) {
    const std::string policy = pin->GetString("mhd", "cgl_firehose_threshold");
    if (policy == "oblique") {
      legacy_threshold = 1.4;
    } else if (policy != "parallel") {
      std::cout << "### FATAL ERROR: <mhd>/cgl_firehose_threshold must be "
                << "oblique or parallel" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (pin->DoesParameterExist("mhd", "firehose_threshold") &&
        pin->GetReal("mhd", "firehose_threshold") != legacy_threshold) {
      std::cout << "### FATAL ERROR: <mhd>/firehose_threshold conflicts with "
                << "cgl_firehose_threshold; remove the legacy parameter" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  return pin->GetOrAddReal("mhd", "firehose_threshold", legacy_threshold);
}

void RequireCGLThreshold(const char *name, const Real value, const Real lower,
                         const bool inclusive) {
  if (!std::isfinite(value) || (inclusive ? value < lower : value <= lower)) {
    std::cout << "### FATAL ERROR: <mhd>/" << name << " must be finite and "
              << (inclusive ? ">= " : "> ") << lower << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

} // namespace

//----------------------------------------------------------------------------------------
// ctor: also calls EOS base class constructor

CGLMHD::CGLMHD(MeshBlockPack *pp, ParameterInput *pin) :
    EquationOfState("mhd", pp, pin) {
  eos_data.is_ideal = true;
  eos_data.is_cgl = true;
  if (!(eos_data.bfloor > 0.0)) {
    std::cout << "### FATAL ERROR: <mhd>/bfloor must be positive for CGL"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (sizeof(Real) == sizeof(float) &&
      eos_data.bfloor*eos_data.bfloor*eos_data.bfloor < FLT_MIN) {
    std::cout << "### FATAL ERROR: CGL bfloor cubed underflows in single precision; "
              << "set <mhd>/bfloor explicitly so bfloor^3 >= FLT_MIN" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  eos_data.mlim = false;
  eos_data.flim = false;
  eos_data.backup_lim = false;
  eos_data.coll = false;
  eos_data.gamma = pin->GetOrAddReal("mhd", "gamma", 5.0/3.0);
  eos_data.nu_coll = 0.0;
  eos_data.lim_coll = 0.0;
  eos_data.firehose_threshold = ParseCGLFirehoseThreshold(pin);
  eos_data.mirror_threshold = pin->GetOrAddReal("mhd", "mirror_threshold", 1.0);
  eos_data.mirror_backup_factor = pin->GetOrAddReal("mhd", "mirror_backup_factor", 2.0);
  eos_data.firehose_backup_factor =
      pin->GetOrAddReal("mhd", "firehose_backup_factor", 1.0);
  eos_data.limiter_backup_nu = pin->GetOrAddReal("mhd", "limiter_backup_nu", 1.0e10);
  RequireCGLThreshold("firehose_threshold", eos_data.firehose_threshold, 0.0, false);
  RequireCGLThreshold("mirror_threshold", eos_data.mirror_threshold, 0.0, false);
  RequireCGLThreshold("mirror_backup_factor", eos_data.mirror_backup_factor, 1.0, true);
  RequireCGLThreshold("firehose_backup_factor", eos_data.firehose_backup_factor,
                       1.0, true);
  RequireCGLThreshold("limiter_backup_nu", eos_data.limiter_backup_nu, 0.0, true);
  eos_data.sigma_max = pin->GetOrAddReal("mhd","sigma_max",(FLT_MAX));

  eos_data.passive = pin->GetOrAddBoolean("mhd", "passive", false);
  eos_data.iso_cs = eos_data.passive ? pin->GetReal("mhd", "iso_sound_speed") : 0.0;
  if (eos_data.passive && (!(eos_data.iso_cs > 0.0) ||
                          !std::isfinite(eos_data.iso_cs))) {
    std::cout << "### FATAL ERROR: passive CGL requires finite positive iso_sound_speed"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  eos_data.mlim = pin->GetOrAddBoolean("mhd", "mirror_limiter", false);
  eos_data.flim = pin->GetOrAddBoolean("mhd", "firehose_limiter", false);
  if (pin->GetOrAddBoolean("mhd", "limiter_hardwall", false)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<mhd>/limiter_hardwall=true is no longer supported; remove it and "
              << "set limiter_nu_coll=1e10 for stiff soft-threshold relaxation once "
              << "per timestep. Set backup_limiters explicitly for backup walls."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  eos_data.backup_lim = pin->GetOrAddBoolean("mhd", "backup_limiters", false);
  if (eos_data.backup_lim && !(eos_data.mlim || eos_data.flim)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<mhd>/backup_limiters requires mirror_limiter or firehose_limiter"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (eos_data.mlim || eos_data.flim) {
    eos_data.lim_coll = pin->GetReal("mhd", "limiter_nu_coll");
    RequireNonnegativeCGLParameter("limiter_nu_coll", eos_data.lim_coll);
    eos_data.coll = true;
  }

  eos_data.nu_coll = pin->GetOrAddReal("mhd", "nu_coll", 0.0);
  RequireNonnegativeCGLParameter("nu_coll", eos_data.nu_coll);
  eos_data.coll = eos_data.coll || (eos_data.nu_coll > 0.0);
}

//----------------------------------------------------------------------------------------
//! \!fn void ConsToPrim()
//! \brief Converts conserved into primitive variables.  Operates over range of cells
//! given in argument list.

void CGLMHD::ConsToPrim(DvceArray5D<Real> &cons, const DvceFaceFld4D<Real> &b,
                          DvceArray5D<Real> &prim, DvceArray5D<Real> &bcc,
                          const bool only_testfloors,
                          const int il, const int iu, const int jl, const int ju,
                          const int kl, const int ku) {
  const int nmhd  = pmy_pack->pmhd->nmhd;
  const int nscal = pmy_pack->pmhd->nscalars;
  const int nmb = pmy_pack->nmb_thispack;
  const EOS_Data eos = eos_data;
  auto fofc_ = pmy_pack->pmhd->fofc;

  const int ni   = (iu - il + 1);
  const int nji  = (ju - jl + 1)*ni;
  const int nkji = (ku - kl + 1)*nji;
  const int nmkji = nmb*nkji;

  int nfloord_=0, nfloore_=0, nfloort_=0;
  Kokkos::parallel_reduce("mhd_c2p",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, int &sumd, int &sume, int &sumt) {
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/ni;
    int i = (idx - m*nkji - k*nji - j*ni) + il;
    j += jl;
    k += kl;

    // load single state conserved variables
    MHDCons1D u;
    u.d  = cons(m,IDN,k,j,i);
    u.mx = cons(m,IM1,k,j,i);
    u.my = cons(m,IM2,k,j,i);
    u.mz = cons(m,IM3,k,j,i);
    u.e  = cons(m,IEN,k,j,i);
    u.mu = cons(m,IAN,k,j,i);

    // load cell-centered fields into conserved state
    // use input CC fields if only testing floors with FOFC
    if (only_testfloors && !eos.passive) {
      u.bx = bcc(m,IBX,k,j,i);
      u.by = bcc(m,IBY,k,j,i);
      u.bz = bcc(m,IBZ,k,j,i);
    // else use simple linear average of face-centered fields
    } else {
      u.bx = 0.5*(b.x1f(m,k,j,i) + b.x1f(m,k,j,i+1));
      u.by = 0.5*(b.x2f(m,k,j,i) + b.x2f(m,k,j+1,i));
      u.bz = 0.5*(b.x3f(m,k,j,i) + b.x3f(m,k+1,j,i));
    }

    // Preserve invalid input evidence: C2P may repair the local energy or A.
    const bool nonfinite_input = only_testfloors && !eos.passive &&
        (!Kokkos::isfinite(u.mu) ||
         !Kokkos::isfinite(u.e - 0.5*(SQR(u.mx) + SQR(u.my) + SQR(u.mz))/u.d
                              - 0.5*(SQR(u.bx) + SQR(u.by) + SQR(u.bz))));

    // call c2p function
    // (inline function in ideal_c2p_mhd.hpp file)
    HydPrim1D w;
    bool dfloor_used=false, efloor_used=false, tfloor_used=false, bfloor_used=false;
    SingleC2P_CGLMHD(u, eos, w, dfloor_used, efloor_used, tfloor_used, bfloor_used);

    // set FOFC flag and quit loop if this function called only to check floors
    if (only_testfloors) {
      if (dfloor_used || (!eos.passive && (efloor_used || tfloor_used || nonfinite_input ||
          !Kokkos::isfinite(w.e) || !Kokkos::isfinite(w.pp) ||
          !(w.e >= eos.pfloor) || !(w.pp >= eos.pfloor)))) {
        fofc_(m,k,j,i) = true;
        if (eos.passive && dfloor_used) cons(m,IDN,k,j,i) = u.d;
        sumd++;  // use dfloor as counter for when either is true
      }
    } else {
      // update counter, reset conserved if floor was hit
      if (dfloor_used) {
        cons(m,IDN,k,j,i) = u.d;
        if (eos.passive) cons(m,IEN,k,j,i) = u.e;
        cons(m,IAN,k,j,i) = u.mu;
        sumd++;
      }
      if (efloor_used) {
        cons(m,IEN,k,j,i) = u.e;
        cons(m,IAN,k,j,i) = u.mu;
        sume++;
      }
      if (bfloor_used) {
        if (eos.passive) cons(m,IEN,k,j,i) = u.e;
        cons(m,IAN,k,j,i) = u.mu;
      }

      //if (tfloor_used) {
      //  cons(m,IEN,k,j,i) = u.e;
      //  sumt++;
      //}

      // store primitive state in 3D array
      prim(m,IDN,k,j,i) = w.d;
      prim(m,IVX,k,j,i) = w.vx;
      prim(m,IVY,k,j,i) = w.vy;
      prim(m,IVZ,k,j,i) = w.vz;
      prim(m,IEN,k,j,i) = w.e;
      prim(m,IPP,k,j,i) = w.pp;
      // store cell-centered fields in 3D array
      bcc(m,IBX,k,j,i) = u.bx;
      bcc(m,IBY,k,j,i) = u.by;
      bcc(m,IBZ,k,j,i) = u.bz;
      // convert scalars (if any), always stored at end of cons and prim arrays.
      for (int n=nmhd; n<(nmhd+nscal); ++n) {
        // apply scalar floor
        if (cons(m,n,k,j,i) < 0.0) {
          cons(m,n,k,j,i) = 0.0;
        }
        prim(m,n,k,j,i) = cons(m,n,k,j,i)/u.d;
      }
    }
  }, Kokkos::Sum<int>(nfloord_), Kokkos::Sum<int>(nfloore_), Kokkos::Sum<int>(nfloort_));

  // store appropriate counters
  if (only_testfloors) {
    pmy_pack->pmesh->ecounter.nfofc += nfloord_;
  } else {
    pmy_pack->pmesh->ecounter.neos_dfloor += nfloord_;
    pmy_pack->pmesh->ecounter.neos_efloor += nfloore_;
    pmy_pack->pmesh->ecounter.neos_tfloor += nfloort_;
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \!fn void CGLMagneticMomentToPrim()
//! \brief Converts conserved variables to primitives while IAN stores p_perp/|B| during
//! a CGL Landau-fluid STS sweep.

void CGLMHD::CGLMagneticMomentToPrim(DvceArray5D<Real> &cons,
                                     const DvceFaceFld4D<Real> &b,
                                     DvceArray5D<Real> &prim,
                                     DvceArray5D<Real> &bcc,
                                     const int il, const int iu,
                                     const int jl, const int ju,
                                     const int kl, const int ku) {
  const int nmhd  = pmy_pack->pmhd->nmhd;
  const int nscal = pmy_pack->pmhd->nscalars;
  const int nmb = pmy_pack->nmb_thispack;
  const EOS_Data eos = eos_data;

  const int ni   = (iu - il + 1);
  const int nji  = (ju - jl + 1)*ni;
  const int nkji = (ku - kl + 1)*nji;
  const int nmkji = nmb*nkji;

  int nfloord_=0, nfloore_=0, nfloort_=0;
  Kokkos::parallel_reduce("mhd_cgl_mub_to_prim",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, int &sumd, int &sume, int &sumt) {
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/ni;
    int i = (idx - m*nkji - k*nji - j*ni) + il;
    j += jl;
    k += kl;

    MHDCons1D u;
    u.d  = cons(m,IDN,k,j,i);
    u.mx = cons(m,IM1,k,j,i);
    u.my = cons(m,IM2,k,j,i);
    u.mz = cons(m,IM3,k,j,i);
    u.e  = cons(m,IEN,k,j,i);
    u.mu = cons(m,IAN,k,j,i);
    u.bx = 0.5*(b.x1f(m,k,j,i) + b.x1f(m,k,j,i+1));
    u.by = 0.5*(b.x2f(m,k,j,i) + b.x2f(m,k,j+1,i));
    u.bz = 0.5*(b.x3f(m,k,j,i) + b.x3f(m,k+1,j,i));

    HydPrim1D w;
    bool dfloor_used=false, efloor_used=false, tfloor_used=false, bfloor_used=false;
    SingleC2P_CGLMHDFromMagneticMoment(u, eos, w, dfloor_used, efloor_used,
                                       tfloor_used, bfloor_used);

    if (dfloor_used) {
      cons(m,IDN,k,j,i) = u.d;
      sumd++;
    }
    if (efloor_used) {
      cons(m,IEN,k,j,i) = u.e;
      cons(m,IAN,k,j,i) = u.mu;
      sume++;
    }
    if (bfloor_used) {
      cons(m,IAN,k,j,i) = u.mu;
    }

    prim(m,IDN,k,j,i) = w.d;
    prim(m,IVX,k,j,i) = w.vx;
    prim(m,IVY,k,j,i) = w.vy;
    prim(m,IVZ,k,j,i) = w.vz;
    prim(m,IEN,k,j,i) = w.e;
    prim(m,IPP,k,j,i) = w.pp;
    bcc(m,IBX,k,j,i) = u.bx;
    bcc(m,IBY,k,j,i) = u.by;
    bcc(m,IBZ,k,j,i) = u.bz;

    for (int n=nmhd; n<(nmhd+nscal); ++n) {
      if (cons(m,n,k,j,i) < 0.0) {
        cons(m,n,k,j,i) = 0.0;
      }
      prim(m,n,k,j,i) = cons(m,n,k,j,i)/u.d;
    }
  }, Kokkos::Sum<int>(nfloord_), Kokkos::Sum<int>(nfloore_), Kokkos::Sum<int>(nfloort_));

  pmy_pack->pmesh->ecounter.neos_dfloor += nfloord_;
  pmy_pack->pmesh->ecounter.neos_efloor += nfloore_;
  pmy_pack->pmesh->ecounter.neos_tfloor += nfloort_;

  return;
}

//----------------------------------------------------------------------------------------
//! \!fn void CGLRefreshPrimFromMagneticMoment()
//! \brief Lightweight primitive refresh for CGL LF STS stages.
//!
//! During a pure CGL Landau-fluid STS sweep, only total energy and the temporary
//! magnetic-moment slot change.  The face-centered and cell-centered magnetic fields are
//! fixed, so this refresh reuses bcc and avoids the full face-to-cell magnetic-field
//! average and scalar primitive conversion done by CGLMagneticMomentToPrim().

void CGLMHD::CGLRefreshPrimFromMagneticMoment(DvceArray5D<Real> &cons,
                                              const DvceArray5D<Real> &bcc,
                                              DvceArray5D<Real> &prim,
                                              const int il, const int iu,
                                              const int jl, const int ju,
                                              const int kl, const int ku) {
  const int nmb = pmy_pack->nmb_thispack;
  const EOS_Data eos = eos_data;

  const int ni   = (iu - il + 1);
  const int nji  = (ju - jl + 1)*ni;
  const int nkji = (ku - kl + 1)*nji;
  const int nmkji = nmb*nkji;

  int nfloord_=0, nfloore_=0, nfloort_=0;
  Kokkos::parallel_reduce("mhd_cgl_lf_refresh_prim",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, int &sumd, int &sume, int &sumt) {
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/ni;
    int i = (idx - m*nkji - k*nji - j*ni) + il;
    j += jl;
    k += kl;

    MHDCons1D u;
    u.d  = cons(m,IDN,k,j,i);
    u.mx = cons(m,IM1,k,j,i);
    u.my = cons(m,IM2,k,j,i);
    u.mz = cons(m,IM3,k,j,i);
    u.e  = cons(m,IEN,k,j,i);
    u.mu = cons(m,IAN,k,j,i);
    u.bx = bcc(m,IBX,k,j,i);
    u.by = bcc(m,IBY,k,j,i);
    u.bz = bcc(m,IBZ,k,j,i);

    HydPrim1D w;
    bool dfloor_used=false, efloor_used=false, tfloor_used=false, bfloor_used=false;
    SingleC2P_CGLMHDFromMagneticMoment(u, eos, w, dfloor_used, efloor_used,
                                       tfloor_used, bfloor_used);

    if (dfloor_used) {
      cons(m,IDN,k,j,i) = u.d;
      sumd++;
    }
    if (efloor_used) {
      cons(m,IEN,k,j,i) = u.e;
      cons(m,IAN,k,j,i) = u.mu;
      sume++;
    }
    if (bfloor_used) {
      cons(m,IAN,k,j,i) = u.mu;
    }
    prim(m,IDN,k,j,i) = w.d;
    prim(m,IVX,k,j,i) = w.vx;
    prim(m,IVY,k,j,i) = w.vy;
    prim(m,IVZ,k,j,i) = w.vz;
    prim(m,IEN,k,j,i) = w.e;
    prim(m,IPP,k,j,i) = w.pp;

    (void) sumt;
  }, Kokkos::Sum<int>(nfloord_), Kokkos::Sum<int>(nfloore_), Kokkos::Sum<int>(nfloort_));

  pmy_pack->pmesh->ecounter.neos_dfloor += nfloord_;
  pmy_pack->pmesh->ecounter.neos_efloor += nfloore_;
  pmy_pack->pmesh->ecounter.neos_tfloor += nfloort_;

  return;
}

//----------------------------------------------------------------------------------------
//! \!fn void PrimToCons()
//! \brief Converts conserved into primitive variables.  Operates over range of cells
//! given in argument list.  Does not change cell- or face-centered magnetic fields.

void CGLMHD::PrimToCons(const DvceArray5D<Real> &prim, const DvceArray5D<Real> &bcc,
                          DvceArray5D<Real> &cons, const int il, const int iu,
                          const int jl, const int ju, const int kl, const int ku) {
  const int nmhd  = pmy_pack->pmhd->nmhd;
  const int nscal = pmy_pack->pmhd->nscalars;
  const int nmb = pmy_pack->nmb_thispack;
  const Real bfloor = eos_data.bfloor;
  const bool passive = eos_data.passive;
  const bool magnetic_moment = pmy_pack->pmhd->cgl_slot_representation ==
      mhd::CGLSlotRepresentation::magnetic_moment;

  par_for("mhd_p2c", DevExeSpace(), 0, (nmb-1), kl, ku, jl, ju, il, iu,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    // load single state primitive variables
    MHDPrim1D w;
    w.d  = prim(m,IDN,k,j,i);
    w.vx = prim(m,IVX,k,j,i);
    w.vy = prim(m,IVY,k,j,i);
    w.vz = prim(m,IVZ,k,j,i);
    w.e  = prim(m,IEN,k,j,i);
    w.pp = prim(m,IPP,k,j,i);

    // load cell-centered fields into primitive state
    w.bx = bcc(m,IBX,k,j,i);
    w.by = bcc(m,IBY,k,j,i);
    w.bz = bcc(m,IBZ,k,j,i);

    // call p2c function
    HydCons1D u;
    if (passive) cgl::PassiveP2C(w, bfloor, u);
    else SingleP2C_CGLMHD(w, bfloor, u);
    // Same representation contract for boundary callbacks and temporary targets.
    if (magnetic_moment) {
      const Real bmag = sqrt(SQR(w.bx) + SQR(w.by) + SQR(w.bz));
      const Real pt = (bmag > bfloor) ? w.pp : TWO_3RDS*w.pp + ONE_3RD*w.e;
      if (passive) u.e = 0.5*w.e + w.pp;
      u.mu = pt/fmax(bmag, bfloor);
    }

    //no need to change pressures here if bfloor was hit as they'll be changed elsewhere

    // store conserved state in 3D array
    cons(m,IDN,k,j,i) = u.d;
    cons(m,IM1,k,j,i) = u.mx;
    cons(m,IM2,k,j,i) = u.my;
    cons(m,IM3,k,j,i) = u.mz;
    cons(m,IEN,k,j,i) = u.e;
    cons(m,IAN,k,j,i) = u.mu;

    // convert scalars (if any), always stored at end of cons and prim arrays.
    for (int n=nmhd; n<(nmhd+nscal); ++n) {
      cons(m,n,k,j,i) = u.d*prim(m,n,k,j,i);
    }
  });

  return;
}

//----------------------------------------------------------------------------------------
//! \!fn void CGLAnisotropyToMagneticMoment()
//! \brief Temporarily reinterpret the CGL IAN slot as mu_B = p_perp/|B| for LF STS.

void CGLMHD::CGLAnisotropyToMagneticMoment(DvceArray5D<Real> &cons,
                                            const DvceArray5D<Real> &bcc,
                                            const int il, const int iu,
                                            const int jl, const int ju,
                                            const int kl, const int ku) {
  const int nmb = pmy_pack->nmb_thispack;
  const Real pfloor = eos_data.pfloor;
  const Real bfloor = eos_data.bfloor;
  const bool passive = eos_data.passive;

  par_for("mhd_cgl_anis_to_mub", DevExeSpace(), 0, (nmb-1), kl, ku, jl, ju, il, iu,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (passive) {
      const Real bmag = sqrt(SQR(bcc(m,IBX,k,j,i)) + SQR(bcc(m,IBY,k,j,i)) +
                             SQR(bcc(m,IBZ,k,j,i)));
      Real pp, pt;
      cgl::PassiveDecode(cons(m,IDN,k,j,i),
                         {cons(m,IEN,k,j,i), cons(m,IAN,k,j,i)},
                         fmax(bmag, bfloor), pp, pt);
      cons(m,IEN,k,j,i) = 0.5*pp + pt;
      cons(m,IAN,k,j,i) = pt/fmax(bmag, bfloor);
      return;
    }
    cons(m,IAN,k,j,i) =
        CGLConservedAnisotropyToMagneticMoment(cons(m,IDN,k,j,i),
                                               cons(m,IM1,k,j,i),
                                               cons(m,IM2,k,j,i),
                                               cons(m,IM3,k,j,i),
                                               cons(m,IEN,k,j,i),
                                               cons(m,IAN,k,j,i),
                                               bcc(m,IBX,k,j,i),
                                               bcc(m,IBY,k,j,i),
                                               bcc(m,IBZ,k,j,i),
                                               bfloor, pfloor);
  });

  return;
}

//----------------------------------------------------------------------------------------
//! \!fn void CGLMagneticMomentToAnisotropy()
//! \brief Convert the temporary LF STS mu_B slot back to conserved anisotropy A.

void CGLMHD::CGLMagneticMomentToAnisotropy(DvceArray5D<Real> &cons,
                                            const DvceArray5D<Real> &bcc,
                                            const int il, const int iu,
                                            const int jl, const int ju,
                                            const int kl, const int ku) {
  const int nmb = pmy_pack->nmb_thispack;
  const Real pfloor = eos_data.pfloor;
  const Real bfloor = eos_data.bfloor;
  const bool passive = eos_data.passive;

  par_for("mhd_cgl_mub_to_anis", DevExeSpace(), 0, (nmb-1), kl, ku, jl, ju, il, iu,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (passive) {
      const Real bmag = sqrt(SQR(bcc(m,IBX,k,j,i)) + SQR(bcc(m,IBY,k,j,i)) +
                             SQR(bcc(m,IBZ,k,j,i)));
      const Real beff = fmax(bmag, bfloor);
      Real pt = cons(m,IAN,k,j,i)*beff;
      Real pp = 2.0*(cons(m,IEN,k,j,i) - pt);
      if (bmag <= bfloor) { pp = TWO_3RDS*cons(m,IEN,k,j,i); pt = pp; }
      const auto q = cgl::PassiveEncodeFloored(cons(m,IDN,k,j,i),
                                              fmax(pp, pfloor), fmax(pt, pfloor),
                                              beff, pfloor);
      cons(m,IEN,k,j,i) = q.j;
      cons(m,IAN,k,j,i) = q.a;
      return;
    }
    Real total_energy = cons(m,IEN,k,j,i);
    Real anisotropy = 0.0;
    CGLMagneticMomentToConservedAnisotropy(cons(m,IDN,k,j,i),
                                           cons(m,IM1,k,j,i),
                                           cons(m,IM2,k,j,i),
                                           cons(m,IM3,k,j,i),
                                           bcc(m,IBX,k,j,i),
                                           bcc(m,IBY,k,j,i),
                                           bcc(m,IBZ,k,j,i),
                                           bfloor, pfloor,
                                           cons(m,IAN,k,j,i),
                                           total_energy, anisotropy);
    cons(m,IEN,k,j,i) = total_energy;
    cons(m,IAN,k,j,i) = anisotropy;
  });

  return;
}

//----------------------------------------------------------------------------------------
//! \!fn void Collisions()
//! \brief Decays pressure anisotropy according to scattering rate. Operates over
//! range of cells given in argument list.

void CGLMHD::Collisions(DvceArray5D<Real> &prim, const DvceArray5D<Real> &bcc,
                          DvceArray5D<Real> &cons, const Real dtc,
                          const CGLCollisionMode mode,
                          const int il, const int iu,
                          const int jl, const int ju, const int kl, const int ku) {
  const int nmhd  = pmy_pack->pmhd->nmhd;
  const int nscal = pmy_pack->pmhd->nscalars;
  const int nmb = pmy_pack->nmb_thispack;
  const EOS_Data eos = eos_data;
  const bool flim = eos_data.flim;
  const bool mlim = eos_data.mlim;
  const Real bfloor = eos_data.bfloor;
  const auto *pcgl_lf = pmy_pack->pmhd->pcgl_lf;
  const bool landau_fluid_active = (pcgl_lf != nullptr);
  const bool backup = landau_fluid_active
                          ? pcgl_lf->effective_backup_limiter
                          : eos_data.backup_lim;

  // TODO(cgl-lf): If the limiter/collision closure becomes nonlocal or
  // gradient-dependent, store nu_eff on the grid here and reconstruct it to LF
  // faces instead of recomputing the present algebraic thresholds in conduction.
  // The update ordering will then need a clear pre/post heat-flux split.

  par_for("mhd_coll", DevExeSpace(), 0, (nmb-1), kl, ku, jl, ju, il, iu,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    // load single state primitive variables
    MHDPrim1D w;
    w.d  = prim(m,IDN,k,j,i);
    w.vx = prim(m,IVX,k,j,i);
    w.vy = prim(m,IVY,k,j,i);
    w.vz = prim(m,IVZ,k,j,i);
    w.e  = prim(m,IPR,k,j,i);
    w.pp = prim(m,IPP,k,j,i);

    // load cell-centered fields into primitive state
    w.bx = bcc(m,IBX,k,j,i);
    w.by = bcc(m,IBY,k,j,i);
    w.bz = bcc(m,IBZ,k,j,i);

    // Leave conserved A untouched when neither rates nor walls change pressures.
    const Real initial_ppar = w.e;
    const Real initial_pperp = w.pp;
    HydCons1D u;
    if (mode == CGLCollisionMode::full) {
      SingleCollRates_CGLMHD(w, eos, dtc, mlim, flim);
    }
    SingleCollWalls_CGLMHD(w, eos, backup);
    if (w.e == initial_ppar && w.pp == initial_pperp) return;
    if (eos.passive) {
      const Real eint = 0.5*initial_ppar + initial_pperp;
      const Real bsqr = SQR(w.bx) + SQR(w.by) + SQR(w.bz);
      const auto q = cgl::PassiveWallEncode(w.d, eint, w.pp - w.e, bsqr,
                                            eos, backup, w.e, w.pp);
      cons(m,IEN,k,j,i) = q.j;
      cons(m,IAN,k,j,i) = q.a;
      prim(m,IPR,k,j,i) = w.e;
      prim(m,IPP,k,j,i) = w.pp;
      return;
    }
    SingleP2C_CGLMHD(w, bfloor, u);
    const Real kinetic = 0.5*(SQR(cons(m,IM1,k,j,i)) + SQR(cons(m,IM2,k,j,i)) +
                              SQR(cons(m,IM3,k,j,i)))/w.d;
    const Real magnetic = 0.5*(SQR(w.bx) + SQR(w.by) + SQR(w.bz));
    const Real eint = cons(m,IEN,k,j,i) - kinetic - magnetic;
    const Real admissible = CGLWallAdmissibleAnisotropy(w, eint, u.mu, eos, backup);
    if (admissible != u.mu) {
      u.mu = admissible;
      CGLRecoverPressuresFromInternalEnergyAndAnisotropy(
          w.d, eint, u.mu, fmax(sqrt(2.0*magnetic), bfloor), w.e, w.pp);
    }

    // Correct conserved anisotropy variable
    cons(m,IAN,k,j,i) = u.mu;

    // Correct pressures
    prim(m,IPR,k,j,i) = w.e;
    prim(m,IPP,k,j,i) = w.pp;
  });

  return;
}
