#ifndef EOS_CGL_AMR_PROJECTION_HPP_
#define EOS_CGL_AMR_PROJECTION_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_amr_projection.hpp
//! \brief Single-cell CGL AMR thermodynamic projection helpers.

#include <limits>

#include "athena.hpp"
#include "eos/eos.hpp"
#include "eos/ideal_c2p_mhd.hpp"

namespace cgl {
namespace amr {

enum SlotRepresentation {
  anisotropy,
  magnetic_moment
};

using RepairMask = unsigned int;

enum Repair : RepairMask {
  kNone                   = 0u,
  kNonfiniteThermo        = 1u << 0,
  kDensityFloor           = 1u << 1,
  kInternalEnergyFloor    = 1u << 2,
  kParallelPressureFloor  = 1u << 3,
  kPerpPressureFloor      = 1u << 4,
  kLowFieldIsotropized    = 1u << 5,
  kFirehoseHardwall       = 1u << 6,
  kMirrorHardwall         = 1u << 7,
  kAnisotropyChanged      = 1u << 8,
  kIntervalEnergyExpanded = 1u << 9,
  kSlopeScaled            = 1u << 10
};

struct ProjectionReport {
  RepairMask repairs = kNone;
  Real U = 0.0;
  Real delta = 0.0;
  Real p_parallel = 0.0;
  Real p_perp = 0.0;
  Real bmag = 0.0;
  Real b_eff = 0.0;
};

KOKKOS_INLINE_FUNCTION
Real InternalEnergyFromPressures(const Real p_parallel, const Real p_perp) {
  return p_perp + 0.5*p_parallel;
}

KOKKOS_INLINE_FUNCTION
Real DeltaFromPressures(const Real p_parallel, const Real p_perp) {
  return p_perp - p_parallel;
}

KOKKOS_INLINE_FUNCTION
void PressuresFromUDelta(const Real U, const Real delta,
                         Real &p_parallel, Real &p_perp) {
  p_parallel = TWO_3RDS*(U - delta);
  p_perp = TWO_3RDS*U + ONE_3RD*delta;
}

KOKKOS_INLINE_FUNCTION
bool Finite3(const Real a, const Real b, const Real c) {
  return Kokkos::isfinite(a) && Kokkos::isfinite(b) && Kokkos::isfinite(c);
}

KOKKOS_INLINE_FUNCTION
void DeltaInterval(const Real U, const Real bsqr, const Real bmag,
                   const EOS_Data &eos, Real &delta_min, Real &delta_max) {
  delta_min = 3.0*eos.pfloor - 2.0*U;
  delta_max = U - 1.5*eos.pfloor;
  if (eos.hardwall_lim && bmag > eos.bfloor) {
    if (eos.flim) {
      delta_min = fmax(delta_min, eos.firehose_threshold*bsqr);
    }
    if (eos.mlim) {
      delta_max = fmin(delta_max, cgl::kMirrorThreshold*bsqr);
    }
  }
}

KOKKOS_INLINE_FUNCTION
bool IsAdmissiblePrimaryState(const Real rho, const Real vx, const Real vy,
                              const Real vz, const Real bx, const Real by,
                              const Real bz, const Real U,
                              const EOS_Data &eos) {
  if (!Finite3(rho, vx, vy) || !Finite3(vz, bx, by) ||
      !Kokkos::isfinite(bz) || !Kokkos::isfinite(U)) {
    return false;
  }
  if (!(rho >= eos.dfloor) || !(U >= 1.5*eos.pfloor)) {
    return false;
  }
  const Real bsqr = SQR(bx) + SQR(by) + SQR(bz);
  const Real bmag = sqrt(bsqr);
  if (!(Kokkos::isfinite(bsqr)) || !(Kokkos::isfinite(bmag))) {
    return false;
  }
  if (bmag <= eos.bfloor) return true;

  Real delta_min, delta_max;
  DeltaInterval(U, bsqr, bmag, eos, delta_min, delta_max);
  const Real tol = 16.0*std::numeric_limits<Real>::epsilon()*
                   fmax(fmax(fabs(delta_min), fabs(delta_max)), fmax(U, bsqr));
  return delta_min <= delta_max + tol;
}

KOKKOS_INLINE_FUNCTION
bool IsAdmissibleUDelta(const Real rho, const Real vx, const Real vy,
                        const Real vz, const Real bx, const Real by,
                        const Real bz, const Real U, const Real delta,
                        const EOS_Data &eos) {
  if (!Finite3(rho, vx, vy) || !Finite3(vz, bx, by) ||
      !Kokkos::isfinite(bz) || !Kokkos::isfinite(U) ||
      !Kokkos::isfinite(delta)) {
    return false;
  }
  if (!(rho >= eos.dfloor) || !(U >= 1.5*eos.pfloor)) {
    return false;
  }
  const Real bsqr = SQR(bx) + SQR(by) + SQR(bz);
  const Real bmag = sqrt(bsqr);
  if (!(Kokkos::isfinite(bsqr)) || !(Kokkos::isfinite(bmag))) {
    return false;
  }
  if (bmag <= eos.bfloor) {
    return fabs(delta) <= 16.0*std::numeric_limits<Real>::epsilon()*fmax(U, 1.0);
  }
  Real delta_min, delta_max;
  DeltaInterval(U, bsqr, bmag, eos, delta_min, delta_max);
  const Real tol = 16.0*std::numeric_limits<Real>::epsilon()*
                   fmax(fmax(fabs(delta_min), fabs(delta_max)), fmax(U, bsqr));
  return delta >= delta_min - tol && delta <= delta_max + tol;
}

KOKKOS_INLINE_FUNCTION
ProjectionReport ProjectUDeltaToCGL(const Real rho_in, const Real vx_in,
                                    const Real vy_in, const Real vz_in,
                                    const Real bx_in, const Real by_in,
                                    const Real bz_in, const Real U_in,
                                    const Real delta_in, const EOS_Data &eos,
                                    const SlotRepresentation slot,
                                    MHDPrim1D &w, HydCons1D &u) {
  ProjectionReport report;
  RepairMask repairs = kNone;

  const Real bx = bx_in;
  const Real by = by_in;
  const Real bz = bz_in;
  const Real bsqr = SQR(bx) + SQR(by) + SQR(bz);
  const Real bmag = sqrt(bsqr);
  if (!Finite3(bx, by, bz) || !Kokkos::isfinite(bsqr) ||
      !Kokkos::isfinite(bmag)) {
    Kokkos::abort(
        "CGL AMR projection requires a finite face-centered magnetic field");
  }

  Real rho = rho_in;
  if (!Kokkos::isfinite(rho)) {
    rho = eos.dfloor;
    repairs |= kNonfiniteThermo | kDensityFloor;
  } else if (rho < eos.dfloor) {
    rho = eos.dfloor;
    repairs |= kDensityFloor;
  }

  Real vx = vx_in;
  Real vy = vy_in;
  Real vz = vz_in;
  if (!Finite3(vx, vy, vz)) {
    if (!Kokkos::isfinite(vx)) vx = 0.0;
    if (!Kokkos::isfinite(vy)) vy = 0.0;
    if (!Kokkos::isfinite(vz)) vz = 0.0;
    repairs |= kNonfiniteThermo;
  }

  Real U = U_in;
  if (!Kokkos::isfinite(U)) {
    U = 1.5*eos.pfloor;
    repairs |= kNonfiniteThermo | kInternalEnergyFloor;
  } else if (U < 1.5*eos.pfloor) {
    U = 1.5*eos.pfloor;
    repairs |= kInternalEnergyFloor;
  }

  Real delta = delta_in;
  if (!Kokkos::isfinite(delta)) {
    delta = 0.0;
    repairs |= kNonfiniteThermo | kAnisotropyChanged;
  }

  const Real b_eff = (bmag > eos.bfloor) ? bmag : eos.bfloor;
  report.bmag = bmag;
  report.b_eff = b_eff;

  Real delta_min, delta_max;
  DeltaInterval(U, bsqr, bmag, eos, delta_min, delta_max);
  if (delta_min > delta_max) {
    const Real U_floor = 1.5*eos.pfloor;
    U = fmax(U, U_floor);
    DeltaInterval(U, bsqr, bmag, eos, delta_min, delta_max);
    if (delta_min > delta_max) {
      delta_min = delta_max = 0.0;
    }
    repairs |= kIntervalEnergyExpanded | kInternalEnergyFloor;
  }

  if (bmag <= eos.bfloor) {
    repairs |= kLowFieldIsotropized;
    if (delta != 0.0) {
      repairs |= kAnisotropyChanged;
    }
    delta = 0.0;
  } else if (delta < delta_min) {
    if (eos.hardwall_lim && eos.flim &&
        eos.firehose_threshold*bsqr >= delta_min &&
        delta < eos.firehose_threshold*bsqr) {
      repairs |= kFirehoseHardwall;
    }
    repairs |= kAnisotropyChanged;
    delta = delta_min;
  } else if (delta > delta_max) {
    if (eos.hardwall_lim && eos.mlim &&
        cgl::kMirrorThreshold*bsqr <= delta_max &&
        delta > cgl::kMirrorThreshold*bsqr) {
      repairs |= kMirrorHardwall;
    }
    repairs |= kAnisotropyChanged;
    delta = delta_max;
  }

  Real p_parallel, p_perp;
  PressuresFromUDelta(U, delta, p_parallel, p_perp);
  if (!Kokkos::isfinite(p_parallel) || p_parallel < eos.pfloor) {
    p_parallel = eos.pfloor;
    repairs |= kParallelPressureFloor;
  }
  if (!Kokkos::isfinite(p_perp) || p_perp < eos.pfloor) {
    p_perp = eos.pfloor;
    repairs |= kPerpPressureFloor;
  }
  U = InternalEnergyFromPressures(p_parallel, p_perp);
  delta = DeltaFromPressures(p_parallel, p_perp);

  w.d = rho;
  w.vx = vx;
  w.vy = vy;
  w.vz = vz;
  w.e = p_parallel;
  w.pp = p_perp;
  w.bx = bx;
  w.by = by;
  w.bz = bz;

  const Real ke = 0.5*rho*(SQR(vx) + SQR(vy) + SQR(vz));
  const Real me = 0.5*bsqr;
  u.d = rho;
  u.mx = rho*vx;
  u.my = rho*vy;
  u.mz = rho*vz;
  u.e = U + ke + me;
  u.mu = (slot == magnetic_moment)
             ? p_perp/b_eff
             : CGLConservedAnisotropy(rho, p_parallel, p_perp, b_eff);

  report.repairs = repairs;
  report.U = U;
  report.delta = delta;
  report.p_parallel = p_parallel;
  report.p_perp = p_perp;
  return report;
}

KOKKOS_INLINE_FUNCTION
ProjectionReport ProjectPrimitiveToCGL(const MHDPrim1D &w_in, const EOS_Data &eos,
                                       const SlotRepresentation slot,
                                       MHDPrim1D &w_out, HydCons1D &u_out) {
  const Real U = InternalEnergyFromPressures(w_in.e, w_in.pp);
  const Real delta = DeltaFromPressures(w_in.e, w_in.pp);
  return ProjectUDeltaToCGL(w_in.d, w_in.vx, w_in.vy, w_in.vz,
                            w_in.bx, w_in.by, w_in.bz, U, delta, eos,
                            slot, w_out, u_out);
}

KOKKOS_INLINE_FUNCTION
ProjectionReport ProjectConservedToCGL(const MHDCons1D &u_in, const EOS_Data &eos,
                                       const SlotRepresentation slot,
                                       MHDPrim1D &w_out, HydCons1D &u_out) {
  Real rho = u_in.d;
  if (!Kokkos::isfinite(rho) || rho < eos.dfloor) rho = eos.dfloor;
  const Real di = 1.0/rho;
  const Real vx = u_in.mx*di;
  const Real vy = u_in.my*di;
  const Real vz = u_in.mz*di;
  const Real bsqr = SQR(u_in.bx) + SQR(u_in.by) + SQR(u_in.bz);
  const Real bmag = sqrt(bsqr);
  const Real b_eff = (bmag > eos.bfloor) ? bmag : eos.bfloor;
  const Real ke = 0.5*di*(SQR(u_in.mx) + SQR(u_in.my) + SQR(u_in.mz));
  const Real me = 0.5*bsqr;
  const Real U = u_in.e - ke - me;

  Real p_parallel, p_perp;
  if (slot == magnetic_moment) {
    if (bmag > eos.bfloor) {
      p_perp = u_in.mu*b_eff;
      p_parallel = 2.0*(U - p_perp);
    } else {
      p_parallel = TWO_3RDS*U;
      p_perp = p_parallel;
    }
  } else if (bmag > eos.bfloor) {
    CGLRecoverPressuresFromInternalEnergyAndAnisotropy(rho, U, u_in.mu, b_eff,
                                                       p_parallel, p_perp);
  } else {
    p_parallel = TWO_3RDS*U;
    p_perp = p_parallel;
  }
  const Real delta = DeltaFromPressures(p_parallel, p_perp);
  return ProjectUDeltaToCGL(rho, vx, vy, vz, u_in.bx, u_in.by, u_in.bz, U,
                            delta, eos, slot, w_out, u_out);
}

KOKKOS_INLINE_FUNCTION
int MaskBit(const RepairMask mask, const Repair bit) {
  return (mask & bit) ? 1 : 0;
}

} // namespace amr
} // namespace cgl

#endif // EOS_CGL_AMR_PROJECTION_HPP_
