#ifndef EOS_CGL_PASSIVE_HPP_
#define EOS_CGL_PASSIVE_HPP_

#include <limits>

#include "eos/isothermal_c2p_mhd.hpp"
#include "eos/cgl_physics.hpp"

namespace cgl {

// Passive thermodynamics carries two material invariants. IEN stores J and
// IAN stores A; neither slot contains the energy of the isothermal flow.
struct PassivePair { Real j, a; };

KOKKOS_INLINE_FUNCTION
PassivePair PassiveSpecificInvariants(const Real rho, const Real ppar,
                                      const Real pperp, const Real bmag) {
  const Real logr = log(rho), logb = log(bmag);
  return {log(ppar) + 2.0*logb - 3.0*logr,
          log(pperp) - log(ppar) + 2.0*logr - 3.0*logb};
}

KOKKOS_INLINE_FUNCTION
PassivePair PassiveEncode(const Real rho, const Real ppar, const Real pperp,
                           const Real bmag) {
  const PassivePair q = PassiveSpecificInvariants(rho, ppar, pperp, bmag);
  return {rho*q.j, rho*q.a};
}

KOKKOS_INLINE_FUNCTION
void PassiveDecode(const Real rho, const PassivePair q, const Real bmag,
                    Real &ppar, Real &pperp) {
  const Real logr = log(rho), logb = log(bmag);
  const Real j = q.j/rho;
  ppar = exp(j + 3.0*logr - 2.0*logb);
  // Preserve exact isotropy even when B is tiny, so roundoff cannot create a
  // spurious fluid-firehose violation at an isotropic low-field state.
  if (q.a == rho*(2.0*logr - 3.0*logb)) {
    pperp = ppar;
  } else {
    pperp = exp(q.a/rho + j + logr + logb);
  }
}

// The two exponentials must recover pressures above their floor. Raising the
// common J coordinate preserves the encoded pressure ratio. This correction
// is used only after an actual repair, never on an unchanged admissible state.
KOKKOS_INLINE_FUNCTION
PassivePair PassiveEncodeFloored(const Real rho, const Real ppar, const Real pperp,
                                 const Real bmag, const Real pfloor) {
  PassivePair q = PassiveEncode(rho, ppar, pperp, bmag);
  if (!Kokkos::isfinite(q.j) || !Kokkos::isfinite(q.a)) {
    Kokkos::abort("Passive CGL invariant encoding is not representable in Real");
  }
  Real pp, pt;
  PassiveDecode(rho, q, bmag, pp, pt);
  if (pp < pfloor || pt < pfloor) {
    const Real shift = rho*fmax(log(pfloor/pp), log(pfloor/pt));
    q.j = Kokkos::nextafter(q.j + shift, std::numeric_limits<Real>::infinity());
    PassiveDecode(rho, q, bmag, pp, pt);
    int tries = 0;
    while (pp < pfloor || pt < pfloor) {
      if (++tries > 32 || !Kokkos::isfinite(q.j) ||
          !Kokkos::isfinite(pp) || !Kokkos::isfinite(pt)) {
        Kokkos::abort("Passive CGL pressure floor cannot be encoded in Real");
      }
      q.j = Kokkos::nextafter(q.j, std::numeric_limits<Real>::infinity());
      PassiveDecode(rho, q, bmag, pp, pt);
    }
  }
  if (!Kokkos::isfinite(q.j) || !Kokkos::isfinite(q.a) ||
      !Kokkos::isfinite(pp) || !Kokkos::isfinite(pt)) {
    Kokkos::abort("Passive CGL pressure floor is not representable in Real");
  }
  return q;
}

KOKKOS_INLINE_FUNCTION
void PassiveC2P(MHDCons1D &u, const EOS_Data &eos, HydPrim1D &w,
                 bool &dfloor_used, bool &efloor_used, bool &bfloor_used) {
  const Real bsqr = SQR(u.bx) + SQR(u.by) + SQR(u.bz);
  const Real bmag = sqrt(bsqr), beff = fmax(bmag, eos.bfloor);
  const Real original_rho = u.d;
  SingleC2P_IsothermalMHD(u, fmax(eos.dfloor, bsqr/eos.sigma_max), w, dfloor_used);
  PassiveDecode(original_rho, {u.e, u.mu}, beff, w.e, w.pp);
  if (!Kokkos::isfinite(w.e) || !(w.e >= eos.pfloor)) {
    w.e = eos.pfloor;
    efloor_used = true;
  }
  if (!Kokkos::isfinite(w.pp) || !(w.pp >= eos.pfloor)) {
    w.pp = eos.pfloor;
    efloor_used = true;
  }
  if (bmag <= eos.bfloor && w.e != w.pp) {
    w.e = TWO_3RDS*w.pp + ONE_3RD*w.e;
    w.pp = w.e;
    bfloor_used = true;
  }
  if (dfloor_used || efloor_used || bfloor_used) {
    const PassivePair q = PassiveEncodeFloored(w.d, w.e, w.pp, beff, eos.pfloor);
    u.e = q.j;
    u.mu = q.a;
    PassiveDecode(w.d, q, beff, w.e, w.pp);
  }
}

KOKKOS_INLINE_FUNCTION
void PassiveP2C(const MHDPrim1D &w, const Real bfloor, HydCons1D &u) {
  // Keep precisely the multiplication sequence of the reference isothermal EOS.
  HydPrim1D flow;
  flow.d = w.d; flow.vx = w.vx; flow.vy = w.vy; flow.vz = w.vz;
  SingleP2C_IsothermalMHD(flow, u);
  const Real bmag = sqrt(SQR(w.bx) + SQR(w.by) + SQR(w.bz));
  Real pp = w.e, pt = w.pp;
  if (bmag <= bfloor) {
    pp = TWO_3RDS*pt + ONE_3RD*pp;
    pt = pp;
  }
  const PassivePair q = PassiveEncode(w.d, pp, pt, fmax(bmag, bfloor));
  if (!Kokkos::isfinite(q.j) || !Kokkos::isfinite(q.a)) {
    Kokkos::abort("Passive CGL initial invariants are not representable in Real");
  }
  u.e = q.j;
  u.mu = q.a;
}

// LF evolves actual U and mu_B temporarily. This avoids subtracting kinetic
// and magnetic energies from a passive thermal quantity at high Mach number.
KOKKOS_INLINE_FUNCTION
void PassiveMomentC2P(MHDCons1D &u, const EOS_Data &eos, HydPrim1D &w,
                       bool &dfloor_used, bool &efloor_used, bool &bfloor_used,
                       const Real bmag) {
  const Real bsqr = SQR(u.bx) + SQR(u.by) + SQR(u.bz);
  SingleC2P_IsothermalMHD(u, fmax(eos.dfloor, bsqr/eos.sigma_max), w, dfloor_used);
  const Real beff = fmax(bmag, eos.bfloor);
  w.pp = u.mu*beff;
  w.e = 2.0*(u.e - w.pp);
  if (bmag <= eos.bfloor) {
    w.e = TWO_3RDS*u.e;
    w.pp = w.e;
    u.mu = w.pp/beff;
    bfloor_used = true;
  }
  if (!Kokkos::isfinite(w.e) || !(w.e >= eos.pfloor)) {
    w.e = eos.pfloor;
    efloor_used = true;
  }
  if (!Kokkos::isfinite(w.pp) || !(w.pp >= eos.pfloor)) {
    w.pp = eos.pfloor;
    efloor_used = true;
  }
  if (efloor_used) {
    u.e = 0.5*w.e + w.pp;
    u.mu = w.pp/beff;
  }
}

// Re-encoding a boundary pressure can round outside the hard wall. Bisect the
// physical anisotropy toward exact isotropy at fixed thermal U until both
// conserved coordinates decode inside. Rates/walls call this only after a change.
KOKKOS_INLINE_FUNCTION
PassivePair PassiveWallEncode(const Real rho, const Real eint, const Real delta,
                               const Real bsqr, const EOS_Data &eos,
                               const bool backup, Real &ppar, Real &pperp) {
  const Real bmag = fmax(sqrt(bsqr), eos.bfloor);
  const Real piso = TWO_3RDS*eint;
  Real trial = delta, outside = delta, inside = 0.0;
  PassivePair q;
  PassivePair best = PassiveEncodeFloored(rho, fmax(piso, eos.pfloor),
                                           fmax(piso, eos.pfloor), bmag, eos.pfloor);
  bool first = true;
  while (true) {
    const Real pp = piso - TWO_3RDS*trial;
    const Real pt = piso + ONE_3RD*trial;
    q = PassiveEncodeFloored(rho, fmax(pp, eos.pfloor), fmax(pt, eos.pfloor),
                             bmag, eos.pfloor);
    PassiveDecode(rho, q, bmag, ppar, pperp);
    if (!HardBoundViolated(pperp - ppar, bsqr, eos, backup)) {
      if (first) return q;
      inside = trial;
      best = q;
    } else {
      outside = trial;
    }
    first = false;
    trial = 0.5*inside + 0.5*outside;
    if (trial == inside || trial == outside) {
      PassiveDecode(rho, best, bmag, ppar, pperp);
      return best;
    }
  }
}

}  // namespace cgl
#endif  // EOS_CGL_PASSIVE_HPP_
