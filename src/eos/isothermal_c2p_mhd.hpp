#ifndef EOS_ISOTHERMAL_C2P_MHD_HPP_
#define EOS_ISOTHERMAL_C2P_MHD_HPP_

//----------------------------------------------------------------------------------------
//! \!fn void SingleC2P_IsothermalMHD()
//! \brief Converts conserved into primitive variables.  Operates over range of cells
//! given in argument list.

KOKKOS_INLINE_FUNCTION
void SingleC2P_IsothermalMHD(MHDCons1D &u, const Real &dfloor_,
                             HydPrim1D &w, bool &dfloor_used) {
  // apply density floor, without changing momentum
  if (u.d < dfloor_) {
    u.d = dfloor_;
    dfloor_used = true;
  }
  w.d = u.d;
  // compute velocities
  Real di = 1.0/u.d;
  w.vx = di*u.mx;
  w.vy = di*u.my;
  w.vz = di*u.mz;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void SingleP2C_IsothermalMHD()
//! \brief Converts single state of primitive variables into conserved variables for
//! non-relativistic MHD with an isothermal EOS.

KOKKOS_INLINE_FUNCTION
void SingleP2C_IsothermalMHD(const HydPrim1D &w, HydCons1D &u) {
  u.d  = w.d;
  u.mx = w.d*w.vx;
  u.my = w.d*w.vy;
  u.mz = w.d*w.vz;
  return;
}

#endif  // EOS_ISOTHERMAL_C2P_MHD_HPP_
