//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file eos.cpp
//! \brief implements constructor and some fns for EquationOfState abstract base class

#include <float.h>
#include <string>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "parameter_input.hpp"
#include "eos/eos.hpp"
#include "eos/cgl_physics.hpp"

//----------------------------------------------------------------------------------------
// EquationOfState constructor

EquationOfState::EquationOfState(std::string bk, MeshBlockPack* pp, ParameterInput *pin) :
    pmy_pack(pp) {
  eos_data.dfloor = pin->GetOrAddReal(bk,"dfloor",(FLT_MIN));
  eos_data.pfloor = pin->GetOrAddReal(bk,"pfloor",(FLT_MIN));
  eos_data.tfloor = pin->GetOrAddReal(bk,"tfloor",(FLT_MIN));
  eos_data.sfloor = pin->GetOrAddReal(bk,"sfloor",(FLT_MIN));
  eos_data.bfloor = pin->GetOrAddReal(bk,"bfloor",sqrt(1024.0*FLT_MIN));
  eos_data.is_cgl = false;
  eos_data.passive = false;
  eos_data.mlim = false;
  eos_data.flim = false;
  eos_data.coll = false;
  eos_data.backup_lim = false;
  eos_data.nu_coll = 0.0;
  eos_data.lim_coll = 0.0;
  eos_data.firehose_threshold = cgl::kFirehoseObliqueThreshold;
}

//----------------------------------------------------------------------------------------
//! \fn void ConsToPrim()
//! \brief No-Op versions of hydro and MHD conservative to primitive functions.
//! Required because each derived class overrides only one.

void EquationOfState::ConsToPrim(DvceArray5D<Real> &cons, DvceArray5D<Real> &prim,
                                 const bool only_testfloors,
                                 const int il, const int iu, const int jl, const int ju,
                                 const int kl, const int ku) {
  Kokkos::abort("NoOp hydro ConsToPrim called.\n"
                "  If using MHD, call MHD version instead.\n"
                "  If using DynGRMHD, use the functions exposed in DynGRMHD instead.\n");
}

void EquationOfState::ConsToPrim(DvceArray5D<Real> &cons, const DvceFaceFld4D<Real> &b,
                                 DvceArray5D<Real> &prim, DvceArray5D<Real> &bcc,
                                 const bool only_testfloors,
                                 const int il, const int iu, const int jl, const int ju,
                                 const int kl, const int ku) {
  Kokkos::abort("NoOp MHD ConsToPrim called.\n"
                "  If using hydro, call hydro version instead.\n"
                "  If using DynGRMHD, use the functions exposed in DynGRMHD instead.\n");
}

//----------------------------------------------------------------------------------------
//! \fn void PrimToCon()
//! \brief No-Op versions of hydro and MHD primitive to conservative functions.
//! Required because each derived class overrides only one.

void EquationOfState::PrimToCons(const DvceArray5D<Real> &prim, DvceArray5D<Real> &cons,
                                 const int il, const int iu, const int jl, const int ju,
                                 const int kl, const int ku) {
  Kokkos::abort("NoOp hydro PrimToCons called.\n"
                "  If using MHD, call MHD version instead.\n"
                "  If using DynGRMHD, use the functions exposed in DynGRMHD instead.\n");
}
void EquationOfState::PrimToCons(const DvceArray5D<Real> &prim,
                                 const DvceArray5D<Real> &bcc, DvceArray5D<Real> &cons,
                                 const int il, const int iu, const int jl, const int ju,
                                 const int kl, const int ku) {
  Kokkos::abort("NoOp MHD PrimToCons called.\n"
                "  If using hydro, call hydro version instead.\n"
                "  If using DynGRMHD, use the functions exposed in DynGRMHD instead.\n");
}

void EquationOfState::Collisions(DvceArray5D<Real> &prim, const DvceArray5D<Real> &bcc,
                                 DvceArray5D<Real> &cons, const Real dtc,
                                 const int il, const int iu,
                                 const int jl, const int ju, const int kl, const int ku) {
  (void) prim;
  (void) bcc;
  (void) cons;
  (void) dtc;
  (void) il;
  (void) iu;
  (void) jl;
  (void) ju;
  (void) kl;
  (void) ku;
}

void EquationOfState::CGLAnisotropyToMagneticMoment(DvceArray5D<Real> &cons,
                                                     const DvceArray5D<Real> &bcc,
                                                     const int il, const int iu,
                                                     const int jl, const int ju,
                                                     const int kl, const int ku) {
  (void) cons;
  (void) bcc;
  (void) il;
  (void) iu;
  (void) jl;
  (void) ju;
  (void) kl;
  (void) ku;
  Kokkos::abort("CGLAnisotropyToMagneticMoment is only defined for CGL MHD.\n");
}

void EquationOfState::CGLMagneticMomentToAnisotropy(DvceArray5D<Real> &cons,
                                                     const DvceArray5D<Real> &bcc,
                                                     const int il, const int iu,
                                                     const int jl, const int ju,
                                                     const int kl, const int ku) {
  (void) cons;
  (void) bcc;
  (void) il;
  (void) iu;
  (void) jl;
  (void) ju;
  (void) kl;
  (void) ku;
  Kokkos::abort("CGLMagneticMomentToAnisotropy is only defined for CGL MHD.\n");
}

void EquationOfState::CGLMagneticMomentToPrim(DvceArray5D<Real> &cons,
                                              const DvceFaceFld4D<Real> &b,
                                              DvceArray5D<Real> &prim,
                                              DvceArray5D<Real> &bcc,
                                              const int il, const int iu,
                                              const int jl, const int ju,
                                              const int kl, const int ku) {
  (void) cons;
  (void) b;
  (void) prim;
  (void) bcc;
  (void) il;
  (void) iu;
  (void) jl;
  (void) ju;
  (void) kl;
  (void) ku;
  Kokkos::abort("CGLMagneticMomentToPrim is only defined for CGL MHD.\n");
}

void EquationOfState::CGLRefreshPrimFromMagneticMoment(DvceArray5D<Real> &cons,
                                                       const DvceArray5D<Real> &bcc,
                                                       DvceArray5D<Real> &prim,
                                                       const int il, const int iu,
                                                       const int jl, const int ju,
                                                       const int kl, const int ku) {
  (void) cons;
  (void) bcc;
  (void) prim;
  (void) il;
  (void) iu;
  (void) jl;
  (void) ju;
  (void) kl;
  (void) ku;
  Kokkos::abort("CGLRefreshPrimFromMagneticMoment is only defined for CGL MHD.\n");
}
