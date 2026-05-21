//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cooling_test.cpp
//! \brief Uniform 1D problem generator for testing the general cooling source term.

#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "pgen/pgen.hpp"
#include "srcterms/cooling_hooks.hpp"

namespace {

struct CoolingTestUserParams {
  Real lambda = 0.0;
  Real gamma = 0.0;
  cooling::UnitSystem lambda_units = cooling::UnitSystem::code;
  cooling::UnitSystem gamma_units = cooling::UnitSystem::code;
};

CoolingTestUserParams user_params;

void CoolingTestFatal(const std::string &message) {
  std::cout << "### FATAL ERROR in cooling_test pgen" << std::endl
            << message << std::endl;
  std::exit(EXIT_FAILURE);
}

cooling::UnitSystem CoolingTestUnitSystem(const std::string &units,
                                          const std::string &name) {
  if (units == "code") return cooling::UnitSystem::code;
  if (units == "cgs") return cooling::UnitSystem::cgs;
  CoolingTestFatal(name + " must be 'code' or 'cgs'.");
  return cooling::UnitSystem::code;
}

struct CoolingTestUserEvaluator {
  Real lambda;
  Real gamma;
  cooling::UnitSystem lambda_units;
  cooling::UnitSystem gamma_units;

  KOKKOS_INLINE_FUNCTION
  cooling::CoolingRates operator()(const cooling::CoolingCellState&) const {
    cooling::CoolingRates rates;
    rates.lambda = lambda;
    rates.gamma = gamma;
    rates.lambda_units = lambda_units;
    rates.gamma_units = gamma_units;
    return rates;
  }
};

void CoolingTestUserCooling(MeshBlockPack *pmbp, const DvceArray5D<Real> &w0,
                            const EOS_Data &eos_data,
                            const cooling::RuntimeData &runtime, const Real bdt,
                            const Real history_bdt, DvceArray5D<Real> &u0,
                            Real &gross_energy, Real &net_energy) {
  CoolingTestUserEvaluator evaluator{user_params.lambda, user_params.gamma,
                                     user_params.lambda_units,
                                     user_params.gamma_units};
  cooling::ApplyCoolingWithEvaluator(pmbp, w0, eos_data, runtime, evaluator, bdt,
                                     history_bdt, u0, gross_energy, net_energy);
}

void CoolingTestUserTimeStep(MeshBlockPack *pmbp, const DvceArray5D<Real> &w0,
                             const EOS_Data &eos_data,
                             const cooling::RuntimeData &runtime,
                             const cooling::TimestepData &timestep,
                             Real &dtnew) {
  CoolingTestUserEvaluator evaluator{user_params.lambda, user_params.gamma,
                                     user_params.lambda_units,
                                     user_params.gamma_units};
  cooling::CoolingNewDtWithEvaluator(pmbp, w0, eos_data, runtime, timestep,
                                     evaluator, dtnew);
}

} // namespace

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::CoolingTest()
//! \brief Initialize a uniform ideal-gas state for isolated cooling-source tests.

void ProblemGenerator::CoolingTest(ParameterInput *pin, const bool restart) {
  if (pin->GetOrAddBoolean("problem", "register_user_cooling", false)) {
    user_params.lambda = pin->GetOrAddReal("problem", "user_lambda", 0.0);
    user_params.gamma = pin->GetOrAddReal("problem", "user_gamma", 0.0);
    user_params.lambda_units = CoolingTestUnitSystem(
        pin->GetOrAddString("problem", "user_lambda_units", "code"),
        "problem/user_lambda_units");
    user_params.gamma_units = CoolingTestUnitSystem(
        pin->GetOrAddString("problem", "user_gamma_units", "code"),
        "problem/user_gamma_units");
    user_cooling_func = CoolingTestUserCooling;
    user_cooling_timestep_func = CoolingTestUserTimeStep;
  }

  // nothing needs to be done on restarts for this pgen
  if (restart) return;

  const Real density = pin->GetOrAddReal("problem", "density", 1.0);
  const Real pressure = pin->GetOrAddReal("problem", "pressure", 1.0);
  const Real vx = pin->GetOrAddReal("problem", "vx", 0.0);
  const Real vy = pin->GetOrAddReal("problem", "vy", 0.0);
  const Real vz = pin->GetOrAddReal("problem", "vz", 0.0);
  const Real scalar0 = pin->GetOrAddReal("problem", "scalar0", 0.0);
  const Real bx = pin->GetOrAddReal("problem", "bx", 0.0);
  const Real by = pin->GetOrAddReal("problem", "by", 0.0);
  const Real bz = pin->GetOrAddReal("problem", "bz", 0.0);

  if (density <= 0.0 || pressure <= 0.0) {
    CoolingTestFatal("problem/density and problem/pressure must be positive.");
  }

  auto &indcs = pmy_mesh_->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;

  if (pmbp->phydro == nullptr && pmbp->pmhd == nullptr) {
    CoolingTestFatal("cooling_test requires either a <hydro> or <mhd> block.");
  }

  if (pmbp->phydro != nullptr) {
    EOS_Data &eos = pmbp->phydro->peos->eos_data;
    if (!eos.is_ideal) CoolingTestFatal("cooling_test requires an ideal-gas EOS.");
    const Real gm1 = eos.gamma - 1.0;
    const Real eint = pressure/gm1;
    const Real ekin = 0.5*density*(vx*vx + vy*vy + vz*vz);
    const int nhydro = pmbp->phydro->nhydro;
    const int nscalars = pmbp->phydro->nscalars;
    auto &u0 = pmbp->phydro->u0;
    auto &w0 = pmbp->phydro->w0;

    par_for("cooling_test_hydro", DevExeSpace(), 0, (pmbp->nmb_thispack - 1),
    ks, ke, js, je, is, ie, KOKKOS_LAMBDA(int m, int k, int j, int i) {
      w0(m, IDN, k, j, i) = density;
      w0(m, IVX, k, j, i) = vx;
      w0(m, IVY, k, j, i) = vy;
      w0(m, IVZ, k, j, i) = vz;
      w0(m, IEN, k, j, i) = eint;
      u0(m, IDN, k, j, i) = density;
      u0(m, IM1, k, j, i) = density*vx;
      u0(m, IM2, k, j, i) = density*vy;
      u0(m, IM3, k, j, i) = density*vz;
      u0(m, IEN, k, j, i) = eint + ekin;
      for (int n = nhydro; n < nhydro + nscalars; ++n) {
        w0(m, n, k, j, i) = scalar0;
        u0(m, n, k, j, i) = density*scalar0;
      }
    });
  }

  if (pmbp->pmhd != nullptr) {
    EOS_Data &eos = pmbp->pmhd->peos->eos_data;
    if (!eos.is_ideal) CoolingTestFatal("cooling_test requires an ideal-gas EOS.");
    const Real gm1 = eos.gamma - 1.0;
    const Real eint = pressure/gm1;
    const Real ekin = 0.5*density*(vx*vx + vy*vy + vz*vz);
    const Real me = 0.5*(bx*bx + by*by + bz*bz);
    const int nmhd = pmbp->pmhd->nmhd;
    const int nscalars = pmbp->pmhd->nscalars;
    auto &u0 = pmbp->pmhd->u0;
    auto &w0 = pmbp->pmhd->w0;
    auto &bcc0 = pmbp->pmhd->bcc0;
    auto &b0 = pmbp->pmhd->b0;

    par_for("cooling_test_mhd", DevExeSpace(), 0, (pmbp->nmb_thispack - 1),
    ks, ke, js, je, is, ie, KOKKOS_LAMBDA(int m, int k, int j, int i) {
      w0(m, IDN, k, j, i) = density;
      w0(m, IVX, k, j, i) = vx;
      w0(m, IVY, k, j, i) = vy;
      w0(m, IVZ, k, j, i) = vz;
      w0(m, IEN, k, j, i) = eint;
      u0(m, IDN, k, j, i) = density;
      u0(m, IM1, k, j, i) = density*vx;
      u0(m, IM2, k, j, i) = density*vy;
      u0(m, IM3, k, j, i) = density*vz;
      u0(m, IEN, k, j, i) = eint + ekin + me;
      bcc0(m, IBX, k, j, i) = bx;
      bcc0(m, IBY, k, j, i) = by;
      bcc0(m, IBZ, k, j, i) = bz;
      b0.x1f(m, k, j, i) = bx;
      b0.x2f(m, k, j, i) = by;
      b0.x3f(m, k, j, i) = bz;
      if (i == ie) b0.x1f(m, k, j, i + 1) = bx;
      if (j == je) b0.x2f(m, k, j + 1, i) = by;
      if (k == ke) b0.x3f(m, k + 1, j, i) = bz;
      for (int n = nmhd; n < nmhd + nscalars; ++n) {
        w0(m, n, k, j, i) = scalar0;
        u0(m, n, k, j, i) = density*scalar0;
      }
    });
  }
}
