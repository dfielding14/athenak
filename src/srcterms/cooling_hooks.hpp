#ifndef SRCTERMS_COOLING_HOOKS_HPP_
#define SRCTERMS_COOLING_HOOKS_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cooling_hooks.hpp
//! \brief Shared templated launchers for pgen-defined cooling evaluators.

#include <cfloat>

#include "athena.hpp"
#include "cooling.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "units/units.hpp"

namespace cooling {

KOKKOS_INLINE_FUNCTION
Real DensityForCooling(DensityKind kind, UnitSystem units, Real rho_code,
                       const RuntimeData &runtime) {
  if (units == UnitSystem::cgs) {
    const Real rho_cgs = rho_code*runtime.density_cgs;
    if (kind == DensityKind::mass_density) return rho_cgs;
    if (kind == DensityKind::number_density) {
      return rho_cgs/(runtime.composition.mu*units::Units::atomic_mass_unit_cgs);
    }
    return rho_cgs*runtime.composition.hydrogen_mass_fraction/
           units::Units::atomic_mass_unit_cgs;
  }
  if (kind == DensityKind::mass_density) return rho_code;
  if (kind == DensityKind::number_density) return rho_code/runtime.composition.mu;
  return rho_code*runtime.composition.hydrogen_mass_fraction;
}

KOKKOS_INLINE_FUNCTION
Real CoolingEdotFromLambda(Real lambda, UnitSystem units, DensityKind density_kind,
                           Real rho_code, const RuntimeData &runtime) {
  const Real q = DensityForCooling(density_kind, units, rho_code, runtime);
  Real edot = q*q*lambda;
  if (units == UnitSystem::cgs) edot *= runtime.edot_cgs_to_code;
  return edot;
}

KOKKOS_INLINE_FUNCTION
Real HeatingEdotFromGamma(Real gamma, UnitSystem units, DensityKind density_kind,
                          Real rho_code, const RuntimeData &runtime) {
  const Real q = DensityForCooling(density_kind, units, rho_code, runtime);
  Real edot = q*gamma;
  if (units == UnitSystem::cgs) edot *= runtime.edot_cgs_to_code;
  return edot;
}

template <typename Evaluator>
void ApplyCoolingWithEvaluator(MeshBlockPack *pmbp, const DvceArray5D<Real> &w0,
                               const EOS_Data &eos_data,
                               const RuntimeData &runtime,
                               const Evaluator evaluator, const Real bdt,
                               DvceArray5D<Real> &u0,
                               Real &gross_energy, Real &net_energy) {
  auto &indcs = pmbp->pmesh->mb_indcs;
  const int is = indcs.is, nx1 = indcs.nx1;
  const int js = indcs.js, nx2 = indcs.nx2;
  const int ks = indcs.ks, nx3 = indcs.nx3;
  const int nmkji = (pmbp->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji = nx2*nx1;
  int nfluid = 0, nscalars = 0;
  if (pmbp->phydro != nullptr) {
    nfluid = pmbp->phydro->nhydro;
    nscalars = pmbp->phydro->nscalars;
  } else if (pmbp->pmhd != nullptr) {
    nfluid = pmbp->pmhd->nmhd;
    nscalars = pmbp->pmhd->nscalars;
  }

  const Real gm1 = eos_data.gamma - 1.0;
  const Real use_e = eos_data.use_e;
  auto &size = pmbp->pmb->mb_size;
  Kokkos::parallel_reduce("user_cooling_source",
                          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, Real &gross_sum, Real &net_sum) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;

    CoolingCellState state;
    state.rho_code = w0(m, IDN, k, j, i);
    state.rho_cgs = state.rho_code*runtime.density_cgs;
    state.nfluid = nfluid;
    state.nscalars = nscalars;
    if (use_e) {
      state.temp_code = w0(m, IEN, k, j, i)/state.rho_code*gm1;
      state.eint_code = w0(m, IEN, k, j, i);
    } else {
      state.temp_code = w0(m, ITM, k, j, i);
      state.eint_code = w0(m, ITM, k, j, i)*state.rho_code/gm1;
    }
    state.temp_cgs = state.temp_code*runtime.temp_cgs;
    if (nscalars > 0) state.scalar0 = w0(m, nfluid, k, j, i);

    const CoolingRates rates = evaluator(state);
    const Real edot_cool = CoolingEdotFromLambda(rates.lambda, rates.lambda_units,
        runtime.cooling_density, state.rho_code, runtime);
    const Real edot_heat = HeatingEdotFromGamma(rates.gamma, rates.gamma_units,
        runtime.heating_density, state.rho_code, runtime);
    const Real edot_net = edot_cool - edot_heat;
    const Real dvol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
    u0(m, IEN, k, j, i) -= bdt*edot_net;
    gross_sum += edot_cool*bdt*dvol;
    net_sum += edot_net*bdt*dvol;
  }, Kokkos::Sum<Real>(gross_energy), Kokkos::Sum<Real>(net_energy));
}

template <typename Evaluator>
void CoolingNewDtWithEvaluator(MeshBlockPack *pmbp, const DvceArray5D<Real> &w0,
                               const EOS_Data &eos_data,
                               const RuntimeData &runtime,
                               const TimestepData &timestep,
                               const Evaluator evaluator, Real &dtnew) {
  auto &indcs = pmbp->pmesh->mb_indcs;
  const int is = indcs.is, nx1 = indcs.nx1;
  const int js = indcs.js, nx2 = indcs.nx2;
  const int ks = indcs.ks, nx3 = indcs.nx3;
  const int nmkji = (pmbp->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji = nx2*nx1;
  int nfluid = 0, nscalars = 0;
  if (pmbp->phydro != nullptr) {
    nfluid = pmbp->phydro->nhydro;
    nscalars = pmbp->phydro->nscalars;
  } else if (pmbp->pmhd != nullptr) {
    nfluid = pmbp->pmhd->nmhd;
    nscalars = pmbp->pmhd->nscalars;
  }

  const Real gm1 = eos_data.gamma - 1.0;
  const Real use_e = eos_data.use_e;
  Kokkos::parallel_reduce("user_cooling_newdt",
                          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, Real &min_dt) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;

    CoolingCellState state;
    state.rho_code = w0(m, IDN, k, j, i);
    state.rho_cgs = state.rho_code*runtime.density_cgs;
    state.nfluid = nfluid;
    state.nscalars = nscalars;
    if (use_e) {
      state.temp_code = w0(m, IEN, k, j, i)/state.rho_code*gm1;
      state.eint_code = w0(m, IEN, k, j, i);
    } else {
      state.temp_code = w0(m, ITM, k, j, i);
      state.eint_code = w0(m, ITM, k, j, i)*state.rho_code/gm1;
    }
    state.temp_cgs = state.temp_code*runtime.temp_cgs;
    if (nscalars > 0) state.scalar0 = w0(m, nfluid, k, j, i);

    if (timestep.use_temperature_bounds) {
      const Real dt_temp = (timestep.temperature_units == UnitSystem::cgs) ?
                           state.temp_cgs : state.temp_code;
      if (dt_temp < timestep.temperature_min || dt_temp > timestep.temperature_max) {
        return;
      }
    }
    if (timestep.use_density_bounds) {
      const Real dt_density = DensityForCooling(timestep.density_kind,
          timestep.density_units, state.rho_code, runtime);
      if (dt_density < timestep.density_min || dt_density > timestep.density_max) {
        return;
      }
    }

    const CoolingRates rates = evaluator(state);
    const Real edot_cool = CoolingEdotFromLambda(rates.lambda, rates.lambda_units,
        runtime.cooling_density, state.rho_code, runtime);
    const Real edot_heat = HeatingEdotFromGamma(rates.gamma, rates.gamma_units,
        runtime.heating_density, state.rho_code, runtime);
    const Real rate = fabs(edot_cool - edot_heat) + FLT_MIN;
    min_dt = fmin(min_dt, timestep.factor*state.eint_code/rate);
  }, Kokkos::Min<Real>(dtnew));
}

} // namespace cooling

#endif // SRCTERMS_COOLING_HOOKS_HPP_
