#ifndef UTILS_TRML_COOLING_HPP_
#define UTILS_TRML_COOLING_HPP_
//========================================================================================
//! \file trml_cooling.hpp
//! \brief Shared cooling/heating helpers for TRML problem generators and diagnostics.
//========================================================================================

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "units/units.hpp"

namespace trml_cooling {

enum HeatingMode {
  kHeatingOff = 0,
  kHeatingDirect = 1,
  kHeatingBalanceCold = 2
};

struct Params {
  Real T_peak = 1.0;
  Real cool_T_min = 0.0;
  Real cool_T_max = std::numeric_limits<Real>::infinity();
  Real lambda_peak = 0.0;
  Real lambda_slope_lo = 0.0;
  Real lambda_slope_hi = 0.0;
  Real lambda_smooth_width = 0.05;
  Real heat_gamma = 0.0;
  Real heat_gamma_T_ref = 1.0;
  Real heat_gamma_T_slope = 0.0;
  Real mass_per_particle = 1.0;
  bool use_dens_ceiling = false;
  Real dens_ceiling = std::numeric_limits<Real>::max();
};

struct Rates {
  Real edot_cool = 0.0;
  Real edot_heat = 0.0;
  Real edot_net = 0.0;
};

struct ReferenceState {
  Real rho_0 = 1.0;
  Real pgas_0 = 1.0;
  Real contrast = 1.0;
  Real gamma = 5.0/3.0;
  Real t_shear = std::numeric_limits<Real>::infinity();
  bool allow_zero_shear = false;
};

struct Setup {
  Params params;
  Real t_cool_0 = std::numeric_limits<Real>::infinity();
  Real xi = std::numeric_limits<Real>::infinity();
  int heating_mode = kHeatingOff;
  bool valid = false;
  bool direct_lambda = false;
};

KOKKOS_INLINE_FUNCTION
Real NumberDensity(const Params &p, const Real dens) {
  return dens/std::max(p.mass_per_particle, static_cast<Real>(1.0e-300));
}

KOKKOS_INLINE_FUNCTION
Real Lambda(const Params &p, const Real temp) {
  if (temp <= 0.0 || p.T_peak <= 0.0 || p.lambda_peak <= 0.0) return 0.0;

  const Real x = std::log(temp/p.T_peak);
  const Real width = std::max(p.lambda_smooth_width, static_cast<Real>(1.0e-12));
  Real ln_shape = 0.0;

  if (x < 0.0) {
    ln_shape = p.lambda_slope_lo*(x + width*(1.0 - std::exp(x/width)));
  } else {
    ln_shape = p.lambda_slope_hi*(x + width*(std::exp(-x/width) - 1.0));
  }
  return p.lambda_peak*std::exp(ln_shape);
}

KOKKOS_INLINE_FUNCTION
Real HeatingShape(const Params &p, const Real temp) {
  if (temp <= 0.0 || p.heat_gamma_T_ref <= 0.0) return 0.0;
  if (p.heat_gamma_T_slope == 0.0) return 1.0;
  return std::pow(temp/p.heat_gamma_T_ref, p.heat_gamma_T_slope);
}

KOKKOS_INLINE_FUNCTION
Real Gamma(const Params &p, const Real temp) {
  if (p.heat_gamma <= 0.0) return 0.0;
  return p.heat_gamma*HeatingShape(p, temp);
}

KOKKOS_INLINE_FUNCTION
Rates Evaluate(const Params &p, const Real dens, const Real temp, const Real /*eint*/) {
  Rates r;
  if (dens <= 0.0 || temp <= 0.0) return r;
  if (temp < p.cool_T_min || temp > p.cool_T_max) return r;
  if (p.use_dens_ceiling && dens > p.dens_ceiling) return r;

  const Real n = NumberDensity(p, dens);
  r.edot_cool = n*n*Lambda(p, temp);
  r.edot_heat = n*Gamma(p, temp);
  r.edot_net = r.edot_cool - r.edot_heat;
  return r;
}

KOKKOS_INLINE_FUNCTION
Real CoolingTimeAtPeak(const Params &p, const ReferenceState &ref) {
  if (p.lambda_peak <= 0.0) return std::numeric_limits<Real>::infinity();
  const Real gm1 = ref.gamma - 1.0;
  const Real dens_peak = ref.pgas_0/p.T_peak;
  const Real eint_peak = ref.pgas_0/gm1;
  const Real n_peak = NumberDensity(p, dens_peak);
  const Real edot_peak = n_peak*n_peak*p.lambda_peak;
  return eint_peak/std::max(edot_peak, static_cast<Real>(1.0e-300));
}

inline bool HasProblem(ParameterInput *pin, const char *name) {
  return pin->DoesParameterExist("problem", name);
}

inline void Fatal(const std::string &message) {
  std::cout << "### FATAL ERROR in TRML cooling setup" << std::endl
            << message << std::endl;
  std::exit(EXIT_FAILURE);
}

inline void Warn(const std::string &message, const bool verbose) {
  if (!verbose) return;
  std::cout << "### WARNING in TRML cooling setup" << std::endl
            << message << std::endl;
}

inline bool GetRealIfPresent(ParameterInput *pin, const char *name, Real &value) {
  if (!HasProblem(pin, name)) return false;
  value = pin->GetReal("problem", name);
  return true;
}

inline std::string HeatingModeName(const int mode) {
  if (mode == kHeatingOff) return "off";
  if (mode == kHeatingDirect) return "direct";
  if (mode == kHeatingBalanceCold) return "balance_cold";
  return "unknown";
}

inline int ParseHeatingMode(const std::string &mode) {
  if (mode == "off" || mode == "none" || mode == "false") return kHeatingOff;
  if (mode == "direct" || mode == "constant") return kHeatingDirect;
  if (mode == "balance_cold" || mode == "balanced_cold") return kHeatingBalanceCold;
  Fatal("Invalid <problem>/heating_mode='" + mode +
        "'. Expected off, direct, or balance_cold.");
  return kHeatingOff;
}

inline Setup ReadInputs(ParameterInput *pin, const units::Units *punit,
                        const ReferenceState &ref, const bool fatal_on_missing,
                        const bool verbose) {
  Setup setup;
  Params &p = setup.params;
  p.mass_per_particle = punit->mu()*punit->atomic_mass_unit();
  const Real T_cold = ref.pgas_0/ref.rho_0/ref.contrast;

  auto missing = [&](const std::string &message) {
    if (fatal_on_missing) Fatal(message);
    setup.valid = false;
    return setup;
  };

  if (HasProblem(pin, "T_peak")) {
    p.T_peak = pin->GetReal("problem", "T_peak");
  } else if (HasProblem(pin, "T_peak_over_T_cold")) {
    p.T_peak = T_cold*pin->GetReal("problem", "T_peak_over_T_cold");
  } else {
    return missing("Missing T_peak or T_peak_over_T_cold.");
  }

  if (HasProblem(pin, "lambda_slope_lo") || HasProblem(pin, "lambda_slope_hi")) {
    if (!HasProblem(pin, "lambda_slope_lo") || !HasProblem(pin, "lambda_slope_hi")) {
      return missing("Specify both lambda_slope_lo and lambda_slope_hi.");
    }
    p.lambda_slope_lo = pin->GetReal("problem", "lambda_slope_lo");
    p.lambda_slope_hi = pin->GetReal("problem", "lambda_slope_hi");
  } else if (HasProblem(pin, "beta_lo") && HasProblem(pin, "beta_hi")) {
    p.lambda_slope_lo = -pin->GetReal("problem", "beta_lo");
    p.lambda_slope_hi = -pin->GetReal("problem", "beta_hi");
    Warn("beta_lo/beta_hi are deprecated; use lambda_slope_lo/lambda_slope_hi. "
         "The legacy beta signs were mapped to actual Lambda(T) slopes.", verbose);
  } else {
    return missing("Missing lambda_slope_lo/lambda_slope_hi or legacy beta_lo/beta_hi.");
  }

  p.lambda_smooth_width = HasProblem(pin, "lambda_smooth_width") ?
      pin->GetReal("problem", "lambda_smooth_width") : 0.05;
  p.heat_gamma_T_ref = HasProblem(pin, "heat_gamma_T_ref") ?
      pin->GetReal("problem", "heat_gamma_T_ref") : T_cold;
  p.heat_gamma_T_slope = HasProblem(pin, "heat_gamma_T_slope") ?
      pin->GetReal("problem", "heat_gamma_T_slope") : 0.0;

  if (HasProblem(pin, "cool_T_min")) {
    p.cool_T_min = pin->GetReal("problem", "cool_T_min");
  } else if (HasProblem(pin, "cooling_below_T_cold") &&
             !pin->GetBoolean("problem", "cooling_below_T_cold")) {
    p.cool_T_min = ref.pgas_0/ref.rho_0/ref.contrast;
    Warn("cooling_below_T_cold=false is deprecated; using cool_T_min=T_cold.", verbose);
  } else {
    p.cool_T_min = 0.0;
  }

  if (HasProblem(pin, "cool_T_max")) {
    p.cool_T_max = pin->GetReal("problem", "cool_T_max");
  } else if (HasProblem(pin, "T_cutoff_over_T_hot")) {
    const Real T_hot = ref.pgas_0/ref.rho_0;
    p.cool_T_max = T_hot*pin->GetReal("problem", "T_cutoff_over_T_hot");
    Warn("T_cutoff_over_T_hot is deprecated; use cool_T_max.", verbose);
  } else {
    p.cool_T_max = std::numeric_limits<Real>::infinity();
  }

  p.use_dens_ceiling = HasProblem(pin, "use_dens_ceiling_cool") ?
      pin->GetBoolean("problem", "use_dens_ceiling_cool") : false;
  if (p.use_dens_ceiling) {
    const Real n0_ceiling = HasProblem(pin, "n0_ceiling_cool") ?
        pin->GetReal("problem", "n0_ceiling_cool") : 20.0;
    p.dens_ceiling = n0_ceiling*p.mass_per_particle;
  }

  if (p.T_peak <= 0.0) Fatal("T_peak must be positive.");
  if (p.cool_T_min < 0.0 || !(p.cool_T_min < p.cool_T_max)) {
    Fatal("Require 0 <= cool_T_min < cool_T_max.");
  }
  const Real window_tol = 64.0*std::numeric_limits<Real>::epsilon()*
      std::max(static_cast<Real>(1.0),
               std::max(std::abs(p.T_peak),
                        std::max(std::abs(p.cool_T_min), std::abs(p.cool_T_max))));
  if (p.T_peak < p.cool_T_min - window_tol || p.T_peak > p.cool_T_max + window_tol) {
    Fatal("T_peak must lie inside the active cooling/heating temperature window.");
  }
  if (p.lambda_smooth_width <= 0.0) Fatal("lambda_smooth_width must be positive.");
  if (p.heat_gamma_T_ref <= 0.0) Fatal("heat_gamma_T_ref must be positive.");
  if (p.lambda_slope_lo < 0.0 || p.lambda_slope_hi > 0.0) {
    Fatal("Require lambda_slope_lo >= 0 and lambda_slope_hi <= 0.");
  }

  setup.direct_lambda = GetRealIfPresent(pin, "lambda_peak", p.lambda_peak);
  if (setup.direct_lambda) {
    if (p.lambda_peak < 0.0) Fatal("lambda_peak must be non-negative.");
    setup.t_cool_0 = CoolingTimeAtPeak(p, ref);
    setup.xi = std::isfinite(ref.t_shear) ? ref.t_shear/setup.t_cool_0 :
               std::numeric_limits<Real>::infinity();
    if (HasProblem(pin, "xi") || HasProblem(pin, "t_cool_0")) {
      Warn("lambda_peak overrides xi/t_cool_0. Effective values are printed for "
           "diagnostics.", verbose);
    }
  } else {
    const bool have_xi = HasProblem(pin, "xi");
    const Real xi_input = have_xi ? pin->GetReal("problem", "xi") : -1.0;
    const Real t_cool_0_input = HasProblem(pin, "t_cool_0") ?
        pin->GetReal("problem", "t_cool_0") : -1.0;

    if (t_cool_0_input > 0.0) {
      setup.t_cool_0 = t_cool_0_input;
    } else if (std::isfinite(ref.t_shear)) {
      if (!have_xi || xi_input <= 0.0) {
        return missing("Require lambda_peak, positive t_cool_0, or positive xi.");
      }
      setup.t_cool_0 = ref.t_shear/xi_input;
    } else {
      return missing("Zero-shear TRML requires lambda_peak or positive t_cool_0.");
    }

    const Real gm1 = ref.gamma - 1.0;
    const Real dens_peak = ref.pgas_0/p.T_peak;
    const Real eint_peak = ref.pgas_0/gm1;
    const Real n_peak = NumberDensity(p, dens_peak);
    p.lambda_peak = eint_peak/(setup.t_cool_0*n_peak*n_peak);
    setup.xi = std::isfinite(ref.t_shear) ? ref.t_shear/setup.t_cool_0 :
               std::numeric_limits<Real>::infinity();
  }

  if (HasProblem(pin, "heating_mode")) {
    setup.heating_mode = ParseHeatingMode(pin->GetString("problem", "heating_mode"));
  } else if (HasProblem(pin, "heat_gamma")) {
    setup.heating_mode = kHeatingDirect;
  } else if (HasProblem(pin, "heating_on")) {
    const bool heating_on = pin->GetBoolean("problem", "heating_on");
    if (!heating_on) {
      setup.heating_mode = kHeatingOff;
      Warn("heating_on=false is deprecated; use heating_mode=off.", verbose);
    } else if (HasProblem(pin, "balanced_heating") &&
               pin->GetBoolean("problem", "balanced_heating")) {
      setup.heating_mode = kHeatingBalanceCold;
      Warn("balanced_heating is deprecated; using heating_mode=balance_cold.", verbose);
    } else {
      return missing("heating_on=true is deprecated. Specify heating_mode=direct "
                     "with heat_gamma, or heating_mode=balance_cold.");
    }
  } else {
    setup.heating_mode = kHeatingOff;
  }

  if (HasProblem(pin, "balanced_heating")) {
    Warn("balanced_heating is deprecated and no longer selects the old pressure-squared "
         "heating law.", verbose);
  }

  if (setup.heating_mode == kHeatingOff) {
    p.heat_gamma = 0.0;
  } else if (setup.heating_mode == kHeatingDirect) {
    p.heat_gamma = HasProblem(pin, "heat_gamma") ? pin->GetReal("problem", "heat_gamma") : 0.0;
    if (p.heat_gamma < 0.0) Fatal("heat_gamma must be non-negative.");
  } else if (setup.heating_mode == kHeatingBalanceCold) {
    if (T_cold < p.cool_T_min || T_cold > p.cool_T_max) {
      Fatal("heating_mode=balance_cold requires T_cold inside [cool_T_min, cool_T_max].");
    }
    const Real rho_cold = ref.rho_0*ref.contrast;
    const Real n_cold = NumberDensity(p, rho_cold);
    const Real heat_shape_cold = HeatingShape(p, T_cold);
    p.heat_gamma = n_cold*Lambda(p, T_cold)/
                   std::max(heat_shape_cold, static_cast<Real>(1.0e-300));
  }

  setup.valid = true;
  return setup;
}

} // namespace trml_cooling

#endif // UTILS_TRML_COOLING_HPP_
