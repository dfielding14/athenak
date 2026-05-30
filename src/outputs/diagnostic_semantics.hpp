//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file diagnostic_semantics.hpp
//! \brief Pure geometry and Newtonian fluid diagnostic formulas for output kernels.

#ifndef OUTPUTS_DIAGNOSTIC_SEMANTICS_HPP_
#define OUTPUTS_DIAGNOSTIC_SEMANTICS_HPP_

#include <cmath>
#include <limits>

#ifdef KOKKOS_INLINE_FUNCTION
#define ATHENAK_DIAGNOSTIC_INLINE KOKKOS_INLINE_FUNCTION
#else
#define ATHENAK_DIAGNOSTIC_INLINE inline
#endif

namespace output_diagnostics {

template <typename T>
ATHENAK_DIAGNOSTIC_INLINE bool IsFinite(T value) {
  constexpr T maximum = std::numeric_limits<T>::max();
  return value == value && value <= maximum && value >= -maximum;
}

template <typename T>
ATHENAK_DIAGNOSTIC_INLINE T ClampUnit(T value) {
  if (value < static_cast<T>(-1.0)) return static_cast<T>(-1.0);
  if (value > static_cast<T>(1.0)) return static_cast<T>(1.0);
  return value;
}

template <typename T>
ATHENAK_DIAGNOSTIC_INLINE T VerticalSign(T z) {
  if (z > static_cast<T>(0.0)) return static_cast<T>(1.0);
  if (z < static_cast<T>(0.0)) return static_cast<T>(-1.0);
  return static_cast<T>(0.0);
}

template <typename T>
struct Geometry {
  T x;
  T y;
  T z;
  T cylindrical_radius;
  T radius;
  T theta;
  T phi;
  T costheta;
  T abscostheta;
  bool valid;
};

template <typename T>
ATHENAK_DIAGNOSTIC_INLINE Geometry<T> BuildGeometry(T x, T y, T z) {
  Geometry<T> geometry{};
  geometry.x = x;
  geometry.y = y;
  geometry.z = z;
  geometry.cylindrical_radius = sqrt(x*x + y*y);
  geometry.radius = sqrt(geometry.cylindrical_radius*geometry.cylindrical_radius + z*z);
  geometry.phi = atan2(y, x);
  if (geometry.phi < static_cast<T>(0.0)) {
    geometry.phi += static_cast<T>(6.283185307179586476925286766559);
  }
  geometry.costheta = geometry.radius > static_cast<T>(0.0)
      ? ClampUnit(z/geometry.radius) : static_cast<T>(1.0);
  geometry.abscostheta = fabs(geometry.costheta);
  geometry.theta = acos(geometry.costheta);
  geometry.valid = IsFinite(x) && IsFinite(y) && IsFinite(z) &&
      IsFinite(geometry.cylindrical_radius) && IsFinite(geometry.radius) &&
      IsFinite(geometry.theta) && IsFinite(geometry.phi) &&
      IsFinite(geometry.costheta) && IsFinite(geometry.abscostheta);
  return geometry;
}

template <typename T>
struct Flow {
  Geometry<T> geometry;
  T vx;
  T vy;
  T vz;
  T radial_velocity;
  T theta_velocity;
  T phi_velocity;
  T cylindrical_radial_velocity;
  T radial_mass_flux;
  T vertical_mass_flux;
  bool valid;
};

template <typename T>
ATHENAK_DIAGNOSTIC_INLINE Flow<T> BuildFlow(T x, T y, T z, T rho,
                                            T mx, T my, T mz) {
  Flow<T> flow{};
  flow.geometry = BuildGeometry(x, y, z);
  flow.valid = flow.geometry.valid && IsFinite(rho) &&
      rho > static_cast<T>(0.0) && IsFinite(mx) && IsFinite(my) && IsFinite(mz);
  if (!flow.valid) return flow;

  flow.vx = mx/rho;
  flow.vy = my/rho;
  flow.vz = mz/rho;
  const T radius = flow.geometry.radius;
  const T cylindrical_radius = flow.geometry.cylindrical_radius;
  flow.radial_velocity = radius > static_cast<T>(0.0)
      ? (flow.vx*x + flow.vy*y + flow.vz*z)/radius : static_cast<T>(0.0);
  flow.cylindrical_radial_velocity = cylindrical_radius > static_cast<T>(0.0)
      ? (flow.vx*x + flow.vy*y)/cylindrical_radius : static_cast<T>(0.0);
  flow.phi_velocity = cylindrical_radius > static_cast<T>(0.0)
      ? (-flow.vx*y + flow.vy*x)/cylindrical_radius : static_cast<T>(0.0);
  flow.theta_velocity = (radius > static_cast<T>(0.0) &&
                         cylindrical_radius > static_cast<T>(0.0))
      ? (z*(flow.vx*x + flow.vy*y)/(radius*cylindrical_radius)
         - flow.vz*cylindrical_radius/radius)
      : static_cast<T>(0.0);
  flow.radial_mass_flux = radius > static_cast<T>(0.0)
      ? (mx*x + my*y + mz*z)/radius : static_cast<T>(0.0);
  flow.vertical_mass_flux = VerticalSign(z)*mz;
  flow.valid = IsFinite(flow.vx) && IsFinite(flow.vy) && IsFinite(flow.vz) &&
      IsFinite(flow.radial_velocity) && IsFinite(flow.theta_velocity) &&
      IsFinite(flow.phi_velocity) && IsFinite(flow.cylindrical_radial_velocity) &&
      IsFinite(flow.radial_mass_flux) && IsFinite(flow.vertical_mass_flux);
  return flow;
}

template <typename T>
struct EnergyFlux {
  Flow<T> flow;
  T kinetic_radial;
  T magnetic_radial;
  T thermal_radial;
  T total_radial;
  T total_vertical;
  bool valid;
};

template <typename T>
ATHENAK_DIAGNOSTIC_INLINE EnergyFlux<T> BuildEnergyFlux(
    T x, T y, T z, T rho, T mx, T my, T mz, T total_energy, T gamma,
    bool is_mhd, T bx, T by, T bz, bool require_total_energy) {
  EnergyFlux<T> energy{};
  energy.flow = BuildFlow(x, y, z, rho, mx, my, mz);
  energy.valid = energy.flow.valid;
  if (!energy.valid) return energy;

  const T velocity_squared = energy.flow.vx*energy.flow.vx +
      energy.flow.vy*energy.flow.vy + energy.flow.vz*energy.flow.vz;
  T magnetic_squared = static_cast<T>(0.0);
  T velocity_dot_magnetic = static_cast<T>(0.0);
  T radial_magnetic = static_cast<T>(0.0);
  if (is_mhd) {
    energy.valid = IsFinite(bx) && IsFinite(by) && IsFinite(bz);
    if (!energy.valid) return energy;
    magnetic_squared = bx*bx + by*by + bz*bz;
    velocity_dot_magnetic = energy.flow.vx*bx + energy.flow.vy*by +
        energy.flow.vz*bz;
    radial_magnetic = energy.flow.geometry.radius > static_cast<T>(0.0)
        ? (bx*x + by*y + bz*z)/energy.flow.geometry.radius
        : static_cast<T>(0.0);
  }
  energy.kinetic_radial =
      static_cast<T>(0.5)*rho*velocity_squared*energy.flow.radial_velocity;
  energy.magnetic_radial =
      magnetic_squared*energy.flow.radial_velocity -
      velocity_dot_magnetic*radial_magnetic;
  if (require_total_energy) {
    energy.valid = IsFinite(total_energy) && IsFinite(gamma);
    if (!energy.valid) return energy;
    const T internal_energy = total_energy - static_cast<T>(0.5)*rho*velocity_squared
        - static_cast<T>(0.5)*magnetic_squared;
    const T enthalpy_plus_ke =
        static_cast<T>(0.5)*rho*velocity_squared + gamma*internal_energy;
    energy.thermal_radial = gamma*internal_energy*energy.flow.radial_velocity;
    energy.total_radial = (enthalpy_plus_ke + magnetic_squared)*
        energy.flow.radial_velocity - velocity_dot_magnetic*radial_magnetic;
    energy.total_vertical = (enthalpy_plus_ke + magnetic_squared)*energy.flow.vz -
        velocity_dot_magnetic*bz;
  }
  energy.valid = IsFinite(energy.kinetic_radial) &&
      IsFinite(energy.magnetic_radial) && IsFinite(energy.thermal_radial) &&
      IsFinite(energy.total_radial) && IsFinite(energy.total_vertical);
  return energy;
}

}  // namespace output_diagnostics

#undef ATHENAK_DIAGNOSTIC_INLINE

#endif  // OUTPUTS_DIAGNOSTIC_SEMANTICS_HPP_
