#ifndef PROFILE_READER_HPP_
#define PROFILE_READER_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file  profile_reader.hpp
//  \brief Allows for reading in profiles for initial conditions and boundary conditions

#include "athena.hpp"
#include <fstream>  // For std::ifstream
#include <string>   // For std::string
#include <sstream>  // For std::istringstream
#include <vector>   // For std::vector

// Forward declaration
class ProfileReaderHost;

// Device-accessible class for profile interpolation
class ProfileReader {
  public:
  // Default constructor
  KOKKOS_INLINE_FUNCTION ProfileReader() {}
  
  // Copy constructor from host reader
  ProfileReader(const ProfileReaderHost& host_reader);

  // Destructor
  KOKKOS_INLINE_FUNCTION
  ~ProfileReader() {
    // Explicitly reset all Views to empty Views
    d_r_vals    = DvceArray1D<Real>();
    d_rho_vals  = DvceArray1D<Real>();
    d_temp_vals = DvceArray1D<Real>();
    d_v_vals    = DvceArray1D<Real>();
  }
  
  // Interpolation functions
  KOKKOS_INLINE_FUNCTION Real GetDensity(Real r) const;
  KOKKOS_INLINE_FUNCTION Real GetTemperature(Real r) const;
  KOKKOS_INLINE_FUNCTION Real GetVelocity(Real r) const;
  KOKKOS_INLINE_FUNCTION Real GetRmin() const;
  
  private:
  // Device views for profile data
  DvceArray1D<Real> d_r_vals;
  DvceArray1D<Real> d_rho_vals;
  DvceArray1D<Real> d_temp_vals;
  DvceArray1D<Real> d_v_vals;
  int num_points;
  
  // Helper for interpolation
  KOKKOS_INLINE_FUNCTION 
  Real Interpolate(const DvceArray1D<Real> &x,
                   const DvceArray1D<Real> &y,
                   Real x_val) const;
};

// Host-only reader that loads the data
class ProfileReaderHost {
public:
  ProfileReaderHost() = default;
  
  // Read profiles from file on host
  void ReadProfiles(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
      throw std::runtime_error("Failed to open profile file: " + filename);
    }
    
    // Read header line and skip it
    std::string header;
    std::getline(file, header);
    
    // Read data
    std::string line;
    while (std::getline(file, line)) {
      std::istringstream iss(line);
      Real r_val, T_val, rho_val, v_val;
      
      if (!(iss >> r_val >> T_val >> rho_val >> v_val)) {
        continue; // Skip invalid lines
      }
      
      h_r_vals.push_back(r_val);
      h_temp_vals.push_back(T_val);
      h_rho_vals.push_back(rho_val);
      h_v_vals.push_back(v_val);
    }
  
    // Check if we have data
    if (h_r_vals.empty()) {
      throw std::runtime_error("No valid data found in profile file");
    }
  }
  
  // Create device-accessible reader
  ProfileReader CreateDeviceReader() const{
    return ProfileReader(*this);
  }
  
  friend class ProfileReader;
  
private:
  // Host vectors for profile data
  std::vector<Real> h_r_vals;
  std::vector<Real> h_rho_vals;
  std::vector<Real> h_temp_vals;
  std::vector<Real> h_v_vals;
};

// Implementation of device reader constructor from host data
ProfileReader::ProfileReader(const ProfileReaderHost& host_reader) {
  // Create device views with appropriate size
  num_points = host_reader.h_r_vals.size();
  
  d_r_vals    = DvceArray1D<Real>("r_vals", num_points);
  d_rho_vals  = DvceArray1D<Real>("rho_vals", num_points);
  d_temp_vals = DvceArray1D<Real>("temp_vals", num_points);
  d_v_vals    = DvceArray1D<Real>("v_vals", num_points);
  
  // Create host mirrors
  auto h_r    = Kokkos::create_mirror_view(d_r_vals);
  auto h_rho  = Kokkos::create_mirror_view(d_rho_vals);
  auto h_temp = Kokkos::create_mirror_view(d_temp_vals);
  auto h_v    = Kokkos::create_mirror_view(d_v_vals);
  
  // Copy data from host vectors to mirrors
  for (int i = 0; i < num_points; i++) {
    h_r(i) = host_reader.h_r_vals[i];
    h_rho(i) = host_reader.h_rho_vals[i];
    h_temp(i) = host_reader.h_temp_vals[i];
    h_v(i) = host_reader.h_v_vals[i];
  }
  
  // Copy from mirrors to device
  Kokkos::deep_copy(d_r_vals, h_r);
  Kokkos::deep_copy(d_rho_vals, h_rho);
  Kokkos::deep_copy(d_temp_vals, h_temp);
  Kokkos::deep_copy(d_v_vals, h_v);
}

// Device-compatible interpolation implementation
KOKKOS_INLINE_FUNCTION
Real ProfileReader::Interpolate(
  const DvceArray1D<Real> &x,
  const DvceArray1D<Real> &y, 
  Real x_val) const 
{
  // Handle out of bounds
  if (x_val <= x(0)) return y(0);
  if (x_val >= x(num_points-1)) return y(num_points-1);
  
  // Find position using binary search
  int left = 0;
  int right = num_points - 1;
  
  while (right - left > 1) {
    int mid = left + (right - left) / 2;
    if (x(mid) > x_val) right = mid;
    else left = mid;
  }
  
  // Linear interpolation
  Real x0 = x(left);
  Real x1 = x(right);
  Real y0 = y(left);
  Real y1 = y(right);
  
  return y0 + (x_val - x0) * (y1 - y0) / (x1 - x0);
}

// Implemented getter functions
KOKKOS_INLINE_FUNCTION
Real ProfileReader::GetDensity(Real r) const {
  return Interpolate(d_r_vals, d_rho_vals, r);
}

KOKKOS_INLINE_FUNCTION
Real ProfileReader::GetTemperature(Real r) const {
  return Interpolate(d_r_vals, d_temp_vals, r);
}

KOKKOS_INLINE_FUNCTION
Real ProfileReader::GetVelocity(Real r) const {
  return Interpolate(d_r_vals, d_v_vals, r);
}

KOKKOS_INLINE_FUNCTION
Real ProfileReader::GetRmin() const {
  // Return the minimum radius from the first value in the r_vals array
  return d_r_vals(0);
}

#endif // PROFILE_READER_HPP_
