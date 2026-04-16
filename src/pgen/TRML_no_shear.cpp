//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file TRML_no_shear.cpp
//! \brief Zero-shear-compatible wrapper around the TRML problem generator.
//!
//! This variant enables <problem>/velocity = 0 while keeping the main TRML
//! implementation shared with src/pgen/TRML.cpp. When run without shear, a positive
//! explicit <problem>/t_cool_0 is required to normalize the cooling/heating source term.

#define TRML_ALLOW_ZERO_SHEAR 1
#include "TRML.cpp"
