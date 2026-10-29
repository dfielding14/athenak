#ifndef MESH_RESTART_PROLONGATION_HPP_
#define MESH_RESTART_PROLONGATION_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================

#include "athena.hpp"

class Mesh;
class ParameterInput;

void ProlongRestartMHD(Mesh *pm, ParameterInput *pin,
                       DvceArray5D<Real> &coarse_cons, DvceFaceFld4D<Real> &coarse_b);

#endif  // MESH_RESTART_PROLONGATION_HPP_
