#ifndef PGEN_RESTART_REFINEMENT_HPP_
#define PGEN_RESTART_REFINEMENT_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================

#include "outputs/io_wrapper.hpp"

class Mesh;
class ParameterInput;

void ReadRefinedRestartMHD(Mesh *pm, ParameterInput *pin, IOWrapper &file,
                          IOWrapperSizeT offset, IOWrapperSizeT block_bytes);

#endif  // PGEN_RESTART_REFINEMENT_HPP_
