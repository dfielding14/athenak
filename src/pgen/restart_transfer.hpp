#ifndef PGEN_RESTART_TRANSFER_HPP_
#define PGEN_RESTART_TRANSFER_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================

#include "outputs/io_wrapper.hpp"

class Mesh;
class ParameterInput;

void ReadRestartMHD(Mesh *pm, ParameterInput *pin, IOWrapper &file,
                    IOWrapperSizeT offset, IOWrapperSizeT block_bytes);

#endif  // PGEN_RESTART_TRANSFER_HPP_
