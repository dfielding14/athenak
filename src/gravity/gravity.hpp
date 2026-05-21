#ifndef GRAVITY_GRAVITY_HPP_
#define GRAVITY_GRAVITY_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gravity.hpp
//! \brief defines Gravity module class

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../mesh/meshblock_pack.hpp"
#include "../parameter_input.hpp"
#include "../multigrid/multigrid.hpp"
#include "../coordinates/coordinates.hpp"
#include "mg_gravity.hpp"

class MeshBlockPack;
class ParameterInput;
class Coordinates;
class Multigrid;
class Driver;
namespace gravity {
class Gravity {
    public:
        Gravity(MeshBlockPack *pmbp, ParameterInput *pin);
        ~Gravity();

        MeshBlockPack* pmy_pack;  // ptr to MeshBlock containing this Field
        DvceArray5D<Real> phi;   // gravitational potential
        Real four_pi_G;
        MGGravityDriver *pmgd = nullptr;
        void SolveStage(Driver *pdriver, int stage, Real dt = 0.0);
};
} // namespace gravity
#endif // GRAVITY_GRAVITY_HPP_
