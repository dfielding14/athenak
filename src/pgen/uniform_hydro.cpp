//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file uniform_hydro.cpp
//! \brief Uniform ideal-gas hydrodynamic state for performance measurements.

#include <cstdlib>
#include <iostream>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "pgen/pgen.hpp"

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::UserProblem()
//! \brief Initialize rho=P=1 and zero velocity in a hydro-only ideal-gas calculation.

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (restart) return;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->phydro == nullptr || pmbp->pmhd != nullptr) {
    std::cout << "### FATAL ERROR in uniform_hydro pgen" << std::endl
              << "uniform_hydro requires a <hydro> block and no <mhd> block."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  auto &eos = pmbp->phydro->peos->eos_data;
  if (!eos.is_ideal) {
    std::cout << "### FATAL ERROR in uniform_hydro pgen" << std::endl
              << "uniform_hydro requires <hydro>/eos=ideal." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  auto &indcs = pmy_mesh_->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nhydro = pmbp->phydro->nhydro;
  const int nscalars = pmbp->phydro->nscalars;
  const Real eint = 1.0/(eos.gamma - 1.0);
  auto &w0 = pmbp->phydro->w0;

  par_for("uniform_hydro", DevExeSpace(), 0, (pmbp->nmb_thispack - 1),
  ks, ke, js, je, is, ie, KOKKOS_LAMBDA(int m, int k, int j, int i) {
    w0(m, IDN, k, j, i) = 1.0;
    w0(m, IVX, k, j, i) = 0.0;
    w0(m, IVY, k, j, i) = 0.0;
    w0(m, IVZ, k, j, i) = 0.0;
    w0(m, IEN, k, j, i) = eint;
    for (int n = nhydro; n < nhydro + nscalars; ++n) {
      w0(m, n, k, j, i) = 0.0;
    }
  });

  auto &u0 = pmbp->phydro->u0;
  pmbp->phydro->peos->PrimToCons(w0, u0, is, ie, js, je, ks, ke);
}
