//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file star_particle_test.cpp
//! \brief Small hydro setup for star-particle formation and accretion tests.

#include <algorithm>
#include <cmath>
#include <iostream>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "parameter_input.hpp"
#include "pgen/pgen.hpp"

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::StarParticleTest()
//! \brief Initializes a uniform gas box with one denser central cell.

void ProblemGenerator::StarParticleTest(ParameterInput *pin, const bool restart) {
  if (restart) return;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->phydro == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "star_particle_test requires a <hydro> block"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  Real base_density = pin->GetOrAddReal("problem", "base_density", 1.0);
  Real peak_density = pin->GetOrAddReal("problem", "peak_density", 10.0);
  Real pressure = pin->GetOrAddReal("problem", "pressure", 1.0);
  Real vx = pin->GetOrAddReal("problem", "vx", 0.0);
  Real vy = pin->GetOrAddReal("problem", "vy", 0.0);
  Real vz = pin->GetOrAddReal("problem", "vz", 0.0);
  Real x1c = pin->GetOrAddReal("problem", "x1c", 0.5);
  Real x2c = pin->GetOrAddReal("problem", "x2c", 0.5);
  Real x3c = pin->GetOrAddReal("problem", "x3c", 0.5);

  auto &indcs = pmy_mesh_->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  auto &size = pmbp->pmb->mb_size;
  auto &w0 = pmbp->phydro->w0;
  int nhydro = pmbp->phydro->nhydro;
  int nscalars = pmbp->phydro->nscalars;
  auto &eos = pmbp->phydro->peos->eos_data;

  par_for("star_particle_test", DevExeSpace(), 0, (pmbp->nmb_thispack-1),
  ks, ke, js, je, is, ie, KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real x1 = CellCenterX(i - is, indcs.nx1, size.d_view(m).x1min,
                          size.d_view(m).x1max);
    Real x2 = CellCenterX(j - js, indcs.nx2, size.d_view(m).x2min,
                          size.d_view(m).x2max);
    Real x3 = CellCenterX(k - ks, indcs.nx3, size.d_view(m).x3min,
                          size.d_view(m).x3max);
    bool in_peak = (fabs(x1 - x1c) <= 0.5*size.d_view(m).dx1 &&
                    fabs(x2 - x2c) <= 0.5*size.d_view(m).dx2 &&
                    fabs(x3 - x3c) <= 0.5*size.d_view(m).dx3);

    w0(m,IDN,k,j,i) = in_peak ? peak_density : base_density;
    w0(m,IVX,k,j,i) = vx;
    w0(m,IVY,k,j,i) = vy;
    w0(m,IVZ,k,j,i) = vz;
    if (nhydro > IEN) {
      w0(m,IEN,k,j,i) = pressure/(eos.gamma - 1.0);
    }
    for (int n=0; n<nscalars; ++n) {
      w0(m,IYF+n,k,j,i) = 0.0;
    }
  });

  pmbp->phydro->peos->PrimToCons(w0, pmbp->phydro->u0, is, ie, js, je, ks, ke);
}
