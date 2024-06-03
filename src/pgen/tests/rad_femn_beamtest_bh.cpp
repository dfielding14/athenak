//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file radiation_femn_beamtest_bh.cpp
//! \brief beam test on a schwarzschild background in isotropic coordinates

// C++ headers
#include <iostream>
#include <sstream>
#include <cmath>        // exp
#include <algorithm>    // max

// AthenaK headers
#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "pgen/pgen.hpp"
#include "radiation_femn/radiation_femn.hpp"
#include "adm/adm.hpp"

void ProblemGenerator::RadiationFEMNBeamtestBH(ParameterInput *pin, const bool restart) {
  if (restart) return;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;

  if (pmbp->pradfemn == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "The 2d beam test in Kerr-Schild problem generator can only be run with radiation-femn, but no "
              << "<radiation-femn> block in input file" << std::endl;
    exit(EXIT_FAILURE);
  }

  if (!pmbp->pmesh->two_d) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "The 2d beam test in Kerr-Schild problem generator can only be run with two dimensions, but parfile"
              << "grid setup is not in 2d" << std::endl;
    exit(EXIT_FAILURE);
  }

  if (pmbp->pradfemn->num_energy_bins != 1) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "The 2d beam test in Kerr-Schild problem generator can only be run with one energy bin!" << std::endl;
    exit(EXIT_FAILURE);
  }

  // capture var pmy_mesh_->mb_indcs;
  auto &indcs = pmy_mesh_->mb_indcs;
  auto &size = pmbp->pmb->mb_size;
  int &is = indcs.is;
  int &ie = indcs.ie;
  int &js = indcs.js;
  int &je = indcs.je;
  int &ks = indcs.ks;
  int &ke = indcs.ke;

  int isg = is - indcs.ng;
  int ieg = ie + indcs.ng;
  int jsg = (indcs.nx2 > 1) ? js - indcs.ng : js;
  int jeg = (indcs.nx2 > 1) ? je + indcs.ng : je;
  int ksg = (indcs.nx3 > 1) ? ks - indcs.ng : ks;
  int keg = (indcs.nx3 > 1) ? ke + indcs.ng : ke;
  int nmb = pmbp->nmb_thispack;
  auto &u_mu_ = pmbp->pradfemn->u_mu;
  adm::ADM::ADM_vars &adm = pmbp->padm->adm;
  Real ADM_mass = pin->GetOrAddReal("adm", "bh_mass", 1.);
  auto metric = pin->GetOrAddString("adm", "metric", "isotropic");

  std::cout << "Choice of metric: " << metric << std::endl;

  if (metric == "kerrschild") {
    par_for("pgen_linetest_metric_initialize_kerrschild", DevExeSpace(), 0, nmb - 1, ksg, keg, jsg, jeg, isg, ieg,
            KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
              Real &x1min = size.d_view(m).x1min;
              Real &x1max = size.d_view(m).x1max;
              int nx1 = indcs.nx1;
              Real x1v = CellCenterX(i - is, nx1, x1min, x1max);

              Real &x2min = size.d_view(m).x2min;
              Real &x2max = size.d_view(m).x2max;
              int nx2 = indcs.nx2;
              Real x2v = CellCenterX(j - js, nx2, x2min, x2max);

              Real r = std::sqrt(std::pow(x2v, 2) + std::pow(x1v, 2));

              for (int a = 0; a < 3; ++a)
                for (int b = a; b < 3; ++b) {
                  adm.g_dd(m, a, b, k, j, i) = (a == b ? 1. : 0.);
                }

              adm.psi4(m, k, j, i) = std::pow(1.0 + 0.5 * ADM_mass / r, 4);

              for (int a = 0; a < 3; ++a)
                for (int b = a; b < 3; ++b) {
                  adm.g_dd(m, a, b, k, j, i) *= adm.psi4(m, k, j, i);
                }

              adm.alpha(m, k, j, i) = (1.0 - 0.5 * ADM_mass / r) / (1.0 + 0.5 * ADM_mass / r);

              u_mu_(m, 0, k, j, i) = sqrt(1. / (adm.alpha(m, k, j, i) * adm.alpha(m, k, j, i)));
              u_mu_(m, 1, k, j, i) = 0.;
              u_mu_(m, 2, k, j, i) = 0.;
              u_mu_(m, 3, k, j, i) = 0.;
            });
  } else {
    par_for("pgen_linetest_metric_initialize_isotropic", DevExeSpace(), 0, nmb - 1, ksg, keg, jsg, jeg, isg, ieg,
            KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
              Real &x1min = size.d_view(m).x1min;
              Real &x1max = size.d_view(m).x1max;
              int nx1 = indcs.nx1;
              Real x1v = CellCenterX(i - is, nx1, x1min, x1max);

              Real &x2min = size.d_view(m).x2min;
              Real &x2max = size.d_view(m).x2max;
              int nx2 = indcs.nx2;
              Real x2v = CellCenterX(j - js, nx2, x2min, x2max);

              Real r = std::sqrt(std::pow(x2v, 2) + std::pow(x1v, 2));
              Real lvec[3] = {x1v / r, x2v / r, 0};

              for (int a = 0; a < 3; ++a)
                for (int b = a; b < 3; ++b) {
                  adm.g_dd(m, a, b, k, j, i) = (a == b ? 1. : 0.);
                  adm.g_dd(m, a, b, k, j, i) += 2. * ADM_mass * lvec[a] * lvec[b] / r;
                }

              adm.psi4(m, k, j, i) = 1.;
              adm.alpha(m, k, j, i) = 1. / sqrt(1.0 + 2. * ADM_mass / r);
              adm.beta_u(m, 0, k, j, i) = 2. * ADM_mass * adm.alpha(m, k, j, i) * adm.alpha(m, k, j, i) * lvec[0] / r;
              adm.beta_u(m, 1, k, j, i) = 2. * ADM_mass * adm.alpha(m, k, j, i) * adm.alpha(m, k, j, i) * lvec[1] / r;
              adm.beta_u(m, 2, k, j, i) = 2. * ADM_mass * adm.alpha(m, k, j, i) * adm.alpha(m, k, j, i) * lvec[2] / r;

              u_mu_(m, 0, k, j, i) = 1. / sqrt(1.0 - 2. * ADM_mass / r);
              u_mu_(m, 1, k, j, i) = 0.;
              u_mu_(m, 2, k, j, i) = 0.;
              u_mu_(m, 3, k, j, i) = 0.;
            });
  }

  return;
}