//========================================================================================
// AthenaK astrophysical fluid dynamics & numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cr_relativistic_pusher_runtime_test.cpp
//! \brief Serial production-path witness for the prescribed-field relativistic CR push.

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "mhd/mhd.hpp"
#include "parameter_input.hpp"
#include "particles/particles.hpp"
#include "particles/relativistic_state.hpp"
#include "pgen/pgen.hpp"

namespace {

[[noreturn]] void Fail(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl
            << message << std::endl;
  std::exit(EXIT_FAILURE);
}

void WriteRuntimeWitness(ParameterInput *pin, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (global_variable::nranks != 1) {
    Fail("CR relativistic prescribed runtime witness is intentionally serial-only");
  }
  if (pmbp == nullptr || pmbp->ppart == nullptr) {
    Fail("CR relativistic prescribed runtime witness requires particles");
  }
  if (pmbp->ppart->nprtcl_thispack != 1) {
    Fail("CR relativistic prescribed runtime witness requires exactly one particle");
  }

  auto pr_h = Kokkos::create_mirror_view_and_copy(
      HostMemSpace(), pmbp->ppart->prtcl_rdata);
  auto pi_h = Kokkos::create_mirror_view_and_copy(
      HostMemSpace(), pmbp->ppart->prtcl_idata);
  particles::relativistic::Vector3 w{
      pr_h(IPWX,0), pr_h(IPWY,0), pr_h(IPWZ,0)};
  Real gamma = 0.0;
  Real kinetic_energy = 0.0;
  if (particles::relativistic::GammaFromW(w, pmbp->ppart->c_model, gamma) !=
          particles::relativistic::StateStatus::success ||
      particles::relativistic::KineticEnergyFromW(
          w, pmbp->ppart->c_model, kinetic_energy) !=
          particles::relativistic::StateStatus::success) {
    Fail("CR relativistic prescribed runtime witness found invalid authoritative w");
  }

  const std::string filename = pin->GetOrAddString(
      "problem", "relativistic_runtime_witness", "cr_relativistic_pusher_runtime.tsv");
  std::ofstream stream(filename);
  if (!stream.is_open()) {
    Fail("could not open CR relativistic prescribed runtime witness output");
  }
  stream << "# cycle\ttime\ttag\tx\ty\tz\tvx\tvy\tvz\twx\twy\twz"
         << "\tcE1\tcE2\tcE3\tb1\tb2\tb3\twork\talpha\tgamma\tkinetic"
         << "\tdx\tdy\tdz\tdb\tlegacy_ipm\n";
  stream << std::setprecision(std::numeric_limits<Real>::max_digits10)
         << pm->ncycle << '\t' << pm->time << '\t' << pi_h(PTAG,0) << '\t'
         << pr_h(IPX,0) << '\t' << pr_h(IPY,0) << '\t' << pr_h(IPZ,0) << '\t'
         << pr_h(IPVX,0) << '\t' << pr_h(IPVY,0) << '\t' << pr_h(IPVZ,0) << '\t'
         << pr_h(IPWX,0) << '\t' << pr_h(IPWY,0) << '\t' << pr_h(IPWZ,0) << '\t'
         << pr_h(IPCEX,0) << '\t' << pr_h(IPCEY,0) << '\t' << pr_h(IPCEZ,0) << '\t'
         << pr_h(IPBX,0) << '\t' << pr_h(IPBY,0) << '\t' << pr_h(IPBZ,0) << '\t'
         << pr_h(IPWORK,0) << '\t' << pr_h(IPALPHA,0) << '\t' << gamma << '\t'
         << kinetic_energy << '\t' << pr_h(IPDX,0) << '\t' << pr_h(IPDY,0) << '\t'
         << pr_h(IPDZ,0) << '\t' << pr_h(IPDB,0) << '\t' << pr_h(IPM,0) << '\n';
}

} // namespace

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  PartRandom(pin, restart);
  pgen_final_func = &WriteRuntimeWitness;
  if (!pin->GetOrAddBoolean("problem", "runtime_poison_live_mhd", false)) {return;}
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &w0 = pmbp->pmhd->w0;
  auto &bcc0 = pmbp->pmhd->bcc0;
  auto &indcs = pmy_mesh_->mb_indcs;
  Kokkos::deep_copy(w0, 0.0);
  Kokkos::deep_copy(bcc0, 9.0);
  par_for("cr_rel_runtime_poison_live_mhd", DevExeSpace(),
          0, (pmbp->nmb_thispack - 1), indcs.ks, indcs.ke,
          indcs.js, indcs.je, indcs.is, indcs.ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    w0(m,IDN,k,j,i) = 1.0;
    w0(m,IVX,k,j,i) = 0.31;
    w0(m,IVY,k,j,i) = -0.27;
    w0(m,IVZ,k,j,i) = 0.19;
    w0(m,IEN,k,j,i) = 0.6;
  });
}
