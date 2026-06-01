//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file pic_paper_smooth_tsc_interface.cpp
//! \brief Deterministic static-SMR paper_smooth TSC interface regression carrier.

#include <iostream>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "outputs/restart_utils.hpp"
#include "particles/particles.hpp"
#include "pgen/pgen.hpp"

namespace {

[[noreturn]] void AbortInterfaceRegression(const char *message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl
            << "pic_paper_smooth_tsc_interface " << message << std::endl;
  restart_utils::AbortOnFatalError();
}

void RequireInterfaceMPISplit(const Mesh *pm) {
  int lower_fine_rank = -1;
  int coarse_rank = -1;
  for (int gid = 0; gid < pm->nmb_total; ++gid) {
    const auto &loc = pm->lloc_eachmb[gid];
    if (loc.level == pm->root_level + 1 &&
        loc.lx1 == 1 && loc.lx2 == 0 && loc.lx3 == 0) {
      lower_fine_rank = pm->rank_eachmb[gid];
    }
    if (loc.level == pm->root_level &&
        loc.lx1 == 1 && loc.lx2 == 0 && loc.lx3 == 0) {
      coarse_rank = pm->rank_eachmb[gid];
    }
  }
  if (lower_fine_rank < 0 || coarse_rank < 0) {
    AbortInterfaceRegression(
        "could not find the lower fine interface block and its coarse neighbor.");
  }
  if (lower_fine_rank == coarse_rank) {
    AbortInterfaceRegression(
        "requires the lower fine interface block and coarse neighbor on different ranks.");
  }
}

int FindLocalParticleMeshBlock(const MeshBlockPack *pmbp, const Real x1,
                               const Real x2, const Real x3) {
  const auto &size = pmbp->pmb->mb_size;
  for (int m = 0; m < pmbp->nmb_thispack; ++m) {
    const auto &block = size.h_view(m);
    if (x1 >= block.x1min && x1 < block.x1max &&
        x2 >= block.x2min && x2 < block.x2max &&
        x3 >= block.x3min && x3 < block.x3max) {
      return m;
    }
  }
  return -1;
}

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::PICPaperSmoothTSCInterface()
//! \brief Add one deterministic particle to a zero-amplitude LinearWave carrier.

void ProblemGenerator::PICPaperSmoothTSCInterface(ParameterInput *pin,
                                                  const bool restart) {
  LinearWave(pin, restart);

  Mesh *pm = pmy_mesh_;
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp == nullptr || pmbp->ppart == nullptr) {
    AbortInterfaceRegression("requires an active <particles> block.");
  }
  if (pmbp->ppart->particle_type != ParticleType::cosmic_ray) {
    AbortInterfaceRegression("requires <particles>/particle_type=cosmic_ray.");
  }

  if (pin->GetOrAddBoolean("problem", "require_interface_mpi_split", false)) {
    RequireInterfaceMPISplit(pm);
  }
  if (restart) return;

  auto *ppart = pmbp->ppart;
  pm->CountParticles();
  if (pm->nprtcl_total != 0) {
    AbortInterfaceRegression("requires <particles>/ppc=0.");
  }

  const Real particle_x = pin->GetReal("problem", "particle_x");
  const Real particle_y = pin->GetOrAddReal("problem", "particle_y", -0.75);
  const Real particle_z = pin->GetOrAddReal("problem", "particle_z", 0.0);
  const int m = FindLocalParticleMeshBlock(pmbp, particle_x, particle_y, particle_z);
  if (m >= 0) {
    HostArray2D<int> h_pi("paper_smooth_tsc_interface_pi", ppart->nidata, 1);
    HostArray2D<Real> h_pr("paper_smooth_tsc_interface_pr", ppart->nrdata, 1);
    for (int n = 0; n < ppart->nidata; ++n) h_pi(n, 0) = 0;
    for (int n = 0; n < ppart->nrdata; ++n) h_pr(n, 0) = 0.0;

    h_pi(PGID, 0) = pmbp->gids + m;
    h_pi(PTAG, 0) = 0;
    h_pi(PSP, 0) = 0;
    h_pi(PCRSOURCE, 0) = static_cast<int>(CRParticleSource::initial);

    h_pr(IPX, 0) = particle_x;
    h_pr(IPY, 0) = particle_y;
    h_pr(IPZ, 0) = particle_z;
    h_pr(IPVX, 0) = 0.0;
    h_pr(IPVY, 0) = 0.0;
    h_pr(IPVZ, 0) = 0.0;
    h_pr(IPM, 0) = 1.0;
    h_pr(IPWT, 0) = 1.0;
    h_pr(IPF0, 0) = 1.0;
    h_pr(IPT_BIRTH, 0) = pm->time;

    Kokkos::resize(ppart->prtcl_idata, ppart->nidata, 1);
    Kokkos::resize(ppart->prtcl_rdata, ppart->nrdata, 1);
    Kokkos::deep_copy(ppart->prtcl_idata, h_pi);
    Kokkos::deep_copy(ppart->prtcl_rdata, h_pr);
    ppart->nprtcl_thispack = 1;
  }

  pm->CountParticles();
  if (pm->nprtcl_total != 1) {
    AbortInterfaceRegression("requires exactly one global manually allocated particle.");
  }
}
