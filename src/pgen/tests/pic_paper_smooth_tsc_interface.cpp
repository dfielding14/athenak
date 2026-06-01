//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file pic_paper_smooth_tsc_interface.cpp
//! \brief Deterministic static-SMR paper_smooth TSC interface regression carrier.

#include <iostream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "mesh/nghbr_index.hpp"
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
        "requires the lower fine interface block and coarse neighbor on "
        "different ranks.");
  }
}

int FindLogicalMeshBlockGID(const Mesh *pm, const int level,
                            const int lx1, const int lx2, const int lx3) {
  for (int gid = 0; gid < pm->nmb_total; ++gid) {
    const auto &loc = pm->lloc_eachmb[gid];
    if (loc.level == level && loc.lx1 == lx1 &&
        loc.lx2 == lx2 && loc.lx3 == lx3) {
      return gid;
    }
  }
  return -1;
}

void RequireRemoteReceiverMPISplit(const Mesh *pm, const MeshBlockPack *pmbp,
                                   const std::string &requirement) {
  if (requirement == "none") return;

  int owner_gid = -1;
  int receiver_gid = -1;
  int slot = -1;
  if (requirement == "g_periodic_x1") {
    owner_gid = FindLogicalMeshBlockGID(pm, pm->root_level + 1, 0, 0, 0);
    receiver_gid = FindLogicalMeshBlockGID(pm, pm->root_level, 1, 0, 0);
    slot = NeighborIndex(-1, 0, 0, 0, 0);
  } else if (requirement == "j_periodic_x1_x3_edge") {
    owner_gid = FindLogicalMeshBlockGID(pm, pm->root_level, 1, 0, 0);
    receiver_gid = FindLogicalMeshBlockGID(pm, pm->root_level + 1, 0, 0, 1);
    slot = NeighborIndex(1, 0, -1, 0, 0);
  } else {
    AbortInterfaceRegression("has an unsupported required remote receiver split.");
  }

  if (owner_gid < 0 || receiver_gid < 0 || slot < 0) {
    AbortInterfaceRegression("could not resolve the required remote receiver split.");
  }
  const int owner_rank = pm->rank_eachmb[owner_gid];
  const int receiver_rank = pm->rank_eachmb[receiver_gid];
  if (owner_rank == receiver_rank) {
    AbortInterfaceRegression("requires the selected receiver on a remote MPI rank.");
  }
  if (owner_rank == global_variable::my_rank) {
    const int owner_m = owner_gid - pmbp->gids;
    if (owner_m < 0 || owner_m >= pmbp->nmb_thispack) {
      AbortInterfaceRegression("could not resolve the local owner MeshBlock.");
    }
    const auto &neighbor = pmbp->pmb->nghbr.h_view(owner_m, slot);
    if (neighbor.gid != receiver_gid || neighbor.rank != receiver_rank) {
      AbortInterfaceRegression("required neighbor slot does not name the receiver.");
    }
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
  if (!pmbp->ppart->UsesPaperVL2Coupling()) {
    AbortInterfaceRegression(
        "requires <particles>/pic_physical_mode=paper_mhd_pic_vl2_tsc.");
  }

  if (pin->GetOrAddBoolean("problem", "require_interface_mpi_split", false)) {
    RequireInterfaceMPISplit(pm);
  }
  RequireRemoteReceiverMPISplit(
      pm, pmbp, pin->GetOrAddString(
          "problem", "required_remote_receiver_mpi_split", "none"));
  if (restart) return;

  auto *ppart = pmbp->ppart;
  pm->CountParticles();
  if (pm->nprtcl_total != 0) {
    AbortInterfaceRegression("requires <particles>/ppc=0.");
  }

  const Real particle_x = pin->GetReal("problem", "particle_x");
  const Real particle_y = pin->GetOrAddReal("problem", "particle_y", -0.75);
  const Real particle_z = pin->GetOrAddReal("problem", "particle_z", 0.0);
  const Real particle_df_weight =
      pin->GetOrAddReal("problem", "particle_df_weight", 0.0);
  const Real particle_vx = pin->GetOrAddReal("problem", "particle_vx", 0.0);
  const Real particle_vy = pin->GetOrAddReal("problem", "particle_vy", 0.0);
  const Real particle_vz = pin->GetOrAddReal("problem", "particle_vz", 0.0);
  const Real particle_dpxdt = pin->GetOrAddReal("problem", "particle_dpxdt", 0.0);
  const Real particle_dpydt = pin->GetOrAddReal("problem", "particle_dpydt", 0.0);
  const Real particle_dpzdt = pin->GetOrAddReal("problem", "particle_dpzdt", 0.0);
  const Real particle_dedt = pin->GetOrAddReal("problem", "particle_dedt", 0.0);
  const Real particle_ebdot = pin->GetOrAddReal("problem", "particle_ebdot", 0.0);
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
    h_pr(IPVX, 0) = particle_vx;
    h_pr(IPVY, 0) = particle_vy;
    h_pr(IPVZ, 0) = particle_vz;
    h_pr(IPM, 0) = 1.0;
    h_pr(IPWT, 0) = 1.0;
    h_pr(IPF0, 0) = 1.0;
    h_pr(IPDFWT, 0) = 0.0;
    if (ppart->UsesDeltaF()) {
      if (particle_df_weight == 1.0) {
        AbortInterfaceRegression("requires problem/particle_df_weight != 1.");
      }
      const Real f0 = particles::PICDeltaFBackgroundValue(
          ppart->pic_deltaf_background, ppart->pic_deltaf_p0,
          ppart->pic_deltaf_kappa, ppart->pic_deltaf_drift_x1,
          ppart->pic_deltaf_drift_x2, ppart->pic_deltaf_drift_x3,
          ppart->pic_deltaf_aniso_x1, ppart->pic_deltaf_aniso_x2,
          ppart->pic_deltaf_aniso_x3, 1.0, 1.0, 1.0,
          particle_vx, particle_vy, particle_vz);
      h_pr(IPF0, 0) = f0/(1.0 - particle_df_weight);
      h_pr(IPDFWT, 0) = particle_df_weight;
    }
    h_pr(IPDPX, 0) = particle_dpxdt;
    h_pr(IPDPY, 0) = particle_dpydt;
    h_pr(IPDPZ, 0) = particle_dpzdt;
    h_pr(IPDE, 0) = particle_dedt;
    h_pr(IPEBDOT, 0) = particle_ebdot;
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
  if (pin->GetOrAddBoolean("problem", "deposit_manual_records_at_startup", false)) {
    ppart->ZeroMoments(nullptr, 2);
    ppart->DepositMoments(nullptr, 2);
  }
}
