//========================================================================================
// Athena++ astrophysical MHD code, Kokkos version
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particle_random.cpp
//! \brief Problem generator that initializes random particle positions and velocities.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>

#include "globals.hpp"
#include "parameter_input.hpp"
#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "particles/particles.hpp"
#include "pgen.hpp"

#include <Kokkos_Random.hpp>

namespace {

KOKKOS_INLINE_FUNCTION
Real WrapUnit(Real x) {
  while (x < 0.0) {x += 1.0;}
  while (x >= 1.0) {x -= 1.0;}
  return x;
}

KOKKOS_INLINE_FUNCTION
Real CenterUnit(Real xmin, Real xmax, Real mesh_min, Real mesh_max) {
  return (0.5*(xmin + xmax) - mesh_min)/(mesh_max - mesh_min);
}

KOKKOS_INLINE_FUNCTION
Real HalfWidthUnit(Real xmin, Real xmax, Real mesh_min, Real mesh_max) {
  return 0.5*(xmax - xmin)/(mesh_max - mesh_min);
}

KOKKOS_INLINE_FUNCTION
Real PeriodicDistance(Real a, Real b) {
  Real d = fabs(a - b);
  return fmin(d, 1.0 - d);
}

KOKKOS_INLINE_FUNCTION
bool OverlapsMovingBox(const RegionSize &sz, const RegionSize &mesh_size,
                       Real c1, Real c2, Real c3, bool multi_d, bool three_d) {
  bool overlaps = PeriodicDistance(
      CenterUnit(sz.x1min, sz.x1max, mesh_size.x1min, mesh_size.x1max), c1)
      <= HalfWidthUnit(sz.x1min, sz.x1max, mesh_size.x1min, mesh_size.x1max) + 0.17;
  if (multi_d) {
    overlaps = overlaps && PeriodicDistance(
        CenterUnit(sz.x2min, sz.x2max, mesh_size.x2min, mesh_size.x2max), c2)
        <= HalfWidthUnit(sz.x2min, sz.x2max, mesh_size.x2min, mesh_size.x2max) + 0.13;
  }
  if (three_d) {
    overlaps = overlaps && PeriodicDistance(
        CenterUnit(sz.x3min, sz.x3max, mesh_size.x3min, mesh_size.x3max), c3)
        <= HalfWidthUnit(sz.x3min, sz.x3max, mesh_size.x3min, mesh_size.x3max) + 0.11;
  }
  return overlaps;
}

void PartRandomRefinementCondition(MeshBlockPack *pmbp) {
  Mesh *pmesh = pmbp->pmesh;
  auto &refine_flag = pmesh->pmr->refine_flag;
  auto &mblev = pmbp->pmb->mb_lev;
  auto &mb_size = pmbp->pmb->mb_size;
  const int nmb = pmbp->nmb_thispack;
  const int mbs = pmesh->gids_eachrank[global_variable::my_rank];
  const int root_level = pmesh->root_level;
  const int target_level = pmesh->max_level;
  const bool multi_d = pmesh->multi_d;
  const bool three_d = pmesh->three_d;
  const RegionSize mesh_size = pmesh->mesh_size;
  const Real phase = static_cast<Real>(pmesh->ncycle);

  par_for("part_random_refinement", DevExeSpace(), 0, nmb-1,
  KOKKOS_LAMBDA(int m) {
    bool refine_region = OverlapsMovingBox(
        mb_size.d_view(m), mesh_size, WrapUnit(0.21 + 0.19*phase),
        WrapUnit(0.37 + 0.11*phase), WrapUnit(0.53 + 0.07*phase),
        multi_d, three_d);
    refine_region = refine_region || OverlapsMovingBox(
        mb_size.d_view(m), mesh_size, WrapUnit(0.74 - 0.13*phase),
        WrapUnit(0.63 + 0.17*phase), WrapUnit(0.28 - 0.09*phase),
        multi_d, three_d);
    const int level = mblev.d_view(m);
    if (refine_region && level < target_level) {
      refine_flag.d_view(m + mbs) = 1;
    } else if ((!refine_region) && level > root_level) {
      refine_flag.d_view(m + mbs) = -1;
    }
  });

  refine_flag.template modify<DevExeSpace>();
  refine_flag.template sync<HostMemSpace>();
}

} // namespace

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::PartRandom()
//! \brief Problem Generator for random particle positions/velocities

void ProblemGenerator::PartRandom(ParameterInput *pin, const bool restart) {
  std::string particle_refinement =
      pin->GetOrAddString("problem","particle_refinement","none");
  if (particle_refinement.compare("moving_boxes") == 0) {
    user_ref_func = PartRandomRefinementCondition;
  } else if (particle_refinement.compare("none") != 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Unknown particle_refinement = '"
              << particle_refinement << "'" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->ppart == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Random particles test requires <particles> block in input file"
              << std::endl;
    exit(EXIT_FAILURE);
  }

  int prtcl_rst_flag = pmbp->ppart->prtcl_rst_flag;
  if (restart && !prtcl_rst_flag) return;

  // capture variables for the kernel
  auto &mbsize = pmbp->pmb->mb_size;
  auto &pr = pmbp->ppart->prtcl_rdata;
  auto &pi = pmbp->ppart->prtcl_idata;
  auto &npart = pmbp->ppart->nprtcl_thispack;
  auto &npart_spec = pmbp->ppart->nprtcl_perspec_thispack;
  auto &nspecies = pmbp->ppart->nspecies;
  auto gids = pmbp->gids;
  auto gide = pmbp->gide;
  auto nmb_thispack = pmbp->nmb_thispack;

  if (prtcl_rst_flag) {
    std::string prst_fname = pin->GetString("problem","prtcl_res_file");
    char rank_dir[20];
    std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d", global_variable::my_rank);
    std::size_t rank_pos = prst_fname.find("rank_00000000");
    if (rank_pos != std::string::npos) {
      prst_fname.replace(rank_pos, sizeof("rank_00000000") - 1, rank_dir);
    }

    FILE* pfile = std::fopen(prst_fname.c_str(),"rb");
    if (pfile == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle restart file '" << prst_fname
                << "' could not be opened" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    Real gen_data[3];
    if (std::fread(gen_data, sizeof(Real), 3, pfile) != 3) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle restart file '" << prst_fname
                << "' has an incomplete header" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    pmbp->ppart->dtnew = gen_data[1];
    pmy_mesh_->time = gen_data[0];
    HostArray1D<Real> host_part_data("host_part_data", 17*npart);
    if (std::fread(host_part_data.data(), sizeof(Real), 17*npart, pfile) !=
        static_cast<std::size_t>(17*npart)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle restart file '" << prst_fname
                << "' has incomplete particle data" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::fclose(pfile);

    DvceArray1D<Real> part_data("part_data", 17*npart);
    Kokkos::deep_copy(part_data, host_part_data);
    par_for("part_restart",DevExeSpace(),0,(npart-1),
    KOKKOS_LAMBDA(const int p) {
      pi(PGID,p) = static_cast<int>(part_data(17*p+0));
      pi(PTAG,p) = static_cast<int>(part_data(17*p+1));
      pi(PSP,p) = static_cast<int>(part_data(17*p+2));
      pr(IPX,p) = part_data(17*p+3);
      pr(IPY,p) = part_data(17*p+4);
      pr(IPZ,p) = part_data(17*p+5);
      pr(IPVX,p) = part_data(17*p+6);
      pr(IPVY,p) = part_data(17*p+7);
      pr(IPVZ,p) = part_data(17*p+8);
      pr(IPM,p) = part_data(17*p+9);
      pr(IPBX,p) = part_data(17*p+10);
      pr(IPBY,p) = part_data(17*p+11);
      pr(IPBZ,p) = part_data(17*p+12);
      pr(IPDX,p) = part_data(17*p+13);
      pr(IPDY,p) = part_data(17*p+14);
      pr(IPDZ,p) = part_data(17*p+15);
      pr(IPDB,p) = part_data(17*p+16);
    });
    if (restart) {
      pmbp->ppart->CheckConsistency("particle restart load");
      return;
    }
  } else {
    Real min_mass = pin->GetOrAddReal("particles","min_mass",1.0);
    Real mass_log_spacing = pin->GetOrAddReal("particles","mass_log_spacing",1.0);
    Real B0x = pin->GetOrAddReal("problem","B0x",0.0);
    Real B0y = pin->GetOrAddReal("problem","B0y",0.0);
    Real B0z = pin->GetOrAddReal("problem","B0z",1.0);

    // initialize particles
    Kokkos::Random_XorShift64_Pool<> rand_pool64(pmbp->gids);
    par_for("part_update",DevExeSpace(),0,(npart-1),
    KOKKOS_LAMBDA(const int p) {
      auto rand_gen = rand_pool64.get_state();  // get random number state this thread
      // choose parent MeshBlock randomly
      int m = static_cast<int>(rand_gen.frand()*(gide - gids + 1.0));
      m = (m < 0) ? 0 : m;
      m = (m > nmb_thispack - 1) ? nmb_thispack - 1 : m;
      pi(PGID,p) = gids + m;
      int spec = (npart_spec > 0) ? p/npart_spec : 0;
      spec = (spec < 0) ? 0 : spec;
      spec = (spec > nspecies - 1) ? nspecies - 1 : spec;
      pi(PSP,p) = spec;

      Real rand = rand_gen.frand();
      pr(IPX,p) = (1. - rand)*mbsize.d_view(m).x1min + rand*mbsize.d_view(m).x1max;
      pr(IPX,p) = fmin(pr(IPX,p),mbsize.d_view(m).x1max);
      pr(IPX,p) = fmax(pr(IPX,p),mbsize.d_view(m).x1min);

      rand = rand_gen.frand();
      pr(IPY,p) = (1. - rand)*mbsize.d_view(m).x2min + rand*mbsize.d_view(m).x2max;
      pr(IPY,p) = fmin(pr(IPY,p),mbsize.d_view(m).x2max);
      pr(IPY,p) = fmax(pr(IPY,p),mbsize.d_view(m).x2min);

      rand = rand_gen.frand();
      pr(IPZ,p) = (1. - rand)*mbsize.d_view(m).x3min + rand*mbsize.d_view(m).x3max;
      pr(IPZ,p) = fmin(pr(IPZ,p),mbsize.d_view(m).x3max);
      pr(IPZ,p) = fmax(pr(IPZ,p),mbsize.d_view(m).x3min);

      pr(IPVX,p) = 2.0*(rand_gen.frand() - 0.5);
      pr(IPVY,p) = 2.0*(rand_gen.frand() - 0.5);
      pr(IPVZ,p) = 2.0*(rand_gen.frand() - 0.5);
      pr(IPM,p) = min_mass*pow(mass_log_spacing, spec);
      pr(IPBX,p) = B0x;
      pr(IPBY,p) = B0y;
      pr(IPBZ,p) = B0z;
      pr(IPDX,p) = 0.0;
      pr(IPDY,p) = 0.0;
      pr(IPDZ,p) = 0.0;
      pr(IPDB,p) = 0.0;

      rand_pool64.free_state(rand_gen);  // free state for use by other threads
    });
  }

  Real B0x = pin->GetOrAddReal("problem","B0x",0.0);
  Real B0y = pin->GetOrAddReal("problem","B0y",0.0);
  Real B0z = pin->GetOrAddReal("problem","B0z",1.0);
  Real cfl_part = pin->GetOrAddReal("particles","cfl_part",0.05);

  if (pmbp->pmhd != nullptr) {
    EOS_Data &eos = pmbp->pmhd->peos->eos_data;
    Real gm1 = eos.gamma - 1.0;
    Real p0 = 1.0/eos.gamma;
    auto &u0 = pmbp->pmhd->u0;
    auto &b0 = pmbp->pmhd->b0;
    auto &bcc0 = pmbp->pmhd->bcc0;
    auto &indcs = pmy_mesh_->mb_indcs;
    int &is = indcs.is, &ie = indcs.ie;
    int &js = indcs.js, &je = indcs.je;
    int &ks = indcs.ks, &ke = indcs.ke;
    bool is_ideal = eos.is_ideal;
    par_for("pgen_mhd", DevExeSpace(), 0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0(m,IDN,k,j,i) = 1.0;
      u0(m,IM1,k,j,i) = 0.0;
      u0(m,IM2,k,j,i) = 0.0;
      u0(m,IM3,k,j,i) = 0.0;
      bcc0(m,IBX,k,j,i) = B0x;
      bcc0(m,IBY,k,j,i) = B0y;
      bcc0(m,IBZ,k,j,i) = B0z;
      b0.x1f(m,k,j,i) = B0x;
      b0.x2f(m,k,j,i) = B0y;
      b0.x3f(m,k,j,i) = B0z;
      if (i==ie) {b0.x1f(m,k,j,i+1) = B0x;}
      if (j==je) {b0.x2f(m,k,j+1,i) = B0y;}
      if (k==ke) {b0.x3f(m,k+1,j,i) = B0z;}
      if (is_ideal) {
        u0(m,IEN,k,j,i) = p0/gm1 + 0.5*(B0x*B0x + B0y*B0y + B0z*B0z);
      }
    });
  }

  // set timestep (which will remain constant for entire run
  // Assumes uniform mesh (no SMR or AMR)
  // Assumes velocities normalized to one, so dt=min(dx)
  Real &dtnew_ = pmbp->ppart->dtnew;
  dtnew_ = std::min(mbsize.h_view(0).dx1, mbsize.h_view(0).dx2);
  dtnew_ = std::min(dtnew_, mbsize.h_view(0).dx3);
  Real min_mass = pin->GetOrAddReal("particles","min_mass",1.0);
  Real bmag = std::sqrt(B0x*B0x + B0y*B0y + B0z*B0z);
  if (bmag > 0.0) {
    dtnew_ = std::min(dtnew_, cfl_part*min_mass/bmag);
  }

  return;
}

#if PART_RANDOM_USER_PROBLEM_ENABLED
void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  PartRandom(pin, restart);
}
#endif
