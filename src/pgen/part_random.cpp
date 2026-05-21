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
#include <sstream>
#include <iostream>

#include "globals.hpp"
#include "parameter_input.hpp"
#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "particles/particles.hpp"

#include <Kokkos_Random.hpp>

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::UserProblem_()
//! \brief Problem Generator for random particle positions/velocities

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
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
    return;
  }

  Real min_mass = pin->GetOrAddReal("particles","min_mass",1.0);
  Real mass_log_spacing = pin->GetOrAddReal("particles","mass_log_spacing",1.0);
  Real B0x = pin->GetOrAddReal("problem","B0x",0.0);
  Real B0y = pin->GetOrAddReal("problem","B0y",0.0);
  Real B0z = pin->GetOrAddReal("problem","B0z",1.0);
  Real cfl_part = pin->GetOrAddReal("particles","cfl_part",0.05);

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
  Real bmag = std::sqrt(B0x*B0x + B0y*B0y + B0z*B0z);
  if (bmag > 0.0) {
    dtnew_ = std::min(dtnew_, cfl_part*min_mass/bmag);
  }

  return;
}
