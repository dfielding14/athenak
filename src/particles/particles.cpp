//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particles.cpp
//! \brief implementation of Particles class constructor and assorted other functions

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "particles.hpp"

namespace particles {
//----------------------------------------------------------------------------------------
// constructor, initializes data structures and parameters

Particles::Particles(MeshBlockPack *ppack, ParameterInput *pin) :
    pmy_pack(ppack) {
  // check this is at least a 2D problem
  if (pmy_pack->pmesh->one_d) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Particle module only works in 2D/3D" <<std::endl;
    std::exit(EXIT_FAILURE);
  }

  // read number of particles per cell, and calculate number of particles this pack
  Real ppc = pin->GetOrAddReal("particles","ppc",1.0);
  nspecies = pin->GetOrAddInteger("particles","nspecies",1);
  if (nspecies < 1) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<particles>/nspecies must be >= 1" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string evolution_t = pin->GetOrAddString("time","evolution","dynamic");
  is_dynamic = (evolution_t.compare("static") != 0);

  // compute number of particles as real number, since ppc can be < 1
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int ncells = indcs.nx1*indcs.nx2*indcs.nx3;
  Real r_npart = ppc*static_cast<Real>((pmy_pack->nmb_thispack)*ncells);
  // then cast to integer
  nprtcl_perspec_thispack = static_cast<int>(r_npart);
  nprtcl_thispack = nprtcl_perspec_thispack*nspecies;

  prtcl_rst_flag = pin->GetOrAddInteger("problem","prtcl_rst_flag",0);
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
    nprtcl_thispack = static_cast<int>(gen_data[2]);
    nprtcl_perspec_thispack = nprtcl_thispack/nspecies;
    std::fclose(pfile);
  }

  // select particle type
  {
    std::string ptype = pin->GetString("particles","particle_type");
    if (ptype.compare("cosmic_ray") == 0) {
      particle_type = ParticleType::cosmic_ray;
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle type = '" << ptype << "' not recognized"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // select pusher algorithm
  {
    std::string ppush = pin->GetString("particles","pusher");
    std::string interp = pin->GetOrAddString("particles","interpolation", "tsc");
    if (ppush.compare("drift") == 0) {
      pusher = ParticlesPusher::drift;
    } else if (ppush.compare("boris") == 0) {
      if (pmy_pack->pmhd == nullptr) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Particle pusher 'boris' requires an <mhd> block"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (interp.compare("tsc") == 0) {
        pusher = ParticlesPusher::boris_tsc;
      } else if (interp.compare("lin") == 0) {
        pusher = ParticlesPusher::boris_lin;
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Interpolation method = '" << interp
                  << "' not recognized" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle pusher must be specified in <particles> block"
                <<std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // set dimensions of particle arrays. Note particles only work in 2D/3D
  if (pmy_pack->pmesh->one_d) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Particles only work in 2D/3D, but 1D problem initialized" <<std::endl;
    std::exit(EXIT_FAILURE);
  }
  switch (particle_type) {
    case ParticleType::cosmic_ray:
      {
        nrdata = 14;
        nidata = 3;
        break;
      }
    default:
      break;
  }
  Kokkos::realloc(prtcl_rdata, nrdata, nprtcl_thispack);
  Kokkos::realloc(prtcl_idata, nidata, nprtcl_thispack);

  // allocate boundary object
  pbval_part = new ParticlesBoundaryValues(this, pin);
}

//----------------------------------------------------------------------------------------
// destructor

Particles::~Particles() {
}

//----------------------------------------------------------------------------------------
// CreateParticleTags()
// Assigns tags to particles (unique integer).  Note that tracked particles are always
// those with tag numbers less than ntrack.

void Particles::CreateParticleTags(ParameterInput *pin) {
  if (prtcl_rst_flag) {return;}

  std::string assign = pin->GetOrAddString("particles","assign_tag","index_order");
  int nprtcl_total_perspec = pmy_pack->pmesh->nprtcl_total/nspecies;

  // tags are assigned sequentially within this rank, starting at 0 with rank=0
  if (assign.compare("index_order") == 0) {
    int tagstart = 0;
    for (int n=1; n<=global_variable::my_rank; ++n) {
      tagstart += pmy_pack->pmesh->nprtcl_eachrank[n-1]/nspecies;
    }

    auto &pi = prtcl_idata;
    int nper = nprtcl_perspec_thispack;
    par_for("ptags",DevExeSpace(),0,(nprtcl_thispack-1),
    KOKKOS_LAMBDA(const int p) {
      int spec = (nper > 0) ? p/nper : 0;
      int p_in_spec = (nper > 0) ? p - spec*nper : p;
      pi(PTAG,p) = spec*nprtcl_total_perspec + tagstart + p_in_spec;
    });

  // tags are assigned sequentially across ranks
  } else if (assign.compare("rank_order") == 0) {
    int myrank = global_variable::my_rank;
    int nranks = global_variable::nranks;
    auto &pi = prtcl_idata;
    int nper = nprtcl_perspec_thispack;
    par_for("ptags",DevExeSpace(),0,(nprtcl_thispack-1),
    KOKKOS_LAMBDA(const int p) {
      int spec = (nper > 0) ? p/nper : 0;
      int p_in_spec = (nper > 0) ? p - spec*nper : p;
      pi(PTAG,p) = spec*nprtcl_total_perspec + myrank + nranks*p_in_spec;
    });

  // tag algorithm not recognized, so quit with error
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Particle tag assignment type = '" << assign << "' not recognized"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

namespace {

void WrapCoordinate(Real &x, Real xmin, Real xmax) {
  Real len = xmax - xmin;
  if (len <= 0.0) {return;}
  while (x < xmin) {x += len;}
  while (x > xmax) {x -= len;}
}

void LogicalLocationToRegion(Mesh *pm, const LogicalLocation &lloc, RegionSize &rs) {
  RegionSize &ms = pm->mesh_size;
  int lev_offset = lloc.level - pm->root_level;

  std::int32_t nmbx1 = pm->nmb_rootx1 << lev_offset;
  rs.x1min = (lloc.lx1 == 0) ? ms.x1min : LeftEdgeX(lloc.lx1, nmbx1,
                                                     ms.x1min, ms.x1max);
  rs.x1max = (lloc.lx1 == nmbx1 - 1) ? ms.x1max : LeftEdgeX(lloc.lx1 + 1, nmbx1,
                                                            ms.x1min, ms.x1max);

  if (pm->multi_d) {
    std::int32_t nmbx2 = pm->nmb_rootx2 << lev_offset;
    rs.x2min = (lloc.lx2 == 0) ? ms.x2min : LeftEdgeX(lloc.lx2, nmbx2,
                                                       ms.x2min, ms.x2max);
    rs.x2max = (lloc.lx2 == nmbx2 - 1) ? ms.x2max : LeftEdgeX(lloc.lx2 + 1, nmbx2,
                                                              ms.x2min, ms.x2max);
  } else {
    rs.x2min = ms.x2min;
    rs.x2max = ms.x2max;
  }

  if (pm->three_d) {
    std::int32_t nmbx3 = pm->nmb_rootx3 << lev_offset;
    rs.x3min = (lloc.lx3 == 0) ? ms.x3min : LeftEdgeX(lloc.lx3, nmbx3,
                                                       ms.x3min, ms.x3max);
    rs.x3max = (lloc.lx3 == nmbx3 - 1) ? ms.x3max : LeftEdgeX(lloc.lx3 + 1, nmbx3,
                                                              ms.x3min, ms.x3max);
  } else {
    rs.x3min = ms.x3min;
    rs.x3max = ms.x3max;
  }
}

bool ContainsParticle(Mesh *pm, const RegionSize &rs, const LogicalLocation &lloc,
                      Real x1, Real x2, Real x3) {
  int lev_offset = lloc.level - pm->root_level;
  std::int32_t nmbx1 = pm->nmb_rootx1 << lev_offset;
  std::int32_t nmbx2 = pm->nmb_rootx2 << lev_offset;
  std::int32_t nmbx3 = pm->nmb_rootx3 << lev_offset;

  bool in_x1 = (x1 >= rs.x1min) && (x1 < rs.x1max ||
               (lloc.lx1 == nmbx1 - 1 && x1 <= rs.x1max));
  bool in_x2 = (!pm->multi_d) || ((x2 >= rs.x2min) && (x2 < rs.x2max ||
               (lloc.lx2 == nmbx2 - 1 && x2 <= rs.x2max)));
  bool in_x3 = (!pm->three_d) || ((x3 >= rs.x3min) && (x3 < rs.x3max ||
               (lloc.lx3 == nmbx3 - 1 && x3 <= rs.x3max)));
  return in_x1 && in_x2 && in_x3;
}

int FindContainingMeshBlock(Mesh *pm, Real x1, Real x2, Real x3) {
  RegionSize rs;
  for (int gid=0; gid<pm->nmb_total; ++gid) {
    LogicalLocation &lloc = pm->lloc_eachmb[gid];
    LogicalLocationToRegion(pm, lloc, rs);
    if (ContainsParticle(pm, rs, lloc, x1, x2, x3)) {return gid;}
  }
  return -1;
}

} // namespace

//----------------------------------------------------------------------------------------
// RemapAfterAMR()
// AMR changes the leaf MeshBlock list and may move blocks between ranks.  Particles are
// remapped by position onto the new leaf list, then the existing particle MPI path moves
// off-rank particles.  This is intentionally separate from face-centered AMR repair.

void Particles::RemapAfterAMR() {
  Mesh *pm = pmy_pack->pmesh;
  if (nprtcl_thispack <= 0) {
    pm->UpdateParticleCounts();
    return;
  }

  auto pr_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_rdata);
  auto pi_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_idata);

  ParticlesBoundaryValues *pb = pbval_part;
  int send_capacity = std::max(1, nprtcl_thispack);
  Kokkos::realloc(pb->sendlist, send_capacity);
  pb->nprtcl_send = 0;
  pb->nprtcl_recv = 0;

  int myrank = global_variable::my_rank;
  for (int p=0; p<nprtcl_thispack; ++p) {
    WrapCoordinate(pr_h(IPX,p), pm->mesh_size.x1min, pm->mesh_size.x1max);
    if (pm->multi_d) {
      WrapCoordinate(pr_h(IPY,p), pm->mesh_size.x2min, pm->mesh_size.x2max);
    }
    if (pm->three_d) {
      WrapCoordinate(pr_h(IPZ,p), pm->mesh_size.x3min, pm->mesh_size.x3max);
    }

    int dest_gid = FindContainingMeshBlock(pm, pr_h(IPX,p), pr_h(IPY,p), pr_h(IPZ,p));
    if (dest_gid < 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Could not remap particle " << p
                << " after AMR; position=(" << pr_h(IPX,p) << ", "
                << pr_h(IPY,p) << ", " << pr_h(IPZ,p) << ")" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    pi_h(PGID,p) = dest_gid;
#if MPI_PARALLEL_ENABLED
    int dest_rank = pm->rank_eachmb[dest_gid];
    if (dest_rank != myrank) {
      int index = pb->nprtcl_send++;
      pb->sendlist.h_view(index).prtcl_indx = p;
      pb->sendlist.h_view(index).dest_gid = dest_gid;
      pb->sendlist.h_view(index).dest_rank = dest_rank;
    }
#endif
  }

  Kokkos::deep_copy(prtcl_rdata, pr_h);
  Kokkos::deep_copy(prtcl_idata, pi_h);
  Kokkos::resize(pb->sendlist, pb->nprtcl_send);
  pb->sendlist.template modify<HostMemSpace>();
  pb->sendlist.template sync<DevExeSpace>();

#if MPI_PARALLEL_ENABLED
  pb->CountSendsAndRecvs();
  pb->InitPrtclRecv();
  pb->PackAndSendPrtcls();
  while (pb->RecvAndUnpackPrtcls() == TaskStatus::incomplete) {}
  pb->ClearPrtclRecv();
  pb->ClearPrtclSend();
#else
  pm->UpdateParticleCounts();
#endif
}

} // namespace particles
