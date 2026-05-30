//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particles.cpp
//! \brief implementation of Particles class constructor and assorted other functions

#include <iostream>
#include <limits>
#include <string>
#include <algorithm>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "particles.hpp"

namespace particles {

namespace {

[[noreturn]] void ParticleInputFatal(const char *file, const int line,
                                     const std::string &field,
                                     const std::string &reason) {
  std::cout << "### FATAL ERROR in " << file << " at line " << line << std::endl
            << "Invalid input " << field << ": " << reason << std::endl;
  std::exit(EXIT_FAILURE);
}

DragParticlesCoupling ParseDragCoupling(ParameterInput *pin, const std::string &field,
                                        const std::string &default_value) {
  std::string value = pin->GetOrAddString("drag_particles", field, default_value);
  if (value.compare("host_cell") == 0 || value.compare("nearest") == 0) {
    return DragParticlesCoupling::host_cell;
  }
  if (value.compare("cloud_in_cell") == 0 || value.compare("cic") == 0) {
    return DragParticlesCoupling::cloud_in_cell;
  }
  ParticleInputFatal(__FILE__, __LINE__, "<drag_particles>/" + field,
                     "expected host_cell or cloud_in_cell");
}

} // namespace

//----------------------------------------------------------------------------------------
// constructor, initializes data structures and parameters

Particles::Particles(MeshBlockPack *ppack, ParameterInput *pin) :
    pmy_pack(ppack),
    dtnew(std::numeric_limits<Real>::max()),
    dtnew_limit(std::numeric_limits<Real>::max()),
    cfl_part(0.5),
    drag_model(DragParticlesModel::none),
    drag_enabled(false),
    drag_backreaction(true),
    drag_include_energy(true),
    drag_orbital_terms(false),
    drag_interpolation(DragParticlesCoupling::host_cell),
    drag_deposition(DragParticlesCoupling::host_cell),
    drag_stopping_time(1.0),
    drag_particle_mass(1.0),
    drag_omega0(1.0),
    drag_qshear(1.5) {
  // check this is at least a 2D problem
  if (pmy_pack->pmesh->one_d) {
    ParticleInputFatal(__FILE__, __LINE__, "<mesh>/nx2,nx3",
                       "particle module only works in 2D/3D");
  }

  // read number of particles per cell, and calculate number of particles this pack
  Real ppc = pin->GetOrAddReal("particles","ppc",1.0);
  if (ppc < 0.0) {
    ParticleInputFatal(__FILE__, __LINE__, "<particles>/ppc", "must be >= 0");
  }
  cfl_part = pin->GetOrAddReal("particles","cfl_part",0.5);
  if (cfl_part <= 0.0 || cfl_part > 1.0) {
    ParticleInputFatal(__FILE__, __LINE__, "<particles>/cfl_part",
                       "must be > 0 and <= 1");
  }

  // compute number of particles as real number, since ppc can be < 1
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int ncells = indcs.nx1*indcs.nx2*indcs.nx3;
  Real r_npart = ppc*static_cast<Real>((pmy_pack->nmb_thispack)*ncells);
  // then cast to integer
  nprtcl_thispack = static_cast<int>(r_npart);

  // select particle type
  {
    std::string ptype = pin->GetString("particles","particle_type");
    if (ptype.compare("cosmic_ray") == 0) {
      particle_type = ParticleType::cosmic_ray;
    } else if (ptype.compare("dust") == 0) {
      particle_type = ParticleType::dust;
    } else {
      ParticleInputFatal(__FILE__, __LINE__, "<particles>/particle_type",
                         "value '" + ptype + "' not recognized");
    }
  }

  // Configure optional drag coupling.  The model lives in a separate block so pgens
  // can select the standard stopping-time form or enroll a fully custom callback.
  bool has_drag_block = pin->DoesBlockExist("drag_particles");
  if (has_drag_block) {
    drag_enabled = pin->GetOrAddBoolean("drag_particles","enabled",true);
    std::string dmodel = pin->GetOrAddString("drag_particles","model","stopping_time");
    Real cfl_drag = pin->GetOrAddReal("drag_particles","cfl_drag",0.25);
    if (cfl_drag <= 0.0) {
      ParticleInputFatal(__FILE__, __LINE__, "<drag_particles>/cfl_drag",
                         "must be > 0");
    }
    if (dmodel.compare("stopping_time") == 0) {
      drag_model = DragParticlesModel::stopping_time;
      drag_stopping_time = pin->GetReal("drag_particles","stopping_time");
      if (drag_stopping_time <= 0.0) {
        ParticleInputFatal(__FILE__, __LINE__, "<drag_particles>/stopping_time",
                           "must be > 0");
      }
      dtnew_limit = std::min(dtnew_limit, cfl_drag*drag_stopping_time);
    } else if (dmodel.compare("user") == 0) {
      drag_model = DragParticlesModel::user;
      if (pin->DoesParameterExist("drag_particles","stopping_time")) {
        drag_stopping_time = pin->GetReal("drag_particles","stopping_time");
        if (drag_stopping_time <= 0.0) {
          ParticleInputFatal(__FILE__, __LINE__, "<drag_particles>/stopping_time",
                             "must be > 0 when provided for model=user");
        }
        dtnew_limit = std::min(dtnew_limit, cfl_drag*drag_stopping_time);
      }
    } else if (dmodel.compare("none") == 0) {
      drag_model = DragParticlesModel::none;
      drag_enabled = false;
    } else {
      ParticleInputFatal(__FILE__, __LINE__, "<drag_particles>/model",
                         "value '" + dmodel + "' not recognized");
    }
    drag_backreaction = pin->GetOrAddBoolean("drag_particles","backreaction",true);
    drag_include_energy = pin->GetOrAddBoolean("drag_particles","include_energy",true);
    drag_interpolation = ParseDragCoupling(pin, "interpolation", "host_cell");
    drag_deposition = ParseDragCoupling(pin, "deposition", "host_cell");
    drag_particle_mass = pin->GetOrAddReal("drag_particles","particle_mass",1.0);
    if (drag_particle_mass <= 0.0) {
      ParticleInputFatal(__FILE__, __LINE__, "<drag_particles>/particle_mass",
                         "must be > 0");
    }
    drag_orbital_terms = pin->GetOrAddBoolean("drag_particles","orbital_terms",false);
    if (drag_orbital_terms && pin->DoesBlockExist("shearing_box")) {
      drag_omega0 = pin->GetOrAddReal("shearing_box","omega0",1.0);
      drag_qshear = pin->GetOrAddReal("shearing_box","qshear",1.5);
    }
    if (drag_orbital_terms) {
      drag_omega0 = pin->GetOrAddReal("drag_particles","omega0",drag_omega0);
      drag_qshear = pin->GetOrAddReal("drag_particles","qshear",drag_qshear);
    }
  }

  // select pusher algorithm
  {
    std::string ppush = pin->GetString("particles","pusher");
    if (ppush.compare("drift") == 0) {
      pusher = ParticlesPusher::drift;
    } else if (ppush.compare("drag") == 0) {
      pusher = ParticlesPusher::drag;
      if (!has_drag_block) {
        ParticleInputFatal(__FILE__, __LINE__, "<particles>/pusher",
                           "drag requires a <drag_particles> input block");
      }
      if ((!drag_enabled) || (drag_model == DragParticlesModel::none)) {
        ParticleInputFatal(__FILE__, __LINE__, "<particles>/pusher",
                           "drag requires an enabled <drag_particles> model");
      }
    } else {
      ParticleInputFatal(__FILE__, __LINE__, "<particles>/pusher",
                         "expected drift or drag");
    }
  }

  // set dimensions of particle arrays. Note particles only work in 2D/3D
  if (pmy_pack->pmesh->one_d) {
    ParticleInputFatal(__FILE__, __LINE__, "<mesh>/nx2,nx3",
                       "particle module only works in 2D/3D");
  }
  switch (particle_type) {
    case ParticleType::cosmic_ray:
      {
        int ndim=4;
        if (pmy_pack->pmesh->three_d) {ndim+=2;}
        nrdata = ndim;
        nidata = 2;
        break;
      }
    case ParticleType::dust:
      {
        nrdata = 7;
        nidata = 2;
        break;
      }
    default:
      break;
  }
  if (drag_enabled) {
    if (particle_type != ParticleType::dust) {
      ParticleInputFatal(__FILE__, __LINE__, "<particles>/particle_type",
                         "drag-coupled particles require dust");
    }
    if (pusher != ParticlesPusher::drag) {
      ParticleInputFatal(__FILE__, __LINE__, "<drag_particles>/enabled",
                         "true requires <particles>/pusher=drag");
    }
    bool has_hydro = (pmy_pack->phydro != nullptr);
    bool has_mhd = (pmy_pack->pmhd != nullptr);
    if (!has_hydro && !has_mhd) {
      ParticleInputFatal(__FILE__, __LINE__, "<drag_particles>/enabled",
                         "true requires exactly one host fluid block: <hydro> or <mhd>");
    }
    if (has_hydro && has_mhd) {
      ParticleInputFatal(__FILE__, __LINE__, "<drag_particles>/enabled",
                         "true currently supports exactly one host fluid block");
    }
    if (!pmy_pack->pmesh->strictly_periodic) {
      ParticleInputFatal(__FILE__, __LINE__, "<mesh>/*_bc",
                         "drag particles currently require periodic boundaries");
    }
    nrdata = std::max(nrdata, IPM + 1);
  }
  dtnew = dtnew_limit;
  Kokkos::realloc(prtcl_rdata, nrdata, nprtcl_thispack);
  Kokkos::realloc(prtcl_idata, nidata, nprtcl_thispack);

  // allocate boundary object
  pbval_part = new ParticlesBoundaryValues(this, pin);
}

//----------------------------------------------------------------------------------------
// destructor

Particles::~Particles() {
  delete pbval_part;
}

//----------------------------------------------------------------------------------------
// SetFixedTimeStepLimit()

void Particles::SetFixedTimeStepLimit(const Real dtlimit) {
  dtnew_limit = std::min(dtnew_limit, dtlimit);
  dtnew = std::min(dtnew, dtnew_limit);
}

//----------------------------------------------------------------------------------------
// UpdateNewTimeStep()
// The particle exchange path communicates only with immediate MeshBlock neighbors.  Keep
// each push below one MeshBlock width so a particle cannot skip over its destination.

void Particles::UpdateNewTimeStep() {
  dtnew = dtnew_limit;
  if (nprtcl_thispack <= 0) {return;}

  auto pr = prtcl_rdata;
  auto pi = prtcl_idata;
  auto &mbsize = pmy_pack->pmb->mb_size;
  int gids = pmy_pack->gids;
  int nmb = pmy_pack->nmb_thispack;
  bool multi_d = pmy_pack->pmesh->multi_d;
  bool three_d = pmy_pack->pmesh->three_d;
  Real crossing_time = std::numeric_limits<Real>::max();

  Kokkos::parallel_reduce("particle_newdt", Kokkos::RangePolicy<DevExeSpace>(0,
      nprtcl_thispack), KOKKOS_LAMBDA(const int p, Real &dtmin) {
    int m = pi(PGID,p) - gids;
    if ((m < 0) || (m >= nmb)) {return;}

    Real speed = fabs(pr(IPVX,p));
    if (speed > 0.0) {
      dtmin = fmin(dtmin, (mbsize.d_view(m).x1max - mbsize.d_view(m).x1min)/speed);
    }
    if (multi_d) {
      speed = fabs(pr(IPVY,p));
      if (speed > 0.0) {
        dtmin = fmin(dtmin, (mbsize.d_view(m).x2max - mbsize.d_view(m).x2min)/speed);
      }
    }
    if (three_d) {
      speed = fabs(pr(IPVZ,p));
      if (speed > 0.0) {
        dtmin = fmin(dtmin, (mbsize.d_view(m).x3max - mbsize.d_view(m).x3min)/speed);
      }
    }
  }, Kokkos::Min<Real>(crossing_time));

  dtnew = std::min(dtnew, cfl_part*crossing_time);
}

//----------------------------------------------------------------------------------------
// CreateParticleTags()
// Assigns tags to particles (unique integer).  Note that tracked particles are always
// those with tag numbers less than ntrack.

void Particles::CreateParticleTags(ParameterInput *pin) {
  std::string assign = pin->GetOrAddString("particles","assign_tag","index_order");

  // tags are assigned sequentially within this rank, starting at 0 with rank=0
  if (assign.compare("index_order") == 0) {
    int tagstart = 0;
    for (int n=1; n<=global_variable::my_rank; ++n) {
      tagstart += pmy_pack->pmesh->nprtcl_eachrank[n-1];
    }

    auto &pi = prtcl_idata;
    par_for("ptags",DevExeSpace(),0,(nprtcl_thispack-1),
    KOKKOS_LAMBDA(const int p) {
      pi(PTAG,p) = tagstart + p;
    });

  // tags are assigned sequentially across ranks
  } else if (assign.compare("rank_order") == 0) {
    int myrank = global_variable::my_rank;
    int nranks = global_variable::nranks;
    auto &pi = prtcl_idata;
    par_for("ptags",DevExeSpace(),0,(nprtcl_thispack-1),
    KOKKOS_LAMBDA(const int p) {
      pi(PTAG,p) = myrank + nranks*p;
    });

  // tag algorithm not recognized, so quit with error
  } else {
    ParticleInputFatal(__FILE__, __LINE__, "<particles>/assign_tag",
                       "value '" + assign + "' not recognized");
  }
}

namespace {

void WrapCoordinate(Real &x, Real xmin, Real xmax) {
  Real len = xmax - xmin;
  if (len <= 0.0) {return;}
  while (x < xmin) {x += len;}
  while (x >= xmax) {x -= len;}
}

} // namespace

//----------------------------------------------------------------------------------------
// RemapAfterAMR()
// AMR rebuilds the leaf MeshBlock list and can move MeshBlocks between ranks.  Particles
// are remapped by position onto the new list, then the normal particle MPI path moves any
// off-rank particles.

void Particles::RemapAfterAMR() {
  Mesh *pm = pmy_pack->pmesh;

  ParticlesBoundaryValues *pb = pbval_part;
  int send_capacity = std::max(1, nprtcl_thispack);
  Kokkos::realloc(pb->sendlist, send_capacity);
  pb->nprtcl_send = 0;
  pb->nprtcl_recv = 0;

  int myrank = global_variable::my_rank;
  if (nprtcl_thispack > 0) {
    auto pr_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_rdata);
    auto pi_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_idata);

    for (int p=0; p<nprtcl_thispack; ++p) {
      WrapCoordinate(pr_h(IPX,p), pm->mesh_size.x1min, pm->mesh_size.x1max);
      if (pm->multi_d) {
        WrapCoordinate(pr_h(IPY,p), pm->mesh_size.x2min, pm->mesh_size.x2max);
      }
      if (pm->three_d) {
        WrapCoordinate(pr_h(IPZ,p), pm->mesh_size.x3min, pm->mesh_size.x3max);
      }

      int dest_gid = pm->FindMeshBlockByPosition(pr_h(IPX,p), pr_h(IPY,p), pr_h(IPZ,p));
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
  }
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
