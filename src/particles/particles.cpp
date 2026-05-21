//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particles.cpp
//! \brief implementation of Particles class constructor and assorted other functions

#include <iostream>
#include <string>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <vector>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "hydro/hydro.hpp"
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

  // select particle type
  {
    std::string ptype = pin->GetString("particles","particle_type");
    if (ptype.compare("cosmic_ray") == 0) {
      particle_type = ParticleType::cosmic_ray;
    } else if (ptype.compare("star") == 0) {
      particle_type = ParticleType::star;
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle type = '" << ptype << "' not recognized"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // read number of particles per cell, and calculate number of particles this pack
  if (particle_type == ParticleType::cosmic_ray) {
    Real ppc = pin->GetOrAddReal("particles","ppc",1.0);

    // compute number of particles as real number, since ppc can be < 1
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    int ncells = indcs.nx1*indcs.nx2*indcs.nx3;
    Real r_npart = ppc*static_cast<Real>((pmy_pack->nmb_thispack)*ncells);
    // then cast to integer
    nprtcl_thispack = static_cast<int>(r_npart);
  } else {
    nprtcl_thispack = 0;
  }

  // select pusher algorithm
  {
    std::string ppush = pin->GetString("particles","pusher");
    if (ppush.compare("drift") == 0) {
      pusher = ParticlesPusher::drift;
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
        int ndim=4;
        if (pmy_pack->pmesh->three_d) {ndim+=2;}
        nrdata = ndim;
        nidata = 2;
        break;
      }
    case ParticleType::star:
      {
        if (!pmy_pack->pmesh->three_d) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "Star particles require a 3D mesh" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        nrdata = 8;
        nidata = 2;

        std::string init_mode = pin->GetOrAddString("particles", "star_init", "none");
        if (init_mode.compare("file") == 0 || init_mode.compare("particle_file") == 0) {
          star_init_from_file = true;
        } else if (init_mode.compare("none") != 0) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<particles>/star_init = '" << init_mode
                    << "' not recognized. Valid choices are none or file." << std::endl;
          std::exit(EXIT_FAILURE);
        }

        star_formation_enabled = pin->GetOrAddBoolean("particles", "star_formation",
                                                      false);
        star_accretion_enabled = pin->GetOrAddBoolean("particles", "star_accretion",
                                                      false);
        if ((star_formation_enabled || star_accretion_enabled) &&
            pmy_pack->phydro == nullptr) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "Star formation/accretion requires a <hydro> block."
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }

        star_formation_density_threshold =
            pin->GetOrAddReal("particles", "star_formation_density_threshold", 0.0);
        star_formation_particle_mass =
            pin->GetOrAddReal("particles", "star_formation_particle_mass", 0.0);
        star_formation_density_floor =
            pin->GetOrAddReal("particles", "star_formation_density_floor", 0.0);
        star_formation_interval =
            pin->GetOrAddInteger("particles", "star_formation_interval", 1);
        star_formation_max_per_cycle =
            pin->GetOrAddInteger("particles", "star_formation_max_per_cycle", -1);
        star_remove_gas_on_formation =
            pin->GetOrAddBoolean("particles", "star_remove_gas_on_formation", true);
        star_accretion_rate =
            pin->GetOrAddReal("particles", "star_accretion_rate", 0.0);
        star_accretion_radius_cells =
            pin->GetOrAddInteger("particles", "star_accretion_radius_cells", 0);
        star_accretion_max_fraction =
            pin->GetOrAddReal("particles", "star_accretion_max_fraction", 0.25);
        star_accretion_density_floor =
            pin->GetOrAddReal("particles", "star_accretion_density_floor",
                              star_formation_density_floor);

        if (star_formation_interval < 1 || star_accretion_radius_cells < 0 ||
            star_accretion_max_fraction < 0.0 || star_accretion_max_fraction > 1.0) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "Invalid star particle formation/accretion parameter."
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (star_formation_enabled &&
            (star_formation_density_threshold <= 0.0 ||
             star_formation_particle_mass <= 0.0)) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "Star formation requires positive "
                    << "star_formation_density_threshold and "
                    << "star_formation_particle_mass." << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (star_accretion_enabled && star_accretion_rate < 0.0) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "star_accretion_rate must be non-negative."
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
        break;
      }
    default:
      break;
  }
  Kokkos::realloc(prtcl_rdata, nrdata, nprtcl_thispack);
  Kokkos::realloc(prtcl_idata, nidata, nprtcl_thispack);

  // allocate boundary object
  pbval_part = new ParticlesBoundaryValues(this, pin);

  if (particle_type == ParticleType::star) {
    if (star_init_from_file) {
      LoadStarsFromFile(pin);
    }
    Real default_dt = std::min(pmy_pack->pmb->mb_size.h_view(0).dx1,
                               pmy_pack->pmb->mb_size.h_view(0).dx2);
    default_dt = std::min(default_dt, pmy_pack->pmb->mb_size.h_view(0).dx3);
    dtnew = pin->GetOrAddReal("particles", "dt", default_dt);
  }
}

//----------------------------------------------------------------------------------------
// destructor

Particles::~Particles() {
  delete pbval_part;
}

//----------------------------------------------------------------------------------------
// FormStars()
// Creates star particles in cells above the configured density threshold.

TaskStatus Particles::FormStars(Driver *pdriver, int stage) {
  if (particle_type != ParticleType::star || !star_formation_enabled) {
    return TaskStatus::complete;
  }
  if ((pmy_pack->pmesh->ncycle % star_formation_interval) != 0) {
    return TaskStatus::complete;
  }

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nhydro_total = pmy_pack->phydro->nhydro + pmy_pack->phydro->nscalars;

  auto uhost = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->phydro->u0);
  auto whost = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->phydro->w0);
  auto &size = pmy_pack->pmb->mb_size;

  struct NewStar {
    std::array<Real, 8> rdata;
    int gid;
    int tag;
  };
  std::vector<NewStar> new_stars;
  Real formed_mass = 0.0;

  auto remove_cell_mass = [&](int m, int k, int j, int i, Real requested_mass) {
    Real rho = uhost(m, IDN, k, j, i);
    Real vol = size.h_view(m).dx1*size.h_view(m).dx2*size.h_view(m).dx3;
    Real available = std::max(static_cast<Real>(0.0),
                              (rho - star_formation_density_floor)*vol);
    Real removed = std::min(requested_mass, available);
    if (removed <= 0.0 || rho <= 0.0) {return 0.0;}

    Real frac = std::min(removed/(rho*vol), static_cast<Real>(1.0));
    for (int n=0; n<nhydro_total; ++n) {
      uhost(m,n,k,j,i) *= (1.0 - frac);
    }
    whost(m,IDN,k,j,i) = uhost(m,IDN,k,j,i);
    if (pmy_pack->phydro->nhydro > IEN) {
      whost(m,IEN,k,j,i) *= (1.0 - frac);
    }
    return removed;
  };

  for (int m=0; m<pmy_pack->nmb_thispack; ++m) {
    bool stop = false;
    for (int k=ks; k<=ke && !stop; ++k) {
      for (int j=js; j<=je && !stop; ++j) {
        for (int i=is; i<=ie; ++i) {
          if (uhost(m,IDN,k,j,i) <= star_formation_density_threshold) {continue;}
          Real particle_mass = star_formation_particle_mass;
          if (star_remove_gas_on_formation) {
            particle_mass = remove_cell_mass(m, k, j, i, star_formation_particle_mass);
          }
          if (particle_mass <= 0.0) {continue;}

          Real x = CellCenterX(i - is, indcs.nx1, size.h_view(m).x1min,
                               size.h_view(m).x1max);
          Real y = CellCenterX(j - js, indcs.nx2, size.h_view(m).x2min,
                               size.h_view(m).x2max);
          Real z = CellCenterX(k - ks, indcs.nx3, size.h_view(m).x3min,
                               size.h_view(m).x3max);
          Real rho = whost(m,IDN,k,j,i);
          Real vx = (rho > 0.0) ? whost(m,IVX,k,j,i) : 0.0;
          Real vy = (rho > 0.0) ? whost(m,IVY,k,j,i) : 0.0;
          Real vz = (rho > 0.0) ? whost(m,IVZ,k,j,i) : 0.0;
          new_stars.push_back({{x, vx, y, vy, z, vz, particle_mass,
                                pmy_pack->pmesh->time}, pmy_pack->gids + m,
                                next_star_tag});
          next_star_tag += global_variable::nranks;
          formed_mass += particle_mass;

          if (star_formation_max_per_cycle >= 0 &&
              static_cast<int>(new_stars.size()) >= star_formation_max_per_cycle) {
            stop = true;
            break;
          }
        }
      }
    }
    if (stop) {break;}
  }

  int local_new_stars = static_cast<int>(new_stars.size());
  int global_new_stars = local_new_stars;
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(&local_new_stars, &global_new_stars, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
#endif

  if (global_new_stars == 0) {
    return TaskStatus::complete;
  }
  if (local_new_stars == 0) {
    UpdateParticleCounts();
    return TaskStatus::complete;
  }

  Kokkos::deep_copy(pmy_pack->phydro->u0, uhost);
  Kokkos::deep_copy(pmy_pack->phydro->w0, whost);

  int old_npart = nprtcl_thispack;
  auto old_rdata = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_rdata);
  auto old_idata = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_idata);

  nprtcl_thispack += local_new_stars;
  Kokkos::realloc(prtcl_rdata, nrdata, nprtcl_thispack);
  Kokkos::realloc(prtcl_idata, nidata, nprtcl_thispack);
  auto rhost = Kokkos::create_mirror_view(prtcl_rdata);
  auto ihost = Kokkos::create_mirror_view(prtcl_idata);

  for (int p=0; p<old_npart; ++p) {
    for (int n=0; n<nrdata; ++n) {rhost(n,p) = old_rdata(n,p);}
    for (int n=0; n<nidata; ++n) {ihost(n,p) = old_idata(n,p);}
  }
  for (int n=0; n<local_new_stars; ++n) {
    int p = old_npart + n;
    for (int q=0; q<nrdata; ++q) {rhost(q,p) = new_stars[n].rdata[q];}
    ihost(PGID,p) = new_stars[n].gid;
    ihost(PTAG,p) = new_stars[n].tag;
  }

  Kokkos::deep_copy(prtcl_rdata, rhost);
  Kokkos::deep_copy(prtcl_idata, ihost);
  star_mass_formed_total += formed_mass;
  UpdateParticleCounts();
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
// AccreteStars()
// Grows star particles by removing gas from the host and neighboring cells.

TaskStatus Particles::AccreteStars(Driver *pdriver, int stage) {
  if (particle_type != ParticleType::star || !star_accretion_enabled ||
      nprtcl_thispack == 0 || star_accretion_rate == 0.0) {
    return TaskStatus::complete;
  }

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nhydro_total = pmy_pack->phydro->nhydro + pmy_pack->phydro->nscalars;
  Real dt = pmy_pack->pmesh->dt;
  Real frac_limit = std::min(star_accretion_max_fraction, star_accretion_rate*dt);
  if (frac_limit <= 0.0) {return TaskStatus::complete;}

  auto rhost = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_rdata);
  auto ihost = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_idata);
  auto uhost = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->phydro->u0);
  auto whost = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->phydro->w0);
  auto &size = pmy_pack->pmb->mb_size;

  Real accreted_mass = 0.0;
  for (int p=0; p<nprtcl_thispack; ++p) {
    int m = ihost(PGID,p) - pmy_pack->gids;
    if (m < 0 || m >= pmy_pack->nmb_thispack) {continue;}
    int ic = static_cast<int>((rhost(IPX,p) - size.h_view(m).x1min)/
                              size.h_view(m).dx1) + is;
    int jc = static_cast<int>((rhost(IPY,p) - size.h_view(m).x2min)/
                              size.h_view(m).dx2) + js;
    int kc = static_cast<int>((rhost(IPZ,p) - size.h_view(m).x3min)/
                              size.h_view(m).dx3) + ks;
    ic = std::max(is, std::min(ie, ic));
    jc = std::max(js, std::min(je, jc));
    kc = std::max(ks, std::min(ke, kc));

    Real particle_gain = 0.0;
    for (int dk=-star_accretion_radius_cells; dk<=star_accretion_radius_cells; ++dk) {
      int k = kc + dk;
      if (k < ks || k > ke) {continue;}
      for (int dj=-star_accretion_radius_cells; dj<=star_accretion_radius_cells; ++dj) {
        int j = jc + dj;
        if (j < js || j > je) {continue;}
        for (int di=-star_accretion_radius_cells; di<=star_accretion_radius_cells; ++di) {
          int i = ic + di;
          if (i < is || i > ie) {continue;}
          Real rho = uhost(m,IDN,k,j,i);
          Real available_rho = std::max(static_cast<Real>(0.0),
                                        rho - star_accretion_density_floor);
          Real delta_rho = frac_limit*available_rho;
          if (delta_rho <= 0.0 || rho <= 0.0) {continue;}
          Real remove_frac = std::min(delta_rho/rho, static_cast<Real>(1.0));
          for (int n=0; n<nhydro_total; ++n) {
            uhost(m,n,k,j,i) *= (1.0 - remove_frac);
          }
          whost(m,IDN,k,j,i) = uhost(m,IDN,k,j,i);
          if (pmy_pack->phydro->nhydro > IEN) {
            whost(m,IEN,k,j,i) *= (1.0 - remove_frac);
          }
          Real vol = size.h_view(m).dx1*size.h_view(m).dx2*size.h_view(m).dx3;
          particle_gain += delta_rho*vol;
        }
      }
    }
    rhost(IPMASS,p) += particle_gain;
    accreted_mass += particle_gain;
  }

  if (accreted_mass > 0.0) {
    Kokkos::deep_copy(prtcl_rdata, rhost);
    Kokkos::deep_copy(pmy_pack->phydro->u0, uhost);
    Kokkos::deep_copy(pmy_pack->phydro->w0, whost);
    star_mass_accreted_total += accreted_mass;
  }
  return TaskStatus::complete;
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
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Particle tag assignment type = '" << assign << "' not recognized"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if (particle_type == ParticleType::star) {
    RefreshNextStarTag();
  }
}

//----------------------------------------------------------------------------------------
// LoadStarsFromFile()
// Initializes star particles from an ASCII table with x y z vx vy vz mass [t_create].

void Particles::LoadStarsFromFile(ParameterInput *pin) {
  std::string particle_file;
  if (pin->DoesParameterExist("particles", "particle_file")) {
    particle_file = pin->GetString("particles", "particle_file");
  } else {
    particle_file = pin->GetString("particles", "star_particle_file");
  }

  std::ifstream infile(particle_file);
  if (!infile) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Unable to open star particle file '" << particle_file
              << "'" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::string columns = pin->GetOrAddString("particles", "particle_file_columns",
                                            "mass_then_time");
  if (columns.compare("mass_then_time") != 0 &&
      columns.compare("time_then_mass") != 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<particles>/particle_file_columns = '" << columns
              << "' not recognized. Valid choices are mass_then_time or time_then_mass."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::vector<std::array<Real, 8>> particles;
  auto &size = pmy_pack->pmb->mb_size;
  const int nmb = pmy_pack->nmb_thispack;

  std::string line;
  while (std::getline(infile, line)) {
    if (line.empty() || line[0] == '#') {continue;}
    std::istringstream iss(line);
    Real x, y, z, vx, vy, vz, col7, col8;
    if (!(iss >> x >> y >> z >> vx >> vy >> vz >> col7)) {continue;}
    Real mass = col7;
    Real tcreate = 0.0;
    if (iss >> col8) {
      if (columns.compare("mass_then_time") == 0) {
        mass = col7;
        tcreate = col8;
      } else {
        tcreate = col7;
        mass = col8;
      }
    }

    for (int m=0; m<nmb; ++m) {
      if (x >= size.h_view(m).x1min && x < size.h_view(m).x1max &&
          y >= size.h_view(m).x2min && y < size.h_view(m).x2max &&
          z >= size.h_view(m).x3min && z < size.h_view(m).x3max) {
        particles.push_back({x, vx, y, vy, z, vz, mass, tcreate});
        break;
      }
    }
  }

  nprtcl_thispack = static_cast<int>(particles.size());
  Kokkos::realloc(prtcl_rdata, nrdata, nprtcl_thispack);
  Kokkos::realloc(prtcl_idata, nidata, nprtcl_thispack);

  auto rhost = Kokkos::create_mirror_view(prtcl_rdata);
  auto ihost = Kokkos::create_mirror_view(prtcl_idata);
  for (int p=0; p<nprtcl_thispack; ++p) {
    int gid = pmy_pack->gids;
    for (int m=0; m<nmb; ++m) {
      if (particles[p][IPX] >= size.h_view(m).x1min &&
          particles[p][IPX] < size.h_view(m).x1max &&
          particles[p][IPY] >= size.h_view(m).x2min &&
          particles[p][IPY] < size.h_view(m).x2max &&
          particles[p][IPZ] >= size.h_view(m).x3min &&
          particles[p][IPZ] < size.h_view(m).x3max) {
        gid = pmy_pack->gids + m;
        break;
      }
    }
    ihost(PGID,p) = gid;
    ihost(PTAG,p) = -1;
    for (int n=0; n<nrdata; ++n) {
      rhost(n,p) = particles[p][n];
    }
  }
  Kokkos::deep_copy(prtcl_rdata, rhost);
  Kokkos::deep_copy(prtcl_idata, ihost);
}

//----------------------------------------------------------------------------------------
// RefreshNextStarTag()
// Finds a globally unique starting tag for dynamically created star particles.

void Particles::RefreshNextStarTag() {
  int local_max = -1;
  if (nprtcl_thispack > 0) {
    auto ihost = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_idata);
    for (int p=0; p<nprtcl_thispack; ++p) {
      local_max = std::max(local_max, ihost(PTAG,p));
    }
  }
  int global_max = local_max;
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif
  next_star_tag = global_max + 1 + global_variable::my_rank;
}

//----------------------------------------------------------------------------------------
// UpdateParticleCounts()
// Updates mesh-level count caches after local particle creation/destruction.

void Particles::UpdateParticleCounts() {
  pmy_pack->pmesh->nprtcl_thisrank = nprtcl_thispack;
  pmy_pack->pmesh->nprtcl_total = 0;
  pmy_pack->pmesh->nprtcl_eachrank[global_variable::my_rank] = nprtcl_thispack;
#if MPI_PARALLEL_ENABLED
  MPI_Allgather(&nprtcl_thispack, 1, MPI_INT, pmy_pack->pmesh->nprtcl_eachrank, 1,
                MPI_INT, MPI_COMM_WORLD);
#endif
  for (int n=0; n<global_variable::nranks; ++n) {
    pmy_pack->pmesh->nprtcl_total += pmy_pack->pmesh->nprtcl_eachrank[n];
  }
}

} // namespace particles
