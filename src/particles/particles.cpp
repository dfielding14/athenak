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
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "hydro/hydro.hpp"
#include "bvals/bvals.hpp"
#include "particle_restart.hpp"
#include "particles.hpp"

namespace particles {
namespace {

void FatalParticleRestartRead(const char *file, int line, const std::string &msg) {
  std::cout << "### FATAL ERROR in " << file << " at line " << line << std::endl
            << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

bool IsOuterMeshBoundary(Real block_max, Real mesh_max) {
  Real scale = std::max(static_cast<Real>(1.0), std::abs(mesh_max));
  Real tol = 100.0*std::numeric_limits<Real>::epsilon()*scale;
  return std::abs(block_max - mesh_max) <= tol;
}

} // namespace
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
    } else if (ppush.compare("gravity") == 0 ||
               ppush.compare("star_gravity") == 0) {
      pusher = ParticlesPusher::gravity;
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
        } else if (init_mode.compare("restart") == 0 ||
                   init_mode.compare("particle_restart") == 0) {
          star_init_from_restart = true;
        } else if (init_mode.compare("none") != 0) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<particles>/star_init = '" << init_mode
                    << "' not recognized. Valid choices are none, file, or restart."
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (pin->DoesParameterExist("particles", "star_prtcl_rst_flag") &&
            pin->GetInteger("particles", "star_prtcl_rst_flag") != 0) {
          star_init_from_restart = true;
          star_init_from_file = false;
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
        star_gravity_conserve_accretion_momentum =
            pin->GetOrAddBoolean("particles", "star_accretion_conserve_momentum",
                                 true);

        if (star_formation_interval < 1 || star_formation_max_per_cycle < -1 ||
            star_accretion_radius_cells < 0 ||
            !std::isfinite(star_formation_density_threshold) ||
            !std::isfinite(star_formation_particle_mass) ||
            !std::isfinite(star_formation_density_floor) ||
            !std::isfinite(star_accretion_rate) ||
            !std::isfinite(star_accretion_density_floor) ||
            !std::isfinite(star_accretion_max_fraction) ||
            star_formation_density_floor < 0.0 ||
            star_accretion_density_floor < 0.0 ||
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
        ParseStarGravity(pin);
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
    if (star_init_from_restart) {
      LoadStarsFromRestart(pin);
    } else if (star_init_from_file) {
      LoadStarsFromFile(pin);
    }
    Real default_dt = std::min(pmy_pack->pmb->mb_size.h_view(0).dx1,
                               pmy_pack->pmb->mb_size.h_view(0).dx2);
    default_dt = std::min(default_dt, pmy_pack->pmb->mb_size.h_view(0).dx3);
    dtnew = pin->GetOrAddReal("particles", "dt", default_dt);
    if (!std::isfinite(dtnew) || dtnew <= 0.0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<particles>/dt must be finite and positive."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    star_gravity_base_timestep = dtnew;
    if (star_gravity_enabled && star_gravity_max_timestep > 0.0) {
      dtnew = std::min(dtnew, star_gravity_max_timestep);
    }
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
      star_accretion_rate == 0.0) {
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

  int local_n = nprtcl_thispack;
  int global_n = local_n;
  std::vector<Real> local_stars(6*local_n, 0.0);
  std::vector<int> local_tags(local_n, -1);
  for (int p=0; p<local_n; ++p) {
    int m = ihost(PGID,p) - pmy_pack->gids;
    if (m < 0 || m >= pmy_pack->nmb_thispack) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Star particle has invalid MeshBlock ownership before "
                << "accretion." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    int ic = static_cast<int>((rhost(IPX,p) - size.h_view(m).x1min)/
                              size.h_view(m).dx1) + is;
    int jc = static_cast<int>((rhost(IPY,p) - size.h_view(m).x2min)/
                              size.h_view(m).dx2) + js;
    int kc = static_cast<int>((rhost(IPZ,p) - size.h_view(m).x3min)/
                              size.h_view(m).dx3) + ks;
    ic = std::max(is, std::min(ie, ic));
    jc = std::max(js, std::min(je, jc));
    kc = std::max(ks, std::min(ke, kc));
    local_stars[6*p] = CellCenterX(ic - is, indcs.nx1, size.h_view(m).x1min,
                                   size.h_view(m).x1max);
    local_stars[6*p + 1] = CellCenterX(jc - js, indcs.nx2, size.h_view(m).x2min,
                                       size.h_view(m).x2max);
    local_stars[6*p + 2] = CellCenterX(kc - ks, indcs.nx3, size.h_view(m).x3min,
                                       size.h_view(m).x3max);
    local_stars[6*p + 3] = size.h_view(m).dx1;
    local_stars[6*p + 4] = size.h_view(m).dx2;
    local_stars[6*p + 5] = size.h_view(m).dx3;
    local_tags[p] = ihost(PTAG,p);
  }

  std::vector<Real> global_stars = local_stars;
  std::vector<int> global_tags = local_tags;
#if MPI_PARALLEL_ENABLED
  std::vector<int> counts(global_variable::nranks, 0);
  MPI_Allgather(&local_n, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
  std::vector<int> real_counts(global_variable::nranks, 0);
  std::vector<int> real_displs(global_variable::nranks, 0);
  std::vector<int> int_displs(global_variable::nranks, 0);
  global_n = 0;
  for (int n=0; n<global_variable::nranks; ++n) {
    int_displs[n] = global_n;
    real_displs[n] = 6*global_n;
    real_counts[n] = 6*counts[n];
    global_n += counts[n];
  }
  global_stars.assign(6*global_n, 0.0);
  MPI_Allgatherv(local_stars.data(), 6*local_n, MPI_ATHENA_REAL,
                 global_stars.data(), real_counts.data(), real_displs.data(),
                 MPI_ATHENA_REAL, MPI_COMM_WORLD);
  global_tags.assign(global_n, -1);
  MPI_Allgatherv(local_tags.data(), local_n, MPI_INT, global_tags.data(),
                 counts.data(), int_displs.data(), MPI_INT, MPI_COMM_WORLD);
#endif
  if (global_n == 0) {return TaskStatus::complete;}

  std::vector<int> order(global_n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(),
            [&](int a, int b) {return global_tags[a] < global_tags[b];});
  std::vector<Real> sorted_stars(6*global_n, 0.0);
  std::map<int, int> tag_to_index;
  for (int q=0; q<global_n; ++q) {
    int source = order[q];
    if (!tag_to_index.emplace(global_tags[source], q).second) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Star particle tags must be unique during accretion."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    for (int n=0; n<6; ++n) {
      sorted_stars[6*q + n] = global_stars[6*source + n];
    }
  }
  global_stars = std::move(sorted_stars);

  auto cell_range = [](Real x, Real radius, Real xmin, Real dx, int start, int end) {
    Real scale = std::max({static_cast<Real>(1.0), std::abs(x), std::abs(xmin)});
    Real tol = 100.0*std::numeric_limits<Real>::epsilon()*scale;
    int lo = static_cast<int>(std::ceil((x - radius - xmin - tol)/dx - 0.5)) + start;
    int hi = static_cast<int>(std::floor((x + radius - xmin + tol)/dx - 0.5)) + start;
    return std::array<int, 2>{std::max(start, lo), std::min(end, hi)};
  };

  std::vector<Real> global_accretion(4*global_n, 0.0);
  Real local_removed_mass = 0.0;
  for (int q=0; q<global_n; ++q) {
    Real x = global_stars[6*q];
    Real y = global_stars[6*q + 1];
    Real z = global_stars[6*q + 2];
    Real rx = star_accretion_radius_cells*global_stars[6*q + 3];
    Real ry = star_accretion_radius_cells*global_stars[6*q + 4];
    Real rz = star_accretion_radius_cells*global_stars[6*q + 5];
    for (int m=0; m<pmy_pack->nmb_thispack; ++m) {
      auto irange = cell_range(x, rx, size.h_view(m).x1min, size.h_view(m).dx1,
                               is, ie);
      auto jrange = cell_range(y, ry, size.h_view(m).x2min, size.h_view(m).dx2,
                               js, je);
      auto krange = cell_range(z, rz, size.h_view(m).x3min, size.h_view(m).dx3,
                               ks, ke);
      for (int k=krange[0]; k<=krange[1]; ++k) {
        for (int j=jrange[0]; j<=jrange[1]; ++j) {
          for (int i=irange[0]; i<=irange[1]; ++i) {
            Real rho = uhost(m,IDN,k,j,i);
            Real available_rho = std::max(static_cast<Real>(0.0),
                                          rho - star_accretion_density_floor);
            Real delta_rho = frac_limit*available_rho;
            if (delta_rho <= 0.0 || rho <= 0.0) {continue;}
            Real remove_frac = std::min(delta_rho/rho, static_cast<Real>(1.0));
            Real vol = size.h_view(m).dx1*size.h_view(m).dx2*size.h_view(m).dx3;
            Real removed_mass = delta_rho*vol;
            global_accretion[4*q] += removed_mass;
            if (star_gravity_conserve_accretion_momentum) {
              global_accretion[4*q + 1] += removed_mass*whost(m,IVX,k,j,i);
              global_accretion[4*q + 2] += removed_mass*whost(m,IVY,k,j,i);
              global_accretion[4*q + 3] += removed_mass*whost(m,IVZ,k,j,i);
            }
            for (int n=0; n<nhydro_total; ++n) {
              uhost(m,n,k,j,i) *= (1.0 - remove_frac);
            }
            whost(m,IDN,k,j,i) = uhost(m,IDN,k,j,i);
            if (pmy_pack->phydro->nhydro > IEN) {
              whost(m,IEN,k,j,i) *= (1.0 - remove_frac);
            }
            local_removed_mass += removed_mass;
          }
        }
      }
    }
  }
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, global_accretion.data(), 4*global_n, MPI_ATHENA_REAL,
                MPI_SUM, MPI_COMM_WORLD);
#endif

  Real accreted_mass = 0.0;
  for (int p=0; p<local_n; ++p) {
    int q = tag_to_index.at(ihost(PTAG,p));
    Real particle_gain = global_accretion[4*q];
    Real old_mass = rhost(IPMASS,p);
    Real old_mom_x = old_mass*rhost(IPVX,p);
    Real old_mom_y = old_mass*rhost(IPVY,p);
    Real old_mom_z = old_mass*rhost(IPVZ,p);
    rhost(IPMASS,p) += particle_gain;
    if (star_gravity_conserve_accretion_momentum && particle_gain > 0.0 &&
        rhost(IPMASS,p) > 0.0) {
      rhost(IPVX,p) = (old_mom_x + global_accretion[4*q + 1])/rhost(IPMASS,p);
      rhost(IPVY,p) = (old_mom_y + global_accretion[4*q + 2])/rhost(IPMASS,p);
      rhost(IPVZ,p) = (old_mom_z + global_accretion[4*q + 3])/rhost(IPMASS,p);
    }
    accreted_mass += particle_gain;
  }

  if (accreted_mass > 0.0) {
    Kokkos::deep_copy(prtcl_rdata, rhost);
    star_mass_accreted_total += accreted_mass;
  }
  if (local_removed_mass > 0.0) {
    Kokkos::deep_copy(pmy_pack->phydro->u0, uhost);
    Kokkos::deep_copy(pmy_pack->phydro->w0, whost);
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
// CreateParticleTags()
// Assigns tags to particles (unique integer).  Note that tracked particles are always
// those with tag numbers less than ntrack.

void Particles::CreateParticleTags(ParameterInput *pin) {
  if (particle_type == ParticleType::star && star_init_from_restart) {
    RefreshNextStarTag();
    return;
  }

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
// FindLocalMeshBlockForPosition()
// Returns the local MeshBlock index that owns a position, or -1 if this rank owns none.

int Particles::FindLocalMeshBlockForPosition(Real x, Real y, Real z) const {
  auto &size = pmy_pack->pmb->mb_size;
  auto &mesh_size = pmy_pack->pmesh->mesh_size;
  for (int m=0; m<pmy_pack->nmb_thispack; ++m) {
    bool in_x = (x >= size.h_view(m).x1min &&
                 (x < size.h_view(m).x1max ||
                  (IsOuterMeshBoundary(size.h_view(m).x1max, mesh_size.x1max) &&
                   x <= size.h_view(m).x1max)));
    bool in_y = (y >= size.h_view(m).x2min &&
                 (y < size.h_view(m).x2max ||
                  (IsOuterMeshBoundary(size.h_view(m).x2max, mesh_size.x2max) &&
                   y <= size.h_view(m).x2max)));
    bool in_z = (z >= size.h_view(m).x3min &&
                 (z < size.h_view(m).x3max ||
                  (IsOuterMeshBoundary(size.h_view(m).x3max, mesh_size.x3max) &&
                   z <= size.h_view(m).x3max)));
    if (in_x && in_y && in_z) {return m;}
  }
  return -1;
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
  int line_number = 0;
  int expected_particles = 0;
  while (std::getline(infile, line)) {
    line_number++;
    std::size_t first = line.find_first_not_of(" \t");
    if (first == std::string::npos || line[first] == '#') {continue;}
    std::istringstream iss(line);
    Real x, y, z, vx, vy, vz, col7, col8;
    if (!(iss >> x >> y >> z >> vx >> vy >> vz >> col7)) {
      std::ostringstream msg;
      msg << "Unable to parse star particle file '" << particle_file
          << "' at line " << line_number << ".";
      FatalParticleRestartRead(__FILE__, __LINE__, msg.str());
    }
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
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z) ||
        !std::isfinite(vx) || !std::isfinite(vy) || !std::isfinite(vz) ||
        !std::isfinite(mass) || !std::isfinite(tcreate) || mass <= 0.0) {
      std::ostringstream msg;
      msg << "Star particle file '" << particle_file << "' has invalid values at line "
          << line_number << ". Mass must be finite and positive.";
      FatalParticleRestartRead(__FILE__, __LINE__, msg.str());
    }
    expected_particles++;

    for (int m=0; m<nmb; ++m) {
      if (x >= size.h_view(m).x1min && x < size.h_view(m).x1max &&
          y >= size.h_view(m).x2min && y < size.h_view(m).x2max &&
          z >= size.h_view(m).x3min && z < size.h_view(m).x3max) {
        particles.push_back({x, vx, y, vy, z, vz, mass, tcreate});
        break;
      }
    }
  }

  int loaded_particles = static_cast<int>(particles.size());
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &loaded_particles, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (loaded_particles != expected_particles) {
    std::ostringstream msg;
    msg << "Star particle file '" << particle_file << "' contains "
        << expected_particles << " particles, but mesh ownership assigned "
        << loaded_particles << ". Check that every initial position is inside the mesh.";
    FatalParticleRestartRead(__FILE__, __LINE__, msg.str());
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
// LoadStarsFromRestart()
// Restores stars from a global sidecar restart and redistributes by position.

void Particles::LoadStarsFromRestart(ParameterInput *pin) {
  std::string restart_file;
  if (pin->DoesParameterExist("particles", "particle_restart_file")) {
    restart_file = pin->GetString("particles", "particle_restart_file");
  } else if (pin->DoesParameterExist("particles", "star_particle_file")) {
    restart_file = pin->GetString("particles", "star_particle_file");
  } else {
    FatalParticleRestartRead(__FILE__, __LINE__,
        "star_init=restart requires <particles>/particle_restart_file.");
  }

  std::ifstream infile(restart_file, std::ios::binary);
  if (!infile) {
    FatalParticleRestartRead(__FILE__, __LINE__,
        "Unable to open particle restart file '" + restart_file + "'.");
  }

  ParticleRestartHeader header;
  infile.read(reinterpret_cast<char*>(&header), sizeof(header));
  if (!infile) {
    FatalParticleRestartRead(__FILE__, __LINE__,
        "Unable to read particle restart header from '" + restart_file + "'.");
  }

  if (std::strncmp(header.magic, kParticleRestartMagic,
                   kParticleRestartMagicSize) != 0 ||
      header.version != kParticleRestartVersion) {
    FatalParticleRestartRead(__FILE__, __LINE__,
        "Particle restart file '" + restart_file + "' has an invalid magic "
        "or version.");
  }
  if (header.header_size != static_cast<int>(sizeof(header)) ||
      header.endian_marker != kParticleRestartEndianMarker) {
    FatalParticleRestartRead(__FILE__, __LINE__,
        "Particle restart file '" + restart_file + "' has an incompatible "
        "header layout or byte order.");
  }
  if (header.real_size != static_cast<int>(sizeof(Real)) ||
      header.int_size != static_cast<int>(sizeof(int))) {
    FatalParticleRestartRead(__FILE__, __LINE__,
        "Particle restart file '" + restart_file + "' was written with incompatible "
        "Real or int sizes.");
  }
  if (header.nrdata != nrdata || header.nidata != nidata) {
    FatalParticleRestartRead(__FILE__, __LINE__,
        "Particle restart file '" + restart_file + "' has particle array dimensions "
        "that do not match this executable/input.");
  }
  if (header.total_particles < 0 || header.nrdata <= 0 || header.nidata <= 0) {
    FatalParticleRestartRead(__FILE__, __LINE__,
        "Particle restart file '" + restart_file + "' has invalid particle counts.");
  }
  if (header.basename_length < 0 ||
      header.basename_length >= kParticleRestartBasenameSize ||
      header.basename[header.basename_length] != '\0' ||
      std::strlen(header.basename) !=
          static_cast<std::size_t>(header.basename_length)) {
    FatalParticleRestartRead(__FILE__, __LINE__,
        "Particle restart file '" + restart_file + "' has invalid basename metadata.");
  }
  if (!std::isfinite(header.time) ||
      !std::isfinite(header.star_mass_formed_total) ||
      !std::isfinite(header.star_mass_accreted_total) ||
      header.star_mass_formed_total < 0.0 ||
      header.star_mass_accreted_total < 0.0 ||
      (header.total_particles == 0 && header.max_tag != -1) ||
      (header.total_particles > 0 && header.max_tag < 0)) {
    FatalParticleRestartRead(__FILE__, __LINE__,
        "Particle restart file '" + restart_file + "' has invalid header metadata.");
  }

  std::uintmax_t expected_size = sizeof(header) +
      static_cast<std::uintmax_t>(header.total_particles)*header.nrdata*sizeof(Real) +
      static_cast<std::uintmax_t>(header.total_particles)*header.nidata*sizeof(int);
  infile.seekg(0, std::ios::end);
  std::streamoff actual_size = infile.tellg();
  if (actual_size < 0 ||
      static_cast<std::uintmax_t>(actual_size) != expected_size) {
    FatalParticleRestartRead(__FILE__, __LINE__,
        "Particle restart file '" + restart_file + "' has an unexpected file size.");
  }
  infile.seekg(sizeof(header), std::ios::beg);

  Real mesh_time = pmy_pack->pmesh->time;
  Real default_tol = 100.0*std::numeric_limits<Real>::epsilon()*
      std::max(static_cast<Real>(1.0),
               std::max(std::abs(mesh_time), std::abs(header.time)));
  Real time_tol = pin->GetOrAddReal("particles", "particle_restart_time_tolerance",
                                    default_tol);
  if (!std::isfinite(time_tol) || time_tol < 0.0) {
    FatalParticleRestartRead(__FILE__, __LINE__,
        "<particles>/particle_restart_time_tolerance must be finite and non-negative.");
  }
  if (header.ncycle != pmy_pack->pmesh->ncycle ||
      std::abs(header.time - mesh_time) > time_tol) {
    std::ostringstream msg;
    msg << "Particle restart file '" << restart_file << "' does not match the mesh "
        << "restart state. mesh(time=" << mesh_time
        << ", cycle=" << pmy_pack->pmesh->ncycle << "), particles(time="
        << header.time << ", cycle=" << header.ncycle << ").";
    FatalParticleRestartRead(__FILE__, __LINE__, msg.str());
  }

  int total_particles = header.total_particles;
  std::vector<Real> all_rdata(static_cast<std::size_t>(total_particles)*nrdata);
  std::vector<int> all_idata(static_cast<std::size_t>(total_particles)*nidata);
  if (total_particles > 0) {
    infile.read(reinterpret_cast<char*>(all_rdata.data()),
                all_rdata.size()*sizeof(Real));
    infile.read(reinterpret_cast<char*>(all_idata.data()),
                all_idata.size()*sizeof(int));
    if (!infile) {
      FatalParticleRestartRead(__FILE__, __LINE__,
          "Particle restart file '" + restart_file + "' ended before all particle data "
          "were read.");
    }
  }
  infile.close();

  std::set<int> restored_tags;
  for (int p=0; p<total_particles; ++p) {
    int tag = all_idata[static_cast<std::size_t>(p)*nidata + PTAG];
    if (tag < 0 || !restored_tags.insert(tag).second) {
      FatalParticleRestartRead(__FILE__, __LINE__,
          "Particle restart file '" + restart_file +
          "' contains invalid or duplicate particle tags.");
    }
  }

  std::vector<Real> local_rdata;
  std::vector<int> local_idata;
  if (total_particles > 0) {
    std::size_t reserve_particles =
        static_cast<std::size_t>(total_particles/global_variable::nranks + 1);
    local_rdata.reserve(reserve_particles*nrdata);
    local_idata.reserve(reserve_particles*nidata);
  }

  int local_max_tag = -1;
  for (int p=0; p<total_particles; ++p) {
    Real x = all_rdata[static_cast<std::size_t>(p)*nrdata + IPX];
    Real y = all_rdata[static_cast<std::size_t>(p)*nrdata + IPY];
    Real z = all_rdata[static_cast<std::size_t>(p)*nrdata + IPZ];
    int local_mb = FindLocalMeshBlockForPosition(x, y, z);
    if (local_mb < 0) {continue;}

    std::size_t rbase = static_cast<std::size_t>(p)*nrdata;
    std::size_t ibase = static_cast<std::size_t>(p)*nidata;
    for (int n=0; n<nrdata; ++n) {
      local_rdata.push_back(all_rdata[rbase + n]);
    }
    for (int n=0; n<nidata; ++n) {
      local_idata.push_back(all_idata[ibase + n]);
    }
    std::size_t local_base = local_idata.size() - nidata;
    local_idata[local_base + PGID] = pmy_pack->gids + local_mb;
    local_max_tag = std::max(local_max_tag, local_idata[local_base + PTAG]);
  }

  nprtcl_thispack = static_cast<int>(local_rdata.size()/nrdata);
  Kokkos::realloc(prtcl_rdata, nrdata, nprtcl_thispack);
  Kokkos::realloc(prtcl_idata, nidata, nprtcl_thispack);
  auto rhost = Kokkos::create_mirror_view(prtcl_rdata);
  auto ihost = Kokkos::create_mirror_view(prtcl_idata);
  for (int p=0; p<nprtcl_thispack; ++p) {
    for (int n=0; n<nrdata; ++n) {
      rhost(n,p) = local_rdata[static_cast<std::size_t>(p)*nrdata + n];
    }
    for (int n=0; n<nidata; ++n) {
      ihost(n,p) = local_idata[static_cast<std::size_t>(p)*nidata + n];
    }
  }
  Kokkos::deep_copy(prtcl_rdata, rhost);
  Kokkos::deep_copy(prtcl_idata, ihost);

  int restored_particles = nprtcl_thispack;
  int restored_max_tag = local_max_tag;
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &restored_particles, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &restored_max_tag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif
  if (restored_particles != total_particles) {
    std::ostringstream msg;
    msg << "Particle restart file '" << restart_file << "' restored "
        << restored_particles << " particles, but its header records "
        << total_particles << ".";
    FatalParticleRestartRead(__FILE__, __LINE__, msg.str());
  }
  if (restored_max_tag != header.max_tag) {
    std::ostringstream msg;
    msg << "Particle restart file '" << restart_file << "' restored max tag "
        << restored_max_tag << ", but its header records " << header.max_tag << ".";
    FatalParticleRestartRead(__FILE__, __LINE__, msg.str());
  }

  star_mass_formed_total =
      (global_variable::my_rank == 0) ? header.star_mass_formed_total : 0.0;
  star_mass_accreted_total =
      (global_variable::my_rank == 0) ? header.star_mass_accreted_total : 0.0;
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
