//========================================================================================
// AthenaK astrophysical fluid dynamics & numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file particles_gravity.cpp
//! \brief Star-particle external and self-gravity support.

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "bvals/bvals.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "parameter_input.hpp"
#include "pgen/pgen.hpp"
#include "particles.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {

struct GravityTreeNode {
  Real cx = 0.0, cy = 0.0, cz = 0.0;
  Real half = 0.0;
  Real mass = 0.0;
  Real cmx = 0.0, cmy = 0.0, cmz = 0.0;
  std::vector<int> particles;
  std::array<int, 8> child = {-1, -1, -1, -1, -1, -1, -1, -1};
};

inline Real SafeInvR3(Real r2, Real eps) {
  Real softened_r2 = r2 + eps*eps;
  if (softened_r2 <= 0.0) {return 0.0;}
  return 1.0/(softened_r2*std::sqrt(softened_r2));
}

inline void MinimumImage(Real &dx, Real length) {
  if (length <= 0.0) {return;}
  if (dx > 0.5*length) {
    dx -= length;
  } else if (dx < -0.5*length) {
    dx += length;
  }
}

inline bool IsLeaf(const GravityTreeNode &node) {
  for (int c=0; c<8; ++c) {
    if (node.child[c] >= 0) {return false;}
  }
  return true;
}

int NewTreeNode(std::vector<GravityTreeNode> &nodes, Real cx, Real cy, Real cz,
                Real half) {
  GravityTreeNode node;
  node.cx = cx;
  node.cy = cy;
  node.cz = cz;
  node.half = half;
  nodes.push_back(node);
  return static_cast<int>(nodes.size()) - 1;
}

int ChildIndex(const GravityTreeNode &node, Real x, Real y, Real z) {
  int idx = 0;
  if (x >= node.cx) {idx += 1;}
  if (y >= node.cy) {idx += 2;}
  if (z >= node.cz) {idx += 4;}
  return idx;
}

int EnsureChild(std::vector<GravityTreeNode> &nodes, int node_index, int child_index) {
  int existing = nodes[node_index].child[child_index];
  if (existing >= 0) {return existing;}

  Real offset = 0.5*nodes[node_index].half;
  Real cx = nodes[node_index].cx + ((child_index & 1) ? offset : -offset);
  Real cy = nodes[node_index].cy + ((child_index & 2) ? offset : -offset);
  Real cz = nodes[node_index].cz + ((child_index & 4) ? offset : -offset);
  int created = NewTreeNode(nodes, cx, cy, cz, offset);
  nodes[node_index].child[child_index] = created;
  return created;
}

void InsertParticle(std::vector<GravityTreeNode> &nodes, int node_index, int particle,
                    const std::vector<Real> &x, const std::vector<Real> &y,
                    const std::vector<Real> &z, int depth) {
  constexpr int max_depth = 48;
  if (!IsLeaf(nodes[node_index])) {
    int child_index =
        ChildIndex(nodes[node_index], x[particle], y[particle], z[particle]);
    InsertParticle(nodes, EnsureChild(nodes, node_index, child_index), particle,
                   x, y, z, depth + 1);
    return;
  }

  if (nodes[node_index].particles.empty() ||
      depth >= max_depth || nodes[node_index].half <= 0.0) {
    nodes[node_index].particles.push_back(particle);
    return;
  }

  std::vector<int> old_particles = std::move(nodes[node_index].particles);
  nodes[node_index].particles.clear();
  for (int old_particle : old_particles) {
    int old_child =
        ChildIndex(nodes[node_index], x[old_particle], y[old_particle], z[old_particle]);
    InsertParticle(nodes, EnsureChild(nodes, node_index, old_child), old_particle,
                   x, y, z, depth + 1);
  }
  int child_index = ChildIndex(nodes[node_index], x[particle], y[particle], z[particle]);
  InsertParticle(nodes, EnsureChild(nodes, node_index, child_index), particle, x, y, z,
                 depth + 1);
}

void BuildNodeMass(std::vector<GravityTreeNode> &nodes, int node_index,
                   const std::vector<Real> &x, const std::vector<Real> &y,
                   const std::vector<Real> &z, const std::vector<Real> &mass) {
  GravityTreeNode &node = nodes[node_index];
  if (IsLeaf(node)) {
    node.mass = 0.0;
    node.cmx = 0.0;
    node.cmy = 0.0;
    node.cmz = 0.0;
    for (int p : node.particles) {
      node.mass += mass[p];
      node.cmx += mass[p]*x[p];
      node.cmy += mass[p]*y[p];
      node.cmz += mass[p]*z[p];
    }
    if (node.mass > 0.0) {
      node.cmx /= node.mass;
      node.cmy /= node.mass;
      node.cmz /= node.mass;
    }
    return;
  }

  node.mass = 0.0;
  node.cmx = 0.0;
  node.cmy = 0.0;
  node.cmz = 0.0;
  for (int c=0; c<8; ++c) {
    int child = node.child[c];
    if (child < 0) {continue;}
    BuildNodeMass(nodes, child, x, y, z, mass);
    node.mass += nodes[child].mass;
    node.cmx += nodes[child].mass*nodes[child].cmx;
    node.cmy += nodes[child].mass*nodes[child].cmy;
    node.cmz += nodes[child].mass*nodes[child].cmz;
  }
  if (node.mass > 0.0) {
    node.cmx /= node.mass;
    node.cmy /= node.mass;
    node.cmz /= node.mass;
  }
}

bool PointInsideNode(const GravityTreeNode &node, Real x, Real y, Real z) {
  return (std::abs(x - node.cx) <= node.half &&
          std::abs(y - node.cy) <= node.half &&
          std::abs(z - node.cz) <= node.half);
}

void AddTreeAcceleration(const std::vector<GravityTreeNode> &nodes, int node_index,
                         int local_tag, Real x, Real y, Real z,
                         const std::vector<Real> &global_x,
                         const std::vector<Real> &global_y,
                         const std::vector<Real> &global_z,
                         const std::vector<Real> &global_mass,
                         const std::vector<int> &global_tag,
                         Real grav_constant, Real softening, Real theta,
                         Real &ax, Real &ay, Real &az) {
  const GravityTreeNode &node = nodes[node_index];
  if (node.mass <= 0.0) {return;}

  if (IsLeaf(node)) {
    for (int q : node.particles) {
      if (global_tag[q] == local_tag) {continue;}
      Real dx = global_x[q] - x;
      Real dy = global_y[q] - y;
      Real dz = global_z[q] - z;
      Real coeff = grav_constant*global_mass[q]*
          SafeInvR3(dx*dx + dy*dy + dz*dz, softening);
      ax += coeff*dx;
      ay += coeff*dy;
      az += coeff*dz;
    }
    return;
  } else if (PointInsideNode(node, x, y, z)) {
    for (int c=0; c<8; ++c) {
      if (node.child[c] >= 0) {
        AddTreeAcceleration(nodes, node.child[c], local_tag, x, y, z, global_x,
                            global_y, global_z, global_mass, global_tag, grav_constant,
                            softening, theta, ax, ay, az);
      }
    }
    return;
  }

  Real dx = node.cmx - x;
  Real dy = node.cmy - y;
  Real dz = node.cmz - z;
  Real r2 = dx*dx + dy*dy + dz*dz;
  Real dist = std::sqrt(r2 + softening*softening);
  if (!IsLeaf(node) && dist > 0.0 && (2.0*node.half/dist) >= theta) {
    for (int c=0; c<8; ++c) {
      if (node.child[c] >= 0) {
        AddTreeAcceleration(nodes, node.child[c], local_tag, x, y, z, global_x,
                            global_y, global_z, global_mass, global_tag, grav_constant,
                            softening, theta, ax, ay, az);
      }
    }
    return;
  }

  Real coeff = grav_constant*node.mass*SafeInvR3(r2, softening);
  ax += coeff*dx;
  ay += coeff*dy;
  az += coeff*dz;
}

} // namespace

namespace particles {

//----------------------------------------------------------------------------------------
//! \fn void Particles::ParseStarGravity()
//! \brief Parse star-particle gravity inputs and validate combinations early.

void Particles::ParseStarGravity(ParameterInput *pin) {
  bool block_exists = pin->DoesBlockExist("star_gravity");
  if (!block_exists) {
    if (pusher == ParticlesPusher::gravity) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "pusher = gravity requires a <star_gravity> block."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    return;
  }

  star_gravity_enabled = pin->GetOrAddBoolean("star_gravity", "enabled", false);
  if (!star_gravity_enabled) {
    if (pusher == ParticlesPusher::gravity) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "pusher = gravity requires <star_gravity>/enabled=true."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    return;
  }
  if (pusher != ParticlesPusher::gravity) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<star_gravity>/enabled=true requires pusher=gravity."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (particle_type != ParticleType::star) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<star_gravity> can only be used with particle_type=star."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  star_gravity_constant = pin->GetReal("star_gravity", "grav_constant");
  star_gravity_self_enabled =
      pin->GetOrAddBoolean("star_gravity", "self_gravity", true);
  star_gravity_softening =
      pin->GetOrAddReal("star_gravity", "softening_length", 0.0);
  star_gravity_timestep_eta =
      pin->GetOrAddReal("star_gravity", "timestep_eta", 0.02);
  star_gravity_max_timestep =
      pin->GetOrAddReal("star_gravity", "max_timestep", -1.0);
  star_gravity_tree_theta =
      pin->GetOrAddReal("star_gravity", "tree_theta", 0.7);
  star_gravity_exact_diagnostics =
      pin->GetOrAddBoolean("star_gravity", "exact_diagnostics", true);

  std::string force = pin->GetOrAddString("star_gravity", "force_method", "direct");
  if (force.compare("direct") == 0) {
    star_gravity_force_method = StarGravityForceMethod::direct;
  } else if (force.compare("tree") == 0) {
    star_gravity_force_method = StarGravityForceMethod::tree;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<star_gravity>/force_method='" << force
              << "' not recognized." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::string integrator = pin->GetOrAddString("star_gravity", "integrator", "kdk");
  if (integrator.compare("kdk") == 0) {
    star_gravity_integrator = StarGravityIntegrator::kdk;
  } else if (integrator.compare("rk4") == 0) {
    star_gravity_integrator = StarGravityIntegrator::rk4;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<star_gravity>/integrator='" << integrator
              << "' not recognized." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::string external =
      pin->GetOrAddString("star_gravity", "external_acceleration", "none");
  if (external.compare("none") == 0) {
    star_gravity_external_mode = StarGravityExternalMode::none;
  } else if (external.compare("constant") == 0) {
    star_gravity_external_mode = StarGravityExternalMode::constant;
  } else if (external.compare("point_mass") == 0 ||
             external.compare("plummer_point_mass") == 0) {
    star_gravity_external_mode = StarGravityExternalMode::point_mass;
  } else if (external.compare("user") == 0) {
    star_gravity_external_mode = StarGravityExternalMode::user;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<star_gravity>/external_acceleration='" << external
              << "' not recognized." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::string timestep =
      pin->GetOrAddString("star_gravity", "timestep_mode", "fixed");
  if (timestep.compare("fixed") == 0) {
    star_gravity_timestep_mode = StarGravityTimestepMode::fixed;
  } else if (timestep.compare("acceleration") == 0) {
    star_gravity_timestep_mode = StarGravityTimestepMode::acceleration;
  } else if (timestep.compare("pair_orbit") == 0) {
    star_gravity_timestep_mode = StarGravityTimestepMode::pair_orbit;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<star_gravity>/timestep_mode='" << timestep
              << "' not recognized." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::string periodic =
      pin->GetOrAddString("star_gravity", "periodic_mode", "domain");
  if (periodic.compare("domain") == 0) {
    star_gravity_periodic_mode = StarGravityPeriodicMode::domain;
  } else if (periodic.compare("minimum_image") == 0) {
    star_gravity_periodic_mode = StarGravityPeriodicMode::minimum_image;
  } else if (periodic.compare("error") == 0) {
    star_gravity_periodic_mode = StarGravityPeriodicMode::error;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<star_gravity>/periodic_mode='" << periodic
              << "' not recognized." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  star_gravity_external_ax =
      pin->GetOrAddReal("star_gravity", "external_accel_x", 0.0);
  star_gravity_external_ay =
      pin->GetOrAddReal("star_gravity", "external_accel_y", 0.0);
  star_gravity_external_az =
      pin->GetOrAddReal("star_gravity", "external_accel_z", 0.0);
  star_gravity_external_mass =
      pin->GetOrAddReal("star_gravity", "external_mass", 0.0);
  star_gravity_external_x =
      pin->GetOrAddReal("star_gravity", "external_x", 0.0);
  star_gravity_external_y =
      pin->GetOrAddReal("star_gravity", "external_y", 0.0);
  star_gravity_external_z =
      pin->GetOrAddReal("star_gravity", "external_z", 0.0);
  star_gravity_external_softening =
      pin->GetOrAddReal("star_gravity", "external_softening_length",
                        star_gravity_softening);

  if (!std::isfinite(star_gravity_constant) ||
      !std::isfinite(star_gravity_softening) ||
      !std::isfinite(star_gravity_timestep_eta) ||
      !std::isfinite(star_gravity_max_timestep) ||
      !std::isfinite(star_gravity_tree_theta) ||
      !std::isfinite(star_gravity_external_ax) ||
      !std::isfinite(star_gravity_external_ay) ||
      !std::isfinite(star_gravity_external_az) ||
      !std::isfinite(star_gravity_external_mass) ||
      !std::isfinite(star_gravity_external_x) ||
      !std::isfinite(star_gravity_external_y) ||
      !std::isfinite(star_gravity_external_z) ||
      !std::isfinite(star_gravity_external_softening) ||
      star_gravity_constant <= 0.0 || star_gravity_softening < 0.0 ||
      star_gravity_timestep_eta <= 0.0 || star_gravity_tree_theta <= 0.0 ||
      star_gravity_external_softening < 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Invalid <star_gravity> parameter.  G, eta, and "
              << "tree_theta must be positive; softening lengths must be non-negative."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (star_gravity_external_mode == StarGravityExternalMode::point_mass &&
      star_gravity_external_mass <= 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "point_mass external acceleration requires "
              << "positive external_mass." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (!star_gravity_self_enabled &&
      star_gravity_external_mode == StarGravityExternalMode::none) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<star_gravity> has no force source: self_gravity=false "
              << "and external_acceleration=none." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (star_gravity_force_method == StarGravityForceMethod::tree &&
      star_gravity_periodic_mode == StarGravityPeriodicMode::minimum_image) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "force_method=tree does not support "
              << "periodic_mode=minimum_image." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (star_gravity_force_method == StarGravityForceMethod::tree &&
      !star_gravity_exact_diagnostics &&
      star_gravity_timestep_mode == StarGravityTimestepMode::pair_orbit) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "timestep_mode=pair_orbit requires "
              << "exact_diagnostics=true." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if (star_gravity_periodic_mode == StarGravityPeriodicMode::error) {
    auto *pm = pmy_pack->pmesh;
    bool periodic_x1 =
        (pm->mesh_bcs[BoundaryFace::inner_x1] == BoundaryFlag::periodic ||
         pm->mesh_bcs[BoundaryFace::outer_x1] == BoundaryFlag::periodic ||
         pm->mesh_bcs[BoundaryFace::inner_x1] == BoundaryFlag::shear_periodic ||
         pm->mesh_bcs[BoundaryFace::outer_x1] == BoundaryFlag::shear_periodic);
    bool periodic_x2 =
        (pm->mesh_bcs[BoundaryFace::inner_x2] == BoundaryFlag::periodic ||
         pm->mesh_bcs[BoundaryFace::outer_x2] == BoundaryFlag::periodic);
    bool periodic_x3 =
        (pm->mesh_bcs[BoundaryFace::inner_x3] == BoundaryFlag::periodic ||
         pm->mesh_bcs[BoundaryFace::outer_x3] == BoundaryFlag::periodic);
    if (periodic_x1 || periodic_x2 || periodic_x3) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "periodic_mode=error forbids periodic mesh boundaries."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void Particles::ComputeStarGravityAcceleration()
//! \brief Compute local accelerations from configured external and star-star gravity.

void Particles::ComputeStarGravityAcceleration(const std::vector<Real> &x,
                                               const std::vector<Real> &y,
                                               const std::vector<Real> &z,
                                               const std::vector<Real> &mass,
                                               const std::vector<int> &tag,
                                               Real eval_time,
                                               std::vector<Real> &ax,
                                               std::vector<Real> &ay,
                                               std::vector<Real> &az,
                                               Real *min_pair,
                                               Real *potential_energy) {
  int local_n = static_cast<int>(x.size());
  ax.assign(local_n, 0.0);
  ay.assign(local_n, 0.0);
  az.assign(local_n, 0.0);
  Real local_min_pair = std::numeric_limits<Real>::max();
  Real local_epot = 0.0;

  if (star_gravity_external_mode == StarGravityExternalMode::constant) {
    for (int p=0; p<local_n; ++p) {
      ax[p] += star_gravity_external_ax;
      ay[p] += star_gravity_external_ay;
      az[p] += star_gravity_external_az;
    }
  } else if (star_gravity_external_mode == StarGravityExternalMode::point_mass) {
    for (int p=0; p<local_n; ++p) {
      Real dx = star_gravity_external_x - x[p];
      Real dy = star_gravity_external_y - y[p];
      Real dz = star_gravity_external_z - z[p];
      Real r2 = dx*dx + dy*dy + dz*dz;
      Real coeff = star_gravity_constant*star_gravity_external_mass*
                   SafeInvR3(r2, star_gravity_external_softening);
      ax[p] += coeff*dx;
      ay[p] += coeff*dy;
      az[p] += coeff*dz;
    }
  }

  std::vector<Real> global_x = x;
  std::vector<Real> global_y = y;
  std::vector<Real> global_z = z;
  std::vector<Real> global_m = mass;
  std::vector<int> global_tag = tag;

#if MPI_PARALLEL_ENABLED
  std::vector<int> counts(global_variable::nranks, 0);
  MPI_Allgather(&local_n, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<int> real_counts(global_variable::nranks, 0);
  std::vector<int> real_displs(global_variable::nranks, 0);
  std::vector<int> int_displs(global_variable::nranks, 0);
  int global_n = 0;
  for (int n=0; n<global_variable::nranks; ++n) {
    int_displs[n] = global_n;
    real_displs[n] = 4*global_n;
    real_counts[n] = 4*counts[n];
    global_n += counts[n];
  }

  std::vector<Real> local_pack(4*local_n, 0.0);
  for (int p=0; p<local_n; ++p) {
    local_pack[4*p] = x[p];
    local_pack[4*p + 1] = y[p];
    local_pack[4*p + 2] = z[p];
    local_pack[4*p + 3] = mass[p];
  }

  std::vector<Real> global_pack(4*global_n, 0.0);
  global_tag.assign(global_n, 0);
  MPI_Allgatherv(local_pack.data(), 4*local_n, MPI_ATHENA_REAL,
                 global_pack.data(), real_counts.data(), real_displs.data(),
                 MPI_ATHENA_REAL, MPI_COMM_WORLD);
  MPI_Allgatherv(const_cast<int*>(tag.data()), local_n, MPI_INT, global_tag.data(),
                 counts.data(), int_displs.data(), MPI_INT, MPI_COMM_WORLD);

  global_x.assign(global_n, 0.0);
  global_y.assign(global_n, 0.0);
  global_z.assign(global_n, 0.0);
  global_m.assign(global_n, 0.0);
  for (int p=0; p<global_n; ++p) {
    global_x[p] = global_pack[4*p];
    global_y[p] = global_pack[4*p + 1];
    global_z[p] = global_pack[4*p + 2];
    global_m[p] = global_pack[4*p + 3];
  }
#endif

  auto &mesh_size = pmy_pack->pmesh->mesh_size;
  Real lx = mesh_size.x1max - mesh_size.x1min;
  Real ly = mesh_size.x2max - mesh_size.x2min;
  Real lz = mesh_size.x3max - mesh_size.x3min;

  if (star_gravity_self_enabled &&
      star_gravity_force_method == StarGravityForceMethod::tree &&
      !global_x.empty()) {
    Real xmin = global_x[0], xmax = global_x[0];
    Real ymin = global_y[0], ymax = global_y[0];
    Real zmin = global_z[0], zmax = global_z[0];
    for (int q=1; q<static_cast<int>(global_x.size()); ++q) {
      xmin = std::min(xmin, global_x[q]);
      xmax = std::max(xmax, global_x[q]);
      ymin = std::min(ymin, global_y[q]);
      ymax = std::max(ymax, global_y[q]);
      zmin = std::min(zmin, global_z[q]);
      zmax = std::max(zmax, global_z[q]);
    }
    Real half = 0.5*std::max({xmax - xmin, ymax - ymin, zmax - zmin});
    half = std::max(half, static_cast<Real>(1.0e-12));
    std::vector<GravityTreeNode> nodes;
    nodes.reserve(2*global_x.size() + 1);
    NewTreeNode(nodes, 0.5*(xmin + xmax), 0.5*(ymin + ymax), 0.5*(zmin + zmax),
                half*(1.0 + 1.0e-12));
    for (int q=0; q<static_cast<int>(global_x.size()); ++q) {
      InsertParticle(nodes, 0, q, global_x, global_y, global_z, 0);
    }
    BuildNodeMass(nodes, 0, global_x, global_y, global_z, global_m);
    for (int p=0; p<local_n; ++p) {
      AddTreeAcceleration(nodes, 0, tag[p], x[p], y[p], z[p], global_x,
                          global_y, global_z, global_m, global_tag,
                          star_gravity_constant, star_gravity_softening,
                          star_gravity_tree_theta, ax[p], ay[p], az[p]);
    }
  }

  if (star_gravity_self_enabled &&
      (star_gravity_force_method == StarGravityForceMethod::direct ||
       min_pair != nullptr || potential_energy != nullptr)) {
    for (int p=0; p<local_n; ++p) {
      for (int q=0; q<static_cast<int>(global_x.size()); ++q) {
        if (tag[p] == global_tag[q]) {continue;}
        Real dx = global_x[q] - x[p];
        Real dy = global_y[q] - y[p];
        Real dz = global_z[q] - z[p];
        if (star_gravity_periodic_mode == StarGravityPeriodicMode::minimum_image) {
          MinimumImage(dx, lx);
          MinimumImage(dy, ly);
          MinimumImage(dz, lz);
        }
        Real r2 = dx*dx + dy*dy + dz*dz;
        if (r2 > 0.0) {
          local_min_pair = std::min(local_min_pair, std::sqrt(r2));
        }
        Real inv_r3 = SafeInvR3(r2, star_gravity_softening);
        if (star_gravity_force_method == StarGravityForceMethod::direct) {
          Real coeff = star_gravity_constant*global_m[q]*inv_r3;
          ax[p] += coeff*dx;
          ay[p] += coeff*dy;
          az[p] += coeff*dz;
        }

        if (tag[p] < global_tag[q]) {
          local_epot -= star_gravity_constant*mass[p]*global_m[q]/
              std::sqrt(r2 + star_gravity_softening*star_gravity_softening);
        }
      }
    }
  }

  if (star_gravity_external_mode == StarGravityExternalMode::user) {
    if (pmy_pack->pmesh->pgen == nullptr ||
        pmy_pack->pmesh->pgen->user_star_particle_accel_func == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "external_acceleration=user requires the pgen to set "
                << "user_star_particle_accel_func." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    pmy_pack->pmesh->pgen->user_star_particle_accel_func(
        pmy_pack->pmesh, this, eval_time, x.data(), y.data(), z.data(), ax.data(),
        ay.data(), az.data(), local_n);
  }

  if (min_pair != nullptr) {*min_pair = local_min_pair;}
  if (potential_energy != nullptr) {*potential_energy = local_epot;}
}

//----------------------------------------------------------------------------------------
//! \fn void Particles::GravityKDKDrift()
//! \brief Kick velocities by half a timestep and drift positions by one timestep.

void Particles::GravityKDKDrift() {
  if (!star_gravity_enabled) {return;}

  Real dt = pmy_pack->pmesh->dt;
  auto rhost = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_rdata);
  auto ihost = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_idata);
  std::vector<Real> x(nprtcl_thispack), y(nprtcl_thispack), z(nprtcl_thispack);
  std::vector<Real> mass(nprtcl_thispack);
  std::vector<int> tag(nprtcl_thispack);
  for (int p=0; p<nprtcl_thispack; ++p) {
    x[p] = rhost(IPX,p);
    y[p] = rhost(IPY,p);
    z[p] = rhost(IPZ,p);
    mass[p] = rhost(IPMASS,p);
    tag[p] = ihost(PTAG,p);
  }

  std::vector<Real> ax, ay, az;
  ComputeStarGravityAcceleration(x, y, z, mass, tag, pmy_pack->pmesh->time,
                                 ax, ay, az);
  for (int p=0; p<nprtcl_thispack; ++p) {
    rhost(IPVX,p) += 0.5*dt*ax[p];
    rhost(IPVY,p) += 0.5*dt*ay[p];
    rhost(IPVZ,p) += 0.5*dt*az[p];
    rhost(IPX,p) += dt*rhost(IPVX,p);
    rhost(IPY,p) += dt*rhost(IPVY,p);
    rhost(IPZ,p) += dt*rhost(IPVZ,p);
  }
  Kokkos::deep_copy(prtcl_rdata, rhost);
}

//----------------------------------------------------------------------------------------
//! \fn void Particles::GravityKDKFinalKick()
//! \brief Complete the second velocity kick after particle migration.

void Particles::GravityKDKFinalKick() {
  if (!star_gravity_enabled) {return;}

  Real dt = pmy_pack->pmesh->dt;
  auto rhost = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_rdata);
  auto ihost = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_idata);
  std::vector<Real> x(nprtcl_thispack), y(nprtcl_thispack), z(nprtcl_thispack);
  std::vector<Real> mass(nprtcl_thispack);
  std::vector<int> tag(nprtcl_thispack);
  for (int p=0; p<nprtcl_thispack; ++p) {
    x[p] = rhost(IPX,p);
    y[p] = rhost(IPY,p);
    z[p] = rhost(IPZ,p);
    mass[p] = rhost(IPMASS,p);
    tag[p] = ihost(PTAG,p);
  }

  std::vector<Real> ax, ay, az;
  ComputeStarGravityAcceleration(x, y, z, mass, tag, pmy_pack->pmesh->time + dt,
                                 ax, ay, az);
  for (int p=0; p<nprtcl_thispack; ++p) {
    rhost(IPVX,p) += 0.5*dt*ax[p];
    rhost(IPVY,p) += 0.5*dt*ay[p];
    rhost(IPVZ,p) += 0.5*dt*az[p];
  }
  Kokkos::deep_copy(prtcl_rdata, rhost);
}

//----------------------------------------------------------------------------------------
//! \fn void Particles::GravityRK4Step()
//! \brief Advance local particles using a fourth-order Runge-Kutta particle update.

void Particles::GravityRK4Step() {
  if (!star_gravity_enabled) {return;}

  Real dt = pmy_pack->pmesh->dt;
  auto rhost = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_rdata);
  auto ihost = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_idata);
  int n = nprtcl_thispack;
  std::vector<Real> x0(n), y0(n), z0(n), vx0(n), vy0(n), vz0(n), mass(n);
  std::vector<int> tag(n);
  for (int p=0; p<n; ++p) {
    x0[p] = rhost(IPX,p);
    y0[p] = rhost(IPY,p);
    z0[p] = rhost(IPZ,p);
    vx0[p] = rhost(IPVX,p);
    vy0[p] = rhost(IPVY,p);
    vz0[p] = rhost(IPVZ,p);
    mass[p] = rhost(IPMASS,p);
    tag[p] = ihost(PTAG,p);
  }

  std::vector<Real> ax1, ay1, az1, ax2, ay2, az2, ax3, ay3, az3, ax4, ay4, az4;
  std::vector<Real> xt(n), yt(n), zt(n), vxt(n), vyt(n), vzt(n);

  Real time = pmy_pack->pmesh->time;
  ComputeStarGravityAcceleration(x0, y0, z0, mass, tag, time, ax1, ay1, az1);
  for (int p=0; p<n; ++p) {
    xt[p] = x0[p] + 0.5*dt*vx0[p];
    yt[p] = y0[p] + 0.5*dt*vy0[p];
    zt[p] = z0[p] + 0.5*dt*vz0[p];
    vxt[p] = vx0[p] + 0.5*dt*ax1[p];
    vyt[p] = vy0[p] + 0.5*dt*ay1[p];
    vzt[p] = vz0[p] + 0.5*dt*az1[p];
  }

  ComputeStarGravityAcceleration(xt, yt, zt, mass, tag, time + 0.5*dt,
                                 ax2, ay2, az2);
  for (int p=0; p<n; ++p) {
    xt[p] = x0[p] + 0.5*dt*vxt[p];
    yt[p] = y0[p] + 0.5*dt*vyt[p];
    zt[p] = z0[p] + 0.5*dt*vzt[p];
    vxt[p] = vx0[p] + 0.5*dt*ax2[p];
    vyt[p] = vy0[p] + 0.5*dt*ay2[p];
    vzt[p] = vz0[p] + 0.5*dt*az2[p];
  }

  ComputeStarGravityAcceleration(xt, yt, zt, mass, tag, time + 0.5*dt,
                                 ax3, ay3, az3);
  for (int p=0; p<n; ++p) {
    xt[p] = x0[p] + dt*vxt[p];
    yt[p] = y0[p] + dt*vyt[p];
    zt[p] = z0[p] + dt*vzt[p];
    vxt[p] = vx0[p] + dt*ax3[p];
    vyt[p] = vy0[p] + dt*ay3[p];
    vzt[p] = vz0[p] + dt*az3[p];
  }

  ComputeStarGravityAcceleration(xt, yt, zt, mass, tag, time + dt,
                                 ax4, ay4, az4);
  for (int p=0; p<n; ++p) {
    rhost(IPX,p) = x0[p] + dt*(vx0[p] + 2.0*(vx0[p] + 0.5*dt*ax1[p]) +
                               2.0*(vx0[p] + 0.5*dt*ax2[p]) +
                               (vx0[p] + dt*ax3[p]))/6.0;
    rhost(IPY,p) = y0[p] + dt*(vy0[p] + 2.0*(vy0[p] + 0.5*dt*ay1[p]) +
                               2.0*(vy0[p] + 0.5*dt*ay2[p]) +
                               (vy0[p] + dt*ay3[p]))/6.0;
    rhost(IPZ,p) = z0[p] + dt*(vz0[p] + 2.0*(vz0[p] + 0.5*dt*az1[p]) +
                               2.0*(vz0[p] + 0.5*dt*az2[p]) +
                               (vz0[p] + dt*az3[p]))/6.0;
    rhost(IPVX,p) = vx0[p] + dt*(ax1[p] + 2.0*ax2[p] + 2.0*ax3[p] + ax4[p])/6.0;
    rhost(IPVY,p) = vy0[p] + dt*(ay1[p] + 2.0*ay2[p] + 2.0*ay3[p] + ay4[p])/6.0;
    rhost(IPVZ,p) = vz0[p] + dt*(az1[p] + 2.0*az2[p] + 2.0*az3[p] + az4[p])/6.0;
  }
  Kokkos::deep_copy(prtcl_rdata, rhost);
}

//----------------------------------------------------------------------------------------
//! \fn void Particles::UpdateGravityDiagnosticsAndTimestep()
//! \brief Update global star-particle gravity diagnostics and dtnew.

void Particles::UpdateGravityDiagnosticsAndTimestep() {
  if (!star_gravity_enabled) {return;}

  auto rhost = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_rdata);
  auto ihost = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_idata);
  int n = nprtcl_thispack;
  std::vector<Real> x(n), y(n), z(n), mass(n);
  std::vector<int> tag(n);
  for (int p=0; p<n; ++p) {
    x[p] = rhost(IPX,p);
    y[p] = rhost(IPY,p);
    z[p] = rhost(IPZ,p);
    mass[p] = rhost(IPMASS,p);
    tag[p] = ihost(PTAG,p);
  }

  std::vector<Real> ax, ay, az;
  Real local_rmin = std::numeric_limits<Real>::max();
  Real local_epot = 0.0;
  Real *rmin = nullptr;
  Real *epot = nullptr;
  if (star_gravity_force_method == StarGravityForceMethod::direct ||
      star_gravity_exact_diagnostics) {
    rmin = &local_rmin;
    epot = &local_epot;
  }
  ComputeStarGravityAcceleration(x, y, z, mass, tag,
                                 pmy_pack->pmesh->time + pmy_pack->pmesh->dt,
                                 ax, ay, az, rmin, epot);

  Real local_mass = 0.0;
  Real local_ekin = 0.0;
  Real local_px = 0.0, local_py = 0.0, local_pz = 0.0;
  Real local_lx = 0.0, local_ly = 0.0, local_lz = 0.0;
  Real local_mx = 0.0, local_my = 0.0, local_mz = 0.0;
  Real local_amax = 0.0;
  for (int p=0; p<n; ++p) {
    Real m = mass[p];
    Real vx = rhost(IPVX,p);
    Real vy = rhost(IPVY,p);
    Real vz = rhost(IPVZ,p);
    local_mass += m;
    local_ekin += 0.5*m*(vx*vx + vy*vy + vz*vz);
    local_px += m*vx;
    local_py += m*vy;
    local_pz += m*vz;
    local_lx += m*(y[p]*vz - z[p]*vy);
    local_ly += m*(z[p]*vx - x[p]*vz);
    local_lz += m*(x[p]*vy - y[p]*vx);
    local_mx += m*x[p];
    local_my += m*y[p];
    local_mz += m*z[p];
    local_amax = std::max(local_amax, std::sqrt(ax[p]*ax[p] + ay[p]*ay[p] +
                                                az[p]*az[p]));
  }

  std::array<Real, 12> global_sums = {
      local_mass, local_ekin, local_epot, local_px, local_py, local_pz,
      local_lx, local_ly, local_lz, local_mx, local_my, local_mz};
  Real global_amax = local_amax;
  Real global_rmin = local_rmin;
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, global_sums.data(), global_sums.size(), MPI_ATHENA_REAL,
                MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &global_amax, 1, MPI_ATHENA_REAL, MPI_MAX,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &global_rmin, 1, MPI_ATHENA_REAL, MPI_MIN,
                MPI_COMM_WORLD);
#endif

  Real global_mass = global_sums[0];
  Real global_ekin = global_sums[1];
  Real global_epot = global_sums[2];
  Real global_px = global_sums[3], global_py = global_sums[4];
  Real global_pz = global_sums[5];
  Real global_lx = global_sums[6], global_ly = global_sums[7];
  Real global_lz = global_sums[8];
  Real global_mx = global_sums[9], global_my = global_sums[10];
  Real global_mz = global_sums[11];
  star_gravity_diag_ekin = global_ekin;
  star_gravity_diag_epot = global_epot;
  star_gravity_diag_px = global_px;
  star_gravity_diag_py = global_py;
  star_gravity_diag_pz = global_pz;
  star_gravity_diag_lx = global_lx;
  star_gravity_diag_ly = global_ly;
  star_gravity_diag_lz = global_lz;
  if (global_mass > 0.0) {
    star_gravity_diag_com_x = global_mx/global_mass;
    star_gravity_diag_com_y = global_my/global_mass;
    star_gravity_diag_com_z = global_mz/global_mass;
  } else {
    star_gravity_diag_com_x = 0.0;
    star_gravity_diag_com_y = 0.0;
    star_gravity_diag_com_z = 0.0;
  }
  star_gravity_diag_amax = global_amax;
  star_gravity_diag_rmin =
      (global_rmin == std::numeric_limits<Real>::max()) ? 0.0 : global_rmin;

  Real gravity_dt = star_gravity_base_timestep;
  if (star_gravity_timestep_mode == StarGravityTimestepMode::fixed) {
    gravity_dt = star_gravity_base_timestep;
  } else if (global_amax > 0.0) {
    Real length = star_gravity_softening;
    if (length <= 0.0 &&
        star_gravity_external_mode == StarGravityExternalMode::point_mass) {
      length = star_gravity_external_softening;
    }
    if (length <= 0.0 && star_gravity_diag_rmin > 0.0) {
      length = star_gravity_diag_rmin;
    }
    if (length > 0.0) {
      gravity_dt = star_gravity_timestep_eta*std::sqrt(length/global_amax);
    }
  }
  if (star_gravity_timestep_mode == StarGravityTimestepMode::pair_orbit &&
      star_gravity_diag_rmin > 0.0 && global_mass > 0.0) {
    Real orbit_dt = star_gravity_timestep_eta*
        std::sqrt((star_gravity_diag_rmin*star_gravity_diag_rmin*
                   star_gravity_diag_rmin)/(star_gravity_constant*global_mass));
    gravity_dt = std::min(gravity_dt, orbit_dt);
  }
  if (star_gravity_max_timestep > 0.0) {
    gravity_dt = std::min(gravity_dt, star_gravity_max_timestep);
  }
  if (gravity_dt > 0.0 && std::isfinite(gravity_dt)) {
    dtnew = gravity_dt;
  }
  star_gravity_diag_dt = dtnew;
}

} // namespace particles
