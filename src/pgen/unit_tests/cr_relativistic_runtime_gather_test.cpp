//========================================================================================
// AthenaK astrophysical fluid dynamics & numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cr_relativistic_runtime_gather_test.cpp
//! \brief Diagnostic-only runtime sampling of live MHD fields at CR particle positions.

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "parameter_input.hpp"
#include "particles/particles.hpp"
#include "particles/trilinear_gather.hpp"
#include "pgen/pgen.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {

enum VelocityProfile {
  kUniformVelocity = 0,
  kAffineVelocity = 1,
  kPeriodicVelocity = 2
};

constexpr Real kTwoPi = 6.283185307179586476925286766559;

struct RuntimeGatherRow {
  int cycle;
  int rank;
  int gid;
  int species;
  int tag;
  Real time;
  Real x1, x2, x3;
  Real u1, u2, u3;
  Real b1, b2, b3;
  Real cE1, cE2, cE3;
  Real cE_dot_b;
};

[[noreturn]] void Fail(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl
            << message << std::endl;
#if MPI_PARALLEL_ENABLED
  int initialized = 0;
  MPI_Initialized(&initialized);
  if (initialized != 0) {MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);}
#endif
  std::exit(EXIT_FAILURE);
}

KOKKOS_INLINE_FUNCTION
particles::ParticleVector3 PrescribedVelocity(
    int profile, Real x1, Real x2, Real x3, Real u1_0, Real u2_0, Real u3_0,
    Real u1_x1, Real u1_x2, Real u1_x3, Real u2_x1, Real u2_x2, Real u2_x3,
    Real u3_x1, Real u3_x2, Real u3_x3, Real uamp, Real uwave) {
  particles::ParticleVector3 u{u1_0, u2_0, u3_0};
  if (profile == kAffineVelocity) {
    u.x += u1_x1*x1 + u1_x2*x2 + u1_x3*x3;
    u.y += u2_x1*x1 + u2_x2*x2 + u2_x3*x3;
    u.z += u3_x1*x1 + u3_x2*x2 + u3_x3*x3;
  } else if (profile == kPeriodicVelocity) {
    const Real k = kTwoPi*uwave;
    u.x += uamp*Kokkos::sin(k*x1)*Kokkos::cos(k*x2);
    u.y += uamp*Kokkos::sin(k*x2)*Kokkos::cos(k*x3);
    u.z += uamp*Kokkos::sin(k*x3)*Kokkos::cos(k*x1);
  }
  return u;
}

int ReadVelocityProfile(ParameterInput *pin) {
  std::string profile = pin->GetOrAddString("problem", "fluid_velocity_profile",
                                            "uniform");
  if (profile.compare("uniform") == 0) {return kUniformVelocity;}
  if (profile.compare("affine") == 0) {return kAffineVelocity;}
  if (profile.compare("periodic") == 0) {return kPeriodicVelocity;}
  Fail("<problem>/fluid_velocity_profile must be 'uniform', 'affine', or 'periodic'");
}

void PlaceBoundaryProbeParticles(ParameterInput *pin, Mesh *pm) {
  std::string layout = pin->GetOrAddString("problem", "runtime_probe_layout", "existing");
  if (layout.compare("existing") == 0) {return;}
  const bool boundary_tags = (layout.compare("boundary_tags") == 0);
  const bool refined_region_tags = (layout.compare("refined_region_tags") == 0);
  if (!boundary_tags && !refined_region_tags) {
    Fail("<problem>/runtime_probe_layout must be 'existing', 'boundary_tags', or "
         "'refined_region_tags'");
  }

  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &pr = pmbp->ppart->prtcl_rdata;
  auto &pi = pmbp->ppart->prtcl_idata;
  const RegionSize mesh_size = pm->mesh_size;
  const int npart = pmbp->ppart->nprtcl_thispack;
  const Real eps1 = 1.0e-10*(mesh_size.x1max - mesh_size.x1min);
  const Real eps2 = 1.0e-10*(mesh_size.x2max - mesh_size.x2min);
  const Real eps3 = 1.0e-10*(mesh_size.x3max - mesh_size.x3min);
  const Real mid1 = 0.5*(mesh_size.x1min + mesh_size.x1max);
  const Real mid2 = 0.5*(mesh_size.x2min + mesh_size.x2max);
  const Real mid3 = 0.5*(mesh_size.x3min + mesh_size.x3max);
  const Real anchor1 = mid1 + 0.17*(mesh_size.x1max - mesh_size.x1min);
  const Real anchor2 = mid2 - 0.13*(mesh_size.x2max - mesh_size.x2min);
  const Real anchor3 = mid3 + 0.11*(mesh_size.x3max - mesh_size.x3min);
  const Real refined_half_width =
      pin->GetOrAddReal("problem", "runtime_refined_half_width", 0.25);
  par_for("cr_rel_runtime_place_boundary_probes", DevExeSpace(), 0, npart - 1,
  KOKKOS_LAMBDA(const int p) {
    pr(IPX,p) = anchor1;
    pr(IPY,p) = anchor2;
    pr(IPZ,p) = anchor3;
    if (boundary_tags) {
      const int pattern = pi(PTAG,p) % 14;
      if (pattern == 0 || pattern == 12) {pr(IPX,p) = mesh_size.x1min + eps1;}
      if (pattern == 1 || pattern == 13) {pr(IPX,p) = mesh_size.x1max - eps1;}
      if (pattern == 2) {pr(IPX,p) = mid1 - eps1;}
      if (pattern == 3) {pr(IPX,p) = mid1 + eps1;}
      if (pattern == 4 || pattern == 12) {pr(IPY,p) = mesh_size.x2min + eps2;}
      if (pattern == 5 || pattern == 13) {pr(IPY,p) = mesh_size.x2max - eps2;}
      if (pattern == 6) {pr(IPY,p) = mid2 - eps2;}
      if (pattern == 7) {pr(IPY,p) = mid2 + eps2;}
      if (pattern == 8 || pattern == 12) {pr(IPZ,p) = mesh_size.x3min + eps3;}
      if (pattern == 9 || pattern == 13) {pr(IPZ,p) = mesh_size.x3max - eps3;}
      if (pattern == 10) {pr(IPZ,p) = mid3 - eps3;}
      if (pattern == 11) {pr(IPZ,p) = mid3 + eps3;}
    } else {
      const int pattern = pi(PTAG,p) % 12;
      if (pattern == 0) {pr(IPX,p) = mid1 - refined_half_width - eps1;}
      if (pattern == 1) {pr(IPX,p) = mid1 - refined_half_width + eps1;}
      if (pattern == 2) {pr(IPX,p) = mid1 + refined_half_width - eps1;}
      if (pattern == 3) {pr(IPX,p) = mid1 + refined_half_width + eps1;}
      if (pattern == 4) {pr(IPY,p) = mid2 - refined_half_width - eps2;}
      if (pattern == 5) {pr(IPY,p) = mid2 - refined_half_width + eps2;}
      if (pattern == 6) {pr(IPY,p) = mid2 + refined_half_width - eps2;}
      if (pattern == 7) {pr(IPY,p) = mid2 + refined_half_width + eps2;}
      if (pattern == 8) {pr(IPZ,p) = mid3 - refined_half_width - eps3;}
      if (pattern == 9) {pr(IPZ,p) = mid3 - refined_half_width + eps3;}
      if (pattern == 10) {pr(IPZ,p) = mid3 + refined_half_width - eps3;}
      if (pattern == 11) {pr(IPZ,p) = mid3 + refined_half_width + eps3;}
    }
  });
  pmbp->ppart->RemapAfterAMR();
}

void SampleRuntimeGather(ParameterInput *pin, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp == nullptr || pmbp->pmhd == nullptr || pmbp->ppart == nullptr) {
    Fail("CR runtime gather diagnostic requires MHD and particle modules");
  }

  auto &pr = pmbp->ppart->prtcl_rdata;
  auto &pi = pmbp->ppart->prtcl_idata;
  auto &w0 = pmbp->pmhd->w0;
  auto &bcc0 = pmbp->pmhd->bcc0;
  auto &mbsize = pmbp->pmb->mb_size;
  auto &indcs = pm->mb_indcs;
  const int npart = pmbp->ppart->nprtcl_thispack;
  const int gids = pmbp->gids;
  const int nmb = pmbp->nmb_thispack;
  const bool multi_d = pm->multi_d;
  const bool three_d = pm->three_d;

  DvceArray2D<Real> sample("cr_rel_runtime_sample", 10, npart);
  DvceArray1D<int> valid("cr_rel_runtime_valid", npart);
  Kokkos::parallel_for("cr_rel_runtime_gather",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, npart),
  KOKKOS_LAMBDA(const int p) {
    const int m = pi(PGID,p) - gids;
    valid(p) = (m >= 0 && m < nmb) ? 1 : 0;
    if (valid(p) == 1) {
      auto fields = particles::GatherNewtonianCellCenteredFields(
          w0, bcc0, mbsize, m, pr(IPX,p), pr(IPY,p), pr(IPZ,p), indcs.is,
          indcs.js, indcs.ks, multi_d, three_d);
      auto cE = particles::ConstructIdealCE(fields.u, fields.b);
      sample(0,p) = fields.u.x;
      sample(1,p) = fields.u.y;
      sample(2,p) = fields.u.z;
      sample(3,p) = fields.b.x;
      sample(4,p) = fields.b.y;
      sample(5,p) = fields.b.z;
      sample(6,p) = cE.x;
      sample(7,p) = cE.y;
      sample(8,p) = cE.z;
      sample(9,p) = cE.x*fields.b.x + cE.y*fields.b.y + cE.z*fields.b.z;
    }
  });

  auto pr_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pr);
  auto pi_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pi);
  auto sample_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), sample);
  auto valid_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), valid);

  std::vector<RuntimeGatherRow> local_rows;
  local_rows.reserve(npart);
  for (int p=0; p<npart; ++p) {
    if (valid_h(p) != 1) {
      Fail("CR runtime gather encountered a particle with nonlocal PGID");
    }
    RuntimeGatherRow row;
    row.cycle = pm->ncycle;
    row.rank = global_variable::my_rank;
    row.gid = pi_h(PGID,p);
    row.species = pi_h(PSP,p);
    row.tag = pi_h(PTAG,p);
    row.time = pm->time;
    row.x1 = pr_h(IPX,p);
    row.x2 = pr_h(IPY,p);
    row.x3 = pr_h(IPZ,p);
    row.u1 = sample_h(0,p);
    row.u2 = sample_h(1,p);
    row.u3 = sample_h(2,p);
    row.b1 = sample_h(3,p);
    row.b2 = sample_h(4,p);
    row.b3 = sample_h(5,p);
    row.cE1 = sample_h(6,p);
    row.cE2 = sample_h(7,p);
    row.cE3 = sample_h(8,p);
    row.cE_dot_b = sample_h(9,p);
    local_rows.push_back(row);
  }

  std::vector<RuntimeGatherRow> global_rows;
  std::vector<int> particle_counts;
  std::vector<int> participants;
#if MPI_PARALLEL_ENABLED
  const int nranks = global_variable::nranks;
  const int local_rank = global_variable::my_rank;
  const int local_bytes = static_cast<int>(local_rows.size()*sizeof(RuntimeGatherRow));
  std::vector<int> byte_counts(nranks, 0);
  std::vector<int> byte_displs(nranks, 0);
  participants.resize(nranks);
  MPI_Gather(&local_rank, 1, MPI_INT, participants.data(), 1, MPI_INT, 0,
             MPI_COMM_WORLD);
  MPI_Gather(&local_bytes, 1, MPI_INT, byte_counts.data(), 1, MPI_INT, 0,
             MPI_COMM_WORLD);
  if (global_variable::my_rank == 0) {
    particle_counts.resize(nranks);
    int total_bytes = 0;
    for (int rank=0; rank<nranks; ++rank) {
      byte_displs[rank] = total_bytes;
      total_bytes += byte_counts[rank];
      particle_counts[rank] =
          byte_counts[rank]/static_cast<int>(sizeof(RuntimeGatherRow));
    }
    if (total_bytes % static_cast<int>(sizeof(RuntimeGatherRow)) != 0) {
      Fail("CR runtime gather received a malformed MPI byte count");
    }
    global_rows.resize(total_bytes/sizeof(RuntimeGatherRow));
  }
  MPI_Gatherv(local_rows.data(), local_bytes, MPI_BYTE, global_rows.data(),
              byte_counts.data(), byte_displs.data(), MPI_BYTE, 0,
              MPI_COMM_WORLD);
#else
  global_rows = local_rows;
  particle_counts = {static_cast<int>(local_rows.size())};
  participants = {0};
#endif

  if (global_variable::my_rank != 0) {return;}
  std::sort(global_rows.begin(), global_rows.end(),
  [](const RuntimeGatherRow &a, const RuntimeGatherRow &b) {
    return std::tie(a.species, a.tag, a.rank, a.gid) <
           std::tie(b.species, b.tag, b.rank, b.gid);
  });

  std::string filename = pin->GetOrAddString(
      "problem", "runtime_gather_tsv", "cr_relativistic_runtime_gather.tsv");
  std::ofstream stream(filename, std::ios::out | std::ios::trunc);
  if (!stream.is_open()) {
    Fail("could not open CR runtime gather diagnostic file '" + filename + "'");
  }
  stream << "# diagnostic-only\tpost-finalize\ttrilinear_cell_centered"
         << "\tcE=-u_cross_B\n";
  stream << "# rank_counts";
  for (std::size_t rank=0; rank<particle_counts.size(); ++rank) {
    stream << '\t' << rank << '=' << particle_counts[rank];
  }
  stream << '\n';
  stream << "# rank_participants\tworld_size=" << participants.size();
  for (const int rank : participants) {
    stream << '\t' << rank;
  }
  stream << '\n';
  stream << "# cycle\ttime\trank\tgid\tspecies\ttag\tx\ty\tz"
         << "\tu1\tu2\tu3\tb1\tb2\tb3\tcE1\tcE2\tcE3\tcE_dot_b\n";
  stream << std::setprecision(std::numeric_limits<Real>::max_digits10);
  for (const auto &row : global_rows) {
    stream << row.cycle << '\t' << row.time << '\t' << row.rank << '\t'
           << row.gid << '\t' << row.species << '\t' << row.tag << '\t'
           << row.x1 << '\t' << row.x2 << '\t' << row.x3 << '\t'
           << row.u1 << '\t' << row.u2 << '\t' << row.u3 << '\t'
           << row.b1 << '\t' << row.b2 << '\t' << row.b3 << '\t'
           << row.cE1 << '\t' << row.cE2 << '\t' << row.cE3 << '\t'
           << row.cE_dot_b << '\n';
  }
}

} // namespace

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  PartRandom(pin, restart);
  pgen_final_func = &SampleRuntimeGather;
  if (restart) {return;}
  PlaceBoundaryProbeParticles(pin, pmy_mesh_);

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp == nullptr || pmbp->pmhd == nullptr) {
    Fail("CR runtime gather diagnostic requires an <mhd> block");
  }

  const int profile = ReadVelocityProfile(pin);
  const Real density = pin->GetOrAddReal("problem", "fluid_density", 1.0);
  const Real pressure = pin->GetOrAddReal("problem", "fluid_pressure", 0.6);
  const Real u1_0 = pin->GetOrAddReal("problem", "u1_0", 0.0);
  const Real u2_0 = pin->GetOrAddReal("problem", "u2_0", 0.0);
  const Real u3_0 = pin->GetOrAddReal("problem", "u3_0", 0.0);
  const Real u1_x1 = pin->GetOrAddReal("problem", "u1_x1", 0.0);
  const Real u1_x2 = pin->GetOrAddReal("problem", "u1_x2", 0.0);
  const Real u1_x3 = pin->GetOrAddReal("problem", "u1_x3", 0.0);
  const Real u2_x1 = pin->GetOrAddReal("problem", "u2_x1", 0.0);
  const Real u2_x2 = pin->GetOrAddReal("problem", "u2_x2", 0.0);
  const Real u2_x3 = pin->GetOrAddReal("problem", "u2_x3", 0.0);
  const Real u3_x1 = pin->GetOrAddReal("problem", "u3_x1", 0.0);
  const Real u3_x2 = pin->GetOrAddReal("problem", "u3_x2", 0.0);
  const Real u3_x3 = pin->GetOrAddReal("problem", "u3_x3", 0.0);
  const Real uamp = pin->GetOrAddReal("problem", "fluid_velocity_amplitude", 0.0);
  const Real uwave = pin->GetOrAddReal("problem", "fluid_velocity_wave_number", 1.0);

  auto &w0 = pmbp->pmhd->w0;
  auto &bcc0 = pmbp->pmhd->bcc0;
  auto &u0 = pmbp->pmhd->u0;
  auto &size = pmbp->pmb->mb_size;
  auto &indcs = pmy_mesh_->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  par_for("cr_rel_runtime_prescribed_velocity", DevExeSpace(),
          0, (pmbp->nmb_thispack - 1), ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real x1 = CellCenterX(i - is, nx1, size.d_view(m).x1min,
                                size.d_view(m).x1max);
    const Real x2 = CellCenterX(j - js, nx2, size.d_view(m).x2min,
                                size.d_view(m).x2max);
    const Real x3 = CellCenterX(k - ks, nx3, size.d_view(m).x3min,
                                size.d_view(m).x3max);
    auto velocity = PrescribedVelocity(
        profile, x1, x2, x3, u1_0, u2_0, u3_0, u1_x1, u1_x2, u1_x3,
        u2_x1, u2_x2, u2_x3, u3_x1, u3_x2, u3_x3, uamp, uwave);
    w0(m,IDN,k,j,i) = density;
    w0(m,IVX,k,j,i) = velocity.x;
    w0(m,IVY,k,j,i) = velocity.y;
    w0(m,IVZ,k,j,i) = velocity.z;
    w0(m,IEN,k,j,i) = pressure;
  });
  pmbp->pmhd->peos->PrimToCons(w0, bcc0, u0, is, ie, js, je, ks, ke);
}
