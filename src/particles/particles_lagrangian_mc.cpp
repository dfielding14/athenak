//========================================================================================
// AthenaK astrophysical fluid dynamics & numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file particles_lagrangian_mc.cpp
//! \brief generic Monte-Carlo mass-flux tracer particles with thermodynamic sampling

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "driver/driver.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "outputs/io_wrapper.hpp"
#include "parameter_input.hpp"
#include "particles.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {

KOKKOS_INLINE_FUNCTION
std::uint64_t SplitMix64(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

KOKKOS_INLINE_FUNCTION
Real HashReal(std::uint64_t x) {
  return static_cast<Real>((SplitMix64(x) >> 11) *
                           (1.0/9007199254740992.0));
}

bool StartsWith(const std::string &value, const std::string &prefix) {
  return value.compare(0, prefix.size(), prefix) == 0;
}

void FatalParticleInput(const std::string &msg) {
  std::cout << "### FATAL ERROR in particles_lagrangian_mc.cpp" << std::endl
            << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

struct CandidateCell {
  int m, k, j, i;
  Real cumulative_end;
};

struct ThermoSample {
  Real rho, pressure, temperature, entropy, eint, v1, v2, v3;
};

RegionSize MeshBlockSizeFromGID(Mesh *pm, int gid) {
  RegionSize mb_size = pm->mesh_size;
  LogicalLocation &lloc = pm->lloc_eachmb[gid];
  std::int32_t lev = lloc.level;

  std::int32_t nmbx1 = pm->nmb_rootx1 << (lev - pm->root_level);
  if (lloc.lx1 == 0) {
    mb_size.x1min = pm->mesh_size.x1min;
  } else {
    mb_size.x1min = LeftEdgeX(lloc.lx1, nmbx1, pm->mesh_size.x1min,
                              pm->mesh_size.x1max);
  }
  if (lloc.lx1 == nmbx1 - 1) {
    mb_size.x1max = pm->mesh_size.x1max;
  } else {
    mb_size.x1max = LeftEdgeX(lloc.lx1 + 1, nmbx1, pm->mesh_size.x1min,
                              pm->mesh_size.x1max);
  }

  if (pm->multi_d) {
    std::int32_t nmbx2 = pm->nmb_rootx2 << (lev - pm->root_level);
    if (lloc.lx2 == 0) {
      mb_size.x2min = pm->mesh_size.x2min;
    } else {
      mb_size.x2min = LeftEdgeX(lloc.lx2, nmbx2, pm->mesh_size.x2min,
                                pm->mesh_size.x2max);
    }
    if (lloc.lx2 == nmbx2 - 1) {
      mb_size.x2max = pm->mesh_size.x2max;
    } else {
      mb_size.x2max = LeftEdgeX(lloc.lx2 + 1, nmbx2, pm->mesh_size.x2min,
                                pm->mesh_size.x2max);
    }
  } else {
    mb_size.x2min = pm->mesh_size.x2min;
    mb_size.x2max = pm->mesh_size.x2max;
  }

  if (pm->three_d) {
    std::int32_t nmbx3 = pm->nmb_rootx3 << (lev - pm->root_level);
    if (lloc.lx3 == 0) {
      mb_size.x3min = pm->mesh_size.x3min;
    } else {
      mb_size.x3min = LeftEdgeX(lloc.lx3, nmbx3, pm->mesh_size.x3min,
                                pm->mesh_size.x3max);
    }
    if (lloc.lx3 == nmbx3 - 1) {
      mb_size.x3max = pm->mesh_size.x3max;
    } else {
      mb_size.x3max = LeftEdgeX(lloc.lx3 + 1, nmbx3, pm->mesh_size.x3min,
                                pm->mesh_size.x3max);
    }
  } else {
    mb_size.x3min = pm->mesh_size.x3min;
    mb_size.x3max = pm->mesh_size.x3max;
  }

  mb_size.dx1 = (mb_size.x1max - mb_size.x1min)/static_cast<Real>(pm->mb_indcs.nx1);
  mb_size.dx2 = (mb_size.x2max - mb_size.x2min)/static_cast<Real>(pm->mb_indcs.nx2);
  mb_size.dx3 = (mb_size.x3max - mb_size.x3min)/static_cast<Real>(pm->mb_indcs.nx3);
  return mb_size;
}

Real WrapOrClamp(Real x, Real xmin, Real xmax, bool periodic) {
  Real width = xmax - xmin;
  if (periodic && width > 0.0) {
    while (x < xmin) x += width;
    while (x >= xmax) x -= width;
    return x;
  }
  if (x < xmin) return xmin;
  if (x >= xmax) return std::nextafter(xmax, xmin);
  return x;
}

bool ContainsPosition(const RegionSize &size, Real x1, Real x2, Real x3) {
  return (x1 >= size.x1min && x1 < size.x1max &&
          x2 >= size.x2min && x2 < size.x2max &&
          x3 >= size.x3min && x3 < size.x3max);
}

int LocateMeshBlockGID(Mesh *pm, Real x1, Real x2, Real x3) {
  for (int gid=0; gid<pm->nmb_total; ++gid) {
    RegionSize size = MeshBlockSizeFromGID(pm, gid);
    if (ContainsPosition(size, x1, x2, x3)) return gid;
  }
  return -1;
}

void SnapToMeshBlockCellCenter(Mesh *pm, int gid, Real &x1, Real &x2, Real &x3) {
  RegionSize size = MeshBlockSizeFromGID(pm, gid);
  auto &indcs = pm->mb_indcs;
  int i = static_cast<int>((x1 - size.x1min)/size.dx1);
  int j = static_cast<int>((x2 - size.x2min)/size.dx2);
  int k = static_cast<int>((x3 - size.x3min)/size.dx3);
  i = std::max(0, std::min(indcs.nx1 - 1, i));
  j = std::max(0, std::min(indcs.nx2 - 1, j));
  k = std::max(0, std::min(indcs.nx3 - 1, k));
  x1 = CellCenterX(i, indcs.nx1, size.x1min, size.x1max);
  x2 = CellCenterX(j, indcs.nx2, size.x2min, size.x2max);
  x3 = CellCenterX(k, indcs.nx3, size.x3min, size.x3max);
}

ThermoSample SampleFromHost(const HostArray5D<Real> &w0, const EOS_Data &eos,
                            int m, int k, int j, int i) {
  ThermoSample out;
  out.rho = w0(m,IDN,k,j,i);
  out.eint = w0(m,IEN,k,j,i);
  out.pressure = eos.IdealGasPressure(out.eint);
  out.temperature = out.pressure/out.rho;
  out.entropy = std::log(out.pressure/std::pow(out.rho, eos.gamma));
  out.v1 = w0(m,IVX,k,j,i);
  out.v2 = w0(m,IVY,k,j,i);
  out.v3 = w0(m,IVZ,k,j,i);
  return out;
}

} // namespace

namespace particles {

//----------------------------------------------------------------------------------------
//! \fn void Particles::ParseTracerSeedSchedules

void Particles::ParseTracerSeedSchedules(ParameterInput *pin) {
  seed_schedules_.clear();
  int fallback_id = 1;

  for (auto it = pin->block.begin(); it != pin->block.end(); ++it) {
    const std::string &block = it->block_name;
    if (!StartsWith(block, "tracer_seed")) continue;

    TracerSeedSchedule sched;
    sched.block_name = block;
    sched.id = pin->GetOrAddInteger(block, "id", fallback_id);
    sched.start_time = pin->GetOrAddReal(block, "start_time", 0.0);
    sched.end_time = pin->GetOrAddReal(block, "end_time", sched.start_time);
    sched.cadence = pin->GetOrAddReal(block, "cadence", -1.0);
    sched.next_time = pin->GetOrAddReal(block, "next_time", sched.start_time);
    sched.event_index = pin->GetOrAddInteger(block, "event_index", 0);
    sched.complete = pin->GetOrAddBoolean(block, "complete", false);
    sched.count_per_event = pin->GetInteger(block, "count_per_event");
    sched.seed = pin->GetOrAddInteger(block, "seed", 0);

    std::string weight = pin->GetOrAddString(block, "weight", "mass");
    if (weight == "mass") {
      sched.weight = TracerSeedWeight::mass;
    } else if (weight == "volume") {
      sched.weight = TracerSeedWeight::volume;
    } else {
      FatalParticleInput(block + ": weight must be mass or volume");
    }

    Mesh *pm = pmy_pack->pmesh;
    sched.x1min = pin->GetOrAddReal(block, "x1min", pm->mesh_size.x1min);
    sched.x1max = pin->GetOrAddReal(block, "x1max", pm->mesh_size.x1max);
    sched.x2min = pin->GetOrAddReal(block, "x2min", pm->mesh_size.x2min);
    sched.x2max = pin->GetOrAddReal(block, "x2max", pm->mesh_size.x2max);
    sched.x3min = pin->GetOrAddReal(block, "x3min", pm->mesh_size.x3min);
    sched.x3max = pin->GetOrAddReal(block, "x3max", pm->mesh_size.x3max);
    sched.center1 = pin->GetOrAddReal(block, "center1", 0.5*(pm->mesh_size.x1min +
                                                              pm->mesh_size.x1max));
    sched.center2 = pin->GetOrAddReal(block, "center2", 0.5*(pm->mesh_size.x2min +
                                                              pm->mesh_size.x2max));
    sched.center3 = pin->GetOrAddReal(block, "center3", 0.5*(pm->mesh_size.x3min +
                                                              pm->mesh_size.x3max));
    sched.radius = pin->GetOrAddReal(block, "radius", 0.0);
    sched.slab_axis = pin->GetOrAddInteger(block, "slab_axis", 1);
    sched.slab_min = pin->GetOrAddReal(block, "slab_min", sched.x1min);
    sched.slab_max = pin->GetOrAddReal(block, "slab_max", sched.x1max);

    std::string region = pin->GetOrAddString(block, "region", "all");
    if (region == "all") {
      sched.region = TracerSeedRegion::all;
    } else if (region == "box") {
      sched.region = TracerSeedRegion::box;
    } else if (region == "sphere") {
      sched.region = TracerSeedRegion::sphere;
      if (sched.radius <= 0.0) FatalParticleInput(block + ": sphere radius must be > 0");
    } else if (region == "slab") {
      sched.region = TracerSeedRegion::slab;
      if (sched.slab_axis < 1 || sched.slab_axis > 3) {
        FatalParticleInput(block + ": slab_axis must be 1, 2, or 3");
      }
    } else {
      FatalParticleInput(block + ": region must be all, box, sphere, or slab");
    }

    if (pin->DoesParameterExist(block, "target")) {
      std::string target = pin->GetString(block, "target");
      if (target == "density") {
        sched.target = TracerSeedTarget::density;
      } else if (target == "temperature") {
        sched.target = TracerSeedTarget::temperature;
      } else if (target == "pressure") {
        sched.target = TracerSeedTarget::pressure;
      } else if (target == "entropy") {
        sched.target = TracerSeedTarget::entropy;
      } else if (StartsWith(target, "scalar")) {
        sched.target = TracerSeedTarget::scalar;
        sched.scalar_index = std::stoi(target.substr(6));
        if (sched.scalar_index < 0 || sched.scalar_index >= GetLagrangianMCScalarCount()) {
          FatalParticleInput(block + ": scalar target index is out of range");
        }
      } else {
        FatalParticleInput(block + ": target must be density, temperature, pressure, "
                           "entropy, or scalarN");
      }
      sched.has_target_min = pin->DoesParameterExist(block, "target_min");
      sched.has_target_max = pin->DoesParameterExist(block, "target_max");
      if (sched.has_target_min) sched.target_min = pin->GetReal(block, "target_min");
      if (sched.has_target_max) sched.target_max = pin->GetReal(block, "target_max");
    }

    if (sched.count_per_event < 0) {
      FatalParticleInput(block + ": count_per_event must be non-negative");
    }
    seed_schedules_.push_back(sched);
    fallback_id++;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void Particles::AppendParticles

void Particles::AppendParticles(const HostArray2D<Real> &new_rdata,
                                const HostArray2D<int> &new_idata, int nnew) {
  if (nnew <= 0) return;

  int nold = nprtcl_thispack;
  HostArray2D<Real> old_r("old_prtcl_rdata", nrdata, nold);
  HostArray2D<int> old_i("old_prtcl_idata", nidata, nold);
  if (nold > 0) {
    Kokkos::deep_copy(old_r, prtcl_rdata);
    Kokkos::deep_copy(old_i, prtcl_idata);
  }

  HostArray2D<Real> merged_r("merged_prtcl_rdata", nrdata, nold + nnew);
  HostArray2D<int> merged_i("merged_prtcl_idata", nidata, nold + nnew);
  for (int p=0; p<nold; ++p) {
    for (int n=0; n<nrdata; ++n) merged_r(n,p) = old_r(n,p);
    for (int n=0; n<nidata; ++n) merged_i(n,p) = old_i(n,p);
  }
  for (int p=0; p<nnew; ++p) {
    for (int n=0; n<nrdata; ++n) merged_r(n,nold+p) = new_rdata(n,p);
    for (int n=0; n<nidata; ++n) merged_i(n,nold+p) = new_idata(n,p);
  }

  nprtcl_thispack = nold + nnew;
  Kokkos::realloc(prtcl_rdata, nrdata, nprtcl_thispack);
  Kokkos::realloc(prtcl_idata, nidata, nprtcl_thispack);
  Kokkos::deep_copy(prtcl_rdata, merged_r);
  Kokkos::deep_copy(prtcl_idata, merged_i);
  pmy_pack->pmesh->UpdateParticleCounts();
}

//----------------------------------------------------------------------------------------
//! \fn void Particles::SeedInitialTracers

void Particles::SeedInitialTracers() {
  if (particle_type != ParticleType::lagrangian_mc) return;
  SeedTracersAtTime(pmy_pack->pmesh->time, true);
  SampleThermodynamics();
}

//----------------------------------------------------------------------------------------
//! \fn void Particles::SeedTracersAtTime

void Particles::SeedTracersAtTime(Real event_time, bool initial_only) {
  if (particle_type != ParticleType::lagrangian_mc) return;

  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb = pmy_pack->nmb_thispack;

  pmy_pack->pmb->mb_size.template sync<HostMemSpace>();
  pmy_pack->pmb->mb_lev.template sync<HostMemSpace>();
  auto h_size = pmy_pack->pmb->mb_size.h_view;
  auto h_lev = pmy_pack->pmb->mb_lev.h_view;

  int nfluid = 0, nscalars = 0;
  EOS_Data eos;
  if (pmy_pack->phydro != nullptr) {
    nfluid = pmy_pack->phydro->nhydro;
    nscalars = pmy_pack->phydro->nscalars;
    eos = pmy_pack->phydro->peos->eos_data;
  } else {
    nfluid = pmy_pack->pmhd->nmhd;
    nscalars = pmy_pack->pmhd->nscalars;
    eos = pmy_pack->pmhd->peos->eos_data;
  }
  int ncells1 = indcs.nx1 + 2*indcs.ng;
  int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*indcs.ng) : 1;
  int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*indcs.ng) : 1;
  int nmb_alloc = std::max(nmb, pm->nmb_maxperrank);
  HostArray5D<Real> h_w0("seed_w0", nmb_alloc, nfluid+nscalars,
                         ncells3, ncells2, ncells1);
  if (pmy_pack->phydro != nullptr) {
    Kokkos::deep_copy(h_w0, pmy_pack->phydro->w0);
  } else {
    Kokkos::deep_copy(h_w0, pmy_pack->pmhd->w0);
  }

  constexpr Real eps = 64.0*std::numeric_limits<Real>::epsilon();
  for (auto &sched : seed_schedules_) {
    if (sched.complete || sched.count_per_event == 0) continue;
    if (initial_only && std::abs(sched.next_time - pm->time) > eps) continue;
    if (!initial_only && sched.next_time > event_time + eps) continue;

    while (!sched.complete && sched.next_time <= event_time + eps) {
      std::vector<CandidateCell> cells;
      Real local_weight = 0.0;

      for (int m=0; m<nmb; ++m) {
        for (int k=ks; k<=ke; ++k) {
          Real x3 = CellCenterX(k-ks, nx3, h_size(m).x3min, h_size(m).x3max);
          for (int j=js; j<=je; ++j) {
            Real x2 = CellCenterX(j-js, nx2, h_size(m).x2min, h_size(m).x2max);
            for (int i=is; i<=ie; ++i) {
              Real x1 = CellCenterX(i-is, nx1, h_size(m).x1min, h_size(m).x1max);

              bool inside = true;
              if (sched.region == TracerSeedRegion::box) {
                inside = (x1 >= sched.x1min && x1 <= sched.x1max &&
                          x2 >= sched.x2min && x2 <= sched.x2max &&
                          x3 >= sched.x3min && x3 <= sched.x3max);
              } else if (sched.region == TracerSeedRegion::sphere) {
                inside = (SQR(x1 - sched.center1) + SQR(x2 - sched.center2) +
                          SQR(x3 - sched.center3)) <= SQR(sched.radius);
              } else if (sched.region == TracerSeedRegion::slab) {
                Real coord = (sched.slab_axis == 1) ? x1 :
                             ((sched.slab_axis == 2) ? x2 : x3);
                inside = (coord >= sched.slab_min && coord <= sched.slab_max);
              }
              if (!inside) continue;

              ThermoSample sample = SampleFromHost(h_w0, eos, m, k, j, i);
              Real target_value = 0.0;
              if (sched.target == TracerSeedTarget::density) {
                target_value = sample.rho;
              } else if (sched.target == TracerSeedTarget::temperature) {
                target_value = sample.temperature;
              } else if (sched.target == TracerSeedTarget::pressure) {
                target_value = sample.pressure;
              } else if (sched.target == TracerSeedTarget::entropy) {
                target_value = sample.entropy;
              } else if (sched.target == TracerSeedTarget::scalar) {
                target_value = h_w0(m,nfluid+sched.scalar_index,k,j,i);
              }
              if (sched.target != TracerSeedTarget::none) {
                if (sched.has_target_min && target_value < sched.target_min) continue;
                if (sched.has_target_max && target_value > sched.target_max) continue;
              }

              Real volume = h_size(m).dx1*h_size(m).dx2*h_size(m).dx3;
              Real weight = (sched.weight == TracerSeedWeight::mass) ?
                            sample.rho*volume : volume;
              if (weight <= 0.0) continue;
              local_weight += weight;
              cells.push_back({m,k,j,i,local_weight});
            }
          }
        }
      }

      Real global_weight = local_weight;
      Real rank_weight_start = 0.0;
#if MPI_PARALLEL_ENABLED
      MPI_Allreduce(MPI_IN_PLACE, &global_weight, 1, MPI_ATHENA_REAL, MPI_SUM,
                    MPI_COMM_WORLD);
      MPI_Exscan(&local_weight, &rank_weight_start, 1, MPI_ATHENA_REAL, MPI_SUM,
                 MPI_COMM_WORLD);
      if (global_variable::my_rank == 0) rank_weight_start = 0.0;
#endif
      if (global_weight <= 0.0) {
        if (global_variable::my_rank == 0) {
          std::cout << "### WARNING: tracer seed schedule " << sched.id
                    << " has no eligible cells at time " << sched.next_time << std::endl;
        }
        sched.complete = (sched.cadence <= 0.0 || sched.next_time >= sched.end_time);
        if (!sched.complete) sched.next_time += sched.cadence;
        continue;
      }

      std::vector<int> chosen;
      std::vector<int> chosen_tag;
      for (int q=0; q<sched.count_per_event; ++q) {
        std::uint64_t key = static_cast<std::uint64_t>(sched.seed)
          ^ (static_cast<std::uint64_t>(sched.id) << 32)
          ^ (static_cast<std::uint64_t>(sched.event_index) << 16)
          ^ static_cast<std::uint64_t>(q);
        Real draw = HashReal(key)*global_weight;
        if (draw < rank_weight_start ||
            draw >= rank_weight_start + local_weight || cells.empty()) {
          continue;
        }
        Real local_draw = draw - rank_weight_start;
        auto it = std::lower_bound(cells.begin(), cells.end(), local_draw,
          [](const CandidateCell &cell, Real value) {
            return cell.cumulative_end < value;
          });
        if (it == cells.end()) it = cells.end() - 1;
        chosen.push_back(static_cast<int>(it - cells.begin()));
        chosen_tag.push_back(static_cast<int>(next_tracer_tag + q));
      }

      int nnew = chosen.size();
      HostArray2D<Real> new_r("new_lagrangian_mc_rdata", nrdata, nnew);
      HostArray2D<int> new_i("new_lagrangian_mc_idata", nidata, nnew);
      for (int p=0; p<nnew; ++p) {
        CandidateCell cell = cells[chosen[p]];
        int m = cell.m, k = cell.k, j = cell.j, i = cell.i;
        ThermoSample sample = SampleFromHost(h_w0, eos, m, k, j, i);
        for (int n=0; n<nrdata; ++n) new_r(n,p) = 0.0;
        for (int n=0; n<nidata; ++n) new_i(n,p) = 0;

        new_i(PGID,p) = pmy_pack->gids + m;
        new_i(PTAG,p) = chosen_tag[p];
        new_i(PLASTMOVE,p) = 0;
        new_i(PLASTLEVEL,p) = h_lev(m);
        new_i(PSEEDID,p) = sched.id;

        new_r(LMCX,p) = CellCenterX(i-is, nx1, h_size(m).x1min, h_size(m).x1max);
        new_r(LMCY,p) = CellCenterX(j-js, nx2, h_size(m).x2min, h_size(m).x2max);
        new_r(LMCZ,p) = CellCenterX(k-ks, nx3, h_size(m).x3min, h_size(m).x3max);
        new_r(LMC_RHO,p) = sample.rho;
        new_r(LMC_PRESSURE,p) = sample.pressure;
        new_r(LMC_TEMPERATURE,p) = sample.temperature;
        new_r(LMC_ENTROPY,p) = sample.entropy;
        new_r(LMC_EINT,p) = sample.eint;
        new_r(LMC_VX,p) = sample.v1;
        new_r(LMC_VY,p) = sample.v2;
        new_r(LMC_VZ,p) = sample.v3;
        new_r(LMC_CREATE_TIME,p) = sched.next_time;
        for (int n=0; n<nscalars; ++n) {
          new_r(LMC_SCALAR0+n,p) = h_w0(m,nfluid+n,k,j,i);
        }
      }
      AppendParticles(new_r, new_i, nnew);

      sched.event_index++;
      next_tracer_tag += sched.count_per_event;
      if (sched.cadence <= 0.0 || sched.next_time + sched.cadence > sched.end_time + eps) {
        sched.complete = true;
      } else {
        sched.next_time += sched.cadence;
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::PushLagrangianMC

TaskStatus Particles::PushLagrangianMC(Driver *pdriver, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  bool &multi_d = pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_pack->pmesh->three_d;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto &mblev = pmy_pack->pmb->mb_lev;
  auto &pi = prtcl_idata;
  auto &pr = prtcl_rdata;
  auto &gids = pmy_pack->gids;
  auto &u1_ = (pmy_pack->phydro != nullptr) ? pmy_pack->phydro->u1 : pmy_pack->pmhd->u1;
  auto &uflxidn_ = (pmy_pack->phydro != nullptr) ?
                   pmy_pack->phydro->uflxidnsaved : pmy_pack->pmhd->uflxidnsaved;
  auto &flx1_ = uflxidn_.x1f;
  auto &flx2_ = uflxidn_.x2f;
  auto &flx3_ = uflxidn_.x3f;
  int ncycle = pmy_pack->pmesh->ncycle;
  std::int64_t rseed = random_seed;
  int nmb = pmy_pack->nmb_thispack;

  par_for("lagrangian_mc_push", DevExeSpace(), 0, (nprtcl_thispack-1),
  KOKKOS_LAMBDA(const int p) {
    if (pi(PLASTMOVE,p) < 0) return;
    int m = pi(PGID,p) - gids;
    if (m < 0 || m >= nmb) {
      pi(PLASTMOVE,p) = -1;
      return;
    }

    int ip = static_cast<int>((pr(LMCX,p) - mbsize.d_view(m).x1min)/
                              mbsize.d_view(m).dx1) + is;
    int jp = js;
    int kp = ks;
    if (multi_d) {
      jp = static_cast<int>((pr(LMCY,p) - mbsize.d_view(m).x2min)/
                            mbsize.d_view(m).dx2) + js;
    }
    if (three_d) {
      kp = static_cast<int>((pr(LMCZ,p) - mbsize.d_view(m).x3min)/
                            mbsize.d_view(m).dx3) + ks;
    }
    int ie = is + indcs.nx1;
    int je = js + indcs.nx2;
    int ke = ks + indcs.nx3;
    if (ip < is || ip >= ie || jp < js || jp >= je || kp < ks || kp >= ke) return;

    Real mass = u1_(m,IDN,kp,jp,ip);
    if (mass <= 0.0) return;
    Real flx1_left  = fmax(-flx1_(m,kp,jp,ip)/mass, 0.0);
    Real flx1_right = fmax( flx1_(m,kp,jp,ip+1)/mass, 0.0);
    Real flx2_left  = multi_d ? fmax(-flx2_(m,kp,jp,ip)/mass, 0.0) : 0.0;
    Real flx2_right = multi_d ? fmax( flx2_(m,kp,jp+1,ip)/mass, 0.0) : 0.0;
    Real flx3_left  = three_d ? fmax(-flx3_(m,kp,jp,ip)/mass, 0.0) : 0.0;
    Real flx3_right = three_d ? fmax( flx3_(m,kp+1,jp,ip)/mass, 0.0) : 0.0;

    std::uint64_t key = static_cast<std::uint64_t>(pi(PTAG,p))*7919ULL
                      + static_cast<std::uint64_t>(ncycle)*104729ULL
                      + static_cast<std::uint64_t>(rseed);
    Real draw = HashReal(key);

    pi(PLASTLEVEL,p) = mblev.d_view(m);
    pi(PLASTMOVE,p) = 32*(ip % 2) + 16*(jp % 2) + 8*(kp % 2);

    if (draw < flx1_left) {
      pr(LMCX,p) -= mbsize.d_view(m).dx1;
      pi(PLASTMOVE,p) += 1;
    } else if (draw < flx1_left + flx1_right) {
      pr(LMCX,p) += mbsize.d_view(m).dx1;
      pi(PLASTMOVE,p) += 2;
    } else if (multi_d && draw < flx1_left + flx1_right + flx2_left) {
      pr(LMCY,p) -= mbsize.d_view(m).dx2;
      pi(PLASTMOVE,p) += 3;
    } else if (multi_d && draw < flx1_left + flx1_right + flx2_left + flx2_right) {
      pr(LMCY,p) += mbsize.d_view(m).dx2;
      pi(PLASTMOVE,p) += 4;
    } else if (three_d && draw < flx1_left + flx1_right + flx2_left + flx2_right +
                                flx3_left) {
      pr(LMCZ,p) -= mbsize.d_view(m).dx3;
      pi(PLASTMOVE,p) += 5;
    } else if (three_d && draw < flx1_left + flx1_right + flx2_left + flx2_right +
                                flx3_left + flx3_right) {
      pr(LMCZ,p) += mbsize.d_view(m).dx3;
      pi(PLASTMOVE,p) += 6;
    }
  });

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::AdjustMeshRefinement

TaskStatus Particles::AdjustMeshRefinement(Driver *pdriver, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  auto &pi = prtcl_idata;
  auto &pr = prtcl_rdata;
  auto &gids = pmy_pack->gids;
  auto &mblev = pmy_pack->pmb->mb_lev;
  bool &multi_d = pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_pack->pmesh->three_d;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto &uflxidn_ = (pmy_pack->phydro != nullptr) ?
                   pmy_pack->phydro->uflxidnsaved : pmy_pack->pmhd->uflxidnsaved;
  auto &flx1_ = uflxidn_.x1f;
  auto &flx2_ = uflxidn_.x2f;
  auto &flx3_ = uflxidn_.x3f;
  int ncycle = pmy_pack->pmesh->ncycle;
  std::int64_t rseed = random_seed;
  int nmb = pmy_pack->nmb_thispack;

  par_for("lagrangian_mc_amr_adjust", DevExeSpace(), 0, (nprtcl_thispack-1),
  KOKKOS_LAMBDA(const int p) {
    if (pi(PLASTMOVE,p) < 0) return;
    int m = pi(PGID,p) - gids;
    if (m < 0 || m >= nmb) {
      pi(PLASTMOVE,p) = -1;
      return;
    }

    int level = mblev.d_view(m);
    int lastlevel = pi(PLASTLEVEL,p);
    int lastmove = pi(PLASTMOVE,p);
    int i_parity = lastmove/32;
    int j_parity = (lastmove % 32)/16;
    int k_parity = (lastmove % 16)/8;
    lastmove = lastmove % 8;

    Real dx1 = mbsize.d_view(m).dx1;
    Real dx2 = multi_d ? mbsize.d_view(m).dx2 : 0.0;
    Real dx3 = three_d ? mbsize.d_view(m).dx3 : 0.0;

    if (level > lastlevel) {
      if (lastmove == 1) {
        pr(LMCX,p) += dx1/2; pr(LMCY,p) -= dx2/2; pr(LMCZ,p) -= dx3/2;
      } else if (lastmove == 2) {
        pr(LMCX,p) -= dx1/2; pr(LMCY,p) -= dx2/2; pr(LMCZ,p) -= dx3/2;
      } else if (multi_d && lastmove == 3) {
        pr(LMCY,p) += dx2/2; pr(LMCX,p) -= dx1/2; pr(LMCZ,p) -= dx3/2;
      } else if (multi_d && lastmove == 4) {
        pr(LMCY,p) -= dx2/2; pr(LMCX,p) -= dx1/2; pr(LMCZ,p) -= dx3/2;
      } else if (three_d && lastmove == 5) {
        pr(LMCZ,p) += dx3/2; pr(LMCX,p) -= dx1/2; pr(LMCY,p) -= dx2/2;
      } else if (three_d && lastmove == 6) {
        pr(LMCZ,p) -= dx3/2; pr(LMCX,p) -= dx1/2; pr(LMCY,p) -= dx2/2;
      }

      int ip = static_cast<int>((pr(LMCX,p) - mbsize.d_view(m).x1min)/dx1) + is;
      int jp = multi_d ? static_cast<int>((pr(LMCY,p) - mbsize.d_view(m).x2min)/dx2) + js
                       : js;
      int kp = three_d ? static_cast<int>((pr(LMCZ,p) - mbsize.d_view(m).x3min)/dx3) + ks
                       : ks;

      Real flx[4] = {0.0, 0.0, 0.0, 0.0};
      if (lastmove == 1) {
        flx[0] = -flx1_(m,kp,jp,ip+1);
        flx[1] = multi_d ? -flx1_(m,kp,jp+1,ip+1) : 0.0;
        flx[2] = three_d ? -flx1_(m,kp+1,jp,ip+1) : 0.0;
        flx[3] = (multi_d && three_d) ? -flx1_(m,kp+1,jp+1,ip+1) : 0.0;
      } else if (lastmove == 2) {
        flx[0] = flx1_(m,kp,jp,ip);
        flx[1] = multi_d ? flx1_(m,kp,jp+1,ip) : 0.0;
        flx[2] = three_d ? flx1_(m,kp+1,jp,ip) : 0.0;
        flx[3] = (multi_d && three_d) ? flx1_(m,kp+1,jp+1,ip) : 0.0;
      } else if (lastmove == 3) {
        flx[0] = -flx2_(m,kp,jp+1,ip);
        flx[1] = -flx2_(m,kp,jp+1,ip+1);
        flx[2] = three_d ? -flx2_(m,kp+1,jp+1,ip) : 0.0;
        flx[3] = three_d ? -flx2_(m,kp+1,jp+1,ip+1) : 0.0;
      } else if (lastmove == 4) {
        flx[0] = flx2_(m,kp,jp,ip);
        flx[1] = flx2_(m,kp,jp,ip+1);
        flx[2] = three_d ? flx2_(m,kp+1,jp,ip) : 0.0;
        flx[3] = three_d ? flx2_(m,kp+1,jp,ip+1) : 0.0;
      } else if (lastmove == 5) {
        flx[0] = -flx3_(m,kp+1,jp,ip);
        flx[1] = -flx3_(m,kp+1,jp,ip+1);
        flx[2] = -flx3_(m,kp+1,jp+1,ip);
        flx[3] = -flx3_(m,kp+1,jp+1,ip+1);
      } else if (lastmove == 6) {
        flx[0] = flx3_(m,kp,jp,ip);
        flx[1] = flx3_(m,kp,jp,ip+1);
        flx[2] = flx3_(m,kp,jp+1,ip);
        flx[3] = flx3_(m,kp,jp+1,ip+1);
      }
      Real total = 0.0;
      for (int n=0; n<4; ++n) {
        flx[n] = fmax(flx[n], 0.0);
        total += flx[n];
      }
      if (total <= 0.0) total = 1.0;
      for (int n=0; n<4; ++n) flx[n] /= total;

      std::uint64_t key = static_cast<std::uint64_t>(pi(PTAG,p))*7919ULL
                        + static_cast<std::uint64_t>(ncycle)*104729ULL
                        + static_cast<std::uint64_t>(rseed) + 1ULL;
      Real draw = HashReal(key);
      int target = 4;
      if (draw < flx[0]) target = 1;
      else if (draw < flx[0] + flx[1]) target = 2;
      else if (draw < flx[0] + flx[1] + flx[2]) target = 3;

      if (lastmove == 1 || lastmove == 2) {
        if (target == 2) pr(LMCY,p) += dx2;
        else if (target == 3) pr(LMCZ,p) += dx3;
        else if (target == 4) {pr(LMCY,p) += dx2; pr(LMCZ,p) += dx3;}
      } else if (lastmove == 3 || lastmove == 4) {
        if (target == 2) pr(LMCX,p) += dx1;
        else if (target == 3) pr(LMCZ,p) += dx3;
        else if (target == 4) {pr(LMCX,p) += dx1; pr(LMCZ,p) += dx3;}
      } else if (lastmove == 5 || lastmove == 6) {
        if (target == 2) pr(LMCX,p) += dx1;
        else if (target == 3) pr(LMCY,p) += dx2;
        else if (target == 4) {pr(LMCX,p) += dx1; pr(LMCY,p) += dx2;}
      }
    } else if (level < lastlevel) {
      pr(LMCX,p) += (i_parity ? -0.25 : 0.25)*dx1;
      if (multi_d) pr(LMCY,p) += (j_parity ? -0.25 : 0.25)*dx2;
      if (three_d) pr(LMCZ,p) += (k_parity ? -0.25 : 0.25)*dx3;
      if (lastmove == 1) pr(LMCX,p) -= dx1/2;
      else if (lastmove == 2) pr(LMCX,p) += dx1/2;
      else if (lastmove == 3) pr(LMCY,p) -= dx2/2;
      else if (lastmove == 4) pr(LMCY,p) += dx2/2;
      else if (lastmove == 5) pr(LMCZ,p) -= dx3/2;
      else if (lastmove == 6) pr(LMCZ,p) += dx3/2;
    }
  });

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void Particles::RemapAfterMeshRefinement

void Particles::RemapAfterMeshRefinement() {
  if (particle_type != ParticleType::lagrangian_mc) return;

  Mesh *pm = pmy_pack->pmesh;
  Kokkos::fence();
  HostArray2D<Real> hr("amr_prtcl_rdata_old", nrdata, nprtcl_thispack);
  HostArray2D<int> hi("amr_prtcl_idata_old", nidata, nprtcl_thispack);
  if (nprtcl_thispack > 0) {
    Kokkos::deep_copy(hr, prtcl_rdata);
    Kokkos::deep_copy(hi, prtcl_idata);
  }

  std::vector<int> dest_rank(nprtcl_thispack, global_variable::my_rank);
  for (int p=0; p<nprtcl_thispack; ++p) {
    hr(LMCX,p) = WrapOrClamp(hr(LMCX,p), pm->mesh_size.x1min, pm->mesh_size.x1max,
                             pm->mesh_bcs[BoundaryFace::inner_x1] == BoundaryFlag::periodic ||
                             pm->mesh_bcs[BoundaryFace::outer_x1] == BoundaryFlag::periodic);
    hr(LMCY,p) = WrapOrClamp(hr(LMCY,p), pm->mesh_size.x2min, pm->mesh_size.x2max,
                             pm->mesh_bcs[BoundaryFace::inner_x2] == BoundaryFlag::periodic ||
                             pm->mesh_bcs[BoundaryFace::outer_x2] == BoundaryFlag::periodic);
    hr(LMCZ,p) = WrapOrClamp(hr(LMCZ,p), pm->mesh_size.x3min, pm->mesh_size.x3max,
                             pm->mesh_bcs[BoundaryFace::inner_x3] == BoundaryFlag::periodic ||
                             pm->mesh_bcs[BoundaryFace::outer_x3] == BoundaryFlag::periodic);

    int gid = LocateMeshBlockGID(pm, hr(LMCX,p), hr(LMCY,p), hr(LMCZ,p));
    if (gid < 0) {
      FatalParticleInput("could not locate tracer particle after mesh refinement");
    }
    SnapToMeshBlockCellCenter(pm, gid, hr(LMCX,p), hr(LMCY,p), hr(LMCZ,p));
    hi(PGID,p) = gid;
    hi(PLASTMOVE,p) = 0;
    hi(PLASTLEVEL,p) = pm->lloc_eachmb[gid].level;
    dest_rank[p] = pm->rank_eachmb[gid];
  }

#if MPI_PARALLEL_ENABLED
  int nranks = global_variable::nranks;
  std::vector<int> send_particles(nranks, 0), recv_particles(nranks, 0);
  for (int p=0; p<nprtcl_thispack; ++p) send_particles[dest_rank[p]]++;
  MPI_Alltoall(send_particles.data(), 1, MPI_INT, recv_particles.data(), 1, MPI_INT,
               MPI_COMM_WORLD);

  std::vector<int> send_real_count(nranks), recv_real_count(nranks);
  std::vector<int> send_int_count(nranks), recv_int_count(nranks);
  std::vector<int> send_real_disp(nranks, 0), recv_real_disp(nranks, 0);
  std::vector<int> send_int_disp(nranks, 0), recv_int_disp(nranks, 0);
  int nrecv = 0;
  for (int r=0; r<nranks; ++r) {
    send_real_count[r] = send_particles[r]*nrdata;
    recv_real_count[r] = recv_particles[r]*nrdata;
    send_int_count[r] = send_particles[r]*nidata;
    recv_int_count[r] = recv_particles[r]*nidata;
    if (r > 0) {
      send_real_disp[r] = send_real_disp[r-1] + send_real_count[r-1];
      recv_real_disp[r] = recv_real_disp[r-1] + recv_real_count[r-1];
      send_int_disp[r] = send_int_disp[r-1] + send_int_count[r-1];
      recv_int_disp[r] = recv_int_disp[r-1] + recv_int_count[r-1];
    }
    nrecv += recv_particles[r];
  }

  std::vector<Real> send_reals(send_real_disp[nranks-1] + send_real_count[nranks-1]);
  std::vector<Real> recv_reals(recv_real_disp[nranks-1] + recv_real_count[nranks-1]);
  std::vector<int> send_ints(send_int_disp[nranks-1] + send_int_count[nranks-1]);
  std::vector<int> recv_ints(recv_int_disp[nranks-1] + recv_int_count[nranks-1]);

  std::vector<int> fill(nranks, 0);
  for (int p=0; p<nprtcl_thispack; ++p) {
    int r = dest_rank[p];
    int q = fill[r]++;
    int rbase = send_real_disp[r] + q*nrdata;
    int ibase = send_int_disp[r] + q*nidata;
    for (int n=0; n<nrdata; ++n) send_reals[rbase+n] = hr(n,p);
    for (int n=0; n<nidata; ++n) send_ints[ibase+n] = hi(n,p);
  }

  MPI_Alltoallv(send_reals.data(), send_real_count.data(), send_real_disp.data(),
                MPI_ATHENA_REAL, recv_reals.data(), recv_real_count.data(),
                recv_real_disp.data(), MPI_ATHENA_REAL, MPI_COMM_WORLD);
  MPI_Alltoallv(send_ints.data(), send_int_count.data(), send_int_disp.data(), MPI_INT,
                recv_ints.data(), recv_int_count.data(), recv_int_disp.data(), MPI_INT,
                MPI_COMM_WORLD);

  nprtcl_thispack = nrecv;
  Kokkos::realloc(prtcl_rdata, nrdata, nprtcl_thispack);
  Kokkos::realloc(prtcl_idata, nidata, nprtcl_thispack);
  HostArray2D<Real> new_r("amr_prtcl_rdata_new", nrdata, nprtcl_thispack);
  HostArray2D<int> new_i("amr_prtcl_idata_new", nidata, nprtcl_thispack);
  for (int p=0; p<nprtcl_thispack; ++p) {
    for (int n=0; n<nrdata; ++n) new_r(n,p) = recv_reals[p*nrdata+n];
    for (int n=0; n<nidata; ++n) new_i(n,p) = recv_ints[p*nidata+n];
  }
#else
  int nlocal = 0;
  for (int p=0; p<nprtcl_thispack; ++p) {
    if (dest_rank[p] == global_variable::my_rank) nlocal++;
  }
  HostArray2D<Real> new_r("amr_prtcl_rdata_new", nrdata, nlocal);
  HostArray2D<int> new_i("amr_prtcl_idata_new", nidata, nlocal);
  int q = 0;
  for (int p=0; p<nprtcl_thispack; ++p) {
    if (dest_rank[p] != global_variable::my_rank) continue;
    for (int n=0; n<nrdata; ++n) new_r(n,q) = hr(n,p);
    for (int n=0; n<nidata; ++n) new_i(n,q) = hi(n,p);
    q++;
  }
  nprtcl_thispack = nlocal;
  Kokkos::realloc(prtcl_rdata, nrdata, nprtcl_thispack);
  Kokkos::realloc(prtcl_idata, nidata, nprtcl_thispack);
#endif

  if (nprtcl_thispack > 0) {
    Kokkos::deep_copy(prtcl_rdata, new_r);
    Kokkos::deep_copy(prtcl_idata, new_i);
  }
  pm->UpdateParticleCounts();
}

//----------------------------------------------------------------------------------------
//! \fn void Particles::SampleThermodynamics

void Particles::SampleThermodynamics() {
  if (particle_type != ParticleType::lagrangian_mc || nprtcl_thispack <= 0) return;

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int nmb = pmy_pack->nmb_thispack;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto &pi = prtcl_idata;
  auto &pr = prtcl_rdata;
  auto &gids = pmy_pack->gids;
  bool &multi_d = pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_pack->pmesh->three_d;
  int nscalars = GetLagrangianMCScalarCount();

  if (pmy_pack->phydro != nullptr) {
    auto w0 = pmy_pack->phydro->w0;
    auto eos = pmy_pack->phydro->peos->eos_data;
    int nfluid = pmy_pack->phydro->nhydro;
    par_for("lagrangian_mc_sample_hydro", DevExeSpace(), 0, nprtcl_thispack-1,
    KOKKOS_LAMBDA(const int p) {
      int m = pi(PGID,p) - gids;
      if (m < 0 || m >= nmb) return;
      int i = static_cast<int>((pr(LMCX,p) - mbsize.d_view(m).x1min)/
                               mbsize.d_view(m).dx1) + is;
      int j = multi_d ? static_cast<int>((pr(LMCY,p) - mbsize.d_view(m).x2min)/
                                       mbsize.d_view(m).dx2) + js : js;
      int k = three_d ? static_cast<int>((pr(LMCZ,p) - mbsize.d_view(m).x3min)/
                                       mbsize.d_view(m).dx3) + ks : ks;
      i = i < is ? is : (i > is + indcs.nx1 - 1 ? is + indcs.nx1 - 1 : i);
      j = j < js ? js : (j > js + indcs.nx2 - 1 ? js + indcs.nx2 - 1 : j);
      k = k < ks ? ks : (k > ks + indcs.nx3 - 1 ? ks + indcs.nx3 - 1 : k);
      Real rho = w0(m,IDN,k,j,i);
      Real eint = w0(m,IEN,k,j,i);
      Real press = eos.IdealGasPressure(eint);
      pr(LMC_RHO,p) = rho;
      pr(LMC_PRESSURE,p) = press;
      pr(LMC_TEMPERATURE,p) = press/rho;
      pr(LMC_ENTROPY,p) = log(press/pow(rho, eos.gamma));
      pr(LMC_EINT,p) = eint;
      pr(LMC_VX,p) = w0(m,IVX,k,j,i);
      pr(LMC_VY,p) = w0(m,IVY,k,j,i);
      pr(LMC_VZ,p) = w0(m,IVZ,k,j,i);
      for (int n=0; n<nscalars; ++n) pr(LMC_SCALAR0+n,p) = w0(m,nfluid+n,k,j,i);
    });
  } else {
    auto w0 = pmy_pack->pmhd->w0;
    auto eos = pmy_pack->pmhd->peos->eos_data;
    int nfluid = pmy_pack->pmhd->nmhd;
    par_for("lagrangian_mc_sample_mhd", DevExeSpace(), 0, nprtcl_thispack-1,
    KOKKOS_LAMBDA(const int p) {
      int m = pi(PGID,p) - gids;
      if (m < 0 || m >= nmb) return;
      int i = static_cast<int>((pr(LMCX,p) - mbsize.d_view(m).x1min)/
                               mbsize.d_view(m).dx1) + is;
      int j = multi_d ? static_cast<int>((pr(LMCY,p) - mbsize.d_view(m).x2min)/
                                       mbsize.d_view(m).dx2) + js : js;
      int k = three_d ? static_cast<int>((pr(LMCZ,p) - mbsize.d_view(m).x3min)/
                                       mbsize.d_view(m).dx3) + ks : ks;
      i = i < is ? is : (i > is + indcs.nx1 - 1 ? is + indcs.nx1 - 1 : i);
      j = j < js ? js : (j > js + indcs.nx2 - 1 ? js + indcs.nx2 - 1 : j);
      k = k < ks ? ks : (k > ks + indcs.nx3 - 1 ? ks + indcs.nx3 - 1 : k);
      Real rho = w0(m,IDN,k,j,i);
      Real eint = w0(m,IEN,k,j,i);
      Real press = eos.IdealGasPressure(eint);
      pr(LMC_RHO,p) = rho;
      pr(LMC_PRESSURE,p) = press;
      pr(LMC_TEMPERATURE,p) = press/rho;
      pr(LMC_ENTROPY,p) = log(press/pow(rho, eos.gamma));
      pr(LMC_EINT,p) = eint;
      pr(LMC_VX,p) = w0(m,IVX,k,j,i);
      pr(LMC_VY,p) = w0(m,IVY,k,j,i);
      pr(LMC_VZ,p) = w0(m,IVZ,k,j,i);
      for (int n=0; n<nscalars; ++n) pr(LMC_SCALAR0+n,p) = w0(m,nfluid+n,k,j,i);
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void Particles::WriteRestartData

void Particles::WriteRestartData(IOWrapper &resfile, bool single_file_per_rank) {
  struct ParticleRestartHeader {
    char magic[16];
    int version;
    int enabled;
    int nrdata;
    int nidata;
    int nlocal;
    int nschedules;
    std::int64_t next_tag;
  };

  ParticleRestartHeader header;
  std::memset(&header, 0, sizeof(header));
  std::strncpy(header.magic, "ATHKPRTCLMC", sizeof(header.magic)-1);
  header.version = 1;
  header.enabled = IsLagrangianMC() ? 1 : 0;
  header.nrdata = nrdata;
  header.nidata = nidata;
  header.nlocal = nprtcl_thispack;
  header.nschedules = static_cast<int>(seed_schedules_.size());
  header.next_tag = next_tracer_tag;

  resfile.Write_any_type(&header, sizeof(header), "byte", single_file_per_rank);
  if (!header.enabled) return;

  std::vector<Real> sched_real(1*header.nschedules);
  std::vector<int> sched_int(3*header.nschedules);
  for (int n=0; n<header.nschedules; ++n) {
    sched_real[n] = seed_schedules_[n].next_time;
    sched_int[3*n] = seed_schedules_[n].id;
    sched_int[3*n+1] = seed_schedules_[n].event_index;
    sched_int[3*n+2] = seed_schedules_[n].complete ? 1 : 0;
  }
  if (header.nschedules > 0) {
    resfile.Write_any_type(sched_real.data(), sched_real.size(), "Real",
                           single_file_per_rank);
    resfile.Write_any_type(sched_int.data(), sched_int.size(), "int", single_file_per_rank);
  }

  HostArray2D<Real> hr("rst_prtcl_rdata", nrdata, nprtcl_thispack);
  HostArray2D<int> hi("rst_prtcl_idata", nidata, nprtcl_thispack);
  if (nprtcl_thispack > 0) {
    Kokkos::deep_copy(hr, prtcl_rdata);
    Kokkos::deep_copy(hi, prtcl_idata);
    resfile.Write_any_type(hr.data(), nrdata*nprtcl_thispack, "Real", single_file_per_rank);
    resfile.Write_any_type(hi.data(), nidata*nprtcl_thispack, "int", single_file_per_rank);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void Particles::ReadRestartData

void Particles::ReadRestartData(IOWrapper &resfile, bool single_file_per_rank) {
  struct ParticleRestartHeader {
    char magic[16];
    int version;
    int enabled;
    int nrdata;
    int nidata;
    int nlocal;
    int nschedules;
    std::int64_t next_tag;
  };

  ParticleRestartHeader header;
  if (resfile.Read_bytes(&header, sizeof(header), 1, single_file_per_rank) != 1) {
    FatalParticleInput("particle restart section is missing or truncated");
  }
  if (std::strncmp(header.magic, "ATHKPRTCLMC", 11) != 0 || header.version != 1) {
    FatalParticleInput("particle restart section has an unrecognized format");
  }
  if (!header.enabled) return;
  if (header.nrdata != nrdata || header.nidata != nidata) {
    FatalParticleInput("particle restart array dimensions do not match this executable");
  }

  std::vector<Real> sched_real(header.nschedules);
  std::vector<int> sched_int(3*header.nschedules);
  if (header.nschedules > 0) {
    if (resfile.Read_Reals(sched_real.data(), sched_real.size(), single_file_per_rank) !=
        sched_real.size()) {
      FatalParticleInput("particle restart schedule real data is truncated");
    }
    if (resfile.Read_bytes(sched_int.data(), sizeof(int), sched_int.size(),
                           single_file_per_rank) != sched_int.size()) {
      FatalParticleInput("particle restart schedule integer data is truncated");
    }
  }

  for (int n=0; n<header.nschedules; ++n) {
    for (auto &sched : seed_schedules_) {
      if (sched.id == sched_int[3*n]) {
        sched.next_time = sched_real[n];
        sched.event_index = sched_int[3*n+1];
        sched.complete = (sched_int[3*n+2] != 0);
      }
    }
  }
  next_tracer_tag = header.next_tag;
  nprtcl_thispack = header.nlocal;
  Kokkos::realloc(prtcl_rdata, nrdata, nprtcl_thispack);
  Kokkos::realloc(prtcl_idata, nidata, nprtcl_thispack);

  HostArray2D<Real> hr("rst_prtcl_rdata_in", nrdata, nprtcl_thispack);
  HostArray2D<int> hi("rst_prtcl_idata_in", nidata, nprtcl_thispack);
  if (nprtcl_thispack > 0) {
    if (resfile.Read_Reals(hr.data(), nrdata*nprtcl_thispack, single_file_per_rank) !=
        static_cast<std::size_t>(nrdata*nprtcl_thispack)) {
      FatalParticleInput("particle restart real data is truncated");
    }
    if (resfile.Read_bytes(hi.data(), sizeof(int), nidata*nprtcl_thispack,
                           single_file_per_rank) !=
        static_cast<std::size_t>(nidata*nprtcl_thispack)) {
      FatalParticleInput("particle restart integer data is truncated");
    }
    Kokkos::deep_copy(prtcl_rdata, hr);
    Kokkos::deep_copy(prtcl_idata, hi);
  }
  pmy_pack->pmesh->UpdateParticleCounts();
}

} // namespace particles
