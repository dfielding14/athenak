//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file pic_migration_restart_ledger.cpp
//! \brief Deterministic mixed particle migration/destruction restart regression.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "particles/particles.hpp"
#include "pgen/pgen.hpp"

namespace {

using particles::Particles;

constexpr int kLedgerSchema = 1;
constexpr int kExpectedCycles = 24;
constexpr int kSplitCycle = 12;
constexpr int kInitialParticles = 192;
constexpr int kExpectedEscapes = 64;

ParameterInput *ledger_pin = nullptr;
std::uint64_t particles_before_cycle = 0;
int committed_cycles = 0;
Real committed_time = 0.0;
int escaped_count = 0;
int boundary_errors = 0;
std::array<Real, Particles::NPIC_BOUNDARY_CONSERVATION> escaped = {};
std::array<Real, Particles::NPIC_BOUNDARY_CONSERVATION> reflected = {};

[[noreturn]] void MigrationFatal(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << message << std::endl;
  std::exit(EXIT_FAILURE);
}

void RequireClose(const std::string &label, const Real measured, const Real expected) {
  const Real scale = std::max(static_cast<Real>(1.0), std::abs(expected));
  if (!std::isfinite(measured) ||
      std::abs(measured - expected) > static_cast<Real>(1.0e-13)*scale) {
    MigrationFatal(label + " does not match the migration regression contract");
  }
}

void StoreLedger() {
  if (ledger_pin == nullptr) MigrationFatal("migration ledger lost ParameterInput");
  ledger_pin->SetInteger("problem", "mig1_ledger_schema", kLedgerSchema);
  ledger_pin->SetBoolean("problem", "mig1_ledger_complete", true);
  ledger_pin->SetInteger("problem", "mig1_committed_cycles", committed_cycles);
  ledger_pin->SetReal("problem", "mig1_committed_time", committed_time);
  ledger_pin->SetInteger("problem", "mig1_escaped_count", escaped_count);
  ledger_pin->SetInteger("problem", "mig1_boundary_errors", boundary_errors);
  constexpr std::array<const char *, Particles::NPIC_BOUNDARY_CONSERVATION> names = {
      "mass", "mom1", "mom2", "mom3", "energy"};
  for (int n = 0; n < Particles::NPIC_BOUNDARY_CONSERVATION; ++n) {
    ledger_pin->SetReal("problem", std::string("mig1_escape_") + names[n],
                        escaped[n]);
    ledger_pin->SetReal("problem", std::string("mig1_reflect_") + names[n],
                        reflected[n]);
  }
}

void LoadLedger(ParameterInput *pin, Mesh *pm) {
  constexpr std::array<const char *, Particles::NPIC_BOUNDARY_CONSERVATION> names = {
      "mass", "mom1", "mom2", "mom3", "energy"};
  if (!pin->DoesParameterExist("problem", "mig1_ledger_schema") ||
      !pin->DoesParameterExist("problem", "mig1_ledger_complete") ||
      pin->GetInteger("problem", "mig1_ledger_schema") != kLedgerSchema ||
      !pin->GetBoolean("problem", "mig1_ledger_complete")) {
    MigrationFatal("migration restart ledger metadata is missing or incomplete");
  }
  committed_cycles = pin->GetInteger("problem", "mig1_committed_cycles");
  committed_time = pin->GetReal("problem", "mig1_committed_time");
  escaped_count = pin->GetInteger("problem", "mig1_escaped_count");
  boundary_errors = pin->GetInteger("problem", "mig1_boundary_errors");
  for (int n = 0; n < Particles::NPIC_BOUNDARY_CONSERVATION; ++n) {
    escaped[n] = pin->GetReal(
        "problem", std::string("mig1_escape_") + names[n]);
    reflected[n] = pin->GetReal(
        "problem", std::string("mig1_reflect_") + names[n]);
  }
  if (committed_cycles != pm->ncycle || committed_time != pm->time ||
      committed_cycles < 0 || escaped_count < 0 || boundary_errors != 0) {
    MigrationFatal("migration restart ledger chronology or values are invalid");
  }
}

std::uint64_t GlobalParticleCount(Mesh *pm) {
  auto *ppart = pm->pmb_pack->ppart;
  std::uint64_t count = static_cast<std::uint64_t>(ppart->nprtcl_thispack);
#if MPI_PARALLEL_ENABLED
  std::uint64_t global_count = 0;
  MPI_Allreduce(&count, &global_count, 1, MPI_UINT64_T, MPI_SUM,
                MPI_COMM_WORLD);
  count = global_count;
#endif
  return count;
}

void MigrationWorkBeforeLoop(Mesh *pm) {
  if (pm == nullptr || pm->pmb_pack == nullptr || pm->pmb_pack->ppart == nullptr) {
    MigrationFatal("migration pre-cycle callback requires particles");
  }
  if (committed_cycles != pm->ncycle || committed_time != pm->time) {
    MigrationFatal("migration pre-cycle ledger chronology is discontinuous");
  }
  pm->pmb_pack->ppart->ResetPICBoundaryConservationDeltas();
  particles_before_cycle = GlobalParticleCount(pm);
}

void MigrationWorkInLoop(Mesh *pm) {
  auto *ppart = pm->pmb_pack->ppart;
  const std::uint64_t particles_after = GlobalParticleCount(pm);
  if (particles_after > particles_before_cycle) {
    MigrationFatal("migration regression unexpectedly created particles");
  }
  const std::uint64_t destroyed = particles_before_cycle - particles_after;
  if (destroyed > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
    MigrationFatal("migration regression destruction count overflowed");
  }

  auto h_escape = Kokkos::create_mirror_view_and_copy(
      HostMemSpace(), ppart->pic_escape_boundary_delta);
  auto h_reflect = Kokkos::create_mirror_view_and_copy(
      HostMemSpace(), ppart->pic_reflecting_boundary_delta);
  auto h_errors = Kokkos::create_mirror_view_and_copy(
      HostMemSpace(), ppart->pic_boundary_conservation_errors);
  std::array<Real, Particles::NPIC_BOUNDARY_CONSERVATION> cycle_escape = {};
  std::array<Real, Particles::NPIC_BOUNDARY_CONSERVATION> cycle_reflect = {};
  for (int n = 0; n < Particles::NPIC_BOUNDARY_CONSERVATION; ++n) {
    cycle_escape[n] = h_escape(n);
    cycle_reflect[n] = h_reflect(n);
  }
  int cycle_errors = h_errors(0);
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, cycle_escape.data(),
                Particles::NPIC_BOUNDARY_CONSERVATION, MPI_ATHENA_REAL, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, cycle_reflect.data(),
                Particles::NPIC_BOUNDARY_CONSERVATION, MPI_ATHENA_REAL, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &cycle_errors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (cycle_errors != 0) {
    MigrationFatal("migration boundary ledger recorded invalid particle state");
  }

  escaped_count += static_cast<int>(destroyed);
  boundary_errors += cycle_errors;
  for (int n = 0; n < Particles::NPIC_BOUNDARY_CONSERVATION; ++n) {
    escaped[n] -= cycle_escape[n];
    reflected[n] += cycle_reflect[n];
  }
  committed_cycles = pm->ncycle + 1;
  committed_time = pm->time + pm->dt;
  StoreLedger();

  if (global_variable::my_rank == 0) {
    const auto old_flags = std::cout.flags();
    const auto old_precision = std::cout.precision();
    std::cout << std::scientific
              << std::setprecision(std::numeric_limits<Real>::max_digits10)
              << "pic_migration_ledger_cycle: cycle=" << committed_cycles
              << " time=" << committed_time
              << " escaped_count=" << destroyed
              << " escape_mass=" << -cycle_escape[Particles::IPIC_BND_MASS]
              << " escape_mom1=" << -cycle_escape[Particles::IPIC_BND_MOM1]
              << " escape_mom2=" << -cycle_escape[Particles::IPIC_BND_MOM2]
              << " escape_mom3=" << -cycle_escape[Particles::IPIC_BND_MOM3]
              << " escape_energy=" << -cycle_escape[Particles::IPIC_BND_ENERGY]
              << " boundary_errors=" << cycle_errors << std::endl;
    std::cout.flags(old_flags);
    std::cout.precision(old_precision);
  }
}

void ValidateCommittedLedger(ParameterInput *pin, Mesh *pm) {
  (void)pin;
  if (committed_cycles != pm->ncycle || committed_time != pm->time ||
      boundary_errors != 0) {
    MigrationFatal("migration checkpoint ledger is not fully committed");
  }
  StoreLedger();
}

void ValidateContract(ParameterInput *pin, Mesh *pm) {
  auto *pmbp = pm->pmb_pack;
  auto *ppart = pmbp->ppart;
  if (global_variable::nranks < 1 || global_variable::nranks > 2 ||
      pmbp->pmhd == nullptr || ppart == nullptr || pm->multilevel || !pm->two_d ||
      pm->mesh_indcs.nx1 != 16 || pm->mesh_indcs.nx2 != 4 ||
      pm->mb_indcs.nx1 != 4 || pm->mb_indcs.nx2 != 4) {
    MigrationFatal(
        "migration regression requires its uniform 16x4, one/two-rank carrier");
  }
  if (pm->mesh_bcs[BoundaryFace::inner_x1] != BoundaryFlag::periodic ||
      pm->mesh_bcs[BoundaryFace::outer_x1] != BoundaryFlag::periodic ||
      pm->mesh_bcs[BoundaryFace::inner_x2] != BoundaryFlag::outflow ||
      pm->mesh_bcs[BoundaryFace::outer_x2] != BoundaryFlag::outflow) {
    MigrationFatal("migration regression boundary contract drifted");
  }
  if (!pmbp->pmhd->peos->eos_data.is_ideal ||
      !ppart->UsesVL2TSCCoupling() || !ppart->pic_boundary_conservation_ledger ||
      !ppart->track_displacement || ppart->nspecies != 3 ||
      pin->GetInteger("time", "nlim") > kExpectedCycles ||
      !pin->GetBoolean("problem", "user_work_in_loop")) {
    MigrationFatal("migration regression MHD-PIC contract drifted");
  }
  RequireClose("<particles>/ppc", pin->GetReal("particles", "ppc"), 3.0);
  RequireClose("<particles>/deposit_qscale",
               pin->GetReal("particles", "deposit_qscale"), 1.0/1024.0);
  RequireClose("<species0>/charge", pin->GetReal("species0", "charge"), 0.0);
  RequireClose("<species1>/charge", pin->GetReal("species1", "charge"), 0.0);
  RequireClose("<species2>/charge", pin->GetReal("species2", "charge"), 0.0);
  RequireClose("<species0>/vx0", pin->GetReal("species0", "vx0"), 4.0);
  RequireClose("<species1>/vx0", pin->GetReal("species1", "vx0"), -4.0);
  RequireClose("<species2>/vy0", pin->GetReal("species2", "vy0"), 1.6);
}

void InitializeUniformMHD(MeshBlockPack *pmbp) {
  auto &indcs = pmbp->pmesh->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  auto &w0 = pmbp->pmhd->w0;
  auto &u0 = pmbp->pmhd->u0;
  auto &b0 = pmbp->pmhd->b0;
  auto &bcc0 = pmbp->pmhd->bcc0;
  par_for("pgen_pic_migration_uniform_mhd", DevExeSpace(),
          0, pmbp->nmb_thispack - 1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    b0.x1f(m, k, j, i) = 0.0;
    b0.x2f(m, k, j, i) = 0.0;
    b0.x3f(m, k, j, i) = 0.0;
    if (i == ie) b0.x1f(m, k, j, i + 1) = 0.0;
    if (j == je) b0.x2f(m, k, j + 1, i) = 0.0;
    if (k == ke) b0.x3f(m, k + 1, j, i) = 0.0;
    w0(m, IDN, k, j, i) = 1.0;
    w0(m, IVX, k, j, i) = 0.0;
    w0(m, IVY, k, j, i) = 0.0;
    w0(m, IVZ, k, j, i) = 0.0;
    w0(m, IPR, k, j, i) = 1.0;
    bcc0(m, IBX, k, j, i) = 0.0;
    bcc0(m, IBY, k, j, i) = 0.0;
    bcc0(m, IBZ, k, j, i) = 0.0;
  });
  pmbp->pmhd->peos->PrimToCons(w0, bcc0, u0, is, ie, js, je, ks, ke);
}

}  // namespace

void ProblemGenerator::PICMigrationRestartLedger(ParameterInput *pin,
                                                  const bool restart) {
  ledger_pin = pin;
  ValidateContract(pin, pmy_mesh_);
  user_work_in_loop = true;
  user_work_before_loop_func = MigrationWorkBeforeLoop;
  user_work_in_loop_func = MigrationWorkInLoop;
  pgen_checkpoint_func = ValidateCommittedLedger;
  pgen_final_func = ValidateCommittedLedger;

  if (restart) {
    LoadLedger(pin, pmy_mesh_);
    return;
  }

  committed_cycles = 0;
  committed_time = pmy_mesh_->time;
  escaped_count = 0;
  boundary_errors = 0;
  escaped.fill(0.0);
  reflected.fill(0.0);
  StoreLedger();
  InitializeUniformMHD(pmy_mesh_->pmb_pack);

  if (GlobalParticleCount(pmy_mesh_) != kInitialParticles) {
    MigrationFatal("migration regression initial particle count drifted");
  }
  if (kSplitCycle >= kExpectedCycles || kExpectedEscapes >= kInitialParticles) {
    MigrationFatal("migration regression compile-time contract is invalid");
  }
}
