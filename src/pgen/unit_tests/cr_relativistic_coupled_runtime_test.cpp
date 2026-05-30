//========================================================================================
// AthenaK astrophysical fluid dynamics & numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cr_relativistic_coupled_runtime_test.cpp
//! \brief Serial witness for experimental frozen-t^n ideal-MHD relativistic CR coupling.

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "parameter_input.hpp"
#include "particles/particles.hpp"
#include "particles/relativistic_state.hpp"
#include "particles/trilinear_gather.hpp"
#include "pgen/pgen.hpp"

namespace {

enum LiveVelocityProfile {
  kUniformVelocity = 0,
  kAffineVelocity = 1,
  kPeriodicVelocity = 2
};

constexpr Real kTwoPi = 6.283185307179586476925286766559;
Real manufactured_acceleration = 0.0;

[[noreturn]] void Fail(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl
            << message << std::endl;
  std::exit(EXIT_FAILURE);
}

KOKKOS_INLINE_FUNCTION
particles::ParticleVector3 LiveVelocity(
    int profile, Real x1, Real x2, Real x3, Real u1, Real u2, Real u3,
    Real u1_x1, Real u1_x2, Real u1_x3, Real u2_x1, Real u2_x2, Real u2_x3,
    Real u3_x1, Real u3_x2, Real u3_x3, Real amplitude, Real wave_number) {
  particles::ParticleVector3 velocity{u1, u2, u3};
  if (profile == kAffineVelocity) {
    velocity.x += u1_x1*x1 + u1_x2*x2 + u1_x3*x3;
    velocity.y += u2_x1*x1 + u2_x2*x2 + u2_x3*x3;
    velocity.z += u3_x1*x1 + u3_x2*x2 + u3_x3*x3;
  } else if (profile == kPeriodicVelocity) {
    const Real k = kTwoPi*wave_number;
    velocity.x += amplitude*Kokkos::sin(k*x1)*Kokkos::cos(k*x2);
    velocity.y += amplitude*Kokkos::sin(k*x2)*Kokkos::cos(k*x3);
    velocity.z += amplitude*Kokkos::sin(k*x3)*Kokkos::cos(k*x1);
  }
  return velocity;
}

int ReadLiveVelocityProfile(ParameterInput *pin) {
  const std::string profile =
      pin->GetOrAddString("problem", "live_velocity_profile", "uniform");
  if (profile.compare("uniform") == 0) {return kUniformVelocity;}
  if (profile.compare("affine") == 0) {return kAffineVelocity;}
  if (profile.compare("periodic") == 0) {return kPeriodicVelocity;}
  Fail("<problem>/live_velocity_profile must be 'uniform', 'affine', or 'periodic'");
}

bool IsCoupledMHDIdeal(MeshBlockPack *pmbp) {
  return pmbp->ppart != nullptr &&
         pmbp->ppart->pusher == ParticlesPusher::relativistic_hc &&
         pmbp->ppart->relativistic_field_source == RelativisticFieldSource::mhd_ideal;
}

void InitializeLiveMHDFields(ParameterInput *pin, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp == nullptr || pmbp->pmhd == nullptr) {
    Fail("CR relativistic coupled runtime witness requires MHD");
  }

  const int profile = ReadLiveVelocityProfile(pin);
  const Real density = pin->GetOrAddReal("problem", "live_density", 1.0);
  const Real pressure = pin->GetOrAddReal("problem", "live_pressure", 0.6);
  const Real u1 = pin->GetOrAddReal("problem", "live_u1", 0.0);
  const Real u2 = pin->GetOrAddReal("problem", "live_u2", 0.0);
  const Real u3 = pin->GetOrAddReal("problem", "live_u3", 0.0);
  const Real u1_x1 = pin->GetOrAddReal("problem", "live_u1_x1", 0.0);
  const Real u1_x2 = pin->GetOrAddReal("problem", "live_u1_x2", 0.0);
  const Real u1_x3 = pin->GetOrAddReal("problem", "live_u1_x3", 0.0);
  const Real u2_x1 = pin->GetOrAddReal("problem", "live_u2_x1", 0.0);
  const Real u2_x2 = pin->GetOrAddReal("problem", "live_u2_x2", 0.0);
  const Real u2_x3 = pin->GetOrAddReal("problem", "live_u2_x3", 0.0);
  const Real u3_x1 = pin->GetOrAddReal("problem", "live_u3_x1", 0.0);
  const Real u3_x2 = pin->GetOrAddReal("problem", "live_u3_x2", 0.0);
  const Real u3_x3 = pin->GetOrAddReal("problem", "live_u3_x3", 0.0);
  const Real amplitude = pin->GetOrAddReal("problem", "live_velocity_amplitude", 0.0);
  const Real wave_number = pin->GetOrAddReal("problem", "live_velocity_wave_number", 1.0);
  const Real b1 = pin->GetOrAddReal("problem", "live_B1", 0.0);
  const Real b2 = pin->GetOrAddReal("problem", "live_B2", 0.0);
  const Real b3 = pin->GetOrAddReal("problem", "live_B3", 0.4);

  auto &w0 = pmbp->pmhd->w0;
  auto &u0 = pmbp->pmhd->u0;
  auto &b0 = pmbp->pmhd->b0;
  auto &bcc0 = pmbp->pmhd->bcc0;
  auto &size = pmbp->pmb->mb_size;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int ng = indcs.ng;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  par_for("cr_rel_coupled_live_mhd", DevExeSpace(), 0, pmbp->nmb_thispack - 1,
          ks - ng, ke + ng, js - ng, je + ng, is - ng, ie + ng,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real x1 = CellCenterX(i - is, nx1, size.d_view(m).x1min,
                                size.d_view(m).x1max);
    const Real x2 = CellCenterX(j - js, nx2, size.d_view(m).x2min,
                                size.d_view(m).x2max);
    const Real x3 = CellCenterX(k - ks, nx3, size.d_view(m).x3min,
                                size.d_view(m).x3max);
    const auto velocity = LiveVelocity(
        profile, x1, x2, x3, u1, u2, u3, u1_x1, u1_x2, u1_x3, u2_x1,
        u2_x2, u2_x3, u3_x1, u3_x2, u3_x3, amplitude, wave_number);
    w0(m,IDN,k,j,i) = density;
    w0(m,IVX,k,j,i) = velocity.x;
    w0(m,IVY,k,j,i) = velocity.y;
    w0(m,IVZ,k,j,i) = velocity.z;
    w0(m,IEN,k,j,i) = pressure;
    bcc0(m,IBX,k,j,i) = b1;
    bcc0(m,IBY,k,j,i) = b2;
    bcc0(m,IBZ,k,j,i) = b3;
    b0.x1f(m,k,j,i) = b1;
    b0.x2f(m,k,j,i) = b2;
    b0.x3f(m,k,j,i) = b3;
    if (i == ie + ng) {b0.x1f(m,k,j,i+1) = b1;}
    if (j == je + ng) {b0.x2f(m,k,j+1,i) = b2;}
    if (k == ke + ng) {b0.x3f(m,k+1,j,i) = b3;}
  });
  pmbp->pmhd->peos->PrimToCons(w0, bcc0, u0, is, ie, js, je, ks, ke);
}

void ApplyManufacturedAcceleration(Mesh *pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &u0 = pmbp->pmhd->u0;
  auto &indcs = pm->mb_indcs;
  const Real acceleration = manufactured_acceleration;
  par_for("cr_rel_coupled_manufactured_acceleration", DevExeSpace(),
          0, pmbp->nmb_thispack - 1, indcs.ks, indcs.ke, indcs.js, indcs.je,
          indcs.is, indcs.ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real density = u0(m,IDN,k,j,i);
    const Real momentum_old = u0(m,IM2,k,j,i);
    const Real momentum_new = momentum_old + density*acceleration*bdt;
    u0(m,IEN,k,j,i) +=
        0.5*(momentum_new*momentum_new - momentum_old*momentum_old)/density;
    u0(m,IM2,k,j,i) = momentum_new;
  });
}

void WritePrepushWitness(ParameterInput *pin, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto pr_h = Kokkos::create_mirror_view_and_copy(
      HostMemSpace(), pmbp->ppart->prtcl_rdata);
  auto pi_h = Kokkos::create_mirror_view_and_copy(
      HostMemSpace(), pmbp->ppart->prtcl_idata);
  const std::string filename = pin->GetOrAddString(
      "problem", "relativistic_coupled_prepush_witness",
      "cr_relativistic_coupled_prepush.tsv");
  std::ofstream stream(filename);
  if (!stream.is_open()) {Fail("could not open CR relativistic coupled prepush output");}
  stream << "# gid\tb1\tb2\tb3\tcE1\tcE2\tcE3\n";
  stream << pi_h(PGID,0) << '\t' << pr_h(IPBX,0) << '\t' << pr_h(IPBY,0) << '\t'
         << pr_h(IPBZ,0) << '\t' << pr_h(IPCEX,0) << '\t' << pr_h(IPCEY,0) << '\t'
         << pr_h(IPCEZ,0) << '\n';
}

void PoisonStoredParticleFields(ParameterInput *pin, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &pr = pmbp->ppart->prtcl_rdata;
  auto &pi = pmbp->ppart->prtcl_idata;
  const int npart = pmbp->ppart->nprtcl_thispack;
  const int invalid_gid = pmbp->gids + pmbp->nmb_thispack;
  const bool poison_gid = pin->GetOrAddBoolean("problem", "runtime_invalid_gid", false);
  const bool poison_position =
      pin->GetOrAddBoolean("problem", "runtime_invalid_stencil_position", false);
  const Real poisoned_x1 =
      pin->GetOrAddReal("problem", "runtime_invalid_stencil_x1", 10.0);
  const Real b1 = pin->GetOrAddReal("problem", "poison_B1", 7.0);
  const Real b2 = pin->GetOrAddReal("problem", "poison_B2", -11.0);
  const Real b3 = pin->GetOrAddReal("problem", "poison_B3", 13.0);
  const Real ce1 = pin->GetOrAddReal("problem", "poison_cE1", -17.0);
  const Real ce2 = pin->GetOrAddReal("problem", "poison_cE2", 19.0);
  const Real ce3 = pin->GetOrAddReal("problem", "poison_cE3", -23.0);
  par_for("cr_rel_coupled_poison_stored_fields", DevExeSpace(), 0, npart - 1,
  KOKKOS_LAMBDA(const int p) {
    pr(IPBX,p) = b1;
    pr(IPBY,p) = b2;
    pr(IPBZ,p) = b3;
    pr(IPCEX,p) = ce1;
    pr(IPCEY,p) = ce2;
    pr(IPCEZ,p) = ce3;
    if (poison_gid) {pi(PGID,p) = invalid_gid;}
    if (poison_position) {pr(IPX,p) = poisoned_x1;}
  });
}

void WriteCoupledRuntimeWitness(ParameterInput *pin, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (global_variable::nranks != 1) {
    Fail("CR relativistic coupled runtime witness is intentionally serial-only");
  }
  if (pmbp == nullptr || pmbp->pmhd == nullptr || pmbp->ppart == nullptr) {
    Fail("CR relativistic coupled runtime witness requires MHD and particles");
  }
  if (pmbp->ppart->nprtcl_thispack != 1) {
    Fail("CR relativistic coupled runtime witness requires exactly one particle");
  }

  auto &pr = pmbp->ppart->prtcl_rdata;
  auto &pi = pmbp->ppart->prtcl_idata;
  auto &w0 = pmbp->pmhd->w0;
  auto &bcc0 = pmbp->pmhd->bcc0;
  auto &mbsize = pmbp->pmb->mb_size;
  auto &indcs = pm->mb_indcs;
  const int gids = pmbp->gids;
  const bool multi_d = pm->multi_d;
  const bool three_d = pm->three_d;
  DvceArray1D<Real> live("cr_rel_coupled_final_live_sample", 10);
  par_for("cr_rel_coupled_final_live_gather", DevExeSpace(), 0, 0,
  KOKKOS_LAMBDA(const int p) {
    const int m = pi(PGID,p) - gids;
    const auto fields = particles::GatherNewtonianCellCenteredFields(
        w0, bcc0, mbsize, m, pr(IPX,p), pr(IPY,p), pr(IPZ,p), indcs.is,
        indcs.js, indcs.ks, multi_d, three_d);
    const auto cE = particles::ConstructIdealCE(fields.u, fields.b);
    live(0) = fields.u.x;
    live(1) = fields.u.y;
    live(2) = fields.u.z;
    live(3) = fields.b.x;
    live(4) = fields.b.y;
    live(5) = fields.b.z;
    live(6) = cE.x;
    live(7) = cE.y;
    live(8) = cE.z;
    live(9) = cE.x*fields.b.x + cE.y*fields.b.y + cE.z*fields.b.z;
  });

  const auto pr_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pr);
  const auto pi_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pi);
  const auto live_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), live);
  const particles::relativistic::Vector3 w{
      pr_h(IPWX,0), pr_h(IPWY,0), pr_h(IPWZ,0)};
  Real gamma = 0.0;
  Real kinetic_energy = 0.0;
  if (particles::relativistic::GammaFromW(w, pmbp->ppart->c_model, gamma) !=
          particles::relativistic::StateStatus::success ||
      particles::relativistic::KineticEnergyFromW(
          w, pmbp->ppart->c_model, kinetic_energy) !=
          particles::relativistic::StateStatus::success) {
    Fail("CR relativistic coupled runtime witness found invalid authoritative w");
  }

  const std::string filename = pin->GetOrAddString(
      "problem", "relativistic_coupled_runtime_witness",
      "cr_relativistic_coupled_runtime.tsv");
  std::ofstream stream(filename);
  if (!stream.is_open()) {
    Fail("could not open CR relativistic coupled runtime witness output");
  }
  stream << "# coupled-runtime\tmhd_ideal\tfrozen_tn\tcE=-u_cross_B\n";
  stream << "# cycle\ttime\tgid\ttag\tx\ty\tz\tvx\tvy\tvz\twx\twy\twz"
         << "\tcE1\tcE2\tcE3\tb1\tb2\tb3\twork\talpha\tgamma\tkinetic"
         << "\tdx\tdy\tdz\tdb\tlegacy_ipm\tfinal_live_u1\tfinal_live_u2"
         << "\tfinal_live_u3\tfinal_live_b1\tfinal_live_b2\tfinal_live_b3"
         << "\tfinal_live_cE1\tfinal_live_cE2\tfinal_live_cE3"
         << "\tfinal_live_cE_dot_b\tparticle_dtnew\tlast_subcycle_steps"
         << "\tlast_push_dt\n";
  stream << std::setprecision(std::numeric_limits<Real>::max_digits10)
         << pm->ncycle << '\t' << pm->time << '\t' << pi_h(PGID,0) << '\t'
         << pi_h(PTAG,0) << '\t'
         << pr_h(IPX,0) << '\t' << pr_h(IPY,0) << '\t' << pr_h(IPZ,0) << '\t'
         << pr_h(IPVX,0) << '\t' << pr_h(IPVY,0) << '\t' << pr_h(IPVZ,0) << '\t'
         << pr_h(IPWX,0) << '\t' << pr_h(IPWY,0) << '\t' << pr_h(IPWZ,0) << '\t'
         << pr_h(IPCEX,0) << '\t' << pr_h(IPCEY,0) << '\t' << pr_h(IPCEZ,0) << '\t'
         << pr_h(IPBX,0) << '\t' << pr_h(IPBY,0) << '\t' << pr_h(IPBZ,0) << '\t'
         << pr_h(IPWORK,0) << '\t' << pr_h(IPALPHA,0) << '\t' << gamma << '\t'
         << kinetic_energy << '\t' << pr_h(IPDX,0) << '\t' << pr_h(IPDY,0) << '\t'
         << pr_h(IPDZ,0) << '\t' << pr_h(IPDB,0) << '\t' << pr_h(IPM,0);
  for (int n=0; n<10; ++n) {stream << '\t' << live_h(n);}
  stream << '\t' << pmbp->ppart->dtnew << '\t'
         << pmbp->ppart->last_subcycle_steps << '\t'
         << pmbp->ppart->last_push_dt << '\n';
}

} // namespace

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  PartRandom(pin, restart);
  pgen_final_func = &WriteCoupledRuntimeWitness;
  if (restart || !IsCoupledMHDIdeal(pmy_mesh_->pmb_pack)) {return;}
  InitializeLiveMHDFields(pin, pmy_mesh_);
  pmy_mesh_->pmb_pack->ppart->MarkRelativisticTimestepBoundDirty();
  const Real dtnew_override =
      pin->GetOrAddReal("problem", "runtime_particle_dtnew_override", -1.0);
  if (dtnew_override > 0.0) {
    pmy_mesh_->pmb_pack->ppart->dtnew = dtnew_override;
    pmy_mesh_->pmb_pack->ppart->relativistic_timestep_bound_dirty = false;
  }
  PoisonStoredParticleFields(pin, pmy_mesh_);
  WritePrepushWitness(pin, pmy_mesh_);
  if (pin->GetOrAddBoolean("problem", "relativistic_coupled_manufactured_test", false)) {
    manufactured_acceleration =
        pin->GetOrAddReal("problem", "manufactured_acceleration", 0.2);
    user_srcs_func = &ApplyManufacturedAcceleration;
  }
}
