//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_landau_fluid.cpp
//! \brief CGL Landau-fluid heat-flux closure and parabolic timestep bound.

#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "mesh/nghbr_index.hpp"
#include "eos/eos.hpp"
#include "eos/cgl_physics.hpp"
#include "diffusion/cgl_landau_fluid.hpp"
#include "diffusion/cgl_landau_fluid_arithmetic.hpp"

namespace {

bool CGLProfileEnvValue(const char *name, bool fallback) {
  const char *value = std::getenv(name);
  if (value == nullptr) {
    return fallback;
  }
  const std::string text(value);
  if (text == "1" || text == "true" || text == "TRUE" ||
      text == "yes" || text == "YES" || text == "on" || text == "ON") {
    return true;
  }
  if (text == "0" || text == "false" || text == "FALSE" ||
      text == "no" || text == "NO" || text == "off" || text == "OFF") {
    return false;
  }
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
            << name << " = '" << text
            << "' is not a boolean value." << std::endl;
  std::exit(EXIT_FAILURE);
}

const char *CGLLFProfileBucketName(CGLLFProfileBucket bucket) {
  switch (bucket) {
    case CGLLFProfileBucket::heat_flux_total:
      return "heat_flux_total";
    case CGLLFProfileBucket::heat_flux_precompute:
      return "heat_flux_precompute";
    case CGLLFProfileBucket::heat_flux_flux1:
      return "heat_flux_flux1";
    case CGLLFProfileBucket::heat_flux_flux1_gradients:
      return "heat_flux_flux1_gradients";
    case CGLLFProfileBucket::heat_flux_flux1_face_state:
      return "heat_flux_flux1_face_state";
    case CGLLFProfileBucket::heat_flux_flux1_closure:
      return "heat_flux_flux1_closure";
    case CGLLFProfileBucket::heat_flux_flux2:
      return "heat_flux_flux2";
    case CGLLFProfileBucket::heat_flux_flux2_gradients:
      return "heat_flux_flux2_gradients";
    case CGLLFProfileBucket::heat_flux_flux2_face_state:
      return "heat_flux_flux2_face_state";
    case CGLLFProfileBucket::heat_flux_flux2_closure:
      return "heat_flux_flux2_closure";
    case CGLLFProfileBucket::heat_flux_flux3:
      return "heat_flux_flux3";
    case CGLLFProfileBucket::heat_flux_flux3_gradients:
      return "heat_flux_flux3_gradients";
    case CGLLFProfileBucket::heat_flux_flux3_face_state:
      return "heat_flux_flux3_face_state";
    case CGLLFProfileBucket::heat_flux_flux3_closure:
      return "heat_flux_flux3_closure";
    case CGLLFProfileBucket::heat_flux_work_diagnostics:
      return "heat_flux_work_diagnostics";
    case CGLLFProfileBucket::timestep_reduction:
      return "timestep_reduction";
    case CGLLFProfileBucket::sweep_begin_conversion:
      return "sweep_begin_conversion";
    case CGLLFProfileBucket::sts_clear_flux:
      return "sts_clear_flux";
    case CGLLFProfileBucket::sts_update_copies:
      return "sts_update_copies";
    case CGLLFProfileBucket::sts_update_kernel:
      return "sts_update_kernel";
    case CGLLFProfileBucket::primitive_refresh:
      return "primitive_refresh";
    case CGLLFProfileBucket::admissibility:
      return "admissibility";
    case CGLLFProfileBucket::sweep_end_conversion:
      return "sweep_end_conversion";
    case CGLLFProfileBucket::post_sweep_collisions:
      return "post_sweep_collisions";
    case CGLLFProfileBucket::parabolic_init_recv:
      return "parabolic_init_recv";
    case CGLLFProfileBucket::parabolic_send_flux:
      return "parabolic_send_flux";
    case CGLLFProfileBucket::parabolic_recv_flux:
      return "parabolic_recv_flux";
    case CGLLFProfileBucket::parabolic_restrict_u:
      return "parabolic_restrict_u";
    case CGLLFProfileBucket::parabolic_send_u:
      return "parabolic_send_u";
    case CGLLFProfileBucket::parabolic_recv_u:
      return "parabolic_recv_u";
    case CGLLFProfileBucket::parabolic_physical_bcs:
      return "parabolic_physical_bcs";
    case CGLLFProfileBucket::parabolic_prolongate:
      return "parabolic_prolongate";
    case CGLLFProfileBucket::count:
      break;
  }
  return "unknown";
}

parabolic::ParabolicIntegratorMode ParseCGLHeatFluxIntegrator(ParameterInput *pin) {
  const std::string integrator =
      pin->GetOrAddString("mhd", "cgl_heat_flux_integrator", "sts");
  if (integrator == "sts") {
    return parabolic::ParabolicIntegratorMode::sts;
  }
  if (integrator == "explicit") {
    return parabolic::ParabolicIntegratorMode::explicit_mode;
  }
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
            << "<mhd>/cgl_heat_flux_integrator = '" << integrator
            << "' is not implemented; valid choices are [sts,explicit]."
            << std::endl;
  std::exit(EXIT_FAILURE);
}

struct CGLLFFaceState {
  Real rho;
  Real ppar;
  Real pperp;
  Real bmag_inv;
  Real bhx;
  Real bhy;
  Real bhz;
  Real bhdir;
  Real cparallel;
  Real lf_k;
  Real nu;
};

KOKKOS_INLINE_FUNCTION
Real ScaledMagneticMagnitude(const Real bx, const Real by, const Real bz) {
  const Real scale = fmax(fabs(bx), fmax(fabs(by), fabs(bz)));
  if (scale == 0.0) {
    return 0.0;
  }
  const Real sx = bx/scale;
  const Real sy = by/scale;
  const Real sz = bz/scale;
  const Real scaled_magnitude = sqrt(sx*sx + sy*sy + sz*sz);
  const Real maximum = std::numeric_limits<Real>::max();
  return (scaled_magnitude > maximum/scale)
             ? maximum
             : scale*scaled_magnitude;
}

KOKKOS_INLINE_FUNCTION
bool BuildCGLLFFaceState(const Real rho_l, const Real rho_r,
                         const Real ppar_l, const Real ppar_r,
                         const Real pperp_l, const Real pperp_r,
                         const Real bx, const Real by, const Real bz, const int dir,
                         const Real lf_k, const bool coeff_local, const Real cparallel0,
                         const bool backup, const EOS_Data &eos,
                         CGLLFFaceState &face) {
  const Real bscale = fmax(fabs(bx), fmax(fabs(by), fabs(bz)));
  if (bscale == 0.0 || lf_k <= 0.0) {
    return false;
  }
  const Real bsx = bx/bscale;
  const Real bsy = by/bscale;
  const Real bsz = bz/bscale;
  const Real bscaled_mag = sqrt(bsx*bsx + bsy*bsy + bsz*bsz);
  const Real maximum = std::numeric_limits<Real>::max();
  const Real bmag = (bscaled_mag > maximum/bscale)
                        ? maximum
                        : bscale*bscaled_mag;
  if (bmag <= eos.bfloor || lf_k <= 0.0) {
    return false;
  }
  face.rho = fmax(static_cast<Real>(0.5)*rho_l +
                  static_cast<Real>(0.5)*rho_r, eos.dfloor);
  face.ppar = fmax(static_cast<Real>(0.5)*ppar_l +
                   static_cast<Real>(0.5)*ppar_r, eos.pfloor);
  face.pperp = fmax(static_cast<Real>(0.5)*pperp_l +
                    static_cast<Real>(0.5)*pperp_r, eos.pfloor);
  face.bmag_inv = (static_cast<Real>(1.0)/bscale)/bscaled_mag;
  face.bhx = bsx/bscaled_mag;
  face.bhy = bsy/bscaled_mag;
  face.bhz = bsz/bscaled_mag;
  face.bhdir = (dir == 0) ? face.bhx : ((dir == 1) ? face.bhy : face.bhz);
  face.cparallel =
      coeff_local ? sqrt(fmax(face.ppar/face.rho, eos.tfloor)) : cparallel0;
  const Real sqrt_max = sqrt(maximum);
  const Real bsqr = (bmag <= sqrt_max) ? bmag*bmag
                                       : maximum;
  const Real nu = fmax(eos.nu_coll, static_cast<Real>(0.0)) +
      cgl::LimiterCollisionRate(face.ppar, face.pperp, bsqr, eos.lim_coll,
                                eos.mlim, eos.flim, eos.firehose_threshold,
                                backup);
  face.lf_k = lf_k;
  face.nu = nu;
  return true;
}

KOKKOS_INLINE_FUNCTION
void CGLLFFlux(const CGLLFFaceState &face, const Real gtpar_x, const Real gtpar_y,
               const Real gtpar_z, const Real gtperp_x, const Real gtperp_y,
               const Real gtperp_z, const Real gb_x, const Real gb_y, const Real gb_z,
               const Real dt_sweep, const Real rkl_weight, Real &eflux,
               Real &muflux, cgl_lf::ScaledValue &weighted_qpar_flux,
               cgl_lf::ScaledValue &weighted_qperp_flux, Real &qpar_ratio,
               Real &qperp_ratio) {
  const Real grad_tpar =
      face.bhx*gtpar_x + face.bhy*gtpar_y + face.bhz*gtpar_z;
  const Real grad_tperp =
      face.bhx*gtperp_x + face.bhy*gtperp_y + face.bhz*gtperp_z;
  const Real grad_b = face.bhx*gb_x + face.bhy*gb_y + face.bhz*gb_z;
  Real signed_qpar_ratio = 0.0;
  signed_qpar_ratio = cgl::ParallelHeatFluxRatio(
      face.cparallel, face.rho, face.ppar, face.lf_k, face.nu, grad_tpar);
  qpar_ratio = fabs(signed_qpar_ratio);
  Real signed_qperp_ratio = 0.0;
  signed_qperp_ratio = cgl::PerpendicularHeatFluxRatio(
      face.cparallel, face.rho, face.ppar, face.pperp, face.bmag_inv,
      face.lf_k, face.nu, grad_tperp, grad_b);
  qperp_ratio = fabs(signed_qperp_ratio);

  weighted_qpar_flux = cgl_lf::LimitedHeatFlux(
      signed_qpar_ratio, face.cparallel, cgl::kSqrtEightOverPi, face.ppar);
  weighted_qpar_flux = cgl_lf::Multiply(weighted_qpar_flux, face.bhdir);
  weighted_qpar_flux = cgl_lf::Multiply(weighted_qpar_flux, dt_sweep);
  weighted_qpar_flux = cgl_lf::Multiply(weighted_qpar_flux, rkl_weight);

  weighted_qperp_flux = cgl_lf::LimitedHeatFlux(
      signed_qperp_ratio, face.cparallel, cgl::kSqrtTwoOverPi, face.pperp);
  weighted_qperp_flux = cgl_lf::Multiply(weighted_qperp_flux, face.bhdir);
  weighted_qperp_flux = cgl_lf::Multiply(weighted_qperp_flux, dt_sweep);
  weighted_qperp_flux = cgl_lf::Multiply(weighted_qperp_flux, rkl_weight);

  // The common stage weight commutes with AMR restriction and remains attached
  // until the final face values are representable.
  eflux = cgl_lf::Materialize(cgl_lf::Add(
      weighted_qperp_flux,
      cgl_lf::Multiply(weighted_qpar_flux, static_cast<Real>(0.5))));
  muflux = cgl_lf::Materialize(
      cgl_lf::Multiply(weighted_qperp_flux, face.bmag_inv));
}

KOKKOS_INLINE_FUNCTION
bool OwnsHeatFluxDiagnosticFace(const int m, const int direction, const int face_index,
                                const int lower, const int upper, const bool multilevel,
                                const int my_level,
                                const DualArray2D<NeighborBlock> &nghbr) {
  if (face_index > lower && face_index <= upper) {
    return true;
  }
  if (face_index != lower && face_index != upper + 1) {
    return false;
  }
  if (!multilevel) {
    return face_index == lower;
  }

  const int side = (face_index == lower) ? -1 : 1;
  for (int n1 = 0; n1 < 2; ++n1) {
    for (int n2 = 0; n2 < 2; ++n2) {
      const int neighbor = (direction == 0) ? NeighborIndex(side, 0, 0, n1, n2)
                         : ((direction == 1) ? NeighborIndex(0, side, 0, n1, n2)
                                             : NeighborIndex(0, 0, side, n1, n2));
      if (nghbr.d_view(m, neighbor).gid < 0) {
        continue;
      }
      const int neighbor_level = nghbr.d_view(m, neighbor).lev;
      if (face_index == lower && neighbor_level > my_level) {
        return false;
      }
      if (face_index == upper + 1 && neighbor_level < my_level) {
        return true;
      }
    }
  }
  return face_index == lower;
}

} // namespace

CGLLandauFluid::CGLLandauFluid(MeshBlockPack *pp, ParameterInput *pin) :
    dtnew(static_cast<Real>(std::numeric_limits<float>::max())),
    lf_k_parallel(0.0),
    lf_coeff_local(true),
    lf_c_parallel0(0.0),
    strict_admissibility(false),
    effective_backup_limiter(false),
    mode(parabolic::ParabolicIntegratorMode::sts),
    pmy_pack(pp),
    tpar_("cgl_lf_tpar", 1, 1, 1, 1),
    tperp_("cgl_lf_tperp", 1, 1, 1, 1),
    bmag_("cgl_lf_bmag", 1, 1, 1, 1) {
  const std::string model = pin->GetString("mhd", "cgl_heat_flux");
  if (model != "landau_fluid") {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<mhd>/cgl_heat_flux = '" << model << "' must be 'landau_fluid'"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  lf_k_parallel = pin->GetReal("mhd", "lf_k_parallel");
  if (lf_k_parallel <= 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<mhd>/lf_k_parallel must be positive." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const std::string coeff_mode =
      pin->GetOrAddString("mhd", "lf_coefficient_mode", "local");
  if (coeff_mode == "background") {
    lf_coeff_local = false;
    lf_c_parallel0 = pin->GetReal("mhd", "lf_c_parallel0");
    if (lf_c_parallel0 <= 0.0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<mhd>/lf_c_parallel0 must be positive." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else if (coeff_mode != "local") {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<mhd>/lf_coefficient_mode = '" << coeff_mode
              << "' must be 'local' or 'background'." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  mode = ParseCGLHeatFluxIntegrator(pin);
  strict_admissibility =
      pin->GetOrAddBoolean("mhd", "cgl_lf_strict_admissibility", false);
  const bool configured_backup =
      pin->GetOrAddBoolean("mhd", "backup_limiters", false);
  const bool instability_limiter_active =
      pin->GetOrAddBoolean("mhd", "mirror_limiter", false) ||
      pin->GetOrAddBoolean("mhd", "firehose_limiter", false);
  effective_backup_limiter = cgl::EffectiveBackupLimiter(
      configured_backup, true, instability_limiter_active,
      strict_admissibility);
  if (effective_backup_limiter && !configured_backup &&
      global_variable::my_rank == 0) {
    std::cout << "CGL Landau-fluid relaxed admissibility enables the emergency "
              << "hard-bound backup limiter; configured backup_limiters remains false."
              << std::endl;
  }
  profile_enabled_ =
      CGLProfileEnvValue("ATHENAK_CGL_LF_PROFILE",
                         pin->GetOrAddBoolean("mhd", "cgl_lf_profile", false));
  profile_detail_enabled_ =
      CGLProfileEnvValue(
          "ATHENAK_CGL_LF_PROFILE_DETAIL",
          pin->GetOrAddBoolean("mhd", "cgl_lf_profile_detail", false));
  profile_detail_enabled_ = profile_enabled_ && profile_detail_enabled_;
  if (profile_enabled_ && global_variable::my_rank == 0) {
    std::cout << "CGL Landau-fluid profiling enabled; timing regions use Kokkos "
              << "fences and are intended for profiling runs only." << std::endl;
    if (profile_detail_enabled_) {
      std::cout << "CGL Landau-fluid detailed profiling enabled; additional "
                << "probe kernels replay directional heat-flux sub-work and "
                << "do not update evolved state." << std::endl;
    }
  }
}

CGLLFProfileRegion::CGLLFProfileRegion(CGLLandauFluid *profile,
                                       CGLLFProfileBucket bucket) :
    profile_(profile),
    bucket_(bucket),
    active_(profile != nullptr && profile->ProfileEnabled()) {
  if (active_) {
    Kokkos::fence();
    timer_.reset();
  }
}

CGLLFProfileRegion::~CGLLFProfileRegion() {
  if (active_) {
    Kokkos::fence();
    profile_->AddProfileTime(bucket_, static_cast<Real>(timer_.seconds()));
  }
}

void CGLLandauFluid::AddProfileTime(CGLLFProfileBucket bucket, Real seconds) {
  if (!profile_enabled_) {
    return;
  }
  const int index = static_cast<int>(bucket);
  if (index < 0 || index >= kCGLLFProfileBucketCount) {
    return;
  }
  profile_seconds_[index] += seconds;
  profile_counts_[index] += 1;
}

void CGLLandauFluid::ReportProfile(const char *context) const {
  if (!profile_enabled_) {
    return;
  }

  Real local_seconds[kCGLLFProfileBucketCount] = {};
  Real local_counts[kCGLLFProfileBucketCount] = {};
  Real sum_seconds[kCGLLFProfileBucketCount] = {};
  Real max_seconds[kCGLLFProfileBucketCount] = {};
  Real sum_counts[kCGLLFProfileBucketCount] = {};
  Real max_counts[kCGLLFProfileBucketCount] = {};
  for (int n = 0; n < kCGLLFProfileBucketCount; ++n) {
    local_seconds[n] = profile_seconds_[n];
    local_counts[n] = static_cast<Real>(profile_counts_[n]);
  }

  int max_nstages = profile_max_nstages_;
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(local_seconds, sum_seconds, kCGLLFProfileBucketCount, MPI_ATHENA_REAL,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(local_seconds, max_seconds, kCGLLFProfileBucketCount, MPI_ATHENA_REAL,
                MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(local_counts, sum_counts, kCGLLFProfileBucketCount, MPI_ATHENA_REAL,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(local_counts, max_counts, kCGLLFProfileBucketCount, MPI_ATHENA_REAL,
                MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&profile_max_nstages_, &max_nstages, 1, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);
#else
  for (int n = 0; n < kCGLLFProfileBucketCount; ++n) {
    sum_seconds[n] = local_seconds[n];
    max_seconds[n] = local_seconds[n];
    sum_counts[n] = local_counts[n];
    max_counts[n] = local_counts[n];
  }
#endif

  if (global_variable::my_rank != 0) {
    return;
  }

  const Mesh *pm = pmy_pack->pmesh;
  const auto &mesh_indcs = pm->mesh_indcs;
  const auto &mb_indcs = pm->mb_indcs;
  const char *slurm_nodes = std::getenv("SLURM_JOB_NUM_NODES");
  const Real inv_nranks =
      (global_variable::nranks > 0) ? static_cast<Real>(1.0/global_variable::nranks)
                                    : static_cast<Real>(1.0);

  std::cout << std::endl
            << "CGL Landau-fluid profile summary ("
            << (context != nullptr ? context : "final") << ")" << std::endl;
  std::cout << "  ranks=" << global_variable::nranks;
  if (slurm_nodes != nullptr) {
    std::cout << " slurm_nodes=" << slurm_nodes;
  }
  std::cout << " meshblocks_total=" << pm->nmb_total
            << " meshblocks_rank0=" << pmy_pack->nmb_thispack
            << " meshblock_cells=" << mb_indcs.nx1 << "x" << mb_indcs.nx2
            << "x" << mb_indcs.nx3
            << " global_cells=" << mesh_indcs.nx1 << "x" << mesh_indcs.nx2
            << "x" << mesh_indcs.nx3 << std::endl;
  std::cout << "  lf_k_parallel=" << lf_k_parallel
            << " coeff_mode=" << (lf_coeff_local ? "local" : "background")
            << " strict_admissibility=" << (strict_admissibility ? "true" : "false")
            << " backup_limiter=" << (effective_backup_limiter ? "true" : "false")
            << " profile_detail=" << (profile_detail_enabled_ ? "true" : "false")
            << " last_nstages=" << profile_last_nstages_
            << " max_nstages=" << max_nstages << std::endl;
  std::cout << "  timing: rank_mean_s is averaged over ranks; rank_max_s is the "
            << "slowest rank. heat_flux_total includes the three directional flux "
            << "regions." << std::endl;
  if (profile_detail_enabled_) {
    std::cout << "  detailed heat-flux buckets are profile-only probe kernels: "
              << "gradients, face_state, and closure replay sub-work without "
              << "updating evolved flux arrays." << std::endl;
  }
  std::cout << std::left << std::setw(34) << "bucket"
            << std::right << std::setw(16) << "rank_mean_s"
            << std::setw(16) << "rank_max_s"
            << std::setw(18) << "calls_rank_mean"
            << std::setw(16) << "calls_rank_max"
            << std::setw(18) << "max_s_per_call" << std::endl;
  for (int n = 0; n < kCGLLFProfileBucketCount; ++n) {
    if (sum_counts[n] <= 0.0) {
      continue;
    }
    const Real calls_mean = sum_counts[n]*inv_nranks;
    const Real seconds_mean = sum_seconds[n]*inv_nranks;
    const Real seconds_per_call =
        (max_counts[n] > 0.0) ? max_seconds[n]/max_counts[n] : 0.0;
    std::cout << std::left << std::setw(34)
              << CGLLFProfileBucketName(static_cast<CGLLFProfileBucket>(n))
              << std::right << std::scientific << std::setprecision(6)
              << std::setw(16) << seconds_mean
              << std::setw(16) << max_seconds[n]
              << std::fixed << std::setprecision(1)
              << std::setw(18) << calls_mean
              << std::setw(16) << max_counts[n]
              << std::scientific << std::setprecision(6)
              << std::setw(18) << seconds_per_call << std::endl;
  }
}

void CGLLandauFluid::AccumulateHeatFluxDiagnostics(const array_sum::GlobalSum &stats) {
  diagnostics.qfaces += static_cast<std::uint64_t>(stats.the_array[0]);
  diagnostics.qpar_cap += static_cast<std::uint64_t>(stats.the_array[1]);
  diagnostics.qpar_cap10 += static_cast<std::uint64_t>(stats.the_array[2]);
  diagnostics.qperp_cap += static_cast<std::uint64_t>(stats.the_array[3]);
  diagnostics.qperp_cap10 += static_cast<std::uint64_t>(stats.the_array[4]);
  stage_qpar_work_ += stats.the_array[5];
  stage_qperp_work_ += stats.the_array[6];
}

void CGLLandauFluid::AddHeatFluxes(const DvceArray5D<Real> &w,
                                   const DvceArray5D<Real> &bcc,
                                   const EOS_Data &eos, Real dt_sweep,
                                   Real rkl_weight, DvceFaceFld5D<Real> &f) {
  CGLLFProfileRegion heat_flux_profile(this, CGLLFProfileBucket::heat_flux_total);
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int ncells1 = indcs.nx1 + 2*indcs.ng;
  const int ncells2 = (indcs.nx2 > 1) ? indcs.nx2 + 2*indcs.ng : 1;
  const int ncells3 = (indcs.nx3 > 1) ? indcs.nx3 + 2*indcs.ng : 1;
  stage_qpar_work_ = 0.0;
  stage_qperp_work_ = 0.0;
  if (tpar_.extent(0) != static_cast<std::size_t>(pmy_pack->nmb_thispack) ||
      tpar_.extent(1) != static_cast<std::size_t>(ncells3) ||
      tpar_.extent(2) != static_cast<std::size_t>(ncells2) ||
      tpar_.extent(3) != static_cast<std::size_t>(ncells1)) {
    Kokkos::realloc(tpar_, pmy_pack->nmb_thispack, ncells3, ncells2, ncells1);
    Kokkos::realloc(tperp_, pmy_pack->nmb_thispack, ncells3, ncells2, ncells1);
    Kokkos::realloc(bmag_, pmy_pack->nmb_thispack, ncells3, ncells2, ncells1);
  }

  auto tpar = tpar_;
  auto tperp = tperp_;
  auto bmag = bmag_;
  {
    CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_precompute);
    par_for("cgl_lf_precompute", DevExeSpace(), 0, nmb1, 0, ncells3 - 1,
            0, ncells2 - 1, 0, ncells1 - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real rho = fmax(w(m,IDN,k,j,i), eos.dfloor);
      tpar(m,k,j,i) = w(m,IPR,k,j,i)/rho;
      tperp(m,k,j,i) = w(m,IPP,k,j,i)/rho;
      bmag(m,k,j,i) = ScaledMagneticMagnitude(
          bcc(m,IBX,k,j,i), bcc(m,IBY,k,j,i), bcc(m,IBZ,k,j,i));
    });
  }

  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  const bool multilevel = pmy_pack->pmesh->multilevel;
  auto size = pmy_pack->pmb->mb_size;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mblev = pmy_pack->pmb->mb_lev;
  const Real lf_k = lf_k_parallel;
  const bool local = lf_coeff_local;
  const Real cpar0 = lf_c_parallel0;
  const bool backup = effective_backup_limiter;
  auto &f1 = f.x1f;
  const int ni1 = ie - is + 2;
  const int nj1 = je - js + 1;
  const int nk1 = ke - ks + 1;
  const int nji1 = nj1*ni1;
  const int nkji1 = nk1*nji1;
  const int nmkji1 = (nmb1 + 1)*nkji1;
  array_sum::GlobalSum qstats1;
  {
    CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux1);
    Kokkos::parallel_reduce("cgl_lf_flux1",
        Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji1),
    KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &qstats) {
    const int m = idx/nkji1;
    const int k = (idx - m*nkji1)/nji1 + ks;
    const int j = (idx - m*nkji1 - (k - ks)*nji1)/ni1 + js;
    const int i = idx - m*nkji1 - (k - ks)*nji1 - (j - js)*ni1 + is;
    Real tx = (tpar(m,k,j,i) - tpar(m,k,j,i-1))/size.d_view(m).dx1;
    Real px = (tperp(m,k,j,i) - tperp(m,k,j,i-1))/size.d_view(m).dx1;
    Real bxg = (bmag(m,k,j,i) - bmag(m,k,j,i-1))/size.d_view(m).dx1;
    Real ty = 0.0, py = 0.0, byg = 0.0, tz = 0.0, pz = 0.0, bzg = 0.0;
    if (multi_d) {
      ty = 0.25*(tpar(m,k,j+1,i) - tpar(m,k,j-1,i) +
                 tpar(m,k,j+1,i-1) - tpar(m,k,j-1,i-1))/size.d_view(m).dx2;
      py = 0.25*(tperp(m,k,j+1,i) - tperp(m,k,j-1,i) +
                 tperp(m,k,j+1,i-1) - tperp(m,k,j-1,i-1))/size.d_view(m).dx2;
      byg = 0.25*(bmag(m,k,j+1,i) - bmag(m,k,j-1,i) +
                  bmag(m,k,j+1,i-1) - bmag(m,k,j-1,i-1))/size.d_view(m).dx2;
    }
    if (three_d) {
      tz = 0.25*(tpar(m,k+1,j,i) - tpar(m,k-1,j,i) +
                 tpar(m,k+1,j,i-1) - tpar(m,k-1,j,i-1))/size.d_view(m).dx3;
      pz = 0.25*(tperp(m,k+1,j,i) - tperp(m,k-1,j,i) +
                 tperp(m,k+1,j,i-1) - tperp(m,k-1,j,i-1))/size.d_view(m).dx3;
      bzg = 0.25*(bmag(m,k+1,j,i) - bmag(m,k-1,j,i) +
                  bmag(m,k+1,j,i-1) - bmag(m,k-1,j,i-1))/size.d_view(m).dx3;
    }
    const Real bx = 0.5*bcc(m,IBX,k,j,i-1) + 0.5*bcc(m,IBX,k,j,i);
    const Real by = 0.5*bcc(m,IBY,k,j,i-1) + 0.5*bcc(m,IBY,k,j,i);
    const Real bz = 0.5*bcc(m,IBZ,k,j,i-1) + 0.5*bcc(m,IBZ,k,j,i);
    CGLLFFaceState face;
    Real eflux = 0.0, muflux = 0.0;
    cgl_lf::ScaledValue weighted_qpar_flux, weighted_qperp_flux;
    Real qpar_ratio = 0.0, qperp_ratio = 0.0;
    if (BuildCGLLFFaceState(w(m,IDN,k,j,i-1), w(m,IDN,k,j,i),
                            w(m,IPR,k,j,i-1), w(m,IPR,k,j,i),
                            w(m,IPP,k,j,i-1), w(m,IPP,k,j,i),
                            bx, by, bz, 0, lf_k, local, cpar0, backup, eos, face)) {
      CGLLFFlux(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                dt_sweep, rkl_weight, eflux, muflux, weighted_qpar_flux,
                weighted_qperp_flux, qpar_ratio, qperp_ratio);
      if (OwnsHeatFluxDiagnosticFace(m, 0, i, is, ie, multilevel,
                                     mblev.d_view(m), nghbr)) {
        qstats.the_array[0] += 1.0;
        if (qpar_ratio > 1.0) qstats.the_array[1] += 1.0;
        if (qpar_ratio > 10.0) qstats.the_array[2] += 1.0;
        if (qperp_ratio > 1.0) qstats.the_array[3] += 1.0;
        if (qperp_ratio > 10.0) qstats.the_array[4] += 1.0;
        const Real area = size.d_view(m).dx2*size.d_view(m).dx3;
        auto qpar_work = cgl_lf::Multiply(weighted_qpar_flux, -area);
        qpar_work = cgl_lf::Multiply(
            qpar_work, tpar(m,k,j,i) - tpar(m,k,j,i-1));
        auto qperp_work = cgl_lf::Multiply(weighted_qperp_flux, -area);
        qperp_work = cgl_lf::Multiply(
            qperp_work, tperp(m,k,j,i) - tperp(m,k,j,i-1));
        qstats.the_array[5] += cgl_lf::Materialize(qpar_work);
        qstats.the_array[6] += cgl_lf::Materialize(qperp_work);
      }
    }
    f1(m,IEN,k,j,i) = eflux;
    f1(m,IAN,k,j,i) = muflux;
    }, Kokkos::Sum<array_sum::GlobalSum>(qstats1));
  }
  AccumulateHeatFluxDiagnostics(qstats1);
  if (profile_detail_enabled_) {
    Real detail = 0.0;
    {
      CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux1_gradients);
      Kokkos::parallel_reduce("cgl_lf_flux1_profile_gradients",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji1),
      KOKKOS_LAMBDA(const int idx, Real &sum) {
      const int m = idx/nkji1;
      const int k = (idx - m*nkji1)/nji1 + ks;
      const int j = (idx - m*nkji1 - (k - ks)*nji1)/ni1 + js;
      const int i = idx - m*nkji1 - (k - ks)*nji1 - (j - js)*ni1 + is;
      const Real tx = (tpar(m,k,j,i) - tpar(m,k,j,i-1))/size.d_view(m).dx1;
      const Real px = (tperp(m,k,j,i) - tperp(m,k,j,i-1))/size.d_view(m).dx1;
      const Real bxg = (bmag(m,k,j,i) - bmag(m,k,j,i-1))/size.d_view(m).dx1;
      Real ty = 0.0, py = 0.0, byg = 0.0, tz = 0.0, pz = 0.0, bzg = 0.0;
      if (multi_d) {
        ty = 0.25*(tpar(m,k,j+1,i) - tpar(m,k,j-1,i) +
                   tpar(m,k,j+1,i-1) - tpar(m,k,j-1,i-1))/size.d_view(m).dx2;
        py = 0.25*(tperp(m,k,j+1,i) - tperp(m,k,j-1,i) +
                   tperp(m,k,j+1,i-1) - tperp(m,k,j-1,i-1))/size.d_view(m).dx2;
        byg = 0.25*(bmag(m,k,j+1,i) - bmag(m,k,j-1,i) +
                    bmag(m,k,j+1,i-1) - bmag(m,k,j-1,i-1))/size.d_view(m).dx2;
      }
      if (three_d) {
        tz = 0.25*(tpar(m,k+1,j,i) - tpar(m,k-1,j,i) +
                   tpar(m,k+1,j,i-1) - tpar(m,k-1,j,i-1))/size.d_view(m).dx3;
        pz = 0.25*(tperp(m,k+1,j,i) - tperp(m,k-1,j,i) +
                   tperp(m,k+1,j,i-1) - tperp(m,k-1,j,i-1))/size.d_view(m).dx3;
        bzg = 0.25*(bmag(m,k+1,j,i) - bmag(m,k-1,j,i) +
                    bmag(m,k+1,j,i-1) - bmag(m,k-1,j,i-1))/size.d_view(m).dx3;
      }
      sum += fabs(tx) + fabs(px) + fabs(bxg) + fabs(ty) + fabs(py) + fabs(byg)
             + fabs(tz) + fabs(pz) + fabs(bzg);
      }, Kokkos::Sum<Real>(detail));
    }
    profile_detail_sink_ += detail;
    detail = 0.0;
    {
      CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux1_face_state);
      Kokkos::parallel_reduce("cgl_lf_flux1_profile_face_state",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji1),
      KOKKOS_LAMBDA(const int idx, Real &sum) {
      const int m = idx/nkji1;
      const int k = (idx - m*nkji1)/nji1 + ks;
      const int j = (idx - m*nkji1 - (k - ks)*nji1)/ni1 + js;
      const int i = idx - m*nkji1 - (k - ks)*nji1 - (j - js)*ni1 + is;
      const Real bx = 0.5*bcc(m,IBX,k,j,i-1) + 0.5*bcc(m,IBX,k,j,i);
      const Real by = 0.5*bcc(m,IBY,k,j,i-1) + 0.5*bcc(m,IBY,k,j,i);
      const Real bz = 0.5*bcc(m,IBZ,k,j,i-1) + 0.5*bcc(m,IBZ,k,j,i);
      CGLLFFaceState face;
      if (BuildCGLLFFaceState(w(m,IDN,k,j,i-1), w(m,IDN,k,j,i),
                              w(m,IPR,k,j,i-1), w(m,IPR,k,j,i),
                              w(m,IPP,k,j,i-1), w(m,IPP,k,j,i),
                              bx, by, bz, 0, lf_k, local, cpar0, backup, eos, face)) {
        sum += face.cparallel + face.bmag_inv + fabs(face.bhdir) + face.nu;
      }
      }, Kokkos::Sum<Real>(detail));
    }
    profile_detail_sink_ += detail;
    detail = 0.0;
    {
      CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux1_closure);
      Kokkos::parallel_reduce("cgl_lf_flux1_profile_closure",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji1),
      KOKKOS_LAMBDA(const int idx, Real &sum) {
      const int m = idx/nkji1;
      const int k = (idx - m*nkji1)/nji1 + ks;
      const int j = (idx - m*nkji1 - (k - ks)*nji1)/ni1 + js;
      const int i = idx - m*nkji1 - (k - ks)*nji1 - (j - js)*ni1 + is;
      const Real tx = (tpar(m,k,j,i) - tpar(m,k,j,i-1))/size.d_view(m).dx1;
      const Real px = (tperp(m,k,j,i) - tperp(m,k,j,i-1))/size.d_view(m).dx1;
      const Real bxg = (bmag(m,k,j,i) - bmag(m,k,j,i-1))/size.d_view(m).dx1;
      Real ty = 0.0, py = 0.0, byg = 0.0, tz = 0.0, pz = 0.0, bzg = 0.0;
      if (multi_d) {
        ty = 0.25*(tpar(m,k,j+1,i) - tpar(m,k,j-1,i) +
                   tpar(m,k,j+1,i-1) - tpar(m,k,j-1,i-1))/size.d_view(m).dx2;
        py = 0.25*(tperp(m,k,j+1,i) - tperp(m,k,j-1,i) +
                   tperp(m,k,j+1,i-1) - tperp(m,k,j-1,i-1))/size.d_view(m).dx2;
        byg = 0.25*(bmag(m,k,j+1,i) - bmag(m,k,j-1,i) +
                    bmag(m,k,j+1,i-1) - bmag(m,k,j-1,i-1))/size.d_view(m).dx2;
      }
      if (three_d) {
        tz = 0.25*(tpar(m,k+1,j,i) - tpar(m,k-1,j,i) +
                   tpar(m,k+1,j,i-1) - tpar(m,k-1,j,i-1))/size.d_view(m).dx3;
        pz = 0.25*(tperp(m,k+1,j,i) - tperp(m,k-1,j,i) +
                   tperp(m,k+1,j,i-1) - tperp(m,k-1,j,i-1))/size.d_view(m).dx3;
        bzg = 0.25*(bmag(m,k+1,j,i) - bmag(m,k-1,j,i) +
                    bmag(m,k+1,j,i-1) - bmag(m,k-1,j,i-1))/size.d_view(m).dx3;
      }
      const Real bx = 0.5*bcc(m,IBX,k,j,i-1) + 0.5*bcc(m,IBX,k,j,i);
      const Real by = 0.5*bcc(m,IBY,k,j,i-1) + 0.5*bcc(m,IBY,k,j,i);
      const Real bz = 0.5*bcc(m,IBZ,k,j,i-1) + 0.5*bcc(m,IBZ,k,j,i);
      CGLLFFaceState face;
      Real eflux = 0.0, muflux = 0.0;
      cgl_lf::ScaledValue weighted_qpar_flux, weighted_qperp_flux;
      Real qpar_ratio = 0.0, qperp_ratio = 0.0;
      if (BuildCGLLFFaceState(w(m,IDN,k,j,i-1), w(m,IDN,k,j,i),
                              w(m,IPR,k,j,i-1), w(m,IPR,k,j,i),
                              w(m,IPP,k,j,i-1), w(m,IPP,k,j,i),
                              bx, by, bz, 0, lf_k, local, cpar0, backup, eos, face)) {
        CGLLFFlux(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                  dt_sweep, rkl_weight, eflux, muflux, weighted_qpar_flux,
                  weighted_qperp_flux, qpar_ratio, qperp_ratio);
        sum += fabs(eflux) + fabs(muflux) + qpar_ratio + qperp_ratio;
      }
      }, Kokkos::Sum<Real>(detail));
    }
    profile_detail_sink_ += detail;
  }
  if (pmy_pack->pmesh->one_d) {
    return;
  }

  auto &f2 = f.x2f;
  const int ni2 = ie - is + 1;
  const int nj2 = je - js + 2;
  const int nk2 = ke - ks + 1;
  const int nji2 = nj2*ni2;
  const int nkji2 = nk2*nji2;
  const int nmkji2 = (nmb1 + 1)*nkji2;
  array_sum::GlobalSum qstats2;
  {
    CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux2);
    Kokkos::parallel_reduce("cgl_lf_flux2",
        Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji2),
    KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &qstats) {
    const int m = idx/nkji2;
    const int k = (idx - m*nkji2)/nji2 + ks;
    const int j = (idx - m*nkji2 - (k - ks)*nji2)/ni2 + js;
    const int i = idx - m*nkji2 - (k - ks)*nji2 - (j - js)*ni2 + is;
    const Real tx = 0.25*(tpar(m,k,j,i+1) - tpar(m,k,j,i-1) +
                          tpar(m,k,j-1,i+1) - tpar(m,k,j-1,i-1))/size.d_view(m).dx1;
    const Real px = 0.25*(tperp(m,k,j,i+1) - tperp(m,k,j,i-1) +
                          tperp(m,k,j-1,i+1) - tperp(m,k,j-1,i-1))/size.d_view(m).dx1;
    const Real bxg = 0.25*(bmag(m,k,j,i+1) - bmag(m,k,j,i-1) +
                           bmag(m,k,j-1,i+1) - bmag(m,k,j-1,i-1))/size.d_view(m).dx1;
    const Real ty = (tpar(m,k,j,i) - tpar(m,k,j-1,i))/size.d_view(m).dx2;
    const Real py = (tperp(m,k,j,i) - tperp(m,k,j-1,i))/size.d_view(m).dx2;
    const Real byg = (bmag(m,k,j,i) - bmag(m,k,j-1,i))/size.d_view(m).dx2;
    Real tz = 0.0, pz = 0.0, bzg = 0.0;
    if (three_d) {
      tz = 0.25*(tpar(m,k+1,j,i) - tpar(m,k-1,j,i) +
                 tpar(m,k+1,j-1,i) - tpar(m,k-1,j-1,i))/size.d_view(m).dx3;
      pz = 0.25*(tperp(m,k+1,j,i) - tperp(m,k-1,j,i) +
                 tperp(m,k+1,j-1,i) - tperp(m,k-1,j-1,i))/size.d_view(m).dx3;
      bzg = 0.25*(bmag(m,k+1,j,i) - bmag(m,k-1,j,i) +
                  bmag(m,k+1,j-1,i) - bmag(m,k-1,j-1,i))/size.d_view(m).dx3;
    }
    const Real bx = 0.5*bcc(m,IBX,k,j-1,i) + 0.5*bcc(m,IBX,k,j,i);
    const Real by = 0.5*bcc(m,IBY,k,j-1,i) + 0.5*bcc(m,IBY,k,j,i);
    const Real bz = 0.5*bcc(m,IBZ,k,j-1,i) + 0.5*bcc(m,IBZ,k,j,i);
    CGLLFFaceState face;
    Real eflux = 0.0, muflux = 0.0;
    cgl_lf::ScaledValue weighted_qpar_flux, weighted_qperp_flux;
    Real qpar_ratio = 0.0, qperp_ratio = 0.0;
    if (BuildCGLLFFaceState(w(m,IDN,k,j-1,i), w(m,IDN,k,j,i),
                            w(m,IPR,k,j-1,i), w(m,IPR,k,j,i),
                            w(m,IPP,k,j-1,i), w(m,IPP,k,j,i),
                            bx, by, bz, 1, lf_k, local, cpar0, backup, eos, face)) {
      CGLLFFlux(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                dt_sweep, rkl_weight, eflux, muflux, weighted_qpar_flux,
                weighted_qperp_flux, qpar_ratio, qperp_ratio);
      if (OwnsHeatFluxDiagnosticFace(m, 1, j, js, je, multilevel,
                                     mblev.d_view(m), nghbr)) {
        qstats.the_array[0] += 1.0;
        if (qpar_ratio > 1.0) qstats.the_array[1] += 1.0;
        if (qpar_ratio > 10.0) qstats.the_array[2] += 1.0;
        if (qperp_ratio > 1.0) qstats.the_array[3] += 1.0;
        if (qperp_ratio > 10.0) qstats.the_array[4] += 1.0;
        const Real area = size.d_view(m).dx1*size.d_view(m).dx3;
        auto qpar_work = cgl_lf::Multiply(weighted_qpar_flux, -area);
        qpar_work = cgl_lf::Multiply(
            qpar_work, tpar(m,k,j,i) - tpar(m,k,j-1,i));
        auto qperp_work = cgl_lf::Multiply(weighted_qperp_flux, -area);
        qperp_work = cgl_lf::Multiply(
            qperp_work, tperp(m,k,j,i) - tperp(m,k,j-1,i));
        qstats.the_array[5] += cgl_lf::Materialize(qpar_work);
        qstats.the_array[6] += cgl_lf::Materialize(qperp_work);
      }
    }
    f2(m,IEN,k,j,i) = eflux;
    f2(m,IAN,k,j,i) = muflux;
    }, Kokkos::Sum<array_sum::GlobalSum>(qstats2));
  }
  AccumulateHeatFluxDiagnostics(qstats2);
  if (profile_detail_enabled_) {
    Real detail = 0.0;
    {
      CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux2_gradients);
      Kokkos::parallel_reduce("cgl_lf_flux2_profile_gradients",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji2),
      KOKKOS_LAMBDA(const int idx, Real &sum) {
      const int m = idx/nkji2;
      const int k = (idx - m*nkji2)/nji2 + ks;
      const int j = (idx - m*nkji2 - (k - ks)*nji2)/ni2 + js;
      const int i = idx - m*nkji2 - (k - ks)*nji2 - (j - js)*ni2 + is;
      const Real tx = 0.25*(tpar(m,k,j,i+1) - tpar(m,k,j,i-1) +
                            tpar(m,k,j-1,i+1) - tpar(m,k,j-1,i-1))
                            /size.d_view(m).dx1;
      const Real px = 0.25*(tperp(m,k,j,i+1) - tperp(m,k,j,i-1) +
                            tperp(m,k,j-1,i+1) - tperp(m,k,j-1,i-1))
                            /size.d_view(m).dx1;
      const Real bxg = 0.25*(bmag(m,k,j,i+1) - bmag(m,k,j,i-1) +
                             bmag(m,k,j-1,i+1) - bmag(m,k,j-1,i-1))
                             /size.d_view(m).dx1;
      const Real ty = (tpar(m,k,j,i) - tpar(m,k,j-1,i))/size.d_view(m).dx2;
      const Real py = (tperp(m,k,j,i) - tperp(m,k,j-1,i))/size.d_view(m).dx2;
      const Real byg = (bmag(m,k,j,i) - bmag(m,k,j-1,i))/size.d_view(m).dx2;
      Real tz = 0.0, pz = 0.0, bzg = 0.0;
      if (three_d) {
        tz = 0.25*(tpar(m,k+1,j,i) - tpar(m,k-1,j,i) +
                   tpar(m,k+1,j-1,i) - tpar(m,k-1,j-1,i))/size.d_view(m).dx3;
        pz = 0.25*(tperp(m,k+1,j,i) - tperp(m,k-1,j,i) +
                   tperp(m,k+1,j-1,i) - tperp(m,k-1,j-1,i))/size.d_view(m).dx3;
        bzg = 0.25*(bmag(m,k+1,j,i) - bmag(m,k-1,j,i) +
                    bmag(m,k+1,j-1,i) - bmag(m,k-1,j-1,i))/size.d_view(m).dx3;
      }
      sum += fabs(tx) + fabs(px) + fabs(bxg) + fabs(ty) + fabs(py) + fabs(byg)
             + fabs(tz) + fabs(pz) + fabs(bzg);
      }, Kokkos::Sum<Real>(detail));
    }
    profile_detail_sink_ += detail;
    detail = 0.0;
    {
      CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux2_face_state);
      Kokkos::parallel_reduce("cgl_lf_flux2_profile_face_state",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji2),
      KOKKOS_LAMBDA(const int idx, Real &sum) {
      const int m = idx/nkji2;
      const int k = (idx - m*nkji2)/nji2 + ks;
      const int j = (idx - m*nkji2 - (k - ks)*nji2)/ni2 + js;
      const int i = idx - m*nkji2 - (k - ks)*nji2 - (j - js)*ni2 + is;
      const Real bx = 0.5*bcc(m,IBX,k,j-1,i) + 0.5*bcc(m,IBX,k,j,i);
      const Real by = 0.5*bcc(m,IBY,k,j-1,i) + 0.5*bcc(m,IBY,k,j,i);
      const Real bz = 0.5*bcc(m,IBZ,k,j-1,i) + 0.5*bcc(m,IBZ,k,j,i);
      CGLLFFaceState face;
      if (BuildCGLLFFaceState(w(m,IDN,k,j-1,i), w(m,IDN,k,j,i),
                              w(m,IPR,k,j-1,i), w(m,IPR,k,j,i),
                              w(m,IPP,k,j-1,i), w(m,IPP,k,j,i),
                              bx, by, bz, 1, lf_k, local, cpar0, backup, eos, face)) {
        sum += face.cparallel + face.bmag_inv + fabs(face.bhdir) + face.nu;
      }
      }, Kokkos::Sum<Real>(detail));
    }
    profile_detail_sink_ += detail;
    detail = 0.0;
    {
      CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux2_closure);
      Kokkos::parallel_reduce("cgl_lf_flux2_profile_closure",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji2),
      KOKKOS_LAMBDA(const int idx, Real &sum) {
      const int m = idx/nkji2;
      const int k = (idx - m*nkji2)/nji2 + ks;
      const int j = (idx - m*nkji2 - (k - ks)*nji2)/ni2 + js;
      const int i = idx - m*nkji2 - (k - ks)*nji2 - (j - js)*ni2 + is;
      const Real tx = 0.25*(tpar(m,k,j,i+1) - tpar(m,k,j,i-1) +
                            tpar(m,k,j-1,i+1) - tpar(m,k,j-1,i-1))
                            /size.d_view(m).dx1;
      const Real px = 0.25*(tperp(m,k,j,i+1) - tperp(m,k,j,i-1) +
                            tperp(m,k,j-1,i+1) - tperp(m,k,j-1,i-1))
                            /size.d_view(m).dx1;
      const Real bxg = 0.25*(bmag(m,k,j,i+1) - bmag(m,k,j,i-1) +
                             bmag(m,k,j-1,i+1) - bmag(m,k,j-1,i-1))
                             /size.d_view(m).dx1;
      const Real ty = (tpar(m,k,j,i) - tpar(m,k,j-1,i))/size.d_view(m).dx2;
      const Real py = (tperp(m,k,j,i) - tperp(m,k,j-1,i))/size.d_view(m).dx2;
      const Real byg = (bmag(m,k,j,i) - bmag(m,k,j-1,i))/size.d_view(m).dx2;
      Real tz = 0.0, pz = 0.0, bzg = 0.0;
      if (three_d) {
        tz = 0.25*(tpar(m,k+1,j,i) - tpar(m,k-1,j,i) +
                   tpar(m,k+1,j-1,i) - tpar(m,k-1,j-1,i))/size.d_view(m).dx3;
        pz = 0.25*(tperp(m,k+1,j,i) - tperp(m,k-1,j,i) +
                   tperp(m,k+1,j-1,i) - tperp(m,k-1,j-1,i))/size.d_view(m).dx3;
        bzg = 0.25*(bmag(m,k+1,j,i) - bmag(m,k-1,j,i) +
                    bmag(m,k+1,j-1,i) - bmag(m,k-1,j-1,i))/size.d_view(m).dx3;
      }
      const Real bx = 0.5*bcc(m,IBX,k,j-1,i) + 0.5*bcc(m,IBX,k,j,i);
      const Real by = 0.5*bcc(m,IBY,k,j-1,i) + 0.5*bcc(m,IBY,k,j,i);
      const Real bz = 0.5*bcc(m,IBZ,k,j-1,i) + 0.5*bcc(m,IBZ,k,j,i);
      CGLLFFaceState face;
      Real eflux = 0.0, muflux = 0.0;
      cgl_lf::ScaledValue weighted_qpar_flux, weighted_qperp_flux;
      Real qpar_ratio = 0.0, qperp_ratio = 0.0;
      if (BuildCGLLFFaceState(w(m,IDN,k,j-1,i), w(m,IDN,k,j,i),
                              w(m,IPR,k,j-1,i), w(m,IPR,k,j,i),
                              w(m,IPP,k,j-1,i), w(m,IPP,k,j,i),
                              bx, by, bz, 1, lf_k, local, cpar0, backup, eos, face)) {
        CGLLFFlux(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                  dt_sweep, rkl_weight, eflux, muflux, weighted_qpar_flux,
                  weighted_qperp_flux, qpar_ratio, qperp_ratio);
        sum += fabs(eflux) + fabs(muflux) + qpar_ratio + qperp_ratio;
      }
      }, Kokkos::Sum<Real>(detail));
    }
    profile_detail_sink_ += detail;
  }
  if (pmy_pack->pmesh->two_d) {
    return;
  }

  auto &f3 = f.x3f;
  const int ni3 = ie - is + 1;
  const int nj3 = je - js + 1;
  const int nk3 = ke - ks + 2;
  const int nji3 = nj3*ni3;
  const int nkji3 = nk3*nji3;
  const int nmkji3 = (nmb1 + 1)*nkji3;
  array_sum::GlobalSum qstats3;
  {
    CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux3);
    Kokkos::parallel_reduce("cgl_lf_flux3",
        Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji3),
    KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &qstats) {
    const int m = idx/nkji3;
    const int k = (idx - m*nkji3)/nji3 + ks;
    const int j = (idx - m*nkji3 - (k - ks)*nji3)/ni3 + js;
    const int i = idx - m*nkji3 - (k - ks)*nji3 - (j - js)*ni3 + is;
    const Real tx = 0.25*(tpar(m,k,j,i+1) - tpar(m,k,j,i-1) +
                          tpar(m,k-1,j,i+1) - tpar(m,k-1,j,i-1))/size.d_view(m).dx1;
    const Real px = 0.25*(tperp(m,k,j,i+1) - tperp(m,k,j,i-1) +
                          tperp(m,k-1,j,i+1) - tperp(m,k-1,j,i-1))/size.d_view(m).dx1;
    const Real bxg = 0.25*(bmag(m,k,j,i+1) - bmag(m,k,j,i-1) +
                           bmag(m,k-1,j,i+1) - bmag(m,k-1,j,i-1))/size.d_view(m).dx1;
    const Real ty = 0.25*(tpar(m,k,j+1,i) - tpar(m,k,j-1,i) +
                          tpar(m,k-1,j+1,i) - tpar(m,k-1,j-1,i))/size.d_view(m).dx2;
    const Real py = 0.25*(tperp(m,k,j+1,i) - tperp(m,k,j-1,i) +
                          tperp(m,k-1,j+1,i) - tperp(m,k-1,j-1,i))/size.d_view(m).dx2;
    const Real byg = 0.25*(bmag(m,k,j+1,i) - bmag(m,k,j-1,i) +
                           bmag(m,k-1,j+1,i) - bmag(m,k-1,j-1,i))/size.d_view(m).dx2;
    const Real tz = (tpar(m,k,j,i) - tpar(m,k-1,j,i))/size.d_view(m).dx3;
    const Real pz = (tperp(m,k,j,i) - tperp(m,k-1,j,i))/size.d_view(m).dx3;
    const Real bzg = (bmag(m,k,j,i) - bmag(m,k-1,j,i))/size.d_view(m).dx3;
    const Real bx = 0.5*bcc(m,IBX,k-1,j,i) + 0.5*bcc(m,IBX,k,j,i);
    const Real by = 0.5*bcc(m,IBY,k-1,j,i) + 0.5*bcc(m,IBY,k,j,i);
    const Real bz = 0.5*bcc(m,IBZ,k-1,j,i) + 0.5*bcc(m,IBZ,k,j,i);
    CGLLFFaceState face;
    Real eflux = 0.0, muflux = 0.0;
    cgl_lf::ScaledValue weighted_qpar_flux, weighted_qperp_flux;
    Real qpar_ratio = 0.0, qperp_ratio = 0.0;
    if (BuildCGLLFFaceState(w(m,IDN,k-1,j,i), w(m,IDN,k,j,i),
                            w(m,IPR,k-1,j,i), w(m,IPR,k,j,i),
                            w(m,IPP,k-1,j,i), w(m,IPP,k,j,i),
                            bx, by, bz, 2, lf_k, local, cpar0, backup, eos, face)) {
      CGLLFFlux(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                dt_sweep, rkl_weight, eflux, muflux, weighted_qpar_flux,
                weighted_qperp_flux, qpar_ratio, qperp_ratio);
      if (OwnsHeatFluxDiagnosticFace(m, 2, k, ks, ke, multilevel,
                                     mblev.d_view(m), nghbr)) {
        qstats.the_array[0] += 1.0;
        if (qpar_ratio > 1.0) qstats.the_array[1] += 1.0;
        if (qpar_ratio > 10.0) qstats.the_array[2] += 1.0;
        if (qperp_ratio > 1.0) qstats.the_array[3] += 1.0;
        if (qperp_ratio > 10.0) qstats.the_array[4] += 1.0;
        const Real area = size.d_view(m).dx1*size.d_view(m).dx2;
        auto qpar_work = cgl_lf::Multiply(weighted_qpar_flux, -area);
        qpar_work = cgl_lf::Multiply(
            qpar_work, tpar(m,k,j,i) - tpar(m,k-1,j,i));
        auto qperp_work = cgl_lf::Multiply(weighted_qperp_flux, -area);
        qperp_work = cgl_lf::Multiply(
            qperp_work, tperp(m,k,j,i) - tperp(m,k-1,j,i));
        qstats.the_array[5] += cgl_lf::Materialize(qpar_work);
        qstats.the_array[6] += cgl_lf::Materialize(qperp_work);
      }
    }
    f3(m,IEN,k,j,i) = eflux;
    f3(m,IAN,k,j,i) = muflux;
    }, Kokkos::Sum<array_sum::GlobalSum>(qstats3));
  }
  AccumulateHeatFluxDiagnostics(qstats3);
  if (profile_detail_enabled_) {
    Real detail = 0.0;
    {
      CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux3_gradients);
      Kokkos::parallel_reduce("cgl_lf_flux3_profile_gradients",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji3),
      KOKKOS_LAMBDA(const int idx, Real &sum) {
      const int m = idx/nkji3;
      const int k = (idx - m*nkji3)/nji3 + ks;
      const int j = (idx - m*nkji3 - (k - ks)*nji3)/ni3 + js;
      const int i = idx - m*nkji3 - (k - ks)*nji3 - (j - js)*ni3 + is;
      const Real tx = 0.25*(tpar(m,k,j,i+1) - tpar(m,k,j,i-1) +
                            tpar(m,k-1,j,i+1) - tpar(m,k-1,j,i-1))
                            /size.d_view(m).dx1;
      const Real px = 0.25*(tperp(m,k,j,i+1) - tperp(m,k,j,i-1) +
                            tperp(m,k-1,j,i+1) - tperp(m,k-1,j,i-1))
                            /size.d_view(m).dx1;
      const Real bxg = 0.25*(bmag(m,k,j,i+1) - bmag(m,k,j,i-1) +
                             bmag(m,k-1,j,i+1) - bmag(m,k-1,j,i-1))
                             /size.d_view(m).dx1;
      const Real ty = 0.25*(tpar(m,k,j+1,i) - tpar(m,k,j-1,i) +
                            tpar(m,k-1,j+1,i) - tpar(m,k-1,j-1,i))
                            /size.d_view(m).dx2;
      const Real py = 0.25*(tperp(m,k,j+1,i) - tperp(m,k,j-1,i) +
                            tperp(m,k-1,j+1,i) - tperp(m,k-1,j-1,i))
                            /size.d_view(m).dx2;
      const Real byg = 0.25*(bmag(m,k,j+1,i) - bmag(m,k,j-1,i) +
                             bmag(m,k-1,j+1,i) - bmag(m,k-1,j-1,i))
                             /size.d_view(m).dx2;
      const Real tz = (tpar(m,k,j,i) - tpar(m,k-1,j,i))/size.d_view(m).dx3;
      const Real pz = (tperp(m,k,j,i) - tperp(m,k-1,j,i))/size.d_view(m).dx3;
      const Real bzg = (bmag(m,k,j,i) - bmag(m,k-1,j,i))/size.d_view(m).dx3;
      sum += fabs(tx) + fabs(px) + fabs(bxg) + fabs(ty) + fabs(py) + fabs(byg)
             + fabs(tz) + fabs(pz) + fabs(bzg);
      }, Kokkos::Sum<Real>(detail));
    }
    profile_detail_sink_ += detail;
    detail = 0.0;
    {
      CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux3_face_state);
      Kokkos::parallel_reduce("cgl_lf_flux3_profile_face_state",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji3),
      KOKKOS_LAMBDA(const int idx, Real &sum) {
      const int m = idx/nkji3;
      const int k = (idx - m*nkji3)/nji3 + ks;
      const int j = (idx - m*nkji3 - (k - ks)*nji3)/ni3 + js;
      const int i = idx - m*nkji3 - (k - ks)*nji3 - (j - js)*ni3 + is;
      const Real bx = 0.5*bcc(m,IBX,k-1,j,i) + 0.5*bcc(m,IBX,k,j,i);
      const Real by = 0.5*bcc(m,IBY,k-1,j,i) + 0.5*bcc(m,IBY,k,j,i);
      const Real bz = 0.5*bcc(m,IBZ,k-1,j,i) + 0.5*bcc(m,IBZ,k,j,i);
      CGLLFFaceState face;
      if (BuildCGLLFFaceState(w(m,IDN,k-1,j,i), w(m,IDN,k,j,i),
                              w(m,IPR,k-1,j,i), w(m,IPR,k,j,i),
                              w(m,IPP,k-1,j,i), w(m,IPP,k,j,i),
                              bx, by, bz, 2, lf_k, local, cpar0, backup, eos, face)) {
        sum += face.cparallel + face.bmag_inv + fabs(face.bhdir) + face.nu;
      }
      }, Kokkos::Sum<Real>(detail));
    }
    profile_detail_sink_ += detail;
    detail = 0.0;
    {
      CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_flux3_closure);
      Kokkos::parallel_reduce("cgl_lf_flux3_profile_closure",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji3),
      KOKKOS_LAMBDA(const int idx, Real &sum) {
      const int m = idx/nkji3;
      const int k = (idx - m*nkji3)/nji3 + ks;
      const int j = (idx - m*nkji3 - (k - ks)*nji3)/ni3 + js;
      const int i = idx - m*nkji3 - (k - ks)*nji3 - (j - js)*ni3 + is;
      const Real tx = 0.25*(tpar(m,k,j,i+1) - tpar(m,k,j,i-1) +
                            tpar(m,k-1,j,i+1) - tpar(m,k-1,j,i-1))
                            /size.d_view(m).dx1;
      const Real px = 0.25*(tperp(m,k,j,i+1) - tperp(m,k,j,i-1) +
                            tperp(m,k-1,j,i+1) - tperp(m,k-1,j,i-1))
                            /size.d_view(m).dx1;
      const Real bxg = 0.25*(bmag(m,k,j,i+1) - bmag(m,k,j,i-1) +
                             bmag(m,k-1,j,i+1) - bmag(m,k-1,j,i-1))
                             /size.d_view(m).dx1;
      const Real ty = 0.25*(tpar(m,k,j+1,i) - tpar(m,k,j-1,i) +
                            tpar(m,k-1,j+1,i) - tpar(m,k-1,j-1,i))
                            /size.d_view(m).dx2;
      const Real py = 0.25*(tperp(m,k,j+1,i) - tperp(m,k,j-1,i) +
                            tperp(m,k-1,j+1,i) - tperp(m,k-1,j-1,i))
                            /size.d_view(m).dx2;
      const Real byg = 0.25*(bmag(m,k,j+1,i) - bmag(m,k,j-1,i) +
                             bmag(m,k-1,j+1,i) - bmag(m,k-1,j-1,i))
                             /size.d_view(m).dx2;
      const Real tz = (tpar(m,k,j,i) - tpar(m,k-1,j,i))/size.d_view(m).dx3;
      const Real pz = (tperp(m,k,j,i) - tperp(m,k-1,j,i))/size.d_view(m).dx3;
      const Real bzg = (bmag(m,k,j,i) - bmag(m,k-1,j,i))/size.d_view(m).dx3;
      const Real bx = 0.5*bcc(m,IBX,k-1,j,i) + 0.5*bcc(m,IBX,k,j,i);
      const Real by = 0.5*bcc(m,IBY,k-1,j,i) + 0.5*bcc(m,IBY,k,j,i);
      const Real bz = 0.5*bcc(m,IBZ,k-1,j,i) + 0.5*bcc(m,IBZ,k,j,i);
      CGLLFFaceState face;
      Real eflux = 0.0, muflux = 0.0;
      cgl_lf::ScaledValue weighted_qpar_flux, weighted_qperp_flux;
      Real qpar_ratio = 0.0, qperp_ratio = 0.0;
      if (BuildCGLLFFaceState(w(m,IDN,k-1,j,i), w(m,IDN,k,j,i),
                              w(m,IPR,k-1,j,i), w(m,IPR,k,j,i),
                              w(m,IPP,k-1,j,i), w(m,IPP,k,j,i),
                              bx, by, bz, 2, lf_k, local, cpar0, backup, eos, face)) {
        CGLLFFlux(face, tx, ty, tz, px, py, pz, bxg, byg, bzg,
                  dt_sweep, rkl_weight, eflux, muflux, weighted_qpar_flux,
                  weighted_qperp_flux, qpar_ratio, qperp_ratio);
        sum += fabs(eflux) + fabs(muflux) + qpar_ratio + qperp_ratio;
      }
      }, Kokkos::Sum<Real>(detail));
    }
    profile_detail_sink_ += detail;
  }
}

void CGLLandauFluid::AdvanceHeatFluxWorkDiagnostics(
    const parabolic::RKL2Coefficients &coeffs, int stage, int nstages) {
  CGLLFProfileRegion profile(this, CGLLFProfileBucket::heat_flux_work_diagnostics);
  if (profile_enabled_) {
    profile_last_nstages_ = nstages;
    if (nstages > profile_max_nstages_) {
      profile_max_nstages_ = nstages;
    }
  }
  if (stage == 1) {
    sweep_qpar_work_ = 0.0;
    sweep_qperp_work_ = 0.0;
    sweep_qpar_work1_ = 0.0;
    sweep_qperp_work1_ = 0.0;
    sweep_qpar_work2_ = 0.0;
    sweep_qperp_work2_ = 0.0;
    sweep_qpar_rhs_ = 0.0;
    sweep_qperp_rhs_ = 0.0;
  }

  const Real qpar_rhs = stage_qpar_work_;
  const Real qperp_rhs = stage_qperp_work_;
  sweep_qpar_work2_ = sweep_qpar_work1_;
  sweep_qperp_work2_ = sweep_qperp_work1_;
  sweep_qpar_work1_ = sweep_qpar_work_;
  sweep_qperp_work1_ = sweep_qperp_work_;
  // The per-sweep diagnostic starts at zero, so the RKL2 S0 term vanishes.
  const Real first_rhs_coeff =
      cgl_lf::CachedRHSCoefficient(coeffs.gammaj_tilde, nstages);
  sweep_qpar_work_ = cgl_lf::WeightedRKL2Update(
      coeffs.muj, sweep_qpar_work1_, coeffs.nuj, sweep_qpar_work2_,
      0.0, 0.0, first_rhs_coeff, sweep_qpar_rhs_, qpar_rhs);
  sweep_qperp_work_ = cgl_lf::WeightedRKL2Update(
      coeffs.muj, sweep_qperp_work1_, coeffs.nuj, sweep_qperp_work2_,
      0.0, 0.0, first_rhs_coeff, sweep_qperp_rhs_, qperp_rhs);
  if (stage == 1) {
    sweep_qpar_rhs_ = qpar_rhs;
    sweep_qperp_rhs_ = qperp_rhs;
  }
  if (stage == nstages) {
    diagnostics.qpar_work += sweep_qpar_work_;
    diagnostics.qperp_work += sweep_qperp_work_;
  }
}

void CGLLandauFluid::AdvancePressureWorkDiagnostics(Real beta_dt, Real gam0, Real gam1,
                                                     int stage, Real pressure_power,
                                                     Real anisotropic_power) {
  if (stage == 1) {
    pressure_work_cycle_start_ = diagnostics.pressure_work;
    anisotropic_pressure_work_cycle_start_ = diagnostics.anisotropic_pressure_work;
  }
  diagnostics.pressure_work =
      gam0*diagnostics.pressure_work + gam1*pressure_work_cycle_start_
      + beta_dt*pressure_power;
  diagnostics.anisotropic_pressure_work =
      gam0*diagnostics.anisotropic_pressure_work
      + gam1*anisotropic_pressure_work_cycle_start_
      + beta_dt*anisotropic_power;
}

void CGLLandauFluid::NewTimeStep(const DvceArray5D<Real> &w, const EOS_Data &eos) {
  CGLLFProfileRegion profile(this, CGLLFProfileBucket::timestep_reduction);
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, nx1 = indcs.nx1;
  const int js = indcs.js, nx2 = indcs.nx2;
  const int ks = indcs.ks, nx3 = indcs.nx3;
  const int nmkji = pmy_pack->nmb_thispack*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji = nx2*nx1;
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  const Real kpar = lf_k_parallel;
  const bool local = lf_coeff_local;
  const Real cpar0 = lf_c_parallel0;
  auto size = pmy_pack->pmb->mb_size;
  dtnew = static_cast<Real>(std::numeric_limits<float>::max());
  Kokkos::parallel_reduce("cgl_lf_newdt", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, Real &min_dt) {
    const int m = idx/nkji;
    int k = (idx - m*nkji)/nji + ks;
    int j = (idx - m*nkji - (k - ks)*nji)/nx1 + js;
    const int i = idx - m*nkji - (k - ks)*nji - (j - js)*nx1 + is;
    const Real rho = fmax(w(m,IDN,k,j,i), eos.dfloor);
    const Real cpar = local ? sqrt(fmax(w(m,IPR,k,j,i)/rho, eos.tfloor)) : cpar0;
    const Real chi = cgl::kSqrtEightOverPi*cpar/kpar;
    if (chi > 0.0) {
      min_dt = fmin(min_dt, SQR(size.d_view(m).dx1)/chi);
      if (multi_d) min_dt = fmin(min_dt, SQR(size.d_view(m).dx2)/chi);
      if (three_d) min_dt = fmin(min_dt, SQR(size.d_view(m).dx3)/chi);
    }
  }, Kokkos::Min<Real>(dtnew));
  const Real fac = three_d ? static_cast<Real>(1.0/6.0) :
                   (multi_d ? static_cast<Real>(0.25) : static_cast<Real>(0.5));
  dtnew *= fac;
}

void CGLLandauFluid::RecordAdmissibility(const DvceArray5D<Real> &u,
                                         const DvceArray5D<Real> &w,
                                         const DvceArray5D<Real> &bcc,
                                         const EOS_Data &eos,
                                         int dfloor_delta, int pfloor_delta,
                                         const char *sweep_name,
                                         int stage, int nstages) {
  CGLLFProfileRegion profile(this, CGLLFProfileBucket::admissibility);
  if (profile_enabled_) {
    profile_last_nstages_ = nstages;
    if (nstages > profile_max_nstages_) {
      profile_max_nstages_ = nstages;
    }
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, nx1 = indcs.nx1;
  const int js = indcs.js, nx2 = indcs.nx2;
  const int ks = indcs.ks, nx3 = indcs.nx3;
  const int nmkji = pmy_pack->nmb_thispack*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji = nx2*nx1;
  int nonfinite = 0, nonpositive = 0, mirror = 0, firehose = 0, hard_bound = 0;
  Kokkos::parallel_reduce("cgl_lf_admissibility",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int idx, int &nbad, int &nnonpos, int &nmirror,
                int &nfirehose, int &nhard) {
    const int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    const int i = idx - m*nkji - k*nji - j*nx1 + is;
    k += ks;
    j += js;
    const Real rho = w(m,IDN,k,j,i);
    const Real ppar = w(m,IPR,k,j,i);
    const Real pperp = w(m,IPP,k,j,i);
    const Real energy = u(m,IEN,k,j,i);
    const Real moment = u(m,IAN,k,j,i);
    if (!Kokkos::isfinite(rho) || !Kokkos::isfinite(ppar) ||
        !Kokkos::isfinite(pperp) || !Kokkos::isfinite(energy) ||
        !Kokkos::isfinite(moment)) {
      ++nbad;
    }
    if (rho <= 0.0 || ppar <= 0.0 || pperp <= 0.0) {
      ++nnonpos;
    }
    const Real bsqr = SQR(bcc(m,IBX,k,j,i)) + SQR(bcc(m,IBY,k,j,i))
                     + SQR(bcc(m,IBZ,k,j,i));
    const Real paniso = pperp - ppar;
    if (eos.mlim && cgl::MirrorLimiterActive(paniso, bsqr)) {
      ++nmirror;
      if (cgl::MirrorHardBoundViolated(paniso, bsqr)) {
        ++nhard;
      }
    }
    if (eos.flim &&
        cgl::FirehoseLimiterActive(paniso, bsqr, eos.firehose_threshold)) {
      ++nfirehose;
      if (cgl::FirehoseHardBoundViolated(paniso, bsqr)) {
        ++nhard;
      }
    }
  }, Kokkos::Sum<int>(nonfinite), Kokkos::Sum<int>(nonpositive),
     Kokkos::Sum<int>(mirror), Kokkos::Sum<int>(firehose),
     Kokkos::Sum<int>(hard_bound));

  diagnostics.nstage += static_cast<std::uint64_t>(nmkji);
  diagnostics.dfloor += static_cast<std::uint64_t>(dfloor_delta);
  diagnostics.pfloor += static_cast<std::uint64_t>(pfloor_delta);
  diagnostics.nonfinite += static_cast<std::uint64_t>(nonfinite);
  diagnostics.nonpositive += static_cast<std::uint64_t>(nonpositive);
  diagnostics.mirror += static_cast<std::uint64_t>(mirror);
  diagnostics.firehose += static_cast<std::uint64_t>(firehose);
  diagnostics.hard_bound += static_cast<std::uint64_t>(hard_bound);

  if (strict_admissibility &&
      (dfloor_delta > 0 || pfloor_delta > 0 || nonfinite > 0 ||
       nonpositive > 0 || hard_bound > 0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "CGL Landau-fluid strict admissibility failed after a split stage: "
              << "sweep=" << sweep_name << " stage=" << stage << "/" << nstages
              << " dfloor=" << dfloor_delta << " pfloor=" << pfloor_delta
              << " nonfinite=" << nonfinite << " nonpositive=" << nonpositive
              << " hard_bound=" << hard_bound << std::endl;
    std::exit(EXIT_FAILURE);
  }
}
