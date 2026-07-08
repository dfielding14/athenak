//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_sts.cpp
//! \brief MHD-owned helpers and task wrappers for super time stepping.

#include <cstdlib>
#include <iostream>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "eos/eos.hpp"
#include "diffusion/conduction.hpp"
#include "diffusion/cgl_landau_fluid.hpp"
#include "diffusion/cgl_landau_fluid_arithmetic.hpp"
#include "diffusion/resistivity.hpp"
#include "diffusion/scalar_diffusion.hpp"
#include "diffusion/viscosity.hpp"
#include "diffusion/hyperviscosity.hpp"
#include "globals.hpp"
#include "mhd.hpp"

namespace {

bool CGLLFTaskTraceEnabled() {
  const char *env = std::getenv("ATHENAK_CGL_LF_TASK_TRACE");
  return env != nullptr && env[0] != '\0' && env[0] != '0';
}

void TraceCGLLFTask(MeshBlockPack *pmbp, const char *task, const char *point,
                    int stage) {
  if (!CGLLFTaskTraceEnabled()) return;
  Kokkos::fence();
  auto *pm = pmbp->pmesh;
  std::cout << "[cgl_lf_task] rank=" << global_variable::my_rank
            << " cycle=" << pm->ncycle
            << " time=" << pm->time
            << " task=" << task
            << " point=" << point
            << " stage=" << stage
            << " nmb=" << pmbp->nmb_thispack
            << std::endl;
}

KOKKOS_INLINE_FUNCTION
bool UpdateSTSMHDVariable(const int n, const bool update_momentum,
                          const bool update_energy, const bool update_cgl_moment,
                          const bool update_scalars,
                          const int nmhd) {
  if (update_momentum && (n == IVX || n == IVY || n == IVZ)) {
    return true;
  }
  if (update_energy && n == IEN) {
    return true;
  }
  if (update_cgl_moment && n == IAN) {
    return true;
  }
  if (update_scalars && n >= nmhd) {
    return true;
  }
  return false;
}

const char *STSSweepName(const Driver *pdrive) {
  if (pdrive->sts.sweep == Driver::STSSweep::pre) return "pre";
  if (pdrive->sts.sweep == Driver::STSSweep::post) return "post";
  return "none";
}

void RefreshCellCenteredBFromFace(const MeshBlockPack *pmbp,
                                  const DvceFaceFld4D<Real> &b,
                                  DvceArray5D<Real> &bcc,
                                  const int il, const int iu,
                                  const int jl, const int ju,
                                  const int kl, const int ku) {
  const int nmb = pmbp->nmb_thispack;
  auto b1 = b.x1f;
  auto b2 = b.x2f;
  auto b3 = b.x3f;
  auto bcc_ = bcc;

  par_for("mhd_refresh_bcc_from_face", DevExeSpace(), 0, nmb - 1,
          kl, ku, jl, ju, il, iu,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    bcc_(m,IBX,k,j,i) = 0.5*(b1(m,k,j,i) + b1(m,k,j,i+1));
    bcc_(m,IBY,k,j,i) = 0.5*(b2(m,k,j,i) + b2(m,k,j+1,i));
    bcc_(m,IBZ,k,j,i) = 0.5*(b3(m,k,j,i) + b3(m,k+1,j,i));
  });
}

} // namespace

namespace mhd {

//! \fn void MHD::AddSelectedDiffusionFluxes()
//! \brief Add only the requested subset of MHD diffusion operators to the live flux
//! scratch array.

void MHD::AddSelectedDiffusionFluxes(DiffusionSelection selection,
                                     Real cgl_dt_sweep, Real cgl_rkl_weight) {
  const bool add_viscosity =
      (selection == DiffusionSelection::explicit_only) ? has_explicit_viscosity
                                                       : has_sts_viscosity;
  const bool add_hyperviscosity =
      (selection == DiffusionSelection::explicit_only) ? has_explicit_hyperviscosity
                                                       : has_sts_hyperviscosity;
  const bool add_conduction =
      (selection == DiffusionSelection::explicit_only) ? has_explicit_conduction
                                                       : has_sts_conduction;
  const bool add_resistivity =
      (selection == DiffusionSelection::explicit_only) ? has_explicit_resistivity
                                                       : has_sts_resistivity;
  const bool add_scalar_diffusion =
      (selection == DiffusionSelection::explicit_only) ? has_explicit_scalar_diffusion
                                                       : has_sts_scalar_diffusion;

  if (add_viscosity && pvisc != nullptr) {
    pvisc->IsotropicViscousFlux(w0, pvisc->nu_iso, peos->eos_data, uflx);
  }
  if (add_hyperviscosity && phypervisc != nullptr) {
    phypervisc->AddHyperViscousFlux(w0, peos->eos_data, uflx);
  }
  if (add_resistivity && presist != nullptr && peos->eos_data.is_ideal) {
    presist->OhmicEnergyFlux(b0, uflx);
  }
  if (add_conduction && pcond != nullptr) {
    pcond->AddHeatFlux(w0, peos->eos_data, uflx);
  }
  if (selection == DiffusionSelection::sts_only && has_cgl_lf_split &&
      pcgl_lf != nullptr) {
    pcgl_lf->AddHeatFluxes(w0, bcc0, peos->eos_data, cgl_dt_sweep,
                           cgl_rkl_weight, uflx);
  }
  if (add_scalar_diffusion && pscalar_diff != nullptr) {
    pscalar_diff->AddScalarDiffusionFlux(w0, nmhd, nscalars, uflx);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void MHD::AddSelectedDiffusionEMF()
//! \brief Add only the requested subset of MHD resistive EMFs to the live electric-field
//! scratch array.

void MHD::AddSelectedDiffusionEMF(DiffusionSelection selection) {
  const bool add_resistivity =
      (selection == DiffusionSelection::explicit_only) ? has_explicit_resistivity
                                                       : has_sts_resistivity;
  if (add_resistivity && presist != nullptr && presist->eta_ohm > 0.0) {
    presist->OhmicEField(b0, efld);
  }
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::ClearSTSFlux()
//! \brief Zero the MHD conserved-flux scratch before one STS stage.

TaskStatus MHD::ClearSTSFlux(Driver *pdrive, int stage) {
  (void) pdrive;
  (void) stage;
  TraceCGLLFTask(pmy_pack, "ClearSTSFlux", "begin", stage);
  CGLLFProfileRegion profile(pcgl_lf, CGLLFProfileBucket::sts_clear_flux);
  const bool lf_only_sts_cell_update =
      has_sts_cgl_lf && pcgl_lf != nullptr && !has_sts_viscosity &&
      !has_sts_hyperviscosity && !has_sts_conduction && !has_sts_resistivity &&
      !has_sts_scalar_diffusion;
  if (lf_only_sts_cell_update) {
    auto flx1 = uflx.x1f;
    auto flx2 = uflx.x2f;
    auto flx3 = uflx.x3f;
    const int nmb = pmy_pack->nmb_thispack;
    const int n31 = static_cast<int>(flx1.extent(2));
    const int n21 = static_cast<int>(flx1.extent(3));
    const int n11 = static_cast<int>(flx1.extent(4));
    const int n32 = static_cast<int>(flx2.extent(2));
    const int n22 = static_cast<int>(flx2.extent(3));
    const int n12 = static_cast<int>(flx2.extent(4));
    const int n33 = static_cast<int>(flx3.extent(2));
    const int n23 = static_cast<int>(flx3.extent(3));
    const int n13 = static_cast<int>(flx3.extent(4));
    par_for("mhd_sts_clear_cgl_lf_flux1", DevExeSpace(), 0, nmb - 1, 0, 1,
            0, n31 - 1, 0, n21 - 1, 0, n11 - 1,
    KOKKOS_LAMBDA(const int m, const int q, const int k, const int j,
                  const int i) {
      const int n = (q == 0) ? IEN : IAN;
      flx1(m,n,k,j,i) = 0.0;
    });
    par_for("mhd_sts_clear_cgl_lf_flux2", DevExeSpace(), 0, nmb - 1, 0, 1,
            0, n32 - 1, 0, n22 - 1, 0, n12 - 1,
    KOKKOS_LAMBDA(const int m, const int q, const int k, const int j,
                  const int i) {
      const int n = (q == 0) ? IEN : IAN;
      flx2(m,n,k,j,i) = 0.0;
    });
    par_for("mhd_sts_clear_cgl_lf_flux3", DevExeSpace(), 0, nmb - 1, 0, 1,
            0, n33 - 1, 0, n23 - 1, 0, n13 - 1,
    KOKKOS_LAMBDA(const int m, const int q, const int k, const int j,
                  const int i) {
      const int n = (q == 0) ? IEN : IAN;
      flx3(m,n,k,j,i) = 0.0;
    });
    TraceCGLLFTask(pmy_pack, "ClearSTSFlux", "end", stage);
    return TaskStatus::complete;
  }
  Kokkos::deep_copy(DevExeSpace(), uflx.x1f, 0.0);
  Kokkos::deep_copy(DevExeSpace(), uflx.x2f, 0.0);
  Kokkos::deep_copy(DevExeSpace(), uflx.x3f, 0.0);
  TraceCGLLFTask(pmy_pack, "ClearSTSFlux", "end", stage);
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::ClearSTSEField()
//! \brief Zero the MHD electric-field scratch before one STS stage.

TaskStatus MHD::ClearSTSEField(Driver *pdrive, int stage) {
  (void) pdrive;
  (void) stage;
  Kokkos::deep_copy(DevExeSpace(), efld.x1e, 0.0);
  Kokkos::deep_copy(DevExeSpace(), efld.x2e, 0.0);
  Kokkos::deep_copy(DevExeSpace(), efld.x3e, 0.0);
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::STSFluxes()
//! \brief Accumulate only the STS-managed MHD diffusion operators that contribute to U.

TaskStatus MHD::STSFluxes(Driver *pdrive, int stage) {
  TraceCGLLFTask(pmy_pack, "STSFluxes", "begin", stage);
  if (!has_any_parabolic_cell_update) {
    TraceCGLLFTask(pmy_pack, "STSFluxes", "skip", stage);
    return TaskStatus::complete;
  }

  AddSelectedDiffusionFluxes(
      DiffusionSelection::sts_only, pdrive->sts.dt_sweep,
      pdrive->sts.coeffs.muj_tilde);
  TraceCGLLFTask(pmy_pack, "STSFluxes", "end", stage);
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::BeginCGLLandauFluidSTSSweep()
//! \brief Switch the CGL extra state from anisotropy to magnetic moment for LF STS.

TaskStatus MHD::BeginCGLLandauFluidSTSSweep(Driver *pdrive, int stage) {
  TraceCGLLFTask(pmy_pack, "BeginCGLLandauFluidSTSSweep", "begin", stage);
  if (!has_cgl_lf_split || !pdrive->sts.enabled || stage != 1) {
    TraceCGLLFTask(pmy_pack, "BeginCGLLandauFluidSTSSweep", "skip", stage);
    return TaskStatus::complete;
  }
  RequireCGLAnisotropyRepresentation("CGL Landau-fluid sweep begin");
  if (diagnose_nonfinite_rk_update) {
    DiagnoseNonfiniteCGLState(stage, "anisotropy-to-magnetic-moment", "pre",
                              STSSweepName(pdrive), "anisotropy", u0);
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int ng = indcs.ng;
  const int n1m1 = indcs.nx1 + 2*ng - 1;
  const int n2m1 = (indcs.nx2 > 1) ? indcs.nx2 + 2*ng - 1 : 0;
  const int n3m1 = (indcs.nx3 > 1) ? indcs.nx3 + 2*ng - 1 : 0;
  {
    CGLLFProfileRegion profile(pcgl_lf, CGLLFProfileBucket::sweep_begin_conversion);
    peos->CGLAnisotropyToMagneticMoment(u0, bcc0, 0, n1m1, 0, n2m1, 0, n3m1);
  }
  cgl_slot_representation = CGLSlotRepresentation::magnetic_moment;
  if (diagnose_nonfinite_rk_update) {
    DiagnoseNonfiniteCGLState(stage, "anisotropy-to-magnetic-moment", "post",
                              STSSweepName(pdrive), "magnetic-moment", u0);
  }
  TraceCGLLFTask(pmy_pack, "BeginCGLLandauFluidSTSSweep", "end", stage);
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::STSEField()
//! \brief Accumulate only the STS-managed resistive EMF contribution.

TaskStatus MHD::STSEField(Driver *pdrive, int stage) {
  (void) pdrive;
  (void) stage;
  if (!has_any_parabolic_field_update) {
    return TaskStatus::complete;
  }

  AddSelectedDiffusionEMF(DiffusionSelection::sts_only);
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::STSUpdateU()
//! \brief Apply one MHD-owned RKL2 STS stage over the enrolled conserved variables.

TaskStatus MHD::STSUpdateU(Driver *pdrive, int stage) {
  TraceCGLLFTask(pmy_pack, "STSUpdateU", "begin", stage);
  if (!has_any_parabolic_cell_update || !(pdrive->sts.enabled)) {
    TraceCGLLFTask(pmy_pack, "STSUpdateU", "skip", stage);
    return TaskStatus::complete;
  }

  if (diagnose_nonfinite_rk_update && has_cgl_lf_split) {
    DiagnoseNonfiniteCGLState(stage, "sts-update-u", "pre", STSSweepName(pdrive),
                              "magnetic-moment", u0);
  }

  const bool lf_only_sts_cell_update =
      has_sts_cgl_lf && pcgl_lf != nullptr && !has_sts_viscosity &&
      !has_sts_hyperviscosity && !has_sts_conduction && !has_sts_resistivity &&
      !has_sts_scalar_diffusion;

  {
    CGLLFProfileRegion profile(pcgl_lf, CGLLFProfileBucket::sts_update_copies);
    if (lf_only_sts_cell_update) {
      const int nmb = pmy_pack->nmb_thispack;
      const int n3 = static_cast<int>(u0.extent(2));
      const int n2 = static_cast<int>(u0.extent(3));
      const int n1 = static_cast<int>(u0.extent(4));
      auto u0_ = u0;
      auto u_sts0_ = u_sts0;
      auto u_sts1_ = u_sts1;
      auto u_sts2_ = u_sts2;
      par_for("mhd_sts_copy_cgl_lf_u", DevExeSpace(), 0, nmb - 1, 0, 1,
              0, n3 - 1, 0, n2 - 1, 0, n1 - 1,
      KOKKOS_LAMBDA(const int m, const int q, const int k, const int j,
                    const int i) {
        const int n = (q == 0) ? IEN : IAN;
        if (stage == 1) {
          u_sts0_(m,n,k,j,i) = u0_(m,n,k,j,i);
        }
        u_sts2_(m,n,k,j,i) = u_sts1_(m,n,k,j,i);
        u_sts1_(m,n,k,j,i) = u0_(m,n,k,j,i);
      });
    } else {
      if (stage == 1) {
        Kokkos::deep_copy(DevExeSpace(), u_sts0, u0);
      }
      Kokkos::deep_copy(DevExeSpace(), u_sts2, u_sts1);
      Kokkos::deep_copy(DevExeSpace(), u_sts1, u0);
    }
  }

  const bool update_momentum = (has_sts_viscosity || has_sts_hyperviscosity);
  const bool update_energy = (has_sts_conduction || has_cgl_lf_split ||
                              ((has_sts_viscosity || has_sts_hyperviscosity ||
                                has_sts_resistivity) &&
                               peos->eos_data.is_ideal));
  const bool update_cgl_moment = has_cgl_lf_split;
  const bool update_scalars = has_sts_scalar_diffusion;
  if (!(update_momentum || update_energy || update_cgl_moment || update_scalars)) {
    TraceCGLLFTask(pmy_pack, "STSUpdateU", "skip", stage);
    return TaskStatus::complete;
  }

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;

  int nmb1 = pmy_pack->nmb_thispack - 1;
  int nvars = nmhd + nscalars;
  int nmhd_vars = nmhd;
  Real dt_sweep = pdrive->sts.dt_sweep;
  auto coeffs = pdrive->sts.coeffs;
  const bool cgl_lf_weighted_flux =
      has_cgl_lf_split && pcgl_lf != nullptr && pcgl_lf->UsesWeightedSTSFlux();
  const Real cgl_first_rkl_weight =
      cgl_lf_weighted_flux
          ? cgl_lf::FirstStageRKLWeight(pdrive->sts.nstages)
          : static_cast<Real>(1.0);
  auto u0_ = u0;
  auto u_sts0_ = u_sts0;
  auto u_sts1_ = u_sts1;
  auto u_sts2_ = u_sts2;
  auto u_sts_rhs_ = u_sts_rhs;
  auto flx1 = uflx.x1f;
  auto flx2 = uflx.x2f;
  auto flx3 = uflx.x3f;
  auto mbsize = pmy_pack->pmb->mb_size;

  int scr_level = 0;
  size_t scr_size = ScrArray1D<Real>::shmem_size(ncells1);
  const int nvars_update = lf_only_sts_cell_update ? 2 : nvars;

  {
    CGLLFProfileRegion profile(pcgl_lf, CGLLFProfileBucket::sts_update_kernel);
    if (lf_only_sts_cell_update) {
      par_for("mhd_sts_update_cgl_lf_u", DevExeSpace(), 0, nmb1, 0, 1,
              ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int q, const int k, const int j,
                    const int i) {
        const int n = (q == 0) ? IEN : IAN;
        Real divf =
            (flx1(m,n,k,j,i+1) - flx1(m,n,k,j,i))/mbsize.d_view(m).dx1;
        if (multi_d) {
          divf +=
              (flx2(m,n,k,j+1,i) - flx2(m,n,k,j,i))/mbsize.d_view(m).dx2;
        }
        if (three_d) {
          divf +=
              (flx3(m,n,k+1,j,i) - flx3(m,n,k,j,i))/mbsize.d_view(m).dx3;
        }
        if (cgl_lf_weighted_flux) {
          // The CGL LF face flux already contains dt_sweep*muj_tilde.
          const Real weighted_rhs = -divf;
          const Real first_rhs_coeff =
              coeffs.gammaj_tilde/cgl_first_rkl_weight;
          u0_(m,n,k,j,i) = cgl_lf::WeightedRKL2Update(
              coeffs.muj, u_sts1_(m,n,k,j,i),
              coeffs.nuj, u_sts2_(m,n,k,j,i),
              1.0 - coeffs.muj - coeffs.nuj, u_sts0_(m,n,k,j,i),
              first_rhs_coeff, u_sts_rhs_(m,n,k,j,i), weighted_rhs);
          if (stage == 1) {
            u_sts_rhs_(m,n,k,j,i) = weighted_rhs;
          }
        } else {
          const Real delta_u = -dt_sweep*divf;
          u0_(m,n,k,j,i) = coeffs.muj*u_sts1_(m,n,k,j,i)
                         + coeffs.nuj*u_sts2_(m,n,k,j,i)
                         + (1.0 - coeffs.muj - coeffs.nuj)*u_sts0_(m,n,k,j,i)
                         + coeffs.gammaj_tilde*u_sts_rhs_(m,n,k,j,i)
                         + coeffs.muj_tilde*delta_u;
          if (stage == 1) {
            u_sts_rhs_(m,n,k,j,i) = delta_u;
          }
        }
      });
    } else {
      par_for_outer("mhd_sts_update_u", DevExeSpace(), scr_size, scr_level, 0, nmb1,
                    0, nvars_update - 1, ks, ke, js, je,
      KOKKOS_LAMBDA(TeamMember_t member, const int m, const int q, const int k,
                    const int j) {
        const int n = lf_only_sts_cell_update ? ((q == 0) ? IEN : IAN) : q;
        if (!UpdateSTSMHDVariable(n, update_momentum, update_energy,
                                  update_cgl_moment, update_scalars,
                                  nmhd_vars)) {
          return;
        }

        ScrArray1D<Real> divf(member.team_scratch(scr_level), ncells1);

        par_for_inner(member, is, ie, [&](const int i) {
          divf(i) = (flx1(m,n,k,j,i+1) - flx1(m,n,k,j,i))/mbsize.d_view(m).dx1;
        });
        member.team_barrier();

        if (multi_d) {
          par_for_inner(member, is, ie, [&](const int i) {
            divf(i) +=
                (flx2(m,n,k,j+1,i) - flx2(m,n,k,j,i))/mbsize.d_view(m).dx2;
          });
          member.team_barrier();
        }

        if (three_d) {
          par_for_inner(member, is, ie, [&](const int i) {
            divf(i) +=
                (flx3(m,n,k+1,j,i) - flx3(m,n,k,j,i))/mbsize.d_view(m).dx3;
          });
          member.team_barrier();
        }

        par_for_inner(member, is, ie, [&](const int i) {
          const bool weighted_cgl_variable =
              cgl_lf_weighted_flux && (n == IEN || n == IAN);
          if (weighted_cgl_variable) {
            // The CGL LF face flux already contains dt_sweep*muj_tilde.
            const Real weighted_rhs = -divf(i);
            const Real first_rhs_coeff =
                coeffs.gammaj_tilde/cgl_first_rkl_weight;
            u0_(m,n,k,j,i) = cgl_lf::WeightedRKL2Update(
                coeffs.muj, u_sts1_(m,n,k,j,i),
                coeffs.nuj, u_sts2_(m,n,k,j,i),
                1.0 - coeffs.muj - coeffs.nuj, u_sts0_(m,n,k,j,i),
                first_rhs_coeff, u_sts_rhs_(m,n,k,j,i), weighted_rhs);
            if (stage == 1) {
              u_sts_rhs_(m,n,k,j,i) = weighted_rhs;
            }
          } else {
            const Real delta_u = -dt_sweep*divf(i);
            u0_(m,n,k,j,i) = coeffs.muj*u_sts1_(m,n,k,j,i)
                           + coeffs.nuj*u_sts2_(m,n,k,j,i)
                           + (1.0 - coeffs.muj - coeffs.nuj)*u_sts0_(m,n,k,j,i)
                           + coeffs.gammaj_tilde*u_sts_rhs_(m,n,k,j,i)
                           + coeffs.muj_tilde*delta_u;
            if (stage == 1) {
              u_sts_rhs_(m,n,k,j,i) = delta_u;
            }
          }
        });
      });
    }
  }

  if (has_cgl_lf_split && pcgl_lf != nullptr) {
    pcgl_lf->AdvanceHeatFluxWorkDiagnostics(
        coeffs, stage, pdrive->sts.nstages);
  }
  if (has_cgl_lf_split && stage == pdrive->sts.nstages) {
    DiagnoseNonfiniteCGLState(stage, "sts-update-u", "post", STSSweepName(pdrive),
                              "magnetic-moment", u0);
  }

  TraceCGLLFTask(pmy_pack, "STSUpdateU", "end", stage);
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::STSUpdateB()
//! \brief Apply one CT-shaped RKL2 STS stage over the enrolled magnetic field.

TaskStatus MHD::STSUpdateB(Driver *pdrive, int stage) {
  if (!has_any_parabolic_field_update || !(pdrive->sts.enabled)) {
    return TaskStatus::complete;
  }

  if (stage == 1) {
    Kokkos::deep_copy(DevExeSpace(), b_sts0.x1f, b0.x1f);
    Kokkos::deep_copy(DevExeSpace(), b_sts0.x2f, b0.x2f);
    Kokkos::deep_copy(DevExeSpace(), b_sts0.x3f, b0.x3f);
  }
  Kokkos::deep_copy(DevExeSpace(), b_sts2.x1f, b_sts1.x1f);
  Kokkos::deep_copy(DevExeSpace(), b_sts2.x2f, b_sts1.x2f);
  Kokkos::deep_copy(DevExeSpace(), b_sts2.x3f, b_sts1.x3f);
  Kokkos::deep_copy(DevExeSpace(), b_sts1.x1f, b0.x1f);
  Kokkos::deep_copy(DevExeSpace(), b_sts1.x2f, b0.x2f);
  Kokkos::deep_copy(DevExeSpace(), b_sts1.x3f, b0.x3f);

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;

  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  Real dt_sweep = pdrive->sts.dt_sweep;
  auto coeffs = pdrive->sts.coeffs;
  auto e1 = efld.x1e;
  auto e2 = efld.x2e;
  auto e3 = efld.x3e;
  auto bx1f = b0.x1f;
  auto bx2f = b0.x2f;
  auto bx3f = b0.x3f;
  auto bx1f_sts0 = b_sts0.x1f;
  auto bx2f_sts0 = b_sts0.x2f;
  auto bx3f_sts0 = b_sts0.x3f;
  auto bx1f_sts1 = b_sts1.x1f;
  auto bx2f_sts1 = b_sts1.x2f;
  auto bx3f_sts1 = b_sts1.x3f;
  auto bx1f_sts2 = b_sts2.x1f;
  auto bx2f_sts2 = b_sts2.x2f;
  auto bx3f_sts2 = b_sts2.x3f;
  auto bx1f_rhs = b_sts_rhs.x1f;
  auto bx2f_rhs = b_sts_rhs.x2f;
  auto bx3f_rhs = b_sts_rhs.x3f;
  auto mbsize = pmy_pack->pmb->mb_size;

  par_for("mhd_sts_update_b1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real delta_b = 0.0;
    if (multi_d) {
      delta_b -= dt_sweep*(e3(m,k,j+1,i) - e3(m,k,j,i))/mbsize.d_view(m).dx2;
      if (three_d) {
        delta_b += dt_sweep*(e2(m,k+1,j,i) - e2(m,k,j,i))/mbsize.d_view(m).dx3;
      }
    }
    bx1f(m,k,j,i) = coeffs.muj*bx1f_sts1(m,k,j,i)
                  + coeffs.nuj*bx1f_sts2(m,k,j,i)
                  + (1.0 - coeffs.muj - coeffs.nuj)*bx1f_sts0(m,k,j,i)
                  + coeffs.gammaj_tilde*bx1f_rhs(m,k,j,i)
                  + coeffs.muj_tilde*delta_b;
    if (stage == 1) {
      bx1f_rhs(m,k,j,i) = delta_b;
    }
  });

  par_for("mhd_sts_update_b2", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real delta_b = dt_sweep*(e3(m,k,j,i+1) - e3(m,k,j,i))/mbsize.d_view(m).dx1;
    if (three_d) {
      delta_b -= dt_sweep*(e1(m,k+1,j,i) - e1(m,k,j,i))/mbsize.d_view(m).dx3;
    }
    bx2f(m,k,j,i) = coeffs.muj*bx2f_sts1(m,k,j,i)
                  + coeffs.nuj*bx2f_sts2(m,k,j,i)
                  + (1.0 - coeffs.muj - coeffs.nuj)*bx2f_sts0(m,k,j,i)
                  + coeffs.gammaj_tilde*bx2f_rhs(m,k,j,i)
                  + coeffs.muj_tilde*delta_b;
    if (stage == 1) {
      bx2f_rhs(m,k,j,i) = delta_b;
    }
  });

  par_for("mhd_sts_update_b3", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real delta_b = -dt_sweep*(e2(m,k,j,i+1) - e2(m,k,j,i))/mbsize.d_view(m).dx1;
    if (multi_d) {
      delta_b += dt_sweep*(e1(m,k,j+1,i) - e1(m,k,j,i))/mbsize.d_view(m).dx2;
    }
    bx3f(m,k,j,i) = coeffs.muj*bx3f_sts1(m,k,j,i)
                  + coeffs.nuj*bx3f_sts2(m,k,j,i)
                  + (1.0 - coeffs.muj - coeffs.nuj)*bx3f_sts0(m,k,j,i)
                  + coeffs.gammaj_tilde*bx3f_rhs(m,k,j,i)
                  + coeffs.muj_tilde*delta_b;
    if (stage == 1) {
      bx3f_rhs(m,k,j,i) = delta_b;
    }
  });

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::CGLLandauFluidPrimitiveRefresh()
//! \brief Recover CGL primitives between LF STS stages while IAN stores magnetic moment.

TaskStatus MHD::CGLLandauFluidPrimitiveRefresh(Driver *pdrive, int stage) {
  TraceCGLLFTask(pmy_pack, "CGLLandauFluidPrimitiveRefresh", "begin", stage);
  RequireCGLMagneticMomentRepresentation("CGL Landau-fluid primitive refresh");
  if (diagnose_nonfinite_rk_update) {
    DiagnoseNonfiniteCGLState(stage, "primitive-refresh", "pre",
                              STSSweepName(pdrive), "magnetic-moment", u0);
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int ng = indcs.ng;
  const int n1m1 = indcs.nx1 + 2*ng - 1;
  const int n2m1 = (indcs.nx2 > 1) ? indcs.nx2 + 2*ng - 1 : 0;
  const int n3m1 = (indcs.nx3 > 1) ? indcs.nx3 + 2*ng - 1 : 0;
  // Intermediate LF stages need only their one-cell stencil.  The final stage must
  // refresh every ghost consumed by the next hyperbolic reconstruction or collision.
  const bool final_stage = (stage == pdrive->sts.nstages);
  const int il = final_stage ? 0 : ((indcs.is > 0) ? indcs.is - 1 : 0);
  const int iu = final_stage ? n1m1 :
      ((indcs.ie < n1m1) ? indcs.ie + 1 : n1m1);
  const int jl = final_stage ? 0 :
      ((indcs.nx2 > 1 && indcs.js > 0) ? indcs.js - 1 : 0);
  const int ju = final_stage ? n2m1 :
      ((indcs.nx2 > 1 && indcs.je < n2m1) ? indcs.je + 1 : n2m1);
  const int kl = final_stage ? 0 :
      ((indcs.nx3 > 1 && indcs.ks > 0) ? indcs.ks - 1 : 0);
  const int ku = final_stage ? n3m1 :
      ((indcs.nx3 > 1 && indcs.ke < n3m1) ? indcs.ke + 1 : n3m1);
  const int dfloor_before = pmy_pack->pmesh->ecounter.neos_dfloor;
  const int pfloor_before = pmy_pack->pmesh->ecounter.neos_efloor;
  {
    CGLLFProfileRegion profile(pcgl_lf, CGLLFProfileBucket::primitive_refresh);
    if (pmy_pack->pmesh->multilevel) {
      RefreshCellCenteredBFromFace(pmy_pack, b0, bcc0, il, iu, jl, ju,
                                   kl, ku);
    }
    peos->CGLRefreshPrimFromMagneticMoment(u0, bcc0, w0, il, iu, jl, ju,
                                           kl, ku);
  }
  pcgl_lf->RecordAdmissibility(
      u0, w0, bcc0, peos->eos_data,
      pmy_pack->pmesh->ecounter.neos_dfloor - dfloor_before,
      pmy_pack->pmesh->ecounter.neos_efloor - pfloor_before,
      pdrive->sts.sweep == Driver::STSSweep::pre ? "pre" : "post",
      stage, pdrive->sts.nstages);
  if (diagnose_nonfinite_rk_update) {
    DiagnoseNonfiniteCGLState(stage, "primitive-refresh", "post",
                              STSSweepName(pdrive), "magnetic-moment", u0);
  }
  TraceCGLLFTask(pmy_pack, "CGLLandauFluidPrimitiveRefresh", "end", stage);
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::EndCGLLandauFluidSTSSweep()
//! \brief Restore conserved anisotropy after the final LF STS stage.

TaskStatus MHD::EndCGLLandauFluidSTSSweep(Driver *pdrive, int stage) {
  TraceCGLLFTask(pmy_pack, "EndCGLLandauFluidSTSSweep", "begin", stage);
  if (!has_cgl_lf_split || !pdrive->sts.enabled || stage != pdrive->sts.nstages) {
    TraceCGLLFTask(pmy_pack, "EndCGLLandauFluidSTSSweep", "skip", stage);
    return TaskStatus::complete;
  }
  RequireCGLMagneticMomentRepresentation("CGL Landau-fluid sweep end");
  if (diagnose_nonfinite_rk_update) {
    DiagnoseNonfiniteCGLState(stage, "magnetic-moment-to-anisotropy", "pre",
                              STSSweepName(pdrive), "magnetic-moment", u0);
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int ng = indcs.ng;
  const int n1m1 = indcs.nx1 + 2*ng - 1;
  const int n2m1 = (indcs.nx2 > 1) ? indcs.nx2 + 2*ng - 1 : 0;
  const int n3m1 = (indcs.nx3 > 1) ? indcs.nx3 + 2*ng - 1 : 0;
  {
    CGLLFProfileRegion profile(pcgl_lf, CGLLFProfileBucket::sweep_end_conversion);
    peos->CGLMagneticMomentToAnisotropy(u0, bcc0, 0, n1m1, 0, n2m1, 0, n3m1);
  }
  cgl_slot_representation = CGLSlotRepresentation::anisotropy;
  if (diagnose_nonfinite_rk_update) {
    DiagnoseNonfiniteCGLState(stage, "magnetic-moment-to-anisotropy", "post",
                              STSSweepName(pdrive), "anisotropy", u0);
  }
  TraceCGLLFTask(pmy_pack, "EndCGLLandauFluidSTSSweep", "end", stage);
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::STSPostSweepCGLCollisions()
//! \brief Apply CGL relaxation after each split LF half-sweep.
//!
//! The pre and post LF sweeps each advance dt_cycle/2. Use that same physical
//! interval here so both source updates together advance exactly one cycle.

TaskStatus MHD::STSPostSweepCGLCollisions(Driver *pdrive, int stage) {
  TraceCGLLFTask(pmy_pack, "STSPostSweepCGLCollisions", "begin", stage);
  if (!has_cgl_lf_split || !peos->eos_data.coll || stage != pdrive->sts.nstages) {
    TraceCGLLFTask(pmy_pack, "STSPostSweepCGLCollisions", "skip", stage);
    return TaskStatus::complete;
  }
  RequireCGLAnisotropyRepresentation("CGL Landau-fluid post-sweep collisions");
  if (diagnose_nonfinite_rk_update) {
    DiagnoseNonfiniteCGLState(stage, "post-sweep-collisions", "pre",
                              STSSweepName(pdrive), "anisotropy", u0);
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int ng = indcs.ng;
  const int n1m1 = indcs.nx1 + 2*ng - 1;
  const int n2m1 = (indcs.nx2 > 1) ? indcs.nx2 + 2*ng - 1 : 0;
  const int n3m1 = (indcs.nx3 > 1) ? indcs.nx3 + 2*ng - 1 : 0;
  {
    CGLLFProfileRegion profile(pcgl_lf, CGLLFProfileBucket::post_sweep_collisions);
    peos->Collisions(w0, bcc0, u0, pdrive->sts.dt_sweep,
                     0, n1m1, 0, n2m1, 0, n3m1);
  }
  if (diagnose_nonfinite_rk_update) {
    DiagnoseNonfiniteCGLState(stage, "post-sweep-collisions", "post",
                              STSSweepName(pdrive), "anisotropy", u0);
  }
  TraceCGLLFTask(pmy_pack, "STSPostSweepCGLCollisions", "end", stage);
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::STSRefreshTimeStep()
//! \brief Refresh MHD-local timestep estimates after the final post sweep stage.

TaskStatus MHD::STSRefreshTimeStep(Driver *pdrive, int stage) {
  TraceCGLLFTask(pmy_pack, "STSRefreshTimeStep", "begin", stage);
  if (!has_any_parabolic_split || !(pdrive->sts.enabled)) {
    TraceCGLLFTask(pmy_pack, "STSRefreshTimeStep", "skip", stage);
    return TaskStatus::complete;
  }
  if (pdrive->sts.sweep == Driver::STSSweep::post && stage == pdrive->sts.nstages) {
    RecomputeTimeStepFromCurrentState(pdrive);
  }
  TraceCGLLFTask(pmy_pack, "STSRefreshTimeStep", "end", stage);
  return TaskStatus::complete;
}

} // namespace mhd
