//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hyperviscous_shear.cpp
//! \brief Analytic verification problem for fourth-derivative velocity damping.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "parameter_input.hpp"
#include "pgen/pgen.hpp"

namespace {

bool create_reference_state = false;

void HyperViscousShearErrors(ParameterInput *pin, Mesh *pm);

void FatalHyperViscousShear(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl << message << std::endl;
  std::exit(EXIT_FAILURE);
}

} // namespace

//----------------------------------------------------------------------------------------
//! \brief Initialize a constant-density transverse Fourier mode with exact damping.

void ProblemGenerator::HyperViscousShear(ParameterInput *pin, const bool restart) {
  if (restart) return;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  const bool is_hydro = (pmbp->phydro != nullptr);
  const bool is_mhd = (pmbp->pmhd != nullptr);
  if (is_hydro == is_mhd) {
    FatalHyperViscousShear(
        "hyperviscous_shear requires exactly one of <hydro> or <mhd>.");
  }
  const std::string block = is_hydro ? "hydro" : "mhd";
  if (!pin->DoesParameterExist(block, "hyperviscosity")) {
    FatalHyperViscousShear(
        "hyperviscous_shear requires a hyperviscosity coefficient.");
  }
  if (pin->GetOrAddBoolean(block, "tdep_viscosity", false)) {
    FatalHyperViscousShear(
        "hyperviscous_shear analytic reference does not support tdep_viscosity.");
  }
  pgen_final_func = HyperViscousShearErrors;

  auto &indcs = pmy_mesh_->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx1 = indcs.nx1;
  const int nmb = pmbp->nmb_thispack;
  auto size = pmbp->pmb->mb_size;
  const Real lx = pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min;
  const Real x1min = pmy_mesh_->mesh_size.x1min;
  const Real amp = pin->GetOrAddReal("problem", "amp", 0.1);
  const Real mode = pin->GetOrAddReal("problem", "mode", 2.0);
  const Real rho = pin->GetOrAddReal("problem", "density", 1.0);
  const Real pressure = pin->GetOrAddReal("problem", "pressure", 1.0);
  const Real kx = 2.0*M_PI*mode/lx;
  const Real nu4 = pin->GetReal(block, "hyperviscosity");
  const Real nu2 = pin->DoesParameterExist(block, "viscosity") ?
      pin->GetReal(block, "viscosity") : 0.0;
  const Real time = create_reference_state ? pmy_mesh_->time : 0.0;
  const Real damping = std::exp(-(nu2*kx*kx + nu4*kx*kx*kx*kx)*time);

  if (is_hydro) {
    EOS_Data &eos = pmbp->phydro->peos->eos_data;
    auto u = create_reference_state ? pmbp->phydro->u1 : pmbp->phydro->u0;
    const Real gm1 = eos.gamma - 1.0;
    par_for("hypervisc_shear_hydro", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      const Real x = CellCenterX(i-is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
      const Real vy = amp*damping*std::sin(kx*(x - x1min));
      u(m,IDN,k,j,i) = rho;
      u(m,IM1,k,j,i) = 0.0;
      u(m,IM2,k,j,i) = rho*vy;
      u(m,IM3,k,j,i) = 0.0;
      if (eos.is_ideal) {
        u(m,IEN,k,j,i) = pressure/gm1 + 0.5*rho*vy*vy;
      }
    });
  } else {
    EOS_Data &eos = pmbp->pmhd->peos->eos_data;
    auto u = create_reference_state ? pmbp->pmhd->u1 : pmbp->pmhd->u0;
    auto &b0 = pmbp->pmhd->b0;
    auto &bcc0 = pmbp->pmhd->bcc0;
    const Real gm1 = eos.gamma - 1.0;
    par_for("hypervisc_shear_mhd", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      const Real x = CellCenterX(i-is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
      const Real vy = amp*damping*std::sin(kx*(x - x1min));
      u(m,IDN,k,j,i) = rho;
      u(m,IM1,k,j,i) = 0.0;
      u(m,IM2,k,j,i) = rho*vy;
      u(m,IM3,k,j,i) = 0.0;
      if (eos.is_ideal) {
        u(m,IEN,k,j,i) = pressure/gm1 + 0.5*rho*vy*vy;
      }
      if (!create_reference_state) {
        b0.x1f(m,k,j,i) = 0.0;
        b0.x2f(m,k,j,i) = 0.0;
        b0.x3f(m,k,j,i) = 0.0;
        if (i == ie) b0.x1f(m,k,j,i+1) = 0.0;
        if (j == je) b0.x2f(m,k,j+1,i) = 0.0;
        if (k == ke) b0.x3f(m,k+1,j,i) = 0.0;
        bcc0(m,IBX,k,j,i) = 0.0;
        bcc0(m,IBY,k,j,i) = 0.0;
        bcc0(m,IBZ,k,j,i) = 0.0;
      }
    });
  }
}

//----------------------------------------------------------------------------------------
//! \brief Compare the evolved transverse momentum with the exact Fourier-mode solution.

namespace {

void HyperViscousShearErrors(ParameterInput *pin, Mesh *pm) {
  create_reference_state = true;
  pm->pgen->HyperViscousShear(pin, false);
  create_reference_state = false;

  MeshBlockPack *pmbp = pm->pmb_pack;
  const bool is_hydro = (pmbp->phydro != nullptr);
  auto u0 = is_hydro ? pmbp->phydro->u0 : pmbp->pmhd->u0;
  auto u1 = is_hydro ? pmbp->phydro->u1 : pmbp->pmhd->u1;
  auto size = pmbp->pmb->mb_size;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nmkji = pmbp->nmb_thispack*nx1*nx2*nx3;
  const int nkji = nx1*nx2*nx3;
  const int nji = nx1*nx2;
  const Real x1min = pm->mesh_size.x1min;
  const Real lx = pm->mesh_size.x1max - pm->mesh_size.x1min;
  const Real mode = pin->GetOrAddReal("problem", "mode", 2.0);
  const Real kx = 2.0*M_PI*mode/lx;
  Real l1 = 0.0, linf = 0.0, projected_amp = 0.0;

  Kokkos::parallel_reduce("hypervisc_shear_l1", Kokkos::RangePolicy<>(DevExeSpace(), 0,
      nmkji), KOKKOS_LAMBDA(const int idx, Real &sum) {
    const int m = idx/nkji;
    const int rem_m = idx - m*nkji;
    const int k0 = rem_m/nji;
    const int rem_k = rem_m - k0*nji;
    const int j0 = rem_k/nx1;
    const int i = rem_k - j0*nx1 + is;
    const int j = j0 + js;
    const int k = k0 + ks;
    const Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
    sum += vol*std::fabs(u0(m,IM2,k,j,i) - u1(m,IM2,k,j,i));
  }, Kokkos::Sum<Real>(l1));
  Kokkos::parallel_reduce("hypervisc_shear_linf", Kokkos::RangePolicy<>(DevExeSpace(), 0,
      nmkji), KOKKOS_LAMBDA(const int idx, Real &max_error) {
    const int m = idx/nkji;
    const int rem_m = idx - m*nkji;
    const int k0 = rem_m/nji;
    const int rem_k = rem_m - k0*nji;
    const int j0 = rem_k/nx1;
    const int i = rem_k - j0*nx1 + is;
    const int j = j0 + js;
    const int k = k0 + ks;
    max_error = fmax(max_error, std::fabs(u0(m,IM2,k,j,i) - u1(m,IM2,k,j,i)));
  }, Kokkos::Max<Real>(linf));
  Kokkos::parallel_reduce("hypervisc_shear_amp", Kokkos::RangePolicy<>(DevExeSpace(), 0,
      nmkji), KOKKOS_LAMBDA(const int idx, Real &sum) {
    const int m = idx/nkji;
    const int rem_m = idx - m*nkji;
    const int k0 = rem_m/nji;
    const int rem_k = rem_m - k0*nji;
    const int j0 = rem_k/nx1;
    const int i = rem_k - j0*nx1 + is;
    const int j = j0 + js;
    const int k = k0 + ks;
    const Real x = CellCenterX(i-is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    const Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
    sum += 2.0*vol*u0(m,IM2,k,j,i)*std::sin(kx*(x - x1min));
  }, Kokkos::Sum<Real>(projected_amp));

#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &l1, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &linf, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &projected_amp, 1, MPI_ATHENA_REAL, MPI_SUM,
                MPI_COMM_WORLD);
#endif

  const Real volume = (pm->mesh_size.x1max - pm->mesh_size.x1min)*
      (pm->mesh_size.x2max - pm->mesh_size.x2min)*
      (pm->mesh_size.x3max - pm->mesh_size.x3min);
  l1 /= volume;
  projected_amp /= volume;

  if (global_variable::my_rank == 0) {
    const std::string fname = pin->GetString("job", "basename") + "-hypervisc-errors.dat";
    FILE *pfile = std::fopen(fname.c_str(), "w");
    if (pfile == nullptr) {
      FatalHyperViscousShear("Could not write hyperviscosity error output.");
    }
    std::fprintf(pfile, "# Nx1 Ncycle L1 Linf ProjectedAmplitude\n");
    std::fprintf(pfile, "%d %d %.16e %.16e %.16e\n",
                 pm->mesh_indcs.nx1, pm->ncycle, l1, linf, projected_amp);
    std::fclose(pfile);
  }
}

} // namespace
