//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file turb.cpp
//  \brief Problem generator for turbulence
#include <cmath>
#include <iostream> // cout

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "particles/particles.hpp"
#include "pgen.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

// User-defined history functions
void TurbulentHistory(HistoryData *pdata, Mesh *pm);


//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::Turb_()
//  \brief Problem Generator for turbulence

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  // Re-enroll the callback on restart without reinitializing the evolved state.
  user_hist_func = TurbulentHistory;
  if (restart) return;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;

  if (pmbp->phydro == nullptr && pmbp->pmhd == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
       << "Turbulence problem generator can only be run with Hydro and/or MHD, but no "
       << "<hydro> or <mhd> block in input file" << std::endl;
    exit(EXIT_FAILURE);
  }

  // capture variables for kernel
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;

  Real cs = pin->GetOrAddReal("eos","iso_sound_speed",1.0);
  Real beta = pin->GetOrAddReal("problem","beta",1.0);

  // Initialize Hydro variables -------------------------------
  if (pmbp->phydro != nullptr) {
    Real d_i = pin->GetOrAddReal("problem","d_i",1.0);
    Real d_n = pin->GetOrAddReal("problem","d_n",1.0);
    auto &u0 = pmbp->phydro->u0;
    EOS_Data &eos = pmbp->phydro->peos->eos_data;
    Real gm1 = eos.gamma - 1.0;
    Real p0 = 1.0/eos.gamma;

    // Set initial conditions
    par_for("pgen_turb", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0(m,IDN,k,j,i) = d_n;
      u0(m,IM1,k,j,i) = 0.0;
      u0(m,IM2,k,j,i) = 0.0;
      u0(m,IM3,k,j,i) = 0.0;
      if (eos.is_ideal) {
        u0(m,IEN,k,j,i) = p0/gm1 +
           0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i)) +
           SQR(u0(m,IM3,k,j,i)))/u0(m,IDN,k,j,i);
      }
    });
  }

  // Initialize MHD variables ---------------------------------
  if (pmbp->pmhd != nullptr) {
    Real d_i = pin->GetOrAddReal("problem","d_i",1.0);
    Real d_n = pin->GetOrAddReal("problem","d_n",1.0);
    int ifield = pin->GetOrAddInteger("problem","ifield",2);
    if (ifield < 0 || ifield > 2) {
      std::cout << "### FATAL ERROR in " << __FILE__
                << " at line " << __LINE__ << std::endl
                << "Invalid <problem>/ifield = " << ifield
                << ", allowed values are 0 (zero field), 1 (zero-net-flux Bz), "
                << "or 2 (uniform Bz)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if (ifield != 0 && beta <= 0.0) {
      std::cout << "### FATAL ERROR in " << __FILE__
                << " at line " << __LINE__ << std::endl
                << "<problem>/beta must be positive when <problem>/ifield is 1 or 2."
                << std::endl;
      exit(EXIT_FAILURE);
    }
    Real B0 = (ifield == 0) ? 0.0 : cs*std::sqrt(2.0*d_i/beta);
    Real x1size = pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min;
    Real kx = 2.0*(M_PI/x1size);
    auto &u0 = pmbp->pmhd->u0;
    auto &b0 = pmbp->pmhd->b0;
    auto &size = pmbp->pmb->mb_size;
    EOS_Data &eos = pmbp->pmhd->peos->eos_data;
    Real gm1 = 0.0;
    Real p0 = 0.0;
    if (eos.is_ideal) {
      gm1 = eos.gamma - 1.0;
      p0 = d_i*SQR(cs)/eos.gamma;
      B0 = (ifield == 0) ? 0.0 : std::sqrt(2.0*p0/beta);
    }

    // Set initial conditions
    par_for("pgen_turb", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0(m,IDN,k,j,i) = d_i;
      u0(m,IM1,k,j,i) = 0.0;
      u0(m,IM2,k,j,i) = 0.0;
      u0(m,IM3,k,j,i) = 0.0;

      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      int nx1 = indcs.nx1;
      Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

      Real bz = 0.0;
      if (ifield == 1) {
        bz = B0*std::sin(kx*x1v);
      } else if (ifield == 2) {
        bz = B0;
      }
      b0.x1f(m,k,j,i) = 0.0;
      b0.x2f(m,k,j,i) = 0.0;
      b0.x3f(m,k,j,i) = bz;
      if (i==ie) {b0.x1f(m,k,j,i+1) = 0.0;}
      if (j==je) {b0.x2f(m,k,j+1,i) = 0.0;}
      if (k==ke) {b0.x3f(m,k+1,j,i) = bz;}

      if (eos.is_ideal) {
        Real bz_cc = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
        u0(m,IEN,k,j,i) = p0/gm1 + 0.5*bz_cc*bz_cc +
           0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i)) +
           SQR(u0(m,IM3,k,j,i)))/u0(m,IDN,k,j,i);
      }
    });
  }

  // Initialize ion-neutral variables -------------------------
  if (pmbp->pionn != nullptr) {
    Real d_i = pin->GetOrAddReal("problem","d_i",1.0);
    Real d_n = pin->GetOrAddReal("problem","d_n",1.0);
    int ifield = pin->GetOrAddInteger("problem","ifield",2);
    if (ifield < 0 || ifield > 2) {
      std::cout << "### FATAL ERROR in " << __FILE__
                << " at line " << __LINE__ << std::endl
                << "Invalid <problem>/ifield = " << ifield
                << ", allowed values are 0 (zero field), 1 (zero-net-flux Bz), "
                << "or 2 (uniform Bz)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if (ifield != 0 && beta <= 0.0) {
      std::cout << "### FATAL ERROR in " << __FILE__
                << " at line " << __LINE__ << std::endl
                << "<problem>/beta must be positive when <problem>/ifield is 1 or 2."
                << std::endl;
      exit(EXIT_FAILURE);
    }
    Real B0 = (ifield == 0) ? 0.0 : cs*std::sqrt(2.0*(d_i+d_n)/beta);
    Real x1size = pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min;
    Real kx = 2.0*(M_PI/x1size);

    // MHD
    auto &u0 = pmbp->pmhd->u0;
    auto &b0 = pmbp->pmhd->b0;
    auto &size = pmbp->pmb->mb_size;
    EOS_Data &eos = pmbp->pmhd->peos->eos_data;
    Real gm1 = eos.gamma - 1.0;
    Real p0 = d_i/eos.gamma; // TODO(@user): multiply by ionized density

    // Set initial conditions
    par_for("pgen_turb_mhd", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0(m,IDN,k,j,i) = d_i;
      u0(m,IM1,k,j,i) = 0.0;
      u0(m,IM2,k,j,i) = 0.0;
      u0(m,IM3,k,j,i) = 0.0;

      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      int nx1 = indcs.nx1;
      Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

      Real bz = 0.0;
      if (ifield == 1) {
        bz = B0*std::sin(kx*x1v);
      } else if (ifield == 2) {
        bz = B0;
      }
      b0.x1f(m,k,j,i) = 0.0;
      b0.x2f(m,k,j,i) = 0.0;
      b0.x3f(m,k,j,i) = bz;
      if (i==ie) {b0.x1f(m,k,j,i+1) = 0.0;}
      if (j==je) {b0.x2f(m,k,j+1,i) = 0.0;}
      if (k==ke) {b0.x3f(m,k+1,j,i) = bz;}

      if (eos.is_ideal) {
        Real bz_cc = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
        u0(m,IEN,k,j,i) = p0/gm1 + 0.5*bz_cc*bz_cc +
           0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i)) +
           SQR(u0(m,IM3,k,j,i)))/u0(m,IDN,k,j,i);
      }
    });
    // Hydro
    auto &u0_ = pmbp->phydro->u0;
    EOS_Data &eos_ = pmbp->phydro->peos->eos_data;
    Real gm1_ = eos_.gamma - 1.0;
    Real p0_ = d_n/eos_.gamma; // TODO(@user): multiply by neutral density

    // Set initial conditions
    par_for("pgen_turb_hydro", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0_(m,IDN,k,j,i) = d_n;
      u0_(m,IM1,k,j,i) = 0.0;
      u0_(m,IM2,k,j,i) = 0.0;
      u0_(m,IM3,k,j,i) = 0.0;
      if (eos_.is_ideal) {
        u0_(m,IEN,k,j,i) = p0_/gm1_ +
            0.5*(SQR(u0_(m,IM1,k,j,i)) + SQR(u0_(m,IM2,k,j,i)) +
            SQR(u0_(m,IM3,k,j,i)))/u0_(m,IDN,k,j,i);
      }
    });
  }

  return;
}


//----------------------------------------------------------------------------------------
// Function for computing dynamo history variables
// The built-in MHD history already records gas mass, momentum, total energy,
// component kinetic energies, and component magnetic energies.  These additional
// reductions supply the volume-weighted velocity/field means, velocity variance,
// discrete div(B), and CR energy/momentum needed by the turbulent-dynamo analysis.
void TurbulentHistory(HistoryData *pdata, Mesh *pm) {
  if (pm == nullptr || pm->pmb_pack == nullptr || pm->pmb_pack->pmhd == nullptr) {
    pdata->nhist = 0;
    return;
  }

  pdata->nhist = 12;
  pdata->label[0] = "vx_vol";
  pdata->label[1] = "vy_vol";
  pdata->label[2] = "vz_vol";
  pdata->label[3] = "v2_vol";
  pdata->label[4] = "Bx_vol";
  pdata->label[5] = "By_vol";
  pdata->label[6] = "Bz_vol";
  pdata->label[7] = "divB_max";
  pdata->label[8] = "cr_Ekin";
  pdata->label[9] = "cr_Px";
  pdata->label[10] = "cr_Py";
  pdata->label[11] = "cr_Pz";

  auto &bcc = pm->pmb_pack->pmhd->bcc0;
  auto &b = pm->pmb_pack->pmhd->b0;
  auto &w0 = pm->pmb_pack->pmhd->w0;
  auto &size = pm->pmb_pack->pmb->mb_size;

  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  const int is = indcs.is;
  const int js = indcs.js;
  const int ks = indcs.ks;
  const int nx1 = indcs.nx1;
  const int nx2 = indcs.nx2;
  const int nx3 = indcs.nx3;
  const int nmkji = pm->pmb_pack->nmb_thispack*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji = nx2*nx1;
  const bool multi_d = pm->multi_d;
  const bool three_d = pm->three_d;

  array_sum::GlobalSum fluid_sums;
  Kokkos::parallel_reduce(
      "TurbulentDynamoHistory",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
      KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &sum) {
        const int m = idx/nkji;
        const int k0 = (idx - m*nkji)/nji;
        const int j0 = (idx - m*nkji - k0*nji)/nx1;
        const int i = idx - m*nkji - k0*nji - j0*nx1 + is;
        const int j = j0 + js;
        const int k = k0 + ks;
        const Real vol =
            size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
        const Real vx = w0(m,IVX,k,j,i);
        const Real vy = w0(m,IVY,k,j,i);
        const Real vz = w0(m,IVZ,k,j,i);

        array_sum::GlobalSum cell;
        cell.the_array[0] = vol*vx;
        cell.the_array[1] = vol*vy;
        cell.the_array[2] = vol*vz;
        cell.the_array[3] = vol*(vx*vx + vy*vy + vz*vz);
        cell.the_array[4] = vol*bcc(m,IBX,k,j,i);
        cell.the_array[5] = vol*bcc(m,IBY,k,j,i);
        cell.the_array[6] = vol*bcc(m,IBZ,k,j,i);
        for (int n = 7; n < NREDUCTION_VARIABLES; ++n) {
          cell.the_array[n] = 0.0;
        }
        sum += cell;
      },
      Kokkos::Sum<array_sum::GlobalSum>(fluid_sums));

  for (int n = 0; n < 7; ++n) {
    pdata->hdata[n] = fluid_sums.the_array[n];
  }

  Real max_divb = 0.0;
  Kokkos::parallel_reduce(
      "TurbulentDynamoMaxDivB",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
      KOKKOS_LAMBDA(const int idx, Real &max_value) {
        const int m = idx/nkji;
        const int k0 = (idx - m*nkji)/nji;
        const int j0 = (idx - m*nkji - k0*nji)/nx1;
        const int i = idx - m*nkji - k0*nji - j0*nx1 + is;
        const int j = j0 + js;
        const int k = k0 + ks;
        Real divb =
            (b.x1f(m,k,j,i+1) - b.x1f(m,k,j,i))/size.d_view(m).dx1;
        if (multi_d) {
          divb +=
              (b.x2f(m,k,j+1,i) - b.x2f(m,k,j,i))/size.d_view(m).dx2;
        }
        if (three_d) {
          divb +=
              (b.x3f(m,k+1,j,i) - b.x3f(m,k,j,i))/size.d_view(m).dx3;
        }
        max_value = fmax(max_value, fabs(divb));
      },
      Kokkos::Max<Real>(max_divb));

#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &max_divb, 1, MPI_ATHENA_REAL, MPI_MAX,
                MPI_COMM_WORLD);
  // HistoryOutput performs a sum reduction after this callback.  Give each rank
  // an equal share of the already-global maximum so the written value remains a max.
  max_divb /= static_cast<Real>(global_variable::nranks);
#endif
  pdata->hdata[7] = max_divb;
  for (int n = 8; n < pdata->nhist; ++n) {
    pdata->hdata[n] = 0.0;
  }

  auto *ppart = pm->pmb_pack->ppart;
  if (ppart != nullptr && ppart->particle_type == ParticleType::cosmic_ray &&
      ppart->nprtcl_thispack > 0) {
    auto &pr = ppart->prtcl_rdata;
    auto &pi = ppart->prtcl_idata;
    auto masses = ppart->species_mass;
    const int nspecies = ppart->nspecies;
    const Real qscale = ppart->deposit_qscale;
    const bool momentum_state = ppart->UsesRelativisticCRState();
    const Real light_speed = ppart->pic_cr_light_speed;

    Real cr_energy = 0.0;
    Real cr_px = 0.0;
    Real cr_py = 0.0;
    Real cr_pz = 0.0;
    Kokkos::parallel_reduce(
        "TurbulentDynamoCRHistory",
        Kokkos::RangePolicy<>(DevExeSpace(), 0, ppart->nprtcl_thispack),
        KOKKOS_LAMBDA(const int p, Real &energy, Real &px, Real &py, Real &pz) {
          const int species = pi(PSP,p);
          if (species < 0 || species >= nspecies) return;
          Real weight = pr(IPWT,p);
          if (weight <= 0.0) weight = 1.0;
          const Real macro_mass = qscale*weight*masses(species);
          energy += macro_mass*particles::CRKineticEnergy(
              momentum_state, light_speed,
              pr(IPVX,p), pr(IPVY,p), pr(IPVZ,p));
          px += macro_mass*pr(IPVX,p);
          py += macro_mass*pr(IPVY,p);
          pz += macro_mass*pr(IPVZ,p);
        },
        Kokkos::Sum<Real>(cr_energy), Kokkos::Sum<Real>(cr_px),
        Kokkos::Sum<Real>(cr_py), Kokkos::Sum<Real>(cr_pz));
    pdata->hdata[8] = cr_energy;
    pdata->hdata[9] = cr_px;
    pdata->hdata[10] = cr_py;
    pdata->hdata[11] = cr_pz;
  }
}
