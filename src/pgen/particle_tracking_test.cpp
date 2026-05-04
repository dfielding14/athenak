//========================================================================================
//! \file particle_tracking_test.cpp
//! \brief Deterministic pgen for validating lagrangian_mc particle tracking output.
//========================================================================================

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "parameter_input.hpp"
#include "particles/particles.hpp"
#include "pgen.hpp"

namespace {

std::string test_mode;
std::string field_profile;
std::string particle_layout;
Real rho0;
Real temp0;
Real scalar0;
Real vx0;
Real vy0;
Real vz0;
Real t_refine;
Real t_derefine;

struct ParticleSeed {
  int gid;
  int level;
  Real x;
  Real y;
  Real z;
};

KOKKOS_INLINE_FUNCTION
Real AnalyticRho(const int i, const int j, const int k, const int level) {
  return 10.0 + static_cast<Real>(i) + 10.0*static_cast<Real>(j) +
         100.0*static_cast<Real>(k) + 1000.0*static_cast<Real>(level);
}

KOKKOS_INLINE_FUNCTION
Real AnalyticTemp(const int i, const int j, const int k) {
  return 1.0 + 0.01*static_cast<Real>(i) + 0.1*static_cast<Real>(j) +
         static_cast<Real>(k);
}

KOKKOS_INLINE_FUNCTION
Real AnalyticScalar(const Real x, const Real y, const Real z) {
  return x + 2.0*y + 3.0*z;
}

void SetHydroFields(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp->phydro == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "particle_tracking_test requires a <hydro> block." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  auto &indcs = pm->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const int nx1 = indcs.nx1;
  const int nx2 = indcs.nx2;
  const int nx3 = indcs.nx3;
  const int nmb = pmbp->nmb_thispack;

  auto &size = pmbp->pmb->mb_size;
  auto &mblev = pmbp->pmb->mb_lev;
  auto &w0 = pmbp->phydro->w0;
  auto &u0 = pmbp->phydro->u0;
  const int nhydro = pmbp->phydro->nhydro;
  const int nscalars = pmbp->phydro->nscalars;
  const Real gamma = pmbp->phydro->peos->eos_data.gamma;
  const Real gm1 = gamma - 1.0;
  const bool analytic = (field_profile == "analytic");
  const Real rho_const = rho0;
  const Real temp_const = temp0;
  const Real scalar_const = scalar0;
  const Real vx = vx0;
  const Real vy = vy0;
  const Real vz = vz0;

  par_for("particle_tracking_fields", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int il = i - is;
    const int jl = j - js;
    const int kl = k - ks;
    const int level = mblev.d_view(m);
    const Real x = CellCenterX(il, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    const Real y = CellCenterX(jl, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    const Real z = CellCenterX(kl, nx3, size.d_view(m).x3min, size.d_view(m).x3max);

    const Real rho = analytic ? AnalyticRho(il, jl, kl, level) : rho_const;
    const Real temp = analytic ? AnalyticTemp(il, jl, kl) : temp_const;
    const Real scalar = analytic ? AnalyticScalar(x, y, z) : scalar_const;
    const Real press = rho*temp;
    const Real eint = press/gm1;
    const Real ekin = 0.5*rho*(vx*vx + vy*vy + vz*vz);

    w0(m,IDN,k,j,i) = rho;
    w0(m,IVX,k,j,i) = vx;
    w0(m,IVY,k,j,i) = vy;
    w0(m,IVZ,k,j,i) = vz;
    w0(m,IEN,k,j,i) = eint;

    u0(m,IDN,k,j,i) = rho;
    u0(m,IM1,k,j,i) = rho*vx;
    u0(m,IM2,k,j,i) = rho*vy;
    u0(m,IM3,k,j,i) = rho*vz;
    u0(m,IEN,k,j,i) = eint + ekin;

    if (nscalars > 0) {
      w0(m,nhydro,k,j,i) = scalar;
      u0(m,nhydro,k,j,i) = rho*scalar;
    }
  });
  Kokkos::fence();
}

bool IncludeParticleCell(const std::string &layout, const int i, const int is, const int ie,
                         const int j, const int js, const int k, const int ks) {
  if (layout == "none") {
    return false;
  }
  if (layout == "cell_centers") {
    return true;
  }
  if (layout == "outer_x1") {
    return i == ie;
  }
  if (layout == "inner_x1") {
    return i == is;
  }
  if (layout == "line_x1") {
    return j == js && k == ks;
  }
  return false;
}

void InitializeParticles(ParameterInput *pin, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp->ppart == nullptr) {
    return;
  }
  if (pmbp->ppart->particle_type != ParticleType::lagrangian_mc) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "particle_tracking_test requires lagrangian_mc particles."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  auto &indcs = pm->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const int nx1 = indcs.nx1;
  const int nx2 = indcs.nx2;
  const int nx3 = indcs.nx3;
  const int max_initial =
      pin->GetOrAddInteger("problem", "max_initial_particles", -1);

  pmbp->pmb->mb_size.template sync<HostMemSpace>();
  pmbp->pmb->mb_lev.template sync<HostMemSpace>();

  std::vector<ParticleSeed> seeds;
  for (int m = 0; m < pmbp->nmb_thispack; ++m) {
    const auto &mb = pmbp->pmb->mb_size.h_view(m);
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie; ++i) {
          if (!IncludeParticleCell(particle_layout, i, is, ie, j, js, k, ks)) {
            continue;
          }
          if (max_initial >= 0 && static_cast<int>(seeds.size()) >= max_initial) {
            break;
          }
          seeds.push_back({pmbp->gids + m,
                           pmbp->pmb->mb_lev.h_view(m),
                           CellCenterX(i - is, nx1, mb.x1min, mb.x1max),
                           CellCenterX(j - js, nx2, mb.x2min, mb.x2max),
                           CellCenterX(k - ks, nx3, mb.x3min, mb.x3max)});
        }
      }
    }
  }

  const int npart = static_cast<int>(seeds.size());
  pmbp->ppart->ReallocateParticles(npart);
  HostArray2D<Real> h_pr("particle_tracking_pr", pmbp->ppart->nrdata, npart);
  HostArray2D<int> h_pi("particle_tracking_pi", pmbp->ppart->nidata, npart);
  for (int p = 0; p < npart; ++p) {
    h_pi(PGID,p) = seeds[p].gid;
    h_pi(PTAG,p) = p;
    h_pi(PLASTMOVE,p) = 0;
    h_pi(PLASTLEVEL,p) = seeds[p].level;
    h_pr(LMCX,p) = seeds[p].x;
    h_pr(LMCY,p) = seeds[p].y;
    h_pr(LMCZ,p) = seeds[p].z;
  }
  Kokkos::deep_copy(pmbp->ppart->prtcl_rdata, h_pr);
  Kokkos::deep_copy(pmbp->ppart->prtcl_idata, h_pi);

  pm->CountParticles();
  pmbp->ppart->CreateParticleTags(pin);
}

void ParticleTrackingSource(Mesh *pm, const Real /*bdt*/) {
  SetHydroFields(pm);
}

void ParticleTrackingWorkInLoop(Mesh *pm) {
  SetHydroFields(pm);
}

void ParticleTrackingRefine(MeshBlockPack *pmbp) {
  Mesh *pm = pmbp->pmesh;
  auto &refine_flag = pm->pmr->refine_flag;
  const int mbs = pm->gids_eachrank[global_variable::my_rank];
  int flag = 0;
  if (t_refine >= 0.0 && pm->time >= t_refine &&
      (t_derefine < 0.0 || pm->time < t_derefine)) {
    flag = 1;
  } else if (t_derefine >= 0.0 && pm->time >= t_derefine) {
    flag = -1;
  }
  for (int m = 0; m < pmbp->nmb_thispack; ++m) {
    refine_flag.h_view(m + mbs) = flag;
  }
  refine_flag.template modify<HostMemSpace>();
  refine_flag.template sync<DevExeSpace>();
}

} // namespace

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  test_mode = pin->GetOrAddString("problem", "test_mode", "static");
  field_profile = pin->GetOrAddString("problem", "field_profile", "analytic");
  particle_layout = pin->GetOrAddString("problem", "particle_layout", "cell_centers");
  rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  temp0 = pin->GetOrAddReal("problem", "temp0", 1.0);
  scalar0 = pin->GetOrAddReal("problem", "scalar0", 0.0);
  vx0 = pin->GetOrAddReal("problem", "vx0", 0.0);
  vy0 = pin->GetOrAddReal("problem", "vy0", 0.0);
  vz0 = pin->GetOrAddReal("problem", "vz0", 0.0);
  t_refine = pin->GetOrAddReal("problem", "t_refine", -1.0);
  t_derefine = pin->GetOrAddReal("problem", "t_derefine", -1.0);

  const bool valid_profile = (field_profile == "analytic" || field_profile == "uniform");
  const bool valid_layout = (particle_layout == "none" ||
                             particle_layout == "cell_centers" ||
                             particle_layout == "outer_x1" ||
                             particle_layout == "inner_x1" ||
                             particle_layout == "line_x1");
  if (!valid_profile || !valid_layout) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Invalid particle_tracking_test profile/layout: field_profile='"
              << field_profile << "', particle_layout='" << particle_layout << "'."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  user_srcs_func = ParticleTrackingSource;
  user_work_in_loop_func = ParticleTrackingWorkInLoop;
  if (pmy_mesh_->adaptive) {
    user_ref_func = ParticleTrackingRefine;
  }

  if (restart) {
    return;
  }

  SetHydroFields(pmy_mesh_);
  InitializeParticles(pin, pmy_mesh_);

  if (global_variable::my_rank == 0) {
    std::cout << "particle_tracking_test initialized mode=" << test_mode
              << " field_profile=" << field_profile
              << " particle_layout=" << particle_layout << std::endl;
  }
}
