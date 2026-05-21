#ifndef PARTICLES_PARTICLES_HPP_
#define PARTICLES_PARTICLES_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particles.hpp
//  \brief definitions for Particles class

#include <cmath>
#include <map>
#include <memory>
#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "tasklist/task_list.hpp"
#include "bvals/bvals.hpp"

// forward declarations

// constants that enumerate ParticlesPusher options
enum class ParticlesPusher {drift, leap_frog, lagrangian_tracer, lagrangian_mc, drag};

// constants that enumerate ParticleTypes
enum class ParticleType {cosmic_ray, dust};

// constants that enumerate drag-particle force models
enum class DragParticlesModel {none, stopping_time, user};

// constants that enumerate particle-fluid interpolation/deposition stencils
enum class DragParticlesCoupling {host_cell, cloud_in_cell};

// Fixed-size interpolation/deposition stencil used by standard drag and pgen callbacks.
// The stencil is block-local and normalized so duplicated/clipped cells still conserve.
struct DragParticleStencil {
  int ncell;
  int i[8], j[8], k[8];
  Real w[8];
};

KOKKOS_INLINE_FUNCTION
int DragClampIndex(const int idx, const int lo, const int hi) {
  return (idx < lo) ? lo : ((idx > hi) ? hi : idx);
}

KOKKOS_INLINE_FUNCTION
void DragStencilAddCell(DragParticleStencil &st, const int i, const int j, const int k,
                        const Real w) {
  if (w == 0.0) {return;}
  for (int n=0; n<st.ncell; ++n) {
    if (st.i[n] == i && st.j[n] == j && st.k[n] == k) {
      st.w[n] += w;
      return;
    }
  }
  int n = st.ncell++;
  st.i[n] = i;
  st.j[n] = j;
  st.k[n] = k;
  st.w[n] = w;
}

KOKKOS_INLINE_FUNCTION
void DragStencilNormalize(DragParticleStencil &st) {
  Real wsum = 0.0;
  for (int n=0; n<st.ncell; ++n) {wsum += st.w[n];}
  if (wsum > 0.0) {
    for (int n=0; n<st.ncell; ++n) {st.w[n] /= wsum;}
  }
}

KOKKOS_INLINE_FUNCTION
void DragStencilAxis(const Real x, const Real xmin, const Real dx, const int is,
                     const int ie, const DragParticlesCoupling method, int idx[2],
                     Real wt[2], int &nidx) {
  if (method == DragParticlesCoupling::host_cell) {
    int i = static_cast<int>((x - xmin)/dx) + is;
    idx[0] = DragClampIndex(i, is, ie);
    wt[0] = 1.0;
    nidx = 1;
  } else {
    Real xi = (x - xmin)/dx + static_cast<Real>(is) - 0.5;
    int i0 = static_cast<int>(floor(xi));
    Real f = xi - static_cast<Real>(i0);
    idx[0] = DragClampIndex(i0, is, ie);
    idx[1] = DragClampIndex(i0 + 1, is, ie);
    wt[0] = 1.0 - f;
    wt[1] = f;
    nidx = 2;
  }
}

KOKKOS_INLINE_FUNCTION
void BuildDragParticleStencil(const DragParticlesCoupling method, const Real x1,
                              const Real x2, const Real x3, const Real x1min,
                              const Real x2min, const Real x3min, const Real dx1,
                              const Real dx2, const Real dx3, const int is,
                              const int ie, const int js, const int je, const int ks,
                              const int ke, const bool multi_d, const bool three_d,
                              DragParticleStencil &st) {
  st.ncell = 0;
  int ii[2], jj[2], kk[2];
  Real wi[2], wj[2], wk[2];
  int ni, nj, nk;
  DragStencilAxis(x1, x1min, dx1, is, ie, method, ii, wi, ni);
  if (multi_d) {
    DragStencilAxis(x2, x2min, dx2, js, je, method, jj, wj, nj);
  } else {
    jj[0] = js;
    wj[0] = 1.0;
    nj = 1;
  }
  if (three_d) {
    DragStencilAxis(x3, x3min, dx3, ks, ke, method, kk, wk, nk);
  } else {
    kk[0] = ks;
    wk[0] = 1.0;
    nk = 1;
  }
  for (int k=0; k<nk; ++k) {
    for (int j=0; j<nj; ++j) {
      for (int i=0; i<ni; ++i) {
        DragStencilAddCell(st, ii[i], jj[j], kk[k], wi[i]*wj[j]*wk[k]);
      }
    }
  }
  DragStencilNormalize(st);
}

template <typename PrimView>
KOKKOS_INLINE_FUNCTION
void InterpolateDragVelocity(const PrimView &w0, const int m,
                             const DragParticleStencil &st, const bool multi_d,
                             const bool three_d, Real &v1, Real &v2, Real &v3) {
  v1 = 0.0;
  v2 = 0.0;
  v3 = 0.0;
  for (int n=0; n<st.ncell; ++n) {
    v1 += st.w[n]*w0(m,IVX,st.k[n],st.j[n],st.i[n]);
    if (multi_d) {v2 += st.w[n]*w0(m,IVY,st.k[n],st.j[n],st.i[n]);}
    if (three_d) {v3 += st.w[n]*w0(m,IVZ,st.k[n],st.j[n],st.i[n]);}
  }
}

template <typename ConsView>
KOKKOS_INLINE_FUNCTION
void DepositDragBackreaction(const ConsView &u0, const int m,
                             const DragParticleStencil &st, const Real inv_vol,
                             const Real dm1, const Real dm2, const Real dm3,
                             const Real denergy, const bool multi_d, const bool three_d,
                             const bool include_energy) {
  for (int n=0; n<st.ncell; ++n) {
    Real coeff = st.w[n]*inv_vol;
    Kokkos::atomic_add(&u0(m,IM1,st.k[n],st.j[n],st.i[n]), coeff*dm1);
    if (multi_d) {
      Kokkos::atomic_add(&u0(m,IM2,st.k[n],st.j[n],st.i[n]), coeff*dm2);
    }
    if (three_d) {
      Kokkos::atomic_add(&u0(m,IM3,st.k[n],st.j[n],st.i[n]), coeff*dm3);
    }
    if (include_energy) {
      Kokkos::atomic_add(&u0(m,IEN,st.k[n],st.j[n],st.i[n]), coeff*denergy);
    }
  }
}

//----------------------------------------------------------------------------------------
//! \struct ParticlesTaskIDs
//  \brief container to hold TaskIDs of all particles tasks

struct ParticlesTaskIDs {
  TaskID push;
  TaskID newgid;
  TaskID count;
  TaskID irecv;
  TaskID sendp;
  TaskID recvp;
  TaskID csend;
  TaskID crecv;
  TaskID drag;
};

namespace particles {

//----------------------------------------------------------------------------------------
//! \class Particles

class Particles {
  friend class ParticlesBoundaryValues;
 public:
  Particles(MeshBlockPack *ppack, ParameterInput *pin);
  ~Particles();

  // data
  ParticleType particle_type;
  int nprtcl_thispack;             // number of particles this MeshBlockPack
  int nrdata, nidata;
//  DvceArray1D<int>  prtcl_gid;     // GID of MeshBlock containing each par
//  DvceArray2D<Real> prtcl_pos;     // positions
//  DvceArray2D<Real> prtcl_vel;     // velocities
  DvceArray2D<Real> prtcl_rdata;   // real number properties each particle (x,v,etc.)
  DvceArray2D<int>  prtcl_idata;   // integer properties each particle (gid, tag, etc.)
  Real dtnew;

  ParticlesPusher pusher;
  DragParticlesModel drag_model;
  bool drag_enabled;
  bool drag_backreaction;
  bool drag_include_energy;
  bool drag_orbital_terms;
  DragParticlesCoupling drag_interpolation;
  DragParticlesCoupling drag_deposition;
  Real drag_stopping_time;
  Real drag_particle_mass;
  Real drag_omega0;
  Real drag_qshear;

  // Boundary communication buffers and functions for particles
  ParticlesBoundaryValues *pbval_part;

  // container to hold names of TaskIDs
  ParticlesTaskIDs id;

  // functions...
  void CreateParticleTags(ParameterInput *pin);
  void AssembleTasks(std::map<std::string, std::shared_ptr<TaskList>> tl);
  TaskStatus Push(Driver *pdriver, int stage);
  TaskStatus NewGID(Driver *pdriver, int stage);
  TaskStatus SendCnt(Driver *pdriver, int stage);
  TaskStatus InitRecv(Driver *pdriver, int stage);
  TaskStatus SendP(Driver *pdriver, int stage);
  TaskStatus RecvP(Driver *pdriver, int stage);
  TaskStatus ClearSend(Driver *pdriver, int stage);
  TaskStatus ClearRecv(Driver *pdriver, int stage);
  TaskStatus ApplyDrag(Driver *pdriver, int stage);
  void ApplyStoppingTimeDrag(const Real bdt);
  void RemapAfterAMR();

 private:
  MeshBlockPack* pmy_pack;  // ptr to MeshBlockPack containing this Particles
};

} // namespace particles
#endif // PARTICLES_PARTICLES_HPP_
