//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file shearing_box_cc.cpp
//! \brief functions to pack/send and recv/unpack boundary values for cell-centered (CC)
//! variables with shearing box boundaries.

#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <utility>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "bvals/bvals.hpp"
#include "shearing_box.hpp"
#include "remap_fluxes.hpp"

namespace {

constexpr int sbox_flux_tag_offset = 8;

int FluxShearCase(const int jr, const int nx2, const int ng) {
  if (jr < ng) return 1;
  if (jr < (nx2-ng)) return 2;
  return 3;
}

int NumFluxBuffers(const int scase) {
  return (scase == 2) ? 2 : 3;
}

int SendFluxJShift(const int n, const int ji, const int l, const int scase) {
  if (scase == 1) {
    return (n==0) ? (ji+l-1) : (l-1-ji);
  } else if (scase == 2) {
    return (n==0) ? (ji+l) : (l-1-ji);
  }
  return (n==0) ? (ji+l) : (l-2-ji);
}

int RecvFluxJShift(const int n, const int ji, const int l, const int scase) {
  return -SendFluxJShift(n,ji,l,scase);
}

void SetFluxJRanges(const int n, const int scase, const int js, const int je,
                    const int ng, const int nx2, const int jr,
                    std::pair<int,int> jsrc[3], std::pair<int,int> jdst[3]) {
  if (scase == 1) {
    if (n==0) {
      jsrc[0] = std::make_pair(js,js+ng-jr);
      jsrc[1] = std::make_pair(js,je+1);
      jsrc[2] = std::make_pair(je-(ng-1)-jr,je+1);
      jdst[0] = std::make_pair(je+1+jr,je+ng+1);
      jdst[1] = std::make_pair(js+jr,je+jr+1);
      jdst[2] = std::make_pair(js-ng,js+jr);
    } else {
      jsrc[0] = std::make_pair(js,js+ng+jr);
      jsrc[1] = std::make_pair(js,je+1);
      jsrc[2] = std::make_pair(je-(ng-1)+jr,je+1);
      jdst[0] = std::make_pair(je+1-jr,je+ng+1);
      jdst[1] = std::make_pair(js-jr,je-jr+1);
      jdst[2] = std::make_pair(js-ng,js-jr);
    }
  } else if (scase == 2) {
    if (n==0) {
      jsrc[0] = std::make_pair(js,je+ng-jr+1);
      jsrc[1] = std::make_pair(je-(ng-1)-jr,je+1);
      jdst[0] = std::make_pair(js+jr,je+ng+1);
      jdst[1] = std::make_pair(js-ng,js+jr);
    } else {
      jsrc[0] = std::make_pair(js,js+ng+jr);
      jsrc[1] = std::make_pair(js-ng+jr,je+1);
      jdst[0] = std::make_pair(je-jr+1,je+ng+1);
      jdst[1] = std::make_pair(js-ng,je-jr+1);
    }
  } else {
    const int dj = nx2-jr;
    if (n==0) {
      jsrc[0] = std::make_pair(js,js+ng+dj);
      jsrc[1] = std::make_pair(js,je+1);
      jsrc[2] = std::make_pair(je-(ng-1)+dj,je+1);
      jdst[0] = std::make_pair(je+1-dj,je+ng+1);
      jdst[1] = std::make_pair(js-dj,je-dj+1);
      jdst[2] = std::make_pair(js-ng,js-dj);
    } else {
      jsrc[0] = std::make_pair(js,js+ng-dj);
      jsrc[1] = std::make_pair(js,je+1);
      jsrc[2] = std::make_pair(je-(ng-1)-dj,je+1);
      jdst[0] = std::make_pair(je+1+dj,je+ng+1);
      jdst[1] = std::make_pair(js+dj,je+dj+1);
      jdst[2] = std::make_pair(js-ng,js+dj);
    }
  }
}

int FluxTagBuffer(const int n, const int l) {
  return sbox_flux_tag_offset + (n<<2) + l;
}

} // namespace

//----------------------------------------------------------------------------------------
// ShearingBoxCC derived class constructor:

ShearingBoxCC::ShearingBoxCC(MeshBlockPack *pp, ParameterInput *pin, int nvar) :
    ShearingBox(pp, pin) {
  // Allocate boundary buffers
  auto &indcs = pp->pmesh->mb_indcs;
  int ncells3 = indcs.nx3 + 2*indcs.ng;
  int ncells2 = indcs.nx2 + 2*indcs.ng;
  int ncells1 = indcs.ng;
  int nmb = std::max(1,std::max(nmb_x1bndry(0),nmb_x1bndry(1)));
  for (int n=0; n<2; ++n) {
    Kokkos::realloc(sendbuf[n].vars,nmb,ncells2,nvar,ncells3,ncells1);
    Kokkos::realloc(recvbuf[n].vars,nmb,ncells2,nvar,ncells3,ncells1);
    Kokkos::realloc(sendbuf[n].flux,nmb,ncells2,nvar,indcs.nx3,1);
    Kokkos::realloc(recvbuf[n].flux,nmb,ncells2,nvar,indcs.nx3,1);
#if MPI_PARALLEL_ENABLED
    sendbuf[n].flux_req = new MPI_Request[3*nmb];
    recvbuf[n].flux_req = new MPI_Request[3*nmb];
    for (int m=0; m<3*nmb; ++m) {
      sendbuf[n].flux_req[m] = MPI_REQUEST_NULL;
      recvbuf[n].flux_req[m] = MPI_REQUEST_NULL;
    }
#endif
  }
}

//----------------------------------------------------------------------------------------
// ShearingBoxCC derived class destructor:

ShearingBoxCC::~ShearingBoxCC() {
#if MPI_PARALLEL_ENABLED
  for (int n=0; n<2; ++n) {
    delete [] sendbuf[n].flux_req;
    delete [] recvbuf[n].flux_req;
  }
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void ShearingBox::PackAndSendCC()
//! \brief Apply shearing sheet BCs to cell-centered variables, including MPI
//! MPI communications. Both the inner_x1 and outer_x1 boundaries are updated.
//! Called on the physics_bcs task after purely periodic BC communication is finished.

TaskStatus ShearingBoxCC::PackAndSendCC(DvceArray5D<Real> &a, ReconstructionMethod rcon) {
  const auto &indcs = pmy_pack->pmesh->mb_indcs;
  const auto &ie = indcs.ie;
  const auto &js = indcs.js, &je = indcs.je;
  const auto &ks = indcs.ks, &ke = indcs.ke;
  const auto &ng = indcs.ng;

  // copy ghost zones at x1-faces into send buffer view
  // apply fractional cell offset to data in send buffers using conservative remap
  const int nvar = a.extent_int(1);  // TODO(@user): 2nd index from L must be NVAR
  const auto &mbsize = pmy_pack->pmb->mb_size;
  int kl=ks, ku=ke;
  if (pmy_pack->pmesh->three_d) {kl -= ng; ku += ng;}
  int nj = indcs.nx2 + 2*ng;
  const int &gids_ = pmy_pack->gids;
  const Real &yshear_ = yshear;
  const auto &x1bndry_mbgid_ = x1bndry_mbgid;
  auto &sbuf = sendbuf;
  int scr_lvl=0;
  size_t scr_size = ScrArray1D<Real>::shmem_size(nj) * 2;
  for (int n=0; n<2; ++n) {
    int nmb1 = nmb_x1bndry(n) - 1;
    par_for_outer("shrcc",DevExeSpace(),scr_size,scr_lvl,0,nmb1,0,(nvar-1),kl,ku,0,(ng-1),
    KOKKOS_LAMBDA(TeamMember_t member,const int m,const int v,const int k,const int i) {
      ScrArray1D<Real> a_(member.team_scratch(scr_lvl), nj); // 1D slice of data
      ScrArray1D<Real> flx(member.team_scratch(scr_lvl), nj); // "flux" at faces
      int mm = x1bndry_mbgid_.d_view(n,m) - gids_;

      // Load scratch array
      if (n==0) {
        par_for_inner(member, 0, nj-1, [&](const int j) {
          a_(j) = a(mm,v,k,j,i);
        });
      } else {
        par_for_inner(member, 0, nj-1, [&](const int j) {
          a_(j) = a(mm,v,k,j,(ie+1)+i);
        });
      }
      member.team_barrier();

      // compute fractional offset
      Real eps = fmod(yshear_,(mbsize.d_view(mm).dx2))/(mbsize.d_view(mm).dx2);
      if (n == 1) {eps *= -1.0;}

      // Compute "fluxes" at shifted cell faces
      switch (rcon) {
        case ReconstructionMethod::dc:
          DC_RemapFlx(member, js, (je+1), eps, a_, flx);
          break;
        case ReconstructionMethod::plm:
          PLM_RemapFlx(member, js, (je+1), eps, a_, flx);
          break;
        case ReconstructionMethod::ppm4:
        case ReconstructionMethod::ppmx:
        case ReconstructionMethod::wenoz:
          PPMX_RemapFlx(member, js, (je+1), eps, a_, flx);
          break;
        default:
          break;
      }
      member.team_barrier();

      // update data in send buffer with fractional shift
      par_for_inner(member, js, je, [&](const int j) {
        sbuf[n].vars(m,j,v,k,i) = a_(j) - (flx(j+1) - flx(j));
      });
    });
  }

  // shift data at x1 boundaries by integer number of cells.
  // Algorithm is broken into three steps: case1/2/3.
  //  * Case1 and case3 are when the integer shift (jr<ng), so that the sending MB
  //    overlaps the ghost cells of the two neighbors, and so requires copy/send
  //    to three separate target MBs.
  //  * Case2 is when the sending MB straddles the boundary between MBs, and so requires
  //    copy/send to only two target MBs.
  // Use deep copy if target MB on same rank, or MPI sends if not
  Kokkos::fence();
  const int &nx2 = indcs.nx2;
  bool no_errors=true;
  for (int n=0; n<2; ++n) {
    for (int m=0; m<nmb_x1bndry(n); ++m) {
      int gid = x1bndry_mbgid.h_view(n,m);
      int mm = gid - pmy_pack->gids;
      // Find integer and fractional number of grids over which offset extends.
      // This assumes every grid has same number of cells in x2-direction!
      int joffset  = static_cast<int>(yshear/(mbsize.h_view(mm).dx2));
      int ji = joffset/nx2;
      int jr = joffset - ji*nx2;

      if (jr < ng) {               //--- CASE 1 (in my nomenclature)
        int tgid, trank;
        std::pair<int,int> jsrc[3],jdst[3];
        if (n==0) {
          jsrc[0] = std::make_pair(js,js+ng-jr);
          jsrc[1] = std::make_pair(js,je+1);
          jsrc[2] = std::make_pair(je-(ng-1)-jr,je+1);
          jdst[0] = std::make_pair(je+1+jr,je+ng+1);
          jdst[1] = std::make_pair(js+jr,je+jr+1);
          jdst[2] = std::make_pair(js-ng,js+jr);
        } else {
          jsrc[0] = std::make_pair(js,js+ng+jr);
          jsrc[1] = std::make_pair(js,je+1);
          jsrc[2] = std::make_pair(je-(ng-1)+jr,je+1);
          jdst[0] = std::make_pair(je+1-jr,je+ng+1);
          jdst[1] = std::make_pair(js-jr,je-jr+1);
          jdst[2] = std::make_pair(js-ng,js-jr);
        }
        // ix1 boundary: send to (target-1) through (target+1)
        // ox1 boundary: send to (target-1) through (target+1)
        for (int l=0; l<3; ++l) {
          int jshift;
          if (n==0) {jshift = ji+l-1;} else {jshift = l-1-ji;} // offset of target
          FindTargetMB(gid,jshift,tgid,trank);
          if (trank == global_variable::my_rank) {
            int tm = TargetIndex(n,tgid);
            using Kokkos::ALL;
            auto src = subview(sendbuf[n].vars,m, jsrc[l],ALL,ALL,ALL);
            auto dst = subview(recvbuf[n].vars,tm,jdst[l],ALL,ALL,ALL);
            deep_copy(DevExeSpace(), dst, src);
#if MPI_PARALLEL_ENABLED
          } else {
            using Kokkos::ALL;
            auto send_ptr = subview(sendbuf[n].vars,m,jsrc[l],ALL,ALL,ALL);
            // create tag using GID of *receiving* MeshBlock
            int lid = tgid - pmy_pack->pmesh->gids_eachrank[trank];
            int tag = CreateBvals_MPI_Tag(lid, ((n<<2) | l));
            int data_size = send_ptr.size();
            int ierr = MPI_Isend(send_ptr.data(), data_size, MPI_ATHENA_REAL, trank, tag,
                                 comm_sbox, &(sendbuf[n].vars_req[3*m + l]));
            if (ierr != MPI_SUCCESS) {no_errors=false;}
#endif
          }
        }
      } else if (jr < (nx2-ng)) {  //--- CASE 2
        int tgid, trank;
        std::pair<int,int> jsrc[2],jdst[2];
        if (n==0) {
          jsrc[0] = std::make_pair(js,je+ng-jr+1);
          jsrc[1] = std::make_pair(je-(ng-1)-jr,je+1);
          jdst[0] = std::make_pair(js+jr,je+ng+1);
          jdst[1] = std::make_pair(js-ng,js+jr);
        } else {
          jsrc[0] = std::make_pair(js,js+ng+jr);
          jsrc[1] = std::make_pair(js-ng+jr,je+1);
          jdst[0] = std::make_pair(je-jr+1,je+ng+1);
          jdst[1] = std::make_pair(js-ng,je-jr+1);
        }
        // ix1 boundary: send to (target  ) through (target+1)
        // ox1 boundary: send to (target-1) through (target  )
        for (int l=0; l<2; ++l) {
          int jshift;
          if (n==0) {jshift = ji+l;} else {jshift = l-1-ji;}
          FindTargetMB(gid,jshift,tgid,trank);
          if (trank == global_variable::my_rank) {
            int tm = TargetIndex(n,tgid);
            using Kokkos::ALL;
            auto src = subview(sendbuf[n].vars,m, jsrc[l],ALL,ALL,ALL);
            auto dst = subview(recvbuf[n].vars,tm,jdst[l],ALL,ALL,ALL);
            deep_copy(DevExeSpace(), dst, src);
#if MPI_PARALLEL_ENABLED
          } else {
            using Kokkos::ALL;
            auto send_ptr = subview(sendbuf[n].vars,m,jsrc[l],ALL,ALL,ALL);
            // create tag using GID of *receiving* MeshBlock
            int lid = tgid - pmy_pack->pmesh->gids_eachrank[trank];
            int tag = CreateBvals_MPI_Tag(lid, ((n<<2) | l));
            int data_size = send_ptr.size();
            int ierr = MPI_Isend(send_ptr.data(), data_size, MPI_ATHENA_REAL, trank, tag,
                                 comm_sbox, &(sendbuf[n].vars_req[3*m + l]));
            if (ierr != MPI_SUCCESS) {no_errors=false;}
#endif
          }
        }
      } else {                     //--- CASE 3
        int tgid, trank;
        std::pair<int,int> jsrc[3],jdst[3];
        if (n==0) {
          jsrc[0] = std::make_pair(js,js+ng+(nx2-jr));
          jsrc[1] = std::make_pair(js,je+1);
          jsrc[2] = std::make_pair(je-(ng-1)+(nx2-jr),je+1);
          jdst[0] = std::make_pair(je+1-(nx2-jr),je+ng+1);
          jdst[1] = std::make_pair(js-(nx2-jr),je-(nx2-jr)+1);
          jdst[2] = std::make_pair(js-ng,js-(nx2-jr));
        } else {
          jsrc[0] = std::make_pair(js,js+ng-(nx2-jr));
          jsrc[1] = std::make_pair(js,je+1);
          jsrc[2] = std::make_pair(je-(ng-1)-(nx2-jr),je+1);
          jdst[0] = std::make_pair(je+1+(nx2-jr),je+ng+1);
          jdst[1] = std::make_pair(js+(nx2-jr),je+(nx2-jr)+1);
          jdst[2] = std::make_pair(js-ng,js+(nx2-jr));
        }
        // ix1 boundary: send to (target  ) through (target+2)
        // ox1 boundary: send to (target-2) through (target  )
        for (int l=0; l<3; ++l) {
          int jshift;
          if (n==0) {jshift = ji+l;} else {jshift = l-2-ji;}
          FindTargetMB(gid,jshift,tgid,trank);
          if (trank == global_variable::my_rank) {
            int tm = TargetIndex(n,tgid);
            using Kokkos::ALL;
            auto src = subview(sendbuf[n].vars,m, jsrc[l],ALL,ALL,ALL);
            auto dst = subview(recvbuf[n].vars,tm,jdst[l],ALL,ALL,ALL);
            deep_copy(DevExeSpace(), dst, src);
#if MPI_PARALLEL_ENABLED
          } else {
            using Kokkos::ALL;
            auto send_ptr = subview(sendbuf[n].vars,m,jsrc[l],ALL,ALL,ALL);
            // create tag using GID of *receiving* MeshBlock
            int lid = tgid - pmy_pack->pmesh->gids_eachrank[trank];
            int tag = CreateBvals_MPI_Tag(lid, ((n<<2) | l));
            int data_size = send_ptr.size();
            int ierr = MPI_Isend(send_ptr.data(), data_size, MPI_ATHENA_REAL, trank, tag,
                                 comm_sbox, &(sendbuf[n].vars_req[3*m + l]));
#endif
          }
        }
      }
    }
  }
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
       << std::endl << "MPI error in posting sends" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \!fn void ShearingBoxCC::RecvAndUnpackCC()
//! \brief Check MPI communication of boundary buffers for CC variables have finished,
//! then copy buffers into ghost zones. Shift has already been performed in
//! PackAndSendCC() function

TaskStatus ShearingBoxCC::RecvAndUnpackCC(DvceArray5D<Real> &a) {
  // create local references for variables in kernel
  const auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int &ng = indcs.ng;
#if MPI_PARALLEL_ENABLED
  //----- STEP 1: check that recv boundary buffer communications have all completed
  const int &nx2 = indcs.nx2;
  bool bflag = false;
  bool no_errors=true;
  for (int n=0; n<2; ++n) {
    for (int m=0; m<nmb_x1bndry(n); ++m) {
      int gid = x1bndry_mbgid.h_view(n,m);
      int mm = gid - pmy_pack->gids;
      // Find integer and fractional number of grids over which offset extends.
      // This assumes every grid has same number of cells in x2-direction!
      int joffset  = static_cast<int>(yshear/(pmy_pack->pmb->mb_size.h_view(mm).dx2));
      int ji = joffset/nx2;
      int jr = joffset - ji*nx2;

      if (jr < ng) {               //--- CASE 1 (in my nomenclature)
        // ix1 boundary: receive from (target+1) through (target-1)
        // ox1 boundary: receive from (target+1) through (target-1)
        for (int l=0; l<3; ++l) {
          int jshift;
          if (n==0) {jshift = -(ji+l-1);} else {jshift = -(l-1-ji);} // offset of sender
          int sgid, srank;
          FindTargetMB(gid,jshift,sgid,srank);
          if (srank != global_variable::my_rank) {
            int test;
            int ierr = MPI_Test(&(recvbuf[n].vars_req[3*m + l]),&test,MPI_STATUS_IGNORE);
            if (ierr != MPI_SUCCESS) {no_errors=false;}
            if (!(static_cast<bool>(test))) {bflag = true;}
          }
        }
      } else if (jr < (nx2-ng)) {  //--- CASE 2
        // ix1 boundary: receive from (target  ) through (target-1)
        // ox1 boundary: receive from (target+1) through (target  )
        for (int l=0; l<2; ++l) {
          int jshift;
          if (n==0) {jshift = -(ji+l);} else {jshift = -(l-1-ji);} // offset of sender
          int sgid, srank;
          FindTargetMB(gid,jshift,sgid,srank);
          if (srank != global_variable::my_rank) {
            int test;
            int ierr = MPI_Test(&(recvbuf[n].vars_req[3*m + l]),&test,MPI_STATUS_IGNORE);
            if (ierr != MPI_SUCCESS) {no_errors=false;}
            if (!(static_cast<bool>(test))) {bflag = true;}
          }
        }
      } else {                     //--- CASE 3
        // ix1 boundary: send to (target  ) through (target+2)
        // ox1 boundary: send to (target-2) through (target  )
        for (int l=0; l<3; ++l) {
          int jshift;
          if (n==0) {jshift = -(ji+l);} else {jshift = -(l-2-ji);} // offset of sender
          int sgid, srank;
          FindTargetMB(gid,jshift,sgid,srank);
          if (srank != global_variable::my_rank) {
            int test;
            int ierr = MPI_Test(&(recvbuf[n].vars_req[3*m + l]),&test,MPI_STATUS_IGNORE);
            if (ierr != MPI_SUCCESS) {no_errors=false;}
            if (!(static_cast<bool>(test))) {bflag = true;}
          }
        }
      }
    }
  }
  // Quit if MPI error detected
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "MPI error in testing non-blocking receives"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // exit if recv boundary buffer communications have not completed
  if (bflag) {return TaskStatus::incomplete;}
#endif

  //----- STEP 2: communications have all completed, so unpack and apply shift
  // copy recv buffer view into ghost zones at x1-faces
  const int nvar = a.extent_int(1);  // TODO(@user): 2nd index from L must be NVAR
  const int &ie = indcs.ie;
  int kl=indcs.ks, ku=indcs.ke;
  if (pmy_pack->pmesh->three_d) {kl -= ng; ku += ng;}
  int nj = indcs.nx2 + 2*ng;
  const int &gids_ = pmy_pack->gids;
  const auto &x1bndry_mbgid_ = x1bndry_mbgid;
  auto &rbuf = recvbuf;
  int scr_lvl=0;
  size_t scr_size = ScrArray1D<Real>::shmem_size(nj) * 3;
  for (int n=0; n<2; ++n) {
    int nmb1 = nmb_x1bndry(n) - 1;
    par_for_outer("shrcc",DevExeSpace(),scr_size,scr_lvl,0,nmb1,0,(nvar-1),kl,ku,0,(ng-1),
    KOKKOS_LAMBDA(TeamMember_t member,const int m,const int v,const int k,const int i) {
      int mm = x1bndry_mbgid_.d_view(n,m) - gids_;
      if (n==0) {
        par_for_inner(member, 0, nj-1, [&](const int j) {
          a(mm,v,k,j,i) = rbuf[n].vars(m,j,v,k,i);
        });
        member.team_barrier();
      } else {
        par_for_inner(member, 0, nj-1, [&](const int j) {
          a(mm,v,k,j,(ie+1)+i) = rbuf[n].vars(m,j,v,k,i);
        });
        member.team_barrier();
      }
    });
  }

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus ShearingBoxCC::InitFluxRecv()
//! \brief Post receives for opposite-face radial flux reconciliation.

TaskStatus ShearingBoxCC::InitFluxRecv() {
#if MPI_PARALLEL_ENABLED
  bool no_errors=true;
  const auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int &js = indcs.js, &je = indcs.je;
  const int &ng = indcs.ng;
  const int &nx2 = indcs.nx2;
  using Kokkos::ALL;
  for (int n=0; n<2; ++n) {
    for (int m=0; m<nmb_x1bndry(n); ++m) {
      int gid = x1bndry_mbgid.h_view(n,m);
      int mm = gid - pmy_pack->gids;
      int joffset = static_cast<int>(
          yshear/(pmy_pack->pmb->mb_size.h_view(mm).dx2));
      int ji = joffset/nx2;
      int jr = joffset - ji*nx2;
      int scase = FluxShearCase(jr,nx2,ng);
      int nbuff = NumFluxBuffers(scase);
      std::pair<int,int> jsrc[3],jdst[3];
      SetFluxJRanges(n,scase,js,je,ng,nx2,jr,jsrc,jdst);
      for (int l=0; l<nbuff; ++l) {
        int jshift = RecvFluxJShift(n,ji,l,scase);
        int sgid, srank;
        FindShearPartnerMB(gid,jshift,sgid,srank);
        if (srank != global_variable::my_rank) {
          auto recv_ptr = subview(recvbuf[n].flux,m,jdst[l],ALL,ALL,ALL);
          int lid = gid - pmy_pack->pmesh->gids_eachrank[global_variable::my_rank];
          int tag = CreateBvals_MPI_Tag(lid,FluxTagBuffer(n,l));
          int ierr = MPI_Irecv(recv_ptr.data(),recv_ptr.size(),MPI_ATHENA_REAL,srank,tag,
                               comm_sbox,&(recvbuf[n].flux_req[3*m+l]));
          if (ierr != MPI_SUCCESS) no_errors=false;
        }
      }
    }
  }
  if (!no_errors) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "MPI error posting shearing flux receives" << std::endl;
    std::exit(EXIT_FAILURE);
  }
#endif
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus ShearingBoxCC::PackAndSendFluxCC()
//! \brief Send each x1 flux to the opposite radial boundary with the integer shear shift.

TaskStatus ShearingBoxCC::PackAndSendFluxCC(DvceFaceFld5D<Real> &flx) {
  const auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int &is = indcs.is, &ie = indcs.ie;
  const int &js = indcs.js, &je = indcs.je;
  const int &ks = indcs.ks, &ke = indcs.ke;
  const int &ng = indcs.ng;
  const int &nx2 = indcs.nx2;
  const int nvar = flx.x1f.extent_int(1);
  const int &gids_ = pmy_pack->gids;
  const auto &x1bndry_mbgid_ = x1bndry_mbgid;
  auto &sbuf = sendbuf;

  for (int n=0; n<2; ++n) {
    const int sn = 1-n;
    const int nmb1 = nmb_x1bndry(sn)-1;
    const int i = (sn==0) ? is : (ie+1);
    par_for("shrflux_pack",DevExeSpace(),0,nmb1,0,nvar-1,ks,ke,js,je,
    KOKKOS_LAMBDA(const int m,const int v,const int k,const int j) {
      int mm = x1bndry_mbgid_.d_view(sn,m) - gids_;
      sbuf[n].flux(m,j,v,k-ks,0) = flx.x1f(mm,v,k,j,i);
    });
  }
  Kokkos::fence();

#if MPI_PARALLEL_ENABLED
  bool no_errors=true;
#endif
  using Kokkos::ALL;
  for (int n=0; n<2; ++n) {
    const int sn = 1-n;
    for (int m=0; m<nmb_x1bndry(sn); ++m) {
      int gid = x1bndry_mbgid.h_view(sn,m);
      int mm = gid - pmy_pack->gids;
      int joffset = static_cast<int>(
          yshear/(pmy_pack->pmb->mb_size.h_view(mm).dx2));
      int ji = joffset/nx2;
      int jr = joffset - ji*nx2;
      int scase = FluxShearCase(jr,nx2,ng);
      int nbuff = NumFluxBuffers(scase);
      std::pair<int,int> jsrc[3],jdst[3];
      SetFluxJRanges(n,scase,js,je,ng,nx2,jr,jsrc,jdst);
      for (int l=0; l<nbuff; ++l) {
        int jshift = SendFluxJShift(n,ji,l,scase);
        int tgid, trank;
        FindShearPartnerMB(gid,jshift,tgid,trank);
        if (trank == global_variable::my_rank) {
          int tm = TargetIndex(n,tgid);
          auto src = subview(sendbuf[n].flux,m,jsrc[l],ALL,ALL,ALL);
          auto dst = subview(recvbuf[n].flux,tm,jdst[l],ALL,ALL,ALL);
          deep_copy(DevExeSpace(),dst,src);
#if MPI_PARALLEL_ENABLED
        } else {
          auto send_ptr = subview(sendbuf[n].flux,m,jsrc[l],ALL,ALL,ALL);
          int lid = tgid - pmy_pack->pmesh->gids_eachrank[trank];
          int tag = CreateBvals_MPI_Tag(lid,FluxTagBuffer(n,l));
          int ierr = MPI_Isend(send_ptr.data(),send_ptr.size(),MPI_ATHENA_REAL,trank,tag,
                               comm_sbox,&(sendbuf[n].flux_req[3*m+l]));
          if (ierr != MPI_SUCCESS) no_errors=false;
#endif
        }
      }
    }
  }
#if MPI_PARALLEL_ENABLED
  if (!no_errors) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "MPI error posting shearing flux sends" << std::endl;
    std::exit(EXIT_FAILURE);
  }
#endif
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus ShearingBoxCC::RecvAndCorrectFluxCC()
//! \brief Remap the opposite-face flux and average it with the local radial flux.

TaskStatus ShearingBoxCC::RecvAndCorrectFluxCC(DvceFaceFld5D<Real> &flx,
                                              ReconstructionMethod rcon) {
  const auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int &is = indcs.is, &ie = indcs.ie;
  const int &js = indcs.js, &je = indcs.je;
  const int &ks = indcs.ks, &ke = indcs.ke;
  const int &ng = indcs.ng;
  const int &nx2 = indcs.nx2;
  const int nvar = flx.x1f.extent_int(1);
#if MPI_PARALLEL_ENABLED
  bool incomplete=false;
  bool no_errors=true;
  for (int n=0; n<2; ++n) {
    for (int m=0; m<nmb_x1bndry(n); ++m) {
      int gid = x1bndry_mbgid.h_view(n,m);
      int mm = gid - pmy_pack->gids;
      int joffset = static_cast<int>(
          yshear/(pmy_pack->pmb->mb_size.h_view(mm).dx2));
      int ji = joffset/nx2;
      int jr = joffset - ji*nx2;
      int scase = FluxShearCase(jr,nx2,ng);
      int nbuff = NumFluxBuffers(scase);
      for (int l=0; l<nbuff; ++l) {
        int jshift = RecvFluxJShift(n,ji,l,scase);
        int sgid, srank;
        FindShearPartnerMB(gid,jshift,sgid,srank);
        if (srank != global_variable::my_rank) {
          int test;
          int ierr = MPI_Test(&(recvbuf[n].flux_req[3*m+l]),&test,MPI_STATUS_IGNORE);
          if (ierr != MPI_SUCCESS) no_errors=false;
          if (!static_cast<bool>(test)) incomplete=true;
        }
      }
    }
  }
  if (!no_errors) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "MPI error testing shearing flux receives" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (incomplete) return TaskStatus::incomplete;
#endif

  const int nj = indcs.nx2 + 2*ng;
  const int &gids_ = pmy_pack->gids;
  const auto &x1bndry_mbgid_ = x1bndry_mbgid;
  const auto &mbsize = pmy_pack->pmb->mb_size;
  const Real &yshear_ = yshear;
  auto &rbuf = recvbuf;
  int scr_lvl=0;
  size_t scr_size = ScrArray1D<Real>::shmem_size(nj)*2;
  for (int n=0; n<2; ++n) {
    int nmb1 = nmb_x1bndry(n)-1;
    par_for_outer("shrflux_remap",DevExeSpace(),scr_size,scr_lvl,0,nmb1,0,nvar-1,
                  ks,ke,KOKKOS_LAMBDA(TeamMember_t member,const int m,const int v,
                                      const int k) {
      ScrArray1D<Real> a_(member.team_scratch(scr_lvl),nj);
      ScrArray1D<Real> flx_(member.team_scratch(scr_lvl),nj);
      int mm = x1bndry_mbgid_.d_view(n,m)-gids_;
      par_for_inner(member,0,nj-1,[&](const int j) {
        a_(j) = rbuf[n].flux(m,j,v,k-ks,0);
      });
      member.team_barrier();
      Real eps = fmod(yshear_,mbsize.d_view(mm).dx2)/mbsize.d_view(mm).dx2;
      if (n==1) eps *= -1.0;
      switch (rcon) {
        case ReconstructionMethod::dc:
          DC_RemapFlx(member,js,je+1,eps,a_,flx_);
          break;
        case ReconstructionMethod::plm:
          PLM_RemapFlx(member,js,je+1,eps,a_,flx_);
          break;
        case ReconstructionMethod::ppm4:
        case ReconstructionMethod::ppmx:
        case ReconstructionMethod::wenoz:
          PPMX_RemapFlx(member,js,je+1,eps,a_,flx_);
          break;
        default:
          break;
      }
      member.team_barrier();
      par_for_inner(member,js,je,[&](const int j) {
        rbuf[n].flux(m,j,v,k-ks,0) = a_(j) - (flx_(j+1)-flx_(j));
      });
    });
  }
  Kokkos::fence();

  for (int n=0; n<2; ++n) {
    int nmb1 = nmb_x1bndry(n)-1;
    int i = (n==0) ? is : (ie+1);
    par_for("shrflux_correct",DevExeSpace(),0,nmb1,0,nvar-1,ks,ke,js,je,
    KOKKOS_LAMBDA(const int m,const int v,const int k,const int j) {
      int mm = x1bndry_mbgid_.d_view(n,m)-gids_;
      flx.x1f(mm,v,k,j,i) =
          0.5*(flx.x1f(mm,v,k,j,i) + rbuf[n].flux(m,j,v,k-ks,0));
    });
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus ShearingBoxCC::ClearFluxRecv()

TaskStatus ShearingBoxCC::ClearFluxRecv() {
#if MPI_PARALLEL_ENABLED
  bool no_errors=true;
  const int &ng = pmy_pack->pmesh->mb_indcs.ng;
  const int &nx2 = pmy_pack->pmesh->mb_indcs.nx2;
  for (int n=0; n<2; ++n) {
    for (int m=0; m<nmb_x1bndry(n); ++m) {
      int gid = x1bndry_mbgid.h_view(n,m);
      int mm = gid - pmy_pack->gids;
      int joffset = static_cast<int>(
          yshear/(pmy_pack->pmb->mb_size.h_view(mm).dx2));
      int ji = joffset/nx2;
      int jr = joffset-ji*nx2;
      int scase = FluxShearCase(jr,nx2,ng);
      int nbuff = NumFluxBuffers(scase);
      for (int l=0; l<nbuff; ++l) {
        int jshift = RecvFluxJShift(n,ji,l,scase);
        int sgid, srank;
        FindShearPartnerMB(gid,jshift,sgid,srank);
        if (srank != global_variable::my_rank) {
          int ierr = MPI_Wait(&(recvbuf[n].flux_req[3*m+l]),MPI_STATUS_IGNORE);
          if (ierr != MPI_SUCCESS) no_errors=false;
        }
      }
    }
  }
  if (!no_errors) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "MPI error clearing shearing flux receives" << std::endl;
    std::exit(EXIT_FAILURE);
  }
#endif
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus ShearingBoxCC::ClearFluxSend()

TaskStatus ShearingBoxCC::ClearFluxSend() {
#if MPI_PARALLEL_ENABLED
  bool no_errors=true;
  const int &ng = pmy_pack->pmesh->mb_indcs.ng;
  const int &nx2 = pmy_pack->pmesh->mb_indcs.nx2;
  for (int n=0; n<2; ++n) {
    const int sn = 1-n;
    for (int m=0; m<nmb_x1bndry(sn); ++m) {
      int gid = x1bndry_mbgid.h_view(sn,m);
      int mm = gid - pmy_pack->gids;
      int joffset = static_cast<int>(
          yshear/(pmy_pack->pmb->mb_size.h_view(mm).dx2));
      int ji = joffset/nx2;
      int jr = joffset-ji*nx2;
      int scase = FluxShearCase(jr,nx2,ng);
      int nbuff = NumFluxBuffers(scase);
      for (int l=0; l<nbuff; ++l) {
        int jshift = SendFluxJShift(n,ji,l,scase);
        int tgid, trank;
        FindShearPartnerMB(gid,jshift,tgid,trank);
        if (trank != global_variable::my_rank) {
          int ierr = MPI_Wait(&(sendbuf[n].flux_req[3*m+l]),MPI_STATUS_IGNORE);
          if (ierr != MPI_SUCCESS) no_errors=false;
        }
      }
    }
  }
  if (!no_errors) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "MPI error clearing shearing flux sends" << std::endl;
    std::exit(EXIT_FAILURE);
  }
#endif
  return TaskStatus::complete;
}
