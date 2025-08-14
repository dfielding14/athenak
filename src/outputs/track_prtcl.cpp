
#include <algorithm>
#include <cstdio>      // fwrite(), fclose(), fopen(), fnprintf(), snprintf()
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>  // mkdir

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "particles/particles.hpp"
#include "outputs.hpp"
#include "mhd/mhd.hpp"

//----------------------------------------------------------------------------------------
// ctor: also calls BaseTypeOutput base class constructor

TrackedParticleOutput::TrackedParticleOutput(ParameterInput *pin, Mesh *pm,
                                             OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  // create new directory for this output. Comments in binary.cpp constructor explain why
  mkdir("trk",0775);
  // allocate arrays
  npout_eachrank.resize(global_variable::nranks);
  ntrack = pin->GetInteger(op.block_name,"nparticles");
  // TODO(@user) improve guess below?
  ntrack_thisrank = ntrack;
}

//----------------------------------------------------------------------------------------
// TrackedParticleOutput::LoadOutputData()
// Copies data for tracked particles on this rank to host outpart array

void TrackedParticleOutput::LoadOutputData(Mesh *pm) {
  particles::Particles *pp = pm->pmb_pack->ppart;

  // Load data for tracked particles on this rank into new device array
  DualArray1D<TrackedParticleData> tracked_prtcl("d_trked",ntrack_thisrank * pp->nspecies);
  int npart = pm->nprtcl_thisrank;
  npart = pm->pmb_pack->ppart->nprtcl_thispack;
  auto &pr = pm->pmb_pack->ppart->prtcl_rdata;
  auto &pi = pm->pmb_pack->ppart->prtcl_idata;
  auto ntrack_ = ntrack;
  auto nspecies_ = pp->nspecies;
  auto nprtcl_total_ = pm->nprtcl_total;
  int species_offset = nprtcl_total_ / nspecies_;

  auto &indcs = pm->mb_indcs;
  const int is = indcs.is;
  const int js = indcs.js;
  const int ks = indcs.ks;
  const int ie = indcs.ie;
  const int je = indcs.je;
  const int ke = indcs.ke;
  const int ngh = indcs.ng;
  auto &mbsize = pm->pmb_pack->pmb->mb_size;
  auto &gids = pm->pmb_pack->gids;
  auto &bcc = pm->pmb_pack->pmhd->bcc0;

  Kokkos::View<int> counter("counter");
  Kokkos::deep_copy(counter, 0);
  par_for("part_trackout",DevExeSpace(),0,(npart-1), KOKKOS_LAMBDA(const int p) {
    if (pi(PTAG,p) % species_offset < ntrack_) {
      int index = Kokkos::atomic_fetch_add(&counter(),1);
      tracked_prtcl.d_view(index).tag = pi(PTAG,p)%species_offset + ntrack_*pi(PSP,p);
      tracked_prtcl.d_view(index).x   = pr(IPX,p);
      tracked_prtcl.d_view(index).y   = pr(IPY,p);
      tracked_prtcl.d_view(index).z   = pr(IPZ,p);
      tracked_prtcl.d_view(index).vx  = pr(IPVX,p);
      tracked_prtcl.d_view(index).vy  = pr(IPVY,p);
      tracked_prtcl.d_view(index).vz  = pr(IPVZ,p);
      tracked_prtcl.d_view(index).Bx  = pr(IPBX,p);
      tracked_prtcl.d_view(index).By  = pr(IPBY,p);
      tracked_prtcl.d_view(index).Bz  = pr(IPBZ,p);

      // Now compute curvature.
      int m = pi(PGID,p) - gids;
      int i = (pr(IPX,p) - mbsize.d_view(m).x1min)/mbsize.d_view(m).dx1 + is;
      int j = (pr(IPY,p) - mbsize.d_view(m).x2min)/mbsize.d_view(m).dx2 + js;
      int k = (pr(IPZ,p) - mbsize.d_view(m).x3min)/mbsize.d_view(m).dx3 + ks;

      Real &bx = bcc(m,IBX,k,j,i);
      Real &by = bcc(m,IBY,k,j,i);
      Real &bz = bcc(m,IBZ,k,j,i);
      Real B_mag_sq = ( bx*bx + by*by + bz*bz);

      // Calculate b_hat vector
      Real b1 = bcc(m,IBX,k,j,i)/sqrt(B_mag_sq);
      Real b2 = bcc(m,IBY,k,j,i)/sqrt(B_mag_sq);
      Real b3 = bcc(m,IBZ,k,j,i)/sqrt(B_mag_sq);

      // calculate b_hat vector at i +/- 1
      Real B_mag_ip1 = sqrt( bcc(m,IBX,k,j,i+1)*bcc(m,IBX,k,j,i+1)
                           + bcc(m,IBY,k,j,i+1)*bcc(m,IBY,k,j,i+1)
                           + bcc(m,IBZ,k,j,i+1)*bcc(m,IBZ,k,j,i+1));
      Real b1_ip1 = bcc(m,IBX,k,j,i+1)/B_mag_ip1;
      Real b2_ip1 = bcc(m,IBY,k,j,i+1)/B_mag_ip1;
      Real b3_ip1 = bcc(m,IBZ,k,j,i+1)/B_mag_ip1;

      Real B_mag_im1 = sqrt( bcc(m,IBX,k,j,i-1)*bcc(m,IBX,k,j,i-1)
                           + bcc(m,IBY,k,j,i-1)*bcc(m,IBY,k,j,i-1)
                           + bcc(m,IBZ,k,j,i-1)*bcc(m,IBZ,k,j,i-1));
      Real b1_im1 = bcc(m,IBX,k,j,i-1)/B_mag_im1;
      Real b2_im1 = bcc(m,IBY,k,j,i-1)/B_mag_im1;
      Real b3_im1 = bcc(m,IBZ,k,j,i-1)/B_mag_im1;

      // calculate b_hat vector at j +/- 1
      Real B_mag_jp1 = sqrt( bcc(m,IBX,k,j+1,i)*bcc(m,IBX,k,j+1,i)
                           + bcc(m,IBY,k,j+1,i)*bcc(m,IBY,k,j+1,i)
                           + bcc(m,IBZ,k,j+1,i)*bcc(m,IBZ,k,j+1,i));
      Real b1_jp1 = bcc(m,IBX,k,j+1,i)/B_mag_jp1;
      Real b2_jp1 = bcc(m,IBY,k,j+1,i)/B_mag_jp1;
      Real b3_jp1 = bcc(m,IBZ,k,j+1,i)/B_mag_jp1;

      Real B_mag_jm1 = sqrt( bcc(m,IBX,k,j-1,i)*bcc(m,IBX,k,j-1,i)
                           + bcc(m,IBY,k,j-1,i)*bcc(m,IBY,k,j-1,i)
                           + bcc(m,IBZ,k,j-1,i)*bcc(m,IBZ,k,j-1,i));
      Real b1_jm1 = bcc(m,IBX,k,j-1,i)/B_mag_jm1;
      Real b2_jm1 = bcc(m,IBY,k,j-1,i)/B_mag_jm1;
      Real b3_jm1 = bcc(m,IBZ,k,j-1,i)/B_mag_jm1;

      // calculate b_hat vector at k +/- 1
      Real B_mag_kp1 = sqrt( bcc(m,IBX,k+1,j,i)*bcc(m,IBX,k+1,j,i)
                           + bcc(m,IBY,k+1,j,i)*bcc(m,IBY,k+1,j,i)
                           + bcc(m,IBZ,k+1,j,i)*bcc(m,IBZ,k+1,j,i));
      Real b1_kp1 = bcc(m,IBX,k+1,j,i)/B_mag_kp1;
      Real b2_kp1 = bcc(m,IBY,k+1,j,i)/B_mag_kp1;
      Real b3_kp1 = bcc(m,IBZ,k+1,j,i)/B_mag_kp1;

      Real B_mag_km1 = sqrt( bcc(m,IBX,k-1,j,i)*bcc(m,IBX,k-1,j,i)
                           + bcc(m,IBY,k-1,j,i)*bcc(m,IBY,k-1,j,i)
                           + bcc(m,IBZ,k-1,j,i)*bcc(m,IBZ,k-1,j,i));
      Real b1_km1 = bcc(m,IBX,k-1,j,i)/B_mag_km1;
      Real b2_km1 = bcc(m,IBY,k-1,j,i)/B_mag_km1;
      Real b3_km1 = bcc(m,IBZ,k-1,j,i)/B_mag_km1;

      // Central differencing of b_hat vector
      Real db1_dx1 = (b1_ip1 - b1_im1)/(2.0*mbsize.d_view(m).dx1);
      Real db2_dx1 = (b2_ip1 - b2_im1)/(2.0*mbsize.d_view(m).dx1);
      Real db3_dx1 = (b3_ip1 - b3_im1)/(2.0*mbsize.d_view(m).dx1);

      Real db1_dx2 = (b1_jp1 - b1_jm1)/(2.0*mbsize.d_view(m).dx2);
      Real db2_dx2 = (b2_jp1 - b2_jm1)/(2.0*mbsize.d_view(m).dx2);
      Real db3_dx2 = (b3_jp1 - b3_jm1)/(2.0*mbsize.d_view(m).dx2);

      Real db1_dx3 = (b1_kp1 - b1_km1)/(2.0*mbsize.d_view(m).dx3);
      Real db2_dx3 = (b2_kp1 - b2_km1)/(2.0*mbsize.d_view(m).dx3);
      Real db3_dx3 = (b3_kp1 - b3_km1)/(2.0*mbsize.d_view(m).dx3);

      // Calculate curvature = |b_hat dot nabla b_hat|
      Real curv1 = b1*db1_dx1 + b2*db1_dx2 + b3*db1_dx3;
      Real curv2 = b1*db2_dx1 + b2*db2_dx2 + b3*db2_dx3;
      Real curv3 = b1*db3_dx1 + b2*db3_dx2 + b3*db3_dx3;

      //Now central difference Bmag vector:
      Real Bmag_dx1 = (B_mag_ip1 - B_mag_im1)/(2.0*mbsize.d_view(m).dx1);
      Real Bmag_dx2 = (B_mag_jp1 - B_mag_jm1)/(2.0*mbsize.d_view(m).dx2);
      Real Bmag_dx3 = (B_mag_kp1 - B_mag_km1)/(2.0*mbsize.d_view(m).dx3);

      tracked_prtcl.d_view(index).K1  = curv1;
      tracked_prtcl.d_view(index).K2  = curv2;
      tracked_prtcl.d_view(index).K3  = curv3;

      tracked_prtcl.d_view(index).dB1  = Bmag_dx1;
      tracked_prtcl.d_view(index).dB2  = Bmag_dx2;
      tracked_prtcl.d_view(index).dB3  = Bmag_dx3;

/*
      // Calculate gradB tensor
      Real dBx_dx = (bcc(m,IBX,k,j,i+1) - bcc(m,IBX,k,j,i-1))/(2.0*mbsize.d_view(m).dx1);
      Real dBx_dy = (bcc(m,IBX,k,j+1,i) - bcc(m,IBX,k,j-1,i))/(2.0*mbsize.d_view(m).dx2);
      Real dBx_dz = (bcc(m,IBX,k+1,j,i) - bcc(m,IBX,k-1,j,i))/(2.0*mbsize.d_view(m).dx3);


      Real dBy_dx = (bcc(m,IBY,k,j,i+1) - bcc(m,IBY,k,j,i-1))/(2.0*mbsize.d_view(m).dx1);
      Real dBy_dy = (bcc(m,IBY,k,j+1,i) - bcc(m,IBY,k,j-1,i))/(2.0*mbsize.d_view(m).dx2);
      Real dBy_dz = (bcc(m,IBY,k+1,j,i) - bcc(m,IBY,k-1,j,i))/(2.0*mbsize.d_view(m).dx3);

      Real dBz_dx = (bcc(m,IBZ,k,j,i+1) - bcc(m,IBZ,k,j,i-1))/(2.0*mbsize.d_view(m).dx1);
      Real dBz_dy = (bcc(m,IBZ,k,j+1,i) - bcc(m,IBZ,k,j-1,i))/(2.0*mbsize.d_view(m).dx2);
      Real dBz_dz = (bcc(m,IBZ,k+1,j,i) - bcc(m,IBZ,k-1,j,i))/(2.0*mbsize.d_view(m).dx3);

      Real BdotGradB_x = (bx * dBx_dx + by * dBx_dy + bz * dBx_dz);
      Real BdotGradB_y = (bx * dBy_dx + by * dBy_dy + bz * dBy_dz);
      Real BdotGradB_z = (bx * dBz_dx + by * dBz_dy + bz * dBz_dz);

      Real Identity_minus_bhat_bhat_xx = 1.0 - bx*bx/B_mag_squared;
      Real Identity_minus_bhat_bhat_xy = 0.0 - bx*by/B_mag_squared;
      Real Identity_minus_bhat_bhat_xz = 0.0 - bx*bz/B_mag_squared;

      Real Identity_minus_bhat_bhat_yx = 0.0 - by*bx/B_mag_squared;
      Real Identity_minus_bhat_bhat_yy = 1.0 - by*by/B_mag_squared;
      Real Identity_minus_bhat_bhat_yz = 0.0 - by*bz/B_mag_squared;

      Real Identity_minus_bhat_bhat_zx = 0.0 - bz*bx/B_mag_squared;
      Real Identity_minus_bhat_bhat_zy = 0.0 - bz*by/B_mag_squared;
      Real Identity_minus_bhat_bhat_zz = 1.0 - bz*bz/B_mag_squared;

      // Calculate curvature which is |(B.gradB).(I - bhat bhat)/B^2|
      Real curv1 = (BdotGradB_x * Identity_minus_bhat_bhat_xx
		      + BdotGradB_y * Identity_minus_bhat_bhat_yx
		      + BdotGradB_z * Identity_minus_bhat_bhat_zx
		      );
      Real curv2 = (BdotGradB_x * Identity_minus_bhat_bhat_xy
		      + BdotGradB_y * Identity_minus_bhat_bhat_yy
		      + BdotGradB_z * Identity_minus_bhat_bhat_zy
		      );
      Real curv3 = (BdotGradB_x * Identity_minus_bhat_bhat_xz
		      + BdotGradB_y * Identity_minus_bhat_bhat_yz
		      + BdotGradB_z * Identity_minus_bhat_bhat_zz
		      );
      Real bdotb = sqrt(curv1*curv1 + curv2*curv2 + curv3*curv3)/B_mag_squared;
      tracked_prtcl.d_view(index).Kmag  = bdotb;
      */

      // Now compute current
      Real j1 = (bcc(m,IBZ,k,j+1,i) - bcc(m,IBZ,k,j-1,i))/(2.0*mbsize.d_view(m).dx2);
      Real j2 = -(bcc(m,IBZ,k,j,i+1) - bcc(m,IBZ,k,j,i-1))/(2.0*mbsize.d_view(m).dx1);
      Real j3 =  (bcc(m,IBY,k,j,i+1) - bcc(m,IBY,k,j,i-1))/(2.0*mbsize.d_view(m).dx1);
      j1 -= (bcc(m,IBY,k+1,j,i) - bcc(m,IBY,k-1,j,i))/(2.0*mbsize.d_view(m).dx3);
      j2 += (bcc(m,IBX,k+1,j,i) - bcc(m,IBX,k-1,j,i))/(2.0*mbsize.d_view(m).dx3);
      j3 -= (bcc(m,IBX,k,j+1,i) - bcc(m,IBX,k,j-1,i))/(2.0*mbsize.d_view(m).dx2);

      tracked_prtcl.d_view(index).jmag  = std::sqrt(j1*j1+j2*j2+j3*j3);
    }
  });
  Kokkos::deep_copy(npout, counter);
  // share number of tracked particles to be output across all ranks
  npout_eachrank[global_variable::my_rank] = npout;
#if MPI_PARALLEL_ENABLED
  MPI_Allgather(&npout, 1, MPI_INT, npout_eachrank.data(), 1, MPI_INT, MPI_COMM_WORLD);
#endif
  tracked_prtcl.resize(npout);
  // sync tracked particle device array with host
  tracked_prtcl.template modify<DevExeSpace>();
  tracked_prtcl.template sync<HostMemSpace>();

  // copy host view into host outpart array
  Kokkos::realloc(outpart, npout);
  Kokkos::deep_copy(outpart, tracked_prtcl.h_view);
}

//----------------------------------------------------------------------------------------
//! \fn void TrackedParticleOutput:::WriteOutputFile(Mesh *pm)
//! \brief Cycles over all tracked particles on this rank and writes ouput data
//! With MPI, all particles are written to the same file.

void TrackedParticleOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {

  if (ncycle_buffer > 1) {
    WriteOutputFileWithBuffer(pm, pin);
    return;
  }
  int big_end = IsBigEndian(); // =1 on big endian machine

  // create filename: "trk/file_basename".trk
  std::string fname;
  fname.assign("trk/");
  fname.append(out_params.file_basename);
  fname.append(".trk");

  // Root process opens/creates file and appends string
  if (global_variable::my_rank == 0) {
    std::stringstream msg;
    msg << std::endl << "# AthenaK tracked particle data at time= " << pm->time
        << "  nranks= " << global_variable::nranks
        << "  cycle=" << pm->ncycle
        << "  ntracked_prtcls=" << ntrack << std::endl;
    FILE *pfile;
    if ((pfile = std::fopen(fname.c_str(),"a")) == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "Output file '" << fname << "' could not be opened" <<std::endl;
      exit(EXIT_FAILURE);
    }
    std::fprintf(pfile,"%s \n",msg.str().c_str());
    std::fclose(pfile);
  }

  // Now all ranks open file and append data
  IOWrapper partfile;
  partfile.Open(fname.c_str(), IOWrapper::FileMode::append);
  std::size_t header_offset = partfile.GetPosition();

  // allocate 1D vector of floats used to convert and output particle data
  float *data = new float[6*npout];
  // Loop over particles, load positions into data[]
  for (int p=0; p<npout; ++p) {
    data[ 6*p   ] = static_cast<float>(outpart(p).x);
    data[(6*p)+1] = static_cast<float>(outpart(p).y);
    data[(6*p)+2] = static_cast<float>(outpart(p).z);
    data[(6*p)+3] = static_cast<float>(outpart(p).vx);
    data[(6*p)+4] = static_cast<float>(outpart(p).vy);
    data[(6*p)+5] = static_cast<float>(outpart(p).vz);
  }

  // calculate local data offset
  std::vector<int> rank_offset(global_variable::nranks, 0);
  int npout_min = pm->nprtcl_eachrank[0];
  for (int n=1; n<global_variable::nranks; ++n) {
    rank_offset[n] = rank_offset[n-1] + npout_eachrank[n-1];
    npout_min = std::min(npout_min, npout_eachrank[n]);
  }

  // Write tracked particle data collectively over minimum shared number of prtcls
  for (int p=0; p<npout_min; ++p) {
    // offset computed assuming tags run 0...(ntrack-1) sequentially
    std::size_t myoffset = header_offset + 6 * outpart(p).tag * sizeof(float);
    // Write particle positions collectively for minimum number of particles across ranks
    if (partfile.Write_any_type_at_all(&(data[6*p]),6,myoffset,"float") != 6) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "particle data not written correctly to tracked particle file"
          << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  // Write particle positions individually for remaining particles on each rank
  for (int p=npout_min; p<npout; ++p) {
    // offset computed assuming tags run 0...(ntrack-1) sequentially
    std::size_t myoffset = header_offset + 6 * outpart(p).tag * sizeof(float);
    // Write particle positions collectively for minimum number of particles across ranks
    if (partfile.Write_any_type_at(&(data[6*p]),6,myoffset,"float") != 6) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "particle data not written correctly to tracked particle file"
          << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  // close the output file and clean up
  partfile.Close();
  delete[] data;

  // increment counters
  if (out_params.last_time < 0.0 || out_params.last_time == pm->time) {
    out_params.last_time = pm->time;
  } else {
    out_params.last_time += out_params.dt;
  }
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void TrackedParticleOutput:::WriteOutputFile(Mesh *pm)
//! \brief Cycles over all tracked particles on this rank and writes ouput data
//! With MPI, all particles are written to the same file.

void TrackedParticleOutput::WriteOutputFileWithBuffer(Mesh *pm, ParameterInput *pin) {
  int big_end = IsBigEndian(); // =1 on big endian machine
  int pp;
  int nspecies_ =  pm->pmb_pack->ppart->nspecies;
  double time_cycle = pm->time;
  // Loop over particles, load positions into data[]
  for (int p=nout_thisrank; p<nout_thisrank + npout; ++p) {
    pp = p - nout_thisrank;
    particle_buffer[(19*p)+0]  = static_cast<float>(icycle_buffer);
    particle_buffer[(19*p)+1]  = static_cast<float>(outpart(pp).tag);
    particle_buffer[(19*p)+2]  = static_cast<float>(time_cycle);
    particle_buffer[(19*p)+3]  = static_cast<float>(outpart(pp).x);
    particle_buffer[(19*p)+4]  = static_cast<float>(outpart(pp).y);
    particle_buffer[(19*p)+5]  = static_cast<float>(outpart(pp).z);
    particle_buffer[(19*p)+6]  = static_cast<float>(outpart(pp).vx);
    particle_buffer[(19*p)+7]  = static_cast<float>(outpart(pp).vy);
    particle_buffer[(19*p)+8]  = static_cast<float>(outpart(pp).vz);
    particle_buffer[(19*p)+9]  = static_cast<float>(outpart(pp).Bx);
    particle_buffer[(19*p)+10]  = static_cast<float>(outpart(pp).By);
    particle_buffer[(19*p)+11] = static_cast<float>(outpart(pp).Bz);
    particle_buffer[(19*p)+12] = static_cast<float>(outpart(pp).K1);
    particle_buffer[(19*p)+13] = static_cast<float>(outpart(pp).K2);
    particle_buffer[(19*p)+14] = static_cast<float>(outpart(pp).K3);
    particle_buffer[(19*p)+15] = static_cast<float>(outpart(pp).dB1);
    particle_buffer[(19*p)+16] = static_cast<float>(outpart(pp).dB2);
    particle_buffer[(19*p)+17] = static_cast<float>(outpart(pp).dB3);
    particle_buffer[(19*p)+18] = static_cast<float>(outpart(pp).jmag);
  }

  nout_thisrank += npout;
  icycle_buffer += 1;

  if (icycle_buffer >= ncycle_buffer){
    if(track_single_file_per_rank){
    std::string fname;
    char rank_dir[20];
    std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d/", global_variable::my_rank);
    fname = std::string("trk/") + std::string(rank_dir) + out_params.file_basename  + ".trk";
      std::stringstream msg;
      msg << std::endl << "# AthenaK tracked particle data at time= " << pm->time
	  << "  rank= " << global_variable::my_rank
          << "  nout= "	<< nout_thisrank
          << "  nranks= " << global_variable::nranks
          << "  cycle=" << pm->ncycle
          << "  buffer cycles" << ncycle_buffer
          << "  ntracked_prtcls=" << ntrack << std::endl;
      FILE *pfile;
      if ((pfile = std::fopen(fname.c_str(),"a")) == nullptr) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "Output file '" << fname << "' could not be opened" <<std::endl;
        exit(EXIT_FAILURE);
      }
      std::fprintf(pfile,"%s \n",msg.str().c_str());
      std::fclose(pfile);

    // Now all ranks open file and append data
    IOWrapper partfile;
    partfile.Open(fname.c_str(), IOWrapper::FileMode::append);
    std::size_t header_offset = partfile.GetPosition();

    std::size_t datasize = sizeof(float);
    int ptag_, icycle_;
    // Write tracked particle data collectively over minimum shared number of prtcls
    for (int p=0; p<nout_thisrank; ++p) {
      // offset computed assuming tags run 0...(ntrack-1) sequentially
      icycle_ = static_cast<int>(particle_buffer[19*p]);
      ptag_ = static_cast<int>(particle_buffer[19*p+1]);
      std::size_t myoffset = header_offset +  18*p* datasize;
      // Write particle positions collectively for minimum number of particles across ranks
      if (partfile.Write_any_type_at_all(&(particle_buffer[19*p+1]),18,myoffset,"float") != 18) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << "particle data not written correctly to tracked particle file"
            << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    // close the output file and clean up
    partfile.Close();
    icycle_buffer = 0;
    nout_thisrank = 0;

    }
    else{
    // create filename: "trk/file_basename".trk
    std::string fname;
    fname.assign("trk/");
    fname.append(out_params.file_basename);
    fname.append(".trk");

    // Root process opens/creates file and appends string
    if (global_variable::my_rank == 0) {
      std::stringstream msg;
      msg << std::endl << "# AthenaK tracked particle data at time= " << pm->time
          << "  nranks= " << global_variable::nranks
          << "  cycle=" << pm->ncycle
	  << "  buffer cycles" << ncycle_buffer
          << "  ntracked_prtcls=" << ntrack << std::endl;
      FILE *pfile;
      if ((pfile = std::fopen(fname.c_str(),"a")) == nullptr) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "Output file '" << fname << "' could not be opened" <<std::endl;
        exit(EXIT_FAILURE);
      }
      std::fprintf(pfile,"%s \n",msg.str().c_str());
      std::fclose(pfile);
    }

    // Now all ranks open file and append data
    IOWrapper partfile;
    partfile.Open(fname.c_str(), IOWrapper::FileMode::append);
    std::size_t header_offset = partfile.GetPosition();

    std::size_t datasize = sizeof(float);
    int ptag_, icycle_;
    // Write tracked particle data collectively over minimum shared number of prtcls
    for (int p=0; p<nout_thisrank; ++p) {
      // offset computed assuming tags run 0...(ntrack-1) sequentially
      icycle_ = static_cast<int>(particle_buffer[19*p]);
      ptag_ = static_cast<int>(particle_buffer[19*p+1]);
      std::size_t myoffset = header_offset + icycle_ * ntrack * nspecies_ * 18 * datasize  +  18*ptag_* datasize;
      // Write particle positions collectively for minimum number of particles across ranks
      if (partfile.Write_any_type_at_all(&(particle_buffer[19*p+1]),18,myoffset,"float") != 18) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << "particle data not written correctly to tracked particle file"
            << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    // close the output file and clean up
    partfile.Close();
    icycle_buffer = 0;
    nout_thisrank = 0;
  }
  }
  // increment counters
  if (out_params.last_time < 0.0 || out_params.last_time == pm->time) {
    out_params.last_time = pm->time;
  } else {
    out_params.last_time += out_params.dt;
  }
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
  return;
}
