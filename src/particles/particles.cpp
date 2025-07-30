//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particles.cpp
//! \brief implementation of Particles class constructor and assorted other functions

#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>
#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "particles.hpp"
#include "units/units.hpp"
#include "utils/sn_scheduler.hpp"

namespace particles {
//----------------------------------------------------------------------------------------
// constructor, initializes data structures and parameters

Particles::Particles(MeshBlockPack *ppack, ParameterInput *pin) :
    pmy_pack(ppack) {
  // check this is at least a 2D problem
  if (pmy_pack->pmesh->one_d) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Particle module only works in 2D/3D" <<std::endl;
    std::exit(EXIT_FAILURE);
  }

  // initialize vector for star particles
  std::vector<std::array<Real, 9>> particle_list;

  // select particle type
  {
    std::string ptype = pin->GetString("particles","particle_type");
    if (ptype.compare("cosmic_ray") == 0) {

      particle_type = ParticleType::cosmic_ray;

      // read number of particles per cell, and calculate number of particles this pack
      Real ppc = pin->GetOrAddReal("particles","ppc",1.0);

      // compute number of particles as real number, since ppc can be < 1
      auto &indcs = pmy_pack->pmesh->mb_indcs;
      int ncells = indcs.nx1*indcs.nx2*indcs.nx3;
      Real r_npart = ppc*static_cast<Real>((pmy_pack->nmb_thispack)*ncells);
      nprtcl_thispack = static_cast<int>(r_npart); // then cast to integer
						   
    } else if (ptype.compare("star") == 0) {

      particle_type = ParticleType::star;
      nprtcl_thispack = 0; // initialize to zero

      // Load particles from file
      std::string particle_file = pin->GetString("particles","star_particle_file");
      std::ifstream infile(particle_file);
      if (!infile) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Unable to open particle file: " << particle_file
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }

      // Read particle positions from file
      int nmb = pmy_pack->nmb_thispack;
      auto &size = pmy_pack->pmb->mb_size;

      // Skip header comments
      while (infile.peek() == '#') {
        infile.ignore(1000, '\n');
      }

      std::array<Real, 8> p;
      while (infile >> p[0] >> p[1] >> p[2] >> p[3] >>
		       p[4] >> p[5] >> p[6] >> p[7]) {
        for (int m=0; m<nmb; ++m) {
          // Loop over particles to see which are in this meshblock
          if (p[0] > size.h_view(m).x1min && p[0] <= size.h_view(m).x1max &&
              p[1] > size.h_view(m).x2min && p[1] <= size.h_view(m).x2max &&
              p[2] > size.h_view(m).x3min && p[2] <= size.h_view(m).x3max) {
            // Add particle to the list if it is within the mesh block bounds
            std::array<Real, 9> new_particle = {p[0], p[1], p[2], 
		                                p[3], p[4], p[5], 
						p[6], p[7], static_cast<Real>(m)};
            particle_list.push_back(new_particle);
          }
        }
      }
      infile.close();
      nprtcl_thispack = particle_list.size(); 

      // Print number of particles loaded
      std::cout << "Loaded " << nprtcl_thispack 
	        << " star particles from file " << particle_file << std::endl;
    }
  }

  // select pusher algorithm
  {
    std::string ppush = pin->GetString("particles","pusher");
    if (ppush.compare("drift") == 0) {
      pusher = ParticlesPusher::drift;
    } else if (ppush.compare("rk4_gravity") == 0) {
      pusher = ParticlesPusher::rk4_gravity;
      // load gravity constants
      r_scale   = pin->GetReal("potential", "r_scale");
      rho_scale = pin->GetReal("potential", "rho_scale");
      m_gal     = pin->GetReal("potential", "mass_gal");
      a_gal     = pin->GetReal("potential", "scale_gal");
      z_gal     = pin->GetReal("potential", "z_gal");
      r_200     = pin->GetReal("potential", "r_200");
      rho_mean  = pin->GetReal("potential", "rho_mean");
      par_grav_dx = pin->GetOrAddReal("particles", "grav_dx", 1e-6);
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle pusher must be specified in <particles> block"
                <<std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // set dimensions of particle arrays. Note particles only work in 2D/3D
  if (pmy_pack->pmesh->one_d) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Particles only work in 2D/3D, but 1D problem initialized" <<std::endl;
    std::exit(EXIT_FAILURE);
  }

  // stars must be 3D
  if (particle_type == ParticleType::star and not pmy_pack->pmesh->three_d) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Star particles only work in 3D" <<std::endl;
    std::exit(EXIT_FAILURE);
  }

  switch (particle_type) {
    case ParticleType::cosmic_ray:
      {
        int ndim=4;
        if (pmy_pack->pmesh->three_d) {ndim+=2;}
        nrdata = ndim;
        nidata = 2;
        break;
      }
    case ParticleType::star:
      {
        nrdata = 9;
        nidata = 3;
        break;
      }
    default:
      break;
  }
  
  Kokkos::realloc(prtcl_rdata, nrdata, nprtcl_thispack);
  Kokkos::realloc(prtcl_idata, nidata, nprtcl_thispack);

  // allocate boundary object
  pbval_part = new ParticlesBoundaryValues(this, pin);

  // Initialize Star particles if particle type is "star"
  {
    std::string ptype = pin->GetString("particles","particle_type");
    if (ptype.compare("star") == 0) {
      // Loop over mesh blocks in this pack
      int nmb = pmy_pack->nmb_thispack;
      auto &size = pmy_pack->pmb->mb_size;
      const int &gids = pmy_pack->gids;
      
      // Copy to device-accessible arrays
      // First create and populate a host view
      HostArray2D<Real> host_pos("host_positions", 9, nprtcl_thispack);
      for (size_t i = 0; i < nprtcl_thispack; ++i) {
	for (size_t j = 0; j < 9; ++j) {
          host_pos(j, i) = particle_list[i][j];
	}
      }

      // Then create the device view and copy data
      auto pos_data = Kokkos::create_mirror_view_and_copy(DevExeSpace(), host_pos);

      auto &pi = prtcl_idata;
      auto &pr = prtcl_rdata;
      int nrdata_ = nrdata;
      Real unit_time = pmy_pack->punit->time_cgs();

      // Initialize particles
      par_for("star_par", DevExeSpace(), 0, nprtcl_thispack-1,
      KOKKOS_LAMBDA(const int p) {   
        int m = static_cast<int>(pos_data(8, p));  
        pi(PGID,p) = gids + m;
	pi(2, p) = 0; // int to track # of SN
        pr(IPX,p)  = pos_data(0, p);
        pr(IPY,p)  = pos_data(1, p);
        pr(IPZ,p)  = pos_data(2, p);
        pr(IPVX,p) = pos_data(3, p);
        pr(IPVY,p) = pos_data(4, p);
        pr(IPVZ,p) = pos_data(5, p);
        pr(6, p) = pos_data(6, p); // time of creation of star particle
        pr(7, p) = pos_data(7, p); // mass of star particle
	pr(8, p) = GetNthSNTime(pr(7,p), pr(6,p), unit_time, 0); // time of next SN

        // Print particle initialization
        // Kokkos::printf("Initialized star particle %d in GID %d at position (%.2f, %.2f, %.2f)\n",
        //          p, gids + m, pos_data(0, p), pos_data(1, p), pos_data(2, p));
      });
    
      dtnew = std::min(size.h_view(0).dx1, size.h_view(0).dx2);
      dtnew = std::min(dtnew, size.h_view(0).dx3);
    }
  }
}

//----------------------------------------------------------------------------------------
// destructor

Particles::~Particles() {
}

//----------------------------------------------------------------------------------------
// CreatePaticleTags()
// Assigns tags to particles (unique integer).  Note that tracked particles are always
// those with tag numbers less than ntrack.

void Particles::CreateParticleTags(ParameterInput *pin) {
  std::string assign = pin->GetOrAddString("particles","assign_tag","index_order");

  // tags are assigned sequentially within this rank, starting at 0 with rank=0
  if (assign.compare("index_order") == 0) {
    int tagstart = 0;
    for (int n=1; n<=global_variable::my_rank; ++n) {
      tagstart += pmy_pack->pmesh->nprtcl_eachrank[n-1];
    }

    auto &pi = prtcl_idata;
    par_for("ptags",DevExeSpace(),0,(nprtcl_thispack-1),
    KOKKOS_LAMBDA(const int p) {
      pi(PTAG,p) = tagstart + p;
    });

  // tags are assigned sequentially across ranks
  } else if (assign.compare("rank_order") == 0) {
    int myrank = global_variable::my_rank;
    int nranks = global_variable::nranks;
    auto &pi = prtcl_idata;
    par_for("ptags",DevExeSpace(),0,(nprtcl_thispack-1),
    KOKKOS_LAMBDA(const int p) {
      pi(PTAG,p) = myrank + nranks*p;
    });

  // tag algorithm not recognized, so quit with error
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Particle tag assinment type = '" << assign << "' not recognized"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

} // namespace particles
